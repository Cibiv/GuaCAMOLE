import argparse
import os
import sys
from datetime import datetime
from time import localtime, strftime

import numpy as np
import pandas as pd

import guacamole.library as lib
from guacamole.bracken import (
    parse_kmer_distributions,
    parse_kraken_report,
    redistribute_reads,
)
from guacamole.distributions import (
    get_genome_length,
    get_sample_dists,
    load_gc_distributions,
)
from guacamole.io import check_excessive_removal, write_output
from guacamole.plotting import save_diagnostic_plots

np.seterr(divide='ignore', invalid='ignore')


def parse_args():
    parser = argparse.ArgumentParser(
        description='Run GuaCAMOLE for GC aware species abundance estimation from metagenomic data'
    )

    parser.add_argument('--kraken_report', metavar='kraken_report', type=str,
                        help='.report file from Kraken2 output')
    parser.add_argument('--kraken_file', metavar='kraken_file', type=str,
                        help='.kraken file from Kraken2 output')
    parser.add_argument('--kraken_db', metavar='kraken_db', type=str,
                        help='Path to Kraken2 Database')
    parser.add_argument('--read_len', metavar='read_len', type=int,
                        help='read length')
    parser.add_argument('--output', metavar='output', type=str,
                        help='name of the output file')
    parser.add_argument('--read_files', metavar='read_files', type=str,
                        help='path(s) to read file(s)', nargs='+')
    parser.add_argument('--threshold', metavar='threshold', type=int,
                        help='Minimum number of reads for a species to estimate its abundance',
                        default=500)
    parser.add_argument('--level', metavar='level', type=str,
                        help='Taxonomy level to be quantified (S=species, G=genus)',
                        default='S')
    parser.add_argument('--length_correction', action='store_true',
                        help='genome size correction')
    parser.add_argument('--plot', action='store_true',
                        help='True if detailed plots should be generated')
    parser.add_argument('--fp_cycles', metavar='fp_cycles', type=int,
                        help='Number of iterations for false positive removal', default=5)
    parser.add_argument('--reg_weight', metavar='reg_weight', type=float,
                        help='Regularization strength [between 0 and 1]', default=0.01)
    parser.add_argument('--fragment_len', metavar='fragment_len', type=int,
                        help='length of the fragment if paired end and known', default=None)
    parser.add_argument('--fasta', action='store_true',
                        help='True if reads are in fasta format, false if fastq')
    parser.add_argument('--quantiles', metavar='quantiles', type=float,
                        help='min and max quantiles for GC distributions',
                        nargs=2, default=[0.05, 0.95])
    parser.add_argument('--verbose', action='store_true',
                        help='Write per-iteration metrics (abundances, efficiencies, residuals) to output files '
                             'for diagnostic visualization')

    return parser.parse_args()


def main():
    args = parse_args()

    report = args.kraken_report
    kraken = args.kraken_file
    kraken_db = args.kraken_db
    read_len = args.read_len
    output = args.output
    threshold = args.threshold
    level = args.level
    length_correction = args.length_correction
    plot = args.plot
    fp_cycles = args.fp_cycles
    reg_weight = args.reg_weight
    fragment_len = args.fragment_len
    fasta = args.fasta
    quantiles = args.quantiles
    verbose = args.verbose

    if len(args.read_files) == 2:
        fastq1, fastq2 = args.read_files
    elif len(args.read_files) == 1:
        fastq1 = args.read_files[0]
        fastq2 = None
    else:
        raise ValueError(f"--read_files expects 1 or 2 files, got {len(args.read_files)}")

    nbin = 100
    kmer_distr = os.path.join(kraken_db, f'database{read_len}mers.kmer_distrib')

    time_start = strftime("%m-%d-%Y %H:%M:%S", localtime())
    sys.stdout.write(f"PROGRAM START TIME: {time_start}\n")

    # ---- Parse kraken report and build taxonomy tree ----
    (
        root_node, all_nodes, lvl_taxids, map2lvl_taxids, _desired_lvl_dct,
        est_reads_dct, _total_reads, _kept_reads, _ignored_reads,
        _n_lvl_total, _n_lvl_est, _n_lvl_del,
    ) = parse_kraken_report(report, level, threshold)

    # ---- Parse kmer distributions ----
    kmer_distr_dict = parse_kmer_distributions(kmer_distr, lvl_taxids, map2lvl_taxids)

    # ---- Load GC distributions (keep raw DF if length_correction is needed) ----
    gc_dists, insert, gc_dists_raw = load_gc_distributions(
        kraken_db, nbin, read_len, fragment_len, return_raw=True
    )

    # ---- Initialize distribution weights ----
    distribution_weights = {taxid: {} for taxid in map2lvl_taxids}

    # ---- Redistribute reads (Bracken-style) ----
    reference_nodes, est_reads_dct = redistribute_reads(
        root_node, level, kmer_distr_dict, map2lvl_taxids, lvl_taxids,
        gc_dists, distribution_weights,
    )

    if len(lvl_taxids) == 0:
        return

    # ---- Build reference distributions ----
    refdists_dict = get_sample_dists(
        sdists=reference_nodes, weights=distribution_weights,
        map2lvl_taxids=map2lvl_taxids, lvl_taxids=lvl_taxids,
        reference=True, gc_dists=gc_dists, kraken_db=kraken_db,
    )

    cols_sorted = sorted(int(col) for col in refdists_dict)
    cols_str = [str(col) for col in cols_sorted]
    refdists = np.column_stack([refdists_dict[col] for col in cols_str])
    ref_normalized = refdists.copy()

    # ---- Create sample distributions ----
    taxids = np.array([int(node) for node in all_nodes], dtype=int)
    taxids = np.sort(np.append(taxids, 1))  # add root node

    sample_dist_file = f'sample_bin_{nbin}.dist'
    if not os.path.exists(sample_dist_file):
        cur_time = strftime("%m-%d-%Y %H:%M:%S", localtime())
        print(f"Starting to create sample distribution at {cur_time}")
        if fragment_len is not None and fastq2 is not None:
            lib.create_sample_dist_avg(
                fastq1=fastq1, fastq2=fastq2, txids=taxids, kraken=kraken, fasta=fasta
            )
        else:
            lib.create_sample_dist(
                fastq1=fastq1, fastq2=fastq2, txids=taxids, kraken=kraken, fasta=fasta
            )
        cur_time = strftime("%m-%d-%Y %H:%M:%S", localtime())
        print(f"Done at {cur_time}")

    sdist = pd.read_csv(sample_dist_file, sep='\t', header=0)
    first_col = sdist.columns[0]
    if first_col.startswith('# '):
        sdist = sdist.rename(columns={first_col: first_col[2:]})

    sample_dists_dict = get_sample_dists(
        sdists=sdist, weights=distribution_weights,
        map2lvl_taxids=map2lvl_taxids, lvl_taxids=lvl_taxids,
        kraken_db=kraken_db,
    )
    s_dists = np.column_stack([sample_dists_dict[col] for col in cols_str])

    # ---- Apply length correction if requested ----
    lengths = None
    if length_correction:
        lengths = get_genome_length(
            kraken_db=kraken_db, nbin=nbin, lvl_taxids=lvl_taxids,
            est_reads_dct=est_reads_dct, map2lvl_taxids=map2lvl_taxids,
            insert=insert, read_len=read_len, gc_dists_raw=gc_dists_raw,
        )
        refdists = refdists * lengths

    # ---- Handle NaN taxids ----
    cols_sorted = np.array(cols_sorted, dtype=int)
    invalid_mask = np.isnan(s_dists.sum(axis=0)) | np.isnan(refdists.sum(axis=0))
    nan_taxids = cols_sorted[invalid_mask]
    if len(nan_taxids) != 0:
        valid_mask = ~invalid_mask
        s_dists = s_dists[:, valid_mask]
        refdists = refdists[:, valid_mask]
        cols_sorted = cols_sorted[valid_mask]
        ref_normalized = ref_normalized[:, valid_mask]
        if length_correction:
            lengths = lengths[valid_mask]

    # ---- Normalize and save intermediate distributions ----
    ref_sum = float(np.sum(refdists))
    if ref_sum > 0:
        refdists = refdists * (float(np.sum(s_dists)) / ref_sum)

    np.savetxt(f'ref_bin_{nbin}_input.dist', refdists)
    np.savetxt(f'sample_bin_{nbin}_input.dist', s_dists)

    # ---- Run GuaCAMOLE core algorithm ----
    ab, taxon_removal_cycle, efficiencies, diagnostics = lib.corrected_abundances(
        s_dists,
        refdists,
        quantiles=quantiles,
        fp_cycles=fp_cycles,
        taxids=cols_sorted,
        plot=False,
        reg_weight=reg_weight,
        verbose=verbose,
        return_diagnostics=True,
    )

    ab[np.isinf(ab)] = np.nan

    # ---- Compute Bracken estimates ----
    ab_br = np.array([float(lvl_taxids[str(taxid)][3]) for taxid in cols_sorted], dtype=float)
    sum_all_reads = float(np.sum(ab_br))
    bracken_ab = ab_br / sum_all_reads if sum_all_reads > 0 else np.zeros_like(ab_br)

    if length_correction:
        corr_ab = bracken_ab / lengths
        bracken_ab = corr_ab / np.sum(corr_ab)
        if len(nan_taxids) != 0:
            print(
                f"ERROR: For taxid(s) {nan_taxids} reference distributions could not be created"
                "(known bug), therefore length"
                " correction is not possible. Please re-run with length_correction set to False."
            )

    # ---- Compute efficiency-based abundances ----
    covered_bins = np.isfinite(efficiencies) & (efficiencies > 0)
    valid_mass = ref_normalized[covered_bins].sum(axis=0)
    raw_avg_eff = ref_normalized[covered_bins].T @ efficiencies[covered_bins]
    avg_eff = np.divide(
        raw_avg_eff,
        valid_mass,
        out=np.full(len(cols_sorted), np.nan, dtype=float),
        where=valid_mass > 0,
    )
    ab_eff_raw = np.divide(
        bracken_ab,
        avg_eff,
        out=np.full(len(cols_sorted), np.nan, dtype=float),
        where=np.isfinite(avg_eff) & (avg_eff > 0),
    )
    eff_sum = float(np.nansum(ab_eff_raw))
    ab_efficiencies = (ab_eff_raw / eff_sum) if eff_sum > 0 else np.full_like(ab_eff_raw, np.nan)
    gc_content = np.argmax(ref_normalized, axis=0).astype(float)

    # ---- Handle NaN taxids in output ----
    if len(nan_taxids) != 0:
        print(
            f"WARNING: For taxid(s) {nan_taxids} reference distributions could not be created, "
            "therefore GuaCAMOLE estimates are set to nan for those!"
        )
        n_nan = len(nan_taxids)
        ab = np.append(ab, np.full(n_nan, np.nan))
        cols_sorted = np.append(cols_sorted, nan_taxids)
        taxon_removal_cycle = np.append(taxon_removal_cycle, np.full(n_nan, np.nan))
        ab_efficiencies = np.append(ab_efficiencies, np.full(n_nan, np.nan))
        gc_content = np.append(gc_content, np.full(n_nan, np.nan))

        ab_br = np.array([float(lvl_taxids[str(taxid)][3]) for taxid in cols_sorted], dtype=float)
        sum_all_reads = float(np.sum(ab_br))
        bracken_ab = ab_br / sum_all_reads if sum_all_reads > 0 else np.zeros_like(ab_br)

    # ---- Build results DataFrame ----
    ab_df = pd.DataFrame({
        'abundance_ls': ab,
        'abundance_br': bracken_ab,
        'taxid': cols_sorted,
        'taxon_removal_cycle': taxon_removal_cycle,
        'GC content': gc_content,
        'abundance_eff': ab_efficiencies,
    })
    ab_df = ab_df.set_index('taxid', drop=False)

    # ---- Check for excessive removal ----
    excessive_removal_detected, efficiencies = check_excessive_removal(ab_df, efficiencies)

    # ---- Write efficiencies ----
    max_eff = float(np.nanmax(efficiencies)) if np.any(np.isfinite(efficiencies)) else np.nan
    if np.isfinite(max_eff) and max_eff > 0:
        efficiencies = efficiencies / max_eff
    np.savetxt('efficiencies.txt', efficiencies)

    # ---- Write output ----
    write_output(output, lvl_taxids, ab_df, sum_all_reads, level, excessive_removal_detected)

    # ---- Generate diagnostic report PDF if requested ----
    if plot:
        pdf_path = save_diagnostic_plots(
            output_path=output,
            diagnostics=diagnostics,
            ab_df=ab_df,
            efficiencies=efficiencies,
            excessive_removal_detected=excessive_removal_detected,
            include_iterations=verbose,
        )
        print(f"Saved diagnostic plots to {pdf_path}")

    # ---- Report timing ----
    time_end = strftime("%m-%d-%Y %H:%M:%S", localtime())
    sys.stdout.write(f"DONE AT: {time_end}\n")
    tdelta = (
        datetime.strptime(time_end, "%m-%d-%Y %H:%M:%S")
        - datetime.strptime(time_start, "%m-%d-%Y %H:%M:%S")
    )
    print(f"DURATION: {tdelta.total_seconds()} seconds")
    np.savetxt('duration_time.txt', np.array([tdelta.total_seconds()]))


if __name__ == "__main__":
    main()
