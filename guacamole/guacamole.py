import os
import sys
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from time import localtime, strftime
from datetime import datetime

import guacamole.library as lib
from guacamole.bracken import parse_kraken_report, parse_kmer_distributions, redistribute_reads
from guacamole.distributions import get_sample_dists, get_genome_length, load_gc_distributions
from guacamole.io import write_output, check_excessive_removal

np.seterr(divide='ignore', invalid='ignore')


def parse_args():
    parser = argparse.ArgumentParser(
        description='Run GuaCAMOLE for GC aware species abundance estimation from metagenomic data')

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
                        nargs=2, default=[0.025, 0.975])

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

    if len(args.read_files) == 2:
        fastq1 = args.read_files[0]
        fastq2 = args.read_files[1]
    elif len(args.read_files) == 1:
        fastq1 = args.read_files[0]
        fastq2 = None

    nbin = 100
    kmer_distr = os.path.join(kraken_db, 'database' + str(read_len) + 'mers.kmer_distrib')

    time_start = strftime("%m-%d-%Y %H:%M:%S", localtime())
    sys.stdout.write("PROGRAM START TIME: " + time_start + '\n')

    # ---- Parse kraken report and build taxonomy tree ----
    (root_node, all_nodes, lvl_taxids, map2lvl_taxids, desired_lvl_dct,
     est_reads_dct, total_reads, kept_reads, ignored_reads,
     n_lvl_total, n_lvl_est, n_lvl_del) = parse_kraken_report(report, level, threshold)

    # ---- Parse kmer distributions ----
    kmer_distr_dict = parse_kmer_distributions(kmer_distr, lvl_taxids, map2lvl_taxids)

    # ---- Load GC distributions ----
    gc_dists, insert = load_gc_distributions(kraken_db, nbin, read_len, fragment_len)

    # ---- Initialize distribution weights ----
    distribution_weights = {}
    for taxid in map2lvl_taxids:
        distribution_weights[taxid] = {}

    # ---- Redistribute reads (Bracken-style) ----
    reference_nodes, est_reads_dct = redistribute_reads(
        root_node, level, kmer_distr_dict, map2lvl_taxids, lvl_taxids,
        gc_dists, distribution_weights
    )

    if len(lvl_taxids) == 0:
        return None

    # ---- Build reference distributions ----
    refdists = get_sample_dists(
        sdists=reference_nodes, weights=distribution_weights,
        map2lvl_taxids=map2lvl_taxids, lvl_taxids=lvl_taxids,
        reference=True, gc_dists=gc_dists, kraken_db=kraken_db
    )

    df = pd.DataFrame(refdists)
    cols_int = [int(col) for col in df.columns]
    cols_sorted = sorted(cols_int)
    cols_str = [str(col) for col in cols_sorted]
    df = df.reindex(cols_str, axis=1)
    refdists = np.array(df)
    ref_normalized = refdists.copy()

    # ---- Create sample distributions ----
    taxids = np.array([int(node) for node in all_nodes])
    taxids = np.append(taxids, 1)  # add root node
    taxids = np.sort(taxids)

    if not os.path.exists('sample_bin_' + str(nbin) + '.dist'):
        time = strftime("%m-%d-%Y %H:%M:%S", localtime())
        print("Starting to create sample distribution at " + time)
        if fragment_len is not None:
            lib.create_sample_dist_avg(
                fastq1=fastq1, fastq2=fastq2, txids=taxids, kraken=kraken, fasta=fasta
            )
        else:
            lib.create_sample_dist(
                fastq1=fastq1, fastq2=fastq2, txids=taxids, kraken=kraken, fasta=fasta
            )
        time = strftime("%m-%d-%Y %H:%M:%S", localtime())
        print("Done at " + time)

    sdist = pd.read_csv('sample_bin_' + str(nbin) + '.dist', sep='\t', header=0)
    newcols = np.array(sdist.columns)
    newcols[0] = newcols[0][2:]
    col_dct = dict(zip(np.array(sdist.columns), newcols))
    sdist.rename(columns=col_dct, inplace=True)

    sample_dists = get_sample_dists(
        sdists=sdist, weights=distribution_weights,
        map2lvl_taxids=map2lvl_taxids, lvl_taxids=lvl_taxids,
        kraken_db=kraken_db
    )

    df = pd.DataFrame(sample_dists)
    cols_int = [int(col) for col in df.columns]
    cols_sorted = sorted(cols_int)
    cols_str = [str(col) for col in cols_sorted]
    df = df.reindex(cols_str, axis=1)
    s_dists = np.array(df)

    # ---- Apply length correction if requested ----
    if length_correction:
        lengths = get_genome_length(
            kraken_db=kraken_db, nbin=nbin, lvl_taxids=lvl_taxids,
            est_reads_dct=est_reads_dct, map2lvl_taxids=map2lvl_taxids,
            insert=insert, read_len=read_len
        )
        refdists = refdists * lengths

    # ---- Handle NaN taxids ----
    nan_taxids = np.array(cols_sorted)[np.isnan(s_dists.sum(0))]
    if len(nan_taxids) != 0:
        non_nan = np.argwhere(~np.isnan(s_dists.sum(0)))
        non_nan = non_nan.reshape((non_nan.shape[0],))
        s_dists = s_dists[:, non_nan]
        refdists = refdists[:, non_nan]
        cols_sorted = np.array(cols_sorted)[non_nan]
        ref_normalized = ref_normalized[:, non_nan]
        if length_correction:
            lengths = lengths[non_nan]

    # ---- Normalize and save intermediate distributions ----
    refdists = refdists / np.sum(refdists)
    refdists = refdists * np.sum(s_dists)

    np.savetxt('ref_bin_' + str(nbin) + '_input.dist', refdists)
    np.savetxt('sample_bin_' + str(nbin) + '_input.dist', s_dists)

    # ---- Run GuaCAMOLE core algorithm ----
    norm, ref, sample = lib.normalize_gc_dists(
        'sample_bin_' + str(nbin) + '_input.dist',
        'ref_bin_' + str(nbin) + '_input.dist',
        quantiles=quantiles
    )

    if plot:
        lib.plot_dist2(norm, line=True)
        plt.savefig('Obs_vs_exp_tresh_' + str(threshold) + '.pdf')
        plt.close()

    ab, taxon_removal_cycle, efficiencies = lib.corrected_abundances(
        'sample_bin_' + str(nbin) + '_input.dist',
        'ref_bin_' + str(nbin) + '_input.dist',
        fp_cycles=fp_cycles,
        taxids=cols_sorted, plot=plot,
        reg_weight=reg_weight
    )

    ab[np.isinf(ab)] = np.nan

    # ---- Compute Bracken estimates ----
    sum_all_reads = 0
    ab_br = []
    for taxid in cols_sorted:
        [name, all_reads, lvl_reads, added_reads] = lvl_taxids[str(taxid)]
        sum_all_reads += float(added_reads)
        ab_br.append(float(added_reads))

    bracken_ab = [reads / sum_all_reads for reads in ab_br]

    if length_correction:
        corr_ab = bracken_ab / lengths
        bracken_ab = corr_ab / np.sum(corr_ab)
        if len(nan_taxids) != 0:
            print(f"ERROR: For taxid(s) {nan_taxids} reference distributions could not be created"
                  f"(known bug), therefore length"
                  f" correction is not possible. Please re-run with length_correction set to False.")

    # ---- Compute efficiency-based abundances ----
    avg_eff = np.matmul(ref_normalized.T, efficiencies.reshape(101, 1)).flatten()
    ab_efficiencies = (bracken_ab / avg_eff) / np.sum(bracken_ab / avg_eff)
    gc_content = np.argmax(ref_normalized, axis=0)

    # ---- Handle NaN taxids in output ----
    if len(nan_taxids) != 0:
        print(f"WARNING: For taxid(s) {nan_taxids} reference distributions could not be created, "
              f"therefore GuaCAMOLE estimates are set to nan for those!")
        ab = np.append(ab, np.repeat(np.nan, len(nan_taxids)))
        cols_sorted = np.append(cols_sorted, nan_taxids)
        taxon_removal_cycle = np.append(taxon_removal_cycle, np.repeat(np.nan, len(nan_taxids)))
        ab_efficiencies = np.append(ab_efficiencies, np.repeat(np.nan, len(nan_taxids)))
        gc_content = np.append(gc_content, np.repeat(np.nan, len(nan_taxids)))

        sum_all_reads = 0
        ab_br = []
        for taxid in cols_sorted:
            [name, all_reads, lvl_reads, added_reads] = lvl_taxids[str(taxid)]
            sum_all_reads += float(added_reads)
            ab_br.append(float(added_reads))

        bracken_ab = [reads / sum_all_reads for reads in ab_br]

    # ---- Build results DataFrame ----
    ab_df = pd.DataFrame({
        'abundance_ls': ab,
        'abundance_br': bracken_ab,
        'taxid': cols_sorted,
        'taxon_removal_cycle': taxon_removal_cycle,
        'GC content': gc_content,
        'abundance_eff': ab_efficiencies
    })
    ab_df.set_index('taxid', inplace=True, drop=False)

    # ---- Check for excessive removal ----
    excessive_removal_detected, efficiencies = check_excessive_removal(ab_df, efficiencies)

    # ---- Write efficiencies ----
    efficiencies = efficiencies / np.max(efficiencies)
    np.savetxt('efficiencies.txt', efficiencies)

    # ---- Write output ----
    write_output(output, lvl_taxids, ab_df, sum_all_reads, level, excessive_removal_detected)

    # ---- Report timing ----
    time_end = strftime("%m-%d-%Y %H:%M:%S", localtime())
    sys.stdout.write("DONE AT: " + time_end + '\n')
    tdelta = (datetime.strptime(time_end, "%m-%d-%Y %H:%M:%S") -
              datetime.strptime(time_start, "%m-%d-%Y %H:%M:%S"))
    print("DURATION: " + str(tdelta.total_seconds()) + " seconds")
    np.savetxt('duration_time.txt', np.array([tdelta.total_seconds()]))


if __name__ == "__main__":
    main()
