import sys

import numpy as np


def write_output(output, lvl_taxids, ab_df, sum_all_reads, level, excessive_removal_detected):
    """Write the GuaCAMOLE output file.

    Parameters
    ----------
    output : str
        Path to the output file.
    lvl_taxids : dict
        Level taxid information with ``[name, all_reads, lvl_reads, est_reads]``,
        where ``est_reads`` is the total post-redistribution read count for the taxon.
    ab_df : pd.DataFrame
        DataFrame with abundance estimates and metadata, indexed by integer ``taxid``.
    sum_all_reads : float
        Total number of redistributed reads across all level taxa.
    level : str
        Taxonomy level (e.g., 'S').
    excessive_removal_detected : bool
        Whether GuaCAMOLE flagged excessive taxa removal.
    """
    with open(output, 'w') as o_file:
        o_file.write(
            'name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\t'
            'added_reads\tnew_est_reads\tfraction_total_reads\t'
            'Bracken_estimate\tGuaCAMOLE_estimate\tGuaCAMOLE_est_eff\t'
            'GC content\tTaxon Removal Cycle\n'
        )

        for taxid, (name, all_reads, _lvl_reads, est_reads) in lvl_taxids.items():
            new_est_reads = float(est_reads)
            new_est_int = round(new_est_reads)
            added_reads = new_est_int - int(all_reads)
            frac_total = (new_est_reads / float(sum_all_reads)) if sum_all_reads > 0 else 0.0

            row = ab_df.loc[int(taxid)]
            bracken_est_val = float(row['abundance_br'])
            taxon_removal_cycle_val = str(row['taxon_removal_cycle'])

            if excessive_removal_detected:
                guacamole_ls_est_val = np.nan
                guacamole_eff_est_val = np.nan
            else:
                guacamole_ls_est_val = float(row['abundance_ls'])
                guacamole_eff_est_val = float(row['abundance_eff'])

            gc_content_val = float(row['GC content'])

            o_file.write(
                f"{name}\t{taxid}\t{level}\t{int(all_reads)}\t{added_reads}\t"
                f"{new_est_int}\t{frac_total:0.5f}\t"
                f"{bracken_est_val:0.5f}\t{guacamole_ls_est_val:0.5f}\t"
                f"{guacamole_eff_est_val:0.5f}\t{gc_content_val:2.0f}\t"
                f"{taxon_removal_cycle_val}\n"
            )


def check_excessive_removal(ab_df, efficiencies):
    """Check if GuaCAMOLE removed an excessive fraction of taxa.

    Parameters
    ----------
    ab_df : pd.DataFrame
        DataFrame with abundance estimates.
    efficiencies : np.ndarray
        Estimated sequencing efficiencies.

    Returns
    -------
    tuple
        ``(excessive_removal_detected: bool, efficiencies: np.ndarray)``
    """
    zero_guacamole_taxa_df = ab_df[ab_df['abundance_ls'] <= 1e-9]
    removed_fraction_bracken = zero_guacamole_taxa_df['abundance_br'].sum()

    excessive_removal_detected = False
    if removed_fraction_bracken > 0.50:
        excessive_removal_detected = True
        efficiencies = np.full_like(efficiencies, np.nan, dtype=float)
        warning_message = (
            f"WARNING: GuaCAMOLE assigned zero abundance to taxa representing "
            f"{removed_fraction_bracken * 100:.2f}% of the Bracken-estimated abundance. "
            "This likely indicates a mismatch between sample and reference genomes or other issues. "
            "GuaCAMOLE abundance and efficiency estimates may be unreliable and are reported as NaN."
        )
        print(warning_message, file=sys.stderr)

    return excessive_removal_detected, efficiencies
