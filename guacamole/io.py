import sys

import numpy as np


def write_output(output, lvl_taxids, ab_df, sum_all_reads, level, excessive_removal_detected):
    """Write the GuaCAMOLE output file.

    Parameters
    ----------
    output : str
        Path to the output file.
    lvl_taxids : dict
        Level taxid information with [name, all_reads, lvl_reads, added_reads].
    ab_df : pd.DataFrame
        DataFrame with abundance estimates and metadata.
    sum_all_reads : float
        Total number of reads across all level taxa.
    level : str
        Taxonomy level (e.g., 'S').
    excessive_removal_detected : bool
        Whether GuaCAMOLE flagged excessive taxa removal.
    """
    o_file = open(output, 'w')
    o_file.write(
        'name\t' + 'taxonomy_id\t' + 'taxonomy_lvl\t' + 'kraken_assigned_reads\t' +
        'added_reads\t' + 'new_est_reads\t' + 'fraction_total_reads\t' +
        'Bracken_estimate\t' + 'GuaCAMOLE_estimate\t' + 'GuaCAMOLE_est_eff\t' +
        'GC content\t' + 'Taxon Removal Cycle\n'
    )

    for taxid in lvl_taxids:
        [name, all_reads, lvl_reads, added_reads] = lvl_taxids[taxid]
        new_all_reads = float(all_reads) + float(added_reads)

        bracken_est_val = float(ab_df.loc[ab_df['taxid'] == int(taxid), 'abundance_br'].iloc[0])
        taxon_removal_cycle_val = str(
            ab_df.loc[ab_df['taxid'] == int(taxid), 'taxon_removal_cycle'].iloc[0]
        )

        if excessive_removal_detected:
            guacamole_ls_est_val = np.nan
            guacamole_eff_est_val = np.nan
        else:
            guacamole_ls_est_val = float(
                ab_df.loc[ab_df['taxid'] == int(taxid), 'abundance_ls'].iloc[0]
            )
            guacamole_eff_est_val = float(
                ab_df.loc[ab_df['taxid'] == int(taxid), 'abundance_eff'].iloc[0]
            )

        gc_content_val = float(ab_df.loc[ab_df['taxid'] == int(taxid), 'GC content'].iloc[0])

        o_file.write(name + '\t')
        o_file.write(taxid + '\t')
        o_file.write(level + '\t')
        o_file.write(str(int(all_reads)) + '\t')
        o_file.write(str(int(new_all_reads) - int(all_reads)) + '\t')
        o_file.write(str(int(new_all_reads)) + '\t')
        o_file.write("%0.5f\t" % (float(int(new_all_reads)) / float(int(sum_all_reads))))
        o_file.write("%0.5f\t" % bracken_est_val)
        o_file.write("%0.5f\t" % guacamole_ls_est_val)
        o_file.write("%0.5f\t" % guacamole_eff_est_val)
        o_file.write("%2.0f\t" % gc_content_val)
        o_file.write(str(taxon_removal_cycle_val) + '\n')
    o_file.close()


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
        (excessive_removal_detected: bool, efficiencies: np.ndarray)
    """
    zero_guacamole_taxa_df = ab_df[ab_df['abundance_ls'] <= 1e-9]
    removed_fraction_bracken = zero_guacamole_taxa_df['abundance_br'].sum()

    excessive_removal_detected = False
    if removed_fraction_bracken > 0.50:
        excessive_removal_detected = True
        efficiencies.fill(np.nan)
        warning_message = (
            f"WARNING: GuaCAMOLE assigned zero abundance to taxa representing "
            f"{removed_fraction_bracken * 100:.2f}% of the Bracken-estimated abundance. "
            "This likely indicates a mismatch between sample and reference genomes or other issues. "
            "GuaCAMOLE abundance and efficiency estimates may be unreliable and are reported as NaN."
        )
        print(warning_message, file=sys.stderr)

    return excessive_removal_detected, efficiencies
