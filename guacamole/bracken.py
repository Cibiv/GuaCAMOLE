"""Taxonomy tree construction, Kraken report parsing, and read redistribution.

The code in this file is adapted from Bracken (Lu et al., 2017,
PeerJ Computer Science 3:e104), available at
https://github.com/jenniferlu717/Bracken under the GPL-3.0 license.
"""

import sys

import numpy as np


class Tree:
    """Tree node used in constructing a taxonomy tree.

    Includes only the taxonomy levels and genomes identified in the Kraken report.
    """

    def __init__(self, name, taxid, level_num, level_id, all_reads, lvl_reads,
                 children=None, parent=None):
        self.name = name
        self.taxid = taxid
        self.level_num = level_num
        self.level_id = level_id
        self.all_reads = all_reads
        self.lvl_reads = lvl_reads
        self.children = []
        self.parent = parent
        if children is not None:
            for child in children:
                self.add_child(child)

    def add_child(self, node):
        assert isinstance(node, Tree)
        self.children.append(node)


def process_kmer_distribution(curr_str, lvl_taxids, map2lvl_taxids):
    """Parse a single line in the kmer distribution file.

    Extracts relevant information for the genomes in this sample.

    Parameters
    ----------
    curr_str : str
        A line from the kmer distribution file.
    lvl_taxids : dict
        Dictionary of taxids at the desired taxonomic level.
    map2lvl_taxids : dict
        Dictionary mapping genome taxids to their level taxid.

    Returns
    -------
    list
        [mapped_taxid, dict] where dict maps genome taxids to lists of fractions.
    """
    split_str = curr_str.strip().split('\t')
    temp_dict = {}
    mapped_taxid = split_str[0]
    for genome_str in split_str[1].split(' '):
        [g_taxid, mkmers, tkmers] = genome_str.split(':')
        mkmers = float(mkmers)
        tkmers = float(tkmers)
        fraction = mkmers / tkmers
        if g_taxid in lvl_taxids or g_taxid in map2lvl_taxids:
            if g_taxid not in temp_dict:
                temp_dict[g_taxid] = [fraction]
            else:
                temp_dict[g_taxid].append(fraction)
    if len(temp_dict) == 0:
        return [-1, {}]
    return [mapped_taxid, temp_dict]


def process_kraken_report(curr_str):
    """Parse a single line in the kraken report.

    Parameters
    ----------
    curr_str : str
        A line from the kraken report file.

    Returns
    -------
    list
        [name, taxid, level_num, level_type, all_reads, level_reads] or empty list.
    """
    split_str = curr_str.strip().split('\t')
    try:
        int(split_str[1])
    except ValueError:
        return []
    all_reads = int(split_str[1])
    level_reads = int(split_str[2])
    level_type = split_str[-3]
    taxid = split_str[-2]
    spaces = 0
    name = split_str[-1]
    for char in name:
        if char == ' ':
            name = name[1:]
            spaces += 1
        else:
            break
    level_num = int(spaces / 2)
    return [name, taxid, level_num, level_type, all_reads, level_reads]


def check_report_file(in_file):
    """Check the format of the kraken report file.

    Parameters
    ----------
    in_file : str
        Path to the report file.

    Returns
    -------
    int
        0 if correct format.
    """
    sys.stderr.write(">> Checking report file: %s\n" % in_file)
    r_file = open(in_file, 'r')
    first_line = r_file.readline()
    if first_line[0] == "C" or first_line[0] == "U":
        sys.stderr.write("\tERROR: Bracken does not use the Kraken default output.\n")
        sys.stderr.write("\t       Bracken requires the Kraken report file (--report option with Kraken)\n")
        exit(1)
    if len(first_line.split("\t")) == 2:
        sys.stderr.write("\tERROR: Bracken is not compatible with mpa-style reports.\n")
        sys.stderr.write("\t       Bracken requires the default Kraken report format\n")
        exit(1)
    r_file.close()
    return 0


def parse_kraken_report(report, level, threshold):
    """Parse the full kraken report and build the taxonomy tree.

    Parameters
    ----------
    report : str
        Path to the kraken report file.
    level : str
        Taxonomy level to quantify (e.g., 'S' for species).
    threshold : int
        Minimum number of reads for a species to be included.

    Returns
    -------
    tuple
        (root_node, all_nodes, lvl_taxids, map2lvl_taxids, desired_lvl_dct,
         est_reads_dct, total_reads, kept_reads, ignored_reads,
         n_lvl_total, n_lvl_est, n_lvl_del)
    """
    main_lvls = ['R', 'K', 'D', 'P', 'C', 'O', 'F', 'G', 'S']
    branch = 0
    if len(level) > 1:
        branch = int(level[1:])
    branch_lvl = main_lvls.index(level[0])

    root_node = None
    prev_node = None
    total_reads = 0
    kept_reads = 0
    ignored_reads = 0
    n_lvl_total = 0
    n_lvl_est = 0
    n_lvl_del = 0
    leaf_nodes = []

    all_nodes = {}
    lvl_taxids = {}
    map2lvl_taxids = {}
    desired_lvl_dct = {}
    est_reads_dct = {}
    last_taxid = -1

    check_report_file(report)

    i_file = open(report, 'r')
    for line in i_file:
        report_vals = process_kraken_report(line)
        if len(report_vals) < 5:
            continue
        [name, taxid, level_num, level_id, all_reads, level_reads] = report_vals
        total_reads += level_reads
        if level_id == 'U':
            continue
        if taxid == '1':
            root_node = Tree(name, taxid, level_num, 'R', all_reads, level_reads)
            prev_node = root_node
            continue
        if level_num != (prev_node.level_num + 1):
            leaf_nodes.append(prev_node)
        while level_num != (prev_node.level_num + 1):
            prev_node = prev_node.parent
        test_branch = 0
        if level_id == '-' or len(level_id) > 1:
            if prev_node.level_id in main_lvls:
                level_id = prev_node.level_id + '1'
                test_branch = 1
            else:
                num = int(prev_node.level_id[-1]) + 1
                test_branch = num
                level_id = prev_node.level_id[:-1] + str(num)
        if level_id == level:
            n_lvl_total += 1
            if all_reads < int(threshold):
                n_lvl_del += 1
                ignored_reads += all_reads
                last_taxid = -1
            else:
                n_lvl_est += 1
                kept_reads += all_reads
                lvl_taxids[taxid] = [name, all_reads, level_reads, 0]
                last_taxid = taxid
                map2lvl_taxids[taxid] = [taxid, level_reads, 0]
                desired_lvl_dct[taxid] = []
        elif (branch > 0 and test_branch > branch):
            if last_taxid != -1:
                map2lvl_taxids[taxid] = [last_taxid, level_reads, 0]
                desired_lvl_dct[last_taxid].append(taxid)
        elif main_lvls.index(level_id[0]) >= branch_lvl:
            if last_taxid != -1:
                map2lvl_taxids[taxid] = [last_taxid, level_reads, 0]
                desired_lvl_dct[last_taxid].append(taxid)
        curr_node = Tree(name, taxid, level_num, level_id, all_reads, level_reads, None, prev_node)
        prev_node.add_child(curr_node)
        all_nodes[curr_node.taxid] = curr_node
        prev_node = curr_node
    i_file.close()
    leaf_nodes.append(prev_node)

    return (root_node, all_nodes, lvl_taxids, map2lvl_taxids, desired_lvl_dct,
            est_reads_dct, total_reads, kept_reads, ignored_reads,
            n_lvl_total, n_lvl_est, n_lvl_del)


def parse_kmer_distributions(kmer_distr, lvl_taxids, map2lvl_taxids):
    """Read and parse the kmer distribution file.

    Parameters
    ----------
    kmer_distr : str
        Path to the kmer distribution file.
    lvl_taxids : dict
        Dictionary of taxids at the desired taxonomic level.
    map2lvl_taxids : dict
        Dictionary mapping genome taxids to their level taxid.

    Returns
    -------
    dict
        Dictionary mapping mapped taxids to their genome fraction dictionaries.
    """
    k_file = open(kmer_distr, 'r')
    kmer_distr_dict = {}
    for line in k_file.readlines()[1:]:
        [mapped_taxid, mapped_taxid_dict] = process_kmer_distribution(line, lvl_taxids, map2lvl_taxids)
        if len(mapped_taxid_dict) == 0:
            continue
        kmer_distr_dict[mapped_taxid] = mapped_taxid_dict
    k_file.close()
    return kmer_distr_dict


def redistribute_reads(root_node, level, kmer_distr_dict, map2lvl_taxids, lvl_taxids,
                       gc_dists, distribution_weights):
    """Redistribute reads from higher taxonomic levels to genomes.

    This implements the Bracken-style probabilistic redistribution,
    while also tracking GC distribution weights and reference node distributions.

    Parameters
    ----------
    root_node : Tree
        Root of the taxonomy tree.
    level : str
        Taxonomy level (e.g., 'S').
    kmer_distr_dict : dict
        Parsed kmer distribution data.
    map2lvl_taxids : dict
        Mapping from genome taxids to level taxids (modified in place).
    lvl_taxids : dict
        Level taxid information (modified in place).
    gc_dists : pd.DataFrame
        GC content distributions per genome.
    distribution_weights : dict
        Weights for distribution (modified in place).

    Returns
    -------
    tuple
        (reference_nodes, est_reads_dct)
    """
    reference_nodes = {}
    est_reads_dct = {}

    curr_nodes = [root_node]
    nondistributed_reads = 0
    distributed_reads = 0

    while len(curr_nodes) > 0:
        curr_node = curr_nodes.pop(0)
        if not isinstance(curr_node, Tree):
            continue
        if curr_node.level_id == level:
            pass
        for child_node in curr_node.children:
            curr_nodes.append(child_node)
        if curr_node.lvl_reads == 0:
            continue
        if curr_node.taxid not in kmer_distr_dict:
            nondistributed_reads += curr_node.lvl_reads
            continue

        distributed_reads += curr_node.lvl_reads
        curr_dict = kmer_distr_dict[curr_node.taxid]
        probability_dict_prelim = {}
        all_genome_reads = 0
        reference_nodes[curr_node.taxid] = 0

        for genome in curr_dict:
            fraction = float(curr_dict[genome][0])
            num_classified_reads = float(map2lvl_taxids[genome][1])
            if genome in kmer_distr_dict and genome in kmer_distr_dict[genome]:
                lvl_fraction = float(kmer_distr_dict[genome][genome][0])
            else:
                lvl_fraction = 1.
            est_genome_reads = num_classified_reads / lvl_fraction
            all_genome_reads += est_genome_reads
            est_reads_dct[genome] = est_genome_reads
            probability_dict_prelim[genome] = [fraction, est_genome_reads]

        if all_genome_reads == 0:
            continue

        total_probability = 0.0
        probability_dict_final = {}
        for genome in probability_dict_prelim:
            [P_R_A, est_g_reads] = probability_dict_prelim[genome]
            P_A = float(est_g_reads) / float(all_genome_reads)
            P_A_R = float(P_R_A) * float(P_A)
            probability_dict_final[genome] = P_A_R
            total_probability += P_A_R

        for genome in probability_dict_final:
            add_fraction = probability_dict_final[genome] / total_probability
            add_reads = add_fraction * float(curr_node.lvl_reads)
            map2lvl_taxids[genome][2] += add_reads
            distribution_weights[genome][curr_node.taxid] = add_reads
            norm_gc_dist = np.array(gc_dists[genome]) / np.sum(np.array(gc_dists[genome]))
            reference_nodes[curr_node.taxid] += add_reads * norm_gc_dist

        reference_nodes[curr_node.taxid] = (
            reference_nodes[curr_node.taxid] / np.sum(reference_nodes[curr_node.taxid])
        )

    # Map reads up to level for all genomes
    for genome in map2lvl_taxids:
        [lvl_taxid, all_reads, add_reads] = map2lvl_taxids[genome]
        lvl_taxids[lvl_taxid][3] += add_reads

    return reference_nodes, est_reads_dct
