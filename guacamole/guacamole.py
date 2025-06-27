# %%
import os, sys, argparse
from time import gmtime, localtime
from time import strftime
from datetime import datetime
import pandas as pd
import numpy as np
import Bio.SeqIO as SeqIO
import guacamole.library as lib
import matplotlib.pyplot as plt
import seaborn as sns
from math import floor
from numpy.random import uniform
from shutil import copyfile
import subprocess
import contextlib
import argparse


np.seterr(divide='ignore', invalid='ignore')

# Tree class
# usage: tree node used in constructing a taxonomy tree
#   including only the taxonomy levels and genomes identified in the Kraken report
class Tree:
    'Tree node.'

    def __init__(self, name, taxid, level_num, level_id, all_reads, lvl_reads, children=None, parent=None):
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

    # process_kmer_distribution


# usage: parses a single line in the kmer distribution file and extracts
# relevant information for the genomes in this sample
# input:
#   - kmer distribution file generaed by generate_kmer_distribution.py.
# returns:
#   - classification taxonomy ID for this line
#   - dictionary of genomes/fractions of the genomes mapping to this classification
def process_kmer_distribution(curr_str, lvl_taxids, map2lvl_taxids):
    split_str = curr_str.strip().split('\t')
    # Parse each genome taxid mapping to m_taxid -- create dictionary instance
    temp_dict = {}
    mapped_taxid = split_str[0]
    for genome_str in split_str[1].split(' '):
        [g_taxid, mkmers, tkmers] = genome_str.split(':')
        mkmers = float(mkmers)
        tkmers = float(tkmers)
        fraction = mkmers / tkmers
        # Only include mappings for genomes within this sample
        if g_taxid in lvl_taxids or g_taxid in map2lvl_taxids:
            if g_taxid not in temp_dict:
                temp_dict[g_taxid] = [fraction]
            else:
                temp_dict[g_taxid].append(fraction)
    # Error check for relevant classifications
    if len(temp_dict) == 0:
        return [-1, {}]
    # Return dictionary
    return [mapped_taxid, temp_dict]


# process_kraken_report
# usage: parses a single line in the kraken report and extracts relevant information
# input: kraken report file with the following tab delimited lines
#   - percent of total reads
#   - number of reads (including at lower levels)
#   - number of reads (only at this level)
#   - taxonomy classification of level
#       (U, - (root), - (cellular org), D, P, C, O, F, G, S)
#   - taxonomy ID (0 = unclassified, 1 = root, 2 = Bacteria...etc)
#   - spaces + name
# returns:
#   - classification/genome name
#   - taxonomy ID for this classification
#   - level for this classification (number)
#   - level name (U, -, D, P, C, O, F, G, S)
#   - all reads classified at this level and below in the tree
#   - reads classified only at this level
def process_kraken_report(curr_str):
    split_str = curr_str.strip().split('\t')
    try:
        int(split_str[1])
    except ValueError:
        return []
    # Extract relevant information
    all_reads = int(split_str[1])
    level_reads = int(split_str[2])
    level_type = split_str[-3]
    taxid = split_str[-2]
    # Get name and spaces
    spaces = 0
    name = split_str[-1]
    for char in name:
        if char == ' ':
            name = name[1:]
            spaces += 1
        else:
            break
            # Determine which level based on number of spaces
    level_num = int(spaces / 2)
    return [name, taxid, level_num, level_type, all_reads, level_reads]


# check_report_file
# usage: checks the format of the report file.
#   - cannot be kraken output file
#   - cannot be mpa style
#   - expect number of columns?
# input: name of input file
# returns: 0 for correct, 1 for incorrect
def check_report_file(in_file):
    sys.stderr.write(">> Checking report file: %s\n" % in_file)
    r_file = open(in_file, 'r')
    first_line = r_file.readline()
    # Test for kraken output file
    if first_line[0] == "C" or first_line[0] == "U":
        sys.stderr.write("\tERROR: Bracken does not use the Kraken default output.\n")
        sys.stderr.write("\t       Bracken requires the Kraken report file (--report option with Kraken)\n")
        exit(1)
        # Test for mpa style
    if len(first_line.split("\t")) == 2:
        sys.stderr.write("\tERROR: Bracken is not compatible with mpa-style reports.\n")
        sys.stderr.write("\t       Bracken requires the default Kraken report format\n")
        exit(1)
    r_file.close()
    return 0


def create_refdist(seqid, kmer_len, index, genome_lib, nbin=100, sampling=50):
    record_dict = SeqIO.index_db(index, genome_lib, 'fasta')
    sequence = record_dict[seqid].seq
    occ = np.zeros(nbin + 1)
    for i in range(len(sequence) - kmer_len + 1)[::sampling]:
        kmer = sequence.upper()[i: i + kmer_len]
        gc = kmer.count("G") + kmer.count("C")
        atgc = kmer.count("A") + kmer.count("T") + gc
        if atgc <= nbin:
            continue
        gc_content = floor((nbin / atgc) * (gc + uniform()))
        occ[gc_content] += 1
    return occ


def get_sample_dists(sdists, weights, map2lvl_taxids, lvl_taxids, gc_dists=None, scale=True, reference=False,
                     kraken_db=None):
    sample_dists = {}
    for genome in map2lvl_taxids:
        genome_dist = 0
        for parent in weights[genome]:
            if parent not in sdists.keys():
                continue
            genome_dist += weights[genome][parent] * sdists[parent]

        lvl_taxid = map2lvl_taxids[genome][0]
        sample_dists[lvl_taxid] = genome_dist + sample_dists[
            lvl_taxid] if lvl_taxid in sample_dists else genome_dist

    nan_txids = [int(txid) for txid in sample_dists if np.sum(sample_dists[txid]) == 0]

    if len(nan_txids) != 0 and reference:
        nodes_dmp = pd.read_csv(os.path.join(kraken_db, 'taxonomy/nodes.dmp'), sep='|', header=None, usecols=[0, 1])
        for txid in nan_txids:
            children = [str(child) for child in np.array(nodes_dmp.loc[nodes_dmp[1] == txid][0])]
            children = [child for child in children if child in gc_dists.columns]
            sample_dists[str(txid)] = np.array(gc_dists[children].mean(1) / np.sum(gc_dists[children].mean(1)))
    if scale:
        for txid in sample_dists:
            sample_dists[txid] = sample_dists[txid] / np.sum(sample_dists[txid])
            if reference:
                continue
            sample_dists[txid] = sample_dists[txid] * lvl_taxids[txid][3]
    return sample_dists


def get_genome_length(kraken_db, nbin, lvl_taxids, map2lvl_taxids, est_reads_dct, read_len, insert=None):
    if insert is None:
        fname = os.path.join(kraken_db, 'gc_bin_' + str(nbin) + '_kmer_' + str(read_len) + '_dist.csv')
    else:
        fname = os.path.join(kraken_db, 'gc_bin_' + str(nbin) + '_kmer_' + str(read_len) +
                             '_insert_' + str(insert) + '_dist.csv')
    gc_dists = pd.read_csv(fname, index_col=0)
    g_seqids = gc_dists.drop('taxid', axis=1).groupby('seqids')
    seq_lengths = g_seqids.sum().sum(axis=1) * 50
    g_taxid = gc_dists.set_index('seqids').sort_index().loc[seq_lengths > 400000, :].groupby('taxid')
    # g_taxid = gc_dists.set_index('seqids').sort_index().loc[seq_lengths > 0, :].groupby('taxid')
    g_taxid = g_taxid.mean().sum(axis=1)
    g_taxid.index = [str(i) for i in g_taxid.index]
    lengths = {txid: 0 for txid in lvl_taxids.keys()}
    weights = {txid: 0 for txid in lvl_taxids.keys()}
    for taxid in [key for key in map2lvl_taxids.keys()]:
        if taxid not in g_taxid.index or taxid not in est_reads_dct.keys():
            continue
        lengths[map2lvl_taxids[taxid][0]] += g_taxid[taxid] * est_reads_dct[taxid]
        weights[map2lvl_taxids[taxid][0]] += est_reads_dct[taxid]
    for taxid in lengths.keys():
        lengths[taxid] = np.sum(lengths[taxid]) / np.sum(weights[taxid])
    sorted_keys = np.sort([int(key) for key in lengths.keys()])
    lengths_out = [lengths[str(taxid)] for taxid in sorted_keys]
    if np.any(np.isnan(lengths_out)):
        nodes_dmp = pd.read_csv(os.path.join(kraken_db, 'taxonomy/nodes.dmp'), sep='|', header=None, usecols=[0, 1])
        nan_taxids = np.where(np.isnan(lengths_out))[0]
        nan_taxids = sorted_keys[nan_taxids]
        for txid in nan_taxids:
            children = nodes_dmp.loc[nodes_dmp[1] == txid][0]
            children_lengths = [g_taxid[str(child)] for child in children if str(child) in g_taxid.index]
            lengths[str(txid)] = np.mean(children_lengths)
    lengths_out = [lengths[str(taxid)] for taxid in sorted_keys]
    if np.any(np.isnan(lengths_out)):
        nan_taxids = np.where(np.isnan(lengths_out))[0]
        nan_taxids = sorted_keys[nan_taxids]
        print("WARNING: COULD NOT ESTIMATE GENOME SIZE FOR TAXIDS: " + str(nan_taxids) + " WILL GIVE THEM MEAN LENGTH")
        for txid in nan_taxids:
            lengths[str(txid)] = np.nanmean(lengths_out)
    lengths_out = [lengths[str(taxid)] for taxid in sorted_keys]
    return lengths_out / np.sum(lengths_out)


def main():
    
    parser = argparse.ArgumentParser(
        description='Run GuaCAMOLE for GC aware species abundance estimation from metagenomic data')

    parser.add_argument('--kraken_report', metavar='kraken_report', type=str, help='.report file from Kraken2 output')
    parser.add_argument('--kraken_file', metavar='kraken_file', type=str, help='.kraken file from Kraken2 output')
    parser.add_argument('--kraken_db', metavar='kraken_db', type=str, help='Path to Kraken2 Database')
    parser.add_argument('--read_len', metavar='read_len', type=int, help='read length')
    parser.add_argument('--output', metavar='output', type=str, help='name of the output file')
    parser.add_argument('--read_files', metavar='read_files', type=str, help='path(s) to read file(s)', nargs='+')
    parser.add_argument('--threshold', metavar='threshold', type=int,
                        help='Minimum number of reads found for a species to estimate its abundance', default=500)
    parser.add_argument('--level', metavar='level', type=str,
                        help='Taxonomy level to be quantified (S=species, G=genus)', default='S')
    parser.add_argument('--length_correction', metavar='length_correction', type=bool,
                        help='genome size correction', default=False)
    parser.add_argument('--plot', metavar='plot', type=bool, help='True if detailed plots should be generated',
                        default=False)
    parser.add_argument('--fp_cycles', metavar='fp_cycles', type=int,
                        help='Number of iterations for false positive removal', default=5)
    parser.add_argument('--reg_weight', metavar='reg_weight', type=float,
                        help='Determines how strong the regularization should be [between 0 and 1]', default=0.01)
    parser.add_argument('--fragment_len', metavar='fragment_len', type=int,
                        help='length of the fragment if paired end and known', default=None)
    parser.add_argument('--fasta', metavar='fasta', type=bool, help='True if reads are in fasta format, false if fastq',
                        default=False)
    parser.add_argument('--quantiles', metavar='quantiles', type=float,
                        help='min and max quantiles of reads that should be used for GC distributions', nargs=2, default=[0.025, 0.975])


    args = parser.parse_args()
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
    seqid2taxid = pd.read_csv(os.path.join(kraken_db, 'seqid2taxid.map'), sep='\t', header=None)
    seqid2taxid.columns = ['seqid', 'taxid']

    time_start = strftime("%m-%d-%Y %H:%M:%S", localtime())
    sys.stdout.write("PROGRAM START TIME: " + time_start + '\n')

    lvl_dict = {'D': 'domains', 'P': 'phylums', 'O': 'orders', 'C': 'classes', 'F': 'families', 'G': 'genuses',
                'S': 'species'}
    if level in lvl_dict:
        abundance_lvl = lvl_dict[level]
    branch = 0
    if len(level) > 1:
        branch = int(level[1:])

    # Initialize variables
    root_node = -1
    prev_node = -1
    main_lvls = ['R', 'K', 'D', 'P', 'C', 'O', 'F', 'G', 'S']
    branch_lvl = main_lvls.index(level[0])
    total_reads = 0
    kept_reads = 0
    ignored_reads = 0
    n_lvl_total = 0
    n_lvl_est = 0
    n_lvl_del = 0
    lvl_nodes = []
    leaf_nodes = []

    all_nodes = {}
    desired_lvl_dct = {}
    est_reads_dct = {}

    # Error Check
    check_report_file(report)

    # Parse kraken report file /and create tree
    i_file = open(report, 'r')
    map2lvl_taxids = {}
    lvl_taxids = {}
    last_taxid = -1
    for line in i_file:
        report_vals = process_kraken_report(line)
        if len(report_vals) < 5:
            continue
        [name, taxid, level_num, level_id, all_reads, level_reads] = report_vals
        total_reads += level_reads
        # Skip unclassified
        if level_id == 'U':
            unclassified_line = line
            u_reads = level_reads
            continue
        # Tree Root
        if taxid == '1':
            root_node = Tree(name, taxid, level_num, 'R', all_reads, level_reads)
            prev_node = root_node
            continue
            # Save leaf nodes
        if level_num != (prev_node.level_num + 1):
            leaf_nodes.append(prev_node)
        # Move to correct parent
        while level_num != (prev_node.level_num + 1):
            prev_node = prev_node.parent
        # Determine correct level ID
        test_branch = 0
        if level_id == '-' or len(level_id) > 1:
            if prev_node.level_id in main_lvls:
                level_id = prev_node.level_id + '1'
                test_branch = 1
            else:
                num = int(prev_node.level_id[-1]) + 1
                test_branch = num
                level_id = prev_node.level_id[:-1] + str(num)
        # Desired level for abundance estimation or below
        if level_id == level:
            n_lvl_total += 1
            # Account for threshold at level
            if all_reads < int(threshold):
                n_lvl_del += 1
                ignored_reads += all_reads
                last_taxid = -1
            else:
                # If level contains enough reads - save for abundance estimation
                n_lvl_est += 1
                kept_reads += all_reads
                lvl_taxids[taxid] = [name, all_reads, level_reads, 0]
                # lvl_taxids[taxid] = [name, all_reads, 0, 0]
                last_taxid = taxid
                # also distribute level reads....
                map2lvl_taxids[taxid] = [taxid, level_reads, 0]
                desired_lvl_dct[taxid] = []
        elif (branch > 0 and test_branch > branch):
            # For all nodes below desired level
            if last_taxid != -1:
                map2lvl_taxids[taxid] = [last_taxid, level_reads, 0]
                desired_lvl_dct[last_taxid].append(taxid)
        elif main_lvls.index(level_id[0]) >= branch_lvl:
            # For all nodes below the desired level
            if last_taxid != -1:
                map2lvl_taxids[taxid] = [last_taxid, level_reads, 0]
                desired_lvl_dct[last_taxid].append(taxid)
        # Add node to tree
        curr_node = Tree(name, taxid, level_num, level_id, all_reads, level_reads, None, prev_node)
        prev_node.add_child(curr_node)
        # Save node
        all_nodes[curr_node.taxid] = curr_node
        prev_node = curr_node
    i_file.close()
    # Add last node
    leaf_nodes.append(prev_node)

    # Read in kmer distribution file
    k_file = open(kmer_distr, 'r')
    kmer_distr_dict = {}
    for line in k_file.readlines()[1:]:
        [mapped_taxid, mapped_taxid_dict] = process_kmer_distribution(line, lvl_taxids, map2lvl_taxids)
        if len(mapped_taxid_dict) == 0:
            continue
        kmer_distr_dict[mapped_taxid] = mapped_taxid_dict
    k_file.close()

    ################################### ADDED BY LAURENZ
    # read in refgen dists
    if fragment_len is not None:
        insert = fragment_len - 2 * read_len
        gc_dists = pd.read_csv(os.path.join(
            kraken_db, 'gc_bin_' + str(nbin) + '_kmer_' + str(read_len) +
                       '_insert_' + str(insert) + '_dist.csv'), index_col=0)
    else:
        insert = None
        gc_dists = pd.read_csv(os.path.join(
            kraken_db, 'gc_bin_' + str(nbin) + '_kmer_' + str(read_len) +
                       '_dist.csv'), index_col=0)

    gc_dists.drop('seqids', axis=1, inplace=True)
    g = gc_dists.groupby('taxid')
    gc_dists = g.sum().T
    gc_dists.columns = [str(i) for i in gc_dists.columns]
    reference_nodes = {}
    distribution_weights = {}
    for taxid in map2lvl_taxids:
        distribution_weights[taxid] = {}
    #################################### END ADDED BY LAURENZ

    # For each PARENT node, distribute level reads to genomes
    curr_nodes = [root_node]
    nondistributed_reads = 0
    distributed_reads = 0
    lvl_reads = 0
    while len(curr_nodes) > 0:
        curr_node = curr_nodes.pop(0)
        # For each child node, add to list of nodes to evaluate
        if not isinstance(curr_node, Tree):
            continue
        # Do not redistribute level reads
        if curr_node.level_id == level:
            pass
            # continue
            # If above level, append
        for child_node in curr_node.children:
            curr_nodes.append(child_node)
            # No reads to distribute
        if curr_node.lvl_reads == 0:
            continue
        # No genomes produce this classification
        if curr_node.taxid not in kmer_distr_dict:
            # print curr_node.name
            nondistributed_reads += curr_node.lvl_reads
            # if curr_node.taxid in distribution_weights:
            #     distribution_weights.pop(curr_node.taxid)
            #     map2lvl_taxids.pop(curr_node.taxid)
            # if curr_node.taxid in lvl_taxids:
            #     lvl_taxids.pop(curr_node.taxid)
            continue
        # Get the dictionary listing all genomes mapping to this node
        distributed_reads += curr_node.lvl_reads
        curr_dict = kmer_distr_dict[curr_node.taxid]
        probability_dict_prelim = {}
        all_genome_reads = 0
        reference_nodes[curr_node.taxid] = 0

        # go through all genomes reads of which map to this node
        for genome in curr_dict:
            # Get the fraction of kmers of the genome expected to map to this node
            fraction = float(curr_dict[genome][0])

            # Determine the number of reads classified by Kraken uniquely for the genome
            # and the fraction of the genome that is unique
            num_classified_reads = float(map2lvl_taxids[genome][1])
            if genome in kmer_distr_dict and genome in kmer_distr_dict[genome]:
                lvl_fraction = float(kmer_distr_dict[genome][genome][0])
            else:
                lvl_fraction = 1.
            # Based on the classified reads and the fraction of unique reads, estimate
            # the true number of reads belonging to this genome in the sample
            est_genome_reads = num_classified_reads / lvl_fraction
            all_genome_reads += est_genome_reads
            # ADDED BY LAURENZ
            est_reads_dct[genome] = est_genome_reads
            # norm_gc_dist = np.array(gc_dists[genome]) / np.sum(np.array(gc_dists[genome]))
            # reference_nodes[curr_node.taxid] += fraction * norm_gc_dist

            # Save values
            probability_dict_prelim[genome] = [fraction, est_genome_reads]

        if all_genome_reads == 0:
            continue

        # Get final probabilities
        # P_R_A = probability that a read is classified at the node given that it belongs to genome A
        # P_A = probability that a randomly selected read belongs to genome A
        # P_A_R = probability that a read belongs to genome A given that its classified at the node
        total_probability = 0.0
        probability_dict_final = {}
        for genome in probability_dict_prelim:
            [P_R_A, est_g_reads] = probability_dict_prelim[genome]
            P_A = float(est_g_reads) / float(all_genome_reads)
            P_A_R = float(P_R_A) * float(P_A)
            probability_dict_final[genome] = P_A_R
            total_probability += P_A_R

        # Find the normalize probabilty and Distribute reads accordingly
        for genome in probability_dict_final:
            add_fraction = probability_dict_final[genome] / total_probability
            add_reads = add_fraction * float(curr_node.lvl_reads)
            map2lvl_taxids[genome][2] += add_reads
            distribution_weights[genome][curr_node.taxid] = add_reads
            norm_gc_dist = np.array(gc_dists[genome]) / np.sum(np.array(gc_dists[genome]))
            reference_nodes[curr_node.taxid] += add_reads * norm_gc_dist

        reference_nodes[curr_node.taxid] = reference_nodes[curr_node.taxid] / np.sum(reference_nodes[curr_node.taxid])

    # For all genomes, map reads up to level
    for genome in map2lvl_taxids:
        [lvl_taxid, all_reads, add_reads] = map2lvl_taxids[genome]
        lvl_taxids[lvl_taxid][3] += add_reads

    if len(lvl_taxids) == 0:
        return None

    # ADDED BY LAURENZ
    refdists = get_sample_dists(sdists=reference_nodes, weights=distribution_weights,
                                map2lvl_taxids=map2lvl_taxids, lvl_taxids=lvl_taxids,
                                reference=True, gc_dists=gc_dists, kraken_db=kraken_db)

    gc_dists_ar = []
    keys = []
    for taxid in refdists.keys():
        keys.append(taxid)
        gc_dists_ar.append([refdists[taxid]][0])

    df = pd.DataFrame(refdists)
    cols_int = [int(col) for col in df.columns]
    cols_sorted = sorted(cols_int)
    cols_str = [str(col) for col in cols_sorted]
    df = df.reindex(cols_str, axis=1)
    refdists = np.array(df)
    ref_normalized = refdists.copy()

    taxids = np.array([int(node) for node in all_nodes])
    # add root node
    taxids = np.append(taxids, 1)
    taxids = np.sort(taxids)
    if not os.path.exists('sample_bin_' + str(nbin) + '.dist'):
        time = strftime("%m-%d-%Y %H:%M:%S", localtime())
        print("Starting to create sample distribution at " + time)
        if fragment_len is not None:
            lib.create_sample_dist_avg(fastq1=fastq1, fastq2=fastq2, txids=taxids, kraken=kraken, fasta=fasta)
        else:
            lib.create_sample_dist(fastq1=fastq1, fastq2=fastq2, txids=taxids, kraken=kraken, fasta=fasta)
        time = strftime("%m-%d-%Y %H:%M:%S", localtime())
        print("Done at " + time)

    sdist = pd.read_csv('sample_bin_' + str(nbin) + '.dist', sep='\t', header=0)
    # remove '#'
    newcols = np.array(sdist.columns)
    newcols[0] = newcols[0][2:]
    col_dct = dict(zip(np.array(sdist.columns), newcols))
    sdist.rename(columns=col_dct, inplace=True)

    sample_dists = get_sample_dists(sdists=sdist, weights=distribution_weights,
                                    map2lvl_taxids=map2lvl_taxids, lvl_taxids=lvl_taxids,
                                    kraken_db=kraken_db)

    df = pd.DataFrame(sample_dists)
    cols_int = [int(col) for col in df.columns]
    cols_sorted = sorted(cols_int)
    cols_str = [str(col) for col in cols_sorted]
    df = df.reindex(cols_str, axis=1)
    s_dists = np.array(df)

    if length_correction:
        lengths = get_genome_length(kraken_db=kraken_db, nbin=nbin, lvl_taxids=lvl_taxids, est_reads_dct=est_reads_dct,
                                    map2lvl_taxids=map2lvl_taxids, insert=insert, read_len=read_len)
        refdists = refdists * lengths

    # sometimes no genome exists even though reads map there, this results in nans in GuaCAMOLE,
    # can still be corrected via efficiencies
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

    refdists = refdists / np.sum(refdists)
    refdists = refdists * np.sum(s_dists)

    np.savetxt('ref_bin_' + str(nbin) + '_input.dist', refdists)
    np.savetxt('sample_bin_' + str(nbin) + '_input.dist', s_dists)

    norm, ref, sample = lib.normalize_gc_dists('sample_bin_' + str(nbin) + '_input.dist',
                                               'ref_bin_' + str(nbin) + '_input.dist',
                                               quantiles=quantiles)
    if plot:
        lib.plot_dist2(norm, line=True)
        plt.savefig('Obs_vs_exp_tresh_' + str(threshold) + '.pdf')
        plt.close()

    ab, taxon_removal_cycle, efficiencies = lib.corrected_abundances('sample_bin_' + str(nbin) + '_input.dist',
                                                                     'ref_bin_' + str(nbin) + '_input.dist',
                                                                     fp_cycles=fp_cycles,
                                                                     taxids=cols_sorted, plot=plot,
                                                                     reg_weight=reg_weight)

    ab[np.isinf(ab)] = np.nan

    # Include Bracken estimates
    # Sum all of the reads for the desired level -- use for fraction of reads
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

    avg_eff = np.matmul(ref_normalized.T, efficiencies.reshape(101, 1)).flatten()
    ab_efficiencies = (bracken_ab / avg_eff) / np.sum(bracken_ab / avg_eff)
    gc_content = np.argmax(ref_normalized, axis=0)

    if len(nan_taxids) != 0:
        print(f"WARNING: For taxid(s) {nan_taxids} reference distributions could not be created, therefore GuaCAMOLE"
              f" estimtates are set to nan for those!")
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


    ab_df = pd.DataFrame({
        'abundance_ls': ab,
        'abundance_br': bracken_ab,
        'taxid': cols_sorted,
        'taxon_removal_cycle': taxon_removal_cycle,
        'GC content': gc_content,
        'abundance_eff': ab_efficiencies
    })
    ab_df.set_index('taxid', inplace=True, drop=False)

    # Identify taxa with zero (or near-zero) GuaCAMOLE abundance
    zero_guacamole_taxa_df = ab_df[ab_df['abundance_ls'] <= 1e-9]

    # Sum their Bracken relative abundances
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
        print(warning_message, file=sys.stderr)  # Print to standard error

    # write efficiencies.txt
    efficiencies = efficiencies / np.max(efficiencies)
    np.savetxt('efficiencies.txt', efficiencies)

    ## write output
    o_file = open(output, 'w')
    o_file.write('name\t' + 'taxonomy_id\t' + 'taxonomy_lvl\t' + 'kraken_assigned_reads\t' + 'added_reads\t' +
                 'new_est_reads\t' + 'fraction_total_reads\t' + 'Bracken_estimate\t' + 'GuaCAMOLE_estimate\t' +
                 'GuaCAMOLE_est_eff\t' + 'GC content\t' + 'Taxon Removal Cycle\n')

    for taxid in lvl_taxids:
        [name, all_reads, lvl_reads, added_reads] = lvl_taxids[taxid]
        # Count up all added reads + all_reads already at the level
        new_all_reads = float(all_reads) + float(added_reads)

        bracken_est_val = float(ab_df.loc[ab_df['taxid'] == int(taxid), 'abundance_br'].iloc[0])
        taxon_removal_cycle_val = str(ab_df.loc[ab_df['taxid'] == int(taxid), 'taxon_removal_cycle'].iloc[0])
        if excessive_removal_detected:
            guacamole_ls_est_val = np.nan
            guacamole_eff_est_val = np.nan
        else:
            # Get GuaCAMOLE estimates (original logic)
            guacamole_ls_est_val = float(ab_df.loc[ab_df['taxid'] == int(taxid), 'abundance_ls'].iloc[0])
            guacamole_eff_est_val = float(ab_df.loc[ab_df['taxid'] == int(taxid), 'abundance_eff'].iloc[0])

        gc_content_val = float(ab_df.loc[ab_df['taxid'] == int(taxid), 'GC content'].iloc[0])

        # Output
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

    time_end = strftime("%m-%d-%Y %H:%M:%S", localtime())
    sys.stdout.write("DONE AT: " + time_end + '\n')
    tdelta = datetime.strptime(time_end, "%m-%d-%Y %H:%M:%S") - datetime.strptime(time_start, "%m-%d-%Y %H:%M:%S")
    print("DURATION: " + str(tdelta.total_seconds()) + " seconds")
    np.savetxt('duration_time.txt', np.array([tdelta.total_seconds()]))


if __name__ == "__main__":
    main()
