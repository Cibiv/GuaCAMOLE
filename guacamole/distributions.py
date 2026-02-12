import os

import numpy as np
import pandas as pd
import Bio.SeqIO as SeqIO
from math import floor
from numpy.random import uniform


def create_refdist(seqid, kmer_len, index, genome_lib, nbin=100, sampling=50):
    """Create a reference GC distribution for a single sequence.

    Parameters
    ----------
    seqid : str
        Sequence identifier.
    kmer_len : int
        Length of kmers (read length).
    index : str
        Path to the SeqIO index database.
    genome_lib : str
        Path to the genome library FASTA.
    nbin : int
        Number of GC bins.
    sampling : int
        Sampling interval along the sequence.

    Returns
    -------
    np.ndarray
        GC content distribution array.
    """
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


def get_sample_dists(sdists, weights, map2lvl_taxids, lvl_taxids, gc_dists=None,
                     scale=True, reference=False, kraken_db=None):
    """Compute sample or reference distributions aggregated to the desired taxonomic level.

    Parameters
    ----------
    sdists : dict or pd.DataFrame
        Per-node sample distributions (or reference node distributions).
    weights : dict
        Distribution weights mapping genomes to node contributions.
    map2lvl_taxids : dict
        Mapping from genome taxids to level taxids.
    lvl_taxids : dict
        Level taxid information.
    gc_dists : pd.DataFrame, optional
        GC distributions for fallback when no reads are available.
    scale : bool
        Whether to scale distributions.
    reference : bool
        Whether computing reference (True) or sample (False) distributions.
    kraken_db : str, optional
        Path to Kraken2 database (needed for fallback lookups).

    Returns
    -------
    dict
        Dictionary mapping level taxids to their aggregated distributions.
    """
    sample_dists = {}
    for genome in map2lvl_taxids:
        genome_dist = 0
        for parent in weights[genome]:
            if parent not in sdists.keys():
                continue
            genome_dist += weights[genome][parent] * sdists[parent]

        lvl_taxid = map2lvl_taxids[genome][0]
        sample_dists[lvl_taxid] = (
            genome_dist + sample_dists[lvl_taxid] if lvl_taxid in sample_dists else genome_dist
        )

    nan_txids = [int(txid) for txid in sample_dists if np.sum(sample_dists[txid]) == 0]

    if len(nan_txids) != 0 and reference:
        nodes_dmp = pd.read_csv(
            os.path.join(kraken_db, 'taxonomy/nodes.dmp'), sep='|', header=None, usecols=[0, 1]
        )
        for txid in nan_txids:
            children = [str(child) for child in np.array(nodes_dmp.loc[nodes_dmp[1] == txid][0])]
            children = [child for child in children if child in gc_dists.columns]
            sample_dists[str(txid)] = np.array(
                gc_dists[children].mean(1) / np.sum(gc_dists[children].mean(1))
            )

    if scale:
        for txid in sample_dists:
            sample_dists[txid] = sample_dists[txid] / np.sum(sample_dists[txid])
            if reference:
                continue
            sample_dists[txid] = sample_dists[txid] * lvl_taxids[txid][3]

    return sample_dists


def get_genome_length(kraken_db, nbin, lvl_taxids, map2lvl_taxids, est_reads_dct,
                      read_len, insert=None):
    """Estimate genome lengths for each taxon.

    Parameters
    ----------
    kraken_db : str
        Path to the Kraken2 database.
    nbin : int
        Number of GC bins.
    lvl_taxids : dict
        Level taxid information.
    map2lvl_taxids : dict
        Mapping from genome taxids to level taxids.
    est_reads_dct : dict
        Estimated read counts per genome.
    read_len : int
        Read length.
    insert : int, optional
        Insert size (fragment_len - 2 * read_len).

    Returns
    -------
    np.ndarray
        Normalized genome length array, ordered by sorted taxid.
    """
    if insert is None:
        fname = os.path.join(
            kraken_db, 'gc_bin_' + str(nbin) + '_kmer_' + str(read_len) + '_dist.csv'
        )
    else:
        fname = os.path.join(
            kraken_db, 'gc_bin_' + str(nbin) + '_kmer_' + str(read_len) +
            '_insert_' + str(insert) + '_dist.csv'
        )
    gc_dists = pd.read_csv(fname, index_col=0)
    g_seqids = gc_dists.drop('taxid', axis=1).groupby('seqids')
    seq_lengths = g_seqids.sum().sum(axis=1) * 50
    g_taxid = gc_dists.set_index('seqids').sort_index().loc[seq_lengths > 400000, :].groupby('taxid')
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
        nodes_dmp = pd.read_csv(
            os.path.join(kraken_db, 'taxonomy/nodes.dmp'), sep='|', header=None, usecols=[0, 1]
        )
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
        print("WARNING: COULD NOT ESTIMATE GENOME SIZE FOR TAXIDS: " +
              str(nan_taxids) + " WILL GIVE THEM MEAN LENGTH")
        for txid in nan_taxids:
            lengths[str(txid)] = np.nanmean(lengths_out)
    lengths_out = [lengths[str(taxid)] for taxid in sorted_keys]
    return lengths_out / np.sum(lengths_out)


def load_gc_distributions(kraken_db, nbin, read_len, fragment_len=None):
    """Load GC content distributions from the Kraken2 database.

    Parameters
    ----------
    kraken_db : str
        Path to the Kraken2 database.
    nbin : int
        Number of GC bins.
    read_len : int
        Read length.
    fragment_len : int, optional
        Fragment length for paired-end data.

    Returns
    -------
    tuple
        (gc_dists DataFrame with taxid columns, insert size or None)
    """
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
    return gc_dists, insert
