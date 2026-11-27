import os

import numpy as np
import pandas as pd
from Bio import SeqIO


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
        GC content distribution array of length ``nbin + 1``.
    """
    record_dict = SeqIO.index_db(index, genome_lib, 'fasta')
    seq_bytes = str(record_dict[seqid].seq).upper().encode('ascii', errors='ignore')
    seq_len = len(seq_bytes)
    if seq_len < kmer_len:
        return np.zeros(nbin + 1, dtype=float)

    arr = np.frombuffer(seq_bytes, dtype=np.uint8)
    is_gc = (arr == ord('G')) | (arr == ord('C'))
    is_valid = is_gc | (arr == ord('A')) | (arr == ord('T'))

    gc_cumsum = np.empty(seq_len + 1, dtype=np.int64)
    gc_cumsum[0] = 0
    np.cumsum(is_gc, out=gc_cumsum[1:])

    valid_cumsum = np.empty(seq_len + 1, dtype=np.int64)
    valid_cumsum[0] = 0
    np.cumsum(is_valid, out=valid_cumsum[1:])

    starts = np.arange(0, seq_len - kmer_len + 1, sampling)
    gc = gc_cumsum[starts + kmer_len] - gc_cumsum[starts]
    atgc = valid_cumsum[starts + kmer_len] - valid_cumsum[starts]

    keep = atgc > nbin
    if not np.any(keep):
        return np.zeros(nbin + 1, dtype=float)

    gc_k = gc[keep]
    atgc_k = atgc[keep]
    bins = np.floor((nbin / atgc_k) * (gc_k + np.random.random(len(gc_k)))).astype(np.int64)
    bins = bins[(bins >= 0) & (bins <= nbin)]
    return np.bincount(bins, minlength=nbin + 1).astype(float)


def get_sample_dists(sdists, weights, map2lvl_taxids, lvl_taxids, gc_dists=None,
                     scale=True, reference=False, kraken_db=None):
    """Compute sample or reference distributions aggregated to the desired taxonomic level.

    Implements Eqs. (2)-(3) for observed read counts (``reference=False``) and
    Eqs. (5)-(6) for genomic GC content distributions (``reference=True``).

    Parameters
    ----------
    sdists : dict or pd.DataFrame
        Per-node sample count distributions (``reference=False``) or
        normalized reference node distributions ``h(i, g)`` (``reference=True``).
    weights : dict
        Distribution weights ``add_reads = M(i) * P(r in G_j | K_i)`` mapping
        genomes to ancestor node contributions.
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
    # Precompute per-node normalized distributions.
    # For reference=False, sdists[parent] contains raw read counts M(parent, g)
    # whose sum is M(parent) (or 2*M(parent) for uncombined paired-end reads).
    # Since weights[genome][parent] = M(parent) * P(r in G_genome | K_parent),
    # normalizing sdists[parent] by its sum yields:
    #   weights[genome][parent] * (M(parent, g) / sum_g M(parent, g))
    #   = P(r in G_genome | K_parent) * M(parent, g)   [Eq. (2)].
    if isinstance(sdists, pd.DataFrame):
        available_keys = set(sdists.columns)
        n_bins = sdists.shape[0]
    else:
        available_keys = set(sdists.keys())
        n_bins = next((len(v) for v in sdists.values() if hasattr(v, '__len__')), 101)

    norm_sdists = {}
    for key in available_keys:
        vec = np.asarray(sdists[key], dtype=float)
        if vec.ndim == 0:
            continue
        vec_sum = float(np.sum(vec))
        if vec_sum > 0:
            norm_sdists[key] = vec if reference else (vec / vec_sum)

    sample_dists = {}
    for genome in map2lvl_taxids:
        genome_dist = np.zeros(n_bins, dtype=float)
        for parent, w in weights[genome].items():
            parent_dist = norm_sdists.get(parent)
            if parent_dist is None:
                continue
            genome_dist += w * parent_dist

        lvl_taxid = map2lvl_taxids[genome][0]
        if lvl_taxid in sample_dists:
            sample_dists[lvl_taxid] = sample_dists[lvl_taxid] + genome_dist
        else:
            sample_dists[lvl_taxid] = genome_dist

    nan_txids = [int(txid) for txid, dist in sample_dists.items() if np.sum(dist) == 0]

    if len(nan_txids) != 0 and reference and kraken_db is not None and gc_dists is not None:
        nodes_path = os.path.join(kraken_db, 'taxonomy/nodes.dmp')
        if os.path.exists(nodes_path):
            nodes_dmp = pd.read_csv(nodes_path, sep='|', header=None, usecols=[0, 1])
            for txid in nan_txids:
                children = [str(child) for child in np.asarray(nodes_dmp.loc[nodes_dmp[1] == txid, 0])]
                children = [child for child in children if child in gc_dists.columns]
                if children:
                    mean_dist = np.asarray(gc_dists[children].mean(axis=1), dtype=float)
                    mean_sum = float(np.sum(mean_dist))
                    if mean_sum > 0:
                        sample_dists[str(txid)] = mean_dist / mean_sum

    if scale:
        for txid, dist in sample_dists.items():
            dist_sum = float(np.sum(dist))
            if dist_sum <= 0:
                sample_dists[txid] = np.full(n_bins, np.nan, dtype=float)
                continue
            scaled = dist / dist_sum
            if not reference:
                scaled = scaled * lvl_taxids[txid][3]
            sample_dists[txid] = scaled

    return sample_dists


def get_genome_length(kraken_db, nbin, lvl_taxids, map2lvl_taxids, est_reads_dct,
                      read_len, insert=None, gc_dists_raw=None):
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
        Estimated read counts per genome (proportional to Bracken priors P(K_j)).
    read_len : int
        Read length.
    insert : int, optional
        Insert size (fragment_len - 2 * read_len).
    gc_dists_raw : pd.DataFrame, optional
        Pre-loaded GC distributions DataFrame (with ``taxid`` and ``seqids`` columns)
        to avoid re-reading the CSV from disk.

    Returns
    -------
    np.ndarray
        Normalized genome length array, ordered by sorted taxid.
    """
    if gc_dists_raw is None:
        if insert is None:
            fname = os.path.join(
                kraken_db, f'gc_bin_{nbin}_kmer_{read_len}_dist.csv'
            )
        else:
            fname = os.path.join(
                kraken_db, f'gc_bin_{nbin}_kmer_{read_len}_insert_{insert}_dist.csv'
            )
        gc_dists_raw = pd.read_csv(fname, index_col=0)

    # Sum across GC bin columns once into a 1D vector instead of grouping 101 columns
    bin_cols = [c for c in gc_dists_raw.columns if c not in ('taxid', 'seqids')]
    row_sums = gc_dists_raw[bin_cols].sum(axis=1)
    seq_df = pd.DataFrame({
        'seqids': gc_dists_raw['seqids'].values,
        'taxid': gc_dists_raw['taxid'].values,
        'windows': row_sums.values,
    })
    seq_totals = seq_df.groupby('seqids')['windows'].sum() * 50
    valid_seqids = seq_totals.index[seq_totals > 400000]
    g_taxid = seq_df.loc[seq_df['seqids'].isin(valid_seqids)].groupby('taxid')['windows'].mean()
    g_taxid.index = [str(i) for i in g_taxid.index]

    lengths = {txid: 0.0 for txid in lvl_taxids}
    weights = {txid: 0.0 for txid in lvl_taxids}
    for taxid in map2lvl_taxids:
        if taxid not in g_taxid.index or taxid not in est_reads_dct:
            continue
        lvl_txid = map2lvl_taxids[taxid][0]
        lengths[lvl_txid] += float(g_taxid[taxid]) * float(est_reads_dct[taxid])
        weights[lvl_txid] += float(est_reads_dct[taxid])

    for taxid, total_len in lengths.items():
        lengths[taxid] = (total_len / weights[taxid]) if weights[taxid] > 0 else np.nan

    sorted_keys = np.sort([int(key) for key in lengths])
    lengths_out = np.array([lengths[str(taxid)] for taxid in sorted_keys], dtype=float)

    if np.any(np.isnan(lengths_out)):
        nodes_path = os.path.join(kraken_db, 'taxonomy/nodes.dmp')
        if os.path.exists(nodes_path):
            nodes_dmp = pd.read_csv(nodes_path, sep='|', header=None, usecols=[0, 1])
            nan_taxids = sorted_keys[np.isnan(lengths_out)]
            for txid in nan_taxids:
                children = nodes_dmp.loc[nodes_dmp[1] == txid, 0]
                children_lengths = [
                    float(g_taxid[str(child)]) for child in children if str(child) in g_taxid.index
                ]
                if children_lengths:
                    lengths[str(txid)] = float(np.mean(children_lengths))
        lengths_out = np.array([lengths[str(taxid)] for taxid in sorted_keys], dtype=float)

    if np.any(np.isnan(lengths_out)):
        nan_taxids = sorted_keys[np.isnan(lengths_out)]
        print(
            f"WARNING: COULD NOT ESTIMATE GENOME SIZE FOR TAXIDS: {nan_taxids} WILL GIVE THEM MEAN LENGTH"
        )
        mean_len = float(np.nanmean(lengths_out)) if not np.all(np.isnan(lengths_out)) else 1.0
        for txid in nan_taxids:
            lengths[str(txid)] = mean_len
        lengths_out = np.array([lengths[str(taxid)] for taxid in sorted_keys], dtype=float)

    return lengths_out / np.sum(lengths_out)


def load_gc_distributions(kraken_db, nbin, read_len, fragment_len=None, return_raw=False):
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
    return_raw : bool
        If True, also return the raw DataFrame (with ``seqids`` and ``taxid``)
        so ``get_genome_length`` can reuse it without reading the CSV again.

    Returns
    -------
    tuple
        ``(gc_dists, insert)`` or ``(gc_dists, insert, gc_dists_raw)`` if ``return_raw=True``.
    """
    if fragment_len is not None:
        insert = fragment_len - 2 * read_len
        fname = os.path.join(
            kraken_db, f'gc_bin_{nbin}_kmer_{read_len}_insert_{insert}_dist.csv'
        )
    else:
        insert = None
        fname = os.path.join(
            kraken_db, f'gc_bin_{nbin}_kmer_{read_len}_dist.csv'
        )

    gc_dists_raw = pd.read_csv(fname, index_col=0)
    gc_dists = gc_dists_raw.drop('seqids', axis=1).groupby('taxid').sum().T
    gc_dists.columns = [str(i) for i in gc_dists.columns]

    if return_raw:
        return gc_dists, insert, gc_dists_raw
    return gc_dists, insert
