import argparse
import os
import sys
from datetime import datetime
from functools import partial
from glob import glob
from multiprocessing import Pool
from time import localtime, strftime

import numpy as np
import pandas as pd
from Bio import SeqIO

_WORKER_RECORD_DICT = None


def _init_worker(paths):
    global _WORKER_RECORD_DICT
    _WORKER_RECORD_DICT = SeqIO.index_db('fasta_index.sql', paths, 'fasta')


def print_time():
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)


def parsing(sequence):
    """Compute cumulative GC and valid-nucleotide (A/C/G/T) counts (case-insensitive)."""
    seq_bytes = str(sequence).upper().encode('ascii', errors='ignore')
    arr = np.frombuffer(seq_bytes, dtype=np.uint8)
    is_gc = (arr == ord('G')) | (arr == ord('C'))
    is_valid = is_gc | (arr == ord('A')) | (arr == ord('T'))

    gc = np.empty(len(arr) + 1, dtype=np.int64)
    gc[0] = 0
    np.cumsum(is_gc, out=gc[1:])

    valid_nucleotides = np.empty(len(arr) + 1, dtype=np.int64)
    valid_nucleotides[0] = 0
    np.cumsum(is_valid, out=valid_nucleotides[1:])

    return gc, valid_nucleotides


def _get_sequence(seqid, paths):
    record_dict = (
        _WORKER_RECORD_DICT
        if _WORKER_RECORD_DICT is not None
        else SeqIO.index_db('fasta_index.sql', paths, 'fasta')
    )
    return record_dict[seqid].seq


def create_refdist(seqid, kmer_len, paths, nbin=100, sampling=50):
    seq = _get_sequence(seqid, paths)
    seq_len = len(seq)
    if seq_len < kmer_len:
        return np.zeros(nbin + 1, dtype=float)

    gc_list, valid_list = parsing(seq)
    starts = np.arange(0, seq_len - kmer_len + 1, sampling)
    gc = gc_list[starts + kmer_len] - gc_list[starts]
    valid = valid_list[starts + kmer_len] - valid_list[starts]

    keep = valid > 0
    if not np.any(keep):
        return np.zeros(nbin + 1, dtype=float)

    gc_k = gc[keep]
    valid_k = valid[keep]
    bins = np.floor((nbin / valid_k) * (gc_k + np.random.random(len(gc_k)))).astype(np.int64)
    bins = bins[(bins >= 0) & (bins <= nbin)]
    return np.bincount(bins, minlength=nbin + 1).astype(float)


def create_refdist_insert(seqid, kmer_len, paths, nbin=100, sampling=50, insert=0):
    seq = _get_sequence(seqid, paths)
    frag_len = 2 * kmer_len + insert
    seq_len = len(seq)
    if seq_len < frag_len:
        return np.zeros(nbin + 1, dtype=float)

    gc_list, valid_list = parsing(seq)
    starts = np.arange(0, seq_len - frag_len + 1, sampling)
    gc_1 = gc_list[starts + kmer_len] - gc_list[starts]
    gc_2 = gc_list[starts + frag_len] - gc_list[starts + kmer_len + insert]
    valid_1 = valid_list[starts + kmer_len] - valid_list[starts]
    valid_2 = valid_list[starts + frag_len] - valid_list[starts + kmer_len + insert]

    valid = valid_1 + valid_2
    keep = valid > 0
    if not np.any(keep):
        return np.zeros(nbin + 1, dtype=float)

    gc_k = (gc_1 + gc_2)[keep]
    valid_k = valid[keep]
    bins = np.floor((nbin / valid_k) * (gc_k + np.random.random(len(gc_k)))).astype(np.int64)
    bins = bins[(bins >= 0) & (bins <= nbin)]
    return np.bincount(bins, minlength=nbin + 1).astype(float)


def main():
    parser = argparse.ArgumentParser(
        description='Create reference distributions with given read length and possibly fragment length to run GuaCAMOLE'
    )

    parser.add_argument('--lib_path', metavar='lib_path', type=str, help='Path to Kraken2 database')
    parser.add_argument('--read_len', metavar='read_len', type=int, help='read length')
    parser.add_argument('--ncores', metavar='ncores', type=int, help='number of threads')
    parser.add_argument('--fragment_len', metavar='fragment_len', type=int, help='fragment length', default=None)
    args = parser.parse_args()
    lib_path = args.lib_path
    read_len = args.read_len
    ncores = args.ncores
    fragment_len = args.fragment_len
    old_gc_dist = None

    insert = (fragment_len - 2 * read_len) if fragment_len is not None else None

    sampling = 50
    nbin = 100
    cur_time = strftime("%m-%d-%Y %H:%M:%S", localtime())
    sys.stdout.write(f"PROGRAM START TIME: {cur_time}\n")
    if old_gc_dist is not None:
        old_df = pd.read_csv(old_gc_dist, index_col=0)

    os.chdir(lib_path)
    seqid2taxid = pd.read_csv('seqid2taxid.map', sep='\t', header=None)
    fasta_paths = sorted(glob('library/*/*.fna'))
    print('Indexing fasta files...')
    record_dict = SeqIO.index_db('fasta_index.sql', fasta_paths, 'fasta')

    print("Indexing complete, querying sequence IDs")
    seqid_to_taxid = dict(zip(seqid2taxid[0], seqid2taxid[1], strict=True))
    seqids = [x for x in seqid2taxid[0].to_numpy() if x in record_dict]

    if old_gc_dist is not None:
        old_dist_seqids = set(old_df['seqids'].to_numpy())
        old_seqids = [x for x in seqids if x in old_dist_seqids]
        seqids = [x for x in seqids if x not in old_dist_seqids]
        print("Existing GC Distribution file provided!")
        print(f"Found GC distributions for {len(old_seqids)} genome sequences in provided file")

    print(f"Computing GC distributions for {len(seqids)} genome sequences")
    if fragment_len is not None:
        worker_fn = partial(
            create_refdist_insert, kmer_len=read_len, nbin=nbin,
            sampling=sampling, paths=fasta_paths, insert=insert,
        )
    else:
        worker_fn = partial(
            create_refdist, kmer_len=read_len, nbin=nbin,
            sampling=sampling, paths=fasta_paths,
        )

    with Pool(ncores, initializer=_init_worker, initargs=(fasta_paths,)) as pool:
        dist_seqids = pool.map(worker_fn, seqids)

    cur_time = strftime("%m-%d-%Y %H:%M:%S", localtime())
    sys.stdout.write(f"PROGRAM END TIME: {cur_time}\n")

    txids = [seqid_to_taxid[s] for s in seqids]
    df = pd.DataFrame(dist_seqids)
    df['taxid'] = np.array(txids)
    df['seqids'] = np.array(seqids)
    df.columns = [str(i) for i in df.columns]
    if old_gc_dist is not None:
        df = pd.concat([old_df.loc[old_df['seqids'].isin(old_seqids), :], df], axis=0)
    if fragment_len is not None:
        df.to_csv(f'gc_bin_{nbin}_kmer_{read_len}_insert_{insert}_dist.csv')
    else:
        df.to_csv(f'gc_bin_{nbin}_kmer_{read_len}_dist.csv')


if __name__ == "__main__":
    main()
