import os, sys
from datetime import datetime
from time import localtime, strftime
import pandas as pd
import numpy as np
from math import floor
from numpy.random import uniform
import Bio.SeqIO as SeqIO
from multiprocessing import Pool
from functools import partial
from glob import glob
import argparse


def print_time():
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)


def parsing(sequence):  # creates list with comulatative GC-content
    gc = [0]
    valid_nucleotides = [0]
    for base in sequence:
        if base in ['G', 'C']:
            gc.append(gc[-1] + 1)
        else:
            gc.append(gc[-1])
        if base in ['G', 'C', 'T', 'A']:
            valid_nucleotides.append(valid_nucleotides[-1] + 1)
        else:
            valid_nucleotides.append(valid_nucleotides[-1])

    return gc, valid_nucleotides


def create_refdist(seqid, kmer_len, paths, nbin=100, sampling=50):
    record_dict = SeqIO.index_db('fasta_index.sql', paths, 'fasta')
    seq = record_dict[seqid].seq
    occ = np.zeros(nbin + 1)
    gc_list, valid_nucleotides_list = parsing(seq)
    for i in range(len(seq) - kmer_len + 1)[::sampling]:
        gc = gc_list[i + kmer_len] - gc_list[i]
        valid_nucleotides = valid_nucleotides_list[i + kmer_len] - valid_nucleotides_list[i]
        if valid_nucleotides == 0:
            continue
        gc_content = floor((nbin / valid_nucleotides) * (gc + uniform()))
        if gc_content > nbin:
            continue
        occ[gc_content] += 1
    return occ


def create_refdist_insert(seqid, kmer_len, paths, nbin=100, sampling=50, insert=0):
    record_dict = SeqIO.index_db('fasta_index.sql', paths, 'fasta')
    seq = record_dict[seqid].seq
    occ = np.zeros(nbin + 1)
    gc_list, valid_nucleotides_list = parsing(seq)
    for i in range(len(seq) - 2*kmer_len - insert + 1)[::sampling]:
        gc_1 = gc_list[i + kmer_len] - gc_list[i]
        gc_2 = gc_list[i + 2*kmer_len + insert] - gc_list[i + kmer_len + insert]
        valid_nucleotides_1 = valid_nucleotides_list[i + kmer_len] - valid_nucleotides_list[i]
        valid_nucleotides_2 = valid_nucleotides_list[i + 2*kmer_len + insert] - valid_nucleotides_list[i + kmer_len + insert]
        valid_nucleotides = valid_nucleotides_1 + valid_nucleotides_2
        if valid_nucleotides == 0:
            continue
        gc = gc_1 + gc_2
        gc_content = floor((nbin / valid_nucleotides) * (gc + uniform()))
        if gc_content > nbin:
            continue
        occ[gc_content] += 1
    return occ


def main():

    parser = argparse.ArgumentParser(
        description='Create reference distributions with given read length and possibly fragment length to run GuaCAMOLE')

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

    if fragment_len is not None:
        insert = fragment_len - 2 * read_len

    sampling = 50
    nbin = 100
    time = strftime("%m-%d-%Y %H:%M:%S", localtime())
    sys.stdout.write("PROGRAM START TIME: " + time + '\n')
    if old_gc_dist is not None:
        old_df = pd.read_csv(old_gc_dist, index_col=0)

    os.chdir(lib_path)
    seqid2taxid = pd.read_csv('seqid2taxid.map', sep='\t', header=None)
    fasta_paths = glob('library/*/*.fna')
    # SeqIO.index_db seems to be sensitive to file order, so make sure it is reproducible
    fasta_paths.sort()
    print('Indexing fasta files...')
    record_dict = SeqIO.index_db('fasta_index.sql', fasta_paths, 'fasta')

    print("Indexing complete, querying sequence IDs")
    seqids = [x for x in np.array(seqid2taxid[0]) if x in record_dict.keys()]

    if old_gc_dist is not None:
        old_dist_seqids = np.array(old_df['seqids'])
        old_seqids = [x for x in seqids if x in old_dist_seqids]
        seqids = [x for x in seqids if x not in old_dist_seqids]
        print("Existing GC Distribution file provided!")
        print("Found GC distributions for " + str(len(old_seqids)) + " genome sequences in provided file")

    print("Computing GC distributions for " + str(len(seqids)) + " genome sequences")
    if fragment_len is not None:
        with Pool(ncores) as pool:
            dist_seqids = pool.map(
                partial(create_refdist_insert, kmer_len=read_len, nbin=nbin, sampling=sampling, paths=fasta_paths,
                        insert=insert),
                seqids
            )
    else:
        with Pool(ncores) as pool:
            dist_seqids = pool.map(
                partial(create_refdist, kmer_len=read_len, nbin=nbin, sampling=sampling, paths=fasta_paths),
                seqids
            )

    time = strftime("%m-%d-%Y %H:%M:%S", localtime())
    sys.stdout.write("PROGRAM END TIME: " + time + '\n')


    in_seqs = seqid2taxid.iloc[:, 0].isin(seqids)
    txids = seqid2taxid.loc[in_seqs, 1]

    df = pd.DataFrame(dist_seqids)
    df['taxid'] = np.array(txids)
    df['seqids'] = np.array(seqids)
    df.columns = [str(i) for i in df.columns]
    if old_gc_dist is not None:
        df = pd.concat([old_df.loc[old_df['seqids'].isin(old_seqids), :], df], axis=0)
    if fragment_len is not None:
        df.to_csv('gc_bin_' + str(nbin) + '_kmer_'+ str(read_len) + '_insert_' + str(insert) + '_dist.csv')
    else:
        df.to_csv('gc_bin_' + str(nbin) + '_kmer_'+ str(read_len) + '_dist.csv')

if __name__ == "__main__":
    main()

