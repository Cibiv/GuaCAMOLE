import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import os
import seaborn as sns
import qpsolvers as qps
import gzip
import scipy
from math import floor
from numpy.random import uniform


def print_time():
    now = datetime.now()
    current_time = now.strftime('%H:%M:%S')
    print("Current Time =", current_time)


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i / inch for i in tupl[0])
    else:
        return tuple(i / inch for i in tupl)


def trim_dist(dist, quantiles=None):
    """
    cuts off bins on extreme ends of gc dist
    :param dist: array with gc dist
    :param quantiles: quantiles within which values are kept
    :return: array with trimmed gc distribution
    """

    if quantiles is None:
        quantiles = [0.05, 0.95]

    dist = np.array(dist)
    cdf = np.cumsum(dist)
    cdf = cdf / cdf[-1]
    lower_bound = np.searchsorted(cdf, min(quantiles))
    upper_bound = np.searchsorted(cdf, max(quantiles))

    new_dist = dist.copy()
    new_dist[:lower_bound] = 0
    new_dist[upper_bound + 1:] = 0

    return new_dist


def normalize_gc_dists(sample_dist, ref_dist, trim=True, quantiles=None):
    """
    normalizes a read gc distribution with an expected distributiom from
    reference sequences
    :param sample_dist: path to the sample distribution (from distribute_reads() )
    :param ref_dist: path to the reference distribution (from get_refgenomes() )
    :param trim: cuts off reads at extreme ends of GC dist (often outliers)
    :param quantiles:
    :return: matrix with the normalized distributions and reference distributions
    """

    sample_dist = np.loadtxt(sample_dist)
    ref_dist = np.loadtxt(ref_dist)
    if trim:
        sample_dist = np.array(pd.DataFrame(sample_dist).apply(trim_dist, result_type='broadcast',
                                                               quantiles=quantiles))
        ref_dist = np.array(pd.DataFrame(ref_dist).apply(trim_dist, result_type='broadcast',
                                                         quantiles=quantiles))
        sample_dist[np.where(ref_dist == 0)] = 0

    norm = sample_dist / ref_dist
    norm[np.isinf(norm)] = 0
    norm[np.isnan(norm)] = 0
    return norm, ref_dist, sample_dist


def get_overlaps(norm_dist):
    """
    return list of list with overlaps of gc distributions of organisms
    :param norm_dist: matrix with normalized gc distributions of organisms
    :return: list of lists containing arrays with pairwise overlaps of organisms
    """
    overlaps = []
    for i in range(norm_dist.shape[1]):
        org1 = np.where(norm_dist[:, i] != 0)[0]
        overlaps.append([])
        for j in range(norm_dist.shape[1]):
            org2 = np.where(norm_dist[:, j] != 0)[0]
            overlaps[i].append([x for x in org1 if x in org2])

    return overlaps


def get_beta(norm_dist, org1, org2, all_overlaps, reg_weight=0.1):
    """
    returns beta_ij for LS matrix
    :param norm_dist: matrix with normalized gc distributions of organisms
    :param ref_dist: reference distributions
    :param org1: first org (index!)
    :param org2: second org (index!)
    :param all_overlaps: list of list as from get_overlaps()
    :return:
    """
    gc = all_overlaps[org1][org2]
    if len(gc) == 0:
        first_term = 0
    else:
        first_term = -np.sum(norm_dist[gc, org1] * norm_dist[gc, org2])

    n_b = (norm_dist != 0).sum(1)
    regularization = 0
    for bin_i in range(norm_dist.shape[0]):
        for bin_j in range(norm_dist.shape[0]):
            if abs(bin_j - bin_i) > 10 or n_b[bin_i] == 0 or n_b[bin_j] == 0:
                continue
            regularization += ((norm_dist[bin_i, org1] * norm_dist[bin_i, org2] / np.square(n_b[bin_i]) +
                                norm_dist[bin_j, org1] * norm_dist[bin_j, org2] / np.square(n_b[bin_j]) -
                                2 * norm_dist[bin_i, org1] * norm_dist[bin_j, org2] / (n_b[bin_i] * n_b[bin_j]))
                               * rbf(bin_i, bin_j))

    return (1 - reg_weight) * first_term / float(norm_dist.shape[1] ** 2) + regularization * reg_weight


def get_alpha(norm_dist, org, all_overlaps, reg_weight=0.1):
    """
    returns alpha_i for LS matrix
    :param norm_dist: matrix with normalized gc distributions of organisms
    :param ref_dist: reference distributions
    :param org: organism index (integer)
    :param all_overlaps: list of list as from get_overlaps()
    :return: alpha_i
    """
    overlaps = all_overlaps[org]
    inner_sum = 0
    for i in range(len(overlaps)):
        if i == org:
            continue
        xsq = np.square(norm_dist[overlaps[i], org])
        inner_sum += np.sum(xsq)

    n_b = (norm_dist != 0).sum(1)
    regularization = 0
    for bin_i in range(norm_dist.shape[0]):
        for bin_j in range(norm_dist.shape[0]):
            if abs(bin_j - bin_i) > 10 or n_b[bin_i] == 0 or n_b[bin_j] == 0:
                continue
            reg_ij = 0
            reg_ij += np.square(norm_dist[bin_i, org]) / np.square(n_b[bin_i])
            reg_ij += np.square(norm_dist[bin_j, org]) / np.square(n_b[bin_j])
            reg_ij -= (2 * norm_dist[bin_i, org] * norm_dist[bin_j, org]) / (n_b[bin_i] * n_b[bin_j])
            regularization += reg_ij * rbf(bin_i, bin_j)

    return (1 - reg_weight) * inner_sum / float(norm_dist.shape[1] ** 2) + regularization * reg_weight


def get_matrix(norm_dist, reg_weight=0.1):
    """
    return matrix for no_mapping pipeline to minimize LS function
    :param norm_dist: matrix with normalized gc distributions of organisms
    :return: matrix
    """
    all_overlaps = get_overlaps(norm_dist=norm_dist)
    a = np.zeros(len(all_overlaps))
    for i in range(len(all_overlaps)):
        a[i] = get_alpha(norm_dist=norm_dist, org=i, reg_weight=reg_weight,
                         all_overlaps=all_overlaps)

    b = np.zeros((len(all_overlaps), len(all_overlaps)))
    for row in range(len(all_overlaps)):
        for col in range(len(all_overlaps)):

            if row == col:
                b[row, col] = a[row]

            elif row > col:
                b[row, col] = b[col, row]

            else:
                b[row, col] = get_beta(norm_dist=norm_dist, reg_weight=reg_weight,
                                       org1=row, org2=col, all_overlaps=all_overlaps)
    return scipy.sparse.csc_matrix(b)


def get_efficiencies(ref_dists, sample_dists, abundances, iteration=0, plot=False):
    np.seterr(divide='ignore')
    lib_size = np.sum(sample_dists)
    r_dists_scaled = abundances * ref_dists / np.sum(abundances * ref_dists)
    ef = sample_dists / (r_dists_scaled * lib_size)
    ef[np.isnan(ef)] = 0
    ef[np.isinf(ef)] = 0
    # weighted mean
    s_dists_norm = np.log(sample_dists)
    s_dists_norm[ef == 0] = 0
    s_dists_norm[np.isinf(s_dists_norm)] = 0
    ef_weights = s_dists_norm / np.sum(s_dists_norm, 1).reshape(s_dists_norm.shape[0], 1)
    w_ef = ef * ef_weights
    filename = 'Efficiencies_it_' + str(iteration)
    w_ef[np.isnan(w_ef)] = 0
    w_ef = w_ef.sum(1)
    if plot:
        plt.figure(figsize=cm2inch(20, 13))
        plt.rcParams.update({'font.size': 16})
        plot_dist2(w_ef.reshape(w_ef.shape[0], 1) / np.max(w_ef))
        plt.ylim(0, 1.1)
        plt.xlabel('GC bin [%]')
        plt.ylabel('Relative Sequencing Efficiency')
        plt.savefig(filename + '.pdf')
        plt.close()

    return w_ef, ef


def get_residuals(ref_dists, sample_dists, abundances, iteration=None, plot=False):
    ef_mean, ef_all = get_efficiencies(ref_dists=ref_dists, sample_dists=sample_dists, abundances=abundances, plot=plot,
                                       iteration=iteration)

    r_dists_scaled = abundances * ref_dists * ef_mean.reshape((ef_mean.shape[0], 1))
    norm_constant = np.sum(abundances * ref_dists * ef_mean.reshape((ef_mean.shape[0], 1)))
    expected = np.sum(sample_dists) * r_dists_scaled / norm_constant
    expected[sample_dists == 0] = 0

    residuals = (sample_dists - expected) / expected
    res_range = np.nanmax(residuals, axis=0) - np.nanmin(residuals, axis=0)

    return res_range, residuals, ef_mean


def corrected_abundances(sample_path, reference_path, quantiles=None, taxids=None, plot=False, reg_weight=0.1,
                         fp_cycles=3):
    """
    return non gc biased abundances
    :param sample_path: path to sample distribution file
    :param reference_path: path to reference distribution file
    :return: lambdas (ordered by txid of organism)
    """
    norm_dist, ref_dist, sample_dist = normalize_gc_dists(sample_path, reference_path, trim=True,
                                                          quantiles=quantiles)

    skipped_taxids = []
    taxids = np.array(taxids)
    cov_rel = np.zeros(norm_dist.shape)
    taxon_removal_cycle = np.repeat(np.nan, len(taxids))
    res_range = np.zeros(len(taxids))
    threshold = 10
    final_threshold = threshold / np.power(2, fp_cycles-1)
    iteration = 1
    rep = 1
    ind_skip = []

    while threshold >= final_threshold or len(ind_skip) != 0:

        ind = np.where(np.isin(taxids, skipped_taxids, invert=True))[0]

        # only re-compute matrix when taxa have been removed or it's initial round
        if len(ind_skip) != 0 or (iteration == 1 and rep == 1):
            m = get_matrix(norm_dist[:, ind], reg_weight=reg_weight)

        a = scipy.sparse.csc_matrix(np.full((m.shape[0],), 1.))
        ab = qps.solve_qp(P=m, q=np.zeros(m.shape[0]), A=a, b=np.array([1.]), solver='cvxopt', verbose=False,
                          abstol=1e-15, reltol=1e-15, lb=np.zeros(m.shape[0]))

        ab = (1 / ab) / np.sum(1 / ab)

        if plot:
            plt.figure(figsize=cm2inch(20, 13))
            plt.rcParams.update({'font.size': 16})
            plot_dist2(norm_dist[:, ind] / ab, line=True)
            plt.xlabel("GC bin (%)")
            plt.ylabel("Obs / (Exp*Abundance)")
            plt.savefig("Obs_vs_Exp_minimized_it_" + str(iteration) + '_rep_' + str(rep) + ".pdf")
            plt.close()

        res_range_new, cov_rel[:, ind], efficiencies = get_residuals(ref_dists=ref_dist[:, ind],
                                                                     sample_dists=sample_dist[:, ind],
                                                                     abundances=ab, plot=plot,
                                                                     iteration=iteration)

        res_range[ind] = res_range_new
        ind_skip = np.where(res_range_new > threshold)[0]
        skipped_taxids.extend(taxids[ind][ind_skip])
        taxon_removal_cycle[ind[ind_skip]] = iteration

        print(
            "Removing taxa " + str(taxids[ind][ind_skip]) + " in cycle " + str(iteration) + " out of " + str(fp_cycles))
        print(f"and with threshold {str(threshold)}")

        residual_df = pd.DataFrame(
            {
                'taxid': taxids[ind],
                'res_range': res_range_new
            }
        )
        residual_df.to_csv('residuals_cycle_' + str(iteration) + '.csv')

        if threshold == final_threshold and len(ind_skip) == 0:
            if plot:
                residual_plot(res_array=cov_rel[:, ind], taxids=taxids[ind], threshold_plot=None)
                plt.savefig('Residuals_final.pdf')
                plt.close()
                os.rename('Efficiencies_it_' + str(iteration) + '.pdf', 'efficiencies_final.pdf')
            break

        rep += 1
        if threshold > final_threshold and len(ind_skip) == 0:
            iteration += 1
            threshold = threshold / 2
            rep = 1

    all_ab = np.zeros(len(taxids))
    all_ab[np.where(np.isin(taxids, skipped_taxids, invert=True))[0]] = ab
    all_ab[np.where(np.isin(taxids, skipped_taxids))[0]] = 0
    return all_ab, taxon_removal_cycle, efficiencies


def plot_dist2(dists, color=None, legend=False, line=False, **kwargs):
    if len(dists.shape) == 1:
        dists = dists.reshape((len(dists), 1))

    dists = dists[:, np.where(np.nansum(dists, axis=0) != 0)[0]]

    dist_df = pd.DataFrame(
        {
            'gc': np.tile(np.arange(0, dists.shape[0]), dists.shape[1]),
            'y': dists.flatten('F'),
            'org': [str(i) for i in np.repeat(np.arange(0, dists.shape[1]), dists.shape[0])]
        }
    )
    if color is None:
        if line:
            return sns.lineplot(x='gc', y='y', hue='org', data=dist_df.loc[dist_df['y'] != 0], legend=False,
                                palette='deep', **kwargs)
        return sns.scatterplot(x='gc', y='y', hue='org', style='org', data=dist_df.loc[dist_df['y'] != 0], legend=False,
                               palette='deep', **kwargs)

    else:
        dist_df['category'] = np.repeat(color, dists.shape[0])
        return sns.lineplot(x='gc', y='y', hue='category', style='org', data=dist_df.loc[dist_df['y'] != 0],
                            legend=legend, **kwargs)


def residual_plot(res_array, taxids, threshold_plot=None):
    gc_bin = []
    res_gc = []
    taxonomy = []
    for i in range(res_array.shape[1]):
        sub_array = res_array[:, i]
        gc_bin.extend(np.where(~np.isnan(sub_array))[0])
        res_gc.extend(sub_array[np.where(~np.isnan(sub_array))[0]])
        taxonomy.extend(np.repeat(taxids[i], len(np.where(~np.isnan(sub_array))[0])))
    res_df = pd.DataFrame(
        {
            'gc': gc_bin,
            'residual': res_gc,
            'taxid': taxonomy
        }
    )
    palette = sns.color_palette(['black'], len(res_df['taxid'].unique()))
    sns.lineplot(x='gc', y='residual', hue='taxid', data=res_df, legend=False, palette=palette, alpha=0.4)
    if threshold_plot is not None:
        plt.hlines(y=threshold_plot, xmin=min(res_df['gc']), xmax=max(res_df['gc']))
    plt.ylabel('Relative Residuals')
    plt.xlabel('GC bin [%]')


def create_sample_dist(fastq1, txids, kraken, fastq2=None, nbin=100, fasta=False):
    denom = 4
    j = 4
    if fasta:
        j = 2
        denom = 2
    i = 0

    if fastq1[-2:] == "gz":
        f1 = gzip.open(fastq1, 'rt')
        if fastq2 is not None:
            f2 = gzip.open(fastq2, 'rt')
    else:
        f1 = open(fastq1, 'r')
        if fastq2 is not None:
            f2 = open(fastq2, 'r')
    occ = np.zeros((nbin + 1, len(txids)), dtype=int)
    idx = dict(zip(txids, range(len(txids))))
    kr = open(kraken, 'r')

    while kr:
        ln_fq1 = f1.readline().strip()
        if fastq2 is not None:
            ln_fq2 = f2.readline().strip()

        if i % j == 0:
            ln_kr = kr.readline()

        if ln_fq1 == "" or ln_kr.strip() == "":
            break
        current_txid = ln_kr.split(sep="\t")[2]

        if int(current_txid) in txids:
            if (i - 1) % denom == 0 or i == 1:
                gc_1 = compute_gc(ln_fq1, nbin=nbin)
                if f2 is not None:
                    gc_2 = compute_gc(ln_fq2, nbin=nbin)
                    if max([gc_1, gc_2]) > nbin:
                        i += 1
                        continue
                    occ[gc_2, idx[int(current_txid)]] += 1
                if gc_1 > nbin:
                    i += 1
                    continue
                occ[gc_1, idx[int(current_txid)]] += 1

        i += 1
    f1.close()
    if f2 is not None:
        f2.close()
    kr.close()

    header = "\t".join(map(str, txids))
    np.savetxt("sample_bin_" + str(nbin) + ".dist", occ,
               fmt='%i', delimiter="\t", header=header)


def compute_gc(seq, nbin=100):
    gc = seq.upper().count('G') + seq.upper().count('C')
    atgc = seq.count('A') + seq.count('T') + gc
    if atgc == 0:
        return 0
    return floor((nbin / atgc) * (gc + uniform()))


def rbf(val1, val2):
    w = np.exp(-abs(val1 - val2))
    return w


def create_sample_dist_avg(fastq1, txids, kraken, fastq2, nbin=100, fasta=False):
    denom = 4
    if fasta:
        denom = 2
    i = 0

    if fastq1[-2:] == "gz":
        f1 = gzip.open(fastq1, 'rt')
        f2 = gzip.open(fastq2, 'rt')
    else:
        f1 = open(fastq1, 'r')
        f2 = open(fastq2, 'r')
    occ = np.zeros((nbin + 1, len(txids)), dtype=int)
    idx = dict(zip(txids, range(len(txids))))
    kr = open(kraken, 'r')

    while kr:
        ln_fq1 = f1.readline().strip()
        ln_fq2 = f2.readline().strip()

        if i % denom == 0:
            ln_kr = kr.readline()

        if ln_fq1 == "" or ln_kr.strip() == "":
            break
        current_txid = ln_kr.split(sep="\t")[2]

        if int(current_txid) in txids:
            if (i - 1) % denom == 0 or i == 1:
                gc = compute_gc(ln_fq1 + ln_fq2, nbin=nbin)
                if gc > nbin:
                    i += 1
                    continue
                occ[gc, idx[int(current_txid)]] += 1
        i += 1
    f1.close()
    f2.close()
    kr.close()

    header = "\t".join(map(str, txids))
    np.savetxt("sample_bin_" + str(nbin) + ".dist", occ,
               fmt='%i', delimiter="\t", header=header)
