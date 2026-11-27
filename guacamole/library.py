import bz2
import gzip
import json
import os
import random
from contextlib import ExitStack
from datetime import datetime
from math import floor

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import qpsolvers as qps
import scipy.sparse
import seaborn as sns


def print_time():
    now = datetime.now()
    current_time = now.strftime('%H:%M:%S')
    print("Current Time =", current_time)


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i / inch for i in tupl[0])
    return tuple(i / inch for i in tupl)


def trim_dist(dist, quantiles=None):
    """Cut off bins on extreme ends of a 1D GC distribution.

    Parameters
    ----------
    dist : array-like
        1D array with GC distribution counts/probabilities.
    quantiles : sequence of float, optional
        Lower and upper quantiles within which values are kept.
        Defaults to ``[0.05, 0.95]``.

    Returns
    -------
    np.ndarray
        Trimmed GC distribution.
    """
    if quantiles is None:
        quantiles = [0.05, 0.95]

    dist = np.asarray(dist, dtype=float)
    new_dist = dist.copy()
    total = np.sum(dist)
    if total <= 0 or np.isnan(total):
        return new_dist

    cdf = np.cumsum(dist) / total
    lower_bound = int(np.searchsorted(cdf, min(quantiles)))
    upper_bound = int(np.searchsorted(cdf, max(quantiles)))

    new_dist[:lower_bound] = 0
    new_dist[upper_bound + 1:] = 0

    return new_dist


def trim_dist_matrix(dist, quantiles=None):
    """Trim every column of a 2D GC distribution matrix to quantile bounds."""
    if quantiles is None:
        quantiles = [0.05, 0.95]

    dist = np.asarray(dist, dtype=float)
    if dist.ndim == 1:
        return trim_dist(dist, quantiles=quantiles)

    q_low = min(quantiles)
    q_high = max(quantiles)
    new_dist = dist.copy()
    totals = np.sum(dist, axis=0)
    for j in range(dist.shape[1]):
        total = totals[j]
        if total <= 0 or np.isnan(total):
            continue
        cdf = np.cumsum(dist[:, j]) / total
        lower_bound = int(np.searchsorted(cdf, q_low))
        upper_bound = int(np.searchsorted(cdf, q_high))
        new_dist[:lower_bound, j] = 0
        new_dist[upper_bound + 1:, j] = 0

    return new_dist


def _load_dist_array(dist_or_path):
    """Load a distribution matrix from a file path or convert an array to 2D float64."""
    arr = np.asarray(dist_or_path, dtype=float) if isinstance(dist_or_path, np.ndarray) else np.loadtxt(dist_or_path)
    arr = np.asarray(arr, dtype=float)
    if arr.ndim == 1:
        arr = arr[:, None]
    return arr


def normalize_gc_dists(sample_dist, ref_dist, trim=True, quantiles=None):
    """Normalize observed read GC distributions by expected reference distributions.

    Parameters
    ----------
    sample_dist : str, os.PathLike, or np.ndarray
        Path to the sample distribution file or 2D array (n_bins, n_taxa).
    ref_dist : str, os.PathLike, or np.ndarray
        Path to the reference distribution file or 2D array (n_bins, n_taxa).
    trim : bool
        If True, trim tails outside ``quantiles``.
    quantiles : sequence of float, optional
        Lower and upper quantiles for trimming.

    Returns
    -------
    tuple of np.ndarray
        ``(norm, ref_dist, sample_dist)`` matrices of shape ``(n_bins, n_taxa)``.
    """
    sample_dist = _load_dist_array(sample_dist)
    ref_dist = _load_dist_array(ref_dist)
    if trim:
        sample_dist = trim_dist_matrix(sample_dist, quantiles=quantiles)
        ref_dist = trim_dist_matrix(ref_dist, quantiles=quantiles)
        sample_dist[ref_dist == 0] = 0

    norm = np.divide(
        sample_dist,
        ref_dist,
        out=np.zeros_like(sample_dist, dtype=float),
        where=ref_dist != 0,
    )
    norm[~np.isfinite(norm)] = 0
    return norm, ref_dist, sample_dist


def get_overlaps(norm_dist):
    """Return pairwise GC-bin overlaps of organisms as a nested list."""
    nonzero = norm_dist != 0
    n_taxa = norm_dist.shape[1]
    overlaps = []
    for i in range(n_taxa):
        row = []
        for j in range(n_taxa):
            row.append(np.flatnonzero(nonzero[:, i] & nonzero[:, j]).tolist())
        overlaps.append(row)
    return overlaps


def rbf(val1, val2):
    return float(np.exp(-abs(val1 - val2)))


def _compute_reg_coefficients(norm_dist, sample_dist=None):
    """Compute matrix C where lambda = C @ x (Eq. 9 of Holcik et al., 2025).

    If ``sample_dist`` is provided, taxon-specific efficiencies in bin g are
    weighted by ``log(O(i, g) + 1)`` and normalized by the sum of weights
    ``n_g = sum_i log(O(i, g) + 1)`` over taxa with nonzero ``norm_dist``.
    Otherwise, falls back to unweighted averaging over active taxa in bin g.
    """
    mask = norm_dist != 0
    if sample_dist is not None:
        weights = np.where(mask, np.log1p(np.maximum(sample_dist, 0.0)), 0.0)
    else:
        weights = mask.astype(float)

    weight_sums = weights.sum(axis=1)
    active = weight_sums > 0
    coeff = np.zeros_like(norm_dist, dtype=float)
    coeff[active] = (weights[active] / weight_sums[active, None]) * norm_dist[active]
    return coeff, active


def get_beta(norm_dist, org1, org2, all_overlaps, reg_weight=0.1, sample_dist=None):
    """Return off-diagonal quadratic-form coefficient Q[org1, org2] (Eqs. 8-10)."""
    gc = all_overlaps[org1][org2]
    first_term = 0.0 if len(gc) == 0 else -2.0 * float(np.sum(norm_dist[gc, org1] * norm_dist[gc, org2]))

    coeff, active = _compute_reg_coefficients(norm_dist, sample_dist=sample_dist)
    n_bins = norm_dist.shape[0]
    idx = np.arange(n_bins)
    diff = np.abs(idx[:, None] - idx[None, :])
    kernel = np.where(active[:, None] & active[None, :], np.exp(-diff), 0.0)
    laplacian = np.diag(kernel.sum(axis=1)) - kernel
    regularization = float(2.0 * (coeff[:, org1] @ laplacian @ coeff[:, org2]))

    n_taxa = norm_dist.shape[1]
    return (1.0 - reg_weight) * first_term / float(n_taxa ** 2) + reg_weight * regularization


def get_alpha(norm_dist, org, all_overlaps, reg_weight=0.1, sample_dist=None):
    """Return diagonal quadratic-form coefficient Q[org, org] (Eqs. 8-10)."""
    overlaps = all_overlaps[org]
    inner_sum = 0.0
    for i, bins in enumerate(overlaps):
        if i == org or len(bins) == 0:
            continue
        inner_sum += float(np.sum(np.square(norm_dist[bins, org])))

    coeff, active = _compute_reg_coefficients(norm_dist, sample_dist=sample_dist)
    n_bins = norm_dist.shape[0]
    idx = np.arange(n_bins)
    diff = np.abs(idx[:, None] - idx[None, :])
    kernel = np.where(active[:, None] & active[None, :], np.exp(-diff), 0.0)
    laplacian = np.diag(kernel.sum(axis=1)) - kernel
    regularization = float(2.0 * (coeff[:, org] @ laplacian @ coeff[:, org]))

    n_taxa = norm_dist.shape[1]
    return (1.0 - reg_weight) * (2.0 * inner_sum) / float(n_taxa ** 2) + reg_weight * regularization


def get_matrix(norm_dist, reg_weight=0.1, sample_dist=None):
    """Build the quadratic programming matrix P for minimizing Eq. (10).

    Constructs ``P`` such that ``0.5 * x^T P x = G_tilde(x)`` where:

        G_tilde(x) = ((1 - r) / n^2) * sum_{i != j} sum_g (R_{g,i} x_i - R_{g,j} x_j)^2
                     + r * sum_{k, l} exp(-|k - l|) * (lambda_k - lambda_l)^2

    with ``lambda = C @ x`` using Eq. (9) weights ``log(O(i, g) + 1)`` when
    ``sample_dist`` is provided.

    Parameters
    ----------
    norm_dist : np.ndarray
        Matrix ``R`` of shape ``(n_bins, n_taxa)`` with observed/expected ratios.
    reg_weight : float
        Regularization hyperparameter ``r`` in ``[0, 1]``.
    sample_dist : np.ndarray, optional
        Matrix ``O`` of shape ``(n_bins, n_taxa)`` with observed read counts
        used for the ``log(O(i, g) + 1)`` weights in Eq. (9).

    Returns
    -------
    scipy.sparse.csc_matrix
        Symmetric positive semidefinite QP matrix ``P`` of shape ``(n_taxa, n_taxa)``.
    """
    norm_dist = np.asarray(norm_dist, dtype=float)
    n_bins, n_taxa = norm_dist.shape
    if n_taxa == 0:
        return scipy.sparse.csc_matrix((0, 0), dtype=float)

    mask = norm_dist != 0
    m_g = mask.sum(axis=1)

    # Pairwise agreement quadratic form G(x) = x^T Q_G x (Eq. 8)
    # Off-diagonal: -2 * sum_g R_{g,i} * R_{g,j}
    # Diagonal:      2 * sum_g R_{g,i}^2 * max(m_g - 1, 0)
    q_g = -2.0 * (norm_dist.T @ norm_dist)
    diag_g = 2.0 * np.sum((norm_dist ** 2) * np.maximum(m_g - 1, 0)[:, None], axis=0)
    np.fill_diagonal(q_g, diag_g)

    # Regularization quadratic form R(x) = sum_{k,l} W_{k,l} (lambda_k - lambda_l)^2
    # = 2 * lambda^T L lambda = x^T (2 * C^T L C) x (Eqs. 9 & 10)
    coeff, active = _compute_reg_coefficients(norm_dist, sample_dist=sample_dist)
    idx = np.arange(n_bins)
    diff = np.abs(idx[:, None] - idx[None, :])
    kernel = np.where(active[:, None] & active[None, :], np.exp(-diff), 0.0)
    laplacian = np.diag(kernel.sum(axis=1)) - kernel
    q_reg = 2.0 * (coeff.T @ laplacian @ coeff)

    # qpsolvers minimizes 0.5 * x^T P x, so P = 2 * Q_total
    q_total = ((1.0 - reg_weight) / float(n_taxa ** 2)) * q_g + reg_weight * q_reg
    p_mat = 2.0 * (0.5 * (q_total + q_total.T))
    return scipy.sparse.csc_matrix(p_mat)


def get_efficiencies(ref_dists, sample_dists, abundances, iteration=0, plot=False):
    """Estimate GC-dependent sequencing efficiencies via Eq. (9).

    Computes the ``log(O(i, g) + 1)``-weighted average of taxon-specific
    efficiencies across active taxa in each GC bin.
    """
    lib_size = float(np.sum(sample_dists))
    ab_ref = abundances * ref_dists
    ab_ref_sum = float(np.sum(ab_ref))
    denom_ef = (ab_ref / ab_ref_sum) * lib_size if ab_ref_sum > 0 else np.zeros_like(ab_ref)
    ef = np.divide(
        sample_dists,
        denom_ef,
        out=np.zeros_like(sample_dists, dtype=float),
        where=denom_ef > 0,
    )
    ef[~np.isfinite(ef)] = 0.0

    # Eq. (9): weighted mean with weights log(O(i, g) + 1)
    weights = np.log1p(np.maximum(sample_dists, 0.0))
    weights[ef == 0] = 0.0
    weight_sums = weights.sum(axis=1, keepdims=True)
    ef_weights = np.divide(
        weights,
        weight_sums,
        out=np.zeros_like(weights, dtype=float),
        where=weight_sums > 0,
    )
    w_ef = np.sum(ef * ef_weights, axis=1)

    if plot:
        filename = f'Efficiencies_it_{iteration}'
        max_ef = float(np.max(w_ef)) if len(w_ef) > 0 else 0.0
        scaled_ef = (w_ef / max_ef) if max_ef > 0 else w_ef
        plt.figure(figsize=cm2inch(20, 13))
        plt.rcParams.update({'font.size': 16})
        plot_dist2(scaled_ef.reshape(scaled_ef.shape[0], 1))
        plt.ylim(0, 1.1)
        plt.xlabel('GC bin [%]')
        plt.ylabel('Relative Sequencing Efficiency')
        plt.savefig(filename + '.pdf')
        plt.close()

    return w_ef, ef


def get_residuals(ref_dists, sample_dists, abundances, iteration=None, plot=False):
    """Compute relative residuals (O - O_hat) / O_hat and per-taxon ranges."""
    ef_mean, _ = get_efficiencies(
        ref_dists=ref_dists,
        sample_dists=sample_dists,
        abundances=abundances,
        plot=plot,
        iteration=iteration,
    )

    r_dists_scaled = abundances * ref_dists * ef_mean[:, None]
    norm_constant = float(np.sum(r_dists_scaled))
    if norm_constant > 0:
        expected = float(np.sum(sample_dists)) * r_dists_scaled / norm_constant
    else:
        expected = np.zeros_like(r_dists_scaled)
    expected[sample_dists == 0] = 0.0

    residuals = np.divide(
        sample_dists - expected,
        expected,
        out=np.full_like(sample_dists, np.nan, dtype=float),
        where=expected > 0,
    )
    # Handle any bins where sample_dists > 0 but expected == 0 as +inf residual
    residuals[(sample_dists > 0) & (expected == 0)] = np.inf

    res_range = np.full(residuals.shape[1], np.nan, dtype=float)
    valid_cols = np.any(~np.isnan(residuals), axis=0)
    if np.any(valid_cols):
        valid_res = residuals[:, valid_cols]
        res_range[valid_cols] = np.nanmax(valid_res, axis=0) - np.nanmin(valid_res, axis=0)

    return res_range, residuals, ef_mean


def _solve_abundances_qp(norm_sub, sample_sub, reg_weight):
    """Solve the GuaCAMOLE QP for a subset of taxa and return normalized abundances."""
    n_taxa = norm_sub.shape[1]
    if n_taxa == 0:
        return np.array([], dtype=float)
    if n_taxa == 1:
        return np.array([1.0], dtype=float)

    m = get_matrix(norm_sub, reg_weight=reg_weight, sample_dist=sample_sub)
    p_scale = float(np.max(np.abs(m.data))) if m.nnz > 0 else 0.0
    if p_scale <= 0.0:
        return np.full(n_taxa, 1.0 / n_taxa, dtype=float)

    # Scale P by its max magnitude so CVXOPT tolerances remain scale-invariant
    p_scaled = m / p_scale
    a_eq = scipy.sparse.csc_matrix(np.ones((1, n_taxa), dtype=float))
    q_zero = np.zeros(n_taxa, dtype=float)
    b_eq = np.array([1.0], dtype=float)
    lb = np.zeros(n_taxa, dtype=float)

    x = qps.solve_qp(
        P=p_scaled, q=q_zero, A=a_eq, b=b_eq, solver='cvxopt', verbose=False,
        abstol=1e-15, reltol=1e-15, lb=lb,
    )
    if x is None:
        # Fallback for rank-deficient / exact-zero-objective PSD systems
        p_reg = p_scaled + scipy.sparse.eye(n_taxa, format='csc') * 1e-12
        x = qps.solve_qp(
            P=p_reg, q=q_zero, A=a_eq, b=b_eq, solver='cvxopt', verbose=False,
            abstol=1e-12, reltol=1e-12, lb=lb,
        )
    if x is None:
        return np.full(n_taxa, np.nan, dtype=float)

    inv_x = np.divide(1.0, x, out=np.full_like(x, np.nan, dtype=float), where=x > 0)
    total_inv = np.nansum(inv_x)
    if total_inv <= 0 or not np.isfinite(total_inv):
        return np.full(n_taxa, np.nan, dtype=float)
    return inv_x / total_inv


def corrected_abundances(sample_path, reference_path, quantiles=None, taxids=None, plot=False, reg_weight=0.1,
                         fp_cycles=3, verbose=False):
    """Estimate GC-bias-corrected species abundances and sequencing efficiencies.

    Parameters
    ----------
    sample_path : str, os.PathLike, or np.ndarray
        Path to sample distribution file or 2D array ``(n_bins, n_taxa)``.
    reference_path : str, os.PathLike, or np.ndarray
        Path to reference distribution file or 2D array ``(n_bins, n_taxa)``.
    quantiles : sequence of float, optional
        Quantile bounds for tail trimming.
    taxids : sequence of int, optional
        Taxonomy IDs corresponding to matrix columns.
    plot : bool
        Generate diagnostic PDF plots if True.
    reg_weight : float
        Regularization strength ``r`` in ``[0, 1]``.
    fp_cycles : int
        Number of false-positive removal threshold cycles.
    verbose : bool
        Write per-iteration JSON files if True.

    Returns
    -------
    tuple
        ``(all_ab, taxon_removal_cycle, efficiencies)``.
    """
    norm_dist, ref_dist, sample_dist = normalize_gc_dists(
        sample_path, reference_path, trim=True, quantiles=quantiles
    )

    n_bins, n_total = norm_dist.shape
    taxids = np.arange(n_total) if taxids is None else np.asarray(taxids)
    skipped_taxids = []
    cov_rel = np.zeros(norm_dist.shape, dtype=float)
    taxon_removal_cycle = np.full(len(taxids), np.nan, dtype=float)
    res_range = np.zeros(len(taxids), dtype=float)

    if fp_cycles <= 0:
        ind = np.arange(len(taxids))
        ab = _solve_abundances_qp(norm_dist, sample_dist, reg_weight)
        _, cov_rel, efficiencies = get_residuals(
            ref_dists=ref_dist,
            sample_dists=sample_dist,
            abundances=ab,
            plot=plot,
            iteration=1,
        )
        if plot:
            residual_plot(res_array=cov_rel, taxids=taxids, threshold_plot=None)
            plt.savefig('Residuals_final.pdf')
            plt.close()
            if os.path.exists('Efficiencies_it_1.pdf'):
                os.replace('Efficiencies_it_1.pdf', 'efficiencies_final.pdf')
        return ab, taxon_removal_cycle, efficiencies

    threshold = 10.0
    final_threshold = threshold / np.power(2.0, fp_cycles - 1)
    iteration = 1
    rep = 1
    ind_skip = []
    need_solve = True
    ab = np.array([], dtype=float)
    efficiencies = np.full(n_bins, np.nan, dtype=float)
    res_range_new = np.array([], dtype=float)
    last_eff_plot_it = None

    while threshold >= final_threshold or len(ind_skip) != 0:
        ind = np.where(np.isin(taxids, skipped_taxids, invert=True))[0]
        if len(ind) == 0:
            ab = np.array([], dtype=float)
            efficiencies = np.full(n_bins, np.nan, dtype=float)
            break

        # Only re-solve QP and recompute efficiencies/residuals when taxa changed or on initial round
        if need_solve:
            ab = _solve_abundances_qp(norm_dist[:, ind], sample_dist[:, ind], reg_weight)

            if plot:
                plt.figure(figsize=cm2inch(20, 13))
                plt.rcParams.update({'font.size': 16})
                plot_dist2(norm_dist[:, ind] / ab, line=True)
                plt.xlabel("GC bin (%)")
                plt.ylabel("Obs / (Exp*Abundance)")
                plt.savefig(f"Obs_vs_Exp_minimized_it_{iteration}_rep_{rep}.pdf")
                plt.close()

            res_range_new, cov_rel[:, ind], efficiencies = get_residuals(
                ref_dists=ref_dist[:, ind],
                sample_dists=sample_dist[:, ind],
                abundances=ab,
                plot=plot,
                iteration=iteration,
            )
            if plot:
                last_eff_plot_it = iteration
            res_range[ind] = res_range_new
            need_solve = False

        ind_skip = np.where(res_range_new > threshold)[0]

        if verbose:
            max_ef = float(np.nanmax(efficiencies)) if len(efficiencies) > 0 else 0.0
            scaled_eff = (efficiencies / max_ef).tolist() if max_ef > 0 else efficiencies.tolist()
            iteration_data = {
                'iteration': iteration,
                'rep': rep,
                'threshold': threshold,
                'taxids': taxids[ind].tolist(),
                'abundances': ab.tolist(),
                'efficiencies': scaled_eff,
                'residuals': cov_rel[:, ind].tolist(),
                'res_range': res_range_new.tolist(),
                'removed_taxids': taxids[ind][ind_skip].tolist() if len(ind_skip) > 0 else [],
            }
            with open(f'iteration_{iteration}_rep_{rep}.json', 'w') as f:
                json.dump(iteration_data, f)

        if len(ind_skip) > 0:
            skipped_taxids.extend(taxids[ind][ind_skip])
            taxon_removal_cycle[ind[ind_skip]] = iteration
            need_solve = True

        print(
            f"Removing taxa {taxids[ind][ind_skip]} in cycle {iteration} out of {fp_cycles}"
        )
        print(f"and with threshold {threshold}")

        if threshold == final_threshold and len(ind_skip) == 0:
            if plot:
                residual_plot(res_array=cov_rel[:, ind], taxids=taxids[ind], threshold_plot=None)
                plt.savefig('Residuals_final.pdf')
                plt.close()
                if last_eff_plot_it is not None:
                    eff_pdf = f'Efficiencies_it_{last_eff_plot_it}.pdf'
                    if os.path.exists(eff_pdf):
                        os.replace(eff_pdf, 'efficiencies_final.pdf')
            break

        rep += 1
        if threshold > final_threshold and len(ind_skip) == 0:
            iteration += 1
            threshold = threshold / 2.0
            rep = 1

    if verbose:
        manifest = {
            'total_iterations': iteration,
            'final_threshold': final_threshold,
            'all_skipped_taxids': [int(t) for t in skipped_taxids],
        }
        with open('iteration_manifest.json', 'w') as f:
            json.dump(manifest, f)

    all_ab = np.zeros(len(taxids), dtype=float)
    kept_mask = np.isin(taxids, skipped_taxids, invert=True)
    if len(ab) == np.sum(kept_mask):
        all_ab[kept_mask] = ab
    return all_ab, taxon_removal_cycle, efficiencies


def plot_dist2(dists, color=None, legend=False, line=False, **kwargs):
    if len(dists.shape) == 1:
        dists = dists.reshape((len(dists), 1))

    dists = dists[:, np.where(np.nansum(dists, axis=0) != 0)[0]]

    dist_df = pd.DataFrame(
        {
            'gc': np.tile(np.arange(0, dists.shape[0]), dists.shape[1]),
            'y': dists.flatten('F'),
            'org': [str(i) for i in np.repeat(np.arange(0, dists.shape[1]), dists.shape[0])],
        }
    )
    if color is None:
        if line:
            return sns.lineplot(
                x='gc', y='y', hue='org', data=dist_df.loc[dist_df['y'] != 0],
                legend=False, palette='deep', **kwargs,
            )
        return sns.scatterplot(
            x='gc', y='y', hue='org', style='org', data=dist_df.loc[dist_df['y'] != 0],
            legend=False, palette='deep', **kwargs,
        )

    dist_df['category'] = np.repeat(color, dists.shape[0])
    return sns.lineplot(
        x='gc', y='y', hue='category', style='org', data=dist_df.loc[dist_df['y'] != 0],
        legend=legend, **kwargs,
    )


def residual_plot(res_array, taxids, threshold_plot=None):
    gc_bin = []
    res_gc = []
    taxonomy = []
    for i in range(res_array.shape[1]):
        sub_array = res_array[:, i]
        valid = np.where(~np.isnan(sub_array))[0]
        gc_bin.extend(valid)
        res_gc.extend(sub_array[valid])
        taxonomy.extend(np.repeat(taxids[i], len(valid)))
    res_df = pd.DataFrame(
        {
            'gc': gc_bin,
            'residual': res_gc,
            'taxid': taxonomy,
        }
    )
    palette = sns.color_palette(['black'], len(res_df['taxid'].unique()))
    sns.lineplot(x='gc', y='residual', hue='taxid', data=res_df, legend=False, palette=palette, alpha=0.4)
    if threshold_plot is not None and len(res_df) > 0:
        plt.hlines(y=threshold_plot, xmin=min(res_df['gc']), xmax=max(res_df['gc']))
    plt.ylabel('Relative Residuals')
    plt.xlabel('GC bin [%]')


def _open_read_file(path):
    """Open plain text, gzip, or bzip2 read files in text mode."""
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    if path.endswith('.bz2'):
        return bz2.open(path, 'rt')
    return open(path)


def compute_gc(seq, nbin=100):
    """Compute dithered GC bin index for a DNA sequence (case-insensitive)."""
    gc = seq.count('G') + seq.count('C') + seq.count('g') + seq.count('c')
    atgc = seq.count('A') + seq.count('T') + seq.count('a') + seq.count('t') + gc
    if atgc == 0:
        return 0
    return floor((nbin / atgc) * (gc + random.random()))


def _compute_gc_pair(seq1, seq2, nbin=100):
    """Compute dithered GC bin index for a concatenated read pair without copying strings."""
    gc = (
        seq1.count('G') + seq1.count('C') + seq1.count('g') + seq1.count('c')
        + seq2.count('G') + seq2.count('C') + seq2.count('g') + seq2.count('c')
    )
    atgc = (
        seq1.count('A') + seq1.count('T') + seq1.count('a') + seq1.count('t')
        + seq2.count('A') + seq2.count('T') + seq2.count('a') + seq2.count('t')
        + gc
    )
    if atgc == 0:
        return 0
    return floor((nbin / atgc) * (gc + random.random()))


def create_sample_dist(fastq1, txids, kraken, fastq2=None, nbin=100, fasta=False):
    """Build per-taxon read GC count histogram (treating paired reads separately)."""
    occ = np.zeros((nbin + 1, len(txids)), dtype=int)
    idx = {str(t): col for col, t in enumerate(txids)}
    skip_after_seq = 0 if fasta else 2

    with ExitStack() as stack:
        f1 = stack.enter_context(_open_read_file(fastq1))
        f2 = stack.enter_context(_open_read_file(fastq2)) if fastq2 is not None else None
        kr = stack.enter_context(open(kraken))

        f1_readline = f1.readline
        f2_readline = f2.readline if f2 is not None else None
        idx_get = idx.get

        for ln_kr in kr:
            header1 = f1_readline()
            if not header1:
                break
            seq1 = f1_readline()
            if not seq1:
                break
            if skip_after_seq:
                f1_readline()
                f1_readline()

            if f2_readline is not None:
                f2_readline()
                seq2 = f2_readline()
                if skip_after_seq:
                    f2_readline()
                    f2_readline()

            parts = ln_kr.split("\t", 3)
            if len(parts) < 3:
                continue
            col = idx_get(parts[2])
            if col is None:
                continue

            gc_1 = compute_gc(seq1, nbin=nbin)
            if f2_readline is not None:
                gc_2 = compute_gc(seq2, nbin=nbin)
                if gc_1 > nbin or gc_2 > nbin:
                    continue
                occ[gc_2, col] += 1
            elif gc_1 > nbin:
                continue
            occ[gc_1, col] += 1

    header = "\t".join(map(str, txids))
    np.savetxt(f"sample_bin_{nbin}.dist", occ, fmt='%i', delimiter="\t", header=header)


def create_sample_dist_avg(fastq1, txids, kraken, fastq2, nbin=100, fasta=False):
    """Build per-taxon fragment GC count histogram by combining paired-end reads."""
    occ = np.zeros((nbin + 1, len(txids)), dtype=int)
    idx = {str(t): col for col, t in enumerate(txids)}
    skip_after_seq = 0 if fasta else 2

    with ExitStack() as stack:
        f1 = stack.enter_context(_open_read_file(fastq1))
        f2 = stack.enter_context(_open_read_file(fastq2))
        kr = stack.enter_context(open(kraken))

        f1_readline = f1.readline
        f2_readline = f2.readline
        idx_get = idx.get

        for ln_kr in kr:
            header1 = f1_readline()
            if not header1:
                break
            seq1 = f1_readline()
            if not seq1:
                break
            if skip_after_seq:
                f1_readline()
                f1_readline()

            f2_readline()
            seq2 = f2_readline()
            if not seq2:
                break
            if skip_after_seq:
                f2_readline()
                f2_readline()

            parts = ln_kr.split("\t", 3)
            if len(parts) < 3:
                continue
            col = idx_get(parts[2])
            if col is None:
                continue

            gc = _compute_gc_pair(seq1, seq2, nbin=nbin)
            if gc > nbin:
                continue
            occ[gc, col] += 1

    header = "\t".join(map(str, txids))
    np.savetxt(f"sample_bin_{nbin}.dist", occ, fmt='%i', delimiter="\t", header=header)
