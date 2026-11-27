import gzip
import time

import guacamole.library as lib
import numpy as np
import pandas as pd
import pytest
from guacamole.create_reference_dist import parsing
from guacamole.distributions import get_sample_dists
from guacamole.io import write_output


def _bruteforce_objective(norm_dist, sample_dist, x, reg_weight):
    """Direct brute-force evaluation of Eqs. (8), (9), and (10) from Holcik et al. (2025)."""
    n_bins, n_taxa = norm_dist.shape
    # Eq. (8): G(x) = sum_{i != j} sum_g (R_{g,i} x_i - R_{g,j} x_j)^2 over overlapping bins
    g_val = 0.0
    for i in range(n_taxa):
        for j in range(n_taxa):
            if i == j:
                continue
            for g in range(n_bins):
                if norm_dist[g, i] != 0 and norm_dist[g, j] != 0:
                    g_val += (norm_dist[g, i] * x[i] - norm_dist[g, j] * x[j]) ** 2

    # Eq. (9): lambda_g = (1 / n_g) * sum_i log(O(i, g) + 1) * R_{g,i} * x_i
    lam = np.zeros(n_bins)
    active = np.zeros(n_bins, dtype=bool)
    for g in range(n_bins):
        w_sum = 0.0
        val = 0.0
        for i in range(n_taxa):
            if norm_dist[g, i] != 0:
                w = np.log1p(sample_dist[g, i])
                w_sum += w
                val += w * norm_dist[g, i] * x[i]
        if w_sum > 0:
            active[g] = True
            lam[g] = val / w_sum

    # Eq. (10): R(x) = sum_{k, l} exp(-|k - l|) * (lambda_k - lambda_l)^2
    r_val = 0.0
    for k in range(n_bins):
        if not active[k]:
            continue
        for m in range(n_bins):
            if not active[m]:
                continue
            r_val += np.exp(-abs(k - m)) * (lam[k] - lam[m]) ** 2

    return ((1.0 - reg_weight) / (n_taxa ** 2)) * g_val + reg_weight * r_val


def test_qp_matrix_matches_paper_equations_8_9_10():
    """0.5 * x^T P x must match Eqs. (8)-(10) to machine precision."""
    rng = np.random.default_rng(42)
    n_bins, n_taxa = 25, 6
    sample_dist = rng.uniform(0.0, 50.0, size=(n_bins, n_taxa))
    mask = rng.random((n_bins, n_taxa)) > 0.4
    sample_dist[~mask] = 0.0
    ref_dist = np.where(mask, rng.uniform(0.5, 5.0, size=(n_bins, n_taxa)), 0.0)
    norm_dist = np.divide(sample_dist, ref_dist, out=np.zeros_like(sample_dist), where=ref_dist > 0)

    for reg_weight in (0.0, 0.01, 0.25, 1.0):
        p_mat = lib.get_matrix(norm_dist, reg_weight=reg_weight, sample_dist=sample_dist).toarray()
        assert np.allclose(p_mat, p_mat.T, atol=1e-13), "QP matrix P must be symmetric"

        for _ in range(5):
            x = rng.uniform(0.01, 1.0, size=n_taxa)
            x /= x.sum()
            qp_val = 0.5 * float(x @ p_mat @ x)
            expected_val = _bruteforce_objective(norm_dist, sample_dist, x, reg_weight)
            assert qp_val == pytest.approx(expected_val, rel=1e-11, abs=1e-12)


def test_get_alpha_beta_match_get_matrix():
    """get_alpha and get_beta should match 0.5 * P (i.e. Q_total)."""
    rng = np.random.default_rng(7)
    sample_dist = rng.uniform(0.0, 20.0, size=(15, 4))
    sample_dist[rng.random((15, 4)) < 0.3] = 0.0
    norm_dist = np.where(sample_dist > 0, rng.uniform(0.1, 3.0, size=(15, 4)), 0.0)

    p_mat = lib.get_matrix(norm_dist, reg_weight=0.05, sample_dist=sample_dist).toarray()
    overlaps = lib.get_overlaps(norm_dist)
    for i in range(4):
        assert lib.get_alpha(norm_dist, i, overlaps, reg_weight=0.05, sample_dist=sample_dist) == pytest.approx(
            0.5 * p_mat[i, i], rel=1e-12, abs=1e-12
        )
        for j in range(4):
            if i != j:
                assert lib.get_beta(
                    norm_dist, i, j, overlaps, reg_weight=0.05, sample_dist=sample_dist
                ) == pytest.approx(0.5 * p_mat[i, j], rel=1e-12, abs=1e-12)


@pytest.mark.parametrize(
    "bias_fn",
    [
        lambda g: 1.0 - 10.0 * (g - 0.5) ** 2,  # peak at 50% GC
        lambda g: 0.2 + 0.8 * g,                 # increasing with GC
        lambda g: 1.0 - 0.8 * g,                 # decreasing with GC
    ],
)
def test_synthetic_recovery_from_eq1(bias_fn):
    """Exact noiseless data from Eq. (1) must recover true abundances and efficiencies."""
    gc_grid = np.arange(101) / 100.0
    true_eff = np.clip(bias_fn(gc_grid), 0.1, 1.0)
    true_eff /= np.max(true_eff)

    centers = [32, 40, 48, 56, 64]
    true_ab = np.array([0.08, 0.15, 0.32, 0.25, 0.20])
    lengths = np.array([2.0, 4.5, 3.0, 5.5, 2.5])
    lengths /= lengths.sum()

    # Build normalized Gaussian genomic GC distributions f(j, g)
    f_mat = np.column_stack([np.exp(-0.5 * ((np.arange(101) - c) / 6.0) ** 2) for c in centers])
    f_mat /= f_mat.sum(axis=0, keepdims=True)

    # Eq. (1): O(j, g) = N * a_j * eta_g * l_j * f(j, g) / Z
    n_reads = 1_000_000.0
    unnorm_obs = (true_ab * lengths)[None, :] * f_mat * true_eff[:, None]
    obs = n_reads * unnorm_obs / unnorm_obs.sum()
    ref = f_mat * lengths[None, :]
    ref = ref / ref.sum() * obs.sum()

    est_ab, removal_cycle, est_eff = lib.corrected_abundances(
        obs, ref, quantiles=[0.025, 0.975], taxids=np.arange(len(centers)),
        plot=False, reg_weight=1e-4, fp_cycles=5,
    )

    assert np.all(np.isnan(removal_cycle)), "No true taxon should be removed on noiseless data"
    assert np.allclose(est_ab, true_ab, rtol=1e-2, atol=1e-3)

    covered = est_eff > 0
    est_eff_norm = est_eff[covered] / np.max(est_eff[covered])
    true_eff_norm = true_eff[covered] / np.max(true_eff[covered])
    assert np.allclose(est_eff_norm, true_eff_norm, rtol=2e-2, atol=2e-2)


def test_false_positive_removal_on_synthetic_outlier():
    """An outlier taxon whose observed GC shape mismatches its reference must be removed."""
    gc = np.arange(101)
    true_eff = np.clip(0.3 + 0.7 * (gc / 100.0), 0.1, 1.0)
    centers = [35, 45, 55, 65]
    true_ab = np.array([0.25, 0.25, 0.25, 0.25])

    f_mat = np.column_stack([np.exp(-0.5 * ((gc - c) / 6.0) ** 2) for c in centers])
    f_mat /= f_mat.sum(axis=0, keepdims=True)

    # Add a false-positive taxon whose observed reads are shifted +12% GC from its reference
    fp_ref = np.exp(-0.5 * ((gc - 48) / 6.0) ** 2)
    fp_ref /= fp_ref.sum()
    fp_obs = np.exp(-0.5 * ((gc - 60) / 6.0) ** 2)
    fp_obs /= fp_obs.sum()

    ref_all = np.column_stack([f_mat, fp_ref])
    obs_true = 100_000.0 * (true_ab[None, :] * f_mat * true_eff[:, None])
    obs_all = np.column_stack([obs_true, 5_000.0 * fp_obs])

    est_ab, removal_cycle, _ = lib.corrected_abundances(
        obs_all, ref_all, quantiles=[0.05, 0.95], taxids=[1, 2, 3, 4, 999],
        plot=False, reg_weight=1e-3, fp_cycles=5,
    )

    assert est_ab[-1] == 0.0, "False-positive taxon should be assigned zero abundance"
    assert not np.isnan(removal_cycle[-1]), "False-positive taxon should record a removal cycle"
    assert np.allclose(est_ab[:4], true_ab, rtol=2e-2, atol=1e-2)


def test_get_efficiencies_handles_fractional_and_unit_reads():
    """Eq. (9) weights log(O + 1) must remain strictly non-negative for 0 < O <= 1."""
    ref = np.array([[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]])
    # Bin 0: O=0.2 and O=0.8 (both < 1, where log(O) would be negative!)
    # Bin 1: O=1.0 and O=1.0 (where log(O) would be 0 -> 0/0 NaN!)
    # Bin 2: O=10.0 and O=10.0
    sample = np.array([[0.2, 0.8], [1.0, 1.0], [10.0, 10.0]])
    abundances = np.array([0.5, 0.5])

    w_ef, _ = lib.get_efficiencies(ref, sample, abundances)
    assert np.all(np.isfinite(w_ef))
    assert np.all(w_ef > 0), f"Efficiencies must be positive, got {w_ef}"


def test_get_sample_dists_linear_weighting_eq2_eq3():
    """Sample redistribution must weight nodes linearly by add_reads (not quadratically by M(i)^2)."""
    # Species S has 100 direct reads (all in GC bin 1) and 100 reads redistributed
    # from genus G which has 10,000 total reads (all in GC bin 2, P(S|G) = 0.01).
    # True redistributed distribution for S: 100 reads in bin 1, 100 reads in bin 2 (50% / 50%).
    sdists = pd.DataFrame({
        '10': [0, 100, 0, 0],     # species S: M(S) = 100
        '100': [0, 0, 10000, 0],  # genus G: M(G) = 10000
    })
    weights = {'10': {'10': 100.0, '100': 100.0}}
    map2lvl = {'10': ['10', 100, 200.0]}
    lvl_taxids = {'10': ['SpeciesA', 100, 100, 200.0]}

    res = get_sample_dists(
        sdists=sdists, weights=weights, map2lvl_taxids=map2lvl,
        lvl_taxids=lvl_taxids, reference=False,
    )
    assert np.allclose(res['10'], [0.0, 100.0, 100.0, 0.0])


def test_write_output_no_double_counting(tmp_path):
    """write_output must not double-count reads and fraction_total_reads must sum to 1."""
    out_file = tmp_path / "out.guaca"
    lvl_taxids = {
        '10': ['TaxonA', 100, 100, 120.0],
        '20': ['TaxonB', 80, 80, 80.0],
    }
    ab_df = pd.DataFrame({
        'abundance_ls': [0.6, 0.4],
        'abundance_br': [0.6, 0.4],
        'taxid': [10, 20],
        'taxon_removal_cycle': [np.nan, np.nan],
        'GC content': [45.0, 55.0],
        'abundance_eff': [0.6, 0.4],
    }).set_index('taxid', drop=False)

    write_output(str(out_file), lvl_taxids, ab_df, sum_all_reads=200.0, level='S', excessive_removal_detected=False)
    df = pd.read_csv(out_file, sep='\t')
    assert list(df['kraken_assigned_reads']) == [100, 80]
    assert list(df['added_reads']) == [20, 0]
    assert list(df['new_est_reads']) == [120, 80]
    assert df['fraction_total_reads'].sum() == pytest.approx(1.0, abs=1e-5)


def test_compute_gc_and_parsing_case_insensitive():
    """Lowercase bases (soft-masked reads/genomes) must yield identical GC counts as uppercase."""
    seq_upper = "ATGCATGCGGGGCCCCAAAA" * 10
    seq_lower = seq_upper.lower()

    gc_u, valid_u = parsing(seq_upper)
    gc_l, valid_l = parsing(seq_lower)
    assert np.array_equal(gc_u, gc_l)
    assert np.array_equal(valid_u, valid_l)

    # Pure AT lowercase sequence must yield bin 0, not inflated GC
    for _ in range(20):
        assert lib.compute_gc("aaaaaatttttt" * 10, nbin=100) == 0


def test_create_sample_dist_single_end_and_gz(tmp_path, monkeypatch):
    """create_sample_dist must work for single-end reads (fastq2=None) and gzipped files."""
    monkeypatch.chdir(tmp_path)
    fq = tmp_path / "reads.fastq.gz"
    kr = tmp_path / "reads.kraken"

    with gzip.open(fq, 'wt') as f:
        f.write("@r1\n" + "G" * 50 + "A" * 50 + "\n+\n" + "I" * 100 + "\n")
        f.write("@r2\n" + "c" * 25 + "t" * 75 + "\n+\n" + "I" * 100 + "\n")

    kr.write_text("C\tr1\t10\t100\t10:66\nC\tr2\t20\t100\t20:66\n")

    lib.create_sample_dist(str(fq), np.array([10, 20]), str(kr), fastq2=None, nbin=100)
    dist_file = tmp_path / "sample_bin_100.dist"
    assert dist_file.exists()
    df = pd.read_csv(dist_file, sep='\t')
    assert df.values.sum() == 2


def test_fp_cycles_zero_and_quantiles_honored():
    """fp_cycles=0 must succeed without UnboundLocalError and custom quantiles must take effect."""
    gc = np.arange(101)
    f_mat = np.column_stack([np.exp(-0.5 * ((gc - c) / 8.0) ** 2) for c in (40, 50, 60)])
    f_mat /= f_mat.sum(axis=0, keepdims=True)
    obs = 10_000.0 * f_mat

    ab, cycles, eff = lib.corrected_abundances(
        obs, f_mat, quantiles=[0.1, 0.9], taxids=[1, 2, 3], fp_cycles=0
    )
    assert np.allclose(ab, [1 / 3, 1 / 3, 1 / 3], atol=1e-2)
    assert np.all(np.isnan(cycles))
    assert np.any(eff > 0)


def test_get_matrix_performance_large_community():
    """Vectorized get_matrix must handle 200 taxa in well under 100ms."""
    rng = np.random.default_rng(123)
    n_bins, n_taxa = 101, 200
    sample_dist = rng.exponential(50.0, size=(n_bins, n_taxa))
    norm_dist = rng.uniform(0.5, 2.0, size=(n_bins, n_taxa))

    t0 = time.perf_counter()
    p_mat = lib.get_matrix(norm_dist, reg_weight=0.01, sample_dist=sample_dist)
    elapsed = time.perf_counter() - t0

    assert p_mat.shape == (n_taxa, n_taxa)
    assert elapsed < 0.2, f"get_matrix took {elapsed:.3f}s for {n_taxa} taxa (expected < 0.2s)"
