"""Publication-quality diagnostic plotting for GuaCAMOLE."""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


def _apply_clean_style():
    """Return a clean matplotlib rcParams dict for publication-ready figures."""
    return {
        'font.family': 'sans-serif',
        'font.sans-serif': ['DejaVu Sans', 'Helvetica', 'Arial'],
        'font.size': 9.5,
        'axes.titlesize': 11,
        'axes.titleweight': 'bold',
        'axes.labelsize': 10,
        'axes.edgecolor': '#333333',
        'axes.linewidth': 0.8,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.labelsize': 8.5,
        'ytick.labelsize': 8.5,
        'xtick.color': '#333333',
        'ytick.color': '#333333',
        'grid.color': '#E0E0E0',
        'grid.linewidth': 0.6,
        'legend.frameon': False,
        'legend.fontsize': 8.5,
    }


def _gc_color_map(gc_values):
    """Map per-taxon genomic GC values (0-100) to perceptually uniform colors."""
    cmap = mpl.colormaps['viridis']
    norm = mpl.colors.Normalize(vmin=25.0, vmax=75.0, clip=True)
    return [cmap(norm(v)) for v in gc_values], cmap, norm


def _plot_curve_lines(ax, x_bins, matrix, colors, alpha=0.65, linewidth=1.1, linestyle='-'):
    """Plot nonzero contiguous segments of each column of ``matrix``."""
    for j in range(matrix.shape[1]):
        y = matrix[:, j]
        valid = np.isfinite(y) & (y != 0)
        if np.sum(valid) < 2:
            continue
        ax.plot(
            x_bins[valid],
            y[valid],
            color=colors[j],
            alpha=alpha,
            linewidth=linewidth,
            linestyle=linestyle,
        )


def _compute_gc_xlim(matrix, pad=2):
    """Compute sensible GC x-axis bounds covering active bins."""
    active = np.any(np.isfinite(matrix) & (matrix != 0), axis=1)
    idx = np.flatnonzero(active)
    if len(idx) == 0:
        return 20, 80
    return max(0, int(idx[0]) - pad), min(matrix.shape[0] - 1, int(idx[-1]) + pad)


def _draw_summary_page(
    diagnostics,
    ab_df,
    efficiencies,
    excessive_removal_detected,
    sample_label,
):
    """Build the 2x2 summary diagnostic figure."""
    norm_dist = diagnostics['norm_dist']
    final_ind = diagnostics['final_ind']
    final_ab = diagnostics['final_ab']
    final_residuals = diagnostics['final_residuals']
    final_threshold = diagnostics['final_threshold']
    taxids = diagnostics['taxids']

    n_bins, n_taxa = norm_dist.shape
    x_bins = np.arange(n_bins)
    kept_mask = np.zeros(n_taxa, dtype=bool)
    kept_mask[final_ind] = True
    removed_ind = np.flatnonzero(~kept_mask)
    n_kept = len(final_ind)
    n_removed = len(removed_ind)

    # Per-taxon GC content (aligned with diagnostics['taxids'])
    gc_by_taxid = ab_df['GC content'].to_dict()
    taxon_gc = np.array([gc_by_taxid.get(int(t), 50.0) for t in taxids], dtype=float)
    colors_all, cmap, gc_norm = _gc_color_map(taxon_gc)
    xmin, xmax = _compute_gc_xlim(norm_dist)

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2), constrained_layout=True)
    ax_a, ax_b, ax_c, ax_d = axes.flat

    # ---- Panel A: Raw Obs/Exp ratios before abundance scaling ----
    if n_removed > 0:
        _plot_curve_lines(
            ax_a, x_bins, norm_dist[:, removed_ind],
            colors=['#D62728'] * n_removed, alpha=0.35, linewidth=0.9, linestyle='--',
        )
    if n_kept > 0:
        kept_colors = [colors_all[i] for i in final_ind]
        _plot_curve_lines(
            ax_a, x_bins, norm_dist[:, final_ind],
            colors=kept_colors, alpha=0.7, linewidth=1.1,
        )
    ax_a.set_title('A. Raw Observed / Expected Ratios (Before Scaling)', loc='left')
    ax_a.set_xlabel('GC content (%)')
    ax_a.set_ylabel('Observed / Expected reads')
    ax_a.set_xlim(xmin, xmax)
    ax_a.grid(True, axis='both', alpha=0.5)

    # Colorbar for genomic GC on Panel A
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=gc_norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax_a, fraction=0.04, pad=0.02)
    cbar.set_label('Taxon GC (%)', fontsize=8.5)
    cbar.ax.tick_params(labelsize=7.5)

    # ---- Panel B: Aligned taxon curves + inferred sequencing efficiency ----
    eff_clean = np.asarray(efficiencies, dtype=float)
    eff_valid = np.isfinite(eff_clean) & (eff_clean > 0)
    eff_max = float(np.nanmax(eff_clean)) if np.any(eff_valid) else np.nan
    eff_norm = (eff_clean / eff_max) if (np.isfinite(eff_max) and eff_max > 0) else np.full_like(eff_clean, np.nan)

    if n_kept > 0 and np.any(eff_valid) and not excessive_removal_detected:
        # Scale taxon curves onto the same [0, 1] relative efficiency scale
        scaled = norm_dist[:, final_ind] / final_ab[None, :]
        # Match global scale to eff_norm via median ratio over active bins
        pos_mask = (scaled > 0) & eff_valid[:, None]
        if np.any(pos_mask):
            ratios = eff_norm[:, None] / np.where(pos_mask, scaled, np.nan)
            scale_factor = float(np.nanmedian(ratios))
            if np.isfinite(scale_factor) and scale_factor > 0:
                scaled_aligned = scaled * scale_factor
                kept_colors = [colors_all[i] for i in final_ind]
                _plot_curve_lines(
                    ax_b, x_bins, scaled_aligned,
                    colors=kept_colors, alpha=0.32, linewidth=0.95,
                )

    # Subtle background bars indicating taxon coverage per GC bin
    if n_kept > 0:
        coverage = (norm_dist[:, final_ind] > 0).sum(axis=1)
        max_cov = int(np.max(coverage))
        if max_cov > 0:
            ax_b.bar(
                x_bins,
                0.15 * (coverage / max_cov),
                width=0.85,
                color='#B0BEC5',
                alpha=0.35,
                zorder=1,
                label=f'Bin support (max {max_cov} taxa)',
            )

    if np.any(eff_valid) and not excessive_removal_detected:
        ax_b.plot(
            x_bins[eff_valid],
            eff_norm[eff_valid],
            color='#111111',
            linewidth=2.2,
            marker='o',
            markersize=3.2,
            zorder=5,
            label='Inferred efficiency',
        )
        peak_gc = int(x_bins[eff_valid][np.argmax(eff_norm[eff_valid])])
        ax_b.annotate(
            f'Peak: {peak_gc}% GC',
            xy=(peak_gc, 1.0),
            xytext=(min(xmax - 12, peak_gc + 6), 1.04),
            fontsize=8.5,
            fontweight='bold',
            color='#111111',
            arrowprops={'arrowstyle': '->', 'color': '#555555', 'lw': 0.8},
        )
    elif excessive_removal_detected:
        ax_b.text(
            0.5, 0.5,
            'Excessive taxon removal (>50%)\nEfficiency unreliable (NaN)',
            transform=ax_b.transAxes, ha='center', va='center',
            fontsize=10, color='#B71C1C', fontweight='bold',
        )

    ax_b.set_title('B. Aligned Taxa & Sequencing Efficiency', loc='left')
    ax_b.set_xlabel('GC content (%)')
    ax_b.set_ylabel('Relative sequencing efficiency')
    ax_b.set_xlim(xmin, xmax)
    ax_b.set_ylim(0, 1.14)
    ax_b.grid(True, axis='both', alpha=0.5)
    ax_b.legend(loc='lower right')

    # ---- Panel C: Relative residuals ----
    ax_c.axhline(0.0, color='#222222', linestyle='--', linewidth=0.9, zorder=1)
    if n_removed > 0:
        for idx_r in removed_ind:
            y = final_residuals[:, idx_r]
            valid = np.isfinite(y)
            if np.sum(valid) >= 2:
                ax_c.plot(
                    x_bins[valid], y[valid],
                    color='#D62728', alpha=0.35, linewidth=0.9, linestyle='--',
                )
        ax_c.plot([], [], color='#D62728', alpha=0.65, linewidth=1.2, linestyle='--',
                  label=f'Removed as FP ({n_removed})')

    if n_kept > 0:
        for idx_k in final_ind:
            y = final_residuals[:, idx_k]
            valid = np.isfinite(y)
            if np.sum(valid) >= 2:
                ax_c.plot(
                    x_bins[valid], y[valid],
                    color='#37474F', alpha=0.42, linewidth=0.95,
                )
        ax_c.plot([], [], color='#37474F', alpha=0.75, linewidth=1.2,
                  label=f'Kept taxa ({n_kept})')

    # Compute sensible y-limits focused on kept taxa (don't let extreme FP outliers flatten kept curves)
    if n_kept > 0:
        kept_res = final_residuals[:, final_ind]
        finite_kept = kept_res[np.isfinite(kept_res)]
        if len(finite_kept) > 0:
            y_abs = max(float(np.max(np.abs(finite_kept))) * 1.25, 0.25)
            ax_c.set_ylim(-y_abs, y_abs)

    if final_threshold is not None and np.isfinite(final_threshold):
        ax_c.text(
            0.02, 0.96,
            f'Final range threshold: {final_threshold:g}',
            transform=ax_c.transAxes, ha='left', va='top',
            fontsize=8.2, color='#444444',
            bbox={'boxstyle': 'round,pad=0.25', 'facecolor': '#F5F5F5', 'edgecolor': '#CCCCCC'},
        )

    ax_c.set_title('C. Model Fit: Relative Residuals', loc='left')
    ax_c.set_xlabel('GC content (%)')
    ax_c.set_ylabel('Relative residual  (Obs - Exp) / Exp')
    ax_c.set_xlim(xmin, xmax)
    ax_c.grid(True, axis='both', alpha=0.5)
    ax_c.legend(loc='upper right')

    # ---- Panel D: Abundance correction (GuaCAMOLE / Bracken) vs GC content ----
    ax_d.axhline(1.0, color='#222222', linestyle='--', linewidth=0.9, zorder=1)
    if not excessive_removal_detected:
        valid_df = ab_df.dropna(subset=['GC content', 'abundance_br']).copy()
        valid_df = valid_df[valid_df['abundance_br'] > 0]

        kept_df = valid_df[valid_df['abundance_ls'] > 1e-9]
        removed_df = valid_df[(valid_df['abundance_ls'] <= 1e-9) & valid_df['abundance_eff'].notna()]

        if len(kept_df) > 0:
            ratio_ls = kept_df['abundance_ls'].to_numpy() / kept_df['abundance_br'].to_numpy()
            sizes = 25.0 + 450.0 * np.sqrt(np.clip(kept_df['abundance_ls'].to_numpy(), 0, 1))
            sc = ax_d.scatter(
                kept_df['GC content'].to_numpy(),
                ratio_ls,
                s=sizes,
                c=kept_df['GC content'].to_numpy(),
                cmap=cmap,
                norm=gc_norm,
                edgecolors='#222222',
                linewidths=0.6,
                alpha=0.88,
                zorder=4,
                label='GuaCAMOLE / Bracken (kept)',
            )
            _ = sc

        if len(removed_df) > 0:
            ratio_eff = removed_df['abundance_eff'].to_numpy() / removed_df['abundance_br'].to_numpy()
            ax_d.scatter(
                removed_df['GC content'].to_numpy(),
                ratio_eff,
                s=45,
                facecolors='none',
                edgecolors='#D62728',
                linewidths=1.1,
                zorder=3,
                label='Eff.-corrected / Bracken (FP)',
            )

        all_ratios = []
        if len(kept_df) > 0:
            all_ratios.extend((kept_df['abundance_ls'] / kept_df['abundance_br']).dropna().tolist())
        if len(removed_df) > 0:
            all_ratios.extend((removed_df['abundance_eff'] / removed_df['abundance_br']).dropna().tolist())
        all_ratios = [r for r in all_ratios if np.isfinite(r) and r > 0]
        if all_ratios:
            ax_d.set_yscale('log', base=2)
            lo = min(min(all_ratios), 0.85) * 0.88
            hi = max(max(all_ratios), 1.15) * 1.14
            ax_d.set_ylim(lo, hi)
            candidate_ticks = [0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 8.0, 16.0]
            ticks = [t for t in candidate_ticks if lo <= t <= hi]
            if len(ticks) < 3:
                ticks = np.geomspace(lo * 1.05, hi * 0.95, 5).tolist()
            ax_d.set_yticks(ticks)
            ax_d.yaxis.set_minor_locator(mpl.ticker.NullLocator())
            ax_d.yaxis.set_major_formatter(
                mpl.ticker.FuncFormatter(lambda v, _pos: f'{v:.2g}x')
            )
    else:
        ax_d.text(
            0.5, 0.5,
            'Excessive taxon removal (>50%)\nCorrection skipped',
            transform=ax_d.transAxes, ha='center', va='center',
            fontsize=10, color='#B71C1C', fontweight='bold',
        )

    ax_d.set_title('D. Abundance Correction vs. Genomic GC', loc='left')
    ax_d.set_xlabel('Genomic GC content (%)')
    ax_d.set_ylabel('Fold change (GuaCAMOLE / Bracken)')
    ax_d.grid(True, axis='both', alpha=0.5)
    if not excessive_removal_detected:
        ax_d.legend(loc='best')

    # Overall figure header
    status_str = (
        f'{n_kept}/{n_taxa} taxa kept'
        if n_removed == 0
        else f'{n_kept}/{n_taxa} taxa kept, {n_removed} removed as FP'
    )
    if excessive_removal_detected:
        status_str += '  [WARNING: >50% Bracken abundance removed]'
    fig.suptitle(f'GuaCAMOLE Diagnostic Summary — {sample_label}  ({status_str})',
                 fontsize=12.5, fontweight='bold')
    return fig


def _draw_iteration_pages(diagnostics, sample_label):
    """Build optional compact pages showing false-positive removal rounds."""
    history = diagnostics.get('history', [])
    # Only include rounds where an optimization occurred and taxa were removed (plus final state)
    interesting = [h for h in history if len(h['removed_idx']) > 0]
    if not interesting:
        return []

    n_bins = diagnostics['norm_dist'].shape[0]
    x_bins = np.arange(n_bins)
    xmin, xmax = _compute_gc_xlim(diagnostics['norm_dist'])

    pages = []
    per_page = 4
    for start in range(0, len(interesting), per_page):
        chunk = interesting[start:start + per_page]
        nrows = len(chunk)
        fig, axes = plt.subplots(
            nrows, 2, figsize=(11.0, 2.7 * nrows + 0.8),
            squeeze=False, constrained_layout=True,
        )

        for row_i, step in enumerate(chunk):
            ax_eff, ax_res = axes[row_i]
            ind = step['ind']
            removed_local = set(step['removed_idx'].tolist())
            eff = step['efficiencies']
            residuals = step['residuals']
            thresh = step['threshold']
            it = step['iteration']
            rep = step['rep']

            # Left: efficiency curve at this step
            valid_ef = np.isfinite(eff) & (eff > 0)
            if np.any(valid_ef):
                ef_norm = eff / np.nanmax(eff)
                ax_eff.plot(x_bins[valid_ef], ef_norm[valid_ef], color='#1565C0', lw=1.8, marker='o', ms=2.5)
            ax_eff.set_title(
                f'Cycle {it} (rep {rep}, T={thresh:g}) — Efficiency ({len(ind)} taxa)',
                loc='left', fontsize=9.5,
            )
            ax_eff.set_ylabel('Rel. efficiency')
            ax_eff.set_xlim(xmin, xmax)
            ax_eff.set_ylim(0, 1.12)
            ax_eff.grid(True, alpha=0.45)

            # Right: residuals highlighting removed outliers in red
            ax_res.axhline(0.0, color='#222222', ls='--', lw=0.8)
            for local_j in range(len(ind)):
                y = residuals[:, local_j]
                valid = np.isfinite(y)
                if np.sum(valid) < 2:
                    continue
                if local_j in removed_local:
                    ax_res.plot(x_bins[valid], y[valid], color='#D62728', alpha=0.8, lw=1.2, zorder=4)
                else:
                    ax_res.plot(x_bins[valid], y[valid], color='#78909C', alpha=0.35, lw=0.8, zorder=2)

            ax_res.set_title(
                f'Residuals — {len(removed_local)} outlier(s) removed (red)',
                loc='left', fontsize=9.5,
            )
            ax_res.set_ylabel('Rel. residual')
            ax_res.set_xlim(xmin, xmax)
            ax_res.grid(True, alpha=0.45)

            if row_i == nrows - 1:
                ax_eff.set_xlabel('GC content (%)')
                ax_res.set_xlabel('GC content (%)')

        fig.suptitle(
            f'GuaCAMOLE False-Positive Removal Progression — {sample_label}',
            fontsize=12, fontweight='bold',
        )
        pages.append(fig)

    return pages


def save_diagnostic_plots(
    output_path,
    diagnostics,
    ab_df,
    efficiencies,
    excessive_removal_detected=False,
    include_iterations=False,
):
    """Generate a single clean diagnostic PDF named ``<output_stem>_plots.pdf``.

    Parameters
    ----------
    output_path : str
        Path to the main GuaCAMOLE output file (used to derive the PDF path and title).
    diagnostics : dict
        Diagnostic dictionary returned by ``corrected_abundances(..., return_diagnostics=True)``.
    ab_df : pd.DataFrame
        Results DataFrame with abundance estimates and GC content.
    efficiencies : np.ndarray
        Final normalized sequencing efficiencies.
    excessive_removal_detected : bool
        Whether GuaCAMOLE flagged excessive taxon removal.
    include_iterations : bool
        If True (e.g. when ``--verbose`` is set), append pages showing each
        false-positive removal round to the same multi-page PDF.

    Returns
    -------
    str
        Path to the saved PDF file.
    """
    stem, _ = os.path.splitext(output_path)
    pdf_path = f'{stem}_plots.pdf'
    sample_label = os.path.basename(stem)

    with mpl.rc_context(_apply_clean_style()), PdfPages(pdf_path) as pdf:
        summary_fig = _draw_summary_page(
            diagnostics=diagnostics,
            ab_df=ab_df,
            efficiencies=efficiencies,
            excessive_removal_detected=excessive_removal_detected,
            sample_label=sample_label,
        )
        pdf.savefig(summary_fig)
        plt.close(summary_fig)

        if include_iterations:
            for page_fig in _draw_iteration_pages(diagnostics, sample_label):
                pdf.savefig(page_fig)
                plt.close(page_fig)

    return pdf_path
