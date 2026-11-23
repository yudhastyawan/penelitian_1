"""
visualize.py
============
Publication-Quality Diagnostic Visualization Suite for Joint Active-Passive
Picking-Free Full f-c Spectrum Inversion via Stein Mixture Variational Inference (SMVI).

Generates three comprehensive diagnostic figures for each benchmark scenario:
1. Waveform & Active-Passive MFJ Broadband Fusion Figure (`*_fig1_waveforms_and_fusion.png`):
   - Active MASW shot gather (with geophone coupling noise) & Passive ambient-noise
     empirical Green's function cross-correlation gather (with colored + azimuthal noise).
   - Individual Active MFJ and Passive MFJ spectrograms, variance-weighted overlap
     stacking curves across 8-12 Hz, and the fused broadband observed spectrum D_obs(f, c).
2. Observed vs. Posterior-Mean Synthetic f-c Spectrograms & Theoretical Modes
   (`*_fig2_fc_spectrum_comparison.png`):
   - Side-by-side comparison of Observed Fused D_obs(f, c), Posterior-Mean Synthetic
     D_syn(f, c | m_bar), and absolute spectral residual |D_obs - D_syn|.
   - Overlays theoretical normal-mode Rayleigh curves (M0, M1, M2 from `disba` / P-SV poles)
     and Leaky/Guided-P waveguide branches, plus SMVI convergence history.
3. Posterior Profiles (Vs, Vp, Poisson's Ratio) & Mixture Marginal Posteriors
   (`*_fig3_posterior_profiles_and_marginals.png`):
   - True vs. inverted 1D depth profiles for Vs(z), Vp(z), and Poisson's ratio nu(z)
     with 95% and 68% credible intervals and individual SMVI Gaussian mixture component means.
   - Mixture 1D & 2D marginal kernel density plots demonstrating uncertainty capture
     in the Low-Velocity Layer (Benchmark A) and Water-Table / Guided-P layer (Benchmark B).
"""

from __future__ import annotations
import os
from typing import Dict, List, Optional, Tuple
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

try:
    from scipy.stats import gaussian_kde
except Exception:

    class gaussian_kde:  # type: ignore[no-redef]
        """Self-contained 1D and 2D Gaussian Kernel Density Estimator fallback."""

        def __init__(self, dataset: np.ndarray, bw_method: float = 0.25) -> None:
            arr = np.asarray(dataset, dtype=np.float64)
            if arr.ndim == 1:
                arr = arr[None, :]
            self.dataset = arr  # shape (d, N)
            self.d, self.n = arr.shape
            stds = np.std(arr, axis=1, keepdims=True) + 1e-12
            self.bandwidth = stds * float(bw_method)

        def __call__(self, points: np.ndarray) -> np.ndarray:
            pts = np.asarray(points, dtype=np.float64)
            if pts.ndim == 1:
                pts = pts[None, :]  # (d, M)
            # Compute scaled squared distance: shape (d, M, N)
            diff = (pts[:, :, None] - self.dataset[:, None, :]) / self.bandwidth[
                :, None, :
            ]
            dist_sq = np.sum(diff**2, axis=0)  # (M, N)
            norm_const = (
                np.prod(self.bandwidth[:, 0])
                * ((2.0 * np.pi) ** (0.5 * self.d))
                * self.n
            )
            return np.sum(np.exp(-0.5 * dist_sq), axis=1) / norm_const


from models import LayeredEarthModel, InversionParameterSpace, evaluate_depth_profiles
from forward_mfj import FullWavefieldMFJOperator, compute_theoretical_dispersion_curves
from fusion import ActivePassiveSpectralFuser

# Publication styling defaults (Magnified for A4 PDF readability)
plt.rcParams.update(
    {
        "font.size": 10.2,
        "axes.titlesize": 11.0,
        "axes.labelsize": 10.2,
        "xtick.labelsize": 9.2,
        "ytick.labelsize": 9.2,
        "legend.fontsize": 8.4,
        "figure.titlesize": 13.0,
        "axes.grid": True,
        "grid.alpha": 0.28,
        "grid.linestyle": "--",
    }
)


def _overlay_dispersion_branches(
    ax: plt.Axes,
    curves: Dict[str, List[Tuple[np.ndarray, np.ndarray]]],
    show_legend: bool = True,
) -> None:
    """Overlay independent disba Rayleigh reference curves only."""
    ray_colors = ["#00ffff", "#76ff03", "#ea80fc"]
    for idx, (f_br, c_br) in enumerate(curves.get("rayleigh_modes", [])):
        col = ray_colors[idx % len(ray_colors)]
        lbl = f"Rayleigh $M_{idx}$ (Normal)" if show_legend else None
        ax.plot(
            f_br,
            c_br,
            color=col,
            linestyle="--",
            linewidth=1.6,
            label=lbl,
            alpha=0.92,
        )


def select_radial_traces(offsets, count=24):
    """Select traces near equally spaced radii, including both aperture ends."""
    r = np.asarray(offsets)
    if r.ndim != 1 or len(r) < 2 or np.any(np.diff(r) <= 0):
        raise ValueError("Trace offsets must be strictly increasing")
    targets = np.linspace(r[0], r[-1], min(count, len(r)))
    indices = np.unique(np.abs(r[:, None] - targets).argmin(axis=0))
    halfwidth = 0.35 * float(np.min(np.diff(r[indices])))
    return indices, halfwidth


def plot_waveforms_and_fusion(
    model: LayeredEarthModel,
    operator: FullWavefieldMFJOperator,
    fuser: ActivePassiveSpectralFuser,
    obs_fusion_dict: Dict[str, np.ndarray],
    t_act: np.ndarray,
    r_act: np.ndarray,
    seis_act: np.ndarray,
    t_pas: np.ndarray,
    r_pas: np.ndarray,
    seis_pas: np.ndarray,
    output_path: str,
) -> None:
    """
    Generate Figure 1: Synthetic Active MASW & Passive EGF Waveforms + Unified MFJ
    Active-Passive Broadband Fusion Diagnostics.
    """
    fig = plt.figure(figsize=(17.2, 10.6), dpi=240)
    gs = GridSpec(2, 3, figure=fig, hspace=0.33, wspace=0.25)

    # --- Panel (a): Active MASW Multi-Offset Gather (Opsi 4: 4m to 114m, Causal Moveout) ---
    ax_a = fig.add_subplot(gs[0, 0])
    max_t_act = min(float(t_act[-1]), max(1.2, float(r_act[-1]) / 100 + 0.15))
    t_mask_a = t_act <= max_t_act
    indices_a, width_a = select_radial_traces(r_act, 26)
    for i in indices_a:
        tr = seis_act[t_mask_a, i]
        tr_norm = tr / (np.max(np.abs(tr)) + 1e-30) * width_a
        ax_a.plot(
            t_act[t_mask_a], r_act[i] + tr_norm, color="#0d47a1", lw=1.0, alpha=0.92
        )
        ax_a.fill_between(
            t_act[t_mask_a],
            r_act[i],
            r_act[i] + tr_norm,
            where=(tr_norm > 0),
            color="#1e88e5",
            alpha=0.38,
        )
    ax_a.set_xlim(0.0, max_t_act)
    ax_a.set_ylim(0.0, r_act[-1] + 6.0)
    ax_a.set_xlabel("Causal Travel Time $t$ (s)")
    ax_a.set_ylabel("Source-Receiver Offset $r$ (m)")
    ax_a.set_title("(a) Active Gather (each trace normalized)", fontweight="bold")

    # --- Panel (b): Passive 2D Nested Array EGF Gather (Opsi 4: 10.9m to 562.6m, Causal Moveout) ---
    ax_b = fig.add_subplot(gs[0, 1])
    max_t_pas = min(float(t_pas[-1]), max(5.0, float(r_pas[-1]) / 100 + 0.2))
    t_mask_p = t_pas <= max_t_pas
    indices_p, width_p = select_radial_traces(r_pas, 24)
    for i in indices_p:
        tr = seis_pas[t_mask_p, i]
        tr_norm = tr / (np.max(np.abs(tr)) + 1e-30) * width_p
        ax_b.plot(
            t_pas[t_mask_p], r_pas[i] + tr_norm, color="#004d40", lw=1.05, alpha=0.92
        )
        ax_b.fill_between(
            t_pas[t_mask_p],
            r_pas[i],
            r_pas[i] + tr_norm,
            where=(tr_norm > 0),
            color="#00897b",
            alpha=0.38,
        )
    ax_b.set_xlim(0.0, max_t_pas)
    ax_b.set_ylim(0.0, r_pas[-1] + 28.0)
    ax_b.set_xlabel("Causal Lag Time $\\tau$ (s)")
    ax_b.set_ylabel("Inter-Station Distance $r_{ij}$ (m)")
    ax_b.set_title("(b) Passive EGFs (each trace normalized)", fontweight="bold")

    # --- Panel (c): Variance-Weighted Fusion Weights Across Overlap Band (8-12 Hz) ---
    ax_c = fig.add_subplot(gs[0, 2])
    f_bb = fuser.freqs_broadband
    ax_c.plot(
        f_bb,
        obs_fusion_dict["w_passive"],
        "o-",
        color="#00897b",
        lw=2.0,
        ms=3.5,
        label="Passive Weight $w_{\\mathrm{pas}}(f)$",
    )
    ax_c.plot(
        f_bb,
        obs_fusion_dict["w_active"],
        "s-",
        color="#3949ab",
        lw=2.0,
        ms=3.5,
        label="Active Weight $w_{\\mathrm{act}}(f)$",
    )
    ax_c.axvspan(
        fuser.config.f_overlap_low,
        fuser.config.f_overlap_high,
        color="#ffd54f",
        alpha=0.32,
        label="8-12 Hz Overlap Transition Band",
    )
    ax_c.set_xlim(f_bb[0], f_bb[-1])
    ax_c.set_ylim(-0.05, 1.08)
    ax_c.set_xlabel("Frequency $f$ (Hz)")
    ax_c.set_ylabel("Shape-variance Fusion Weight")
    ax_c.set_title("(c) Shape-Based Active-Passive Fusion Weights", fontweight="bold")
    ax_c.legend(loc="center right", framealpha=0.9)

    c_grid = fuser.c_grid
    # --- Panel (d): Passive MFJ Spectrogram (0.5 - 12 Hz) ---
    ax_d = fig.add_subplot(gs[1, 0])
    im_d = ax_d.pcolormesh(
        operator.freqs_passive,
        c_grid,
        obs_fusion_dict["spec_passive_raw"].T,
        shading="auto",
        cmap="magma",
        vmin=0.0,
        vmax=1.0,
    )
    ax_d.axvspan(8.0, 12.0, color="cyan", alpha=0.14, label="8-12 Hz Overlap")
    ax_d.set_xlabel("Frequency $f$ (Hz)")
    ax_d.set_ylabel("Phase Velocity $c$ (m/s)")
    ax_d.set_title("(d) Passive MFJ $f$-$c$ Spectrum (0.5-12 Hz)", fontweight="bold")
    plt.colorbar(im_d, ax=ax_d, pad=0.02, label="Normalized MFJ Amplitude")

    # --- Panel (e): Active MFJ Spectrogram (8 - 55 Hz) ---
    ax_e = fig.add_subplot(gs[1, 1])
    im_e = ax_e.pcolormesh(
        operator.freqs_active,
        c_grid,
        obs_fusion_dict["spec_active_raw"].T,
        shading="auto",
        cmap="magma",
        vmin=0.0,
        vmax=1.0,
    )
    ax_e.axvspan(8.0, 12.0, color="cyan", alpha=0.14, label="8-12 Hz Overlap")
    alias_c = 2 * operator.active_geom.dx * operator.freqs_active
    ax_e.plot(operator.freqs_active, alias_c, "w--", lw=1.0, label="c = 2 dx f")
    ax_e.fill_between(
        operator.freqs_active,
        c_grid[0],
        np.maximum(c_grid[0], alias_c),
        color="white",
        alpha=0.16,
    )
    ax_e.legend(loc="upper right", fontsize=7)
    ax_e.set_ylim(c_grid[0], c_grid[-1])

    ax_e.set_xlabel("Frequency $f$ (Hz)")
    ax_e.set_ylabel("Phase Velocity $c$ (m/s)")
    ax_e.set_title("(e) Active MASW MFJ $f$-$c$ Spectrum (8-55 Hz)", fontweight="bold")
    plt.colorbar(im_e, ax=ax_e, pad=0.02, label="Normalized MFJ Amplitude")

    # --- Panel (f): Fused Broadband Observed Spectrum D_obs(f, c) with Theoretical Modes ---
    ax_f = fig.add_subplot(gs[1, 2])
    im_f = ax_f.pcolormesh(
        f_bb,
        c_grid,
        obs_fusion_dict["d_fused"].T,
        shading="auto",
        cmap="magma",
        vmin=0.0,
        vmax=1.0,
    )
    curves = compute_theoretical_dispersion_curves(model, f_bb, c_grid)
    _overlay_dispersion_branches(ax_f, curves, show_legend=True)
    ax_f.axvline(8.0, color="white", ls=":", lw=1.1, alpha=0.8)
    ax_f.axvline(12.0, color="white", ls=":", lw=1.1, alpha=0.8)
    ax_f.set_xlabel("Frequency $f$ (Hz)")
    ax_f.set_ylabel("Phase Velocity $c$ (m/s)")
    ax_f.set_title(
        "(f) Fused Broadband Spectrum $D_{\\mathrm{obs}}(f, c)$ (1-55 Hz)",
        fontweight="bold",
    )
    ax_f.legend(
        loc="upper right", framealpha=0.85, facecolor="#121212", labelcolor="white"
    )
    plt.colorbar(im_f, ax=ax_f, pad=0.02, label="Fused Spectral Energy")

    fig.suptitle(
        f"{model.name}\nUnified Active-Passive Modified Frequency-Bessel (MFJ) Processing & Broadband Fusion",
        fontweight="bold",
        y=0.98,
    )
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def plot_fc_spectrum_comparison(
    true_model: LayeredEarthModel,
    fuser: ActivePassiveSpectralFuser,
    d_obs: np.ndarray,
    smvi_results: Dict[str, np.ndarray],
    output_path: str,
) -> None:
    """
    Generate Figure 2: Observed vs. Posterior-Mean Synthetic Fused f-c Spectrograms
    with independent Rayleigh reference curves and optimization history.
    """
    fig = plt.figure(figsize=(15.5, 9.0), dpi=220)
    gs = GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.24)

    f_bb = fuser.freqs_broadband
    c_grid = fuser.c_grid
    predictive = "d_posterior_predictive_mean" in smvi_results
    d_syn = smvi_results.get(
        "d_posterior_predictive_mean", smvi_results["d_syn_posterior_mean"]
    )
    curves_true = compute_theoretical_dispersion_curves(true_model, f_bb, c_grid)

    # --- Panel (a): Observed Fused Broadband f-c Spectrogram ---
    ax_a = fig.add_subplot(gs[0, 0])
    im_a = ax_a.pcolormesh(
        f_bb, c_grid, d_obs.T, shading="auto", cmap="turbo", vmin=0.0, vmax=1.0
    )
    _overlay_dispersion_branches(ax_a, curves_true, show_legend=True)
    ax_a.axvspan(8.0, 12.0, color="white", alpha=0.10)
    ax_a.set_xlabel("Frequency $f$ (Hz)")
    ax_a.set_ylabel("Phase Velocity $c$ (m/s)")
    ax_a.set_title(
        "(a) Observed Broadband Fused Spectrum $D_{\\mathrm{obs}}(f, c)$",
        fontweight="bold",
    )
    ax_a.legend(
        loc="upper right", framealpha=0.85, facecolor="#1a1a1a", labelcolor="white"
    )
    plt.colorbar(im_a, ax=ax_a, pad=0.02, label="Normalized Spectral Energy")

    # --- Panel (b): Posterior-Mean Synthetic Fused Broadband f-c Spectrogram ---
    ax_b = fig.add_subplot(gs[0, 1])
    im_b = ax_b.pcolormesh(
        f_bb, c_grid, d_syn.T, shading="auto", cmap="turbo", vmin=0.0, vmax=1.0
    )
    _overlay_dispersion_branches(ax_b, curves_true, show_legend=True)
    ax_b.axvspan(8.0, 12.0, color="white", alpha=0.10)
    ax_b.set_xlabel("Frequency $f$ (Hz)")
    ax_b.set_ylabel("Phase Velocity $c$ (m/s)")
    ax_b.set_title(
        (
            "(b) Mean predicted spectrum across mixture samples"
            if predictive
            else "(b) Spectrum evaluated at the mean physical model"
        ),
        fontweight="bold",
    )
    ax_b.legend(
        loc="upper right", framealpha=0.85, facecolor="#1a1a1a", labelcolor="white"
    )
    plt.colorbar(im_b, ax=ax_b, pad=0.02, label="Normalized Spectral Energy")

    # --- Panel (c): Spectral Residual Map |D_obs - D_syn| ---
    ax_c = fig.add_subplot(gs[1, 0])
    residual = np.abs(d_obs - d_syn)
    im_c = ax_c.pcolormesh(
        f_bb, c_grid, residual.T, shading="auto", cmap="inferno", vmin=0.0, vmax=0.6
    )
    _overlay_dispersion_branches(ax_c, curves_true, show_legend=False)
    ax_c.set_xlabel("Frequency $f$ (Hz)")
    ax_c.set_ylabel("Phase Velocity $c$ (m/s)")
    mean_res = float(np.mean(residual))
    ax_c.set_title(
        f"(c) Absolute Spectral Misfit Residual $|D_{{\\mathrm{{obs}}}} - D_{{\\mathrm{{syn}}}}|$ (Mean = {mean_res:.3f})",
        fontweight="bold",
    )
    plt.colorbar(im_c, ax=ax_c, pad=0.02, label="Absolute Spectral Residual")

    # --- Panel (d): SMVI Convergence Diagnostics (NCC, Wasserstein W2, and ELBO) ---
    ax_d = fig.add_subplot(gs[1, 1])
    iters = np.arange(1, len(smvi_results["history_misfit_min"]) + 1)
    ax_d.plot(
        iters,
        smvi_results["history_misfit_mean"],
        "--",
        color="#5c6bc0",
        lw=1.6,
        label="Mixture MC-average misfit",
    )
    ax_d.plot(
        iters,
        smvi_results["history_misfit_min"],
        "-",
        color="#d32f2f",
        lw=2.2,
        label="Minimum component MC-average misfit",
    )
    ax_d.plot(
        iters,
        smvi_results["history_ncc"],
        ":",
        color="#00796b",
        lw=1.8,
        label="Frequency-Slice $1 - \\mathrm{NCC}$ Misfit",
    )
    ax_d.set_xlabel("SMVI Iteration")
    ax_d.set_ylabel("Picking-Free Spectral Misfit")
    ax_d.set_title(
        "(d) SMVI Optimization History\n(Convergence Unconfirmed)", fontweight="bold"
    )

    ax_d2 = ax_d.twinx()
    ax_d2.plot(
        iters,
        smvi_results["history_w2"],
        "-.",
        color="#ef6c00",
        lw=1.7,
        label="1D Wasserstein $W_2$ RMS (m/s)",
    )
    ax_d2.set_ylabel("1D Wasserstein $W_2$ Distance (m/s)", color="#ef6c00")
    ax_d2.tick_params(axis="y", labelcolor="#ef6c00")
    ax_d2.grid(False)

    lines1, labels1 = ax_d.get_legend_handles_labels()
    lines2, labels2 = ax_d2.get_legend_handles_labels()
    ax_d.legend(lines1 + lines2, labels1 + labels2, loc="upper right", framealpha=0.9)

    fig.suptitle(
        f"{true_model.name}\nPicking-Free Numerical P-SV Spectrum Inversion Fit",
        fontweight="bold",
        y=0.98,
    )
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def plot_posterior_profiles_and_marginals(
    true_model: LayeredEarthModel,
    param_space: InversionParameterSpace,
    smvi_results: Dict[str, np.ndarray],
    output_path: str,
) -> None:
    """
    Generate Figure 3: True vs. Inverted Vs, Vp, and Poisson's Ratio Depth Profiles
    with 95% Credible Intervals + Mixture Marginal Posterior Density Plots.
    """
    fig = plt.figure(figsize=(16.0, 9.6), dpi=220)
    gs = GridSpec(2, 3, figure=fig, hspace=0.34, wspace=0.28)

    z_max = 32.0
    z_grid = np.linspace(0.0, z_max, 320)
    samples_phys = smvi_results["posterior_samples_phys"]
    prof_stats = evaluate_depth_profiles(samples_phys, param_space, z_grid)
    true_prof = evaluate_depth_profiles(
        param_space.true_vector[None, :], param_space, z_grid
    )
    comp_means_phys = smvi_results["comp_means_phys"]
    comp_prof = evaluate_depth_profiles(comp_means_phys, param_space, z_grid)
    weights = smvi_results["weights"]

    # --- Panels (a, b, c): 1D Depth Profiles for Vs(z), Vp(z), and Poisson's Ratio nu(z) ---
    panel_configs = [
        ("vs", "$V_s$ (m/s)", "(a) Shear Velocity $V_s(z)$ Profile", "#1565c0"),
        ("vp", "$V_p$ (m/s)", "(b) Compressional Velocity $V_p(z)$ Profile", "#00695c"),
        (
            "nu",
            "Poisson's Ratio $\\nu$",
            "(c) Dynamic Poisson's Ratio $\\nu(z)$ Profile",
            "#6a1b9a",
        ),
    ]

    for idx, (key, xlabel, title, color) in enumerate(panel_configs):
        ax = fig.add_subplot(gs[0, idx])
        # 95% and 68% credible intervals
        ax.fill_betweenx(
            z_grid,
            prof_stats[f"{key}_p025"],
            prof_stats[f"{key}_p975"],
            color=color,
            alpha=0.18,
            label="95% mixture quantiles (uncalibrated)",
        )
        ax.fill_betweenx(
            z_grid,
            prof_stats[f"{key}_p16"],
            prof_stats[f"{key}_p84"],
            color=color,
            alpha=0.32,
            label="68% mixture quantiles",
        )
        # Plot individual SMVI Gaussian Mixture component means
        for k in range(len(weights)):
            lbl = "SMVI Mixture Components" if k == 0 else None
            ax.plot(
                comp_prof[f"{key}_samples"][k],
                z_grid,
                color="#ff8f00",
                ls=":",
                lw=1.2,
                alpha=0.75,
                label=lbl,
            )
        # Posterior mean & True profile
        ax.plot(
            prof_stats[f"{key}_mean"],
            z_grid,
            color=color,
            lw=2.3,
            label="SMVI Posterior Mean",
        )
        ax.plot(
            true_prof[f"{key}_mean"],
            z_grid,
            color="#d32f2f",
            ls="--",
            lw=2.2,
            label="True Benchmark Profile",
        )
        ax.set_ylim(z_max, 0.0)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Depth $z$ (m)")
        ax.set_title(title, fontweight="bold")
        ax.legend(loc="lower left", framealpha=0.88, fontsize=7.8)

    # --- Panels (d, e, f): Mixture Marginal & Joint Posterior Densities ---
    # Target Layer 2 (the Low-Velocity Layer in Benchmark A, or the Shallow Water-Table
    # Saturated Aquifer Layer in Benchmark B)
    is_bench_a = true_model.code == "benchmark_a"
    vs2_samples = samples_phys[:, 1]
    vp2_samples = samples_phys[:, 5]
    h1_samples = samples_phys[:, 8]
    h2_samples = samples_phys[:, 9]
    r2 = (vp2_samples / vs2_samples) ** 2
    nu2_samples = (0.5 * r2 - 1.0) / (r2 - 1.0)

    # Panel (d): Marginal KDE of Layer 2 Vs (LVL or Saturated Layer) showing individual mixture components
    ax_d = fig.add_subplot(gs[1, 0])
    x_vs2 = np.linspace(param_space.vs_min[1], param_space.vs_max[1], 250)
    kde_vs2 = gaussian_kde(vs2_samples, bw_method=0.22)
    y_vs2 = kde_vs2(x_vs2)
    ax_d.fill_between(x_vs2, 0, y_vs2, color="#1565c0", alpha=0.25)
    ax_d.plot(x_vs2, y_vs2, color="#0d47a1", lw=2.2, label="SMVI Mixture Posterior KDE")

    # Overlay individual Gaussian mixture component densities in physical space
    for k in range(len(weights)):
        comp_mask = smvi_results["posterior_comp_labels"] == k
        if np.sum(comp_mask) > 8:
            try:
                kde_k = gaussian_kde(vs2_samples[comp_mask], bw_method=0.35)
                lbl_k = "Weighted component KDE" if k == 0 else None
                ax_d.plot(
                    x_vs2,
                    weights[k] * kde_k(x_vs2),
                    "--",
                    color="#ef6c00",
                    lw=1.3,
                    alpha=0.85,
                    label=lbl_k,
                )
            except Exception:
                pass

    ax_d.axvline(
        true_model.vs[1],
        color="#d32f2f",
        ls="--",
        lw=2.2,
        label=f"True $V_{{s2}} = {true_model.vs[1]:.0f}$ m/s",
    )
    ax_d.axvline(
        np.mean(vs2_samples),
        color="#0d47a1",
        ls=":",
        lw=2.0,
        label=f"Posterior Mean = {np.mean(vs2_samples):.1f} m/s",
    )
    ax_d.set_xlabel("Layer 2 Shear Velocity $V_{s2}$ (m/s)")
    ax_d.set_ylabel("Marginal Posterior Density")
    title_d = (
        "(d) Mixture LVL $V_{s2}$ Posterior Density"
        if is_bench_a
        else "(d) Saturated Layer $V_{s2}$ Posterior Density"
    )
    ax_d.set_title(title_d, fontweight="bold")
    ax_d.legend(loc="upper right", framealpha=0.9)

    # Panel (e): Marginal KDE of Layer 2 Vp / Poisson's Ratio (Water-Table Guided-P sensitivity)
    ax_e = fig.add_subplot(gs[1, 1])
    x_vp2 = np.linspace(param_space.vp_min[1], param_space.vp_max[1], 250)
    kde_vp2 = gaussian_kde(vp2_samples, bw_method=0.22)
    y_vp2 = kde_vp2(x_vp2)
    ax_e.fill_between(x_vp2, 0, y_vp2, color="#00695c", alpha=0.25)
    ax_e.plot(x_vp2, y_vp2, color="#004d40", lw=2.2, label="SMVI Mixture Posterior KDE")

    for k in range(len(weights)):
        comp_mask = smvi_results["posterior_comp_labels"] == k
        if np.sum(comp_mask) > 8:
            try:
                kde_k = gaussian_kde(vp2_samples[comp_mask], bw_method=0.35)
                lbl_k = "Weighted component KDE" if k == 0 else None
                ax_e.plot(
                    x_vp2,
                    weights[k] * kde_k(x_vp2),
                    "--",
                    color="#ef6c00",
                    lw=1.3,
                    alpha=0.85,
                    label=lbl_k,
                )
            except Exception:
                pass

    ax_e.axvline(
        true_model.vp[1],
        color="#d32f2f",
        ls="--",
        lw=2.2,
        label=f"True $V_{{p2}} = {true_model.vp[1]:.0f}$ m/s",
    )
    ax_e.axvline(
        np.mean(vp2_samples),
        color="#004d40",
        ls=":",
        lw=2.0,
        label=f"Posterior Mean = {np.mean(vp2_samples):.1f} m/s",
    )
    ax_e.set_xlabel("Layer 2 Compressional Velocity $V_{p2}$ (m/s)")
    ax_e.set_ylabel("Marginal Posterior Density")
    title_e = (
        f"(e) LVL $V_{{p2}}$ Posterior ($\\bar{{\\nu}}_2 = {np.mean(nu2_samples):.3f}$)"
        if is_bench_a
        else f"(e) Water-Table $V_{{p2}}$ Posterior ($\\bar{{\\nu}}_2 = {np.mean(nu2_samples):.3f}$)"
    )
    ax_e.set_title(title_e, fontweight="bold")
    ax_e.legend(loc="upper right", framealpha=0.9)

    # Panel (f): 2D Joint Posterior Density Contour showing Mixture Trade-Off
    ax_f = fig.add_subplot(gs[1, 2])
    if is_bench_a:
        # LVL Velocity Vs_2 vs LVL Thickness h_2 trade-off
        x_var, y_var = vs2_samples, h2_samples
        x_true, y_true = true_model.vs[1], true_model.thickness[1]
        x_comp, y_comp = comp_means_phys[:, 1], comp_means_phys[:, 9]
        x_lbl = "LVL Shear Velocity $V_{s2}$ (m/s)"
        y_lbl = "LVL Thickness $h_2$ (m)"
        title_f = "(f) Joint Mixture Posterior: LVL $V_{s2}$ vs. $h_2$"
    else:
        # Water-Table Depth h_1 vs Saturated Aquifer Vp_2 trade-off
        x_var, y_var = vp2_samples, h1_samples
        x_true, y_true = true_model.vp[1], true_model.thickness[0]
        x_comp, y_comp = comp_means_phys[:, 5], comp_means_phys[:, 8]
        x_lbl = "Saturated Aquifer $V_{p2}$ (m/s)"
        y_lbl = "Water-Table Depth $h_1$ (m)"
        title_f = "(f) Joint Mixture Posterior: Water-Table $V_{p2}$ vs. $h_1$"

    xg = np.linspace(np.min(x_var) * 0.94, np.max(x_var) * 1.06, 90)
    yg = np.linspace(np.min(y_var) * 0.92, np.max(y_var) * 1.08, 90)
    XX, YY = np.meshgrid(xg, yg)
    kde_2d = gaussian_kde(np.vstack([x_var, y_var]), bw_method=0.26)
    ZZ = kde_2d(np.vstack([XX.ravel(), YY.ravel()])).reshape(XX.shape)

    ax_f.contourf(XX, YY, ZZ, levels=14, cmap="YlGnBu", alpha=0.85)
    ax_f.contour(XX, YY, ZZ, levels=6, colors="#0d47a1", linewidths=0.7, alpha=0.6)
    ax_f.scatter(x_var, y_var, s=8, c="#263238", alpha=0.25, label="Posterior Samples")
    ax_f.scatter(
        x_comp,
        y_comp,
        s=75,
        c="#ff8f00",
        edgecolors="black",
        marker="^",
        zorder=5,
        label="SMVI Mixture Centers $\\boldsymbol{\\mu}_k$",
    )
    ax_f.scatter(
        [x_true],
        [y_true],
        s=130,
        c="#d32f2f",
        edgecolors="white",
        marker="*",
        zorder=6,
        label="True Benchmark Parameter",
    )
    ax_f.set_xlabel(x_lbl)
    ax_f.set_ylabel(y_lbl)
    ax_f.set_title(title_f, fontweight="bold")
    ax_f.legend(loc="upper right", framealpha=0.9, fontsize=7.5)

    fig.suptitle(
        f"{true_model.name}\nSMVI mixture profiles and samples (uncalibrated)",
        fontweight="bold",
        y=0.98,
    )
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
