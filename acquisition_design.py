"""Measured geometry diagnostics and actual transform comparisons. No prescribed winner."""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Tuple
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


@dataclass
class FieldArrayDesign:
    """Represents a 2D passive station geometry + active MASW spread."""

    name: str
    short_label: str
    passive_xy: np.ndarray  # shape (N_stations, 2) in meters
    active_x: np.ndarray  # shape (N_channels,) in meters
    shot_offsets: np.ndarray
    color: str
    is_recommended: bool = False

    @property
    def n_stations(self) -> int:
        return self.passive_xy.shape[0]

    def pairwise_distances_and_azimuths(self) -> Tuple[np.ndarray, np.ndarray]:
        """Compute all N*(N-1)/2 inter-station radial distances r_ij and azimuths theta_ij."""
        n = self.n_stations
        dists = []
        azims = []
        for i in range(n):
            for j in range(i + 1, n):
                dx = self.passive_xy[j, 0] - self.passive_xy[i, 0]
                dy = self.passive_xy[j, 1] - self.passive_xy[i, 1]
                r = np.hypot(dx, dy)
                if r > 0.5:
                    dists.append(r)
                    azims.append(np.arctan2(dy, dx))
        return np.array(dists, dtype=np.float64), np.array(azims, dtype=np.float64)


def build_four_candidate_field_designs() -> List[FieldArrayDesign]:
    """Construct the 4 candidate field acquisition designs (each using 25 passive sensors + 48 active geophones)."""
    active_x = np.arange(48, dtype=np.float64) * 2.0  # 0 to 94 m

    # --- Option 1: Classical 1D Collinear Array (25 stations @ 25 m spacing = 600 m) ---
    xy_opt1 = np.zeros((25, 2), dtype=np.float64)
    xy_opt1[:, 0] = np.arange(25, dtype=np.float64) * 25.0

    # --- Option 2: Orthogonal 2D L-Array (13 stations along X @ 25 m, 12 along Y @ 25 m) ---
    xy_opt2 = np.zeros((25, 2), dtype=np.float64)
    xy_opt2[0:13, 0] = np.arange(13, dtype=np.float64) * 25.0
    xy_opt2[13:25, 1] = np.arange(1, 13, dtype=np.float64) * 25.0

    # --- Option 3: Concentric Circular / Nested Triangle SPAC Rings (1 center + 3 rings of 8) ---
    xy_opt3 = np.zeros((25, 2), dtype=np.float64)
    radii = [25.0, 95.0, 280.0]
    idx = 1
    for r_ring in radii:
        angles = np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False)
        for ang in angles:
            xy_opt3[idx, 0] = r_ring * np.cos(ang)
            xy_opt3[idx, 1] = r_ring * np.sin(ang)
            idx += 1

    # Option 4: 9 inner stations and 16 logarithmic spiral stations.
    xy_opt4 = np.zeros((25, 2), dtype=np.float64)
    core_x = np.linspace(4.0, 88.0, 9)
    core_y = np.array([0.0, 4.0, -4.0, 7.0, -7.0, 5.0, -5.0, 3.0, 0.0])
    xy_opt4[0:9, 0] = core_x
    xy_opt4[0:9, 1] = core_y

    golden_angle = np.pi * (3.0 - np.sqrt(5.0))  # 137.50776 deg
    radii_spiral = np.geomspace(52.0, 320.0, 16)
    for k in range(16):
        ang = k * golden_angle + 0.35
        r_s = radii_spiral[k]
        xy_opt4[9 + k, 0] = 45.0 + r_s * np.cos(ang)
        xy_opt4[9 + k, 1] = r_s * np.sin(ang)

    return [
        FieldArrayDesign(
            name="Option 1: 1D Collinear Array (25 stn @ 25 m + 1-Shot MASW)",
            short_label="Opt 1: 1D Collinear",
            passive_xy=xy_opt1,
            active_x=active_x,
            shot_offsets=np.array([4.0]),
            color="#d32f2f",
            is_recommended=False,
        ),
        FieldArrayDesign(
            name="Option 2: 2D Orthogonal L-Array (300x300 m + Dual-Shot MASW)",
            short_label="Opt 2: 2D L-Array",
            passive_xy=xy_opt2,
            active_x=active_x,
            shot_offsets=np.array([4.0, 10.0]),
            color="#f57c00",
            is_recommended=False,
        ),
        FieldArrayDesign(
            name="Option 3: Concentric Circular Rings (SPAC/MAM + MASW)",
            short_label="Opt 3: Circular SPAC",
            passive_xy=xy_opt3,
            active_x=active_x,
            shot_offsets=np.array([4.0, 10.0]),
            color="#7b1fa2",
            is_recommended=False,
        ),
        FieldArrayDesign(
            name="Option 4: Hybrid Nested Log-Linear 2D + Multi-Offset Spine",
            short_label="Opt 4: Hybrid Nested 2D",
            passive_xy=xy_opt4,
            active_x=active_x,
            shot_offsets=np.array([-20.0, -10.0, -4.0, 4.0, 10.0, 20.0]),
            color="#00796b",
            is_recommended=False,
        ),
    ]


def evaluate_acquisition_metrics(designs):
    results = {}
    for d in designs:
        distances, angles = d.pairwise_distances_and_azimuths()
        counts, _ = np.histogram(distances, bins=np.linspace(5.0, 585.0, 30))
        results[d.short_label] = {
            "station_count": d.n_stations,
            "pair_count": len(distances),
            "raw_min_distance_m": float(distances.min()),
            "max_distance_m": float(distances.max()),
            "populated_radial_bins_pct": float(100 * np.mean(counts > 0)),
            "second_azimuthal_moment": float(abs(np.mean(np.exp(2j * angles)))),
            "minimum_spacing_nyquist_proxy_hz_at_190mps": float(
                190 / (2 * distances.min())
            ),
            "nyquist_proxy_caveat": "Minimum-spacing proxy only; does not establish alias-free bandwidth of an irregular array.",
        }
    return results


def make_operator(design, freqs_act, freqs_pas, c_grid, quadrature_density=1.0):
    from forward_mfj import (
        ActiveSurveyGeometry,
        PassiveSurveyGeometry,
        FullWavefieldMFJOperator,
    )

    active = ActiveSurveyGeometry(near_offsets=tuple(design.shot_offsets))
    passive = PassiveSurveyGeometry(
        n_stations=design.n_stations,
        coordinates=design.passive_xy,
        use_option4_2d_nested=False,
        n_radial_bins=0,
    )
    return FullWavefieldMFJOperator(
        active, passive, freqs_act, freqs_pas, c_grid, quadrature_density
    )


def plot_acquisition_design_comparison(output_path):
    designs = build_four_candidate_field_designs()
    metrics = evaluate_acquisition_metrics(designs)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout="constrained")
    for d in designs:
        r, angles = d.pairwise_distances_and_azimuths()
        axes[0, 0].scatter(*d.passive_xy.T, s=12, color=d.color, label=d.short_label)
        counts, edges = np.histogram(r, bins=np.linspace(0, 600, 31))
        axes[0, 1].plot(
            (edges[1:] + edges[:-1]) / 2, counts, color=d.color, label=d.short_label
        )
        # True normalized 2D plane-wave array response along kx, ky=0.
        k = np.linspace(-0.5, 0.5, 500)
        arf = abs(np.exp(1j * k[:, None] * d.passive_xy[None, :, 0]).mean(axis=1)) ** 2
        axes[1, 0].plot(k, arf, color=d.color, label=d.short_label)
        az = np.linspace(0, np.pi, 181)
        axes[1, 1].plot(
            np.rad2deg(az),
            abs(np.cos(2 * (angles[:, None] - az)).mean(axis=0)),
            color=d.color,
        )
    titles = [
        "Station coordinates",
        "Pair-distance histogram",
        "Array response along kx (ky=0)",
        "Directional second moment (not velocity bias)",
    ]
    labels = [
        ("X (m)", "Y (m)"),
        ("Distance (m)", "Number of pairs"),
        ("kx (rad/m)", "Normalized power"),
        ("Azimuth (degrees)", "Absolute second moment"),
    ]
    for ax, title, (x, y) in zip(axes.flat, titles, labels):
        ax.set(title=title, xlabel=x, ylabel=y)
        ax.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=7)
    axes[0, 0].set_aspect("equal")
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
    return metrics


def transform_spectrum(operator, u, mode, kind="mfj"):
    from forward_mfj import _mfj_hankel_transform_numba
    from scipy.special import j0

    if kind == "mfj":
        return operator.compute_mfj_spectrum(u, mode, use_ccfj_if_available=False)
    active = mode == "active"
    f = operator.freqs_active if active else operator.freqs_passive
    r = operator.r_act if active else operator.r_pas
    dr = operator.dr_act if active else operator.dr_pas
    win = operator.win_act if active else operator.win_pas
    kr = (
        2 * np.pi * f[:, None, None] * r[None, None, :] / operator.c_grid[None, :, None]
    )
    kernel = np.exp(-1j * kr) if kind == "plane_wave" else j0(kr).astype(complex)
    return _mfj_hankel_transform_numba(u, r, dr, win, np.ascontiguousarray(kernel))


def plot_dispersion_image_acquisition_comparison(
    model, c_grid, output_path, freqs_act=None, freqs_pas=None, quadrature_density=1.0
):
    from fusion import ActivePassiveSpectralFuser

    freqs_act = np.linspace(8, 56, 33) if freqs_act is None else freqs_act
    freqs_pas = np.linspace(0.5, 12, 24) if freqs_pas is None else freqs_pas
    designs = build_four_candidate_field_designs()
    panels = []
    report = {}
    for i, design in enumerate(designs):
        op = make_operator(design, freqs_act, freqs_pas, c_grid, quadrature_density)
        fuser = ActivePassiveSpectralFuser(op)
        ua, _ = op.synthesize_frequency_gather(
            model, "active", True, np.random.default_rng(2029)
        )
        up, _ = op.synthesize_frequency_gather(
            model, "passive", True, np.random.default_rng(2030)
        )
        # Reference is noise-free data from the same design, never the observed image itself.
        clean_a, _ = op.synthesize_frequency_gather(model, "active", False)
        clean_p, _ = op.synthesize_frequency_gather(model, "passive", False)
        kinds = ["plane_wave", "standing_bessel", "mfj"] if i == 0 else ["mfj"]
        for kind in kinds:
            ref = fuser.fuse_spectra(
                transform_spectrum(op, clean_a, "active", kind),
                transform_spectrum(op, clean_p, "passive", kind),
            )["d_fused"]
            spec = fuser.fuse_spectra(
                transform_spectrum(op, ua, "active", kind),
                transform_spectrum(op, up, "passive", kind),
            )["d_fused"]
            title = design.short_label + " / " + kind
            ov = (fuser.freqs_broadband >= 8) & (fuser.freqs_broadband <= 12)
            ncc = float(
                np.sum(spec[ov] * ref[ov])
                / (np.linalg.norm(spec[ov]) * np.linalg.norm(ref[ov]) + 1e-30)
            )
            report[title] = {
                "overlap_ncc_against_noise_free_same_geometry": ncc,
                "spectral_mae_against_noise_free_same_geometry": float(
                    np.mean(abs(spec - ref))
                ),
                "reference": f"noise-free {kind} of the same geometry; noise consistency only, not method ranking",
            }
            panels.append((title, spec, fuser.freqs_broadband))
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), layout="constrained")
    for ax, (title, spec, freqs) in zip(axes.flat, panels):
        im = ax.pcolormesh(
            freqs, c_grid, spec.T, shading="auto", vmin=0, vmax=1, cmap="turbo"
        )
        ax.set(title=title, xlabel="Frequency (Hz)", ylabel="Phase velocity (m/s)")
        ax.title.set_fontsize(9)
        fig.colorbar(im, ax=ax)
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
    return report
