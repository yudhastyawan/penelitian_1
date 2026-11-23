"""
fusion.py
=========
Active-Passive MFJ Spectral Normalization and Variance-Weighted Broadband Fusion Module.

Native-frequency slice normalization and heuristic spectral-shape weighting.
No interpolation adds new data frequencies. The variance-named arrays contain
shape proxies, not calibrated noise variances or optimal statistical weights.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Optional, Tuple
import numpy as np

try:
    from numba import njit
except ImportError:

    def njit(*args, **kwargs):
        def decorator(func):
            return func

        return decorator


from models import LayeredEarthModel
from forward_mfj import FullWavefieldMFJOperator


@dataclass
class BroadbandFusionConfig:
    """
    Configuration for active-passive f-c spectrum fusion.

    Attributes
    ----------
    f_min_global : float
        Minimum frequency of the fused broadband spectrum (Hz), default 1.0 Hz.
    f_max_global : float
        Maximum frequency of the fused broadband spectrum (Hz), default 55.0 Hz.
    f_overlap_low : float
        Lower boundary of the active-passive overlap transition band (Hz), default 8.0 Hz.
    f_overlap_high : float
        Upper boundary of the active-passive overlap transition band (Hz), default 12.0 Hz.
    norm_power : float
        Sharpness exponent applied to normalized frequency slices (default 1.5) to
        enhance modal ridges while preserving secondary higher and guided-P branches.
    """

    f_min_global: float = 1.0
    f_max_global: float = 55.0
    f_overlap_low: float = 8.0
    f_overlap_high: float = 12.0
    norm_power: float = 1.4


@njit(cache=True, fastmath=True)
def normalize_frequency_slices_numba(
    spec_fc: np.ndarray,
    power: float = 1.0,
) -> np.ndarray:
    """
    Normalize each frequency slice of an f-c matrix to [0, 1] with optional
    contrast exponent `power`.
    """
    n_f, n_c = spec_fc.shape
    out = np.zeros((n_f, n_c), dtype=np.float64)
    for n in range(n_f):
        s_min = spec_fc[n, 0]
        s_max = spec_fc[n, 0]
        for q in range(1, n_c):
            val = spec_fc[n, q]
            if val < s_min:
                s_min = val
            if val > s_max:
                s_max = val
        span = s_max - s_min
        if span > 1e-30:
            inv_span = 1.0 / span
            for q in range(n_c):
                norm_val = (spec_fc[n, q] - s_min) * inv_span
                if power != 1.0:
                    out[n, q] = norm_val**power
                else:
                    out[n, q] = norm_val
    return out


@njit(cache=True, fastmath=True)
def estimate_slice_variance_numba(spec_fc: np.ndarray) -> np.ndarray:
    """
    Estimate a spectral-shape variance proxy (not a calibrated noise variance)
    based on background floor median relative to peak modal energy.
    Lower background floor -> higher SNR -> smaller variance sigma^2(f_n).
    """
    n_f, n_c = spec_fc.shape
    var_est = np.zeros(n_f, dtype=np.float64)
    for n in range(n_f):
        row = spec_fc[n, :]
        mean_val = np.mean(row)
        max_val = np.max(row) + 1e-12
        # Ratio of mean energy to peak energy measures diffuse noise floor
        floor_ratio = mean_val / max_val
        var_est[n] = max(floor_ratio * floor_ratio, 1e-4)
    return var_est


class ActivePassiveSpectralFuser:
    """
    Fuses active MASW f-c spectra and passive ambient-noise cross-correlation f-c
    spectra across a shared broadband frequency grid with variance-weighted stacking
    in the 8-12 Hz transition zone.
    """

    def __init__(
        self,
        operator: FullWavefieldMFJOperator,
        config: Optional[BroadbandFusionConfig] = None,
    ) -> None:
        self.operator = operator
        self.config = config or BroadbandFusionConfig()

        # Evaluate only frequencies actually synthesized. Interpolation of sparse
        # spectra can invent or erase ridges and must not be used as new data.
        fa, fp = operator.freqs_active, operator.freqs_passive
        lo, hi = self.config.f_overlap_low, self.config.f_overlap_high
        # Inside the transition both channels must exist at the SAME native
        # frequency. A raw union would alternate pure-passive and joint slices
        # when grids differ, introducing artificial frequency striping.
        both = np.array(
            [
                f
                for f in fa
                if lo < f < hi and np.any(np.isclose(fp, f, rtol=0.0, atol=1e-9))
            ]
        )
        native = np.unique(np.concatenate((fp[fp <= lo], both, fa[fa >= hi])))
        self.freqs_broadband = np.ascontiguousarray(
            native[
                (native >= self.config.f_min_global)
                & (native <= self.config.f_max_global)
            ]
        )
        if len(self.freqs_broadband) < 1:
            raise ValueError("No native frequencies in the requested fusion band")
        self.c_grid = operator.c_grid

        # Precompute interpolation indices/weights from active and passive frequency grids
        self._setup_broadband_interpolation()

    def _setup_broadband_interpolation(self) -> None:
        """Precompute cosine-squared transition weights across the 8-12 Hz overlap band."""
        f_bb = self.freqs_broadband
        f_lo = self.config.f_overlap_low
        f_hi = self.config.f_overlap_high

        # Prior aperture confidence alpha(f) in [0, 1]:
        # alpha = 0 for f < 8 Hz (pure passive),
        # alpha = sin^2(pi/2 * (f - 8)/(12 - 8)) in [8, 12] Hz,
        # alpha = 1 for f > 12 Hz (pure active).
        alpha = np.zeros_like(f_bb)
        in_overlap = (f_bb >= f_lo) & (f_bb <= f_hi)
        above_overlap = f_bb > f_hi
        phase = 0.5 * np.pi * (f_bb[in_overlap] - f_lo) / max(f_hi - f_lo, 1e-6)
        alpha[in_overlap] = np.sin(phase) ** 2
        alpha[above_overlap] = 1.0
        self.alpha_aperture = alpha

    def _interp_spectrum_to_broadband(
        self, spec_src: np.ndarray, freqs_src: np.ndarray
    ) -> np.ndarray:
        """Copy matching native frequency slices onto the broadband union grid."""
        # Historical name retained for compatibility: this now copies exact
        # native slices only; missing frequencies remain unavailable.
        out = np.zeros((len(self.freqs_broadband), spec_src.shape[1]))
        for i, f in enumerate(self.freqs_broadband):
            matches = np.flatnonzero(np.isclose(freqs_src, f, rtol=0.0, atol=1e-9))
            if matches.size:
                out[i] = spec_src[matches[0]]
        return out

    def fuse_spectra(
        self,
        spec_active: np.ndarray,
        spec_passive: np.ndarray,
    ) -> Dict[str, np.ndarray]:
        """
        Normalize and fuse active and passive f-c matrices via inverse-variance
        spectral stacking across the overlapping 8-12 Hz transition band.

        Parameters
        ----------
        spec_active : np.ndarray
            Active MFJ spectrum on `operator.freqs_active`, shape (N_f_act, N_c).
        spec_passive : np.ndarray
            Passive MFJ spectrum on `operator.freqs_passive`, shape (N_f_pas, N_c).

        Returns
        -------
        result : dict
            Contains:
            - 'd_fused': Broadband fused f-c spectrum D(f, c), shape (N_f_bb, N_c)
            - 'd_active_bb': Normalized active spectrum on broadband grid
            - 'd_passive_bb': Normalized passive spectrum on broadband grid
            - 'w_active': Frequency-dependent heuristic active weights w_act(f)
            - 'w_passive': Frequency-dependent heuristic passive weights w_pas(f)
            - 'sigma2_active': Active spectral-shape proxy sigma^2_act(f)
            - 'sigma2_passive': Passive spectral-shape proxy sigma^2_pas(f)
        """
        # 1. Normalize individual active and passive spectra slice-by-slice
        norm_act = normalize_frequency_slices_numba(
            spec_active, power=self.config.norm_power
        )
        norm_pas = normalize_frequency_slices_numba(
            spec_passive, power=self.config.norm_power
        )

        # 2. Copy matching native slices to the broadband frequency union
        d_act_bb = self._interp_spectrum_to_broadband(
            norm_act, self.operator.freqs_active
        )
        d_pas_bb = self._interp_spectrum_to_broadband(
            norm_pas, self.operator.freqs_passive
        )

        # Re-normalize slices inside their native support bands
        d_act_bb = normalize_frequency_slices_numba(d_act_bb, power=1.0)
        d_pas_bb = normalize_frequency_slices_numba(d_pas_bb, power=1.0)

        # 3. Compute uncalibrated spectral-shape proxies
        var_act = estimate_slice_variance_numba(d_act_bb)
        var_pas = estimate_slice_variance_numba(d_pas_bb)

        # 4. Compute variance-weighted stacking coefficients in overlap band [8, 12] Hz
        alpha = self.alpha_aperture
        available_a = np.any(
            np.isclose(
                self.freqs_broadband[:, None],
                self.operator.freqs_active[None, :],
                rtol=0.0,
                atol=1e-9,
            ),
            axis=1,
        )
        available_p = np.any(
            np.isclose(
                self.freqs_broadband[:, None],
                self.operator.freqs_passive[None, :],
                rtol=0.0,
                atol=1e-9,
            ),
            axis=1,
        )
        alpha = np.where(
            available_a & ~available_p,
            1.0,
            np.where(available_p & ~available_a, 0.0, alpha),
        )
        prec_act = available_a * alpha / np.maximum(var_act, 1e-5)
        prec_pas = available_p * (1.0 - alpha) / np.maximum(var_pas, 1e-5)
        denom = prec_act + prec_pas + 1e-30

        w_act = prec_act / denom
        w_pas = prec_pas / denom

        # 5. Stack and re-normalize each broadband frequency slice
        stacked = w_act[:, None] * d_act_bb + w_pas[:, None] * d_pas_bb
        d_fused = normalize_frequency_slices_numba(stacked, power=1.0)

        return {
            "d_fused": d_fused,
            "d_active_bb": d_act_bb,
            "d_passive_bb": d_pas_bb,
            "w_active": w_act,
            "w_passive": w_pas,
            "sigma2_active": var_act,
            "sigma2_passive": var_pas,
        }

    def generate_fused_spectrum_for_model(
        self,
        model: LayeredEarthModel,
        add_noise: bool = False,
        rng: Optional[np.random.Generator] = None,
    ) -> Dict[str, np.ndarray]:
        """
        End-to-end convenience method: forward-model active & passive wavefields,
        apply the unified MFJ transform, and fuse into D(f, c).
        """
        spec_act, spec_pas, u_act, u_pas = self.operator.forward_active_passive_spectra(
            model, add_noise=add_noise, rng=rng
        )
        fused_dict = self.fuse_spectra(spec_act, spec_pas)
        fused_dict["spec_active_raw"] = spec_act
        fused_dict["spec_passive_raw"] = spec_pas
        fused_dict["u_active"] = u_act
        fused_dict["u_passive"] = u_pas
        return fused_dict
