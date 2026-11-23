"""P-SV surface compliance, numerical Sommerfeld synthesis and phase-only MFJ.

The legacy heuristic mode/residue approximation has been removed. Numerical
quadrature requires convergence checks; field noise is a phenomenological model.
External libraries are optional references and are not silently substituted.
"""

from __future__ import annotations
from dataclasses import dataclass
import shutil
import subprocess
import tempfile
from typing import Dict, List, Optional, Tuple
import numpy as np

try:
    from numba import njit, prange

    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False

    def njit(*args, **kwargs):
        def decorator(func):
            return func

        return decorator

    prange = range


@njit(cache=True, fastmath=True)
def _bessel_j0_y0_tensor_numba(kr: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Numba-accelerated evaluation of cylindrical Bessel functions J_0(x) and Y_0(x)
    over a 3D argument tensor `kr` of shape (N_f, N_c, N_r).
    Uses exact power series with harmonic numbers for x <= 12.0 and Hankel's
    high-order asymptotic expansion for x > 12.0.
    """
    n_f, n_c, n_r = kr.shape
    j0_out = np.zeros((n_f, n_c, n_r), dtype=np.float64)
    y0_out = np.zeros((n_f, n_c, n_r), dtype=np.float64)
    gamma_e = 0.5772156649015328606
    two_over_pi = 0.6366197723675813431
    pi_over_4 = 0.7853981633974483096

    for n in range(n_f):
        for q in range(n_c):
            for m in range(n_r):
                x = kr[n, q, m]
                if x < 1e-8:
                    x = 1e-8
                if x <= 12.0:
                    half_x2 = 0.25 * x * x
                    term = 1.0
                    j_sum = 1.0
                    y_sum = 0.0
                    harm = 0.0
                    for k_idx in range(1, 34):
                        term = -term * half_x2 / (k_idx * k_idx)
                        harm += 1.0 / k_idx
                        j_sum += term
                        y_sum -= term * harm
                        if abs(term) < 1e-16:
                            break
                    j0_val = j_sum
                    y0_val = two_over_pi * (
                        (np.log(0.5 * x) + gamma_e) * j0_val + y_sum
                    )
                else:
                    inv_x = 1.0 / x
                    inv_x2 = inv_x * inv_x
                    # Hankel asymptotic expansion coefficients P_0(x), Q_0(x)
                    p0 = 1.0 - inv_x2 * (
                        0.0703125
                        - inv_x2 * (0.112152099609375 - inv_x2 * 0.5725014209747314)
                    )
                    q0 = (
                        -0.125
                        * inv_x
                        * (1.0 - inv_x2 * (0.234375 - inv_x2 * 0.64599609375))
                    )
                    amp = np.sqrt(two_over_pi * inv_x)
                    chi = x - pi_over_4
                    c_chi = np.cos(chi)
                    s_chi = np.sin(chi)
                    j0_val = amp * (p0 * c_chi - q0 * s_chi)
                    y0_val = amp * (p0 * s_chi + q0 * c_chi)

                j0_out[n, q, m] = j0_val
                y0_out[n, q, m] = y0_val

    return j0_out, y0_out


try:
    import ccfj

    CCFJ_AVAILABLE = True
except ImportError:
    CCFJ_AVAILABLE = False

try:
    import pyprop8

    PYPROP8_AVAILABLE = True
except ImportError:
    PYPROP8_AVAILABLE = False

try:
    from disba import PhaseDispersion

    DISBA_AVAILABLE = True
except ImportError:
    DISBA_AVAILABLE = False

from models import LayeredEarthModel

# =============================================================================
# 1. Acquisition Geometry Specifications (Active MASW + Passive Ambient Noise)
# =============================================================================


@dataclass
class ActiveSurveyGeometry:
    """
    Active-source MASW acquisition geometry (Option 4 field geometry):
    - 48-channel linear vertical/3C geophone array (4.5 Hz)
    - 2.0 m receiver spacing (94.0 m spread: x in [0, 94] m)
    - 6 bidirectional multi-offset shots at r_0 in (-20, -10, -4, +4, +10, +20) m
      from the array ends, providing dense unique source-receiver radial offsets
      from 4.0 m to 114.0 m.
    - 25.0 Hz Ricker wavelet vertical point force covering 8 - 60 Hz.
    """

    n_channels: int = 48
    dx: float = 2.0
    near_offsets: Tuple[float, ...] = (-20.0, -10.0, -4.0, 4.0, 10.0, 20.0)
    f_center: float = 25.0
    f_min: float = 8.0
    f_max: float = 60.0
    dt: float = 0.0015
    nt: int = 1024

    @property
    def offsets_primary(self) -> np.ndarray:
        """Primary shot gather offsets (4 m near-offset, 4 m to 98 m)."""
        return 4.0 + np.arange(self.n_channels, dtype=np.float64) * self.dx

    @property
    def offsets_secondary(self) -> np.ndarray:
        """Secondary shot gather offsets (10 m near-offset, 10 m to 104 m)."""
        return 10.0 + np.arange(self.n_channels, dtype=np.float64) * self.dx

    @property
    def combined_unique_offsets(self) -> np.ndarray:
        """
        Combined sorted unique radial source-receiver offsets across all 6 bidirectional
        shots of Opsi 4 (±4 m, ±10 m, ±20 m -> r in [4.0, 114.0] m).
        """
        rec_x = np.arange(self.n_channels, dtype=np.float64) * self.dx
        spread_len = rec_x[-1]
        all_r_list = []
        for off in self.near_offsets:
            absoff = abs(float(off))
            if off < 0:
                shot_x = -absoff
            else:
                shot_x = spread_len + absoff
            all_r_list.append(np.abs(rec_x - shot_x))
        all_r = np.concatenate(all_r_list)
        return np.unique(np.round(all_r, 4))


@dataclass
class PassiveSurveyGeometry:
    """
    Passive-source ambient noise cross-correlation acquisition geometry.
    Directly coupled by default to Opsi 4 (Hybrid Nested 2D Golden-Angle Log-Spiral):
    - 25 nodal 3C stations: 9 dense inner core stations co-located along the 94 m active
      spine (dx ~ 10.5 m with 2D crossline stagger) + 16 outer 2D golden-angle logarithmic
      spiral stations extending out to 562.6 m maximum inter-station aperture.
    - Retains all actual inter-station distances by default (300 station pairs).
      Optional radial binning is explicit and may discard short-distance resolution.
    """

    n_stations: int = 25
    dx_station: float = 25.0
    f_min: float = 0.5
    f_max: float = 12.0
    dt: float = 0.015
    nt: int = 1024
    use_option4_2d_nested: bool = True
    n_radial_bins: int = 0
    coordinates: Optional[np.ndarray] = None

    @property
    def station_xy(self) -> np.ndarray:
        """Return 2D coordinates (N_stations, 2) for Opsi 4 Hybrid Nested 2D or 1D fallback."""
        if self.coordinates is not None:
            xy = np.asarray(self.coordinates, dtype=float)
            if xy.shape != (self.n_stations, 2):
                raise ValueError("coordinates must have shape (n_stations, 2)")
            return xy.copy()
        if self.use_option4_2d_nested and self.n_stations < 10:
            raise ValueError("nested geometry requires at least 10 stations")
        xy = np.zeros((self.n_stations, 2), dtype=np.float64)
        if not self.use_option4_2d_nested:
            xy[:, 0] = np.arange(self.n_stations, dtype=np.float64) * self.dx_station
            return xy

        # Option 4: 9 dense inner core stations + 16 golden-angle 2D logarithmic spiral stations
        core_x = np.linspace(4.0, 88.0, 9)
        core_y = np.array([0.0, 4.0, -4.0, 7.0, -7.0, 5.0, -5.0, 3.0, 0.0])
        xy[0:9, 0] = core_x
        xy[0:9, 1] = core_y

        golden_angle = np.pi * (3.0 - np.sqrt(5.0))
        radii_spiral = np.geomspace(52.0, 320.0, self.n_stations - 9)
        for k in range(self.n_stations - 9):
            ang = k * golden_angle + 0.35
            r_s = radii_spiral[k]
            xy[9 + k, 0] = 45.0 + r_s * np.cos(ang)
            xy[9 + k, 1] = r_s * np.sin(ang)
        return xy

    def _compute_binned_2d_pairs(
        self, dominant_noise_azimuth_deg: float = 25.0
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Compute all N*(N-1)/2 inter-station 2D distances r_ij and azimuths theta_ij,
        and stack them into continuous radial distance bins with their exact 2D
        azimuthal directional bias factor <cos(2*(theta_ij - theta_0))>.
        """
        xy = self.station_xy
        n = self.n_stations
        dists = []
        azims = []
        for i in range(n):
            for j in range(i + 1, n):
                dx = xy[j, 0] - xy[i, 0]
                dy = xy[j, 1] - xy[i, 1]
                r = np.hypot(dx, dy)
                if r > 1.0:
                    dists.append(r)
                    azims.append(np.arctan2(dy, dx))
        dists_arr = np.array(dists, dtype=np.float64)
        azims_arr = np.array(azims, dtype=np.float64)
        th0 = np.deg2rad(dominant_noise_azimuth_deg)

        if not self.use_option4_2d_nested and self.coordinates is None:
            unique_r = np.arange(1, n, dtype=np.float64) * self.dx_station
            weights = np.arange(n - 1, 0, -1, dtype=np.float64)
            az_factor = np.full_like(unique_r, np.cos(2.0 * th0))
            return unique_r, weights, az_factor

        if self.n_radial_bins <= 0:
            # Group coincident distances only; retain actual short lags and local azimuths.
            rounded = np.round(dists_arr, 8)
            unique, inverse, counts = np.unique(
                rounded, return_inverse=True, return_counts=True
            )
            factors = (
                np.bincount(inverse, weights=np.cos(2 * (azims_arr - th0))) / counts
            )
            return unique, counts.astype(float), factors

        # Bin the 300 2D station pairs into n_radial_bins continuous radial lags
        r_min, r_max = float(np.min(dists_arr)), float(np.max(dists_arr))
        edges = np.linspace(r_min * 0.99, r_max * 1.01, self.n_radial_bins + 1)
        r_centers = []
        r_counts = []
        r_az_bias = []

        for b in range(self.n_radial_bins):
            mask = (dists_arr >= edges[b]) & (dists_arr < edges[b + 1])
            cnt = int(np.sum(mask))
            if cnt > 0:
                r_centers.append(float(np.mean(dists_arr[mask])))
                r_counts.append(float(cnt))
                # Blended local + array-wide 2D azimuthal cancellation
                local_cos2 = float(np.mean(np.cos(2.0 * (azims_arr[mask] - th0))))
                r_az_bias.append(local_cos2)

        return (
            np.array(r_centers, dtype=np.float64),
            np.array(r_counts, dtype=np.float64),
            np.array(r_az_bias, dtype=np.float64),
        )

    @property
    def interstation_offsets(self) -> np.ndarray:
        """Radial distances r_m (10.9 m to 562.6 m for Opsi 4)."""
        r_c, _, _ = self._compute_binned_2d_pairs()
        return r_c

    @property
    def pair_redundancy_weights(self) -> np.ndarray:
        """Number of station pairs contributing to each radial bin."""
        _, cnt, _ = self._compute_binned_2d_pairs()
        return cnt

    @property
    def azimuthal_anisotropy_factors(self) -> np.ndarray:
        """Residual 2D azimuthal directional factor per radial bin."""
        _, _, az_f = self._compute_binned_2d_pairs()
        return az_f


# =============================================================================
# 2. Low-Level Numba Full-Wavefield P-SV Layered Green's Function Kernel
#    (Structured with flat typed memory layouts for Rust / C++ FFI porting)
# =============================================================================


@njit(cache=True, fastmath=True)
def _invert_2x2_numba(
    a00: complex, a01: complex, a10: complex, a11: complex
) -> Tuple[complex, complex, complex, complex]:
    """Closed-form 2x2 complex matrix inverse with regularization guard."""
    det = a00 * a11 - a01 * a10
    if abs(det) < 1e-30:
        det = det + (1e-30 + 1e-30j)
    inv_det = 1.0 / det
    return a11 * inv_det, -a01 * inv_det, -a10 * inv_det, a00 * inv_det


@njit(cache=True, fastmath=True)
def _psv_surface_compliance_single_wk(
    omega_c: complex,
    k_val: float,
    h_arr: np.ndarray,
    vs_arr: np.ndarray,
    vp_arr: np.ndarray,
    rho_arr: np.ndarray,
    qs_arr: np.ndarray,
    qp_arr: np.ndarray,
) -> complex:
    """
    Evaluate the exact 1D layered P-SV surface vertical dynamic compliance G_zz(omega, k)
    using the decaying-exponential
    Kennett Generalized Reflection / Stiffness matrix recursion (only decaying exponentials).

    All inputs use flat contiguous float64 arrays suitable for porting
    to a C++/Rust `extern "C"` function.
    """
    n_layers = len(vs_arr)
    k = k_val
    k2 = k * k
    w2 = omega_c * omega_c

    # --- Start at bottom half-space (layer index n_layers - 1) ---
    idx = n_layers - 1
    vs_c = vs_arr[idx] * (1.0 - 0.5j / qs_arr[idx])
    vp_c = vp_arr[idx] * (1.0 - 0.5j / qp_arr[idx])
    mu = rho_arr[idx] * vs_c * vs_c

    nu_p = np.sqrt(k2 - w2 / (vp_c * vp_c))
    nu_s = np.sqrt(k2 - w2 / (vs_c * vs_c))
    if nu_p.real < 0.0:
        nu_p = -nu_p
    if nu_s.real < 0.0:
        nu_s = -nu_s

    gamma = 2.0 * k2 - w2 / (vs_c * vs_c)

    # Down-going partition matrices in half-space: K_N = T_d * U_d^{-1}
    ud00, ud01 = k + 0.0j, nu_s
    ud10, ud11 = nu_p, k + 0.0j
    td00, td01 = -2.0 * mu * k * nu_p, -mu * gamma
    td10, td11 = -mu * gamma, -2.0 * mu * k * nu_s

    iud00, iud01, iud10, iud11 = _invert_2x2_numba(ud00, ud01, ud10, ud11)
    k00 = td00 * iud00 + td01 * iud10
    k01 = td00 * iud01 + td01 * iud11
    k10 = td10 * iud00 + td11 * iud10
    k11 = td10 * iud01 + td11 * iud11

    # --- Upward recursion through finite layers j = n_layers - 2 down to 0 ---
    for j in range(n_layers - 2, -1, -1):
        h = h_arr[j]
        vs_c = vs_arr[j] * (1.0 - 0.5j / qs_arr[j])
        vp_c = vp_arr[j] * (1.0 - 0.5j / qp_arr[j])
        mu = rho_arr[j] * vs_c * vs_c

        nu_p = np.sqrt(k2 - w2 / (vp_c * vp_c))
        nu_s = np.sqrt(k2 - w2 / (vs_c * vs_c))
        if nu_p.real < 0.0:
            nu_p = -nu_p
        if nu_s.real < 0.0:
            nu_s = -nu_s

        gamma = 2.0 * k2 - w2 / (vs_c * vs_c)

        ud00, ud01 = k + 0.0j, nu_s
        ud10, ud11 = nu_p, k + 0.0j
        uu00, uu01 = k + 0.0j, -nu_s
        uu10, uu11 = -nu_p, k + 0.0j

        td00, td01 = -2.0 * mu * k * nu_p, -mu * gamma
        td10, td11 = -mu * gamma, -2.0 * mu * k * nu_s
        tu00, tu01 = 2.0 * mu * k * nu_p, -mu * gamma
        tu10, tu11 = -mu * gamma, 2.0 * mu * k * nu_s

        # A = T_u - K_{j+1} * U_u
        a00 = tu00 - (k00 * uu00 + k01 * uu10)
        a01 = tu01 - (k00 * uu01 + k01 * uu11)
        a10 = tu10 - (k10 * uu00 + k11 * uu10)
        a11 = tu11 - (k10 * uu01 + k11 * uu11)

        # B = T_d - K_{j+1} * U_d
        b00 = td00 - (k00 * ud00 + k01 * ud10)
        b01 = td01 - (k00 * ud01 + k01 * ud11)
        b10 = td10 - (k10 * ud00 + k11 * ud10)
        b11 = td11 - (k10 * ud01 + k11 * ud11)

        ia00, ia01, ia10, ia11 = _invert_2x2_numba(a00, a01, a10, a11)

        # R_j = -A^{-1} * B
        r00 = -(ia00 * b00 + ia01 * b10)
        r01 = -(ia00 * b01 + ia01 * b11)
        r10 = -(ia10 * b00 + ia11 * b10)
        r11 = -(ia10 * b01 + ia11 * b11)

        # Decaying phase factors E_j = diag(ep, es)
        ep = np.exp(-nu_p * h)
        es = np.exp(-nu_s * h)

        # R_tilde = E_j * R_j * E_j
        rt00 = ep * r00 * ep
        rt01 = ep * r01 * es
        rt10 = es * r10 * ep
        rt11 = es * r11 * es

        # M_T = T_d + T_u * R_tilde
        mt00 = td00 + (tu00 * rt00 + tu01 * rt10)
        mt01 = td01 + (tu00 * rt01 + tu01 * rt11)
        mt10 = td10 + (tu10 * rt00 + tu11 * rt10)
        mt11 = td11 + (tu10 * rt01 + tu11 * rt11)

        # M_U = U_d + U_u * R_tilde
        mu00 = ud00 + (uu00 * rt00 + uu01 * rt10)
        mu01 = ud01 + (uu00 * rt01 + uu01 * rt11)
        mu10 = ud10 + (uu10 * rt00 + uu11 * rt10)
        mu11 = ud11 + (uu10 * rt01 + uu11 * rt11)

        imu00, imu01, imu10, imu11 = _invert_2x2_numba(mu00, mu01, mu10, mu11)

        # K_j = M_T * M_U^{-1}
        k00 = mt00 * imu00 + mt01 * imu10
        k01 = mt00 * imu01 + mt01 * imu11
        k10 = mt10 * imu00 + mt11 * imu10
        k11 = mt10 * imu01 + mt11 * imu11

    # Surface dynamic compliance G(0) = -K_1^{-1}
    ik00, ik01, ik10, ik11 = _invert_2x2_numba(k00, k01, k10, k11)
    g_zz = -ik11

    return g_zz


@njit(cache=True, fastmath=True)
def _psv_surface_compliance_grid_numba(
    freqs: np.ndarray,
    k_grid: np.ndarray,
    eta_damp: float,
    h_arr: np.ndarray,
    vs_arr: np.ndarray,
    vp_arr: np.ndarray,
    rho_arr: np.ndarray,
    qs_arr: np.ndarray,
    qp_arr: np.ndarray,
) -> np.ndarray:
    """
    Compute complex P-SV surface vertical dynamic compliance G_zz(omega_n, k_l)
    over a 2D frequency-wavenumber grid of shape (N_f, N_k).
    """
    n_f = len(freqs)
    n_k = len(k_grid)
    g_zz_mat = np.zeros((n_f, n_k), dtype=np.complex128)

    for n in range(n_f):
        omega_c = 2.0 * np.pi * freqs[n] + 1j * eta_damp
        for l in range(n_k):
            g_zz_mat[n, l] = _psv_surface_compliance_single_wk(
                omega_c,
                k_grid[l],
                h_arr,
                vs_arr,
                vp_arr,
                rho_arr,
                qs_arr,
                qp_arr,
            )
    return g_zz_mat


# =============================================================================
# 3. Modular Forward Adapters (PyProp8, CPS hspec96, and Built-In Numba Solver)
# =============================================================================


class PyProp8ForwardAdapter:
    """
    Modular adapter for `PyProp8` (O'Toole & Woodhouse reflectivity/propagator matrix
    library) when available in the Python environment.
    """

    @staticmethod
    def is_available() -> bool:
        return PYPROP8_AVAILABLE

    @staticmethod
    def compute_seismograms(
        model: LayeredEarthModel,
        offsets_m: np.ndarray,
        nt: int,
        dt: float,
        f_center: float = 25.0,
    ) -> Optional[np.ndarray]:
        if not PYPROP8_AVAILABLE:
            raise ImportError("Install pyprop8 for the independent elastic adapter")
        # SI units throughout; PyProp8 uses consistent caller-supplied units.
        # This adapter is elastic (Q is not included); do not compare its
        # amplitudes to attenuating native output without matching moduli.
        layers = [
            (
                model.thickness[i] if i < model.n_layers - 1 else np.inf,
                model.vp[i],
                model.vs[i],
                model.rho[i],
            )
            for i in range(model.n_layers)
        ]
        earth = pyprop8.LayeredStructureModel(layers)
        source = pyprop8.PointSource(
            0.0, 0.0, 1e-3, np.zeros((3, 3)), np.array([[0.0], [0.0], [1.0]]), 0.0
        )
        receivers = pyprop8.ListOfReceivers(
            np.asarray(offsets_m),
            np.zeros_like(offsets_m),
            depth=0.0,
            geometry="cartesian",
        )

        def ricker(omega):
            x = (omega / (2 * np.pi * f_center)) ** 2
            return x * np.exp(-x)

        _, seis = pyprop8.compute_seismograms(
            earth,
            source,
            receivers,
            nt,
            dt,
            source_time_function=ricker,
            show_progress=False,
            squeeze_outputs=False,
            stencil_kwargs={
                "kmin": 0.0,
                "kmax": max(8 / min(offsets_m), 6 * np.pi * f_center / min(model.vs)),
                "nk": 8000,
            },
        )
        return np.asarray(seis[0, :, 2, :].T, dtype=float)


class CPSHspec96Adapter:
    """
    Modular adapter for Herrmann's Computer Programs in Seismology (`hspec96` / `hpulse96`)
    wavenumber integration suite if installed on system PATH.
    """

    @staticmethod
    def is_available() -> bool:
        return (
            shutil.which("hspec96") is not None and shutil.which("hpulse96") is not None
        )

    @staticmethod
    def write_model96(model: LayeredEarthModel, filepath: str) -> None:
        """Write a CPS MODEL.01 format file."""
        lines = [
            "MODEL.01",
            model.name,
            "ISOTROPIC",
            "KGS",
            "FLAT EARTH",
            "1-D",
            "CONSTANT VELOCITY",
            "LINE08",
            "LINE09",
            "LINE10",
            "LINE11",
            "H(KM) VP(KM/S) VS(KM/S) RHO(GM/CC) QP QS ETAP ETAS FREFP FREFS",
        ]
        for i in range(model.n_layers):
            h_km = model.thickness[i] / 1000.0 if i < model.n_layers - 1 else 0.0
            vp_kms = model.vp[i] / 1000.0
            vs_kms = model.vs[i] / 1000.0
            rho_gcc = model.rho[i] / 1000.0
            lines.append(
                f"{h_km:.4f} {vp_kms:.4f} {vs_kms:.4f} {rho_gcc:.4f} "
                f"{model.qp[i]:.1f} {model.qs[i]:.1f} 0.0 0.0 1.0 1.0"
            )
        with open(filepath, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")


# =============================================================================
# 4. Phase-only one-way Hankel MFJ transform
# =============================================================================


@njit(cache=True, fastmath=True)
def _mfj_hankel_transform_numba(
    u_freq_r: np.ndarray,
    r_offsets: np.ndarray,
    dr_weights: np.ndarray,
    window_r: np.ndarray,
    h0_out_kernel: np.ndarray,
) -> np.ndarray:
    """
    Numba-accelerated discrete Modified Frequency-Bessel (MFJ) Transform using
    one-way outgoing Hankel kernel H_0^(2)(omega * r / c) = J_0(kr) - i * Y_0(kr).

    Parameters
    ----------
    u_freq_r : np.ndarray (complex128)
        Complex wavefield spectrum U(f_n, r_m), shape (N_f, N_r).
    r_offsets : np.ndarray (float64)
        Radial offset distances r_m in meters, shape (N_r,).
    dr_weights : np.ndarray (float64)
        Trapezoidal radial integration weights Delta_r_m, shape (N_r,).
    window_r : np.ndarray (float64)
        Spatial Hanning tapering window w(r_m), shape (N_r,).
    h0_out_kernel : np.ndarray (complex128)
        Precomputed one-way outgoing Hankel tensor H_0^(2)(2*pi*f_n * r_m / c_q),
        shape (N_f, N_c, N_r).

    Returns
    -------
    fc_spec : np.ndarray (float64)
        Normalized MFJ dispersion spectrogram I_MFJ(f_n, c_q), shape (N_f, N_c).
    """
    n_f, n_r = u_freq_r.shape
    n_c = h0_out_kernel.shape[1]
    fc_spec = np.zeros((n_f, n_c), dtype=np.float64)

    for n in range(n_f):
        # 1. Compute max amplitude across traces at frequency f_n for spectral whitening
        max_amp = 0.0
        for m in range(n_r):
            amp = abs(u_freq_r[n, m])
            if amp > max_amp:
                max_amp = amp
        reg_eps = 0.02 * max_amp + 1e-30

        # 2. Pre-weight trace spectrum with whitening, spatial Hanning taper, and cylindrical r * dr
        u_weighted = np.zeros(n_r, dtype=np.complex128)
        for m in range(n_r):
            amp = abs(u_freq_r[n, m])
            whitened = u_freq_r[n, m] / (amp + reg_eps)
            # Weight by sqrt(r) * dr so cylindrical spreading 1/sqrt(r) is balanced across array
            u_weighted[m] = (
                window_r[m] * whitened * np.sqrt(r_offsets[m]) * dr_weights[m]
            )

        # 3. Project onto one-way Hankel basis across phase velocities c_q
        slice_max = 0.0
        for q in range(n_c):
            acc = 0.0 + 0.0j
            for m in range(n_r):
                acc += u_weighted[m] * h0_out_kernel[n, q, m]
            val = abs(acc)
            fc_spec[n, q] = val
            if val > slice_max:
                slice_max = val

        # 4. Normalize frequency slice to [0, 1]
        if slice_max > 1e-30:
            inv_max = 1.0 / slice_max
            for q in range(n_c):
                fc_spec[n, q] *= inv_max

    return fc_spec


def compute_trapezoidal_weights(r_offsets: np.ndarray) -> np.ndarray:
    """Compute 1D trapezoidal quadrature weights Delta_r_m for arbitrary sorted offsets."""
    n = len(r_offsets)
    w = np.zeros(n, dtype=np.float64)
    if n == 1:
        w[0] = 1.0
        return w
    w[0] = 0.5 * (r_offsets[1] - r_offsets[0])
    w[-1] = 0.5 * (r_offsets[-1] - r_offsets[-2])
    for i in range(1, n - 1):
        w[i] = 0.5 * (r_offsets[i + 1] - r_offsets[i - 1])
    return w


def build_modified_hanning_taper(
    r_offsets: np.ndarray, edge_fraction: float = 0.25
) -> np.ndarray:
    """
    Construct a cosine-tapered (Tukey/Hanning) spatial window w(r) across the receiver
    spread to suppress finite-aperture truncation ringing while preserving array resolution.
    """
    r_min = r_offsets[0]
    r_max = r_offsets[-1]
    span = max(r_max - r_min, 1e-6)
    u = (r_offsets - r_min) / span
    w = np.ones_like(u)
    left = u < edge_fraction
    right = u > (1.0 - edge_fraction)
    w[left] = 0.5 * (1.0 - np.cos(np.pi * u[left] / edge_fraction))
    w[right] = 0.5 * (1.0 - np.cos(np.pi * (1.0 - u[right]) / edge_fraction))
    return w


def precompute_hankel_kernels(
    freqs: np.ndarray,
    c_grid: np.ndarray,
    r_offsets: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Precompute forward incoming/outgoing Hankel integration kernel and inverse
    Modified Frequency-Bessel (MFJ) one-way Hankel kernel on (f, c, r) grids.

    Returns
    -------
    h0_fwd : np.ndarray (complex128)
        Forward cylindrical wavefield synthesis kernel H_0^(1)(2*pi*f * r / c),
        shape (N_f, N_c, N_r).
    h0_mfj : np.ndarray (complex128)
        Inverse MFJ imaging kernel H_0^(2)(2*pi*f * r / c), shape (N_f, N_c, N_r).
    """
    # Compute kr argument tensor of shape (N_f, N_c, N_r)
    omega = 2.0 * np.pi * freqs[:, None, None]
    c = c_grid[None, :, None]
    r = r_offsets[None, None, :]
    kr = np.ascontiguousarray(np.maximum(omega * r / c, 1e-4), dtype=np.float64)

    j0_kr, y0_kr = _bessel_j0_y0_tensor_numba(kr)

    # Normalize the one-way cylindrical Hankel kernels H_0^(1)(kr) and H_0^(2)(kr)
    # by their cylindrical modulus sqrt(J_0^2(kr) + Y_0^2(kr)) so the exact cylindrical
    # near-field and far-field Hankel phase is preserved while completely eliminating
    # the Y_0(kr -> 0) logarithmic amplitude blow-up at high phase velocities c.
    modulus_kr = np.sqrt(j0_kr * j0_kr + y0_kr * y0_kr) + 1e-30
    j0_unit = j0_kr / modulus_kr
    y0_unit = y0_kr / modulus_kr

    h0_fwd = (j0_unit + 1j * y0_unit).astype(np.complex128)
    h0_mfj = (j0_unit - 1j * y0_unit).astype(np.complex128)
    return np.ascontiguousarray(h0_fwd), np.ascontiguousarray(h0_mfj)


# =============================================================================
# 5. Unified Full-Wavefield Forward + MFJ Operator Class
# =============================================================================


class FullWavefieldMFJOperator:
    """
    Fast, precomputed Full-Wavefield P-SV Green's Function + Modified Frequency-Bessel (MFJ)
    forward operator for both synthetic dataset generation (with phenomenological active/passive
    noise) and picking-free SMVI inversion.
    """

    def __init__(
        self,
        active_geom,
        passive_geom,
        freqs_active,
        freqs_passive,
        c_grid,
        quadrature_density=1.0,
    ):
        self.active_geom, self.passive_geom = active_geom, passive_geom
        self.freqs_active = np.ascontiguousarray(freqs_active, dtype=float)
        self.freqs_passive = np.ascontiguousarray(freqs_passive, dtype=float)
        self.c_grid = np.ascontiguousarray(c_grid, dtype=float)
        self.quadrature_density = float(quadrature_density)
        if self.quadrature_density <= 0:
            raise ValueError("quadrature_density must be positive")
        self.r_act = active_geom.combined_unique_offsets
        self.r_pas, self.pair_counts_pas, self.az_factors_pas = (
            passive_geom._compute_binned_2d_pairs()
        )
        self.dr_act = compute_trapezoidal_weights(self.r_act)
        self.dr_pas = compute_trapezoidal_weights(self.r_pas)
        self.win_act = build_modified_hanning_taper(self.r_act, 0.22)
        self.win_pas = build_modified_hanning_taper(self.r_pas, 0.18)
        _, self.h0_mfj_act = precompute_hankel_kernels(
            self.freqs_active, self.c_grid, self.r_act
        )
        _, self.h0_mfj_pas = precompute_hankel_kernels(
            self.freqs_passive, self.c_grid, self.r_pas
        )

    def quadrature_grid(self, frequency, offsets, cache=True):
        """Model-independent quadrature; shared by native and AD backends."""
        from scipy.special import j0

        if not hasattr(self, "_quadrature_cache"):
            self._quadrature_cache = {}
        key = (float(frequency), tuple(offsets), self.quadrature_density)
        if cache and key in self._quadrature_cache:
            return self._quadrature_cache[key]
        omega = 2 * np.pi * frequency
        ksplit = omega / 5000.0
        kmax = max(3 * omega / max(0.5 * self.c_grid[0], 40.0), 8 / np.min(offsets))
        k = np.concatenate(
            (
                np.linspace(1e-10, ksplit, 33)[:-1],
                np.geomspace(
                    ksplit, kmax, max(64, int(1024 * self.quadrature_density))
                ),
            )
        )
        weights = compute_trapezoidal_weights(k)
        basis = np.ascontiguousarray(
            j0(k[:, None] * offsets[None, :]) * weights[:, None]
        )
        if cache:
            self._quadrature_cache[key] = k, basis
        return k, basis

    def _green_gather(self, model, freqs, offsets, cache=True):
        """Sommerfeld J0 integral of the complete P-SV compliance.

        Integrate a static-subtracted kernel on a frequency-dependent composite linear/logarithmic k grid.
        Restore the analytic Boussinesq 1/r term. No imposed modal branches are used.
        quadrature_density must be increased to check numerical convergence.
        """
        h, vs, vp, rho, qs, qp = model.to_flat_arrays()
        u = np.empty((len(freqs), len(offsets)), complex)
        mu = rho[0] * (vs[0] * (1 - 0.5j / qs[0])) ** 2
        ratio2 = (vp[0] * (1 - 0.5j / qp[0]) / (vs[0] * (1 - 0.5j / qs[0]))) ** 2
        nu = (0.5 * ratio2 - 1) / (ratio2 - 1)
        static_coeff = (1 - nu) / mu
        for n, f in enumerate(freqs):
            eta = 0.005 * 2 * np.pi * f
            k, basis = self.quadrature_grid(f, offsets, cache)
            g = _psv_surface_compliance_grid_numba(
                np.array([f]), k, eta, h, vs, vp, rho, qs, qp
            )[0]
            residual = g * k - static_coeff
            u[n] = residual @ basis + static_coeff / offsets
        return u

    def synthesize_frequency_gather(
        self, model, mode="active", add_noise=False, rng=None
    ):
        if mode not in ("active", "passive"):
            raise ValueError("mode must be active or passive")
        rng = np.random.default_rng(42) if rng is None else rng
        active = mode == "active"
        freqs = self.freqs_active if active else self.freqs_passive
        offsets = self.r_act if active else self.r_pas
        if active:
            x = (freqs / self.active_geom.f_center) ** 2
            source = x * np.exp(-x)
        else:
            source = 1 / (1 + (freqs / 6) ** 1.5)
        u = self._green_gather(model, freqs, offsets) * source[:, None]
        if add_noise:
            if active:
                u *= rng.lognormal(0, 0.35, (1, len(offsets)))
                u *= np.exp(1j * rng.normal(0, 0.08, (1, len(offsets))))
                noise_scale = 0.08
            else:
                # Explicit phenomenological directional perturbation, not an ambient-source solver.
                u *= np.exp(
                    1j
                    * 0.045
                    * self.az_factors_pas[None, :]
                    * 2
                    * np.pi
                    * freqs[:, None]
                    * offsets[None, :]
                    / 260
                )
                noise_scale = (
                    0.045 / np.sqrt(np.maximum(self.pair_counts_pas, 1))[None, :]
                )
            rms = np.sqrt(np.mean(np.abs(u) ** 2, axis=1, keepdims=True))
            u += (
                noise_scale
                * rms
                * (rng.normal(size=u.shape) + 1j * rng.normal(size=u.shape))
                / np.sqrt(2)
            )
        return np.ascontiguousarray(u), None

    def synthesize_time_gather(self, model, mode="active", add_noise=True, rng=None):
        """Evaluate the Green function directly on the FFT grid; conjugate for irfft."""
        rng = np.random.default_rng(42) if rng is None else rng
        active = mode == "active"
        geom = self.active_geom if active else self.passive_geom
        offsets = self.r_act if active else self.r_pas
        nt = max(2048, geom.nt)
        freqs = np.fft.rfftfreq(nt, geom.dt)
        mask = (freqs >= geom.f_min) & (freqs <= geom.f_max)
        fb = freqs[mask]
        if active:
            x = (fb / geom.f_center) ** 2
            source = x * np.exp(-x)
        else:
            source = 1 / (1 + (fb / 6) ** 1.5)
        ub = self._green_gather(model, fb, offsets, cache=False) * source[:, None]
        if add_noise:
            rms = np.sqrt(np.mean(np.abs(ub) ** 2))
            ub += (
                0.04
                * rms
                * (rng.normal(size=ub.shape) + 1j * rng.normal(size=ub.shape))
                / np.sqrt(2)
            )
        spec = np.zeros((len(freqs), len(offsets)), complex)
        t0 = 0.045 if active else 0.12
        spec[mask] = np.conj(ub) * np.exp(-2j * np.pi * fb[:, None] * t0)
        return np.arange(nt) * geom.dt, offsets, np.fft.irfft(spec, n=nt, axis=0)

    def compute_mfj_spectrum(
        self,
        u_fr: np.ndarray,
        mode: str = "active",
        use_ccfj_if_available: bool = True,
    ) -> np.ndarray:
        """
        Transform a frequency-offset wavefield gather U(f_n, r_m) into the
        frequency-phase velocity (f-c) domain using the native phase-only one-way Hankel MFJ transform.
        The legacy use_ccfj_if_available argument is ignored; no backend switching occurs.
        """
        if mode == "active":
            freqs = self.freqs_active
            r_off = self.r_act
            dr_w = self.dr_act
            win_r = self.win_act
            h0_mfj = self.h0_mfj_act
        else:
            freqs = self.freqs_passive
            r_off = self.r_pas
            dr_w = self.dr_pas
            win_r = self.win_pas
            h0_mfj = self.h0_mfj_pas

        return _mfj_hankel_transform_numba(u_fr, r_off, dr_w, win_r, h0_mfj)

    def forward_active_passive_spectra(
        self,
        model: LayeredEarthModel,
        add_noise: bool = False,
        rng: Optional[np.random.Generator] = None,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Run the complete full-wavefield forward modeling + MFJ transform for both
        active and passive surveys on a candidate Earth model.

        Returns
        -------
        spec_act : np.ndarray (float64)
            Normalized active MFJ f-c spectrogram, shape (N_f_act, N_c).
        spec_pas : np.ndarray (float64)
            Normalized passive MFJ f-c spectrogram, shape (N_f_pas, N_c).
        u_act : np.ndarray (complex128)
            Active frequency-offset wavefield gather, shape (N_f_act, N_r_act).
        u_pas : np.ndarray (complex128)
            Passive frequency-offset EGF gather, shape (N_f_pas, N_r_pas).
        """
        u_act, _ = self.synthesize_frequency_gather(
            model, mode="active", add_noise=add_noise, rng=rng
        )
        u_pas, _ = self.synthesize_frequency_gather(
            model, mode="passive", add_noise=add_noise, rng=rng
        )
        spec_act = self.compute_mfj_spectrum(u_act, mode="active")
        spec_pas = self.compute_mfj_spectrum(u_pas, mode="passive")
        return spec_act, spec_pas, u_act, u_pas


# =============================================================================
# 6. Independent theoretical Rayleigh dispersion (disba)
# =============================================================================


def compute_theoretical_dispersion_curves(
    model: LayeredEarthModel,
    freqs: np.ndarray,
    c_grid: np.ndarray,
    n_rayleigh_modes: int = 3,
) -> Dict[str, List[Tuple[np.ndarray, np.ndarray]]]:
    """
    Obtain independent Rayleigh-mode curves from disba if available.
    No heuristic fallback or unvalidated Guided-P reference is returned.

    Returns
    -------
    curves : dict
        Dictionary with keys:
        - 'rayleigh_modes': list of (f_branch, c_branch) tuples for M0, M1, M2
        - 'guided_p_modes': list of (f_branch, c_branch) tuples for Guided-P branches
    """
    rayleigh_modes: List[Tuple[np.ndarray, np.ndarray]] = []
    guided_p_modes: List[Tuple[np.ndarray, np.ndarray]] = []

    # 1. Try disba for exact normal-mode Rayleigh curves
    if DISBA_AVAILABLE:
        try:
            # disba expects thickness, vp, vs, rho in km, km/s, g/cm^3
            h_km = model.thickness.copy() / 1000.0
            h_km[-1] = 0.0
            vp_kms = model.vp / 1000.0
            vs_kms = model.vs / 1000.0
            rho_gcc = model.rho / 1000.0
            periods = 1.0 / np.sort(freqs)[::-1]
            pd = PhaseDispersion(h_km, vp_kms, vs_kms, rho_gcc)
            for m_idx in range(n_rayleigh_modes):
                try:
                    c_res = pd(periods, mode=m_idx, wave="rayleigh")
                    if len(c_res.period) > 2:
                        f_br = (1.0 / c_res.period)[::-1]
                        c_br = (c_res.velocity * 1000.0)[::-1]
                        valid = (c_br >= c_grid[0]) & (c_br <= c_grid[-1])
                        if np.sum(valid) > 2:
                            rayleigh_modes.append((f_br[valid], c_br[valid]))
                except Exception:
                    continue
        except Exception:
            pass

    # No heuristic branches are advertised as theoretical eigenmodes. Guided/leaky
    # pole extraction remains unvalidated and therefore has no reference overlay.
    return {"rayleigh_modes": rayleigh_modes, "guided_p_modes": guided_p_modes}
