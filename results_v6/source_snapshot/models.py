"""
models.py
=========
Synthetic 1D Layered Earth Models and Parameterization Utilities for
Joint Active-Passive Picking-Free Full f-c Spectrum Inversion.

Defines a low-Vs-layer model and a shallow Vp-contrast model for controlled
synthetic tests. Mode jumping, guided-P identification and simultaneous Vp/Vs
resolution are research questions, not properties established by these fixtures.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import numpy as np


@dataclass
class LayeredEarthModel:
    """
    1D Horizontally Layered Isotropic Viscoelastic Earth Model.

    Attributes
    ----------
    name : str
        Descriptive name of the benchmark scenario.
    code : str
        Short identifier ('benchmark_a' or 'benchmark_b').
    thickness : np.ndarray
        Layer thicknesses in meters, shape (N,). The last entry (half-space) is 0.0 or np.inf.
    vs : np.ndarray
        Shear-wave velocity (Vs) in m/s, shape (N,).
    vp : np.ndarray
        Compressional-wave velocity (Vp) in m/s, shape (N,).
    rho : np.ndarray
        Bulk density in kg/m^3, shape (N,).
    qs : np.ndarray
        Shear quality factor Qs (dimensionless), shape (N,).
    qp : np.ndarray
        Compressional quality factor Qp (dimensionless), shape (N,).
    description : str
        Physical and geological description of the model.
    """

    name: str
    code: str
    thickness: np.ndarray
    vs: np.ndarray
    vp: np.ndarray
    rho: np.ndarray
    qs: np.ndarray
    qp: np.ndarray
    description: str = ""

    @property
    def n_layers(self) -> int:
        return len(self.vs)

    @property
    def poisson_ratio(self) -> np.ndarray:
        """Compute dynamic Poisson's ratio nu from Vp/Vs ratio."""
        r2 = (self.vp / self.vs) ** 2
        return (0.5 * r2 - 1.0) / (r2 - 1.0)

    @property
    def vp_vs_ratio(self) -> np.ndarray:
        return self.vp / self.vs

    @property
    def interface_depths(self) -> np.ndarray:
        """Cumulative interface depths z_1, z_2, ..., z_{N-1} in meters."""
        return np.cumsum(self.thickness[:-1])

    def to_flat_arrays(
        self,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Return contiguous float64 flat arrays suitable for Numba / Rust / C++ FFI kernels.
        """
        h = np.ascontiguousarray(self.thickness, dtype=np.float64)
        vs = np.ascontiguousarray(self.vs, dtype=np.float64)
        vp = np.ascontiguousarray(self.vp, dtype=np.float64)
        rho = np.ascontiguousarray(self.rho, dtype=np.float64)
        qs = np.ascontiguousarray(self.qs, dtype=np.float64)
        qp = np.ascontiguousarray(self.qp, dtype=np.float64)
        return h, vs, vp, rho, qs, qp

    def summary_table(self) -> str:
        """Return a formatted ASCII summary table of the layered model."""
        lines = [
            f"=== {self.name} ({self.code}) ===",
            f"{'Layer':<6} {'h (m)':<10} {'Depth (m)':<12} {'Vs (m/s)':<10} "
            f"{'Vp (m/s)':<10} {'Vp/Vs':<8} {'Poisson':<9} {'Rho (kg/m3)':<12}",
            "-" * 82,
        ]
        z_top = 0.0
        nu = self.poisson_ratio
        for i in range(self.n_layers):
            h_val = self.thickness[i]
            if i < self.n_layers - 1 and np.isfinite(h_val) and h_val > 0:
                h_str = f"{h_val:.2f}"
                z_bot = z_top + h_val
                z_str = f"{z_top:.1f}-{z_bot:.1f}"
                z_top = z_bot
            else:
                h_str = "inf"
                z_str = f">{z_top:.1f}"
            lines.append(
                f"{i+1:<6} {h_str:<10} {z_str:<12} {self.vs[i]:<10.1f} "
                f"{self.vp[i]:<10.1f} {self.vp[i]/self.vs[i]:<8.2f} "
                f"{nu[i]:<9.4f} {self.rho[i]:<12.1f}"
            )
        return "\n".join(lines)


@dataclass
class InversionParameterSpace:
    """
    Defines the prior box bounds and bijective logit/sigmoid transforms for the
    11-dimensional inversion parameter vector:
        m = [Vs_1, Vs_2, Vs_3, Vs_4, Vp_1, Vp_2, Vp_3, Vp_4, h_1, h_2, h_3]
    where density rho and Qs/Qp are fixed nuisance profiles specified by the
    synthetic experiment. Vp uses a conditional lower bound of 1.45 Vs.
    """

    reference_model: LayeredEarthModel
    vs_min: np.ndarray
    vs_max: np.ndarray
    vp_min: np.ndarray
    vp_max: np.ndarray
    h_min: np.ndarray
    h_max: np.ndarray
    param_names: List[str] = field(init=False)

    def __post_init__(self) -> None:
        n = self.reference_model.n_layers
        self.param_names = (
            [f"Vs_{i+1}" for i in range(n)]
            + [f"Vp_{i+1}" for i in range(n)]
            + [f"h_{i+1}" for i in range(n - 1)]
        )

    @property
    def dim(self) -> int:
        return len(self.param_names)

    @property
    def lower_bounds(self) -> np.ndarray:
        return np.concatenate([self.vs_min, self.vp_min, self.h_min]).astype(np.float64)

    @property
    def upper_bounds(self) -> np.ndarray:
        return np.concatenate([self.vs_max, self.vp_max, self.h_max]).astype(np.float64)

    @property
    def true_vector(self) -> np.ndarray:
        return self.model_to_vector(self.reference_model)

    def model_to_vector(self, model: LayeredEarthModel) -> np.ndarray:
        """Pack a LayeredEarthModel into a 1D parameter vector m in R^d."""
        return np.concatenate(
            [
                model.vs,
                model.vp,
                model.thickness[:-1],
            ]
        ).astype(np.float64)

    def vector_to_model(
        self, m_vec: np.ndarray, name_suffix: str = "sample"
    ) -> LayeredEarthModel:
        """Unpack a 1D parameter vector m in R^d into a valid LayeredEarthModel."""
        n = self.reference_model.n_layers
        m_vec = np.asarray(m_vec, dtype=float)
        if m_vec.shape != (self.dim,) or not np.all(np.isfinite(m_vec)):
            raise ValueError("invalid model vector")
        if np.any(m_vec < self.lower_bounds) or np.any(m_vec > self.upper_bounds):
            raise ValueError("model outside prior bounds")
        vs, vp = m_vec[:n].copy(), m_vec[n : 2 * n].copy()
        if np.any(vp < 1.45 * vs):
            raise ValueError("Vp must be at least 1.45 Vs")
        thickness = np.append(m_vec[2 * n :], 0.0)
        # Fixed nuisance density and Q are explicit assumptions of the synthetic test.
        rho = self.reference_model.rho.copy()

        return LayeredEarthModel(
            name=f"{self.reference_model.name} ({name_suffix})",
            code=self.reference_model.code,
            thickness=thickness,
            vs=vs,
            vp=vp,
            rho=rho,
            qs=self.reference_model.qs.copy(),
            qp=self.reference_model.qp.copy(),
        )

    def _conditional_lower(self, physical):
        low = np.broadcast_to(self.lower_bounds, np.shape(physical)).copy()
        n = self.reference_model.n_layers
        low[..., n : 2 * n] = np.maximum(low[..., n : 2 * n], 1.45 * physical[..., :n])
        return low

    def physical_to_latent(self, m_phys):
        m = np.asarray(m_phys, dtype=float)
        low = self._conditional_lower(m)
        u = (m - low) / (self.upper_bounds - low)
        if np.any((u <= 0) | (u >= 1)):
            raise ValueError("logit requires interior physical parameters")
        return np.log(u) - np.log1p(-u)

    def latent_to_physical(self, z_latent):
        from scipy.special import expit

        z = np.asarray(z_latent, dtype=float)
        sig = expit(z)
        m = self.lower_bounds + (self.upper_bounds - self.lower_bounds) * sig
        low = self._conditional_lower(m)
        n = self.reference_model.n_layers
        m[..., n : 2 * n] = (
            low[..., n : 2 * n]
            + (self.upper_bounds[n : 2 * n] - low[..., n : 2 * n]) * sig[..., n : 2 * n]
        )
        return m

    def latent_jacobian_diag(self, z_latent):
        from scipy.special import expit

        z = np.asarray(z_latent, dtype=float)
        m = self.latent_to_physical(z)
        return (self.upper_bounds - self._conditional_lower(m)) * expit(z) * expit(-z)

    def latent_log_abs_det_jacobian(self, z_latent):
        z = np.asarray(z_latent, dtype=float)
        m = self.latent_to_physical(z)
        widths = self.upper_bounds - self._conditional_lower(m)
        return np.sum(
            np.log(widths) - np.logaddexp(0, -z) - np.logaddexp(0, z), axis=-1
        )


def get_benchmark_a_lvl() -> Tuple[LayeredEarthModel, InversionParameterSpace]:
    """
    Construct Benchmark A: 4-Layer Onshore Sedimentary/Pyroclastic Profile
    with a Trapped Soft Low-Velocity Layer (LVL) in Layer 2.

    Physics Goal:
    - Vs drops sharply from 260 m/s in Layer 1 to 150 m/s in Layer 2 (4-10 m depth)
      before rising to 390 m/s in Layer 3 and 680 m/s in the bedrock half-space.
    - Modal energy exchange is to be measured, not imposed in the model definition.
    """
    model = LayeredEarthModel(
        name="Benchmark A: Low-Velocity Layer (Mode-Jumping & Osculation Test)",
        code="benchmark_a",
        thickness=np.array([4.0, 6.0, 14.0, 0.0], dtype=np.float64),
        vs=np.array([260.0, 150.0, 390.0, 680.0], dtype=np.float64),
        vp=np.array([540.0, 480.0, 880.0, 1480.0], dtype=np.float64),
        rho=np.array([1780.0, 1620.0, 1960.0, 2180.0], dtype=np.float64),
        qs=np.array([25.0, 18.0, 35.0, 55.0], dtype=np.float64),
        qp=np.array([50.0, 40.0, 75.0, 110.0], dtype=np.float64),
        description=(
            "4-layer onshore sedimentary/pyroclastic profile containing a trapped soft "
            "Low-Velocity Layer (Vs=150 m/s at 4-10 m depth) beneath a stiffer surface crust "
            "(Vs=260 m/s); used to investigate possible modal complexity."
        ),
    )

    param_space = InversionParameterSpace(
        reference_model=model,
        vs_min=np.array([170.0, 95.0, 270.0, 500.0], dtype=np.float64),
        vs_max=np.array([350.0, 240.0, 510.0, 860.0], dtype=np.float64),
        vp_min=np.array([380.0, 320.0, 620.0, 1100.0], dtype=np.float64),
        vp_max=np.array([760.0, 720.0, 1200.0, 1900.0], dtype=np.float64),
        h_min=np.array([2.2, 3.5, 9.0], dtype=np.float64),
        h_max=np.array([6.5, 9.5, 20.0], dtype=np.float64),
    )
    return model, param_space


def get_benchmark_b_water_table() -> Tuple[LayeredEarthModel, InversionParameterSpace]:
    """
    Construct Benchmark B: 4-Layer Saturated Coastal Alluvial Profile
    with a Shallow Groundwater Table (Guided-P Sensitivity Test).

    Physics Goal:
    - Vs increases monotonically and smoothly with depth (175 -> 235 -> 350 -> 590 m/s).
    - Vp jumps sharply across the groundwater table at z = 3.5 m from 550 m/s
      (unsaturated vadose zone, nu = 0.444) to 1650 m/s (saturated aquifer,
      nu = 0.490 > 0.45). Identification of guided/leaky modes requires further validation.
    """
    model = LayeredEarthModel(
        name="Benchmark B: Shallow Groundwater Table (Guided-P Sensitivity Test)",
        code="benchmark_b",
        thickness=np.array([3.5, 6.5, 14.0, 0.0], dtype=np.float64),
        vs=np.array([175.0, 235.0, 350.0, 590.0], dtype=np.float64),
        vp=np.array([550.0, 1650.0, 1780.0, 1980.0], dtype=np.float64),
        rho=np.array([1660.0, 1980.0, 2060.0, 2200.0], dtype=np.float64),
        qs=np.array([22.0, 28.0, 38.0, 55.0], dtype=np.float64),
        qp=np.array([45.0, 85.0, 95.0, 120.0], dtype=np.float64),
        description=(
            "4-layer saturated coastal alluvial profile where Vp jumps sharply at the "
            "water table (z=3.5 m, Vp: 550 -> 1650 m/s, Poisson's ratio = 0.490) while "
            "Vs increases smoothly, exciting prominent leaky/guided-P modes."
        ),
    )

    param_space = InversionParameterSpace(
        reference_model=model,
        vs_min=np.array([115.0, 160.0, 250.0, 440.0], dtype=np.float64),
        vs_max=np.array([240.0, 320.0, 460.0, 750.0], dtype=np.float64),
        vp_min=np.array([360.0, 1200.0, 1350.0, 1550.0], dtype=np.float64),
        vp_max=np.array([820.0, 2050.0, 2200.0, 2450.0], dtype=np.float64),
        h_min=np.array([2.0, 3.8, 9.0], dtype=np.float64),
        h_max=np.array([5.8, 9.8, 20.0], dtype=np.float64),
    )
    return model, param_space


def evaluate_depth_profiles(
    models_vec: np.ndarray,
    param_space: InversionParameterSpace,
    z_grid: np.ndarray,
) -> Dict[str, np.ndarray]:
    """
    Evaluate an ensemble of 1D layered models (N_samples, dim) onto a regular
    depth grid z_grid (N_z,) to compute posterior mean, median, and 95% credible
    intervals for Vs(z), Vp(z), and Poisson's ratio nu(z).

    Parameters
    ----------
    models_vec : np.ndarray
        Array of physical model vectors, shape (N_samples, dim).
    param_space : InversionParameterSpace
        Inversion parameter space object.
    z_grid : np.ndarray
        1D array of depth points in meters, shape (N_z,).

    Returns
    -------
    profiles : dict
        Dictionary containing arrays of shape (N_samples, N_z) for 'vs', 'vp', 'nu',
        plus summary statistics ('vs_mean', 'vs_p025', 'vs_p975', etc.).
    """
    n_samples = models_vec.shape[0]
    n_z = len(z_grid)
    n_layers = param_space.reference_model.n_layers

    vs_prof = np.zeros((n_samples, n_z), dtype=np.float64)
    vp_prof = np.zeros((n_samples, n_z), dtype=np.float64)
    nu_prof = np.zeros((n_samples, n_z), dtype=np.float64)

    for s in range(n_samples):
        m = param_space.vector_to_model(models_vec[s])
        z_int = m.interface_depths
        nu = m.poisson_ratio
        # Determine layer index for each depth point in z_grid
        layer_idx = np.searchsorted(z_int, z_grid, side="right")
        layer_idx = np.clip(layer_idx, 0, n_layers - 1)
        vs_prof[s, :] = m.vs[layer_idx]
        vp_prof[s, :] = m.vp[layer_idx]
        nu_prof[s, :] = nu[layer_idx]

    out: Dict[str, np.ndarray] = {
        "z_grid": z_grid,
        "vs_samples": vs_prof,
        "vp_samples": vp_prof,
        "nu_samples": nu_prof,
    }
    for key, arr in [("vs", vs_prof), ("vp", vp_prof), ("nu", nu_prof)]:
        out[f"{key}_mean"] = np.mean(arr, axis=0)
        out[f"{key}_median"] = np.median(arr, axis=0)
        out[f"{key}_p025"] = np.percentile(arr, 2.5, axis=0)
        out[f"{key}_p975"] = np.percentile(arr, 97.5, axis=0)
        out[f"{key}_p16"] = np.percentile(arr, 16.0, axis=0)
        out[f"{key}_p84"] = np.percentile(arr, 84.0, axis=0)

    return out


if __name__ == "__main__":
    mod_a, space_a = get_benchmark_a_lvl()
    mod_b, space_b = get_benchmark_b_water_table()
    print(mod_a.summary_table())
    print()
    print(mod_b.summary_table())
