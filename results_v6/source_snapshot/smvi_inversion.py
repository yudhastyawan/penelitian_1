"""Picking-free inference with exact discrete W2 and consistent full covariance.

SMVI uses equal-weight components and Stein transport on their complete free
parameter vectors. Seismic posterior calibration remains unvalidated. K=1 uses a Gaussian reparameterized ELBO baseline.
Gradients of the target use central finite differences; no autodiff is claimed.
"""

from dataclasses import dataclass
import time
import numpy as np
from scipy.special import logsumexp
from scipy.linalg import solve_triangular
from models import LayeredEarthModel, InversionParameterSpace
from fusion import ActivePassiveSpectralFuser

try:
    from numba import njit
except ImportError:

    def njit(*args, **kwargs):
        return lambda fn: fn


from importlib.util import find_spec

TORCH_AVAILABLE = find_spec("torch") is not None


@njit(cache=True)
def discrete_wasserstein2(p, q, positions):
    """Exact monotone transport cost between two discrete probability masses."""
    if len(p) != len(q) or len(p) != len(positions) or len(p) < 2:
        raise ValueError("transport arrays must have matching lengths >= 2")
    for idx in range(len(p)):
        if (
            not np.isfinite(p[idx])
            or not np.isfinite(q[idx])
            or p[idx] < 0
            or q[idx] < 0
        ):
            raise ValueError("transport masses must be finite and nonnegative")
        if not np.isfinite(positions[idx]) or (
            idx > 0 and positions[idx] <= positions[idx - 1]
        ):
            raise ValueError(
                "transport positions must be finite and strictly increasing"
            )
    total_p, total_q = np.sum(p), np.sum(q)
    if (
        not np.isfinite(total_p)
        or not np.isfinite(total_q)
        or total_p <= 0
        or total_q <= 0
    ):
        raise ValueError("transport masses must have finite positive totals")
    i = j = 0
    rp, rq = p[0] / total_p, q[0] / total_q
    cost = 0.0
    while i < len(p) and j < len(q):
        mass = min(rp, rq)
        cost += mass * (positions[i] - positions[j]) ** 2
        rp -= mass
        rq -= mass
        if rp <= 1e-15:
            i += 1
            if i < len(p):
                rp = p[i] / total_p
        if rq <= 1e-15:
            j += 1
            if j < len(q):
                rq = q[j] / total_q
    return cost


@njit(cache=True)
def compute_spectral_misfit_numba(
    d_obs,
    d_syn,
    c_grid,
    freq_weights,
    wasserstein_weight=0.35,
    sharpen_beta=2.0,
    freqs=None,
    passive_cutoff_hz=12.0,
    passive_velocity_max=680.0,
    wasserstein_scale_mps=200.0,
):
    """Slice NCC plus exact discrete W2, with an explicit frequency/velocity mask."""
    n_f, n_c = d_obs.shape
    if d_syn.shape != d_obs.shape or len(c_grid) != n_c or len(freq_weights) != n_f:
        raise ValueError("incompatible spectrum, velocity or frequency-weight shape")
    if not np.isfinite(wasserstein_scale_mps) or wasserstein_scale_mps <= 0:
        raise ValueError("positive finite Wasserstein scale required")
    if not np.all(np.isfinite(d_obs)) or not np.all(np.isfinite(d_syn)):
        raise ValueError("spectra must be finite")
    if not np.all(np.isfinite(freq_weights)) or np.any(freq_weights < 0):
        raise ValueError("frequency weights must be finite and nonnegative")
    if n_c < 2 or not np.all(np.isfinite(c_grid)) or np.any(np.diff(c_grid) <= 0):
        raise ValueError("velocity grid must be finite and strictly increasing")
    if freqs is not None and (len(freqs) != n_f or not np.all(np.isfinite(freqs))):
        raise ValueError("incompatible or nonfinite frequency grid")
    ncc_acc = w2_acc = wsum = 0.0
    for n in range(n_f):
        if freq_weights[n] <= 0:
            continue
        limit = n_c
        if freqs is not None and freqs[n] <= passive_cutoff_hz:
            limit = max(2, np.searchsorted(c_grid, passive_velocity_max, side="right"))
            limit = min(limit, n_c)
        obs, syn = d_obs[n, :limit], d_syn[n, :limit]
        no, ns = np.sqrt(np.sum(obs**2)), np.sqrt(np.sum(syn**2))
        if no <= 1e-30 and ns <= 1e-30:
            continue
        ncc = np.sum(obs * syn) / max(no * ns, 1e-30)
        # Explicit floor regularization, with no threshold-dependent hidden ridge removal.
        po = np.maximum(obs, 0.0) ** sharpen_beta + 1e-12
        pq = np.maximum(syn, 0.0) ** sharpen_beta + 1e-12
        w2 = discrete_wasserstein2(po, pq, c_grid[:limit])
        wf = freq_weights[n]
        ncc_acc += wf * (1 - ncc)
        w2_acc += wf * w2
        wsum += wf
    if wsum == 0:
        return 0.0, 0.0, 0.0
    jn, jw = ncc_acc / wsum, w2_acc / wsum
    return (
        (1 - wasserstein_weight) * jn
        + wasserstein_weight * jw / wasserstein_scale_mps**2,
        jn,
        np.sqrt(jw),
    )


@dataclass
class SMVIConfig:
    n_components: int = 6
    n_iterations: int = 55
    n_mc_samples_per_comp: int = 2
    lr_mean: float = 0.04
    lr_chol: float = 0.015
    sigma_likelihood: float = 0.08
    wasserstein_weight: float = 0.35
    repulsion_strength: float = 0.2
    entropy_weight: float = 1.0
    n_posterior_draws: int = 600
    random_seed: int = 2026
    passive_cutoff_hz: float = 12.0
    passive_velocity_max: float = 680.0
    gradient_backend: str = "finite_difference"
    initialization_refinement_steps: int = 0
    # None resolves once from the FULL parent grid, before any velocity mask.
    wasserstein_scale_mps: float | None = None
    predictive_draws: int = 128
    predictive_max_draws: int = 512
    predictive_bootstrap_repeats: int = 32
    predictive_mc_tolerance: float = 0.01


class PickingFreeSpectralLikelihood:
    def __init__(
        self, fuser, param_space, d_obs, config, velocity_max=None, channel="fused"
    ):
        self.fuser, self.param_space, self.config = fuser, param_space, config
        c = np.asarray(fuser.c_grid)
        f = np.asarray(fuser.freqs_broadband)
        if (
            c.ndim != 1
            or len(c) < 2
            or not np.all(np.isfinite(c))
            or np.any(c <= 0)
            or np.any(np.diff(c) <= 0)
        ):
            raise ValueError("positive finite increasing velocity grid required")
        if (
            f.ndim != 1
            or len(f) < 2
            or np.any(f <= 0)
            or not np.all(np.isfinite(f))
            or np.any(np.diff(f) <= 0)
        ):
            raise ValueError(
                "finite increasing frequency grid with at least two nodes required"
            )
        if np.iscomplexobj(d_obs):
            raise ValueError("observed spectra must be real-valued")
        self.d_obs = np.array(d_obs, dtype=float, order="C", copy=True)
        if (
            self.d_obs.shape != (len(f), len(c))
            or not np.all(np.isfinite(self.d_obs))
            or np.any(self.d_obs < 0)
        ):
            raise ValueError(
                "observations must be finite nonnegative spectra on the declared grids"
            )
        if not np.isfinite(config.sigma_likelihood) or config.sigma_likelihood <= 0:
            raise ValueError("positive finite likelihood sigma required")
        if (
            not np.isfinite(config.wasserstein_weight)
            or not 0 <= config.wasserstein_weight <= 1
        ):
            raise ValueError("Wasserstein weight must be in [0,1]")
        self.wasserstein_scale_mps = float(
            config.wasserstein_scale_mps
            if config.wasserstein_scale_mps is not None
            else 0.15 * (c[-1] - c[0])
        )
        if (
            not np.isfinite(self.wasserstein_scale_mps)
            or self.wasserstein_scale_mps <= 0
        ):
            raise ValueError("positive finite Wasserstein scale required")
        self.eval_count = 0
        self.velocity_max = velocity_max
        if velocity_max is not None and (
            not np.isfinite(velocity_max) or np.sum(c <= velocity_max) < 2
        ):
            raise ValueError("velocity window must contain at least two grid nodes")
        if channel not in ("fused", "active_bb", "passive_bb"):
            raise ValueError("invalid spectrum channel")
        self.channel = channel
        self._channel_support = np.ones(len(f), bool)
        if channel != "fused":
            source = (
                fuser.operator.freqs_active
                if channel == "active_bb"
                else fuser.operator.freqs_passive
            )
            self._channel_support = np.any(
                np.isclose(f[:, None], source[None, :], rtol=0.0, atol=1e-9), axis=1
            )
        self.set_frequency_weights(np.ones(len(f)))

    def set_frequency_weights(self, importance):
        """Integrate each retained contiguous band on its own support/endpoints."""
        from forward_mfj import compute_trapezoidal_weights

        importance = np.asarray(importance, dtype=float)
        f = self.fuser.freqs_broadband
        if (
            importance.shape != f.shape
            or not np.all(np.isfinite(importance))
            or np.any(importance < 0)
        ):
            raise ValueError(
                "frequency importance must be finite, nonnegative and match frequency nodes"
            )
        selected = (importance > 0) & self._channel_support
        edges = np.diff(np.r_[False, selected, False].astype(int))
        weights = np.zeros(len(f))
        for start, end in zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)):
            if end - start < 2:
                raise ValueError(
                    "each retained frequency band needs at least two nodes"
                )
            weights[start:end] = (
                compute_trapezoidal_weights(f[start:end]) * importance[start:end]
            )
        if not np.any(weights > 0):
            raise ValueError("at least one nonempty frequency band required")
        self.freq_weights = weights

    def spectral_metrics(self, syn):
        mask = (
            np.ones(len(self.fuser.c_grid), bool)
            if self.velocity_max is None
            else self.fuser.c_grid <= self.velocity_max
        )
        return compute_spectral_misfit_numba(
            self.d_obs[:, mask],
            syn[:, mask],
            self.fuser.c_grid[mask],
            self.freq_weights,
            self.config.wasserstein_weight,
            2.0,
            self.fuser.freqs_broadband,
            self.config.passive_cutoff_hz,
            self.config.passive_velocity_max,
            self.wasserstein_scale_mps,
        )

    def evaluate_physical_model(self, m_phys):
        self.eval_count += 1
        model = self.param_space.vector_to_model(m_phys)
        syn = self.fuser.generate_fused_spectrum_for_model(model, add_noise=False)[
            "d_" + self.channel
        ]
        jt, jn, jw = self.spectral_metrics(syn)
        return -0.5 * jt / self.config.sigma_likelihood**2, jt, jn, jw, syn

    def log_posterior_latent(self, z_latent):
        if not np.all(np.isfinite(z_latent)):
            return -np.inf, np.inf, np.inf, np.inf
        try:
            m = self.param_space.latent_to_physical(z_latent)
            ll, jt, jn, jw, _ = self.evaluate_physical_model(m)
        except ValueError:
            return -np.inf, np.inf, np.inf, np.inf
        # Uniform prior on the admissible physical domain; exact triangular Jacobian.
        return (
            ll + float(self.param_space.latent_log_abs_det_jacobian(z_latent)),
            jt,
            jn,
            jw,
        )

    def score_and_post_latent(self, z_latent, fd_eps=0.02):
        if self.config.gradient_backend == "autodiff":
            if not hasattr(self, "_ad_target"):
                from autodiff_likelihood import AutodiffSpectralLikelihood

                self._ad_target = AutodiffSpectralLikelihood(
                    self.fuser,
                    self.param_space,
                    self.d_obs,
                    self.config,
                    self.velocity_max,
                    self.channel,
                )
            self._ad_target.freq_weights = self.freq_weights.copy()
            self.eval_count += 1
            return self._ad_target.score_and_post_latent(z_latent)
        if self.config.gradient_backend != "finite_difference":
            raise ValueError("Unknown gradient backend")
        lp, jt, jn, jw = self.log_posterior_latent(z_latent)
        grad = np.empty(len(z_latent))
        for i in range(len(z_latent)):
            dz = np.zeros(len(z_latent))
            dz[i] = fd_eps
            grad[i] = (
                self.log_posterior_latent(z_latent + dz)[0]
                - self.log_posterior_latent(z_latent - dz)[0]
            ) / (2 * fd_eps)
        if not np.all(np.isfinite(grad)):
            raise FloatingPointError("nonfinite target gradient")
        return grad, lp, jt, jn, jw


class SteinMixtureVariationalInference:
    def __init__(self, likelihood, config=None):
        self.likelihood = likelihood
        self.param_space = likelihood.param_space
        self.config = config or likelihood.config
        self.dim, self.K = self.param_space.dim, self.config.n_components
        if min(self.K, self.config.n_iterations, self.config.n_mc_samples_per_comp) < 1:
            raise ValueError("positive components, iterations and MC count required")
        if self.config.entropy_weight != 1:
            raise ValueError("ELBO requires entropy_weight=1")
        self.rng = np.random.default_rng(self.config.random_seed)
        self.weights = np.ones(self.K) / self.K
        self.means_z = np.zeros((self.K, self.dim))
        self.log_stds_z = np.full((self.K, self.dim), -0.7)
        self.offdiag_z = np.zeros((self.K, self.dim, self.dim))

    def _cholesky(self, k):
        return np.tril(self.offdiag_z[k], -1) + np.diag(np.exp(self.log_stds_z[k]))

    def _initialize_mixture_components(self):
        count = max(24, self.K * 8)
        u = np.empty((count, self.dim))
        for j in range(self.dim):
            u[:, j] = (
                self.rng.permutation(count) + self.rng.uniform(0.01, 0.99, count)
            ) / count
        # Latin hypercube in sigmoid coordinates, independent of true model values.
        candidates = np.log(u) - np.log1p(-u)
        scores = np.array(
            [self.likelihood.log_posterior_latent(z)[0] for z in candidates]
        )
        selected = []
        order = np.argsort(scores)[::-1]
        for idx in order:
            if (
                not selected
                or min(
                    np.linalg.norm(candidates[idx] - candidates[j]) for j in selected
                )
                > 0.5
            ):
                selected.append(idx)
            if len(selected) == self.K:
                break
        selected += [i for i in order if i not in selected][: self.K - len(selected)]
        self.means_z = candidates[selected].copy()
        self.initialization_records = []
        if self.config.initialization_refinement_steps:
            from scipy.optimize import minimize

            for k in range(self.K):

                def objective(z):
                    grad, lp, *_ = self.likelihood.score_and_post_latent(z)
                    return -lp, -grad

                opt = minimize(
                    objective,
                    self.means_z[k],
                    jac=True,
                    method="L-BFGS-B",
                    options={
                        "maxiter": self.config.initialization_refinement_steps,
                        "maxls": 12,
                    },
                )
                if np.isfinite(opt.fun):
                    self.means_z[k] = opt.x
                self.initialization_records.append(
                    {
                        "success": bool(opt.success),
                        "iterations": int(opt.nit),
                        "objective": float(opt.fun),
                    }
                )
        return selected

    def _mixture_details(self, z):
        logs, scores = [], []
        for k in range(self.K):
            L = self._cholesky(k)
            y = solve_triangular(L, z - self.means_z[k], lower=True)
            logs.append(
                np.log(self.weights[k])
                - 0.5 * self.dim * np.log(2 * np.pi)
                - np.log(np.diag(L)).sum()
                - 0.5 * y @ y
            )
            scores.append(-solve_triangular(L.T, y, lower=False))
        lq = float(logsumexp(logs))
        responsibilities = np.exp(np.array(logs) - lq)
        return lq, responsibilities @ np.asarray(scores), responsibilities

    def _compute_mixture_log_prob_and_score(self, z):
        lq, score, _ = self._mixture_details(z)
        return lq, score

    def _stein_kernel_update(self, particles, forces):
        """SMVI transport over ALL free component parameters, equal weights.

        Theta=(mu, log diagonal L, strict lower L). With unweighted component
        ELBO gradients, attraction and entropy repulsion both carry 1/K.
        RBF bandwidth h=median squared distance/log(K+1).
        """
        if self.K == 1:
            return forces.copy()
        diff = particles[:, None, :] - particles[None, :, :]
        ds = np.sum(diff**2, axis=-1)
        positive = ds[np.triu_indices(self.K, 1)]
        h = max(float(np.median(positive)) / np.log(self.K + 1), 1e-8)
        kernel = np.exp(-ds / h)
        attraction = kernel.T @ forces / self.K
        repulsion = -2 * np.sum(kernel[:, :, None] * diff, axis=0) / (h * self.K)
        return attraction + self.config.repulsion_strength * repulsion

    def run(self, verbose=True):
        start = time.time()
        self._initialize_mixture_components()
        hist = {
            key: []
            for key in ["elbo", "misfit_mean", "misfit_min", "ncc", "w2", "means_phys"]
        }
        adam = {}

        def update(name, value, gradient, rate, iteration):
            if name not in adam:
                adam[name] = (np.zeros_like(value), np.zeros_like(value))
            m, v = adam[name]
            gradient = np.clip(gradient, -100.0, 100.0)
            m[:] = 0.9 * m + 0.1 * gradient
            v[:] = 0.999 * v + 0.001 * gradient**2
            return value + rate * (m / (1 - 0.9**iteration)) / (
                np.sqrt(v / (1 - 0.999**iteration)) + 1e-8
            )

        for it in range(1, self.config.n_iterations + 1):
            gm = np.zeros_like(self.means_z)
            gl = np.zeros_like(self.offdiag_z)
            elbos = np.zeros(self.K)
            metrics = np.zeros((self.K, 3))
            eps_shared = self.rng.normal(
                size=(self.config.n_mc_samples_per_comp, self.dim)
            )
            for k in range(self.K):
                L = self._cholesky(k)
                for eps in eps_shared:
                    z = self.means_z[k] + L @ eps
                    g, lp, jt, jn, jw = self.likelihood.score_and_post_latent(z)
                    lq, gq = self._compute_mixture_log_prob_and_score(z)
                    residual = g - gq
                    gm[k] += residual / self.config.n_mc_samples_per_comp
                    gl[k] += np.outer(residual, eps) / self.config.n_mc_samples_per_comp
                    elbos[k] += (lp - lq) / self.config.n_mc_samples_per_comp
                    metrics[k] += (
                        np.array([jt, jn, jw]) / self.config.n_mc_samples_per_comp
                    )
            gdiag = np.diagonal(gl, axis1=1, axis2=2) * np.exp(self.log_stds_z)
            lower = np.tril_indices(self.dim, -1)
            particles = np.concatenate(
                (self.means_z, self.log_stds_z, self.offdiag_z[:, lower[0], lower[1]]),
                axis=1,
            )
            forces = np.concatenate((gm, gdiag, gl[:, lower[0], lower[1]]), axis=1)
            transport = self._stein_kernel_update(particles, forces)
            self.means_z = update(
                "means", self.means_z, transport[:, : self.dim], self.config.lr_mean, it
            )
            self.log_stds_z = update(
                "diag",
                self.log_stds_z,
                transport[:, self.dim : 2 * self.dim],
                self.config.lr_chol,
                it,
            )
            off = update(
                "offdiag",
                self.offdiag_z[:, lower[0], lower[1]],
                transport[:, 2 * self.dim :],
                self.config.lr_chol,
                it,
            )
            self.offdiag_z[:, lower[0], lower[1]] = off
            self.log_stds_z = np.clip(self.log_stds_z, -6.0, 2.0)
            # The SMVI mixture has equal weights; no ad hoc weight dynamics.
            self.weights[:] = 1 / self.K
            best = np.argmin(metrics[:, 0])
            for key, val in [
                ("elbo", self.weights @ elbos),
                ("misfit_mean", self.weights @ metrics[:, 0]),
                ("misfit_min", metrics[best, 0]),
                ("ncc", metrics[best, 1]),
                ("w2", metrics[best, 2]),
                ("means_phys", self.param_space.latent_to_physical(self.means_z)),
            ]:
                hist[key].append(val)
            if verbose and (it == 1 or it % 10 == 0 or it == self.config.n_iterations):
                print(
                    f'  [SMVI {it}/{self.config.n_iterations}] ELBO={hist["elbo"][-1]:.3f}, misfit={hist["misfit_min"][-1]:.4f}',
                    flush=True,
                )
        zs, ms, labels = self.sample_posterior(self.config.n_posterior_draws)
        center_metrics = [self.likelihood.log_posterior_latent(z) for z in self.means_z]
        # A true global MAP is not inferred from a finite set of component centers.
        best = int(np.argmin([v[1] for v in center_metrics]))
        best_phys = self.param_space.latent_to_physical(self.means_z[best])
        mean = ms.mean(axis=0)
        _, jm, nm, wm, ds_mean = self.likelihood.evaluate_physical_model(mean)
        _, jb, nb, wb, ds_best = self.likelihood.evaluate_physical_model(best_phys)
        result = {
            "weights": self.weights.copy(),
            "means_z": self.means_z.copy(),
            "log_stds_z": self.log_stds_z.copy(),
            "cholesky_z": np.array([self._cholesky(k) for k in range(self.K)]),
            "comp_means_phys": self.param_space.latent_to_physical(self.means_z),
            "comp_misfits": np.array([v[1] for v in center_metrics]),
            "posterior_samples_z": zs,
            "posterior_samples_phys": ms,
            "posterior_comp_labels": labels,
            "posterior_mean_phys": mean,
            "best_component_phys": best_phys,
            "d_syn_posterior_mean": ds_mean,
            "d_syn_best_component": ds_best,
            "mean_model_metrics": {
                "misfit": jm,
                "ncc_similarity": 1 - nm,
                "wasserstein_w2_mps": wm,
            },
            "best_component_metrics": {
                "misfit": jb,
                "ncc_similarity": 1 - nb,
                "wasserstein_w2_mps": wb,
            },
            "elapsed_sec": np.array([time.time() - start]),
            "n_forward_evals": np.array([self.likelihood.eval_count]),
        }
        if self.config.predictive_draws > 0 and hasattr(self.likelihood, "fuser"):
            from posterior_prediction import evaluate_mixture_prediction

            result.update(
                evaluate_mixture_prediction(
                    self.likelihood,
                    self.means_z,
                    np.array([self._cholesky(k) for k in range(self.K)]),
                    self.weights,
                    self.config,
                )
            )
        result["optimization_diagnostics"] = {
            "converged": False,
            "status": "requires_independent_restarts_and_budget_study",
            "initialization_refinement": self.initialization_records,
        }
        result["elapsed_sec"] = np.array([time.time() - start])
        result["n_forward_evals"] = np.array([self.likelihood.eval_count])
        result.update({"history_" + k: np.asarray(v) for k, v in hist.items()})
        return result

    def sample_posterior(self, n_samples=600):
        labels = self.rng.choice(self.K, n_samples, p=self.weights)
        z = np.array(
            [
                self.means_z[k] + self._cholesky(k) @ self.rng.normal(size=self.dim)
                for k in labels
            ]
        )
        return z, self.param_space.latent_to_physical(z), labels


class AffineInvariantEnsembleMCMC:
    def __init__(
        self,
        likelihood,
        n_walkers=24,
        n_steps=1000,
        burn_in=300,
        stretch_a=2.0,
        seed=2026,
    ):
        self.likelihood = likelihood
        self.param_space = likelihood.param_space
        self.dim = self.param_space.dim
        self.n_walkers = max(n_walkers, 2 * self.dim + 2)
        self.n_steps = n_steps
        self.burn_in = burn_in
        self.stretch_a = stretch_a
        self.rng = np.random.default_rng(seed)
        if not 0 <= burn_in < n_steps or stretch_a <= 1:
            raise ValueError("invalid burn-in or stretch factor")

    def run(self, init_centers_z=None):
        start = time.time()
        if init_centers_z is None:
            z = self.rng.normal(0, 1, (self.n_walkers, self.dim))
        else:
            z = np.array(
                [
                    init_centers_z[i % len(init_centers_z)]
                    + self.rng.normal(0, 0.2, self.dim)
                    for i in range(self.n_walkers)
                ]
            )
        lp = np.array([self.likelihood.log_posterior_latent(v)[0] for v in z])
        chain = np.empty((self.n_steps, self.n_walkers, self.dim))
        accepted = 0
        for step in range(self.n_steps):
            for w in range(self.n_walkers):
                j = self.rng.integers(self.n_walkers - 1)
                j += j >= w
                stretch = (
                    (self.stretch_a - 1) * self.rng.uniform() + 1
                ) ** 2 / self.stretch_a
                proposal = z[j] + stretch * (z[w] - z[j])
                newlp = self.likelihood.log_posterior_latent(proposal)[0]
                if (
                    np.isfinite(newlp)
                    and np.log(self.rng.uniform())
                    < (self.dim - 1) * np.log(stretch) + newlp - lp[w]
                ):
                    z[w] = proposal
                    lp[w] = newlp
                    accepted += 1
            chain[step] = z
        post = chain[self.burn_in :]
        # Per-walker ACF, averaged over walkers. Short estimates are explicitly untrusted.
        taus = []
        for d in range(self.dim):
            x = post[:, :, d] - post[:, :, d].mean(axis=0)
            n = len(x)
            ft = np.fft.rfft(x, n=2 * n, axis=0)
            ac = np.fft.irfft(ft * np.conj(ft), axis=0)[:n].mean(axis=1)
            ac = ac / max(ac[0], 1e-30)
            tau = 1.0
            for lag in range(1, n):
                if ac[lag] < 0:
                    break
                tau += 2 * ac[lag]
                if lag > 5 * tau:
                    break
            taus.append(max(tau, 1.0))
        tau = np.array(taus)
        reliable = bool(len(post) >= 50 * np.max(tau))
        samples = post.reshape(-1, self.dim)
        return {
            "samples_z": samples,
            "samples_phys": self.param_space.latent_to_physical(samples),
            "chain_z": chain,
            "acceptance_rate": np.array([accepted / (self.n_steps * self.n_walkers)]),
            "autocorrelation_time_steps": tau,
            "ess": len(samples) / tau,
            "convergence_status": (
                "acf_length_sufficient_unconfirmed"
                if reliable
                else "insufficient_chain_length"
            ),
            "converged": False,  # independent starts/mode exploration still need examination
            "elapsed_sec": np.array([time.time() - start]),
            "n_evals": np.array([self.n_walkers * (self.n_steps + 1)]),
        }
