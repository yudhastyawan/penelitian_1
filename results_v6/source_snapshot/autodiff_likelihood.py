"""CPU float64 automatic derivatives of the existing numerical forward model.

This is a vectorized translation of the Numba recursion, not an independent
physics reference. Regression tests compare spectra, target values and gradients
against the Numba implementation. No approximate/surrogate dispersion is used.
"""

import numpy as np
import torch
from torch.nn.functional import softplus
from smvi_inversion import PickingFreeSpectralLikelihood


def _inverse(a, b, c, d):
    det = a * d - b * c
    det = torch.where(abs(det) < 1e-30, det + (1e-30 + 1e-30j), det)
    return d / det, -b / det, -c / det, a / det


def compliance(omega_c, k_val, h_arr, vs_arr, vp_arr, rho_arr, qs_arr, qp_arr):
    n_layers = len(vs_arr)
    k = k_val
    k2 = k * k
    w2 = omega_c * omega_c

    # --- Start at bottom half-space (layer index n_layers - 1) ---
    idx = n_layers - 1
    vs_c = vs_arr[idx] * (1.0 - 0.5j / qs_arr[idx])
    vp_c = vp_arr[idx] * (1.0 - 0.5j / qp_arr[idx])
    mu = rho_arr[idx] * vs_c * vs_c

    nu_p = torch.sqrt(k2 - w2 / (vp_c * vp_c))
    nu_s = torch.sqrt(k2 - w2 / (vs_c * vs_c))
    nu_p = torch.where(nu_p.real < 0, -nu_p, nu_p)
    nu_s = torch.where(nu_s.real < 0, -nu_s, nu_s)

    gamma = 2.0 * k2 - w2 / (vs_c * vs_c)

    # Down-going partition matrices in half-space: K_N = T_d * U_d^{-1}
    ud00, ud01 = k + 0.0j, nu_s
    ud10, ud11 = nu_p, k + 0.0j
    td00, td01 = -2.0 * mu * k * nu_p, -mu * gamma
    td10, td11 = -mu * gamma, -2.0 * mu * k * nu_s

    iud00, iud01, iud10, iud11 = _inverse(ud00, ud01, ud10, ud11)
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

        nu_p = torch.sqrt(k2 - w2 / (vp_c * vp_c))
        nu_s = torch.sqrt(k2 - w2 / (vs_c * vs_c))
        nu_p = torch.where(nu_p.real < 0, -nu_p, nu_p)
        nu_s = torch.where(nu_s.real < 0, -nu_s, nu_s)

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

        ia00, ia01, ia10, ia11 = _inverse(a00, a01, a10, a11)

        # R_j = -A^{-1} * B
        r00 = -(ia00 * b00 + ia01 * b10)
        r01 = -(ia00 * b01 + ia01 * b11)
        r10 = -(ia10 * b00 + ia11 * b10)
        r11 = -(ia10 * b01 + ia11 * b11)

        # Decaying phase factors E_j = diag(ep, es)
        ep = torch.exp(-nu_p * h)
        es = torch.exp(-nu_s * h)

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

        imu00, imu01, imu10, imu11 = _inverse(mu00, mu01, mu10, mu11)

        # K_j = M_T * M_U^{-1}
        k00 = mt00 * imu00 + mt01 * imu10
        k01 = mt00 * imu01 + mt01 * imu11
        k10 = mt10 * imu00 + mt11 * imu10
        k11 = mt10 * imu01 + mt11 * imu11

    # Surface dynamic compliance G(0) = -K_1^{-1}
    ik00, ik01, ik10, ik11 = _inverse(k00, k01, k10, k11)
    g_zz = -ik11

    return g_zz


def normalize(s, power=1.0):
    low = s.amin(dim=1, keepdim=True)
    span = s.amax(dim=1, keepdim=True) - low
    norm = (s - low) / span.clamp_min(1e-30)
    return norm.clamp_min(0.0) ** power


def wasserstein2(p, q, c):
    """Exact discrete quantile transport; piecewise derivatives through CDFs."""
    p = p / p.sum()
    q = q / q.sum()
    cp = p.cumsum(0)
    cq = q.cumsum(0)
    # Exclude terminal sums (mathematically exactly 1) to prevent roundoff
    # exceeding the quantile domain and producing negative interval lengths.
    points = torch.sort(
        torch.cat(
            (
                torch.zeros(1, dtype=p.dtype),
                cp[:-1],
                cq[:-1],
                torch.ones(1, dtype=p.dtype),
            )
        )
    ).values
    mid = (points[1:] + points[:-1]) / 2
    ip = torch.searchsorted(cp.contiguous(), mid.contiguous()).clamp_max(len(c) - 1)
    iq = torch.searchsorted(cq.contiguous(), mid.contiguous()).clamp_max(len(c) - 1)
    return ((points[1:] - points[:-1]) * (c[ip] - c[iq]) ** 2).sum()


class AutodiffSpectralLikelihood(PickingFreeSpectralLikelihood):
    """Same target and quadrature as PickingFreeSpectralLikelihood, AD scores."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        torch.set_num_threads(1)
        op = self.fuser.operator
        self._tensor = lambda a: torch.as_tensor(
            np.asarray(a).copy(), dtype=torch.float64
        )
        self._obs = self._tensor(self.d_obs)
        self._c = self._tensor(op.c_grid)
        self._low = self._tensor(self.param_space.lower_bounds)
        self._high = self._tensor(self.param_space.upper_bounds)
        model = self.param_space.reference_model
        self._rho, self._qs, self._qp = map(
            self._tensor, (model.rho, model.qs, model.qp)
        )
        self._grids = {}
        for mode, fs, r, dr, win, kernel in [
            ("active", op.freqs_active, op.r_act, op.dr_act, op.win_act, op.h0_mfj_act),
            (
                "passive",
                op.freqs_passive,
                op.r_pas,
                op.dr_pas,
                op.win_pas,
                op.h0_mfj_pas,
            ),
        ]:
            nodes = [op.quadrature_grid(f, r) for f in fs]
            omega = 2 * np.pi * fs
            source = (
                (fs / op.active_geom.f_center) ** 2
                * np.exp(-((fs / op.active_geom.f_center) ** 2))
                if mode == "active"
                else 1 / (1 + (fs / 6) ** 1.5)
            )
            indices = np.array(
                [
                    (
                        np.flatnonzero(np.isclose(fs, f, rtol=0, atol=1e-9))[0]
                        if np.any(np.isclose(fs, f, rtol=0, atol=1e-9))
                        else -1
                    )
                    for f in self.fuser.freqs_broadband
                ]
            )
            self._grids[mode] = {
                "k": self._tensor(np.stack([x[0] for x in nodes])),
                "basis": self._tensor(np.stack([x[1] for x in nodes])),
                "omega": self._tensor(omega[:, None]) * (1 + 0.005j),
                "r": self._tensor(r),
                "source": self._tensor(source),
                "weights": self._tensor(np.sqrt(r) * dr * win),
                "kernel": torch.as_tensor(kernel.copy()),
                "indices": torch.as_tensor(indices.clip(0)),
                "available": self._tensor(indices >= 0),
            }
        self._alpha = self._tensor(self.fuser.alpha_aperture)

    def _physical(self, z):
        n = len(self._rho)
        sig = torch.sigmoid(z)
        vs = self._low[:n] + (self._high[:n] - self._low[:n]) * sig[:n]
        vlo = torch.maximum(self._low[n : 2 * n], 1.45 * vs)
        vp = vlo + (self._high[n : 2 * n] - vlo) * sig[n : 2 * n]
        h = (
            self._low[2 * n :]
            + (self._high[2 * n :] - self._low[2 * n :]) * sig[2 * n :]
        )
        widths = torch.cat(
            (
                self._high[:n] - self._low[:n],
                self._high[n : 2 * n] - vlo,
                self._high[2 * n :] - self._low[2 * n :],
            )
        )
        jac = (widths.log() - softplus(-z) - softplus(z)).sum()
        return vs, vp, torch.cat((h, torch.zeros(1, dtype=z.dtype))), jac

    def _spectrum(self, vs, vp, h, mode):
        g = self._grids[mode]
        k = g["k"]
        response = compliance(g["omega"], k, h, vs, vp, self._rho, self._qs, self._qp)
        mu = self._rho[0] * (vs[0] * (1 - 0.5j / self._qs[0])) ** 2
        ratio = (
            vp[0] * (1 - 0.5j / self._qp[0]) / (vs[0] * (1 - 0.5j / self._qs[0]))
        ) ** 2
        a = (1 - (0.5 * ratio - 1) / (ratio - 1)) / mu
        res = response * k - a
        # Real Bessel bases avoid allocating a second complex copy of the cache.
        u = (
            torch.bmm(res.real[:, None, :], g["basis"]).squeeze(1)
            + 1j * torch.bmm(res.imag[:, None, :], g["basis"]).squeeze(1)
            + a / g["r"]
        )
        u = u * g["source"][:, None]
        amp = abs(u)
        u = u / (amp + 0.02 * amp.amax(dim=1, keepdim=True) + 1e-30)
        image = abs(torch.bmm(g["kernel"], (u * g["weights"])[:, :, None]).squeeze(-1))
        image = image / image.amax(dim=1, keepdim=True).clamp_min(1e-30)
        image = normalize(image, self.fuser.config.norm_power)
        return normalize(image[g["indices"]] * g["available"][:, None]), g["available"]

    def spectrum_and_target(self, z):
        vs, vp, h, jac = self._physical(z)
        a, has_a = self._spectrum(vs, vp, h, "active")
        p, has_p = self._spectrum(vs, vp, h, "passive")
        alpha = torch.where(
            (has_a > 0) & (has_p == 0),
            1.0,
            torch.where((has_p > 0) & (has_a == 0), 0.0, self._alpha),
        )
        va = (a.mean(dim=1) / (a.amax(dim=1) + 1e-12)) ** 2
        vp_ = (p.mean(dim=1) / (p.amax(dim=1) + 1e-12)) ** 2
        pa = has_a * alpha / va.clamp_min(1e-4)
        pp = has_p * (1 - alpha) / vp_.clamp_min(1e-4)
        ds = normalize((pa[:, None] * a + pp[:, None] * p) / (pa + pp + 1e-30)[:, None])
        if self.channel == "active_bb":
            ds = a
        if self.channel == "passive_bb":
            ds = p
        jn = []
        jw = []
        weights = []
        mask = (
            np.ones(len(self._c), bool)
            if self.velocity_max is None
            else self.fuser.c_grid <= self.velocity_max
        )
        c = self._c[mask]
        observed = self._obs[:, mask]
        synthetic = ds[:, mask]
        for i, f in enumerate(self.fuser.freqs_broadband):
            if self.freq_weights[i] <= 0:
                continue
            limit = len(c)
            if f <= self.config.passive_cutoff_hz:
                limit = min(
                    limit,
                    max(
                        2,
                        int(
                            np.searchsorted(
                                c.numpy(),
                                self.config.passive_velocity_max,
                                side="right",
                            )
                        ),
                    ),
                )
            o = observed[i, :limit]
            s = synthetic[i, :limit]
            no = torch.linalg.vector_norm(o)
            ns = torch.linalg.vector_norm(s)
            if float(no.detach()) <= 1e-30 and float(ns.detach()) <= 1e-30:
                continue
            jn.append(1 - (o * s).sum() / (no * ns).clamp_min(1e-30))
            jw.append(
                wasserstein2(
                    o.clamp_min(0) ** 2 + 1e-12, s.clamp_min(0) ** 2 + 1e-12, c[:limit]
                )
            )
            weights.append(self.freq_weights[i])
        w = self._tensor(weights)
        jn = (torch.stack(jn) * w).sum() / w.sum()
        jw = (torch.stack(jw) * w).sum() / w.sum()
        loss = (
            (1 - self.config.wasserstein_weight) * jn
            + self.config.wasserstein_weight * jw / self.wasserstein_scale_mps**2
        )
        return (
            ds,
            -0.5 * loss / self.config.sigma_likelihood**2 + jac,
            loss,
            jn,
            torch.sqrt(jw),
        )

    def score_and_post_latent(self, z_latent, fd_eps=None):
        self.eval_count += 1
        z = self._tensor(z_latent).requires_grad_(True)
        ds, lp, jt, jn, jw = self.spectrum_and_target(z)
        (grad,) = torch.autograd.grad(lp, z)
        if not bool(torch.isfinite(grad).all()):
            raise FloatingPointError("nonfinite AD target gradient")
        return (
            grad.detach().numpy(),
            float(lp.detach()),
            float(jt.detach()),
            float(jn.detach()),
            float(jw.detach()),
        )
