# -*- coding: utf-8 -*-
"""
Model C core:
Physics model + small-sample calibration + GPR residual correction.

This version intentionally keeps the original rheology equations as much as possible,
but replaces the hand-written binary loop with:
    1) logarithmic scan for a valid bracket;
    2) Brent root finding if a root exists;
    3) bounded minimization / boundary fallback if no root exists.

The physical model is calibrated by df and k_phi:
    phi_p_eff = k_phi * phi_cb

Then a Gaussian Process learns the remaining residual:
    log10(eta_exp) = log10(eta_physics) + GP(log10(gamma))

With only one calibration formulation, the GP residual is assumed to be a
systematic shear-rate-dependent correction shared by nearby formulations.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Tuple

import numpy as np
from scipy.optimize import brentq, differential_evolution, minimize_scalar, least_squares
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, RBF

# ----------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------
PI = 3.141592653
kb = 1.38e-23
TT = 298.15
hp = 6.626e-34
v0 = 4e15
eps0 = 8.854e-12


# ----------------------------------------------------------------------
# Original auxiliary functions
# ----------------------------------------------------------------------
def cal_dlvo(
    aa=75e-9,
    eps_cb=2.78,
    eps_nmp=32.2,
    n_cb=2.0,
    n_nmp=1.47,
    sigma=6e-9,
    nb=12,
    lmd=1,
    zeta=10e-3,
    verbose=False,
):
    """Vectorized extended-DLVO calculation."""
    ah1 = (3 / 4) * kb * TT * ((eps_cb - eps_nmp) / (eps_cb + eps_nmp)) ** 2
    ah2 = (
        3
        / (16 * math.sqrt(2))
        * hp
        * v0
        * (n_cb**2 - n_nmp**2) ** 2
        / (n_cb**2 + n_nmp**2) ** (3 / 2)
    )
    Ah = ah1 + ah2

    dx = 2e-12
    hh = np.arange(1, 5001, dtype=float) * dx
    rc = 2 + hh / aa

    vdw1 = 2 * aa**2 / (hh * (4 * aa + hh))
    vdw2 = 2 * aa**2 / (2 * aa + hh) ** 2
    vdw3 = np.log((hh * (4 * aa + hh)) / ((2 * aa + hh) ** 2))
    V_vdw = -(Ah / 6) * (vdw1 + vdw2 + vdw3)

    V_ele = (4 * PI * eps0 * eps_nmp * zeta**2 * aa**2) / (hh + 2 * aa)

    born1 = 4 * Ah * (sigma / aa) ** (nb - 6)
    born2 = math.factorial(nb - 6) / math.factorial(nb - 2)
    born3 = lmd / ((lmd + 1) * (rc - 1 - lmd) ** (nb - 5))
    V_born = born1 * born2 * born3

    V_total = V_vdw + V_ele + V_born

    if verbose:
        for i in range(0, len(hh), 100):
            print(
                f"hh={hh[i]/1e-9:.3f} nm | "
                f"vdW={V_vdw[i]/(kb*TT):.3f} kBT | "
                f"ele={V_ele[i]/(kb*TT):.3f} kBT | "
                f"born={V_born[i]/(kb*TT):.3f} kBT | "
                f"total={V_total[i]/(kb*TT):.3f} kBT"
            )

    return hh, V_vdw, V_born, V_ele, V_total


def cal_param(list_x, list_y, dx=2e-12, verbose=False):
    """
    Keep the definition close to the original implementation so that the
    Figure-5C reproduction remains directly comparable.
    """
    x = np.asarray(list_x, dtype=float)
    y = np.asarray(list_y, dtype=float)

    dy = np.gradient(y, x)
    d2y = np.gradient(dy, x)

    idx_uc = int(np.argmin(y))
    Uc = -float(y[idx_uc])
    hc = float(x[idx_uc])

    idx_fc = int(np.argmax(dy))
    Fc = float(dy[idx_fc])

    if not np.isfinite(Fc) or Fc <= 0:
        raise ValueError(f"Invalid Fc={Fc}. Check the DLVO curve.")

    k0 = float(d2y[idx_uc] * hc / Fc)
    if not np.isfinite(k0) or k0 <= 0:
        raise ValueError(f"Invalid k0={k0}. Check the DLVO derivative.")

    if verbose:
        print(
            f"Uc={Uc/(kb*TT):.3f} kBT | hc={hc/1e-9:.3f} nm | "
            f"Fc={Fc/1e-12:.3f} pN | k0={k0:.3f}"
        )

    return Uc, hc, Fc, k0


def cal_miu(w_pvdf, eta_sol=1.65e-3, mode="5130"):
    if mode in ["5130", "SOLEF5130", "Solef5130"]:
        k_eta = 2.11e6
    elif mode in ["6020", "SOLEF6020", "Solef6020"]:
        k_eta = 1.18e6
    elif mode in ["1015", "SOLEF1015", "Solef1015"]:
        k_eta = 4.95e6
    else:
        raise ValueError(
            f"Unsupported mode: {mode}. Supported modes: 5130, 6020, 1015."
        )
    eta_sp = k_eta * w_pvdf**3
    return eta_sp * eta_sol + eta_sol


# ----------------------------------------------------------------------
# Robust physical rheology solver
# ----------------------------------------------------------------------
@dataclass
class RheologyPhysicsParameters:
    hc: float
    uc: float
    fc: float
    k0: float
    aa: float = 75e-9
    dl: float = 1.55
    miu: float = 3.6e-3
    phi_m: float = 0.6
    zz: float = 5.0
    phi_ref: float = 1.2
    q_min_ratio: float = 1.0
    q_max_ratio: float = 30.0


def _validate_physics(p: RheologyPhysicsParameters, df: float, phi_p: float):
    if p.hc <= 0 or p.uc <= 0 or p.fc <= 0 or p.k0 <= 0 or p.aa <= 0:
        raise ValueError("hc, uc, fc, k0 and aa must all be > 0.")
    if p.miu <= 0:
        raise ValueError("miu must be > 0.")
    if not (1.0 < df < 3.0):
        raise ValueError(f"df={df} invalid; require 1 < df < 3.")
    if p.dl <= 0:
        raise ValueError("dl must be > 0.")
    if not (0 < p.phi_m < 1):
        raise ValueError("phi_m must be between 0 and 1.")
    if phi_p <= 0:
        raise ValueError("effective phi_p must be > 0.")


def _state_at_q(
    q: float,
    gamma: float,
    phi_p: float,
    df: float,
    p: RheologyPhysicsParameters,
) -> Dict[str, float]:
    """Original rheology equations evaluated at one q and gamma."""
    phi_int = (q / p.aa) ** (df - 3.0)
    phi_a = phi_p / phi_int

    # Keep the original clipping rule for strict comparability with old code.
    phi_a_cap = p.phi_m * 0.95 / (p.phi_ref**3)
    phi_a = min(phi_a, phi_a_cap)

    tau_0 = 6.0 * PI * p.miu * p.aa**3 / (kb * TT)

    log_u = math.log(p.uc / (kb * TT))
    log_tau_0 = math.log(tau_0)
    log_k0 = math.log(p.k0)
    log_aa_hc = math.log(p.aa / p.hc)
    log_qa = math.log(q / p.aa)

    # tau_si from the original tau_sic expression
    log_tau_sic = (
        -log_tau_0
        + 0.5 * log_k0
        + math.log((p.uc * p.aa) / (kb * TT * p.hc))
        - p.zz * p.uc / (kb * TT)
    )
    log_tau_si = -log_tau_sic

    # log-domain time constants
    log_tau_sc = (
        -log_tau_0
        + 0.5 * log_k0
        + 1.5 * log_u
        + 2.0 * log_aa_hc
        + math.log(p.aa / q)
        - p.uc / (kb * TT)
    )
    log_tau_s = -log_tau_sc

    log_tau_rc = (
        0.5 * log_k0
        + math.log(gamma)
        + log_aa_hc
        - p.dl * log_qa
    )
    log_tau_r = -log_tau_rc
    log_tau_sr = log_tau_0 + log_qa
    log_tau_rs = log_tau_si - p.dl * log_qa

    log_r_plus_rs = np.logaddexp(log_tau_r, log_tau_rs)
    log_num = log_tau_s + log_r_plus_rs
    log_den = np.logaddexp(
        log_tau_s + log_tau_rs,
        log_r_plus_rs + log_tau_sr,
    )
    log_alpha = float(log_num - log_den)
    alpha = float(np.exp(np.clip(log_alpha, -700, 700)))

    base = 1.0 - (phi_a / p.phi_m) * p.phi_ref**3
    base = max(base, 1e-12)
    eta1 = p.miu * base ** (-2.0)

    exponent = (5.0 - 2.0 * df - p.dl) / (3.0 - df)

    # eta2 in log-domain
    log_eta2 = (
        log_alpha
        + math.log(p.fc)
        + math.log(p.hc)
        - 0.5 * log_k0
        + 2.0 * math.log(phi_a)
        - 3.0 * math.log(p.aa)
        - math.log(gamma)
        + exponent * math.log(phi_p / phi_a)
    )

    log_eta = float(np.logaddexp(math.log(eta1), log_eta2))
    eta = float(np.exp(np.clip(log_eta, -700, 700)))

    log_gamma_c = (
        math.log(2.0 * p.fc * p.aa)
        - math.log(5.0 * PI)
        - log_eta
        - 3.0 * math.log(q)
    )
    gamma_c = float(np.exp(np.clip(log_gamma_c, -700, 700)))
    residual = float(log_gamma_c - math.log(gamma))

    return {
        "eta": eta,
        "eta1": eta1,
        "eta2": float(np.exp(np.clip(log_eta2, -700, 700))),
        "phi_a": phi_a,
        "alpha": alpha,
        "gamma_c": gamma_c,
        "residual": residual,
        "q": q,
        "q_ratio": q / p.aa,
    }


def solve_rheology_curve(
    gamma_array: Iterable[float],
    phi_cb: float,
    df: float,
    k_phi: float,
    physics: RheologyPhysicsParameters,
    n_scan: int = 120,
) -> Tuple[np.ndarray, np.ndarray, list]:
    """
    Predict a rheology curve for raw CB volume fraction phi_cb.

    effective particle fraction used by the original model:
        phi_p = k_phi * phi_cb
    """
    gamma_array = np.asarray(gamma_array, dtype=float)
    if np.any(gamma_array <= 0):
        raise ValueError("All shear rates must be > 0.")

    phi_p = float(k_phi * phi_cb)
    _validate_physics(physics, df, phi_p)

    q_grid = physics.aa * np.geomspace(
        physics.q_min_ratio, physics.q_max_ratio, n_scan
    )

    eta_out = []
    details = []

    for gamma in gamma_array:
        f_grid = np.empty_like(q_grid)

        for j, q in enumerate(q_grid):
            f_grid[j] = _state_at_q(q, gamma, phi_p, df, physics)["residual"]

        finite = np.isfinite(f_grid)
        bracket = None

        for j in range(len(q_grid) - 1):
            if not (finite[j] and finite[j + 1]):
                continue
            if f_grid[j] == 0:
                bracket = (q_grid[j], q_grid[j])
                break
            if f_grid[j] * f_grid[j + 1] < 0:
                bracket = (q_grid[j], q_grid[j + 1])
                break

        if bracket is not None:
            q1, q2 = bracket
            if q1 == q2:
                q_star = q1
            else:
                def root_fun(q):
                    return _state_at_q(q, gamma, phi_p, df, physics)["residual"]

                q_star = brentq(
                    root_fun,
                    q1,
                    q2,
                    xtol=physics.aa * 1e-10,
                    rtol=1e-10,
                    maxiter=100,
                )
            status = "root"
        else:
            # No mathematical root in the allowed q-domain.
            # Find the physically closest self-consistent state.
            def objective(log_q):
                q = math.exp(log_q)
                r = _state_at_q(q, gamma, phi_p, df, physics)["residual"]
                return r * r if np.isfinite(r) else 1e100

            res = minimize_scalar(
                objective,
                bounds=(
                    math.log(physics.aa * physics.q_min_ratio),
                    math.log(physics.aa * physics.q_max_ratio),
                ),
                method="bounded",
                options={"xatol": 1e-12},
            )
            q_star = math.exp(res.x)

            # Explicitly compare both physical boundaries.
            candidates = [
                physics.aa * physics.q_min_ratio,
                q_star,
                physics.aa * physics.q_max_ratio,
            ]
            candidate_states = [
                _state_at_q(q, gamma, phi_p, df, physics) for q in candidates
            ]
            idx = int(np.argmin([abs(s["residual"]) for s in candidate_states]))
            q_star = candidates[idx]

            if idx == 0:
                status = "lower_bound"
            elif idx == 2:
                status = "upper_bound"
            else:
                status = "no_root"

        state = _state_at_q(q_star, gamma, phi_p, df, physics)
        state["status"] = status
        eta_out.append(state["eta"])
        details.append(state)

    return gamma_array, np.asarray(eta_out, dtype=float), details


# ----------------------------------------------------------------------
# Model C
# ----------------------------------------------------------------------
class ModelC:
    """
    Physics + calibration + Gaussian-process residual.

    One-formulation demo:
      - calibrate df and k_phi using one experimental curve;
      - GP learns shear-rate-dependent model discrepancy;
      - correction is transferred to another formulation.

    This is intentionally conservative: with only one formulation used for
    calibration, the GP is NOT allowed to learn composition dependence.
    Composition dependence comes only from the physical model.
    """

    def __init__(
        self,
        physics: RheologyPhysicsParameters,
        df_prior: float = 1.8,
        k_phi_prior: float = 3.20345,
        df_bounds=(1.20, 2.80),
        k_phi_bounds=(0.5, 8.0),
        regularization=0.01,
        random_state=42,
    ):
        self.physics = physics
        self.df_prior = float(df_prior)
        self.k_phi_prior = float(k_phi_prior)
        self.df_bounds = tuple(df_bounds)
        self.k_phi_bounds = tuple(k_phi_bounds)
        self.regularization = float(regularization)
        self.random_state = int(random_state)

        self.df_ = None
        self.k_phi_ = None
        self.gp_ = None
        self.calibration_phi_cb_ = None
        self.train_gamma_ = None
        self.train_eta_ = None

    @staticmethod
    def _clean_curve(gamma, eta):
        gamma = np.asarray(gamma, dtype=float)
        eta = np.asarray(eta, dtype=float)
        mask = (
            np.isfinite(gamma)
            & np.isfinite(eta)
            & (gamma > 0)
            & (eta > 0)
        )
        gamma = gamma[mask]
        eta = eta[mask]
        order = np.argsort(gamma)
        return gamma[order], eta[order]

    def _physics_log_residual(self, theta, gamma, eta_exp, phi_cb):
        df, k_phi = theta
        _, eta_phy, _ = solve_rheology_curve(
            gamma,
            phi_cb=phi_cb,
            df=df,
            k_phi=k_phi,
            physics=self.physics,
        )

        residual = np.log10(eta_phy) - np.log10(eta_exp)

        # Weak regularization toward the original model.
        # It improves identifiability when only one formulation is available.
        reg = self.regularization * np.array([
            (df - self.df_prior) / 0.5,
            (math.log(k_phi) - math.log(self.k_phi_prior)) / 0.8,
        ])
        return np.r_[residual, reg]

    def calibrate(self, phi_cb, gamma_exp, eta_exp, verbose=True):
        gamma_exp, eta_exp = self._clean_curve(gamma_exp, eta_exp)

        if len(gamma_exp) < 4:
            raise ValueError("At least 4 valid rheology points are required.")

        bounds = [self.df_bounds, self.k_phi_bounds]

        # Global search first, to avoid poor local minima.
        def scalar_loss(theta):
            r = self._physics_log_residual(theta, gamma_exp, eta_exp, phi_cb)
            return float(np.mean(r**2))

        de = differential_evolution(
            scalar_loss,
            bounds=bounds,
            seed=self.random_state,
            polish=False,
            popsize=10,
            maxiter=35,
            tol=1e-6,
            updating="immediate",
        )

        # Local least-squares refinement.
        ls = least_squares(
            lambda t: self._physics_log_residual(t, gamma_exp, eta_exp, phi_cb),
            x0=de.x,
            bounds=(
                [self.df_bounds[0], self.k_phi_bounds[0]],
                [self.df_bounds[1], self.k_phi_bounds[1]],
            ),
            xtol=1e-10,
            ftol=1e-10,
            gtol=1e-10,
            max_nfev=300,
        )

        self.df_, self.k_phi_ = map(float, ls.x)
        self.calibration_phi_cb_ = float(phi_cb)
        self.train_gamma_ = gamma_exp.copy()
        self.train_eta_ = eta_exp.copy()

        # Physics prediction after parameter calibration.
        _, eta_phy, _ = solve_rheology_curve(
            gamma_exp,
            phi_cb=phi_cb,
            df=self.df_,
            k_phi=self.k_phi_,
            physics=self.physics,
        )

        # Residual GP in log-space.
        # With one formulation, only gamma is used as GP input.
        X = np.log10(gamma_exp).reshape(-1, 1)
        y = np.log10(eta_exp) - np.log10(eta_phy)

        # Near-interpolating GP: tiny alpha, no white-noise term.
        # RBF length scale is optimized but constrained to avoid wild oscillation.
        kernel = ConstantKernel(
            1.0, (1e-5, 1e3)
        ) * RBF(
            length_scale=0.6,
            length_scale_bounds=(0.03, 5.0),
        )

        self.gp_ = GaussianProcessRegressor(
            kernel=kernel,
            alpha=1e-12,
            normalize_y=False,
            n_restarts_optimizer=8,
            random_state=self.random_state,
        )
        self.gp_.fit(X, y)

        # Diagnostics
        pred = self.predict(phi_cb, gamma_exp, return_std=False)
        rmse_log = float(
            np.sqrt(np.mean((np.log10(pred["eta_model_c"]) - np.log10(eta_exp)) ** 2))
        )
        mape = float(np.mean(np.abs(pred["eta_model_c"] / eta_exp - 1.0)) * 100)

        if verbose:
            print("=" * 72)
            print("Model C calibration")
            print("=" * 72)
            print(f"Calibration phi_CB = {phi_cb:.6f}")
            print(f"df*              = {self.df_:.6f}")
            print(f"k_phi*           = {self.k_phi_:.6f}")
            print(f"GP kernel        = {self.gp_.kernel_}")
            print(f"Training logRMSE = {rmse_log:.6e}")
            print(f"Training MAPE    = {mape:.6e} %")
            print("=" * 72)

        return {
            "df": self.df_,
            "k_phi": self.k_phi_,
            "training_log_rmse": rmse_log,
            "training_mape_percent": mape,
        }

    def predict(self, phi_cb, gamma, return_std=True):
        if self.df_ is None or self.k_phi_ is None or self.gp_ is None:
            raise RuntimeError("Call calibrate() before predict().")

        gamma = np.asarray(gamma, dtype=float)

        gamma, eta_phy, details = solve_rheology_curve(
            gamma,
            phi_cb=phi_cb,
            df=self.df_,
            k_phi=self.k_phi_,
            physics=self.physics,
        )

        X = np.log10(gamma).reshape(-1, 1)

        if return_std:
            delta_mean, delta_std = self.gp_.predict(X, return_std=True)
        else:
            delta_mean = self.gp_.predict(X)
            delta_std = np.zeros_like(delta_mean)

        log_eta_c = np.log10(eta_phy) + delta_mean
        eta_c = 10.0**log_eta_c

        # GP uncertainty in log10-space -> multiplicative band.
        eta_low = 10.0 ** (log_eta_c - 1.96 * delta_std)
        eta_high = 10.0 ** (log_eta_c + 1.96 * delta_std)

        return {
            "gamma": gamma,
            "eta_physics": eta_phy,
            "eta_model_c": eta_c,
            "eta_low95": eta_low,
            "eta_high95": eta_high,
            "gp_delta_log10": delta_mean,
            "gp_std_log10": delta_std,
            "details": details,
        }
