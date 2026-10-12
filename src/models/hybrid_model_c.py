"""Physics baseline plus a composition/shear-rate Gaussian-process residual."""

from __future__ import annotations

import copy
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, RBF, WhiteKernel
from sklearn.preprocessing import StandardScaler

from .microrheology import (
    RheologyPhysicsParameters,
    cal_dlvo,
    cal_miu,
    cal_param,
    solve_rheology_curve,
)


def build_base_physics(config: dict) -> RheologyPhysicsParameters:
    material = config["material"]
    w_pvdf = material["w_pvdf"]
    aa = material["primary_radius_m"]
    sigma = (6.0 - 4.0 * w_pvdf**0.15) * 1e-9
    dl = 1.3 + w_pvdf**0.3
    hh, _, _, _, total = cal_dlvo(aa=aa, sigma=sigma)
    uc, hc, fc, k0 = cal_param(hh, total)
    return RheologyPhysicsParameters(
        hc=hc,
        uc=uc,
        fc=fc,
        k0=k0,
        aa=aa,
        dl=dl,
        miu=cal_miu(w_pvdf, mode=material["pvdf_grade"]),
        phi_m=material["maximum_packing_fraction"],
        zz=5.0,
        phi_ref=1.2,
        q_min_ratio=1.0,
        q_max_ratio=30.0,
    )


def _downsample(gamma: np.ndarray, eta: np.ndarray, n: int):
    if len(gamma) <= n:
        return gamma, eta
    index = np.unique(np.round(np.linspace(0, len(gamma) - 1, n)).astype(int))
    return gamma[index], eta[index]


@dataclass
class FitSummary:
    df: float
    k_phi: float
    objective_mse: float


class HybridModelC:
    def __init__(self, physics: RheologyPhysicsParameters, calibration: dict):
        self.physics = copy.deepcopy(physics)
        self.calibration = dict(calibration)
        self.df_ = None
        self.k_phi_ = None
        self.residual_scaler_ = StandardScaler()
        self.residual_gpr_ = None

    def _physics_residual(self, theta, curves):
        df, k_phi = theta
        residuals = []
        n_points = self.calibration["points_per_curve"]
        n_scan = self.calibration["q_scan_calibration"]
        for curve in curves:
            gamma, eta = _downsample(curve["gamma"], curve["eta"], n_points)
            try:
                _, eta_phy, _ = solve_rheology_curve(
                    gamma, curve["phi_cb"], df, k_phi, self.physics, n_scan=n_scan
                )
                residuals.extend((np.log10(eta_phy) - np.log10(eta)) / np.sqrt(len(eta)))
            except (ValueError, OverflowError, FloatingPointError):
                return np.full(sum(min(len(c["gamma"]), n_points) for c in curves) + 2, 1e3)
        residuals.extend([0.02 * (df - 1.8) / 0.5, 0.02 * np.log(k_phi / 3.20345)])
        return np.asarray(residuals)

    def fit(self, curves: list) -> FitSummary:
        train = [curve for curve in curves if curve["split"] == "train"]
        bounds = (
            [self.calibration["df_bounds"][0], self.calibration["k_phi_bounds"][0]],
            [self.calibration["df_bounds"][1], self.calibration["k_phi_bounds"][1]],
        )
        starts = ([1.8, 3.20345], [1.6, 2.0], [2.0, 4.0], [2.2, 6.0], [1.4, 1.0])
        best = None
        for start in starts:
            result = least_squares(
                lambda theta: self._physics_residual(theta, train),
                x0=np.asarray(start), bounds=bounds, max_nfev=80,
                ftol=3e-6, xtol=3e-6, gtol=3e-6,
            )
            mse = float(np.mean(self._physics_residual(result.x, train) ** 2))
            if best is None or mse < best[0]:
                best = (mse, result.x)
        assert best is not None
        self.df_, self.k_phi_ = map(float, best[1])

        features, targets = [], []
        for curve in train:
            _, eta_phy, _ = solve_rheology_curve(
                curve["gamma"], curve["phi_cb"], self.df_, self.k_phi_, self.physics,
                n_scan=self.calibration["q_scan_prediction"],
            )
            features.extend(np.column_stack([
                np.full(len(curve["gamma"]), curve["phi_cb"]),
                np.log10(curve["gamma"]),
            ]))
            targets.extend(np.log10(curve["eta"]) - np.log10(eta_phy))
        x_scaled = self.residual_scaler_.fit_transform(np.asarray(features))
        kernel = (
            ConstantKernel(1.0, (1e-5, 1e3))
            * RBF([2.0, 1.0], [(0.5, 30.0), (0.05, 30.0)])
            + WhiteKernel(1e-5, (1e-9, 5e-2))
        )
        self.residual_gpr_ = GaussianProcessRegressor(
            kernel=kernel, alpha=1e-9, normalize_y=True, n_restarts_optimizer=5,
            random_state=self.calibration["random_state"],
        )
        self.residual_gpr_.fit(x_scaled, np.asarray(targets))
        return FitSummary(self.df_, self.k_phi_, best[0])

    def predict(self, phi_cb: float, gamma, return_std: bool = True) -> dict:
        if self.df_ is None or self.k_phi_ is None or self.residual_gpr_ is None:
            raise RuntimeError("Fit the model before prediction.")
        gamma = np.asarray(gamma, dtype=float)
        _, eta_phy, details = solve_rheology_curve(
            gamma, phi_cb, self.df_, self.k_phi_, self.physics,
            n_scan=self.calibration["q_scan_prediction"],
        )
        features = np.column_stack([np.full(len(gamma), phi_cb), np.log10(gamma)])
        scaled = self.residual_scaler_.transform(features)
        if return_std:
            delta, std = self.residual_gpr_.predict(scaled, return_std=True)
        else:
            delta, std = self.residual_gpr_.predict(scaled), np.zeros(len(gamma))
        log_eta = np.log10(eta_phy) + delta
        return {
            "eta": 10.0**log_eta,
            "eta_physics": eta_phy,
            "delta_log10": delta,
            "std_log10": std,
            "eta_low95": 10.0 ** (log_eta - 1.96 * std),
            "eta_high95": 10.0 ** (log_eta + 1.96 * std),
            "details": details,
        }
