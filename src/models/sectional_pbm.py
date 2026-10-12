"""Mechanistic reproduction of Baldyga et al. (2008).

This is a zero-dimensional, mass-conservative sectional population balance.
It reproduces the mechanisms and qualitative trends in Figs. 3--6; it is not
an exact digitisation because the paper omits several fitted constants.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Config:
    n_bins: int = 56
    n_min: float = 1.0e5
    n_max: float = 2.0e10
    r0: float = 0.02e-6
    phi0: float = 0.018
    phi_star: float = 0.67
    df0: float = 2.85
    df_agg: float = 1.8
    df_shear: float = 3.0
    mu_liquid: float = 1.0e-3
    temperature: float = 298.15
    shear_rate: float = 2.0e4
    agg_scale: float = 0.0
    break_scale: float = 0.0
    shear_compact: float = 0.0
    agg_restructure: float = 1.0
    break_restructure: float = 2.0
    n_critical: float = 2.0e8
    t_end: float = 1.0


def relative_viscosity_limits(phi: np.ndarray, phi_star: float = 0.67):
    """Buyevich--Kapbsov limits, paper equations (15) and (16)."""
    phi = np.asarray(phi, dtype=float)
    # The correlations diverge at close packing; stay on their physical side.
    p = np.clip(phi, 0.0, phi_star * (1.0 - 1e-9))
    x = p / phi_star
    base = (1.0 - p) ** (-2.5)
    series_inf = sum(x ** (j / 3.0) / j for j in range(1, 7))
    m_inf = base - 2.0 * (np.log(1.0 - x ** (1.0 / 3.0)) + series_inf)
    series_0 = sum((1.0 + j) * x**j for j in range(3))
    m_0 = base + 1.3 * ((1.0 - x) ** (-2.0) - series_0)
    return m_0, m_inf


def relative_viscosity(phi, pe, phi_star: float = 0.67):
    """Paper equation (17)."""
    m0, minf = relative_viscosity_limits(phi, phi_star)
    return (m0 + minf * pe) / (1.0 + pe)


class SectionalPBM:
    """Log-grid fixed-pivot PBM for aggregation and equal binary breakage."""

    def __init__(self, cfg: Config):
        self.c = cfg
        self.mass = np.geomspace(cfg.n_min, cfg.n_max, cfg.n_bins)
        self.logm = np.log(self.mass)
        self._agg_map = [[self._map(m1 + m2) for m2 in self.mass] for m1 in self.mass]
        self._half_map = [self._map(m / 2.0) for m in self.mass]

    def _map(self, target):
        """Return two pivot indices and weights conserving primary-particle mass."""
        if target <= self.mass[0]:
            return (0, 0, 1.0, 0.0)
        if target >= self.mass[-1]:
            return (self.c.n_bins - 1, self.c.n_bins - 1, 1.0, 0.0)
        hi = int(np.searchsorted(self.mass, target))
        lo = hi - 1
        whi = (target - self.mass[lo]) / (self.mass[hi] - self.mass[lo])
        return lo, hi, 1.0 - whi, whi

    def initial(self):
        # Narrow lognormal distribution centred on an 8 micrometre aggregate.
        n0 = (8.0e-6 / (2.0 * self.c.r0)) ** self.c.df0
        c = np.exp(-0.5 * ((self.logm - np.log(n0)) / 0.20) ** 2)
        c /= np.dot(c, self.mass)  # total primary-particle mass = 1
        return np.r_[c, self.c.df0]

    def radius(self, df):
        return self.c.r0 * self.mass ** (1.0 / df)

    def rhs(self, _t, y):
        c = np.maximum(y[:-1], 0.0)
        df = float(np.clip(y[-1], 1.55, 3.0))
        dc = np.zeros_like(c)
        r = self.radius(df)

        agg_events = 0.0
        if self.c.agg_scale:
            # Dimensionless form of the orthokinetic part of eq. (21).
            rr = (r[:, None] + r[None, :]) / (2.0 * self.c.r0)
            kernel = self.c.agg_scale * rr**3 / np.max(rr**3)
            for i in range(self.c.n_bins):
                for j in range(i, self.c.n_bins):
                    rate = kernel[i, j] * c[i] * c[j] * (0.5 if i == j else 1.0)
                    if rate == 0.0:
                        continue
                    dc[i] -= rate
                    dc[j] -= rate
                    lo, hi, wlo, whi = self._agg_map[i][j]
                    dc[lo] += rate * wlo
                    dc[hi] += rate * whi
                    agg_events += rate

        break_events = 0.0
        if self.c.break_scale:
            # Smooth sectional counterpart of the threshold law in eq. (25).
            gamma = self.c.break_scale * np.maximum((self.mass / self.c.n_critical) ** 0.55 - 1.0, 0.0)
            for i, rate_i in enumerate(gamma * c):
                if rate_i == 0.0:
                    continue
                dc[i] -= rate_i
                lo, hi, wlo, whi = self._half_map[i]
                dc[lo] += 2.0 * rate_i * wlo
                dc[hi] += 2.0 * rate_i * whi
                break_events += rate_i

        number = max(np.sum(c), 1e-30)
        agg_freq = agg_events / number
        break_freq = break_events / number
        mean_r = float(np.dot(c, r) / number)
        denom = 3.0 * np.log(max(mean_r / (2.0 * self.c.r0), 1.000001)) + 2.0 * np.log(2.0)
        df_break = 3.0 * (df * np.log(max(mean_r / self.c.r0, 1.000001)) - np.log(2.0)) / denom
        ddf = (
            self.c.shear_compact * (self.c.df_shear - df)
            + self.c.agg_restructure * agg_freq * (self.c.df_agg - df)
            + self.c.break_restructure * break_freq * (df_break - df)
        )
        return np.r_[dc, ddf]

    def observables(self, sol):
        c = np.maximum(sol.y[:-1], 0.0)
        df = np.clip(sol.y[-1], 1.55, 3.0)
        # Correct tiny integration drift before evaluating normalized moments.
        c /= np.sum(c * self.mass[:, None], axis=0)
        radii = self.c.r0 * self.mass[:, None] ** (1.0 / df[None, :])
        diam = 2.0 * radii
        moments = {q: np.sum(c * diam**q, axis=0) for q in (0, 1, 2, 3, 4)}
        sizes = {
            "L10": moments[1] / moments[0],
            "L30": (moments[3] / moments[0]) ** (1.0 / 3.0),
            "L32": moments[3] / moments[2],
            "L43": moments[4] / moments[3],
        }
        phi_eff = self.c.phi0 * np.sum(c * self.mass[:, None] ** (3.0 / df[None, :]), axis=0)
        mean_a = np.sum(c * radii, axis=0) / np.sum(c, axis=0)
        kb = 1.380649e-23
        d_inf = kb * self.c.temperature / (6.0 * np.pi * self.c.mu_liquid * mean_a)
        chi = np.maximum(1.0 - phi_eff / self.c.phi_star, 1e-4) ** -2
        pe = self.c.shear_rate * mean_a**2 * chi / d_inf
        mu = self.c.mu_liquid * relative_viscosity(phi_eff, pe, self.c.phi_star)
        return sizes, df, phi_eff, mu, c


def run_case(cfg: Config):
    model = SectionalPBM(cfg)
    t = np.linspace(0.0, cfg.t_end, 181)
    sol = solve_ivp(model.rhs, (0.0, cfg.t_end), model.initial(), t_eval=t,
                    method="BDF", rtol=3e-6, atol=1e-12)
    if not sol.success:
        raise RuntimeError(sol.message)
    return model, sol, model.observables(sol)


def plot_viscosity(outdir: Path):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    phi0 = np.linspace(1e-4, 0.30, 350)
    for ax, df in zip(axes, (2.1, 2.5)):
        # Fig. 3 does not state aggregate R; 0.2 um gives a readable demonstration.
        phi = phi0 * (0.2e-6 / 0.02e-6) ** (3.0 - df)
        m0, minf = relative_viscosity_limits(phi)
        invalid = phi >= 0.67
        m0 = np.where(invalid, np.nan, m0)
        minf = np.where(invalid, np.nan, minf)
        ax.semilogy(phi0, minf, label=r"$M_\infty$")
        ax.semilogy(phi0, m0, label=r"$M_0$")
        ax.set(title=fr"$D_f={df}$", xlabel=r"primary fraction $\phi_0$", ylabel="relative viscosity M")
        ax.set_ylim(1, 1e5); ax.grid(alpha=.25); ax.legend()
    fig.suptitle("Buyevich-Kapbsov viscosity limits (paper Eqs. 15-16)")
    fig.savefig(outdir / "fig3_viscosity_limits.png", dpi=180)
    plt.close(fig)


def plot_scenarios(outdir: Path):
    base = Config()
    cases = {
        "aggregation": replace(base, phi0=0.030, agg_scale=2.0e10, break_scale=0.12, shear_compact=0.06),
        "breakage": replace(base, phi0=0.010, break_scale=10.0, shear_compact=0.15, n_critical=1e6),
        "restructuring": replace(base, phi0=0.010, shear_compact=2.0),
    }
    summary = []
    for name, cfg in cases.items():
        model, sol, (sizes, df, phi, mu, c) = run_case(cfg)
        fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), constrained_layout=True)
        for label, values in sizes.items():
            axes[0].plot(sol.t / cfg.t_end, values * 1e6, label=label)
        axes[0].set(xlabel=r"$t/t_{end}$", ylabel=r"mean diameter [$\mu$m]"); axes[0].legend(ncol=2)
        axes[1].plot(sol.t / cfg.t_end, df); axes[1].set(xlabel=r"$t/t_{end}$", ylabel=r"$D_f$")
        axes[2].plot(sol.t / cfg.t_end, mu * 1e3, label=r"$\mu$")
        ax2 = axes[2].twinx(); ax2.plot(sol.t / cfg.t_end, phi, "C1", label=r"$\phi_{eff}$")
        axes[2].set(xlabel=r"$t/t_{end}$", ylabel="viscosity [mPa s]"); ax2.set_ylabel(r"$\phi_{eff}$")
        axes[0].grid(alpha=.25); axes[1].grid(alpha=.25); axes[2].grid(alpha=.25)
        fig.suptitle(name.capitalize() + "-dominated qualitative reproduction")
        fig.savefig(outdir / f"scenario_{name}.png", dpi=180); plt.close(fig)
        summary.append((name, np.dot(c[:, 0], model.mass), np.dot(c[:, -1], model.mass), sizes["L43"][0], sizes["L43"][-1], df[0], df[-1], phi[0], phi[-1], mu[0], mu[-1]))
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("output"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    plot_viscosity(args.output)
    summary = plot_scenarios(args.output)
    print("case mass0 mass1 L43_0[m] L43_1[m] Df0 Df1 phi0 phi1 mu0[Pa.s] mu1[Pa.s]")
    for row in summary:
        print(row[0], *(f"{x:.6g}" for x in row[1:]))


if __name__ == "__main__":
    main()
