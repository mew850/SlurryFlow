from __future__ import annotations

import json
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd


def load_config(path: Union[str, Path]) -> dict:
    with Path(path).open("r", encoding="utf-8") as stream:
        return json.load(stream)


def load_rheology_curves(config: dict, project_root: Union[str, Path] = ".") -> list:
    root = Path(project_root)
    data_dir = root / config["data_dir"]
    scale = 1e-3 if config.get("viscosity_input_unit") == "mPa.s" else 1.0
    curves = []
    for item in config["samples"]:
        frame = pd.read_csv(data_dir / item["file"])
        gamma = pd.to_numeric(frame["Shear Rate"], errors="coerce").to_numpy(float)
        eta = pd.to_numeric(frame["Viscosity"], errors="coerce").to_numpy(float) * scale
        valid = np.isfinite(gamma) & np.isfinite(eta) & (gamma > 0) & (eta > 0)
        order = np.argsort(gamma[valid])
        curves.append({**item, "gamma": gamma[valid][order], "eta": eta[valid][order]})
    return curves


def validate_curves(curves: list) -> None:
    if not curves:
        raise ValueError("No rheology curves were loaded.")
    phis = [curve["phi_cb"] for curve in curves]
    if len(set(phis)) != len(phis):
        raise ValueError("Each curve must have a unique phi_cb.")
    if not any(curve["split"] == "train" for curve in curves):
        raise ValueError("At least one training curve is required.")
