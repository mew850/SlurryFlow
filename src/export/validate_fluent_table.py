from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np

from .fluent_table import interpolate_table


def validate(model_path, table_path, output_path, n_points=40, seed=42):
    with model_path.open("rb") as stream:
        model = pickle.load(stream)
    table = np.load(table_path)
    phi, log_gamma, log_eta = table["phi"], table["log_gamma"], table["log_eta"]
    rng = np.random.default_rng(seed)
    errors = []
    for _ in range(n_points):
        p = rng.uniform(phi[0], phi[-1])
        lg = rng.uniform(log_gamma[0], log_gamma[-1])
        exact = model.predict(p, [10.0**lg], False)["eta"][0]
        interpolated = interpolate_table(p, 10.0**lg, phi, log_gamma, log_eta)
        errors.append(abs(np.log10(interpolated / exact)))
    result = {
        "n_random_points": n_points,
        "max_abs_log10_error": float(np.max(errors)),
        "mean_abs_log10_error": float(np.mean(errors)),
        "max_multiplicative_error": float(10.0**np.max(errors) - 1.0),
        "viscosity_min_Pa_s": float(10.0**np.min(log_eta)),
        "viscosity_max_Pa_s": float(10.0**np.max(log_eta)),
    }
    output_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=Path("artifacts/model_c/model_c.pkl"))
    parser.add_argument("--table", type=Path, default=Path("artifacts/fluent_table/model_c_table.npz"))
    parser.add_argument("--output", type=Path, default=Path("artifacts/fluent_table/validation.json"))
    args = parser.parse_args()
    result = validate(args.model, args.table, args.output)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
