from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.data.rheology import load_config, load_rheology_curves, validate_curves
from src.models.hybrid_model_c import HybridModelC, build_base_physics
from src.utils.metrics import curve_metrics


def train(config_path: Path, project_root: Path, artifact_dir: Path):
    config = load_config(config_path)
    curves = load_rheology_curves(config, project_root)
    validate_curves(curves)
    model = HybridModelC(build_base_physics(config), config["calibration"])
    summary = model.fit(curves)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    with (artifact_dir / "model_c.pkl").open("wb") as stream:
        pickle.dump(model, stream)

    rows = []
    figure, axis = plt.subplots(figsize=(10, 7), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.05, 0.95, len(curves)))
    for color, curve in zip(colors, curves):
        prediction = model.predict(curve["phi_cb"], curve["gamma"])
        rows.append({
            "sample": curve["sample"], "phi_cb": curve["phi_cb"], "split": curve["split"],
            **curve_metrics(curve["eta"], prediction["eta"]),
        })
        marker = "o" if curve["split"] == "train" else "s"
        axis.loglog(curve["gamma"], curve["eta"], marker, ms=4, ls="none", color=color)
        dense = np.geomspace(curve["gamma"].min(), curve["gamma"].max(), 150)
        axis.loglog(dense, model.predict(curve["phi_cb"], dense, False)["eta"], color=color,
                    label="{} phi={:.4f} ({})".format(curve["sample"], curve["phi_cb"], curve["split"]))
    axis.set(xlabel="Shear rate [1/s]", ylabel="Viscosity [Pa s]")
    axis.grid(True, which="both", alpha=0.25); axis.legend(fontsize=8, ncol=2)
    figure.savefig(artifact_dir / "model_c_validation.png", dpi=180); plt.close(figure)

    pd.DataFrame(rows).to_csv(artifact_dir / "validation_metrics.csv", index=False)
    metadata = {
        "df": summary.df, "k_phi": summary.k_phi, "objective_mse": summary.objective_mse,
        "residual_kernel": str(model.residual_gpr_.kernel_),
        "train_samples": [c["sample"] for c in curves if c["split"] == "train"],
        "test_samples": [c["sample"] for c in curves if c["split"] == "test"],
    }
    (artifact_dir / "model_c_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return model, config, curves


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("configs/model_c.json"))
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--artifact-dir", type=Path, default=Path("artifacts/model_c"))
    args = parser.parse_args()
    train(args.config, args.project_root, args.artifact_dir)


if __name__ == "__main__":
    main()
