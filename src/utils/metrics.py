import numpy as np


def curve_metrics(y_true, y_pred):
    y_true, y_pred = np.asarray(y_true), np.asarray(y_pred)
    error = np.log10(y_pred) - np.log10(y_true)
    centered = np.log10(y_true) - np.mean(np.log10(y_true))
    return {
        "log_mae": float(np.mean(np.abs(error))),
        "log_rmse": float(np.sqrt(np.mean(error**2))),
        "mape_percent": float(np.mean(np.abs(y_pred / y_true - 1.0)) * 100.0),
        "r2_log": float(1.0 - np.sum(error**2) / max(np.sum(centered**2), 1e-30)),
    }
