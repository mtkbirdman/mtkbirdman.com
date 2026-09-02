"""Fit the canonical parabolic drag polar from the archived VSPAERO polar columns."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "aero_reference" / "G103A_drag_polar.csv"
DEFAULT_OUTPUT = ROOT / "aero_reference" / "G103A_drag_polar_fit.json"


def load_points(path: Path):
    with path.open("r", encoding="utf-8", newline="") as file:
        rows = list(csv.DictReader(file))
    alpha = np.asarray([float(row["AoA_deg"]) for row in rows], dtype=float)
    cl = np.asarray([float(row["CLwtot"]) for row in rows], dtype=float)
    cd = np.asarray([float(row["CDwtot"]) for row in rows], dtype=float)
    return alpha, cl, cd


def fit(alpha, cl, cd):
    design = np.column_stack((np.ones_like(cl), cl**2))
    cd0, k = np.linalg.lstsq(design, cd, rcond=None)[0]
    prediction = design @ np.asarray([cd0, k])
    residual = cd - prediction
    ss_res = float(np.sum(residual**2))
    ss_tot = float(np.sum((cd - np.mean(cd))**2))
    return {
        "model": "CD=CD0+K*CL^2",
        "source_original_file": "G103A.polar",
        "source_repository_file": "aero_reference/G103A_drag_polar.csv",
        "source_columns": ["CLwtot", "CDwtot"],
        "fit_method": "unweighted least squares of CDwtot against [1, CLwtot^2]",
        "alpha_range_deg": [float(np.min(alpha)), float(np.max(alpha))],
        "point_count": int(len(alpha)),
        "CD0": float(cd0),
        "K": float(k),
        "r_squared": 1.0 - ss_res / ss_tot,
        "rmse": float(np.sqrt(np.mean(residual**2))),
        "max_abs_residual": float(np.max(np.abs(residual))),
        "points": [
            {
                "AoA_deg": float(a),
                "CLwtot": float(c_l),
                "CDwtot": float(c_d),
                "CD_fit": float(c_fit),
                "residual": float(r),
            }
            for a, c_l, c_d, c_fit, r in zip(alpha, cl, cd, prediction, residual)
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    alpha, cl, cd = load_points(args.input)
    result = fit(alpha, cl, cd)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("CD0", "K", "r_squared", "rmse", "max_abs_residual")}, indent=2))


if __name__ == "__main__":
    main()
