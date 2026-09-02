"""Regenerate Club dynamic-FEM launch-power-ramp sensitivity and early dt convergence."""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from winch_launch_simulator.simulation import load_case, run_simulation
FEM = ROOT / "cases" / "g103a_club_dynamic_winch_fem.json"


def event(result, code):
    return next((item for item in result["events"] if item["code"] == code), None)


def summarize_early(case, result, ramp_duration_s=None):
    h = result["history"]
    t = np.asarray(h["time_s"], dtype=float)
    liftoff = event(result, "LIFTOFF")
    if liftoff is None:
        raise RuntimeError("LIFTOFF not found in early-transient run")
    i = int(np.argmin(np.abs(t - liftoff["time"])))
    row = {
        "dt_s": case["simulation"]["dt_s"],
        "liftoff_time_s": float(liftoff["time"]),
        "liftoff_distance_m": float(h["x_m"][i] - h["x_m"][0]),
        "liftoff_eas_m_s": float(h["eas_m_s"][i]),
        "liftoff_ground_speed_m_s": float(h["ground_speed_x_m_s"][i]),
        "mean_ground_acceleration_m_s2": float(h["ground_speed_x_m_s"][i] / liftoff["time"]),
        "max_tension_to_end_n": float(np.max(np.asarray(h["cable_tension_n"], dtype=float))),
        "max_eas_to_end_m_s": float(np.max(np.asarray(h["eas_m_s"], dtype=float))),
        "warnings": [item["code"] for item in result["events"] if item["severity"] == "warning"],
        "stops": [item["code"] for item in result["events"] if item["severity"] == "stop"],
    }
    if ramp_duration_s is not None:
        row["ramp_duration_s"] = ramp_duration_s
    return row


def early_case(dt_s=None, ramp_s=None):
    case = copy.deepcopy(load_case(FEM))
    case["simulation"]["max_time_s"] = 6.0
    case["simulation"]["post_release_duration_s"] = 0.0
    if dt_s is not None:
        case["simulation"]["dt_s"] = dt_s
    if ramp_s is not None:
        case["winch"]["control"]["launch_power_ramp_duration_s"] = ramp_s
    return case


def main():
    ramp_rows = []
    for ramp in (2.5, 3.0, 3.5, 4.0, 4.5):
        case = early_case(ramp_s=ramp)
        row = summarize_early(case, run_simulation(case), ramp)
        ramp_rows.append(row)
        print(json.dumps({"section": "dynamic_fem_ramp", **row}), flush=True)

    dt_rows = []
    for dt in (0.005, 0.0025, 0.00125):
        case = early_case(dt_s=dt)
        row = summarize_early(case, run_simulation(case))
        dt_rows.append(row)
        print(json.dumps({"section": "dynamic_fem_dt", **row}), flush=True)

    record = {
        "version": "v1.3.0",
        "architecture_change": "elastic Secant removed; Club throttle-ramp sensitivity applies only to Dynamic winch + FEM",
        "purpose": "Club dynamic-FEM launch-power-ramp sensitivity and early-transient dt/dt2/dt4 convergence",
        "canonical_ramp_duration_s": 3.5,
        "dynamic_fem_ramp_sweep": ramp_rows,
        "dynamic_fem_early_transient_dt_convergence": dt_rows,
        "scope_note": "The ideal-Secant pair uses a common prescribed-tension excitation and has no throttle-ramp parameter.",
        "source_version": "v1.3.0",
    }
    out = ROOT / "validation" / "reference" / "launch_power_ramp.json"
    out.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(out)


if __name__ == "__main__":
    main()
