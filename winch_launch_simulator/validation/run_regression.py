"""Re-run the five canonical cases and write the current regression record."""
from __future__ import annotations

import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np

from winch_launch_simulator.simulation import load_case, run_simulation

ROOT = Path(__file__).resolve().parents[1]
CASE_PATHS = [
    ROOT / "cases" / "g103a_prescribed_tension.json",
    ROOT / "cases" / "g103a_dynamic_winch.json",
    ROOT / "cases" / "g103a_club_sop_dynamic_winch.json",
    ROOT / "cases" / "g103a_dynamic_winch_fem.json",
    ROOT / "cases" / "g103a_club_sop_dynamic_winch_fem.json",
]


def event(result, code):
    return next((item for item in result["events"] if item["code"] == code), None)


def summarize(case, result, runtime_s):
    history = result["history"]
    time_axis = np.asarray(history["time_s"], dtype=float)
    liftoff = event(result, "LIFTOFF")
    release = event(result, "PILOT_COMMAND_RELEASE")
    weak_link = event(result, "WEAK_LINK_BREAK")
    terminal = release or weak_link

    liftoff_index = 0 if liftoff is None else int(np.argmin(np.abs(time_axis - liftoff["time"])))
    terminal_index = len(time_axis) - 1 if terminal is None else int(np.argmin(np.abs(time_axis - terminal["time"])))
    airborne = slice(liftoff_index, terminal_index + 1)

    return {
        "case_id": case["case_id"],
        "dt_s": case["simulation"]["dt_s"],
        "runtime_s": runtime_s,
        "normal_release": release is not None,
        "release_time_s": None if release is None else float(release["time"]),
        "release_altitude_m": None if release is None else float(history["agl_m"][terminal_index]),
        "release_eas_m_s": None if release is None else float(history["eas_m_s"][terminal_index]),
        "max_eas_m_s": float(np.max(np.asarray(history["eas_m_s"], dtype=float)[airborne])),
        "max_tension_n": float(np.max(np.asarray(history["cable_tension_n"], dtype=float)[airborne])),
        "min_stall_margin_m_s": float(np.min(np.asarray(history["stall_margin_m_s"], dtype=float)[airborne])),
        "max_abs_pitch_rate_deg_s": float(np.max(np.abs(np.asarray(history["pitch_rate_deg_s"], dtype=float)[airborne]))),
        "max_pitch_deg": float(np.max(np.asarray(history["theta_deg"], dtype=float)[airborne])),
        "warnings": [item["code"] for item in result["events"] if item["severity"] == "warning"],
        "stops": [item["code"] for item in result["events"] if item["severity"] == "stop"],
        "validation": list(result["validation"]),
        "event_codes": [item["code"] for item in result["events"]],
        "initialization": result["metadata"]["initialization"],
        "final_fem_node_count": result["metadata"]["final_fem_node_count"],
        "fem_nodes_removed": result["metadata"]["fem_nodes_removed"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "validation" / "reference" / "canonical_regression.json",
    )
    args = parser.parse_args()

    rows = []
    for path in CASE_PATHS:
        case = copy.deepcopy(load_case(path))
        case["simulation"]["post_release_duration_s"] = 0.0
        started = time.perf_counter()
        result = run_simulation(case)
        row = summarize(case, result, time.perf_counter() - started)
        rows.append(row)
        print(json.dumps(row, ensure_ascii=False), flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps({"version": "v1.0", "cases": rows}, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(args.output, flush=True)


if __name__ == "__main__":
    main()
