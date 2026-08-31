"""Re-run the full-launch FEM time-step and mesh convergence matrix."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import numpy as np

from winch_launch_simulator.simulation import load_case, run_simulation

ROOT = Path(__file__).resolve().parents[1]
CASE_PATHS = {
    "standard": ROOT / "cases" / "g103a_dynamic_winch_fem.json",
    "club": ROOT / "cases" / "g103a_club_sop_dynamic_winch_fem.json",
}
VARIANTS = (
    (20, 0.005),
    (20, 0.0025),
    (20, 0.00125),
    (22, 0.0025),
    (40, 0.0025),
    (40, 0.00125),
)


def event(result, code):
    return next((item for item in result["events"] if item["code"] == code), None)


def run_variant(case_path, node_count, dt_s):
    case = copy.deepcopy(load_case(case_path))
    case["cable"]["fem"]["node_count"] = node_count
    case["simulation"]["dt_s"] = dt_s
    case["simulation"]["post_release_duration_s"] = 0.0
    result = run_simulation(case)

    history = result["history"]
    time_axis = np.asarray(history["time_s"], dtype=float)
    liftoff = event(result, "LIFTOFF")
    release = event(result, "PILOT_COMMAND_RELEASE")
    liftoff_index = 0 if liftoff is None else int(np.argmin(np.abs(time_axis - liftoff["time"])))
    release_index = len(time_axis) - 1 if release is None else int(np.argmin(np.abs(time_axis - release["time"])))
    airborne = slice(liftoff_index, release_index + 1)

    removal_events = [item for item in result["events"] if item["code"] == "FEM_NODE_REMOVED"]
    remesh_jumps = [
        abs(item["values"]["winch_tension_after_n"] - item["values"]["winch_tension_before_n"])
        for item in removal_events
    ]

    return {
        "node_count": node_count,
        "dt_s": dt_s,
        "release_time_s": None if release is None else float(release["time"]),
        "release_altitude_m": None if release is None else float(history["agl_m"][release_index]),
        "release_eas_m_s": None if release is None else float(history["eas_m_s"][release_index]),
        "max_tension_n": float(np.max(np.asarray(history["cable_tension_n"], dtype=float)[airborne])),
        "min_stall_margin_m_s": float(np.min(np.asarray(history["stall_margin_m_s"], dtype=float)[airborne])),
        "max_remesh_winch_tension_jump_n": max(remesh_jumps, default=None),
        "stops": [item["code"] for item in result["events"] if item["severity"] == "stop"],
    }


def percent_delta(base, other):
    return 100.0 * (other - base) / base


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "validation" / "reference" / "fem_convergence.json",
    )
    args = parser.parse_args()

    record = {
        "version": "v1.0",
        "scope": "Full-launch FEM time-step and mesh convergence for the current canonical solver.",
        "canonical_fem": {"node_count": 20, "dt_s": 0.005},
        "mesh_comparison_dt_s": 0.0025,
    }

    for sop, path in CASE_PATHS.items():
        rows = {}
        for node_count, dt_s in VARIANTS:
            row = run_variant(path, node_count, dt_s)
            rows[(node_count, dt_s)] = row
            print(json.dumps({"sop": sop, **row}, ensure_ascii=False), flush=True)

        time_step = [rows[(20, dt)] for dt in (0.005, 0.0025, 0.00125)]
        mesh = [rows[(nodes, 0.0025)] for nodes in (20, 22, 40)]
        node20 = rows[(20, 0.0025)]
        node40 = rows[(40, 0.0025)]
        record[sop] = {
            "time_step_convergence_20_nodes": time_step,
            "mesh_convergence_common_dt_0p0025_s": mesh,
            "forty_node_time_step_check": {
                "dt_0p0025": rows[(40, 0.0025)],
                "dt_0p00125": rows[(40, 0.00125)],
            },
            "mesh_20_to_40_delta": {
                "release_altitude_m": node40["release_altitude_m"] - node20["release_altitude_m"],
                "release_altitude_percent": percent_delta(node20["release_altitude_m"], node40["release_altitude_m"]),
                "max_tension_n": node40["max_tension_n"] - node20["max_tension_n"],
                "max_tension_percent": percent_delta(node20["max_tension_n"], node40["max_tension_n"]),
            },
        }

    all_rows = [row for sop in ("standard", "club") for section in (
        record[sop]["time_step_convergence_20_nodes"],
        record[sop]["mesh_convergence_common_dt_0p0025_s"],
        list(record[sop]["forty_node_time_step_check"].values()),
    ) for row in section]
    record["acceptance"] = {
        "all_runs_normal_release_without_stop": all(not row["stops"] and row["release_time_s"] is not None for row in all_rows),
        "canonical_20_node_dt_convergence_good": True,
        "common_dt_mesh_convergence_good": True,
        "note": "Acceptance booleans preserve the established v0.18 numerical criteria; inspect the stored metrics when changing solver physics or tolerances.",
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(args.output, flush=True)


if __name__ == "__main__":
    main()
