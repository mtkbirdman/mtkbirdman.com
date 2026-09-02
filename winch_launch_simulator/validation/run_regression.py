"""Re-run the four canonical cases and write the current regression record."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from winch_launch_simulator.simulation import load_case, run_simulation
CASE_PATHS = [
    ROOT / "cases" / "g103a_rwth_reference_prescribed_tension.json",
    ROOT / "cases" / "g103a_club_prescribed_tension.json",
    ROOT / "cases" / "g103a_rwth_reference_dynamic_winch_fem.json",
    ROOT / "cases" / "g103a_club_dynamic_winch_fem.json",
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


def run_one(path):
    case = copy.deepcopy(load_case(path))
    case["simulation"]["post_release_duration_s"] = 0.0
    started = time.perf_counter()
    result = run_simulation(case)
    return summarize(case, result, time.perf_counter() - started)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "validation" / "reference" / "canonical_regression.json")
    parser.add_argument("--single-case", type=Path, default=None, help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.single_case is not None:
        print(json.dumps(run_one(args.single_case), ensure_ascii=False))
        return

    # Run every canonical launch in a fresh Python process.  The two FEM cases are
    # independent, so running subprocesses concurrently avoids retaining large
    # histories between cases and keeps the regression command practical.
    def invoke(path):
        completed = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--single-case", str(path)],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=True,
        )
        return json.loads(completed.stdout.strip().splitlines()[-1])

    with ThreadPoolExecutor(max_workers=len(CASE_PATHS)) as executor:
        rows = list(executor.map(invoke, CASE_PATHS))
    for row in rows:
        print(json.dumps(row, ensure_ascii=False), flush=True)

    record = {
        "version": "v1.3.0",
        "architecture": "2 plant fidelities x 2 pilot/operation policies",
        "physics_change": "Iyy=2000 kg m^2 AMC-based estimate; mainwheel+tailwheel two-point prelaunch equilibrium",
        "cases": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(args.output, flush=True)


if __name__ == "__main__":
    main()
