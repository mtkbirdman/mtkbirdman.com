"""Regenerate FEM numerical checks for v1.3.0.

The default run is intentionally targeted and reproducible in ordinary development:
- RWTH Dynamic-FEM full-launch dt / dt/2 comparison at 20 nodes.
- Club Dynamic-FEM early-transient dt / dt/2 / dt/4 results are referenced from
  ``validation/reference/launch_power_ramp.json``.

The FEM equations and spatial discretization are unchanged from v1.1.0.  The archived
20/22/40-node full matrix therefore remains plant-discretization evidence, but its old
policy-dependent trajectory values are not relabelled as v1.3.0 results.

Use ``--full-matrix`` to run a new, computationally expensive full-launch spatial/time
matrix for both v1.3.0 policies.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from winch_launch_simulator.simulation import load_case, run_simulation
CASE_PATHS = {
    "rwth_reference": ROOT / "cases" / "g103a_rwth_reference_dynamic_winch_fem.json",
    "club": ROOT / "cases" / "g103a_club_dynamic_winch_fem.json",
}
FULL_VARIANTS = (
    (20, 0.005),
    (20, 0.0025),
    (20, 0.00125),
    (22, 0.0025),
    (40, 0.0025),
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


def run_variants_isolated(jobs):
    """Run each long FEM case in a fresh Python subprocess.

    Recycling the interpreter after every full-launch case keeps memory and runtime
    predictable without changing solver physics.
    """
    script = Path(__file__).resolve()

    def invoke(job):
        case_path, node_count, dt_s = job
        completed = subprocess.run(
            [
                sys.executable,
                str(script),
                "--worker-case",
                str(case_path),
                "--worker-nodes",
                str(node_count),
                "--worker-dt",
                str(dt_s),
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=ROOT,
        )
        return json.loads(completed.stdout.strip().splitlines()[-1])

    with ThreadPoolExecutor(max_workers=min(2, len(jobs))) as executor:
        return list(executor.map(invoke, jobs))

def targeted_record():
    rwth_jobs = [(str(CASE_PATHS["rwth_reference"]), 20, dt_s) for dt_s in (0.005, 0.0025)]
    rwth_rows = run_variants_isolated(rwth_jobs)
    for row in rwth_rows:
        print(json.dumps({"section": "rwth_full_launch_dt", **row}), flush=True)

    ramp_path = ROOT / "validation" / "reference" / "launch_power_ramp.json"
    ramp = json.loads(ramp_path.read_text(encoding="utf-8"))
    club_dt_rows = ramp["dynamic_fem_early_transient_dt_convergence"]

    return {
        "version": "v1.3.0",
        "scope": (
            "Targeted numerical checks after the plant/policy architecture change. "
            "FEM plant equations are unchanged from v1.1.0; the full v1.1.0 mesh matrix "
            "is archived separately."
        ),
        "canonical_fem": {"node_count": 20, "dt_s": 0.005},
        "rwth_reference_full_launch_time_step": rwth_rows,
        "club_early_transient_time_step": {
            "reference_file": "validation/reference/launch_power_ramp.json",
            "window_s": 6.0,
            "dt_values_s": [row["dt_s"] for row in club_dt_rows],
            "rows": club_dt_rows,
            "note": (
                "Covers the complete 3.5 s launch-power ramp and liftoff; liftoff metrics "
                "are converged to substantially below operational/model uncertainty."
            ),
        },
        "archived_full_mesh_reference": "validation/reference/fem_convergence_v1.1.0.json",
        "interpretation": (
            "The v1.1.0 full 20/22/40-node study remains evidence for the unchanged FEM "
            "plant discretization. It is not relabelled as a v1.3.0 current-policy trajectory result."
        ),
        "acceptance": {
            "rwth_full_launch_normal_release_without_stop": all(
                row["release_time_s"] is not None and not row["stops"] for row in rwth_rows
            ),
            "club_early_transient_without_stop": all(not row["stops"] for row in club_dt_rows),
        },
    }


def full_matrix_record():
    record = {
        "version": "v1.3.0",
        "scope": "Computationally expensive full-launch time-step and mesh matrix for both current FEM policies.",
        "canonical_fem": {"node_count": 20, "dt_s": 0.005},
        "mesh_comparison_dt_s": 0.0025,
    }
    for policy, path in CASE_PATHS.items():
        jobs = [(str(path), node_count, dt_s) for node_count, dt_s in FULL_VARIANTS]
        computed = run_variants_isolated(jobs)
        rows = {}
        for (node_count, dt_s), row in zip(FULL_VARIANTS, computed):
            rows[(node_count, dt_s)] = row
            print(json.dumps({"policy": policy, **row}), flush=True)
        node20 = rows[(20, 0.0025)]
        node40 = rows[(40, 0.0025)]
        record[policy] = {
            "time_step_convergence_20_nodes": [rows[(20, dt)] for dt in (0.005, 0.0025, 0.00125)],
            "mesh_convergence_common_dt_0p0025_s": [rows[(nodes, 0.0025)] for nodes in (20, 22, 40)],
            "mesh_20_to_40_delta": {
                "release_altitude_m": node40["release_altitude_m"] - node20["release_altitude_m"],
                "release_altitude_percent": percent_delta(node20["release_altitude_m"], node40["release_altitude_m"]),
                "max_tension_n": node40["max_tension_n"] - node20["max_tension_n"],
                "max_tension_percent": percent_delta(node20["max_tension_n"], node40["max_tension_n"]),
            },
        }
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full-matrix", action="store_true")
    parser.add_argument("--worker-case", type=Path, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--worker-nodes", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--worker-dt", type=float, default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "validation" / "reference" / "fem_convergence.json",
    )
    args = parser.parse_args()

    if args.worker_case is not None:
        if args.worker_nodes is None or args.worker_dt is None:
            parser.error("worker mode requires --worker-nodes and --worker-dt")
        print(json.dumps(run_variant(args.worker_case, args.worker_nodes, args.worker_dt)))
        return

    record = full_matrix_record() if args.full_matrix else targeted_record()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(args.output, flush=True)


if __name__ == "__main__":
    main()
