"""Regenerate the canonical Club Vw feedback-path reference case."""
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
CASE = ROOT / "cases" / "g103a_club_dynamic_winch_fem.json"
OUT = ROOT / "validation" / "reference" / "vw_feedback.json"


def event(result, code):
    return next((item for item in result["events"] if item["code"] == code), None)


def main():
    case = copy.deepcopy(load_case(CASE))
    case["atmosphere"]["wind_x_m_s"] = -4.5
    result = run_simulation(case)
    h = result["history"]
    agl = np.asarray(h["agl_m"], dtype=float)
    eas = np.asarray(h["eas_m_s"], dtype=float)
    stall = np.asarray(h["stall_margin_m_s"], dtype=float)
    pilot_phase = np.asarray(h["pilot_phase"], dtype=object)

    request = event(result, "WINCH_POWER_REDUCTION_REQUEST")
    release = event(result, "PILOT_COMMAND_RELEASE")
    reduced_indices = [
        i for i, value in enumerate(h["club_speed_reduction_active"]) if bool(value)
    ]
    reduced_time = None if not reduced_indices else float(h["time_s"][reduced_indices[0]])
    airborne = pilot_phase != "GROUND_ROLL_TAILWHEEL"
    below = agl < 60.96
    above = ~below

    record = {
        "version": "v1.3.0",
        "physics_change": "Iyy=2000 kg m^2 from AMC 22.481-derived estimate; mainwheel+tailwheel two-point prelaunch equilibrium",
        "case_id": case["case_id"],
        "wind_x_m_s": -4.5,
        "request_time_s": None if request is None else float(request["time"]),
        "reduced_throttle_active_time_s": reduced_time,
        "operator_delay_s": case["winch"]["control"]["speed_protection"]["operator_response_delay_s"],
        "normal_release": release is not None,
        "release_time_s": None if release is None else float(release["time"]),
        "max_eas_below_200ft_m_s": float(np.max(eas[below])),
        "max_eas_at_or_above_200ft_m_s": float(np.max(eas[above])),
        "vw_m_s": case["limits"]["max_winch_eas_m_s"],
        "min_airborne_stall_margin_m_s": float(np.min(stall[airborne])),
        "warnings": [item["code"] for item in result["events"] if item["severity"] == "warning"],
        "stops": [item["code"] for item in result["events"] if item["severity"] == "stop"],
        "source_version": "v1.3.0",
    }
    OUT.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
