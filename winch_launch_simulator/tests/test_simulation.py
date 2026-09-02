import copy
import math

import numpy as np
import pytest

from winch_launch_simulator.initialization import initialize_launch_state
from winch_launch_simulator.models.controllers import club_pilot_control, rwth_reference_pilot_control
from winch_launch_simulator.models.winch import full_throttle_drive_torque
from winch_launch_simulator.simulation import TAUT, detect_step_event, evaluate_dynamics, load_case, run_simulation


def event_state(**updates):
    state = {
        "cable_tension_n": 1000.0,
        "eas_m_s": 25.0,
        "alpha_deg": 2.0,
        "stall_margin_m_s": 5.0,
        "cable_ground_angle_deg": 20.0,
        "cable_los_ground_angle_deg": 20.0,
        "pitch_rate_deg_s": 0.0,
        "gamma_ground_deg": 0.0,
        "pilot_phase": "RWTH_AIRSPEED_CONTROL",
        "pilot_release_preparation_progress": 0.0,
        "pilot_release_pitch_error_deg": 0.0,
        "pilot_power_reduction_requested": False,
        "winch_terminal_fade_progress": 0.0,
        "aircraft_ground_signed_depth_m": -0.05,
        "agl_m": 100.0,
    }
    state.update(updates)
    return state


def airborne_pilot_state(**updates):
    state = {
        "liftoff_time_s": 0.0,
        "rotation_start_time_s": None,
        "rotation_start_pitch_rad": None,
        "rwth_safety_altitude_time_s": 1.0,
        "eas_rate_m_s2": 0.0,
        "power_reduction_request_time_s": None,
        "terminal_fade_progress": 0.0,
    }
    state.update(updates)
    return state


def test_load_case_requires_positive_parabolic_drag_coefficients(tmp_path):
    import json
    source = load_case("cases/g103a_rwth_reference_prescribed_tension.json")
    for field in ("CD0", "K"):
        case = copy.deepcopy(source)
        case["aerodynamics"][field] = 0.0
        path = tmp_path / f"bad_{field}.json"
        path.write_text(json.dumps(case), encoding="utf-8")
        with pytest.raises(ValueError, match=field):
            load_case(path)


def test_ground_start_initializer_returns_static_ideal_secant_state():
    case = load_case("cases/g103a_rwth_reference_prescribed_tension.json")
    state, metadata = initialize_launch_state(case)
    assert state.shape == (9,)
    assert np.isfinite(state).all()
    assert np.allclose(state[2:4], 0.0)
    assert metadata["ideal_secant"] is True
    assert metadata["initial_reel_in_m"] == 0.0


def test_rwth_reference_is_trim_only_before_safety_altitude():
    case = load_case("cases/g103a_rwth_reference_prescribed_tension.json")
    air = {"density_kg_m3": 1.225, "airspeed_m_s": 20.0}
    result = rwth_reference_pilot_control(case, 2.0, air, 0.2, {"rwth_safety_altitude_time_s": None, "terminal_fade_progress": 0.0})
    assert result["phase"] == "RWTH_TRIM"
    assert math.degrees(result["elevator_rad"]) == pytest.approx(-3.0)


def test_rwth_reference_airspeed_feedback_has_correct_direction_after_fade():
    case = load_case("cases/g103a_rwth_reference_prescribed_tension.json")
    state = {"rwth_safety_altitude_time_s": 0.0, "terminal_fade_progress": 0.0}
    slow = rwth_reference_pilot_control(case, 10.0, {"density_kg_m3": 1.225, "airspeed_m_s": 25.0}, 0.0, state)
    fast = rwth_reference_pilot_control(case, 10.0, {"density_kg_m3": 1.225, "airspeed_m_s": 32.0}, 0.0, state)
    # Positive elevator is nose-down in the simulator; slow flight commands more nose-down.
    assert slow["elevator_rad"] > fast["elevator_rad"]


def test_club_ground_roll_adds_up_elevator_when_tailwheel_unloads():
    case = load_case("cases/g103a_club_prescribed_tension.json")
    air = {"density_kg_m3": 1.225, "airspeed_m_s": 10.0, "eas_m_s": 10.0, "gamma_ground_rad": 0.0}
    state = {"liftoff_time_s": None, "rotation_start_time_s": None, "terminal_fade_progress": 0.0}
    loaded = club_pilot_control(case, 0.0, air, 0.0, 0.0, 0.0, 2500.0, 600.0, state)
    unloaded = club_pilot_control(case, 0.0, air, 0.0, 0.0, 0.0, 2500.0, 0.0, state)
    assert unloaded["elevator_rad"] < loaded["elevator_rad"]


def test_g72_g73_load_limit_has_priority_over_weak_link_when_both_cross():
    case = copy.deepcopy(load_case("cases/g103a_rwth_reference_prescribed_tension.json"))
    case["limits"]["weak_link_break_load_n"] = 12000.0
    event = detect_step_event(event_state(cable_tension_n=11000.0), event_state(cable_tension_n=12100.0), case, TAUT, airborne_pilot_state())
    assert event["code"] == "G72_G73_LOAD_EXCEEDED"


def test_liftoff_crossing_is_detected_from_last_wheel_clearance():
    case = load_case("cases/g103a_rwth_reference_prescribed_tension.json")
    event = detect_step_event(
        event_state(aircraft_ground_signed_depth_m=0.002),
        event_state(aircraft_ground_signed_depth_m=-0.002),
        case,
        TAUT,
        {"liftoff_time_s": None},
    )
    assert event["code"] == "LIFTOFF"
    assert event["fraction"] == 0.5


def test_stall_and_aero_validity_events_are_disabled_before_liftoff():
    case = load_case("cases/g103a_rwth_reference_prescribed_tension.json")
    event = detect_step_event(
        event_state(alpha_deg=13.0, stall_margin_m_s=1.0),
        event_state(alpha_deg=15.0, stall_margin_m_s=-1.0),
        case,
        TAUT,
        {"liftoff_time_s": None},
    )
    assert event is None


def test_rwth_release_requires_terminal_fade_low_tension_and_los_angle():
    case = load_case("cases/g103a_rwth_reference_prescribed_tension.json")
    release = case["pilot"]["release"]
    previous = event_state(cable_tension_n=600.0, cable_los_ground_angle_deg=77.0, winch_terminal_fade_progress=0.99)
    current = event_state(cable_tension_n=400.0, cable_los_ground_angle_deg=78.0, winch_terminal_fade_progress=1.0)
    event = detect_step_event(previous, current, case, TAUT, airborne_pilot_state())
    assert current["cable_tension_n"] <= release["max_tension_n"]
    assert event["code"] == "PILOT_COMMAND_RELEASE"


def test_club_vw_becomes_hard_constraint_at_200ft_crossing():
    case = load_case("cases/g103a_club_prescribed_tension.json")
    vw = case["limits"]["max_winch_eas_m_s"]
    previous = event_state(agl_m=60.90, eas_m_s=vw + 0.1, pilot_phase="POST_LIFTOFF_HOLD_DOWN")
    current = event_state(agl_m=61.02, eas_m_s=vw + 0.1, pilot_phase="PROGRESSIVE_ROTATION")
    event = detect_step_event(previous, current, case, TAUT, airborne_pilot_state())
    assert event["code"] == "MAX_WINCH_SPEED_EXCEEDED"
    assert event["severity"] == "stop"


def test_dynamic_winch_control_modes_are_policy_specific():
    rwth = load_case("cases/g103a_rwth_reference_dynamic_winch_fem.json")
    club = load_case("cases/g103a_club_dynamic_winch_fem.json")
    assert rwth["winch"]["control"]["mode"] == "tension_control"
    assert club["winch"]["control"]["mode"] == "club_throttle"


def test_club_launch_power_ramp_starts_at_pretension_equilibrium():
    case = load_case("cases/g103a_club_dynamic_winch_fem.json")
    state, initialization = initialize_launch_state(case)
    runtime = {"initial_reel_in_m": initialization["initial_reel_in_m"], "removed_node_count": 0, "initial_node_count": case["cable"]["fem"]["node_count"]}
    control = case["winch"]["control"]
    pretension = case["launch_setup"]["pretension_n"]
    winch = case["winch"]
    initial_throttle = (pretension * winch["drum_radius_m"] + winch["loss_torque_nm"]) / full_throttle_drive_torque(winch, 0.0)
    pilot_state = airborne_pilot_state(liftoff_time_s=None, rwth_safety_altitude_time_s=None)
    _, diagnostics = evaluate_dynamics(0.0, state, case, TAUT, pilot_state, runtime)
    assert diagnostics["winch_throttle_command"] == pytest.approx(initial_throttle)
    assert diagnostics["winch_launch_power_ramp_progress"] == 0.0
    assert math.isnan(diagnostics["target_tension_n"])


def test_prescribed_pair_full_launches_release_without_stop():
    for path in ["cases/g103a_rwth_reference_prescribed_tension.json", "cases/g103a_club_prescribed_tension.json"]:
        result = run_simulation(load_case(path))
        assert any(event["code"] == "PILOT_COMMAND_RELEASE" for event in result["events"])
        assert not [event for event in result["events"] if event["severity"] == "stop"]


def test_club_prescribed_has_positive_canonical_stall_margin():
    result = run_simulation(load_case("cases/g103a_club_prescribed_tension.json"))
    release_time = next(event["time"] for event in result["events"] if event["code"] == "PILOT_COMMAND_RELEASE")
    airborne = [
        margin
        for t, margin, agl in zip(
            result["history"]["time_s"],
            result["history"]["stall_margin_m_s"],
            result["history"]["agl_m"],
        )
        if agl > 1.0 and t <= release_time
    ]
    assert min(airborne) >= 1.5
