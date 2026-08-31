import copy
import math

import numpy as np
import pytest

from winch_launch_simulator.initialization import initialize_launch_state
from winch_launch_simulator.models.winch import full_throttle_drive_torque
from winch_launch_simulator.models.controllers import club_pilot_control, standard_pilot_control
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
        "pilot_phase": "INITIAL_CLIMB",
        "pilot_release_preparation_progress": 0.0,
        "pilot_release_pitch_error_deg": 0.0,
        "pilot_power_reduction_requested": False,
        "winch_terminal_fade_progress": 0.0,
        "aircraft_ground_signed_depth_m": -0.05,
        "agl_m": 100.0,
    }
    state.update(updates)
    return state


def airborne_pilot_state():
    return {
        "liftoff_time_s": 0.0,
        "rotation_start_time_s": None,
        "rotation_start_pitch_rad": None,
        "eas_rate_m_s2": 0.0,
        "power_reduction_request_time_s": None,
        "terminal_fade_progress": 0.0,
    }


def test_ground_start_initializer_returns_static_aircraft_state():
    case = load_case("cases/g103a_prescribed_tension.json")
    state, metadata = initialize_launch_state(case)
    assert state.shape == (9,)
    assert np.isfinite(state).all()
    assert np.allclose(state[2:4], 0.0)
    assert state[5] == 0.0
    assert metadata["ground_start"] is True
    assert metadata["target_pretension_n"] == 2500.0
    assert metadata["initial_cg_height_m"] > 0.0


def test_ground_start_produces_a_liftoff_event_before_rotation():
    case = copy.deepcopy(load_case("cases/g103a_dynamic_winch.json"))
    case["simulation"]["max_time_s"] = 12.0
    result = run_simulation(case)
    liftoff = [event for event in result["events"] if event["code"] == "LIFTOFF"]
    assert len(liftoff) == 1
    assert result["metadata"]["liftoff_time_s"] == liftoff[0]["time"]
    assert result["history"]["pilot_phase"][0] == "GROUND_ROLL"
    assert "INITIAL_CLIMB" in result["history"]["pilot_phase"]


def test_g72_g73_load_limit_has_priority_over_weak_link_when_both_cross():
    case = copy.deepcopy(load_case("cases/g103a_prescribed_tension.json"))
    case["limits"]["weak_link_break_load_n"] = 12000.0
    event = detect_step_event(
        event_state(cable_tension_n=11000.0),
        event_state(cable_tension_n=12100.0),
        case,
        TAUT,
        airborne_pilot_state(),
    )
    assert event["code"] == "G72_G73_LOAD_EXCEEDED"


def test_liftoff_crossing_is_detected_from_the_last_wheel_clearance():
    case = load_case("cases/g103a_prescribed_tension.json")
    event = detect_step_event(
        event_state(aircraft_ground_signed_depth_m=0.002, pilot_phase="GROUND_ROLL"),
        event_state(aircraft_ground_signed_depth_m=-0.002, pilot_phase="GROUND_ROLL"),
        case,
        TAUT,
        {"liftoff_time_s": None},
    )
    assert event["code"] == "LIFTOFF"
    assert event["fraction"] == 0.5


def test_stall_and_aero_validity_events_are_disabled_before_liftoff():
    case = load_case("cases/g103a_prescribed_tension.json")
    event = detect_step_event(
        event_state(alpha_deg=13.0, stall_margin_m_s=1.0, pilot_phase="GROUND_ROLL"),
        event_state(alpha_deg=15.0, stall_margin_m_s=-1.0, pilot_phase="GROUND_ROLL"),
        case,
        TAUT,
        {"liftoff_time_s": None},
    )
    assert event is None


def test_standard_takeoff_attitude_feedback_has_the_correct_direction():
    case = load_case("cases/g103a_prescribed_tension.json")
    air = {"density_kg_m3": 1.225, "airspeed_m_s": 25.0}
    target = math.radians(case["pilot"]["initial_climb"]["pitch_target_deg"])
    state = airborne_pilot_state()
    nose_high = standard_pilot_control(case, 0.0, air, target + math.radians(2), 0.0, 0.0, 0.0, state)
    nose_low = standard_pilot_control(case, 0.0, air, target - math.radians(2), 0.0, 0.0, 0.0, state)
    assert nose_high["elevator_rad"] > nose_low["elevator_rad"]


def test_club_ground_roll_adds_up_elevator_when_tailwheel_unloads():
    case = load_case("cases/g103a_club_sop_dynamic_winch.json")
    air = {"density_kg_m3": 1.225, "airspeed_m_s": 10.0, "eas_m_s": 10.0, "gamma_ground_rad": 0.0}
    state = {"liftoff_time_s": None, "rotation_start_time_s": None, "terminal_fade_progress": 0.0}
    loaded = club_pilot_control(case, 0.0, air, 0.0, 0.0, 0.0, 2500.0, 600.0, state)
    unloaded = club_pilot_control(case, 0.0, air, 0.0, 0.0, 0.0, 2500.0, 0.0, state)
    assert unloaded["elevator_rad"] < loaded["elevator_rad"]


def test_club_rotation_starts_at_200ft():
    case = copy.deepcopy(load_case("cases/g103a_club_sop_dynamic_winch.json"))
    case["limits"]["max_winch_eas_m_s"] = 100.0
    result = run_simulation(case)
    history = result["history"]
    index = next(i for i, phase in enumerate(history["pilot_phase"]) if phase == "PROGRESSIVE_ROTATION")
    assert history["agl_m"][index] >= 60.96 - 0.2


def test_club_vw_becomes_hard_constraint_at_200ft_crossing():
    case = load_case("cases/g103a_club_sop_dynamic_winch.json")
    vw = case["limits"]["max_winch_eas_m_s"]
    previous = event_state(agl_m=60.90, eas_m_s=vw + 0.1, pilot_phase="POST_LIFTOFF_HOLD_DOWN")
    current = event_state(agl_m=61.02, eas_m_s=vw + 0.1, pilot_phase="PROGRESSIVE_ROTATION")
    event = detect_step_event(previous, current, case, TAUT, airborne_pilot_state())
    assert event["code"] == "MAX_WINCH_SPEED_EXCEEDED"
    assert event["severity"] == "stop"
    assert 0.0 < event["fraction"] < 1.0


def test_standard_release_requires_low_tension_angle_pitch_rate_and_pitch_attitude():
    case = load_case("cases/g103a_prescribed_tension.json")
    release = case["pilot"]["release"]
    previous = event_state(
        cable_tension_n=release["max_tension_n"] + 50,
        cable_ground_angle_deg=release["min_ground_angle_deg"] - 0.2,
        pitch_rate_deg_s=release["max_abs_pitch_rate_deg_s"] + 0.2,
        pilot_phase="RELEASE_PREPARATION",
        pilot_release_preparation_progress=0.99,
        pilot_release_pitch_error_deg=release["max_pitch_error_deg"] + 0.2,
    )
    current = event_state(
        cable_tension_n=release["max_tension_n"] - 10,
        cable_ground_angle_deg=release["min_ground_angle_deg"] + 0.2,
        pitch_rate_deg_s=release["max_abs_pitch_rate_deg_s"] - 0.2,
        pilot_phase="RELEASE_PREPARATION",
        pilot_release_preparation_progress=1.0,
        pilot_release_pitch_error_deg=release["max_pitch_error_deg"] - 0.2,
    )
    event = detect_step_event(previous, current, case, TAUT, airborne_pilot_state())
    assert event["code"] == "PILOT_COMMAND_RELEASE"


def test_dynamic_winch_control_modes_remain_distinct():
    standard = load_case("cases/g103a_dynamic_winch.json")
    club = load_case("cases/g103a_club_sop_dynamic_winch.json")
    assert standard["winch"]["control"]["mode"] == "tension_control"
    assert club["winch"]["control"]["mode"] == "club_throttle"


def test_club_launch_power_ramp_starts_at_pretension_drive_and_ignores_liftoff_state():
    case = load_case("cases/g103a_club_sop_dynamic_winch.json")
    state, initialization = initialize_launch_state(case)
    runtime = {"initial_reel_in_m": initialization["initial_reel_in_m"]}
    control = case["winch"]["control"]
    pretension = case["launch_setup"]["pretension_n"]
    winch = case["winch"]
    initial_throttle = (
        pretension * winch["drum_radius_m"] + winch["loss_torque_nm"]
    ) / full_throttle_drive_torque(winch, 0.0)

    ground_state = {
        "liftoff_time_s": None,
        "rotation_start_time_s": None,
        "power_reduction_request_time_s": None,
        "terminal_fade_progress": 0.0,
    }
    _, start = evaluate_dynamics(0.0, state, case, TAUT, ground_state, runtime)
    assert start["winch_throttle_command"] == pytest.approx(initial_throttle)
    assert start["winch_launch_power_ramp_progress"] == 0.0
    assert math.isnan(start["target_tension_n"])

    half_time = 0.5 * control["launch_power_ramp_duration_s"]
    _, mid_ground = evaluate_dynamics(half_time, state, case, TAUT, ground_state, runtime)
    expected_mid = 0.5 * (initial_throttle + control["launch_throttle"])
    assert mid_ground["winch_throttle_command"] == pytest.approx(expected_mid)
    assert mid_ground["winch_launch_power_ramp_progress"] == pytest.approx(0.5)

    airborne_state = dict(ground_state, liftoff_time_s=1.0)
    _, mid_airborne = evaluate_dynamics(half_time, state, case, TAUT, airborne_state, runtime)
    assert mid_airborne["winch_throttle_command"] == pytest.approx(
        mid_ground["winch_throttle_command"]
    )

    _, full = evaluate_dynamics(
        control["launch_power_ramp_duration_s"], state, case, TAUT, airborne_state, runtime
    )
    assert full["winch_throttle_command"] == control["launch_throttle"]
    assert full["winch_launch_power_ramp_progress"] == 1.0


def test_club_headwind_can_request_and_activate_winch_power_reduction():
    case = copy.deepcopy(load_case("cases/g103a_club_sop_dynamic_winch.json"))
    case["atmosphere"]["wind_x_m_s"] = -4.5
    result = run_simulation(case)
    request = [event for event in result["events"] if event["code"] == "WINCH_POWER_REDUCTION_REQUEST"]
    assert request
    active = [
        t for t, value in zip(result["history"]["time_s"], result["history"]["club_speed_reduction_active"])
        if value
    ]
    assert active
    delay = case["winch"]["control"]["speed_protection"]["operator_response_delay_s"]
    assert active[0] >= request[0]["time"] + delay - case["simulation"]["dt_s"]
    assert max(e for e, a in zip(result["history"]["eas_m_s"], result["history"]["agl_m"]) if a >= 60.96) < case["limits"]["max_winch_eas_m_s"]
    assert any(event["code"] == "PILOT_COMMAND_RELEASE" for event in result["events"])


def test_dynamic_drive_torque_never_exceeds_shared_full_throttle_curve():
    for path in ["cases/g103a_dynamic_winch.json", "cases/g103a_club_sop_dynamic_winch.json"]:
        result = run_simulation(load_case(path))
        for drive, available in zip(result["history"]["drive_torque_nm"], result["history"]["winch_full_throttle_torque_nm"]):
            assert drive <= available + 1.0e-9


def test_post_release_terminal_time_is_exact(project_root):
    case = load_case(project_root / "cases" / "g103a_dynamic_winch.json")
    case["simulation"]["post_release_duration_s"] = 0.137
    result = run_simulation(case)
    release = next(event for event in result["events"] if event["code"] == "PILOT_COMMAND_RELEASE")
    assert result["metadata"]["final_time_s"] == pytest.approx(release["time"] + 0.137, abs=1.0e-12)
