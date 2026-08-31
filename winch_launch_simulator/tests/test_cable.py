import math

from winch_launch_simulator.models.cable import hook_winch_geometry, secant_elastic_tension
from winch_launch_simulator.models.winch import scheduled_tension
from winch_launch_simulator.simulation import load_case


def test_hook_winch_geometry_points_from_hook_to_winch():
    hook = {"x_from_cg_m": 0.0, "z_from_cg_m": 0.0}
    winch = {"x_m": 100.0, "z_m": 0.0}
    geometry = hook_winch_geometry(0.0, -10.0, 0.0, 0.0, 20.0, 0.0, hook, winch)
    assert geometry["unit_x_ground"] > 0.0
    assert geometry["unit_z_ground"] > 0.0
    assert math.isclose(geometry["length_m"], math.hypot(100.0, 10.0), rel_tol=1e-12)
    assert geometry["length_rate_m_s"] < 0.0


def test_secant_elastic_tension_uses_unstrained_length_and_derived_reel_in():
    case = load_case("cases/g103a_dynamic_winch.json")
    cable = case["cable"]
    winch = case["winch"]
    reel_in = 0.5
    deployed = cable["unstrained_length_m"] - reel_in
    tension = secant_elastic_tension(
        cable,
        winch,
        deployed + 1.0,
        0.0,
        0.0,
        0.0,
        reel_in,
    )
    assert tension == cable["stiffness_n_m"]


def test_tension_schedule_holds_pretension_until_liftoff_then_uses_rotation_rate():
    settings = {
        "initial_rise_rate_n_s": 10.0,
        "rotation_rise_rate_n_s": 100.0,
        "max_tension_n": 5000.0,
        "release_prepare_angle_deg": 80.0,
        "release_angle_deg": 85.0,
        "release_target_n": 400.0,
    }
    pretension = 1000.0

    on_ground = scheduled_tension(settings, pretension, 5.0, 0.0, None, None)
    before_rotation = scheduled_tension(settings, pretension, 7.0, 0.0, 2.0, None)
    after_rotation = scheduled_tension(settings, pretension, 9.0, 0.0, 2.0, 7.0)

    assert on_ground == 1000.0
    assert before_rotation == 1050.0
    assert after_rotation == 1250.0
