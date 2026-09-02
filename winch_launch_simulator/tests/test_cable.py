import math

import pytest

from winch_launch_simulator.models.cable import hook_winch_geometry
from winch_launch_simulator.models.winch import prescribed_tension, rwth_target_tension


def test_hook_winch_geometry_points_from_hook_to_winch():
    hook = {"x_from_cg_m": 0.0, "z_from_cg_m": 0.0}
    winch = {"x_m": 100.0, "z_m": 0.0}
    geometry = hook_winch_geometry(0.0, -10.0, 0.0, 0.0, 20.0, 0.0, hook, winch)
    assert geometry["unit_x_ground"] > 0.0
    assert geometry["unit_z_ground"] > 0.0
    assert math.isclose(geometry["length_m"], math.hypot(100.0, 10.0), rel_tol=1e-12)
    assert geometry["length_rate_m_s"] < 0.0


def test_common_prescribed_tension_is_liftoff_relative_and_linear_between_points():
    settings = {
        "time_points_after_liftoff": [
            {"time_s": 0.0, "tension_n": 2500.0},
            {"time_s": 1.5, "tension_n": 1600.0},
            {"time_s": 18.0, "tension_n": 4750.0},
        ],
        "terminal_fade": {
            "start_ground_angle_deg": 60.0,
            "end_ground_angle_deg": 75.0,
            "final_tension_n": 400.0,
        },
    }
    pretension = 2500.0
    on_ground, progress = prescribed_tension(settings, pretension, 10.0, 0.0, None)
    at_liftoff, _ = prescribed_tension(settings, pretension, 10.0, 0.0, 10.0)
    midpoint, _ = prescribed_tension(settings, pretension, 10.75, 0.0, 10.0)
    assert on_ground == pytest.approx(pretension)
    assert at_liftoff == pytest.approx(2500.0)
    assert midpoint == pytest.approx(2050.0)
    assert progress == 0.0


def test_common_prescribed_tension_terminal_fade_uses_cable_angle():
    settings = {
        "time_points_after_liftoff": [
            {"time_s": 0.0, "tension_n": 2500.0},
            {"time_s": 1.5, "tension_n": 1600.0},
            {"time_s": 18.0, "tension_n": 4750.0},
        ],
        "terminal_fade": {
            "start_ground_angle_deg": 60.0,
            "end_ground_angle_deg": 75.0,
            "final_tension_n": 400.0,
        },
    }
    target, progress = prescribed_tension(
        settings,
        2500.0,
        30.0,
        math.radians(67.5),
        10.0,
    )
    assert progress == pytest.approx(0.5)
    assert target == pytest.approx((4750.0 + 400.0) / 2.0)


def test_rwth_target_tension_uses_source_structure_with_g103a_adapted_cap():
    settings = {
        "rise_rate_n_s": 500.0,
        "max_tension_n": 5000.0,
        "terminal_fade": {
            "start_ground_angle_deg": 65.0,
            "end_ground_angle_deg": 75.0,
            "final_tension_n": 400.0,
        },
    }
    before_liftoff, _ = rwth_target_tension(settings, 2500.0, 5.0, 0.0, None)
    rising, _ = rwth_target_tension(settings, 2500.0, 7.0, 0.0, 5.0)
    capped, _ = rwth_target_tension(settings, 2500.0, 20.0, 0.0, 5.0)
    assert before_liftoff == pytest.approx(2500.0)
    assert rising == pytest.approx(3500.0)
    assert capped == pytest.approx(5000.0)
