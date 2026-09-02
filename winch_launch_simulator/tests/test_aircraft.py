import math

from winch_launch_simulator.models.aircraft import air_data, aerodynamic_forces_and_moment
from winch_launch_simulator.simulation import load_case


def test_no_wind_air_data_matches_body_state():
    air = air_data(20.0, 2.0, math.radians(5.0), 0.0, 0.0, 1.225)
    assert math.isclose(air["airspeed_m_s"], math.hypot(20.0, 2.0), rel_tol=1e-12)
    assert math.isclose(air["alpha_rad"], math.atan2(2.0, 20.0), rel_tol=1e-12)




def test_drag_coefficient_uses_parabolic_polar_at_instantaneous_cl():
    case = load_case("cases/g103a_rwth_reference_prescribed_tension.json")
    air = {
        "airspeed_m_s": 30.0,
        "alpha_rad": math.radians(4.0),
    }
    result = aerodynamic_forces_and_moment(
        case, air, math.radians(2.0), math.radians(-3.0), 1.225
    )
    aero = case["aerodynamics"]
    expected_cd = aero["CD0"] + aero["K"] * result["CL"] ** 2
    assert math.isclose(result["CD"], expected_cd, rel_tol=1e-12, abs_tol=1e-15)
    assert result["CD"] >= aero["CD0"]

def test_positive_elevator_command_gives_more_nose_down_cm():
    case = load_case("cases/g103a_rwth_reference_prescribed_tension.json")
    air = {
        "airspeed_m_s": 30.0,
        "alpha_rad": math.radians(2.0),
    }
    negative = aerodynamic_forces_and_moment(case, air, 0.0, math.radians(-1.0), 1.225)
    positive = aerodynamic_forces_and_moment(case, air, 0.0, math.radians(+1.0), 1.225)
    assert positive["Cm"] < negative["Cm"]


def test_wheel_ground_contact_is_unilateral_and_has_no_attractive_normal_force():
    from winch_launch_simulator.models.aircraft import wheel_ground_forces_and_moment

    case = load_case("cases/g103a_rwth_reference_dynamic_winch_fem.json")
    # Far above ground: no tire contact.
    clear = wheel_ground_forces_and_moment(case, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0)
    assert clear["contact_count"] == 0
    assert clear["normal_force_n"] == 0.0

    # Lower the aircraft enough to penetrate the tire penalty surface.
    contact = wheel_ground_forces_and_moment(case, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0)
    assert contact["contact_count"] > 0
    assert contact["normal_force_n"] > 0.0
    assert contact["force_ground_n"][1] < 0.0
