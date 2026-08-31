import math

from winch_launch_simulator.models.winch import drum_derivatives, full_throttle_drive_torque
from winch_launch_simulator.simulation import load_case


def test_dynamic_winch_energy_identity():
    case = load_case("cases/g103a_dynamic_winch.json")
    tension = 4000.0
    drive = 2500.0
    omega = 40.0
    _, omega_dot = drum_derivatives(case, tension, drive, omega)
    winch = case["winch"]
    left = winch["equivalent_inertia_kg_m2"] * omega * omega_dot
    right = drive * omega - tension * winch["drum_radius_m"] * omega - winch["loss_torque_nm"] * omega
    assert math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-9)


def test_full_throttle_drive_torque_interpolates_shared_curve():
    winch = load_case("cases/g103a_dynamic_winch.json")["winch"]
    assert full_throttle_drive_torque(winch, 25.0) == 3200.0
    assert full_throttle_drive_torque(winch, 40.0) == 2300.0


def test_standard_and_club_cases_use_the_same_dynamic_winch_hardware():
    standard = load_case("cases/g103a_dynamic_winch.json")
    club = load_case("cases/g103a_club_sop_dynamic_winch.json")
    for key in ["x_m", "z_m", "drum_radius_m", "equivalent_inertia_kg_m2", "loss_torque_nm", "drive"]:
        assert standard["winch"][key] == club["winch"][key]
    for key in ["stiffness_n_m", "damping_n_s_m", "unstrained_length_m"]:
        assert standard["cable"][key] == club["cable"][key]
