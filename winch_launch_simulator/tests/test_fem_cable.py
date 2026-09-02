import copy
import math

import numpy as np
import pytest

from winch_launch_simulator.initialization import initialize_launch_state
from winch_launch_simulator.models.fem_cable import (
    fem_cable_dynamics,
    fem_element_drag_forces,
    fem_reference_length,
    fem_winch_boundary_tension,
    fem_winch_reference_length,
    pack_fem_state,
    unpack_fem_state,
)
from winch_launch_simulator.simulation import BROKEN, RELEASED, TAUT, _apply_fem_node_removal, evaluate_dynamics, load_case, run_simulation


def make_cable(node_count=2, total_length=None):
    return {
        "model": "fem",
        "unstrained_length_m": float(total_length if total_length is not None else node_count),
        "fem": {
            "node_count": node_count,
            "diameter_m": 0.005,
            "linear_density_kg_m": 1.88e-2,
            "youngs_modulus_pa": 4.025e10,
            "strain_rate_stress_pa_s": 1.56e7,
            "include_weight": False,
            "include_aerodynamic_drag": False,
            "reeling": {"enabled": False},
        },
    }


def area(cable):
    return math.pi * cable["fem"]["diameter_m"] ** 2 / 4.0


def test_unstretched_and_compressed_fem_elements_carry_no_tension():
    cable = make_cable(node_count=2)
    for x in (1.0, 0.9):
        result = fem_cable_dynamics(
            cable,
            hook_position=np.array([0.0, 0.0]),
            hook_velocity=np.zeros(2),
            free_node_positions=np.array([[x, 0.0]]),
            free_node_velocities=np.zeros((1, 2)),
            winch_position=np.array([2.0, 0.0]),
            winch_tension_n=0.0,
        )
        assert result["hook_tension_n"] == pytest.approx(0.0, abs=1e-12)


def test_fem_elastic_and_damping_terms_match_constitutive_law():
    cable = make_cable(node_count=2)
    strain = 0.01
    elastic = fem_cable_dynamics(
        cable,
        np.array([0.0, 0.0]),
        np.zeros(2),
        np.array([[1.0 + strain, 0.0]]),
        np.zeros((1, 2)),
        np.array([2.0 + strain, 0.0]),
        0.0,
    )
    assert elastic["hook_tension_n"] == pytest.approx(
        area(cable) * cable["fem"]["youngs_modulus_pa"] * strain, rel=1e-12
    )

    rate = 0.2
    damped = fem_cable_dynamics(
        cable,
        np.array([0.0, 0.0]),
        np.zeros(2),
        np.array([[1.0, 0.0]]),
        np.array([[rate, 0.0]]),
        np.array([2.0, 0.0]),
        0.0,
    )
    assert damped["hook_tension_n"] == pytest.approx(
        area(cable) * cable["fem"]["strain_rate_stress_pa_s"] * rate, rel=1e-12
    )


def test_equal_internal_and_winch_tension_gives_zero_free_node_acceleration():
    cable = make_cable(node_count=2)
    strain = 0.01
    tension = area(cable) * cable["fem"]["youngs_modulus_pa"] * strain
    result = fem_cable_dynamics(
        cable,
        np.array([0.0, 0.0]),
        np.zeros(2),
        np.array([[1.0 + strain, 0.0]]),
        np.zeros((1, 2)),
        np.array([2.0 + strain, 0.0]),
        tension,
    )
    np.testing.assert_allclose(result["free_node_acceleration_m_s2"], 0.0, atol=1e-9)


def test_cable_weight_produces_gravity_acceleration_on_free_node():
    cable = make_cable(node_count=2)
    cable["fem"]["include_weight"] = True
    result = fem_cable_dynamics(
        cable,
        np.array([0.0, 0.0]),
        np.zeros(2),
        np.array([[1.0, 0.0]]),
        np.zeros((1, 2)),
        np.array([2.0, 0.0]),
        0.0,
    )
    np.testing.assert_allclose(result["free_node_acceleration_m_s2"], [[0.0, 9.80665]], atol=1e-12)


def test_fem_ground_contact_is_unilateral_and_pushes_upward():
    case = load_case("cases/g103a_rwth_reference_dynamic_winch_fem.json")
    cable = make_cable(node_count=2)
    cable["fem"]["include_weight"] = False
    result = fem_cable_dynamics(
        cable,
        np.array([0.0, -1.0]),
        np.zeros(2),
        np.array([[1.0, 0.001]]),
        np.zeros((1, 2)),
        np.array([2.0, 0.0]),
        0.0,
        ground=case["ground"],
    )
    assert result["ground_contact_node_count"] == 1
    assert result["max_ground_penetration_m"] > 0.0
    assert result["free_node_acceleration_m_s2"][0, 1] < 0.0


def test_ground_contact_damping_is_zero_before_penetration():
    case = load_case("cases/g103a_rwth_reference_dynamic_winch_fem.json")
    cable = make_cable(node_count=2)
    result = fem_cable_dynamics(
        cable,
        np.array([0.0, -1.0]),
        np.zeros(2),
        np.array([[1.0, -0.01]]),
        np.array([[0.0, 5.0]]),
        np.array([2.0, 0.0]),
        0.0,
        ground=case["ground"],
    )
    assert result["ground_contact_node_count"] == 0
    assert result["ground_normal_force_n"] == pytest.approx(0.0, abs=0.0)
    assert result["ground_friction_force_n"] == pytest.approx(0.0, abs=0.0)


def test_ground_contact_force_is_approximately_mesh_independent():
    ground = load_case("cases/g103a_rwth_reference_dynamic_winch_fem.json")["ground"]
    totals = []
    for node_count in (20, 40):
        cable = make_cable(node_count=node_count, total_length=1000.0)
        positions = np.column_stack(
            (
                np.linspace(1000.0 / node_count, 1000.0 * (node_count - 1) / node_count, node_count - 1),
                np.full(node_count - 1, 0.001),
            )
        )
        result = fem_cable_dynamics(
            cable,
            np.array([0.0, -1.0]),
            np.zeros(2),
            positions,
            np.zeros_like(positions),
            np.array([1000.0, 0.0]),
            0.0,
            ground=ground,
        )
        totals.append(result["ground_normal_force_n"])
    assert totals[1] == pytest.approx(totals[0], rel=0.04)


def test_wlo_element_drag_is_zero_at_zero_relative_wind_and_opposes_motion():
    velocity = np.array([[4.0, -2.0], [0.0, 10.0]])
    force = fem_element_drag_forces(
        0.005,
        1.3,
        np.array([[10.0, 0.0], [10.0, 0.0]]),
        np.array([10.0, 10.0]),
        velocity,
        np.array([1.225, 1.225]),
        np.array([[4.0, -2.0], [0.0, 0.0]]),
    )
    np.testing.assert_allclose(force[0], 0.0, atol=0.0)
    assert np.dot(force[1], velocity[1]) < 0.0


def test_winch_reference_length_includes_initial_reel_in_and_remesh_continuity():
    cable = make_cable(node_count=20, total_length=1000.0)
    before = fem_winch_reference_length(
        cable, drum_angle_rad=40.0, removed_node_count=0,
        initial_reel_in_m=3.0, drum_radius_m=0.5
    )
    after = fem_winch_reference_length(
        cable, drum_angle_rad=40.0, removed_node_count=1,
        initial_reel_in_m=3.0, drum_radius_m=0.5
    )
    assert before == pytest.approx(27.0)
    assert after - before == pytest.approx(fem_reference_length(cable))


def test_winch_boundary_tension_matches_elastic_strain():
    cable = make_cable(node_count=2)
    strain = 0.01
    result = fem_winch_boundary_tension(
        cable,
        np.array([0.0, 0.0]),
        np.zeros(2),
        np.array([1.0 + strain, 0.0]),
        1.0,
        0.0,
    )
    assert result["tension_n"] == pytest.approx(
        area(cable) * cable["fem"]["youngs_modulus_pa"] * strain, rel=1e-12
    )


def test_canonical_fem_initializer_solves_pretension_and_ground_equilibrium():
    case = load_case("cases/g103a_rwth_reference_dynamic_winch_fem.json")
    state, metadata = initialize_launch_state(case)
    positions, velocities = unpack_fem_state(case["cable"], state[9:])

    assert state.size == 85
    assert positions.shape == (19, 2)
    np.testing.assert_allclose(velocities, 0.0, atol=0.0)
    assert metadata["initial_winch_tension_n"] == pytest.approx(2500.0, abs=2.0)
    assert metadata["initial_hook_tension_n"] == pytest.approx(2500.0, abs=2.0)
    assert metadata["fem_static_max_free_node_force_n"] < 0.1
    assert metadata["fem_static_max_ground_penetration_m"] < 0.003
    assert 3.0 < metadata["initial_reel_in_m"] < 5.0


def test_prelaunch_equilibrium_removes_the_old_large_initial_tension_pulse():
    case = copy.deepcopy(load_case("cases/g103a_rwth_reference_dynamic_winch_fem.json"))
    case["simulation"]["max_time_s"] = 0.2
    result = run_simulation(case)
    tension = result["history"]["hook_tension_n"]
    assert min(tension) > 2000.0
    assert max(tension) < 3000.0


@pytest.mark.parametrize("connection_state", [RELEASED, BROKEN])
def test_released_or_broken_fem_cable_applies_no_aircraft_cable_force(connection_state):
    case = load_case("cases/g103a_rwth_reference_dynamic_winch_fem.json")
    state, metadata = initialize_launch_state(case)
    runtime = {
        "initial_reel_in_m": metadata["initial_reel_in_m"],
        "removed_node_count": 0,
        "initial_node_count": case["cable"]["fem"]["node_count"],
    }
    derivative, diagnostics = evaluate_dynamics(
        0.0,
        state,
        case,
        connection_state,
        {"liftoff_time_s": None, "rotation_start_time_s": None, "eas_rate_m_s2": 0.0},
        runtime,
    )
    assert diagnostics["cable_tension_n"] == 0.0
    assert diagnostics["hook_moment_nm"] == 0.0
    assert diagnostics["cable_power_to_aircraft_w"] == 0.0
    np.testing.assert_allclose(derivative[9:], 0.0, atol=0.0)


def test_fem_node_removal_preserves_remaining_states_reference_length_and_mass():
    case = load_case("cases/g103a_rwth_reference_dynamic_winch_fem.json")
    state, metadata = initialize_launch_state(case)
    positions, velocities = unpack_fem_state(case["cable"], state[9:])
    kept_position = positions[0].copy()
    kept_velocity = velocities[0].copy()

    positions[-1] = np.array([1095.0, -0.0025])
    state[9 : 9 + positions.size] = positions.ravel()
    runtime = {
        "initial_reel_in_m": metadata["initial_reel_in_m"],
        "removed_node_count": 0,
        "initial_node_count": case["cable"]["fem"]["node_count"],
    }
    pilot_state = {"liftoff_time_s": None, "rotation_start_time_s": None, "eas_rate_m_s2": 0.0}
    _, diagnostics = evaluate_dynamics(0.0, state, case, TAUT, pilot_state, runtime)
    new_state, new_diagnostics, events = _apply_fem_node_removal(
        0.0, state, case, TAUT, pilot_state, runtime, diagnostics
    )
    new_positions, new_velocities = unpack_fem_state(case["cable"], new_state[9:])

    assert len(events) == 1
    values = events[0]["values"]
    assert values["free_reference_length_after_m"] == pytest.approx(values["free_reference_length_before_m"], abs=1e-12)
    assert values["free_cable_mass_after_kg"] == pytest.approx(values["free_cable_mass_before_kg"], abs=1e-12)
    np.testing.assert_allclose(new_positions[0], kept_position, atol=0.0)
    np.testing.assert_allclose(new_velocities[0], kept_velocity, atol=0.0)
    assert new_diagnostics["fem_active_node_count"] == 19
    assert new_state.size == state.size - 4


def test_repository_fem_case_reports_ground_contact_diagnostics():
    case = copy.deepcopy(load_case("cases/g103a_rwth_reference_dynamic_winch_fem.json"))
    case["simulation"]["max_time_s"] = 0.05
    result = run_simulation(case)
    history = result["history"]
    assert history["fem_ground_contact_node_count"][0] > 0
    assert history["fem_max_ground_penetration_m"][0] < 0.003
    assert history["fem_cable_weight_n"][0] > 0.0


def test_ground_start_fem_has_no_large_artificial_initial_tension_collapse(project_root):
    case = load_case(project_root / "cases" / "g103a_rwth_reference_dynamic_winch_fem.json")
    case["simulation"]["max_time_s"] = 0.2
    result = run_simulation(case)
    tension = result["history"]["cable_tension_n"]
    assert tension[0] == pytest.approx(2500.0, abs=2.0)
    assert min(tension) > 2300.0
    assert max(tension) < 2600.0
