import json
import math

from winch_launch_simulator.simulation import load_case, run_simulation


def canonical_paths(project_root):
    return [
        project_root / "cases" / "g103a_rwth_reference_prescribed_tension.json",
        project_root / "cases" / "g103a_club_prescribed_tension.json",
        project_root / "cases" / "g103a_rwth_reference_dynamic_winch_fem.json",
        project_root / "cases" / "g103a_club_dynamic_winch_fem.json",
    ]


def test_canonical_ground_start_cases_are_executable(project_root):
    for path in canonical_paths(project_root):
        case = load_case(path)
        case["simulation"]["max_time_s"] = 0.05
        result = run_simulation(case)
        expected = "GROUND_ROLL_TAILWHEEL" if case["pilot"]["sop"] == "club" else "RWTH_TRIM"
        assert result["history"]["pilot_phase"][0] == expected
        assert not [event for event in result["events"] if event["severity"] == "stop"]


def test_case_tree_is_two_plants_by_two_pilot_policies(project_root):
    paths = canonical_paths(project_root)
    assert all(path.exists() for path in paths)
    assert len(list((project_root / "cases").glob("*.json"))) == 4
    assert not any("dynamic_winch.json" == path.name for path in paths)


def test_no_removed_standard_or_elastic_secant_case_remains(project_root):
    for path in canonical_paths(project_root):
        case = load_case(path)
        assert case["pilot"]["sop"] != "standard"
        if case["cable"]["model"] == "ideal_secant":
            assert case["winch"]["model"] == "prescribed_tension"
            assert "stiffness_n_m" not in case["cable"]
            assert "damping_n_s_m" not in case["cable"]
        else:
            assert case["cable"]["model"] == "fem"
            assert case["winch"]["model"] == "dynamic"


def test_prescribed_pair_uses_identical_tension_profile(project_root):
    rwth = load_case(project_root / "cases" / "g103a_rwth_reference_prescribed_tension.json")
    club = load_case(project_root / "cases" / "g103a_club_prescribed_tension.json")
    assert rwth["winch"]["tension_profile"] == club["winch"]["tension_profile"]
    assert rwth["cable"] == club["cable"] == {"model": "ideal_secant"}


def test_each_pilot_policy_is_identical_across_plant_fidelity(project_root):
    rwth_low = load_case(project_root / "cases" / "g103a_rwth_reference_prescribed_tension.json")
    rwth_high = load_case(project_root / "cases" / "g103a_rwth_reference_dynamic_winch_fem.json")
    club_low = load_case(project_root / "cases" / "g103a_club_prescribed_tension.json")
    club_high = load_case(project_root / "cases" / "g103a_club_dynamic_winch_fem.json")
    assert rwth_low["pilot"] == rwth_high["pilot"]
    assert club_low["pilot"] == club_high["pilot"]


def test_rwth_reference_structure_matches_source_derived_features(project_root):
    case = load_case(project_root / "cases" / "g103a_rwth_reference_dynamic_winch_fem.json")
    pilot = case["pilot"]
    assert pilot["sop"] == "rwth_reference"
    assert pilot["trim_elevator_deg"] == -3.0
    assert pilot["safety_altitude_m"] == 50.0
    assert pilot["target_eas_m_s"] == 30.0
    assert case["winch"]["control"]["mode"] == "tension_control"
    schedule = case["winch"]["control"]["target_tension_schedule"]
    assert schedule["rise_rate_n_s"] == 500.0
    assert schedule["max_tension_n"] == 5000.0


def test_club_structure_is_shared_and_independent_of_removed_standard(project_root):
    case = load_case(project_root / "cases" / "g103a_club_dynamic_winch_fem.json")
    pilot = case["pilot"]
    assert pilot["sop"] == "club"
    assert pilot["hold_down"]["until_agl_m"] == 60.96
    assert pilot["rotation"]["pitch_target_deg"] == 45.0
    assert pilot["full_climb"]["pitch_target_deg"] == 45.0
    assert pilot["full_climb"]["pitch_kp"] == 0.25
    assert pilot["full_climb"]["elevator_bias_deg"] == -2.0
    assert case["winch"]["control"]["mode"] == "club_throttle"


def test_dynamic_fem_pair_uses_same_hardware_and_cable(project_root):
    rwth = load_case(project_root / "cases" / "g103a_rwth_reference_dynamic_winch_fem.json")
    club = load_case(project_root / "cases" / "g103a_club_dynamic_winch_fem.json")
    for key in ["x_m", "z_m", "drum_radius_m", "equivalent_inertia_kg_m2", "loss_torque_nm", "drive"]:
        assert rwth["winch"][key] == club["winch"][key]
    assert rwth["cable"] == club["cable"]


def test_dimensioned_g103a_geometry_is_common(project_root):
    case = load_case(canonical_paths(project_root)[0])
    assert case["hook"]["x_from_cg_m"] == 0.535
    assert case["hook"]["z_from_cg_m"] == 0.366
    contacts = {item["name"]: item for item in case["ground"]["aircraft"]["contacts"]}
    assert contacts["nosewheel"]["x_from_cg_m"] == 2.029
    assert contacts["mainwheel"]["x_from_cg_m"] == -0.169


def test_canonical_aerodynamic_provenance_is_wake15(project_root):
    for path in canonical_paths(project_root):
        case = load_case(path)
        state = case["aerodynamics"]["state"]
        assert state["aerodynamic_dataset_id"] == "g103a_stab_wake15_alpha2"
        assert state["reference_reynolds"] == 4400000.0
        assert state["wake_num_iter"] == 15
        assert state["num_wake_nodes"] == 64
        assert state["fixed_wake_flag"] is False
        assert state["drag_model"] == "CD=CD0+K*CL^2"
        assert case["aerodynamics"]["CD0"] > 0.0
        assert case["aerodynamics"]["K"] > 0.0


def test_drag_polar_fit_record_matches_canonical_cases(project_root):
    fit = json.loads((project_root / "aero_reference" / "G103A_drag_polar_fit.json").read_text(encoding="utf-8"))
    assert fit["model"] == "CD=CD0+K*CL^2"
    assert fit["r_squared"] > 0.9999
    for path in canonical_paths(project_root):
        case = load_case(path)
        assert math.isclose(case["aerodynamics"]["CD0"], fit["CD0"], abs_tol=1e-15)
        assert math.isclose(case["aerodynamics"]["K"], fit["K"], abs_tol=1e-15)


def test_provenance_has_no_removed_model_as_current_baseline(project_root):
    for path in canonical_paths(project_root):
        case = load_case(path)
        text = "\n".join(case["provenance"]["analysis_baseline"])
        assert "stiffness=18181" not in text
        assert "Standard SOP" not in text
        assert case["provenance"]["validation_reference"] == "validation/reference/canonical_regression.json"


def test_canonical_pitch_inertia_uses_amc_22481_rounded_value(project_root):
    for path in canonical_paths(project_root):
        case = load_case(path)
        assert case["aircraft"]["iyy_kg_m2"] == 2000.0


def test_ground_initializer_starts_on_mainwheel_and_tailwheel_only(project_root):
    from winch_launch_simulator.initialization import initialize_launch_state

    case = load_case(project_root / "cases" / "g103a_club_dynamic_winch_fem.json")
    _, metadata = initialize_launch_state(case)
    assert metadata["ground_support_mode"] == "mainwheel_tailwheel_two_point"
    assert metadata["mainwheel_normal_force_n"] > 0.0
    assert metadata["tailwheel_normal_force_n"] > 0.0
    assert metadata["nosewheel_normal_force_n"] == 0.0
    assert metadata["nosewheel_clearance_m"] > 0.0
