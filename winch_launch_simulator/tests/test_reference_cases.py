from winch_launch_simulator.simulation import load_case, run_simulation


def canonical_paths(project_root):
    return [
        project_root / "cases" / "g103a_prescribed_tension.json",
        project_root / "cases" / "g103a_dynamic_winch.json",
        project_root / "cases" / "g103a_club_sop_dynamic_winch.json",
        project_root / "cases" / "g103a_dynamic_winch_fem.json",
        project_root / "cases" / "g103a_club_sop_dynamic_winch_fem.json",
    ]


def test_canonical_ground_start_cases_are_executable(project_root):
    for path in canonical_paths(project_root):
        case = load_case(path)
        case["simulation"]["max_time_s"] = 0.05
        result = run_simulation(case)
        expected = "GROUND_ROLL_TAILWHEEL" if case["pilot"]["sop"] == "club" else "GROUND_ROLL"
        assert result["history"]["pilot_phase"][0] == expected
        assert not [event for event in result["events"] if event["severity"] == "stop"]


def test_standard_and_club_pilot_models_are_intentionally_distinct(project_root):
    standard = load_case(project_root / "cases" / "g103a_dynamic_winch.json")
    club = load_case(project_root / "cases" / "g103a_club_sop_dynamic_winch.json")
    assert standard["pilot"]["sop"] == "standard"
    assert club["pilot"]["sop"] == "club"
    assert "target_eas_m_s" in standard["pilot"]["full_climb"]
    assert "target_eas_m_s" not in club["pilot"]["full_climb"]
    assert club["pilot"]["hold_down"]["until_agl_m"] == 60.96
    assert club["pilot"]["full_climb"]["pitch_target_deg"] == 45.0


def test_case_tree_matches_three_level_fidelity_hierarchy(project_root):
    assert not (project_root / "cases" / "optimized").exists()
    assert not (project_root / "cases" / "fem" / "g103a_prescribed_tension_fem.json").exists()
    assert all(path.exists() for path in canonical_paths(project_root))


def test_case_schema_uses_1100m_club_launch_setup_and_derived_prelaunch_reel_in(project_root):
    secant = load_case(project_root / "cases" / "g103a_dynamic_winch.json")
    fem = load_case(project_root / "cases" / "g103a_dynamic_winch_fem.json")
    for case in (secant, fem):
        assert case["launch_setup"]["mode"] == "ground_start"
        assert case["launch_setup"]["pretension_n"] == 2500.0
        assert case["winch"]["x_m"] == 1100.0
        assert case["cable"]["unstrained_length_m"] == 1100.0
        assert "initial_reel_in_m" not in case["cable"]
        assert "initial_drum_rate_rad_s" not in case["winch"]


def test_measured_g103a_hook_and_wheel_center_geometry_is_canonical(project_root):
    case = load_case(project_root / "cases" / "g103a_dynamic_winch.json")
    assert case["hook"]["x_from_cg_m"] == 0.535
    assert case["hook"]["z_from_cg_m"] == 0.366
    contacts = {item["name"]: item for item in case["ground"]["aircraft"]["contacts"]}
    assert contacts["nosewheel"]["x_from_cg_m"] == 2.029
    assert contacts["nosewheel"]["z_from_cg_m"] == 0.508
    assert contacts["mainwheel"]["x_from_cg_m"] == -0.169
    assert contacts["mainwheel"]["z_from_cg_m"] == 0.366


def test_v018_canonical_control_structure(project_root):
    prescribed = load_case(project_root / "cases" / "g103a_prescribed_tension.json")
    secant_standard = load_case(project_root / "cases" / "g103a_dynamic_winch.json")
    secant_club = load_case(project_root / "cases" / "g103a_club_sop_dynamic_winch.json")
    fem_standard = load_case(project_root / "cases" / "g103a_dynamic_winch_fem.json")
    fem_club = load_case(project_root / "cases" / "g103a_club_sop_dynamic_winch_fem.json")

    assert prescribed["pilot"]["sop"] == "standard"
    assert prescribed["pilot"]["ground_roll_elevator_deg"] == -5.2
    assert prescribed["winch"]["tension_schedule"]["rotation_rise_rate_n_s"] == 90.0

    assert secant_standard["pilot"]["sop"] == "standard"
    assert secant_standard["pilot"]["pitch_rate_gain_s"] == 2.4
    assert secant_standard["pilot"]["full_climb"]["target_eas_m_s"] == 27.8
    assert "eas_ki" not in secant_standard["pilot"]["full_climb"]

    assert fem_standard["pilot"]["sop"] == "standard"
    assert fem_standard["winch"]["control"]["target_tension_schedule"]["max_tension_n"] == 5000.0

    for club in (secant_club, fem_club):
        assert club["pilot"]["sop"] == "club"
        assert club["pilot"]["ground_roll"]["tailwheel_load_cue_n"] == 500.0
        assert club["pilot"]["hold_down"]["low_tension_cue_n"] == 4000.0
        assert club["pilot"]["rotation"]["pitch_target_deg"] == 45.0
        assert club["pilot"]["full_climb"]["pitch_target_deg"] == 45.0
        assert club["pilot"]["full_climb"]["pitch_kp"] == 0.5
        assert club["winch"]["control"]["mode"] == "club_throttle"
        assert club["winch"]["control"]["launch_throttle"] == 1.0
        assert club["winch"]["control"]["launch_power_ramp_duration_s"] == 3.5
        assert "ground_roll_tension_hold" not in club["winch"]["control"]
        assert club["winch"]["control"]["terminal_fade"]["start_ground_angle_deg"] == 60.0
        assert club["winch"]["control"]["terminal_fade"]["end_ground_angle_deg"] == 75.0

    for case in (prescribed, secant_standard, secant_club, fem_standard, fem_club):
        specification = case["provenance"]["specification"]
        assert specification == "docs/implementation_spec.md"
        assert (project_root / specification).exists()


def test_v018_club_secant_and_fem_share_the_final_common_sop(project_root):
    secant = load_case(project_root / "cases" / "g103a_club_sop_dynamic_winch.json")
    fem = load_case(project_root / "cases" / "g103a_club_sop_dynamic_winch_fem.json")
    for case in (secant, fem):
        assert case["pilot"]["full_climb"]["pitch_kp"] == 0.5
        assert case["pilot"]["full_climb"]["elevator_bias_deg"] == -4.0
        assert case["pilot"]["full_climb"]["power_reduction_request_margin_m_s"] == 0.9
        assert case["winch"]["control"]["launch_throttle"] == 1.0


def test_v018_case_provenance_contains_no_stale_club_or_1000m_baseline(project_root):
    stale = ("1000 m", "2700 N", "12 deg", "throttle=0.95", "copied from")
    for path in canonical_paths(project_root):
        case = load_case(path)
        text = "\n".join(case["provenance"]["analysis_baseline"])
        assert not any(token in text for token in stale)
        validation_reference = case["provenance"]["validation_reference"]
        assert validation_reference == "validation/reference/canonical_regression.json"
        assert (project_root / validation_reference).exists()
