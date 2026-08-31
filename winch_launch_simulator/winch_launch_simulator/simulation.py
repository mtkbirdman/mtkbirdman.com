import copy
import csv
import json
import math
from pathlib import Path

import numpy as np

from .initialization import initialize_launch_state
from .models.aircraft import air_data, aerodynamic_forces_and_moment, wheel_ground_forces_and_moment
from .models.atmosphere import G, atmosphere
from .models.cable import hook_winch_geometry, secant_elastic_tension
from .models.controllers import club_pilot_control, standard_pilot_control, tension_control_drive_torque
from .models.fem_cable import (
    fem_cable_dynamics,
    fem_node_count_from_state,
    fem_reference_length,
    fem_winch_boundary_tension,
    fem_winch_reference_length,
    pack_fem_state,
    unpack_fem_state,
)
from .models.winch import club_winch_throttle, drum_derivatives, full_throttle_drive_torque, scheduled_tension

TAUT = 1
SLACK = 2
RELEASED = 3
BROKEN = 4

BASE_STATE_NAMES = (
    "x_m",
    "z_m",
    "u_m_s",
    "w_m_s",
    "theta_rad",
    "q_rad_s",
    "drum_angle_rad",
    "drum_rate_rad_s",
    "winch_integral",
)

def state_names(case, state=None):
    """Return state names, including the currently active variable-size FEM suffix."""
    names = list(BASE_STATE_NAMES)
    if case["cable"]["model"] != "fem":
        return tuple(names)

    if state is None:
        node_count = case["cable"]["fem"]["node_count"]
    else:
        node_count = fem_node_count_from_state(np.asarray(state)[len(BASE_STATE_NAMES) :])

    free_node_count = node_count - 1
    for index in range(free_node_count):
        names.extend((f"cable_node_{index + 2}_x_m", f"cable_node_{index + 2}_z_m"))
    for index in range(free_node_count):
        names.extend(
            (
                f"cable_node_{index + 2}_vx_m_s",
                f"cable_node_{index + 2}_vz_m_s",
            )
        )
    return tuple(names)

def load_case(path):
    """Load one JSON case and validate the physical choices needed by the solver."""
    with open(path, "r", encoding="utf-8") as file:
        case = json.load(file)

    required = {
        "case_id",
        "simulation",
        "launch_setup",
        "ground",
        "aircraft",
        "aerodynamics",
        "hook",
        "cable",
        "winch",
        "pilot",
        "atmosphere",
        "limits",
        "output",
        "provenance",
    }
    missing = sorted(required - set(case))
    if missing:
        raise ValueError(f"Missing case sections: {', '.join(missing)}")

    if case["simulation"]["dt_s"] <= 0.0:
        raise ValueError("simulation.dt_s must be positive.")
    if case["aircraft"]["mass_kg"] <= 0.0 or case["aircraft"]["iyy_kg_m2"] <= 0.0:
        raise ValueError("aircraft mass and pitch inertia must be positive.")
    if case["aerodynamics"]["reference_area_m2"] <= 0.0:
        raise ValueError("aerodynamics.reference_area_m2 must be positive.")

    launch_setup = case["launch_setup"]
    if launch_setup.get("mode") != "ground_start":
        raise ValueError("launch_setup.mode must be 'ground_start'.")
    if launch_setup.get("pretension_n", 0.0) <= 0.0:
        raise ValueError("launch_setup.pretension_n must be positive.")

    ground = case["ground"]
    aircraft_ground = ground.get("aircraft", {})
    contacts = aircraft_ground.get("contacts", [])
    if {item.get("name") for item in contacts} != {"nosewheel", "mainwheel", "tailwheel"}:
        raise ValueError("ground.aircraft.contacts must define nosewheel, mainwheel and tailwheel.")
    if aircraft_ground.get("normal_stiffness_n_m", 0.0) <= 0.0:
        raise ValueError("ground.aircraft.normal_stiffness_n_m must be positive.")
    if aircraft_ground.get("normal_damping_n_s_m", 0.0) < 0.0:
        raise ValueError("ground.aircraft.normal_damping_n_s_m must be non-negative.")
    if aircraft_ground.get("rolling_resistance_coefficient", -1.0) < 0.0:
        raise ValueError("ground.aircraft.rolling_resistance_coefficient must be non-negative.")
    if any(item.get("radius_m", 0.0) <= 0.0 for item in contacts):
        raise ValueError("All wheel radii must be positive.")

    cable = case["cable"]
    winch = case["winch"]
    if cable.get("model") not in {"secant", "fem"}:
        raise ValueError("cable.model must be 'secant' or 'fem'.")
    if winch.get("model") not in {"prescribed_tension", "dynamic"}:
        raise ValueError("winch.model must be 'prescribed_tension' or 'dynamic'.")

    if cable["model"] == "secant" and winch["model"] == "dynamic":
        for key in ("stiffness_n_m", "damping_n_s_m", "unstrained_length_m"):
            if cable.get(key, 0.0) <= 0.0:
                raise ValueError(f"cable.{key} must be positive for a dynamic Secant case.")

    if cable["model"] == "fem":
        if cable.get("unstrained_length_m", 0.0) <= 0.0:
            raise ValueError("cable.unstrained_length_m must be positive for a FEM case.")
        fem = cable.get("fem")
        if fem is None:
            raise ValueError("cable.fem is required when cable.model='fem'.")
        required_fem = {
            "node_count",
            "diameter_m",
            "linear_density_kg_m",
            "youngs_modulus_pa",
            "strain_rate_stress_pa_s",
            "include_weight",
            "include_aerodynamic_drag",
            "reeling",
        }
        missing_fem = sorted(required_fem - set(fem))
        if missing_fem:
            raise ValueError(f"Missing cable.fem fields: {', '.join(missing_fem)}")
        if not isinstance(fem["node_count"], int) or fem["node_count"] < 2:
            raise ValueError("cable.fem.node_count must be an integer of at least 2.")
        for key in (
            "diameter_m",
            "linear_density_kg_m",
            "youngs_modulus_pa",
            "strain_rate_stress_pa_s",
        ):
            if fem[key] <= 0.0:
                raise ValueError(f"cable.fem.{key} must be positive.")
        if not isinstance(fem["include_weight"], bool) or not isinstance(
            fem["include_aerodynamic_drag"], bool
        ):
            raise ValueError("FEM weight and aerodynamic-drag switches must be boolean.")
        if fem["include_aerodynamic_drag"]:
            if fem.get("drag_coefficient", 0.0) <= 0.0:
                raise ValueError(
                    "cable.fem.drag_coefficient must be positive when aerodynamic drag is enabled."
                )
            if fem.get("aerodynamic_model") != "williams_lansdorp_ockels_2007_cylinder":
                raise ValueError(
                    "cable.fem.aerodynamic_model must be "
                    "'williams_lansdorp_ockels_2007_cylinder'."
                )

        cable_ground = ground.get("cable", {})
        if cable_ground.get("normal_stiffness_n_m2", 0.0) <= 0.0:
            raise ValueError("ground.cable.normal_stiffness_n_m2 must be positive for FEM cases.")
        if cable_ground.get("normal_damping_n_s_m2", 0.0) < 0.0:
            raise ValueError("ground.cable.normal_damping_n_s_m2 must be non-negative.")
        if cable_ground.get("sliding_friction_coefficient", -1.0) < 0.0:
            raise ValueError("ground.cable.sliding_friction_coefficient must be non-negative.")

        reeling = fem["reeling"]
        if not isinstance(reeling.get("enabled"), bool):
            raise ValueError("cable.fem.reeling.enabled must be boolean.")
        if reeling["enabled"]:
            factor = reeling.get("node_removal_length_factor")
            minimum = reeling.get("minimum_node_count")
            if factor is None or not 0.0 < factor < 1.0:
                raise ValueError(
                    "cable.fem.reeling.node_removal_length_factor must be between 0 and 1."
                )
            if not isinstance(minimum, int) or not 2 <= minimum <= fem["node_count"]:
                raise ValueError(
                    "cable.fem.reeling.minimum_node_count must be an integer from 2 through node_count."
                )
        if winch["model"] == "dynamic" and not reeling["enabled"]:
            raise ValueError(
                "cable.model='fem' with winch.model='dynamic' requires "
                "cable.fem.reeling.enabled=true."
            )

    if winch["model"] == "dynamic":
        control = winch["control"]
        if winch["drive"]["model"] != "equivalent_full_throttle_torque_curve":
            raise ValueError(
                "winch.drive.model must be 'equivalent_full_throttle_torque_curve'."
            )
        if control["mode"] not in {"tension_control", "club_throttle"}:
            raise ValueError(
                "winch.control.mode must be 'tension_control' or 'club_throttle'."
            )
        curve = winch["drive"]["full_throttle_torque_curve"]
        if len(curve) < 2:
            raise ValueError("winch.drive.full_throttle_torque_curve requires at least two points.")
        rates = [point["drum_rate_rad_s"] for point in curve]
        torques = [point["drive_torque_nm"] for point in curve]
        if any(current <= previous for previous, current in zip(rates, rates[1:])):
            raise ValueError("winch full-throttle curve drum rates must increase strictly.")
        if rates[0] < 0.0 or any(torque < 0.0 for torque in torques):
            raise ValueError("winch full-throttle curve rates and torques must be non-negative.")
        if control["mode"] == "club_throttle":
            if not 0.0 <= control["launch_throttle"] <= 1.0:
                raise ValueError("winch.control.launch_throttle must be between 0 and 1.")
            if control["launch_power_ramp_duration_s"] <= 0.0:
                raise ValueError("winch.control.launch_power_ramp_duration_s must be positive.")
            zero_speed_full_torque = full_throttle_drive_torque(winch, 0.0)
            pretension_drive_torque = (
                launch_setup["pretension_n"] * winch["drum_radius_m"]
                + winch["loss_torque_nm"]
            )
            if zero_speed_full_torque <= 0.0 or pretension_drive_torque > zero_speed_full_torque:
                raise ValueError(
                    "Club pretension equilibrium must be supportable by the zero-speed drive torque."
                )
            speed = control["speed_protection"]
            if speed["operator_response_delay_s"] < 0.0 or not 0.0 <= speed["reduced_throttle"] <= 1.0:
                raise ValueError("Invalid Club winch speed-protection settings.")
            terminal = control["terminal_fade"]
            if terminal["end_ground_angle_deg"] <= terminal["start_ground_angle_deg"]:
                raise ValueError("Club terminal-fade end angle must exceed its start angle.")
            if not 0.0 <= terminal["final_throttle"] <= 1.0:
                raise ValueError("Club terminal-fade final_throttle must be between 0 and 1.")

    pilot = case["pilot"]
    if pilot.get("sop") not in {"standard", "club"}:
        raise ValueError("pilot.sop must be 'standard' or 'club'.")

    if pilot["sop"] == "standard":
        initial = pilot["initial_climb"]
        if initial["pitch_kp"] < 0.0:
            raise ValueError("pilot.initial_climb.pitch_kp must be non-negative.")
        rotation = pilot["rotation"]
        if rotation["min_duration_s"] <= 0.0 or rotation["start_min_agl_m"] < 0.0:
            raise ValueError("Invalid standard-pilot rotation settings.")
        if rotation["full_climb_min_agl_m"] <= 0.0:
            raise ValueError("pilot.rotation.full_climb_min_agl_m must be positive.")
        release_preparation = pilot["release_preparation"]
        if release_preparation["end_ground_angle_deg"] <= release_preparation["start_ground_angle_deg"]:
            raise ValueError("Standard release-preparation end angle must exceed its start angle.")
    else:
        if winch["model"] != "dynamic" or winch["control"]["mode"] != "club_throttle":
            raise ValueError("pilot.sop='club' requires dynamic winch control mode 'club_throttle'.")
        if pilot["ground_roll"]["tailwheel_load_cue_n"] <= 0.0:
            raise ValueError("Club tailwheel_load_cue_n must be positive.")
        hold = pilot["hold_down"]
        if hold["until_agl_m"] <= 0.0 or hold["low_tension_cue_n"] <= 0.0:
            raise ValueError("Invalid Club hold-down settings.")
        if pilot["rotation"]["duration_s"] <= 0.0:
            raise ValueError("Club rotation duration_s must be positive.")
        full = pilot["full_climb"]
        if full["pitch_kp"] < 0.0 or full["speed_protection_margin_m_s"] < 0.0:
            raise ValueError("Invalid Club full-climb settings.")
        if full["power_reduction_request_margin_m_s"] < 0.0:
            raise ValueError("Club power-reduction request margin must be non-negative.")

    if case["limits"]["initial_phase_max_winch_speed_factor"] < 1.0:
        raise ValueError("limits.initial_phase_max_winch_speed_factor must be at least 1.0.")

    return case

def _fem_winch_boundary_state(case, fem_positions, fem_velocities, drum_angle, drum_rate, cable_runtime):
    """Return variable-length final-element geometry and tension for a dynamic winch."""
    winch = case["winch"]
    removed = int(cable_runtime.get("removed_node_count", 0))
    final_reference_length = fem_winch_reference_length(
        case["cable"],
        drum_angle,
        removed,
        initial_reel_in_m=float(cable_runtime.get("initial_reel_in_m", 0.0)),
        drum_radius_m=winch["drum_radius_m"],
    )
    final_reference_rate = -winch["drum_radius_m"] * drum_rate
    boundary = fem_winch_boundary_tension(
        case["cable"],
        fem_positions[-1],
        fem_velocities[-1],
        np.array([winch["x_m"], winch.get("z_m", 0.0)]),
        final_reference_length,
        final_reference_rate,
    )
    boundary["reference_length_m"] = final_reference_length
    boundary["reference_length_rate_m_s"] = final_reference_rate
    return boundary

def evaluate_dynamics(
    time_s,
    state,
    case,
    connection_state,
    pilot_state,
    cable_runtime=None,
    *,
    collect_diagnostics=True,
):
    """Evaluate one complete coupled aircraft/cable/winch state in physical order."""
    if cable_runtime is None:
        cable_runtime = {}
        if case["cable"]["model"] == "fem":
            cable_runtime = {
                "removed_node_count": 0,
                "initial_node_count": case["cable"]["fem"]["node_count"],
            }

    (
        x,
        z,
        u,
        w,
        theta,
        pitch_rate,
        drum_angle,
        drum_rate,
        winch_integral,
    ) = state[: len(BASE_STATE_NAMES)]

    altitude = -z
    agl = case["ground"]["z_m"] - z

    # 1. Atmosphere and aerodynamic state.
    rho, wind_x, wind_z = atmosphere(case, altitude)
    air = air_data(u, w, theta, wind_x, wind_z, rho)
    air["density_kg_m3"] = rho

    # 2. Rigid-body kinematics and hook/winch geometry.
    x_dot = u * math.cos(theta) + w * math.sin(theta)
    z_dot = -u * math.sin(theta) + w * math.cos(theta)
    cable_geometry = hook_winch_geometry(
        x,
        z,
        theta,
        pitch_rate,
        x_dot,
        z_dot,
        case["hook"],
        case["winch"],
    )
    cable_los_ground_angle_rad = cable_geometry["ground_angle_rad"]

    cable_model = case["cable"]["model"]
    fem_positions = None
    fem_velocities = None
    fem_result = None
    if cable_model == "fem":
        fem_positions, fem_velocities = unpack_fem_state(
            case["cable"], state[len(BASE_STATE_NAMES) :]
        )

        first_segment = fem_positions[0] - np.array(
            [cable_geometry["x_hook_m"], cable_geometry["z_hook_m"]]
        )
        first_segment_length = float(np.linalg.norm(first_segment))
        if first_segment_length <= 0.0:
            raise RuntimeError("FEM cable hook-side element length became non-positive.")
        first_segment_unit = first_segment / first_segment_length
        cable_geometry["ground_angle_rad"] = math.atan2(
            abs(first_segment_unit[1]), abs(first_segment_unit[0])
        )

    # 3. Cable tension and winch dynamics.
    winch = case["winch"]
    target_tension = math.nan
    throttle_command = math.nan
    full_throttle_torque = math.nan
    launch_power_ramp_progress = math.nan
    club_speed_reduction_active = False
    terminal_fade_progress_raw = 0.0
    terminal_fade_progress = float(pilot_state.get("terminal_fade_progress", 0.0))
    drive_torque = 0.0
    winch_integral_dot = 0.0
    drum_angle_dot = 0.0
    drum_rate_dot = 0.0
    winch_tension = 0.0
    fem_boundary = None

    if connection_state in {RELEASED, BROKEN}:
        tension = 0.0
        if winch["model"] == "dynamic":
            full_throttle_torque = full_throttle_drive_torque(winch, drum_rate)
            if winch["control"]["mode"] == "club_throttle":
                throttle_command = 0.0
            drum_angle_dot, drum_rate_dot = drum_derivatives(case, 0.0, 0.0, drum_rate)
    elif winch["model"] == "prescribed_tension":
        target_tension = scheduled_tension(
            winch["tension_schedule"],
            case["launch_setup"]["pretension_n"],
            time_s,
            cable_geometry["ground_angle_rad"],
            pilot_state.get("liftoff_time_s"),
            pilot_state.get("rotation_start_time_s"),
        )
        winch_tension = target_tension
        tension = target_tension
    else:
        if cable_model == "fem":
            fem_boundary = _fem_winch_boundary_state(
                case,
                fem_positions,
                fem_velocities,
                drum_angle,
                drum_rate,
                cable_runtime,
            )
            winch_tension = fem_boundary["tension_n"]
            tension = winch_tension
        else:
            tension = secant_elastic_tension(
                case["cable"],
                case["winch"],
                cable_geometry["length_m"],
                cable_geometry["length_rate_m_s"],
                drum_angle,
                drum_rate,
                float(cable_runtime.get("initial_reel_in_m", 0.0)),
            )
            winch_tension = tension

        control = winch["control"]
        full_throttle_torque = full_throttle_drive_torque(winch, drum_rate)

        if control["mode"] == "tension_control":
            target_tension = scheduled_tension(
                control["target_tension_schedule"],
                case["launch_setup"]["pretension_n"],
                time_s,
                cable_geometry["ground_angle_rad"],
                pilot_state.get("liftoff_time_s"),
                pilot_state.get("rotation_start_time_s"),
            )
            drive_torque, winch_integral_dot, _ = tension_control_drive_torque(
                winch,
                control,
                target_tension,
                tension,
                winch_integral,
                full_throttle_torque,
            )
        else:  # club_throttle
            (
                throttle_command,
                launch_power_ramp_progress,
                club_speed_reduction_active,
                terminal_fade_progress_raw,
                terminal_fade_progress,
            ) = club_winch_throttle(
                winch,
                case["launch_setup"]["pretension_n"],
                time_s,
                cable_los_ground_angle_rad,
                pilot_state,
            )
            drive_torque = throttle_command * full_throttle_torque

        drum_angle_dot, drum_rate_dot = drum_derivatives(
            case, winch_tension, drive_torque, drum_rate
        )

    if cable_model == "fem":
        active_node_count = fem_positions.shape[0] + 1
        removed_node_count = int(cable_runtime.get("removed_node_count", 0))
        reference_length = fem_reference_length(case["cable"])

        if connection_state in {RELEASED, BROKEN}:
            force_x_ground = 0.0
            force_z_ground = 0.0
            tension = 0.0
            winch_tension = 0.0
            fem_position_dot = np.zeros_like(fem_positions)
            fem_velocity_dot = np.zeros_like(fem_velocities)
            cable_length = cable_geometry["length_m"]
            cable_length_rate = cable_geometry["length_rate_m_s"]
            fem_max_element_tension = 0.0
            fem_max_strain = math.nan
            fem_min_segment_length = math.nan
            fem_hook_resultant_force = 0.0
            fem_winch_reference_length_m = math.nan
            fem_winch_span_length_m = math.nan
            fem_free_reference_length_m = math.nan
            fem_free_cable_mass_kg = math.nan
            fem_cable_weight_n = 0.0
            fem_cable_drag_n = 0.0
            fem_cable_drag_power_w = 0.0
            fem_ground_contact_node_count = 0
            fem_max_ground_penetration_m = 0.0
            fem_ground_normal_force_n = 0.0
            fem_ground_friction_force_n = 0.0
            fem_ground_friction_power_w = 0.0
        else:
            if winch["model"] == "dynamic":
                if fem_boundary is None:
                    fem_boundary = _fem_winch_boundary_state(
                        case,
                        fem_positions,
                        fem_velocities,
                        drum_angle,
                        drum_rate,
                        cable_runtime,
                    )
                final_reference_length = fem_boundary["reference_length_m"]
                final_reference_rate = fem_boundary["reference_length_rate_m_s"]
                winch_tension = fem_boundary["tension_n"]
            else:
                final_reference_length = reference_length
                final_reference_rate = 0.0

            hook_position = np.array(
                [cable_geometry["x_hook_m"], cable_geometry["z_hook_m"]]
            )
            hook_velocity = np.array(
                [
                    cable_geometry["x_hook_dot_m_s"],
                    cable_geometry["z_hook_dot_m_s"],
                ]
            )
            winch_position = np.array([winch["x_m"], winch.get("z_m", 0.0)])

            segment_density = None
            segment_wind = None
            if case["cable"]["fem"]["include_aerodynamic_drag"]:
                cable_points = np.vstack((hook_position, fem_positions, winch_position))
                segment_centers = 0.5 * (cable_points[:-1] + cable_points[1:])
                (
                    segment_density,
                    segment_wind_x,
                    segment_wind_z,
                ) = atmosphere(case, -segment_centers[:, 1])
                segment_wind = np.column_stack((segment_wind_x, segment_wind_z))

            fem_result = fem_cable_dynamics(
                case["cable"],
                hook_position,
                hook_velocity,
                fem_positions,
                fem_velocities,
                winch_position,
                winch_tension,
                final_reference_length_m=final_reference_length,
                final_reference_length_rate_m_s=final_reference_rate,
                segment_air_density_kg_m3=segment_density,
                segment_wind_ground_m_s=segment_wind,
                gravity_m_s2=G,
                ground=case["ground"],
            )
            force_x_ground, force_z_ground = fem_result["hook_force_ground_n"]
            tension = fem_result["hook_tension_n"]
            local_unit = fem_result["hook_segment_unit_ground"]
            cable_geometry["ground_angle_rad"] = math.atan2(
                abs(local_unit[1]), abs(local_unit[0])
            )
            fem_position_dot = fem_velocities
            fem_velocity_dot = fem_result["free_node_acceleration_m_s2"]
            cable_length = fem_result["current_length_m"]
            cable_length_rate = fem_result["current_length_rate_m_s"]
            fem_max_element_tension = fem_result["max_element_tension_n"]
            fem_max_strain = fem_result["max_strain"]
            fem_min_segment_length = fem_result["min_segment_length_m"]
            fem_hook_resultant_force = fem_result["hook_resultant_force_n"]
            fem_winch_reference_length_m = fem_result["winch_reference_length_m"]
            fem_winch_span_length_m = fem_result["winch_span_length_m"]
            fem_free_reference_length_m = fem_result["free_reference_length_m"]
            fem_free_cable_mass_kg = fem_result["free_cable_mass_kg"]
            fem_cable_weight_n = fem_result["cable_weight_n"]
            fem_cable_drag_n = fem_result["cable_drag_n"]
            fem_cable_drag_power_w = fem_result["cable_drag_power_w"]
            fem_ground_contact_node_count = fem_result["ground_contact_node_count"]
            fem_max_ground_penetration_m = fem_result["max_ground_penetration_m"]
            fem_ground_normal_force_n = fem_result["ground_normal_force_n"]
            fem_ground_friction_force_n = fem_result["ground_friction_force_n"]
            fem_ground_friction_power_w = fem_result["ground_friction_power_w"]
    else:
        force_x_ground = tension * cable_geometry["unit_x_ground"]
        force_z_ground = tension * cable_geometry["unit_z_ground"]
        cable_length = cable_geometry["length_m"]
        cable_length_rate = cable_geometry["length_rate_m_s"]
        fem_max_element_tension = math.nan
        fem_max_strain = math.nan
        fem_min_segment_length = math.nan
        fem_hook_resultant_force = math.nan
        active_node_count = 0
        removed_node_count = 0
        reference_length = math.nan
        fem_winch_reference_length_m = math.nan
        fem_winch_span_length_m = math.nan
        fem_free_reference_length_m = math.nan
        fem_free_cable_mass_kg = math.nan
        fem_cable_weight_n = math.nan
        fem_cable_drag_n = math.nan
        fem_cable_drag_power_w = math.nan
        fem_ground_contact_node_count = 0
        fem_max_ground_penetration_m = math.nan
        fem_ground_normal_force_n = math.nan
        fem_ground_friction_force_n = math.nan
        fem_ground_friction_power_w = math.nan

    # 4. Wheel/ground reaction. Club ground-roll control uses the actual
    # tailwheel normal load; the contacts themselves remain unconditional.
    ground_result = wheel_ground_forces_and_moment(
        case, x, z, theta, x_dot, z_dot, pitch_rate
    )
    ground_force_x, ground_force_z = ground_result["force_ground_n"]
    ground_x_body = math.cos(theta) * ground_force_x - math.sin(theta) * ground_force_z
    ground_z_body = math.sin(theta) * ground_force_x + math.cos(theta) * ground_force_z

    # 5. Pilot command. Standard and Club SOPs are deliberately separate because
    # their control objectives differ, not merely their parameter values.
    if case["pilot"]["sop"] == "standard":
        pilot = standard_pilot_control(
            case,
            time_s,
            air,
            theta,
            pitch_rate,
            agl,
            cable_geometry["ground_angle_rad"],
            pilot_state,
        )
    else:
        pilot = club_pilot_control(
            case,
            time_s,
            air,
            theta,
            pitch_rate,
            agl,
            tension,
            ground_result["contacts"]["tailwheel"]["normal_force_n"],
            pilot_state,
        )
    elevator = pilot["elevator_rad"]

    # 6. Aerodynamic forces and pitching moment.
    aero = aerodynamic_forces_and_moment(case, air, pitch_rate, elevator, rho)

    # 7. Cable force and hook moment.
    cable_x_body = math.cos(theta) * force_x_ground - math.sin(theta) * force_z_ground
    cable_z_body = math.sin(theta) * force_x_ground + math.cos(theta) * force_z_ground
    hook_moment = (
        case["hook"]["z_from_cg_m"] * cable_x_body
        - case["hook"]["x_from_cg_m"] * cable_z_body
    )

    # 8. Total forces, total moment and state derivatives.
    total_x = aero["x_force_n"] + cable_x_body + ground_x_body
    total_z = aero["z_force_n"] + cable_z_body + ground_z_body
    total_moment = aero["pitch_moment_nm"] + hook_moment + ground_result["pitch_moment_nm"]
    mass = case["aircraft"]["mass_kg"]
    iyy = case["aircraft"]["iyy_kg_m2"]

    base_derivative = np.array(
        [
            x_dot,
            z_dot,
            total_x / mass - pitch_rate * w - G * math.sin(theta),
            total_z / mass + pitch_rate * u + G * math.cos(theta),
            pitch_rate,
            total_moment / iyy,
            drum_angle_dot,
            drum_rate_dot,
            winch_integral_dot,
        ],
        dtype=float,
    )
    if cable_model == "fem":
        derivative = np.concatenate(
            (base_derivative, fem_position_dot.ravel(), fem_velocity_dot.ravel())
        )
    else:
        derivative = base_derivative

    if not collect_diagnostics:
        return derivative, None

    n_aero = aero["lift_n"] / (mass * G)
    stall_speed = case["limits"]["stall_speed_1g_eas_m_s"] * math.sqrt(max(n_aero, 0.0))
    stall_margin = air["eas_m_s"] - stall_speed
    max_winch_speed = case["limits"]["max_winch_eas_m_s"]
    max_winch_speed_margin = max_winch_speed - air["eas_m_s"]
    if case["pilot"]["sop"] == "club":
        early_speed_allowance = agl < case["pilot"]["hold_down"]["until_agl_m"]
    else:
        early_speed_allowance = pilot["phase"] in {"GROUND_ROLL", "INITIAL_CLIMB", "PROGRESSIVE_ROTATION"}
    active_speed_limit = (
        max_winch_speed * case["limits"]["initial_phase_max_winch_speed_factor"]
        if early_speed_allowance
        else max_winch_speed
    )
    active_speed_stop_margin = active_speed_limit - air["eas_m_s"]
    weak_link_margin = case["limits"]["weak_link_break_load_n"] - tension
    release_load_margin = case["hook"]["max_cable_load_n"] - tension
    cable_power = (
        force_x_ground * cable_geometry["x_hook_dot_m_s"]
        + force_z_ground * cable_geometry["z_hook_dot_m_s"]
    )

    if case["winch"]["model"] == "dynamic":
        drum_energy = 0.5 * case["winch"]["equivalent_inertia_kg_m2"] * drum_rate**2
        winch_balance_rhs = (
            drive_torque * drum_rate
            - winch_tension * case["winch"]["drum_radius_m"] * drum_rate
            - case["winch"]["loss_torque_nm"] * drum_rate
        )
    else:
        drum_energy = math.nan
        winch_balance_rhs = math.nan

    diagnostics = {
        "time_s": time_s,
        "x_m": x,
        "relative_altitude_m": altitude,
        "agl_m": agl,
        "u_m_s": u,
        "w_m_s": w,
        "ground_speed_x_m_s": air["ground_speed_x_m_s"],
        "climb_rate_m_s": air["climb_rate_m_s"],
        "airspeed_m_s": air["airspeed_m_s"],
        "eas_m_s": air["eas_m_s"],
        "gamma_air_deg": math.degrees(air["gamma_air_rad"]),
        "gamma_ground_deg": math.degrees(air["gamma_ground_rad"]),
        "alpha_deg": math.degrees(air["alpha_rad"]),
        "theta_deg": math.degrees(theta),
        "pitch_rate_deg_s": math.degrees(pitch_rate),
        "elevator_deg": math.degrees(elevator),
        "pilot_speed_error": pilot["speed_error"],
        "pilot_phase": pilot["phase"],
        "pilot_rotation_progress": pilot["rotation_progress"],
        "pilot_release_preparation_progress": pilot["release_preparation_progress"],
        "pilot_initial_climb_mode": (
            "takeoff_attitude" if case["pilot"]["sop"] == "standard" else "taut_cable_hold_down"
        ),
        "pilot_initial_climb_tension_deficit": pilot["initial_climb_tension_deficit"],
        "pilot_pitch_target_deg": pilot["pitch_target_deg"],
        "pilot_speed_protection_pitch_increment_deg": pilot["speed_protection_pitch_increment_deg"],
        "pilot_power_reduction_requested": pilot["power_reduction_requested"],
        "pilot_release_pitch_reference_deg": (
            case["pilot"]["release"]["pitch_reference_deg"]
            if case["pilot"]["sop"] == "standard"
            else math.nan
        ),
        "pilot_release_pitch_error_deg": (
            math.degrees(theta) - case["pilot"]["release"]["pitch_reference_deg"]
            if case["pilot"]["sop"] == "standard"
            else math.nan
        ),
        "pilot_eas_rate_m_s2": pilot_state.get("eas_rate_m_s2", 0.0),
        "CL": aero["CL"],
        "CD": aero["CD"],
        "Cm": aero["Cm"],
        "dynamic_pressure_pa": aero["dynamic_pressure_pa"],
        "lift_n": aero["lift_n"],
        "drag_n": aero["drag_n"],
        "aero_pitch_moment_nm": aero["pitch_moment_nm"],
        "aero_load_factor": n_aero,
        "aircraft_ground_contact_count": ground_result["contact_count"],
        "aircraft_max_ground_penetration_m": ground_result["max_penetration_m"],
        "aircraft_ground_signed_depth_m": ground_result["max_signed_contact_depth_m"],
        "aircraft_ground_normal_force_n": ground_result["normal_force_n"],
        "aircraft_rolling_resistance_n": ground_result["rolling_force_n"],
        "aircraft_ground_dissipation_w": ground_result["rolling_power_w"],
        "nosewheel_normal_force_n": ground_result["contacts"]["nosewheel"]["normal_force_n"],
        "mainwheel_normal_force_n": ground_result["contacts"]["mainwheel"]["normal_force_n"],
        "tailwheel_normal_force_n": ground_result["contacts"]["tailwheel"]["normal_force_n"],
        "hook_x_m": cable_geometry["x_hook_m"],
        "hook_z_m": cable_geometry["z_hook_m"],
        "winch_x_m": case["winch"]["x_m"],
        "winch_z_m": case["winch"].get("z_m", 0.0),
        "cable_model": cable_model,
        "cable_length_m": cable_length,
        "cable_length_rate_m_s": cable_length_rate,
        "cable_ground_angle_deg": math.degrees(cable_geometry["ground_angle_rad"]),
        "cable_los_ground_angle_deg": math.degrees(cable_los_ground_angle_rad),
        "target_tension_n": target_tension,
        "cable_tension_n": tension,
        "hook_tension_n": tension,
        "winch_tension_n": winch_tension if cable_model == "fem" else tension,
        "fem_hook_resultant_force_n": fem_hook_resultant_force,
        "fem_max_element_tension_n": fem_max_element_tension,
        "fem_max_strain": fem_max_strain,
        "fem_min_segment_length_m": fem_min_segment_length,
        "fem_active_node_count": active_node_count if cable_model == "fem" else math.nan,
        "fem_removed_node_count": removed_node_count if cable_model == "fem" else math.nan,
        "fem_reference_element_length_m": reference_length,
        "fem_winch_reference_length_m": fem_winch_reference_length_m,
        "fem_winch_span_length_m": fem_winch_span_length_m,
        "fem_free_reference_length_m": fem_free_reference_length_m,
        "fem_free_cable_mass_kg": fem_free_cable_mass_kg,
        "fem_cable_weight_n": fem_cable_weight_n,
        "fem_cable_drag_n": fem_cable_drag_n,
        "fem_cable_drag_power_w": fem_cable_drag_power_w,
        "fem_ground_contact_node_count": fem_ground_contact_node_count,
        "fem_max_ground_penetration_m": fem_max_ground_penetration_m,
        "fem_ground_normal_force_n": fem_ground_normal_force_n,
        "fem_ground_friction_force_n": fem_ground_friction_force_n,
        "fem_ground_friction_power_w": fem_ground_friction_power_w,
        "winch_control_mode": (
            winch["control"]["mode"] if winch["model"] == "dynamic" else "prescribed_tension"
        ),
        "winch_throttle_command": throttle_command,
        "winch_launch_power_ramp_progress": launch_power_ramp_progress,
        "club_speed_reduction_active": club_speed_reduction_active,
        "winch_terminal_fade_progress_raw": terminal_fade_progress_raw,
        "winch_terminal_fade_progress": terminal_fade_progress,
        "winch_full_throttle_torque_nm": full_throttle_torque,
        "hook_moment_nm": hook_moment,
        "cable_power_to_aircraft_w": cable_power,
        "weak_link_margin_n": weak_link_margin,
        "release_load_margin_n": release_load_margin,
        "max_winch_speed_margin_m_s": max_winch_speed_margin,
        "active_speed_stop_margin_m_s": active_speed_stop_margin,
        "stall_margin_m_s": stall_margin,
        "drum_angle_rad": drum_angle if case["winch"]["model"] == "dynamic" else math.nan,
        "drum_rate_rad_s": drum_rate if case["winch"]["model"] == "dynamic" else math.nan,
        "drive_torque_nm": drive_torque if case["winch"]["model"] == "dynamic" else math.nan,
        "drum_energy_j": drum_energy,
        "winch_energy_balance_rhs_w": winch_balance_rhs,
        "connection_state": connection_state,
        "density_kg_m3": rho,
        "wind_x_m_s": wind_x,
        "wind_z_m_s": wind_z,
    }

    if cable_model == "fem":
        initial_free_node_count = case["cable"]["fem"]["node_count"] - 1
        active_free_node_count = len(fem_positions)
        for index in range(initial_free_node_count):
            node_number = index + 2
            if index < active_free_node_count:
                diagnostics[f"cable_node_{node_number}_x_m"] = float(fem_positions[index, 0])
                diagnostics[f"cable_node_{node_number}_z_m"] = float(fem_positions[index, 1])
            else:
                diagnostics[f"cable_node_{node_number}_x_m"] = math.nan
                diagnostics[f"cable_node_{node_number}_z_m"] = math.nan

    return derivative, diagnostics

def rk4_step(time_s, state, dt_s, case, connection_state, pilot_state, cable_runtime=None):
    """Advance one fixed RK4 step with discrete cable and pilot phase state held fixed."""
    if cable_runtime is None:
        cable_runtime = {}
        if case["cable"]["model"] == "fem":
            cable_runtime = {
                "removed_node_count": 0,
                "initial_node_count": case["cable"]["fem"]["node_count"],
            }
    k1, _ = evaluate_dynamics(
        time_s,
        state,
        case,
        connection_state,
        pilot_state,
        cable_runtime,
        collect_diagnostics=False,
    )
    k2, _ = evaluate_dynamics(
        time_s + 0.5 * dt_s,
        state + 0.5 * dt_s * k1,
        case,
        connection_state,
        pilot_state,
        cable_runtime,
        collect_diagnostics=False,
    )
    k3, _ = evaluate_dynamics(
        time_s + 0.5 * dt_s,
        state + 0.5 * dt_s * k2,
        case,
        connection_state,
        pilot_state,
        cable_runtime,
        collect_diagnostics=False,
    )
    k4, _ = evaluate_dynamics(
        time_s + dt_s,
        state + dt_s * k3,
        case,
        connection_state,
        pilot_state,
        cable_runtime,
        collect_diagnostics=False,
    )
    return state + dt_s * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0

def _apply_fem_node_removal(
    time_s,
    state,
    case,
    connection_state,
    pilot_state,
    cable_runtime,
    diagnostics,
):
    """Apply discrete winch-side FEM remeshing after an accepted integration step.

    During each RK4 step the state dimension is fixed.  At the accepted step end,
    the winch-side mass point is removed when the physical final span falls below
    k*L_ref.  All remaining mass-point positions and velocities are preserved.
    """
    if case["cable"]["model"] != "fem":
        return state, diagnostics, []
    if case["winch"]["model"] != "dynamic":
        return state, diagnostics, []
    if connection_state not in {TAUT, SLACK}:
        return state, diagnostics, []

    reeling = case["cable"]["fem"]["reeling"]
    if not reeling["enabled"]:
        return state, diagnostics, []

    reference_length = fem_reference_length(case["cable"])
    threshold = reeling["node_removal_length_factor"] * reference_length
    minimum_node_count = reeling["minimum_node_count"]
    removal_events = []

    while (
        diagnostics["fem_active_node_count"] > minimum_node_count
        and diagnostics["fem_winch_span_length_m"] < threshold
    ):
        before_node_count = int(diagnostics["fem_active_node_count"])
        before_span = float(diagnostics["fem_winch_span_length_m"])
        before_tension = float(diagnostics["winch_tension_n"])
        before_cable_length = float(diagnostics["cable_length_m"])
        before_reference_length = float(diagnostics["fem_free_reference_length_m"])
        before_cable_mass = float(diagnostics["fem_free_cable_mass_kg"])
        positions, velocities = unpack_fem_state(
            case["cable"], state[len(BASE_STATE_NAMES) :]
        )
        removed_position = positions[-1].copy()
        removed_velocity = velocities[-1].copy()

        positions = positions[:-1]
        velocities = velocities[:-1]
        state = np.concatenate(
            (
                state[: len(BASE_STATE_NAMES)],
                pack_fem_state(positions, velocities),
            )
        )
        cable_runtime["removed_node_count"] = int(
            cable_runtime.get("removed_node_count", 0)
        ) + 1

        _, diagnostics = evaluate_dynamics(
            time_s,
            state,
            case,
            connection_state,
            pilot_state,
            cable_runtime,
        )
        removal_events.append(
            {
                "time": time_s,
                "code": "FEM_NODE_REMOVED",
                "severity": "info",
                "message": "Winch-side FEM mass point removed during reel-in remeshing.",
                "state_snapshot": dict(zip(state_names(case, state), state.tolist())),
                "values": {
                    "node_count_before": before_node_count,
                    "node_count_after": int(diagnostics["fem_active_node_count"]),
                    "removed_node_count": int(cable_runtime["removed_node_count"]),
                    "node_removal_threshold_m": threshold,
                    "winch_span_before_m": before_span,
                    "winch_span_after_m": float(diagnostics["fem_winch_span_length_m"]),
                    "winch_tension_before_n": before_tension,
                    "winch_tension_after_n": float(diagnostics["winch_tension_n"]),
                    "cable_length_before_m": before_cable_length,
                    "cable_length_after_m": float(diagnostics["cable_length_m"]),
                    "free_reference_length_before_m": before_reference_length,
                    "free_reference_length_after_m": float(diagnostics["fem_free_reference_length_m"]),
                    "free_cable_mass_before_kg": before_cable_mass,
                    "free_cable_mass_after_kg": float(diagnostics["fem_free_cable_mass_kg"]),
                    "removed_node_position_m": removed_position.tolist(),
                    "removed_node_velocity_m_s": removed_velocity.tolist(),
                },
            }
        )

    return state, diagnostics, removal_events

def _crossing_fraction(previous, current, threshold, direction):
    tolerance = 1.0e-7 * max(1.0, abs(threshold))
    if direction == "up":
        if previous < threshold - tolerance and current >= threshold:
            return (threshold - previous) / (current - previous)
    elif direction == "down":
        if previous > threshold + tolerance and current <= threshold:
            return (previous - threshold) / (previous - current)
    else:
        raise ValueError("direction must be 'up' or 'down'.")
    return None

def detect_step_event(previous, current, case, connection_state, pilot_state):
    """Return the highest-priority event crossed during one candidate RK4 step."""
    if connection_state not in {TAUT, SLACK}:
        return None

    candidates = []

    if pilot_state.get("liftoff_time_s") is None:
        fraction = _crossing_fraction(
            previous["aircraft_ground_signed_depth_m"],
            current["aircraft_ground_signed_depth_m"],
            0.0,
            "down",
        )
        if fraction is not None:
            candidates.append((1, fraction, "LIFTOFF", "transition"))

    fraction = _crossing_fraction(
        previous["cable_tension_n"], current["cable_tension_n"], case["hook"]["max_cable_load_n"], "up"
    )
    if fraction is not None:
        candidates.append((2, fraction, "G72_G73_LOAD_EXCEEDED", "stop"))

    fraction = _crossing_fraction(
        previous["cable_tension_n"], current["cable_tension_n"], case["limits"]["weak_link_break_load_n"], "up"
    )
    if fraction is not None:
        candidates.append((3, fraction, "WEAK_LINK_BREAK", "transition"))

    limits = case["limits"]
    vw = limits["max_winch_eas_m_s"]
    early_stop = vw * limits["initial_phase_max_winch_speed_factor"]

    if case["pilot"]["sop"] == "club":
        boundary = case["pilot"]["hold_down"]["until_agl_m"]
        altitude_fraction = _crossing_fraction(previous["agl_m"], current["agl_m"], boundary, "up")

        vw_fraction = _crossing_fraction(previous["eas_m_s"], current["eas_m_s"], vw, "up")
        if vw_fraction is not None:
            before_boundary = (
                previous["agl_m"] < boundary
                and (altitude_fraction is None or vw_fraction < altitude_fraction)
            )
            severity = "warning" if before_boundary else (
                "stop" if limits.get("stop_on_overspeed", True) else "warning"
            )
            candidates.append((4, vw_fraction, "MAX_WINCH_SPEED_EXCEEDED", severity))

        early_stop_fraction = _crossing_fraction(
            previous["eas_m_s"], current["eas_m_s"], early_stop, "up"
        )
        if early_stop_fraction is not None:
            before_boundary = (
                previous["agl_m"] < boundary
                and (altitude_fraction is None or early_stop_fraction < altitude_fraction)
            )
            if before_boundary:
                severity = "stop" if limits.get("stop_on_overspeed", True) else "warning"
                candidates.append(
                    (4, early_stop_fraction, "INITIAL_PHASE_OVERSPEED_STOP_LIMIT_EXCEEDED", severity)
                )

        # The Club SOP permits Vw exceedance only below 200 ft.  If the aircraft
        # reaches 200 ft already above Vw, the hard constraint starts exactly at
        # the altitude crossing rather than at the next integration sample.
        if altitude_fraction is not None:
            eas_at_boundary = previous["eas_m_s"] + altitude_fraction * (
                current["eas_m_s"] - previous["eas_m_s"]
            )
            if eas_at_boundary > vw:
                severity = "stop" if limits.get("stop_on_overspeed", True) else "warning"
                candidates.append((4, altitude_fraction, "MAX_WINCH_SPEED_EXCEEDED", severity))
    else:
        early_phases = {"GROUND_ROLL", "INITIAL_CLIMB", "PROGRESSIVE_ROTATION"}
        previous_early = previous["pilot_phase"] in early_phases
        current_early = current["pilot_phase"] in early_phases
        if previous_early and current_early:
            fraction = _crossing_fraction(previous["eas_m_s"], current["eas_m_s"], vw, "up")
            if fraction is not None:
                candidates.append((4, fraction, "MAX_WINCH_SPEED_EXCEEDED", "warning"))
            fraction = _crossing_fraction(previous["eas_m_s"], current["eas_m_s"], early_stop, "up")
            if fraction is not None:
                severity = "stop" if limits.get("stop_on_overspeed", True) else "warning"
                candidates.append((4, fraction, "INITIAL_PHASE_OVERSPEED_STOP_LIMIT_EXCEEDED", severity))
        else:
            fraction = _crossing_fraction(previous["eas_m_s"], current["eas_m_s"], vw, "up")
            if fraction is None and previous_early and not current_early and current["eas_m_s"] > vw:
                fraction = 1.0
            if fraction is not None:
                severity = "stop" if limits.get("stop_on_overspeed", True) else "warning"
                candidates.append((4, fraction, "MAX_WINCH_SPEED_EXCEEDED", severity))

    if pilot_state.get("liftoff_time_s") is not None:
        alpha_high = _crossing_fraction(
            previous["alpha_deg"], current["alpha_deg"], limits["aero_alpha_max_deg"], "up"
        )
        alpha_low = _crossing_fraction(
            previous["alpha_deg"], current["alpha_deg"], limits["aero_alpha_min_deg"], "down"
        )
        alpha_fraction = alpha_high if alpha_high is not None else alpha_low
        if alpha_fraction is not None:
            severity = "stop" if limits.get("stop_on_aero_validity", True) else "warning"
            candidates.append((4, alpha_fraction, "AERO_MODEL_VALIDITY_EXCEEDED", severity))

        fraction = _crossing_fraction(
            previous["stall_margin_m_s"], current["stall_margin_m_s"], 0.0, "down"
        )
        if fraction is not None:
            severity = "stop" if limits.get("stop_on_stall", True) else "warning"
            candidates.append((4, fraction, "STALL_MARGIN_EXCEEDED", severity))

    if (
        case["pilot"]["sop"] == "club"
        and pilot_state.get("power_reduction_request_time_s") is None
        and current.get("pilot_power_reduction_requested", False)
    ):
        threshold = vw - case["pilot"]["full_climb"]["power_reduction_request_margin_m_s"]
        fraction = _crossing_fraction(previous["eas_m_s"], current["eas_m_s"], threshold, "up")
        if fraction is None:
            fraction = 1.0
        candidates.append((5, fraction, "WINCH_POWER_REDUCTION_REQUEST", "transition"))

    release = case["pilot"]["release"]
    if case["pilot"]["sop"] == "standard":
        ready = (
            current["pilot_phase"] == "RELEASE_PREPARATION"
            and current["pilot_release_preparation_progress"] >= 1.0
            and current["cable_tension_n"] <= release["max_tension_n"]
            and current["cable_ground_angle_deg"] >= release["min_ground_angle_deg"]
            and abs(current["pitch_rate_deg_s"]) <= release["max_abs_pitch_rate_deg_s"]
            and abs(current["pilot_release_pitch_error_deg"]) <= release["max_pitch_error_deg"]
        )
        if ready:
            fractions = []
            checks = [
                (previous["cable_tension_n"], current["cable_tension_n"], release["max_tension_n"], "down"),
                (previous["cable_ground_angle_deg"], current["cable_ground_angle_deg"], release["min_ground_angle_deg"], "up"),
                (previous["pilot_release_preparation_progress"], current["pilot_release_preparation_progress"], 1.0, "up"),
                (abs(previous["pitch_rate_deg_s"]), abs(current["pitch_rate_deg_s"]), release["max_abs_pitch_rate_deg_s"], "down"),
                (abs(previous["pilot_release_pitch_error_deg"]), abs(current["pilot_release_pitch_error_deg"]), release["max_pitch_error_deg"], "down"),
            ]
            for a, b, threshold, direction in checks:
                already = a <= threshold if direction == "down" else a >= threshold
                fractions.append(0.0 if already else (_crossing_fraction(a, b, threshold, direction) or 1.0))
            candidates.append((6, max(fractions), "PILOT_COMMAND_RELEASE", "transition"))
    else:
        ready = (
            current["pilot_phase"] == "LEVEL_OFF"
            and current["winch_terminal_fade_progress"] >= 1.0
            and current["cable_tension_n"] <= release["max_tension_n"]
            and current["cable_los_ground_angle_deg"] >= release["min_los_ground_angle_deg"]
            and abs(current["pitch_rate_deg_s"]) <= release["max_abs_pitch_rate_deg_s"]
            and abs(current["gamma_ground_deg"]) <= release["max_abs_gamma_ground_deg"]
        )
        if ready:
            candidates.append((6, 1.0, "PILOT_COMMAND_RELEASE", "transition"))

    if not candidates:
        return None

    candidates.sort(key=lambda item: (item[0], item[1]))
    priority, fraction, code, severity = candidates[0]
    return {
        "priority": priority,
        "fraction": min(max(fraction, 0.0), 1.0),
        "code": code,
        "severity": severity,
    }

def run_simulation(case):
    """Run one winch-launch case and return history, events, metadata and case snapshot."""
    state, initialization = initialize_launch_state(case)
    initial_state_dimension = len(state)
    cable_runtime = {"initial_reel_in_m": initialization.get("initial_reel_in_m", 0.0)}
    if case["cable"]["model"] == "fem":
        cable_runtime.update({
            "removed_node_count": 0,
            "initial_node_count": case["cable"]["fem"]["node_count"],
        })
    dt = case["simulation"]["dt_s"]
    max_time = case["simulation"]["max_time_s"]
    post_release_duration = case["simulation"].get("post_release_duration_s", 2.0)
    time_s = 0.0
    connection_state = TAUT
    release_time = None
    history = []
    events = []
    validation = []

    pilot_state = {
        "liftoff_time_s": None,
        "rotation_start_time_s": None,
        "rotation_start_pitch_rad": None,
        "eas_rate_m_s2": 0.0,
        "power_reduction_request_time_s": None,
        "terminal_fade_progress": 0.0,
    }

    _, diagnostics = evaluate_dynamics(
        time_s, state, case, connection_state, pilot_state, cable_runtime
    )
    history.append(diagnostics)

    while time_s < max_time:
        remaining_time = max_time - time_s
        if release_time is not None:
            remaining_time = min(
                remaining_time, release_time + post_release_duration - time_s
            )
        if remaining_time <= 1.0e-12:
            break
        step_dt = min(dt, remaining_time)
        try:
            candidate_state = rk4_step(
                time_s,
                state,
                step_dt,
                case,
                connection_state,
                pilot_state,
                cable_runtime,
            )
        except RuntimeError as error:
            events.append(
                {
                    "time": time_s,
                    "code": "NUMERICAL_FAILURE",
                    "severity": "stop",
                    "message": str(error),
                    "state_snapshot": dict(zip(state_names(case, state), state.tolist())),
                    "values": {},
                }
            )
            validation.append(f"Simulation stopped: {error}")
            break
        candidate_time = time_s + step_dt

        if not np.all(np.isfinite(candidate_state)):
            events.append(
                {
                    "time": time_s,
                    "code": "NUMERICAL_FAILURE",
                    "severity": "stop",
                    "message": "Non-finite state generated by the integrator.",
                    "state_snapshot": dict(zip(state_names(case, state), state.tolist())),
                    "values": {},
                }
            )
            validation.append("Simulation stopped because a non-finite state was generated.")
            break

        _, candidate_diagnostics = evaluate_dynamics(
            candidate_time,
            candidate_state,
            case,
            connection_state,
            pilot_state,
            cable_runtime,
        )
        candidate_eas_rate = (
            candidate_diagnostics["eas_m_s"] - diagnostics["eas_m_s"]
        ) / (candidate_time - time_s)
        candidate_diagnostics["pilot_eas_rate_m_s2"] = candidate_eas_rate

        event = detect_step_event(
            diagnostics, candidate_diagnostics, case, connection_state, pilot_state
        )

        if event is None:
            next_time = candidate_time
            next_state = candidate_state
            next_diagnostics = candidate_diagnostics
            next_eas_rate = candidate_eas_rate
        else:
            fraction = event["fraction"]
            next_time = time_s + fraction * (candidate_time - time_s)
            next_state = state + fraction * (candidate_state - state)
            _, next_diagnostics = evaluate_dynamics(
                next_time,
                next_state,
                case,
                connection_state,
                pilot_state,
                cable_runtime,
            )
            if next_time > time_s:
                next_eas_rate = (
                    next_diagnostics["eas_m_s"] - diagnostics["eas_m_s"]
                ) / (next_time - time_s)
            else:
                next_eas_rate = pilot_state["eas_rate_m_s2"]
            next_diagnostics["pilot_eas_rate_m_s2"] = next_eas_rate

        state = next_state
        time_s = next_time
        diagnostics = next_diagnostics
        pilot_state["eas_rate_m_s2"] = next_eas_rate

        if case["pilot"]["sop"] == "club":
            pilot_state["terminal_fade_progress"] = max(
                pilot_state.get("terminal_fade_progress", 0.0),
                diagnostics.get("winch_terminal_fade_progress", 0.0),
            )

        rotation = case["pilot"]["rotation"]
        if (
            pilot_state["rotation_start_time_s"] is None
            and pilot_state.get("liftoff_time_s") is not None
            and connection_state in {TAUT, SLACK}
        ):
            if case["pilot"]["sop"] == "standard":
                rotation_ready = (
                    diagnostics["eas_m_s"] >= rotation["start_eas_m_s"]
                    and next_eas_rate >= rotation["min_eas_rate_m_s2"]
                    and diagnostics["agl_m"] >= rotation["start_min_agl_m"]
                )
            else:
                rotation_ready = diagnostics["agl_m"] >= case["pilot"]["hold_down"]["until_agl_m"]
            if rotation_ready:
                pilot_state["rotation_start_time_s"] = time_s
                pilot_state["rotation_start_pitch_rad"] = float(state[4])
                _, diagnostics = evaluate_dynamics(
                    time_s, state, case, connection_state, pilot_state, cable_runtime
                )
                diagnostics["pilot_eas_rate_m_s2"] = next_eas_rate

        if event is None:
            if connection_state in {TAUT, SLACK}:
                connection_state = TAUT if diagnostics["cable_tension_n"] > 0.0 else SLACK
                diagnostics["connection_state"] = connection_state

            state, diagnostics, removal_events = _apply_fem_node_removal(
                time_s,
                state,
                case,
                connection_state,
                pilot_state,
                cable_runtime,
                diagnostics,
            )
            diagnostics["pilot_eas_rate_m_s2"] = next_eas_rate
            diagnostics["connection_state"] = connection_state
            events.extend(removal_events)
            history.append(diagnostics)
        else:
            if event["code"] == "LIFTOFF":
                pilot_state["liftoff_time_s"] = time_s
                _, diagnostics = evaluate_dynamics(
                    time_s, state, case, connection_state, pilot_state, cable_runtime
                )
                diagnostics["pilot_eas_rate_m_s2"] = next_eas_rate
            elif event["code"] == "WINCH_POWER_REDUCTION_REQUEST":
                pilot_state["power_reduction_request_time_s"] = time_s
                _, diagnostics = evaluate_dynamics(
                    time_s, state, case, connection_state, pilot_state, cable_runtime
                )
                diagnostics["pilot_eas_rate_m_s2"] = next_eas_rate

            history.append(diagnostics)
            events.append(
                {
                    "time": time_s,
                    "code": event["code"],
                    "severity": event["severity"],
                    "message": event["code"].replace("_", " ").title(),
                    "state_snapshot": dict(zip(state_names(case, state), state.tolist())),
                    "values": {
                        "cable_tension_n": diagnostics["cable_tension_n"],
                        "cable_ground_angle_deg": diagnostics["cable_ground_angle_deg"],
                        "cable_los_ground_angle_deg": diagnostics["cable_los_ground_angle_deg"],
                        "eas_m_s": diagnostics["eas_m_s"],
                        "alpha_deg": diagnostics["alpha_deg"],
                        "pitch_rate_deg_s": diagnostics["pitch_rate_deg_s"],
                        "gamma_ground_deg": diagnostics["gamma_ground_deg"],
                        "pilot_phase": diagnostics["pilot_phase"],
                        "winch_terminal_fade_progress": diagnostics["winch_terminal_fade_progress"],
                        "aircraft_ground_contact_count": diagnostics["aircraft_ground_contact_count"],
                        "aircraft_ground_signed_depth_m": diagnostics["aircraft_ground_signed_depth_m"],
                    },
                }
            )

            if event["code"] == "LIFTOFF":
                pass
            elif event["code"] == "WEAK_LINK_BREAK":
                connection_state = BROKEN
                release_time = time_s
            elif event["code"] == "PILOT_COMMAND_RELEASE":
                connection_state = RELEASED
                release_time = time_s
            elif event["severity"] == "stop":
                break

    result_history = {key: [row[key] for row in history] for key in history[0]}
    cable_model = case["cable"]["model"]
    final_node_count = (
        int(history[-1]["fem_active_node_count"]) if cable_model == "fem" else None
    )
    metadata = {
        "case_id": case["case_id"],
        "samples": len(history),
        "final_time_s": history[-1]["time_s"],
        "final_altitude_m": history[-1]["relative_altitude_m"],
        "winch_model": case["winch"]["model"],
        "cable_model": cable_model,
        "state_dimension": len(state),
        "initial_state_dimension": initial_state_dimension,
        "final_state_dimension": len(state),
        "initial_fem_node_count": (
            case["cable"]["fem"]["node_count"] if cable_model == "fem" else None
        ),
        "final_fem_node_count": final_node_count,
        "fem_nodes_removed": (
            int(cable_runtime.get("removed_node_count", 0)) if cable_model == "fem" else 0
        ),
        "fem_reeling_enabled": (
            bool(case["cable"]["fem"]["reeling"]["enabled"])
            if cable_model == "fem"
            else False
        ),
        "winch_control": (
            case["winch"]["control"]["mode"]
            if case["winch"]["model"] == "dynamic"
            else "prescribed_tension"
        ),
        "pilot_model": (
            "ground roll -> takeoff attitude -> progressive rotation -> EAS full climb -> release preparation"
            if case["pilot"]["sop"] == "standard"
            else "tailwheel hold -> taut-cable hold-down to 200 ft -> 45 deg rotation/full climb -> level off"
        ),
        "initialization": copy.deepcopy(initialization),
        "liftoff_time_s": pilot_state.get("liftoff_time_s"),
        "aerodynamic_dataset": copy.deepcopy(case["aerodynamics"]["state"]),
        "automatic_release_enabled": bool(case["hook"].get("automatic_release_enabled", False)),
        "automatic_release_note": (
            "G72/G73 source values are stored, but the angle event is disabled until the installation-axis mapping is confirmed."
            if not case["hook"].get("automatic_release_enabled", False)
            else "Automatic-release axis mapping supplied by the case."
        ),
    }

    return {
        "history": result_history,
        "events": events,
        "metadata": metadata,
        "validation": validation,
        "case_snapshot": copy.deepcopy(case),
    }

def save_result(result, output_root="outputs"):
    """Write the result formats required by the specification."""
    case_id = result["metadata"]["case_id"]
    directory = Path(output_root) / case_id
    directory.mkdir(parents=True, exist_ok=True)

    history = result["history"]
    columns = list(history)
    row_count = len(history[columns[0]])
    with open(directory / f"{case_id}_history.csv", "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=columns)
        writer.writeheader()
        for index in range(row_count):
            writer.writerow({column: history[column][index] for column in columns})

    for suffix, value in (
        ("events", result["events"]),
        ("case_snapshot", result["case_snapshot"]),
        ("validation", result["validation"]),
    ):
        with open(directory / f"{case_id}_{suffix}.json", "w", encoding="utf-8") as file:
            json.dump(value, file, ensure_ascii=False, indent=2, allow_nan=True)

    return directory
