import math

from .atmosphere import RHO0


GROUND_ROLL = "GROUND_ROLL"
INITIAL_CLIMB = "INITIAL_CLIMB"
PROGRESSIVE_ROTATION = "PROGRESSIVE_ROTATION"
FULL_CLIMB = "FULL_CLIMB"
RELEASE_PREPARATION = "RELEASE_PREPARATION"

GROUND_ROLL_TAILWHEEL = "GROUND_ROLL_TAILWHEEL"
POST_LIFTOFF_HOLD_DOWN = "POST_LIFTOFF_HOLD_DOWN"
LEVEL_OFF = "LEVEL_OFF"


def standard_pilot_control(
    case,
    time_s,
    air,
    theta_rad,
    pitch_rate_rad_s,
    agl_m,
    cable_ground_angle_rad,
    pilot_state,
):
    """Standard longitudinal pilot model, read in launch-phase order."""
    pilot = case["pilot"]
    lower = math.radians(pilot["elevator_min_deg"])
    upper = math.radians(pilot["elevator_max_deg"])
    pitch_rate_gain = pilot["pitch_rate_gain_s"]

    initial = pilot["initial_climb"]
    initial_pitch = math.radians(initial["pitch_target_deg"])
    initial_command = (
        math.radians(initial["elevator_bias_deg"])
        + initial["pitch_kp"] * (theta_rad - initial_pitch)
        + pitch_rate_gain * pitch_rate_rad_s
    )

    full = pilot["full_climb"]
    target_eas = full["target_eas_m_s"]
    speed_error = (
        RHO0 * target_eas**2
        - air["density_kg_m3"] * air["airspeed_m_s"] ** 2
    ) / (RHO0 * target_eas**2)
    full_command = (
        math.radians(full["elevator_bias_deg"])
        + full["eas_kp"] * speed_error
        + pitch_rate_gain * pitch_rate_rad_s
    )

    rotation_start_time = pilot_state.get("rotation_start_time_s")
    rotation_progress = 0.0
    release_progress = 0.0

    if pilot_state.get("liftoff_time_s") is None:
        phase = GROUND_ROLL
        if "ground_roll_elevator_deg" in pilot:
            unsaturated = (
                math.radians(pilot["ground_roll_elevator_deg"])
                + pitch_rate_gain * pitch_rate_rad_s
            )
        else:
            unsaturated = initial_command
    elif rotation_start_time is None:
        phase = INITIAL_CLIMB
        unsaturated = initial_command
    else:
        rotation = pilot["rotation"]
        rotation_progress = min(
            max((time_s - rotation_start_time) / rotation["min_duration_s"], 0.0),
            1.0,
        )
        full_climb_ready = (
            rotation_progress >= 1.0
            and agl_m >= rotation["full_climb_min_agl_m"]
        )

        if not full_climb_ready:
            phase = PROGRESSIVE_ROTATION
            blend = 0.5 * (1.0 - math.cos(math.pi * rotation_progress))
            rotation_pitch = math.radians(rotation["pitch_target_deg"])
            rotation_pitch_command = (
                rotation["pitch_kp"]
                * (theta_rad - rotation_pitch)
                * 16.0
                * rotation_progress**2
                * (1.0 - rotation_progress) ** 2
            )
            unsaturated = (
                initial_command
                + (full_command - initial_command) * blend
                + rotation_pitch_command
            )
        else:
            release = pilot["release_preparation"]
            cable_angle_deg = math.degrees(cable_ground_angle_rad)
            if cable_angle_deg < release["start_ground_angle_deg"]:
                phase = FULL_CLIMB
                unsaturated = full_command
            else:
                phase = RELEASE_PREPARATION
                release_progress = min(
                    max(
                        (cable_angle_deg - release["start_ground_angle_deg"])
                        / (
                            release["end_ground_angle_deg"]
                            - release["start_ground_angle_deg"]
                        ),
                        0.0,
                    ),
                    1.0,
                )
                blend = 0.5 * (1.0 - math.cos(math.pi * release_progress))
                release_pitch = math.radians(release["pitch_target_deg"])
                release_command = (
                    math.radians(release["elevator_bias_deg"])
                    + release["pitch_kp"] * (theta_rad - release_pitch)
                    + pitch_rate_gain * pitch_rate_rad_s
                )
                unsaturated = full_command + (release_command - full_command) * blend

    return {
        "elevator_rad": min(max(unsaturated, lower), upper),
        "speed_error": speed_error,
        "phase": phase,
        "rotation_progress": rotation_progress,
        "release_preparation_progress": release_progress,
        "initial_climb_tension_deficit": 0.0,
        "pitch_target_deg": math.nan,
        "speed_protection_pitch_increment_deg": 0.0,
        "power_reduction_requested": False,
    }


def club_pilot_control(
    case,
    time_s,
    air,
    theta_rad,
    pitch_rate_rad_s,
    agl_m,
    cable_tension_n,
    tailwheel_normal_force_n,
    pilot_state,
):
    """Club SOP: tailwheel hold, taut-cable hold-down, 45 deg climb, then level off."""
    pilot = case["pilot"]
    lower = math.radians(pilot["elevator_min_deg"])
    upper = math.radians(pilot["elevator_max_deg"])
    pitch_rate_gain = pilot["pitch_rate_gain_s"]

    ground_roll = pilot["ground_roll"]
    tailwheel_cue = ground_roll["tailwheel_load_cue_n"]
    tailwheel_deficit = max(0.0, (tailwheel_cue - tailwheel_normal_force_n) / tailwheel_cue)
    ground_command = (
        math.radians(ground_roll["elevator_bias_deg"])
        - math.radians(ground_roll["tailwheel_load_gain_deg"]) * tailwheel_deficit
        + pitch_rate_gain * pitch_rate_rad_s
    )

    hold_down = pilot["hold_down"]
    tension_cue = hold_down["low_tension_cue_n"]
    tension_deficit = max(0.0, (tension_cue - cable_tension_n) / tension_cue)
    hold_command = (
        math.radians(hold_down["forward_elevator_deg"])
        - math.radians(hold_down["low_tension_relief_gain_deg"]) * tension_deficit
        + pitch_rate_gain * pitch_rate_rad_s
    )

    full = pilot["full_climb"]
    speed_protection = max(
        0.0,
        air["eas_m_s"]
        - (case["limits"]["max_winch_eas_m_s"] - full["speed_protection_margin_m_s"]),
    )
    speed_pitch_increment_deg = min(
        full["max_speed_protection_pitch_increment_deg"],
        full["speed_protection_pitch_gain_deg_per_m_s"] * speed_protection,
    )
    full_pitch_target_deg = full["pitch_target_deg"] + speed_pitch_increment_deg
    full_command = (
        math.radians(full["elevator_bias_deg"])
        + full["pitch_kp"] * (theta_rad - math.radians(full_pitch_target_deg))
        + pitch_rate_gain * pitch_rate_rad_s
    )

    power_request = (
        air["eas_m_s"]
        >= case["limits"]["max_winch_eas_m_s"] - full["power_reduction_request_margin_m_s"]
    )

    rotation_progress = 0.0
    level_off_progress = float(pilot_state.get("terminal_fade_progress", 0.0))
    pitch_target_deg = math.nan

    if pilot_state.get("liftoff_time_s") is None:
        phase = GROUND_ROLL_TAILWHEEL
        unsaturated = ground_command
        power_request = False
    elif agl_m < hold_down["until_agl_m"] and pilot_state.get("rotation_start_time_s") is None:
        phase = POST_LIFTOFF_HOLD_DOWN
        unsaturated = hold_command
        power_request = False
    elif level_off_progress > 0.0:
        phase = LEVEL_OFF
        level_off = pilot["level_off"]
        gamma_error = air["gamma_ground_rad"] - math.radians(level_off["target_gamma_ground_deg"])
        level_command = (
            math.radians(level_off["elevator_bias_deg"])
            + level_off["gamma_kp"] * gamma_error
            + pitch_rate_gain * pitch_rate_rad_s
        )
        blend = 0.5 * (1.0 - math.cos(math.pi * min(max(level_off_progress, 0.0), 1.0)))
        unsaturated = full_command + (level_command - full_command) * blend
        pitch_target_deg = full_pitch_target_deg
        power_request = False
    else:
        rotation_start_time = pilot_state.get("rotation_start_time_s")
        if rotation_start_time is None:
            phase = POST_LIFTOFF_HOLD_DOWN
            unsaturated = hold_command
            power_request = False
        else:
            rotation = pilot["rotation"]
            rotation_progress = min(
                max((time_s - rotation_start_time) / rotation["duration_s"], 0.0),
                1.0,
            )
            if rotation_progress < 1.0:
                phase = PROGRESSIVE_ROTATION
                blend = 0.5 * (1.0 - math.cos(math.pi * rotation_progress))
                start_pitch = pilot_state.get("rotation_start_pitch_rad", theta_rad)
                rotation_target = start_pitch + blend * (
                    math.radians(rotation["pitch_target_deg"]) - start_pitch
                )
                pitch_target_deg = math.degrees(rotation_target)
                pitch_bump = (
                    rotation["pitch_kp"]
                    * (theta_rad - math.radians(rotation["pitch_target_deg"]))
                    * 16.0
                    * rotation_progress**2
                    * (1.0 - rotation_progress) ** 2
                )
                unsaturated = hold_command + (full_command - hold_command) * blend + pitch_bump
                power_request = False
            else:
                phase = FULL_CLIMB
                unsaturated = full_command
                pitch_target_deg = full_pitch_target_deg

    return {
        "elevator_rad": min(max(unsaturated, lower), upper),
        "speed_error": math.nan,
        "phase": phase,
        "rotation_progress": rotation_progress,
        "release_preparation_progress": level_off_progress,
        "initial_climb_tension_deficit": tension_deficit,
        "pitch_target_deg": pitch_target_deg,
        "speed_protection_pitch_increment_deg": speed_pitch_increment_deg,
        "power_reduction_requested": bool(power_request and phase == FULL_CLIMB),
    }


def tension_control_drive_torque(
    winch,
    control,
    target_tension_n,
    current_tension_n,
    integral_state,
    available_drive_torque_nm,
):
    """Return drive torque requested by tension control, limited by current capability."""
    error = target_tension_n - current_tension_n
    torque = (
        winch["loss_torque_nm"]
        + target_tension_n * winch["drum_radius_m"]
        + control["kp_nm_per_n"] * error
        + control["ki_nm_per_n_s"] * integral_state
    )
    torque = min(max(torque, 0.0), available_drive_torque_nm)

    integral_dot = error
    if torque >= available_drive_torque_nm and error > 0.0:
        integral_dot = 0.0
    elif torque <= 0.0 and error < 0.0:
        integral_dot = 0.0

    return torque, integral_dot, error
