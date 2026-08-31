import math

import numpy as np


def full_throttle_drive_torque(winch, drum_rate_rad_s):
    """Interpolate the equivalent full-throttle drive torque at the drum shaft."""
    curve = winch["drive"]["full_throttle_torque_curve"]
    rates = [point["drum_rate_rad_s"] for point in curve]
    torques = [point["drive_torque_nm"] for point in curve]
    return float(np.interp(drum_rate_rad_s, rates, torques))


def scheduled_tension(settings, pretension_n, time_s, ground_angle_rad, liftoff_time_s, rotation_start_time_s):
    """Return the prescribed/controlled tension target for the standard launch sequence."""
    if liftoff_time_s is None:
        target = pretension_n
    elif rotation_start_time_s is None or time_s <= rotation_start_time_s:
        target = pretension_n + settings["initial_rise_rate_n_s"] * max(
            0.0, time_s - liftoff_time_s
        )
    else:
        tension_at_rotation = pretension_n + settings["initial_rise_rate_n_s"] * max(
            0.0, rotation_start_time_s - liftoff_time_s
        )
        target = tension_at_rotation + settings["rotation_rise_rate_n_s"] * (
            time_s - rotation_start_time_s
        )

    target = min(target, settings["max_tension_n"])

    start = math.radians(settings["release_prepare_angle_deg"])
    end = math.radians(settings["release_angle_deg"])
    if ground_angle_rad > start:
        fraction = min(1.0, (ground_angle_rad - start) / (end - start))
        target += fraction * (settings["release_target_n"] - target)

    return max(0.0, target)


def club_winch_throttle(winch, pretension_n, time_s, cable_los_ground_angle_rad, pilot_state):
    """Return the continuous Club launch command, speed reduction and terminal fade."""
    control = winch["control"]

    # t=0 is ALL OUT with pretension already established.  Start from the
    # equivalent command that balances that pretension at zero drum speed,
    # then ramp to the selected launch command without using LIFTOFF as a switch.
    pretension_drive_torque = (
        pretension_n * winch["drum_radius_m"] + winch["loss_torque_nm"]
    )
    zero_speed_full_torque = full_throttle_drive_torque(winch, 0.0)
    initial_throttle = pretension_drive_torque / zero_speed_full_torque
    ramp_progress = min(max(time_s / control["launch_power_ramp_duration_s"], 0.0), 1.0)
    throttle = initial_throttle + ramp_progress * (
        control["launch_throttle"] - initial_throttle
    )

    request_time = pilot_state.get("power_reduction_request_time_s")
    speed_reduction_active = False
    if request_time is not None:
        speed = control["speed_protection"]
        if time_s >= request_time + speed["operator_response_delay_s"]:
            throttle = min(throttle, speed["reduced_throttle"])
            speed_reduction_active = True

    terminal = control["terminal_fade"]
    angle_deg = math.degrees(cable_los_ground_angle_rad)
    raw_progress = min(
        max(
            (angle_deg - terminal["start_ground_angle_deg"])
            / (terminal["end_ground_angle_deg"] - terminal["start_ground_angle_deg"]),
            0.0,
        ),
        1.0,
    )
    progress = max(float(pilot_state.get("terminal_fade_progress", 0.0)), raw_progress)
    throttle += progress * (terminal["final_throttle"] - throttle)

    return (
        min(max(throttle, 0.0), 1.0),
        ramp_progress,
        speed_reduction_active,
        raw_progress,
        progress,
    )


def drum_derivatives(case, tension_n, drive_torque_nm, drum_rate_rad_s):
    """Return drum angle/rate derivatives for the dynamic winch."""
    winch = case["winch"]
    angular_acceleration = (
        drive_torque_nm
        - winch["loss_torque_nm"]
        - tension_n * winch["drum_radius_m"]
    ) / winch["equivalent_inertia_kg_m2"]
    return drum_rate_rad_s, angular_acceleration
