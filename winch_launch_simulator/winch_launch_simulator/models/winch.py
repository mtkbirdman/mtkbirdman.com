import math

import numpy as np


def full_throttle_drive_torque(winch, drum_rate_rad_s):
    """Interpolate the equivalent full-throttle drive torque at the drum shaft."""
    curve = winch["drive"]["full_throttle_torque_curve"]
    rates = [point["drum_rate_rad_s"] for point in curve]
    torques = [point["drive_torque_nm"] for point in curve]
    return float(np.interp(drum_rate_rad_s, rates, torques))


def _angle_fade(value, ground_angle_rad, start_deg, end_deg, final_value):
    start = math.radians(start_deg)
    end = math.radians(end_deg)
    if ground_angle_rad <= start:
        return float(value), 0.0
    fraction = min(max((ground_angle_rad - start) / (end - start), 0.0), 1.0)
    return float(value + fraction * (final_value - value)), float(fraction)


def prescribed_tension(settings, pretension_n, time_s, ground_angle_rad, liftoff_time_s):
    """Common open-loop tension excitation used by both ideal-Secant pilot cases.

    The profile is an analysis input rather than a winch-operator SOP.  It is defined
    relative to each case's LIFTOFF so the RWTH-reference and Club-pilot cases use the same prescribed forcing law.
    """
    points = settings["time_points_after_liftoff"]
    times = [float(point["time_s"]) for point in points]
    values = [float(point["tension_n"]) for point in points]
    if liftoff_time_s is None:
        elapsed = 0.0
        target = float(pretension_n)
    else:
        elapsed = max(0.0, time_s - liftoff_time_s)
        target = float(np.interp(elapsed, times, values))
    terminal = settings["terminal_fade"]
    target, progress = _angle_fade(
        target,
        ground_angle_rad,
        terminal["start_ground_angle_deg"],
        terminal["end_ground_angle_deg"],
        terminal["final_tension_n"],
    )
    if abs(values[0] - pretension_n) > 1.0e-9:
        raise ValueError("First prescribed-tension point must equal launch pretension.")
    return max(0.0, target), progress


def rwth_target_tension(settings, pretension_n, time_s, ground_angle_rad, liftoff_time_s):
    """G103A-adapted RWTH reference force-control target.

    Source-derived structure: moderate initial force, approximately linear rise after
    liftoff, force cap, and terminal reduction with cable angle.  The numerical cap is
    a G103A analysis setting rather than the ASK-21 value from the RWTH reference.
    """
    if liftoff_time_s is None:
        target = pretension_n
    else:
        target = pretension_n + settings["rise_rate_n_s"] * max(0.0, time_s - liftoff_time_s)
    target = min(target, settings["max_tension_n"])
    terminal = settings["terminal_fade"]
    target, progress = _angle_fade(
        target,
        ground_angle_rad,
        terminal["start_ground_angle_deg"],
        terminal["end_ground_angle_deg"],
        terminal["final_tension_n"],
    )
    return max(0.0, target), progress


def club_winch_throttle(winch, pretension_n, time_s, cable_los_ground_angle_rad, pilot_state):
    """Return the continuous Club launch command, speed reduction and terminal fade."""
    control = winch["control"]

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
