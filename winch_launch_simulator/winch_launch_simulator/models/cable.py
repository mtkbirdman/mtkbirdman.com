import math


def hook_winch_geometry(x, z, theta, pitch_rate, x_dot, z_dot, hook, winch):
    """Return tow-hook kinematics and the straight hook-to-winch geometry."""
    x_h_body = hook["x_from_cg_m"]
    z_h_body = hook["z_from_cg_m"]
    cos_theta = math.cos(theta)
    sin_theta = math.sin(theta)

    x_hook = x + cos_theta * x_h_body + sin_theta * z_h_body
    z_hook = z - sin_theta * x_h_body + cos_theta * z_h_body

    x_hook_dot = x_dot + pitch_rate * (-sin_theta * x_h_body + cos_theta * z_h_body)
    z_hook_dot = z_dot + pitch_rate * (-cos_theta * x_h_body - sin_theta * z_h_body)

    dx = winch["x_m"] - x_hook
    dz = winch.get("z_m", 0.0) - z_hook
    length = math.hypot(dx, dz)
    if length <= 0.0:
        raise RuntimeError("Hook-to-winch distance became non-positive.")

    unit_x = dx / length
    unit_z = dz / length
    length_rate = -(unit_x * x_hook_dot + unit_z * z_hook_dot)
    ground_angle = math.atan2(abs(dz), abs(dx))

    return {
        "x_hook_m": x_hook,
        "z_hook_m": z_hook,
        "x_hook_dot_m_s": x_hook_dot,
        "z_hook_dot_m_s": z_hook_dot,
        "length_m": length,
        "length_rate_m_s": length_rate,
        "unit_x_ground": unit_x,
        "unit_z_ground": unit_z,
        "ground_angle_rad": ground_angle,
    }


def secant_elastic_tension(
    cable,
    winch,
    cable_length_m,
    cable_length_rate_m_s,
    drum_angle_rad,
    drum_rate_rad_s,
    initial_reel_in_m,
):
    """Return tension from the Secant model's equivalent spring and damper."""
    deployed_natural_length = (
        cable["unstrained_length_m"]
        - initial_reel_in_m
        - winch["drum_radius_m"] * drum_angle_rad
    )
    extension = cable_length_m - deployed_natural_length
    extension_rate = cable_length_rate_m_s + winch["drum_radius_m"] * drum_rate_rad_s
    tension = cable["stiffness_n_m"] * extension + cable["damping_n_s_m"] * extension_rate
    return max(0.0, tension)
