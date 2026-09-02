import math

from .atmosphere import RHO0


def air_data(u, w, theta, wind_x, wind_z, rho):
    """Convert ground-referenced body-axis velocity to aerodynamic quantities."""
    u_air = u - wind_x * math.cos(theta) + wind_z * math.sin(theta)
    w_air = w - wind_x * math.sin(theta) - wind_z * math.cos(theta)
    airspeed = math.hypot(u_air, w_air)
    alpha = math.atan2(w_air, u_air)
    gamma_air = theta - alpha

    ground_x = u * math.cos(theta) + w * math.sin(theta)
    climb_rate = u * math.sin(theta) - w * math.cos(theta)
    gamma_ground = math.atan2(climb_rate, ground_x)
    eas = airspeed * math.sqrt(rho / RHO0)

    return {
        "u_air_m_s": u_air,
        "w_air_m_s": w_air,
        "airspeed_m_s": airspeed,
        "eas_m_s": eas,
        "alpha_rad": alpha,
        "gamma_air_rad": gamma_air,
        "gamma_ground_rad": gamma_ground,
        "ground_speed_x_m_s": ground_x,
        "climb_rate_m_s": climb_rate,
    }


def aerodynamic_forces_and_moment(case, air, pitch_rate_rad_s, elevator_rad, rho):
    """Evaluate the canonical G103A longitudinal aerodynamic model.

    CL and Cm use the local-linear VSPAERO stability derivatives.  CD uses
    the parabolic static-polar fit CD = CD0 + K * CL**2 evaluated at the
    instantaneous CL.
    """
    aero = case["aerodynamics"]
    airspeed = max(air["airspeed_m_s"], 1.0e-6)
    alpha = air["alpha_rad"]
    alpha_ref = math.radians(aero["reference_alpha_deg"])
    chord = aero["reference_chord_m"]

    cl = (
        aero["CL_ref"]
        + aero["CL_alpha_per_rad"] * (alpha - alpha_ref)
        + aero["CL_qhat"] * pitch_rate_rad_s * chord / (2.0 * airspeed)
        + aero["CL_elevator_per_rad"] * elevator_rad
    )
    cd = aero["CD0"] + aero["K"] * cl**2
    cm = (
        aero["Cm_ref"]
        + aero["Cm_alpha_per_rad"] * (alpha - alpha_ref)
        + aero["Cm_qhat"] * pitch_rate_rad_s * chord / (2.0 * airspeed)
        + aero["Cm_elevator_per_rad"] * elevator_rad
    )

    dynamic_pressure = 0.5 * rho * airspeed**2
    lift = dynamic_pressure * aero["reference_area_m2"] * cl
    drag = dynamic_pressure * aero["reference_area_m2"] * cd
    moment = dynamic_pressure * aero["reference_area_m2"] * chord * cm

    x_force = -drag * math.cos(alpha) + lift * math.sin(alpha)
    z_force = -drag * math.sin(alpha) - lift * math.cos(alpha)

    return {
        "dynamic_pressure_pa": dynamic_pressure,
        "CL": cl,
        "CD": cd,
        "Cm": cm,
        "lift_n": lift,
        "drag_n": drag,
        "x_force_n": x_force,
        "z_force_n": z_force,
        "pitch_moment_nm": moment,
    }


def wheel_ground_forces_and_moment(case, x_m, z_m, theta_rad, x_dot_m_s, z_dot_m_s, pitch_rate_rad_s):
    """Return the summed 2-D wheel/ground reactions in ground axes.

    Wheel coordinates are body-fixed wheel-center coordinates.  The ground is the
    flat plane ``z = ground.z_m`` with +z downward.  Each tire uses a unilateral
    spring/damper normal reaction and a smooth rolling-resistance force.  Contact
    transitions are therefore consequences of the rigid-body motion rather than
    hard-coded launch phases.
    """
    ground = case["ground"]
    settings = ground["aircraft"]
    ground_z = ground["z_m"]
    stiffness = settings["normal_stiffness_n_m"]
    damping = settings["normal_damping_n_s_m"]
    rolling = settings["rolling_resistance_coefficient"]
    velocity_scale = settings["velocity_regularization_m_s"]

    cos_theta = math.cos(theta_rad)
    sin_theta = math.sin(theta_rad)
    force_x = 0.0
    force_z = 0.0
    moment = 0.0
    contact_count = 0
    max_penetration = 0.0
    normal_total = 0.0
    rolling_total = 0.0
    rolling_power = 0.0
    contact_values = {}
    max_signed_depth = -math.inf

    for contact in settings["contacts"]:
        x_body = contact["x_from_cg_m"]
        z_body = contact["z_from_cg_m"]
        radius = contact["radius_m"]

        center_x_rel = cos_theta * x_body + sin_theta * z_body
        center_z_rel = -sin_theta * x_body + cos_theta * z_body
        center_x = x_m + center_x_rel
        center_z = z_m + center_z_rel
        center_x_dot = x_dot_m_s + pitch_rate_rad_s * (
            -sin_theta * x_body + cos_theta * z_body
        )
        center_z_dot = z_dot_m_s + pitch_rate_rad_s * (
            -cos_theta * x_body - sin_theta * z_body
        )

        signed_depth = center_z + radius - ground_z
        max_signed_depth = max(max_signed_depth, signed_depth)
        penetration = max(0.0, signed_depth)
        normal = 0.0
        rolling_force = 0.0
        if penetration > 0.0:
            normal = max(0.0, stiffness * penetration + damping * center_z_dot)
            if normal > 0.0:
                contact_count += 1
                rolling_force = -rolling * normal * math.tanh(center_x_dot / velocity_scale)

        # The tire contact point lies one radius vertically below the wheel center.
        contact_x_rel = center_x_rel
        contact_z_rel = center_z_rel + radius
        force_x += rolling_force
        force_z -= normal
        moment += contact_z_rel * rolling_force - contact_x_rel * (-normal)
        normal_total += normal
        rolling_total += rolling_force
        rolling_power += rolling_force * center_x_dot
        max_penetration = max(max_penetration, penetration)
        contact_values[contact["name"]] = {
            "center_x_m": center_x,
            "center_z_m": center_z,
            "penetration_m": penetration,
            "normal_force_n": normal,
            "rolling_force_n": rolling_force,
            "center_x_velocity_m_s": center_x_dot,
            "center_z_velocity_m_s": center_z_dot,
        }

    return {
        "force_ground_n": (force_x, force_z),
        "pitch_moment_nm": moment,
        "contact_count": contact_count,
        "max_penetration_m": max_penetration,
        "normal_force_n": normal_total,
        "rolling_force_n": rolling_total,
        "rolling_power_w": rolling_power,
        "max_signed_contact_depth_m": max_signed_depth,
        "contacts": contact_values,
    }
