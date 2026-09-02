import math

import numpy as np
from scipy.optimize import least_squares, minimize

from .models.aircraft import wheel_ground_forces_and_moment
from .models.atmosphere import G
from .models.cable import hook_winch_geometry
from .models.fem_cable import (
    fem_cable_dynamics,
    fem_reference_length,
    fem_winch_boundary_tension,
    pack_fem_state,
)

BASE_STATE_SIZE = 9

def _contact_kinematics(contact, z_m, theta_rad, ground_z):
    """Return static 2-D wheel-center/contact geometry for one tire."""
    x_body = contact["x_from_cg_m"]
    z_body = contact["z_from_cg_m"]
    radius = contact["radius_m"]
    cos_theta = math.cos(theta_rad)
    sin_theta = math.sin(theta_rad)
    center_x_rel = cos_theta * x_body + sin_theta * z_body
    center_z_rel = -sin_theta * x_body + cos_theta * z_body
    signed_depth = z_m + center_z_rel + radius - ground_z
    return {
        "center_x_rel_m": center_x_rel,
        "center_z_rel_m": center_z_rel,
        "contact_x_rel_m": center_x_rel,
        "signed_depth_m": signed_depth,
    }


def _ground_pose(case, pretension_n):
    """Solve a mainwheel+tailwheel pre-launch static equilibrium.

    The launch starts from the practical tail-down two-point support condition:
    mainwheel and tailwheel are active unilateral contacts, while the nosewheel
    must remain clear of the ground.  The compliant tire law is retained during
    the static solve so that the t=0 state is continuous with the time-domain
    ground-contact model rather than starting from zero tire reaction.
    """
    x_m = case["launch_setup"].get("x_m", 0.0)
    mass = case["aircraft"]["mass_kg"]
    hook = case["hook"]
    winch = case["winch"]
    ground = case["ground"]
    ground_z = ground["z_m"]
    stiffness = ground["aircraft"]["normal_stiffness_n_m"]
    contacts = {item["name"]: item for item in ground["aircraft"]["contacts"]}
    main = contacts["mainwheel"]
    tail = contacts["tailwheel"]
    nose = contacts["nosewheel"]

    # First obtain the rigid two-point attitude as a robust initial guess.
    def rigid_residual(values):
        z_value, theta_value = values
        return np.array(
            [
                _contact_kinematics(main, z_value, theta_value, ground_z)["signed_depth_m"],
                _contact_kinematics(tail, z_value, theta_value, ground_z)["signed_depth_m"],
            ]
        )

    rigid_solution = least_squares(
        rigid_residual,
        np.array([-0.60, math.radians(5.5)]),
        bounds=(np.array([-1.5, math.radians(-5.0)]), np.array([0.2, math.radians(15.0)])),
        xtol=1.0e-13,
        ftol=1.0e-13,
        gtol=1.0e-13,
        max_nfev=500,
    )
    if not rigid_solution.success:
        raise RuntimeError("Rigid mainwheel/tailwheel ground attitude did not converge.")

    def residual(values):
        z_value, theta_value = values
        main_geom = _contact_kinematics(main, z_value, theta_value, ground_z)
        tail_geom = _contact_kinematics(tail, z_value, theta_value, ground_z)
        main_normal = stiffness * max(main_geom["signed_depth_m"], 0.0)
        tail_normal = stiffness * max(tail_geom["signed_depth_m"], 0.0)

        geometry = hook_winch_geometry(
            x_m, z_value, theta_value, 0.0, 0.0, 0.0, hook, winch
        )
        force_x = pretension_n * geometry["unit_x_ground"]
        force_z = pretension_n * geometry["unit_z_ground"]

        cos_theta = math.cos(theta_value)
        sin_theta = math.sin(theta_value)
        hook_x_rel = cos_theta * hook["x_from_cg_m"] + sin_theta * hook["z_from_cg_m"]
        hook_z_rel = -sin_theta * hook["x_from_cg_m"] + cos_theta * hook["z_from_cg_m"]
        cable_moment = hook_z_rel * force_x - hook_x_rel * force_z

        vertical_residual = -main_normal - tail_normal + force_z + mass * G
        pitch_residual = (
            main_geom["contact_x_rel_m"] * main_normal
            + tail_geom["contact_x_rel_m"] * tail_normal
            + cable_moment
        )
        scale = mass * G
        return np.array([vertical_residual / scale, pitch_residual / scale])

    # A slight tail penetration is expected because the tail reaction is small at
    # the canonical 2500 N pretension.  Starting from the rigid support attitude
    # selects the physically intended two-contact equilibrium rather than the
    # alternative mainwheel-only equilibrium of the unilateral contact system.
    z0, theta0 = rigid_solution.x
    equilibrium_solution = least_squares(
        residual,
        np.array([z0 + 0.010, theta0]),
        bounds=(np.array([-1.5, math.radians(3.5)]), np.array([0.2, math.radians(8.0)])),
        xtol=1.0e-13,
        ftol=1.0e-13,
        gtol=1.0e-13,
        max_nfev=1000,
    )
    if (
        not equilibrium_solution.success
        or np.linalg.norm(residual(equilibrium_solution.x), ord=np.inf) > 1.0e-7
    ):
        raise RuntimeError("Pre-launch mainwheel/tailwheel equilibrium did not converge.")

    z_m, theta_rad = map(float, equilibrium_solution.x)
    ground_result = wheel_ground_forces_and_moment(case, x_m, z_m, theta_rad, 0.0, 0.0, 0.0)
    main_result = ground_result["contacts"]["mainwheel"]
    tail_result = ground_result["contacts"]["tailwheel"]
    nose_result = ground_result["contacts"]["nosewheel"]

    if main_result["normal_force_n"] <= 0.0 or tail_result["normal_force_n"] <= 0.0:
        raise RuntimeError("Requested two-point launch attitude does not support both mainwheel and tailwheel.")
    if nose_result["penetration_m"] > 1.0e-10:
        raise RuntimeError("Requested two-point launch attitude places the nosewheel on the ground.")

    # Positive clearance means the bottom of the tire lies above the ground plane.
    nose_geom = _contact_kinematics(nose, z_m, theta_rad, ground_z)
    nose_clearance = max(0.0, -nose_geom["signed_depth_m"])
    metadata = {
        "ground_support_mode": "mainwheel_tailwheel_two_point",
        "rigid_two_point_pitch_deg": math.degrees(float(rigid_solution.x[1])),
        "mainwheel_normal_force_n": float(main_result["normal_force_n"]),
        "tailwheel_normal_force_n": float(tail_result["normal_force_n"]),
        "nosewheel_normal_force_n": float(nose_result["normal_force_n"]),
        "mainwheel_penetration_m": float(main_result["penetration_m"]),
        "tailwheel_penetration_m": float(tail_result["penetration_m"]),
        "nosewheel_clearance_m": float(nose_clearance),
    }
    return z_m, theta_rad, metadata

def _relax_fem_prelaunch(case, hook_position, winch_position, target_tension_n):
    """Find a gravity/ground equilibrium FEM shape and the reel-in giving target winch tension."""
    cable = case["cable"]
    fem = cable["fem"]
    ground = case["ground"]
    node_count = fem["node_count"]
    free_count = node_count - 1
    reference_length = fem_reference_length(cable)
    total_unstrained_length = cable["unstrained_length_m"]
    area = math.pi * fem["diameter_m"] ** 2 / 4.0
    axial_rigidity = fem["youngs_modulus_pa"] * area
    linear_density = fem["linear_density_kg_m"]

    fractions = np.arange(1, node_count, dtype=float) / node_count
    positions = hook_position + fractions[:, None] * (winch_position - hook_position)
    positions[:, 1] = ground["z_m"] - 0.5 * fem["diameter_m"]

    geometric_length = float(np.linalg.norm(hook_position - positions[0]))
    if free_count > 1:
        geometric_length += float(np.sum(np.linalg.norm(positions[1:] - positions[:-1], axis=1)))
    geometric_length += float(np.linalg.norm(winch_position - positions[-1]))
    target_strain = target_tension_n / axial_rigidity
    reel_in = max(0.0, total_unstrained_length - geometric_length / (1.0 + target_strain))
    reel_in = min(reel_in, 0.8 * reference_length)

    ground_cable = ground["cable"]
    ground_z = ground["z_m"]
    cable_radius = 0.5 * fem["diameter_m"]
    stiffness_per_length = ground_cable["normal_stiffness_n_m2"]

    nominal_x = positions[:, 0].copy()
    x_half_width = 0.45 * reference_length
    bounds = []
    for x_value in nominal_x:
        bounds.append((max(hook_position[0] + 1.0e-4, x_value - x_half_width), min(winch_position[0] - 1.0e-4, x_value + x_half_width)))
        bounds.append((ground_z - 5.0, ground_z + 0.05))

    last_result = None
    last_positions = positions.copy()

    for _ in range(8):
        final_reference = reference_length - reel_in
        if final_reference <= 0.2 * reference_length:
            raise RuntimeError("Pre-launch FEM reel-in made the final reference element too short.")

        internal_refs = np.full(node_count - 1, reference_length, dtype=float)
        span_refs = np.concatenate((internal_refs, [final_reference]))
        node_mass = np.empty(node_count, dtype=float)
        node_mass[0] = 0.5 * linear_density * reference_length
        if node_count > 2:
            node_mass[1:-1] = linear_density * reference_length
        node_mass[-1] = 0.5 * linear_density * (reference_length + final_reference)
        tributary_length = node_mass / linear_density

        def energy_and_gradient(flat):
            free_positions = np.asarray(flat, dtype=float).reshape(free_count, 2)
            points = np.vstack((hook_position, free_positions, winch_position))
            vectors = points[1:] - points[:-1]
            lengths = np.linalg.norm(vectors, axis=1)
            if np.any(lengths <= 1.0e-8):
                return 1.0e30, np.zeros_like(flat)

            extensions = np.maximum(lengths - span_refs, 0.0)
            elastic_energy = float(np.sum(0.5 * axial_rigidity / span_refs * extensions**2))
            gravity_energy = float(-np.sum(node_mass[1:] * G * free_positions[:, 1]))
            penetrations = np.maximum(0.0, free_positions[:, 1] + cable_radius - ground_z)
            ground_energy = float(
                np.sum(0.5 * stiffness_per_length * tributary_length[1:] * penetrations**2)
            )
            energy = elastic_energy + gravity_energy + ground_energy

            zero_velocities = np.zeros_like(free_positions)
            boundary = fem_winch_boundary_tension(
                cable,
                free_positions[-1],
                zero_velocities[-1],
                winch_position,
                final_reference,
                0.0,
            )
            result = fem_cable_dynamics(
                cable,
                hook_position,
                np.zeros(2),
                free_positions,
                zero_velocities,
                winch_position,
                boundary["tension_n"],
                final_reference_length_m=final_reference,
                final_reference_length_rate_m_s=0.0,
                include_weight=True,
                include_aerodynamic_drag=False,
                gravity_m_s2=G,
                ground=ground,
            )
            free_force = result["free_node_acceleration_m_s2"] * node_mass[1:, None]
            gradient = -free_force.ravel()
            return energy, gradient

        optimized = minimize(
            energy_and_gradient,
            last_positions.ravel(),
            method="L-BFGS-B",
            jac=True,
            bounds=bounds,
            options={"maxiter": 1000, "ftol": 1.0e-12, "gtol": 1.0e-5, "maxls": 40},
        )
        if not optimized.success and np.linalg.norm(optimized.jac, ord=np.inf) > 0.5:
            raise RuntimeError(f"Pre-launch FEM equilibrium did not converge: {optimized.message}")

        last_positions = optimized.x.reshape(free_count, 2)
        zero_velocities = np.zeros_like(last_positions)
        boundary = fem_winch_boundary_tension(
            cable,
            last_positions[-1],
            zero_velocities[-1],
            winch_position,
            final_reference,
            0.0,
        )
        last_result = fem_cable_dynamics(
            cable,
            hook_position,
            np.zeros(2),
            last_positions,
            zero_velocities,
            winch_position,
            boundary["tension_n"],
            final_reference_length_m=final_reference,
            final_reference_length_rate_m_s=0.0,
            include_weight=True,
            include_aerodynamic_drag=False,
            gravity_m_s2=G,
            ground=ground,
        )
        tension_error = target_tension_n - boundary["tension_n"]
        if abs(tension_error) <= 2.0:
            break
        # The whole cable has an equivalent axial stiffness close to EA/L.  This
        # update is deliberately simple; each iteration re-solves the full static shape.
        reel_in += tension_error / (axial_rigidity / total_unstrained_length)
        reel_in = min(max(reel_in, 0.0), 0.8 * reference_length)
    else:
        raise RuntimeError("Pre-launch FEM pretension solve did not converge.")

    max_force_residual = 0.0
    if last_result is not None:
        max_force_residual = float(
            np.max(np.linalg.norm(last_result["free_node_acceleration_m_s2"] * node_mass[1:, None], axis=1))
        )

    return last_positions, float(reel_in), {
        "initial_winch_tension_n": float(boundary["tension_n"]),
        "initial_hook_tension_n": float(last_result["hook_tension_n"]),
        "fem_static_max_free_node_force_n": max_force_residual,
        "fem_static_max_ground_penetration_m": float(last_result["max_ground_penetration_m"]),
    }

def initialize_launch_state(case):
    """Build the physical t=0 state and derived pre-launch setup for a ground-start case."""
    setup = case["launch_setup"]
    if setup["mode"] != "ground_start":
        raise ValueError("Only launch_setup.mode='ground_start' is supported by the canonical cases.")

    x_m = setup.get("x_m", 0.0)
    pretension = setup["pretension_n"]
    winch = case["winch"]

    # Start from the practical mainwheel+tailwheel two-point support state.
    z_m, theta_rad, ground_metadata = _ground_pose(case, pretension)

    geometry = hook_winch_geometry(x_m, z_m, theta_rad, 0.0, 0.0, 0.0, case["hook"], winch)
    base_state = np.array(
        [x_m, z_m, 0.0, 0.0, theta_rad, 0.0, 0.0, 0.0, 0.0],
        dtype=float,
    )

    initial_reel_in = 0.0
    initialization = {
        "ground_start": True,
        "initial_cg_height_m": -z_m + case["ground"]["z_m"],
        "initial_pitch_deg": math.degrees(theta_rad),
        "target_pretension_n": pretension,
        **ground_metadata,
    }

    if case["cable"]["model"] == "ideal_secant":
        initialization["initial_reel_in_m"] = 0.0
        initialization["ideal_secant"] = True
        return base_state, initialization

    hook_position = np.array([geometry["x_hook_m"], geometry["z_hook_m"]])
    winch_position = np.array([winch["x_m"], winch.get("z_m", case["ground"]["z_m"])])
    free_positions, initial_reel_in, fem_metadata = _relax_fem_prelaunch(
        case, hook_position, winch_position, pretension
    )
    initialization.update(fem_metadata)
    initialization["initial_reel_in_m"] = float(initial_reel_in)
    fem_state = pack_fem_state(free_positions, np.zeros_like(free_positions))
    return np.concatenate((base_state, fem_state)), initialization
