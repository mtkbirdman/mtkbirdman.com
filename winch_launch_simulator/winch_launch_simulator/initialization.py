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

def _ground_pose(case, cable_force_ground_n):
    """Solve the vertical/pitch pre-launch equilibrium of the three-wheel airframe."""
    x_m = case["launch_setup"].get("x_m", 0.0)
    force_x, force_z = cable_force_ground_n
    mass = case["aircraft"]["mass_kg"]
    hook = case["hook"]

    # The measured G103A geometry gives a nose/main ground attitude near 1.4 deg.
    theta_guess = math.radians(case["launch_setup"].get("ground_pitch_guess_deg", 1.4))
    main = next(
        item for item in case["ground"]["aircraft"]["contacts"] if item["name"] == "mainwheel"
    )
    z_guess = (
        main["x_from_cg_m"] * math.sin(theta_guess)
        - main["z_from_cg_m"] * math.cos(theta_guess)
        - main["radius_m"]
        + case["ground"]["z_m"]
    )

    def residual(values):
        z_m, theta_rad = values
        ground = wheel_ground_forces_and_moment(
            case, x_m, z_m, theta_rad, 0.0, 0.0, 0.0
        )
        cos_theta = math.cos(theta_rad)
        sin_theta = math.sin(theta_rad)
        hook_x_rel = cos_theta * hook["x_from_cg_m"] + sin_theta * hook["z_from_cg_m"]
        hook_z_rel = -sin_theta * hook["x_from_cg_m"] + cos_theta * hook["z_from_cg_m"]
        cable_moment = hook_z_rel * force_x - hook_x_rel * force_z
        vertical_residual = ground["force_ground_n"][1] + force_z + mass * G
        pitch_residual = ground["pitch_moment_nm"] + cable_moment
        scale = mass * G
        return (vertical_residual / scale, pitch_residual / scale)

    solution = least_squares(
        residual,
        np.array([z_guess, theta_guess]),
        bounds=(np.array([-1.5, math.radians(-10.0)]), np.array([0.2, math.radians(15.0)])),
        xtol=1.0e-12,
        ftol=1.0e-12,
        gtol=1.0e-12,
        max_nfev=500,
    )
    if not solution.success or np.linalg.norm(residual(solution.x), ord=np.inf) > 1.0e-6:
        raise RuntimeError("Pre-launch aircraft ground equilibrium did not converge.")
    return float(solution.x[0]), float(solution.x[1])

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

    # First place the aircraft using the nominal pretension along the hook/winch line.
    z_guess = -0.59
    theta_guess = math.radians(setup.get("ground_pitch_guess_deg", 1.4))
    geometry = hook_winch_geometry(x_m, z_guess, theta_guess, 0.0, 0.0, 0.0, case["hook"], winch)
    nominal_force = (
        pretension * geometry["unit_x_ground"],
        pretension * geometry["unit_z_ground"],
    )
    z_m, theta_rad = _ground_pose(case, nominal_force)

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
    }

    if case["cable"]["model"] == "secant":
        if winch["model"] == "dynamic":
            effective_natural_length = geometry["length_m"] - pretension / case["cable"]["stiffness_n_m"]
            initial_reel_in = case["cable"]["unstrained_length_m"] - effective_natural_length
            if initial_reel_in < 0.0:
                raise RuntimeError("Secant pre-launch geometry requires a cable longer than its unstrained length.")
        initialization["initial_reel_in_m"] = float(initial_reel_in)
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
