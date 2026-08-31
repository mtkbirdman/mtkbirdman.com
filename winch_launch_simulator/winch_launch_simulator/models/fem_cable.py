import math

import numpy as np


_EPS = 1.0e-12


def fem_reference_length(cable):
    """Return the nominal unstrained FEM element length at launch start."""
    return cable["unstrained_length_m"] / cable["fem"]["node_count"]


def fem_node_count_from_state(state):
    """Infer the active FEM mass-point count from a FEM-only state suffix."""
    state = np.asarray(state, dtype=float)
    if state.size % 4 != 0:
        raise ValueError(
            "FEM cable state suffix must contain 4 scalars per free node."
        )
    return state.size // 4 + 1


def pack_fem_state(free_node_positions, free_node_velocities):
    """Pack free-node positions and velocities into the FEM state suffix."""
    positions = np.asarray(free_node_positions, dtype=float)
    velocities = np.asarray(free_node_velocities, dtype=float)
    if positions.shape != velocities.shape or positions.ndim != 2 or positions.shape[1] != 2:
        raise ValueError("FEM node positions and velocities must both have shape (N, 2).")
    return np.concatenate((positions.ravel(), velocities.ravel()))


def unpack_fem_state(cable, state, node_count=None):
    """Return free-node positions and velocities from the FEM state suffix."""
    state = np.asarray(state, dtype=float)
    if node_count is None:
        node_count = fem_node_count_from_state(state)

    free_node_count = node_count - 1
    coordinate_count = free_node_count * 2
    expected = 2 * coordinate_count
    if state.size != expected:
        raise ValueError(
            f"FEM cable state requires {expected} values, received {state.size}."
        )

    positions = state[:coordinate_count].reshape(free_node_count, 2)
    velocities = state[coordinate_count:].reshape(free_node_count, 2)
    return positions, velocities


def fem_winch_reference_length(
    cable,
    drum_angle_rad,
    removed_node_count,
    initial_reel_in_m=0.0,
    drum_radius_m=0.0,
):
    """Return the unstrained length of the active winch-side FEM element.

    Positive drum angle is the reel-in direction. Each removed FEM mass point adds
    one nominal element length back to the active final element, preserving the
    remaining free natural cable length across a remesh.
    """
    reference_length = fem_reference_length(cable)
    reeled_length = initial_reel_in_m + drum_radius_m * drum_angle_rad
    return reference_length * (removed_node_count + 1) - reeled_length


def fem_winch_boundary_tension(
    cable,
    last_node_position,
    last_node_velocity,
    winch_position,
    final_reference_length_m,
    final_reference_length_rate_m_s,
):
    """Evaluate the elastic/damped tension in the variable-length winch element."""
    fem = cable["fem"]
    last_node_position = np.asarray(last_node_position, dtype=float)
    last_node_velocity = np.asarray(last_node_velocity, dtype=float)
    winch_position = np.asarray(winch_position, dtype=float)

    span_vector = winch_position - last_node_position
    span_length = float(np.linalg.norm(span_vector))
    if span_length <= 0.0:
        raise RuntimeError("FEM cable winch-side span length became non-positive.")
    if final_reference_length_m <= _EPS:
        raise RuntimeError("FEM cable winch-side reference length became non-positive.")

    span_unit = span_vector / span_length
    span_length_rate = -float(np.dot(span_unit, last_node_velocity))
    strain = span_length / final_reference_length_m - 1.0
    strain_rate = (
        span_length_rate / final_reference_length_m
        - span_length
        * final_reference_length_rate_m_s
        / final_reference_length_m**2
    )

    area = math.pi * fem["diameter_m"] ** 2 / 4.0
    tension = area * (
        fem["youngs_modulus_pa"] * strain
        + fem["strain_rate_stress_pa_s"] * strain_rate
    )

    return {
        "tension_n": max(float(tension), 0.0),
        "span_length_m": span_length,
        "span_length_rate_m_s": span_length_rate,
        "span_unit_ground": span_unit,
        "strain": float(strain),
        "strain_rate_s_inv": float(strain_rate),
    }


def fem_element_drag_forces(
    diameter_m,
    drag_coefficient,
    element_vectors,
    element_lengths_m,
    element_center_velocities,
    air_density_kg_m3,
    wind_ground_m_s,
):
    """Return WLO07 aerodynamic drag for all active cylindrical cable elements."""
    element_vectors = np.asarray(element_vectors, dtype=float)
    element_lengths = np.asarray(element_lengths_m, dtype=float)
    center_velocity = np.asarray(element_center_velocities, dtype=float)
    density = np.asarray(air_density_kg_m3, dtype=float)
    wind = np.asarray(wind_ground_m_s, dtype=float)

    element_unit = element_vectors / element_lengths[:, None]
    relative_velocity = center_velocity - wind
    speed = np.linalg.norm(relative_velocity, axis=1)

    relative_unit = np.zeros_like(relative_velocity)
    moving = speed > _EPS
    relative_unit[moving] = relative_velocity[moving] / speed[moving, None]

    tangential_speed_signed = np.sum(relative_velocity * element_unit, axis=1)
    tangential_velocity = tangential_speed_signed[:, None] * element_unit
    normal_velocity = relative_velocity - tangential_velocity

    tangential_speed = np.linalg.norm(tangential_velocity, axis=1)
    normal_speed = np.linalg.norm(normal_velocity, axis=1)
    tangential_unit = np.zeros_like(tangential_velocity)
    normal_unit = np.zeros_like(normal_velocity)
    tangential = tangential_speed > _EPS
    normal = normal_speed > _EPS
    tangential_unit[tangential] = (
        tangential_velocity[tangential] / tangential_speed[tangential, None]
    )
    normal_unit[normal] = normal_velocity[normal] / normal_speed[normal, None]

    cos_eta = np.clip(
        np.abs(np.sum(relative_unit * element_unit, axis=1)),
        0.0,
        1.0,
    )
    eta = np.arccos(cos_eta)
    normal_loading = (
        0.5
        - 0.1 * np.cos(eta)
        + 0.1 * np.sin(eta)
        - 0.4 * np.cos(2.0 * eta)
        - 0.011 * np.sin(2.0 * eta)
    )
    tangential_loading = 0.01 * (
        2.008
        - 0.3858 * eta
        + 1.9159 * eta**2
        - 4.16147 * eta**3
        + 3.5064 * eta**4
        - 1.187299 * eta**5
    )

    direction = (
        normal_loading[:, None] * normal_unit
        + tangential_loading[:, None] * tangential_unit
    )
    magnitude = (
        0.5
        * density
        * float(drag_coefficient)
        * float(diameter_m)
        * element_lengths
        * speed**2
    )
    return -magnitude[:, None] * direction


def _node_masses(linear_density_kg_m, internal_reference_length_m, final_reference_length_m, node_count):
    """Half-lump adjacent element masses to each active cable mass point."""
    masses = np.empty(node_count, dtype=float)
    if node_count == 1:
        masses[0] = 0.5 * linear_density_kg_m * final_reference_length_m
        return masses

    masses[0] = 0.5 * linear_density_kg_m * internal_reference_length_m
    if node_count > 2:
        masses[1:-1] = linear_density_kg_m * internal_reference_length_m
    masses[-1] = 0.5 * linear_density_kg_m * (
        internal_reference_length_m + final_reference_length_m
    )
    return masses


def fem_cable_dynamics(
    cable,
    hook_position,
    hook_velocity,
    free_node_positions,
    free_node_velocities,
    winch_position,
    winch_tension_n,
    *,
    final_reference_length_m=None,
    final_reference_length_rate_m_s=0.0,
    segment_air_density_kg_m3=None,
    segment_wind_ground_m_s=None,
    include_weight=None,
    include_aerodynamic_drag=None,
    gravity_m_s2=9.80665,
    ground=None,
):
    """Evaluate the elastic lumped-mass FEM cable model.

    The active cable consists of mass points from the tow hook toward the winch and
    a final span from the last mass point to the fixed winch location. Internal
    elements use damped Hooke elasticity. The final span tension is supplied as the
    winch boundary condition; for a dynamic winch it is normally obtained from
    ``fem_winch_boundary_tension`` using the variable final reference length.

    Cable gravity, WLO07 cylindrical aerodynamic drag, and unilateral flat-ground
    contact are optional case-level completions. Element drag is evaluated at
    element centers and half-lumped to adjacent mass points. The half-load at the
    fixed winch endpoint is a boundary reaction and is not integrated as a free
    cable mass.
    """
    fem = cable["fem"]

    hook_position = np.asarray(hook_position, dtype=float)
    hook_velocity = np.asarray(hook_velocity, dtype=float)
    free_node_positions = np.asarray(free_node_positions, dtype=float)
    free_node_velocities = np.asarray(free_node_velocities, dtype=float)
    winch_position = np.asarray(winch_position, dtype=float)

    if free_node_positions.ndim != 2 or free_node_positions.shape[1:] != (2,):
        raise ValueError("free_node_positions must have shape (N, 2).")
    if free_node_velocities.shape != free_node_positions.shape:
        raise ValueError("free_node_velocities must match free_node_positions shape.")

    node_count = free_node_positions.shape[0] + 1
    if node_count < 2:
        raise ValueError("The active FEM cable requires at least two mass points.")

    reference_length = fem_reference_length(cable)
    if final_reference_length_m is None:
        final_reference_length_m = reference_length
    final_reference_length_m = float(final_reference_length_m)
    if final_reference_length_m <= _EPS:
        raise RuntimeError("FEM cable winch-side reference length became non-positive.")

    if include_weight is None:
        include_weight = fem["include_weight"]
    if include_aerodynamic_drag is None:
        include_aerodynamic_drag = fem["include_aerodynamic_drag"]

    positions = np.vstack((hook_position, free_node_positions))
    velocities = np.vstack((hook_velocity, free_node_velocities))

    # n active mass points imply n spans: n-1 elastic internal spans plus the
    # winch-side boundary span.
    internal_vector = positions[1:] - positions[:-1]
    internal_length = np.linalg.norm(internal_vector, axis=1)
    if np.any(internal_length <= 0.0):
        raise RuntimeError("FEM cable element length became non-positive.")
    internal_unit = internal_vector / internal_length[:, None]

    relative_velocity = velocities[1:] - velocities[:-1]
    internal_length_rate = np.sum(relative_velocity * internal_unit, axis=1)
    internal_strain = internal_length / reference_length - 1.0

    area = math.pi * fem["diameter_m"] ** 2 / 4.0
    internal_tension = area * (
        fem["youngs_modulus_pa"] * internal_strain
        + fem["strain_rate_stress_pa_s"]
        * internal_length_rate
        / reference_length
    )
    # physical_completion: a flexible tow cable transmits tension, not compression.
    internal_tension = np.maximum(internal_tension, 0.0)

    last_to_winch = winch_position - positions[-1]
    last_span_length = float(np.linalg.norm(last_to_winch))
    if last_span_length <= 0.0:
        raise RuntimeError("FEM cable winch-side span length became non-positive.")
    last_span_unit = last_to_winch / last_span_length
    last_span_length_rate = -float(np.dot(last_span_unit, velocities[-1]))

    winch_tension_n = max(float(winch_tension_n), 0.0)

    # Internal axial forces.
    node_force = np.zeros_like(positions)
    internal_force = internal_tension[:, None] * internal_unit
    node_force[:-1] += internal_force
    node_force[1:] -= internal_force
    node_force[-1] += winch_tension_n * last_span_unit

    # Build all n physical spans for external-force calculations.
    span_start = positions
    span_end = np.vstack((positions[1:], winch_position))
    span_vector = span_end - span_start
    span_length = np.linalg.norm(span_vector, axis=1)
    span_velocity_start = velocities
    span_velocity_end = np.vstack((velocities[1:], np.zeros(2, dtype=float)))
    span_center_velocity = 0.5 * (span_velocity_start + span_velocity_end)

    segment_drag_force = np.zeros((node_count, 2), dtype=float)
    if include_aerodynamic_drag:
        if segment_air_density_kg_m3 is None or segment_wind_ground_m_s is None:
            raise ValueError(
                "Cable aerodynamic drag requires segment density and wind arrays."
            )
        density = np.asarray(segment_air_density_kg_m3, dtype=float)
        wind = np.asarray(segment_wind_ground_m_s, dtype=float)
        if density.shape != (node_count,):
            raise ValueError(
                f"segment_air_density_kg_m3 must have shape {(node_count,)}."
            )
        if wind.shape != (node_count, 2):
            raise ValueError(
                f"segment_wind_ground_m_s must have shape {(node_count, 2)}."
            )

        segment_drag_force = fem_element_drag_forces(
            fem["diameter_m"],
            fem["drag_coefficient"],
            span_vector,
            span_length,
            span_center_velocity,
            density,
            wind,
        )

        # Half-lump each cylinder force to its two endpoints. The winch-end half
        # is a fixed-boundary reaction and therefore not part of node_force.
        half_drag = 0.5 * segment_drag_force
        node_force += half_drag
        node_force[1:] += half_drag[:-1]

    node_mass = _node_masses(
        fem["linear_density_kg_m"],
        reference_length,
        final_reference_length_m,
        node_count,
    )
    if np.any(node_mass <= 0.0):
        raise RuntimeError("FEM cable lumped mass became non-positive.")

    if include_weight:
        node_force[:, 1] += node_mass * float(gravity_m_s2)

    ground_normal_total = 0.0
    ground_friction_total = 0.0
    ground_friction_power = 0.0
    ground_contact_count = 0
    max_ground_penetration = 0.0
    if ground is not None:
        settings = ground["cable"]
        ground_z = ground["z_m"]
        radius = 0.5 * fem["diameter_m"]
        stiffness_per_length = settings["normal_stiffness_n_m2"]
        damping_per_length = settings["normal_damping_n_s_m2"]
        friction_coefficient = settings["sliding_friction_coefficient"]
        velocity_scale = settings["velocity_regularization_m_s"]

        # node_mass / linear_density is the tributary unstrained cable length
        # represented by each lumped mass. Scaling the contact law with this
        # length keeps the ground support approximately mesh-independent.
        free_position = positions[1:]
        free_velocity = velocities[1:]
        penetration = np.maximum(0.0, free_position[:, 1] + radius - ground_z)
        tributary_length = node_mass[1:] / fem["linear_density_kg_m"]
        normal_force = np.zeros_like(penetration)
        penetrating = penetration > 0.0
        normal_force[penetrating] = np.maximum(
            0.0,
            tributary_length[penetrating]
            * (
                stiffness_per_length * penetration[penetrating]
                + damping_per_length * free_velocity[penetrating, 1]
            ),
        )
        friction_force = (
            -friction_coefficient
            * normal_force
            * np.tanh(free_velocity[:, 0] / velocity_scale)
        )
        node_force[1:, 0] += friction_force
        node_force[1:, 1] -= normal_force

        active_contact = normal_force > 0.0
        ground_normal_total = float(np.sum(normal_force))
        ground_friction_total = float(np.sum(friction_force))
        ground_friction_power = float(np.sum(friction_force * free_velocity[:, 0]))
        ground_contact_count = int(np.count_nonzero(active_contact))
        max_ground_penetration = (
            float(np.max(penetration[active_contact])) if ground_contact_count else 0.0
        )

    # Santel's first mass point is kinematically the hook. The total force lumped
    # there is passed to the aircraft. The scalar hook tension remains the axial
    # tension in the first cable element for weak-link/release diagnostics.
    hook_force = node_force[0].copy()
    hook_tension = float(internal_tension[0])
    hook_segment_unit = internal_unit[0].copy()

    free_node_acceleration = node_force[1:] / node_mass[1:, None]

    current_length = float(np.sum(span_length))
    current_length_rate = float(
        np.sum(internal_length_rate) + last_span_length_rate
    )
    free_reference_length = float(
        reference_length * (node_count - 1) + final_reference_length_m
    )

    final_strain = last_span_length / final_reference_length_m - 1.0
    max_strain = float(max(np.max(internal_strain), final_strain))
    max_element_tension = float(max(np.max(internal_tension), winch_tension_n))

    total_drag_vector = np.sum(segment_drag_force, axis=0)
    drag_power = float(np.sum(segment_drag_force * span_center_velocity))

    return {
        "hook_force_ground_n": hook_force,
        "hook_resultant_force_n": float(np.linalg.norm(hook_force)),
        "hook_tension_n": hook_tension,
        "hook_unit_ground": hook_segment_unit,
        "hook_segment_unit_ground": hook_segment_unit,
        "winch_tension_n": winch_tension_n,
        "free_node_acceleration_m_s2": free_node_acceleration,
        "current_length_m": current_length,
        "current_length_rate_m_s": current_length_rate,
        "free_reference_length_m": free_reference_length,
        "reference_element_length_m": reference_length,
        "winch_reference_length_m": final_reference_length_m,
        "winch_span_length_m": last_span_length,
        "winch_span_length_rate_m_s": last_span_length_rate,
        "max_element_tension_n": max_element_tension,
        "max_strain": max_strain,
        "min_segment_length_m": float(np.min(span_length)),
        "active_node_count": int(node_count),
        "free_cable_mass_kg": float(fem["linear_density_kg_m"] * free_reference_length),
        "integrated_node_mass_kg": float(np.sum(node_mass)),
        "cable_weight_n": float(np.sum(node_mass) * gravity_m_s2) if include_weight else 0.0,
        "cable_drag_force_ground_n": total_drag_vector,
        "cable_drag_n": float(np.linalg.norm(total_drag_vector)),
        "cable_drag_power_w": drag_power,
        "ground_contact_node_count": int(ground_contact_count),
        "max_ground_penetration_m": float(max_ground_penetration),
        "ground_normal_force_n": float(ground_normal_total),
        "ground_friction_force_n": float(ground_friction_total),
        "ground_friction_power_w": float(ground_friction_power),
    }
