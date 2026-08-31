import matplotlib.pyplot as plt
import numpy as np

from winch_launch_simulator.simulation import load_case, run_simulation
from winch_launch_simulator.visualization import make_winch_launch_animation, plot_summary


def _cable_collection(axes):
    return next(
        collection
        for collection in axes["trajectory"].collections
        if collection.get_label() == "tow cable"
    )


def test_summary_plot_has_four_axes():
    case = load_case("cases/g103a_prescribed_tension.json")
    case["simulation"]["max_time_s"] = 0.05
    result = run_simulation(case)
    figure = plot_summary(result)
    assert len(figure.axes) == 4
    plt.close(figure)


def test_animation_frame_count_tracks_simulation_time_not_integrator_samples():
    case = load_case("cases/g103a_prescribed_tension.json")
    case["simulation"]["max_time_s"] = 2.0
    result = run_simulation(case)
    fps = 20

    figure, ani, axes = make_winch_launch_animation(
        result["history"], fps=fps, mode="overview"
    )
    frame_count = len(list(ani.new_frame_seq()))
    duration = result["history"]["time_s"][-1] - result["history"]["time_s"][0]

    assert "trajectory" in axes
    assert frame_count == max(1, int(round(duration * fps)))
    assert abs(frame_count / fps - duration) <= 0.5 / fps + 1.0e-12
    assert frame_count < len(result["history"]["time_s"])

    ani._draw_was_started = True
    plt.close(figure)


def _pixels_per_metre(axis):
    # transData is the final data-to-display transform used by the renderer.
    # Measuring it after draw() verifies the actual pixel scale rather than
    # inferring it from requested Figure/Axes geometry.
    axis.figure.canvas.draw()
    origin = axis.transData.transform((0.0, 0.0))
    x_point = axis.transData.transform((100.0, 0.0))
    z_point = axis.transData.transform((0.0, 100.0))
    return (
        abs(x_point[0] - origin[0]) / 100.0,
        abs(z_point[1] - origin[1]) / 100.0,
    )


def test_animation_keeps_equal_x_z_physical_scale_in_follow_and_overview_modes():
    case = load_case("cases/g103a_dynamic_winch.json")
    case["simulation"]["max_time_s"] = 0.05
    result = run_simulation(case)

    for mode in ("follow", "overview"):
        figure, ani, axes = make_winch_launch_animation(
            result["history"], fps=20, mode=mode
        )
        ani._func(0)
        x_pixels_per_metre, z_pixels_per_metre = _pixels_per_metre(
            axes["trajectory"]
        )
        np.testing.assert_allclose(
            x_pixels_per_metre,
            z_pixels_per_metre,
            rtol=0.0,
            atol=1.0e-9,
        )

        ani._draw_was_started = True
        plt.close(figure)


def test_overview_keeps_ground_near_bottom_at_equal_rendered_scale():
    history = {
        "time_s": [0.0, 1.0],
        "x_m": [0.0, 1100.0],
        "relative_altitude_m": [0.0, 350.0],
        "hook_x_m": [0.0, 1100.0],
        "hook_z_m": [0.0, -350.0],
        "winch_x_m": [1100.0, 1100.0],
        "winch_z_m": [0.0, 0.0],
        "connection_state": [1, 1],
        "cable_model": ["secant", "secant"],
    }

    figure, ani, axes = make_winch_launch_animation(
        history, fps=20, mode="overview"
    )
    ani._func(0)
    axis = axes["trajectory"]

    assert axis.get_ylim()[0] == -20.0
    assert axis.get_ylim()[1] < 400.0
    assert figure.get_figwidth() > figure.get_figheight()

    x_pixels_per_metre, z_pixels_per_metre = _pixels_per_metre(axis)
    np.testing.assert_allclose(
        x_pixels_per_metre,
        z_pixels_per_metre,
        rtol=0.0,
        atol=1.0e-9,
    )

    ani._draw_was_started = True
    plt.close(figure)


def test_secant_animation_uses_one_straight_cable_segment():
    case = load_case("cases/g103a_dynamic_winch.json")
    case["simulation"]["max_time_s"] = 0.05
    result = run_simulation(case)

    figure, ani, axes = make_winch_launch_animation(
        result["history"], fps=20, mode="overview"
    )
    ani._func(0)
    segments = _cable_collection(axes).get_segments()

    assert len(segments) == 1
    np.testing.assert_allclose(
        segments[0][0],
        [result["history"]["hook_x_m"][0], -result["history"]["hook_z_m"][0]],
    )
    np.testing.assert_allclose(
        segments[0][1],
        [result["history"]["winch_x_m"][0], -result["history"]["winch_z_m"][0]],
    )

    ani._draw_was_started = True
    plt.close(figure)


def test_fem_animation_draws_each_active_element_from_saved_node_positions():
    case = load_case("cases/g103a_dynamic_winch_fem.json")
    case["simulation"]["max_time_s"] = 0.01
    result = run_simulation(case)
    history = result["history"]

    assert "cable_node_2_x_m" in history
    assert "cable_node_20_z_m" in history
    assert np.isfinite(history["cable_node_2_x_m"][0])
    assert np.isfinite(history["cable_node_20_z_m"][0])

    figure, ani, axes = make_winch_launch_animation(history, fps=20, mode="overview")
    ani._func(0)
    segments = _cable_collection(axes).get_segments()

    assert len(segments) == int(history["fem_active_node_count"][0])

    hook_x = history["hook_x_m"][0]
    hook_altitude = -history["hook_z_m"][0]
    winch_x = history["winch_x_m"][0]
    winch_altitude = -history["winch_z_m"][0]
    node_x = np.array([history[f"cable_node_{n}_x_m"][0] for n in range(2, 21)])
    node_altitude = -np.array([history[f"cable_node_{n}_z_m"][0] for n in range(2, 21)])
    los_altitude = (
        hook_altitude
        + (winch_altitude - hook_altitude)
        * (node_x - hook_x)
        / (winch_x - hook_x)
    )
    assert np.max(np.abs(node_altitude - los_altitude)) > 0.01

    ani._draw_was_started = True
    plt.close(figure)


def test_fem_animation_drops_removed_nodes_and_hides_released_cable():
    history = {
        "time_s": [0.0, 1.0, 2.0, 3.0],
        "x_m": [0.0, 1.0, 2.0, 3.0],
        "relative_altitude_m": [1.0, 2.0, 3.0, 3.0],
        "hook_x_m": [0.0, 1.0, 2.0, 3.0],
        "hook_z_m": [-1.0, -2.0, -3.0, -3.0],
        "winch_x_m": [10.0, 10.0, 10.0, 10.0],
        "winch_z_m": [0.0, 0.0, 0.0, 0.0],
        "connection_state": [1, 1, 3, 3],
        "cable_model": ["fem", "fem", "fem", "fem"],
        "fem_active_node_count": [3, 2, 2, 2],
        "cable_node_2_x_m": [3.0, 5.0, 5.0, 5.0],
        "cable_node_2_z_m": [0.2, 0.1, 0.1, 0.1],
        "cable_node_3_x_m": [7.0, np.nan, np.nan, np.nan],
        "cable_node_3_z_m": [0.4, np.nan, np.nan, np.nan],
    }

    figure, ani, axes = make_winch_launch_animation(history, fps=1, mode="overview")
    cable = _cable_collection(axes)

    ani._func(0)
    assert len(cable.get_segments()) == 3

    ani._func(1)
    assert len(cable.get_segments()) == 2

    ani._func(2)
    assert len(cable.get_segments()) == 0

    ani._draw_was_started = True
    plt.close(figure)
