import importlib
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter, writers
from matplotlib.collections import LineCollection

def plot_summary(result):
    """Create the standard compact static review plot from one simulation history."""
    history = result["history"]
    time = history["time_s"]

    figure, axes = plt.subplots(4, 1, figsize=(10, 12))

    axes[0].plot(history["x_m"], history["relative_altitude_m"])
    axes[0].set_xlabel("x [m]")
    axes[0].set_ylabel("relative altitude [m]")
    axes[0].grid(True)

    axes[1].plot(time, history["eas_m_s"], label="EAS")
    axes[1].plot(time, history["airspeed_m_s"], label="TAS")
    axes[1].set_ylabel("speed [m/s]")
    axes[1].legend()
    axes[1].grid(True)

    axes[2].plot(time, history["theta_deg"], label="theta")
    axes[2].plot(time, history["gamma_air_deg"], label="gamma_air")
    axes[2].plot(time, history["alpha_deg"], label="alpha")
    axes[2].set_ylabel("angle [deg]")
    axes[2].legend()
    axes[2].grid(True)

    target_tension = np.asarray(history["target_tension_n"], dtype=float)
    if np.isfinite(target_tension).any():
        axes[3].plot(time, target_tension, label="target")
    axes[3].plot(time, history["cable_tension_n"], label="hook")
    axes[3].set_xlabel("time [s]")
    axes[3].set_ylabel("tension [N]")
    axes[3].legend()
    axes[3].grid(True)

    figure.suptitle(result["metadata"]["case_id"])
    figure.tight_layout()
    return figure

def make_winch_launch_animation(history, fps=30, mode="follow"):
    """Create a 1x-real-time, equal-x/z-scale animation from calculated history."""
    if fps <= 0:
        raise ValueError("fps must be positive.")

    time = np.asarray(history["time_s"], dtype=float)
    x = np.asarray(history["x_m"], dtype=float)
    altitude = np.asarray(history["relative_altitude_m"], dtype=float)
    hook_x = np.asarray(history["hook_x_m"], dtype=float)
    hook_altitude = -np.asarray(history["hook_z_m"], dtype=float)
    winch_x = np.asarray(history["winch_x_m"], dtype=float)
    winch_altitude = -np.asarray(history["winch_z_m"], dtype=float)
    connection_state = np.asarray(history["connection_state"])
    cable_model = str(history["cable_model"][0])

    fem_node_numbers = sorted(
        int(key.removeprefix("cable_node_").removesuffix("_x_m"))
        for key in history
        if key.startswith("cable_node_") and key.endswith("_x_m")
    )

    duration = time[-1] - time[0]
    frame_count = max(1, int(round(duration * fps)))
    frame_times = time[0] + np.arange(frame_count, dtype=float) / fps
    frame_indices = np.searchsorted(
        time,
        frame_times,
        side="right",
    ) - 1
    frame_indices = np.clip(frame_indices, 0, len(time) - 1)
    frame_indices[0] = 0

    # Let Matplotlib enforce the x-z physical scale directly.  The data
    # limits below describe only the useful launch region; equal aspect is
    # responsible for making 1 m in x the same rendered length as 1 m in z.
    # This avoids coupling the physical scale to hand-calculated Figure/Axes
    # ratios.
    if mode == "overview":
        figure_width_in = 14.0
        figure_height_in = 5.5
    else:
        figure_width_in = 10.0
        figure_height_in = 7.0

    figure = plt.figure(figsize=(figure_width_in, figure_height_in))
    axis = figure.add_axes([0.09, 0.11, 0.88, 0.84])
    axis.set_aspect("equal", adjustable="box")

    path_line, = axis.plot([], [])
    aircraft_point, = axis.plot([], [], marker="o")
    cable_segments = LineCollection([], linewidths=1.5, label="tow cable")
    axis.add_collection(cable_segments)
    time_text = axis.text(0.02, 0.95, "", transform=axis.transAxes)
    axis.set_xlabel("x [m]")
    axis.set_ylabel("relative altitude [m]")
    axis.grid(True)

    # Overview limits describe the useful physical region only.  Keep the
    # ordinary ground-side limit at -20 m instead of inventing hundreds of
    # metres of underground space to obtain a desired aspect ratio.
    overview_x_min = min(
        0.0,
        float(np.min(x)),
        float(np.min(hook_x)),
        float(np.min(winch_x)),
    )
    overview_x_max = max(
        float(np.max(x)),
        float(np.max(hook_x)),
        float(np.max(winch_x)),
    )
    overview_h_data_min = min(
        float(np.min(altitude)),
        float(np.min(hook_altitude)),
        float(np.min(winch_altitude)),
    )
    overview_h_max = max(
        float(np.max(altitude)),
        float(np.max(hook_altitude)),
        float(np.max(winch_altitude)),
    )

    x_margin = 0.05 * max(overview_x_max - overview_x_min, 1.0)
    h_margin = 0.05 * max(overview_h_max - overview_h_data_min, 1.0)
    overview_xlim = (
        overview_x_min - x_margin,
        overview_x_max + x_margin,
    )
    overview_h_min = min(-20.0, overview_h_data_min)
    if overview_h_data_min < -20.0:
        overview_h_min -= h_margin
    overview_ylim = (
        overview_h_min,
        overview_h_max + h_margin,
    )

    # Follow mode keeps a fixed 400 m x-span and 280 m z-span.  Because the
    # ratio of those spans is constant, equal aspect also keeps the Axes box
    # stable while the window follows the aircraft.
    follow_x_before_m = 150.0
    follow_x_after_m = 250.0
    follow_h_span = 280.0
    follow_fraction_below_aircraft = 100.0 / follow_h_span

    def update(frame_number):
        index = frame_indices[frame_number]
        path_line.set_data(x[: index + 1], altitude[: index + 1])
        aircraft_point.set_data([x[index]], [altitude[index]])

        if connection_state[index] in (1, 2):
            points = [[hook_x[index], hook_altitude[index]]]

            if cable_model == "fem":
                active_node_count = int(history["fem_active_node_count"][index])
                for node_number in fem_node_numbers:
                    if node_number > active_node_count:
                        break
                    node_x = float(history[f"cable_node_{node_number}_x_m"][index])
                    node_z = float(history[f"cable_node_{node_number}_z_m"][index])
                    if np.isfinite(node_x) and np.isfinite(node_z):
                        points.append([node_x, -node_z])

            points.append([winch_x[index], winch_altitude[index]])
            points = np.asarray(points, dtype=float)
            cable_segments.set_segments(
                np.stack((points[:-1], points[1:]), axis=1)
            )
        else:
            cable_segments.set_segments([])

        if mode == "follow":
            axis.set_xlim(
                x[index] - follow_x_before_m,
                x[index] + follow_x_after_m,
            )

            follow_h_min = (
                altitude[index]
                - follow_fraction_below_aircraft * follow_h_span
            )
            if follow_h_min < -20.0:
                follow_h_min = -20.0
            axis.set_ylim(follow_h_min, follow_h_min + follow_h_span)
        else:
            axis.set_xlim(*overview_xlim)
            axis.set_ylim(*overview_ylim)

        time_text.set_text(f"t = {time[index]:.1f} s")
        return path_line, aircraft_point, cable_segments, time_text

    ani = FuncAnimation(
        figure,
        update,
        frames=len(frame_indices),
        interval=1000.0 / fps,
        blit=False,
        repeat=True,
    )

    axes = {
        "trajectory": axis,
    }
    return figure, ani, axes

def configure_ffmpeg(auto_install=True):
    """Configure Matplotlib to use FFmpeg, matching the animation notebooks."""
    if writers.is_available("ffmpeg"):
        return plt.rcParams.get("animation.ffmpeg_path", "ffmpeg")

    try:
        import imageio_ffmpeg
    except ImportError:
        if not auto_install:
            return None

        print(
            "FFmpegが見つからないため、imageio-ffmpeg を"
            "現在のPython環境へインストールします。"
        )
        try:
            subprocess.check_call(
                [sys.executable, "-m", "pip", "install", "-q", "imageio-ffmpeg"]
            )
            importlib.invalidate_caches()
            import imageio_ffmpeg
        except Exception as exc:
            print(f"imageio-ffmpeg の準備に失敗しました: {exc}")
            return None

    ffmpeg_path = imageio_ffmpeg.get_ffmpeg_exe()
    plt.rcParams["animation.ffmpeg_path"] = ffmpeg_path

    if writers.is_available("ffmpeg"):
        print(f"FFmpeg: {ffmpeg_path}")
        return ffmpeg_path

    return None

def save_mp4(
    ani,
    path="winch_launch_animation.mp4",
    fps=30,
    dpi=120,
):
    """Save an animation as MP4."""
    if configure_ffmpeg(auto_install=True) is None:
        raise RuntimeError(
            "MP4保存に必要なFFmpegを利用できません。"
            "FFmpegまたはimageio-ffmpegを利用できる状態にしてください。"
        )

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    writer = FFMpegWriter(fps=fps)
    ani.save(path, writer=writer, dpi=dpi)
    return path

def save_gif(
    ani,
    path="winch_launch_animation.gif",
    fps=30,
    dpi=100,
):
    """Save an animation as GIF."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    writer = PillowWriter(fps=fps)
    ani.save(path, writer=writer, dpi=dpi)
    return path
