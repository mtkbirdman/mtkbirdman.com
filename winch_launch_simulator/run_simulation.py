import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from winch_launch_simulator.simulation import load_case, run_simulation, save_result
from winch_launch_simulator.visualization import (
    make_winch_launch_animation,
    plot_summary,
    save_gif,
    save_mp4,
)

def main():
    parser = argparse.ArgumentParser(description="Run one winch-launch simulation case.")
    parser.add_argument("case", help="Path to a JSON case file.")
    parser.add_argument("--output-root", default="outputs")
    parser.add_argument("--animation", action="store_true", help="Save MP4 animation.")
    parser.add_argument("--gif", action="store_true", help="Also save GIF animation.")
    parser.add_argument(
        "--fps",
        type=int,
        default=30,
        help="Animation frames per simulated real-time second.",
    )
    parser.add_argument("--overview", action="store_true", help="Use overview animation framing.")
    args = parser.parse_args()

    case = load_case(args.case)
    result = run_simulation(case)
    output_directory = save_result(result, args.output_root)

    figure = plot_summary(result)
    figure.savefig(Path(output_directory) / f"{case['case_id']}_summary.png", dpi=150)
    plt.close(figure)

    if args.animation or args.gif or case["output"].get("save_animation", False):
        figure, ani, _ = make_winch_launch_animation(
            result["history"],
            fps=args.fps,
            mode="overview" if args.overview else "follow",
        )

        if args.animation or case["output"].get("save_animation", False):
            save_mp4(
                ani,
                path=Path(output_directory) / f"{case['case_id']}_winch_launch.mp4",
                fps=args.fps,
                dpi=120,
            )

        if args.gif:
            save_gif(
                ani,
                path=Path(output_directory) / f"{case['case_id']}_winch_launch.gif",
                fps=args.fps,
                dpi=100,
            )

        plt.close(figure)

    print(f"case_id      : {case['case_id']}")
    print(f"samples      : {result['metadata']['samples']}")
    print(f"final_time_s : {result['metadata']['final_time_s']:.3f}")
    print(f"final_alt_m  : {result['metadata']['final_altitude_m']:.3f}")
    if result["events"]:
        print("events:")
        for event in result["events"]:
            print(f"  {event['time']:.3f} s  {event['code']}  [{event['severity']}]")
    else:
        print("events       : none")

if __name__ == "__main__":
    main()
