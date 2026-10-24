"""Representative heatmaps of marker distance to the socket wall.

Fig. 4A / S4A (EMCN+ vessels) and Fig. 4G / S6A (Osx+ cells). For each
sample, an absolute-distance and a max-normalized heatmap are written.

Usage (from the repository root):
    python code/image_analysis/distance_heatmaps.py --marker emcn --day d7 --samples 7p11_2
    python code/image_analysis/distance_heatmaps.py --marker osx --day d3 --samples 3p11_1
    python code/image_analysis/distance_heatmaps.py --marker osx --day d7      # all samples
"""

from __future__ import annotations

import argparse
from pathlib import Path

from socket_margin import DATA_DIR, RESULTS_DIR, osx_socket_id
from socket_margin.plots import osx_heatmaps, vessel_heatmaps

SUFFIX = {"emcn": "_threshold.png", "osx": "_osx.png"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--marker", choices=sorted(SUFFIX), required=True)
    parser.add_argument("--day", choices=["d3", "d7"], required=True)
    parser.add_argument("--samples", nargs="*", help="sample IDs, e.g. 3p11_1 (default: all)")
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--out-dir", type=Path, default=RESULTS_DIR / "heatmaps")
    args = parser.parse_args()

    mask_dir = args.data_dir / args.marker / args.day
    socket_dir = args.data_dir / "emcn" / args.day
    suffix = SUFFIX[args.marker]
    samples = args.samples or sorted(p.name[: -len(suffix)] for p in mask_dir.glob(f"*{suffix}"))
    out_dir = args.out_dir / args.marker / args.day
    out_dir.mkdir(parents=True, exist_ok=True)

    for sample_id in samples:
        mask_path = mask_dir / f"{sample_id}{suffix}"
        if args.marker == "emcn":
            vessel_heatmaps(socket_dir / f"{sample_id}.png", mask_path, out_dir)
        else:
            osx_heatmaps(socket_dir / f"{osx_socket_id(sample_id)}.png", mask_path, out_dir)


if __name__ == "__main__":
    main()
