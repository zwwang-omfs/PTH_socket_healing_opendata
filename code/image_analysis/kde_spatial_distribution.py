"""Group KDE curves of marker position within the extraction socket.

Fig. 4B / S4B (EMCN+ vessels) and Fig. 4H / S6B (Osx+ cells).

For every image, each marker pixel inside the socket ROI is assigned its
distance to the socket wall. A Gaussian KDE is computed per image, and the
curves are averaged within each group (shaded band: SEM across images).
Two versions are produced:

* normalized: rank of the distance among all ROI pixels of that image (0 = wall,
  1 = socket centre), so that sockets of different size are comparable;
* absolute: the raw distance.

Usage (from the repository root):
    python code/image_analysis/kde_spatial_distribution.py --marker emcn
    python code/image_analysis/kde_spatial_distribution.py --marker osx
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from socket_margin import DATA_DIR, RESULTS_DIR, osx_socket_id
from socket_margin.core import collect_marker_pixels, load_osx_mask, load_vessel_mask
from socket_margin.density import group_mean_density
from socket_margin.plots import GROUP_ORDER, plot_group_kde

MARKERS = {
    # sample_per_image: maximum pixels drawn per image (EMCN: effectively all).
    "emcn": dict(mask_dir="emcn", suffix="_threshold.png", load=load_vessel_mask,
                 sample_per_image=500000, y_label="EMCN+ Vessel Density"),
    "osx": dict(mask_dir="osx", suffix="_osx.png", load=load_osx_mask,
                sample_per_image=50000, y_label="Osx+ Cell Density"),
}


def absolute_axis(marker: str, pixels: pd.DataFrame):
    """Grid and plotted column for the absolute-distance KDE.

    EMCN: column ``d_raw_um``; the upper limit is the 99.5th percentile of the
    pixel distances. Osx: column ``d_raw`` with a fixed upper limit of 300.
    """
    if marker == "osx":
        upper = 300.0
        return np.linspace(0, upper, 300), "d_raw", upper, "KDE of raw distance"
    raw = pd.to_numeric(pixels["d_raw"], errors="coerce")
    raw = raw[np.isfinite(raw) & (raw > 0)]
    upper = float(np.percentile(raw, 99.5)) if raw.size > 50 else float(raw.max())
    return np.linspace(0, upper, 300), "d_raw_um", upper, "KDE of raw vessel distance"


def run(marker: str, day: str, data_dir: Path, out_dir: Path) -> None:
    cfg = MARKERS[marker]
    socket_dir = data_dir / "emcn" / day
    socket_id = osx_socket_id if marker == "osx" else (lambda sid: sid)
    pixels = collect_marker_pixels(
        data_dir / cfg["mask_dir"] / day,
        cfg["suffix"],
        lambda sid: socket_dir / f"{socket_id(sid)}.png",
        cfg["load"],
        sample_per_image=cfg["sample_per_image"],
    )
    prefix = day[1]  # 'd3' -> '3'
    groups = [f"{prefix}{g}" for g in GROUP_ORDER]
    out_dir.mkdir(parents=True, exist_ok=True)

    pixels.groupby(["group", "filename"]).size().rename("pixels").reset_index().to_csv(
        out_dir / f"{marker}_{day}_pixels_per_image.csv", index=False)

    grid = np.linspace(0.0, 1.0, 200)
    curves = group_mean_density(pixels, "d_cdf", grid, groups, clip_range=(0.0, 1.0))
    curves.to_csv(out_dir / f"{marker}_{day}_kde_normalized.csv", index=False)
    plot_group_kde(curves, prefix, "Normalized Distance", cfg["y_label"],
                   "KDE of normalized distance", out_dir / f"{marker}_{day}_kde_normalized")

    grid, column, upper, title = absolute_axis(marker, pixels)
    curves = group_mean_density(pixels, column, grid, groups, clip_range=(0.0, upper))
    curves.to_csv(out_dir / f"{marker}_{day}_kde_absolute.csv", index=False)
    plot_group_kde(curves, prefix, "Distance (${\\mu}$m)", cfg["y_label"], title,
                   out_dir / f"{marker}_{day}_kde_absolute")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--marker", choices=sorted(MARKERS), required=True)
    parser.add_argument("--days", nargs="+", default=["d3", "d7"], choices=["d3", "d7"])
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--out-dir", type=Path, default=RESULTS_DIR / "kde")
    args = parser.parse_args()
    for day in args.days:
        run(args.marker, day, args.data_dir, args.out_dir)


if __name__ == "__main__":
    main()
