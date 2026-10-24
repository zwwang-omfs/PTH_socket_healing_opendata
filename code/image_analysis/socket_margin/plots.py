"""Figure helpers: group KDE curves and representative distance heatmaps."""

from __future__ import annotations

from pathlib import Path

import cv2
import matplotlib.pyplot as plt
import numpy as np

from .core import (
    ROI_OSX_HEATMAP,
    ROI_VESSEL_HEATMAP,
    UM_PER_PIXEL,
    distance_feature_maps,
    load_osx_mask,
    load_vessel_mask_median,
    read_gray,
    socket_outline_contours,
    socket_roi,
)

# Group suffix -> (legend label, colour, line style), as in Fig. 4B/H.
GROUP_STYLE = {
    "c": ("Ctrl", "#1f77b4", "-"),
    "cy": ("Cy", "#ff7f0e", "--"),
    "cyp": ("Cy + iPTH", "#2ca02c", "-."),
    "p": ("iPTH", "#d62728", ":"),
}
GROUP_ORDER = ["c", "cy", "cyp", "p"]

SCALE_BAR_UM = 500


def _save(fig, out_stem, formats, **kwargs):
    for fmt in formats:
        path = Path(f"{out_stem}.{fmt}")
        fig.savefig(path, **kwargs)
        print(f"Saved {path}")


# ---------------------------------------------------------------------------
# KDE curves
# ---------------------------------------------------------------------------

def plot_group_kde(curves, day_prefix, x_label, y_label, title, out_stem,
                   formats=("pdf", "svg"), font_size=20):
    """Overlay group mean density curves with +/- SEM bands."""
    with plt.rc_context({"font.family": "Arial", "font.size": font_size}):
        fig = plt.figure(figsize=(10, 5))
        for suffix in GROUP_ORDER:
            sub = curves[curves["group"] == f"{day_prefix}{suffix}"]
            if sub.empty:
                continue
            label, colour, style = GROUP_STYLE[suffix]
            x, mean, sem = sub["x"].values, sub["mean"].values, sub["sem"].values
            plt.plot(x, mean, color=colour, label=label, linestyle=style)
            plt.fill_between(x, mean - sem, mean + sem, color=colour, alpha=0.25, linewidth=0)
        plt.xlim(curves["x"].min(), curves["x"].max())
        plt.ylim(bottom=0)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        plt.legend()
        plt.tight_layout()
        _save(fig, out_stem, formats)
        plt.close(fig)


# ---------------------------------------------------------------------------
# Heatmaps
# ---------------------------------------------------------------------------

def _draw_outline_and_scale_bar(ax, socket_path):
    for contour in socket_outline_contours(socket_path):
        points = contour.squeeze()
        if points.ndim == 2 and len(points) > 1:
            ax.plot(points[:, 0], points[:, 1], color="#000000", linewidth=1)
    length_px = SCALE_BAR_UM / UM_PER_PIXEL
    ax.plot([50, 50 + length_px], [50, 50], color="black", linewidth=3, solid_capstyle="butt")
    ax.text(50 + length_px / 2, 70, f"{SCALE_BAR_UM} μm", fontsize=10, ha="center", va="top", color="black")


def vessel_heatmaps(socket_path, vessel_path, out_dir, font_size=12):
    """EMCN+ vessel heatmaps: absolute (pixels) and max-normalized distance.

    Vessel pixels are dilated (5 x 5) for visibility only. Outputs
    ``<id>_absolute_distance_heatmap.pdf`` and
    ``<id>_max_normalized_distance_heatmap.pdf``.
    """
    socket_img = read_gray(socket_path)
    roi, distance = socket_roi(socket_path, **ROI_VESSEL_HEATMAP)
    vessel = load_vessel_mask_median(vessel_path, shape=roi.shape)
    vessel = cv2.dilate(vessel, np.ones((5, 5), np.uint8), iterations=1)
    inside = (vessel > 0) & (roi > 0)
    if not np.any(inside):
        raise ValueError(f"No vessel pixels inside the socket ROI: {vessel_path}")
    heat_abs = np.where(inside, distance, np.nan)
    heat_norm = np.where(inside, distance_feature_maps(roi, distance)["d_norm_max"], np.nan)

    h, w = socket_img.shape
    extent = [0, w, h, 0]
    sample_id = Path(socket_path).stem
    panels = [
        (heat_abs, (0, 300), "Vessel distance heatmap (absolute)",
         "Distance to socket wall (${\\mu}$m)", "absolute_distance_heatmap"),
        (heat_norm, (0, 1), "Vessel distance heatmap (max normalized)",
         "Normalized distance (max)", "max_normalized_distance_heatmap"),
    ]
    with plt.rc_context({"font.size": font_size}):
        for heat, (vmin, vmax), title, cbar_label, suffix in panels:
            fig, ax = plt.subplots(1, 1, figsize=(w / 100, h / 100), dpi=100)
            ax.set_facecolor("white")
            fig.patch.set_facecolor("white")
            ax.imshow(np.full_like(socket_img, 255), cmap="gray", vmin=0, vmax=255, extent=extent, origin="upper")
            image = ax.imshow(np.ma.masked_invalid(heat), cmap="viridis", alpha=1.0, extent=extent,
                              origin="upper", vmin=vmin, vmax=vmax)
            ax.contour(vessel, levels=[0.5], colors="white", linewidths=1.5, extent=extent, origin="upper")
            _draw_outline_and_scale_bar(ax, socket_path)
            ax.set_xlim(0, w)
            ax.set_ylim(h, 0)
            ax.set_title(title, fontsize=font_size)
            ax.set_axis_off()
            cbar = plt.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label(cbar_label, fontsize=12, fontfamily="Arial")
            cbar.ax.tick_params(labelsize=10)
            for label in cbar.ax.get_yticklabels():
                label.set_fontfamily("Arial")
            plt.tight_layout()
            _save(fig, Path(out_dir) / f"{sample_id}_{suffix}", ["pdf"], bbox_inches="tight", dpi=300)
            plt.close(fig)
    return {"absolute": heat_abs, "normalized": heat_norm, "vessel_mask": vessel, "roi": roi}


def osx_heatmaps(socket_path, osx_path, out_dir, font_size=12):
    """Osx+ cell heatmaps: absolute (pixels) and max-normalized distance.

    Cell pixels are dilated (5 x 5) for visibility only. Outputs
    ``<id>_osx_distance_heatmap.{pdf,svg}`` and
    ``<id>_osx_distance_heatmap_normalized.{pdf,svg}``.
    """
    socket_img = read_gray(socket_path)
    roi, distance, _ = socket_roi(socket_path, return_wall=True, **ROI_OSX_HEATMAP)
    cells = load_osx_mask(osx_path, shape=roi.shape)
    cells = cv2.dilate(cells, np.ones((5, 5), np.uint8), iterations=1)
    inside = (cells > 0) & (roi > 0)
    if not np.any(inside):
        raise ValueError(f"No Osx+ pixels inside the socket ROI: {osx_path}")
    max_distance = float(np.max(distance[roi > 0]))
    heat_abs = np.where(inside, distance, np.nan)
    heat_norm = np.where(inside, distance / max_distance, np.nan)

    stem = Path(osx_path).stem
    panels = [
        (heat_abs, (0, 300), "Vessel distance heatmap vs. socket wall",
         "Distance to socket wall (${\\mu}$m)", f"{stem}_distance_heatmap"),
        (heat_norm, (0, 1), f"Normalized distance heatmap (max={max_distance:.1f} μm)",
         "Normalized distance (0-1)", f"{stem}_distance_heatmap_normalized"),
    ]
    with plt.rc_context({"font.size": font_size}):
        for heat, (vmin, vmax), title, cbar_label, name in panels:
            fig, ax = plt.subplots(figsize=(6, 6))
            ax.set_facecolor("white")
            fig.patch.set_facecolor("white")
            ax.imshow(np.full_like(socket_img, 255), cmap="gray", vmin=0, vmax=255)
            image = ax.imshow(np.ma.masked_invalid(heat), cmap="viridis", alpha=1, vmin=vmin, vmax=vmax)
            ax.contour(cells, levels=[0.5], colors="white", linewidths=1)
            _draw_outline_and_scale_bar(ax, socket_path)
            ax.set_title(title, fontsize=font_size)
            ax.set_axis_off()
            cbar = plt.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label(cbar_label, fontsize=font_size)
            plt.tight_layout()
            _save(fig, Path(out_dir) / name, ["pdf", "svg"], bbox_inches="tight", dpi=300)
            plt.close(fig)
    values = heat_norm[inside]
    print(f"{stem}: max distance {max_distance:.2f}; normalized distance "
          f"min {values.min():.3f}, max {values.max():.3f}, mean {values.mean():.3f}, median {np.median(values):.3f}")
    return {"absolute": heat_abs, "normalized": heat_norm, "cell_mask": cells, "roi": roi,
            "max_distance": max_distance}
