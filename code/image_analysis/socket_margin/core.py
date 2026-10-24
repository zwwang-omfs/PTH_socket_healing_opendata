"""Socket ROI reconstruction, mask loading and per-pixel distance features.

The extraction-socket outline (a thin traced line) is closed at its opening, the
enclosed area is taken as the region of interest (ROI), and every ROI pixel is
assigned its Euclidean distance to the socket wall. Marker pixels (EMCN+ vessels
or Osx+ cells) inside the ROI are then described by that distance, either in
pixels or after per-image normalization.
"""

from __future__ import annotations

import re
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
from skimage.morphology import skeletonize

UM_PER_PIXEL = 1.243

# ROI variants used in the manuscript. The KDE analyses (Fig. 4B/H, S4B, S6B)
# used the closed, convex-hull ROI; the representative heatmaps used the
# variants below, exactly as they were generated for Fig. 4A/G, S4A and S6A.
ROI_KDE = dict(dilate_wall=True, close_wall=True, convex_hull=True)
ROI_VESSEL_HEATMAP = dict(dilate_wall=False, close_wall=False, convex_hull=False)
ROI_OSX_HEATMAP = dict(dilate_wall=True, close_wall=True, convex_hull=False)

GROUP_PATTERN = r"(\d\D*?)(?=\d)"


def read_gray(path: str | Path) -> np.ndarray:
    """Read an image as 8-bit grayscale; works with non-ASCII paths on Windows."""
    data = np.fromfile(str(path), dtype=np.uint8)
    img = cv2.imdecode(data, cv2.IMREAD_GRAYSCALE)
    if img is None:
        raise ValueError(f"Cannot read image: {path}")
    return img


def sample_group(sample_id: str) -> str | None:
    """Return the day+group prefix of a sample ID, e.g. '3cyp11_2' -> '3cyp'."""
    match = re.search(GROUP_PATTERN, sample_id)
    return match.group(1) if match else None


# ---------------------------------------------------------------------------
# Socket ROI
# ---------------------------------------------------------------------------

def _find_skeleton_endpoints(skeleton01: np.ndarray) -> list[tuple[int, int]]:
    """Skeleton pixels with exactly one 8-connected neighbour, as (x, y)."""
    kernel = np.array([[1, 1, 1], [1, 10, 1], [1, 1, 1]], np.int32)
    neighbours = cv2.filter2D(skeleton01, -1, kernel, borderType=cv2.BORDER_CONSTANT)
    ys, xs = np.where(neighbours == 11)
    return [(int(x), int(y)) for x, y in zip(xs, ys)]


def _estimate_wall_thickness_px(wall255: np.ndarray) -> int:
    inverted = (255 - wall255).astype(np.uint8)
    distance = cv2.distanceTransform(inverted, cv2.DIST_L2, 5)
    edges = cv2.Canny(wall255, 50, 150)
    values = distance[edges > 0]
    if values.size == 0:
        return 3
    return int(np.clip(2 * np.median(values), 2, 12))


def _pick_opening_endpoints(endpoints, wall01):
    """Choose the endpoint pair that is closest together without crossing the wall."""
    if len(endpoints) < 2:
        return None
    wall = wall01.astype(bool)
    canvas = np.zeros(wall01.shape, np.uint8)
    best, best_cost = None, 1e18
    for i in range(len(endpoints)):
        for j in range(i + 1, len(endpoints)):
            p1, p2 = endpoints[i], endpoints[j]
            canvas.fill(0)
            cv2.line(canvas, p1, p2, 1, 1)
            crossings = np.count_nonzero((canvas > 0) & wall)
            cost = np.hypot(p1[0] - p2[0], p1[1] - p2[1]) + 1e6 * crossings
            if cost < best_cost:
                best_cost, best = cost, (p1, p2)
    return best


def socket_roi(socket_image_path, *, dilate_wall=True, close_wall=True,
               convex_hull=True, return_wall=False):
    """Reconstruct the socket ROI from a traced outline image.

    Returns the ROI (0/1 uint8) and the distance-to-wall map in pixels
    (float32); with ``return_wall`` also the closed wall mask (0/1 uint8).
    """
    socket_img = read_gray(socket_image_path)
    _, wall = cv2.threshold(socket_img, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    if cv2.countNonZero(wall) > wall.size // 2:
        wall = cv2.bitwise_not(wall)
    if dilate_wall:
        wall = cv2.dilate(wall, None, iterations=1)

    # Close the socket opening by joining the two free ends of the outline.
    skeleton = skeletonize((wall > 0).astype(bool)).astype(np.uint8)
    pair = _pick_opening_endpoints(_find_skeleton_endpoints(skeleton), (wall > 0).astype(np.uint8))
    if pair is not None:
        thickness = _estimate_wall_thickness_px(wall)
        cv2.line(wall, pair[0], pair[1], 255, thickness=max(2, thickness))

    if close_wall:
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
        wall = cv2.morphologyEx(wall, cv2.MORPH_CLOSE, kernel)

    # Everything not reachable from the image corner is inside the socket.
    free_space = (wall == 0).astype(np.uint8)
    h, w = free_space.shape
    flood = free_space * 255
    cv2.floodFill(flood, np.zeros((h + 2, w + 2), np.uint8), (0, 0), 128)
    outside = (flood == 128).astype(np.uint8)
    roi = (free_space & (1 - outside)).astype(np.uint8)

    if convex_hull:
        contours, _ = cv2.findContours(roi, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if contours:
            hull = cv2.convexHull(np.vstack(contours))
            roi.fill(0)
            cv2.drawContours(roi, [hull], -1, 1, thickness=-1)

    distance = cv2.distanceTransform(roi * 255, cv2.DIST_L2, 5).astype(np.float32)
    if return_wall:
        return roi, distance, (wall > 0).astype(np.uint8)
    return roi, distance


def socket_outline_contours(socket_image_path):
    """Contours of the traced socket outline, for drawing on figures."""
    socket_img = read_gray(socket_image_path)
    _, bw = cv2.threshold(socket_img, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    if cv2.countNonZero(bw) > bw.size // 2:
        bw = cv2.bitwise_not(bw)
    outline = (cv2.bitwise_not(bw) < 255).astype(np.uint8)
    contours, _ = cv2.findContours(outline, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)
    return contours


# ---------------------------------------------------------------------------
# Marker masks
# ---------------------------------------------------------------------------

def _fit_shape(mask01, shape):
    if shape is not None and mask01.shape != shape:
        mask01 = cv2.resize(mask01, (shape[1], shape[0]), interpolation=cv2.INTER_NEAREST)
    return mask01.astype(np.uint8)


def load_vessel_mask(path, shape=None, kernel_size=3):
    """EMCN+ vessels are black (0) in the threshold image.

    Small black specks are removed by dilating the white background.
    """
    mask01 = (read_gray(path) == 0).astype(np.uint8)
    ksize = int(kernel_size) if kernel_size and kernel_size > 0 else 3
    if ksize % 2 == 0:
        ksize += 1
    background = cv2.dilate(1 - mask01, np.ones((ksize, ksize), np.uint8), iterations=1)
    return _fit_shape(1 - background, shape)


def load_vessel_mask_median(path, shape=None, ksize=5):
    """Vessel mask denoised with a median filter (used for the vessel heatmaps)."""
    mask01 = cv2.medianBlur((read_gray(path) == 0).astype(np.uint8), ksize)
    return _fit_shape(mask01, shape)


def load_osx_mask(path, shape=None):
    """Osx+ cells are white (255) in the threshold image."""
    return _fit_shape((read_gray(path) == 255).astype(np.uint8), shape)


# ---------------------------------------------------------------------------
# Per-pixel distance features
# ---------------------------------------------------------------------------

def distance_feature_maps(roi01, distance):
    """Raw distance plus three per-image normalizations, NaN outside the ROI.

    d_norm_max    distance / maximum distance in the ROI
    d_cdf         rank of the distance within all ROI pixels (0-1)
    d_over_median distance / median distance in the ROI
    """
    roi = roi01 > 0
    d = np.where(roi, distance, np.nan).astype(float)
    values = d[roi]
    out = {"d_raw": d}
    if values.size == 0:
        nan = d * 0 + np.nan
        out.update(d_norm_max=nan, d_cdf=nan, d_over_median=nan,
                   edge_half_flag=np.zeros_like(d, dtype=np.uint8))
        return out
    d_max = float(np.nanmax(values)) if np.isfinite(values).any() else np.nan
    out["d_norm_max"] = d / d_max if d_max and np.isfinite(d_max) and d_max > 0 else d * 0 + np.nan
    ranks = np.searchsorted(np.sort(values), d, side="right")
    out["d_cdf"] = np.where(roi, ranks / max(1, values.size), np.nan).astype(float)
    median = float(np.nanmedian(values)) if np.isfinite(values).any() else np.nan
    out["d_over_median"] = d / median if median and np.isfinite(median) and median > 0 else d * 0 + np.nan
    edge = (d < median).astype(np.uint8)
    edge[~roi] = 0
    out["edge_half_flag"] = edge
    return out


def marker_pixel_table(socket_path, mask_path, load_mask, sample=None, rng_seed=0):
    """One row per marker pixel inside the ROI, with its distance features.

    If ``sample`` is set and the image has more marker pixels, a reproducible
    random subset of that size is drawn.
    """
    roi, distance = socket_roi(socket_path, **ROI_KDE)
    marker = (load_mask(mask_path, shape=roi.shape) > 0) & (roi > 0)
    features = distance_feature_maps(roi, distance)
    yy, xx = np.where(marker)
    if sample is not None and sample > 0 and yy.size > sample:
        keep = np.random.default_rng(rng_seed).choice(np.arange(yy.size), size=sample, replace=False)
        yy, xx = yy[keep], xx[keep]
    d_raw = features["d_raw"][yy, xx]
    return pd.DataFrame({
        "y": yy.astype(int),
        "x": xx.astype(int),
        "d_raw": d_raw,
        "d_raw_um": d_raw * UM_PER_PIXEL,
        "d_norm_max": features["d_norm_max"][yy, xx],
        "d_cdf": features["d_cdf"][yy, xx],
        "d_over_median": features["d_over_median"][yy, xx],
        "edge_half_flag": features["edge_half_flag"][yy, xx].astype(np.uint8),
    })


def collect_marker_pixels(mask_dir, mask_suffix, socket_path_for, load_mask,
                          sample_per_image=None, rng_seed=0):
    """Pixel table for every ``*<mask_suffix>`` image in ``mask_dir``.

    ``socket_path_for(sample_id)`` returns the matching socket outline image.
    Images are processed in sorted order and image *i* uses seed
    ``rng_seed + i``, so subsampling is reproducible.
    """
    mask_dir = Path(mask_dir)
    if not mask_dir.is_dir():
        raise FileNotFoundError(f"Folder not found: {mask_dir}")
    mask_files = sorted(f.name for f in mask_dir.iterdir() if f.name.endswith(mask_suffix))
    tables = []
    for i, name in enumerate(mask_files):
        sample_id = name[: -len(mask_suffix)]
        socket_path = Path(socket_path_for(sample_id))
        if not socket_path.exists():
            print(f"Skipping {name}: socket outline not found ({socket_path})")
            continue
        table = marker_pixel_table(socket_path, mask_dir / name, load_mask,
                                   sample=sample_per_image, rng_seed=rng_seed + i)
        if table.empty:
            continue
        table.insert(0, "group", sample_group(sample_id))
        table.insert(1, "filename", sample_id)
        tables.append(table)
    if not tables:
        return pd.DataFrame(columns=["group", "filename", "y", "x"])
    pixels = pd.concat(tables, ignore_index=True)
    print(f"{mask_dir}: {pixels['filename'].nunique()} images, "
          f"{pixels['group'].nunique()} groups, {len(pixels)} pixels")
    return pixels
