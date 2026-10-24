"""Per-sample kernel density estimates and group mean curves (mean +/- SEM)."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


def sample_density(values, grid, clip_range=(0.0, 1.0), bw_method=None):
    """Gaussian KDE of one sample on ``grid``, rescaled to unit area.

    Values outside ``clip_range`` are discarded before estimation.
    """
    values = np.asarray(values)
    values = values[np.isfinite(values)]
    if clip_range is not None:
        lo, hi = clip_range
        values = values[(values >= lo) & (values <= hi)]
    if values.size == 0:
        return np.zeros_like(grid, dtype=float)
    density = gaussian_kde(values, bw_method=bw_method)(grid)
    area = float(np.trapezoid(density, grid))
    if area > 0:
        density = density / area
    return density.astype(float)


def group_mean_density(pixels, value_col, grid, groups, clip_range=(0.0, 1.0),
                       group_col="group", sample_col="filename"):
    """Average the per-sample densities within each group.

    Each image (sample) contributes one curve, so the SEM reflects
    between-sample variation rather than the number of pixels.
    """
    present = set(pixels[group_col].dropna().unique())
    rows = []
    for group in [g for g in groups if g in present]:
        curves = [
            sample_density(sub[value_col].values, grid, clip_range=clip_range)
            for _, sub in pixels[pixels[group_col] == group].groupby(sample_col)
        ]
        if not curves:
            continue
        stacked = np.vstack(curves)
        mean = stacked.mean(axis=0)
        n = stacked.shape[0]
        sem = stacked.std(axis=0, ddof=1) / np.sqrt(n) if n > 1 else np.zeros_like(mean)
        rows.append(pd.DataFrame({"group": group, "x": grid, "mean": mean, "sem": sem, "n": n}))
    if not rows:
        return pd.DataFrame(columns=["group", "x", "mean", "sem", "n"])
    return pd.concat(rows, ignore_index=True)
