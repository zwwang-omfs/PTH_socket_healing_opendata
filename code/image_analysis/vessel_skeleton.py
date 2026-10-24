"""Initial automatic centrelines of EMCN+ vessels (Fig. S5, first step).

For every ``<id>_threshold.png`` (vessels black on white) an
``<id>_skeleton.png`` is written with skimage.morphology.skeletonize(), and an
``<id>_dilated.png`` with cv2.dilate() (white-background dilation, which removes
small black specks). As described in Fig. S5, the automatic skeleton was the
starting point for the semi-manual centreline used for vessel length and
branch-point quantification; that manual correction is not scripted.

Usage (from the repository root):
    python code/image_analysis/vessel_skeleton.py data/images/emcn/d3 --out-dir results/image_analysis/skeleton/d3
"""

from __future__ import annotations

import argparse
from pathlib import Path

import cv2
import numpy as np
from skimage.morphology import skeletonize
from skimage.util import img_as_ubyte

from socket_margin.core import read_gray


def vessel_skeleton(threshold_img: np.ndarray) -> np.ndarray:
    """Skeleton of the black vessels, drawn as black lines on white."""
    _, vessels = cv2.threshold(threshold_img, 127, 255, cv2.THRESH_BINARY_INV)
    return 255 - img_as_ubyte(skeletonize(vessels > 0))


def dilate_background(threshold_img: np.ndarray, kernel_size: int = 1, iterations: int = 1) -> np.ndarray:
    kernel = np.ones((kernel_size, kernel_size), np.uint8)
    return cv2.dilate(threshold_img, kernel, iterations=iterations)


def write_png(path: Path, image: np.ndarray) -> None:
    ok, buffer = cv2.imencode(".png", image)
    if not ok:
        raise RuntimeError(f"Could not encode {path}")
    buffer.tofile(str(path))
    print(f"Saved {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("folder", type=Path, help="folder containing *_threshold.png")
    parser.add_argument("--out-dir", type=Path, default=None, help="default: the input folder")
    parser.add_argument("--dilate-kernel", type=int, default=1)
    args = parser.parse_args()

    out_dir = args.out_dir or args.folder
    out_dir.mkdir(parents=True, exist_ok=True)
    for path in sorted(args.folder.glob("*_threshold.png")):
        image = read_gray(path)
        sample_id = path.name[: -len("_threshold.png")]
        write_png(out_dir / f"{sample_id}_skeleton.png", vessel_skeleton(image))
        write_png(out_dir / f"{sample_id}_dilated.png", dilate_background(image, args.dilate_kernel))


if __name__ == "__main__":
    main()
