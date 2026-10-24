"""Spatial distribution of EMCN+ vessels and Osx+ cells relative to the socket wall."""

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = REPO_ROOT / "data" / "images"
RESULTS_DIR = REPO_ROOT / "results" / "image_analysis"


def osx_socket_id(sample_id: str) -> str:
    """Socket-outline ID for an Osx image.

    Day-7 iPTH Osx images are named ``7p<mouse>_<root>`` whereas the matching
    socket outlines are named ``7p1<mouse>_<root>``; all other IDs are shared.
    """
    return re.sub(r"^7p(\d)_", r"7p1\1_", sample_id)
