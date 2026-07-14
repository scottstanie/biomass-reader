#!/usr/bin/env python
"""Make compact, consistently scaled figures from a completed Dolphin run.

The interferogram panel removes each pair's circular-mean phase before display.
That preserves its fringe pattern while giving every panel the same phase origin
instead of assigning an arbitrary red/blue offset to individual pairs.

Example
-------
python scripts/plot_dolphin_diagnostics.py \
    /data/biomass/gslc_10stack/dolphin
"""

from __future__ import annotations

import argparse
import math
import re
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import rasterio
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap
from rasterio.enums import Resampling

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", message="Failed to write rc setting")
    import ultraplot as uplt

OIL_SLICK_NODES = [
    (0.02, 0.01, 0.06),
    (0.22, 0.00, 0.48),
    (0.02, 0.08, 0.82),
    (0.00, 0.58, 0.88),
    (0.00, 0.72, 0.32),
    (0.72, 0.92, 0.00),
    (1.00, 0.82, 0.00),
    (1.00, 0.40, 0.00),
    (0.88, 0.02, 0.18),
    (0.52, 0.00, 0.38),
    (0.02, 0.01, 0.06),
]
OIL_SLICK_CMAP = LinearSegmentedColormap.from_list("oil_slick", OIL_SLICK_NODES, 512)

PAIR = re.compile(r"(\d{8})_(\d{8})\.int(?:\.cor)?\.tif$")
CLOSURE = re.compile(r"closure_phase_(\d{8})_(\d{8})_(\d{8})\.tif$")


def _read_display(path: Path, max_dimension: int = 600) -> tuple[np.ndarray, float]:
    """Read a raster at a compact display resolution using area averaging."""
    with rasterio.open(path) as dataset:
        scale = max(1, math.ceil(max(dataset.shape) / max_dimension))
        shape = (
            math.ceil(dataset.height / scale),
            math.ceil(dataset.width / scale),
        )
        return (
            dataset.read(1, out_shape=shape, resampling=Resampling.average),
            dataset.nodata,
        )


def _pair_files(dolphin_dir: Path, suffix: str) -> list[Path]:
    """Find final Dolphin pair products, ordered by their acquisition dates."""
    paths = sorted((dolphin_dir / "interferograms").glob(f"*.int{suffix}.tif"))
    if not paths:
        raise FileNotFoundError(
            f"no final *{suffix}.tif products in {dolphin_dir / 'interferograms'}"
        )
    if any(PAIR.search(path.name) is None for path in paths):
        raise ValueError("could not parse dates from one or more Dolphin pair products")
    return paths


def _closure_files(dolphin_dir: Path) -> list[Path]:
    """Find final closure phases; fall back to phase-linking intermediates."""
    paths = sorted((dolphin_dir / "interferograms").glob("closure_phase_*.tif"))
    paths = [path for path in paths if CLOSURE.search(path.name)]
    if not paths:
        pattern = "phase_linking/linked_phase/*/closure_phases/*.tif"
        paths = sorted(dolphin_dir.glob(pattern))
    if not paths:
        raise FileNotFoundError("no closure_phase_<date>_<date>_<date>.tif found")
    return paths


def _pair_label(path: Path) -> str:
    """Make a short panel label and temporal separation from a pair filename."""
    match = PAIR.search(path.name)
    assert match is not None
    start, end = (datetime.strptime(value, "%Y%m%d") for value in match.groups())
    return f"{start:%b %d}–{end:%b %d}  ({(end - start).days} d)"


def _closure_label(path: Path) -> str:
    """Make a short label for an adjacent three-date closure phase."""
    match = CLOSURE.search(path.name)
    assert match is not None
    return "–".join(
        datetime.strptime(value, "%Y%m%d").strftime("%b %d")
        for value in match.groups()
    )


def _circular_center(phase: np.ndarray, valid: np.ndarray) -> tuple[np.ndarray, float]:
    """Shift a wrapped phase raster to zero circular mean."""
    center = np.angle(np.mean(np.exp(1j * phase[valid]))) if valid.any() else np.nan
    return np.angle(np.exp(1j * (phase - center))), float(center)


def _valid_data(data: np.ndarray, nodata: float) -> np.ndarray:
    """Mask the explicit Dolphin nodata value, including complex zero."""
    valid = np.isfinite(data)
    if nodata is not None:
        valid &= np.abs(data - nodata) > 1.0e-12
    return valid


def _phase_cmap():
    """Use sarlet's oil-slick colors, leaving nodata transparent."""
    cmap = OIL_SLICK_CMAP.copy()
    cmap.set_bad((1, 1, 1, 0))
    return cmap


def _coherence_cmap():
    """Use a clear sequential coherence map with transparent nodata."""
    cmap = colormaps["viridis"].copy()
    cmap.set_bad((1, 1, 1, 0))
    return cmap


def _panel_grid(count: int, *, max_columns: int = 5) -> tuple[int, int]:
    """Use a compact row-major small-multiples layout."""
    columns = min(count, max_columns)
    return math.ceil(count / columns), columns


def _format_panels(axes: np.ndarray, count: int) -> list:
    """Flatten axes and hide unused cells without changing the grid geometry."""
    flat = list(np.asarray(axes).flat)
    for axis in flat[count:]:
        axis.axis("off")
    return flat[:count]


def _save_pair_panel(
    paths: list[Path],
    output: Path,
    *,
    coherence: bool,
) -> Path:
    """Render all final interferograms or coherences with one shared scale."""
    rows, columns = _panel_grid(len(paths))
    fig, axes = uplt.subplots(
        nrows=rows,
        ncols=columns,
        refwidth="2.2in",
        share=False,
        abc="a)",
    )
    last_image = None
    for axis, path in zip(_format_panels(axes, len(paths)), paths, strict=True):
        data, nodata = _read_display(path)
        valid = _valid_data(data, nodata)
        if coherence:
            image = axis.imshow(
                np.ma.masked_where(~valid, data),
                cmap=_coherence_cmap(),
                vmin=0,
                vmax=1,
            )
            subtitle = _pair_label(path)
        else:
            phase, center = _circular_center(np.angle(data), valid)
            image = axis.imshow(
                np.ma.masked_where(~valid, phase),
                cmap=_phase_cmap(),
                vmin=-np.pi,
                vmax=np.pi,
            )
            subtitle = f"{_pair_label(path)}  μ={center:+.2f}"
        axis.format(title=subtitle, xticks=[], yticks=[])
        last_image = image

    assert last_image is not None
    label = "coherence" if coherence else "phase relative to circular mean (rad)"
    fig.colorbar(last_image, loc="r", label=label)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    uplt.close(fig)
    return output


def _save_closure_panel(paths: list[Path], output: Path) -> Path:
    """Render adjacent-triplet closure phase with a physically common zero."""
    rows, columns = _panel_grid(len(paths))
    fig, axes = uplt.subplots(
        nrows=rows,
        ncols=columns,
        refwidth="2.2in",
        share=False,
        abc="a)",
    )
    last_image = None
    for axis, path in zip(_format_panels(axes, len(paths)), paths, strict=True):
        phase, nodata = _read_display(path)
        valid = _valid_data(phase, nodata)
        image = axis.imshow(
            np.ma.masked_where(~valid, phase),
            cmap=_phase_cmap(),
            vmin=-np.pi,
            vmax=np.pi,
        )
        axis.format(title=_closure_label(path), xticks=[], yticks=[])
        last_image = image

    assert last_image is not None
    fig.colorbar(last_image, loc="r", label="closure phase (rad)")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    uplt.close(fig)
    return output


def run(dolphin_dir: Path, output_dir: Path | None = None) -> list[Path]:
    """Render the three requested Dolphin diagnostic panels."""
    if not dolphin_dir.is_dir():
        raise NotADirectoryError(dolphin_dir)
    output_dir = output_dir or dolphin_dir.parent / "diagnostics"
    interferograms = _pair_files(dolphin_dir, "")
    coherences = _pair_files(dolphin_dir, ".cor")
    if {path.stem.removesuffix(".int") for path in interferograms} != {
        path.stem.removesuffix(".int.cor") for path in coherences
    }:
        raise ValueError("interferogram and coherence pair sets do not match")
    closures = _closure_files(dolphin_dir)
    return [
        _save_pair_panel(
            interferograms, output_dir / "dolphin_interferograms.png", coherence=False
        ),
        _save_pair_panel(
            coherences, output_dir / "dolphin_coherences.png", coherence=True
        ),
        _save_closure_panel(closures, output_dir / "dolphin_closure_phases.png"),
    ]


def main() -> None:
    """Parse a completed Dolphin run directory and write diagnostic panels."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dolphin_dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    for output in run(**vars(parser.parse_args())):
        print(output)


if __name__ == "__main__":
    main()
