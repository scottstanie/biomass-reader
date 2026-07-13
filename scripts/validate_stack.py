#!/usr/bin/env python
"""Quantify alignment, coherence, and closure phase for a GSLC stack."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import rasterio


def _looks_sum(array: np.ndarray, looks: int) -> np.ndarray:
    rows = array.shape[0] // looks * looks
    cols = array.shape[1] // looks * looks
    return np.nansum(
        array[:rows, :cols].reshape(rows // looks, looks, cols // looks, looks),
        axis=(1, 3),
    )


def _pair(first: np.ndarray, second: np.ndarray, looks: int):
    valid = np.isfinite(first) & np.isfinite(second)
    cross = _looks_sum(np.where(valid, first * second.conj(), np.nan), looks)
    power_first = _looks_sum(np.where(valid, np.abs(first) ** 2, np.nan), looks)
    power_second = _looks_sum(np.where(valid, np.abs(second) ** 2, np.nan), looks)
    with np.errstate(divide="ignore", invalid="ignore"):
        coherence = np.abs(cross) / np.sqrt(power_first * power_second)
    coherence[(power_first == 0) | (power_second == 0)] = np.nan
    return cross, np.clip(coherence, 0, 1)


def run(
    paths: list[Path],
    output_dir: Path,
    looks: int = 6,
    dolphin_dir: Path | None = None,
) -> dict:
    """Validate a three-or-more acquisition stack and write JSON/PNG results."""
    import matplotlib.pyplot as plt

    paths = sorted(paths)
    if len(paths) < 3:
        raise ValueError("at least three GSLCs are required for closure validation")
    arrays = []
    metadata = []
    for path in paths:
        with rasterio.open(path) as dataset:
            arrays.append(dataset.read(1))
            metadata.append(
                {
                    "path": str(path.resolve()),
                    "shape": dataset.shape,
                    "crs": dataset.crs.to_string(),
                    "transform": tuple(dataset.transform)[:6],
                    "nodata": str(dataset.nodata),
                }
            )
    aligned = all(
        item["shape"] == metadata[0]["shape"]
        and item["crs"] == metadata[0]["crs"]
        and item["transform"] == metadata[0]["transform"]
        for item in metadata
    )
    if not aligned:
        raise ValueError("GSLC rasters are not on an identical grid")

    pairs = {}
    products = {}
    for first_idx in range(len(arrays) - 1):
        for second_idx in range(first_idx + 1, len(arrays)):
            key = f"{paths[first_idx].stem[-8:]}_{paths[second_idx].stem[-8:]}"
            cross, coherence = _pair(
                arrays[first_idx], arrays[second_idx], looks
            )
            pairs[(first_idx, second_idx)] = (cross, coherence)
            finite = np.isfinite(coherence)
            products[key] = {
                "valid_fraction": float(np.mean(finite)),
                "mean_coherence": float(np.nanmean(coherence)),
                "median_coherence": float(np.nanmedian(coherence)),
                "fraction_coherence_ge_0.5": float(np.mean(coherence[finite] >= 0.5)),
            }

    closure = None
    closure_browse = None
    closure_stats = None
    if len(arrays) >= 3:
        ifg01, coh01 = pairs[(0, 1)]
        ifg12, coh12 = pairs[(1, 2)]
        ifg02, coh02 = pairs[(0, 2)]
        closure = np.angle(ifg01 * ifg12 * ifg02.conj())
        min_coherence = np.minimum(np.minimum(coh01, coh12), coh02)
        closure_browse = np.where(min_coherence >= 0.3, closure, np.nan)
        high = np.isfinite(closure) & (min_coherence >= 0.5)
        closure_stats = {
            "high_coherence_pixel_count": int(np.count_nonzero(high)),
            "median_absolute_phase_rad": float(np.median(np.abs(closure[high]))),
            "p90_absolute_phase_rad": float(np.percentile(np.abs(closure[high]), 90)),
        }

    unwrap_stats = {}
    if dolphin_dir is not None:
        for path in sorted((dolphin_dir / "unwrapped").glob("*.unw.conncomp.tif")):
            with rasterio.open(path) as dataset:
                labels = dataset.read(1)
                nodata = dataset.nodata
            valid = (
                labels != nodata
                if nodata is not None
                else np.ones_like(labels, dtype=bool)
            )
            components = np.unique(labels[valid & (labels > 0)])
            unwrap_stats[path.name] = {
                "valid_fraction": float(np.mean(valid)),
                "labeled_fraction_of_valid": float(
                    np.mean(labels[valid] > 0) if np.any(valid) else 0.0
                ),
                "component_count": int(components.size),
            }

    result = {
        "aligned": aligned,
        "looks": looks,
        "rasters": metadata,
        "pairs": products,
        "closure": closure_stats,
        "whirlwind": unwrap_stats,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "validation.json").write_text(json.dumps(result, indent=2) + "\n")

    first_pair = pairs[(0, 1)]
    amplitude = 10 * np.log10(_looks_sum(np.abs(arrays[0]) ** 2, looks) + 1e-7)
    figure, axes = plt.subplots(1, 4, figsize=(18, 5), constrained_layout=True)
    axes[0].imshow(
        amplitude,
        cmap="gray",
        vmin=np.nanpercentile(amplitude, 2),
        vmax=np.nanpercentile(amplitude, 98),
    )
    axes[0].set_title("Reference amplitude")
    axes[1].imshow(np.angle(first_pair[0]), cmap="twilight", vmin=-np.pi, vmax=np.pi)
    axes[1].set_title("First-pair wrapped phase")
    image = axes[2].imshow(first_pair[1], cmap="magma", vmin=0, vmax=1)
    axes[2].set_title("First-pair coherence")
    figure.colorbar(image, ax=axes[2], shrink=0.7)
    axes[3].imshow(closure_browse, cmap="twilight", vmin=-np.pi, vmax=np.pi)
    axes[3].set_title("Three-date closure phase")
    for axis in axes:
        axis.set_axis_off()
    figure.savefig(output_dir / "validation.png", dpi=120)
    plt.close(figure)
    return result


def main() -> None:
    """Run the validation CLI."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gslcs", type=Path, nargs="+")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--looks", type=int, default=6)
    parser.add_argument("--dolphin-dir", type=Path)
    args = parser.parse_args()
    print(
        json.dumps(
            run(args.gslcs, args.output_dir, args.looks, args.dolphin_dir), indent=2
        )
    )


if __name__ == "__main__":
    main()
