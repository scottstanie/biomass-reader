#!/usr/bin/env python
"""Explore temporal and polarimetric phase observables in BIOMASS GSLCs.

This is a research diagnostic, not part of the production reader. All products
must be on the same geogrid and follow ``biomass_<POL>_<YYYYMMDD>.tif`` naming.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import rasterio


def _looks_sum(array: np.ndarray, looks: int) -> np.ndarray:
    rows = array.shape[0] // looks * looks
    cols = array.shape[1] // looks * looks
    reshaped = array[:rows, :cols].reshape(
        rows // looks, looks, cols // looks, looks
    )
    return np.nansum(reshaped, axis=(1, 3))


def _cross(first_path: Path, second_path: Path, looks: int):
    with rasterio.open(first_path) as dataset:
        first = dataset.read(1)
    with rasterio.open(second_path) as dataset:
        second = dataset.read(1)
    valid = np.isfinite(first) & np.isfinite(second)
    cross = _looks_sum(np.where(valid, first * second.conj(), np.nan), looks)
    first_power = _looks_sum(np.where(valid, np.abs(first) ** 2, np.nan), looks)
    second_power = _looks_sum(np.where(valid, np.abs(second) ** 2, np.nan), looks)
    del first, second
    with np.errstate(divide="ignore", invalid="ignore"):
        coherence = np.abs(cross) / np.sqrt(first_power * second_power)
    coherence[(first_power == 0) | (second_power == 0)] = np.nan
    return cross, np.clip(coherence, 0, 1)


def _phase_stats(phase: np.ndarray, mask: np.ndarray) -> dict:
    values = phase[mask]
    resultant = np.abs(np.mean(np.exp(1j * values))) if values.size else np.nan
    return {
        "pixel_count": int(values.size),
        "circular_mean_rad": float(np.angle(np.mean(np.exp(1j * values)))),
        "circular_resultant": float(resultant),
        "median_absolute_phase_rad": float(np.median(np.abs(values))),
    }


def run(gslc_dirs: list[Path], output_dir: Path, looks: int = 6) -> dict:
    """Compute polarimetric phase diagnostics and write JSON/PNG outputs."""
    import matplotlib.pyplot as plt

    paths: dict[str, dict[str, Path]] = {}
    for directory in gslc_dirs:
        matches = sorted(directory.glob("biomass_??_????????.tif"))
        polarizations = {path.stem.split("_")[1] for path in matches}
        if len(polarizations) != 1:
            raise ValueError(f"expected one polarization in {directory}")
        polarization = polarizations.pop()
        if polarization in paths:
            raise ValueError(f"duplicate directory for {polarization}")
        paths[polarization] = {path.stem[-8:]: path for path in matches}
    required = {"HH", "HV", "VH", "VV"}
    if set(paths) != required:
        raise ValueError(
            f"expected directories for {sorted(required)}, got {sorted(paths)}"
        )
    dates = sorted(set.intersection(*(set(items) for items in paths.values())))
    if len(dates) < 3:
        raise ValueError("need at least three common dates for HH, HV, VH, and VV")

    temporal: dict[str, dict[tuple[int, int], tuple[np.ndarray, np.ndarray]]] = {}
    result: dict = {"dates": dates, "looks": looks, "temporal": {}, "polarimetric": {}}
    for polarization in paths:
        temporal[polarization] = {}
        result["temporal"][polarization] = {}
        for first in range(len(dates) - 1):
            for second in range(first + 1, len(dates)):
                cross, coherence = _cross(
                    paths[polarization][dates[first]],
                    paths[polarization][dates[second]],
                    looks,
                )
                temporal[polarization][(first, second)] = (cross, coherence)
                finite = np.isfinite(coherence)
                result["temporal"][polarization][
                    f"{dates[first]}_{dates[second]}"
                ] = {
                    "median_coherence": float(np.nanmedian(coherence)),
                    "fraction_coherence_ge_0.5": float(
                        np.mean(coherence[finite] >= 0.5)
                    ),
                }

    polar_cross = {}
    for first_pol, second_pol in (("HH", "VV"), ("HV", "VH")):
        key = f"{first_pol}_{second_pol}"
        result["polarimetric"][key] = {}
        for date in dates:
            cross, coherence = _cross(
                paths[first_pol][date], paths[second_pol][date], looks
            )
            polar_cross[(key, date)] = (cross, coherence)
            high = np.isfinite(coherence) & (coherence >= 0.5)
            result["polarimetric"][key][date] = {
                "median_coherence": float(np.nanmedian(coherence)),
                **_phase_stats(np.angle(cross), high),
            }

    differentials = {}
    differential_masks = {}
    result["temporal_differential_phase"] = {}
    for first in range(len(dates) - 1):
        for second in range(first + 1, len(dates)):
            key = f"{dates[first]}_{dates[second]}"
            hh, hh_coh = temporal["HH"][(first, second)]
            vv, vv_coh = temporal["VV"][(first, second)]
            differential = np.angle(hh * vv.conj())
            high = np.isfinite(differential) & (hh_coh >= 0.5) & (vv_coh >= 0.5)
            differentials[(first, second)] = differential
            differential_masks[(first, second)] = (
                np.isfinite(differential) & (hh_coh >= 0.3) & (vv_coh >= 0.3)
            )
            result["temporal_differential_phase"][key] = _phase_stats(
                differential, high
            )

    result["closure"] = {}
    closures = {}
    closure_browse_masks = {}
    for polarization in paths:
        ifg01, coh01 = temporal[polarization][(0, 1)]
        ifg12, coh12 = temporal[polarization][(1, 2)]
        ifg02, coh02 = temporal[polarization][(0, 2)]
        closure = np.angle(ifg01 * ifg12 * ifg02.conj())
        high = (
            np.isfinite(closure)
            & (coh01 >= 0.5)
            & (coh12 >= 0.5)
            & (coh02 >= 0.5)
        )
        closures[polarization] = closure
        closure_browse_masks[polarization] = (
            np.isfinite(closure)
            & (coh01 >= 0.3)
            & (coh12 >= 0.3)
            & (coh02 >= 0.3)
        )
        result["closure"][polarization] = _phase_stats(closure, high)

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "polarimetric_phase.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    figure, axes = plt.subplots(2, 4, figsize=(18, 9), constrained_layout=True)
    hhvv = polar_cross[("HH_VV", dates[0])]
    hvvh = polar_cross[("HV_VH", dates[0])]
    differential_closure = np.angle(
        np.exp(1j * (closures["HH"] - closures["VV"]))
    )
    panels = [
        (np.where(hhvv[1] >= 0.3, np.angle(hhvv[0]), np.nan), "HH-VV phase"),
        (hhvv[1], "HH-VV coherence"),
        (np.where(hvvh[1] >= 0.3, np.angle(hvvh[0]), np.nan), "HV-VH phase"),
        (hvvh[1], "HV-VH coherence"),
        (
            np.where(differential_masks[(0, 1)], differentials[(0, 1)], np.nan),
            "HH-VV temporal phase difference",
        ),
        (
            np.where(closure_browse_masks["HH"], closures["HH"], np.nan),
            "HH closure",
        ),
        (
            np.where(closure_browse_masks["HV"], closures["HV"], np.nan),
            "HV closure",
        ),
        (
            np.where(
                closure_browse_masks["HH"] & closure_browse_masks["VV"],
                differential_closure,
                np.nan,
            ),
            "HH-VV differential closure",
        ),
    ]
    for axis, (data, title) in zip(axes.flat, panels, strict=True):
        if "coherence" in title.lower():
            axis.imshow(data, cmap="magma", vmin=0, vmax=1)
        else:
            finite = np.abs(data[np.isfinite(data)])
            limit = min(np.pi, max(0.05, float(np.percentile(finite, 98))))
            axis.imshow(data, cmap="twilight", vmin=-limit, vmax=limit)
        axis.set_title(title)
        axis.set_axis_off()
    figure.savefig(output_dir / "polarimetric_phase.png", dpi=120)
    plt.close(figure)
    return result


def main() -> None:
    """Run the exploratory polarimetric CLI."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gslc_dirs", type=Path, nargs=4)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--looks", type=int, default=6)
    args = parser.parse_args()
    print(json.dumps(run(args.gslc_dirs, args.output_dir, args.looks), indent=2))


if __name__ == "__main__":
    main()
