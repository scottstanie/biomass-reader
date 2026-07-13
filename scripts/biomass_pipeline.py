#!/usr/bin/env python
"""Build a BIOMASS GSLC stack and optionally run Dolphin/Whirlwind.

Data acquisition is intentionally independent (see ``download_maap.py``). This
command is restartable and operates entirely on already extracted products.

Example
-------
python scripts/biomass_pipeline.py \
    --products /data/BIO_S1_SCS__1S_* \
    --dem /data/dem_utm.tif --work /data/validation \
    --polarization HH --spacing 30 --run-dolphin
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import time
from pathlib import Path

import rasterio

from biomass_reader import BiomassSlc
from biomass_reader.geocode import (
    geocode_slc,
    make_shared_geogrid,
    write_stack_provenance,
)


def _make_dolphin_config(
    gslcs: list[Path], work_dir: Path, wavelength: float, path: Path
) -> Path:
    """Create a validated Dolphin configuration using Whirlwind."""
    try:
        from dolphin.workflows.config import DisplacementWorkflow
    except ImportError as error:
        raise RuntimeError("Dolphin is required for --run-dolphin") from error

    config = DisplacementWorkflow.model_validate(
        {
            "cslc_file_list": gslcs,
            "work_directory": work_dir,
            "input_options": {"wavelength": wavelength},
            "output_options": {"strides": {"x": 6, "y": 6}},
            "phase_linking": {
                "ministack_size": len(gslcs),
                "half_window": {"x": 3, "y": 3},
                "write_closure_phase": len(gslcs) >= 3,
            },
            "interferogram_network": {"max_bandwidth": min(3, len(gslcs) - 1)},
            "unwrap_options": {
                "run_unwrap": True,
                "unwrap_method": "whirlwind",
                "n_parallel_jobs": 1,
            },
            "timeseries_options": {
                "run_inversion": len(gslcs) >= 3,
                "run_velocity": False,
            },
        }
    )
    config.to_yaml(path, with_comments=False)
    return path


def _stack_key(product_id: str) -> tuple[str, str]:
    match = re.search(r"_(T\d{3})_(F\d{3})_", product_id)
    if match is None:
        raise ValueError(f"cannot determine track/frame from {product_id!r}")
    return match.group(1), match.group(2)


def _matches_output(path: Path, geogrid, slc: BiomassSlc) -> bool:
    if not path.exists():
        return False
    with rasterio.open(path) as dataset:
        expected_transform = (
            geogrid.spacing_x,
            0.0,
            geogrid.start_x,
            0.0,
            geogrid.spacing_y,
            geogrid.start_y,
        )
        return (
            dataset.shape == (geogrid.length, geogrid.width)
            and dataset.crs is not None
            and dataset.crs.to_epsg() == geogrid.epsg
            and tuple(dataset.transform)[:6] == expected_transform
            and dataset.nodata is not None
            and dataset.tags().get("BIOMASS_PRODUCT_ID") == slc.product_id
            and dataset.tags().get("POLARIZATION") == slc.polarization
        )


def run(args: argparse.Namespace) -> list[Path]:
    """Run the local-product GSLC workflow."""
    products = sorted({path.resolve() for path in args.products})
    if len(products) < 2:
        raise ValueError("provide at least two extracted repeat-pass products")
    slcs = sorted(
        (BiomassSlc.from_dir(path, args.polarization) for path in products),
        key=lambda slc: slc.sensing_start,
    )
    if len({slc.swath for slc in slcs}) != 1:
        raise ValueError("all products must use the same swath")
    if len({_stack_key(slc.product_id) for slc in slcs}) != 1:
        raise ValueError("all products must use the same track and frame")

    args.work.mkdir(parents=True, exist_ok=True)
    gslc_dir = args.work / "gslc"
    gslc_dir.mkdir(exist_ok=True)
    import isce3

    dem = isce3.io.Raster(str(args.dem))
    geogrid = make_shared_geogrid(
        slcs, dem.get_epsg(), args.spacing, extent=args.extent
    )
    print(
        f"shared {args.extent} geogrid: {geogrid.length} x {geogrid.width}, "
        f"EPSG:{geogrid.epsg}, {args.spacing:g} m"
    )

    outputs: list[Path] = []
    for slc in slcs:
        date = slc.sensing_start.strftime("%Y%m%d")
        output = gslc_dir / f"biomass_{slc.polarization}_{date}.tif"
        if args.resume and _matches_output(output, geogrid, slc):
            print(f"reusing {output}")
        else:
            started = time.monotonic()
            geocode_slc(slc, args.dem, geogrid, output, flatten=True)
            print(f"wrote {output} in {time.monotonic() - started:.1f} s")
        outputs.append(output)

    write_stack_provenance(
        args.work / "stack.json",
        slcs,
        args.dem,
        geogrid,
        args.extent,
        flatten=True,
    )
    if args.run_dolphin:
        config = _make_dolphin_config(
            outputs,
            args.work / "dolphin",
            slcs[0].wavelength,
            args.work / "dolphin.yaml",
        )
        subprocess.run(
            [sys.executable, "-m", "dolphin", "run", str(config)], check=True
        )
    return outputs


def main() -> None:
    """Parse command-line arguments and run the workflow."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--products", type=Path, nargs="+", required=True)
    parser.add_argument("--dem", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--polarization", default="HH")
    parser.add_argument("--spacing", type=float, default=30.0)
    parser.add_argument("--extent", choices=("union", "intersection"), default="union")
    parser.add_argument("--resume", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--run-dolphin", action="store_true")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
