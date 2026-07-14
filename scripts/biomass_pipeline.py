#!/usr/bin/env python
"""Build a BIOMASS GSLC stack and optionally run Dolphin/Whirlwind.

Data acquisition is intentionally independent (see ``download_maap.py``). This
command is restartable and operates entirely on already extracted products.

Example
-------
python scripts/biomass_pipeline.py \
    --products /data/BIO_S1_SCS__1S_* \
    --dem /data/dem_utm.tif --work /data/validation \
    --polarization HH --posting 30 30 --run-dolphin
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import time
from dataclasses import asdict
from pathlib import Path

import rasterio

from biomass_reader import BiomassSlc, parse_correction_status
from biomass_reader.geocode import (
    geocode_slc,
    make_shared_geogrid,
    native_posting,
    write_stack_provenance,
)
from biomass_reader.ionosphere import ionosphere_lut_summary


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
            "output_options": {"strides": {"x": 3, "y": 3}},
            "phase_linking": {
                "ministack_size": len(gslcs),
                "half_window": {"x": 7, "y": 7},
                "write_closure_phase": len(gslcs) >= 3,
                "nearest_n_coherence": 3,
            },
            "interferogram_network": {"max_bandwidth": 3},
            "unwrap_options": {
                "run_unwrap": True,
                "unwrap_method": "whirlwind",
                "n_parallel_jobs": 1,
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


def _matches_output(
    path: Path, geogrid, slc: BiomassSlc, ionosphere_record: dict
) -> bool:
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
            and dataset.tags().get("IONOSPHERE_POLICY") == ionosphere_record["policy"]
        )


def _ionosphere_provenance(slcs: list[BiomassSlc], policy: str) -> list[dict]:
    """Validate the L1 correction state and record the shipped LUT schema."""
    records: list[dict] = []
    required_layers = {"phaseScreen", "rangeShifts"}
    for slc in slcs:
        status = parse_correction_status(slc.product.annotation_xml)
        lut = ionosphere_lut_summary(slc.product.lut_nc)
        available_layers = set(lut["variables"])
        missing_layers = sorted(required_layers - available_layers)
        l1_applied = status.ionospheric_phase_screen and status.ionospheric_group_delay
        if policy == "require-l1-applied" and (not l1_applied or missing_layers):
            raise ValueError(
                f"{slc.product_id} is not safe for the default GSLC policy: "
                f"phase-screen/group-delay applied={l1_applied}, "
                f"missing ionosphere LUT layers={missing_layers}. "
                "Use --ionosphere-policy allow-unapplied only when you will "
                "apply and validate a residual correction separately."
            )
        records.append(
            {
                "product": slc.product_id,
                "policy": (
                    "l1_phase_and_group_delay_already_applied"
                    if l1_applied
                    else "unapplied_l1_ionosphere_allowed"
                ),
                "correction_status": asdict(status),
                "lut": lut,
            }
        )
    return records


def _write_ionosphere_tags(path: Path, record: dict) -> None:
    """Attach the correction policy to an already-written GSLC."""
    status = record["correction_status"]
    with rasterio.open(path, "r+") as dataset:
        dataset.update_tags(
            IONOSPHERE_POLICY=record["policy"],
            IONOSPHERIC_PHASE_SCREEN_CORRECTED=str(
                status["ionospheric_phase_screen"]
            ).lower(),
            IONOSPHERIC_GROUP_DELAY_CORRECTED=str(
                status["ionospheric_group_delay"]
            ).lower(),
            IONOSPHERE_LUT_VARIABLES=",".join(record["lut"]["variables"]),
        )


def _resolve_posting(
    args: argparse.Namespace, slcs: list[BiomassSlc], epsg: int
) -> tuple[tuple[float, float], str]:
    """Choose the requested explicit, legacy, or native UTM posting."""
    if args.native:
        return native_posting(slcs, epsg), "native"
    if args.posting is not None:
        return tuple(args.posting), "explicit"
    if args.spacing is not None:
        return (args.spacing, args.spacing), "legacy-spacing"
    return (30.0, 30.0), "default"


def _plot_dolphin_diagnostics(dolphin_dir: Path) -> None:
    """Make the standard diagnostic panels after a successful Dolphin run."""
    try:
        from plot_dolphin_diagnostics import run as plot
    except ImportError as error:
        raise RuntimeError(
            "default Dolphin plotting requires biomass-reader[plot]; "
            "install it or pass --no-plot"
        ) from error
    for path in plot(dolphin_dir):
        print(f"wrote {path}")


def run(args: argparse.Namespace, *, min_products: int = 2) -> list[Path]:
    """Run the local-product GSLC workflow."""
    products = sorted({path.resolve() for path in args.products})
    if len(products) < min_products:
        raise ValueError(f"provide at least {min_products} extracted product(s)")
    slcs = sorted(
        (BiomassSlc.from_dir(path, args.polarization) for path in products),
        key=lambda slc: slc.sensing_start,
    )
    if len({slc.swath for slc in slcs}) != 1:
        raise ValueError("all products must use the same swath")
    if len({_stack_key(slc.product_id) for slc in slcs}) != 1:
        raise ValueError("all products must use the same track and frame")
    ionosphere = _ionosphere_provenance(slcs, args.ionosphere_policy)

    args.work.mkdir(parents=True, exist_ok=True)
    gslc_dir = args.work / "gslc"
    gslc_dir.mkdir(exist_ok=True)
    import isce3

    dem = isce3.io.Raster(str(args.dem))
    posting, posting_mode = _resolve_posting(args, slcs, dem.get_epsg())
    geogrid = make_shared_geogrid(
        slcs, dem.get_epsg(), posting, extent=args.extent
    )
    print(
        f"shared {args.extent} geogrid: {geogrid.length} x {geogrid.width}, "
        f"EPSG:{geogrid.epsg}, {posting[0]:g} x {posting[1]:g} m "
        f"({posting_mode})"
    )

    outputs: list[Path] = []
    for slc, ionosphere_record in zip(slcs, ionosphere, strict=True):
        date = slc.sensing_start.strftime("%Y%m%d")
        output = gslc_dir / f"biomass_{slc.polarization}_{date}.tif"
        if args.resume and _matches_output(output, geogrid, slc, ionosphere_record):
            print(f"reusing {output}")
        else:
            started = time.monotonic()
            geocode_slc(slc, args.dem, geogrid, output, flatten=True)
            _write_ionosphere_tags(output, ionosphere_record)
            print(f"wrote {output} in {time.monotonic() - started:.1f} s")
        outputs.append(output)

    write_stack_provenance(
        args.work / "stack.json",
        slcs,
        args.dem,
        geogrid,
        args.extent,
        flatten=True,
        ionosphere=ionosphere,
        posting_mode=posting_mode,
    )
    if getattr(args, "run_dolphin", False):
        config = _make_dolphin_config(
            outputs,
            args.work / "dolphin",
            slcs[0].wavelength,
            args.work / "dolphin.yaml",
        )
        subprocess.run(
            [sys.executable, "-m", "dolphin", "run", str(config)], check=True
        )
        if args.plot:
            _plot_dolphin_diagnostics(args.work / "dolphin")
    return outputs


def make_parser(
    *, include_dolphin: bool = True, description: str | None = None
) -> argparse.ArgumentParser:
    """Build the shared GSLC-workflow parser."""
    parser = argparse.ArgumentParser(description=description or __doc__)
    parser.add_argument("--products", type=Path, nargs="+", required=True)
    parser.add_argument("--dem", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--polarization", default="HH")
    posting = parser.add_mutually_exclusive_group()
    posting.add_argument(
        "--posting",
        type=float,
        nargs=2,
        metavar=("X", "Y"),
        help="output easting and northing posting in metres",
    )
    posting.add_argument(
        "--native",
        action="store_true",
        help="choose conservative 10 / 2**n UTM postings from source GCP sampling",
    )
    posting.add_argument(
        "--spacing",
        type=float,
        help="deprecated square-posting alias; use --posting X Y",
    )
    parser.add_argument("--extent", choices=("union", "intersection"), default="union")
    parser.add_argument("--resume", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument(
        "--ionosphere-policy",
        choices=("require-l1-applied", "allow-unapplied"),
        default="require-l1-applied",
        help=(
            "default verifies that L1 already applied phase-screen and group-delay "
            "corrections; it never applies the shipped LUT a second time"
        ),
    )
    if include_dolphin:
        parser.add_argument("--run-dolphin", action="store_true")
        parser.add_argument(
            "--plot",
            action=argparse.BooleanOptionalAction,
            default=True,
            help="write Dolphin diagnostic panels after a successful run (default)",
        )
    return parser


def main() -> None:
    """Parse command-line arguments and run the repeat-pass workflow."""
    parser = make_parser()
    run(parser.parse_args())


if __name__ == "__main__":
    main()
