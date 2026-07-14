"""Command-line interface for biomass-reader.

Examples
--------
Print scene metadata::

    biomass-reader info /path/to/BIO_S1_SCS__1S_..._DU1SS4

Write a complex SLC GeoTIFF for one polarization::

    biomass-reader to-geotiff /path/to/product --polarization HH --out hh.tif
"""

from __future__ import annotations

import argparse
from pathlib import Path

from .slc import BiomassSlc


def _print_info(args: argparse.Namespace) -> None:
    """Print radar-geometry metadata for one BIOMASS SCS product."""
    slc = BiomassSlc.from_dir(args.product_dir, args.polarization)
    m = slc._meta
    print(f"product         : {slc.identifier}")
    print(f"mission phase   : {m.mission_phase}   swath: {m.swath}")
    print(f"polarizations   : {', '.join(m.polarizations)}")
    print(f"look side       : {slc.look_side}")
    print(f"sensing start   : {slc.sensing_start.isoformat()}")
    print(f"shape (az, rg)  : {slc.shape}")
    freq_mhz = slc.center_frequency / 1e6
    print(f"wavelength [m]  : {slc.wavelength:.4f}  ({freq_mhz:.2f} MHz)")
    print(f"starting range  : {m.starting_range:.1f} m")
    print(f"range spacing   : {m.range_pixel_spacing:.3f} m")
    print(f"prf [Hz]        : {m.prf:.3f}")
    print(f"doppler ests    : {len(m.doppler_estimates)}")
    print(f"bounds (WSEN)   : {tuple(round(b, 3) for b in slc.bounds)}")


def _to_geotiff(args: argparse.Namespace) -> None:
    """Write a complex SLC GeoTIFF in radar coordinates."""
    slc = BiomassSlc.from_dir(args.product_dir, args.polarization)
    path = slc.to_geotiff(args.out)
    print(f"wrote {path} {slc.shape} ({slc.polarization})")


def main(argv: list[str] | None = None) -> None:
    """Entry point for the ``biomass-reader`` command."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    info = commands.add_parser("info", help="print product radar geometry")
    info.add_argument("product_dir", type=Path)
    info.add_argument("--polarization", default="HH")
    info.set_defaults(handler=_print_info)

    geotiff = commands.add_parser("to-geotiff", help="write a complex GeoTIFF")
    geotiff.add_argument("product_dir", type=Path)
    geotiff.add_argument("--out", type=Path, required=True)
    geotiff.add_argument("--polarization", default="HH")
    geotiff.set_defaults(handler=_to_geotiff)

    args = parser.parse_args(argv)
    args.handler(args)


if __name__ == "__main__":
    main()
