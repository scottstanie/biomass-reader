"""Command-line interface for biomass-reader.

Examples
--------
Print scene metadata::

    biomass-reader info /path/to/BIO_S1_SCS__1S_..._DU1SS4

Write a complex SLC GeoTIFF for one polarization::

    biomass-reader to-geotiff /path/to/product --polarization HH --out hh.tif
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import tyro

from .slc import BiomassSlc


@dataclass
class Info:
    """Print radar-geometry metadata for a BIOMASS SCS product."""

    product_dir: Path
    """Path to the extracted SCS product directory."""
    polarization: str = "HH"
    """Polarization to report (HH, HV, VH, VV)."""

    def run(self) -> None:  # noqa: D102
        slc = BiomassSlc.from_dir(self.product_dir, self.polarization)
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


@dataclass
class ToGeotiff:
    """Write a complex SLC GeoTIFF (radar coordinates) for one polarization."""

    product_dir: Path
    """Path to the extracted SCS product directory."""
    out: Path
    """Output GeoTIFF path."""
    polarization: str = "HH"
    """Polarization to extract (HH, HV, VH, VV)."""

    def run(self) -> None:  # noqa: D102
        slc = BiomassSlc.from_dir(self.product_dir, self.polarization)
        path = slc.to_geotiff(self.out)
        print(f"wrote {path} {slc.shape} ({slc.polarization})")


def main() -> None:
    """Entry point for the ``biomass-reader`` command."""
    command = tyro.cli(Info | ToGeotiff)
    command.run()


if __name__ == "__main__":
    main()
