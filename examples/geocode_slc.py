#!/usr/bin/env python
"""Geocode one BIOMASS L1a SCS scene into a GSLC GeoTIFF.

This is the minimal end-to-end workflow: parse the product with
``biomass-reader``, then hand the isce3 objects to
``isce3.geocode.geocode_slc`` -- the same primitive COMPASS uses for
Sentinel-1 and the NISAR GSLC SAS uses for NISAR.

Usage
-----
    python examples/geocode_slc.py PRODUCT_DIR DEM.tif OUT_GSLC.tif \
        --polarization HH --spacing 50

Notes
-----
This does *not* yet apply the ionosphere correction (dominant at P-band); that
is a follow-on step that consumes ``biomass_reader.ionosphere``.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def geocode_scene(
    product_dir: Path,
    dem_file: Path,
    out_file: Path,
    polarization: str = "HH",
    spacing: float = 50.0,
) -> Path:
    """Geocode one polarization of a BIOMASS SCS scene to a GSLC GeoTIFF."""
    import isce3

    from biomass_reader import BiomassSlc
    from biomass_reader.geocode import geocode_slc, make_shared_geogrid

    slc = BiomassSlc.from_dir(product_dir, polarization)
    dem_raster = isce3.io.Raster(str(dem_file))
    epsg = dem_raster.get_epsg()
    geogrid = make_shared_geogrid([slc], epsg, spacing)
    geocode_slc(slc, dem_file, geogrid, out_file)
    print(f"wrote {out_file}  {(geogrid.length, geogrid.width)}  EPSG:{epsg}")
    return out_file


def main() -> None:  # noqa: D103
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("product_dir", type=Path)
    p.add_argument("dem_file", type=Path)
    p.add_argument("out_file", type=Path)
    p.add_argument("--polarization", default="HH")
    p.add_argument("--spacing", type=float, default=50.0)
    args = p.parse_args()
    geocode_scene(
        args.product_dir,
        args.dem_file,
        args.out_file,
        args.polarization,
        args.spacing,
    )


if __name__ == "__main__":
    main()
