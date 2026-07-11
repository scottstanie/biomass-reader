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

import numpy as np


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

    slc = BiomassSlc.from_dir(product_dir, polarization)
    radar_grid = slc.radar_grid
    orbit = slc.orbit
    native_doppler = slc.doppler

    dem_raster = isce3.io.Raster(str(dem_file))
    epsg = dem_raster.get_epsg()
    proj = isce3.core.make_projection(epsg)
    ellipsoid = proj.ellipsoid

    # Geogrid covering the scene footprint at the requested ground spacing.
    geogrid = isce3.product.bbox_to_geogrid(
        radar_grid, orbit, native_doppler, spacing, -spacing, epsg
    )

    # Read the full complex scene (radar coordinates).
    rdr_data = slc.read_complex()
    geo_data = np.zeros((geogrid.length, geogrid.width), dtype=np.complex64)

    isce3.geocode.geocode_slc(
        geo_data_blocks=geo_data,
        rdr_data_blocks=rdr_data,
        dem_raster=dem_raster,
        radargrid=radar_grid,
        geogrid=geogrid,
        orbit=orbit,
        native_doppler=native_doppler,
        image_grid_doppler=isce3.core.LUT2d(),
        ellipsoid=ellipsoid,
        threshold_geo2rdr=1.0e-8,
        num_iter_geo2rdr=25,
        flatten=True,
    )

    _write_geotiff(geo_data, geogrid, epsg, out_file)
    print(f"wrote {out_file}  {geo_data.shape}  EPSG:{epsg}")
    return out_file


def _write_geotiff(data, geogrid, epsg, out_file: Path) -> None:
    from osgeo import gdal, osr

    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(
        str(out_file),
        geogrid.width,
        geogrid.length,
        1,
        gdal.GDT_CFloat32,
        options=["COMPRESS=DEFLATE", "TILED=YES"],
    )
    ds.SetGeoTransform(
        (geogrid.start_x, geogrid.spacing_x, 0, geogrid.start_y, 0, geogrid.spacing_y)
    )
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(data)
    ds.FlushCache()


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
