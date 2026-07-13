"""Phase-preserving geocoding helpers for BIOMASS SLC stacks."""

from __future__ import annotations

import json
import math
from collections.abc import Sequence
from pathlib import Path
from typing import Literal

import numpy as np
import rasterio
from affine import Affine

from .slc import BiomassSlc


def make_shared_geogrid(
    slcs: Sequence[BiomassSlc],
    epsg: int,
    spacing: float,
    extent: Literal["union", "intersection"] = "union",
):
    """Create one snapped geogrid covering a stack of BIOMASS acquisitions."""
    import isce3

    if not slcs:
        raise ValueError("at least one SLC is required")
    if spacing <= 0:
        raise ValueError("spacing must be positive")

    grids = [
        isce3.product.bbox_to_geogrid(
            slc.radar_grid,
            slc.orbit,
            slc.doppler,
            spacing,
            -spacing,
            epsg,
        )
        for slc in slcs
    ]
    xmins = [grid.start_x for grid in grids]
    xmaxs = [grid.start_x + grid.width * spacing for grid in grids]
    ymaxs = [grid.start_y for grid in grids]
    ymins = [grid.start_y - grid.length * spacing for grid in grids]

    if extent == "union":
        xmin, xmax = min(xmins), max(xmaxs)
        ymin, ymax = min(ymins), max(ymaxs)
    elif extent == "intersection":
        xmin, xmax = max(xmins), min(xmaxs)
        ymin, ymax = max(ymins), min(ymaxs)
        if xmin >= xmax or ymin >= ymax:
            raise ValueError("the SLC geogrids do not overlap")
    else:
        raise ValueError(f"unknown extent policy: {extent!r}")

    # Snap to a projection-fixed lattice so repeated runs and subset stacks
    # produce pixel-aligned grids.
    xmin = math.floor(xmin / spacing) * spacing
    xmax = math.ceil(xmax / spacing) * spacing
    ymin = math.floor(ymin / spacing) * spacing
    ymax = math.ceil(ymax / spacing) * spacing
    width = int(round((xmax - xmin) / spacing))
    length = int(round((ymax - ymin) / spacing))
    return isce3.product.GeoGridParameters(
        xmin, ymax, spacing, -spacing, width, length, epsg
    )


def geocode_slc(
    slc: BiomassSlc,
    dem_file: str | Path,
    geogrid,
    output_file: str | Path,
    *,
    flatten: bool = True,
) -> Path:
    """Geocode one SLC to a supplied shared geogrid and retain invalid pixels."""
    import isce3

    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    dem = isce3.io.Raster(str(dem_file))
    ellipsoid = isce3.core.make_projection(geogrid.epsg).ellipsoid
    invalid = np.complex64(np.nan + 1j * np.nan)
    output = np.full((geogrid.length, geogrid.width), invalid, dtype=np.complex64)

    isce3.geocode.geocode_slc(
        geo_data_blocks=[output],
        rdr_data_blocks=[slc.read_complex()],
        dem_raster=dem,
        radargrid=slc.radar_grid,
        geogrid=geogrid,
        orbit=slc.orbit,
        native_doppler=slc.doppler,
        image_grid_doppler=isce3.core.LUT2d(),
        ellipsoid=ellipsoid,
        threshold_geo2rdr=1.0e-8,
        num_iter_geo2rdr=25,
        flatten=flatten,
        invalid_value=invalid,
    )
    _write_complex_geotiff(output, geogrid, output_file, slc, flatten)
    return output_file


def _write_complex_geotiff(
    data: np.ndarray,
    geogrid,
    output_file: Path,
    slc: BiomassSlc,
    flatten: bool,
) -> None:
    transform = Affine(
        geogrid.spacing_x,
        0.0,
        geogrid.start_x,
        0.0,
        geogrid.spacing_y,
        geogrid.start_y,
    )
    with rasterio.open(
        output_file,
        "w",
        driver="GTiff",
        width=geogrid.width,
        height=geogrid.length,
        count=1,
        dtype="complex64",
        crs=f"EPSG:{geogrid.epsg}",
        transform=transform,
        nodata=np.nan,
        tiled=True,
        compress="deflate",
    ) as dataset:
        dataset.write(data, 1)
        dataset.set_band_description(1, slc.polarization)
        dataset.update_tags(
            BIOMASS_PRODUCT_ID=slc.product_id,
            POLARIZATION=slc.polarization,
            WAVELENGTH_METERS=str(slc.wavelength),
            FLATTENED=str(flatten).lower(),
            NATIVE_DOPPLER="annotation geometryDCPolynomial",
        )


def write_stack_provenance(
    path: str | Path,
    slcs: Sequence[BiomassSlc],
    dem_file: str | Path,
    geogrid,
    extent: str,
    flatten: bool,
) -> Path:
    """Write a machine-readable record of the GSLC stack geometry and inputs."""
    path = Path(path)
    record = {
        "products": [slc.product_id for slc in slcs],
        "polarization": slcs[0].polarization,
        "dem": str(Path(dem_file).resolve()),
        "extent_policy": extent,
        "flatten": flatten,
        "wavelength_m": slcs[0].wavelength,
        "geogrid": {
            "epsg": geogrid.epsg,
            "start_x": geogrid.start_x,
            "start_y": geogrid.start_y,
            "spacing_x": geogrid.spacing_x,
            "spacing_y": geogrid.spacing_y,
            "width": geogrid.width,
            "length": geogrid.length,
        },
    }
    path.write_text(json.dumps(record, indent=2) + "\n")
    return path
