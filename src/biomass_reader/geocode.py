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
from rasterio.crs import CRS
from rasterio.warp import transform

from .slc import BiomassSlc


def _posting(posting: float | tuple[float, float]) -> tuple[float, float]:
    """Normalize a scalar legacy spacing or an ``(x, y)`` posting."""
    if isinstance(posting, float | int):
        posting_x = posting_y = float(posting)
    else:
        posting_x, posting_y = (float(value) for value in posting)
    if posting_x <= 0 or posting_y <= 0:
        raise ValueError("posting values must be positive")
    return posting_x, posting_y


def projected_sampling(slc: BiomassSlc, epsg: int) -> tuple[float, float]:
    """Estimate native GCP sampling along projected easting and northing.

    The producer VRT's ground-control points form a radar-to-map affine fit.
    For each map coordinate, the larger contribution from the azimuth/range
    lattice is the limiting native sampling scale. This handles BIOMASS's
    rotated swath rather than assuming range is easting and azimuth northing.
    """
    target_crs = CRS.from_epsg(epsg)
    if not target_crs.is_projected:
        raise ValueError(f"native posting requires a projected CRS, not EPSG:{epsg}")
    with rasterio.open(slc.slc_path) as dataset:
        gcps, source_crs = dataset.gcps
    if source_crs is None or len(gcps) < 3:
        raise ValueError(f"{slc.product_id} has insufficient VRT ground-control points")

    x, y = transform(
        source_crs,
        target_crs,
        [gcp.x for gcp in gcps],
        [gcp.y for gcp in gcps],
    )
    design = np.column_stack(
        (
            np.ones(len(gcps)),
            [gcp.row for gcp in gcps],
            [gcp.col for gcp in gcps],
        )
    )
    x_coefficients = np.linalg.lstsq(design, x, rcond=None)[0]
    y_coefficients = np.linalg.lstsq(design, y, rcond=None)[0]
    jacobian = np.array([x_coefficients[1:], y_coefficients[1:]])
    sampling = np.max(np.abs(jacobian), axis=1)
    return float(sampling[0]), float(sampling[1])


def native_posting(slcs: Sequence[BiomassSlc], epsg: int) -> tuple[float, float]:
    """Choose conservative posting from a standard grid series.

    For native sampling < 10m: uses the ``10 / 2**n`` series (10, 5, 2.5, 1.25, ...).
    For native sampling >= 10m: uses 10m increments (10, 20, 30, 40, ..., 100).

    Each output-axis posting is the largest value from the appropriate series that
    is not coarser than the least-well-sampled source scene in that direction.

    Examples:
        - 3m native (Sentinel-1 range) → 2.5m posting (10/2^n series)
        - 14m native (Sentinel-1 azimuth) → 10m posting (10/2^n series)
        - 46m native (BIOMASS range) → 40m posting (10m increment series, 1.15x upsampling)
        - 7m native (BIOMASS azimuth) → 5m posting (10/2^n series, 1.4x upsampling)

    This avoids excessive upsampling (>5x) that would occur if the 10/2^n series
    were used for coarse-resolution sensors. BIOMASS L1A native ground resolution
    is typically 6-10m (azimuth) and 25-50m (range, geometry-dependent), with
    range resolution often exceeding 10m at high latitudes or far range.
    """
    if not slcs:
        raise ValueError("at least one SLC is required")
    sampling = np.min([projected_sampling(slc, epsg) for slc in slcs], axis=0)

    def standard_posting(value: float) -> float:
        """Select posting from standard series based on native sampling."""
        if value < 10.0:
            # Fine resolution: use 10/2^n series (10, 5, 2.5, 1.25, ...)
            candidate = 10.0
            while candidate > value:
                candidate /= 2.0
            return candidate
        else:
            # Coarse resolution: use 10m increments (10, 20, 30, ..., 100)
            # Round DOWN to nearest 10m to ensure posting ≤ native (avoid coarser than native)
            # Capped at 100m for very coarse sensors
            return min(math.floor(value / 10.0) * 10.0, 100.0)

    return float(standard_posting(sampling[0])), float(standard_posting(sampling[1]))


def make_shared_geogrid(
    slcs: Sequence[BiomassSlc],
    epsg: int,
    posting: float | tuple[float, float],
    extent: Literal["union", "intersection"] = "union",
):
    """Create one snapped geogrid covering a stack of BIOMASS acquisitions."""
    import isce3

    if not slcs:
        raise ValueError("at least one SLC is required")
    spacing_x, spacing_y = _posting(posting)

    grids = [
        isce3.product.bbox_to_geogrid(
            slc.radar_grid,
            slc.orbit,
            slc.doppler,
            spacing_x,
            -spacing_y,
            epsg,
        )
        for slc in slcs
    ]
    xmins = [grid.start_x for grid in grids]
    xmaxs = [grid.start_x + grid.width * spacing_x for grid in grids]
    ymaxs = [grid.start_y for grid in grids]
    ymins = [grid.start_y - grid.length * spacing_y for grid in grids]

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
    xmin = math.floor(xmin / spacing_x) * spacing_x
    xmax = math.ceil(xmax / spacing_x) * spacing_x
    ymin = math.floor(ymin / spacing_y) * spacing_y
    ymax = math.ceil(ymax / spacing_y) * spacing_y
    width = int(round((xmax - xmin) / spacing_x))
    length = int(round((ymax - ymin) / spacing_y))
    return isce3.product.GeoGridParameters(
        xmin, ymax, spacing_x, -spacing_y, width, length, epsg
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
    ionosphere: list[dict] | None = None,
    posting_mode: str = "explicit",
) -> Path:
    """Write a machine-readable record of the GSLC stack geometry and inputs."""
    path = Path(path)
    record = {
        "products": [slc.product_id for slc in slcs],
        "polarization": slcs[0].polarization,
        "dem": str(Path(dem_file).resolve()),
        "extent_policy": extent,
        "posting_mode": posting_mode,
        "flatten": flatten,
        "wavelength_m": slcs[0].wavelength,
        "ionosphere": ionosphere,
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
