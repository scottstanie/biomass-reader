"""Focused quality checks against metadata shipped with BIOMASS products."""

from __future__ import annotations

import numpy as np
import rasterio

from .slc import BiomassSlc


def gcp_radar_coordinate_errors(slc: BiomassSlc) -> np.ndarray:
    """Compare product VRT GCPs with isce3 zero-Doppler geo2rdr geometry.

    Returns
    -------
    numpy.ndarray
        Array shaped ``(number_of_gcps, 2)`` containing errors in
        ``(azimuth_line, range_sample)`` pixels.
    """
    import isce3

    with rasterio.open(slc.slc_path) as dataset:
        gcps, _ = dataset.gcps
    if not gcps:
        raise ValueError(f"product VRT contains no GCPs: {slc.slc_path}")

    radar_grid = slc.radar_grid
    errors = []
    for gcp in gcps:
        azimuth_time, slant_range = isce3.geometry.geo2rdr(
            [np.deg2rad(gcp.x), np.deg2rad(gcp.y), gcp.z],
            orbit=slc.orbit,
            doppler=isce3.core.LUT2d(),
            wavelength=slc.wavelength,
            side=slc.look_side,
        )
        line = (
            azimuth_time - radar_grid.sensing_start
        ) / radar_grid.az_time_interval
        sample = (
            slant_range - radar_grid.starting_range
        ) / radar_grid.range_pixel_spacing
        errors.append((line - gcp.row, sample - gcp.col))
    return np.asarray(errors)
