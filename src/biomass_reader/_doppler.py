"""Build the native Doppler-centroid ``LUT2d`` for a BIOMASS SCS scene.

Unlike NISAR (zero-Doppler focused), BIOMASS SCS carries a non-zero Doppler
centroid described by per-azimuth polynomials in slant-range time. This mirrors
the Sentinel-1 situation, so ``geocode_slc`` must be given a real
``native_doppler`` LUT rather than an empty one.
"""

from __future__ import annotations

import numpy as np

from ._annotation import AnnotationMetadata
from ._constants import SPEED_OF_LIGHT
from ._orbit import reference_epoch


def build_doppler_lut(meta: AnnotationMetadata, n_range: int = 40):
    """Build a native Doppler-centroid ``LUT2d`` on a uniform (az, range) grid.

    The BIOMASS Doppler estimates are polynomials in two-way slant-range time,
    sampled at a handful of azimuth times. This evaluates them onto a uniform
    grid (slant range in meters vs. azimuth time in seconds-since-ref-epoch) so
    the result can be handed directly to ``isce3.geocode.geocode_slc``.

    Parameters
    ----------
    meta
        Parsed annotation metadata.
    n_range
        Number of slant-range samples in the LUT grid.

    Returns
    -------
    isce3.core.LUT2d
        Doppler centroid [Hz] as a function of (slant range, azimuth time),
        sharing the scene reference epoch. ``bounds_error`` is disabled so
        geocoding can query slightly outside the fitted extent.
    """
    import isce3

    ref = reference_epoch(meta.sensing_start)

    # Slant-range grid [m] -> two-way slant-range time [s] for polynomial eval.
    r0 = meta.starting_range
    r1 = r0 + meta.range_pixel_spacing * (meta.number_of_samples - 1)
    ranges = np.linspace(r0, r1, n_range)
    slant_times = 2.0 * ranges / SPEED_OF_LIGHT

    # Azimuth grid [s since ref epoch] at the estimate times.
    est = meta.doppler_estimates
    az_times = np.array([(e.azimuth_time - ref).total_seconds() for e in est])

    rows = []
    for e in est:
        # np.polyval wants highest-power-first; coefficients are ascending.
        rows.append(np.polyval(e.coefficients[::-1], slant_times - e.t0))
    data = np.asarray(rows, dtype=np.float64)  # shape (n_az, n_range)

    # A single azimuth estimate cannot define a grid; duplicate it across the
    # scene so the LUT spans the full azimuth extent.
    if data.shape[0] == 1:
        stop = (meta.sensing_stop - ref).total_seconds()
        az_times = np.array([az_times[0], stop])
        data = np.vstack([data, data])

    dy = (az_times[-1] - az_times[0]) / (len(az_times) - 1)
    dx = (ranges[-1] - ranges[0]) / (len(ranges) - 1)
    return isce3.core.LUT2d(
        float(ranges[0]),
        float(az_times[0]),
        float(dx),
        float(dy),
        data,
        "bilinear",
        False,
    )
