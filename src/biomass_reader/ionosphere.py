"""Access the BIOMASS annotation LUT NetCDF, including ionosphere correction.

At P-band (435 MHz) the ionosphere is the dominant error source: dispersive
range delay and Faraday rotation scale as roughly ``(f_C / f_P)**2 ~ 63x`` worse
than Sentinel-1 C-band. The L1a SCS annotation therefore ships an
``ionosphereCorrection`` LUT group (alongside ``geometry``, ``radiometry``,
``denoising``, ``rfiMitigation``). This module exposes those groups as xarray
datasets so a downstream correction step can apply them.

This is intentionally thin: it surfaces the data, it does not (yet) apply a
correction model. The applying step belongs in the geocoding workflow.
"""

from __future__ import annotations

from pathlib import Path

import xarray as xr

LUT_GROUPS = (
    "rfiMitigation",
    "ionosphereCorrection",
    "denoising",
    "geometry",
    "radiometry",
)


def open_lut_group(lut_nc: str | Path, group: str) -> xr.Dataset:
    """Open one group of the annotation LUT NetCDF as an xarray dataset.

    Parameters
    ----------
    lut_nc
        Path to the ``annotation/*.nc`` LUT file.
    group
        Group name, one of :data:`LUT_GROUPS`.

    Returns
    -------
    xarray.Dataset
    """
    assert group in LUT_GROUPS, f"unknown LUT group {group!r}; expected {LUT_GROUPS}"
    return xr.open_dataset(lut_nc, group=group)


def open_ionosphere(lut_nc: str | Path) -> xr.Dataset:
    """Open the ``ionosphereCorrection`` LUT group.

    Parameters
    ----------
    lut_nc
        Path to the ``annotation/*.nc`` LUT file.

    Returns
    -------
    xarray.Dataset
        Ionosphere-correction variables on the annotation LUT grid.
    """
    return open_lut_group(lut_nc, "ionosphereCorrection")
