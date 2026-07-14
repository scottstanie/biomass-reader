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
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
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
    try:
        import xarray as xr
    except ImportError as error:
        raise ImportError(
            "opening annotation LUT data requires biomass-reader[ionosphere]"
        ) from error
    assert group in LUT_GROUPS, f"unknown LUT group {group!r}; expected {LUT_GROUPS}"
    dataset = xr.open_dataset(lut_nc, group=group)
    with xr.open_dataset(lut_nc) as root:
        coordinates = {
            dim: root[dim].load()
            for dim in dataset.dims
            if dim in root.variables and root[dim].dims == (dim,)
        }
    return dataset.assign_coords(coordinates)


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


def ionosphere_lut_summary(lut_nc: str | Path) -> dict[str, Any]:
    """Return small, JSON-safe metadata for the ionosphere correction LUT.

    This deliberately reads coordinates and attributes, not the full LUT.  It
    is useful for recording exactly which product-supplied ionosphere layers
    were available when a GSLC was made.
    """
    try:
        from netCDF4 import Dataset
    except ImportError as error:
        raise ImportError(
            "inspecting annotation LUT metadata requires biomass-reader[ionosphere]"
        ) from error
    with Dataset(lut_nc) as root:
        group = root.groups["ionosphereCorrection"]
        variables = {
            name: {
                "dims": list(variable.dimensions),
                "units": getattr(variable, "units", None),
            }
            for name, variable in group.variables.items()
        }
        dimensions = {
            dimension
            for variable in group.variables.values()
            for dimension in variable.dimensions
        }
        coordinates = [
            dimension
            for dimension in dimensions
            if dimension in root.variables
            and root.variables[dimension].dimensions == (dimension,)
        ]
    return {"variables": variables, "coordinates": coordinates}
