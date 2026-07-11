"""Read BIOMASS L1a SCS complex measurement data.

The SCS measurement is stored as two 4-band GeoTIFFs -- amplitude (``*abs*``)
and phase (``*phase*``) -- with bands ordered ``HH, HV, VH, VV``. The complex
SLC for a polarization is reconstructed as ``amplitude * exp(1j * phase)``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window

from ._constants import POLARIZATION_ORDER


def _band_index(polarization: str) -> int:
    pol = polarization.upper()
    assert pol in POLARIZATION_ORDER, (
        f"unknown polarization {polarization!r}; expected one of {POLARIZATION_ORDER}"
    )
    return POLARIZATION_ORDER.index(pol) + 1  # rasterio bands are 1-indexed


def read_complex(
    abs_path: str | Path,
    phase_path: str | Path,
    polarization: str,
    window: Window | None = None,
) -> np.ndarray:
    """Read one polarization as a complex64 array.

    Parameters
    ----------
    abs_path, phase_path
        Paths to the amplitude and phase measurement GeoTIFFs.
    polarization
        One of ``HH``, ``HV``, ``VH``, ``VV``.
    window
        Optional rasterio window to read a spatial subset.

    Returns
    -------
    numpy.ndarray
        Complex64 SLC array (azimuth, range).
    """
    band = _band_index(polarization)
    with rasterio.open(abs_path) as src:
        amp = src.read(band, window=window)
    with rasterio.open(phase_path) as src:
        phase = src.read(band, window=window)
    return (amp * np.exp(1j * phase)).astype(np.complex64)


def write_geotiff(data: np.ndarray, out_path: str | Path) -> Path:
    """Write a complex SLC array to a CFloat32 GeoTIFF (radar coordinates).

    No geotransform is attached: the array is in slant-range/azimuth radar
    geometry. Geocoding is a downstream step (``isce3.geocode.geocode_slc``).
    """
    out_path = Path(out_path)
    profile = {
        "driver": "GTiff",
        "height": data.shape[0],
        "width": data.shape[1],
        "count": 1,
        "dtype": "complex64",
        "compress": "deflate",
        "tiled": True,
    }
    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(data, 1)
    return out_path
