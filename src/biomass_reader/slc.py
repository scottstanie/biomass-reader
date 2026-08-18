"""``BiomassSlc``: a BIOMASS L1a SCS scene as an isce3-ready SLC.

Unlike ``s1reader.Sentinel1BurstSlc`` (a frozen dataclass with ~70 fields),
the isce3 objects are built lazily from a small parsed
:class:`~biomass_reader._annotation.AnnotationMetadata`.
"""

from __future__ import annotations

from datetime import datetime
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from . import _measurement
from ._annotation import AnnotationMetadata, parse_annotation
from ._doppler import build_doppler_lut
from ._orbit import load_orbit, reference_epoch
from .product import BiomassProduct

if TYPE_CHECKING:
    import isce3


class BiomassSlc:
    """One polarization of a BIOMASS L1a SCS scene, as an isce3-ready SLC.

    Parameters
    ----------
    product
        Resolved product file paths.
    metadata
        Parsed annotation metadata.
    polarization
        Selected polarization (``HH``, ``HV``, ``VH``, ``VV``).
    """

    def __init__(
        self,
        product: BiomassProduct,
        metadata: AnnotationMetadata,
        polarization: str,
    ) -> None:
        self._product = product
        self._meta = metadata
        pol = polarization.upper()
        assert pol in metadata.polarizations, (
            f"polarization {pol!r} not in product polarizations "
            f"{metadata.polarizations}"
        )
        self._polarization = pol

    @classmethod
    def from_dir(cls, product_dir: str | Path, polarization: str = "HH") -> BiomassSlc:
        """Load a ``BiomassSlc`` from an extracted SCS product directory.

        Parameters
        ----------
        product_dir
            Path to the SCS product directory.
        polarization
            Polarization to select. Default ``HH``.

        Returns
        -------
        BiomassSlc
        """
        product = BiomassProduct.from_dir(product_dir)
        metadata = parse_annotation(product.annotation_xml, product.mph)
        return cls(product, metadata, polarization)

    # ------------------------------------------------------------------ #
    # SLC protocol
    # ------------------------------------------------------------------ #
    @cached_property
    def radar_grid(self) -> isce3.product.RadarGridParameters:
        """isce3 radar grid for the full-resolution SCS scene."""
        import isce3

        m = self._meta
        ref_dt = reference_epoch(m.sensing_start)
        ref = isce3.core.DateTime(ref_dt)
        sensing_start = (m.sensing_start - ref_dt).total_seconds()
        return isce3.product.RadarGridParameters(
            sensing_start,
            m.wavelength,
            m.prf,
            m.starting_range,
            m.range_pixel_spacing,
            self.look_side,
            m.number_of_lines,
            m.number_of_samples,
            ref,
        )

    @cached_property
    def orbit(self) -> isce3.core.Orbit:
        """isce3 orbit from the navigation file, sharing the scene ref epoch."""
        return load_orbit(self._product.orbit_xml, self._meta.sensing_start)

    @cached_property
    def doppler(self) -> isce3.core.LUT2d:
        """Native Doppler-centroid LUT (non-zero, like Sentinel-1)."""
        return build_doppler_lut(self._meta)

    @property
    def azimuth_carrier(self) -> isce3.core.Poly2d:
        """Azimuth carrier polynomial (zero: BIOMASS is stripmap, not TOPS)."""
        import isce3

        return isce3.core.Poly2d(np.zeros((1, 1)))

    @property
    def wavelength(self) -> float:
        """Radar wavelength [m]."""
        return self._meta.wavelength

    @property
    def center_frequency(self) -> float:
        """Radar carrier frequency [Hz]."""
        return self._meta.carrier_frequency

    @property
    def sensing_start(self) -> datetime:
        """Azimuth time of the first line (UTC)."""
        return self._meta.sensing_start

    @property
    def shape(self) -> tuple[int, int]:
        """Scene shape (azimuth lines, range samples)."""
        return (self._meta.number_of_lines, self._meta.number_of_samples)

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        """Geographic bounds (west, south, east, north) from the footprint."""
        return self._meta.footprint.bounds

    @cached_property
    def center_target(self) -> np.ndarray:
        """Approximate the scene-center ground target in ECEF coordinates.

        BIOMASS does not store a distinguished center-pixel ECEF coordinate in
        the main annotation.
        The footprint centroid at zero ellipsoidal height.
        """
        import isce3

        center = self._meta.footprint.centroid
        xyz = isce3.core.Ellipsoid().lon_lat_to_xyz(
            [np.deg2rad(center.x), np.deg2rad(center.y), 0.0]
        )
        return np.asarray(xyz, dtype=np.float64)

    @property
    def identifier(self) -> str:
        """Unique identifier: ``<product_id>_<pol>``."""
        return f"{self._product.product_id}_{self._polarization}"

    @property
    def product(self) -> BiomassProduct:
        """Resolved files belonging to the source L1a product."""
        return self._product

    @property
    def product_id(self) -> str:
        """Source L1a product identifier."""
        return self._product.product_id

    @property
    def swath(self) -> str:
        """BIOMASS acquisition swath identifier."""
        return self._meta.swath

    @property
    def look_side(self) -> str:
        """Look direction (``left`` for BIOMASS)."""
        return self._meta.look_side

    @property
    def polarization(self) -> str:
        """Selected polarization."""
        return self._polarization

    @property
    def pixel_spacing(self) -> tuple[float, float]:
        """Pixel spacing (azimuth_m, range_m)."""
        return (self._meta.azimuth_pixel_spacing, self._meta.range_pixel_spacing)

    @property
    def slc_path(self) -> Path:
        """Path to the product-supplied four-band complex VRT."""
        return self._product.measurement_vrt

    # ------------------------------------------------------------------ #
    # Data access
    # ------------------------------------------------------------------ #
    def read_complex(self, window=None) -> np.ndarray:
        """Read the complex SLC (``amplitude * exp(1j*phase)``) for this pol."""
        return _measurement.read_complex(
            self._product.measurement_abs,
            self._product.measurement_phase,
            self._polarization,
            window=window,
        )

    def to_geotiff(self, out_path: str | Path) -> Path:
        """Write the complex SLC to a CFloat32 GeoTIFF (radar coordinates)."""
        return _measurement.write_geotiff(self.read_complex(), out_path)

    def __repr__(self) -> str:
        return (
            f"BiomassSlc({self._product.product_id!r}, "
            f"pol={self._polarization!r}, shape={self.shape})"
        )
