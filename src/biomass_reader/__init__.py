"""isce3-native reader for ESA BIOMASS Level-1a SCS products.

BIOMASS L1a SCS is a Sentinel-1-like SAFE product (slant-range complex data +
annotation XML + orbit file). This package parses it into isce3 objects so the
existing isce3 / dolphin / sweets machinery can geocode and interfere BIOMASS
scenes:

>>> from biomass_reader import BiomassSlc
>>> slc = BiomassSlc.from_dir("BIO_S1_SCS__1S_..._DU1SS4", polarization="HH")
>>> radar_grid = slc.radar_grid   # isce3.product.RadarGridParameters
>>> orbit = slc.orbit             # isce3.core.Orbit
>>> doppler = slc.doppler         # isce3.core.LUT2d (native centroid)
"""

from __future__ import annotations

from ._annotation import AnnotationMetadata, DopplerEstimate, parse_annotation
from .product import BiomassProduct
from .slc import BiomassSlc

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "0.0.0+unknown"

__all__ = [
    "AnnotationMetadata",
    "BiomassProduct",
    "BiomassSlc",
    "DopplerEstimate",
    "__version__",
    "parse_annotation",
]
