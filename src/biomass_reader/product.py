"""Locate the files inside a BIOMASS L1a SCS product directory.

The SCS product is a Sentinel-1-like folder tree::

    <PRODUCT>.SAFE-style dir/
      <bio...>.xml               # MPH (main product header)
      measurement/
        *abs*.tiff               # 4-band amplitude (HH, HV, VH, VV)
        *phase*.tiff             # 4-band phase
      annotation/
        *.xml                    # main annotation (radar geometry)
        *.nc                     # LUT NetCDF (iono / geometry / radiometry ...)
        navigation/
          *orb*.xml              # orbit state vectors
          *att*.xml              # attitude
      preview/ , schema/
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


def _first(directory: Path, pattern: str) -> Path:
    matches = sorted(directory.glob(pattern))
    assert matches, f"no file matching {pattern!r} in {directory}"
    return matches[0]


@dataclass(frozen=True)
class BiomassProduct:
    """Resolved paths to the components of one BIOMASS L1a SCS product."""

    path: Path
    mph: Path
    annotation_xml: Path
    lut_nc: Path
    orbit_xml: Path
    attitude_xml: Path
    measurement_abs: Path
    measurement_phase: Path

    @classmethod
    def from_dir(cls, path: str | Path) -> BiomassProduct:
        """Resolve product file paths from a product directory.

        Parameters
        ----------
        path
            Path to the extracted SCS product directory.

        Returns
        -------
        BiomassProduct
            Resolved component paths.
        """
        path = Path(path)
        assert path.is_dir(), f"not a directory: {path}"
        measurement = path / "measurement"
        annotation = path / "annotation"
        navigation = annotation / "navigation"
        return cls(
            path=path,
            mph=_first(path, "bio*.xml"),
            annotation_xml=_first(annotation, "*.xml"),
            lut_nc=_first(annotation, "*.nc"),
            orbit_xml=_first(navigation, "*orb*.xml"),
            attitude_xml=_first(navigation, "*att*.xml"),
            measurement_abs=_first(measurement, "*abs*.tiff"),
            measurement_phase=_first(measurement, "*phase*.tiff"),
        )

    @property
    def product_id(self) -> str:
        """The product identifier (directory name)."""
        return self.path.name
