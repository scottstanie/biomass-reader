"""Parse the BIOMASS L1a SCS main annotation XML into radar-geometry metadata.

The element paths below are derived from the official BIOMASS Processing Suite
schemas (``bio-l1ab-main-annotation.xsd`` / ``bio-l1-annotations.xsd``). Parsing
is namespace-agnostic (matched on ``local-name()``) so it is robust to the
concrete namespace prefixes a given processor baseline emits.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
from dateutil.parser import isoparse
from lxml import etree
from shapely.geometry import Polygon

from ._constants import SPEED_OF_LIGHT


@dataclass(frozen=True)
class DopplerEstimate:
    """A single Doppler-centroid estimate at one azimuth time.

    The centroid is a polynomial in ``(slant_range_time - t0)``.

    Attributes
    ----------
    azimuth_time
        Azimuth time of the estimate (UTC).
    t0
        Reference two-way slant-range time of the polynomial [s].
    coefficients
        Polynomial coefficients ``[c0, c1, ...]`` (ascending power), giving
        Doppler centroid [Hz] as ``sum_k c_k * (tau - t0)**k``.
    """

    azimuth_time: datetime
    t0: float
    coefficients: np.ndarray


@dataclass(frozen=True)
class AnnotationMetadata:
    """Radar-geometry metadata parsed from a BIOMASS L1a SCS annotation.

    This is the minimal set needed to build an ``isce3`` radar grid, Doppler
    LUT, and to interpret the measurement raster. Orbit state vectors live in a
    separate navigation file (see :mod:`biomass_reader._orbit`).
    """

    sensing_start: datetime
    sensing_stop: datetime
    azimuth_time_interval: float
    range_pixel_spacing: float
    azimuth_pixel_spacing: float
    first_sample_slant_range_time: float
    number_of_samples: int
    number_of_lines: int
    carrier_frequency: float
    look_side: str
    polarizations: tuple[str, ...]
    mission_phase: str
    swath: str
    footprint: Polygon
    doppler_estimates: tuple[DopplerEstimate, ...]

    @property
    def starting_range(self) -> float:
        """One-way starting slant range [m] of the first sample."""
        return 0.5 * SPEED_OF_LIGHT * self.first_sample_slant_range_time

    @property
    def wavelength(self) -> float:
        """Radar wavelength [m] from the carrier frequency."""
        return SPEED_OF_LIGHT / self.carrier_frequency

    @property
    def prf(self) -> float:
        """Effective azimuth PRF [Hz] of the SCS grid."""
        return 1.0 / self.azimuth_time_interval


def _findall(root: etree._Element, name: str) -> list[etree._Element]:
    """Return all descendants whose local (namespace-stripped) tag is ``name``."""
    return root.xpath(".//*[local-name()=$n]", n=name)


def _text(root: etree._Element, name: str) -> str:
    els = _findall(root, name)
    assert els, f"annotation missing required element <{name}>"
    txt = els[0].text
    assert txt is not None, f"annotation element <{name}> is empty"
    return txt.strip()


def _float(root: etree._Element, name: str) -> float:
    return float(_text(root, name))


def _int(root: etree._Element, name: str) -> int:
    return int(_text(root, name))


def _time(root: etree._Element, name: str) -> datetime:
    return _parse_time(_text(root, name))


def _parse_time(value: str) -> datetime:
    """Parse a BIOMASS ``timeType`` value, tolerating a ``UTC=`` prefix."""
    value = value.strip()
    for prefix in ("UTC=", "TAI=", "UT1="):
        if value.startswith(prefix):
            value = value[len(prefix) :]
    return isoparse(value)


def _parse_footprint(root: etree._Element) -> Polygon:
    """Parse ``sarImage/footprint`` (flat ``lat lon lat lon ...``) to a polygon."""
    coords = [float(v) for v in _text(root, "footprint").split()]
    assert len(coords) % 2 == 0, "footprint must have an even number of values"
    points = [(coords[i + 1], coords[i]) for i in range(0, len(coords), 2)]
    return Polygon(points)


def _doppler_coeffs(est: etree._Element) -> np.ndarray:
    """Return the usable Doppler polynomial coefficients for one estimate.

    Prefers ``combinedDCPolynomial`` (data-derived centroid); real products
    often leave it empty and populate only ``geometryDCPolynomial``
    (geometry-only), so fall back to that.
    """
    for name in ("combinedDCPolynomial", "geometryDCPolynomial"):
        els = est.xpath(".//*[local-name()=$n]", n=name)
        if els and els[0].text and els[0].text.strip():
            return np.array([float(v) for v in els[0].text.split()])
    raise AssertionError("dcEstimate has no usable Doppler polynomial")


def _parse_doppler(root: etree._Element) -> tuple[DopplerEstimate, ...]:
    """Parse the ``dopplerParameters/dcEstimateList`` into estimates."""
    estimates = [
        DopplerEstimate(
            azimuth_time=_time(est, "azimuthTime"),
            t0=_float(est, "t0"),
            coefficients=_doppler_coeffs(est),
        )
        for est in _findall(root, "dcEstimate")
    ]
    assert estimates, "annotation contains no Doppler estimates"
    return tuple(estimates)


def parse_annotation(path: str | Path, mph_path: str | Path) -> AnnotationMetadata:
    """Parse a BIOMASS L1a SCS main annotation XML file.

    Parameters
    ----------
    path
        Path to the ``annotation/*.xml`` main annotation file.  mph_path
        Path to the product MPH (``bio*.xml`` at the product root). The look
        direction is carried in the MPH (``sar:antennaLookDirection``), not the
        annotation.

    Returns
    -------
    AnnotationMetadata
        Parsed radar-geometry metadata.
    """
    root = etree.parse(str(path)).getroot()
    mph = etree.parse(str(mph_path)).getroot()

    look_els = mph.xpath(".//*[local-name()='antennaLookDirection']")
    assert look_els and look_els[0].text, "MPH missing antennaLookDirection"
    look_side = look_els[0].text.strip().lower()

    polarizations = tuple(
        el.text.strip() for el in _findall(root, "polarisation") if el.text
    )

    return AnnotationMetadata(
        sensing_start=_time(root, "firstLineAzimuthTime"),
        sensing_stop=_time(root, "lastLineAzimuthTime"),
        azimuth_time_interval=_float(root, "azimuthTimeInterval"),
        range_pixel_spacing=_float(root, "rangePixelSpacing"),
        azimuth_pixel_spacing=_float(root, "azimuthPixelSpacing"),
        first_sample_slant_range_time=_float(root, "firstSampleSlantRangeTime"),
        number_of_samples=_int(root, "numberOfSamples"),
        number_of_lines=_int(root, "numberOfLines"),
        carrier_frequency=_float(root, "radarCarrierFrequency"),
        look_side=look_side,
        polarizations=polarizations,
        mission_phase=_text(root, "missionPhaseID"),
        swath=_text(root, "swath"),
        footprint=_parse_footprint(root),
        doppler_estimates=_parse_doppler(root),
    )
