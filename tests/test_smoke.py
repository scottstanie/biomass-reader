"""Smoke tests that do not require a downloaded product.

Product-level tests (parsing a real annotation, geocoding a scene) are added
once a sample SCS product is available in the test-data area.
"""

from __future__ import annotations

import numpy as np

import biomass_reader
from biomass_reader._annotation import _parse_time


def test_imports():
    """Top-level API is importable."""
    assert hasattr(biomass_reader, "BiomassSlc")
    assert hasattr(biomass_reader, "parse_annotation")


def test_parse_time_strips_prefix():
    """``timeType`` parsing tolerates the CFI ``UTC=`` prefix."""
    a = _parse_time("UTC=2026-04-23T10:04:58.533000")
    b = _parse_time("2026-04-23T10:04:58.533000")
    assert a == b
    assert a.year == 2026 and a.microsecond == 533000


def test_polarization_band_index():
    """Measurement band order maps HH,HV,VH,VV -> 1,2,3,4."""
    from biomass_reader._measurement import _band_index

    assert (_band_index("HH"), _band_index("VV")) == (1, 4)


def test_wavelength_from_carrier():
    """Wavelength derives from carrier frequency (P-band ~ 0.69 m)."""
    from biomass_reader._constants import NOMINAL_CARRIER_FREQUENCY_HZ, SPEED_OF_LIGHT

    lam = SPEED_OF_LIGHT / NOMINAL_CARRIER_FREQUENCY_HZ
    assert np.isclose(lam, 0.689, atol=0.01)
