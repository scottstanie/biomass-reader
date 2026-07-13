"""Regression tests against locally staged ESA BIOMASS products."""

from __future__ import annotations

import importlib.util

import numpy as np
import pytest
import rasterio
from rasterio.windows import Window

from biomass_reader import (
    BiomassSlc,
    gcp_radar_coordinate_errors,
    parse_correction_status,
)
from biomass_reader.ionosphere import open_ionosphere
from biomass_reader.product import BiomassProduct


@pytest.mark.parametrize(
    "polarization,band", [("HH", 1), ("HV", 2), ("VH", 3), ("VV", 4)]
)
def test_complex_samples_match_product_vrt(product_dirs, polarization, band):
    """Amplitude/phase reconstruction agrees with the producer's VRT."""
    product_dir = product_dirs[0]
    slc = BiomassSlc.from_dir(product_dir, polarization)
    window = Window(100, 100, 16, 16)
    actual = slc.read_complex(window=window)
    with rasterio.open(slc.slc_path) as dataset:
        expected = dataset.read(band, window=window)
        assert dataset.descriptions[band - 1] == polarization
    np.testing.assert_allclose(actual, expected, rtol=1e-6, atol=1e-7)


def test_real_product_geometry_and_corrections(product_dirs):
    """All staged repeat acquisitions expose consistent radar metadata."""
    slcs = [BiomassSlc.from_dir(path) for path in product_dirs]
    assert all(slc.shape == (20_767, 1_367) for slc in slcs)
    assert all(slc.look_side == "left" for slc in slcs)
    assert all(np.isclose(slc.center_frequency, 435e6) for slc in slcs)
    assert all(slc.center_target.shape == (3,) for slc in slcs)
    assert all(np.isfinite(slc.center_target).all() for slc in slcs)
    assert all(slc.radar_grid.length == slc.shape[0] for slc in slcs)
    assert all(slc.radar_grid.width == slc.shape[1] for slc in slcs)
    assert all(slc.orbit.size > 2 for slc in slcs)

    for slc in slcs:
        status = parse_correction_status(slc.product.annotation_xml)
        assert status.polarimetric
        assert status.faraday_rotation
        assert status.ionospheric_phase_screen
        assert status.ionospheric_group_delay
        assert status.range_spreading_loss
        assert status.faraday_rotation_applied
        assert not status.autofocus_requested


def test_lut_coordinates_are_attached(product_dirs):
    """Group variables retain the coordinate axes stored at NetCDF root."""
    product = BiomassProduct.from_dir(product_dirs[0])
    dataset = open_ionosphere(product.lut_nc)
    assert dataset.phaseScreen.attrs["units"] == "rad"
    assert "relativeAzimuthTimeSLC" in dataset.coords
    assert "slantRangeTimeSLC" in dataset.coords
    assert dataset.phaseScreen.dims == (
        "relativeAzimuthTimeSLC",
        "slantRangeTimeSLC",
    )


def test_zero_doppler_geometry_matches_product_gcps(product_dirs):
    """Orbit/radar-grid geometry reproduces the producer's VRT GCPs."""
    for product_dir in product_dirs:
        errors = np.abs(gcp_radar_coordinate_errors(BiomassSlc.from_dir(product_dir)))
        assert np.percentile(errors[:, 0], 95) < 0.5
        assert np.percentile(errors[:, 1], 95) < 0.001


@pytest.mark.skipif(
    importlib.util.find_spec("sarlet") is None,
    reason="sarlet is not installed",
)
def test_satisfies_sarlet_slc_protocol(product_dirs):
    """The reader conforms to sarlet's runtime-checkable SLC protocol."""
    from sarlet._types import SLC

    assert isinstance(BiomassSlc.from_dir(product_dirs[0]), SLC)
