"""Shared pytest fixtures for optional BIOMASS real-product tests."""

from __future__ import annotations

import os
from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def product_dirs() -> list[Path]:
    """Return real SCS products selected by ``BIOMASS_TEST_DATA``."""
    value = os.environ.get("BIOMASS_TEST_DATA")
    if value is None:
        pytest.skip("set BIOMASS_TEST_DATA to run real-product tests")
    root = Path(value).expanduser()
    products = sorted(root.glob("BIO_S1_SCS__1S_*"))
    products = [path for path in products if path.is_dir()]
    if not products:
        pytest.fail(f"no extracted SCS products found below {root}")
    return products
