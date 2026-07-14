#!/usr/bin/env python
"""Geocode one or more BIOMASS L1a SCS products to complex GSLC GeoTIFFs.

This is the GSLC-only counterpart to ``biomass_pipeline.py``: no Dolphin,
interferograms, unwrapping, or time-series inversion are run.  By default it
verifies that the source product's ionospheric phase-screen and group-delay
corrections were already applied and records the correction/LUT provenance in
both each GSLC and ``stack.json``.

Example
-------
python scripts/geocode_gslc.py \
    --products /data/BIO_S1_SCS__1S_* --dem /data/dem_utm.tif \
    --work /data/biomass/gslc --polarization HH --posting 30 30
"""

from __future__ import annotations

from biomass_pipeline import make_parser, run


def main() -> None:
    """Parse GSLC-only arguments and permit a single source product."""
    parser = make_parser(include_dolphin=False, description=__doc__)
    run(parser.parse_args(), min_products=1)


if __name__ == "__main__":
    main()
