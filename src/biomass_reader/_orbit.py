"""Build an ``isce3.core.Orbit`` from a BIOMASS navigation (orbit) XML file.

BIOMASS ships CFI-style orbit files (``annotation/navigation/*orb*.xml``) whose
state vectors live under ``Data_Block/List_of_OSVs/OSV`` with ECEF position
``X/Y/Z`` [m] and velocity ``VX/VY/VZ`` [m/s], the same layout as the
Sentinel-1 precise-orbit EOFs.
"""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

from lxml import etree

# Reference epoch offset: place the orbit reference epoch this far before the
# scene start, matching the isce3 / s1-reader convention so radar grid, orbit,
# and Doppler LUT all share one reference epoch.
_REF_EPOCH_OFFSET = timedelta(days=2)


def reference_epoch(sensing_start: datetime) -> datetime:
    """Return the shared orbit/radar-grid reference epoch for a scene."""
    return sensing_start - _REF_EPOCH_OFFSET


def load_orbit(path: str | Path, sensing_start: datetime):
    """Build an ``isce3.core.Orbit`` from a BIOMASS navigation XML file.

    Parameters
    ----------
    path
        Path to the ``*orb*.xml`` orbit file.
    sensing_start
        Scene sensing start, used to set the shared reference epoch.

    Returns
    -------
    isce3.core.Orbit
        Orbit with ECEF state vectors and reference epoch
        ``sensing_start - 2 days``.
    """
    import isce3

    root = etree.parse(str(path)).getroot()
    osvs = root.xpath(".//*[local-name()='OSV']")
    assert osvs, f"no <OSV> state vectors found in {path}"

    def _val(osv: etree._Element, name: str) -> str:
        els = osv.xpath(".//*[local-name()=$n]", n=name)
        assert els, f"OSV missing <{name}>"
        return els[0].text.strip()

    def _utc(osv: etree._Element):
        # Parse straight to isce3 DateTime to keep sub-microsecond precision
        # (python datetime truncates to microseconds).
        s = _val(osv, "UTC")
        for prefix in ("UTC=", "TAI=", "UT1="):
            s = s.removeprefix(prefix)
        return isce3.core.DateTime(s)

    times = [_utc(osv) for osv in osvs]
    positions = [
        [float(_val(o, "X")), float(_val(o, "Y")), float(_val(o, "Z"))] for o in osvs
    ]
    velocities = [
        [float(_val(o, "VX")), float(_val(o, "VY")), float(_val(o, "VZ"))] for o in osvs
    ]

    # BIOMASS orbit timestamps carry ~1 microsecond jitter, which isce3's Orbit
    # rejects (it requires exactly-uniform spacing). Snap to a uniform grid; at
    # ~7.5 km/s a 1 microsecond shift is millimeter-level and negligible.
    t0 = times[0]
    span = times[-1] - t0
    dt = span.total_seconds() / (len(times) - 1)
    deviations = [abs((t - t0).total_seconds() - i * dt) for i, t in enumerate(times)]
    assert max(deviations) < 1e-3, (
        f"orbit state vectors are not uniformly sampled "
        f"(max deviation {max(deviations):.6f} s from {dt:.6f} s spacing)"
    )

    state_vectors = [
        isce3.core.StateVector(t0 + isce3.core.TimeDelta(i * dt), pos, vel)
        for i, (pos, vel) in enumerate(zip(positions, velocities, strict=True))
    ]
    ref_epoch = isce3.core.DateTime(reference_epoch(sensing_start))
    return isce3.core.Orbit(state_vectors, ref_epoch)
