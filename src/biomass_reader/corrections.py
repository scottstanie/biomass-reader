"""Report corrections recorded in a BIOMASS L1a SCS annotation.

The L1 annotation distinguishes requested/performed corrections from auxiliary
LUTs retained for provenance and downstream stack processing.  In particular,
the presence of ``phaseScreen`` or ``faradayRotation`` in the NetCDF file does
not mean that an L1a consumer should apply it again.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from lxml import etree


@dataclass(frozen=True)
class CorrectionStatus:
    """Correction flags reported by one BIOMASS L1a SCS product."""

    polarimetric: bool
    faraday_rotation: bool
    ionospheric_phase_screen: bool
    ionospheric_group_delay: bool
    range_spreading_loss: bool
    autofocus_requested: bool
    faraday_rotation_applied: bool
    autofocus_shifts_applied: bool


def _bool(root: etree._Element, name: str) -> bool:
    elements = root.xpath(".//*[local-name()=$name]", name=name)
    if not elements or elements[0].text is None:
        raise ValueError(f"annotation missing required correction flag <{name}>")
    value = elements[0].text.strip().lower()
    if value not in {"true", "false"}:
        raise ValueError(f"invalid boolean value for <{name}>: {value!r}")
    return value == "true"


def parse_correction_status(annotation_xml: str | Path) -> CorrectionStatus:
    """Parse correction status from a BIOMASS main annotation XML file."""
    root = etree.parse(str(annotation_xml)).getroot()
    return CorrectionStatus(
        polarimetric=_bool(root, "polarimetricCorrectionFlag"),
        faraday_rotation=_bool(root, "faradayRotationCorrectionFlag"),
        ionospheric_phase_screen=_bool(
            root, "ionosphericPhaseScreenCorrectionFlag"
        ),
        ionospheric_group_delay=_bool(root, "groupDelayCorrectionFlag"),
        range_spreading_loss=_bool(root, "rangeSpreadingLossCompensationFlag"),
        autofocus_requested=_bool(root, "autofocusFlag"),
        faraday_rotation_applied=_bool(root, "faradayRotationCorrectionApplied"),
        autofocus_shifts_applied=_bool(root, "autofocusShiftsApplied"),
    )
