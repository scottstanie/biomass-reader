# BIOMASS L1a SCS correction inventory

This matrix distinguishes corrections already performed by the L1 processor
from auxiliary estimates retained in the product. It was checked against the
three T007/F004 products from 2026-04-23 through 2026-04-29, their shipped XSDs,
and the BIOMASS Processing Suite stack processor.

| Effect | Product evidence | Status in test products | Reader action |
|---|---|---:|---|
| Polarimetric calibration | `polarimetricCorrectionFlag` | applied | Do not apply again |
| Faraday rotation | flag plus `faradayRotationCorrectionApplied`; NetCDF angle and uncertainty LUTs | applied | Expose LUT for QA/provenance; do not rotate again |
| Ionospheric phase screen | `ionosphericPhaseScreenCorrectionFlag`; `phaseScreen` in radians | applied | Expose LUT; do not subtract it again |
| Ionospheric group delay | `groupDelayCorrectionFlag`; `rangeShifts` in metres | applied | Expose LUT; do not shift the L1a again |
| Ionospheric azimuth shift | `azimuthShifts` in metres | supplied; no independent applied flag | Treat as provenance pending BPS semantics; do not apply by default |
| Autofocus | `autofocusFlag`, `autofocusShiftsApplied` | not requested/applied | No action |
| Range-spreading loss | `rangeSpreadingLossCompensationFlag` | applied | Do not rescale again |
| Radiometry | `sigmaNought`, `gammaNought` LUTs | normalization LUTs supplied | Not applied to complex GSLCs by this reader |
| Denoising | per-polarization `denoising*` LUTs | noise estimates supplied | Not applied; retain complex measurements |
| RFI | per-polarization frequency masks on RAW axes | diagnostic/intermediate masks supplied | Not directly applicable to focused SLC pixels |
| Geometry | latitude, longitude, height, incidence/elevation and terrain-slope LUTs | supplied | QA/reference; isce3 geocoding uses orbit, radar grid, and user DEM |
| Orbit | navigation XML state vectors | consumed | Build `isce3.core.Orbit`; preserve the product orbit |
| Attitude | navigation XML | supplied | Not currently consumed by `geocode_slc`; document rather than pretend support |
| Topographic flattening | not an L1 correction | performed during GSLC geocoding | `flatten=True` with the user DEM |
| Solid-Earth tide/troposphere | no L1a correction field found | not applied here | Candidate downstream interferometric corrections, not reader functions |

The main practical conclusion is that the ionosphere LUT must not be blindly
applied to these L1a samples. BPS uses the phase-screen and range-shift layers
when assembling and calibrating stacks, but the L1 annotation explicitly says
the corresponding phase and group-delay corrections were performed. Any
residual-ionosphere method should estimate a residual from the stack and label
it separately from the L1 correction.

`biomass_reader.ionosphere.open_lut_group` attaches the root NetCDF coordinate
axes to group variables. `biomass_reader.corrections.parse_correction_status`
reports the applied/requested flags. Both `scripts/biomass_pipeline.py` and
the GSLC-only `scripts/geocode_gslc.py` now check this by default, require the
`phaseScreen` and `rangeShifts` LUT layers, and record the flags/LUT schema in
`stack.json` and GeoTIFF tags. `--ionosphere-policy allow-unapplied` exists for
future residual-correction experiments, but does not apply a correction itself.
