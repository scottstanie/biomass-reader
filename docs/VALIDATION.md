# T007/F004 validation report

## Inputs and environment

Three BPS 4.4.2 `S1_SCS__1S` products were processed:

- 2026-04-23, repeat cycle C01;
- 2026-04-26, repeat cycle C02;
- 2026-04-29, repeat cycle C03.

The run used Python 3.12, isce3 0.25.8, Dolphin
0.42.3.post1.dev77, Whirlwind 0.5.0, EPSG:32620, a 30 m output
grid, topographic flattening, and 6 by 6 looks for the independent validation
statistics.

```bash
export BIOMASS_DATA=/Volumes/WD_BLACK_SN7100_4TB/Documents/Learning/sweets-testing/biomass
export PYTHONPATH=src

python scripts/biomass_pipeline.py \
  --products "$BIOMASS_DATA"/BIO_S1_SCS__1S_*/ \
  --dem "$BIOMASS_DATA"/dem_biomass_utm.tif \
  --work "$BIOMASS_DATA"/validation-v2 \
  --polarization HH --posting 30 30 --run-dolphin

python scripts/validate_stack.py \
  "$BIOMASS_DATA"/validation-v2/gslc/*.tif \
  --output-dir "$BIOMASS_DATA"/validation-v2/diagnostics \
  --dolphin-dir "$BIOMASS_DATA"/validation-v2/dolphin --looks 6

python scripts/polarimetric_phase.py \
  "$BIOMASS_DATA"/validation-v2/gslc \
  "$BIOMASS_DATA"/validation-v2-HV/gslc \
  "$BIOMASS_DATA"/validation-v2-VH/gslc \
  "$BIOMASS_DATA"/validation-v2-VV/gslc \
  --output-dir "$BIOMASS_DATA"/validation-v2/polarimetric --looks 6
```

## Geometry and samples

- All four reconstructed complex channels agree with the BPS VRT within
  float32 numerical tolerance.
- The 48 VRT GCPs per product agree with zero-Doppler isce3 geo2rdr at the 95th
  percentile to 0.38 azimuth pixel or better and 0.00034 range pixel or better.
- All GSLCs are 5,123 by 3,638 pixels with the same transform, CRS, and NaN
  nodata declaration.
- The geogrid is a projection-snapped union: origin (285210 m, 468750 m),
  spacing (30 m, -30 m), EPSG:32620.

The available DEM does not cover the entire conservative western geogrid bound;
isce3 reports this explicitly. Those pixels remain NaN. A production example
should fetch a slightly larger DEM, but this does not affect the common valid
footprint used for the statistics below.

## Interferometric results

| Pair | Valid fraction | Median coherence | Fraction of valid pixels ≥ 0.5 |
|---|---:|---:|---:|
| 2026-04-23–26 | 0.5108 | 0.4499 | 0.3890 |
| 2026-04-23–29 | 0.5036 | 0.2392 | 0.0717 |
| 2026-04-26–29 | 0.5121 | 0.4506 | 0.3888 |

Among 15,576 looked pixels where all three pair coherences are at least 0.5,
the median absolute closure phase is 0.1017 rad and the 90th percentile is
0.2687 rad.

Dolphin plus Whirlwind completed in 28 seconds with 1.64 GB reported peak
memory. The fraction of valid pixels assigned to a connected component was
0.9997, 0.9813, and 0.9996 for the three pairs, respectively. Each retained one
connected component. The metric excludes the declared uint16 nodata value
65535.

## Full-polarimetric diagnostic

HH, HV, VH, and VV were independently geocoded to the same grid. HV/VH
reciprocity is strong: median coherence is approximately 0.991 and median
absolute phase difference is 0.021–0.022 rad. HH/VV median coherence is about
0.415 and its spatial circular mean phase remains near -0.06 rad across the
three dates. Interpretation and proposed external validation are in
[`POLARIMETRIC_RESEARCH.md`](POLARIMETRIC_RESEARCH.md).

Generated results are retained under `validation-v2` beside the source data,
including `stack.json`, `dolphin.yaml`, GSLCs, interferograms, unwrapped rasters,
`diagnostics/validation.json`, `diagnostics/validation.png`, and the
polarimetric JSON/PNG.
