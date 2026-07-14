# biomass-reader

An isce3-native reader for **ESA BIOMASS Level-1a SCS** (single-look complex,
slant-range) products.

BIOMASS L1a SCS is structurally a *Sentinel-1-like SAFE product*: slant-range
complex data + an annotation XML + a CFI orbit file, in
`measurement/ annotation/ preview/ schema/` folders. This package plays the role
`s1-reader` plays for Sentinel-1 — it parses the product annotation into isce3
objects so the existing isce3 / dolphin / sweets stack can geocode and interfere
BIOMASS scenes with no bespoke SAR code.

```python
from biomass_reader import BiomassSlc

slc = BiomassSlc.from_dir("BIO_S1_SCS__1S_..._DU1SS4", polarization="HH")

radar_grid = slc.radar_grid    # isce3.product.RadarGridParameters
orbit      = slc.orbit         # isce3.core.Orbit
doppler    = slc.doppler       # isce3.core.LUT2d  (native centroid, non-zero)
data       = slc.read_complex()  # complex64 (azimuth, range)
```

`BiomassSlc` conforms to the `SLC` protocol used across sarlet/sweets
(`radar_grid`, `orbit`, `doppler`, `wavelength`, `shape`, `bounds`, ...), so it
drops straight into `isce3.geocode.geocode_slc`. See
[`examples/geocode_slc.py`](examples/geocode_slc.py) for a single-scene GSLC
driver.

## Why not just extend s1-reader / a light sarlet sensor?

- A light sarlet-style sensor (e.g. `sarlet.sensors.nisar`) *wraps* something
  that already emits isce3 objects. BIOMASS has no such producer yet — that is
  exactly the gap this package fills. Once it exists, a thin
  `sarlet.sensors.biomass` / sweets adapter on top is trivial.
- We keep s1-reader's *idea* (product annotation -> isce3 objects) but not its
  API: no 70-field frozen dataclass. `BiomassSlc` holds a small parsed
  `AnnotationMetadata` and builds the isce3 objects lazily in properties.

## Key BIOMASS-vs-Sentinel-1 differences handled here

| Aspect | BIOMASS | Handling |
|---|---|---|
| Band | P-band, 435 MHz (lambda ~= 0.69 m) | `carrier_frequency` from annotation |
| Look side | **left** | `LookSide` from `antennaLookDirection` |
| Pixels | 4-band amplitude + 4-band phase GeoTIFFs (HH,HV,VH,VV) | reconstruct `amp * exp(1j*phase)` |
| Doppler | non-zero centroid (like S1) | `LUT2d` from `dcEstimate` polynomials |
| Ionosphere | dominant at P-band; correction flags and LUT shipped | flags and labeled LUT coordinates exposed; never reapplied implicitly |
| Granule | full-swath scene (not bursts) | one `BiomassSlc` per scene+pol |

## Status

**Pre-alpha, validated prototype.** Annotation units, orbit layout, complex
sample reconstruction, polarization order, Doppler construction, sarlet
protocol compatibility, stack geocoding, Dolphin, and Whirlwind have been
exercised with three BPS 4.4.2 repeat-pass products from track T007/frame F004.
See the [architecture review](docs/ARCHITECTURE.md) and
[correction inventory](docs/CORRECTIONS.md) for verified behavior and remaining
limitations. The separate [polarimetric research note](docs/POLARIMETRIC_RESEARCH.md)
lays out the moisture/vegetation hypotheses and the evidence needed to test them.
Exact commands and measured results are recorded in the
[T007/F004 validation report](docs/VALIDATION.md).

## Repeat-pass GSLC workflow

Download is an independent, optional stage. Credentials are read from a local
file and are never accepted as command-line values:

```bash
python scripts/download_maap.py \
  --bbox -66.5 2.5 -64.5 4.5 --start 2026-04-22 --end 2026-04-30 \
  --track T007 --frame F004 --dest /data/biomass --max 3 --unzip
```

Build a shared union-grid stack from existing products and optionally run
Dolphin with the Whirlwind unwrapper:

```bash
python scripts/biomass_pipeline.py \
  --products /data/biomass/BIO_S1_SCS__1S_*/ \
  --dem /data/dem_utm.tif --work /data/biomass/validation \
  --polarization HH --posting 30 30 --run-dolphin

python scripts/validate_stack.py /data/biomass/validation/gslc/*.tif \
  --output-dir /data/biomass/validation/diagnostics --looks 6
```

`--posting X Y` accepts independent UTM easting/northing postings. `--native`
chooses conservative native postings from the source GCP lattice, snapped to
the `10 / 2**n` series (the staged BIOMASS scenes resolve to 10 m x 5 m).
Dolphin diagnostic panels are written by default after a successful
`--run-dolphin`; use `--no-plot` to skip them.

Every GSLC carries NaN nodata and basic provenance tags; `stack.json` records
the common grid, DEM, source products, wavelength, flattening choice, and L1
ionosphere correction/LUT state. The default checks that the L1a ionospheric
phase-screen and group-delay corrections were already applied; it does **not**
apply the LUT a second time.

For GSLCs only (including a single source product), use:

```bash
python scripts/geocode_gslc.py \
  --products /data/biomass/BIO_S1_SCS__1S_*/ \
  --dem /data/dem_utm.tif --work /data/biomass/gslc \
  --polarization HH --posting 30 30
```

The base install is reader-only. Install optional capabilities explicitly:

```bash
pip install '.[download,ionosphere,plot]'
pixi run -e pipeline python scripts/biomass_pipeline.py --help
```

## Tests

```bash
pytest
BIOMASS_TEST_DATA=/data/biomass pytest
```

The second form enables windowed regression tests against local full products;
large ESA data are not copied into the repository.

## Install

```bash
pixi install
pixi run biomass-reader --help
```

## Data access

BIOMASS L1a SCS products are distributed via the **ESA-MAAP** STAC catalog
(`https://catalog.maap.eo.esa.int/catalogue/`), collection `BiomassLevel1a`,
`productType='S1_SCS__1S'`. Catalog *search* is public; *download* needs a MAAP
Bearer token (offline token from the MAAP portal). This is a different identity
system from `earth.esa.int/eogateway`.

## License

BSD-3-Clause OR Apache-2.0.
