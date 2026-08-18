# biomass-reader

An [ISCE3](https://github.com/isce-framework/isce3) reader for [**ESA BIOMASS Level-1a SCS**](https://earth.esa.int/eogateway/missions/biomass) (single-look complex, slant-range) products.

BIOMASS L1a SCS is a Sentinel-1-like SAFE product with metadata formatted as XML annotation files.
This package takes a similar approach to [`s1-reader`](https://github.com/isce-framework/s1-reader) (parsing Sentinel-1 metadata for use in ISCE3).
This allows the creation of Geocoded SLCS (GSLCS) readable by [dolphin](https://github.com/isce-framework/dolphin) /  [sweets](https://github.com/isce-framework/sweets), or for simplifying interferogram creation.

## Example Usage

```python
from biomass_reader import BiomassSlc

slc = BiomassSlc.from_dir("BIO_S1_SCS__1S_..._DU1SS4", polarization="HH")

radar_grid = slc.radar_grid    # isce3.product.RadarGridParameters
orbit      = slc.orbit         # isce3.core.Orbit
doppler    = slc.doppler       # isce3.core.LUT2d  (native centroid, non-zero)
data       = slc.read_complex()  # complex64 (azimuth, range)
```

`BiomassSlc` creates an  `SLC` protocol with attributes like `radar_grid`, `orbit`, `doppler`, `wavelength`, `shape`, `bounds`.
This can be fed into [`isce3.geocode.geocode_slc`](https://isce-framework.github.io/isce3/api/python/isce3/geocode/geocode_slc.html).
See [`examples/geocode_slc.py`](examples/geocode_slc.py) for a single-scene GSLC driver.

## BIOMASS-vs-Sentinel-1 differences handled here

| Aspect     | BIOMASS                                           | Handling                                                    |
| ---------- | ------------------------------------------------- | ----------------------------------------------------------- |
| Band       | P-band, 435 MHz (lambda ~= 0.69 m)                | `carrier_frequency` from annotation                         |
| Look side  | **left**                                          | `LookSide` from `antennaLookDirection`                      |
| Doppler    | non-zero centroid (like S1)                       | `LUT2d` from `dcEstimate` polynomials                       |
| Ionosphere | strong at P-band, but provided as correction LUTs | flags and labeled LUT coordinates exposed but not reapplied |
| Granule    | full-swath scene, as opposed to S1 bursts         | one `BiomassSlc` per scene+pol                              |

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

`--posting X Y` accepts independent UTM easting/northing postings. `--native` chooses conservative native postings from the source GCP lattice, snapped to the `10 / 2**n` series (the staged BIOMASS scenes resolve to 10 m x 5 m).

Dolphin diagnostic panels are written by default after a successful `--run-dolphin`; use `--no-plot` to skip them.

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

## Reader Status

Pre-alpha prototype.
Annotation units, orbit layout, complex sample reconstruction, polarization order, Doppler construction, stack geocoding, Dolphin, and Whirlwind have been exercised with three BPS 4.4.2 repeat-pass products from track T007/frame F004.

## Tests

```bash
pytest
BIOMASS_TEST_DATA=/data/biomass pytest
```

The second form enables windowed regression tests against local full products; large ESA data are not copied into the repository.

## Install

```bash
pixi install
pixi run biomass-reader --help
```

## Data access

BIOMASS L1a SCS products are distributed via the **ESA-MAAP** STAC catalog (`https://catalog.maap.eo.esa.int/catalogue/`), collection `BiomassLevel1a`, `productType='S1_SCS__1S'`. Catalog search is public but download needs a MAAP Bearer token (offline token from the MAAP portal). 

## License

BSD-3-Clause OR Apache-2.0.
