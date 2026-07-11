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
| Ionosphere | dominant at P-band; `ionosphereCorrection` LUT shipped | surfaced via `biomass_reader.ionosphere` |
| Granule | full-swath scene (not bursts) | one `BiomassSlc` per scene+pol |

## Status

**Pre-alpha scaffold.** The annotation/orbit element paths are derived from the
official BIOMASS Processing Suite XSDs
(`bio-l1ab-main-annotation.xsd`, `bio-l1-annotations.xsd`) and parsed
namespace-agnostically. They still need to be **validated against a real
downloaded SCS product** — see `VALIDATION` checklist below. The Doppler-LUT
grid construction and orbit-file layout in particular should be confirmed
end-to-end by geocoding one scene.

### VALIDATION checklist (do once first product is on disk)

- [ ] Confirm annotation element names/units match a real `annotation/*.xml`
      (esp. `firstSampleSlantRangeTime` units, `radarCarrierFrequency`,
      `antennaLookDirection` value).
- [ ] Confirm orbit file layout (`OSV` / `X,Y,Z,VX,VY,VZ`, `UTC=` prefix).
- [ ] Confirm measurement band order is `HH, HV, VH, VV`.
- [ ] Geocode one scene (`examples/geocode_slc.py`) and eyeball against the
      product `preview/` quicklook + KML footprint.

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
