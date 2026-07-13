# Architecture and correctness review

## Scope

`biomass-reader` is a sensor adapter, not a deformation-processing framework.
It resolves one BIOMASS L1a SCS product, reconstructs its complex scattering
channels, and presents the geometry as isce3 and sarlet objects. Download,
stack geocoding, Dolphin, and scientific experiments remain thin workflows on
top.

## Public surface

- `BiomassProduct`: deterministic paths to MPH, annotation, navigation, LUT,
  complex VRT, and source amplitude/phase rasters.
- `BiomassSlc`: one polarization implementing sarlet's runtime `SLC` protocol.
- `parse_annotation`: the minimal radar-grid and Doppler metadata.
- `parse_correction_status`: correction flags, separate from radar geometry.
- `open_lut_group`: labeled access to auxiliary LUTs with their coordinates.
- `make_shared_geogrid` and `geocode_slc`: phase-preserving stack primitives.

The product-supplied complex VRT is the canonical `slc_path`. Direct
amplitude/phase reconstruction is retained for efficient windowed reads and is
regression-tested against all four VRT bands.

## Verified assumptions

The three staged T007/F004 products verify:

- shape 20,767 azimuth lines by 1,367 range samples;
- left-looking geometry and 435 MHz carrier frequency;
- HH, HV, VH, VV band order;
- `amplitude * exp(1j * phase)` agrees with the BPS VRT;
- CFI orbit files produce isce3 orbits spanning sensing time;
- Doppler estimates use two-way slant-range time and geometry polynomials when
  the combined polynomial is empty;
- zero-Doppler geo2rdr reproduces the 48 product VRT GCPs with less than 0.5
  azimuth pixel and 0.001 range pixel error at the 95th percentile;
- the reader satisfies sarlet's current runtime-checkable `SLC` protocol.

## Workflow choices

The shared geogrid uses the union of every acquisition's isce3-derived extent
and snaps it to a projection-fixed pixel lattice. Intersection is available as
an explicit alternative. Outputs retain NaN nodata, CRS, transform,
polarization, wavelength, product ID, and flattening provenance.

The annotation explicitly defines image line times as zero-Doppler times, so
`image_grid_doppler` is zero during geocoding. The non-zero annotation Doppler
is still supplied as isce3's `native_doppler`; it is used for SLC baseband/phase
handling and bounds, not as the output image-grid geometry. Comparing GCPs with
the native Doppler as if it defined the image grid produces an expected apparent
azimuth offset and is not a timing correction.

The DEM should cover the valid scene footprint plus geo2rdr height margin. The
current demonstration DEM does not cover the entire conservative union bounding
box, so isce3 reports a western-limit warning; those output pixels remain
invalid. This is visible rather than converted into plausible zero-valued data.

## Remaining limitations

- Attitude is resolved but not used by the current isce3 primitive.
- XML lookup is namespace-independent, but processor-baseline variation beyond
  the tested BPS 4.4.2 products needs additional fixtures.
- No residual ionosphere estimator is implemented.
- A compact redistributable real-data fixture is still needed for CI; local
  full-product tests are enabled through `BIOMASS_TEST_DATA`.
- GSLC is currently a complex GeoTIFF plus provenance JSON, not a formal ESA or
  NISAR product specification.
