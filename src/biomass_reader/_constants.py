"""Physical and mission constants for BIOMASS."""

from __future__ import annotations

# Speed of light in vacuum [m/s].
SPEED_OF_LIGHT = 299_792_458.0

# BIOMASS P-band carrier: 435 MHz, 6 MHz bandwidth, fully polarimetric.
# The exact carrier is read per-product from `instrumentParameters`; this is a
# fallback / sanity-check value only.
NOMINAL_CARRIER_FREQUENCY_HZ = 435e6

# Measurement band order in the L1a SCS abs/phase GeoTIFFs.
POLARIZATION_ORDER = ("HH", "HV", "VH", "VV")
