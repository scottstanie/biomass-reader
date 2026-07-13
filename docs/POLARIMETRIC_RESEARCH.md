# Separating moisture and vegetation effects: a testable phase program

This is a research plan, not a correction algorithm. Three acquisitions are
enough to expose useful observables, but not enough to identify soil moisture,
vegetation structure, motion, and propagation independently.

## Measurement model

For polarization (p) and acquisition (t), write the focused complex return
as a coherent sum of surface, double-bounce, and vegetation contributions:

\[
S_p(t)=\sum_k A_{p,k}(t)\exp\{i\phi_{p,k}(t)\}.
\]

The multilooked temporal interferogram is

\[
I_p(t_a,t_b)=\langle S_p(t_a)S_p^*(t_b)\rangle.
\]

Its phase contains acquisition geometry and displacement, residual propagation,
and the polarization-dependent change in the effective scattering phase center.
A co-polar differential observable,

\[
D_{HH,VV}^{ab}=\arg\{I_{HH}^{ab}(I_{VV}^{ab})^*\},
\]

cancels phase terms that are genuinely common to HH and VV. It does **not**
automatically isolate moisture: different mixtures of surface, double-bounce,
and volume scattering can move the HH and VV phase centers differently.
Polarization-dependent residual ionosphere and calibration errors are also
possible, even though the tested L1a products report that Faraday rotation,
phase screen, and group delay were corrected.

For three dates, the closure phase

\[
C_p^{123}=\arg\{I_p^{12}I_p^{23}(I_p^{13})^*\}
\]

cancels phase that can be represented as a scalar per acquisition. Non-zero
closure therefore diagnoses non-separable scattering evolution and estimator
noise; it is not uniquely a deformation, soil-moisture, or vegetation signal.
The soil-scattering model of [De Zan et al. (2014)](https://doi.org/10.1109/TGRS.2013.2241069)
specifically predicts that moisture-dependent subsurface scattering can break
triplet phase consistency. Classical PolInSAR provides the broader framework
for polarization-dependent phase centers
([Cloude and Papathanassiou, 1998](https://doi.org/10.1109/36.718859)).

## First T007/F004 observations

The exploratory script uses 6 by 6 spatial looks on the 30 m GSLCs and gates
phase summaries at coherence 0.5.

| Observable | 2026-04-23 | 2026-04-26 | 2026-04-29 |
|---|---:|---:|---:|
| median HV/VH coherence | 0.9909 | 0.9911 | 0.9909 |
| median absolute HV-VH phase | 0.0222 rad | 0.0215 rad | 0.0209 rad |
| median HH/VV coherence | 0.4140 | 0.4150 | 0.4152 |
| circular mean HH-VV phase | -0.0660 rad | -0.0617 rad | -0.0644 rad |

HV/VH agreement is a strong reciprocity and processing sanity check. HH/VV is
less coherent, as expected when the channels weight scattering mechanisms
differently, but its spatial circular mean is stable over these dates.

For the 23–26, 23–29, and 26–29 April pairs, the median absolute high-coherence
HH–VV temporal differential phases are 0.176, 0.193, and 0.178 rad. Median
absolute closure phase is 0.095–0.102 rad across the four channels. These values
are large enough to study but are not source attribution.

The result is consistent with published warnings that polarization-dependent
interferometric scattering phase can be driven by both soil-moisture structure
and vegetation changes. Airborne observations found agricultural phase
diversity related to wet biomass in some cases rather than soil moisture
([Brancato et al., 2021](https://doi.org/10.1029/2020EA001445)). P-band
polarimetric work likewise separates surface, double-bounce, and vegetation
volume contributions rather than assigning a single channel to a single cause
([Alemohammad et al., 2018](https://doi.org/10.1016/j.rse.2018.02.032)).

## Proposed discrimination experiment

1. Build at least a season-long time series. Three dates provide one closure
   triplet but cannot establish repeatable response functions.
2. Form the calibrated scattering vector
   \([S_{HH},\sqrt{2}S_{HV},S_{VV}]^T\), after verifying HV/VH reciprocity.
   Estimate covariance/coherency matrices using a physically stated window.
3. Stratify pixels by polarimetric mechanism (surface, double bounce, volume),
   forest/non-forest class, incidence angle, slope, and temporal coherence.
4. Test temporal HH–VV differential phase and per-channel closure against rain,
   modeled/in-situ soil moisture, vegetation water content, wind, and phenology.
5. Fit nested models rather than a single phase threshold:

   - acquisition effects: residual common and polarization-dependent ramps;
   - surface state: moisture/dielectric and roughness interaction;
   - vegetation state: water content, motion, and structural change;
   - stable spatial effects: topography, incidence angle, and scattering class.

6. Use held-out dates and sites. A moisture interpretation is credible only if
   it follows independent moisture observations within surface-dominated pixels
   and transfers across dates after rainfall, vegetation, and geometry controls.
7. Treat closure as a model-checking observable. If a proposed acquisition-wise
   correction removes pair phase but leaves structured closure, its scattering
   model is incomplete.

The companion `scripts/polarimetric_phase.py` writes the observables and browse
image. It intentionally performs no inversion and labels differential phases as
diagnostics rather than moisture products.
