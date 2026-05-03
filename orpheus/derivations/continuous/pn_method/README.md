# `pn_method/` — Spherical-harmonics expansion (P_N method) — reserved

**Status:** empty / placeholder. Implementation deferred.

## Intent

Modal expansion of the angular flux in Legendre / spherical-harmonic
polynomials. The historical alternative to the discrete-ordinates
S_N family. Multi-region sphere k_eff and flux shapes — the natural
cross-check target for `trajectory_resolvent/` multi-region sphere
results that ORPHEUS already verifies against Garcia 2021 truth values.

Having `pn_method/` reserved as a folder (even empty) makes the
reference asymmetry explicit: today ORPHEUS uses Garcia 2021 as a
*truth set without a method-of-record*; once this folder is populated
the cross-check becomes a structurally-independent two-method
agreement.

## Canonical references (LOCAL in `scratch/literature/`)

* **Garcia, R.D.M., Siewert, C.E. (2017).** P_N method for the sphere.
  Local PDF in `scratch/literature/`.
* **Garcia, R.D.M. et al. (2019).** Spherical-shell P_N.
  Local PDF in `scratch/literature/`.
* **Garcia, R.D.M. et al. (2021).** Multi-region sphere P_N. The
  primary truth set ORPHEUS already cross-checks against in
  `tests/derivations/test_trajectory_resolvent_garcia2021.py`.
  Local PDF in `scratch/literature/`.
* **Davison, B. (1957).** *Neutron Transport Theory*. Oxford, chs 9–11.
  The textbook foundational reference for the P_N method.
