# `galerkin_sn_hybrid/` — Hybrid collocation-Galerkin-S_N — reserved (lower priority)

**Status:** empty / placeholder. Implementation deferred. Lower
priority — research-grade with narrow application.

## Intent

Morel's collocation-Galerkin-S_N angular discretisation: a hybrid that
sits between modal expansion (P_N) and discrete-ordinates (S_N),
keeping discrete-ordinate sweeps as the spatial mechanism while using
Galerkin projection in the angular variable. Useful for problems where
neither P_N nor S_N alone behaves well (highly anisotropic scattering
combined with strong streaming, where ray effects in S_N and ringing
in P_N both bite).

## Canonical reference (LOCAL in `scratch/literature/`)

* **Morel, J.E. (1989).** "A Hybrid Collocation-Galerkin-Sn Method for
  Solving the Boltzmann Transport Equation."
  Local PDF: `scratch/literature/Morel1989.pdf`.
