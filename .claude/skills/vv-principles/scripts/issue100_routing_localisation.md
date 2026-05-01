# Issue #100 — Falsification is not localisation (Class-B normalisation routing)

**The bug.** The rank-N closure for Class-B (sphere) hollow
geometries routed mode-0 through the legacy `compute_P_esc` (no
surface Jacobian, "isotropic escape probability") while modes
n ≥ 1 used the canonical `compute_P_esc_mode` (with the
`(ρ_max/R)²` surface-to-observer Jacobian). The two primitives
live in *different normalisation spaces*. The mismatch is
partially absorbed in 1-region geometry where the rank-1 Mark
closure is calibrated to the legacy mode-0; in multi-region with
a strong-scatterer outer shell, the calibration breaks. Sphere
1G/2R rank-2 came in at +57 % k_eff; rank-3..8 plateaued at +66 %
rather than converging.

**Evidence that existed before it was caught.** Rank-1 Mark agreed
with reference to acceptable precision on every tested geometry.
The rank-N Marshak machinery passed every internal consistency
check — mode-0 alone was correct; mode-n alone was correct; the
rank-N operator was symmetric, the closure matrix invertible, the
eigenvalue iteration convergent. Each component, tested in
isolation, was verified.

**Why that evidence didn't catch it.** Single-mode tests cannot
expose mode-routing bugs — by definition, only one routing path
is exercised. Monte Carlo parity revealed the *symptom* (wrong
k_eff) but could not localise the *cause*: an external reference
tells you the answer is wrong, not which term is wrong. The
investigation required an internal algebraic invariant — the
Schur-equivalence identity (the rank-N closure must reduce to the
G·R·P matrix at unit rank, and to the Mark closure on a
homogeneous infinite medium) — to localise the defect to the
mode-0 vs mode-n routing fork. The 8-probe cascade narrowed the
bug only when Probe G replaced the legacy mode-0 with the
canonical form and observed the systematic shift propagate
through every multi-region case.

**What evidence class would have caught it.** When an
investigation must distinguish "the answer is wrong" from "this
specific term is wrong," external falsification (MC, analytical
reference) is necessary but **NEVER** sufficient — it provides a
discrepancy, not a location. **Instead**, the diagnostic chain
**MUST** include an *internal algebraic invariant* that
constrains the structure of the operator: Schur equivalence
under rank reduction, conservation under specific limits,
reduction-to-known-closure on degenerate geometry. Single-mode
tests must be augmented by mode-coupling tests (rank-N
parametrised over rank, with the rank-1 reduction asserted as
an invariant). Without an algebraic-invariant probe, MC
localisation is structurally insufficient for closure-routing
bugs. See vv-principles reference.md §1 (structural independence
applied to localisation, not just falsification).

**References.** Issue #100 (open); Issue #103 (open); Issue #132
(open follow-up); numerics-investigator agent memory
`issue_100_class_b_mr_mg.md`; probe cascade
`scratch/derivations/diagnostics/diag_class_b_rank_n_probe_*.py`
(probes B–H).
