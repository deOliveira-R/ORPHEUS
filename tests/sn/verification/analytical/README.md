# SN L1 analytical-reference test suite

L1 verification tests that consume **structurally-independent
analytical references** from `orpheus/derivations/continuous/` and
assert the SN solver reproduces them.

This directory replaces the bulk of full-solver eigenvalue tests in
the legacy SN test files. The migration was sequenced per the SN
test-architecture matrix (a since-retired design doc; the tagging
conventions now live in `docs/theory/verification/harness.rst`)
and `.claude/plans/sn_reshape.md` Phase 2.

## Why analytical-reference

Each test here imports a `ContinuousReferenceSolution` from
`orpheus.derivations.continuous.analytical.homogeneous` (1G/2G/4G
homogeneous medium k_inf and spectrum), or from one of the other
production research-grade reference modules:

* `orpheus.derivations.continuous.fn_method` — slab / sphere F_N
  expansion (closed-form analytic continuation).
* `orpheus.derivations.continuous.singular_eigenfunction` — Case
  singular-eigenfunction expansion for slab.
* `orpheus.derivations.continuous.trajectory_resolvent` — Sanchez
  2002 closed-form MoC (research-grade accuracy for heterogeneous /
  curvilinear).

The references are computed by structurally-orthogonal mathematics
(linear algebra of cross-section matrices for homogeneous; spectral
expansion for F_N; singular-integral methods for Case;
trajectory-domain closed-form for `trajectory_resolvent`). A passing
test cross-checks the SN discretisation chain end-to-end against a
mathematically-independent ground.

## Memory + runtime budget

Each test runs **one** SN solver invocation on a small mesh
(typically ≤ 20 cells, 1 group ≤ 4) with claim-driven quadrature
(GL-1D N=8 for slab/sphere; ProductQuadrature(2x4) for cylinder).
Per-test peak ψ-array footprint is well under 100 MB; runtimes are
<1 s. The tier-1 default (`pytest -m "(l0 or l1) and not slow"`)
fits comfortably in any developer environment.

## Markers

* `pytestmark = pytest.mark.l1` — every file marked at module scope.
* `@pytest.mark.verifies("<label>")` — links to Sphinx equation
  labels in `docs/theory/`. Nexus consumes these to build
  test ↔ equation edges.
* No `@pytest.mark.slow` — anything heavy enough to want `slow`
  belongs in `tests/sn/l2_integration/` (planned for after the
  reshape).
