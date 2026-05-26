---
name: phase1-verification-plan
description: Verification plan for Phase 1 of the moment-space-and-layering refactor (ERR-039 fix via SphericalHarmonicBasis / SphericalHarmonicSpace / MomentProjection / .T-vs-.H correction). Pin tests for legacy contracts, new tests for the new conventions, L1 MMS gates that must stay green, P2/P3 forward-compatibility, P3.1 import-linter spec, and an ordered gates checklist.
metadata:
  type: project
---

# Phase 1 — Verification Plan (Moment-Space & Layering)

**Branch**: `refactor/moment-space-and-layering` (worktree:
`.claude/worktrees/moment-space-and-layering/`)
**Plan**: `.claude/plans/moment_space_and_layering_plan.md`
**ERR**: ERR-039 (root cause); plan §1 enumerates the four operators
(`S₀`, `Πᵀ`, `Π.H`, `R`) collapsed into the `(2ℓ+1)` literal.

## Scope and ordering

Phase 1 has steps P1.0 (recon, done) → P1.1 (basis) → P1.2 (space) →
P1.3 (projection rewire) → P1.4 (`.T` / `.H` correction) → P1.5 (new
V&V tests) → P1.6 (docstring + retire `assert_galerkin_idempotency`) →
P1.7 (retire `sn/solver.py` inline `(2ℓ+1)`). This plan ladders the
test work to those steps. Every step **MUST** preserve the gates from
all prior steps; the "Gates checklist" at the end gives the exact
invocation sequence the implementer runs after each P1.x.

## CRITICAL — claim layer and pillar gating (per `vv-principles`)

Per `vv-principles` §"Hierarchical claim taxonomy" and §"The three
pillars of verification", every row below declares its claim layer and
pillar. Phase 1 is a representation-and-adjoint refactor — it touches
**no eigenvalue logic, no flux-shape solution, no spatial
discretisation**. The claim layers are therefore:

- **Algebraic-identity claims** (the equation `R = (2ℓ+1)·S₀`,
  `Π·R = 4π·I`, `Π.H = g_C·S₀`, etc.) — verified against
  **closed-form** SymPy / hand-derived identities. Pillar: closed-form
  (Branch 1A per `algebra-of-record`). Structural independence:
  the references in the V&V tests are derived from the SH addition
  theorem and the Galerkin adjoint identity — both upstream of any
  ORPHEUS code path under test.
- **Bit-identity invariance claims** (the SN production path remains
  bit-identical, or drifts only within the FP-non-associativity bound
  per `vv-principles` §"Bit-identity vs principled-equivalence") —
  verified against frozen snapshots (`tests/sn/regression/snapshots/`)
  and the L1 MMS gates which carry analytical references.
- **NO eigenvalue claims are made by Phase 1.** The aniso scattering
  source builder (`ScatteringOperator.build_aniso_source`) does not
  change — its consumers (L1 MMS-aniso, eigenvalue runs in the
  regression snapshots) keep their existing references.

MMS gates (P1.7-affected via `_build_rhs_cartesian` consumers AND the
L1 anisotropic scattering tests) are convergence-order + flux-shape
claims — pillar: MMS (Branch 1C). They do NOT prove eigenvalues; their
purpose here is to gate the operator's spatial / angular discretisation
through the rewire.

Cardinal-rule check (`vv-principles` §"1-group degeneracy"): every
multi-group test below includes ≥2G. The L1 MMS-aniso gate uses 1G but
is a flux-shape MMS claim, NOT an eigenvalue claim — degeneracy does
not apply.

---

## Section A — Legacy tests that pin behaviour through the refactor

The four legacy test surfaces are (1) the projection-operator unit
tests, (2) the SH evaluation tests, (3) the SN regression snapshots,
(4) the SN scattering tests. Their behaviour through P1.1–P1.7 is
specified below at file:line granularity.

### A.1 `tests/numerics/test_projection_operators.py` (410 lines)

Verified against the worktree copy as of 62994ad. Class-by-class
disposition:

| Lines | Class / Test | Disposition through Phase 1 |
|---|---|---|
| 47-57 | `TestABCs` (3 tests) | **LEGACY, ADAPT**. P1.3 renames `HarmonicMomentProjection → MomentProjection`. Imports + class names update; assertions unchanged (subclassing remains). |
| 65-117 | `TestHarmonicMomentProjectionShapes` (3 tests) | **LEGACY, ADAPT names** only. The `apply` einsum (line 340 of `projection.py`) is unchanged by P1.3/P1.4 — only `apply_transpose` changes. Assertions stand bit-identical. |
| 112-117 | `TestHarmonicMomentProjectionShapeValidation` (1 test) | **LEGACY, ADAPT names**. The `__post_init__` validator is unchanged. |
| 121-136 | `TestHarmonicMomentProjectionFromMeasure` (2 tests) | **LEGACY, ADAPT names**. P1.3 keeps `from_measure` working bit-identical (the constructor classmethod is unchanged). The new `from_spherical_harmonic_space` classmethod is exercised by NEW tests in §B. |
| 144-175 | `TestHarmonicMomentReconstructionShapes` (3 tests) | **LEGACY, PRESERVED bit-identical**. P1.3 keeps `from_Y` as a back-compat shim (plan §P1.3); the `apply` einsum (line 475 of `projection.py`) is **deliberately not modified** (Anti-recommendation 3, plan §2.1). Test at line 173 (`two_l_plus_one = 2.0 * np.arange(L + 1) + 1.0`) is the spec; it stays. |
| **194-231** | `TestGalerkinIdempotencyOnLebedev.test_pi_R_is_identity_on_band_limited` (the `Π·R = 4π·I` identity, parametrised over `(order=7,L=2)`, `(13,3)`, `(17,4)`) | **LEGACY, CRITICAL — MUST CONTINUE PASSING bit-identical.** This is the genuine Galerkin identity; the P1.5 sibling will assert it again via the new `SphericalHarmonicSpace` API but this test pins the existing (Π, R) pair. Plan §P1.6 explicitly preserves it. |
| **234-302** | `TestApplyTransposeIsWWeightedAdjoint` (`@catches("ERR-039")`, 2 tests at 248-271 and 273-302) | **LEGACY, REWRITE under P1.4.** Currently asserts `apply_transpose` returns bare `S₀` (no `w_n`). After P1.4, `apply_transpose` (or its `.T` successor) returns `Πᵀ = w_n · S₀`. The W-weighted adjoint identity in the first test (LHS `⟨Πψ,c⟩_C` line 264, RHS `⟨ψ, Π* c⟩_V` line 269) MUST still hold — but the RHS computation MUST drop the EXTERNAL `measure.weights` factor that line 269 currently applies, because the operator now CARRIES `w_n` internally. **Specific edits:** |
| | | Line 268: change `M.apply_transpose(c_masked)` (a) keep call site under the new `.T` semantics; the returned array now equals `w_n · S₀(c_masked)` (shape `(N,)`); (b) line 269 RHS becomes `float(np.sum(psi * Pi_T_c))` — drop the `measure.weights *` factor since `Pi_T_c` already carries it. The identity LHS=RHS still pins ERR-039. |
| | | Lines 273-302 (`test_apply_transpose_no_2l_plus_1_factor`): semantics shift from "bare `S₀`" to "`w_n · S₀`". The "differs from `R`" assertion still holds (`R = (2ℓ+1)·S₀` ≠ `w_n · S₀`). The expected RHS at line 299 changes from `np.einsum("nlm,lm->n", M.Y, c_masked)` to `np.einsum("n,nlm,lm->n", M.weights, M.Y, c_masked)`. Test name updates to `test_T_returns_w_weighted_transpose`. |
| | | Decorator: keep `@pytest.mark.catches("ERR-039")` + `@pytest.mark.l1`. The test class moves to NEW file `tests/numerics/test_spherical_harmonic_space.py` (P1.5) as the canonical ERR-039 catching site, OR the rewrite stays here and the new file adds the broader space-based coverage. **Recommendation**: rewrite in-place under the new name; the new P1.5 file covers the space-level invariants. |
| 305-350 | `TestGalerkinAdjointPairing.test_adjoint_pairing_matches` | **LEGACY, ADAPT.** Same shift as above: the RHS at line 349 currently writes `np.sum(psi * measure.weights * Rc_no_factor)` because the prose constructs `Rc_no_factor = np.einsum("nlm,lm->n", Y, c_masked)` — i.e. the OLD `apply_transpose` semantics. Under P1.4 this is `M.T.apply(c_masked)` (which now includes `w_n`), so drop the explicit `measure.weights *` from line 349. The mathematics is invariant; only the wiring updates. |
| 358-381 | `TestAssertGalerkinIdempotencyMethod.test_method_signals_violation` | **DELETED** under P1.6 (Resolution #4 in plan §0.2; CC.5 in plan §1.5). The method `assert_galerkin_idempotency` is itself retired; the test exists only to exercise it. Both go in the same P1.6 commit per `lessons-L20` (retirement = test migration). |
| 389-409 | `TestCapabilities` (2 tests) + `test_petrov_galerkin_is_abstract_no_ship_concrete` | **LEGACY, ADAPT names**. `CAP_APPLY_TRANSPOSE` is still advertised by `MomentProjection` (the operator now correctly carries `w_n`); the bare `from_Y` reconstruction still advertises `CAP_APPLY` only. |

**Aggregate disposition**:
- Pin (bit-identical): TestABCs (renames), TestHarmonicMomentProjectionShapes, TestHarmonicMomentProjectionShapeValidation, TestHarmonicMomentProjectionFromMeasure, TestHarmonicMomentReconstructionShapes, TestGalerkinIdempotencyOnLebedev, TestCapabilities, test_petrov_galerkin_is_abstract_no_ship_concrete.
- Rewrite (semantics shift under P1.4): TestApplyTransposeIsWWeightedAdjoint (2 tests), TestGalerkinAdjointPairing (1 test).
- Delete (retirement under P1.6): TestAssertGalerkinIdempotencyMethod (1 test).

### A.2 `tests/numerics/test_spherical_harmonics.py`

The function `evaluate_real_sh` relocates from
`orpheus/numerics/spherical_harmonics.py` to a **method** on
`SphericalHarmonicBasis` (plan §P1.1). The relocation is mechanical;
the math is identical.

| Test class | Disposition |
|---|---|
| `TestEvaluateRealSHBasic` (L0, lines 30-54) | **LEGACY, IMPORT-PATH CHANGE ONLY**. P1.1 ships a back-compat shim at `orpheus/numerics/spherical_harmonics.py` that re-exports `evaluate_real_sh` (plan §P1.1: "thin re-export shim of `SphericalHarmonicBasis.evaluate`"). The existing imports (`from orpheus.numerics.spherical_harmonics import evaluate_real_sh`) keep working through Phase 1; the shim deletes in P3.2 along with these imports. |
| `TestAdditionTheorem` (L1, lines 57+) | **LEGACY, IMPORT-PATH CHANGE ONLY**. Same shim story. |

**Required action in P1.1**: re-verify all imports of `evaluate_real_sh` continue working — the four sn-quadrature delegators (`numerics/quadrature/directional.py:75, 393` plus the four sn-quadrature wrappers per plan §P1.1) AND these two test classes.

### A.3 `tests/sn/regression/snapshots/` (15+ .npz files)

The frozen regression suite is the **bit-identity contract** for the
production SN path. Phase 1 changes ONLY (a) the projection's
`apply_transpose` semantics (P1.4) — zero production callers per memo
§A1.a; (b) the SH evaluation import path (P1.1) — production path
unchanged; (c) the `_build_rhs_cartesian` inline `(2ℓ+1)` migration
(P1.7) — production path computational changes possible.

| Step | Bit-identity contract |
|---|---|
| P1.1–P1.6 | **STRICT bit-identity required.** Slab cases `slab_*` keep `rtol=1e-12, atol=1e-13`; curvilinear cases `sphere_*`, `cyl_*` keep `rtol=5e-6, atol=1e-7` (the existing iteration-floor relaxation; see `test_dd_regression.py:78-80`). Reason: production code paths (the `ScatteringOperator.build_aniso_source` pipeline at `sn/scattering.py:609-657` AND the SI / Krylov outer loop) are not touched in P1.1–P1.6. P1.4 affects ONLY `apply_transpose`, which has zero production callers (plan §P1.4). |
| **P1.7** | **PRINCIPLED-EQUIVALENCE acceptable per `vv-principles` §"Bit-identity vs principled-equivalence"** (the three criteria). The migration replaces the hand-rolled triple-loop `(2*l+1) * (sig_s[mid][l].T @ SUM) / sum_w` (line 930) with `ReconstructionOperator @ MomentProjection @ ...` einsum chain. The reduction tree changes; ULP drift expected. The three criteria below MUST be satisfied for the drift to be accepted: |

**P1.7 bit-identity criteria (per `vv-principles` §"Bit-identity vs principled-equivalence")**:

1. **Principled at every step.** The new chain `R @ Λ @ M @ ψ` produces named intermediates: `moments = M(ψ)` is the harmonic-moment field, `scattered = Λ(moments)` is the per-ℓ-block-diagonal scaling output, `q_aniso = R(scattered) / W` is the per-ordinate aniso source. Each is a reactor-physics quantity (HarmonicMomentField; cf. `sn/harmonic_moment_field.py`). The legacy `fiL[ix, iy, :, l, l+m]` array and `qS` inner-loop accumulator at `sn/solver.py:887-930` are unnamed scratch intermediates. **Criterion 1: SATISFIED** by the migration's nature.
2. **Structurally-independent reference**: each affected snapshot must agree with the analytical limit where one exists. The homogeneous-reflective `*_homogeneous_*` snapshots converge to `k_inf = νΣ_f/Σ_a` exactly (cited in `test_dd_regression.py:46`); the heterogeneous snapshots are pinned by the L1 MMS-aniso gate (`tests/sn/test_mms_aniso.py`) which uses MMS as its independent reference. **Action**: re-run the snapshot generator post-P1.7 and verify each `*_homogeneous_*` snapshot's `keff` still matches `k_inf` to within `outer_iters × ULP`; verify each heterogeneous snapshot's flux remains within `rtol=5e-6` (the existing curvilinear bound — slabs may stay at `rtol=1e-12` if the drift is truly FP-only). **Criterion 2: PIN VIA EXISTING REFERENCES** — no new references needed.
3. **Drift dimensionally explainable**: the new chain has reduction depth `O(N · (L+1)²)` (one einsum) vs the legacy triple loop's `O(N · (L+1)² · ng)` reduction. Drift bound per step 3 of the criteria: `reduction_depth × ULP × condition_number ≤ outer_iters × ULP_floor`. For slabs (well-conditioned), drift should remain ≤ 1e-13; for curvilinear (worse-conditioned), drift should remain ≤ existing 5e-6 floor. **Action**: if any snapshot exceeds these bounds post-P1.7, investigate as a real bug (not a tolerance drift) — the bound is the contract.

**P1.7 snapshot relaxation policy**: do NOT loosen `rtol` or `atol` in `test_dd_regression.py`. If the drift exceeds the existing bound, the migration is wrong. If it stays within, no change needed. Document the verification in the P1.7 commit message: "Verified post-migration drift ≤ existing snapshot bound per `vv-principles` §Bit-identity criterion 3."

### A.4 `tests/sn/test_scattering_operator.py` and related

The single production consumer of the `(M, Λ, R)` triple is
`ScatteringOperator.build_aniso_source` (`sn/scattering.py:525-657`).
Phase 1 changes **only** the construction site of `R` — its `apply`
einsum is bit-identical (Anti-recommendation 3). Therefore:

| Test file | Disposition |
|---|---|
| `tests/sn/test_scattering_operator.py` | **STRICT bit-identical** through P1.1–P1.6. Production `build_aniso_source` calls `HarmonicMomentReconstruction.from_Y(Y)` (line 633); this constructor (line 461 of `projection.py`) is unchanged by P1.3 (kept as back-compat shim). Even if the new `from_spherical_harmonic_space` classmethod sources `two_l_plus_one` from the space, the values are identical and the einsum is identical. |
| `tests/sn/test_legendre_moment_scattering.py` | **STRICT bit-identical** (Λ operator is untouched by Phase 1). |
| `tests/sn/test_operators_apply_typed.py` | **STRICT bit-identical** (FD-matvec / probe-tests; bypass the type layer entirely per `sn/scattering.py:652-654`). |
| `tests/sn/test_harmonic_moment_field.py` | **STRICT bit-identical** (P1 does not touch `HarmonicMomentField`'s typed-field algebra; relocation happens in P3.3 not P1). |

---

## Section B — New tests under P1.5

Test file: `tests/numerics/test_spherical_harmonic_space.py`. ALL tests
`@pytest.mark.catches("ERR-039")` + `@pytest.mark.l1` +
`@pytest.mark.verifies(<label>)` where the label resolves to a
`.. math:: :label:` block added to
`docs/theory/spherical_harmonics.rst` in P1.6.

### B.1 The five plan-§P1.5 identities

The plan §P1.5 enumerates five identities. Each becomes a single
pytest function. Tolerances chosen per the existing
`test_projection_operators.py` precedent (rtol=1e-10 to 1e-15
depending on whether quadrature exactness or einsum-direct).

```python
"""Tests for SphericalHarmonicSpace and the metric / adjoint algebra.

ERR-039 root-cause invariants: the four operators (S0, Pi^T, Pi.H, R)
collapsed into the (2l+1) literal are now distinct typed objects
with executable identities.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.basis import SphericalHarmonicBasis
from orpheus.numerics.spaces import SphericalHarmonicSpace
from orpheus.numerics.projection import (
    MomentProjection,
    ReconstructionOperator,  # or HarmonicMomentReconstruction
)
from orpheus.numerics.quadrature import lebedev_sphere


@pytest.fixture
def lebedev_L_pair():
    """Pair a Lebedev rule with an L that the rule integrates exactly.

    The rule of degree `deg` integrates `Y_l Y_l'` exactly for
    `l + l' <= deg`. For L=3, need deg >= 6 → order=13.
    """
    measure = lebedev_sphere(13)
    L = 3
    return measure, L


# ─────────────────────────────────────────────────────────────────────
# B.1 — the five P1.5 identities (one test each)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.catches("ERR-039")
@pytest.mark.verifies("sh-space-metric")
def test_space_inner_product_weights_equal_4pi_over_2l_plus_1():
    """SphericalHarmonicSpace.from_L(L).inner_product_weights[l] == 4pi/(2l+1).

    The metric g_C lives in exactly one place — on the space.
    """
    L = 4
    space = SphericalHarmonicSpace.from_L(L)
    expected_per_ell = 4.0 * np.pi / (2.0 * np.arange(L + 1) + 1.0)
    # The (L+1, 2L+1) array carries 4pi/(2l+1) in the 2l+1 valid slots
    # of row l and 0 in the |m|>l padding.
    weights = space.inner_product_weights
    assert weights.shape == (L + 1, 2 * L + 1)
    for ell in range(L + 1):
        # Row ell: first 2*ell+1 slots equal expected_per_ell[ell],
        # remaining slots equal 0.
        np.testing.assert_allclose(
            weights[ell, : 2 * ell + 1],
            expected_per_ell[ell],
            rtol=1e-15,
        )
        np.testing.assert_array_equal(
            weights[ell, 2 * ell + 1 :], 0.0,
        )


@pytest.mark.l1
@pytest.mark.catches("ERR-039")
@pytest.mark.verifies("sh-quadrature-gram")
def test_basis_mass_matrix_against_lebedev(lebedev_L_pair):
    """SphericalHarmonicBasis.mass_matrix(lebedev_measure) ≈ MomentMassMatrix.

    Pins the SH convention against a structurally-independent
    discrete-orthogonality computation:
        sum_n w_n Y_l^m(Omega_n) Y_l'^m'(Omega_n) ≈ (4pi/(2l+1)) δ_ll' δ_mm'
    for l + l' <= deg.
    """
    measure, L = lebedev_L_pair
    basis = SphericalHarmonicBasis(L=L)
    M_mass = basis.mass_matrix(measure)  # MomentMassMatrix
    # Diagonal entries (l, m) equal 4pi / (2l+1)
    expected_per_ell = 4.0 * np.pi / (2.0 * np.arange(L + 1) + 1.0)
    diag = M_mass.diagonal()  # whatever API the MomentMassMatrix exposes
    for ell in range(L + 1):
        np.testing.assert_allclose(
            diag[ell, : 2 * ell + 1],
            expected_per_ell[ell],
            rtol=1e-12,  # quadrature exactness + FP roundoff
        )
    # Off-diagonal (l, l') for l != l' equal 0 to quadrature precision
    # (verified via einsum on Y).
    Y = basis.evaluate(measure.nodes)  # (N, L+1, 2L+1)
    # Build the full Gram and check off-diagonal block
    G = np.einsum("n,nlm,nij->lmij", measure.weights, Y, Y)
    for ell in range(L + 1):
        for ell_p in range(L + 1):
            if ell == ell_p:
                continue
            block = G[ell, : 2 * ell + 1, ell_p, : 2 * ell_p + 1]
            np.testing.assert_allclose(block, 0.0, atol=1e-12)


@pytest.mark.l1
@pytest.mark.catches("ERR-039")
@pytest.mark.verifies("sh-addition-theorem-reconstruction")
def test_R_equals_2l_plus_1_times_S0(lebedev_L_pair):
    """R.apply(c) == (2l+1) * S_0(c) for random band-limited c.

    The addition-theorem reconstruction R is g_C^{-1} · S_0 up to 4pi.
    This pins the ERR-039 identity: R is NOT the Hilbert adjoint.
    """
    measure, L = lebedev_L_pair
    space = SphericalHarmonicSpace.from_L(L)
    basis = SphericalHarmonicBasis(L=L)
    R = ReconstructionOperator.from_spherical_harmonic_space(
        space, basis.evaluate(measure.nodes),
    )
    rng = np.random.default_rng(seed=2026)
    c = rng.standard_normal((L + 1, 2 * L + 1))
    # Mask |m| > l entries to zero
    for ell in range(L + 1):
        c[ell, 2 * ell + 1 :] = 0.0

    actual = R.apply(c)
    # S_0 is the naked synthesis (basis method, plan §P1.4)
    S0_c = basis.synthesize(c, measure.nodes)  # einsum("nlm,lm->n", Y, c)
    two_l_plus_one = 2.0 * np.arange(L + 1) + 1.0
    expected = np.einsum("nlm,l,lm->n", basis.evaluate(measure.nodes),
                         two_l_plus_one, c)
    np.testing.assert_allclose(actual, expected, rtol=1e-14)
    # Also pin: actual == (2l+1) · S_0(c), per-ell
    # (sum over m for each l of Y_lm * c_lm, times 2l+1)
    # Cross-check via direct construction
    per_ell_synthesis = np.zeros_like(actual)
    Y = basis.evaluate(measure.nodes)
    for ell in range(L + 1):
        contrib = np.einsum("nm,m->n", Y[:, ell, :], c[ell, :])
        per_ell_synthesis += (2 * ell + 1) * contrib
    np.testing.assert_allclose(actual, per_ell_synthesis, rtol=1e-14)


@pytest.mark.l1
@pytest.mark.catches("ERR-039")
@pytest.mark.verifies("galerkin-idempotency-4pi")
def test_pi_R_is_4pi_identity_on_band_limited(lebedev_L_pair):
    """Pi · R = 4pi · I on the band-limited coefficient space.

    Sister of tests/numerics/test_projection_operators.py:194-231
    (TestGalerkinIdempotencyOnLebedev) — this one constructs the
    operators through the new SphericalHarmonicSpace API; the legacy
    test pins the same identity through the legacy (Π, R) pair.
    Both MUST agree to quadrature precision.
    """
    measure, L = lebedev_L_pair
    space = SphericalHarmonicSpace.from_L(L)
    basis = SphericalHarmonicBasis(L=L)
    Y = basis.evaluate(measure.nodes)
    M = MomentProjection.from_spherical_harmonic_space(
        space, weights=measure.weights, Y=Y,
    )
    R = ReconstructionOperator.from_spherical_harmonic_space(space, Y)

    rng = np.random.default_rng(seed=42)
    c = rng.standard_normal((L + 1, 2 * L + 1))
    for ell in range(L + 1):
        c[ell, 2 * ell + 1 :] = 0.0

    out = M.apply(R.apply(c))
    np.testing.assert_allclose(out, 4.0 * np.pi * c,
                               rtol=1e-10, atol=1e-12)


@pytest.mark.l1
@pytest.mark.catches("ERR-039")
@pytest.mark.verifies("hilbert-adjoint-equals-metric-times-S0")
def test_H_equals_g_C_times_S0(lebedev_L_pair):
    r"""M.H computed generically equals g_C ⊙ S_0(c).

    The .H machinery composes g_V^{-1} ∘ M.T ∘ g_W. With g_V the
    metric on the coefficient space (g_C, the MomentMassMatrix) and
    g_W the angular metric (quadrature weights), the result on
    coefficient input c is (4pi/(2l+1)) · sum_n Y_lm(n) c_lm — the
    "metric-aware" synthesis.
    """
    measure, L = lebedev_L_pair
    space = SphericalHarmonicSpace.from_L(L)
    basis = SphericalHarmonicBasis(L=L)
    Y = basis.evaluate(measure.nodes)
    M = MomentProjection.from_spherical_harmonic_space(
        space, weights=measure.weights, Y=Y,
    )

    rng = np.random.default_rng(seed=99)
    c = rng.standard_normal((L + 1, 2 * L + 1))
    for ell in range(L + 1):
        c[ell, 2 * ell + 1 :] = 0.0

    actual = M.H.apply(c)  # the generic adjoint composition
    # Expected: g_C · S_0(c) with g_C[l] = 4pi/(2l+1)
    g_C_per_ell = 4.0 * np.pi / (2.0 * np.arange(L + 1) + 1.0)
    c_scaled = c.copy()
    for ell in range(L + 1):
        c_scaled[ell, :] *= g_C_per_ell[ell]
    expected = np.einsum("nlm,lm->n", Y, c_scaled)
    np.testing.assert_allclose(actual, expected, rtol=1e-12)
```

### B.2 Additional new tests required for P1.1–P1.4 correctness

The five identities above pin the ERR-039 math. The following four
tests pin the constructor / API surface that P1.1–P1.4 introduces.

```python
# ─────────────────────────────────────────────────────────────────────
# B.2 — constructor / API surface tests
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
@pytest.mark.verifies("sh-basis-mass-matrix-quadrature-independence")
@pytest.mark.parametrize("quadrature_name", ["lebedev", "level_symmetric", "product_gauss"])
def test_mass_matrix_under_multiple_quadratures(quadrature_name):
    """SphericalHarmonicBasis.mass_matrix is exact-by-construction
    against any quadrature that integrates Y_l Y_l' for l + l' <= 2L.

    Probe: for L=2, deg ≥ 4 quadrature gives g_C = diag(4pi/(2l+1));
    deg < 4 gives a discretisation error proportional to truncation.
    Test pins both regimes:
      - exact regime: |M_mass - g_C| < 1e-12
      - undersampled regime: error grows as expected
    """
    from orpheus.numerics.quadrature import (
        lebedev_sphere, level_symmetric_sn, product_gauss_sphere,
    )
    L = 2
    basis = SphericalHarmonicBasis(L=L)

    # Exact-degree rule per family
    if quadrature_name == "lebedev":
        measure = lebedev_sphere(7)  # deg=7 >= 2L=4
    elif quadrature_name == "level_symmetric":
        measure = level_symmetric_sn(8)  # LS_8 has deg ≥ 5
    else:  # product_gauss
        measure = product_gauss_sphere(n_polar=8, n_azimuthal=16)

    M_mass = basis.mass_matrix(measure)
    expected_per_ell = 4.0 * np.pi / (2.0 * np.arange(L + 1) + 1.0)
    diag = M_mass.diagonal()
    for ell in range(L + 1):
        np.testing.assert_allclose(
            diag[ell, : 2 * ell + 1], expected_per_ell[ell],
            rtol=1e-10,
            err_msg=f"quad={quadrature_name}, ell={ell}",
        )


@pytest.mark.l0
@pytest.mark.verifies("moment-projection-codomain-type")
def test_moment_projection_codomain_is_spherical_harmonic_space():
    """MomentProjection.codomain returns a SphericalHarmonicSpace of
    matching L. Type-level guarantee that .H propagation works.
    """
    from orpheus.numerics.spaces import SphericalHarmonicSpace
    measure = lebedev_sphere(7)
    L = 2
    basis = SphericalHarmonicBasis(L=L)
    space = SphericalHarmonicSpace.from_L(L)
    M = MomentProjection.from_spherical_harmonic_space(
        space, weights=measure.weights, Y=basis.evaluate(measure.nodes),
    )
    cod = M.codomain
    assert isinstance(cod, SphericalHarmonicSpace)
    assert cod.L == L
    assert cod.shape == (L + 1, 2 * L + 1)
    # Equality is (name, shape) — same L gives same codomain
    assert cod == SphericalHarmonicSpace.from_L(L)


@pytest.mark.l0
@pytest.mark.verifies("from-spherical-harmonic-space-roundtrip")
def test_from_spherical_harmonic_space_roundtrip():
    """Construction via the new classmethod gives a (M, R) pair that
    is bit-identical to construction via the legacy (weights, Y)
    constructors — the migration is API-only, no numerical change.
    """
    measure = lebedev_sphere(13)
    L = 3
    basis = SphericalHarmonicBasis(L=L)
    space = SphericalHarmonicSpace.from_L(L)
    Y = basis.evaluate(measure.nodes)

    # New API
    M_new = MomentProjection.from_spherical_harmonic_space(
        space, weights=measure.weights, Y=Y,
    )
    R_new = ReconstructionOperator.from_spherical_harmonic_space(space, Y)

    # Legacy API (preserved as shim per plan §P1.3)
    from orpheus.numerics.projection import (
        HarmonicMomentProjection, HarmonicMomentReconstruction,
    )
    M_legacy = HarmonicMomentProjection(
        weights=measure.weights, Y=Y, L=L,
    )
    R_legacy = HarmonicMomentReconstruction.from_Y(Y)

    rng = np.random.default_rng(seed=12345)
    psi = rng.standard_normal(measure.n_points)
    c = rng.standard_normal((L + 1, 2 * L + 1))
    for ell in range(L + 1):
        c[ell, 2 * ell + 1 :] = 0.0

    np.testing.assert_array_equal(M_new.apply(psi), M_legacy.apply(psi))
    np.testing.assert_array_equal(R_new.apply(c), R_legacy.apply(c))


@pytest.mark.l1
@pytest.mark.catches("ERR-039")
@pytest.mark.verifies("T-vs-H-representation-vs-hilbert")
def test_T_is_w_weighted_representation_transpose(lebedev_L_pair):
    """The two adjoint identities both hold under the Euclidean and
    W-weighted inner products on V respectively:

        Euclidean on V:    <Pi psi, c>_C = <psi, Pi.T c>_V_Euclidean
        W-weighted on V:   <Pi psi, c>_C = <psi, Pi.H c>_V_W

    .T is the representation transpose = w_n · S_0.
    .H is the Hilbert adjoint = g_C · S_0.

    Cf. test_projection_operators.py:248-271 (legacy form); this is
    the corrected form under P1.4. The legacy test rewrites in §A.1.
    """
    measure, L = lebedev_L_pair
    space = SphericalHarmonicSpace.from_L(L)
    basis = SphericalHarmonicBasis(L=L)
    Y = basis.evaluate(measure.nodes)
    M = MomentProjection.from_spherical_harmonic_space(
        space, weights=measure.weights, Y=Y,
    )

    rng = np.random.default_rng(seed=7)
    psi = rng.standard_normal(measure.n_points)
    c = rng.standard_normal((L + 1, 2 * L + 1))
    for ell in range(L + 1):
        c[ell, 2 * ell + 1 :] = 0.0

    # <Pi psi, c>_C — Euclidean on coefficient space
    M_psi = M.apply(psi)
    lhs = float(np.sum(M_psi * c))

    # <psi, Pi.T c>_V_Euclidean — Euclidean on V; Pi.T already
    # carries w_n
    T_c = M.T.apply(c)
    rhs_euclidean = float(np.sum(psi * T_c))
    np.testing.assert_allclose(lhs, rhs_euclidean, rtol=1e-12, atol=1e-14)

    # <Pi psi, c>_C with C in its native metric is non-trivial; the
    # canonical identity is the Euclidean-on-C form above. The
    # W-weighted-on-V form uses Pi.H:
    #   <Pi psi, c>_C_native = <psi, Pi.H c>_V_W_native
    # where C_native carries g_C. Verify the metric-aware form:
    g_C_per_ell = 4.0 * np.pi / (2.0 * np.arange(L + 1) + 1.0)
    c_in_C_metric = c.copy()
    for ell in range(L + 1):
        c_in_C_metric[ell, :] *= g_C_per_ell[ell]
    lhs_metric = float(np.sum(M_psi * c_in_C_metric))

    H_c = M.H.apply(c)
    # <psi, H_c>_V_W = sum_n w_n psi_n (H c)_n
    rhs_metric = float(np.sum(measure.weights * psi * H_c))
    np.testing.assert_allclose(lhs_metric, rhs_metric, rtol=1e-12,
                               atol=1e-14)
```

### B.3 Forward-compatibility test that paves the way for Phase 2

```python
@pytest.mark.l1
@pytest.mark.verifies("space-equality-by-name-shape")
def test_spherical_harmonic_space_equality_by_name_shape():
    """Two SphericalHarmonicSpace instances with equal L compare equal.

    Per plan §P1.2: shape == (L+1, 2L+1); __eq__/__hash__ inherit
    (name, shape) from FunctionSpace. Forward-compatible with Phase
    2's DualSpace/TensorProductSpace tests where Space identity is
    load-bearing.
    """
    from orpheus.numerics.spaces import SphericalHarmonicSpace
    a = SphericalHarmonicSpace.from_L(3)
    b = SphericalHarmonicSpace.from_L(3)
    c = SphericalHarmonicSpace.from_L(4)
    assert a == b
    assert hash(a) == hash(b)
    assert a != c
```

---

## Section C — L1 MMS gates that must continue passing

The single production consumer of the `(M, Λ, R)` pipeline is
`ScatteringOperator.build_aniso_source` (`sn/scattering.py:609-657`).
Phase 1 changes neither the algorithm nor the einsum trees in this
path; the bit-identity contract holds.

| Test file | What it exercises | Disposition |
|---|---|---|
| `tests/sn/test_mms_aniso.py::test_sn_p1_aniso_mms_converges_second_order` | P1 anisotropic MMS slab; `solve_sn_fixed_source` with `scattering_order=1`. Calls `_build_aniso_scattering` → `ScatteringOperator.build_aniso_source` → the (M, Λ, R) pipeline. | **STRICT bit-identical** through P1.1–P1.6. After P1.7 (the `_build_rhs_cartesian` migration), MUST still show O(h²) convergence; the per-iterate flux may drift within `outer_iters × ULP` but the convergence rate is invariant. Decorator stack unchanged: `@pytest.mark.l1` + `@pytest.mark.verifies("transport-cartesian", "dd-cartesian-1d", "dd-slab", "pn-scatter", "sn-mms-p1-psi", "sn-mms-p1-qext")`. |
| `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py` | Curvilinear MMS aniso convergence (sphere / cylinder; uses the same R-Λ-M pipeline). | **STRICT bit-identical** through P1.1–P1.6. Post-P1.7: convergence rate invariant; per-iterate drift within iteration floor. |
| `tests/sn/test_scattering_operator.py` | Direct unit tests of `ScatteringOperator.build_aniso_source` and the (M, Λ, R) pipeline composition. | **STRICT bit-identical** through P1.1–P1.7 (the pipeline composition is unchanged; only `_build_rhs_cartesian`'s parallel hand-rolled form retires). |
| `tests/sn/test_legendre_moment_scattering.py` | Direct unit tests of `LegendreMomentScattering` (Λ). | **STRICT bit-identical** (Λ is untouched). |
| `tests/sn/test_harmonic_moment_field.py` | Typed `HarmonicMomentField` tests. | **STRICT bit-identical** (P3.3 not P1 relocates this class). |
| `tests/sn/test_phase_c_mms.py` | Phase-C MMS tests (curvilinear M-M sweep). | **STRICT bit-identical**. |

**Invariance claim**: P1.1–P1.6 touch ONLY the rename, the new
basis/space classes, the `.T`/`.H` semantics on `MomentProjection`,
and the deletion of `assert_galerkin_idempotency`. None of these
intersect with the production `build_aniso_source` einsum chain. P1.7
touches `_build_rhs_cartesian` only — the parallel hand-rolled form
predating the (M, Λ, R) pipeline. The pipeline itself is invariant
through all of Phase 1.

**Failure-mode check (per `vv-principles` §"6 AI failure modes"):**
- Sign flip in `(2ℓ+1)` migration (mode 1) — caught by the L1 MMS
  convergence gate AND the snapshot suite (`sphere_2g_p1_aniso_dd_n20.npz`).
- Convention drift between `MomentMassMatrix` and the inline literal
  (mode 6) — caught by §B.1 `test_R_equals_2l_plus_1_times_S0` and
  the snapshot suite.
- Missing `w_n` in `.T` (mode 3) — caught by §B.2
  `test_T_is_w_weighted_representation_transpose` and the rewritten
  §A.1 `TestApplyTransposeIsWWeightedAdjoint`.
- Variable swap `M.T` ↔ `M.H` (mode 2) — caught by §B.1
  `test_H_equals_g_C_times_S0` AND the type-level `codomain` check
  (§B.2 `test_moment_projection_codomain_is_spherical_harmonic_space`).

---

## Section D — Phase 2 invariants spec (preview, naming pattern)

Phase 2 is out of scope for this session but the P1.5 test file
**establishes the naming pattern** that Phase 2 inherits. Names and
decorator stacks for the Phase 2 tests:

| Test | Where | Decorator stack |
|---|---|---|
| `test_dual_of_dual_is_identity` | `tests/numerics/test_dual_space.py` | `@pytest.mark.l1` + `@pytest.mark.verifies("dual-dual-identity")` |
| `test_composite_adjoint_distributes` | `tests/numerics/test_operator_adjoint_algebra.py` | `@pytest.mark.l1` + `@pytest.mark.verifies("composite-adjoint-distributivity")` (asserts `(A @ B).H == B.H @ A.H` on non-self-adjoint pairs via `assert_adjoint_consistency`) |
| `test_assert_adjoint_consistency_W_inner_product` | `tests/numerics/test_operator_adjoint_algebra.py` | `@pytest.mark.l1` + `@pytest.mark.verifies("hilbert-adjoint-pairing-identity")` (the canonical ⟨A x, y⟩_W = ⟨x, A.H y⟩_V over random samples) |
| `test_tensor_product_shape` | `tests/numerics/test_tensor_product_space.py` | `@pytest.mark.l0` + `@pytest.mark.verifies("tensor-product-shape-concatenation")` (`(V * W).shape == V.shape + W.shape`) |
| `test_direct_sum_shape_along_coupled_axis` | `tests/numerics/test_direct_sum_space.py` | `@pytest.mark.l0` + `@pytest.mark.verifies("direct-sum-coupled-axis-sum")` |
| `test_sum_of_tensor_products_streaming_form` | `tests/numerics/test_sum_of_tensor_products_operator.py` | `@pytest.mark.l1` + `@pytest.mark.verifies("streaming-as-sum-of-tensor-products")` — **issue #172 anchor** |
| `test_tensor_product_operator_adjoint_distributivity` | `tests/numerics/test_tensor_product_operator.py` | `@pytest.mark.l1` + `@pytest.mark.verifies("tensor-operator-adjoint-distributivity")` — **issue #173 anchor** |

**Naming convention** (inherited by Phase 2):
- File name: `test_<concept>.py` mirrors the production module name
  (`spaces/spherical_harmonic_space.py` → `test_spherical_harmonic_space.py`;
  `spaces/dual_space.py` → `test_dual_space.py`).
- Test function: `test_<claim>` where `<claim>` is the math identity
  in lowercase-underscore form.
- Decorators: ordered as `@pytest.mark.l{0,1}`,
  `@pytest.mark.catches("ERR-NNN")` (if applicable),
  `@pytest.mark.verifies("<sphinx-label>")`.

**Forward-compatibility action in P1.5**: when registering the
Sphinx labels under `docs/theory/spherical_harmonics.rst`, use a
single prefix `sh-` (e.g. `sh-space-metric`, `sh-quadrature-gram`)
so Phase 2's labels (`dual-`, `tensor-`, `sum-`) sit alongside in
the same namespace.

---

## Section E — P3.1 import-linter test spec

Test file: `tests/test_layer_imports.py`. The test enforces the
import contract defined in plan §P3.0 (layer table). It is written
in P3.1 (FIRST step of Phase 3) before any module moves; its initial
failures are the migration to-do list.

### Structure

```python
"""Import-linter: enforce the L1 / L2 / L3 layer contract.

Per plan §P3.0:
  L1 (numerics/)       — math primitives; knows no transport
  L2 (transport/)      — transport vocabulary; method-agnostic
  L3 (sn/, pn/, moc/,  — one method's machinery; method-specific
      cp/, mc/, diffusion/, kinetics/)

Imports flow only L3 -> L2 -> L1 (and L3 -> L3 within the same
package). Forbidden edges raise the test.
"""
from __future__ import annotations

import ast
import pathlib
from collections.abc import Iterable

import pytest

ORPHEUS_ROOT = pathlib.Path(__file__).parent.parent / "orpheus"

# Layer assignment by top-level package directory
L1_PACKAGES: frozenset[str] = frozenset({"numerics"})
L2_PACKAGES: frozenset[str] = frozenset({"transport"})
L3_PACKAGES: frozenset[str] = frozenset({
    "sn", "pn", "moc", "cp", "mc", "diffusion", "kinetics",
})
INPUT_PACKAGES: frozenset[str] = frozenset({"geometry", "data"})

# Forbidden edges (source_package → forbidden_target_packages)
FORBIDDEN_EDGES: dict[str, frozenset[str]] = {
    # L1 cannot import L2 or L3
    "numerics": L2_PACKAGES | L3_PACKAGES,
    # L2 cannot import L3
    "transport": L3_PACKAGES,
    # L3 packages cannot import each other (each method is sealed)
    "sn":         L3_PACKAGES - {"sn"},
    "pn":         L3_PACKAGES - {"pn"},
    "moc":        L3_PACKAGES - {"moc"},
    "cp":         L3_PACKAGES - {"cp"},
    "mc":         L3_PACKAGES - {"mc"},
    "diffusion":  L3_PACKAGES - {"diffusion"},
    "kinetics":   L3_PACKAGES - {"kinetics"},
}

# Whitelisted exception edges (plan §P3.1):
#   - numerics/eigenvalue.py legacy shim, allowed until P3.4 retires it
WHITELIST: frozenset[tuple[str, str]] = frozenset({
    # (source_module_relative_path, target_top_level_package)
    # e.g. ("numerics/eigenvalue.py", "sn") — until retired in P3.4
})


def _iter_python_modules(root: pathlib.Path) -> Iterable[pathlib.Path]:
    """Yield every .py file under orpheus/, excluding __pycache__."""
    for p in root.rglob("*.py"):
        if "__pycache__" in p.parts:
            continue
        yield p


def _top_level_package(rel_path: pathlib.Path) -> str:
    """Return the L1/L2/L3 package name for a file under orpheus/."""
    return rel_path.parts[0]


def _imports_of(module_path: pathlib.Path) -> list[tuple[str, bool]]:
    """Parse a module's imports.

    Returns list of (full_module_name, is_type_checking) tuples.
    `from orpheus.X.Y import Z` → ("orpheus.X.Y", False).
    `import orpheus.X.Y as Z`  → ("orpheus.X.Y", False).
    Imports inside `if TYPE_CHECKING:` blocks are marked is_type_checking=True.
    """
    src = module_path.read_text()
    tree = ast.parse(src, filename=str(module_path))
    results: list[tuple[str, bool]] = []

    def _visit(node: ast.AST, in_type_checking: bool) -> None:
        if isinstance(node, ast.If):
            # Detect `if TYPE_CHECKING:` block
            is_tc_guard = (
                isinstance(node.test, ast.Name)
                and node.test.id == "TYPE_CHECKING"
            )
            for child in node.body:
                _visit(child, in_type_checking or is_tc_guard)
            for child in node.orelse:
                _visit(child, in_type_checking)
            return
        if isinstance(node, ast.ImportFrom) and node.module:
            results.append((node.module, in_type_checking))
        elif isinstance(node, ast.Import):
            for alias in node.names:
                results.append((alias.name, in_type_checking))
        for child in ast.iter_child_nodes(node):
            _visit(child, in_type_checking)

    _visit(tree, in_type_checking=False)
    return results


def _check_module(module_path: pathlib.Path) -> list[str]:
    """Return list of violation messages for a single module."""
    rel = module_path.relative_to(ORPHEUS_ROOT)
    src_pkg = _top_level_package(rel)
    if src_pkg not in FORBIDDEN_EDGES:
        return []  # not a tracked source layer (e.g. tests/_harness)

    forbidden = FORBIDDEN_EDGES[src_pkg]
    violations: list[str] = []
    for module_name, is_tc in _imports_of(module_path):
        if not module_name.startswith("orpheus."):
            continue
        tgt_pkg = module_name.split(".")[1]
        if tgt_pkg not in forbidden:
            continue
        # TYPE_CHECKING imports of L3 into L1/L2 are allowed
        if is_tc and src_pkg in (L1_PACKAGES | L2_PACKAGES):
            continue
        # Whitelist
        rel_str = str(rel).replace("\\", "/")
        if (rel_str, tgt_pkg) in WHITELIST:
            continue
        violations.append(
            f"{rel} imports {module_name} "
            f"(forbidden: {src_pkg} → {tgt_pkg})"
        )
    return violations


# Parametrise per module so a failure isolates the offender
_ALL_MODULES = sorted(_iter_python_modules(ORPHEUS_ROOT))


@pytest.mark.foundation
@pytest.mark.parametrize("module_path", _ALL_MODULES, ids=str)
def test_no_forbidden_imports(module_path: pathlib.Path) -> None:
    violations = _check_module(module_path)
    assert not violations, "\n".join(violations)
```

### Behaviour through Phase 3

- **P3.1 (this test lands)**: expect failures for every misplaced
  module — those failures ARE the migration to-do list. The first
  P3.2 / P3.3 commits make them go green one module at a time.
- **`numerics/eigenvalue.py`**: add to `WHITELIST` until P3.4
  retires it.
- **`AngularFlux.from_flat_with_traces`** (the L2-leaks-to-L3 case
  per plan CC.4): the P3.3 design splits the class. The L2 base is
  clean; the SN-specific `from_flat_with_traces` lives in
  `sn/angular_flux_b1pp_adapter.py` (L3) — no whitelist needed once
  the split lands.
- **`TYPE_CHECKING` guards**: the test's `_visit` function detects
  `if TYPE_CHECKING:` blocks and skips the forbidden-edge check for
  L3-type-only imports in L1/L2. Sample exemption: `transport/`
  modules may type-annotate with `from orpheus.sn import SNMesh`
  inside a TYPE_CHECKING block (string-only annotation).

### Decorator and tagging

- `@pytest.mark.foundation` (per `vv-principles` §"V&V level
  taxonomy") — this is a software invariant, NOT a verifies-an-
  equation claim. The harness recognises the tag.
- NO `@pytest.mark.verifies(...)` — there is no Sphinx equation
  label.
- Parametrised per-module so a failing test's name names the
  offender.

---

## Gates checklist — ordered invocations the implementer runs

Run after EACH P1.x step. Every gate MUST be green before advancing.

### After P1.0 (recon — already complete)
- `pytest -q tests/numerics/test_projection_operators.py tests/numerics/test_spherical_harmonics.py`
- Establishes the green baseline.

### After P1.1 (`SphericalHarmonicBasis`)
- `pytest -q tests/numerics/test_spherical_harmonics.py` — STRICT bit-identical (shim works).
- `pytest -q tests/numerics/test_projection_operators.py` — STRICT bit-identical (no projection change yet).
- `pytest -q tests/sn/test_scattering_operator.py` — STRICT bit-identical (production unchanged).
- `pytest -q tests/sn/regression/test_dd_regression.py -m regression` — STRICT bit-identical (existing tolerances).

### After P1.2 (`SphericalHarmonicSpace`, `MomentMassMatrix`)
- All P1.1 gates green.
- New: `pytest -q tests/numerics/test_spherical_harmonic_space.py::test_space_inner_product_weights_equal_4pi_over_2l_plus_1`
- New: `pytest -q tests/numerics/test_spherical_harmonic_space.py::test_basis_mass_matrix_against_lebedev`
- New: `pytest -q tests/numerics/test_spherical_harmonic_space.py::test_spherical_harmonic_space_equality_by_name_shape`

### After P1.3 (`MomentProjection` / `ReconstructionOperator` rewire)
- All P1.2 gates green.
- `pytest -q tests/numerics/test_projection_operators.py` — STRICT bit-identical EXCEPT class-name imports adapt (see §A.1; `TestABCs`, `TestHarmonicMomentProjection*` rename to `TestMomentProjection*`).
- New: `pytest -q tests/numerics/test_spherical_harmonic_space.py::test_moment_projection_codomain_is_spherical_harmonic_space`
- New: `pytest -q tests/numerics/test_spherical_harmonic_space.py::test_from_spherical_harmonic_space_roundtrip`

### After P1.4 (`.T` / `.H` correction)
- All P1.3 gates green.
- **Rewritten under P1.4**: `pytest -q tests/numerics/test_projection_operators.py::TestApplyTransposeIsWWeightedAdjoint` (now `TestTReturnsWWeightedTranspose`).
- **Rewritten under P1.4**: `pytest -q tests/numerics/test_projection_operators.py::TestGalerkinAdjointPairing` (RHS drops external `w_n`).
- New: `pytest -q tests/numerics/test_spherical_harmonic_space.py::test_T_is_w_weighted_representation_transpose`
- New: `pytest -q tests/numerics/test_spherical_harmonic_space.py::test_R_equals_2l_plus_1_times_S0`
- New: `pytest -q tests/numerics/test_spherical_harmonic_space.py::test_pi_R_is_4pi_identity_on_band_limited`
- New: `pytest -q tests/numerics/test_spherical_harmonic_space.py::test_H_equals_g_C_times_S0`
- `pytest -q tests/sn/test_scattering_operator.py` — STRICT bit-identical (production calls `R.apply`, untouched).
- `pytest -q tests/sn/regression/test_dd_regression.py -m regression` — STRICT bit-identical.

### After P1.5 (the new test file as a whole)
- All P1.4 gates green; this step is the canonical "did P1.5 land?" gate.
- `pytest -q tests/numerics/test_spherical_harmonic_space.py` — entire file green (5 plan-identity + 4 API-surface + 1 forward-compat = 10 tests).
- `pytest -q -k catches_ERR_039` (or `-m "catches('ERR-039')"`) — confirms the ERR-039 catalog test set has been expanded.

### After P1.6 (docstring surgery + `assert_galerkin_idempotency` retired)
- All P1.5 gates green.
- **Deleted**: `pytest -q tests/numerics/test_projection_operators.py::TestAssertGalerkinIdempotencyMethod` — MUST fail to collect (the class no longer exists).
- `pytest -q tests/numerics/test_projection_operators.py` — strict green minus the deleted class.
- `python -m tests._harness.audit` — verify the V&V matrix carries the new `sh-*` equation labels and no orphan `catches("ERR-039")` entries remain unmapped.
- `sphinx-build -W docs docs/_build/html` — clean (the new `sh-space-metric` etc. labels were added; the warning boxes at `projection.py:402-422` and `sn/scattering.py:555-560` collapsed to one line each).

### After P1.7 (retire `sn/solver.py:930` inline `(2ℓ+1)`)
- All P1.6 gates green.
- `pytest -q tests/sn/test_scattering_operator.py` — STRICT bit-identical (the `R · Λ · M` pipeline at `sn/scattering.py:609-657` is the canonical site; unchanged).
- `pytest -q tests/sn/test_mms_aniso.py` — convergence rate O(h²) preserved.
- `pytest -q tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py` — convergence rate preserved.
- `pytest -q tests/sn/regression/test_dd_regression.py -m regression` — slab `*_dd_n*` snapshots STRICT bit-identical at `rtol=1e-12`; curvilinear `sphere_*`, `cyl_*` snapshots within existing `rtol=5e-6` floor. If any snapshot breaches its existing bound, the migration is wrong (see §A.3 "P1.7 bit-identity criteria"). DO NOT loosen the tolerance.
- For each `*_homogeneous_*` snapshot, verify `result.keff` matches `k_inf = νΣ_f/Σ_a` to within `outer_iters × ULP` — confirms criterion 2 (structurally-independent analytical reference).
- Document the snapshot verification in the P1.7 commit message per the §A.3 policy.
- **Final gate before Phase 1 close**: full suite — `pytest -q` — green.

### Phase-1-close audit
- `python -m tests._harness.audit` — V&V matrix fully populated for ERR-039; no orphan equation labels under `docs/theory/spherical_harmonics.rst`.
- `grep -rn "(2 \* l + 1)\|(2\*l+1)\|two_l_plus_one" orpheus/` — exactly TWO sources remain: `SphericalHarmonicSpace.addition_theorem_factor` (the canonical) and `ReconstructionOperator.two_l_plus_one` field initialised from it (the cache). The `_build_rhs_cartesian` site at `sn/solver.py:862, 884-898, 930` is gone. The test fixture at `tests/numerics/test_projection_operators.py:173, 293` is the test-as-spec form per plan §P1.3 and stays.
- `grep -rn "assert_galerkin_idempotency" orpheus/ tests/` — zero hits.
- `sphinx-build -W docs docs/_build/html` — clean.

---

## Self-improvement entries

### New failure mode → `vv-principles` skill update (NONE for this plan)

This plan covers existing failure modes (modes 1, 2, 3, 6 in
`vv-principles` §"6 AI failure modes"). No new failure mode is
introduced. No skill table update required.

### Plan-rejection counter-examples (none yet)

If qa or the user rejects any row in this plan, log a one-paragraph
counter-example under
`.claude/agent-memory/test-architect/feedback_*.md` per the
test-architect AGENT.md self-improvement directive.

---

## Pointers

- Plan: `/Users/rodrigo/git/nuclear/ORPHEUS/.claude/worktrees/moment-space-and-layering/.claude/plans/moment_space_and_layering_plan.md`
- Legacy test (the canonical pin): `/Users/rodrigo/git/nuclear/ORPHEUS/.claude/worktrees/moment-space-and-layering/tests/numerics/test_projection_operators.py`
- SN regression snapshot driver: `/Users/rodrigo/git/nuclear/ORPHEUS/.claude/worktrees/moment-space-and-layering/tests/sn/regression/test_dd_regression.py`
- The production consumer: `/Users/rodrigo/git/nuclear/ORPHEUS/.claude/worktrees/moment-space-and-layering/orpheus/sn/scattering.py` (lines 525-657)
- The inline-`(2ℓ+1)` retirement site: `/Users/rodrigo/git/nuclear/ORPHEUS/.claude/worktrees/moment-space-and-layering/orpheus/sn/solver.py` (lines 862, 884-898, 930)
- Skills: `vv-principles`, `coding-elegance` Pattern 7 (definition-site normalisation), `numerical-bug-signatures` Signature 3 (transpose convention drift)
- ERR-039: `/Users/rodrigo/git/nuclear/ORPHEUS/.claude/skills/vv-principles/error_catalog.md`
