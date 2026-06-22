---
name: issue-196-phase-g-step1-verification-gates
description: Verification-gate design for Phase G Step 1 (Issue #196) — type-system promotion of DiamondDifference → SNCellOperator(LinearOperator) and MorelMontryAngularSweep → AngularRedistribution(LinearOperator). Seven gates, with explicit Phase F twin-path defense; sketches three pytest skeleton modules; no production-code change in scope.
metadata:
  type: project
  branch: refactor/sn-operator-algebra
  step: Phase G Step 1
  authored: 2026-05-12
---

# Phase G Step 1 — Verification-gate design

**Purpose.** Phase G Step 1 promotes two existing classes to
`LinearOperator` subclasses with NO mathematical change:

- `DiamondDifference` → `SNCellOperator(LinearOperator)`. Adds
  `apply(cell_avg) → residual` (cell-level discretised operator
  residual `L_cell · ψ_avg − q`) on top of the existing
  `solve(visit, total_xs, source, upstream_state) → CellResult` which
  IS the current `DiamondDifference.update` body. Declares
  `capabilities = frozenset({CAP_APPLY, CAP_SOLVE})`.
- `MorelMontryAngularSweep` → `AngularRedistribution(LinearOperator)`.
  Wraps the M-M angular recurrence `__call__` and composes the
  Carlson seed (`carlson_inward_sweep_from_source` / `CarlsonInwardSweep`)
  as a private inner `LinearOperator`. Declares
  `capabilities = frozenset({CAP_APPLY, CAP_SOLVE})` (the latter
  exposes inverse of the M-M recurrence — algebraically the
  recurrence is upper-triangular in ordinate index so it inverts by
  back-substitution).

This is a **pure type-system promotion**. The plan goal is
bit-identity (`np.array_equal`) on all paths the existing tests
exercise.

**Failure-mode profile this defends against.** Phase F bug: Phase D
fixed the apply-path Carlson seed; the SI/sweep path kept the broken
`np.zeros` seed for a month; the regression suite couldn't see it
because Gate 1.1 only ran the apply path AND the curvilinear
snapshots were SI-generated under the wrong seed (so bit-identity
preserved the bug). Step 1's gates MUST close this class of bug by
construction.

---

## CRITICAL: Gate matrix overview

| # | Gate | Level | Catches |
|---|---|---|---|
| 1 | SNCellOperator.solve ≡ DiamondDifference.update (bit-identity) | foundation | Mode 6 (convention drift across promotion) |
| 2 | apply ∘ solve = q (Protocol round-trip) | foundation | Modes 1–4 in the cell algebra |
| 3 | Capability surface (CAP_APPLY, CAP_SOLVE, MissingCapability) | foundation | Pattern 4 (illegal-state-unrepresentable) |
| 4 | Curvilinear + Cartesian coverage of gates 1–3 | foundation | Mode 7 (test-design simplification — slab nulls angular redistribution) |
| 5 | **Phase F twin-path defense** — apply-path + sweep-path SAME SNCellOperator on sphere_2g_3reg | l0 | ERR-026 manifestation #6 + #7 by construction |
| 6 | AngularRedistribution algebraic identities (linearity, flat-flux closure, Carlson seed equivalence) | foundation | Modes 2, 3 in M-M recurrence |
| 7 | Regression-snapshot invariance (all 11 snapshots bit-identical at rtol=1e-12) | regression | Mode 5 (index/wiring drift in promotion) |

---

## Gate 1 — `SNCellOperator.solve ≡ DiamondDifference.update` (bit-identity)

**Claim layer.** Software-invariant claim — this is a foundation-tag
gate, NOT an L0/L1/L2 claim about physics. The promotion is pure
plumbing; the gate pins that plumbing didn't lose any algebra.

**Pillar.** Closed-form. The reference is `DiamondDifference.update`
itself — a sibling code path computing the same scalar formula. This
is **NOT structural independence** (both paths are in the same
codebase); the gate is a regression contract against drift during
the promotion, not a verification claim about correctness.

**Why this gate exists.** The Phase F bug shipped under apparently-
green tests because the bit-identity invariant ran only the
apply-path. Step 1 explicitly parametrizes over BOTH consumers — the
apply-path consumer (`transport_operator_matvec_spherical`) AND the
SI-path consumer (`_sweep_1d_spherical`) — with the SAME parameter
sweep on each, so a divergence between the two consumers under the
promotion cannot hide.

**Parametrize matrix:**

| geometry | n_groups | regions | quadrature | direction (sign of μ) | cell_idx |
|---|---|---|---|---|---|
| slab | {1, 2} | {homog, 3-region} | GL-4 | {+1, −1} | {0, interior, last} |
| sphere | {1, 2} | {homog, 3-region} | GL-4, GL-8 | {+1, −1} | {0 (pole), interior, last (outer face)} |
| cylinder | {1, 2} | {homog, 3-region} | LS-4, product | {+1, −1} | {0 (axis), interior, last} |

For each combination, build the synthetic `(visit, total_xs, source,
upstream_state)` packet and assert:

```python
assert np.array_equal(
    SNCellOperator().solve(visit, total_xs, source, upstream).cell_average_flux,
    DiamondDifference().update(visit, total_xs, source, upstream).cell_average_flux,
)
```

Plus the same on `outgoing_spatial_flux` and `outgoing_angular_state`.

**Coverage.** This subsumes the existing `tests/sn/spatial/test_diamond.py`
gates (28 tests) but parametrizes the SAME inputs through the new
operator API in addition. Net: 28 existing + ~36 new = ~64.

---

## Gate 2 — Apply-vs-solve consistency (Protocol round-trip)

**Claim layer.** Software-invariant — this is the contract claim
`apply ∘ solve = id` that the LinearOperator capability set
advertises.

**Pillar.** Closed-form algebraic identity. For a `LinearOperator`
`L` with both `CAP_APPLY` and `CAP_SOLVE`:

```
L.apply(L.solve(q)) ≡ q   (up to FP rounding, rtol=1e-12)
```

**Why this gate exists.** Without round-trip, the apply and solve
methods could disagree on what the operator means — the same defect
shape as the Phase F twin-path bug, but at the type-system level
rather than at the call-site level.

**Critical design choice — source MUST be non-trivial.** A round-
trip on `q = 0` passes trivially because both `apply` and `solve`
preserve zero. The gate MUST cover non-trivial sources to exercise
the algebraic structure. Specifically:

- Non-uniform `source` (per-group, per-cell variation).
- Non-trivial `upstream_state` (positive `spatial_upstream`,
  positive `angular_upstream` for curvilinear).
- Heterogeneous `total_xs` (per-group, multi-region).

**Parametrize matrix:**

| geometry | n_groups | source | upstream_state | total_xs |
|---|---|---|---|---|
| slab | {1, 2} | {Q ∝ Σ_t (flat-ψ source), Q random, Q linear in r} | {ψ_in = 0, ψ_in = 0.5} | {Σ_t homog, Σ_t 3-region} |
| sphere | {1, 2} | same three | {(ψ_s, ψ_a) = (0,0), (0.5, 0.3)} | same two |
| cylinder | {1, 2} | same three | same two | same two |

**Assertion shape:**

```python
result = SNCellOperator().solve(visit, total_xs, source, upstream)
# apply: compute the cell-level residual at the just-solved cell_avg
residual = SNCellOperator().apply(
    cell_avg=result.cell_average_flux,
    visit=visit,
    total_xs=total_xs,
    upstream=upstream,
)
# L_cell · ψ_avg − q at the converged cell average ≡ 0
np.testing.assert_allclose(residual, source, rtol=1e-12)
```

Note: the `apply` method's signature is the design judgment for
the method-implementer — see Open Questions §1 below.

---

## Gate 3 — Capability surface

**Claim layer.** Software-invariant — every capability the class
declares (`CAP_APPLY`, `CAP_SOLVE`) MUST have a passing test exercising
the matching method. Capabilities the class does NOT declare MUST
raise `MissingCapability`.

**Pillar.** Closed-form — the capability set IS the spec.

**Why this gate exists.** Per coding-elegance Pattern 4
(illegal-states-unrepresentable), the capability set is the project's
mechanism for catching "calling the wrong thing" at construction
time. A class that declares `CAP_APPLY` but whose `.apply` is broken,
or that fails to declare a capability the consumer needs, is a bug
class the test must close.

**Tests required:**

```python
def test_sncell_operator_capabilities():
    op = SNCellOperator()
    assert CAP_APPLY in op.capabilities
    assert CAP_SOLVE in op.capabilities
    assert CAP_APPLY_TRANSPOSE not in op.capabilities  # NOT in Step 1 scope

def test_sncell_operator_missing_capability_on_unsupported():
    op = SNCellOperator()
    # Attempting to access a capability not in the set MUST raise
    # MissingCapability via the composer or via direct check
    with pytest.raises(MissingCapability):
        # e.g., compose with an operator requiring CAP_APPLY_TRANSPOSE
        _ = compose_requiring_adjoint(op)
```

For `AngularRedistribution`: same shape with the M-M-specific
capability set (`CAP_APPLY` + `CAP_SOLVE`; `CAP_APPLY_TRANSPOSE`
out of scope at Step 1).

---

## Gate 4 — Curvilinear + Cartesian coverage

**Claim layer.** Software-invariant — the gates above MUST exercise
all three geometries (slab, sphere, cylinder).

**Pillar.** Closed-form. Mirrors `vv-principles` Mode 7 (MMS
simplification bias) for verification-coverage design.

**Why this gate exists.** Phase F's bug only manifested in
curvilinear (slab DD doesn't have angular redistribution; the M-M
recurrence is null). A test suite that runs only on slab cannot see
ERR-026. Step 1's gate matrix MUST hit all three geometries for
`SNCellOperator`, and BOTH curvilinear geometries for
`AngularRedistribution`.

**Concrete checklist** (per gate):

- Gate 1 (bit-identity): slab ✓, sphere ✓, cylinder ✓ (already in
  parametrize matrix).
- Gate 2 (round-trip): slab ✓, sphere ✓, cylinder ✓.
- Gate 3 (capability): geometry-agnostic; one test per class.
- Gate 5 (twin-path defense): sphere ✓ (canonical defect site),
  cylinder ✓ (for completeness; the α-dome telescoping hid ERR-026
  here pre-Phase-F).
- Gate 6 (M-M identities): sphere ✓, cylinder ✓; slab does not apply
  (no angular redistribution).

---

## Gate 5 — Phase F twin-path defense (sphere_2g_3reg sentinel)

**Claim layer.** L0 software-invariant claim. Tagged `@pytest.mark.l0`
because it pins a discretised solver invariant against a
heterogeneous multi-region multi-group probe (the canonical L0
diagnostic configuration per `vv-principles`).

**Pillar.** Closed-form algebraic identity at the cell level. The
claim is: "the apply-path consumer and the sweep-path consumer of
`SNCellOperator` produce IDENTICAL per-cell, per-ordinate ψ at the
same fixed point". The reference is the `SNCellOperator` itself
(same code path, two consumers); the claim is the **consistency** of
the two consumers, NOT correctness of either against an external
ground truth.

**Why this gate exists.** This is the explicit gate the Phase G plan
requires: "Design at least ONE test that would have caught the Phase
F bug if Step 1 had shipped in late Phase D." The Phase F twin-path
bug had `sf[0]/sf[1] = 0.522` on SI vs `1.029` on Krylov on
`sphere_2g_3reg` n=40 (Phase F memo §"Empirical evidence"). Step 1's
gate captures that asymmetry as a regression sentinel against the
ratio diverging.

**CRITICAL — what this gate proves and does NOT prove.** This is a
**twin-path consistency** test, not an absolute-correctness test.
Two consumers can agree to machine precision on a wrong fixed point
(the canonical hidden-twin failure mode). At Step 1, the test only
proves the two consumers route through the SAME promoted
`SNCellOperator`. The absolute correctness chain comes from the
existing `tests/sn/test_phase_c_crosscheck.py::test_phase_d_trajectory_resolvent_crosscheck`
+ MMS convergence tests, which are NOT in Step 1 scope.

**Gate design.**

```python
@pytest.mark.l0
@pytest.mark.verifies("dd-curvilinear-scalar")
@pytest.mark.catches("ERR-026")
@pytest.mark.parametrize("n_cells", [40])  # the Phase F diagnostic mesh
@pytest.mark.parametrize("geometry", ["sphere"])  # cylinder under separate ID
def test_apply_path_and_sweep_path_consume_same_sncell_operator(
    n_cells, geometry,
):
    """Phase F twin-path defense — apply-matvec and SI-sweep agree
    at the per-cell, per-ordinate level on sphere_2g_3reg.

    The Phase F manifestation #6 bug had sf[0]/sf[1] = 0.522 (SI) vs
    1.029 (Krylov) on this same probe. Step 1 closes this BY
    CONSTRUCTION: both call sites consume the SAME SNCellOperator
    instance with the SAME closure config, so divergence at the
    per-cell level cannot exist — they execute the same code.

    This test pins that invariant. If the two paths produce
    different per-cell, per-ordinate ψ at cell 0 of sphere_2g_3reg,
    the SNCellOperator promotion failed to unify the closure.
    """
    # Build sphere_2g_3reg n=40 (the Phase F canonical probe).
    materials = _materials_2g_3reg()
    mesh = _spherical_mesh(n_cells=n_cells, regions=3)
    quad = GaussLegendre1D.create(8)

    # Solve via SI (sweep path) to converged fixed point
    psi_si = solve_sn(materials, mesh, quad, inner_solver="source_iteration", ...)

    # Solve via Krylov (apply path) to converged fixed point
    psi_kr = solve_sn(materials, mesh, quad, inner_solver="krylov", ...)

    # Both fixed points consume the SAME SNCellOperator (post-Step 1).
    # Per-cell, per-ordinate ψ at cell 0 (the Phase F defect site) MUST agree.
    np.testing.assert_allclose(
        psi_si.psi[:, :, 0],  # (ng, M) at i=0
        psi_kr.psi[:, :, 0],
        rtol=1e-10,
        err_msg=(
            "Phase F twin-path bug regression: apply-path and SI-path "
            "produce different per-cell ψ at i=0 of sphere_2g_3reg. "
            "Step 1's SNCellOperator promotion should have unified the "
            "cell-update closure across both consumers."
        ),
    )

    # Concrete Phase F regression sentinel: sf[0]/sf[1] should NOT
    # equal 0.522 (the broken SI fixed point) and should NOT diverge
    # from the Krylov ratio.
    sf_si = (quad.weights @ psi_si.psi)  # (ng, nx)
    sf_kr = (quad.weights @ psi_kr.psi)
    ratio_si = sf_si[0, 0] / sf_si[0, 1]
    ratio_kr = sf_kr[0, 0] / sf_kr[0, 1]
    np.testing.assert_allclose(
        ratio_si, ratio_kr, rtol=1e-8,
        err_msg=(
            f"sf[0]/sf[1] differs across paths: SI={ratio_si:.4f}, "
            f"Krylov={ratio_kr:.4f}. The Phase F bug had this ratio "
            f"at 0.522 (SI) vs 1.029 (Krylov)."
        ),
    )
```

**CAVEAT (open question for method-implementer).** At Step 1, the
SNCellOperator is promoted but the call sites still live in the
twin functions `transport_operator_matvec_spherical` and
`_sweep_1d_spherical`. Whether Step 1 alone CLOSES manifestation #6
depends on whether both functions are wired to consume the same
SNCellOperator at Step 1 or only at Step 2. Per the plan,
**Step 2** is where the call-site unification happens; Step 1 is
type-system promotion only. **The test above is therefore xfail at
Step 1** (with `strict=False` since this is API-dependent), and
transitions to xpass at Step 2 when the call sites unify. Mark with
`@pytest.mark.xfail(reason="closes at Phase G Step 2 — call-site unification", strict=False)`
at Step 1.

This is the ERR-026 manifestation #7 closure marker that the plan
goal #5 demands.

---

## Gate 6 — `AngularRedistribution` algebraic identities

**Claim layer.** Foundation software-invariant — these are the
algebraic properties of the M-M recurrence (linearity, flat-flux
closure) and the Carlson seed (apply-vs-source-driven equivalence).

**Pillar.** Closed-form algebraic identities derived in
`docs/theory/discrete_ordinates.rst` §Phase F (`hebert-3-432` /
`hebert-3-434` / `hebert-3-435` / `phase-f-carlson-seed-source-driven`
labels) and the M-M weights identity (`mm-weights`).

**Why this gate exists.** The 57 existing tests at
`tests/sn/spatial/test_sweep_vs_apply_consistency.py` pin these
identities at the FUNCTION level (`carlson_inward_sweep_from_source`
helper, `CarlsonInwardSweep` strategy class). Step 1 lifts these to
the OPERATOR level — the same identities must hold when
`AngularRedistribution` consumes the Carlson seed as a composed
inner `LinearOperator`.

**Identities to pin** (mirroring the 57 existing tests):

1. **Linearity in `Q_bar`** (Carlson seed). `seed(αQ₁ + βQ₂) =
   α·seed(Q₁) + β·seed(Q₂)`.
2. **Linearity in `bc_outer_value`** (Carlson seed). Same shape on
   the BC argument.
3. **Apply-path ↔ source-driven equivalence** (Carlson seed). The
   apply-path `CarlsonInwardSweep.__call__` and the source-driven
   helper `carlson_inward_sweep_from_source` produce bit-identical
   output on matched inputs.
4. **Flat-flux closure** (M-M recurrence, Σw = 2 Hébert convention).
   On flat ψ_const with Q_bar = Σ_t · ψ_const, the recurrence's
   seed is φ̄_i = ψ_const at every cell (the Phase F Gate 1.6
   identity at `tests/sn/test_phase_c_gates.py:551+`).
5. **Linearity of `AngularRedistribution.apply` in input**. The M-M
   recurrence is affine in its scalar-flux input; the operator's
   `.apply` is linear (the affine intercept is captured in the
   Carlson seed contribution, treated as a separate term).

**Parametrize matrix:**

| geometry | n_groups | ψ_const | σ_t | bc_outer |
|---|---|---|---|---|
| sphere | {1, 2} | {0.5, 1.0, 2.5} | {0.1, 0.5, 1.5} | {0.0, 0.3, 1.0} |
| cylinder | {1, 2} | same | same | same |

**Reuse rather than rewrite.** The 57 existing tests already cover
this matrix at the function level. Step 1's tests EXTEND the existing
file (`tests/sn/spatial/test_sweep_vs_apply_consistency.py`) with
new tests that exercise the OPERATOR API
(`AngularRedistribution(...).apply(...)` / `.solve(...)`) on the same
matrix. Net add: ~15 new tests.

---

## Gate 7 — Regression-snapshot invariance

**Claim layer.** Regression contract — all 11 frozen snapshots stay
bit-identical (`assert_allclose(rtol=1e-12, atol=1e-13)`) after the
Step 1 promotion.

**Pillar.** Bit-identity contract per `vv-principles` §"Bit-identity
vs principled-equivalence". Step 1 is **pure type-system promotion**;
NO floating-point reduction tree changes; NO algorithmic changes; NO
closure-config changes. Bit-identity MUST hold throughout.

**Conditions under which a snapshot break would be principled at
Step 1.** None. Step 1 is the wrong step for a principled break.
Step 2 (closure unification) is where the SI fixed point may shift
to match the Krylov fixed point on heterogeneous curvilinear
snapshots; that shift, if it happens, is the principled break per
plan goal #7 with the structurally-independent reference being the
MMS convergence test at converged-to value.

**At Step 1, the contract is: zero snapshot regenerations.** If any
of the 11 snapshots break during Step 1, the promotion is wrong;
investigate before regenerating.

**Snapshots in scope:**

```
slab_2g_3reg_dd_n40.npz
slab_2g_homogeneous_dd_n20.npz
slab_2g_p1_aniso_dd_n20.npz
slab_fixed_source_dd_n20.npz
sphere_2g_3reg_dd_n40.npz
sphere_2g_homogeneous_dd_n20.npz
sphere_2g_p1_aniso_dd_n20.npz
cyl_1g_homogeneous_LS4_dd_n20.npz
cyl_1g_homogeneous_product_dd_n20.npz
cyl_2g_3reg_LS4_dd_n40.npz
2d_1g_LS4_dd_15x15.npz
```

Gate: `pytest tests/sn/regression/test_dd_regression.py` MUST pass
unchanged.

---

## Test specifications (file-level layout)

Three new/extended pytest modules. **Method-implementer fills the
bodies; this memo defines the structure.**

### File 1 — `tests/sn/spatial/test_sncell_operator.py` (NEW)

Module-level concerns: bit-identity + apply-vs-solve + capability
surface + curvilinear coverage for `SNCellOperator`.

Test classes (sketch):

- `TestSNCellOperatorCapabilities` — Gate 3 for SNCellOperator.
- `TestSNCellOperatorBitIdenticalSlab` — Gate 1 + Gate 4 (slab branch).
- `TestSNCellOperatorBitIdenticalSphere` — Gate 1 + Gate 4 (sphere branch).
- `TestSNCellOperatorBitIdenticalCylinder` — Gate 1 + Gate 4 (cylinder).
- `TestSNCellOperatorBitIdenticalCylindricalDegenerate` — Gate 1 +
  Gate 4 (degenerate branch, abs_mu < 1e-15).
- `TestSNCellOperatorApplySolveRoundTrip` — Gate 2 across all
  geometries.

### File 2 — `tests/sn/spatial/test_angular_redistribution.py` (NEW)

Module-level concerns: M-M algebraic identities + Carlson seed
equivalence + flat-flux closure + capability surface for
`AngularRedistribution`.

Test classes (sketch):

- `TestAngularRedistributionCapabilities` — Gate 3 for
  AngularRedistribution.
- `TestAngularRedistributionApplyLinearity` — Gate 6 #5.
- `TestAngularRedistributionFlatFluxClosure` — Gate 6 #4 (sphere +
  cylinder).
- `TestAngularRedistributionCarlsonSeedEquivalence` — Gate 6 #3
  (operator-level instead of function-level).
- `TestAngularRedistributionApplySolveRoundTrip` — Gate 2 for
  AngularRedistribution.

### File 3 — Extension to `tests/sn/spatial/test_sweep_vs_apply_consistency.py`

Add tests that exercise the OPERATOR API on the same parametrize
matrix as the existing 57 function-level tests:

- `test_sncell_operator_apply_sweep_equivalence_flat_psi` — Gate 5
  preview at flat-ψ level (round-trip the Phase F probe through
  `SNCellOperator` from both consumers).
- `test_angular_redistribution_apply_sweep_equivalence_flat_psi` —
  same for `AngularRedistribution`.

### File 4 — Phase F twin-path defense

Either add to `tests/sn/test_phase_c_gates.py` (as Gate 1.7) or to
`tests/sn/spatial/test_sweep_vs_apply_consistency.py`. The plan's
list of suggested locations names the latter; prefer that for
operator-level cohesion.

The test is `test_apply_path_and_sweep_path_consume_same_sncell_operator`
described in detail in Gate 5 above. **Mark xfail-strict=False at
Step 1** (the call-site unification happens at Step 2); the marker
self-removes at Step 2 when the test xpasses.

---

## Test-design failure-mode self-audit (Mode 7 — MMS simplification bias)

For each gate, declare which terms the test's ansatz **activates**
and which it **nulls**. If a test's ansatz nulls a load-bearing
term, flag it and add an angularly-non-trivial companion.

| Gate | Activates | Nulls | Risk |
|---|---|---|---|
| 1 (bit-identity) | All terms (synthetic per-cell inputs are arbitrary) | None — parametrize over geometry × source × ψ_in covers full per-cell algebra | LOW — synthetic-input parametrization avoids ansatz bias. |
| 2 (apply-solve round-trip) | Streaming, collision, source coupling. **Slab branch nulls angular redistribution (no upstream_state.angular)**. | Angular redistribution on slab branch. | MEDIUM — but slab is the simplest branch by design; the curvilinear branches in the same gate activate angular redistribution. |
| 3 (capability) | None — pure type-system check. | All physics. | LOW — this is a software contract gate, no physics in scope. |
| 4 (geometry coverage) | All branches when matrix is fully populated. | — | LOW — coverage IS the deliverable. |
| 5 (twin-path defense, sphere_2g_3reg) | **Heterogeneous (3-region) + multi-group (2G) + curvilinear** — full Mode 7 antidote. Activates streaming, collision, scattering, angular redistribution, BC, multi-group coupling. | None of the load-bearing terms. | **LOW** — the Phase F probe is the canonical Mode 7 antidote. |
| 6 (M-M identities) | Flat-ψ closure activates Σw = 2 normalization, seed-recurrence linearity. **The flat-ψ ansatz NULLS angular variation** (the most aggressive form of redistribution stress) by construction. | High-angular-mode redistribution. | **MEDIUM — angularly-trivial ansatz**. The 57 existing function-level tests already use the flat-ψ ansatz; reusing it at the operator level inherits the same blind spot. **MITIGATION**: gate 5 (heterogeneous multi-group sphere_2g_3reg) is the angularly-non-trivial companion. |
| 7 (regression snapshots) | All terms exercised by the 11 snapshots. | None at Step 1 (no closure change). | LOW. |

**Mode 7 net assessment.** The gate matrix passes the Mode 7 audit
because Gate 5's heterogeneous multi-group probe is the angularly-
non-trivial companion to Gate 6's flat-ψ ansatz. No gate is
unaccompanied.

**Action item for method-implementer.** Do NOT design any
additional `AngularRedistribution` test with ONLY a flat-ψ ansatz
without also providing the heterogeneous-MR companion. The 57
existing tests are exempt from this rule (their scope is
function-level Carlson-seed math, where flat-ψ IS the canonical
identity); new operator-level tests must respect the companion
discipline.

---

## V&V harness tagging discipline

Per `tests/_harness/` discipline and
`docs/testing/architecture.rst`:

| Test | Marker(s) | Rationale |
|---|---|---|
| Gate 1 (bit-identity) | `@pytest.mark.foundation` | Software invariant — not an equation-level claim |
| Gate 2 (round-trip) | `@pytest.mark.foundation` | Software invariant per LinearOperator Protocol |
| Gate 3 (capability) | `@pytest.mark.foundation` | Software invariant per Pattern 4 |
| Gate 4 (geometry coverage) | (parametrize discipline, no extra marker) | Coverage is an artifact of Gate 1–3 matrices |
| Gate 5 (twin-path defense, sphere_2g_3reg) | `@pytest.mark.l0` + `@pytest.mark.verifies("dd-curvilinear-scalar")` + `@pytest.mark.catches("ERR-026")` + `@pytest.mark.xfail(strict=False)` at Step 1 | L0 because it pins a discretised solver invariant on the canonical multi-group heterogeneous curvilinear probe; xfail at Step 1 because the call-site unification happens at Step 2. |
| Gate 6 (M-M identities) | `@pytest.mark.foundation` + `@pytest.mark.verifies(...)` for the algebraic-identity labels (`phase-f-carlson-seed-source-driven`, `phase-f-q-bar-twin-forms`, `mm-weights`) | Foundation because identities are software contracts; `verifies` because they appear in Sphinx theory page. |
| Gate 7 (regression) | (already tagged `@pytest.mark.regression` via module-level `pytestmark`) | Existing gate; unchanged. |

**Note on `@pytest.mark.verifies` for Step 1**: at Step 1 there are
no NEW equation labels in `docs/theory/`. The existing labels
(`dd-curvilinear-scalar`, `dd-mm-closure-constants`, `phase-f-*`,
`mm-weights`, `hebert-3-432`...) cover the Step 1 claims. Phase G
narrative additions to Sphinx theory pages are a separate
deliverable (per plan §"Sphinx documentation campaign"), out of
scope here.

---

## Open questions for the method-implementer

1. **`SNCellOperator.apply` signature.** The plan says:

   > `SNCellOperator.apply(cell_avg) → residual` computes
   > `L_cell · ψ_avg − q` (the discretised operator residual at the
   > cell level).

   But the cell-level residual ALSO depends on the upstream state
   and the cell's streaming terms (i.e. the `visit` packet) and the
   source `q`. The signature is therefore actually:

   ```python
   def apply(
       self, cell_avg: np.ndarray,
       *, visit: CellVisit, total_xs: np.ndarray,
       upstream_state: UpstreamState, source: np.ndarray,
   ) -> np.ndarray:
       ...
   ```

   OR alternatively the operator gets bound to the per-cell context
   at construction time:

   ```python
   bound_op = SNCellOperator(visit=visit, total_xs=total_xs, ...)
   residual = bound_op.apply(cell_avg)
   ```

   This is a design judgment that affects the round-trip gate's
   test shape (Gate 2). Recommend the second pattern for elegance
   (Pattern 5: build the right primitive; the per-cell bound
   operator is the primitive) but accept the first if the
   method-implementer determines the bound-operator pattern leaks
   into per-iteration allocation cost.

2. **`SNCellOperator.solve` vs `DiamondDifference.update` — naming
   and signature drift.** The plan says `solve(spatial_upstream,
   angular_upstream, source) → cell_avg`. But `DiamondDifference.update`
   takes `(visit, total_xs, source, upstream_state) → CellResult`.
   Does `SNCellOperator.solve` return a scalar `cell_avg`, or the
   full `CellResult` (with `outgoing_spatial_flux` and
   `outgoing_angular_state`)? Recommend **return CellResult** so
   the bit-identity contract (Gate 1) is a direct field-by-field
   `np.array_equal`. The method-implementer can decide; the gate
   design accommodates either via field-level vs whole-result
   assertion.

3. **`AngularRedistribution.apply(fi) → fi'` semantics.** The plan
   says:

   > `apply(fi) → fi'` is the M-M recurrence transform; the Carlson
   > seed becomes a private inner `LinearOperator`.

   The M-M recurrence is fundamentally a per-level loop with the
   Carlson seed as initial condition. Is the `apply` consuming
   `fi` = the scalar flux moment φ_0 of shape `(ng, nx)`, or the
   per-ordinate ψ of shape `(ng, M, nx)`? Recommend `(ng, nx)`
   φ_0 since that's the M-M recurrence's natural input (the seed
   already folds ψ → φ_0). The method-implementer can decide based
   on how the call-site consumers want to compose.

4. **`AngularRedistribution.solve` semantics.** The plan says
   `capabilities = frozenset({CAP_APPLY, CAP_SOLVE})`. The M-M
   recurrence's "solve" is back-substitution from the input scalar
   flux to the per-ordinate ψ field. Is this within scope at Step
   1, or should the capability declaration drop `CAP_SOLVE` for
   Step 1 (and add it later when the back-substitution is
   implemented)? Recommend **drop CAP_SOLVE at Step 1** to minimize
   scope; Gate 3's capability test becomes one-line simpler. The
   plan said `{CAP_APPLY, CAP_SOLVE}` but if the back-substitution
   isn't implementing anything load-bearing, dropping `CAP_SOLVE`
   is principled per Pattern 4 (illegal-states-unrepresentable —
   advertise only what works).

5. **MissingCapability test trigger** (Gate 3). To raise
   `MissingCapability`, the test needs to compose the operator with
   another operator requiring an unsupported capability. Is there
   an existing canonical fixture (e.g.
   `compose_requiring_adjoint(op)`) in
   `tests/numerics/test_operator.py` that Step 1's tests can reuse?
   If yes, refer to it. If no, the method-implementer needs to
   construct one; alternatively, simply assert
   `CAP_APPLY_TRANSPOSE not in op.capabilities` and let the
   existing composer-level tests cover the raising behavior.

---

## Test-architect risks identified

1. **Risk: bit-identity contract on Gate 1 may fail at FP ULP if the
   `apply` method's residual computation re-orders FP operations.**
   The current `DiamondDifference.update` has a documented FP-
   operation-order contract (per the module docstring lines 173–184).
   The new `SNCellOperator.apply` is NEW code with no historical
   operation order; if its construction picks a different reduction
   tree from `update`, the gates pass but the `apply ∘ solve = id`
   gate may fail at ULP level on heterogeneous inputs.

   **Mitigation**: Gate 2 asserts `rtol=1e-12` (NOT `np.array_equal`),
   which absorbs ULP-level FP-non-associativity. The `apply` method
   is allowed to choose its own reduction tree; the round-trip
   identity holds at rtol=1e-12.

2. **Risk: Gate 5 (twin-path defense, sphere_2g_3reg) is currently
   xfail at Step 1.** If the marker self-removal (per plan goal #6
   on Phase E sentinel) doesn't trigger at Step 2 because
   manifestation #7 isn't fully closed by Step 2, the gate becomes
   a quiet failure rather than a loud regression.

   **Mitigation**: use `strict=False` xfail (not strict). This is
   API-dependent (the call-site unification happens at Step 2, not
   Step 1) per the v&v-tagging idiom (see
   `feedback_vv_tagging.md`). When the gate xpasses at Step 2, the
   method-implementer removes the marker as part of the Step 2
   closeout; if the gate stays xfail, the Step 2 closeout memo
   MUST document why.

3. **Risk: Gate 6's flat-ψ ansatz inherits the Mode 7 blind spot
   from the existing 57 tests.** If the new operator-level tests
   only run the flat-ψ probe, they cannot see ERR-026 manifestation
   #7 (which requires heterogeneous multi-group). Discussed in the
   Mode 7 audit above.

   **Mitigation**: Gate 5 is the angularly-non-trivial companion.
   Method-implementer MUST NOT ship Gate 6 without Gate 5 also
   shipping.

4. **Risk: Regression-snapshot Gate 7 passes on a stale
   snapshot if a snapshot has been regenerated under the wrong
   (Phase B / pre-Phase-F-broken) state.** The plan goal #7 says
   "11 frozen regression snapshots pass at `assert_allclose(rtol=1e-12)`
   throughout the migration (bit-identity preserved unless principled
   per `vv-principles`)". The snapshots are AS-OF Phase F closeout
   `b0cc1b1`; six were regenerated by Phase F. If any of those
   regenerated values is wrong relative to the structurally-
   independent Variant-α reference, Step 1's bit-identity gate
   preserves the bug bit-identically (exactly the Phase F failure
   mode for ERR-026 in pre-Phase-F regression snapshots).

   **Mitigation**: Gate 7 is a regression contract, NOT a
   correctness gate. The absolute-correctness chain runs through
   `tests/sn/test_phase_c_crosscheck.py::test_phase_d_trajectory_resolvent_crosscheck`
   (Variant α reference) — that test MUST also pass at Step 1.
   The method-implementer's closeout memo for Step 1 must confirm
   both: (a) all 11 snapshots bit-identical, (b) Gate 4.2 cross-
   check still green at Phase E rtols.

---

## Cross-step deliverables (out of scope for this memo)

Step 1's tests prepare the gate but DO NOT cover:

- **Manifestation #7 closure verification** — happens at Step 2 via
  the call-site unification; Gate 5 transitions xfail → xpass.
- **`.H` adjoint reciprocity** — Step 5 deliverable; new test
  module `tests/sn/test_adjoint.py` is in the plan's Step 5 scope.
- **MMS L1 convergence** — already pinned by existing
  `tests/sn/l1_analytical/` suite; Step 1 doesn't change the
  numerical algebra, so existing convergence rates are preserved
  by construction.

---

## Pointers

- Plan: `.claude/plans/issue_196_phase_g_four_operator_unification.md`
  (Step 1 starts at line 76; pre-step charter at line 80).
- Reconciliation: `.claude/agent-memory/general-purpose/phase_g_four_operator_architecture_reconciliation.md`
  §2.1 (SNCellOperator verdict — IN L+C as SweepRepresentation
  strategy), §2.3 (AngularRedistribution verdict — IN L as
  curvilinear streaming closure detail).
- Phase F closeout: `.claude/agent-memory/method-implementer/issue_168_phase_f_closeout.md`
  — failure-mode profile this gate matrix defends against.
- Existing tests: `tests/sn/spatial/test_diamond.py` (28
  bit-identity tests, baseline for Gate 1's extension);
  `tests/sn/spatial/test_sweep_vs_apply_consistency.py` (57
  function-level Carlson + twin-path tests, baseline for Gate 6's
  operator-level extension); `tests/sn/test_phase_c_gates.py`
  Gate 1.6 (the canonical M-M flat-flux closure at sweep level).
- LinearOperator Protocol: `orpheus/numerics/operator.py:115-200`
  (CAP_APPLY, CAP_SOLVE, MissingCapability, Protocol contract).
- Existing source: `orpheus/sn/spatial/diamond.py` (DiamondDifference,
  3-branch geometry-polymorphic update), `orpheus/sn/spatial/pole_angular_closure.py`
  (MorelMontryAngularSweep, lines 466+),
  `orpheus/sn/spatial/psi_half_angle_seed.py`
  (CarlsonInwardSweep + carlson_inward_sweep_from_source helper).

---

## Linked memories

- `[[feedback-vv-tagging]]` — module-level pytestmark vs per-test
  decorators; xfail strict=False discipline.
- `[[feedback-diagnostic-promotion]]` — three foundation classes for
  promoted tests; this memo applies that pattern to type-system
  promotion (foundation tests for bit-identity + capability
  surface; L0 tests for the heterogeneous probe).
- `[[feedback-cross-method-protocol]]` — pairwise agreement tolerance
  rule (`max(tol_a, tol_b)`); applied at Gate 5 (`rtol=1e-10` for
  apply-path vs SI-path agreement, both at fixed-point precision).
