---
name: Phase D Gate 1.1 sphere MMS diagnosis — Carlson seed lives in the M-M angular recurrence, NOT the WDD spatial sweep
description: Numerics-investigator diagnosis for Issue #168 Phase D Step 2. Gate 1.1 sphere MMS failure under MorelMontryAngularSweep is the M-M angular recurrence's hardcoded ψ_{1/2,i}=0 seed; the Hébert §3.9.4 (3.432)-(3.435) Carlson coupled-pole sweep should seed the ANGULAR recurrence, not the WDD spatial pole-face IC.
type: project
---

# Phase D Gate 1.1 sphere MMS diagnosis — closeout

**Branch**: `refactor/sn-operator-algebra` 2026-05-12. Phase D Step 2
of the plan at `/home/vscode/.claude/plans/structured-booping-parrot.md`.

**Diagnostic script**: `tests/sn/diagnostics/gate_1_1_sphere_mms_failure.py`
— self-contained CLI probe; run with `python tests/sn/diagnostics/gate_1_1_sphere_mms_failure.py`.

**Headline finding**: the Phase D hypothesis ("Carlson coupled-pole
inward μ = −1 sweep closes Gate 1.1") is empirically CONFIRMED on the
flat-ψ probe, BUT the **point of injection differs from what the
literature memo + Phase D plan implied**. The literature memo's §7
implementation note routes the inward-sweep result into the WDD
**spatial** pole-face IC (`psi_face_in` at the outward sweep's first
cell). That intervention CHANGES NOTHING for flat ψ. The intervention
that closes the residual is seeding the M-M **angular** recurrence's
half-angle face flux `ψ_{1/2,i,g}` with the inward sweep result —
i.e., the `psi_half_left = 0` hardcode at
`orpheus/sn/spatial/pole_angular_closure.py:411` is the offending
line, NOT `psi_face_in = fi[..., 0]` at
`orpheus/sn/operator.py:738`.

This is a corrective amendment to the literature memo's §7 architectural
sketch, not a falsification of the Hébert §3.9.4 math. The math is
correct; the production code's misplacement of "the Carlson seed" is
the bug.

## 1. Failure profile

**Problem setup**: nx=10, R=2.0, GL4 quadrature, sphere, reflective
BC, homogeneous Σ_t. Packed flat-ψ probe `psi = ones(n_unknowns)`.
Expected: `L · ψ = Σ_t · ψ` (pure collision residue, redistribution
and streaming cancel per ordinate on flat ψ).

**Baseline residuals (Phase C-shipped pole-angular closure)**:

| Closure | Σ_t = 0 | Σ_t = 0.5 |
|---|---|---|
| `LegacyTauSymmetricInterpolation` | **1.78e-15 PASS** | **1.78e-15 PASS** |
| `BaileyFlatFluxRedist`             | **1.78e-15 PASS** | **1.78e-15 PASS** |
| `MorelMontryAngularSweep`          | **1.89e+01 FAIL** | **1.89e+01 FAIL** |

The M-M failure is independent of Σ_t — the residual is the same
1.89e+01 at both values because it lives in the redistribution
term (Σ_t cancels in `r = L·ψ − Σ_t·ψ`).

**Failure profile (M-M, Σ_t = 0.5)**:

* **Per-ordinate**: all 4 GL4 ordinates show similar magnitude residual
  (12.9, 18.9, 18.3, 11.9). The largest residuals are on the
  small-|μ| ordinates (n=1, 2) where the redistribution coefficient
  ΔA/w is largest.
* **Per-cell**: residuals peak at i=0 (the pole cell, max|r| = 18.9)
  and decay roughly as `1/r_c` as r grows. At the outer cell i=9 the
  residual is still 1.28. The bug is NOT localised — it pollutes
  every cell, but with magnitude controlled by the geometry factor
  ΔA/w (which is larger near r=0).
* **Per-(ordinate, cell) field** (g=0, shape (4, 10)):
  ```
  [[ 12.9   5.5   3.4   2.4   1.9   1.6   1.3   1.1   1.0   0.  ]   ← μ = -0.86
   [-18.9  -8.1  -5.0  -3.6  -2.8  -2.3  -1.9  -1.7  -1.5   0.  ]   ← μ = -0.34
   [ 18.3   7.9   4.8   3.5   2.7   2.2   1.9   1.6   1.4   1.3]   ← μ = +0.34
   [-11.9  -5.1  -3.1  -2.2  -1.7  -1.4  -1.2  -1.1  -0.9  -0.8]]  ← μ = +0.86
  ```
  Sign pattern: alternates per-ordinate (+, −, +, −). Magnitude
  decreases with r. The two outermost-cell entries in the inward
  ordinates are exact-zero — because the eq_map excludes those
  positions (BC fixes them) and the residual field is built only
  on unknowns. No 0 entries in the outward ordinates' outermost
  cell (those ARE unknowns).

**Reading the fingerprint**: this is NOT spatial-streaming-localised
(if it were, residuals would be concentrated at one boundary). It is
NOT an iteration-scheme bug (it shows up on a single matvec, not
across iterations). The sign-alternation per ordinate paired with the
1/r_c falloff is the signature of a **per-ordinate constant offset in
the redistribution term** that propagates across cells via the
α-cascade — exactly the prediction of a wrong `ψ_{1/2,i}` seed for
the M-M angular recurrence.

## 2. Carlson intervention sweep (the empirical evidence)

Four interventions tested against the M-M failing configuration:

| Intervention | What it changes | max\|residual\| |
|---|---|---|
| `[A]` Carlson seed for outward-sweep `psi_face_in` (literature memo §7) | `orpheus/sn/operator.py:738` `psi_face_in = phi_aux[g,0]` | **1.89e+01 FAIL** (unchanged) |
| `[B]` Carlson seed for M-M half-angle ψ_{1/2,i} | `pole_angular_closure.py:411` `psi_half_left = phi_aux` | **1.78e-15 PASS** |
| `[C]` BOTH `[A]` + `[B]` | both sites | **1.78e-15 PASS** (no extra effect) |
| `[D]` M-M half-angle ψ_{1/2,i} = ψ_cell of ordinate 0 (broadcast) | `pole_angular_closure.py:411` `psi_half_left = ψ_cell` | **1.78e-15 PASS** |

**Reading the results**:

* `[A]` confirms the WDD spatial pole-face IC is NOT what's wrong.
  The Phase C `psi_face_in = fi[:, outgoing_mask, 0, 0]` Lewis-Miller
  seed is already the cell-centre value, which on flat ψ already
  EQUALS the Carlson inward-sweep output (both = 1.0). Replacing it
  changes nothing.
* `[B]` is the canonical Carlson intervention: feeding `phi_aux[g,i]`
  (the FULL spatial profile from the inward μ = −1 sweep) into the
  M-M recurrence's `psi_half_left` seed closes the residual to
  machine precision.
* `[C]` is bit-identical to `[B]` — confirms `[A]` is a no-op on this
  probe.
* `[D]` is a sanity-check: for flat ψ where `phi_aux[g,i] ≡ ψ_cell`,
  any seed equal to ψ_cell works. The flat-ψ probe **cannot
  distinguish** Carlson seed [B] from naive cell-centre seed [D].

## 3. Structural-independence cross-check (vacuum BC)

The flat-ψ-on-reflective probe is a degenerate fixed point: `[B]`,
`[C]`, `[D]` all coincide because the inward sweep returns ψ_cell
exactly. To prove the Carlson seed is canonical (not merely
coincidentally correct), I added a vacuum-BC probe where
`phi_aux ≠ ψ_cell`.

**Vacuum BC inward sweep result** (R=2.0, nx=10, Σ_t=0.5):
```
phi_aux = [0.613, 0.572, 0.527, 0.478, 0.423, 0.362, 0.295,
           0.220, 0.138, 0.048]
```
vs `ψ_cell = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]` — the seeds differ by
up to 0.95 in absolute value.

**Vacuum-BC residuals on flat ψ** (NOT a fixed point of the operator
— residual is non-zero by design; we compare relative magnitudes):

| Seed | max\|residual\| |
|---|---|
| Baseline (ψ_{1/2}=0) | 23.98 |
| `[B]` Carlson inward sweep | 12.40 |
| `[D]` ψ_cell broadcast | 12.92 |

Seed `[B]` is **strictly better** than seed `[D]` on the vacuum-BC
probe, **and** the residual fields differ by max-abs 7.31 (i.e. they
are structurally distinct outputs). This is the structural-
independence evidence: the Carlson seed `[B]` is **not** a
coincidental match for the naive seed `[D]`; it is mathematically
distinct and quantitatively superior.

The remaining 12.40 residual on vacuum BC is expected — flat ψ = 1 is
NOT a fixed point of the operator under vacuum BC (the BC removes
inflow, so the cell-centre values cannot stay at 1.0 in the actual
solution). The flat-ψ probe under vacuum BC is a probe of the
*operator action*, not its converged solution.

## 4. Diagnosis: where the bug lives in production code

Production file: `orpheus/sn/spatial/pole_angular_closure.py`
Function: `_mm_weighted_angular_recurrence_single_level`
Line: 411

Production code (verbatim):

```python
psi_half_left = np.zeros((ng, nx), dtype=psi_level.dtype)
for m in range(M):
    tau_m = tau_level[m]
    psi_half_right = (
        psi_level[:, m, :] - (1.0 - tau_m) * psi_half_left
    ) / tau_m
    redist[:, m, :] = (
        dAw_level[:, m].reshape(1, nx)
        * (alpha_level[m + 1] * psi_half_right
           - alpha_level[m] * psi_half_left)
        / volume.reshape(1, nx)
    )
    psi_half_left = psi_half_right
```

The literal `psi_half_left = np.zeros(...)` is the wrong Carlson seed.
Per Hébert §3.9.4 Eqs. (3.432)-(3.435), the correct seed is the
inward μ = −1 spatial sweep's cell-centred output `φ̄_{1/2,i}` —
which IS a function of (Σ_t, source moments, BC), NOT identically
zero.

**Note on the production docstring**: line 52-58 of
`pole_angular_closure.py` already documents the right answer:

> initialised at the Carlson zero-weight starting direction
> μ = −1 (Hébert Eqs. 3.432-3.435 give the source-driven sweep that
> fixes φ_{1/2,i} for the **inverse** problem; for the **forward
> apply** matvec we adopt φ_{1/2,i} = 0, the unique choice that makes
> the recursion's seed consistent with α_{1/2} = 0 and that the sweep
> converges to under fixed-point iteration).

This is the design rationale that Phase B baked into the code. The
empirical evidence in this diagnosis shows that rationale is **wrong
for the apply matvec**: the M-M recurrence requires `φ_{1/2,i}` to
equal the source-driven inward-sweep output, not 0, for the operator's
flat-ψ fixed point to satisfy `L·ψ = Σ_t·ψ`. The `0` seed makes the
operator inconsistent with its sweep partner — the converged sweep
result has `φ_{1/2,i}` ≠ 0 (it's the value the sweep iterates to from
the BC trace law), and the apply matvec must use that same value to
be a faithful linearisation of the discrete operator.

## 5. Failure-mode characterisation per `vv-principles`

This bug is **Mode 3 — Missing factor / wrong term initialization**.
The hardcoded `psi_half_left = 0` is the equivalent of a missing
`ΔA/w` factor — both are "wrong term initialization with
1-group/homogeneous-flat invisibility".

Why it survived Phase C tests:
- Gate 1.1 with `LegacyTauSymmetricInterpolation` and `BaileyFlatFluxRedist`
  pass because both closures use the cell-centre value as the
  half-angle face flux (Legacy collapses τ-symmetrically on flat ψ,
  Bailey collapses by definition). For flat ψ the cell-centre ≡
  φ̄_{1/2,i} ≡ 1, so the half-angle "face flux" is consistent with
  the spatial profile.
- Cylindrical Gate 1.1 with M-M closure passes because each
  μ-level has its own α-dome that telescopes back to zero at the
  level's edges. The pole-face seed `psi_half_left = 0` cancels
  cleanly in the cylindrical case via the level-by-level α-dome
  cancellation (literature memo §8 open question 3 hypothesis 2 is
  empirically confirmed: the cylindrical structure absorbs the wrong
  seed; the sphere does not).
- Phase B's L1 flat-flux-identity test passed because (per the
  pole_angular_closure docstring §"Three strategies") that test pins
  the equivalence of the three closures **on flat ψ**, where the
  per-cell DD angular recurrence's ψ_{1/2,i} = 0 seed plus the
  α-cascade plus the redistribution combine in a way that makes the
  three closures give *equal* (but, in M-M's case, also wrong)
  output for flat ψ on a particular boundary. The test isn't checking
  that any of them satisfy `L·ψ = Σ_t·ψ` — only that they agree with
  each other.

## 6. Confidence in the Carlson approach

**Confidence: HIGH** that the Carlson-coupled-pole inward sweep is the
right structural fix, BUT with one architectural amendment to the
Phase D plan:

* The `CarlsonCoupledPole` strategy should expose its output as the
  **half-angle face flux seed for the M-M angular recurrence**, not
  the spatial pole-face IC for the outward WDD sweep.
* The Protocol shape proposed in the Phase D plan
  (`PoleFaceInitialCondition` consumed at `operator.py:734-740`) is
  the WRONG injection point per intervention `[A]`. The right
  injection point is `_mm_weighted_angular_recurrence_single_level`'s
  `psi_half_left` initialisation.
* The shape contract should be `(ng, nx)` per the literature memo
  §7 output shape correction — NOT `(ng, n_outgoing)`. This is
  important: the M-M recurrence consumes a per-cell value, not a
  per-outgoing-ordinate value.

**Architectural recommendation**: the `CarlsonCoupledPole` strategy
should live INSIDE `_mm_weighted_angular_recurrence_single_level` (or
as a hook the recurrence calls), NOT as a separate Protocol consumed
by the matvec. The `PoleAngularClosure` Protocol is the right home
for the seed; the new `PoleFaceInitialCondition` Protocol the Phase D
plan proposes is a misnomer — it should be `MMHalfAngleInitialCondition`
or similar.

For cleanliness: the Phase D plan's two-Protocol architecture (Option
B) can still be honoured if `MorelMontryAngularSweep.__call__`
accepts an optional `psi_half_seed` argument and delegates it to
`_mm_weighted_angular_recurrence_single_level`. The
`PoleFaceInitialCondition` Protocol then carries the same
`(ng, nx)` output but the matvec routes its output INTO the
`pole_angular_closure` call as the seed kwarg, not into `psi_face_in`.

## 7. Side observations for the method-implementer

1. **The literature memo §7's `(ng, n_outgoing)` shape is wrong** —
   per the empirical evidence above, the seed shape is `(ng, nx)`,
   one value per radial cell. The memo's own §4 ("Pole-face shape vs
   cell-centre values") flags this concern; the diagnostic confirms
   the cell-centred profile is the load-bearing shape.

2. **Caching opportunity is real**. The inward sweep's output
   `phi_aux[g, i]` depends on `(Σ_t, Q_moments, bc_outer, mesh)`. In
   the M-M recurrence's calling context (one call per matvec, called
   from `transport_operator_matvec_spherical`), the inputs are
   essentially constant across Krylov inner iterations: Σ_t, mesh, BC
   are problem-fixed; `Q_moments` are the source moments which change
   per outer iteration but NOT per Krylov inner. Cache the result at
   the SNStreamingOperator level keyed on outer-iteration state.

3. **The M-M recurrence currently uses `psi_half_left = 0` for the
   CYLINDRICAL case too**, but cylindrical Gate 1.1 passes empirically.
   This is because the cylindrical `_mm_weighted_angular_recurrence_single_level`
   is called per μ-level, and each level's α-dome ends at α=0 by
   antisymmetry, so the wrong `psi_half_left = 0` cancels at the
   level boundary. The bug exists structurally on both geometries
   but its observable consequence is only sphere-side. The Phase D
   fix should still update the cylindrical call sites for
   structural correctness (architectural alignment with the canonical
   form), but cylindrical tests will be a regression-stability check,
   not a new pass.

4. **The fact that intervention `[D]` (ψ_cell broadcast) ALSO works
   on flat-ψ reflective should be in the method-implementer's V&V
   plan as an explicit FALSIFICATION probe**. On a vacuum-BC or
   non-flat-ψ MMS the [B] vs [D] seeds give different residuals
   (max-abs difference 7.31 on the diagnostic's vacuum probe). The
   production unit test for `CarlsonCoupledPole` MUST use a
   non-degenerate probe (vacuum BC, or non-flat ψ, or both) so it
   distinguishes the canonical Carlson seed from the naive ψ_cell
   broadcast. A flat-ψ-reflective test will not catch a future
   regression that replaces `phi_aux` with `ψ_cell.broadcast` —
   they are coincidentally equivalent on that probe.

5. **The Phase D plan's Step 3c instruction** (replace
   `psi_face_in = fi[:, outgoing_mask, i0, 0].copy()` with the new
   strategy) is INCORRECT per the diagnostic evidence. The site that
   needs the strategy is the M-M angular recurrence at
   `pole_angular_closure.py:411`, NOT the matvec site at
   `operator.py:738`. The Phase D plan should be amended to reflect
   this.

## 8. Acceptance gate this diagnostic establishes

The diagnostic SHOULD be promoted in two forms after Phase D ships:

* **Promote-as-pytest** the failure-mode signature: a per-ordinate
  flat-flux test under M-M closure on reflective sphere that asserts
  `max|residual| ≤ 1e-12`. This is the Gate 1.1 MMS xfail-strict
  marker removal (Phase D plan §5a). It already exists as a
  parametrised xfail in `tests/sn/test_phase_c_gates.py`; Phase D
  removes the xfail marker.

* **Promote-as-pytest** the structural-independence safeguard: a
  vacuum-BC test that pins `[B] − [D] ≠ 0` (i.e. the Carlson seed
  produces a residual field that is structurally distinct from a
  naive ψ_cell-broadcast seed). This is a NEW gate the diagnostic
  identifies; the Phase D plan should add it as a regression test in
  `tests/sn/spatial/test_pole_face_initial_condition.py` (or
  whatever the new test module is named).

## 9. Phase D plan deviations summary

| Plan claim | Diagnostic finding |
|---|---|
| Inject seed at `operator.py:734-740` `psi_face_in` | Wrong site — `[A]` no-op on flat ψ |
| New Protocol named `PoleFaceInitialCondition` | Wrong name — should be `MMHalfAngleInitialCondition` (or merge into `PoleAngularClosure`'s call signature) |
| Output shape `(ng, n_outgoing)` | Wrong shape — should be `(ng, nx)` per cell, not per ordinate |
| Phase D Step 3c modifies the matvec | The matvec is fine; the modification should be in `_mm_weighted_angular_recurrence_single_level` |

All four are CORRECTIONS to the plan, not falsifications of the
core Hébert §3.9.4 Carlson hypothesis. The literature memo's math
is correct; the plan's *application* of that math to the production
code architecture mis-identified the injection point.

## Pointers

* **Diagnostic script**: `tests/sn/diagnostics/gate_1_1_sphere_mms_failure.py`
* **Production code site**:
  `orpheus/sn/spatial/pole_angular_closure.py:411`
  (`psi_half_left = np.zeros(...)`)
* **Existing Phase C empirical test**:
  `tests/sn/test_phase_c_gates.py::test_apply_curvilinear_per_ordinate_flat_flux_residual`
  parametrised over `mms` × `sphere` × Σ_t — currently xfail(strict=False).
* **Literature memo this corrects**:
  `.claude/agent-memory/literature-researcher/phase_d_carlson_coupled_pole.md`
  §7 (Implementation notes — shape and injection point)
* **Phase D plan this amends**:
  `/home/vscode/.claude/plans/structured-booping-parrot.md` §Step 3c
  (matvec injection-point)
* **Phase C closeout this builds on**:
  `.claude/agent-memory/method-implementer/issue_168_phase_c_closeout.md`
* **`vv-principles` reference**:
  Mode 3 (Missing factor / wrong term initialization) +
  Anti-pattern 4 (Homogeneous-only verification) — the bug survived
  Phase B's L1 flat-flux-identity test because that test compared
  closures against each other on flat ψ rather than against the
  closed-form `L·ψ = Σ_t·ψ` identity.
