---
name: issue-247-moment-source-gate-spec
description: "#247 L1 verification spec — close the slope-SOURCE half of the LM-1989 trap in 2-D Cartesian LD SN (Mode-10 gap): quadrature-projected moment external+boundary source gate, mutation controls, projection-correctness sub-gate, bit-identity guard. Branch refactor/sn-foundation-cleanup. PRE-IMPL."
metadata:
  type: project
---

# #247 — 2-D Cartesian LD slope-SOURCE gate (the Mode-10 closeout)

Issue #247, L1 VERIFICATION. Branch `refactor/sn-foundation-cleanup`.
Host `.venv/bin/python -O`. I design the gate; production change (widen
`solve_sn_fixed_source` to a moment-resolved-bulk typed union + a moment
boundary trace) is owned by method-implementer. Extends
[[d5b_s3_inc_c_moment_iterate_verification]] (the FP-completion that landed
S3) + [[issue_158_ld_spatial_verification]] + [[feedback_vv_tagging]] +
[[feedback_regression_tolerance_design]]. The Mode-7 honest-scope ledger
this closes is in `SN2DCartesianLDStressMMSCase` (`mms/sn.py:1132-1161`).

## 0. THE PROBED RESULT THAT SHAPES EVERYTHING (read first)

⭐⭐ The slope-SOURCE sign error is genuinely an **O(h²)-SMALL forcing
that rides ON TOP of the O(h²) discretization floor** — confirmed live.
Average-moment L2 vs `phi_exact` under the x-slope-SOURCE sign flip:

```
AVG-moment err vs phi_exact, nc=[16,32,64,128]
  correct slope : [8.18e-3, 1.99e-3, 4.86e-4, 1.18e-4]  orders [2.04, 2.03, 2.05]
  x-slope FLIP  : [1.17e-2, 2.96e-3, 7.38e-4, 1.81e-4]  orders [1.98, 2.01, 2.03]
```

Both converge at CLEAN O(h²) (~2.0). The flipped error is ~**1.4–1.5×**
larger at EVERY mesh and the ratio is roughly CONSTANT across refinement
(1.43, 1.49, 1.52, 1.54). ⇒ **A pure convergence-ORDER leg on the
converged flux is BLIND to the sign** (order stays 2 both ways), AND **a
fixed-mesh VALUE-BAND on the average scalar flux is too fragile** (1.4×
gap — you'd need a tolerance band tighter than the factor 1.5, which the
O(h²) discretization error itself eats). This is Mode-10 in pure form;
the brief's CRUX warning is fully confirmed.

CONSEQUENCE for the gate design: the teeth do NOT come from the converged
flux. They come from TWO places where the sign flip is O(1):
1. **The consumed SOURCE moment vector** at the RHS level (structural,
   array-level — flip the production slope-source-row symbol → the source
   the kernel receives flips → O(1) red). This is the PRIMARY mutation
   control.
2. **A single-sweep per-moment value-band** vs the manufactured ψ moments
   with EXACT inflow + EXACT scattering held (one sweep isolates the
   source at O(1)): probed avg-moment ratio 2.3×, x-slope ratio 2.7× at
   nc=24 — sharper than the converged-flux 1.4×, still modest.

Per-moment ratios under the x-slope-SOURCE flip (CONVERGED, nc=32, vs
manufactured A moments): avg 3.6×, x-slope 1.6×, xy 1.0×. The AVG moment
vs `phi_exact` at a FIXED mesh is the BEST converged signal (3.6× when
referenced to the manufactured A, because the flipped average error grows
relative to the shrinking-correct error) — but it is STILL a value-band,
not a teeth-by-construction structural gate. Lead with the structural
RHS-moment mutation control; back it with the value-band.

## 1. THE PROJECTION NORMALIZATION CRUX (locked, probed)

The brief's CRUX: the projection slope-moment normalization MUST match
the production UBLD closure's moment definition. NAILED:

- The UBLD mass matrix is `M = diag(h, θh)` per axis (`_ubld.py:mass_1d`,
  `θ=1/3`), the L2-Gram of the **Legendre moment basis {P₀=1, P₁=ξ}** on
  `ξ∈[-1,1]` (`⟨P₀,P₀⟩=h`, `⟨P₁,P₁⟩=θh=h/3`).
- The cell kernel consumes a per-volume moment vector `S_moments`
  (`linear_discontinuous.py:495-500`) and forms `R_source = M @ S_moments`.
  Symbolic confirm (`ld_ubld.py:derive_d1_reduction_to_production`):
  `R_prod = [Q̄·h, θ·Q̂·h]`. So **`M` applies the volume/θ weighting; the
  projection supplies the BARE per-volume Legendre coefficients**:
  - `q̄ = (1/V)∫q` = cell average → slot 0 (`AVERAGE_MOMENT`).
  - `q̂_a = ⟨q,P₁(ξ_a)⟩/⟨P₁,P₁⟩` = the Legendre P₁ coefficient on axis a
    (for a linear `q=a+bx`: `q̂ = b·h/2`). NO `θ`, NO `h`, NO V in the
    projected coefficient — the kernel's `M` adds them.
- **d=2 Kronecker moment order = `[ψ̄, ψ̂_y, ψ̂_x, ψ̂_xy]`** (axis 0 = x
  OUTER, axis 1 = y INNER; slot order `[(0,0),(0,1),(1,0),(1,1)]`).
  Probed live: `np.diag(M) == kron(diag(h_x,θh_x), diag(h_y,θh_y))`.
  So the x-slope is **slot 2**, y-slope slot 1, xy slot 3.
- **GLOBAL frame, not sweep frame.** `_CellSolve.cell` reframes the source
  global→sweep via `octant_moment_frame_signs` (`sweep_graph.py:934`). The
  projection supplies GLOBAL-frame coefficients (the natural Legendre
  coeffs of `q(x,y)` in global x/y); production handles the per-octant
  slope-sign flip. Confirmed empirically: feeding global-frame `Qm` to
  `rep.sweep` gives correct O(h²) convergence.

The projection (q_nodes=6 Gauss) reproduces the cell AVERAGE in slot 0 to
7.7e-9 (vs a q_nodes=12 reference) and emits non-trivial slope moments
(y-slope max 0.18, x-slope 0.10, xy 0.05 on this case). The structural
independence (L11): the projection computes `∫q·P₁` by `leggauss` directly
— it does NOT call `_lift_external_source_to_moments`, does NOT call any
LD cell-op. It descends only from the manufactured closed-form `Q^ext`
(or, for the sub-gate, a hand-laid polynomial) + `numpy.polynomial`.

⚠ The crux is NOT harder than stated, with ONE caveat: `external_source`
(the existing flat producer) evaluates Q at the cell CENTRE; the
projection computes the cell AVERAGE. They differ by O(h²) (slot-0 ratio
~0.93 at nc=8). This is CORRECT (projection does genuine averaging) but
means the gate must NOT cross-check projection-slot-0 against the existing
flat `external_source` for equality — only against an independent fine-quad
cell-average reference (the sub-gate does exactly this).

## 2. DELIVERABLE 1 — the per-moment value-band gate

**Field compared**: per-moment scalar-flux moments `φ̂` (at least slot 0
AVERAGE and slot 2 x-slope; recommend ALL four) vs the manufactured
tensor-Legendre moments of `φ=A` (project A by the SAME quadrature
projector — `φ̄=A`-average, `φ̂_x=(h_x/2)∂_xA`-coeff, etc.).

**Norm**: volume-weighted L2 PER MOMENT (`_l2_2d(phi_m[...,k]-ref[...,k],
volumes)`), one error per moment slot per group. (NOTE: `volume_weighted_l2`
named in the brief does NOT exist on THIS branch — `tests/sn/_test_helpers.py`
has no such symbol; the #249 hoist was on the #236 branch per memory. Use
the local `_l2_2d` already in `test_mms_ld_2d.py`, or hoist it now.)

**Meshes / structure**: BOTH a fixed-mesh value-band AND an O(h²)
per-moment convergence leg, because neither alone is sufficient:
- **Per-moment O(h²) CONVERGENCE leg** (nc=[16,32,64]): each moment's L2
  error vs the manufactured-A moment must show order → 2 (finest > 1.85,
  all > 1.7). This proves the slope rows are CONSISTENT (the
  slope-UNKNOWN + the source together produce a 2nd-order moment). It is
  NECESSARY but (per §0) NOT sufficient for the SIGN.
- **Fixed-mesh value-band** (nc=32, the avg moment vs `phi_exact`):
  `1e-9 < err < BAND` where BAND is calibrated tight enough that the
  flipped-source error (≈2.7× the correct at nc=32 when referenced to A
  growing) is OUT of band but the correct is IN. ⚠ This band is FRAGILE
  (the gap is only ~1.5–3.6×); it is a SECONDARY guard. The PRIMARY
  sign-catcher is the structural mutation control (§3 / Deliverable 2).

**Sign sensitivity**: the value-band ALONE is too weak (Mode-10). Mark the
per-moment convergence leg `@verifies("ld-cartesian-2d")` (it upgrades the
honest scope) but DOCUMENT that its sign-teeth come from the paired
mutation control, not the band.

## 3. DELIVERABLE 2 — the mutation control (anti-pattern #11)

The PRIMARY sign-catcher. Two distinct mutations (EXTERNAL Q̂ vs SCATTERING
Σ_s·φ̂), each positive (correct → green) AND negative (flipped → red),
runtime monkeypatch reverted in `finally`. The teeth are STRUCTURAL: assert
the consumed source MOMENT VECTOR (the RHS the kernel receives) matches the
projected reference, where a sign flip is O(1) (not the O(h²) converged
flux). Recommended teeth = compare the `S_moments`/`R_source` the LD cell
op assembles for the x-slope row, OR (simpler, equally sharp) assert the
single-sweep per-ordinate x-slope MOMENT vs manufactured ψ at a fixed mesh
goes red (ratio ≥ 2 probed).

### Mutation table

| # | Source | Symbol to flip (where) | Sign-flip mutation | NEW gate must | EXISTING scalar gate (`test_ld_2d_stress_converges_second_order`) |
|---|--------|------------------------|--------------------|---------------|---------------------------------------------------------------------|
| M1 | EXTERNAL Q̂ (x-slope) | the production moment-source LIFT (the NEW `_lift_external_source_to_moments` successor that passes slope rows through) — monkeypatch it to negate the x-slope slot (slot 2) of the lifted moment source | `lifted[..., 2] *= -1` inside the lift, reverted in `finally` | go RED (single-sweep x-slope moment vs manufactured ψ ratio ≥2; OR converged avg-vs-A band breached) | stays GREEN — it feeds the FLAT source (slope rows already 0), so a flip of an already-zero slope row is a no-op. THIS asymmetry IS the Mode-10 gap being closed. |
| M2 | EXTERNAL Q̂ (y-slope) | same lift, slot 1 | `lifted[..., 1] *= -1` | go RED | stays GREEN (flat → no-op) |
| M3 | EXTERNAL Q̂ (xy) | same lift, slot 3 | `lifted[..., 3] *= -1` | go RED (weaker — xy ratio ~1.0 converged, so use the single-sweep moment teeth or the RHS-vector structural assert) | stays GREEN |
| M4 | SCATTERING Σ_s·φ̂ (slope) | `ScatteringOperator._combine_iso_aniso` / the iso accumulator that scatters `Σ_s⊗I` over EVERY spatial moment (`scattering.py:527-544`) — monkeypatch to negate the slope rows (slots 1,2,3) of the iso source | `iso.values[..., 1:] *= -1` (post-accumulation), reverted in `finally` | go RED | stays GREEN — the scalar gate's converged flux is only ~1.4× sensitive; the scattering-slope flip needs the per-moment single-sweep teeth |

⚠ M1–M3 (external) are only CONSUMED once the production widens the lift to
pass slope rows. Until then they are xfail-strict (the lift zeros them → the
flip is a no-op → both green → the test that asserts "flipped reddens" fails
→ xfail). M4 (scattering) is consumed TODAY (the S3 iterate path) — but its
teeth still need the per-moment single-sweep gate (the converged flux is
sub-floor sensitive). DISTINGUISH the two clearly: M1–M3 verify the NEW
external consumption; M4 verifies the EXISTING scattering consumption was
never sign-blind. Both are required by the brief.

PREFERRED structural teeth (sharper than the value-band, O(1) by
construction, recommended for M1/M4): build the projected moment source,
feed it to ONE `rep.sweep` with EXACT manufactured-ψ inflow + EXACT
scattering `Σ_s·φ̂` (from manufactured A moments) held, and assert the
per-ordinate x-slope OUTPUT moment `ang[...,2]` matches the manufactured
ψ x-slope moment to `rtol` — then the monkeypatch flip breaks it by ratio
≥2 (probed 2.7×). This isolates the source at O(1) and does not wait for
SI convergence (fast + sharp).

## 4. DELIVERABLE 3 — projection-correctness sub-gate (L11, foundation)

A `@pytest.mark.foundation` test (NO `verifies` — it's a software
invariant on the projector, not a solver claim). Project a KNOWN bilinear
polynomial `q = a00 + a10·x + a01·y + a11·x·y` on a cell `[xL,xR]×[yL,yR]`
and assert the four tensor-Legendre coefficients EXACTLY (machine
precision, q_nodes≥2 integrates a bilinear exactly):

```
q̄    = a00 + a10·xc + a01·yc + a11·xc·yc          (slot 0)
q̂_y  = (hy/2)·(a01 + a11·xc)                       (slot 1)
q̂_x  = (hx/2)·(a10 + a11·yc)                       (slot 2)
q̂_xy = (hx/2)·(hy/2)·a11                           (slot 3)
```

Probed: all four match to machine precision. The hand-derived coefficients
are the structurally-independent reference. The convention (`(h/2)·∂`
coefficient on `{1,ξ}`, `ξ∈[-1,1]`) MATCHES `mass_1d`'s `diag(h,θh)`
Legendre basis — this is what makes the projection apples-to-apples with
production.

**HOW IT AVOIDS TAUTOLOGY**: (a) the projector code uses only
`numpy.polynomial.legendre.leggauss` + the hand-derived integral, NEVER
`_lift_external_source_to_moments` nor any `_ubld`/`LinearDiscontinuous`
method; (b) the reference is hand-laid polynomial algebra, not a production
echo. Branch 1 (`derive_2d_cartesian_ld_stress_mms`) emits only the
per-ordinate `Q_closed`, NOT a symbolic spatial moment — so the
"numeric-matches-symbolic-moment" leg the brief offers as optional is NOT
available here (Branch 1 has no symbolic moment emitter). The polynomial
hand-check is the L11 ground. (If a symbolic moment emitter is added to
Branch 1 later, add the numeric==symbolic leg then — flagged for
archivist/method-implementer.)

A SECOND sub-gate leg: project the manufactured A and assert slot 0 ==
cell-average computed by an INDEPENDENT fine quadrature (q_nodes=12)
to ~1e-8 — pins that the projector's average IS the cell average, not the
cell-centre value (the §1 caveat: the existing flat `external_source` uses
cell-centre, which differs by O(h²) — do NOT cross-check against it).

## 5. DELIVERABLE 4 — bit-identity guard (DriftWarning)

The DEFAULT flat-source path (slopes zeroed) MUST stay bit-identical after
the typed-union widening. Production today: `_lift_external_source_to_moments`
(`solver.py:1896`) lifts a flat `(N,ng,*spatial)` onto slot 0, zeros the
rest; DD/Step (per_axis==1) returns input unchanged. The widening adds a
moment-resolved-bulk BRANCH; the flat branch must not move.

Gate: the existing strict DriftWarning regression gate. Run
```
pytest tests/sn/sweep/core tests/sn/solve \
  -W "error::tests.sn.regression._regression_assert.DriftWarning"
```
(the canonical strict-clean gate from [[d5_nd_polymorphism_verification]]
D5a.3 — baseline DD strict 505p/1s/4xf). Plus the existing LD bit-id
two-paths `test_ld_2d_two_paths_ffw_equals_mfw` (random vacuum source) +
`test_ld_2d_stress_two_paths_ffw_equals_mfw` MUST stay GREEN (the flat +
boundary-average path unchanged). ALSO add a NEW explicit DD bit-id leg:
solve a DD 2-D fixed-source with a flat bulk through the WIDENED entry and
assert `array_equal` vs the pre-widening snapshot (DD per_axis==1 → the
moment branch is never taken → byte-identical). The flat-LD path
(slope-zeroed) must also stay bit-identical: assert that passing a flat
`(N,ng,nx,ny)` bulk through the widened entry gives the SAME result as
today (the slope rows are still zeroed unless a moment-resolved bulk is
explicitly supplied — the honest default).

## 6. DELIVERABLE 5 — legacy-behavior pin

**Probed live**: passing a moment-resolved bulk `(N,ng,nx,ny,2^d)` to
`solve_sn_fixed_source` TODAY raises
`ValueError: fixed-source bulk shape (24,2,8,6,4) does not match
(N, ng, *spatial) = (24,2,8,6)` (`solver.py:1876`).

**Grepped**: NO existing test asserts this `ValueError` (searched
`tests/` for "fixed-source bulk shape" / "does not match (N, ng" — zero
hits in SN). ⇒ There is NO legacy ValueError pin to migrate. The
production change RELAXES this check (typed union: flat bulk OR
moment-resolved bulk OR `TimedFullField`); the only migration is that the
relaxed check must STILL reject a genuinely-wrong shape (e.g.
`(N,ng,nx,ny,5)` where 5≠2^d, or `(N,ng,nx)` on a 2-D mesh). Recommend a
NEW negative pin: assert the widened entry STILL raises on a
trailing-moment-axis length ≠ `2^d` (the moment-vector contract), so the
relaxation doesn't swallow real shape bugs.

## 7. DELIVERABLE 6 — scope of legs (external Q̂ vs boundary face-slope)

**Recommend TWO legs (separate gates), not one.** The external Q̂ and the
boundary transverse-face-slope are independent production paths with
different magnitudes and different consumption sites:

- **Leg A — external Q̂** (the dominant, primary): mutations M1–M4 above.
  The bulk slope source. This is what `_lift_external_source_to_moments`
  must widen.
- **Leg B — boundary transverse-face-slope** (secondary): the manufactured
  inflow ψ varies ALONG each face (transverse). `_inflow_to_moments`
  (`loss_representation.py:357`) zeros the `2^{d-1}` transverse face-slope.
  Probed: the transverse face-slope is ~12% of the face-average on this
  case (a per-cell-along-face O(h) quantity). Separate gate: project the
  manufactured inflow onto the transverse face Legendre basis, feed the
  moment-resolved boundary trace, and assert the near-boundary cells'
  flux improves / the converged value-band tightens vs the
  face-average-only trace. Mutation: flip the transverse-face-slope sign.

Keep them SEPARATE because: (a) Leg B needs a DIFFERENT production change
(moment-resolved boundary trace, not the bulk lift); (b) Leg B's signal is
smaller (boundary-localized) — bundling would let Leg A's stronger signal
mask a Leg B regression; (c) the honest-scope ledger names them as two
distinct deferred items. xfail-strict each until its production path lands.

## 8. -O safety (Mode 8) + file location

The canonical invocation is `python -O` (confirmed: pytest warns "assert
statements are not executed"). The NEW gate file MUST use `np.testing.*` /
`pytest.fail` only — NO bare `assert` (it would be a false green under -O).
Mirror the existing `test_mms_ld_2d.py` style.

**File**: extend `tests/sn/verification/mms/test_mms_ld_2d.py` (the D5b LD
2-D home) — the new gates are the S4 closeout of the same case. The
projection-correctness foundation sub-gate could also live in
`tests/sn/spatial/test_linear_discontinuous.py` (kernel-level) but the
moment-SOURCE consumption is end-to-end → keep with the MMS gates.

## 9. Honest-scope docstring update (after production lands)

`SN2DCartesianLDStressMMSCase` (`mms/sn.py:1132-1161`) and the existing
`test_ld_2d_stress_converges_second_order` docstring both say the
slope-SOURCE half is DEFERRED (S4). When the production widening + these
gates land, UPDATE both to say the slope-SOURCE sign IS now VERIFIED (the
mutation control reddens). This is a `verifies`-tagging update
([[feedback_vv_tagging]]): the new per-moment gate carries
`@verifies("ld-cartesian-2d")`; the mutation controls carry
`@catches(...)` ONLY if a real production bug surfaces (then mint the ERR;
next free is ERR-063). Until then no ERR (Mode-10 is a proactive-gap close,
not a caught bug).

## 10. Baseline (live @ HEAD, branch refactor/sn-foundation-cleanup)

```
pytest tests/sn/verification/mms/test_mms_ld_2d.py -q -m "not slow"
  → 5 passed, 2 deselected  (35s)
```
ERR ceiling = ERR-062 (the pure-z moment-broadcast, already gated here).
`LinearDiscontinuous.theta == 1/3`, `spatial_basis_per_axis == 2`.
NEVER run all `tests/sn` (#212 hang at
`test_keff_slab::test_heterogeneous_absolute_keff`); scope to
`tests/sn/verification/mms` + `tests/derivations`.

## 11. Self-improvement (Mode-10 reinforced)

This task is a SECOND clean instance of Mode-10 (activated-but-unconstrained
term): the slope-SOURCE rows are genuinely consumed (S3 scattering iterate)
yet a sign flip leaves the converged flux at O(h²) and ~1.4×. The vv Mode-10
row already exists and names THIS exact instance (#240 D5b-S4) — NO new row
needed; the spec REINFORCES the existing entry's recipe (mutation-verify
the activated term is also constrained; if the converged-value mutation is
sub-floor, add an O(1)-isolating companion). The companion here = the
structural RHS-moment mutation control + the single-sweep per-moment teeth.
This is the canonical resolution of a Mode-10 gap; no skill edit required.
