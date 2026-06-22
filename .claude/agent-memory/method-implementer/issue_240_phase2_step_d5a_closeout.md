---
name: issue-240-phase2-step-d5a-closeout
description: #240 Phase 2 Step D5a — fold the 2-D ScanMarch's inline DiamondDifference onto the scheme coefficient model (cartesian_scan_coefficients + reflect_scan_coefficients + residual_kernel_batch); principled ~1-ULP re-baseline; scheme owns the 2/w; zero inline DD
metadata:
  type: project
---

# #240 Phase 2 Step D5a — ScanMarch onto the coefficient model (#239) — CLOSEOUT

Branch `feature/sn-space-angle-tier2`, 2026-06-16, **NOT committed** (main
agent runs elegance-enforcer + qa on the diff, then commits).

**Why:** ScanMarch (the 2-D Cartesian row-march loss-representation) hard-coded
DiamondDifference's `alpha/beta/D_row` math inline, so it ONLY worked for DD and
no scheme could ride the row-march polymorphically. D1-D4 homed the generic
advection–reaction coefficient ops onto `DiscretizationSchemeBase`. D5a routes
ScanMarch's two interior kernels THROUGH those ops, making the row-march
scheme-generic over every `transverse_coupling_is_facewise` closure (DD today;
Step rides it free once Step exists). **PRINCIPLED ~1-ULP RE-BASELINE** of
existing DD math (precedent #158-B1): converged VALUES unchanged, only the
arithmetic path re-associates.

**How to apply:** when a sweep/matvec body inlines a scheme's discretization
constants (the diamond `2`, the blend `0.5`), the fold is to add a scheme
**coefficient producer** that owns those constants and have the body consume the
generic reconstruction ops. The litmus (coding-elegance Pattern 2 / the
coefficient-model litmus): *if an explicit-matrix representation would have to
re-derive a calculation the sweep does, that calculation belongs on the scheme.*

## The mechanism (the fold)

Two new scheme methods (NOT in the `DiscretizationScheme` Protocol — they are
OPT-IN capabilities guarded by the `transverse_coupling_is_facewise` trait,
carried on `DiscretizationSchemeBase` with raising defaults, exactly like
`affine_scan_coefficients` / `cell_kernel_batch`):

- `DiscretizationSchemeBase.cartesian_scan_coefficients(*, s_scan, s_transverse, sig_t)`
  → `(a, inverse_denom, w, transverse_couplings)`. The multi-D scan-march
  analogue of `affine_scan_coefficients` (the 1-D/curvilinear chain producer):
  `s_scan`/`s_transverse` are RAW down-face streaming `g = |μ|/Δ`; the scheme
  applies its diamond `2`. DD override: `S = Σ_t + 2g_x + 2g_y`,
  `a = 2(2g_x)·inverse_denom − 1`, `w = ½`, `transverse_couplings = (2g_y,)`.
- `DiscretizationSchemeBase.reflect_scan_coefficients(psi_bar)` → `(α, β)`. The
  apply-direction x-face reflection scan: `α = −(1−w)/w`, `β = ψ̄/w`. DD override:
  `α = −1`, `β = 2ψ̄`. (The recurrence form of `outgoing_face_from_average`.)

**SOLVE (`ScanMarch._sweep_interior`):** per y-row, `cartesian_scan_coefficients`
→ `(a, inverse_denom, w, (c_y,))`; effective source `QV_eff = Q + c_y·ψ_y_in`;
`beta = source_emission(QV_eff, inverse_denom, w)`; `alpha = a`; `_scanmarch_row`
closes via `cell_average(in_x, out_x, w)` + `outgoing_face_from_average(ψ̄, ψ_y, w)`.

**APPLY (`ScanMarch._loss_action_interior`):** per y-row, reconstruct in_x off
the probe via `reflect_scan_coefficients` reflection scan; the per-cell residual
+ the transverse-y outflow ride `scheme.residual_kernel_batch(psi_bar, psi_in=(in_x, ψ_y_in), s_axes=(g_x, g_y), sigt, Q=0)`
(the ÷V uniform matvec kernel every facewise closure shares — the apply twin of
the scan was ALREADY established as `residual_kernel_batch` territory in
`diamond.py:531-543`; D5a just routes ScanMarch through it).

`_scanmarch_row` (scan.py) gained a `w` parameter (was hardcoded `0.5`); now
`cell_average`/`outgoing_face_from_average` consume it.

## Files touched

- `orpheus/sn/spatial/scheme.py` — NEW `cartesian_scan_coefficients` +
  `reflect_scan_coefficients` raising-default methods on
  `DiscretizationSchemeBase` (after `affine_scan_coefficients`, ~`:937-1066`);
  a NOTE block in the `DiscretizationScheme` Protocol (~`:520-533`) recording
  WHY the scan-family methods are base-only opt-in capabilities, NOT Protocol
  members (the Protocol declares ONLY `update`/`residual`).
- `orpheus/sn/spatial/diamond.py` — NEW DD overrides `cartesian_scan_coefficients`
  (~`:522-574`) + `reflect_scan_coefficients` (~`:576-595`).
- `orpheus/sn/loss_representation.py` — folded `ScanMarch._sweep_interior`
  (~`:1328-1356`) + `ScanMarch._loss_action_interior` (~`:1433-1471`); both now
  read `self.mesh.scheme`, zero inline `2.0*`/`Diamond`/`0.5`.
- `orpheus/sn/spatial/scan.py` — `_scanmarch_row` (~`:312-352`) gained `w` param;
  closes via the generic staticmethods (was inline `0.5*(in+out)` + `outgoing_face_from_average(...,0.5)`).
- `tests/sn/operators/test_streaming_operator.py` — EXTENDED
  `TestT4bPreT4RegressionSnapshot` with `_assert_cart2d_arm` + 3 cart2d apply
  arms (D5a.2; `~:875-960`).
- `tests/sn/solve/test_affine_carve_bit_identity.py` — re-baselined the 4 2-D
  GOLDEN sha hashes + regen-history comment (D5a entry, 2026-06-16).
- `tests/sn/_fixtures/wave_t_t4/pre_t4_snapshots.npz` — regenerated 3 cart2d
  `*_apply_bulk` keys to post-D5a value (boundary byte-identical; all 47 keys
  preserved).
- `tests/sn/_data/bc_extraction_2d_baseline/vacuum_bulk_2d_seed{0,1,2}.npy` —
  regenerated via `--capture-baseline` (2-D vacuum matvec ~1-ULP re-baseline).
- `docs/theory/loss_representations.rst` — NEW stub
  `.. _loss-rep-scanmarch-coefficient-model:` with `.. todo::` (archivist #240
  D6) + STALE flag on `loss-rep-scanmarch-solve`/`-apply` (those eqs still show
  the inline-DD form) + `:mod:`/`:meth:` cross-refs.
- `docs/verification/matrix.rst` — auto-regen (foundation 3093→3097, +4 tests).

## Gate results (verbatim summaries)

- **D5a.1** (ScanMarch ≡ FullFieldWavefront oracle, nulp): `11 passed`.
  Existing parametrize ALREADY covers the Mode-9 degeneracy-break — the
  `(12, 7, 6, 2, "reflective")` tuple is non-square / ng=2 / specular (reflective),
  plus `(5, 9, 4, 4, "reflective")` non-square 4g. No ADD needed.
- **D5a.2** (cart2d matvec frozen reference + 3 slab arms): `6 passed`. The
  pre-fold cart2d matvec was VERIFIED bit-identical (0 ULP) to the frozen
  pre-#240 snapshot, so the frozen IS the correct pre-D5a reference; the fold
  re-associates ~1 ULP (relΔ ~3e-14) which the `reduction_depth=nx` nulp metric
  AMPLIFIES at near-zero cancellation cells (1 entry @ 256 ULP at |ψ|~0.014,
  absΔ~4.4e-16 = 1 ULP of the O(1) field). **Regenerated the cart2d `*_apply_bulk`
  snapshots to the post-D5a value** (same-build → 0 ULP henceforth, exactly like
  the slab arms), documented as the D5a re-baseline. Boundary byte-identical.
- **D5a.3** (negative control, strict DriftWarning escalation):
  `513 passed / 1 skipped / 4 xfailed` — IDENTICAL pre vs post fold. The 1-D
  slab sha (`si_slab_2g_het`) UNCHANGED (the 1-D CumprodScan scan path is
  byte-identical — the fold did not leak into the 1-D solve). The 2-D arms
  (`si_2d_p1_aniso_het` + `krylov_2d_p1_aniso_het`) RE-BASELINED — see below.
- **Stay-green anchors:** 2-D MMS (DD value reference) + kinf homogeneous (≥2G):
  `39 passed, 3 xfailed`. ScanMarch end-to-end + 2-D windowing: `8 passed`.
- **Route-around** (`tests/sn/operators spatial sweep/core sweep/cartesian_2d`):
  `1044 passed` (+3 new cart2d arms), `7 failed` — ALL 7 confirmed PRE-EXISTING
  on the clean baseline (sphere 1-D matvec SPH #206-family ×3, `Face 'ymin'
  requires genuine mu_y` ×2, sphere curvilinear apply snapshots ×2).
- **Audit** exit 0. **Sphinx `-W -E`** exit 1 on the PRE-EXISTING `paramref`
  role error in `orpheus/geometry/mesh.py` (confirmed identical on baseline —
  NOT a D5a regression; my `loss-rep-scanmarch-coefficient-model` stub builds
  warning-free).

## THE LOAD-BEARING FINDING — the spec's D5a.3 byte-identity premise was OFF

The spec (`d5_nd_polymorphism_verification.md` §D5a.3) claimed "the SOLVE-path
scan snapshots stay byte-identical (`2·QV·inv ≡ QV·inv/0.5`)" and treated a RED
solve snapshot as a fold-leaked-into-solve BUG. **That premise did not hold for
the 2-D SOLVE.** Reason: the pre-D5a `_sweep_interior` was the ONE remaining
place using `÷D_row` **division** (`α = 2sx2/D_row`, `β = 2(...)/D_row`), NOT the
coefficient-model `×inverse_denom` **reciprocal** that the 1-D CumprodScan
already rides. The spec's `≡` argument assumed `×inv` was already in use. Folding
the 2-D solve onto the coefficient model switches `÷D → ×(1/D)`, a genuine
~1-ULP FP-association change — the SAME sanctioned re-baseline the 1-D path took
at #158 Inc B. So `si_2d_p1_aniso_het` (SI, solve fold) and
`krylov_2d_p1_aniso_het` (apply fold through GMRES) BOTH re-baselined; the 1-D
slab byte-pin stayed exact (the true negative control).

**Verified the re-baseline against the vv-principles 3-criteria** (not just
old-vs-new proximity):
1. Principled: named coefficient-model ops, inline DD removed.
2. Structurally-independent reference: (a) post-fold 2-D SI `.solve` φ ≡ Krylov
   `.apply` φ to **4.19e-12 rel** (DIFFERENT code paths converging to the same
   fixed point), AND (b) the D5a.1 ScanMarch≡FullFieldWavefront oracle pins the
   value to analytical `k_inf=1.875` / `φ=Q/Σ_t`.
3. FP-non-associativity: SI = iter_count × ULP; Krylov = GMRES_tol × ULP.

This is consistent with the affine-carve file's OWN history — it already
re-baselined the krylov-2D hashes at #240 Phase 2 Step B for the same apply
re-association. D5a just extends that to the SI-2D arm (the solve fold) too.

## Zero-inline-DD audit (D5a.4 — confirmed by audit, xfail placeholder SKIPPED)

`ScanMarch._sweep_interior` + `._loss_action_interior` + `_scanmarch_row`:
CLEAN — no inline `2.0*`, no `Diamond` reference, no hardcoded `0.5` in CODE
(only in docstrings as `:meth:` cross-refs and explanatory `w = ½` notes). The
diamond `2` and blend `w` live entirely in the scheme. Polymorphic-readiness is
satisfied by construction (Step rides the row-march with zero further ScanMarch
change once `Step.cartesian_scan_coefficients`/`reflect_scan_coefficients` land).

## Design decisions

- **`cartesian_scan_coefficients` is a SEPARATE method from `affine_scan_coefficients`,
  not a generalization.** `affine_scan_coefficients` is fundamentally 1-D-chain-
  shaped (curvilinear geometry args `abs_mu/A_down/A_total/dA_w/c_out/V`, NO
  transverse slot). The row-march needs the transverse `2g_y` in the diagonal +
  the scan-axis RAW `g`. A d-generic tuple signature (`s_transverse` as a tuple)
  keeps it ready for d=3 (the row would sum multiple transverse couplings).
- **Both methods are base-only opt-in (NOT Protocol).** First attempt added them
  to the `DiscretizationScheme` Protocol → broke 2 `@runtime_checkable`
  conformance tests (synthetic strategies with only `update`/`residual` failed
  `isinstance` on 3.12+ member-presence checks). FIXED by moving them to the base
  with raising defaults + a Protocol NOTE block. This mirrors `affine_scan_coefficients`
  / `cell_kernel_batch` (also base-only). The Protocol stays minimal: the
  universal `update`/`residual` contract only.
- **APPLY routes through `residual_kernel_batch` (reuse), SOLVE needs a new
  scan producer.** The apply direction has a KNOWN probe ψ̄ → per-cell residual
  is independent → the existing uniform ÷V matvec kernel applies directly. The
  solve direction is x-coupled (cell i inflow = cell i−1 outflow) → cannot use
  the independent-batch kernel → needs the closed-form affine scan coefficients.

## Owed (DISPATCH_REQUEST to archivist queued)

- The rich-narrative expansion of `loss-rep-scanmarch-coefficient-model` (#240
  D6): rewrite `loss-rep-scanmarch-solve`/`-apply` eqs in coefficient-model
  form, re-derive the byte/ULP posture. The stub carries the full brief + the
  STALE flag + the 3-criteria re-baseline evidence.

## Lesson → coding-elegance / algebra-of-record

The "byte-identical scan" claim in a coefficient-model docstring is only true
WHERE the consumer already uses the `×inverse_denom` reciprocal form. A consumer
still on `÷denom` division (a leftover inline path) re-baselines ~1 ULP when it
joins the coefficient model — that's a FEATURE (the model's single source of
truth), not a regression, but it WILL trip a strict byte/sha gate. **When folding
a leftover-inline path onto an established coefficient model, expect a ~1-ULP
re-baseline of any strict byte gate downstream, and verify it against the model's
structurally-independent reference (here: the FullFieldWavefront oracle + the
SI≡Krylov cross-check), not against old-vs-new proximity.** The 1-D vs 2-D
division/reciprocal asymmetry is the kind of convention drift a Pattern-7
crosswalk would surface at design time.
