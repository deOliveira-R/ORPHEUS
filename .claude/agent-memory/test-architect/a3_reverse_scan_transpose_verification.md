---
name: a3-reverse-scan-transpose-verification
description: Reusable recipe — verifying the #280 reverse-SCAN transpose-solve (SweepOperator.apply_transpose) + the 2.5a apply-loop unification, WITH the assembled-Mᵀ oracle, in the post-#226 typed-predicate operator world. Cited by the A3/#280 spec update.
metadata:
  type: project
---

# Verifying the #280 reverse-scan transpose-solve (Phase 2.5, post-assembly world)

Full worked method behind the A3/#280 gate-spec extension (design-time,
2026-07-04). The compressed lesson is `lessons.md` L17; this file keeps
the reusable RECIPE. Cited by `.claude/plans/archive/a3_solve_transpose_verification.md`.

**Why:** Phase 2.5 of the stencil-assembly campaign = GitHub #280 —
unify the SN loss-operator sweep machinery into ONE orientation×kernel
walk; `sweep_transpose` + a shared coherent adjoint frame emerge. The
existing a3 gate spec (G1–G6) predates THREE things: the #226 taxonomy
merge (frozenset→typed-predicate + reified inverse operators), Phase
2's assembly mode, and the Phase 2a module relocation. This is the
reconciliation.

**How to apply:** any time the reverse-DAG transpose-solve, the 2.5a
apply-loop unification, or the assembled-Mᵀ oracle is touched, start here.

## 0. The operator-terminology DRIFT (surface this FIRST — the a3 spec is stale)

The a3 spec's `InvertibleOperator.solve_transpose` + `CAP_SOLVE_TRANSPOSE`
+ `_AdjointOperator.solve` + `CAP_SOLVE` vocabulary is **RETIRED** by the
#226 taxonomy merge (memory L13/L15 — the frozenset→typed-predicate
terminal step landed). The CURRENT #280 surface (verified in code
2026-07-04):

| a3-spec name (retired) | CURRENT surface (post-#226) | where |
|---|---|---|
| `InvertibleOperator.solve_transpose` (the reverse-scan primitive) | **`SweepOperator.apply_transpose`** — `SweepOperator = (L+C).inverse()`; the reverse-scan IS the inverse operator's transpose matvec `(L+C)⁻ᵀb` | `sweep_operator.py` (delegates like `.apply`→`inner.solve`) |
| `CAP_SOLVE_TRANSPOSE` string tag | **`SweepOperator.is_adjointable` → True** (runtime property; today inherits base `False` from `InverseWrapMixin`) | `operator.py:538` |
| `_AdjointOperator.solve` + `CAP_SOLVE` | **NOTHING NEW** — the metric transpose-solve `A.inverse().H.apply(b)=G⁻¹·SweepOperator.apply_transpose(G·b)` falls out of the EXISTING `_AdjointOperator.apply` (already routes `inner.apply_transpose`) FOR FREE once `SweepOperator.apply_transpose` exists | `operator.py:1127` |
| the `.H.solve` coherence | **`_AdjointOperator.inverse()` + `is_invertible`** closing `A.H.inverse() ≡ A.inverse().H` (predicate honesty: `inner.is_invertible and inner.inverse().is_adjointable`) | #280 comment |

`InverseWrapMixin` docstring already names #280: *"the adjoint-inverse is
the #280 family, deferred (free for the iterative branch, a reverse-DAG
`sweep_transpose` for the direct sweep)"*. **Design gates against the
CURRENT surface; do NOT resurrect the retired string tags.**

## 1. Current paths (post Phase-2a relocation — re-derive by grep, these drift)

- **Production** `orpheus/sn/loss_representation/__init__.py`: `loss_action_transpose`
  (1-D adj matvec) `:2777`; `_apply_walk` (1-D fwd matvec) `:2426`; `_run`
  (1-D solve SCAN) `:2994`; `reversed(cells)` `:2903`; `mirror=quad.reflection_index("x")`
  `:2833`; `angular_adjoint` `:2943`. Multi-D Cartesian adjoint still `raise
  NotImplementedError` `:2822`.
- **Kernel** `orpheus/transport/spatial/`: `cell_balance_for_streaming`
  `cell_balance.py:120`; `_cartesian_streaming_diagonal` `diamond.py:363`;
  `residual_kernel_batch` `diamond.py:494`(DD)/`linear_discontinuous.py:646`(LD).
- **2.5 rename estate** `orpheus/sn/spatial/` = {pole_angular_closure,
  psi_half_angle_seed, sweep_cache, scan, pairing} (R9).
- **Assembly** `orpheus/sn/loss_representation/assembly.py` — `assemble_ordinate_blocks`,
  `ordinate_walk_order`; Cartesian-only.
- **Canary tests** (post-relocation): `tests/sn/sweep/core/test_sweep_cache.py`
  (was `tests/sn/spatial/`); `tests/sn/operators/{test_g_adjoint_reciprocity,
  test_invertible_operator,test_one_octant_walk,test_capability_survival,
  test_removal_form_matvec_sweep}.py`; `tests/sn/sweep/core/{test_cell_kernel_batch,
  test_phase_c_gates}.py`; `tests/sn/sweep/test_assembly_mode.py`.

## 2. The 2.5a apply-loop unification — bit-identity of a RELOCATION (both orientations)

2.5a folds `{_apply_walk (fwd matvec), loss_action_transpose (adj matvec)}`
into ONE orientation-parametrized per-cell loop frame, **bit-identical for
BOTH orientations**.

**Canaries that pin each today** (both 0-ULP `array_equal`):
- fwd `_apply_walk`: `test_removal_form_matvec_sweep::test_invertible_apply_is_M_of_C_sigma_bit_identical`
  (slab/sphere/cyl/cart2d 2G) + `test_invertible_operator` Q/Σ_t + `test_cell_kernel_batch`.
- adj `loss_action_transpose`: `test_removal_form_matvec_sweep::test_invertible_apply_transpose_is_M_transpose_of_C_sigma_bit_identical`
  (slab/sphere/cyl 2G, NOT cart2d — deferral) + `test_phase_c_gates::test_apply_apply_transpose_reciprocity_under_sweep_frame`
  (sphere/cyl, rel<1e-12) + `test_g_adjoint_reciprocity` (slab/sphere/cyl 1g+2g).

**THE TRAP — these `array_equal` canaries are SELF-REFERENTIAL.** The
reference is the SAME `loss_action[_transpose]` on a FRESH `StreamingOperator(σ_r)`
instance. They pin "the composite override OWNS the matvec (routes through
the representation) vs the leaf-sum leak" — a reduction-tree discriminator,
NOT a value oracle. Under a RELOCATION that moves BOTH `op.apply`'s path AND
the reference's path (SAME relocated code), they move TOGETHER and stay
GREEN even if values shifted. NECESSARY, NOT SUFFICIENT for the relocation.
**The sufficient addition:** a FROZEN pre-carve baseline snapshot
(`tests/sn/regression/_regression_assert.py::assert_regression`,
`--capture-baseline`) of the fwd+adj 1-D matvec VALUES (slab/sphere/cyl,
2G), re-compared post-carve at `array_equal` — structurally independent of
the relocated code because captured BEFORE it. This is the genuine 0-ULP
relocation proof.

**Structural pin — there is NO 1-D walk spy today** (unlike multi-D's
`test_one_octant_walk.py`: runtime `_interior_walk` spy + AST tripwire on
`{is_solve,is_apply,is_matvec}`). 2.5a CREATES the shared 1-D loop frame, so
it needs its OWN pin:
1. **Runtime spy** (Mode-11-sharpened WRAP, autouse/monkeypatch counter):
   BOTH orientations route through the ONE shared loop method; count > 0.
2. **AST tripwire** — orientation SIBLING of the `is_solve` rule: forbid
   `is_adjoint`/`is_forward`/`is_transpose`/`is_reverse` as real identifiers
   (Name/Attribute/arg) in the shared frame; demand an orientation-carrying
   OBJECT (like `_ApplyOperands`/`_SolveOperands`/`_SweepEmit`). Verified:
   none of those identifiers exist in the module today — the tripwire starts clean.

## 3. The assembled-Mᵀ oracle — a NEW structurally-independent reference (Cartesian only)

Phase 2's assembly gives `M` = per-(ordinate,group) sparse block of `L+C`,
walk-order LOWER-triangular (`permuted = block.as_matrix()[ix_(order,order)]`).
The existing fwd G2 is `solve_triangular(permuted, q[order], lower=True) ≡
(L+C).solve(q)` (LAPACK `dtrtrs` ≡ sweep; the #284 discharge). The
**transpose-solve oracle** for `sweep_transpose`:
```
x_ref = solve_triangular(permuted.T, b[order], lower=False)  # upper-tri back-sub
# then un-permute: sweep_transpose(b)[order] ≡ x_ref
```
`SparseAssembledOperator.apply_transpose` (exact CSR `M.T @ x`,
`is_adjointable`=True) is the matvec analog.

**L11 classification (honest):** the assembled blocks are extracted from
`residual_kernel_batch` UNIT PROBES — INDEPENDENT of the reverse-walk code
(the `apply_transpose`/`sweep_transpose` bodies), but SHARING the forward
KERNEL. So:
- **Catches (that dense-apply G2 doesn't):** wrong TRANSPOSED-SCAN-COEFFICIENT
  (a'/b' for the reversed affine recurrence) — LAPACK `trtrs` shares NO code
  with the scan; PLUS the EXACT triangularity certificate; PLUS discharges the
  "reverse-substitution IS sweep_transpose" contract at object level (the
  transpose analog of #284).
- **Misses (that dense-apply G2 catches):** everything CURVILINEAR. Assembly is
  Cartesian-ONLY (`SNMesh.streaming` gate; `test_curvilinear_refuses_the_cartesian_walk`).
  The sphere μ-reversal/mirror bug lives ONLY in curvilinear → **the SPHERE
  leg leans on G2's dense-apply oracle `np.linalg.solve(M.T,b)` (M from fwd
  `apply` on basis vectors), which stays the KEYSTONE.**
- **Jointly blind (G1+G2+assembled):** a bug in the SHARED forward kernel
  (`residual_kernel_batch`/`cell_balance_for_streaming` denom) that moves
  `apply`, `apply_transpose`, AND `sweep_transpose` TOGETHER — all three
  references build from that one kernel. Closed by the FORWARD-correctness
  ground (`test_removal_form` array_equal + `test_invertible_operator` Q/Σ_t +
  reciprocity) + the one-source teeth (`test_assembly_mode::test_teeth_shared_kernel_sign_flip`).
  Mode-11-adjacent (twin-path/coverage), NOT Mode-12.
- **LD caveat:** the 1-D adjoint may be DD/scalar-ONLY (`loss_action_transpose`
  allocates `psi_bar=(ng,N,nx)` no moment tail — fragmentation-map landmine #2).
  If so, `sweep_transpose` inherits DD-only, and the LD assembled-Mᵀ oracle is
  MOOT for 2.5b (no LD reverse to compare). VERIFY the LD-adjoint scope first.

**Worth a gate row?** YES on Cartesian (slab + 2-D); it earns L2-cross-check
status. **Reserve `keff`/scalar functionals — Mode-12:** NEVER credit
`sweep_transpose` on `k*` or a norm (`eig(Kᵀ)=eig(K)` is adjoint-blind by
construction; memory L16 "NEVER keff(asm)≡keff(apply)"). Object-level (full
field) gates ONLY.

## 4. Reverse-SCAN failure modes the a3 reverse-LOOP spec did NOT cover

The a3 spec assumed a reverse-LOOP (it was written when `loss_action_transpose`
— a per-cell reverse loop — was the reference). 2.5b builds the reverse-SCAN
coherent with `_run` (Blelloch prefix scan). NEW modes:

| Failure mode | Catcher |
|---|---|
| Wrong transposed scan coefficients (a'/b' of the reversed affine recurrence — scan ≠ loop) | G2 dense-Mᵀ + assembled-Mᵀ (LAPACK trtrs), non-symmetric het+≥2G config |
| Two-denom seam: reverse-scan rides ÷V `residual_kernel_batch` (apply) instead of ×V `affine_scan_coefficients` (the `_run` form) — same slab denom two ways (#242) | G1 round-trip (denom mismatch breaks I) + G2 |
| Curvilinear per-μ-level ordinate-coupled scan REVERSAL + `angular_adjoint` threading + Carlson mirror-seed cross-sweep | G1/G2 SPHERE leg (slab NULLS — degenerate curvilinear) — unchanged keystone |
| Relocation moves op.apply AND its self-referential reference together | FROZEN pre-carve baseline snapshot (§2), NOT the array_equal self-reference |
| Orientation as boolean flag | AST tripwire `is_adjoint`/`is_forward`/`is_transpose` + orientation OBJECT (§2) |
| sweep_transpose credited on k*/norm | Mode-12: object-level full-field gate only |

## 5. Scaffold order (pre-carve canary vs post-carve proof)

- **S0 PRE-CARVE (harden the surface the carve rebuilds):** G3 full-loss
  `(L+C−S−B)` reciprocity — extend `test_g_adjoint_reciprocity` from `L+C−B`
  to `L+C−S−B` (S† landed `15185e5`; asymmetric SigS ≥2G; per-group one-hot φ,
  vv L27). Gives 2.5a's adjoint bit-identity a COMPOSITE-level canary. + the
  FROZEN baseline snapshot (§2).
- **2.5a apply-loop unification (bit-id RELOCATION):** frozen-baseline array_equal
  (both orientations) + existing removal-form array_equal + NEW 1-D spy + AST
  tripwire. Forward canaries stay green.
- **2.5b sweep_transpose reverse-scan (NEW values → re-baseline NOT bit-id):**
  G1 round-trip + G2 dense-Mᵀ (slab+SPHERE, keystone) + assembled-Mᵀ
  (Cartesian slab+2D) + the §4 mutations RED.
- **2.5c operator wiring (CURRENT surface):** `SweepOperator.apply_transpose`
  routes to `sweep_transpose`; `SweepOperator.is_adjointable`→True;
  `_AdjointOperator.inverse()`+`is_invertible` (coherence). Gates: Mode-11
  routing sentinel (`A.inverse().H.apply` executes `SweepOperator.apply_transpose`,
  not fwd apply) + value round-trip via forward-matvec reciprocity (the #280
  comment's Gate-1, never calls the transpose path) + swap-law
  `A.H.inverse().apply(b)≡A.inverse().H.apply(b)` + predicate flips
  (`SweepOperator.is_adjointable` True; `InverseOperator`/`GreenOperator` STAY
  False — no consumer) + `assert_type(A.H.inverse(), LinearOperator[D,C])`.

## 6. #282-fix gate plan (CONDITIONAL — scope not user-ruled)

If fix route (a) lands inside 2.5 (carry ψ(·,μ=−1) per μ-level as explicit
AUGMENTED state → the spherical walk becomes genuinely triangular over the
augmented system):
- The existing characterization `test_282_characterization_spherical_seed_is_a_back_edge`
  is a loud-flip: it asserts `above > 1e-12*scale`; when the back edge vanishes
  it REDS with an actionable "rewrite as spherical G2" message (L16 loud-flip
  contract).
- **Successor:** a spherical G2 gate (triangularity + LAPACK ≡ sweep on the
  AUGMENTED system). Augmented structure: the μ=−1 zero-weight Carlson
  starting-direction seed rows (Hébert 3.432–3.435) sit FIRST in each μ-level
  block; the per-ordinate cell blocks march after in μ-increasing order; the
  current BACK edge (ordinate→seed at mirror) becomes a within-augmented-state
  FORWARD edge (seed→first-ordinate inflow via the α-dome recursion).
- **Teeth:** swap the augmented seed's coupling direction (feed the LAST
  ordinate) → reintroduces the back edge → triangularity leg reds; OR a sign
  flip in the Hébert 3.432–3.435 starting-direction solve → LAPACK≡sweep leg reds.
- Curvilinear assembly is out of 2b scope precisely BECAUSE of this back edge;
  the fix is what UNBLOCKS a curvilinear assembled-Mᵀ oracle.

## 7. The CYLINDER arm (#280 2.5b-cyl, task #28) — three extra structures, ONE mandatory config

The cyl reverse-scan (`_run_transpose` cyl branch, retiring the guard at
`__init__.py:4480-4490`) transposes the forward `_run` cyl path, which has
THREE structures the sphere arm lacks: (1) MULTI-LEVEL M-M thread
(`quad.level_indices`, each independent — `psi_angle_bar` re-inits per level,
sphere is single-level `[None]`); (2) the m0 DIRECT-SEED FOLD transpose (the
2.5b-cyl-fwd landmark: non-carrying, `seed_cot=None`/`m_seed=None`, `−mm_a_in·
out_bar` routes to `psi_avg_chain_p_bar` NOT `psi_a_in_chain_bar`, no
`ang_coeff·s_bar` feedback, folded `(c_out−c_in)` coeffs); (3) DEGENERATE
pure-azimuthal ords (μ_r=η=0) via the reversed per-cell `dag_walk`+`scheme.update`
path — the ERR-066 catcher family.

**THE MANDATORY CONFIG — `Quadrature.product(n_mu=4, n_phi=8)` cyl (empirically confirmed 2026-07-05):**

| structure | product cyl | level_symmetric cyl |
|---|---|---|
| degenerate ords (μ_r=0) | **8** (ords 2,6,10,…) | **0** |
| seed-ord M-M `mm_a_in_coeff[m0]` | **1.0 (LIVE)** | **0.0 (DEAD)** |
| seed fold c_in[m0] | 0.14–0.48 (live) | 0.0 (no-op) |
| levels / carrying | 4 / non-carrying | 4 / non-carrying |

`level_symmetric` NULLS **both** of the two hardest transpose terms (degenerate
+ seed-fold M-M). §0.6 in the flesh — every prior cyl gate rode `level_symmetric`
and was BLIND (only the `cyl_product_2g` reciprocity row caught the 2.5a ERR-066
adjoint-drop). **Product cyl is mandatory; LS is the everything-nulled CONTROL
(a mutation that reds product MUST stay green on LS — that asymmetry IS the
config-blindness proof).** Both mesh helpers already exist in
`tests/sn/sweep/test_cyl_direct_seed_fold.py` (`_cyl_product_mesh`,
`_cyl_level_symmetric_mesh`).

### Gate roster (extend `tests/sn/operators/test_loss_transpose_solve.py` — the a3 spec's `tests/sn/sweep/` path is stale)

Add `"cyl_product"`, `"cyl_ls"` to `_MESHES` (cyl `space is None` → G1/G2 helpers
reduce to the slab bulk-only branch out-of-the-box), plus new gates:

- **G1 round-trip (bulk):** `solve_transpose∘apply_transpose = I` + dual. Fast
  smoke vs the independent `loss_action_transpose` (adj matvec, cyl-verified
  post-ERR-066). NECESSARY-not-sufficient.
- **G2 dense-Mᵀ oracle (bulk — THE structural keystone):** `solve_transpose(b)
  ≡ np.linalg.solve(M.T, b)`, `M = _probe_augmented_matrix_one_group` (column-
  probes the FORWARD apply — a THIRD independent reference; M invertible for cyl,
  the forward gate `np.linalg.solve(M,q)` proves it). Cyl `space is None` →
  `_source_carried_mask` all-True (NO seed-outflow-corner mask, unlike sphere —
  the forward gate compares the full bulk). Catches (a)-(d) on the bulk.
- **G3 full-field reciprocity (bulk⊕boundary — NEW, covers the #284 slot):**
  `⟨A.solve(q),p⟩ = ⟨q,A.solve_transpose(p)⟩` summing ALL FullField slots.
  **Validated 1.8e-16 on the landed sphere; bulk-only form fails at 42% → the
  boundary carries load-bearing IP mass, so a boundary-cotangent (`m_boundary`
  passthrough) bug has O(1) teeth here that G1/G2-bulk MISS.** Structurally
  independent via the independently-verified forward `solve` (forward cyl gate).
  Add parametrized over ALL meshes (retroactively hardens slab/sphere too).
- **G4 m_seed=None + R12a-refuse contract:** cyl `_run_transpose` returns
  `m_seed is None` (`starting_direction is None` on the output); passing a
  `seed_cot` raises `ValueError` (the `:4389` elif). Mirrors the forward R12a.
- **G5 config-activation sentinel (foundation):** assert product-cyl `n_degenerate
  > 0` AND `mm_a_in_coeff[m0] != 0`; LS-cyl `n_degenerate == 0` AND
  `mm_a_in_coeff[m0] == 0`. Pins the mandatory config's activation so a future
  quad change can't silently decay it into blindness (§0.6 / Mode-7 declaration).

**assembled-Mᵀ (LAPACK `solve_triangular`): NOT available for cyl — Cartesian-ONLY**
(`SparseAssembledOperator`/`ordinate_walk_order` gate on `SNMesh.streaming`; memo
§3, `test_curvilinear_refuses_the_cartesian_walk`). The dense-Mᵀ G2 is the cyl's
keystone in its place — exactly the sphere's situation. LD moment-tail transpose
= the #280 kernel-pair deferral (the `:4398` `cell_moment_count>1` guard raises,
DD/scalar-only, same as slab/sphere).

### Mutation → gate matrix (post-carve teeth-verification; each RED under `-O`; ALL value gates np.testing/pytest.fail)

Revert every mutation by **manual re-edit / in-process monkeypatch — NEVER `git
checkout`** (the `_run_transpose` file is UNCOMMITTED during the surgical carve;
process-discipline.md). All are code edits to the new cyl branch:

| # | mutation | reds | config-asymmetry |
|---|---|---|---|
| **(a)** seed-ord M-M misroute: `−mm_a_in·out_bar` → `psi_a_in_chain_bar` (the sphere/bulk arm) instead of `psi_avg_chain_p_bar` | G2, G3 (+G1) | RED product / **GREEN LS** (mm_a_in_coeff[m0]=0 nulls it) |
| **(a2)** seed-ord keeps the `ang_coeff·s_bar` feedback (ang_contrib was 0) OR unfolded coeffs not the `(c_out−c_in)` folded pair | G2, G3 | RED product / GREEN LS (dead fold) |
| **(b) LOAD-BEARING** degenerate-ord DROP: skip the `if geom.is_degenerate: <reverse dag_walk>; continue` branch (ERR-066 reincarnation) | G2, G3 (+G1) | **RED product / GREEN LS (0 degenerate) — the ERR-066 config-blindness, run BOTH to demonstrate** |
| **(c)** multi-level thread leak: hoist `psi_angle_bar=zeros` OUT of the per-level loop | G2, G3 (+G1) | RED product & LS (both 4-level) |
| **(d)** coupled-pole mirror reverse: `pole_outflow_bar[global_n]` not `[mirror[global_n]]` (or drop the `pole_outflow_bar` re-init per level) | G2, G3 | RED product |
| **(e)** m_seed contract: return a non-None `StartingDirectionSourceSink`, OR drop the `:4389` R12a refuse | G4 | pytest.raises / assert-None |

G2/G3 are the mutation keystones (structurally independent — forward apply /
forward solve). G1 is corroborating smoke. Mutation (b) MUST be run on BOTH
configs to positively exhibit the product-RED/LS-GREEN asymmetry (the evidence
that product is the mandatory config, not a preference).

### Scaffold order (cyl 2.5b)

- **PRE-CARVE:** the forward references are ALREADY landed+green (forward solve =
  `test_cyl_direct_seed_fold.py`; adj matvec = `test_removal_form_matvec_sweep`
  cyl). No new value to freeze (2.5b is NEW values, re-baseline NOT bit-id).
  Optional loud-flip: a `pytest.raises(NotImplementedError)` char-test on the
  current cyl guard, removed at carve time (churn — skip unless the main agent
  wants the L16 characterization). Confirm slab/sphere keystone + forward-cyl
  GREEN (regression baseline).
- **POST-CARVE (the proof):** G1+G2 (product+LS via `_MESHES`) + G3 (all meshes)
  + G4 + G5 GREEN; then the (a)-(e) mutation matrix RED-verified, (a)/(b) showing
  the product-RED/LS-GREEN asymmetry.
