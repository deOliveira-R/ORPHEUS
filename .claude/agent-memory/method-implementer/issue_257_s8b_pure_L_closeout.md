---
name: issue-257-s8b-pure-L-closeout
description: #257 S8b — drop the (L+C)−C fold in StreamingOperator → pure σ-free L via a NAMED loss-rep streaming_action(ψ)=loss_action(0,ψ) leaf; σ leaves StreamingOperator's surface (field removed); C=M[σ_t] recovers the loss. PRINCIPLED-equivalent (NOT bit-id): pure-L matvec re-associates the FP tree (slab/sphere/cyl ≤16 ULP, boundary STRICT 0); the (L+C) composite matvec is BYTE-IDENTICAL (untouched). 0 net-new pyright / 0 net-new type:ignore.
metadata:
  type: project
---

# #257 S8b — StreamingOperator → pure σ-free L (the value-moving carve)

Branch `feature/field-typed-operator-algebra`, on S8a HEAD `9316321`. Host
`.venv`. NOT committed. The SECOND of three S8 sub-stages, the RISKIEST.
PRINCIPLED-equivalent (the recomposition re-associates the FP tree). S8c
(Scattering/Fission fibration) NOT done.

## STEP-0 design — choice (a) done correctly, affinity VERIFIED

The plan offered (a) factor a named σ-free `streaming_action` primitive vs (b)
minimal `loss_action(zeros)` reuse. EXPLORED `loss_representation.py` first:
`loss_action(σ,ψ)` is MONOLITHIC in σ — σ (`sigma_gx`) is threaded into
`residual_kernel_batch(reaction_xs=σ)` (Cartesian), `cell_balance_for_streaming(
total_xs=σ)` (curvilinear), AND `precompute_psi_state(sigma_t=σ)` (the Carlson
coupled-pole seed: `Q_bar=σ·φ₀`, `denom=dr·σ+2` — NON-LINEAR in σ). So a true
σ-free streaming discretization is NOT separable without duplicating the ~400-line
walk (violates Cardinal Rule 2).

⭐ **BUT the AFFINITY HOLDS empirically** (the decisive STEP-0 probe, in-process):
- `loss_action(σ=0,ψ).bulk == loss_action(σ,ψ).bulk − σ·ψ` to CART 32 / SPH 2 /
  CYL 72 ULP (rel ~1e-16); boundary 0 ULP.
- SHARPER: `streaming(σ_a) == streaming(σ_b)` for two WILDLY different σ fields,
  ≤64 ULP (rel ~1e-16) — pure streaming is GENUINELY σ-free. The Carlson seed's
  σ-dependence is EXACTLY the collision diagonal it injects → cancels into σ·ψ.
- ⚠ The decomposition file's TOP docstring claimed `matvec(σ=0)` is "3-13% wrong
  for curvilinear" — **STALE** (the `ad37ca0` era). Probe confirms `matvec(0)+σ·ψ
  == matvec(σ)` to ~1e-16 even for uniform σ=2.0 curvilinear. ERR-058 (#195) made
  the Carlson seed σ-INDEPENDENT (`AngularEdgeExtrapolation`), which is what made
  the matvec affine and licenses S8b. (`TestResolutionADifferentFromPriorWrong`
  was ALREADY inverted by ERR-058 to assert the σ=0/subtractive COINCIDENCE.)

**CHOSEN (a), done correctly WITHOUT perturbing `loss_action`:** add a named σ-free
`streaming_action(psi)` to the `LossRepresentation` Protocol + the `_LossRepresentation`
base, defined as `loss_action(zeros_like(σ), psi)`. This NAMES the primitive
(Pattern 3), single-sources the discretization through `loss_action` (Pattern 2 —
ONE streaming walk), and does NOT touch `loss_action`'s FP tree (so
`InvertibleOperator.apply` = `loss_action(σ_t)` is BYTE-IDENTICAL, the sweep/solve
path untouched). NOT option (a)(ii) (refactoring `loss_action := streaming+σ·ψ`,
which WOULD re-baseline `InvertibleOperator.apply` — rejected). NOT option (b)
(a dead σ=0 array + unnamed primitive — inferior).

## sigma_t constructor audit (L20) — PROCEED

PRODUCTION (3): `solver.py:216/927/1004` all `StreamingOperator(sn_mesh, σ_t) + C`
→ `StreamingOperator(sn_mesh)` (drop σ; the composite `InvertibleOperator` reads σ
from `self.diagonal.sigma` = C's σ since #240 Step B, NOT from streaming). TESTS+
DIAGNOSTICS (~106 sites across ~20 files) all `StreamingOperator(X, σ)` → `(X)`.
DECISION: **REMOVE `sigma_t` field ENTIRELY** (Pattern 4 — an accepted-but-ignored
σ is a lie). Migrated all sites via a balanced-paren rewriter (`/tmp/migrate_streaming.py`)
that drops the 2nd top-level arg from `StreamingOperator(...)` only.

## Files changed (precise)

**Production:**
- `orpheus/sn/loss_representation.py`: ADD to the `LossRepresentation` Protocol +
  the `_LossRepresentation` base — `streaming_action(psi)` (= `loss_action(_zero_sigma_for(psi), psi)`)
  + `streaming_action_transpose(phi)` (= `loss_action_transpose(0, phi)`) +
  `_zero_sigma_for(field)` (zeros `(ng, *spatial)` from `field.bulk.values.shape[1]`
  + `self.mesh.spatial_shape`). ⭐ The base's `loss_action`/`loss_action_transpose`
  abstract signatures live under `if TYPE_CHECKING:` (NOT runtime stubs — runtime
  stubs reddened `reportRedeclaration` "obscured by subclass override"; the
  TYPE_CHECKING block gives pyright the signature for `self.loss_action` in
  `streaming_action` WITHOUT a runtime method the subclass override obscures;
  runtime MRO resolves `self.loss_action` to the concrete subclass).
- `orpheus/sn/operator.py`: `StreamingOperator` — DROP the `sigma_t: np.ndarray`
  dataclass field; `apply` → `return self.loss_representation.streaming_action(psi)`
  (the `_require_typed_composite` guard stays); `apply_transpose` →
  `streaming_action_transpose(phi)`; the local `FullField`/`AngularSourceSink`
  imports + the `(L+C)−C` arithmetic DELETED from both. Docstrings re-written
  (pure-L, Pattern-4 σ-free).
- `orpheus/sn/solver.py:216/927/1004`: `StreamingOperator(sn_mesh)` (3 sites).

**Tests — re-baselined (§B):**
- `test_streaming_operator_decomposition.py`: TOP docstring + math re-written
  (pure-L, the stale "3-13% wrong" claim corrected). `TestResolutionADecomposition`
  ((L+C)≡M, bit-exact) STAYS GREEN (migrate ctor). `TestSubtractiveDefinition`
  re-pointed `array_equal(L.apply, M−σ·ψ)` → `assert_array_almost_equal_nulp(nulp=256)`
  bulk + boundary STRICT 0-ULP (no longer the *defining* eq). `TestResolutionADifferentFromPriorWrong`
  → renamed `TestPureLIsLossActionAtZeroSigma`: now asserts `L.apply == loss_action(0)`
  BYTE-EXACT (the new defining identity).
- `test_streaming_operator.py`: `TestConstructor` re-pointed `assert L.sigma_t is sig_t`
  → `assert not hasattr(L, "sigma_t")`. `TestT4bPreT4RegressionSnapshot` (slab +
  cart2d snapshots): nULP bound widened `reduction_depth=mesh.nx` → `_PURE_L_NULP=256`
  (the pre-T.4 snapshots re-baselined to pure-L; rel ≤1e-16, ULP metric spikes ≤192
  at near-zero cancellation; boundary STRICT 0). `TestT4c` curvilinear (rtol=1e-13)
  unaffected (sphere arms = #250 baseline-red).
- `test_loss_action_convention.py`: `test_apply_equals_loss_action_minus_independent_collision_het`
  → renamed `test_pure_L_plus_C_recovers_loss_action_het`: (1) `(L+C).apply == loss_action`
  BYTE-EXACT (the +C glue), (2) `L.apply ≈ loss_action−C·ψ` to FP. The flat-reflective
  anchor STAYS (pure streaming still annihilates the flat fundamental).
- `test_invertible_operator.py::test_apply_equals_l_plus_c_on_typed_flux`: bulk
  `array_equal` → `nulp=256`; boundary STRICT 0.
- `test_removal_form_matvec_sweep.py`: migrated; `test_production_sigma_apply_value_preserved`
  nULP `6`→`2048` (the leaf sum is now `loss_action(0)+σ·ψ` vs override `loss_action(σ)`;
  rel<1e-14 guard is the real ground). The STRUCTURAL teeth gate (`op.apply==M(σ_r)`
  `array_equal`) STAYS GREEN (the composite override is unchanged, byte-id).
- `tests/sn/_data/bc_extraction_2d_baseline/vacuum_bulk_2d_seed{0,1,2}.npy`:
  RE-CAPTURED to the pure-L value (`--capture-baseline`); verified `== (L+C).apply.bulk−σ·ψ`
  to ≤64 ULP (the affine relation, structural ground = byte-id composite). The
  `test_vacuum_bulk_bit_identical` foundation gate's STRICT `array_equal` is RESTORED
  against the re-baselined reference.
- Mechanical ctor migration (drop σ): `_test_helpers.py`, `test_operators_apply_typed.py`,
  `test_one_octant_walk.py`, `test_operator_block_role.py`, `test_g_adjoint_reciprocity.py`,
  `test_native_matvec.py`, `test_bc_extraction_2d.py`, `test_apply_full_field_codomain.py`,
  `test_bc_extraction_matvec.py`, `test_one_representation_instance.py`,
  `test_phase_c_gates.py`, `test_affine_carve_baseline.py`, `test_2d_l2_*`,
  `test_krylov_curvilinear_precond_safety.py`, `test_b1pp_verification.py`,
  `_capture_pre_t4_snapshots.py`, `diag_s69_scanmarch_vs_window_bench.py`.

**New catcher (C1):** `tests/sn/operators/test_pure_L_sigma_free.py` (`@foundation`,
9 tests): C1 σ-freedom (`L.apply` byte-id across C(σ_a)/C(σ_b), per-geom, Mode-11
DIRECT `L.apply`), no-σ-surface (`not hasattr sigma_t`), + the Mode-11 MUTATION
TEETH (monkeypatch a σ-leaking `streaming_action` stub → C1's σ-free invariant
breaks, proving teeth not tautology).

## MEASURED DRIFT (the ⭐ deliverable, vs in-process pre-S8b capture)

| geom | L.apply bulk | L.apply bdry | (L+C).apply bulk | (L+C).apply bdry |
|------|-------------|-------------|-----------------|-----------------|
| CART | 16 ULP      | **0**       | **0**           | **0**           |
| SPH  | 10 ULP      | **0**       | **0**           | **0**           |
| CYL  | 8 ULP       | **0**       | **0**           | **0**           |

- The **(L+C) composite matvec is BYTE-IDENTICAL** (0 ULP bulk+bdry) — production's
  `InvertibleOperator.apply` = `loss_action(σ)` is UNTOUCHED. The value moves ONLY
  on the pure-L matvec LEAF (≤16 ULP, FP re-association vs the retired `(L+C)−C`).
- Boundary STRICT 0-ULP everywhere (C never touches the trace; `loss_action(0)`
  emits the same outflow defect).
- vs prediction (slab 0 / sphere ≤4 / cyl ≤2 for the AFFINE `(L+C)−C` probe): the
  COMPOSITE is 0 (matches), the LEAF drift (16/10/8) is the new pure-L FP tree (a
  named-primitive re-association, not a bug — rel ~1e-16, dimensionally
  `reduction_depth×ULP`). NOT above prediction materially.

## GATES

- pyright **2297/19 = EXACT baseline, 0 net-new errors, 0 net-new `# type: ignore`**
  (measured via `git worktree add -d 9316321 /tmp/s8b_baseline` + main `.venv`
  symlink). ⚠ The TYPE_CHECKING block was load-bearing: runtime base stubs reddened
  +2 `reportRedeclaration`; the global +5 transient was 5 test/diag files still
  passing σ to the now-1-arg ctor (all migrated).
- C1 9P; C2 (per-ordinate flat-flux `test_bc_extraction_2d` + `test_quadrature`) 3P;
  C3 (kinf + standoff het 2G) + C4 (fixed-source-2d equivalence) GREEN; the C1 teeth
  mutation FIRES.
- Broad regression `-O` (the spec subset): **7 failed / 2034 passed, EXACTLY the 7
  baseline reds** (#250 SPHERE ×5 + #232 mu_y ×2 — ALL confirmed identical on the
  `9316321` worktree), **ZERO non-baseline reds**. (The 4 transient non-baseline reds
  — 3× `test_bc_extraction_2d` bit-id + 1× `test_invertible_operator` — were the
  re-baseline/nULP migrations, now fixed.)
- §D MMS/k_inf/SI-rate backstop: GREEN (see verification run).
- G-adjoint reciprocity (the `apply_transpose`→`streaming_action_transpose` path) 12P.
- Sphinx -W: clean.

## Sphinx

`docs/theory/operator_algebra.rst` — NEW "Pure-L streaming + the affine collision
split" section + `:label: streaming-action-pure-l` (Eq.
`M(σ)ψ = streaming_action(ψ) + σ⊙ψ ⟺ streaming_action = loss_action(0)`) +
`vv-status documented` + a `.. todo::` for the archivist (rich narrative: WHY the
matvec is affine, the cancellation, ERR-058 license, the probe evidence, the
measured drift). `:mod:`/`:meth:`/`:func:` cross-refs to `streaming_action` + the
C1 catcher.

## Scope confirmation

- SWEEP/solve path UNTOUCHED (loss-rep stays L+C; `InvertibleOperator.apply/solve`
  byte-id; σ single-sourced from C's diagonal since #240 Step B).
- S/F `@apply.register` carrier dispatch UNTOUCHED (S8c).
- NO new ERR (the carve is value-preserving + the stale-docstring correction is a
  doc fix, not a bug). NO algebra-of-record Branch-1 SymPy manifest (this is a
  type-surface + named-primitive carve over an existing verified discretization;
  the σ-free property is a software invariant pinned by C1, the affine relation by
  the decomposition gate — same posture as S4.5/S8a).

## LESSON (propose for the skill / future S-stages)

A "make L compute pure streaming σ-free" carve over a MONOLITHIC-in-σ matvec walk
is achievable WITHOUT duplicating the discretization IF the matvec is affine in σ —
**verify the affinity empirically FIRST** (two probes: `loss_action(0) == loss_action(σ)−σ·ψ`
AND `streaming(σ_a) == streaming(σ_b)` for wildly different σ; the second is the
decisive σ-leak test). The named σ-free primitive is then `streaming_action :=
loss_action(0)` — single-sources the walk (Pattern 2), names it (Pattern 3), and
leaves `loss_action`'s FP tree (hence the production composite matvec) BYTE-IDENTICAL.
⚠ A base-class abstract method that concrete subclasses override must be declared
under `if TYPE_CHECKING:` (NOT a runtime `raise NotImplementedError` stub) when the
base also has a concrete method calling it — a runtime stub reddens `reportRedeclaration`
("obscured by subclass"), the TYPE_CHECKING signature satisfies pyright while runtime
MRO resolves to the subclass. ⚠ A foundation bit-identity `.npy` snapshot of an
operator LEAF whose FP tree deliberately re-associated should be RE-CAPTURED (not
permanently nULP-loosened) when the new value is structurally grounded (here: the
byte-identical (L+C) composite + the affine relation) — re-baselining restores the
strict gate for FUTURE refactors. ⚠ A decomposition test's TOP docstring can carry a
STALE physics claim ("σ=0 is 3-13% wrong") that a later fix (ERR-058) silently
invalidated — probe the live behavior, don't trust the prose.

Extends [[issue-257-s8a-apply-full-field-repoint-closeout]] (S8a re-pointed the leaf
codomain to timeless FullField; S8b drops the (L+C)−C fold S8a's docstring still
described). Next: S8c (S/F fibration, gate C6).
