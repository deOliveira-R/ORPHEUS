---
name: issue-240-phase2-step-d2-homing
description: #240 Phase 2 D2 review — the 3 affine reconstruction ops homed from affine_closure.py free-func module onto CellUpdateBase staticmethods; affine_closure.py git-rm deleted; LD w-fold collapsed. PASS clean.
metadata:
  type: project
---

#240 Phase 2 Step D2 (branch `feature/sn-space-angle-tier2`, UNCOMMITTED at review): the architectural
HOMING step my D1 review predicted. The 3 generic advection–reaction reconstruction ops
(`source_emission`, `cell_average`, `outgoing_face_from_average`) MOVE from the dangling free-function
module `spatial/affine_closure.py` → `@staticmethod`s on `CellUpdateBase` (the DiscretizationScheme Base);
`affine_closure.py` `git rm`-deleted. PLUS the D1-carry-forward LD `w`-fold landed.

**VERDICT: PASS clean** (1 latent NIT, do-not-block). The homing is CORRECT + COMPLETE.

**Why PASS (durable facts I verified — trust on re-review):**
- **`@staticmethod` is RIGHT.** All 3 ops are pure functions of (face fluxes / source, blend weight `w`) —
  zero instance state, zero scheme-specific constant. No `self`/`cls`. `CellUpdateBase` is the right home
  (lit-confirmed generic advection–reaction; diffusion scheme will subclass + inherit) — not a dedicated
  mixin (Pattern-6: only one consumer-family today; a mixin would be premature abstraction-over-nothing).
- **Homing is COMPLETE.** Grep `affine_closure` across `orpheus/`+`tests/` → ZERO live imports. Only survivors:
  (a) ONE test COMMENT in test_cell_kernel_batch.py:280 documenting the homing as HISTORY (correct, not a
  ref); (b) doc prose in discrete_ordinates.rst :1509 + matrix.rst :338 naming `test_affine_closure` (the
  TEST module, which KEEPS its name — repointed internally); (c) `docs/_build/html` STALE artifacts
  (regenerate next Sphinx build, not source). The `.. todo::` D6 stub UPDATED to name D2 + the homing.
- **Three caller idioms all CORRECT for their context:**
  - scheme-internal (diamond.py/lin_disc.py): `self.<op>` — the scheme IS-A CellUpdateBase. Right.
  - loss-rep (loss_representation.py :1455/:2637/:2654/:2807/:2827): `self.mesh.cell_update.<op>` — reads
    "the mesh's scheme's reconstruction". VERIFIED the dominant idiom in that file (`mesh.cell_update`
    appears ~12× : .is_affine_scannable :712, cell_update= :901/:963/:1104/:1162, .residual_kernel_batch
    :2070). NOT a smell — it is the established scheme-instance access path. Right.
  - scan.py (free-function module, no instance): `CellUpdateBase.<op>` (class-level staticmethod). Right —
    scan.py has no `self`; calling on the class is the only honest spelling. (TYPE_CHECKING-free import OK.)
- **Rich docstring PRESERVED (Cardinal Rule 3).** The full `affine_closure` module docstring
  (coefficient-model derivation, convex-blend universality, ×V/÷V convention, byte-id-vs-principled note)
  relocated VERBATIM into the `CellUpdateBase` class docstring §"Reconstruction ops"; the 3 per-staticmethod
  docstrings carry over verbatim. Nothing lost. The `:func:`→`:meth:` role-flip is correct (free func→method).
- **LD `w`-fold = GENUINE single-sourcing, byte-identical.** `_kernel_terms` return changed
  `(eff_denom, rhs, d2, g_over_theta, in0)` → `(eff_denom, rhs, w, in0)`; the line
  `w = 1.0/(1.0 + g_over_theta/d2)` moved VERBATIM from BOTH cell_kernel_batch + residual_kernel_batch
  bodies INTO `_kernel_terms` (computed once, both arms consume). Byte-id: same operands, same op order,
  hoisted not rewritten. Verified `d2`/`g_over_theta` are NOT used downstream in either method after the
  fold (dropping them from the return is safe). This RESOLVES my D1 NIT-1 EXACTLY at the predicted
  collapse destination (`_kernel_terms` returns `w`). The future `w(Σ)` Péclet hazard is now closed.
- **Annotation widen** `w: np.ndarray` → `w: "float | np.ndarray"` (stringized, no `from __future__`
  needed on this file but harmless) is CORRECT: DD passes literal `0.5`, LD/coeff-model pass arrays.
  String-only change, no runtime effect.
- **Test migration FAITHFUL.** test_affine_closure.py keeps its 3 fns + bodies verbatim, re-binds
  `cell_average = CellUpdateBase.cell_average` (import-repoint only). kernel-sha goldens REGENERATED
  (source change `op(...)`→`self.op(...)` → new source/AST hash, byte-identical NUMERICS — the strict
  DriftWarning gate stays green per method-implementer 505/1/4 + full 1083).
- **The −1 total-test drop (5250→5249) is BENIGN + EXPECTED.** `test_layer_imports` auto-discovers modules
  via `_iter_python_modules(ORPHEUS_ROOT)` and parametrizes 1 test/module; deleting affine_closure.py removed
  exactly 1 auto-param (271→270 in that file). NO real verification test lost. (Recurring tell: a module
  `git rm` legitimately drops the layer-import count by 1 — don't flag it as coverage loss.)

**NIT (latent, do-not-block):** the diff comment at lin_disc.py :441-444 (now inside `_kernel_terms`)
keeps the LD reconstruction-form prose `(1+k)ψ̄ − k·ψ_in` as "derived ONCE here so both kernel arms share
it" — accurate. No dangling. Watch on D6: when `w(Σ)` Péclet/κ-blend lands, it MUST land in `_kernel_terms`
(the single source) — the fold is the guard that makes one-not-the-other drift unrepresentable.

Carry-forward: D3 renames `CellUpdateBase`→`DiscretizationScheme` + relocates assembly; the 3 staticmethods
ride along. The diffusion-scheme subclass (the homing's whole point) is the first NEW consumer — re-check the
inheritance reads naturally then. All D-series LD pins STILL ride Q̂=0 (the LM-1989 slope-source sign-trap
remains the recurring untested corner; flag again at the slope-Q̂ admission).
