---
name: issue-158-matvec-via-kernel-trait
description: #158 Inc B follow-up — matvec_via_kernel capability flag (the matvec analog of is_affine_scannable); DD-byte-id-via-flag PASS
metadata:
  type: project
---

# #158 Increment-B follow-up — `matvec_via_kernel` capability flag (PASS, 2 nits)

Verdict on the focused re-review AFTER my Inc-B coefficient-model PASS-WITH-NITS (all
those NITs fixed). Branch `feature/sn-space-angle-tier2`, working tree (UNCOMMITTED at
review time; #158 lives on `2b56348` + the 4-file working diff
`loss_representation.py`/`cell_update.py`/`linear_discontinuous.py`/`diamond.py`;
`affine_closure.py`+`cell_balance.py` already committed).

**Trigger for the change:** the version I (and qa) reviewed routed the 1-D **Cartesian
matvec** UNIFORMLY through `cell_update.residual_kernel_batch` for BOTH DD and LD. The ÷V
`residual_kernel_batch` re-associates vs the old ×V `cell_balance_for_streaming` → DD
drifts ~1 ULP on **non-power-of-2 cell widths** → 11 STRICT `slab_*_apply_bit_identical`
snapshots broke (`test_streaming_operator.py`). So my earlier "DD byte-id" finding held
ONLY on the gate's power-of-2 meshes — a real hole the uniform routing opened.

**The fix = a capability trait, NOT a re-baseline:**
- `cell_update.py:566` `CellUpdateBase.matvec_via_kernel: ClassVar[bool] = False` (opt-out
  default; the matvec analog of `is_affine_scannable`).
- `linear_discontinuous.py:282` LD overrides `= True` (LD's matvec IS its
  `residual_kernel_batch` — LD has NO separate `cell_balance` density form).
- `loss_representation.py:2020` `_apply_walk` Cartesian branch:
  `if curvature=="cartesian" and sn_mesh.cell_update.matvec_via_kernel:` → LD rides the ÷V
  kernel; **DD (`False`) falls through to its byte-identical `cell_balance_for_streaming`
  density path** (the #206 Phase-C carve; ALSO DD's curvilinear M-M matvec arm).

## RULINGS (all 3 of the user's questions → PASS)

1. **NOT a twin path → legitimate polymorphic dispatch.** The two arms are GENUINELY
   different per-scheme cell math: DD = `(denom·ψ̄ − numer)/V` diamond `out=2ψ̄−in` via
   `cell_balance_for_streaming` (which ALSO carries the curvilinear M-M
   `angular_denom_term`/`angular_numer_upstream` — geometry-SPANNING, can't be a pure
   `(a,inverse_denom,w)` coeff); LD = ÷V Schur `S·ψ̄−rhs`, `out=ψ̄+(g/θ)(ψ̄−in)/d2`. SSOT
   intact PER SCHEME: LD arm shares `_kernel_terms` with `cell_kernel_batch`
   (linear_discontinuous.py:443/468); DD arm shares `cell_balance_for_streaming` with DD
   `residual`. NO new fold either side → NOT the phased-carve twin-DELIVERY smell
   (institutional #1/#2). Selector reads a trait the scheme declares about ITSELF, not a
   caller-injected mode → Pattern 1/5, NOT anti-pattern-3 boolean dispatch.

2. **Capability flag was the RIGHT call vs re-baselining 11 snapshots.** `cell_balance` is
   STRUCTURALLY REQUIRED regardless (it's DD's curvilinear matvec AND the scan-cache denom
   form, #206 Phase-C) → routing DD through it keeps DD on its ONE native bit-id form, NOT
   preserving a superseded path. Re-baselining would (a) lose DD-hot-matvec regression
   resolution, (b) COUPLE DD goldens to LD-convention choices (a future LD edit nudges DD's
   snapshot). Flag DECOUPLES. anti-pattern-17: the `-W error::DriftWarning` gate IS the
   contract; the flag DEFENDS it, doesn't relax it. Default `False`=opt-out mirrors
   `is_affine_scannable=False` + `cell_kernel_batch` raise-by-default → "DD via kernel" is
   UNREPRESENTABLE (DD never declares the trait) = Pattern 4.

3. **Docstrings accurate; stale "kernel for DD and LD alike" claim GONE everywhere.**
   `affine_closure.py:43-56` apply note correct (apply NOT a generic op; LD declares flag;
   DD keeps `cell_balance`); the `.. note::` :58-69 = textbook honest split (DD `w=½`
   byte-id since power-of-2 scaling commutes w/ IEEE; principled-equiv `w≠½`=LD only) AND
   explicitly says DD snapshots do NOT re-baseline (:64-66) → CLOSES my old Inc-B NIT-2.
   `_apply_walk` docstring :1894-1905 + inline :2020-2049 + `cell_update.py:566-581` all
   precise (non-pow2 ~1ULP rationale, LD-never-reaches-`cell_balance`, geometry-blind
   curvilinear caveat).

**PINS verified (branch is NOT a latent path):** LD arm pinned END-TO-END by
`test_mms_ld_slab.py::test_sn_1d_slab_ld_mms_krylov_matches_si` (:120 — Krylov inner drives
LD `residual_kernel_batch` via `_apply_walk` ≡ SI sweep, L14 matvec≡sweep). DD byte-id arm
pinned by `test_streaming_operator.py` `slab_*_apply_bit_identical` (:801-849) under
DriftWarning-error. Both `if`-arms have a live witness. Carry-forward (NOT new): both LD
pins ride Q̂=0 flat slice → slope-source sign-trap (LM-1989 §1.4/§6) un-probed at matvec
until Inc-C (the standing #158 carry; `linear_discontinuous.py:401-406` states it).

## NITS (non-blocking)
- NIT-1 (do-now cosmetic, Smell-#7): `cell_update.py:580-581` "geometry-blind to this flag"
  is loose — the flag is only CONSULTED on the Cartesian branch; curvilinear DD bypasses it
  structurally (the gate lives inside `if curvature=="cartesian"`). Reword to "—that arm is
  on the non-Cartesian branch, which never consults this flag." Won't mislead future
  maintainer into thinking the flag has curvilinear semantics. Zero code impact.
- NIT-2 (record, Inc-C seam): when curvilinear LD publishes (`affine_scan_coefficients`
  slab-only guard linear_discontinuous.py:543), `matvec_via_kernel` becomes
  geometry-entangled (LD-Cart-matvec=÷V kernel, LD-curv-matvec needs M-M `cell_balance` or
  a curv kernel). bool ClassVar fine for slab-only present (LD raises on curv at construction
  → illegal state already unrepresentable). Flag the moment the trait may need to become
  geometry-aware. Belongs to existing #158 curvilinear-LD scope, NO new issue.

APPROVED for commit; no blocking conditions. Idiom catalogue: `matvec_via_kernel` joins
`is_affine_scannable` as the per-leaf MATVEC capability trait — a scheme declares whether
its apply IS its `residual_kernel_batch` (kernel) or a specialised density form
(`cell_balance`). Future Step/LLD declares accordingly w/ zero `_apply_walk` change.
