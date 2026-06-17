---
name: issue-240-d5b-s3-a-consumer-closeout
description: #240/#38/#37 D5b-S3-A — the φ̂-iterate CONSUMER fold (select SpatialMomentSpace + thread Σ_s·φ̂ through forward/SI). STOPPED at a genuine architectural fork the crosswalk does NOT resolve — surfaced + characterized, NO code written. The fork the prior in-flight closeout did not reach because it was BLOCKED on the typed-field space (now resolved by S3-A0): WHERE the slope-scattering source Σ_s·φ̂ enters the d=1 paths. BOTH 1-D production paths (SI scan + Krylov matvec) lack a (2,)-moment SOURCE seam, and `-S` composes element-wise on the Schur-REDUCED scalar matvec output, so the crosswalk's "thread via D1ClosedForm.kernel_rhs/schur_xV s_hat" cannot be wired without a decision on the operator-composition layer. Returned a concise question to the main agent.
metadata:
  type: project
---

# #240 D5b-S3-A — Inc C φ̂-iterate CONSUMER: STOPPED at an architectural fork

**Branch** `feature/sn-space-angle-tier2` (S3-A0 done + committed `d313d16`/`96dfc96`).
**NOT committed; NO code written this session.** Host env, `.venv/bin/python -O`.
Canonical `python -O -m pytest`; NEVER all `tests/sn` (#212).

## STATUS — pre-reads DONE, full seam verification DONE, STOPPED at a fork

I read all 6 mandatory pre-reads in full and verified every seam the crosswalk
lists against the LIVE tree. The typed-field-space BLOCKER from the prior
in-flight closeout (`issue_240_d5b_s3_a_inc_c_closeout.md` §BLOCKER) IS resolved
by S3-A0 (the `spatial_moments` factory param + `_compose_spatial_moments` +
`find_factor` are all live in `_bases.py`/`harmonic_moment_field.py`; verified by
Read). The `Σ_s ⊗ I` scattering einsum lift in `material_xs_field.py` is committed
(verified — NOT to be redone). The cell-emit, both iterate carriers, and the
scattering accumulator widening are all mechanically straightforward GIVEN a
resolved source-entry decision.

But there is a genuine architectural fork the crosswalk's seam map does NOT
resolve, at the EXACT place the prior dispatch never reached (it died on the
field-space BLOCKER one layer earlier). Wiring it by guessing would ship either
a correctness bug behind a green gate (the precise Mode-9-misclassification the
test-architect memo §0/§4 warns about) OR a shape crash. Per the brief's
explicit instruction ("If you hit a genuine architectural fork the crosswalk
does not resolve, STOP and return a concise question rather than guessing"), I
stopped and surfaced it.

## THE FORK — where does the slope-scattering source `Σ_s·φ̂` enter the d=1 paths?

The crosswalk seam #3 (d=1) says: *"`D1ClosedForm.kernel_rhs` hard-codes Q̂=0;
`schur_xV` already has an `s_hat` arg but `_kernel_terms` calls `kernel_rhs`
FLAT. Thread the scattering-slope into the d=1 scan source so the slope row gets
`Σ_s·φ̂`."* Verifying the LIVE routing shows this is under-specified — the d=1
LD slab uses TWO distinct production paths, and NEITHER has a `(2,)`-moment
SOURCE seam, and the operator-composition layer that subtracts `−S` does not
route a slope source into the cell system:

### Fact 1 — the GATE 1 leg 1a tripwire (#37) uses KRYLOV, not SI.
`tests/sn/verification/mms/test_mms_ld_slab.py::test_ld_thick_diffusive_limit_xfail`
(`:257-265`) drives BOTH DD and LD with `inner_solver="krylov"`. So the 1-D leg's
SUT is the d=1 Krylov **matvec**, not the SI scan. (The brief's "1-D #37 via the
existing d=1 krylov-slab matvec" matches this.)

### Fact 2 — the d=1 Krylov matvec is the Schur-REDUCED SCALAR residual; `−S` is element-wise at the OperatorSum level.
- `InvertibleOperator.apply(ψ)` (`operator.py:1354`) = `loss_representation.loss_action(self.sigma, ψ)` = `(L+C)ψ`, the FULL within-group loss. For 1-D it routes `CumprodScan.loss_action` → `_OneDimScanWalk.loss_action` → `_apply_walk` → per cell `scheme.residual_kernel_batch(... Q_cells=_MATVEC_ZERO_SOURCE ...)` (`loss_representation.py:2165-2173`). d=1 `residual_kernel_batch` (`linear_discontinuous.py:656-662`) calls `_kernel_terms` → `D1ClosedForm.kernel_rhs(Q_cells, in0)` → returns a SCALAR per-cell residual `eff_denom·ψ̄ − rhs`. The slope `ψ̂` is Schur-ELIMINATED inside `eff_denom`/`rhs` — there is no slope ROW in the matvec output; the probe `psi_bar` is scalar `(ng, K)`.
- The full Krylov operator `A = (L+C) − S` is an `OperatorSum`: `S.apply(ψ)` is subtracted ELEMENT-WISE on the bulk values (the `OperatorSum.apply` leaf sum; `InvertibleOperator.apply` only OWNS `(L+C)`, see `operator.py:1339`). `S.apply(AngularFlux)` (`scattering.py:1067`) → `AngularSourceSink` per-ordinate.

⇒ If the iterate `ψ` carries φ̂ (the `(N,ng,nx,2)` AngularFlux), then `S.apply(ψ)`
produces a `(N,ng,nx,2)` per-ordinate source (the `Σ_s ⊗ I` lift is moment-
agnostic — verified). But `(L+C)ψ` from the d=1 Schur-reduced matvec is SCALAR
`(N,ng,nx)` — element-wise `(L+C)ψ − S.apply(ψ)` is a SHAPE MISMATCH (or, if
`(L+C)ψ` is also widened to `(N,ng,nx,2)`, then the d=1 matvec must produce a
GENUINE slope-moment ROW, i.e. STOP Schur-reducing — a different code path than
the crosswalk's "fold `s_hat` into `kernel_rhs`").

### Fact 3 — the 1-D SI scan path has NO `(2,)`-moment source seam at all.
The within-group SI driver (`_within_group_si`) on a 1-D mesh runs
`CumprodScan` → `_OneDimScanWalk.sweep` → `_run` (`loss_representation.py:2556`).
The scan source is `b = scheme.source_emission(QV_chain, inverse_denom, w)`
(`:2733`/`:2903`) where `QV_chain` is the per-ordinate SCALAR cell-average source
`Q·V` `(K,ng,nx)`. The scan affine recurrence `ψ_out = a·ψ_in + b` carries ONLY
the cell-average source `b`; the slope `ψ̂` is reconstructed locally from the
`(a,inverse_denom,w)` coefficients (via `scan_xV`, which has NO `s_hat` arg —
unlike `schur_xV`). So the SI scan has NO place to inject `Σ_s·φ̂`. The crosswalk
pointed at `kernel_rhs`/`schur_xV`, but the SI scan uses NEITHER — it uses
`affine_scan_coefficients` → `scan_xV`.

### Fact 4 — the 2-D leg (SI) is cleaner, but shares the same `−S` composition question.
The 2-D thick-box analog (GATE 1 leg 1b, NEW) MUST run via SI (the d≥2 Krylov
matvec raise `_CellResidual.cell:929` is S3-B, out of scope). The 2-D SI path is
the windowed `_MomentWindowedResolvent` → `solve_moments` → `moment_buf`, and
the d≥2 `_ubld_system` `R = M·S_moments` (`linear_discontinuous.py:531-536`)
DOES have a genuine `(2^d,ng)` source slot (currently scalar-lifted to slot 0).
Here the source seam is real: thread `(2^d,ng)` `Q_cells` carrying `Σ_s·φ̂` in the
slope rows. BUT: the SI loop builds `rhs = q_ext + S.apply(ψ) + B.apply(ψ)`
(`iteration.py:503-505`) and passes it to `L.solve(rhs)` (the sweep). So the
slope-scattering source arrives as part of `Q` THROUGH the sweep — which means
`Q_cells` must carry the `(2^d,ng)` moment source, AND the sweep's
`_SolveOperands.Q` carrier must grow the trailing `2^d` axis (the crosswalk's
"gate the widening on LD-d≥2" — DD must stay byte-id). This is mechanically
do-able, but only AFTER deciding whether the d=1 Krylov path is unified with it
or stays separate (Fact 2/3).

## THE QUESTION (returned to the main agent verbatim in the assistant message)

For the d=1 (1-D #37) path, which of the following is the intended architecture
for routing `Σ_s·φ̂` into the LD cell balance?

- **Option A — full `(2,)`-moment d=1 matvec.** Stop Schur-reducing the d=1
  matvec; make `(L+C)ψ` return the genuine `2`-moment residual (average ROW +
  slope ROW), so the iterate carries φ̂ end-to-end and `−S.apply(ψ)` subtracts
  the `(2,)`-moment scattering source row-for-row. Unifies d=1 with d≥2 (both
  carry the full `2^d` moment vector through the matvec). COST: the d=1 matvec
  leaves the L16 closed-form scalar fast path; the Krylov flat vector doubles;
  AND the 1-D SI scan (`_run`/`scan_xV`) must ALSO be widened to a `(2,)`-moment
  scan OR the 1-D SI path must switch off the scan onto the wavefront/kernel for
  LD when φ̂ is active. (The 1-D SI #37 leg is via Krylov so the scan widening is
  not strictly required for GATE 1, but GATE 2's `test_sn_1d_slab_ld_mms_*`
  exercises the SI scan AND the two-paths gate `CumprodScan ≡ FullFieldWavefront`
  must stay green — so the scan path's source convention has to remain consistent
  with whatever the wavefront does.)

- **Option B — keep the Schur-reduced scalar d=1 matvec; fold `Σ_s·φ̂` into the
  cell SOURCE (not the OperatorSum `−S`).** The crosswalk's literal reading: the
  slope-scattering source enters `kernel_rhs`'s `eff_source` via the `s_hat`
  fold (exactly as `schur_xV` already does for the per-cell `update`/`residual`).
  This keeps the d=1 matvec SCALAR and the L16 fast path. BUT it requires the
  slope-scattering source to reach the SWEEP/MATVEC's `Q_cells` as a `(2,)`-moment
  source — which means `−S` must NOT be a post-matvec element-wise subtraction
  for the slope component; the slope row of `S.apply(ψ)` must be routed INTO
  `Q_cells` (the cell-system RHS), while the average row stays the usual `−S`.
  That splits `S.apply(ψ)` across two composition sites (average → OperatorSum
  `−S`; slope → sweep `Q_cells`), which breaks the clean `A = (L+C) − S`
  OperatorSum algebra (Pattern 1 "read as the math" + the InvertibleOperator
  `loss_action(self.sigma)` single-source). For the SI scan path this is even
  more invasive (the scan has no moment-source slot — Fact 3).

- **Option C — something else** (e.g. the φ̂ slope is a PURELY LOCAL per-cell
  reconstruction that never travels between sweeps, and the diffusion limit is
  recovered by the slope-SCATTERING source being assembled inside the cell
  system from the iterate's average moment only — but that contradicts the
  crosswalk's "carry φ̂ between sweeps" and the test-architect's "the FP CHANGES
  because the slope source `Σ_s·φ̂` couples the slope rows").

My read: the crosswalk + the durable plan + the S3-A0 closeout all point at
Option A's SPIRIT (φ̂ is a first-class typed factor that travels between sweeps;
the `SpatialMomentSpace` was minted exactly so the iterate can HOLD it) — but the
crosswalk's seam #3 WORDING ("thread via `kernel_rhs`/`schur_xV` `s_hat`") is
Option B's mechanism, and the two are incompatible at the OperatorSum `−S`
composition layer. The d=1 Krylov matvec being Schur-reduced-scalar is the crux:
it cannot subtract a `(2,)`-moment `S.apply(ψ)` element-wise without either
(A) widening its own output to `(2,)`-moments or (B) splitting `−S` across two
sites. This is the decision I need before writing any code, because it dictates
(i) whether the d=1 matvec stays scalar or becomes `(2,)`-moment, (ii) whether
the Krylov flat vector / `to_flat` grows, (iii) whether `_SolveOperands.Q` /
`_ApplyOperands.probe` widen, (iv) whether the 1-D SI scan path is widened or
bypassed for LD-with-φ̂, and (v) how `OperatorSum.apply` composes `(L+C) − S`
for a moment-carrying LD iterate.

## What I did NOT touch (NO code written)
- No production code, no test code, no Sphinx, no `git add`.
- `material_xs_field.py` Σ_s⊗I lift: confirmed committed (S3-A0) — LEFT.
- The d≥2 Krylov matvec raise (`_CellResidual.cell:929`): S3-B — LEFT.

## Seams VERIFIED live (the crosswalk's map is accurate EXCEPT the d=1 source entry)
- `_CellSolve.cell` (`sweep_graph.py:858-896`): the φ̂ emit gate `len(s_axes) > 1
  and spatial_basis_per_axis > 1` drops to `psi_avg[..., AVERAGE_MOMENT]` (`:884`).
  Mechanically straightforward to widen the moment-mode emit
  `"nlm,ngd,n->lmgd"` → `"nlm,ngdp,n->lmgdp"` once the iterate carries `p`.
- `ScatteringOperator.apply` (`scattering.py:1050/1067/1093/1002`): moment-agnostic
  post-S3-A0 (the `_assemble_per_ordinate_source` accumulator
  `ScalarSourceSink.zeros_on(mesh)` `:514` needs `spatial_moments=` selected — the
  prior closeout's accumulator-broadcast blocker; the S3-A0 param makes this a
  one-line select once the iterate width is known).
- `integrate_angular` (`angular_flux.py:131`) `"n,ng...->g..."` — moment-agnostic.
- `D1ClosedForm.kernel_rhs` (`_ubld.py:324`) hard-codes Q̂=0; `schur_xV` (`:343`)
  has `s_hat`; `scan_xV` (`:383`) has NO `s_hat` (the SI scan gap, Fact 3).
- `_ubld_system` (`linear_discontinuous.py:531-536`) `R = M·S_moments` HAS the
  genuine `(2^d,ng)` slot (d≥2 source seam — the clean one).
- The un-windowed iterate construction: `_within_group_si` (`solver.py:672`),
  cold-start `_windowed_cold_start` (`:607`, windowed only — the un-windowed
  full-angular cold-start is `initial_guess=None` → `_zeros_like(q_ext)` in
  `SourceIteration.solve:481`), `_maybe_window` (`:575`).

## LESSON (candidate for the crosswalk / coding-elegance Pattern 7)
The S3-A0 closeout resolved the typed-field SPACE blocker (where φ̂ LIVES), but a
SECOND, distinct fork sits one layer deeper: where φ̂'s scattering source ENTERS
the matvec/sweep, given that `A = (L+C) − S` is an OperatorSum that subtracts `S`
ELEMENT-WISE on the matvec OUTPUT, while the LD cell balance needs the slope
source folded INTO the cell SYSTEM (the Schur `eff_source`). A carve that adds a
moment axis to an iterate must map BOTH (a) where the axis LIVES (the field space
— S3-A0) AND (b) how the operator ALGEBRA composes the axis-carrying source
across `(L+C) − S` (this fork). The crosswalk mapped (a) and the d≥2 cell-system
source slot, but not the d=1 `−S`-composition layer where the Schur-reduced
scalar matvec and the moment-carrying scattering source meet. Add an
OPERATOR-COMPOSITION row to the Pattern-7 crosswalk whenever a carve adds an axis
to a quantity that flows through an OperatorSum (`A = Σ ops`): "for each op, does
the axis-carrying source compose element-wise on the op's OUTPUT, or must it fold
INTO the op's internal system?"
