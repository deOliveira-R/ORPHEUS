# B — Operator machinery (Part III-A) + Phase 2 (splitting) + Phase 3 (pencil): reconciliation vs HEAD

Reconciled 2026-08-19 against HEAD `7aae9bf1` (main; `orpheus/` + `tests/` clean, zero
uncommitted edits in the audited tree). Method: direct reads + greps + one measured
suite run; Nexus not needed for this slice (all questions were symbol/text-shaped).
Every `[M]` below was produced this session; commands stated inline or implied by the
file:line. Line numbers are at HEAD `7aae9bf1` — re-derive after any merge.

## 0. Headlines — the three findings that reshape the campaign plan

**H1 — The plan has an older, TRACKED, user-ruled sibling covering the same ground, and
never cites it.** `[M]` the report v2 is **untracked** (`git ls-files --error-unmatch`
fails); `.claude/plans/operator_strategy_realization_campaign.md` (2026-07-28, tracked)
is the "Operator · Splitting · Realization" campaign whose P3/P4/P5/P6 phases occupy
exactly the report's Phase 2 + Phase 3 territory, with **user rulings** (2026-07-28/29)
the report does not know. `[M]` `grep -in "strategy|realization campaign"` over the
report: zero cross-references. Campaign state: **P0 COMPLETE + MERGED; P1 NOT started**
— `[M]` today: `python -O -m pytest tests/sn/architecture/ -q` → **105 passed,
21 xfailed, 1.66 s**, bit-matching the campaign plan's own 2026-08-13 checkpoint.

**H2 — Phase 3's proposed mechanism is contradicted by a standing user constraint.**
Report 3.1: "Split `KEigenvalue` → `CriticalityProblem` + `PowerIteration`". Campaign
P4 (the tracked plan, §"⬢ P4"): `GeneralizedEigenPencil(lhs, rhs, spectral_parameter)`
is an **ADDITIONAL Layer-2a object** and "**must NOT displace `power_iteration`'s late
binding**, which is strictly more general (it admits CP BiCGSTAB, diffusion's direct
inverse, the homogeneous dense solve — none of which has an `(L,S,F)` factorization)"
(constraint L23a). The report's *goal* (pencil as object; α/shift as posing rows)
survives; its *means* (a `PowerIteration` solver class replacing the loop) is refuted
by the campaign's own constraint. `[M]` neither `Pencil` nor `GeneralizedEigenPencil`
nor `CriticalityProblem` exists anywhere (`grep -rn` over `orpheus/ tests/`: 0 hits).

**H3 — The #340 convergence campaign (2026-08-09..13) rebuilt the exact seam Phases
2/3/C would wire into, after the plan was written.** `EigenvalueSolver.converged() →
bool` was **replaced** by `measure_stopping_criteria(...) -> tuple[StoppingCriterion,...]`
(eigenvalue.py:106, ⛔ note dated 2026-08-09, #340 N2b); `power_iteration` now returns
`PowerIterationOutcome` (frozen dataclass, eigenvalue.py:304) carrying an
`IterationRecord` tree (per-criterion trajectories + per-inner children + typed
`IterationBudget`). Any Phase-2.3/2.4/3.1/C design that says "add a reading to the
loop" now has a named, solver-agnostic channel — and any design drafted against the
old `(keff, history, flux)` triple is stale.

---

## 1. Part III-A — claim-by-claim verdicts

| # | Plan claim (III-A) | Verdict | `[M]` evidence |
|---|---|---|---|
| A1 | Four-member inverse family + `InverseWrapMixin` | **CONFIRMED** | `InverseWrapMixin` operator.py:2123 (abstract `apply(x, /, *, initial_guess=None)`; solve=forward; involution `inverse()→inner` by object identity). Members: `InverseOperator` operator.py:2219; `SweepOperator` sn/operators/sweep_operator.py:83 (`InverseWrapMixin["SweepInvertible"]`, union at :80 = `StreamingCollisionOperator \| ScheduledInvertibleOperator`); `GreenOperator` numerics/green_operator.py:192 (`InverseWrapMixin[OperatorSum]`); `MatrixInverseOperator` numerics/matrix_inverse_operator.py:95 (eager LU; `MatrixTooLarge` gate 4096; "is_invertible deliberately not consulted — reads values, not structure"). No renames since the plan. |
| A2 | Structure-keyed `.inverse()` | **CONFIRMED** | Leaves→`InverseOperator` (also `OperatorProduct.inverse()→InverseOperator(self)`, #285); `OperatorSum.inverse()`→`GreenOperator` (operator.py:1522); `StreamingCollisionOperator.inverse()`→`SweepOperator` (MRO shadow; streaming.py class :445 declares it "the SOLE invertible OperatorSum"); diffusion/homogeneous construct `MatrixInverseOperator` explicitly where the Green diverges (diffusion/operators.py:117 comment). |
| A3 | Seeded-apply as structure (#285 closed) | **CONFIRMED** | Abstract in the mixin (operator.py:~2190); `SupportsSeededApply` Protocol iteration.py:187; `seeded_inverse()` iteration.py:260. |
| A4 | `B.split(schedule)` → `(B_lower, B_upper)` | **CONFIRMED, sharpened** | sn/operators/boundary.py:819 `def split(self, schedule) -> "BoundarySplit"`; `BoundarySplit(NamedTuple)` :1257 with `.lower`/`.upper` = `SNMaskedBoundaryOperator` ("named pair so the two construction sites cannot be swapped silently"). Landed `cc293ef3` 2026-07-01 (pre-plan). Docstring now carries the **post-plan #341 sharpening**: "a *splitting*, NOT a 'regular splitting' — see #341 and :ref:`sn-boundary-gs-not-regular`". |
| A5 | `ScheduledInvertibleOperator` reified, machine-precision round-trip | **CONFIRMED — but slated to DISAPPEAR** | sn/operators/scheduled_invertible.py:96 — `OperatorSum[FullField, FullField, StreamingCollisionOperator, ScaledOperator[..., SNMaskedBoundaryOperator]]` = `M=(L+C−B_lower)`, schedule-triangular solve. ⚠ **User ruling 2026-07-28** (campaign §P5): "expected to DISAPPEAR with P5 — retire, don't rename" once the fold becomes partition+assignment. Do not design new machinery onto this type. |
| A6 | "N is an unreified `*gains` bag — no object owns the (M, N) pair" | **CHANGED (imprecise even at write time)** | Two-level truth. (a) Driver level: CONFIRMED — `SourceIteration.__init__(A_inv, *gains, max_iter=1000, tol=1e-8, corrector=None, budget_name="max_iter")` iteration.py:611 — still variadic, no Splitting type, no ρ. (b) Record level: `WithinGroupSystem` (sn/coupled_system.py:326) **already owns the triple** as fields `loss: CoupledOperator`, `space: CoupledSpace`, `implicit_operator: CoupledOperator\|StreamingCollisionOperator`, `explicit_gains: tuple[LinearOperator,...]` — and it predates the plan (campaign P0 marker 2026-07-28 names it "SAFE TO TOUCH"). The real defect is **R7**: `_select_si_splitting` (sn/solver.py:987) **re-derives** the G-S pair instead of consuming the record — "the splitting is a TWIN" (test_stage_separation.py:83, strict-xfail pinned). So the honest statement: *the pair has a passive record owner; no object owns it as a SPLITTING (with schedule/ρ/contract), and the G-S arm bypasses the record.* |
| A7 | Frame hierarchy `FrameBase → PetrovGalerkinFrame → GalerkinFrame` (#268); `S = RΛM` documented | **CONFIRMED** | numerics/frame.py:114/:334/:374 (+ `_FrameAnalysis`:412, `_FrameReconstruction`:446). |
| A8 | "missing: the kernel/operator split around it and the tightness gate" | **CONFIRMED still missing** | `[M]` `grep -rn "tightness\|frame_bounds\|tight_frame\|is_tight"` over frame.py + scattering.py + tests/numerics/test_frame*: **0 hits**. No `.on(V)` binding exists (`grep "def on("`: only `GeneratingMeasure.on`, unrelated). Apply-time dispatch on S persists — #331 measures `ScatteringOperator.apply` raising "unsupported input type … expected TimedFullField, ScalarFlux, AngularFlux, or HarmonicMomentFlux". |
| A9 | `sweep_acyclicity` with SCC computation (1-D) | **CONFIRMED — but it is a DERIVATIONS oracle, not production** | `orpheus/derivations/discrete/sn/sweep_acyclicity.py` (landed `d9b092c2` 2026-07-28): `TraceSlot`:119, `TraceDigraph`:135 (`strongly_connected_components()`, `is_acyclic`, `cyclic_components()`), `build_slab_trace_digraph`:285 (all 5 law kinds incl. periodic cross-face), `derive_slab_trace_acyclicity`:374. Gate: tests/sn/sweep/test_sweep_acyclicity.py. ⚠ #324: "production never consumes it" — the schedule uses the `permutes_ordinates` heuristic (over-lags acyclic configs; cannot see genuine cycles). No closing-edge RANK anywhere (`grep rank` in the module: 0). |
| A10 | Default splitting recipe implicit (left-spine head = M) | **CONFIRMED** | `OperatorSum.is_invertible` operator.py:1495 (`return self.a.is_invertible`; docstring names "the CANONICAL-ORDERING contract"); `GreenOperator` derives the split via `_left_spine_terms` green_operator.py:154 + `_negated`:178. ⚠ #374 (OPEN) measures the cost: `(I + (−1)·I).is_invertible == True` on the exact zero map — the predicate answers "admits a Green preconditioner", not "is invertible". |
| A11 | `KEigenvalue` fuses Problem and Solver | **CONFIRMED, with a correction to the framing** | iteration.py:1192. It does NOT own a parallel loop (docstring: "delegates the outer loop to power_iteration — the SINGLE power-iteration loop"; it is ONE implementer of `EigenvalueSolver` beside the 5 families). The fusion is *posing + boundary-realization in one class*: `__init__(A, S, F, *, max_outer, keff_tol, flux_tol, max_inner=None, inner_tol, eigenvalue_method)` builds `A.inverse()` (posing) AND realizes `compute_fission_source`/`solve_fixed_source`/`compute_keff`/`measure_stopping_criteria` (solver boundary). Post-plan additions (#340): `max_inner=None` derives the budget from `inner_tol` via `default_iteration_budget` (convergence.py:326); `inner_records` accumulation (⛔ note: "until 2026-08-09 this discarded the inner trajectory INSIDE a shared numerics primitive"). |
| A12 | `EigenvalueSolver` Protocol speaks fission vocabulary | **CONFIRMED** | eigenvalue.py:106 — `compute_fission_source` ("Q_f = χ·(νΣ_f·φ)/k_eff"), `compute_keff`. **CHANGED around it**: `converged()` replaced by `measure_stopping_criteria` (2026-08-09); new refinement Protocols `RecordingSolver`:222 (`inner_records`; all 5 shipped families declare it since #340 N4) and `ProductionRateSolver`:270 (ERR-052 normalisation). |
| A13 | LossRepresentation: "deleted; contents split by answer/cost rule" | **NOT LANDED as stated — evolved differently** | `LossRepresentation` is now a **Protocol** at sn/loss_representation/__init__.py:251 inside a loss_representation PACKAGE; `StreamingCollisionOperator.loss_representation` (streaming.py:~640) delegates to the streaming leaf's cached instance (S6.5/#222 single-source). The σ-free carve happened instead (#257 S8b: `L.apply` IS `loss_action(0, ψ)`). The deletion never occurred; the answer/cost separation partially realized by other means. ERR-026 family retired via the curvilinear campaign (separate route). |
| A14 | ERR-039: "land the in-flight fix with the tightness gate" | **STALE at write time** | Catalog (error_catalog.rst:3107): ERR-039 "**CAUGHT IN QA REVIEW 2026-05-10, fixed in same commit as introduction**" — 3 months before the plan. No tightness gate exists (A8). The plan's "in-flight fix" has no referent at HEAD. |
| A15 | `domain=None` escape hatch verified live; closes as quotient instance | **CONFIRMED still live** | `LinearOperator.domain/codomain` return `Optional[FunctionSpace]` default None (operator.py:~662-687, "pre-date Issue 9.6… composability check is skipped"); R1 xfail reason pins `homogeneous/solver.py:193` `F.domain is None` (measured). The quotient-instance closure is unbuilt. |
| A16 | Guide drift: `BlockOperator` (v3) vs `CoupledOperator` (code) | **CONFIRMED unchanged** | `CoupledOperator` remains the code name (numerics/coupled_system.py). `TrackSweepOperator`: 0 hits (correction #10 holds). |

## 2. Phase 2 rows — state at HEAD

| Row | Verdict | Evidence / who owns it now |
|---|---|---|
| 2.1 Reify `Splitting(M, *gains)` | **NOT LANDED**; superseded-in-shape by campaign **P3** | No `Splitting` type (`grep -rni "splitting"` production: prose only). Campaign P3 ruling (user, 2026-07-29): shape symmetry — loss/M/N take the SAME partition; seedless arm's `explicit_gains` becomes a `CoupledOperator` gain grid; "M⁻¹N becomes formable ⟹ ρ is probeable as an operator expression". P3's measured table: carrying arm already complies (`gains=(N,)`, one 2×2 grid); seedless arm is the anomaly (bare `(S, B)` pair). |
| 2.2 `SplittingRecipe` | **NOT LANDED**; campaign **P5** exit = "a strategy is a value; `_select_si_splitting` and the str flags retire" | Same goal, different vocabulary. R7 strict-xfail (test_stage_separation.py:83) is the landing gate already in the tree. |
| 2.3 `Splitting.spectral_radius()` | **NOT LANDED** — but a *different* ρ landed 2026-08-14 | `5def63b0`: per-closure **face-transmission damping spectrum** with `spectral_radius: float\|None` field (transport/spatial/scheme.py:311, `=max(radii)` :416), consumed by `loss_kernel_gauge.py:398-409` refusal messages. This is the #341/#343 undamped-sawtooth machinery (B's octant action has d−1 eigenvalues exactly −1), NOT ρ(M⁻¹N). |
| 2.4 A priori ρ_SI (3-D form, correction #14) | **NOT LANDED** | No `arctan` ρ formula in production (grep). #14 still unverified. |
| 2.5 SCC closing-edge rank in `sweep_acyclicity` | **NOT LANDED**; issue-tracked | SCC exists derivations-side (A9); rank nowhere. #324 (wire SCC into the B_lower/B_upper decision) + #320 (wire-or-retire the named invariants) are the fuller statements; campaign P5.3 says "Schedule = SCC of the implicit set; generalize sweep_acyclicity.py off the trace digraph (#320)". |
| 2.6 `MonodromyOperator` (promote `TrackMonodromy`) | **NOT LANDED**; ⚠ premise misnamed | `TrackMonodromy` **never existed in code** — `[M]` `git log -S TrackMonodromy` → only `f20b6c98` (the grand-report v3 doc); it is an MoC sheaf-layer CONCEPT there (grand report :2119), not a boundary-closure object. The real item is **#300** (OPEN, blocked on #299): close the boundary SCC via Woodbury on the rank-1 B — with the measured CP evidence that `cp/solver.py:389-397`'s `P_inf` construction IS Sherman-Morrison longhand. #300 is strictly the fuller statement of 2.6. |
| 2.7 Mesh diagnostics (τ>2 warns) | **NOT AUDITED HERE** (out of slice; no `max.*tau` machinery seen in numerics core) | — |
| 2.8 Quarter-symmetry acyclicity gate | **NOT LANDED** | #324's cost-1 ("it over-lags": acyclic configs like `reflective\|vacuum` iterate today) is the live statement; nothing takes a one-pass path. |
| 2.9 CP conditioning fingerprint | **NOT AUDITED HERE** (CP out of slice) | — |
| 2.10 `partition_along` + block extraction `A_ij = π_i A ι_j` | **NOT LANDED**; campaign **P5.1/P5.2** is the ruled design | P5: `Partition` with reconstruction identity `Σ J_i R_i = I` as its gate; "a **second constructor** on the existing block interface: derive `A_ij = R_i A J_j` from `(A, partition)` — *not* a parallel type; one interface, two constructors. (This resolves the fork between the predecessor plan and #296.)" Same formula as 2.10; the campaign adds the anti-fork ruling. |
| 2.11 Reify `BlockDigraph` | **NOT LANDED** | No such type (grep 0). The §I.10 middle layer remains unowned; campaign P5.3's "Schedule = SCC of the implicit set" is the nearest ruled step but does not name a BlockDigraph object. |

**Post-plan Phase-2-relevant facts the plan lacks** (all `[M]`): #341 (2026-08-09,
doc label `sn-boundary-gs-not-regular` in docs/theory/methods/sn/) — the boundary G-S
is NOT a regular splitting; Varga's ρ_GS ≤ ρ_J is **not guaranteed** (d=2/d=3 sign
inversion measured) ⟹ any 2.3/2.4 machinery must not assume regular-splitting
comparison theorems. #343 (2026-08-09) — the octant ORDER is an unowned rate lever:
n_GS/n_Jacobi spans [0.771, 1.892] over the 25 achievable d=3 fold patterns; shipped
lexicographic order is 24th of 25. #373 (2026-08-14) — `M_GS⁻¹` is a **subspace
inverse**: full-space `‖M M⁻¹ − I‖ = 3.25e-01`, defect exactly the (inflow-rows ×
outflow-cols) block; exact on `range(q+Nψ)`; a #200-style Krylov preconditioner built
from it is outside its contract. Any reified `Splitting` type must state this
subspace contract or close it.

## 3. Phase 3 rows — state at HEAD

| Row | Verdict | Evidence |
|---|---|---|
| 3.1 Split `KEigenvalue` → `CriticalityProblem` + `PowerIteration` | **NOT LANDED; means REFUTED by campaign P4/L23a** (see H2) | No such classes (grep 0). What DID land at the "neutral solver boundary": `measure_stopping_criteria` + `PowerIterationOutcome` + `RecordingSolver` (#340) — the reading-based seam CW brackets would use (report 3.1's v2 note) now exists in exactly the right shape, under different names. Vocabulary still fission (A12). |
| 3.2 `Pencil(A, Π).at(σ)` | **NOT LANDED**; campaign P4.1 names it `GeneralizedEigenPencil(lhs, rhs, spectral_parameter)` with contract `domain(A)==domain(M)` "**now non-vacuous, because P1 made domains non-Optional**" — i.e. **P4 is gated on P1**, which has not started. `[M]` 0 hits. Correction #11's twin (`LaplaceResolvent`/`AlphaEigenproblem`): still never minted (0 hits) — the anti-twin ruling holds vacuously. |
| 3.3 α-aware guard branch in `StreamingCollisionOperator.__init__` | **premise CONFIRMED; branch NOT LANDED** | The strict guard exists: streaming.py:597 `if not np.all(diagonal.coefficient.values > 0): raise ValueError(...)` ("min(sigma)=…; if sigma_r … dipping negative, the multi-group set is physically inconsistent"). No α-aware text. ⚠ New context the plan lacks: `StreamingCollisionOperator` is now an `OperatorSum` subclass (streaming.py:445, `L + C` dispatch one-directional per #261) with a mesh-identity invariant (`streaming.sn_mesh is diagonal.coefficient.mesh`) — a shift that rebinds the diagonal must preserve BOTH guards. |
| 3.4 α shift-invert | **NOT LANDED** | Consistent with my prior O-4 recon (pose_alpha / TimeMassOperator: 0 hits). |
| 3.5 Well-formedness assertion | **NOT LANDED** | No cycle-completeness assertion; #320's `assert_cycles_are_declared` remains named-and-unbuilt (sweep_graph.py:49 docstring only). |

## 4. Part I tree-facing premises (§I.1, §I.4, §I.6, §I.10)

- §I.1 boundary-trace row (M=(L+C−B_lower), N=B_upper): **LANDED and wired** —
  `SweepSchedule.gauss_seidel` → `B.split` → `ScheduledInvertibleOperator` →
  `_select_si_splitting` (solver.py:987); the split docstring states exactly the plan's
  algebra ("M = (L+C) − parts.lower and gains = (S, parts.upper); Jacobi ⟹ B_lower=0").
  The plan knew this ("half of the splitting stage exists") — CONFIRMED.
- §I.4 "`StreamingCollisionOperator.__init__` validates σ>0 strictly": **CONFIRMED**
  (streaming.py:597; see 3.3).
- §I.6 rank table: the *analysis* stands; the *tree* implements only the "split"
  column (boundary G-S). Nothing computes rank; nothing closes low-rank SCCs. Rows map
  to issues: white/albedo/periodic close → #300; acyclic-but-iterated rows → #324;
  reflective|reflective split → landed + refined by #341/#343/#373.
- §I.6's "`sweep_acyclicity` already computes the components; the rank computation is
  one addition": **MISLEADING as written** — the components are computed in a
  *derivations* module no production code consumes (#324). "One addition" understates:
  the wiring decision (#320's "production type vs derivations oracle") is open.
- §I.10 three layers: unchanged as analysis; **no layer is reified** (no Partition
  type, no BlockDigraph, schedule = `SweepSchedule` enum-like with the #343 unowned
  order). Campaign P5 owns the reification path.

## 5. Corrections ledger spot-checks (#2, #6–#12, #14, #15/#17 where tree-facing)

| # | Verdict at HEAD |
|---|---|
| 2 (dense-triangular exists) | Analysis-only; no tree contradiction. Consistent with time-implicit design notes. |
| 6 (CausalFactor rename retracted; object deleted) | **Deletion never happened** — see A13. The retraction's premise ("being deleted") is itself stale. |
| 7 (pipeline order) | Analysis-only; campaign pipeline (§1's stage list) agrees posing precedes splitting-instantiation. |
| 8 (shift rides as gain) | Analysis-only; no shift machinery exists to check. |
| 9 (reserve "resolvent" for the λ-family) | **PARTIALLY APPLIED**: the module layer is named `green_operator` (right), but iteration.py prose still calls A "the INVERTIBLE resolvent operand" (:15,:29,:34) and `KEigenvalue.solve_fixed_source` docstring says "Resolvent ``A_loss⁻¹ q``" (:1412). The misnomer survives in 6+ docstring sites. |
| 10 (`TrackSweepOperator` rejected) | **HOLDS** — 0 hits. |
| 11 (Pencil unification) | Twin never minted; the unification target also not built (0 hits both). |
| 12 (M symbol) | Prose discipline only; not checkable. Note the tree adds a 5th M: `MultiplicationOperator` `M[σ]` (streaming docstrings) — the report's reservation list is already short one. |
| 14 (ρ_SI 3-D unverified) | **CONFIRMED still unverified**; blocks 2.4 as stated. |
| 15 (torsor overturned; "prose must be swept") | **NOT APPLIED — and the tree ENFORCES the overturned ontology.** `[M]` `orpheus/numerics/coupled_system.py:108,239,292` still state "the affine flux torsor…"; `orpheus/transport/fields/_flux_role.py` implements the torsor ALGEBRA as arithmetic law (`flux ⊕ displacement → flux; flux ⊕ flux → TypeError`, :97). No `cone.rst`, no cone predicate, no battery (D7 unstarted; `grep Collatz\|Birkhoff` in docs/theory/foundations: 0). The cone ruling is paper-only. **Deep dive = sibling D's slice**; flagged here because #331 (displacement acceptance) sits directly on Phase-2's Krylov path. |

## 6. What the plan does not know (in scope for Phases 2/3/C)

1. **The #340 convergence-contract rebuild** (2026-08-09..13, `d9b027d7`..`0f5ca91c`):
   `measure_stopping_criteria` replaces `converged` (eigenvalue.py:106, ⛔ 2026-08-09);
   `MINIMUM_OUTER_ITERATIONS` single-sourced (was transcribed in all 5 realizers);
   `PowerIterationOutcome` (eigenvalue.py:304) replaces the lossy triple —
   deliberately NOT tuple-unpackable; `converged` DERIVED from the record, never
   stored (#342); `IterationRecord` (convergence.py:910) with `converged`:1066,
   `exhausted_budget`:1109, `fully_converged`:1141 (tree verdict) + `children` per
   inner solve; `RecordingSolver`/`ProductionRateSolver` refinements; budgets derived
   from tolerance (`default_iteration_budget`, convergence.py:326, `ffcce8df`);
   SourceIteration's stop is now the **ρ-honest equation residual**
   `‖rhs_{n−1} − rhs_n‖/‖q_ext‖` (supersedes the iterate-increment norm; iteration.py
   :~500-530) with an end-of-solve certificate (`_certify_within_group_exit`).
2. **ERR-079 → `IterationBudget`** (`0f5ca91c`, 2026-08-13): typed triple
   `limit`/`name`/`iterations_per_unit` (convergence.py:827-829; class :691) — the
   Krylov restart×Arnoldi exchange (`callbacks == maxiter*restart`, measured 12/12
   rows) is now a named conversion. **Residue for Phase 2/3**: any new
   Splitting/Pencil driver must state its budget's unit exchange through this type;
   `budget_name` threads through `SourceIteration`/`power_iteration` signatures.
   That is the whole residue — the defect class itself is closed.
3. **The DSA `corrector` seam** on `SourceIteration` (iteration.py:611,
   `corrector: LinearOperator | None`) — the accelerated two-step with FP-invariance
   contract. A reified `Splitting` must decide whether the corrector is part of the
   recipe or the driver.
4. **#341/#343/#373 splitting theory + measurements** (§2 above) — not-regular
   splitting, unowned octant order (2.5× lever), subspace-inverse contract.
5. **`5def63b0`** (2026-08-14): closure-damping spectrum with per-face
   `spectral_radius` (transport/spatial/scheme.py:311) + `loss_kernel_gauge` refusals.
6. **The assembly axis** (pre-plan `83a0db7b` 2026-07-04, absent from III-A's
   "already built" list): `SupportsAssembly`/`assemblable()`/`MissingAssembly`
   (operator.py:1156/:1212/:508), `OperatorSum.assemble` = `[A]+[B]` sparse sum,
   `SparseAssembledOperator` (assembled_operator.py:77). A third structural axis
   beside inverse/adjoint that Phase-2's "materialize the monodromy" and
   `MatrixInverseOperator` route through.
7. **#280 adjoint-inverse swap law**: `_AdjointOperator.is_invertible`/`inverse()`
   implement `(A*)⁻¹=(A⁻¹)*` by object identity (operator.py:~1330-1370);
   `SweepOperator` carries `apply_transpose`/`is_adjointable` (reverse-scan, #280
   2.5c). The III-A audit predates none of this but does not mention it.
8. **`IntegralKernelOperator` became a category Protocol** (transport/operators/
   integral_kernel_operator.py:164, `Protocol[V]`, the §5.6 locality criterion) — the
   1.5 data/arrow split remains unbuilt, but the object 1.5 proposes to split is no
   longer a concrete operator claiming kernelhood; re-scope 1.5 before executing.

## 7. Issue-overlap map (issue ↔ plan item; which is fuller)

| Issue | Plan item | Fuller statement |
|---|---|---|
| #330 (spaces MANDATORY; space owns shape+traversal; OPEN, 3 comments) | Phase-1 gating + §I.8 identity doctrine; the campaign-P1 issue-home | **Issue** — it carries the user ruling (2026-08-04), the measured Defect-1 (None ⟹ `.H` silently Euclidean), and three tree-wide remainder items (TraceRestrictionOperator lengths blocked on mandatory binding; per-face iteration primitive deferred at 1 consumer; space-declared ordinate axis). The plan's §I.8 is the fuller THEORY; #330 the fuller WORK statement. |
| #331 (leaves disagree on displacement space; OPEN) | §I.8 "error in the domain W"; touches 2.x Krylov | **Issue** — measured per-leaf table (L OK; S, B TypeError); the plan never states the displacement question. Needs a user decision (does domain include V?). |
| #374 (`OperatorSum.is_invertible` wrong question; OPEN) | A10 / 2.1-2.2 | **Issue** — measured zero-map witness + the preferred fix (split the name). A reified Splitting would make the predicate's real question ("admits Green preconditioner") a property of the recipe, dissolving the misname. |
| #375 (`A.H` dead end; OPEN) | Adjoint machinery (III-A periphery) | **Issue** — verified at HEAD: raise stub operator.py:~1316-1321, `is_adjointable` honestly False (not a mis-claim); the 4-line composition `(A*)ᵀ = G_C·A·G_D⁻¹` is the fix. Orthogonal to Phases 2/3. |
| #296 (reify block-operator; OPEN) | 2.10/2.11 | **Split**: #296 is fuller on consumption modes (assemble/apply/solve) + the strategy table (sweep=block-triangular, SI=block-Jacobi/angle, multigroup=block-GS/energy); the plan §I.10 is fuller on the three-layer ownership; **campaign P5.2 rules the fork** (one interface, two constructors — cites #296 verbatim). |
| #300 (close boundary SCC closed-form, Woodbury; OPEN, blocked on #299) | 2.6 | **Issue** — has the CP Sherman-Morrison-longhand evidence (cp/solver.py:389-397) and the rank table. The plan's `TrackMonodromy` premise is misnamed (§2 row 2.6). |
| #320 (creates_sweep_cycle write-only; OPEN) | 2.5 periphery / 3.5 | **Issue** — the flag was RETIRED (`ddc1ee10` 2026-07-30); `assert_cycles_are_declared` named-and-unbuilt; definition-of-done includes the production-type-vs-oracle decision. |
| #324 (wire SCC into the splitting; OPEN) | 2.5 | **Issue** — names the production heuristic (`permutes_ordinates`), its two costs (over-lag; cycle-blind), and the wiring proposal. |
| #343 (octant order unowned rate lever; OPEN) | §I.10 schedule layer / 2.2 | **Issue** — measured 25-pattern spectrum + the separating predicate `max_a L_a > Σ_{b≠a} L_b`; the plan does not know the order matters 2.5×. |
| #373 (M_GS⁻¹ subspace inverse; OPEN) | 2.1/2.2 contract | **Issue** — block-resolved defect table; plan unaware. |
| #366 (two per-level `converged` surfaces; OPEN) | Phase C periphery | **Issue** — post-plan (#340-era); names `PowerIterationOutcome.converged` (zero production readers) and `IterationHistory.converged`. Relevant to where CW brackets report. |

## 8. Refuted / stale premises — one-line structural reasons

1. **"Split KEigenvalue → CriticalityProblem + PowerIteration"** (3.1): refuted as a
   mechanism — `power_iteration`'s late binding is strictly more general than any
   `(L,S,F)` factorization, so the pencil must be additive, not a replacement
   (campaign P4 constraint L23a; user-ruled).
2. **"No object owns the (M, N) pair"** (III-A): the pair had a record owner
   (`WithinGroupSystem`) before the plan was written; the true defect is R7 — the
   driver re-derives instead of consuming — which is a different repair (consume the
   record / make the partition a value) than "give N an owner".
3. **"Promote `TrackMonodromy`"** (2.6): nothing to promote — the name only ever
   existed in the grand-report v3 as an MoC sheaf concept; the boundary-closure work
   item is #300 (blocked on #299).
4. **"ERR-039 in-flight fix"** (III-A): fixed 2026-05-10 per the catalog — no
   in-flight referent existed at plan-writing time; only the tightness GATE remains
   real (and remains absent).
5. **"LossRepresentation: deleted"** (III-A slate): not deleted — reshaped into a
   Protocol + package with σ-free streaming; re-derive the answer/cost residue from
   the current shape before re-planning D5/1.3.
6. **"sweep_acyclicity … the rank computation is one addition"** (§I.6): the SCC
   object is a derivations oracle with zero production consumers (#324), so the
   addition is preceded by an unmade wiring decision (#320).
7. **Correction #15's sweep instruction**: unexecuted, and inverted in force — the
   torsor is now load-bearing ARITHMETIC (`_flux_role.py`), so "sweep the prose"
   understates the change to a type-system carve (sibling D's territory).
8. **The 21-xfail todo list** (P1): `[M]` today — R1: `test_leaf_space_annotation_is_
   not_optional[S,L,F,C,B]` (5) + `test_model_generic_leaf_declares_a_space[{F,C}×
   {2g,4g}]` (4); R2: `test_leaf_without_a_space_refuses_construction[C,S,F]` (3);
   R6: `test_wrong_carrier_refusal…[B-{slab,cylinder,sphere,cart2d}×{alien,
   coupled_field}]` (8, all naming SNBoundaryOperator's unguarded read at
   boundary.py:343); R7: `test_driver_consumes_the_records_own_splitting
   [cart2d-gauss_seidel]` (1). Marker definitions: test_monomorphic_leaves.py:251/
   :263/:276; test_stage_separation.py:83.

## 9. For siblings

- **A (space layer)**: #330's three remainder items (comment 2026-08-07) are the
  live space-side blockers; `_agreed_space` (operator.py:338) is the shared
  commutative space-resolution law both `OperatorSum` and `TensorProductOperator`
  now route through (ERR-076 fix — eager at `__init__` AND recomputed in the
  `domain`/`codomain` properties, operator.py:~2858-2910; catalog entry
  error_catalog.rst:5597 with the mutation-battery counts 18/8/7/3).
- **C (kernel binding)**: apply-time dispatch on S/F persists (#331's measured
  TypeError lists the accepted carrier types); no `.on(V)`; `IntegralKernelOperator`
  is now a category Protocol (§6.8) — re-scope 1.5 against that.
- **D (flux ontology)**: correction #15 is unexecuted and the torsor algebra is
  enforced arithmetic (`_flux_role.py:97`); `coupled_system.py:108,239,292` carry the
  doctrine the plan overturned. #331 is the displacement-space fork.
