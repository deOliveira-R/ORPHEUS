---
name: coupled-block-boundary-unweld-rulings
description: Reusable design-review rulings for a per-system boundary block (B = Σ_x ι_x∘B_x∘π_x) — the precedent DSA/multiphysics will copy; from the SN 1b un-weld review.
metadata:
  type: project
---

Reviewed the SN boundary un-weld (Step 1b, coupled-block campaign, branch
`refactor/sn-walk-unification`, 2026-07-07): the welded `SNBoundaryOperator`
(trace + ψ½ ray hand-fused in one `_apply_faces`) split into siblings
`B_a` (`SNBoundaryOperator`, trace) + `B_b` (`RadialCharacteristicBoundaryOperator`,
ray corner) over ONE `FullFieldSpace`, composed `B = B_a + B_b` via `OperatorSum`.
Explicitly PRECEDENT-SETTING (DSA / multiphysics copy this shape). These rulings
are durable and reusable when the NEXT per-system boundary block lands.

## The rulings that HELD (credit these; don't re-litigate)

- **Present-zero, never None, on blocks a sub-operator doesn't own** — CORRECT, and
  doubly so: (1) FORCED by `FullField._combine_radial_characteristic` (raises on a
  seeded⊕unseeded presence mismatch, so the `OperatorSum.apply` fold
  `B_a.apply(ψ)+B_b.apply(ψ)` would raise if B_a emitted a None ray); (2)
  algebraically it IS the direct-sum embedding `ι_x∘B_x∘π_x` — a block operator is
  the ZERO MAP off its block. "present" ⟺ the block EXISTS in the space; VALUE-zero
  ⟺ this operator doesn't ACT there — two orthogonal facts, correctly separated
  (Pattern 4). The rejected "emit None + special-case the sum" conflates them. When
  reviewing the next block op: verify each sub-operator emits present-zero on the
  complement, and that a gate pins the disjoint reconstruction byte-for-byte.

- **The `isinstance(B, PlainType)` narrowing before a grading verb (`.split()`) is
  HONEST, not a type-smell** — when the triple returns the COMPOSED `B`
  (`PlainType | OperatorSum`) and only ONE specialization path needs to decompose
  it. It faithfully expresses "a solver-internal grading lives on ONE system's block,
  never the composite" (RULING P1). The two alternatives are WORSE: `OperatorSum.split`
  delegating VIOLATES the ruling (grading on the composite); "triple returns the parts
  separately" pushes `B_a+B_b` composition to N consumer sites (Pattern-2 regression).
  The narrowing can never fire in correct operation (its guard condition —
  multi-D-Cartesian — is disjoint from seed-carrying), so it's a defensive
  impossible-state assertion with a ruling-citing message. PASS it.

- **Role grid respected** — every `.apply`/matvec output is a SOURCE/SINK
  (`AngularSourceSink` / `AngularBoundarySourceSink` / `RadialCharacteristicSourceSink`),
  the input trace is parsed as a FLUX. Correct per the institutional role grid.

- **Twin-delivery single-sourced at the operator CORE** — two routes (the
  `OperatorSum.apply` matvec fold vs the in-place `reflect_inflow_inplace` +
  `reflect_corner_inplace` double-call in `_reflect_outflow_into_inflow`) both bottom
  out at `_reflect_trace` (B_a) / `_reflect_corner` (B_b). Two routes = ACCEPTABLE
  (institutional Pattern 1); the standing remedy ("drivers take the whole
  L+C−S−F−B") is tracked. Watch for a THIRD route.

## The recurring FINDINGS to check on the NEXT per-system block

1. **Role-parse-guard SYMMETRY (L-010).** B_a guards `isinstance(psi.boundary,
   AngularBoundaryFlux)` (its slot is role-erased to the base per #289); B_b did NOT
   guard `isinstance(psi.radial_characteristic, RadialCharacteristicFlux)` despite the
   ray slot being identically erased (`RadialCharacteristicField | None`) and a
   flux-role type existing. A per-system boundary block that reads an erased composite
   slot should assert its FLUX role symmetrically — the template must carry the guard
   on EVERY block or DSA inherits the gap. (Milder when the block reads raw values vs
   feeds a realized law operator — note the nuance, still flag.)

2. **"Ruled-kinds" set single-source.** `is_adjointable` returned `kind in
   ("vacuum","reflective")` while `_reflect_corner` re-encoded the SAME set as its
   if/elif/else. Coextensive TODAY ⟹ NIT, but the collapse trigger (rule a 3rd
   law kind) is ON THE ROADMAP (the loud-defer `else` names white/albedo/periodic
   "2.5d plan-of-record"). Bug habitat: implement kind #3 in one site not the other →
   `apply_transpose` false-advertised or unreachable-for-a-supported-kind. Demand a
   `_RULED_CORNER_KINDS` constant both read.

3. **Dead `type: ignore` on the moved/duplicated method** — see lessons L-012
   sharpening. Run the `reportUnnecessaryTypeIgnoreComment` check.

4. **In-diff stale docstring on the delivery helper** — the un-weld updated the
   inline comments but left the `_reflect_outflow_into_inflow` DOCSTRING BODY asserting
   the pre-un-weld "the SAME single SNBoundaryOperator" contract (now B_a + B_b). And
   the module docstring's "SNBoundaryOperator is the A_ss leaf" is now the TRACE-A_ss
   leaf. anti-#11 / Cardinal Rule 3. Reconcile the carve's own file.

5. **Naming asymmetry (judgment call).** B_a kept the generic `SNBoundaryOperator`;
   B_b got the descriptive `RadialCharacteristicBoundaryOperator`. Docstrings frame
   them co-equal ("System A / System B") but the names frame B_a primary / B_b
   auxiliary — which may be HONEST (B_b IS subsidiary: block-triangular, coupled via
   lagged scattering) or may undercut the precedent (DSA adds a 3rd sibling — which is
   "the" SN boundary?). Raise it; let the user rule.

## Architectural opportunities surfaced (file as issues if not tracked)

- **Mesh-identity guard duplicated ~8–9 sites** (streaming ×5, windowing,
  scheduled_invertible, boundary ×2) — the verbatim-modulo-message
  `if psi.bulk.mesh is not mesh: raise ValueError("… mesh-identity invariant …")`.
  B_b's copy FOLLOWS convention (not a new twin), so the fix is a tree-wide
  `_require_mesh_identity(mesh, psi, op_name)` helper, NOT a this-carve blocker.
- **`bc[face].kind == "reflective"` stringly-typed dispatch** (the #186 shim `.kind`
  tag) is a codebase-wide convention → enum (anti-#4). Check if already tracked.

## Step 2 — the A_BA Schur fold un-weld (reviewed 2026-07-08, verdict PASS/ship)

Sibling carve: the ψ½ bulk→ray fold `A_BA = Reconstruction ∘ Emission`, hand-rolled
in 5 places, extracted into ONE `RadialCharacteristicReconstruction(LinearOperator)`
(NEW, `transport/operators/radial_characteristic_reconstruction.py`) + a single-source
weights helper `_radial_characteristic_reconstruction_weights(n_moments, sign)` shared
by the forward `fold_moments_to_radial_characteristic` AND the new `…_transpose`. All
6 design decisions ACCEPTED; zero VIOLATIONS. Durable rulings (reusable for the DSA /
multiphysics A_BA copy):

- **Home follows the DEPENDENCY, not the concept — and it is VERIFIABLE.** A_BA lives
  in `transport/` (not beside its conceptual sibling A_BB in `sn/`) because it produces
  a transport `RadialCharacteristicSourceSink`, is consumed by transport S/F, uses the
  numerics fold, and needs sn only for the SNMesh TYPE. The load-bearing proof is a
  grep: `from orpheus.sn` in `transport/` must be 100% under `TYPE_CHECKING` (verified —
  all 9 are TYPE_CHECKING SNMesh; the carve added 0 runtime sn imports). A_BB stays in
  sn because it wraps the sn sweep engine. Conceptual-sibling-split-across-layers is NOT
  a smell when the `RadialCharacteristic*` name family keeps both greppable.
- **The birth-FACTORY and the OPERATOR are correctly NOT unified (the decision-2 crux).**
  `RadialCharacteristicSourceSink.from_angular_source` (per-level Legendre ANALYSIS of a
  per-ordinate `(N,ng,nx)` source, then fold) and `Reconstruction.apply` (broadcast fold
  of a GIVEN `(n_moments,ng,nx)` moment) share the fold MATH (single-sourced at the
  helper) but are different operations for different callers computing DIFFERENT additive
  terms of q½ (fixed-source/cold-start rhs vs the lagged S/F gain). Unifying would lean on
  the R12a 1-level coincidence (broadcast==per-level only at 1 level) — fragile, correctly
  rejected (Pattern 6). Residual: the 2-line `for sign: cells_view[...]=fold(...)` scatter
  scaffold is duplicated between them — NIT, coextensive-today, collapse trigger = a 3rd
  divergent fold-scatter consumer (or a metric-weight re-home to the cells_view write).
- **`n_moments` FIXED at construction is Pattern 4, not rigidity.** It makes `apply` /
  `apply_transpose` a consistent adjoint pair by construction (the adjoint identity needs
  both ends to agree on the moment dimension); inferring per-call would let apply(n=5) and
  transpose(n=1) desync. default=1 = the isotropic production reach.
- **Broadcast-forward / sum-transpose is the correct adjoint of a copy, and the docstring
  SAYS SO.** apply broadcasts one moment to all (level,sign); apply_transpose sums the
  per-(level,sign) cotangents. Reads like the math.
- **`domain=None` (bare-ndarray moment domain) mirrors K_iso `IsotropicScattering`
  exactly** (space-anonymous, is_adjointable=True, realized transpose) — but A_BA is MORE
  typed (codomain IS the ray space; only the moment input is bare, and `apply` shape-guards
  it). Acceptable deferral; mint a MomentSpace at the 3rd moment-domain consumer.
- **Forward realized here while A_BB defers its forward — JUSTIFIED asymmetry.** A_BA.apply
  is a self-contained fold; A_BB.apply is woven into `(L+C).apply` (extracting standalone =
  twin path), so A_BB raises NotImplementedError until step 4. Both operators DOCUMENT the
  choice with the twin-path rationale.
- **The Pattern-2/7 WIN to credit:** the S-adjoint's hand-rolled `0.5` is retired to the
  shared `…_transpose` (weights read from the ONE helper, so fwd/adj `(2ℓ+1)/2·(±1)^ℓ`
  can't drift). Retirement COMPLETE — all 3 hand-rolled seed loops deleted (grep: no
  `sd_values`/`cells_bar_sum` survive; fold has exactly 3 call sites).

Recurring FINDINGS to check on the next A_BA-family carve:
1. **Line-number cross-refs to now-DELETED loops go stale** (test `_ba_oldloop_reference`
   cites `scattering.py:1592-1596 / fission.py:524-530`; the fold moved 408→447). L-004
   nit — prefer symbol-refs / "reproduced here" over drifting line-cites. Historically
   framed ⟹ follow-up, not MUST-FIX.
2. **Theory-page (docs/theory/) prose for A_BA is absent** (only the build-artifact
   graph.json mentions it). In-code module/class docstrings are exemplary; the Cardinal
   Rule 3 gap is the operator-algebra page section — confirm it's on the step-4 doc list.
3. **apply_transpose reads a role-erased cotangent (mesh-only guard, no Flux-role assert)**
   — CONSISTENT with A_BB `solve_transpose` (values-only Euclidean transpose, the milder
   L-010 case), so NOT a carve-introduced asymmetry. Note, don't block.
4. **Per-matvec construction** of `RadialCharacteristicReconstruction(mesh)` in each seed
   arm — cheap (3 attrs + None-check), reified as a held block at step 4. Micro, not arch.

## Step 3 — A_AB = RadialCharacteristicSeeding (reviewed 2026-07-08, verdict PASS-WITH-NITS)

The ray→bulk seed injection (`RadialCharacteristicSeeding`, added BESIDE A_BB in
`sn/operators/radial_characteristic.py`). CELL-LOCAL ANGULAR (seed at cell i feeds cell
i's ordinate recurrence, NO spatial coupling), so BOTH directions realize HERE as thin
WRAPs of the single-sourced M-M closure (`precompute_psi_state`/`cell_contribution`/
`angular_adjoint`), unlike A_BB whose forward matvec is spatially woven and defers. Code
is SHIP-QUALITY, zero VIOLATIONS. Durable rulings:

- **The WRAP-vs-twin line is RIGHT and identical to A_BB's.** The kernel is single-sourced
  (both A_AB and the in-sweep call the SAME closure methods — verified via the Mode-11
  class-patch spy). Only the `∓numer/V` residual-PLACEMENT orchestration is duplicated;
  coextensiveness is bit-pinned (`array_equal`, 0-ULP) by
  `test_apply_matches_the_in_sweep_seed_injection`. So NIT/CONCERN, never VIOLATION.
- **Isolate-by-linearity is elegant, not a smell.** `apply` zeroes the bulk
  (`precompute_psi_state(np.zeros, radial_characteristic=seed)`) so A_AA's redistribution
  drops out; `apply_transpose` discards `psi_ang_bar`, keeps `seed_cells_bar`. σ-INDEPENDENT
  (bulk=0 ⇒ collision/streaming drop out) — the ctor takes `sn_mesh` ONLY (contrast A_BB's
  σ_t). The isolation premise is POSITIVELY asserted (reference identical at two σ slopes).
- **Codomain typed as the ROLE-ERASED base** (`LinearOperator[RadialCharacteristicField,
  AngularField]`; apply returns the `AngularSourceSink` SUBTYPE — `issubclass` True, so
  covariantly sound). Operators compose on role-agnostic SPACES; the role rides the field
  instance. Consistent with A_BB's single-param `[RadialCharacteristicField]`.
- **The structurally-INDEPENDENT correctness gate handles the L13/Mode-11 blind spot
  EXEMPLARILY.** Bit-identity gates route both sides through the same kernel (blind to a
  shared-method bug — acknowledged in-docstring). `test_euclidean_adjoint_consistency`
  compares `apply` (precompute+cell_contribution) vs `apply_transpose` (angular_adjoint) —
  two separately-implemented duals — + a `cell_contribution`-sign-flip tooth. Credit this.
- **Twin docstring cited the in-sweep by MODULE ref, not line-cite** — dodged the L-004
  drift trap (contrast A_BB which line-cites `4104-4119`). Prefer this.

Recurring FINDINGS (the 4 do-now NITs; check on the A_AA / step-4 carve):
1. **A matvec sub-block (A_AB) mints a NEW tracked-transient twin DISTINCT from the
   resolvent twins.** The RETIREMENT LIST tracked only the A_BB `(L+C).solve` twin
   (`loss_representation:4104-4119`) + `sweep_transpose` (`:4621-4649`). A_AB's twin lives
   in the *matvec* entry point — `(L+C).apply` forward at `loss_representation:3168-3186`
   (`m_full=(denom·ψ−numer)/V`, seed part) + `(L+C).apply_transpose` at `~3475-3590`
   (`numer_bar=-ob/V` → `seed_cells_bar`). The list did NOT enumerate these ⇒ the step-4
   audit could retire the solve-side inline blocks and MISS the apply-side one. **Add
   retirement-list entries keyed by the DIFFERENT production entry point (apply vs solve).**
   This is the headline: a tracked-transient twin whose tracking is incomplete is one edit
   from an untracked twin. (coding-standards "retire as you go" / anti-#11.)
2. **Docstring twin-tracking ASYMMETRY vs A_BB.** A_BB's module docstring has a full
   "Tracked transient twin" section (names the in-sweep block + retire-at-step-4/5 trigger);
   A_AB documents the KERNEL single-source but does NOT name `loss_representation:3168-3186`
   as its placement twin with the trigger. Make it symmetric.
3. **Autodoc REACH gap (Cardinal Rule 3).** `docs/api/discrete_ordinates.rst` `.. automodule`s
   `streaming` + `boundary` but NOT `radial_characteristic` — so the exemplary A_BB + A_AB
   docstrings never reach the rendered API page. PRE-EXISTING since step 1c (A_BB shares it);
   A_AB compounds it; one-line `.. automodule:: orpheus.sn.operators.radial_characteristic`
   surfaces both. (Theory-PROSE for the block operators still deferred to step-4, consistent
   with A_BA finding 2.)
4. **Family-inventory docstring drift.** `sn/operators/__init__.py`'s docstring inventories
   the SN operators (L, L+C, (L+C)⁻¹, B) but omits the ψ½ System-B block family (A_BB, A_AB,
   B_b). Pre-existing; cheap to complete now that the family is being extended.
Role-guard (mesh-only `_check_mesh`, no Flux-role assert) = CONSISTENT with A_BB/A_BA milder
L-010 case; not a carve-introduced asymmetry. No bare `assert` in the gates (vv Mode 8 OK).

## Step 4a — the `SystemRole {A, B, COUPLED}` two-system role lattice (reviewed 2026-07-08, verdict PASS-WITH-NITS)

Foundation carve: a NEW parallel role axis `SystemRole` (numerics/operator.py) — the COARSE
two-system partition (A=transport bulk⊕trace / B=ψ½ ray / COUPLED=off-diagonal), ORTHOGONAL to
the existing `BlockRole {BULK,FULL,BOUNDARY}` (within-System-A refinement). Clones block_role's
mechanism EXACTLY: `SystemRole` enum + `_join_system_roles` (same-role-stays/else-COUPLED, None
propagates) + `system_role: Optional = None` on the Protocol + all 3 composers propagate. Stamps
the four ψ½ blocks (A_BB→B, B_b→B, A_AB→COUPLED, A_BA→COUPLED); every model-generic/System-A op
stays None. Zero VIOLATIONS. Durable rulings (reusable for the NEXT parallel-role-axis, incl. 4d/4e
and any DSA/multiphysics augmentation axis):

- **A parallel enum (NOT a merged member) is STRUCTURALLY FORCED when one operator occupies a
  non-None point on BOTH axes.** The witness here: B_b carries `block_role=BOUNDARY` AND
  `system_role=B` simultaneously (boundary.py:537+540). A single merged enum can't represent a
  2-D classification coordinate; the only single-enum alternative is a combinatorial product
  (`BOUNDARY_A`/`BOUNDARY_B`/…) = the flat-enum-encoding-a-product anti-pattern. Two attributes
  for two orthogonal axes is Pattern 4, not ceremony. This is the decisive PASS argument — check
  for a both-axes operator before accepting/rejecting a parallel-vs-merged enum decision.

- **Cloning the EXACT mechanism at TWO axes is RIGHT; the generic "role-lattice machinery" is
  premature AND typing-regressive.** The two joins ARE the same lattice (non-empty subsets of a
  2-set, identical modulo top FULL↔COUPLED) ⟹ `_join_system_roles` is a coextensive twin of
  `_join_block_roles`. But (a) no forced divergent edit (join law is math; both docstrings share
  the None=conservative-unknown intent) ⟹ NIT not VIOLATION; (b) the generic `_join_role_lattice(
  a,b,*,top)` loses domain-naming (Beck rule 2 > rule 4) and doesn't touch the deeper per-composer
  duplication; (c) the abstraction that WOULD single-source everything — a RoleAxis registry with
  `setattr(self, axis.attr, …)` reflection — is invisible to pyright, REGRESSING the campaign's
  ratchet discipline. Collapse trigger = the THIRD role axis (L-002), not the second. Verdict:
  clone, cross-ref the twin, defer the generic.

- **`None` on a model-generic/System-A leaf is HONEST, not an overload.** None uniformly means "no
  INTRINSIC two-system membership." C/S/F/L/B_a genuinely have none — their System-A membership is
  CONTEXTUAL (assigned by the SN assembly), not a property of a model-generic op. Consequence to
  carry to 4d: a raw `L+C−S−B` sum joins to None (NOT SystemRole.A), so the assembled CoupledOperator
  (COUPLED) AND its A_AA block (A) must be EXPLICITLY-stamped wrappers, never join-derived.

- **A defined-but-unstamped enum member (`SystemRole.A`) is the honest third of a complete
  partition, NOT dead** — when it is live in the join law + tested + has a scheduled stamping one
  step away (4d A_AA block). Fails the L-007 self-conceding-trim on both legs (has a consumer, has a
  scheduled stamp). Define the complete partition once; shipping `{B,COUPLED}` and editing the enum
  at 4d is the higher-churn path.

The recurring NIT family = "the new axis is an INCOMPLETE MIRROR of its established sibling"
(symmetry-in-code / L-010). Check every one on the next parallel-axis carve:
1. **Missing exact-member pin.** BlockRole pins `{r.name for r in BlockRole}=={"BULK","FULL",
   "BOUNDARY"}` (test_operator_protocols.py:58); SystemRole had NO analogue — yet `_join_system_roles`
   hard-codes the 2-atoms+top lattice (`else COUPLED`), so a 4th member silently corrupts the join.
   The parallel-enum decision LEANS on BlockRole's pin but failed to replicate it. Do-now one-liner.
2. **Untested composer leg.** The established `TestComposerRoleDerivation` tests all 3 composers
   (incl. `2.0*L`, `-C` for ScaledOperator); the new propagation test covered only OperatorSum +
   _AdjointOperator, leaving the added ScaledOperator line (operator.py:1757) ungated. Reachability
   proven by the sibling suite. Do-now one-liner (`(2.0*a_ab).system_role is COUPLED`).
3. **One-directional twin cross-ref.** `_join_system_roles`→`_join_block_roles` ref present (good);
   the reverse + a stated collapse trigger absent. Docstring-only.

Credit (got right): composer symmetry COMPLETE (all 3 propagate both roles — L-010 latent class
closed); isinstance markers correctly NOT cloned (Pattern 6 — no consumer yet, natural home at 4d);
test pins the private join algebra directly (better None-propagation coverage than block_role's
composer-only tests); no broken `CoupledOperator` xref (plain literal, not `:class:`).

## Step 4b — A_BB forward realization via SHARED-KERNEL EXTRACTION (reviewed 2026-07-08, verdict PASS-WITH-NITS)

Completed A_BB (`RadialCharacteristicOperator`): realized `apply` (μ∂_r+σ_t matvec),
`apply_transpose`, `inverse()`; flipped `is_invertible`/`is_adjointable`→True. Per USER
RULING "extract the shared kernel NOW over a tracked twin," the forward was NOT reimplemented
standalone — three shared functions added to `sn/spatial/psi_half_angle_seed.py`
(`radial_characteristic_residual_march` verbatim-relocated from the walk's deleted
`_seed_residual_march`; `radial_characteristic_forward_residual` orchestration;
`…_transpose` PURE A_BB^T). Both `A_BB.apply` AND the walk `_seed_rows_forward` now call the
ONE body. Code is SHIP-QUALITY, ZERO VIOLATIONS. This is the PRECEDENT for the DSA/multiphysics
A_BB-forward copy. Durable rulings (all 5 crux decisions PASS):

- **Extraction HOME = the leaf both consumers already import, BESIDE the inverse march
  (symmetry-in-code).** The forward kernels live in `psi_half_angle_seed.py` next to their
  inverse `carlson_inward_sweep_from_source` (forward march ↔ inverse march = neighbors), NOT
  in the operator. Decisive argument: putting the kernel in the operator would force the walk
  (`loss_representation`) to import the HEAVY operator module (LinearOperator/InverseOperator
  machinery) for a pure array function — the leaf is the lightest shared home. `space` param is
  duck-typed (TYPE_CHECKING-only `RadialCharacteristicSpace` import; runtime calls
  `.levels/.cells_view/.corner_view`) — correct sn/spatial→numerics layering, 0 runtime import.
  Same "home follows the DEPENDENCY" law as the step-2 A_BA ruling, applied to a pure kernel.
- **The transpose SPLIT across two places is the correct BLOCK decomposition, not a smell.**
  Walk `_seed_rows_transpose` = shared PURE A_BB^T (`…_forward_residual_transpose`) + the A_AB^T
  `seed_cells_bar` coupling added in the walk; `A_BB.apply_transpose` wraps ONLY the pure part
  (verified: no `seed_cells_bar` in the operator). A_BB^T and A_AB^T are DIFFERENT operators;
  each single-sourced; the walk's full System-B-row transpose = A_BB^T + A_AB^T composed
  (Pattern 2 + 5). The split IMPROVES A_AB^T isolation (its `seed_cells_bar` add is now visually
  separated, ready for its own step-3-tracked retirement). NOT a new twin — `seed_cells_bar` is
  the SAME tracked A_AB matvec-side twin from step-3 finding #1.
- **Keeping the thin walk methods is JUSTIFIED — for two DIFFERENT reasons.**
  `_seed_rows_forward` (one-line delegation) is load-bearing: `test_282_direct_seed_fixed_point.py:373`
  `mp.setattr(_OneDimScanWalk,"_seed_rows_forward",_flipped)` is a live mutation-teeth seam;
  retiring the method forces rewiring that test to the shared kernel, which would ALSO patch the
  operator's apply (broader blast, conflates two callers). `_seed_rows_transpose` is NOT
  vestigial either — it's genuine composition (shared A_BB^T + inline A_AB^T), not a pure
  delegation, so it earns a name regardless of the seam. Its docstring DOCUMENTS the keep-reason
  (anti-#11 done right).
- **`apply` returning the CONCRETE `RadialCharacteristicSourceSink` (base declares
  `-> Codomain`=`RadialCharacteristicField`) is covariant NARROWING, verified honest.**
  `RadialCharacteristicSourceSink(RadialCharacteristicField)` IS a subclass ⟹ a narrower return
  is LSP-sound; pyright 0 errors on the file (Leg-3 live-tree check). Role grid respected:
  apply/apply_transpose output = SOURCE/SINK (Aψ), solve returns FLUX, solve_transpose returns
  SourceSink. Identical to the A_AB (step-3) precedent.
- **`inverse()` via GENERIC `InverseOperator` is the correct taxonomy — because A_BB.solve is a
  DIRECT, seed-INDEPENDENT exact march.** It matches the InverseOperator contract (`apply` does
  `del initial_guess; return inner.solve`), NOT the iterative-SEEDED `SweepOperator` (which
  threads `initial_guess` into the curvilinear Carlson closure). Minting a distinguished type
  needs a distinguishing invariant a CONSUMER exploits — none exists (Pattern 6 / type-vs-property).
  The pre-guard `if not self.is_invertible: raise NotInvertible(...)  # pragma: no cover` FOLLOWS
  the `MultiplicationOperator.inverse` family convention exactly (verified) — dead here only because
  A_BB's is_invertible is hardcoded True (vs the value-dependent leaves), honestly pragma-marked.
  Do NOT flag the guard: family-uniform spelling is Beck-rule-2.

Credit to reinforce: the Mode-11 anti-twin spy (`_install_forward_spy`) patches BOTH namespaces
(`_rc_mod` + `_lr_mod`) → the both-callers-one-body single-source proof; the shared-kernel
docstring NAMES both callers with `:meth:` refs and the code ACTUALLY calls it (verified by spy —
the INVERSE of my recurring L-001 catch, credit it); test uses ψ0=solve(q0) consistent subspace
for solve∘apply=id (author understands the +1-outflow-corner free-datum null structure).

The ONE do-now NIT (the recurring in-diff-staleness sub-pattern — RECORD IT):
- **Module docstring UPDATED, CLASS docstring left STALE.** `radial_characteristic.py:176-179`
  (the class docstring, NOT in the diff) still says "the step-1c realization scope … the forward
  `apply`/`inverse`/`is_invertible` land in step 4" — directly contradicts the code 10 lines below
  AND the just-rewritten module docstring ("realizes BOTH directions"). anti-#11 false-invariant
  doc / Cardinal Rule 3 / L-004 but IN the carve's own file ⟹ do-now, not follow-up. **Recurring
  tell for the NEXT block-op carve: when a carve rewrites the MODULE docstring's scope section,
  grep the CLASS + method docstrings of the same file for the OLD scope language — the module/class
  docstring pair drifts because they restate the same scope contract in two places** (itself a mild
  Pattern-2: the "realization scope" is stated in both; the class one should just `:ref:` the module
  one, never restate). Bug habitat: a maintainer reads the class docstring (shown first), believes
  apply is deferred/NotImplementedError, and re-adds a guard or writes a twin "to match the doc."

## Step 4c — THE LIFT: S/F→pure bulk, A_BA=`RadialCharacteristicEmission` own gain (reviewed 2026-07-08, verdict CONCERNS/ship-after-retirement)

The scatter LIFT: model-generic S/F stopped hand-rolling the ψ½ ray; the bulk→ray coupling
posed as `A_BA = RadialCharacteristicEmission (Fold ∘ K ∘ integrate)`, applied by the SI/Krylov
driver as its OWN lagged gain (the #208 B-from-S pattern extended). ZERO elegance-VIOLATIONS
(3-leg bar); all 5 design decisions SOUND. Durable rulings:

- **A_BA as `FullField→FullField` gain is FORCED, not a smell — verified at the fold.** The driver
  folds gains via `for g in self.gains: rhs += g.apply` (SI, iteration.py:620) / `out -= g.apply`
  (Krylov, :870) over the FullField iterate — so a peer gain MUST be FullField→FullField. Injecting
  A_BA as a 3rd gain IS the honest `(L+C−S−A_BA−B)ψ`. Mirrors B_b exactly (present-zero bulk/boundary
  = the ι_B∘core∘π_A embedding; reads only `psi.bulk`, its true block domain). Does NOT muddy the
  (B,A) identity — the un-embedded core is trivially recoverable. Correctly typed for its CURRENT
  (gain) consumer per Pattern 6; re-type at 4d when the CoupledOperator exists.
- **LATENT 4d reconciliation (Arch-Opp, not a blocker):** the four blocks now split across TWO typing
  conventions — A_BB (`[RadialCharacteristicField]` ray→ray) + A_AB (`[Ray, Angular]` sub-space) are
  BLOCK-typed; A_BA + B_b are FullField-embedded. 4d's 2×2 assembly must reconcile. Flag so 4d isn't
  surprised; forcing block-typing now = premature abstraction toward a non-existent consumer.
- **`_lagged_gains(S,B,mesh)` at the gain-assembly seam (NOT widening `_within_group_triple`) is the
  RIGHT clean-before-extend.** Single-sources the A_BA injection (all 4 driver sites route through it:
  2 Krylov direct, 2 SI via `_select_si_resolvent`→`_lagged_gains`). Variable-arity gains (2 seedless /
  3 seed-carrying) CANNOT be a fixed 4-tuple, and the G-S arm needs S,B separately for `.split()` — so
  the triple rightly stays `(L+C,S,B)` and A_BA layers on. Avoids the 99-ref rename. **Recurring do-now
  NIT:** the peer-gain add left the gain-ENUMERATING docstrings/comments stale — solver.py:194-195
  (`_within_group_triple` matvec `(L+C−S−B)` + rhs `q_ext+S+B`), :376/382-384 (`_within_group_krylov`),
  :645 (`_within_group_si`), :1394/:1564 (inline rhs) all omit A_BA on carrying meshes (anti-#11 / L-004,
  blast lands OUTSIDE the diff). Latent habitat: a NEW driver site unpacks the triple, forgets
  `_lagged_gains`, silently drops A_BA on the sphere.
- **`emission_kernel` genericity SOUND (user-ruled machinery + necessary DI) but the DOCSTRING
  overstates fission as a live consumer (L-001).** The operator composes an injected kernel (`Fold ∘ K
  ∘ integrate`) — DI is necessary, and K_iso vs the fission dyad are two real emission kernels, so keeping
  it kernel-agnostic is honest "build the machinery." BUT the param docstring (radial_characteristic.py:
  1143-1144 + class-doc :1086/:1097) says "pass `fission_op.kernel` for the fission coupling" — the plan
  (commit 2, SAME session) routes fission's PRODUCTION coupling AROUND this operator: a direct
  `Reconstruction(fission_source[None])` fold at the OUTER q_ext seam, because fission's K∘integrate is
  pre-computed as the eigenvalue outer source (HAZARD 5 — S=within-group-lagged, F=outer). So
  `Fold∘K∘integrate` is the SCATTER shape; fission uses `Fold` alone elsewhere. The "2nd consumer" is a
  phantom for THIS operator. Bug habitat: a maintainer wiring `RadialCharacteristicEmission(mesh,F.kernel)`
  into fission DOUBLE-APPLIES K∘integrate (outer source already has it). Fix = soften the doc: "F.kernel
  smoke-verified as machinery; fission's production coupling rides the outer q_ext seam." CONCERN, not
  VIOLATION (no divergence today).
- **Present-zero ray on the S/F drops = CORRECT (re-affirms step-1b).** FORCED by the composite
  mixed-presence law (`FullField.__add__` raises on mixed ray presence, so `S.apply+A_BA.apply` needs both
  present) + IS the direct-sum embedding. No cleaner spelling.
- **Pattern-7 WIN to credit:** the moved S-adj transpose swapped `self.weights`→`mesh.quad.weights`
  (both = `quadrature.weights`, verified) — the more canonical source since A_BA holds the mesh, not S.
  Docstrings read like the math (`apply`/`apply_transpose` bodies are literally `Fold∘K∘integrate` and its
  reverse, named intermediates φ0/q0/ray); `is_invertible inherits False` verified accurate (base:621);
  my step-3 autodoc-reach finding was ADDRESSED (discrete_ordinates.rst:137 automodules the sn module).

**THE headline (retirement incomplete — MUST complete before commit):** the migration moved
`RadialCharacteristicReconstruction` transport→sn but the OLD file
`orpheus/transport/operators/radial_characteristic_reconstruction.py` STILL EXISTS as a byte-faithful
twin (same class/`__init__`/apply/transpose) — zero production importers, one TEST importer
(`test_psi_half_coupling.py:74`). The plan's own commit-1 list says "retire the transport
radial_characteristic_reconstruction.py". Two identical classes coexist = anti-#1 twin habitat (a future
fold-weight edit lands on one). NOT an elegance-VIOLATION (coextensive today → NIT) but an aggressive-
retirement FLOOR breach + plan-mandated commit-1 deliverable ⟹ hard Approval Condition: delete the old
file + coordinate the `_rcr_mod` test repoint (test agent) so they land together.

## Step B.2a — the N-system `CoupledField`/`CoupledSpace`/`CoupledOperator` substrate (reviewed 2026-07-10, verdict PASS-WITH-NITS)

The semantics-agnostic numerics-layer block machinery (`orpheus/numerics/coupled_system.py`, NEW,
uncommitted on `refactor/sn-walk-unification`): typed N×N block grid over present-`None` sparsity,
the block vector, the direct-sum space. Chosen over flat `OperatorSum` + present-zero padding (padding
keeps wrong multiplications representable — user-decisive). ZERO VIOLATIONS; 39/39 green under -O;
pyright 0 on the file (ratchet unaffected — new-file-at-0). Durable rulings (reusable for B.2b–d instance
wire AND the DSA/multiphysics 2×2 copy):

- **THE headline twin question — `CoupledField` is NOT a duplicated `Composite` algebra.** The ψ½
  instance is `CoupledField[FullField, RadialCharacteristicComposite]` where `FullField` IS-A
  `Composite` — so `CoupledField` NESTS `Composite`, it is not parallel to it. Different arity (N-tuple
  vs `Composite`'s fixed interior⊕boundary + `_recombine`-extensible 3rd block), different layer
  (inter-system vs intra-system), and the actual arithmetic + affine-torsor role law single-sources on
  the LEAVES — both containers merely FAN OUT `op` across sub-parts (neither computes any arithmetic).
  Unifying them MULTIPLIES concepts (a generic over block-count-and-namedness), failing the
  concept-count test. PASS. Collapse trigger: a THIRD direct-sum fan-out container, or a refactor making
  `Composite` itself N-ary. **When reviewing the next block-vector: check whether it NESTS or PARALLELS
  the existing composite — nesting is not a twin.**

- **`eq=False` on an ndarray-bearing block vector is CORRECT and follows the `Field` convention, NOT a
  divergence to reconcile.** `numerics/field.py:142` is itself `@dataclass(frozen=True, eq=False)`.
  `CoupledField(eq=False)` is MORE defensive than transport `Composite` (eq ON) — `Composite`'s auto-eq
  is safe ONLY because it is narrowly typed to identity-eq `Field` leaves; a general aggregate accepting
  any `SystemField` (some members carry array-valued eq — the test toys do) needs `eq=False` or `==`
  raises "truth value ambiguous". The three ndarray-frozen-dataclass eq strategies in the tree are all
  legitimate for their kind: SPACES → explicit `(name,shape)` identity eq; STATE VECTORS → `eq=False`
  (no identity tuple exists); `Composite` → auto-eq (safe by narrow leaf typing).

- **The block Hilbert adjoint comes FREE and it VERIFIES (Mode-12).** `CoupledOperator` implements ONLY
  the Euclidean transposed grid (`apply_transpose`) + `CoupledSpace` member-wise metrics; `A.H` is the
  generic `_AdjointOperator` doing `G_V⁺ ⊙ apply_transpose(G_W ⊙ y)` via `codomain.apply_metric` /
  `domain.apply_inverse_metric` (operator.py:1216). ZERO adjoint code in the block operator. Verified
  live: the M-ADJ-metric tooth monkeypatches out the conjugation → reciprocity reds O(1). This is the
  Pattern-1∘Pattern-2 template every future block operator must copy — a hand-rolled "Euclidean block
  .H" is the ERR-067 reopening. Guard-narrow `apply_transpose`/`assemble` (per-block `adjointable`/
  `assemblable` TypeGuards, `all(present)` predicates, block-naming raises) faithfully mirror
  `OperatorSum`; `assemble` via `sparse.block_array` with the all-None-row/col construction guard
  DOUBLING as block_array's size-inference precondition (Pattern 4).

- **No `_recombine`-style subclass hook on `CoupledField` is the RIGHT defer-until-2, not a Liskov
  trap.** No second `CoupledField` flavor exists (the instance uses `CoupledField` directly; semantics
  ride the members). `replace(self, systems=...)` returns `type(self)` → Liskov-safe-by-default for a
  hypothetical subclass; the hook lands when a second flavor needs transform semantics (mirroring
  `Composite._recombine`). NB the endo/non-endo split is correctly observed: metric methods use
  `replace(x,...)` (endomorphism, preserve x's type); `apply`/`apply_transpose` construct fresh
  `CoupledField` (domain→codomain, possibly different). `reduce(lambda a,b:a+b, contributions)`
  correctly avoids `sum()`'s int-0 seed footgun.

The ONE recurring NIT (L-001 again — my #1 catch): **`_is_system_field` / `_member_from_flat` open-code
iteration.py's `_is_ravellable` / `_unravel_like` that their OWN docstrings NAME** ("the _is_ravellable
pair", "mirroring _unravel_like"). Coextensive today → NIT not VIOLATION. But the fix is NOT
"import the private now" — `iteration.py`'s `_is_ravellable`/`_unravel_like`/`_ravel` are ad-hoc
`_`-privates the codebase is MIGRATING AWAY from (vector.py promotes the arithmetic half to the named
`runtime_checkable Vector` protocol, explicitly excluding to_flat/copy). So per L-002 (two coextensive
spellings acceptable; collapse at the THIRD divergent consumer OR when the concept gets a named home):
the ravellable-PAIR member concept wants its own named Protocol the way `Vector` got one for the
dunders — `SystemField` is the first step toward that home. Bug habitat: a ravellable-contract change
lands on iteration.py and the CoupledField member-gate/reconstruction silently diverges from the Krylov
boundary that actually consumes CoupledField (where ERR-053 sizing lives). Architectural opportunity to
track: converge `SystemField` + `_is_ravellable` + `Vector` (three overlapping member-contract concepts
in `numerics`) — Cardinal Rule 2. Doc-accuracy micro-nit: SystemField's "narrower than any honest
Protocol spelling" overstates (Vector honestly spells the simple `Self+Self→Self` dunders and the
members satisfy it); the accurate reason for duck-typing is the affine-torsor signatures
(flux+displacement→flux) + the `Any`-lambda `_map_binary` posture.

## Step B.2d d1 — the block-native driver + `A = M − N` system record (reviewed 2026-07-11, verdict PASS-WITH-NITS)

The de-twin that collapsed the TWO within-group construction sites (`_within_group_triple`
in solver.py ⊕ `build_coupled_system` in coupled_system.py — the "sanctioned transient twin"
of B.2b/B.2c) into ONE: `build_within_group_system(sn_mesh, mat_xs, *, scattering_op=None,
scattering_order=0) -> WithinGroupSystem(loss, space, resolvent, gains)` — the named Hackbusch
regular splitting `A = M − N`. `build_coupled_system` reduced to a loss-VIEW delegation. The
SI/Krylov drivers go block-native: iterate = `CoupledField[ψ_A 3-block, ψ_B]` on carrying
meshes. M = `CoupledInvertibleOperator` (the fused (L+C) walk behind a pack/split bridge);
M⁻¹ = `CoupledSweepOperator` (`InverseWrapMixin` sibling). The B.2b FullField-gain adapters
(`_RayEmissionFullFieldGain`/`_RayBoundaryFullFieldGain`) RETIRED. ZERO code VIOLATIONS;
pyright 0/0 on all 6 files; 0 `# type: ignore` / 0 bare `assert` added. Durable rulings:

- **The `A = M − N` record single-sources at the PIECE level, not the operator level — and
  that is the correct Pattern-2 shape.** `loss` grid `[[A_AA,A_AB],[-emission,A_BB]]`, `M`
  (`CoupledInvertibleOperator(LC)`), and `N` (`CoupledOperator[[S+B_a,∅],[emission,B_b]]`) are
  THREE compositions of the SAME shared objects (one `LC`, `S`, `B_a`, `emission`, `B_b`,
  `A_BB`). Verified the algebra closes: `A_AA=LC−(S+B_a)=M_AA−N_AA`; `A_AB=Seeding−0`;
  `−emission=0−emission`; `A_BB−B_b=A_BB−B_b`. Primitives compose; the three views can't
  drift because they reference one object set. NOT a twin. (The M-vs-grid ACTION equivalence
  — fused joint recurrence vs explicit block sum — is a pre-existing principled-equiv row
  ~5.5e-16 pinned by G-c2.3, NOT a construction twin; don't conflate.)

- **The pack/split bridge is ONE systematic adapter, not a twin.** `CoupledInvertibleOperator`'s
  four surfaces each do `_split_fused_state(self.fused.<op>(_fuse_coupled_state(x)))` — pack→
  act→split, uniform across apply/solve/apply_transpose/solve_transpose. `fused` (public, read
  by the wraps-predicate gate) is high-signal: "the fused (L+C) walk," distinct from
  `resolvent`(=M) and the mixin's `inner`(=M, on the M⁻¹ side). The distinct wrap-names are
  HONEST (they wrap different things) — do not flag as inconsistency.

- **`CoupledSweepOperator` is a faithful `InverseWrapMixin` sibling of `SweepOperator`** —
  identical back-half (`del initial_guess  # direct exact inverse (#282/2.5d); return
  self.inner.solve(rhs)`), same ctor-guard-is-the-TYPE pattern. The `Coupled` prefix keeps the
  family greppable (`grep InvertibleOperator`/`SweepOperator` finds both). Credit; copy for DSA.

- **`[Any]` on the driver returns (`SourceIteration[Any]`/`KrylovAcceleration[Any]`, widened
  from `[FullField]`) is INVARIANCE-FORCED, not a dodge.** V is invariant-bound to `Vector`;
  the coupled/seedless heteromorphism (CoupledField vs TimedFullField) is not statically
  discriminable from `WithinGroupSystem.resolvent`'s `Coupled…|Invertible…` union; consumers
  re-narrow with loud isinstance parses. Matches the C2 `default=Any`-under-invariance
  precedent + L-009 honest-cast-with-loud-renarrow. pyright 0/0 confirms. PASS it.

- **Transitional hygiene is the EXEMPLAR to hold up** (anti-#11 done right): EVERY transient
  names its retirement point — the 4 bridge helpers + the 4 birth-seam `if coupled:
  _split_fused_state` comments → d2; the `CoupledInvertibleOperator` internal bridge → 4e;
  deleted symbols left navigational tombstones pointing to the successor. Illegal-states loud:
  `_system_b_member` REFUSES a live-ray fused state (the dead-slot double-count hazard); the
  seedless SI arm parses `S`'s type; windowed-reconstruction asserts non-coupled structurally.

- **Numerics NOT contaminated with ψ½ vocab.** `_flux_displacement_leaf` recurses on
  `systems[0]` — numerics-native direct-sum vocabulary, duck-typed, no transport import (the
  "System A" mention is an EXAMPLE). The GMRES exact-breakdown carve-out CODE is a pure
  `residual_history[-1]==0.0` predicate (generic, permanent, correct). Clean.

The recurring FINDINGS (check on d2/d3 and the DSA copy):
1. **THE headline — deletion blast radius = 13 broken `:func:`_within_group_triple` doc xrefs**
   (operator_algebra.rst:1750,5025,5133,5237,5329,5722,6736; discrete_ordinates.rst:12233,
   13820,14746,14820,14842,17125). Mix of `:func:` (broken xref, no -W warning per Cardinal
   Rule 3 / aggressive-retirement 3rd search) and ``literal`` (false present-tense prose). The
   present-tense ones ("the within-group builder `_within_group_triple`", "returns the (L+C,S,B)
   triple") are higher-severity — a maintainer navigating them re-ADDS the retired symbol,
   re-opening the twin. Tracked for d3 archivist per the plan spine, but ENUMERATE them or fix
   in-commit (find-replace `_within_group_triple`→`build_within_group_system` + tense). This is
   THE recurring doc-gate finding of every retirement carve in this campaign (step-4c headline,
   #226 C5) — expect it again on d2 (the FullField/mixed-presence-law retirement).
2. **The exact-breakdown guard's PERMANENCE rides a TRANSIENT trigger.** iteration.py GMRES
   carve-out is a correct GENERAL guard, but its comment cites only the ψ½ dead-pad as the
   trigger ("the padding retires at the d2 eviction"). Post-d2 the pad is gone; the guard risks
   becoming untriggered-and-untested, and a d2 maintainer may retire it WITH the pad. Anchor it
   with a minimal general singular-consistent GMRES test, or reword so the general invariant is
   load-bearing and the ψ½ pad is "first caller, not the reason." (A permanent guard whose only
   demonstrated trigger is transitional = a sub-species of the fuller-view-oracle question.)
3. **`_fuse_coupled_state` seam asymmetry (L-010, mild, transient).** It funnels System A through
   `_system_a_member` but inlines `RadialCharacteristicComposite.require_member(systems[1])`,
   bypassing the sibling `_system_b_member` seam whose stated purpose is "one funnel per system."
   Defensible (require_member = clean non-optional) but mixed; funnel both or inline both.
4. **`build_coupled_system` docstring overstates current consumers (L-001 shape on consumers).**
   New delegation doc says present-tense "residual evaluation, assembly, the DSA substrate take
   this pair" — grep = ZERO production callers (only tests). They are LATENT/planned. Future-tense
   the prose; keep the function (the loss-view IS the right forward shape).
5. **Pre-existing LC-spelling TRIPLICATION persists** (`build_within_group_system` builder +
   `self.L` @ solver.py:849 init + :929 rebind_sigma_t), documented at operator_algebra.rst:1748.
   NOT introduced/worsened here; tracked for d3 estate. The 3rd spelling IS the L-002 collapse
   trigger — confirm d3 collapses the two `self.L` legs onto the builder (or a shared `_build_lc`).

## Step B.2d d2 — THE ATOMIC RAY EVICTION (reviewed 2026-07-11, verdict PASS-WITH-NITS)

The ψ½ ray EVICTED from the SN composite: `FullField` collapsed to a pure 2-block
`Composite[BulkField, BoundaryField]` (inherits the WHOLE algebra — 6 dunders, flat protocol,
zeros, copy, _recombine — from the base; adds ONLY the concrete-locus `__post_init__` guards;
~250 lines of ψ½ machinery + `_seed_space_for`/`radial_characteristic_space` field DELETED).
System B = its own `RadialCharacteristicComposite`, coupled via `WithinGroupSystem`. The six
walk signatures re-typed to EXPLICIT leaf kwargs; the transitional `_split_fused_state`/
`_fuse_coupled_state` bridge DELETED, replaced by `_require_coupled_pair` (parse) +
`_coupled_flux_state`/`_coupled_source_state` (native birth). ZERO VIOLATIONS. The d1
transitional-hygiene promises (bridge dissolves at d2, ray eviction) all DISCHARGED.

Durable rulings (reusable for the DSA/multiphysics 2×2 eviction — the whole campaign is
precedent-setting):

- **The typed-leg protocol is `input=role-erased-base` / `output=concrete-role-class`, and it
  is CONSISTENT across all 8 signatures — verify it, don't flag the apparent base/concrete
  asymmetry.** Which leg is base LOOKS inconsistent (loss_action: flux=base, source=concrete;
  solve: source=base, flux=concrete) but is CORRECT: the READ leg is the permissive base, the
  FILLED buffer is the concrete role class matching the surface's OUTPUT role. matvec reads flux
  → writes SourceSink; solve reads source → writes Flux; transposes read cotangent → write
  SourceSink. I nearly flagged a false inconsistency — the input/output role SWAP between
  matvec and solve is what swaps which leg is base. Role-honest, not a NIT. CREDIT.

- **The (A,A) block action's zero-substitution is FORWARD-ONLY and single-sourced; the transpose
  DISCARDS, and that asymmetry is symmetry-in-MATH.** M = [[L+C, +Seeding],[0, A_BB]]. Forward
  y_A = (L+C)x_A + Seeding·x_B, so the (A,A) action MUST feed a zero seed (zero x_B) — ONE
  spelling at `loss_representation:3158` (`seed_field = RadialCharacteristicFlux.zeros_on` when
  no leg + carrying). Transpose x̄_A = (L+C)ᵀ ȳ_A is INDEPENDENT of seed_cot ((B,A)=0), so it just
  skips the pullback (`if chi_seed is not None and seed_cot_out is not None`). No second zero
  spelling leaked (all other `RadialCharacteristicFlux.zeros_on` are OUTPUT-buffer allocs, a
  different concept). When reviewing the next eviction: confirm the off-diagonal ZEROING lands on
  the input-read side of the forward and is a DISCARD (not a substitution) on the transpose.

- **Two guard families, correctly separated: the R12a biconditional (`_require`/`_refuse`,
  single-sourced helpers) vs the both-or-neither PAIR-XOR (inline ×2).** matvec allows both-or-
  neither (the (A,A) block action is a legal no-legs call); solve REQUIRES both (no ray-decoupled
  inverse spelling — the joint M⁻¹ is the only inverse). This matvec-vs-solve distinction is
  PRINCIPLED and documented. The XOR `(a is None) != (b is None)` is the one inline duplicate
  (NIT-2 below).

- **`is_adjointable` third factor = eager-refuse-at-predicate (illegal-states-unrepresentable).**
  `SweepOperator.is_adjointable` gained `and self.inner.sn_mesh.radial_characteristic_space is
  None`: a carrying-mesh fused wrap's leg-less `apply_transpose(b)` can't thread seed_cot, so the
  honest joint home is `CoupledSweepOperator` — refuse EAGERLY, not mid-apply. Safe attr chain
  (guarded by the `isinstance(self.inner, InvertibleOperator)` conjunct first).

- **The refactor DISSOLVED the `_system_b_member` live-ray refusal — correct (the guard the type
  now makes unnecessary).** Post-eviction a live-ray ψ_A is unrepresentable, so the B.2c double-
  count guard branch is DELETED, not kept. "The type system is the guard" (memo R3). This is the
  positive form of illegal-states-unrepresentable — CREDIT, don't ask for the guard back.

- **`boundary_vs_interior_split` identity √(b²+i²)=‖r_A‖ is now EXACT** (pre-eviction 3-block
  silently excluded seed rows). The split residual pair (Interior/Boundary, mirroring the B.2b
  SourceSink split) dissolves the retired unified leaf's corner-units deviation — each locus leaf
  declares its own honest units (ANGULAR_RATE interior / ANGULAR_FLUX corner). CREDIT.

Retirement COMPLETE at code level (verified live): 0 orphaned `.radial_characteristic` block
reads (every hit is a module-path `:class:`/`:meth:` ref, an import, or the NEW `Solution.
radial_characteristic` member); 0 refs to the deleted `_split_fused_state`/`_fuse_coupled_state`/
`_seed_space_for`/`_zero_radial_characteristic_like`/`_radial_characteristic_scaled`/
`radial_characteristic_space=`; all 3 SN `Solution(...)` sites pass the new biconditional member;
doc source updated (theory torsor-table row → split pair; verification matrix regenerated). The
d1-predicted doc-xref blast radius did NOT recur in SOURCE docs (only stale `_build/` artifacts +
one past-tense prose mention — both benign).

THE finding to carry forward (CONCERN, d2-INTRODUCED, the recurring eviction trap):
- **`evaluate_residual`'s BARE arm has no carrying-mesh refusal — re-opens Mode-12 blindness at
  the INPUT boundary.** Pre-eviction the 3-block FullField ALWAYS carried the seed (mesh-keyed
  presence), so the seed residual was always computed. Post-eviction a caller can pass a bare
  System-A `FullField` on a CARRYING mesh and silently get a System-A-only residual (r_B dropped)
  — exactly the Mode-12 (b) blindness the interior-residual file's OWN docstring says it exists to
  prevent. Not wrong today (only caller is the coupled-arm test), but the docstring names the
  future DSA #2 consumer. Fix mirrors the campaign's own posture: the bare arm refuses when
  `q_ext.interior.mesh.radial_characteristic_space is not None` ("pass the coupled pair"). This is
  the eviction's signature footgun — an evicted concept becomes silently OMISSIBLE at a consumer
  that used to get it for free. CHECK THIS on every future block eviction.

Two do-now NITs (coextensive-today, collapse triggers named):
- **NIT (4× shell in `CoupledInvertibleOperator`):** apply/solve/apply_transpose/solve_transpose
  each do parse(`_require_coupled_pair`)→alloc-buffer→`self.fused.<surface>(...)`→repack
  (`CoupledField(systems=(result_a, RadialCharacteristicComposite.from_unified(buf)))`). Parse +
  repack are byte-identical ×4; the MIDDLE genuinely varies (kwarg names flux/source vs
  seed_cot/_out, buffer type, IN/OUT role) so a callback-helper trades dup for indirection (Beck
  rule 2 ≥ rule 4 — a wash). ACCEPTABLE-as-is. Collapse trigger: a metric-weighted ψ_B re-home
  or a `.copy()` landing on the repack → then fold a `_pack_b(result_a, buf)` (repack line only).
- **NIT (both-or-neither XOR ×2):** the pair-guard is inline in `loss_action` + `loss_action_
  transpose`; the R12a `_require`/`_refuse` guard family is ALREADY single-sourced in the same
  module, so the XOR is the third member that should join it (`_require_leg_pair`). Cheap.
- **Micro-NIT (greppability):** forward legs `radial_characteristic_{flux,source}` vs transpose
  `seed_cot{,_out}` — a grep for the family prefix misses the transpose legs. Both "seed" and
  "radial_characteristic" are established ψ½ vocab, so low-severity.

## Step B.2d d3 — the LC-collapse ESTATE + the B_aᵀ vacuum transpose fix (reviewed 2026-07-11, verdict PASS-WITH-NITS)

Two production changes on the uncommitted diff: (1) the LC-triplication collapse — the L-002
trigger I flagged at d1 finding #5 — LANDED: `build_streaming_collision(sn_mesh, mat_xs)`
(coupled_system.py:559-580) is now THE one `StreamingOperator + MultiplicationOperator[σ_t]`
spelling, called by `build_within_group_system` + both `SNSolver` legs; grep-confirmed the ONLY
inline LC in the tree. (2) the `SNBoundaryOperator._reflect_trace` transpose fix. ZERO
VIOLATIONS; both ship-quality. Durable, reusable rulings (DSA/multiphysics copies these):

- **LC-collapse home = the ASSEMBLY layer, not the operator-vocabulary layer.** `build_streaming_collision`
  correctly lives beside the `build_*` family (its only consumers) in coupled_system.py, NOT in
  streaming.py beside `StreamingOperator`/`InvertibleOperator`. Decisive: streaming.py is the operator-CLASS
  vocabulary and has NO `MaterialXSField` dependency; hosting a `(SNMesh, MaterialXSField)`-consuming builder
  there inverts the layering. "Home follows the dependency" — coupled_system.py already carries every runtime
  dep. Name is family-consistent via the `build_` prefix; correctly NOT `_system` (it builds a FACTOR/resolvent
  core, not a system). Pure extract-function (3 prior spellings byte-identical), so rebind live-read preserved.

- **THE reusable transpose ruling: `(P_sel ∘ law)ᵀ = lawᵀ ∘ P_sel` — a projected-forward operator's
  transpose restricts the INPUT, it does NOT project the OUTPUT.** The bug: forward `B_face = P_inflow ∘ law`
  (output-projected onto the consistency row); the old transpose spelled `P_outflow ∘ lawᵀ` (OUTPUT projection
  of lawᵀ) — the transpose ONLY when `law` is an inflow↔outflow PERMUTATION (reflective/albedo). For a
  DIAGONAL-mask law (vacuum = zero-on-inflow ⊕ id-on-outflow) it extracted the harmless id-on-outflow block
  into a spurious `+1` where the forward is the ZERO map. **The tell (institutional): every permutation-law
  fixture stays GREEN over the bug — permutation laws coincide bit-identically under both spellings, so only a
  DIAGONAL/non-permutation law (vacuum) exposes it.** When reviewing any projected-forward operator's transpose,
  demand the input-restriction spelling and demand a NON-permutation-law gate (a dense `Bᵀ==B.T` on a vacuum/
  albedo face), never trust the reflective reciprocity gate. Fix byte-preserves the `apply` branch (only the
  buggy transpose arm changed); the unified `sel` serves whole-face AND row-restricted (`rows`) cases.

- **The per-face `adjointable(law)` TypeGuard-raise is Pattern-4 + anti-#19 done right, and it is
  single-sourced with the composite `is_adjointable`.** `adjointable()` (numerics/operator.py:1125) is a
  `TypeGuard` returning `op.is_adjointable`; it narrows `law` so `law.apply_transpose(masked)` types WITHOUT a
  `# type: ignore`. The raise is that TypeGuard's mandatory else-arm, bottoming out at the SAME `law.is_adjointable`
  the composite `all(law.is_adjointable ...)` aggregates — not a second spelling. Lazy-at-apply placement is
  family-consistent with B_b (`_reflect_corner` raises lazily too, boundary.py:621). Comment honestly hedges
  "unreachable via .H's eager gate, but loud if a caller bypasses the predicate" — accurate (the Euclidean
  `apply_transpose` has NO eager gate at its entry; the per-face raise IS the primary refusal there).

- **THE do-now NIT (recurring): a newly-added fail-loud `MissingAdjoint`/refusal raise ships WITHOUT a
  bite-test even when the fixture already exists.** Verified across the WHOLE live `tests/` tree (L-011): no
  test hits `SNBoundaryOperator.apply_transpose`'s raise, though `test_adjointability_drops_when_a_face_lacks_it`
  (test_sn_boundary_operator.py:226) already builds the exact `_BWithStubFace` non-adjointable fixture — a
  ~3-line `pytest.raises(MissingAdjoint)` extension. A guard whose loudness is unproven is one refactor from
  silent. **On every block-op carve that adds a fail-loud raise: grep the live tree for a negative gate hitting
  THAT method (not just the `is_adjointable` predicate), and check the fixture isn't already sitting there.**

- **NOTE — a rewired per-face test that mirrors the production formula becomes a PLUMBING pin, not a LAW
  oracle.** The rewired `test_sn_boundary_operator.py:203` assertion now computes `bc.apply_transpose(masked)`
  byte-identically to production ⟹ coextensive, catches only face-routing plumbing across `_CASES`. That's FINE
  because the defining-law oracle correctly moved to `test_psi_half_coupling.py` (`dense(Bᵀ)==dense(B).T` +
  `dense(Bᵀ)==0` on vacuum — structurally independent). When a fix rewrites a per-face test to mirror the fixed
  formula, verify the STRUCTURALLY-INDEPENDENT law oracle exists elsewhere (dense round-trip), and say which test
  is the oracle of record so a maintainer doesn't weaken it.

- **CREDIT (reinforce): the GMRES exact-breakdown anchor closed my d1 finding #2** — the permanent guard whose
  only prod trigger was the transient ψ½ pad now has a general singular-consistent (`diag(1,0)`) anchor +
  version-independent stub branch-pins + a nonzero-tail TEETH arm. This is the fuller-view/permanent-guard
  resolution done right: anchor the general invariant so the transitional caller is "first caller, not the reason."

## Phase C 4e — THE WALK UN-WEAVE (solve-leg extraction, reviewed 2026-07-11, verdict CONCERNS/ship-after-1-doc-fix)

The fused (L+C) 1-D walk's welded two-leg ψ½ ray-solve orchestration DELETED; the walk now
routes System B through the NAMED resolvent `RadialCharacteristicOperator.solve` (e2 @ `98fe2e36`),
after e1 (`63702e7`) retired the unified `RadialCharacteristicField` leaf family + unified space +
`from_unified`/`to_unified` bridge (3 modules deleted) and e1b (`ea7f919c`) reminted the freed name
onto the System-B composite. ZERO VIOLATIONS. The two durable, reusable rulings (the NEXT solve-leg
extraction — DSA/multiphysics — copies these):

- **The solve-orchestration extraction template (credit — the INVERSE of my L-001 catch).** The walk
  WRAPS the named resolvent, it does not re-implement: `RadialCharacteristicOperator.solve`
  (radial_characteristic.py:524,529) itself calls the single-sourced kernel
  `carlson_inward_sweep_from_source` — the exact two-leg march (inward + outward-on-reversed) the walk
  used to inline. Verify the extraction with THREE greps: (1) kernel refs in the walk module → 0
  (`carlson` 8→0 here); (2) the resolvent's `.solve` body actually calls the kernel (wrap, not fork);
  (3) the transpose mirrors (`.solve_transpose`→`carlson_inward_sweep_transpose`). The Mode-11 sentinel
  re-aims from the walk-namespace kernel onto the resolvent's `.solve` (CLASS-level wrap) + an S2
  tripwire asserts zero kernel-name text in the walk source (an inline march creeping back reds before
  any value gate goes blind). The augmented-cotangent transpose (`seed_cot.copy()` → per-level
  `+= psi_angle_bar` on the coupled-feed slot → ONE `solve_transpose` after the loop) is the clean
  adjoint of "System B is also a direct output": total cotangent = direct-output leg (`seed_cot`) +
  fused A_ABᵀ pullback; `.copy()` is defensive, not hidden state; strictly cleaner than the retired
  bare-`np.zeros_like`-buffer hand-assembly (Pattern 3 typed-intermediate).

- **The recurring CONCERN = operator-CONSTRUCTION triplication (institutional Pattern 1, coextensive-
  today).** Single-sourcing the solve ORCHESTRATION leaves the operator CONSTRUCTION `RadialChar…Operator(
  mesh, σ_t)` at N sites: driver `A_BB` (coupled_system.py:680, the loss-grid block, WITH `- B_b`) +
  walk forward (:4207) + walk transpose (:4775). Class + `.solve`/`.apply` single-sourced (the win);
  construction triplicated. Verify value-coextensiveness (here `_run`'s `sig_t` == `mat_xs.
  total_cross_section_field.values`) ⇒ CONCERN not VIOLATION. The driver's block is a DIFFERENT
  composition (`- B_b`) so it can't be threaded as-is; minting-inside-the-walk is the RIGHT shape when
  (i) the op is cheap (direct march, no factorization), (ii) function-local import is a real cycle-break
  (`operators.streaming → loss_representation`), (iii) threading touches the 6-signature `_run` chain.
  Remedy = reciprocal cross-ref + tracked trigger (walk→driver ref present; driver→walk absent), NOT
  premature unification (Pattern 6). Bug habitat: a future σ_t transform / new ctor arg lands on the
  residual-path block but not the production-solve mints. The 2nd route (loss-grid `A_BB`) may be LATENT
  (residual/DSA #2 consumer) — that keeps it a step-note, not a hot blocker.

- **THE do-now (recurring tell, do-now because in the carve's own file): the MODULE docstring rewritten
  clean, the OPERATOR CLASS docstring left STALE.** `radial_characteristic_field.py` correctly says
  "the unified layout + bridge are RETIRED", but `RadialCharacteristicOperator`'s class docstring
  (radial_characteristic.py:213-217) still asserts present-tense "bridge to the unified ψ½ layout" — the
  e1b diff proves the clause rode the mechanical `Composite→Field` rename UNCHANGED (a leftover e1 patched
  with a contradicting parenthetical but never rewrote). Same tell as step-4b (:176-179). On every un-weave:
  when the carve rewrites one docstring's scope section, grep the SAME operator's class + method docstrings
  for the OLD scope verb ("bridge to the unified …"). Bug habitat: a maintainer re-adds the retired bridge
  "to match the doc." Retirement is otherwise COMPLETE at code level (0 symbol refs to every deleted class/
  name/module); the only doc-orphans are historically-framed broken `:class:` xrefs to the deleted
  `RadialCharacteristicSourceSink` (2 surviving source-sink files) + the theory Key Facts `RadialChar…Flux`
  xref — de-linkify in the archivist pass (blast landed OUTSIDE the 52-file diff, per L4).
