# Operator-algebra carve — *inverse-as-operator* (grand-report alignment)

> **⚠ SUPERSEDED IN PART (2026-06-22).** Phase 1 (the domain/codomain two-param
> split) is now **Phase A of `.claude/plans/typed_carrier_grid_carve.md`**, where it
> finally has a real consumer (the non-endomorphic projection/reconstruction edges of
> the typed-carrier grid) — it is no longer speculative. The **foundation** (RC1
> `block_role` + RC3 `_as_boundary`, originally implied here) **landed as `272caa9`**.
> This plan's **Phases 2–5 (inverse-as-operator: `.inverse()` returns an operator,
> `.solve`→`.inverse().apply`, retire `CAP_SOLVE`/`MissingCapability`) REMAIN the
> follow-on**, built on the grid carve's Phase A. Do the grid carve FIRST, then resume
> here from Phase 2.

**A DETOUR from the #226 pyright burndown.** This is a deliberate architectural
carve, motivated by elegance — NOT a compliance pass. The #226 B4 "reds" are the
*symptom*; this plan fixes the *cause* and the reds clear as a consequence. There
is an explicit **EXIT back to the pyright burndown** in §8.

> **Cold-pickup contract:** a fresh post-compaction session should be able to run
> this end-to-end. Read §0 first, then §1–§5. This is **surgical / main-agent-direct
> + user-steered** work (`.claude/rules/delegation.md`): the main agent writes the
> code with the user steering; **do NOT dispatch `method-implementer`**.
> `test-architect` / `explorer` / `qa` / `elegance-enforcer` remain available.

---

## 0. Pickup-cold orientation (read first)

- **What this is:** move the SN operator algebra from *capabilities-as-strings +
  `solve`-as-gated-method* toward the grand report's **inverse-as-operator** design,
  where `solve` and `matvec` are two views of a dual operator pair and the
  sweep-vs-Krylov choice is *which inverse operator* `.inverse()` returns.
- **Why now:** the user rejected a "fix the 18 reds" framing and asked to use the
  type errors as a lever to improve the architecture so it "plays well with Python."
  The grand report (`.claude/plans/neutron_transport_grand_report_v3.md`) already
  specifies the target; current code is a partial realization.
- **Environment:** Host → use `.venv/bin/python`. Temp files in `$CLAUDE_JOB_DIR/tmp`.
  Run the **SessionStart protocol** (vv-principles, subagent-handoff-protocol,
  coding-elegance skills; `.claude/lessons.md`; nexus `session_briefing`;
  `docs/development.rst`).
- **Branch:** start a dedicated branch `refactor/operator-inverse-algebra` off
  `main`. (At time of writing main is `c22b6b5`; trust git, not this hash —
  `git log --oneline -1`.) Per-phase **ff-merge + push** when green (user-authorized
  cadence). `main` always green.
- **Oracle & constraints:** `npx --no-install pyright --outputjson orpheus/` CLI is
  the error oracle — the streamed `<new-diagnostics>` LAGS edits, never trust it over
  a CLI run. **NO `# type: ignore`** (hard user constraint — fix the real type).
  Canonical tests: `python -O -m pytest`. Re-tighten the ratchet after each cluster
  (`python -m tests._harness.pyright_ratchet --update`; gate `tests/test_pyright_ratchet.py`).
- **This carve is BEHAVIOR-PRESERVING.** Every numerical output (matvec, `.H`, solve,
  keff, MMS) stays **bit-identical** (or principled-equivalent only where FP-reassoc
  is forced). Any non-bit change in the matvec≡sweep twins is a BUG, not progress.

### Reading list for the cold pickup
1. This plan (§1–§8).
2. `.claude/plans/issue_226_b4_operator_generics_verification.md` — the test-architect
   verification net (contract-pinning tests, risk map RC1–RC4, 3 gaps, 5 mutation recipes).
3. `.claude/agent-memory/explorer/issue_226_operator_generics_map.md` — the structural map.
4. Grand report **§1** (operator algebra `K = A_loss.inverse() @ F`, ~lines 55–110),
   **§3.1** (capability Protocols, line 157), **§3.3** (`LinearOperator(ABC)` taxonomy,
   lines 275–332), **§22** (resolvent, ~4377–4415).

---

## 0.5 Report-aligned decisions (LOCKED 2026-06-29, user-steered)

The grand report's operator-algebra sections (§1, §3.1, §5.7, §6.3, §21–22, §30–31, §38)
were read to guide this carve; it validates the direction wholesale and locks these:

1. **Inverse-operator naming = HYBRID.** `.inverse()` returns the report's math-canonical
   concretes — **`SweepOperator`** (structured/triangular direct inverse of `(L+C)`; the
   sweep IS the inverse application, §5.7), **`GreenOperator`** (general resolvent,
   Krylov-applied, §5.7/§22), future **`ResolventOperator`** (parametrized `(sT+A).inverse()`,
   α-work). The **greppable family marker is a `SupportsInverse` Protocol** (`def inverse(self)
   -> LinearOperator`), matching the report's §3.1 `Supports*` family (`SupportsApply`,
   `SupportsAdjoint`); `grep SupportsInverse` finds every inverse-producer. (Name `SupportsInverse`
   preferred over `Invertible` to avoid collision with the concrete `InvertibleOperator` = `(L+C)`,
   and to sit in the `Supports*` family. SUPERSEDES the plan body's `SweepInverseOperator`/
   `KrylovInverseOperator`.)
2. **Capability retirement is SOLVE-AXIS-ONLY, Protocol-based (§3.1 + test-architect
   CORRECTION).** Phase 4 retires the **`CAP_SOLVE` axis** → the structural `SupportsInverse`
   Protocol; invertibility = `isinstance(op, SupportsInverse)` / having `.inverse()`. **KEEP
   `CAP_APPLY_TRANSPOSE` as a runtime frozenset SoT** — it gates the `OperatorSum/Product.
   apply_transpose` "propagates iff both" law (operator.py:762-763/841-842) that the LANDED
   #276 S† rests on; making `apply_transpose` universal-on-base would let a half-transposable
   sum forge a wrong adjoint (breaks `test_isotropic_scattering.py:219`, `test_capability_
   survival.py:139`, the S† bit-id twin). Delete the `solve` AND `apply_transpose` `# type:
   ignore`s by **declaring** both methods on the Protocol (optional-presence, RC4 spelling)
   while the transpose RUNTIME gate stays the frozenset. The honest solve residue (value-
   dependent singular: `DiagonalOperator`/`MultiplicationOperator` zero-coefficient,
   `min|f|=0`) → an **eager RUNTIME RAISE inside `.inverse()`** (raise, not bare-assert, for
   `-O`; positive+negative gated), NOT a surviving solve frozenset. (§6.3 resolved; static
   `ScaledOperator(0.0)` is already a ctor `ValueError`.)
   **`CAP_SOLVE` consumer surface (Phase 3, wider than first mapped):** SourceIteration +
   `_GaussSeidelResolvent`/`_MomentWindowedResolvent` (the real `L` on curvilinear paths) +
   KrylovAcceleration preconditioner-selection (`_has(L, CAP_SOLVE)`, iteration.py:712) +
   DiagonalOperator/MultiplicationOperator/every composer. Full table: verification spec §1a.
3. **`solve(A, b)` free-fn KEPT as call-site sugar** ≡ `A.inverse() @ b` (§6.4 fork → keep;
   §21.2 power iteration uses it alongside `loss.inverse() @ F`).
4. **`.H` vs `.T` (§6.3) — already landed by the adjoint campaign (now on main).** `.T` =
   Euclidean representation transpose (`loss_action_transpose`); `.H` = Hilbert adjoint
   (`G⁻¹AᵀG` wrapper). The coherence identity **`A.H.inverse() == A.inverse().H`**
   (adjoint-of-inverse = inverse-of-adjoint) IS #280's transpose-solve — the inverse carve
   SUBSUMES A3, no `CAP_SOLVE_TRANSPOSE`. (test-architect gate #4 pins it.)
5. **Operator ≠ context (§38).** `SweepOperator` is the *semantic* inverse; the sweep
   coefficient cache is its application *context*, not part of the semantic object.
6. **§6.5 fork → `SupportsInverse` Protocol (structural) + concrete `SweepOperator`/
   `GreenOperator`.** The concrete `InvertibleOperator` (= `(L+C)`, `OperatorSum["FullField"]`)
   STAYS concrete AND satisfies `SupportsInverse` structurally; `.inverse()` on it returns a
   `SweepOperator`. Gate the no-inverse residue STRUCTURALLY (`hasattr(op,"inverse")`, not
   `isinstance(InvertibleOperator)`) so it survives either shape.
7. **Phase 2 ships `SweepOperator` ONLY; `GreenOperator` DEFERRED (§4c new fork, Pattern 6).**
   `(L+C).inverse() → SweepOperator` has the real consumer (SourceIteration's inner solve).
   The general-sum inverse `A_loss.inverse() = (L+C−S).inverse()` (= `GreenOperator`, Krylov-
   on-sum) has NO current consumer (the k-eigenvalue path is power-iteration-with-SI-inner,
   not a literal `A_loss.inverse()` operator) → defer until a consumer (α-resolvent / Krylov
   eigenvalue). Do NOT ship `OperatorSum.inverse()` now.
8. **`A.H.inverse()` (the transpose-solve) is #280, NOT this carve.** The carve makes it
   CONSTRUCTIBLE (the `.H`/`.inverse()` algebra composes) but does NOT implement the reverse-DAG
   `sweep_transpose`; `(L+C).H` correctly has no `.inverse()` until #280 (redesigned onto
   `A.H.inverse()`). The adjoint-coherence gate (spec §5) rides `xfail(strict=False, "#280")`
   until then — the algebra is proven coherent, the transpose-solve impl stays paused per the
   "Full carve, redesign #280" scope.

> Where the plan body below says `SweepInverseOperator`/`KrylovInverseOperator`, read
> `SweepOperator`/`GreenOperator`; where it says "retire `CAP_SOLVE` via isinstance/Protocol",
> the Protocol is `SupportsInverse`.

**`SweepOperator` location:** its own file `orpheus/sn/operators/sweep_operator.py` (user, 2026-06-29)
— out of `streaming.py` until a generic `numerics` home earns itself.

---

## 0.6 Unified inverse + adjoint materialization — NO twin-path (user directive 2026-06-29)

**The hazard the user blocked.** test-architect's "retire `CAP_SOLVE`, KEEP `CAP_APPLY_TRANSPOSE`"
is a correct LOCAL fix (don't break the landed S† "propagate iff both") but, left as-is, splits the
architecture into a **twin-path**: invertibility queried via the `SupportsInverse` Protocol / `.inverse()`,
adjointability via the `CAP_APPLY_TRANSPOSE` string. Two mechanisms for the identical structural
question — the stringly-typed anti-pattern surviving on one axis while the other modernizes (Cardinal
Rule 2). FORBIDDEN: the carve must materialize the adjoint via the SAME direction as the inverse.

**The realization that dissolves it: the "propagate iff both" LAW already lives in the method bodies,
NOT the frozenset.** `OperatorSum.apply_transpose(x)` IS `a.apply_transpose(x) + b.apply_transpose(x)`
(operator.py:787) — the adjoint of a sum is the sum of adjoints, by construction; a summand without a
transpose makes the expression raise. Identically `(AB)⁻¹ = B⁻¹A⁻¹` IS `OperatorProduct.solve`
(operator.py:858). Both axes already propagate structurally. The `frozenset[str]` is only a redundant,
stringly-typed, axis-conflating ADVERTISEMENT ("will it work?" pre-check). It retires for BOTH axes via
ONE mechanism:

| Concern | Solve axis | Adjoint axis |
|---|---|---|
| The operator | `.inverse() -> SweepOperator/GreenOperator` | `.H -> _AdjointOperator` |
| Static contract | `SupportsInverse(Protocol)` | `SupportsAdjoint(Protocol)` (report §3.1) |
| Composition law | `(AB).inverse()=B⁻¹@A⁻¹` (method body) | `(A+B).H=A.H+B.H`, `(AB).H=B.H@A.H` (method body, already) |
| "iff both" QUERY | derived `@property is_invertible` (recursive `a.is_invertible and b.is_invertible`) | derived `@property is_adjointable` (recursive) |
| Undefined leaf | eager raise `NotInvertible` at `.inverse()` | raise `MissingAdjoint` at `.H` |

`apply` is universal (the base method — no CAP). So the **entire `capabilities: frozenset[str]` retires**:
CAP_APPLY → the base `apply`; CAP_SOLVE → `SupportsInverse` + `is_invertible`; CAP_APPLY_TRANSPOSE →
`SupportsAdjoint` + `is_adjointable`. Per-axis, typed, structurally COMPUTED (not a stored set that can
drift), each paired with its operator-returning method. "Propagate iff both" = the derived predicate's
recursive `and`, NOT a frozenset intersection. EAFP (try `.H`/`.inverse()`, catch the typed exception) is
the alternative for consumers that don't pre-query; the genuine pre-query consumer is
KrylovAcceleration's preconditioner-selection (`_has(L, CAP_SOLVE)` → `L.is_invertible`).

**Safe for the landed #276 S†:** the reciprocity LAW is the method body (`a.H + b.H`), untouched; only the
ADVERTISEMENT moves from `CAP_APPLY_TRANSPOSE ∈ caps` to `is_adjointable`. Gates: the S† bit-id twin +
`test_g_adjoint_reciprocity` (+ wrong-metric L11 control) stay green, AND a new gate pins
`is_adjointable == (old CAP_APPLY_TRANSPOSE membership)` for every operator (advertisement faithful) with
an M-mutation (drop a summand's adjoint → composite `is_adjointable` False, the sum's `.H` raises).

**Sub-decision (settle with test-architect / elegance-enforcer):** the runtime QUERY form — a derived
`@property is_invertible`/`is_adjointable` (instance-accurate on composites, my lean) vs. `isinstance(op,
SupportsInverse/SupportsAdjoint)` + EAFP (Protocol presence is class-uniform, so not instance-accurate on
composites — insufficient alone for the pre-query consumer). The Protocols carry the STATIC contract
either way; the derived predicate carries the RUNTIME instance truth.

---

## 0.7 IMPLEMENTATION STATUS — cold-pickup resume point (2026-06-29)

> **⚠ SUPERSEDED IN PART (2026-07-01):** the P3+ phases below are RE-SCOPED by the
> taxonomy design `.claude/plans/operator_machinery_taxonomy.md` (§12 phase order;
> §13 name-earning invariants; §14 the sphere-seed #282 + CP #283 rulings). The
> committed P2a/P2b/P2c foundation stands unchanged. Resume from the TAXONOMY plan,
> not from the P3 block below — the duck-typed resolvents dissolve into configured
> SweepOperators, and SourceIteration consumes an already-inverse operator via
> `.apply` (no per-consumer `.solve→.inverse().apply` rollout).

> A fresh post-compaction session resumes HERE. Branch `refactor/inverse-as-operator` off
> `main` @ `f903686` (**adjoint-inclusive** — the #276 K_iso / F† / S† campaign is now ON
> main; trust git: `git log --oneline -5`). Surgical, **main-agent-direct, user-steered, NO
> method-implementer**. Per-phase ff-commit WHEN THE USER ASKS; `main` always green. Oracle
> `npx pyright orpheus/` = **410** (count-neutral target; the carve drives it DOWN at P4).
> Canonical gate `.venv/bin/python -O -m pytest <paths> -m "not slow" -q -p no:xdist -p
> no:cacheprovider --timeout=300`. NO `# type: ignore`; NEVER `git checkout/restore/stash` on
> uncommitted files.

### DONE — Phase 2a, the predicate foundation (committed, verified)
- `32314cc` chore: ratchet re-baseline 420→410 (the adjoint merge's −10, never re-tightened).
- `c1126c5` feat: **Phase 2a** in `orpheus/numerics/operator.py` — ADDITIVE, count-neutral,
  faithfulness-gated. Frozenset path UNTOUCHED, nothing retired.
  - `SupportsInverse`/`SupportsAdjoint` Protocols (NON-`runtime_checkable` — static contract
    only; in `__all__`).
  - base `LinearOperator.is_invertible`/`is_adjointable` `@property` (default `False`, after
    `codomain`).
  - composer recursion: `OperatorSum` (is_adjointable=a∧b; is_invertible inherits False),
    `OperatorProduct` (both=a∧b), `ScaledOperator` (pass-through, α≠0 by ctor), `IdentityOperator`
    (both True), `_AdjointOperator` (inherits False/False — correct, no CAP_SOLVE/CAP_APPLY_TRANSPOSE).
  - leaves: `Zero`/`IncomingOrdinateMaskTensor`/`PeriodicWrap` (adjointable-only),
    `Permutation` (both), `TensorProduct`/`SumOfTensorProducts` (recurse `self.ops`/`self.summands`),
    `Diagonal` (is_invertible=`np.all(coeff!=0)`, adjointable=True), `RankOne` (NONE — {apply} only).
  - gate `tests/numerics/test_operator_capability_predicates.py`: `is_X≡CAP_X∈caps` for every op
    + asymmetry fixtures (`ZeroOperator` adj-not-inv, `_ApplyOnly` neither) + recursive pins. 12 green.

### DONE — Phase 2b, the transport + SN leaf predicates (committed `d36507b`, verified)
ADDITIVE; frozenset path UNTOUCHED. CLI pyright **410 → 409** (a −1 WIN, not just neutral),
faithfulness-gated, all canaries green. The comprehensive `grep -rn capabilities orpheus/`
surfaced advertisers the original P2b note MISSED — the real scope:
- **transport energy ops** (all `{APPLY,APPLY_TRANSPOSE}` → `is_adjointable=True`, is_invertible
  inherits False): `scattering.py` `LegendreMomentScattering`/`N2NMomentOperator`/`ScatteringOperator`,
  `fission.py` `FissionOperator`, `isotropic_scattering.py` `IsotropicScattering`/`IsotropicN2N`.
- **`multiplication_operator.py` `MultiplicationOperator`** → **DELEGATE** to `self.engine`
  (`is_invertible=engine.is_invertible` value-dependent min|f|>0; `is_adjointable=engine.is_adjointable`).
- **numerics faces**: `projection.py` `AnalysisOperator` ABC → `is_adjointable=True` (`_FrameAnalysis`
  inherits); `frame.py` `_FrameReconstruction` → `is_adjointable=True` (OVERRIDE — its ABC
  `ReconstructionOperator` advertises only `{APPLY}`, so the predicate follows the FACE's caps, not the ABC).
- **sn**: `streaming.py` `StreamingOperator` → `is_adjointable=True`; `InvertibleOperator` →
  `is_invertible=True` (is_adjointable inherits OperatorSum a∧b); `operators/boundary.py`
  `SNBoundaryOperator` → `is_adjointable = bool(laws) and all(law.is_adjointable for law in laws)`
  (mirrors the caps intersection rule).
- **BUG FOUND + FIXED** (the comprehensive grep's payoff): `geometry/boundary/_bound_compat.py`
  `_BoundBoundaryOperator` delegated `capabilities` to `self.inner` but NOT the predicates → it
  inherited base `False`, which would have made the `SNBoundaryOperator` aggregator UNFAITHFUL
  (caps via inner = True, is_adjointable = base False). Added `is_invertible`/`is_adjointable`
  delegation. Also typed `boundary.py:_face_laws -> dict[str, LinearOperator]` (was `object`) — this
  is the −1 pyright win (resolved the pre-existing `object`-attribute-access at the `.capabilities` read).
- **gate** (single-source helper `tests/_harness/predicates.py::assert_capability_faithful`, explicit
  `raise` so it fires under `-O` — a helper module is NOT pytest-assert-rewritten, Mode-8 trap):
  retrofit `tests/numerics/test_operator_capability_predicates.py` to import it; ADD frame-face
  faithfulness to `tests/numerics/test_frame.py`; ADD `TestPredicateFaithfulness` (transport energy +
  SN streaming/boundary + the `_BoundBoundaryOperator` wrappers + the value-dependent zero-coeff
  asymmetry) to `tests/sn/operators/test_capability_survival.py`. 42 + 132 green.
- **DEFERRED to P3 (correct phase boundary):** the `sn/solver.py` duck-typed resolvents
  (`_GaussSeidelResolvent`/`_MomentWindowedResolvent`) advertise `CAP_SOLVE` but are NOT
  `LinearOperator` subclasses (no base predicate to inherit). Their `is_invertible` declaration lands
  WITH their consumer migration (`SourceIteration` guard `_has(L,CAP_SOLVE)`→`L.is_invertible`) and its
  pinning gate (`test_keff_curvilinear`), so the predicate + consumer + gate land together at P3.

### DONE — Phase 2c, the SweepOperator (the inverse OPERATOR) (uncommitted, verified)
CLI pyright **409** (count-neutral), keystone gate green + **teeth proven** (M-EQ in-process).
- NEW FILE `orpheus/sn/operators/sweep_operator.py`: `SweepOperator(LinearOperator["FullField"])`
  (plain wrapper, like `_GaussSeidelResolvent` — avoids dataclass hash on the mutable inner).
  `apply(rhs, *, initial_guess=None)` → `self.inner.solve(rhs, initial_guess=...)` (BIT-IDENTICAL);
  `domain=inner.codomain`/`codomain=inner.domain` (the inverse SWAPS; equal here — endomorphic);
  `capabilities={CAP_APPLY}`; `.H`/`.inverse()`/`is_*` DEFER → inherit base False (adjoint-inverse =
  #280). Docstring carries the operator-vs-context (§38) + sweep-vs-Krylov (§5.7) rationale.
- `InvertibleOperator.inverse()` → `return SweepOperator(self)` (runtime LATE import to break the
  SweepOperator↔InvertibleOperator cycle) + a `TYPE_CHECKING` import of `SweepOperator` so the
  `-> "SweepOperator"` annotation resolves (the +1 pyright the late-import-alone caused). `__init__`
  docstring bullet added for the new operator.
- KEYSTONE gate `tests/sn/operators/test_inverse_operator_equivalence.py` (`foundation`,
  `verifies("inverse-as-operator")`): `inverse().apply(b) == solve(b)` `assert_array_equal` on
  bulk+boundary, slab+sphere+cyl ≥2G/vacuum; the curvilinear seed on sphere+cyl PAIRED with a Mode-11
  in-process `solve`-spy (the EXACT seed object threads through); + a SweepOperator-structure pin.
  **M-EQ teeth**: monkeypatching `apply`→forward `inner.apply` reddens the value gate (proven). 7 green.
- **DEFERRED to P5** (composition return-type algebra, spec Gate 2.2): `OperatorProduct.inverse()`
  (`(AB)⁻¹=B⁻¹A⁻¹`) + `ScaledOperator.inverse()` (`(αL)⁻¹=(1/α)L⁻¹`). No Phase-3 consumer calls
  `.solve` on those composites (only InvertibleOperator + the resolvents), so they ride P5.

**P3 — consumer migration** (`.solve`→`.inverse().apply`; `_has(CAP_X)`→`.is_X`; full §1a table in
the spec):
- `numerics/iteration.py`: `:435` SourceIteration guard `_has(L,CAP_SOLVE)`→`L.is_invertible`;
  `:712` KrylovAccel preconditioner `elif _has(L,CAP_SOLVE)`→`L.is_invertible` (the genuine pre-query
  consumer); `:467/468` `L.solve(...)`→`L.inverse().apply(...)` + DELETE the `_solve_accepts_seed`
  inspect probe.
- `numerics/operator.py` composer transpose reads {646,669,762,841,900,1333,1356,1437,1448,1601,1819}
  `_has(op,CAP_APPLY_TRANSPOSE)`→`op.is_adjointable`; `sn/operators/boundary.py:141`.
- `sn/solver.py` `_GaussSeidelResolvent`/`_MomentWindowedResolvent` (the real curvilinear `L`) +
  `base_resolvent.solve(:2347)`.

**P4 — FULL frozenset retirement, BOTH axes (gated — CHECKPOINT THE USER FIRST)**:
- delete `capabilities` frozenset (all leaves+composers), both `_has` (operator.py:594,
  iteration.py:154), `MissingCapability` solve+transpose, `.solve()` methods, the ~10
  `# type: ignore[attr-defined]`. CAP_APPLY→base `apply`; KEEP a one-line `hasattr(op,"apply")`
  eager composition guard (§4e, test-architect rec).
- eager `.H` raise: `_AdjointOperator.apply` MissingCapability (operator.py:669) → raise eagerly at
  `.H`/`adjoint()`; :669 dead. BEHAVIOR CHANGE (lazy→eager) — re-grep tests for lazy-`.H`-raise.
- **THE TRAP (test-architect): literal-string test preconditions** `"apply_transpose" not in
  A.capabilities` at `test_g_adjoint_reciprocity.py:210` + `test_removal_form_matvec_sweep.py:308`
  (the very S† canaries) + `test_scattering_operator:176`/`test_fission_operator:79`/
  `test_bc_universal_invariants:419`/`test_boundary:439` → rewire to `not A.is_adjointable`. **"S†
  stays green" REQUIRES this rewire.** (`tests/sn/solve/conftest.py:4 cap("solve")` = vv marker — LEAVE.)
- the S† reciprocity LAW is the method body (`a.H+b.H`), UNTOUCHED; gates stay green: the S† bit-id
  twin + `test_g_adjoint_reciprocity` + wrong-metric L11 control; M-mutations per spec §7.

**P5 — docs + #280 redesign**: operator-algebra theory + SN dev-history changelog
(`discrete_ordinates.rst`: inverse-as-operator, why CAP retired, sweep-vs-Krylov as inverse kinds);
update GH **#280** to drop `CAP_SOLVE_TRANSPOSE` → `A.H.inverse()`; **EXIT** to the pyright burndown
(`pyright_signal_cleanup.md`, mark B4 DONE-via-carve).

### Pointers
- **Verification spec (gate authority):** `.claude/plans/issue_226_inverse_operator_verification.md`
  (§2 keystone equivalence, §2.3 faithfulness, §3b eager-raise, §4 retirement delete/rewire/keep map,
  §5 adjoint-coherence `A.H.inverse()` = `xfail("#280")`, §7 M-mutations, §8 zero-scaling residue).
  Phase-0/1 (landed): `issue_226_b4_operator_generics_verification.md`.
- **test-architect instances** (resume via SendMessage if needed): `a484bc0c04f169b42` (solve-axis net),
  `a8e56b7736037102b` (both-axes extension).
- **Tasks:** #132 (P2, in_progress), #133 (P3), #134 (P4 full retirement), #135 (P5).
- **Decisions locked** in §0.5 (naming/Protocol/solve-sugar) + §0.6 (the unified no-twin-path mechanism).

---

## 1. The architectural thesis (the WHY)

The current operator algebra reimplements **at runtime** three things the type system
expresses natively, and the grand report fixes all three with one idea — *the inverse
is an operator*:

1. **Capabilities as a `frozenset[str]`** (`CAP_APPLY/CAP_SOLVE/CAP_APPLY_TRANSPOSE`)
   + `MissingCapability` is stringly-typed dispatch (a `coding-elegance` anti-pattern),
   and it is *why* the `# type: ignore`s exist (`solve`/`apply_transpose` aren't on the
   base type).
2. **`CAP_SOLVE` conflates "is invertible" with "has a *direct* (sweep) solver."**
   `(L+C−S)` is mathematically invertible (Krylov) but has no sweep, so it is denied
   `CAP_SOLVE` — the "sum drops solve" law is an artifact of this conflation.
3. **Even the *shape* of `solve` is introspected** (`inspect.signature(...).parameters`
   in `SourceIteration` to detect seedable vs plain solve).

**Grand-report resolution (the three load-bearing pieces):**

- **The inverse is an operator** (§1, line 76): `A_loss = L + C − S; K = A_loss.inverse() @ F`.
  No `solve` method, no `CAP_SOLVE`. "Solve" = **apply the inverse operator**. So
  `A @ psi` (forward view) and `A.inverse() @ q` (inverse view) are the *two views of
  one operator* — duals like `A` and `A.H`. (`InvertibleOperator` becomes "a
  `LinearOperator` that knows its inverse," forward view `apply`, inverse view
  `inverse().apply`.)
- **Sweep vs Krylov = which inverse operator `.inverse()` returns** (§3.3: `SweepOperator`,
  `GreenOperator`, `ResolventOperator` are all `LinearOperator` subclasses). Structured/
  triangular → returns a `SweepInverseOperator` (direct forward substitution); general →
  returns a `KrylovInverseOperator` (GMRES on `A.apply`). Same interface (`apply`),
  different algorithm. This dissolves the "sum drops solve" conflation: every invertible
  operator can produce an inverse operator; "has a sweep" is just *which kind*.
- **Domain/codomain split fixes variance + rehomes non-endomorphic operators** (§3.1,
  line 170: `SupportsApply(Protocol[X, Y]): def apply(self, x: X) -> Y`). Separate
  domain/codomain type params (domain **contravariant**, codomain **covariant**) instead
  of the single invariant `apply(x: V) -> V`. This clears the `TensorProductOperator →
  LinearOperator[Unknown]` invariance rejection AND makes `ProjectionOperator`
  (rich→coeffs) / `ReconstructionOperator` (coeffs→rich) **first-class** instead of the
  quarantined outliers they are today.

**The honest residue (one real type boundary):** singular maps (projections, masks,
zero) genuinely have no inverse — they are `LinearOperator` but NOT `InvertibleOperator`
(no `.inverse()`). Plus the value-dependent `ScaledOperator(0·op)` singularity (a runtime
fact, not a type). This is the only place a runtime check legitimately survives.

---

## 2. Current state (the FROM) — exact anchors

`orpheus/numerics/operator.py`:
- `class MissingCapability(TypeError)` L229.
- `class LinearOperator(Protocol[V])` L252 (declares `capabilities: frozenset[str]` ≈L289, `apply`).
- `class LinearOperatorMixin(Generic[V])` L344 (`capabilities` instance attr ≈L378; `block_role: ClassVar[Optional[BlockRole]]` ≈L398).
- `_AdjointOperator(LinearOperatorMixin[V])` L591; `apply_transpose` L663; sets `self.block_role` ≈L623.
- `OperatorSum(LinearOperatorMixin[V])` L671; `__init__` L694; caps derivation L727–731 ("solve does NOT propagate"); `block_role = _join_block_roles(...)` L736; `apply` L750; `apply_transpose` L753 **with `# type: ignore[attr-defined]` L754**.
- `OperatorProduct(LinearOperatorMixin[V])` L757; `solve` L821 (product propagates solve).
- `ScaledOperator(LinearOperatorMixin[V])` L830; `solve` L880; sets `self.block_role` ≈L867.
- Leaves with class-attr `capabilities` + per-leaf `solve`/`apply_transpose`: ≈L896/903/906/956/968/1026/1059/1062/1120/1161/1201/1210/1310/1323/1402/1636/1640.
- `class TensorProductOperator(LinearOperatorMixin)` L1214 — **BARE, no `[V]` → `Unknown`** (the variance-break source). `IncomingSourceOperator` similar.

`orpheus/numerics/iteration.py` (`SourceIteration`):
- gains `CAP_APPLY` guard L428–434; `L` `CAP_SOLVE` guard L435–443; `self.L = L` L445;
  `_solve_accepts_seed = "initial_guess" in inspect.signature(self.L.solve).parameters` L457–462;
  `_solve_with_seed` L464–468. **`.solve` reds: L459/467/468.**
- (Adjacent, OUT of this carve's scope — leave for B5/other: L228 scipy LinearOperator call, L527/568 union, L825 `tol`, L1029 return.)

`orpheus/sn/operator.py`:
- `class InvertibleOperator(OperatorSum["FullField"])` **L631** — "the one operator":
  the `(L+C)` sum with the sweep `solve`. This is what the end state generalizes.
  Parametrized with the `["FullField"]` **string proxy** (layering: `numerics` can't
  import `transport.TimedFullField`).
- `L922` singledispatch `apply` arms (scattering/fission union).

`orpheus/sn/boundary_realizer.py`:
- `_as_boundary` does `op.block_role = BlockRole.BOUNDARY` → `:125` "Cannot assign to
  attribute block_role for class LinearOperator[Unknown]". `:176/206/207/223/224/240/252`
  `TensorProductOperator → LinearOperator[Unknown]` (reportArgumentType); `:260`
  `IncomingSourceOperator → LinearOperator[Unknown]` (reportReturnType).

`orpheus/sn/angular_operator.py` L102/270, `orpheus/sn/boundary_operator.py` L137 —
`capabilities` ClassVar-overrides (`reportIncompatibleVariableOverride`).

`orpheus/numerics/vector.py:151`: `V = TypeVar("V", bound=Vector)` (the bound already
landed — do NOT re-recommend an unbounded V).

**B4 cluster size:** 18 reds in scope (per the verification plan); ~10 adjacent reds in
these files are B5/other and MUST NOT be touched or regressed.

---

## 3. Target state (the TO) — the grand-report algebra

```python
# numerics/operator.py — the base
class LinearOperator(Protocol[D_contra, C_co]):     # domain contravariant, codomain covariant
    def apply(self, x: D_contra) -> C_co: ...
    @property
    def H(self) -> "LinearOperator[C_co, D_contra]": ...   # adjoint operator (swapped)

Endomorphism = LinearOperator[V, V]                 # readability alias for the common case

class InvertibleOperator(LinearOperator[D, C], Protocol):
    def inverse(self) -> LinearOperator[C, D]: ...   # the inverse OPERATOR (swapped)

# solve(A, b) free fn == A.inverse() @ b  (grand report line 4631 / line 76)
```

- `apply` (matvec) and `.H` (adjoint) are **universal** on `LinearOperator`.
- `.inverse()` is the distinguishing method of `InvertibleOperator`; it returns an
  operator (`SweepInverseOperator` for structured/triangular, `KrylovInverseOperator`
  for general). `solve` the *method* is **retired**; `A.inverse() @ b` (or `solve(A,b)`)
  is the inverse view.
- `capabilities: frozenset[str]` + `CAP_SOLVE` + `MissingCapability` for the solve axis
  are **retired** → invertibility is `isinstance(op, InvertibleOperator)` / the
  `.inverse()`-bearing subtype. Apply/transpose universal ⇒ the `# type: ignore`s on
  `apply_transpose` delete (not suppress). Keep a minimal runtime check ONLY for the
  value-dependent singular edge (zero-scaling) if it genuinely exists — investigate.
- Composition algebra encoded in `@`/`+`/`.inverse()` **return types**
  (`InvertibleOperator @ InvertibleOperator -> InvertibleOperator`; a sum is invertible
  via Krylov so `.inverse()` exists there too).
- Non-endomorphic `ProjectionOperator`/`ReconstructionOperator` become first-class
  `LinearOperator[Rich, Coeffs]` / `[Coeffs, Rich]`.

The acceptance test (coding-elegance master standard): the algebra reads like the math —
`K = A_loss.inverse() @ F`, `psi = (L + C).inverse() @ q`, `psi_dagger = A_loss.H.inverse() @ R`.

---

## 4. Verification net (must stay green at every phase)

Full inventory in `.claude/plans/issue_226_b4_operator_generics_verification.md`. The
load-bearing gates (PURE re-typing ⇒ **zero** tests reclassify to principled-equivalent):

- **Bit-identity ground:** `tests/sn/operators/test_removal_form_matvec_sweep.py::test_invertible_apply_is_M_of_C_sigma_bit_identical` (`assert_array_equal`); `TensorProduct` `test_apply_against_einsum_reference`.
- **Capability closure:** `tests/numerics/test_operator.py` (`test_sum_solve_does_not_propagate`, `test_product_solve_{propagates,drops}`, `test_scaled_preserves_all_capabilities`, 3 `MissingCapability` negatives); `tests/numerics/test_tensor_product_operator.py::TestTensorProductCapabilities`.
- **block_role join/partition:** `tests/sn/operators/test_operator_block_role.py` (esp. `test_sum_join_is_union_of_blocks`: `BULK⊔BOUNDARY=FULL`; `test_realized_bc_advertises_boundary_operator` — the `_as_boundary` stamp).
- **`.H` G-adjoint:** `test_operator.py::test_adjoint_{euclidean,weighted}_identity`; `tests/sn/operators/test_g_adjoint_reciprocity.py` (incl. wrong-metric negative control).
- **End-to-end:** `tests/sn/eigenvalue/test_keff_curvilinear.py` (L1 k_inf, L2 2G-heterogeneous), MMS slab/2d, `test_prescribed_inflow_consistency`.

**3 gaps to add as tests BEFORE the carve (Phase 0):**
- GAP-1: `TensorProductOperator`/`IncomingSourceOperator`/`SumOfTensorProductsOperator` default `block_role is None` (→ `test_operator_protocols.py::TestGenericOperatorsAreUnclassified`).
- GAP-2: `IncomingSourceOperator` per-leaf `capabilities ==` pin + a `(L+C−B).capabilities` composite-survival pin for the property-based `SNBoundaryOperator`.
- GAP-3: an apply-only leaf stays non-solvable after the base gains a `solve`/`inverse` surface (`CAP_SOLVE ∉ caps` AND `SourceIteration(apply_only)` raises) — the analogue under the new design is "an apply-only / singular operator has no `.inverse()`".

**Mutation recipes M1–M5** (in the verification plan): inject a `block_role` flip / drop a
leaf capability / swap an `OperatorSum` operand etc.; each MUST turn its gate RED under
`-O` (bare `assert` is stripped — gates use `np.testing.*`/`pytest.fail`).

---

## 5. Migration order (the HOW) — incremental, each phase behavior-preserving + gated

> Sequence chosen so the *typing foundation* lands first (low-risk, clears most reds),
> then the inverse-operator is introduced *alongside* the old path and proven equivalent
> before anything is retired. ff-merge+push each phase when green.

### Phase 0 — Safety net (additive; autonomous-safe)
- Confirm every §4 twin green at branch HEAD. Add GAP-1/2/3 tests. Run M1–M5, confirm RED.
- **Gate:** net green; mutations RED. This is the precondition; do not start Phase 1 without it.

### Phase 1 — Domain/codomain typevar split (the report's `[X,Y]` shape)
- `LinearOperator(Protocol[V])` → `LinearOperator(Protocol[D_contra, C_co])`,
  `apply(self, x: D) -> C`; `LinearOperatorMixin(Generic[D, C])`. Add `Endomorphism[V] =
  LinearOperator[V, V]` alias for readability.
- Endomorphic operators parametrize `[V, V]`; bare `TensorProductOperator` →
  `LinearOperatorMixin[np.ndarray, np.ndarray]` (or `Endomorphism[np.ndarray]`);
  Projection/Reconstruction get their true `[Rich, Coeffs]`/`[Coeffs, Rich]`.
- `.H` typed `-> LinearOperator[C, D]` (swapped).
- **Clears:** the `TensorProductOperator → Unknown` variance reds (boundary_realizer 8) +
  `:125`. **Zero runtime.** **Gate:** all twins green; CLI shows those reds gone.
- DECISION (see §6): whether to rehome Projection/Reconstruction now (bigger) or alias
  them endomorphically and defer.

### Phase 2 — Introduce `.inverse()` as a first-class operator (coexists with `.solve()`)
- Add inverse-operator types `SweepInverseOperator` (wraps the existing sweep — for
  `InvertibleOperator`, `sn/operator.py:631`) and `KrylovInverseOperator` (GMRES on
  `A.apply`), both `LinearOperator[C, D]` whose `.apply(b, *, initial_guess=None)` returns
  `A⁻¹b`.
- Add `inverse()` to `InvertibleOperator` (returns `SweepInverseOperator` wrapping the
  current sweep — **bit-identical** to today's `.solve`). Add `inverse()` to the general
  invertible composites (`OperatorSum`/`OperatorProduct`) returning `KrylovInverseOperator`
  (the "sum is Krylov-invertible" insight).
- **Gate:** NEW test `A.inverse().apply(b) == A.solve(b)` **bit-identical** for every op
  that has both; existing twins green. Nothing removed yet.

### Phase 3 — Migrate consumers `.solve()` → `.inverse().apply()`
- `SourceIteration`: annotate `L: InvertibleOperator`; `self.L.solve(...)` →
  `self._L_inv = L.inverse(); self._L_inv.apply(rhs, initial_guess=...)`. The
  `_solve_accepts_seed` **inspect probe deletes** (the inverse operator's `apply` has
  `initial_guess` canonically). The `CAP_SOLVE` guard → `isinstance(L, InvertibleOperator)`
  / the annotation.
- Eigenvalue/Krylov consumers (`KrylovAccelerator`, `loss.inverse() @ production` — already
  the report's target form, grand report line 4159). Migrate one consumer at a time.
- **Clears:** iteration.py `.solve` reds (L459/467/468). **Gate:** keff/MMS/twins green per
  consumer.

### Phase 4 — Retire `.solve()` + `CAP_SOLVE`/`MissingCapability` (solve axis)
- Once no consumer calls `.solve()`, **delete** the `.solve()` methods (aggressive
  retirement; rewire their behavioral tests to `.inverse().apply`). Retire `CAP_SOLVE` and
  the solve-axis `MissingCapability` guards.
- Decide `apply`/`apply_transpose` capability strings: apply is universal; put
  `apply_transpose`/`.H` on the base ⇒ **delete the `# type: ignore`s** (operator.py:754
  etc.). Keep a frozenset ONLY if the zero-scaling value-dependent edge needs it.
- Migrate the capability-closure tests to the new shape (an apply-only/singular op has no
  `.inverse()`). **Gate:** twins green; `# type: ignore`s deleted; capability reds clear.

### Phase 5 — Composition algebra in return types + cleanup + docs
- `@`/`+`/`.inverse()` return types reflect invertibility. Final pyright sweep: B4 cluster
  → 0; re-baseline ratchet. Remove dead capability scaffolding. Consider tightening the
  `["FullField"]` string proxy only if layering allows (else leave — §7 guard).
- **Sphinx (Cardinal Rule 3):** this is a real architectural milestone — update the
  operator-algebra theory/docs + the SN theory page "Development history" changelog
  (`docs/theory/discrete_ordinates.rst`): the inverse-as-operator design, why
  capabilities-as-strings was retired, the dual-pair `apply`/`inverse().apply`,
  sweep-vs-Krylov as inverse-operator kinds. Rebuild Sphinx clean.

---

## 6. Decision points / forks (the user steers these)

1. **Projection/Reconstruction rehome (Phase 1):** bring them in as first-class
   `[Rich, Coeffs]` now (architecturally right per the report; bigger blast radius) vs
   alias endomorphically and defer. *Lean: rehome — it's the report's design and removes
   the quarantine.*
2. **How far in one step:** land Phases 1–5 as one branch with per-phase gates (recommended)
   vs ship Phase 1 (typing foundation) and pause.
3. **Does any capability `frozenset` survive?** Only for value-dependent singularity
   (`ScaledOperator(0·op)`). Investigate whether that edge is real in production; if not,
   retire the frozenset entirely.
4. **`solve(A, b)` free function:** keep as sugar for `A.inverse() @ b` (grand report has
   it, line 4631) vs operator-form only. *Lean: keep — it reads well at call sites.*
5. **`InvertibleOperator` shape:** it is currently `OperatorSum["FullField"]`
   (`sn/operator.py:631`). Decide whether it stays a concrete subclass or becomes the
   `SupportsInverse` Protocol the composites satisfy structurally.

---

## 7. Guard flags — do NOT cross without explicit intent
- The `["FullField"]` **string forward-ref proxy** — tightening to real `TimedFullField`
  forces a runtime/`TYPE_CHECKING` import across the `numerics`↔`transport` layer
  boundary. Leave unless Phase 5 explicitly tackles layering.
- `singledispatchmethod` `apply` over a union (`ScatteringOperator`/`FissionOperator`,
  sn/operator.py:922) — annotate the dispatch base to the within-group V, keep arms
  concrete; don't force one TypeVar over the union.
- Adjacent B5/other reds in `iteration.py`/`solver.py`/`scattering.py` — NOT this carve;
  do not touch (they belong to the pyright burndown's later clusters).

---

## 8. EXIT — back to the pyright burndown

**Trigger:** Phase 5 complete — B4 cluster (18 reds) cleared *by architecture*, several
`# type: ignore`s deleted, Sphinx rebuilt, branch merged to `main`.

**Steps to resume the burndown:**
1. `python -m tests._harness.pyright_ratchet --update` (orpheus/ will have dropped from
   **420**; B4's 18 + the deleted ignores). Commit the baseline; confirm
   `tests/test_pyright_ratchet.py` green.
2. Open **`.claude/plans/pyright_signal_cleanup.md`** → "### Resume here (next)". Mark
   **B4 = DONE via the operator-inverse-algebra carve** (not patched — re-architected).
   Update the orpheus/ count.
3. Resume the remaining clusters, in order of value:
   - **B5 — union dispatch** (`scattering.py` ~13: `ndarray | ScalarSourceSink | None`,
     `.mesh`/`.values`/`.spatial_moments_per_axis` on raw ndarray; `solver.py` BulkField vs
     AngularField). Needs source-build flow understanding; moderate risk.
   - **units.py (5)** — pint-stub-fighting (`PlainUnit | Unknown` vs `Unit`); third-party
     category, low value — defer with B6.
   - **B6 + Workstream C** — derivations SymPy (~228) + non-production trees policy.
   - **Workstream A** — LSP-noise hook (user-gated; touches commit-protected `.claude/hooks/`).
4. Operating rules unchanged: CLI is the oracle; **NO `# type: ignore`**; re-tighten the
   ratchet after each cluster; per-cluster ff-merge + push.

**Open, unrelated:** **#250** (stale SPHERE bit-identity snapshots from the W1 clamp fix
`b2d8a6d` — a V&V re-baseline, not pyright). Independent of this carve.

---

## 9. References
- **Grand report:** `.claude/plans/neutron_transport_grand_report_v3.md` — §1 (line 76
  `K = A_loss.inverse() @ F`), §3.1 (line 157 capability Protocols; line 170
  `SupportsApply[X,Y]`), §3.3 (lines 275–332 `LinearOperator(ABC)` + taxonomy incl.
  `GreenOperator`/`SweepOperator`/`ProjectionOperator`), §22 (~4377–4415 resolvent),
  line 4159 (`loss.inverse() @ production`), line 4631 (`solve(A.H, ...)`).
- **B4 verification plan:** `.claude/plans/issue_226_b4_operator_generics_verification.md`.
- **Operator-generics design map:** `.claude/agent-memory/explorer/issue_226_operator_generics_map.md`.
- **Pyright burndown (exit target):** `.claude/plans/pyright_signal_cleanup.md`.
- **Skills/rules:** `coding-elegance` (Pattern 1 dunder algebra, Pattern 2 single source,
  Pattern 4 illegal-states-unrepresentable, Pattern 5 build-the-primitive), `vv-principles`
  (bit-identity vs principled; Mode 8 `-O` strip; Mode 11 gate-never-executes),
  `.claude/rules/delegation.md` (surgical = main-agent-direct, user-steered),
  `.claude/rules/coding-standards.md` (clean-before-extend, aggressive retirement).
