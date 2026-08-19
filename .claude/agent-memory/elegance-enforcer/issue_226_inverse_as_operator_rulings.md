---
name: issue-226-inverse-as-operator-rulings
description: Durable design-review rulings from #226 taxonomy §12 step-1 (.inverse() on every advertiser) AND step-2 (reify M + windowed product); reusable for the Green/Matrix inverse siblings and future forward composites
metadata:
  type: project
---

# #226 inverse-as-operator — step-1 review rulings (reviewed 2026-07-01)

Branch `refactor/inverse-as-operator` (uncommitted step-1 delivery over HEAD `69ed531`).
Plan: `.claude/plans/archive/operator_machinery_taxonomy.md` (§13 name-earning invariants,
§12 phase order, §9 keystone). **Merge-status goes stale — reconcile against git
(`git merge-base --is-ancestor`) before treating any of this as landed.**

**Why:** these are the reusable rulings for when the sibling inverse types
(`GreenOperator` §12 step 4, `MatrixInverseOperator` §12 step 5) land — the
review should apply the SAME criteria, and the back-half twin below collapses at
that point.

**How to apply:**

- **DELEGATION-not-reciprocal is the right inverse for value-bearing leaves.**
  `DiagonalOperator.inverse()` / `MultiplicationOperator.inverse()` return
  `InverseOperator(self)` whose `apply` DELEGATES to the leaf's `solve` (division),
  NOT `DiagonalOperator(1/c)`. Three reasons, strongest last: (1) no division-vs-
  reciprocal ULP twin while `solve` coexists; (2) no units-dishonest `1/Σ`
  "reciprocal cross-section" (a mean free path — anti-pattern #13); (3) the
  reciprocal breaks I2 involution at ULP (`1/(1/c) ≠ c`), delegation makes
  `(A⁻¹)⁻¹ is A` hold by OBJECT IDENTITY, exact. This deviates from the §9 plan
  text ("reciprocal") deliberately and correctly.

- **The faithfulness keystone FORCES `solve = inner.apply` on inverse objects.**
  `is_invertible ≡ CAP_SOLVE∈caps` (`tests/_harness/predicates.py`). Wanting I2
  (`is_invertible=True` so you can invert an inverse) ⟹ must carry `CAP_SOLVE` ⟹
  must have a real `solve`. `A⁻¹.solve(b) = A.apply(b)` is mathematically exact
  (solving the inverse IS applying the forward), not filler. Retires at P4 with
  the whole `.solve` surface.

- **RECURRING TWIN (collapse trigger set): the wrap-delegate back-half.**
  `SweepOperator` (sn/operators/sweep_operator.py) and `InverseOperator`
  (numerics/operator.py) carry a byte-identical back-half — `capabilities`,
  `domain→inner.codomain`, `codomain→inner.domain`, `solve→inner.apply`,
  `is_invertible→True`, `inverse()→inner`, repr shape — differing only in `apply`
  (Sweep threads `initial_guess`) + the `__init__` guard + carrier typing.
  Coextensive today (pinned by `test_i1_inverse_apply_is_solve_bit_identical`) ⟹
  NOT a VIOLATION, a forward-looking CONCERN. **Collapse trigger:** when
  `GreenOperator`/`MatrixInverseOperator` (the 3rd/4th wrap-delegate siblings)
  land, extract the back-half into a shared mixin/base (leaves override
  `apply`+`repr`); do NOT hand-re-derive it a 3rd time (keystone-correctness
  habitat). Taxonomy §3/§8 "siblings, no ABC" is about not sharing an INVARIANT —
  a mechanism mixin is a different axis, not precluded. Remedy NOW = reciprocal
  cross-ref comments + tracked trigger, not premature unification (Pattern 6).

- **STRUCTURAL inverses (Product/Scaled/Tensor/Permutation/Identity) route
  through their own constructors** — `(AB)⁻¹=B⁻¹A⁻¹` re-runs the capability-closure
  + composability guard on the inverse (Pattern 4 ∩ 2). Domain algebra coheres:
  the forward guard `a.domain==b.codomain` is identical to the inverse guard
  `b⁻¹.domain==a⁻¹.codomain`. "Bit-identical to solve" holds by induction; the
  gate uses a NON-commuting `D@P` so the factor-order swap has teeth.

- **The base `LinearOperator` Protocol carries NO `solve` and NO `inverse`** (only
  `apply`/`H`/`is_invertible`/`is_adjointable`). So `_SolveBackedLeaf` (module-
  private) is genuinely needed for `inner.solve` to type-check, and it correctly
  stays private (solve is a dying surface — never mint a public `SupportsSolve`).

- **The `_BoundBoundaryOperator` shim gap (pre-existing, keystone-armed):** it
  forwards `capabilities` (⟹ advertises `CAP_SOLVE` when inner is a reflective
  Permutation-based law) but has NO `solve` method ⟹ the standard
  `if CAP_SOLVE in caps: op.solve()` guard passes then AttributeErrors. Latent
  (no consumer, §17), self-resolves at P4. When reviewing shim promise-completion,
  check ALL forwarded caps have a backing method, not just the one being added.

---

# #226 step-2 review rulings (reviewed 2026-07-01, over base `e7115a2`)

Step 2 = "Reify M + windowed product" (§17 W1/W2). Uncommitted working tree on
`refactor/inverse-as-operator` (step 1 now COMMITTED as `e7115a2`). Verdict:
**BLOCK — but the ONLY hard gate is retirement-completeness (doc blast radius);
the code + tests are ship-quality.** Reusable rulings:

- **W2 reify-M is honest and correct.** `ScheduledInvertibleOperator(OperatorSum)`
  = M = `(LC, ScaledOperator(-1, B_lower))` — `apply` is the inherited leaf sum
  (no override), `inverse()→SweepOperator(self)`, `solve→invertible._solve_timed_
  full_field(schedule=…, reflect=lower.reflect_rows_inplace)`. This dissolves the
  old `_GaussSeidelResolvent` whose `apply=(L+C)` / `solve=(L+C−B_lower)⁻¹` were
  inverses of DIFFERENT operators (round-trip defect 2.667). PASS.

- **The SweepInvertible union widening is the RIGHT move (avoids a 3rd wrap-
  delegate class).** `SweepOperator.inner: Union[InvertibleOperator,
  ScheduledInvertibleOperator]` — the existing inverse wrapper ABSORBS the new
  forward type instead of minting a parallel `ScheduledSweepOperator`. Exactly the
  step-1 "don't hand-re-derive the back-half" ruling, applied preventively.

- **NEW twin (shape #1, forward side) — CONCERN not VIOLATION.** Step 1 flagged the
  INVERSE-side wrap-delegate back-half (SweepOperator↔InverseOperator, shape #2).
  Step 2 creates a 2nd member of the FORWARD-side back-half (shape #1):
  `InvertibleOperator` + `ScheduledInvertibleOperator` now BOTH carry
  `is_invertible→True` + `inverse()→SweepOperator(self)` + `solve()→
  self._solve_timed_full_field(rhs, initial_guess=…)` (3 trivial methods,
  byte-coextensive; only `_solve_timed_full_field` bodies differ). Per L2 (collapse
  on the 3rd, not the 2nd) + Pattern 6 (a mixin ADDS a concept to save 3 trivial
  lines ⟹ concept-count goes UP) ⟹ DEFER. Remedy NOW = reciprocal cross-ref
  (ScheduledInvertibleOperator→InvertibleOperator already exists one-way; add the
  back-ref) + tracked collapse trigger. NB shape #1's trigger differs from shape
  #2's: Green/Matrix are INVERSE-side (shape #2), so shape #1 only grows if a new
  FORWARD composite appears.

- **"Two reflect verbs" = genuinely different semantics, single-sourced core =
  PASS (not a twin).** `SNBoundaryOperator.reflect_inflow_inplace` (whole-face
  ASSIGNMENT `ψ.inflow ← B·ψ`, for the SI direct loop + reconstruction) vs
  `SNMaskedBoundaryOperator.reflect_rows_inplace` (row-masked ADDITIVE
  `bf[rows] += (B·bf)[rows]`, the inhomogeneous forward-subst row `z_in=y_row+
  (Bz)_row`). Both route through the ONE `_reflect_trace`. This is TWO WRITE-BACK
  POLICIES over one core, not two implementations. The additive-vs-overwrite is THE
  load-bearing correctness distinction (the old whole-face overwrite dropped y_row →
  O(1)-wrong as an inverse). Reusable shape: when you see two `*_inplace` verbs,
  check (a) shared core and (b) whether the write semantics genuinely differ before
  crying twin.

- **Windowing-as-composition (W1) = the codomain-changing-emit smell resolved
  right.** `WindowedSweep(OperatorProduct) = P @ A.inverse()`, P =
  `BulkAnalysisOperator` (frame.analysis on bulk ⊕ Id on trace) — a block
  coisometry with caps `{apply}`, so `OperatorProduct`'s capability closure gives
  `is_invertible=False` + no `CAP_SOLVE` FOR FREE (honest "no round-trip promise").
  `WindowedSweep.apply` OVERRIDES with the fused substrate moment-emit; the
  inherited `OperatorProduct.apply` body (`a.apply(b.apply(x))` = `P.apply(A⁻¹.
  apply(rhs))`) IS the deforestation oracle, kept per the fuller-view exception.
  `codomain=None` is the HONEST "moment-bulk FunctionSpace not minted yet" (faking
  the angular space would be a LIE) — Pattern-6 deferral, tracked. Retires the
  public `solve_moments` (codomain-change-as-a-config).

- **The sweep-realized-inverse SOURCE-SUBSPACE restriction is family-wide, not
  step-2-specific.** The sweep `shed` OVERWRITES outflow-definition rows, so
  `M.solve` (and `(L+C).solve`) realize the true inverse only on `{y: outflow-rows
  = 0}` (all production rhs). Deferring a raise-guard is consistent (a guard would
  break the step-1 keystone gate that feeds `A.apply(random x)` AND would have to be
  applied to the whole family). WATCH: when the taxonomy's `M⁻¹A` preconditioner
  algebra lands (grander-scheme verdict), a residual `r=b−Ax` has nonzero outflow
  rows ⟹ silently projected. The family needs a louder contract THEN.

## Step-2 findings delivered (for continuity if a step-2 re-review happens)

1. **MUST-FIX (the block): retirement blast radius — stale docs.** The carve deleted
   `_GaussSeidelResolvent` + `solve_moments` (valid at `e7115a2`) but left
   PRESENT-TENSE dead-contract refs across 3 theory pages: `discrete_ordinates.rst`
   :11594 ("the resolvent … IS `_GaussSeidelResolvent`") + :11786 (the OLD factory
   return `(_GaussSeidelResolvent(L+C,B,schedule),(S,))` — now `((L+C)-B_lower,
   (S,B_upper))`) + :14407/:14437/:17833; `operator_algebra.rst` :7418/:7482/:7704-5/
   :7850/:7869-72 (a whole section on the retired solve_moments/resolvent story);
   `cross_section_data.rst` :277. Sphinx has NO `-W`/nitpicky ⟹ silent rot (Rule 3
   theory page). operator_algebra.rst needs a real rewrite, not find-replace.
2. **SHOULD-FIX: `reflect_rows_inplace` docstring over-claims** "exact on an
   ARBITRARY rhs (the W2-round-trip keystone)" — the keystone round-trips a
   CONSISTENT state, and exactness is only on `{outflow-rows=0}` (per the test
   module docstring). "Arbitrary" is true only for the LOWER-INFLOW rows; reword so
   the two docstrings are coextensive.
3. **SHOULD-FIX: intra-carve typed-honesty inconsistency.** `windowing.py:201/206`
   narrows `self.a`/`self.b` with `# type: ignore[return-value]`; the SIBLING
   `scheduled_invertible.py` (same carve) uses `cast(...)`. Unify on `cast`
   (anti-pattern #19 — declare the narrowing the isinstance guard proved).
4. **NOTE: transitional `_MomentWindowedResolvent` two-hat** (`apply`→base forward
   ANGULAR, `solve`→product MOMENTS) transiently re-introduces the codomain-mismatch
   the carve dissolves — OK because step 3 deletes it; verify step 3 lands.
5. **NOTE (test-arch): windowed×G-S corner** (`P @ M⁻¹`, the DEFAULT for 2-D
   Cartesian eigenvalue) pinned only at integration/ℓ=0 (FP-equivalence), not at the
   unit/ℓ≥1 deforestation level (that gate uses `P @ (L+C)⁻¹`). Acceptable: the
   moment kernel and the boundary schedule are ORTHOGONAL and each is independently
   pinned — but a corner worth a dedicated pin.

---

# #226 step-3 review rulings (reviewed 2026-07-02, over base `cc293ef`)

Step 3 = "Solver builds the inverse operator; SI/Krylov apply it — dissolves the
duck-typed resolvents + consumer migration in one move." Uncommitted working tree
on `refactor/inverse-as-operator` (step 2 now COMMITTED `cc293ef`). Verdict:
**SHIP-WITH-FIXES — code+tests ship-quality; the sole MUST-FIX-before-MERGE is the
doc-retirement blast radius (deferrable to archivist, but scope must EXPAND).**
Mirrors the step-1/step-2 shape (code ships, doc-retirement is the gate).

- **The step-3 mandate is DELIVERED elegantly.** `SourceIteration` now
  `self.L_inv.apply(rhs, initial_guess=psi_prev)` UNCONDITIONALLY — the
  `inspect.signature` probe + `_solve_with_seed`/`_solve_accepts_seed` DELETED (the
  duck-typed resolvent dissolved into the typed inverse). `_MomentWindowedResolvent`
  DELETED (my step-2-demanded two-hat retirement). `_maybe_window` literally spells
  `BulkAnalysisOperator(...) @ sweep` = `P @ A⁻¹`. Master-standard PASS.

- **RULING — the seeded-apply contract is CONVENTION, not STRUCTURE, across the
  operator family.** Only `InverseOperator`/`SweepOperator`/`WindowedSweep` (explicit
  `initial_guess` kwarg) + `ScaledOperator` (transparent `**kwextra` forward) accept
  the seed; `OperatorProduct`/`IdentityOperator`/`PermutationOperator`/
  `TensorProductOperator`/base-leaf apply (`apply(self,x,/)`) do NOT. So KEigenvalue's
  `cast("SupportsSeededApply[Any]", cast("SupportsInverse", self.L).inverse())` is
  runtime-safe ONLY for the REACHABLE set (SN loss → `SweepOperator`; L0 test →
  `InverseOperator`). The accompanying comment "**every inverse the family returns
  carries the canonical seeded-apply signature**" is a FALSE UNIVERSAL (anti-#11
  false-invariant comment) — CONCERN not VIOLATION (latent: no current consumer
  composes a structural-inverse `L` into KEigenvalue; but `is_invertible` passes for
  `A@B`, then TypeError at the first inner `apply`). SHOULD-FIX: scope the comment to
  the reachable inverses. **Reusable when Green/Matrix land:** they are the SEED-LESS
  inverse kind → they MUST `del initial_guess` like `InverseOperator` to keep the
  family contract; verify their `.apply` accepts the kwarg.

- **`del <arg>  # <reason>` IS the codebase idiom** for a Protocol-mandated
  accepted-and-ignored arg (precedent `sn/spatial/psi_half_angle_seed.py:405
  del context  # unused — Protocol-shape compatibility only`; also carlvik/spectrum).
  `del initial_guess` on `InverseOperator.apply`/`WindowedSweep.apply` is CORRECT and
  consistent — PASS (scrutiny a). Not `_ = x` (that idiom is reserved for touching a
  property for its side-effect, e.g. `_ = self.mat_xs.absorption_cross_section`).

- **The double-cast is the posing-layer shape; rule-of-two helper candidate.**
  `_seeded_inverse(L) -> SupportsSeededApply[Any]` would serve BOTH cast sites (both
  in iteration.py: KEigenvalue `L_inv` + Krylov preconditioner `.apply`) and read
  like the domain ("build the seeded inverse"). SHOULD-FIX (defensible to keep — casts
  are honest per L-009). Fold the scoped seeded-apply rationale (ruling above) into
  the helper docstring. NB solver.py needs NO cast (base_resolvent is concretely
  `InvertibleOperator`/`ScheduledInvertibleOperator`) — the 2 cast sites are exactly
  KEigenvalue + Krylov.

- **`SupportsSeededApply(Protocol[V])` name/placement/bound all PASS** (scrutiny c):
  mirrors `SupportsInverse` (static-only, not runtime_checkable), lives in iteration.py
  (the consumer's contract), `V bound=Vector`. `rhs` is positional-only (`/`) in the
  Protocol → each implementor names its first param in its own vocab (`rhs` sweeps /
  `x` generic inverse) without breaking conformance (widening pos-only→pos-or-kw is
  safe). Deliberate, not drift.

- **L→L_inv migration COMPLETE** (scrutiny d PASS): all 3 prod (`iteration.py:1050`
  KEigenvalue._inner, `solver.py:539`) + all test `SourceIteration(...)` sites pass a
  pre-built inverse; no forward-as-inverse latent bug; no external reader of the
  renamed attr. Krylov correctly keeps `self.L` (FORWARD — its matvec applies it).

## Step-3 findings delivered

1. **MUST-FIX-before-MERGE: doc-retirement blast radius (deferrable to archivist,
   scope EXPAND).** 6 broken `:class:`/`:meth:` cross-refs to the deleted
   `_MomentWindowedResolvent` in `operator_algebra.rst` (7258, 7422, 7483, 7493, 7932,
   8065) — Sphinx has NO `-n`/`-W` ⟹ silent plain-text rot (Rule 3); 7258 is
   present-tense-FALSE ("behind the … adapter, is therefore never even constructed").
   **The brief UNDERCOUNTS:** the step-2 debt (`_GaussSeidelResolvent`/`solve_moments`)
   is STILL OPEN — `operator_algebra.rst` (10) + `discrete_ordinates.rst` (6). Those
   are mostly literals (render fine) ⟹ lower severity, but present-tense prose claims
   among them are MUST-FIX. operator_algebra.rst "needs a real rewrite, not
   find-replace" (step-2 note stands). Discriminator applied: `:class:`/`:meth:`
   cross-ref to a deleted symbol = broken+silent = higher-severity than a literal.
2. **SHOULD-FIX: KEigenvalue universal seeded-apply claim (see ruling above)** — scope
   the comment; file/track the structural-vs-convention question for the taxonomy.
3. **SHOULD-FIX: extract `_seeded_inverse` helper (see ruling above).**
4. **PASS callouts:** two-hat retirement done; test-migration textbook
   (`requires_solve_on_L` → `invertibility_obligation_lives_at_the_inverse_builder` +
   `keigenvalue_requires_invertible_L` + `ApplyOnlyStep` end-to-end; GAP-3 migrated);
   my step-2-demanded `P@M⁻¹` deforestation pin landed with all 3 legs; seed spy is
   route-invariant/`-O`-proof/values-equal; pyright ratchet 152→148 (clean reduction).
   Verified green: 82 gates in 0.38 s; ratchet delta confirmed.

---

# #226 step-4 review rulings (reviewed 2026-07-02, over base `1ab7429`)

Step 4 = "GreenOperator + OperatorSum.inverse() + InverseWrapMixin extraction" (§12
step 4). Uncommitted working tree (step 3 now COMMITTED `1ab7429`). Verdict:
**SHIP-WITH-FIXES — ZERO code defects; all 3 fixes are doc-sync.** The step-1
predicted collapse EXECUTED correctly. Reusable rulings for step 5
(`MatrixInverseOperator`):

- **The mixin extraction is the step-1 prediction, executed right.** `InverseWrapMixin
  (Generic[_ForwardT], metaclass=ABCMeta)` abstracts the byte-identical back-half
  (`capabilities` ClassVar `{apply,solve}`, domain↔codomain swap, `solve→inner.apply`,
  `is_invertible→True`, `inverse()→inner` by object identity, `__init__`). The 3 siblings
  (`InverseOperator`/`SweepOperator`/`GreenOperator`) keep EXACTLY guard+apply+repr — the
  cut my step-1 ruling specified. Fired at the 3rd sibling per defer-≥2/extract-at-3
  (Pattern 6). No metaclass conflict (`_ProtocolMeta` ⊃ `ABCMeta`). `capabilities` as a
  ClassVar is a single-source IMPROVEMENT (value is instance-invariant) — and NOT the
  dataclass-field trap (mixin is a plain ABCMeta class, not a dataclass). PASS.

- **The #285 STRUCTURAL resolution is real and load-bearing.** The mixin's abstract
  `apply(x, /, *, initial_guess=None)` means ABCMeta blocks a missing override AND pyright
  rejects a kwarg-dropping one — a new sibling CANNOT forget the seed keyword. This
  UPGRADES the step-3 "seeded-apply is CONVENTION not STRUCTURE" ruling: it is now
  STRUCTURE for the wrap-delegate family; the residue is `OperatorProduct.inverse` (a
  composition, apply is `(self, x, /)` — no kwarg), closes at step 5.

- **`_SolveBackedLeaf → _InvertibleForward` rename = accuracy fix, not churn.** Green's
  `inner` is an `OperatorSum` (NOT a leaf), so "Leaf" was wrong; "InvertibleForward"
  generalizes. Bound requires `solve` (needed by `InverseOperator.apply→inner.solve`; a
  union contract Green satisfies but doesn't use — transitional, retires at P4 with the
  whole `.solve` surface). 0 residue of the old name.

- **Refinement loop = honest wrap (attack #1).** The loop relocates ONLY the stopping test
  (driver's ρ-blind increment → true residual `‖(A−B)ψ−q‖/‖q‖`); the splitting math stays
  in `SourceIteration`. Proven 3 ways: bit-identical-to-hand-SI (`assert_array_equal`,
  0-ULP — a re-impl can't be bit-identical); the §18.A near-critical gate (budget=600
  BETWEEN increment-stop ~460 & refinement-close ~920 ⟹ increment-mutant silently wrong,
  honest Green raises; budget=5000 = positive control the promise is DRIVEN not check-only,
  which would falsely raise ∀ρ>1/2). `steps += max(len(history),1)` counts total Richardson
  steps; the `max(…,1)` guards the empty-history (`max_iter=0`) infinite loop. Residual
  matvec correctly EXCLUDED from the "splitting-step budget" (a check is not a step). Re-seed
  preserves the Neumann-partial-sum identity (Richardson is memoryless).

- **`OperatorSum.is_invertible`=leading-term is honest (attack #2).** `is_invertible=True`
  means "inverse OPERATOR constructible," NOT "inversion guaranteed" — the faithfulness-
  keystone semantics carried to sums; divergence raises `ConvergenceFailure` LOUDLY
  (Cardinal Rule 1; the C+L trap pins it). NO prod consumer of `OperatorSum.solve` (grep
  verified: only `SourceIteration.solve`/`krylov.solve` in `orpheus/`) ⟹ `solve=fresh-Green
  -per-call` is off-hot-path, transitional, retires P4. `_left_spine_terms` uses
  `while type(node) is OperatorSum` (EXACT type, not isinstance) — load-bearing: stops at
  fused `InvertibleOperator` leaves, preserving the MRO shadow (`(L+C).inverse()→Sweep`,
  `a_loss.inverse()→Green`).

- **The forward-side shape-#1 twin stays correctly DEFERRED.** `OperatorSum` now also carries
  `is_invertible/inverse/solve`, but its bodies DIVERGE from InvertibleOperator/
  ScheduledInvertibleOperator (`inverse()→Green` not `Sweep`; `solve→inverse().apply`) — so
  it's honest polymorphism, NOT a 3rd byte-identical twin. Shape #1 still = 2 members;
  trigger (a new FORWARD composite) not met. My step-2 ruling holds.

- **NEW smell flavor (add to the tells): FORWARD-staleness — a docstring calls a
  JUST-DELIVERED thing "deferred / no consumer yet / open."** The retirement-side L-001 is
  "doc names a primitive the code doesn't call"; this is its mirror — the carve DELIVERS X
  but a doc still frames X as pending (sweep_operator.py:12 "Green deferred, no consumer
  yet"; operator_algebra.rst:8285 "#285 … the open decision … lands when GreenOperator do").
  Grep on any additive carve: `grep -niE "<newthing>.{0,40}(defer|no consumer|not yet|
  pattern 6|future|open decision)"`. SHOULD-FIX (in-diff code docstring = fix inline).

## Step-4 findings delivered

1. **SHOULD-FIX inline: sweep_operator.py:12** "a Krylov GreenOperator (deferred — Pattern 6,
   no consumer yet)" — delivered now; also "Krylov" contradicts Green's own naming (Krylov =
   realization axis, `KrylovInverseOperator` was rejected). → "iterative/Green, delivered".
2. **SHOULD-FIX (archivist): operator_algebra.rst:8285-8297** frames #285 open/future when the
   carve resolved it structural (iteration.py docstrings already updated — theory drifted).
3. **FOLLOW-UP (Rule 3, archivist): Green theory-page section** — only 2 incidental mentions
   (788, 8293); the inline module docstring is Rule-3-quality and suffices short-term.
4. **PASS callouts:** prior doc-retirement gate CLOSED (no broken `:class:` xrefs; historical
   framing; `discrete_ordinates.rst:11604` names the live public entry); my step-3
   SHOULD-FIXes landed (`_seeded_inverse` extracted single-source; KEigenvalue comment
   scoped). Verified: 77 green/0.31s (`-O` serial); pyright 0-new (green_operator/operator/
   sweep_operator/streaming all 0; the 7 iteration.py errors are pre-existing untouched lines).

---

# #226 step-5 review rulings (reviewed 2026-07-02, over base `9333305`)

Step 5 = "`MatrixInverseOperator` + `LinearOperator.as_matrix()`" (§12 step 5 — the
TERMINAL taxonomy step). Uncommitted working tree + 1 new file
(`orpheus/numerics/matrix_inverse_operator.py`). Verdict: **SHIP-WITH-FIXES — ZERO code
defects; the ONE SHOULD-FIX is a docstring-vs-truth drift the diff itself introduced.** The
step-1/step-4 predictions all EXECUTED right. Reusable rulings for the earn-the-class
consumers (CP `[P]` §14b; homogeneous full-operator spelling task #138):

- **The 4th sibling joined `InverseWrapMixin` EXACTLY as predicted** (step-1/step-4): keeps
  only guard+apply+repr; `del initial_guess  # exact direct inverse — M-direct IS
  seed-independence` (the seed-less contract step-3 mandated for Green/Matrix). The
  `_InvertibleForward`→`_WrappedForward` bound RELAXATION is an ACCURACY win *earned by the
  4th sibling*: MatrixInverseOperator provably never reads `inner.solve`/`inner.is_invertible`
  (it inverts the MATERIALIZATION — values not structure), exposing the mixin's TRUE minimal
  bound (domain/codomain/apply). `_InvertibleForward` survives as InverseOperator's narrowing
  (+is_invertible/solve). Not speculative — a real consumer forced the split. PASS.

- **`OperatorProduct.inverse()→InverseOperator(self)` closes the #285 kwarg residue** step-4
  named. Bit-identical VERIFIED against live tree: `OperatorProduct.solve` (operator.py:1174)
  = `b.solve(a.solve(q))`, so `InverseOperator(product).apply(q)` = the old reversed-product
  spelling exactly, now + the seeded-apply contract + object-identity involution. The guard is
  satisfied (inverse() checks `is_invertible` before wrapping). PASS.

- **`_resolve_basis_shape` is GENUINELY single-sourced** (the brief's drift worry = ZERO): the
  sibling imports the FUNCTION OBJECT by reference (not a copied rule); both base `as_matrix`
  and the ctor call the ONE impl; the ctor's re-pass through `as_matrix` is IDEMPOTENT on
  explicit input. Textbook Pattern 7 + Pattern 2.

- **NO memory held twice** (brief's hypothetical does NOT occur): ctor's `matrix` is a LOCAL,
  GC'd after `lu_factor`; only `self._lu` + `self._basis_shape` retained (the irreducible
  state). Optional NIT: `lu_factor(matrix, overwrite_a=True)` is SAFE (matrix is a fresh local)
  and shaves the transient 2×n² construction peak → 1×n².

- **Illegal-states unrepresentable AT CONSTRUCTION** (Cardinal Rule 1, never return a
  non-inverse): non-square→ValueError, too-large→MatrixTooLarge, exact-singular (`==0.0` U
  pivot)→LinAlgError. NEAR-singular deliberately NOT refused (priced into the M-direct κ(A)
  bound) — documented + pinned. scipy's LinAlgWarning silenced locally (catch_warnings
  restores) + replaced by the loud construction error. Every invariant the docstrings claim is
  test-pinned (seed bit-identity, machine·cond-grain w/ Green contrast, involution-object-id,
  rectangular-honest, C-order, `as_matrix ≡ retired _as_dense` oracle) — my L1 catch does NOT
  fire.

- **`as_matrix` override is LSP-clean, not ceremony:** it honors the base ERROR contract
  (inconsistent basis_shape→ValueError, tighter max_dimension→MatrixTooLarge), so it's
  substitutable on the failure paths too. `MatrixTooLarge(RuntimeError)` as a resource-effect
  (not TypeError/ValueError) with deliberately-NO `is_materializable` predicate (a pure
  resource precheck, unlike structural is_invertible/is_adjointable) is a principled
  illegal-states/resource distinction (§17 A2).

- **RETIREMENT CLEAN — a first for this campaign** (steps 2/3 both BLOCKED on doc-retirement
  blast radius; step 5 does NOT). `_as_dense` deleted: NO dangling import (test-migration
  done — the behavioral char-test rehomed to `tests/numerics/test_matrix_inverse_operator.py`
  INLINES the reference loop as a drift-guard oracle, not an import of the dead symbol); NO
  source `.rst` cross-ref (the `docs/_build/*` `_as_dense` hits are regenerated build output);
  all surviving mentions are HISTORICAL-framed ("né _as_dense", "retired _as_dense"). When
  reviewing a leaf→primitive promotion, this is the model: promote+retire+rehome-the-char-test-
  as-oracle in one move.

- **`dense_per_material` honesty fix REMOVES a false consumer claim** (elegance WIN): the old
  docstring said it was "the as_dense consumption mode for the LHS-fold solvers
  (diffusion/homogeneous)" — verified FALSE (grep: only test transpose-oracles consume it; no
  prod/diffusion consumer). New docstring correctly frames it as the STORAGE-SIDE oracle
  (L11 structurally-independent-of-`apply`'s-einsum). Direction of the fix (remove a false
  claim) is exactly right.

## Step-5 findings delivered

1. **SHOULD-FIX (docstring-vs-truth, INTRODUCED by this diff): `InverseOperator` class docstring
   is now an undercount.** operator.py:~1599-1611 says it "serves the value-bearing LEAVES
   (Diagonal/Multiplication) whose inverse action is an exact pointwise division" + the whole
   "reciprocal twin / 1/Σ mean-free-path" section is leaf-specific — but this diff made
   `OperatorProduct.inverse()` return `InverseOperator(product)`, a NON-leaf whose inverse
   action is `b.solve(a.solve(q))`, NOT a division. Code is already general+correct
   (`apply→inner.solve`); only the narrative narrows. Bug habitat: a future narrowing
   assertion (`isinstance(inner,(Diagonal,Multiplication))`) or a reciprocal-materialization
   "optimization" trusting "always division" breaks the composite case. Fix: broaden to "generic
   wrapper for any solve-backed forward lacking a more specific named inverse — value-bearing
   leaves AND invertible composites (OperatorProduct)"; scope the reciprocal-twin note to the
   leaf case.
2. **NIT (optional): two-kinds taxonomy placement.** wrap-delegate-vs-algebra-closed contrast
   lives only in `OperatorProduct.inverse()`. Stated ONCE = correct (don't duplicate onto
   Scaled/Perm/Tensor/Identity); but the family root (`InverseWrapMixin` docstring) is a more
   discoverable home. No bug habitat.
3. **NIT (optional, perf): `lu_factor(..., overwrite_a=True)`** — see memory-held-once ruling.
4. **FOLLOW-UP (Rule 3, archivist): operator_algebra.rst step-5 section** — the `as_matrix`
   functor (Op→Mat), MatrixInverseOperator (M-row), MatrixTooLarge (§17 A2), #285 closure. Inline
   module docstrings are Rule-3-quality, suffice short-term (same as Green @ step 4). Sphinx
   rebuild (standard merge step) regenerates the stale `docs/_build/*` `_as_dense` HTML.
5. **PASS callouts** above. Verified state (per brief, not re-run): 117 gates green; 14/14
   mutations red; pyright ratchet 148 (0-new); tier numerics+sn+transport 2998/0; homogeneous
   k∞ byte-identical. MatrixInverseOperator has NO prod consumer yet (tests only) — honest
   latent-consumer framing (direct_eigenvalue inline dense solve / CP `[P]` / task #138), same
   ruling GreenOperator shipped under @ step 4.

---

# #226 step-6 review rulings (reviewed 2026-07-02, over `f4919b1` — carve P4, THE TERMINAL retirement)

Step 6 = "retire the `capabilities` frozenset, both axes" (W1+W2 landed; the elegance
review IS the W3 gate). Verdict: **PASS — ZERO MUST-FIX code defects, ZERO VIOLATIONS.**
The cleanest carve of the campaign (steps 2/3 BLOCKED on doc-retirement; step 6's only
gate is the SAME recurring doc rot, now a listed W3 archivist deliverable). Reusable
rulings + the durable findings:

- **Retirement completeness is a 2-grep PASS.** `grep -rnE "CAP_|MissingCapability|
  \.capabilities" orpheus/ tests/` returns EXACTLY 6 deliberate docstring-history mentions
  (operator.py:289/:466/:483 "successor to CAP_X" — historical-framed, allowed; predicates.py
  module docstring; 2 test docstrings) and ZERO functional reads. The frozenset twin did NOT
  survive anywhere. When reviewing a symbol-retirement, this grep + tense-discrimination
  (historical "successor to X" OK / present-tense "advertises X" = rot) is the whole gate.

- **The getattr-consistency judgment (brief-requested, the sharp finding).** After a
  predicate is PROMOTED from a capability-set read to a base-Protocol property (default
  `False`), a `getattr(operand, "<predicate>", False)` read is DEAD-relative-to-type. The
  discriminator honest-vs-smell is NOT "could the operand be duck-typed?" — it is **does a
  SIBLING read the SAME operand directly?** `OperatorSum.is_invertible`:1151 uses `getattr`
  but `OperatorSum.is_adjointable`:1124 reads `self.a.is_adjointable` DIRECT, and EVERY other
  composite (`OperatorProduct` 1286/1328, `Scaled` 1394/1415, `TensorProduct` 2120/2146,
  `SNBoundaryOperator`, `_BoundBoundaryOperator` 116/120) reads DIRECT → the getattr is the
  lone odd-one-out, its `False` fallback provably dead (suite green ⟹ no operand lacks the
  attr), so it is a POLISH/NIT (coextensive today), NOT a VIOLATION. Latent bug habitat: base
  makes the predicate abstract → direct reads RED under pyright, this getattr silently keeps
  `False`, masking a missing impl asymmetrically (loud on adjoint axis, silent on invert).
  The ctor-guard getattr (KEigenvalue:1090, GreenOperator:256) are SOFTER — a posing-layer
  boundary-parse (malformed arg → not-known-invertible → refuse) is defensible; the smell is
  only that they are a THIRD spelling of the question the carve otherwise canonicalized with
  the `invertible()` bridge (used at `_seeded_inverse`:263, SAME file). Reusable: "one
  question, three spellings (bridge / direct / getattr-default)" is a Pattern-7 consistency
  NIT, sharpest where a sibling proves the getattr dead.

- **The RankOne-bug-class recurrence check is MANDATORY on any advertisement-retype.** The
  carve's own caught bug: RankOne lost its ONLY adjointability advertisement (per-instance
  caps deleted, no `is_adjointable` override) → F† dead behind a factor guard. So for EVERY
  class that removed `capabilities = {…CAP_APPLY_TRANSPOSE…}`, VERIFY a matching
  `is_adjointable` override exists (directly OR inherited from a role base). Confirmed clean
  here: scattering/fission/isotropic each have `def is_adjointable`; frame faces inherit
  `AnalysisOperator.is_adjointable=True` (role-level) / `_FrameReconstruction` overrides
  (reconstruction role is apply-only base-False). The teeth: `test_capability_survival.py`
  enumerates every advertiser + pins its `(is_invertible, is_adjointable, inverse_contract)`
  triple via the keystone, mutation-verified (M-ADJ-EAGER). Grep `def is_adjointable` in each
  retyped file; do NOT trust the docstring's `is_adjointable=True` claim (L-001).

- **Keystone v2 is the model for a coexistence-scaffold retirement.** `assert_capability_
  faithful` (pinned `is_X ≡ CAP_X ∈ caps` — could only run WHILE the frozenset existed)
  DELETED with a rationale docstring ("its job, licensing the deletion, is done"); replaced
  by the permanent `assert_inverse_adjoint_contract` referencing NO frozenset, single-sourced
  across both enumerations (numerics + SN/transport). Three legs: three-valued inverse
  contract (INVERTIBLE/VALUE_RAISE/STRUCTURAL_ABSENT via `hasattr` — the runtime shadow of
  declared-vs-not) + eager-`.H` return-or-raise + **bridge-consistency leg** (`_inv_bridge(op)
  == op.is_invertible`) which is what gives the one-line TypeGuards their teeth. Explicit
  `raise` not `assert` (un-collected helper, `-O`-strips asserts, Mode-8). This is how a
  licensing-scaffold retires: delete it WITH the thing it licensed, keep its permanent
  successor.

- **`_seeded_inverse` two-kinds dispatch = type-as-structure done right.** Keyed by
  `isinstance(inv, InverseWrapMixin)` (`_wrap_delegate_member` TypeGuard) — STRUCTURAL family
  membership, NOT a signature/`inspect` probe: mixin members carry the seeded `apply` by
  construction (#285) → returned as-is; algebra-closed inverses (Perm/Identity/Scaled/TP) are
  plain forwards → wrapped in the accept-and-drop `_SeededExactApply`. The TypeGuard is
  genuinely load-bearing (without it `return inv` fails to narrow to `SupportsSeededApply`).
  Fixes the pre-existing crash on algebra-closed preconditioner heads at the ONE home.

## Step-6 findings delivered

1. **POLISH-1 (do-now, 1 line): `OperatorSum.is_invertible` operator.py:1151** — the lone
   getattr in a direct-read family (see ruling). Fix `return self.a.is_invertible` to match
   `is_adjointable`:1124.
2. **POLISH-2 (optional, consistency): KEigenvalue iteration.py:1090 + GreenOperator
   green_operator.py:256** — the two ctor-guard getattr, a 3rd spelling of the invertibility
   question; prefer direct `x.is_invertible` (base guarantees it) or the `invertible()` bridge.
3. **FILE-AS-ISSUE (W3 archivist Sphinx pass — the SAME doc gate as steps 2/3): 13 BROKEN
   `:class:MissingCapability` xrefs** (operator_algebra.rst 9 / galerkin_projection.rst 2 /
   discrete_ordinates.rst 2 — deleted symbol → silent plain-text rot, no `-W`) + stale
   CAP_/`capabilities`/`OperatorSum.solve` prose across 6 theory pages. operator_algebra.rst's
   whole "honest CAP_SOLVE gate" section (~:1325-1348, the MultiplicationOperator story) needs
   a REWRITE to the NotInvertible/`min|f|>0` surface, NOT find-replace.
4. **NOTE (pre-existing, not carve-introduced, optional co-target): `# type: ignore[override]`
   at _bound_compat.py:141** (block_role property-over-attribute) — the carve just DELETED a
   sibling `# type: ignore[override]` on the retired `capabilities` property 2 lines up; the
   user rejects `# type: ignore` (#255 tracks). Fix opportunistically or leave to #255.
5. **PASS callouts:** repr getattr:824-825 is the ONE defensible getattr (a repr must never
   raise) — keep it; §44.D honesty fix landed (`_InvertibleForward.solve` docstring = "permanent
   face"); §44.F resolved (repr axis-flags); my step-5 SHOULD-FIX landed (InverseOperator
   docstring broadened LEAVES+COMPOSITES); the step-1 keystone-armed `_BoundBoundaryOperator`
   caps-without-solve gap is CLOSED (no solve stub, forwards inverse()/apply_transpose).

---

# #226 step-5b review rulings (reviewed 2026-07-02, over `cb62310` — task #138, the FIRST prod consumer)

Step 5b = "wire the homogeneous solver onto MatrixInverseOperator: `K =
MatrixInverseOperator(loss) @ production`, `dominant_eigenpair(K.as_matrix())`; extract the
Perron–Frobenius validation into public `dominant_eigenpair`; fold the sign convention into
private `_sign_normalised`." Uncommitted tree. Verdict: **SHIP-WITH-FIXES — ZERO code defects;
the 2 SHOULD-FIX are BOTH doc-staleness (one in-file, one blast-radius-outside-diff), 1 POLISH.**
The cleanest consumer-wiring of the campaign. Reusable rulings:

- **The dominant_eigenpair extraction is a twin-in-waiting DISSOLVED (Pattern 2 win).** The
  argmax-real / complex-reject / real-cast / sign-normalise body was ABOUT to be copied for the
  K-path; instead it single-sources into public `dominant_eigenpair(M)` and `direct_eigenvalue`
  delegates via `dominant_eigenpair(np.linalg.solve(A,F))`. Byte-verified: the extraction is
  identical to the old inline (only additions = boundary `asarray(dtype=float)` + squareness
  guard, both no-ops on the solve(A,F) path → direct_eigenvalue stays bit-identical). Two ENTRIES
  (posed (A,F) engine vs operator-algebra K-path), ONE kernel = honest polymorphism, NOT a twin
  (mirrors step-3 "two entries one kernel = PASS"). One-home pinned by two monkeypatch
  neuter-proofs (complex axis + sign axis).

- **`_sign_normalised` folds a REAL byte-identical twin.** `direct_eigenvalue`'s
  `if phi.sum()<0.0: phi=-phi` and RQI's `v if v.sum()>=0.0 else -v` are the SAME map (flip iff
  sum<0, incl. the sum==0 corner → both keep v). Fold is exact. Pattern 7 (normalise at the
  definition site). Correctly PRIVATE (module-underscore) — a convention detail with no
  cross-module consumer; exposing it would invite re-application (anti-Pattern-7). `dominant_eigenpair`
  correctly PUBLIC (cross-module consumer = homogeneous solver + own test class).

- **The two `basis_shape=(ng,1)` are HONEST per-call, NOT duplication (scrutiny c PASS).** Verified
  against live tree: `loss` (MultiplicationOperator/from_mesh) and `production` (FissionOperator)
  are meshless → `space=None`/`full_field_space=None` → `domain=None`; `OperatorProduct.domain`
  = `self.b.domain` = production.domain = None. So BOTH `MatrixInverseOperator(loss, basis_shape=)`
  AND `K.as_matrix(basis_shape=)` structurally CANNOT self-derive — `_resolve_basis_shape` raises
  without explicit. Same value from ONE source (`ng=mix.ng`) → no drift habitat. It is the SECOND
  witness to the DELIBERATE space-anonymous-operator design (Mult.domain docstring: "mesh-free
  default"); if a meshless FunctionSpace is ever minted (CP `[P]`) both args collapse — an
  architectural NOTE, not a fix. NB: `OperatorProduct` does NOT override `as_matrix` (uses base
  apply-to-basis); `MatrixInverseOperator.as_matrix` DOES override (needs no basis_shape — reads
  stored LU dim), but K is the PRODUCT so it routes base.

- **`dominant_eigenpair` squareness-ValueError + float-cast BELONG there (scrutiny a PASS).** It is
  now a PUBLIC primitive with DIRECT callers (K-path + tests) → boundary parse (Pattern 4/8) is
  correct, not ceremony. The squareness guard gives a domain ValueError vs LAPACK's opaque
  LinAlgError, pinned by `test_non_square_raises` (has teeth). Redundant when reached via
  direct_eigenvalue (solve(A,F) always square) — the honest cost of guarding a distinct public
  surface. `-> "OperatorSum"` annotation HONEST (scrutiny b PASS): `__sub__` returns
  `OperatorSum[Domain,Codomain]`, so `collision - k_iso` genuinely is one at runtime.

- **NEW reusable tell — CONSUMER-MIGRATION doc blast radius lands on the ABANDONED PRODUCER's
  docstring.** When a carve moves consumer C off producer P onto producer Q, P's OWN docstring, if
  it named C as its exemplar consumer, goes stale OUTSIDE the diff. Here: `direct_eigenvalue`
  docstring (eigenvalue.py:348) still cites `(solve_homogeneous_infinite)` as its use case, but the
  diff severed that link (solver now calls `dominant_eigenpair` via the K-path). Grep the abandoned
  producer's docstring for the migrated consumer's name. A specific direction of L-004.

- **NEW reusable tell — pre-impl scaffolding retired in the BODY but the module docstring left
  naming the never-shipped speculative CLASS.** test_eigenvalue.py:1-26 headlines the file as
  verifying `DirectEigenvalue` (a CamelCase class that NEVER shipped — the primitive shipped as the
  free function `direct_eigenvalue`); the diff removed the `getattr(_eig,"DirectEigenvalue",None)`
  pre-impl tolerance but left the header. Bug habitat: a contributor greps `DirectEigenvalue`, finds
  nothing in orpheus/, and re-adds a class "to match the docstring," re-opening the class-vs-function
  question the carve closed. On any pre-impl→shipped transition, sweep the test module docstring for
  the speculative symbol spelling.

# #226 step-6 TAIL / carve-P5 slice-1 review (reviewed 2026-07-02, over `394d8c0`)

P5 slice-1 = "composition return types": ONE prod change (`ScaledOperator.inverse()`
annotation `"ScaledOperator"`→`"ScaledOperator[Codomain, Domain]"` — the lone
parameter-DROPPING `.inverse()` among the nine advertisers, now carrying the swap
like its siblings; honest-to-runtime since the body returns
`ScaledOperator(1/α, self.op.inverse())` where `op.inverse()` is already `[C,D]`)
+ a new pyright-only pin bank `_composition_algebra_return_type_static_pins` (13
params). Verdict: **SHIP — ZERO code defects, ZERO VIOLATIONS; 1 cosmetic NIT.**
Reusable rulings:

- **The concrete-method-annotation pin vs the TypeGuard-bridge-Protocol-face pin
  are NOT a twin even when they assert the SAME type today.** `two_space.inverse()`
  (535, inside `if invertible(two_space)`) pins the `SupportsInverse.inverse()`
  bridge face (operator.py:872) reached via the TypeGuard; `a_sum.inverse()` (620,
  direct call on concrete `OperatorSum[S,A]`) pins the concrete `OperatorSum.inverse()`
  annotation (operator.py:1153). Both = `LinearOperator[AngularFlux, ScalarFlux]`
  TODAY. Coextensiveness check (my VIOLATION-standard leg 2, exonerating direction):
  tighten `OperatorSum.inverse()`→`GreenOperator` (plausible future) ⟹ 620 reds, 535
  does NOT (bridge face unchanged); tighten the `SupportsInverse` Protocol / the
  `invertible()` TypeGuard ⟹ 535 reds, 620 does NOT. DISTINCT red-surfaces ⟹ two
  different laws sharing a value ⟹ NOT duplication. Reusable discriminator for any
  pin-bank twin flag.

- **Bare-class `assert_type` pins (perm/tensor/product/diag→their class) are NOT
  vacuous.** `assert_type` demands EXACT type match, so widening `PermutationOperator.
  inverse()`'s annotation to `LinearOperator` reds `assert_type(perm.inverse(),
  PermutationOperator)`. Smaller red-surface than a parametrized pin (no type-args to
  drift), but real teeth + earns its line as documentation-of-the-law. Same for the
  `_ForwardT`-binding pins (`green.inverse()→OperatorSum`, `matrix_inv.inverse()→
  LinearOperator`): they guard the mixin's per-sibling `InverseWrapMixin[X]` binding —
  change `X` and the pin reds.

- **13-param fixture-free pin fn is ACCEPTABLE (not a split trigger).** The file's
  idiom is pyright-only `_..._static_pins` fns where every operand is a param (4 and 7
  precedents); each param is self-documenting via its inline type annotation, 3
  `# ──`-divider sub-sections segment the body, and ONE banner tells ONE coherent
  "composition-algebra return types" story — splitting into constructors-bank vs
  inverse-kinds-bank would fracture the banner. Don't demand a split just for arity.

- **Banner/docstring claim-correctness ALL verified against the live tree** (my L1
  discipline): "two inverse kinds ... canonical statement on `OperatorProduct.inverse`"
  ✓ (the WRAP-DELEGATE vs ALGEBRA-CLOSED contrast is at operator.py:1310-1315);
  "loose family face by design" ✓ (matches `OperatorSum.inverse` docstring 1157-1162);
  "#280 seam — no `A.H.inverse()` pin because `_AdjointOperator` deliberately lacks
  `inverse()`" ✓ (confirmed `_AdjointOperator` has only apply/apply_transpose, AND the
  `InverseWrapMixin` docstring 1623-1626 records the adjoint-inverse as the deferred
  #280 family). All 19 pins match their live annotations.

- **The 1 NIT (cosmetic, non-blocking): pin-fn docstring "Two-space operands ...
  throughout" is literally false for 7/13 params** (endo is endo; identity single-space;
  perm/tensor/diag/green/matrix_inv bare). NOT a bug habitat — `assert_type` reds on any
  mismatch regardless of operand arity, so the imprecision gives no false comfort that
  could admit a swap bug. Offered a one-line reword; did NOT inflate to a fold.

## Step-5b findings delivered

1. **SHOULD-FIX (in-file doc): test_eigenvalue.py:1-26** stale `DirectEigenvalue` ×4 + omits
   `dominant_eigenpair` (the primitive this diff centers with a new TestDominantEigenpair class).
2. **SHOULD-FIX (blast-radius doc): eigenvalue.py:348** `direct_eigenvalue` cites
   `(solve_homogeneous_infinite)` as consumer — the diff severed that link.
3. **POLISH (≤5-line, risk-free): solver.py:196-203** trim "one LU backsolve per column" —
   falsifiable by a future batched `OperatorProduct.as_matrix` override (the diff's OWN
   apply-liveness test contemplates exactly that count change). First comment block (183-191)
   EARNS its lines (loss.inverse()→Green strategy rationale, Rule 3).
4. **PASS callouts:** first-production-consumer claim VERIFIED true (`solver.py:194` is the only
   `MatrixInverseOperator(` construction in orpheus/); old `_assemble_loss_matrix` name fully
   retired (rename, not leave-both); rename matrix→operator follows the return-type change;
   Mode-11 spy hard-`_require` retires the silent-skip hazard; K-resolvent gate procedurally
   independent (`np.linalg.solve` vs operator algebra); principled re-baseline P4-D→5b with all 3
   vv criteria stated. Not re-run (brief: 848 green, ratchet 145).
