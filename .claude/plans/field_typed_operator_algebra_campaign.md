# Field-typed operator algebra — campaign plan (GH #256)

**STATUS: RE-SCOPED into #257 (2026-06-19). Steps 1+3 are the FOUNDATION; the
rest folds into the larger CoefficientField/operator-as-promotion thread.**

⭐⭐ **RE-SCOPE (user architecture review):** the user pushed back that the bare
`Generic[V]` is unexplanatory ("if type hinting doesn't help, it's compliance
theatrics") and that cross-sections should be fields. This surfaced that the
grand report §5.5–5.7 ALREADY specifies the target: **`CoefficientField`** (Field
sibling) promoted to **`MultiplicationOperator`** (`C = M[σ_t]`), `IntegralKernelOperator`
(S/F from Kernels), the Operator/Kernel/Functional suffix law, and the honest
**`TransportState`** bound (refines `Vector`; `np.ndarray` does NOT satisfy it) that
makes `apply(x: TransportState) -> TransportState` read like the domain. Filed as
**#257** (the new home). Frame hardened: `.claude/agent-memory/cross-domain-attacker/coefficient_field_promotion_frames.md`
(multiplier-algebra embedding M: L^∞→B(L²); coefficient space = cone/simplex not
torsor; fission = `M_χ∘ProductionRateFunctional∘M_νΣf`; rename §4158 `.multiplication_operator()`
to avoid the §5.7 class collision). **#256 steps 1 (`cfb651b`) + 3 (`41a92cb`) STAND as
the foundation; steps 4/5/6 (scipy single-source, Functional=§5.6 suffix, BoundaryMomentField)
re-home into #257's implementation plan.** NEXT = explorer blast-radius census → #257 plan.

---

**(historical) STATUS: EXECUTING on branch `feature/field-typed-operator-algebra`.**

⭐ **PICKUP DECISIONS (2026-06-19, settled with user):**
- **Branch base:** merged the 15 unpushed foundation+pyright commits to local
  `main` (ff to `05fa1ef`, NOT pushed — those 7 baseline reds below were already
  on main), branched `feature/field-typed-operator-algebra` off it.
- **Fork #2 → DEFER P3.6** (step 2): lean on the structural `Vector` bound alone
  (TimedFullField satisfies it today). Steps 1/3/4/5 land + clear the pyright
  errors without the `Field` promotion; revisit P3.6 when kinetics `flux⊕precursors`
  gives the genuine 2nd direct-sum instance.
- **Fork #1 → BOTH (general moment-tail)** for step 6: #251 already ships the
  transverse-spatial tail, so adding the angular-moment need is the 2nd instance
  that justifies a general moment-tail layout (not premature).
- **Fork #3 (S/F singledispatch) → decide at step 3** with `scattering.py`/`fission.py` open.

⭐ **PROGRESS:**
- **Step 1 ✅ DONE — `cfb651b`** `refactor(numerics)`. `orpheus/numerics/vector.py`:
  structural `Vector` Protocol + `V` TypeVar, exported. ⚠ **Two corrections vs the
  Q1 memo's recommended `{__add__,__sub__,__mul__}`** (found by auditing the actual
  call sites): the algebra does scalar-mul ONLY as `scalar*v` (`__rmul__`, never
  `v*scalar`), AND divides the carrier (`__truediv__`) at `iteration.py:1021` `Fψ/k`
  + `eigenvalue.py:226` `ψ/p`. Honest contract = `{__add__,__sub__,__rmul__,__truediv__}`.
  Gate `tests/numerics/test_vector_protocol.py` (8 pass `-O`, pyright clean,
  elegance+qa reviewed — qa mutation-verified the negatives).
- **Step 2 (P3.6) — DEFERRED** per fork #2.
- **Step 3 — IN PROGRESS. ⚠ DESIGN RECONCILIATION (type-forced, verified against pyright error sites):**
  the Q1 memo's "single method-scoped `V`, NOT per-operator generics" is type-theoretically
  UNTENABLE. `InvertibleOperator.apply` overriding `OperatorSum.apply` (`operator.py:858/886`)
  and `OperatorSum.__init__(a,b: LinearOperator)` accepting `StreamingOperator`/`CollisionOperator`
  (`operator.py:813`) can ONLY clear if `OperatorSum`+`InvertibleOperator` share a CLASS-level
  `V=TimedFullField` (a method-scoped `V` narrows contravariantly → still incompatible). So the
  design = **class-level `Generic[V]`/`Protocol[V]` (the [[issue-226-operator-generics-map]] §7
  mechanism) WITH `V` BOUND to `Vector`** (step-1's contribution — answers the Q1 memo's
  "unbounded ceremony" objection: a BOUNDED generic names the contract AND clears the overrides,
  and the existing `# type: ignore[arg-type]` on every Mixin dunder get REMOVED).
  ✅ **Step 3 DONE — `41a92cb`** `refactor(numerics)`. Class-level `Generic[V]`/`Protocol[V]`
  bound to `Vector`; SN ops bind `[TimedFullField]`; `apply(self, x: V, /)` positional-only;
  drivers `Generic[V]` with `q_ext: V`. Operator-algebra pyright fileset **92→74 (−18)** — the
  6 `reportIncompatibleMethodOverride` + the whole q_ext/TFF `reportArgumentType` cascade GONE;
  10 dunder `# type: ignore[arg-type]` REMOVED, zero new ignores. Bit-identical (exactly the 7
  baseline reds, 1739 passed). `vector.py` amended: dunders `-> Vector` → `-> Self` (PEP 673;
  needed for `TimedFullField` to satisfy the bound under the generic). A `TYPE_CHECKING`-only
  `apply` stub on the Mixin lets `self` satisfy `LinearOperator[V]` in the dunders (elegance
  measured it strictly cleaner than per-dunder self-typing — self-typing re-opens the cascade).
  ⚠ **Fork #3 RESOLVED → S/F left PLAIN, NOT bound.** Binding+annotating the singledispatch base
  poisoned the `ScalarFlux→ScalarSourceSink` arm; the 2 `reportIncompatibleMethodOverride` at
  `fission.py:223`/`scattering.py:962` are a PRE-EXISTING pyright `singledispatchmethod[NoReturn]`
  modelling gap (invariant at HEAD), orthogonal to #256.
  ⚠ **Step 3↔4 COUPLING (discovered):** threading `q_ext: V` (which clears the solver.py cascade)
  is INSEPARABLE from the ravel boundary — the residual iteration.py driver-internal widening
  (psi from the untyped `_zeros_like`/`_unravel_like` is `Unknown`) is closed by **step 4** (the
  `Ravellable` protocol + `as_scipy_linop_typed`). NIT-1 (typing `_solve_with_seed`) was tried and
  empirically INCREASES errors (74→77, psi widens upstream from `_zeros_like`) → deferred to step 4.
- **NEXT = Step 4** (single-source scipy serialization + `Ravellable` protocol). This ALSO completes
  the iteration.py driver generic (clears the step-3 ravel-boundary widening). Design = §1 D6 +
  Q6 in the cross-domain-attacker memo + §4 in the explorer code-map.

⚠ **BASELINE REDS (route around in every gate):** `main` carries 7 pre-existing
reds — #250 (5 stale SPHERE bit-identity snapshots) + #232 (2 `mu_y` 2-D-mesh +
1-D-quad test-construction) — plus #212 (`continuous_get` hang at
`test_keff_slab::test_heterogeneous_absolute_keff`). The campaign introduces none
of these; deselect them, don't run all `tests/sn`.

---

Read this file + the two design memos in full before touching code:
- `.claude/agent-memory/cross-domain-attacker/issue_226_container_algebra_design.md` (the structural design — Q1–Q7 with recommendations)
- `.claude/agent-memory/explorer/issue_226_field_typed_algebra_map.md` (the exhaustive code map — census, file:lines, gap analysis, blast radius)
- `.claude/plans/issue_226_operator_generics.md` (the SUPERSEDED `Generic[V]` alternative — kept for its layering analysis only)

Line numbers below were current 2026-06-19 on branch `refactor/pyright-signal`
(base `e5f2b1c`); re-confirm at pickup (`stale: file changed` warnings are
expected).

---

## 0. Thesis (why this campaign exists)

The operator algebra's "ndarray operator strays" are **architecture gaps**, not a
genuine second vector population. `TimedFullField` (`transport/timed_full_field.py:123`)
is ALREADY the cross-method generic container — `bulk: BulkField` ⊕
`boundary: BoundaryField` ⊕ history — and `HarmonicMomentField` IS a `BulkField`,
so the container can carry MOMENTS today (only a `BoundaryMomentField` is missing).
With the container as the vector type:
- `MomentProjection` (`projection.py:241`, `domain`=angular, `codomain`=`SphericalHarmonicSpace`)
  is container-ENDOMORPHIC: `TimedFullField→TimedFullField`, the inner `.bulk` flips
  `AngularFlux↔HarmonicMomentField` at runtime. The apparent V→W **dissolves**.
- Flat `np.ndarray` is a **serialization boundary** (scipy Krylov), NOT a vector
  space in the algebra.

This SUPERSEDES the `Generic[V]`-over-heterogeneous-types proposal. The ~30
operator-algebra pyright errors under #226 are fixed here as a byproduct of the
architecture being correct — NOT typed-around.

## 1. Design decisions (settled)

**D1 — Universal vector type = a structural `Vector` Protocol in `numerics`.**
New `orpheus/numerics/vector.py`:
```python
@runtime_checkable
class Vector(Protocol):
    def __add__(self, other: "Vector") -> "Vector": ...
    def __sub__(self, other: "Vector") -> "Vector": ...
    def __mul__(self, scalar: float) -> "Vector": ...
V = TypeVar("V", bound=Vector)
```
`apply(self, x: V) -> V` replaces `apply(x: np.ndarray) -> np.ndarray`
(`operator.py:317`). `np.ndarray`, every `Field` leaf, and `TimedFullField` satisfy
`Vector` via their existing dunders — this just NAMES the duck-type the codebase
already relies on (`_is_ravellable`, `iteration.py:178`). REJECT per-operator
`Generic[F]` (forces every composer generic; cannot express "same container,
different inner leaf").

**D2 — The algebra is endomorphic; three carve-outs.**
`(L+C−S−F−B)` is `V→V` (forced by the iteration `ψ=L⁻¹(Σgᵢψ+q)` — everything
addable ⇒ one container). Carve-outs:
- **Functionals** `V→scalar` (keff/production estimators, `inner_product`/`norm`):
  a SEPARATE `Functional` Protocol `evaluate(x: V) -> float`. NOT `LinearOperator`s
  (capability sets + composer closure are meaningless for a covector).
- **Trace ι\*** `full→boundary`: the internal off-diagonal `A_bs`/`A_sb` block of
  `L`'s FULL block, NOT a public leaf. The boundary law `B` is `V_boundary→V_boundary`
  realised as the `A_ss` block (identity on bulk) ⇒ a `TimedFullField→TimedFullField`
  endomorphism. Already `BlockRole.BULK/FULL/BOUNDARY`-encoded. (Wave O O.4a.2 — keep
  the ι\*/B-extraction DECOUPLED from this campaign; it is a separate behavioral change.)
- **Projection/reconstruction**: same container, different inner leaf — covered by the
  single `apply(x:V)->V` signature (this is WHY D1 beats per-op generics).

**D3 — The 7 matrix primitives stay FLAT + private.**
`PermutationOperator`(961), `IncomingOrdinateMaskTensor`(1065), `PeriodicWrapOperator`(1155),
`DiagonalOperator`(1417), `RankOneOperator`(1535), `TensorProductOperator`(1203),
`SumOfTensorProductsOperator`(1324), + `IdentityOperator`(877)/`ZeroOperator`(899)
(all `numerics/operator.py`). Each acts on ONE tagged numpy axis of a field's
`.values` (ordinate `axis=0` / group axis) and is already `& IdentityOperator`-folded
into a `TensorProductOperator`. **Scope rule (the campaign's dividing line):**
*composite state ⇒ container-typed `LinearOperator`; single tensor axis of one leaf's
`.values` ⇒ flat ndarray primitive used INSIDE a typed operator's apply.* They stay
`np.ndarray`-typed (which satisfies the `Vector` bound) — do NOT container-ify them
(that interposes dataclass alloc on the zero-copy hot path). `ZeroOperator.codomain_zero`
(947–955) is already this pattern in miniature.

**D4 — `BoundaryMomentField` + `MomentFullField`.**
`BoundaryMomentField` = a `BoundaryField` subclass (gets mesh/TraceSpace/`face_view`/
factories free; `_bases.py:481`) storing a moment-resolved trace — the boundary-locus
sibling of `HarmonicMomentField`. Storage lever = `geometry.py:boundary_face_layout`
appending a per-face moment tail (the #251 Leg B mechanism; `trace_space.py` UNCHANGED;
metric/selectors broadcast). `MomentFullField` = a named FACTORY `moment_full_field(φ_lm, b_lm)`,
NOT a subclass — `TimedFullField` stays the ONE container (a subclass re-introduces the
{flux,moment}×{SN,CP,MoC} grid, Cardinal Rule 2). Cross-class safety already holds via
the leaf-class gate (`tff_moment + tff_flux` raises — different `.bulk` classes).

**D5 — P3.6 keystone: write `DirectSumSpace`, promote `TimedFullField`→`Field`.**
⚠ `DirectSumSpace` DOES NOT EXIST (named only in docstrings/plans). BUT `FullFieldSpace`
(`numerics/spaces/full_field_space.py`, the bulk⊕trace composite with per-block G-metric,
already built for the adjoint path, referenced `operator.py:351,356`) is its SEED.
P3.6 = (a) install `FullFieldSpace`(=`DirectSumSpace` after the cross-method rename) as
the container's `.space`; (b) promote `TimedFullField` to `Field` via **STRUCTURED
storage** (keep `.bulk`/`.boundary` attributes; relax `Field.values` to `Vector`) — NOT
flat storage (that couples storage to serialization). The container's delegate dunders
already match `Field`'s, so promotion DELETES them. TimedFullField currently LACKS
`values`/`space`/`UNITS`/`l2`/`inner_product` that `Field` demands; route `l2`/`inner_product`
through the FullFieldSpace direct-sum metric. Design doc: `.claude/plans/depth_b_field_on_function_space.md` §3.8.

**D6 — scipy as typed serialization (single-source).**
Today: TWO ravel paths (`as_scipy_linop` `operator.py:1661` with ZERO prod callers +
the inline `_ravel`/`_unravel_like` closure feeding `spla.LinearOperator` at
`iteration.py:755`). Merge into `as_scipy_linop_typed(op, template)`: serialization lives
INSIDE matvec (`flat→Vector→op.apply→flat` via `template.from_flat/to_flat`), invisible to
scipy; `n = template.to_flat().size` (kills the `as_scipy_linop` shape-guess footgun
`1679–1683`). GMRES/BiCGSTAB still see a plain `(n,n)` matvec. Drivers DROP the inline
closure ⇒ operators never touch ndarrays.

## 2. OPEN FORKS — decide at pickup (do not start step 6 until #1 is named)

1. **`BoundaryMomentField` axis.** Per-face ANGULAR-moment trace (boundary twin of
   `HarmonicMomentField`, for windowing the trace into SH moments — the DSA/P_N need)
   vs per-face TRANSVERSE-SPATIAL-moment trace (#251 Leg B face-slope). DIFFERENT axes.
   Name which (likely the angular-moment one for the moment-projection closure; #251 is
   a separate spatial concern).
2. **Two-instance rule for `DirectSumSpace`.** Is moment `bulk⊕boundary` the SECOND
   direct-sum instance (so writing `DirectSumSpace` now is justified), or must kinetics
   `flux⊕precursors` land first? (`feedback_unify_after_two_instances`.) If the latter,
   D5 step 2 waits and the Protocol retype (step 3) leans on the structural `Vector`
   bound alone (TimedFullField satisfies it TODAY, pre-promotion).
3. **S/F `singledispatchmethod` apply** (`scattering.py`, `fission.py` — the biggest
   typing wart; their true (Vin,Vout) is a UNION TFF/AngularFlux/ScalarFlux/ndarray).
   Keep the multi-arm union or split into typed siblings. Decide before step 3 touches them.

## 3. Execution sequence (steps 1–5 BIT-IDENTICAL typing; step 6 the only behavioral change)

Each step: implement → `elegance-enforcer` + `qa` review → CLI pyright re-measure +
regression gate → commit (`feat`/`refactor` + a `chore(claude)` records commit).

1. **`Vector` Protocol** (`numerics/vector.py`). Pure addition, zero runtime risk.
   Gate: import smoke; `isinstance(np.zeros(3), Vector)` / `isinstance(timed_full_field, Vector)` True.
2. **P3.6** (gated on fork #2): install `FullFieldSpace` as `TimedFullField.space`;
   promote to `Field` (structured storage). Gate: `test_g_adjoint_reciprocity` ≤1e-13
   survives; `from_flat(to_flat(x))==x` 0 ULP; delegate-dunder deletion is behavior-neutral.
3. **Protocol retype `apply(x: V) -> V`** across `operator.py` (Protocol + Mixin +
   OperatorSum + composers) + leaf parameterization (SN ops, projection/reconstruction).
   Clears the 7 `reportIncompatibleMethodOverride` + the `q_ext: np.ndarray` cascade
   (only 2 prod sites: `iteration.py:450,697`). Decide fork #3 here.
4. **Single-source the scipy adapter** (`as_scipy_linop_typed`); drop the inline closure.
   Gate: SN Krylov keff unchanged to solver tol; `from_flat(to_flat)` 0 ULP.
5. **`Functional` category** for keff/production/`inner_product` (split out of `LinearOperator`).
   Negative gate: a `Functional` is NOT a `LinearOperator` (no `apply`/capabilities).
6. **`BoundaryMomentField` + moment-state closure** (gated on fork #1) — the ONLY
   return-shape change. ⭐ PROACTIVELY dispatch **test-architect** FIRST (crosses
   scalar↔angular↔moment — the Smell-16 hazard; the `solve_moments`→`(moment_buf, None)`
   boundary-drop is the gap this closes). Gate: moment tensor byte-identical pre/post;
   the new boundary block provably == the old `None` (zero-equivalent).

## 4. Verification (Cardinal Rule 1 — this is a TYPING/architecture refactor; numerics MUST NOT move)

- Full SN regression snapshot + MMS L1 gates GREEN at **0 ULP** through steps 1–5
  (`-W "error::tests.sn.regression._regression_assert.DriftWarning"`).
- `from_flat(to_flat(x)) == x` to 0 ULP (now also for the moment container).
- `test_g_adjoint_reciprocity` ≤ 1e-13 survives P3.6.
- Negative type-gates: `tff_flux + tff_moment` raises; `apply` returns the same container;
  a flat primitive rejects a container; a `Functional` is not a `LinearOperator`.
- Step 6 only: moment tensor byte-identical; new boundary block == old `None`.
- Route-arounds: NEVER all `tests/sn` (#212 `continuous_get` hang at
  `test_keff_slab::test_heterogeneous_absolute_keff`). Use `tests/sn/operators
  tests/sn/spatial tests/sn/sweep/core tests/sn/solve` + the relevant MMS.

## 5. Blast radius (from the explorer map)

- `q_ext: np.ndarray` — 2 prod sites (`iteration.py:450` SourceIteration, `:697` Krylov);
  retype to `Vector`/ravellable, NOT `Field` (TimedFullField isn't a `Field` pre-P3.6).
- Operator annotation surface = the §1 census in the explorer memo (exhaustive file:lines)
  + composer family + the iteration drivers' `L`/`*gains` + estimators.
- Π/R consumers (prod): `sn/scattering.py`, `sn/solver.py`, `sn/operator.py`,
  `sn/loss_representation.py` (windowed-moment path — HIGHEST-value target),
  `numerics/spaces/spherical_harmonic_space.py`, `numerics/basis/spherical_harmonic_basis.py`,
  `transport/fields/harmonic_moment_field.py`.

## 6. Related issues

- **#208** (CLOSED) — landed the substrate: block-role typing, `FullFieldSpace`, G-adjoint,
  typed Source/Residual leaves, box-7 `evaluate_residual` (`sn/solver.py:225`). Build on it.
- **#2** (OPEN, DSA) — consumes the typed residual; downstream beneficiary.
- **#200** (OPEN, block-inverse Krylov precond) — lands at the `M_scipy` precond closure
  (`iteration.py:760-769`); downstream.
- **#226** (OPEN, pyright) — the operator-algebra errors are tracked-deferred to THIS campaign.
- **#251** — the trace-moment-axis storage lever (shared with fork #1).
- **TransportState** = renamed to `TimedFullField` (2026-05-28); no separate class.

## 7. Execution discipline

New branch off `main` (e.g. `feature/field-typed-operator-algebra`) — NOT the
`refactor/pyright-signal` branch (that is the orthogonal pyright pass). EXPLICIT paths
only; never `.claude/skills/*`, `CLAUDE.md`/`.claude/rules`/`.claude/hooks`, `docs/_build/`,
or the forbidden untracked set. Forbid `git checkout/restore/stash/reset` on tracked paths
in every sub-agent brief (L28). Commit only when asked; push needs an explicit ask. Trailer:
`Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`. After step 6 the
archivist documents the container algebra + the moment-boundary closure; rebuild Sphinx clean.
