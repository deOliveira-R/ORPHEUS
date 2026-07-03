# Verification plan — #226 B4 *inverse-as-operator* carve (Phases 2–5)

**Sibling to** `.claude/plans/issue_226_b4_operator_generics_verification.md` (the
Phase 0/1 spec — its GAP-1/2/3, M1–M5, §1 contract-pinning inventory are LANDED
on main, `ab75d02`/`272caa9`). This spec **EXTENDS** that net for Phases 2–5; it
does NOT re-derive the Phase 0/1 contract-pinning inventory — read that file's §1
for the leaf capability / block_role / `.H`-reciprocity / matvec-twin tables. The
carve itself is `.claude/plans/operator_inverse_algebra_carve.md` (§3 target, §5
phases, §6 forks).

**Carve scope (Phases 2–5, verbatim from the carve plan).**
- **P2** `.inverse()` returns an OPERATOR: `SweepInverseOperator` (wraps the
  existing `(L+C)` sweep, for `InvertibleOperator` at
  `orpheus/sn/operators/streaming.py:509`), `KrylovInverseOperator` (GMRES on
  `A.apply`, for general sums) — both `LinearOperator[C, D]` whose
  `.apply(b, *, initial_guess=None)` returns `A⁻¹b`. Coexists with `.solve()`.
- **P3** migrate consumers `.solve()` → `.inverse().apply()`; DELETE the
  `_solve_accepts_seed` `inspect.signature` probe (iteration.py:457–462).
- **P4** RETIRE the ENTIRE `capabilities` frozenset — BOTH axes (`CAP_APPLY` +
  `CAP_SOLVE` + `CAP_APPLY_TRANSPOSE`) + `MissingCapability`; delete the
  `.solve`/`.apply_transpose` `# type: ignore`s. (Scope widened — see §0′; the
  carve-plan's solve-axis-only sketch is superseded by §0.6.)
- **P5** composition algebra in `@`/`+`/`.inverse()` return types + docs.

**Claim layer (§1.5 gate).** This is a *behaviour-preserving* carve — every
numeric output stays **bit-identical** (no FP-reassociation expected). So the
overwhelming majority of gates are **foundation** software-invariant claims
(capability sets, operator-algebra laws, the matvec≡sweep / inverse≡solve twins)
plus the **structural-identity** adjoint-reciprocity gate. The value-bearing
reference pillar is **closed-form** (`k_inf = νΣ_f/Σ_a`, `φ = Q/Σ_t`, dense
`np.linalg.solve`) and **structural identity** (G-reciprocity
`⟨Aψ,φ⟩_G = ⟨ψ, A.Hφ⟩_G`). The L1 keff/MMS net rides ON TOP as the end-to-end
safety. **No gate in this carve makes an eigenvalue claim on an MMS reference**
(anti-#2; MMS proves flux-shape/order only — the keff claims ride the closed-form
k_inf anchor).

**The discipline that governs this carve (§0.5).** A behaviour-preserving carve
is "done" only when, for EVERY gate, a named mutation reddens it under the
canonical invocation. A green run proves nothing; the mutation proves the teeth.

**Canonical invocation** (Host, Python 3.14 — xdist UNSTABLE on tests/sn +
tests/numerics, SERIAL only):

```
.venv/bin/python -O -m pytest <paths> -m "not slow" -q -rfE --timeout=300 \
    -p no:xdist -p no:cacheprovider
```

**Constraints.** Behaviour-preserving / bit-identical; **NO `# type: ignore`**
(principled typing only — anti-#19); CLI `npx pyright orpheus/` is the oracle.

---

## 0′. ⚠ SCOPE CHANGE (user-steered 2026-06-29) — retire the ENTIRE frozenset, BOTH axes

**This SUPERSEDES §1b ("KEEP `CAP_APPLY_TRANSPOSE`") and §10 fork 1.** §1b was the
correct *local* fix (don't break the landed S† "propagate iff both"), but leaving
`CAP_APPLY_TRANSPOSE` as a string while `solve` moves to a Protocol is a **twin-path**
(Cardinal Rule 2): two mechanisms for the identical structural question. The user
FORBADE it. The locked mechanism (`operator_inverse_algebra_carve.md` §0.6) retires
`capabilities: frozenset[str]` ENTIRELY — BOTH axes — via ONE shape:

| Concern | Solve axis | Adjoint axis |
|---|---|---|
| operator-returning method | `.inverse() -> SweepOperator` | `.H -> _AdjointOperator` (exists) |
| static contract Protocol | `SupportsInverse` | `SupportsAdjoint` |
| runtime instance QUERY | derived `@property is_invertible` | derived `@property is_adjointable` |
| undefined leaf | eager raise `NotInvertible` at `.inverse()` | eager raise `MissingAdjoint` at `.H` |
| "propagate iff both" | METHOD BODY (`(AB)⁻¹=B⁻¹A⁻¹`, op.py:858) + recursive `is_invertible and` | METHOD BODY (`(A+B)ᵀ=Aᵀ+Bᵀ`, op.py:787) + recursive `is_adjointable and` |

`apply` is the universal base method ⇒ **`CAP_APPLY` retires too.** VERIFIED on HEAD:
the "iff both" LAW lives in the method bodies (`OperatorSum.apply_transpose` op.py:787-788
= `a.apply_transpose(x)+b.apply_transpose(x)`; `OperatorProduct.solve` :858-860;
`.apply_transpose` :862-864; `ScaledOperator.apply_transpose` :921-922) — the frozenset
is a redundant, stringly-typed ADVERTISEMENT, not the law.

**S† is SAFE iff the faithfulness keystone gates it (§2.3).** The reciprocity LAW is the
method body (`a.H + b.H`), UNTOUCHED; only the advertisement moves (`CAP_APPLY_TRANSPOSE
∈ caps` → `is_adjointable`). Deleting the frozenset does NOT perturb #276's S† PROVIDED
§2.3 proves `is_adjointable ≡ CAP_APPLY_TRANSPOSE-membership` for EVERY operator BEFORE
the frozenset is deleted. That faithfulness invariant is the load-bearing keystone (§2.3,
DELTA #5).

**Runtime-query verdict (settles the §0.6 sub-decision): derived `@property`, NOT
`isinstance(SupportsInverse/SupportsAdjoint)`.** See §1d.3 — Protocol `isinstance` is
class-uniform (blind to half-adjointable composites AND the value-dependent singular
edge); the Protocols earn their place as the STATIC contract ONLY. The new adjoint-axis
deltas are: §1d (consumer surface + verdict), §2.3 (faithfulness keystone), §3b
(consumer migration), §4e (full retirement + S†-stays-green gate migration), §7
(M-ADJ-PROP / M-ADJ-FORGE / M-ADJ-EAGER / M-INV-FORGE).

---

## 0. Orientation — the CURRENT tree (verified at branch HEAD `refactor/inverse-as-operator`, off `f903686`)

Facts the Phase 0/1 spec and carve plan PRE-DATE (main absorbed #276's adjoint
campaign + the sn/ reorg + the domain/codomain split — verify before trusting any
older line anchor):

1. **Phase 1 (domain/codomain split) is LANDED.** `LinearOperator(Protocol[Domain,
   Codomain])` (operator.py:283), PEP-696 default `LinearOperator[V] ≡
   LinearOperator[V, V]`. `.H`/`.adjoint()` return `_AdjointOperator` (op.py:615).
   The carve plan's "Phase 1" is DONE; this spec is Phases 2–5 only.
2. **`.inverse()` / `SweepInverseOperator` / `KrylovInverseOperator` DO NOT EXIST
   yet** (grep-confirmed — only a forward-looking comment in
   `transport/operators/scattering.py:1271`). P2 introduces them; the branch has
   **no commits yet** (clean `git diff main...HEAD`).
3. **`.H` is the METRIC Hilbert adjoint, NOT the Euclidean transpose**
   (op.py:668–695): `A.H.apply(y) = G_V⁺ ⊙ apply_transpose(G_W ⊙ y)`. `A.H` (=
   `_AdjointOperator`) advertises `{CAP_APPLY}` **iff** the inner advertises
   `CAP_APPLY_TRANSPOSE`; it does NOT propagate solve, and its own
   `apply_transpose` raises `NotImplementedError`. **This is load-bearing for
   gate #4** — `A.H.inverse()` must implement the swap `(A.H)⁻¹ = (A⁻¹).H`, and
   the transpose-solve flows through `loss_action_transpose` + the G-metric, NOT
   a plain Euclidean `Aᵀ` (lessons L12: `.H` ≠ Euclidean transpose).
4. **`InvertibleOperator` = `(L+C)`** (streaming.py:509, `OperatorSum["FullField"]`),
   caps `{CAP_APPLY, CAP_APPLY_TRANSPOSE, CAP_SOLVE}` (CAP_SOLVE added in
   `__init__` at :665). `apply` → `loss_representation.loss_action(σ,ψ)`;
   `apply_transpose` → `loss_action_transpose(σ,φ)` (the LANDED S† matvec);
   `solve(rhs, *, initial_guess=None)` → the WDD sweep with the Carlson seed
   (:761–825). **`A.inverse()` wraps `A.solve`; `A.H.inverse()` rides
   `A.apply_transpose = loss_action_transpose`.**
5. **The SourceIteration `L` is NOT always an `InvertibleOperator`.** The
   curvilinear / windowed paths pass a `_GaussSeidelResolvent` (sn/solver.py:338,
   advertises `{CAP_APPLY, CAP_SOLVE}`) or `_MomentWindowedResolvent` (:486),
   consumed via `base_resolvent.solve(rhs, initial_guess=psi)` (:2347). **P2's
   `.inverse()` and P3's migration MUST cover the resolvents, not just
   `InvertibleOperator`** (see §3, §4).

---

## 1. ⚠ RETIREMENT-BLOCKER ANALYSIS — the SOLVE axis is wide; the TRANSPOSE axis must NOT be retired (MAIN_AGENT_AUGMENTATION)

**This section is the answer to "does anything on main hard-depend on `CAP_SOLVE`
as a string beyond the `.solve`/SourceIteration path, and would Phase 4 break
LANDED adjoint code?" Read it BEFORE Phase 4.**

### 1a. The `CAP_SOLVE` consumer surface is much wider than the carve plan's sketch

The carve plan §4 frames P4 as "retire `.solve` + `CAP_SOLVE` (the solve axis)"
as if it were SourceIteration + `InvertibleOperator`. It is NOT. `CAP_SOLVE` is
woven through (grep-verified, current tree):

| Consumer | Site | What it does with `CAP_SOLVE` |
|---|---|---|
| `OperatorProduct` | operator.py:839–840, `solve` :858 | `(AB)⁻¹=B⁻¹A⁻¹` propagation iff BOTH |
| `ScaledOperator` | :898–899, `solve` :917 | `(αL)⁻¹=(1/α)L⁻¹` survival |
| `IdentityOperator` | :933, `solve` :940 | trivial |
| `TensorProductOperator` | :1331, `solve` :1368 | caps **intersection**; solve iff all |
| `DiagonalOperator` | :1369 spectrum law, `solve` :1376 | **CAP_SOLVE iff `min|f|>0`** (value-dep) |
| `MultiplicationOperator` (C engine) | `transport/operators/multiplication_operator.py:196`, `solve` :369 | inherits DiagonalOperator's value-dep CAP_SOLVE |
| `_GaussSeidelResolvent` | `sn/solver.py:367,384` | advertises `{CAP_APPLY, CAP_SOLVE}` |
| `_MomentWindowedResolvent` | `sn/solver.py:544` | mirrors base incl. CAP_SOLVE |
| `SourceIteration` | iteration.py:435, `L.solve` :467/468 | construction gate `_has(L, CAP_SOLVE)` |
| `KrylovAcceleration` | iteration.py:712, `L.solve` :713 (`# type: ignore`) | preconditioner SELECTION (`elif _has(L, CAP_SOLVE)`) |
| `numerics/frame.py` | :74 | analysis-face advisory |

**Implication for the migration map (§4):** retiring `CAP_SOLVE` is principled —
"is this operator invertible?" becomes a STRUCTURAL fact (the operator either
has a `.inverse()` method / satisfies `SupportsInverse`, or it doesn't) — but the
blast radius is the WHOLE table, not two sites. Every one of those `_has(op,
CAP_SOLVE)` reads must migrate to `hasattr(op, "inverse")` / `isinstance(op,
SupportsInverse)`, and every `.solve` body to `.inverse().apply`. This is a P4
scope the carve plan under-counts — **flag it to the user before P4 starts.**

### 1b. The SOLVE axis and the TRANSPOSE axis are DIFFERENT — do NOT retire `CAP_APPLY_TRANSPOSE`

The carve plan §4/§1 says "apply/**transpose** universal ⇒ the `# type: ignore`s
on `apply_transpose` delete (not suppress)" and "keep a frozenset ONLY if the
zero-scaling edge needs it" — i.e. it contemplates making `apply_transpose`
universal-on-base and retiring the whole frozenset. **That would break LANDED
#276 adjoint code. The retirement-blocker is `CAP_APPLY_TRANSPOSE`, and it is a
HARD blocker.**

`CAP_APPLY_TRANSPOSE` is NOT trivially universal — `apply_transpose` is a real
partition exactly like solve. Operators that genuinely lack a transpose:
`_AdjointOperator.apply_transpose` RAISES `NotImplementedError` (op.py:697); a
one-directional operator has no Euclidean transpose. The capability gates the
propagation laws that #276's "S† falls out free" rests on (lessons L12):

- `OperatorSum.apply_transpose` propagates **iff both** have `CAP_APPLY_TRANSPOSE`
  (op.py:762–763). `(A+B)ᵀ=Aᵀ+Bᵀ` is only valid when both are transposable.
- `OperatorProduct.apply_transpose` **iff both** (`(AB)ᵀ=BᵀAᵀ`, :841–842).
- `_AdjointOperator.apply` GATES on `_has(inner, CAP_APPLY_TRANSPOSE)` (:669).

**LANDED tests that assert `CAP_APPLY_TRANSPOSE ∈ caps AND CAP_SOLVE ∉ caps`** —
these break if you retire `CAP_APPLY_TRANSPOSE`, and they are the #276 substrate:

- `tests/sn/operators/test_isotropic_scattering.py::test_apply_and_transpose_not_solve`
  (:219), `::test_euclidean_reciprocity` (:120),
  `::test_iso_transpose_matches_dense_per_material` (:133) — K_iso S† crosscheck.
- `tests/sn/operators/test_frame_conjugate_carve.py::test_lambda_apply_and_transpose_not_solve`
  (:170).
- `tests/sn/operators/test_scattering_operator.py::test_capability_set_apply_and_transpose`
  (:170).
- `tests/sn/operators/test_capability_survival.py::test_composite_transpose_follows_closure_law`
  (:139) — the `(L+C) − B` transpose closure (depends on the "iff both" law).
- `tests/sn/operators/test_removal_form_matvec_sweep.py::test_invertible_apply_transpose_is_M_transpose_…_bit_identical`
  (:290) — the S† matvec bit-identity (slab/sphere/cyl).

**RECOMMENDED P4 design (the principled spelling, applied to BOTH axes
separately):**

- **SOLVE axis → RETIRE.** Delete `.solve` methods, `CAP_SOLVE`, the solve-axis
  `MissingCapability` guards. Invertibility becomes the structural `SupportsInverse`
  /`.inverse()`-bearing property. Delete the `.solve` `# type: ignore`s
  (op.py:860, 919, 1376; iteration.py:713) because the call sites move to
  `.inverse().apply`.
- **TRANSPOSE axis → KEEP `CAP_APPLY_TRANSPOSE`.** To delete the `apply_transpose`
  `# type: ignore[attr-defined]`s (op.py:692, 788, 864, 922, 1365) WITHOUT
  retiring the capability: **declare `apply_transpose` on the `LinearOperator`
  Protocol as an optional-presence method** (the same RC4 principled spelling the
  Phase 0/1 spec used for the static contract) while KEEPING `CAP_APPLY_TRANSPOSE`
  as the runtime SoT for the propagation laws. The static type-check passes (the
  method is declared on the base → no `[attr-defined]`), and the "iff both"
  closure laws stay intact. **Do NOT make `apply_transpose` a working
  universal-on-base method** — that dissolves the "iff both" law and silently
  green-lights a wrong adjoint on a half-transposable sum (lessons L12: a
  forgotten/blanket `CAP_APPLY_TRANSPOSE` SILENTLY drops/forges the capability).

Gate this split with §6's canaries (the S† bit-identity + reciprocity + the
`*_not_solve` capability pins) plus M-T (§7).

### 1c. The §6 forks — recommendations (do any change the gate design?)

- **§6.3 (surviving frozenset for zero-scaling).** The STATIC `ScaledOperator(0.0,
  op)` edge is ALREADY handled: the constructor RAISES `ValueError` ("use
  ZeroOperator explicitly", op.py:883–891). The genuine **value-dependent**
  singular edge is `DiagonalOperator`/`MultiplicationOperator` with a zero
  coefficient (`CAP_SOLVE iff min|f|>0`, op.py:1369). **Recommendation:** retire
  the `CAP_SOLVE` frozenset; the singular check SURVIVES as a runtime guard
  *inside* `.inverse()` (raise a clear singular error when `min|f|==0`), NOT as a
  frozenset. **This DOES materially shape gate #7** — it is a positive+negative
  pair on `DiagonalOperator/MultiplicationOperator.inverse()` (§8).
- **§6.4 (keep `solve(A, b)` free fn).** KEEP as sugar for `A.inverse() @ b`
  (grand report line 4631; reads at call sites). Adds ONE smoke gate (`solve(A,b)
  == A.inverse().apply(b)`); otherwise no gate-design change.
- **§6.5 (`InvertibleOperator` concrete vs `SupportsInverse` Protocol).** Does NOT
  change the gates **if gate #3' tests the STRUCTURAL `.inverse()`-bearing
  property** (`hasattr(op, "inverse")` / `isinstance(op, SupportsInverse)`), NOT
  `isinstance(op, InvertibleOperator)`. **Recommendation:** write gate #3'
  structurally so it survives whichever the carve picks; lean concrete-for-now
  (smaller blast radius), Protocol when a 2nd inverse-bearing leaf family appears.

---

## 1d. The ADJOINT-axis consumer surface + the runtime-query verdict (DELTA #3 + augmentation)

### 1d.1 — The `CAP_APPLY_TRANSPOSE` query sites (the parallel to §1a's `CAP_SOLVE` table)

Grep-verified on HEAD `f903686`. **Unlike the SOLVE axis** (§1a — woven through
KrylovAcceleration + the sn/solver resolvents + every composer), the ADJOINT axis
is read **ONLY by composers**. There is **NO external pre-query consumer** — no
KrylovAccel-analogue; `iteration.py`'s `_has` reads are all `CAP_APPLY`/`CAP_SOLVE`
(422/429/435/687/694/712/1003), none transpose. This makes the adjoint migration
SMALLER than the solve axis: there is no preconditioner-selection site to rewire.

| Consumer (production) | Site | Read today | → migrates to |
|---|---|---|---|
| `_AdjointOperator.__init__` caps-swap | op.py:646 | `_has(inner, CAP_APPLY_TRANSPOSE)` → adds `CAP_APPLY` to the adjoint | adjoint exists iff `inner.is_adjointable`; else `.H` eager-raises `MissingAdjoint` |
| `_AdjointOperator.apply` GATE | op.py:669 | `if not _has(inner, CAP_APPLY_TRANSPOSE): raise MissingCapability` | UNREACHABLE under the eager `.H` raise (adjoint only built when adjointable) → DELETE the gate |
| `OperatorSum.__init__` | op.py:762-763 | `_has(a,T) and _has(b,T)` → `caps.add` | recursive `is_adjointable = a.is_adjointable and b.is_adjointable` |
| `OperatorProduct.__init__` | op.py:839/841-842 | `_has(a,·) and _has(b,·)` (solve + transpose) | recursive `is_invertible and` / `is_adjointable and` |
| `ScaledOperator.__init__` | op.py:898/900-901 | `_has(op, ·)` → `survivors.add` | `is_invertible = op.is_invertible`; `is_adjointable = op.is_adjointable` (α≠0) |
| `TensorProductOperator.__init__` | op.py:1333 | `all(_has(op,cap) for op in ops)` | `all(f.is_adjointable for f in factors)` |
| `SumOfTensorProductsOperator.__init__` | op.py:1437 | `all(_has(s,cap) for s in summands)` | recursive `all(...)` |
| `TensorProduct.apply_transpose` GUARD | op.py:1356 | `if CAP_APPLY_TRANSPOSE not in self.capabilities: raise` | `if not self.is_adjointable: raise MissingAdjoint` |
| `SumOfTensorProducts.apply_transpose` GUARD | op.py:1448 | same | same |
| `IncomingSourceOperator.__init__` | op.py:1601, 1819 | caps build incl. transpose/solve | recursive property |
| `SNBoundaryOperator.__init__` | sn/operators/boundary.py:141-143 | `all(CAP_APPLY_TRANSPOSE in law.capabilities for law in laws)` → `caps.add` | `is_adjointable = all(law.is_adjointable for law in laws)` |

**The leaf ADVERTISERS** (the `is_adjointable` ground truth, NOT consumers) — the
class-attr / `default_factory` frozensets become a leaf `is_adjointable` property
(all carry `CAP_APPLY_TRANSPOSE` → return `True`):
`transport/operators/{scattering,fission,isotropic_scattering,multiplication_operator}.py`,
`numerics/{frame,projection}.py`, `sn/operators/{streaming,boundary}.py`, op.py
leaves (Identity/Zero/Dense/dyad/mask). **Two ASYMMETRY fixtures prove the axes are
independent:** `ZeroOperator` is `{CAP_APPLY, CAP_APPLY_TRANSPOSE}` (op.py:213) →
`is_adjointable=True` BUT `is_invertible=False`; `MultiplicationOperator(zero-coeff)`
is the value-dependent twin (adjointable, NOT invertible — §8). The faithfulness
keystone (§2.3) MUST enumerate both or it is blind to the axis SEPARATION.

### 1d.2 — ⚠ LITERAL-STRING reads a `CAP_APPLY_TRANSPOSE`-constant grep MISSES (augmentation flag)

The named #276 gates read the cap as a **bare string** `"apply_transpose"` /
`"solve"`, NOT the constant — the migration audit MUST grep BOTH. Worse, **two are
precondition guards INSIDE the very S† canaries the user named as the "S† stays
green" evidence** — they `AttributeError` (no `.capabilities`) the instant the
frozenset is deleted, unless rewired. "S† stays green" is therefore NOT automatic:

| Gate | Site | Reads | Rewire (REQUIRED for the canary to survive P4) |
|---|---|---|---|
| `test_g_adjoint_reciprocity_full_block` | :210 | `if "apply_transpose" not in A.capabilities: pytest.fail` | `if not A.is_adjointable:` |
| `test_removal_form_matvec_sweep` (S† twin) | :308 | `if "apply_transpose" not in op.capabilities:` | `if not op.is_adjointable:` |
| `test_scattering_operator` | :176 | `"solve" not in op.capabilities` | `not op.is_invertible` |
| `test_fission_operator` | :79 | `"solve" not in op.capabilities` | `not op.is_invertible` |
| `test_bc_universal_invariants` | :419-420 | `"solve"`/`"apply_transpose" not in op.capabilities` | `not is_invertible` / `not is_adjointable` |
| `test_boundary` | :439 | `"apply_transpose" in realized.capabilities` | `realized.is_adjointable` |

(`tests/sn/solve/conftest.py:4` `cap("solve")` is the vv-LEVEL marker, NOT a
frozenset read — LEAVE it.) **No production consumer reads the transpose cap as a
literal string** beyond the composers (grep-confirmed) — the literal-string surface
is test-only, but it is the Phase-3/4 migration trap that the constant-grep hides.

### 1d.3 — VERDICT on the runtime-query form (augmentation sub-decision): derived `@property`, CONFIRMED

The user's lean is correct; `isinstance(op, SupportsInverse/SupportsAdjoint)` is
REFUTED as the *runtime* query, for TWO independent reasons:

1. **Class-uniform on composites.** A `runtime_checkable` Protocol `isinstance`
   checks METHOD PRESENCE at the class level (`hasattr`). `OperatorSum` the CLASS
   defines `apply_transpose` (op.py:787) — so `isinstance(any_sum, SupportsAdjoint)`
   is `True` for EVERY `OperatorSum`, including the half-adjointable
   `matrix_full + matrix_apply_only` whose `apply_transpose` RAISES at call time (a
   summand lacks it). isinstance cannot distinguish adjointable from non-adjointable
   composites — exactly the distinction the pre-query consumer must make.
2. **Blind to the value-dependent edge (§8).** `MultiplicationOperator` the CLASS
   has `solve`/`inverse` — so `isinstance(M_zero, SupportsInverse)` is `True` even
   when `min|f|==0` makes it singular. Only the derived `M.is_invertible` property
   (reading the DATA, `min|f|>0`) is instance-accurate.

⇒ The **derived `@property is_invertible`/`is_adjointable`** (recursive on
composites, data-reading on leaves) is MANDATORY for any pre-query branch. On the
SOLVE axis the genuine pre-query consumer is **KrylovAcceleration's
preconditioner-selection** (`elif _has(L, CAP_SOLVE)`, iteration.py:712 → `elif
L.is_invertible`). The ADJOINT axis has **NO** such consumer (§1d.1) — its
"pre-query" is purely the composer recursion + the eager `.H` raise; the recursive
property is built REGARDLESS because the composers need it to derive their own
`is_adjointable`. EAFP (`try: A.H … except MissingAdjoint`) suffices for a consumer
that just wants to USE the adjoint; the property is for consumers that BRANCH first.

**`runtime_checkable` Protocols STILL earn their place — as the STATIC contract
ONLY.** `SupportsInverse`/`SupportsAdjoint` let pyright verify an operator passed
where one is expected declares `.inverse()`/`apply_transpose`, and give consumers a
clean type to annotate (`def precondition(L: SupportsInverse)`). They are the static
half; the derived property is the runtime-instance half — complementary (matches the
§0.6 table), NOT competing. Do NOT use `isinstance` for a runtime branch.

---

## 2. PHASE 2 KEYSTONE — the equivalence gate (`A.inverse().apply(b) == A.solve(b)`)

**THE proof that `.inverse()` is the same operator as today's `.solve`, in
operator form.** Bit-identity inheritance: `SweepInverseOperator.apply` wraps
`A.solve`, so equality is bit-exact *by construction* — the gate proves the
wrapping introduced ZERO change. (Necessary-NOT-sufficient on its own — it proves
"unchanged", not "was right"; it RIDES the existing closed-form anchors
`test_removal_form_kinf_independent_reference_2g` + keff curvilinear. Lessons L2
two-anchor template.)

**Claim layer:** foundation (bit-identity twin). **Reference:** the legacy
`A.solve` output (inheritance) + the closed-form anchors named above (independence).

**New file:** `tests/sn/operators/test_inverse_operator_equivalence.py`,
`@pytest.mark.foundation`, `@pytest.mark.verifies("inverse-as-operator")` (add the
label as a P5 doc anchor in `discrete_ordinates.rst`).

### Gate 2.1 — InvertibleOperator `(L+C)`, every geometry incl. the curvilinear seed

For each `case ∈ {slab_2g, sphere_2g, cyl_2g}` (≥2G heterogeneous, vacuum BC,
non-zero random rhs `b` — reuse the `_slab/_sphere/_cyl` builders from
`test_removal_form_matvec_sweep.py:147`):

```python
A = L + C                      # InvertibleOperator
b = _random_fullfield(case)    # fixed-seed het ≥2G, non-zero (flat nulls redistribution)
np.testing.assert_array_equal(            # BIT-IDENTICAL — not allclose
    A.inverse().apply(b),
    A.solve(b),
    err_msg=f"{case}: inverse().apply diverged from solve",
)
```

**The curvilinear Carlson coupled-pole seed MUST be exercised** (§0.6 — slab is
BLIND to the seed-threading bug; only sphere/cyl activate the angular
redistribution recursion that reads the seed). Add a seeded variant on
**sphere + cyl** (NOT slab):

```python
seed = _offset_seed(case)      # a FullField seed deliberately ≠ the default Carlson inward sweep
np.testing.assert_array_equal(
    A.inverse().apply(b, initial_guess=seed),
    A.solve(b, initial_guess=seed),
    err_msg=f"{case}: inverse().apply dropped the initial_guess seed",
)
```

**Robustness against a seed-insensitive solve (Mode-11 belt-and-braces).** If the
sphere/cyl solve happens to converge to a seed-independent fixed point, the value
gate above cannot see a dropped seed. PAIR it with an in-process wrap pinning the
PATH (lessons L4/L12):

```python
def test_inverse_threads_initial_guess_to_solve(monkeypatch):
    A = L + C
    captured = {}
    real_solve = A.solve
    def spy(rhs, *, initial_guess=None):
        captured["seed"] = initial_guess
        return real_solve(rhs, initial_guess=initial_guess)
    monkeypatch.setattr(A, "solve", spy)
    seed = _offset_seed("sphere_2g")
    A.inverse().apply(b, initial_guess=seed)
    assert captured["seed"] is seed   # the EXACT object threaded through
```

### Gate 2.2 — the composites (`OperatorProduct`, `ScaledOperator`)

For every op that has BOTH `.solve` and `.inverse()` after P2:

```python
for A in (matrix_full @ matrix_full,           # OperatorProduct (CAP_SOLVE iff both)
          2.0 * matrix_full):                  # ScaledOperator ((αL)⁻¹=(1/α)L⁻¹)
    np.testing.assert_array_equal(A.inverse().apply(b), A.solve(b))
```

Reuse the `matrix_full` fixture (`tests/numerics/test_operator.py:…`). These pin
that `OperatorProduct.inverse()` reverses order (`(AB)⁻¹=B⁻¹A⁻¹`) and
`ScaledOperator.inverse()` divides — the composition-algebra return-type law (P5)
verified at the VALUE level, bit-identical to the legacy `.solve`.

**Mutation (M-EQ):** see §7. Route `SweepInverseOperator.apply` to the FORWARD
`A.apply` instead of `A.solve` → gate 2.1 reds (the forward matvec ≠ the inverse).

### Gate 2.3 — THE FAITHFULNESS KEYSTONE (DELTA #1 + #5; the load-bearing pre-retirement gate)

**Without this gate, Phase 4's frozenset deletion is UNGUARDED.** It proves the
derived predicates mirror the frozenset EXACTLY, for EVERY operator, while BOTH
still exist (Phase 2–3 coexistence) — THEN Phase 4 deletes the frozenset trusting
the equivalence. It is a **SCAFFOLD gate**: it can only run while `op.capabilities`
exists, so it is DELETED in Phase 4 with the frozenset (its job done). The PERMANENT
successors are the recursive-composition pins (below, which never reference the
frozenset) + the rewritten closure-law tests (§4e) + the S† method canaries (§6).

**Claim layer:** foundation (software-invariant — the advertisement-faithfulness
contract). **Reference:** the frozenset itself, for the duration of coexistence — a
pure new-≡-old equivalence (necessary by construction, and SUFFICIENT here because
the frozenset is the LANDED, already-verified SoT this carve preserves bit-for-bit).

**New file:** `tests/numerics/test_capability_faithfulness.py`, `@pytest.mark.foundation`.

**(a) The two-axis faithfulness invariant — enumerate EVERY operator** (leaves +
every composite, slab + curvilinear). Reuse `matrix_full`/`matrix_apply_only`
(test_operator.py:87/96), the SN leaf factories (`IsotropicScattering`,
`ScatteringOperator`, `FissionOperator`, `MultiplicationOperator` both nonzero AND
zero-coeff), and the `_slab/_sphere/_cyl` loss builders:

```python
for op in ALL_OPERATORS:        # Identity, Zero, Dense, dyad, mask, S_iso, S, F,
                                # C(nonzero), C(zero-coeff), L+C, L+C-B, A@B, 2*A,
                                # A.H, TensorProduct, SumOfTP, B_law, ...
    assert op.is_invertible  == (CAP_SOLVE           in op.capabilities), repr(op)
    assert op.is_adjointable == (CAP_APPLY_TRANSPOSE in op.capabilities), repr(op)
```

**(b) The recursive-composition pins (PERMANENT — survive Phase 4).** These relate
the COMPOSITE property to its OPERANDS, never referencing the frozenset, so they
remain the canonical "propagate iff both" gate after the scaffold (a) is deleted:

```python
assert (a + b).is_adjointable == (a.is_adjointable and b.is_adjointable)
assert (a + b).is_invertible  is False                  # a sum is not DIRECTLY invertible (§4c defer)
assert (a @ b).is_invertible  == (a.is_invertible  and b.is_invertible)
assert (a @ b).is_adjointable == (a.is_adjointable and b.is_adjointable)
assert (2.0 * a).is_invertible  == a.is_invertible      # ScaledOperator, α≠0
assert (2.0 * a).is_adjointable == a.is_adjointable
# _AdjointOperator: A.H exists iff a.is_adjointable; (A.H).is_adjointable is False (A.H.H deferred, op.py:697)
```

**Config-blindness (§0.6) — the activating fixtures.** The enumeration MUST include
BOTH a value-dependent **adjointable-but-NOT-invertible** leaf (`ZeroOperator`;
`MultiplicationOperator(true-zero-coeff)` — §8) AND a **half-adjointable composite**
(`matrix_full + matrix_apply_only`, `is_adjointable=False`). A buggy
`is_adjointable` that simply returned `is_invertible` (or vice versa) passes an
all-symmetric (both-True / both-False) suite SILENTLY; only the asymmetric fixtures
break the coincidence. Verify the zero coefficient is a TRUE zero, not `1e-300`
(which `min|f|>0` passes — same trap as §8).

**Teeth (DELTA #5):** M-ADJ-PROP, M-ADJ-FORGE, M-INV-FORGE (§7) all red THIS gate
during coexistence — it is the single net that catches both a broken recursion AND
a forged property. Every other Phase-4 deletion rides on it staying green first.

---

## 3. PHASE 3 — per-consumer migration gates (the existing net pins the behaviour)

P3 rewires `.solve()` consumers to `.inverse().apply()` and DELETES the
`_solve_accepts_seed` `inspect.signature` probe (iteration.py:457–462). No NEW
gates — the migration is pinned by the EXISTING end-to-end net, which must stay
green per consumer. The discipline: migrate ONE consumer, run its pinning gate,
repeat.

| Consumer (migration target) | Pinning gate (must stay green through the rewire) |
|---|---|
| `SourceIteration` (`self.L.solve` → `self._L_inv.apply`) | `tests/sn/solve/test_si_single_primitive_contract.py`; `tests/numerics/test_iteration.py::test_source_iteration_recovers_direct_solve` (:98), `::test_source_iteration_with_explicit_solve_realisation` (:164) |
| …its `L` is a **resolvent** in production | `tests/sn/eigenvalue/test_keff_curvilinear.py` (L1 k_inf het ≥2G; the `_GaussSeidelResolvent` path), MMS slab+2d (`tests/sn/verification/mms/test_mms_ld_slab.py`, `test_mms_ld_2d.py`) |
| `KrylovAcceleration` (preconditioner `L.solve` → `L.inverse().apply`) | `tests/numerics/test_iteration.py::test_krylov_acceleration_recovers_direct_solve` (:200); deletes the iteration.py:713 `# type: ignore` |
| direct `base_resolvent.solve(rhs, initial_guess=)` (sn/solver.py:2347) | `test_keff_curvilinear.py` + `test_prescribed_inflow_consistency.py` |

### 3a. The inspect-probe deletion is provably safe

`grep -rn "_solve_accepts_seed\|_solve_with_seed" tests/` → **no test references
the probe** (verified — a private dispatch detail; the lone `inspect.signature`
hit in tests is `test_compute_psi_half_per_level.py:77`, an UNRELATED
curvilinear-sweep signature check, NOT the SI solve-probe). The probe exists
ONLY because some `L.solve` accept `initial_guess` and some don't; under the new
design the inverse operator's `apply(b, *, initial_guess=None)` takes it
canonically (carve plan §3). **Safety gate:** the EXISTING
`test_source_iteration_with_explicit_solve_realisation` (:164) exercises a custom
`L` whose `.solve`/`.inverse().apply` realises a specific inverse — confirm it
stays green after the probe deletes. **Mutation (M-PROBE):** in the migrated
`SourceIteration`, drop the `initial_guess=psi_prev` kwarg from the
`.inverse().apply` call → the curvilinear `test_keff_curvilinear` Carlson-seeded
path must red (the seed no longer threads). If it stays green, the seed path is
not exercised by the suite → ADD a seeded curvilinear SI case BEFORE deleting the
probe.

### 3b. Adjoint-axis consumer migration (DELTA #3 — composer caps-derivations + the eager `.H` raise)

The §1d.1 sites migrate together — they are the same `_has(op, CAP_APPLY_TRANSPOSE)`
recursion as the solve axis. Discipline mirrors §3: rewire one composer, run its
pinning gate, repeat. The pinning gates are the EXISTING closure-law tests (migrated
in §4e) + the S† canaries (§6, unchanged):

| Migration | Pinning gate (stays green through the rewire) |
|---|---|
| `OperatorSum/Product/Scaled` caps-derivation → recursive `is_adjointable`/`is_invertible` | `test_operator.py` closure pins (§4e), `test_capability_survival::test_composite_transpose_follows_closure_law` |
| `_AdjointOperator.apply` gate (op.py:669) → DELETE (eager `.H` raise makes it unreachable) | `test_g_adjoint_reciprocity_full_block` (after its :210 guard rewires to `is_adjointable`) |
| `SNBoundaryOperator` caps (boundary.py:141) → `all(law.is_adjointable …)` | `test_capability_survival`, `test_boundary.py:439` (rewired) |
| `TensorProduct`/`SumOfTP` caps + `apply_transpose` guards (op.py:1333/1356/1437/1448) | `test_tensor_product_operator.py`, `test_outer_dyad.py`, `test_incoming_ordinate_mask_tensor.py` (§4e) |

**⚠ The eager-`.H`-raise behavior change (LOCKED §0.6: "raise MissingAdjoint at
`.H`").** Today `A.H` on a NON-adjointable A SUCCEEDS — it constructs an
`_AdjointOperator` (op.py:642-657, unconditional) and the raise is LAZY at `.apply`
(the `_has(inner, CAP_APPLY_TRANSPOSE)` gate, op.py:669). The locked design raises
EAGERLY at `.H` (checking `self.is_adjointable` in `.adjoint()`/`.H` before
constructing the wrapper), so the op.py:669 gate becomes dead and DELETES. This is a
genuine behavior change — the raise SITE moves from `.apply` to `.H`.

**Migration audit:** grep `tests/` for any assertion that `A.H` / `A.H.apply` RAISES
(the lazy form `wrapper = A.H; pytest.raises(...): wrapper.apply(x)` migrates to
`pytest.raises(MissingAdjoint): A.H`). On HEAD the only `.capabilities`-transpose
reads in the reciprocity/removal gates are `pytest.fail` PRECONDITIONS (they skip,
they don't test the raise — §1d.2), so no lazy-raise assertion is currently broken —
but the implementer MUST re-grep AFTER Phase 2 lands the eager raise, because it is a
real semantic change. **Mutation (M-ADJ-EAGER, §7):** build `(matrix_full +
matrix_apply_only).H` → MUST raise `MissingAdjoint` eagerly; a returned wrapper means
the eager guard is missing.

---

## 4. PHASE 4 — the retirement MIGRATION MAP (delete / rewire / keep)

Apply the project **retirement-means-test-migration** rule (coding-standards.md):
behavioral contract → REWIRE to `.inverse()`; API-smoke (CAP_SOLVE exists) →
DELETE; characterization → KEEP. The NEW invariant replacing "lacks `CAP_SOLVE`"
is **"a singular / apply-only operator has no `.inverse()`"** (structural:
`not hasattr(op, "inverse")` / `not isinstance(op, SupportsInverse)`).

### 4a. The migration table — every test touching `CAP_SOLVE`/`MissingCapability`/`.solve`

| Test (file::name) | Today asserts | Disposition | New shape |
|---|---|---|---|
| `test_operator.py::test_sum_solve_does_not_propagate` (:132) | `CAP_SOLVE ∉ (A+B).caps` | **REWIRE** | `not hasattr(A+B, "inverse")` — a sum is not directly invertible (becomes Krylov-invertible only if the carve adds `KrylovInverseOperator` to sums; see 4c) |
| `::test_product_solve_propagates_with_both` (:158) | product CAP_SOLVE iff both | **REWIRE** | `(A@B).inverse()` exists iff both have `.inverse()`; `.inverse().apply` reverses order |
| `::test_product_solve_drops_when_one_lacks` (:163) | product drops CAP_SOLVE | **REWIRE** | `not hasattr(matrix_full @ matrix_apply_only, "inverse")` |
| `::test_scaled_preserves_all_capabilities` (:178) | scaled keeps CAP_SOLVE | **REWIRE** | `(2*A).inverse()` exists iff `A.inverse()` exists |
| `::test_zero_lacks_solve` (:204) | `CAP_SOLVE ∉ Zero.caps` | **REWIRE** | `not hasattr(ZeroOperator(...), "inverse")` |
| `::test_identity_full_capabilities` (:199) | Identity has CAP_SOLVE | **REWIRE** | `IdentityOperator().inverse()` exists and is the identity |
| `::test_product_solve_reverses_order` (:287) | `(AB)⁻¹=B⁻¹A⁻¹` value | **REWIRE** | same value via `.inverse().apply` |
| `::test_scaled_solve_divides` (:301) | `(αL)⁻¹=(1/α)L⁻¹` value | **REWIRE** | same value via `.inverse().apply` |
| `test_iteration.py::test_source_iteration_requires_solve_on_L` (:493) | `MissingCapability(match="solve")` when L lacks CAP_SOLVE | **REWIRE** | the NEW GAP-3' shape (§4b): apply-only L → SourceIteration raises (no `.inverse()`) |
| `test_capability_survival.py::test_l_plus_c_is_invertible_with_solve` (:127) | `CAP_SOLVE ∈ (L+C).caps` | **REWIRE** | `(L+C).inverse()` exists / `isinstance(L+C, SupportsInverse)` |
| `::test_composite_keeps_apply_drops_solve` (:134) | `CAP_SOLVE ∉ (L+C−B).caps` | **REWIRE** | `not hasattr(L+C−B, "inverse")` |
| `test_diagonal_operator.py` (:117 `MissingCapability` "non-zero coefficient") | zero-coeff Diagonal lacks solve | **REWIRE** | zero-coeff `DiagonalOperator.inverse()` RAISES (the value-dependent edge, §8) |
| `test_multiplication_operator.py` (:377 `MissingCapability`) | singular M lacks solve | **REWIRE** | singular `MultiplicationOperator.inverse()` RAISES |
| `test_operator.py::test_*_rejects_missing_apply_at_composition` (:248–266) | `MissingCapability(match="apply")` | **KEEP** | the APPLY-axis guard is NOT retired (only the SOLVE axis is) |
| `test_iteration.py::test_*_requires apply` (:355,452,466,480) | `MissingCapability(match="apply")` | **KEEP** | apply-axis guard survives |
| `test_tensor_product_operator.py::TestTensorProductRequiresApply` (:124) | TP ctor `MissingCapability` apply | **KEEP** | apply-axis |
| `test_frame.py` (:425 `MissingCapability` "DENSE") | unbuilt dense solve (#275) | **KEEP / REWIRE** | if the frame analysis face migrates to `.inverse()`, rewire the match; else keep |
| `test_isotropic_scattering.py::test_apply_and_transpose_not_solve` (:219) | caps == {apply, transpose} (∌ solve) | **REWIRE the solve half, KEEP transpose** | `{apply, apply_transpose} ⊆ caps AND not hasattr(S, "inverse")` — KEEP `CAP_APPLY_TRANSPOSE ∈ caps` (§1b) |
| `test_frame_conjugate_carve.py::test_lambda_apply_and_transpose_not_solve` (:170) | same | **REWIRE solve half, KEEP transpose** | same |
| `test_scattering_operator.py::test_capability_set_apply_and_transpose` (:170) | apply+transpose | **KEEP** | `CAP_APPLY_TRANSPOSE` survives |

**Pure API-smoke (CAP_SOLVE-exists) tests → DELETE.** Audit for any test whose
SOLE assertion is "the string `CAP_SOLVE` is importable / present in a recognised
set" with no behavioral content — delete with the symbol. (None found in the
inventory above; all `CAP_SOLVE` assertions are behavioral closure-law pins →
rewire. State NONE-FOUND explicitly so a later reader doesn't hunt.)

**Characterization tests → KEEP under `tests/<module>/characterization/`.** None
in the `CAP_SOLVE` inventory are characterization (all are foundation
closure-laws). The matvec≡sweep / inverse≡solve twins stay foundation.

### 4b. GAP-3' — the NEW invariant (replaces Phase 0/1's GAP-3)

Phase 0/1 GAP-3 = "an apply-only leaf stays non-solvable after the base gains a
solve surface" (`CAP_SOLVE ∉ caps` AND `SourceIteration(apply_only)` raises). The
P4 successor:

```python
# tests/numerics/test_operator.py — GAP-3'
def test_apply_only_operator_has_no_inverse(matrix_apply_only):
    """P4: a singular / apply-only operator has NO .inverse() — invertibility
    is a STRUCTURAL property (no CAP_SOLVE frozenset). Structural so it survives
    the §6.5 concrete-vs-Protocol fork."""
    assert not hasattr(matrix_apply_only, "inverse")          # or: not isinstance(.., SupportsInverse)
    with pytest.raises((MissingCapability, AttributeError, TypeError)):
        SourceIteration(matrix_apply_only)        # L not invertible → raises at construction
```

(Replace `MissingCapability` in the `pytest.raises` with whatever P4 names the
not-invertible error — see 4d. Match on a STRING the carve commits to, e.g.
`match="invert"`.)

### 4c. The "sum is Krylov-invertible" decision shapes the sum rows

The carve plan §1 floats "a sum is invertible via Krylov so `.inverse()` exists
there too" (P5 composition algebra). This is a FORK:

- **If `OperatorSum.inverse()` is added** (returns `KrylovInverseOperator(self,
  …)`): then `test_sum_solve_does_not_propagate` REWIRES to "the sum's
  `.inverse()` is a `KrylovInverseOperator`, NOT a `SweepInverseOperator`" — i.e.
  the sum IS invertible, just by a different KIND. Gate it: `(L+C−B).inverse()`
  exists and `type(...)==KrylovInverseOperator`, and `.apply(b)` solves
  `(L+C−B)x=b` to GMRES tol (cross-checked vs a dense `np.linalg.solve` on a small
  het ≥2G case — structurally independent).
- **If NOT** (sums stay non-invertible until a consumer needs it — Pattern 6
  defer): `test_sum_solve_does_not_propagate` REWIRES to `not hasattr(sum,
  "inverse")`.

**Recommendation:** DEFER `OperatorSum.inverse()` (no current consumer — the
within-group loss `(L+C−S−B)` is solved by the resolvent/Krylov DRIVER, not by an
operator `.inverse()`; the production `K = A_loss.inverse() @ F` form is grand-
report-future). Write the sum rows as `not hasattr(sum, "inverse")` now; the
KrylovInverseOperator-on-sum gate lands with its first consumer. **Flag this fork
to the user** — it decides whether P5's "a sum is Krylov-invertible" is built now
or deferred.

### 4d. M4' — prove the migrated guard bites

Phase 0/1's M4 neutered the `_has(L, CAP_SOLVE)` gate. Its P4 successor monkeypatches
the NEW structural guard:

```python
# in-process: make every operator look invertible
import orpheus.numerics.iteration as it
monkeypatch.setattr(it, "_is_invertible", lambda op: True)   # the new SupportsInverse probe
# Run GAP-3' + test_source_iteration_requires_solve_on_L → MUST RED
# (the pytest.raises no longer raises — the apply-only L is wrongly accepted)
```

If P4 keeps a `MissingCapability`-style raise inside `.inverse()` for the
value-dependent edge, M4' ALSO covers §8's guard.

### 4e. FULL retirement, BOTH axes — S† STAYS GREEN (DELTA #2; supersedes §4a's "KEEP transpose" rows)

Phase 4 now deletes `capabilities: frozenset[str]` ENTIRELY — `CAP_APPLY` +
`CAP_SOLVE` + `CAP_APPLY_TRANSPOSE` + the solve/transpose `MissingCapability` guards
— plus ALL `# type: ignore`s on `solve`/`apply_transpose` (the methods are declared
on the base / leaves; universal `apply` needs no cap). The faithfulness keystone
(§2.3 scaffold) is DELETED here too.

**The S†-STAYS-GREEN proof (DELTA #2) — the named gates split into TWO kinds.** The
reciprocity LAW is the METHOD BODY (`a.H + b.H`, op.py:787; `loss_action_transpose`),
UNTOUCHED by the retirement. So the gates that EXERCISE the method stay green; the
gates that ASSERT frozenset MEMBERSHIP rewire:

| Gate | What it touches | Disposition |
|---|---|---|
| `test_removal_form_matvec_sweep::test_invertible_apply_transpose_…_bit_identical` (:290) | exercises `op.apply_transpose` (METHOD) + a `"apply_transpose" not in op.capabilities` precondition at **:308** | **STAYS** (`assert_array_equal` S† bit-id) — REWIRE its :308 guard → `not op.is_adjointable` (else `AttributeError` post-deletion) |
| `test_g_adjoint_reciprocity_full_block` (:200) | exercises `A.H.apply` (METHOD) + a `"apply_transpose" not in A.capabilities` precondition at **:210** | **STAYS** (`rel<1e-12`) — REWIRE :210 → `not A.is_adjointable` |
| `test_wrong_trace_metric_breaks_reciprocity` (:268, L11 control) | exercises the METHOD only | **STAYS UNCHANGED** — MUST still RED on the wrong metric |
| `test_isotropic_scattering::test_apply_and_transpose_not_solve` (:219) | ASSERTS `CAP_APPLY_TRANSPOSE in caps`, `CAP_SOLVE not in caps` | **REWIRE** → `op.is_adjointable` + `not op.is_invertible` / `not hasattr(op,"inverse")`; drop `CAP_APPLY in caps` (apply universal) |
| `test_capability_survival::test_composite_transpose_follows_closure_law` (:139) | ASSERTS `(CAP_APPLY_TRANSPOSE in composite.caps) == all(… for op in (L+C,B))` | **REWIRE** → `composite.is_adjointable == ((L+C).is_adjointable and B.is_adjointable)` (the recursive pin, §2.3b) |
| `test_fission_operator::…not_solve` (:69/:79), `test_scattering_operator::…apply_and_transpose` (:170/:176), `test_legendre_moment_scattering::…not_solve` (:80) | ASSERT membership (literal strings, §1d.2) | **REWIRE** → `is_adjointable` / `not is_invertible` |

**⚠ The first two are the very canaries the user named as "S† stays green"
evidence** — and BOTH carry a frozenset-string PRECONDITION (reciprocity:210,
removal:308) that breaks on deletion. "S† stays green" is NOT automatic: it REQUIRES
rewiring those two preconditions to `is_adjointable` AS PART OF the retirement. Miss
them and the canaries `AttributeError` instead of asserting — a false RED that masks
whether S† actually held.

**The `test_operator.py` closure pins** (the main membership table) migrate per:

| Test (file::name, §-anchor in test_operator.py) | Today | Rewire |
|---|---|---|
| `test_sum_transpose_propagates_with_both` (:138) | `CAP_APPLY_TRANSPOSE in s.caps` | `s.is_adjointable` (s = full+full) |
| `test_sum_transpose_drops_when_one_lacks` (:143) | `CAP_APPLY_TRANSPOSE not in s.caps` | `not s.is_adjointable` (s = full+apply_only — the half-adjointable pin; M-ADJ-PROP target) |
| `test_product_transpose_propagates_with_both` (:168) | `CAP_APPLY_TRANSPOSE in p.caps` | `p.is_adjointable` |
| `test_scaled_preserves_all_capabilities` (:178) | `s.caps == full.caps` | `s.is_invertible==full.is_invertible and s.is_adjointable==full.is_adjointable` |
| `test_scaled_apply_only_stays_apply_only` (:184) | `s.caps == {CAP_APPLY}` | `not s.is_invertible and not s.is_adjointable` |
| `test_identity_full_capabilities` (:199) | `== {APPLY,SOLVE,TRANSPOSE}` | `Identity().is_invertible and Identity().is_adjointable` |
| `test_zero_lacks_solve` (:204) | `SOLVE∉, APPLY∈, TRANSPOSE∈` | `not Zero().is_invertible and Zero().is_adjointable` (the asymmetry fixture) |
| `test_apply_only_operator_is_not_solvable` (:216, GAP-3) | `CAP_SOLVE∉` + SI raises | GAP-3' (§4b): `not is_invertible` + SI raises |
| `test_multiplication_operator.py` (:359/:369) M0 vs M_zero | `M0.caps=={A,T,SOLVE}`, `M_zero.caps=={A,T}` | `M0.is_invertible and M0.is_adjointable`; `M_zero.is_adjointable and not M_zero.is_invertible` (the value-dep asymmetry) |

**`apply`-universal note.** `test_*_apply_propagates` (:127,:153) + the
`MissingCapability(match="apply")` composition guards (:248-267) — with `CAP_APPLY`
retired and `apply` universal-on-base, an operand genuinely lacking `apply` is caught
by the base Protocol / an `AttributeError`, not a frozenset check. **Fork for the
user:** keep a minimal `hasattr(op,"apply")` guard at composition raising a typed
error (eager, matches the Wave-A "fail at composition not call time" contract — just
no longer reading a frozenset) vs let duck-typing handle it. Recommend KEEP the eager
guard; migrate these tests to `pytest.raises((MissingCapability, TypeError,
AttributeError), match="apply")`. The `NoApplyOperator` fixture (test_operator.py:238)
stays as the negative case.

---

## 5. THE ADJOINT-COHERENCE GATE (#4 — load-bearing; pre-builds #280's substrate)

**Why load-bearing:** the carve SUBSUMES the paused adjoint A3/#280 — `A.H.inverse()`
IS the transpose-solve (= `sweep_transpose`). This gate pins that the inverse
carve's `.H` stays coherent with the LANDED #276 adjoint machinery
(`loss_action_transpose`, `test_g_adjoint_reciprocity`), so #280 can later swap
the Krylov-on-transpose for a direct transpose-sweep with the gate already
guarding the value.

**The design decision (recommended).** Implement `A.H.inverse()` NOW via the swap
`(A.H)⁻¹ = (A⁻¹).H`, with the transpose-solve riding the LANDED
`A.apply_transpose = loss_action_transpose`:
- `_AdjointOperator.inverse()` returns `self.inner.inverse().H` (the swap, by
  construction).
- `A.inverse().H.apply(b) = G_V⁺ ⊙ A.inverse().apply_transpose(G_W ⊙ b)`, where
  `SweepInverseOperator.apply_transpose(b)` solves `Aᵀx=b` via **GMRES on
  `A.apply_transpose`** (= `loss_action_transpose`, already on main). #280 later
  replaces that Krylov-on-transpose with the direct `sweep_transpose` — SAME
  value, gate stays green.

**If the carve DEFERS the transpose-solve to #280** (`A.H.inverse()` raises): mark
the gate `@pytest.mark.xfail(strict=False, reason="#280 sweep_transpose / A.H.inverse()
transpose-solve not yet landed")` so it FLIPS to xpass when #280 lands (lessons
L4 xfail-strict-False discipline — NOT strict=True, which a stale failure would
satisfy for the wrong reason).

**Claim layer:** structural-identity (operator-algebra law) + foundation
self-consistency. **Reference:** the G-reciprocity bilinear identity (closed-form,
structurally independent of the transpose-solve algorithm) + the LANDED
`test_g_adjoint_reciprocity` metric.

**New file:** `tests/sn/operators/test_inverse_adjoint_coherence.py`,
`@pytest.mark.foundation`. Cases: slab_2g, sphere_2g, cyl_2g (≥2G het; the
curvilinear metric is the non-degenerate one — slab's trace metric is simpler).

### Gate 5.1 — the LOAD-BEARING reciprocity pin (structurally independent)

For `x = A.H.inverse().apply(b)` (so `A.H x = b`), verify via the reciprocity
identity WITHOUT recomputing the transpose-solve the same way — use the FORWARD
matvec `A.apply` (the sweep-form `loss_action`) + the G-metric:

```python
# ⟨A ψ, x⟩_G == ⟨ψ, b⟩_G  for random ψ   ⟺   A.H x = b
A = L + C
x = A.H.inverse().apply(b)
for psi in random_fullfields(n=3, case=case):
    lhs = g_inner(A.apply(psi), x, space=A.codomain)     # forward matvec + metric
    rhs = g_inner(psi, b, space=A.domain)
    np.testing.assert_allclose(lhs, rhs, rtol=1e-9,
        err_msg=f"{case}: A.H.inverse() is not the metric transpose-solve")
```

`A.apply` (forward sweep-form matvec) is structurally INDEPENDENT of how `x` was
computed (transpose-solve via `loss_action_transpose`/Krylov), and the `g_inner`
reuses the EXACT metric `test_g_adjoint_reciprocity` validates. This is the
genuine cross-check — NOT a round-trip tautology.

### Gate 5.2 — the swap-algebra value law

```python
np.testing.assert_allclose(
    A.H.inverse().apply(b),
    A.inverse().H.apply(b),
    rtol=1e-12,                      # both go through G-metric FP → allclose, not array_equal
    err_msg=f"{case}: (A.H)⁻¹ ≠ (A⁻¹).H — swap algebra broken")
```

(If implemented as the swap, this is near-tautological — keep it anyway: it
GUARDS a future refactor that gives `_AdjointOperator.inverse()` an independent
body, and documents the P5 composition-algebra return-type law.)

### Gate 5.3 — consistency with the LANDED transpose matvec (round-trip)

```python
# A.H.apply uses loss_action_transpose (landed); A.H.inverse() is a different algorithm
np.testing.assert_allclose(A.H.apply(A.H.inverse().apply(b)), b, rtol=1e-9,
    err_msg=f"{case}: A.H ∘ (A.H)⁻¹ ≠ I")
```

This round-trips the inverse against the **landed `loss_action_transpose`** (the
forward `A.H.apply`), confirming the inverse-carve's `.H` is consistent with the
#276 work. Necessary-not-sufficient alone (both arms could share a metric bug) →
its independence comes from being PAIRED with gate 5.1 (which never calls
`A.H.apply`).

**Mutations (§7 M-ADJ):** (a) make `_AdjointOperator.inverse()` return
`self.inner.inverse()` WITHOUT the `.H` (drop the swap) → 5.1 + 5.2 red; (b) drop
the `G_V⁺`/`G_W` metric application in the transpose-solve path (use the Euclidean
`Aᵀ` instead of the metric `.H`) → 5.1 reds on sphere/cyl (the non-trivial trace
metric), stays green on a Euclidean-metric slab (the lessons-L12 `.H`≠Euclidean
discriminator — curvilinear ACTIVATES the metric).

---

## 6. BIT-IDENTITY CANARIES THAT MUST NOT MOVE (deliverable 5)

A behaviour-preserving carve moves NO value. These LANDED gates stay
`assert_array_equal` / unchanged-tolerance; ANY drift is a carve BUG, not FP
noise (§1.5 — there is no FP-reassociation in this carve). Name them so the carve
re-runs them per phase:

| Canary | File::name | Stays |
|---|---|---|
| Removal-form matvec twin | `test_removal_form_matvec_sweep.py::test_invertible_apply_is_M_of_C_sigma_bit_identical` (:239) | `assert_array_equal` (slab/sphere/cyl/2d) |
| Removal-form **transpose** matvec twin (S†) | `…::test_invertible_apply_transpose_is_M_transpose_of_C_sigma_bit_identical` (:290) | `assert_array_equal` (slab/sphere/cyl) — the TRANSPOSE-axis canary; reds if P4 perturbs `apply_transpose`/`CAP_APPLY_TRANSPOSE` |
| matvec≡sweep inverse-twin | `…::test_removal_form_matvec_sweep_roundtrip` (:350, sphere excluded) | `assert_allclose` (unchanged tol) |
| k_inf independent ref | `…::test_removal_form_kinf_independent_reference_2g` (:455) | closed-form anchor (unchanged) — gate 2.1's independence source |
| G-adjoint reciprocity | `test_g_adjoint_reciprocity.py::test_g_adjoint_reciprocity_full_block` (:200) | `assert_allclose rtol≤1e-10` |
| …wrong-metric L11 control | `…::test_wrong_trace_metric_breaks_reciprocity` (:268) | MUST STILL RED on the wrong metric (negative control) |
| Scattering S† dense crosscheck | `test_isotropic_scattering.py::test_iso_transpose_matches_dense_per_material` (:133) | per-material dense transpose (unchanged) |
| Scattering Euclidean reciprocity | `…::test_euclidean_reciprocity` (:120) | unchanged |
| keff L1/L2 net | `test_keff_curvilinear.py::test_homogeneous_exact` (≥2G), `TestL2…::test_2g_heterogeneous_fuel_moderator` | tolerances UNCHANGED (a re-typing/rewire cannot move keff) |
| MMS net | `test_mms_ld_slab.py`, `test_mms_ld_2d.py` | flux-shape/order, unchanged |
| Static apply-contract | `test_operators_apply_typed.py` (C6 `assert_type` block, :360–371) | EXTEND, do not weaken (see §6a) |

### 6a. The static `assert_type` gate is the home for the `# type: ignore` deletions

P4 deletes 8 `# type: ignore`s in operator.py + 1 in iteration.py. The PROOF the
deletions are principled (anti-#19, not suppression) is a pyright-checked
`assert_type` block. EXTEND `tests/sn/operators/test_operators_apply_typed.py`
(one home per apply-contract, lessons L10) with the inverse/transpose contracts:

```python
# Phase-4 static contract: apply_transpose is declared on the base (no [attr-defined])
assert_type(A.apply_transpose(b), FullField)          # was a # type: ignore site
assert_type(A.inverse().apply(b), FullField)          # SweepInverseOperator.apply
assert_type(A.H.inverse(), LinearOperator)            # the swap return type
```

**Prove its teeth (CLI-pyright):** mutate one declared overload return (e.g. make
`SweepInverseOperator.apply` return `np.ndarray`) and confirm
`reportAssertTypeFailure` fires under `npx pyright orpheus/ tests/sn/operators/`
(lessons L10). The runtime gate is `pytest`; the static gate is `npx pyright` —
BOTH must be in the per-phase checklist.

---

## 7. MUTATION RECIPES — the inverse axis (prove every NEW gate bites, under `-O`)

In-process monkeypatch ONLY (the carve holds uncommitted edits — NEVER `git
checkout`/`restore`/`stash`, which reverts to HEAD and destroys them;
process-discipline.md). Run each, expect RED under `python -O -m pytest`, restore.

| ID | Target gate | Mutation | Expected RED |
|---|---|---|---|
| **M-EQ** | 2.1 (equivalence keystone) | Wrap `SweepInverseOperator.apply` to route the FORWARD `A.apply` (matvec) instead of `A.solve` | `assert_array_equal` reds (forward ≠ inverse) on every geometry |
| **M-SEED** | 2.1 seeded + 2.x wrap | Make `SweepInverseOperator.apply` ignore `initial_guess` (`apply(b)` only) | sphere/cyl seeded value gate reds IF seed-sensitive; the Mode-11 wrap (`captured["seed"] is seed`) reds UNCONDITIONALLY |
| **M-PROBE** | 3a (probe deletion safety) | In migrated `SourceIteration`, drop `initial_guess=psi_prev` from `.inverse().apply` | `test_keff_curvilinear` Carlson-seeded path reds (seed no longer threads) |
| **M-ADJ-swap** | 5.1 + 5.2 | `_AdjointOperator.inverse()` returns `self.inner.inverse()` (drop the `.H`) | reciprocity 5.1 + swap 5.2 red |
| **M-ADJ-metric** | 5.1 (curvilinear) | Transpose-solve path uses Euclidean `Aᵀ` (skip `G_V⁺`/`G_W`) | 5.1 reds on **sphere/cyl** (non-trivial trace metric), GREEN on slab — the `.H`≠Euclidean discriminator (lessons L12) |
| **M-NOINV** | GAP-3' (4b) | Give an apply-only matrix a no-op `.inverse` attr | GAP-3' reds (`hasattr` True; SourceIteration wrongly accepts) |
| **M4'** | migrated guard (4d) | `monkeypatch _is_invertible → True` | GAP-3' + `test_source_iteration_requires_solve_on_L` red (raise no longer fires) |
| **M-SING** | §8 (value-dependent edge) | `DiagonalOperator(coeff_with_a_zero).inverse()` made to NOT raise | §8 negative gate reds |
| **M-T** | §1b — COEXISTENCE-era proof | Drop `CAP_APPLY_TRANSPOSE` from `OperatorSum`'s `caps.add` (op.py:763) | `test_composite_transpose_follows_closure_law` + S† twin red — proves the frozenset is not inert. **Superseded post-retirement by M-ADJ-PROP** (op.py:763 is gone; the recursive property replaces it) |
| **M-ADJ-PROP** (DELTA #4) | 2.3 keystone + migrated drop-test (§4e :143) + eager `.H` | Break the recursive `and`: `OperatorSum.is_adjointable = a.is_adjointable` (drop `and b.is_adjointable`) | `not (full+apply_only).is_adjointable` reds (wrongly True); 2.3 reds (`True==False`); `(full+apply_only).H` stops eager-raising — proves "propagate iff both" bites on the ADJOINT axis |
| **M-ADJ-FORGE** (DELTA #4) | 2.3 keystone + downstream `.H` | Monkeypatch the half-adjointable `full+apply_only` `is_adjointable → True` (force the lie) | 2.3 reds (property True vs frozenset False, coexistence). Post-retirement: the forge bypasses the eager `.H` raise → `(full+apply_only).H.apply(v)` **AttributeErrors** on the missing `apply_transpose` (a RAISE, not a wrong value — a half-adjointable summand has no transpose to forge) → any `.H`-exercising gate reds — proves `is_adjointable` is LOAD-BEARING |
| **M-ADJ-EAGER** | 3b eager-`.H` raise | `.adjoint()`/`.H` constructs the wrapper WITHOUT the `is_adjointable` precheck | `(full+apply_only).H` returns a wrapper instead of raising `MissingAdjoint` → the eager-raise gate reds |
| **M-INV-FORGE** | 2.3 keystone (solve-axis twin) | `ZeroOperator.is_invertible → True` (force; the value-dep / non-invertible axis) | 2.3 reds on the Zero row (`True` vs `CAP_SOLVE∉`=False) — the solve-axis faithfulness twin, symmetric to M-ADJ-FORGE |

**M-ADJ-FORGE precision note.** The user's framing ("forges a wrong adjoint → a
reciprocity gate reds") resolves to a RAISE, not a wrong VALUE: forcing
`is_adjointable=True` on a *half*-adjointable sum cannot fabricate a transpose for the
summand that lacks one, so the downstream `.H.apply` AttributeErrors (which pytest
counts as a red/error). A wrong-VALUE reciprocity red would require a summand with a
*present-but-wrong* `apply_transpose` — that is a different (value) mutation, covered
by the S† leaf-transpose gates (lessons L12), not by a capability forge. The
capability forge's PRIMARY catcher is the §2.3 faithfulness gate; the downstream raise
is the "why it is load-bearing" consequence.

**Mode-8 hygiene:** every gate above uses `np.testing.assert_*` / `pytest.raises`
/ `assert_type` / collected-module `assert` (pytest-rewritten — fires under `-O`,
lessons L4). NO new gate relies on a PRODUCTION-code bare `assert`. The eager
`NotInvertible`/`MissingAdjoint` raises and the §8 singular guard MUST `raise` (not
`assert`) so they survive `-O` in production. The §2.3 faithfulness `assert`s live in
a collected `tests/` module → rewritten → fire under `-O`.

---

## 8. THE HONEST RESIDUE — the value-dependent singular edge (§6.3 fork; deliverable 7)

**Investigation result.** Two distinct "singular" edges exist; only ONE is
value-dependent and needs a surviving runtime check:

1. **Static `ScaledOperator(0.0, op)`** — ALREADY handled: the constructor raises
   `ValueError` ("use ZeroOperator explicitly", op.py:883–891). The zero scalar is
   known at construction → no runtime check needed; this is NOT the residue. The
   existing `test_scaled_rejects…`/`ValueError` pin stays (KEEP).
2. **Value-dependent `DiagonalOperator(coeff)` / `MultiplicationOperator(f)` with a
   zero entry** — THE genuine residue. `CAP_SOLVE` is advertised iff `min|f|>0`
   (the spectrum law, op.py:1369; the collision diagonal `C` flows through
   `MultiplicationOperator`'s engine). The SAME class is invertible or not
   depending on its DATA — so invertibility CANNOT be a pure type
   (`isinstance(InvertibleOperator)`); a runtime check survives.

**Recommendation (gates the carve must ship).** Retire the `CAP_SOLVE` frozenset;
move the singular check INTO `.inverse()` — `DiagonalOperator.inverse()` /
`MultiplicationOperator.inverse()` raise a clear singular error (the P4-named
not-invertible exception) when `min|f|==0`, and return a working inverse operator
when `min|f|>0`. Positive + negative pair (anti-#11 — a contract-validation
method needs BOTH a must-not-raise and a must-raise):

```python
# tests/numerics/test_diagonal_operator.py  (+ mirror in test_multiplication_operator.py)
def test_nonsingular_diagonal_inverse_applies():
    """POSITIVE: a nonzero-coefficient diagonal IS invertible (min|f|>0)."""
    D = DiagonalOperator(np.array([2.0, 4.0, 0.5]), axis=0)
    b = rng.standard_normal(3)
    np.testing.assert_allclose(D.inverse().apply(b), b / np.array([2.0, 4.0, 0.5]))

def test_singular_diagonal_has_no_working_inverse():
    """NEGATIVE: a zero entry makes the operator non-invertible — .inverse()
    RAISES (the value-dependent edge; min|f|==0). The check is a runtime RAISE
    (NOT a bare assert — survives -O in production)."""
    D = DiagonalOperator(np.array([2.0, 0.0, 0.5]), axis=0)
    with pytest.raises((MissingCapability, ValueError), match="singular|invert|non-zero"):
        D.inverse()                     # or D.inverse().apply(b) if the check is lazy
```

**Eager vs lazy fork (flag to the user):** raise at `.inverse()` construction
(eager — matches the Wave-A "fail at composition not call time" philosophy) vs at
`.inverse().apply` (lazy). **Recommend eager** (consistent with `MissingCapability`
being raised at composition today). The negative gate's `pytest.raises` target
moves accordingly; M-SING (§7) mutates whichever the carve picks.

**Config-blindness note (§0.6):** a diagonal whose coefficient is uniformly
non-zero NULLS this edge — the negative gate MUST use a coefficient with an actual
zero entry (the activating fixture), not a near-zero `1e-300` (which `min|f|>0`
passes). Verify the zero is a TRUE zero at design time.

---

## 9. PRE-CARVE / PER-PHASE CHECKLIST (the net)

**Before P2 (additive):**
1. Confirm the §6 canaries GREEN at branch HEAD (the bit-identity ground inherited
   from main). Capture the pyright baseline `npx pyright --outputjson orpheus/`
   (the 18 B4 reds + the 8+1 `# type: ignore` sites the carve drives to 0/deleted).
2. The Phase 0/1 net (GAP-1/2/3 + M1–M5) is already on main — re-run M1–M5
   in-process, confirm still RED (regression guard the reorg/#276 didn't rot them).

**P2:** write §2 (2.1 all-geometry + seeded sphere/cyl + Mode-11 wrap; 2.2
composites; **2.3 the faithfulness keystone — BOTH axes, every operator, +recursive
pins, +the two asymmetry fixtures**). Add the derived `is_invertible`/`is_adjointable`
properties + `SupportsInverse`/`SupportsAdjoint` Protocols + eager
`NotInvertible`/`MissingAdjoint` (coexisting with the frozenset). Run M-EQ, M-SEED,
M-ADJ-PROP, M-ADJ-FORGE, M-INV-FORGE → RED (2.3 catches the last three). Canaries
GREEN.

**P3:** migrate one consumer at a time; its §3/§3b pinning gate GREEN per step.
SOLVE axis: delete the inspect probe (M-PROBE → RED); rewire KrylovAccel
`_has(L,CAP_SOLVE)` → `L.is_invertible`. ADJOINT axis (§3b): rewire the composer
caps-derivations → recursive properties; switch `.H` to the eager raise (M-ADJ-EAGER
→ RED); delete the op.py:669 lazy gate. Delete the migrated `# type: ignore`s on
solve/apply_transpose; `npx pyright` → those reds gone, no new red.

**P4 (FULL retirement, BOTH axes — §4e):** rewire the §4a + §4e tables
(delete/rewire/keep), INCLUDING the two literal-string canary preconditions
(reciprocity:210, removal:308 → `is_adjointable`) and the `test_operator.py` closure
pins; land GAP-3' + M4' + §8 positive/negative. **DELETE the entire
`capabilities: frozenset` (CAP_APPLY + CAP_SOLVE + CAP_APPLY_TRANSPOSE +
`MissingCapability` solve/transpose guards)** and the §2.3 faithfulness scaffold (its
job done). Delete ALL solve/apply_transpose `# type: ignore`s; EXTEND the
`assert_type` block (§6a, add the `apply_transpose`/`inverse` contracts) + prove
`reportAssertTypeFailure` teeth under CLI-pyright. Canaries (S† transpose twin +
reciprocity + L11 wrong-metric control) GREEN. **The permanent post-retirement
adjoint gate is the §2.3b recursive pins + the migrated closure-law tests + the S†
method canaries** (none reference the frozenset).

**P5:** composition-algebra return types (the sum-inverse fork 4c — defer unless a
consumer); gate 5.x adjoint-coherence (or its xfail if the transpose-solve is
deferred to #280). Sphinx (Cardinal Rule 3): the `inverse-as-operator` design + the
UNIFIED both-axes retirement (`is_invertible`/`is_adjointable` derived from the method
bodies, NO twin-path — §0′/§0.6) + sweep-vs-Krylov as inverse-kinds, into
`discrete_ordinates.rst` "Development history" + the operator-algebra theory page.
Re-baseline the pyright ratchet.

**Done when:** for EVERY gate (2.1, 2.2, **2.3**, 3a, **3b**, GAP-3', M4',
5.1/5.2/5.3, §8±), a named mutation reddens it under `-O` (§7); the §2.3 faithfulness
keystone proved `is_invertible ≡ CAP_SOLVE` AND `is_adjointable ≡
CAP_APPLY_TRANSPOSE` for EVERY operator BEFORE the frozenset deleted; B4 cluster
pyright → 0; **the ENTIRE `capabilities` frozenset RETIRED (both axes), the
solve/apply_transpose `# type: ignore`s DELETED**; every §6 canary bit-identical
(incl. the two rewired literal-string preconditions); the recursive `is_X` pins
GREEN as the permanent successors; Sphinx clean.

---

## 10. Forks the user must steer (surfaced for the surgical-carve checkpoints)

1. **§1b — ✅ RESOLVED by §0′ (user-steered 2026-06-29): RETIRE the whole frozenset,
   BOTH axes.** §1b's "KEEP CAP_APPLY_TRANSPOSE" was the correct *local* fix but a
   twin-path (Cardinal Rule 2). The retirement is SAFE because the "iff both" LAW
   lives in the METHOD BODIES (`a.H+b.H`, op.py:787; `(AB)⁻¹=B⁻¹A⁻¹`, op.py:858),
   NOT the frozenset — the recursive `is_adjointable`/`is_invertible` properties
   replace the advertisement, and the §2.3 faithfulness keystone (proving the
   property ≡ the old frozenset for EVERY operator, before deletion) is the gate
   that licenses it. M-T becomes M-ADJ-PROP. No retirement-blocker remains.
2. **§4c — build `OperatorSum.inverse()` (KrylovInverseOperator-on-sum) now vs
   defer.** Recommend DEFER (no consumer; Pattern 6). Decides the sum-row shape in
   §4a + whether P5's "sum is Krylov-invertible" ships now.
3. **§5 — implement `A.H.inverse()` now (swap + Krylov-on-`loss_action_transpose`)
   vs xfail-defer to #280.** Recommend implement-now (cheap; rides landed
   machinery; makes gate #4 live and genuinely pre-builds #280's substrate). xfail
   (strict=False) is the honest fallback if deferred.
4. **§8 — eager (`.inverse()` construction) vs lazy (`.inverse().apply`) singular
   raise.** Recommend eager (matches Wave-A fail-at-composition).
5. **§6.5 — `InvertibleOperator` concrete subclass vs `SupportsInverse` Protocol.**
   Gate #3' is written structurally to survive either; recommend concrete-for-now.

---

# PART II — the RE-SCOPED §12 steps 1–3 (2026-07-01 extension)

> **Why this Part exists.** Parts §0–§10 above were written for the carve
> plan's ORIGINAL P2–P5 order. The taxonomy plan
> (`operator_machinery_taxonomy.md`) was re-evaluated 2026-07-01 (§17) and
> re-scoped into six steps (§12). This Part adds concrete gate specs for
> **steps 1–3** (the ones landing NOW, main-agent-direct). Parts §0–§10
> REMAIN VALID and map onto the re-scoped order as follows — read this
> crosswalk before treating any older §ref as current:

| §0–§10 section | Re-scoped home | Status |
|---|---|---|
| §2.1/§2.2 equivalence keystone (`InvertibleOperator` + `OperatorProduct`/`ScaledOperator`) | the `InvertibleOperator` slice of **Step 1** | LANDED (`test_inverse_operator_equivalence.py` green @ `69ed531`) |
| §2.3 faithfulness keystone + §4 retirement map + §4e both-axes | **Step 6** (frozenset retirement — the LAST step) | valid, deferred |
| §3/§3b consumer migration | **Step 3** (SOLVE axis) + Step 6 (adjoint axis) | Step 3 = §14 prereq below |
| §5 adjoint-coherence (`A.H.inverse()`) | Step 6 / #280 | valid, deferred |
| §6 canaries, §7 mutation bank, §8 singular edge | ALL steps (§8 = Step 1 negative controls) | valid, extended below |

**⚠ TWO of the brief's Step-1 premises are SUPERSEDED by the live
implementation (L10 refute-the-premise — verified by reading the code being
written concurrently, `operator.py`/`multiplication_operator.py` @ `69ed531`+).**
The implementer did NOT build "reciprocal Diagonals". A new generic
`InverseOperator` (`operator.py:1224`) DELEGATES to the leaf's `.solve`
(the division `b/c`), explicitly rejecting the reciprocal twin for two
reasons in its own docstring (`:1241–1248`): (a) `(1/c)·b ≠ b/c` in FP (a
drift), and (b) `1/Σ` is a *mean free path* — a units-dishonest named
quantity (coding-elegance Pattern 3). Consequences that **rewrite the gate
tolerances** below:

1. **`inverse().apply(b) == solve(b)` is BIT-IDENTICAL for ALL NINE**, not
   ~1 ULP — the design delegates to / composes bit-id-preserving structural
   inverses. (The brief's "reciprocal diagonals are ~1-ULP" premise is moot.)
2. **`(A⁻¹)⁻¹ = A` holds by OBJECT IDENTITY** for the `InverseOperator`-backed
   leaves (`InverseOperator.inverse()` returns `self.inner`, `:1293`), NOT the
   brief's "action-equality at ~ULP". Object identity only degrades to
   action-nulp for the recomputed-scalar composite (`ScaledOperator`, `1/(1/α)`).

The **round-trip** `inverse().apply(A.apply(x)) ≈ x` is still per-type nulp
(the two directions use *different* arithmetic — multiply vs divide — so the
round-trip is NOT bit-exact even when equivalence-to-solve is). This split
(equivalence-to-solve = bit-id; round-trip = nulp) is the load-bearing
distinction of §12 below.

---

## 12. STEP 1 — deliver `.inverse()` on the NINE `is_invertible` advertisers

**Decision (user, locked): DELIVER-ALL, no retirement.** All nine advertisers
gain a working `.inverse()`. The `is_invertible` docstring already PUBLISHES
the `.inverse()` contract (`operator.py:387`) — Step 1 makes every advertiser
honour it.

**Claim layer:** foundation (software-invariant — operator-algebra laws +
bit-id inheritance). **Reference pillar:** closed-form (the analytic inverse
of each structure: `q/c`, `take(q, π⁻¹)`, `q`, `B⁻¹A⁻¹`) + bit-id inheritance
from the LANDED `.solve`. NO eigenvalue/MMS claim here (anti-#2).

### 12.0 The live design + per-type tolerance table (the load-bearing table)

Verified by reading the concurrent implementation (line anchors DRIFT — the
file is being edited live; target BEHAVIOUR, not the line):

| Advertiser | `.inverse()` returns | `apply` realization | `inv.apply == solve`? | **I1 round-trip tol** | `(A⁻¹)⁻¹` form |
|---|---|---|---|---|---|
| `DiagonalOperator` (`op.py`) | `InverseOperator(self)` | `self.inner.solve(x)` = `x/c` | **bit-id** | `nulp(2)` (×then÷) | **`is A`** (object id) |
| `MultiplicationOperator` (`transport/…`) | `InverseOperator(self)` | `inner.solve` = `q/f` (FullField) | **bit-id** | `nulp(2)` | **`is A`** |
| `OperatorProduct` | `OperatorProduct(b.inverse(), a.inverse())` | `b⁻¹(a⁻¹(q))` | **bit-id** (induction) | `nulp(tol_a+tol_b)` | new obj, action-eq |
| `ScaledOperator` | `ScaledOperator(1/α, op.inverse())` | `(1/α)·op⁻¹(q)` | **bit-id** given op bit-id | `nulp(tol_op+2)` | action-`nulp` (`1/(1/α)`) |
| `IdentityOperator` | `IdentityOperator()` | `x` | bit-id | `array_equal` (0 ULP) | action-exact |
| `PermutationOperator` | `PermutationOperator(inverse_perm, axis)` | `take(x, π⁻¹)` | **bit-id** | `array_equal` (0 ULP) | perm-array-equal |
| `TensorProductOperator` | `TensorProduct(f.inverse() ∀ f)` | chained factor inverses | **bit-id** | `nulp(Σ tol_f)` | new obj, action-eq |
| `InvertibleOperator` (`sn/…`) | `SweepOperator(self)` | `self.solve` (WDD sweep) | **bit-id** (§2.1, LANDED) | `rtol=1e-10` (sweep; **sphere excluded**, #282) | deferred (#280) |
| `InverseOperator` (generic) | `self.inner` (object id) | `inner.solve` | — | — | `is inner` |
| `_BoundBoundaryOperator` shim (`geometry/…`) | **Step 1 adds** → `self.inner.inverse()` | inner's | inherits inner | inherits inner | inherits inner |

**Tolerance rule (task A.1 — NEVER a blanket loose tolerance):** the round-trip
tolerance IS the reduction depth of its FP path — *count the arithmetic ops*.
Zero-arithmetic gathers/passthroughs (`Identity`, `Permutation`) are
`array_equal`; one multiply + one divide (`Diagonal`, `Multiplication`) is
`nulp(2)`; composites inherit the sum; the sweep is `rtol=1e-10` (curvilinear
conditioning) with **sphere excluded until #282** (matches the LANDED
`test_removal_form_matvec_sweep_roundtrip` sphere-exclusion — a sphere
round-trip would false-RED on the seed-lag defect, not the carve).

### 12.1 Gate I1 — the universal round-trip (the plan's keystone, §9/§12)

`inv.apply(A.apply(x)) ≈ x` AND `A.apply(inv.apply(q)) ≈ q`, parametrized over
EVERY `(A, A.inverse())` pair, per-type tolerance from §12.0. This is the
STRONGEST Step-1 gate: it needs NO legacy `.solve` (survives the Step-6 solve
retirement) and proves the inverse inverts *its own* forward.

**Structural independence (L2 — round-trip alone is necessary-NOT-sufficient;
both arms share the coefficient).** PAIR each round-trip with a **closed-form
value anchor** (the analytic inverse, structurally independent of the
composition): `D.inverse().apply(q) == q/c`; `P.inverse().apply(q) ==
np.take(q, π⁻¹)`; `I.inverse().apply(q) == q`. The closed form is the "was
right"; the round-trip is the "is self-consistent".

**New file:** `tests/numerics/test_inverse_universal.py`,
`@pytest.mark.foundation`. Parametrize a builder registry of `(name, A,
closed_form_inverse, roundtrip_tol)` rows spanning the numerics family
(`Diagonal`, `Identity`, `Permutation`, `Scaled∘Diagonal`,
`Diagonal @ Permutation` — **non-commuting**, `TensorProduct(Diagonal,
Permutation)`). Fixed-seed non-trivial inputs (a `Diagonal` with **non-unit,
non-uniform** `c` — a uniform or all-ones `c` NULLS the ×÷ round-trip drift
AND the reciprocal mutation, §0.6).

```python
@pytest.mark.parametrize("name,A,closed,tol", _INVERSE_ROWS, ids=lambda r: r)
def test_i1_roundtrip_and_closed_form(name, A, closed, tol):
    x = _nontrivial_input(A, seed=226)
    inv = A.inverse()
    _assert_close(inv.apply(A.apply(x)), x, tol, f"{name}: A⁻¹(Ax) ≠ x")
    _assert_close(A.apply(inv.apply(x)), x, tol, f"{name}: A(A⁻¹x) ≠ x")
    _assert_close(inv.apply(x), closed(x), tol,          # STRUCTURAL anchor
                  f"{name}: A⁻¹ ≠ its closed form")
```

`_assert_close` dispatches `array_equal` (tol==0) vs
`assert_array_almost_equal_nulp(nulp=tol)` vs `assert_allclose(rtol=tol)` per
the row. **The `MultiplicationOperator` I1 lives in
`tests/transport/test_multiplication_operator.py`** (its carrier is a
`FullField`; round-trip on `.bulk.values`, needs an `SNMesh`). **The
`InvertibleOperator` I1 is the EXISTING** `test_removal_form_matvec_sweep_roundtrip`
(§6 canary, sphere excluded) — reference it, do NOT re-add.

### 12.2 Gate I2 — functoriality

`(αA)⁻¹=(1/α)A⁻¹`, `(AB)⁻¹=B⁻¹A⁻¹`, `(A⁻¹)⁻¹≈A` — as ACTION equalities:

```python
def test_i2_scaled():          # (αL)⁻¹ = (1/α)L⁻¹
    _assert_close((α*A).inverse().apply(q),
                  ((1.0/α) * A.inverse().apply(q)), nulp=2, ...)
def test_i2_product():         # (AB)⁻¹ = B⁻¹A⁻¹  — NON-COMMUTING A,B
    _assert_close((A@B).inverse().apply(q),
                  B.inverse().apply(A.inverse().apply(q)), ...)
def test_i2_involution():      # (A⁻¹)⁻¹
    # PER-TYPE assertion form (§12.0):
    assert D.inverse().inverse() is D          # InverseOperator-backed: OBJECT ID
    np.testing.assert_array_equal(P.inverse().inverse().perm, P.perm)   # perm exact
    _assert_close((sc.inverse().inverse()).apply(x), sc.apply(x), nulp=4)  # Scaled: action
```

**⚠ Correction to the brief (task A.2):** for the `InverseOperator`-backed
leaves (`Diagonal`/`Multiplication`) `(A⁻¹)⁻¹ is A` by OBJECT IDENTITY — assert
`is`, NOT action-nulp. The brief's "must be action-equality at ~ULP" is the
recomputed-reciprocal design that was NOT chosen; it applies ONLY to
`ScaledOperator` (`1/(1/α)`). **`(A⁻¹)⁻¹` is NOT testable on
`InvertibleOperator`** (its inverse `SweepOperator.inverse()` is deferred to
#280 — `sweep_operator.py:26`); scope I2-involution to the numerics family.

**The `(AB)⁻¹` row MUST use NON-COMMUTING factors** (`Diagonal @ Permutation`),
else the order-swap mutation M-INV-ORDER (§12.6) is INVISIBLE (config-blindness
§0.6 — a `Diagonal @ Diagonal` product commutes, so `B⁻¹A⁻¹ == A⁻¹B⁻¹` and the
swap cannot red).

### 12.3 Gate I3 — faithfulness continuity of the inverse OBJECT (extend §2.3)

The NEW inverse objects must themselves satisfy the faithfulness keystone
(`is_invertible ≡ CAP_SOLVE`, `is_adjointable ≡ CAP_APPLY_TRANSPOSE`). **Extend
`tests/numerics/test_operator_capability_predicates.py`** (the existing
faithfulness home — it already has `assert_capability_faithful`):

```python
def test_inverse_objects_are_faithful():
    D = DiagonalOperator(_C)                     # nonsingular
    _assert_faithful(D.inverse())                # InverseOperator: {APPLY,SOLVE}, is_inv True, is_adj False
    _assert_faithful((D @ PermutationOperator(...)).inverse())   # OperatorProduct inverse
    _assert_faithful(ScaledOperator(2.0, D).inverse())
    # InverseOperator's own (A⁻¹)⁻¹ = leaf, faithful:
    assert D.inverse().inverse() is D
```

`InverseOperator` advertises `{CAP_APPLY, CAP_SOLVE}`, `is_invertible=True`,
`is_adjointable=False` (`:1261/1290`) → faithful. `SweepOperator` advertises
`{CAP_APPLY}` only → faithful (both False). The existing
`test_operator_capability_predicates.py` predicate net already pins the
FORWARD advertisers; this row adds the INVERSE objects.

### 12.4 Negative controls (§8 is the Diagonal/Mult instance; extend to composites)

Every `.inverse()` is EAGER-guarded (raises at `.inverse()`, matching §8's
recommendation, as implemented): `DiagonalOperator`/`MultiplicationOperator`
raise `MissingCapability` on a true-zero coefficient; `OperatorProduct`/
`ScaledOperator`/`TensorProductOperator` raise when a factor/operand is
non-invertible; `InverseOperator.__init__` raises on a non-invertible inner.
anti-#11 — one POSITIVE + one NEGATIVE per type:

```python
def test_singular_diagonal_inverse_raises():          # NEGATIVE (extends §8)
    with pytest.raises(MissingCapability, match="non-zero|invert|zero entry"):
        DiagonalOperator(np.array([2.0, 0.0, 0.5])).inverse()
def test_product_with_singular_factor_inverse_raises():
    with pytest.raises(MissingCapability, match="invert|both"):
        (DiagonalOperator(_CZ) @ IdentityOperator()).inverse()
def test_nonsingular_diagonal_inverse_applies():      # POSITIVE
    np.testing.assert_array_equal(
        DiagonalOperator(_C).inverse().apply(b), b / _C)   # == solve (bit-id)
```

**§0.6 activating fixture:** the singular coefficient MUST be a TRUE zero
(`0.0`), NOT `1e-300` (which `np.all(c != 0)` passes — same trap §8). **α=0 is
NOT this edge:** `ScaledOperator(0.0, …)` is ctor-rejected (`ValueError`,
`operator.py:996–1004`) — the `ScaledOperator` negative control is a
non-invertible OPERAND, not a zero scalar. The `raise` (not `assert`) survives
`-O` (Mode-8).

### 12.5 The `_BoundBoundaryOperator` shim (the NINTH advertiser)

The shim forwards `capabilities`/`is_invertible`/`is_adjointable`/`apply`/
`apply_transpose`/`block_role` to `self.inner` (`_bound_compat.py:104–134`) but
has NO `.inverse()` yet. Step 1 adds `def inverse(self): return
self.inner.inverse()`. Gate (a geometry/boundary test home,
`tests/sn/…/test_boundary*.py` or `tests/geometry/`):

```python
def test_bound_shim_forwards_inverse():
    inner = <a realized invertible boundary law>
    shim = _BoundBoundaryOperator(inner, kind="reflective")
    # forwards to the SAME inner inverse (object-equal action):
    np.testing.assert_array_equal(
        shim.inverse().apply(q), inner.inverse().apply(q))
    # a non-invertible inner → the shim's .inverse() propagates the raise:
    with pytest.raises(MissingCapability):
        _BoundBoundaryOperator(<apply-only inner>).inverse()
```

### 12.6 Mutation list — each REDs a NAMED gate (under `-O`, monkeypatch only)

| ID | Target gate | Mutation | Expected RED | Config that ACTIVATES it |
|---|---|---|---|---|
| **M-INV-APPLY** | I1 (§12.1) | `InverseOperator.apply` delegates to `inner.apply` (not `inner.solve`) | round-trip `A⁻¹(Ax)=c²x≠x` reds | non-unit `c` (uniform `c=1` blind) |
| **M-INV-ORDER** | I2 product (§12.2) | `OperatorProduct.inverse` returns `a.inverse() @ b.inverse()` (un-reversed) | `(AB)⁻¹` reds | **non-commuting** `Diagonal @ Permutation` (commuting blind) |
| **M-INV-SCALE** | I2 scaled (§12.2) | `ScaledOperator.inverse` uses `α` (not `1/α`) | round-trip `α²·op⁻¹op·x≠x` reds | `α≠1` (α=1 blind) |
| **M-INV-PERM** | I1 permutation (§12.1) | `PermutationOperator.inverse` returns `PermutationOperator(self.perm)` (not `inverse_perm`) | round-trip reds | **non-involution** perm (3-cycle `[1,2,0]`; involution `[1,0,2]` blind) |
| **M-INV-SHIM** | §12.5 | shim `.inverse()` forwards to a WRONG inner (or `self`) | shim forward-equality reds | any invertible inner |
| **M-INV-NOGUARD** | §12.4 | drop the singular guard in `Diagonal.inverse` / `InverseOperator.__init__` | negative control stops raising | TRUE-zero coefficient |
| **M-INV-FAITH** | I3 (§12.3) | `InverseOperator.is_invertible → False` (lie) | `_assert_faithful` reds (`False` vs `CAP_SOLVE∈caps`) | any |

**Every mutation names its activating config** — the mutation-blindness
discipline (§0.6) applies to the MUTATIONS, not just the SUT: a `Diagonal @
Diagonal` (commuting) product cannot red M-INV-ORDER; an involution
permutation cannot red M-INV-PERM; a unit `c` cannot red M-INV-APPLY. The
parametrization registry (§12.1) MUST carry the activating variant of each.

### 12.7 Homes + markers (Step 1)

| Home | Contents | Marker |
|---|---|---|
| **NEW** `tests/numerics/test_inverse_universal.py` | I1 (numerics family) + I2 functoriality | `foundation`, not slow |
| EXTEND `tests/numerics/test_operator_capability_predicates.py` | I3 inverse-object faithfulness | `foundation` |
| EXTEND `tests/numerics/test_diagonal_operator.py` | §12.4 singular-diagonal ± pair | `foundation` |
| EXTEND `tests/transport/test_multiplication_operator.py` | MultiplicationOperator I1 (FullField) + singular pair | `foundation` |
| EXTEND `tests/numerics/test_permutation_operator.py` / `test_tensor_product_operator.py` | per-type I1/I2 if not covered by the universal registry | `foundation` |
| NEW row in a `tests/sn/…/test_boundary*.py` (or `tests/geometry/`) | §12.5 shim forwarding ± | `foundation` |
| EXISTING `tests/sn/operators/test_inverse_operator_equivalence.py` | `InvertibleOperator` equivalence + seed (LANDED) | `foundation` |

All `foundation`, none slow, all `np.testing.*`/`pytest.raises` (fire under
`-O`, Mode-8). **NO** `verifies(<physics-eq>)` — these are software invariants
(L9; foundation MUST NOT carry a physics label).

---

## 13. STEP 2 — reify `M = (L+C−B_lower)` + windowed = `P @ A.inverse()`

**Claim layer:** foundation (machine-precision round-trip; regular-splitting
identity) + a Mode-9-safe fixed-point equivalence. **Reference:** closed-form
(the reified `M.apply` — a genuine forward to round-trip against) + the
existing deforestation oracle. NO eigenvalue/MMS claim.

> **✅ IMPLEMENTED @ `cc293ef` (2026-07-01) — deviations & findings recorded (L10):**
>
> - **Domain finding (re-shapes §13.1's x):** the sweep substrate's shed
>   OVERWRITES the outflow-definition rows (`z.out = streamed`), so the walk
>   realizes `M⁻¹` exactly on the SOURCE SUBSPACE `{y : y.out-rows = 0}` —
>   every production rhs (`q + Sψ + B_upperψ` write bulk/inflow only) —
>   whose M-preimage is the TRACE-CONSISTENT states (`x.out =
>   streamed(x.bulk)`; solve outputs). §13.1 therefore round-trips a
>   CONSISTENT x (built as `LC.solve(random)`) and asserts machine precision
>   on bulk AND trace (STRONGER than the bulk-only falsifier). Discriminator
>   teeth verified intact on consistent x (the confused pairing still REDs
>   O(1) — its apply lacks the subtraction). Family-wide pre-existing
>   property (also `(L+C).solve`); the contract decision for
>   non-source consumers (`M⁻¹A` preconditioning) is **issue #284**.
>   An additive shed (exact on ALL of V) was REJECTED: persistent-buffer
>   consumers (direct SI loop, reconstruction) reuse out-slots across sweeps.
> - **The in-sweep reflect became the row-masked ADDITIVE
>   `SNMaskedBoundaryOperator.reflect_rows_inplace`** (`bf[rows] +=
>   (B·bf)[rows]` — the inhomogeneous row `z_in = y_row + (Bz)_row`); the
>   old whole-face OVERWRITE dropped `y_row` (benign in production where it
>   is 0 on reflective faces, O(1)-wrong as an inverse) and stamped fresh
>   values onto upper rows the splitting defines as LAGGED. Production
>   bit-identity check: the DD regression snapshot drift (6920 ULP /
>   9.81e-13) is IDENTICAL pre/post carve — pre-existing, not this change.
> - **M-SPLIT is UNREPRESENTABLE as spec'd** — the masked `B_lower`
>   single-sources the row split for BOTH `M.apply` and the reflect, so
>   "mask ≠ walk fold" cannot be constructed (prevention-by-construction).
>   Replaced by the pair (one per gate, in-tree permanent mutation tests):
>   **M-SPLIT-DIR** (flip the split direction in `lower_inflow_rows` →
>   flipped rows are read before their reflect → §13.1 REDs) and
>   **M-SPLIT-PART** (doctor one half post-split → §13.2 REDs).
> - **Measured:** W2 round-trip 5.2e-16 bulk / 4.4e-16 trace (was 2.667);
>   split bit-exact; W2-FP G-S≡Jacobi ≤ 50·1e-12 (LS-S4 het vacuum-x/
>   reflective-y box); W1 fused≡deforested 1.8e-16 (bound 4·N·eps=2.1e-14).
> - **§13.4 factor 3 (static codomain):** `FullField → moment-bulk` is NOT
>   statically expressible (the composite is not bulk-generic); pinned as
>   `assert_type(product.apply(rhs), TimedFullField)`
>   (`_windowed_product_static_typing_pins`) + a RUNTIME
>   `isinstance(out.bulk, HarmonicMomentFlux)` in the renamed gate. Minting
>   a bulk-generic composite is a separate typing carve.
> - **Extra pin gained:** `test_prescribed_inflow_consistency` (migrated)
>   now three-way-pins J ≡ GS ≡ K ≤ 1e-11 THROUGH the reified pipeline
>   with a non-zero prescribed-inflow boundary source.
> - Elegance review: SHIP (doc blast radius cleared by archivist in the
>   same pass; SF docstring-scope + cast + twin cross-ref applied).
>   Enforcer note for step 3: a dedicated `P @ M⁻¹` deforestation pin
>   would close the windowed×G-S corner — fold into the step-3 gate design.

### 13.0 Findings that RE-SHAPE the gate configs (L10 — verified against `sn/solver.py`)

- **F-A: the G-S resolvent is d-GENERIC (2-D AND 3-D Cartesian), NOT "2-D
  only".** The brief/§17 phrase "2-D Cartesian, the G-S resolvent's only
  production domain" is imprecise: `_GaussSeidelResolvent`'s docstring
  (`solver.py:358–364`) calls the "2-D ONLY" claim "stale Phase-3 narration";
  `_select_si_resolvent` excludes ONLY 1-D (`and not sn_mesh.is_1d`, `:708`).
  The load-bearing fact for the gate: **G-S needs a multi-D mesh with ≥2
  reflective faces (non-empty `B_lower`); 2-D Cartesian is the cheapest
  activating config**, and 3-D is a valid (already #225-value-gated) domain —
  do NOT write a gate that ASSUMES 2-D-exclusivity.
- **F-B: the WINDOWED path IS 2-D-Cartesian ONLY** (`_maybe_window`, `:620`,
  `if sn_mesh.is_cartesian and sn_mesh.ndim == 2`). Genuinely 2-D-only
  (distinct from F-A) — curvilinear stays full-angular (the Carlson seed the
  moment tensor cannot carry, `:514–523`). The `P @ A.inverse()` deforestation
  gate config MUST be 2-D Cartesian.
- **F-C: the frame identity is `M∘R = 4π·I`, NOT `M∘R = I`.** The W1
  resolution's "coisometry `M∘R=I`" is loose — the LANDED
  `test_frame.py::test_pi_R_is_4pi_identity_through_the_frame` (`:171`) pins
  `analysis ∘ reconstruction = 4π·I` (the no-prefactor SH convention).
  Asserting `=I` would BE the ERR-051 bug (`vv` anti-#11). **Cite the existing
  `4π·I` gate; do NOT mint a `=I` coisometry gate.**

### 13.1 Gate W2-round-trip — the runnable keystone (machine-precision on `M`)

REIFY `M.apply = (L+C).apply − B_lower.apply` (`B_lower` = the existing
boundary `LinearOperator` masked to the schedule's strictly-lower octant
pairs; `B` has no octant-diagonal so `B = B_lower + B_upper` is exact). Then
`M.inverse()` is a genuine `SweepOperator` on `M`:

```python
# 2-D Cartesian, ≥2 reflective faces (non-empty B_lower), ≥2G het
M = (L + C) - B_lower                    # reified forward
x = _random_fullfield_2d(seed=...)
np.testing.assert_allclose(M.inverse().apply(M.apply(x)), x, rtol=1e-10,
    err_msg="M⁻¹(Mx) ≠ x — the reified G-S resolvent does not invert its OWN forward")
```

This is what the confused `_GaussSeidelResolvent` CANNOT pass: §17 falsifier-3
ran it — the clean `(L+C).inverse()` round-trips at **3.8e-16**, the legacy
`_GaussSeidelResolvent` round-trip vs its (vestigial) forward is **2.667, O(1)
RED** (defect = `(L+C−B_lower)⁻¹B_lower x`). The reification is DONE when the
2.667 becomes machine-precision.

### 13.2 Gate W2-split — split exactness `B == B_lower + B_upper`

```python
for face in trace.faces:
    np.testing.assert_array_equal(
        B.apply(psi).face_view(face),
        (B_lower.apply(psi) + B_upper.apply(psi)).face_view(face),
        err_msg="B ≠ B_lower + B_upper — the octant split is not exact")
```

Bit-identical (the split is a partition of octant pairs, no arithmetic
change). **Mutation M-SPLIT:** drop one octant pair from `B_lower`'s mask →
reds (and W2-round-trip also reds — `M` is then the wrong forward).

### 13.3 Gate W2-FP — fixed-point equivalence vs Jacobi (Mode-9-SAFE)

The converged fixed point MUST be splitting-invariant (`vv` Mode 9). **The
config MUST break the degenerate coincidence:** a DIAGONAL cubature
(`lebedev`/`level_symmetric` — shared faces) on a heterogeneous 2-D box, NOT
an axis-aligned `product` quad (which makes octant-G-S accidentally exact,
ERR-056). Assert the G-S-resolvent converged flux == the Jacobi (plain-SI)
converged flux to `SAFETY × conv_tol` (L7):

```python
phi_gs  = solve_via(M_resolvent, ...)      # (L+C−B_lower)⁻¹ G-S splitting
phi_jac = solve_via(plain_LC, ...)         # Jacobi (B external)
np.testing.assert_allclose(phi_gs, phi_jac, rtol=SAFETY*conv_tol, ...)
```

**⚠ the converged-SI-equivalence fallback ALONE is Mode-9-DEGENERATE** (§17
W2): the fixed point is splitting-invariant — it cannot even distinguish G-S
from Jacobi — so W2-FP is NECESSARY but insufficient; **W2-round-trip (§13.1)
is the load-bearing correctness gate**. ρ(M⁻¹N)<1 is a convergence
*certificate*, not correctness — spec it as an optional rate check, not a
correctness gate.

### 13.4 Gate W1 — the windowed product `P @ A.inverse()` (the deforestation rename)

**The existing oracle to rename (task B.4 — FOUND):**
`tests/sn/solve/test_2d_anisotropic_windowing.py::test_2d_windowed_moments_in_sweep_equal_post_projection`
(leg (d), `:256`). Today its SUT is `LC.solve_moments(rhs, frame)` and its
ORACLE is `frame.analysis.apply(LC.solve(rhs).bulk.values)`, bound
`4·N·eps` (scale-relative, principled-equiv — NOT bit-id, NOT element-wise
nulp; the ℓ≥1 block is ~10⁻³× the ℓ=0). **Step 2 migrates the SUT to the
product**: `(P @ A.inverse()).apply(rhs)` where `P` = the scattering frame's
`analysis` face on the bulk ⊕ identity on the trace (`FullField→MomentField`);
the ORACLE stays the post-projection form. The fusion (never materialize ψ)
stays a SUBSTRATE moment-emit (`_SweepEmit` MOMENT mode) — the gate proves the
FUSED product ≡ the post-projection oracle:

```python
# RENAMED: test_windowed_product_equals_post_projection (the deforestation gate)
sut    = (P @ A.inverse()).apply(rhs).bulk.values          # fused, moments
oracle = frame.analysis.apply(A.inverse().apply(rhs).bulk.values)   # deforested
rel_drift = np.abs(sut - oracle).max() / np.abs(oracle).max()
assert rel_drift <= 4 * quad.N * eps        # principled-equiv, scale-relative
```

Three honest-contract FACTORS (§17 W1 — do NOT collapse into one):
1. **fusion-correctness** = the renamed gate above (windowed ≡ post-projection).
2. **coisometry** = the EXISTING `test_pi_R_is_4pi_identity_through_the_frame`
   (`M∘R = 4π·I`, F-C — cite it, do NOT re-mint as `=I`).
3. **codomain typing** = a pyright-visible `assert_type((P @ A.inverse()).codomain-carrier, MomentField)` extend of the `assert_type` block (§6a) — the `FullField→MomentField` type change is a static fact.

**The anti-degeneracy leg is MANDATORY and already present** (leg-2, `:327` —
ℓ≥1 must carry signal, else the moment comparison is vacuous). Keep it. The
config is anisotropic P1 het 2-D (Mode-7 — an isotropic/P0 config is BLIND to
a dropped ℓ≥1 moment, L1/§0.6). **A "moment-proxy residual"
`M[(L+C)Rφ]−Mq` is REJECTED** (§17 W1 — category-confused: `Rφ≠ψ`; it tests a
different P_N-reduced system and would RED a perfect fused sweep).

### 13.5 Pinning gates — the existing behaviour MUST stay bit-identical

Step 2 dissolves the resolvents but keeps their converged behaviour. These
LANDED gates pin it (grep `tests/` for `gauss_seidel`/`window`):

| Pin | File | Stays |
|---|---|---|
| DAG window ≡ full-field (STORAGE oracle, distinct from W1) | `tests/sn/sweep/core/test_sweep_graph_window_equivalence.py` | `array_equal` (bit-id) — substrate-level, Step 2 must NOT perturb it |
| windowed-moments ≡ post-projection | `test_2d_anisotropic_windowing.py::…post_projection` | the W1 gate (§13.4) — migrated to the product |
| G-S / windowed keff+flux end-to-end | `tests/sn/eigenvalue/test_keff_2d.py`, `test_d3_admission.py`, `test_2d_anisotropic_windowing.py` (a,b,c) | tolerances UNCHANGED (a dissolution cannot move keff/flux) |
| SI-≡-Krylov cross-check (structural anchor) | `test_2d_anisotropic_windowing.py::…si_krylov_agree` | `rtol=1e-6` (Krylov stays full-angular → genuine independent ref) |

**Mutation M-DEFOREST:** make the fused product drop the ℓ≥1 accumulation
(emit ℓ=0 only) → the W1 gate + leg (a) SI≡Krylov red on the P1 config (a P0
config stays green — the Mode-7 witness that the anisotropic config is
load-bearing).

### 13.6 Homes (Step 2)

| Home | Contents |
|---|---|
| NEW `tests/sn/…/test_gauss_seidel_reification.py` (or extend `test_solver_components.py`) | W2-round-trip (§13.1), W2-split (§13.2), W2-FP (§13.3) |
| RENAME/migrate in `test_2d_anisotropic_windowing.py` | W1 deforestation gate (§13.4) — SUT → `P @ A.inverse()` |
| EXTEND `test_operators_apply_typed.py` | the `FullField→MomentField` codomain `assert_type` |

W2 gates `foundation`; the windowing gates keep their existing `l1`+`verifies`
markers (they ARE equation-level — anisotropic scattering + 2-D closure).

---

## 14. STEP 3 PREREQUISITE — the fast-lane Mode-11 seed-threading spy

> **✅ IMPLEMENTED @ `1ab7429` (2026-07-02) — step 3 complete (spy + rewire
> + consumer migrations + the A-rename).**
> The spy landed as specced (`tests/sn/solve/test_seed_threading_spy.py`,
> GL-8 het-2G 10-cell VACUUM sphere via `curvilinear_two_region_mesh` +
> mixtures A/B, `max_inner=6` + `inner_tol=0.0` → exactly 6 recorded
> solves, no reconstruction call on the unwindowed path; markers
> foundation+sentinel+cap("solve")) and was GREEN pre-rewire (0.13 s) and
> post-rewire — the §14.1 route-invariance claim held exactly.
> **Teeth:** M-SEED-DROP/ZERO/STALE (throwaway in-process patches on the
> then-extant `_solve_with_seed`) all RED pre-rewire; **M-PROBE**
> (post-rewire realization: sever the `initial_guess` thread at
> `SweepOperator.apply`) REDs the same `seed DROPPED` guard in <1 s — the
> §17 F4 delta closed. §14.5 row 2 delivered:
> `test_si_krylov_eigenvalue_equivalence_sphere` tagged
> `@catches("ERR-026", "M-SEED-DROP")` (M-SEED-VALUE |Δk|≈3.46e-2 stays
> its @slow evidence). §14.3 honored: NO fast sphere value gate.
> **§3 consumer migrations executed** (SourceIteration → `A_inv.apply`
> uniform threading; probe + `_solve_with_seed` + ctor CAP_SOLVE gate
> deleted per §3a; Krylov fallback → `_seeded_inverse(A).apply`, its
> type:ignore deleted; the sn/solver direct site →
> `.inverse().apply(rhs, initial_guess=...)`), all §3 pinning gates green
> incl. the migrated single-primitive contract (SweepOperator + `.inner`
> identity). **Beyond spec:** the enforcer's windowed×G-S corner pin
> (`test_2d_windowed_product_over_gauss_seidel_M_equals_post_projection`,
> `P @ M⁻¹` fused ≡ deforested ≤4·N·eps + B_lower non-degeneracy guard)
> and **the A-rename** (numerics params/attrs/messages/equations
> L→A per the project convention; `SupportsSeededApply` + NEW
> `_seeded_inverse(A)` single-source the narrowing; #285 files the
> structural-vs-convention decision into steps 4–5).
> Tier 2981/0; ratchet 148; Sphinx `-W` clean. The spec text below is the
> as-designed record.

**This gate lands BEFORE the Step-3 driver rewire** (the §12-step-3
prerequisite, §17 F4). It pins the CONTRACT (the per-iterate seed threads),
so it is GREEN on the CURRENT driver and STAYS GREEN through the rewire (it
does not pin the route).

### 14.0 The §17 F4 finding (why the fast net must be the spy)

The keystone `inverse().apply ≡ solve` (§2.1) pins only the WRAPPER
(bit-identical by construction). The Step-3 migration rewires the LOOP surface
— `SourceIteration._solve_with_seed` (`iteration.py:464`) threads `psi_prev`
as `initial_guess` on each inner step (`self.L.solve(rhs,
initial_guess=psi_prev)`, gated by the `_solve_accepts_seed` inspect-probe
`:453–462` that Step 3 DELETES). Empirically (§17 F4): a simulated seed-drop
on het-2G SPHERE gives **|Δk|=3.46e-2**, yet under the canonical `-m "not
slow"` run only a fragile 1G monotone margin reddens — the strong het-2G
sphere value catcher is `@pytest.mark.slow`-DESELECTED, and cylinder's
seeded-value gates are VACUOUS (the per-level α-dome telescopes → the seed
cancels exactly, §14a). So the load-bearing FAST guard is a Mode-11 path-spy.

### 14.1 The spy — pin the per-iterate threading contract

Wrap the route-invariant delegate `InvertibleOperator.solve` (pre-Step-3 the
driver calls `resolvent.solve` = this for the 1-D sphere; post-Step-3 it calls
`L_inv.apply` = `SweepOperator.apply` which DELEGATES to `inner.solve` = this —
so wrapping here survives the rewire). Assert the k-th call receives the
(k−1)-th ITERATE (not `None`, not zeros, not the rhs):

```python
# tests/sn/solve/test_seed_threading_spy.py  (NEW)
def test_si_threads_previous_iterate_as_seed(monkeypatch):
    sn, A, ... = _small_het_2g_sphere()           # §14.2
    calls = []                                     # (initial_guess, return) per solve
    real = InvertibleOperator.solve
    def spy(self, rhs, *, initial_guess=None):
        out = real(self, rhs, initial_guess=initial_guess)
        calls.append((initial_guess, out))
        return out
    monkeypatch.setattr(InvertibleOperator, "solve", spy)

    solve_sn_fixed_source(..., inner_solver="source_iteration", max_inner=6, ...)

    assert len(calls) >= 2, "SI did not iterate — no threading to observe"
    for k in range(1, len(calls)):
        seed_k, _ = calls[k]
        assert seed_k is not None, f"call {k}: seed DROPPED (initial_guess=None)"
        prev_return = calls[k-1][1]
        np.testing.assert_array_equal(          # the seed IS the previous iterate
            seed_k.bulk.values, prev_return.bulk.values,
            err_msg=f"call {k}: seed ≠ previous iterate (zeros/rhs substituted?)")
```

`foundation` + `sentinel`-style (one cheap always-on canary). `-O`-proof
(`np.testing.*`/`assert` in a collected test → pytest-rewritten, fires under
`-O`; but prefer `pytest.fail`/`np.testing` for the `is not None` guard to be
belt-and-braces). It pins the CONTRACT (values-equal, robust to a defensive
copy — do NOT assert `is`, which a copy would false-RED).

### 14.2 The config — seconds, sphere, no NaN

Small het-2G SPHERE where the threading is load-bearing (§17 #6: the sphere
seed-lag is where a drop MOVES the value). **Use GL quadrature (NOT
level-symmetric)** — LS-S4/16-cell fixed-source SI → NaN (§17 #6); GL avoids
the μ=−1 pole NaN. ~8–12 cells, `max_inner=6` (we need ≥2 inner iterations to
observe threading, NOT convergence), moderate `c` (het-2G A/B mixtures). Runs
in << 1 s. Reuse the `_sphere` builder from `test_removal_form_matvec_sweep.py`
with a GL quad override.

### 14.3 The value-gate companion decision (justified: SPY-ALONE + keep the @slow value gate)

**A fast (non-slow) sphere |Δk| value-gate does NOT make sense — the spy alone
is the right fast net.** Reasoning (Mode-10/11 logic):

- The seed-drop's VALUE signal (|Δk|=3.46e-2) requires the sphere's lagged
  Carlson closure ACTIVE *and* a physical config — which is GL-S16/40-cell =
  SLOW (§17 #6). A small/cheap sphere either NaNs (LS-S4/16-cell) or nulls the
  sensitivity → a fast sphere value-gate that ACTIVATES the drop does not
  exist.
- Cylinder value-gates are VACUOUS for seed-drops (telescoping, §14a); a
  cylinder companion would be a false-green.
- A fast value-gate below the activation floor is a Mode-10 sub-floor gate —
  **worse than none** (it manufactures false confidence). AGENT.md §5: "value
  gates below the floor are worse than none."

**Therefore:** the SPY (§14.1) is the fast net (Mode-11, pins the PATH in
seconds); the strong VALUE catcher STAYS `@pytest.mark.slow` (the existing
het-2G sphere eigenvalue) — mark it explicitly `catches` the seed-drop and
keep it in the full (`not "not slow"`) gate. Do NOT add a fast sphere
value-gate.

### 14.4 Mutations

| ID | Target | Mutation | Expected RED |
|---|---|---|---|
| **M-SEED-DROP** | §14.1 spy | patch `_solve_with_seed` → `self.L.solve(rhs)` (no `initial_guess`) | spy `seed_k is not None` reds |
| **M-SEED-ZERO** | §14.1 spy | thread `BoundaryFlux.zeros`/`rhs` instead of `psi_prev` | spy values-equal reds |
| **M-SEED-STALE** | §14.1 spy | thread `calls[0]` (the cold seed) every step (not the fresh iterate) | spy values-equal reds at k≥2 |
| **M-SEED-VALUE** (the @slow companion) | the kept @slow het-2G sphere eig | same drop, run the SLOW gate | |Δk|≈3.46e-2 reds — proves the slow value gate is the true value catcher |

**M-SEED-DROP is the §17 F4 witness**: it must red the FAST spy (seconds),
where today it reddens only the `@slow`-deselected value gate. That delta IS
the reason this prerequisite exists.

### 14.5 Homes + markers (Step 3 prereq)

| Home | Contents | Marker |
|---|---|---|
| NEW `tests/sn/solve/test_seed_threading_spy.py` | the §14.1 spy | `foundation` + `sentinel`, **not slow** |
| EXISTING het-2G sphere eigenvalue gate (find via `grep -rn "slow" tests/sn/eigenvalue/`) | the kept value catcher; add `@catches` for the seed-drop | `slow` |

The spy pins the CONTRACT (Step 3 keeps it green through the rewire); the
§3/§3b consumer-migration pinning gates (already specced) pin the end-to-end.
**M-PROBE (§7) is the sibling**: it drops `initial_guess` in the migrated
`SourceIteration` — the spy is its fast, always-on complement.

---

## 15. SUMMARY — each gate, its home, its teeth, and the land-order

**Implementation order (test-first; the harness is editing production LIVE, so
per operator: author the gate → confirm RED under its mutation → confirm GREEN
on the landed method).**

| # | Gate | Home | Mutation (teeth) | Lands relative to production |
|---|---|---|---|---|
| **STEP 1** | | | | |
| 1.1 | I1 round-trip + closed-form (numerics family) | NEW `test_inverse_universal.py` | M-INV-APPLY (non-unit c) | WITH each `.inverse()` |
| 1.2 | I2 functoriality (α, product, involution) | `test_inverse_universal.py` | M-INV-ORDER (non-commuting), M-INV-SCALE (α≠1) | WITH `.inverse()` |
| 1.3 | I3 inverse-object faithfulness | EXTEND `test_operator_capability_predicates.py` | M-INV-FAITH | WITH `.inverse()` |
| 1.4 | Negative controls (singular/composite) | EXTEND `test_diagonal_operator.py` / `test_multiplication_operator.py` | M-INV-NOGUARD (true zero) | WITH the eager guard |
| 1.5 | Permutation I1 (non-involution) | `test_inverse_universal.py` / `test_permutation_operator.py` | M-INV-PERM (3-cycle) | WITH `PermutationOperator.inverse` |
| 1.6 | MultiplicationOperator I1 (FullField) | EXTEND `test_multiplication_operator.py` | M-INV-APPLY | WITH `MultiplicationOperator.inverse` |
| 1.7 | Shim forwarding ± | `tests/sn/…/test_boundary*.py` | M-INV-SHIM | WITH the shim `.inverse()` |
| **STEP 2** | | | | |
| 2.1 | W2-round-trip (machine-precision on M) | NEW `test_gauss_seidel_reification.py` | M-SPLIT | WITH the `M` reification |
| 2.2 | W2-split exactness | ↑ | M-SPLIT | WITH the `B_lower` mask |
| 2.3 | W2-FP Mode-9-safe (diagonal cubature) | ↑ | (rides M-SPLIT) | AFTER 2.1 |
| 2.4 | W1 deforestation (product ≡ post-projection) | RENAME in `test_2d_anisotropic_windowing.py` | M-DEFOREST (drop ℓ≥1) | WITH the `P @ A.inverse()` product |
| 2.5 | coisometry `M∘R=4π·I` (cite, not re-mint) | EXISTING `test_frame.py::test_pi_R_is_4pi_identity…` | (landed) | — |
| **STEP 3 PREREQ** | | | | |
| 3.0 | seed-threading spy (fast, Mode-11) | NEW `test_seed_threading_spy.py` | M-SEED-DROP/ZERO/STALE | **BEFORE** the Step-3 driver rewire |
| 3.1 | @slow het-2G sphere value catcher | EXISTING eig gate + `@catches` | M-SEED-VALUE (|Δk|≈3.46e-2) | kept @slow |

**The keystone threading all steps** (§9/§12): `inv.apply(A.apply(x)) ≈ x` for
every `(A, A.inverse())` pair — REDs the confused resolvent (2.667, §13.1),
passes the clean inverse (3.8e-16). **Every gate above names the mutation that
reddens it AND the config that activates that mutation** (§0.5 gate-liveness +
§0.6 config-blindness — the standing discipline).

---

# PART III — STEP 4: `GreenOperator` + `OperatorSum.inverse()` + the wrap-delegate mixin (2026-07-02)

> **Why this Part exists.** Steps 1–3 (PART II §12–§14) LANDED the
> exact/direct inverse family (`InverseOperator`, `SweepOperator`, the
> nine advertisers, the reified `M`, the seed-threading spy). Step 4
> (taxonomy §12 step 4, §17 W3) adds the FIRST *iterative* inverse —
> `GreenOperator` = `A_loss.inverse()` wrapping `SourceIteration` — and
> the `OperatorSum` invertibility contract change that routes to it, plus
> the wrap-delegate MIXIN extraction the two landed twins' docstrings
> already name as their collapse TRIGGER (operator.py `InverseOperator`
> docstring; sweep_operator.py module docstring). This Part is the
> verification net for that step, main-agent-direct (surgical). It EXTENDS
> PART II; the §0.5 gate-liveness / §0.6 config-blindness standing
> disciplines and the §7/§12.6 mutation bank remain in force.

**The four production cuts this Part gates** (brief, verified against the
live `operator.py`/`iteration.py`/`sweep_operator.py` @ HEAD — line anchors
DRIFT, target behaviour):

1. **Mixin extraction** — the wrap-delegate back-half (`capabilities =
   {CAP_APPLY, CAP_SOLVE}` / domain↔codomain swap / `solve→inner.apply` /
   `is_invertible=True` / `inverse()→inner`) shared byte-identically by
   `InverseOperator` + `SweepOperator` extracts to a shared mechanism
   mixin at the 3rd sibling. The mixin declares the ABSTRACT canonical
   seeded-apply signature `apply(x, /, *, initial_guess=None)` — pyright
   LSP-enforces the kwarg on every sibling override (#285's structural
   resolution). Siblings keep: ctor guard, `apply` body, `__repr__`.
   **Behavior-preserving** (existing suites stay green unchanged).
2. **`GreenOperator(inner: OperatorSum, *, max_iter=1000, tol=1e-8)`** —
   the 3rd sibling, `orpheus/numerics/green_operator.py`. Split derivation:
   flatten the sum's LEFT SPINE into terms (exact-type `OperatorSum` nodes
   ONLY — subclasses like `InvertibleOperator` are leaf-terms); the LEADING
   term must be `is_invertible` (else `MissingCapability` naming the #261
   canonical ordering) → preconditioner `_seeded_inverse(leading)`; the
   rest → NEGATED gains (`ScaledOperator(-1,t)`, with the `ScaledOperator(
   -1,X)→X` unwrap so `A−S` yields gains `[S]`). Driver:
   `SourceIteration(precond, *gains, max_iter, tol)` (Green WRAPS the
   driver — taxonomy §11.2, never re-implements). `apply(q, *,
   initial_guess=None)` threads the seed as the driver's START and
   **raises `ConvergenceFailure`** (new, lives in green_operator.py) on
   non-convergence. Back-half via the mixin.
3. **`OperatorSum` contract change** — `is_invertible → self.a.is_invertible`
   (the recursion makes the left-spine leading term the designated
   preconditioner); conditional `CAP_SOLVE` in `__init__` (lockstep with
   `is_invertible`); `solve(b) = self.inverse().apply(b)`; `inverse()` →
   LATE import of green_operator, `return GreenOperator(self)` (W3
   placement; mirrors the streaming.py `InvertibleOperator.inverse()` →
   `SweepOperator` precedent).
4. **`_seeded_inverse` docstring scope-widening** — the family now
   structurally carries the kwarg via the mixin; residue (composed
   `OperatorProduct.inverse()` still lacks it) = #285.

**Claim layer (§1.5 gate).** Unlike steps 1–3 (behaviour-PRESERVING,
bit-id inheritance), **`GreenOperator` is the FIRST inverse in the family
with NO legacy `.solve` to inherit from** — the sum was NOT invertible
before, so there is no bit-identity twin to ride (lessons L13's
"equivalence keystone = bit-id INHERITANCE" does NOT apply here). Green's
correctness rests entirely on **structural-independence** anchors
(closed-form dense-LU + the G-Neumann multiple-scattering expansion), NOT
on inheritance. Every Green value gate is therefore a genuine
flux-shape/foundation claim against a structurally-independent reference —
**no eigenvalue claim, no MMS reference** (anti-#2). The `OperatorSum`
contract change (cut 3) IS behaviour-preserving for the value of every
EXISTING solve path (it only turns previously-non-invertible sums INTO
Green-invertible ones); its gates are foundation software-invariant
(faithfulness) + the migration of two frozen `is_invertible is False`
pins that flip by design.

---

## ✅ IMPLEMENTED @ `9333305` (2026-07-02) — as-built record + THREE deltas from this net

All four cuts + the §19–§24 gate set landed (main-agent-direct); **all 14
§25 mutations verified** — 12 bite their named gates, 2 designed-green
controls hold (the always-wrap `_negated` no-op control — proving the
unwrap is pure deforestation — and the M-GRN-SEED blindness proof: the
landed §14 spy + the L0 anchor stay GREEN under the seed drop while the
§23 spy reds). M-MIXIN-KWARG verified via CLI pyright on a cp-backup
file mutation (`reportIncompatibleMethodOverride` fires; file restored by
`cp`, never git). Pyright ratchet EXACTLY at baseline 148; Sphinx `-W`
EXIT 0. Homes as §26.1 (`test_green_operator.py` 11 gates,
`test_green_operator_sn.py` 5, migrations + extensions as mapped).

**Delta 1 — §18.A implemented as a REFINEMENT LOOP, not a check.** At
increment-stop the true residual is ρ/(1−ρ)·tol, so the recommended
check-only raise would falsely fire for EVERY split with ρ > 1/2 —
including §21.2's own het-slab gate. `Green.apply` therefore DRIVES the
promise: after each increment-stopped driver call it measures the true
residual and, if unmet with budget remaining, re-seeds the driver with
its own iterate (`steps` accumulates `len(history)` across calls against
ONE total `max_iter`). The driver remains the sole iteration engine —
the loop is tolerance bookkeeping only. The §18.A falsifier still
passes: the near-critical gate's budget is calibrated BETWEEN the
increment-stop (~460 steps at ρ=0.99/tol=1e-4) and the refinement close
(~920), so honest-Green raises where the increment-mutant returns
silently; the generous-budget leg is the refinement's positive control
(promise MET at true-residual grain). Two adversarial FP shapes found by
the divergence tooth are additionally caught in `apply` (documented
in-code): increment=nan, and the stopping-DENOMINATOR overflow false
convergence (`res = finite/inf = 0.0` onto a ~1e154 iterate).

**Delta 2 — §21.2's dense-LU anchor replaced by a trace-consistent
MANUFACTURED anchor.** The sketch's apply-to-basis LU would disagree
with the sweep-realized Green on arbitrary trace rows (#284: sweep
inverses are exact on the SOURCE subspace only). The landed gate
manufactures at the operator level: `x_tc = (L+C)⁻¹(random)` is
trace-consistent, `q = A_loss·x_tc` lies in the source subspace, and the
exact solution IS `x_tc` — same independence class (the anchor shares
only the forward operands, never the iteration), exact, no trace caveat.
The Green-vs-hand-SI equivalence is pinned at the DRIVER grain
(`green._driver.solve ≡ SourceIteration(sweep, S).solve`, bit-identical)
because `green.apply` legitimately runs refinement steps past the hand
loop.

**Delta 3 — M-MIXIN-CODOM's domain/codomain leg is
exercised-but-unconstrained (honest scope).** Every wrapped forward
today is endomorphic (`FullField→FullField`) or space-less
(`domain=None` leaves), so a swapped-the-wrong-way mixin
`domain`/`codomain` is value-invisible — no committed gate can red it
until a non-endomorphic invertible exists. The `solve→inner.solve`
half of that mutation IS constrained (the new back-half-solve anchor
reds it — M-MIXIN-SOLVE bites). Declared here per the Mode-10 third
state; revisit when a `Codomain ≠ Domain` invertible lands.

Also landed in passing: `operator.py`'s module docstring still carried
the pre-A-convention `(L−S−F)ψ=q` — fixed to `(A−S−F)ψ=q` (the code
sibling of the archivist's same-day theory-page L→A sweep).

---

## 18. Two DESIGN OBJECTIONS to the planned cuts (surface BEFORE writing gates)

The brief invites design objections. Two, each with its falsifier. **A is
a firm correctness recommendation; B is acceptable-with-documentation.**

### 18.A (FIRM) — `ConvergenceFailure` MUST test the TRUE residual, not the driver's increment

`SourceIteration.solve` returns `(psi, residual_history)` and does NOT
raise on `max_iter` exhaustion (iteration.py — the loop `break`s only when
`res < tol`, else returns the final iterate); its `res` is the **iterate
INCREMENT** `‖ψₙ−ψₙ₋₁‖ / ‖ψₙ‖`. The brief says Green raises "when the final
relative residual ≥ tol". **If "residual" is read as the driver's
increment `residual_history[-1]`, the exception does NOT fulfil its stated
purpose.** By numerical-bug-signatures **Signature 9** (ρ-blind stopping)
the increment UNDERSTATES the true error by `1/(1−ρ)`: as the scattering
ratio `c → 1` (near-critical / thick-reflective) the driver reports
`increment < tol` while `ψ` is still far from `A⁻¹q`. So an
increment-check would let a NON-converged iterate pass silently — exactly
the wrong answer the exception exists to prevent (Cardinal Rule 1).

**Recommendation:** `GreenOperator.apply` computes and checks the **true
relative residual** `‖inner.apply(ψ) − q‖ / ‖q‖ ≥ tol` (one extra matvec,
cheap relative to the iterative solve), NOT `residual_history[-1]`.

**Falsifier / gate (M-GRN-INCREMENT, §25):** a synthetic split `A = D − N`
with `ρ(D⁻¹N) = 0.99` (near-critical) whose exact answer is known
(`np.linalg.solve`). At a modest `max_iter`, the increment drops below
`tol` while the true residual is `~100×` larger. A Green that checks the
increment returns the wrong iterate SILENTLY (gate reds — it should have
raised); a Green that checks the true residual raises `ConvergenceFailure`.
This same gate makes the exception's teeth testable.

### 18.B (ACCEPTABLE w/ doc) — `is_invertible = self.a.is_invertible` is "leading-term-preconditionable", not "mathematically invertible"

The recursion conflates two different predicates. For `(−S) + A` (leading
term spelled non-invertible), the OPERATOR `A − S` IS mathematically
invertible (it is the within-group loss — Green-invertible when spelled
`A − S`), yet `((−S)+A).is_invertible = (−S).is_invertible = False`. So the
predicate reports `False` for a mathematically-invertible operator — a
narrowing of the north-star meaning of `is_invertible`.

**Why acceptable:** (a) the **#261 canonical-ordering rule** already makes
operand ORDER semantically load-bearing (`L+C` fuses, `C+L` does not), so
`is_invertible` is honestly read as "invertible via the leading-term
preconditioner AT THIS SPELLING"; (b) the §17 consumer audit found **ZERO**
production consumers of `is_invertible` outside composer propagation, so no
consumer relies on the spelling-independent meaning TODAY; (c) the refusal
is LOUD — `(−S)+A).inverse()` raises `MissingCapability` naming the
canonical ordering (`A − S`), not a silent wrong answer.

**Recommendation:** the `OperatorSum.is_invertible` docstring MUST state it
is "leading-term-preconditionable at this operand order (§261)", NOT
"mathematically invertible"; and the `(−S)+A` refusal message is
contract-PINNED (§22, §25 M-GRN-ORDER) so the semantic is tested, not
merely documented.

**Falsifier (latent):** a FUTURE consumer that branches on `is_invertible`
meaning spelling-independent mathematical invertibility would mis-handle
`(−S)+A`. None exists today; if one appears, this predicate is the site to
revisit. Filed as a note on the #261 canonical-ordering ruling.

---

## 19. STEP 4a — the wrap-delegate MIXIN extraction (behaviour-preserving)

**Claim layer:** foundation (a pure refactor — the back-half moves, no
value/behaviour changes). **Reference:** the existing landed suites
(inheritance — they must stay bit-identical) + the pyright static contract.

The mixin carries: `capabilities = {CAP_APPLY, CAP_SOLVE}`, `domain`
(=`inner.codomain`), `codomain` (=`inner.domain`), `solve(b) =
inner.apply(b)`, `is_invertible = True`, `inverse() = self.inner`, and the
**abstract** `apply(x, /, *, initial_guess=None)`. Siblings
(`InverseOperator`, `SweepOperator`, and now `GreenOperator`) keep only:
ctor guard (`self.inner = …`), the concrete `apply` body, `__repr__`.

### 19.1 The no-op proof — existing suites stay GREEN UNCHANGED

The extraction is behaviour-preserving iff every landed gate on the two
existing siblings passes with NO edit and NO re-baseline:

| Sibling | Suite (must stay green, unchanged) |
|---|---|
| `InverseOperator` | `test_inverse_universal.py` (I1/I2/faithfulness, §12), `test_diagonal_operator.py` + `test_multiplication_operator.py` singular ± pairs |
| `SweepOperator` | `test_inverse_operator_equivalence.py` (§2.1 equivalence + seed), `test_removal_form_matvec_sweep_roundtrip` (§6 canary, sphere excluded), the seed-threading spy (§14) |

No new runtime gate is needed for the extraction ITSELF (the above ARE its
regression net); the NEW gate is the STATIC contract (§19.2) + the mixin
mutation bank (§19.3).

### 19.2 The static gate — the abstract seeded-apply signature (the #285 resolution)

EXTEND `tests/sn/operators/test_operators_apply_typed.py` (§6a — one home
per apply-contract). The mixin's abstract `apply(x, /, *,
initial_guess=None)` is what LSP-enforces the kwarg on every sibling; make
that a pyright-visible fact:

```python
# every sibling accepts the canonical seeded-apply signature
assert_type(InverseOperator(D).apply(x), <carrier>)
assert_type(SweepOperator(lc).apply(rhs, initial_guess=seed), TimedFullField)
assert_type(GreenOperator(a_loss).apply(q, initial_guess=seed), <carrier>)
si: SupportsSeededApply = GreenOperator(a_loss)   # structural conformance (no reportAssignmentType)
```

**Prove teeth (CLI-pyright, lessons L10):** mutate ONE sibling's `apply`
override to drop the kwarg (`def apply(self, x, /)`) → `npx pyright
orpheus/` reports `reportIncompatibleMethodOverride` against the mixin's
abstract signature (M-MIXIN-KWARG, §25). This is the structural proof that
"every sibling carries the kwarg" is ENFORCED, not hoped — the #285
structural-vs-convention decision resolves to STRUCTURAL (the mixin
declares it; pyright checks it).

### 19.3 Mixin mutations — the shared back-half reds ALL THREE siblings

The whole POINT of the mixin (Pattern 2: single source of truth) is that a
back-half bug becomes UN-spellable-per-sibling. Prove it: a mixin mutation
must red gates on every sibling, not one.

- **M-MIXIN-CODOM** — mixin swaps `domain`/`codomain` the wrong way (or
  `solve → inner.solve` instead of `inner.apply`) → the involution /
  round-trip / equivalence gates red on InverseOperator AND SweepOperator
  AND GreenOperator.
- **M-MIXIN-INV** — mixin `inverse()` returns `self` not `self.inner` →
  `(A⁻¹)⁻¹ is A` (I2 involution) reds on all three.

(These live in §25's consolidated table; run them after 4b lands so all
three siblings exist.)

---

## 20. STEP 4b — `GreenOperator` joins the universal contract (I1/I2 + the tolerance extension)

**Claim layer:** foundation + flux-shape (the value against a
structurally-independent dense LU). **Reference:** closed-form
(`np.linalg.solve` of the materialized sum) + the G-Neumann expansion
(§21). **New file:** `tests/numerics/test_green_operator.py`,
`@pytest.mark.foundation`.

### 20.0 The §12.0 per-type tolerance table — the Green row (task A.1 extension)

`GreenOperator` is ITERATIVE, so its round-trip tolerance is DRIVER-TOL,
NOT nulp (lessons L7 — iterative → `SAFETY × conv_tol`). This is the one
row that breaks the §12.0 "count the arithmetic ops → nulp" rule, because
the FP path is an ITERATION, not a fixed reduction depth:

| Advertiser | `.inverse()` returns | `apply` realization | `inv.apply == solve`? | **I1 round-trip tol** | `(A⁻¹)⁻¹` form |
|---|---|---|---|---|---|
| `GreenOperator` (`numerics/green_operator.py`) | `OperatorSum.inverse()` → `GreenOperator(self)` for a Green-invertible sum | `SourceIteration(precond, *gains).solve(q, initial_guess)` | **N/A — tautological** (no legacy solve; `OperatorSum.solve := inverse().apply`) | `rtol = SAFETY(10) × tol` (DRIVER-TOL; L7 `ρ/(1−ρ)` headroom) | **`is inner`** (object id, via mixin) |

**⚠ The `inverse().apply ≡ solve` gate is a TAUTOLOGY for Green /
OperatorSum — do NOT include it as evidence.** In step 1 that gate was
bit-id INHERITANCE (the inverse WRAPPED a pre-existing `.solve`). Here
`OperatorSum.solve` is DEFINED as `self.inverse().apply(b)` (cut 3) — so
`sum.solve(b) == sum.inverse().apply(b)` proves NOTHING (it is the
definition). The evidence that Green is correct is the **dense-LU anchor**
(§21), never a solve-equivalence. State this explicitly so a later reader
does not add a green tautology and mistake it for coverage.

### 20.1 Gate G-I1 — the universal round-trip at DRIVER TOL + closed-form anchor

`Green.apply(A.apply(x)) ≈ x` AND `A.apply(Green.apply(q)) ≈ q`, at `rtol =
SAFETY × tol`, PAIRED with the closed-form anchor (L2 — round-trip alone is
self-consistency; both arms share the operands):

```python
# L0 dense split A = D - N, ρ(D⁻¹N) < 1 (convergent), NON-symmetric, N≠0
A = D - N                       # OperatorSum(a=D, b=ScaledOperator(-1,N))
green = A.inverse()             # GreenOperator
q = _nontrivial(seed=226)
np.testing.assert_allclose(green.apply(A.apply(x)), x, rtol=SAFETY*green.tol)
np.testing.assert_allclose(A.apply(green.apply(q)), q, rtol=SAFETY*green.tol)
np.testing.assert_allclose(                                   # STRUCTURAL anchor
    green.apply(q), np.linalg.solve(A_dense, q), rtol=SAFETY*green.tol,
    err_msg="Green ≠ dense LU of the materialized sum")
```

**Config-blindness (§0.6):** the split MUST have `N ≠ 0` and be
CONVERGENT (`ρ(D⁻¹N) < 1`) but NON-TRIVIAL (a zero `N` makes Green ≡ `D⁻¹`
and nulls every gain mutation, §25). Use a `D` diagonal, `N` off-diagonal,
`ρ ≈ 0.5`.

### 20.2 Gate G-I2 — functoriality (involution by object identity)

`(G)⁻¹ → the sum, by OBJECT IDENTITY` (mixin `inverse()→inner`):

```python
A = D - N
assert A.inverse().inverse() is A          # OBJECT identity (not action-nulp)
```

`(αA)⁻¹`/`(AB)⁻¹` functoriality is NOT re-tested on Green here — those are
the `ScaledOperator`/`OperatorProduct` composites of step 1 (already
gated); a Green nested in a product/scale is a step-5 composition concern
(defer). The Green-specific I2 is the involution-by-identity only.

### 20.3 Gate G-I3 — the Green OBJECT is faithful (extend §12.3)

EXTEND `tests/numerics/test_operator_capability_predicates.py`
(`test_inverse_objects_are_faithful`): `GreenOperator` advertises
`{CAP_APPLY, CAP_SOLVE}` (via mixin), `is_invertible=True`,
`is_adjointable=False` (adjoint-inverse is #280) → `_assert_faithful`
passes on both axes.

```python
green = (D - N).inverse()
_assert_faithful(green)                     # is_inv≡CAP_SOLVE, is_adj≡CAP_APPLY_TRANSPOSE
assert green.inverse() is (D - N-equivalent inner)   # involution, object id
```

---

## 21. STEP 4c — the NAME-EARNING G-invariants (the taxonomy §13 bar)

The user's bar (taxonomy §13): "to be NAMED `GreenOperator` it must be
TESTED to have Green-operator invariants" — a distinguishing property a
bare `InverseOperator` does NOT automatically have. Rulings on which of the
four candidate invariants are LOAD-BEARING vs redundant-with-anchor:

| Invariant | Ruling | Rationale |
|---|---|---|
| **G-Neumann** | **LOAD-BEARING** (L0 + L1) | The distinguishing invariant — a generic `A⁻¹` has NO splitting to satisfy. No other gate proves the multiple-scattering expansion structure. |
| **forward dense-LU anchor** | **LOAD-BEARING** (L0 + L1) | The correctness "was right" (structurally independent of the iteration). G-Neumann's convergence TARGET. |
| **G-reciprocity** | **INCLUDE (name-earning, lighter)** | The physical Green's reciprocity theorem; earns the name (§13). NOT eliminated: it is the CHEAP proof (no 2nd dense oracle) that the split derivation is correct for a DIFFERENT operand config (the transpose) — catches a leading-term-selection asymmetry. |
| **G-kernel** (`G.apply(δ_j)` = col j of A⁻¹) | **FOLD into the anchor** | Redundant as a standalone — it IS the forward anchor evaluated on unit vectors. Make the §20.1/§21.1 anchor input set INCLUDE the `δ_j` basis so the "column of the inverse = unit-point-source flux" reading is exercised without a separate gate. |

### 21.1 Gate G-Neumann-L0 — the splitting structure (§17 falsifier-5, run)

The partial sums of the preconditioned Neumann series converge geometrically
to `Green.apply(q)` AND to the independent dense LU:

```python
# A = P - N (P = leading invertible term, N = gain); ρ(P⁻¹N) = c < 1
P_inv, N = D.inverse(), N_op
partial = P_inv.apply(q)                      # k=0 term
acc = partial.copy()
ratios = []
for k in range(1, K):
    partial = P_inv.apply(N.apply(partial))   # (P⁻¹N)^k P⁻¹ q
    ratios.append(norm(partial)/norm(prev))   # geometric: → c
    acc = acc + partial
np.testing.assert_allclose(acc, green.apply(q), rtol=SAFETY*green.tol)       # Neumann → Green
np.testing.assert_allclose(acc, np.linalg.solve(A_dense, q), rtol=1e-10)     # → dense LU (independent)
assert abs(np.mean(ratios[-3:]) - c) < 0.05   # decay ratio ≈ the scattering ratio (0.4995≈0.5 in §17)
```

This is the §17 falsifier-5 shape, RUN (ratio 0.4995 ≈ c=0.5, agree
4.2e-13 by k=40). It is the LOAD-BEARING name-earner. Include the `δ_j`
basis in the `q` set (G-kernel fold): `green.apply(δ_j) == inv(A_dense)[:,j]`.

### 21.2 Gate G-Neumann-L1 — SN `A − S`, het ≥2G VACUUM slab (Mode-9 discipline)

The L1 instance on the real SN operators. **The config MUST be a
heterogeneous ≥2G VACUUM slab, NOT an isotropic reflective box** (§0.6 /
Mode-9): a reflective box with flat flux NULLS the scattering-redistribution
the Neumann series expands (the multiple-scattering terms are trivial); the
vacuum slab's flux gradient ACTIVATES the streaming↔scattering coupling so
the series is non-trivial. ≥2G is mandatory for the scattering-matrix
convention (Mode-6 — a 1G `Σ_s` is symmetric and blind to `SigSᵀ` drift).

```python
# tests/sn/operators/test_green_operator_sn.py — @pytest.mark.foundation (some @slow)
A_loss = (L + C) - S            # het 2G vacuum slab; leading (L+C) = InvertibleOperator (leaf-term)
green  = A_loss.inverse()       # GreenOperator; precond = sweep((L+C)), gains = [S]
# G-Neumann: Σ_k ((L+C)⁻¹ S)^k (L+C)⁻¹ q  →  green.apply(q), ratio ≈ c
# anchor: dense LU of the MATERIALIZED (L+C-S) (apply-to-basis) — structurally independent
A_dense = _as_dense(A_loss)     # apply-to-basis; small slab
np.testing.assert_allclose(green.apply(q), np.linalg.solve(A_dense, q_flat).reshape(...), rtol=SAFETY*green.tol)
```

**Green-vs-hand-SI equivalence (the auto-split-reproduces-canonical-SI
pin).** The flatten's whole job is to reproduce the CANONICAL within-group
source iteration `sweep(q + Sψ + Bψ)` as ONE driver loop. Pin it
(bit-id-inheritance analog — necessary-NOT-sufficient, so it RIDES the
dense-LU anchor for "was right"):

```python
green_result, _ = green.apply(q), ...
hand = SourceIteration(_seeded_inverse(L+C), S, max_iter=..., tol=green.tol).solve(q)[0]
np.testing.assert_array_equal(green.apply(q).bulk.values, hand.bulk.values)   # auto-split ≡ hand-split
```

### 21.3 Gate G-reciprocity-L0 — the Euclidean transpose (no `.H` minting)

Build `Gᵀ` MANUALLY as the Green over the transposed operands (NOT via
`Green.H` — that is #280-adjacent). `A = P − N` NON-symmetric; `Aᵀ = Pᵀ −
Nᵀ` built from the transposed operands; check the reciprocity bilinear
identity in the EUCLIDEAN inner product:

```python
A   = P - N                     # non-symmetric dense split
A_T = P.T_op - N.T_op           # transposed-operand sum (manual, not A.H)
gA, gAT = A.inverse(), A_T.inverse()
for phi1, phi2 in random_pairs:
    lhs = np.dot(phi2, gA.apply(phi1))          # ⟨φ₂, G φ₁⟩
    rhs = np.dot(gAT.apply(phi2), phi1)          # ⟨Gᵀ φ₂, φ₁⟩
    np.testing.assert_allclose(lhs, rhs, rtol=SAFETY*gA.tol)
```

**Division of labour (why this is not redundant):** the forward anchor
(§21.1) pins `gA` against `LU(A)` — but NOTHING else pins `gAT` (the
transposed-operand split). G-reciprocity pins `gAT`'s correctness via the
identity, WITHOUT a second dense oracle. **Config-blindness:** `A` MUST be
NON-symmetric (`A ≠ Aᵀ`), else `gAT == gA` and the reciprocity is vacuous
(§0.6). This catches a leading-term-selection ASYMMETRY bug (Green flattens
`Pᵀ − Nᵀ` picking the wrong leading term).

---

## 22. STEP 4c — the ORDERING-RULING teeth + `ConvergenceFailure`

The #261 order-dependent-strategy ruling has THREE edges + one exception
class. All four need teeth.

### 22.1 The four edges

| Spelling | Structure | `.inverse()` | apply behaviour | Gate |
|---|---|---|---|---|
| `L + C` | `InvertibleOperator` (MRO shadows the sum) | `SweepOperator` (override) | direct WDD sweep | `type((L+C).inverse()) is SweepOperator` (NOT Green) — the MRO shadow holds |
| `A_loss = (L+C) − S − B` | plain sum, leading `(L+C)` invertible | `GreenOperator` | converges → does NOT raise | §21.2 + `not raises` |
| `C + L` | plain sum, leading `C` invertible | `GreenOperator` | collision-preconditioned Richardson DIVERGES → raises | §22.2 |
| `(−S) + A` | plain sum, leading `(−S)` NOT invertible | REFUSES at construction | `MissingCapability` naming #261 | §22.3 |

### 22.2 Gate — the divergent split raises LOUDLY (`C+L` + the L0 analog)

The convergent case does NOT raise; the divergent case DOES.

```python
# L0 analog: D + K with ρ(D⁻¹K) > 1 (streaming-dominant transport shape)
div = D + K_big            # OperatorSum(a=D invertible, b=K), ρ(D⁻¹K) > 1
with pytest.raises(ConvergenceFailure, match="converge|residual"):
    div.inverse().apply(q)
# convergent control — MUST NOT raise (a bare pytest.raises with no control is Mode-10 blind)
conv = D - N_small         # ρ < 1
conv.inverse().apply(q)    # no raise
# L1: C + L (transport) — collision-preconditioned Richardson diverges
with pytest.raises(ConvergenceFailure):
    (C + L).inverse().apply(q)   # het 2G; C = MultiplicationOperator leading term
```

`ConvergenceFailure` is a `raise` (NOT a bare assert — survives `-O`,
Mode-8). Its check is the TRUE residual (§18.A), not the increment.

### 22.3 Gate — leading-non-invertible refuses at CONSTRUCTION with the canonical-ordering message

```python
with pytest.raises(MissingCapability, match="canonical|ordering|leading|invert"):
    ((-1.0 * S) + A).inverse()     # leading (−S) not invertible → refuse; DON'T build a Green
```

**Config-blindness:** the leading term MUST be genuinely non-invertible
(`S` = scattering, apply-only). The message match pins the §18.B semantic
(the refusal names the canonical spelling `A − S`).

---

## 23. Seed-threading — RULING: Green needs its OWN Mode-11 pin (the landed spy is BLIND to it)

**The question (brief item 4):** is `Green.apply(q, initial_guess=x0)`
threading into the driver covered by the landed seed-threading spy (§14),
or does it need its own pin?

**RULING: Green needs its OWN Mode-11 pin.** The landed spy
(`test_seed_threading_spy.py`) wraps `InvertibleOperator.solve` and pins
the DRIVER's per-iterate threading `driver → precond.apply(rhs,
initial_guess=psi_prev)` — a level BELOW Green's driver-start. When
`Green.apply(q, initial_guess=x0)` runs it threads `x0` as the START of its
INTERNAL `SourceIteration.solve(q, initial_guess=x0)` — a DIFFERENT
threading. If `Green.apply` DROPS `x0` (calls `driver.solve(q)` bare), the
internal per-iterate threading STILL works (each iterate seeds the next),
so the landed spy stays GREEN. The drop is also VALUE-INVISIBLE — a warm
start changes only the inner convergence RATE, not the converged value
(Green converges to `A_loss⁻¹q` regardless of start), so a value gate
cannot see it either (Mode-11 territory, exactly like the §14 sphere case).

**The gate (in `test_green_operator.py`):** an in-process wrap on the
INTERNAL driver's `solve`, asserting the FIRST call receives `x0` (values-
equal, robust to a defensive copy — do NOT assert `is`, L4/§14.1):

```python
def test_green_threads_initial_guess_to_driver_start(monkeypatch):
    A = D - N
    green = A.inverse()
    seen = []
    real = SourceIteration.solve
    def spy(self, q_ext, initial_guess=None):
        seen.append(initial_guess); return real(self, q_ext, initial_guess=initial_guess)
    monkeypatch.setattr(SourceIteration, "solve", spy)
    green.apply(q, initial_guess=x0)
    assert seen and seen[0] is not None, "Green dropped initial_guess to the driver start"
    np.testing.assert_array_equal(seen[0], x0)   # the EXACT start threaded through
```

Mutation **M-GRN-SEED** (§25): `Green.apply` calls `driver.solve(q)` w/o
`initial_guess` → this spy reds; the landed §14 spy + every value gate stay
GREEN (proving the landed spy is blind to Green's drop). This is the
`GreenOperator` analog of the §17 F4 finding — the fast, always-on Mode-11
net that the composition of existing gates does NOT provide.

---

## 24. STEP 4c — the `OperatorSum` contract change + the MIGRATION MAP

The contract change (`is_invertible → a.is_invertible`; conditional
`CAP_SOLVE`; `solve`; `inverse()→Green`) FLIPS two frozen `is_invertible
is False` pins from `False` to `True` BY DESIGN — the sum with an invertible
LEADING term now HAS a general inverse (the Green). This is the headline
migration the brief names.

### 24.1 The blast-radius audit (retirement-audit discipline — grep, don't guess)

`grep -rn "is_invertible is False" tests/` (WITHOUT filtering — an
adjoint-line filter HIDES `predicates.py:94`). Classify EACH hit: a LEAF /
PRODUCT / SCALED op stays `False` (unchanged rules); a PLAIN `OperatorSum`
with an invertible LEADING term FLIPS. Verified inventory:

| Site | Op | Verdict |
|---|---|---|
| `test_operator.py` `not matrix_apply_only.is_invertible` | apply-only leaf | **STAYS** False |
| `test_frame.py:69` `face.is_invertible is False` | frame face leaf | **STAYS** |
| `test_operator_capability_predicates.py:112` `p_singular` | `Diag(_CZ) @ Id` — OperatorPRODUCT (singular factor) | **STAYS** (`a and b` unchanged) |
| `test_capability_survival.py:247` `C_singular` | singular `MultiplicationOperator` leaf | **STAYS** |
| **`test_operator_capability_predicates.py:94`** `s_both = ident + Diag` | plain sum, leading `ident` invertible | **FLIPS → True** |
| **`test_capability_survival.py:276`** `(L+C) − B` | plain sum, leading `(L+C)` invertible | **FLIPS → True** (headline) |

### 24.2 The two migrating pins (REWIRE — the `is_invertible is False` premise is now FALSE by design)

| Pin (file::name) | Today asserts | New shape |
|---|---|---|
| `test_capability_survival.py::test_invertible_operator_is_the_sole_invertible_sum` (:266–276) | `(L+C)` is the SOLE invertible sum; `(L+C)−B).is_invertible is False  # no general (A+B)⁻¹` | **The premise is retired.** `(L+C)−B).is_invertible is True`; `type((L+C-B).inverse()) is GreenOperator`; the comment `# no general (A+B)⁻¹` DELETES (a sum with an invertible leading term HAS a general inverse: the Green). Rename the test → `test_sum_with_invertible_leading_term_is_green_invertible`. KEEP the adjoint-axis assertions unchanged. |
| `test_operator_capability_predicates.py::test_sum_predicates_recursive_and_faithful` (:90–98) | `s_both.is_invertible is False`; docstring "never invertible (no general sum inverse)" | `s_both.is_invertible is True` (leading `ident` invertible); add the recursive pin `s_both.is_invertible == s_both.a.is_invertible`; docstring → "invertible iff the LEADING term is (§261 canonical order)". `s_half = ident + apply_only` also flips to `is_invertible True` (leading `ident` invertible) — its `.inverse()` is a Green that may or may not converge; keep the `is_adjointable is False` pin. |

### 24.3 The faithfulness keystone STAYS GREEN — but only if `CAP_SOLVE` is added in LOCKSTEP

`tests/_harness/predicates.py::assert_capability_faithful` (`is_invertible
≡ CAP_SOLVE∈caps`) is UNCHANGED (it is the definition). It stays green iff
cut 3 adds `CAP_SOLVE` to `OperatorSum.__init__` under the EXACT SAME
condition as `is_invertible`:

```python
# OperatorSum.__init__  (both read self.a.is_invertible — single source)
if self.a.is_invertible:
    caps.add(CAP_SOLVE)
# OperatorSum.is_invertible property:  return self.a.is_invertible
```

**The two asymmetric fixtures that PIN the `a`-recursion (config-blindness
§0.6 — a both-invertible sum is BLIND to a/b/and/or).** ADD to the
faithfulness enumeration:

```python
inv_then_ao = IdentityOperator() + _ApplyOnly()   # a inv, b not → is_inv True,  CAP_SOLVE∈caps
ao_then_inv = _ApplyOnly() + IdentityOperator()   # a not, b inv → is_inv False, CAP_SOLVE∉caps
assert inv_then_ao.is_invertible is True  and ao_then_inv.is_invertible is False
_assert_faithful(inv_then_ao); _assert_faithful(ao_then_inv)
```

The `(True, False)` pattern for `(inv+ao, ao+inv)` UNIQUELY identifies the
`a`-rule (b-rule→(F,T); and→(F,F); or→(T,T)). Both fixtures are mandatory.

### 24.4 The dispatch / MRO gates (cut 3 routing)

```python
assert type((L + C).inverse()) is SweepOperator      # InvertibleOperator override shadows the sum (MRO)
assert type((C + L).inverse()) is GreenOperator       # plain sum → Green
assert type(((L+C) - S).inverse()) is GreenOperator   # A_loss → Green
```

---

## 25. MUTATION RECIPES — Step 4 (each REDs a NAMED gate, under `-O`, monkeypatch only)

In-process monkeypatch ONLY (uncommitted edits — NEVER `git checkout`,
process-discipline.md). Every mutation names its ACTIVATING config (§0.6
applies to mutations).

| ID | Target gate | Mutation | Expected RED | Activating config |
|---|---|---|---|---|
| **M-MIXIN-KWARG** | 19.2 static | drop `initial_guess` from one sibling's `apply` override | `npx pyright` `reportIncompatibleMethodOverride` vs the mixin's abstract sig | any sibling |
| **M-MIXIN-CODOM** | 19.3 / involution | mixin swaps `domain`/`codomain` OR `solve→inner.solve` | round-trip/involution red on ALL THREE siblings | any |
| **M-MIXIN-INV** | 19.3 / G-I2 | mixin `inverse()` returns `self` not `self.inner` | `(A⁻¹)⁻¹ is A` red on all three | any |
| **M-GRN-SIGN** | G-Neumann + anchor (§21) | gain-sign flip: use `+t` as the gain (drop the `-1` negation) | driver solves `(P+N)x=q` → anchor + Neumann red | `N ≠ 0` (zero gain blind) |
| **M-GRN-UNWRAP** | anchor (§21) | negation-unwrap dropped: gain stays `ScaledOperator(-1,X)` → driver adds `−Xψ` (double sign) | anchor reds (wrong fixed point) | `N ≠ 0` |
| **M-GRN-SWAP** | anchor + Neumann | precond/gain swap: `_seeded_inverse(a gain)` as preconditioner | anchor reds / precond raises if the gain is non-invertible | leading term ≠ a gain |
| **M-GRN-FLATTEN** | flatten-exactness (§21.2/§24.4) | flatten THROUGH subclasses (treat `InvertibleOperator` as exact `OperatorSum`) | `A_loss.inverse()` REFUSES (leading `L` non-invertible) where it should build a Green | `A_loss=(L+C)−S−B` |
| **M-GRN-TOL** | ConvergenceFailure (§22.2) | delete the residual-vs-tol check in `Green.apply` | the divergent split (`C+L` / `D+K` ρ>1) returns garbage SILENTLY — the `pytest.raises` gate reds (no raise) | `ρ(P⁻¹N) > 1` |
| **M-GRN-INCREMENT** | ConvergenceFailure (§18.A) | check `residual_history[-1]` (increment) instead of the true `‖Ax−q‖/‖q‖` | a `c=0.99` split passes silently (increment<tol, true≫tol) — gate reds | `c ≥ 0.99` near-critical |
| **M-GRN-SEED** | seed spy (§23) | `Green.apply` calls `driver.solve(q)` w/o `initial_guess` | the Green Mode-11 spy reds; landed §14 spy + value gates stay green (proves the landed spy is blind) | any (1st driver call) |
| **M-GRN-ORDER** | refusal (§22.3) | `((−S)+A).inverse()` builds a Green anyway (ignore leading-non-invertible) | the canonical-ordering refusal gate reds (no raise) | leading term non-invertible |
| **M-SUM-INV-B** | is_invertible recursion (§24.3) | `OperatorSum.is_invertible = self.b.is_invertible` (or `a and b`) | `(inv+ao).is_invertible`/`(ao+inv)` flip wrong; faithfulness reds | BOTH asymmetric fixtures |
| **M-SUM-CAPDRIFT** | faithfulness lockstep (§24.3) | conditional `CAP_SOLVE` uses `_has(a) and _has(b)` (not `a.is_invertible`) | `assert_capability_faithful` reds (`is_invertible ≠ CAP_SOLVE∈caps`) | `inv+ao` (a inv, b not) |

**Mode-8 hygiene:** every gate uses `np.testing.*` / `pytest.raises` /
`assert_type` / collected-module `assert` (pytest-rewritten → fires under
`-O`). `ConvergenceFailure` and the leading-non-invertible refusal MUST
`raise` (not `assert`) so they survive `-O` in production (§7).

---

## 26. HOMES, MARKERS, and LAND-ORDER

### 26.1 Homes + markers

| Home | Contents | Marker |
|---|---|---|
| **NEW** `tests/numerics/test_green_operator.py` | G-I1 (driver-tol round-trip + dense-LU anchor + δ_j), G-I2 involution, G-Neumann-L0, G-reciprocity-L0, divergent `D+K` `ConvergenceFailure` (+convergent control + true-residual `c=0.99`), Green seed spy (§23) | `foundation`, not slow |
| **NEW** `tests/sn/operators/test_green_operator_sn.py` | G-Neumann-L1 (het 2G VACUUM slab), Green-vs-hand-SI equivalence, `(L+C)−S−B`/`(L+C)−B` invertibility flip, `(C+L)` divergent `ConvergenceFailure`, `(L+C).inverse()` still `SweepOperator` (MRO), `(−S)+A` refusal | `foundation`; the het-2G convergence cases `@slow` |
| EXTEND `tests/numerics/test_inverse_universal.py` | dispatch rows (`type((C+L).inverse()) is GreenOperator`) | `foundation` |
| EXTEND `tests/numerics/test_operator_capability_predicates.py` | G-I3 Green-object faithfulness; the two asymmetric sum fixtures (§24.3) | `foundation` |
| EXTEND `tests/sn/operators/test_operators_apply_typed.py` | the mixin abstract seeded-apply `assert_type` + `SupportsSeededApply` conformance (§19.2) | `foundation` (static) |
| MIGRATE `tests/sn/operators/test_capability_survival.py` (:266–276) | the headline flip → `test_sum_with_invertible_leading_term_is_green_invertible` (§24.2) | (keep existing) |
| MIGRATE `tests/numerics/test_operator_capability_predicates.py` (:90–98) | `s_both.is_invertible` flip + recursive pin (§24.2) | (keep existing) |

**NO** `verifies(<physics-eq>)` — all Step-4 gates are software invariants
(foundation; L9 — foundation MUST NOT carry a physics label). The L1 SN
cases are foundation (they verify the operator-algebra CONSTRUCTION, not a
theory equation), not `l1`+`verifies`.

### 26.2 Land-order (test-first — author gate → RED under its mutation → GREEN on landed method)

| # | Cut | Gates that land WITH it | Teeth |
|---|---|---|---|
| **4a** | mixin extraction (refactor the 2 existing siblings) | 19.1 no-op (existing suites green UNCHANGED), 19.2 static contract | M-MIXIN-KWARG (pyright) |
| **4b** | `GreenOperator` class (constructed DIRECTLY in tests — before routing) | G-I1, G-I2, G-Neumann-L0, G-reciprocity-L0, `ConvergenceFailure` (D+K + true-residual), Green seed spy; the mixin now has 3 users (M-MIXIN-CODOM/INV bite all three) | M-GRN-SIGN/UNWRAP/SWAP/TOL/INCREMENT/SEED |
| **4c** | `OperatorSum` contract (is_invertible=a.is_invertible + conditional CAP_SOLVE + solve + inverse()→Green) | §24.2 pin migration, §24.3 faithfulness re-verify + asymmetric fixtures, §24.4 dispatch/MRO, §22.3 `(−S)+A` refusal | M-SUM-INV-B, M-SUM-CAPDRIFT, M-GRN-ORDER, M-GRN-FLATTEN |
| **4d** | L1 SN wiring (`A_loss = (L+C)−S(−B)`) | G-Neumann-L1 (vacuum slab), Green-vs-hand-SI, `(C+L)` divergent end-to-end | M-GRN-FLATTEN (het), M-GRN-TOL (`C+L`) |

**Done when:** for EVERY Step-4 gate a named mutation reds it under `-O`;
`GreenOperator` earns its name (G-Neumann + G-reciprocity green,
mutation-proven); the two `is_invertible is False` pins migrated and the
faithfulness keystone re-verified green with the asymmetric fixtures; the
mixin extraction proven a no-op (existing suites unchanged) + its abstract
seeded-apply signature pyright-enforced; `ConvergenceFailure` fires on the
divergent split (true-residual, §18.A) and NOT on the convergent control.

---

# PART IV — STEP 5: `MatrixInverseOperator` + `as_matrix()` + `OperatorSum`/`OperatorProduct` closure (#285) (2026-07-02)

> **✅ IMPLEMENTED @ `d82fa77` (2026-07-02, `refactor/inverse-as-operator`;
> closes #285).** 18 gates in the NEW
> `tests/numerics/test_matrix_inverse_operator.py` + the §30.10/§31
> universal extensions + M-I3 faithfulness + the 4th-sibling static pins;
> **14/14 mutations verified** (13 in-process via the job-tmp
> `step5teeth/run_teeth.py` bank + M-MINV-KWARG via CLI pyright
> `reportIncompatibleMethodOverride`, cp-backup/restore); ratchet exactly
> 148; tier 2998/0; homogeneous k∞ pins byte-green. **As-built deltas
> (5, none weakening):** (1) the §32.2 HOMOG-EQUIV gate is HOSTED in the
> numerics file, not `test_homogeneous.py` — that file's file-level
> `l1 + verifies(...)` marks would have written FALSE equation-TESTS
> edges for a pure software gate (the registry warned on the conflict);
> the gate documents the home delta in its docstring. (2) §30.4
> MINV-DIRECT gained a NON-symmetric second leg (`D − N`, nilpotent) so
> M-MINV-LUTRANS's activating config matches its target gate — the
> spec's own §34 config column ("NON-symmetric inner; symmetric blind")
> contradicted its target column (§30.1/§30.4, whose diag fixture is
> symmetric). (3) the §29.2 per-call edges use tighten/lift on a
> 5-vector + a default-gate refusal at (4097,) — the same three edges
> with no 200 MB (5000,)-lift materialization. (4) the §30.6 non-square
> match pins the ctor's domain-language fragment ("SQUARE
> materialization"), not a loose "square" — scipy's own "expected square
> matrix" would otherwise satisfy the match and M-MINV-NOGUARD-SQ would
> be toothless. (5) the §30.3 Green-contrast sanity leg is `1e3·green.tol`
> (the spec's 10× under-budgeted the G·A column-error cond
> amplification); the DISTINGUISHING leg (resid > machine·cond) is
> unchanged. §31.4's negative already existed as
> `test_product_with_singular_factor_inverse_raises` — not duplicated.

> **Why this Part exists.** Steps 1–4 (PART II–III) landed the exact
> (`InverseOperator`/`SweepOperator`), reified (`M`) and iterative
> (`GreenOperator`) inverse families. Step 5 (taxonomy §12 step 5, §13
> M-row, §17 A2) adds the **materialising** family: `as_matrix()` as a
> `LinearOperator` BASE method (the promoted `_as_dense` apply-to-basis
> loop) and `MatrixInverseOperator` — the 4th `InverseWrapMixin` sibling
> that EAGERLY factorises `inner.as_matrix()` and inverts by `lu_solve`.
> It also CLOSES **#285**: `OperatorProduct.inverse()` stops returning a
> raw reversed product and returns `InverseOperator(self)` (the generic
> family member), so a composed inverse finally carries the canonical
> seeded-apply kwarg. Finally it RETIRES `_as_dense`, rewiring the
> homogeneous solver's two call sites to `op.as_matrix(basis_shape=(ng,1))`.
> This Part is the verification net for that step, main-agent-direct
> (surgical). It EXTENDS PARTS II–III; the §0.5 gate-liveness / §0.6
> config-blindness standing disciplines and the §7/§12.6/§25 mutation
> banks remain in force. Line anchors DRIFT (the file is being edited
> LIVE — `MatrixTooLarge` op.py:290 + `_resolve_basis_shape` op.py:310
> are already present as uncommitted scaffold on `9333305`); target
> BEHAVIOUR.

**The production cuts this Part gates** (brief; verified against the live
`operator.py`/`homogeneous/solver.py`/`isotropic_scattering.py`):

1. **`as_matrix()` on `LinearOperator`** (base method): signature
   `as_matrix(*, basis_shape: tuple[int,...] | None = None,
   max_dimension: int = 4096) -> np.ndarray`. The promoted `_as_dense`
   loop: per flat index `i` (C-order over `basis_shape`) build `e_i`,
   `op.apply(e_i)`, ravel (C-order), `column_stack` → shape
   `(prod(out_shape), prod(basis_shape))` (RECTANGULAR allowed — the
   output dim emerges from `apply`). Basis-shape resolution via the
   **single-source** `_resolve_basis_shape` (op.py:310, shared with the
   `MatrixInverseOperator` ctor — Pattern 2): explicit `basis_shape`
   wins; else `op.domain.shape`; else a `ValueError` naming BOTH
   remedies. Refuses `MatrixTooLarge(RuntimeError)` when
   `prod(basis_shape) > max_dimension` (§17 A2 — a RESOURCE effect on a
   TOTAL functor; NO `is_materializable` predicate). Return is a DENSE
   ndarray; the sparse-triangular streaming override (sweep_graph.py:66)
   is DEFERRED with its 3-D consumer.
2. **Mixin bound relaxation** — `_ForwardT`'s bound `_InvertibleForward →
   LinearOperator` (the back-half consumes only domain/codomain/apply).
   Pure widening; the 125+ landed family gates stay green UNCHANGED +
   pyright clean.
3. **`MatrixInverseOperator`** (NEW `orpheus/numerics/matrix_inverse_operator.py`,
   green_operator-style leaf): 4th `InverseWrapMixin` sibling. Ctor
   `(inner, *, basis_shape=None, max_dimension=4096)` EAGERLY materialises
   `inner.as_matrix(...)` (MatrixTooLarge propagates; non-square → raise;
   exact-singular → scipy `lu_factor` raises `LinAlgError` naturally).
   **Does NOT require `inner.is_invertible`** — the matrix reads VALUES
   not STRUCTURE (the §3 strategy-override witness: it inverts what Green
   REFUSES). `apply(q, /, *, initial_guess=None)` = accept-and-ignore seed
   (`del initial_guess` — M-direct IS seed-independence), `lu_solve`,
   reshape. `as_matrix()` OVERRIDE = one batched `lu_solve(lu, eye(n))`.
   Back-half ({apply, solve}, domain↔codomain swap, `solve=inner.apply`,
   `is_invertible=True`, `inverse()→inner` by object id) from the mixin.
4. **#285 closure** — `OperatorProduct.inverse()` → `InverseOperator(self)`
   (was `OperatorProduct(b.inverse(), a.inverse())`). `apply(q) =
   product.solve(q) = b.solve(a.solve(q))` — bit-identical composition;
   seed accept-and-ignore; involution STRENGTHENS to object identity;
   non-invertible product still raises `MissingCapability`.
   Permutation/Identity/Scaled/TensorProduct inverses stay ALGEBRA-CLOSED
   (a perm's inverse IS a perm — the OTHER kind of inverse; documented,
   not gated for seeded-apply).
5. **Homogeneous rewire** — `_as_dense` RETIRES; `_assemble_loss_matrix`
   (loss) + fission call sites → `op.as_matrix(basis_shape=(ng,1))`.
   `direct_eigenvalue` keeps its ndarray boundary. Output BIT-IDENTICAL
   (same basis columns, order, eig call). The full `K =
   MatrixInverseOperator(loss) @ F` spelling is #138 — OUT of scope.
6. **`dense_per_material` docstring honesty** (isotropic_scattering.py:149,
   :219) — doc-only: it has ZERO production consumers; it is the
   storage-side (`sig_s_legendre(mid)[0].T`) STRUCTURALLY-INDEPENDENT
   oracle the tests use, NOT a diffusion/homogeneous consumption path
   (homogeneous used `_as_dense`). No behaviour change.

**Claim layer (§1.5 gate).** `as_matrix()` value gates are **closed-form**
(hand-built matrices for `Diagonal`/`Permutation`; `dense_per_material` as
the structurally-independent storage oracle for the energy leaves).
`MatrixInverseOperator` is a **DIRECT dense inverse** — unlike Green
(iterative, driver-tol), its tolerances are **machine·cond** grain (NOT
`SAFETY×conv_tol`, NOT nulp): I1/M-materialise/M-direct assert against a
closed-form dense reference at LU precision. The `#285` product change and
the homogeneous rewire are **behaviour-preserving / bit-identical**
(foundation software-invariant + closed-form k∞ anchor). **No gate here
makes an eigenvalue claim on an MMS reference** (anti-#2); the homogeneous
k∞ rides its LANDED SymPy `case.k_inf` anchor (`test_kinf_exact`, 1e-12).

---

## 27. DESIGN PRESSURE-TEST FINDINGS (surface BEFORE writing gates)

The settled design is not to be redesigned — but three points are
under-specified or superseded by the LIVE tree. Flag, gate around, do NOT
silently paper over.

### 27.A (SHARPEN — the gate wording changes) — M-materialise's "a matrix-free inverse CANNOT satisfy it (no `as_matrix()`)" is SUPERSEDED

The taxonomy §13 M-row earns the name on: "**M-materialise** —
`Ainv.as_matrix() @ A.as_matrix() ≈ I` both ways — a matrix-free inverse
CANNOT satisfy it (no `as_matrix()`)." That parenthetical is now FALSE:
cut 1 makes `as_matrix()` a **universal base method** — a `GreenOperator`
inherits it (apply-to-basis), so `green.as_matrix() @ A.as_matrix() ≈ I`
ALSO holds (each Green column is `green.apply(e_j) ≈ A⁻¹e_j`). So
"no `as_matrix()`" no longer distinguishes `MatrixInverseOperator`.

**What genuinely earns the name is the PRECISION GRAIN, not the method's
existence.** `MatrixInverseOperator.as_matrix()` is `lu_solve(lu, eye)` —
the factorisation's own inverse at **machine·cond** precision; a Green's
inherited base `as_matrix()` is `n` iterative solves at **driver-tol**.
So the name-earning invariant is:

- **M-materialise gated at MACHINE precision** (`atol = K·eps·cond`, NOT
  `SAFETY×tol`): a Green passes it only to driver-tol (~1e-8), so the
  MACHINE-grain assertion is what a direct dense inverse satisfies and an
  iterative one does not. THIS is the distinguishing test.
- **M-direct's machine·cond residual** (the brief's own grain) is the
  companion — `‖A(A⁻¹q)−q‖ ≤ K·eps·cond·‖q‖`, un-achievable by an
  iterative inverse (whose residual floors at its tol).

**Recommendation:** gate M-materialise + M-direct at machine·cond, and add
one explicit CONTRAST assertion — the same `A`, a `GreenOperator(A)`, does
NOT meet the machine-grain residual (it meets only driver-tol) — so the
name-earner is proven DISTINGUISHING, not merely satisfied. (This is the
§21 "a bare invariant that a generic inverse also satisfies is not
name-earning" discipline, applied to the superseded framing.)

### 27.B (HONEST-SCOPE — the witness fixture changes) — the `(−S)+(L+C)` VALUE-vs-STRUCTURE witness is FullField, OUT of the ndarray `as_matrix` scope

The brief's strategy-override witness is `(−S)+(L+C)`: Green REFUSES it
(leading `−S` non-invertible) but `MatrixInverseOperator` inverts it.
BUT the real SN `(−S)+(L+C)` is a **FullField-carrier** operator, and
`as_matrix()`'s apply-to-basis builds **bare-ndarray** `e_i` — FullField
carriers are OUT of scope (honest-scope §35). Streaming `(L+C)` is also
"never materialised" (taxonomy §5). So `MatrixInverseOperator((−S)+(L+C))`
cannot materialise the FullField sum.

**Recommendation:** the witness uses an **ndarray-carrier structural
analog** that reproduces the EXACT refusal/inversion asymmetry: a synthetic
`neg_S_ao + D` where `S_ao` is an apply-only leaf (`_ApplyOnly`-style, a
fixed dense action, `is_invertible=False`) and `D` a leading… no — spell
it `(−S_ao) + D` so the LEFT-SPINE HEAD is `ScaledOperator(-1, S_ao)`
(`is_invertible = S_ao.is_invertible = False`). Then `((−S_ao)+D).inverse()`
→ Green REFUSES (canonical-ordering `MissingCapability`), while
`MatrixInverseOperator((−S_ao)+D)` materialises `D − S_ao` (both have
`apply`) and inverts it — anchored vs `np.linalg.solve(D_dense − S_ao_dense,
q)`. The STRUCTURE (leading-non-invertible) is what the witness proves; the
physics is incidental. State in-gate that the FullField `(−S)+(L+C)` is the
motivation, the ndarray split the realisable proof.

### 27.C (CONFIRM — two distinct error classes) — the None-domain raise is `ValueError`, the too-big raise is `MatrixTooLarge`

The live `_resolve_basis_shape` (op.py:329) raises **`ValueError`** when an
operator has NO domain AND no explicit `basis_shape` (a request that is
ill-POSED, not resource-refused); `as_matrix` raises **`MatrixTooLarge`**
(RuntimeError) when `prod(basis_shape) > max_dimension` (a well-posed
request that is resource-REFUSED). These are DIFFERENT classes with
DIFFERENT semantics (§17 A2 rationale is baked into the `MatrixTooLarge`
docstring). The boundary gates MUST pin BOTH separately, and the `pytest.
raises` match strings MUST discriminate — a gate that catches
`(ValueError, MatrixTooLarge)` loosely would pass under a bug that
collapsed the two. This is a Pattern-4 illegal-states check: "un-materialisable
as posed" ≠ "too big to materialise here".

---

## 28. `as_matrix()` L0 CORRECTNESS — the base method (closed-form + storage oracle)

**Claim layer:** foundation + closed-form. **Reference:** hand-built
matrices (`np.diag`, permutation-matrix gather) for the structured leaves;
the STRUCTURALLY-INDEPENDENT `dense_per_material` storage transpose for the
energy leaves. **New file:** `tests/numerics/test_matrix_inverse_operator.py`,
`@pytest.mark.foundation`, `not slow`.

### 28.1 Gate ASM-STRUCTURED — exact against a hand-built matrix

For the two structured leaves whose dense form is closed-form, `as_matrix`
is BIT-EXACT (apply-to-basis on `e_i` gathers a single column with no
accumulation → `array_equal`):

```python
c = np.array([2.0, 5.0, 0.5, 3.0])          # non-uniform, non-±1 (config-blind §0.6)
np.testing.assert_array_equal(
    DiagonalOperator(c).as_matrix(basis_shape=(4,)), np.diag(c),
    err_msg="as_matrix(Diagonal) ≠ diag(c)")
perm = np.array([1, 2, 3, 0])               # 4-cycle: NON-symmetric (transpose-blind guard)
np.testing.assert_array_equal(
    PermutationOperator(perm).as_matrix(basis_shape=(4,)), np.eye(4)[:, perm],
    err_msg="as_matrix(Permutation) ≠ the gather matrix")
```

**Config-blindness (§0.6):** `Diagonal` MUST be non-uniform (a uniform `c`
is blind to a per-column scaling bug) and the `Permutation` MUST be a
NON-involution / NON-symmetric cycle (a symmetric matrix is blind to
M-ASM-TRANSPOSE — the row-vs-column assembly swap). The permutation-matrix
convention (`np.eye(N)[:, perm]` vs `[perm]`) is itself the thing the
transpose mutation flips — pin the exact one the apply produces.

### 28.2 Gate ASM-ORACLE — the energy leaves vs `dense_per_material` (structural independence)

**The independence argument (spell it in-gate).** `as_matrix(basis_shape=
(ng,1))` drives `apply` → the C kernel `apply_p0_in_scatter` (a per-column
matvec accumulation over `g'`); `dense_per_material()[mid]` reads
`sig_s_legendre(mid)[0].T` — a direct storage TRANSPOSE-copy. Neither
computes the other: one is the apply path, one is the storage side. They
agree ONLY if the operator's apply is faithful to its stored cross-sections
— exactly the L0 claim. (This is the §35 note the `dense_per_material`
docstring-honesty fix records: `dense_per_material` is the ORACLE side, not
a production consumer.)

```python
mat_xs = MaterialMesh.from_materials({0: mix_2g_asymmetric}).material_xs_field()
for op in (IsotropicScattering(mat_xs), IsotropicN2N(mat_xs)):
    got = op.as_matrix(basis_shape=(mix.ng, 1))
    ref = op.dense_per_material()[0]         # the single meshless material
    np.testing.assert_allclose(got, ref, rtol=0, atol=4 * np.finfo(float).eps * abs(ref).max(),
        err_msg=f"{type(op).__name__}: apply-to-basis ≠ storage transpose")
# the OperatorSum path — C − K_iso materialises as the difference of the two
loss = MultiplicationOperator.from_mesh(mat_xs.total_cross_section_field, mat_xs.mesh) \
       - (IsotropicScattering(mat_xs) + IsotropicN2N(mat_xs))
ref_sum = np.diag(sig_t) - (sig_s0 + 2.0 * sig_2).T      # fused storage oracle
np.testing.assert_allclose(loss.as_matrix(basis_shape=(mix.ng, 1)), ref_sum, atol=1e-12)
```

**Config-blindness (§0.6 / Mode-6):** the mixture MUST be **≥2G with an
ASYMMETRIC `SigS`** (`Σ_s[g→g'] ≠ Σ_s[g'→g]`) — a 1G or symmetric-`SigS`
mixture makes `M = Mᵀ`, so the transpose-assembly mutation (M-ASM-TRANSPOSE)
is INVISIBLE and the `.T` in `dense_per_material` is un-pinned (Mode-6, the
`SigSᵀ` convention-drift class). Reuse the `homo_2eg_n2n` asymmetric case
(it already de-vacuums the n2n term — see `test_homogeneous.py:57`).

### 28.3 Gate ASM-COLUMN — the C-order flat-index ↔ column-j convention

The apply-to-basis loop's load-bearing convention: flat index `i`
(C-order over `basis_shape`) → column `j = i`, with the output raveled
C-order into the rows. Pin it with a NON-symmetric operator on a
MULTI-element basis so a C↔F order swap in EITHER the basis enumeration or
the output ravel reddens:

```python
# a test-local (ng,1)-carrier operator whose (i,j) matrix entry is DISTINCT
# per cell (M[i,j] = 10*i + j) so any index transposition is O(1)-visible
A = _IndexStampOperator(shape=(3, 1))        # apply(e_ij) column = the stamped col
np.testing.assert_array_equal(
    A.as_matrix(basis_shape=(3, 1)), _expected_stamp((3, 1)),
    err_msg="as_matrix column/row convention (C-order) is transposed")
```

### 28.4 Gate ASM-RECTANGULAR — the output dimension emerges from `apply`

A test-local operator mapping `(2,) → (3,)` (a fixed `(3,2)` action)
materialises to a `(3,2)` matrix — proving `as_matrix` does NOT assume
square and reads the output shape from `apply` itself:

```python
M = _RNG.standard_normal((3, 2))
A = _DenseActionOperator(M)                   # apply(x) = M @ x, (2,)→(3,)
got = A.as_matrix(basis_shape=(2,))
assert got.shape == (3, 2), f"rectangular as_matrix shape {got.shape} ≠ (3,2)"
np.testing.assert_allclose(got, M, atol=1e-14, err_msg="rectangular as_matrix ≠ M")
```

---

## 29. BASIS-SHAPE RESOLUTION + `MatrixTooLarge` BOUNDARY

**Claim layer:** foundation (the resolution rule + the resource gate).
**Reference:** the `_resolve_basis_shape` single-source contract (§27.C).

### 29.1 Gate ASM-RESOLVE — the three resolution arms

```python
def test_basis_shape_explicit_wins():
    # explicit basis_shape overrides even a domain-carrying operator
    A = _SpaceCarrying(shape=(4,))           # domain.shape == (4,)
    assert A.as_matrix(basis_shape=(4,)).shape == (4, 4)

def test_basis_shape_domain_default():
    # domain-carrying operator, NO explicit basis_shape → uses domain.shape
    A = _SpaceCarrying(shape=(3,))
    assert A.as_matrix().shape == (3, 3)     # resolved from A.domain.shape

def test_basis_shape_none_domain_raises_valueerror():
    # NEITHER domain NOR basis_shape → ValueError naming BOTH remedies (§27.C)
    with pytest.raises(ValueError, match="basis_shape|domain"):
        DiagonalOperator(np.array([1.0, 2.0])).as_matrix()   # bare ndarray, domain=None
```

**Config-blindness:** the None-domain fixture MUST be a genuine
domain-`None` leaf (a bare-ndarray `DiagonalOperator` — verified
`domain is None`; ALL numerics leaves are `domain=None`), NOT one that
happens to carry a space. The domain-default fixture MUST carry a `domain`
whose `.shape` differs from any hard-coded default so the resolution is
proven to READ the space.

**Honest-scope (verified):** EVERY ndarray-carrier numerics leaf is
`domain=None` (Diagonal/Permutation/Identity/energy leaves all return
`None`), so the **domain-default resolution arm has NO production
exerciser in step 5** — it is reached only by the test-local
`_SpaceCarrying` fixture (which returns a real `FunctionSpace` with a
`.shape`). The arm becomes production-live when a domain-carrying
ndarray operator appears, or when the FullField `as_matrix` carve lands
(both future). The RULE (`_resolve_basis_shape`) is tested faithfully by
the test-local op; do NOT force a production domain-carrier now.

### 29.2 Gate ASM-TOOLARGE — the resource boundary (at / above / per-call)

The size gate is `prod(basis_shape) > max_dimension` (strict `>` — AT the
threshold PASSES). Pin all three edges + the message quality:

```python
def test_at_threshold_materialises():                 # designed-GREEN control
    A = _SpaceCarrying(shape=(4,))
    A.as_matrix(basis_shape=(4,), max_dimension=4)     # prod==4==max → OK, no raise

def test_one_above_threshold_raises_matrix_too_large():
    A = _SpaceCarrying(shape=(5,))
    with pytest.raises(MatrixTooLarge, match="dimension|4|materiali"):
        A.as_matrix(basis_shape=(5,), max_dimension=4)  # prod 5 > 4 → refuse

def test_per_call_max_dimension_override():
    A = _SpaceCarrying(shape=(5000,))                   # exceeds the DEFAULT 4096
    with pytest.raises(MatrixTooLarge):
        A.as_matrix(basis_shape=(5000,))                # default gate refuses
    A.as_matrix(basis_shape=(5000,), max_dimension=5000)  # per-call lift → OK
```

**Mode-8:** `pytest.raises` (a context manager, not a bare `assert`) fires
under `-O`. **Boundary discipline (§0.6):** the AT-threshold designed-green
control is MANDATORY — a lone `one_above` raises-gate cannot distinguish `>`
from `>=` (an off-by-one gate mutation, M-ASM-GATE-OFFBYONE). The
`ValueError` (§29.1) and `MatrixTooLarge` (§29.2) matches MUST be
class-discriminating (§27.C) — never a shared `(ValueError, RuntimeError)`.

---

## 30. `MatrixInverseOperator` — the 4th sibling (universal contract + M-invariants)

**Claim layer:** foundation + closed-form. **Reference:** `np.linalg.solve`
/ `np.linalg.inv` of the materialised inner (structurally independent of the
LU path? — see the note), at **machine·cond** grain. **Home:**
`tests/numerics/test_matrix_inverse_operator.py`.

### 30.0 The per-type tolerance row (extends §12.0 / §20.0)

| Advertiser | `.inverse()` from | `apply` realization | round-trip tol | `(A⁻¹)⁻¹` form |
|---|---|---|---|---|
| `MatrixInverseOperator` | DIRECT construction (NOT a factory `.inverse()` in step 5 — #138/CP defer) | `lu_solve(lu, q_flat).reshape(basis_shape)` | `atol = K·eps·cond(A)` (**machine·cond**, NOT nulp, NOT driver-tol) | **`is inner`** (object id, via mixin) |

**⚠ Independence caveat (state it).** The M-materialise / I1 anchor
compares `MatrixInverseOperator(A)` against `np.linalg.solve(A_dense, q)`.
When `A_dense = A.as_matrix()` and the anchor is `np.linalg.solve(A.as_matrix(),
q)`, BOTH ride the SAME `as_matrix()` — so that comparison is a LU-vs-LU
self-consistency (necessary, not sufficient). The STRUCTURALLY-INDEPENDENT
anchor is a **hand-built** dense reference (for a `DiagonalOperator` inner:
`A_hand = np.diag(c)`, `inv = np.diag(1/c)`) whose matrix does NOT come from
`as_matrix`. Every §30 value gate MUST anchor on the hand-built reference,
not on `np.linalg.solve(A.as_matrix(), ·)` alone (§27.A precision-grain +
L2 two-anchor).

### 30.1 Gate MINV-I1 — round-trip both ways + hand-built closed-form anchor

```python
c = np.array([2.0, 5.0, 0.5, 3.0])           # non-uniform (cond = 10)
A = DiagonalOperator(c)
minv = MatrixInverseOperator(A, basis_shape=(4,))
q = _RNG.standard_normal(4)
cond = c.max() / c.min()
tol = 32 * np.finfo(float).eps * cond
np.testing.assert_allclose(minv.apply(A.apply(q)), q, atol=tol,
    err_msg="MINV: A⁻¹(Ax) ≠ x")
np.testing.assert_allclose(A.apply(minv.apply(q)), q, atol=tol,
    err_msg="MINV: A(A⁻¹x) ≠ x")
np.testing.assert_allclose(minv.apply(q), q / c, atol=tol,   # HAND-built closed form
    err_msg="MINV ≠ its closed-form inverse (diag(1/c))")
```

### 30.2 Gate MINV-I2 — involution by OBJECT IDENTITY

```python
A = DiagonalOperator(c)
minv = MatrixInverseOperator(A)
assert minv.inverse() is A          # mixin inverse()→inner, object identity
```

### 30.3 Gate MINV-MATERIALISE — the name-earner at MACHINE precision + the Green CONTRAST (§27.A)

```python
A = DiagonalOperator(c)
minv = MatrixInverseOperator(A)
n = 4
Ainv_mat, A_mat = minv.as_matrix(), A.as_matrix(basis_shape=(4,))
mach = 64 * np.finfo(float).eps * (c.max() / c.min())
np.testing.assert_allclose(Ainv_mat @ A_mat, np.eye(n), atol=mach,   # both ways
    err_msg="M-materialise: A⁻¹·A ≠ I at machine precision")
np.testing.assert_allclose(A_mat @ Ainv_mat, np.eye(n), atol=mach)
np.testing.assert_allclose(Ainv_mat, np.diag(1.0 / c), atol=mach,    # HAND anchor
    err_msg="MINV.as_matrix() ≠ diag(1/c)")
# DISTINGUISHING contrast (§27.A): the iterative Green does NOT meet the
# machine grain — a generic inverse satisfying M-materialise only to
# driver-tol is NOT a MatrixInverseOperator.
green = (DiagonalOperator(c) - ScaledOperator(0.3, PermutationOperator(np.array([1,2,3,0])))).inverse()
green_mat = green.as_matrix(basis_shape=(4,))           # inherited base apply-to-basis
A2_mat = (DiagonalOperator(c) - ScaledOperator(0.3, PermutationOperator(np.array([1,2,3,0])))).as_matrix(basis_shape=(4,))
resid = np.abs(green_mat @ A2_mat - np.eye(4)).max()
assert resid > mach, "the Green as_matrix met MACHINE grain — the name-earner is not distinguishing"
assert resid < 10 * green.tol, "sanity: the Green still meets driver-tol"
```

### 30.4 Gate MINV-DIRECT — machine·cond residual + seed BIT-IDENTITY (M-direct = seed-independence)

```python
A = DiagonalOperator(c)
minv = MatrixInverseOperator(A)
q = _RNG.standard_normal(4)
x = minv.apply(q)
resid = np.linalg.norm(A.apply(x) - q) / np.linalg.norm(q)
assert resid <= 32 * np.finfo(float).eps * (c.max()/c.min()), \
    f"M-direct residual {resid:.2e} exceeds machine·cond"
# seed-independence — the SAME q under a garbage seed gives a BIT-IDENTICAL result
junk = _RNG.standard_normal(4) * 1e6
np.testing.assert_array_equal(minv.apply(q), minv.apply(q, initial_guess=junk),
    err_msg="MINV.apply consumed the seed — M-direct is NOT seed-independent")
```

**Config-blindness (§0.6):** the seed MUST be a non-zero, non-`q` junk
vector (a zero seed or `seed==q` cannot red M-MINV-SEED-CONSUME — a
consuming apply that warm-started from `q` or `0` would be invisible).
`array_equal` (bit-identity, not `allclose`): the CORRECT `del initial_guess`
gives byte-identical output; a consuming mutation drifts.

### 30.5 Gate MINV-BACKHALF — `solve` = the FORWARD matvec (anchored, NOT tautological)

The mixin's `solve` on the inverse object IS the forward action. Anchor
against the HAND-built forward, never via `minv.inverse().apply` (which is a
definition-tautology, §20.2):

```python
A = DiagonalOperator(c)
minv = MatrixInverseOperator(A)
x = _RNG.standard_normal(4)
np.testing.assert_allclose(minv.solve(x), c * x, rtol=1e-14,   # forward = c ⊙ x
    err_msg="mixin solve ≠ the forward matvec")
```

### 30.6 Gate MINV-GUARDS — ctor refuses non-square, singular, too-big (anti-#11 ± pairs)

```python
def test_minv_non_square_raises():
    A = _DenseActionOperator(_RNG.standard_normal((3, 2)))    # (2,)→(3,), rectangular
    with pytest.raises((ValueError, np.linalg.LinAlgError), match="square"):
        MatrixInverseOperator(A, basis_shape=(2,))

def test_minv_singular_raises_linalgerror():                  # TRUE zero, §0.6
    A = DiagonalOperator(np.array([2.0, 0.0, 0.5]))           # exact-singular matrix
    with pytest.raises(np.linalg.LinAlgError):
        MatrixInverseOperator(A, basis_shape=(3,))            # lu_factor raises naturally

def test_minv_too_large_propagates():
    A = _SpaceCarrying(shape=(5000,))
    with pytest.raises(MatrixTooLarge):
        MatrixInverseOperator(A, basis_shape=(5000,))         # ctor as_matrix refuses

def test_minv_nonsingular_constructs():                       # POSITIVE control
    MatrixInverseOperator(DiagonalOperator(c))                # no raise
```

**Config-blindness:** the singular fixture MUST be a TRUE `0.0` coefficient
(a `1e-300` passes `lu_factor` — same trap as §8/§12.4). **Note the
DIVISION of guards:** non-square is a `MatrixInverseOperator` ctor check;
singular is `scipy.lu_factor`'s NATURAL `LinAlgError` (the design does NOT
add a redundant singular guard — the brief); too-big propagates from
`as_matrix`. Three distinct raise SITES, three distinct classes.

### 30.7 Gate MINV-WITNESS — VALUE-vs-STRUCTURE: inverts what Green REFUSES (§27.B)

```python
def test_matrix_inverts_what_green_refuses():
    """The §3 strategy-override witness (ndarray analog of (−S)+(L+C), §27.B):
    a sum whose LEADING term is non-invertible — Green REFUSES at construction,
    the dense strategy INVERTS the materialised matrix."""
    S_ao = _DenseActionOperator(_off_diagonal_nilpotent(4))   # apply-only, is_invertible=False
    D = DiagonalOperator(c)
    sum_op = (-1.0 * S_ao) + D          # left-spine head = ScaledOperator(-1, S_ao), NOT invertible
    with pytest.raises(MissingCapability, match="canonical ordering|invert"):
        sum_op.inverse()                # Green REFUSES (leading non-invertible)
    minv = MatrixInverseOperator(sum_op, basis_shape=(4,))    # materialise D − S_ao, invert
    q = _RNG.standard_normal(4)
    A_hand = np.diag(c) - _off_diagonal_nilpotent(4)          # structurally independent
    np.testing.assert_allclose(minv.apply(q), np.linalg.solve(A_hand, q), atol=1e-12,
        err_msg="MatrixInverse did not invert the leading-non-invertible sum")
```

**Config-blindness:** `D − S_ao` MUST be genuinely non-singular as a MATRIX
(pick `S_ao` strictly off-diagonal / nilpotent so `det(D − S_ao) ≈ ∏c ≠ 0`),
AND `S_ao.is_invertible` MUST be `False` (else Green would NOT refuse and the
asymmetry is vacuous). This is the gate that proves the DESIGN CLAIM "the
matrix realization reads VALUES not STRUCTURE".

### 30.8 Gate MINV-FAITH — the 4th inverse object is faithful (extend §12.3 home)

EXTEND `tests/numerics/test_operator_capability_predicates.py::
test_inverse_objects_are_faithful`:

```python
minv = MatrixInverseOperator(DiagonalOperator(_C))
_assert_faithful(minv)                       # {APPLY, SOLVE}, is_inv True, is_adj False
assert minv.inverse() is DiagonalOperator... # the inner, object identity
```

### 30.9 Gate MINV-STATIC — the 4th-sibling seeded-apply static pins (extend §19.2)

EXTEND `tests/sn/operators/test_operators_apply_typed.py::
_inverse_family_seeded_apply_static_pins` with a `matrix:
MatrixInverseOperator` parameter:

```python
_ = matrix.apply(arr, initial_guess=arr)     # accepts the canonical kwarg (mixin abstract)
d: SupportsSeededApply = matrix              # structural conformance (no reportAssignmentType)
del d
```

**Teeth (M-MINV-KWARG, CLI-pyright, §34):** drop `initial_guess` from
`MatrixInverseOperator.apply` → `npx pyright orpheus/` reports
`reportIncompatibleMethodOverride` vs the mixin's abstract sig. Same
structural proof as M-MIXIN-KWARG (§19.2), now on the 4th sibling.

### 30.10 Gate MINV-REGISTRY — participation in the universal family (extend `test_inverse_universal.py`)

`MatrixInverseOperator` is NOT a factory `.inverse()` dispatch target in
step 5 (direct construction only — the `A.inverse() → MatrixInverseOperator`
routing is #138/CP). So it joins the universal net as a DIRECTLY-constructed
registry row, NOT a `type(A.inverse()) is MatrixInverseOperator` dispatch
pin:

```python
# EXTEND test_inverse_universal.py — a direct-construction row
def test_matrix_inverse_operator_universal_roundtrip():
    A = DiagonalOperator(_C7)
    minv = MatrixInverseOperator(A)
    np.testing.assert_allclose(minv.apply(A.apply(_X7)), _X7, atol=1e-12)
    np.testing.assert_allclose(minv.apply(_X7), _X7 / _C7, atol=1e-12)   # closed form
    assert minv.inverse() is A                                           # involution
```

State in-gate: the `.inverse()`-FACTORY row (`type(...) is
MatrixInverseOperator`) lands with #138 — do NOT add a phantom dispatch pin
now.

---

## 31. #285 — `OperatorProduct.inverse()` → `InverseOperator(self)`

**Claim layer:** foundation (behaviour-preserving VALUE; a NEW seeded-apply
surface + a STRENGTHENED involution). **Reference:** the composition value
`b.solve(a.solve(q))` (bit-identity inheritance) + object identity.
**Home:** `tests/numerics/test_inverse_universal.py` (extend the product
rows) + `tests/numerics/test_operator.py` (the product-inverse pins).

### 31.1 The migration inventory (grep-verified — the surface is SMALL)

`.solve` on `OperatorProduct` is UNCHANGED (`b.solve(a.solve(q))`); #285
touches ONLY `.inverse()`. Grep result:

| Existing pin | Today | Under #285 | Disposition |
|---|---|---|---|
| `test_operator.py::test_product_solve_reverses_order` (:293) | `p.solve(vector)` value | `.solve` UNCHANGED | **STAYS GREEN** (no touch) |
| `test_operator.py::test_dunder_matmul_returns_product` (:367) | `isinstance(A@A, OperatorProduct)` | forward `@` construction | **STAYS** (not `.inverse()`) |
| `test_diagonal_operator.py::test_product_of_two_diagonals` (:133) | `isinstance(D1@D2, OperatorProduct)` | forward `@` | **STAYS** |
| `test_inverse_universal.py::test_i2_product_reversed` (:127) | VALUE `(D@P).inverse().apply(q) == P⁻¹(D⁻¹q)` | `InverseOperator(D@P).apply(q) = (D@P).solve(q) = P.solve(D.solve(q))` — SAME value, bit-id | **STAYS GREEN** — add a companion structure pin (§31.3) |
| `test_operator_capability_predicates.py::test_inverse_objects_are_faithful` (:156) | `_assert_faithful((D@P).inverse())` | now `InverseOperator(D@P)` — {APPLY,SOLVE}, faithful | **STAYS GREEN** |

**No test asserts `type((A@B).inverse()) is OperatorProduct` or compares to
`OperatorProduct(b.inverse(), a.inverse())`** (grep-confirmed) — so the
return-type change breaks NO existing structure pin. State NONE-FOUND so a
later reader doesn't hunt.

### 31.2 Gate PROD-285-REPRO — the seeded call is now ACCEPTED (the #285 witness)

```python
def test_product_inverse_accepts_seeded_apply():
    """#285: a composed inverse now carries the canonical seeded-apply kwarg.
    TypeError today (raw OperatorProduct.apply is positional-only); accepted
    after (InverseOperator's mixin apply(x, /, *, initial_guess=None))."""
    D, P = DiagonalOperator(_C7), PermutationOperator(_P7)
    inv = (D @ P).inverse()
    q = _RNG.standard_normal(7)
    out = inv.apply(q, initial_guess=q)      # <-- TypeError pre-#285, accepted post
    np.testing.assert_array_equal(out, inv.apply(q))   # seed accept-and-ignore
```

### 31.3 Gate PROD-285-VALUE + INVOLUTION — bit-id composition path + strengthened involution

```python
def test_product_inverse_value_and_involution():
    D, P = DiagonalOperator(_C7), PermutationOperator(_P7)
    prod = D @ P
    inv = prod.inverse()
    assert type(inv) is InverseOperator          # NEW: the generic member, not a raw product
    q = _RNG.standard_normal(7)
    np.testing.assert_array_equal(               # bit-id composition: b.solve(a.solve(q))
        inv.apply(q), P.solve(D.solve(q)),
        err_msg="#285 product-inverse value ≠ b.solve(a.solve(q))")
    assert inv.inverse() is prod                 # involution STRENGTHENS to object identity
```

**Config-blindness (§0.6):** the factors MUST be NON-COMMUTING (`Diagonal @
Permutation`) — a commuting `D1 @ D2` cannot red M-PROD-FACTORORDER (the
`b.solve`/`a.solve` order-swap), exactly as §12.2 required.

### 31.4 Gate PROD-285-NEGATIVE + ALGEBRA-CLOSED preservation

```python
def test_non_invertible_product_still_raises():
    with pytest.raises(MissingCapability, match="invertible|both"):
        (DiagonalOperator(_CZ) @ IdentityOperator()).inverse()   # a factor singular
def test_algebra_closed_inverses_unchanged():
    # a perm's inverse IS a perm (the OTHER kind — first-class forward, NOT a wrap);
    # #285 does NOT route these through InverseOperator (documented, not gated for seed)
    assert type(PermutationOperator(_P7).inverse()) is PermutationOperator
    assert type(IdentityOperator().inverse()) is IdentityOperator
    assert type(ScaledOperator(2.0, DiagonalOperator(_C7)).inverse()) is ScaledOperator
```

---

## 32. HOMOGENEOUS REWIRE — `_as_dense` RETIRES → `as_matrix(basis_shape=(ng,1))`

**Claim layer:** foundation (bit-identical output) + closed-form (the
LANDED k∞ anchor). **Reference:** the EXISTING homogeneous gates (they
already pin k∞ + flux + rates TIGHTLY) + an inlined local `_as_dense`-loop
oracle for the equivalence.

### 32.1 The existing pins are STRONG — lean on them (inventory)

The homogeneous suite ALREADY pins the values the rewire must preserve
(`tests/homogeneous/test_homogeneous.py`):

| Existing gate | Pins | Rewire disposition |
|---|---|---|
| `test_kinf_exact` (:62) | k∞ vs SymPy `case.k_inf`, 1e-12, 4 cases | **STAYS GREEN** (structurally-independent anchor; MUST NOT move) |
| `test_post_solve_production_rate_is_100` (:80) | flux normalisation | **STAYS GREEN** |
| `test_assemble_loss_matrix_matches_fused_oracle` (:118) | `_assemble_loss_matrix` ≈ fused `diag(Σt)−(Σs0+2Σ2)ᵀ`, 1e-12 | **CALL-SITE MIGRATES** (it imports `_assemble_loss_matrix`, which drops `_as_dense`) — value pin unchanged |
| `test_kinf_is_the_direct_eigenvalue_of_the_assembled_pair` (:257) | solver k∞ == `direct_eigenvalue(A,F)[0]` BYTE-equal; imports `_as_dense` directly (:276) | **CALL-SITE MIGRATES** (`_as_dense` import → `.as_matrix(basis_shape=(ng,1))`); byte-identity pin unchanged |
| `test_rates_via_integrated_reaction_rate_are_bit_identical` (:299) | rate rerouting bit-id | **STAYS GREEN** |

**So the homogeneous bit-identity is ALREADY gated** — the rewire's job is
to keep `test_kinf_exact` + `test_kinf_is_the_direct_eigenvalue...` GREEN
byte-for-byte. Do NOT add a redundant fresh k∞ pin; migrate the two
`_as_dense`-importing call sites and let the LANDED pins guard the values.

### 32.2 Gate HOMOG-EQUIV — `as_matrix(basis_shape=(ng,1))` ≡ the retired `_as_dense` (local oracle)

The equivalence gate INLINES the retired loop as a local oracle (the
fuller-view-oracle discipline — the retired code becomes a test-local
reference, then the production path is proven identical):

```python
def test_as_matrix_equals_retired_as_dense_loop():
    """The promoted as_matrix reproduces the retired _as_dense apply-to-basis
    loop BYTE-for-byte (the retirement is a pure relocation, not a re-derivation)."""
    mat_xs = MaterialMesh.from_materials({0: mix_2g_n2n}).material_xs_field()
    ng = mix.ng
    loss = MultiplicationOperator.from_mesh(mat_xs.total_cross_section_field, mat_xs.mesh) \
           - (IsotropicScattering(mat_xs) + IsotropicN2N(mat_xs))
    # local oracle = the retired loop, verbatim
    cols = []
    for i in range(ng):
        e_i = np.zeros((ng, 1)); e_i[i, 0] = 1.0
        cols.append(np.asarray(loss.apply(e_i)).ravel())
    oracle = np.column_stack(cols)
    np.testing.assert_array_equal(          # BYTE-identical (same loop, promoted)
        loss.as_matrix(basis_shape=(ng, 1)), oracle,
        err_msg="as_matrix diverged from the retired _as_dense loop")
```

This gate lives with the retirement; delete it if/when the local-oracle
duplication is judged pure noise (per the retirement-audit — but keep it
through the merge cycle so a later `as_matrix` refactor cannot silently
move the homogeneous basis columns).

### 32.3 The `dense_per_material` docstring-honesty fix (doc-only, no runtime gate)

Cut 6 is doc-only — `dense_per_material` (isotropic_scattering.py:150-157,
:220) claims to be "the `as_dense` consumption mode for the LHS-fold solvers
(diffusion / homogeneous)". That is FALSE (homogeneous used `_as_dense`;
grep-confirmed ZERO production consumers of `dense_per_material`). The
rewrite states its TRUE role: the STRUCTURALLY-INDEPENDENT storage-side
oracle (`sig_s_legendre(mid)[0].T`) that the §28.2 `as_matrix` gate
cross-checks against. **No runtime gate** (behaviour unchanged); the §28.2
ASM-ORACLE gate IS the thing that keeps `dense_per_material` honest (it uses
it as the oracle, so a drift reddens ASM-ORACLE). Verify via the Sphinx `-W`
build (the docstring is a `:math:` block) + a grep in the closeout proving
no production caller was mis-retired.

---

## 33. MIXIN BOUND RELAXATION — the no-op proof

**Claim layer:** foundation (pure widening — no value/behaviour change).
**Reference:** the LANDED family suites (inheritance) + pyright.

Cut 2 relaxes `_ForwardT`'s bound `_InvertibleForward → LinearOperator`.
This is a PURE WIDENING (the back-half consumes only
domain/codomain/apply); the three landed siblings (`InverseOperator`,
`SweepOperator`, `GreenOperator`) are unaffected. The no-op is proven by
the EXISTING family gates staying green UNCHANGED:

| Regression net (stays green, UNCHANGED) | Proves |
|---|---|
| `test_inverse_universal.py` (I1/I2/faithfulness) | leaf inverses unchanged |
| `test_green_operator.py` (all §20–§23 gates) | Green unchanged |
| `test_inverse_operator_equivalence.py` + the seed spy | Sweep unchanged |
| `test_operators_apply_typed.py::_inverse_family_seeded_apply_static_pins` | the abstract sig still LSP-enforced |

**Pyright expectation:** the ratchet stays EXACTLY at baseline **148**
(numerics 21 / sn 102 / transport 25) — the relaxation ADDS no error (it
widens an accepted type) and `MatrixInverseOperator` + `as_matrix` introduce
NO new `# type: ignore` (anti-#19; `cast` OK per the `_seeded_inverse`
precedent). A ratchet move is a RED. Verify via CLI `npx pyright orpheus/`
(the oracle), not the LSP.

---

## 34. MUTATION BANK — Step 5 (each REDs a NAMED gate, under `-O`, monkeypatch only)

In-process monkeypatch ONLY (uncommitted edits — NEVER `git checkout`,
process-discipline.md; §27 notes the file is LIVE-edited). Every mutation
names its ACTIVATING config (§0.6 applies to mutations).

| ID | Target gate | Mutation | Expected RED | Activating config |
|---|---|---|---|---|
| **M-ASM-TRANSPOSE** | §28.1/§28.2 | `as_matrix` assembles row-major (`np.vstack`) instead of `column_stack` (transpose the matrix) | Diagonal-exact + dense_per_material ORACLE red | NON-symmetric op (asymmetric SigS / 4-cycle perm; symmetric BLIND) |
| **M-ASM-RAVEL** | §28.3/§28.4 | ravel output in F-order (`ravel("F")`) instead of C-order | column-convention + rectangular red | multi-element basis (`(ng,1)`/`(2,)`; scalar `(n,)` blind) |
| **M-ASM-GATE-DELETE** | §29.2 | delete the `prod > max_dimension` check | one-above gate reds (no `MatrixTooLarge`) | `prod > max_dimension` |
| **M-ASM-GATE-OFFBYONE** | §29.2 | `>=` instead of `>` in the size gate | AT-threshold designed-green control reds (raises when it should pass) | `prod == max_dimension` EXACTLY |
| **M-ASM-RESOLVE** | §29.1 | `_resolve_basis_shape` ignores explicit `basis_shape`, always uses `domain.shape` | explicit-wins gate reds (wrong shape / None-domain raise) | explicit `basis_shape` ≠ `domain.shape` |
| **M-MINV-LUTRANS** | §30.1/§30.4 | `lu_solve(lu, q, trans=1)` (solves Aᵀx=q) | I1 + M-direct red | NON-symmetric inner (symmetric blind) |
| **M-MINV-SEED-CONSUME** | §30.4 | `apply` warm-starts / adds `initial_guess` (drop `del initial_guess`) | seed bit-identity red | non-zero junk seed ≠ q |
| **M-MINV-ASM-FORWARD** | §30.3 | `as_matrix()` OVERRIDE returns `inner.as_matrix()` (forward) not `lu_solve(lu, eye)` | M-materialise red (A·A=A²≠I) | A ≠ A⁻¹ (non-±1 `c`) |
| **M-MINV-NOGUARD-SQ** | §30.6 | drop the non-square ctor guard | non-square gate: wrong/no clean raise | rectangular inner |
| **M-MINV-KWARG** | §30.9 static | drop `initial_guess` from `MatrixInverseOperator.apply` | `npx pyright` `reportIncompatibleMethodOverride` vs mixin abstract | any (CLI-pyright) |
| **M-MINV-WITNESS** | §30.7 | `MatrixInverseOperator` ctor ADDS an `inner.is_invertible` guard (reads STRUCTURE) | the witness gate reds (refuses the leading-non-invertible sum it SHOULD invert) | `(−S_ao)+D` (leading non-invertible, matrix non-singular) |
| **M-PROD-RAWWRAP** | §31.2/§31.3 | `OperatorProduct.inverse()` returns the raw `OperatorProduct(b.inverse(), a.inverse())` (revert #285) | seeded-call repro reds (TypeError, positional-only) + involution-object-id reds | seeded call / involution |
| **M-PROD-FACTORORDER** | §31.3 | `OperatorProduct.solve` → `a.solve(b.solve(q))` (a/b swap) | product-inverse value red (`≠ b.solve(a.solve(q))`) | NON-commuting `Diagonal @ Permutation` |
| **M-HOMOG-ASDENSE** | §32.2 | (rides M-ASM-TRANSPOSE/RAVEL) — `as_matrix` diverges from the retired loop | HOMOG-EQUIV red + `test_kinf_is_the_direct_eigenvalue...` byte-id red | homo_2eg_n2n |

**Designed-GREEN controls (mark explicitly):**
- §29.2 `test_at_threshold_materialises` — the `prod==max` PASS (pairs with
  one-above; M-ASM-GATE-OFFBYONE reds THIS, proving `>` not `>=`).
- §30.6 `test_minv_nonsingular_constructs` — the POSITIVE ctor control.
- §30.3 the Green-contrast `resid < 10·green.tol` sanity — proving the Green
  IS a working inverse (just to driver-tol), so the machine-grain gap is the
  distinguishing signal, not a broken Green.
- §31.2 `np.testing.assert_array_equal(out, inv.apply(q))` — the seed IS
  accept-and-ignore (a designed-green: the product inverse ignores the seed
  by delegation, exactly like the exact leaves).

**Mode-8 hygiene:** every gate uses `np.testing.*` / `pytest.raises` /
`assert_type` / collected-module `assert` (pytest-rewritten → fires under
`-O`). The ctor guards RAISE (not `assert`) so they survive `-O`.

---

## 35. HOMES, MARKERS, LAND-ORDER, HONEST-SCOPE

### 35.1 Homes + markers

| Home | Contents | Marker |
|---|---|---|
| **NEW** `tests/numerics/test_matrix_inverse_operator.py` | §28 as_matrix L0 (structured/oracle/column/rectangular) + test-local `_IndexStampOperator`/`_DenseActionOperator`/`_SpaceCarrying`; §29 basis-resolution + MatrixTooLarge; §30 MINV I1/I2/materialise/direct/backhalf/guards/witness | `foundation`, not slow |
| EXTEND `tests/numerics/test_inverse_universal.py` | §30.10 MINV direct-construction registry row; §31 the #285 product-inverse repro/value/involution/negative/algebra-closed | `foundation` |
| EXTEND `tests/numerics/test_operator_capability_predicates.py` | §30.8 MINV-object faithfulness | `foundation` |
| EXTEND `tests/sn/operators/test_operators_apply_typed.py` | §30.9 the 4th-sibling seeded-apply `assert_type` + `SupportsSeededApply` conformance | `foundation` (static) |
| MIGRATE `tests/homogeneous/test_homogeneous.py` | §32.1 the two `_as_dense`-importing call sites → `.as_matrix(basis_shape=(ng,1))`; §32.2 ADD the HOMOG-EQUIV local-oracle gate | (keep existing l1/verifies) |

**NO** `verifies(<physics-eq>)` on the Step-5 gates — all software
invariants (foundation; L9). The homogeneous file KEEPS its landed
`l1`+`verifies` markers (the k∞ chain is equation-level); the migrated call
sites do NOT change those markers.

### 35.2 Land-order (test-first — author gate → RED under its mutation → GREEN on landed method)

| # | Cut | Gates that land WITH it | Teeth |
|---|---|---|---|
| **5a** | `as_matrix()` base method + `_resolve_basis_shape` + `MatrixTooLarge` | §28 (structured/oracle/column/rectangular), §29 (resolve/toolarge) | M-ASM-TRANSPOSE/RAVEL/GATE-DELETE/GATE-OFFBYONE/RESOLVE |
| **5b** | mixin bound relaxation | §33 no-op (existing suites green UNCHANGED) + pyright 148 | (ratchet move = RED) |
| **5c** | `MatrixInverseOperator` (direct construction) | §30 MINV I1/I2/materialise/direct/backhalf/guards/witness/faith/static/registry | M-MINV-LUTRANS/SEED-CONSUME/ASM-FORWARD/NOGUARD-SQ/KWARG/WITNESS |
| **5d** | `OperatorProduct.inverse()` → `InverseOperator(self)` (#285) | §31 repro/value/involution/negative/algebra-closed | M-PROD-RAWWRAP, M-PROD-FACTORORDER |
| **5e** | homogeneous rewire + `_as_dense` retirement + `dense_per_material` doc | §32 call-site migration + HOMOG-EQUIV; existing k∞ pins stay byte-green | M-HOMOG-ASDENSE |

### 35.3 Honest-scope declarations (Mode-10 third state + out-of-scope)

- **`as_matrix` on TYPED-CARRIER (FullField) operators is OUT of scope.**
  The apply-to-basis loop builds bare-ndarray `e_i`; a FullField-carrier
  operator's `apply` needs a FullField basis vector. Step 5 gates
  NDARRAY-carrier operators only (Diagonal/Permutation/energy leaves/dense
  test operators). The FullField `(L+C)`/`(−S)+(L+C)` materialisation is a
  future carve (its 3-D sparse-triangular sibling is already deferred,
  taxonomy §5 / sweep_graph.py:66). The §30.7 witness uses the ndarray
  analog for exactly this reason (§27.B).
- **`MatrixInverseOperator` domain/codomain SWAP is exercised-but-unconstrained
  (Mode-10, inherited from Green delta 3).** Every current inner is
  square + endomorphic (`FullField→FullField`) or space-less
  (`domain=None`), so a mixin `domain`/`codomain` swapped-the-wrong-way is
  VALUE-INVISIBLE — no committed gate reds it until a `Codomain ≠ Domain`
  square invertible exists. The `solve→inner.solve` half of that mixin
  mutation IS constrained (§30.5 back-half anchor reds it). Declared per the
  Mode-10 third state; revisit at a non-endomorphic invertible.
- **`MatrixInverseOperator` is NOT a `.inverse()` factory-dispatch target in
  step 5** — direct construction only. The `A.inverse() →
  MatrixInverseOperator` routing (and the `K = MatrixInverseOperator(loss) @
  F` normal form) land with #138 / CP. Do NOT add a `type(A.inverse()) is
  MatrixInverseOperator` dispatch pin now (it would be a phantom).

**Done when:** for EVERY Step-5 gate a named mutation reds it under `-O`
(§34); `as_matrix` is exact vs hand-built + the structurally-independent
`dense_per_material` oracle, with the C-order column convention pinned on a
NON-symmetric op; the `MatrixTooLarge`/`ValueError` boundary is
class-discriminated with the AT-threshold designed-green control;
`MatrixInverseOperator` earns its name at MACHINE·cond grain (M-materialise
+ M-direct, with the Green driver-tol CONTRAST proving it distinguishing)
and its seed-independence is byte-pinned; the §30.7 witness proves it
inverts the leading-non-invertible sum Green refuses; #285's product inverse
accepts the seeded kwarg (repro flips TypeError→accepted) with the
strengthened object-identity involution and the value bit-identical to
`b.solve(a.solve(q))`; the homogeneous k∞/flux pins stay BYTE-green through
the `_as_dense` retirement with the HOMOG-EQUIV local-oracle guarding the
relocation; the mixin relaxation is a proven no-op (existing suites
unchanged, pyright 148); `dense_per_material`'s docstring tells the truth.

---

# PART V — STEP 6: the frozenset RETIREMENT (both axes) — Design "C + B" (2026-07-02)

> **Why this Part exists.** Steps 1–5 (PARTS II–IV) LANDED the whole
> inverse family (`InverseOperator`/`SweepOperator`/`GreenOperator`/
> `MatrixInverseOperator`), the derived `is_invertible`/`is_adjointable`
> properties, the static `SupportsInverse`/`SupportsAdjoint` Protocols,
> and the coexistence-era faithfulness scaffold (`tests/_harness/predicates.py`
> `assert_capability_faithful`, called from the numerics + SN/transport
> enumerations). Step 6 (taxonomy §12 step 6; carve §0.6; spec §0′/§0.6)
> RETIRES the `capabilities: frozenset[str]` — CAP_APPLY + CAP_SOLVE +
> CAP_APPLY_TRANSPOSE + `MissingCapability` + `_has` — for BOTH axes at
> once, under the user-locked "C + B" design (2026-07-02). This Part is
> the verification net for that retirement, main-agent-direct (surgical).
> It EXTENDS PARTS I–IV; the §0.5 gate-liveness / §0.6 config-blindness
> standing disciplines and the §7/§12.6/§25/§34 mutation banks remain in
> force. It **supersedes** §4a/§4e (whose dispositions predate steps 4–5)
> and §2.3 (whose scaffold DELETES here — §36 is its permanent successor).
> Line anchors are VERIFIED against the live tree @ HEAD `d82fa77`
> (`git diff main...HEAD` clean beyond steps 1–5); they DRIFT under the
> live carve — target BEHAVIOUR.

**The production state this Part gates (verified @ `d82fa77`).** The
derived properties and static Protocols already EXIST and coexist with
the frozenset (steps 0–5 landed them). Step 6 does NOT introduce them —
it deletes the redundant advertisement they duplicate. The nine cuts:

1. **Mint `NotInvertible(TypeError)` + `MissingAdjoint(TypeError)`**
   (operator.py) → REPLACE `MissingCapability` (op.py:269). TypeError
   parentage for except-clause continuity; **verified ZERO production
   `except MissingCapability` clauses** (grep on `d82fa77`), so no
   `except` breaks. (Test-side `pytest.raises(MissingCapability)` — 48
   lines / 11 files — migrate by axis; §37.)
2. **Free narrowing TypeGuards** `invertible(op) -> TypeGuard[SupportsInverse]`
   / `adjointable(op) -> TypeGuard[SupportsAdjoint]` (PEP-647, NOT TypeIs —
   value-dependent, so the False branch makes no static claim); free
   functions because a property cannot narrow `self`.
3. **Promote `SupportsInverse`/`SupportsAdjoint`** (op.py:808/822) to
   EXTEND `LinearOperator[Domain, Codomain]` (the narrowed branch retains
   `apply`/`H`/…). `SupportsInverse` declares `inverse()`; `SupportsAdjoint`
   declares `apply_transpose(...)`.
4. **`.adjoint()`/`.H` (op.py:644/674) gains the EAGER gate** —
   `if not adjointable(self): raise MissingAdjoint(...)` then build the
   wrapper; the narrowing types `_AdjointOperator(self)`'s ctor so the
   op.py:934 ignore dies. The LAZY gate inside `_AdjointOperator.apply`
   (op.py:911–916) DELETES. **Raise site moves `.apply` → `.H`** (a real
   behavior change — §38, M-ADJ-EAGER).
5. **`inverse()` is NOT declared on the `LinearOperator` base** —
   structural non-invertibility (Zero, masks, RankOne, streaming `L`,
   energy S/F leaves, `B`) = method ABSENT (pyright error at call sites,
   `AttributeError` only if forced). Value-dependent non-invertibility
   (zero-coeff Diagonal/Mult, sum-with-non-invertible-head) = method
   EXISTS, raises `NotInvertible` eagerly.
6. **Composer law bodies guard-narrow** instead of `# type: ignore` —
   `OperatorSum.apply_transpose` (op.py:1040), `OperatorProduct.apply_transpose`
   (:1180), `ScaledOperator.apply_transpose` (:1287), `TensorProduct`
   (:2040), `SumOfTP` (:2153), `RankOne` (:2571 region): `if not
   adjointable(self.a) or not adjointable(self.b): raise MissingAdjoint(...)`.
   All 4 `cast(SupportsInverse, …)` sites (op.py:1310 Scaled.inverse,
   :2076 TP.inverse per-factor loop, iteration.py:231 [QUOTED string],
   _bound_compat.py:143) → guard+narrow.
7. **`solve` pruned to native realizations** (Design B): DELETE `solve`
   on `OperatorSum` (:1097, "transitional" line executed — the inverse
   action IS the GreenOperator), `IdentityOperator` (:1334),
   `PermutationOperator` (:1753), `ScaledOperator` (:1282),
   `TensorProductOperator` (:2043), `_BoundBoundaryOperator` (_bound_compat.py:161).
   KEEP `solve` on `DiagonalOperator` (:2397, gate→NotInvertible),
   `MultiplicationOperator` (:387/:231→NotInvertible), `OperatorProduct`
   (:1174, RE-ROUTED — §40), the sweep composites (`InvertibleOperator`
   streaming.py:844, `_ScheduledInvertibleOperator` :184), the mixin
   (`InverseWrapMixin.solve` :1569 = inner.apply — the inverse-family
   un-invert face).
8. **Retire the frozenset + tags + `_has`** — `capabilities` field
   (op.py:378) + every leaf/composite class-attr / `default_factory` /
   `field` / `ClassVar` / `@property` set + `CAP_APPLY`/`CAP_SOLVE`/
   `CAP_APPLY_TRANSPOSE` (:136-138) + `_has` (op.py:836 AND the copy
   iteration.py:204) + `MissingCapability` + all solve/apply_transpose
   `# type: ignore`s. Composition apply-guards KEEP but raise plain
   `TypeError` (SAME messages) — Sum :977-986, Product :1130-1139, Scaled
   :1243-1247, TP :1998-2001; SI/Krylov apply-guards (iteration.py:507-523/
   :757-768) → `callable(getattr(op,"apply",None))` raising TypeError
   (`match="apply"` preserved). The repr caps-read (op.py:787) needs a
   successor decision — FLAGGED (§44), not decided.
9. **`tests/_harness/predicates.py` DELETES** with the frozenset; §36 is
   its permanent successor.

**Claim layer (§1.5 gate).** Step 6 is **behaviour-preserving** — no
value moves (the "iff both" laws already live in the method bodies;
§0′). Every gate is **foundation software-invariant** (the two-axis
faithfulness contract, the re-route bit-identity) + rides the LANDED
closed-form anchors (k∞, k_inf-reference) for independence. **No gate
here makes an eigenvalue or MMS claim** (anti-#2). The one genuine
VALUE assertion is the `OperatorProduct.solve` re-route (§40) — gated
bit-identity (direct rows) / iterative-tol (the Green-factor row)
against the pre-carve baseline + a structurally-independent dense
anchor.

---

## 36. KEYSTONE v2 — the PERMANENT two-axis faithfulness successor (replaces the §2.3 scaffold)

**The §2.3 scaffold (`assert_capability_faithful`, `is_X ≡ CAP_X∈caps`)
DELETES with the frozenset — it can only run while `op.capabilities`
exists.** Its permanent successor asserts the inverse/adjoint contract
DIRECTLY, referencing NO frozenset, so it survives the retirement and
becomes the standing net. It is the load-bearing gate of this Part:
every other cut rides on it staying green.

**Claim layer:** foundation (the structural inverse/adjoint contract).
**Reference:** the contract itself (no old-≡-new — the frozenset is
gone). Independence for the VALUE that `is_invertible`/`is_adjointable`
report rides the LANDED I1 round-trip (`test_inverse_universal.py`) +
the closed-form k_inf/keff anchors — a wrong predicate that still lets
`.inverse()` return would be caught by I1's round-trip, not here; §36
pins the CONTRACT (returns-vs-raises-vs-absent), not the inverse value.

**Home.** The scaffold's two enumeration callers MIGRATE (they are not
deleted — coding-standards retirement-means-test-migration): strip every
`_assert_faithful(op)` call (it reads the deleted `capabilities`), KEEP
the direct axis assertions already present, and route through a NEW
shared helper `tests/_harness/predicates.py::assert_inverse_adjoint_contract`
(the ONE-place contract, Cardinal Rule 2 — the file is REWRITTEN, not
deleted). The helper is `-O`-safe (explicit `raise`, not `assert` — it
is an un-collected helper module; the §2.3 predicates.py docstring
already documents this Mode-8 rule).

- `tests/numerics/test_operator_capability_predicates.py` — numerics
  leaves/composites (the `_LEAVES` list :62, the recursive-sum/product
  tests, the inverse-OBJECT faithfulness :145).
- `tests/sn/operators/test_capability_survival.py` — the SN/transport
  advertisers (`TestPredicateFaithfulness` :203, the boundary-law-wrapper
  delegation :268, the `:277-278` literal-CAP reads).

### 36.1 The shared helper (rewritten predicates.py) — the four legs

```python
# tests/_harness/predicates.py  (REWRITTEN — no frozenset read; -O-safe raise)
_STRUCTURAL_ABSENT = "structural-absent"   # inverse() NOT declared → AttributeError only if forced
_VALUE_RAISE       = "value-raise"         # inverse() declared, raises NotInvertible when singular
_INVERTIBLE        = "invertible"          # inverse() returns

def assert_inverse_adjoint_contract(op, *, invertible, adjointable, inverse_contract):
    """The PERMANENT successor to assert_capability_faithful (spec §36).

    invertible : bool          expected op.is_invertible
    adjointable: bool          expected op.is_adjointable
    inverse_contract: str      _INVERTIBLE | _VALUE_RAISE | _STRUCTURAL_ABSENT
    Explicit raise (-O-safe — un-collected helper, vv Mode-8)."""
    from orpheus.numerics.operator import (
        NotInvertible, MissingAdjoint, invertible as _inv_guard, adjointable as _adj_guard)
    # (a) inverse axis — the structural/value split IS part of the contract
    if op.is_invertible != invertible:
        raise AssertionError(f"{op!r}: is_invertible={op.is_invertible} != {invertible}")
    if invertible:
        op.inverse()                                    # MUST return (I1 covers the value)
    elif inverse_contract == _VALUE_RAISE:
        if not hasattr(op, "inverse"):
            raise AssertionError(f"{op!r}: value-dep type must DECLARE inverse()")
        try:
            op.inverse(); raise AssertionError(f"{op!r}: singular inverse() did not raise")
        except NotInvertible:
            pass
    else:                                               # _STRUCTURAL_ABSENT
        if hasattr(op, "inverse"):
            raise AssertionError(f"{op!r}: structural non-invertible must NOT declare inverse()")
    # (b) adjoint axis — uniformly eager-return-or-raise (.H is on the base)
    if op.is_adjointable != adjointable:
        raise AssertionError(f"{op!r}: is_adjointable={op.is_adjointable} != {adjointable}")
    if adjointable:
        op.H                                            # eager: returns the wrapper, no raise
    else:
        try:
            op.H; raise AssertionError(f"{op!r}: non-adjointable .H did not eager-raise")
        except MissingAdjoint:
            pass
    # (c) bridge-consistency — pins the one-line TypeGuard bodies against drift
    if _inv_guard(op) != op.is_invertible:
        raise AssertionError(f"{op!r}: invertible() bridge drifted from is_invertible")
    if _adj_guard(op) != op.is_adjointable:
        raise AssertionError(f"{op!r}: adjointable() bridge drifted from is_adjointable")
```

**Leg (a) — the structural/value split enumeration (the contract's
teeth).** Every operator is tagged `_INVERTIBLE` / `_VALUE_RAISE` /
`_STRUCTURAL_ABSENT`; the tag is CONTRACT, not incidental. From the
verified tree:

| Contract | Operators |
|---|---|
| `_VALUE_RAISE` (declares `inverse()`, raises `NotInvertible` when singular) | `DiagonalOperator`, `MultiplicationOperator`, `OperatorSum` (non-invertible head), `OperatorProduct` (non-invertible factor), `ScaledOperator` (non-invertible operand), `TensorProductOperator` (non-invertible factor) |
| `_STRUCTURAL_ABSENT` (`not hasattr(op,"inverse")`) | `ZeroOperator`, `_ApplyOnly`, `IncomingOrdinateMaskTensor`, `PeriodicWrapOperator`, `RankOneOperator`, `StreamingOperator` (`L`), `IsotropicScattering`, `IsotropicN2N`, `ScatteringOperator`, `LegendreMomentScattering`, `N2NMomentOperator`, `FissionOperator`, `SNBoundaryOperator` (`B`), `ReconstructionOperator`, `_AdjointOperator` (`A.H` — `A.H.inverse()` #280-deferred) |
| `_INVERTIBLE` (declares `inverse()`, returns) | `IdentityOperator`, `PermutationOperator` (algebra-closed), `DiagonalOperator`/`Mult` (min\|f\|>0), invertible `OperatorSum`/`Product`/`Scaled`/`TP`, the inverse family (`InverseOperator`/`SweepOperator`/`GreenOperator`/`MatrixInverseOperator`), `_BoundBoundaryOperator` (reflective law) |

**Config-blindness (§0.6) — the two asymmetry fixtures are MANDATORY**
(else blind to the axis SEPARATION, exactly as the §2.3 scaffold
required): (i) `ZeroOperator` — `is_adjointable=True`, `is_invertible=False`,
`_STRUCTURAL_ABSENT`; (ii) `MultiplicationOperator`/`DiagonalOperator`
with a **TRUE zero** coefficient (not `1e-300`) — `is_adjointable=True`,
`is_invertible=False`, `_VALUE_RAISE`. And the half-adjointable composite
`matrix_full + matrix_apply_only` — `is_adjointable=False`, `is_invertible=True`
(leading term invertible). Without all three, a buggy `is_adjointable`
that returned `is_invertible` (or vice versa) passes silently.

**Leg (d) — the recursive-composition pins survive VERBATIM.** These
already live in the enumeration callers (`test_sum_predicates_recursive_and_faithful`,
`test_product_predicates_recursive_and_faithful`) and never reference the
frozenset — only the `_assert_faithful(op)` CALLS inside them strip out,
replaced by `assert_inverse_adjoint_contract(...)`:

```python
assert (a + b).is_adjointable == (a.is_adjointable and b.is_adjointable)
assert (a + b).is_invertible  == a.is_invertible          # leading-term rule (#261; step-4 landed)
assert (a @ b).is_invertible  == (a.is_invertible  and b.is_invertible)
assert (a @ b).is_adjointable == (a.is_adjointable and b.is_adjointable)
assert (2.0 * a).is_invertible  == a.is_invertible
assert (2.0 * a).is_adjointable == a.is_adjointable
# _AdjointOperator: A.H.is_invertible is False AND A.H.is_adjointable is False (#280 deferred)
```

**Teeth:** M-ADJ-PROP, M-INV-FORGE, M-BRIDGE (§41) all red THIS gate.

---

## 37. MIGRATION-MAP UPDATES — supersede §4a/§4e; the WIDER inventory the 2026-06-29 map missed

**Pressure-test (c) finding (I re-grepped `capabilities` across `tests/`
myself, as the brief demanded).** The §4a/§4e tables predate steps 4–5
and are INCOMPLETE: **127 reflective `.capabilities` reads across 33
files**, plus NEW step-2/4/5 files the map never saw. The migration is
NOT a fixed table — it is a MECHANICAL RULE + a completeness re-grep gate.

### 37.1 The mechanical rewrite rule (apply to every reflective read)

| Today | → Rewire |
|---|---|
| `CAP_SOLVE in op.capabilities` | `op.is_invertible` |
| `CAP_SOLVE not in op.capabilities` | `not op.is_invertible` |
| `CAP_APPLY_TRANSPOSE in op.capabilities` | `op.is_adjointable` |
| `CAP_APPLY_TRANSPOSE not in op.capabilities` | `not op.is_adjointable` |
| `CAP_APPLY in op.capabilities` | DROP (apply universal) or `callable(getattr(op,"apply",None))` |
| `op.capabilities == frozenset({CAP_APPLY})` (strict eq) | `not op.is_invertible and not op.is_adjointable` |
| `op.capabilities == frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})` | `op.is_adjointable and not op.is_invertible` |
| `op.capabilities == frozenset({CAP_APPLY, CAP_SOLVE, CAP_APPLY_TRANSPOSE})` | `op.is_invertible and op.is_adjointable` |
| **bare-string** `"solve"/"apply_transpose" in op.capabilities` (§1d.2) | `op.is_invertible` / `op.is_adjointable` |
| `pytest.raises(MissingCapability, match="invert\|solve\|singular")` | `pytest.raises(NotInvertible, match="invert")` |
| `pytest.raises(MissingCapability, match="apply_transpose\|adjoint")` | `pytest.raises(MissingAdjoint, match="adjoint")` |
| `pytest.raises(MissingCapability, match="apply")` (composition apply-guard) | `pytest.raises(TypeError, match="apply")` |

The **strict-equality** reads (which caught a *spurious added* cap —
`test_capability_survival.py:106-108` IncomingSource, `test_tensor_product_operator.py:99-123`)
map to the axis-predicate CONJUNCTION; the "spurious cap" concern
becomes "spurious `is_invertible`/`is_adjointable` True", still caught.

### 37.2 The COMPLETE file inventory (33 files — the map's §4a/§4e caught ~12)

**Newly-surfaced (steps 2/4/5 added them; the 2026-06-29 map predates
them) — MUST be in the migration:** `test_matrix_inverse_operator.py` (3),
`test_inverse_operator_equivalence.py` (1), `test_scattering_adjoint.py`
(11), `test_fission_adjoint.py` (2), `test_invertible_operator.py` (4),
`test_streaming_operator.py` (5), `test_periodic_wrap_operator.py` (1),
`test_incoming_ordinate_mask_tensor.py` (1), `test_collision_operator.py`
(1), `test_angular_average_operator.py` (1), `test_estimators_as_functionals.py`
(1). **Enumerated in §4a/§4e (verify anchors drifted):** `test_operator.py`
(17), `test_capability_survival.py` (11), `test_tensor_product_operator.py`
(11), `test_outer_dyad.py` (8), `test_multiplication_operator.py` (7),
`test_diagonal_operator.py` (5), `test_sn_boundary_operator.py` (4),
`test_legendre_moment_scattering.py` (4), `test_isotropic_scattering.py`
(4), `test_fission_operator.py` (4), `test_scattering_operator.py` (3),
`test_g_adjoint_reciprocity.py` (2, incl. the :210 precondition),
`test_frame_conjugate_carve.py` (2), `test_iteration.py` (2),
`test_frame.py` (2), `test_bc_universal_invariants.py` (3),
`test_bound_compat.py` (2), `test_removal_form_matvec_sweep.py` (1, the
:308 precondition), `test_boundary.py` (1), `test_prescribed_inflow_consistency.py`
(1), `test_permutation_operator.py` (1), and `tests/_harness/predicates.py`
(the scaffold — DELETES/rewrites, §36).

### 37.3 The COMPLETENESS GATE (this is the safety net for pressure-test c)

The migration is DONE only when **`grep -rn "\.capabilities\|MissingCapability\|CAP_SOLVE\|CAP_APPLY_TRANSPOSE\|CAP_APPLY\b\|_has(" tests/ orpheus/`
returns ZERO hits** (except the git-history / plan-file cross-refs). Run
it as the final W2 gate; a single surviving reflective read
`AttributeError`s the instant the frozenset deletes (the frozenset is
gone → `getattr(op,"capabilities")` raises). **Mutation M-LITERAL-STRING**
(§41): leave ONE bare-string precondition unrewired
(`"apply_transpose" not in op.capabilities`, the §1d.2 canary trap at
`test_g_adjoint_reciprocity:210` / `test_removal_form_matvec_sweep:308`)
→ that canary `AttributeError`s post-deletion (a false RED masking
whether S† held). The constant-grep MISSES these (bare string) — the
`.capabilities`-grep catches them; both greps are mandatory.

### 37.4 GAP-3′ gets the TWO-LEGGED shape (structural-absent vs value-raise)

§4b's GAP-3′ (`not hasattr(op,"inverse")` + SI raises) is now TWO-legged
— the structural/value split (§36 leg a) is the contract:

```python
def test_apply_only_operator_has_no_inverse(matrix_apply_only):
    """P6: an APPLY-ONLY / STRUCTURAL non-invertible has NO inverse() method."""
    assert not hasattr(matrix_apply_only, "inverse")          # structural-absent
    with pytest.raises(TypeError, match="apply"):             # SI apply-axis guard (design B)
        SourceIteration(matrix_apply_only)     # NOTE: SI now gates on apply, not invertibility

def test_singular_value_dep_operator_inverse_raises():
    """P6: a VALUE-DEP non-invertible DECLARES inverse(), which raises."""
    D = DiagonalOperator(np.array([2.0, 0.0, 0.5]))           # a TRUE zero (not 1e-300)
    assert hasattr(D, "inverse")
    with pytest.raises(NotInvertible, match="invert|singular"):
        D.inverse()
```

**⚠ Design B changes the SI gate semantics.** SI no longer gates on
invertibility (the operator arrives PRE-inverted — iteration.py:507-523
already gates CAP_APPLY only, step 3). So the apply-only rejection is an
APPLY-axis `TypeError(match="apply")`, NOT a `NotInvertible`. M4′ (§41)
neuters this apply-guard.

### 37.5 §4a/§4e ROWS SUPERSEDED by steps 4–5 (write the corrected dispositions)

| §4a/§4e row | §4-era shape (OBSOLETE) | Step-6 corrected shape |
|---|---|---|
| `test_sum_solve_does_not_propagate` (:132) | `not hasattr(A+B, "inverse")` — "sum not invertible" | **OBSOLETE** (step 4: sums DO invert via Green). New: `type((A+B).inverse()) is GreenOperator` (leading-term invertible) AND `not hasattr(A+B, "solve")` (Design B deletes `OperatorSum.solve`). The `is_invertible is False` premise RETIRED (§24.2 landed it). |
| §4c fork (sum Krylov-invertible?) | open fork | **RESOLVED** — sums invert via `GreenOperator` (step 4). No `KrylovInverseOperator`-on-sum question remains. |
| `test_scaled_solve_divides` (:301) | value via `.solve` | `ScaledOperator.solve` DELETED (Design B) → rewire to `(2*A).inverse().apply(b)` (= `(1/α)·A⁻¹b`, bit-id). |
| `test_product_solve_reverses_order` (:287/:293) | value via `.solve` | `OperatorProduct.solve` KEPT-rerouted → the value pin STAYS on `.solve` AND gains the §40 re-route equivalence gate + sentinel. |
| GAP-3′ (§4b) | one-legged `not hasattr` | TWO-legged (§37.4). |

### 37.6 The solve-DELETION blast radius (Design B — grep `\.solve(` in tests + prod)

Deleting `solve` on Identity/Permutation/Scaled/TP/OperatorSum/`_BoundBoundaryOperator`:

- **`OperatorSum.solve` (:1097)** — `InvertibleOperator` (streaming.py:844)
  and `_ScheduledInvertibleOperator` (:184) OVERRIDE it (their own direct
  sweep) → their MRO is unaffected; only GENERIC-sum callers break.
  Migrate `(L+C−B).solve(b)` → `.inverse().apply(b)` (= `GreenOperator(...).apply`).
  Step-4 tests exercising `Sum.solve` (`test_green_operator.py` back-half
  anchor, `test_capability_survival.py:154-163`) → rewire to `.inverse().apply`.
- **`PermutationOperator.solve` (:1753)** = `Pᵀb` (integer gather) —
  behavioral callers rewire to `.inverse().apply` (bit-id, §40); API-smoke
  (`CAP_SOLVE`-exists / `test_permutation_operator.py` caps read) → delete.
- **`ScaledOperator.solve` / `IdentityOperator.solve` / `TensorProductOperator.solve`**
  — behavioral → `.inverse().apply`; API-smoke → delete.
- **`_BoundBoundaryOperator.solve` (_bound_compat.py:161)** — the
  `test_bound_compat.py` solve-forward test → rewire to `.inverse().apply`
  (the shim's inverse forwards the inner law's inverse) OR delete if it
  was only a caps-existence smoke.

Per coding-standards retirement-means-test-migration: behavioral → rewire;
API-smoke → delete; state NONE-FOUND for any category with zero members.

---

## 38. THE EAGER-`.H` MIGRATION AUDIT (§3b executed) — the raise site moves `.apply` → `.H`

**The behavior change.** Today `A.H` on a NON-adjointable `A` SUCCEEDS
(constructs `_AdjointOperator` unconditionally, op.py:644) and the raise
is LAZY at `.apply` (op.py:911-916). Design C raises EAGERLY at `.H`
(`if not adjointable(self): raise MissingAdjoint` in `.adjoint()`), so
the op.py:911-916 lazy gate is dead and DELETES.

**Migration audit (grep tests/ for the lazy-raise pattern).** On HEAD the
only `.capabilities`-transpose reads in the reciprocity/removal gates are
`pytest.fail` PRECONDITIONS (`test_g_adjoint_reciprocity:210`,
`test_removal_form_matvec_sweep:308` — §1d.2/§37.3), NOT lazy-raise
assertions. Grep for BOTH shapes at carve time (the tree edits live):

1. `wrapper = A.H; with pytest.raises(...): wrapper.apply(x)` (lazy) →
   `with pytest.raises(MissingAdjoint): A.H` (eager). Enumerate + rewire.
2. `A.H` on a non-adjointable operand relied upon to SUCCEED constructing
   the wrapper (any test that builds `A.H` then never calls `.apply`, for
   a non-adjointable `A`) → now raises eagerly → those constructions
   break. Grep `\.H\b` in tests, filter for non-adjointable operands.

**Gate (the M-ADJ-EAGER target).** The eager-raise pin (NEW, home
`test_operator_capability_predicates.py` or the keystone leg b):

```python
def test_half_adjointable_sum_H_raises_eagerly():
    s = matrix_full + matrix_apply_only          # is_adjointable False
    with pytest.raises(MissingAdjoint, match="adjoint"):
        s.H                                        # EAGER — not lazy at .apply
```

**Mutation M-ADJ-EAGER (§41):** `.adjoint()` builds the wrapper WITHOUT
the `adjointable(self)` precheck → `(full+apply_only).H` returns a wrapper
instead of raising → this gate reds.

---

## 39. STATIC-LAYER GATES — pin the narrowing like step-5's kwarg tooth (CLI pyright is the ORACLE)

The frozenset retirement deletes **10 `# type: ignore`s** (op.py:934,
1040, 1176, 1180, 1284, 1287, 2040, 2051 + _bound_compat.py:110) and **4
`cast(SupportsInverse,…)`** (op.py:1310, 2076, iteration.py:231,
_bound_compat.py:143). The PROOF the deletions are principled (anti-#19,
not suppression) is CLI-pyright-checked. Three teeth:

### 39.1 Guard-delete pyright teeth (M-GUARD-DELETE-PYRIGHT)

Delete a guard-narrow → the narrowing is gone → CLI pyright RED. Name the
expected diagnostic per site:

| Delete | Expected `npx pyright orpheus/` diagnostic |
|---|---|
| the `if not invertible(self.op): raise` in `ScaledOperator.inverse` (op.py:1310 region) | `reportAttributeAccessIssue` on `self.op.inverse()` (base `LinearOperator` has no `inverse`) |
| the guard-narrow in `OperatorProduct.solve` re-route (§40) | `reportAttributeAccessIssue` on `self.a.inverse()` / `self.b.inverse()` |
| the `if not adjointable(...)` in a composer `apply_transpose` (e.g. `OperatorSum` op.py:1040) | `reportAttributeAccessIssue` on `self.a.apply_transpose(x)` |

### 39.2 The ratchet — re-baseline DOWN in-commit

Deleting 10 ignores + 4 casts REDUCES the numerics pyright count. The
ratchet (`tests/_harness/pyright_ratchet.py`, `pyright_baseline.json`,
current **148** = numerics 21 / sn 102 / transport 25) MUST be
re-baselined DOWN in the same commit (`python -m tests._harness.pyright_ratchet
--update` after the carve). A ratchet that STAYS at 148 means the carve
either added a new red (regression) or the ignores masked already-clean
lines — investigate. Expected: numerics drops by the count of net-cleared
reds; assert the EXACT new floor (lessons L10 — never trust the count,
assert it). Verify via CLI `npx --no-install pyright orpheus/`, NOT the
LSP (advisory noise, #226).

### 39.3 The static-pin file — extend `test_operators_apply_typed.py` (one home, lessons L10)

EXTEND `_inverse_family_seeded_apply_static_pins` (:447) / the C6 block
(:355) with the narrowing contracts, pyright-checked:

```python
# member-access narrowing type-checks with NO cast (the Design-C probe)
def _narrowing_static_pins(op: LinearOperator, b) -> None:
    if invertible(op):
        assert_type(op.inverse(), LinearOperator)   # narrowed → SupportsInverse.inverse
    if adjointable(op):
        assert_type(op.apply_transpose(b), object)  # narrowed → SupportsAdjoint.apply_transpose
```

**`Zero().inverse()` unspellable** — a SEPARATE pyright-checked file
(NOT runtime; the expression is a static error, so it cannot live in a
collected test that would `AttributeError` at import). Home: a
`# pyright: strict` fixture file with an expected-error comment, run
under `npx pyright`:

```python
# tests/sn/operators/_static_zero_inverse_unspellable.py  (pyright-only, NOT collected)
ZeroOperator().inverse()   # pyright: reportAttributeAccessIssue  (inverse NOT on base)
```

**Prove the teeth (CLI pyright).** Mutate one declared overload return
(make `SweepInverseOperator.apply` return `np.ndarray`) → `reportAssertTypeFailure`
fires under `npx pyright orpheus/ tests/sn/operators/`. M-MINV-KWARG /
M-MIXIN-KWARG precedent: the pyright tooth is verified on a `cp`-backup
mutation (never `git checkout` — the file holds uncommitted edits).

> **✅ §39.3 IMPLEMENTED (2026-07-02, W3) — as-built record + one delta:**
> 1. `_typeguard_bridge_narrowing_static_pins` appended to
>    `test_operators_apply_typed.py` — three legs per axis: member access
>    RESOLVES; the narrowed signature carries the honest carrier SWAP
>    (pinned on a `LinearOperator[ScalarFlux, AngularFlux]` two-space
>    operand — an endomorphism cannot see a swap bug); the narrowed
>    operand still IS a `LinearOperator` (the Protocol-promotion payoff,
>    §44.E). No `cast` anywhere.
> 2. **DELTA — the unspellable fixture is INVERTED-POLARITY, not a
>    standing red.** `tests/` sits in the CLI-checked tree (the pyproject
>    `ignore` list covers only orpheus subpackages), so the spec's bare
>    expected-error file would have parked a PERMANENT red in every
>    full-tree run — the very noise #226 kills. As-built:
>    `_static_zero_inverse_unspellable.py` puts the forbidden expression
>    inside a never-called function under a rule-scoped
>    `# pyright: ignore[reportAttributeAccessIssue]` with file-level
>    `reportUnnecessaryTypeIgnoreComment=error` — GREEN while the
>    unspellability holds (the ignore is necessary), RED the moment
>    `inverse()` appears on Zero/base (the ignore turns unnecessary).
>    The ignore comment is the ASSERTION mechanism, not suppression
>    (delete it → red today; property breaks → red then). ⚠ Brushes the
>    user's no-`# type: ignore` line — flagged for veto; swap to a
>    config-excluded must-red fixture is a 5-minute change if rejected.
> 3. **Teeth (both cp-backup on `operator.py`, CLI-verified):**
>    M-BRIDGE-ANNOT — widening `invertible()`'s return annotation to
>    `bool` reds the pins file 3→8 (the `SupportsInverse` assignment +
>    both `.inverse()` accesses un-narrow); M-STATIC-UNSPELLABLE —
>    declaring `inverse()` on `ZeroOperator` reds the fixture 0→1
>    (`reportUnnecessaryTypeIgnoreComment`). Both restored byte-identical.
> 4. The pins file carries 3 PRE-EXISTING CLI errors (`BC.reflective`
>    dynamic-enum reads at :82/:93 — the #254/#255 declare-dynamic-attrs
>    family, tests-tree debt outside the ratchet scope; untouched).
>    Ratchet after all W3 edits: EXACTLY 145 (21/99/25), exit 0.

---

## 40. THE `OperatorProduct.solve` RE-ROUTE — the one body rewrite (Mode-11 sentinel + bit-id baseline + dense anchor)

**The rewrite (Design B).** `OperatorProduct.solve` (op.py:1174) changes
`self.b.solve(self.a.solve(b_vec))` →
`if not invertible(self.a) or not invertible(self.b): raise NotInvertible(...);
return self.b.inverse().apply(self.a.inverse().apply(b_vec))`. This kills
the op.py:1176 `# type: ignore` (the guard narrows `a`/`b` to `SupportsInverse`)
AND FIXES the deleted-solve factors: a `Permutation`/`Scaled`/`Sum`/`TP`
factor (whose `.solve` Design B deleted) now works via its `.inverse().apply`.

**Why it EXECUTES (pressure-test a — I verified the call graph).**
`OperatorProduct.inverse()` (#285, op.py:1189) returns `InverseOperator(self)`;
`InverseOperator.apply` (op.py:1656) = `self.inner.solve` = `OperatorProduct.solve`
(re-routed). So `(A@B).inverse().apply(q)` DOES execute the re-route. But
a value gate comparing `.inverse().apply` to a reference SPELLED the same
way is tautological — the gate MUST (i) construct MIXED, NON-COMMUTING,
solve-deleted-kind factors, (ii) SENTINEL-wrap `OperatorProduct.solve` to
prove execution (Mode-11), (iii) anchor VALUE on a structurally-independent
dense reference, (iv) prove INHERITANCE vs the pre-carve baseline.

### 40.1 Gate REROUTE — the four legs

```python
def test_product_solve_reroute(monkeypatch):
    """§40: OperatorProduct.solve re-routes through .inverse().apply,
    bit-identically, and EXECUTES (Mode-11 sentinel)."""
    D = DiagonalOperator(_C7)                    # kept-solve leaf
    P = PermutationOperator(_P7_3cycle)          # solve-DELETED, algebra-closed, NON-involution
    prod = D @ P                                 # NON-COMMUTING (config-blindness §0.6)
    q = _RNG.standard_normal(7)
    # (ii) SENTINEL — prove OperatorProduct.solve is entered (else Mode-11 vacuous)
    calls = []
    real = OperatorProduct.solve
    monkeypatch.setattr(OperatorProduct, "solve",
                        lambda self, b: calls.append(1) or real(self, b))
    out = prod.inverse().apply(q)                # routes InverseOperator → OperatorProduct.solve
    assert calls, "OperatorProduct.solve NEVER executed — gate is Mode-11 vacuous"
    # (iii) INDEPENDENCE — dense closed-form (structurally independent of the re-route)
    dense = np.linalg.solve(prod.as_matrix(basis_shape=(7,)), q)
    np.testing.assert_allclose(out, dense, rtol=1e-12)
    # (iv) INHERITANCE — vs the pre-carve --capture-baseline snapshot (direct → array_equal)
    assert_regression("product_reroute_DP", out)     # bit-id (integer gather + division)
```

### 40.2 The factor-kind matrix + the honest tolerance split (pressure-test b resolved)

The re-route's bit-identity claim holds per factor kind — I traced each:

| Factor kind | old `factor.solve(x)` | new `factor.inverse().apply(x)` | tolerance |
|---|---|---|---|
| `DiagonalOperator` | `x/diag` | `InverseOperator(D).apply` = `D.solve` = `x/diag` | **array_equal** (same call) |
| `PermutationOperator` | `Pᵀx` gather (:1753) | inverse-perm `.apply` gather | **array_equal** (integer gather) |
| `ScaledOperator(α,op)` | `(1/α)·op.solve(x)` (:1284) | `ScaledOperator(1/α,op.inv).apply` = `(1/α)·op.solve(x)` | **array_equal** (identical float path) |
| `IdentityOperator` | `x` | `Identity.apply` = `x` | **array_equal** |
| `OperatorSum` (invertible head) | `GreenOperator(sum).apply(q)` (:1104) | `GreenOperator(sum).apply(q)` | **`assert_regression` (SAFETY×driver_tol)** — iterative (below) |

**The Green-factor row (pressure-test b — I decided it HONESTLY).**
Green-vs-Green IS mathematically bit-identical (both spellings construct
`GreenOperator(sum)` DETERMINISTICALLY — same left-spine split, same
preconditioner, same negated gains — and call `.apply(q, initial_guess=None)`,
a deterministic sweep+vector sequence). BUT it routes through the
ITERATIVE driver + the §18.A refinement loop, and lessons L7 warns
iterated results drift 1000s of ULP CROSS-RUN/CROSS-MACHINE (FP jitter in
the driver's internal reductions). So `array_equal` is the WRONG gate for
the sum-factor row — it would false-RED on a re-run. **Decision:** gate
the four DIRECT factor rows at `assert_array_equal` (bit-id, no reduction
change); gate the SUM-factor row at `assert_regression` (SAFETY×driver_tol,
lessons L7 iterative), layered under `-W error::DriftWarning` so a
genuine value change (not FP jitter) still escalates. Do NOT claim
`array_equal` on the iterative row — that is the L7 "iterated snapshot as
bit-id gate" trap.

**Config-blindness (§0.6):** the factors MUST be NON-COMMUTING
(`Diagonal @ Permutation`, and a 3-cycle permutation — non-involution)
for M-PROD-REROUTE (the factor-order swap) to red; a commuting `D1 @ D2`
is BLIND (the designed-green control, §41).

---

## 41. MUTATION MATRIX — Step 6 (post-retirement forms; each REDs a NAMED gate, under `-O`, monkeypatch only)

In-process monkeypatch ONLY (uncommitted edits — NEVER `git checkout`/
`restore`/`stash`, process-discipline.md; the tree is LIVE-edited).
Every mutation names its ACTIVATING config (§0.6).

| ID | Target gate | Mutation | Expected RED | Activating config |
|---|---|---|---|---|
| **M-ADJ-PROP** | §36 leg d + eager-`.H` (§38) | `OperatorSum.is_adjointable = a.is_adjointable` (drop `and b.is_adjointable`) | `(full+apply_only).is_adjointable` wrongly True → §36 leg d reds; `(full+apply_only).H` stops eager-raising → §38 reds | **half-adjointable** `full+apply_only` (both-adjointable sum is BLIND) |
| **M-ADJ-FORGE** | downstream `.H.apply` (reciprocity canary §42) | monkeypatch the half-adjointable sum's class `is_adjointable → True` | the eager `.H` now BUILDS a wrapper; `(full+apply_only).H.apply(v)` **AttributeErrors** (the apply_only summand has no `apply_transpose`) → the `.H`-exercising gate reds (a RAISE, not a wrong value — brief's precision note) | half-adjointable sum |
| **M-INV-FORGE** | §36 leg a (structural/absence) | `ZeroOperator.is_invertible → True` (force) | §36's `if invertible: op.inverse()` branch → `Zero.inverse()` **AttributeError** (Zero has no `inverse` — structural-absent) → leg a reds | `ZeroOperator` (the asymmetry fixture) |
| **M-ADJ-EAGER** | §38 eager-raise pin | `.adjoint()` builds the wrapper WITHOUT the `adjointable(self)` precheck | `(full+apply_only).H` returns a wrapper instead of raising `MissingAdjoint` → §38 reds | half-adjointable sum |
| **M-PROD-REROUTE** | §40 REROUTE (value + sentinel) | swap the re-route order: `self.a.inverse().apply(self.b.inverse().apply(b_vec))` | dense-anchor + baseline reds (value `≠ P⁻¹(D⁻¹q)`); sentinel STILL fires (execution unaffected — proves the sentinel is not the value catch) | **NON-COMMUTING** `Diagonal @ 3-cycle-Permutation` |
| **M-BRIDGE** | §36 leg c (bridge-consistency) | `invertible(op)` body → `return True` | leg c reds on any non-invertible op (`invertible(Zero)`=True ≠ `Zero.is_invertible`=False) | `ZeroOperator` |
| **M-LITERAL-STRING** | §37.3 + §42 canaries | leave `test_g_adjoint_reciprocity:210` precondition as `"apply_transpose" not in A.capabilities` (unrewired) | the canary **AttributeError**s post-deletion (no `.capabilities`) — a false RED masking whether S† held (the §1d.2 trap) | frozenset deleted |
| **M4′** | §37.4 GAP-3′ SI guard | neuter the SI apply-guard (`callable(getattr(op,"apply",None))` → always True) | the apply-only rejection (`SourceIteration(matrix_apply_only)`) no longer raises → GAP-3′ reds | apply-less L |
| **M-GUARD-DELETE-PYRIGHT** | §39.1 static teeth | delete a guard-narrow (Scaled.inverse / Product.solve re-route / a composer apply_transpose) | CLI `npx pyright` reds (`reportAttributeAccessIssue` — base has no `inverse`/`apply_transpose`) | any (CLI pyright, cp-backup) |

**Designed-GREEN control (mark explicitly):** **M-PROD-COMMUTING** — run
M-PROD-REROUTE's factor-order swap on a COMMUTING `D1 @ D2` → the value
gate does NOT red (commuting factors make the swap value-invisible). This
proves the §0.6 config-blindness discipline: the re-route value gate is
BLIND without non-commuting factors, so the non-commuting fixture is
load-bearing, not decoration.

> **✅ W3 BANK EXECUTED (2026-07-02) — every mutation RED-verified under
> `-O`.** Teeth at job-tmp `step6teeth/` (in-process plugin monkeypatch
> for the five class/bridge mutations; cp-backup byte-copy restore for
> the three file mutations — never git; every restore verified
> byte-identical + clean `git status`). Verdicts:
>
> | ID | Verdict |
> |---|---|
> | M-ADJ-EAGER | RED @ the W1 boundary (`m_adj_eager.py`) |
> | M-ADJ-PROP | RED — 4 named gates (keystone `[Sum-half-adjointable]`, the §38 eager pin, `test_sum_predicates_recursive_and_faithful`, `test_sum_transpose_drops_when_one_lacks`); both-adjointable rows stayed green (blind, as designed) |
> | M-ADJ-FORGE | RED — 3 gates (keystone `[Sum-leading-not]`+`[Sum-half-adjointable]`, eager pin) + inline precision demo: the forged wrapper's `.apply` raised **`MissingAdjoint`** (the as-built guard-narrow upgraded the spec's predicted raw `AttributeError` — still a RAISE, never a wrong value) |
> | M-INV-FORGE | RED — keystone `[Zero]` (leg-a predicate-mismatch fires before the `inverse()` branch; the structural-absent `hasattr` leg backstops) |
> | M-BRIDGE | RED — EXACTLY the 10 `is_invertible=False` rows at leg (c), all invertible rows green (perfect selectivity) |
> | M-PROD-REROUTE | RED — both §40 gates; the sentinel STILL fired (execution proof ≠ value catch); dense anchor + baseline were the catchers |
> | M-PROD-COMMUTING | **GREEN** (designed control) — the swapped body on `D1 @ D2` passes the rtol=1e-12 dense anchor; claim scoped to the VALUE gate (the bit-id baseline row may ULP-drift — not the control's leg) |
> | M-LITERAL-STRING | RED (false-red demo) — all 5 canary cases `AttributeError: 'OperatorSum' object has no attribute 'capabilities'` IN THE PRECONDITION, before any reciprocity math (the §1d.2 trap made concrete) |
> | M4′ | RED — `test_source_iteration_requires_apply_on_A_inv` fails "DID NOT RAISE TypeError" under the neutered guard |
> | M-GUARD-DELETE-PYRIGHT | RED — `operator.py` CLI errors 0 → 1 (`reportAttributeAccessIssue`: "Attribute inverse is unknown" on `self.op.inverse()` in `ScaledOperator.inverse`) → restored 0 |

**Mode-8 hygiene:** every gate uses `np.testing.*` / `pytest.raises` /
`assert_type` / collected-module `assert` (pytest-rewritten → fires under
`-O`). The eager `NotInvertible`/`MissingAdjoint` raises, the composition
`TypeError` apply-guards, and the value-dep singular guards all `raise`
(not `assert`) so they survive `-O` in production. The §36 helper
`raise`s (un-collected module, §36 note).

---

## 42. S†-STAYS-GREEN + the CANARY SET (§4e executed) — the no-regression floor

The reciprocity LAW is the METHOD BODY (`a.H + b.H`, op.py:1040;
`loss_action_transpose`), UNTOUCHED by the retirement — only the
advertisement moves. So the S†-exercising gates stay green; the
membership-asserting gates rewire (§37.1). **"S† stays green" is NOT
automatic** — it REQUIRES the two literal-string PRECONDITION rewires:

| Canary | Precondition rewire (REQUIRED — else AttributeError post-deletion) |
|---|---|
| `test_g_adjoint_reciprocity_full_block` (:200) | `:210` `if "apply_transpose" not in A.capabilities: pytest.fail` → `if not A.is_adjointable:` |
| `test_removal_form_matvec_sweep` S† twin (:290) | `:308` `if "apply_transpose" not in op.capabilities:` → `if not op.is_adjointable:` |

**The floor (stays green / stays RED, unchanged tolerance):**

| Canary | Disposition |
|---|---|
| `test_invertible_apply_transpose_…_bit_identical` (:290, slab/sphere/cyl) | STAYS `assert_array_equal` (S† bit-id) + `:308` rewired |
| `test_g_adjoint_reciprocity_full_block` (:200) | STAYS `rtol≤1e-10` + `:210` rewired |
| `test_wrong_trace_metric_breaks_reciprocity` (:268, L11 control) | STAYS UNCHANGED — MUST STILL RED on the wrong metric (the negative control that proves the reciprocity gate has teeth; §0.6 the metric-vs-Euclidean discriminator activates on sphere/cyl) |
| `test_iso_transpose_matches_dense_per_material` (:133), `test_euclidean_reciprocity` (:120) | unchanged (S† dense crosscheck) |
| keff L1/L2 net (`test_keff_curvilinear.py::test_homogeneous_exact` ≥2G, `test_2g_heterogeneous_fuel_moderator`) | tolerances UNCHANGED (a re-typing/retirement cannot move keff) |
| MMS net (`test_mms_ld_slab.py`, `test_mms_ld_2d.py` — not-slow slices) | flux-shape/order, unchanged |
| the ~117 landed step-4/5 family gates (`test_green_operator*.py`, `test_matrix_inverse_operator.py`, `test_inverse_universal.py`, `test_inverse_operator_equivalence.py`) | the no-regression floor — stay green UNCHANGED |
| the 0-ULP scattering S† canary (`test_isotropic_scattering.py`) | `assert_array_equal` unchanged |

---

## 43. LAND-ORDER CHECKLIST (§9-style) — W1 additive → W2 atomic retirement → W3 teeth

> **✅ W1+W2 IMPLEMENTED @ `f4919b1` (2026-07-02) — as-built record.**
> Both waves landed in one commit (63 files, +1765/−1425); full tier `-O`
> serial **3853/0**; §37.3 completeness re-grep ZERO; ratchet
> re-baselined **148 → 145** (numerics 21 / sn 99 / transport 25).
> Deltas vs this Part's letter:
> 1. **W1 parentage bridge:** `NotInvertible` was born subclassing
>    `MissingCapability` and the 8 `.inverse()`-guard raises swapped in
>    W1 (not W2) — resolving §43's internal tension (keystone leg (a)
>    green at the W1 boundary vs the 26 `pytest.raises(MissingCapability)`
>    lines). W2 re-parented to `TypeError`. Final state = the approved
>    design.
> 2. **Two production bugs found and fixed root-cause:**
>    (a) `_seeded_inverse` crashed on ALGEBRA-CLOSED preconditioner heads
>    (pre-existing since step 4; the §40 sum-row fixture exposed it) —
>    fixed with the `_wrap_delegate_member` TypeGuard + the
>    accept-and-drop `_SeededExactApply` adapter (iteration.py), pinned
>    end-to-end by the baseline's `sum_green` row; (b) `RankOneOperator`
>    lost its ONLY adjointability advertisement with the per-instance
>    caps deletion (**F† dead** behind the TP factor guard) — both
>    migration agents converged independently, pinned as strict-xfails;
>    the `is_adjointable` override landed, pins flipped, markers deleted.
> 3. **`assert_regression` did not exist** — §40 leg (iv) realized as a
>    PRE-carve `.npz` snapshot (`tests/numerics/data/
>    step6_product_reroute_baseline.npz`, captured from the OLD
>    `b.solve(a.solve(q))` spelling before the rewrite; seed 2860);
>    direct rows `assert_array_equal`, Green row `assert_allclose`
>    rtol=1e-8. Gates live in `tests/numerics/test_operator.py`.
> 4. **Two files sat OUTSIDE the §37.2 inventory** and were caught by the
>    tree-wide collection probe + the tier: `test_frame_conjugate_carve.py`
>    (strict-eq caps read) and `test_integral_kernel_category.py`
>    (surface-tuple `hasattr(op, "capabilities")` — a reflective read the
>    CAP-grep could not see). Both migrated.
> 5. The eager-`.H` §38 audit found ZERO lazy-raise patterns in tests —
>    the migration bucket was instead 4 W1 fixture-honesty completions
>    (`MatrixOperator`/`_SpacedMatrixOperator` predicates + a guarded
>    `inverse()`) and 2 stale pre-step-4 "sum never solves" pins rewired.
> 6. A transient +3 numerics pyright regression (swapping the dyad
>    transpose's isinstance for a predicate read LOST the narrowing) was
>    caught by §39.2's exact-floor discipline — the isinstance IS the
>    narrowing; kept in the body, the predicate delegates to the same fact.
> W3 (below) remains: the §41 bank (M-ADJ-EAGER already RED-verified at
> the W1 boundary), §39.3 static pins + the Zero-unspellable file, docs,
> enforcer.
>
> **✅ W3 COMMITTED @ `cb62310` (2026-07-02) — step 6 COMPLETE.** The
> full §41 bank RED-verified (verdict table at the §41 stamp), §39.3
> pins + inverted-polarity fixture + both CLI teeth (§39.3 stamp),
> archivist three-layer docs pass (+ the corrected pre-existing false
> S/F-adjointability claim, empirically confirmed), enforcer PASS with
> the 2 getattr-flattening polish items applied, 2 main-agent
> cross-check doc catches. Tier 3853/0; ratchet EXACTLY 145; Sphinx
> `-E -W` exit 0. The step-3 inverter-narrative docs staleness is
> tracked on #226 for the P5 docs slice.
>
> **✅ CARVE P5 (step-6 tail) COMPLETE — COMMITTED @ `70da74f`
> (2026-07-02; 3 files, +263/−58; task #135).** The §9 P5 charter executed:
> (1) composition-algebra return types — `ScaledOperator.inverse()` →
> `"ScaledOperator[Codomain, Domain]"` (the ONE parameter-dropping
> annotation among the nine surfaces) + the pyright-only
> `_composition_algebra_return_type_static_pins` bank in
> `test_operators_apply_typed.py` (M-SCALED-BARE teeth red-verified via
> `reportAssertTypeFailure`; enforcer SHIP zero-defects, twin-path
> acquitted, 1 NIT folded); the sum-inverse fork 4c needed NO action
> (resolved at steps 4/5: Sum→Green, `(L+C)`-subclass MRO-shadows→Sweep).
> (2) The §5 adjoint-coherence gates are NOT landed as xfail — an
> AttributeError-xfail verifies nothing; gates 5.1–5.3 + M-ADJ-swap/
> M-ADJ-metric + the double-swap static pin are INLINED self-contained
> in the #280 redesign comment (4871342146), which re-grounds #280 onto
> `A.H.inverse()` = the swap `(A.H)⁻¹ = (A⁻¹).H` +
> `SweepOperator.apply_transpose` (Euclidean reverse-scan; metric rides
> `_AdjointOperator`). (3) Docs: the inverter section rewritten →
> "Choosing the A⁻¹ realisation (the inverse-operator family)" (anchor
> `choosing-inverse-realisation`), tracking note deleted, Dev-history P5
> row added; Sphinx `-E -W` exit 0; #226 tracker comment RESOLVED
> (4871409512). (4) Ratchet unchanged at EXACTLY 145 — nothing to
> re-baseline (the carve cleared no reds; the tightening was
> ratchet-neutral by design). Full record: taxonomy plan §12 "6-tail".

### W1 — ADDITIVE (coexists with the frozenset; no reflective read breaks yet)

- Mint `NotInvertible`/`MissingAdjoint` (TypeError subclasses); add the
  free `invertible`/`adjointable` TypeGuards; promote
  `SupportsInverse`/`SupportsAdjoint` to extend `LinearOperator[D,C]`.
- Flip `.adjoint()`/`.H` EAGER (raising the new `MissingAdjoint`) +
  migrate ONLY the `.H` raise-site tests (§38) — the value-dep `.inverse()`
  guards STILL raise `MissingCapability` in W1 (unchanged; the full
  MissingCapability→NotInvertible swap is W2, else it breaks the 48
  `pytest.raises(MissingCapability)` lines mid-flight).
- **Gates at the W1 boundary:** §36 legs a/b/c run GREEN (they read the
  NEW predicates/guards, not the frozenset); the §2.3 scaffold STILL green
  (frozenset present); M-ADJ-EAGER reds; the `.H` migration audit clean;
  existing suites GREEN (additive); pyright ADDS no red (Protocol promotion
  + TypeGuards are accepted types — verify, ratchet still 148).

### W2 — THE RETIREMENT (atomic — the frozenset deletion breaks every reflective read at once)

- Delete the `capabilities` frozenset (all leaves + composers + the base
  field), `CAP_*` tags, `_has` (op.py:836 + iteration.py:204),
  `MissingCapability`, the §2.3 scaffold (`predicates.py` rewritten to §36).
- Swap value-dep `.inverse()` guards → `NotInvertible`; composition
  apply-guards → plain `TypeError` (SAME messages); adjoint guards →
  `MissingAdjoint`. Re-route `OperatorProduct.solve` (§40). Delete the 5
  algebra-type `.solve`s + `_BoundBoundaryOperator.solve` (§37.6). Delete
  the 10 `# type: ignore`s + 4 casts (§39).
- **FULL test migration** (§37): all 127 reflective reads (33 files, incl.
  the newly-surfaced step-2/4/5 files), all 48 `MissingCapability` lines by
  axis, the two literal-string canary preconditions (:210/:308), the §36
  keystone v2, GAP-3′ two-legged.
- **Gates at the W2 boundary:** §36 is the SOLE faithfulness gate (green);
  §40 re-route (baseline + dense anchor + sentinel) green; every §42
  canary green (incl. the two rewired preconditions + the L11 wrong-metric
  control STILL red); the **COMPLETENESS re-grep (§37.3) returns ZERO**
  hits; M-LITERAL-STRING confirmed (leaving one unrewired → AttributeError).

### W3 — TEETH + TIER + RATCHET + DOCS

- Run the §41 mutation bank: every mutation reds its named gate under
  `-O` (in-process monkeypatch); the M-PROD-COMMUTING control stays green;
  M-GUARD-DELETE-PYRIGHT via `cp`-backup CLI pyright.
- Extend the static-pin file (§39.3) + the `Zero().inverse()`-unspellable
  pyright file; prove `reportAssertTypeFailure` / `reportAttributeAccessIssue`
  teeth under `npx pyright`.
- **Re-baseline the ratchet DOWN** in-commit (§39.2) — assert the EXACT
  new floor (the 10 ignores + 4 casts clear; numerics count drops).
- Sphinx (Cardinal Rule 3): the both-axes retirement (`is_invertible`/
  `is_adjointable` from method bodies, NO twin-path — §0′/§0.6) into
  `discrete_ordinates.rst` "Development history" + the operator-algebra
  theory page. Sphinx `-W` EXIT 0.

**Done when:** for EVERY §36/§37/§38/§39/§40 gate a named §41 mutation
reds it under `-O`; §36 keystone v2 (the two-axis structural/value
contract) is the SOLE, PERMANENT faithfulness gate, referencing NO
frozenset; the completeness re-grep returns ZERO; the entire
`capabilities` frozenset + `MissingCapability` + `_has` RETIRED (both
axes); the 10 `# type: ignore`s + 4 casts DELETED; the ratchet
re-baselined DOWN and the EXACT new floor asserted; every §42 canary
bit-identical (incl. the two rewired literal-string preconditions and the
L11 wrong-metric control still red); Sphinx clean.

---

## 44. PRESSURE-TEST FINDINGS + HONEST-SCOPE + DESIGN OBJECTIONS

Surfaced BEFORE writing gates (the brief invited them; my PART IV
pressure-tests caught 3 real issues).

### 44.A (Mode-11, RESOLVED with a sentinel) — the re-routed `Product.solve` IS on the gate's call graph, but a value gate spelled the same way is tautological

`(A@B).inverse().apply(q)` → `InverseOperator(A@B).apply` = `(A@B).solve`
(re-routed) — VERIFIED the call graph (op.py:1189→InverseOperator,
:1656→inner.solve). So the re-route IS executed. But the value reference
MUST be structurally independent (dense `np.linalg.solve(as_matrix, q)`),
NOT `b.inverse().apply(a.inverse().apply(q))` (which re-spells the
production). §40 pins execution with a `OperatorProduct.solve` sentinel
(counter>0) AND anchors value on the dense reference AND inherits from the
pre-carve baseline. Without the sentinel the gate is Mode-11-vacuous if a
future refactor routes `.inverse().apply` around `.solve`.

### 44.B (pressure-test b, DECIDED HONESTLY) — Green-vs-Green is math-bit-identical but must NOT be gated `array_equal`

The sum-factor product row (`(Sum_invertible @ D).solve`) routes through
`GreenOperator(sum).apply` on BOTH the old and new spelling — mathematically
bit-identical (deterministic split + deterministic seed=None apply). BUT
it is ITERATIVE, and lessons L7 warns iterated results drift 1000s of ULP
cross-run/cross-machine. **Gate it at `assert_regression` (SAFETY×driver_tol),
NOT `array_equal`** — layered under `-W error::DriftWarning` so a genuine
value change still escalates. The four DIRECT factor rows (Diagonal/Perm/
Scaled/Identity) ARE `array_equal` (integer gather / identical float path,
no reduction change). Do NOT over-claim `array_equal` on the iterative row.

### 44.C (pressure-test c, RESOLVED) — the test surface is WIDER than the §4a/§4e map (I re-grepped)

**127 reflective `.capabilities` reads across 33 files** — the 2026-06-29
map enumerated ~12; steps 2/4/5 added ELEVEN files it never saw
(`test_matrix_inverse_operator`, `test_inverse_operator_equivalence`,
`test_scattering_adjoint`, `test_fission_adjoint`, `test_invertible_operator`,
`test_streaming_operator`, `test_periodic_wrap_operator`,
`test_incoming_ordinate_mask_tensor`, `test_collision_operator`,
`test_angular_average_operator`, `test_estimators_as_functionals`). §37
gives the MECHANICAL rule + the full inventory + the **completeness re-grep
gate** (§37.3) — the migration is done only when the grep returns ZERO.
This is the safety net a fixed table cannot provide.

### 44.D (HONEST-SCOPE — `_InvertibleForward.solve` docstring is now FALSE; UPDATE, do not execute)

`_InvertibleForward` (op.py:1463-1480, `InverseOperator`'s narrowing
bound) declares `solve` (:1480) with a docstring: *"the `solve` spelling
retires at P4, shrinking this contract with it."* That line is now WRONG:
`InverseOperator` wraps ONLY `DiagonalOperator` / `MultiplicationOperator`
/ `OperatorProduct` — and Design B KEEPS `solve` on all three (Product
re-routed). So `_InvertibleForward.solve` does NOT shrink; the contract
STAYS. **Fix the docstring** (state the leaf `solve` is KEPT as the
InverseOperator un-invert face), do not execute the shrink. Recorded as a
deliberate design finding, not a cut.

### 44.E (HONEST-SCOPE — TypeGuard REPLACES not INTERSECTS; guards belong at LinearOperator-typed sites)

The `invertible`/`adjointable` TypeGuards narrow a `LinearOperator` to
`SupportsInverse`/`SupportsAdjoint`. They belong at sites where the static
type is `LinearOperator` (the composer bodies, consumer pre-query). They
do NOT compose with an already-narrowed type (a site already typed
`SupportsInverse` needs no guard). Note it so the implementer does not
sprinkle redundant guards.

### 44.F (FLAG, do NOT decide — the repr caps-read successor, op.py:787)

`__repr__` reads `caps = sorted(getattr(self, "capabilities", frozenset()))`
(op.py:787). Post-retirement `getattr(..., frozenset())` returns the empty
set → EVERY operator reprs `caps=[]` (a degraded, misleading repr). The
successor is a design choice (show `is_invertible`/`is_adjointable`
flags, or drop the field) — FLAGGED for the user's steer, NOT decided
here. A gate is cheap (`repr(op)` contains the axis flags) once the
choice is made; do not gate the degraded interim.

### 44.G (verify-at-carve — the cast inventory + the SI-gate semantics)

Two items to RE-verify against the live tree at carve time (they drift):
(1) the 4 `cast("SupportsInverse", …)` sites (op.py:1310, 2076,
iteration.py:231 — QUOTED string, a bare-`cast(SupportsInverse` grep
MISSES it — _bound_compat.py:143); (2) the SI-apply-guard already gates
CAP_APPLY-only (step 3) — GAP-3′'s apply-only rejection is a
`TypeError(match="apply")`, NOT a `NotInvertible` (§37.4). A brief that
reads "SI requires invertible L" is STALE (step 3 pre-inverted the L).

---

## 11. References
- **Phase 0/1 spec (LANDED net this extends):**
  `.claude/plans/issue_226_b4_operator_generics_verification.md`.
- **The carve:** `.claude/plans/operator_inverse_algebra_carve.md` (§3 target, §5
  phases, §6 forks).
- **Map:** `.claude/agent-memory/explorer/issue_226_operator_generics_map.md` (§8
  B4 cluster; line anchors DRIFTED — re-derive via Nexus).
- **Landed #276 adjoint substrate (gate #4 / §1b):** `loss_action_transpose`
  (`sn/loss_representation/__init__.py`), `A.H` (op.py:543–702, the metric
  adjoint), `test_g_adjoint_reciprocity.py`, `test_isotropic_scattering.py` (S†),
  `test_removal_form_matvec_sweep.py` (the apply/apply_transpose bit-id twins).
- **Skills/rules:** `vv-principles` (anti-#2 MMS≠eigenvalue, anti-#19 no type:ignore,
  Mode-8 `-O`-strip, Mode-11 gate-never-executes, bit-identity vs principled),
  `coding-elegance` (Pattern 4 illegal-states, anti-#3 stringly-typed dispatch →
  the CAP_SOLVE retirement rationale), `.claude/rules/coding-standards.md`
  (retirement-means-test-migration), `.claude/rules/process-discipline.md`
  (monkeypatch-not-checkout for mutation under uncommitted edits).
- **test-architect lessons:** L2 (two-anchor: bit-id inheritance + closed-form),
  L4 (gate teeth / xfail-strict-False / call-count), L10 (refute-the-premise;
  assert_type static gate + CLI-pyright teeth), L12 (`.H` metric ≠ Euclidean
  transpose; CAP_APPLY_TRANSPOSE propagation; Mode-11 leaf-wrap).
