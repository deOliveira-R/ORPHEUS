---
name: issue-257-s5-functional-category-closeout
description: #257 S5 — ADD the §5.6 Functional(Protocol[V_contra,R_co]) category (numerics L1) + the concrete ProductionRateFunctional (transport L2, the F=M_χ∘ProductionRate∘M_νΣf middle); BIT-IDENTICAL additive stage, RankOne NOT rewired (S6); 0 net-new pyright, 0 type:ignore
metadata:
  type: project
---

# #257 S5 — Functional category (§5.6 suffix law) closeout

Branch `feature/field-typed-operator-algebra`, built on S4.5 HEAD `b404ae1`
(NOT committed — main agent commits after elegance + qa review). PURELY
ADDITIVE / BIT-IDENTICAL stage: ADD two new types + cross-check; NO rewire of
the fission operator (S6), NO touch to SN/CP eigenvalue production paths.

## Deliverables (the 7-item manifest)

1. **Branch-1 / numerics floor** — `orpheus/numerics/functional.py` (NEW, L1):
   `@runtime_checkable Functional(Protocol[V_contra, R_co])` with ONE method
   `evaluate(self, x: V_contra, /) -> R_co`. The §5.6 suffix law (Operator=
   field→field, Kernel=integrated-against-measure, **Functional=field→scalar**).
   Co-vector companion of `Vector`, SIBLING of `LinearOperator` (NOT subclass —
   deliberately ZERO operator surface: no apply/capabilities/solve/
   apply_transpose/H/domain/codomain). `__all__=["Functional","V_contra","R_co"]`.
   Exported from `numerics/__init__.py` alongside `Vector`/`V`.
2. **(no separate Branch-1 SymPy module)** — S5 mints NO new closed-form/
   semi-analytical reference. The §5.6 functional is a TYPE-SURFACE carve +
   a numerics-primitive delegation (the group contraction already exists as
   `RankOneOperator.inner`); the structurally-independent correctness ref is
   the test-side hand-derived double-loop (`hand_derived_production_density`,
   explicit Python accumulation, NO numpy reduction — L11 honoured). Same
   posture as S3a/S3b/#247/#251: NO algebra-of-record Branch-1/L1-derivation
   manifest owed.
3. **Branch-2 / transport leaf** — `orpheus/transport/production_rate_functional.py`
   (NEW, L2 — NOT sn/): `@dataclass(frozen=True) ProductionRateFunctional` with
   `nu_sigma_f: CrossSectionField` (the kw constructor the helper assumes) +
   `evaluate(phi: Field | NDArray) -> np.ndarray` =
   `(self.nu_sigma_f.values * phi_values).sum(axis=0, keepdims=True)`. Placed in
   transport/ because it imports ONLY transport types (CrossSectionField + the
   flux carrier) → the shared SN/CP/MoC leaf. The §5.6 frame: this is the MIDDLE
   of `F = M_χ ∘ ProductionRateFunctional ∘ M_νΣf` (Frame 3); S6 composes it.
   BIT-IDENTICAL to `RankOneOperator.apply`'s `inner` (right=νΣf, axis=0,
   keepdims). NO volume measure (density, not integral), NO χ re-broadcast.
4. **Foundation test gate (category)** — `tests/transport/test_functional_category.py`
   (test-architect-authored, PRE-LANDED). My code turns its 9 PRE-IMPL skips
   into live passes.
5. **L1 cross-check (correctness + equivalence)** —
   `tests/transport/test_production_rate_functional.py` (test-architect, PRE-
   LANDED). B.1 = hand-loop correctness (struct-indep ref, the L1 cross-check
   by pillar = the explicit double-loop, NOT the RankOne it'll replace), B.2 =
   RankOne `inner` bit-identity (DEMARCATED equivalence, de-risks S6), B.3 =
   no-measure Mode-3 guard. DEVIATION from the literal
   `tests/derivations/test_<name>_xverif*.py` folder convention: the
   verification spec deliberately lives in `tests/transport/` (the spec note
   is explicit; the cross-check ref is a test-side hand-loop, not a
   derivations/ Branch-1 module), so the xverif IS B.1 in this file. Justified.
6. **Estimators-as-functionals gate** — `tests/numerics/test_estimators_as_functionals.py`
   (test-architect, PRE-LANDED). C.1 = estimator bit-identity (PASSED today,
   unchanged arithmetic) + C.2 = category-honesty. The 1 remaining SKIP is the
   OPTIONAL `ProductionFunctional` bound wrapper leg — NOT shipped, per the
   test-architect recommendation (the lighter-touch design keeps the bare
   `(L,S,F,ψ)` estimator callables; they are NOT `Functional[V]` because they
   consume the operator triple, not a lone field). The category simply now
   NAMES what their field→scalar core is.
7. **Sphinx stub** — `docs/theory/operator_algebra.rst` new `.. _functional-category:`
   section: `:label: production-rate-functional` (Eq. 8) + a `functional-category`
   anchor + `:eq:`/`:mod:`/`:class:`/`:meth:`/`:func:` cross-refs + 2
   `.. vv-status: ... documented` + a `.. todo:: Archivist expansion needed`.
   DID NOT write the rich narrative (archivist's deliverable).
8. **Test-helper probe** (Deliverable 3) — `tests/transport/_functional_helpers.py`
   `require_production_rate_functional` now probes
   `"orpheus.transport.production_rate_functional"` FIRST (the transport home),
   sn/ candidates kept as fallback. THE ONLY test-file edit (the helper says
   "change this ONE function").
9. **Closeout memo** — this file.
10. **archivist DISPATCH_REQUEST** emitted (followup:false).

## THE 2 PYRIGHT HINGES (structural, 0 `# type: ignore`)

- **(1) Protocol typevar variance.** The literal brief signature
  `class Functional(Protocol[V, R])` reusing the IMPORTED INVARIANT
  `V` (from vector.py) emits 2 `reportInvalidTypeVarUse` WARNINGS
  ("V should be contravariant", "R should be covariant") — because a
  functional CONSUMES V (input-only) and PRODUCES R (output-only), unlike the
  operator endomorphism `apply(x:V)->V` which uses V BOTH ways (→ invariant,
  no warning). FIX (structural, NOT suppression): declare
  `V_contra = TypeVar("V_contra", bound=Vector, contravariant=True)` +
  `R_co = TypeVar("R_co", covariant=True)`. This is the CORRECT variance for a
  functional and the co-vector mirror of vector.py's invariant `V`. Renamed the
  public typevars accordingly (brief's `R` → `R_co`; added `V_contra`).
- **(2) `.values` on `Vector`.** The brief's reference body
  `phi.values if hasattr(phi, "values") else np.asarray(phi)` trips
  `reportAttributeAccessIssue` ("Cannot access attribute values for class
  Vector") — `Vector` Protocol has NO `.values`, and pyright does NOT narrow it
  via `hasattr`. FIX (structural): annotate `phi: Field | NDArray` (Field is L1
  numerics, has `.values: NDArray`; the flux carriers are Field subclasses,
  ScalarFlux IS-A Field) and narrow via `isinstance(phi, Field)` not `hasattr`.
  Reads honestly ("a typed flux field or a raw array"). Conforms to
  `Functional[Field, np.ndarray]`; nothing in S5 assigns it to a typed
  `Functional[Vector,...]` slot so the contravariant-narrowing is not flagged.
  Same concrete-carrier-annotation pattern `MultiplicationOperator.apply`
  (`psi: "TimedFullField"`, reads `psi.bulk.values`) uses.

## GATES (all green, `python -O` canonical; `-O` is the INTERPRETER flag)

- **pyright (CLI, the trusted oracle): 2307 errors / 19 warnings = EXACT
  baseline @ `b404ae1`, 0 net-new, 0 net-new `# type: ignore`.** ⚠ NOTE: the
  brief quoted baseline 2295/19 — STALE (from the S4.5 worktree count); the
  ACTUAL `b404ae1` host-tree baseline is 2307/19 (verified BEFORE writing any
  code). Both new files 0/0 standalone. The +13 delta the brief's 2295 would
  imply is NOT mine — it is the pre-existing #226 host-tree count.
- **S5 suite: 21 passed / 1 skipped** (was 3p/19s PRE-IMPL). The 1 skip = the
  OPTIONAL wrapper (intended). MUTATION-verified: `sum(axis=0→1)` REDS 6 gates
  (both hand-loop legs + RankOne equivalence + both shape/measure gates +
  discriminator), reverted clean.
- **Bit-identity regression: 92 passed / 2 xfailed / 0 FAILED** on
  test_fission_operator + test_multiplication_operator + test_kinf_homogeneous
  (≥2G, NOT 1G) + test_invertible_operator. RankOne untouched → fission
  byte-identical by construction (purely additive 2 new files + 2 export lines +
  1 helper probe line — touches NO SN/CP/operator production code). Baseline
  reds = 7 (#250 SPHERE ×5 + #232 mu_y ×2) live in OTHER files I routed around
  (never all tests/sn, #212 hangs); my subset introduces ZERO new reds.
- **`import orpheus` clean; layer imports 281 passed** (`-O`) + the 2 new-module
  rows pass WITHOUT `-O` (asserts firing, Mode-8). functional.py imports only
  typing + numerics.vector (L1-clean); production_rate_functional.py imports
  numpy + numerics.field (L1) + transport CrossSectionField (L2) — NO L3 edge.
- **Sphinx clean** — full build, only pre-existing `SyntaxWarning`s in 2 test
  files (invalid escape sequences in test_projection_operators / 
  test_fission_operator, triggered by the Nexus graph builder importing test
  modules — NOT mine, NOT rST). NO new rST warning, NO broken xref, NO dup
  label. Both labels in HTML; Nexus graph picked up `math:equation:
  production-rate-functional` (degree 5) + the `ProductionRateFunctional`
  class/method nodes.

## Scope discipline (NOT done, per brief)

- RankOne / FissionOperator NOT rewired (S6, the principled-not-bit-identical
  composition `F = M_χ ∘ ProductionRateFunctional ∘ M_νΣf`).
- `iteration.py` estimators (`_default_production_estimator`/`_default_keff_
  estimator`) UNTOUCHED — bit-identical (C.1 pins). NOT wrapped in a Functional
  object (they take the 4-arg `(L,S,F,ψ)` operator triple, not a lone field →
  NOT `Functional[V]`; the category just NAMES their core).
- SN/CP `compute_production_rate`/`compute_group_production_rate` UNTOUCHED
  (a different volume-weighted group-last reduction; collapsing the 3-way
  duplication is deferred).

## vv classification + NO ERR

- Claim layers: A = CATEGORY/structural (Protocol membership, ref = the
  Protocol defs, non-self-referential via pos+neg+discriminator). B = FLUX-SHAPE/
  value (struct-indep hand double-loop ref; RankOne-equivalence is the DEMARCATED
  equivalence leg). C.1 = regression/value (bit-id inherited). ZERO eigenvalue/
  MMS/convergence claims (those ride the EXISTING test_iteration KEigenvalue gates).
- NO new vv failure mode. NO ERR (additive, nothing caught — the wrong-axis
  mutation is a PROACTIVE gate-teeth check, not a found bug; next free ERR-063
  stays reserved).

## LESSON

A co-vector Protocol `evaluate(x: V) -> R` has a DIFFERENT variance profile than
the operator endomorphism `apply(x: V) -> V`: V is consumed-only (→ declare a
`contravariant=True` input typevar) and R is produced-only (→ `covariant=True`),
where the operator's V is invariant-by-dual-use (and so emits NO variance
warning). Reusing vector.py's invariant `V` in a one-directional Protocol
position reddens pyright with `reportInvalidTypeVarUse` WARNINGS — fix by
declaring properly-variant typevars (the co-vector mirror), NOT by suppression.
And: a `.values` access on the abstract `Vector` carrier (no such member)
narrows via `isinstance(phi, Field)` not `hasattr` — annotate the concrete
`Field | NDArray` (the established `MultiplicationOperator.apply(psi:
TimedFullField)` concrete-carrier pattern), and pyright narrows cleanly.
Extends [[issue-257-s3b-multiplication-operator-closeout]] (the §5.7 sibling
promotion; same additive-bit-identical + 0-net-new-pyright framing; S3b promotes
the coefficient to an OPERATOR, S5 names the field→scalar FUNCTIONAL).
