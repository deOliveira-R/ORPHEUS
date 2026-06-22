---
name: issue-257-s1-coefficient-fields
description: PASS-WITH-NITS verdict on #257 S1 (CoefficientField + CoefficientRole marker mixin + CrossSectionField cone + SpectrumField simplex) on branch feature/field-typed-operator-algebra; the cone/simplex mirror of flux/displacement, invariant-as-single-gate, role-marker-as-taxonomy.
metadata:
  type: project
---

# #257 S1 — CoefficientField + operator-as-promotion (Stage 1: fields only)

PASS-WITH-NITS, `feature/field-typed-operator-algebra`, pre-commit pure-addition
(2 new typed leaves + role mixin + 2 unit constants + intrinsic-property test;
NO existing prod changed). Grand report §5.5–5.7 promotion: cross sections are
`CoefficientField`s, operators are those fields PROMOTED (S3, not S1). Frame memo
`.claude/agent-memory/cross-domain-attacker/coefficient_field_promotion_frames.md`
Frame 2 (coefficient = cone/simplex, NOT flux torsor).

## The four judgment calls — all KEEP (correct)

1. **`CoefficientRole` dunder-EMPTY marker mixin = KEEP, NOT Pattern-6.** The
   2nd concrete consumer exists IN THE SAME DIFF (both leaves carry it) → rule-of-two
   satisfied at landing. Earned on 2 LIVE-IN-S1 grounds (not the deferred `*`):
   (a) taxonomy — `issubclass(leaf, CoefficientRole)` is the discrimination axis
   the test asserts (`not issubclass(leaf, FluxRole)`), MRO-slot mirrors FluxRole
   (mixed in BEFORE storage base); (b) named home for S3's `M_f@M_g=M_{f·g}`
   field·field `*` WITH a tracked trigger (issue+plan cited → discharges anti-#11).
   FOLD would re-create attribute/stringly dispatch (frame Smell #16) + the S3 `*`
   arrives homeless → twin-path. LESSON: a dunder-empty marker is the RARE case
   where a long docstring on an empty body is ELEGANT (content = "absence of the
   gate", cannot self-document via code) — not a smell.

2. **`CrossSectionField` does NOT enforce σ≥0 (cone-as-property) = KEEP; the
   asymmetry vs SpectrumField IS principled.** Discriminator = the consuming
   algebra's domain. Multiplier embedding `M:L^∞→B(L²)` is LINEAR ON THE FULL
   SIGNED vector space (homomorphism law `M_{af+bg}=aM_f+bM_g` becomes the S3
   suite); `Σ−Σ′` perturbation is a legit signed element. Enforcing σ≥0 would make
   inherited `__sub__` SOMETIMES-RAISE-ON-VALUES (worse illegal-states than
   "cone = tested property"). Cone enforced at the DATA boundary (XS reader), not
   every intermediate (anti-#8). Asymmetry principled: signed Σ is MEANINGFUL,
   off-simplex χ is NEVER meaningful (always Σ=2) → different meaningfulness →
   different enforcement locus. Pattern-4 applied WITH DISCRIMINATION not uniformly.

3. **`χ+χ` rejected via EMERGENT `__post_init__`, not explicit `__add__` = KEEP
   (the invariant IS the gate).** VERIFIED chain: `Field.__add__`→`_check_partner`
   (passes)→`replace()` on frozen dc RE-RUNS `__post_init__`→simplex check Σ=2→
   ValueError("simplex"). MORE elegant than FluxRole's explicit `__add__` TypeError:
   FluxRole MUST hand-write the gate (flux+flux is type-coherent-but-void, no
   construction invariant catches it); the simplex has ONE source of truth (the
   `__post_init__` invariant) that EVERY illegal-χ path routes through. An explicit
   `__add__` override = a 2nd statement of the simplex law = Pattern-2 dup (relax
   the tol in one, the other goes stale). Discoverability discharged: class
   docstring + the ValueError message both explain in domain language. STANDING
   PATTERN: when a construction invariant exists, prefer emergent-gate-via-`replace`
   over a hand-written dunder override — it single-sources the law.

4. **`SpectrumField.mix` parallels `FluxRole.affine_combination` = KEEP, do NOT
   extract.** Shared blend SKELETON (empty-check / float-coeffs / Σλ-tol /
   `_check_partner` loop / `sum(start=coeffs[0]*first.values)` seed / `replace`)
   but the CONSTRAINT defines the op: affine_combination = Σλ=1 SIGNED-λ-OK
   (relaxation ω>1 has λ<0, legal on affine subspace); mix = Σλ=1 AND λ≥0 (CONVEX,
   λ<0 pushes off simplex). Extracting `_blend(pairs, *, require_convex: bool)` =
   anti-#3 BOOLEAN-FLAG dispatch (the flag is the habitat for "wrong flag → silent
   affine blend where convex needed → off-simplex χ, ERR-039 class"). The two also
   live on OPPOSITE role families (FluxRole vs CoefficientRole-family) → a shared
   helper forces a cross-family dep the whole design keeps complementary-but-separate.
   COLLAPSE TRIGGER = a THIRD blend with the SAME constraint, NOT the 2nd with a
   DIFFERENT one (Metz: dup cheaper than wrong abstraction).

## Nits (do-now in same commit)
- **NIT-1 (Cardinal-Rule-3 doc accuracy)**: `mix` docstring (`:142`) says "analogue
  of affine_combination with the extra nonnegativity restriction" — UNDERSELLS:
  affine→convex CHANGES the operation's identity, and that is WHY they're not
  single-sourced. An "it's basically X" comment is the documented invite to the
  J4 wrong-refactor. Add a sentence: constraints define different ops, do not unify.
- **NIT-2 (twin-arithmetic test gap)**: `isinstance(_chi, Vector)` is True
  (`@runtime_checkable` checks METHOD PRESENCE only — all 4 dunders inherited) but
  `χ+χ`/`χ−χ` RAISE (simplex-gated). Green "is a Vector" next to a type whose `+`
  raises = latent contradiction; the day someone wires χ into a generic Vector
  accumulator it type-checks + passes isinstance + raises deep in the loop. Add a
  1-line comment: Vector conformance is structural-presence-only for SpectrumField,
  additive ops are invariant-gated (xref `test_addition_leaves_simplex`).
- **NIT-3 (KNOWN ITEM — CONFIRMED CORRECT)**: units.py:183 `CROSS_SECTION_UNITS =
  UREG.cm**-1` trips pint `**`-stub PlainUnit→Unit. grep-VERIFIED units.py carries
  ZERO `# type: ignore` — the 4 pre-existing constants (:156/:162/:166/:170) ship
  the SAME error bare. Matching the convention (5 identical bare) is RIGHT; casting
  only :183 = 2 spellings of one situation in one file (worse Pattern-7 drift);
  casting all 5 = scope-creep on a pure-add PR. → FILE a uniform-cleanup issue
  (module:numerics, type:improvement), don't patch.

## Minor (no action)
- NIT-4: `_SIMPLEX_TOL` (named, `~ng·eps` rationale, single-sourced for both simplex
  + Σλ checks) is EXEMPLARY. Sub-nit: same `1e-12` literal is bare in
  `affine_combination:135` (the affine Σλ check) → 2 spellings of partition-of-unity
  tol across the 2 role families. NOT do-now (independently defensible, predates PR);
  real twin on a 3rd Σλ=1 check.
- NIT-5: both leaves correctly mirror ScalarFlux's
  `@dataclass(frozen=True, eq=False, kw_only=True, repr=False)`, inherit units-aware
  `Field.__repr__` + no-content-eq. Consistent.

## Arch-opp (S3 watch, not now)
Frame memo leans to `CoefficientField(ScalarField)` ABC holding the commutative-algebra
machinery WHEN field·field `*` lands. S1 correctly minted role-marker-only (defers the
base decision). S3 reviewer: decide whether `*` lives on `CoefficientRole` or a new
`CoefficientField` base. No action S1.
