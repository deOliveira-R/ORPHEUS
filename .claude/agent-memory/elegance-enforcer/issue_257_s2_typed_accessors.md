---
name: issue-257-s2-typed-accessors
description: #257 S2 PASS-WITH-NITS — typed CrossSectionField accessors on MaterialXSField; symmetric-build vs Pattern-6 adjudication, zero-copy lens pinned by `is`-test
metadata:
  type: project
---

# #257 S2 — MaterialXSField typed CrossSectionField accessors

PASS-WITH-NITS (`feature/field-typed-operator-algebra`, HEAD `505e1b7`, pre-commit, pure-addition,
bit-identical). Diff = 1 import (`material_xs_field.py:93`, runtime not TYPE_CHECKING) + 3 `@property`
accessors (`total_cross_section_field` / `absorption_cross_section_field` / `fission_production_field` at :398-408 — ASYMMETRIC: the fission lens drops the `_cross_section` infix because νΣf is a production quantity, verified vs source @S5 review; the earlier "fission_production_cross_section_field" spelling here was WRONG), each
`CrossSectionField.from_mesh(<raw cached view>, self.mesh)`. The FIELD-SIDE half of the S3 promotion
`C = M[σ_t]`. Test `tests/sn/test_material_xs_field_typed.py` 10✓ in 0.48s. No new pyright errors;
2 pre-existing `cells_by_material` errors untouched (out of diff range).

**Why:** S2 of the #257 coefficient-field campaign (plan `.claude/plans/issue_257_coefficient_field_promotion.md`);
documented as "typed accessors ALONGSIDE raw ndarray views, existing consumers untouched, bit-identical."

**How to apply:** the 4 judgment-call verdicts + the reusable adjudication below.

## The 4 judgment calls — verdicts

- **J4 Layering (PASS/KEEP):** L3→L2 runtime import is layer-legal AND cycle-free — VERIFIED by fresh
  interpreter in BOTH import orders (exit 0). The `transport.fields → sn` back-edge exists ONLY under
  `TYPE_CHECKING` (`_bases.py:79`, the `SNMesh` annotation). Runtime (not TYPE_CHECKING) is CORRECT
  because `CrossSectionField` is USED at runtime in the accessor body (`.from_mesh(...)`), not merely
  annotated — TYPE_CHECKING would NameError on first access. Clean, no deferred-import band-aid
  (contrast #245's numerics→sn cycle break).

- **J3 Construct-on-access vs cache (PASS/KEEP do-NOT-cache):** the EXPENSIVE step (per-cell
  `assemble_cell_xs`) is ALREADY cached behind the raw view (`_ensure_cell_views`, run-once); the
  accessor wraps that cached array zero-copy (`from_mesh:338` passes `values=values` straight to the
  frozen-dataclass ctor, NO copy — STRUCTURAL, verified). Caching the FIELD object would mint a 2nd
  mutable cache slot needing lockstep invalidation with `_sig_t_cell` → a state-coherence obligation
  (two-spelling-of-one-value divergence habitat). Per-call construction REMOVES that obligation. S3's
  `MultiplicationOperator` reads the field once at operator-build time (not per matvec) → no hot
  pressure. STANDING PATTERN: a frozen-dataclass lens built-on-access over an already-cached array is
  MORE elegant than a derived cache — it cannot go stale.

- **J1 `*_field` suffix naming (PASS/KEEP this stage):** the two accessors return DIFFERENT types
  (`np.ndarray` vs `CrossSectionField`) on the same σ → MUST differ in name while coexisting; `_field`
  is the conventional self-documenting transitional sibling. CAVEAT (record, no action): the suffix is
  correct BECAUSE temporary — the terminal elegant name is UNSUFFIXED `total_cross_section` returning
  `CrossSectionField`, raw retiring or becoming the exception. The S-stage that flips the primary path
  must RENAME, not accumulate a 3rd spelling (anti-#11 removal trigger lives in the plan, S3+).

- **J2 Typing σ_a with NO near-term consumer (CONCERN→KEEP-with-rationale-NIT):** the ONE axis where
  shortest-elegance and symmetry-elegance pull apart. σ_t→S3, νΣf→S6, σ_a→NO consumer through S9.
  ⭐ THE ADJUDICATION (reusable): Pattern-6 defer-until-evidence applies to ABSTRACTIONS, not to a 3rd
  INSTANCE of an already-proven pattern. The wrap-pattern (`<raw view> → CrossSectionField.from_mesh`)
  has its rule-of-two met IN-DIFF (σ_t + νΣf both have consumers); σ_a is a 3rd instance of that proven
  pattern, NOT a 3rd speculative abstraction. "Symmetry-in-math ⟹ symmetry-in-code": σ_t/σ_a/νΣf are
  co-equal macroscopic XS (the S1 leaf docstring already enumerates the family); an asymmetric
  "two-typed-one-raw" surface is ITSELF a smell that invites the future σ_a consumer to open-code the
  wrap inline (Pattern-2 habitat) or pollute its diff with the accessor. KEEP closes the habitat; DEFER
  opens one. Cost of KEEP = 3 lines that can't diverge (all 3 pinned by the SAME parametrized `is`-test
  via the `ACCESSORS` table → σ_a is NOT untested dead weight). NIT (do-now): the comment block
  (:378-387) justifies the typed-lens but NOT why-all-three; add one sentence so a future reader does
  not DRY-delete σ_a (anti-#11 applied to a build-ahead accessor — the symmetric build is fine, the
  REASON must be in the code).

## STANDING TELLS / reusable lessons

- ⭐ **Symmetric-build of a co-equal family is NOT a Pattern-6 violation when the wrap-pattern's
  rule-of-two is already met in-diff.** Discriminator: is the Nth member a new ABSTRACTION (defer) or
  the Nth INSTANCE of a proven pattern over a domain-symmetric set (build, closes the open-code
  habitat)? Demand the rationale-comment so the build-ahead member survives a future DRY-sweep.
- **`.values is raw_view` (`is`, NOT `==`/`array_equal`) is the correct gate for a zero-copy LENS
  claim.** An `==` test passes even if a copy creeps in and silently breaks the "no second
  representation" invariant. The team uses this correctly here — the byte-identity pin a twin-adjacent
  carve demands. Negative control `not isinstance(raw, CrossSectionField)` proves pure-additive.
- **The team's #1 institutional smell (phased carve leaves twin DELIVERY single-sourced only at the
  operator) does NOT fire on an additive typed-LENS stage** — S2 introduces no second COMPUTATION of
  any XS; every accessor wraps the ONE cached ndarray. Distinguish a typed-lens-over-cache (no twin)
  from a second realization route (the dangerous twin).
- **Comment shorthand "(1/cm, the cone)" risks reading as a per-field ≥0 GUARANTEE** — the S1 leaf is
  explicit the cone is a PROPERTY not a type invariant (signed Σ−Σ′ is valid). Micro-nit; soften to
  "physical XS, so cone-valued." (Links the #257 S1 J2 σ≥0-not-enforced ruling.)
