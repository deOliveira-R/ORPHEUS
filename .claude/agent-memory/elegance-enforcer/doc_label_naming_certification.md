---
name: doc-label-naming-certification
description: Method + durable rulings for certifying descriptive equation :label: naming on the #231 theory corpus (greppability/one-concept-one-spelling discipline)
metadata:
  type: project
---

Reusable when asked to certify NAMING of a batch of equation `:label:` options on
`docs/theory/` (the #231 corpus campaign; ratified disciplines = kebab/no-`eq:`
prefix, descriptive-not-positional, page-family prefix + one word-order per family,
project-global self-locating namespace, internal-consistency = one concept one
spelling).

**Method that worked (fast, high-coverage):**
1. Ground-truth the diff: `git diff -U0 -- docs/theory | grep '^+ *:label:'` must be
   the ONLY additions; `grep '^-'` must be empty (pure-additive claim). Cross-check
   the `docs/verification/matrix.rst` orphan count delta = N-new-labels exactly (it
   corroborates zero label collisions: 108→243 for 135 adds).
2. Master corpus: `grep -rnE "^\s+:label:" docs/theory` → `file :: name`. Prefix
   histogram (`sed 's/-.*//' | uniq -c`) reveals the page families.
3. Collision scan = grep the corpus for CONCEPT TOKENS (balance, source, scatter,
   multiplication, adjoint, transpose, si-, cell-balance, keff-update…) and read the
   equations behind any same-token pair. This is the highest-value check; parallel
   batches can't see each other.
4. Section-anchor `-eq` rule: a label `foo-eq` is JUSTIFIED iff a `.. _foo:` anchor
   exists on the page (grep `_foo:`); it disambiguates eq-from-section. Verify BOTH
   members of a 1d/2d pair carry the same rationale.
5. Mechanics via diff hunks (`git diff -U3`): `:label:` sits immediately under
   `.. math::`, indent tracks directive depth (list-nested = 5-6 spaces), aligned
   `&\;=\;` block = ONE label. Python one-liner catches double-labels / misplacement.

**Durable ruling — the "one domain term, two legitimate concepts" pattern (→ SHOULD-
CONSIDER, not MUST-FIX):** `mg-multiplication-operator` (K=A⁻¹F, neutron-multiplication)
vs pre-existing `multiplication-operator-{embedding,action}` (M[f], the multiplier
algebra) share the "multiplication-operator" token for DIFFERENT operators. All three
VIOLATION legs are present (bug-habitat = grep-ambiguation / wrong `:eq:` target;
not-coextensive = genuinely two operators; live-verified). BUT it is a SHOULD-CONSIDER
because: (a) the prose overload pre-exists corpus-wide (~5 pages call K "the
multiplication operator" — domain-inherent, not backfill-introduced); (b) the label
strings are distinct & each disambiguated (`mg-` geo-prefix vs `-embedding/-action`);
(c) the fix (`mg-keff-operator`) TRADES one inconsistency (token overload) for another
(label-vs-surrounding-prose divergence). When a rename swaps one inconsistency for
another and the term is entrenched domain vocab, it's a terminology-owner judgment
call, not a mechanical fix I certify as required.

**Non-collisions that LOOK like collisions (checked, cleared):** `self-collision-
probability-slab` (P_ii, the probability ∈[0,1]) vs `self-slab` (r_ii, reduced CP with
the Σt·t term) — different quantities, new name is the more-explicit one.
`region-areas-pincell` (absolute areas cm², MoC) vs `pin-cell-volume-fractions`
(dimensionless fractions, homogenization) — different quantities/pages (only nit:
pincell vs pin-cell hyphenation). `si-` = source-iteration consistently whether prefix
(SN pages) or infix (`inverse-driver-si-update`, inverse-family page) — word position
differs by page-family, which the discipline ALLOWS.

**References batches (peierls-*/fn-*/wm72-*) are the gold standard** — dense
`<geometry>-<concept>` / `<region>-<concept>` families, consistent word order,
paper-`-eqNN` completions (join existing `-eqNN` on literature-mirror pages; rule these
WITH paper-number context, NOT as bare positional). Bare-dominant ratios on derivation-
mirror pages (2 labeled / 31 bare) are by-design, not under-labeling.

**Cross-page governing-eq parallel is CORRECT, not a twin path:** the transport eq is
labeled on 6 pages (each a different form: FN moment-space `fn-method-bte`, Peierls 3-D
`peierls-transport-equation-3d`, MoC `transport-equation`, methods-index
`methods-local-transport-equation`, …). Same physics, different per-method writing →
each page self-contained. Only wrinkle = the concept token splits `bte` (spectral-
reference cluster: galerkin_spectral, fn_method) vs `transport-equation` (elsewhere);
the split pre-dates the backfill and `fn-method-bte` correctly mirrors its sibling
`galerkin-spectral-bte` → SHOULD-CONSIDER at most.
