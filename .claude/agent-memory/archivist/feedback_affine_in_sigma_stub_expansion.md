---
name: affine-in-sigma-stub-expansion
description: #240 Phase 2 Step B — expand the loss-rep-removal-form-matvec stub into rich narrative. A principled-equivalence (NOT bug) carve: composite owns its matvec via loss_action(self.sigma); affine-in-σ leaf-sum is value-correct-by-coincidence; removal form C(σ_r) makes the two-source distinction observable.
metadata:
  type: feedback
---

#240 Phase 2 Step B stub→rich expansion in
`docs/theory/loss_representations.rst` (`.. _loss-rep-removal-form-matvec:`).
A **principled-equivalence carve**, not a bug fix — the apply value never
shipped wrong. The narrative's job is to explain WHY a value-correct path was
still wrong-by-construction, and why the carve is prophylactic (no production
caller of the removal form yet — #200).

**Why:** the deepest source by far was the TEST FILE
(`tests/sn/operators/test_removal_form_matvec_sweep.py`), not the closeout
memo. Its module docstring carries the full "PREMISE CORRECTION" — the
affine-in-σ algebra, the four gate categories (a teeth / b value-ground /
c production-σ / d negative), the vv-principles compliance, the Mode-9
gate-design-degeneracy nuance. The closeout memo gave the re-baseline ULP
table + blast-radius split; the test docstring gave the math + the gate
taxonomy. The production docstrings (`InvertibleOperator.apply` lines
1321-1366, `StreamingOperator.apply` 775-841) gave the prose seed verbatim.

**How to apply (the shape that worked — 5 H3 subsections under the H2):**

1. **Lead with the framing**: the composite has its own matvec too; getting it
   to single-source σ from where solve does is the carve. State up front it is
   principled-equivalence (no wrong value shipped).
2. **"The σ parameter is now explicit"** — the protocol signature flip
   `(operator, psi)`→`(sigma, psi)`, symmetric with `sweep(Q, sig_t)`; caller
   single-sources. Frame as the one-instance theorem extended to the diagonal.
3. **"The affine-in-σ structure"** — the LOAD-BEARING derivation. Two NEW eq
   labels: `loss-rep-affine-cell` (M(σ)ψ = denom·ψ̄ − numer, denom = Σs_a + σ,
   ψ̄ KNOWN in apply) → `loss-rep-affine` (M(σ)ψ = streaming_action + σ·ψ).
   Then the CONTRAST: sweep has ψ̄ = numer/denom → σ in the DENOMINATOR →
   non-affine. The asymmetry is WHY the round-trip can't catch a σ-routing
   error in apply alone. Then `loss-rep-leaf-sum`: L.apply gives
   streaming_action (Res-A subtracts σ_t), C.apply gives σ_r·ψ, sum = M(σ_r)ψ
   — right value, WRONG source.
4. **"The two-source hazard"** — a 2-row table (solve sources σ from C;
   inherited apply sources from L). Agree only because production builds L,C
   from same σ_t. Show the override code-block (loss_action(self.sigma)).
   Pattern 2 single-source. Note multi-D adjoint = deferral RAISE not silent.
5. **"The removal form makes it observable"** — `loss-rep-removal-sigma`
   (σ_r = σ_t − Σ_s0). NO production caller yet (#200 fold not wired); σ_r-sweep
   is NOT a correct accelerator (#215, 46-56% aniso errors); carve is
   prophylactic — trap only if a future L-leaf becomes non-affine in σ.
6. **"Verification"** — the four gates (a/b/c/d) lifted from the test file,
   mapped to the #240-Step-A value-preserving-re-association shape + the
   vv 3 criteria. The blast-radius split (APPLY snapshots re-baseline ≤5 ULP,
   SWEEP/SOLVE stay bit-identical; Krylov-2D φ ≡ SI-2D φ to 3.9e-12) IS the
   structural-independence evidence. k_inf eigenvalue ground = xfail until #200.

**Disciplines confirmed:**

- **Carve SHARPENS an existing eq-label, no new :label: for the claim.** The
  teeth gate `verifies("loss-rep-resolution-a")` — the same Resolution-A glue,
  now single-sourcing σ. I added 4 DERIVATION/representational eq-labels
  (affine-cell, affine, leaf-sum, removal-sigma) but NONE is a new solver
  claim — they are the algebra OF the existing claim. (Did not vv-status tag
  them — they sit inside a verified narrative chain, derivation-category like
  scanmarch-solve/apply which are untagged on the same page.)
- **[AdamsLarsen2002]_ is defined in discrete_ordinates.rst** → cite as plain
  text "Adams & Larsen 2002" inline (NOT `[AdamsLarsen2002]_`) to avoid the
  cross-doc duplicate-citation warning. The page's own References table already
  uses plain "Adams & Larsen (2002)".
- **This entire page's :meth:/:class: code-xrefs render PLAIN TEXT** (0 resolve
  to hyperlinks — confirmed by grepping the built HTML for `<a href>` anchors:
  even pre-existing `StreamingOperator.apply` renders as plain `…apply()`).
  `orpheus.sn.operator` is not member-autodoc'd with resolvable anchors here,
  and adding automodule would dup-label (operator.py docstrings carry
  `:math::label:`). USE the page's existing `:meth:`~orpheus.sn.operator.…``
  path syntax for consistency — do NOT half-surface 1-2 leaves.
- **Build gate**: `-E -W --keep-going` baseline = 1 (the standing
  `Mesh1D.from_geometry :paramref:` ERROR). Post-edit identical = 1. The
  acceptance gate is count-unchanged, NOT zero. "build finished with problems,
  1 warning" refers solely to that pre-existing item.
- Em-dash titles: underline length by `len(title)` code points (verified all
  5 H3 underlines ≥ title len). H2 `----` parent → H3 `~~~~` is valid here.

**Quality self-assessment (1-5):**

| Dimension | Score | Note |
|-----------|:---:|------|
| Derivation depth | 5 | full affine-in-σ derivation, sweep-vs-apply contrast, leaf-sum collapse, all intermediate steps |
| Cross-references | 4 | every :eq:/:ref: resolves; :meth:/:class: render plain-text by page convention (not my regression) |
| Numerical evidence | 5 | ULP re-baseline table by geometry (2/2/4/5), 3.9e-12 Krylov≡SI cross-check, gate-by-gate |
| Failed approaches | 5 | the leak's value-correctness, the round-trip gate's Mode-9 degeneracy, sphere round-trip divergence (ERR-058), bitexact-spec re-baseline |
| Code traceability | 5 | every gate named, production overrides + StreamingOperator.apply + CollisionOperator cited |
| Derivation source | 5 | algebra from test-file docstring + production docstrings (the algebra-of-record for an apply-direction carve is the operator code, not a SymPy module) |

Weakest = cross-references (4), but only because the page-wide plain-text
code-xref convention is out of scope for a stub expansion (surfacing the whole
operator package is its own architectural-docs task).
