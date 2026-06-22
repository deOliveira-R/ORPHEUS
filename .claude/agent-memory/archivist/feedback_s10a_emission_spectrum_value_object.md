---
name: s10a-emission-spectrum-value-object
description: #257 S10a — expand the EmissionSpectrum value-object stub (cross_section_data.rst) into rich narrative. A DATA-LAYER value-object (not a solver/operator carve): probability-simplex law, source-not-per-cell enforcement, PRODUCTION-keying (νΣf>0) vs fissionability (Σf>0), behavior-neutral precursor, S10b convex-combination seam. SOFTWARE-INVARIANT not solver claim → foundation gates no verifies.
metadata:
  type: feedback
---

Expanded the S10a stub "The Fission Emission Spectrum as a Validated
Value-Object" (label `emission-spectrum-simplex-law`) in
`docs/theory/cross_section_data.rst` — removed the `.. todo:: Archivist`
and replaced with 6 `-`-level (H2) subsections. NOT a close-out, NOT a
retirement; a DATA-LAYER value-object documentation task. NOT committed.

## What was distinctive about this task

- **Data-layer value-object, not a solver/operator carve.** Most recent
  archivist work was SN operator-algebra carves (apply-solve, decomp-leaf,
  taxonomy). This is the `orpheus/data/` probability-simplex newtype. The
  invariant is a SOFTWARE invariant (a probability-distribution law), NOT
  a theory-page solver equation → the gates are `@pytest.mark.foundation`
  with NO `verifies(...)` edge. The eq-labels I minted are documented-only
  structural/representational identities (the simplex law, the
  fission-source consumption form, the χ_mix convex-combination), not
  flux/eigenvalue claims.

- **The stub already carried the keying correction.** The brief warned a
  2026-06-21 review re-keyed the guard SigF→PRODUCTION; the method-implementer
  had ALREADY updated the stub's label-math + TODO point (5) to production-
  keying. My job was the rich-narrative expansion, NOT re-keying. VERIFY-THE-
  STUB-FIRST paid off (it matched the closeout memo's REFINEMENT block).

- **The closeout memo's `## ⭐ S10a REFINEMENT` block was the richest prose
  seed** — it carried the production-keying derivation, the billiard-fixture
  resolution (SigF=0 honest, hack RETIRED), the two-predicate divergence,
  the convex-combination property, the tol probe (U_235 sum=1.0000000000000002,
  2.22e-16), the DD-byte-id inertness proof. The production DOCSTRINGS
  (emission_spectrum.py module + enforce_emission_spectrum + the two
  __post_init__ comments + Isotope.is_fissile/is_producing docstrings) were
  the SECOND seed — they already articulate WHY production-keying. Read both,
  cross-check against LIVE code (mixture.py:211-214 first-fissile χ = the
  S10b seam; billiard fixture line 106-107 SigF=zeros/SigP=copy).

## The 6 subsections (the shape that worked)

1. **Why the law is enforced at the source (and never per-cell)** — the
   single-source `enforce_emission_spectrum` helper (Pattern 2, the two
   __post_init__ bodies collapsed); WHY-not-per-cell (broadcast χ[None,:]
   legitimately holds χ=0 in non-fissile cells → per-cell validator can't
   distinguish genuine null from malformed → discriminating predicate
   (is_producing) lives at the material); the value-object-vs-#263-SpectrumField
   boundary (validated newtype, NOT flux-algebra member, no __add__/morphism).
2. **The tolerance rationale: atol=1e-12** — 3-row band table (real GENDF FP
   residual 2.2e-16 MUST-pass / 1e-12 threshold / 1e-6 physical-error MUST-fail);
   rtol=0 rationale; assert_null EXACT-zero (constructed not summed); the
   two-clause independence ([1.2,-0.2] Σ=1 but negative).
3. **Why the law keys on production, not fissionability** — THE LOAD-BEARING
   subsection. Mint eq-label `emission-spectrum-fission-source`
   (q^fiss_g = χ_g Σ νΣ_f,g' φ_g') = the WHY (χ never used alone, always ×νΣ_f).
   Two-predicate table (is_fissile=Σf>0 "can fission?" → compute_macro_xs
   fissile_indices; is_producing=νΣf>0 "does it emit?" → the χ law). COINCIDE
   for real compute_macro_xs (νΣf>0⟺Σf>0, proven TestRealGendfConstructs
   U_235/O_016/H_001), DIVERGE for synthetic SigP>0/SigF=0 (trajectory-resolvent
   billiard). The hack-RETIRED narrative (earlier draft keyed is_fissile +
   patched billiard SigF stand-in; re-keying retired it, fixture now SigF=0
   honest) + the isolation leg test_non_producing_nonzero_raises (sigF>0 nubar=0).
   GENERALIZED lesson: key a guard-at-source on the predicate matching WHY the
   guarded value EXISTS.
4. **The behavior-neutral precursor** — dead non-fissile χ in synthetic
   xs-library zeroed; inert because χ·νΣf≡0 for non-producing; PROVEN (not
   asserted) byte-identical by test_dd_regression.py (no snapshot moved = χ
   was inert); bit-identity-vs-principled-equivalence framing (regression
   contract IS the right gate when the change alters nothing); the L20 blast-
   radius scope lesson (every DIRECT constructor not just factory path; AST scan).
5. **The S10b seam** — mixture.py:211-214 first-fissile χ verbatim = the seam;
   mint eq-label `emission-spectrum-chi-mix` (production-weighted convex combo);
   THE decisive property: convex combination of simplices IS a simplex → χ_mix
   passes assert_normalized VERBATIM (no new law) — the value-object designed
   so S10b inherits its floor free.
6. **Verification chain** — 2-row gate table (intrinsic + container) + the
   real-GENDF no-red proof + Mode-8 (-O safe, pytest.raises/np.testing/_require
   not bare assert) + DD-regression third leg.

## Cross-ref reality (matches AGENT.md standing convention, NOT a new finding)

- `:eq:` (emission-spectrum-fission-source/-chi-mix/-simplex) + intra-doc
  `:ref:` (-production-keying/-verification) + cross-doc `:ref:`
  (synthetic-xs-library → verification.rst:108) ALL resolve to real hrefs.
  Grep-gated every `:eq:`/`:ref:` target before building.
- Code-xrefs SPLIT by automodule coverage: `Mixture`/`compute_macro_xs`
  resolve (api/data.rst automodules macro_xs.mixture) → real href;
  `EmissionSpectrum`/`enforce_emission_spectrum`/`Isotope`/`micro_xs.isotope`
  render PLAIN STYLED TEXT (no warning, non-nitpicky config) because
  emission_spectrum + micro_xs.isotope are NOT automodule'd. This is EXACTLY
  how the pre-existing stub's refs already rendered — consistent, NOT half-
  surfacing. DEFERRED the automodule surfacing (its own arch task; +
  emission_spectrum.py module docstring has a `.. math::` block → automodule
  would risk duplicate-content). Did NOT half-surface 1-2 leaves.

## Build gate

- FORCED `-W -E` rebuild. EXIT=1 from the LONE standing baseline (the
  `mesh.py` Mesh1D.from_geometry `:paramref:` ERROR, needs sphinx-paramlinks,
  out of scope). WARNING/ERROR/CRITICAL grep = that ONE ERROR, unchanged.
  The LD `verifies(...) — skipping` lines are INFO severity (don't appear in
  the W/E/C grep). ZERO `emission` mentions in the build log → my edit added
  ZERO new notices. Baseline count UNCHANGED (1).
- Em-dash/code-point underline check: 1 H2 underline was 1 char too long
  (title 58 code-points, underline 59) → fixed via python len()-check on the
  title vs `^-{5,}$` underline (NOT wc -c). χ in the title is 1 code-point.

## Quality self-score (rubric)

- Derivation-depth 5 (full WHY chains: per-cell-discards-the-predicate
  argument, production-vs-fissionability from the consumption form, convex-
  combination-is-a-simplex proof).
- Cross-references 5 (every gate file + both __post_init__ + Mixture/
  compute_macro_xs hyperlinked + synthetic-xs-library cross-doc).
- Numerical-evidence 5 (the 3-row tol band table with the real 2.22e-16
  probe, the DD-byte-id inertness proof, U_235/O_016/H_001 GENDF legs).
- Failed-approaches 5 (the is_fissile-keyed earlier draft + billiard SigF
  stand-in hack documented as RETIRED-history).
- Code-traceability 5 (3 minted eq-labels + the mixture.py:211-214 seam line
  quoted verbatim + every cited test by path).
- Derivation-source N/A (software invariant verified by foundation gates +
  byte-id regression — the correct source for a probability-simplex law;
  no SymPy module owed).

Extends [[quadrature-invariant-subsection]] (same documented-vv-status
software-invariant-not-solver-claim pattern). The method-implementer
closeout this expands: `issue_257_s10a_emission_spectrum_closeout.md`.
