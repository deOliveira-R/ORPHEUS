---
name: si-rate-recovery-docs
description: Pattern for documenting an SI within-group splitting + a partial-feature rate recovery (Phase 3 boundary Gauss-Seidel) in discrete_ordinates.rst — orphan-label resolution, honest-scope discipline (boundary-only B vs scattering c-mode), the failed σ_r-fold trap, ERR-056 shared-face rule, and the Case–Zweifel-vs-within-group c nuance
metadata:
  type: feedback
---

# Documenting the Phase 3 SI boundary-Gauss-Seidel rate recovery (commit 1cbb383)

Task: add a `:label: si-spectral-rate` block (resolves an orphan
`@pytest.mark.verifies`) AND the full within-group SI story to
`docs/theory/discrete_ordinates.rst`. Shape that worked — a new H2
subsection `_si-within-group-splitting:` inserted INTO the existing
`sn-iteration-primitives` H1 section, between the `SourceIteration:
discrete fixed-point realisation` and `KEigenvalue: outer power
iteration` subsections (the existing SI subsection already had an
UNLABELED ρ=c sentence — I added a forward `:ref:`+`:eq:` pointer into
the new section rather than relabel its loose paragraph).

**Section anatomy (9 parts, ~440 LoC):**
1. Lead paragraph naming the 3 code anchors (`sweep_schedule` module,
   `_GaussSeidelResolvent`, `solve_sn_fixed_source` `inner_schedule=`).
2. `.. admonition:: Key Facts (SI rate) :class: tip` — 5 bullets:
   ρ_J=c iteration matrix `(L+C)⁻¹(S+B)`; boundary-G-S ~0.86–0.92×
   (NOT c²-halving); scattering c-mode = Krylov/DSA only; fixed-point
   invariant; 1-D falls back to Jacobi.
3. **Four-operator within-group eq** `si-within-group-operator-eq`
   `(L+C-S-B)ψ=q` with per-operator bullet glossary (S=Σ_s0·P_iso
   isotropic-projection, B=trace-only).
4. **Splitting** `si-jacobi-fixed-point` `(L+C)⁻¹(Sψ_n+Bψ_n+q)`.
5. **`si-spectral-rate`** (THE mandatory label): ρ_J=c=max Σ_s/Σ_t AND
   the n_Jacobi≈log(ε)/log(c) right-hand identity in ONE eq block; then
   a `.. note::` on the Case–Zweifel-vs-within-group c nuance (the test
   uses `Mixture.scattering_ratio` = (Σ_s+νΣ_f)/Σ_t which is SLIGHTLY
   LARGER than the within-group Σ_s0/Σ_t in the eq — flag it explicitly,
   cite the 0.6–1.2 band + 655-vs-728 measured/predicted ratio 0.90).
6. **Jacobi vs G-S** 2-row schedule list-table + selection prose
   (`_select_si_resolvent` → `(GaussSeidelResolvent,(S,))` vs `(L+C,(S,B))`).
7. `.. warning::` HONEST SCOPE (load-bearing per user): boundary-only,
   ~0.86–0.92× NOT c²-halving, scattering=Krylov(302 vs 860)/DSA(#2).
8. **σ_r-fold trap** (#215): `si-sigma-r-fold-mismatch` eq (diagonal
   removal σ_r·I ≠ isotropic-projection Σ_s0·P_iso) + a 3-row failure
   table (reflective exact / anisotropic 46–56% / "exact" variant
   DIVERGES Σ_s0/σ_r≈39) + the foldable_part/residual_part = DSA-input-
   not-sweep-fold conclusion.
9. **ERR-056** shared-face rule: axis-aligned (product, 1 face/octant)
   immune vs diagonal (lebedev/level_symmetric, 2 faces/octant) needs
   LAST-outflowing-group reflect; white/vacuum excluded.
10. **Numerical evidence** 4-row list-table (box 697→641 0.92×; slab
    no-op 655=655; vacuum 128=128 negative control; Jacobi-vs-Krylov
    2.85×) + eigenvalue-path total_inner_iterations seam (SI 371 /
    Krylov 310) + DSA-spike 8–21× comparison.

**Orphan-label resolution mechanics:** the test file already carried
`@pytest.mark.verifies("si-spectral-rate")` on an ORPHAN (label didn't
exist). Defining `:label: si-spectral-rate` ONCE resolves it. Verified
via `grep -c "^   :label: si-spectral-rate$"` == 1 + the built HTML
anchor `equation-si-spectral-rate` present. The doctree is the source
of truth for orphan-resolution, not a Nexus query (graph.db absent).

**vv-status / vocabulary discipline (Directive 5):** the rate target
is a CLOSED-FORM property of the cross sections (ρ=c), structurally
independent — NOT another ORPHEUS solver, NOT MMS (MMS does not prove
RATES against an eigenvalue, same boundary as "MMS does not prove
eigenvalues"). Said so explicitly in the note. The fixed-point-
invariance-under-splitting framing cites `vv-principles` Mode 9
(synthetic-acceleration verified to SAME fixed point). The σ_r-trap
cites anti-pattern #4 (homogeneous/isotropic blind to angular
structure) — the isotropic box can't see the 46–56% error.

**Reference hygiene:** `[Pautz2002]` was NOT defined — added the
`.. [Pautz2002]` bib entry next to `[AdamsLarsen2002]` AND extended the
AdamsLarsen entry with the §II (ρ=c) + §VI (KBA parallel sweep) section
pointers the new prose cites. (Undefined-citation-warns + the
`grep '^\.\. \[Key\]'`-before-citing rule are in AGENT.md.)

**Code-source anchors used (all rich docstrings — the algebra-of-record
seeds):** `sweep_schedule.py` module docstring (Jacobi/G-S split,
B_lower/B_upper cyclic edges, mesh-time-derived); `_GaussSeidelResolvent`
+ `_select_si_resolvent` docstrings (the seed-then-overwrite loop +
resolvent/gains tuple); `scattering.py:907-922` σ_r LATENT-TRAP comment
(verbatim 46–56% / Σ_s0/σ_r≈39 numbers); error_catalog ERR-056 (the
shared-face post-mortem); the test file's PREMISE-CORRECTION docstring
(the 0.86–0.92× re-scoping + measured counts). The test docstrings are
the BEST evidence source — they carry the measured numbers + the
honest-scope correction in prose.

(Build-gate, venv/worktree, and Nexus-unavailable discipline are in
AGENT.md "Build-Gating & Cross-Ref Reality".)
