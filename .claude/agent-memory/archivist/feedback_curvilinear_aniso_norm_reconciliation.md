---
name: curvilinear-aniso-norm-reconciliation
description: W5 docs for the #229+#9+#233 curvilinear-aniso program — a multi-error norm-reconciliation close-out (one "floor" = 3 distinct errors split by volume-L2-vs-L∞), a τ-clamp mis-citation DESIGN finding, a NEW ERR-059 mint for an inherent WONTFIX limitation, and a conflated-bibliography-key flag. Sequel to the ERR-058/#195/#196 family in discrete_ordinates.rst.
metadata:
  type: feedback
---

W5 of `fix/curvilinear-aniso-pole-and-clamp`. Wrote: NEW H1 section
`sn-curvilinear-aniso-norm-reconciliation` in `discrete_ordinates.rst`
(inserted AFTER the #196 close-out, before "Krylov inner solver");
ERR-059 + τ-clamp finding + ERR-026-surviving-manifestation note in
`error_catalog.md`; Deliverable-3 doc surfaces (`reduced_operator.py`
docstring, `structured_geometry.rst`); `@catches("ERR-059")` ×4 on the
W2 gate. Baseline + final both = **1 warning** (pre-existing `paramref`).

## The load-bearing patterns (reuse these)

1. **Norm-reconciliation as the headline.** When ONE apparent "floor"
   resolves into N distinct errors because two investigations measured
   DIFFERENT norms (volume-weighted L2 `√Σ V·diff²` vs pointwise/L∞),
   LEAD with the norm difference, not the errors. The 3-error table
   (#/error/dominant-norm/quadrature-scaling/status) IS the section's
   spine. The killer fact: a pole-cell error at ONE cell of V~h³
   contributes √V~h^1.5 → ~h^2.5 to L2 → subdominant → INVISIBLE to the
   production L2/keff gates → needs an L∞/per-cell probe to surface.
   This is a "norm gotcha" Key Facts bullet.

2. **A clamp/patch mis-citation is a correctness-of-DESIGN finding, NOT
   an ERR-NNN bug instance.** It gets a labeled SUBSECTION in
   error_catalog.md (no `## ERR-NNN` header → not counted in the ERR
   total, no catching test required, audit doesn't flag it as MISSING).
   Triple-confirm structure: (1) literature (the patch is mis-cited; the
   real source recommends the UNpatched form), (2) the guarded-against
   pathology never arises on physical fields (every activation spurious),
   (3) stability holds without it. Record the "mixed accuracy signature"
   gotcha: removing a fortuitous-cancellation floor can RAISE a fine
   floor while cleaning the coarse rate — the lower floor was an
   accident, not a gain.

3. **Mint a NEW ERR-NNN for an INHERENT WONTFIX limitation** (not just
   for fixed bugs). ERR-059 = pole-cell O(h), failure mode #5
   (closure/index at a coordinate singularity), Status "DOCUMENTED
   INHERENT LIMITATION (WONTFIX-for-DD)". The tracking issue STAYS OPEN
   (tracks the future higher-order scheme). When you mint it, ADD the
   `@catches("ERR-NNN")` to the characterization gate (brief required it)
   — including the GUARANTEE tests (they pin the L2-invisibility that is
   part of the mechanism) not just the characterization tests. Then flip
   the gate's docstring "minted in W5 / @catches will be added" TODO to
   present-tense. Verify via `tests._harness.audit` exit 0 + ERR-NNN NOT
   in the MISSING list.

4. **ERR-026 "surviving manifestation" note.** When an old wrong-FP
   family is CLOSED but a DISTINCT inherent defect survives, add a
   "**Surviving manifestation (distinct root, NOT the X class)**" block
   to the old entry's final status, pointing at the new ERR + explaining
   WHY it's a different root (here: ERR-026=wrong fixed point CLOSED;
   ERR-059=spatial ORDER at the singular cell, FP now correct). Also flip
   the stale "4 xfail markers remain, tracked at #229" → "REMOVED by W3,
   retuned to assert-what-is-true".

5. **Two-unrelated-paths framing for "anisotropic".** In curvilinear SN,
   "anisotropic" = Path-(I) geometric α-dome redistribution (P0-only,
   #229) XOR Path-(II) Legendre P1+ SCATTERING (geometry-agnostic, #9).
   Document the L0 "operator-admits" trick (feed known aniso ψ, isolate
   S₁.apply−S₀.apply per-ordinate vs SH-table-indep hand-ref — NOT
   weight-summed, α-dome telescopes) + L1 directional eigenvalue (P1
   forward-peaked LOWERS keff via enhanced leakage in a VACUUM sphere;
   leakage-monotone R-small>R-large>0 is the structural neg-control). 1g
   is LEGITIMATE for an operator/flux-shape claim (the Cardinal 1g bar is
   EIGENVALUE-only).

## RST/build specifics that bit (or would have)

- **Title-underline-too-short** on em-dash (—, 1 codepoint/3 bytes) and
  `:math:`-bearing H4 titles. Compute `len(title)` in python3, NOT wc -c.
  Got 4 warnings, fixed all by exact-length `-` underlines.
- **Duplicate-citation `[Hebert2009]`** — already defined in
  collision_probability.rst. Citations resolve CROSS-DOC, so DON'T
  redefine it; put the rich §3.9.4 detail in PROSE instead and rely on
  the cross-doc bib entry. (I initially added a 2nd `[Hebert2009]` def →
  warning → removed it.) NEW keys (`[BaileyMorelChang2010]`,
  `[WuXieFischer1999]`) were fine to add (page-unique).
- **`matrix.rst` is AUTO-GENERATED** ("Do not edit by hand"). A label
  reference there (`morel-montry-clamp`) regenerates from the registry;
  since the label still exists + stays `vv-status documented`, the row is
  correct and needs NO hand edit. The brief's "update matrix.rst:741"
  was satisfied by NOT touching it (it auto-regenerates).
- **Duplicate `:label:` across .py-docstring + .rst.** `reduced_operator.py`
  is NOT automodule'd → its docstring math is NOT rendered → the
  `:label: morel-montry-clamp` there did NOT emit a warning (latent
  smell only). Collapse = make the RST the single rendered DEFINITION;
  de-LABEL the .py block (turn it into a plain math block + `:eq:` +
  `:ref:` cross-link). This is the safe way to collapse a .py↔.rst dup
  label without an automodule.

## FLAG raised (reported to main agent)

The `[Bailey2009]` bibliography entry in discrete_ordinates.rst CONFLATES
two papers: it carries the Bailey-Adams-Yang-Zika *title* (diffusion-FE,
unrelated) but attributes the curvilinear SN Eq.50/74 to it. The real
source is Bailey-Morel-Chang 2010 (Eq.43 = the exact-on-linear M-M
weight). I added a corrected `[BaileyMorelChang2010]` entry + a
`.. warning::` on `[Bailey2009]` recommending body cites migrate. Did NOT
blanket-rewrite the 6+ `[Bailey2009]_` body cites (out of W5 scope;
flagged for a follow-up). The .py docstring's References already had this
correction noted — the bib entry just hadn't caught up.

Cross-ref: [[feedback_err058_success_closeout_supersedes_phase_chain]]
(the #195 close-out this builds on), [[feedback_issue_196_eigenvalue_verification_closeout]]
(the #196 sequel), [[feedback_retirement_docs]] (the retraction-tombstone
+ flip-stale-narrative discipline used on the Phase C/D sections).
