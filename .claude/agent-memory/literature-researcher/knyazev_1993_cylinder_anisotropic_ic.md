---
name: Knyazev 1993 (NOT Knyazev-Selivanov 2014) — cylinder linearly-anisotropic kernel
description: The canonical "Knyazev cylinder Ki_{2+k}" reference is single-author A.P. Knyazev 1993 in Atomic Energy. "Knyazev-Selivanov 2014" was a phantom citation hallucinated in a wish-list memo and propagated forward.
type: reference
---

## Verified citation

**A. P. Knyazev (1993).** "Solution of the transport equation in
integral form in a one-dimensional cylindrical geometry with linearly
anisotropic scattering." *Atomic Energy* (Plenum/Springer translation
of *Atomnaya Énergiya*), **Vol. 74, No. 5, pp. 368-374**, May 1993.
DOI: `10.1007/BF00844623`.

Russian original: *Atomnaya Énergiya* 74(5), 403-410, May 1993.
Author at Institute of Technical and Experimental Physics. Sole author.
Paywalled (Springer). Not OA. No arXiv / HAL / Zenodo deposit.

## What the paper provides

Closed-form 1-D quadrature evaluator for the homogeneous
linearly-anisotropic-scattering cylinder kernel, in Bickley-Naylor
functions of order 2+k. Structurally the cylinder analog of Sanchez
1986 sphere Eq. (A6). The Ki_{3+k_m+k_n} expansion that Phase 4 cyl
matrix-Galerkin and Issue #132 cylinder Hébert closure both depend on
comes from this paper.

## The "Knyazev-Selivanov 2014" phantom

The following claims do NOT correspond to any catalogued publication:

- "Knyazev-Selivanov 2014 JCP, stable cylindrical Ki_n series" —
  no JCP paper, no Knyazev-Selivanov co-authorship anywhere in
  CrossRef / OpenAlex / Semantic Scholar / OSTI / web search.
- "Knyazev, B. A. & Selivanov, A. N. (2014). *Bickley-Naylor functions
  and their integrals*, Nauka, Moscow / English translation" —
  no such book in WorldCat / Russian State Library / shadow indices.

Both phantoms originated in
`peierls_greens_extensions_lit.md` (literature-researcher wish-list)
and propagated to `sanchez_chandrasekhar_three_meanings.md` and
`/workspaces/ORPHEUS/.claude/scratch/sanchez_chandrasekhar_gap.md`.
The author surname "Knyazev" + a guessed year/co-author/journal was
treated by downstream sub-agents as a real reference. It is not.

## How to apply

When the user (or a sub-agent's brief) cites "Knyazev-Selivanov 2014"
for cylinder Bickley-Naylor / Ki_{2+k} / 3-D cylinder collision
probability machinery, redirect to the actual 1993 single-author
paper. Use DOI `10.1007/BF00844623`. The 1993 paper is paywalled but
accessible via Springer institutional subscription; the Russian
original is in *Atomnaya Énergiya* 74(5), 403-410.

Search queries that surface this paper reliably:
- CrossRef `journal_search("1063-4258", query="Knyazev cylindrical")`
- CrossRef `get_work("10.1007/BF00844623")`

Search queries that DO NOT find it (and indicate paywall depth):
- Semantic Scholar by DOI returned title only, no abstract.
- OpenAlex returned `cited_by_count: 0` (the paper is invisible to
  the open citation graph — typical for early-1990s Russian-translation
  journal papers).

## Cross-references in tree

- `/workspaces/ORPHEUS/.claude/agent-memory/numerics-investigator/issue_132_cylinder_hebert.md` —
  cites the 1993 paper correctly with the right DOI. **Use this as
  the trusted form.**
- `/workspaces/ORPHEUS/.claude/scratch/knyazev_selivanov_lookup.md` —
  full lookup brief written 2026-05-03 with negative-finding evidence
  for the 2014 phantom.
