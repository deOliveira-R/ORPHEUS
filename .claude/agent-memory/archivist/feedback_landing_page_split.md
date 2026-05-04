---
name: Theory-section discrete-vs-reference split via landing pages
description: Pattern for restructuring a flat theory toctree into a tree of landing pages without `git mv`-ing files
type: feedback
---

When the theory/ section grows beyond ~10 pages and starts mixing
fundamentally different artefacts (production solver theory vs
verification reference theory), the cheap restructuring move is
**landing pages, NOT folder reorganization**. The 2026-05-03 split
of docs/theory/ into infrastructure / discrete / reference branches
(commit `e90029a`) used this pattern.

**Why:** A flat 14-page Theory toctree was conflating two audiences
(end users hacking the SN solver vs verification engineers writing
reference solutions). The user's brief proposed `git mv` into
`theory/discrete/`, `theory/reference/`, `theory/infrastructure/`
subfolders. But moving files at this scale risks breaking 49 active
`:doc:` references across the corpus, perturbs Nexus indexing, and
forces a one-shot atomic restructure that's hard to validate.

**How to apply:**
1. Inventory `:doc:` references first (`grep -rn ":doc:" docs/`).
   The decision metric: if move count × avg-cross-refs is high,
   landing pages win.
2. Add a top-level `theory/index.rst` that toctree-includes the
   sub-landing pages (one per branch).
3. Each sub-landing page (`theory/transport_methods.rst` and
   `theory/reference_solvers.rst` in this case) has its own toctree
   listing the family's pages.
4. Top-level `docs/index.rst` collapses to ONE entry: `theory/index`.
5. Existing `:doc:\`peierls_nystrom\`` references continue resolving
   without modification — paths haven't changed.

**The landing page IS the architectural document** — it teaches
the reader the taxonomy, links to all members of the branch, and
defines what each branch's V&V claim is. The `reference_solvers.rst`
landing page in this commit is the canonical instance: it documents
the three-meanings Green's-function taxonomy + per-solver pillar
classification table that previously lived only in `.claude/scratch/`
literature memos.

**Reserve folder stubs** are cheap insurance: when adding a new
empty package (`spectral_resolvent/`, `pn_method/`, etc.), drop a
stub theory page with just Key Facts + canonical references +
`:label:`. Future code lifting into that folder gains a working
cross-reference graph for free.

**Build constraint:** Sphinx -W must remain clean across the
restructure. The pre-existing Nexus "verifies has no matching
equation" info-lines (~17 in this project) are NOT warnings (they're
info-level Nexus extension output) — count them in the baseline and
ensure post-edit count is unchanged.
