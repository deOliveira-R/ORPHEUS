---
name: canonical-axis-convention-ssot-section
description: Recipe for a single-source-of-truth convention section documenting a just-landed boundary-normalisation flip (an axis ordering enforced once at data ingest), plus the orthogonal-axis terminology disambiguation that must ride with it.
metadata:
  type: feedback
---

When a refactor lands a **canonical axis-ordering convention enforced
once at a foreign-data ingest boundary** (here: energy group 0 = FASTEST,
`eg` descending, `SigS` downscatter upper-triangular — normalised from
NJOY's opposite order by `_to_canonical_group_order` in `gendf.py`), the
documentation deliverable is a **single-source-of-truth (SSOT)
convention section** on the data-foundation page, plus a set of
`:ref:`-back pointers from every page that already states the target
convention. The shape that worked (`cross_section_data.rst`, branch
`refactor/group-0-fastest-convention`):

**The SSOT section (4 sub-parts, in order):**
1. **`.. note::` banner declaring it the single source of truth** — "every
   page defers here". Mint a clean linkable label (`.. _canonical-group-convention:`)
   so other pages `:ref:` it. Grep `_<label>:` repo-wide first (it was free).
2. **The convention stated precisely** as a bullet list — index meaning,
   array-ordering (with the interval each index spans), and the matrix
   triangle structure (downscatter upper / upscatter thermal-only-lower).
3. **An `.. important::` asserting ORTHOGONALITY to the sibling axis.**
   THE TRAP: a "convention" page for a DIFFERENT axis already exists
   (`index_convention.rst` = the array STORAGE LAYOUT `(N,ng,nx,ny)` =
   axis ORDERING within arrays). That is NOT the same as the energy-index
   ENERGY ordering. State the distinction explicitly ("which physical
   group each index labels" vs "which array axis the group sits on") and
   `:ref:` the sibling for the orthogonal concern — do NOT let a reader
   conflate them, and do NOT over-claim the sibling page documents
   something it doesn't (I first wrote ":ref: the SigS convention note in
   the SN index page" — but that page covers layout-transpose, NOT the
   `SigS^T` in-scatter transpose; reworded to inline the `Σ_s^T φ` fact
   and point at the sibling only for the layout-axis distinction).
4. **Enforcement sub-section** naming the single boundary normalisation
   (`:func:` the ingest reverser), enumerating exactly which arrays flip
   (and which do NOT — the non-energy-indexed `sig0`), and a downstream
   "order-transparent" statement (caches written post-flip).
5. **Rationale sub-section** — the WHY, ordered by weight. Here:
   (1) convergence (the flipped source was the LONE outlier; everything
   else was already canonical — name them), (2) physics-identical (the
   in-scatter is a FULL CONTRACTION, no energy-sweep ⟹ permutation-
   invariant; back it with the real L2 permutation gate — keff invariant,
   flux reversed), (3) the natural sweep direction. The contraction
   equation gets a `:label:` + `.. vv-status: <label> documented`
   (representational group-coupling definition, NOT a solver claim;
   rationale comment names the permutation gate as the verifiable
   content — L-004).

**The orthogonal-axis TERMINOLOGY note (rides WITH the SSOT, same page):**
A landed convention about axis A routinely collides with an overloaded
word that ALSO names a concept on axis B. Here "group" = energy bin
(axis A) vs **octant group** = a set of ordinate directions sharing a
sweep direction (axis B, the SN angular sweep). Document the distinction
as its own labelled section so neither is ever read as the other:
- a 2-row `.. list-table::` (axis | "group" means | where it appears);
- ground the angular sense against the LIVE algebra-of-record docstring
  (`_GaussSeidelResolvent` / `_select_si_resolvent` in `sn/solver.py`):
  it's an **octant-group/BOUNDARY** Gauss-Seidel folding the reflective
  `B` into a multi-D wavefront sweep, FP-identical to Jacobi (vv Mode 9),
  NOT a within-group scattering fold and NOT an energy sweep;
- a `.. note::` stating the absent concept EXPLICITLY ("an energy-group
  Gauss-Seidel — a downscatter cascade — does NOT exist in ORPHEUS
  today") and tying it back to rationale-point-2 (this absence is WHY the
  ordering is free to be a pure relabel). `:eq:` the contraction label.

**Extraction-internals vs stored-data index trap (deliverable 3 class):**
a constant that gates extraction (`_IG_THRESH = 95`, the thermal cutoff)
may index the **NATIVE foreign order** because it runs BEFORE the ingest
flip (verify the call order in the live code: `_init_scattering` inside
`_build_isotope`, THEN `_to_canonical_group_order` in `convert_gxs`). The
fix is NOT to renumber it — it's to ADD an `.. important::` flagging
"this boundary is in the native index; after the flip it's at the
high-index end ≈ G−95" and tag every `g ≤ 95` / `g > 95` in the section
"(native NJOY index)". Anchor the section so the SSOT can `:ref:` it.

**The flag-don't-fix satisfied-claims sweep (deliverable 5 class):**
the ~8 pages that ALREADY stated the target convention are now SATISFIED
— add a `:ref:` to the SSOT at each (the dedicated "ORPHEUS vs Sood
convention" sections and the "Indexing:" blocks are the natural spots),
and FLAG (don't rewrite) subtle mismatches: `homogeneous.rst` uses
1-based math labels (group 1 = fast) while the code index is 0-based —
resolved with a one-clause "the code index is 0-based; the 1-based
labels below are a presentation choice", NOT a renumber. A page calling
itself "project-wide single source of truth" for the same convention
(`peierls_nystrom.rst` §sig_s) gets a forward `:ref:` subordinating it to
the new SSOT, heading-rename left as a flagged larger call.

**Verification:** baseline `-E -W` (was EXIT=0, zero W/E/C) → all 7
cross-doc `:ref:`s resolved to real `href`s in HTML (grep
`href=".*#<label>"` per consuming page — the `-W`-blind plain-text check,
L-002); intra-doc `:ref:`/`:eq:` + the cross-doc sibling ref all resolved;
the SN code-xrefs render plain-text by page convention (OK — module paths
verified accurate against live `grep "class …"`). The auto-regen V&V
matrix (L-008, never hand-edit) correctly filed the new label under
documented-only AND surfaced the two real backing tests
(`test_gendf_canonical_order` foundation, `test_group_permutation_invariance`
L2) — read the L2 test's docstring to confirm the rationale's "keff
invariant, flux reversed" claim matches verbatim and is NOT over-claimed
as solver correctness (it's an L2 SYMMETRY claim).

Quality self-assessment (Directive 3): Derivation depth 4 (the convention
+ contraction + enforcement are fully stated; no SymPy derivation needed —
it's a convention, not a derived result), Cross-references 5 (every page
linked both ways, sibling-axis orthogonality made explicit), Numerical
evidence 4 (the permutation gate cited with its exact assertions; no
convergence table because nothing converges differently — a flip is a
relabel, so "numerical evidence" is structurally the invariance gate, not
a table — same structural-absence note as L-010's rubric tail), Failed
approaches N/A (convention is settled, brief said don't re-litigate),
Code traceability 5 (every claim `:func:`/`:ref:`-linked to live code +
the ingest reverser + the backing tests), Derivation source N/A.
