---
name: solver-replacement-campaign-closeout
description: Recipe for the P8 docs close-out when a whole solver is REPLACED (a legacy island → the operator-algebra family) across a theory page + verification page + sibling forward-ref page, incl. a mis-named-law re-attribution and a brief-vs-live registry-mechanism catch.
metadata:
  type: feedback
---

The P8 documentation phase of a campaign that REPLACES an entire legacy
solver (a MATLAB-port "island" → the shared operator-algebra family) is
a 3-surface job with a distinctive shape. Instance: #290 diffusion
integration (`diffusion_1d.rst` overhaul + `verification.rst`
retired-tier note + `operator_algebra.rst` DSA forward-ref
reconciliation). The moves, reusable:

**1. Theory-page overhaul = Key-Facts operator-family reading + a NEW
"Operator-family architecture (production)" §, keeping the continuous
references intact.** The island page carried classical-form + continuous
references only. Add the production layer WITHOUT rewriting the
derivations: (a) Key Facts gains the `(L+C−S−B)ψ = Fψ/k` reading (new
labeled eq + removal-as-C−S-theorem), (b) a new h1 § documents the
scalar composite / L,B leaves / shared C,S,F arms / resolvent-choice
rationale / mesh+Protocol layering / DSA seam / object-level gates, (c)
the continuous derivations + MMS § stay byte-untouched (they are the
references, not the production path). The rejected-alternative
(Green-splitting resolvent diverges for elliptic — ρ≥1, not the
nilpotent ρ=0 of a triangular sweep) is load-bearing Sphinx-as-brain
content — a `.. warning::` block, stated as the campaign's empirical
finding, NOT an over-claimed ρ→1-therefore-diverges mechanism (ρ→1 is
slow convergence, not divergence — Cardinal Rule 1 on the mechanism).

**2. Mis-named-law re-attribution: rename the LAW word across the page,
NEVER the math or the labels.** #290 ruling 3: the island called a
hard-Dirichlet φ=0 wall "vacuum"; the faithful vacuum is Marshak J⁻=0,
and φ=0 is its own honestly-named `zero_flux` law. The analytic sine
references satisfy φ=0 ⟹ they are ZERO-FLUX references. Sweep every
"vacuum" prose mention (incl. the section heading + resize its
underline + the `(vacuum left/right, N equations)` annotations + the
investigation-history dead-ends) to "zero-flux", KEEP every eq `:label:`
and all math verbatim, and add ONE `.. important::` block stating the
distinction + that the faithful vacuum has property-gates today + the
analytic Robin reference is a follow-up (#293). A renamed helper
(`_solve_2region_vacuum`→`_zero_flux_eigenvalue`) is a dead `:func:` —
repoint (grep gate; it is private ⟹ plain-text anyway).

**3. Investigation-history reframe: split reference-solution dead ends
(LIVE) from the retired-island solver dead end (MOOT).** A page's
"abandoned approaches" § usually mixes two machines. Dead ends about the
CONTINUOUS REFERENCE (transfer-matrix / mode-basis) are LIVE — keep them
present-tense as "why the reference is built this way". The dead end
about the RETIRED solver (hardcoded outer-tol / BiCGSTAB) is MOOT for
the modern exact-LU path — de-role its `:class:`/:func:` to LITERALS
(``DiffusionSolver`` names "refer to the island, now deleted"), past-
tense it, and state which conclusion transfers (the general
check-your-tolerances lesson) vs which is moot (no inner iteration
exists to mis-tolerance). Fix the pre-existing "two dead ends" count
+ wrong Cardinal-Rule cite while you are in the paragraph.

**4. Sibling forward-ref page: flip "seam expected" → "seam now real",
keep the still-unbuilt consumer unbuilt.** `operator_algebra.rst` had
DSA passages: "A diffusion operator `A_diff = L+C−S` was EXPECTED".
Update the operator to the landed spelling (`L+C−S−B`), add
`:doc:`/theory/X`` + `:mod:` cross-refs, and note the operator now
EXISTS (#290 P4) — but the CONSUMER (DSA #2, `as_dsa_source`) is still
unbuilt. Do NOT over-flip: "the accelerator itself is still unbuilt;
its consumable operator now exists." Same for the "no-(A,S,F)-triple"
resolvent taxonomy: diffusion NOW has the L/C/S/F family but still
belongs to the monolithic-resolvent camp because it has NO SWEEP (its
`A⁻¹` is a direct LU of the fused A, not an `(A−S)⁻¹` iterated over a
sweep) — reword "no-triple" → "monolithic-resolvent (no-sweep)", not a
blanket flip.

**5. The brief can mis-state the CURRENT mechanism — verify the registry
against live code (L-001 face).** The brief said "the registry is now
all-analytical/eager" for `verification.rst`. LIVE `reference_values.py`
still has a lazy tier — but a DIFFERENT one: `_CONTINUOUS_BUILDERS`
name→thunk (Issue #212, for O(minutes) Peierls references), NOT the
retired Richardson `_LAZY_LOADERS`/production-solver side-channel. Fold
the RETIRED mechanism (Richardson + `_LAZY_LOADERS` + circular-import
rationale) into a `.. note:: Historical`, and describe the CURRENT
tier (build-cost laziness) accurately. Read the retrieval fn
(`continuous_get`, `_build_continuous_registry`) — do not transcribe
the brief's characterization.

**6. Development-history § = the discrete_ordinates exemplar, per-phase
rows.** Newest-first `.. list-table::` (Phase | milestone what+why |
commit-hash), one row per phase + a `pre-#290` island-era row. Unmerged
branch: state "the branch is unmerged at the time of writing (P8)".

Gates that MATTER here beyond the standard `-E -W`/audit: (a) the
`:noindex:`-whole-package plain-text discovery ([[lessons]] L-002
sharpening — import-verify every symbol, the api page governs the
link); (b) preserve the two verifies-target labels (`diffusion-coefficient`
in `tests/data/`, `diffusion-mms` in `tests/diffusion/`) — grep tree-wide
for verifies(), not just the module's own test dir. See also
[[feedback-capstone-root-cause-ruling]] (structural-WHY retrofit) and the
AGENT.md Close-Out Narrative Arc (the FALSIFICATION variant; this one is
a SUCCESS-story replacement, so no tombstones — just island→family
tense-flips).
