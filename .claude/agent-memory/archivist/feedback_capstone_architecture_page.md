---
name: capstone-architecture-page
description: Writing a NEW capstone theory page for a multi-phase refactor that COMPLETED (a layer ABOVE existing per-method pages) — section shape, cross-ref-not-duplicate discipline, measurement+decision section, verbatim user-decision admonition
metadata:
  type: feedback
---

Writing the S5.5 capstone for the SN loss-representation architecture
(issue #222 sweep-strategy carve + S6 re-layering, COMPLETE). NEW page
`docs/theory/loss_representations.rst`. Sibling of
[[feedback_post_wave_cleanup_docs]] (close-out arc) but distinct: this
is a CAPSTONE for a FINISHED multi-phase campaign, not a single
follow-up close-out.

**Rule: a capstone documents the LAYER, not the methods.** The page sits
ABOVE existing per-method pages. discrete_ordinates already owned the
WDD cell algebra + cumprod recurrence + anti-diagonal wavefront + the
S6.4(e) 3-layer walk/level-op/kernel stack (`sweep-dispatch-relayering`);
operator_algebra owned the Wave-O algebra + matvec≡sweep + cochain. The
capstone's UNIQUE content = the REPRESENTATION layer: the N schedules as
ONE lower-triangular operator, the selection SSOT (supports/default_for/
registry), the one-walk + one-instance type-fact theorems, the
measurement+decision. CROSS-REFERENCE everything else via :ref:/:eq:/
:doc: — NEVER re-derive. A capstone that re-derives the cell math is a
duplicate, not a capstone.

**Section shape that worked (9 H1s):** Key Facts (admonition) → native
mathematical frame (the unifying abstraction — here: lower-triangular,
SOLVE=fwd-subst, APPLY=row-action, Resolution-A glue) → the N variants
(list-table + one H2 per variant) → selection (the SSOT, 3-consumer
table) → the structural theorems (one-walk/one-instance, each with its
discriminating spy + WHY-it's-not-a-tautology) → the
measurement+decision → verification architecture (bit-id vs principled-
equiv, the load-bearing taxonomy) → history/rationale (carve arc table +
rejected alternatives + succession + deferred extension point) →
Literature → See also.

**The measurement+decision section is the WHY of the default** and the
highest-value content (Cardinal Rule 3). Shape: the benchmark table
verbatim (cite the diagnostic script path + the #222 comment id) →
end-to-end number → 3 NUMBERED findings (each explains a column, e.g.
"no memory edge at d=2", "advantage narrows with angular order", "end-
to-end dilution") → the decision in a VERBATIM `.. admonition::
(verbatim) :class: important` block (the user's exact words are
archival-load-bearing — they justify KEEPING the non-default peer) →
the mechanical consequence (one-line registry reorder; what stayed
unchanged = blast-radius).

**Verification section maps 1:1 to vv-principles.** Two granularities:
ACROSS schedules = principled-equivalence (nulp/abs-tol oracle + solver-
tol Mode-9 FP-invariance; cite G4.a/G4.b/G6 + the ERR-056 d=2 honest-
limitation + the context-manager-not-fixture forcing + the Fork-B2
polarity inversion); WITHIN a schedule = bit-identity (window≡full
array_equal + the kernel sha256 source-of-record + the affine converged-
bytes golden with the regenerate-in-commit discipline). Add the MMS-
reaches-flux-shape-not-eigenvalue note (closed-form k_inf is the
structural-independence anchor; SI≡Krylov is necessary-not-sufficient
twin agreement). Read the ACTUAL test docstrings — they carry the
nuance (e.g. scan_march_equivalence uses abs-tol NOT nulp because near-
zero shed amplifies ULP; the end-to-end test inverted forcing polarity
at the flip).

**Production-code staleness is OUT OF SCOPE.** The module docstring's
phase-status block + the registry docstring were stale future-tense
("S2: MATVEC side", "ScanMarch is OPT-IN Fork B1 default unchanged") —
the code had advanced past them. Document the FINAL state from the live
class bodies + registry ORDER + plan STATUS history; do NOT transcribe
the stale prose, do NOT touch production code on a docs task. (The
`_sweep_1d_cumprod` ref in discrete_ordinates:1368 is also stale —
dissolved — but not my file to fix either.)

**Cross-ref grep-gate BEFORE building** (the :ref:/:eq: do warn intra-
doc; code-xrefs render plain-text silently in this non-nitpicky repo).
Three greps: (1) every :ref: I cite → `^\.\. _NAME:` exists; (2) every
:eq: I cite → `:label: NAME` exists; (3) my OWN new :label:s are unique.
The `\b` boundary in grep matches `loss-rep-scanmarch` AND
`loss-rep-scanmarch-solve` — count DEFINITIONS not boundary-matches to
confirm uniqueness. Code-xrefs verified by having read the live source
this session (ground truth), NOT by the warning count.

**Baseline this session = 11 (NOT the 1 in AGENT.md).** Worktree
`sn-nd-layout` baseline: 1 paramref ERROR + 2 homogeneous plot WARNINGs
+ 8 verification.rst CRITICAL include-path = 11, EXIT=0. The AGENT.md
"baseline is 1 warning" is the MAIN-checkout figure; a worktree with
un-built generated includes carries more. Always establish the baseline
with a `-E` build in THE WORKTREE before editing; gate on count-AND-set
unchanged, not on AGENT.md's number. Grep WARNING:|ERROR:|CRITICAL:
(three severities) + the summary line + EXIT.

Quality self-assessment (1-5): Derivation depth 4 (frame + Resolution-A
+ scanmarch coeffs derived; cell math cross-ref'd not re-derived, correct
for a capstone); Cross-references 5 (every class/func/meth/mod linked,
all :ref:/:eq: grep-verified); Numerical evidence 5 (the full Fork-B2
benchmark table + end-to-end + the 3 findings); Failed approaches 5
(rejected-alternatives subsection + the carve-arc table + succession +
deferred ExplicitMatrix); Code traceability 5 (every equation tied to a
:func:/:meth:); Derivation source 3 (no SymPy derivations/ script — the
math is discretization/scheduling, sourced from the diamond/scan
docstrings + Blelloch/L&M/Hébert literature, appropriate here).
