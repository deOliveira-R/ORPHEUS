---
name: capstone-completion-status-reaudit
description: Recipe for the COMPLETION (final) phase of a multi-phase campaign's capstone/overview page — the first pass is re-auditing every ship-state STATUS claim against what merged since the page was written, BEFORE adding new content; plus documenting a designed-but-UNBUILT sibling type as a seam (literal not :class:, tied to the real guard) and the success-capstone "documented-future seam" honest-status block. Distinct from [[feedback-capstone-architecture-page]] (NEW page for a layer) and [[feedback-capstone-root-cause-ruling]] (retrofit a theorem).
metadata:
  type: feedback
---

# Capstone-completion = STATUS RE-AUDIT first, new content second

When the task is the **completion / final phase** of a multi-phase
campaign's **capstone** (the overview page that ties the incrementally-
documented pieces together — e.g. the frame-projection campaign's P7
on `galerkin_projection.rst`), the load-bearing work is NOT the new
prose. It is **re-auditing every ship-state STATUS claim** the page
made when it was written at an EARLY phase, because those claims stale
SILENTLY as later phases land concrete consumers and the prose still
"reads fine".

**Why:** a capstone written at phase N says things like "X is the
single concrete frame shipping today", "no concrete Y ships yet", "no Y
numerical evidence ships yet", "Y will inherit the tests when its
test_basis lands", and a consumer-table **Live/Pending column**. When
N+1/N+2 ship the concrete Y (here: forward `Solution.homogenize` (P3) /
`Solution.condense` (P5) shipped as concrete `PetrovGalerkinFrame`
consumers), EVERY such claim is now false — but nothing warns (prose
staleness is invisible to `-W`; L-002). The reader trusts the capstone
as the navigational layer, so a stale "Pending" there mis-routes a
fresh session worse than a stale leaf page.

**How to apply (the re-audit pass, FIRST — before drafting new
sections):**

1. **Grep the live code for each consumer verb the page frames as
   future** (`grep "def homogenize\|def condense"`). If it exists,
   every "ships today is only X" / "no concrete Y yet" / "Y will
   inherit tests" claim is stale → flip to "Y shipped (Pn)", cross-ref
   the consumer's derivation + its L0 gates.
2. **Walk the consumer table row by row.** The Live/Pending column is
   the densest staleness. Split a stale "Y — Pending" into the SHIPPED
   forward row (Live, Pn) and the still-pending NON-degenerate row
   (the adjoint-weighted / higher-rank case — Pending, blocked on Z).
3. **Audit the Key Facts + the closing of each discipline section +
   the Numerical-evidence section** — these carry the "single concrete
   / no evidence yet" phrasings. Re-grep them specifically.
4. Only THEN add the new synthesis content. The new content (a
   composed-verbs section, a missing primitive) is the EASY half; the
   stale-status correction is the half that actually completes the
   capstone and the one a "just add a section" reading misses.

# A designed-but-UNBUILT sibling type is a SEAM, not a class

The brief named a hierarchy `{GalerkinFrame, LeastSquaresFrame}` — but
only `GalerkinFrame` is a class (`grep __all__`); the "least-squares"
case is the `GramStructure.DENSE` → `MissingCapability` REFUSAL in
`FrameBase.project`/`.gram` (gated to an issue). A sibling doc had a
**dangling `:class:LeastSquaresFrame`** (renders plain-text, no `-W`
warning — L-002) pointing at the capstone, which never documented it.

**Resolution (trust the code's real seam over the brief's class name):**
- Document the unbuilt sibling as a **named seam**: literal
  ````LeastSquaresFrame```` (NOT `:class:`, which falsely implies an
  importable class), in a labelled subsection, marked "designed but not
  built".
- **Tie it to the ACTUAL code mechanism** — the enum value, the guard
  that raises, the deferral issue — so it's grounded, not vapor. The
  3-way `GramStructure` gate (DIAGONAL / PARTITION_OF_UNITY → built
  row-sum probe; DENSE → refused dense solve) IS the seam's real shape;
  a `.. list-table::` of "structure × normalisation × built?" carries it.
- **Repoint the dangling ref** to the new anchor and change its
  `:class:` → literal. Verify in built HTML: 0 `py-class` xref spans
  for the name, ≥1 plain-literal span (the L-002 distinction that
  proves the dangling ref is gone, not just relocated).

# The success-capstone "documented-future seam" honest-status block

The success analogue of the close-out's status banner (AGENT.md
Close-Out arc). For a theory that is fully derived but NOT implemented
(here P6 adjoint-weighted homogenization, blocked on the adjoint flux
φ* from a sibling campaign):

- A dedicated anchored section + an `.. important:: Status — theory
  documented, implementation NOT built (Pn)` block stating: (a) what is
  NOT wired (φ* not in `homogenize`/`condense`; both run the forward
  φ*=φ degenerate), (b) what it's blocked on (the upstream campaign +
  the exact operator it needs — `(L+C−S)ᵀψ*=q*`), (c) a sentinel: "do
  **not** read the forward homogenisation's green rate gates as evidence
  for the adjoint-weighted case" (the L-010 vocabulary: green gates pin
  the φ*=φ degenerate, not the lift).
- **Frame it as a one-line `test_basis` swap** (φ→φ*) that the
  type-carries-the-discipline architecture buys — that's the payoff of
  the whole campaign, so the seam doubles as the architecture's
  vindication, not just a TODO.
- **Generalises a shipped degenerate:** the eigenvalue-consistent
  bilinear ⟨φ*,Σφ⟩ generalises the φ†=1 degenerate
  `IntegratedReactionRate(Σx)=∫⟨Σx,φ⟩dV` (already minted) — name the
  degenerate so the lift is "replace the implicit φ†=1 with a real
  adjoint", concrete not abstract.
- **Sharpen sibling prose that says the future case "is built in a
  later phase"** (implies scheduled) → "deferred, blocked on X"
  (honest). **Bidirectional wiring:** capstone seam ↔ the derivation
  page's bilinear `:eq:` label + its vv-status rationale (sharpen the
  comment too — comments take no `:ref:`, prose does).
- The bilinear label stays `.. vv-status: <label> documented` (L-004,
  literature-transcribed / representational — NOT a solver claim); do
  NOT claim it tested.

**Quality self-assessment (Directive 3):** Derivation-depth 4 (added
the adjoint Σ_R=∫φ*Σφ/∫φ*φ + the φ†=1-degenerate generalisation + the
3 composed-verb compositions, but CROSS-REF'd the full bilinear
derivation to `discrete_ordinates` rather than re-deriving — correct
capstone discipline, synthesis not duplication); Cross-refs 5 (3 new
anchors + 6 cross-doc refs + the bidirectional seam wiring, every one
grep-gated AND built-HTML-href-verified; the dangling
`:class:LeastSquaresFrame` killed, confirmed 0 py-class spans);
Numerical-evidence 3 (pointed at the shipped forward L0 gates
`tests.sn.test_homogenization`; no NEW tables — a synthesis capstone
adds none, structurally absent per the rubric note, not a deficit —
WEAKEST); Failed-approaches 5 (preserved the #268 two-reversal history,
the metric-fold-breaks-under-adjoint argument, the
LeastSquaresFrame-not-needed-at-rank-0 note); Code-traceability 5
(every verb/class tied to live code — `FrameBase.conjugate/project/gram`,
`GramStructure`, `MissingCapability`, `Solution.homogenize/condense`,
`IntegratedReactionRate`, all grep-confirmed); Derivation-source 4 (math
from the live `frame.py`/`basis/base.py` docstrings + the
`discrete_ordinates` derivation, read first per L-005; the bilinear is
Hébert §6/§13 literature-transcribed — no SymPy script applies to a
projection identity).

See [[lessons]] L-001 (the brief's hierarchy `{Galerkin, LeastSquares}`
was the DESIGN, only Galerkin built — verify against `__all__`, trust
code), L-002 (dangling `:class:` renders plain-text; built-HTML
py-class-span audit is the real gate), L-004 (bilinear label stays
vv-status:documented), L-007 (preserve the #268 reversal history /
honest seam), L-010 (forward green gates are NOT evidence for the
adjoint lift). Siblings: [[feedback-capstone-architecture-page]],
[[feedback-capstone-root-cause-ruling]],
[[petrov-galerkin-homogenization-reframe]].
