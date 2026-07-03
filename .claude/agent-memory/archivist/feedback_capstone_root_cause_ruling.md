---
name: capstone-root-cause-ruling
description: Recipe for archiving a CAPSTONE architectural RULING that supplies the structural WHY (a theorem) behind a design split the existing pages only ASSERTED axis-by-axis — distinct from [[feedback-capstone-architecture-page]] (a NEW page for a LAYER) and from the close-out arc (a falsification). Instance: the Funk-Hecke/Schur eigenbasis-ownership ruling that unifies the Galerkin-vs-Petrov-Galerkin split.
metadata:
  type: feedback
---

# Capstone root-cause ruling — archiving the WHY behind an already-asserted split

The rule: when a literature + structural investigation produces a
CAPSTONE result that explains the **root cause** of a design decision
the docs already *state* (here: "angular scattering is Galerkin,
energy/spatial are Petrov-Galerkin" — asserted axis-by-axis on
`galerkin_projection.rst` before the WHY existed), the archival move
is to **retrofit the theorem into the existing pages**, NOT to write a
new page. The ruling is a *derivation of a previously-asserted fact*,
so it lands as a new SECTION in the page that owns the asserted fact,
plus thin cross-referenced echoes in the sibling pages.

**Why:** the value of a root-cause ruling is that it converts a set of
axis-by-axis conventions into ONE structural principle ("an operator
owns its frame iff the frame is its eigenbasis"). A future session
that reads the principle can re-derive every axis's discipline instead
of memorising a table. The principle is the campaign-wide prize — give
it its own subsection with a `.. list-table::` reading it across every
axis, and state it as a biconditional in bold.

**How to apply** — the 6-move shape (instance: the eigenbasis ruling,
`galerkin_projection.rst` §`frame-eigenbasis-ownership`, #268 family):

1. **Lead with the principle as a bold biconditional**, before any
   proof. ("An operator owns its frame **iff** the frame is its
   eigenbasis.") The reader must grasp the unifying claim before the
   machinery.
2. **The structural leg = the theorem written out.** When the
   production code already has the operator factorisation (here
   `S = R∘Λ∘M`), SHOW it is a named theorem (spectral theorem
   `A=UΣU*`: M=U*, Λ=Σ, R=U). The existing addition-theorem/kernel
   identities become the *spectral resolution* — re-frame, don't
   re-derive. Cite the theorem source (Funk–Hecke → Müller 1966;
   Schur's lemma for the block structure + the irrep-dimension
   weights).
3. **The asymmetry that fixes ownership.** A "X owns the frame" claim
   is only load-bearing if a SIBLING operator does NOT get the same
   basis. Document the contrast in a `.. list-table::` (scattering
   diagonalised / streaming block-tridiagonal via Clebsch–Gordan).
   Without the asymmetry the ownership is a coin toss.
4. **Literature corroboration as a NEGATIVE-SPACE table.** The
   strongest evidence for "M exists only because of scattering" is a
   paper that REMOVES M (Ahrens 2014 LDO). Build a `.. list-table::`
   citing every reference with eq/section number; state explicitly
   "N references, ZERO cross-validation against any non-scattering
   use" (the success-story analogue of the close-out's
   "ZERO cross-validation" line).
5. **The unifying principle subsection** = the campaign prize. One
   `.. list-table::` (axis × symmetry × eigenbasis? × discipline),
   then the forced consequence (no symmetry ⟹ no eigenbasis ⟹
   solution-weighted ⟹ Petrov-Galerkin on the test side). Include the
   negative consequences too (here: "fission does not own an angular
   frame" — its axis is energy, no eigenbasis).
6. **The relocation tripwire as a forward-looking `.. _label:`
   section.** A constructor-ownership claim ("scattering owns the
   frame") is true UNTIL a 2nd consumer with an independent parameter
   arrives. Document the exact flip condition + that the neutral
   factory already exists and anticipates it. Give it its own anchor
   so a future session landing on the 2nd-consumer task finds it.

**vv-status discipline (L-004):** every eq-label minted for a theorem
transcription (zonal-kernel definition, Funk–Hecke eigenvalue formula,
spectral-theorem factorisation, Pℓ recurrence) is `.. vv-status:
<label> documented` with a rationale comment naming (a) the literature
source and (b) the bit-identity gate that pins the *implementing* code
(here the 0-ULP windowed-vs-full crosscheck). These are structural
transcriptions, NOT solver claims — they land in "Documented-only".

**Cross-ref-not-duplicate (the sibling echoes):** the SH page and the
operator-algebra page each get a SHORT section/paragraph stating the
eigenbasis fact in their own vocabulary + a `:ref:` to the deep
treatment. Do NOT duplicate the proof. The SH page answers "why these
functions"; the operator-algebra page re-reads its `frame.conjugate`
2-cell as spectral; both point to the one full derivation.

**Build traps hit this session (all mine, all fixed pre-report):**
- The em-dash `—` in a section title is 1 code point; my underline
  came out 1 short → "Title underline too short". Size with
  `len(title)` in python (L-009).
- `**NN**(4)` in a list-table cell → "Inline strong start-string
  without end-string": a closing `**` immediately followed by `(` is
  illegal (an *opening* delimiter is not in the allowed
  after-closing-`**` set; `:` IS allowed, which is why `**84**:33`
  did not warn). Put a space: `**77** (4)`.
- An intra-doc `:ref:` to a not-yet-anchored `.. note::` DOES warn
  under `-W` (L-002) — add the `.. _label:` to the note in the SAME
  edit. Cross-doc `:ref:` to the new section resolves silently when
  the target exists.
- Wrote `:meth:` for a CLASS name (`SNSolver`) — renders plain-text,
  no warning (L-002). Grep-gated every code-xref; fixed to `:class:`.

**Branch reality:** this landed on the worktree branch
`refactor/operator-inverse-algebra` where the doc-debt baseline is
`-W` CLEAN (EXIT=0, zero warnings) — NOT the 1-warning `mesh.py`
baseline some notes cite. Always run the `-E` baseline first to learn
the real acceptance floor for the branch you are on.
