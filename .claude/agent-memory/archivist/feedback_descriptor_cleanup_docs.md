---
name: BC descriptor-cleanup docs (Issue #186 / B3 + β2)
description: Pattern for documenting an architectural cleanup that retires a same-page interim (Option A → pure-descriptor). Predecessor as anti-pattern; type-system enforcement framing.
type: feedback
---

Pattern for documenting an architectural close-out that retires a
previously-documented interim approach **on the same page** — where
the predecessor was the chosen path of the immediately-preceding
session (here Issue #176 / C176.3's Option A) and the new path
moves to a strictly purer architecture.

**Why:** This is the third architectural wave in the BC trace-law
refactor (Wave 12 → Issue #188+#176 → Issue #186). Each predecessor
explicitly chose a tradeoff that the next wave overrode. The docs
must preserve the design history (rejected paths are load-bearing
intellectual content) without leaving stale "this is the chosen
path" claims that confuse readers about the current state.

**How to apply:** When closing a same-page interim:

1. **Replace the canonical-section anchor**, do NOT just edit in
   place. The pre-cleanup `bc-option-a-signatures` was the chosen
   path's documentation; in the cleanup it becomes a retired form.
   Delete the entire section, replace with a new section under a
   forward-facing anchor (`bc-trace-law-descriptor-model`). The
   predecessor's content moves into a "What was tried and rejected"
   subsection of the new section — full retrospective with
   rationale, not just a one-line dismissal.
2. **Cross-ref-rewrite-on-retire.** Every `:ref:` to the deleted
   anchor (in this same doc and in sister docs) must be flipped to
   point at the new anchor with retrospective framing. Three sites
   needed flipping: `boundary_conditions.rst:1885` (own-page
   forward pointer), `discrete_ordinates.rst:2543` and `:2691`
   (sister-doc references). Building with `-W` catches every
   undefined-label hit.
3. **Predecessor enumeration in numerical order, not chronological.**
   The "What was tried" subsection should list each predecessor
   approach as Option A / β1 / etc with the SAME enumeration the
   plan file used. Match the plan vocabulary verbatim so a reader
   landing on the section from a plan link finds the same names.
4. **Anti-pattern catalog rewrite, not delete.** The previous wave's
   anti-pattern entry "Option B (rejected)" is now the chosen path;
   flip it to "Option A (chosen-then-retired)". Add a NEW entry for
   the new failure mode that the cleanup enables ("calling apply on
   a raw BC descriptor — static type error"). Use the existing
   numbered list; bump existing entries' numbers if needed for
   logical ordering.
5. **Type-system-enforcement framing.** When the cleanup converts a
   convention into a static type contract (here: "law is not a
   callable" was a docstring claim → after #186 it's a static type
   error), this framing IS the load-bearing improvement. Use it
   liberally — "enforced by the type system, not by convention" is
   the one-line summary that captures why the cleanup landed.
6. **`bc-rank-n-algebra` rewrite around the two-type-families
   distinction.** When a section's name is unchanged but its body
   describes a fundamentally different algebra (Wave-0 operator
   algebra over realized leaves → descriptor-tree algebra over
   unrealized laws), don't rename the anchor — rewrite the body
   and preserve the section title. The anchor is load-bearing for
   cross-doc references; rewriting the body sharpens the content
   without breaking links. The new body needs a list-table
   contrasting "descriptor tree" vs "operator tree" (node types,
   `apply` presence, when built) — this contrast IS the
   architectural payoff.

**Predecessor retrospective subsection size** — multi-paragraph
when the predecessor was a deliberate design landing (Option A
got ~80 LoC of retrospective; β1 got ~30 LoC). NOT a one-line
dismissal: the user needs to know *why* each approach was
chosen-then-retired, not just that it was.

**Numerical evidence** — when foundation/L1 tests pin the new
contract (here `test_law_composition.py` has 18 tests), cite the
exact count in the V&V status admonition. Always grep the test
file for `def test_` to confirm the count before claiming it.

**Module-docstring depth-of-coverage rule** — the package-level
`__init__.py` docstring is the package "executive summary"; rewrite
to the new architecture and add a forward-pointer to the theory
page section that carries the full retrospective. Add a
descriptor-tree-composition-module subsection (briefly describing
`_composition.py`'s purpose + the three exports). Cross-references
at the bottom must list `_composition.py` alongside `_base.py` /
`SNBoundaryRealizer` / `realize_recursively`.

**Wave-E retrospective treatment** — historical Wave-E Round 3
prose referencing `BoundaryOperator.apply` is now wrong about the
current contract but right about the architectural conclusion.
Don't delete; add a parenthetical note: "(Post Issue #186 the
law itself is a pure descriptor; the realiser produces the
callable. See :ref:`bc-trace-law-descriptor-model`.)" — preserves
the historical narrative while pointing forward.

**Acceptance gate** — Sphinx `-W` build warning count unchanged
from pre-edit (here 9 baseline pre-existing warnings: paramref
unknown role, two title-underlines, six unknown skill docs). NEVER
chase a "zero warnings" goal when the baseline is non-zero — the
baseline warnings are pre-existing project-state issues that are
out-of-scope for the cleanup.

**Role-name pitfalls** — `:pydata:` and `:pyattr:` are NOT valid
RST roles in this build (they require `sphinx-paramref` or a
similar third-party extension). Use plain ``literal`` for ClassVar
names and `:attr:` for instance attributes; the build will catch
the wrong role with an explicit "Unknown interpreted text role"
error.
