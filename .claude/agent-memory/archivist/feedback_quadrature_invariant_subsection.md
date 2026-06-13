---
name: quadrature-invariant-subsection
description: Documenting a quadrature-table structural invariant a closure seed silently relies on (#193 re-scope into the ERR-058 coupled-pole section); two-independent-code-sites "emergent invariant" framing + Mode-7-at-operator-internals warning + foundation-gate-as-tripwire
metadata:
  type: feedback
---

Documenting an implementation-level invariant that a documented physics
seed RELIES ON (Issue #193, re-scoped, added into the ERR-058 close-out
section of `discrete_ordinates.rst`).

**The shape that worked** (a focused H4 subsection placed ADJACENT to the
math block whose seed it underwrites — manifestation (a)'s coupled-pole
continuity eq, NOT at the end of the close-out):

1. **Lead with the concrete code line the invariant lives in** — quote
   the exact expression (`outflow_at_inner.T[quad.reflection_index("x")]`)
   across all twin sites (both matvec twins + SI sweep twin), so the
   reader knows the single index that is load-bearing.
2. **State the invariant crisply** in prose AND as a named eq-label
   (the 3 equalities: μ_x sign-flip, μ_y/μ_z held → intra-level). The
   eq-label is `documented` vv-status (structural/representational
   identity, NOT a solver claim — pinned by a foundation gate, not a
   flux/eigenvalue). Add the `.. vv-status:` DIRECTIVE + rationale
   comment naming the foundation gate as the verifiable content.
3. **"Why the physics demands it"** — the WHY (pole continuity is at
   FIXED axial direction; only the radial cosine reverses at r=0).
4. **"Why it holds by construction today"** — the EMERGENT-INVARIANT
   framing: TWO INDEPENDENT code sites conspire (the reflection table
   flips only μ_x + holds μ_z; the level is keyed on μ_z). Neither site
   enforces it alone → it's a property of their agreement. This framing
   is what makes the gotcha legible.
5. **`.. warning::` = the silent-corruption surface** — frame as
   "Mode-7 blind spot at the OPERATOR-INTERNALS level" (vv-principles:
   the same flat-field blindness that hid ERR-058 itself; on flat ψ the
   mirror value == own value so every flat-flux gate stays green +
   scalar/particle-balance blind via α-dome telescoping). Name the exact
   future change that would break it (a new cubature reflection table; a
   refactor of `_compute_sphere_reflection_partners`).
6. **"Why it now has its own foundation gate"** — frame the gate as the
   TRIPWIRE that turns silent→loud: it pins the invariant DIRECTLY as a
   quadrature-table property (not through any flux/eigenvalue), so it
   fires BEFORE any solver runs. Cite both foundation tests; note the
   `verifies(...)` on the first ties the table-invariant to the
   continuity eq it underwrites.
7. **`.. note:: Re-scope of Issue #N`** — one honest paragraph: the
   ORIGINAL target, why it DISSOLVED (a keystone deletion made the old
   concern vacuous), where the load-bearing invariant MOVED to.

**Cross-ref gotcha caught this session:** I drafted a `:ref:` to an
invented anchor (`sn-bulk-boundary-split`) for the Wave-O B-split. It
did NOT exist. Grep-gate (`grep -rn "_<label>:" docs/`) caught it BEFORE
the build (intra-doc dangling `:ref:` WOULD warn under -W; cross-doc
renders plain-text SILENTLY — neither is acceptable on correctness
grounds). Fix: repoint at the REAL canonical anchor for the concept
(`bc-extraction` in `operator_algebra.rst:1639`, the bc-extraction
record). ALWAYS grep `^\.\. _<label>:` for every `:ref:` target before
building; verify cross-doc resolution in the built HTML
(`href="operator_algebra.html#bc-extraction"`).

**Verify-before-citing discipline that paid off:** every code claim in
the brief was checked against live source (operator.py:493,976;
loss_representation.py:2057/2110/2161; directional.py:176;
rules_sphere.py:213; rules_product.py:127,133). The "cylinder level
indexed by a non-x cosine" claim was the load-bearing one — confirmed
`level_indices` groups on |μ_z| (sphere/LS) or fixes μ_z per level
(product), NEVER on μ_x. The test file's exact function names + markers
were verified before being cited as `:func:` targets.

**Build gate:** baseline `-E` = 1 pre-existing warning (Mesh1D.from_geometry
:paramref:). Post-edit `-E` = identical set (EXIT=0, diff IDENTICAL).
Test-module `:func:` refs (`tests.sn.sweep.curvilinear...`) render
plain-text by existing-page convention (tests not automodule'd) —
consistent with the surrounding verification-chain section; SAFE, no
warning.

**Quality self-score** (rubric): derivation-depth 5 (full physics WHY +
the construction-proof via two independent sites + the linear-in-μ
exactness inherited from the parent block), cross-references 5 (every
twin site + both gate functions + the re-scope's dissolved concern
anchor), numerical-evidence N/A (this is a structural invariant, not a
numerical result — the foundation gate IS the evidence, cited), failed-
approaches 5 (the dissolved #193 original target documented as honest
history), code-traceability 5 (eq-label + all three twin sites named),
derivation-source N/A (no SymPy — quadrature-table invariant verified by
a foundation test, the correct source for a software invariant).
