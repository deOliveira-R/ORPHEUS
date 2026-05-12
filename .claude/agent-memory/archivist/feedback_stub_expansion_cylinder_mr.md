---
name: Stub-expansion cylinder MR Phase 1b (Cardinal Rule 3 rich narrative from method-implementer stubs)
description: Pattern for archivist expansion of method-implementer-shipped stub labels into rich Sphinx narrative, using closeout-memo + verification-plan + Branch-1/Branch-2 sources as the synthesis seed
type: feedback
---

# Cylinder MR Phase 1b stub-to-narrative expansion (2026-05-12)

Replaced 8 `.. todo::` stub labels in `docs/theory/trajectory_resolvent.rst`
`.. _peierls-greens-cylinder-mr:` section with full math + prose +
numerical-evidence tables, mirroring the sphere MR section's depth/rigor.
~880 LoC growth, zero new Sphinx warnings.

**Why:** Cardinal Rule 3 (Sphinx is the LLM's brain) requires every theory
page to carry full derivations, design rationale (not just what — WHY),
gotchas, and numerical evidence. Method-implementer ships stubs (per
`algebra-of-record` § "Sphinx stub vs rich narrative"); archivist expands.

**How to apply:**

1. **Read sources in this order**: closeout memo (carries the gate-results
   table + structural-risk falsification record) → verification plan
   (carries the math claim hierarchy + V&V level matrix per gate +
   anti-patterns avoided) → Branch-1 SymPy module (carries the canonical
   algebra-of-record identities the rich narrative cites verbatim) →
   Branch-2 production code (carries the implementation details
   the prose references via `:func:`/`:class:`) → tests (carries the
   numerical-evidence numbers via achieved-vs-target). Then the EXISTING
   sphere MR section for style/depth template.

2. **Per-label expansion shape** (one section per `:label:`):

   - **Restate the math** with full equation block (preserved verbatim from
     the stub for label-resolution stability).
   - **`:vv-status: <label> documented`** or `verified` comment block
     immediately after the `.. math::` directive, matching the pattern in
     the existing sphere MR sections (line 1393).
   - **Why-this-equation prose**: 2-4 paragraphs explaining the geometric
     intuition + the structural distinction from the sphere analog.
   - **Algebraic-ancestor cross-ref**: the SymPy `derive_*` function +
     its foundation test, as a `:func:` link.
   - **Implementation prose**: which Branch-2 function + which oracle
     primitives realise the equation, as `:func:`/`:class:` links.
   - **Numerical-evidence table** when the equation carries a gate
     result: target column + achieved column.
   - **Gotcha / structural-risk falsification**: for each risk flagged in
     the verification plan §6.X, document whether it materialised and what
     the gate evidence shows (e.g., "tangential grazing-ray 0/0
     cancellation generalises because both numerator and denominator
     carry the same axial Jacobian — Gate 3 hits 3.5e-12 at the K=1
     reduction, which would be impossible with NaN").

3. **Key Facts admonition at the top** mirrors the sphere MR / Peierls
   pattern: 5-7 bullets that frontload the load-bearing claims (phase-
   space, segmentation strategy, Cardinal Rule 6 enforcement,
   structural-independence anchor, V&V limitations).

4. **Gate-results table at the section head** (per Cardinal Rule 3
   "numerical evidence" mandate) using `.. list-table::` with columns:
   Gate / Type / Pillar / Target / Achieved. Pulled verbatim from the
   closeout memo. Captures the verification chain in 12 rows.

5. **V&V limitations section** explicit and named — the verification plan
   defines Gate 5 (literature MR benchmark) as TENTATIVE; the closeout
   confirms NO benchmark found. Archivist documents this as a known V&V
   limitation with the candidates inspected and the future-search
   suggestion (Russian/Japanese archives, Sanchez-Pomraning 1989). This
   is the load-bearing transparency piece that survives in Sphinx after
   the closeout memo is archived.

6. **L0 bugs caught subsection** records the development bug-hunt
   outcome: in this case "no solver bugs; two test-data-shape errors
   in the test code, fixed inline; no ERR-NNN catalog entries". This is
   the algebra-of-record proof point — Branch-1 SymPy identities
   derived first means Branch-2 production code lands clean.

**RST gotchas hit this session:**

- **`:ref:`<skill-name>`` does NOT resolve.** Skills (`vv-principles`,
  `numerical-bug-signatures`, `algebra-of-record`) are `.md` files under
  `.claude/skills/`, not Sphinx documents. Cite them as plain prose:
  `the \`\`vv-principles\`\` skill § "..."`. The pre-existing
  `cross_method.rst` warning floor (6 `unknown document: '/skills/...'`
  warnings) is the canonical instance — the warning-count gate is
  unchanged from baseline, do NOT try to "fix" those (they were
  pre-existing).
- **`:func:` cross-ref target lives where the FUNCTION is defined**,
  not where it's imported. `apply_variant_alpha_closure` is defined in
  `variant_alpha_core.py`, exported via `chord_oracle.py` — the Sphinx
  cross-ref MUST use `...variant_alpha_core.apply_variant_alpha_closure`
  even though oracle callers import it from `chord_oracle`. Sphinx
  resolves by definition-site, not import-site.
- **`:doc:`/development`` works for top-level docs but NOT for the
  `:doc:`/skills/X`` invocations** in `cross_method.rst` (those are the
  unresolved warnings; the skill is not in the Sphinx toctree).
- **Pre-existing 9-warning baseline floor** must be diffed pre/post-edit
  (per `feedback_closeout_docs.md` "warning-count diff as the
  acceptance gate"). For this session: pre-edit 9 warnings (all in
  unrelated files), post-edit 9 warnings (same set). Gate satisfied
  even though `-W` exits non-zero. NEVER attempt to suppress pre-
  existing warnings as part of an unrelated docs edit.

**vv-status discipline (matches existing sphere MR section):**

- `:vv-status: <label> documented` for definitional / notation
  introductions (the trajectory parametrisation, the optical-depth
  definition, the bounce-sum closure shape).
- `:vv-status: <label> verified` for testable claims with numerical
  evidence (the homogeneous-limit reduction passing Gate 1, the WM-72
  cross-check passing Gate 2, the k_∞ recovery passing Gate 3, the
  interface continuity passing Gate 4, the quadrature convergence
  passing Gate 6).
- The rationale comment block precedes the `:vv-status:` line and uses
  the `.. (vv-status rationale) <category>:` prefix matching the
  sphere MR pattern at line 1393.

**Cross-link discipline:**

- Cite `:eq:` for equations on the same page (uses Sphinx equation
  resolver, e.g. `:eq:`peierls-greens-mr-trajectory-segments``).
- Cite `:ref:` for section labels on OTHER pages (e.g.
  `:ref:`theory-sood-registry``).
- Cite `:func:`/`:class:`/`:mod:` for Python symbols, full dotted path
  with `~` prefix to show only the leaf name.
- Cite `:doc:` for cross-document path references (e.g.
  `:doc:`/development``).

Self-rating: 5/5 derivation depth (full piecewise-τ algebra, V_α1_cyl_mr
algebraic ancestor, 0/0 cancellation analysis); 5/5 cross-references
(every `:label:` linked to its `derive_*` SymPy function and to its
test gate); 5/5 numerical evidence (12-row gate-results table at section
head + 3 per-gate evidence tables); 5/5 failed approaches (V&V Gate 5
limitation explicit, Gate 4 spline-across-jump floor documented, K=1
2G asymmetric SigS as ERR-002 anti-Mode-#6 explicit); 5/5 code
traceability (every equation linked to `_cylinder_trajectory_segments_2d`
or `MultiRegionCylinderChordOracle` or `apply_variant_alpha_closure`);
5/5 derivation source (every algebraic claim cited the originating
`derive_*_cylinder_mr*` SymPy module).
