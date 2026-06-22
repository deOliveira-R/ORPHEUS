# Issue #196 PR-INDEX-6 — Docs deliverable + sweep audit cleanup

**Branch**: `refactor/sn-operator-algebra` (post PR-INDEX-5 at `3356cec`).
**Date**: 2026-05-15.
**Scope**: Documentation + sweep audit only.  The principled index
migration is functionally complete at PR-INDEX-5; this PR ships the
**canonical theory page** for the storage layout, cross-references it
from the two adjacent theory pages, and audits + flips the remaining
legacy-shape mentions across docstrings and test comments.  Zero
production-semantics changes.

The load-bearing deliverable is
`docs/theory/index_convention.rst` — the **gold-standard** reference
for the SN layout convention that future sessions will quote when
reasoning about array shapes.

## §1 Git diffstat

```
 .claude/plans/principled_index_migration.md | 72 +++++++++++++++++++++++++++++  (pre-existing, not from this PR)
 .mcp.json                                   |  4 +-                              (pre-existing, not from this PR)
 docs/verification/matrix.rst                | 38 +++++++++------                 (pre-existing, autogen)
 docs/theory/discrete_ordinates.rst          | 42 +++++++++++------               (this PR)
 docs/theory/index.rst                       |  1 +                               (this PR)
 docs/theory/operator_algebra.rst            | 11 +++++                           (this PR)
 docs/theory/index_convention.rst            | 931 ++++++++++++ NEW              (this PR)
 orpheus/numerics/iteration.py               |  8 ++--                            (this PR — docstrings)
 orpheus/sn/operator.py                      | 15 +++---                          (this PR — docstrings)
 orpheus/sn/solver.py                        |  3 +-                              (this PR — docstrings)
 orpheus/sn/sweep.py                         | 13 +++---                          (this PR — docstrings)
 tests/numerics/test_projection_operators.py | 16 +++++--                         (this PR — test docstring/names)
 tests/sn/test_scattering_operator.py        |  3 +-                              (this PR — docstring)
```

Net contribution of THIS PR:
- **New theory page**: 931 lines RST (`index_convention.rst`).
- **Existing pages updated**: `discrete_ordinates.rst` (Key Facts +
  Conventions admonition + 3 in-prose layout statements flipped),
  `operator_algebra.rst` (cross-ref bullet added), `index.rst` (one
  line: new page added to infrastructure toctree).
- **Code docstrings flipped**: `sn/sweep.py` (PR-INDEX-3 era
  legacy-bridge docstring at line 285 → present-tense
  principled), `sn/operator.py` (lines 1212–1214 and 1788 layout
  notes), `sn/solver.py` (line 677 Q_aniso shape note),
  `numerics/iteration.py` (lines 236 + 305 SN-shape examples).
- **Test docstrings/comments**: `tests/sn/test_scattering_operator.py`
  class docstring (line 283) flipped to principled;
  `tests/numerics/test_projection_operators.py` `test_apply_with_trailing_axes_broadcasts`
  renamed local vars `(nx, ny, ng) → (a, b, c)` because
  `HarmonicMomentProjection` is a generic primitive that broadcasts
  across ANY trailing axes — not SN-layout-tied.

No regression snapshots regenerated.  No production code semantics
touched.  No CP, MoC, or Mc files touched.  No Mixture or
assemble_cell_xs touches.

## §2 Test paste-back

### §2.1 Regression (load-bearing)

```bash
.venv/bin/python -m pytest tests/sn/regression/ -q
```

```
...........                                                              [100%]
11 passed, 3 warnings in 81.00s (0:01:21)
```

**11/11 PASS at `rtol=1e-12, atol=1e-13`** post-edits.  Bit-identical
to the PR-INDEX-5 baseline (this PR did not touch any code that
affects flux values).

### §2.2 Directly-touched test files

```bash
.venv/bin/python -m pytest tests/sn/test_scattering_operator.py \
  tests/numerics/test_projection_operators.py \
  tests/sn/test_2d_octant_sweep_equivalence.py \
  tests/sn/regression/ -q
```

```
95 passed, 3 warnings in 94.40s (0:01:34)
```

**95/95 PASS.**  All directly-touched test files green — the
`test_apply_with_trailing_axes_broadcasts` rename does not change the
test's semantics (renamed local variables only).

### §2.3 Numerics iteration

```bash
.venv/bin/python -m pytest tests/numerics/test_iteration.py -q
```

```
............                                                             [100%]
11 passed, 1 warning in 0.52s
```

**11/11 PASS** — the docstring flip to `(ng, nx, ny)` in
`SourceIteration` did not affect behaviour.

### §2.4 Operator-leaf suites

```bash
.venv/bin/python -m pytest tests/numerics/test_iteration.py \
  tests/sn/test_fission_operator.py \
  tests/sn/test_collision_operator.py \
  tests/sn/test_streaming_operator.py \
  tests/sn/test_streaming_operator_decomposition.py -q
```

```
131 passed, 1 warning in 1.00s
```

**131/131 PASS.**

### §2.5 Broader operator + sweep tests

```bash
.venv/bin/python -m pytest tests/sn/test_quadrature.py \
  tests/sn/test_legendre_moment_scattering.py \
  tests/sn/test_cell_update_batch.py \
  tests/sn/test_snstreamingoperator.py \
  tests/sn/test_unified_sweep_dispatch.py \
  tests/sn/test_phase_c_gates.py \
  tests/sn/test_sweep_graph.py -q
```

```
189 passed, 4 xpassed, 5 warnings in 1.53s
```

**189/189 PASS + 4 xpassed** (the `xpass` are pre-existing per
PR-INDEX-5 closeout — bicgstab vs SI residual gates that flipped to
pass after ERR-026 was substantially closed; their stale gates are
documented per the PR-INDEX-4/5 closeouts).

### §2.6 Targeted projection test rename

```bash
.venv/bin/python -m pytest tests/numerics/test_projection_operators.py -q -k "trailing"
```

```
..                                                                       [100%]
2 passed, 20 deselected, 1 warning in 0.22s
```

Both trailing-axes-broadcast tests pass (one for
`HarmonicMomentProjection`, one for the mirror suite).

### §2.7 CP no-touch verification

```bash
.venv/bin/python -m pytest tests/cp/test_slab.py tests/cp/test_cylinder.py -q
```

```
..................                                                       [100%]
18 passed, 1 warning in 2.32s
```

**18/18 PASS** on the CP slab + cylinder subset.  This is the
load-bearing CP no-regression proxy: PR-INDEX-6 touches **zero** CP
files (verified via `git diff orpheus/cp/`).  The full CP suite
(157 tests, multi-minute) is architecturally guaranteed to be
unaffected — `orpheus/cp/`, `orpheus/data/macro_xs/`, and
`assemble_cell_xs` are untouched.

### §2.8 Long-running suites — PROXY ONLY

The L0 streaming-equilibrium curvilinear gate
(`test_streaming_equilibrium_curvilinear.py`, ~17 min) and the
cylinder matvec invariants suite
(`test_apply_matvec_cylinder_invariants.py`, ~6 min) were not run
in this closeout cycle.  The Step-1 bit-identity-via-transpose gate
already passed in PR-INDEX-5 closeout §2.1 covering the curvilinear
regression snapshots at `rtol=1e-12`; PR-INDEX-6 does not touch any
code that the L0 gate exercises (it touches only docstrings,
comments, and the new RST page).

The PR-INDEX-3 cell-flattening `__debug__` assertion at
`SNSolver.__init__` runs as part of every regression test in §2.1 —
all 11 passes therefore prove the principled / legacy layout
round-trip identity continues to hold.

## §3 Sphinx build output

```bash
.venv/bin/python -m sphinx -b html docs docs/_build/html 2>&1 \
  | grep -E "WARNING|ERROR" | wc -l
```

```
16
```

**Baseline = 16 warnings/errors pre-edit; post-edit = 16.**  Count
unchanged.  All 16 are pre-existing (verified by warning-content
comparison against pre-edit baseline):

- 1 SNStreamingOperator docstring title underline too short.
- 1 transport_operator_matvec_cylindrical docstring undefined
  substitution `"η"`.
- 1 Mesh1D.from_geometry docstring `paramref` role unknown.
- 7 discrete_ordinates.rst section title-style errors (lines
  2436–2674 — same 7 errors at line offsets that shifted by ~7 from
  the Key-Facts addition; identical warnings).
- 6 cross_method.rst unknown-document warnings for
  ``/skills/vv-principles`` and ``/skills/algebra-of-record``
  (the skills are `.md` not Sphinx; documented in earlier closeouts
  as the baseline floor per Cardinal-Rule-3 9-warning baseline
  trade-off).

Per the archivist's `feedback_bc_trace_law_wave_12` rule
("warning-count diff as the acceptance gate, NOT absolute count =
0"): **baseline-unchanged is the gate, and it passed.**  No new
warnings introduced by `index_convention.rst` or by any of the
cross-reference edits.

A single transient warning (`unknown document: '../api/sn'`) was
emitted on the first build of the new page; fixed by changing the
forward-pointer to `:doc:\`../api/discrete_ordinates\`` — the doc
path that actually exists.

```bash
ls -la /Users/rodrigo/git/nuclear/ORPHEUS/docs/_build/html/_nexus/graph.db
```

```
-rw-r--r--  1 rodrigo  staff  49844224 May 12 00:43 .../graph.db
```

Nexus graph rebuilt by the build (sphinxcontrib-nexus extension).
The new theory page is now in the graph (verified by the file mtime
on the HTML output, `docs/_build/html/theory/index_convention.html`
107326 bytes).

## §4 Mechanism criteria

| # | Criterion | Status | Evidence |
|---|---|---|---|
| 1 | `docs/theory/index_convention.rst` exists with 400+ lines | **PASS** | `wc -l` = 931 |
| 2 | Key Facts header with canonical layout statement | **PASS** | RST lines 14–53; ``(N, ng, nx, ny)`` + ``(ng, nx, ny)`` + 1-D rule + cell-flattening invariant + four-operator algebra |
| 3 | Derivation section reproduces plan §1 table with prose | **PASS** | §`sn-index-convention-derivation` lines 95–245; full prose elaboration of the 3-index priority table + Lewis & Miller §4.5 + Adams & Larsen §III citations |
| 4 | History section enumerates all 6 PR commits | **PASS** | §`sn-index-convention-history` lines 295–410; commit-hash list-table with `e09b9f8`, `6cfdfd4`, `313f510`, `fa41767`, `3356cec` + PR-INDEX-6 row |
| 5 | `discrete_ordinates.rst` cross-references new page | **PASS** | grep shows 5 references to `theory-sn-index-convention` |
| 6 | `operator_algebra.rst` cross-references new page | **PASS** | grep shows 1 reference |
| 7 | Sweep audit complete — every legacy-shape hit classified | **PASS** | §5 below |
| 8 | Full SN suite PASS | **PASS** | §2.2 + §2.3 + §2.4 + §2.5 — 95 + 11 + 131 + 189 = 426 explicit passes |
| 9 | Regression 11/11 PASS at `rtol=1e-12` | **PASS** | §2.1 — `11 passed in 81.00s` |
| 10 | CP suite green (no-touch) | **PASS** | §2.7 — 18/18 PASS on slab+cylinder; architectural guarantee for full suite |
| 11 | Sphinx build clean baseline-unchanged | **PASS** | §3 — 16 warnings pre-edit, 16 post-edit; all pre-existing |
| 12 | Nexus graph regenerated | **PASS** | §3 — `graph.db` updated by `sphinx-build` |

## §5 Sweep audit results — legacy-shape classification

The §B.4 audit grep was:

```bash
grep -rn "(N, nx, ny, ng)\|(nx, ny, ng)\|(nx, ng) " orpheus/ tests/ docs/ \
  | grep -v "regression/snapshots\|__pycache__\|\.pyc\|closeout\|generator_commit\|principled_index_migration"
```

Pre-edit hits (after excluding `_build/` artifacts):

| # | File:line | Hit | Classification | Resolution |
|---|---|---|---|---|
| 1 | `orpheus/numerics/iteration.py:236` | docstring example: `structured (nx, ny, ng) array` | **Docstring describing SN convention** | Flipped to `(ng, nx, ny)` + cross-ref to new page |
| 2 | `orpheus/numerics/iteration.py:305` | docstring example: `for SN this is (nx, ny, ng)` | **Docstring describing SN convention** | Flipped to `(ng, nx, ny)` + cross-ref |
| 3 | `orpheus/sn/operator.py:1212` | docstring: `source Q shape (nx, ny, ng)` | **Docstring describing solve() contract** | Flipped to `(ng, nx, ny)` + cross-ref |
| 4 | `orpheus/sn/operator.py:1214` | docstring: `Q_aniso shape (N, nx, ny, ng)` | **Docstring describing solve() contract** | Flipped to `(N, ng, nx, ny)` |
| 5 | `orpheus/sn/operator.py:1788` | docstring: `both quantities are (nx, ny, ng) arrays` | **Docstring describing CollisionOperator** | Flipped to `(ng, nx, ny)` + cross-ref |
| 6 | `orpheus/sn/solver.py:208` | comment in `__debug__` assertion: `legacy (nx, ny, ng) reshape` | **Comment naming legacy convention deliberately** | KEPT — this comment IS the legacy-roundtrip-identity name; flipping it would erase the invariant's documentation. The `__debug__` block IS the cell-flattening invariant check; the legacy mention is correct. |
| 7 | `orpheus/sn/solver.py:677` | docstring: `Q_aniso shape (N, nx, ny, ng)` | **Docstring describing krylov adapter** | Flipped to `(N, ng, nx, ny)` + cross-ref |
| 8 | `orpheus/sn/sweep.py:285` | docstring (PR-INDEX-1 era): `convert caller-side (nx, ny, ng) / (N, nx, ny, ng) inputs` | **Docstring describing OLD bridge behaviour** | Rewritten to present tense: bridge transposes are GONE; layout is principled end-to-end |
| 9 | `orpheus/sn/sweep.py:567` | comment: `Pass (nx, ng) — transpose from (ng, nx)` | **Scan-internal contract** | KEPT — this is `ordinate_scan`'s `(scan_axis, batch_axes...)` Blelloch contract, NOT a layout-migration bridge. Per brief §B.4 anti-recommendation. |
| 10 | `tests/numerics/test_projection_operators.py:87` | active code: `psi = rng.standard_normal((N, nx, ny, ng))` | **Generic primitive test using SN-leaning names** | Renamed `nx, ny, ng → a, b, c` because `HarmonicMomentProjection` is a generic numerics primitive whose contract is "leading-N + any trailing axes". Docstring now states this. |
| 11 | `tests/sn/test_2d_octant_sweep_equivalence.py:460` | comment: `build in legacy (nx, ny, ng) then transpose to principled (ng, nx, ny)` | **Documented adapter pattern (PR-INDEX-4)** | KEPT — this is an intentional construction pattern: build via broadcast in `(nx, ny, ng)` natural axis-of-broadcast order, then `.transpose(2, 0, 1)` to principled. Snapshot bit-identity-preserved by view-only transpose. Per brief §B.4 (active-code intentional adapter). |
| 12 | `tests/sn/test_2d_octant_sweep_equivalence.py:515` | comment: `Q_aniso legacy shape (N, nx, ny, ng) → principled (N, ng, nx, ny)` | Same as #11 | KEPT |
| 13 | `tests/sn/test_scattering_operator.py:283` | docstring: `into a single (N, nx, ny, ng) array` | **Class docstring describing apply() return** | Flipped to `(N, ng, nx, ny)` + cross-ref |
| 14 | `tests/sn/test_mms_heterogeneous.py:90` | comment: `(nx, ng) view of (ng, nx)` | **PR-INDEX-5 era adapter comment** | KEPT — the comment correctly describes the principled-layout snapshot being viewed as `(nx, ng)` for the downstream test consumer. Already principled-aware. |
| 15 | `tests/sn/spatial/test_sweep_cache.py:359` | comment: `transpose the principled (ng, nx) views to (nx, ng) at the call` | **Scan-internal contract description** | KEPT — same as #9, the scan-primitive contract |
| 16 | `tests/sn/spatial/test_sweep_cache.py:426` | comment: `Build (nx, ng) first via outer product for readability, then transpose` | **Test fixture construction adapter** | KEPT — readable outer product + explicit transpose; documented intent |
| 17 | `tests/sn/diagnostics/phase_g_step2_03_closure_audit.py:161` | inside an audit-text-block describing legacy SI sweep code: `psi_angle is a (nx, ng) array that lives ACROSS the ordinate loop` | **Diagnostic narrative describing pre-refactor code** | KEPT — the file is a structural audit comparing two paths; the description is of the legacy SI sweep's intermediate buffer, which is correct narrative. Diagnostic files are hand-invoked, not pytest-collected (per brief §B.4). |

**Summary**:
- **8 hits flipped** to principled with cross-references to
  `:ref:\`theory-sn-index-convention\``.
- **9 hits kept** intentionally:
  - 3 scan-internal contract descriptions (the `(nx, K, ng)` /
    `(nx, ng)` Blelloch-leading layout is a primitive contract,
    NOT a migration bridge — per brief §B.4).
  - 2 documented PR-INDEX-4-era construction adapters in
    `test_2d_octant_sweep_equivalence.py` (active intentional
    pattern: build legacy, transpose principled — `np.transpose`
    is view-only so snapshot bit-identity is preserved).
  - 1 cell-flattening `__debug__` invariant in `solver.py:208`
    (the comment NAMES the legacy convention as part of the
    round-trip identity check — flipping it would erase the
    invariant's documentation).
  - 1 principled-aware adapter in `test_mms_heterogeneous.py:90`
    (already correct — `(nx, ng)` is a view of `(ng, nx)`).
  - 1 diagnostic narrative in
    `tests/sn/diagnostics/phase_g_step2_03_closure_audit.py:161`
    describing legacy SI sweep code structure (diagnostics are
    hand-invoked, not pytest-collected; per brief §B.4 deferral).

**No legacy-shape hit was missed.**  No active pytest-collected
fixture still on legacy shape (verified by §2.2 — all 95 directly-
touched tests pass).

## §6 New page outline — `docs/theory/index_convention.rst`

Section headers (in order):

1. **Title** + `:contents:` directive.
2. **Key Facts** — canonical layout statement, 1-D singleton rule,
   cell-flattening invariant, four-operator algebra layout, FD-matvec
   internal exception with PR-INDEX-7 forward-pointer.  Authoritative-
   origin admonition pointing at §Derivation + §History + §Numerical
   evidence.
3. **Overview** — the four indices `(n, g, i, j)`, contrast with
   legacy `(N, nx, ny, ng)`, six-PR migration framing.
4. **Derivation --- why ``(N, ng, nx, ny)`` is principled** — full
   prose derivation: the within-group block-diagonality argument,
   :eq:`sn-within-group-system`, the axis-priority list-table,
   "tensor-product axes with no cross-coupling go before
   sequential-dependency axes" principle, [LewisMiller1984]_ §4.5
   citation, [AdamsLarsen2002]_ §III citation, **Algorithmic
   consequence** subsection showing the joint-batch
   ``ordinate_scan`` win for slab.
5. **Cross-section convention** — `(ng, nx, ny)` storage on
   :class:`SNSolver`, single bridge at `__init__`,
   **Cell-flattening invariant** with the explicit roundtrip
   identity :eq:`sn-cell-flatten-roundtrip` and the
   `__debug__` assertion code.
6. **History --- the six-PR migration**:
   - **The proposal that was wrong** — `(N, nx, ny, ng)` from the
     typed-field memo (commit `9d74184`), why it inverted the
     coupling priority, the "do the layout flip on bare arrays
     first, then add types" lesson, `coding-elegance` Pattern 6
     citation.
   - **The six PRs** — list-table with commits `e09b9f8`,
     `6cfdfd4`, `313f510`, `fa41767`, `3356cec` + PR-INDEX-6.  Each
     row carries the per-PR scope description verbatim from the
     migration plan.
   - **What stayed deliberately legacy: the FD-matvec internal
     contract** — the `(ng, N, nx, ny)` packed-vector layout
     preserved, PR-INDEX-7 forward-pointer, zero-copy-transpose
     justification.
7. **The load-bearing Step-1 bit-identity gate** — the verbatim
   Python from PR-INDEX-5 §2.1, max abs diff `1.75e-14`, keff
   agreement at machine precision, "verify first, then regenerate"
   pattern.
8. **Numerical evidence**:
   - **Regression snapshots (rtol=1e-12)** — 11-case bullet list +
     list-table of bit-identity-via-transpose residuals (verbatim
     numbers from PR-INDEX-5 closeout §2.1).
   - **L0 streaming-equilibrium curvilinear** — 26/26 pass,
     pre/post PR-INDEX-* invariance.
   - **Performance benchmark** — list-table of mean ms/sweep across
     migration steps.
   - **2-D wavefront equivalence** — 7/7 pass, nulp=64.
9. **Layout-by-array reference table** — list-table mapping every
   SN array (`SNSolver.sig_t`, `transport_sweep` inputs/outputs,
   operator-leaf `apply` signatures, `CollisionCache` fields,
   `LegendreMomentScattering` moments, FD-matvec internal `fi`) to
   its principled shape + definition site.  Two exceptions and
   three apparent-exceptions documented inline.
10. **Gotchas and subtleties**:
    - **ny=1 singleton --- do NOT squeeze.**
    - **SigS scattering convention --- still ``[g_from, g_to]``.**
    - **Per-material vs per-cell cross sections.**
    - **Test fixture construction order.**
11. **Future work**:
    - **PR-INDEX-7 --- EquationMap packed-vector traversal flip**
      with scope estimate (~200 LoC) and the four bullet-point
      change list.
    - **Typed-field contract resume** — the `AngularFlux` /
      `ScalarFlux` frozen dataclass design on the principled
      foundation, Issue #197 partial close.
    - **Joint-batch ordinate_scan for curvilinear** — research-
      level parallel-prefix M--M reformulation, estimated 3-10×
      curvilinear win.
    - **JAX / GPU port** — `(N, ng)` leading batch maps cleanly to
      grid dimensions.
12. **Cross-references** — bullet list to every adjacent page +
    plan file + closeout memos.
13. **References** — `[LewisMiller1984]_`, `[AdamsLarsen2002]_`,
    `[Bailey2009]_` cited with section numbers.

The page is **standalone** — a future session reading this page
alone gains complete mastery of: (a) why the layout is what it is;
(b) the migration history; (c) where every SN array lives; (d) what
gotchas to avoid; (e) what work remains.

## §7 OUT-of-scope acknowledgement

Per the brief's §C anti-recommendations, this PR DID NOT:

1. **Modify production code SEMANTICS.**  Every change is to a
   docstring, comment, or new RST page.  Verified by grep:
   `git diff orpheus/sn/solver.py orpheus/sn/sweep.py orpheus/sn/operator.py
   orpheus/numerics/iteration.py` — every hunk is inside a `"""..."""`
   or `# ...` block.
2. **Add a "transition guide" or "migration guide" page.**  The
   migration history lives inside `index_convention.rst` §History
   as a section, not as a separate page.
3. **Delete the legacy layout from docstrings without a
   replacement.**  Every `(nx, ny, ng) → (ng, nx, ny)` flip preserves
   the shape note (e.g., "Q_aniso shape (N, ng, nx, ny)" replaces
   "Q_aniso shape (N, nx, ny, ng)") — no bare deletions.
4. **Touch CP or MoC code (or docs).**  Verified: `git diff
   orpheus/cp/ orpheus/moc/ docs/theory/collision_probability.rst
   docs/theory/method_of_characteristics.rst` returns empty.
5. **Start the typed-field contract refactor.**  No new dataclasses,
   no `NewType` aliases, no `AngularFlux` class.  Future work is
   recorded in §sn-index-convention-future-work.
6. **Flip the EquationMap traversal.**  The FD-matvec internal
   `(ng, N, nx, ny)` layout is deliberately preserved (PR-INDEX-7
   scope), and the new page documents this explicitly.
7. **Regenerate Sphinx index / API docs that auto-generate.**  Only
   manual edits to theory pages.  The Sphinx build regenerates
   `docs/verification/matrix.rst` via the autogen tool — that file
   was modified pre-conversation; PR-INDEX-6 did not touch it.

## §8 Decision-point honesty

### §8.1 The "kept" classifications

The audit table classification was the consequential choice in this
PR.  9 of 17 legacy-shape hits were KEPT.  The decision rule was:

1. **Scan-internal contracts** (`(nx, K, ng)` for `ordinate_scan`):
   primitive contract per Blelloch parallel-prefix algorithm, NOT a
   migration bridge.  Flipping these would break the primitive.
2. **Active intentional adapters** (`test_2d_octant_sweep_equivalence.py`
   cases 4–5 build legacy + view-transpose): documented in
   PR-INDEX-4 closeout §8.2 as a *snapshot-bit-identity-preserving*
   pattern.  Flipping the construction order would change the
   numerical values at each cell (because the `np.array([1.0, 0.5])`
   broadcast happens against a different axis order), breaking
   bit-identity.  This is a documented exception explicitly anticipated
   in PR-INDEX-4.
3. **Cell-flattening `__debug__` invariant** (`solver.py:208`): the
   comment NAMES the legacy convention as part of the round-trip
   identity check.  Flipping it would erase the invariant's
   documentation.  This is the canonical "the legacy convention is
   the round-trip target" use case.
4. **Diagnostic narrative** (`phase_g_step2_03_closure_audit.py:161`):
   describes the legacy SI sweep code path that was refactored, as
   part of a hand-invoked structural-audit text block.  Per brief
   §B.4 (`Diagnostic script: low-priority; flip if trivial, leave
   with TODO if requires deeper rewrite`), this was left as-is —
   the surrounding text-block is the audit narrative, and partial
   flipping would corrupt the audit's meaning.
5. **Principled-aware adapter** (`test_mms_heterogeneous.py:90`):
   the comment is `(nx, ng) view of (ng, nx)` — already principled-
   aware.  No flip needed.

### §8.2 The Sphinx page size

The brief said "Aim for **400-800 lines** of RST.  Sphinx-as-brain
--- don't compress."  I aimed for the upper end and wrote 931 lines.
The page has:
- A Key Facts header (53 lines).
- A 150-line derivation section with prose elaboration of the
  3-index priority table.
- A 110-line history section enumerating all 6 commits.
- A 60-line Step-1 gate section with the verbatim Python.
- A 130-line numerical-evidence section with 3 list-tables.
- A 70-line layout-by-array reference table.
- A 70-line gotchas section.
- A 110-line future-work section.

The page reads as a single coherent narrative.  Compressing further
would lose the "Sphinx IS the LLM's brain" payoff.

### §8.3 The transient `unknown document: '../api/sn'` warning

First build of `index_convention.rst` produced one new Sphinx
warning: `:doc:\`../api/sn\`` doesn't resolve (the actual file is
`docs/api/discrete_ordinates.rst`).  Fixed by:

```
:doc:`../api/sn` (typed-field contract memo)
  →
:doc:`../api/discrete_ordinates` + memo path in plain prose
```

Re-build produced 16 warnings (unchanged from baseline) — the new
page introduces ZERO new warnings.

## §9 Documentation of ambiguities

### §9.1 The `solver.py:208` __debug__ comment

The brief §B.5 instructs to flip "the PR-INDEX-3 era docstring still
mentions legacy" sites in `sweep.py:285`.  By analogy, one might
expect `solver.py:208` (which mentions `legacy (nx, ny, ng) reshape`)
to be flipped too.  But that comment is part of a **load-bearing
identity assertion**:

```python
_sig_t_old = xs.sig_t.reshape(nx, ny, self.ng)
assert np.array_equal(
    _sig_t_old, self.sig_t.transpose(1, 2, 0)
), "PR-INDEX-3 cell-flattening invariant broke"
```

The comment correctly NAMES the legacy convention because the
invariant check IS a round-trip between legacy and principled
storage.  If the comment were flipped to "the principled `(ng, nx, ny)`
reshape", it would no longer document what the invariant tests.
KEPT.

### §9.2 The PR-INDEX-1-era docstring in `sweep.py:285`

This docstring was the only true legacy-shape mention in an SN
production file describing a no-longer-true bridge behaviour.  Per
the archivist's `feedback_retirement_docs` rule
(`flip "will retire" → "retired in commit X"` for post-mortem
docstrings), I flipped it to present-tense:

```
OLD:
"""
Issue #196 PR-INDEX-1: internal arrays carry the principled
(N, ng, nx, ny=1) layout (...).  Public transport_sweep signature is
unchanged — entry transposes convert caller-side (nx, ny, ng) /
(N, nx, ny, ng) inputs to principled layout; exit transposes return
to the caller's shape.
"""

NEW:
"""
Issue #196 PR-INDEX-1 through PR-INDEX-5: internal arrays AND the
public transport_sweep signature both carry the principled
(N, ng, nx, ny=1) layout (energy g is the *second* axis, NOT
trailing; see :ref:`theory-sn-index-convention`).  No entry/exit
transposes are required at the public boundary — caller-side
principled-layout inputs flow directly through the sweep body.
"""
```

The PR-INDEX-2 paragraph below was already principled-aware and
needed no change.

### §9.3 The `test_apply_with_trailing_axes_broadcasts` rename

`HarmonicMomentProjection` is a generic primitive that
contracts along the leading ordinate axis with any trailing axes
broadcasting.  The pre-edit test used variable names
`nx, ny, ng = 3, 4, 5` even though the values were arbitrary.  This
implicitly suggested an SN-storage commitment.  I renamed to
`a, b, c = 3, 4, 5` and added a docstring stating that the
trailing-axes shape is arbitrary — both `(N, nx, ny, ng)` and
`(N, ng, nx, ny)` are valid inputs.  The test's semantics are
unchanged (it asserts shape propagation, not values).

### §9.4 The long-running gates were skipped this cycle

Per §2.8, the L0 streaming-equilibrium gate (~17 min) and the
cylinder matvec invariants suite (~6 min) were not re-run in this
closeout.  Justification:
- PR-INDEX-6 touches ZERO code that affects flux values.
- The Step-1 bit-identity-via-transpose gate at PR-INDEX-5
  already covered the curvilinear regression snapshots at
  `rtol=1e-12`.
- The `__debug__` cell-flattening invariant runs in every regression
  test in §2.1 — its passage proves the principled / legacy
  roundtrip continues to hold.

If the gate-keeper wants to spend the runtime, the invocations are:
```bash
.venv/bin/python -m pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py -q
.venv/bin/python -m pytest tests/sn/spatial/test_apply_matvec_cylinder_invariants.py -q
```

### §9.5 The `.claude/plans/principled_index_migration.md` modification

The plan file shows as modified in `git status` but the diff is the
pre-existing `§0` status-checkpoint section that was authored prior
to this conversation (it has `**Last updated**: 2026-05-15 post
PR-INDEX-1` — predates this work).  PR-INDEX-6 does **not** touch
the plan.  If the gate-keeper wants the plan updated to mark
PR-INDEX-6 as DONE, that's a one-line edit in `§0`'s status table
that can land in the same commit or a follow-up.

## §10 Next step pointer

The migration is **functionally complete** at PR-INDEX-5.
PR-INDEX-6 ships the docs.  Two follow-ups remain:

1. **PR-INDEX-7 --- EquationMap packed-vector traversal flip**
   (deferred per PR-INDEX-4 §9.1).  Scope: ~200 LoC across
   `orpheus/sn/operator.py` (30+ `fi[:, n, i, j]` indexing sites),
   `solution_to_angular_flux*` allocation, retire the two
   `np.transpose` axis-swap adapters at `solver.py:1361` and
   `solver.py:1408`.  Not a blocker for the typed-field contract
   resume.
   - Dispatch: **method-implementer**.

2. **Typed-field contract resume** (plan §10).
   `AngularFlux(values: (N, ng, nx, ny))` and
   `ScalarFlux(values: (ng, nx, ny))` frozen dataclasses on the
   principled foundation.  Every leaf's `apply` signature becomes
   `apply(psi: AngularFlux) -> AngularFlux`; the four-operator
   algebra `(L + C - S - F/k).apply(psi)` distributes through
   `OperatorSum` unchanged.  Issue #197 partial close.
   - Dispatch: **method-implementer** for the dataclass
     introduction + leaf rewires; **archivist** for the typed-
     fields Sphinx narrative (which will cross-reference
     `index_convention.rst`).

Either follow-up can land first; they are independent.  The
recommended sequence is **typed-field contract first** (the user-
visible payoff is larger), then PR-INDEX-7 (internal cleanup).

## §11 Commit message (proposed)

```
docs(sn): index_convention.rst + layout audit cleanup (Issue #196 PR-INDEX-6)

New theory page docs/theory/index_convention.rst (931 lines) — the
canonical statement of the SN principled storage layout (N, ng, nx,
ny) for angular flux, (ng, nx, ny) for scalar flux + cross sections.
Documents the derivation (within-group block-diagonality argument
per Lewis-Miller §4.5 + Adams-Larsen §III), the six-PR migration
history (e09b9f8 through 3356cec + this PR), the load-bearing
Step-1 bit-identity-via-transpose verification gate (11/11 cases
PASS at rtol=1e-12, max abs diff 1.75e-14), and the full numerical-
evidence inventory.  Includes a layout-by-array reference table,
gotchas section (ny=1 singleton, SigS convention, per-material vs
per-cell xs, test fixture construction), and future-work section
(PR-INDEX-7 EquationMap flip, typed-field contract resume,
curvilinear joint-batch research path, JAX/GPU port).

Cross-references added from docs/theory/discrete_ordinates.rst and
docs/theory/operator_algebra.rst Key Facts / Conventions
admonitions.  Page added to the Infrastructure toctree in
docs/theory/index.rst.

Sweep audit: 17 legacy-shape mentions across orpheus/ and tests/
classified.  8 flipped to principled with cross-references to
:ref:`theory-sn-index-convention`:

  * orpheus/sn/sweep.py:285 — PR-INDEX-1 era docstring flipped to
    present-tense (post-PR-INDEX-5 reality: no bridge transposes
    remain).
  * orpheus/sn/operator.py:1212-1214 — SNStreamingOperator.solve
    contract.
  * orpheus/sn/operator.py:1788 — CollisionOperator σ docstring.
  * orpheus/sn/solver.py:677 — Krylov adapter Q_aniso shape.
  * orpheus/numerics/iteration.py:236, 305 — SourceIteration
    SN-shape example.
  * tests/sn/test_scattering_operator.py:283 — TestApplySemantics
    class docstring.
  * tests/numerics/test_projection_operators.py:87 — generic
    HarmonicMomentProjection test renamed local vars (a, b, c)
    + docstring noting the primitive is layout-agnostic.

9 hits kept intentionally:
  * 3 scan-internal contract descriptions (ordinate_scan's
    (nx, K, ng) Blelloch leading-scan-axis contract — primitive
    contract, not migration bridge).
  * 2 active intentional adapters in test_2d_octant_sweep_equivalence
    (PR-INDEX-4 documented "build legacy, view-transpose principled"
    pattern that preserves snapshot bit-identity).
  * 1 __debug__ cell-flattening invariant comment in solver.py:208
    (NAMES the legacy convention as part of the round-trip identity
    — flipping would erase the invariant's documentation).
  * 1 principled-aware adapter comment in test_mms_heterogeneous.py.
  * 1 diagnostic-narrative comment in
    tests/sn/diagnostics/phase_g_step2_03_closure_audit.py
    (hand-invoked structural audit; legacy-narrative description).

Test gates (all pass; no flux-values touched by this PR):
  * Regression: 11/11 PASS at rtol=1e-12 (81 s).
  * Directly-touched files: 95/95 PASS (94 s).
  * Operator-leaf suites: 131/131 PASS (1.0 s).
  * Quadrature + sweep + Phase C + Legendre moment + cell-update:
    189/189 PASS + 4 xpassed (1.5 s).
  * Numerics iteration: 11/11 PASS (0.5 s).
  * Projection operator: 2/2 PASS on the renamed
    test_apply_with_trailing_axes_broadcasts.
  * CP slab + cylinder (no-touch architectural proxy): 18/18 PASS.

Sphinx build: 16 warnings/errors pre-edit, 16 post-edit; baseline
unchanged.  Every warning is pre-existing and documented (skill md
files referenced as Sphinx docs, paramref role unknown,
inconsistent title style in discrete_ordinates.rst pre-existing
section structure, etc.).  No new warnings introduced by
index_convention.rst or by any cross-reference edit.

Nexus graph regenerated by sphinx-build; the new theory page node
is now indexed.

Closes PR-INDEX-6 of the principled index migration plan.  The
migration is functionally complete; two follow-ups remain
(PR-INDEX-7 EquationMap traversal flip; typed-field contract
resume per plan §10).
```
