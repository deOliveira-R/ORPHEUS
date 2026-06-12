---
name: face-name-carve docs (C4 / Issue #220)
description: Documenting a STORAGE-KEYING carve that collapses N hand-listed per-geometry branches to one crosswalk + one face-name-keyed dict, retiring named BC attributes (bc_xmin/bc_left/...) for a dict. Two-surface migration (INPUT-mesh bc_* KEEP vs SNMesh-resolved bc_* MIGRATE), structural-absence-vs-null sharpening, latent-dN-bug-closed-by-construction (NO ERR entry — type unconstructible), historical-code-block retirement-note (keep verbatim + note), absorbing a now-retired feature section into the carve record.
type: feedback
---

# Face-name carve docs (C4 / Issue #220)

The SN N-D layout campaign's C4 carve: a STORAGE-KEYING refactor.
`FaceLabel.face_name` crosswalk single-sources the `(axis, endpoint)
→ "{axis}{min|max}"` rendering; `boundary_face_layout` and
`SNMesh.bc` both collapse from 3-branch per-geometry hand-lists to
one comprehension over `face_labels`; named BC attributes
(`bc_xmin`/`bc_left`/`bc_ymin`...) retired for a face-name-keyed dict.

## The decisive doc moves (reusable for any storage-keying carve)

1. **Two-surface bc_* migration discipline — the load-bearing
   triage.** A grep for `bc_left`/`bc_xmin`/... hits TWO distinct
   surfaces that look identical:
   - **INPUT-mesh dataclass fields** (`Mesh1D.bc_left`,
     `Mesh2D.bc_xmin`, `bc_right=BC("vacuum")` constructor kwargs) —
     the user-facing geometry DECLARATION. **UNCHANGED → KEEP.** All
     of `docs/api/geometry.rst`'s hits and the
     `discrete_ordinates.rst` "Stage 1 — Geometry declaration" prose
     are this class.
   - **SNMesh RESOLVED-attribute spellings** (`sn_mesh.bc_left`,
     "stored as bc_xmin", `bc_right.apply`, `sn_mesh.bc_left ==
     "vacuum"`) — the solver's resolved BC surface. **MIGRATE** to
     `sn_mesh.bc["xmin"]` dict-keyed spellings.
   Disambiguate by reading 2 lines of context per hit: is it
   declaring (`Mesh1D carries...`, `bc_left=BC(...)`) or consuming a
   resolved op (`.apply`, `== "vacuum"`, "stored as")? Build a
   KEEP/MIGRATE table BEFORE editing.

2. **Curvilinear `bc_right`/`bc_left` → `bc["xmax"]`/`bc["xmin"]`,
   not `bc["right"]`.** A sphere/cylinder's outer radius keys
   `"xmax"` (the historical curvilinear convention:
   `_ENDPOINT_SUFFIX["outer"] = "max"`). So a Phase-D SPHERE test's
   `sn_mesh.bc_right.apply` migrates to `sn_mesh.bc["xmax"].apply` —
   add a parenthetical ("a sphere's `"outer"` endpoint renders as
   `"xmax"`") at the migration site so the reader isn't confused why
   a "right" boundary is keyed "xmax".

3. **Structural-absence-vs-null sharpening = a 3-column before/after
   table.** The carve's most subtle payoff: a pole-BC that was
   `bc_left = bc_xmin = None` (a named attr that EXISTS and HOLDS
   None) becomes a dict key that DOESN'T EXIST (KeyError on access).
   Document with a `.. list-table::` (Aspect × Pre-C4-None ×
   Post-C4-absence): the `bc` surface, the access result, the
   failure mode of a buggy consumer. This IS the
   illegal-states-unrepresentable (Pattern 4) load-bearing content.

4. **Latent d=N bug closed by construction → its own section + a
   "NO ERR entry" `.. note::`.** A hand-listed `"y" if face in (...)
   else "x"` axis dispatch is correct at d≤2 by string coincidence
   but a z-face would silently build the wrong reflection partner
   (vv Mode-9 class). The carve derives the axis from the label's
   own `AXIS_NAMES[axis_index]`, correct at any d. CRITICAL framing:
   because the triggering type (`Mesh3D`) is UNCONSTRUCTIBLE today,
   NO production bug ever shipped → NO ERR entry. The error catalog
   records SHIPPED L0-caught bugs; a defect closed-by-construction
   before its triggering type exists is documented in the theory
   page, NOT `error_catalog.md`. Cite the d=2 observable proxy test
   (with its non-vacuity guard) as the structural pin that makes the
   d=3 extension correct.

5. **Absorbing a now-retired feature section INTO the carve record.**
   The pre-existing `bc-1d-y-placeholder-design` section documented
   the exact y-placeholder C4 retires. Don't just delete it — the
   carve record's "What C4 retired" subsection ABSORBS it: the
   pre-C4 "why the placeholders were once safe" rationale moves into
   the retirement bullet ("...preserved in the
   :ref:`bc-curvilinear-realizer-unification` history above; C4
   removes the need for the rationale by removing the
   placeholders"). The anchor `bc-1d-y-placeholder-design` was
   removed → grep ALL of docs/+orpheus/+.claude/plans/ for the old
   anchor and repoint every LIVE `:ref:` (one in
   discrete_ordinates.rst:5719 → repointed to
   `bc-face-name-carve-what-retired`). The `docs/_build/html/` and
   `_nexus/graph.json` hits are stale ARTIFACTS — ignore.

6. **Historical code-block retirement-note (keep verbatim + note).**
   The Wave-2 "L7 trap" code blocks spell BOTH a retired 2-arg
   `.apply(psi, quad)` AND a retired `sn_mesh.bc_xmin` attribute.
   Keep the blocks verbatim (they document the L7-trap STRUCTURE,
   independent of spelling) + add ONE `.. note::` after them naming
   both retirements with forward `:ref:`s. Same pattern for the
   prose "this was `sn_mesh.bc_left == "vacuum"` pre-C4" inline
   historical flips.

7. **Section-vs-Closure ownership when inserting an H1 mid-arc.** I
   inserted a new H1 carve section between an existing H1
   ("curvilinear realizer unification") and its `Closure` H2. The
   `----` Closure then visually attached to MY H1. Fix: retitle the
   orphaned Closure to "Predecessor closure — ..." and add a short
   "C4 closure" H2 of my own ABOVE it, so each closure clearly owns
   its arc. Always check what H2/section follows your insertion
   point — a trailing same-level section silently re-parents.

## Cross-ref reality on this carve

- **`orpheus.sn.axis` is NOT automodule'd anywhere** → my
  `:class:`~orpheus.sn.axis.FaceLabel`` / `:func:`...face_labels``
  refs render PLAIN TEXT (degrade to `FaceLabel`/`face_labels` —
  clean, readable). `axis.py` has NO `:label:` docstrings so it
  WOULD be automodule-safe, but surfacing it is a separate
  architectural-docs task (would surface FaceLabel/AxisMesh/
  RadialAxisMesh/Axis1D + all pure shape fns while sibling sn leaves
  sweep_graph/sweep_schedule stay unsurfaced = inconsistent
  half-surfacing). DEFER + note the gap. `orpheus.sn.geometry` IS
  automodule'd → `:class:`SNMesh`` / `:attr:`SNMesh.bc`` /
  `:meth:`_resolve_one`` resolve fine.
- **My invented anchors were wrong twice — grep BEFORE citing.**
  `:ref:`Mode-9 <vv-mode-9>`` (no `.. _vv-mode-9:` defined →
  intra-doc dangling, `-W`-caught) → plain prose "vv-principles
  Mode-9 class". `:ref:`sn-curvilinear-pole`` (invented) → the real
  anchor is `sn-pole-angular-closure-protocol`. ALWAYS `grep
  "^\.\. _<anchor>:"` before a cross-doc `:ref:`.
- **Em-dash title underline: 1 code point, measure with
  `len(title)`.** My C4 title with `—` was 77 code points; I sized
  the underline 76 → "Title underline too short" (the ONE new
  warning, count 11→12). Fixed to 77. Scan ALL new section titles
  with a python `len()` loop, not `wc -c`.

## Build gate (this worktree)

Baseline = **11 warnings** (forced `-E`): 8 `_generated/*.rst`
InputErrors + 8 verification.rst include CRITICALs + 2 `.h5`
FileNotFound + 2 homogeneous plotting + 1 mesh.py paramref ERROR.
ALL env-artifact, NONE touch SN theory pages. Gate = post-edit
warning SET byte-identical to baseline (verified by line-normalized
`sed 's/:[0-9]+:/:L:/'` + sort + uniq) AND grep for
`boundary_conditions|discrete_ordinates|face-name|Title underline`
in the log returns EMPTY, EXIT=0. The Sphinx summary line ("build
succeeded, N warnings") is authoritative over a raw grep count (which
double-counts multi-line CRITICAL/InputError entries — I saw 21 from
grep vs the true 11 from the summary).

## Quality scores (this task)

| Dimension | Score | Note |
|---|---|---|
| Derivation depth | 5 | full face_labels→{layout,bc} chain + crosswalk code + structural-absence table + d=3 bug section |
| Cross-references | 4 | every SNMesh/realizer symbol linked; sn.axis refs plain-text (unsurfaced module, deferred) |
| Numerical evidence | 4 | structural carve — no numbers change; bit-identity-by-inheritance + sha256-golden + 1455-test green cited; L0 pin inventory enumerated |
| Failed approaches | 5 | string-vs-FaceLabel key decision + WHY; latent d=3 bug; degenerate-placeholder no-consumer retirement |
| Code traceability | 5 | every claim → :meth:/:attr: + the two consumer migrations named |
| Derivation source | n/a | structural/architecture, no SymPy derivation needed |
