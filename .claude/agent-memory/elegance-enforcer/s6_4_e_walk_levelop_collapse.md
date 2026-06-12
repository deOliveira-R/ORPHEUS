---
name: s6-4-e-walk-levelop-collapse
description: S6.4(e) review — 4 graph walks (apply/residual ×full/windowed) → 2 storage walks × _CellSolve/_CellResidual level-op; CONCERNS (no FAIL)
metadata:
  type: project
---

# S6.4(e): the graph direction×storage product collapse — REVIEW VERDICT

The fifth sub-step of issue #222's S6 re-layering (branch `worktree-sn-nd-layout`,
HEAD `3438131` at review). Reviewed uncommitted working-tree changes 2026-06-11.

**Verdict: PASS-with-CONCERNS (no FAIL, no VIOLATION).** The carve is correct and
elegant at the load-bearing axes; the concerns are documentation-correctness +
one illegal-states-representable gap, not structural bugs.

**Why:** the architectural core — collapsing the `apply`/`residual`×`full`/`windowed`
2×2 method product into `walk_full`/`walk_windowed` (storage) × `_CellSolve`/
`_CellResidual` (direction-as-OBJECT, frozen dataclass with one `.cell` method) — is
exactly the right shape. Direction is never a boolean flag (Pattern 3 satisfied);
the storage walk owns the level loop + cochain, the level op owns cell-algebra
dispatch + per-level emit; the emit expressions relocated VERBATIM (einsum strings +
left-fold order preserved → bit-identity intact, pinned by window≡full
`assert_array_equal` + sha256 kernel-source pin). `diamond.py` retired to pure cell
algebra (storage adapters DELETED with honest in-file note); `SweepCellSlice` retired
from `cell_update.py`; the kernel-pair extension point (`cell_kernel_batch`/
`residual_kernel_batch` NotImplementedError stubs) is complete + honest for future
Step/LD/EC. Verified 155 green under `-O` (kernel + one-octant-walk + 3 graph files).

**How to apply (the 3 concerns, for the next reviewer / a follow-up):**

1. **CONCERN — `_CellSolve` is missing the `_SweepEmit.__post_init__` mode guard
   (illegal-states-representable).** `_CellSolve` carries the SAME angular-XOR-moment
   DI fields as `_SweepEmit` (`weights`, `angular_flux`, `scalar_flux`, `moment_buf`,
   `Y`) and the SAME `if moment_buf is None:` branch in `.cell`, but has NO
   `__post_init__`. SMOKE-VERIFIED two illegal builds succeed: (a) all output buffers
   `None`, (b) `moment_buf` set + `Y_octant=None`. Both fail only DEEP in the walk
   (`None.__iadd__` AttributeError / einsum on None) far from the construction site.
   `_SweepEmit` one level up rejects both at `__post_init__`. Production is safe TODAY
   (the 2 prod `_CellSolve` sites derive their mode from the already-guarded `emit`),
   but the type is the override/extension point future Step/LD kernels construct
   directly. **Bug habitat:** a future direct constructor (or a test) half-wires the
   mode and gets an opaque deep-walk crash. **Required:** add the exactly-one-mode
   `__post_init__` to `_CellSolve` (mirror `_SweepEmit` lines 424-431) so the illegal
   half-wired output is unrepresentable at construction. Cheap, local, closes the gap.

2. **CONCERN (Cardinal Rule 3 doc-correctness) — stale `:meth:` refs to the 4 RETIRED
   graph methods survive in `.py` docstrings/comments** even though the archivist
   already repointed them in `operator_algebra.rst` (clean). 8 production sites:
   `loss_representation.py` lines 10, 972, 981, 1069, 1105, 1126; `sweep.py:855`;
   `operator.py:101` — all `SweepDependencyGraph.apply`/`.residual`. Plus test
   docstrings/comments: `test_sweep_graph.py:297,360,392`,
   `test_sweep_graph_nd_admission.py:390,391,416,450`,
   `test_sweep_graph_window_equivalence.py:4,5,218,240` (`apply_windowed`/
   `residual_windowed`/`update_batch`/`graph.apply`). `conf.py` has NO `nitpicky`, so
   these warn-not-fail — but they are dangling xrefs to deleted methods. **Required:**
   repoint to `walk_full`/`walk_windowed` (× the level op where the direction matters).

3. **NOTE (not a finding) — the moment-mode `if moment_buf is None:` branch now appears
   in THREE interior emit sites:** `_CellSolve.cell` (per-level, sweep_graph.py),
   `ScanMarch._sweep_interior` (per-row), `_SweepEmit.pure_z` (whole-octant). This is
   institutional-memory pattern #2 ("twin emit appears N times, verified identical,
   acceptable-for-now") — ACCEPTABLE: each is the SAME einsum at a different granularity
   (the granularity IS the kernel's distinguishing concern), value-identical, pinned by
   the window≡full + moment oracles. STOPS being acceptable at a FOURTH site or an edit
   landing on one-not-all. Collapse destination if it grows: a shared per-granularity
   emit primitive the three reference.

**Test-migration audit (test_cell_update_batch.py → test_cell_kernel_batch.py):
FAITHFUL.** Old 9 tests mapped: 6 preserved (closed-form, bit-id→now per-ordinate-loop
ref, isotropic-Q, 2× default-raises, residual closed/affine, round-trip), 2 DISSOLVED
as documented (negative-octant variant = gather INDICES = storage; residual-without-probe
= now signature-enforced via required `psi_bar` kwarg), `TestFaceFluxScatter` (2)
genuinely MOVED to `test_sweep_graph.py` (verified: walk drives `walk_full`, scatter
checked vs per-cell loop, all 4 octants). Two NEW gates added: sha256 kernel-source pin
(the genuine FP-reduction-tree-of-record exception — left-fold non-associativity) +
`test_graph_exposes_two_walks_not_four` retirement assertion. Nothing dropped beyond the
2 documented dissolutions.

**`walk_full` inline gather/scatter vs `FullFieldWavefront._octant_face_cochain`
addressing — NOT a duplication VIOLATION.** `walk_full` inlines the per-axis advanced
-index gather/scatter (relocated VERBATIM from the retired DiamondDifference adapters via
the relocated module-fn `_cell_face_selector`); `_octant_face_cochain`/`_edge_outflow` in
loss_representation handle the octant-buffer ALLOCATION + in-edge SEED + out-edge EXTRACT
(a different concern — boundary embedding, not interior face addressing). They share the
`2 + axis` lattice-position idiom but address different things (interior level faces vs
domain edge slots). Bounded, single-purpose each. No collapse demanded.
