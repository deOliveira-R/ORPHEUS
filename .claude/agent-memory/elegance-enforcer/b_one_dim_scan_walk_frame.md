---
name: b-one-dim-scan-walk-frame
description: #206 Phase B — extract _OneDimScanWalk (1-D analogue of _OctantWalk) frozen mesh-holder frame; PASS-with-nits; relocation bit-faithful but 3 dangling :func: rst xrefs + PEP8 method spacing
metadata:
  type: project
---

# #206 Phase B — `_OneDimScanWalk` frame extraction (the 1-D analogue of `_OctantWalk`)

PASS-WITH-NITS (no commit-blockers). The SOLVE-only frame extraction following the A1/ScanMarch-2D
face-closure single-sourcing ([[a1_diamond_face_closure_seam]], [[scanmarch_2d_diamond_closure_routing]]).
Working-tree diff; bit-id gates ran separately (per dispatch constraint — I did `git diff`/`grep`/`Read` only).

**What B did:** extracted `@dataclass(frozen=True) class _OneDimScanWalk` (sole field `mesh: SNMesh`,
no `__init__`) and relocated the former module-level free helpers `_sweep_1d_unified`(→`.sweep`),
`_ensure_geom_cache`/`_ensure_coll_cache`(→methods), `_run_1d_sweep`(→`._run`) INTO it. DETERMINISTIC
transform: re-indent +4, whole-word `sn_mesh`→`self.mesh`, the 3 `.sweep` cross-calls→self-methods;
the ~380-line `_run` body otherwise VERBATIM → bit-identical. 2 prod callers (`CumprodScan.sweep:726`,
`ScanMarch.sweep` 1-D branch `:1246`) route through `_OneDimScanWalk(self.mesh).sweep(...)`.
`_initial_guess_values` (pure container-extractor, no mesh) KEPT module-level (`:1675`). B2 gate
`TestOneDimScanWalkFrame` (frozen+sole-`mesh`-field + the AST `is_solve`/`is_apply`/`is_matvec` tripwire).

**AXIS-1 — abstraction carries its weight NOW (PASS).** NOT premature Pattern-6 scaffolding. There are
TWO live consumers TODAY (`CumprodScan` + `ScanMarch`-1D both call `.sweep`) — unify-after-two is SATISFIED,
exactly as `_OctantWalk` earns its keep with 2 consumers. The frame is the genuine 1-D twin of the landed
`_OctantWalk` (both `@dataclass(frozen=True)`, sole field `mesh:"SNMesh"`, own their sweep body — verified
`:444-445` vs `:1719-1735`). Forwarding the matvec to Phase C is the CORRECT defer (no apply-kernel param
built ahead of its consumer — Pattern-6 respected; the injection point is DOCUMENTED in the class docstring
`:1730-1732` "MATVEC attaches in Phase C as the per-ordinate apply-kernel, mirroring `_OctantWalk`'s
cell-kernel injection", NOT pre-built). The frame establishes the home that Phase C's apply-kernel needs;
extracting it as a relocation-only step (zero behaviour change, bit-id) is the textbook clean-before-extend
move ([[clean-before-extend]]). Does NOT need to wait for Phase C.

**AXIS-2 — relocation fidelity (PASS, exhaustive).** Whole-word transform COMPLETE + correct: `awk` over
the class body (1712-2005) = ZERO surviving `sn_mesh` inside the frame; all surviving `sn_mesh` in the file
are OUTSIDE (`_compute_decomposition` 617-678, `transport_sweep`/`_sweep_scheduled` 1519+, `_sweep_2d_*`
2220+, module docstring 17-18). No `sn_mesh` in a docstring/string was wrongly flipped (the only old-name
mentions inside the frame are the relocation-provenance note `:1726-1727` naming the retired helpers — correct).
Pre-existing local-vs-attribute `cell_update` inconsistency FAITHFULLY preserved (local `cell_update =
self.mesh.cell_update` `:1885` used by `.update()` `:2144`, but `.cell_average_from_faces` reached via
`self.mesh.cell_update.` directly `:2019,2184`) — pre-A1, NOT introduced; correct to leave (fixing risks
the bit-id gate; Phase-C-cleanup candidate). Not an alias smell (whole-word, not `sn_mesh = self.mesh` alias).

**AXIS-3 — `_initial_guess_values` kept module-level (PASS).** Correct call. It is a pure container-agnostic
extractor (duck-types `.bulk`/`.values`, no mesh, no self-state); making it a `@staticmethod` on the frame
would be FALSE cohesion (it is consumed by `_run` AND referenced by `transport_sweep`'s docstring `:1604`;
it is a free utility of the iterate-container, not of the 1-D walk). Mirrors the design directive. Leave it.

**AXIS-4 — geometry (`is_slab`) fork, SOLVE-only (PASS).** Appropriate for Phase B. The `is_slab` fork lives
INSIDE `_run` (the 1-D scan's only variation IS geometry); `.sweep` is the solve entry; NO apply-kernel
param exists yet (correct — Phase C). Clean seam: the frame has exactly 4 methods, zero premature
apply/matvec scaffolding (grep-verified — `apply`/`matvec` appear ONLY in comments: the bare-sweep note +
Pattern-2 cross-refs to the matvec's shared Carlson seed). Pattern-2 routing to the SHARED `psi_half_seed`
strategy preserved (`:2050,2122` — sweep/matvec stay ONE discrete system, not re-inlined).

**AXIS-5 — B2 tripwire + frozen mesh-holder assertion (PASS, sufficient).** `TestOneDimScanWalkFrame` is a
faithful clone of the established `test_octant_walk_is_kernel_parameterized_not_boolean`
(`tests/sn/operators/test_one_octant_walk.py:84-119`): identical `smells={"is_solve","is_apply","is_matvec"}`,
same AST walk over Name/Attribute/arg/keyword identifiers, `-O`-safe `pytest.fail`. `test_frame_resolves_as_
frozen_mesh_holder` pins frozen + `fields==["mesh"]` — the structural mirror of `_OctantWalk`. Correct
+ sufficient for Phase B; it is the anti-degradation guard that stops the Phase-C matvec attachment from
becoming a bool flag (the exact 1-D mirror of the octant tripwire).

**AXIS-6 — no elegance regression (PASS) bar the nits below.** No alias smell. `:func:`→`:meth:` updates
inside the 3 touched .py files (operator.py:1116, scan.py:41, sweep_cache.py:76 + loss_representation
internal docstrings 8/696/1632) are accurate and point at the real new home.

## NITS (non-blocking)

**NIT-1 (CONCERN — doc drift, Cardinal Rule 3 / anti-pattern #11): 3 LIVE `:func:` rst xrefs to the
RETIRED `_sweep_1d_unified` were MISSED.** The diff updated the `:func:`→`:meth:` refs in the 3 prod .py
files but NOT the theory docs. SURVIVING dangling cross-reference ROLES:
- `docs/theory/loss_representations.rst:273` — `:func:`~orpheus.sn.loss_representation._sweep_1d_unified``
- `docs/theory/index_convention.rst:1176` — `:func:`_sweep_1d_unified``
- `docs/theory/index_convention.rst:1422` — `:func:`_sweep_1d_unified``
These are LIVE Sphinx roles (not historical literals) → now name a function that no longer exists. The
function was never `autofunction::`'d so they likely already soft-dangled (Sphinx warns, builds unless
nitpicky — `conf.py` has no `nitpicky`), but the carve makes them SEMANTICALLY stale. Fix → `:meth:`~orpheus.
sn.loss_representation._OneDimScanWalk.sweep``. DISCRIMINATOR (carried from [[c5_2_phantom_shim_retirement]]/
[[c5_3_geometry_blind_trace]]): a `:func:` ROLE that dangles = FIX; a literal double-backtick mention in
HISTORICAL prose = LEAVE. The LEAVE bucket (do NOT touch): `index_convention.rst:405` (PR-INDEX changelog
history row, frozen record), `tests/.../test_boundary_face_layout.py:18`, `test_sweep_graph_nd_admission.py:20`,
`test_unified_sweep_dispatch.py:403` (plain-text thing-under-test docstrings, non-breaking), `test_streaming_
operator.py:1095` (historical T.5-scope prose). Recommend sweeping NIT-1 in THIS carve (no-issues-for-inline-
fixes; the A/B series owns the #206 doc surface) — it is a 3-line edit, archivist-class.

**NIT-2 (cosmetic — PEP8): ZERO blank lines between the 4 frame methods.** `sweep`→`_ensure_geom_cache`
(`:1797`→`:1798`), `_ensure_geom_cache`→`_ensure_coll_cache`, `_ensure_coll_cache`→`_run` all abut with no
blank line. PEP8 + the template want ONE (`_OctantWalk` `:524`→`:526` has it). Artifact of the mechanical
re-indent (the free fns had 2 blank lines; collapsed to 0). 3-line fix; not a correctness/structure issue.

Files: `loss_representation.py` (1719-2210 the frame; 726/1246 callers; 1675 `_initial_guess_values` kept;
8/696/1632 internal docstring refs updated), `operator.py:1116`, `scan.py:41`, `sweep_cache.py:76`
(`:func:`→`:meth:` updated), `tests/sn/sweep/core/test_unified_sweep_dispatch.py:396-457` (B2 gate).
MISSED docs: `loss_representations.rst:273`, `index_convention.rst:1176,1422`.

Related: [[a1_diamond_face_closure_seam]] + [[scanmarch_2d_diamond_closure_routing]] (the #206 A-series
face-closure single-sourcing this builds on), [[s6_4_c_dag_ownership_move]] (the `_OctantWalk` frame this is
the 1-D analogue of), [[clean-before-extend]] (the relocation-before-Phase-C-extension discipline).
