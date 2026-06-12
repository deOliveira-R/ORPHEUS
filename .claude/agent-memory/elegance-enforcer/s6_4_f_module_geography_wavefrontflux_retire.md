---
name: s6-4-f-module-geography-wavefrontflux-retire
description: S6.4(f) #222 final sub-step — sweep.py dissolved into loss_representation + scan.py; WavefrontFlux/InteriorFaceSpace retired. PASS-clean pure-relocation review.
metadata:
  type: project
---

# S6.4(f) — module geography + WavefrontFlux retirement (issue #222, FINAL sub-step)

**Verdict: PASS (clean).** Pure relocation/retirement, zero arithmetic change. No CONCERNs
worth blocking. This closes the S6.4 (c)→(d)→(e)→(f) carve series ([[s6-4-dag-ownership-move]],
[[s6-4-e-walk-levelop-collapse]]).

**Why:** the previous (c)/(d)/(e) sub-steps each raised stale-docstring/role CONCERNs (3 stale
WavefrontFlux orchestrator docstrings at (d); 8 prod + ~11 test stale `:meth:` refs at (e)). (f)
SWEEPS those up as a side-effect of the relocation — the de-roling discipline (live→`:func:` new
home, dead→literal) is applied wholesale and is the cleanest of the five.

## What (f) did (all verified)
- `orpheus/sn/sweep.py` DELETED. Bodies relocated VERBATIM (byte-diffed against `git show
  HEAD:orpheus/sn/sweep.py`): `transport_sweep`/`_unwrap_source`/`_initial_guess_values`/
  `_sweep_1d_unified`(+`_ensure_geom_cache`/`_ensure_coll_cache`/`_run_1d_sweep`)/
  `_sweep_2d_wavefront` → `loss_representation.py` under an ORCHESTRATION banner (carries the
  Hébert/Bailey-Morel-Chang/Lewis-Miller/Blelloch refs). `_x_scan_faces`/`_scanmarch_row` →
  `spatial/scan.py`. RENAME `_sweep_2d_scheduled`→`_sweep_scheduled` (buffers went d-generic at
  (d), the `2d` lied).
- Only deltas vs originals: (1) the rename in defs + `:func:` docstring refs, (2) lazy-import
  drops (`from .loss_representation import default_for/_OctantWalk/_SolveOperands/_SweepEmit/
  MovingFrontierWindow` — all now in-module; the historical sweep↔loss_representation cycle GONE),
  (3) relocation banners, (4) 3 dropped section banners whose rationale RELOCATED to richer homes
  (the LD/EC/Step CellUpdate-Protocol open-problem note → `docs/theory/discrete_ordinates.rst`
  §2241-2266; the scan-march principled-equiv note → `ScanMarch` class docstring). No rationale lost.
- `WavefrontFlux` + `InteriorFaceSpace` git-rm'd (`transport/fields/wavefront_flux.py`,
  `numerics/spaces/interior_face_space.py`, `tests/transport/fields/test_wavefront_flux.py`).

## The load-bearing seam (spy-import re-points) — VERIFIED CORRECT
`operator.py` resolves `transport_sweep` via a **function-local** `from .loss_representation import
transport_sweep` (line ~2094); the 3 spies in `test_invertible_operator.py` patch
`"orpheus.sn.loss_representation.transport_sweep"` — patch target == call-time resolution site, so
the spy is observed. The `ordinate_scan` monkeypatch-by-assignment spy
(`test_ordinate_scan_joint_batch.py`) reassigns `loss_representation.ordinate_scan` (the
module-global imported at loss_representation.py:135); the 1-D sweep call sites (`_run_1d_sweep`
@1966) resolve through THAT global → observed. (The scan.py-internal `ordinate_scan` call @295 is
the 2-D scan-march path, not exercised by the slab-only spy — non-issue.)

## Recurring TELL this carve handled WELL (note for future relocations)
De-roling discipline on a module dissolve: a `:func:`orpheus.sn.sweep.X`` role survives a delete
only if X relocated. `_sweep_1d_cumprod`/`_solve_recurrence`/`_sweep_1d_spherical`/
`_sweep_1d_cylindrical` were LONG-DEAD pre-session symbols (count 0 in HEAD sweep.py) — correctly
de-roled to plain literals "(the dissolved `sweep.py`)" NOT re-roled (re-roling → broken Sphinx
ref). Live relocated symbols got proper `:func:`~orpheus.sn.loss_representation.X`` roles. Grep
gate to run on any module dissolve: `:func:`[~]?orpheus\.<mod>\.` with TRAILING DOT (catches the
member-role form, skips the bare `:mod:` survivors of sibling modules like sweep_graph/
sweep_schedule/sweep_scratch).

## Retirement succession notes — ADEQUATE (every named symbol verified real)
`transport/fields/__init__.py` + `sweep_graph.py` `_FrontierPlan` docstring carry the succession:
the C¹_int concept lives on as `_MovingFrontier` (sweep_graph.py:252, rolling (d−1)-frontier
seed/shed) + `FullFieldWavefront._octant_face_cochain` (loss_representation.py:985) + `_edge_outflow`
(1014); whole-trace exchange = `_OctantWalk` (450). All four symbols grep-confirmed present. The
retitled theory `wavefront-flux-cochain` section (parallel archivist) holds the history.

## loss_representation.py size (review point 2) — NOT a concern
~2394 lines = selector (Compatibility/IncompatibleRepresentation/LossRepresentation/default_for)
+ 4 representations (CumprodScan/MovingFrontierWindow/FullFieldWavefront/ScanMarch) + walk frames
(_OctantWalk/_DAGWavefront) + ORCHESTRATION (transport_sweep + schedule loop + 1-D body). USER
explicitly approved "sweep.py disappears since its functions have been entirely taken by
loss_representation.py". `__all__` complete+honest (8 public symbols + LOSS_REPRESENTATIONS;
nothing private leaks; transport_sweep annotated "relocated from the dissolved sweep.py").

Import smoke test green (worktree on PYTHONPATH, main venv): all relocated symbols resolve,
`_sweep_2d_scheduled` retired, `orpheus.sn.sweep`→ModuleNotFoundError, both retired types absent
from package surfaces.
