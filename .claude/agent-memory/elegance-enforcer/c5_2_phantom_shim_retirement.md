---
name: c5-2-phantom-shim-retirement
description: C5.2 (#225) phantom-shim retirement review — ny/dy/dx RETIRED keep-nx asymmetry PASS-with-nits; twin-collapse dividend; stale mesh.ny docstrings in sibling transport layer
metadata:
  type: project
---

C5.2 (#225, commit `a28798d`, branch worktree-sn-nd-layout) — retires the
phantom-bearing legacy spatial metadata on SNMesh. Reviewed PASS-with-nits.
Sibling to [[c5-1-axis-primary-snmesh]] (C5.2 folds in C5.1's NIT 1 cross-ref
comment + the volumes docstring). The keystone-after-keystone before 3-D.

**The retirement cut (the load-bearing asymmetry, REINFORCE):** retire LIARS +
cheap DUPLICATES, keep HONEST high-traffic SUGAR. `SNMesh.ny`/`.dy` LIED at d=1
(phantom `1`/`[1.0]` — the #214 class) and underspecify at d≥3 → RETIRED
(AttributeError). `SNMesh.dx` was a pure duplicate of `axis_widths[0]` → RETIRED.
`SNMesh.nx` SURVIVES as documented `spatial_shape[0]` sugar — honest at any d AND
23 production reads dominated by the 1-D `iter_cells` direction logic
(geometry.py:1216-1378), 1-D bare sweep, sweep cache, pole closure. The asymmetry
is principled, NOT arbitrary: nx has a broad legitimate 1-D base, the others are
liars/dupes. `BulkField.nx/.ny` outright retired (an (nx,ny)-keyed field read
SILENTLY TRUNCATES a 3-D tensor — live d=3 landmine, correctly killed NOW).
`MaterialXSField`/`ScatteringOperator` nx/ny → ONE rank-generic `spatial_shape`
read-through each.

**The mechanism (consistent with C5.1):** `_axis_widths → public SNMesh.axis_widths`
is a stored INSTANCE ATTRIBUTE set in `_init_core`, consistent with its siblings
(`self.coord`, `self.nx`, `self.quad`, `self.reduced`, `self._volumes`) — SNMesh is
a plain class storing derived metadata as attrs, `@property` reserved for the
read-through facade (`spatial_shape`, `face_labels`, `volumes`, `volume_measure`).
axis_widths-as-public-attr is the established convention, not a slip. No consumer
mutates it.

**THE ELEGANCE DIVIDEND (the headline win):** the retirement let a TEST DELETE A
TWIN BRANCH. `test_unified_sweep_dispatch` went from explicit
`if sn_mesh.reduced is not None: (ng,nx) else (ng,nx,ny)` Σ_t builder → ONE
rank-generic `np.ones((ng, *spatial))`. Pattern-2 applied to the test itself,
enabled by the retirement. This is the "retire-pays-a-dividend" signature: a
retirement that COLLAPSES a path elsewhere, not just renames.

**The 3 phantom-PINNING test rewrites (all rank-honest, faithful migration):**
(1) krylov restart `N*ng*nx*ny` → `n_cells=int(np.prod(spatial_shape)); N*ng*n_cells`
— the MODEL rewrite (named intermediate `n_cells`, Pattern-3, correct at any rank;
the old form only worked because `ny==1` on the 1-D/curvilinear mesh). (2)
sweep-regression `sn_mesh.ny==1` assertion DELETED + `range(sn_mesh.ny)` →
`range(spatial_shape[1])` — BUT kept `mesh.dx[i]`/`mesh.dy[j]` (the **Mesh2D
DATACLASS**, a legitimate INDEPENDENT ORACLE, NOT the retired SNMesh.dx/.dy — the
test verifies the SNMesh streaming stencil against the dataclass widths; keeping
it is the RIGHT Pattern-2 exception). (3) slab metadata tuples rank-honest. Field
`nx`/`ny` property tests → RETIREMENT NEGATIVES (pin positive `spatial_shape==(3,5)`
AND negative `pytest.raises(AttributeError)` on `.nx`/`.ny` — inventory-pin
discipline). New C5-G13 `volume_measure` delegate-byte-identical pin. All
np.testing-only (Mode-8 -O-safe).

**volume_measure-as-pure-delegation is HONEST not premature.** Pure read-through
to the dataclass's own `volume_measure` SSOT while the adapter exists; the
`self.mesh is None` (d≥3 axis-native) arm EXPLICITLY deferred to C5.5 in the
docstring (not a hidden stub). `self.mesh` never None in present regime (line 249)
→ cannot NPE. SN-side consumers (solver.py:949,988) now read `sn_mesh.volume_measure`
not `sn_mesh.mesh.volume_measure` — closes the "mesh adapter is construction detail
not data path" seam + removes a `.mesh` reach-through from the keff hot path
(established benefit, defer-abstraction does NOT block this — 1 consumer pair, real
benefit). pole_angular_closure repointed `sn_mesh.mesh.coord → sn_mesh.coord` (C5.1
axis-derived).

**THE TWO NITS (both CONCERN-or-below, neither blocks):**
1. **(CONCERN, Pattern-7 / s6_4_f de-roling)** Six transport-layer docstrings name
the RETIRED attr `mesh.ny` in `from_mesh` derivation prose (`angular_flux.py:73`,
`scalar_flux.py:117`, `angular_residual.py:99`, `scalar_residual.py:78-79`,
`angular_source_sink.py:87`, `scalar_source_sink.py:98-99`). These PREDATE the
commit (last touched D-H.1a) but THIS commit makes them stale (pre-C5.2 `mesh.ny`
resolved via shim; post-C5.2 AttributeErrors). Production paths CLEAN (`from_mesh`
reads `spatial_shape`, verified) → CONCERN not VIOLATION (misleads a reader, breaks
no run). DISCRIMINATOR (per [[s6_9_fork_b2_default_flip]] stale-doc discipline):
literal `mesh.ny`/`mesh.nx` ATTRIBUTE refs = STALE → fix to `*mesh.spatial_shape`;
bare `(N,ng,nx,ny)` LAYOUT-axis-name placeholders = FINE, leave alone. The carve
that newly invalidates a ref owns sweeping it (de-roling discipline) — it swept
`_bases.py` but not the six siblings in the same layer.
2. **(cosmetic)** `_ref_compute_keff`/fission rank-1 mix `range(sn_mesh.nx)` (sugar)
+ `spatial_shape[1]` within one loop nest — a FAITHFUL consequence of the keep-nx
asymmetry (mirrors the prod surface), reads unevenly vs sibling
`(nx,ny),ng = ...spatial_shape, ...ng`. No bug; pure legibility nit.

**Migration idiom catalogue (uniformly high quality, for future N-D carve review):**
`*X.spatial_shape` splat (BEST — any-d), `(nx,ny),ng = X.spatial_shape, X.ng`
(nested unpack preserving ng), `nx,ny = X.spatial_shape` (clean d=2 unpack, self-
documents the 2-D assertion), `spatial_shape[1]` (single-axis), `np.prod(spatial_
shape)` (total cells). The `nx,ny = X.spatial_shape` form is GOOD not bad — it
ValueErrors at d=3 which is the honest 2-D-only assertion for a 2-D-specific test.
