---
name: Issue #168 Phase C — sweep-frame matvec rewrite + BoundaryFaceFlux retirement + empirical default-flip deferral
description: Phase C of Issue #168 ships on `refactor/sn-operator-algebra` 2026-05-12. Sweep-frame apply matvec rewrite — WDD spatial closure + BC trace law at boundary edge + retired Phase A BoundaryFaceFlux Protocol. Empirical Gate 1.1 outcome → MorelMontryAngularSweep default flip DEFERRED to Phase D. ERR-026 stays at PARTIAL CLOSURE (spatial-closure aligned; pole-face refinement Phase D scope).
type: project
---

# Issue #168 Phase C — closeout

`refactor/sn-operator-algebra` 2026-05-12. The architectural alignment
that Phase B identified as the load-bearing precondition for the
curvilinear-default flip. Resumed per `.claude/plans/issue_168_phase_c.md`
after two post-pause unblockers (trajectory_resolvent cylinder MR
Phase 1b coverage + cross-domain-attacker's identification of the
§16A.3 affine-BC-contract architectural inconsistency).

## Why this memo exists

Future implementers (Phase D) need the Phase C surface area pinned:
what shipped, what the empirical Gate 1.1 outcome was, and which
follow-up the deferred default flip depends on. **Why it matters:**
ERR-026 stays at PARTIAL CLOSURE through Phase C. The 4 xfail-strict
markers on curvilinear MMS tests stay xfail. Phase D is the closing
fix and depends on:

- The Phase C architectural alignment (sweep-frame matvec + BC trace
  contract honoured + retired BoundaryFaceFlux Protocol) which IS
  the load-bearing infrastructure;
- A new pole-face spatial closure refinement that makes
  `MorelMontryAngularSweep` work cleanly with the WDD recurrence at
  the spherical pole (the failure mode Gate 1.1 surfaced).

**How to apply:** Before scoping Phase D, reread this memo + the
Phase B closeout. The architectural changes Phase C ships are
final — Phase D extends them with one additional pole-face closure
strategy, NOT a rewrite.

## What shipped

### Commits (5 atomic per the plan §3 sequencing)

1. `eae6f05` — `feat(sn): iter_cells_by_direction + unknowns_at_cell_for_mask helpers (Issue #168 Phase C)`
2. `16e33d7` — `test(sn): Phase C gate skeleton — linearity, reciprocity, flat-flux probe (Issue #168 Phase C)`
3. `3fd1302` — `feat(sn): sweep-frame apply matvec + retire BoundaryFaceFlux Protocol (Issue #168 Phase C)`
4. `3a95e3b` — `test(sn): Phase C MMS + cross-check gates + empirical xfail markers (Issue #168 Phase C)`
5. `d445a8f` — `docs(sn,bc): Phase C Sphinx + ERR-026 closure narrative (Issue #168 Phase C)`

Plus `7497cec` — `test(sn): regenerate 6th curvilinear snapshot + sn_reshape Wave H Phase C row`.

### New APIs

* **`SNMesh.iter_cells_by_direction(direction_sign[, mu_level_idx])`**
  — companion to `iter_cell_visits` that yields cells in DAG-topological
  order keyed by the bulk sweep direction sign. Vectorised across the
  ordinate subsets (`outgoing_mask = quad.mu_x > +eps` /
  `incoming_mask = quad.mu_x < -eps`) without needing per-ordinate
  iteration. Cylindrical case takes a per-level `mu_level_idx`.
  Foundation tests pin bit-identity to `iter_cell_visits(ordinate_idx=
  representative_n)` across sphere / slab / cylindrical for every
  non-degenerate ordinate.
* **`EquationMap.unknowns_at_cell_for_mask(cell_idx, ordinate_mask)`**
  — precomputed inverse lookup `(cell, ordinate) → k`. Lazily builds
  an `(nx, N) int` table with `-1` sentinels; O(1) call-time per
  (cell, mask) combination.

### Sweep-frame matvec rewrite

`transport_operator_matvec_spherical` and `_cylindrical` rewritten
per the plan §2 pseudocode:

1. **Outward sweep** (μ ≥ 0): for each cell in `iter_cells_by_direction(+1)`,
   compute WDD diamond `ψ_face_out = 2·ψ_cell − ψ_face_in`, accumulate
   streaming + redistribution + collision per ordinate subset.
2. **BC trace law applied ONCE** at the boundary edge:
   `bc_outer.apply(outflow_at_boundary)` consumes the WDD-propagated
   outflow face vector and returns the inflow face vector for the
   inward sweep.
3. **Inward sweep** (μ < 0): for each cell in `iter_cells_by_direction(-1)`,
   same WDD recurrence consuming the BC-determined `psi_face_in`.

**Pole-face initial condition** (the deviation from the plan's `0`
pseudocode): `ψ_face_in(pole) = ψ_cell[0]` (Lewis-Miller §4.5). For
true flat ψ this gives `ψ_face_out = 2·ψ_cell − ψ_cell = ψ_cell` at
every cell, preserving the per-ordinate flat-flux invariant. Naive
`ψ_face_in(pole) = 0` (per the plan's pseudocode) caused oscillating
face fluxes `0, 2c, 0, 2c, …` that broke flat-flux balance on all
three pole-closure strategies.

**Vectorisation across ordinates**: per user constraint 1, the sweep
body has NO `for n in range(quad.N)` loop. Cell-level streaming is
computed as `mu_out[None, :] * (A[i+1] * psi_face_out - A[i] *
psi_face_in) / V[i]` on the (ng, n_out) subset; the cylindrical
case adds an outer per-level loop (the level-DAG topology is
intrinsic to cylindrical's azimuthal-level structure).

### `solution_to_angular_flux_*` simplification

Return value changed from `(fi, boundary_face_flux)` (Phase A) to
just `fi`. The companion `boundary_face_flux` array retired with
the Protocol. The inward-at-outer-boundary cell-centre slots
`fi[:, n_inward, -1, 0]` are now filled with the **reflected-partner
cell-centre value** as an analytical extension — the eq_map excludes
these from unknowns (BC determines them), but the WDD recurrence on
flat ψ requires the cell-centre to be consistent so the per-ordinate
flat-flux invariant holds.

### BoundaryFaceFlux Protocol RETIRED

Phase A's Protocol + concrete strategies + tests + SNMesh field
RETIRED entirely:

- `orpheus/sn/spatial/boundary_face_flux.py` DELETED (415 LOC)
- `tests/sn/spatial/test_boundary_face_flux.py` DELETED (21 foundation
  tests)
- `SNMesh.boundary_face_flux` constructor field + attribute REMOVED
- `boundary_face_flux_closure` kwarg REMOVED from
  `transport_operator_matvec_spherical` and `_cylindrical`
- `orpheus/sn/spatial/__init__.py` exports cleaned up

The boundary-face closure is now INSIDE the sweep-frame WDD
propagation chain — the BC trace law owns the boundary edge per
§16A.3.

### Test updates

`tests/sn/test_snstreamingoperator.py` — 5 Phase-A-specific tests
retired or rewritten (full diff in commit `3fd1302`):

* `test_solution_to_angular_flux_spherical_returns_tuple` →
  `test_solution_to_angular_flux_spherical_returns_single_array`
  (pins Phase C single-fi return).
* `test_solution_to_angular_flux_spherical_fi_is_pure_cell_center_vacuum`
  → `test_solution_to_angular_flux_spherical_inward_extension`
  (pins reflected-partner fill at i=N-1).
* `test_apply_spherical_outer_face_uses_dd_extrapolation` →
  `test_apply_spherical_wdd_recurrence_deterministic` (bit-stable
  via `np.array_equal`).
* `test_apply_spherical_constant_flux_yields_zero_collisionless` →
  `test_apply_spherical_constant_flux_yields_zero_collisionless_reflective`
  (REFLECTIVE BC variant — Phase C invariant holds; vacuum BC no
  longer preserves flat ψ structurally).
* `test_apply_spherical_vacuum_bc_constant_flux_no_corruption` →
  `test_apply_spherical_vacuum_bc_residual_is_consistent` (linearity
  + BC physically removes inflow contract).
* `test_apply_cylindrical_outer_face_uses_dd_extrapolation` →
  `test_apply_cylindrical_wdd_recurrence_deterministic`.
* `test_snmesh_accepts_boundary_face_flux_kwarg` →
  `test_snmesh_no_longer_accepts_boundary_face_flux_kwarg` (field
  retirement contract pinned).
* `test_apply_spherical_constant_flux_under_morel_montry_canonical_form`
  rewritten to use a reflective sphere (angular-integrated invariant
  preserved by α-telescoping at interior cells under MMS + reflective
  BC).

### New test files

* `tests/sn/test_iter_cells_by_direction.py` — 11 foundation tests
  for the new APIs.
* `tests/sn/test_phase_c_gates.py` — Gates 1.1 (parametrised over 3
  pole closures × 2 Σ_t × 2 geometries), 1.2 (apply determinism),
  1.3 (apply ↔ apply_transpose reciprocity), 1.4 (linearity, the
  precondition), 1.5 (BC trace contract).
* `tests/sn/test_phase_c_mms.py` — Gates 3.1, 3.2, 3.3 (MMS spatial
  + angular convergence; 3.1 + 3.2 marked xfail).
* `tests/sn/test_phase_c_crosscheck.py` — Gate 4.1 (k_∞ recovery,
  PASSES) + Gate 4.2 (trajectory_resolvent crosscheck, SKIP for
  Phase D).

### Sphinx narrative

`docs/theory/discrete_ordinates.rst` gains:

* Phase A subsection rewritten to note the Protocol retirement.
* New Phase C subsection with 6 `:label:` anchors:
  - `phase-c-sweep-frame-matvec` (introduction)
  - `phase-c-apply-sweep-equivalence`
  - `bc-trace-contract-respected-by-matvec`
  - `sn-mms-spherical-aniso-spatial-convergence`
  - `sn-mms-cylindrical-aniso-spatial-convergence`
  - `sn-curvilinear-homogeneous-kinf-recovery`
  - `sn-curvilinear-trajectory-resolvent-crosscheck`

`docs/theory/boundary_conditions.rst` §16A.3 gains a note pinning
the Phase C architectural fix: the SN apply matvec now honours the
affine-BC-trace contract.

`.claude/skills/vv-principles/error_catalog.md` ERR-026 entry
extended with "What Wave H Phase C added" subsection.

### Regression snapshots

6 curvilinear snapshots regenerated under the Phase C operator:

- `cyl_1g_homogeneous_LS4_dd_n20`
- `cyl_1g_homogeneous_product_dd_n20`
- `cyl_2g_3reg_LS4_dd_n40`
- `sphere_2g_3reg_dd_n40`
- `sphere_2g_homogeneous_dd_n20`
- `sphere_2g_p1_aniso_dd_n20`

5 Cartesian snapshots stay bit-identical (Gate 2.1 PASS).

## Architectural deviations from the plan

### Deviation 1: pole-face initial condition is `psi_cell[0]`, not 0

The plan §2 pseudocode reads `psi_face_in = np.zeros((ng, n_out))`
with the comment "Pole face: ψ_face = 0 by symmetry (also multiplied
by A[0]=0)". Empirically, this initial condition combined with the
WDD recurrence on flat ψ produces oscillating face fluxes
`0, 2c, 0, 2c, …` that break the per-ordinate flat-flux invariant on
**all three pole-closure strategies** (not just MMS).

The shipped code uses `psi_face_in(pole) = psi_cell[0]` (Lewis-
Miller §4.5 "Carlson starting direction"). For true flat ψ this
yields `psi_face_out = 2·psi_cell - psi_cell = psi_cell` at every
cell, preserving the invariant. Cylindrical analogue identical.

This is documented in code comments + the closeout. The plan's
pseudocode was a simplification; the working code's pole-face
initialisation is the Carlson seed standard.

### Deviation 2: empirical Gate 1.1 outcome → default flip DEFERRED

Per the user's explicit constraint 7, the empirical Gate 1.1 probe
with `MorelMontryAngularSweep` determines the default flip. Results
on sphere + cylinder with all 4 (Σ_t × pole-closure) combinations:

| Geometry | Pole closure | Σ_t = 0 | Σ_t = 0.5 |
|----------|--------------|---------|-----------|
| Sphere   | Legacy       | PASS    | PASS      |
| Sphere   | BFF          | PASS    | PASS      |
| Sphere   | MMS          | **FAIL**| **FAIL**  |
| Cyl      | Legacy       | PASS    | PASS      |
| Cyl      | BFF          | PASS    | PASS      |
| Cyl      | MMS          | PASS    | PASS      |

Sphere MMS fails; cylindrical MMS passes. The asymmetry is the
spherical pole-face WDD initial condition interacting with the
canonical Hébert §3.9.4 angular recurrence (the recurrence's
half-angle face flux at the pole is sensitive to the initial
condition in a way the cylindrical per-level α-telescoping
absorbs). Per user constraint 7, **default flip DEFERRED to Phase D**.

Cylindrical MMS Gate 1.1 PASS is a strong signal for the eventual
Phase D fix — the cylindrical structure has the per-level α-dome
telescoping that the spherical case lacks at the pole. Phase D's
pole-face spatial closure refinement is the architectural addition
needed.

### Deviation 3: Gate 4.2 placeholder, not full implementation

The plan §5 Gate 4.2 calls for trajectory_resolvent cross-check on
5 snapshots with rtol ≤ 1e-9 (relaxable to 5e-4 if SN nx-
discretisation dominates). Phase C ships a SKIPPED placeholder
because:

1. The cross-check's premise (SN flux shape converges to the
   continuous Green's-function reference) requires ERR-026 to be
   CLOSED on flux shape — which Phase C explicitly does not deliver
   (default flip deferred).
2. The full agreement test is structurally complex (interpolation
   onto SN cell centres + normalisation matching), which is its own
   work item separate from the architectural alignment.

The placeholder pins the structurally-independent reference for
Phase D: the bare entry points
`solve_greens_function_{cylinder,sphere_mr,cylinder_mr}` documented
in the test docstring.

### Deviation 4: ERR-026 xfail-strict markers STAY xfail

The plan §3 step 11 says "Markers REMOVED (not relaxed to
strict=False)" — but this was conditional on Gate 1.1 passing. Since
Gate 1.1 failed under MMS on sphere, the underlying ERR-026 flux-
shape bug remains and the four xfail-strict markers correctly stay
xfail. Their reason strings stay informative for Phase D.

## Verification gates — summary

| Gate | Description | Status |
|------|-------------|--------|
| 1.1  | Per-ordinate flat-flux residual | PASS for Legacy + BFF on both geometries; xfail(strict=False) for MMS on sphere; PASS for MMS on cylinder |
| 1.2  | apply determinism via np.array_equal | PASS |
| 1.3  | apply ↔ apply_transpose reciprocity to rtol=1e-12 | PASS |
| 1.4  | apply linearity to rtol=1e-13 (PRECONDITION) | PASS |
| 1.5  | BC trace contract — apply(0)=0 under vacuum + reflective | PASS |
| 2.1  | 5 Cartesian regression snapshots bit-identical (rtol=1e-12) | PASS |
| 2.2  | Phase B 28 foundation tests | PASS |
| 2.3  | Phase B 5 L1 flat-flux-identity tests | PASS |
| 2.4  | 21 Phase A boundary_face_flux tests DELETED | DONE |
| 3.1  | Spatial MMS spherical anisotropic ansatz | xfail (ERR-026 PARTIAL) |
| 3.2  | Spatial MMS cylindrical anisotropic ansatz | xfail (ERR-026 PARTIAL) |
| 3.3  | Angular convergence at fixed nx | PASS |
| 4.1  | k_∞ recovery 2G reflective sphere | PASS (rtol < 5e-4) |
| 4.2  | trajectory_resolvent cross-check | SKIP (Phase D placeholder) |

L1 analytical suite (`tests/sn/l1_analytical/`): 20 passed, 2 xfailed
(ERR-026 tripwires).

## What remains for Phase D

1. **Pole-face spatial closure refinement.** The WDD recurrence at
   the spherical pole with `psi_face_in(pole) = psi_cell[0]` plus
   `MorelMontryAngularSweep` angular closure produces inconsistent
   per-ordinate residuals. The fix is likely a Carlson-style
   coupled pole sweep where outward ordinates' pole-face
   initialisation is determined by inward ordinates' pole-face
   propagation (the symmetry condition at r=0 in continuous form).
2. **Default flip**: once 1. is in place, flip
   `SNMesh.pole_angular_closure` default
   `LegacyTauSymmetricInterpolation` → `MorelMontryAngularSweep`
   and the curvilinear `solve_sn_fixed_source` default
   `"source_iteration"` → `"krylov"`.
3. **Marker removal**: the 4 xfail-strict ERR-026 tripwires
   (sphere + cylinder × iso + aniso MMS) come off. ERR-026 status:
   PARTIAL CLOSURE → CLOSED.
4. **Snapshot regeneration**: regenerate the 6 curvilinear
   snapshots under the Phase D operator + cross-check each via
   `trajectory_resolvent` at rtol ≤ 5e-4.
5. **Gate 4.2 full implementation**: replace the SKIP placeholder
   with the actual SN-vs-trajectory_resolvent flux-shape agreement
   test on snapshots 1, 2, 4, 5, 6.
6. **error_catalog.md update**: ERR-026 status flip PARTIAL →
   CLOSED + Verification section pointing to Phase D's L1 +
   trajectory_resolvent cross-check evidence.

## LOC delta

Production code:

* `orpheus/sn/operator.py`: +210 / -200 (net ~+10; matvec rewrite
  + simplification of solution_to_angular_flux_spherical;
  EquationMap.unknowns_at_cell_for_mask helper added).
* `orpheus/sn/geometry.py`: +140 / -20 (iter_cells_by_direction
  helper; `boundary_face_flux` field removed).
* `orpheus/sn/solver.py`: -10 (tuple-unpack site simplifications).
* `orpheus/sn/spatial/__init__.py`: -10 (BoundaryFaceFlux exports).
* `orpheus/sn/spatial/boundary_face_flux.py`: -415 (DELETED).

Test code:

* `tests/sn/test_iter_cells_by_direction.py`: +260 NEW.
* `tests/sn/test_phase_c_gates.py`: +370 NEW.
* `tests/sn/test_phase_c_mms.py`: +135 NEW.
* `tests/sn/test_phase_c_crosscheck.py`: +130 NEW.
* `tests/sn/test_snstreamingoperator.py`: rewriting net ~-50.
* `tests/sn/spatial/test_boundary_face_flux.py`: -232 (DELETED).

Documentation:

* `docs/theory/discrete_ordinates.rst`: +130 / -90 (Phase A
  rewrite + Phase C subsection with 6 labels).
* `docs/theory/boundary_conditions.rst`: +15 (Phase C §16A.3 note).
* `.claude/skills/vv-principles/error_catalog.md`: +90 (Phase C
  ERR-026 narrative extension).
* `.claude/plans/sn_reshape.md`: +30 / -15 (Wave H Phase C row).

Snapshots: 6 curvilinear regenerated (binary files).

## Files touched (Phase C summary)

NEW production:
- (none)

NEW tests:
- `tests/sn/test_iter_cells_by_direction.py` (260 lines, 11 foundation tests)
- `tests/sn/test_phase_c_gates.py` (370 lines, ~30 tests across 5 gate sets)
- `tests/sn/test_phase_c_mms.py` (135 lines, 3 tests)
- `tests/sn/test_phase_c_crosscheck.py` (130 lines, 2 tests)

MODIFIED:
- `orpheus/sn/operator.py` (matvec rewrite, signature changes,
  EquationMap.unknowns_at_cell_for_mask)
- `orpheus/sn/geometry.py` (iter_cells_by_direction; field
  retirement)
- `orpheus/sn/solver.py` (tuple-unpack sites)
- `orpheus/sn/spatial/__init__.py` (BoundaryFaceFlux exports removed)
- `tests/sn/test_snstreamingoperator.py` (5 tests rewritten)
- `docs/theory/discrete_ordinates.rst` (Phase A subsection rewrite
  + Phase C subsection + 6 new :label: anchors)
- `docs/theory/boundary_conditions.rst` (§16A.3 Phase C note)
- `docs/verification/matrix.rst` (auto-regenerated)
- `.claude/skills/vv-principles/error_catalog.md` (ERR-026 Phase C
  extension)
- `.claude/plans/sn_reshape.md` (Wave H Phase C row flipped to DONE)

DELETED:
- `orpheus/sn/spatial/boundary_face_flux.py` (415 LOC)
- `tests/sn/spatial/test_boundary_face_flux.py` (232 LOC, 21 tests)

REGENERATED:
- 6 curvilinear regression snapshots

## Sub-agent dispatch (archivist)

The Phase C Sphinx subsection is **stub-extended**: it has the 6
`:label:` anchors, the architectural narrative, the verification-
gate cross-references, and the empirical-Gate-1.1-outcome rationale.
The rich-narrative expansion (full mathematical derivation of the
WDD sweep-frame matvec, the §16A.3 contract enforcement, the
pole-face Carlson starting-direction convention) is the archivist's
deliverable.

Per CLAUDE.md the dispatch is OWED, not optional. Drafting brief
here for the user to dispatch when ready:

> **DISPATCH_REQUEST to archivist** (for the user to invoke):
>
> Expand the Phase C subsection at
> `docs/theory/discrete_ordinates.rst` under the
> ``phase-c-sweep-frame-matvec`` label into a rich narrative
> following Cardinal Rule 3 maximum-effort standard. Source
> artefacts:
> - The 5 Phase C commits on `refactor/sn-operator-algebra`
>   (`eae6f05`, `16e33d7`, `3fd1302`, `3a95e3b`, `d445a8f`).
> - `orpheus/sn/operator.py` (the rewritten matvecs + their
>   docstrings — the algebra and the BC trace contract
>   enforcement).
> - `orpheus/sn/geometry.py` `iter_cells_by_direction` (the new
>   API and its foundation tests).
> - This closeout memo
>   (`.claude/agent-memory/method-implementer/issue_168_phase_c_closeout.md`).
> - The cross-domain-attacker memo at
>   `.claude/agent-memory/cross-domain-attacker/issue_168_phase_c_sweep_frame.md`
>   (the structural-frame analysis that informed the rewrite).
> - The Phase B closeout
>   (`.claude/agent-memory/method-implementer/issue_168_phase_b_closeout.md`)
>   for the spatial-closure-alignment narrative.
> - `docs/theory/boundary_conditions.rst` §16A.3 + §16A.5 for the
>   affine-BC-trace contract.
>
> Output: rich narrative covering: (a) the architectural shift
> from arithmetic-average spatial closure (Phase A) to WDD
> sweep-frame propagation (Phase C); (b) the §16A.3 contract
> enforcement and how the sweep-frame matvec honours it
> structurally; (c) the pole-face Carlson starting-direction
> convention and why `psi_face_in(pole) = psi_cell[0]` (not 0)
> preserves the per-ordinate flat-flux invariant; (d) the BoundaryFaceFlux
> Protocol retirement and the architectural reasoning ("two paths
> to the same operator → unify after the second instance"); (e)
> the empirical Gate 1.1 finding and what it reveals about the
> spherical-pole-vs-cylindrical-level structural asymmetry; (f)
> the Phase D scope (pole-face spatial closure refinement).

## Related work

- `.claude/plans/issue_168_phase_c.md` — the Phase C plan.
- `.claude/agent-memory/cross-domain-attacker/issue_168_phase_c_sweep_frame.md`
  — the structural-frame analysis that informed the rewrite.
- `.claude/agent-memory/method-implementer/issue_168_phase_a_closeout.md`
  — Phase A closeout.
- `.claude/agent-memory/method-implementer/issue_168_phase_b_closeout.md`
  — Phase B closeout (PoleAngularClosure Protocol + Hébert §3.9.4).
- ERR-026 in `.claude/skills/vv-principles/error_catalog.md` —
  PARTIAL CLOSURE status through Phase C.
- `scratch/literature/Hebert(2009)Chapter3.pdf` — Hébert §3.9.4
  primary source.
- `.claude/agent-memory/literature-researcher/sphere_sn_pole_closure_canonical.md`
  — Hébert canonical reference + Bailey citation correction
  (corrected across the codebase in Phase B).
