# The octant face-transmission map — an algebra of record for diamond, step and linear discontinuous

Opened 2026-09-23. Status: ruled in outline, not started. The ontology of the LD comparison is still being searched (`plan-authoring` §0), so this is a living plan, and the design discussion with the user comes before code.

## The instruction (the user, verbatim)

*"The face-transmission proposal you made is accepted, and we shall turn it into an excellent algebra of record improving and extending it (potentially showing how Step and LD do the job comparatively)."* (2026-09-22, ruling R5 of `.claude/plans/vv_suite_layout.md`.) Then, 2026-09-23: *"we will do that after you sort what you need to do with the performance test and prepare to compact context before the algebra of record development of the case."*

## What exists today `[M]` 2026-09-23

- `derivations/sn_dd_face_transmission.py`, the last file under the repository-root `derivations/`. For one Cartesian cell in d dimensions, the diamond closure `psi_out = 2 psi_cell - psi_in` per axis, with `psi_cell` eliminated through the cell balance, gives the d x d map from an octant's d inflow faces to its d outflow faces: `Sigma_DD = (2/D) 1 w^T - I`, with `w_a = 2 |mu_a| A_a` and `D = Sigma_t V + sum_b w_b`. It proves symbolically, per d:
  - Claim 1: `spec(Sigma_DD) = {1 - 2 Sigma_t V / D} u {-1}^(d-1)`;
  - Claim 2 (the control): step closure (`psi_out = psi_cell`) gives `Sigma_step = (1/D) 1 w^T`, with `spec = {(D - Sigma_t V)/D} u {0}^(d-1)`.
  
  It also carries a characteristic-polynomial cross-check (d <= 4) and a numeric spot-check.
- What rests on it:
  - the (d-1)-fold `-1` is the undamped face sawtooth (w-weighted zero-average faces, invisible to collision) that costs 1631 sweeps on the d = 3 all-reflective box, and it is why the default iteration budget targets rho = 0.986 (`orpheus/numerics/convergence.py:289`);
  - negative eigenvalues void Varga's regular-splitting hypothesis, so boundary Gauss-Seidel is not guaranteed faster than Jacobi (#341; the tree said "regular splitting" in 11 places until 2026-08-09);
  - step sends the mode to 0, so the `-1` belongs to diamond's `-psi_in` term, not to transport.
- Cited by: `docs/theory/methods/sn/cartesian_multid.rst` (the label `dd-face-transmission-spectrum` near line 4521, the section `sn-boundary-gs-not-regular` near 4443, and lines 75, 4795, 4801); `docs/theory/foundations/discretization.rst:1270`; `orpheus/numerics/convergence.py:289`; the docstring of `tests/gates/sn/solve/test_reflective_si_iteration_budget.py`.
- **Its defects:**
  - Nothing runs it: the label's V&V status is `documented` (`.. vv-status: dd-face-transmission-spectrum documented`), with no `verifies` edge.
  - Its checks are bare `assert`s, which `python -O` strips.
  - Both its proof and its numeric check build the matrix by hand. Nothing compares against the map the production sweep applies: a transcription gap its own docstring names.
- **The import gate:** `tests/gates/derivations/test_diagnostics_resolve_their_imports.py` scans `derivations/` and names this file as the known member of that root. Moving the file retires that root: drop `DERIVATIONS_ROOT` and its control member in the same commit.
- **The sibling that sets the pattern:** `orpheus/derivations/discrete/sn/sweep_acyclicity.py` ("the algebra of record for SN sweep acyclicity"), gated by `tests/gates/sn/sweep/test_sweep_acyclicity.py`.
- **Production closures:** `DiamondDifference` (`orpheus/transport/spatial/diamond.py:165`) and `LinearDiscontinuous` (`orpheus/transport/spatial/linear_discontinuous.py:245`). There is no step closure in production: step exists only as this derivation's control. LD in d >= 2 is the bilinear (UBLD) closure, with 2^d moments per cell (Branch 1: `orpheus/derivations/discrete/sn/ld_ubld.py`). It rides the DAG wavefront, not the scan-march (#38, #240 D5b).

## Rulings

- **R5** (the user, 2026-09-22): move it to `orpheus/derivations/discrete/sn/` as an algebra of record, improved and extended, possibly comparing step and LD with diamond. The proposal it accepted:
  1. the move, `assert`s becoming raises, and the module returning the matrices and spectra it proves;
  2. a test at `l0` with `verifies("dd-face-transmission-spectrum")`: the symbolic claims for d = 2 to 4, plus a leg that extracts the face-to-face map from the production sweep (unit inflow on each face of one cell, one octant) and asserts it equals `Sigma_DD`; a diamond-to-step mutation must redden it;
  3. repoint the page and `convergence.py` citations, and retire the repository-root `derivations/`.
- Workflow `[R]`: an algebra-of-record carve. Load the `algebra-of-record` skill. The main agent writes the derivation with the user steering; the test-architect designs the gates, with `rests_on` edges; qa and the elegance-enforcer review; the archivist writes the page.

## Open questions, for the discussion at resume

1. **LD's face map is not d x d.** A bilinear cell's face carries 2^(d-1) moments (its average and its transverse slopes), so the octant map acts on `d * 2^(d-1)` face moments, not d face averages. What is the right object to compare across closures: the full face-moment map; its restriction to face averages (which needs a closure of the transverse moments); or its spectrum on the zero-average subspace that carries diamond's sawtooth? `[HYPOTHESIS]` The comparative question the page needs answered is "does closure X carry an undamped mode at zero leakage?", which is a statement about the spectral radius and the unit-modulus eigenvalues of each map.
2. **Step is not a production closure.** Is it a control in the derivation only (as today), or a reason to build it? `[R]` Control only, unless the user rules otherwise.
3. **The production-extraction leg:**
   - which operator surface exposes the single-cell, single-octant inflow-to-outflow map (the sweep through `StreamingCollisionOperator.solve` on a one-cell all-reflective mesh, or the boundary operator `B` restricted to one octant);
   - whether LD's map can be extracted the same way;
   - the premise to measure first: that a one-cell mesh is admissible.
4. **What the extension claims.** For LD: its spectrum in closed form (SymPy), or only numerically for chosen d? The UBLD cell system is `2^d x 2^d` per cell, so a closed form may be heavy for d = 3.

## ⏸ COMPACTION POINT — 2026-09-23

Resume here. Nothing is started. First, read this file and the `algebra-of-record` skill, then `derivations/sn_dd_face_transmission.py`, then `orpheus/derivations/discrete/sn/sweep_acyclicity.py` with its test, as the pattern. Bring questions 1 to 4 to the user before any code: this is ontology search, and the LD comparison's object is not yet named.
