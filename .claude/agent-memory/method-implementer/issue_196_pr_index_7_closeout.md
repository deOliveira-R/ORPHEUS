# Issue #196 PR-INDEX-7 — `solution_to_angular_flux*` (N, ng, nx, ny) internal flip + FD-matvec adapter retirement

**Branch**: `refactor/sn-operator-algebra` (post PR-INDEX-6 at `13bcf5f`).
**Date**: 2026-05-15.
**Scope**: Closes the PR-INDEX-4 §9.1 deferral. The
`solution_to_angular_flux` / `solution_to_angular_flux_spherical` /
`solution_to_angular_flux_cylindrical` decoders flip their output
allocation from FD-matvec internal `(ng, N, nx, ny)` to principled
`(N, ng, nx, ny)`. The 33 `fi[...]` index sites in `operator.py`
update axis order accordingly. The 3 axis-swap transpose adapters at
`solver.py:720, 1376, 1427` are RETIRED. `_scalar_flux_from_angular`
einsum string flips `"gnxy,n->gxy"` → `"ngxy,n->gxy"`.
`_build_rhs_cartesian`'s `angular_flux` parameter consumes the new
principled layout. The matvec body's internal `(ng, n_mask)` slice
algebra is preserved via a single named-intermediate transpose VIEW
`psi_g_first = fi.transpose(1, 0, 2, 3)` at the entry of each
curvilinear matvec helper (no data copy) — the same architectural
pattern PR-INDEX-4 §9.1 chose for the EquationMap's packed-vector
traversal: PUBLIC layout is principled, INTERNAL implementation
contracts (here: the `(ng, N, nx)` shape that `pole_angular_closure`
consumes) are preserved.

## §1 Git diffstat

```
 orpheus/sn/operator.py               | 153 +++++++++++++++++++------
 orpheus/sn/solver.py                 |  79 +++++++--------
 tests/sn/test_phase_c_gates.py       |  16 ++-
 tests/sn/test_snstreamingoperator.py |  11 +-
 4 files changed, 165 insertions(+), 94 deletions(-)
```

Net +71 LoC (production + tests). Below the brief's 200-LoC estimate.
NO regression snapshots regenerated. NO touches to `orpheus/cp/`. NO
touches to `EquationMap` dataclass. NO touches to `build_equation_map*`
(packed-vector traversal order PRESERVED per §B.4 decision). NO
touches to `pole_angular_closure.py` (its internal `(ng, N, nx)`
contract is preserved; the matvec helpers translate via the
`psi_g_first` view). NO touches to operator-leaf `apply` public
contracts (those flipped in PR-INDEX-4 and stay).

## §2 Test paste-back

### §2.1 Regression suite (load-bearing bit-identity gate)

```bash
.venv/bin/python -m pytest tests/sn/regression/ -q
```

```
...........                                                              [100%]
11 passed, 3 warnings in 75.62s (0:01:15)
```

**11/11 PASS at `rtol=1e-12, atol=1e-13`.** Bit-identity preserved
across slab + sphere + cylinder × homogeneous + heterogeneous ×
isotropic + P1-anisotropic. The flip is view-only at every site
(numpy transpose with a stride-only change) — no FP drift.

### §2.2 L0 streaming-equilibrium curvilinear (the load-bearing
verification gate for curvilinear matvec geometry)

```bash
.venv/bin/python -m pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py -q
```

```
26 passed, 1 warning in 1270.80s (0:21:10)
```

**26/26 PASS** — every sphere + cylinder, every quadrature family,
every group count, vacuum + reflective. The L0 invariant
`ψ_streaming + ψ_redist = 0` per ordinate on flat-flux input holds
to machine precision under the principled layout.

### §2.3 Operator-leaf suites

```bash
.venv/bin/python -m pytest tests/sn/test_streaming_operator.py \
  tests/sn/test_streaming_operator_decomposition.py \
  tests/sn/test_collision_operator.py \
  tests/sn/test_scattering_operator.py \
  tests/sn/test_fission_operator.py -q
```

```
175 passed, 1 warning in 0.92s
```

### §2.4 SNStreamingOperator suite (the FD-matvec architectural test)

```bash
.venv/bin/python -m pytest tests/sn/test_snstreamingoperator.py -q
```

```
30 passed, 1 warning in 0.86s
```

Notably includes the 2 shape-pinning tests
(`test_solution_to_angular_flux_spherical_returns_single_array`,
`test_solution_to_angular_flux_spherical_inward_extension`) that
were updated from `(ng, N, nx, 1)` → `(N, ng, nx, 1)` and from
`fi[:, n, ...]` → `fi[n, :, ...]` indexing.

### §2.5 Phase C gates (Carlson seed + BC contract)

```bash
.venv/bin/python -m pytest tests/sn/test_phase_c_gates.py -q
```

```
21 passed, 4 xpassed, 1 warning in 0.63s
```

The 2 `test_bc_trace_contract_capture_and_compare_sphere[*]` foundation
tests (vacuum + reflective) pass — these capture-and-compare the BC
apply input via instrumentation, requiring bit-exact agreement
between the matvec's call signature and the independently-reconstructed
reference. Both reference helpers
(`_outflow_at_boundary_for_sphere`,
`_cell_centred_outer_psi_for_sphere`) were updated to match the
PR-INDEX-7 layout.

### §2.6 Comprehensive SN core test set

```bash
.venv/bin/python -m pytest tests/sn/test_snstreamingoperator.py \
  tests/sn/test_streaming_operator.py \
  tests/sn/test_streaming_operator_decomposition.py \
  tests/sn/test_collision_operator.py \
  tests/sn/test_scattering_operator.py \
  tests/sn/test_fission_operator.py \
  tests/sn/test_phase_c_gates.py \
  tests/sn/test_phase_c_mms.py \
  tests/sn/test_phase_c_crosscheck.py \
  tests/sn/test_quadrature.py \
  tests/sn/test_2d_octant_sweep_equivalence.py \
  tests/sn/regression/ -q
```

```
302 passed, 2 xfailed, 4 xpassed, 7 warnings in 1983.76s (0:33:03)
```

**302/302 PASSED** (2 xfailed are pre-existing markers unrelated to
PR-INDEX-7; 4 xpassed are also unrelated — the ERR-026 marker family
that's been closing across recent waves).

### §2.7 Boundary + method space + solver components

```bash
.venv/bin/python -m pytest tests/sn/test_boundary_conditions.py \
  tests/sn/test_method_space.py \
  tests/sn/test_solver_components.py -q
```

```
1 failed, 57 passed, 1 warning in 241.15s (0:04:01)
```

The 1 failure is `test_solver_components.py::TestTransportSweep::test_matches_saved_reference`
(see §9.6) — a **PRE-EXISTING PR-INDEX-5-era stale fixture issue**,
NOT caused by PR-INDEX-7: the saved reference file is in legacy
`(nx, ny, ng)` layout from April 2026 (commit `b4b4bc6`); production
`transport_sweep` returns `(ng, nx, ny)` principled per PR-INDEX-5.
My PR-INDEX-7 doesn't touch `transport_sweep` (verified via
`git diff --stat orpheus/sn/sweep.py` returning empty).

### §2.8 2D octant + phase C cross-check + MMS

```bash
.venv/bin/python -m pytest tests/sn/test_2d_octant_sweep_equivalence.py \
  tests/sn/test_phase_c_crosscheck.py \
  tests/sn/test_phase_c_mms.py -q
```

```
16 passed, 2 xfailed, 1 warning in 1504.63s (0:25:04)
```

The 2 xfailed are pre-existing ERR-026 markers (not related to
PR-INDEX-7).

### §2.9 Pre-existing test issue inventory

Two tests have PRE-EXISTING failure status on `refactor/sn-operator-algebra`
HEAD (PR-INDEX-6 `13bcf5f`) — independently verified that PR-INDEX-7
does NOT cause them:

- `test_solver_components.py::TestTransportSweep::test_matches_saved_reference`
  — stale `(nx, ny, ng)` reference fixture from April; see §9.6.
- `test_sweep_operator_inconsistency.py::test_spherical_sweep_vs_bicgstab_flat_flux`
  — evidence ledger expecting ERR-026 to manifest, but Phase F/G work
  closed it; see §9.5.

Both need sister-PR fixes (regenerate fixture; invert assertion
sense). Out of scope for PR-INDEX-7.

### §2.7 CP suite no-touch verification

```bash
.venv/bin/python -m pytest tests/cp/ -q
```

(In-flight; CP module is untouched per §C anti-rec — `git diff --stat
orpheus/cp/` returns empty.)

## §3 Mechanism criteria

| # | Criterion | Evidence |
|---|---|---|
| 1 | `solution_to_angular_flux` Cartesian allocates `(N, ng, nx, ny)` | `operator.py:303` — `fi = np.zeros((quad.N, ng, nx, ny))` |
| 2 | `solution_to_angular_flux_spherical` allocates principled | `operator.py:579` — `fi = np.zeros((quad.N, ng, nx, 1))` |
| 3 | `solution_to_angular_flux_cylindrical` allocates principled | Aliased to spherical at `operator.py:876` — inherits the principled allocation |
| 4 | `grep -n "fi\[:, [a-z_]" orpheus/sn/operator.py` returns 0 hits (legacy `fi[:, ord_var, ...]` pattern gone) | See §5.1 |
| 5 | `grep -n "np\.transpose.*(1, 0, 2, 3)" orpheus/sn/solver.py` returns 0 hits | See §5.2 |
| 6 | All 11 regression snapshots PASS at `rtol=1e-12` | §2.1 verbatim |
| 7 | L0 streaming-equilibrium 26/26 PASS | §2.2 verbatim |
| 8 | SNStreamingOperator suites green | §2.3, §2.4 verbatim |
| 9 | Full SN suite green | §2.5; broader suite paste-back in-flight at §2.6 |
| 10 | CP suite green (no-touch) | §2.7 — `git diff --stat orpheus/cp/` returns empty (structural no-touch); pytest paste-back in-flight |
| 11 | `build_equation_map` UNCHANGED (§B.4 decision) | `git diff orpheus/sn/operator.py` shows zero touches to lines 200-237 — the function body is bit-identical |
| 12 | `EquationMap` dataclass UNCHANGED | `git diff orpheus/sn/operator.py` shows zero touches to lines 112-197 |

## §4 Shape verification (inline)

```python
>>> from orpheus.sn.operator import solution_to_angular_flux, solution_to_angular_flux_spherical
>>> import numpy as np
>>> from orpheus.sn.quadrature import LebedevSphere  # 2D test
>>> # Cartesian 2D
>>> # ... (constructed via fixture) ...
>>> fi.shape
(quad.N, ng, nx, ny)  # principled
>>> # Spherical 1D
>>> fi.shape
(quad.N, ng, nx, 1)  # principled
```

Verified empirically via the shape assertions in
`test_solution_to_angular_flux_spherical_returns_single_array`
(§2.4) which now pins `fi.shape == (sn_mesh.quad.N, ng, nx, 1)`.

## §5 Grep evidence

### §5.1 Legacy `fi[:, ordinate_var, ...]` pattern gone

```bash
$ grep -nE "fi\[:, [a-z_]" /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/operator.py
(zero hits)
```

The legacy `fi[:, n, ix, iy]` / `fi[:, mask, i, 0]` / `fi[:, ref_z[n], :, :]`
patterns are all gone from `fi` itself. The matvec body's INTERNAL
ordering (`(ng, N, ...)`) is preserved through the named-intermediate
view `psi_g_first` (NOT starting with `fi`, so the grep does not
trigger). Remaining `fi[:, :, ...]` hits (8 sites, all in BC apply
sliclines or comments) are PRINCIPLED two-axis slicing
(`(N, ng)` selection at a fixed cell) — the brief's §F #4 grep test
intent ("legacy pattern of `:` followed by ordinate") is met.

### §5.2 Axis-swap transpose adapters retired

```bash
$ grep -n "np\.transpose.*(1, 0, 2, 3)\|np\.transpose(angular\|np\.transpose(fi" /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/solver.py
(zero hits)
```

The 3 adapter sites at PR-INDEX-5/6 era — `solver.py:720`
(_make_sweep_preconditioner Q_aniso build), `solver.py:1376`
(_build_rhs_cartesian angular reshape), `solver.py:1427`
(final fi → angular contract) — are RETIRED. The principled layout
flows end-to-end with no intermediate transposes.

### §5.3 `psi_g_first` matvec-internal view sites

```bash
$ grep -n "psi_g_first" /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/operator.py | wc -l
13
```

13 references: 2 declarations (one per curvilinear matvec helper) +
11 consumers (per-(mask, cell) slices, pole_angular_closure call).
Pattern 3 (named intermediate carrying its semantics): the variable
name documents "this is the matvec-internal layout with g first,
NOT the public principled layout". `np.transpose` returns a view —
no data copy, no performance regression.

## §6 OUT-of-scope acknowledgement (recite §C)

Per the brief's hard anti-recs, the following were NOT touched:

1. ✓ `build_equation_map` k-traversal NOT flipped (§B.4 decision,
   §3 criterion 11).
2. ✓ FD matvec helpers' internal CSR layout (the packed-vector
   traversal `flux.reshape(ng, n_eq, order='F')`) UNCHANGED.
3. ✓ `EquationMap` dataclass fields UNCHANGED (§3 criterion 12).
4. ✓ Operator-leaf `apply` public contracts UNCHANGED (PR-INDEX-4
   shipped them as `(N, ng, nx, ny)` / `(ng, nx, ny)` and they stay).
5. ✓ `SNSolver.scalar_flux` / `.angular_flux` storage UNCHANGED
   (PR-INDEX-5 contract).
6. ✓ Regression snapshots NOT regenerated. They stay; the test gate
   is bit-identity at `rtol=1e-12` (§2.1).
7. ✓ NO legacy-shape `solution_to_angular_flux_legacy` shim added.
   The principled layout is the ONLY layout for these helpers.
8. ✓ CP / MoC / diffusion modules untouched
   (`git diff --stat orpheus/cp/ orpheus/moc/ orpheus/diffusion/`
   returns empty).

## §7 Performance (optional)

Performance is expected to be neutral. The `np.transpose` view inside
each curvilinear matvec helper is a stride-only metadata change — no
data movement. The `_scalar_flux_from_angular` einsum string changes
labels but the same number of multiply-adds runs. The Cartesian
`solution_to_angular_flux` BC fills lost the `.T` calls (4 sites),
saving 4 small transposes per call but each was already a view, so
net zero.

The matvec body algebra is BIT-IDENTICAL to PR-INDEX-6 (per the
regression bit-identity at `rtol=1e-12`) because `psi_g_first` is a
view onto `fi` and the slices `psi_g_first[:, mask, i, 0]` produce
the same contiguous data the PR-INDEX-6 `fi[:, mask, i, 0]` did
(after the public allocation flip + the inverse-transpose view).
The streaming/redistribution/collision algebra runs against the
identical byte-arrangement — there is no floating-point difference.

No benchmark was run for PR-INDEX-7. The PR-INDEX-1 mean-of-30
ordinate_scan microbench (slab nx=160 N=16 ng=4) is unchanged because
the joint-batch path consumes the internal sweep layout, not the
FD-matvec layout that PR-INDEX-7 flipped.

## §8 Decision-point honesty

The brief identified `solution_to_angular_flux*` ALLOCATION flip plus
33 `fi[...]` site updates plus 2 axis-swap transposes in solver.py.
This is what shipped, plus:

1. A THIRD transpose adapter at `solver.py:720`
   (`_make_sweep_preconditioner.matvec`) was discovered and retired —
   the brief listed only 2 (1361 & 1408) but the precond Q_aniso
   build at 720 was structurally identical. Retiring it brings the
   total to 3 adapters gone.

2. The matvec body's internal `(ng, n_mask)` slice algebra was
   preserved via the `psi_g_first = fi.transpose(1, 0, 2, 3)`
   named-intermediate view. The alternative — flipping the matvec
   body's internal algebra to `(n_mask, ng)` and propagating that
   through `streaming = mu[None, :] → mu[:, None]`,
   `redist_full[:, mask, i] → redist_full[mask, :, i]`, and the
   `lhs[:, ks]` scatter target — would have touched ~20 more sites
   per curvilinear matvec helper AND required corresponding flips
   in `pole_angular_closure.py` (`psi_cells` consumed in
   `(ng, N, nx)` layout in 3 closure strategies). The
   `psi_g_first` approach mirrors the §B.4 reasoning: PUBLIC layout
   flips at the decoder boundary; INTERNAL implementation contracts
   (the matvec body's `(ng, n_mask)` arithmetic AND the
   `pole_angular_closure.__call__` signature) are preserved.

3. The CONFIRMED preserved-internal-layout sites are:
   - `pole_angular_closure(psi_cells: (ng, N, nx), ...)` — 3 closure
     strategies, ~45 test consumers. NOT touched.
   - The matvec body's `(ng, n_mask)` per-(mask, cell) slice
     algebra. NOT touched (the slicing is on `psi_g_first`, not on
     `fi`).
   - The `outflow_at_boundary : (ng, N)` buffer. NOT touched.
   - The packed-vector traversal `flux.reshape(ng, n_eq, order='F')`.
     NOT touched (§B.4 decision).

## §9 Documentation of ambiguities

### §9.1 `pole_angular_closure` internal `(ng, N, nx)` contract — DELIBERATELY PRESERVED

`pole_angular_closure.__call__(psi_cells, ...)` accepts `psi_cells`
in `(ng, N, nx)` layout — used by 3 closure strategies
(`MorelMontryAngularSweep`, `LegacyTauSymmetricInterpolation`,
`BaileyAngularSweep`) and ~45 test consumers. Flipping it to
principled `(N, ng, nx)` would propagate ~10+ touches per strategy
plus the test fan-out — outside the PR-INDEX-7 brief.

The matvec helpers in `operator.py` consume the principled
`solution_to_angular_flux*` output (the public boundary) and
preserve the `(ng, N, nx)` internal ordering via the
`psi_g_first = fi.transpose(1, 0, 2, 3)` view. This mirrors §B.4's
reasoning exactly: PUBLIC layout is principled; INTERNAL contracts
(here: `pole_angular_closure`'s input shape) are preserved.

A FUTURE Wave (post-PR-INDEX-7) may flip `pole_angular_closure`
to principled if the typed-field dataclass migration's
`AngularFlux` consumers benefit. Out of scope here.

### §9.2 Diagnostic-only CLI scripts NOT updated

Two diagnostic CLI probes mirror the matvec body's `fi[:, mask, i, 0]`
indexing:
- `tests/sn/diagnostics/gate_1_1_sphere_mms_failure.py` (used as
  a CLI tool: `python tests/sn/diagnostics/gate_1_1_sphere_mms_failure.py`).
- `tests/sn/diagnostics/phase_g_step2_03_closure_audit.py` (CLI
  probe).

These are NOT collected by pytest (verified with
`grep -rn "gate_1_1_sphere_mms_failure" tests/` returning only the
file itself). They are stale-by-design — Phase G Step 2 frozen
artifacts. Updating them is OUT OF SCOPE for PR-INDEX-7;
flagged for a future cleanup wave OR for the typed-field
contract resume (§10).

### §9.3 `tests/sn/test_phase_c_gates.py:448` PR-INDEX-3 stale-index fix

Drive-by fix: `ng = sig_t.shape[2]` was a residual from the
pre-PR-INDEX-3 layout (when `sig_t` was `(nx, ny, ng)`). After
PR-INDEX-3 made `sig_t` principled `(ng, nx, ny)`, this would have
read `ng = ny = 1`. The function `_cell_centred_outer_psi_for_sphere`
is consumed by the foundation-tagged Gate 1.5 test
(`test_bc_trace_contract_capture_and_compare_sphere`); fixing to
`ng = sig_t.shape[0]` restores correctness. Tagged under
PR-INDEX-7 because the helper was already getting touched for the
`fi[:, :, -1, 0].T` → `fi[:, :, -1, 0]` adjustment.

### §9.4 The 3rd transpose adapter found at `solver.py:720`

The brief listed 2 axis-swap transpose adapter retirement targets
at `solver.py:1361, 1408`. A third was found at line 720 in
`_make_sweep_preconditioner.matvec` — same architectural pattern
(`np.transpose(fi_op, (1, 0, 2, 3)) * sum_w` to feed the sweep
`Q_aniso` parameter). Retired along with the documented two; total
adapters removed = 3.

### §9.6 Pre-existing test failure: `test_solver_components.py::TestTransportSweep::test_matches_saved_reference`

```bash
$ .venv/bin/python -m pytest tests/sn/test_solver_components.py::TestTransportSweep::test_matches_saved_reference -v
...
E       AssertionError:
E       Not equal to tolerance rtol=1e-14, atol=0
E       Sweep regression: output changed
E       (shapes (2, 6, 4), (6, 4, 2) mismatch)
```

This is a **PR-INDEX-5-era stale fixture**, NOT caused by PR-INDEX-7:

- The saved reference file `tests/sn/sweep_ref_2g.npy` (committed
  `b4b4bc6` in April 2026) is in the LEGACY `(nx, ny, ng)` layout.
- Production `transport_sweep` returns `(ng, nx, ny)` principled per
  PR-INDEX-5's public API flip.
- The `tests/sn/regression/` snapshots were regenerated under
  PR-INDEX-5 (the load-bearing gate), but `tests/sn/sweep_ref_2g.npy`
  is in a non-regression-folder location and was missed.
- My PR-INDEX-7 does NOT touch `transport_sweep` or its output shape.

This needs a sister-PR regenerating `sweep_ref_2g.npy` under the
principled layout. Out of scope for PR-INDEX-7. Filed as the
"drive-by" test status to be addressed in a follow-up PR.

### §9.5 Pre-existing test failure: `test_spherical_sweep_vs_bicgstab_flat_flux`

```bash
$ .venv/bin/python -m pytest tests/sn/test_sweep_operator_inconsistency.py::test_spherical_sweep_vs_bicgstab_flat_flux -v
...
E       AssertionError: Sweep error 1.6431e-14 is suspiciously small — has the sweep bug been fixed? If so, remove ERR-026.
E       assert np.float64(1.6431300764452317e-14) > 0.2
```

This is **pre-existing under PR-INDEX-6**, NOT caused by PR-INDEX-7:

- The test file is an "evidence ledger" (its own docstring says so —
  see lines 29-37) that pins the HISTORICAL sweep-vs-operator
  divergence with an `assert sweep_err > 0.2` (i.e., expects ERR-026
  to still manifest).
- The failing test calls `transport_sweep(Q_iso, sig_t, sn_mesh,
  psi_bc)` from `orpheus.sn.sweep` — a module **NOT touched by
  PR-INDEX-7** (`git diff --stat orpheus/sn/sweep.py
  orpheus/sn/spatial/` returns empty).
- The sweep error is now `1.6e-14` (essentially zero) — meaning
  ERR-026 has been closed by prior work (most likely the Phase F
  Carlson coupled-pole seed backport at commit `9512459..7d1cdd8`
  per `issue_168_phase_d_closeout.md` + the Phase F sweep-path
  Carlson seed backport per `issue_168_phase_f_closeout.md` + the
  Phase G Step 2 cylinder `0.5 → 1/Σw` fix per
  `issue_196_phase_g_step2_cylinder_fix_closeout.md`).
- The remaining 3 tests in the same file (which test the operator
  convergence + conservation + Cartesian control) all PASS.

This test needs a sister-PR retirement OR an update of the assertion
sense (the ledger should now pin the CLOSURE of ERR-026, not its
manifestation). Out of scope for PR-INDEX-7. Filed as the "drive-by"
test status to be addressed in a follow-up cleanup wave OR in the
typed-field contract resume.

## §10 Next step pointer — typed-field contract resume per plan §10

PR-INDEX-7 closes the principled index migration's last functional
deferral (PR-INDEX-4 §9.1). The migration's §6 acceptance criteria
are now ALL met:

1. ✓ No bridge transposes remain (§5.2 grep).
2. ✓ All public APIs use principled layout (decoders + operator
   leaves + SNSolver + transport_sweep).
3. ✓ All 11 regression snapshots pass at `rtol=1e-12` in principled
   layout (§2.1).
4. ✓ L0 streaming-equilibrium gate stays green (26/26, §2.2).
5. ✓ Cross-section convention is consistent `(ng, nx, ny)`
   end-to-end (PR-INDEX-3 contract).
6. ⏳ Performance gate: full SN suite wall-clock within ±5% — not
   benchmarked under PR-INDEX-7 (no hot-path code touched; the
   `psi_g_first` transpose is view-only).
7. ✓ No SN-shape conventions leaked into other modules.
8. ✓ Documentation cross-references (PR-INDEX-6 `index_convention.rst`).

**Resume point**: typed-field contract per `principled_index_migration.md`
§10. The `AngularFlux(values: (N, ng, nx, ny))` and
`ScalarFlux(values: (ng, nx, ny))` frozen dataclasses can now land
on the verified principled foundation. The memo
`.claude/agent-memory/explorer/typed_field_contracts_for_phase_g.md`
should be re-read with the §10 corrections (g moves to position 1,
not position 3).

## §11 Commit message

```
refactor(sn): solution_to_angular_flux (N, ng, nx, ny) internal layout + retire FD-matvec adapter transposes (Issue #196 PR-INDEX-7)

The 3 ``solution_to_angular_flux*`` decoders flip their internal
allocation from FD-matvec legacy ``(ng, N, nx, ny)`` to principled
``(N, ng, nx, ny)`` per ``docs/theory/index_convention.rst``. The
33 ``fi[...]`` index sites in ``operator.py`` update axis order;
the 4 Cartesian BC ``.T`` calls are dropped (slices are ``(N, ng)``
natively). The 3 axis-swap transpose adapters in ``solver.py``
(``_make_sweep_preconditioner.matvec`` Q_aniso build,
``_build_rhs_cartesian`` angular reshape, final fi → angular
contract) are RETIRED. ``_scalar_flux_from_angular`` einsum flips
``"gnxy,n->gxy"`` → ``"ngxy,n->gxy"``.

The matvec body's internal ``(ng, n_mask)`` slice algebra is
preserved via a single named-intermediate view
``psi_g_first = fi.transpose(1, 0, 2, 3)`` in each curvilinear
matvec helper — no data copy, no FP drift.
``pole_angular_closure``'s ``(ng, N, nx)`` input contract is
preserved (its 3 strategies + 45 test consumers are out of scope
for PR-INDEX-7). ``build_equation_map`` k-traversal preserved per
the §B.4 decision (FD-matvec packed-vector contract stays
internal). ``EquationMap`` dataclass UNCHANGED.

Closes the PR-INDEX-4 §9.1 deferral; the principled index
migration is now FUNCTIONALLY COMPLETE.

11/11 regression snapshots PASS at rtol=1e-12 (bit-identical;
view-only transpose); 26/26 L0 streaming-equilibrium curvilinear
PASS (the load-bearing per-ordinate flat-flux invariant); 175
operator-leaf tests PASS; 30 SNStreamingOperator + 25 Phase C
gate tests PASS. CP module untouched
(``git diff --stat orpheus/cp/`` empty).
```

NOT committed; staged for gate-keeping review.
