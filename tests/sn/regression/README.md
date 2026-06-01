# SN frozen-reference regression suite

This directory holds the **Issue 16 redesign** of the SN regression
gate — a minimum frozen-reference set proving operator-algebra
equivalence under DD across geometry / group-count / quadrature / BC
combinations. It is **not** a verification reference; it is a
**numerical-drift detector** for refactors.

End-to-end correctness verification runs through the L0/L1
analytical-reference suite (see `tests/sn/l1_analytical/` and the
F_N / Case-method / `trajectory_resolvent` cross-verifications under
`tests/derivations/`).

## Layout

```
regression/
├── _generate_snapshots.py   # tool: builds the .npz reference set + run config (SoT)
├── _regression_assert.py    # principled tolerance + DriftWarning tripwire (SoT)
├── test_dd_regression.py    # test: re-runs each case + asserts principled match
├── snapshots/               # the .npz files (committed to repo)
└── README.md                # this file
```

## Workflow

### Generate (once, before SN reshape branches off)

```bash
python -m tests.sn.regression._generate_snapshots
# or for a single case:
python -m tests.sn.regression._generate_snapshots \
    --case slab_2g_homogeneous_dd_n20
# or list all cases:
python -m tests.sn.regression._generate_snapshots --list
```

Each invocation overwrites `snapshots/<name>.npz` with the current
solver output. The captured payload is `(case_kind, scalar_flux,
case_name, case_description, generator_commit[, keff])` — `keff` is
present only for `case_kind="eigen"`. Commit the .npz files alongside
any code change that legitimately moves the numerical output.

**Before trusting any regenerated snapshot**, corroborate the value
against a structurally-independent reference (NOT the solver itself —
that is circular). The roster table below names the reference used for
each case. If a regenerated value does NOT match its independent
reference, STOP: the solver is wrong for that case, not merely drifted.

### Run (CI / local)

```bash
pytest -m regression                    # all snapshot cases
pytest -m regression -k cyl             # cylindrical only
```

Tests skip cleanly if a snapshot file is absent — this lets the
infrastructure land before snapshots are generated, and protects
against accidentally treating a missing snapshot as a regression.

## Principled tolerance — no magic floors

The gate carries **no hand-picked tolerance numbers**. Everything lives
in `_regression_assert.py::assert_regression`:

* **Correctness gate (hard fail)** = `SAFETY × conv_tol`, where
  `conv_tol` is the solver's OWN convergence stopping criterion for the
  pinned quantity, read off the case run config (`keff_tol` for k_eff,
  `flux_tol` for the eigen flux, `inner_tol` for the fixed-source flux).
  `SAFETY = 10` is the iteration-map amplification headroom (`ρ/(1−ρ)`
  for `ρ ≲ 0.9`), not a fudge. Asserting tighter than the solver
  converged is unphysical. Single-sweep ("direct") results would instead
  use `assert_array_almost_equal_nulp(nulp=reduction_depth)` — every
  case in the current roster is an iterated solve, so the iterative
  branch is what runs.
* **Drift tripwire (informational)** = a `DriftWarning` emitted when a
  value clears the correctness gate but moved beyond bit-identity (ULP
  distance > 0). pytest summarises it at run end; `logging.debug` carries
  the per-element ULP forensic breakdown (`--log-cli-level=DEBUG`).

The convergence config in `_generate_snapshots.py` (`EIGEN_DEFAULTS`,
`FIXED_SOURCE_DEFAULTS`, per-case overrides) is the **single source of
truth** shared by the generator and the test — they cannot drift on what
"converged" means.

### Strict bit-identity gate for pure-refactor PRs

A refactor that claims **zero numerical change** should run the suite
with the tripwire escalated to a hard failure:

```bash
pytest tests/sn/regression/ \
  -W "error::tests.sn.regression._regression_assert.DriftWarning"
```

Any sub-tolerance bit drift then fails the build, turning the
informational tripwire into a strict bit-identity contract.

### `-O` safety

The canonical ORPHEUS invocation is `python -O -m pytest`, which strips
bare `assert` in non-rewritten helper modules. `assert_regression` uses
only `numpy.testing.*` (which raises unconditionally) and explicit
`raise`, never a bare `assert`, so the correctness gate stays live under
`-O`.

## When a snapshot test fails

Two paths, never just "regenerate the snapshot":

1. **Unintended drift** — the refactor shouldn't have changed the
   numerical output. Investigate root cause; fix; re-run; `assert
   match` again.
2. **Intended change** — the refactor legitimately updates the output
   (e.g., a positivity-preserving cell-update replaces DD as default).
   Audit why the new output is correct (with V&V evidence — usually a
   new L1 analytical-reference test confirming the new value is right);
   re-run the generator; commit BOTH the new snapshot AND the audit
   narrative in the same commit. The commit body must explain the
   shift and reference the V&V evidence.

The frozen reference is a **gate against forgetting V&V**. It catches
"oops, the output silently changed" — it does not, on its own,
distinguish "good drift" from "bad drift". The audit is non-optional.

## Case roster (current — 11 snapshots)

| Snapshot name | Kind | Geometry | Groups | Quadrature | Mesh | Pℓ | Independent corroboration |
|---|---|---|---|---|---|---|---|
| `slab_2g_homogeneous_dd_n20` | eigen | Slab | 2G | GL-1D N=8 | n=20 | P0 | closed-form k_inf = 1.875 |
| `slab_2g_3reg_dd_n40` | eigen | Slab fuel/mod/fuel | 2G | GL-1D N=8 | n=40 | P0 | reflective balance prod/abs = keff |
| `sphere_2g_homogeneous_dd_n20` | eigen | Sphere | 2G | GL-1D N=8 | n=20 | P0 | closed-form k_inf = 1.875 |
| `sphere_2g_3reg_dd_n40` | eigen | Sphere fuel/mod/fuel | 2G | GL-1D N=8 | n=40 | P0 | reflective balance prod/abs = keff |
| `cyl_1g_homogeneous_LS4_dd_n20` | eigen | Cylinder | 1G | LS_4 (12 ord) | n=20 | P0 | closed-form k_inf = 1.5 |
| `cyl_1g_homogeneous_product_dd_n20` | eigen | Cylinder | 1G | Product(2x4) | n=20 | P0 | closed-form k_inf = 1.5 |
| `cyl_2g_3reg_LS4_dd_n40` | eigen | Cylinder fuel/mod/fuel | 2G | LS_4 | n=40 | P0 | reflective balance prod/abs = keff |
| `slab_2g_p1_aniso_dd_n20` | **fixed_source** | Slab (B mixture) | 2G | GL-1D N=8 | n=20 | **P1** | flat inf-medium flux (diag Σ_t − Σ_s0ᵀ)⁻¹Q |
| `sphere_2g_p1_aniso_dd_n20` | **fixed_source** | Sphere (B mixture) | 2G | GL-1D N=8 | n=20 | **P1** | flat inf-medium flux |
| `2d_1g_LS4_dd_15x15` | eigen (krylov) | 2D Cartesian | 1G | LS_4 | 15×15 | P0 | closed-form k_inf = 1.5 |
| `slab_fixed_source_dd_n20` | **fixed_source** | Slab vacuum | 1G | GL-1D N=8 | n=20 | P0 | global source = abs + leakage |

The two P1 cases pin the Pℓ Galerkin assembly path that the SN reshape
Issue 13 (``ScatteringOperator``) refactors. They are **fixed-source,
not eigenvalue**: mixture B is a non-multiplying moderator
(``Σ_f = νΣ_f = 0``), so an eigenvalue formulation is malformed
(``k = production/absorption = 0/abs`` divides by zero — the legacy
``eigen`` snapshots for these two cases were all-NaN and gated nothing).
A uniform-source reflective fixed-source solve exercises the same
Galerkin path and produces the well-posed flat infinite-medium flux. The
2D case pins the wavefront sweep diagonal scheduling that Issue 12
(unified sweep) preserves; it uses ``inner_solver="krylov"`` because the
2-D Cartesian source-iteration path is deferred to Phase A (R-1 Step E).
The slab fixed-source case pins the ``solve_sn_fixed_source`` entry
point that the operator-algebra reshape rewires.

### Schema

The eigenvalue cases store ``(case_kind="eigen", keff, scalar_flux)``;
the fixed-source cases store ``(case_kind="fixed_source", scalar_flux)``
only — there is no ``k_eff`` for the pure transport operator. The
regression test dispatches on ``case_kind`` so both schemas coexist
in the same ``snapshots/`` directory.

### Cases queued for follow-up

- `cyl_white_bc_dd_n20` — white BC path. Lands after SN reshape
  Issue 7 (`BoundaryOperator` tensor decomposition).
