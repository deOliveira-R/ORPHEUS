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
├── _generate_snapshots.py   # tool: builds the .npz reference set
├── test_dd_regression.py    # test: re-runs each case + asserts match
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
solver output. The captured payload is `(keff, scalar_flux, case_name,
case_description, generator_commit)`. Commit the .npz files alongside
any code change that legitimately moves the numerical output.

### Run (CI / local)

```bash
pytest -m regression                    # all snapshot cases
pytest -m regression -k cyl             # cylindrical only
```

Tests skip cleanly if a snapshot file is absent — this lets the
infrastructure land before snapshots are generated, and protects
against accidentally treating a missing snapshot as a regression.

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

## Case roster (current)

| Snapshot name | Geometry | Groups | BC | Quadrature | Mesh |
|---|---|---|---|---|---|
| `slab_2g_homogeneous_dd_n20` | Slab | 2G | reflective/reflective | GL-1D N=8 | n=20 |
| `slab_2g_3reg_dd_n40` | Slab fuel/mod/fuel | 2G | reflective/reflective | GL-1D N=8 | n=40 |
| `sphere_2g_homogeneous_dd_n20` | Sphere | 2G | refl/refl (k_inf path) | GL-1D N=8 | n=20 |
| `sphere_2g_3reg_dd_n40` | Sphere fuel/mod/fuel | 2G | refl/refl | GL-1D N=8 | n=40 |
| `cyl_1g_homogeneous_LS4_dd_n20` | Cylinder | 1G | refl/refl | LS_4 (12 ord) | n=20 |
| `cyl_1g_homogeneous_product_dd_n20` | Cylinder | 1G | refl/refl | ProductQuadrature(2x4) | n=20 |
| `cyl_2g_3reg_LS4_dd_n40` | Cylinder fuel/mod/fuel | 2G | refl/refl | LS_4 | n=40 |

### Cases queued for follow-up

- `slab_p1_aniso_dd_n20` / `sphere_p1_aniso_dd_n20` — Pℓ scattering
  paths. Land after the SN reshape Pℓ-Galerkin Issue 13.
- `2d_1g_LS4_dd_15x15` — 2D wavefront sweep. Land alongside the 2D
  quadrature audit.
- `slab_fixed_source_dd_n20` — fixed-source path via
  `solve_sn_fixed_source`.
- `cyl_white_bc_dd_n20` — white BC path. Lands after SN reshape
  Issue 7 (`ResolvedBC` tensor decomposition).
