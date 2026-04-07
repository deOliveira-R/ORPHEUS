# QA State — ORPHEUS

Mutable state for the QA agent. Updated after each session.
Last updated: 2026-04-06.


## Running Tests

```bash
.venv/bin/python -m pytest tests/ -v -k "not slow"   # ~500 tests, ~2 min
.venv/bin/python -m pytest tests/ -v                  # +slow (~10 min)
.venv/bin/python -m sphinx -b html docs docs/_build/html
```


## Test Classification

| File | ~Count | Level | What it covers |
|------|--------|-------|----------------|
| `test_cp_slab.py` | 9 | L1 | Semi-analytical eigenvalue (E₃) |
| `test_cp_cylinder.py` | 9 | L1 | Semi-analytical eigenvalue (Ki₄) |
| `test_cp_sphere.py` | 9 | L1 | Semi-analytical eigenvalue (exp) |
| `test_cp_properties.py` | 20 | L0 | Row sums, reciprocity, non-negativity |
| `test_cp_diagnostics.py` | 30 | L0 | Matrix structure, optical limits |
| `test_cp_verification.py` | 40 | L0–L2 | Upscatter, (n,2n), convergence rate |
| `test_sn_1d.py` | 20 | L1–L2 | Homogeneous exact, Richardson convergence |
| `test_sn_cylindrical.py` | 25 | L0–L2 | Curvilinear sweep, heterogeneous |
| `test_sn_spherical.py` | 25 | L0–L2 | Curvilinear sweep, heterogeneous |
| `test_sn_quadrature.py` | 28 | L0 | Weight sums, symmetry, α dome |
| `test_sn_solver_components.py` | 35 | L0–L1 | Sweep regression, normalization |
| `test_sn_properties.py` | 10 | L0 | Flux symmetry, particle balance |
| `test_sn_sweep_regression.py` | 5 | L0 | Bitwise sweep output |
| `test_moc.py` | 10 | L1 | Homogeneous exact, heterogeneous |
| `test_moc_verification.py` | 30 | L0–L2 | Weight formula, (n,2n), convergence |
| `test_moc_properties.py` | 20 | L0 | Particle balance, flux positivity |
| `test_moc_quadrature.py` | 15 | L0 | TY polar weights, azimuthal uniformity |
| `test_moc_ray_tracing.py` | 20 | L0 | Ray-circle intersection, reflective links |
| `test_monte_carlo.py` | 20 | L1 | z-score vs analytical, heterogeneous |
| `test_mc_properties.py` | 25 | L0 | Geometry protocol, 1G deterministic |
| `test_mc_gaps.py` | 30 | L0 | Gap region handling, tracking |
| `test_mc_convergence.py` | 5 | L2 | 1/√N convergence |
| `test_mc_cross_verification.py` | 5 | L4 | MC vs CP (benchmarking) |
| `test_homogeneous.py` | 3 | L1 | Matrix eigenvalue |
| `test_diffusion.py` | 5 | L1–L2 | Buckling eigenvalue, Richardson |
| `test_diffusion_properties.py` | 5 | L0 | Flux symmetry, positivity |
| `test_convergence.py` | 1 | L4 | SN vs CP (benchmarking) |
| `test_geometry.py` | 15 | L0 | Mesh volumes, coordinate systems |


## Known Gaps

- **MMS tests** (TS-20260403-002): no manufactured-solution tests for
  SN, MOC, or diffusion.  Biggest L1 gap.
- **L3 validation**: no ICSBEP/IRPhE comparison.  Aspirational for
  an educational code.
- **Richardson caching** (TS-20260403-001): heterogeneous references
  recomputed every run (~15 min).
- **L0 markers** (TS-20260405-002): no `@pytest.mark.l0` tagging yet.


## ORPHEUS-Specific Bug History

### Cylindrical DD (ERR-014/015)

Wrong α recursion + missing ΔA/w.  Passed: homogeneous exact (all
groups), particle balance, conservation, non-negativity.  Failed ONLY
on heterogeneous convergence (keff diverged: 1.15→0.90→0.52).

**QA rule:** curvilinear transport needs fixed-source flat-flux
diagnostic (uniform Q, Σ_t=1, 40+ cells).  Check φ ≈ Q/Σ_t AND
bounded flux range.

### MOC Weight Bug (ERR-019)

Weight formula `4π·ω_a·ω_p·t_s·sin(θ_p)` — the 4π and sin(θ_p)
cancel in keff for homogeneous.  Only heterogeneous multi-region
exposes the error.

### Cross-Section Convention

`Mixture.SigS[l][g_from, g_to]` — source uses transpose:
`Q = SigS^T @ phi`.  Swap gives wrong multi-group, correct 1-group.

### Quadrature Weight Sums

GL: `sum(w) = 2`.  Lebedev/LS/Product: `sum(w) = 4π`.  Wrong sum
cancels in keff — catch with `φ = Q/Σ_t` normalization test.

### The α Dome

Curvilinear α coefficients: non-negative dome (0→peak→0).  Negative
α → NaN.  Test: `assert np.all(alpha >= -1e-14)`.
