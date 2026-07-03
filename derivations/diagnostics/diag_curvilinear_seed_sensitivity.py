"""Diagnostic: is the curvilinear (L+C).solve a seed-DEPENDENT lagged inverse?

TAXONOMY-RULING investigation (#226). READ+DIAGNOSE ONLY.

Hypothesis under test: the 1-D curvilinear (sphere/cyl) SN sweep reads the
M-M half-angle starting seed psi_{1/2} from the PREVIOUS ITERATE, making
(L+C).solve(rhs, initial_guess=...) seed-DEPENDENT -- a true triangular
inverse ONLY at the SI fixed point, not a machine-precision one-shot inverse.

Probes:
  1. seed sensitivity: cold vs two random seeds (bulk rel-diff)
  2. residual r = A.apply(A.solve(b)) - b, normalized
  3. per-ordinate localization of the residual
  4. slab CONTROL (expect bit-exact seed-independence + machine residual)
  5. fixed-point iteration psi <- A.solve(b, ig=psi): does r -> machine?

Run: .venv/bin/python -O diag_curvilinear_seed_sensitivity.py
"""
from __future__ import annotations

import sys
import numpy as np

# Repo root on sys.path so tests.* and orpheus.* import.
REPO = "/Users/rodrigo/git/nuclear/ORPHEUS"
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from orpheus.geometry import (
    BC, CoordSystem, Mesh1D, Region, RegionMesh, StructuredGeometry,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import InvertibleOperator, StreamingOperator
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import placeholder_materials


# ── Geometry builders (mirror tests/sn/operators/test_removal_form_matvec_sweep.py)
def _slab(nx=6, n_ord=4, ng=2, bc="vacuum"):
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC(bc), BC(bc)),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _sphere(nx=6, ng=2, bc="vacuum"):
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC(bc),),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _cyl(nx=6, ng=2, bc="vacuum"):
    geom = StructuredGeometry(
        geometry="CYL",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC(bc),),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _removal_sigmas(sn, *, seed):
    rng = np.random.default_rng([seed, 240])
    sig_t = rng.uniform(1.0, 3.0, size=(sn.ng, *sn.spatial_shape))
    sig_s0 = rng.uniform(0.2, 0.8, size=(sn.ng, *sn.spatial_shape))
    sig_r = sig_t - sig_s0
    assert np.all(sig_r > 0.0)
    return sig_t, sig_r


def _random_state(sn, *, seed):
    rng = np.random.default_rng([seed, 7])
    state = TimedFullField.zeros(bulk=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn)
    state.bulk.values[...] = rng.standard_normal(state.bulk.values.shape)
    for face in state.boundary.layout.faces:
        fv = state.boundary.face_view(face)
        fv[...] = rng.standard_normal(fv.shape)
    return state


def _zero_boundary(field):
    for face in field.boundary.layout.faces:
        field.boundary.face_view(face)[...] = 0.0


def _build(sn, *, seed):
    sig_t, sig_r = _removal_sigmas(sn, seed=seed)
    op = InvertibleOperator(
        StreamingOperator(sn), MultiplicationOperator.from_mesh(sig_r, sn),
    )
    return op


# ═══════════════════════════════════════════════════════════════════
# Probe 1 + 2 + 3 : seed sensitivity, residual, per-ordinate localization
# ═══════════════════════════════════════════════════════════════════
def probe_seed_and_residual(name, sn, seed):
    op = _build(sn, seed=seed)

    # Purely volumetric vacuum source (zero the inflow trace).
    b = _random_state(sn, seed=seed + 2)
    _zero_boundary(b)

    # Two DIFFERENT random seeds + cold start.
    X1 = _random_state(sn, seed=seed + 100)
    X2 = _random_state(sn, seed=seed + 200)

    psi_none = op.solve(b)
    psi_X1 = op.solve(b, initial_guess=X1)
    psi_X2 = op.solve(b, initial_guess=X2)

    vN = psi_none.bulk.values
    v1 = psi_X1.bulk.values
    v2 = psi_X2.bulk.values

    def reldiff(a, c):
        denom = max(np.abs(a).max(), 1e-300)
        return float(np.abs(a - c).max() / denom)

    d_12 = reldiff(v1, v2)
    d_1none = reldiff(v1, vN)

    # Residual r = A.apply(A.solve(b)) - b on the COLD solve.
    q_back = op.apply(psi_none)
    r = q_back.bulk.values - b.bulk.values
    r_inf = float(np.abs(r).max())
    b_inf = float(np.abs(b.bulk.values).max())
    r_rel = r_inf / max(b_inf, 1e-300)

    # Per-ordinate residual localization (max over cells+groups per ordinate n).
    # r shape (N, ng, nx, ny). Reduce all but axis 0.
    r_per_ord = np.abs(r).reshape(r.shape[0], -1).max(axis=1)

    # Per-ordinate seed sensitivity (|v1 - v2| per ordinate).
    dd = np.abs(v1 - v2).reshape(v1.shape[0], -1).max(axis=1)

    mu = sn.quad.mu_x
    print(f"\n=== {name} (N={sn.quad.N}, nx={sn.nx}, ng={sn.ng}) ===")
    print(f"  seed sensitivity  max|psi_X1 - psi_X2| / max|psi_X1| = {d_12:.3e}")
    print(f"  seed vs cold      max|psi_X1 - psi_none|/ max|psi_X1| = {d_1none:.3e}")
    print(f"  residual  ||A.apply(A.solve(b)) - b||_inf / ||b||_inf = {r_rel:.3e}"
          f"   (abs {r_inf:.3e})")
    print(f"  per-ordinate residual (n: mu_n -> |r|_inf, |seedΔ|):")
    for n in range(sn.quad.N):
        print(f"     n={n:2d}  mu={mu[n]:+.4f}   |r|={r_per_ord[n]:.3e}   "
              f"|seedΔ|={dd[n]:.3e}")
    return dict(d_12=d_12, d_1none=d_1none, r_rel=r_rel)


# ═══════════════════════════════════════════════════════════════════
# Probe 5 : fixed-point iteration psi <- A.solve(b, ig=psi)
# ═══════════════════════════════════════════════════════════════════
def probe_fixed_point(name, sn, seed, n_iter=40):
    op = _build(sn, seed=seed)
    b = _random_state(sn, seed=seed + 2)
    _zero_boundary(b)

    print(f"\n=== FIXED POINT: {name} ===")
    psi = op.solve(b)  # cold start
    prev = None
    for k in range(n_iter):
        r = op.apply(psi).bulk.values - b.bulk.values
        r_rel = float(np.abs(r).max()) / max(np.abs(b.bulk.values).max(), 1e-300)
        incr = (
            float(np.abs(psi.bulk.values - prev).max())
            if prev is not None else float("nan")
        )
        if k < 8 or k % 5 == 0 or k == n_iter - 1:
            print(f"   iter {k:2d}  ||r||_rel={r_rel:.3e}   ||Δψ||_inf={incr:.3e}")
        prev = psi.bulk.values.copy()
        psi = op.solve(b, initial_guess=psi)
    return r_rel


if __name__ == "__main__":
    results = {}
    for name, sn, seed in [
        ("SLAB   (control)", _slab(ng=2, bc="vacuum"), 11),
        ("CYLINDER", _cyl(ng=2, bc="vacuum"), 22),
        ("SPHERE", _sphere(ng=2, bc="vacuum"), 33),
    ]:
        results[name] = probe_seed_and_residual(name, sn, seed)

    for name, sn, seed in [
        ("SLAB   (control)", _slab(ng=2, bc="vacuum"), 11),
        ("CYLINDER", _cyl(ng=2, bc="vacuum"), 22),
        ("SPHERE", _sphere(ng=2, bc="vacuum"), 33),
    ]:
        probe_fixed_point(name, sn, seed)

    print("\n\n================ SUMMARY ================")
    for name, r in results.items():
        print(f"{name:20s}  seedΔ(X1,X2)={r['d_12']:.2e}  "
              f"residual_rel={r['r_rel']:.2e}")
