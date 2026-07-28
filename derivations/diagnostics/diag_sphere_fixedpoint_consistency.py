"""Diagnostic: is the sphere (L+C).solve EXACT when seeded AT the fixed point?

Decides whether the sphere seed-dependence is a *lagged inverse of the CORRECT
operator* (fixed point correct; only the seed is lagged) vs a *wrong operator*.

Probe 6 (fixed-point consistency): pick arbitrary psi0; b = A.apply(psi0);
   psi_back = A.solve(b, initial_guess=psi0). If psi_back == psi0 to machine
   precision, the sweep IS the exact inverse WHEN SEEDED AT THE FIXED POINT ->
   the operator's fixed point is CORRECT and the seed-lag is the ONLY defect.
   (This is the solve o apply = identity test WITH the correct seed -- the
   test-suite round-trip uses the WRONG (cold) seed, hence its divergence.)

Probe 7 (production sigma + physical source): repeat probe 1/2 with
   sigma_C == sigma_t (production) and a smooth POSITIVE source, to rule out
   the divergence being a removal-form / random-source artifact.

Probe 8 (does real sphere SN converge?): a genuine fixed-source SI solve, to
   confirm production sphere is NOT broken -- the bare-operator seed-iteration
   divergence is specific to iterating solve-o-apply on the raw operator.

Run: .venv/bin/python -O diag_sphere_fixedpoint_consistency.py
"""
from __future__ import annotations

import sys
import numpy as np

REPO = "/Users/rodrigo/git/nuclear/ORPHEUS"
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from orpheus.geometry import (
    BC, Mesh1D, Region, RegionMesh, StructuredGeometry,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import StreamingCollisionOperator, StreamingOperator
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import placeholder_materials


def _sphere(nx=6, ng=2, bc="vacuum"):
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC(bc),),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _random_state(sn, *, seed):
    rng = np.random.default_rng([seed, 7])
    state = TimedFullField.zeros(bulk=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn)
    state.bulk.values[...] = rng.standard_normal(state.bulk.values.shape)
    for face in state.boundary.layout.faces:
        state.boundary.face_view(face)[...] = rng.standard_normal(
            state.boundary.face_view(face).shape
        )
    return state


def _relerr(a, c):
    return float(np.abs(a - c).max() / max(np.abs(a).max(), 1e-300))


# ═══════════════════════════════════════════════════════════════════
# Probe 6 : fixed-point consistency (solve o apply = identity WITH correct seed)
# ═══════════════════════════════════════════════════════════════════
def probe_fixed_point_consistency(sig_mode):
    sn = _sphere(ng=2, bc="vacuum")
    rng = np.random.default_rng([7, 240])
    sig_t = rng.uniform(1.0, 3.0, size=(sn.ng, *sn.spatial_shape))
    if sig_mode == "removal":
        sig_r = sig_t - rng.uniform(0.2, 0.8, size=sig_t.shape)
    else:  # production
        sig_r = sig_t
    op = StreamingCollisionOperator(
        StreamingOperator(sn), MultiplicationOperator.from_mesh(sig_r, sn),
    )

    psi0 = _random_state(sn, seed=55)
    # b = A.apply(psi0). apply returns a timeless FullField; wrap into a
    # TimedFullField source (what the SI/Krylov driver does).
    b_src = op.apply(psi0)
    b = TimedFullField(
        bulk=b_src.bulk, boundary=b_src.boundary,
        _history=(), history_depth=psi0.history_depth,
    )

    # WRONG seed (cold): the test-suite round-trip.
    psi_cold = op.solve(b)
    # CORRECT seed (psi0 = the fixed point for this b).
    psi_fp = op.solve(b, initial_guess=psi0)

    e_cold = _relerr(psi0.bulk.values, psi_cold.bulk.values)
    e_fp = _relerr(psi0.bulk.values, psi_fp.bulk.values)
    print(f"\n=== Probe 6  sphere fixed-point consistency  [{sig_mode} sigma] ===")
    print(f"   solve(apply(psi0))            [cold seed ]  rel err vs psi0 = {e_cold:.3e}")
    print(f"   solve(apply(psi0), ig=psi0)   [correct sd]  rel err vs psi0 = {e_fp:.3e}")
    print("   -> if the 'correct seed' err is machine-precision, the sweep is the")
    print("      EXACT inverse at the fixed point: fixed point CORRECT, seed-lag only.")
    return e_cold, e_fp


# ═══════════════════════════════════════════════════════════════════
# Probe 8 : does a genuine sphere fixed-source SI converge? (physical source)
# ═══════════════════════════════════════════════════════════════════
def probe_real_si(with_scatter):
    """Run a hand SI: psi <- (L+C).solve(q_ext + Sigma_s0 * phi), threading the
    previous psi as the seed. Physical: sigma_t > sigma_s0 > 0, positive iso
    source. Confirms production sphere converges (or, with pure absorber,
    that the seed-lag alone still converges under a physical source)."""
    sn = _sphere(nx=10, ng=2, bc="vacuum")
    ng, nx = sn.ng, sn.nx
    sig_t = np.full((ng, nx), 1.0)
    sig_s0 = np.full((ng, nx), 0.6 if with_scatter else 0.0)
    sig_r = sig_t - sig_s0            # removal diagonal for (L+C)
    op = StreamingCollisionOperator(
        StreamingOperator(sn), MultiplicationOperator.from_mesh(sig_r, sn),
    )

    # Isotropic unit source, per-ordinate magnitude Q/W (producer applies /W).
    W = float(sn.quad.weights.sum())
    q_iso = 1.0
    from orpheus.transport.source_sinks import AngularSourceSink

    def build_rhs(phi_scalar):
        # within-group scatter source Sigma_s0 * phi, isotropic -> /W per ord.
        src_scalar = q_iso + sig_s0 * phi_scalar          # (ng, nx)
        per_ord = np.broadcast_to(
            (src_scalar / W)[None, :, :], (sn.quad.N, ng, nx),
        ).copy()
        rhs = TimedFullField.zeros(
            bulk=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn,
        )
        rhs.bulk.values[...] = per_ord[..., None] if rhs.bulk.values.ndim == 4 else per_ord
        return rhs

    phi = np.zeros((ng, nx))
    psi = None
    print(f"\n=== Probe 8  genuine sphere SI  (with_scatter={with_scatter}) ===")
    for k in range(60):
        rhs = build_rhs(phi)
        psi = op.solve(rhs, initial_guess=psi)
        # scalar flux from the sweep's returned angular flux
        w = sn.quad.weights
        av = psi.bulk.values
        av = av[..., 0] if av.ndim == 4 else av           # (N, ng, nx)
        phi_new = np.einsum("n,ngx->gx", w, av)
        dphi = float(np.abs(phi_new - phi).max() / max(np.abs(phi_new).max(), 1e-300))
        phi = phi_new
        if k < 6 or k % 10 == 0 or k == 59:
            print(f"   SI iter {k:2d}  max phi={phi.max():.4e}  ||dphi||_rel={dphi:.3e}")
        if dphi < 1e-11:
            print(f"   CONVERGED at iter {k} (||dphi||_rel={dphi:.2e})")
            break
    else:
        print("   did NOT converge in 60 iters")
    return phi


if __name__ == "__main__":
    probe_fixed_point_consistency("removal")
    probe_fixed_point_consistency("production")
    probe_real_si(with_scatter=False)
    probe_real_si(with_scatter=True)
