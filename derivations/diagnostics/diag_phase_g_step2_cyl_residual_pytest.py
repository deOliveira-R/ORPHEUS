"""Diagnostic (PROMOTABLE): cylinder apply-matvec preserves flat ψ.

Created by numerics-investigator on 2026-05-13.
If this test catches a real bug, promote to tests/sn/test_apply_matvec_cylinder_invariants.py

The cylinder transport apply-matvec L on ψ_flat (constant) MUST yield
Σ_t·ψ_flat at every (ordinate, cell): streaming and angular redistribution
must telescope to zero on flat-flux input (the Bailey α-dome invariant +
Pomraning structural-singularity preservation). Any deviation is a
structural defect in the curvilinear operator.

Pre-fix (production 2026-05-13, MorelMontryAngularSweep + sphere Σw=2
hardcoded Carlson seed factor): residual = 4-17 (mesh-dependent, grows
with refinement) → DIVERGENT. Catches ERR-048 manifestation #3.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.quadrature import ProductQuadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operator import (
    build_equation_map_cylindrical,
    transport_operator_matvec_cylindrical,
)


@pytest.mark.parametrize("n_cells", [10, 20, 40])
@pytest.mark.parametrize("n_mu", [2, 4])
@pytest.mark.parametrize("n_phi", [2, 4])
def test_cylinder_apply_matvec_preserves_flat_psi(n_cells, n_mu, n_phi):
    """Cylinder apply-matvec must satisfy L·ψ_flat = Σ_t·ψ_flat exactly."""
    fuel = get_mixture("B", "1g")
    R = 2.0
    geom = StructuredGeometry(
        geometry="CYL",
        regions=(Region(mat_id=0, outer_thickness_cm=R),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(
        geom, region_meshes=(RegionMesh(n_cells=n_cells),),
    )
    quad = ProductQuadrature.create(n_mu=n_mu, n_phi=n_phi)
    sn_mesh = SNMesh(mesh, quad)
    nx = n_cells
    ng = 1
    sig_t = np.full((nx, 1, ng), 2.0)
    psi_flat = 10.0 / quad.weights.sum()
    eq_map = build_equation_map_cylindrical(nx, quad, ng)
    sol_flat = np.full(eq_map.n_unknowns, psi_flat)

    lhs = transport_operator_matvec_cylindrical(
        sol_flat, eq_map, quad, sig_t, nx, ng,
        sn_mesh.reduced.face_areas, sn_mesh.volumes,
        sn_mesh.reduced.alpha_per_level,
        sn_mesh.reduced.redist_dAw_per_level,
        sn_mesh.reduced.tau_mm_per_level,
        sn_mesh=sn_mesh,
        pole_angular_closure=sn_mesh.pole_angular_closure,
    )

    expected = 2.0 * psi_flat  # = Σ_t · ψ_flat
    residual = np.max(np.abs(lhs - expected))
    assert residual < 1e-10, (
        f"Cylinder apply-matvec residual on flat ψ: {residual:.6e}, "
        f"expected < 1e-10. Configuration n_cells={n_cells}, n_mu={n_mu}, "
        f"n_phi={n_phi}, Σw={quad.weights.sum():.6f}, ψ_flat={psi_flat:.6f}. "
        f"This is the L0 per-ordinate flat-flux invariant — Bailey's α-dome "
        f"telescoping + collision = Σ_t·ψ. Pre-fix the residual grows with "
        f"mesh refinement (Signature 1 fingerprint, ERR-048 manifestation #3)."
    )
