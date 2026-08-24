r"""ERR-058 catcher — the curvilinear operator ADMITS the manufactured solution.

Promoted 2026-06-12 from ``derivations/diagnostics/diag_195_probe3_residual_audit.py``
(numerics-investigator, Issue #195) after the ERR-058 closure-seed fix.

The decisive structural check behind the #195 diagnosis: build the
typed within-group triple ``(L+C, S, B)`` exactly as
``solve_sn_fixed_source`` does, apply it to the EXACT manufactured
angular flux :math:`\psi_{\text{ref},n} = A(r)/W`, and demand the
**per-ordinate** balance residual

.. math::

   r_{n,g,i} \;=\; \bigl[(L + C - S - B)\,\psi_{\text{ref}}\bigr]_{n,g,i}
                   \;-\; Q^{\text{ext}}_{n,g,i}

decay under spatial refinement in the volume-weighted norm.

Why PER-ORDINATE and why VOLUME-WEIGHTED
-----------------------------------------

* **Per-ordinate**: the M-M angular redistribution α-dome telescopes
  under the angular weight sum, so a SCALAR (weight-summed) residual is
  structurally blind to half-angle-thread defects — exactly how
  ERR-058 b (the Carlson proxy-source seed, per-ordinate residual
  O(10)) hid behind a scalar residual of O(h²) during the #195
  investigation (vv-principles anti-pattern #8, instantiated in a
  diagnostic).
* **Volume-weighted**: the sphere's pole-adjacent cells legitimately
  carry a bounded non-decaying pointwise residual (the closure
  truncation meets the :math:`\Delta A/V \sim 1/h` geometry factor on
  cells whose volume vanishes as :math:`r^2\,dr`); the solution-error
  ladders prove it harmless.  The volume-weighted norm is the
  stability-relevant measure: sphere decays at order ≈ 1.5, cylinder
  at order ≈ 2.0 (measured 2026-06-12 post-fix).

Bug-era values this gate must catch (pre-ERR-058, measured): the
volume-weighted per-ordinate residual was O(1e-1)-class
(per-ordinate pointwise up to ±55 at the pole, ±13 in the bulk) —
3+ orders of magnitude above the asserted bounds.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import (
    build_cylindrical_mms_case,
    build_spherical_mms_case,
)
from orpheus.sn.solver import (
    SNSolver,
    _as_sn_mesh,
    _build_fixed_source_rhs,
)
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import radial_characteristic_edge_seed


def _vol_weighted_per_ordinate_residual(case, nc: int) -> float:
    """RMS volume-weighted per-ordinate residual of ψ_ref under (L+C−S−B)."""
    mesh = case.build_mesh(nc)
    Q = case.external_source(mesh)
    sn_mesh = _as_sn_mesh(
        mesh, case.quadrature, case.materials, "vacuum", mat_map=None,
    )
    solver = SNSolver(
        sn_mesh, inner_solver="source_iteration",
        scattering_order=0, max_inner=10, inner_tol=1e-13,
    )
    q_ext = _build_fixed_source_rhs(Q, sn_mesh)
    # B.2d: on a carrying mesh the rhs builder returns the coupled pair —
    # this probe reads System A's member (the q½ member is A_BB's rhs).
    if hasattr(q_ext, "systems"):
        q_ext = q_ext.systems[0]
    # B.2d: the triple retired into build_within_group_system; this fused
    # 3-block probe reads the production surfaces directly. B = B_a alone is
    # bit-identical here: on vacuum cases the ray-corner B_b term is exactly
    # zero, and B_a pads the ray slot present-zero like the retired composite.
    from orpheus.sn.coupled_system import build_streaming_collision
    from orpheus.sn.operators.boundary import SNBoundaryOperator
    LC = build_streaming_collision(solver.sn_mesh, solver.mat_xs)
    S = solver.scattering_op
    B = SNBoundaryOperator(solver.sn_mesh)

    # ψ_ref,n = A(r)/W — isotropic per ordinate, zero boundary (vacuum).
    A = case.phi_exact(mesh.centers)
    W = float(sn_mesh.quad.weights.sum())
    vals = np.zeros((sn_mesh.quad.N, sn_mesh.ng, *sn_mesh.spatial_shape))
    vals[:, 0, :] = (A / W)[None, :]
    zero = TimedFullField.zeros(
        interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn_mesh.full_field_space,
    )
    psi_ref = TimedFullField(
        interior=AngularFlux.from_mesh(vals, sn_mesh), boundary=zero.boundary,
    )
    # #282 route (a) → B.2d: the CONSISTENT edge-extrapolated ψ½ seed of the
    # MMS trial (the trial's own μ = −1 starting datum) rides the walk's
    # EXPLICIT flux leg, so LC.apply reproduces the pre-route-(a) operator
    # action and the residual decays as before; no legs on non-carrying.
    seed_leg = radial_characteristic_edge_seed(vals, sn_mesh)
    if seed_leg is None:
        lc_out = LC.apply(psi_ref)
    else:
        # step 6: the joint row-A action rides THE GRID (presence structural).
        from orpheus.numerics.coupled_system import CoupledField

        from tests.sn._test_helpers import joint_m_grid

        grid, _space = joint_m_grid(sn_mesh, LC)
        lc_out = grid.apply(
            CoupledField(systems=(psi_ref, seed_leg)),
        ).systems[0]

    rv = (
        lc_out.interior.values
        - S.apply(psi_ref).interior.values
        - B.apply(psi_ref).interior.values
        - q_ext.interior.values
    )
    V = sn_mesh.volumes
    return float(np.sqrt(np.einsum("x,ngx->", V, rv**2) / V.sum()))


@pytest.mark.l1
@pytest.mark.catches("ERR-058")
@pytest.mark.parametrize(
    "builder, name, min_order, abs_bound_80",
    [
        # Measured post-fix (2026-06-12): sphere 1.97e-3 → 9.7e-4 →
        # 3.4e-4 (orders ~1.5 — the pole-adjacent O(1) band under the
        # r²dr volume weight); cylinder 5.50e-5 → 1.37e-5 → 3.44e-6
        # (orders 2.0).  Bug-era values were O(1e-1)-class.
        pytest.param(
            build_spherical_mms_case, "sphere", 1.3, 2.0e-3, id="sphere",
        ),
        pytest.param(
            build_cylindrical_mms_case, "cylinder", 1.8, 3.0e-5, id="cylinder",
        ),
    ],
)
def test_operator_admits_isotropic_mms_per_ordinate(
    builder, name, min_order, abs_bound_80,
):
    """(L+C−S−B)·ψ_ref − Q decays per-ordinate (volume-weighted)."""
    case = builder()
    r40 = _vol_weighted_per_ordinate_residual(case, 40)
    r80 = _vol_weighted_per_ordinate_residual(case, 80)
    order = float(np.log2(r40 / r80))
    assert order > min_order, (
        f"{name}: volume-weighted per-ordinate residual does not decay "
        f"(r40={r40:.3e}, r80={r80:.3e}, order={order:.2f}) — the "
        f"discrete operator no longer admits the manufactured solution "
        f"per ordinate (ERR-058 closure-seed regression class)."
    )
    assert r80 < abs_bound_80, (
        f"{name}: r80={r80:.3e} exceeds {abs_bound_80:.1e} — bug-era "
        f"values were O(1e-1); a re-floored closure seed lands here."
    )
