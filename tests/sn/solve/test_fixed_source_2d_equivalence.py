r"""Phase 1 (gap) — 2-D Cartesian fixed-source Krylov un-gate verification.

NEW-1 of the Phase 1 verification plan
(``.claude/agent-memory/test-architect/phase1_solve_unification_verification_plan.md``).

The Phase 1 "gap" sub-step (2026-06-04) deleted the
``NotImplementedError`` guard that made
:func:`~orpheus.sn.solver._solve_fixed_source_krylov` raise on 2-D
Cartesian meshes.  This file is the POSITIVE correctness pin that
REPLACES the retired raise-pin
(``test_fixed_source_g1.py::TestTwoDCartesianRaises``): instead of
pinning "2-D Krylov raises", we pin "2-D Krylov is CORRECT".

The un-gate is safe because the 2-D fixed-source Krylov is the structural
twin of two already-verified paths — the live 2-D *eigenvalue* Krylov
inner (:meth:`SNSolver._solve_krylov`) and the geometry-agnostic 2-D
*fixed-source SI* path — sharing the identical operator triple
(:func:`_within_group_triple`) and Krylov driver
(:func:`_within_group_krylov`), differing only in ``q_ext``.

Two legs (vv-principles §three-pillars):

* **Leg 1 — closed-form Q/Σ_t (the correctness anchor).** On a
  homogeneous, fully-reflective box with a uniform source and NO
  scattering (``placeholder_materials`` → ``Σ_s = 0``), the converged
  per-ordinate flux is the analytic streaming equilibrium
  ``ψ_n = q_n / Σ_t = (q_iso / W) / Σ_t`` (Signature 4).  This is
  structurally independent of any ORPHEUS solver — the infinite-medium
  balance.  ≥2G (per-group independent absorber) so the operand is
  non-degenerate.
* **Leg 2 — SI≡Krylov twin (necessary, not sufficient).** On a ≥2G
  HETEROGENEOUS NON-FLAT case (fuel | moderator split, vacuum-x), the
  two inners converge to the same solution.  Anchored by Leg 1's
  closed form; the degenerate-flux guard (``max/min > 1.2``) proves the
  redistribution terms are exercised (vv-principles anti-pattern #4).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem
from orpheus.geometry.mesh import Mesh2D
from orpheus.sn import solve_sn_fixed_source
from orpheus.sn.geometry import SNMesh
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.source_sinks import AngularSourceSink
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.l1


# ── Leg 1: closed-form Q/Σ_t streaming equilibrium ───────────────────


@pytest.mark.verifies("transport-cartesian")
def test_2d_homogeneous_reflective_krylov_hits_q_over_sigma_t() -> None:
    r"""2-D fixed-source Krylov converges to the analytic ``q_n / Σ_t``.

    Homogeneous, fully reflective, uniform source, ``Σ_s = 0`` ⟹ flat
    flux at the infinite-medium balance ``ψ_n = q_iso / (W · Σ_t)``
    (Σ_t = 1 here).  ≥2G (each group an independent absorber).  The
    closed form is structurally independent of the solver (L11).
    """
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, 5),
        edges_y=np.linspace(0.0, 2.0, 5),
        mat_map=np.zeros((4, 4), dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_xmin=BC("reflective"), bc_xmax=BC("reflective"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    # Genuine 2-D Cartesian reflective box ⟹ O_h symmetry. The SN-canonical
    # level-symmetric set is the right tool; Lebedev is the SO(3) moment
    # cubature (over-quadrature: O_h N=110 doe=17). level_symmetric(4) is the
    # SAME O_h group at N=24 doe=3. The asserted value q_iso/(W·Σ_t) is the
    # closed-form streaming equilibrium — structurally quadrature-independent
    # (W = Σ_n w_n = 4π for every rule; each ordinate equilibrates to q_n/Σ_t),
    # so the swap is exact. Verified: relerr 1.05e-12 ≤ rtol 1e-10.
    quad = Quadrature.level_symmetric(sn_order=4)
    materials = placeholder_materials(ng=2)  # Σ_t = 1, Σ_s = 0, ≥2G
    sn_mesh = SNMesh(mesh, quad, materials)
    sum_w = float(quad.weights.sum())

    q_iso = 0.5
    Q_iso = np.full((sn_mesh.ng, sn_mesh.nx, sn_mesh.ny), q_iso)
    src = AngularSourceSink.from_isotropic(Q_iso, sn_mesh)

    result = solve_sn_fixed_source(
        materials=materials, mesh=mesh, quadrature=quad,
        external_source=src.values,
        inner_solver="krylov",
        max_inner=300, inner_tol=1e-12,
    )
    assert result.history.converged, "2-D fixed-source Krylov did not converge"

    # Per-ordinate ψ uniform across (N, ng, nx, ny) at the fixed point;
    # Σ_t = 1.0 ⟹ ψ_n = q_iso / W.  Exact analytic limit → only FP noise.
    per_ord = result.angular_flux.bulk.values
    expected_per_ord = q_iso / sum_w
    np.testing.assert_allclose(
        per_ord, expected_per_ord, rtol=1e-10,
        err_msg=(
            "2-D fixed-source Krylov per-ord ψ does not match the "
            "closed-form streaming equilibrium q_iso/(W·Σ_t); a spike "
            "would indicate a missing volume/ΔA factor (failure mode 3)."
        ),
    )


# ── Leg 2: SI ≡ Krylov on a heterogeneous non-flat case ──────────────


@pytest.mark.verifies("transport-cartesian")
def test_2d_heterogeneous_si_krylov_equivalence() -> None:
    r"""2-D fixed-source SI and Krylov converge to the SAME solution.

    Fuel (A) | moderator (B) split across x with vacuum-x / reflective-y
    and a uniform source ⟹ a non-flat flux.  Both inners solve the
    identical ``(L + C − S − B) ψ = q`` operator equation, so they agree
    to iterative-solver tolerance.  Anchored by Leg 1's closed form.
    """
    mat_map = np.zeros((8, 4), dtype=int)
    mat_map[:4, :] = 2  # fuel half
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 4.0, 9),
        edges_y=np.linspace(0.0, 2.0, 5),
        mat_map=mat_map,
        coord=CoordSystem.CARTESIAN,
        bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    # SI-vs-Krylov EQUIVALENCE on a genuine 2-D Cartesian (8×4, fuel|mod,
    # vacuum-x/reflective-y) heterogeneous case ⟹ O_h symmetry. Both inners
    # solve the identical (L+C−S−B)ψ=q operator on the SAME quadrature, so the
    # equivalence holds for ANY quadrature; swap both sides together.
    # level_symmetric(4) (O_h, N=24 doe=3) replaces Lebedev (O_h, N=110 doe=17).
    # Verified: non-flat guard max/min=1.729 (>1.2 fires), SI≡Krylov to 1.2e-9
    # (within rtol 1e-6 / atol 1e-8). 6.4s → ~1.4s.
    quad = Quadrature.level_symmetric(sn_order=4)
    materials = {2: get_mixture("A", "2g"), 0: get_mixture("B", "2g")}
    sn_mesh = SNMesh(mesh, quad, materials)

    Q_iso = np.ones((sn_mesh.ng, sn_mesh.nx, sn_mesh.ny))
    src = AngularSourceSink.from_isotropic(Q_iso, sn_mesh)

    common = dict(
        materials=materials, mesh=mesh, quadrature=quad,
        external_source=src.values, max_inner=500, inner_tol=1e-10,
    )
    sol_si = solve_sn_fixed_source(inner_solver="source_iteration", **common)
    sol_kry = solve_sn_fixed_source(inner_solver="krylov", **common)

    phi_si = sol_si.scalar_flux.values
    phi_kry = sol_kry.scalar_flux.values

    # Degenerate-flux guard: the profile must be genuinely non-flat so the
    # redistribution / streaming terms are exercised (anti-pattern #4).
    prof = phi_si[0]
    assert prof.max() / prof.min() > 1.2, (
        f"flux profile is too flat (max/min = {prof.max() / prof.min():.3f}); "
        f"the heterogeneous case must exercise redistribution."
    )

    # Both inners solve the identical operator equation → same solution
    # (magnitude AND shape) to iterative-solver tolerance.
    np.testing.assert_allclose(
        phi_kry, phi_si, rtol=1e-6, atol=1e-8,
        err_msg=(
            "2-D fixed-source SI and Krylov disagree beyond iterative "
            "tolerance — they share the operator triple, so divergence "
            "signals a convention/sign drift between the two drivers."
        ),
    )
