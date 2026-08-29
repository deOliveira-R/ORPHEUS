r"""L1 MMS — Linear-Discontinuous (LD) SN on a 1-D slab, via the affine scan.

LD is the first non-Diamond-Difference spatial scheme (#158).  Since Increment B
it is **affine-scannable** (the Schur complement eliminates the slope, leaving
the single-upstream recurrence ``ψ_out = a·ψ_in + b``), so a 1-D slab LD mesh
rides the fast DAG-free ``CumprodScan`` via the coefficient model — it supplies
``(a, inverse_denom, w)`` through ``affine_scan_coefficients`` and the generic
base reconstruction staticmethods (``source_emission`` / ``cell_average`` /
``outgoing_face_from_average`` on ``DiscretizationSchemeBase``) do the rest.  The
polymorphic ``FullFieldWavefront`` DAG
oracle (Increment A, via ``cell_kernel_batch``) remains the verification
reference — the two-paths gate pins ``CumprodScan``-LD ≡ ``FullFieldWavefront``-LD.
This file pins:

* **Routing** — a 1-D slab LD mesh dispatches to ``CumprodScan`` by
  ``default_for`` (and a DD mesh still picks ``CumprodScan`` too — both affine,
  no selection regression).
* **Spatial order** — LD converges :math:`\mathcal{O}(h^2)` on the manufactured
  slab problem, run end-to-end through
  ``solve_sn_fixed_source(..., scheme=LinearDiscontinuous())``.
* **Two-paths (#158 Inc B)** — ``CumprodScan``-LD ≡ ``FullFieldWavefront``-LD
  (principled-equivalent at nULP; the ×V scan vs ÷V kernel conventions agree).
* **Matvec ≡ sweep** — the Krylov apply (LD's group-2 ``residual_kernel_batch``
  via ``_apply_walk``) matches the SI sweep.

The per-cell exactness-on-linears (the strong structurally-independent
correctness oracle) lives in ``tests/transport/spatial/test_linear_discontinuous.py``.

Related: ``orpheus.transport.spatial.linear_discontinuous`` (the occupant);
``orpheus.transport.spatial.scheme.DiscretizationSchemeBase`` (the generic coefficient-model
reconstruction staticmethods);
``orpheus.sn.loss_representation.{CumprodScan, FullFieldWavefront}``;
``.claude/plans/mellow-swinging-breeze.md`` (Increment B).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import build_1d_slab_mms_case
from orpheus.sn import solve_sn_fixed_source
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.loss_representation import (
    CumprodScan,
    FullFieldWavefront,
    default_for,
)
from orpheus.transport.spatial import LinearDiscontinuous

from tests.sn._test_helpers import volume_weighted_l2


@pytest.mark.foundation
def test_ld_slab_mesh_routes_to_cumprod_scan() -> None:
    """A 1-D slab LD mesh dispatches to the fast affine scan (#158 Inc B); DD
    stays on the scan too (both affine — no selection regression)."""
    case = build_1d_slab_mms_case()
    mesh = case.build_mesh(16)
    ld_mesh = SNMesh(mesh, case.quadrature, case.materials,
                     scheme=LinearDiscontinuous())
    dd_mesh = SNMesh(mesh, case.quadrature, case.materials)
    if not isinstance(default_for(ld_mesh, ld_mesh.scheme, ld_mesh.angular_closure), CumprodScan):
        pytest.fail(
            "LD slab mesh routed to "
            f"{type(default_for(ld_mesh, ld_mesh.scheme, ld_mesh.angular_closure)).__name__}, expected CumprodScan "
            "(LD is affine-scannable since Increment B)"
        )
    if not isinstance(default_for(dd_mesh, dd_mesh.scheme, dd_mesh.angular_closure), CumprodScan):
        pytest.fail(
            "DD slab mesh no longer routes to CumprodScan "
            f"({type(default_for(dd_mesh, dd_mesh.scheme, dd_mesh.angular_closure)).__name__}) — selection regression"
        )


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("ld-ubld-d1-reduction", "transport-cartesian")
def test_sn_1d_slab_ld_mms_converges_second_order() -> None:
    r"""LD SN on a 1-D slab shows measured :math:`\mathcal{O}(h^2)` end-to-end.

    Drives the manufactured slab problem through
    ``solve_sn_fixed_source(scheme=LinearDiscontinuous())`` — which routes
    to ``FullFieldWavefront`` calling LD's ``cell_kernel_batch`` — at four mesh
    refinements.  Asserts the order converges to 2 (the finest pair > 1.95; the
    pre-asymptotic coarse pair > 1.85) AND the converged value (the magnitude
    band against the manufactured ``phi_exact`` — rate alone is necessary, never
    sufficient, vv-principles §5).
    """
    case = build_1d_slab_mms_case()
    n_cells = [20, 40, 80, 160]
    errors = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        Q = case.external_source(mesh)
        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, Q,
            boundary_condition="vacuum", max_inner=500, inner_tol=1e-13,
            scheme=LinearDiscontinuous(),
        )
        phi_num = result.scalar_flux.values[0, :]
        phi_ref = case.phi_exact(mesh.centers)
        errors.append(volume_weighted_l2(phi_num, phi_ref, mesh.widths))

    errors = np.asarray(errors)
    orders = np.log2(errors[:-1] / errors[1:])
    # O(h²) converging to 2: finest pair clears 1.95, all pairs clear 1.85
    # (the smooth sin ansatz is mildly pre-asymptotic on the coarse pair; the
    #  strong exactness-on-linears oracle is the per-cell gate-1 companion).
    assert orders[-1] > 1.95, f"LD finest order {orders[-1]:.3f} ≯ 1.95 ({orders})"
    assert np.all(orders > 1.85), f"LD order dropped below 1.85: {orders}"
    assert 1e-9 < errors[-1] < 1e-2, (
        f"LD finest-mesh error {errors[-1]:.3e} outside the value band — "
        "rate may be O(h²) to the WRONG limit (vv-principles §5)"
    )


@pytest.mark.l1
@pytest.mark.verifies("ld-ubld-d1-reduction")
def test_sn_1d_slab_ld_mms_krylov_matches_si() -> None:
    r"""LD forward matvec ≡ SI sweep (the L14 matvec≡sweep leg).

    Since Increment B LD routes to ``CumprodScan``: the Krylov inner exercises
    LD's ``residual_kernel_batch`` via the scan walk's ``_apply_walk`` cartesian
    arm (the apply twin of the SI sweep's group-3 ``affine_scan_coefficients`` +
    generic ``cell_average``); both must converge to the SAME manufactured flux.
    """
    case = build_1d_slab_mms_case()
    mesh = case.build_mesh(80)
    Q = case.external_source(mesh)
    si = solve_sn_fixed_source(
        case.materials, mesh, case.quadrature, Q, boundary_condition="vacuum",
        max_inner=500, inner_tol=1e-13, scheme=LinearDiscontinuous(),
    )
    kry = solve_sn_fixed_source(
        case.materials, mesh, case.quadrature, Q, boundary_condition="vacuum",
        max_inner=500, inner_tol=1e-13, inner_solver="krylov",
        scheme=LinearDiscontinuous(),
    )
    np.testing.assert_allclose(
        kry.scalar_flux.values[0], si.scalar_flux.values[0],
        rtol=1e-9, atol=1e-11,
    )


@pytest.mark.l1
@pytest.mark.verifies("ld-ubld-d1-reduction")
def test_ld_two_paths_scan_equals_dag_oracle() -> None:
    r"""#158 Inc B headline — ``CumprodScan``-LD ≡ ``FullFieldWavefront``-LD.

    The fast affine scan (``CumprodScan``, group-3 ``affine_scan_coefficients``
    in the ×V convention) and the polymorphic DAG oracle
    (``FullFieldWavefront``, group-2 ``cell_kernel_batch`` in the ÷V ``g=|μ|/Δ``
    convention) solve the SAME LD discrete system via two different schedules
    AND two different scalings.  They are **principled-equivalent** (NOT
    byte-identical: ``S_scan = V·S_kernel`` + distinct reduction trees), so the
    two-paths cross-check pins LD's affine coefficients against the trusted
    Increment-A kernel — the system-level analogue of the unit group2≡group3
    gate, and the slope-row sign-trap catcher at the schedule level.

    STRESSING config (Mode-9 — a flat/homogeneous two-paths gate would be
    FP-coincidentally exact): heterogeneous per-cell-per-group ``σ_t`` (≥2G) +
    a non-flat random per-ordinate source.  One sweep each, on a fresh boundary
    (the sweep mutates it in place).
    """
    from orpheus.geometry import BC, CoordSystem, Mesh1D
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
    from orpheus.derivations.continuous.mms.sn import _make_2g_asymmetric_mixture

    nx = 24
    materials = {0: _make_2g_asymmetric_mixture(
        np.array([1.0, 1.5]), np.array([[0.3, 0.2], [0.0, 0.6]]),
    )}
    mesh = Mesh1D(
        edges=np.linspace(0.0, 5.0, nx + 1), mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN, bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(8)
    ld_mesh = SNMesh(mesh, quad, materials, scheme=LinearDiscontinuous())
    N, ng = quad.N, ld_mesh.ng

    rng = np.random.default_rng(20260614)
    sig_t = rng.uniform(0.3, 3.0, size=(ng, nx))          # het, ≥2G
    Q = rng.standard_normal((N, ng, nx))                  # non-flat per-ordinate

    _, phi_scan = CumprodScan.pose(ld_mesh).sweep(
        Q, sig_t, AngularBoundaryFlux.zeros(ld_mesh.angular_trace),
    )
    _, phi_dag = FullFieldWavefront.pose(ld_mesh).sweep(
        Q, sig_t, AngularBoundaryFlux.zeros(ld_mesh.angular_trace),
    )
    # Principled-equivalent (×V scan vs ÷V kernel): tight nULP-scale band, well
    # below any algorithmic-difference signature (a sign-trap would be O(1)).
    np.testing.assert_allclose(phi_scan, phi_dag, rtol=1e-12, atol=1e-13)


@pytest.mark.foundation
def test_ld_curvilinear_solve_fails_fast() -> None:
    """A 1-D *curvilinear* LD solve must fail FAST, not silently run
    DD-shaped curvature math.

    The claim is the end-to-end one: such a mesh matches
    ``CumprodScan.supports`` (``is_1d and is_affine_scannable``), so without
    a refusal somewhere it would sweep with the wrong math.

    ⚠ **The refusal's PROVENANCE moved on 2026-08-26 and this test moved with
    it.**  It used to assert ``match="slab/Cartesian"`` — the LD *scan
    closure* guard (#158 Inc B), which fires at the ``CollisionCache`` build
    in ``SNSolver.__init__``.  P1 item 9 added an EARLIER guard on the moment
    MASS: an LD curvilinear space cannot even be posed, because the true
    ``M/V`` there is cell-dependent and non-diagonal.  That is the deeper
    cause — you cannot pose the space, let alone scan it — so the earlier,
    more specific refusal is the right one to surface.

    ⭐ The scan-closure guard is NOT orphaned by that move: this test was its
    only witness, so a DIRECT one was written in the same commit —
    ``test_ld_scan_closure_refuses_non_neutral_curvature`` below.  That is
    strictly better evidence, because it exercises the guard's actual
    predicate (a *value* signal — P4.9a: the assembled ``angular_denom_term``
    non-neutral) instead of a path that merely happened to reach it.
    """
    from orpheus.geometry import BC, CoordSystem, Mesh1D
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.derivations.continuous.mms.sn import _make_1g_mixture

    nx = 8
    materials = {0: _make_1g_mixture(1.0, 0.5)}
    sphere = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1), mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL, bc_left=BC("reflective"), bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(8)
    Q = np.ones((quad.N, 1, nx)) / quad.weights.sum()
    with pytest.raises(NotImplementedError, match="no moment mass on a spherical"):
        solve_sn_fixed_source(
            materials, sphere, quad, Q, boundary_condition="vacuum",
            scheme=LinearDiscontinuous(),
        )


@pytest.mark.foundation
def test_ld_scan_closure_refuses_non_neutral_curvature() -> None:
    """The #158 scan-closure guard, pinned at its OWN predicate.

    Written 2026-08-26 because the end-to-end test above stopped being this
    guard's witness when an earlier refusal was added, and it was the only
    one (`[M]` ``grep "slab/Cartesian" tests/`` returned exactly one hit).  A
    guard nothing can redden is a coverage claim an audit will trust —
    ``vv-principles`` #17.

    The guard keys on a **value** signal, not a chart tag: slab carries an
    exactly-zero assembled ``angular_denom_term`` (the identity closure's
    zero constants over ``slab_streaming``'s neutral element) and
    curvilinear does not.  So it is reachable directly, with no mesh at
    all — which no earlier guard can preempt.

    ⚠ P4.9a re-key, stated so nobody reads more than it checks: the former
    two-signal predicate (``dA_w`` / ``c_out`` separately) collapsed to the
    ASSEMBLED product, so a ``dA_w ≠ 0`` with ``c_out ≡ 0`` configuration
    (an identity closure hand-mounted on a curvilinear mesh) no longer trips
    THIS guard — it is refused upstream by the walk admission
    (``supports_curvilinear=False``), the primary guard.

    Both legs, per ``vv-principles`` #11: the neutral input is ACCEPTED and
    the non-neutral input is REFUSED.  Without the positive leg this is
    satisfied by a guard that rejects everything.
    """
    ld = LinearDiscontinuous()
    kw = dict(
        abs_mu=np.array([0.5]),
        A_down=np.array([[1.0]]),
        A_total=np.array([[2.0]]),
        V=np.array([[0.25]]),
        reaction_xs=np.array([[[1.3]]]),
    )

    # POSITIVE: the slab's neutral element is admitted.
    a, inv, w = ld.affine_scan_coefficients(
        **kw, angular_denom_term=np.array([[0.0]]),
    )
    assert np.all(np.isfinite(a)) and np.all(np.isfinite(inv))
    assert np.all(np.isfinite(w))

    # NEGATIVE: a non-neutral assembled contribution is refused.
    with pytest.raises(NotImplementedError, match="slab/Cartesian"):
        ld.affine_scan_coefficients(
            **kw, angular_denom_term=np.array([[0.21]]),
        )


@pytest.mark.l1
@pytest.mark.verifies(
    "ld-ubld-d1-reduction", "ld-ubld-slope-angular-reduction",
    "ld-ubld-octant-moment-frame-signs",
)
@pytest.mark.catches("ERR-061")
def test_ld_thick_diffusive_limit() -> None:
    r"""Thick diffusive slab: LD recovers the diffusion limit (≈ DD interior).

    On an optically thick, scattering-dominated mesh (σ_t·h = 10/cell, c→1) DD
    reproduces the interior diffusion solution (LMM-1987 Table I); the FULL LD
    with the threaded slope source Σ_s·φ̂ (Increment C) does too — LD ≈ DD at the
    midpoint at the COARSE mesh (NOT only under refinement).

    Increment C threads the spatial-slope iterate φ̂.  #240 D5b-S3 fixed the
    slope-moment FRAME (ERR-061): the per-ordinate slope ψ̂_n is stored in the
    GLOBAL-x frame so the angular reduction φ̂ = Σ_n w_n ψ̂_n (which feeds the
    isotropic scattering slope source Σ_s·φ̂) does NOT cancel forward against
    backward ordinates.  Pre-fix LD was 38.9% off DD at nx=4 (the flat-source
    signature persisting through the under-driven slope); post-fix it is < 5%.

    Reference: the continuous diffusion solution (LMM-1987 §; DD-as-anchor is the
    diffusion proxy — DD's 1-D diffusion-limit consistency is the
    separately-published ground, vv §three-pillars).  The 2G-het group-coupled
    companion is ``test_ld_thick_diffusive_limit_2g`` (Mode-6 — the slope source
    Σ_s^T·φ̂ couples groups; a 1G-only gate is a degeneracy guard failure).
    """
    from orpheus.derivations.continuous.mms.sn import _make_1g_mixture
    from orpheus.geometry import BC, CoordSystem, Mesh1D
    from orpheus.numerics.quadrature import Quadrature

    nx, length, sigma_t, c = 4, 1.0, 40.0, 0.99   # σ_t·h = 10/cell (thick), c→1
    materials = {0: _make_1g_mixture(sigma_t, c * sigma_t)}
    mesh = Mesh1D(
        edges=np.linspace(0.0, length, nx + 1), mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN, bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(8)
    Q = np.ones((quad.N, 1, nx)) / quad.weights.sum()   # uniform iso source
    dd = solve_sn_fixed_source(
        materials, mesh, quad, Q, boundary_condition="vacuum",
        inner_solver="krylov", max_inner=2000, inner_tol=1e-10,
    )
    ld = solve_sn_fixed_source(
        materials, mesh, quad, Q, boundary_condition="vacuum",
        inner_solver="krylov", max_inner=2000, inner_tol=1e-10,
        scheme=LinearDiscontinuous(),
    )
    mid = nx // 2
    dd_mid = float(dd.scalar_flux.values[0, mid])
    ld_mid = float(ld.scalar_flux.values[0, mid])
    rel = abs(ld_mid - dd_mid) / abs(dd_mid)
    # Mode-8: a function-call assertion (fires under -O), NOT a bare assert.
    np.testing.assert_array_less(
        rel, 0.05,
        err_msg=(
            f"LD midpoint {ld_mid:.4f} vs DD {dd_mid:.4f} (rel {rel:.1%}) — "
            "the LD diffusion limit needs the global-frame slope moment "
            "(#240 D5b-S3 / ERR-061; the scattering slope source Σ_s·φ̂)"
        ),
    )


@pytest.mark.l1
@pytest.mark.verifies(
    "ld-ubld-d1-reduction", "ld-ubld-slope-angular-reduction",
    "ld-ubld-octant-moment-frame-signs",
)
@pytest.mark.catches("ERR-061")
def test_ld_thick_diffusive_limit_2g() -> None:
    r"""2G-het companion of the thick-diffusion limit (Mode-6 — group coupling).

    The slope-scattering source ``Σ_s^T·φ̂`` is group-coupled; a 1G test cannot
    see a ``Σ_s^T`` transpose / convention bug on the slope rows (1G Σ_s is
    scalar).  This drives a 2-group asymmetric-downscatter thick slab and asserts
    LD ≈ DD at the coarse mesh in BOTH groups — the non-degenerate companion to
    the 1G gate (the `vv-principles` 1-group-degeneracy rule / H1; #240 D5b-S3 GATE 5).
    """
    from orpheus.derivations.continuous.mms.sn import _make_2g_asymmetric_mixture
    from orpheus.geometry import BC, CoordSystem, Mesh1D
    from orpheus.numerics.quadrature import Quadrature

    nx, length = 4, 1.0
    # Thick, scattering-dominated, ASYMMETRIC 2G downscatter (g0→g1 off-diagonal).
    sig_t = np.array([40.0, 40.0])
    sig_s = np.array([[30.0, 9.6], [0.0, 39.6]])   # c≈0.99 each group, asym
    materials = {0: _make_2g_asymmetric_mixture(sig_t, sig_s)}
    mesh = Mesh1D(
        edges=np.linspace(0.0, length, nx + 1), mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN, bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(8)
    Q = np.ones((quad.N, 2, nx)) / quad.weights.sum()
    common = dict(boundary_condition="vacuum", inner_solver="krylov",
                  max_inner=4000, inner_tol=1e-10)
    dd = solve_sn_fixed_source(materials, mesh, quad, Q, **common)
    ld = solve_sn_fixed_source(materials, mesh, quad, Q,
                               scheme=LinearDiscontinuous(), **common)
    mid = nx // 2
    for g in range(2):
        dd_mid = float(dd.scalar_flux.values[g, mid])
        ld_mid = float(ld.scalar_flux.values[g, mid])
        rel = abs(ld_mid - dd_mid) / abs(dd_mid)
        np.testing.assert_array_less(
            rel, 0.06,
            err_msg=(
                f"group {g}: LD {ld_mid:.4f} vs DD {dd_mid:.4f} (rel {rel:.1%}) "
                "— 2G group-coupled slope source (#240 D5b-S3 / ERR-061, Mode-6)"
            ),
        )
