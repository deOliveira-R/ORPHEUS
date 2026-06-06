r"""Phase 4 / O.2b 4.6 — NON-VACUUM prescribed-inflow MMS convergence (L1).

These rows verify the SN spatial + angular + **prescribed-inflow**
discretisation on a fixed-source problem whose exact angular flux is
known in closed form AND is NON-ZERO at the boundary.  The entire
existing MMS catalog is vacuum-automatic (every ansatz vanishes at the
faces, so :math:`\gamma_-\psi \equiv 0`); 4.6 fills the
:math:`q.\text{boundary} \neq 0` gap.

These rows drive the PUBLIC composite-source API:
``solve_sn_fixed_source`` accepts the composite ``q = q_bulk ⊕ q_∂`` (a
:class:`~orpheus.transport.timed_full_field.TimedFullField`), so the
non-vacuum prescribed inflow is supplied as
``case.fixed_source(sn).boundary`` — a
:class:`~orpheus.transport.source_sinks.BoundarySourceSink` built by the
ergonomic ``BoundarySourceSink.prescribed_inflow`` generator. The
affine-BC inhomogeneous term ``q`` is consumed by ``(L+C).solve`` as the
sweep inflow seed.  (These rows originally drove the ``(L+C),S,B``
operator-triple bypass because ``solve_sn_fixed_source`` hardcoded vacuum
``q.boundary``; once the public API gained composite-source support they
were migrated onto it — retirement = test migration.)

**The load-bearing assertion is the converged VALUE, not the rate.**
A dropped ``q.boundary`` (a silent vacuum solve) converges cleanly at
O(h²) to a DIFFERENT, boundary-ZERO flux — passing the rate gate.  Only
the value-vs-:math:`A(x)` (a0>0, non-zero at the boundary) assertion
sees the dropped inflow (``vv-principles`` §5 / H4: rate is
necessary-not-sufficient).

Mode-7 activates/nulls (mandatory declaration):
- **Slab ACTIVATES**: streaming :math:`\mu(A'+\mu B')/W`, within-group
  scatter c<1, :math:`\mu^2` second-moment coupling, NON-VACUUM
  :math:`\gamma_-\psi` (the novelty), 2G ``SigSᵀ`` group transfer (T2).
- **Slab NULLS**: curvilinear angular redistribution (→ the sphere
  companion ``test_mms_prescribed_inflow_sphere_activates_redistribution``
  is mandatory; NEVER ship slab-only — ERR-026 territory) + fission
  (non-fissile mixture; MMS proves no eigenvalue).

See:
- ``.claude/skills/vv-principles/SKILL.md`` (failure mode #7, H4).
- ``.claude/plans/phase4_o2b_46_nonvacuum_mms_plan.md`` (decision E — bypass).
- :mod:`orpheus.derivations.continuous.mms.sn` (the Branch-1/Branch-2 cases).
- ``docs/theory/discrete_ordinates.rst`` labels ``sn-mms-nonvacuum-*``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import (
    build_slab_2g_nonvacuum_mms_case,
    build_slab_nonvacuum_mms_case,
    build_sphere_nonvacuum_mms_case,
)
from orpheus.sn import solve_sn_fixed_source
from orpheus.sn.geometry import SNMesh


# ── shared solve helper (the PUBLIC composite-source API) ────────────


def _l2_error(phi_num: np.ndarray, phi_ref: np.ndarray, measure: np.ndarray) -> float:
    """Measure-weighted discrete L2 norm of the flux error.

    Slab uses ``mesh.widths``; curvilinear uses ``mesh.volumes`` (the
    cross-domain-attacker H-vol note: the L2 norm MUST be volume-weighted
    for the curvilinear order to be measured correctly).
    """
    diff = phi_num - phi_ref
    return float(np.sqrt(np.sum(measure * diff * diff)))


def _solve(case, n_cells: int, *, g: int = 0):
    r"""Solve the non-vacuum fixed-source problem via the PUBLIC composite
    API: :func:`~orpheus.sn.solve_sn_fixed_source` consuming
    ``case.fixed_source(sn)`` — the bulk ⊕ prescribed-inflow
    :class:`~orpheus.transport.timed_full_field.TimedFullField`.

    This is the migration off the operator-triple bypass (the public API
    now carries non-vacuum sources via the composite RHS; retirement =
    test migration). The slab uses SI (1-D Jacobi); the sphere uses the
    curvilinear Krylov default — both converge to the same fixed point.
    Returns ``(mesh, sn, psi, phi_g)`` where ``phi_g`` is the scalar flux
    on the radial slice (ny=0).
    """
    mesh = case.build_mesh(n_cells)
    materials = (
        case.build_materials(mesh)
        if hasattr(case, "build_materials")
        else case.materials
    )
    sn = SNMesh(mesh, case.quadrature, materials)
    result = solve_sn_fixed_source(
        materials, mesh, case.quadrature, case.fixed_source(sn),
        max_inner=1000, inner_tol=1e-13,
    )
    psi = result.angular_flux
    phi = result.scalar_flux.values[g, :, 0]   # (ng, nx, ny) → (nx,)
    return mesh, sn, psi, phi


def _require(condition: bool, message: str) -> None:
    """``-O``-survivable assertion (bare ``assert`` is stripped under -O)."""
    if not condition:
        raise AssertionError(message)


# ═══════════════════════════════════════════════════════════════════════
# T1 / T2 — slab non-vacuum convergence (1g + 2g asymmetric Σs)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.parametrize("case_kind", ["1g", "2g_asymmetric_S"])
@pytest.mark.verifies(
    "transport-cartesian",
    "dd-cartesian-1d",
    "dd-slab",
    "sn-mms-nonvacuum-psi",
    "sn-mms-nonvacuum-qext",
)
def test_mms_prescribed_inflow_slab_converges_second_order(case_kind: str):
    r"""Diamond-difference SN on a NON-VACUUM slab shows O(h²) AND
    converges to the imposed boundary-non-zero flux :math:`A(x)`.

    Three assertions per group (rate is necessary-not-sufficient):

    1. **Rate** ``orders > 1.9`` — DD design order on a smooth ansatz.
    2. **Converged VALUE** (the structural-independence anchor): the
       finest-mesh :math:`\phi_{\rm num}` matches :math:`A(x)` — which is
       NON-zero at the boundary (a0>0).  A dropped ``q.boundary`` (silent
       vacuum) converges O(h²) to the boundary-ZERO limit, so only this
       value check sees it.
    3. **Inflow-honoured spot-check**: the solved trace slot equals the
       imposed :math:`\gamma_-\psi = (A + \mu B)/W`.

    Measured (probe, nc∈{20,40,80}): 1g orders [2.04, 2.01], finest
    L2err ~1.2e-3, pointwise max|φ-A| ~8e-5; 2g g0 [2.04, 2.01],
    g1 [2.05, 2.01].  T2 additionally exercises the ERR-002 ``SigSᵀ``
    group-transfer hazard via the asymmetric downscatter-only Σs.
    """
    if case_kind == "1g":
        case = build_slab_nonvacuum_mms_case()
        n_groups = 1
    elif case_kind == "2g_asymmetric_S":
        case = build_slab_2g_nonvacuum_mms_case()
        n_groups = 2
    else:  # pragma: no cover - guarded by parametrize
        raise AssertionError(case_kind)

    n_cells = [20, 40, 80]
    W = float(case.quadrature.weights.sum())
    mu = case.quadrature.mu_x

    for g in range(n_groups):
        errors = []
        finest = None
        for nc in n_cells:
            mesh, sn, psi, phi = _solve(case, nc, g=g)
            ref = case.phi_exact(mesh.centers, g)
            errors.append(_l2_error(phi, ref, mesh.widths))
            finest = (mesh, sn, psi, phi, ref)
        errors = np.asarray(errors)
        orders = np.log2(errors[:-1] / errors[1:])

        # (1) Rate.
        _require(
            np.all(orders > 1.9),
            f"[{case_kind} g={g}] expected O(h²), got orders={orders} "
            f"from errors={errors}",
        )

        # (2) Converged VALUE (catches a dropped q.boundary — vacuum
        # converges O(h²) to the WRONG, boundary-zero limit).
        mesh_f, sn_f, psi_f, phi_f, ref_f = finest
        np.testing.assert_allclose(
            phi_f, ref_f, rtol=5e-3, atol=5e-3,
            err_msg=(
                f"[{case_kind} g={g}] non-vacuum φ converged to wrong "
                f"limit — prescribed inflow not honoured"
            ),
        )
        # The reference IS non-zero at the boundary (a0>0) — pin it so
        # the value check above is genuinely discriminating.
        _require(
            abs(float(ref_f[0])) > 0.1 and abs(float(ref_f[-1])) > 0.1,
            f"[{case_kind} g={g}] reference flux is ~vacuum at the faces "
            f"(ref[0]={ref_f[0]}, ref[-1]={ref_f[-1]}) — non-vacuum-ness lost",
        )

        # (3) Inflow-honoured spot-check on the finest mesh.
        for face, x_face in (("xmin", 0.0), ("xmax", case.slab_length)):
            inflow = sn_f.trace.inflow_indices_for_face(face)
            got = psi_f.boundary.face_view(face)[inflow, g]
            want = (case.A(x_face, g) + mu[inflow] * case.B(x_face, g)) / W
            np.testing.assert_allclose(
                got, want, rtol=1e-9, atol=1e-12,
                err_msg=f"[{case_kind} g={g}] γ₋ψ slot not honoured on {face}",
            )


# ═══════════════════════════════════════════════════════════════════════
# T3 — sphere non-vacuum (activates curvilinear redistribution); xfail #195
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.catches("ERR-026")
@pytest.mark.xfail(
    strict=True,
    reason=(
        "ERR-026 PARTIAL / #195 curvilinear DD pre-asymptotic — see "
        "tests/sn/verification/mms/test_mms_curvilinear.py. The non-vacuum "
        "+ redistribution machinery is verified (rate + value), but the "
        "underlying curvilinear DD spatial convergence it rides on is the "
        "OPEN #195 gate (the absolute-magnitude band fails until the "
        "pole-face spatial closure aligns with Hébert §3.9.4). Marker comes "
        "off when #195 closes."
    ),
)
@pytest.mark.verifies(
    "transport-spherical",
    "sn-mms-nonvacuum-sph-psi",
    "sn-mms-nonvacuum-sph-qext",
)
def test_mms_prescribed_inflow_sphere_activates_redistribution():
    r"""Curvilinear angular redistribution :math:`(1-\mu^2)/r\,\partial_\mu\psi
    = (1-\mu^2)B(r)/r` under NON-VACUUM inflow at r=R (Mode-7 mandatory
    companion to the slab rows).

    Activates: radial streaming, within-group scatter (c<1), :math:`\mu^2`
    second moment, the redistribution term (the ERR-026 path the slab
    nulls), and the prescribed non-vacuum inflow at the outer face
    (A(R)>0).  Pole r=0 is handled by the symmetry BC (not a face;
    B(0)=0 keeps the redistribution regular — HAZARD H1).

    Ships ``xfail(strict)`` on #195: the curvilinear DD spatial
    convergence is pre-asymptotic at these meshes (the absolute-magnitude
    band fails despite the Phase D Carlson seed restoring the O(h²) rate
    on the empirical probe).  This MUST collect + xfail strictly (NOT
    xpass — an xpass means #195 closed and the marker must come off).
    """
    case = build_sphere_nonvacuum_mms_case()
    n_cells = [20, 40, 80]
    W = float(case.quadrature.weights.sum())
    mu = case.quadrature.mu_x

    errors = []
    finest = None
    for nc in n_cells:
        mesh, sn, psi, phi = _solve(case, nc, g=0)
        ref = case.phi_exact(mesh.centers)
        # volume-weighted L2 (curvilinear — H-vol note).
        errors.append(_l2_error(phi, ref, mesh.volumes))
        finest = (mesh, sn, psi, phi, ref)
    errors = np.asarray(errors)
    orders = np.log2(errors[:-1] / errors[1:])

    # Rate gate (the #195 magnitude trap fires the xfail; if the rate
    # also regresses the xfail still strict-passes).
    _require(
        np.all(orders[-2:] > 1.9),
        f"sphere non-vacuum: expected O(h²), got orders={orders}",
    )

    # The #195 magnitude band — this is the gate that FAILS pre-asymptotic
    # at these meshes (mirrors test_mms_curvilinear.py). The xfail rides
    # on this assertion.
    mesh_f, sn_f, psi_f, phi_f, ref_f = finest
    _require(
        1e-8 < errors[-1] < 1e-3,
        f"sphere non-vacuum: finest L2 error {errors[-1]:.3e} outside the "
        f"#195 in-band [1e-8, 1e-3] (pre-asymptotic).",
    )
    np.testing.assert_allclose(
        phi_f, ref_f, rtol=2e-2, atol=2e-2,
        err_msg="sphere non-vacuum φ converged to wrong limit",
    )


# ═══════════════════════════════════════════════════════════════════════
# T3g — GREEN sphere structural check (non-convergence-dependent)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_sphere_nonvacuum_inflow_honoured_and_redistribution_live():
    r"""Sphere non-vacuum structural check that exercises the non-vacuum
    inflow + redistribution path NOW (NOT convergence-dependent — green
    today, complementing the xfail T3).

    Two claims:
    1. The prescribed inflow at r=R is honoured per inflow ordinate
       (:math:`\gamma_-\psi = (A(R) + \mu B(R))/W`, A(R)>0 non-vacuum).
    2. The redistribution source :math:`(1-\mu^2)B(r)/r` is non-zero on
       the manufactured source (the ERR-026 term is live under the 4.6
       ansatz — B(r)≠0 on the open interval, with B(0)=0 pole-regular).
    """
    case = build_sphere_nonvacuum_mms_case()
    mesh, sn, psi, phi = _solve(case, 40, g=0)
    W = float(case.quadrature.weights.sum())
    mu = case.quadrature.mu_x

    # (1) inflow honoured at r=R (the sphere's only boundary face).
    face = "xmax"
    inflow = sn.trace.inflow_indices_for_face(face)
    got = psi.boundary.face_view(face)[inflow, 0]
    R = case.radius
    want = (case.A(np.array([R]))[0] + mu[inflow] * case.B(np.array([R]))[0]) / W
    np.testing.assert_allclose(
        got, want, rtol=1e-9, atol=1e-12,
        err_msg="sphere γ₋ψ slot at r=R not honoured",
    )
    # ... and A(R) is genuinely non-vacuum.
    _require(
        abs(float(case.A(np.array([R]))[0])) > 0.1,
        "sphere reference A(R) ~ vacuum — non-vacuum-ness lost",
    )

    # (2) the redistribution source (1-μ²)B(r)/r is non-zero on the mesh
    # interior (B(0)=0 pole-regular; B(r)≠0 elsewhere → ERR-026 path live).
    r = mesh.centers
    B_r = case.B(r)                                # (nx,)
    redistribution = (1.0 - mu[:, None] ** 2) * (B_r / r)[None, :]  # (N, nx)
    _require(
        np.abs(redistribution).max() > 1e-6,
        "redistribution source (1-μ²)B/r is ~zero — ERR-026 path nulled",
    )
