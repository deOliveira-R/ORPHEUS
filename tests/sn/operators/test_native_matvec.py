r"""Foundation contract pins for the path-forward
:func:`~orpheus.sn.operator._transport_operator_matvec_unified`.

R-1 Step 4 Step G0 (2026-05-22) — the unified matvec switched to a
native :class:`AngularFlux` ↔ :class:`AngularFlux` signature.  The
face-state I/O is direct ``(N, ng)`` :class:`BoundaryFlux` access (no
``EquationMap.face_outer_ordinate`` slot map; the inflow / outflow
ordinate masks derive from ``quad.mu_x`` direction signs).

This file pins the structural contract the path-forward matvec must
satisfy.  Per V&V plan §G0.2 (``.claude/plans/r1_step4_g_verification_plan.md``)
the seven foundation pins are:

1. Zero input → zero output (linearity sentinel).
2. Uniform ψ on homogeneous reflective medium → ``(L+C)·ψ = σ_t·ψ``
   (the streaming-equilibrium / per-ordinate flat-flux invariant —
   the ERR-026 / ERR-006 / ERR-007 canonical diagnostic; Signature
   1 of ``numerical-bug-signatures``).
3. Linearity in ψ — ``M(αψ + βφ) = αM(ψ) + βM(φ)`` at FP-ULP.
4. Output shape ``(N, ng, nx, ny)`` for cells + ``(N, ng)`` for
   ``boundary.xmax_face`` (and slab-inner) — the path-forward
   typed contract.
5. Face residual zero at non-outflow ordinates (no equation there;
   inflow ords get their value from the BC, not from the iterate).
6. Quadrature-derived outflow masks at the outer face match the
   legacy ``face_outer_ordinate`` slot positions exactly.  Pin via
   ``build_equation_map_with_traces`` cross-check.
7. 2-D Cartesian raises ``NotImplementedError`` (mirrors the
   existing 2-D guard pre-Step-G; Phase A absorbs).

Cross-references:

* ``.claude/plans/r1_step4_g_convention_crosswalk.md`` Axis 3 + 6.
* ``.claude/plans/r1_step4_g_dependency_audit.md`` CORRECTION
  section (Path forward).
* ``.claude/lessons.md`` L18 (Pattern 7 producer-side normalisation;
  the face-mask derivation from ``quad`` is the equivalent of
  Pattern 7 applied to the slot-map elimination).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn.geometry import SNMesh
from tests.sn._test_helpers import _LC_matvec
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.source_sinks import AngularSourceSink
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField
from orpheus.numerics.quadrature import Quadrature
from tests.sn._test_helpers import placeholder_materials


pytestmark = pytest.mark.foundation


# ── Mesh fixtures ────────────────────────────────────────────────────


def _slab_mesh(nx: int = 4, length: float = 1.0, ng: int = 1) -> SNMesh:
    """1-D slab with vacuum BCs + GL N=4."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, length, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _sphere_mesh(nx: int = 4, radius: float = 1.0, ng: int = 1) -> SNMesh:
    """1-D sphere with reflective inner / vacuum outer + GL N=4."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _cylinder_mesh(nx: int = 4, radius: float = 1.0, ng: int = 1) -> SNMesh:
    """1-D cylinder with reflective inner / vacuum outer + LS-4."""
    mesh = Mesh1D(
        edges=np.linspace(0.01, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


GEOMETRIES = [
    ("slab", _slab_mesh),
    ("sphere", _sphere_mesh),
    ("cylinder", _cylinder_mesh),
]


def _zero_flux(sn_mesh: SNMesh) -> TimedFullField:
    """Construct a zero :class:`TimedFullField` on ``sn_mesh``."""
    return TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh)


def _uniform_flux(sn_mesh: SNMesh, value: float = 1.0) -> TimedFullField:
    """Construct a uniform-ψ :class:`TimedFullField` with face state matching.

    The boundary face state is set to ``value`` on every face the geometry
    carries, preserving the pre-D-H.2-C4c semantic where bulk-uniform
    implies boundary-at-the-value (the flat-flux invariant input).
    """
    state = TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh)
    state.bulk.values[:] = value
    for face in ("xmin", "xmax"):
        if face in state.boundary.layout.faces:
            state.boundary.face_view(face)[:] = value
    return state


# ── Pin 1: zero input → zero output ─────────────────────────────────


class TestZeroInputZeroOutput:
    """The matvec is linear; zero ψ + zero boundary → zero output."""

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_zero_input_zero_output(self, name, builder) -> None:
        sn_mesh = builder()
        sigma_t = np.full((sn_mesh.ng, sn_mesh.nx, sn_mesh.ny), 2.0)
        result = _LC_matvec(
            _zero_flux(sn_mesh), sigma_t,
        )
        np.testing.assert_array_equal(
            result.bulk.values, np.zeros_like(result.bulk.values),
        )
        # Face residual at outflow positions = (WDD-propagated 0) - (stored 0)
        # = 0; at inflow positions = 0 by default.  Whole array zero.
        xmax_view = result.boundary.face_view("xmax")
        np.testing.assert_array_equal(xmax_view, np.zeros_like(xmax_view))
        if "xmin" in result.boundary.layout.faces:
            xmin_view = result.boundary.face_view("xmin")
            np.testing.assert_array_equal(xmin_view, np.zeros_like(xmin_view))


# ── Pin 2: uniform ψ on homogeneous reflective → σ_t·ψ ──────────────


class TestUniformFluxSigmaT:
    r"""On a homogeneous reflective medium with uniform face state, the
    streaming term cancels exactly (per-ordinate flat-flux invariant)
    and the matvec returns ``(L+C)·ψ = σ_t·ψ`` per ordinate.

    This is the CANONICAL ERR-026 / ERR-006 / ERR-007 diagnostic — the
    Pomraning isotropy test extended to all geometries.  Pre-Phase G
    Step 2 the curvilinear SI sweep failed this with ~22% L0 error
    (ERR-048 manifestation).  The path-forward matvec inherits the
    geometry-agnostic body and the M-M Carlson seed structure that
    closed ERR-048.
    """

    @pytest.mark.verifies("transport-cartesian", "hebert-3-432")
    @pytest.mark.parametrize("name,builder", [
        ("slab_reflective", lambda: _make_reflective_slab()),
        ("sphere_reflective", lambda: _make_reflective_sphere()),
        ("cylinder_reflective", lambda: _make_reflective_cylinder()),
    ])
    def test_uniform_flux_on_homogeneous_reflective_gives_sigma_t(
        self, name, builder,
    ) -> None:
        sn_mesh = builder()
        ng = sn_mesh.ng
        sigma_t_val = 2.0
        sigma_t = np.full((ng, sn_mesh.nx, sn_mesh.ny), sigma_t_val)
        result = _LC_matvec(
            _uniform_flux(sn_mesh, value=1.0), sigma_t,
        )
        # Per-ordinate cell action: (L+C)·1 = σ_t·1 = 2.0.  Flat-flux
        # invariant holds for every ordinate, every cell.
        np.testing.assert_allclose(
            result.bulk.values, sigma_t_val, rtol=1e-12, atol=1e-13,
            err_msg=(
                f"{name}: uniform-ψ flat-flux invariant violated; max "
                f"deviation = {np.abs(result.bulk.values - sigma_t_val).max():.3e}"
            ),
        )


def _make_reflective_slab(nx: int = 4, length: float = 1.0) -> SNMesh:
    """Reflective slab (homogeneous infinite slab) — flat-flux invariant
    needs reflective BC so the outflow at the boundary reflects back as
    inflow, nulling the net streaming on uniform ψ."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, length, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials())


def _make_reflective_sphere(nx: int = 4, radius: float = 1.0) -> SNMesh:
    """Reflective outer sphere — uniform ψ flat-flux invariant requires
    the outer BC to reflect the outflow back as inflow (so the net
    radial streaming on uniform ψ cancels).  Pole at r=0 always
    handles regularity via M-M Carlson seed."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),  # unused — pole regularity
        bc_right=BC("reflective"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials())


def _make_reflective_cylinder(nx: int = 4, radius: float = 1.0) -> SNMesh:
    """Reflective outer cylinder — same flat-flux logic as sphere."""
    mesh = Mesh1D(
        edges=np.linspace(0.01, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials())


# ── Pin 3: linearity ────────────────────────────────────────────────


class TestLinearity:
    """``M(αψ + βφ) = αM(ψ) + βM(φ)`` at FP-ULP."""

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_linearity(self, name, builder) -> None:
        sn_mesh = builder()
        ng = sn_mesh.ng
        N, nx, ny = sn_mesh.quad.N, sn_mesh.nx, sn_mesh.ny
        sigma_t = np.full((ng, nx, ny), 0.5)

        rng = np.random.default_rng(seed=42)

        def _random_state() -> TimedFullField:
            state = TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh)
            state.bulk.values[:] = rng.standard_normal((N, ng, nx, ny))
            state.boundary.face_view("xmax")[:] = rng.standard_normal((N, ng))
            if "xmin" in state.boundary.layout.faces:
                state.boundary.face_view("xmin")[:] = rng.standard_normal((N, ng))
            return state

        # ψ + φ with random face state too — linearity must hold
        # across cell + boundary slots.
        psi = _random_state()
        phi = _random_state()

        alpha, beta = 1.7, -0.3

        # M(αψ + βφ)
        sum_psi = TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh)
        sum_psi.bulk.values[:] = alpha * psi.bulk.values + beta * phi.bulk.values
        sum_psi.boundary.face_view("xmax")[:] = (
            alpha * psi.boundary.face_view("xmax")
            + beta * phi.boundary.face_view("xmax")
        )
        if "xmin" in sum_psi.boundary.layout.faces:
            sum_psi.boundary.face_view("xmin")[:] = (
                alpha * psi.boundary.face_view("xmin")
                + beta * phi.boundary.face_view("xmin")
            )
        m_sum = _LC_matvec(sum_psi, sigma_t)

        # αM(ψ) + βM(φ)
        m_psi = _LC_matvec(psi, sigma_t)
        m_phi = _LC_matvec(phi, sigma_t)

        np.testing.assert_allclose(
            m_sum.bulk.values,
            alpha * m_psi.bulk.values + beta * m_phi.bulk.values,
            rtol=1e-12, atol=1e-13,
            err_msg=f"{name}: linearity violated on cell slot",
        )
        np.testing.assert_allclose(
            m_sum.boundary.face_view("xmax"),
            alpha * m_psi.boundary.face_view("xmax")
            + beta * m_phi.boundary.face_view("xmax"),
            rtol=1e-12, atol=1e-13,
            err_msg=f"{name}: linearity violated on xmax_face",
        )


# ── Pin 4: output shape contract ────────────────────────────────────


class TestOutputShape:
    """The matvec returns ``TimedFullField`` with the path-forward shapes."""

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_output_shape_matches_input(self, name, builder) -> None:
        sn_mesh = builder()
        sigma_t = np.full((sn_mesh.ng, sn_mesh.nx, sn_mesh.ny), 1.0)
        result = _LC_matvec(
            _zero_flux(sn_mesh), sigma_t,
        )
        # Composite carrier; the (L+C).apply output bulk is an
        # AngularSourceSink (rate-density source/sink, not a flux — B.5.2).
        assert isinstance(result, TimedFullField)
        assert isinstance(result.bulk, AngularSourceSink)
        assert isinstance(result.boundary, BoundaryFlux)
        # Cell values: (N, ng, nx, ny).
        assert result.bulk.values.shape == (
            sn_mesh.quad.N, sn_mesh.ng, sn_mesh.nx, sn_mesh.ny,
        )
        # Outer face: (N, ng) for every geometry.
        assert result.boundary.face_view("xmax").shape == (
            sn_mesh.quad.N, sn_mesh.ng,
        )
        # Inner face: (N, ng) for slab; absent for curvilinear.
        curv = getattr(sn_mesh, "curvature", None)
        if curv is None:
            assert "xmin" in result.boundary.layout.faces
            assert result.boundary.face_view("xmin").shape == (
                sn_mesh.quad.N, sn_mesh.ng,
            )
        else:
            assert "xmin" not in result.boundary.layout.faces


# ── Pin 5: face residual zero at non-outflow ────────────────────────


class TestFaceResidualMask:
    r"""The matvec writes the face residual ONLY at outflow ordinates;
    inflow ordinates at the same face get a zero residual (no equation
    — their value comes from the BC, not from the iterate).

    Outflow at the outer face: ``quad.mu_x > 0``.  Inflow: ``mu_x < 0``.
    """

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_outer_face_residual_zero_at_inflow_ordinates(
        self, name, builder,
    ) -> None:
        sn_mesh = builder()
        ng = sn_mesh.ng
        sigma_t = np.full((ng, sn_mesh.nx, sn_mesh.ny), 1.0)

        # Random ψ — face residual at inflow ords must stay zero
        # regardless of input data (no equation there).
        rng = np.random.default_rng(seed=11)
        psi = TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh)
        psi.bulk.values[:] = rng.standard_normal(
            (sn_mesh.quad.N, ng, sn_mesh.nx, sn_mesh.ny),
        )
        psi.boundary.face_view("xmax")[:] = rng.standard_normal(
            (sn_mesh.quad.N, ng),
        )
        if "xmin" in psi.boundary.layout.faces:
            psi.boundary.face_view("xmin")[:] = rng.standard_normal(
                (sn_mesh.quad.N, ng),
            )

        result = _LC_matvec(psi, sigma_t)

        mu_x = sn_mesh.quad.mu_x
        eps = 1e-15
        inflow_outer = mu_x <= -eps  # μ_x < 0 = inflow at outer face
        np.testing.assert_array_equal(
            result.boundary.face_view("xmax")[inflow_outer, :],
            np.zeros((np.sum(inflow_outer), ng)),
            err_msg=(
                f"{name}: outer-face residual is non-zero at inflow "
                f"ordinates (μ_x < 0).  These have no equation; the "
                f"residual at these positions MUST stay zero."
            ),
        )

    def test_slab_inner_face_residual_zero_at_inflow_ordinates(
        self,
    ) -> None:
        """Slab xmin face: outflow at μ < 0; inflow at μ > 0."""
        sn_mesh = _make_reflective_slab()
        ng = sn_mesh.ng
        sigma_t = np.full((ng, sn_mesh.nx, sn_mesh.ny), 1.0)

        rng = np.random.default_rng(seed=22)
        psi = TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh)
        psi.bulk.values[:] = rng.standard_normal(
            (sn_mesh.quad.N, ng, sn_mesh.nx, sn_mesh.ny),
        )
        psi.boundary.face_view("xmax")[:] = rng.standard_normal(
            (sn_mesh.quad.N, ng),
        )
        psi.boundary.face_view("xmin")[:] = rng.standard_normal(
            (sn_mesh.quad.N, ng),
        )

        result = _LC_matvec(psi, sigma_t)

        mu_x = sn_mesh.quad.mu_x
        eps = 1e-15
        # Inflow at xmin face: μ_x > 0 (right-going).
        inflow_inner = mu_x >= +eps
        np.testing.assert_array_equal(
            result.boundary.face_view("xmin")[inflow_inner, :],
            np.zeros((np.sum(inflow_inner), ng)),
            err_msg=(
                "slab inner-face residual is non-zero at inflow "
                "ordinates (μ_x > 0).  Inflow at xmin has no equation; "
                "residual MUST stay zero."
            ),
        )


# ── Pin 7: 2-D Cartesian raises NotImplementedError ─────────────────
#
# (Pin 6 — ``TestQuadDerivedMaskEqualsLegacySlotMap`` — retired with
# D-J 2026-05-30.  It pinned the equivalence ``quad.mu_x > 0`` ≡
# ``eq_map.face_outer_ordinate``.  The eq_map side of that equation no
# longer exists; the production code now uses the quad-derived mask
# directly.  The equivalence is no longer a gate, it is the definition.)


class TestTwoDCartesianRaises:
    r"""Mirrors the existing 2-D guard pre-Step-G — the unified matvec
    only handles 1-D slab + curvilinear; 2-D Cartesian absorbs in
    Phase A.  Failure here would mean a silent 2-D regression.
    """

    def test_two_d_cartesian_routes_through_apply_2d_cartesian(self) -> None:
        """Wave T post-T.5 matvec retirement: the legacy 2-D guard
        in ``_transport_operator_matvec_unified`` (which raised
        ``NotImplementedError``) is GONE; 2-D Cartesian routes through
        :meth:`StreamingOperator._apply_2d_cartesian` instead, which
        successfully computes the matvec.  The new contract: 2-D
        Cartesian (L+C).apply returns a valid TimedFullField result
        (no exception)."""
        from orpheus.geometry.mesh import Mesh2D
        from orpheus.transport.timed_full_field import TimedFullField
        # Need ny > 1 for the 2-D path to fire.
        mesh = Mesh2D(
            edges_x=np.array([0.0, 1.0, 2.0]),
            edges_y=np.array([0.0, 1.0, 2.0]),
            mat_map=np.zeros((2, 2), dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
            bc_ymin=BC("vacuum"), bc_ymax=BC("vacuum"),
        )
        quad = Quadrature.gauss_legendre(n_ordinates=4)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())
        sigma_t = np.full((sn_mesh.ng, sn_mesh.nx, sn_mesh.ny), 1.0)
        psi = TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh)
        result = _LC_matvec(psi, sigma_t)
        assert isinstance(result, TimedFullField)
        # On zero input the matvec is zero (linearity sentinel).
        np.testing.assert_array_equal(
            result.bulk.values, np.zeros_like(result.bulk.values),
        )


# ── Sentinel: rejects non-AngularFlux input ─────────────────────────


class TestTypeContract:
    """Path-forward signature accepts :class:`TimedFullField` only;
    bare ndarrays (and the bulk-only L2 AngularFlux) are rejected with
    TypeError.  D-H.2-C4c flipped the kernel signature from the
    bulk-only AngularFlux to the L2 composite carrier."""

    def test_rejects_bare_ndarray(self) -> None:
        """The (L+C).apply operator-algebra path rejects bare ndarrays
        via `StreamingOperator.apply`'s TimedFullField type check
        (post-D-I.3d).  Wave T post-T.5 (matvec retirement) routes
        the check through the public operator-algebra surface."""
        from orpheus.sn.operator import (
            CollisionOperator,
            StreamingOperator,
        )
        sn_mesh = _slab_mesh()
        sigma_t = np.full((sn_mesh.ng, sn_mesh.nx, sn_mesh.ny), 1.0)
        psi_bare = np.zeros((sn_mesh.quad.N, sn_mesh.ng, sn_mesh.nx, sn_mesh.ny))
        L_op = StreamingOperator(sn_mesh, sigma_t)
        C_op = CollisionOperator(sn_mesh, sigma_t)
        with pytest.raises(TypeError, match="TimedFullField"):
            (L_op + C_op).apply(psi_bare)
