r"""Foundation contract pins for the path-forward ``(L + C)`` matvec.

R-1 Step 4 Step G0 (2026-05-22) — the unified matvec switched to a
native typed-field ↔ typed-field signature; the module-level helper
``_transport_operator_matvec_unified`` it landed as was DELETED at the
Wave T T.5 matvec retirement.  The kernel now lives on the loss
representation (:meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk`,
the fused ``(L+C)ψ`` walk) and is reached through the public operator
algebra :meth:`~orpheus.sn.operators.streaming.StreamingCollisionOperator.apply`;
every gate below drives it through the ``_LC_matvec`` shim.  The
face-state I/O is direct ``(N, ng)``
:class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
access (no ``EquationMap.face_outer_ordinate`` slot map; the inflow /
outflow ordinate masks derive from ``quad.mu_x`` direction signs).

This file pins the structural contract the path-forward matvec must
satisfy.  The V&V plan §G0.2
(``.claude/plans/r1_step4_g_verification_plan.md``) listed seven
foundation pins; five survive as gates and two retired with the
machinery they cross-checked:

1. Zero input → zero output (linearity sentinel).
2. Uniform ψ on homogeneous reflective medium → ``(L+C)·ψ = σ_t·ψ``
   (the streaming-equilibrium / per-ordinate flat-flux invariant —
   the ERR-026 / ERR-006 / ERR-007 canonical diagnostic; Signature
   1 of ``numerical-bug-signatures``).
3. Linearity in ψ — ``M(αψ + βφ) = αM(ψ) + βM(φ)`` at FP-ULP.
4. Output shape ``(N, ng, *spatial)`` for cells + ``(N, ng)`` per
   boundary face (``boundary.face_view("xmax")``, and the slab inner
   face) — the path-forward typed contract.
5. Face residual, per face class.  The partition is **THREE-way**
   (inflow ⊔ outflow ⊔ tangential), so this pin has three clauses and
   they have independent fates:

   * **outflow** — the outflow-definition defect
     ``streamed − ψ.outflow``.
   * **inflow** — the ``I·ψ.inflow`` IDENTITY, *not* zero.  (Wave O
     #208 O.4a.2 inverted this clause: pre-carve the matvec re-applied
     the BC internally and emitted a zero residual there; post-carve
     the trace inflow is an explicit unknown driven by the boundary
     consistency residual ``ψ.inflow − B·ψ.outflow``.)
   * **tangential** (``μ_r = 0``) — an exact structural ZERO, and
     still true.  ⚠ It is NOT asserted in this file, and cannot be:
     all three fixtures here (GL-4 slab, GL-4 sphere, LS-4 cylinder)
     carry **zero** tangential ordinates — even-order Gauss-Legendre
     is the one production family with none.  A test added here would
     pass vacuously.  The tree's sole catcher is
     ``test_sweep_inverse_identity.py`` on its ``cyl_product``
     fixture, which now fails loudly if that fixture ever loses its
     tangential ordinates.

   Historical note: this pin used to read "residual zero at
   **non-outflow** ordinates", which conflated the last two clauses —
   "not outflow" is not "inflow".  When Wave O inverted the inflow
   half, a 2026-08-03 correction rewrote the sentence around it and
   dropped the tangential half entirely; measurement restored it.
6. *(RETIRED — D-J 2026-05-30.)*  ``TestQuadDerivedMaskEqualsLegacySlotMap``
   pinned ``quad.mu_x > 0`` ≡ ``eq_map.face_outer_ordinate``; the
   ``eq_map`` side no longer exists, so the equivalence is the
   definition rather than a gate.
7. *(INVERTED — Wave T post-T.5.)*  The 2-D Cartesian
   ``NotImplementedError`` guard is GONE; 2-D Cartesian ``(L+C).apply``
   drives the representation's ``loss_action`` and returns a valid
   result.  ``TestTwoDCartesianRaises`` now pins the positive contract.

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
from orpheus.sn.mesh.augmented_mesh import SNMesh
from tests.sn._test_helpers import _LC_matvec
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.source_sinks import AngularSourceSink, AngularBoundarySourceSink
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.radial_characteristic_field import RadialCharacteristicField
from orpheus.transport.full_field import FullField
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
    """1-D cylinder with reflective inner / vacuum outer + the folded rule."""
    mesh = Mesh1D(
        edges=np.linspace(0.01, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.folded_product(n_mu=4, n_phi=8)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


GEOMETRIES = [
    ("slab", _slab_mesh),
    ("sphere", _sphere_mesh),
    ("cylinder", _cylinder_mesh),
]


def _zero_flux(sn_mesh: SNMesh) -> TimedFullField:
    """Construct a zero :class:`TimedFullField` on ``sn_mesh``."""
    # #282 route (a): pass the seed leaf UNIFORMLY — the R12a predicate
    # allocates the present-but-ZERO ψ½ block iff the mesh carries levels
    # (sphere yes; slab/cyl stay None).  apply(0)=0 holds on every block.
    return TimedFullField.zeros(
        interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh,
    )


def _uniform_flux(sn_mesh: SNMesh, value: float = 1.0) -> TimedFullField:
    """Construct a uniform-ψ :class:`TimedFullField` with face state matching.

    The boundary face state is set to ``value`` on every face the geometry
    carries, preserving the pre-D-H.2-C4c semantic where bulk-uniform
    implies boundary-at-the-value (the flat-flux invariant input).
    """
    state = TimedFullField.zeros(
        interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh,
    )
    state.interior.values[:] = value
    for face in ("xmin", "xmax"):
        if face in state.boundary.layout.faces:
            state.boundary.face_view(face)[:] = value
    # #282 route (a) → B.2d: the CONSISTENT ψ½ seed for a per-ordinate flat
    # field is the same constant ``value`` (constant-preserving), riding the
    # walk's EXPLICIT flux leg on a carrying mesh — so the pole march
    # reproduces ``value`` at every cell and the flat-flux invariant
    # ``(L+C)·ψ = σ_t·ψ`` holds.  A zero seed leg would inject a wrong pole
    # datum and break the value assertion (rule 2).
    seed = None
    if sn_mesh.radial_characteristic_field_space is not None:
        seed = RadialCharacteristicField.from_mesh(sn_mesh)
        seed.interior.values[:] = value
        seed.boundary.values[:] = value
    return state, seed


# ── Pin 1: zero input → zero output ─────────────────────────────────


class TestZeroInputZeroOutput:
    """The matvec is linear; zero ψ + zero boundary → zero output."""

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_zero_input_zero_output(self, name, builder) -> None:
        sn_mesh = builder()
        sigma_t = np.full((sn_mesh.ng, *sn_mesh.spatial_shape), 2.0)
        result = _LC_matvec(
            _zero_flux(sn_mesh), sigma_t,
        )
        np.testing.assert_array_equal(
            result.interior.values, np.zeros_like(result.interior.values),
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
        sigma_t = np.full((ng, *sn_mesh.spatial_shape), sigma_t_val)
        state, seed = _uniform_flux(sn_mesh, value=1.0)
        result = _LC_matvec(
            state, sigma_t, radial_characteristic_flux=seed,
        )
        # Per-ordinate cell action: (L+C)·1 = σ_t·1 = 2.0.  Flat-flux
        # invariant holds for every ordinate, every cell.
        np.testing.assert_allclose(
            result.interior.values, sigma_t_val, rtol=1e-12, atol=1e-13,
            err_msg=(
                f"{name}: uniform-ψ flat-flux invariant violated; max "
                f"deviation = {np.abs(result.interior.values - sigma_t_val).max():.3e}"
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
    quad = Quadrature.folded_product(n_mu=4, n_phi=8)
    return SNMesh(mesh, quad, placeholder_materials())


# ── Pin 3: linearity ────────────────────────────────────────────────


class TestLinearity:
    """``M(αψ + βφ) = αM(ψ) + βM(φ)`` at FP-ULP."""

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_linearity(self, name, builder) -> None:
        sn_mesh = builder()
        ng = sn_mesh.ng
        N = sn_mesh.quad.N
        sigma_t = np.full((ng, *sn_mesh.spatial_shape), 0.5)

        rng = np.random.default_rng(seed=42)

        def _random_state() -> TimedFullField:
            state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh)
            state.interior.values[:] = rng.standard_normal((N, ng, *sn_mesh.spatial_shape))
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
        sum_psi = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh)
        sum_psi.interior.values[:] = alpha * psi.interior.values + beta * phi.interior.values
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
            m_sum.interior.values,
            alpha * m_psi.interior.values + beta * m_phi.interior.values,
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
        sigma_t = np.full((sn_mesh.ng, *sn_mesh.spatial_shape), 1.0)
        result = _LC_matvec(
            _zero_flux(sn_mesh), sigma_t,
        )
        # Composite carrier; the (L+C).apply output bulk AND boundary are
        # the source/sink role leaves (AngularSourceSink / AngularBoundarySourceSink)
        # — the operator output is Aψ (a rate-density / flux source), NOT a
        # residual.  B.5.2: a residual only arises from from_balance(Aψ, b),
        # never straight off an operator output (the boundary mirrors the bulk).
        # #257 S8a — the matvec leaf is a base arrow ``FullField -> FullField``
        # (the comonad lives on the driver); the output is the TIMELESS FullField.
        assert isinstance(result, FullField)
        assert not isinstance(result, TimedFullField)
        assert isinstance(result.interior, AngularSourceSink)
        assert isinstance(result.boundary, AngularBoundarySourceSink)
        # Cell values: (N, ng, *spatial).
        assert result.interior.values.shape == (
            sn_mesh.quad.N, sn_mesh.ng, *sn_mesh.spatial_shape,
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
    r"""The matvec writes the outflow-definition defect ``streamed −
    ψ.outflow`` at the outflow ordinates, and the IDENTITY ``ψ.inflow`` at
    the inflow ordinates.

    Wave O #208 O.4a.2 (BC extraction): pre-carve the inflow ordinates got
    a zero residual ("no equation — value comes from the BC"); the keystone
    re-applied the BC inside the matvec.  Post-carve ``L_full`` reads
    ``ψ.boundary.inflow`` as a GIVEN and EMITS the ``I·ψ.inflow`` identity
    on the inflow slots — the trace inflow became an explicit unknown,
    driven by the boundary consistency residual ``ψ.inflow − B·ψ.outflow``
    (the ``−B`` sibling supplies the ``−B·ψ.outflow`` term in
    ``(L+C−S−F−B)``).  So the inflow-ordinate output now EQUALS the input
    ``ψ.inflow`` (identity), not zero.

    Outflow at the outer face: ``quad.mu_x > 0``.  Inflow: ``mu_x < 0``.
    """

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_outer_face_inflow_slots_carry_the_identity(
        self, name, builder,
    ) -> None:
        sn_mesh = builder()
        ng = sn_mesh.ng
        sigma_t = np.full((ng, *sn_mesh.spatial_shape), 1.0)

        # Random ψ — the inflow-ordinate output is the I·ψ.inflow identity
        # row (the consistency-residual diagonal), so it tracks the input.
        rng = np.random.default_rng(seed=11)
        psi = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh)
        psi.interior.values[:] = rng.standard_normal(
            (sn_mesh.quad.N, ng, *sn_mesh.spatial_shape),
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
            psi.boundary.face_view("xmax")[inflow_outer, :],
            err_msg=(
                f"{name}: outer-face inflow-ordinate output (μ_x < 0) is "
                f"NOT the input ψ.inflow identity.  Post-extraction L_full "
                f"emits the I·ψ.inflow consistency-row diagonal there; the "
                f"reflective coupling is the sibling −B, not an intra-matvec "
                f"re-apply."
            ),
        )

    def test_slab_inner_face_inflow_slots_carry_the_identity(
        self,
    ) -> None:
        """Slab xmin face: outflow at μ < 0; inflow at μ > 0 (identity row)."""
        sn_mesh = _make_reflective_slab()
        ng = sn_mesh.ng
        sigma_t = np.full((ng, *sn_mesh.spatial_shape), 1.0)

        rng = np.random.default_rng(seed=22)
        psi = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh)
        psi.interior.values[:] = rng.standard_normal(
            (sn_mesh.quad.N, ng, *sn_mesh.spatial_shape),
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
            psi.boundary.face_view("xmin")[inflow_inner, :],
            err_msg=(
                "slab inner-face inflow-ordinate output (μ_x > 0) is NOT "
                "the input ψ.inflow identity.  Post-extraction the inflow "
                "slot carries the I·ψ.inflow consistency-row diagonal."
            ),
        )


# ── Pin 7: 2-D Cartesian returns a result (the guard is GONE) ───────
#
# (Pin 6 — ``TestQuadDerivedMaskEqualsLegacySlotMap`` — retired with
# D-J 2026-05-30.  It pinned the equivalence ``quad.mu_x > 0`` ≡
# ``eq_map.face_outer_ordinate``.  The eq_map side of that equation no
# longer exists; the production code now uses the quad-derived mask
# directly.  The equivalence is no longer a gate, it is the definition.)


class TestTwoDCartesianRaises:
    r"""2-D Cartesian coverage of the ``(L+C)`` matvec.

    The class NAME is historical: pre-Step-G the unified matvec handled
    only 1-D slab + curvilinear and raised ``NotImplementedError`` on
    2-D Cartesian, and this class pinned that guard.  Wave T post-T.5
    removed the guard — 2-D Cartesian drives the representation's
    ``loss_action`` and returns a result — so the pin INVERTED into the
    positive contract asserted below.  Failure here would still mean a
    silent 2-D regression.
    """

    def test_two_d_cartesian_loss_action_returns_result(self) -> None:
        """Wave T post-T.5 matvec retirement: the legacy 2-D guard
        in ``_transport_operator_matvec_unified`` (which raised
        ``NotImplementedError``) is GONE; 2-D Cartesian (L+C).apply now
        drives the representation's ``loss_action`` (the matvec walk that
        since S6.3 lives on the loss representation, off the operator;
        ``ScanMarch`` default since S6.9), which successfully computes
        the matvec.  The new
        contract: 2-D Cartesian (L+C).apply returns a valid
        TimedFullField result (no exception)."""
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
        # A genuine-2-D quadrature is required: the y-faces need ordinates
        # with non-zero mu_y. The 1-D gauss_legendre set has mu_y=0 for
        # every ordinate, which the trace-space guard correctly rejects.
        quad = Quadrature.level_symmetric(sn_order=4)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())
        sigma_t = np.full((sn_mesh.ng, *sn_mesh.spatial_shape), 1.0)
        psi = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh)
        result = _LC_matvec(psi, sigma_t)
        # #257 S8a — base arrow output is the TIMELESS FullField.
        assert isinstance(result, FullField)
        assert not isinstance(result, TimedFullField)
        # On zero input the matvec is zero (linearity sentinel).
        np.testing.assert_array_equal(
            result.interior.values, np.zeros_like(result.interior.values),
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
        from orpheus.sn.operators.streaming import (
            StreamingOperator,
        )
        from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
        sn_mesh = _slab_mesh()
        sigma_t = np.full((sn_mesh.ng, *sn_mesh.spatial_shape), 1.0)
        psi_bare = np.zeros((sn_mesh.quad.N, sn_mesh.ng, *sn_mesh.spatial_shape))
        L_op = StreamingOperator(sn_mesh)
        C_op = MultiplicationOperator.from_mesh(sigma_t, sn_mesh)
        with pytest.raises(TypeError, match="TimedFullField"):
            (L_op + C_op).apply(psi_bare)
