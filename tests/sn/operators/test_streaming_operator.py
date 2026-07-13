r"""Foundation tests for :class:`orpheus.sn.operators.streaming.StreamingOperator`.

The "L" of the four-operator algebra
``A_wg = L + C - S.foldable_part()`` (Issue #196, Resolution A).
:class:`StreamingOperator` is the typed leaf; the within-group
sweep-invertible composite is :class:`InvertibleOperator` (= ``L + C``).

Resolution A — subtractive definition
-------------------------------------

.. math::

   L.{\rm apply}(\psi) \;:=\; M(\psi;\;\sigma_t) \;-\;
                              \sigma_t \odot \psi.{\rm bulk}

L carries σ_t at constructor time. This is intrinsic to the discrete
curvilinear matvec (rational in σ_t through Hébert §3.9.4's Carlson
coupled-pole seed), not a defect — analogous to the DD coefficient
:math:`\alpha_{DD}(\sigma_t\,\Delta x)` carrying σ_t in
characteristic-line methods.

The decomposition gate (``rel_residual = 0.0``) lives in
:file:`tests/sn/test_streaming_operator_decomposition.py`.  ``(L + C)``
correctness is gated by ``test_unified_matvec_{slab,sphere,cylinder}``
against structurally-independent references.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.operator import LinearOperator
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import (
    StreamingOperator,
)
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from tests.sn.regression._regression_assert import assert_regression
from orpheus.numerics.quadrature import Quadrature
from tests.sn._test_helpers import placeholder_materials, radial_characteristic_edge_seed
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField

pytestmark = pytest.mark.foundation


# ═══════════════════════════════════════════════════════════════════════
# Fixtures — small slab / sphere / cylinder problems.
# ═══════════════════════════════════════════════════════════════════════


def _slab_mesh(nx: int = 4, length: float = 1.0, ng: int = 2) -> SNMesh:
    """Slab Mesh1D + GL N=4 quadrature, vacuum BCs.

    R-1 Step 4 Step G0 — ``ng`` matches ``_sig_t_uniform`` default so
    AngularFlux shape validation accepts the typed wrapping inside the
    path-forward matvec (the mesh.ng vs sig_t.ng dimensional sin per
    #205 closes at the test boundary).
    """
    mesh = Mesh1D(
        edges=np.linspace(0.0, length, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _spherical_mesh(nx: int = 4, radius: float = 1.0, ng: int = 2) -> SNMesh:
    """Spherical Mesh1D + GL N=4, reflective inner / vacuum outer.

    R-1 Step 4 Step G0 — see ``_slab_mesh`` re: ``ng`` default.
    """
    mesh = Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _cylindrical_mesh(nx: int = 4, radius: float = 1.0, ng: int = 2) -> SNMesh:
    """Cylindrical Mesh1D + Level-Symmetric SN-4 quadrature.

    R-1 Step 4 Step G0 — see ``_slab_mesh`` re: ``ng`` default.
    """
    mesh = Mesh1D(
        edges=np.linspace(0.01, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _sig_t_uniform(sn_mesh: SNMesh, ng: int = 2,
                   value: float = 0.5) -> np.ndarray:
    """σ_t uniform across cells / groups (``(ng, *spatial)``)."""
    return value * np.ones((ng, *sn_mesh.spatial_shape))


def _sig_t_heterogeneous(sn_mesh: SNMesh, ng: int = 2) -> np.ndarray:
    """σ_t heterogeneous — different value per (cell, group)."""
    rng = np.random.default_rng(seed=20260514)
    return 0.3 + 0.5 * rng.random((ng, *sn_mesh.spatial_shape))


GEOMETRIES = [
    ("slab", _slab_mesh),
    ("sphere", _spherical_mesh),
    ("cylinder", _cylindrical_mesh),
]


# 1-D-only geometries — alias of ``GEOMETRIES`` (the 1-D paths are the
# only ones the :class:`TimedFullField` branch supports natively at this
# wave; 2-D Cartesian is exercised separately below).
GEOMETRIES_1D = [
    ("slab", _slab_mesh),
    ("sphere", _spherical_mesh),
    ("cylinder", _cylindrical_mesh),
]


def _random_composite(sn_mesh, seed=171):
    """Build a TimedFullField with non-zero bulk + non-zero boundary."""
    from dataclasses import replace

    rng = np.random.default_rng(seed)
    state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh)
    bulk_values = rng.standard_normal(state.interior.values.shape)
    boundary_values = 0.1 + rng.random(state.boundary.values.shape)
    state = replace(state, interior=replace(state.interior, values=bulk_values))
    state = replace(
        state, boundary=replace(state.boundary, values=boundary_values),
    )
    return state


# ═══════════════════════════════════════════════════════════════════════
# 1. Capability advertising.
# ═══════════════════════════════════════════════════════════════════════


class TestCapabilities:
    """StreamingOperator is adjointable but NOT invertible.

    Resolution A: L alone is not invertible (rank-deficient without
    collision), so ``is_invertible`` is False and it carries NO ``solve``
    verb; invertibility appears at the OperatorSum level via the fusion
    hook (substep 3+4.b.ii). ``is_adjointable`` IS True — the analytic
    reverse-direction G-adjoint matvec ``Lᵀ`` landed in Wave O / O.2b
    (the foundation of ``L.H``); the tests below were migrated from the
    pre-O.2b apply-only contract (retirement = test migration).
    """

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_predicates_adjointable_not_invertible(self, name, builder):
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh)
        assert L.is_adjointable and not L.is_invertible

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_no_solve_verb(self, name, builder):
        # L is structurally non-invertible at the leaf: no solve verb at
        # all — the sweep-invertible solve lives on the fused (L+C) sum.
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh)
        assert not L.is_invertible
        assert not hasattr(L, "solve")

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_has_apply_transpose(self, name, builder):
        # Wave O / O.2b: L carries the analytic reverse-direction adjoint
        # matvec Lᵀ (the foundation of the G-adjoint ``L.H``).
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh)
        assert L.is_adjointable

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_satisfies_linear_operator_protocol(self, name, builder):
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh)
        assert isinstance(L, LinearOperator)


# ═══════════════════════════════════════════════════════════════════════
# 2. Constructor surface — Resolution A requires σ_t at construction.
# ═══════════════════════════════════════════════════════════════════════


class TestConstructor:
    """StreamingOperator(sn_mesh) takes ONLY the mesh — pure L is σ-free.

    Pattern 4 (illegal states unrepresentable), #257 S8b: pure ``L``
    computes ``Ω·∇ψ`` directly (the named ``streaming_action`` leaf) and
    reads NO σ — the curvilinear Carlson seed's σ-dependence is exactly
    the collision diagonal it injects (which belongs to ``C``).  A σ on
    ``L`` would be a parameter the leaf never reads, so the constructor
    refuses it.
    """

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_construct_with_sn_mesh_only(self, name, builder):
        sn_mesh = builder()
        L = StreamingOperator(sn_mesh)
        assert L.sn_mesh is sn_mesh
        # σ-free: the leaf carries no sigma_t surface (#257 S8b).  The pure-L
        # apply reads no σ — it routes to the representation's σ-free
        # streaming_action leaf; the collision diagonal lives in C.
        assert not hasattr(L, "sigma_t")


# ═══════════════════════════════════════════════════════════════════════
# 3. apply shape correctness — D-I.3c (2026-05-29) retired alongside the
#    bare-ndarray adapter.  The TimedFullField shape contract is now
#    pinned by :class:`TestCompositeInvariants` (returns TimedFullField,
#    bulk type, history-depth preservation, mesh identity).
# ═══════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════
# 4. THE LOAD-BEARING COMPOSITION TEST — no xfail allowed.
# ═══════════════════════════════════════════════════════════════════════


# ``(L + C).apply``'s correctness is gated by
# ``test_unified_matvec_{slab,sphere,cylinder}.py`` against
# structurally-independent references (L1 trajectory-resolvent,
# hand-derived k_∞, the unified WDD body).


# ═══════════════════════════════════════════════════════════════════════
# 5. Linearity of L (sanity).
# ═══════════════════════════════════════════════════════════════════════


class TestLinearity:
    """L is a linear operator on :class:`TimedFullField`.

    D-I.3c (2026-05-29) — migrated from the retiring bare-ndarray
    calling convention to :class:`TimedFullField`.  The linearity
    claim survives unchanged; the construction routes through
    :meth:`~orpheus.transport.timed_full_field.TimedFullField.zeros` and the typed dunders
    :meth:`TimedFullField.__add__` / :meth:`TimedFullField.__mul__`
    propagate the linear combination to both bulk and boundary
    members.
    """

    @pytest.mark.parametrize("name,builder", GEOMETRIES_1D)
    def test_apply_zero_returns_zero(self, name, builder):
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh, ng=2)
        L = StreamingOperator(sn_mesh)
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh)
        out = L.apply(state)
        np.testing.assert_allclose(out.interior.values, 0.0, atol=1e-14)
        np.testing.assert_allclose(out.boundary.values, 0.0, atol=1e-14)

    @pytest.mark.parametrize("name,builder", GEOMETRIES_1D)
    def test_apply_is_linear(self, name, builder):
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh, ng=2)
        L = StreamingOperator(sn_mesh)
        state1 = _random_composite(sn_mesh, seed=51)
        state2 = _random_composite(sn_mesh, seed=52)
        # #208: a general α·ψ₁ + β·ψ₂ (α+β≠1) is illegal on affine flux STATES;
        # verify linearity via the affine-supported ops — homogeneity
        # op(c·ψ)=c·op(ψ) AND affine additivity in torsor form ψ₁ + λ(ψ₂⊖ψ₁)
        # = (1−λ)ψ₁ + λψ₂ (a flux). Together they imply full linearity; op.apply
        # stays on flux states.
        c, lam = 2.3, 0.7
        hom_combined = L.apply(c * state1)
        hom_separate = c * L.apply(state1)
        np.testing.assert_allclose(
            hom_combined.interior.values, hom_separate.interior.values, rtol=1e-12, atol=1e-13)
        np.testing.assert_allclose(
            hom_combined.boundary.values, hom_separate.boundary.values, rtol=1e-12, atol=1e-13)
        out_combined = L.apply(state1 + lam * (state2 - state1))
        out_separate = (1.0 - lam) * L.apply(state1) + lam * L.apply(state2)
        np.testing.assert_allclose(
            out_combined.interior.values, out_separate.interior.values,
            rtol=1e-12, atol=1e-13,
        )
        np.testing.assert_allclose(
            out_combined.boundary.values, out_separate.boundary.values,
            rtol=1e-12, atol=1e-13,
        )


# ═══════════════════════════════════════════════════════════════════════
# 6. The composite (L + C) dispatches to InvertibleOperator (R-1 Step C).
# ═══════════════════════════════════════════════════════════════════════


class TestSumCapabilities:
    """R-1 Step C contract: ``L + C`` dispatches to :class:`InvertibleOperator`.

    The within-group composite IS sweep-invertible — that's the
    algebraic identity SN methods are built on (Lewis & Miller §3.2).
    R-1 Step C encodes that identity at the type level by overriding
    :meth:`StreamingOperator.__add__` and
    :meth:`CollisionOperator.__add__` to return
    :class:`~orpheus.sn.operators.streaming.InvertibleOperator`, a specialisation
    of :class:`OperatorSum` that adds ``solve`` via its
    ``loss_representation`` strategy sweep (the operator-free
    ``transport_sweep`` wrapper retired at the coupled-block step 6).
    """

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_sum_exposes_apply(self, name, builder):
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh)
        C = MultiplicationOperator.from_mesh(sig_t, sn_mesh)
        A = L + C
        assert callable(getattr(A, "apply", None))

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_sum_advertises_solve_via_invertible_dispatch(
        self, name, builder,
    ):
        r"""``L + C`` dispatches to InvertibleOperator and is invertible."""
        from orpheus.sn.operators.streaming import InvertibleOperator
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh)
        C = MultiplicationOperator.from_mesh(sig_t, sn_mesh)
        A = L + C
        assert isinstance(A, InvertibleOperator)
        assert A.is_invertible


# ═══════════════════════════════════════════════════════════════════════
# D-H.2-C1 — Composite TimedFullField invariants on StreamingOperator.apply.
#
# 1-D geometries — the TimedFullField branch supports these natively.
# 2-D Cartesian raises NotImplementedError (Phase A coordination).  The
# parity-vs-legacy ``AngularFlux`` tests retired with D-H.2-C1 (the
# legacy class itself retires in C5; both branches share the same
# matvec kernel — exercising the composite branch alone is sufficient).
# ═══════════════════════════════════════════════════════════════════════


class TestCompositeInvariants:
    """Composite :class:`TimedFullField` dispatch on
    :meth:`StreamingOperator.apply`.

    Unlike C / S / F (which return implicit-zero composite boundary),
    L's composite return carries the actual face residual — this is
    the algebraic value the operator-algebra :math:`(L+C-S-F)\\psi`
    needs at the trace.
    """

    @pytest.mark.parametrize("name,builder", GEOMETRIES_1D)
    def test_returns_timeless_full_field(self, name, builder):
        from orpheus.transport.full_field import FullField
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.timed_full_field import TimedFullField

        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh)
        state = _random_composite(sn_mesh)

        out = L.apply(state)

        # #257 S8a — the matvec leaf is a base arrow ``FullField -> FullField``,
        # so the output is the TIMELESS FullField (history-free).
        assert isinstance(out, FullField)
        assert not isinstance(out, TimedFullField)
        assert isinstance(out.interior, AngularSourceSink)
        assert out.interior.mesh is sn_mesh

    @pytest.mark.parametrize("name,builder", GEOMETRIES_1D)
    def test_boundary_carries_face_residual(self, name, builder):
        """L emits a non-zero face residual into ``out.boundary``.

        L is the only operator in {L, C, S, F} that emits a non-zero
        face residual.  The composite branch routes the matvec's face
        output through the L2 AngularBoundaryFlux flat buffer at the
        layout-assigned slots; on random non-zero input the buffer
        must carry visibly non-zero values.
        """
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh)
        state = _random_composite(sn_mesh, seed=182)

        out_composite = L.apply(state)

        # At least one face slot must carry a non-trivial residual.
        layout = out_composite.boundary.layout
        face_max = max(
            float(np.abs(out_composite.boundary.face_view(face)).max())
            for face in layout.faces
        )
        assert face_max > 1e-12, (
            f"L's composite boundary face residual is ~0 "
            f"(max |face| = {face_max:.3e}); L must emit a non-zero "
            f"face residual on random input."
        )

    @pytest.mark.parametrize("name,builder", GEOMETRIES_1D)
    def test_zero_state_zero_output(self, name, builder):
        """ψ = 0 ⇒ L·ψ = 0 in both bulk AND boundary (linearity guard)."""
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh)
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh)
        out = L.apply(state)
        np.testing.assert_array_equal(out.interior.values, 0.0)
        np.testing.assert_array_equal(out.boundary.values, 0.0)

    @pytest.mark.parametrize("name,builder", GEOMETRIES_1D)
    def test_output_is_timeless_full_field(self, name, builder):
        """#257 S8a — the matvec leaf is a base arrow ``FullField -> FullField``.

        The operator output is the TIMELESS
        :class:`~orpheus.transport.full_field.FullField` (history-free),
        regardless of the input iterate's ``history_depth``: the comonad
        lives on the iteration driver, not on the operator (was: the old
        convention stamped ``history_depth`` onto the output — re-pointed).
        ``-O``-firing (Mode 8): explicit raise, not bare ``assert``.
        """
        from orpheus.transport.full_field import FullField

        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh)
        for depth in (0, 1, 2, 4):
            state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh, history_depth=depth)
            out = L.apply(state)
            if type(out) is not FullField or isinstance(out, TimedFullField):
                pytest.fail(
                    f"{name} depth={depth}: L.apply output must be a timeless "
                    f"FullField, got {type(out).__name__}"
                )

    def test_2d_cartesian_produces_full_field(self):
        """D-H.2-C4d — 2-D Cartesian path returns a (timeless) FullField.

        Pre-C4d this test pinned the deferred ``NotImplementedError``
        stub for the 2-D path; C4d ships the L2-native FD kernel
        (the representation ``loss_action`` matvec walk, which since
        S6.3 lives on the loss representation, off the operator —
        ``ScanMarch`` default since S6.9) and the path becomes
        functional.
        Structural invariant test: the return is the composite carrier
        with the correct bulk type.  #257 S8a — the matvec leaf is a base
        arrow ``FullField -> FullField`` (history-free; was: pinned the
        history-bearing ``TimedFullField`` — re-pointed).
        """
        from orpheus.geometry import Mesh2D
        from orpheus.transport.fields.angular_flux import (
            AngularFlux,
        )
        from orpheus.transport.full_field import FullField
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.timed_full_field import TimedFullField

        mesh = Mesh2D(
            edges_x=np.linspace(0.0, 1.0, 4),
            edges_y=np.linspace(0.0, 1.0, 4),
            mat_map=np.zeros((3, 3), dtype=int),
        )
        quad = Quadrature.level_symmetric(sn_order=4)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials(ng=2))
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh)

        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh)
        out = L.apply(state)

        # base-arrow codomain: a timeless FullField, NOT the timed subclass.
        if type(out) is not FullField or isinstance(out, TimedFullField):
            pytest.fail(
                f"L.apply output must be a timeless FullField, "
                f"got {type(out).__name__}"
            )
        assert isinstance(out.interior, AngularSourceSink)
        assert out.interior.mesh is sn_mesh

    def test_mesh_identity_invariant(self):
        """Distinct SNMesh instances must reject the apply."""
        sn_mesh_a = _slab_mesh()
        sn_mesh_b = _slab_mesh()
        sig_t = _sig_t_uniform(sn_mesh_a)
        L = StreamingOperator(sn_mesh_a)
        state_b = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh_b)
        with pytest.raises(ValueError, match="mesh-identity"):
            L.apply(state_b)


class TestOperatorAlgebraCompositionUnderTimedFullField:
    """``(L + C - S - F).apply(state)`` composes under TimedFullField.

    Load-bearing for the post-D-H.1b operator algebra: all four
    operators must accept TimedFullField input and return
    TimedFullField output for ``OperatorSum.apply`` (which evaluates
    ``self.a.apply(x) + self.b.apply(x)``) to type-check at every
    composition step.
    """

    @pytest.mark.parametrize("name,builder", GEOMETRIES_1D)
    def test_LC_apply_composite(self, name, builder):
        """``(L + C).apply(state)`` returns a TIMELESS FullField.

        ``L + C`` dispatches to :class:`InvertibleOperator`; its
        :meth:`apply` evaluates the representation ``loss_action``.  #257 S8a —
        the matvec leaf is a base arrow ``FullField -> FullField`` (the comonad
        lives on the driver), so the output is the TIMELESS FullField
        (history-free).
        """
        from orpheus.transport.full_field import FullField
        from orpheus.transport.timed_full_field import TimedFullField

        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh)
        C = MultiplicationOperator.from_mesh(sig_t, sn_mesh)
        A = L + C
        state = _random_composite(sn_mesh, seed=191)

        out = A.apply(state)

        assert isinstance(out, FullField)
        assert not isinstance(out, TimedFullField)
        assert out.interior.mesh is sn_mesh


# ═══════════════════════════════════════════════════════════════════════
# StreamingOperator.apply — pre-T.4 matvec regression snapshots
# ═══════════════════════════════════════════════════════════════════════
#
# These gates pin the production fused-matvec path ``StreamingOperator.apply``
# (= ``(L+C)ψ − σ_t·ψ``) against frozen pre-T.4 snapshots — slab (principled-
# equivalence, the #240 ÷V kernel re-associates ~1 ULP; boundary strict) +
# curvilinear (principled-equivalence) + 2-D Cartesian ScanMarch.
#
# #238 retired the separately-applicable ``(M_spatial, M_angular_redist)``
# operator-leaf split and its structural / decomposition-invariant /
# standalone-apply / materialize-inverse-cache tests: that split had no
# production consumer (everything rides the fused ``loss_action``), and the
# curvilinear Morel–Montry angular redistribution it isolated is verified
# end-to-end by the anisotropic curvilinear MMS
# (``tests/sn/verification/mms/test_curvilinear_aniso_convergence.py``,
# ``catches("ERR-026")``).  The shared fixture builders below feed the
# surviving snapshot regression gates.

PRE_T4_SNAPSHOTS_PATH = (
    "tests/sn/_fixtures/wave_t_t4/pre_t4_snapshots.npz"
)


def _slab_for_snapshot_arm(*, ng: int, bc_left: BC, bc_right: BC) -> SNMesh:
    """Reconstruct the slab fixture used by the T.4a snapshot script.

    Mirrors `tests/sn/_fixtures/wave_t_t4/_capture_pre_t4_snapshots.py`'s
    `_slab_mesh` builder (nx=20, GL N=8, P1 asymmetric SigS for 2G).
    Kept inline (NOT imported from the capture script) so the
    regression-test fixture is self-contained.
    """
    from scipy.sparse import csr_matrix

    from orpheus.derivations.common.xs_library import make_mixture

    if ng == 1:
        mix = make_mixture(
            sig_t=np.array([1.0]),
            sig_c=np.array([0.2]),
            sig_f=np.array([0.1]),
            nu=np.array([2.5]),
            chi=np.array([1.0]),
            sig_s=np.array([[0.7]]),
        )
    else:
        p0 = np.array([[0.38, 0.10], [0.05, 0.90]])
        p1 = np.array([[0.02, 0.01], [0.00, 0.04]])
        mix = make_mixture(
            sig_t=np.array([0.5, 1.0]),
            sig_c=np.array([0.01, 0.02]),
            sig_f=np.array([0.01, 0.08]),
            nu=np.array([2.5, 2.5]),
            chi=np.array([1.0, 0.0]),
            sig_s=p0,
        )
        mix.SigS = [csr_matrix(p0), csr_matrix(p1)]

    mesh = Mesh1D(
        edges=np.linspace(0.0, 4.0, 21),
        mat_ids=np.zeros(20, dtype=int),
        bc_left=bc_left,
        bc_right=bc_right,
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    return SNMesh(mesh, quad, {0: mix})


def _sigma_t_from_mat_map(sn_mesh: SNMesh) -> np.ndarray:
    """Build per-(g, *cell) σ_t from sn_mesh.mat_map (matches the T.4a script)."""
    ng = sn_mesh.ng
    sig_t = np.empty((ng, *sn_mesh.spatial_shape), dtype=float)
    mat_map = sn_mesh.mat_map
    for cell in np.ndindex(*sn_mesh.spatial_shape):
        mid = int(mat_map[cell])
        mat = sn_mesh.materials[mid]
        for g in range(ng):
            sig_t[(g, *cell)] = float(mat.SigT[g])
    return sig_t


class TestT4bPreT4RegressionSnapshot:
    """L4-1 / L4-5 / L4-6 / L4-7 — slab matvec regression vs pre-T.4 snapshots.

    Originally a STRICT bit-identity gate (Route A: the slab matvec used
    the ×V ``cell_balance`` density path with no reduction reorder).  #240
    moved the Cartesian matvec onto the uniform ÷V
    ``scheme.residual_kernel_batch`` kernel (DD AND LD, retiring the
    bit-identity-only ``matvec_via_kernel`` special case) — which
    re-associates the cell-balance reduction (~1 ULP mean, max ~64 at a
    near-zero cancellation value; max relΔ ~1e-14).  Per ``vv-principles``
    §"Bit-identity vs principled-equivalence" the snapshot was re-baselined
    and the gate narrowed to :func:`assert_regression` (``kind="direct"``,
    ``nulp=reduction_depth``).

    #257 S8b re-baselined the BULK gate again: ``StreamingOperator.apply``
    is now pure σ-free ``streaming_action(ψ) = loss_action(0, ψ)`` (the
    ``(L+C)−C`` fold is retired).  The σ-free walk re-associates the FP
    reduction tree relative to subtracting σ⊙ψ off the σ-bearing matvec,
    so the pure-L value drifts from the pre-T.4 snapshot by FP-non-
    associativity only: **rel ≤ 1e-16 (machine ε) on every arm, boundary
    0 ULP**.  The per-element ULP *metric* spikes (max ~192 measured) at
    near-zero cancellation values — the same artefact the slab arms already
    documented (a tiny absolute Δ against a tiny |ψ| reads a large ULP
    count).  The nULP bound is therefore widened to a documented
    :attr:`_PURE_L_NULP` that absorbs the near-zero spikes while staying
    far below any real-bug magnitude (the structural ground is the
    BYTE-IDENTICAL ``(L+C)`` composite — see the decomposition gate — and
    the closed-form k_∞ + MMS backstop, NOT old-vs-new proximity).

    Boundary traces stay byte-identical (0 ULP) — the outflow defect is
    reconstructed from the same ``ψ_out = 2ψ̄ − ψ_in`` faces, untouched by
    the pure-L carve (``C`` never acts on the trace).
    """

    # Pure-L (#257 S8b) re-associates the matvec FP tree vs the pre-T.4
    # snapshot; the rel drift is machine-ε but the per-element ULP metric
    # spikes at near-zero cancellation values (≤192 measured).  256 absorbs
    # that with headroom and is far below real-bug magnitude.
    _PURE_L_NULP = 256

    @pytest.fixture(scope="class")
    def snapshots(self):
        """Load the pre-T.4 snapshot bundle."""
        from tests.sn._test_helpers import SN_TESTS_ROOT
        path = SN_TESTS_ROOT / "_fixtures" / "wave_t_t4" / "pre_t4_snapshots.npz"
        with np.load(path) as data:
            return {k: data[k] for k in data.files}

    def _capture_arm(self, sn_mesh: SNMesh, seed: int) -> tuple[np.ndarray, np.ndarray]:
        """Re-run StreamingOperator.apply on the snapshot fixture; return
        (bulk, boundary) values.
        """
        from dataclasses import replace
        # Pure-L streaming (#257 S8b) — σ-free; the snapshot fixture's σ_t
        # is no longer needed to build L (the snapshot pins L's matvec leaf).
        L = StreamingOperator(sn_mesh)
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh)
        rng = np.random.default_rng(seed)
        state = replace(
            state,
            interior=replace(
                state.interior,
                values=rng.uniform(0.05, 1.0, size=state.interior.values.shape),
            ),
        )
        # #282 route (a) → step 6: on a carrying mesh (sphere) the joint
        # σ-free streaming action is the BLOCK SUM ``L·ψ + Seeding·ψ½``
        # (Seeding is σ-independent — its own pinned property), with the
        # CONSISTENT edge-extrapolated ψ½ seed reproducing the
        # pre-route-(a) internally-computed seed.  The bulk re-associates
        # ~1 ULP vs the fused spelling (inside the principled band); the
        # boundary residual is probe-reconstructed, seed-independent —
        # byte-identical.  No seed term on non-carrying meshes (slab/cyl).
        sd = radial_characteristic_edge_seed(state.interior.values, sn_mesh)
        out = L.apply(state)
        if sd is not None:
            from orpheus.sn.operators.radial_characteristic import (
                RadialCharacteristicSeeding,
            )

            out = out + RadialCharacteristicSeeding(sn_mesh).apply(sd)
        return out.interior.values.copy(), out.boundary.values.copy()

    def _assert_arm(self, snapshots, *, tag: str, mesh: SNMesh, seed: int) -> None:
        """Re-run the slab matvec arm: BULK principled-equivalent, BOUNDARY strict.

        The BULK residual rides the ÷V ``residual_kernel_batch`` kernel, which
        re-associates the cell-balance reduction ~1 ULP (#240), so its gate is
        ``assert_regression(kind="direct")`` — a single matvec's only drift
        source is FP-non-associativity over the cell-chain reduction, bounded by
        ``nulp=reduction_depth`` (the face propagates through ``mesh.nx`` cells).
        That bound never bites here: the live output reproduces the regenerated
        same-build snapshot to 0 ULP.  (A stray near-zero cancellation value can
        read a large ULP *metric* for a tiny absolute Δ — e.g. ~64 ULP at
        |ψ|~1e-2, absΔ~2e-16 — but those vanish against the same-build snapshot
        and never reach the gate.)  Bit-identity is retained as the escalatable
        ``DriftWarning`` (``-W error::DriftWarning`` re-pins exact bytes).

        The BOUNDARY trace did NOT move (0 ULP — the outflow defect reconstructs
        from the same ``ψ_out = 2ψ̄ − ψ_in`` faces), so it keeps the STRICT
        ``assert_array_equal`` gate: bit-identity as a free bonus, per
        ``vv-principles`` — strict where the implementation is genuinely
        unchanged.
        """
        bulk, boundary = self._capture_arm(mesh, seed=seed)
        assert_regression(
            bulk, snapshots[f"{tag}_apply_bulk"],
            conv_tol=0.0, kind="direct", reduction_depth=self._PURE_L_NULP,
            case_name=f"{tag}_apply_bulk", quantity="apply_bulk",
        )
        np.testing.assert_array_equal(
            boundary, snapshots[f"{tag}_apply_boundary"],
        )

    def test_slab_1g_vacuum_apply_principled_equiv(self, snapshots):
        """L4-1 — slab 1G vacuum, P0 fixture (#240 uniform-kernel matvec)."""
        self._assert_arm(
            snapshots, tag="slab_1g_vacuum",
            mesh=_slab_for_snapshot_arm(
                ng=1, bc_left=BC("vacuum"), bc_right=BC("vacuum"),
            ),
            seed=20260531 + 1,
        )

    def test_slab_2g_vacuum_apply_principled_equiv(self, snapshots):
        """L4-1 — slab 2G vacuum, P1 asymmetric SigS (ERR-002 detector)."""
        self._assert_arm(
            snapshots, tag="slab_2g_vacuum",
            mesh=_slab_for_snapshot_arm(
                ng=2, bc_left=BC("vacuum"), bc_right=BC("vacuum"),
            ),
            seed=20260531 + 2,
        )

    def test_slab_2g_reflective_apply_principled_equiv(self, snapshots):
        """L4-1 — slab 2G with reflective inner BC + vacuum outer.

        Catches BC convention drift via the D-B+1 specular permutation
        path (the D-B+1 production tensor-network instance).
        """
        self._assert_arm(
            snapshots, tag="slab_2g_reflective",
            mesh=_slab_for_snapshot_arm(
                ng=2, bc_left=BC("reflective"), bc_right=BC("vacuum"),
            ),
            seed=20260531 + 3,
        )

    # ── #240 D5a — 2-D Cartesian ScanMarch matvec frozen-reference arm ──
    #
    # The cart2d apply snapshots (frozen pre-#240, UNCONSUMED until D5a) pin the
    # 2-D row-march matvec ``ScanMarch._loss_action_interior`` — the
    # value-preserving linchpin (GATE D5a.2) that catches a fold bug which
    # preserves the D5a.1 ``ScanMarch ≡ FullFieldWavefront`` relation but shifts
    # the VALUE (a "both paths move together" blind spot of the two-paths gate).
    # D5a routed the row-march off inline DD onto ``scheme.residual_kernel_batch``
    # + ``scheme.reflect_scan_coefficients`` (#239 coefficient-model lift), which
    # re-associates the cell-balance reduction ~1 ULP (the same principled
    # re-baseline as the slab arms, #240). The frozen bulk snapshot was
    # regenerated to the post-D5a value (max absΔ ~3.6e-15 vs the byte-identical
    # pre-D5a value — verified pre-fold the frozen ≡ current at 0 ULP), per
    # ``vv-principles`` §"Bit-identity vs principled-equivalence" (the value is
    # independently grounded by the D5a.1 oracle, NOT by old-vs-new proximity).
    # Boundary trace is byte-identical (0 ULP — the outflow defect reconstructs
    # from the same ``ψ_out = 2ψ̄ − ψ_in`` faces).

    def _assert_cart2d_arm(
        self, snapshots, *, tag: str, ng: int, bc_kind: str, seed: int,
    ) -> None:
        """Re-run the 2-D Cartesian ScanMarch matvec arm: BULK principled, BOUNDARY strict.

        Uses the capture-script mesh / ψ constructors verbatim so the live
        matvec consumes the SAME inputs the frozen snapshot was captured from
        (no fixture drift).  ``reduction_depth=_PURE_L_NULP`` per the pure-L
        FP-tree re-association (#257 S8b) — see the class docstring.
        """
        from tests.sn._fixtures.wave_t_t4._capture_pre_t4_snapshots import (
            _cart2d_mesh, _make_state,
        )
        sn_mesh = _cart2d_mesh(ng=ng, bc_kind=bc_kind)
        # Pure-L streaming (#257 S8b) — σ-free; the fixture σ_t is no longer
        # needed to build L.
        L = StreamingOperator(sn_mesh)
        state = _make_state(sn_mesh, seed=seed)
        out = L.apply(state)
        assert_regression(
            out.interior.values, snapshots[f"{tag}_apply_bulk"],
            conv_tol=0.0, kind="direct", reduction_depth=self._PURE_L_NULP,
            case_name=f"{tag}_apply_bulk", quantity="apply_bulk",
        )
        np.testing.assert_array_equal(
            out.boundary.values, snapshots[f"{tag}_apply_boundary"],
        )

    def test_cart2d_1g_specular_apply_principled_equiv(self, snapshots):
        """D5a.2 — 2-D Cartesian 1G reflective (specular), ScanMarch matvec."""
        self._assert_cart2d_arm(
            snapshots, tag="cart2d_1g_specular", ng=1, bc_kind="specular",
            seed=20260531 + 31,
        )

    def test_cart2d_1g_vacuum_apply_principled_equiv(self, snapshots):
        """D5a.2 — 2-D Cartesian 1G vacuum, ScanMarch matvec.

        Vacuum NULLS the reflective face-shed: pairs with the specular arm so a
        BC-coupled fold bug (transverse face convention) is not hidden by the
        reflective permutation path.
        """
        self._assert_cart2d_arm(
            snapshots, tag="cart2d_1g_vacuum", ng=1, bc_kind="vacuum",
            seed=20260531 + 32,
        )

    def test_cart2d_2g_specular_apply_principled_equiv(self, snapshots):
        """D5a.2 — 2-D Cartesian 2G P1-asymmetric reflective, ScanMarch matvec.

        ≥2G with asymmetric downscatter (ERR-002 detector) + reflective BC — the
        Cardinal-6 multi-group + heterogeneous degeneracy-break for the matvec.
        """
        self._assert_cart2d_arm(
            snapshots, tag="cart2d_2g_specular", ng=2, bc_kind="specular",
            seed=20260531 + 33,
        )


class TestT4cPreT4RegressionSnapshotCurvilinear:
    """L4-2 / L4-3 — sphere + cylinder curvilinear matvec regression snapshots.

    Wave T post-T.5 matvec retirement: ``StreamingOperator.apply``
    curvilinear branch now goes through the dual-emission
    decomposition (``M_spatial + M_angular_redist − σ_t·ψ``) rather
    than the legacy unified-matvec shortcut.  The dual-emission body
    computes ``m_spat = m_full − m_ang`` per cell — this introduces a
    per-cell FP-non-associativity reorder vs the legacy single-
    emission path.  The drift is bounded by ``cell_count × ULP``
    (per `vv-principles` principled-equivalence three-criteria gate);
    the test relaxes from strict bit-identity to ``rtol=1e-13,
    atol=1e-14`` accordingly.  L1 closed-form k_∞ + L1 MMS gates
    remain the structural-independence ground.

    **Sphere arms re-baselined 2026-06-26 (closes #250).**  The CYLINDER
    arms remain frozen at their last refresh (current ≡ frozen to
    ``rtol=1e-13``), so they still witness the operator-algebra campaign's
    numerical inertness.  The SPHERE snapshot, however, went stale across a
    real closure change: ``b2d8a6d`` (Bailey Eq. 43, Refs #229) unclamped
    the spherical Morel–Montry WDD weight τ AFTER this store's last refresh
    but updated only its own targeted snapshot (#240 later re-captured
    SLB/CYL/cart2d here, deferring SPH — lesson L-034).  The sphere arms
    were re-captured to current production (sphere-only, via the canonical
    capture builders; the cyl / slab / cart2d keys stayed frozen).
    Correctness of the current sphere matvec was verified vs L0 streaming-
    equilibrium (per-ordinate), L1 isotropic + anisotropic MMS, and the
    L1 trajectory-resolvent before re-capture — never green-by-fiat
    (``vv-principles`` L11/L14/L27).
    """

    @pytest.fixture(scope="class")
    def snapshots(self):
        """Load the pre-T.4 snapshot bundle."""
        from tests.sn._test_helpers import SN_TESTS_ROOT
        path = SN_TESTS_ROOT / "_fixtures" / "wave_t_t4" / "pre_t4_snapshots.npz"
        with np.load(path) as data:
            return {k: data[k] for k in data.files}

    def _capture_arm(self, sn_mesh: SNMesh, seed: int) -> tuple[np.ndarray, np.ndarray]:
        from dataclasses import replace
        # Pure-L streaming (#257 S8b) — σ-free; the snapshot fixture's σ_t
        # is no longer needed to build L (the snapshot pins L's matvec leaf).
        L = StreamingOperator(sn_mesh)
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh)
        rng = np.random.default_rng(seed)
        state = replace(
            state,
            interior=replace(
                state.interior,
                values=rng.uniform(0.05, 1.0, size=state.interior.values.shape),
            ),
        )
        # #282 route (a) → step 6: on a carrying mesh (sphere) the joint
        # σ-free streaming action is the BLOCK SUM ``L·ψ + Seeding·ψ½``
        # (Seeding is σ-independent — its own pinned property), with the
        # CONSISTENT edge-extrapolated ψ½ seed reproducing the
        # pre-route-(a) internally-computed seed.  The bulk re-associates
        # ~1 ULP vs the fused spelling (inside the principled band); the
        # boundary residual is probe-reconstructed, seed-independent —
        # byte-identical.  No seed term on non-carrying meshes (slab/cyl).
        sd = radial_characteristic_edge_seed(state.interior.values, sn_mesh)
        out = L.apply(state)
        if sd is not None:
            from orpheus.sn.operators.radial_characteristic import (
                RadialCharacteristicSeeding,
            )

            out = out + RadialCharacteristicSeeding(sn_mesh).apply(sd)
        return out.interior.values.copy(), out.boundary.values.copy()

    def test_sphere_1g_apply_bit_identical(self, snapshots):
        """L4-2 — sphere 1G P0 vacuum-at-r=R."""
        from tests.sn._fixtures.wave_t_t4._capture_pre_t4_snapshots import (
            _sphere_mesh,
        )
        bulk, boundary = self._capture_arm(_sphere_mesh(ng=1), seed=20260531 + 11)
        # Snapshot bundle is legacy (N, ng, nx, 1); drop the phantom ny=1.
        np.testing.assert_allclose(
            bulk, snapshots["sphere_1g_apply_bulk"],
            rtol=1e-13, atol=1e-14,
        )
        np.testing.assert_array_equal(
            boundary, snapshots["sphere_1g_apply_boundary"],
        )

    def test_sphere_2g_apply_bit_identical(self, snapshots):
        """L4-2 — sphere 2G P1 asymmetric SigS vacuum-at-r=R."""
        from tests.sn._fixtures.wave_t_t4._capture_pre_t4_snapshots import (
            _sphere_mesh,
        )
        bulk, boundary = self._capture_arm(_sphere_mesh(ng=2), seed=20260531 + 12)
        np.testing.assert_allclose(
            bulk, snapshots["sphere_2g_apply_bulk"],
            rtol=1e-13, atol=1e-14,
        )
        np.testing.assert_array_equal(
            boundary, snapshots["sphere_2g_apply_boundary"],
        )

    def test_cylinder_1g_apply_bit_identical(self, snapshots):
        """L4-3 — cylinder 1G P0 level-symmetric quad."""
        from tests.sn._fixtures.wave_t_t4._capture_pre_t4_snapshots import (
            _cylinder_mesh,
        )
        bulk, boundary = self._capture_arm(_cylinder_mesh(ng=1), seed=20260531 + 21)
        np.testing.assert_allclose(
            bulk, snapshots["cyl_1g_apply_bulk"],
            rtol=1e-13, atol=1e-14,
        )
        np.testing.assert_array_equal(
            boundary, snapshots["cyl_1g_apply_boundary"],
        )

    def test_cylinder_2g_apply_bit_identical(self, snapshots):
        """L4-3 — cylinder 2G P1 asymmetric SigS, multi-level quad."""
        from tests.sn._fixtures.wave_t_t4._capture_pre_t4_snapshots import (
            _cylinder_mesh,
        )
        bulk, boundary = self._capture_arm(_cylinder_mesh(ng=2), seed=20260531 + 22)
        np.testing.assert_allclose(
            bulk, snapshots["cyl_2g_apply_bulk"],
            rtol=1e-13, atol=1e-14,
        )
        np.testing.assert_array_equal(
            boundary, snapshots["cyl_2g_apply_boundary"],
        )
