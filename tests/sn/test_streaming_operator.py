r"""Foundation tests for :class:`orpheus.sn.operator.StreamingOperator`.

Phase G Step 3+4.b.i (Issue #196, Resolution A). The "L" of the
four-operator algebra ``A_wg = L + C - S.foldable_part()``. This
substep lands it as an ADDITIVE leaf alongside the legacy
:class:`SNStreamingOperator` (which bundles ``L + C`` with σ_t);
no production code consumes the leaf yet (substep 3+4.c wires
:class:`SNSolver`).

Resolution A — subtractive definition
-------------------------------------

.. math::

   L.{\rm apply}(\psi) \;:=\; M(\psi;\;\sigma_t) \;-\;
                              \sigma_t \odot \psi_{\rm packed}

L carries σ_t at constructor time. This is intrinsic to the discrete
curvilinear matvec (rational in σ_t through Hébert §3.9.4's Carlson
coupled-pole seed), not a defect — analogous to the DD coefficient
:math:`\alpha_{DD}(\sigma_t\,\Delta x)` carrying σ_t in
characteristic-line methods.

Load-bearing test (no xfail allowed)
------------------------------------

The composition equivalence
``(L + C(σ_t)).apply(ψ) ≈ SNStreamingOperator(σ_t).apply(ψ)`` MUST
hold bit-exact (``rtol=1e-12``) on ALL THREE geometries: slab,
sphere, cylinder. With both homogeneous AND heterogeneous σ_t. By
Resolution A this is guaranteed by construction.

The decomposition itself (``rel_residual = 0.0``) lives in
:file:`tests/sn/test_streaming_operator_decomposition.py`.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    CAP_SOLVE,
    LinearOperator,
)
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import (
    CollisionOperator,
    StreamingOperator,
)
from orpheus.numerics.quadrature import Quadrature
from tests.sn._test_helpers import placeholder_materials

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
    """σ_t uniform across cells / groups (PR-INDEX-3: ``(ng, nx, ny)``)."""
    return value * np.ones((ng, sn_mesh.nx, sn_mesh.ny))


def _sig_t_heterogeneous(sn_mesh: SNMesh, ng: int = 2) -> np.ndarray:
    """σ_t heterogeneous — different value per (cell, group)."""
    nx, ny = sn_mesh.nx, sn_mesh.ny
    rng = np.random.default_rng(seed=20260514)
    return 0.3 + 0.5 * rng.random((ng, nx, ny))


def _packed_psi(sn_mesh: SNMesh, ng: int, seed: int = 42) -> np.ndarray:
    """Random packed-vector ψ matching :class:`StreamingOperator`'s
    eq_map size.

    PR-TYPED-6.5 Phase 3b — :class:`StreamingOperator` now uses the
    B1'' face-aware packed layout for 1-D (cell-centres + outer-face
    + inner-face for slab); :class:`SNStreamingOperator` retains the
    legacy compressed layout (no face slots, curvilinear inward-at-
    outer cell-centres BC-resolved).  ``_packed_psi`` sizes for the
    new leaf; tests that compare against the legacy bundle need
    their own legacy-sized random vector (or accept that B1'' fixes
    the bug and the two layouts no longer overlap — see
    :class:`TestCompositionEquivalence` xfail markers).
    """
    sig_t = _sig_t_uniform(sn_mesh, ng=ng)
    L = StreamingOperator(sn_mesh, sig_t)
    rng = np.random.default_rng(seed)
    return rng.standard_normal(L.n_unknowns)


GEOMETRIES = [
    ("slab", _slab_mesh),
    ("sphere", _spherical_mesh),
    ("cylinder", _cylindrical_mesh),
]


# ═══════════════════════════════════════════════════════════════════════
# 1. Capability advertising.
# ═══════════════════════════════════════════════════════════════════════


class TestCapabilities:
    """StreamingOperator advertises {apply} only — no solve, no apply_T.

    Resolution A: L alone is not invertible (rank-deficient without
    collision). solve appears at the OperatorSum level via the
    fusion hook (substep 3+4.b.ii).
    """

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_capabilities_apply_only(self, name, builder):
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        assert L.capabilities == frozenset({CAP_APPLY})

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_no_solve_capability(self, name, builder):
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        assert CAP_SOLVE not in L.capabilities

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_no_apply_transpose_capability(self, name, builder):
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        assert CAP_APPLY_TRANSPOSE not in L.capabilities

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_satisfies_linear_operator_protocol(self, name, builder):
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        assert isinstance(L, LinearOperator)


# ═══════════════════════════════════════════════════════════════════════
# 2. Constructor surface — Resolution A requires σ_t at construction.
# ═══════════════════════════════════════════════════════════════════════


class TestConstructor:
    """StreamingOperator(sn_mesh, sigma_t) takes σ_t as a required arg.

    Pattern 4 (illegal states unrepresentable): the discrete L is
    intrinsically σ_t-coupled (Hébert §3.9.4 Carlson seed); a
    constructor without σ_t would silently produce a different
    operator. Resolution A's contract requires σ_t at construction.
    """

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_construct_with_sn_mesh_and_sigma_t(self, name, builder):
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        assert L.sn_mesh is sn_mesh
        assert L.sigma_t is sig_t


# ═══════════════════════════════════════════════════════════════════════
# 3. apply shape correctness.
# ═══════════════════════════════════════════════════════════════════════


class TestApplyShape:
    """apply preserves packed-vector shape across geometries."""

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_apply_preserves_shape(self, name, builder):
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh, ng=2)
        L = StreamingOperator(sn_mesh, sig_t)
        psi = _packed_psi(sn_mesh, ng=2)
        out = L.apply(psi)
        assert out.shape == psi.shape

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_apply_returns_packed_vector(self, name, builder):
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh, ng=2)
        L = StreamingOperator(sn_mesh, sig_t)
        psi = _packed_psi(sn_mesh, ng=2)
        out = L.apply(psi)
        assert out.ndim == 1


# ═══════════════════════════════════════════════════════════════════════
# 4. THE LOAD-BEARING COMPOSITION TEST — no xfail allowed.
# ═══════════════════════════════════════════════════════════════════════


# D-K.5 (2026-05-29): ``TestCompositionEquivalence`` retired together
# with :class:`SNStreamingOperator`.  The three xfail-strict tests
# (``test_uniform_sigma_t_homogeneous``, ``test_heterogeneous_sigma_t``,
# ``test_via_operator_sum_algebra``) pinned
# ``(L + C).apply(ψ) ≈ SNStreamingOperator(σ_t).apply(ψ)``.  The
# comparison was documented as retiring "when SNStreamingOperator goes
# away at PR-TYPED-6c Step 7"; that step is D-K.  ``(L + C).apply``'s
# correctness is now gated by ``test_unified_matvec_{slab,sphere,
# cylinder}.py`` against structurally-independent references
# (L1 trajectory-resolvent, hand-derived k_∞, the unified WDD body
# itself).


# ═══════════════════════════════════════════════════════════════════════
# 5. Linearity of L (sanity).
# ═══════════════════════════════════════════════════════════════════════


class TestLinearity:
    """L is a linear operator on packed vectors."""

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_apply_zero_returns_zero(self, name, builder):
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh, ng=2)
        L = StreamingOperator(sn_mesh, sig_t)
        psi = _packed_psi(sn_mesh, ng=2)
        zero = np.zeros_like(psi)
        out = L.apply(zero)
        np.testing.assert_allclose(out, 0.0, atol=1e-14)

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_apply_is_linear(self, name, builder):
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh, ng=2)
        L = StreamingOperator(sn_mesh, sig_t)
        psi1 = _packed_psi(sn_mesh, ng=2, seed=51)
        psi2 = _packed_psi(sn_mesh, ng=2, seed=52)
        alpha = 2.3
        beta = -0.7
        out_combined = L.apply(alpha * psi1 + beta * psi2)
        out_separate = alpha * L.apply(psi1) + beta * L.apply(psi2)
        np.testing.assert_allclose(
            out_combined, out_separate, rtol=1e-12, atol=1e-13,
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
    :class:`~orpheus.sn.operator.InvertibleOperator`, a specialisation
    of :class:`OperatorSum` that adds ``solve`` via
    :func:`~orpheus.sn.sweep.transport_sweep`.
    """

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_sum_advertises_apply(self, name, builder):
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        C = CollisionOperator(sn_mesh, sig_t)
        A = L + C
        assert CAP_APPLY in A.capabilities

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_sum_advertises_solve_via_invertible_dispatch(
        self, name, builder,
    ):
        r"""``L + C`` dispatches to InvertibleOperator and advertises solve."""
        from orpheus.sn.operator import InvertibleOperator
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        C = CollisionOperator(sn_mesh, sig_t)
        A = L + C
        assert isinstance(A, InvertibleOperator)
        assert CAP_SOLVE in A.capabilities


# ═══════════════════════════════════════════════════════════════════════
# D-H.2-C1 — Composite TimedFullField invariants on StreamingOperator.apply.
#
# 1-D geometries — the TimedFullField branch supports these natively.
# 2-D Cartesian raises NotImplementedError (Phase A coordination).  The
# parity-vs-legacy ``AngularFlux`` tests retired with D-H.2-C1 (the
# legacy class itself retires in C5; both branches share the same
# matvec kernel — exercising the composite branch alone is sufficient).
# ═══════════════════════════════════════════════════════════════════════


GEOMETRIES_1D = [
    ("slab", _slab_mesh),
    ("sphere", _spherical_mesh),
    ("cylinder", _cylindrical_mesh),
]


def _random_composite(sn_mesh, seed=171):
    """Build a TimedFullField with non-zero bulk + non-zero boundary."""
    from dataclasses import replace

    rng = np.random.default_rng(seed)
    state = sn_mesh.zeros_timed_full_field()
    bulk_values = rng.standard_normal(state.bulk.values.shape)
    boundary_values = 0.1 + rng.random(state.boundary.values.shape)
    state = replace(state, bulk=replace(state.bulk, values=bulk_values))
    state = replace(
        state, boundary=replace(state.boundary, values=boundary_values),
    )
    return state


class TestCompositeInvariants:
    """Composite :class:`TimedFullField` dispatch on
    :meth:`StreamingOperator.apply`.

    Unlike C / S / F (which return implicit-zero composite boundary),
    L's composite return carries the actual face residual — this is
    the algebraic value the operator-algebra :math:`(L+C-S-F)\\psi`
    needs at the trace.
    """

    @pytest.mark.parametrize("name,builder", GEOMETRIES_1D)
    def test_returns_timed_full_field(self, name, builder):
        from orpheus.transport.fields.angular_flux import (
            AngularFlux,
        )
        from orpheus.transport.timed_full_field import TimedFullField

        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        state = _random_composite(sn_mesh)

        out = L.apply(state)

        assert isinstance(out, TimedFullField)
        assert isinstance(out.bulk, AngularFlux)
        assert out.bulk.mesh is sn_mesh
        assert out.history_depth == state.history_depth
        assert out._history == ()

    @pytest.mark.parametrize("name,builder", GEOMETRIES_1D)
    def test_boundary_carries_face_residual(self, name, builder):
        """L emits a non-zero face residual into ``out.boundary``.

        L is the only operator in {L, C, S, F} that emits a non-zero
        face residual.  The composite branch routes the matvec's face
        output through the L2 BoundaryFlux flat buffer at the
        layout-assigned slots; on random non-zero input the buffer
        must carry visibly non-zero values.
        """
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
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
        L = StreamingOperator(sn_mesh, sig_t)
        state = sn_mesh.zeros_timed_full_field()
        out = L.apply(state)
        np.testing.assert_array_equal(out.bulk.values, 0.0)
        np.testing.assert_array_equal(out.boundary.values, 0.0)

    @pytest.mark.parametrize("name,builder", GEOMETRIES_1D)
    def test_history_depth_preserved(self, name, builder):
        """Composite return preserves input ``history_depth`` capacity."""
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        for depth in (0, 1, 2, 4):
            state = sn_mesh.zeros_timed_full_field(history_depth=depth)
            out = L.apply(state)
            assert out.history_depth == depth

    def test_2d_cartesian_produces_timed_full_field(self):
        """D-H.2-C4d — 2-D Cartesian path returns a TimedFullField.

        Pre-C4d this test pinned the deferred ``NotImplementedError``
        stub for the 2-D path; C4d ships the L2-native FD kernel
        ``_apply_2d_cartesian`` and the path becomes functional.
        Structural invariant test: the return is the composite carrier
        with the correct bulk type.
        """
        from orpheus.geometry import Mesh2D
        from orpheus.transport.fields.angular_flux import (
            AngularFlux,
        )
        from orpheus.transport.timed_full_field import TimedFullField

        mesh = Mesh2D(
            edges_x=np.linspace(0.0, 1.0, 4),
            edges_y=np.linspace(0.0, 1.0, 4),
            mat_map=np.zeros((3, 3), dtype=int),
        )
        quad = Quadrature.level_symmetric(sn_order=4)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials(ng=2))
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)

        state = sn_mesh.zeros_timed_full_field()
        out = L.apply(state)

        assert isinstance(out, TimedFullField)
        assert isinstance(out.bulk, AngularFlux)
        assert out.bulk.mesh is sn_mesh
        assert out.history_depth == state.history_depth

    def test_mesh_identity_invariant(self):
        """Distinct SNMesh instances must reject the apply."""
        sn_mesh_a = _slab_mesh()
        sn_mesh_b = _slab_mesh()
        sig_t = _sig_t_uniform(sn_mesh_a)
        L = StreamingOperator(sn_mesh_a, sig_t)
        state_b = sn_mesh_b.zeros_timed_full_field()
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
        """``(L + C).apply(state)`` returns a TimedFullField.

        ``L + C`` dispatches to :class:`InvertibleOperator` (subclass
        of :class:`OperatorSum`); its inherited ``.apply`` evaluates
        ``L.apply(state) + C.apply(state)`` and the addition succeeds
        only if BOTH return TimedFullField.
        """
        from orpheus.transport.timed_full_field import TimedFullField

        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        C = CollisionOperator(sn_mesh, sig_t)
        A = L + C
        state = _random_composite(sn_mesh, seed=191)

        out = A.apply(state)

        assert isinstance(out, TimedFullField)
        assert out.bulk.mesh is sn_mesh
        assert out.history_depth == state.history_depth

