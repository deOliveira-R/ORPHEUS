r"""Foundation tests for :class:`orpheus.sn.operator.StreamingOperator`.

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
from tests.sn._test_helpers import _LC_matvec
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
    state = sn_mesh.zeros_timed_full_field()
    bulk_values = rng.standard_normal(state.bulk.values.shape)
    boundary_values = 0.1 + rng.random(state.boundary.values.shape)
    state = replace(state, bulk=replace(state.bulk, values=bulk_values))
    state = replace(
        state, boundary=replace(state.boundary, values=boundary_values),
    )
    return state


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
    :meth:`SNMesh.zeros_timed_full_field` and the typed dunders
    :meth:`TimedFullField.__add__` / :meth:`TimedFullField.__mul__`
    propagate the linear combination to both bulk and boundary
    members.
    """

    @pytest.mark.parametrize("name,builder", GEOMETRIES_1D)
    def test_apply_zero_returns_zero(self, name, builder):
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh, ng=2)
        L = StreamingOperator(sn_mesh, sig_t)
        state = sn_mesh.zeros_timed_full_field()
        out = L.apply(state)
        np.testing.assert_allclose(out.bulk.values, 0.0, atol=1e-14)
        np.testing.assert_allclose(out.boundary.values, 0.0, atol=1e-14)

    @pytest.mark.parametrize("name,builder", GEOMETRIES_1D)
    def test_apply_is_linear(self, name, builder):
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh, ng=2)
        L = StreamingOperator(sn_mesh, sig_t)
        state1 = _random_composite(sn_mesh, seed=51)
        state2 = _random_composite(sn_mesh, seed=52)
        alpha = 2.3
        beta = -0.7
        out_combined = L.apply(alpha * state1 + beta * state2)
        out_separate = alpha * L.apply(state1) + beta * L.apply(state2)
        np.testing.assert_allclose(
            out_combined.bulk.values, out_separate.bulk.values,
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


# ═══════════════════════════════════════════════════════════════════════
# Wave T T.4b — M_spatial / M_angular_redist decomposition tests
# ═══════════════════════════════════════════════════════════════════════
#
# Per `.claude/agent-memory/test-architect/wave_t_t4_streaming_verification_spec.md`
# §4 — the algebraic-identity tests for the new
# `StreamingOperator.M_spatial` + `M_angular_redist` decomposition.
#
# Spec deviation from §4 (resolved via user-endorsed design follow-up
# captured in the T.4b commit message):
#
# * L2-1 / L2-2 / L2-7 were specified for a clean SOTP form
#   `(D_axis & Ω_axis & I_g)` with one summand per spatial axis.  The
#   honest decomposition has 2 bespoke per-direction-sign summands (the
#   WDD recurrence couples spatial + angular axes via sweep direction;
#   not a clean tensor-product factor — same MA-Q1 critique applies as
#   for curvilinear M_angular_redist).
# * Revised: L2-1 = M_spatial is `_MSpatialOperatorSum` (an OperatorSum
#   subclass) with 2 `_SpatialSweepDirection` summands (forward μ_x > 0,
#   backward μ_x < 0).  L2-2 = each summand is a bespoke LinearOperator
#   leaf, NOT a `TensorProductOperator`.  L2-7 = direction_sign values
#   are disjoint (+1 / -1).
# * T.4b scope: SLAB ONLY.  Curvilinear (sphere, cylinder) goes through
#   the legacy path until T.4c lands the bespoke AngularRedistributionOperator
#   leaf.  2-D Cartesian stays procedural per Q1 hybrid.

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
    """Build per-(g, i, j) σ_t from sn_mesh.mat_map (matches the T.4a script)."""
    ng = sn_mesh.ng
    nx, ny = sn_mesh.nx, sn_mesh.ny
    sig_t = np.empty((ng, nx, ny), dtype=float)
    mat_map = sn_mesh.mat_map
    for ix in range(nx):
        for iy in range(ny):
            mid = int(mat_map[ix, iy])
            mat = sn_mesh.materials[mid]
            for g in range(ng):
                sig_t[g, ix, iy] = float(mat.SigT[g])
    return sig_t


class TestT4bMSpatialStructure:
    """L2-1 / L2-3 / L2-5 / L2-6 / L2-7 — type-level structural tests.

    Per the T.4 verification spec §4, with deviation noted above
    (per-direction summands, NOT 3-factor TP).
    """

    def test_M_spatial_is_operator_sum_with_two_per_direction_summands_slab(self):
        """L2-1 (revised) — slab M_spatial is `_MSpatialOperatorSum`
        with 2 per-direction `_SpatialSweepDirection` summands.
        """
        from orpheus.sn.operator import (
            _MSpatialOperatorSum,
            _SpatialSweepDirection,
        )
        from orpheus.numerics.operator import OperatorSum

        sn_mesh = _slab_mesh()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)

        M_spat = L.M_spatial
        assert isinstance(M_spat, _MSpatialOperatorSum)
        assert isinstance(M_spat, OperatorSum)
        # Default `OperatorSum` is binary — `a` is forward, `b` is
        # backward.  Wave O will read these as the per-direction
        # algebra exposure.
        assert isinstance(M_spat.a, _SpatialSweepDirection)
        assert isinstance(M_spat.b, _SpatialSweepDirection)

    def test_M_spatial_summand_direction_signs_disjoint(self):
        """L2-7 (revised) — the two per-direction summands have
        disjoint direction_sign values (one +1, one -1).
        """
        from orpheus.sn.operator import _SpatialSweepDirection

        sn_mesh = _slab_mesh()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        M_spat = L.M_spatial

        signs = {M_spat.a.direction_sign, M_spat.b.direction_sign}
        assert signs == {+1, -1}

    def test_M_angular_redist_slab_is_zero_operator(self):
        """L2-3 — slab `M_angular_redist` is `ZeroOperator`.

        Slab has no curvilinear angular redistribution (the M-M
        half-grid degenerates; ``IdentityAngularClosure`` returns zero
        contributions per the cell-balance algebra).
        """
        from orpheus.numerics.operator import ZeroOperator

        sn_mesh = _slab_mesh()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        assert isinstance(L.M_angular_redist, ZeroOperator)

    def test_M_angular_redist_curvilinear_is_bespoke_leaf(self):
        """L2-4 (T.4c) — curvilinear ``M_angular_redist`` is the bespoke
        :class:`AngularRedistributionOperator` leaf.

        Per the spec Q2 = (iii) bespoke leaf — the M-M half-grid
        recurrence is sequentially coupled along the angular axis;
        wrapping it as a 3-factor TP would false-assert separability
        the recurrence doesn't support.  The leaf carries
        ``{CAP_APPLY}`` and standalone applies per-cell M-M
        ``cell_contribution`` algebra.
        """
        from orpheus.sn.operator import AngularRedistributionOperator

        for builder in (_spherical_mesh, _cylindrical_mesh):
            sn_mesh = builder()
            sig_t = _sig_t_uniform(sn_mesh)
            L = StreamingOperator(sn_mesh, sig_t)
            M_ang = L.M_angular_redist
            assert isinstance(M_ang, AngularRedistributionOperator)
            assert M_ang.capabilities == frozenset({CAP_APPLY})
            # The leaf carries the same sn_mesh and sigma_t as the
            # parent operator (no copy; shared reference).
            assert M_ang.sn_mesh is sn_mesh
            assert M_ang.sigma_t is sig_t

    def test_M_spatial_capabilities_apply_only(self):
        """L2-6 — M_spatial advertises `{CAP_APPLY}` only.

        OperatorSum's capability closure: apply propagates iff both
        operands have apply (they do); solve does NOT propagate; the
        transpose closure is unset because per-direction leaves carry
        `{CAP_APPLY}` only (no `apply_transpose` yet — adjoint Krylov
        is plan §4 item 4, deferred).
        """
        sn_mesh = _slab_mesh()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        assert L.M_spatial.capabilities == frozenset({CAP_APPLY})

    def test_M_spatial_cached_property_returns_same_instance(self):
        """The summand identities are stable across re-access — a
        :func:`cached_property` invariant.  Important for Wave O
        / adjoint sites that may hold references to summands across
        invocations.
        """
        sn_mesh = _slab_mesh()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        assert L.M_spatial is L.M_spatial


class TestT4bAlgebraDecompositionInvariantSlab:
    """A-1 — slab algebra-decomposition invariant.

    `L.apply(ψ).bulk + σ_t·ψ.bulk == M_spatial.apply(ψ).bulk + M_angular_redist.apply(ψ).bulk`

    For slab, M_angular_redist = ZeroOperator, so the RHS reduces to
    M_spatial.apply(ψ).bulk alone.  This pins the Q3 (γ) subtraction
    contract: σ_t·ψ subtraction at apply boundary; algebra unwinds
    bit-exact.
    """

    def test_slab_apply_decomposition_invariant(self):
        """L = M_spatial + M_angular_redist − σ_t·ψ — bit-exact for slab."""
        sn_mesh = _slab_for_snapshot_arm(
            ng=2, bc_left=BC("vacuum"), bc_right=BC("vacuum"),
        )
        sig_t = _sigma_t_from_mat_map(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)

        rng = np.random.default_rng(20260531 + 41)
        state = sn_mesh.zeros_timed_full_field()
        from dataclasses import replace
        state = replace(
            state,
            bulk=replace(
                state.bulk,
                values=rng.uniform(0.05, 1.0, size=state.bulk.values.shape),
            ),
        )

        L_out = L.apply(state)
        M_spat_out = L.M_spatial.apply(state)
        # M_angular_redist.apply on slab returns the zero operator's
        # output; the algebra invariant unwinds via M_spatial alone.

        expected_bulk = (
            M_spat_out.bulk.values
            - sig_t[None, :, :, :] * state.bulk.values
        )
        np.testing.assert_array_equal(L_out.bulk.values, expected_bulk)
        # Boundary residual is carried by M_spatial alone (M_angular_redist
        # does not write to boundary; Resolution A subtraction is bulk-only).
        np.testing.assert_array_equal(
            L_out.boundary.values, M_spat_out.boundary.values,
        )


class TestT4bPreT4RegressionSnapshot:
    """L4-1 / L4-5 / L4-6 / L4-7 — bit-identity vs pre-T.4 snapshots.

    Per the T.4 verification spec §6 R3 (reduction-order drift) — the
    slab path's M_spatial.apply delegates to `_transport_operator_matvec_unified`
    bit-exact (no reduction reorder), so Route A strict bit-identity
    applies.

    Curvilinear (sphere, cylinder) arms are also pinned here because
    T.4b leaves curvilinear ON the legacy path — these snapshots
    must stay bit-identical.  T.4c will revisit them.
    """

    @pytest.fixture(scope="class")
    def snapshots(self):
        """Load the pre-T.4 snapshot bundle."""
        import os
        path = os.path.join(
            os.path.dirname(__file__), "_fixtures", "wave_t_t4",
            "pre_t4_snapshots.npz",
        )
        with np.load(path) as data:
            return {k: data[k] for k in data.files}

    def _capture_arm(self, sn_mesh: SNMesh, seed: int) -> tuple[np.ndarray, np.ndarray]:
        """Re-run StreamingOperator.apply on the snapshot fixture; return
        (bulk, boundary) values.
        """
        from dataclasses import replace
        sig_t = _sigma_t_from_mat_map(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        state = sn_mesh.zeros_timed_full_field()
        rng = np.random.default_rng(seed)
        state = replace(
            state,
            bulk=replace(
                state.bulk,
                values=rng.uniform(0.05, 1.0, size=state.bulk.values.shape),
            ),
        )
        out = L.apply(state)
        return out.bulk.values.copy(), out.boundary.values.copy()

    def test_slab_1g_vacuum_apply_bit_identical(self, snapshots):
        """L4-1 — slab 1G vacuum, P0 fixture."""
        bulk, boundary = self._capture_arm(
            _slab_for_snapshot_arm(
                ng=1, bc_left=BC("vacuum"), bc_right=BC("vacuum"),
            ),
            seed=20260531 + 1,
        )
        np.testing.assert_array_equal(
            bulk, snapshots["slab_1g_vacuum_apply_bulk"],
        )
        np.testing.assert_array_equal(
            boundary, snapshots["slab_1g_vacuum_apply_boundary"],
        )

    def test_slab_2g_vacuum_apply_bit_identical(self, snapshots):
        """L4-1 — slab 2G vacuum, P1 asymmetric SigS (ERR-002 detector)."""
        bulk, boundary = self._capture_arm(
            _slab_for_snapshot_arm(
                ng=2, bc_left=BC("vacuum"), bc_right=BC("vacuum"),
            ),
            seed=20260531 + 2,
        )
        np.testing.assert_array_equal(
            bulk, snapshots["slab_2g_vacuum_apply_bulk"],
        )
        np.testing.assert_array_equal(
            boundary, snapshots["slab_2g_vacuum_apply_boundary"],
        )

    def test_slab_2g_reflective_apply_bit_identical(self, snapshots):
        """L4-1 — slab 2G with reflective inner BC + vacuum outer.

        Catches BC convention drift via the D-B+1 specular permutation
        path (the D-B+1 production tensor-network instance).
        """
        bulk, boundary = self._capture_arm(
            _slab_for_snapshot_arm(
                ng=2, bc_left=BC("reflective"), bc_right=BC("vacuum"),
            ),
            seed=20260531 + 3,
        )
        np.testing.assert_array_equal(
            bulk, snapshots["slab_2g_reflective_apply_bulk"],
        )
        np.testing.assert_array_equal(
            boundary, snapshots["slab_2g_reflective_apply_boundary"],
        )


class TestT4cAlgebraDecompositionInvariantCurvilinear:
    """A-2 / A-3 — curvilinear algebra-decomposition invariant.

    `(L+C).apply(ψ) == M_spatial.apply(ψ) + M_angular_redist.apply(ψ)`
    at principled-equivalence ULP precision (per `vv-principles` §
    "Bit-identity vs principled-equivalence" three-criteria gate —
    M_spatial is computed via subtraction `(L+C) - M_ang` in
    `_MSpatialOperatorSum.apply` for curvilinear, so the equivalence
    differs from bit-identity at FP-non-associativity ULP).
    """

    def _build_curvilinear_fixture(self, builder, seed: int):
        """Construct a curvilinear sn_mesh + sig_t + state fixture."""
        from dataclasses import replace
        sn_mesh = builder(ng=2)
        sig_t = _sigma_t_from_mat_map(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        state = sn_mesh.zeros_timed_full_field()
        rng = np.random.default_rng(seed)
        state = replace(
            state,
            bulk=replace(
                state.bulk,
                values=rng.uniform(0.05, 1.0, size=state.bulk.values.shape),
            ),
        )
        return L, state

    def test_sphere_LpC_equals_M_spatial_plus_M_angular_redist(self):
        """A-2 — sphere: `(L+C)·ψ == M_spat·ψ + M_ang·ψ` (ULP-clean)."""
        from tests.sn._fixtures.wave_t_t4._capture_pre_t4_snapshots import (
            _sphere_mesh,
        )
        L, state = self._build_curvilinear_fixture(_sphere_mesh, seed=20260531 + 42)

        unified = _LC_matvec(state, L.sigma_t)
        M_spat = L.M_spatial.apply(state)
        M_ang = L.M_angular_redist.apply(state)

        # `(L+C) == M_spatial + M_angular_redist` at principled-
        # equivalence ULP — the subtraction `M_spatial = (L+C) - M_ang`
        # inside `_MSpatialOperatorSum.apply` re-introduces FP-non-
        # associativity drift bounded by ~reduction_depth × ULP.
        np.testing.assert_allclose(
            M_spat.bulk.values + M_ang.bulk.values,
            unified.bulk.values,
            rtol=1e-14, atol=1e-14,
        )
        # Boundary: M_angular_redist writes zero → M_spatial.boundary
        # equals the unified boundary exactly (bit-identical).
        np.testing.assert_array_equal(
            M_spat.boundary.values, unified.boundary.values,
        )

    def test_cylinder_LpC_equals_M_spatial_plus_M_angular_redist(self):
        """A-3 — cylinder: same invariant on multi-level fixture."""
        from tests.sn._fixtures.wave_t_t4._capture_pre_t4_snapshots import (
            _cylinder_mesh,
        )
        L, state = self._build_curvilinear_fixture(_cylinder_mesh, seed=20260531 + 43)

        unified = _LC_matvec(state, L.sigma_t)
        M_spat = L.M_spatial.apply(state)
        M_ang = L.M_angular_redist.apply(state)

        np.testing.assert_allclose(
            M_spat.bulk.values + M_ang.bulk.values,
            unified.bulk.values,
            rtol=1e-14, atol=1e-14,
        )
        np.testing.assert_array_equal(
            M_spat.boundary.values, unified.boundary.values,
        )

    def test_curvilinear_M_angular_redist_writes_no_boundary(self):
        """MA-Q4 — `M_angular_redist.apply` writes ZERO to boundary
        for sphere AND cylinder.  The angular redistribution is an
        interior-cell operation per the cell-balance algebra; spatial
        face residuals belong to `M_spatial` alone (Wave O typing
        consequence: M_angular_redist is a BulkOperator, M_spatial
        is a FullOperator).
        """
        from tests.sn._fixtures.wave_t_t4._capture_pre_t4_snapshots import (
            _sphere_mesh, _cylinder_mesh,
        )
        for builder, seed in (
            (_sphere_mesh, 20260531 + 81),
            (_cylinder_mesh, 20260531 + 82),
        ):
            L, state = self._build_curvilinear_fixture(builder, seed=seed)
            M_ang = L.M_angular_redist.apply(state)
            np.testing.assert_array_equal(
                M_ang.boundary.values, np.zeros_like(M_ang.boundary.values),
            )

    def test_curvilinear_M_angular_redist_linearity(self):
        """A-8 — linearity of `M_angular_redist.apply` for curvilinear.

        `M_ang(α·ψ + β·φ) = α·M_ang(ψ) + β·M_ang(φ)`.  Verified at
        FP precision because M-M's `cell_contribution` is linear in
        `psi_state`, and `precompute_psi_state` is linear in ψ
        (Carlson coupled-pole seed is linear; M-M closure is
        ``is_linear = True``).
        """
        from dataclasses import replace
        from tests.sn._fixtures.wave_t_t4._capture_pre_t4_snapshots import (
            _sphere_mesh,
        )

        sn_mesh = _sphere_mesh(ng=2)
        sig_t = _sigma_t_from_mat_map(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        rng = np.random.default_rng(20260531 + 91)

        # Two random ψ states.
        psi_shape = (sn_mesh.quad.N, sn_mesh.ng, sn_mesh.nx, sn_mesh.ny)
        psi_a = rng.uniform(0.05, 1.0, size=psi_shape)
        psi_b = rng.uniform(0.05, 1.0, size=psi_shape)
        alpha, beta = 0.37, -0.91

        def _make(values):
            state = sn_mesh.zeros_timed_full_field()
            return replace(state, bulk=replace(state.bulk, values=values))

        out_a = L.M_angular_redist.apply(_make(psi_a)).bulk.values
        out_b = L.M_angular_redist.apply(_make(psi_b)).bulk.values
        out_lin = L.M_angular_redist.apply(_make(alpha * psi_a + beta * psi_b)).bulk.values

        expected = alpha * out_a + beta * out_b
        np.testing.assert_allclose(out_lin, expected, rtol=1e-13, atol=1e-14)


class TestT4cPreT4RegressionSnapshotCurvilinear:
    """L4-2 / L4-3 — sphere + cylinder principled-equivalence vs pre-T.4 snapshot.

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
    """

    @pytest.fixture(scope="class")
    def snapshots(self):
        """Load the pre-T.4 snapshot bundle."""
        import os
        path = os.path.join(
            os.path.dirname(__file__), "_fixtures", "wave_t_t4",
            "pre_t4_snapshots.npz",
        )
        with np.load(path) as data:
            return {k: data[k] for k in data.files}

    def _capture_arm(self, sn_mesh: SNMesh, seed: int) -> tuple[np.ndarray, np.ndarray]:
        from dataclasses import replace
        sig_t = _sigma_t_from_mat_map(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        state = sn_mesh.zeros_timed_full_field()
        rng = np.random.default_rng(seed)
        state = replace(
            state,
            bulk=replace(
                state.bulk,
                values=rng.uniform(0.05, 1.0, size=state.bulk.values.shape),
            ),
        )
        out = L.apply(state)
        return out.bulk.values.copy(), out.boundary.values.copy()

    def test_sphere_1g_apply_bit_identical(self, snapshots):
        """L4-2 — sphere 1G P0 vacuum-at-r=R."""
        from tests.sn._fixtures.wave_t_t4._capture_pre_t4_snapshots import (
            _sphere_mesh,
        )
        bulk, boundary = self._capture_arm(_sphere_mesh(ng=1), seed=20260531 + 11)
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


class TestT5MaterializeInverseCache:
    """T.5 close-out — `M_spatial.materialize_inverse_cache(sigma_t)` API surface.

    Exposes the `CollisionCache` from M_spatial's natural angle
    (Pattern 2 dual-view of `CollisionCache.from_geometry`).  The
    method delegates to the existing factory; no new cache logic is
    introduced.  This test pins the API and verifies the cache fields
    match what `CollisionCache.from_geometry` would produce
    independently — same factory, same numerics.

    Future leverage (post-T.5 cache-unification micro-wave, NOT in
    T.5 scope): `_ensure_coll_cache` in `sn/sweep.py` could route
    through this method as the canonical cache-construction path,
    making M_spatial the single source of truth for its own inverse
    cache.  Today both pathways co-exist; this method exposes the
    operator-side angle so future refactoring has a name to migrate
    toward.
    """

    def test_M_spatial_materialize_inverse_cache_returns_collision_cache(self):
        """`M_spatial.materialize_inverse_cache()` returns a CollisionCache
        with the same fields as the canonical `from_geometry` factory."""
        from orpheus.sn.spatial.sweep_cache import (
            CollisionCache,
            GeometryCoefficients,
        )

        sn_mesh = _slab_mesh()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)

        # Build the cache via M_spatial.
        cache_via_op = L.M_spatial.materialize_inverse_cache()
        assert isinstance(cache_via_op, CollisionCache)

        # Compare to the canonical from_geometry path.
        geom = GeometryCoefficients.from_mesh_and_quad(sn_mesh)
        cache_canonical = CollisionCache.from_geometry(geom, sig_t[:, :, 0])

        # All three cache fields bit-identical (both paths call the
        # SAME `from_geometry` factory; this verifies the delegation
        # is wired correctly).
        np.testing.assert_array_equal(
            cache_via_op.inverse_denom, cache_canonical.inverse_denom,
        )
        np.testing.assert_array_equal(
            cache_via_op.a_attenuation, cache_canonical.a_attenuation,
        )
        np.testing.assert_array_equal(
            cache_via_op.cumprod_a, cache_canonical.cumprod_a,
        )


class TestT4dApply2DCartesianSourceHashPin:
    """A2D-1 — defensive source-hash pin on ``_apply_2d_cartesian``.

    Per the T.4 verification spec §1 Q1 = (c) HYBRID: T.4 lifts the 1-D
    matvec into the ``M_spatial + M_angular_redist`` decomposition, but
    LEAVES the 2-D Cartesian path procedural (cell-centre proxy + upwind
    FD via ``_compute_gradients``).  The decision to keep 2-D procedural
    rests on the architectural payload of the cell-centre-proxy ↔
    face-view-as-trace rewire (10% k_∞ drift comment at
    ``sn/operator.py:862-868``) and ``[[feedback_unify_after_two_instances]]``.

    Without a defensive pin, future author drift on the 2-D Cartesian
    path could silently change the behavior — the M_spatial /
    M_angular_redist properties don't touch it; the L1 MMS 2-D suite
    only catches a subset of failure modes.  This source-hash pin
    surfaces ANY intentional edit (the developer MUST update the pin)
    while preserving full freedom to refactor (the pin is one line).

    When updating ``_apply_2d_cartesian``:

    1. Run this test to see the new hash in the failure message.
    2. Update ``EXPECTED_SHA256`` with the new value.
    3. Commit the test + production change together.

    The hash invariant is INTENTIONALLY brittle — it's a tripwire, not
    a soft constraint.
    """

    # The SHA-256 hash of ``inspect.getsource(StreamingOperator._apply_2d_cartesian)``
    # as of T.4d (commit landing this test).  See class docstring for
    # the update procedure when the 2-D Cartesian path is intentionally
    # refactored.
    EXPECTED_SHA256: str = (
        "7a5dfee89499f246334c505330d86e4b30aeec8c2f33975c3306b862353d114e"
    )

    def test_apply_2d_cartesian_source_hash_unchanged(self):
        """Defensive: any source-level edit to ``_apply_2d_cartesian``
        MUST be deliberate — update the hash with the new value
        printed in the failure message.
        """
        import hashlib
        import inspect

        src = inspect.getsource(StreamingOperator._apply_2d_cartesian)
        actual = hashlib.sha256(src.encode("utf-8")).hexdigest()
        if actual != self.EXPECTED_SHA256:
            raise AssertionError(
                "A2D-1: `StreamingOperator._apply_2d_cartesian` source code "
                "has changed.  This is the T.4 hybrid-Q1 defensive pin — "
                "any intentional edit to the 2-D Cartesian path MUST "
                "update `TestT4dApply2DCartesianSourceHashPin.EXPECTED_SHA256` "
                "in this test file.\n"
                f"  EXPECTED: {self.EXPECTED_SHA256}\n"
                f"  ACTUAL:   {actual}\n"
                "If the change is intentional, copy the ACTUAL value above "
                "into EXPECTED_SHA256 and commit alongside the production "
                "change.  If the change is unintentional, revert it."
            )


class TestT4bMSpatialStandaloneApply:
    """Standalone per-direction `_SpatialSweepDirection.apply` slow-path
    invariants.  Used by future Wave-O / adjoint composition where the
    per-direction algebra is exercised in isolation.

    For T.4b: verify that the standalone summands SUM to the unified
    M_spatial.apply output bit-exact (algebraic-equivalence gate).
    """

    def test_sum_of_per_direction_standalone_applies_equals_unified(self):
        """`M_x_forward.apply(ψ) + M_x_backward.apply(ψ) == M_spatial.apply(ψ)`.

        The standalone per-direction apply is slow (each redoes the
        full bidirectional sweep then masks) but algebraically
        consistent with the orchestrated apply (the OperatorSum
        algebraic identity holds at value level).

        Verified on slab to keep the test fast; the per-direction
        decomposition's correctness is independent of geometry.
        """
        from dataclasses import replace

        sn_mesh = _slab_for_snapshot_arm(
            ng=2, bc_left=BC("vacuum"), bc_right=BC("vacuum"),
        )
        sig_t = _sigma_t_from_mat_map(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        state = sn_mesh.zeros_timed_full_field()
        rng = np.random.default_rng(20260531 + 73)
        state = replace(
            state,
            bulk=replace(
                state.bulk,
                values=rng.uniform(0.05, 1.0, size=state.bulk.values.shape),
            ),
        )

        M_spat = L.M_spatial
        forward_out = M_spat.a.apply(state)
        backward_out = M_spat.b.apply(state)
        sum_bulk = forward_out.bulk.values + backward_out.bulk.values
        sum_boundary = forward_out.boundary.values + backward_out.boundary.values

        unified_out = M_spat.apply(state)
        np.testing.assert_array_equal(sum_bulk, unified_out.bulk.values)
        np.testing.assert_array_equal(
            sum_boundary, unified_out.boundary.values,
        )

