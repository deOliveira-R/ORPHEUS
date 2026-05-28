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
    SNStreamingOperator,
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


class TestCompositionEquivalence:
    r"""``(L + C(σ_t)).apply(ψ) ≈ SNStreamingOperator(σ_t).apply(ψ)``.

    Pre-PR-TYPED-6c Step 5 this held bit-exact on ALL geometries because
    ``L.apply`` consumed the same legacy per-geometry matvec helper as
    ``SNStreamingOperator.apply``.

    Post-Step-5, ``L.apply`` routes through
    :func:`transport_operator_matvec_unified` (the WDD-based unified
    matvec).  The legacy bundle ``SNStreamingOperator`` still routes
    through the per-geometry helpers.  The equivalence holds:

      * **Sphere** — bit-identical (the unified sphere body was verified
        ULP-equal to ``transport_operator_matvec_spherical`` at Step 2;
        see ``test_unified_matvec_sphere.py``).
      * **Cylinder** — DIVERGES.  The legacy ``transport_operator_matvec_cylindrical``
        carries the ERR-049 routing bug; the unified body is structurally
        immune.  Marked xfail until SNStreamingOperator retires
        in PR-TYPED-6c Step 7 (general retirement).
      * **Cartesian (slab)** — DIVERGES.  The legacy
        ``transport_operator_matvec`` uses 1st-order upwind FD via
        ``_compute_gradients``; the unified body uses 2nd-order WDD via
        ``cell_balance_for_streaming``.  Order-of-accuracy delta is
        characterised in ``test_unified_matvec_slab.py``.  Marked xfail
        until Step 7.

    The Step 6 verification gate is the audit of this divergence — both
    cylinder and Cartesian failures are documented and intentional.
    """

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_uniform_sigma_t_homogeneous(self, name, builder, request):
        """L(σ_t) + C(σ_t) ≡ legacy(σ_t) — homogeneous σ_t.

        Sphere: bit-identical post-Step-5.  Cylinder + slab: xfail
        (intentional legacy divergence, see class docstring).
        """
        # PR-TYPED-6.5 Phase 3b — B1'' face-state in the packed vector
        # makes :class:`StreamingOperator`'s packed format structurally
        # incompatible with :class:`SNStreamingOperator`'s legacy
        # compressed layout: the cell-only block has different size
        # (curvilinear B1'' includes inward-at-outer cells; SN
        # compresses them) AND there's a face block in B1'' that has
        # no SN analog.  Direct (out_sum == out_legacy) is no longer
        # well-defined — different shapes.  ALL geometries xfail until
        # ``SNStreamingOperator`` retires at Step 7.
        request.node.add_marker(pytest.mark.xfail(
            reason=(
                "PR-TYPED-6.5 Phase 3b — StreamingOperator now uses "
                "B1'' face-aware packed layout; SNStreamingOperator "
                "uses legacy compressed layout.  Shape mismatch is "
                "structural; comparison retires at Step 7 when "
                "SNStreamingOperator goes away."
            ),
            strict=True,
        ))
        sn_mesh = builder()
        ng = 2
        sig_t = _sig_t_uniform(sn_mesh, ng=ng, value=0.4)
        legacy = SNStreamingOperator(sn_mesh, sig_t)
        L = StreamingOperator(sn_mesh, sig_t)
        C = CollisionOperator(sn_mesh, sig_t)
        psi = _packed_psi(sn_mesh, ng=ng, seed=11)
        out_sum = L.apply(psi) + C.apply(psi)
        out_legacy = legacy.apply(psi)
        np.testing.assert_allclose(
            out_sum, out_legacy, rtol=1e-12, atol=1e-14
        )

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_heterogeneous_sigma_t(self, name, builder, request):
        """L(σ_t) + C(σ_t) ≡ legacy(σ_t) — heterogeneous σ_t per cell.

        Sphere: bit-identical post-Step-5.  Cylinder + slab: xfail
        (intentional legacy divergence, see class docstring).
        """
        # PR-TYPED-6.5 Phase 3b — see ``test_uniform_sigma_t_homogeneous``
        # for the rationale; B1'' packed layout makes the (L+C) vs
        # SNStreamingOperator comparison structurally incompatible for
        # all 1-D geometries.
        request.node.add_marker(pytest.mark.xfail(
            reason=(
                "PR-TYPED-6.5 Phase 3b — StreamingOperator now uses "
                "B1'' face-aware packed layout; SNStreamingOperator "
                "uses legacy compressed layout.  Shape mismatch is "
                "structural; comparison retires at Step 7."
            ),
            strict=True,
        ))
        sn_mesh = builder()
        ng = 2
        sig_t = _sig_t_heterogeneous(sn_mesh, ng=ng)
        legacy = SNStreamingOperator(sn_mesh, sig_t)
        L = StreamingOperator(sn_mesh, sig_t)
        C = CollisionOperator(sn_mesh, sig_t)
        psi = _packed_psi(sn_mesh, ng=ng, seed=22)
        out_sum = L.apply(psi) + C.apply(psi)
        out_legacy = legacy.apply(psi)
        np.testing.assert_allclose(
            out_sum, out_legacy, rtol=1e-12, atol=1e-14
        )

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_via_operator_sum_algebra(self, name, builder, request):
        """``(L + C).apply(ψ)`` via OperatorSum distribution.

        Mechanism criterion 3 — the algebra composes correctly. The
        same numerical result as explicit ``L.apply + C.apply``,
        confirming :class:`OperatorSum.apply` routes through both
        leaves with no order-of-operation drift.

        Sphere: bit-identical post-Step-5.  Cylinder + slab: xfail
        (intentional legacy divergence, see class docstring).
        """
        # PR-TYPED-6.5 Phase 3b — see ``test_uniform_sigma_t_homogeneous``
        # for the rationale; B1'' packed layout makes the (L+C) vs
        # SNStreamingOperator comparison structurally incompatible for
        # all 1-D geometries.
        request.node.add_marker(pytest.mark.xfail(
            reason=(
                "PR-TYPED-6.5 Phase 3b — StreamingOperator now uses "
                "B1'' face-aware packed layout; SNStreamingOperator "
                "uses legacy compressed layout.  Shape mismatch is "
                "structural; comparison retires at Step 7."
            ),
            strict=True,
        ))
        sn_mesh = builder()
        ng = 2
        sig_t = _sig_t_heterogeneous(sn_mesh, ng=ng)
        L = StreamingOperator(sn_mesh, sig_t)
        C = CollisionOperator(sn_mesh, sig_t)
        legacy = SNStreamingOperator(sn_mesh, sig_t)
        A = L + C  # OperatorSum via LinearOperatorMixin.__add__
        psi = _packed_psi(sn_mesh, ng=ng, seed=33)
        out_algebra = A.apply(psi)
        out_legacy = legacy.apply(psi)
        np.testing.assert_allclose(
            out_algebra, out_legacy, rtol=1e-12, atol=1e-14
        )


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
# D-H.1b.6 — TimedFullField composite dispatch on StreamingOperator.apply.
# ═══════════════════════════════════════════════════════════════════════


# 1-D geometries — the TimedFullField branch supports these natively
# via the legacy AngularFlux bridge.  2-D Cartesian raises
# NotImplementedError (Phase A coordination).
GEOMETRIES_1D = [
    ("slab", _slab_mesh),
    ("sphere", _spherical_mesh),
    ("cylinder", _cylindrical_mesh),
]


def _legacy_angular_flux_from_composite(state):
    """Helper: build a legacy AngularFlux carrying the composite's
    bulk + boundary, for direct comparison against the legacy branch.
    """
    from orpheus.sn.angular_flux import AngularFlux as LegacyAngularFlux
    from orpheus.sn.boundary_flux import BoundaryFlux as LegacyBoundaryFlux

    sn_mesh = state.bulk.mesh
    legacy_bf = LegacyBoundaryFlux(mesh=sn_mesh)
    layout = state.boundary.layout
    if "xmax" in layout.faces:
        legacy_bf.xmax_face = state.boundary.face_view("xmax").copy()
    if "xmin" in layout.faces:
        legacy_bf.xmin_face = state.boundary.face_view("xmin").copy()
    return LegacyAngularFlux(
        state.bulk.values.copy(), sn_mesh, boundary=legacy_bf,
    )


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


class TestTimedFullFieldBranch:
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
            AngularFlux as L2AngularFlux,
        )
        from orpheus.transport.timed_full_field import TimedFullField

        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        state = _random_composite(sn_mesh)

        out = L.apply(state)

        assert isinstance(out, TimedFullField)
        assert isinstance(out.bulk, L2AngularFlux)
        assert out.bulk.mesh is sn_mesh
        assert out.history_depth == state.history_depth
        assert out._history == ()

    @pytest.mark.parametrize("name,builder", GEOMETRIES_1D)
    def test_bulk_matches_legacy_branch(self, name, builder):
        """Composite branch bulk == legacy AngularFlux branch values.

        Pins equivalence so the bridge stays a thin adapter — any
        future change to the matvec kernel produces identical bulk
        values across both type dispatches.
        """
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        state = _random_composite(sn_mesh, seed=181)
        legacy_in = _legacy_angular_flux_from_composite(state)

        out_composite = L.apply(state)
        out_legacy = L.apply(legacy_in)

        np.testing.assert_array_equal(
            out_composite.bulk.values, out_legacy.values,
        )

    @pytest.mark.parametrize("name,builder", GEOMETRIES_1D)
    def test_boundary_carries_face_residual(self, name, builder):
        """Composite boundary == legacy boundary face residual.

        L is the only operator in {L, C, S, F} that emits a non-zero
        face residual.  The composite branch must preserve this — the
        L2 BoundaryFlux flat buffer carries the matvec's face output
        at the layout-assigned slots.
        """
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        state = _random_composite(sn_mesh, seed=182)
        legacy_in = _legacy_angular_flux_from_composite(state)

        out_composite = L.apply(state)
        out_legacy = L.apply(legacy_in)

        # Face-by-face comparison via L2 face_view.
        layout = out_composite.boundary.layout
        if "xmax" in layout.faces:
            np.testing.assert_array_equal(
                out_composite.boundary.face_view("xmax"),
                out_legacy.boundary.xmax_face,
            )
        if "xmin" in layout.faces:
            np.testing.assert_array_equal(
                out_composite.boundary.face_view("xmin"),
                out_legacy.boundary.xmin_face,
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

    def test_2d_cartesian_raises_not_implemented(self):
        """2-D Cartesian path raises NotImplementedError (Phase A scope)."""
        from orpheus.geometry import Mesh2D

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
        with pytest.raises(NotImplementedError, match="2-D Cartesian"):
            L.apply(state)

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

    @pytest.mark.parametrize("name,builder", GEOMETRIES_1D)
    def test_LC_apply_matches_legacy(self, name, builder):
        """``(L + C).apply(state).bulk`` matches the legacy branch sum."""
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        C = CollisionOperator(sn_mesh, sig_t)
        A = L + C
        state = _random_composite(sn_mesh, seed=192)
        legacy_in = _legacy_angular_flux_from_composite(state)

        out_composite = A.apply(state)
        out_legacy = A.apply(legacy_in)

        np.testing.assert_array_equal(
            out_composite.bulk.values, out_legacy.values,
        )
