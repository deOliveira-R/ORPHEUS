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
from orpheus.sn.quadrature import GaussLegendre1D, LevelSymmetricSN

pytestmark = pytest.mark.foundation


# ═══════════════════════════════════════════════════════════════════════
# Fixtures — small slab / sphere / cylinder problems.
# ═══════════════════════════════════════════════════════════════════════


def _slab_mesh(nx: int = 4, length: float = 1.0) -> SNMesh:
    """Slab Mesh1D + GL N=4 quadrature, vacuum BCs."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, length, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = GaussLegendre1D.create(n_ordinates=4)
    return SNMesh(mesh, quad)


def _spherical_mesh(nx: int = 4, radius: float = 1.0) -> SNMesh:
    """Spherical Mesh1D + GL N=4, reflective inner / vacuum outer."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = GaussLegendre1D.create(n_ordinates=4)
    return SNMesh(mesh, quad)


def _cylindrical_mesh(nx: int = 4, radius: float = 1.0) -> SNMesh:
    """Cylindrical Mesh1D + Level-Symmetric SN-4 quadrature."""
    mesh = Mesh1D(
        edges=np.linspace(0.01, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = LevelSymmetricSN.create(sn_order=4)
    return SNMesh(mesh, quad)


def _sig_t_uniform(sn_mesh: SNMesh, ng: int = 2,
                   value: float = 0.5) -> np.ndarray:
    """σ_t uniform across cells / groups."""
    return value * np.ones((sn_mesh.nx, sn_mesh.ny, ng))


def _sig_t_heterogeneous(sn_mesh: SNMesh, ng: int = 2) -> np.ndarray:
    """σ_t heterogeneous — different value per (cell, group)."""
    nx, ny = sn_mesh.nx, sn_mesh.ny
    rng = np.random.default_rng(seed=20260514)
    return 0.3 + 0.5 * rng.random((nx, ny, ng))


def _packed_psi(sn_mesh: SNMesh, ng: int, seed: int = 42) -> np.ndarray:
    """Random packed-vector ψ matching the geometry's eq_map size."""
    sig_t = _sig_t_uniform(sn_mesh, ng=ng)
    legacy = SNStreamingOperator(sn_mesh, sig_t)
    rng = np.random.default_rng(seed)
    return rng.standard_normal(legacy.n_unknowns)


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

    The mechanism-criterion 3 contract. By Resolution A this holds
    **bit-exact** on ALL THREE geometries:

    .. math::

        L.{\rm apply}(\psi) + C.{\rm apply}(\psi)
        \;=\; [M(\psi;\sigma_t) - \sigma_t\odot\psi] + [\sigma_t\odot\psi]
        \;=\; M(\psi;\sigma_t)

    No xfail. The decomposition is algebraic, not an approximation.
    """

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_uniform_sigma_t_homogeneous(self, name, builder):
        """L(σ_t) + C(σ_t) ≡ legacy(σ_t) — homogeneous σ_t."""
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
    def test_heterogeneous_sigma_t(self, name, builder):
        """L(σ_t) + C(σ_t) ≡ legacy(σ_t) — heterogeneous σ_t per cell."""
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
    def test_via_operator_sum_algebra(self, name, builder):
        """``(L + C).apply(ψ)`` via OperatorSum distribution.

        Mechanism criterion 3 — the algebra composes correctly. The
        same numerical result as explicit ``L.apply + C.apply``,
        confirming :class:`OperatorSum.apply` routes through both
        leaves with no order-of-operation drift.
        """
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
# 6. The OperatorSum (L + C) advertises apply but NOT solve.
# ═══════════════════════════════════════════════════════════════════════


class TestSumCapabilities:
    """Pre-3+4.b.ii contract: (L + C).solve must NOT be available.

    The within-group fusion hook (substep 3+4.b.ii) will add a
    SN-aware solve at the OperatorSum level. Before that lands, the
    sum cannot invert itself — :class:`OperatorSum` does not
    propagate solve in general (Sherman-Morrison applies only under
    low-rank structure).
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
    def test_sum_does_not_advertise_solve(self, name, builder):
        sn_mesh = builder()
        sig_t = _sig_t_uniform(sn_mesh)
        L = StreamingOperator(sn_mesh, sig_t)
        C = CollisionOperator(sn_mesh, sig_t)
        A = L + C
        assert CAP_SOLVE not in A.capabilities
