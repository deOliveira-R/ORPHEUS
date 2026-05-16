"""Foundation tests for :class:`orpheus.sn.operator.CollisionOperator`.

Phase G Step 3+4.b.i (Issue #196). The "C" of the four-operator
algebra ``A_wg = L + C - S.foldable_part()``.

``CollisionOperator`` is the simplest leaf: diagonal in position,
group, and direction. ``apply`` is σ · ψ; ``solve`` is q / σ;
``apply_transpose`` equals ``apply`` (self-adjoint). All three
capabilities are analytic.

Convention-agnostic σ: the same operator class accepts either the
full ``σ_t`` (total cross-section) or the within-group removal
cross-section ``σ_r = σ_t − Σ_{s,0}^{g→g}``. No ``is_removal`` flag —
the operator's action is identical for either input.
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
from orpheus.sn.operator import CollisionOperator, SNStreamingOperator
from orpheus.sn.quadrature import GaussLegendre1D, LevelSymmetricSN
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation


# ═══════════════════════════════════════════════════════════════════════
# Fixtures (matching test_streaming_operator.py).
# ═══════════════════════════════════════════════════════════════════════


def _slab_mesh(nx: int = 4, length: float = 1.0) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, length, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = GaussLegendre1D.create(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials())


def _spherical_mesh(nx: int = 4, radius: float = 1.0) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = GaussLegendre1D.create(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials())


def _cylindrical_mesh(nx: int = 4, radius: float = 1.0) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.01, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = LevelSymmetricSN.create(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials())


def _packed_psi(sn_mesh: SNMesh, sigma: np.ndarray,
                seed: int = 42) -> np.ndarray:
    """Random packed-vector ψ matching the geometry's eq_map size."""
    legacy = SNStreamingOperator(sn_mesh, sigma)
    rng = np.random.default_rng(seed)
    return rng.standard_normal(legacy.n_unknowns)


def _sigma_total(sn_mesh: SNMesh, ng: int = 2) -> np.ndarray:
    """Random per-cell per-group cross-section, bounded away from 0.

    PR-INDEX-3: ``(ng, nx, ny)`` principled layout.
    """
    nx, ny = sn_mesh.nx, sn_mesh.ny
    rng = np.random.default_rng(seed=20260514)
    return 0.3 + 0.5 * rng.random((ng, nx, ny))


def _sigma_removal(sn_mesh: SNMesh, ng: int = 2) -> np.ndarray:
    """Synthetic σ_r — same shape, smaller magnitude. Handled identically."""
    nx, ny = sn_mesh.nx, sn_mesh.ny
    rng = np.random.default_rng(seed=20260515)
    return 0.1 + 0.3 * rng.random((ng, nx, ny))


GEOMETRIES = [
    ("slab", _slab_mesh),
    ("sphere", _spherical_mesh),
    ("cylinder", _cylindrical_mesh),
]


# ═══════════════════════════════════════════════════════════════════════
# 1. Capability advertising — all three capabilities are analytic.
# ═══════════════════════════════════════════════════════════════════════


class TestCapabilities:
    """CollisionOperator advertises apply + solve + apply_transpose."""

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_all_three_capabilities(self, name, builder):
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = CollisionOperator(sn_mesh, sigma)
        assert C.capabilities == frozenset(
            {CAP_APPLY, CAP_SOLVE, CAP_APPLY_TRANSPOSE}
        )

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_satisfies_linear_operator_protocol(self, name, builder):
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = CollisionOperator(sn_mesh, sigma)
        assert isinstance(C, LinearOperator)


# ═══════════════════════════════════════════════════════════════════════
# 2. apply — σ · ψ element-wise.
# ═══════════════════════════════════════════════════════════════════════


class TestApply:
    """apply preserves shape and computes σ · ψ slot-by-slot."""

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_apply_preserves_shape(self, name, builder):
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = CollisionOperator(sn_mesh, sigma)
        psi = _packed_psi(sn_mesh, sigma)
        out = C.apply(psi)
        assert out.shape == psi.shape

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_apply_zero_returns_zero(self, name, builder):
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = CollisionOperator(sn_mesh, sigma)
        psi = _packed_psi(sn_mesh, sigma)
        zero = np.zeros_like(psi)
        out = C.apply(zero)
        np.testing.assert_array_equal(out, np.zeros_like(psi))

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_apply_constant_sigma_uniform(self, name, builder):
        """With σ constant scalar c, apply(ψ) == c · ψ slot-by-slot."""
        sn_mesh = builder()
        nx, ny = sn_mesh.nx, sn_mesh.ny
        c = 0.4
        sigma = c * np.ones((2, nx, ny))  # PR-INDEX-3: (ng, nx, ny)
        C = CollisionOperator(sn_mesh, sigma)
        psi = _packed_psi(sn_mesh, sigma, seed=99)
        out = C.apply(psi)
        np.testing.assert_allclose(out, c * psi, rtol=1e-14, atol=1e-15)

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_apply_is_linear(self, name, builder):
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = CollisionOperator(sn_mesh, sigma)
        psi1 = _packed_psi(sn_mesh, sigma, seed=51)
        psi2 = _packed_psi(sn_mesh, sigma, seed=52)
        alpha = 2.3
        beta = -0.7
        out_combined = C.apply(alpha * psi1 + beta * psi2)
        out_separate = alpha * C.apply(psi1) + beta * C.apply(psi2)
        np.testing.assert_allclose(
            out_combined, out_separate, rtol=1e-14, atol=1e-15,
        )


# ═══════════════════════════════════════════════════════════════════════
# 3. solve — q / σ element-wise.
# ═══════════════════════════════════════════════════════════════════════


class TestSolve:
    """solve computes q / σ element-wise."""

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_solve_constant_sigma_uniform(self, name, builder):
        """With σ constant scalar c, solve(q) == q / c."""
        sn_mesh = builder()
        nx, ny = sn_mesh.nx, sn_mesh.ny
        c = 0.4
        sigma = c * np.ones((2, nx, ny))  # PR-INDEX-3: (ng, nx, ny)
        C = CollisionOperator(sn_mesh, sigma)
        q = _packed_psi(sn_mesh, sigma, seed=88)
        out = C.solve(q)
        np.testing.assert_allclose(out, q / c, rtol=1e-14, atol=1e-15)

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_solve_preserves_shape(self, name, builder):
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = CollisionOperator(sn_mesh, sigma)
        q = _packed_psi(sn_mesh, sigma)
        out = C.solve(q)
        assert out.shape == q.shape

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_solve_equals_division(self, name, builder):
        """solve(b) == b / sigma_packed at rtol=1e-14."""
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = CollisionOperator(sn_mesh, sigma)
        b = _packed_psi(sn_mesh, sigma, seed=200)
        out = C.solve(b)
        # Compute sigma_packed independently via the same gather pattern.
        from orpheus.sn.operator import (
            build_equation_map,
            build_equation_map_cylindrical,
            build_equation_map_spherical,
        )
        nx, ny = sn_mesh.nx, sn_mesh.ny
        quad = sn_mesh.quad
        curv = getattr(sn_mesh, "curvature", None)
        ng = int(sigma.shape[0])  # PR-INDEX-3: ng at axis 0
        if curv == "spherical":
            eq_map = build_equation_map_spherical(nx, quad, ng)
        elif curv == "cylindrical":
            eq_map = build_equation_map_cylindrical(nx, quad, ng)
        else:
            eq_map = build_equation_map(nx, ny, quad, ng)
        # PR-INDEX-3: principled (ng, nx, ny) — advanced index gives (ng, n_eq).
        sigma_packed = sigma[:, eq_map.ix, eq_map.iy].ravel(order='F')
        np.testing.assert_allclose(out, b / sigma_packed, rtol=1e-14)


# ═══════════════════════════════════════════════════════════════════════
# 4. apply ∘ solve = identity (mechanism criterion 6).
# ═══════════════════════════════════════════════════════════════════════


class TestApplySolveIdentity:
    """apply(solve(q)) ≈ q — the operator's inverse contract."""

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_apply_inverts_solve(self, name, builder):
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = CollisionOperator(sn_mesh, sigma)
        q = _packed_psi(sn_mesh, sigma, seed=77)
        round_trip = C.apply(C.solve(q))
        np.testing.assert_allclose(round_trip, q, rtol=1e-12, atol=1e-14)

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_solve_inverts_apply(self, name, builder):
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = CollisionOperator(sn_mesh, sigma)
        psi = _packed_psi(sn_mesh, sigma, seed=78)
        round_trip = C.solve(C.apply(psi))
        np.testing.assert_allclose(round_trip, psi, rtol=1e-12, atol=1e-14)


# ═══════════════════════════════════════════════════════════════════════
# 5. apply_transpose == apply (self-adjoint).
# ═══════════════════════════════════════════════════════════════════════


class TestApplyTranspose:
    """C is self-adjoint: C^T · ψ == C · ψ."""

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_apply_transpose_equals_apply(self, name, builder):
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = CollisionOperator(sn_mesh, sigma)
        psi = _packed_psi(sn_mesh, sigma, seed=4242)
        out_apply = C.apply(psi)
        out_transpose = C.apply_transpose(psi)
        # Bit-exact: apply_transpose delegates to apply.
        np.testing.assert_array_equal(out_transpose, out_apply)


# ═══════════════════════════════════════════════════════════════════════
# 6. Convention-agnostic σ — σ_t and σ_r are interchangeable.
# ═══════════════════════════════════════════════════════════════════════


class TestSigmaInterpretation:
    """The same CollisionOperator works for σ_t AND σ_r without API change.

    Pattern 4: the operator does not encode an interpretation of its
    σ — both are valid; the consumer decides which one to pass.
    """

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_works_with_total_cross_section(self, name, builder):
        sn_mesh = builder()
        sigma_t = _sigma_total(sn_mesh)
        C = CollisionOperator(sn_mesh, sigma_t)
        psi = _packed_psi(sn_mesh, sigma_t, seed=10)
        out = C.apply(psi)
        assert np.any(np.abs(out) > 0)

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_works_with_removal_cross_section(self, name, builder):
        sn_mesh = builder()
        sigma_r = _sigma_removal(sn_mesh)
        C = CollisionOperator(sn_mesh, sigma_r)
        psi = _packed_psi(sn_mesh, sigma_r, seed=10)
        out = C.apply(psi)
        assert np.any(np.abs(out) > 0)

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_no_is_removal_kwarg(self, name, builder):
        """CollisionOperator does not accept an is_removal kwarg.

        Anti-recommendation 8 from the brief: the type does not encode
        an interpretation it doesn't act on.
        """
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        with pytest.raises(TypeError):
            CollisionOperator(  # type: ignore[call-arg]
                sn_mesh, sigma=sigma, is_removal=False,
            )


# ═══════════════════════════════════════════════════════════════════════
# 7. Sigma layout sanity — packed gather matches per-cell read.
# ═══════════════════════════════════════════════════════════════════════


class TestSigmaLayout:
    """The packed-vector gather σ[ix(k), iy(k), g] is correct.

    σ non-zero only at one cell → out is non-zero only at that
    cell's packed slots.
    """

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_localised_sigma_localised_output(self, name, builder):
        sn_mesh = builder()
        ng = 2
        nx, ny = sn_mesh.nx, sn_mesh.ny
        ix_target = nx // 2
        iy_target = 0
        sigma = np.zeros((ng, nx, ny))  # PR-INDEX-3: (ng, nx, ny)
        sigma[:, ix_target, iy_target] = 1.0
        C = CollisionOperator(sn_mesh, sigma)
        sigma_ref = 0.5 * np.ones((ng, nx, ny))  # for sizing
        psi = _packed_psi(sn_mesh, sigma_ref, seed=33)
        out = C.apply(psi)
        from orpheus.sn.operator import (
            build_equation_map,
            build_equation_map_cylindrical,
            build_equation_map_spherical,
        )
        quad = sn_mesh.quad
        curv = getattr(sn_mesh, "curvature", None)
        if curv == "spherical":
            eq_map = build_equation_map_spherical(nx, quad, ng)
        elif curv == "cylindrical":
            eq_map = build_equation_map_cylindrical(nx, quad, ng)
        else:
            eq_map = build_equation_map(nx, ny, quad, ng)
        out_grid = out.reshape(ng, eq_map.n_eq, order='F')
        mask_target = (eq_map.ix == ix_target) & (eq_map.iy == iy_target)
        np.testing.assert_array_equal(
            out_grid[:, ~mask_target], 0.0,
        )
        psi_grid = psi.reshape(ng, eq_map.n_eq, order='F')
        np.testing.assert_allclose(
            out_grid[:, mask_target],
            1.0 * psi_grid[:, mask_target],
            rtol=1e-14, atol=1e-15,
        )
