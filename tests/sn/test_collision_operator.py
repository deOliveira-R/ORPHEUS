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

R-1 Step F (2026-05-19) — these tests migrated from the legacy
packed-1D contract (``rng.standard_normal(legacy.n_unknowns)``) to
the typed :class:`~orpheus.sn.angular_flux.AngularFlux` contract.
The legacy ``SNStreamingOperator.n_unknowns`` was sized to the
pre-B1'' eq_map (cell-only); the new operator-algebra contract uses
``build_equation_map_with_traces`` which sizes ``n_unknowns`` to
include the B1'' face block.  Operating on typed ``AngularFlux``
dissolves the size-mismatch concern: arithmetic, shape, and inverse
contracts all read directly off the typed surface.
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
from orpheus.sn.angular_flux import AngularFlux
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import CollisionOperator
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
    return SNMesh(mesh, quad, placeholder_materials(ng=2))


def _spherical_mesh(nx: int = 4, radius: float = 1.0) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = GaussLegendre1D.create(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=2))


def _cylindrical_mesh(nx: int = 4, radius: float = 1.0) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.01, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = LevelSymmetricSN.create(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=2))


def _random_angular_flux(
    sn_mesh: SNMesh, ng: int = 2, seed: int = 42,
) -> AngularFlux:
    """Random :class:`AngularFlux` shaped ``(N, ng, nx, ny)``.

    R-1 Step F migration target: replaces the legacy ``_packed_psi``
    helper that built bare-ndarray packed vectors of the
    ``SNStreamingOperator.n_unknowns`` size.  The typed surface ties
    the random field to the SN phase-space carrier directly — no
    eq_map dispatch required at the test level.
    """
    rng = np.random.default_rng(seed)
    N = sn_mesh.quad.N
    nx, ny = sn_mesh.nx, sn_mesh.ny
    return AngularFlux(
        rng.standard_normal((N, ng, nx, ny)), sn_mesh,
    )


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
        psi = _random_angular_flux(sn_mesh)
        out = C.apply(psi)
        assert isinstance(out, AngularFlux)
        assert out.values.shape == psi.values.shape

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_apply_zero_returns_zero(self, name, builder):
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = CollisionOperator(sn_mesh, sigma)
        N = sn_mesh.quad.N
        nx, ny = sn_mesh.nx, sn_mesh.ny
        zero = AngularFlux(np.zeros((N, sigma.shape[0], nx, ny)), sn_mesh)
        out = C.apply(zero)
        np.testing.assert_array_equal(out.values, 0.0)

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_apply_constant_sigma_uniform(self, name, builder):
        """With σ constant scalar c, apply(ψ) == c · ψ slot-by-slot."""
        sn_mesh = builder()
        nx, ny = sn_mesh.nx, sn_mesh.ny
        c = 0.4
        sigma = c * np.ones((2, nx, ny))  # PR-INDEX-3: (ng, nx, ny)
        C = CollisionOperator(sn_mesh, sigma)
        psi = _random_angular_flux(sn_mesh, seed=99)
        out = C.apply(psi)
        np.testing.assert_allclose(
            out.values, c * psi.values, rtol=1e-14, atol=1e-15,
        )

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_apply_is_linear(self, name, builder):
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = CollisionOperator(sn_mesh, sigma)
        psi1 = _random_angular_flux(sn_mesh, seed=51)
        psi2 = _random_angular_flux(sn_mesh, seed=52)
        alpha = 2.3
        beta = -0.7
        out_combined = C.apply(alpha * psi1 + beta * psi2)
        out_separate = alpha * C.apply(psi1) + beta * C.apply(psi2)
        np.testing.assert_allclose(
            out_combined.values, out_separate.values,
            rtol=1e-14, atol=1e-15,
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
        q = _random_angular_flux(sn_mesh, seed=88)
        out = C.solve(q)
        np.testing.assert_allclose(
            out.values, q.values / c, rtol=1e-14, atol=1e-15,
        )

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_solve_preserves_shape(self, name, builder):
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = CollisionOperator(sn_mesh, sigma)
        q = _random_angular_flux(sn_mesh)
        out = C.solve(q)
        assert isinstance(out, AngularFlux)
        assert out.values.shape == q.values.shape

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_solve_equals_division(self, name, builder):
        """solve(q) == q.values / sigma (broadcast over ordinates)."""
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = CollisionOperator(sn_mesh, sigma)
        q = _random_angular_flux(sn_mesh, seed=200)
        out = C.solve(q)
        # σ has shape (ng, nx, ny); broadcasts over the ordinate axis.
        np.testing.assert_allclose(
            out.values, q.values / sigma[None, :, :, :], rtol=1e-14,
        )


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
        q = _random_angular_flux(sn_mesh, seed=77)
        round_trip = C.apply(C.solve(q))
        np.testing.assert_allclose(
            round_trip.values, q.values, rtol=1e-12, atol=1e-14,
        )

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_solve_inverts_apply(self, name, builder):
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = CollisionOperator(sn_mesh, sigma)
        psi = _random_angular_flux(sn_mesh, seed=78)
        round_trip = C.solve(C.apply(psi))
        np.testing.assert_allclose(
            round_trip.values, psi.values, rtol=1e-12, atol=1e-14,
        )


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
        psi = _random_angular_flux(sn_mesh, seed=4242)
        out_apply = C.apply(psi)
        out_transpose = C.apply_transpose(psi)
        # Bit-exact: apply_transpose delegates to apply.
        np.testing.assert_array_equal(out_transpose.values, out_apply.values)


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
        psi = _random_angular_flux(sn_mesh, seed=10)
        out = C.apply(psi)
        assert np.any(np.abs(out.values) > 0)

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_works_with_removal_cross_section(self, name, builder):
        sn_mesh = builder()
        sigma_r = _sigma_removal(sn_mesh)
        C = CollisionOperator(sn_mesh, sigma_r)
        psi = _random_angular_flux(sn_mesh, seed=10)
        out = C.apply(psi)
        assert np.any(np.abs(out.values) > 0)

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
# 7. Sigma layout sanity — broadcast over the ordinate axis is correct.
# ═══════════════════════════════════════════════════════════════════════


class TestSigmaLayout:
    """The σ field broadcasts correctly across the ordinate axis.

    σ non-zero only at one cell → out is non-zero only at that
    cell's ``(ng, ix, iy)`` slot, replicated across every ordinate.
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
        psi = _random_angular_flux(sn_mesh, ng=ng, seed=33)
        out = C.apply(psi)
        # Build a mask shaped (nx, ny) selecting only the target cell.
        cell_mask = np.zeros((nx, ny), dtype=bool)
        cell_mask[ix_target, iy_target] = True
        # Output at NON-target cells is zero.
        np.testing.assert_array_equal(
            out.values[:, :, ~cell_mask], 0.0,
        )
        # Output at the TARGET cell equals 1.0 · psi at that cell
        # (sigma == 1.0 at the target).
        np.testing.assert_allclose(
            out.values[:, :, cell_mask],
            psi.values[:, :, cell_mask],
            rtol=1e-14, atol=1e-15,
        )
