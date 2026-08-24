"""Foundation tests for the collision multiplier :math:`M[\\sigma_t]`
(a :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`).

Phase G Step 3+4.b.i (Issue #196). The "C" of the four-operator
algebra ``A_wg = L + C - S.foldable_part()``. #261 retired the former
``CollisionOperator`` thin subclass — the collision instance is now a
plain ``MultiplicationOperator`` carrying ``σ_t`` as its coefficient.

The collision multiplier ``C = M[σ_t]`` is the simplest leaf: diagonal
in position, group, and direction. ``apply`` is σ · ψ; ``solve`` is
q / σ; ``apply_transpose`` equals ``apply`` (self-adjoint). Both
structural predicates (``is_invertible``, ``is_adjointable``) are
analytic — True for a positive σ.

Convention-agnostic σ: the same operator class accepts either the
full ``σ_t`` (total cross-section) or the within-group removal
cross-section ``σ_r = σ_t − Σ_{s,0}^{g→g}``. No ``is_removal`` flag —
the operator's action is identical for either input.

D-H.2-C1 (2026-05-28) — these tests migrated from the legacy
``orpheus.sn.angular_flux.AngularFlux`` input contract to the
composite :class:`~orpheus.transport.timed_full_field.TimedFullField`
carrier (bulk = L2 :class:`~orpheus.transport.fields.angular_flux.AngularFlux`).
The legacy class (a distinct pre-D-G design: bulk values + conflated
boundary buffer + history tuple) was DELETED in D-H.2-C5 phase 2
(commit ``d8843ba9``); rewriting fixtures in-place preserved the
breadth of coverage (per-σ-shape, per-geometry, linearity) while
exercising the composite branch of every operator method.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.operator import LinearOperator
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.source_sinks import AngularSourceSink
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import placeholder_materials
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux

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
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=2))


def _spherical_mesh(nx: int = 4, radius: float = 1.0) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=2))


def _cylindrical_mesh(nx: int = 4, radius: float = 1.0) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.01, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.folded_product(n_mu=4, n_phi=8)
    return SNMesh(mesh, quad, placeholder_materials(ng=2))


def _random_state(
    sn_mesh: SNMesh, ng: int = 2, seed: int = 42,
) -> TimedFullField:
    """Random :class:`TimedFullField` whose bulk has shape ``(N, ng, *spatial)``.

    D-H.2-C1: the composite carrier replaced the legacy
    ``orpheus.sn.angular_flux.AngularFlux`` input.  Bulk values
    are random; boundary is zero (the collision multiplier ``M[σ_t]`` is
    bulk-only — boundary is structurally implicit-zero per Option β3 /
    Issue #208).
    """
    rng = np.random.default_rng(seed)
    N = sn_mesh.quad.N
    state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn_mesh.full_field_space)
    return replace(
        state,
        interior=replace(
            state.interior, values=rng.standard_normal((N, ng, *sn_mesh.spatial_shape)),
        ),
    )


def _sigma_total(sn_mesh: SNMesh, ng: int = 2) -> np.ndarray:
    """Random per-cell per-group cross-section, bounded away from 0.

    PR-INDEX-3: ``(ng, *spatial)`` principled layout.
    """
    rng = np.random.default_rng(seed=20260514)
    return 0.3 + 0.5 * rng.random((ng, *sn_mesh.spatial_shape))


def _sigma_removal(sn_mesh: SNMesh, ng: int = 2) -> np.ndarray:
    """Synthetic σ_r — same shape, smaller magnitude. Handled identically."""
    rng = np.random.default_rng(seed=20260515)
    return 0.1 + 0.3 * rng.random((ng, *sn_mesh.spatial_shape))


GEOMETRIES = [
    ("slab", _slab_mesh),
    ("sphere", _spherical_mesh),
    ("cylinder", _cylindrical_mesh),
]


# ═══════════════════════════════════════════════════════════════════════
# 1. Predicate advertising — both structural axes are analytic.
# ═══════════════════════════════════════════════════════════════════════


class TestPredicates:
    """The collision multiplier ``M[σ_t]`` is invertible AND adjointable (σ > 0)."""

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_invertible_and_adjointable(self, name, builder):
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = MultiplicationOperator.from_mesh(sigma, sn_mesh)
        assert C.is_invertible and C.is_adjointable

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_satisfies_linear_operator_protocol(self, name, builder):
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = MultiplicationOperator.from_mesh(sigma, sn_mesh)
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
        C = MultiplicationOperator.from_mesh(sigma, sn_mesh)
        psi = _random_state(sn_mesh)
        out = C.apply(psi)
        assert isinstance(out, FullField)  # #257 S8a: timeless codomain (base arrow)
        assert isinstance(out.interior, AngularSourceSink)
        assert out.interior.values.shape == psi.interior.values.shape

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_apply_zero_returns_zero(self, name, builder):
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = MultiplicationOperator.from_mesh(sigma, sn_mesh)
        zero = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn_mesh.full_field_space)
        out = C.apply(zero)
        np.testing.assert_array_equal(out.interior.values, 0.0)

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_apply_constant_sigma_uniform(self, name, builder):
        """With σ constant scalar c, apply(ψ) == c · ψ slot-by-slot."""
        sn_mesh = builder()
        c = 0.4
        sigma = c * np.ones((2, *sn_mesh.spatial_shape))  # PR-INDEX-3: (ng, *spatial)
        C = MultiplicationOperator.from_mesh(sigma, sn_mesh)
        psi = _random_state(sn_mesh, seed=99)
        out = C.apply(psi)
        np.testing.assert_allclose(
            out.interior.values, c * psi.interior.values, rtol=1e-14, atol=1e-15,
        )

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_apply_is_linear(self, name, builder):
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = MultiplicationOperator.from_mesh(sigma, sn_mesh)
        psi1 = _random_state(sn_mesh, seed=51)
        psi2 = _random_state(sn_mesh, seed=52)
        # Linearity, stated directly (campaign 1 CS3 — flux lives in V, so
        # ψ₁ + ψ₂ is legal): homogeneity op(c·ψ) = c·op(ψ) AND additivity
        # op(ψ₁+ψ₂) = op(ψ₁)+op(ψ₂). Additivity alone reds an affine op
        # (the retired blend spelling op(ψ₁+λ(ψ₂−ψ₁)) could not — affine
        # maps preserve affine combinations; see the sharpness argument in
        # tests/sn/operators/test_declared_law_is_linear.py).
        c = 2.3
        np.testing.assert_allclose(
            C.apply(c * psi1).interior.values, (c * C.apply(psi1)).interior.values,
            rtol=1e-14, atol=1e-15,
        )
        out_sum = C.apply(psi1 + psi2)
        out_separate = C.apply(psi1) + C.apply(psi2)
        np.testing.assert_allclose(
            out_sum.interior.values, out_separate.interior.values,
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
        c = 0.4
        sigma = c * np.ones((2, *sn_mesh.spatial_shape))  # PR-INDEX-3: (ng, *spatial)
        C = MultiplicationOperator.from_mesh(sigma, sn_mesh)
        q = _random_state(sn_mesh, seed=88)
        out = C.solve(q)
        np.testing.assert_allclose(
            out.interior.values, q.interior.values / c, rtol=1e-14, atol=1e-15,
        )

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_solve_preserves_shape(self, name, builder):
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = MultiplicationOperator.from_mesh(sigma, sn_mesh)
        q = _random_state(sn_mesh)
        out = C.solve(q)
        assert isinstance(out, FullField)  # #257 S8a: timeless codomain (base arrow)
        assert isinstance(out.interior, AngularFlux)
        assert out.interior.values.shape == q.interior.values.shape

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_solve_equals_division(self, name, builder):
        """solve(q) == q.interior.values / sigma (broadcast over ordinates)."""
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = MultiplicationOperator.from_mesh(sigma, sn_mesh)
        q = _random_state(sn_mesh, seed=200)
        out = C.solve(q)
        # σ has shape (ng, *spatial); broadcasts over the ordinate axis.
        np.testing.assert_allclose(
            out.interior.values, q.interior.values / sigma[None], rtol=1e-14,
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
        C = MultiplicationOperator.from_mesh(sigma, sn_mesh)
        q = _random_state(sn_mesh, seed=77)
        round_trip = C.apply(C.solve(q))
        np.testing.assert_allclose(
            round_trip.interior.values, q.interior.values, rtol=1e-12, atol=1e-14,
        )

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_solve_inverts_apply(self, name, builder):
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = MultiplicationOperator.from_mesh(sigma, sn_mesh)
        psi = _random_state(sn_mesh, seed=78)
        round_trip = C.solve(C.apply(psi))
        np.testing.assert_allclose(
            round_trip.interior.values, psi.interior.values, rtol=1e-12, atol=1e-14,
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
        C = MultiplicationOperator.from_mesh(sigma, sn_mesh)
        psi = _random_state(sn_mesh, seed=4242)
        out_apply = C.apply(psi)
        out_transpose = C.apply_transpose(psi)
        # Bit-exact: apply_transpose delegates to apply.
        np.testing.assert_array_equal(
            out_transpose.interior.values, out_apply.interior.values,
        )


# ═══════════════════════════════════════════════════════════════════════
# 6. Convention-agnostic σ — σ_t and σ_r are interchangeable.
# ═══════════════════════════════════════════════════════════════════════


class TestSigmaInterpretation:
    """The same collision multiplier ``M[σ_t]`` works for σ_t AND σ_r without API change.

    Pattern 4: the operator does not encode an interpretation of its
    σ — both are valid; the consumer decides which one to pass.
    """

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_works_with_total_cross_section(self, name, builder):
        sn_mesh = builder()
        sigma_t = _sigma_total(sn_mesh)
        C = MultiplicationOperator.from_mesh(sigma_t, sn_mesh)
        psi = _random_state(sn_mesh, seed=10)
        out = C.apply(psi)
        assert np.any(np.abs(out.interior.values) > 0)

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_works_with_removal_cross_section(self, name, builder):
        sn_mesh = builder()
        sigma_r = _sigma_removal(sn_mesh)
        C = MultiplicationOperator.from_mesh(sigma_r, sn_mesh)
        psi = _random_state(sn_mesh, seed=10)
        out = C.apply(psi)
        assert np.any(np.abs(out.interior.values) > 0)

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
        nx = sn_mesh.nx
        ix_target = nx // 2
        sigma = np.zeros((ng, *sn_mesh.spatial_shape))  # PR-INDEX-3: (ng, *spatial)
        sigma[:, ix_target] = 1.0
        C = MultiplicationOperator.from_mesh(sigma, sn_mesh)
        psi = _random_state(sn_mesh, ng=ng, seed=33)
        out = C.apply(psi)
        # Build a mask shaped (*spatial,) selecting only the target cell.
        cell_mask = np.zeros(sn_mesh.spatial_shape, dtype=bool)
        cell_mask[ix_target] = True
        # Output at NON-target cells is zero.
        np.testing.assert_array_equal(
            out.interior.values[:, :, ~cell_mask], 0.0,
        )
        # Output at the TARGET cell equals 1.0 · psi at that cell
        # (sigma == 1.0 at the target).
        np.testing.assert_allclose(
            out.interior.values[:, :, cell_mask],
            psi.interior.values[:, :, cell_mask],
            rtol=1e-14, atol=1e-15,
        )


# ═══════════════════════════════════════════════════════════════════════
# 8. Composite-specific invariants — boundary + history_depth.
#
# Apply / solve / apply_transpose algebra is already covered above on
# the :class:`TimedFullField` carrier.  The tests below pin the
# composite-only contract: collision is bulk-only (Option β3 / Issue
# #208), so the output boundary is the implicit-zero
# :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`, and
# the history-depth metadata threads through unchanged.
# ═══════════════════════════════════════════════════════════════════════


class TestCompositeInvariants:
    """Composite-specific invariants beyond the bulk algebra above."""

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_apply_implicit_zero_boundary(self, name, builder):
        """Collision is bulk-only — boundary member is all zeros."""
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = MultiplicationOperator.from_mesh(sigma, sn_mesh)
        state = _random_state(sn_mesh, seed=172)

        out = C.apply(state)

        np.testing.assert_array_equal(out.boundary.values, 0.0)

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_solve_implicit_zero_boundary(self, name, builder):
        """Pseudoinverse leaves the rank-deficient face block zero."""
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = MultiplicationOperator.from_mesh(sigma, sn_mesh)
        state = _random_state(sn_mesh, seed=175)

        out = C.solve(state)

        np.testing.assert_array_equal(out.boundary.values, 0.0)

    @pytest.mark.parametrize("name,builder", GEOMETRIES)
    def test_output_is_timeless_full_field(self, name, builder):
        """#257 S8a — the matvec leaf is a base arrow ``FullField -> FullField``.

        The operator output is the TIMELESS
        :class:`~orpheus.transport.full_field.FullField` (history-free),
        regardless of the input iterate's ``history_depth``: the comonad
        lives on the iteration driver, not on the operator (was: the old
        convention stamped ``history_depth`` onto the output — re-pointed).
        """
        sn_mesh = builder()
        sigma = _sigma_total(sn_mesh)
        C = MultiplicationOperator.from_mesh(sigma, sn_mesh)
        for depth in (0, 1, 2, 4):
            state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn_mesh.full_field_space, history_depth=depth)
            out_apply = C.apply(state)
            out_solve = C.solve(state)
            # base-arrow codomain: a timeless FullField, NOT the timed subclass.
            # ``-O``-firing (Mode 8): explicit raises, not bare ``assert``.
            if type(out_apply) is not FullField or isinstance(out_apply, TimedFullField):
                pytest.fail(
                    f"{name} depth={depth}: C.apply output must be a timeless "
                    f"FullField, got {type(out_apply).__name__}"
                )
            if type(out_solve) is not FullField or isinstance(out_solve, TimedFullField):
                pytest.fail(
                    f"{name} depth={depth}: C.solve output must be a timeless "
                    f"FullField, got {type(out_solve).__name__}"
                )
