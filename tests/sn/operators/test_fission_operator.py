"""Foundation tests for :class:`orpheus.transport.operators.fission.FissionOperator`.

Round 1.2 of Wave D of the SN reshape campaign (Issue #162). The
operator carries the same math :class:`SNSolver` used to expose under
``compute_fission_source``; that method is now a thin delegator. These
tests pin the lifted math at the operator level and the convention
that :meth:`FissionOperator.apply` returns :math:`F\\,\\phi` *without*
the :math:`1/k` division (the caller divides).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture, make_mixture
from orpheus.geometry import Mesh2D
from orpheus.numerics.operator import LinearOperator
from orpheus.transport.operators.fission import FissionOperator
from orpheus.transport.reaction_rate_functional import ReactionRateFunctional
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import SNSolver
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.timed_full_field import TimedFullField

pytestmark = pytest.mark.foundation


def _uniform_2d(nx, ny, delta, mat_map):
    return Mesh2D(
        edges_x=np.linspace(0, nx * delta, nx + 1),
        edges_y=np.linspace(0, ny * delta, ny + 1),
        mat_map=np.asarray(mat_map, dtype=int),
    )


@pytest.fixture
def solver_2g():
    fuel = get_mixture("A", "2g")
    mod = get_mixture("B", "2g")
    materials = {2: fuel, 0: mod}

    nx, ny = 6, 4
    delta = 0.2
    mat = np.zeros((nx, ny), dtype=int)
    mat[:3, :] = 2
    mat[3:, :] = 0

    mesh = _uniform_2d(nx, ny, delta, mat)
    quad = Quadrature.lebedev(order=17)
    sn_mesh = SNMesh(mesh, quad, materials)
    return SNSolver(sn_mesh)


# ──────────────────────────────────────────────────────────────────────
# Protocol contract
# ──────────────────────────────────────────────────────────────────────


class TestProtocolCompliance:
    """FissionOperator must satisfy the LinearOperator Protocol."""

    def test_implements_linear_operator(self, solver_2g):
        assert isinstance(solver_2g.fission_op, LinearOperator)

    def test_predicates_adjointable_not_invertible(self, solver_2g):
        """``is_adjointable`` True, ``is_invertible`` False — rank-1 in energy.

        F is non-invertible (rank-1 in energy has no useful inverse), but it
        HAS a transpose: the adjoint fission F† = |νΣf⟩⟨χ| (campaign #276),
        the χ↔νΣf dyad swap.
        """
        op = solver_2g.fission_op
        assert op.is_adjointable and not op.is_invertible


# ──────────────────────────────────────────────────────────────────────
# Bit-identical extraction (the load-bearing claim)
# ──────────────────────────────────────────────────────────────────────


class TestBitIdenticalExtraction:
    """FissionOperator math must match the legacy compute_fission_source."""

    def test_apply_matches_legacy_compute_at_k_eq_1(self, solver_2g):
        r"""apply(φ) = compute_fission_source(φ, k=1.0) within nulp=4.

        Wave T step T.2 relaxed this from ``np.array_equal`` to
        ``assert_array_almost_equal_nulp(nulp=4)`` per the
        ``vv-principles`` §"Bit-identity vs principled-equivalence"
        three-criteria gate:

        1. **Principled at every step** — the new path routes through
           the typed 2-factor :class:`TensorProductOperator` kernel
           (``outer(chi, ReactionRateFunctional(νΣf)) & IdentityOperator()``);
           every intermediate is a named physical quantity (per-cell
           production rate, emission spectrum × rate).
        2. **Structurally-independent reference** — the homogeneous-
           reflective :math:`k_\infty = \nu\Sigma_f / \Sigma_a`
           analytical limit is pinned by
           ``tests/sn/l1_analytical/test_kinf_homogeneous.py``.
        3. **FP-non-associativity, dimensionally explainable** —
           drift ≤ ``reduction_depth × ULP`` (here ng=2, drift ≈ 1
           ULP); reduction tree differs between ``np.einsum`` and
           ``(right * x).sum(axis=0)`` even though both compute the
           same algebraic sum.

        Pre-T.2 the apply body inlined the einsum and matched the
        legacy bit-identically.  Post-T.2 the inlined math lives in
        :class:`~orpheus.numerics.operator.RankOneOperator.apply`
        (single source of truth); the reduction order is whatever
        numpy chooses for ``(a * b).sum(axis=0)``.
        """
        np.random.seed(42)
        (nx, ny), ng = solver_2g.sn_mesh.spatial_shape, solver_2g.ng
        # PR-INDEX-4: FissionOperator.apply consumes / returns principled
        # (ng, nx, ny).  Use that shape directly.
        phi = np.random.rand(ng, nx, ny) + 0.1

        out_op = solver_2g.fission_op.apply(phi)
        # Reference: hand-coded version of the legacy method (no division by k).
        # All operands principled (ng, nx, ny).
        fission_rate = np.einsum("gxy,gxy->xy", solver_2g.mat_xs.fission_production, phi)
        expected = solver_2g.mat_xs.emission_spectrum * fission_rate[None, :, :]

        # Wave T step T.2: nulp=4 relaxation (see docstring).
        np.testing.assert_array_almost_equal_nulp(out_op, expected, nulp=4)

    def test_delegator_matches_apply_with_k(self, solver_2g):
        """SNSolver.compute_fission_source(φ, k) = fission_op.apply(φ) / k.

        Issue #196 PR-INDEX-5: both the delegator and the operator
        consume / return principled ``(ng, nx, ny)``.  No bridges.

        Wave T step T.2: the delegator and the operator both route
        through the same :attr:`FissionOperator.kernel` (single source
        of truth), so bit-identity holds between them — only the
        legacy-vs-RankOneOperator reduction-order drift (above test)
        requires nulp relaxation.
        """
        np.random.seed(7)
        (nx, ny), ng = solver_2g.sn_mesh.spatial_shape, solver_2g.ng
        phi = np.random.rand(ng, nx, ny) + 0.1

        for k in [1.0, 0.93, 1.27, 0.5]:
            out_via_delegator = solver_2g.compute_fission_source(phi, k)
            out_via_operator = solver_2g.fission_op.apply(phi) / k
            np.testing.assert_array_equal(out_via_delegator, out_via_operator)


# ──────────────────────────────────────────────────────────────────────
# Algebraic identities
# ──────────────────────────────────────────────────────────────────────


class TestRank1EnergyStructure:
    """F = χ ⊗ νΣ_f is rank-1 in energy per cell."""

    def test_zero_flux_zero_source(self, solver_2g):
        """φ = 0 => F·φ = 0 (linearity guard)."""
        (nx, ny), ng = solver_2g.sn_mesh.spatial_shape, solver_2g.ng
        phi = np.zeros((ng, nx, ny))
        out = solver_2g.fission_op.apply(phi)
        np.testing.assert_array_equal(out, np.zeros_like(phi))

    def test_apply_linearity(self, solver_2g):
        """F·(αφ_1 + βφ_2) = αF·φ_1 + βF·φ_2."""
        (nx, ny), ng = solver_2g.sn_mesh.spatial_shape, solver_2g.ng
        np.random.seed(13)
        phi1 = np.random.rand(ng, nx, ny) + 0.1
        phi2 = np.random.rand(ng, nx, ny) + 0.1
        alpha, beta = 2.5, -1.7

        lhs = solver_2g.fission_op.apply(alpha * phi1 + beta * phi2)
        rhs = (alpha * solver_2g.fission_op.apply(phi1)
               + beta * solver_2g.fission_op.apply(phi2))
        np.testing.assert_allclose(lhs, rhs, rtol=1e-12, atol=1e-13)

    def test_constant_flux_uniform_material(self):
        """In a uniform-material homogeneous problem, φ_g = c =>
        (Fφ)_g = χ_g · c · Σ_g' νΣ_{f,g'}."""
        mix = make_mixture(
            sig_t=np.array([0.5, 1.0]),
            sig_c=np.array([0.01, 0.02]),
            sig_f=np.array([0.01, 0.08]),
            nu=np.array([2.5, 2.5]),
            chi=np.array([1.0, 0.0]),
            sig_s=np.array([[0.38, 0.10], [0.00, 0.90]]),
        )
        nx, ny = 3, 3
        mesh = _uniform_2d(nx, ny, 0.5, np.zeros((nx, ny), dtype=int))
        quad = Quadrature.lebedev(order=17)
        solver = SNSolver(SNMesh(mesh, quad, {0: mix}))

        c = 1.5
        # PR-INDEX-4: principled (ng, nx, ny).
        phi = c * np.ones((solver.ng, nx, ny))
        out = solver.fission_op.apply(phi)

        # Hand-computed: per-cell fission rate = c · Σ_g νΣ_f[g] = c · Σ_g SigP[g].
        # (SigP is the production cross-section νΣ_f on the Mixture.)
        nuSigF_total = float(np.sum(mix.SigP))
        rate = c * nuSigF_total
        for ix in range(nx):
            for iy in range(ny):
                expected = mix.chi * rate
                np.testing.assert_allclose(
                    out[:, ix, iy], expected, rtol=1e-14,
                    err_msg=f"FissionOperator.apply mismatch at cell ({ix}, {iy})",
                )

    def test_chi_normalisation(self, solver_2g):
        """Σ_g χ_g = 1 per material (sanity check on the per-cell χ array)."""
        # Note: for materials with no fission, chi may be zero; we only check
        # the cells whose mixture has nonzero νΣ_f.
        nx, ny = solver_2g.sn_mesh.spatial_shape
        # PR-INDEX-3: solver.mat_xs.emission_spectrum / solver.mat_xs.fission_production are (ng, nx, ny).
        for mid, (ix_arr, iy_arr) in solver_2g.mat_xs.cells_by_material.items():
            for ix, iy in zip(ix_arr, iy_arr):
                chi_cell = solver_2g.mat_xs.emission_spectrum[:, ix, iy]
                if np.sum(solver_2g.mat_xs.fission_production[:, ix, iy]) > 1e-15:
                    np.testing.assert_allclose(chi_cell.sum(), 1.0, rtol=1e-12)


class TestKDivisionConvention:
    """The 1/k division is at the SOLVER level, not the operator level."""

    def test_apply_does_not_divide_by_k(self, solver_2g):
        """apply(φ) is independent of any eigenvalue — pure linear action."""
        np.random.seed(99)
        (nx, ny), ng = solver_2g.sn_mesh.spatial_shape, solver_2g.ng
        # PR-INDEX-4: principled (ng, nx, ny).
        phi = np.random.rand(ng, nx, ny) + 0.1

        # Construct the operator twice with no k handed in — apply
        # should not depend on any external state.
        out1 = solver_2g.fission_op.apply(phi)
        out2 = solver_2g.fission_op.apply(phi.copy())
        np.testing.assert_array_equal(out1, out2)

    def test_compute_fission_source_does_divide_by_k(self, solver_2g):
        """SNSolver.compute_fission_source(φ, k) = apply(φ) / k.

        Issue #196 PR-INDEX-5: ``phi`` principled ``(ng, nx, ny)``.
        """
        np.random.seed(101)
        (nx, ny), ng = solver_2g.sn_mesh.spatial_shape, solver_2g.ng
        phi = np.random.rand(ng, nx, ny) + 0.1

        out_k_one = solver_2g.compute_fission_source(phi, 1.0)
        out_k_two = solver_2g.compute_fission_source(phi, 2.0)
        np.testing.assert_allclose(out_k_one, 2.0 * out_k_two, rtol=1e-14)


# ──────────────────────────────────────────────────────────────────────
# D-H.2-C1 — Composite TimedFullField invariants (bulk-only fission).
#
# Per Option β3 (Issue #208), the FissionOperator has no boundary
# action — the emission spectrum lives only on cell-centred volumes.
# The composite branch returns a :class:`TimedFullField` with the
# fission per-ordinate source in the bulk member and an implicit-zero
# :class:`AngularBoundaryFlux` in the boundary member.  The parity test
# vs. legacy AngularFlux retired with D-H.2 (legacy branch deleted in
# C5); the algebra is exercised on the composite carrier directly.
# ──────────────────────────────────────────────────────────────────────


class TestCompositeInvariants:
    """Composite :class:`TimedFullField` variant: bulk-only fission."""

    def test_returns_timed_full_field(self, solver_2g):
        """Composite input → composite output (dispatch contract)."""
        from dataclasses import replace

        from orpheus.transport.fields.angular_flux import (
            AngularFlux,
        )
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.timed_full_field import TimedFullField

        sn_mesh = solver_2g.sn_mesh
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn_mesh.full_field_space)
        # Seed bulk with a deterministic non-zero per-ordinate ψ.
        np.random.seed(31)
        bulk_values = np.random.rand(*state.interior.values.shape) + 0.1
        state = replace(state, interior=replace(state.interior, values=bulk_values))

        out = solver_2g.fission_op.apply(state)

        # #257 S8a — the matvec leaf is a base arrow ``FullField -> FullField``,
        # so the output is the TIMELESS FullField (history-free).
        assert isinstance(out, FullField)
        assert not isinstance(out, TimedFullField)
        assert isinstance(out.interior, AngularSourceSink)
        assert out.interior.space is sn_mesh.angular_bulk_space

    def test_implicit_zero_boundary(self, solver_2g):
        """Fission has no boundary action — boundary member is all zeros."""
        from dataclasses import replace

        sn_mesh = solver_2g.sn_mesh
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn_mesh.full_field_space)
        np.random.seed(32)
        bulk_values = np.random.rand(*state.interior.values.shape) + 0.1
        state = replace(state, interior=replace(state.interior, values=bulk_values))

        out = solver_2g.fission_op.apply(state)

        # Implicit-zero boundary: every flat-buffer entry is exactly 0.
        # Option β3 / Wave O Issue #208 — the bulk-only nature is encoded in the
        # type: ``block_role = BlockRole.BULK`` (the :class:`BulkOperator` marker,
        # shipped in Wave O O.1).
        np.testing.assert_array_equal(out.boundary.values, 0.0)

    def test_zero_bulk_zero_output(self, solver_2g):
        """ψ = 0 ⇒ F·ψ = 0 (linearity guard at composite layer)."""
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, space=solver_2g.sn_mesh.full_field_space)
        out = solver_2g.fission_op.apply(state)
        np.testing.assert_array_equal(out.interior.values, 0.0)
        np.testing.assert_array_equal(out.boundary.values, 0.0)

    def test_output_is_timeless_full_field(self, solver_2g):
        """#257 S8a — the matvec leaf is a base arrow ``FullField -> FullField``.

        The output is the TIMELESS FullField (history-free) regardless of the
        input iterate's ``history_depth``: the comonad lives on the iteration
        driver, not the operator (was: the old convention stamped
        ``history_depth`` onto the output — re-pointed).
        """
        sn_mesh = solver_2g.sn_mesh
        for depth in (0, 1, 2, 4):
            state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn_mesh.full_field_space, history_depth=depth)
            out = solver_2g.fission_op.apply(state)
            assert isinstance(out, FullField)
            assert not isinstance(out, TimedFullField)


# ──────────────────────────────────────────────────────────────────────
# Wave T step T.2 — RankOneOperator TP kernel verification gate.
#
# T.2 lifts the inlined ``np.einsum("gxy,gxy->xy", sig_p, phi) *
# chi[None, :, :]`` body into the typed 2-factor TP kernel
# ``outer(chi, ReactionRateFunctional(νΣf)) & IdentityOperator()``.  The
# `kernel` property is the type-visible §16A.10 ``B = G_patch ⊗ K_omega
# ⊗ K_g`` decomposition for fission.  This class pins the kernel
# structure so future refactors don't silently undo the T.2 lift.
# ──────────────────────────────────────────────────────────────────────


class TestRankOneTensorProductKernel:
    """Wave T step T.2 — fission kernel is a 2-factor TensorProductOperator
    wrapping :class:`~orpheus.numerics.operator.RankOneOperator`.
    """

    def test_kernel_is_2_factor_tensor_product(self, solver_2g):
        """``fission_op.kernel`` returns ``RankOneOperator & IdentityOperator``.

        The 2-factor TP shape mirrors D-B+1's specular BC pattern and
        T.1's vacuum / white / albedo / periodic lifts: a leading
        operator that does the actual algebraic work on the
        non-trivial axis (here ``axis=0``, the group axis) plus an
        :class:`IdentityOperator` factor that advertises the spatial-
        axis broadcast.
        """
        from orpheus.numerics.operator import (
            IdentityOperator,
            RankOneOperator,
            TensorProductOperator,
        )

        kernel = solver_2g.fission_op.kernel
        assert isinstance(kernel, TensorProductOperator)
        assert len(kernel.ops) == 2
        assert isinstance(kernel.ops[0], RankOneOperator)
        assert isinstance(kernel.ops[1], IdentityOperator)

    def test_kernel_reconstruction_is_chi_functional_is_production_rate(self, solver_2g):
        """The RankOneOperator dyad: reconstruction = χ, functional = ⟨νΣ_f| (axis 0).

        Post-refactor the rank-1 op is the dyad ``outer(χ, ReactionRateFunctional(νΣf))``:
        its ``reconstruction`` column is the emission spectrum χ (bound by
        REFERENCE — ``mat_xs`` shares the same numpy buffer across calls, so
        ``is`` holds), and its ``functional`` row co-vector is the
        production-rate ``ReactionRateFunctional`` (the §5.6 contraction,
        ``axis=0`` over groups, ``weight == νΣ_f``).
        """
        kernel = solver_2g.fission_op.kernel
        rank_one = kernel.ops[0]
        # ``reconstruction`` binds the χ buffer by reference (not a copy).
        assert rank_one.reconstruction is solver_2g.mat_xs.emission_spectrum
        # The row co-vector is the production-rate reaction-rate functional.
        assert isinstance(rank_one.functional, ReactionRateFunctional)
        assert rank_one.functional.axis == 0
        np.testing.assert_array_equal(
            np.asarray(rank_one.functional.weight),
            solver_2g.mat_xs.fission_production,
        )

    def test_kernel_apply_matches_apply_dispatch(self, solver_2g):
        """``kernel.apply(phi.values)`` equals ``fission_op.apply(phi).values``
        for the bare-ndarray and ScalarFlux dispatch arms — both go
        through the same kernel.  Single source of truth for the
        rank-1 outer-product math.
        """
        np.random.seed(57)
        nx, ny, ng = (
            solver_2g.sn_mesh.nx,
            solver_2g.sn_mesh.spatial_shape[1],
            solver_2g.ng,
        )
        phi_arr = np.random.rand(ng, nx, ny) + 0.1

        out_via_apply = solver_2g.fission_op.apply(phi_arr)
        out_via_kernel = solver_2g.fission_op.kernel.apply(phi_arr)
        # Same code path — bit-identical.
        np.testing.assert_array_equal(out_via_apply, out_via_kernel)

    def test_kernel_predicates_adjointable_not_invertible(self, solver_2g):
        """Kernel is adjointable but NOT invertible — the TP closure laws.

        Predicate closure: :class:`RankOneOperator` (with the
        :class:`ReactionRateFunctional` row, an ``InnerProductFunctional``) is
        adjointable-not-invertible (campaign #276 — the dual dyad);
        :class:`IdentityOperator` is both. The TP conjunction keeps
        adjointability (both factors have it) and drops invertibility (only
        Identity has it) — the §15 rank-1 fission structure has no useful
        inverse, but it DOES transpose (F† = |νΣf⟩⟨χ|).
        """
        kernel = solver_2g.fission_op.kernel
        assert kernel.is_adjointable and not kernel.is_invertible
