"""Foundation tests for the Issue #197 PR-TYPED-5 :class:`Solution` type.

Pins the structural contract of :class:`Solution` +
:class:`IterationHistory` + :class:`SolutionDiff`:

* shape consistency at construction (delegates to the typed fields);
* mesh-binding identity across angular_flux / scalar_flux / boundary_flux;
* eigenvalue vs fixed-source discrimination via
  :meth:`Solution.is_eigenvalue` / :meth:`Solution.is_fixed_source`;
* :meth:`IterationHistory.dominance_ratio` for a recorded trajectory;
* :meth:`Solution.compare` field-by-field summary;
* :meth:`Solution.reaction_rate_density` :math:`\\sigma \\cdot \\phi` math.

These are foundation tests — software invariants of the typed return
contract, not physics claims about a solver, so they carry
``@pytest.mark.foundation`` per the V&V harness convention.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn.geometry import SNMesh
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.timed_full_field import TimedFullField
from orpheus.sn.solution import IterationHistory, Solution, SolutionDiff

from tests.sn._test_helpers import placeholder_materials
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux

pytestmark = pytest.mark.foundation


# ── Fixtures ─────────────────────────────────────────────────────────


def _slab_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    """Build a small slab :class:`SNMesh` for unit testing."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _make_fluxes(sn_mesh: SNMesh, fill: float = 1.0):
    """Build (state, scalar, boundary) for the given mesh.

    D-H.1b — :class:`Solution.angular_flux` is now a
    :class:`~orpheus.transport.timed_full_field.TimedFullField`
    composite. The triple return-tuple is preserved for fixture-shape
    stability; the first slot is now a TimedFullField (bulk angular
    flux + boundary trace + history), the third slot is the
    composite's boundary trace.
    """
    state = TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh)
    state.bulk.values[:] = fill
    state.boundary.values[:] = fill
    phi = ScalarFlux.from_mesh(np.full((sn_mesh.ng, *sn_mesh.spatial_shape), fill), sn_mesh)
    return state, phi, state.boundary


# ════════════════════════════════════════════════════════════════════
# IterationHistory
# ════════════════════════════════════════════════════════════════════


class TestIterationHistory:
    def test_default_empty(self) -> None:
        """An IterationHistory with no entries has no dominance ratio."""
        h = IterationHistory()
        assert h.keff_history == ()
        assert h.flux_residuals == ()
        assert h.n_inner is None
        assert h.n_outer is None
        assert h.converged is True
        assert h.dominance_ratio() is None
        assert h.latest_keff() is None
        assert h.latest_residual() is None

    def test_dominance_ratio_single_entry(self) -> None:
        h = IterationHistory(keff_history=(1.234,))
        assert h.dominance_ratio() is None  # need ≥ 2

    def test_dominance_ratio_three_iters(self) -> None:
        h = IterationHistory(keff_history=(1.0, 1.1, 1.10005))
        # |1.10005 - 1.1| / |1.1| ≈ 4.545e-5
        ratio = h.dominance_ratio()
        assert ratio is not None
        np.testing.assert_allclose(ratio, 5e-5 / 1.1, rtol=1e-10)

    def test_dominance_ratio_zero_prev(self) -> None:
        """A zero penultimate keff returns None to avoid divide-by-zero."""
        h = IterationHistory(keff_history=(0.0, 1.0))
        assert h.dominance_ratio() is None

    def test_latest_keff(self) -> None:
        h = IterationHistory(keff_history=(1.0, 1.1, 1.05))
        assert h.latest_keff() == 1.05

    def test_latest_residual(self) -> None:
        h = IterationHistory(flux_residuals=(1e-3, 1e-6, 1e-9))
        assert h.latest_residual() == 1e-9

    def test_frozen(self) -> None:
        """IterationHistory is frozen — fields cannot be reassigned."""
        h = IterationHistory(keff_history=(1.0,))
        with pytest.raises((AttributeError, Exception)):
            h.keff_history = (1.0, 1.1)  # type: ignore[misc]

    def test_keff_history_is_tuple(self) -> None:
        """The trajectory MUST be a tuple — frozen-with-mutable-list anti-pattern."""
        h = IterationHistory(keff_history=(1.0, 1.1))
        assert isinstance(h.keff_history, tuple)
        # Pinning the contract: a list passed in survives but the field
        # is documented as a tuple; consumers MUST treat it as immutable.


# ════════════════════════════════════════════════════════════════════
# Solution — construction + mesh-identity contract
# ════════════════════════════════════════════════════════════════════


class TestSolutionConstruction:
    def test_construct_fixed_source(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol = Solution(
            angular_flux=psi,
            scalar_flux=phi,
            mesh=m,
        )
        assert sol.angular_flux is psi
        assert sol.scalar_flux is phi
        assert sol.boundary_flux is bf
        assert sol.mesh is m
        assert sol.keff is None
        assert sol.history is None

    def test_construct_eigenvalue(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        h = IterationHistory(keff_history=(1.0, 1.05), n_outer=2)
        sol = Solution(
            angular_flux=psi, scalar_flux=phi,
            mesh=m, keff=1.05, history=h,
        )
        assert sol.keff == 1.05
        assert sol.history is h

    def test_mesh_identity_angular_flux(self) -> None:
        """A foreign-mesh TimedFullField must be rejected at construction."""
        m1 = _slab_mesh()
        m2 = _slab_mesh()  # distinct instance, same shape
        state_foreign = TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=m2)
        _, phi, bf = _make_fluxes(m1)
        with pytest.raises(ValueError, match="angular_flux.bulk.mesh"):
            Solution(
                angular_flux=state_foreign, scalar_flux=phi, mesh=m1,
            )

    def test_mesh_identity_scalar_flux(self) -> None:
        m1 = _slab_mesh()
        m2 = _slab_mesh()
        phi_foreign = ScalarFlux.from_mesh(np.zeros((m2.ng, *m2.spatial_shape)), m2)
        psi, _, bf = _make_fluxes(m1)
        with pytest.raises(ValueError, match="scalar_flux.mesh"):
            Solution(
                angular_flux=psi, scalar_flux=phi_foreign, mesh=m1,
            )

    def test_boundary_flux_delegates_to_angular_flux(self) -> None:
        """D-H.1b — ``sol.boundary_flux is sol.angular_flux.boundary``.

        :attr:`Solution.boundary_flux` is a property that returns
        :attr:`TimedFullField.boundary` — the composite's owned
        boundary trace.
        """
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol = Solution(angular_flux=psi, scalar_flux=phi, mesh=m)
        assert sol.boundary_flux is psi.boundary
        assert sol.boundary_flux is bf  # bf = state.boundary from _make_fluxes

    def test_frozen(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol = Solution(angular_flux=psi, scalar_flux=phi, mesh=m)
        with pytest.raises((AttributeError, Exception)):
            sol.keff = 1.0  # type: ignore[misc]


# ════════════════════════════════════════════════════════════════════
# Solution — discrimination
# ════════════════════════════════════════════════════════════════════


class TestSolutionDiscrimination:
    def test_fixed_source_no_keff(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol = Solution(angular_flux=psi, scalar_flux=phi, mesh=m)
        assert sol.is_fixed_source() is True
        assert sol.is_eigenvalue() is False

    def test_eigenvalue_has_keff(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol = Solution(
            angular_flux=psi, scalar_flux=phi,
            mesh=m, keff=1.234,
        )
        assert sol.is_eigenvalue() is True
        assert sol.is_fixed_source() is False

    def test_discrimination_methods_are_mutually_exclusive(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        for keff in (None, 1.0, 0.5, 2.0):
            sol = Solution(
                angular_flux=psi, scalar_flux=phi,
                mesh=m, keff=keff,
            )
            assert sol.is_eigenvalue() != sol.is_fixed_source()


# ════════════════════════════════════════════════════════════════════
# Solution — convergence diagnostics (delegating to history)
# ════════════════════════════════════════════════════════════════════


class TestSolutionDiagnostics:
    def test_dominance_ratio_three_iterations(self) -> None:
        """Solution.dominance_ratio delegates to IterationHistory."""
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        h = IterationHistory(keff_history=(1.0, 1.1, 1.10005))
        sol = Solution(angular_flux=psi, scalar_flux=phi,
                       mesh=m, keff=1.10005, history=h)
        ratio = sol.dominance_ratio()
        assert ratio is not None
        np.testing.assert_allclose(ratio, 5e-5 / 1.1, rtol=1e-10)

    def test_dominance_ratio_no_history(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol = Solution(angular_flux=psi, scalar_flux=phi,
                       mesh=m)
        assert sol.dominance_ratio() is None

    def test_converged_no_history(self) -> None:
        """A solution without history is assumed converged."""
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol = Solution(angular_flux=psi, scalar_flux=phi,
                       mesh=m)
        assert sol.converged() is True

    def test_converged_with_history(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        h_yes = IterationHistory(converged=True)
        h_no = IterationHistory(converged=False)
        sol_yes = Solution(
            angular_flux=psi, scalar_flux=phi,
            mesh=m, history=h_yes,
        )
        sol_no = Solution(
            angular_flux=psi, scalar_flux=phi,
            mesh=m, history=h_no,
        )
        assert sol_yes.converged() is True
        assert sol_no.converged() is False

    def test_keff_history_list(self) -> None:
        """The legacy keff_history accessor returns a plain list."""
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        h = IterationHistory(keff_history=(1.0, 1.05, 1.04))
        sol = Solution(angular_flux=psi, scalar_flux=phi,
                       mesh=m, keff=1.04, history=h)
        lst = sol.keff_history_list()
        assert isinstance(lst, list)
        assert lst == [1.0, 1.05, 1.04]

    def test_keff_history_list_no_history(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol = Solution(angular_flux=psi, scalar_flux=phi,
                       mesh=m)
        assert sol.keff_history_list() == []


# ════════════════════════════════════════════════════════════════════
# Solution — reaction-rate accessor
# ════════════════════════════════════════════════════════════════════


class TestReactionRate:
    def test_reaction_rate_density_shape(self) -> None:
        """σ · φ has the same shape as φ (per-cell rate density)."""
        m = _slab_mesh(nx=4, ng=2)
        psi, _, bf = _make_fluxes(m)
        # Non-trivial flux + non-trivial XS to exercise the einsum
        phi_values = np.arange(
            m.ng * m.nx, dtype=float,
        ).reshape(m.ng, *m.spatial_shape) + 1.0
        phi = ScalarFlux.from_mesh(phi_values, m)
        sol = Solution(angular_flux=psi, scalar_flux=phi,
                       mesh=m)

        xs = np.full((m.ng, *m.spatial_shape), 0.5)
        rate = sol.reaction_rate_density(xs)
        assert rate.shape == (m.ng, *m.spatial_shape)

    def test_reaction_rate_density_math(self) -> None:
        r"""σ · φ at each cell — the named math reads as the formula."""
        m = _slab_mesh(nx=3, ng=2)
        psi, _, bf = _make_fluxes(m)
        phi_values = np.array([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
        ])  # (ng=2, nx=3)
        phi = ScalarFlux.from_mesh(phi_values, m)
        sol = Solution(angular_flux=psi, scalar_flux=phi,
                       mesh=m)

        # σ constant per group → rate[g,x] = σ[g] · φ[g,x]
        xs = np.full((m.ng, *m.spatial_shape), 0.0)
        xs[0, :] = 0.5
        xs[1, :] = 0.25
        rate = sol.reaction_rate_density(xs)

        np.testing.assert_allclose(
            rate[0, :], [0.5, 1.0, 1.5], rtol=1e-15,
        )
        np.testing.assert_allclose(
            rate[1, :], [1.0, 1.25, 1.5], rtol=1e-15,
        )

    def test_reaction_rate_density_zero_flux(self) -> None:
        """σ · 0 = 0 everywhere."""
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m, fill=0.0)
        sol = Solution(angular_flux=psi, scalar_flux=phi,
                       mesh=m)
        xs = np.full((m.ng, *m.spatial_shape), 0.5)
        rate = sol.reaction_rate_density(xs)
        np.testing.assert_array_equal(rate, np.zeros_like(rate))


# ════════════════════════════════════════════════════════════════════
# Solution — compare
# ════════════════════════════════════════════════════════════════════


class TestSolutionCompare:
    def test_compare_identical_within_tolerance(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol_a = Solution(angular_flux=psi, scalar_flux=phi,
                         mesh=m, keff=1.0)
        sol_b = Solution(angular_flux=psi, scalar_flux=phi,
                         mesh=m, keff=1.0)
        diff = sol_a.compare(sol_b, rtol=1e-12)
        assert isinstance(diff, SolutionDiff)
        assert diff.keff_abs == 0.0
        assert diff.angular_flux_linf == 0.0
        assert diff.scalar_flux_linf == 0.0
        assert diff.within_tolerance is True

    def test_compare_different_keff(self) -> None:
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol_a = Solution(angular_flux=psi, scalar_flux=phi,
                         mesh=m, keff=1.0)
        sol_b = Solution(angular_flux=psi, scalar_flux=phi,
                         mesh=m, keff=1.001)
        diff = sol_a.compare(sol_b, rtol=1e-6)
        assert diff.keff_abs is not None
        np.testing.assert_allclose(diff.keff_abs, 0.001, rtol=1e-12)
        assert diff.within_tolerance is False  # 0.001 > 1e-6 * 1.001

    def test_compare_different_flux(self) -> None:
        m = _slab_mesh()
        psi, phi_a, bf = _make_fluxes(m, fill=1.0)
        _, phi_b, _ = _make_fluxes(m, fill=1.1)
        sol_a = Solution(angular_flux=psi, scalar_flux=phi_a,
                         mesh=m)
        sol_b = Solution(angular_flux=psi, scalar_flux=phi_b,
                         mesh=m)
        diff = sol_a.compare(sol_b, rtol=1e-12)
        # phi_a - phi_b = -0.1, L∞ = 0.1
        np.testing.assert_allclose(diff.scalar_flux_linf, 0.1, rtol=1e-12)
        assert diff.within_tolerance is False

    def test_compare_keff_only_one_eigenvalue(self) -> None:
        """When one is fixed-source, keff_abs is None."""
        m = _slab_mesh()
        psi, phi, bf = _make_fluxes(m)
        sol_a = Solution(angular_flux=psi, scalar_flux=phi,
                         mesh=m, keff=1.0)
        sol_b = Solution(angular_flux=psi, scalar_flux=phi,
                         mesh=m)  # fixed-source
        diff = sol_a.compare(sol_b)
        assert diff.keff_abs is None

    def test_compare_cross_mesh_rejected(self) -> None:
        m1 = _slab_mesh()
        m2 = _slab_mesh()
        psi_a, phi_a, bf_a = _make_fluxes(m1)
        psi_b, phi_b, bf_b = _make_fluxes(m2)
        sol_a = Solution(angular_flux=psi_a, scalar_flux=phi_a, mesh=m1)
        sol_b = Solution(angular_flux=psi_b, scalar_flux=phi_b, mesh=m2)
        with pytest.raises(ValueError, match="meshes differ"):
            sol_a.compare(sol_b)
