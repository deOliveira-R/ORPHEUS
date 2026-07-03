r"""The modern diffusion k-eigenvalue solver — engine + solution gates (#290 P5).

Three gate families over :class:`orpheus.diffusion.DiffusionSolver` /
:func:`orpheus.diffusion.solve_diffusion_1d`:

* **cross-engine** — the protocol-driven
  :func:`~orpheus.numerics.eigenvalue.power_iteration` path ≡
  :func:`~orpheus.numerics.eigenvalue.direct_eigenvalue` on the
  materialized ``(A, F)`` pair, ``|Δk| < 1e-10`` + the full composite
  eigenvector, across the albedo family (the campaign's committed
  engine gate);
* **L2 infinite-medium** — reflective diffusion k ≡ homogeneous
  :math:`k_\infty` (structurally independent loss posings sharing the
  verified K_iso kernels; the homogeneous side is anchored to a SymPy
  closed form). Includes the **3-group asymmetric case** — the
  discriminator the legacy island's hardcoded 2-group ``[::-1]`` flip
  trick structurally cannot represent;
* **solution semantics** — per-law trace identities at the converged
  mode (vacuum :math:`J^-=0`, albedo :math:`J^-=\alpha J^+`,
  reflective :math:`J_{\rm net}=0`, zero-flux :math:`\phi_\Gamma=0`),
  the asymmetric-BC face-mapping pin, the integrated balance identity
  ``P/k = absorption + leakage``, k-monotonicity in the albedo, and
  the ProductionRateSolver contract (#270 arm).

Conventions: ``.claude/plans/diffusion_crosswalk.md``.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.diffusion import DiffusionSolver, solve_diffusion_1d
from orpheus.geometry import BC
from orpheus.geometry.mesh import Mesh1D, Mesh2D
from orpheus.homogeneous import solve_homogeneous_infinite
from orpheus.numerics.eigenvalue import ProductionRateSolver, direct_eigenvalue
from orpheus.numerics.flat_operator import FlattenedOperator
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.reaction_rate_functional import IntegratedReactionRate

# ═══════════════════════════════════════════════════════════════════════
# Materials — heterogeneous asymmetric 2G (the P4 anti-degeneracy tables)
# and a fully-asymmetric 3G table (distinct χ spread, skip-scatter 1→3)
# whose scatter structure has NO 2-group flip representation.
# ═══════════════════════════════════════════════════════════════════════

_SIG_T_A = np.array([0.2181, 0.7850])
_SIG_S_A = np.array([[0.1900, 0.0160], [0.0, 0.4200]])
_SIG_F_A = np.array([0.0024, 0.0489])

_SIG_T_B = np.array([0.3416, 0.9431])
_SIG_S_B = np.array([[0.1000, 0.0020], [0.0, 0.0500]])

_SIG_T_3G = np.array([0.30, 0.55, 1.10])
_SIG_S_3G = np.array([
    [0.10, 0.05, 0.01],
    [0.00, 0.30, 0.08],
    [0.00, 0.00, 0.60],
])
_SIG_F_3G = np.array([0.002, 0.010, 0.120])
_NU_3G = np.array([2.9, 2.5, 2.4])
_CHI_3G = np.array([0.7, 0.3, 0.0])


def _fuel_2g():
    return make_mixture(
        sig_t=_SIG_T_A,
        sig_c=_SIG_T_A - _SIG_F_A - _SIG_S_A.sum(axis=1),
        sig_f=_SIG_F_A, nu=np.array([2.54, 2.47]),
        chi=np.array([1.0, 0.0]), sig_s=_SIG_S_A,
    )


def _reflector_2g():
    return make_mixture(
        sig_t=_SIG_T_B,
        sig_c=_SIG_T_B - _SIG_S_B.sum(axis=1),
        sig_f=np.zeros(2), nu=np.zeros(2),
        chi=np.zeros(2), sig_s=_SIG_S_B,
    )


def _fuel_3g():
    return make_mixture(
        sig_t=_SIG_T_3G,
        sig_c=_SIG_T_3G - _SIG_F_3G - _SIG_S_3G.sum(axis=1),
        sig_f=_SIG_F_3G, nu=_NU_3G,
        chi=_CHI_3G, sig_s=_SIG_S_3G,
    )


_EDGES_HET = np.array([0.0, 0.5, 1.5, 3.0, 5.0])   # non-uniform
_MAT_IDS_HET = np.array([0, 1, 1, 0])


def _het_mesh(bc_left: BC, bc_right: BC) -> Mesh1D:
    return Mesh1D(
        edges=_EDGES_HET, mat_ids=_MAT_IDS_HET,
        bc_left=bc_left, bc_right=bc_right,
    )


def _het_materials():
    return {0: _fuel_2g(), 1: _reflector_2g()}


def _solve_tight(materials, mesh, **kwargs):
    return solve_diffusion_1d(
        materials, mesh,
        keff_tol=1e-13, flux_tol=1e-11, max_outer=5000, **kwargs,
    )


# ═══════════════════════════════════════════════════════════════════════
# Cross-engine gate: power_iteration ≡ direct_eigenvalue at 1e-10
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestEngineCrossGate:
    """The two engines solve the SAME posed (A, F) — k and the full
    composite eigenvector must coincide (the campaign's 1e-10 gate)."""

    @pytest.mark.parametrize(
        "bc_left, bc_right",
        [
            (BC("zero_flux"), BC("zero_flux")),
            (BC("reflective"), BC("zero_flux")),
            (BC("reflective"), BC("albedo", {"albedo": 0.3})),
            (BC("vacuum"), BC("vacuum")),
        ],
        ids=["zeroflux", "refl-zeroflux", "refl-albedo", "marshak-vacuum"],
    )
    def test_power_matches_direct(self, bc_left, bc_right):
        materials = _het_materials()
        mesh = _het_mesh(bc_left, bc_right)

        result = _solve_tight(materials, mesh)

        solver = DiffusionSolver(MaterialMesh(mesh, materials))
        A_mat = FlattenedOperator(solver.loss, solver.template).as_matrix()
        F_mat = FlattenedOperator(solver.fission, solver.template).as_matrix()
        k_direct, phi_direct = direct_eigenvalue(A_mat, F_mat)

        assert abs(result.keff - k_direct) < 1e-10

        # The full composite eigenvector (bulk ⊕ trace), both normalized
        # to unit production rate — the fundamental mode is unique.
        flat_power = np.asarray(result.flux.to_flat())
        production_direct = solver.compute_production_rate(phi_direct)
        phi_direct_normed = phi_direct / production_direct
        np.testing.assert_allclose(
            flat_power, phi_direct_normed,
            rtol=1e-6, atol=1e-8 * float(np.abs(phi_direct_normed).max()),
        )


# ═══════════════════════════════════════════════════════════════════════
# L2 infinite-medium: reflective diffusion ≡ homogeneous k∞
# (independent loss posings — L+C−S−B on the composite vs C−K_iso
# meshless — sharing the verified K_iso kernels; #276 anchored the
# homogeneous side to a SymPy closed form.)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l2
class TestInfiniteMedium:
    @pytest.mark.parametrize(
        "mixture",
        [pytest.param(_fuel_2g, id="2g"), pytest.param(_fuel_3g, id="3g")],
    )
    def test_reflective_slab_reproduces_k_infinity(self, mixture):
        mix = mixture()
        ng = mix.ng
        # Non-uniform mesh — the constant mode is annihilated by L on ANY
        # spacing (the P4 gate), so k∞ must be reproduced exactly.
        mesh = Mesh1D(
            edges=np.array([0.0, 0.7, 1.5, 2.2]),
            mat_ids=np.zeros(3, dtype=int),
            bc_left=BC("reflective"), bc_right=BC("reflective"),
        )
        result = _solve_tight({0: mix}, mesh)
        reference = solve_homogeneous_infinite(mix)

        assert abs(result.keff - reference.k_inf) < 1e-11

        # The mode is spatially flat with the homogeneous spectrum.
        phi = result.flux.bulk.values                     # (ng, 3)
        assert float(np.ptp(phi, axis=1).max()) < 1e-10 * float(phi.max())
        spectrum = phi[:, 0] / phi[:, 0].sum()
        expected = reference.flux / reference.flux.sum()
        np.testing.assert_allclose(spectrum, expected.reshape(ng), rtol=1e-9)

    def test_three_group_is_structurally_beyond_the_island(self):
        """The 3G discriminator: the modern path is ng-generic by
        construction (K_iso); the island's ``sig_s[::-1]`` flip trick is
        2-group-only. This solve exercising ng=3 end-to-end (with
        skip-scatter 1→3 and a split χ) IS the structural kill of the
        flip trick — no 2G rearrangement reproduces it."""
        mix = _fuel_3g()
        mesh = Mesh1D(
            edges=np.linspace(0.0, 40.0, 9),
            mat_ids=np.zeros(8, dtype=int),
            bc_left=BC("zero_flux"), bc_right=BC("zero_flux"),
        )
        result = _solve_tight({0: mix}, mesh)
        # Leaky ⟹ strictly below k∞; positive fundamental mode.
        assert 0.0 < result.keff < solve_homogeneous_infinite(mix).k_inf
        assert np.all(result.flux.bulk.values > 0.0)


# ═══════════════════════════════════════════════════════════════════════
# Solution semantics: trace identities, balance, ordering, protocol
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestSolutionTraceSemantics:
    """The converged composite carries its boundary law — typed,
    inspectable, and exact to LU roundoff (the inflow rows of ``Aψ`` are
    constraint rows with zero source)."""

    def _trace_scale(self, trace) -> float:
        return float(np.abs(trace.values).max())

    def test_vacuum_means_zero_incoming(self):
        result = _solve_tight(
            _het_materials(), _het_mesh(BC("vacuum"), BC("vacuum")),
        )
        trace = result.flux.boundary
        atol = 1e-12 * self._trace_scale(trace)
        for face in ("xmin", "xmax"):
            np.testing.assert_allclose(
                trace.inflow_view(face), 0.0, atol=atol,
            )

    def test_albedo_law_holds_at_solution(self):
        alpha = 0.3
        result = _solve_tight(
            _het_materials(),
            _het_mesh(BC("reflective"), BC("albedo", {"albedo": alpha})),
        )
        trace = result.flux.boundary
        np.testing.assert_allclose(
            trace.inflow_view("xmax"),
            alpha * trace.outflow_view("xmax"),
            rtol=1e-10,
        )

    def test_reflective_has_zero_net_current(self):
        result = _solve_tight(
            _het_materials(), _het_mesh(BC("reflective"), BC("reflective")),
        )
        trace = result.flux.boundary
        atol = 1e-12 * self._trace_scale(trace)
        for face in ("xmin", "xmax"):
            np.testing.assert_allclose(
                trace.net_current(face), 0.0, atol=atol,
            )

    def test_zero_flux_law_zeroes_the_boundary_scalar_flux(self):
        result = _solve_tight(
            _het_materials(), _het_mesh(BC("zero_flux"), BC("zero_flux")),
        )
        trace = result.flux.boundary
        atol = 1e-12 * self._trace_scale(trace)
        for face in ("xmin", "xmax"):
            np.testing.assert_allclose(
                trace.p1_boundary_scalar_flux(face), 0.0, atol=atol,
            )

    def test_asymmetric_bcs_pin_the_face_mapping(self):
        """Reflective LEFT + zero-flux RIGHT: both laws must land on
        their own faces (a face-mapping swap satisfies neither)."""
        result = _solve_tight(
            _het_materials(), _het_mesh(BC("reflective"), BC("zero_flux")),
        )
        trace = result.flux.boundary
        atol = 1e-12 * self._trace_scale(trace)
        np.testing.assert_allclose(trace.net_current("xmin"), 0.0, atol=atol)
        np.testing.assert_allclose(
            trace.p1_boundary_scalar_flux("xmax"), 0.0, atol=atol,
        )
        # The flux sags toward the Dirichlet sink: higher at the
        # reflective edge than at the zero-flux edge, every group.
        phi = result.flux.bulk.values
        assert np.all(phi[:, 0] > phi[:, -1])

    def test_current_profile_is_reconstructed_at_every_face(self):
        result = _solve_tight(
            _het_materials(), _het_mesh(BC("reflective"), BC("zero_flux")),
        )
        ng, nx = result.flux.bulk.values.shape
        assert result.current.shape == (ng, nx + 1)
        trace = result.flux.boundary
        # Boundary slots are the axis-signed trace net currents …
        np.testing.assert_array_equal(
            result.current[:, 0], -trace.net_current("xmin"),
        )
        np.testing.assert_array_equal(
            result.current[:, -1], +trace.net_current("xmax"),
        )
        # … and the reflective face carries none, while the sink drains
        # rightward (positive axis-signed current at xmax).
        atol = 1e-12 * float(np.abs(result.current).max())
        np.testing.assert_allclose(result.current[:, 0], 0.0, atol=atol)
        assert np.all(result.current[:, -1] > 0.0)


@pytest.mark.foundation
class TestBalanceAndOrdering:
    def test_integrated_balance_identity(self):
        """P/k = absorption + leakage at the converged mode (the
        column-sum theorem + divergence telescoping, end-to-end)."""
        materials = _het_materials()
        mesh = _het_mesh(BC("vacuum"), BC("albedo", {"albedo": 0.5}))
        result = _solve_tight(materials, mesh)

        material_mesh = result.mesh
        mat_xs = material_mesh.material_xs_field()
        phi = result.flux.bulk.values
        production = IntegratedReactionRate(
            mat_xs.fission_production_field
        ).evaluate(phi)
        absorption = IntegratedReactionRate(
            mat_xs.absorption_cross_section_field
        ).evaluate(phi)

        trace = result.flux.boundary
        areas = material_mesh.areas
        leakage = float(
            areas[0] * trace.net_current("xmin").sum()
            + areas[-1] * trace.net_current("xmax").sum()
        )

        lhs = production / result.keff
        rhs = absorption + leakage
        assert abs(lhs - rhs) < 1e-10 * abs(rhs)

    def test_keff_is_monotone_in_the_albedo(self):
        """𝒜 orders the eigenvalue: zero-flux < vacuum < albedo(0.6)
        < reflective (more return current ⟹ less leakage ⟹ larger k)."""
        materials = _het_materials()

        def k(bc):
            return _solve_tight(materials, _het_mesh(bc, bc)).keff

        k_zero_flux = k(BC("zero_flux"))
        k_vacuum = k(BC("vacuum"))
        k_albedo = k(BC("albedo", {"albedo": 0.6}))
        k_reflective = k(BC("reflective"))
        assert k_zero_flux < k_vacuum < k_albedo < k_reflective

    def test_returned_mode_is_production_normalized(self):
        """The #270 contract: power_iteration renormalizes through
        compute_production_rate, so the returned mode carries
        ∫νΣf·φ dV = 1 (the island's ``fi /= max`` hack is retired)."""
        result = _solve_tight(
            _het_materials(), _het_mesh(BC("zero_flux"), BC("zero_flux")),
        )
        solver = DiffusionSolver(result.mesh)
        production = solver.compute_production_rate(
            np.asarray(result.flux.to_flat())
        )
        assert abs(production - 1.0) < 1e-12


@pytest.mark.foundation
class TestProtocolAndRefusals:
    def test_solver_satisfies_the_production_rate_contract(self):
        solver = DiffusionSolver(
            MaterialMesh(
                _het_mesh(BC("reflective"), BC("vacuum")), _het_materials(),
            )
        )
        assert isinstance(solver, ProductionRateSolver)

    def test_unsupported_bc_kind_is_refused_with_the_supported_list(self):
        # "white" is DELIBERATELY absent: at P1 it coincides with
        # reflective (the P3 realizer's coincidence note).
        mesh = _het_mesh(BC("white"), BC("vacuum"))
        with pytest.raises(ValueError, match="'white'.*Supported.*albedo"):
            DiffusionSolver(MaterialMesh(mesh, _het_materials()))

    def test_albedo_without_parameter_is_refused(self):
        mesh = _het_mesh(BC("albedo"), BC("vacuum"))
        with pytest.raises(ValueError, match="albedo.*parameter"):
            DiffusionSolver(MaterialMesh(mesh, _het_materials()))

    def test_multi_dimensional_mesh_is_refused(self):
        mesh2d = Mesh2D(
            edges_x=np.array([0.0, 1.0, 2.0]),
            edges_y=np.array([0.0, 1.0]),
            mat_map=np.zeros((2, 1), dtype=int),
        )
        with pytest.raises(ValueError, match="1-D"):
            DiffusionSolver(MaterialMesh(mesh2d, {0: _fuel_2g()}))

    def test_undeclared_bcs_default_to_reflective(self):
        """The infinite-lattice convention (the SN default, mirrored)."""
        materials = _het_materials()
        bare = Mesh1D(edges=_EDGES_HET, mat_ids=_MAT_IDS_HET)
        explicit = _het_mesh(BC("reflective"), BC("reflective"))
        k_bare = _solve_tight(materials, bare).keff
        k_explicit = _solve_tight(materials, explicit).keff
        assert abs(k_bare - k_explicit) < 1e-13
