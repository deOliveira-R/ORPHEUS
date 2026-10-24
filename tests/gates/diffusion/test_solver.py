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

from orpheus.transport.fields.cross_section_field import CrossSectionField

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.diffusion import DiffusionMesh, DiffusionSolver, solve_diffusion_1d
from orpheus.diffusion import ZeroFluxCondition
from orpheus.geometry import StructuredGeometry
from orpheus.geometry.deck import MirrorDeck
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import declare_sides
from orpheus.mesh import CellEdges, CellsByCount, Mesh1D, Mesher
from orpheus.homogeneous import solve_homogeneous_infinite
from orpheus.numerics.eigenvalue import ProductionRateSolver, direct_eigenvalue
from orpheus.numerics.flat_operator import FlattenedOperator
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


def _het_mesh(left, right) -> tuple[Mesh1D, dict]:
    """Fuel [0, 0.5] | reflector [0.5, 3] (two cells, at 1.5) | fuel [3, 5]:
    non-uniform cells, one ``CellEdges`` per material region. Each side is a deck or a condition;
    returns the mesh and the ``{tag: condition}`` its tags stand for (rule R1)."""
    (role_left, role_right), conditions = declare_sides(left, right)
    geometry = StructuredGeometry.slab(
        (0.0, 0.5, 3.0, 5.0), (0, 1, 0), left=role_left, right=role_right,
    )
    mesh = Mesher(geometry).partition((
        CellEdges(np.array([0.0, 0.5])),
        CellEdges(np.array([0.5, 1.5, 3.0])),
        CellEdges(np.array([3.0, 5.0])),
    )).mesh
    return mesh, conditions


def _het_materials():
    return {0: _fuel_2g(), 1: _reflector_2g()}


def _solve_tight(materials, mesh, conditions=None, **kwargs):
    return solve_diffusion_1d(
        materials, mesh,
        keff_tol=1e-13, flux_tol=1e-11, max_outer=5000,
        boundary_conditions=conditions if conditions is not None else {}, **kwargs,
    )


# ═══════════════════════════════════════════════════════════════════════
# Cross-engine gate: power_iteration ≡ direct_eigenvalue at 1e-10
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestEngineCrossGate:
    """The two engines solve the SAME posed (A, F) — k and the full
    composite eigenvector must coincide (the campaign's 1e-10 gate)."""

    @pytest.mark.parametrize(
        "left, right",
        [
            (ZeroFluxCondition(), ZeroFluxCondition()),
            (MirrorDeck(), ZeroFluxCondition()),
            (MirrorDeck(), WallCondition.white(0.3)),
            (WallCondition.vacuum(), WallCondition.vacuum()),
        ],
        ids=["zeroflux", "refl-zeroflux", "refl-albedo", "marshak-vacuum"],
    )
    def test_power_matches_direct(self, left, right):
        materials = _het_materials()
        mesh, conditions = _het_mesh(left, right)

        result = _solve_tight(materials, mesh, conditions)

        solver = DiffusionSolver(DiffusionMesh(mesh, materials, boundary_conditions=conditions))
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
        mesh = Mesher(
            StructuredGeometry.from_homogeneous(2.2, MirrorDeck()),
        ).partition(CellEdges(np.array([0.0, 0.7, 1.5, 2.2]))).mesh
        result = _solve_tight({0: mix}, mesh)
        reference = solve_homogeneous_infinite(mix)

        assert abs(result.keff - reference.k_inf) < 1e-11

        # The mode is spatially flat with the homogeneous spectrum.
        phi = result.flux.interior.values                     # (ng, 3)
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
        mesh = Mesher(
            StructuredGeometry.from_homogeneous(40.0, 0),
        ).partition(CellsByCount.uniform_width(8)).mesh
        result = _solve_tight({0: mix}, mesh, {0: ZeroFluxCondition()})
        # Leaky ⟹ strictly below k∞; positive fundamental mode.
        assert 0.0 < result.keff < solve_homogeneous_infinite(mix).k_inf
        assert np.all(result.flux.interior.values > 0.0)


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
            _het_materials(), *_het_mesh(WallCondition.vacuum(), WallCondition.vacuum()),
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
            *_het_mesh(MirrorDeck(), WallCondition.white(alpha)),
        )
        trace = result.flux.boundary
        np.testing.assert_allclose(
            trace.inflow_view("xmax"),
            alpha * trace.outflow_view("xmax"),
            rtol=1e-10,
        )

    def test_reflective_has_zero_net_current(self):
        result = _solve_tight(
            _het_materials(), *_het_mesh(MirrorDeck(), MirrorDeck()),
        )
        trace = result.flux.boundary
        atol = 1e-12 * self._trace_scale(trace)
        for face in ("xmin", "xmax"):
            np.testing.assert_allclose(
                trace.net_current(face), 0.0, atol=atol,
            )

    def test_zero_flux_law_zeroes_the_boundary_scalar_flux(self):
        result = _solve_tight(
            _het_materials(), *_het_mesh(ZeroFluxCondition(), ZeroFluxCondition()),
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
            _het_materials(), *_het_mesh(MirrorDeck(), ZeroFluxCondition()),
        )
        trace = result.flux.boundary
        atol = 1e-12 * self._trace_scale(trace)
        np.testing.assert_allclose(trace.net_current("xmin"), 0.0, atol=atol)
        np.testing.assert_allclose(
            trace.p1_boundary_scalar_flux("xmax"), 0.0, atol=atol,
        )
        # The flux sags toward the Dirichlet sink: higher at the
        # reflective edge than at the zero-flux edge, every group.
        phi = result.flux.interior.values
        assert np.all(phi[:, 0] > phi[:, -1])

    def test_current_profile_is_reconstructed_at_every_face(self):
        result = _solve_tight(
            _het_materials(), *_het_mesh(MirrorDeck(), ZeroFluxCondition()),
        )
        ng, nx = result.flux.interior.values.shape
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
        mesh, conditions = _het_mesh(WallCondition.vacuum(), WallCondition.white(0.5))
        result = _solve_tight(materials, mesh, conditions)

        material_mesh = result.mesh
        mat_xs = material_mesh.mat_xs
        phi = result.flux.interior.values
        production = IntegratedReactionRate(
            CrossSectionField.fission_production_of(mat_xs)
        ).evaluate(phi)
        absorption = IntegratedReactionRate(
            CrossSectionField.absorption_of(mat_xs)
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
            return _solve_tight(materials, *_het_mesh(bc, bc)).keff

        k_zero_flux = k(ZeroFluxCondition())
        k_vacuum = k(WallCondition.vacuum())
        k_albedo = k(WallCondition.white(0.6))
        k_reflective = k(MirrorDeck())
        assert k_zero_flux < k_vacuum < k_albedo < k_reflective

    def test_returned_mode_is_production_normalized(self):
        """The #270 contract: power_iteration renormalizes through
        compute_production_rate, so the returned mode carries
        ∫νΣf·φ dV = 1 (the island's ``fi /= max`` hack is retired)."""
        result = _solve_tight(
            _het_materials(), *_het_mesh(ZeroFluxCondition(), ZeroFluxCondition()),
        )
        solver = DiffusionSolver(result.mesh)
        production = solver.compute_production_rate(
            np.asarray(result.flux.to_flat())
        )
        assert abs(production - 1.0) < 1e-12


@pytest.mark.foundation
class TestProtocolAndDefaults:
    def test_solver_satisfies_the_production_rate_contract(self):
        mesh, conditions = _het_mesh(MirrorDeck(), WallCondition.vacuum())
        solver = DiffusionSolver(
            DiffusionMesh(mesh, _het_materials(), boundary_conditions=conditions)
        )
        assert isinstance(solver, ProductionRateSolver)

    # (BC-tag / albedo-parameter / multi-D refusals moved to the PHASE
    # SPACE at #290 P7a — they fire at DiffusionMesh construction and
    # are gated in test_augmented_mesh.py.)
