"""Tests for the SN boundary condition infrastructure.

Verifies the declaration-to-realization path: the roles declared on the geometry (a deck, or a tag whose condition is declared at
the solve), their resolution at SNProblem construction through the sealed realizer table, and correct behavior in sweeps.
(Wave 2 S2b: the string registry ``BOUNDARY_OPERATOR_REGISTRY`` is retired; its successor is gated in ``test_sn_boundary_realizer.py``.)
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh2D, Mesher
from orpheus.sn.boundary.realized_operators import RankOneReemission, ScaledOrdinatePermutation, ZeroMap
from orpheus.sn.problem import SNProblem
from orpheus.numerics.quadrature import Quadrature
from tests.gates.sn._test_helpers import placeholder_materials
from orpheus.geometry.deck import MirrorDeck
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for

#: Exterior conditions of the tags this file's meshes carry (wave 2: a face's role is a deck or an int tag; the condition a tag
#: stands for is declared at the solve, keyword ``boundary_conditions``).
_CONDITIONS = {0: WallCondition.vacuum(), 1: WallCondition.white(1.0), 2: WallCondition.specular(1.0)}


# SN boundary-condition infrastructure: structural invariants of the
# SNProblem BC wiring (no theory-page :label:). Foundation, not a physics
# equation gate. (Was a V&V orphan before the taxonomy reorg forced a
# marker — see .claude/plans/archive/sn_test_taxonomy.md.)
pytestmark = pytest.mark.foundation


@pytest.fixture
def quad():
    return Quadrature.gauss_legendre(4)


@pytest.fixture
def slab_mesh():
    return Mesher(StructuredGeometry.slab(
        (0.0, 5.0), (0,), left=MirrorDeck(), right=MirrorDeck(),
    )).partition(CellsByCount.uniform_width(10)).mesh


# ═══════════════════════════════════════════════════════════════════════
# Registry tests: RETIRED (wave 2 S2b)
# ═══════════════════════════════════════════════════════════════════════
#
# ``TestSNBCRegistry`` (``test_registry_keys``, ``test_registry_docstrings``, ``test_registry_programmatic_query``) asserted that
# ``SNProblem.BOUNDARY_OPERATOR_REGISTRY`` held "vacuum" and "reflective" and that every factory carried a docstring "used as a
# description for UI query". The registry was a string-keyed table of factories; it is replaced by the sealed
# ``RealizerTable`` ``SN_REALIZERS`` keyed by UNIT TYPE, whose membership and door are gated in ``test_sn_boundary_realizer.py``.
# The UI-description claim has no successor: nothing queried it (0 consumers in the tree at the retirement, ``grep`` of
# ``BOUNDARY_OPERATOR_REGISTRY`` over ``orpheus/``, ``tests/`` and ``docs/``).


# ═══════════════════════════════════════════════════════════════════════
# Resolution tests
# ═══════════════════════════════════════════════════════════════════════

class TestSNBCResolution:
    """BC resolution at SNProblem construction time."""

    def test_declared_reflective(self, slab_mesh, quad):
        """The mirror deck declared on both faces resolves to the unit ordinate permutation."""
        sn = SNProblem(slab_mesh, quad, placeholder_materials(), boundary_conditions=conditions_for(slab_mesh, _CONDITIONS))
        for face in ("xmin", "xmax"):
            assert isinstance(sn.bc[face], ScaledOrdinatePermutation) and sn.bc[face].alpha == 1.0

    def test_explicit_vacuum(self, slab_mesh, quad):
        mesh = replace(slab_mesh, face_roles={"xmin": 0, "xmax": 0})
        sn = SNProblem(mesh, quad, placeholder_materials(), boundary_conditions=conditions_for(mesh, _CONDITIONS))
        assert isinstance(sn.bc["xmin"], ZeroMap)
        assert isinstance(sn.bc["xmax"], ZeroMap)

    def test_mixed_bcs(self, slab_mesh, quad):
        mesh = replace(slab_mesh, face_roles={"xmin": MirrorDeck(), "xmax": 0})
        sn = SNProblem(mesh, quad, placeholder_materials(), boundary_conditions=conditions_for(mesh, _CONDITIONS))
        assert isinstance(sn.bc["xmin"], ScaledOrdinatePermutation)
        assert isinstance(sn.bc["xmax"], ZeroMap)

    def test_white_is_realized_now(self, slab_mesh, quad):
        """Meaning change at wave 2 S2b. ``test_unknown_bc_raises`` pinned the REFUSAL of a white face ("does not support ...
        'white'"): the string registry listed two kinds. The diffuse chain is registered in the realizer table, so a white wall on a
        Cartesian face realizes to the rank-one re-emission; what SN still refuses is the diffuse chain on a RADIAL outer surface
        (issue #653, below)."""
        mesh = replace(slab_mesh, face_roles={"xmin": 1, "xmax": MirrorDeck()})
        sn = SNProblem(mesh, quad, placeholder_materials(), boundary_conditions=conditions_for(mesh, _CONDITIONS))
        assert isinstance(sn.bc["xmin"], RankOneReemission)

    def test_a_diffuse_wall_on_a_radial_outer_surface_is_refused(self, quad):
        """Re-posed from ``test_error_lists_supported`` (the refusal listed the registry's supported kinds when a periodic tag was
        met; the half-declared translation pair is now a geometry-level ``DeckError``, gated in ``test_snmesh_realizer_wiring.py``).
        The refusal SN keeps is the scope boundary of issue #653: a diffuse (white) wall at a radial outer surface."""
        from orpheus.transport.boundary_conditions import BoundaryRefused

        geom = StructuredGeometry.sphere((0.0, 1.0), (0,), outer=1)
        mesh = Mesher(geom).partition(CellsByCount.uniform_width(5)).mesh
        with pytest.raises(BoundaryRefused, match="corner action is ruled only for a permuting return or a zero return") as refused:
            SNProblem(mesh, quad, placeholder_materials(), boundary_conditions=conditions_for(mesh, _CONDITIONS))
        assert refused.value.admission.issue == 653

    def test_2d_mesh_resolution(self):
        """2-D Cartesian BC resolution. Needs a genuine-2-D quadrature:
        the y-faces require ordinates with non-zero mu_y, which the shared
        1-D ``gauss_legendre`` fixture lacks (the trace-space guard at
        ``angular_trace_space.py`` correctly rejects it). ``level_symmetric``
        carries genuine mu_y."""
        quad = Quadrature.level_symmetric(sn_order=4)
        mesh = Mesh2D(
            edges_x=np.linspace(0, 2, 3),
            edges_y=np.linspace(0, 2, 3),
            mat_map=np.zeros((2, 2), dtype=int),
            face_roles={
                "xmin": MirrorDeck(),
                "xmax": 0,
                "ymin": MirrorDeck(),
                "ymax": 0,
            },
        )
        sn = SNProblem(mesh, quad, placeholder_materials(), boundary_conditions=conditions_for(mesh, _CONDITIONS))
        assert isinstance(sn.bc["xmin"], ScaledOrdinatePermutation)
        assert isinstance(sn.bc["xmax"], ZeroMap)
        assert isinstance(sn.bc["ymin"], ScaledOrdinatePermutation)
        assert isinstance(sn.bc["ymax"], ZeroMap)

    def test_curvilinear_vacuum_resolves(self, quad):
        """Spherical/cylindrical accept vacuum; sweep enforces zero-incoming
        at the outer face. Vacuum support added in commits 655e3e5 / 37c5bbf
        (the curvilinear-only-reflective gate was removed and the inward
        sweep now branches on ``is_vacuum_outer``)."""
        geom = StructuredGeometry.sphere((0.0, 1.0), (0,), outer=0)
        mesh = Mesher(geom).partition(CellsByCount.uniform_width(5)).mesh
        sn = SNProblem(mesh, quad, placeholder_materials(), boundary_conditions=conditions_for(mesh, _CONDITIONS))
        assert isinstance(sn.bc["xmax"], ZeroMap)


# ═══════════════════════════════════════════════════════════════════════
# Sweep behavior tests
# ═══════════════════════════════════════════════════════════════════════

def _err052_fixture():
    """The ERR-052 configuration: a 2 cm 2-group slab, reflective vs vacuum.

    Returns ``(materials, mesh_refl, mesh_vac, cell_volumes, sigma_p)``.
    ``sigma_p`` is the mixture's production cross-section
    :math:`\\nu\\Sigma_f` and ``cell_volumes`` the slab cell widths — the
    two ingredients of the hand-computed production rate, taken from the
    MIXTURE and the MESH, never from the solver.
    """
    import numpy as _np
    from orpheus.derivations.reference_values import get

    case = get("sn_slab_2eg_1rg")
    mix = next(iter(case.materials.values()))
    materials = {0: mix}
    mesh_refl = Mesher(StructuredGeometry.slab(
        (0.0, 2.0), (0,), left=MirrorDeck(), right=MirrorDeck(),
    )).partition(CellsByCount.uniform_volume(20)).mesh
    mesh_vac = replace(mesh_refl, face_roles={"xmin": 0, "xmax": 0})
    # The fixture mixture is (n,2n)-free, so total production IS fission
    # production — asserted here so the hand formula below stays honest
    # if the reference case ever gains a Σ₂ channel.
    assert mix.Sig2[0].nnz == 0, (
        "the hand-computed production rate below omits the (n,2n) "
        "emission term; the fixture mixture has grown a Σ₂ channel"
    )
    return (
        materials, mesh_refl, mesh_vac,
        _np.diff(mesh_refl.edges), _np.asarray(mix.SigP),
    )


def _hand_production_rate(result, cell_volumes, sigma_p) -> float:
    r"""``∫_V Σ_g νΣ_{f,g} φ_g dV`` from the SOLUTION's scalar flux.

    Deliberately hand-assembled from ``result.scalar_flux.values``
    ``(ng, nx)``, the mesh cell widths and the mixture's ``SigP`` — NOT
    via :meth:`SNSolver.compute_production_rate`, which is the very
    routine the power iteration divides by. Structural independence
    (``vv-principles`` L11): a reference computed by the normaliser
    would make the gate below a tautology.
    """
    import numpy as _np

    phi = _np.asarray(result.scalar_flux.values)
    assert phi.shape == (sigma_p.size, cell_volumes.size), phi.shape
    return float((sigma_p[:, None] * phi * cell_volumes[None, :]).sum())


class TestSNBCSweepBehavior:
    """Verify that resolved BCs produce correct sweep behavior."""

    def test_vacuum_keff_lower_than_reflective(self, quad):
        """Vacuum BC loses neutrons → lower keff than reflective.

        B0.3 NOTE — this test used to carry ``@catches("ERR-052")``.
        **It never caught ERR-052** (measured): re-introducing the exact
        documented bug moves the instrumentation decisively
        (renormalisation calls 6 → 0, ``|φ|max`` 7.60 → 0.61) and this
        assertion stays GREEN, because it is an ORDERING comparison with
        a ~10× margin (1.875 > 0.164) and the ordering is a
        leakage-physics fact that survives a badly-scaled iterate.

        It was a true catcher when written; the *fixture* drifted out of
        the regime. The marker now lives on
        :meth:`test_power_iteration_renormalises_to_unit_production_rate`
        below, which asserts on the mechanism instead. This row keeps its
        own honest claim — the leakage ordering — and no coverage claim
        it cannot back (``vv-principles``: a ``catches`` marker is a
        COVERAGE CLAIM, not a topic tag, and it has a SHELF LIFE).
        """
        from orpheus.sn.solver import solve_sn

        materials, mesh_refl, mesh_vac, _, _ = _err052_fixture()
        result_refl = solve_sn(materials, mesh_refl, quad, boundary_conditions=conditions_for(mesh_refl, _CONDITIONS))
        result_vac = solve_sn(materials, mesh_vac, quad, boundary_conditions=conditions_for(mesh_vac, _CONDITIONS))

        # Reflective has higher keff (no leakage vs leakage)
        assert result_refl.outcome.keff > result_vac.outcome.keff

    @pytest.mark.catches("ERR-052")
    @pytest.mark.parametrize("tol", [1e-7, 1e-12])
    def test_power_iteration_renormalises_to_unit_production_rate(
        self, quad, tol,
    ):
        r"""ERR-052 catcher, re-posed onto the MECHANISM (B0.3).

        ERR-052 is *power iteration without per-step flux
        renormalisation*. The production fix
        (``orpheus/numerics/eigenvalue.py``, the
        ``isinstance(solver, ProductionRateSolver)`` branch) divides the
        iterate by its production rate every outer step, and its
        DOCUMENTED consequence is an output convention:

        .. math::

            P(\phi) \;=\; \int_V \sum_g \nu\Sigma_{f,g}\,\phi_g \,dV
            \;=\; 1 .

        That is the instrumented quantity the bug moves, so it is what
        this test asserts — with **no margin** (``rtol=1e-12``) and, by
        construction, **no regime dependence**: :math:`P(\phi)=1` holds
        at every outer count and every tolerance, so this marker cannot
        go blind the way its predecessor did when the fixture drifted to
        6 outers.

        MUTATION-VERIFIED (the ``catches`` recipe). Re-introducing the
        exact documented bug in-process — making the
        ``ProductionRateSolver`` narrowing never match, i.e. the pre-fix
        un-normalised trajectory verbatim — moves the measured quantity:

        =============  ==========  ==========================
        leg            unmutated   ERR-052 re-introduced
        =============  ==========  ==========================
        reflective     1.0         0.84375     (16 % low)
        vacuum         1.0         0.0739524   **13.5× low**
        =============  ==========  ==========================

        Both legs red at ``rtol=1e-12``.

        HONEST SCOPE — a finding recorded so the next reader does not
        chase it: the catalogue's *denormalised-FP underflow after
        ~30-60 outers* is **NOT reachable in this fixture at any
        depth**. Measured on the vacuum leg with the bug re-introduced,
        ``|φ|max`` sits at 5.6552e-01 after 6, 24 and 32 outers — the
        un-normalised iterate reaches a stable magnitude rather than
        decaying without bound, because once ``k`` has converged the
        ``F·φ/k`` step is scale-neutral. So no tightening of the
        tolerance could ever have rescued the old ordering assertion;
        only moving to the mechanism does. The ``tol`` parametrisation
        (6 outers vs 32) is kept as the regime-independence proof, not
        as a route to the underflow.

        The BOTH-BC coverage is deliberate and is the catalogue's own
        lesson ("every BC type needs its own eigenvalue test at
        multi-group"): reflective is the easy supercritical case, vacuum
        the subcritical one where the un-normalised iterate decays.
        """
        import numpy as np

        from orpheus.sn.solver import solve_sn

        materials, mesh_refl, mesh_vac, volumes, sigma_p = _err052_fixture()
        for label, mesh in (("reflective", mesh_refl), ("vacuum", mesh_vac)):
            result = solve_sn(
                materials, mesh, quad,
                keff_tol=tol, flux_tol=tol, max_outer=200, boundary_conditions=conditions_for(mesh, _CONDITIONS),
            )
            production = _hand_production_rate(result, volumes, sigma_p)
            np.testing.assert_allclose(
                production, 1.0, rtol=1e-12, atol=0.0,
                err_msg=(
                    f"{label} BC, tol={tol:g}: the converged flux does "
                    f"NOT carry the ERR-052 unit-production-rate "
                    f"convention — P(phi)={production!r} != 1. The "
                    f"per-step renormalisation in power_iteration is "
                    f"missing or was bypassed (check the "
                    f"ProductionRateSolver narrowing)."
                ),
            )
