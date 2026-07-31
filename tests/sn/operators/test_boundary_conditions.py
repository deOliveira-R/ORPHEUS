"""Tests for the SN boundary condition infrastructure.

Verifies the BOUNDARY_OPERATOR_REGISTRY pattern: declaration on geometry, resolution at
SNMesh construction, and correct behavior in sweeps.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Mesh2D, CoordSystem
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.numerics.quadrature import Quadrature
from tests.sn._test_helpers import placeholder_materials

# SN boundary-condition infrastructure: structural invariants of the
# SNMesh BC wiring (no theory-page :label:). Foundation, not a physics
# equation gate. (Was a V&V orphan before the taxonomy reorg forced a
# marker — see .claude/plans/sn_test_taxonomy.md.)
pytestmark = pytest.mark.foundation


@pytest.fixture
def quad():
    return Quadrature.gauss_legendre(4)


@pytest.fixture
def slab_mesh():
    return Mesh1D(edges=np.linspace(0, 5, 11), mat_ids=np.zeros(10, dtype=int))


# ═══════════════════════════════════════════════════════════════════════
# Registry tests
# ═══════════════════════════════════════════════════════════════════════

class TestSNBCRegistry:
    """BOUNDARY_OPERATOR_REGISTRY is the single source of truth: resolves, validates, advertises."""

    def test_registry_keys(self):
        assert "vacuum" in SNMesh.BOUNDARY_OPERATOR_REGISTRY
        assert "reflective" in SNMesh.BOUNDARY_OPERATOR_REGISTRY

    def test_registry_docstrings(self):
        """Every factory has a docstring (used as description for UI query)."""
        for kind, factory in SNMesh.BOUNDARY_OPERATOR_REGISTRY.items():
            assert factory.__doc__ is not None, f"BC factory '{kind}' has no docstring"

    def test_registry_programmatic_query(self):
        """Descriptions are queryable via factory docstrings."""
        descriptions = {k: v.__doc__ for k, v in SNMesh.BOUNDARY_OPERATOR_REGISTRY.items()}
        assert "vacuum" in descriptions
        assert "reflective" in descriptions


# ═══════════════════════════════════════════════════════════════════════
# Resolution tests
# ═══════════════════════════════════════════════════════════════════════

class TestSNBCResolution:
    """BC resolution at SNMesh construction time."""

    def test_default_is_reflective(self, slab_mesh, quad):
        """None on mesh resolves to 'reflective' (eigenvalue default)."""
        sn = SNMesh(slab_mesh, quad, placeholder_materials())
        assert sn.bc["xmin"] == "reflective"
        assert sn.bc["xmax"] == "reflective"

    def test_explicit_vacuum(self, slab_mesh, quad):
        mesh = Mesh1D(
            edges=slab_mesh.edges, mat_ids=slab_mesh.mat_ids,
            bc_left=BC.vacuum, bc_right=BC.vacuum,
        )
        sn = SNMesh(mesh, quad, placeholder_materials())
        assert sn.bc["xmin"] == "vacuum"
        assert sn.bc["xmax"] == "vacuum"

    def test_mixed_bcs(self, slab_mesh, quad):
        mesh = Mesh1D(
            edges=slab_mesh.edges, mat_ids=slab_mesh.mat_ids,
            bc_left=BC.reflective, bc_right=BC.vacuum,
        )
        sn = SNMesh(mesh, quad, placeholder_materials())
        assert sn.bc["xmin"] == "reflective"
        assert sn.bc["xmax"] == "vacuum"

    def test_unknown_bc_raises(self, slab_mesh, quad):
        mesh = Mesh1D(
            edges=slab_mesh.edges, mat_ids=slab_mesh.mat_ids,
            bc_left=BC("white"),
        )
        with pytest.raises(ValueError, match="does not support.*'white'"):
            SNMesh(mesh, quad, placeholder_materials())

    def test_error_lists_supported(self, slab_mesh, quad):
        mesh = Mesh1D(
            edges=slab_mesh.edges, mat_ids=slab_mesh.mat_ids,
            bc_left=BC("periodic"),
        )
        with pytest.raises(ValueError, match="'reflective'.*'vacuum'"):
            SNMesh(mesh, quad, placeholder_materials())

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
            bc_xmin=BC.reflective, bc_xmax=BC.vacuum,
            bc_ymin=BC.reflective, bc_ymax=BC.vacuum,
        )
        sn = SNMesh(mesh, quad, placeholder_materials())
        assert sn.bc["xmin"] == "reflective"
        assert sn.bc["xmax"] == "vacuum"
        assert sn.bc["ymin"] == "reflective"
        assert sn.bc["ymax"] == "vacuum"

    def test_curvilinear_vacuum_resolves(self, quad):
        """Spherical/cylindrical accept vacuum; sweep enforces zero-incoming
        at the outer face. Vacuum support added in commits 655e3e5 / 37c5bbf
        (the curvilinear-only-reflective gate was removed and the inward
        sweep now branches on ``is_vacuum_outer``)."""
        mesh = Mesh1D(
            edges=np.linspace(0.1, 1.0, 6), mat_ids=np.zeros(5, dtype=int),
            coord=CoordSystem.SPHERICAL,
            bc_right=BC.vacuum,
        )
        sn = SNMesh(mesh, quad, placeholder_materials())
        assert sn.bc["xmax"] == "vacuum"


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
    from orpheus.geometry import (
        Mesh1D as _Mesh1D,
        Region,
        RegionMesh,
        StructuredGeometry,
    )

    case = get("sn_slab_2eg_1rg")
    mix = next(iter(case.materials.values()))
    materials = {0: mix}
    mesh_refl = _Mesh1D.from_geometry(
        StructuredGeometry(
            geometry="SLB",
            regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
            bcs=(BC.reflective, BC.reflective),
        ),
        region_meshes=(RegionMesh(n_cells=20),),
    )
    mesh_vac = Mesh1D(
        edges=mesh_refl.edges, mat_ids=mesh_refl.mat_ids,
        bc_left=BC.vacuum, bc_right=BC.vacuum,
    )
    # The fixture mixture is (n,2n)-free, so total production IS fission
    # production — asserted here so the hand formula below stays honest
    # if the reference case ever gains a Σ₂ channel.
    assert mix.Sig2.nnz == 0, (
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
        result_refl = solve_sn(materials, mesh_refl, quad)
        result_vac = solve_sn(materials, mesh_vac, quad)

        # Reflective has higher keff (no leakage vs leakage)
        assert result_refl.keff > result_vac.keff

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
                keff_tol=tol, flux_tol=tol, max_outer=200,
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
