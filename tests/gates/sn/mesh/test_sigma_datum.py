r"""The total cross section is a DATUM of the Problem — consumers campaign step 2,
sub-step C3a (O-5/O-6, 2026-09-13).

``MaterialMesh.sigma_t_cell`` is the per-cell :math:`\sigma_t` in the principled
``(ng, *spatial)`` layout: derived from the materials by default through the ONE
per-cell assembler (:func:`~orpheus.data.macro_xs.cell_xs.assemble_cell_xs`) and
REPLACED by :meth:`~orpheus.transport.mesh.material_mesh.MaterialMesh.with_cross_sections`
— a NEW Problem (another ``_identity_key``) over the SAME phase space (the same
``_contractibility_key``, so ``same_phase_space`` holds and fields pair).  The
hub's :attr:`mat_xs` is ONE :class:`MaterialXSField` per Problem, and its
``total_cross_section`` view IS the datum.  The σ-free geometry table (Stratum 1)
is shared by identity across σ-variants (the intern keys on contractibility) and
lives as long as a hub that used it.

Gate specification: ``scratch/_consumers/planning/test_architect_c3_delta.md``
§A.1–A.3, §B.3 (rows a1–a4, b1–b3, c1–c2, e1/e2/e4/e5 and the lifetime row).
Every ``[M]`` below is the test-architect's, taken pre-carve at ``93225c65``.
"""
from __future__ import annotations

import gc
import weakref

import numpy as np
import pytest

from orpheus.data.macro_xs.cell_xs import assemble_cell_xs
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry.mesh import BC, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.problem import SNProblem
from orpheus.transport.mesh.material_mesh import MaterialMesh

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    """``-O``-safe assertion (the canonical runner strips ``assert``)."""
    if not cond:
        raise AssertionError(msg)


def _mesh() -> Mesh1D:
    return Mesh1D(
        edges=np.linspace(0.0, 2.0, 9),
        mat_ids=np.array([0] * 4 + [1] * 4, dtype=int),
        bc_left=BC("reflective"), bc_right=BC("vacuum"),
    )


def _mats():
    return {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}


def _bare() -> MaterialMesh:
    return MaterialMesh(_mesh(), _mats())


def _sn() -> SNProblem:
    return SNProblem(_mesh(), Quadrature.gauss_legendre(4), _mats())


_TIERS = pytest.mark.parametrize("build", [_bare, _sn], ids=["MaterialMesh", "SNProblem"])


# ═══════════════════════════════════════════════════════════════════════
# A.1 — the σ datum enters ``_identity_key`` ONLY
# ═══════════════════════════════════════════════════════════════════════


class TestLawTheSigmaDatumIsIdentityOnly:
    """The datum moves the IDENTITY and not the PHASE SPACE — on both tiers."""

    @_TIERS
    def test_law_a_sigma_variant_is_a_DIFFERENT_problem(self, build) -> None:
        hub = build()
        hub_b = hub.with_cross_sections(3.0 * hub.sigma_t_cell)
        _require(hub_b != hub, "a σ-variant Problem must have another identity")
        _require(hash(hub_b) != hash(hub), "…and another hash (the key hashes the datum's bytes)")

    @_TIERS
    def test_law_a_sigma_variant_shares_the_PHASE_SPACE(self, build) -> None:
        hub = build()
        hub_b = hub.with_cross_sections(3.0 * hub.sigma_t_cell)
        _require(hub_b.same_phase_space(hub), (
            "the datum leaked into _contractibility_key — a depletion step is "
            "the SAME phase space (fields must pair; the geometry table must be shared)"
        ))

    @_TIERS
    def test_law_re_declaring_the_SAME_sigma_is_the_same_problem(self, build) -> None:
        """The NORMALISATION theorem: the ``total_cross_section`` view (a
        non-contiguous ``.T.reshape`` of the assembler's output) re-declared
        as the datum is the SAME Problem — same key, same hash."""
        hub = build()
        hub_c = hub.with_cross_sections(hub.mat_xs.total_cross_section)
        _require(hub_c == hub and hash(hub_c) == hash(hub), (
            "re-declaring the hub's own σ_t must be the same Problem — the datum "
            "is normalised (C-contiguous float64, -0.0 → +0.0) before it is keyed"
        ))

    @_TIERS
    def test_law_the_stored_datum_is_READ_ONLY(self, build) -> None:
        """The key is cached from the datum's bytes — an in-place write would
        desynchronise identity from content, so the datum refuses it."""
        hub = build()
        _require(hub.sigma_t_cell.flags.writeable is False, "the datum must be read-only")
        with pytest.raises(ValueError, match="read-only"):
            hub.sigma_t_cell[...] = 0.0

    def test_law_a_wrong_shape_is_REFUSED(self) -> None:
        hub = _sn()
        with pytest.raises(ValueError, match="sigma_t_cell must be"):
            hub.with_cross_sections(np.ones((hub.ng, hub.nx + 1)))


# ═══════════════════════════════════════════════════════════════════════
# The round-trip law the retired ``if __debug__`` assert used to state
# (``sn-cell-flatten-roundtrip``) — a REAL witness now, not an ``assert``
# the canonical ``-O`` runner strips.
# ═══════════════════════════════════════════════════════════════════════


class TestLawTheRoundTrip:
    @pytest.mark.verifies("sn-cell-flatten-roundtrip")
    @_TIERS
    def test_law_the_datum_is_the_assembled_cell_sigma_t(self, build) -> None:
        hub = build()
        xs = assemble_cell_xs(hub.materials, hub.mat_map)
        expected = xs.sig_t.T.reshape(hub.ng, *hub.spatial_shape)
        _require(np.array_equal(hub.sigma_t_cell, expected), (
            "the default datum must be the assembler's (N_cells, ng) output "
            "transposed into the principled (ng, *spatial) layout — bit-exactly"
        ))
        _require(np.array_equal(hub.mat_xs.total_cross_section, expected), (
            "the field's σ_t view must read the datum"
        ))


# ═══════════════════════════════════════════════════════════════════════
# A.2 / A.3 — ONE field per Problem; a σ-variant's field reads the override
# ═══════════════════════════════════════════════════════════════════════


class TestLawOneFieldPerProblem:
    def test_law_one_field_per_problem(self) -> None:
        """The COUNT row: every production consumer of ONE hub gets the same
        object — ``[M]`` pre-carve the two mints per entry
        (``solver.py`` + the hub's ``fission``) were 2 distinct objects."""
        from orpheus.sn.solver import solve_sn
        from orpheus.transport.mesh.material_xs_field import MaterialXSField

        hub = _sn()
        _require(hub.mat_xs is hub.mat_xs, "mat_xs is a cached_property")
        mints: list[int] = []
        original = MaterialXSField.from_mesh.__func__

        def spy(cls, mesh):  # noqa: ANN001, ANN202
            mints.append(id(mesh))
            return original(cls, mesh)

        MaterialXSField.from_mesh = classmethod(spy)  # type: ignore[method-assign]
        try:
            mats, mesh, quad = _mats(), _mesh(), Quadrature.gauss_legendre(4)
            solve_sn(mats, mesh, quad, max_outer=3)
        finally:
            MaterialXSField.from_mesh = classmethod(original)  # type: ignore[method-assign]
        _require(mints, "POSITIVE CONTROL: the spy must observe at least one mint")
        _require(len(mints) == 1, (
            f"ONE MaterialXSField per Problem: a k-solve minted {len(mints)} "
            "(pre-C3a: 2 — the solver's and the hub's fission)"
        ))

    def test_law_two_sigma_variants_have_DIFFERENT_fields(self) -> None:
        hub = _sn()
        hub_b = hub.with_cross_sections(3.0 * hub.sigma_t_cell)
        _require(hub_b.mat_xs is not hub.mat_xs, "a σ-variant owns its own field")
        _require(not np.array_equal(hub_b.mat_xs.total_cross_section, hub.mat_xs.total_cross_section),
                 "…whose σ_t view differs")

    def test_law_the_override_moves_sigma_t_and_NOTHING_else(self) -> None:
        """The honest read-set: σ_t and the DERIVED diffusion coefficient move
        (RULED 2026-09-13 fork 4 (a): ``D = 1/(3(σ_t − Σ_s1-outflow))`` per cell
        follows the datum, C3b-2); absorption, fission production and the
        emission spectrum re-derive from the MATERIALS and are unmoved."""
        hub = _sn()
        hub_b = hub.with_cross_sections(3.0 * hub.sigma_t_cell)
        _require(np.array_equal(hub_b.mat_xs.total_cross_section, 3.0 * hub.mat_xs.total_cross_section),
                 "σ_t: the field reads the hub's OWN datum (a scalar multiply is bit-exact)")
        _require(not np.array_equal(hub_b.mat_xs.diffusion_coefficient, hub.mat_xs.diffusion_coefficient),
                 "D follows the overridden σ_t (the derived coefficient moves)")
        for view in ("absorption_cross_section", "fission_production", "emission_spectrum"):
            _require(np.array_equal(getattr(hub_b.mat_xs, view), getattr(hub.mat_xs, view)),
                     f"{view} must be unmoved by a σ_t override (it derives from the materials)")

    def test_law_the_operators_posed_over_the_variant_MOVE(self) -> None:
        """Activation: the datum is USED — the k-solve over the ×3 σ-variant differs.

        ⚠ Spelled at the SOLVE level on purpose: the raw-data entry
        (``solve_sn(materials, mesh, quad)``) rebuilds a hub from the materials
        and cannot see a σ-variant, and ``np.asarray`` of a coupled field is an
        OBJECT array on which ``array_equal`` is always False — ``[M]`` a first
        draft compared applies that way and was green under the
        `with_cross_sections`-returns-`self` arm (a tautology, caught by the
        C3a battery)."""
        from orpheus.numerics.eigenvalue import power_iteration
        from orpheus.sn.solver import SNSolver

        hub = _sn()
        hub_b = hub.with_cross_sections(3.0 * hub.sigma_t_cell)
        ks = []
        for h in (hub, hub_b):
            solver = SNSolver(h, inner_solver="source_iteration", keff_tol=1e-9, flux_tol=1e-8)
            ks.append(float(power_iteration(solver, max_iter=300, budget_name="max_outer").keff))
        _require(abs(ks[1] - ks[0]) / abs(ks[0]) > 1e-3, f"a ×3 σ-variant must move k: {ks[0]!r} vs {ks[1]!r}")


# ═══════════════════════════════════════════════════════════════════════
# B.3 — the geometry intern under σ-variants
# ═══════════════════════════════════════════════════════════════════════


class TestLawTheGeometryInternUnderSigmaVariants:
    def test_law_two_sigma_variants_SHARE_the_geometry_table(self) -> None:
        from orpheus.sn.loss_representation import geometry_cache_for

        hub = _sn()
        hub_b = hub.with_cross_sections(3.0 * hub.sigma_t_cell)
        g, g_b = geometry_cache_for(hub, hub.angular_closure), geometry_cache_for(hub_b, hub_b.angular_closure)
        _require(g_b is g, (
            "the σ-free Stratum-1 table is the PHASE SPACE's: a σ-variant must "
            "read the same object (the intern keys on contractibility, not identity)"
        ))

    def test_law_the_table_is_built_ONCE_per_solve(self) -> None:
        from orpheus.sn.loss_representation import _GEOM_CACHE_INTERN
        from orpheus.sn.solver import solve_sn
        from orpheus.sn.sweep.cache import StreamingCoefficientCache

        builds: list[int] = []
        original = StreamingCoefficientCache.from_mesh_and_quad.__func__

        def spy(cls, mesh):  # noqa: ANN001, ANN202
            builds.append(1)
            return original(cls, mesh)

        _GEOM_CACHE_INTERN.clear()
        StreamingCoefficientCache.from_mesh_and_quad = classmethod(spy)  # type: ignore[method-assign]
        try:
            solve_sn(_mats(), _mesh(), Quadrature.gauss_legendre(4), max_outer=3)
        finally:
            StreamingCoefficientCache.from_mesh_and_quad = classmethod(original)  # type: ignore[method-assign]
        _require(builds, "POSITIVE CONTROL: the spy must observe at least one build")
        _require(len(builds) == 1, (
            f"Stratum 1 built {len(builds)}x in one solve — a weak-valued intern "
            "with no holder rebuilds per sweep ([M] 550 per slab eigen solve)"
        ))

    def test_law_the_interned_VALUE_retains_no_hub(self) -> None:
        from orpheus.sn.loss_representation import geometry_cache_for

        hub = _sn()
        table = geometry_cache_for(hub, hub.angular_closure)
        for name in getattr(type(table), "__slots__", ()):
            if name == "__weakref__":
                continue
            value = getattr(table, name)
            _require(value is not hub and not isinstance(value, (SNProblem, Mesh1D)),
                     f"the interned table must hold no hub reference (field {name!r})")

    def test_law_the_interned_table_dies_with_its_last_holder(self) -> None:
        """The LIFETIME claim the re-key changes — a MUST-STAY-GREEN pin:
        the table is weak-referenceable and evicts when every holder is gone."""
        from orpheus.sn.loss_representation import _GEOM_CACHE_INTERN, geometry_cache_for

        _GEOM_CACHE_INTERN.clear()
        hub = _sn()
        table = geometry_cache_for(hub, hub.angular_closure)
        ref = weakref.ref(table)
        _require(len(_GEOM_CACHE_INTERN) == 1, "one live table")
        del table, hub
        gc.collect()
        _require(ref() is None, "the table must die with its last holder (the hub)")
        _require(len(_GEOM_CACHE_INTERN) == 0, "…and leave the intern empty")

    def test_law_the_intern_is_BOUNDED_by_distinct_phase_spaces(self) -> None:
        from orpheus.sn.loss_representation import _GEOM_CACHE_INTERN, geometry_cache_for

        _GEOM_CACHE_INTERN.clear()
        hub = _sn()
        hubs = [hub.with_cross_sections((1.0 + 0.1 * k) * hub.sigma_t_cell) for k in range(5)]
        tables = [geometry_cache_for(h, h.angular_closure) for h in [hub, *hubs]]
        _require(len(_GEOM_CACHE_INTERN) == 1, (
            f"6 σ-variants of ONE phase space must be ONE intern entry; got {len(_GEOM_CACHE_INTERN)}"
        ))
        _require(all(t is tables[0] for t in tables), "…all reading one table")
