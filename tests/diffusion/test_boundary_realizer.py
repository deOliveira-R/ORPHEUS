r"""DiffusionBoundaryRealizer — the albedo-family realization (#290 P3).

This file FLIPPED from the Wave-5 stub-invariant tests (the realizer
raised ``NotImplementedError`` by design) to positive realization
tests when the functional body shipped (#290 P3, closes #182). The
identity pins (``method_name``; :class:`BoundaryRealizer` Protocol
conformance — the successor of the registry-lookup pin after the
#290 P7b registry dissolution) open the file; everything else pins
the now-real physics:

* the law → 𝒜 table of :math:`J^- = \mathcal{A}\,J^+` (#290 ruling 2),
* the ruling-3 semantics — vacuum IS Marshak zero-incoming
  (:math:`J^- = 0`), zero-flux IS the Dirichlet idealization
  (:math:`\phi_\Gamma = 2(J^+ + J^-) = 0` under P1, 𝒜 = −1),
* the structure-keyed collapses (Zero / Identity / 𝒜·Identity)
  mirroring the SN realizer's albedo branch,
* the ``BlockRole.BOUNDARY`` stamp on every realized operator,
* the two refusals (periodic: trace-block wrap, not a per-face
  albedo; prescribed inflow: affine source, not an operator),
* rank-N composition through the descriptor-tree walker with the
  diffusion realizer at the leaves (the #290 P3 walker
  generalization), and
* the :class:`DiffusionMethodSpace` construction surface (mirror of
  ``SNMethodSpace.minimal`` / ``for_face``).

Convention contract: ``.claude/plans/diffusion_crosswalk.md`` (the
BC-family table these tests transcribe — vv L11: production derives
the map through ``_partial_current_albedo``, the tests hand-list it).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.diffusion.augmented_mesh import DiffusionMesh
from orpheus.diffusion.boundary_realizer import DiffusionBoundaryRealizer
from orpheus.diffusion.method_space import DiffusionMethodSpace
from orpheus.geometry.boundary import (
    AlbedoBoundary,
    BoundaryError,
    BoundaryRealizer,
    PeriodicBoundary,
    PrescribedInflow,
    ReflectiveBoundary,
    VacuumInflow,
    WhiteBoundary,
    ZeroFluxBoundary,
    realize_recursively,
)
from orpheus.numerics.operator import (
    BlockRole,
    BoundaryOperator,
    IdentityOperator,
    ZeroOperator,
)
from orpheus.geometry.mesh import Mesh1D

# V&V LEVELS — B0.3 REPAIR (2026-07-30). This file used to carry a
# module-level ``pytestmark = [pytest.mark.foundation]`` blanket while
# ``TestRulingThreeSemantics`` also carried ``@pytest.mark.l0``. The
# two stack, so all FOUR tests in that class resolved with conflicting
# level markers ``['foundation', 'l0']`` and tripped
# ``PytestUnknownMarkWarning`` out of ``tests/conftest.py:120``
# ("using 'l0'" — the physics ladder dominates by design).
#
# The blanket is therefore gone and each class states its own level:
# ``foundation`` for the structural / registration / dispatch pins (no
# theory-page ``:label:``), ``l0`` for the ruling-3 semantics class,
# whose four tests ARE term-level physics claims against the P1
# partial-current dictionary (J⁻ = 0 for vacuum, φ_Γ = 2(J⁺+J⁻) = 0
# for zero-flux, J⁺ − J⁻ = 0 for reflective, J⁻ ≥ 0 for the physical
# sub-family). Resolution is unchanged for every test — only the
# conflict is gone.


@pytest.fixture
def realizer() -> DiffusionBoundaryRealizer:
    return DiffusionBoundaryRealizer()


@pytest.fixture
def minimal_ms() -> DiffusionMethodSpace:
    return DiffusionMethodSpace.minimal()


@pytest.fixture
def outflow_probe() -> np.ndarray:
    """A J⁺ face-slot half: shape ``(ng, *face_spatial)`` = (2, 3),
    strictly positive (a physical outflow partial current)."""
    return np.array([[1.7, 0.3, 2.1], [0.2, 0.9, 0.4]])


# ═══════════════════════════════════════════════════════════════════════
# Registration — the pins that survive from the stub era, now positive
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestRegistration:
    def test_method_name_attribute(self):
        assert DiffusionBoundaryRealizer.method_name == "diffusion"

    def test_conforms_to_the_boundary_realizer_protocol(self):
        """The structural pin that replaced the registry-lookup pin
        when #290 P7b dissolved ``BoundaryRealizerRegistry``: the
        realizer satisfies the :class:`BoundaryRealizer` Protocol the
        walker and the method-mesh hook dispatch through."""
        assert isinstance(DiffusionBoundaryRealizer(), BoundaryRealizer)

    def test_realize_is_functional(self, realizer, minimal_ms):
        """The stub-era ``NotImplementedError`` contract is GONE —
        realizing a law returns a working 1-arg operator (#182)."""
        op = realizer.realize(VacuumInflow(), minimal_ms)
        assert op.apply(np.ones((2,))) is not None


# ═══════════════════════════════════════════════════════════════════════
# The albedo table — law → 𝒜 in J⁻ = 𝒜·J⁺ (the crosswalk BC-family
# table hand-transcribed; production derives it via
# ``_partial_current_albedo`` — vv L11 structural independence)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestAlbedoTable:
    @pytest.mark.parametrize(
        ("law", "expected_albedo"),
        [
            pytest.param(VacuumInflow(), 0.0, id="vacuum->0"),
            pytest.param(
                ReflectiveBoundary(axis="x"), 1.0, id="reflective->1",
            ),
            pytest.param(
                WhiteBoundary(axis="x", outward_sign=+1), 1.0,
                id="white->1",
            ),
            pytest.param(
                AlbedoBoundary(albedo=0.37), 0.37, id="albedo(0.37)",
            ),
            pytest.param(ZeroFluxBoundary(), -1.0, id="zero_flux->-1"),
            pytest.param(
                ReflectiveBoundary(axis="x", albedo=0.7), 0.7,
                id="reflective(0.7)",
            ),
            pytest.param(
                WhiteBoundary(axis="x", outward_sign=+1, albedo=0.8),
                0.8, id="white(0.8)",
            ),
        ],
    )
    def test_realized_action_is_albedo_multiple(
        self, realizer, minimal_ms, outflow_probe, law, expected_albedo,
    ):
        """``op.apply(J⁺) == 𝒜·J⁺`` exactly (a single scalar multiply
        — bit-identical to the hand-computed product)."""
        op = realizer.realize(law, minimal_ms)
        np.testing.assert_array_equal(
            op.apply(outflow_probe), expected_albedo * outflow_probe,
        )

    def test_white_coincides_with_reflective_at_p1(
        self, realizer, minimal_ms, outflow_probe,
    ):
        """The documented P1 coincidence: specular vs Lambertian
        return differ only in angular redistribution, which the
        half-range ℓ=0 moments integrate out — identical realized
        action for equal return amplitude α."""
        specular = realizer.realize(
            ReflectiveBoundary(axis="x", albedo=0.6), minimal_ms,
        )
        lambertian = realizer.realize(
            WhiteBoundary(axis="x", outward_sign=+1, albedo=0.6),
            minimal_ms,
        )
        np.testing.assert_array_equal(
            specular.apply(outflow_probe), lambertian.apply(outflow_probe),
        )

    @pytest.mark.parametrize(
        ("law", "expected_type"),
        [
            pytest.param(VacuumInflow(), ZeroOperator, id="vacuum->Zero"),
            pytest.param(
                AlbedoBoundary(albedo=0.0), ZeroOperator, id="albedo(0)->Zero",
            ),
            pytest.param(
                ReflectiveBoundary(axis="x"), IdentityOperator,
                id="reflective->Identity",
            ),
            pytest.param(
                WhiteBoundary(axis="x", outward_sign=+1), IdentityOperator,
                id="white->Identity",
            ),
            pytest.param(
                AlbedoBoundary(albedo=1.0), IdentityOperator,
                id="albedo(1)->Identity",
            ),
        ],
    )
    def test_structure_keyed_collapse(
        self, realizer, minimal_ms, law, expected_type,
    ):
        """𝒜 ∈ {0, 1} collapses to the structural primitive (mirror of
        the SN realizer's AlbedoBoundary branch)."""
        assert isinstance(realizer.realize(law, minimal_ms), expected_type)

    def test_every_realized_operator_is_boundary_block(
        self, realizer, minimal_ms,
    ):
        """The A_ss stamp (#208 / Wave O): every realized law carries
        ``BlockRole.BOUNDARY`` and satisfies the isinstance marker."""
        laws = [
            VacuumInflow(),
            ReflectiveBoundary(axis="x"),
            WhiteBoundary(axis="x", outward_sign=+1),
            AlbedoBoundary(albedo=0.37),
            ZeroFluxBoundary(),
        ]
        for law in laws:
            op = realizer.realize(law, minimal_ms)
            assert op.block_role is BlockRole.BOUNDARY, type(law).__name__
            assert isinstance(op, BoundaryOperator), type(law).__name__


# ═══════════════════════════════════════════════════════════════════════
# Ruling-3 semantics — the honest-naming laws mean what they say
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l0
class TestRulingThreeSemantics:
    def test_vacuum_means_zero_incoming_current(
        self, realizer, minimal_ms, outflow_probe,
    ):
        """Vacuum IS the Marshak condition: J⁻ = 0 exactly, whatever
        leaks out (#290 ruling 3 — NOT the zero-flux Dirichlet the
        legacy island mis-registered under this name)."""
        inflow = realizer.realize(VacuumInflow(), minimal_ms).apply(
            outflow_probe,
        )
        np.testing.assert_array_equal(inflow, np.zeros_like(outflow_probe))

    def test_zero_flux_means_zero_boundary_scalar_flux(
        self, realizer, minimal_ms, outflow_probe,
    ):
        """Zero-flux delivers its defining property under the P1
        dictionary: φ_Γ = 2(J⁺ + J⁻) = 0 EXACTLY (𝒜 = −1 ⟹
        J⁻ = −J⁺ ⟹ the pair sums to zero bit-exactly)."""
        inflow = realizer.realize(ZeroFluxBoundary(), minimal_ms).apply(
            outflow_probe,
        )
        boundary_scalar_flux = 2.0 * (outflow_probe + inflow)
        np.testing.assert_array_equal(
            boundary_scalar_flux, np.zeros_like(outflow_probe),
        )

    def test_reflective_means_zero_net_current(
        self, realizer, minimal_ms, outflow_probe,
    ):
        """Reflective: J = J⁺ − J⁻ = 0 exactly (perfect return)."""
        inflow = realizer.realize(
            ReflectiveBoundary(axis="x"), minimal_ms,
        ).apply(outflow_probe)
        net_current = outflow_probe - inflow
        np.testing.assert_array_equal(
            net_current, np.zeros_like(outflow_probe),
        )

    def test_positivity_is_a_physical_law_property(
        self, realizer, minimal_ms, outflow_probe,
    ):
        """J⁻ ≥ 0 for every PHYSICAL law (𝒜 ∈ [0, 1]) given J⁺ ≥ 0;
        zero-flux (𝒜 = −1) deliberately violates it — positivity is a
        property of the physical sub-family, not a type invariant
        (crosswalk BC-family table)."""
        physical = [
            VacuumInflow(),
            ReflectiveBoundary(axis="x"),
            WhiteBoundary(axis="x", outward_sign=+1),
            AlbedoBoundary(albedo=0.37),
        ]
        for law in physical:
            inflow = realizer.realize(law, minimal_ms).apply(outflow_probe)
            assert np.all(inflow >= 0.0), type(law).__name__
        dirichlet_inflow = realizer.realize(
            ZeroFluxBoundary(), minimal_ms,
        ).apply(outflow_probe)
        assert np.all(dirichlet_inflow < 0.0)


# ═══════════════════════════════════════════════════════════════════════
# Refusals — the laws that are NOT a per-face rank-1 albedo
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestRefusals:
    def test_periodic_refused(self, realizer, minimal_ms):
        """Periodic couples a face to its OPPOSITE face — a
        trace-block permutation, not a per-face albedo. Refused
        loudly rather than realized as the wrong (identity) thing."""
        with pytest.raises(BoundaryError, match="OPPOSITE face"):
            realizer.realize(PeriodicBoundary(), minimal_ms)

    def test_prescribed_inflow_refused(self, realizer, minimal_ms):
        """Prescribed inflow is the rank-0 AFFINE law — a boundary
        SOURCE, not a linear boundary operator B (the same
        operator/source split SN keeps un-stamped)."""
        with pytest.raises(BoundaryError, match="AFFINE"):
            realizer.realize(PrescribedInflow(), minimal_ms)

    def test_unknown_law_type_refused_with_case_list(
        self, realizer, minimal_ms,
    ):
        """A law outside the dispatch table names itself and the
        realizable cases (mirror of the SN fallthrough contract)."""

        class _NotALaw:
            pass

        with pytest.raises(BoundaryError, match="_NotALaw"):
            realizer.realize(_NotALaw(), minimal_ms)  # type: ignore[arg-type]


# ═══════════════════════════════════════════════════════════════════════
# Rank-N composition — the walker generalized over leaf realizers
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestComposition:
    def test_law_sum_realizes_to_operator_sum_action(
        self, realizer, minimal_ms, outflow_probe,
    ):
        """A Marshak-style mix realizes leaf-by-leaf through the
        diffusion realizer: (0.3·reflective + 0.7·albedo(0.5)) acts
        as the same weighted sum of realized leaf actions."""
        tree = (
            0.3 * ReflectiveBoundary(axis="x")
            + 0.7 * AlbedoBoundary(albedo=0.5)
        )
        composed = realize_recursively(tree, minimal_ms, realizer=realizer)
        expected = 0.3 * outflow_probe + 0.7 * (0.5 * outflow_probe)
        np.testing.assert_array_equal(composed.apply(outflow_probe), expected)

    def test_scaled_zero_flux_leaf(
        self, realizer, minimal_ms, outflow_probe,
    ):
        """A LawScaled root over the new law: 0.5·zero_flux acts as
        0.5·(−J⁺)."""
        composed = realize_recursively(
            0.5 * ZeroFluxBoundary(), minimal_ms, realizer=realizer,
        )
        np.testing.assert_array_equal(
            composed.apply(outflow_probe),
            0.5 * (-1.0 * outflow_probe),
        )

    def test_composition_refusal_propagates(self, realizer, minimal_ms):
        """A refused leaf inside a composition tree refuses the whole
        realization — no silent partial tree."""
        tree = 0.5 * PeriodicBoundary() + 0.5 * VacuumInflow()
        with pytest.raises(BoundaryError, match="OPPOSITE face"):
            realize_recursively(tree, minimal_ms, realizer=realizer)


# ═══════════════════════════════════════════════════════════════════════
# DiffusionMethodSpace — the mirrored construction surface
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestDiffusionMethodSpace:
    @staticmethod
    def _mesh() -> DiffusionMesh:
        mesh1d = Mesh1D(np.linspace(0.0, 10.0, 5), np.zeros(4, dtype=int))
        return DiffusionMesh(mesh1d, {0: get_mixture("A", "2g")})

    def test_minimal_is_metadata_free(self):
        ms = DiffusionMethodSpace.minimal()
        assert ms.mesh is None and ms.face is None and ms.trace is None

    def test_for_face_reads_the_mesh_trace(self):
        """#290 P7a: for_face reads the trace OFF the mesh — the single
        source; the stored trace IS the mesh's cached scalar trace."""
        mesh = self._mesh()
        ms = DiffusionMethodSpace.for_face(mesh=mesh, face="xmax")
        assert ms.face == "xmax"
        assert ms.trace is mesh.scalar_trace
        assert ms.mesh is mesh

    def test_for_face_validates_face_against_trace_inventory(self):
        """A method space for a nonexistent face is unrepresentable —
        the mirror of SN's for_face deriving (and thereby validating)
        per-face indices from its trace."""
        with pytest.raises(ValueError, match="ymax"):
            DiffusionMethodSpace.for_face(mesh=self._mesh(), face="ymax")

    def test_frozen(self):
        ms = DiffusionMethodSpace.minimal()
        with pytest.raises(AttributeError):
            ms.face = "xmin"  # type: ignore[misc]
