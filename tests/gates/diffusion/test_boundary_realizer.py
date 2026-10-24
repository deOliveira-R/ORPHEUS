r"""The diffusion realizer table — the albedo-family realization (#290 P3; the table of wave 2).

This file FLIPPED from the Wave-5 stub-invariant tests to positive realization tests when the functional body shipped (#290 P3,
closes #182), and was re-posed on the typed vocabulary by the wave-2 carve (slice S3): the realizer is the table
:data:`~orpheus.diffusion.boundary_realizer.DIFFUSION_REALIZERS` (the class ``DiffusionBoundaryRealizer`` and its
``_partial_current_albedo`` are retired), and a declaration is a ``WallCondition`` (or the method-only ``ZeroFluxCondition``, or a
deck). The structural facts of the table (the six registered classes, the tiers, the pin rows, the functor's sum) are the gates of
``test_diffusion_realizer_table.py`` and ``tests/gates/transport/test_admission_matrix.py``; this file keeps the PHYSICS the
realized operators carry, read through the door a problem uses:

* the condition → 𝒜 table of :math:`J^- = \mathcal{A}\,J^+` (#290 ruling 2),
* the ruling-3 semantics — vacuum IS Marshak zero-incoming (:math:`J^- = 0`), zero-flux IS the Dirichlet idealization
  (:math:`\phi_\Gamma = 2(J^+ + J^-) = 0` under P1, 𝒜 = −1),
* the structure-keyed collapses (Zero / Identity / 𝒜·Identity) mirroring the SN realizer's scalar branch,
* the ``BlockRole.BOUNDARY`` stamp on every realized operator,
* the refusals (a translation deck: trace-block wrap, not a per-face albedo; a boundary source: affine, not an operator),
* the sum of a mixed wall (the functor), and
* the :class:`DiffusionMethodSpace` construction surface (mirror of ``SNMethodSpace.minimal`` / ``for_face``).

Retired with the walker (the descriptor-tree ``realize_recursively``) and the law ladder: the composition-tree tests (the sum of
the chains' scalars is the functor, gated in ``test_diffusion_realizer_table.py``), the "unknown law type" fallthrough (a
declaration that is not a deck, a condition or an inflow is a ``TypeError`` at the table's door, pinned below) and the
"unpopulated factors" law (a ``BoundaryCondition`` subclass that omits an accessor cannot be constructed, pinned below).

Convention contract: ``.claude/plans/diffusion_crosswalk.md`` (the BC-family table these tests transcribe — vv L11: production
derives the map through the realizer classes, the tests hand-list it).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.diffusion.augmented_mesh import DiffusionMesh
from orpheus.diffusion.boundary_realizer import DIFFUSION_REALIZERS
from orpheus.diffusion.method_space import DiffusionMethodSpace
from orpheus.diffusion.zero_flux import ZeroFluxCondition
from orpheus.geometry import CoordSystem, StructuredGeometry
from orpheus.geometry.deck import MirrorDeck, TranslationDeck
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.operator import BlockRole, BoundaryOperator, IdentityOperator, ZeroOperator
from orpheus.transport.boundary_conditions import BoundaryCondition, BoundaryRefused, FaceSite, RefusalAxis, WallCondition

# V&V LEVELS — B0.3 REPAIR (2026-07-30). This file used to carry a module-level ``pytestmark = [pytest.mark.foundation]`` blanket
# while ``TestRulingThreeSemantics`` also carried ``@pytest.mark.l0``. The two stack, so all FOUR tests in that class resolved with
# conflicting level markers ``['foundation', 'l0']`` and tripped ``PytestUnknownMarkWarning``. The blanket is therefore gone and each
# class states its own level: ``foundation`` for the structural pins, ``l0`` for the ruling-3 semantics class, whose tests ARE
# term-level physics claims against the P1 partial-current dictionary (J⁻ = 0 for vacuum, φ_Γ = 2(J⁺+J⁻) = 0 for zero-flux,
# J⁺ − J⁻ = 0 for a mirror, J⁻ ≥ 0 for the physical sub-family).

_SITE = FaceSite("xmax", CoordSystem.CARTESIAN, 1)


@pytest.fixture
def minimal_ms() -> DiffusionMethodSpace:
    return DiffusionMethodSpace.minimal()


@pytest.fixture
def outflow_probe() -> np.ndarray:
    """A J⁺ face-slot half: shape ``(ng, *face_spatial)`` = (2, 3), strictly positive (a physical outflow partial current)."""
    return np.array([[1.7, 0.3, 2.1], [0.2, 0.9, 0.4]])


@pytest.fixture
def realize(minimal_ms):
    """``realize(declared)``: the table's operator for a declaration on a Cartesian face (the door a problem uses)."""
    return lambda declared: DIFFUSION_REALIZERS.realize(declared, _SITE, minimal_ms)


# ═══════════════════════════════════════════════════════════════════════
# Registration — the pins that survive from the stub era, now positive
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestRegistration:
    def test_method_name_attribute(self):
        assert DIFFUSION_REALIZERS.method == "diffusion"

    def test_the_mesh_holds_the_table_it_realizes_through(self):
        """The structural successor of the registry-lookup pin (#290 P7b dissolved ``BoundaryRealizerRegistry``): you hold the
        method-mesh, therefore you hold its table."""
        assert DiffusionMesh.boundary_realizers is DIFFUSION_REALIZERS

    def test_realize_is_functional(self, realize):
        """The stub-era ``NotImplementedError`` contract is GONE — realizing a condition returns a working 1-arg operator (#182)."""
        op = realize(WallCondition.vacuum())
        assert op.apply(np.ones((2,))) is not None


# ═══════════════════════════════════════════════════════════════════════
# The albedo table — condition → 𝒜 in J⁻ = 𝒜·J⁺ (the crosswalk BC-family
# table hand-transcribed; production derives it through the realizers —
# vv L11 structural independence)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestAlbedoTable:
    @pytest.mark.parametrize(
        ("declared", "expected_albedo"),
        [
            pytest.param(WallCondition.vacuum(), 0.0, id="vacuum->0"),
            pytest.param(MirrorDeck(), 1.0, id="mirror->1"),
            pytest.param(WallCondition.white(1.0), 1.0, id="white->1"),
            pytest.param(WallCondition.white(0.37), 0.37, id="albedo(0.37)"),
            pytest.param(ZeroFluxCondition(), -1.0, id="zero_flux->-1"),
            pytest.param(WallCondition.specular(0.7), 0.7, id="specular(0.7)"),
            pytest.param(WallCondition.white(0.8), 0.8, id="white(0.8)"),
        ],
    )
    def test_realized_action_is_albedo_multiple(self, realize, outflow_probe, declared, expected_albedo):
        """``op.apply(J⁺) == 𝒜·J⁺`` exactly (a single scalar multiply — bit-identical to the hand-computed product)."""
        np.testing.assert_array_equal(realize(declared).apply(outflow_probe), expected_albedo * outflow_probe)

    def test_white_coincides_with_specular_at_p1(self, realize, outflow_probe):
        """The documented P1 coincidence (the tier-2 witness): specular vs Lambertian return differ only in angular redistribution,
        which the half-range ℓ=0 moments integrate out — identical realized action for equal return amplitude α."""
        np.testing.assert_array_equal(
            realize(WallCondition.specular(0.6)).apply(outflow_probe), realize(WallCondition.white(0.6)).apply(outflow_probe),
        )

    @pytest.mark.parametrize(
        ("declared", "expected_type"),
        [
            pytest.param(WallCondition.vacuum(), ZeroOperator, id="vacuum->Zero"),
            pytest.param(WallCondition.white(0.0), ZeroOperator, id="albedo(0)->Zero"),
            pytest.param(MirrorDeck(), IdentityOperator, id="mirror->Identity"),
            pytest.param(WallCondition.white(1.0), IdentityOperator, id="white->Identity"),
            pytest.param(WallCondition.specular(1.0), IdentityOperator, id="specular(1)->Identity"),
        ],
    )
    def test_structure_keyed_collapse(self, realize, declared, expected_type):
        """𝒜 ∈ {0, 1} collapses to the structural primitive (mirror of the SN realizer's scalar branch)."""
        assert isinstance(realize(declared), expected_type)

    def test_every_realized_operator_is_boundary_block(self, realize):
        """The A_ss stamp (#208 / Wave O): every realized condition carries ``BlockRole.BOUNDARY`` and satisfies the isinstance marker."""
        for declared in (WallCondition.vacuum(), MirrorDeck(), WallCondition.white(1.0), WallCondition.white(0.37), ZeroFluxCondition()):
            op = realize(declared)
            assert op.block_role is BlockRole.BOUNDARY, type(declared).__name__
            assert isinstance(op, BoundaryOperator), type(declared).__name__


# ═══════════════════════════════════════════════════════════════════════
# Ruling-3 semantics — the honest-naming conditions mean what they say
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l0
class TestRulingThreeSemantics:
    def test_vacuum_means_zero_incoming_current(self, realize, outflow_probe):
        """Vacuum IS the Marshak condition: J⁻ = 0 exactly, whatever leaks out (#290 ruling 3 — NOT the zero-flux Dirichlet the
        legacy island mis-registered under this name)."""
        inflow = realize(WallCondition.vacuum()).apply(outflow_probe)
        np.testing.assert_array_equal(inflow, np.zeros_like(outflow_probe))

    def test_zero_flux_means_zero_boundary_scalar_flux(self, realize, outflow_probe):
        """Zero-flux delivers its defining property under the P1 dictionary: φ_Γ = 2(J⁺ + J⁻) = 0 EXACTLY (𝒜 = −1 ⟹ J⁻ = −J⁺ ⟹ the
        pair sums to zero bit-exactly)."""
        inflow = realize(ZeroFluxCondition()).apply(outflow_probe)
        boundary_scalar_flux = 2.0 * (outflow_probe + inflow)
        np.testing.assert_array_equal(boundary_scalar_flux, np.zeros_like(outflow_probe))

    def test_a_mirror_means_zero_net_current(self, realize, outflow_probe):
        """A mirror: J = J⁺ − J⁻ = 0 exactly (perfect return)."""
        inflow = realize(MirrorDeck()).apply(outflow_probe)
        np.testing.assert_array_equal(outflow_probe - inflow, np.zeros_like(outflow_probe))

    def test_positivity_is_a_physical_property(self, realize, outflow_probe):
        """J⁻ ≥ 0 for every PHYSICAL condition (𝒜 ∈ [0, 1]) given J⁺ ≥ 0; zero-flux (𝒜 = −1) deliberately violates it — positivity is
        a property of the physical sub-family, not a type invariant (crosswalk BC-family table)."""
        for declared in (WallCondition.vacuum(), MirrorDeck(), WallCondition.white(1.0), WallCondition.white(0.37)):
            assert np.all(realize(declared).apply(outflow_probe) >= 0.0), type(declared).__name__
        assert np.all(realize(ZeroFluxCondition()).apply(outflow_probe) < 0.0)


# ═══════════════════════════════════════════════════════════════════════
# Refusals — the declarations that are NOT a per-face rank-1 albedo
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestRefusals:
    def test_a_translation_deck_is_refused(self, realize):
        """A translation deck couples a face to its OPPOSITE face — a trace-block permutation, not a per-face albedo. Refused
        loudly (``TOPOLOGY``) rather than realized as the wrong (identity) thing."""
        with pytest.raises(BoundaryRefused, match="opposite face") as info:
            realize(TranslationDeck())
        assert info.value.admission.axis is RefusalAxis.TOPOLOGY

    def test_a_source_is_refused(self, realize):
        """A boundary source is the AFFINE part of the trace map — a boundary SOURCE, not a linear boundary operator B (the same
        operator/source split SN keeps). Refused as plumbing, #656."""
        with pytest.raises(BoundaryRefused) as info:
            realize(WallCondition.uniform(1.0))
        assert info.value.admission.axis is RefusalAxis.PLUMBING and info.value.admission.issue == 656

    def test_a_non_declaration_is_refused_at_the_door(self, minimal_ms):
        """What the "unknown law type" leg of the law ladder pinned (mirror of the SN fallthrough contract): an object that is not a
        deck, a condition or an inflow is a ``TypeError`` at the table's door, naming what it takes, before any realizer runs."""

        class _NotACondition:
            pass

        with pytest.raises(TypeError, match="DeckTransformation, a BoundaryCondition or an InflowTerms"):
            DIFFUSION_REALIZERS.realize(_NotACondition(), _SITE, minimal_ms)  # type: ignore[arg-type]

    def test_a_condition_without_its_accessors_cannot_be_constructed(self):
        """The successor of the "unpopulated factors" leg: a :class:`BoundaryCondition` subclass shipping without its three
        accessors (``return_fraction``, ``source``, ``face_local``) is abstract, so it is refused at construction by name and no
        table is ever reached."""

        class _ConditionWithoutAccessors(BoundaryCondition):
            pass

        with pytest.raises(TypeError, match="_ConditionWithoutAccessors"):
            _ConditionWithoutAccessors()  # type: ignore[abstract]


# ═══════════════════════════════════════════════════════════════════════
# The sum of a mixed wall — the functor on the response network
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestComposition:
    def test_a_mixed_wall_realizes_to_the_sum_of_the_scalars_action(self, realize, outflow_probe):
        """A Marshak-style mix realizes chain by chain: 0.3·specular + 0.5·white acts as the sum of the two realized chain actions,
        𝒜 = 0.3 + 0.5 (the functor, plan law L-L; its own gate pins the order and the single-chain images bit for bit)."""
        mixed = realize(WallCondition.from_amplitudes(0.3, 0.5)).apply(outflow_probe)
        np.testing.assert_array_almost_equal_nulp(mixed, 0.3 * outflow_probe + 0.5 * outflow_probe, nulp=2)

    def test_a_scaled_zero_flux_does_not_exist(self):
        """The ``0.5 * ZeroFluxBoundary()`` leaf of the old walker is unspellable: a method-only condition defines no arithmetic
        (law L-F), so the scaled Dirichlet condition is Python's own ``TypeError``."""
        with pytest.raises(TypeError):
            0.5 * ZeroFluxCondition()  # type: ignore[operator]


# ═══════════════════════════════════════════════════════════════════════
# DiffusionMethodSpace — the mirrored construction surface
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestDiffusionMethodSpace:
    @staticmethod
    def _mesh() -> DiffusionMesh:
        geometry = StructuredGeometry.from_homogeneous(10.0, MirrorDeck())
        mesh1d = Mesher(geometry).partition(CellsByCount.uniform_width(4)).mesh
        return DiffusionMesh(mesh1d, {0: get_mixture("A", "2g")})

    def test_minimal_is_metadata_free(self):
        ms = DiffusionMethodSpace.minimal()
        assert ms.face is None and ms.trace is None

    def test_for_face_stores_the_mesh_trace(self):
        """#290 P7a: the stored trace IS the trace passed in (the mesh's
        cached scalar trace, the single source), and the method space
        holds no mesh."""
        mesh = self._mesh()
        ms = DiffusionMethodSpace.for_face(trace=mesh.scalar_trace, face="xmax")
        assert ms.face == "xmax"
        assert ms.trace is mesh.scalar_trace
        assert not hasattr(ms, "mesh")

    def test_for_face_validates_face_against_trace_inventory(self):
        """A method space for a nonexistent face is unrepresentable —
        the mirror of SN's for_face deriving (and thereby validating)
        per-face indices from its trace."""
        with pytest.raises(ValueError, match="ymax"):
            DiffusionMethodSpace.for_face(trace=self._mesh().scalar_trace, face="ymax")

    def test_frozen(self):
        ms = DiffusionMethodSpace.minimal()
        with pytest.raises(AttributeError):
            ms.face = "xmin"  # type: ignore[misc]
