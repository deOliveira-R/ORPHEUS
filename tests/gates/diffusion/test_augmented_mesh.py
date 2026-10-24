r"""DiffusionMesh construction laws — the diffusion phase space (#290 P7a).

Intrinsic-property gates for the method-mesh type (project standard:
every math-bearing type ships a test of its DEFINING laws). The
defining laws of :class:`~orpheus.diffusion.augmented_mesh.DiffusionMesh`
are its CONSTRUCTION invariants — what makes an instance a solvable
diffusion phase space rather than bare data:

* **realized boundary conditions at construction** — ``mesh.bc`` maps each
  boundary face to the realized albedo operator :math:`J^- =
  \mathcal{A} J^+` (vacuum = Marshak :math:`\mathcal{A}=0`; zero-flux =
  the Dirichlet idealization :math:`\mathcal{A}=-1`; a mirror deck
  :math:`\mathcal{A}=1`), joined to the faces by the binding of the
  declared ``boundary_conditions``;
* **operator coverage ≡ face coverage** — ``bc`` and the trace derive from
  the ONE ``face_labels`` inventory (the structural invariant that
  retired the P4 boundary operator's coverage validation);
* **the composite carrier identity** — ``full_field_space`` is built
  over THIS mesh's cached ``scalar_trace`` (operator/field identity
  guards ride on it);
* **promotion** — ``from_material_mesh`` re-derives the data block
  bit-identically; the only extra input is the keyword
  ``boundary_conditions`` (the tags live on the axes' roles, the
  conditions they name are the method's input);
* **admission refusals** — multi-D, a source (#656), a translation deck,
  an unbound tag (the mesh-less refusal retired at the CS4c coda with the
  fabricated carrier it refused): each fires AT CONSTRUCTION, so "a
  diffusion phase space with unresolved/unrealizable BCs" is
  unrepresentable (the refusals moved here from the P5 solver at P7a).
  The two refusals the string registry made (``white`` absent from the
  admitted tags, an ``albedo`` tag with no parameter) are not
  refusals of the typed vocabulary: ``white(a)`` is admitted at tier 2 and
  an albedo without its number is not constructible.

Foundation tier — software invariants on the phase-space type; the
downstream physics (stencil, k, MMS) is gated in ``test_operators.py``
/ ``test_solver.py`` / ``test_mms.py``.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.diffusion import DiffusionMesh
from orpheus.diffusion import ZeroFluxCondition
from orpheus.geometry import StructuredGeometry
from orpheus.geometry.deck import MirrorDeck, TranslationDeck
from orpheus.mesh import CellsByCount, Mesh1D, Mesh2D, Mesher
from orpheus.transport.boundary_conditions import BoundaryRefused, RefusalAxis, SubmarkovViolationError, UnboundTagError, WallCondition
from orpheus.transport.mesh.material_mesh import MaterialMesh
from tests.gates._exterior import declare_sides

pytestmark = [pytest.mark.foundation]

_MATS = {0: get_mixture("A", "2g")}


def _mesh1d(left: Any = MirrorDeck(), right: Any = MirrorDeck()) -> tuple[Mesh1D, dict]:
    """The slab and the ``{tag: condition}`` of its tags: a side is a deck or a condition (``declare_sides`` mints the tags, rule R1)."""
    (role_left, role_right), conditions = declare_sides(left, right)
    geometry = StructuredGeometry.slab((0.0, 10.0), (0,), left=role_left, right=role_right)
    return Mesher(geometry).partition(CellsByCount.uniform_width(4)).mesh, conditions


def _diffusion_mesh(left: Any = MirrorDeck(), right: Any = MirrorDeck()) -> DiffusionMesh:
    mesh, conditions = _mesh1d(left, right)
    return DiffusionMesh(mesh, _MATS, boundary_conditions=conditions)


_PROBE = np.array([1.5, 2.5])   # a J⁺ face slot, (ng,)


class TestRealizedBoundaryConditions:
    def test_bc_realized_per_face_with_ruling3_semantics(self):
        """The condition → 𝒜 table THROUGH the mesh path: vacuum J⁻=0,
        zero-flux J⁻=−J⁺ — realized operators, not stored tags."""
        dm = _diffusion_mesh(WallCondition.vacuum(), ZeroFluxCondition())
        np.testing.assert_array_equal(
            np.asarray(dm.bc["xmin"].apply(_PROBE)), 0.0,
        )
        np.testing.assert_array_equal(
            np.asarray(dm.bc["xmax"].apply(_PROBE)), -_PROBE,
        )

    def test_mirror_and_albedo_actions(self):
        dm = _diffusion_mesh(MirrorDeck(), WallCondition.white(0.3))
        np.testing.assert_array_equal(
            np.asarray(dm.bc["xmin"].apply(_PROBE)), _PROBE,
        )
        np.testing.assert_allclose(
            np.asarray(dm.bc["xmax"].apply(_PROBE)), 0.3 * _PROBE,
            rtol=1e-15,
        )

    def test_operator_coverage_equals_face_coverage(self):
        """bc and the trace derive from the ONE face_labels inventory —
        the structural invariant that made the P4 boundary operator's
        coverage validation unrepresentable-state ceremony."""
        dm = _diffusion_mesh(WallCondition.vacuum(), MirrorDeck())
        assert set(dm.bc) == set(dm.scalar_trace.face_names)

    def test_sphere_pole_is_not_a_bc_face(self):
        """Curvilinear: the r=0 pole is a regularity condition, not a
        face — bc and trace both carry only the outer face."""
        sphere = Mesher(
            StructuredGeometry.sphere((0.0, 5.0), (0,), outer=0),
        ).partition(CellsByCount.uniform_width(3)).mesh
        # A curved face has no mirror deck: the perfect mirror is the specular wall tag (rule R2).
        dm = DiffusionMesh(sphere, _MATS, boundary_conditions={0: WallCondition.specular(1.0)})
        assert dm.scalar_trace.face_names == ("xmax",)
        assert set(dm.bc) == {"xmax"}


class TestSpaces:
    def test_composite_carrier_rides_the_cached_trace(self):
        dm = _diffusion_mesh()
        ffs = dm.full_field_space
        assert ffs.trace_space is dm.scalar_trace
        assert ffs is dm.full_field_space          # cached — one identity
        # CS4b S4: the interior IS the carrier's cached axis-built mint
        # (the Q2 three-family unification's diffusion leg) — and its
        # metric ACTION is cell volumes broadcast over the group axis
        # (asserted on the action, not the storage: axis-built spaces
        # carry per-axis weights, not one dense array).
        assert ffs.interior_space is dm.bulk_space
        x = np.arange(1.0, 1.0 + float(np.prod(ffs.interior_space.shape))).reshape(
            ffs.interior_space.shape
        )
        np.testing.assert_array_equal(
            ffs.interior_space.apply_metric(x),
            x * np.asarray(dm.volumes)[None, :],
        )


class TestPromotion:
    def test_from_material_mesh_is_a_bit_identical_data_join(self):
        """The promotion re-derives the data block from the carrier's
        own axes/mesh/mat_map/materials — values bit-identical; the
        conditions the axes' tags stand for come in by keyword."""
        mesh, conditions = _mesh1d(WallCondition.vacuum(), MirrorDeck())
        mm = MaterialMesh(mesh, _MATS)
        dm = DiffusionMesh.from_material_mesh(mm, boundary_conditions=conditions)
        assert isinstance(dm, MaterialMesh)
        assert dm.materials is mm.materials
        assert dm.mesh is mm.mesh
        assert dm.axes == mm.axes
        np.testing.assert_array_equal(dm.mat_map, mm.mat_map)
        np.testing.assert_array_equal(dm.volumes, mm.volumes)
        assert dm.ng == mm.ng
        # The conditions realized: vacuum on xmin (J⁻ = 0).
        np.testing.assert_array_equal(
            np.asarray(dm.bc["xmin"].apply(_PROBE)), 0.0,
        )


class TestAdmissionRefusals:
    def test_multi_d_mesh_is_refused_at_construction(self):
        mesh2d = Mesh2D(
            edges_x=np.array([0.0, 1.0, 2.0]),
            edges_y=np.array([0.0, 1.0]),
            mat_map=np.zeros((2, 1), dtype=int),
            face_roles={"xmin": MirrorDeck(), "xmax": MirrorDeck(), "ymin": MirrorDeck(), "ymax": MirrorDeck()},
        )
        with pytest.raises(ValueError, match="1-D"):
            DiffusionMesh(mesh2d, _MATS)

    def test_multi_d_promotion_is_refused(self):
        mesh2d = Mesh2D(
            edges_x=np.array([0.0, 1.0, 2.0]),
            edges_y=np.array([0.0, 1.0]),
            mat_map=np.zeros((2, 1), dtype=int),
            face_roles={"xmin": MirrorDeck(), "xmax": MirrorDeck(), "ymin": MirrorDeck(), "ymax": MirrorDeck()},
        )
        with pytest.raises(ValueError, match="1-D"):
            DiffusionMesh.from_material_mesh(MaterialMesh(mesh2d, _MATS))

    def test_white_is_admitted_at_tier_two_as_the_same_operator_as_specular(self):
        """Re-posed from the retired "unsupported tag" refusal: the string registry omitted ``white`` ("at P1 it coincides with
        reflective"); the typed table ADMITS it, stating the coincidence as the tier (``Realized(False)``): ``white(0.5)`` and
        ``specular(0.5)`` realize one operator, 𝒜 = 0.5."""
        white = _diffusion_mesh(WallCondition.white(0.5), WallCondition.vacuum())
        specular = _diffusion_mesh(WallCondition.specular(0.5), WallCondition.vacuum())
        np.testing.assert_array_equal(np.asarray(white.bc["xmin"].apply(_PROBE)), 0.5 * _PROBE)
        np.testing.assert_array_equal(np.asarray(white.bc["xmin"].apply(_PROBE)), np.asarray(specular.bc["xmin"].apply(_PROBE)))

    def test_an_albedo_that_is_not_a_physical_number_is_not_constructible(self):
        """Re-posed from the retired "albedo without parameter" refusal: the albedo is the argument of the facade and is parsed at
        construction, so a missing or non-numeric one cannot be spelled (``None``: ``TypeError``; NaN: ``ValueError``), and neither
        can one above one (``SubmarkovViolationError``, the cone parse: #265's albedo item)."""
        with pytest.raises(TypeError):
            WallCondition.white(None)  # type: ignore[arg-type]
        with pytest.raises(ValueError, match="NaN"):
            WallCondition.white(float("nan"))
        with pytest.raises(SubmarkovViolationError):
            WallCondition.white(5.0)

    def test_a_source_is_refused_as_plumbing_at_construction(self):
        """A boundary source has no channel in diffusion yet: ``Refused(PLUMBING, #656)`` at construction, naming the face."""
        with pytest.raises(BoundaryRefused) as info:
            _diffusion_mesh(WallCondition.uniform(1.0), WallCondition.vacuum())
        assert info.value.admission.axis is RefusalAxis.PLUMBING and info.value.admission.issue == 656 and info.value.face == "xmin"

    def test_a_translation_deck_is_refused_as_topological_at_construction(self):
        """The periodic slab couples two faces, which a per-face scalar cannot hold: ``Refused(TOPOLOGY)`` at construction."""
        with pytest.raises(BoundaryRefused) as info:
            _diffusion_mesh(TranslationDeck(), TranslationDeck())
        assert info.value.admission.axis is RefusalAxis.TOPOLOGY

    def test_an_unbound_tag_is_refused_at_construction(self):
        """No default wall: a face tagged with no condition declared is ``UnboundTagError``."""
        mesh, _ = _mesh1d(WallCondition.vacuum(), MirrorDeck())
        with pytest.raises(UnboundTagError):
            DiffusionMesh(mesh, _MATS)
