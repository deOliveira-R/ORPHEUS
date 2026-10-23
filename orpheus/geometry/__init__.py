"""The geometry layer: shapes, coordinate systems and face roles, with no discretisation.

* :class:`StructuredGeometry`: the 1-D layered geometry (a coordinate
  system, breakpoints, one material id per interval and one face role per
  boundary point), a pure shape and boundary description with no cell
  counts. Reference solvers consume it directly. A role is a deck or an
  ``int`` tag; the geometry names no response.
* :class:`CoordSystem` and the volume and area formulas of each coordinate
  system, and the :class:`AngularChart` each declares for a direction.
* :mod:`~orpheus.geometry.transformation`: rigid motions and permutations.
* :mod:`~orpheus.geometry.symmetry`: the closed subgroups of :math:`O(3)`,
  realized; :mod:`~orpheus.geometry.roots_of_unity`, the exact realization of
  :math:`C_n` they generate.
* :mod:`~orpheus.geometry.manifold`: the point sets the groups act on (the
  :class:`~orpheus.geometry.manifold.Manifold` closed sum, its G-spaces and
  orbit spaces), which numerics' measures are supported on.
* :mod:`~orpheus.geometry.deck`: the decks (how a domain continues across a face), a closed sum of
  :class:`MirrorDeck` and :class:`TranslationDeck`; :mod:`~orpheus.geometry.colouring`: the stabiliser of a
  layered geometry's fundamental domain.
* :mod:`~orpheus.geometry.faces`: the axis letters and the boundary-face
  naming bijection.

The last three moved down from :mod:`orpheus.numerics` at layer DAG WP-B
(2026-10-22): numerics imports the points and the group action, never the
reverse.

The mesh is an overlay on the geometry and lives in its own package,
:mod:`orpheus.mesh`, which imports this one and never the reverse:
a :class:`~orpheus.mesh.mesher.Mesher` discretises a
:class:`StructuredGeometry`.
"""

from .extent import GeometryExtent
from .colouring import ColouringSymmetry, FaceIsometry, Symmetries, UncomputedSymmetries, colouring_symmetry
from .deck import (
    DeckError,
    DeckTransformation,
    FacePairing,
    FaceRole,
    IsometryClass,
    MirrorDeck,
    TranslationDeck,
    admitted_decks,
    check_pairings_close,
    mirror_motion,
    pairing_of,
    parse_face_role,
    translation_motion,
)
from .coord import AngularChart, CoordSystem, MeasureCoordinate, compute_areas_1d, compute_volumes_1d, compute_volumes_2d
from .structured_geometry import StructuredGeometry
from .transformation import (
    NotAFinitePointGroupError,
    Permutation,
    RigidMotion,
    close_group,
)

__all__ = [
    "GeometryExtent",
    "AngularChart",
    "ColouringSymmetry",
    "CoordSystem",
    "DeckError",
    "DeckTransformation",
    "FaceIsometry",
    "FacePairing",
    "FaceRole",
    "IsometryClass",
    "MirrorDeck",
    "Symmetries",
    "TranslationDeck",
    "UncomputedSymmetries",
    "MeasureCoordinate",
    "NotAFinitePointGroupError",
    "Permutation",
    "RigidMotion",
    "StructuredGeometry",
    "admitted_decks",
    "check_pairings_close",
    "close_group",
    "colouring_symmetry",
    "compute_areas_1d",
    "compute_volumes_1d",
    "compute_volumes_2d",
    "mirror_motion",
    "pairing_of",
    "parse_face_role",
    "translation_motion",
]
