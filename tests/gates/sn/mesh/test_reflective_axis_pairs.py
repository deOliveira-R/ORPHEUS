r"""``SNProblem.reflective_axis_pairs`` — the geometry half of the gauge predicate.

An undamped face mode (the closure half —
:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.face_transmission_spectrum`)
only closes into a null vector of :math:`A = L + C - S - B` if it can
return to itself, which needs a **closed** reflective loop.  `[M]` #344
measured the requirement as **two** such axes: at ``d = 2`` a single
vacuum face collapses ``dim ker A`` from 12 to 0.

⭐ Counts PAIRS, not faces.  :func:`test_a_MIXED_axis_contributes_nothing`
is the load-bearing row: counting faces would call ``xmin=R, xmax=V`` a
half-reflective axis and let a mode "close" through a leg that actually
escapes.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.data.macro_xs.synthetic import get_mixture
from orpheus.geometry.deck import MirrorDeck
from orpheus.mesh import Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.problem import SNProblem
from orpheus.sn.solver import _as_problem
from orpheus.mesh import AxisMesh
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for

#: The roles this file is written in: a mirror deck is the reflective face, the tag 0 (declared the vacuum wall) the vacuum one.
_R = MirrorDeck()
_V = 0
_CONDITIONS = {0: WallCondition.vacuum()}
_EDGES = np.linspace(0.0, 1.0, 4)


def _sn_mesh_2d(conditions=None, **face_roles) -> SNProblem:
    mesh = Mesh2D(
        edges_x=_EDGES, edges_y=_EDGES,
        mat_map=np.zeros((3, 3), dtype=int),
        face_roles=face_roles,
    )
    return SNProblem(
        mesh, Quadrature.level_symmetric(4), {0: get_mixture("B", "2g")},
        boundary_conditions=conditions_for(mesh, _CONDITIONS if conditions is None else conditions),
    )


@pytest.mark.foundation
def test_an_all_reflective_box_counts_every_axis():
    assert _sn_mesh_2d(
        xmin=_R, xmax=_R, ymin=_R, ymax=_R,
    ).reflective_axis_pairs == 2


@pytest.mark.foundation
def test_one_vacuum_face_breaks_its_axis():
    """`[M]` #344: this is the configuration where ``dim ker A`` drops 12 → 0."""
    assert _sn_mesh_2d(
        xmin=_R, xmax=_R, ymin=_R, ymax=_V,
    ).reflective_axis_pairs == 1


@pytest.mark.foundation
def test_a_MIXED_axis_contributes_nothing():
    """``xmin=R, xmax=V`` is worth ZERO, not a half.

    A mode travelling out through ``xmax`` escapes; there is no return
    leg, so the x axis cannot host a closed loop no matter how
    reflective its other face is.
    """
    mesh = _sn_mesh_2d(xmin=_R, xmax=_V, ymin=_R, ymax=_R)
    assert mesh.reflective_axis_pairs == 1

    reflective_faces = sum(
        1 for face in mesh.bc
        if isinstance(mesh.binding.decks.get(face), MirrorDeck)
    )
    assert reflective_faces == 3


@pytest.mark.foundation
def test_TWO_mixed_axes_are_worth_zero_not_one():
    """⭐ The row that actually separates counting PAIRS from counting FACES.

    ⚠ Mutation-found.  The obvious mixed-axis fixture above does NOT
    discriminate: with ``xmin=R, xmax=V, ymin=R, ymax=R`` a
    faces-counting implementation (``reflective_faces // 2 = 3 // 2``)
    returns 1 — the right answer for the wrong reason.  `[M]` a
    ``faces // 2`` mutation left all six of the original rows GREEN.

    Both axes mixed is where the two rules part company: 2 reflective
    faces, ``2 // 2 = 1`` by the face rule, but **0** closed loops —
    every path out of this box escapes on its first bounce.
    """
    mesh = _sn_mesh_2d(xmin=_R, xmax=_V, ymin=_R, ymax=_V)

    reflective_faces = sum(
        1 for face in mesh.bc
        if isinstance(mesh.binding.decks.get(face), MirrorDeck)
    )
    assert reflective_faces == 2, "fixture drifted"
    assert reflective_faces // 2 == 1, (
        "fixture no longer discriminates: a faces-counting implementation "
        "must disagree with the pair count here, or this row is inert"
    )
    assert mesh.reflective_axis_pairs == 0


@pytest.mark.foundation
def test_an_all_vacuum_box_has_none():
    assert _sn_mesh_2d(
        xmin=_V, xmax=_V, ymin=_V, ymax=_V,
    ).reflective_axis_pairs == 0


@pytest.mark.foundation
def test_a_unit_specular_WALL_closes_its_axis_exactly_as_a_mirror_deck_does():
    """The predicate reads the REALIZED operator of each face, not the declaration it came from (S2b amend; first red of the review).

    A mirror deck and a wall tagged with the unit specular condition realize the same ``ScaledOrdinatePermutation`` at ``alpha = 1``,
    hence the same loss operator and the same kernel. ``[M]`` at the reviewed commit (LS4, 3 x 3 box) the deck box returned 2 closed
    axes and ``dim ker A = 12`` and the wall box returned 0 and no gauge, with ``array_equal(A_deck, A_wall)`` True: the gauge freedom
    was decided by how the wall was declared, so the trace of a wall box was returned un-gauged with a tangential current of
    7.4e-2 on ``ymin`` and ``ymax`` where the exact answer is 0. A PARTIAL specular wall returns a fraction and closes nothing; a
    unit specular wall facing a vacuum face is a mixed axis and is worth zero, as a deck would be.
    """
    from orpheus.sn.operators.loss_kernel_gauge import gauge_freedom

    walls = _sn_mesh_2d({0: WallCondition.specular(1.0)}, xmin=0, xmax=0, ymin=0, ymax=0)
    decks = _sn_mesh_2d(xmin=_R, xmax=_R, ymin=_R, ymax=_R)
    assert walls.reflective_axes == decks.reflective_axes == (0, 1)
    assert gauge_freedom(walls).present and gauge_freedom(decks).present
    assert walls.loss_kernel_gauge.dimension == decks.loss_kernel_gauge.dimension == 12

    partial = _sn_mesh_2d({0: WallCondition.specular(0.9)}, xmin=0, xmax=0, ymin=0, ymax=0)
    assert partial.reflective_axis_pairs == 0

    mixed = _sn_mesh_2d(
        {0: WallCondition.specular(1.0), 1: WallCondition.vacuum()}, xmin=0, xmax=1, ymin=0, ymax=0,
    )
    assert mixed.reflective_axes == (1,)


@pytest.mark.foundation
def test_it_counts_at_three_dimensions_too():
    """Derived from the face inventory, so a third axis needs no edit."""
    quadrature = Quadrature.level_symmetric(4)
    materials = {0: get_mixture("B", "2g")}

    def _axis(low, high):
        return AxisMesh(edges=_EDGES, role_low=low, role_high=high)

    all_reflective = (_axis(_R, _R),) * 3
    problem = _as_problem(all_reflective, quadrature, materials)
    assert problem.ndim == 3
    assert problem.reflective_axis_pairs == 3

    one_vacuum_pair = (_axis(_V, _V),) + all_reflective[1:]
    assert _as_problem(
        one_vacuum_pair, quadrature, materials, boundary_conditions=conditions_for(one_vacuum_pair, _CONDITIONS),
    ).reflective_axis_pairs == 2
