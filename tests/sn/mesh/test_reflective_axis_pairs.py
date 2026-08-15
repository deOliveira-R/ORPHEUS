r"""``SNMesh.reflective_axis_pairs`` — the geometry half of the gauge predicate.

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

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.solver import _as_sn_mesh
from orpheus.transport.mesh.axis import AxisMesh

_R = BC("reflective")
_V = BC("vacuum")
_EDGES = np.linspace(0.0, 1.0, 4)


def _sn_mesh_2d(**boundary_conditions) -> SNMesh:
    mesh = Mesh2D(
        edges_x=_EDGES, edges_y=_EDGES,
        mat_map=np.zeros((3, 3), dtype=int),
        **boundary_conditions,
    )
    return SNMesh(
        mesh, Quadrature.level_symmetric(4), {0: get_mixture("B", "2g")},
    )


@pytest.mark.foundation
def test_an_all_reflective_box_counts_every_axis():
    assert _sn_mesh_2d(
        bc_xmin=_R, bc_xmax=_R, bc_ymin=_R, bc_ymax=_R,
    ).reflective_axis_pairs == 2


@pytest.mark.foundation
def test_one_vacuum_face_breaks_its_axis():
    """`[M]` #344: this is the configuration where ``dim ker A`` drops 12 → 0."""
    assert _sn_mesh_2d(
        bc_xmin=_R, bc_xmax=_R, bc_ymin=_R, bc_ymax=_V,
    ).reflective_axis_pairs == 1


@pytest.mark.foundation
def test_a_MIXED_axis_contributes_nothing():
    """``xmin=R, xmax=V`` is worth ZERO, not a half.

    A mode travelling out through ``xmax`` escapes; there is no return
    leg, so the x axis cannot host a closed loop no matter how
    reflective its other face is.
    """
    mesh = _sn_mesh_2d(bc_xmin=_R, bc_xmax=_V, bc_ymin=_R, bc_ymax=_R)
    assert mesh.reflective_axis_pairs == 1

    reflective_faces = sum(
        1 for face in mesh.bc
        if type(mesh.bc[face].law).__name__ == "ReflectiveBoundary"
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
    mesh = _sn_mesh_2d(bc_xmin=_R, bc_xmax=_V, bc_ymin=_R, bc_ymax=_V)

    reflective_faces = sum(
        1 for face in mesh.bc
        if type(mesh.bc[face].law).__name__ == "ReflectiveBoundary"
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
        bc_xmin=_V, bc_xmax=_V, bc_ymin=_V, bc_ymax=_V,
    ).reflective_axis_pairs == 0


@pytest.mark.foundation
def test_a_BARE_mesh_is_all_reflective_and_the_predicate_sees_it():
    """⚠ The silent default, and why the exposure is broad.

    ``resolve_boundary_conditions`` fills unset faces with
    ``BC("reflective")``, and neither ``solve_sn`` nor
    ``solve_sn_adjoint`` even takes a ``boundary_condition`` argument —
    so the canonical ``k_inf`` lattice, written with zero BC arguments,
    is in the singular class.  The predicate must read the REALIZED law,
    not the tag a caller passed.
    """
    mesh = Mesh2D(
        edges_x=_EDGES, edges_y=_EDGES, mat_map=np.zeros((3, 3), dtype=int),
    )
    sn_mesh = SNMesh(
        mesh, Quadrature.level_symmetric(4), {0: get_mixture("B", "2g")},
    )
    assert sn_mesh.reflective_axis_pairs == sn_mesh.ndim == 2


@pytest.mark.foundation
def test_it_counts_at_three_dimensions_too():
    """Derived from the face inventory, so a third axis needs no edit."""
    quadrature = Quadrature.level_symmetric(4)
    materials = {0: get_mixture("B", "2g")}

    def _axis(low, high):
        return AxisMesh(edges=_EDGES, bc_low=low, bc_high=high)

    all_reflective = (_axis(_R, _R),) * 3
    sn_mesh = _as_sn_mesh(all_reflective, quadrature, materials)
    assert sn_mesh.ndim == 3
    assert sn_mesh.reflective_axis_pairs == 3

    one_vacuum_pair = (_axis(_V, _V),) + all_reflective[1:]
    assert _as_sn_mesh(
        one_vacuum_pair, quadrature, materials,
    ).reflective_axis_pairs == 2
