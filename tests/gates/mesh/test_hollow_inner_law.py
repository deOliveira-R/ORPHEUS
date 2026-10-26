"""A face role is a deck or a tag: the mesh refuses ``None`` as the inner role of a hollow curvilinear body.

Until P1 step 3b an undeclared (``None``) inner law of a hollow cylinder or sphere was legal, and the RECORD here showed that it
gave the same SN flux, bit for bit, as a declared reflective one. ``None`` is no longer a role: the bare ``Mesh1D`` constructor
refuses it on every face (``TypeError``), so the only spelling of that cavity is the inner tag bound to
``WallCondition.specular(1.0)``. The row executes only the mesh.

This file held SN's and diffusion's refusals of a declared inner condition (#511, the P1 specification's S3.9) until the layer-DAG
re-filing (Phase 1b, 2026-10-26). Each refusal is one method's claim, so those rows moved, bodies unchanged:
``tests/gates/sn/architecture/test_sn_hollow_inner_law.py`` (SN's refusal and the void-cavity verification of the admitted unit
specular wall) and ``tests/gates/diffusion/test_diffusion_hollow_inner_law.py``, with their shared body in
``tests/gates/transport/_hollow_inner_helpers.py``. ``_direct_hollow_mesh`` is also defined in the SN file, which builds the same
bare mesh with a legal role (a mesh gate may not import a transport helper).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import CoordSystem
from orpheus.mesh import Mesh1D

pytestmark = pytest.mark.foundation

_CURVILINEAR = (CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL)
#: The outer face of the hollow body is the tag 0 (rule R6).
_OUTER = 0


def _direct_hollow_mesh(coord: CoordSystem, inner_role) -> Mesh1D:
    """The same body built by the bare constructor, with no geometry."""
    edges = np.linspace(0.5, 2.0, 9)
    return Mesh1D(
        coord=coord, edges=edges, volumes=coord.measure(edges),
        region_ids=np.zeros(8, int), region_materials=(0,), face_roles={"xmin": inner_role, "xmax": _OUTER},
    )


@pytest.mark.parametrize("coord", _CURVILINEAR, ids=lambda c: c.name.lower())
def test_an_undeclared_inner_law_is_the_reflective_one(coord):
    """RE-POSED at P1 step 3b. Until then this was a RECORD that an
    undeclared (``None``) inner law and a declared reflective one gave the
    same flux bit for bit. ``None`` is no longer a role: the mesh refuses it,
    so the only way to spell that cavity is the inner tag bound to
    ``WallCondition.specular(1.0)``, whose answer
    ``tests/gates/sn/architecture/test_sn_hollow_inner_law.py::test_a_reflective_inner_law_is_a_void_cavity``
    verifies."""
    with pytest.raises(TypeError, match="a face role is a MirrorDeck, a TranslationDeck or an int tag, got None"):
        _direct_hollow_mesh(coord, None)
