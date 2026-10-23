r"""SNProblem materials parameter — Issue #197 PR-TYPED-0 acceptance tests.

Pins the §F mechanism criteria for the PR:

* Criterion 2: ``SNProblem(mesh, quad)`` without ``materials`` raises TypeError.
* Criterion 3: ``SNProblem.ng`` property exists + returns the uniform ng.
* Criterion 4: heterogeneous-``ng`` materials raise InconsistentMaterialsError.
* Criterion 5: ``mat_map`` referencing a missing material id raises ValueError.
* Plus: empty materials dict raises ValueError; ``.materials`` accessor
  matches what was passed.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByMaxWidth, EqualWidth, Mesh1D, Mesher
from orpheus.data.materials import InconsistentMaterialsError
from orpheus.sn.problem import SNProblem
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for

#: Exterior conditions of the tags this file's meshes carry (wave 2: a face's role is a deck or an int tag; the condition a tag
#: stands for is declared at the solve, keyword ``boundary_conditions``).
_CONDITIONS = {0: WallCondition.vacuum()}


pytestmark = pytest.mark.foundation


def _mix(ng: int) -> Mixture:
    z = np.zeros(ng)
    z_mat = csr_matrix(np.zeros((ng, ng)))
    return Mixture(
        SigC=z.copy(), SigL=z.copy(), SigF=z.copy(),
        SigP=z.copy(), SigT=np.ones(ng),
        SigS=[z_mat], Sig2=[z_mat], chi=z.copy(),
    )


def _slab_mesh(
    breakpoints: tuple[float, ...] = (0.0, 1.0), mat_ids: tuple[int, ...] = (0,),
) -> Mesh1D:
    """A vacuum slab over ``[0, 1]`` in cells of width 1/4, one material per interval."""
    return Mesher(StructuredGeometry.slab(
        breakpoints, mat_ids, left=0, right=0,
    )).partition(CellsByMaxWidth(0.25, EqualWidth())).mesh


def test_materials_required_positional_arg() -> None:
    """Criterion 2: SNProblem(mesh, quad) without materials raises TypeError."""
    mesh = _slab_mesh()
    quad = Quadrature.gauss_legendre(4)
    with pytest.raises(TypeError, match="materials"):
        SNProblem(mesh, quad, boundary_conditions=conditions_for(mesh, _CONDITIONS))


def test_ng_property_returns_uniform_ng() -> None:
    """Criterion 3: SNProblem.ng returns the materials' uniform ng."""
    mesh = _slab_mesh()
    quad = Quadrature.gauss_legendre(4)
    problem = SNProblem(mesh, quad, {0: _mix(ng=2)}, boundary_conditions=conditions_for(mesh, _CONDITIONS))
    assert problem.ng == 2

    sn_mesh4 = SNProblem(mesh, quad, {0: _mix(ng=4)}, boundary_conditions=conditions_for(mesh, _CONDITIONS))
    assert sn_mesh4.ng == 4


def test_inconsistent_ng_raises_inconsistent_materials_error() -> None:
    """Criterion 4: mismatched ng across materials raises InconsistentMaterialsError."""
    mesh = _slab_mesh((0.0, 0.5, 1.0), (0, 1))
    quad = Quadrature.gauss_legendre(4)
    materials = {0: _mix(ng=2), 1: _mix(ng=4)}
    with pytest.raises(InconsistentMaterialsError, match="uniform ng"):
        SNProblem(mesh, quad, materials, boundary_conditions=conditions_for(mesh, _CONDITIONS))


def test_missing_material_id_raises_value_error() -> None:
    """Criterion 5: ``mat_map`` referencing an id missing from materials raises."""
    mesh = _slab_mesh((0.0, 0.25, 0.75, 1.0), (0, 1, 0))
    quad = Quadrature.gauss_legendre(4)
    # Only id 0 is present; id 1 used in mat_ids triggers the validation.
    with pytest.raises(ValueError, match=r"material ids \[1\]"):
        SNProblem(mesh, quad, {0: _mix(ng=2)}, boundary_conditions=conditions_for(mesh, _CONDITIONS))


def test_empty_materials_raises_value_error() -> None:
    mesh = _slab_mesh()
    quad = Quadrature.gauss_legendre(4)
    with pytest.raises(ValueError, match="non-empty materials"):
        SNProblem(mesh, quad, {}, boundary_conditions=conditions_for(mesh, _CONDITIONS))


def test_materials_attribute_is_dict_passed() -> None:
    """The ``.materials`` attribute carries the dict from the constructor."""
    mesh = _slab_mesh()
    quad = Quadrature.gauss_legendre(4)
    materials = {0: _mix(ng=2)}
    problem = SNProblem(mesh, quad, materials, boundary_conditions=conditions_for(mesh, _CONDITIONS))
    # Un-weld arc (R20/R21): parsed into the stage-1 declaration at the
    # boundary — entries identical, the mapping itself no longer aliased.
    from orpheus.data.materials import Materials
    assert isinstance(problem.materials, Materials)
    assert all(problem.materials[k] is materials[k] for k in materials)


def test_inconsistent_materials_error_is_value_error() -> None:
    """InconsistentMaterialsError is a ValueError subclass for ``except`` ergonomics."""
    assert issubclass(InconsistentMaterialsError, ValueError)
