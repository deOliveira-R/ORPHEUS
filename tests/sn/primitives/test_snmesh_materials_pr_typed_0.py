r"""SNMesh materials parameter — Issue #197 PR-TYPED-0 acceptance tests.

Pins the §F mechanism criteria for the PR:

* Criterion 2: ``SNMesh(mesh, quad)`` without ``materials`` raises TypeError.
* Criterion 3: ``SNMesh.ng`` property exists + returns the uniform ng.
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
from orpheus.geometry import BC, Mesh1D
from orpheus.sn.mesh.augmented_mesh import InconsistentMaterialsError, SNMesh
from orpheus.numerics.quadrature import Quadrature

pytestmark = pytest.mark.foundation


def _mix(ng: int) -> Mixture:
    z = np.zeros(ng)
    z_mat = csr_matrix(np.zeros((ng, ng)))
    return Mixture(
        SigC=z.copy(), SigL=z.copy(), SigF=z.copy(),
        SigP=z.copy(), SigT=np.ones(ng),
        SigS=[z_mat], Sig2=z_mat, chi=z.copy(),
    )


def _slab_mesh(mat_ids=None) -> Mesh1D:
    if mat_ids is None:
        mat_ids = np.zeros(4, dtype=int)
    return Mesh1D(
        edges=np.linspace(0.0, 1.0, len(mat_ids) + 1),
        mat_ids=np.asarray(mat_ids, dtype=int),
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )


def test_materials_required_positional_arg() -> None:
    """Criterion 2: SNMesh(mesh, quad) without materials raises TypeError."""
    mesh = _slab_mesh()
    quad = Quadrature.gauss_legendre(4)
    with pytest.raises(TypeError, match="materials"):
        SNMesh(mesh, quad)


def test_ng_property_returns_uniform_ng() -> None:
    """Criterion 3: SNMesh.ng returns the materials' uniform ng."""
    mesh = _slab_mesh()
    quad = Quadrature.gauss_legendre(4)
    sn_mesh = SNMesh(mesh, quad, {0: _mix(ng=2)})
    assert sn_mesh.ng == 2

    sn_mesh4 = SNMesh(mesh, quad, {0: _mix(ng=4)})
    assert sn_mesh4.ng == 4


def test_inconsistent_ng_raises_inconsistent_materials_error() -> None:
    """Criterion 4: mismatched ng across materials raises InconsistentMaterialsError."""
    mesh = _slab_mesh(mat_ids=[0, 0, 1, 1])
    quad = Quadrature.gauss_legendre(4)
    materials = {0: _mix(ng=2), 1: _mix(ng=4)}
    with pytest.raises(InconsistentMaterialsError, match="uniform ng"):
        SNMesh(mesh, quad, materials)


def test_missing_material_id_raises_value_error() -> None:
    """Criterion 5: ``mat_map`` referencing an id missing from materials raises."""
    mesh = _slab_mesh(mat_ids=[0, 1, 1, 0])
    quad = Quadrature.gauss_legendre(4)
    # Only id 0 is present; id 1 used in mat_ids triggers the validation.
    with pytest.raises(ValueError, match=r"material ids \[1\]"):
        SNMesh(mesh, quad, {0: _mix(ng=2)})


def test_empty_materials_raises_value_error() -> None:
    mesh = _slab_mesh()
    quad = Quadrature.gauss_legendre(4)
    with pytest.raises(ValueError, match="non-empty materials"):
        SNMesh(mesh, quad, {})


def test_materials_attribute_is_dict_passed() -> None:
    """The ``.materials`` attribute carries the dict from the constructor."""
    mesh = _slab_mesh()
    quad = Quadrature.gauss_legendre(4)
    materials = {0: _mix(ng=2)}
    sn_mesh = SNMesh(mesh, quad, materials)
    # Un-weld arc (R20/R21): parsed into the stage-1 declaration at the
    # boundary — entries identical, the mapping itself no longer aliased.
    from orpheus.data.materials import Materials
    assert isinstance(sn_mesh.materials, Materials)
    assert all(sn_mesh.materials[k] is materials[k] for k in materials)


def test_inconsistent_materials_error_is_value_error() -> None:
    """InconsistentMaterialsError is a ValueError subclass for ``except`` ergonomics."""
    assert issubclass(InconsistentMaterialsError, ValueError)
