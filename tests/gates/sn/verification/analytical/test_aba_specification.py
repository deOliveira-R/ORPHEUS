r"""The A|B|A cross-check problem is posed once, and the SN fixtures pose that problem (#405 P2 step 7b.2, gate R7b2.1).

:mod:`._aba_reference` defines the A|B|A specification. The cross-check's SN
fixtures spell the same body in three places (the snapshot generator's
``_sphere_3region``/``_cylinder_3region``, the standoff and unified files'
cylinder meshes) and its materials in three ways (``get_mixture``, the
retiring helper's ``aba_xs_2g`` arrays, the ``_make_2g_mixture`` rebuilds).
Nothing asserted that they pose one problem (X4); these rows do: the
geometry's breakpoints, region materials and outer law, and the isotropic
transport data (σ_t, the P0 transfer matrix, νΣ_f, χ, bit for bit).

First red on ``65b0de93``: the helper ``_aba_reference`` did not exist (this
change writes it). Declared finding the rows pin: ``get_mixture("A"|"B",
"2g")`` carries a P1 moment the cross-check never solves (both sides are
isotropic), so the specification's materials are the P0 rebuilds.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem
from tests.gates.sn.verification.analytical._aba_reference import (
    ABA_MATERIAL_IDS,
    ABA_RADII,
    aba_materials,
    aba_specification,
)

pytestmark = pytest.mark.foundation

_COORDS = (CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL)


def _transport_data(mixture):
    return (np.asarray(mixture.SigT), mixture.SigS[0].toarray(), np.asarray(mixture.SigP), np.asarray(mixture.chi))


@pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
def test_r7b2_1_the_specification(coord) -> None:
    spec = aba_specification(coord)
    assert spec.geometry.breakpoints == (0.0, *ABA_RADII)
    assert tuple(spec.geometry.mat_ids) == ABA_MATERIAL_IDS
    assert spec.geometry.boundaries == (BC.reflective,)
    assert spec.n_groups == 2 and spec.n_regions == 3
    assert all(len(m.SigS) == 1 for m in spec.materials.values()), "the specification's materials are isotropic"
    assert aba_specification(coord) is spec  # one definition, built once


def test_r7b2_1_the_library_mixtures_carry_a_p1_moment_the_problem_drops() -> None:
    """Activation of the isotropic rebuild: the library's moderator carries a P1 moment (mean cosine 0.6), so a
    specification built from ``get_mixture`` would pose another problem."""
    assert len(get_mixture("B", "2g").SigS) == 2 and np.any(get_mixture("B", "2g").SigS[1].toarray() != 0.0)


@pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
def test_r7b2_1_the_snapshot_fixtures_mesh_the_specifications_body(coord) -> None:
    from tests.gates.sn.regression import _generate_snapshots as generator

    config = (generator._sphere_3region("2g", 40) if coord is CoordSystem.SPHERICAL
              else generator._cylinder_3region("2g", 40, "folded_4x8"))
    mesh, spec = config["mesh"], aba_specification(coord)
    assert mesh.coord is coord
    assert set(spec.geometry.breakpoints) <= set(float(e) for e in mesh.edges)
    assert mesh.region_materials == tuple(spec.geometry.mat_ids)
    assert mesh.outer_law == BC.reflective
    for mid, mixture in config["materials"].items():
        for a, b in zip(_transport_data(mixture), _transport_data(spec.materials[mid])):
            np.testing.assert_array_equal(a, b)


def test_r7b2_1_the_cylinder_rows_mesh_the_specifications_body() -> None:
    """The standoff and unified files' cylinder meshes and their ``_make_2g_mixture`` rebuilds."""
    from tests.gates.sn.sweep.curvilinear.test_unified_matvec_cylinder import _build_mr_cylinder_mesh
    from tests.gates.sn.verification.analytical.test_l1_standoff_slab_cylinder import _build_cyl_mesh

    spec = aba_specification(CoordSystem.CYLINDRICAL)
    for mesh, materials in (_build_mr_cylinder_mesh(40), _build_cyl_mesh(40)):
        assert set(spec.geometry.breakpoints) <= set(float(e) for e in mesh.edges)
        np.testing.assert_array_equal(mesh.mat_ids[[0, mesh.N // 2, -1]], [0, 1, 0])
        for mid in set(mesh.mat_ids.tolist()):
            for a, b in zip(_transport_data(materials[mid]), _transport_data(aba_materials()[mid])):
                np.testing.assert_array_equal(a, b)
