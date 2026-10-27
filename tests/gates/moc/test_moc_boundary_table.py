r"""MoC's boundary table: ``MOC_REALIZERS``, the preset's declaration and the carry (wave 2, WF2.3 slice S4; plan sections 2.8, 4.2, 6.4).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 4.2 (the MoC column), 4.3 delta 9, 6.4. MoC is NOT improved on its own terms
(``CLAUDE.md``, "The direction of development"): each entry is today's behaviour re-keyed by type, a harmonisation onto the shared
machinery. The MoC rows of the shared admission pin carry the grid; THIS file carries the table's identity, the preset and the carry.

What MoC computes, and why the specular wall is APPROXIMATED
============================================================

``MOCMesh`` reads a solid cylindrical ``Mesh1D`` as the Wigner-Seitz cylinder of a square pin cell: region 0 is a disk, the pitch is
``edges[-1] * sqrt(pi)`` and the rays are traced across the SQUARE, linked by specular reflection at its flat faces (the registry's one
kind, ``reflective``, applied to every track). A specular wall of weight 1 on the declared cylinder is therefore realized on a different
shape: ``Approximated.of("an equal-area square cell with reflecting flat faces, not the declared cylinder", 547)`` (#547, section 4.2;
the name is pinned character for character by the shared admission pin). Every other wall is ``Refused(PLUMBING, #654)``: the registry
held nothing else (``BC_REGISTRY = {"reflective": ...}``). A slab, a hollow body and any non-cylinder are refused by MoC's own guard
(#514, ``SCOPE-BOUNDARY[guard]``), not by admission, and the migrated ``tests/gates/moc/test_moc_dropped_laws_are_refused.py::TestMOC`` rows
keep them.

The ladder
==========

::

    test_the_moc_table_is_sealed_and_the_old_registry_is_gone    foundation: the table exists; BC_REGISTRY, bc_kind and the factory do not
        ^
    test_the_preset_declares_a_unit_specular_wall                the preset's cylinder carries tag 0 and solve_moc binds it to specular(1.0)
    test_an_explicit_condition_overrides_the_preset              the binding is the caller's when given (a vacuum is refused, #654)
        ^
    test_the_carry_of_the_moc_solve                              bit-equal to the HEAD capture (cylinder, mirror, 5 + 5 cells)

First reds at the base tip ``d4da4328``: ``ImportError`` from ``orpheus.moc.geometry`` and ``orpheus.moc.solver`` (the retired ``BC``) for
every row but the fixture row (file read only). Mutations for the reviewers: M-bc-kind-restored (``bc_kind`` and ``BC_REGISTRY`` put back);
M-preset-vacuum (``default_pin_cell_mesh``'s tag bound to ``vacuum()``: the preset row refuses at build); M-links-off (the track linking
skipped: the carry row moves by O(1)).
"""

from __future__ import annotations

import math
from functools import cache
from typing import Any

import numpy as np
import pytest

from orpheus.moc import geometry, solver
from tests.gates.transport import _carry_cases as carry
from tests.gates.moc import _moc_s4_helpers as MOC
from tests.gates.transport import _wave2_s4_helpers as H

_THIS = "tests/gates/moc/test_moc_boundary_table.py"
_TABLE = f"{_THIS}::test_the_moc_table_is_sealed_and_the_old_registry_is_gone"
_FIXTURE_TEST = "tests/gates/sn/architecture/test_wave2_carry.py::test_the_fixture_is_the_head_capture"


@pytest.mark.foundation
def test_the_moc_table_is_sealed_and_the_old_registry_is_gone() -> None:
    r"""``MOC_REALIZERS`` is a sealed ``RealizerTable`` of method ``"moc"`` on the non-negative-flux cone whose unregistered-unit door cites
    #654; ``MOCMesh.BC_REGISTRY``, ``MOCMesh._resolve_bc``, the instance attribute ``bc_kind`` (0 readers, ``_census_g24.out``) and the
    module function ``_moc_bc_reflective`` are gone in the same commit.

    Red input at the base tip: ``orpheus.moc.geometry`` cannot be imported (the retired ``BC``); once it can, the table does not exist
    while ``BC_REGISTRY`` does. Mutation: M-bc-kind-restored.
    """
    B = H.bc()
    table = MOC.moc_table()
    assert isinstance(table, B.RealizerTable)
    assert table.method == "moc" and table.state_cone is B.StateCone.NONNEGATIVE_FLUX and table.unregistered_issue == 654
    table.catalogue()
    assert not hasattr(geometry.MOCMesh, "BC_REGISTRY") and not hasattr(geometry.MOCMesh, "_resolve_bc")
    assert not hasattr(geometry, "_moc_bc_reflective"), "the registry's factory is retired with the registry"
    built = MOC.moc_head(MOC.moc_mesh1d(), {0: B.WallCondition.specular(1.0)})
    assert not hasattr(built, "bc_kind"), "bc_kind has 0 readers and retires (section 2.8)"


@pytest.mark.l1
@pytest.mark.rests_on(_TABLE)
def test_the_preset_declares_a_unit_specular_wall() -> None:
    r"""``solve_moc``'s default mesh (``default_pin_cell_mesh``) is the Wigner-Seitz radii cylinder with its outer face the TAG 0 (section
    6.4: ``moc/solver.py:97`` declares ``cylinder(..., outer=0)``), and a call with no mesh and no ``boundary_conditions`` binds that tag
    to ``WallCondition.specular(1.0)`` itself: the preset is a specular wall, the one wall MoC realizes, so ``solve_moc(materials)``
    still runs. The mesh read alone carries the tag and no condition (a tagged carrier must say what its tags mean).

    The solve is coarse (4 azimuthal angles, 1 polar, ray spacing 0.4): the row reads that the call is admitted and returns a finite
    positive ``k``, not a value.

    Red input at the base tip: ``ImportError``. Mutation: M-preset-vacuum.
    """
    mesh = solver.default_pin_cell_mesh()
    assert dict(mesh.face_roles) == {"xmax": 0}, f"the preset's roles are {dict(mesh.face_roles)}, not the tag 0 on the outer face"
    from orpheus.data.macro_xs.synthetic import get_mixture

    materials = {2: get_mixture("A", "2g"), 1: get_mixture("B", "2g"), 0: get_mixture("B", "2g")}
    with H.quiet():
        result = solver.solve_moc(materials, n_azi=4, n_polar=1, ray_spacing=0.4, max_outer=3)
    assert math.isfinite(result.keff) and result.keff > 0.0


@pytest.mark.l1
@pytest.mark.rests_on(_TABLE)
def test_an_explicit_condition_overrides_the_preset() -> None:
    r"""The condition a caller passes wins over the preset's: a vacuum on the preset's tag is ``Refused(PLUMBING, #654)`` (MoC links only
    specular tracks), so the same call that succeeds with the preset is refused with the caller's declaration. Control for the preset row:
    it shows the preset is a DEFAULT and not the only binding the solve reads.

    Red input: ``ImportError`` at the base tip. Mutation: a ``solve_moc`` that ignores the keyword (the vacuum row builds).
    """
    B = H.bc()
    from orpheus.data.macro_xs.synthetic import get_mixture

    materials = {2: get_mixture("A", "2g"), 1: get_mixture("B", "2g"), 0: get_mixture("B", "2g")}
    with pytest.raises(B.BoundaryRefused) as info, H.quiet():
        solver.solve_moc(materials, n_azi=4, n_polar=1, ray_spacing=0.4, max_outer=3, boundary_conditions={0: B.WallCondition.vacuum()})
    assert info.value.admission.axis is B.RefusalAxis.PLUMBING and info.value.admission.issue == 654


@cache
def _fixture() -> Any:
    assert carry.CARRY_FIXTURE.exists(), f"{carry.CARRY_FIXTURE} is missing: the carry fixture is the instrument of this file"
    return np.load(carry.CARRY_FIXTURE)


@pytest.mark.foundation
@pytest.mark.rests_on(_FIXTURE_TEST)
def test_the_fixture_holds_the_moc_row() -> None:
    r"""The fixture holds ``keff_moc/cylinder/mirror`` captured at ``5d04c296`` and it is physical: MoC's specular cylinder (an equal-area
    square cell) is within 5e-3 of CP's white-closed cylinder on the same body (both close the leakage on the same equal-area geometry by
    different closures; ``0.73586`` against ``0.73491`` at capture), so a flat or foreign column fails here. This is a code-to-code
    SANITY of the capture (L4, no correctness claim), the control that the number was captured on the body the gate rebuilds.

    Red input: none today; it reddens on a truncated fixture or one captured from another tree.
    """
    fixture = _fixture()
    assert str(fixture["meta/commit"]) == carry.HEAD_COMMIT
    moc, cp_white = H.fixture_k("moc", "cylinder", ("mirror",)), H.fixture_k("cp", "cylinder", ("white", 1.0))
    assert abs(moc - cp_white) < 5e-3, f"MoC {moc!r} against CP white {cp_white!r}"


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_the_fixture_holds_the_moc_row")
def test_the_carry_of_the_moc_solve() -> None:
    r"""``solve_moc`` in the new spelling (the cylinder's outer face a tag bound to ``WallCondition.specular(1.0)``, rule R2) equals the
    HEAD ``k`` bit for bit: cylinder ``(0, 1, 2)``, materials ``(1, 0)``, 5 + 5 cells, ``n_azi=8, n_polar=2, ray_spacing=0.1``.

    Red input: ``ImportError`` at the base tip. Mutations: M-links-off; ``specular(1.0)`` realized with another return weight.
    """
    k = MOC.moc_solve(MOC.moc_mesh1d(), {0: H.bc().WallCondition.specular(1.0)}).keff
    expected = H.fixture_k("moc", "cylinder", ("mirror",))
    assert k == expected, f"k = {k!r}, the HEAD k = {expected!r}"
