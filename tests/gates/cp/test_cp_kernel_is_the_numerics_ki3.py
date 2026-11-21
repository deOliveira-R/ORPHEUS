r"""G5b of the unit ``insulation`` (#641 item 1): the cylinder ``CPMesh`` binds the kernel ``orpheus.numerics.bickley.ki3``.

The binding leg lives here and not in ``tests/gates/derivations/test_cp_ki3_copy.py``: a ``tests/gates/derivations/`` file may not import
``orpheus.cp`` (an L3 package). The cylinder ``CPMesh`` (built as in ``test_cp_chord_spelling``, with ``{0: WallCondition.vacuum()}``) binds
``ki3`` (``mesh._kernel is ki3``), ``mesh._kernel_zero`` equals ``ki3(np.zeros(n_y))`` bit for bit, and the sphere binds the exponential kernel.

First red at the base tip: ``ModuleNotFoundError`` for ``orpheus.numerics.bickley`` (the kernel is the reference's ``_ki3_mp`` there).
Mutation: ``monkeypatch.setattr(orpheus.cp.solver, "ki3", lambda tau: np.exp(-tau))`` before building the cylinder ``CPMesh`` (the production
module reads the name at call time): the identity row and the zero-row equality redden. The sphere row is the control that is NOT
moved by it.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.transport.boundary_conditions import WallCondition
from tests.gates.cp import _cp_s4_helpers as CP


def _mesh(chart: str):
    return CP.cp_head(CP.cp_mesh(chart), {0: WallCondition.vacuum()})


@pytest.mark.l1
def test_the_cylinder_binds_the_numerics_bickley_kernel() -> None:
    from orpheus.numerics.bickley import ki3

    mesh = _mesh("cylinder")
    assert mesh._kernel is ki3
    assert np.array_equal(mesh._kernel_zero, ki3(np.zeros(len(mesh._y_pts))))


@pytest.mark.l1
def test_the_sphere_binds_the_exponential_kernel() -> None:
    """The sphere's kernel is ``exp(-tau)``: ``_kernel_zero`` is ones and the kernel agrees with ``exp`` on a grid. Green at the base tip."""
    mesh = _mesh("sphere")
    tau = np.linspace(0.0, 5.0, 11)
    assert np.array_equal(mesh._kernel(tau), np.exp(-tau)) and np.array_equal(mesh._kernel_zero, np.ones(len(mesh._y_pts)))


@pytest.mark.l1
def test_the_binding_row_reddens_when_the_kernel_name_is_replaced(monkeypatch: pytest.MonkeyPatch) -> None:
    """The mutation: ``orpheus.cp.solver.ki3`` replaced by ``exp(-tau)`` before the cylinder is built; the cylinder then binds a kernel that
    is not ``numerics.bickley.ki3``, and its zero row is ones instead of ``Ki_3(0)``."""
    from orpheus.cp import solver

    from orpheus.numerics.bickley import ki3

    monkeypatch.setattr(solver, "ki3", lambda tau: np.exp(-tau))
    mesh = _mesh("cylinder")
    assert mesh._kernel is not ki3
    assert not np.array_equal(mesh._kernel_zero, ki3(np.zeros(len(mesh._y_pts))))
