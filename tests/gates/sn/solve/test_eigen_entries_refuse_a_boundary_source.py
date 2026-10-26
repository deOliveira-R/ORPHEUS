"""An eigenvalue question has no source term, so a declared boundary source is refused.

``solve_sn`` and ``solve_sn_adjoint`` pose the homogeneous pencil
:math:`A\\psi = \\mu F\\psi`. A face declared with a ``WallCondition.uniform(q)``
carries a source :math:`q \\neq 0`, which that question does not contain.

The first red (qa review of ERR-094, 2026-09-30): before the refusal both
entries answered the source-free problem and dropped the declared :math:`q`
in silence; ``[M]`` k was bitwise the vacuum k for q = 1 and q = 50. A
source belongs to a source question (``solve_sn_fixed_source``,
``solve_sn_multiplying_source``).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.data.macro_xs.synthetic import get_mixture
from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn, solve_sn_adjoint
from orpheus.transport.boundary_conditions import WallCondition


_MATERIALS = {0: get_mixture("A", "2g")}
_QUADRATURE = Quadrature.gauss_legendre(8)


def _slab():
    """The slab whose left face carries tag 0 and whose right face carries tag 1; the left condition is the test's variable."""
    geometry = StructuredGeometry.slab((0.0, 2.0), (0,), left=0, right=1)
    return Mesher(geometry).partition(CellsByCount.uniform_width(8)).mesh


@pytest.mark.l0
@pytest.mark.catches("ERR-095")
@pytest.mark.parametrize("entry", [solve_sn, solve_sn_adjoint], ids=lambda f: f.__name__)
def test_a_declared_boundary_source_is_refused(entry) -> None:
    with pytest.raises(ValueError, match=r"eigenvalue question.*\['xmin'\].*prescribed inflow"):
        entry(_MATERIALS, _slab(), _QUADRATURE, boundary_conditions={0: WallCondition.uniform(1.0), 1: WallCondition.vacuum()})


@pytest.mark.l0
@pytest.mark.parametrize("entry", [solve_sn, solve_sn_adjoint], ids=lambda f: f.__name__)
def test_a_source_free_declaration_is_admitted(entry) -> None:
    """The control: the same slab with a vacuum left face solves."""
    result = entry(_MATERIALS, _slab(), _QUADRATURE, boundary_conditions={0: WallCondition.vacuum(), 1: WallCondition.vacuum()})
    assert np.isfinite(result.outcome.keff)
