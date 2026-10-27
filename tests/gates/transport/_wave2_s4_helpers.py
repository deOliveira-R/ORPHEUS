r"""The method-agnostic builders of the S4 gates on the CP, MoC and MC boundary tables (wave 2, WF2.3 slice S4).

Not a test module (no ``test_`` prefix). Spec: ``.claude/plans/wave2_boundary_carve.md`` sections 2.8 (the three tables), 4.2 (the CP, MoC
and MC columns), 6.0 (the carry fixture) and 6.4 (the seven done-when gates).

Each method's builders live beside its gates since 2026-10-26: ``tests/gates/cp/_cp_s4_helpers.py``, ``tests/gates/moc/_moc_s4_helpers.py``
and ``tests/gates/mc/_mc_s4_helpers.py`` (the tables' accessors, the carried meshes, solves and heads). They were here until layer DAG
ruling 27 turned their string imports into statements, which made this transport helper a gate importing four method packages; the
gates that span the three tables moved to ``tests/gates/cross_method/``. What stays is what every method's gates share:

* the faces' :class:`~orpheus.transport.boundary_conditions.FaceSite` of the section 4.2 site classes, CP's slab split by face
  (:func:`site`: ``CART_XMIN``, ``CART_XMAX``, ``RAD_OUTER``, ``RAD_INNER``), because CP's table is keyed on ``(site_class, face)``;
* :func:`white_one`, :func:`quiet`, and :func:`fixture_k`, the frozen HEAD eigenvalue of the carry fixture for one carried row (the
  carried problems are ``tests/gates/transport/_carry_cases.py``'s, with the SAME charts, breakpoints, materials, cell counts, solver
  settings and seed the HEAD capture script used, because a bit-equal ``k`` is meaningful only for the same problem);
* the shared chart builder and admission summary of ``tests.gates.diffusion._wave2_helpers`` (re-exported: one builder, Pattern 2).
"""

from __future__ import annotations

import contextlib
import io
import types
from typing import Any

import numpy as np

from orpheus.transport import boundary_conditions
from tests.gates.diffusion._wave2_helpers import chart_mesh, summary  # noqa: F401  (re-exported: the one builder, the one summary)
from tests.gates.transport import _carry_cases as carry


def bc() -> types.ModuleType:
    """``orpheus.transport.boundary_conditions``."""
    return boundary_conditions


def site(label: str) -> Any:
    """The sites of section 4.2 as ``FaceSite`` values. ``CART_XMIN`` and ``CART_XMAX`` are the Cartesian faces CP tells apart."""
    from orpheus.geometry.coord import CoordSystem

    B = bc()
    match label:
        case "CART_XMIN":
            return B.FaceSite("xmin", CoordSystem.CARTESIAN, 1)
        case "CART_XMAX" | "CART":
            return B.FaceSite("xmax", CoordSystem.CARTESIAN, 1)
        case "RAD_OUTER":
            return B.FaceSite("xmax", CoordSystem.CYLINDRICAL, 1)
        case "RAD_INNER":
            return B.FaceSite("xmin", CoordSystem.CYLINDRICAL, 1)
    raise ValueError(label)


def white_one() -> Any:
    """``WallCondition.white(1.0)``: the diffuse chain at unit weight, the one weight CP realizes (section 4.3 delta 3)."""
    return bc().WallCondition.white(1.0)


@contextlib.contextmanager
def quiet():
    """The solvers print a banner; the gates read numbers, not stdout."""
    with contextlib.redirect_stdout(io.StringIO()):
        yield


def fixture_k(method: str, chart: str, declaration: tuple) -> float:
    """The frozen HEAD eigenvalue of ``carry_head.npz`` for one carried row."""
    fixture = np.load(carry.CARRY_FIXTURE)
    return float(fixture[carry.other_key(method, chart, declaration)])
