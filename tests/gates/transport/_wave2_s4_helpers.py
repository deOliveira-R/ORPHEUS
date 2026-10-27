r"""Shared builders of the S4 gates on the CP, MoC and MC boundary tables (wave 2, WF2.3 slice S4).

Not a test module (no ``test_`` prefix). Spec: ``.claude/plans/wave2_boundary_carve.md`` sections 2.8 (the three tables), 4.2 (the CP, MoC
and MC columns), 6.0 (the carry fixture) and 6.4 (the seven done-when gates).

How this module reads the production modules
============================================

The cp, moc, mc and diffusion modules are imported at the top, as MODULES, and the accessors below read their attributes at call
time, so a monkeypatch of a module attribute reaches every gate (layer DAG ruling 27: an import whose target is known when the test
is written is an import statement, which the test-layer gate sees). Until 2026-10-26 each accessor imported its module by string,
because at the base tip of slice S4 every cp, moc and mc import raised ``ImportError`` (each solver module imported the retired
``BC``); that reason ended with the slice. The statements make visible what the strings hid: a transport helper importing four
method packages, listed in ``OPEN_TEST_LAYER_VIOLATIONS`` (``tests/gates/test_test_layers.py``).

Where the tables live (the contract the gates read)
===================================================

Section 6.4 lists the production files of the slice as ``orpheus/cp/solver.py``, ``orpheus/mc/solver.py`` and
``orpheus/moc/{geometry.py, solver.py}``, and section 2.8 says the three tables are ``CP_REALIZERS``, ``MOC_REALIZERS`` and
``MC_REALIZERS``, "each in its method's package": the tables replace the three ``BC_REGISTRY`` dicts, which live today in
``orpheus.cp.solver``, ``orpheus.moc.geometry`` and ``orpheus.mc.solver``. The accessors :func:`cp_table`, :func:`moc_table` and
:func:`mc_table` read those modules. A builder who homes a table in another module moves the three accessors (and says so in its
report): the gates read nothing else of the location.

What it builds
==============

* the faces' :class:`~orpheus.transport.boundary_conditions.FaceSite` of the section 4.2 site classes, CP's slab split by face
  (:func:`site`: ``CART_XMIN``, ``CART_XMAX``, ``RAD_OUTER``, ``RAD_INNER``), because CP's table is keyed on ``(site_class, face)``;
* the carried problems of ``tests/gates/transport/_carry_cases.py`` in the NEW spelling for CP (:func:`cp_mesh`, :func:`cp_solve`),
  MoC (:func:`moc_mesh`, :func:`moc_solve`) and MC (:func:`mc_mesh`, :func:`mc_solve`): the SAME charts, breakpoints, materials, cell
  counts, solver settings and seed the HEAD capture script used (``_capture_carry_head._capture_other_methods`` (deleted in WF2.4; recover it with ``git show 197549f6:tests/gates/transport/_capture_carry_head.py``)), because a bit-equal
  ``k`` is meaningful only for the same problem. The charts are the shared ones of ``tests.gates.diffusion._wave2_helpers.chart_mesh``
  (one builder, Pattern 2; that module imports no diffusion module at the top);
"""

from __future__ import annotations

import contextlib
import io
import types
from typing import Any

import numpy as np

from orpheus.cp import solver as cp_solver
from orpheus.diffusion import zero_flux as diffusion_zero_flux
from orpheus.mc import solver as mc_solver
from orpheus.moc import geometry as moc_geometry
from orpheus.moc import quadrature as moc_quadrature
from orpheus.moc import solver as moc_solver
from orpheus.transport import boundary_conditions
from tests.gates.diffusion._wave2_helpers import chart_mesh, summary  # noqa: F401  (re-exported: the one builder, the one summary)
from tests.gates.transport import _carry_cases as carry

#: The reason text the MoC table's one named approximation carries (section 4.2, the MoC column; pinned exactly).
MOC_DOMAIN_NAME = "an equal-area square cell with reflecting flat faces, not the declared cylinder"
#: A fragment MC's TOPOLOGY refusal must carry (section 4.3 delta 9: the torus, and the issue it is tracked under).
MC_TORUS_FRAGMENT = "torus"


def bc() -> types.ModuleType:
    """``orpheus.transport.boundary_conditions``."""
    return boundary_conditions


def cp_table() -> Any:
    """``CP_REALIZERS``, the sealed CP realizer table (``orpheus/cp/solver.py``)."""
    return cp_solver.CP_REALIZERS


def moc_table() -> Any:
    """``MOC_REALIZERS``, the sealed MoC realizer table (``orpheus/moc/geometry.py``)."""
    return moc_geometry.MOC_REALIZERS


def mc_table() -> Any:
    """``MC_REALIZERS``, the sealed MC realizer table (``orpheus/mc/solver.py``)."""
    return mc_solver.MC_REALIZERS


def zero_flux() -> Any:
    """A ``ZeroFluxCondition()`` (diffusion's method-only condition, which no CP, MoC or MC table registers)."""
    return diffusion_zero_flux.ZeroFluxCondition()


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


# ── CP ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────


def cp_mesh(chart: str, *, outer: Any = 0, left: Any = None) -> Any:
    """The carried CP body of ``chart``: a slab has a mirror at the left face (CP's kernel images one) and ``outer`` on the right; a
    solid radial body has ``outer`` on its one face. ``left`` overrides the slab's left role (a tag or a deck)."""
    from orpheus.geometry.deck import MirrorDeck

    if chart == "slab":
        return chart_mesh(chart, left=left if left is not None else MirrorDeck(), right=outer)
    return chart_mesh(chart, outer=outer)


def cp_solve(mesh: Any, conditions: dict[int, Any] | None = None, params: Any = None) -> Any:
    """``solve_cp`` at ``CPParams()`` (what the HEAD capture ran) with ``boundary_conditions=`` iff given; stdout silenced."""
    kwargs: dict[str, Any] = {} if conditions is None else {"boundary_conditions": conditions}
    with quiet():
        return cp_solver.solve_cp(carry.materials(), mesh, params if params is not None else cp_solver.CPParams(), **kwargs)


def cp_head(mesh: Any, conditions: dict[int, Any] | None = None) -> Any:
    """A ``CPMesh`` over ``mesh`` with ``boundary_conditions=`` iff given (construction only: the table is read here)."""
    return cp_solver.CPMesh(mesh) if conditions is None else cp_solver.CPMesh(mesh, boundary_conditions=conditions)


# ── MoC ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────


def moc_mesh1d(*, outer: Any = 0) -> Any:
    """The carried MoC cylinder: breakpoints ``(0, 1, 2)``, materials ``(1, 0)``, ``MOC_CELLS`` cells, ``outer`` on its one face."""
    return chart_mesh("cylinder", outer=outer, cells=carry.MOC_CELLS)


def moc_solve(mesh: Any, conditions: dict[int, Any] | None = None) -> Any:
    """``solve_moc`` at the carried ray settings with ``boundary_conditions=`` iff given; stdout silenced."""
    kwargs: dict[str, Any] = {} if conditions is None else {"boundary_conditions": conditions}
    with quiet():
        return moc_solver.solve_moc(carry.materials(), mesh, **carry.MOC_ARGS, **kwargs)


def moc_head(mesh: Any, conditions: dict[int, Any] | None = None) -> Any:
    """A ``MOCMesh`` over ``mesh`` (4 azimuthal angles, one polar) with ``boundary_conditions=`` iff given."""
    quadrature = moc_quadrature.MOCQuadrature.create(n_azi=4, n_polar=1)
    return moc_geometry.MOCMesh(mesh, quadrature) if conditions is None else moc_geometry.MOCMesh(mesh, quadrature, boundary_conditions=conditions)


# ── MC ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────


def mc_slab_mesh1d(*, left: Any, right: Any) -> Any:
    """The carried MC slab: breakpoints ``(0, 1.5, 3)``, materials ``(2, 0)`` (fuel, moderator), 6 + 6 cells; the faces' roles as given."""
    return chart_mesh("slab", left=left, right=right, mat_ids=(2, 0))


def mc_materials() -> dict[int, Any]:
    """``{2: fuel, 0: moderator}``: the carried mixtures under MC's ids (the HEAD capture's ``mc_mats``)."""
    mats = carry.materials()
    return {2: mats[1], 0: mats[0]}


def mc_head(mesh: Any, pitch: float, conditions: dict[int, Any] | None = None) -> Any:
    """An ``MCMesh`` over ``mesh`` with ``boundary_conditions=`` iff given."""
    return mc_solver.MCMesh(mesh, pitch=pitch) if conditions is None else mc_solver.MCMesh(mesh, pitch=pitch, boundary_conditions=conditions)


def mc_solve(geometry: Any) -> Any:
    """``solve_monte_carlo`` with the carried parameters (fixed seed), the geometry given; stdout silenced."""
    with quiet():
        return mc_solver.solve_monte_carlo(mc_materials(), mc_solver.MCParams(**carry.MC_PARAMS, geometry=geometry))


def fixture_k(method: str, chart: str, declaration: tuple) -> float:
    """The frozen HEAD eigenvalue of ``carry_head.npz`` for one carried row."""
    fixture = np.load(carry.CARRY_FIXTURE)
    return float(fixture[carry.other_key(method, chart, declaration)])
