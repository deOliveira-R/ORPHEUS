r"""A vacuum-bounded CP eigenvalue must not exceed the white-bounded one on the same body (#546, a defect left visible by ruling D4).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 6.4 ("The #546 strict xfail, spelled out") and section 4.3 delta 2.

The claim and its status
========================

``CPSolver.compute_keff`` divides production by net removal with no leakage term, a premise true only under a white wall, so a
vacuum-bounded eigenvalue lies ABOVE the white-bounded one on the same body, which physics forbids (a leaky wall cannot raise
``k``): ``k_vacuum <= k_white``. The user's ruling of 2026-10-10 leaves CP's pencil to a later campaign, and decision D4 of the twelfth
exchange forbids masking the defect: the carve admits CP times the empty network (vacuum) as ``Realized(faithful=True)`` (section 4.2),
so the wall's admission is silent about the defect, and THIS file is where it stays visible. It states the physics a reader can check
without CP's internals (the sibling ``test_cp_first_collision.py::test_the_vacuum_k_is_the_eigenvalue_of_the_cp_pencil`` compares CP
with its own dense pencil and stays beside it).

The ladder
==========

::

    test_the_control_sn_satisfies_the_same_inequality[...]            REFERENCE: the same comparison through SN holds on the same body
        ^                                                              (slab: vacuum vs white; curved: vacuum vs the unit specular wall,
        |                                                               because SN refuses the white wall on a radial face, #653)
    test_a_vacuum_bounded_cp_eigenvalue_does_not_exceed_the_white_...  RECORD of the defect: strict xfail on 3 of 3 charts
    test_the_vacuum_above_white_rows_are_strict_xfails                 foundation: the three marks are strict, raise AssertionError, cite #546

Deviation from the literal code of section 6.4, and why (``instrument-doctrine`` X1)
====================================================================================

Section 6.4 spells the marks ``xfail(strict=True, reason=...)``. A bare ``xfail`` swallows ANY failure of the body, including the
``ImportError`` the body raises at the base tip (every ``orpheus.cp`` import raises, slice S1 removed ``BC``) and the ``TypeError`` of
the missing ``boundary_conditions`` keyword before the carve's keyword exists: the three rows would read XFAILED at the base tip for a
reason that has nothing to do with #546 (``vv-principles`` mode 8, arm 4: a misattributed strict xfail). The marks here carry
``raises=AssertionError``, so only the inequality failing is the expected failure, and the rows are FAILED (red) at the base tip and
XFAILED after the carve. The companion foundation test asserts ``raises`` as well as ``strict`` and the issue number.

The numbers (``[M]``, section 6.4: ``scratch/layer_dag/wave2_spec/_probe_546_xfail.py`` at ``5d04c296``): slab 1.4003162 against
1.2308420 (gap 0.1695), cylinder 0.9354010 against 0.7393875 (gap 0.1960), sphere 0.5854991 against 0.4140867 (gap 0.1714), on the 2G
``A|B`` body, uniform_volume(4) + (4), ``CPParams(keff_tol=1e-11)``.

First reds at the base tip ``d4da4328`` (the command is in the report): the three xfail rows FAIL with ``ModuleNotFoundError``/
``ImportError`` (``orpheus.cp.solver`` imports the retired ``BC``) instead of XFAIL; the control rows pass (SN is carved); the
foundation row passes (it reads the marks only). Mutation for the reviewers: ``compute_keff`` given the leakage term (the three
rows XPASS, strictness reddens them: the day the defect is fixed the marks are deleted).

Since 2026-10-26 ``orpheus.cp.solver`` is imported at the top (layer DAG ruling 27: an import whose target is known when the test is
written is an import statement), so a module that cannot be imported is a collection error, itself a red; ``raises=AssertionError``
still keeps every other failure of the body from reading as the #546 record. The #546 body (:func:`_mesh_546`, :func:`_mixtures_546`)
lives here, its one consumer: the 2G ``A|B`` body ``cp_body(coord, (0.0, 0.5, 1.0), (2, 0), outer)`` meshed by ``uniform_volume(4)``
per region, as ``tests/gates/cp/test_cp_first_collision.py:_vacuum_mesh`` builds it.
"""

from __future__ import annotations

from typing import Any

import pytest

from orpheus.cp import solver as cp_solver
from orpheus.data.macro_xs.synthetic import get_mixture
from orpheus.geometry import CoordSystem
from orpheus.mesh import CellsByCount, Mesher
from tests.gates.cp._cp_geometry import cp_body
from tests.gates.transport import _wave2_s4_helpers as H

_CHARTS = ("slab", "cylinder", "sphere")
_REASON = (
    "#546: CPSolver.compute_keff divides production by net removal with no leakage term, so a vacuum-bounded "
    "eigenvalue lies above the white-bounded one on the same body (measured 2026-10-22 on the 2G A|B body, cp_body, "
    "uniform_volume(4)+(4), CPParams(keff_tol=1e-11): slab 1.4003162 vs 1.2308420 (gap 0.1695), cylinder 0.9354010 vs "
    "0.7393875 (gap 0.1960), sphere 0.5854991 vs 0.4140867 (gap 0.1714)); physics requires k_vacuum <= k_white. "
    "Not fixed by the user's ruling of 2026-10-10 (CP's pencil is a later campaign); D4 of the twelfth exchange: no masking."
)
_VACUUM_ABOVE_WHITE = [
    pytest.param(chart, id=chart, marks=pytest.mark.xfail(strict=True, raises=AssertionError, reason=_REASON)) for chart in _CHARTS
]
_THIS = "tests/gates/cp/test_cp_vacuum_not_above_white.py"


def _mesh_546(chart: str, outer: Any = 0) -> Any:
    """The 2G ``A|B`` body of the #546 rows: ``cp_body(coord, (0.0, 0.5, 1.0), (2, 0), outer)``, 4 + 4 equal-volume cells."""
    coord = {"slab": CoordSystem.CARTESIAN, "cylinder": CoordSystem.CYLINDRICAL, "sphere": CoordSystem.SPHERICAL}[chart]
    return Mesher(cp_body(coord, (0.0, 0.5, 1.0), (2, 0), outer)).partition((CellsByCount.uniform_volume(4), CellsByCount.uniform_volume(4))).mesh


def _mixtures_546() -> dict[int, Any]:
    """``{2: A, 0: B}`` of the 2-group library, the materials of the #546 body."""
    return {2: get_mixture("A", "2g"), 0: get_mixture("B", "2g")}


def _k(chart: str, condition) -> float:
    """``solve_cp``'s k on the #546 body (``cp_body``, 2G ``A|B``, 4 + 4 cells) with the outer tag bound to ``condition``."""
    params = cp_solver.CPParams(keff_tol=1e-11)
    with H.quiet():
        return cp_solver.solve_cp(_mixtures_546(), _mesh_546(chart, 0), params, boundary_conditions={0: condition}).keff


@pytest.mark.l1
@pytest.mark.parametrize("chart", _VACUUM_ABOVE_WHITE)
def test_a_vacuum_bounded_cp_eigenvalue_does_not_exceed_the_white_bounded_one(chart: str) -> None:
    r"""``k_vacuum <= k_white`` on the 2G ``A|B`` body, through CP. Strict xfail on 3 of 3 charts (the defect is #546).

    Red input at the base tip: the body raises ``ImportError`` (not ``AssertionError``), so ``raises=AssertionError`` makes the row
    FAIL rather than XFAIL. After the carve the inequality fails by the measured gaps and the row XFAILs; when #546 is fixed it XPASSes
    and strictness reddens it. Mutation: the measured defect fixed (leakage in ``compute_keff``) in a copy.
    """
    k_vac, k_white = _k(chart, H.bc().WallCondition.vacuum()), _k(chart, H.white_one())
    assert k_vac <= k_white, f"{chart}: k_vacuum {k_vac!r} exceeds k_white {k_white!r} by {k_vac - k_white:.6e}"


def _sn_k(chart: str, condition) -> float:
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.solver import solve_sn

    quadrature = Quadrature.folded_product(n_mu=4, n_phi=8) if chart == "cylinder" else Quadrature.gauss_legendre(n_ordinates=8)
    return solve_sn(_mixtures_546(), _mesh_546(chart, 0), quadrature, boundary_conditions={0: condition}).outcome.keff


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_the_vacuum_above_white_rows_are_strict_xfails")
@pytest.mark.parametrize("chart", _CHARTS)
def test_the_control_sn_satisfies_the_same_inequality(chart: str) -> None:
    r"""The assertion CAN pass: the same comparison through SN holds on the same body. Slab: ``k_vacuum < k_white`` (``[M]`` 2026-10-23
    at ``d4da4328``: 0.15554 against 1.22974). The radial rows compare vacuum with the unit specular wall (0.04111 against 0.74726 on the
    cylinder, 0.01867 against 0.42896 on the sphere) because SN refuses the white wall on a radial outer face
    (``Refused(PLUMBING, #653)``, section 4.2): section 6.4's text says "all three charts" for the white wall and is refuted by that
    cell, so the control is the nearest wall SN realizes there. The pair is a REFERENCE (a different method, a different carrier), so it
    shares no code with CP's ``compute_keff``.

    Red input: none at the base tip (SN is carved). The control reddens on an SN whose vacuum ``k`` exceeds its return-wall ``k``
    (a vacuum realized as a mirror: both ``k`` equal, ``<`` fails).
    """
    W = H.bc().WallCondition
    k_vac = _sn_k(chart, W.vacuum())
    k_return = _sn_k(chart, W.white(1.0) if chart == "slab" else W.specular(1.0))
    assert k_vac < k_return, f"SN {chart}: k_vacuum {k_vac!r} is not below the return-wall k {k_return!r}"


@pytest.mark.foundation
def test_the_vacuum_above_white_rows_are_strict_xfails() -> None:
    r"""The three rows carry one ``xfail`` each, ``strict=True``, ``raises=AssertionError``, the reason naming #546: a mark that lost
    ``strict`` inside ``pytest.param(marks=...)`` would turn the repair into a silent XPASS, and one without ``raises`` would
    swallow an ``ImportError`` as the expected failure (the deviation from section 6.4's literal code, stated above).

    Red input: none (the marks are read, not run); it reddens on a mark that lost any of the three properties.
    """
    assert len(_VACUUM_ABOVE_WHITE) == 3
    for param in _VACUUM_ABOVE_WHITE:
        xfails = [m for m in param.marks if getattr(m, "name", None) == "xfail"]
        assert len(xfails) == 1 and xfails[0].kwargs.get("strict") is True, param.id
        assert xfails[0].kwargs.get("raises") is AssertionError, param.id
        assert "#546" in xfails[0].kwargs.get("reason", ""), param.id
