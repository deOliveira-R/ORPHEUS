r"""The carry of wave 2 for diffusion: the twelve eigenvalues the carve must not move (law L-K; plan sections 6.0, 6.3 gate 5, 8).

``tests/gates/transport/carry_head.npz`` was captured FROM THE HEAD TREE (``main`` ``5d04c296``, a detached worktree, before any
edit) by ``tests/gates/transport/_capture_carry_head.py``: ``solve_diffusion_1d`` (default tolerances) on the three one-dimensional
charts, a slab ``(0, 1.5, 3)`` and a cylinder and a sphere ``(0, 1, 2)``, each of two regions (materials ``1`` then ``0``) with six
cells apiece and the two-group library, under four declarations of the outer wall: vacuum, the perfect mirror, the bare albedo
``0.5`` and zero flux (12 rows, key ``keff_diffusion/<chart>/<declaration>``). The same problems, spelled in the NEW vocabulary
(``tests/gates/diffusion/_wave2_helpers.build``), must give the same ``k`` BIT-EQUAL (``==`` on the float; a non-bit difference is a
finding to investigate, never a tolerance to widen): the slab mirror is a ``MirrorDeck`` on each face and the curved mirror a tag
bound to ``WallCondition.specular(1.0)`` (rule R2, Appendix B's meaning change); the bare albedo is
``WallCondition.white(0.5)`` (R4: the shape is invisible to diffusion, the only method that admitted the bare tag); zero flux is
``ZeroFluxCondition()`` on a tag.

What each row can and cannot see
================================

* The solve reads every layer between the declaration and ``k``: the condition, the binding, the table, the realized scalar, the
  leakage operator's face closure and the power iteration. A move of any scalar (``A``), of the face closure, or of the realized
  operator's type moves a ``k`` here; the four declarations hit ``A = 0, 1, 0.5, -1``.
* Blind: a defect that is the same for every ``A`` (the bulk operator), or that is zero on a two-region two-group body only if the
  body were homogeneous (it is not), so redistribution is in view.

Beside the carry, the end-to-end legs that the carry cannot reach, each a structural statement and not a record:

* the **tier-2 witness through the solve**: ``white(a)``, ``specular(a)`` and (``a = 1``) the mirror deck give the SAME ``k`` bit for
  bit (diffusion cannot tell them apart: they realize one operator);
* the **functor law through the solve**: a mixed wall ``from_amplitudes(s, d)`` with ``s + d != 1`` gives the ``k`` of
  ``white(s + d)`` to ``10 x keff_tol`` (iterative: ten times the solver's own tolerance, ``SAFETY(10)`` of
  ``feedback_regression_tolerance_design``), while the ``k`` of ``white(s)`` and of ``white(d)`` each differ from it by orders more
  (a mutation that realizes one chain only, or rescales, is seen).

First red at the base tip: ``ImportError`` from ``import orpheus.diffusion`` (the carve's own import chain), for every row that
executes diffusion; the instrument row (``test_the_fixture_holds_the_diffusion_rows``) reads the file only and is green today (it
pins the instrument, not the carve).

Mutations for the reviewers: M-zero-flux-positive (``+I``: reddens the three ``zero_flux`` rows); M-white-as-specular-weight-one
(the diffuse chain realized as ``1.0 * I``: reddens the ``albedo_0.5`` rows); M-mirror-deck-as-zero (``DiffusionDecks`` realizing
the mirror as ``ZeroOperator``: reddens the slab mirror row and the tier-2 witness); M-rescale and M-first-chain-only (the functor
law leg).
"""

from __future__ import annotations

from functools import cache
from typing import Any

import numpy as np
import pytest

from tests.gates.diffusion import _wave2_helpers as H
from tests.gates.transport import _carry_cases as carry

_THIS = "tests/gates/diffusion/test_diffusion_boundary_carry.py"
_FIXTURE_TEST = "tests/gates/sn/architecture/test_wave2_carry.py::test_the_fixture_is_the_head_capture"
_CASES = [(chart, d) for chart in carry.CHARTS for d in carry.DIFFUSION_DECLARATIONS]
_IDS = [f"{chart}-{carry.tag(d)}" for chart, d in _CASES]


@cache
def _fixture() -> Any:
    assert carry.CARRY_FIXTURE.exists(), f"{carry.CARRY_FIXTURE} is missing: the carry fixture is the instrument of this file"
    return np.load(carry.CARRY_FIXTURE)


@pytest.mark.foundation
@pytest.mark.rests_on(_FIXTURE_TEST)
def test_the_fixture_holds_the_diffusion_rows() -> None:
    r"""The fixture holds the 12 ``keff_diffusion/*`` rows (3 charts x 4 declarations) captured at ``5d04c296``, and they are
    physical: on each chart ``k`` is ordered vacuum < albedo 0.5 < mirror (leakage falls as the return rises), and zero flux, the
    Dirichlet idealization ``J^- = -J^+``, leaks the most (the lowest ``k``). The ordering is the positive control of the
    capture: a fixture of a flat or inverted column fails here.

    Red input: none today (the capture is complete); it reddens on a truncated fixture or one captured from another tree.
    """
    fixture = _fixture()
    assert str(fixture["meta/commit"]) == carry.HEAD_COMMIT
    keys = [carry.other_key("diffusion", chart, d) for chart, d in _CASES]
    assert len(keys) == 12 and all(k in fixture.files for k in keys), [k for k in keys if k not in fixture.files]
    for chart in carry.CHARTS:
        k = {d[0]: float(fixture[carry.other_key("diffusion", chart, d)]) for d in carry.DIFFUSION_DECLARATIONS}
        assert k["zero_flux"] < k["vacuum"] < k["albedo"] < k["mirror"], f"{chart}: {k}"


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_the_fixture_holds_the_diffusion_rows")
@pytest.mark.parametrize("chart,declaration", _CASES, ids=_IDS)
def test_the_carry_of_the_diffusion_solves(chart: str, declaration: tuple) -> None:
    r"""``solve_diffusion_1d`` in the new spelling equals the HEAD ``k`` bit for bit, on the chart and declaration of the row.

    Red input: ``ImportError`` at the base tip (every row); the zero-flux rows additionally need ``orpheus.diffusion.zero_flux``.
    Mutations: M-zero-flux-positive, M-white-as-specular-weight-one, M-mirror-deck-as-zero.
    """
    mesh, conditions = H.build(chart, declaration)
    k = H.solve(mesh, conditions).keff
    expected = _fixture()[carry.other_key("diffusion", chart, declaration)]
    assert k == float(expected), f"{chart} / {carry.tag(declaration)}: k = {k!r}, the HEAD k = {float(expected)!r}"


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_the_fixture_holds_the_diffusion_rows")
@pytest.mark.parametrize("chart", list(carry.CHARTS))
def test_white_specular_and_the_mirror_are_one_image_through_the_solve(chart: str) -> None:
    r"""The tier-2 witness, end to end: diffusion annihilates the shape, so ``white(a)`` and ``specular(a)`` solve to the same
    ``k`` bit for bit (they realize one operator), and at ``a = 1`` both equal the mirror (a ``MirrorDeck`` on the slab faces, the
    unit specular tag on the curved one). The ``k`` of ``a = 0.5`` differs from the mirror's (the control: the row can fail).

    Red input: ``ImportError`` at the base tip. Mutation: M-white-as-specular-weight-one (the white row reads 1.0), a diffuse
    realizer that reads the shape.
    """
    W = H.bc().WallCondition
    cases = {name: H.build(chart, d) for name, d in (("mirror", ("mirror",)), ("albedo_half", ("albedo", 0.5)))}
    mirror = H.solve(*cases["mirror"]).keff
    roles, _ = H.roles_and_conditions(chart, None)
    for name, wall in (("white(1)", W.white(1.0)), ("specular(1)", W.specular(1.0))):
        mesh = H.chart_mesh(chart, **roles)
        assert H.solve(mesh, {0: wall}).keff == mirror, f"{chart}: {name} is not the mirror"
    white = H.solve(H.chart_mesh(chart, **roles), {0: W.white(0.5)}).keff
    specular = H.solve(H.chart_mesh(chart, **roles), {0: W.specular(0.5)}).keff
    assert white == specular, f"{chart}: white(0.5) and specular(0.5) differ: {white!r} vs {specular!r}"
    assert white != mirror and white == H.solve(*cases["albedo_half"]).keff


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_white_specular_and_the_mirror_are_one_image_through_the_solve")
def test_a_mixed_wall_solves_as_the_sum_of_scalars() -> None:
    r"""The functor law through the solve, on the slab: ``from_amplitudes(0.3, 0.5)`` gives the ``k`` of ``white(0.8)`` within
    ``10 x keff_tol`` (``1e-9``: the solver's default ``keff_tol`` is ``1e-10``), while the ``k`` of ``white(0.3)`` and of
    ``white(0.5)`` each lie more than ``1e-3`` away (``s + d != 1`` and a rescale by the total is seen; a single chain is seen).
    ``specular(0.3) + white(0.5)`` built in the other order gives the same ``k`` bit for bit.

    Red input: ``ImportError`` at the base tip. Mutations: M-rescale, M-first-chain-only.
    """
    W = H.bc().WallCondition
    keff_tol = 1e-10
    roles, _ = H.roles_and_conditions("slab", None)

    def k_of(wall: Any) -> float:
        return H.solve(H.chart_mesh("slab", **roles), {0: wall}).keff

    mixed, target = k_of(W.from_amplitudes(0.3, 0.5)), k_of(W.white(0.8))
    assert abs(mixed - target) <= 10 * keff_tol, f"mixed {mixed!r} vs white(0.8) {target!r}"
    for single in (0.3, 0.5):
        assert abs(k_of(W.white(single)) - target) > 1e-3, f"white({single}) is not distinguishable from white(0.8): the row cannot see a dropped chain"
    assert k_of(W.specular(0.3) + W.white(0.5)) == k_of(W.white(0.5) + W.specular(0.3))
