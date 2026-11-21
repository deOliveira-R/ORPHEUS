r"""G3 of the unit ``insulation`` (#641 item 1): the composite Gauss-Legendre rule CP draws over the impact parameter.

``orpheus.cp.solver._composite_gauss_legendre(partition, n_per_panel)`` (``partition.breakpoints`` the panel ends) is the DIRECT SUM of the panels' Gauss-Legendre measures, each
pushed forward along the inclusion of its panel into ``Interval(b_0, b_last)`` and returned as a ``DiscreteMeasure``. Its laws:

* (b) size: ``nodes.shape == weights.shape == ((len(b) - 1) * n,)`` and ``support == Interval(b_0, b_last)``;
* (c) order: the nodes are strictly increasing and panel ``k`` lies in ``[b_k, b_{k+1}]``;
* (d) exactness PER PANEL, the law the rule exists for: for a function that is a polynomial of degree at most ``2n - 1`` on each panel,
  ``weights @ f(nodes)`` equals the integral over ``[b_0, b_last]`` to rounding. The witness is the kinked ``|x - 0.4|^3``: a rule that is
  exact per panel integrates it exactly (the kink is a breakpoint); one global rule of the same size does not (error 5.2e-9);
* (e) agreement with the reference rule, ``orpheus.derivations.common.quadrature.composite_gauss_legendre`` (numpy's ``leggauss`` mapped
  by float arithmetic, an independent construction), at ``rtol = 1e-10``: not bit-equal (the reference's small end weights carry about
  1e-12 relative error; the production rule is correctly rounded);
* (f) total in ``n``: the odd ``n`` legs. A panel such as ``[0.2, 0.4]`` has an exact rational centre node on a float64 rounding tie at
  odd ``n``; it raised ``ArithmeticError`` before the unit ``quadrature-ties`` and builds after it.
* (a) the refusals of a malformed partition belong to ``ConcentricPartition`` (the geometry type's own law), not to the helper; the two
  refused legs pin the geometry law and are not a witness for the CP input boundary, which ``Mesh1D`` parses first.

First red at the base tip: every row that calls the helper reds with ``AttributeError`` (the helper does not exist; it is read at call time so
each row reports on its own). The laws are written ONCE, in ``_law_violations``, and the mutant rows below apply that checker to three
mutants built from numpy alone (no production helper), so each arm reddening is a fact of the checker and not of the landing:
M-drop-last-panel, M-wrong-n (63 points per panel where 64 were asked), M-one-global-rule (one ``LEGENDRE.on(b_0, b_last).gauss(n k)``:
reddens (d) by 5.2e-9 against 0.0, and (c)).
"""

from __future__ import annotations

import types

import numpy as np
import pytest

from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.chart import Chart
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.manifold import Interval
from orpheus.numerics.generating_measure import LEGENDRE

_EDGES = (0.0, 0.2, 0.4, 0.55, 0.63)
_KINK = 0.4


def _breakpoints() -> tuple[float, ...]:
    """The edges through the partition, as ``CPMesh`` hands them to the helper (the one parse)."""
    return ConcentricPartition(Chart(CoordSystem.CYLINDRICAL), _EDGES).breakpoints


def _helper(breakpoints: tuple[float, ...], n: int):
    from orpheus.cp import solver

    partition = ConcentricPartition(Chart(CoordSystem.CYLINDRICAL), _EDGES)
    assert partition.breakpoints == breakpoints
    return solver._composite_gauss_legendre(partition, n)


def _exact_kinked_integral() -> float:
    """``integral_0^0.63 |x - 0.4|^3 dx = 0.4^4 / 4 + 0.23^4 / 4``, from the closed form."""
    return _KINK**4 / 4 + (_EDGES[-1] - _KINK) ** 4 / 4


def _law_violations(rule, breakpoints: tuple[float, ...], n: int) -> list[str]:
    """Laws (b), (c), (d) of a rule given by ``.nodes``, ``.weights``, ``.support``; the empty list is a pass."""
    found: list[str] = []
    panels = len(breakpoints) - 1
    x, w = np.asarray(rule.nodes), np.asarray(rule.weights)
    if not (x.shape == w.shape == (panels * n,)):
        return [f"(b) size: nodes {x.shape}, weights {w.shape}, expected {(panels * n,)}"]
    if rule.support != Interval(breakpoints[0], breakpoints[-1]):
        found.append(f"(b) support {rule.support!r}")
    if not np.all(np.diff(x) > 0):
        found.append("(c) the nodes are not strictly increasing")
    for k in range(panels):
        block = x[k * n : (k + 1) * n]
        if not (breakpoints[k] <= block.min() and block.max() <= breakpoints[k + 1]):
            found.append(f"(c) panel {k} leaves [{breakpoints[k]}, {breakpoints[k + 1]}]")
    error = abs(float(w @ np.abs(x - _KINK) ** 3) - _exact_kinked_integral())
    if not error <= 1e-15:
        found.append(f"(d) per-panel exactness on |x - {_KINK}|^3: error {error!r}")
    return found


def _reference(breakpoints: tuple[float, ...], n: int):
    from orpheus.derivations.common.quadrature import composite_gauss_legendre

    return composite_gauss_legendre(list(breakpoints), n)


@pytest.mark.l1
@pytest.mark.parametrize("n", [8, 64, 3, 63], ids=["n8", "n64", "n3-odd", "n63-odd"])
def test_the_composite_rule_obeys_its_laws_and_agrees_with_the_reference(n: int) -> None:
    r"""Laws (b), (c), (d) and (e) on the four-panel edges, at two even and two odd ``n`` (law (f): the odd rows are the panels
    ``[0.2, 0.4]`` and ``[0.55, 0.63]``, which raised ``ArithmeticError`` on a rounding tie before ``quadrature-ties``).

    Red input at the base tip: ``AttributeError`` (no helper). The odd rows additionally read ``ArithmeticError`` on the tree before
    ``quadrature-ties``; that unit is integrated at this base, so they are green once the helper exists.
    Mutations: M-drop-last-panel, M-wrong-n, M-one-global-rule (see the mutant rows).
    """
    b = _breakpoints()
    rule = _helper(b, n)
    assert _law_violations(rule, b, n) == []
    ref = _reference(b, n)
    np.testing.assert_allclose(rule.nodes, ref.pts, rtol=1e-10, atol=0.0)
    np.testing.assert_allclose(rule.weights, ref.wts, rtol=1e-10, atol=0.0)


@pytest.mark.l1
def test_the_composite_rule_is_the_concatenation_of_its_panels_bit_for_bit() -> None:
    """The direct sum equals the hand concatenation of the panels' correctly rounded rules (``array_equal``: a pushforward along the
    identity inclusion and a sum reorder no addition). Red input at the base tip: ``AttributeError``. Mutation: M-drop-last-panel."""
    b = _breakpoints()
    rule = _helper(b, 8)
    panels = [LEGENDRE.on(lo, hi).gauss(8) for lo, hi in zip(b[:-1], b[1:])]
    assert np.array_equal(rule.nodes, np.concatenate([p.nodes for p in panels]))
    assert np.array_equal(rule.weights, np.concatenate([p.weights for p in panels]))


def _mutant(nodes: list[np.ndarray], weights: list[np.ndarray], support: Interval) -> types.SimpleNamespace:
    return types.SimpleNamespace(nodes=np.concatenate(nodes), weights=np.concatenate(weights), support=support)


def _mutants() -> dict[str, tuple[types.SimpleNamespace, int]]:
    b = _breakpoints()
    n = 8
    whole = Interval(b[0], b[-1])
    panels = [LEGENDRE.on(lo, hi).gauss(n) for lo, hi in zip(b[:-1], b[1:])]
    nodes, weights = [p.nodes for p in panels], [p.weights for p in panels]
    dropped = _mutant(nodes[:-1], weights[:-1], whole)
    wrong_n_panels = [LEGENDRE.on(lo, hi).gauss(n - 1) for lo, hi in zip(b[:-1], b[1:])]
    wrong_n = _mutant([p.nodes for p in wrong_n_panels], [p.weights for p in wrong_n_panels], whole)
    one = LEGENDRE.on(b[0], b[-1]).gauss(n * (len(b) - 1))
    global_rule = types.SimpleNamespace(nodes=one.nodes, weights=one.weights, support=whole)
    honest = _mutant(nodes, weights, whole)
    return {"honest": (honest, n), "drop-last-panel": (dropped, n), "wrong-n": (wrong_n, n), "one-global-rule": (global_rule, n)}


@pytest.mark.foundation
@pytest.mark.parametrize("arm", ["drop-last-panel", "wrong-n", "one-global-rule"])
def test_the_laws_redden_on_each_mutant_of_the_rule(arm: str) -> None:
    r"""Each mutant is built from numpy and ``LEGENDRE`` alone, never from the production helper. Positive control: the ``honest``
    rule (the panels concatenated) passes the same checker. ``one-global-rule`` is read by law (d), error 5.2e-9 against 0.0, and also
    by law (c) (its nodes cross the panel edges).
    """
    honest, n = _mutants()["honest"]
    assert _law_violations(honest, _breakpoints(), n) == []
    mutant, n = _mutants()[arm]
    violations = _law_violations(mutant, _breakpoints(), n)
    assert violations, f"the checker is blind to the mutant {arm}"
    if arm == "one-global-rule":
        assert any(v.startswith("(d)") for v in violations), violations


@pytest.mark.foundation
@pytest.mark.parametrize("edges", [(0.0,), (0.0, 0.5, 0.5)], ids=["one-breakpoint", "not-increasing"])
def test_a_malformed_partition_is_refused_by_the_geometry_type(edges: tuple[float, ...]) -> None:
    r"""Law (a): the helper makes no check of its own; the refusals are ``ConcentricPartition``'s (``ValueError``, fragments
    ``strictly increasing`` and ``at least two``). Mutation: a ``monkeypatch`` of ``ConcentricPartition.__post_init__`` that skips the
    check makes these rows read no error (the next test); the rows pin the geometry law and are not a witness for the CP input boundary."""
    with pytest.raises(ValueError, match="strictly increasing"):
        ConcentricPartition(Chart(CoordSystem.CYLINDRICAL), edges)


@pytest.mark.foundation
def test_the_partition_refusal_row_reddens_when_the_check_is_skipped(monkeypatch: pytest.MonkeyPatch) -> None:
    def skipped(self: ConcentricPartition) -> None:
        object.__setattr__(self, "breakpoints", tuple(self.breakpoints))

    monkeypatch.setattr(ConcentricPartition, "__post_init__", skipped)
    ConcentricPartition(Chart(CoordSystem.CYLINDRICAL), (0.0, 0.5, 0.5))  # accepted: the refusal row above would read red
