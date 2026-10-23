r"""The colouring symmetry of a layered geometry and the one coordinate-agreement convention
(wave 2 of the layer DAG, WF2.2 unit L1; ``.claude/plans/wave2_boundary_carve.md`` sections 2.2, 5.1;
``scratch/layer_dag/wave2_spec/measure_float_rule.md``).

Claims, with layer and kind:

* gate 7, ``test_colouring_symmetry_float_rule``: REFERENCE (the measured battery of ``measure_float_rule.md``:
  13 of 13 symmetric slabs accepted, 6 of 6 asymmetric refused, by ``64 * ulp(max(|b_0|, |b_n|))``); the activation
  is the interior breakpoint, so every row carries the palindromic id triple ``(0, 1, 0)`` (with one id the run
  merge deletes the interior breakpoints and the float rule has nothing to read);
* gate 8, ``test_colouring_symmetry_reads_runs_and_roles``: THEOREM (the group is a function of the colouring, not
  of the cell partition: refinement independence) and exact integer comparison of material ids and tags;
* gate 9, ``test_one_agreement_convention``: a STRUCTURE gate (X4: one definition per quantity), plus the
  convention's own value rows.

Like ``test_deck.py`` this module reads the production modules at run time (``colouring``, ``chart``), so the red
of each row is its own.

The calling convention the rows use, and the one question it leaves to the build (it is in this unit's NEEDS):
``colouring_symmetry(coord, ends, breaks, ids, role_low, role_high)`` with ``breaks`` the full ``b_0 .. b_n`` (the
ends included) and ``ends == (b_0, b_n)``.

Mutations (each on a COPY of ``colouring.py`` or ``chart.py``):

* gate 7: three mutants, each caught by a named row: the exact comparison (refused: the ``(0, .1, .2, .3)``
  rows), ``k * eps * (b - a)`` (refused: ``offset 1e6``), an absolute ``1e-9`` band (accepted: the
  ``1e-8``-scaled asymmetric rows); and a constant of 2 or 8 (refused: ``64 ulp - 4``) or of 1e6 (accepted:
  ``64 ulp + 6``);
* gate 8: compare material ids with ``np.isclose`` (the ``100000, 0, 100001`` row); do not merge adjacent equal ids
  (the refinement rows);
* gate 9: keep ``_LEVEL_AGREEMENT`` beside the new name, spell ``64 * np.spacing`` inside ``colouring.py``, make
  ``coordinates_agree`` scalar-only, give ``colouring`` a second definition.
"""

from __future__ import annotations

import ast
import importlib
import re
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from orpheus.geometry import CoordSystem
from orpheus.values.content import content_digest
from tests.gates._content_identity_helpers import require

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/geometry/test_colouring_symmetry.py"
_ROOT = Path(__file__).resolve().parents[3]
_CHORD_LEVEL_REFUSAL = "tests/gates/geometry/test_chord_level.py::test_a_level_that_disagrees_with_the_impact_parameter_is_refused"


class _Lazy:
    def __init__(self, name: str) -> None:
        self._name = name

    def __getattr__(self, attr: str) -> Any:
        return getattr(importlib.import_module(self._name), attr)


colouring: Any = _Lazy("orpheus.geometry.colouring")
chart: Any = _Lazy("orpheus.geometry.chart")
deck: Any = _Lazy("orpheus.geometry.deck")

_CART = CoordSystem.CARTESIAN
_PALINDROME = (0, 1, 0)


def _group(breaks, ids, low: Any = 1, high: Any = 1, coord: CoordSystem = _CART):
    """``colouring_symmetry`` on the RAW partition ``breaks`` / ``ids`` (the ends are the first and last break)."""
    return colouring.colouring_symmetry(coord, (breaks[0], breaks[-1]), tuple(breaks), tuple(ids), low, high)


def _mirrors(breaks, ids, low: Any = 1, high: Any = 1) -> bool:
    """``_group(...)`` is the mirror of ITS OWN partition's ends."""
    return _is_mirror(_group(breaks, ids, low, high), breaks)


def _is_mirror(group, breaks) -> bool:
    """The group is ``{e, x -> b_0 + b_n - x}``: identity first, one reflection swapping the two end faces, AT ITS
    POSITION: the motion sends ``b_0`` to ``b_n`` and ``b_n`` to ``b_0`` (``breaks`` are the partition's, ends included)."""
    els = group.elements
    if len(els) != 2 or not els[0].is_identity or els[1].is_identity:
        return False
    if dict(els[1].faces) != {"xmin": "xmax", "xmax": "xmin"} or not np.array_equal(els[1].motion.linear, [[-1.0]]):
        return False
    ends = np.array([[breaks[0]], [breaks[-1]]])
    return bool(np.allclose(els[1].motion.on_points(ends), ends[::-1], rtol=0.0, atol=1e-12 * max(abs(breaks[0]), abs(breaks[-1]))))


def _is_trivial(group) -> bool:
    return len(group.elements) == 1 and group.elements[0].is_identity


# ── gate 7: the float rule ───────────────────────────────────────────────────

_THIRD = (0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0)
_SYMMETRIC = {
    "(0,.1,.2,.3)": (0.0, 0.1, 0.2, 0.3),
    "(0,1,2,3)": (0.0, 1.0, 2.0, 3.0),
    "(0,1/3,2/3,1)": _THIRD,
    "(0,.1,.2,.3)*1e-8": tuple(b * 1e-8 for b in (0.0, 0.1, 0.2, 0.3)),
    "(0,1,2,3)*1e-8": tuple(b * 1e-8 for b in (0.0, 1.0, 2.0, 3.0)),
    "(0,1/3,2/3,1)*1e-8": tuple(b * 1e-8 for b in _THIRD),
    "offset 1e6": (1e6, 1e6 + 0.1, 1e6 + 0.2, 1e6 + 0.3),
    "negative (-.3,-.2,-.1,0)": (-0.3, -0.2, -0.1, 0.0),
    "unequal widths (0,.1,.3,.4)": (0.0, 0.1, 0.3, 0.4),
}
_ASYMMETRIC = {
    "(0,.1,.2,.30000001)": (0.0, 0.1, 0.2, 0.30000001),
    "(0,1,2.0000001,3)": (0.0, 1.0, 2.0000001, 3.0),
    "(0,.1,.2,.30000001)*1e-8": tuple(b * 1e-8 for b in (0.0, 0.1, 0.2, 0.30000001)),
    "(0,1,2.0000001,3)*1e-8": tuple(b * 1e-8 for b in (0.0, 1.0, 2.0000001, 3.0)),
    "(0,.1,.25,.3)": (0.0, 0.1, 0.25, 0.3),
}


@pytest.mark.rests_on(_CHORD_LEVEL_REFUSAL)
@pytest.mark.parametrize("name", list(_SYMMETRIC), ids=list(_SYMMETRIC))
def test_colouring_symmetry_float_rule_accepts_the_symmetric(name: str) -> None:
    """A symmetric slab is the mirror, however its breakpoints round (13 of 13 in the measurement).

    Reddens on: the exact comparison (every row but ``(0,1,2,3)``, the only row that does not activate the rule:
    ``(0,.1,.2,.3)``, ``(0,1/3,2/3,1)``, the three scalings, the negative slab and ``unequal widths`` all round
    inexactly), the ``k * eps * (b - a)`` mutant (``offset 1e6``), and a mirror whose position is not
    ``b_0 + b_n`` (``offset=0.0`` or ``offset=total``: ``_is_mirror`` reads the motion at the ends).
    """
    b = _SYMMETRIC[name]
    require(_mirrors(b, _PALINDROME), f"{name}: a symmetric slab {b} was not read as the mirror")


@pytest.mark.rests_on(_CHORD_LEVEL_REFUSAL)
@pytest.mark.parametrize("name", list(_ASYMMETRIC), ids=list(_ASYMMETRIC))
def test_colouring_symmetry_float_rule_refuses_the_asymmetric(name: str) -> None:
    """An asymmetric slab is the identity alone (6 of 6 in the measurement; the smallest real asymmetry is 4.5e7 ulp).

    Reddens on: the absolute ``1e-9`` band (the ``*1e-8`` rows), a tolerance wider than 64 ulp.
    """
    b = _ASYMMETRIC[name]
    require(_is_trivial(_group(b, _PALINDROME)), f"{name}: the asymmetric slab {b} was read as symmetric")


@pytest.mark.rests_on(_HERE + "::test_colouring_symmetry_float_rule_accepts_the_symmetric", _HERE + "::test_one_agreement_convention")
@pytest.mark.parametrize("ulps,accepted", [(60, True), (70, False)], ids=["64 ulp - 4", "64 ulp + 6"])
def test_colouring_symmetry_float_rule_constant_is_64_ulp_at_the_end_magnitude(ulps: int, accepted: bool) -> None:
    """The constant and the scale are pinned from both sides: the slab ``(0, 1, 2 + d, 3)`` with ``d`` a multiple of
    ``spacing(3)`` has mirror error ``d`` at breakpoints 1 and 2, and is the mirror iff ``d <= 64 * spacing(3)``
    (the scale ``max(|b_0|, |b_n|) = 3``).

    Reddens on: any constant below 60 or above 70 (2, 8, 1e6), and a scale taken as the interval's own width.
    """
    d = ulps * np.spacing(3.0)
    b = (0.0, 1.0, 2.0 + float(d), 3.0)
    require(2.0 + float(d) - 2.0 == d, "the perturbation is exactly representable at 2")
    require(_mirrors(b, _PALINDROME) == accepted, f"d = {ulps} ulp(3): expected accepted = {accepted}")


# ── gate 8: runs and roles ───────────────────────────────────────────────────

_B = (0.0, 1.0, 2.0, 3.0)


@pytest.mark.rests_on(_HERE + "::test_colouring_symmetry_float_rule_accepts_the_symmetric", _HERE + "::test_colouring_symmetry_float_rule_refuses_the_asymmetric")
def test_colouring_symmetry_reads_runs_and_roles() -> None:
    """Material runs compare EXACTLY (integer ids), and the end roles must be the same tag or the same deck type.

    Reddens on: a copy comparing material ids with a tolerance (``100000, 0, 100001``); a copy that treats two
    different tags as equal; a copy that lets a deck stand for a tag.
    """
    mirror, translation = deck.MirrorDeck(), deck.TranslationDeck()
    require(_mirrors(_B, (0, 1, 0)), "ids (0,1,0) are symmetric")
    require(_is_trivial(_group(_B, (0, 1, 2))), "ids (0,1,2) are not")
    require(_is_trivial(_group(_B, (100000, 0, 100001))), "ids (100000, 0, 100001): a tolerance on integers reads them equal")
    require(_mirrors(_B, (100000, 0, 100000)), "ids (100000, 0, 100000) are symmetric")

    for low, high, symmetric, why in (
        (1, 1, True, "the same tag"),
        (mirror, mirror, True, "two mirror decks"),
        (translation, translation, True, "two translation decks"),
        (0, 0, True, "the same tag 0 (a falsy tag is a tag)"),
        (1, 2, False, "two different tags"),
        (1, mirror, False, "a tag against a deck"),
        (translation, 0, False, "a deck against a tag"),
        (mirror, translation, False, "a mirror against a translation"),
    ):
        got = _group(_B, (0, 1, 0), low, high)
        require(_is_mirror(got, _B) if symmetric else _is_trivial(got), f"{why}: expected symmetric = {symmetric}")

    # radial charts: G_c absorbs the rotation and W is trivial, so the identity only, whatever the breaks
    for coord in (CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL):
        got = _group((0.5, 1.0, 1.5), (3, 3), 0, 0, coord)
        require(_is_trivial(got), f"{coord}: a radial chart returned {len(got.elements)} elements")

    # the pessimistic typed answer where no group machinery exists is a value
    u = colouring.UncomputedSymmetries("Mesh2D")
    require(u.carrier == "Mesh2D" and u == colouring.UncomputedSymmetries("Mesh2D"), "UncomputedSymmetries is a value")
    require(u != colouring.UncomputedSymmetries("Mesh3D"), "the carrier is part of the value")
    require(content_digest(u) != content_digest(colouring.UncomputedSymmetries("Mesh3D")), "and of its digest")
    require(not isinstance(u, colouring.ColouringSymmetry), "an uncomputed answer is not a computed group")

    # identity first, always
    require(_group(_B, (0, 1, 2)).elements[0].is_identity and _group(_B, (0, 1, 0)).elements[0].is_identity, "identity first")


@pytest.mark.rests_on(_HERE + "::test_colouring_symmetry_reads_runs_and_roles")
def test_colouring_symmetry_is_independent_of_the_cell_partition() -> None:
    """The refinement control (C8): a palindromic slab ``(0, 1, 2)`` of one material, given as two intervals, as a
    non-uniform cell partition with repeated ids and as a refined one, is the SAME group, because the function
    merges equal adjacent ids into runs itself. A mesh-derived group (cell edges as run breaks) would read the
    non-uniform partition as asymmetric: the attack of the thirteenth exchange (math 3).

    Reddens on: a copy that does not merge adjacent equal ids.
    """
    forms = {
        "two intervals": ((0.0, 1.0, 2.0), (0, 0)),
        "non-uniform cells": ((0.0, 0.25, 1.0, 2.0), (0, 0, 0)),
        "refined": ((0.0, 0.5, 1.0, 1.5, 2.0), (0, 0, 0, 0)),
        "one interval": ((0.0, 2.0), (0,)),
    }
    groups = {name: _group(b, ids) for name, (b, ids) in forms.items()}
    for name, g in groups.items():
        require(_is_mirror(g, forms[name][0]), f"{name}: a one-material slab is the mirror")
        require(g == groups["two intervals"], f"{name}: the group depends on the cell partition")

    runs = {
        "three runs": ((0.0, 1.0, 2.0, 3.0), (0, 1, 0)),
        "outer runs refined, uniformly": ((0.0, 0.5, 1.0, 2.0, 2.5, 3.0), (0, 0, 1, 0, 0)),
        "outer runs refined, asymmetrically": ((0.0, 0.25, 1.0, 2.0, 3.0), (0, 0, 1, 0)),
        "the middle run split": ((0.0, 1.0, 1.4, 2.0, 3.0), (0, 1, 1, 0)),
    }
    base = _group(*runs["three runs"])
    require(_is_mirror(base, runs["three runs"][0]), "the palindromic runs are the mirror")
    for name, (b, ids) in runs.items():
        require(_group(b, ids) == base, f"{name}: the group depends on the cell partition")
    # and a refinement that IS a different colouring is a different group
    require(_is_trivial(_group((0.0, 1.0, 1.5, 3.0), (0, 1, 0))), "moving a run break breaks the symmetry")


# ── gate 9: ONE agreement convention ─────────────────────────────────────────


def _geometry_sources() -> dict[Path, str]:
    paths = sorted((_ROOT / "orpheus" / "geometry").rglob("*.py"))
    require(len(paths) > 20 and (_ROOT / "orpheus/geometry/chart.py") in paths, f"the scan read {len(paths)} files")
    return {p: p.read_text(encoding="utf-8") for p in paths}


def _assignments_of(tree: ast.AST, name: str) -> list[ast.AST]:
    return [
        node
        for node in ast.walk(tree)
        if (isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name for t in node.targets))
        or (isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and node.target.id == name)
    ]


def _calls(node: ast.AST, name: str) -> list[ast.Call]:
    return [
        c
        for c in ast.walk(node)
        if isinstance(c, ast.Call) and ((isinstance(c.func, ast.Name) and c.func.id == name) or (isinstance(c.func, ast.Attribute) and c.func.attr == name))
    ]


def _spells_the_convention(node: ast.AST) -> list[str]:
    """The ways a body spells the convention itself: a ``spacing`` call, or the literal ``64``."""
    found = [f"spacing@{c.lineno}" for c in _calls(node, "spacing")]
    found += [f"64@{n.lineno}" for n in ast.walk(node) if isinstance(n, ast.Constant) and n.value in (64, 64.0) and not isinstance(n.value, bool)]
    return found


@pytest.mark.rests_on(_CHORD_LEVEL_REFUSAL)
def test_one_agreement_convention() -> None:
    """Two coordinates of one geometry AGREE when they differ by at most ``COORDINATE_AGREEMENT_ULPS`` (64) ulp at a
    stated scale: ONE constant, ONE function, read by the chord check and by the colouring rule.

    AST over every ``.py`` under ``orpheus/geometry/`` (exclusions: none; the definition site is included and is
    the point). Population: exactly one assignment of the constant (to ``64.0``), exactly one definition of
    ``coordinates_agree``, no ``_LEVEL_AGREEMENT`` anywhere (code, docstring or comment), the chord check
    (``RadialImage.at_level``) and ``colouring.py`` each CALL it and spell neither ``spacing`` nor ``64``.

    Reddens on: the base tip itself (``git grep -n _LEVEL_AGREEMENT`` reads 3 sites); a copy keeping the old name
    beside the new one; a copy of ``colouring.py`` spelling ``64 * np.spacing(...)``; a second definition.
    """
    sources = _geometry_sources()
    names = {p: ast.parse(s) for p, s in sources.items()}

    old = {str(p.relative_to(_ROOT)): [m.start() for m in re.finditer(r"_LEVEL_AGREEMENT", s)] for p, s in sources.items()}
    require(not any(old.values()), f"_LEVEL_AGREEMENT survives: { {k: len(v) for k, v in old.items() if v} }")

    constant = [(p, a) for p, t in names.items() for a in _assignments_of(t, "COORDINATE_AGREEMENT_ULPS")]
    require(len(constant) == 1, f"{len(constant)} assignments of COORDINATE_AGREEMENT_ULPS, expected 1")
    path, node = constant[0]
    require(path.name == "chart.py", f"the constant is defined in {path.name}, not chart.py")
    value = node.value if isinstance(node, ast.Assign | ast.AnnAssign) else None
    require(isinstance(value, ast.Constant) and value.value == 64.0 and isinstance(value.value, float), "the constant is 64.0")

    definitions = [(p, n) for p, t in names.items() for n in ast.walk(t) if isinstance(n, ast.FunctionDef) and n.name == "coordinates_agree"]
    require(len(definitions) == 1 and definitions[0][0].name == "chart.py", f"coordinates_agree defined at {[str(p) for p, _ in definitions]}")

    chart_tree = names[_ROOT / "orpheus/geometry/chart.py"]
    at_level = [n for n in ast.walk(chart_tree) if isinstance(n, ast.FunctionDef) and n.name == "at_level"]
    require(len(at_level) == 1, "RadialImage.at_level exists once")
    require(len(_calls(at_level[0], "coordinates_agree")) == 1, "the chord check calls coordinates_agree once")
    require(not _spells_the_convention(at_level[0]), f"the chord check spells the convention itself: {_spells_the_convention(at_level[0])}")

    colouring_tree = names[_ROOT / "orpheus/geometry/colouring.py"]
    require(len(_calls(colouring_tree, "coordinates_agree")) >= 1, "colouring.py never calls coordinates_agree")
    require(not _spells_the_convention(colouring_tree), f"colouring.py spells the convention itself: {_spells_the_convention(colouring_tree)}")

    require(colouring.coordinates_agree is chart.coordinates_agree, "colouring holds a second coordinates_agree")
    require(chart.COORDINATE_AGREEMENT_ULPS == 64.0, "the constant")
    require(not hasattr(colouring, "COORDINATE_AGREEMENT_ULPS"), "colouring holds a second address of the constant")
    for name in ("COORDINATE_AGREEMENT_ULPS", "coordinates_agree"):
        require(name in chart.__all__ and name not in colouring.__all__, f"{name}: one public address, chart's")


@pytest.mark.rests_on(_HERE + "::test_one_agreement_convention")
def test_the_agreement_convention_is_elementwise_and_scale_dependent() -> None:
    """``coordinates_agree(a, b, *, scale) := |a - b| <= 64 * spacing(scale)``, elementwise (the chord check passes an
    array ``scale``), on a scalar row and an array row.

    Reddens on: a copy that is scalar-only (the array row), a copy that ignores ``scale`` (the two scale rows
    disagree), ``<`` for ``<=`` (the edge row), a copy dividing by the scale (relative, not ulp).
    """
    agree = chart.coordinates_agree
    eps1 = np.spacing(1.0)
    require(bool(agree(1.0, 1.0 + 64 * eps1, scale=1.0)), "the edge: exactly 64 ulp agrees")
    require(not bool(agree(1.0, 1.0 + 65 * eps1, scale=1.0)), "65 ulp does not")
    require(not bool(agree(0.0, 1e-10, scale=1.0)), "1e-10 at scale 1 does not agree")
    require(bool(agree(0.0, 1e-10, scale=1e6)), "1e-10 at scale 1e6 agrees (spacing(1e6) = 1.2e-10)")
    require(bool(agree(2.5, 2.5, scale=0.0)), "equal values agree at any scale")

    a = np.array([1.0, 1e6, 3.0, 5.0])
    b = a + np.array([10 * np.spacing(1.0), 10 * np.spacing(1e6), 100 * np.spacing(3.0), 0.0])
    got = np.asarray(agree(a, b, scale=a))
    require(got.shape == (4,) and got.dtype == bool, f"the array row returned {got.shape} {got.dtype}")
    require(got.tolist() == [True, True, False, True], f"the array row read {got.tolist()}")
    scales = np.array([1.0, 1e6, 3.0, 5.0])
    require(np.asarray(agree(a, b, scale=scales)).tolist() == got.tolist(), "an array scale is read elementwise")
    require(np.asarray(agree(a, b, scale=1e6)).tolist() == [True, True, True, True], "a scalar scale broadcasts")


# ── the partition is parsed, the value's invariants hold ─────────────────────


@pytest.mark.rests_on(_HERE + "::test_colouring_symmetry_reads_runs_and_roles")
def test_colouring_symmetry_parses_the_partition_before_the_chart() -> None:
    """The raw partition is refused when it is not a partition, on every chart, before the chart is read.

    One negative leg per rule: the length (one more break than ids), no interval, unsorted breaks, a zero-width
    interval, ends that disagree with the first and last break; on a Cartesian AND a radial chart (the radial early
    return used to skip the check).

    Reddens on: a copy without the length check, a copy parsing after the radial return, a copy admitting
    non-increasing breaks, a copy that trusts ``ends`` over the breaks.
    """
    refused = {
        "too few breaks": ((0.0, 1.0), (0, 1, 0), (0.0, 1.0)),
        "too many breaks": ((0.0, 1.0, 2.0, 3.0, 4.0), (0, 1, 0), (0.0, 4.0)),
        "no interval": ((), (), (0.0, 0.0)),
        "one break, no ids": ((1.0,), (), (1.0, 1.0)),
        "unsorted": ((0.0, 2.0, 1.0, 3.0), (0, 1, 0), (0.0, 3.0)),
        "zero-width interval": ((0.0, 1.0, 1.0, 3.0), (0, 1, 0), (0.0, 3.0)),
        "ends disagree with the breaks": ((0.0, 1.0, 2.0, 3.0), (0, 1, 0), (0.0, 5.0)),
    }
    for coord in (_CART, CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL):
        for name, (breaks, ids, ends) in refused.items():
            with pytest.raises(ValueError):
                colouring.colouring_symmetry(coord, ends, breaks, ids, 1, 1)
                raise AssertionError(f"{coord} / {name}: the partition was accepted")  # noqa: B904 (a red, not a ValueError)
    require(_mirrors((0.0, 1.0, 2.0, 3.0), (0, 1, 0)), "the positive control: a well-formed partition is read")


@pytest.mark.rests_on(_HERE + "::test_colouring_symmetry_reads_runs_and_roles")
def test_the_colouring_values_state_their_invariants() -> None:
    """``FaceIsometry`` is a permutation of its own faces and ``is_identity`` reads BOTH fields; a
    ``ColouringSymmetry`` is non-empty with the identity first.

    Reddens on: a ``FaceIsometry.is_identity`` that reads the motion only (the swapped-faces row), a copy without
    either ``__post_init__``.
    """
    from orpheus.geometry.transformation import RigidMotion
    from orpheus.values.content import FrozenMapping

    swapped = FrozenMapping({"xmin": "xmax", "xmax": "xmin"})
    fixed = FrozenMapping({"xmin": "xmin", "xmax": "xmax"})
    require(not colouring.FaceIsometry(RigidMotion.identity(1), swapped).is_identity, "an identity motion with swapped faces is not the identity")
    require(colouring.FaceIsometry(RigidMotion.identity(1), fixed).is_identity, "the identity motion with fixed faces is")
    with pytest.raises(ValueError):
        colouring.FaceIsometry(RigidMotion.identity(1), FrozenMapping({"xmin": "xmax", "xmax": "xmax"}))
    with pytest.raises(ValueError):
        colouring.ColouringSymmetry(())
    with pytest.raises(ValueError):
        colouring.ColouringSymmetry((colouring.FaceIsometry(RigidMotion.identity(1), swapped),))
    require(not bool(chart.coordinates_agree(-1.0, -1.0 - 1e-3, scale=-1.0)), "and still refuses a real difference")
    require(bool(chart.coordinates_agree(-1.0, -1.0, scale=-1.0)), "a negative scale is its magnitude (np.spacing of a negative is negative)")
