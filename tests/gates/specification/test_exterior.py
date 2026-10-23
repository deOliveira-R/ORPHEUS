r"""The reference's wall value ``WallAmplitudes(specular, diffuse)`` (wave 2, unit L5).

Specified by the test-architect (2026-10-23, ``.claude/plans/wave2_boundary_carve.md``
5.5, decisions D2 and D3). The reference specification states a wall by the
cone coordinates ``(s, d)`` of its linear part, as PLAIN DATA: two reals in
[0, 1] whose sum is at most 1 and which are not both positive. A value
outside the bounds is a parse refusal (``ValueError``; ``TypeError`` for a
non-real, the ``values.scalars`` convention); a combined specular-and-diffuse
wall is a scope boundary that issue #650 lifts (``NotImplementedError``
naming #650). The bounds have ONE owner: the message fragment
``_FRAGMENT`` is defined once (law L-G); until slice S5 re-points
``derivations/continuous/characteristic/walls.py::Wall`` onto this class, the
old definition there is the population S5 retires (the strict ``xfail`` row
flips at S5).
"""

from __future__ import annotations

import ast
import dataclasses
from pathlib import Path

import pytest

from tests.gates._content_identity_helpers import require

pytestmark = pytest.mark.foundation

_ROOT = Path(__file__).resolve().parents[3]
_EXTERIOR = _ROOT / "orpheus" / "specification" / "exterior.py"
_FRAGMENT = "a wall is specular or diffuse"
_BOUNDS = "not a physical wall"
_FRAGMENTS = (_FRAGMENT, _BOUNDS)
_FORBIDDEN = ("orpheus.transport", "orpheus.derivations", "orpheus.geometry.boundary")


def _wall_amplitudes():
    from orpheus.specification.exterior import WallAmplitudes

    return WallAmplitudes


@pytest.mark.parametrize("specular,diffuse", [(0, 0), (0.4, 0), (0, 0.4), (1, 0), (0, 1), (0.0, 1.0)])
def test_wall_amplitudes_a_stated_wall_constructs(specular, diffuse) -> None:
    """Vacuum, a partial specular wall, a partial diffuse wall, the mirror and the white
    wall construct. Red at base: the import error. Mutation (a copy):
    ``__post_init__`` refuses ``(1, 0)``."""
    wall = _wall_amplitudes()(specular, diffuse)
    require(isinstance(wall.specular, float) and isinstance(wall.diffuse, float), "the coordinates are canonical floats")
    require((wall.specular, wall.diffuse) == (float(specular), float(diffuse)), f"read back {(wall.specular, wall.diffuse)}")


@pytest.mark.parametrize(
    "specular,diffuse",
    [(0.6, 0.6), (1.2, 0), (0, 1.2), (-0.1, 0), (0, -0.1), (1.0, 1.0), (float("nan"), 0), (0, float("nan"))],
)
def test_wall_amplitudes_a_value_outside_the_bounds_is_a_parse_refusal(specular, diffuse) -> None:
    """Outside [0, 1], a sum above 1 (the bounds are checked BEFORE the combined-wall
    scope boundary, so ``(0.6, 0.6)`` and ``(1, 1)`` are ``ValueError``), and NaN.
    Mutation: a copy that drops the sum bound reddens the ``(0.6, 0.6)`` row (it would
    raise ``NotImplementedError``); a copy that drops the lower bound reddens ``-0.1``."""
    cls = _wall_amplitudes()
    # NaN is refused upstream by ``parse_real`` (its own message); the bounds rows carry the pinned fragment.
    match = None if "nan" in (str(specular), str(diffuse)) else _BOUNDS
    with pytest.raises(ValueError, match=match):
        cls(specular, diffuse)


@pytest.mark.parametrize("bad", [True, False, "0.5", None, 0.5 + 0.1j])
def test_wall_amplitudes_a_non_real_is_refused_as_a_type(bad) -> None:
    """``bool`` is not a number here (``values.scalars.parse_real``): ``True`` must not
    be read as the mirror. Both coordinates are checked."""
    cls = _wall_amplitudes()
    with pytest.raises(TypeError):
        cls(bad, 0.0)
    with pytest.raises(TypeError):
        cls(0.0, bad)


@pytest.mark.parametrize("specular,diffuse", [(0.5, 0.5), (0.25, 0.75), (1e-300, 1e-300), (0.5, 1e-9)])
def test_wall_amplitudes_a_combined_wall_is_a_scope_boundary_naming_650(specular, diffuse) -> None:
    """Both coordinates positive (sum within the bound) is a wall the closed vocabulary
    does not state yet: ``NotImplementedError`` naming #650, NOT ``ValueError`` (the
    refusal is lifted by widening, not by a better input). ``(1e-300, 1e-300)`` and
    ``(0.5, 1e-9)`` make the 'both positive' predicate bite at the extremes (a
    ``> 1e-9`` threshold passes the first). Mutation: a copy accepting ``(0.5, 0.5)``
    reddens all rows."""
    cls = _wall_amplitudes()
    with pytest.raises(NotImplementedError, match="#650"):
        cls(specular, diffuse)
    with pytest.raises(NotImplementedError, match=_FRAGMENT):
        cls(specular, diffuse)


def test_wall_amplitudes_is_a_frozen_content_identity_with_canonical_parts() -> None:
    """``ContentIdentity``, ``frozen=True``, ``eq=False`` (the S5.9 population gate reds on
    ``eq=True``), the parts exactly ``(specular, diffuse)``, and ``-0.0`` folded."""
    from orpheus.values.content import ContentIdentity

    cls = _wall_amplitudes()
    require(issubclass(cls, ContentIdentity), "WallAmplitudes is a ContentIdentity")
    params = cls.__dataclass_params__  # type: ignore[attr-defined]
    require(params.frozen and not params.eq, f"frozen={params.frozen}, eq={params.eq}")
    require([f.name for f in dataclasses.fields(cls)] == ["specular", "diffuse"], "the fields are (specular, diffuse)")
    wall = cls(-0.0, 0.5)
    require(str(wall.specular) == "0.0", f"-0.0 is folded, got {wall.specular!r}")
    with pytest.raises(dataclasses.FrozenInstanceError):
        wall.specular = 0.1  # pyright: ignore[reportAttributeAccessIssue]


def _defines_fragment(source: str, fragment: str = _FRAGMENT) -> int:
    """The number of string literals in ``source`` containing ``fragment`` (AST, never a line grep)."""
    return sum(
        isinstance(n, ast.Constant) and isinstance(n.value, str) and fragment in n.value
        for n in ast.walk(ast.parse(source))
    )


def test_the_wall_bounds_message_is_owned_once_in_the_new_class() -> None:
    """Law L-G at L5: ``exterior.py`` defines the combined-wall message fragment in exactly
    one literal. Positive control: the counter finds the fragment in a planted source and
    counts two literals as two. Red at base: ``exterior.py`` does not exist."""
    require(_defines_fragment(f"x = '{_FRAGMENT} (planted)'") == 1, "control: one planted literal counts 1")
    require(_defines_fragment(f"x = '{_FRAGMENT}'\ny = '{_FRAGMENT}'") == 2, "control: two literals count 2")
    for fragment in _FRAGMENTS:
        n = _defines_fragment(_EXTERIOR.read_text(), fragment)
        require(n == 1, f"{_EXTERIOR.name} defines {n} literals containing {fragment!r}, not 1")


@pytest.mark.xfail(
    strict=True, raises=AssertionError,
    reason="slice S5 re-points Wall onto WallAmplitudes; until then walls.py:111-125 still defines the fragment (law L-G)",
)
def test_the_wall_bounds_message_is_defined_in_one_place_across_orpheus() -> None:
    """Flips at S5: after the re-point exactly one file under ``orpheus/`` (``exterior.py``)
    defines the fragment. Strictness is asserted by introspection below."""
    for fragment in _FRAGMENTS:
        owners = sorted(
            str(p.relative_to(_ROOT))
            for p in (_ROOT / "orpheus").rglob("*.py")
            if _defines_fragment(p.read_text(), fragment) > 0
        )
        require(owners == ["orpheus/specification/exterior.py"], f"{fragment!r} is defined in {owners}")


def test_the_one_place_row_is_a_strict_xfail_with_a_reason() -> None:
    """A marker moved into ``pytest.param(marks=...)`` loses strictness: assert it on the row itself."""
    marks = [m for m in test_the_wall_bounds_message_is_defined_in_one_place_across_orpheus.pytestmark if m.name == "xfail"]
    require(len(marks) == 1 and marks[0].kwargs.get("strict") is True and bool(marks[0].kwargs.get("reason")), "strict xfail with reason")


def _forbidden_imports(source: str) -> list[str]:
    found: list[str] = []
    for node in ast.walk(ast.parse(source)):
        names: list[str] = []
        if isinstance(node, ast.Import):
            names = [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                names = ["<relative>"]
            else:
                base = node.module or ""
                names = [base] + [f"{base}.{a.name}" for a in node.names]
        found += [n for n in names if n == "<relative>" or any(n == f or n.startswith(f + ".") for f in _FORBIDDEN)]
    return found


def test_the_value_imports_nothing_across() -> None:
    """``exterior.py`` imports nothing from transport, derivations or geometry.boundary,
    and no relative import (which could reach them). Positive controls: the scanner finds
    each forbidden prefix in each import spelling, including ``from pkg import module``.
    Mutation: a planted ``from orpheus.transport import x`` on a copy."""
    for planted in (
        "import orpheus.transport.core",
        "from orpheus.derivations import continuous",
        "from orpheus.geometry import boundary",
        "from orpheus.geometry.boundary import _base",
        "from . import sibling",
    ):
        require(bool(_forbidden_imports(planted)), f"control: {planted!r} not found")
    require(_forbidden_imports("from orpheus.values.scalars import parse_real\nimport math") == [], "control: a clean source is clean")
    found = _forbidden_imports(_EXTERIOR.read_text())
    require(found == [], f"{_EXTERIOR.name} imports {found}")
