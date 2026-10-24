r"""The reference pairing loop and the plain-data edge it crosses (wave 2, S5; plan 2.9 and gate 8 of 6.5; decision C10).

Specified by the test-architect (2026-10-23, ``.claude/plans/wave2_boundary_carve.md`` 6.5). A reference states its exterior as
plain data (``WallAmplitudes(specular, diffuse)``, the cone coordinates of a wall's linear part); production states it as
``WallCondition`` values; the translation is ONE test-side loop, :func:`tests.gates._reference_pairing.conditions_of`, built on
the primitives-taking class method ``WallCondition.from_amplitudes(specular, diffuse)``. This module gates that loop and the
two places the edge is spelled in the references' own plain data:

* **The loop** (``test_the_reference_pairing_loop_*``): ``conditions_of({1: WallAmplitudes(0.3, 0), 2: WallAmplitudes(0, 0.7)})``
  is a ``BoundaryConditions`` of ``WallCondition.specular(0.3)`` and ``WallCondition.white(0.7)`` with the key set preserved
  exactly. Four legs: the positive (the stated pair), the key set (a tag in, the same tag out, the empty exterior in, the empty
  table out), the vertices of the cone (``(0, 0)`` is the vacuum condition, ``(1, 0)`` the mirror, ``(0, 1)`` the white wall),
  and the NEGATIVE (the swapped pair is not equal to it: an equality that cannot fail certifies nothing, ``vv-principles``
  anti-pattern 11). The loop is NOT a tautology (``instrument-doctrine`` X4): ``WallAmplitudes`` is read twice, by the reference
  (into ``Wall``) and by ``from_amplitudes`` (into a production condition), and the pairing evidence is the comparison of the
  reference SOLUTION with the production SOLUTION (``test_reference_carry.py`` for the reference side; the method pairing tests
  for production); no row here asserts that ``from_amplitudes(w)`` "realizes the same albedo as" ``Wall(w)``.
* **The field name** (``test_problem_spec_*``, an AST read of ``continuous_reference.py`` so a red here is the field and not an
  import chain): ``ProblemSpec`` has no field named ``boundary_conditions``. That is the production keyword's exact name, kept
  across a string-valued plain-data edge, and a reader meeting ``boundary_conditions`` on a reference could take it for the
  production keyword it is not (decision C10). It is ``boundary_kinds``.
* **The vocabulary** (``test_the_reference_boundary_kind_*``): the ``Literal`` alias is ``ReferenceBoundaryKind`` (the name
  ``BoundaryCondition`` is the production abstract base of ``orpheus.transport.boundary_conditions``, and two things under one
  name is exactly what the carve is retiring), its value set contains ``"zero_flux"`` (the latent type lie of
  ``cases/diffusion.py``, which the alias never listed), and every literal the derivations pass as a ``boundary_kinds``
  (or, at the base, ``boundary_conditions``) value is a member of it. The last row is a census (population stated in its own
  assertion message, its filter proven by a planted source), and it reads RED at the base for a second, separate reason: the
  derivations also pass ``"white_rank2"`` (``peierls_nystrom/cylinder.py``, ``sphere.py``), a value the alias never listed
  either; the build adds it to the alias or states the kind differently, and the row says which of the two it must be.

Declared red at the base tip (``de354496``) and the input that reddens each row (``instrument-doctrine`` X1): the field-name
row reads ``boundary_conditions`` in the ``ProblemSpec`` body; the alias rows find no ``ReferenceBoundaryKind`` assignment; the
census finds ``zero_flux`` and ``white_rank2`` outside the alias. The loop rows are GREEN at the base: ``from_amplitudes`` and
``BoundaryConditions`` exist since S2a, ``WallAmplitudes`` since L5, and ``conditions_of`` is created in this commit (test-side
infrastructure, plan 6.5), so they are the live control that the new edge is crossed correctly. Mutations the reviewers apply:
(1) in ``conditions_of`` swap ``w.specular`` and ``w.diffuse``: the positive and the vertices rows redden; (2) drop one tag
(``list(exterior)[1:]``): the key-set row reddens; (3) build every condition from ``WallCondition.vacuum()``: positive,
vertices and negative rows redden; (4) name the field ``boundary_kinds`` but keep ``boundary_conditions`` as a second field
(a leftover): the field-name row reddens, which is why it asserts the field's ABSENCE and the new field's PRESENCE.
"""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Any

import pytest

from tests.gates._content_identity_helpers import require

pytestmark = pytest.mark.foundation

_ROOT = Path(__file__).resolve().parents[3]
_CONTINUOUS_REFERENCE = _ROOT / "orpheus" / "derivations" / "common" / "continuous_reference.py"
_HERE = "tests/gates/derivations/test_reference_pairing_loop.py::"
_EXTERIOR = "tests/gates/specification/test_exterior.py::"


# ── the loop ─────────────────────────────────────────────────────────────


def _amplitudes(specular: float, diffuse: float):
    from orpheus.specification.exterior import WallAmplitudes

    return WallAmplitudes(specular, diffuse)


@pytest.mark.rests_on(_EXTERIOR + "test_wall_amplitudes_a_stated_wall_constructs")
def test_the_reference_pairing_loop_maps_each_tag_to_the_condition_its_amplitudes_state() -> None:
    """``{1: (0.3, 0), 2: (0, 0.7)}`` becomes the specular 0.3 and white 0.7 conditions, under the same tags.

    The expected value is spelled through the production facades (``specular``, ``white``), not through ``from_amplitudes``:
    the facade and the primitive are two spellings the production side keeps equal, and the loop is correct iff it lands on
    the one the facade names. First red (the mutation of the docstring, (1)): ``specular`` and ``diffuse`` swapped reads
    ``white(0.3)`` and ``specular(0.7)``."""
    from orpheus.transport.boundary_conditions import BoundaryConditions, WallCondition
    from tests.gates._reference_pairing import conditions_of

    exterior = {1: _amplitudes(0.3, 0.0), 2: _amplitudes(0.0, 0.7)}
    built = conditions_of(exterior)
    require(isinstance(built, BoundaryConditions), f"conditions_of returns a BoundaryConditions, got {type(built).__name__}")
    require(built[1] == WallCondition.specular(0.3), f"tag 1 reads {built[1]!r}, not specular(0.3)")
    require(built[2] == WallCondition.white(0.7), f"tag 2 reads {built[2]!r}, not white(0.7)")
    require(built == BoundaryConditions({1: WallCondition.specular(0.3), 2: WallCondition.white(0.7)}), "the table is the stated one")


@pytest.mark.rests_on(_HERE + "test_the_reference_pairing_loop_maps_each_tag_to_the_condition_its_amplitudes_state")
def test_the_reference_pairing_loop_preserves_the_key_set_exactly() -> None:
    """A tag in, the same tag out; the empty exterior (a deck-only geometry) gives the empty table; non-contiguous tags survive.

    The production exact cover refuses a condition for a tag no face carries and a face with no condition, so a loop that
    renumbered, dropped or invented a tag would move the refusal from the specification to the solve. First red: ``list(exterior)[1:]``
    (the key row), a ``range(len(exterior))`` renumbering (the non-contiguous tags ``{3, 7}`` come back ``{0, 1}``)."""
    from tests.gates._reference_pairing import conditions_of

    for tags in ((), (1, 2), (3, 7), (0,), (5, 0, 9)):
        exterior = {t: _amplitudes(0.0, 0.0) for t in tags}
        built = conditions_of(exterior)
        require(sorted(built.keys()) == sorted(tags), f"tags {tags} come back as {sorted(built.keys())}")


@pytest.mark.rests_on(_HERE + "test_the_reference_pairing_loop_maps_each_tag_to_the_condition_its_amplitudes_state")
def test_the_reference_pairing_loop_reads_the_vertices_of_the_cone_as_their_conditions() -> None:
    """``(0, 0)`` is the vacuum condition, ``(1, 0)`` the mirror (``specular(1)``), ``(0, 1)`` the white wall (``white(1)``).

    The three vertices of the triangle of sub-Markov walls, each against its facade; a loop that treated ``(1, 0)`` as
    ``vacuum`` (an amplitude test of ``> 1`` for ``>= 1``) reddens the mirror vertex; one that dropped a zero reddens vacuum."""
    from orpheus.transport.boundary_conditions import WallCondition
    from tests.gates._reference_pairing import conditions_of

    expected = {
        (0.0, 0.0): WallCondition.vacuum(),
        (1.0, 0.0): WallCondition.specular(1.0),
        (0.0, 1.0): WallCondition.white(1.0),
    }
    built = conditions_of({k: _amplitudes(*pair) for k, pair in enumerate(expected)})
    for k, (pair, condition) in enumerate(expected.items()):
        require(built[k] == condition, f"the vertex {pair} reads {built[k]!r}, not {condition!r}")


@pytest.mark.rests_on(_HERE + "test_the_reference_pairing_loop_maps_each_tag_to_the_condition_its_amplitudes_state")
def test_the_reference_pairing_loop_distinguishes_a_swapped_pair() -> None:
    """The NEGATIVE leg: ``(0.3, 0)`` and ``(0, 0.3)`` are different conditions, and a different amplitude is a different condition.

    Without it the positive row passes any ``conditions_of`` for which ``WallCondition`` equality is coarse (a ``__eq__`` that
    ignores the chain), which is a second reading of the same value (``vv-principles`` anti-pattern 19: cite the reading under
    the deliberately wrong structure)."""
    from tests.gates._reference_pairing import conditions_of

    specular = conditions_of({0: _amplitudes(0.3, 0.0)})
    diffuse = conditions_of({0: _amplitudes(0.0, 0.3)})
    other = conditions_of({0: _amplitudes(0.4, 0.0)})
    require(specular != diffuse, "a 0.3 specular wall and a 0.3 diffuse wall are one table")
    require(specular != other, "a 0.3 and a 0.4 specular wall are one table")
    require(specular == conditions_of({0: _amplitudes(0.3, 0.0)}), "the same stated wall is not equal to itself")


# ── the field name and the vocabulary, read from the source ──────────────


def _tree() -> ast.Module:
    return ast.parse(_CONTINUOUS_REFERENCE.read_text())


def _class_fields(tree: ast.Module, name: str) -> list[str]:
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == name:
            return [n.target.id for n in node.body if isinstance(n, ast.AnnAssign) and isinstance(n.target, ast.Name)]
    raise AssertionError(f"no class {name} in {_CONTINUOUS_REFERENCE.name}")


def _literal_members(tree: ast.Module, alias: str) -> list[str] | None:
    """The string members of ``alias = Literal[...]`` at module level, or None when it is not assigned."""
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == alias for t in node.targets):
            value = node.value
            if isinstance(value, ast.Subscript):
                inner = value.slice
                elements = inner.elts if isinstance(inner, ast.Tuple) else [inner]
                return [e.value for e in elements if isinstance(e, ast.Constant) and isinstance(e.value, str)]
    return None


def test_problem_spec_has_no_field_named_boundary_conditions_and_has_boundary_kinds() -> None:
    """The field is ``boundary_kinds``: its absence is checked beside the new field's presence (a leftover second field reds it).

    Positive control: the AST reader finds the known field ``geometry_type`` in the same class, and finds no field in an
    absent class name (the reader raises instead of reading an empty list). Red at base: the body declares
    ``boundary_conditions``."""
    tree = _tree()
    fields = _class_fields(tree, "ProblemSpec")
    require("geometry_type" in fields, f"control: the reader found no known field in {fields}")
    with pytest.raises(AssertionError, match="no class"):
        _class_fields(tree, "NoSuchClass")
    require("boundary_conditions" not in fields, f"ProblemSpec still has the field boundary_conditions: {fields}")
    require("boundary_kinds" in fields, f"ProblemSpec has no field boundary_kinds: {fields}")


@pytest.mark.rests_on(_HERE + "test_problem_spec_has_no_field_named_boundary_conditions_and_has_boundary_kinds")
def test_the_reference_boundary_kind_alias_is_renamed_and_lists_zero_flux() -> None:
    """``ReferenceBoundaryKind`` is assigned, ``"zero_flux"`` is one of its values, and ``BoundaryCondition`` is no longer assigned here.

    ``[M]`` the alias value set at the base is ``vacuum, reflective, white, marshak, periodic`` (5 members): ``"zero_flux"``
    is used by ``cases/diffusion.py`` (2 sites) and was never listed, a type lie. Red at base: no ``ReferenceBoundaryKind``."""
    tree = _tree()
    require(_literal_members(tree, "OperatorForm") is not None, "control: the reader finds the known alias OperatorForm")
    members = _literal_members(tree, "ReferenceBoundaryKind")
    require(members is not None, "continuous_reference.py does not assign ReferenceBoundaryKind")
    assert members is not None
    require("zero_flux" in members, f"ReferenceBoundaryKind lacks 'zero_flux': {members}")
    require(len(members) == len(set(members)), f"a value is listed twice in {members}")
    require(_literal_members(tree, "BoundaryCondition") is None, "the alias BoundaryCondition is still assigned in continuous_reference.py")


def _kind_values(source: str, keywords: tuple[str, ...]) -> tuple[list[str], int]:
    """The string literals that are DICT VALUES of a keyword argument named in ``keywords``, and the count of non-literal values."""
    literals: list[str] = []
    other = 0
    for node in ast.walk(ast.parse(source)):
        if not isinstance(node, ast.keyword) or node.arg not in keywords or not isinstance(node.value, ast.Dict):
            continue
        for value in node.value.values:
            if isinstance(value, ast.Constant) and isinstance(value.value, str):
                literals.append(value.value)
            else:
                other += 1
    return literals, other


@pytest.mark.rests_on(_HERE + "test_the_reference_boundary_kind_alias_is_renamed_and_lists_zero_flux")
def test_the_reference_boundary_kind_covers_every_literal_the_derivations_pass() -> None:
    """Every string literal passed as a boundary-kind dict value anywhere under ``orpheus/derivations`` is an alias member.

    The census reads the keyword ``boundary_kinds`` (and ``boundary_conditions``, the pre-rename spelling, so the row also
    covers the base). Its filter is proven first: a planted source with ``boundary_kinds={"a": "x"}`` yields ``["x"]`` and one
    with a non-literal value counts it in the second return. The population (files, literal values, non-literal values) is in
    the failure message (``instrument-doctrine`` X2). Red at base: ``zero_flux`` (2 sites, ``cases/diffusion.py``) and
    ``white_rank2`` (``peierls_nystrom/cylinder.py``, ``sphere.py``) are not members; the build either lists ``white_rank2`` or
    states that kind another way, and this row follows it."""
    kw = ("boundary_kinds", "boundary_conditions")
    planted = _kind_values('f(boundary_kinds={"a": "x", "b": value})', kw)
    require(planted == (["x"], 1), f"control: the filter reads {planted}")
    require(_kind_values('f(other={"a": "x"})', kw) == ([], 0), "control: an unrelated keyword is not read")
    # The pre-rename spelling is the fallback so that, at the base, this row reads the VALUES and not the missing name
    # (the alias row above owns the name): its own red is then ``zero_flux`` and ``white_rank2``.
    members = _literal_members(_tree(), "ReferenceBoundaryKind") or _literal_members(_tree(), "BoundaryCondition")
    require(members is not None, "continuous_reference.py assigns neither ReferenceBoundaryKind nor BoundaryCondition")
    assert members is not None
    seen: dict[str, set[str]] = {}
    files = 0
    non_literal = 0
    for path in sorted((_ROOT / "orpheus" / "derivations").rglob("*.py")):
        files += 1
        literals, other = _kind_values(path.read_text(), kw)
        non_literal += other
        for value in literals:
            seen.setdefault(value, set()).add(str(path.relative_to(_ROOT)))
    require(bool(seen), "control: the census found no literal kind at all in orpheus/derivations")
    outside: dict[str, Any] = {v: sorted(where) for v, where in seen.items() if v not in members}
    require(
        not outside,
        f"kinds passed by the derivations and absent from ReferenceBoundaryKind {members}: {outside} "
        f"(population: {files} files under orpheus/derivations, {sum(len(v) for v in seen.values())} file-value pairs, "
        f"{non_literal} non-literal values not read)",
    )
