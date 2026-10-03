r"""The emission channels, the cell coefficient and the group-count rule (#405 P1 step 8: S8.6, S8.7, S8.9).

Specified by the test-architect (2026-10-02, ``.claude/plans/reference_p1_spec.md``
§1.8). The module is ``orpheus/data/cells.py`` (``Channel``, ``CellCoefficient``);
the group-count rule and ``InconsistentMaterialsError`` live in
``orpheus/data/materials.py`` (ruling 7 of 2026-10-02).

A CELL is ``(material id, Channel)``; ``Channel`` is the closed set of the three
emission channels (ruling 2). A ``CellCoefficient`` is the direction that
scales a set of cells together. It holds explicit ``cells`` and, as a second
field, ``channels_in_every_material`` (spelled ``CellCoefficient.every(*channels)``:
"every material of the problem that carries the channel"), both content, and
``resolve(materials)`` returns the explicit non-zero cells (ruling 5: a
material carries a cell when the cell is non-zero).
"""

from __future__ import annotations

import ast
import enum
from pathlib import Path

import numpy as np
import pytest

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.data.materials import InconsistentMaterialsError, Materials
from orpheus.numerics.content import ContentIdentity
from tests.gates._content_identity_helpers import require
from tests.gates.specification._fixtures import fission_only, fuel, moderator, n2n_only, slab2

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/data/test_cells.py"
_ROOT = Path(__file__).resolve().parents[3]
F, S, N2 = Channel.FISSION_EMISSION, Channel.SCATTERING_EMISSION, Channel.N2N_EMISSION


# ═════════════════════════════════════════════════════════════════════════════
# S8.6 — the channels: a closed set, and the carriage predicate
# ═════════════════════════════════════════════════════════════════════════════


def test_s8_6_the_channels_are_the_three_emissions() -> None:
    """Closed: an ``Enum`` of exactly three members, none an int (an ``IntEnum``
    member would compare equal to its value, and two channels to two ints)."""
    require(issubclass(Channel, enum.Enum) and not issubclass(Channel, int), "Channel is not a plain Enum")
    names = [c.name for c in Channel]
    require(names == ["FISSION_EMISSION", "SCATTERING_EMISSION", "N2N_EMISSION"], f"the channels are {names}")


_CARRIES = {
    # mixture builder -> the channels it carries (non-zero cells), by construction
    fuel: {F, S, N2},
    moderator: {S},
    n2n_only: {N2},
    fission_only: {F},
}


@pytest.mark.parametrize("build", list(_CARRIES), ids=lambda b: b.__name__)
@pytest.mark.parametrize("channel", list(Channel), ids=lambda c: c.name)
def test_s8_6_the_carriage_predicate(build, channel: Channel) -> None:
    """``channel.is_carried_by(mixture)`` is True iff the cell is non-zero: fission
    emission reads ``SigP`` (the producing predicate), scattering any ``SigS``
    block, (n,2n) any ``Sig2`` block. The single-channel fixtures separate the three."""
    expected = channel in _CARRIES[build]
    require(channel.is_carried_by(build()) is expected, f"{channel.name} on {build.__name__}: expected {expected}")


def test_s8_6_a_high_order_block_alone_is_carried() -> None:
    """A stack whose P0 block is zero and whose P1 block is not carries the channel
    (the cell is the whole Legendre stack; a P0-only reader drops it)."""
    from dataclasses import replace

    from scipy.sparse import csr_matrix

    m = moderator()
    p1_only = replace(m, SigS=(csr_matrix((2, 2)), csr_matrix(np.array([[0.01, 0.0], [0.0, 0.02]]))),
                      SigT=m.SigC + m.SigL + m.SigF)
    require(Channel.SCATTERING_EMISSION.is_carried_by(p1_only), "a P1-only scattering stack read as zero")


# ═════════════════════════════════════════════════════════════════════════════
# S8.7 — the cell coefficient: a value, its quantifier, its resolution
# ═════════════════════════════════════════════════════════════════════════════


def test_s8_7_a_cell_coefficient_is_a_content_value() -> None:
    a = CellCoefficient({(0, F), (1, S)})
    require(isinstance(a, ContentIdentity), "not a ContentIdentity")
    require(a == CellCoefficient([(1, S), (0, F), (0, F)]), "order or a duplicate cell is content")  # pyright: ignore[reportArgumentType]  # the coercion is the subject
    require(a == CellCoefficient({(np.int64(0), F), (1, S)}), "a numpy id is another id")  # pyright: ignore[reportArgumentType]  # the coercion is the subject
    require(a != CellCoefficient({(0, F)}) and len({a, CellCoefficient({(0, F)})}) == 2, "a cell is not content")
    require(isinstance(a.cells, frozenset), f"cells is a {type(a.cells).__name__}")


@pytest.mark.parametrize(
    "cells,error,fragment",
    [
        pytest.param(set(), ValueError, r"at least one cell", id="empty"),
        pytest.param({(0,)}, TypeError, r"a \(material id, Channel\) pair", id="not-a-pair"),
        pytest.param({(0, "fission")}, TypeError, r"a channel is a Channel, got a str", id="channel-str"),
        pytest.param({(True, F)}, TypeError, r"a material id is an int, got bool", id="bool-id"),  # the shared parse_integer's wording (#559)
        pytest.param({(0.0, F)}, TypeError, r"a material id is an int, got float", id="float-id"),
    ],
)
def test_s8_7_construction_refusals(cells, error, fragment: str) -> None:
    with pytest.raises(error, match=fragment):
        CellCoefficient(cells)


def test_s8_7_a_channel_in_every_material_is_a_channel() -> None:
    with pytest.raises(TypeError, match=r"a channel is a Channel, got a str"):
        CellCoefficient(channels_in_every_material=["fission emission"])  # pyright: ignore[reportArgumentType]  # the refusal is the subject


def test_s8_7_the_two_fields() -> None:
    """``every`` sets ``channels_in_every_material`` and no cell; explicit cells set
    no channel; both are frozen sets, and either alone is a direction."""
    every = CellCoefficient.every(F, S)
    require(every.cells == frozenset() and every.channels_in_every_material == frozenset({F, S}), f"every is {every!r}")
    explicit = CellCoefficient({(0, F)})
    require(explicit.channels_in_every_material == frozenset(), f"explicit is {explicit!r}")
    mixed = CellCoefficient({(0, F)}, {S})
    require(isinstance(mixed.cells, frozenset) and isinstance(mixed.channels_in_every_material, frozenset), "not frozen")
    require(mixed != explicit and mixed != CellCoefficient.every(S), "a field is not content")


def test_s8_7_every_is_content_with_one_meaning() -> None:
    """``every(*channels)`` is a value: equal for equal channel sets in any
    order, unequal to every explicit set, and refused with no channel."""
    require(CellCoefficient.every(F, S) == CellCoefficient.every(S, F, F), "every's channel order is content")
    require(CellCoefficient.every(F) != CellCoefficient({(0, F)}), "the quantifier equals an explicit cell set")
    require(CellCoefficient.every(F) != CellCoefficient.every(S), "every ignores its channel")
    with pytest.raises(ValueError, match="at least one channel"):
        CellCoefficient.every()


_TWO = Materials({0: fuel(), 1: moderator(), 2: n2n_only()})


@pytest.mark.parametrize(
    "key,resolved",
    [
        pytest.param(CellCoefficient.every(F), {(0, F)}, id="every-F"),
        pytest.param(CellCoefficient.every(S), {(0, S), (1, S)}, id="every-S"),
        pytest.param(CellCoefficient.every(N2), {(0, N2), (2, N2)}, id="every-N2"),
        pytest.param(CellCoefficient.every(F, S, N2), {(0, F), (0, S), (1, S), (0, N2), (2, N2)}, id="every-emission"),
        pytest.param(CellCoefficient({(1, F), (0, F)}), {(0, F)}, id="explicit-zero-cell-dropped"),
        pytest.param(CellCoefficient({(2, N2), (1, S)}), {(2, N2), (1, S)}, id="explicit-kept"),
    ],
)
def test_s8_7_resolution_by_channel(key: CellCoefficient, resolved: set) -> None:
    """``resolve`` returns the explicit NON-ZERO cells, a ``CellCoefficient`` with int ids only."""
    got = key.resolve(_TWO)
    require(got == CellCoefficient(resolved), f"{key} resolved to {got.cells}")
    require(all(isinstance(m, int) and not isinstance(m, bool) for m, _ in got.cells), f"a non-int id survived {got.cells}")
    require(got.channels_in_every_material == frozenset(), f"a channel in every material survived {got!r}")


def test_s8_7_resolution_is_idempotent() -> None:
    for key in (CellCoefficient.every(F, S, N2), CellCoefficient({(1, F), (0, F), (2, N2)})):
        once = key.resolve(_TWO)
        require(once.resolve(_TWO) == once, f"{key}: resolve is not idempotent")


@pytest.mark.parametrize(
    "key,fragment",
    [
        pytest.param(CellCoefficient({(7, F)}), r"material 7 is not among the problem's materials \(ids: \[0, 1, 2\]\)", id="undeclared"),
        pytest.param(CellCoefficient({(0, F), (7, S)}), r"material 7 is not among the problem's materials", id="undeclared-beside-a-cell"),
        pytest.param(CellCoefficient({(9, F)}, {F}),
                     r"material 9 is not among the problem's materials", id="undeclared-beside-every"),
        pytest.param(CellCoefficient({(1, F)}), r"zero direction", id="one-zero-cell"),
        pytest.param(CellCoefficient({(1, F), (2, S)}), r"zero direction", id="all-cells-zero"),
    ],
)
def test_s8_7_resolution_refusals(key: CellCoefficient, fragment: str) -> None:
    with pytest.raises(ValueError, match=fragment):
        key.resolve(_TWO)


def test_s8_7_every_on_a_declaration_carrying_none_is_a_zero_direction() -> None:
    with pytest.raises(ValueError, match="zero direction"):
        CellCoefficient.every(F).resolve(Materials({0: moderator(), 1: n2n_only()}))


# ═════════════════════════════════════════════════════════════════════════════
# S8.9 — the group-count rule: one home in data, every declared material
# ═════════════════════════════════════════════════════════════════════════════


def test_s8_9_the_error_lives_in_data() -> None:
    require(InconsistentMaterialsError.__module__ == "orpheus.data.materials",
            f"InconsistentMaterialsError is defined in {InconsistentMaterialsError.__module__}")
    require(issubclass(InconsistentMaterialsError, ValueError), "no longer a ValueError")


def _imports_of(name: str) -> list[tuple[str, int, str]]:
    hits = []
    for base in ("orpheus", "tests"):
        for path in sorted((_ROOT / base).rglob("*.py")):
            if "__pycache__" in path.parts:
                continue
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            for node in ast.walk(tree):
                if isinstance(node, ast.ImportFrom) and any(a.name == name for a in node.names):
                    hits.append((str(path.relative_to(_ROOT)), node.lineno, node.module or "." * node.level))
    return hits


def test_s8_9_every_importer_reads_the_one_home_with_no_shim() -> None:
    """An AST census over ``orpheus/`` and ``tests/``: every ``from X import
    InconsistentMaterialsError`` names ``orpheus.data.materials``, and no
    package re-exports it (no shim, ruling 7). ``[M]`` first red at ``e2e62cde``:
    4 importers, 0 of them in data (``orpheus/sn/problem.py:46``,
    ``orpheus/transport/mesh/__init__.py:17``, two test modules)."""
    hits = _imports_of("InconsistentMaterialsError")
    require(len(hits) >= 1, "activation: no importer found (the census reads nothing)")
    stray = [h for h in hits if h[2] != "orpheus.data.materials"]
    require(not stray, f"{len(stray)} of {len(hits)} importers read another home: {stray}")


def test_s8_9_the_rule_reads_every_declared_material() -> None:
    """A declaration with a spectator of another group count is refused, by the
    data rule and by ``MaterialMesh`` (which calls it), with the pinned fragment."""
    from orpheus.mesh import CellsByCount, Mesher
    from orpheus.transport.mesh import MaterialMesh

    declaration = Materials({0: fuel(), 5: fuel(ng=3)})
    with pytest.raises(InconsistentMaterialsError, match="uniform ng"):
        declaration.uniform_group_count()
    mesh = Mesher(slab2((0, 0))).partition(CellsByCount.uniform_width(2)).mesh
    with pytest.raises(InconsistentMaterialsError, match="uniform ng"):
        MaterialMesh(mesh, declaration)
    require(Materials({0: fuel(), 5: moderator()}).uniform_group_count() == 2, "a uniform declaration refused")


def test_s8_9_material_mesh_and_the_specification_call_the_one_rule(monkeypatch: pytest.MonkeyPatch) -> None:
    """The ROUTE gate (Pattern 2): replace the data rule with a decoy that raises
    a sentinel, and require both consumers to raise the sentinel. A consumer that
    keeps its own copy of the check stays green under the decoy and reds here.
    Names the helper it targets: ``Materials.uniform_group_count``."""
    from orpheus.data.cells import CellCoefficient
    from orpheus.mesh import CellsByCount, Mesher
    from orpheus.numerics.question import Eigen
    from orpheus.specification import GeometrySpecification, InfiniteMediumSpecification
    from orpheus.transport.mesh import MaterialMesh

    class Sentinel(Exception):
        pass

    calls: list[str] = []

    def decoy(self: Materials) -> int:
        calls.append("called")
        raise Sentinel

    monkeypatch.setattr(Materials, "uniform_group_count", decoy)
    mesh = Mesher(slab2((0, 0))).partition(CellsByCount.uniform_width(2)).mesh
    with pytest.raises(Sentinel):
        MaterialMesh(mesh, {0: fuel()}).ng
    with pytest.raises(Sentinel):
        InfiniteMediumSpecification(0, fuel(), Eigen(CellCoefficient.every(F)))
    with pytest.raises(Sentinel):
        GeometrySpecification(Materials({0: fuel()}), slab2((0, 0)), Eigen(CellCoefficient.every(F)))
    require(len(calls) >= 3, f"activation: the decoy ran {len(calls)} times")
