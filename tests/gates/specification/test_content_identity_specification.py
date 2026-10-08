r"""Content identity of the specification and its keys (#405 P1 step 8: S8.2, S8.5).

Specified by the test-architect (2026-10-02, ``.claude/plans/reference_p1_spec.md``
§1.8). Its ``ROSTER`` (``Specification``, ``CellCoefficient``,
``GeometryExtent``) joins the union in
``tests/gates/numerics/test_content_identity.py``, whose ``_PACKAGES`` gains
``"orpheus.specification"`` in the same change. ``[M]`` on the prototype, the
roster union removed: with the package walked, S5.9 reds naming all three
classes; without the walk it names only ``CellCoefficient`` and
``GeometryExtent``, so ``Specification`` escapes the population gate.

The specification's digest is the digest of its CANONICAL content (S8.10):
two spellings of one direction are one key; two spellings of one FUNCTION
(a ``RegionwiseConstant`` and the ``Symbolic`` it lowers to) are two keys,
a miss and never a wrong hit.
"""

from __future__ import annotations

import os
import subprocess
import sys
import textwrap
from dataclasses import replace
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import sympy as sp

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.data.materials import Materials
from orpheus.geometry import BC, GeometryExtent, StructuredGeometry
from orpheus.numerics.content import FrozenMapping, content_digest
from orpheus.numerics.mesh_free_function import RegionwiseConstant, Symbolic
from orpheus.numerics.question import Eigen, FixedSource, Fundamental, Nearest, Response
from orpheus.specification import GeometrySpecification, InfiniteMediumSpecification
from tests.gates._content_identity_helpers import (
    Entry,
    check_equal_pair,
    check_perturbation,
    check_pickle,
    check_population,
    leg,
    pair_ids,
    param_id,
    perturbation_ids,
    require,
)
from tests.gates.specification._fixtures import fuel, moderator, slab2

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/specification/test_content_identity_specification.py"
_CI = "tests/gates/numerics/test_content_identity.py"
_QCI = "tests/gates/numerics/test_content_identity_question.py"
_SPEC = "tests/gates/specification/test_specification.py"
_ROOT = Path(__file__).resolve().parents[3]

F, S, N2 = Channel.FISSION_EMISSION, Channel.SCATTERING_EMISSION, Channel.N2N_EMISSION
_TABLE = np.array([[1.5, 0.25], [0.0, 3.0]])


def _ulp(x: float) -> float:
    return float(np.nextafter(x, np.inf))


def _materials() -> Materials:
    return Materials({0: fuel(), 1: moderator()})


def _fuel_ulp():
    m = fuel()
    sig_c = m.SigC.copy()
    sig_c[0] = _ulp(sig_c[0])
    return replace(m, SigC=sig_c)


def _absorber():
    """A pure absorber: carries no emission channel, so it moves no canonical key."""
    from tests.gates.specification._fixtures import _mixture

    return _mixture(fission=False, scattering=False, n2n=False)


def _table(values: np.ndarray = _TABLE) -> RegionwiseConstant:
    return RegionwiseConstant(values.copy())


def _table_ulp() -> RegionwiseConstant:
    moved = _TABLE.copy()
    moved[1, 1] = _ulp(moved[1, 1])
    return RegionwiseConstant(moved)


def _point(offset: float = 0.25, key=None):
    return {key if key is not None else CellCoefficient.every(S): offset}


def _spec(materials=None, geometry=None, question=None) -> GeometrySpecification:
    return GeometrySpecification(
        materials=_materials() if materials is None else materials,
        geometry=slab2() if geometry is None else geometry,
        question=FixedSource(_table(), _point()) if question is None else question,
    )


def _medium(material_id: Any = 0, mixture: Any = None, question: Any = None) -> InfiniteMediumSpecification:
    return InfiniteMediumSpecification(
        material_id, fuel() if mixture is None else mixture,
        Response(_table(_TABLE[:1]), _point()) if question is None else question,
    )


_INFINITE_MEDIUM = Entry(
    cls=InfiniteMediumSpecification, base=_medium, parts=("material_id", "mixture", "question"),
    perturb={
        # The material id moves the resolved point key too (its cells name the id).
        "material_id": (leg("another id", lambda: _medium(material_id=3), "question"),),
        "mixture": (leg("one coefficient by one ULP", lambda: _medium(mixture=_fuel_ulp())),),
        "question": (
            leg("the detector's second coefficient by one ULP", lambda: _medium(question=Response(_table(np.array([[_TABLE[0, 0], _ulp(_TABLE[0, 1])]])), _point()))),
            leg("the point's offset by one ULP", lambda: _medium(question=Response(_table(_TABLE[:1]), _point(_ulp(0.25))))),
            leg("the role", lambda: _medium(question=FixedSource(_table(_TABLE[:1]), _point()))),
            leg("k", lambda: _medium(question=Eigen(CellCoefficient.every(F)))),
        ),
    },
    pairs=(
        ("two builds", _medium, _medium),
        ("an int id vs a numpy id", _medium, lambda: _medium(material_id=np.int64(0))),
        ("every(S) vs its explicit cell", _medium,
         lambda: _medium(question=Response(_table(_TABLE[:1]), _point(key=CellCoefficient({(0, S)}))))),
    ),
)

_SPECIFICATION = Entry(
    cls=GeometrySpecification, base=_spec, parts=("materials", "geometry", "question"),
    perturb={
        "materials": (
            leg("one coefficient by one ULP", lambda: _spec(materials=Materials({0: _fuel_ulp(), 1: moderator()}))),
        ),
        "geometry": (
            leg("a breakpoint by one ULP", lambda: _spec(
                geometry=StructuredGeometry.slab((0.0, _ulp(1.0), 3.0), (0, 1), left=BC.reflective, right=BC.vacuum))),
            leg("a boundary law", lambda: _spec(
                geometry=StructuredGeometry.slab((0.0, 1.0, 3.0), (0, 1), left=BC.vacuum, right=BC.vacuum))),
        ),
        "question": (
            leg("one source coefficient by one ULP", lambda: _spec(question=FixedSource(_table_ulp(), _point()))),
            leg("the point's offset by one ULP", lambda: _spec(question=FixedSource(_table(), _point(_ulp(0.25))))),
            leg("the point key's channel", lambda: _spec(question=FixedSource(_table(), _point(key=CellCoefficient.every(N2))))),
            leg("the point key's kind", lambda: _spec(question=FixedSource(_table(), _point(key=GeometryExtent(0))))),
            leg("the role (a response of the same function)", lambda: _spec(question=Response(_table(), _point()))),
            leg("k", lambda: _spec(question=Eigen(CellCoefficient.every(F)))),
        ),
    },
    pairs=(
        ("two builds", _spec, _spec),
        # The orchestrator's ruling 3 on the step-8 NEEDS (2026-10-02): the
        # specification keeps only the materials its geometry assigns, so a
        # spectator, carrying the channel or not, is not part of the key.
        ("a spectator declared", _spec, lambda: _spec(materials=Materials({0: fuel(), 1: moderator(), 9: moderator()}))),
        ("a non-carrying spectator declared", _spec, lambda: _spec(materials=Materials({0: fuel(), 1: moderator(), 9: _absorber()}))),
        ("every(S) vs its explicit cells", _spec,
         lambda: _spec(question=FixedSource(_table(), _point(key=CellCoefficient({(0, S), (1, S)}))))),
        ("an explicit zero cell vs none",
         lambda: _spec(question=Eigen(CellCoefficient({(0, F)}))),
         lambda: _spec(question=Eigen(CellCoefficient({(0, F), (1, F)})))),
        ("the default point vs the empty point",
         lambda: _spec(question=Eigen(CellCoefficient.every(F))),
         lambda: _spec(question=Eigen(CellCoefficient.every(F), {}, Fundamental()))),
    ),
)

def _cc(cells: Iterable[tuple[int, Channel]] = frozenset({(0, F), (1, S)}),
        every: Iterable[Channel] = frozenset({N2})) -> CellCoefficient:
    return CellCoefficient(cells, every)


_CELL_COEFFICIENT = Entry(
    cls=CellCoefficient, base=_cc, parts=("cells", "channels_in_every_material"),
    perturb={
        "cells": (
            leg("a cell's channel", lambda: _cc(cells={(0, F), (1, N2)})),
            leg("a cell's material", lambda: _cc(cells={(0, F), (2, S)})),
            leg("a cell added", lambda: _cc(cells={(0, F), (1, S), (1, N2)})),
            leg("no cell", lambda: _cc(cells=frozenset())),
        ),
        "channels_in_every_material": (
            leg("another channel", lambda: _cc(every={S})),
            leg("a channel added", lambda: _cc(every={N2, F})),
            leg("none", lambda: _cc(every=frozenset())),
        ),
    },
    pairs=(("two builds", _cc, lambda: CellCoefficient([(1, S), (0, F), (0, F)], [N2, N2])),  # pyright: ignore[reportArgumentType]  # the coercion is the subject
           ("every in any channel order", lambda: CellCoefficient.every(F, S), lambda: CellCoefficient.every(S, F)),
           ("every(...) is the field", lambda: CellCoefficient.every(F, S), lambda: CellCoefficient(channels_in_every_material={F, S}))),
)

_GEOMETRY_EXTENT = Entry(
    cls=GeometryExtent, base=lambda: GeometryExtent(1), parts=("interval",),
    perturb={"interval": (leg("another interval", lambda: GeometryExtent(0)),)},
    pairs=(("an int vs a numpy int", lambda: GeometryExtent(1), lambda: GeometryExtent(np.int64(1))),),  # pyright: ignore[reportArgumentType]  # the coercion is the subject
)

ROSTER: tuple[Entry, ...] = (_SPECIFICATION, _INFINITE_MEDIUM, _CELL_COEFFICIENT, _GEOMETRY_EXTENT)


# ── S8.2: every part is content; equal content is one value ──────────────────


@pytest.mark.rests_on(f"{_CI}::test_s5_9_the_population_is_the_roster")
@pytest.mark.parametrize("entry", ROSTER, ids=lambda e: e.id)
def test_s8_2_population(entry: Entry) -> None:
    check_population(entry)


@pytest.mark.rests_on(f"{_HERE}::test_s8_2_population", f"{_QCI}::test_s7_9_population")
@pytest.mark.parametrize(
    "entry,part,the_leg",
    [pytest.param(e, p, l, id=param_id(e.id, p, l[0])) for e, p, l in perturbation_ids(ROSTER)],
)
def test_s8_2_a_moved_part_moves_the_digest(entry: Entry, part: str, the_leg) -> None:
    check_perturbation(entry, part, the_leg)


@pytest.mark.rests_on(f"{_SPEC}::test_s8_10_every_resolves_to_the_explicit_non_zero_cells")
@pytest.mark.parametrize(
    "entry,pair", [pytest.param(e, p, id=param_id(e.id, p[0])) for e, p in pair_ids(ROSTER)]
)
def test_s8_2_equal_content_is_one_value(entry: Entry, pair) -> None:
    check_equal_pair(entry, pair)


@pytest.mark.parametrize("entry", ROSTER, ids=lambda e: e.id)
def test_s8_2_pickle_round_trip(entry: Entry) -> None:
    check_pickle(entry)


# ── S8.2: one function in two spellings is two keys ──────────────────────────


def test_s8_2_a_table_and_the_symbolic_it_lowers_to_are_one_rhs_and_two_keys() -> None:
    """A ``RegionwiseConstant`` ``Q`` and the ``Symbolic`` ``q = Q`` written as a
    ``Piecewise`` in r on the same geometry define the same right-hand side,
    checked in Branch 1 at points strictly inside each region (never at a
    breakpoint, §0 ruling 4); their specifications' digests DIFFER: two
    spellings of one function are two cache keys (a miss, never a wrong hit;
    census Q5.6: ``ContentIdentity.__eq__`` is ``NotImplemented`` across types)."""
    r = Symbolic.r
    geometry = slab2()
    (a, b), (_, c) = geometry.intervals
    per_group = [sp.Piecewise((float(_TABLE[0, g]), r < b), (float(_TABLE[1, g]), True)) for g in range(2)]
    symbolic = Symbolic.of(*per_group)
    for x in (0.5 * (a + b), 0.5 * (b + c), a + 1e-3 * (b - a), c - 1e-3 * (c - b)):
        region = 0 if x < b else 1
        values = [float(e.subs(r, x)) for e in symbolic.expressions]
        require(values == list(_TABLE[region]), f"x = {x}: the Symbolic reads {values}, the table {_TABLE[region]}")
    by_table = _spec(question=FixedSource(_table()))
    by_symbol = _spec(question=FixedSource(symbolic))
    require(by_table.content_digest != by_symbol.content_digest and by_table != by_symbol,
            "two spellings of one function share a key")


# ── S8.2: the same digests in every process ──────────────────────────────────

_SEED_SCRIPT = textwrap.dedent(
    """
    import numpy as np
    import orpheus
    from orpheus.data.cells import CellCoefficient, Channel
    from orpheus.data.materials import Materials
    from orpheus.geometry import GeometryExtent
    from orpheus.numerics.mesh_free_function import RegionwiseConstant
    from orpheus.numerics.question import Eigen, FixedSource
    from orpheus.specification import GeometrySpecification, InfiniteMediumSpecification
    from tests.gates.specification._fixtures import fuel, moderator, slab2

    print("FILE", orpheus.__file__)
    t = RegionwiseConstant(np.array([[1.5, 0.25], [0.0, 3.0]]))
    F, S = Channel.FISSION_EMISSION, Channel.SCATTERING_EMISSION
    values = {
        "k infinite medium": InfiniteMediumSpecification(0, fuel(), Eigen(CellCoefficient.every(F))),
        "fixed source slab": GeometrySpecification(materials=Materials({0: fuel(), 1: moderator()}), geometry=slab2(),
                                           question=FixedSource(t, {CellCoefficient.every(S): 0.25})),
        "critical extent": GeometrySpecification(materials=Materials({0: fuel(), 1: moderator()}), geometry=slab2(),
                                         question=Eigen(GeometryExtent(1))),
        "a cell coefficient": CellCoefficient.every(F, S),
    }
    for name, value in values.items():
        print("DIGEST", name, value.content_digest.hex())
        print("HASH", name, hash(value))
    print("CONTROL", hash("a salted str hash"))
    """
)


def _run_seeded(seed: str) -> list[str]:
    env = {**os.environ, "PYTHONHASHSEED": seed, "PYTHONPATH": str(_ROOT)}
    out = subprocess.run([sys.executable, "-O", "-c", _SEED_SCRIPT], cwd=_ROOT, env=env,
                         capture_output=True, text=True, timeout=300)
    require(out.returncode == 0, f"seed {seed}: {out.stderr[-3000:]}")
    return out.stdout.splitlines()


@pytest.mark.rests_on(f"{_QCI}::test_s7_11_digests_and_hashes_are_seed_stable")
def test_s8_2_digests_and_hashes_are_seed_stable() -> None:
    """``PYTHONHASHSEED`` 1 and 2 print the same four digests and hashes; the
    ``str`` CONTROL must differ (the harness sees salting, X1). The cell set is a
    ``frozenset`` of enum-bearing tuples, whose ITERATION order is seed-dependent:
    the leg that would expose an encoder reading a set in iteration order."""
    one, two = _run_seeded("1"), _run_seeded("2")
    files = [line.split(" ", 1)[1] for line in one if line.startswith("FILE ")]
    require(files and files[0].startswith(str(_ROOT)), f"the subprocess imported {files} (L22)")
    pick = lambda lines, kind: [l for l in lines if l.startswith(kind + " ")]  # noqa: E731
    for kind in ("DIGEST", "HASH"):
        require(len(pick(one, kind)) == 4, f"{kind}: {len(pick(one, kind))} lines")
        require(pick(one, kind) == pick(two, kind), f"{kind} lines differ across seeds")
    require(pick(one, "CONTROL") != pick(two, "CONTROL"), "control: str hashes did not differ")


# ── S8.5: RECORD, the bytes ──────────────────────────────────────────────────

_PIN_AFTER_LANDING = "pin after landing"
# Re-pinned 2026-10-02 (the test-architect): the layer became the type
# (``InfiniteMediumSpecification`` / ``GeometrySpecification``, user ruling "The
# infinite medium is the point in phase space") and ``CellCoefficient`` gained
# ``channels_in_every_material``; both schema tags moved, so both pins moved.
_PINNED: dict[str, str] = {
    # Re-pinned 2026-10-08: Eigen gained its declared gauge (the user's ruling), which the specification's
    # canonical question writes in (the default production, resolved); the fixed-source pin did not move.
    "k infinite medium": "6eb99bfc997eca9114e9b62a3d59124067d28374eefd88ba7ed1bf4349528007",
    "fixed source slab": "737d57afd4df1d4f9a63d0bcf963acbe7404bb8102fbbf70fab6f00d3bdd6366",
}


def _fingerprint_values() -> dict[str, object]:
    t = RegionwiseConstant(_TABLE.copy())
    return {
        "k infinite medium": InfiniteMediumSpecification(0, fuel(), Eigen(CellCoefficient.every(F))),
        "fixed source slab": _spec(question=FixedSource(t, _point())),
    }


@pytest.mark.rests_on(f"{_HERE}::test_s8_2_digests_and_hashes_are_seed_stable")
@pytest.mark.rests_on(f"{_QCI}::test_s7_12_the_digest_bytes_are_recorded")
def test_s8_5_the_digest_bytes_are_recorded() -> None:
    """RECORD (designed to red on any edit of the specification's schema, the
    keys' schema, the canonical form or the encoder): the producer fingerprint the
    P3 cache keys rest on. Re-pin only with the reason, in the commit message."""
    got = {name: content_digest(value).hex() for name, value in _fingerprint_values().items()}
    if any(pin == _PIN_AFTER_LANDING for pin in _PINNED.values()):
        pytest.fail(f"pin after landing: replace the placeholders of _PINNED with {got}")
    for name, pin in _PINNED.items():
        require(got[name] == pin, f"{name}: the digest bytes moved ({got[name]} != {pin}): every "
                                  f"cache key holding a specification is invalidated; re-pin only with the reason")
