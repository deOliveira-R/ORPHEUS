r"""Content identity of the question values (#405 P1 step 7, S7.9-S7.12).

DRAFT (test-architect, 2026-10-02). Lands as
``tests/gates/numerics/test_content_identity_question.py``; its ``ROSTER`` joins
the union in ``tests/gates/numerics/test_content_identity.py`` (S5.6's route gate reads
it; S5.9 walks ``orpheus.numerics`` and reds while the five classes have no
roster entry, ``[M]`` 2026-10-02 on the prototype).

Every part of every value moves the digest (the population is the type's
``content_parts``, never a hand list); the role is the TYPE, so a fixed
source and a response holding the SAME function at the SAME point are two
values with two digests; the digests are the same in every process, and a
RECORD row pins their bytes.
"""

from __future__ import annotations

import os
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest

from orpheus.numerics.content import FrozenMapping, content_digest
from orpheus.numerics.mesh_free_function import RegionwiseConstant, Symbolic
from orpheus.numerics.question import Eigen, FixedSource, Fundamental, Nearest, Response
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

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/numerics/test_content_identity_question.py"
_CI = "tests/gates/numerics/test_content_identity.py"
_MF = "tests/gates/numerics/test_content_identity_mesh_free.py"
_ROOT = Path(__file__).resolve().parents[3]

_TABLE = np.array([[1.5, 0.25], [0.0, 3.0]])
_POINT = {"fission-emission": -1.0, "boron": 0.125}


def _table(values: np.ndarray = _TABLE) -> RegionwiseConstant:
    return RegionwiseConstant(values.copy())


def _ulp_table() -> RegionwiseConstant:
    moved = _TABLE.copy()
    moved[1, 1] = np.nextafter(moved[1, 1], np.inf)
    return RegionwiseConstant(moved)


def _symbolic() -> Symbolic:
    return Symbolic.of(1 + Symbolic.mu, Symbolic.r**2)


def _point(**changes: float) -> dict[str, float]:
    return {**_POINT, **changes}


def _point_ulp() -> dict[str, float]:
    return _point(boron=float(np.nextafter(0.125, 1.0)))


def _point_legs(build):
    """The point legs every kind shares: one offset by one ULP, a key added, a
    key relabelled (the same offsets under another key), and the str key
    against its int twin (``"1" != 1``, the encoder's rule)."""
    return (
        leg("one offset by one ULP", lambda: build(point=_point_ulp())),
        leg("a key added", lambda: build(point=_point(xenon=0.5))),
        leg("a key relabelled", lambda: build(point={"fission-emission": -1.0, "boron-10": 0.125})),
        leg("the empty point", lambda: build(point={})),
    )


def _eigen(**kw) -> Eigen:
    return Eigen(kw.pop("parameter", "fission-emission"), **{"point": _POINT, "mode": Nearest(0.5), **kw})


def _fixed(**kw) -> FixedSource:
    return FixedSource(kw.pop("source", _table()), **{"point": _POINT, **kw})


def _response(**kw) -> Response:
    return Response(kw.pop("detector", _table()), **{"point": _POINT, **kw})


_EIGEN = Entry(
    cls=Eigen, base=_eigen, parts=("parameter", "point", "mode"),
    perturb={
        "parameter": (
            leg("another key", lambda: _eigen(parameter="boron")),
            leg("a tuple key", lambda: _eigen(parameter=("cell", 1))),
        ),
        "point": _point_legs(_eigen),
        "mode": (
            leg("fundamental", lambda: _eigen(mode=Fundamental())),
            leg("tau by one ULP", lambda: _eigen(mode=Nearest(float(np.nextafter(0.5, 1.0))))),
        ),
    },
    pairs=(
        ("two builds", _eigen, _eigen),
        ("the point's insertion order", _eigen,
         lambda: _eigen(point=dict(reversed(list(_POINT.items()))))),
        ("a dict vs the FrozenMapping", _eigen, lambda: _eigen(point=FrozenMapping(_POINT))),
        ("-0.0 vs 0.0 offset", lambda: _eigen(point=_point(boron=-0.0)), lambda: _eigen(point=_point(boron=0.0))),
        ("an int offset vs its float twin", lambda: _eigen(point=_point(boron=1)), lambda: _eigen(point=_point(boron=1.0))),
        ("tau 1 vs 1.0", lambda: _eigen(mode=Nearest(1)), lambda: _eigen(mode=Nearest(1.0))),
    ),
)

_FIXED = Entry(
    cls=FixedSource, base=_fixed, parts=("source", "point"),
    perturb={
        "source": (
            leg("one coefficient by one ULP", lambda: _fixed(source=_ulp_table())),
            leg("a Symbolic of the same group count", lambda: _fixed(source=_symbolic())),
        ),
        "point": _point_legs(_fixed),
    },
    pairs=(("two builds", _fixed, _fixed),
           ("a dict vs the FrozenMapping", _fixed, lambda: _fixed(point=FrozenMapping(_POINT)))),
)

_RESPONSE = Entry(
    cls=Response, base=_response, parts=("detector", "point"),
    perturb={
        "detector": (
            leg("one coefficient by one ULP", lambda: _response(detector=_ulp_table())),
            leg("a Symbolic of the same group count", lambda: _response(detector=_symbolic())),
        ),
        "point": _point_legs(_response),
    },
    pairs=(("two builds", _response, _response),),
)

_FUNDAMENTAL = Entry(cls=Fundamental, base=Fundamental, parts=(), perturb={},
                     pairs=(("two builds", Fundamental, Fundamental),))

_NEAREST = Entry(
    cls=Nearest, base=lambda: Nearest(0.5), parts=("tau",),
    perturb={"tau": (leg("one ULP", lambda: Nearest(float(np.nextafter(0.5, 1.0)))),
                     leg("the sign", lambda: Nearest(-0.5)))},
    pairs=(("1 vs 1.0", lambda: Nearest(1), lambda: Nearest(1.0)),
           ("-0.0 vs 0.0", lambda: Nearest(-0.0), lambda: Nearest(0.0))),
)

ROSTER: tuple[Entry, ...] = (_EIGEN, _FIXED, _RESPONSE, _FUNDAMENTAL, _NEAREST)


# ── S7.9: every part is content ─────────────────────────────────────────────


@pytest.mark.rests_on(f"{_CI}::test_s5_9_the_population_is_the_roster")
@pytest.mark.parametrize("entry", ROSTER, ids=lambda e: e.id)
def test_s7_9_population(entry: Entry) -> None:
    check_population(entry)


@pytest.mark.rests_on(f"{_MF}::test_s6_9_population")
@pytest.mark.parametrize(
    "entry,part,the_leg",
    [pytest.param(e, p, l, id=param_id(e.id, p, l[0])) for e, p, l in perturbation_ids(ROSTER)],
)
def test_s7_9_a_moved_part_moves_the_digest(entry: Entry, part: str, the_leg) -> None:
    check_perturbation(entry, part, the_leg)


@pytest.mark.parametrize(
    "entry,pair", [pytest.param(e, p, id=param_id(e.id, p[0])) for e, p in pair_ids(ROSTER)]
)
def test_s7_9_equal_content_is_one_value(entry: Entry, pair) -> None:
    check_equal_pair(entry, pair)


@pytest.mark.parametrize("entry", ROSTER, ids=lambda e: e.id)
def test_s7_9_pickle_round_trip(entry: Entry) -> None:
    check_pickle(entry)


# ── S7.10: the role is the type ─────────────────────────────────────────────


@pytest.mark.rests_on(f"{_HERE}::test_s7_9_population")
@pytest.mark.parametrize("function", [_table, _symbolic], ids=["regionwise", "symbolic"])
@pytest.mark.parametrize("point", [{}, _POINT], ids=["physical_point", "offset_point"])
def test_s7_10_a_source_and_a_detector_of_one_function_are_two_questions(function, point) -> None:
    """``FixedSource(q, p)`` and ``Response(q, p)`` hold the SAME function at
    the SAME point and are two values: unequal both ways, two digests, two
    members of a set; and neither is equal to the function it holds. A cache
    can therefore never serve a flux for an importance (ruling 3)."""
    q = function()
    source, detector = FixedSource(q, point), Response(q, point)
    require(source.source == detector.detector, "activation: the two do not hold one function")
    require(not (source == detector) and not (detector == source), "== merged the roles")
    require(source.content_digest != detector.content_digest, "one digest for two roles")
    require(len({source, detector}) == 2, "a set merged the roles")
    require(source != q and detector != q, "a question equals its datum")


# ── S7.11: the same digests in every process ─────────────────────────────────

_SEED_SCRIPT = textwrap.dedent(
    """
    import numpy as np
    import orpheus
    from orpheus.numerics.mesh_free_function import RegionwiseConstant, Symbolic
    from orpheus.numerics.question import Eigen, FixedSource, Fundamental, Nearest, Response

    print("FILE", orpheus.__file__)
    t = RegionwiseConstant(np.array([[1.5, 0.25], [0.0, 3.0]]))
    s = Symbolic.of(1 + Symbolic.mu, Symbolic.r**2)
    values = {
        "eigen": Eigen("fission-emission"),
        "eigen offset nearest": Eigen(("cell", 3), {"boron": 0.125, "fission-emission": -1.0}, Nearest(0.5)),
        "fixed source table": FixedSource(t, {"fission-emission": -1.0}),
        "fixed source symbolic": FixedSource(s),
        "response table": Response(t),
        "fundamental": Fundamental(),
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


@pytest.mark.rests_on(f"{_CI}::test_s5_1_digests_and_hashes_are_seed_stable")
def test_s7_11_digests_and_hashes_are_seed_stable() -> None:
    """``PYTHONHASHSEED`` 1 and 2 print the same six digests and hashes; the
    ``str`` CONTROL hash must differ (the harness sees salting, X1). A ``str``
    parameter key is the leg that would expose a ``hash()``-based encoder."""
    one, two = _run_seeded("1"), _run_seeded("2")
    files = [line.split(" ", 1)[1] for line in one if line.startswith("FILE ")]
    require(files and files[0].startswith(str(_ROOT)), f"the subprocess imported {files} (L22)")
    pick = lambda lines, kind: [l for l in lines if l.startswith(kind + " ")]  # noqa: E731
    for kind in ("DIGEST", "HASH"):
        require(len(pick(one, kind)) == 6, f"{kind}: {len(pick(one, kind))} lines")
        require(pick(one, kind) == pick(two, kind), f"{kind} lines differ across seeds")
    require(pick(one, "CONTROL") != pick(two, "CONTROL"), "control: str hashes did not differ")


# ── S7.12: RECORD, the bytes ─────────────────────────────────────────────────

_PIN_AFTER_LANDING = "pin after landing"
_PINNED: dict[str, str] = {
    "eigen at the physical point": "fd81fd3c2a705a94de9aec7372aed78cdefa090681e073755d5d37785dcbbbd8",
    "eigen offset nearest": "44f0ee840797ce961c39a43cd66c2a7b7b57eb7bc7d27629cb71e4392d33552e",
    "fixed source table": "5278ab77706178e70dfe4afb0d329df10dc1916eb7b6e29b3ce381e344e40fae",
    "response table": "18e27725dbf70e77ae723d64b3fd4f42c766e591bd74e81afcfcdd0099ede33b",
}


def _fingerprint_values() -> dict[str, object]:
    t = RegionwiseConstant(np.array([[1.5, 0.25], [0.0, 3.0]]))
    return {
        "eigen at the physical point": Eigen("fission-emission"),
        "eigen offset nearest": Eigen(("cell", 3), {"boron": 0.125, "fission-emission": -1.0}, Nearest(0.5)),
        "fixed source table": FixedSource(t, {"fission-emission": -1.0}),
        "response table": Response(t),
    }


@pytest.mark.rests_on(f"{_HERE}::test_s7_11_digests_and_hashes_are_seed_stable")
@pytest.mark.rests_on(f"{_CI}::test_s5_7_the_digest_bytes_are_recorded")
def test_s7_12_the_digest_bytes_are_recorded() -> None:
    """RECORD (designed to red on any edit of the question schema or the
    encoder): the producer fingerprint the P3 cache keys rest on. Re-pin only
    with the reason, in the commit message."""
    got = {name: content_digest(value).hex() for name, value in _fingerprint_values().items()}
    if any(pin == _PIN_AFTER_LANDING for pin in _PINNED.values()):
        pytest.fail(f"pin after landing: replace the placeholders of _PINNED with {got}")
    for name, pin in _PINNED.items():
        require(got[name] == pin, f"{name}: the digest bytes moved ({got[name]} != {pin}): every "
                                  f"cache key holding a question is invalidated; re-pin only with the reason")
