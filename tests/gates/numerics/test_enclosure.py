r"""The enclosure of a real number and its quotient (#405 P2 step 1, R1.1-R1.9).

Specified by the test-architect (2026-10-03, ``.claude/plans/reference_p2_spec.md``
§1.1). ``Enclosure(value, bound)`` claims that an exact real number lies in
the closed interval :math:`[v - b, v + b]`, read in EXACT real arithmetic on
the two doubles; a reference's reading of its own exact answer is one, and an
exact value is the enclosure with bound 0. The quotient of two enclosures
encloses every quotient of their members (the containment law), so a
reference's bound propagates through a ``Ratio`` observable without a second
definition of the bound.

The laws are checked in exact rational arithmetic (``fractions.Fraction``):
no row evaluates a bound in floating point, so no row can be the victim of
the rounding it gates (``tests/_harness/float_bounds.py``, the same stance).

Every class is resolved on its module at run time (``_E()``), never captured
at collection, so a mutation battery that rebinds it reaches every row
(test-architect lessons §2, ``L100``).
"""

from __future__ import annotations

import importlib
import math
import os
import pickle
import random
import subprocess
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import pytest

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

_HERE = "tests/gates/numerics/test_enclosure.py"
_ROOT = Path(__file__).resolve().parents[3]
_MODULE = "orpheus.numerics.enclosure"


def _E() -> Any:
    """The class, read off its module now (a rebinding arm reaches it)."""
    return importlib.import_module(_MODULE).Enclosure


def _hull(e: Any) -> tuple[Fraction, Fraction]:
    """The enclosure's interval in exact arithmetic on its two doubles."""
    v, b = Fraction(e.value), Fraction(e.bound)
    return v - b, v + b


def _quotient_corners(a: Any, b: Any) -> list[Fraction]:
    """The exact quotients at the four corners of the operand box; with zero
    outside the denominator's interval, x / y is monotone in each variable on
    the box, so these four are the extremes of every quotient of members."""
    (xa, xb), (ya, yb) = _hull(a), _hull(b)
    return [x / y for x in (xa, xb) for y in (ya, yb)]


def _exact_half_width(a: Any, b: Any, centre: float) -> Fraction:
    """The smallest bound about ``centre`` that encloses every quotient, exactly."""
    return max(abs(q - Fraction(centre)) for q in _quotient_corners(a, b))


# ── the population ──────────────────────────────────────────────────────────

_SEED = 20261003
_N_DRAWN = 3000
#: Relative widths spanning an exact value, rounding-sized, small, and large
#: enclosures; the denominator stops below 1 so zero stays outside it.
_REL_NUM = (0.0, 1e-16, 1e-8, 1e-2, 0.5, 0.99)
_REL_DEN = (0.0, 1e-16, 1e-8, 1e-2, 0.5, 0.9)


def _draw(rng: random.Random, rels: tuple[float, ...]) -> tuple[float, float]:
    value = rng.choice((-1.0, 1.0)) * 10.0 ** rng.uniform(-30, 30) * rng.uniform(1, 10)
    return value, abs(value) * rng.choice(rels) * rng.uniform(0, 1)


#: Hand-picked rows each population draw is unlikely to hit: exact operands
#: whose quotient is not a double (outward rounding is then the whole bound),
#: an exact quotient, every sign combination with wide enclosures (the
#: corner choice), and a denominator one part in 10^15 from touching zero.
_SPECIAL: tuple[tuple[str, tuple[float, float], tuple[float, float]], ...] = (
    ("one-third-exact", (1.0, 0.0), (3.0, 0.0)),
    ("two-sevenths-exact", (2.0, 0.0), (7.0, 0.0)),
    ("one-half-exact", (1.0, 0.0), (2.0, 0.0)),
    ("pos-pos-wide", (3.0, 2.5), (2.0, 1.5)),
    ("pos-neg-wide", (3.0, 2.5), (-2.0, 1.5)),
    ("neg-pos-wide", (-3.0, 2.5), (2.0, 1.5)),
    ("neg-neg-wide", (-3.0, 2.5), (-2.0, 1.5)),
    ("numerator-straddles-zero", (0.5, 2.0), (4.0, 1.0)),
    ("denominator-near-zero", (1.0, 0.0), (1.0, 1.0 - 1e-15)),
    ("tiny-over-huge", (1e-300, 1e-310), (1e300, 1e290)),
)


def _population() -> list[tuple[str, tuple[float, float], tuple[float, float]]]:
    rng = random.Random(_SEED)
    drawn = [(f"draw{i}", _draw(rng, _REL_NUM), _draw(rng, _REL_DEN)) for i in range(_N_DRAWN)]
    return list(_SPECIAL) + drawn


def _quotients() -> list[tuple[str, Any, Any, Any]]:
    """Every population row whose quotient constructs, with the operands."""
    E = _E()
    rows = []
    for name, (va, ba), (vb, bb) in _population():
        a, b = E(va, ba), E(vb, bb)
        try:
            rows.append((name, a, b, a / b))
        except (ZeroDivisionError, ValueError):
            continue
    require(len(rows) >= 2900, f"activation: only {len(rows)} of {len(_SPECIAL) + _N_DRAWN} rows formed a quotient")
    return rows


# ═════════════════════════════════════════════════════════════════════════════
# R1.1 — admission
# ═════════════════════════════════════════════════════════════════════════════

_REFUSED: tuple[tuple[str, tuple[Any, Any], type[BaseException], str], ...] = (
    ("nan-value", (math.nan, 0.0), ValueError, "Enclosure: the value is NaN"),
    ("inf-value", (math.inf, 0.0), ValueError, "Enclosure: the value is infinite"),
    ("nan-bound", (1.0, math.nan), ValueError, "Enclosure: the bound is NaN"),
    ("inf-bound", (1.0, math.inf), ValueError, "Enclosure: the bound is infinite"),
    ("negative-bound", (1.0, -1e-300), ValueError, "the bound is a distance"),
    ("complex-value", (1j, 0.0), TypeError, "Enclosure: the value must be a real number"),
    ("str-bound", (1.0, "0"), TypeError, "Enclosure: the bound must be a real number"),
    ("none-value", (None, 0.0), TypeError, "Enclosure: the value must be a real number"),
    ("bool-bound", (1.0, True), TypeError, "Enclosure: the bound must be a real number"),
)


@pytest.mark.parametrize("args, error, fragment", [r[1:] for r in _REFUSED], ids=[r[0] for r in _REFUSED])
def test_r1_1_admission_refuses(args: tuple[Any, Any], error: type[BaseException], fragment: str) -> None:
    """A value and a bound are finite reals and the bound is non-negative;
    each refusal names the field (the message is the gate)."""
    with pytest.raises(error, match=fragment):
        _E()(*args)


def test_r1_1_admission_admits() -> None:
    """Positive legs (anti-#11): an exact value, an integer pair (stored as
    doubles), a negative value with a wide bound, and ``-0.0`` as a bound
    (canonical ``+0.0``: a distance has no sign)."""
    E = _E()
    exact = E(2.0, 0.0)
    require((exact.value, exact.bound) == (2.0, 0.0), f"exact: {exact!r}")
    ints = E(3, 1)
    require(type(ints.value) is float and type(ints.bound) is float, f"integers stored as {type(ints.value)}, {type(ints.bound)}")
    require(E(-3.5, 2.0).value == -3.5, "a negative value is admitted")
    signed_zero = E(1.0, -0.0)
    require(math.copysign(1.0, signed_zero.bound) == 1.0, f"-0.0 bound stored as {signed_zero.bound!r}")


# ═════════════════════════════════════════════════════════════════════════════
# R1.2 — the containment law of the quotient (THEOREM)
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.rests_on(f"{_HERE}::test_r1_1_admission_admits")
def test_r1_2_the_quotient_encloses_every_quotient() -> None:
    r"""For every population row: every exact quotient :math:`x / y`, with
    :math:`x` in the numerator's interval and :math:`y` in the denominator's,
    lies in the quotient's interval (checked at the four corners, which bound
    the quotient of members because :math:`x/y` is monotone in each variable
    on a box that excludes :math:`y = 0`). Population: the 10 special rows and
    3000 seeded draws (seed 20261003); the row prints how many formed."""
    rows = _quotients()
    failures = [
        f"{name}: {a!r} / {b!r} = {r!r} misses {float(_exact_half_width(a, b, r.value) - Fraction(r.bound)):.3e}"
        for name, a, b, r in rows
        if _exact_half_width(a, b, r.value) > Fraction(r.bound)
    ]
    print(f"R1.2: {len(rows)} quotients checked")
    require(not failures, f"{len(failures)} of {len(rows)} quotients fail to enclose:\n" + "\n".join(failures[:10]))


@pytest.mark.rests_on(f"{_HERE}::test_r1_2_the_quotient_encloses_every_quotient")
def test_r1_2_an_exact_quotient_that_is_not_a_double_has_a_positive_bound() -> None:
    """``Enclosure(1, 0) / Enclosure(3, 0)``: one third is not a double, so the
    quotient's bound is positive, and it contains one third exactly. The row
    the outward rounding alone carries (a rounded-to-nearest bound is 0)."""
    E = _E()
    third = E(1.0, 0.0) / E(3.0, 0.0)
    require(third.bound > 0.0, f"1/3 is not a double, so its exact enclosure has a positive bound; got {third!r}")
    require(abs(Fraction(1, 3) - Fraction(third.value)) <= Fraction(third.bound), f"{third!r} does not contain 1/3")


# ═════════════════════════════════════════════════════════════════════════════
# R1.3 — the centre is the quotient of the values; R1.4 — tightness
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.rests_on(f"{_HERE}::test_r1_2_the_quotient_encloses_every_quotient")
def test_r1_3_the_centre_is_the_quotient_of_the_values() -> None:
    """The quotient's value is ``a.value / b.value`` bit for bit: a ratio
    observable reads the ratio of its readings, and the bound carries the
    rest (a centre at the hull's midpoint would move the reading with the
    operands' bounds)."""
    bad = [f"{n}: {r.value!r} != {a.value / b.value!r}" for n, a, b, r in _quotients() if r.value != a.value / b.value]
    require(not bad, f"{len(bad)} centres moved:\n" + "\n".join(bad[:10]))


#: The derived slack of an outward-rounded quotient, in units of the last
#: place of the largest exact corner |q|: each operand endpoint is rounded
#: once and stepped outward once (at most 1.5 ulp of itself, 3u relative), so
#: a corner quotient is perturbed by at most 6u relative, rounded (u) and
#: stepped outward (2u): 9u|q| <= 9 ulp(q); the half-width's subtraction is
#: rounded and stepped up, 1.5 ulp(q) more. ``[M]`` 2026-10-03 on the
#: prototype over this population: 7.84.
_SLACK_ULP = Fraction(21, 2)


@pytest.mark.rests_on(f"{_HERE}::test_r1_3_the_centre_is_the_quotient_of_the_values")
def test_r1_4_the_bound_is_the_hull_up_to_its_own_rounding() -> None:
    """The bound exceeds the exact half-width about the centre by at most
    10.5 ulp of the largest corner: an enclosure that is merely valid (a bound
    of 1e300 encloses everything) makes every verification floor fail, so the
    law has two sides."""
    bad = []
    for name, a, b, r in _quotients():
        corners = _quotient_corners(a, b)
        unit = Fraction(math.ulp(float(max(abs(q) for q in corners))))
        excess = (Fraction(r.bound) - _exact_half_width(a, b, r.value)) / unit
        if excess > _SLACK_ULP:
            bad.append(f"{name}: excess {float(excess):.2f} ulp")
    require(not bad, f"{len(bad)} bounds looser than their own rounding:\n" + "\n".join(bad[:10]))


# ═════════════════════════════════════════════════════════════════════════════
# R1.5 — the refusals of the quotient
# ═════════════════════════════════════════════════════════════════════════════

_ZERO_DENOMINATORS = (
    ("exact-zero", (0.0, 0.0)),
    ("touches-zero-from-above", (1.0, 1.0)),
    ("straddles-zero", (1.0, 2.0)),
    ("straddles-zero-negative", (-1.0, 1.5)),
)


@pytest.mark.parametrize("den", [d for _, d in _ZERO_DENOMINATORS], ids=[n for n, _ in _ZERO_DENOMINATORS])
def test_r1_5_a_denominator_holding_zero_is_refused(den: tuple[float, float]) -> None:
    """``ZeroDivisionError`` naming the reason: a quotient over an interval
    holding zero is unbounded, so there is no enclosure to return."""
    E = _E()
    with pytest.raises(ZeroDivisionError, match="contains zero"):
        E(1.0, 0.0) / E(*den)


def test_r1_5_an_overflowing_quotient_is_refused() -> None:
    """``1e300 / 1e-300`` is beyond a double: the quotient is refused, never an
    enclosure with an infinite bound (an infinite bound encloses nothing a
    floor can use)."""
    E = _E()
    with pytest.raises(ValueError, match="infinite"):
        E(1e300, 0.0) / E(1e-300, 0.0)


@pytest.mark.parametrize("other", [3.0, 3, None], ids=["float", "int", "none"])
def test_r1_5_an_exact_operand_is_spelled_as_an_enclosure(other: Any) -> None:
    """A bare number is not an enclosure: dividing by one is a ``TypeError``
    both ways, so an exact value is spelled ``Enclosure(v, 0)`` once and a
    float can never stand in for a bounded reading."""
    E = _E()
    with pytest.raises(TypeError):
        E(1.0, 0.0) / other
    with pytest.raises(TypeError):
        other / E(1.0, 0.0)


# ═════════════════════════════════════════════════════════════════════════════
# R1.6 — the fields are the claim; R1.7 — content identity
# ═════════════════════════════════════════════════════════════════════════════


def test_r1_6_the_fields_are_the_value_and_the_bound() -> None:
    """``dataclasses.fields`` is exactly ``(value, bound)``: the interval's
    ends are derived, never stored beside the centre (one source, X4)."""
    import dataclasses

    names = tuple(f.name for f in dataclasses.fields(_E()))
    require(names == ("value", "bound"), f"fields {names}")


def _enc(value: float = 1.25, bound: float = 0.5) -> Any:
    return _E()(value, bound)


ROSTER: tuple[Entry, ...] = (
    Entry(
        cls=importlib.import_module(_MODULE).Enclosure,
        base=_enc,
        parts=("value", "bound"),
        perturb={
            "value": (
                leg("one ULP", lambda: _enc(value=math.nextafter(1.25, 2.0))),
                leg("sign", lambda: _enc(value=-1.25)),
            ),
            "bound": (
                leg("one ULP", lambda: _enc(bound=math.nextafter(0.5, 1.0))),
                leg("exact", lambda: _enc(bound=0.0)),
            ),
        },
        pairs=(
            ("int and float", lambda: _E()(2, 1), lambda: _E()(2.0, 1.0)),
            ("bound -0.0 and 0.0", lambda: _E()(2.0, -0.0), lambda: _E()(2.0, 0.0)),
            ("value -0.0 and 0.0", lambda: _E()(-0.0, 1.0), lambda: _E()(0.0, 1.0)),
        ),
    ),
)


@pytest.mark.parametrize("entry", ROSTER, ids=[e.id for e in ROSTER])
def test_r1_7_the_population_is_the_parts(entry: Entry) -> None:
    check_population(entry)


@pytest.mark.parametrize(
    "entry, part, the_leg",
    perturbation_ids(ROSTER),
    ids=[param_id(e.id, p, lg[0]) for e, p, lg in perturbation_ids(ROSTER)],
)
def test_r1_7_each_part_moves_the_digest(entry: Entry, part: str, the_leg: Any) -> None:
    check_perturbation(entry, part, the_leg)


@pytest.mark.parametrize("entry, pair", pair_ids(ROSTER), ids=[param_id(e.id, p[0]) for e, p in pair_ids(ROSTER)])
def test_r1_7_equal_content_is_one_value(entry: Entry, pair: Any) -> None:
    check_equal_pair(entry, pair)


@pytest.mark.parametrize("entry", ROSTER, ids=[e.id for e in ROSTER])
def test_r1_7_pickle(entry: Entry) -> None:
    check_pickle(entry)


# ═════════════════════════════════════════════════════════════════════════════
# R1.8 — the layer
# ═════════════════════════════════════════════════════════════════════════════

_LAYER_SCRIPT = """
import sys
import orpheus.numerics.enclosure as m
e = m.Enclosure(1.0, 0.0) / m.Enclosure(3.0, 0.5)
print(m.__file__)
print(sorted({k.split('.')[1] for k in sys.modules if k.startswith('orpheus.')}))
print('sympy' in sys.modules)
"""


def test_r1_8_the_module_is_numerics() -> None:
    """(a) A fresh interpreter that imports the module and forms one quotient
    loads only the ``orpheus`` sub-packages a cold ``import orpheus.numerics``
    loads (``{geometry, numerics}``, ``[M]`` 2026-10-03) and no SymPy. (b) By
    AST, every ``orpheus`` import in the module (module-level, local,
    ``TYPE_CHECKING``) is under ``orpheus.numerics``; activation: the import
    of ``orpheus.numerics.content`` is seen."""
    import ast

    env = {**os.environ, "PYTHONPATH": str(_ROOT)}
    out = subprocess.run([sys.executable, "-O", "-c", _LAYER_SCRIPT], cwd=_ROOT, env=env, capture_output=True, text=True)
    require(out.returncode == 0, out.stderr)
    file, packages, sympy_loaded = out.stdout.strip().splitlines()
    require(file.startswith(str(_ROOT)), f"the subprocess imported {file} (L22)")
    require(packages == "['geometry', 'numerics']", f"loaded orpheus packages {packages}")
    require(sympy_loaded == "False", "the module loads SymPy")
    module_file = importlib.import_module(_MODULE).__file__
    require(module_file is not None, "the module has no file")
    ours: list[str] = []
    for node in ast.walk(ast.parse(Path(str(module_file)).read_text())):
        if isinstance(node, ast.ImportFrom) and node.module:
            ours.append(node.module)
        elif isinstance(node, ast.Import):
            ours.extend(alias.name for alias in node.names)
    ours = [m for m in ours if m.startswith("orpheus")]
    require("orpheus.numerics.content" in ours, f"activation: the content import is not seen in {ours}")
    outside = [m for m in ours if not m.startswith("orpheus.numerics")]
    require(not outside, f"imports outside numerics: {outside}")
