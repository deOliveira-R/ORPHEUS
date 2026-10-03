r"""The readings, by role (#405 P2 step 1, R1.10-R1.19).

Specified by the test-architect (2026-10-03, ``.claude/plans/reference_p2_spec.md``
§1.1), on the user's G3 ruling and the orchestrator's Q1/Q4 rulings of
2026-10-03:

* a REFERENCE's reading is a guarantee about the exact answer of its own
  equation, ``ReferenceReading = Enclosure | Printed``: an
  :class:`~orpheus.numerics.enclosure.Enclosure`, or a :class:`Printed` value,
  the cited author's claim, never recomputed;
* a PRODUCTION reading is a self-report verification puts on trial,
  ``ProductionReading = Measured`` (``Estimated`` joins it when a Monte Carlo
  consumer exists); the two sums are disjoint, so an ``Estimated`` (or any
  production) reference is unspellable;
* ``Printed(text, citation)``: the printed decimal string is the one source;
  its value (the correctly rounded double) and its half unit in the last
  printed digit are DERIVED, and the citation carries the place the value is
  printed (a locator).

The package ``orpheus.reference`` is new: an input-tier package above the
specification, registered in ``tests/gates/test_layer_imports.py``. Until it
is registered the import linter SKIPS it (``_check_source`` returns ``[]``
for a package it does not know), which R1.18 reads as a red.
"""

from __future__ import annotations

import ast
import dataclasses
import importlib
import math
import os
import subprocess
import sys
import typing
from decimal import Decimal
from pathlib import Path
from typing import Any

import pytest

from orpheus.data.citation import Citation
from orpheus.numerics.content import content_digest
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

_HERE = "tests/gates/reference/test_readings.py"
_ENC = "tests/gates/numerics/test_enclosure.py"
_ROOT = Path(__file__).resolve().parents[3]
_READING = "orpheus.reference.reading"
_CITE = Citation("SoodForsterParsons2003", "Table 10")


def _reading() -> Any:
    return importlib.import_module(_READING)


def _P(text: Any, citation: Any = _CITE) -> Any:
    return _reading().Printed(text, citation)


# ═════════════════════════════════════════════════════════════════════════════
# R1.10 — the printed text is canonical: one claim, one value
# ═════════════════════════════════════════════════════════════════════════════

#: Spellings of ONE claim (the same number printed to the same last digit).
_ONE_CLAIM = (
    ("exponent-zero", "1.0E0", "1.0"),
    ("leading-plus", "+1.0", "1.0"),
    ("scaled-mantissa", "10E-1", "1.0"),
    ("exponent-form", "1.234E-03", "0.001234"),
    ("lowercase-e", "1.234e-3", "0.001234"),
    ("negative-zero", "-0.000", "0.000"),
    ("surrounding-space", " 0.6123 ", "0.6123"),
)
#: Two claims each (a precision, a last-digit place, a sign).
_TWO_CLAIMS = (
    ("trailing-zero", "1.0", "1.00"),
    ("ten-and-1E+1", "10", "1E+1"),
    ("sign", "1.0", "-1.0"),
    ("last-digit", "0.99996", "0.99997"),
)


@pytest.mark.parametrize("a, b", [r[1:] for r in _ONE_CLAIM], ids=[r[0] for r in _ONE_CLAIM])
def test_r1_10_spellings_of_one_claim_are_one_value(a: str, b: str) -> None:
    """Equal both ways, one hash, one digest, one member of a set."""
    pa, pb = _P(a), _P(b)
    require(pa == pb and pb == pa, f"{a!r} and {b!r} print one claim, got {pa.text!r} and {pb.text!r}")
    require(hash(pa) == hash(pb) and len({pa, pb}) == 1, f"{a!r} and {b!r}: hash or set separates them")
    require(content_digest(pa) == content_digest(pb), f"{a!r} and {b!r}: two digests")


@pytest.mark.parametrize("a, b", [r[1:] for r in _TWO_CLAIMS], ids=[r[0] for r in _TWO_CLAIMS])
def test_r1_10_two_claims_are_two_values(a: str, b: str) -> None:
    """A trailing zero IS a claim (``1.0`` and ``1.00`` promise different
    precision), so the text is never reduced to its number."""
    pa, pb = _P(a), _P(b)
    require(pa != pb and len({pa, pb}) == 2, f"{a!r} and {b!r} are two claims, merged as {pa.text!r}")
    require(content_digest(pa) != content_digest(pb), f"{a!r} and {b!r}: one digest")


# ═════════════════════════════════════════════════════════════════════════════
# R1.11 — the value; R1.12 — the half unit (REFERENCE: a hand-written table)
# ═════════════════════════════════════════════════════════════════════════════

#: (printed text, its half unit in the last printed digit), written by hand
#: from the digits: the independent route to the exponent arithmetic.
_HALF_UNITS = (
    ("1.00000", "0.000005"),
    ("0.99996", "0.000005"),
    ("0.6123", "0.00005"),
    ("1.234E-03", "0.0000005"),
    ("12", "0.5"),
    ("1E+1", "5"),
    ("0", "0.5"),
    ("-2.5", "0.05"),
    ("6.02214076E+23", "5E+14"),
)


@pytest.mark.parametrize("text", [t for t, _ in _HALF_UNITS])
def test_r1_11_the_value_is_the_correctly_rounded_double(text: str) -> None:
    """``value`` is a ``float`` equal to Python's correctly rounded parse of
    the text (``float(str)`` rounds to nearest)."""
    p = _P(text)
    require(type(p.value) is float, f"{text!r}: value is a {type(p.value).__name__}")
    require(p.value == float(text), f"{text!r}: value {p.value!r} != {float(text)!r}")


@pytest.mark.parametrize("text, half", _HALF_UNITS, ids=[t for t, _ in _HALF_UNITS])
def test_r1_12_the_half_unit_is_exact(text: str, half: str) -> None:
    """``half_unit`` is a ``Decimal`` equal, exactly, to half a unit in the
    last printed digit: the author's precision claim, not recomputed."""
    p = _P(text)
    require(isinstance(p.half_unit, Decimal), f"{text!r}: half_unit is a {type(p.half_unit).__name__}")
    require(p.half_unit == Decimal(half), f"{text!r}: half_unit {p.half_unit} != {half}")


# ═════════════════════════════════════════════════════════════════════════════
# R1.13 — refusals of the text; R1.14 — the citation carries a place
# ═════════════════════════════════════════════════════════════════════════════

_BAD_TEXT = (
    ("empty", "", ValueError, "not a printed decimal number"),
    ("comma", "1,0", ValueError, "not a printed decimal number"),
    ("word", "abc", ValueError, "not a printed decimal number"),
    ("two-points", "1.0.0", ValueError, "not a printed decimal number"),
    ("hex-float", "0x1p3", ValueError, "not a printed decimal number"),
    ("nan", "nan", ValueError, "not a finite printed number"),
    ("NaN", "NaN", ValueError, "not a finite printed number"),
    ("sNaN", "sNaN", ValueError, "not a finite printed number"),
    ("inf", "inf", ValueError, "not a finite printed number"),
    ("Infinity", "-Infinity", ValueError, "not a finite printed number"),
    ("overflow", "1e400", ValueError, "beyond the range of a double"),
    ("underflow", "1e-400", ValueError, "beyond the range of a double"),
    ("grouping-underscore", "1_0.0", ValueError, "not a printed decimal number"),
    ("a-float", 1.0, TypeError, "a printed value is its decimal text"),
    ("a-decimal", Decimal("1.0"), TypeError, "a printed value is its decimal text"),
    ("none", None, TypeError, "a printed value is its decimal text"),
)


@pytest.mark.parametrize("text, error, fragment", [r[1:] for r in _BAD_TEXT], ids=[r[0] for r in _BAD_TEXT])
def test_r1_13_a_non_decimal_text_is_refused(text: Any, error: type[BaseException], fragment: str) -> None:
    with pytest.raises(error, match=fragment):
        _P(text)


def test_r1_14_the_citation_names_a_place() -> None:
    """A printed value is printed at a place: a citation with no locator is
    refused (``ValueError`` naming the locator), a bare key string is not a
    ``Citation`` (``TypeError``); the positive leg keeps the citation."""
    with pytest.raises(ValueError, match="locator"):
        _P("1.0", Citation("SoodForsterParsons2003"))
    with pytest.raises(TypeError, match="Citation"):
        _P("1.0", "SoodForsterParsons2003")
    require(_P("1.0").citation == _CITE, "the citation is kept")


# ═════════════════════════════════════════════════════════════════════════════
# R1.15 — one source: the fields are the text and the citation
# ═════════════════════════════════════════════════════════════════════════════


def test_r1_15_the_text_is_the_one_source() -> None:
    """``dataclasses.fields(Printed)`` is exactly ``(text, citation)``: no
    stored ``value``, ``digits`` or ``half_unit`` beside the text, so the
    value and its precision cannot disagree (X4)."""
    names = tuple(f.name for f in dataclasses.fields(_reading().Printed))
    require(names == ("text", "citation"), f"fields {names}")


_PRINTED_BASE = "0.99996"


def _printed(text: str = _PRINTED_BASE, citation: Citation = _CITE) -> Any:
    return _P(text, citation)


ROSTER: tuple[Entry, ...] = (
    Entry(
        cls=importlib.import_module(_READING).Printed,
        base=_printed,
        parts=("text", "citation"),
        perturb={
            "text": (
                leg("last digit", lambda: _printed("0.99997")),
                leg("one more digit", lambda: _printed("0.999960")),
            ),
            "citation": (
                leg("another place", lambda: _printed(citation=Citation("SoodForsterParsons2003", "Table 11"))),
                leg("another work", lambda: _printed(citation=Citation("SoodLA13511_1999", "Table 10"))),
            ),
        },
        pairs=(("plus-sign", lambda: _printed("+0.99996"), lambda: _printed("0.99996")),),
    ),
)


@pytest.mark.parametrize("entry", ROSTER, ids=[e.id for e in ROSTER])
def test_r1_15_the_population_is_the_parts(entry: Entry) -> None:
    check_population(entry)


@pytest.mark.parametrize(
    "entry, part, the_leg",
    perturbation_ids(ROSTER),
    ids=[param_id(e.id, p, lg[0]) for e, p, lg in perturbation_ids(ROSTER)],
)
def test_r1_15_each_part_moves_the_digest(entry: Entry, part: str, the_leg: Any) -> None:
    check_perturbation(entry, part, the_leg)


@pytest.mark.parametrize("entry, pair", pair_ids(ROSTER), ids=[param_id(e.id, p[0]) for e, p in pair_ids(ROSTER)])
def test_r1_15_equal_content_is_one_value(entry: Entry, pair: Any) -> None:
    check_equal_pair(entry, pair)


@pytest.mark.parametrize("entry", ROSTER, ids=[e.id for e in ROSTER])
def test_r1_15_pickle(entry: Entry) -> None:
    check_pickle(entry)


# ═════════════════════════════════════════════════════════════════════════════
# R1.16 — the two sums are disjoint; R1.17 — Estimated is not built
# ═════════════════════════════════════════════════════════════════════════════


def _members(alias: Any) -> set[type]:
    """The classes an alias names: a union's arguments, or the one class."""
    args = typing.get_args(alias)
    return set(args) if args else {alias}


@pytest.mark.rests_on(f"{_ENC}::test_r1_6_the_fields_are_the_value_and_the_bound")
def test_r1_16_the_reading_sums_are_disjoint_by_role() -> None:
    """``ReferenceReading`` names exactly ``{Enclosure, Printed}``;
    ``ProductionReading`` names exactly ``{Measured}``, and that ``Measured``
    IS ``orpheus.numerics.outcome.Measured`` (one definition of a measured
    number); the two sets share no class, and no member of either is a
    subclass of a member of the other, so a production reading cannot be
    passed where a reference reading is required. Instance legs: an
    ``Enclosure`` and a ``Printed`` are not production readings, a
    ``Measured`` is not a reference reading."""
    from orpheus.numerics import outcome
    from orpheus.numerics.enclosure import Enclosure

    reading = _reading()
    ref, prod = _members(reading.ReferenceReading), _members(outcome.ProductionReading)
    require(ref == {Enclosure, reading.Printed}, f"ReferenceReading names {sorted(c.__qualname__ for c in ref)}")
    require(prod == {outcome.Measured}, f"ProductionReading names {sorted(c.__qualname__ for c in prod)}")
    crossed = [(r.__qualname__, p.__qualname__) for r in ref for p in prod if issubclass(r, p) or issubclass(p, r)]
    require(not crossed, f"a reference reading and a production reading are related by subclassing: {crossed}")
    require(not isinstance(Enclosure(1.0, 0.0), outcome.ProductionReading), "an Enclosure is a production reading")
    require(not isinstance(_P("1.0"), outcome.ProductionReading), "a Printed is a production reading")
    require(not isinstance(outcome.Measured(1.0), reading.ReferenceReading), "a Measured is a reference reading")


def _class_definitions() -> tuple[int, dict[str, list[str]]]:
    """Every class definition under ``orpheus/``, by name, with its files."""
    files = sorted((_ROOT / "orpheus").rglob("*.py"))
    found: dict[str, list[str]] = {}
    for path in files:
        for node in ast.walk(ast.parse(path.read_text(), filename=str(path))):
            if isinstance(node, ast.ClassDef):
                found.setdefault(node.name, []).append(str(path.relative_to(_ROOT)))
    return len(files), found


def test_r1_17_estimated_is_not_built_and_measured_is_defined_once() -> None:
    """By AST over every ``.py`` under ``orpheus/`` (input count printed, more
    than 300 asserted): no class ``Estimated`` is defined (the user's G3
    ruling: it waits for a Monte Carlo consumer, so this row is edited ON
    PURPOSE when it lands), and ``Measured`` is defined exactly once, in
    ``orpheus/numerics/outcome.py`` (positive control: the census finds it)."""
    n_files, classes = _class_definitions()
    print(f"R1.17: {n_files} files parsed")
    require(n_files > 300, f"activation: only {n_files} files parsed")
    require(classes.get("Measured") == ["orpheus/numerics/outcome.py"], f"Measured defined at {classes.get('Measured')}")
    require("Estimated" not in classes, f"Estimated is defined at {classes.get('Estimated')}: edit R1.16 and R1.17 with it")


# ═════════════════════════════════════════════════════════════════════════════
# R1.18 — the reference package's layer
# ═════════════════════════════════════════════════════════════════════════════

#: The packages ``orpheus.reference`` may import (the orchestrator's Q4 ruling).
_REFERENCE_MAY_IMPORT = {"reference", "specification", "numerics", "data", "geometry"}
#: The packages it must not import, each of which the linter must list.
_REFERENCE_MUST_NOT = {"mesh", "transport", "derivations", "sn", "diffusion", "homogeneous", "cp", "moc", "mc"}
#: The packages below it, which must not import it.
_BELOW_REFERENCE = ("numerics", "geometry", "data", "mesh", "specification")


def test_r1_18_the_package_is_registered_with_the_linter() -> None:
    """``tests/gates/test_layer_imports.py`` knows the package: its forbidden
    set holds every package in ``_REFERENCE_MUST_NOT``; every package below it
    is forbidden to import it; ``derivations`` (which writes reference
    solutions) is not; ``orpheus.reference`` is a cold entry point. An
    unregistered package is SKIPPED by the linter, so this row is what makes
    its zero violations mean something."""
    from tests.gates import test_layer_imports as lint

    edges = lint.FORBIDDEN_EDGES
    require("reference" in edges, "the linter does not know the package: its imports are never checked")
    missing = sorted(_REFERENCE_MUST_NOT - edges["reference"])
    require(not missing, f"the linter lets the reference package import {missing}")
    unguarded = [p for p in _BELOW_REFERENCE if "reference" not in edges[p]]
    require(not unguarded, f"packages below the reference package may import it: {unguarded}")
    require("reference" not in edges["derivations"], "derivations may not import the reference package")
    marks = getattr(lint.test_entry_point_imports_in_a_fresh_interpreter, "pytestmark")
    entries = next(m for m in marks if m.name == "parametrize").args[1]
    require("orpheus.reference" in entries, "orpheus.reference is not a cold entry point")


_LAYER_SCRIPT = """
import sys
import orpheus.reference.reading as m
from orpheus.data.citation import Citation
m.Printed('1.0', Citation('SoodForsterParsons2003', 'Table 10'))
print(m.__file__)
print(sorted({k.split('.')[1] for k in sys.modules if k.startswith('orpheus.')}))
"""


def test_r1_18_a_cold_import_loads_only_the_layers_below() -> None:
    """A fresh interpreter importing the module and building one ``Printed``
    loads only packages in ``_REFERENCE_MAY_IMPORT`` (``[M]`` 2026-10-03 on
    the prototype: data, geometry, numerics, reference)."""
    env = {**os.environ, "PYTHONPATH": str(_ROOT)}
    out = subprocess.run([sys.executable, "-O", "-c", _LAYER_SCRIPT], cwd=_ROOT, env=env, capture_output=True, text=True)
    require(out.returncode == 0, out.stderr)
    file, packages = out.stdout.strip().splitlines()
    require(file.startswith(str(_ROOT)), f"the subprocess imported {file} (L22)")
    loaded = set(ast.literal_eval(packages))
    require("reference" in loaded, f"activation: {loaded}")
    require(loaded <= _REFERENCE_MAY_IMPORT, f"loads {sorted(loaded - _REFERENCE_MAY_IMPORT)}")


# ═════════════════════════════════════════════════════════════════════════════
# R1.19 — the digests in every process
# ═════════════════════════════════════════════════════════════════════════════

_SEED_SCRIPT = """
from orpheus.numerics.enclosure import Enclosure
from orpheus.reference.reading import Printed
from orpheus.data.citation import Citation
from orpheus.numerics.content import content_digest
import orpheus.reference.reading as m
print(m.__file__)
for v in (Enclosure(1.25, 0.5), Enclosure(1.0, 0.0) / Enclosure(3.0, 0.0),
          Printed('0.99996', Citation('SoodForsterParsons2003', 'Table 10'))):
    print(content_digest(v).hex(), hash(v))
print('control', hash('salted'))
"""


def _seed_run(seed: str) -> list[str]:
    env = {**os.environ, "PYTHONHASHSEED": seed, "PYTHONPATH": str(_ROOT)}
    out = subprocess.run([sys.executable, "-O", "-c", _SEED_SCRIPT], cwd=_ROOT, env=env, capture_output=True, text=True)
    require(out.returncode == 0, out.stderr)
    return out.stdout.strip().splitlines()


@pytest.mark.rests_on("tests/gates/numerics/test_content_identity.py::test_s5_1_digests_and_hashes_are_seed_stable")
def test_r1_19_digests_and_hashes_are_seed_stable() -> None:
    """``PYTHONHASHSEED`` 1 and 2 print the same three digests and hashes (a
    ``str`` part, the printed text, is where a salted hash would leak); the
    control line, a ``str`` hash, must differ (the seeds took effect)."""
    one, two = _seed_run("1"), _seed_run("2")
    require(one[0].startswith(str(_ROOT)), f"the subprocess imported {one[0]} (L22)")
    require(one[1:-1] == two[1:-1], f"digests or hashes moved with the seed:\n{one}\n{two}")
    require(one[-1] != two[-1], "control: the str hash did not move, so the seeds did not take effect")


def test_r1_12_the_half_unit_ignores_the_callers_decimal_context() -> None:
    """The half unit is exact under ANY decimal context: a caller's narrowed
    context (here ``Emin=-99``, ``Emax=99``) must neither flush ``1.00E-200``'s
    half unit to zero nor overflow on ``1.00E+150`` (qa, 2026-10-03: the
    context-dependent ``Decimal(5).scaleb`` returned ``0E-126`` and raised)."""
    import decimal

    with decimal.localcontext() as context:
        context.Emin, context.Emax = -99, 99
        small, large = _P("1.00E-200").half_unit, _P("1.00E+150").half_unit
    require(small == Decimal("5E-203"), f"half unit of 1.00E-200 under a narrowed context: {small}")
    require(large == Decimal("5E+147"), f"half unit of 1.00E+150 under a narrowed context: {large}")
