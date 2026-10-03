r"""The published solution, the standing, the one Withdrawal and the printed enclosure (#405 P2 step 5, R5w, R5e, R5p).

Specified by the test-architect (2026-10-03, ``.claude/plans/reference_p2_spec.md``
§1.5), on the user's rulings of 2026-10-03 and the orchestrator's API ruling:

* ``Withdrawal(reason, issue)`` has ONE home, ``orpheus/reference/withdrawal.py``;
  P0's lock in ``derivations/common/withdrawal.py`` imports it; no re-export
  shim; every module that binds the name binds that one object;
* ``Standing = Current | Withdrawal``, one sum for a publication and a
  certificate;
* ``PublishedSolution(specification, printed, standing)``: each printed value
  carries its own citation (``Printed``); ``read`` returns the printed value,
  refuses an unprinted observable (``NotPrinted``) and refuses every read of a
  withdrawn publication; construction refuses an observable the
  specification cannot pose, through the ONE admission function
  ``admit_observable(observable, specification)`` the step-6 reader reuses;
* ``Printed.enclosure()``: the printed claim as an ``Enclosure``, rounded
  outward (the review note C4: one verb, never an ``isinstance`` branch per
  consumer).
"""

from __future__ import annotations

import ast
import math
import os
import subprocess
import sys
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from typing import Any

import pytest

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
from tests.gates.reference import _step5 as s5

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/reference/test_published.py"
_ROOT = Path(__file__).resolve().parents[3]


# ═════════════════════════════════════════════════════════════════════════════
# R5w — the one Withdrawal
# ═════════════════════════════════════════════════════════════════════════════


def test_r5w_1_withdrawal_is_defined_once_in_the_reference_package() -> None:
    """By AST over every ``.py`` under ``orpheus/`` and ``tests/`` (input count
    printed, more than 1000 asserted): exactly one class ``Withdrawal`` is
    defined, in ``orpheus/reference/withdrawal.py``; positive control: the
    census finds ``Printed`` in ``orpheus/reference/reading.py``."""
    files = sorted((_ROOT / "orpheus").rglob("*.py")) + sorted((_ROOT / "tests").rglob("*.py"))
    print(f"R5w.1: {len(files)} files parsed")
    require(len(files) > 1000, f"activation: {len(files)} files parsed")
    defined: dict[str, list[str]] = {}
    for path in files:
        for node in ast.walk(ast.parse(path.read_text(), filename=str(path))):
            if isinstance(node, ast.ClassDef) and node.name in ("Withdrawal", "Printed"):
                defined.setdefault(node.name, []).append(str(path.relative_to(_ROOT)))
    require(defined.get("Printed") == ["orpheus/reference/reading.py"], f"control: Printed at {defined.get('Printed')}")
    require(defined.get("Withdrawal") == ["orpheus/reference/withdrawal.py"], f"Withdrawal defined at {defined.get('Withdrawal')}")


_BINDERS = (
    "orpheus.derivations.common.withdrawal",
    "orpheus.derivations.continuous.peierls_nystrom",
    "tests._harness.registry",
    "tests.conftest",
)


@pytest.mark.rests_on(f"{_HERE}::test_r5w_1_withdrawal_is_defined_once_in_the_reference_package")
def test_r5w_2_every_binding_is_the_one_object() -> None:
    """Every module that binds the name ``Withdrawal`` (the lock, the
    withdrawn Peierls family, the test registry, the root conftest) binds the
    reference package's class, by identity: one type, so a reference
    certificate's ``Withdrawal`` standing and P0's lock agree by construction.
    The lock's own machinery (``withdrawn_generator``, ``GeneratorWithdrawn``)
    stays in derivations."""
    import importlib

    one = s5.c(s5.WITHDRAWAL, "Withdrawal")
    for name in _BINDERS:
        bound = getattr(importlib.import_module(name), "Withdrawal", None)
        require(bound is one, f"{name}.Withdrawal is {bound!r}, not {one!r}")
    lock = importlib.import_module("orpheus.derivations.common.withdrawal")
    require(hasattr(lock, "withdrawn_generator") and hasattr(lock, "GeneratorWithdrawn"), "the lock moved with the type")


def test_r5w_3_the_law_is_unchanged() -> None:
    """Positive and negative legs of the moved law (its full gate is
    ``tests/gates/test_withdrawal.py``, re-pointed): a reason and a positive
    issue construct; an empty reason, a ``bool`` issue and a zero issue are
    refused naming the field."""
    W = s5.c(s5.WITHDRAWAL, "Withdrawal")
    w = W("a reason", 12)
    require((w.reason, w.issue) == ("a reason", 12), f"{w!r}")
    for args, fragment in ((("", 1), "reason"), (("r", True), "issue"), (("r", 0), "issue")):
        with pytest.raises(ValueError, match=fragment):
            W(*args)


# ═════════════════════════════════════════════════════════════════════════════
# R5e — the printed claim as an enclosure
# ═════════════════════════════════════════════════════════════════════════════

_TEXTS = ("1.00000", "0.99996", "0.6123", "1.234E-03", "12", "1E+1", "0", "-2.5", "6.02214076E+23", "0.1")


@pytest.mark.parametrize("text", _TEXTS)
def test_r5e_the_printed_interval_lies_inside_its_enclosure(text: str) -> None:
    r"""The printed interval :math:`[t - h, t + h]` (the text ``t`` and its half
    unit ``h``, in EXACT decimal arithmetic) lies inside ``enclosure()``'s
    interval (exact rationals on its two doubles), and the enclosure is tight:
    each end exceeds the printed one by less than 2 ulp of it (one rounding,
    one outward step). The centre is the printed value's double."""
    p = s5.printed(text)
    e = p.enclosure()
    require(type(e) is s5.c(s5.ENCLOSURE, "Enclosure"), f"enclosure() returned a {type(e).__name__}")
    t, h = Fraction(Decimal(p.text)), Fraction(p.half_unit)
    lo, hi = Fraction(e.value) - Fraction(e.bound), Fraction(e.value) + Fraction(e.bound)
    require(lo <= t - h and t + h <= hi, f"{text}: [{float(t - h)}, {float(t + h)}] not inside [{float(lo)}, {float(hi)}]")
    for end, edge in ((lo, t - h), (hi, t + h)):
        slack = 2 * Fraction(math.ulp(float(edge))) + 2 * Fraction(math.ulp(e.value))
        require(abs(end - edge) < slack, f"{text}: an end is loose by {float(abs(end - edge)):.3e}")
    require(e.value == p.value, f"{text}: centre {e.value!r} != printed value {p.value!r}")


# ═════════════════════════════════════════════════════════════════════════════
# R5p — the published solution
# ═════════════════════════════════════════════════════════════════════════════


def test_r5p_1_read_returns_the_printed_value_and_refuses_an_unprinted_one() -> None:
    """``read`` returns the very ``Printed`` stored for the observable; a second,
    equal-content spelling of the observable reads the same value; an
    observable the publication does not print raises ``NotPrinted`` (a
    ``LookupError``) matching "does not print"."""
    pub = s5.published()
    require(pub.read(s5.eigenvalue()) == s5.printed("1.00000"), "the eigenvalue reads another value")
    require(pub.read(s5.flux()) == s5.printed("0.6123", "Table 11"), "an equal-content flux integral reads another value")
    not_printed = s5.c(s5.PUBLISHED, "NotPrinted")
    require(issubclass(not_printed, LookupError), "NotPrinted is not a LookupError")
    with pytest.raises(not_printed, match="does not print"):
        pub.read(s5.flux(((2.0, 0.5), (0.0, 2.0))))
    with pytest.raises(not_printed, match="does not print"):
        pub.read(s5.point(0.5, 1))


def test_r5p_2_a_withdrawn_publication_refuses_every_read() -> None:
    """Built with a ``Withdrawal`` standing (an erratum), every read refuses,
    naming the issue; the same publication ``Current`` reads (the positive
    leg)."""
    pub = s5.published(standing=s5.withdrawal(issue=4242))
    with pytest.raises(Exception, match="4242"):
        pub.read(s5.eigenvalue())
    s5.published().read(s5.eigenvalue())


_BAD = (
    ("specification-not-a-specification", lambda: s5.published(spec="slab"), TypeError, "specification"),
    ("prints-nothing", lambda: s5.published(printed_map={}), ValueError, "prints nothing"),
    ("key-not-an-observable", lambda: s5.published(printed_map={"k": s5.printed("1.0")}), TypeError, "observable"),
    ("value-not-printed", lambda: s5.published(printed_map={s5.eigenvalue(): 1.0}), TypeError, "Printed"),
    ("value-an-enclosure", lambda: s5.published(printed_map={s5.eigenvalue(): s5.enclosure(1.0, 0.0)}), TypeError, "Printed"),
    ("standing-a-string", lambda: s5.published(standing="current"), TypeError, "standing"),
)


@pytest.mark.parametrize("build, error, fragment", [b[1:] for b in _BAD], ids=[b[0] for b in _BAD])
def test_r5p_3_admission(build: Any, error: type[BaseException], fragment: str) -> None:
    with pytest.raises(error, match=fragment):
        build()


_UNPOSABLE = (
    ("eigenvalue-of-a-source-question", lambda: s5.published(spec=s5.source_slab(), printed_map={s5.eigenvalue(): s5.printed("1.0")})),
    ("point-value-on-the-infinite-medium", lambda: s5.published(spec=s5.eigen_medium(), printed_map={s5.point(0.5, 0): s5.printed("1.0")})),
    ("weight-with-three-regions", lambda: s5.published(printed_map={s5.flux(((1.0, 0.0), (1.0, 0.0), (1.0, 0.0))): s5.printed("1.0")})),
    ("weight-with-one-group", lambda: s5.published(printed_map={s5.flux(((1.0,), (1.0,))): s5.printed("1.0")})),
    ("point-value-group-out-of-range", lambda: s5.published(printed_map={s5.point(0.5, 2): s5.printed("1.0")})),
    ("point-value-outside-the-geometry", lambda: s5.published(printed_map={s5.point(3.5, 0): s5.printed("1.0")})),
)


@pytest.mark.parametrize("build", [b for _, b in _UNPOSABLE], ids=[n for n, _ in _UNPOSABLE])
def test_r5p_4_an_unposable_observable_is_refused(build: Any) -> None:
    """A publication cannot print what its specification cannot pose. The
    refusal is ``admit_observable``'s (rows R5p.5 prove the route)."""
    with pytest.raises((TypeError, ValueError)):
        build()


def test_r5p_4_the_posable_observables_construct() -> None:
    """Positive legs: the eigenvalue on an eigen slab and on the infinite
    medium; a group-0 point inside the slab; a fitting weight on a source
    slab."""
    s5.published()
    s5.published(spec=s5.eigen_medium(), printed_map={s5.eigenvalue(): s5.printed("1.40")})
    s5.published(printed_map={s5.point(2.5, 1): s5.printed("0.25")})
    s5.published(spec=s5.source_slab(), printed_map={s5.flux(): s5.printed("3.1")})


@pytest.mark.rests_on(f"{_HERE}::test_r5p_4_the_posable_observables_construct")
def test_r5p_5_the_admission_is_the_specifications_one_function(monkeypatch: pytest.MonkeyPatch) -> None:
    """ROUTE gate (X4: one admission). Rebind ``admit_observable`` in every
    module that binds it to a decoy that raises a marked error: building a
    POSABLE publication must now raise the decoy's error (the publication
    calls the function), and an unposable one too (no second admission
    refuses first). Activation: the decoy's call count is positive."""
    import importlib

    module_name, name = s5.ADMISSION
    honest = getattr(importlib.import_module(module_name), name)
    calls: list[Any] = []

    def decoy(*args: Any, **kwargs: Any) -> None:
        calls.append(args)
        raise ValueError("DECOY-ADMISSION")

    rebound = 0
    for mod in list(sys.modules.values()):
        if mod is not None and getattr(mod, name, None) is honest:
            monkeypatch.setattr(mod, name, decoy)
            rebound += 1
    require(rebound >= 1, "activation: nothing bound the admission function")
    with pytest.raises(ValueError, match="DECOY-ADMISSION"):
        s5.published()
    with pytest.raises(ValueError, match="DECOY-ADMISSION"):
        s5.published(spec=s5.eigen_medium(), printed_map={s5.point(0.5, 0): s5.printed("1.0")})
    require(len(calls) >= 2, f"activation: the decoy ran {len(calls)} times")


# ═════════════════════════════════════════════════════════════════════════════
# R5p.6 — content identity (Current and PublishedSolution)
# ═════════════════════════════════════════════════════════════════════════════

ROSTER: tuple[Entry, ...] = (
    Entry(cls=s5.c("orpheus.reference.withdrawal", "Withdrawal"), base=s5.withdrawal, parts=("reason", "issue"),
          perturb={"reason": (leg("another reason", lambda: s5.c("orpheus.reference.withdrawal", "Withdrawal")("another reason", 999)),),
                   "issue": (leg("another issue", lambda: s5.c("orpheus.reference.withdrawal", "Withdrawal")("the publication was found unfit as a whole", 1000)),)},
          pairs=(("two builds", s5.withdrawal, s5.withdrawal),)),
    Entry(cls=s5.c(s5.PUBLISHED, "Current"), base=s5.current, parts=(), perturb={},
          pairs=(("two builds", s5.current, s5.current),)),
    Entry(
        cls=s5.c(s5.PUBLISHED, "PublishedSolution"),
        base=s5.published,
        parts=("specification", "printed", "standing"),
        perturb={
            "specification": (leg("another material in region 1", lambda: s5.published(spec=s5.eigen_slab(s5.sf.n2n_only()))),),
            "printed": (
                leg("a last printed digit", lambda: s5.published(printed_map={s5.eigenvalue(): s5.printed("1.00001"), s5.flux(): s5.printed("0.6123", "Table 11")})),
                leg("another citation locator", lambda: s5.published(printed_map={s5.eigenvalue(): s5.printed("1.00000", "Table 12"), s5.flux(): s5.printed("0.6123", "Table 11")})),
                leg("another observable", lambda: s5.published(printed_map={s5.eigenvalue(): s5.printed("1.00000"), s5.point(0.5, 0): s5.printed("0.6123", "Table 11")})),
            ),
            "standing": (leg("withdrawn", lambda: s5.published(standing=s5.withdrawal())),),
        },
        pairs=(("the printed mapping's insertion order", s5.published,
                lambda: s5.published(printed_map={s5.flux(): s5.printed("0.6123", "Table 11"), s5.eigenvalue(): s5.printed("1.00000")})),),
    ),
)


@pytest.mark.parametrize("entry", ROSTER, ids=[e.id for e in ROSTER])
def test_r5p_6_the_population_is_the_parts(entry: Entry) -> None:
    check_population(entry)


@pytest.mark.parametrize("entry, part, the_leg", perturbation_ids(ROSTER),
                         ids=[param_id(e.id, p, lg[0]) for e, p, lg in perturbation_ids(ROSTER)])
def test_r5p_6_each_part_moves_the_digest(entry: Entry, part: str, the_leg: Any) -> None:
    check_perturbation(entry, part, the_leg)


@pytest.mark.parametrize("entry, pair", pair_ids(ROSTER), ids=[param_id(e.id, p[0]) for e, p in pair_ids(ROSTER)])
def test_r5p_6_equal_content_is_one_value(entry: Entry, pair: Any) -> None:
    check_equal_pair(entry, pair)


@pytest.mark.parametrize("entry", ROSTER, ids=[e.id for e in ROSTER])
def test_r5p_6_pickle(entry: Entry) -> None:
    check_pickle(entry)


def test_r5p_6_a_withdrawn_and_a_current_publication_are_two_values() -> None:
    a, b = s5.published(), s5.published(standing=s5.withdrawal())
    require(a != b and content_digest(a) != content_digest(b) and len({a, b}) == 2, "the standing is not content")


# ═════════════════════════════════════════════════════════════════════════════
# R5.12 — the layer
# ═════════════════════════════════════════════════════════════════════════════

_LAYER_SCRIPT = """
import sys
import orpheus.reference.published, orpheus.reference.certificate, orpheus.reference.withdrawal
print(orpheus.reference.published.__file__)
print(sorted({k.split('.')[1] for k in sys.modules if k.startswith('orpheus.')}))
"""


def test_r5_12_a_cold_import_loads_only_the_layers_below() -> None:
    """Importing the three step-5 modules loads only
    ``{reference, specification, numerics, data, geometry}``: no derivations
    (the reference package may not import it; P0's lock imports the
    reference package, the other way)."""
    env = {**os.environ, "PYTHONPATH": str(_ROOT)}
    out = subprocess.run([sys.executable, "-O", "-c", _LAYER_SCRIPT], cwd=_ROOT, env=env, capture_output=True, text=True)
    require(out.returncode == 0, out.stderr)
    file, packages = out.stdout.strip().splitlines()
    require(file.startswith(str(_ROOT)), f"the subprocess imported {file} (L22)")
    loaded = set(ast.literal_eval(packages))
    require("reference" in loaded, f"activation: {loaded}")
    allowed = {"reference", "specification", "numerics", "data", "geometry"}
    require(loaded <= allowed, f"loads {sorted(loaded - allowed)}")
