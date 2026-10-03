r"""The reference solution and its reading (#405 P2 step 6, R6.1-R6.9).

Specified by the test-architect (2026-10-03, ``.claude/plans/reference_p2_spec.md``
§1.6), on the user's rulings of 2026-10-03 (a reference's reading is its
natural extension, one answer; no stored representation in P2; certification
decoupled from caching) and qa F2 of step 5:

* ``ReferenceSolution(specification, derivation, certificate)``: the
  derivation (a ``Derivation``, ``establish(observable) -> Exact |
  DerivedBound``: the natural extension evaluated WITH its derived bound, on
  demand, for any admissible observable) is a plain reference; the
  certificate is a ``ReferenceCertificate`` or ``None``;
* construction enforces what a certificate cannot: every claimed observable is
  posable on the specification (``admit_observable``), every anchor prints for
  the SAME specification, and no claim is on a ``Ratio`` (a ratio's reading is
  the quotient of its operands' readings, one definition);
* ``read(observable)``: ``admit_observable`` first; a ratio by
  ``Enclosure.__truediv__`` of its operands' readings; an eigenvalue or a
  linear observable is ESTABLISHED by the derivation and its enclosure
  returned (a derived enclosure is a guarantee on its own, G3); when the
  certificate claims the observable, the claim's enclosure and the
  established one must share a point (refused, "disagrees", otherwise);
* an observable the derivation cannot bound raises ``NotCertified`` from
  ``establish``; what an uncertified reading returns instead is the user's open
  question, and that row is marked to be re-posed.
"""

from __future__ import annotations

import ast
import math
import os
import subprocess
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import pytest

from tests.gates._content_identity_helpers import require
from tests.gates.reference import _step5 as s5
from tests.gates.reference import _step6 as s6
from tests.gates.specification import _fixtures as sf

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/reference/test_reference_solution.py"
_ROOT = Path(__file__).resolve().parents[3]


def _contains(enclosure: Any, exact: Fraction) -> bool:
    return abs(Fraction(enclosure.value) - exact) <= Fraction(enclosure.bound)


# ═════════════════════════════════════════════════════════════════════════════
# R6.1 — construction
# ═════════════════════════════════════════════════════════════════════════════


def test_r6_1_the_fields_and_their_types() -> None:
    """The three fields, all required (no default on the certificate: a family
    without a derived bound passes ``None`` and says so at the call); the
    test-side table derivation satisfies the ``Derivation`` protocol; a
    certificate of ``None`` constructs."""
    import dataclasses
    import inspect

    cls = s5.c(s6.SOLUTION, "ReferenceSolution")
    names = tuple(f.name for f in dataclasses.fields(cls))
    require(names == ("specification", "derivation", "certificate"), f"fields {names}")
    defaulted = [p.name for p in inspect.signature(cls).parameters.values() if p.default is not inspect.Parameter.empty]
    require(not defaulted, f"defaulted fields {defaulted}")
    spec, derivation, cert = s6.medium_table()
    require(isinstance(derivation, s6.derivation_protocol()), "the table derivation is not a Derivation")
    s6.reference(spec, derivation, cert)
    s6.reference(spec, derivation, None)


_BAD = (
    ("specification-a-string", lambda d, c: s6.reference("medium", d, c), TypeError, "specification"),
    ("derivation-without-evaluate", lambda d, c: s6.reference(s5.eigen_medium(), object(), c), TypeError, "derivation"),
    ("certificate-a-mapping", lambda d, c: s6.reference(s5.eigen_medium(), d, {}), TypeError, "certificate"),
)


@pytest.mark.parametrize("build, error, fragment", [b[1:] for b in _BAD], ids=[b[0] for b in _BAD])
def test_r6_1_admission(build: Any, error: type[BaseException], fragment: str) -> None:
    _, derivation, cert = s6.medium_table()
    with pytest.raises(error, match=fragment):
        build(derivation, cert)


def test_r6_1_a_claim_the_specification_cannot_pose_is_refused() -> None:
    """qa F2: a certificate claiming a ``PointValue`` cannot be held by an
    infinite-medium reference (the medium has no position); the refusal is
    ``admit_observable``'s (R6.2's route covers the function)."""
    _, derivation, _ = s6.medium_table()
    cert = s5.certificate({s5.point(0.5, 0): s6.exact_claim(Fraction(1, 2))})
    with pytest.raises(ValueError, match="position"):
        s6.reference(s5.eigen_medium(), derivation, cert)


def test_r6_1_an_anchor_of_another_specification_is_refused() -> None:
    """qa F2, ``[M]`` before this step: a certificate claiming k on the medium,
    anchored by a publication printing for a DIFFERENT specification (the
    slab), read ``Valid``. The holder refuses it, naming the specification;
    the same anchor printed for the holder's own specification is admitted."""
    spec, derivation, _ = s6.medium_table()
    foreign = s5.published(spec=s5.eigen_slab(), printed_map={s5.eigenvalue(): s5.printed("1.40000")})
    own = s5.published(spec=spec, printed_map={s5.eigenvalue(): s5.printed("1.40000")})
    claims = {s5.eigenvalue(): s6.exact_claim(Fraction(7, 5))}
    with pytest.raises(ValueError, match="specification"):
        s6.reference(spec, derivation, s5.certificate(claims, [s5.corroboration(foreign)]))
    s6.reference(spec, derivation, s5.certificate(claims, [s5.corroboration(own)]))


def test_r6_1_a_ratio_claim_is_refused() -> None:
    """PROPOSAL (§1.6): a ratio's reading is the quotient of its operands'
    readings, so a claim on a ``Ratio`` would be a second definition of one
    number; the holder refuses it, naming the quotient."""
    spec, derivation, _ = s6.medium_table()
    ratio = s5.c(s5.OBSERVABLE, "Ratio")(s6.group_flux(0), s6.group_flux(1))
    cert = s5.certificate({ratio: s6.exact_claim(Fraction(390, 77))})
    with pytest.raises(ValueError, match="quotient"):
        s6.reference(spec, derivation, cert)


# ═════════════════════════════════════════════════════════════════════════════
# R6.2 — read admits first, through the one function
# ═════════════════════════════════════════════════════════════════════════════

_UNPOSABLE = (
    ("point-value-on-the-medium", lambda: s5.point(0.5, 0)),
    ("weight-with-two-regions", lambda: s6.group_flux(0, regions=2)),
    ("weight-with-three-groups", lambda: s6.group_flux(0, groups=3)),
)


@pytest.mark.parametrize("make", [m for _, m in _UNPOSABLE], ids=[n for n, _ in _UNPOSABLE])
def test_r6_2_an_unposable_observable_is_refused_before_the_derivation_runs(make: Any) -> None:
    spec, derivation, cert = s6.medium_table()
    ref = s6.reference(spec, derivation, cert)
    derivation.calls.clear()
    with pytest.raises(ValueError):
        ref.read(make())
    require(not derivation.calls, f"the derivation ran on an unposable observable: {derivation.calls}")


def test_r6_2_read_admits_through_the_specifications_one_function(monkeypatch: pytest.MonkeyPatch) -> None:
    """ROUTE gate: with ``admit_observable`` rebound in every binding to a decoy,
    a posable read raises the decoy's error (the reader calls it) and so does
    an unposable one (no second admission in the reader refuses first)."""
    import importlib

    spec, derivation, cert = s6.medium_table()
    ref = s6.reference(spec, derivation, cert)
    module_name, name = s5.ADMISSION
    honest = getattr(importlib.import_module(module_name), name)

    def decoy(*args: Any, **kwargs: Any) -> None:
        raise ValueError("DECOY-ADMISSION")

    rebound = 0
    for mod in list(sys.modules.values()):
        if mod is not None and getattr(mod, name, None) is honest:
            monkeypatch.setattr(mod, name, decoy)
            rebound += 1
    require(rebound >= 1, "activation: nothing bound admit_observable")
    for observable in (s5.eigenvalue(), s5.point(0.5, 0)):
        with pytest.raises(ValueError, match="DECOY-ADMISSION"):
            ref.read(observable)


# ═════════════════════════════════════════════════════════════════════════════
# R6.3 — a read is the derivation's established enclosure, checked against a claim
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("which", ["eigenvalue", "group-0", "group-1"])
def test_r6_3_a_read_is_the_derivations_established_enclosure(which: str) -> None:
    """The reading equals ``derivation.establish(observable).enclosure()``, and
    the derivation was asked for exactly that observable, on demand."""
    spec, derivation, cert = s6.medium_table()
    ref = s6.reference(spec, derivation, cert)
    observable = {"eigenvalue": s5.eigenvalue(), "group-0": s6.group_flux(0), "group-1": s6.group_flux(1)}[which]
    derivation.calls.clear()
    reading = ref.read(observable)
    require(derivation.calls == [observable], f"{which}: the derivation was called with {derivation.calls}")
    require(reading == derivation.establish(observable).enclosure(), f"{which}: read {reading!r}")


def test_r6_3_the_reading_is_the_established_enclosure_not_the_claims() -> None:
    """When the claim (exact 7/5) and the established enclosure (a derived bound
    1.4 +- 1e-6) agree, the reading is the ESTABLISHED one: the certificate
    corroborates, it does not replace the derivation's answer."""
    spec, derivation, cert = s6.medium_table()
    derived = s5.derived(1.4, 1e-6, "a residual bound")
    derivation.table[s5.eigenvalue()] = derived
    reading = s6.reference(spec, derivation, cert).read(s5.eigenvalue())
    require(reading == derived.enclosure(), f"read {reading!r}, established {derived.enclosure()!r}")


def test_r6_3_an_unclaimed_observable_reads_its_established_enclosure() -> None:
    """G2's point: an observable the certificate never declared (a weight
    (1, 1), a test mesh's cell average in general) reads the derivation's
    enclosure, with and without a certificate."""
    spec, derivation, cert = s6.medium_table()
    both = s5.flux(((1.0, 1.0),))
    derivation.table[both] = Fraction(30, 7) + Fraction(11, 13)
    for certificate in (cert, None):
        reading = s6.reference(spec, derivation, certificate).read(both)
        require(_contains(reading, Fraction(30, 7) + Fraction(11, 13)), f"{reading!r}")


def test_r6_3_a_derivation_that_disagrees_with_its_claim_is_refused() -> None:
    """A derivation establishing k = 7/5 + 1/10^9 exactly, against an exact claim
    7/5: the two enclosures share no point (the step-5 family law), so the read
    is refused, naming the disagreement: the certificate is about another
    answer than the one the reference computes."""
    spec, derivation, cert = s6.medium_table()
    derivation.table[s5.eigenvalue()] = Fraction(7, 5) + Fraction(1, 10**9)
    ref = s6.reference(spec, derivation, cert)
    with pytest.raises(ValueError, match="disagree"):
        ref.read(s5.eigenvalue())


@pytest.mark.parametrize("standing", ["invalid", "withdrawn"])
def test_r6_3_a_non_valid_certificate_still_reads(standing: str) -> None:
    """Ruled 2026-10-03: ``read`` answers on an ``Invalid`` or ``Withdrawn``
    certificate (only the ``VerificationCertificate`` refuses it)."""
    spec, derivation, _ = s6.medium_table()
    if standing == "invalid":
        cert = s5.certificate({s5.eigenvalue(): s5.claim(1e-30, s6.exact_of(Fraction(7, 5)))})
        require(s5.state_kind(cert) == "Invalid", f"activation: {cert.state!r}")
    else:
        cert = s5.certificate({s5.eigenvalue(): s6.exact_claim(Fraction(7, 5))}, standing=s5.withdrawal())
    reading = s6.reference(spec, derivation, cert).read(s5.eigenvalue())
    require(_contains(reading, Fraction(7, 5)), f"{reading!r}")


# ═════════════════════════════════════════════════════════════════════════════
# R6.4 — a ratio is the quotient of its operands' readings
# ═════════════════════════════════════════════════════════════════════════════


def test_r6_4_a_ratio_reads_through_the_one_quotient(monkeypatch: pytest.MonkeyPatch) -> None:
    """ROUTE gate (X4): with ``Enclosure.__truediv__`` wrapped by a counter, a
    ratio's reading equals ``read(a) / read(b)`` and the counter ran once;
    the quotient's enclosure contains the exact ratio (30/7)/(11/13)."""
    Enclosure = s5.c(s5.ENCLOSURE, "Enclosure")
    honest = Enclosure.__truediv__
    calls: list[Any] = []

    def counting(self: Any, other: Any) -> Any:
        calls.append((self, other))
        return honest(self, other)

    spec, derivation, cert = s6.medium_table()
    ref = s6.reference(spec, derivation, cert)
    a, b = ref.read(s6.group_flux(0)), ref.read(s6.group_flux(1))
    monkeypatch.setattr(Enclosure, "__truediv__", counting)
    reading = ref.read(s5.c(s5.OBSERVABLE, "Ratio")(s6.group_flux(0), s6.group_flux(1)))
    require(len(calls) == 1, f"the quotient ran {len(calls)} times")
    require(reading == honest(a, b), f"{reading!r} != {honest(a, b)!r}")
    require(_contains(reading, Fraction(30, 7) / Fraction(11, 13)), f"{reading!r} misses the exact ratio")


# ═════════════════════════════════════════════════════════════════════════════
# R6.5 — no derivable bound: refused today (RE-POSED when the user rules)
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("certified", [True, False], ids=["with-a-certificate", "without"])
def test_r6_5_an_observable_the_derivation_cannot_bound_is_refused_today(certified: bool) -> None:
    """RE-POSED WHEN THE USER RULES ON THE UNCERTIFIED READING (#405 P2, the
    open question of 2026-10-03). A family with no derived bound for an
    observable raises ``NotCertified`` (a ``LookupError``) from ``establish``,
    naming the missing certification, and ``read`` lets it through: no float
    without a bound is returned."""
    spec, derivation, cert = s6.medium_table()
    ref = s6.reference(spec, derivation, cert if certified else None)
    error = s6.not_certified()
    require(issubclass(error, LookupError), "NotCertified is not a LookupError")
    with pytest.raises(error, match="certif"):
        ref.read(s5.flux(((0.5, 2.0),)))


# ═════════════════════════════════════════════════════════════════════════════
# R6.6 — the first instance: the exact infinite medium
# ═════════════════════════════════════════════════════════════════════════════

_MIXTURES = (("fuel-2g-upscatter-n2n", sf.fuel), ("fission-only-2g", sf.fission_only), ("fuel-4g", lambda: sf.fuel(4)))


def _medium(mixture: Any) -> Any:
    from orpheus.data.cells import CellCoefficient, Channel
    from orpheus.numerics.question import Eigen
    from orpheus.specification import InfiniteMediumSpecification

    return InfiniteMediumSpecification(0, mixture, Eigen(CellCoefficient.every(Channel.FISSION_EMISSION)))


@pytest.mark.parametrize("make", [m for _, m in _MIXTURES], ids=[n for n, _ in _MIXTURES])
def test_r6_6_the_exact_infinite_medium_reads_the_exact_answer(make: Any) -> None:
    """REFERENCE (exact rationals). For each mixture (2 and 4 groups; anti-#3: a
    1-group k is flux-shape blind), the factory's reference is certified
    ``Valid``; its ``Eigenvalue`` reading contains the exact k_inf of
    ``exact_infinite_medium_of`` (read in the k chart, NEEDS 1) with a bound at
    most one ulp; every group's indicator flux integral contains the exact
    group flux under the module's gauge; and the ratio of the first two groups
    contains the exact ratio, gauge-free; an UNCLAIMED weight (g + 1)/4 reads
    the exact weighted sum (the derivation establishes every flux integral of a
    regionwise-constant weight exactly)."""
    mixture = make()
    exact = s6.exact_medium_of(mixture)
    ref = s6.exact_medium_reference(_medium(mixture))
    require(ref.certificate is not None and s5.state_kind(ref.certificate) == "Valid", f"{ref.certificate!r}")
    k = ref.read(s5.eigenvalue())
    require(_contains(k, exact.k_inf) and k.bound <= math.ulp(k.value), f"k reading {k!r} against {float(exact.k_inf)!r}")
    groups = len(exact.flux)
    for g in range(groups):
        phi = ref.read(s6.group_flux(g, groups=groups))
        require(_contains(phi, exact.flux[g]), f"group {g}: {phi!r} against {float(exact.flux[g])!r}")
    ratio = ref.read(s5.c(s5.OBSERVABLE, "Ratio")(s6.group_flux(0, groups=groups), s6.group_flux(1, groups=groups)))
    require(_contains(ratio, exact.flux[0] / exact.flux[1]), f"ratio {ratio!r}")
    weights = [(g + 1) / 4 for g in range(groups)]  # dyadic, so the exact weighted sum is the reference
    mixed = ref.read(s5.flux((tuple(weights),)))
    require(_contains(mixed, sum((Fraction(w) * f for w, f in zip(weights, exact.flux)), Fraction(0))), f"weighted {mixed!r}")


def test_r6_6_the_factory_refuses_a_question_it_does_not_answer() -> None:
    """The exact family answers the k question on the infinite medium only: a
    geometry specification, and a fixed-source question, are refused."""
    with pytest.raises((TypeError, ValueError)):
        s6.exact_medium_reference(s5.eigen_slab())
    from orpheus.numerics.mesh_free_function import RegionwiseConstant
    from orpheus.numerics.question import FixedSource
    from orpheus.specification import InfiniteMediumSpecification
    import numpy as np

    source = InfiniteMediumSpecification(0, sf.fuel(), FixedSource(RegionwiseConstant(np.ones((1, 2)))))
    with pytest.raises((TypeError, ValueError)):
        s6.exact_medium_reference(source)


# ═════════════════════════════════════════════════════════════════════════════
# R6.7 — the layer
# ═════════════════════════════════════════════════════════════════════════════

_LAYER_SCRIPT = """
import sys
import orpheus.reference.solution as m
print(m.__file__)
print(sorted({k.split('.')[1] for k in sys.modules if k.startswith('orpheus.')}))
"""


def test_r6_7_the_reader_imports_no_derivation() -> None:
    """A cold import of ``orpheus.reference.solution`` loads only
    ``{reference, specification, numerics, data, geometry}``: the derivations
    implement the protocol and import the reference package, never the
    reverse."""
    env = {**os.environ, "PYTHONPATH": str(_ROOT)}
    out = subprocess.run([sys.executable, "-O", "-c", _LAYER_SCRIPT], cwd=_ROOT, env=env, capture_output=True, text=True)
    require(out.returncode == 0, out.stderr)
    file, packages = out.stdout.strip().splitlines()
    require(file.startswith(str(_ROOT)), f"the subprocess imported {file} (L22)")
    loaded = set(ast.literal_eval(packages))
    require("reference" in loaded, f"activation: {loaded}")
    allowed = {"reference", "specification", "numerics", "data", "geometry"}
    require(loaded <= allowed, f"loads {sorted(loaded - allowed)}")
