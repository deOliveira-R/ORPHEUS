r"""The response as a tensor network of typed factors, and the closed family of chains (wave 2, slice S2a, gate 9).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 2.4 ("``network.py``", decisions 48 and 49) and law L-M.

A chain is a sequence of FACTORS whose tiers compose from ``"outflow"`` to ``"inflow"`` (the tier marker is
numerics' ``TraceRole``, widened by one member in leaf L4: no second enumeration of the tiers). The family of
chains is a CLOSED sum, ``Chain = SpecularChain | DiffuseChain``: a ``typing.Union`` has no creation hook, so a
check "when the class is created" could cover only the chains the module lists, and an open family would have been
unenforceable (decision 48). Hence the three structural gates below, each with a planted positive control, because a
gate over a list that finds nothing reads exactly like a clean tree (X1).

The claims (``vv-principles``):

* ``test_chains_compose_and_the_family_is_closed`` (law L-M): STRUCTURE. The decidable composition parse
  ``chain_composes`` accepts both shipped chains and refuses a planted wrong-order tuple and a tuple ending on
  ``"current"``; it runs at import of ``network.py`` (a copy whose ``DiffuseChain`` swaps its factors raises
  ``TypeError`` at import, executed in a subprocess so the planted classes never enter this process's registries);
  ``network.py`` defines exactly the chain classes ``get_args(Chain)`` lists; every ``match`` naming a chain class
  ends in ``case _: assert_never(...)``.
* ``test_the_unit_return_of_both_chains``: THEOREM on the declared closed forms (law L-A's rider gates the typed
  image in slice S2b): each chain declares ``unit_return = 1.0``, ``ResponseNetwork.return_fraction`` of ``(s, d)`` is
  ``s + d`` (``math.fsum``), and ``cone_point`` is its inverse.
* ``test_the_response_network_is_a_canonical_form``: the coordinates ARE the canonical form (no leaf collapse, no
  stored partial sum): weights keyed by chain instances, finite reals, ``-0.0`` canonical, a zero weight not stored,
  order not content, ``terms`` in ``CHAIN_ORDER``.

First reds, at the base tip ``4ec95bc6``: ``ModuleNotFoundError: No module named
'orpheus.transport.boundary_conditions'`` on every row.

Mutations (each on a COPY of ``network.py``; the reviewers' battery), named in the law's row:

* swap the factors of ``DiffuseChain`` (the import-time row, ``TypeError``, and ``chain_composes`` on the shipped
  chain);
* add a third chain class carrying ``factors`` without editing the alias (the closure scan names it);
* drop a ``case`` or the ``case _: assert_never(chain)`` arm of ``cone_point`` (the match scan);
* make ``unit_return`` ``0.5`` for ``DiffuseChain``.

Stated limit: ``math.fsum`` and ``+`` agree on TWO terms (a sum of two doubles is already correctly rounded), so with
the shipped two-chain family no mutation between them is observable by any gate; the rows pin the correctly rounded
value, which binds any implementation, and the choice of ``fsum`` matters only for a third chain (the family is closed).
"""

from __future__ import annotations

import math
import subprocess
import sys
import textwrap
from typing import get_args

import pytest

from orpheus.numerics.spaces.angular_trace_space import TraceRole
from orpheus.values.content import FrozenMapping, content_digest
from tests.gates._content_identity_helpers import require
from tests.gates.transport._boundary_conditions_helpers import (
    PACKAGE_DIR,
    ROOT,
    bc,
    chain_classes_in,
    cone_grid,
    matches_naming,
    weights,
)

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/transport/test_response_network.py"
_CHAINS = _HERE + "::test_chains_compose_and_the_family_is_closed"

_CHAIN_NAMES = frozenset({"SpecularChain", "DiffuseChain"})

_PLANTED_THIRD_CHAIN = '''
from typing import ClassVar
class SpecularChain:
    factors: ClassVar = (MirrorPermutation,)
class DiffuseChain:
    factors: ClassVar = (CurrentRetraction, IsotropicSection)
class SquaredChain:
    factors: ClassVar = (MirrorPermutation, MirrorPermutation)
class NotAChain:
    domain: ClassVar = "outflow"
'''
_PLANTED_MATCH_BAD = '''
def f(chain):
    match chain:
        case SpecularChain():
            return 1
        case DiffuseChain():
            return 2
'''
_PLANTED_MATCH_RAISE = '''
def f(chain):
    match chain:
        case SpecularChain():
            return 1
        case _:
            raise ValueError("unregistered")
'''
_PLANTED_MATCH_GOOD = '''
def f(chain):
    match chain:
        case SpecularChain():
            return 1
        case DiffuseChain():
            return 2
        case _:
            assert_never(chain)
'''
_PLANTED_MATCH_OUTSIDE = '''
def f(thing):
    match thing:
        case int():
            return 1
        case _:
            raise TypeError(thing)
'''


def _exec_network(source: str) -> str:
    """Execute a source text of ``network.py`` in a FRESH interpreter (so a planted class never enters this
    process's ``ContentIdentity`` or dataclass registries), as a module of the package; print ``OK`` or the
    exception's class name and message. The interpreter imports the tree under test (``PYTHONPATH``)."""
    script = textwrap.dedent(
        """
        import sys, types
        import orpheus
        from tests.gates.transport._boundary_conditions_helpers import load
        load()  # the real package, and the parent plumbing the shim owns (planted classes then share this process only)
        src = sys.stdin.read()
        name = "orpheus.transport.boundary_conditions._planted_network"
        module = types.ModuleType(name)
        module.__package__ = "orpheus.transport.boundary_conditions"
        sys.modules[name] = module  # dataclasses reads the defining module for string annotations
        try:
            exec(compile(src, "planted_network.py", "exec"), module.__dict__)
        except BaseException as err:
            print("RAISED", type(err).__name__, err)
        else:
            print("OK")
        print("FILE", orpheus.__file__)
        """
    )
    import os
    import tempfile

    with tempfile.TemporaryDirectory() as cwd:
        out = subprocess.run(
            [sys.executable, "-O", "-c", script], input=source, cwd=cwd, capture_output=True, text=True, timeout=300,
            env={**os.environ, "PYTHONPATH": str(ROOT)},
        )
    require(out.returncode == 0, f"the planted-module subprocess failed:\n{out.stderr[-2000:]}")
    lines = out.stdout.splitlines()
    files = [line.split(" ", 1)[1] for line in lines if line.startswith("FILE ")]
    require(files and files[0].startswith(str(ROOT)), f"the subprocess imported {files}, not the tree under test")
    return lines[0]


# ── gate 9: chains compose and the family is closed (law L-M) ────────────────


def test_chains_compose_and_the_family_is_closed() -> None:
    """Law L-M, three structural facts and one behaviour.

    * ``chain_composes`` accepts ``(MirrorPermutation,)`` and ``(CurrentRetraction, IsotropicSection)`` and refuses
      ``(IsotropicSection, CurrentRetraction)``, ``(CurrentRetraction,)`` and ``(IsotropicSection,)`` (the planted
      wrong-order tuple and the tuples ending on ``"current"`` or starting on it), a doubled mirror and a mirror
      followed by a section (tiers do not meet).
    * ``get_args(TraceRole)`` is the four tiers of one face's ladder, and every tier a factor names is among them.
    * ``network.py`` DEFINES exactly the chain classes ``get_args(Chain)`` lists (classes with a ``factors`` class
      attribute); the scanner's controls: a planted source with a third chain class reports it, one with the
      exact set reports nothing.
    * every ``match`` naming a chain class ends in a guard-free ``case _: assert_never(...)``; population at least
      one (``cone_point``'s match); controls: a match without the arm and a match whose catch-all raises are
      violations, the correct match is clean, a match naming no chain is outside the population.
    * at import, a copy of ``network.py`` whose ``DiffuseChain`` swaps its factors raises ``TypeError`` (an explicit
      raise: the suite runs under ``python -O``); the unmodified source executes without error (the control).

    Reddens on: the five mutations of the module docstring. At the base tip the first read of ``bc`` is the red.
    """
    spec, diff = bc.SpecularChain, bc.DiffuseChain
    require(set(get_args(bc.Chain)) == {spec, diff} and len(get_args(bc.Chain)) == 2, f"Chain is {get_args(bc.Chain)}")
    require(tuple(type(c) for c in bc.CHAIN_ORDER) == (spec, diff), f"CHAIN_ORDER is {bc.CHAIN_ORDER}")
    mirror, retraction, section = bc.MirrorPermutation, bc.CurrentRetraction, bc.IsotropicSection
    require(set(get_args(bc.Factor)) == {mirror, retraction, section}, f"Factor is {get_args(bc.Factor)}")

    for chain in get_args(bc.Chain):
        require(bc.chain_composes(chain.factors), f"{chain.__name__}.factors do not compose")
    require(spec.factors == (mirror,) and diff.factors == (retraction, section), "the shipped chains are s P and d (E o Rt)")
    require(bc.chain_composes((mirror,)) and bc.chain_composes((retraction, section)), "a shipped chain is refused")
    for planted in ((section, retraction), (retraction,), (section,), (mirror, mirror), (mirror, section), (retraction, retraction)):
        require(not bc.chain_composes(planted), f"chain_composes accepted the planted tuple {[f.__name__ for f in planted]}")

    tiers = {"full", "outflow", "current", "inflow"}
    require(set(get_args(TraceRole)) == tiers, f"TraceRole has {get_args(TraceRole)}: the four tiers of one face's ladder")
    for factor in get_args(bc.Factor):
        require({factor.domain, factor.codomain} <= set(get_args(TraceRole)), f"{factor.__name__} names a tier TraceRole lacks")
    require((mirror.domain, mirror.codomain) == ("outflow", "inflow"), "P: outflow -> inflow")
    require((retraction.domain, retraction.codomain) == ("outflow", "current"), "Rt: outflow -> current")
    require((section.domain, section.codomain) == ("current", "inflow"), "E: current -> inflow")

    # the closure scan, with its controls (X1)
    require(chain_classes_in(_PLANTED_THIRD_CHAIN) == {"SpecularChain", "DiffuseChain", "SquaredChain"},
            "control: the closure scanner does not see a planted third chain")
    require(chain_classes_in(_PLANTED_MATCH_GOOD) == set(), "control: the closure scanner reports a source with no chain class")
    source = (PACKAGE_DIR / "network.py").read_text(encoding="utf-8")
    listed = {c.__name__ for c in get_args(bc.Chain)}
    require(chain_classes_in(source) == listed, f"network.py defines the chain classes {sorted(chain_classes_in(source))}, Chain lists {sorted(listed)}")

    # the exhaustive-match scan, with its controls
    require(matches_naming(_PLANTED_MATCH_BAD, _CHAIN_NAMES) == (1, ["f:3"]), "control: a match without the arm is not a violation")
    require(matches_naming(_PLANTED_MATCH_RAISE, _CHAIN_NAMES) == (1, ["f:3"]), "control: a catch-all that raises is not a violation")
    require(matches_naming(_PLANTED_MATCH_GOOD, _CHAIN_NAMES) == (1, []), "control: the correct match is flagged")
    require(matches_naming(_PLANTED_MATCH_OUTSIDE, _CHAIN_NAMES) == (0, []), "control: a match naming no chain is in the population")
    population, violations, sites = 0, [], []
    files = sorted(PACKAGE_DIR.glob("*.py"))
    require(len(files) == 10, f"the scan read {len(files)} package files, not the 10 of the table")
    for path in files:
        n, bad = matches_naming(path.read_text(encoding="utf-8"), _CHAIN_NAMES)
        population += n
        violations += [f"{path.name}:{b}" for b in bad]
        if n:
            sites.append(path.name)
    require(population >= 1 and "network.py" in sites, f"cone_point's match over the chain was not found; sites: {sites}")
    require(not violations, f"a match over a chain lacks `case _: assert_never(...)`: {violations}")

    # the composition check runs at import: a planted swap raises TypeError; the unmodified source does not
    swapped = source.replace("(CurrentRetraction, IsotropicSection)", "(IsotropicSection, CurrentRetraction)")
    require(swapped != source, "the mutant could not be built: network.py does not spell DiffuseChain's factors as (CurrentRetraction, IsotropicSection)")
    require(_exec_network(source) == "OK", "control: the unmodified network.py does not execute cleanly as a module")
    require(_exec_network(swapped).startswith("RAISED TypeError"), "a DiffuseChain with swapped factors imports without a TypeError")


# ── the unit return, the chart and its inverse ───────────────────────────────


@pytest.mark.rests_on(_CHAINS)
def test_the_unit_return_of_both_chains() -> None:
    """Each chain declares ``unit_return == 1.0``; ``ResponseNetwork.return_fraction`` of the point ``(s, d)`` is
    ``math.fsum([s, d])`` (correctly rounded) over the whole cone grid; ``cone_point`` is its inverse on the grid
    (one ``match``, total because the closed sum IS the pair); the empty network returns 0.0 and sits at the origin;
    ``0.9, 0.1`` sums to exactly 1.0 (the float rule of the cone sum).

    Reddens on: ``DiffuseChain.unit_return = 0.5``; a ``cone_point`` that
    returns the weights in insertion order.
    """
    require(bc.SpecularChain.unit_return == 1.0 and bc.DiffuseChain.unit_return == 1.0, "both shipped chains declare unit_return 1.0")
    grid = cone_grid()
    require(len(grid) == 48, f"the grid has {len(grid)} points, not the 48 (45 lattice + 3 non-dyadic) it was written with")
    for s, d in grid:
        net = weights(specular=s, diffuse=d)
        require(net.return_fraction == math.fsum([s, d]), f"return_fraction at {(s, d)} is {net.return_fraction!r}, fsum {math.fsum([s, d])!r}")
        require(net.cone_point == (s, d), f"cone_point at {(s, d)} reads {net.cone_point}")
    require(weights().return_fraction == 0.0 and weights().cone_point == (0.0, 0.0), "the empty network returns nothing")
    require(weights(specular=0.9, diffuse=0.1).return_fraction == 1.0, "0.9 + 0.1 is correctly rounded to 1.0")
    require(weights(specular=0.5).cone_point == (0.5, 0.0) and weights(diffuse=0.5).cone_point == (0.0, 0.5), "a single chain reads on its own axis")


# ── the canonical form ───────────────────────────────────────────────────────


@pytest.mark.rests_on(_CHAINS)
def test_the_response_network_is_a_canonical_form() -> None:
    """The weights are parsed (chain keys, finite reals, ``-0.0`` canonical, a zero weight not stored), their
    order is not content, ``terms`` follow ``CHAIN_ORDER``, and ``+`` and scalar ``*`` act chain by chain with the
    zero-weight drop reaching the result. Signed weights are admitted here (the cone is the condition's, so a
    displacement can hold them). The marker ``ZeroNetwork`` is a value with identity by type.

    Reddens on: a network that stores ``0.0`` or ``-0.0`` (a zero weight appears in ``terms`` and the digest
    differs from the empty network's); one that keeps the caller's order (``terms`` or the digest reorder); one that
    accepts a non-chain key, ``NaN``, an infinity or a ``bool`` weight.
    """
    forward = weights(specular=0.5, diffuse=0.25)
    backward = bc.ResponseNetwork(FrozenMapping({bc.DiffuseChain(): 0.25, bc.SpecularChain(): 0.5}))
    require(forward == backward and content_digest(forward) == content_digest(backward), "the order of the weights is content")
    for net in (forward, backward):
        require(net.terms == ((0.5, bc.SpecularChain()), (0.25, bc.DiffuseChain())), f"terms are {net.terms}, not in CHAIN_ORDER")
    require(isinstance(forward.weights, FrozenMapping), "the weights are a FrozenMapping")
    empty = bc.ResponseNetwork(FrozenMapping({}))
    for zero in (weights(specular=0.0), weights(specular=-0.0, diffuse=0.0), weights(specular=0.0, diffuse=0.5) + weights(diffuse=0.5) * -1.0):
        require(zero == empty and content_digest(zero) == content_digest(empty), f"a zero weight is stored: {zero!r}")
        require(zero.terms == () and len(zero.weights) == 0, "the empty network has terms")

    for bad in (FrozenMapping({bc.MirrorPermutation(): 1.0}), FrozenMapping({"specular": 1.0}), FrozenMapping({bc.SpecularChain: 1.0})):
        with pytest.raises(TypeError):
            bc.ResponseNetwork(bad)
    for bad_weight in (float("nan"), float("inf"), True, "0.5"):
        with pytest.raises((TypeError, ValueError)):
            bc.ResponseNetwork(FrozenMapping({bc.SpecularChain(): bad_weight}))
    signed = weights(specular=-0.5, diffuse=2.0)
    require(signed.cone_point == (-0.5, 2.0), "a signed network is admitted")

    require((forward + forward).cone_point == (1.0, 0.5), "+ adds chain by chain")
    require((forward + weights(diffuse=0.5)).cone_point == (0.5, 0.75), "+ merges the chains the operands share")
    require(forward + forward * -1.0 == empty, "x + (-1) x is the empty network")
    require((forward * 2.0).cone_point == (1.0, 0.5), "scalar * scales every weight")
    require(forward * 0.0 == empty and forward * -0.0 == empty, "a zero scalar drops every weight")
    require(bc.ZeroNetwork() == bc.ZeroNetwork() and hash(bc.ZeroNetwork()) == hash(bc.ZeroNetwork()), "ZeroNetwork is a value by type")
    with pytest.raises(TypeError):
        forward + 1.0  # type: ignore[operator]  # a number is not a network


def test_the_canonical_order_is_the_alias() -> None:
    """``CHAIN_ORDER`` is DERIVED from ``get_args(Chain)``: it holds exactly the chain types of the closed sum, once
    each, in the alias's order, so a chain added to the alias is ordered (and its weight kept) with no second edit.

    Reddens on: a hand-written ``CHAIN_ORDER`` (the mutation re-spells it; the next chain added to the alias would be
    dropped by ``ResponseNetwork`` silently).
    """
    members = get_args(bc.Chain)
    require(len(members) >= 2, f"control: the closed sum has {len(members)} members")
    require(tuple(type(chain) for chain in bc.CHAIN_ORDER) == members, f"CHAIN_ORDER is {[type(c).__name__ for c in bc.CHAIN_ORDER]}, the alias {[m.__name__ for m in members]}")
    source = (PACKAGE_DIR / "network.py").read_text()
    require("CHAIN_ORDER: tuple[Chain, ...] = tuple(" in source and "get_args(Chain)" in source.split("CHAIN_ORDER: tuple[Chain, ...] =")[1].splitlines()[0], "CHAIN_ORDER is not read from the alias at its definition")
