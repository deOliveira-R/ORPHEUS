r"""A test-local realizer table for the admission mechanism (wave 2, slice S2a, gate 13).

Not a test module (no ``test_`` prefix). The mechanism (``RealizerTable``, the three realizer bases, the lattice) is
exercised on a table the gates build themselves, never on a method's: the real tables are slices S2b, S3 and S4, and
their rows are appended to ``test_admission_matrix.py``.

The default table, ``make_table()``, reproduces the SN columns of section 4.2 as DATA (hand-typed here, never read
from a production table): the specular chain is admitted everywhere except on the hollow inner face, where it is
admitted iff its weight is 1.0; the diffuse chain is refused on a radial outer face (issue 653) and on the inner face
(issue 511); the zero network is refused on the inner face; the source is refused on both radial sites; the deck
realizer admits the mirror and refuses the translation on its topology. Every realizer records its calls in
:class:`Spies`, so a gate can name which unit was consulted with which weight, site and space.

The image type :class:`Img` is the test-local ``OpT``: a sum of named terms with ``+`` (merging a name that
recurs), so that the functor law (the image of a sum is the sum of the images) is assertable with equality.
"""

from __future__ import annotations

import dataclasses
import types
from collections.abc import Callable
from typing import Any

from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.deck import MirrorDeck, TranslationDeck
from tests.gates.transport._boundary_conditions_helpers import bc

#: The opaque "method space" every realize call receives (the table never reads it).
SPACE = object()


def site(kind: str) -> Any:
    """The three geometric sites of section 4.2 as :class:`FaceSite` values."""
    match kind:
        case "CART":
            return bc.FaceSite("xmax", CoordSystem.CARTESIAN, 1, False)
        case "RAD_OUTER":
            return bc.FaceSite("xmax", CoordSystem.CYLINDRICAL, 1, False)
        case "RAD_INNER":
            return bc.FaceSite("xmin", CoordSystem.CYLINDRICAL, 1, True)
    raise ValueError(kind)


SITES = ("CART", "RAD_OUTER", "RAD_INNER")


@dataclasses.dataclass(frozen=True)
class Img:
    """A sum of ``(name, weight)`` terms; ``+`` adds the weights of a name that recurs, in first-appearance order."""

    terms: tuple[tuple[str, float], ...]

    def __add__(self, other: "Img") -> "Img":
        if not isinstance(other, Img):
            return NotImplemented
        acc: dict[str, float] = dict(self.terms)
        for name, weight in other.terms:
            acc[name] = acc.get(name, 0.0) + weight
        return Img(tuple(acc.items()))

    def scaled(self, c: float) -> "Img":
        return Img(tuple((name, c * weight) for name, weight in self.terms))


@dataclasses.dataclass
class Spies:
    admit: list[tuple] = dataclasses.field(default_factory=list)
    realize: list[tuple] = dataclasses.field(default_factory=list)
    sample: list[tuple] = dataclasses.field(default_factory=list)
    deck: list[tuple] = dataclasses.field(default_factory=list)

    def consulted(self, name: str) -> int:
        return sum(1 for entry in self.admit if entry[0] == name)


def make_table(
    *,
    method: str = "mech",
    state_cone: Any = None,
    seal: bool = True,
    with_zero: bool = True,
    with_source: bool = True,
    with_deck: bool = True,
    zero_admit: Callable[[Any], Any] | None = None,
    source_admit: Callable[[Any], Any] | None = None,
    sample_override: Callable[[Any, Any, Any], Any] | None = None,
) -> types.SimpleNamespace:
    """A table (sealed unless ``seal=False``) holding the SN-like test rows; returns the table, the spies and the
    realizer classes (``Spec``, ``Diff``, ``Zero``, ``Src``, ``DeckR``; ``None`` for one left out)."""
    cone = state_cone if state_cone is not None else bc.StateCone.NONNEGATIVE_FLUX
    table = bc.RealizerTable(method, state_cone=cone)
    spies = Spies()
    R, Refused, Axis = bc.Realized, bc.Refused, bc.RefusalAxis
    inner = lambda issue: Refused(Axis.PLUMBING, "the hollow inner face", issue)  # noqa: E731
    outer = lambda: Refused(Axis.PLUMBING, "the radial corner action is not built", 653)  # noqa: E731

    class Spec(bc.ResponseRealizer, table=table, realizes=bc.SpecularChain):
        status = R(True)

        def admit(self, unit, weight, site_, space=None):
            spies.admit.append(("spec", unit, weight, site_, space))
            if site_.site_class is bc.SiteClass.RAD_INNER:
                return R(True) if weight == 1.0 else inner(511)
            return self.status

        def realize(self, unit, weight, site_, space):
            spies.realize.append(("spec", unit, weight, site_, space))
            return Img((("spec", weight),))

    class Diff(bc.ResponseRealizer, table=table, realizes=bc.DiffuseChain):
        status = R(True)

        def admit(self, unit, weight, site_, space=None):
            spies.admit.append(("diff", unit, weight, site_, space))
            match site_.site_class:
                case bc.SiteClass.CART:
                    return self.status
                case bc.SiteClass.RAD_OUTER:
                    return outer()
                case bc.SiteClass.RAD_INNER:
                    return inner(511)

        def realize(self, unit, weight, site_, space):
            spies.realize.append(("diff", unit, weight, site_, space))
            return Img((("diff", weight),))

    Zero = Src = DeckR = None

    if with_zero:
        class Zero(bc.ResponseRealizer, table=table, realizes=bc.ZeroNetwork):  # type: ignore[no-redef]
            status = R(True)

            def admit(self, unit, weight, site_, space=None):
                spies.admit.append(("zero", unit, weight, site_, space))
                if zero_admit is not None:
                    return zero_admit(site_)
                return inner(511) if site_.site_class is bc.SiteClass.RAD_INNER else self.status

            def realize(self, unit, weight, site_, space):
                spies.realize.append(("zero", unit, weight, site_, space))
                return Img((("zero", 1.0),))

    if with_source:
        class Src(bc.SourceRealizer, table=table):  # type: ignore[no-redef]
            status = R(True)

            def admit(self, terms, site_, space=None):
                spies.admit.append(("source", terms, 1.0, site_, space))
                if source_admit is not None:
                    return source_admit(site_)
                match site_.site_class:
                    case bc.SiteClass.CART:
                        return self.status
                    case bc.SiteClass.RAD_OUTER:
                        return outer()
                    case bc.SiteClass.RAD_INNER:
                        return inner(511)

            def sample(self, terms, site_, space):
                spies.sample.append((terms, site_, space))
                if sample_override is not None:
                    return sample_override(terms, site_, space)
                return super().sample(terms, site_, space)

    if with_deck:
        class DeckR(bc.DeckRealizer, table=table):  # type: ignore[no-redef]
            def admit(self, deck, pairing, site_, space=None):
                spies.deck.append(("admit", deck, pairing, site_, space))
                match deck:
                    case MirrorDeck():
                        return R(True)
                    case TranslationDeck():
                        return Refused(Axis.TOPOLOGY, "the operator cannot couple two faces", None)
                    case _:
                        from typing import assert_never

                        assert_never(deck)

            def realize(self, deck, pairing, site_, space):
                spies.deck.append(("realize", deck, pairing, site_, space))
                return Img(((f"deck:{type(deck).__name__}", 1.0),))

    if seal:
        table.seal()
    return types.SimpleNamespace(table=table, spies=spies, Spec=Spec, Diff=Diff, Zero=Zero, Src=Src, DeckR=DeckR)


def values() -> dict[str, Any]:
    """The value classes of section 4.2 as declarations (labels are the matrix rows' names)."""
    W = bc.WallCondition
    return {
        "vacuum": W.vacuum(),
        "specular 1.0": W.specular(1.0),
        "specular 0.5": W.specular(0.5),
        "white 0.5": W.white(0.5),
        "both 0.3 0.5": W.from_amplitudes(0.3, 0.5),
        "inflow 1.0": W.uniform(1.0),
        "specular 0.5 over inflow 1.0": W.specular(0.5) + W.uniform(1.0),
    }
