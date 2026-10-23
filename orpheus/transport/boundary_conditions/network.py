r"""The response of a wall as a tensor network of typed factors (``network.py`` of the conditions package).

**What this module is.** A wall that ends the domain closes the transport equation with the affine trace map

.. math:: \gamma_-\psi|_f = R\,\gamma_+\psi|_f + q ,

and this module holds the LINEAR part :math:`R`, never :math:`q` (that is :mod:`.inflow`). :math:`R` is a network of
typed factor tensors on the rungs of one face's trace ladder, the ladder of
:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`: the outflow half-trace :math:`\Gamma_+`
(tier ``"outflow"``), the partial current :math:`S` (``"current"``) and the inflow half-trace :math:`\Gamma_-`
(``"inflow"``). The tier marker is numerics' own :data:`~orpheus.numerics.spaces.angular_trace_space.TraceRole`: no
second enumeration of the tiers exists here (decision 49), and "embedding" is not an operator name in this project.

**The factors** (each a field-less frozen dataclass: a factor is a TYPE, its realization is the method's):

* :class:`MirrorPermutation` :math:`P`: :math:`\Gamma_+ \to \Gamma_-`, :math:`(P\phi)(x,\Omega) =
  \phi(x, \Omega - 2(\Omega\cdot n_f)n_f)`, the ordinate permutation of the face's mirror plane;
* :class:`CurrentRetraction` :math:`Rt`: :math:`\Gamma_+ \to S`, integration against the trace measure
  :math:`|\Omega\cdot n_f|\,w`;
* :class:`IsotropicSection` :math:`E`: :math:`S \to \Gamma_-`, :math:`EJ = J/M_-` with :math:`M_-` the INFLOW mass
  :math:`\sum_{\Omega\cdot n_f<0}|\Omega\cdot n_f|\,w`. :math:`Rt_{in}\circ E = \mathrm{id}`: the split pair of
  :mod:`orpheus.numerics.axis_marginal` (``AxisRetractionOperator``, ``AxisSectionOperator``).

**The chains.** Two chains compose from :math:`\Gamma_+` to :math:`\Gamma_-`: the specular chain :math:`P` and the
diffuse chain :math:`E\circ Rt` (listed in the order the factors are applied). They form the CLOSED sum
``Chain = SpecularChain | DiffuseChain`` (decision 48): a ``typing.Union`` has no creation hook, so an open family would
be unenforceable, and the family is instead closed the way the decks are: every place that reads a chain does an
exhaustive ``match`` ending in ``assert_never`` (here, :attr:`ResponseNetwork.cone_point`), and a third chain is an
edit of the alias, after which every such match reddens until it handles the new member. A chain declares its
``unit_return`` (the return fraction of the chain at weight one) and the left-eigenvector rider of the theory page
(``R^T (w|mu|)_- = rho (w|mu|)_+``) witnesses it on the realized image; both shipped chains declare ``1.0`` (the split
pair, and ``P`` preserving the partial current). A method-only kernel (the neutron-guide wall, an anisotropic SN
kernel, zero flux) is another :class:`~.condition.BoundaryCondition` type and never a chain.

**Decidable composition (law L-M).** :func:`chain_composes` is the predicate "the factors' tiers compose from
``"outflow"`` to ``"inflow"``". It is applied to every member of the closed sum when this module is imported, so a
chain such as ``(IsotropicSection, CurrentRetraction)`` or one that ends on ``"current"`` is unspellable in the shipped
module; the check raises ``TypeError`` explicitly (not ``assert``: the suite runs under ``python -O``). The admission
test for a NEW chain is the two-test criterion of ``orpheus/geometry/boundary/_factors.py`` (multiplicativity, then
"is it a quotient?"): a mirror-like chain is constitutive (it belongs to :math:`R`), an averaging chain is a kernel,
neither is a deck.

**The network.** :class:`ResponseNetwork` is :math:`R = \sum_c w_c\,c` as a canonical weight map. The canonical form
IS the coordinates: weights are finite reals (``-0.0`` is ``0.0``), a zero weight is not stored, and the entry order is
not content (a :class:`~orpheus.values.content.FrozenMapping` keyed by the field-less chains is order-free, and the
mapping is stored in :data:`CHAIN_ORDER` so every reader sees one order). There is no collapse of a one-term sum to a
leaf, so no stored value holds an intermediate rounding that a second addition rounds again. The weights ARE the point
:math:`(s, d)` of the convex cone :math:`K = \{s, d \ge 0,\ s + d \le 1\}` that :attr:`ResponseNetwork.cone_point`
reads, by one exhaustive match: the closed sum is the pair, so the chart is total and not lossy.

**Why a marker for the empty network.** Realization is a functor on the network and the image of the empty sum is the
method's zero map, which is method-specific (SN's ``ZeroMap``, diffusion's zero scalar). :class:`ZeroNetwork` is the key
under which a method registers it. ``WallCondition.specular(0.0)``, ``.white(0.0)`` and ``.vacuum()`` have the same
empty network, hence one value and one digest.

**Gotchas.**

* ``weights`` is keyed by chain INSTANCES (field-less singletons) while a method's realizer table is keyed by chain
  TYPES: the functor looks up ``type(chain)`` and carries the weight beside it (type-keyed dispatch, instance-carried
  data), which stays correct the day a chain gains a parameter.
* Summing the images of several chains is done in :data:`CHAIN_ORDER` (the canonical order), so a method's
  realization of a two-chain wall does not depend on how the wall was built.
* :attr:`ResponseNetwork.return_fraction` is ``math.fsum``, correctly rounded: ``fsum([0.9, 0.1]) <= 1.0`` is True,
  so an intended unit wall written as ``0.9`` plus ``0.1`` is inside the cone, and ``fsum([1.0, 2**-52]) <= 1.0`` is
  False, so a true excess of one double step is outside it.
* ``ResponseNetwork`` parses only its own keys and reals. The CONE (non-negative weights, return fraction at most one)
  is the condition's: a network with negative weights is the linear part of a displacement.

Layer: ``transport`` (L2). Imports ``numerics`` (the tier marker) and ``values`` only.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import ClassVar, TypeAlias, assert_never, get_args

from orpheus.numerics.spaces.angular_trace_space import TraceRole
from orpheus.values.content import ContentIdentity, FrozenMapping
from orpheus.values.scalars import parse_finite_real

__all__ = [
    "CHAIN_ORDER",
    "Chain",
    "CurrentRetraction",
    "DiffuseChain",
    "Factor",
    "IsotropicSection",
    "MirrorPermutation",
    "ResponseNetwork",
    "SpecularChain",
    "ZeroNetwork",
    "chain_composes",
]


@dataclass(frozen=True)
class MirrorPermutation:
    r""":math:`P`, the ordinate permutation of the face's mirror plane, :math:`\Gamma_+ \to \Gamma_-`.

    Face-local: it holds no axis and no outward sign, because the realizer reads the face's normal from the face name.
    """

    domain: ClassVar[TraceRole] = "outflow"
    codomain: ClassVar[TraceRole] = "inflow"


@dataclass(frozen=True)
class CurrentRetraction:
    r""":math:`Rt`, the retraction of the split pair: :math:`\Gamma_+ \to S`, integration against :math:`|\Omega\cdot n_f|\,w`."""

    domain: ClassVar[TraceRole] = "outflow"
    codomain: ClassVar[TraceRole] = "current"


@dataclass(frozen=True)
class IsotropicSection:
    r""":math:`E`, the section of the split pair: :math:`S \to \Gamma_-`, :math:`EJ = J/M_-` (divides by the INFLOW mass)."""

    domain: ClassVar[TraceRole] = "current"
    codomain: ClassVar[TraceRole] = "inflow"


Factor: TypeAlias = MirrorPermutation | CurrentRetraction | IsotropicSection


@dataclass(frozen=True)
class SpecularChain:
    r""":math:`s\,P`, the specular chain; ``factors`` are listed in the order they are applied."""

    factors: ClassVar[tuple[type[Factor], ...]] = (MirrorPermutation,)
    unit_return: ClassVar[float] = 1.0


@dataclass(frozen=True)
class DiffuseChain:
    r""":math:`d\,(E\circ Rt)`, the diffuse (Lambertian) chain; ``factors`` are listed in the order they are applied."""

    factors: ClassVar[tuple[type[Factor], ...]] = (CurrentRetraction, IsotropicSection)
    unit_return: ClassVar[float] = 1.0


#: The chains of the shared family, a CLOSED sum: editing this alias is the only way to add a chain.
Chain: TypeAlias = SpecularChain | DiffuseChain

#: The canonical order of the chains (the order of the alias): the order of :attr:`ResponseNetwork.terms` and of the sum of
#: realized images. DERIVED from the alias, never a second enumeration of it (a chain added to the alias is ordered with it).
CHAIN_ORDER: tuple[Chain, ...] = tuple(chain() for chain in get_args(Chain))


def chain_composes(factors: tuple[type[Factor], ...]) -> bool:
    r"""Whether ``factors``, applied in order, compose from the ``"outflow"`` rung to the ``"inflow"`` rung.

    ``factors[0].domain == "outflow"``, ``factors[i].codomain == factors[i+1].domain`` and
    ``factors[-1].codomain == "inflow"``. The empty tuple does not compose (the empty network is :class:`ZeroNetwork`,
    not an empty chain).
    """
    if not factors or factors[0].domain != "outflow" or factors[-1].codomain != "inflow":
        return False
    return all(a.codomain == b.domain for a, b in zip(factors, factors[1:], strict=False))


def _refuse_a_chain_that_does_not_compose() -> None:
    """Run once at import: every member of the closed sum ``Chain`` composes from outflow to inflow."""
    for chain in get_args(Chain):
        if not chain_composes(chain.factors):
            raise TypeError(f"{chain.__name__}: its factors {[f.__name__ for f in chain.factors]} do not compose from outflow to inflow")


_refuse_a_chain_that_does_not_compose()


@dataclass(frozen=True)
class ZeroNetwork:
    r"""The marker of the empty network (:math:`R = 0`): the key under which a method registers its zero image."""


@dataclass(frozen=True, eq=False)
class ResponseNetwork(ContentIdentity):
    r"""The linear part :math:`R = \sum_c w_c\,c` of a wall, a weighted sum of chains in canonical form.

    ``weights`` is parsed on construction: every key is a member of the closed sum :data:`Chain`, every weight a finite
    real (``-0.0`` becomes ``0.0``), a zero weight is dropped, and the mapping is stored in :data:`CHAIN_ORDER`. The
    cone is NOT checked here (see the module docstring).
    """

    weights: FrozenMapping[Chain, float]

    def __post_init__(self) -> None:
        parsed: dict[Chain, float] = {}
        for chain, weight in dict(self.weights).items():
            if not isinstance(chain, get_args(Chain)):
                raise TypeError(f"ResponseNetwork: a key is a chain (one of {[c.__name__ for c in get_args(Chain)]}), got {chain!r}")
            weight = parse_finite_real(weight, f"ResponseNetwork weight of {type(chain).__name__}")
            if weight != 0.0:
                parsed[chain] = weight
        object.__setattr__(self, "weights", FrozenMapping({c: parsed[c] for c in CHAIN_ORDER if c in parsed}))

    @property
    def cone_point(self) -> tuple[float, float]:
        r"""The coordinates :math:`(s, d)`; ``0.0`` for a chain that is absent. One exhaustive match over the closed sum."""
        specular = diffuse = 0.0
        for chain, weight in self.weights.items():
            match chain:
                case SpecularChain():
                    specular = weight
                case DiffuseChain():
                    diffuse = weight
                case _:
                    assert_never(chain)
        return (specular, diffuse)

    @property
    def return_fraction(self) -> float:
        r""":math:`\rho = \sum_c w_c\,\mathrm{unit\_return}(c)`, correctly rounded (``math.fsum``)."""
        return math.fsum(weight * chain.unit_return for chain, weight in self.weights.items())

    @property
    def terms(self) -> tuple[tuple[float, Chain], ...]:
        """The ``(weight, chain)`` pairs in :data:`CHAIN_ORDER`; the realizers' functor reads this."""
        return tuple((self.weights[chain], chain) for chain in CHAIN_ORDER if chain in self.weights)

    def __add__(self, other: object) -> ResponseNetwork:
        """The network whose weight on each chain is the sum of the two (the sum of the maps)."""
        if not isinstance(other, ResponseNetwork):
            return NotImplemented
        summed = dict(self.weights)
        for chain, weight in other.weights.items():
            summed[chain] = summed.get(chain, 0.0) + weight
        return ResponseNetwork(FrozenMapping(summed))

    def __mul__(self, scalar: object) -> ResponseNetwork:
        """Every weight multiplied by a real scalar (a ``bool`` is not one)."""
        if isinstance(scalar, bool) or not isinstance(scalar, int | float):
            return NotImplemented
        return ResponseNetwork(FrozenMapping({chain: scalar * weight for chain, weight in self.weights.items()}))
