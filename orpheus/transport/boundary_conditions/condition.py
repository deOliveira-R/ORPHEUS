r"""The exterior conditions of a face (``condition.py`` of the conditions package: the walls and their arithmetic).

**The mathematics.** A face that ends the domain closes the transport equation with the affine trace map
:math:`\gamma_-\psi|_f = R\,\gamma_+\psi|_f + q`: the linear part :math:`R` is a :class:`~.network.ResponseNetwork`
(a point :math:`(s, d)` of the sub-Markov simplex :math:`K = \{s, d \ge 0,\ s + d \le 1\}`, the weights of the specular chain
:math:`P` and the diffuse chain :math:`E\circ Rt`), and :math:`q \ge 0` is an :class:`~.inflow.InflowTerms`.

**The types.**

* :class:`BoundaryCondition`: the abstract base of every exterior condition, the shared family and the method-only
  ones (``ZeroFluxCondition`` is the one that ships; an anisotropic SN kernel would be another). It is a ``ContentIdentity``. Its contract is three accessors, no storage:
  :attr:`~BoundaryCondition.return_fraction` (:math:`\rho = J^-/J^+` of the linear part, the DECLARATION's closed form;
  a zero-flux condition reports ``-1.0``), :attr:`~BoundaryCondition.source` (the affine part :math:`q`, read by
  ``Binding.sources`` for EVERY condition type with no ``isinstance`` on a leaf) and
  :attr:`~BoundaryCondition.face_local` (stated in the face's own frame, so every isometry carries it from face to
  face). :meth:`~BoundaryCondition.invariant_under` is concrete: ``face_local or isometry.is_identity``.
* :class:`WallCondition`: the ONE shared exterior type, the affine pair ``(network, inflow)``. Vacuum, white, specular
  and inflow are thin facade classmethods over it, each returning an EQUAL ``WallCondition`` and never a second type:
  ``WallCondition.specular(0.0)``, ``.white(0.0)`` and ``.vacuum()`` are one value with one digest, because
  :math:`R = 0` is one point of :math:`K` however it is reached.

**Admissibility is a VALUE predicate (decision 50).** ``_cone_refusal(network, inflow)`` is the ONE definition of "inside"
(the network in the simplex :math:`K`, the inflow in the cone :math:`q \ge 0`; the identifier keeps the older name "cone"): every network weight is non-negative (``BoundaryResponseNotPositiveError``), the return fraction
(``ResponseNetwork.return_fraction``, the correctly rounded ``math.fsum``) is at most one (``SubmarkovViolationError``), and every inflow weight is non-negative
(``BoundaryResponseNotPositiveError``). ``__post_init__`` raises what it returns, so a ``WallCondition`` outside them
(and a facade outside it: ``specular(1.2)``, ``white(-0.1)``, a direct two-chain ``(0.6, 0.6)``) does not exist. The
arithmetic (``+``, ``-``, unary ``-``, scalar ``*``, shared with :class:`WallDisplacement` through
``_PairArithmetic``) computes the pair of the result and exits through ONE function, ``_combine``: a ``WallCondition``
iff the pair lies in the cone, else a ``WallDisplacement`` (a sensitivity, never bindable). The sign of a scalar and
the presence of ``-`` decide nothing; what cannot be bound is a value outside :math:`K`, however it was spelled, and
``x + h * delta`` for finite differences re-enters :math:`K` without re-typing. There is no ``__truediv__`` (no
consumer; division would leave the dyadic coordinates), and a method-only condition or a deck defines no operator at
all (``ZeroFluxCondition() + x`` and ``0.5 * MirrorDeck()`` are Python's own ``TypeError``, law L-F).

**What a displacement is.** The cone :math:`K` of walls (``s, d >= 0``, ``s + d <= 1``, every inflow weight ``>= 0``)
is a convex set; a combination of conditions that leaves it (``white(0.3) - specular(0.5)``, ``0.6 * specular +
0.6 * white``, a negative multiple) is a point of the AFFINE HULL of :math:`K` times the source space: a direction in
which a parameter may move, the sensitivity coordinates of a perturbation study, not a wall. :class:`WallDisplacement`
is that value. It is NOT a :class:`BoundaryCondition`, so it is never bindable: ``BoundaryConditions``
refuses it with a ``TypeError`` naming the type, and that refusal is the one narrowing site (no consumer branches on
the sort).

**Why a separate type and not a flag.** The cone is a value predicate (decision 50): the result of arithmetic is a
``WallCondition`` iff its coordinates lie in :math:`K`. Weights of a displacement are signed, which a ``WallCondition``
cannot hold, so the two cannot share a type without a runtime flag that every consumer would have to read.

**Same operators, one exit.** A displacement defines ``+``, ``-``, unary ``-`` and scalar ``*`` through the same
``_PairArithmetic`` as the condition, exiting through ``_combine``: ``WallCondition + WallDisplacement`` is the
condition of the summed pair iff the sum re-enters :math:`K` (a finite-difference step and its inverse returns to a
condition), else a displacement.

**Gotchas of the displacement.**

* :meth:`WallDisplacement.of` ALWAYS builds a displacement, even when the combination lies in :math:`K` (then it is
  still unbindable and unequal to the condition of the same coordinates, since the types differ). It is the
  constructor tests use when the operand signs must not decide the type.
* The network and the inflow of a displacement hold the same drop rule as a condition's (a zero weight is not
  stored, ``-0.0`` is ``0.0``).

**Face-local, no orientation field anywhere.** ``WallCondition`` and the factors hold no ``axis`` and no
``outward_sign``: the realizer reads the face's normal from the face name, so one condition object is correct on every
face of its tag, which a tag-keyed binding requires. ``face_local`` implies invariance under the colouring group only
in ONE dimension: on a Cartesian slab the group acts on faces by the mirror :math:`x \mapsto a + b - x`, which has no
tangential coordinate to reverse; in higher dimension a mirror fixing a face reverses or permutes its tangent frame, so
a ``FaceFrameInflow`` reading a tangential coordinate would not be invariant (not live: a ``Mesh2D`` carries
``UncomputedSymmetries``, where ``bind`` refuses every non-face-local condition on a multi-face tag).

**Gotchas.**

* **The field and the facade have two names.** The dataclass field ``inflow`` is the stored :class:`~.inflow.InflowTerms`
  (decision 47: a field named ``source`` under the base's abstract ``source`` property could not be instantiated), and
  the facade of a uniform source is ``WallCondition.uniform(value)``, named after the coordinate
  ``InflowTerms.uniform``. The specification spelled the facade ``inflow(value)``, a name a class attribute cannot share
  with the field (it was installed after the dataclass was built, and a type checker saw only the field, so every
  consumer was red); the rename was a finding of the S2a review.
* ``__radd__`` and ``__rsub__`` are not defined: the left operand of every pairing the arithmetic admits is itself a
  ``_PairArithmetic`` (a condition or a displacement), whose ``__add__`` / ``__sub__`` answers, so a reflected method
  is unreachable (measured: 0 calls across the three pairings). ``sum([...])`` from the integer ``0`` is therefore
  not supported.
* ``albedo = 0`` is admitted and IS the vacuum (the albedo is a parameter coordinate, and the retired
  ``AlbedoBoundary(0, ...)`` spelling was the same vacuum): the empty network is registered as ``ZeroNetwork`` by every method.
* At the exact boundary of the cone two associations of one sum may differ by a rounding, one a condition and one a
  displacement (the binding door then refuses the latter); nothing in the tree builds a wall at :math:`\rho` within a
  few ulp of one other than by this arithmetic.

Theory: :ref:`bc-architecture-of-record`.

Layer: ``transport`` (L2). Imports ``numerics`` and ``values`` only; ``FaceIsometry`` is imported under
``TYPE_CHECKING``.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Self, assert_never

from orpheus.transport.boundary_conditions.errors import BoundaryError, BoundaryResponseNotPositiveError, SubmarkovViolationError
from orpheus.transport.boundary_conditions.inflow import InflowProfile, InflowTerms
from orpheus.transport.boundary_conditions.network import CHAIN_ORDER, Chain, DiffuseChain, ResponseNetwork, SpecularChain
from orpheus.values.content import ContentIdentity, FrozenMapping
from orpheus.values.scalars import parse_finite_real

if TYPE_CHECKING:
    from orpheus.geometry.colouring import FaceIsometry

__all__ = ["BoundaryCondition", "WallCondition", "WallDisplacement"]


class BoundaryCondition(ContentIdentity, ABC):
    r"""The base of every exterior condition: three accessors, no storage.

    A leaf spells its own storage (``WallCondition`` stores ``inflow`` and exposes it as ``source``).
    """

    @property
    @abstractmethod
    def return_fraction(self) -> float:
        r""":math:`\rho = J^-/J^+` of the linear part: the declaration's closed form (zero flux reports ``-1.0``)."""

    @property
    @abstractmethod
    def source(self) -> InflowTerms:
        r"""The affine part :math:`q`; :meth:`InflowTerms.empty` for a condition with none."""

    @property
    @abstractmethod
    def face_local(self) -> bool:
        """Stated in the face's own frame, so every isometry carries it from face to face."""

    def invariant_under(self, isometry: FaceIsometry) -> bool:
        """Whether the condition is unchanged by ``isometry``: ``face_local or isometry.is_identity``."""
        return self.face_local or isometry.is_identity


def _cone_refusal(network: ResponseNetwork, inflow: InflowTerms) -> BoundaryError | None:
    r"""The ONE definition of "outside the cone :math:`K` x (positive source cone)": the refusal, or None inside."""
    for _, weight in network.weights.items():
        if weight < 0.0:
            return BoundaryResponseNotPositiveError(f"a response weight {weight} is negative: not a physical wall", network)
    if network.return_fraction > 1.0:
        return SubmarkovViolationError(f"the return fraction exceeds 1: not a physical wall (cone point {network.cone_point})", network)
    if inflow.uniform < 0.0 or any(weight < 0.0 for weight, _ in inflow.sampled):
        return BoundaryResponseNotPositiveError("an inflow weight is negative: q is in the positive cone", inflow)
    return None


def _combine(network: ResponseNetwork, inflow: InflowTerms) -> WallCondition | WallDisplacement:
    """The ONE exit of the arithmetic: a condition iff the pair lies in the cone, else a displacement."""
    if _cone_refusal(network, inflow) is None:
        return WallCondition(network, inflow)
    return WallDisplacement(network, inflow)


def _pair_of(operand: object) -> tuple[ResponseNetwork, InflowTerms] | None:
    """``(network, inflow)`` of a condition or a displacement; None for anything else (the operator then answers NotImplemented)."""
    if isinstance(operand, WallCondition | WallDisplacement):
        return operand.network, operand.inflow
    return None


class _PairArithmetic:
    """The affine arithmetic of a ``(network, inflow)`` pair: ``+``, ``-``, unary ``-`` and scalar ``*``, ONE definition.

    Shared by :class:`WallCondition` and :class:`WallDisplacement`, so a chain of arithmetic needs no
    narrowing; every operator exits through :func:`_combine`. ``a - b`` adds ``b`` scaled by ``-1.0``, which negates
    each coordinate exactly, so it equals the coordinate-wise difference.
    """

    network: ResponseNetwork
    inflow: InflowTerms

    def __add__(self, other: object) -> WallCondition | WallDisplacement:
        pair = _pair_of(other)
        if pair is None:
            return NotImplemented
        return _combine(self.network + pair[0], self.inflow + pair[1])

    def __sub__(self, other: object) -> WallCondition | WallDisplacement:
        pair = _pair_of(other)
        if pair is None:
            return NotImplemented
        return _combine(self.network + pair[0] * -1.0, self.inflow + pair[1] * -1.0)

    def __neg__(self) -> WallCondition | WallDisplacement:
        return _combine(self.network * -1.0, self.inflow * -1.0)

    def __mul__(self, scalar: object) -> WallCondition | WallDisplacement:
        if isinstance(scalar, bool) or not isinstance(scalar, int | float):
            return NotImplemented
        return _combine(self.network * scalar, self.inflow * scalar)

    __rmul__ = __mul__


def _amplitude(chain: Chain, specular: float, diffuse: float) -> float:
    """The coordinate of ``chain`` in the pair ``(s, d)``: one exhaustive match, the section of :attr:`ResponseNetwork.cone_point`."""
    match chain:
        case SpecularChain():
            return specular
        case DiffuseChain():
            return diffuse
        case _:
            assert_never(chain)


@dataclass(frozen=True, eq=False)
class WallCondition(_PairArithmetic, BoundaryCondition):
    r"""The ONE shared exterior type: the affine pair :math:`(R, q)` of a wall.

    ``network`` is parsed against the cone (every weight non-negative, ``return_fraction <= 1``) and ``inflow``
    against the positive source cone; a pair outside raises the cone's refusal. The facades build the usual walls.
    """

    network: ResponseNetwork
    inflow: InflowTerms = field(default_factory=InflowTerms.empty)

    def __post_init__(self) -> None:
        refusal = _cone_refusal(self.network, self.inflow)
        if refusal is not None:
            raise refusal

    @property
    def source(self) -> InflowTerms:
        """The affine part :math:`q` (the stored ``inflow``)."""
        return self.inflow

    @property
    def return_fraction(self) -> float:
        r""":math:`\rho` of the network."""
        return self.network.return_fraction

    @property
    def face_local(self) -> bool:
        """The network is face-local by construction; the source is face-local iff no term is in the laboratory frame."""
        return self.inflow.face_local

    @property
    def cone_point(self) -> tuple[float, float]:
        r"""The coordinates :math:`(s, d)` of the network, so a test never rebuilds the chart by ``isinstance``."""
        return self.network.cone_point

    @classmethod
    def vacuum(cls) -> Self:
        """The zero network with no source."""
        return cls(ResponseNetwork(FrozenMapping({})))

    @classmethod
    def specular(cls, albedo: float = 1.0) -> Self:
        r""":math:`\mathrm{albedo}\cdot P`; ``specular(0.0)`` is the vacuum."""
        return cls.from_amplitudes(parse_finite_real(albedo, "albedo"), 0.0)

    @classmethod
    def white(cls, albedo: float = 1.0) -> Self:
        r""":math:`\mathrm{albedo}\cdot(E\circ Rt)`, the Lambertian wall; the first positional argument is the albedo."""
        return cls.from_amplitudes(0.0, parse_finite_real(albedo, "albedo"))

    @classmethod
    def uniform(cls, value: float) -> Self:
        """The zero network over a uniform ``q = value`` (the coordinate ``InflowTerms.uniform``)."""
        return cls(ResponseNetwork(FrozenMapping({})), InflowTerms(parse_finite_real(value, "inflow"), ()))

    @classmethod
    def framed(cls, profile: InflowProfile, weight: float = 1.0) -> Self:
        """The zero network over the sampled source ``weight * profile``."""
        return cls(ResponseNetwork(FrozenMapping({})), InflowTerms(0.0, ((weight, profile),)))

    @classmethod
    def from_amplitudes(cls, specular: float, diffuse: float) -> Self:
        r"""The plain-data edge: the condition at the cone point :math:`(s, d)`; the section of :attr:`cone_point`.

        Takes the two primitives and never a ``WallAmplitudes`` (references are insulated from production).
        """
        s = parse_finite_real(specular, "specular amplitude")
        d = parse_finite_real(diffuse, "diffuse amplitude")
        return cls(ResponseNetwork(FrozenMapping({chain: _amplitude(chain, s, d) for chain in CHAIN_ORDER})))


@dataclass(frozen=True, eq=False)
class WallDisplacement(_PairArithmetic, ContentIdentity):
    """A combination outside the cone: a point of the affine hull of the wall cone times the source space. Never bindable."""

    network: ResponseNetwork
    inflow: InflowTerms

    @classmethod
    def of(cls, *terms: tuple[float, WallCondition]) -> WallDisplacement:
        """The signed combination ``sum_i c_i x_i`` of conditions, ALWAYS a displacement."""
        network = ResponseNetwork(FrozenMapping({}))
        inflow = InflowTerms.empty()
        for coefficient, condition in terms:
            network = network + condition.network * coefficient
            inflow = inflow + condition.inflow * coefficient
        return cls(network, inflow)
