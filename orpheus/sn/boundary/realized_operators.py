r"""The typed realized boundary operators of SN: the images of the response network's chains (``realized_operators.py``).

**The mathematics.** A face of an SN problem is closed by the affine trace map :math:`\gamma_-\psi|_f = R\,\gamma_+\psi|_{f'} + q`
(:mod:`orpheus.transport.boundary_conditions`), whose linear part :math:`R` is a weighted sum of factor chains on the rungs of
the face's trace ladder, outflow half-trace :math:`\Gamma_+`, partial current :math:`S`, inflow half-trace :math:`\Gamma_-`
(:mod:`orpheus.numerics.spaces.angular_trace_space`). Realization is a FUNCTOR on that network: the SN table
(:data:`orpheus.sn.boundary.realizer.SN_REALIZERS`) sends each chain to an operator :math:`\Gamma_+ \to \Gamma_-`, and the image
of a sum is the sum of the images.

* the empty network (:math:`0`) goes to the zero map :math:`\Gamma_+\to\Gamma_-`, the type :class:`ZeroMap`;
* :math:`s\,P` (``SpecularChain``) goes to :math:`s\,P_\pi` with :math:`\pi` certified, the type :class:`ScaledOrdinatePermutation`;
* :math:`d\,(E\circ R_t)` (``DiffuseChain``) goes to :math:`d\,(\text{emission}\circ\text{contraction})`, the type :class:`RankOneReemission`;
* a sum of chains goes to the sum of the images, the type :class:`RealizedSum`.

* :math:`P_\pi` is the ordinate permutation of a rigid motion: for the mirror of the face's own plane it sends the ordinate
  :math:`\Omega` to its reflection :math:`\Omega - 2(\Omega\cdot n_f)\,n_f` (a wall's constitutive specular return, and the
  mirror DECK, which is the same matrix at unit weight); for the translation wrap of a periodic pair it is the identity
  relabelling between two different index sets, :math:`\Gamma_+(f') \to \Gamma_-(f)`.
* The Lambertian image is the factored rank-one kernel. :math:`R_t\phi = \sum_{\Omega\cdot n_f>0} w\,|\Omega\cdot n_f|\,\phi`
  contracts :math:`\Gamma_+` to the partial current, and :math:`E J = J/M_-` re-emits it uniformly on :math:`\Gamma_-`, with
  :math:`M_- = \sum_{\Omega\cdot n_f<0} w\,|\Omega\cdot n_f|` the INFLOW mass (the divisor rule, below).

**The facts are properties of the type.** The five questions the SN consumers ask of a realized face (does it return
everything, :attr:`~RealizedBoundary.return_fraction`; does it relabel ordinates, :attr:`~RealizedBoundary.permutes_ordinates`;
whose outflow does it consume, :attr:`~RealizedBoundary.domain_face`; is it the zero map, :attr:`~RealizedBoundary.is_vacuum`;
and what operator is it, :attr:`~RealizedBoundary.inner`) used to be answered by asking the LAW the operator came from
(``law_permutes_ordinates``, ``response_kernel.is_zero``, ``.amplitude``, ``.geometry_map.domain_face``, a ``kind`` string). The
realizer that builds the operator already knows each fact, so each is DERIVED from the type's parameters, never stored beside
them. A type holds ONLY its parameters ``(alpha, ..., face, trace)``; the operator is built from spaces READ OFF the trace
ladder (``trace.outflow_space(f)``, ``trace.current_space(f)``, ``trace.inflow_space(f)``), so a lying ``alpha`` is unspellable
and ``dataclasses.replace(op, alpha=...)`` REBUILDS the operator rather than editing a fact. The types are frozen and compare
and hash by their parameters (the type, the face, the partner, ``alpha`` and the permutation's index bytes), not by a stored
operator: a ``LinearOperator`` has no content identity, and a trace space compares by ``(name, shape)`` only, so two traces of
different quadratures with the same total size would otherwise be indistinguishable. The old shim's ``str`` equality arm (the
retired ``op == "reflective"``, #553) is not carried.

**One body per builder (decision 45).** This module OWNS the one body of each operator builder. The realizer's four builders
(``_narrowed_zero_operator``, ``_deck_kernel``, ``_checked_angular_average`` and ``_specular_crossed_diagnosis``) retired in the
commit that added this file, and the only helper shared by the constructors is :func:`_attenuated_kernel_operator` (scale and
role stamp). The import edge is ``realizer -> realized_operators`` only; this module imports nothing from ``realizer.py``.

**The divisor rule.** ``E`` divides by the INFLOW mass (a plain ``numpy`` sum of the inflow half-space's metric), which is what
makes ``Rt_in o E = id`` and the return fraction exactly ``alpha``. The shipped operator before this carve divided by the
OUTFLOW sum, so a Lambertian on a quadrature with :math:`M_+ \neq M_-` had return fraction :math:`d\,M_-/M_+`, silently.
``[M]`` ``tests/gates/sn/architecture/test_realized_boundary.py::test_the_mass_census_inflow_equals_outflow_on_every_registry_rule``:
:math:`M_+ = M_-` bitwise on 70 of 70 (registry rule, face) pairs (35 rules, 2 faces), so the carry is bit-identical wherever the
two masses are bit-equal, which is every shipped rule (plan section 8, row R-11).

**Rider (b) of the WF2.2 record, decided (S2b): no ``OrdinatePairing`` carrier.** The plan asked whether
``ScaledOrdinatePermutation`` should take a frozen ``OrdinatePairing(permutation, motion, axis)`` whose only constructor is
``ordinate_pairing``, so that an uncertified ``quadrature.ordinate_permutation(motion)`` could not be realized. Decision: not
adopted in this slice, for three structural reasons. (1) ``ordinate_pairing`` certifies an AXIS REFLECTION and refuses every other
motion, so the translation wrap's identity permutation (the other consumer of this type) has no constructor under the proposed
rule, and the carrier would need a second constructor, which is no longer "only". (2) Python cannot make a frozen dataclass's
constructor private, so "only constructor" would be a convention and not a structure, which is exactly the prose that X3 of the
instrument doctrine refuses. (3) The production paths that build a mirror image call ``ordinate_pairing`` and nothing else
(:class:`~orpheus.sn.boundary.realizer.SNDecks`, ``SNSpecularChain``), and the translation wrap is certified at construction by
the crossing check of :func:`_crossing` (the quotient identification :math:`\Gamma_+(f') \equiv \Gamma_-(f)`, executed rather
than asserted). The strict ``xfail`` of ``test_an_uncertified_permutation_cannot_be_realized`` stays, with its RECORD row, as the
measured statement of the gap. Trigger to revisit: a second realizer (CP, MoC or MC mirror) consuming a permutation, or #657.

**Rider (c): the realization functor checks the tiers.** Each leaf type names the chain of the closed sum it realizes
(:attr:`ScaledOrdinatePermutation.chain`, :attr:`RankOneReemission.chain`) and refuses, at construction, a trace whose
``outflow`` and ``inflow`` half-spaces do not carry the tiers the chain's first and last factors declare
(``chain.factors[0].domain == outflow.role``, ``chain.factors[-1].codomain == inflow.role``). The check reads
``AngularFaceTraceSpace.role`` and never dispatches on it. The declared composition of
:mod:`orpheus.transport.boundary_conditions.network` and the realized composition of this module therefore meet at one check
instead of being two independent spellings of "the rungs compose".

Gotchas
=======

* ``apply`` takes the face's outflow half-trace of the DOMAIN face (:attr:`~RealizedBoundary.domain_face`), which is the
  partner for a translation wrap and the face itself otherwise; ``block_role`` is :attr:`~orpheus.numerics.operator.BlockRole.BOUNDARY`.
* A zero scale is the :class:`ZeroMap`, refused as an ``alpha`` here (``alpha`` lies in ``(0, 1]``); a super-Markov or non-finite
  scale is refused at construction.
* ``RealizedSum`` is the image of a sum of chains only: it is parsed (at least two terms, none a :class:`ZeroMap`, one face and
  one domain face for all), and ``a + b`` of two realized faces builds it. ``sum(...)`` is compensated since Python 3.12 and
  :attr:`~RealizedSum.return_fraction` uses ``math.fsum`` on purpose (correctly rounded, as the declared wall's is).
* The ``inner`` operator is built on first use and cached; it is not part of the type's identity, and the type compares and hashes
  without ever building it. A pairing that does not cross (ERR-045) therefore raises at the first use of ``inner``, not at
  construction; the SN table's ``admit`` runs the same crossing check up front (``ordinate_pairing`` for a mirror, :func:`_crossing`
  for a translation wrap), so a declaration reaching ``realize`` through a problem has crossed already.

References
----------
* ``.claude/plans/wave2_boundary_carve.md`` sections 2.6 and 3 (laws L-A, L-E, L-L) and 6.2 (S2b gates 14 to 18, 20).
* ``docs/theory/foundations/boundary_conditions.rst`` (``bc-tensor-decomposition``, ``bc-method-realizability``).
"""

from __future__ import annotations

import math
import operator
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from functools import cached_property, reduce
from typing import TYPE_CHECKING, ClassVar, Optional

import numpy as np

from orpheus.numerics.elementary_operators import PermutationOperator
from orpheus.numerics.operator import (
    IdentityOperator,
    LinearOperator,
    MissingAdjoint,
    NotInvertible,
    ZeroMorphism,
    adjointable,
    invertible,
)
from orpheus.sn.boundary.angular import IsotropicEmissionOperator, PartialCurrentOperator
from orpheus.transport.boundary_conditions.network import Chain, DiffuseChain, SpecularChain
from orpheus.transport.boundary_conditions.realization import stamp_boundary_role

if TYPE_CHECKING:
    from orpheus.geometry.transformation import Permutation
    from orpheus.numerics.space import FunctionSpace
    from orpheus.numerics.spaces.angular_trace_space import AngularFaceTraceSpace, AngularTraceSpace

__all__ = [
    "RankOneReemission",
    "RealizedBoundary",
    "RealizedSum",
    "ScaledOrdinatePermutation",
    "ZeroMap",
]

#: The ordinate axis of every per-face half-trace space: LEADING (the half-trace spaces state it as "row order is the contract").
_ORDINATE_AXIS: int = 0


def _attenuated_kernel_operator(kernel: LinearOperator, alpha: float) -> LinearOperator:
    r"""``alpha * (kernel (x) K_g)``, stamped :attr:`~orpheus.numerics.operator.BlockRole.BOUNDARY`: the one shared body.

    The section-16A.10 decomposition of a boundary operator is ``B = G_patch (x) K_omega (x) K_g``; every leaf this module builds
    is that shape with the group factor an identity and the amplitude out front. The amplitude regimes are the algebra's, not
    branches: :math:`1 \cdot X` IS :math:`X` (the one-term collapse of the formal linear combination), so ``alpha = 1`` is the bare
    tensor product, and the role lives on the LEAF (the user's ruling, 2026-10-20): the kernel is stamped, and ``alpha * kernel``
    derives ``BOUNDARY`` by the join.
    """
    return float(alpha) * stamp_boundary_role(kernel & IdentityOperator())


def _parse_alpha(alpha: object, owner: str) -> float:
    r"""``alpha`` as a float in ``(0, 1]``: a zero scale is the zero map, above one is super-Markov, and NaN is not a scale."""
    if isinstance(alpha, bool) or not isinstance(alpha, (int, float, np.integer, np.floating)):
        raise TypeError(f"{owner}: alpha is a real scale in (0, 1], got {alpha!r}")
    value = float(alpha)
    if not math.isfinite(value) or not 0.0 < value <= 1.0:
        raise ValueError(
            f"{owner}: alpha must lie in (0, 1] (0 is the ZeroMap, above 1 is super-Markov, NaN and infinity are not scales), got {value!r}"
        )
    return value


def _crossing(permutation: "Permutation", outflow: "AngularFaceTraceSpace", inflow: "AngularFaceTraceSpace", owner: str) -> np.ndarray:
    r"""The local gather of a pairing :math:`\Gamma_+(f') \to \Gamma_-(f)`: row :math:`j` reads the ordinate :math:`\pi^{-1}(\mathrm{inflow}[j])`.

    The quotient equivariance :math:`\psi(g x, g\Omega) = \psi(x, \Omega)` evaluated on the trace gives
    :math:`\gamma_-\psi|_f(\Omega) = \gamma_+\psi|_{f'}(h^{-1}\Omega)`, with :math:`h` the deck element carrying the domain face onto
    the installation face. With :math:`\pi` the ordinate permutation :math:`h` induces (convention :math:`h(\Omega_i) = \Omega_{\pi(i)}`),
    row :math:`j` therefore gathers ordinate :math:`\pi^{-1}(\mathrm{inflow}[j])` at that ordinate's position INSIDE :math:`\Gamma_+(f')`,
    which is what the half-space's ``to_local`` computes (the space owns its row order; on a slab the mirror REVERSES order, so a
    hand-written ``arange`` would be wrong). The inverse is load-bearing only for a latent case (a sector rotation, where
    :math:`\pi \neq \pi^{-1}`); every shipped pairing is an involution or the identity.

    The ERR-045 invariant (a pairing maps inflow to OUTFLOW, never to itself) is enforced by ``to_local``: a gathered ordinate
    outside :math:`\Gamma_+(f')` raises. A pairing is a BIJECTION between the two half-traces, so two sets of different size have
    no realization.

    Raises
    ------
    ValueError
        The half-traces differ in size, or the pairing sends an inflow ordinate outside :math:`\Gamma_+` (ERR-045).
    """
    if inflow.shape[0] != outflow.shape[0]:
        raise ValueError(
            f"{owner}: |Gamma_-| = {inflow.shape[0]} but |Gamma_+| = {outflow.shape[0]}. A pairing is a BIJECTION between the two "
            f"half-traces, so a pair that differs in size has no realization (the quadrature's tangential band swallowed ordinates asymmetrically)."
        )
    try:
        return outflow.to_local(permutation.inverse().indices[inflow.ordinate_indices])
    except ValueError as error:
        raise ValueError(
            f"{owner}: the pairing does not map Gamma_-({inflow.name}) into Gamma_+({outflow.name}) (ERR-045: an inflow ordinate "
            f"whose partner is not outflow). Underlying: {error}"
        ) from error


def _tiers_meet(chain: type, outflow: "AngularFaceTraceSpace", inflow: "AngularFaceTraceSpace", owner: str) -> None:
    r"""The declared tiers of ``chain`` meet the tiers of the spaces the image derives from (a check, never a dispatch on ``role``)."""
    first, last = chain.factors[0].domain, chain.factors[-1].codomain  # type: ignore[attr-defined]
    if first != outflow.role or last != inflow.role:
        raise ValueError(
            f"{owner}: the chain {chain.__name__} runs {first!r} to {last!r} but the spaces it is realized on carry the tiers "
            f"{outflow.role!r} to {inflow.role!r}"
        )


class RealizedBoundary(LinearOperator, ABC):
    r"""A realized face operator :math:`\Gamma_+ \to \Gamma_-`, with the facts its consumers read as properties of its TYPE.

    The abstract base of the four images of the response network; replaces ``_BoundBoundaryOperator`` (a realized operator paired
    with the law it came from). The algebra object delegates: :meth:`apply`, :meth:`apply_transpose`, :attr:`domain`,
    :attr:`codomain`, :attr:`block_role` and the structural predicates are those of :attr:`inner`, the operator the type builds
    from its parameters. The ``face`` attribute (a dataclass field of the leaves, a derived property of the sum) is the face whose
    INFLOW the operator fills.

    Subclasses are frozen dataclasses (``eq=False``) holding their parameters only and defining ``__eq__``/``__hash__`` over
    them; the abstract members below are properties derived from those parameters.
    """

    face: str

    @property
    @abstractmethod
    def domain_face(self) -> str:
        """The face whose OUTFLOW half-trace this operator consumes (the partner for a translation wrap, else the face itself)."""

    @property
    @abstractmethod
    def return_fraction(self) -> float:
        r""":math:`J^-/J^+` of the realized map on this face's trace (the partial-current return ratio)."""

    @property
    @abstractmethod
    def permutes_ordinates(self) -> bool:
        """The operator is a scaled ordinate permutation that is not the identity (a mirror: True; a translation wrap: False)."""

    @property
    def is_vacuum(self) -> bool:
        """The realized map is the zero map, ``R = 0``: derived from the TYPE (the question is answered by the realized map)."""
        return isinstance(self, ZeroMap)

    @property
    @abstractmethod
    def inner(self) -> LinearOperator:
        """The operator, DERIVED from the parameters (never a constructor argument)."""

    # ── the algebra object delegates to the derived operator ─────────────────────────────────────────

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        return getattr(self.inner, "domain", None)

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return getattr(self.inner, "codomain", None)

    @property
    def is_invertible(self) -> bool:
        return self.inner.is_invertible

    @property
    def is_adjointable(self) -> bool:
        return self.inner.is_adjointable

    def inverse(self) -> LinearOperator:
        r"""The derived operator's inverse (a permutation image has one; the projecting images raise the inner operator's own guard)."""
        if not invertible(self.inner):
            raise NotInvertible(f"{type(self).__name__}.inverse: the realized {type(self.inner).__name__} is not invertible.")
        return self.inner.inverse()

    @property
    def block_role(self):  # type: ignore[override]
        return getattr(self.inner, "block_role", None)

    def apply(self, psi: np.ndarray) -> np.ndarray:
        return self.inner.apply(psi)

    def apply_transpose(self, psi: np.ndarray) -> np.ndarray:
        if not adjointable(self.inner):
            raise MissingAdjoint(f"{type(self).__name__}.apply_transpose: the realized {type(self.inner).__name__} is not adjointable.")
        return self.inner.apply_transpose(psi)

    def __add__(self, other: object):  # type: ignore[override]
        r"""The image of a sum of chains: ``a + b`` of two realized faces is a :class:`RealizedSum` (flattened, in order)."""
        if isinstance(other, RealizedBoundary):
            return RealizedSum(terms=(*_terms_of(self), *_terms_of(other)))
        return super().__add__(other)  # type: ignore[misc]

    # ── identity: the parameters, never the stored operator ──────────────────────────────────────────

    @abstractmethod
    def _key(self) -> tuple:
        """The parameters that are the identity of the type."""

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, RealizedBoundary):
            return NotImplemented
        return type(self) is type(other) and self._key() == other._key()

    def __hash__(self) -> int:
        return hash((type(self).__name__, self._key()))


def _terms_of(op: RealizedBoundary) -> tuple[RealizedBoundary, ...]:
    return op.terms if isinstance(op, RealizedSum) else (op,)


@dataclass(frozen=True, eq=False)
class ZeroMap(RealizedBoundary):
    r"""The image of the empty network, :math:`R = 0`: the zero map :math:`\Gamma_+(f) \to \Gamma_-(f)`.

    Vacuum's whole content is :math:`R = 0`; a prescribed inflow's linear part is zero too (its content is the affine :math:`q`,
    which travels the boundary-source channel); and any wall of zero weight is this object. A zero map between two DIFFERENT
    spaces must emit the zero of the space it lands in, not an echo of the one it came from: the bound
    :class:`~orpheus.numerics.operator.ZeroMorphism` derives both zeros from its declared pair.

    Parameters
    ----------
    face
        The face whose inflow it fills (and whose outflow it consumes).
    trace
        The :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` of the problem; the spaces are read off its ladder.
    """

    face: str
    trace: "AngularTraceSpace" = field(repr=False)

    @property
    def domain_face(self) -> str:
        return self.face

    @property
    def return_fraction(self) -> float:
        return 0.0

    @property
    def permutes_ordinates(self) -> bool:
        return False

    @cached_property
    def inner(self) -> LinearOperator:
        return stamp_boundary_role(ZeroMorphism(domain=self.trace.outflow_space(self.face), codomain=self.trace.inflow_space(self.face)))

    def _key(self) -> tuple:
        return (self.face,)


@dataclass(frozen=True, eq=False)
class ScaledOrdinatePermutation(RealizedBoundary):
    r"""The image of ``alpha * SpecularChain`` (and of a deck): :math:`R = \alpha\,P_\pi`, :math:`\Gamma_+(f') \to \Gamma_-(f)`.

    Row :math:`j` of the image reads the preimage of the :math:`j`-th inflow ordinate under :math:`\pi`, at that ordinate's
    position inside :math:`\Gamma_+(f')` (:func:`_crossing`). A symmetry has no amplitude, so a mirror DECK is this type at
    ``alpha = 1`` and a partial specular wall is the same matrix scaled by its weight.

    Parameters
    ----------
    alpha
        The scale, in ``(0, 1]``. It IS the return fraction: the permutation preserves the partial current (a measure-preserving
        bijection, certified by ``ordinate_pairing`` for a mirror), so :math:`J^-(R\phi) = \alpha J^+(\phi)`.
    permutation
        The ordinate permutation :math:`\pi` (``pi.indices[i]`` is the image of ordinate ``i``): the certified one a mirror pairing
        returns, or the identity for a translation wrap.
    face
        The face whose inflow it fills.
    partner
        The face whose outflow it consumes: ``face`` for a mirror, ``face_opposite(face)`` for a translation wrap.
    trace
        The problem's :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`.

    Raises
    ------
    ValueError
        At construction: ``alpha`` outside ``(0, 1]``, or the spaces do not carry the tiers of the specular chain. At first use of
        :attr:`inner`: the half-traces differ in size or the pairing does not map inflow into outflow (ERR-045).
    """

    alpha: float
    permutation: "Permutation" = field(repr=False)
    face: str
    partner: str
    trace: "AngularTraceSpace" = field(repr=False)

    chain: ClassVar[type] = SpecularChain

    def __post_init__(self) -> None:
        object.__setattr__(self, "alpha", _parse_alpha(self.alpha, "ScaledOrdinatePermutation"))
        _tiers_meet(self.chain, self.trace.outflow_space(self.partner), self.trace.inflow_space(self.face), "ScaledOrdinatePermutation")

    @property
    def domain_face(self) -> str:
        return self.partner

    @property
    def return_fraction(self) -> float:
        return self.alpha

    @property
    def permutes_ordinates(self) -> bool:
        return not np.array_equal(self.permutation.indices, np.arange(self.permutation.n))

    @cached_property
    def inner(self) -> LinearOperator:
        outflow, inflow = self.trace.outflow_space(self.partner), self.trace.inflow_space(self.face)
        local = _crossing(self.permutation, outflow, inflow, f"ScaledOrdinatePermutation({self.partner!r} -> {self.face!r})")
        return _attenuated_kernel_operator(PermutationOperator(local, axis=_ORDINATE_AXIS, domain=outflow, codomain=inflow), self.alpha)

    def _key(self) -> tuple:
        return (self.face, self.partner, self.alpha, self.permutation.indices.tobytes())


@dataclass(frozen=True, eq=False)
class RankOneReemission(RealizedBoundary):
    r"""The image of ``alpha * DiffuseChain``: :math:`R = \alpha\,(E \circ R_t)`, the factored Lambertian :math:`\Gamma_+(f) \to \Gamma_-(f)`.

    The two-link chain :math:`\Gamma_+ \xrightarrow{R_t} S \xrightarrow{E} \Gamma_-`: the contraction against the trace measure
    :math:`w|\Omega\cdot n_f|` (the outgoing partial current :math:`J^+`, a named intermediate rather than an anonymous local) and
    the uniform re-emission divided by the INFLOW mass. Neither the norm, nor the contraction weights, nor the inflow extent is a
    parameter: all three are READ from the half-spaces of the trace ladder (the ladder's metric is the retraction's covector), so a
    free ``norm`` that disagrees with its half-space is unspellable. The hemisphere IS ``trace.outflow_space(face)``, derived from
    the face name once, so the orientation cross-check of the retired builder has nothing left to compare.

    Parameters
    ----------
    alpha
        The scale, in ``(0, 1]``: the return fraction, because ``Rt_in o E = id`` on the inflow rung.
    face
        The face whose inflow it fills (and whose outflow it consumes).
    trace
        The problem's :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`.
    """

    alpha: float
    face: str
    trace: "AngularTraceSpace" = field(repr=False)

    chain: ClassVar[type] = DiffuseChain

    def __post_init__(self) -> None:
        object.__setattr__(self, "alpha", _parse_alpha(self.alpha, "RankOneReemission"))
        _tiers_meet(self.chain, self.trace.outflow_space(self.face), self.trace.inflow_space(self.face), "RankOneReemission")
        if self.trace.inflow_space(self.face).shape[0] == 0:
            raise ValueError(f"RankOneReemission: face {self.face!r} has no inflow ordinates, so the re-emitted flux has no Gamma_- to land on")

    @property
    def domain_face(self) -> str:
        return self.face

    @property
    def return_fraction(self) -> float:
        return self.alpha

    @property
    def permutes_ordinates(self) -> bool:
        return False

    @cached_property
    def inner(self) -> LinearOperator:
        outflow, current, inflow = self.trace.outflow_space(self.face), self.trace.current_space(self.face), self.trace.inflow_space(self.face)
        contraction = PartialCurrentOperator(outflow.inner_product_weights, domain=outflow, codomain=current)
        emission = IsotropicEmissionOperator(float(np.sum(inflow.inner_product_weights)), int(inflow.shape[0]), domain=current, codomain=inflow)
        return _attenuated_kernel_operator(emission @ contraction, self.alpha)

    def _key(self) -> tuple:
        return (self.face, self.alpha)


@dataclass(frozen=True, eq=False)
class RealizedSum(RealizedBoundary):
    r"""The image of a sum of chains: the sum of the images, :math:`R = \sum_k R_k` (a mixed wall).

    PARSED: at least two terms (one term is the term), no :class:`ZeroMap` among them (a zero term is the empty network's image
    and is never summed), and one face and one domain face for all (a wall's chains all fill the same inflow from the same
    outflow). Built by ``a + b`` of realized faces, in the canonical chain order of the declared wall.

    Parameters
    ----------
    terms
        The images of the chains, in order.
    """

    terms: tuple[RealizedBoundary, ...]

    def __post_init__(self) -> None:
        terms = tuple(self.terms)
        if len(terms) < 2:
            raise ValueError(f"RealizedSum: a sum of chains has at least two terms, got {len(terms)}")
        if not all(isinstance(term, RealizedBoundary) for term in terms):
            raise TypeError(f"RealizedSum: every term is a RealizedBoundary, got {[type(term).__name__ for term in terms]}")
        if any(isinstance(term, ZeroMap) for term in terms):
            raise ValueError("RealizedSum: a ZeroMap term is the image of the empty network and is never summed")
        if len({term.face for term in terms}) != 1 or len({term.domain_face for term in terms}) != 1:
            raise ValueError(f"RealizedSum: the terms fill different faces or consume different outflows: {[(t.face, t.domain_face) for t in terms]}")
        object.__setattr__(self, "terms", terms)

    @property
    def face(self) -> str:  # type: ignore[override]
        return self.terms[0].face

    @property
    def domain_face(self) -> str:
        return self.terms[0].domain_face

    @property
    def return_fraction(self) -> float:
        return math.fsum(term.return_fraction for term in self.terms)

    @property
    def permutes_ordinates(self) -> bool:
        return False

    @cached_property
    def inner(self) -> LinearOperator:
        return reduce(operator.add, (term.inner for term in self.terms))

    def _key(self) -> tuple:
        return self.terms
