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

**The facts are properties of the type.** The questions the SN consumers ask of a realized face (does it return
everything, :attr:`~RealizedBoundary.return_fraction`; does it close its axis into a mirror loop, :attr:`~RealizedBoundary.closes_axis`; does it relabel ordinates, :attr:`~RealizedBoundary.permutes_ordinates`;
whose outflow does it consume, :attr:`~RealizedBoundary.domain_face`; is it the zero map, :attr:`~RealizedBoundary.is_vacuum`;
and what operator is it, :attr:`~RealizedBoundary.inner`) used to be answered by asking the LAW the operator came from
(``law_permutes_ordinates``, ``response_kernel.is_zero``, ``.amplitude``, ``.geometry_map.domain_face``, a ``kind`` string). The
realizer that builds the operator already knows each fact, so each is DERIVED from the type's parameters, never stored beside
them. A type holds ONLY its parameters ``(alpha, ..., face, trace)``; the operator is built from spaces READ OFF the trace
ladder (``trace.outflow_space(f)``, ``trace.current_space(f)``, ``trace.inflow_space(f)``), so a lying ``alpha`` is unspellable
and ``dataclasses.replace(op, alpha=...)`` REBUILDS the operator rather than editing a fact. The types are frozen and compare
and hash by their CONTENT: the type, the face, the partner, ``alpha``, the permutation's index bytes AND the content of the two
half-trace spaces the operator is derived from (their ordinate indices and metric weights), not by a stored operator (a
``LinearOperator`` has no content identity) and not by the trace spaces' own equality (a trace space compares by ``(name, shape)``
only, so two traces of different quadratures with the same total size are equal as spaces yet carry different operators). Two
operators on different discretizations are therefore unequal and hash apart. The old shim's ``str`` equality arm (the retired
``op == "reflective"``, #553) is not carried.

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

**Rider (b) of the WF2.2 record, decided (S2b amend): no ``OrdinatePairing`` carrier, because the type checks the facts itself.**
The plan asked whether ``ScaledOrdinatePermutation`` should take a frozen ``OrdinatePairing(permutation, motion, axis)`` whose only
constructor is ``ordinate_pairing``, so that an uncertified ``quadrature.ordinate_permutation(motion)`` could not be realized
(Pattern 4). The reasons against, in the order they were found. (1) A certified and a raw permutation of the same closed rule are
the same array, so the carrier would certify a PROVENANCE and not a property; what must hold is a property of the pair
(permutation, trace). (2) The properties the facts of this type rest on are exactly two, and the constructor now checks both on the
trace it is given: the pairing crosses (ERR-045, :func:`_crossing`) and it carries the outflow trace measure onto the inflow trace
measure, which is what makes ``return_fraction == alpha``. ``[M]`` before the amend a mirror composed with a shift constructed,
reported ``return_fraction`` 1.0 and measured 1.198 on a non-uniform outflow; it is refused at construction now
(``test_a_pairing_that_does_not_cross_or_moves_measure_is_refused_at_construction``). The third check of ``ordinate_pairing``
(involution) is irrelevant to the facts: the gather uses the inverse. (3) The translation wrap's identity permutation has no
constructor under the proposed rule, which would need a second one, and Python cannot make a frozen dataclass's constructor private,
so "only constructor" would be a convention (X3 of the instrument doctrine). The uncertified permutation that survives the two
checks is one that is correct on this trace, which is the claim the type makes. Trigger to revisit: a realizer whose pairing is not
read off a trace ladder (a second consumer in CP, MoC or MC), or #657.

**Rider (c): the tiers are checked once, by the declared network.** The WF2.2 record asked that the line "the realization functor
checks the tiers" land only with the check. The check is :func:`orpheus.transport.boundary_conditions.network.chain_composes`,
run at import over every member of the closed sum ``Chain`` (``factors[0].domain == "outflow"``, consecutive factors meet,
``factors[-1].codomain == "inflow"``), and the trace ladder's rungs are constants (``trace.outflow_space(f)`` is always the outflow
tier, ``trace.inflow_space(f)`` always the inflow tier). A second check inside the leaf constructors could not fail on any
reachable input (``[M]`` S2b review: turning it into a no-op reddened 0 of 667 tests), so it was deleted rather than kept as a
second spelling of ``chain_composes``; the attribute :attr:`ScaledOrdinatePermutation.chain` names the chain of the closed sum a
leaf realizes, and the declared and the realized composition meet because both read the ladder's two fixed rungs.

Gotchas
=======

* ``apply`` takes the face's outflow half-trace of the DOMAIN face (:attr:`~RealizedBoundary.domain_face`), which is the
  partner for a translation wrap and the face itself otherwise; ``block_role`` is :attr:`~orpheus.numerics.operator.BlockRole.BOUNDARY`.
* A zero scale is the :class:`ZeroMap`, refused as an ``alpha`` here (``alpha`` lies in ``(0, 1]``); a super-Markov or non-finite
  scale is refused at construction.
* ``return_fraction == alpha`` holds only if the permutation preserves the trace measure, so :class:`ScaledOrdinatePermutation`
  checks, AT CONSTRUCTION, that the pairing crosses (ERR-045) and that the inflow measure equals the gathered outflow measure
  (relative tolerance ``1e-12``, the certifier's). An uncertified permutation that does cross but does not preserve the measure
  is refused instead of reporting a false return fraction.
* ``RealizedSum`` is the image of a sum of chains only: it is parsed (at least two leaf terms, none a :class:`ZeroMap`, one face,
  one domain face and one trace for all, no chain type twice, ``math.fsum`` of the return fractions at most 1: the same interval
  as :func:`_parse_alpha`), and ``a + b`` of two realized faces builds it. The zero map is the additive identity: ``ZeroMap + X``
  is ``X`` and a sum that is left with one term collapses to it. :attr:`~RealizedSum.return_fraction` uses ``math.fsum`` on
  purpose (correctly rounded, as the declared wall's is).
* The ``inner`` operator is built on first use and cached; it is not part of the type's identity. The pairing's local gather is
  computed (and cached) at construction, so a pairing that does not cross never constructs.

References
----------
* ``.claude/plans/wave2_boundary_carve.md`` sections 2.6 and 3 (laws L-A, L-E, L-L) and 6.2 (S2b gates 14 to 18, 20).
* ``docs/theory/foundations/boundary_conditions.rst`` (:ref:`bc-architecture-of-record`; ``bc-tensor-decomposition``, ``bc-method-realizability``).
"""

from __future__ import annotations

import math
import operator
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from functools import cached_property, reduce
from typing import TYPE_CHECKING, ClassVar, Optional, cast

import numpy as np

from orpheus.numerics.elementary_operators import PermutationOperator
from orpheus.numerics.quadrature.ordinate_pairing import MEASURE_RTOL
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
from orpheus.values.scalars import parse_finite_real

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
    r"""``alpha`` as a float in ``(0, 1]``: a zero scale is the zero map, above one is super-Markov, and NaN is not a scale.

    The real, finite scalar is the one parser's (:func:`~orpheus.values.scalars.parse_finite_real`: a bool or a non-real is a
    ``TypeError``, NaN and infinity are keyed ``ValueError`` refusals, owned once); only the interval is this type's.
    """
    value = parse_finite_real(alpha, f"{owner}.alpha")
    if not 0.0 < value <= 1.0:
        raise ValueError(f"{owner}: alpha must lie in (0, 1] (0 is the ZeroMap, above 1 is super-Markov), got {value!r}")
    return value


def _space_content(space: "AngularFaceTraceSpace") -> tuple:
    r"""The content of a half-trace space: its name, shape, ordinate indices and metric weights, as bytes (a space's own ``==`` sees only ``(name, shape)``)."""
    weights = space.inner_product_weights
    return (space.name, space.shape, space.ordinate_indices.tobytes(), b"" if weights is None else np.asarray(weights).tobytes())


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
    def closes_axis(self) -> bool:
        r"""The realized map is a mirror loop on its own face: an ordinate permutation (not the identity) that returns everything.

        Derived from the REALIZED operator, never from how the wall was declared: a wall tagged with a unit specular condition and a
        mirror deck realize the same :class:`ScaledOrdinatePermutation` at ``alpha = 1``, hence the same operator, hence the same
        kernel of the loss operator. The gauge freedom of the loss operator and the proven mirror row of the diffusion-synthetic
        low-order system are both questions about this one fact. A face that permutes ordinates onto ANOTHER face's outflow (none is
        constructed today) or that returns only a fraction does not close its axis.
        """
        return self.permutes_ordinates and self.domain_face == self.face and self.return_fraction == 1.0

    @property
    def is_vacuum(self) -> bool:
        """The realized map is the zero map, ``R = 0``: derived from the TYPE (the question is answered by the realized map)."""
        return isinstance(self, ZeroMap)

    @cached_property
    def inner(self) -> LinearOperator:
        """The operator, DERIVED from the parameters (never a constructor argument) and built once, at first use."""
        return self._derive_inner()

    @abstractmethod
    def _derive_inner(self) -> LinearOperator:
        """The one body of the derived operator: each realized type states it, :attr:`inner` caches it."""

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
        r"""The image of a sum of chains: ``a + b`` of two realized faces is a :class:`RealizedSum` (flattened, in order).

        The zero map is the identity: its terms drop out, a sum left with one term IS that term, and a sum of zero maps is the zero map.
        Every operand, a zero map included, fills the same face on the same inflow space: a zero map on another face or another
        discretization is the same refusal as a sum across faces.
        """
        if not isinstance(other, RealizedBoundary):
            return super().__add__(other)  # type: ignore[misc]
        operands = (*_leaves(self), *_leaves(other))
        if len({term.face for term in operands}) != 1:
            raise ValueError(f"RealizedBoundary.__add__: the terms fill different faces: {[term.face for term in operands]}")
        if len({_space_content(term.trace.inflow_space(term.face)) for term in operands}) != 1:
            raise ValueError("RealizedBoundary.__add__: the terms land on different inflow spaces (different discretizations of the face)")
        survivors = tuple(term for term in operands if not isinstance(term, ZeroMap))
        if not survivors:
            return self
        return survivors[0] if len(survivors) == 1 else RealizedSum(terms=survivors)

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


def _leaves(op: RealizedBoundary) -> "tuple[ZeroMap | ScaledOrdinatePermutation | RankOneReemission, ...]":
    r"""The leaf images of ``op``: a :class:`RealizedSum`'s terms, or ``op`` itself (the four images are the closed family)."""
    return cast("tuple[ZeroMap | ScaledOrdinatePermutation | RankOneReemission, ...]", op.terms if isinstance(op, RealizedSum) else (op,))


def _spaces_key(op: "ZeroMap | ScaledOrdinatePermutation | RankOneReemission") -> tuple:
    r"""The content of the two half-trace spaces a leaf operator is derived from (its domain and codomain)."""
    return (_space_content(op.trace.outflow_space(op.domain_face)), _space_content(op.trace.inflow_space(op.face)))


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

    def __post_init__(self) -> None:
        # Both half-traces are read at construction, so a face the trace does not hold is refused here (ValueError "Unknown face"),
        # as the sibling ScaledOrdinatePermutation refuses through its gather, and not at the first ``.domain`` read.
        self.trace.outflow_space(self.face)
        self.trace.inflow_space(self.face)

    @property
    def domain_face(self) -> str:
        return self.face

    @property
    def return_fraction(self) -> float:
        return 0.0

    @property
    def permutes_ordinates(self) -> bool:
        return False

    def _derive_inner(self) -> LinearOperator:
        return stamp_boundary_role(ZeroMorphism(domain=self.trace.outflow_space(self.face), codomain=self.trace.inflow_space(self.face)))

    def _key(self) -> tuple:
        return (self.face, _spaces_key(self))


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
        At construction: ``alpha`` outside ``(0, 1]``, the half-traces differ in size, the pairing does not map inflow into outflow
        (ERR-045), or the pairing does not carry the outflow measure onto the inflow measure.
    """

    alpha: float
    permutation: "Permutation" = field(repr=False)
    face: str
    partner: str
    trace: "AngularTraceSpace" = field(repr=False)

    chain: ClassVar[type[Chain]] = SpecularChain

    def __post_init__(self) -> None:
        object.__setattr__(self, "alpha", _parse_alpha(self.alpha, "ScaledOrdinatePermutation"))
        self._gather  # the pairing crosses and preserves the measure, or the type does not exist

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
    def _gather(self) -> np.ndarray:
        r"""The local gather of the pairing, computed once at construction: it crosses (ERR-045) and preserves the trace measure.

        ``return_fraction == alpha`` rests on the second check: the permutation carries the outflow measure onto the inflow measure
        (row :math:`j` reads an ordinate whose metric equals the inflow's :math:`j`-th), so :math:`J^-(R\phi) = \alpha J^+(\phi)`
        exactly. An uncertified permutation that crosses but moves measure (a mirror composed with a shift) reports a false return
        fraction and is refused here.
        """
        owner = f"ScaledOrdinatePermutation({self.partner!r} -> {self.face!r})"
        outflow, inflow = self.trace.outflow_space(self.partner), self.trace.inflow_space(self.face)
        local = _crossing(self.permutation, outflow, inflow, owner)
        if outflow.inner_product_weights is not None and inflow.inner_product_weights is not None:
            carried, wanted = np.asarray(outflow.inner_product_weights)[local], np.asarray(inflow.inner_product_weights)
            if not np.allclose(carried, wanted, rtol=MEASURE_RTOL, atol=0.0):
                raise ValueError(
                    f"{owner}: the pairing does not carry the outflow trace measure onto the inflow trace measure (worst relative gap "
                    f"{float(np.max(np.abs(carried - wanted) / np.abs(wanted))):.3e}), so its return fraction is not alpha = {self.alpha}."
                )
        return local

    def _derive_inner(self) -> LinearOperator:
        outflow, inflow = self.trace.outflow_space(self.partner), self.trace.inflow_space(self.face)
        return _attenuated_kernel_operator(PermutationOperator(self._gather, axis=_ORDINATE_AXIS, domain=outflow, codomain=inflow), self.alpha)

    def _key(self) -> tuple:
        return (self.face, self.partner, self.alpha, self.permutation.indices.tobytes(), _spaces_key(self))


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

    chain: ClassVar[type[Chain]] = DiffuseChain

    def __post_init__(self) -> None:
        object.__setattr__(self, "alpha", _parse_alpha(self.alpha, "RankOneReemission"))
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

    def _derive_inner(self) -> LinearOperator:
        outflow, current, inflow = self.trace.outflow_space(self.face), self.trace.current_space(self.face), self.trace.inflow_space(self.face)
        out_weights, in_weights = outflow.inner_product_weights, inflow.inner_product_weights
        assert out_weights is not None and in_weights is not None  # narrowing only: a half-trace built from a quadrature always carries its metric
        contraction = PartialCurrentOperator(out_weights, domain=outflow, codomain=current)
        emission = IsotropicEmissionOperator(float(np.sum(in_weights)), int(inflow.shape[0]), domain=current, codomain=inflow)
        return _attenuated_kernel_operator(emission @ contraction, self.alpha)

    def _key(self) -> tuple:
        return (self.face, self.alpha, _spaces_key(self))


@dataclass(frozen=True, eq=False)
class RealizedSum(RealizedBoundary):
    r"""The image of a sum of chains: the sum of the images, :math:`R = \sum_k R_k` (a mixed wall).

    PARSED: at least two leaf terms (one term is the term, a zero term is the identity and a sum is flattened: ``a + b`` does all
    three), one face, one domain face and one trace for all (a wall's chains all fill the same inflow from the same outflow on one
    discretization), each chain type at most once, and the return fractions summing to at most 1 (the declared wall's cone). Built by ``a + b`` of realized faces, in the canonical chain order of the declared wall.

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
        if not all(isinstance(term, RealizedBoundary) and not isinstance(term, RealizedSum) for term in terms):
            raise TypeError(f"RealizedSum: every term is a leaf image of a chain (a sum is flattened by ``a + b``), got {[type(term).__name__ for term in terms]}")
        if any(isinstance(term, ZeroMap) for term in terms):
            raise ValueError("RealizedSum: a ZeroMap term is the image of the empty network, the additive identity, and is never summed (``a + b`` drops it)")
        leaves = cast("tuple[ScaledOrdinatePermutation | RankOneReemission, ...]", terms)
        if len({term.face for term in leaves}) != 1 or len({term.domain_face for term in leaves}) != 1:
            raise ValueError(f"RealizedSum: the terms fill different faces or consume different outflows: {[(t.face, t.domain_face) for t in leaves]}")
        if len({_spaces_key(term) for term in leaves}) != 1:
            raise ValueError("RealizedSum: the terms are realized on different trace spaces (different discretizations of the face)")
        if len({type(term) for term in leaves}) != len(leaves):
            raise ValueError(f"RealizedSum: a chain appears twice among {[type(term).__name__ for term in leaves]}; a wall has one weight per chain")
        total = math.fsum(term.return_fraction for term in leaves)
        if total > 1.0:
            raise ValueError(f"RealizedSum: the return fractions sum to {total!r} > 1 (super-Markov), outside the same interval as alpha")
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

    def _derive_inner(self) -> LinearOperator:
        return reduce(operator.add, (term.inner for term in self.terms))

    def _key(self) -> tuple:
        return self.terms
