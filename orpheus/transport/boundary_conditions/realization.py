r"""The per-method realizer tables: admission and realization through ONE door (``realization.py`` of the conditions package).

**The mathematics.** Realization is a partial map from declarations to operators on a method's trace space
(:mod:`.admission`); for the shared family it is a FUNCTOR on the response network: the image of the empty network is the
method's zero map, the image of one chain is that chain's realizer's, and the image of a sum is the sum of the images
(summed in ``CHAIN_ORDER``, so the realization does not depend on how the wall was built). The affine part :math:`q`
realizes separately, into the boundary-source channel, as an array.

**The table.** A :class:`RealizerTable` belongs to ONE method and is keyed by TYPE, filled by realizers:

* the chain types :class:`~.network.SpecularChain` and :class:`~.network.DiffuseChain` (exactly the members of the closed
  sum ``Chain``), the marker :class:`~.network.ZeroNetwork` (the key of the zero image), the marker
  :class:`~.inflow.InflowTerms` (the source channel; its entry is a :class:`SourceRealizer`, the others
  :class:`ResponseRealizer`) and the method-only condition types (a ``ZeroFluxCondition``);
* ``WallCondition`` itself is NOT a key: it is realized generically by the functor, so a method never writes a branch
  over the shared family's cone point;
* decks have exactly one :class:`DeckRealizer` per table (an exhaustive match over ``DeckTransformation``).

**Registration is a class keyword plus ``seal()``.** ``class SNSpecularChain(ResponseRealizer, table=SN_REALIZERS,
realizes=SpecularChain)`` registers at class creation and refuses a duplicate, a sealed table and a key that is not a
legal unit. The table and every realizer that registers into it live in ONE module per method, and that module calls
``table.seal()`` on its last line: anyone holding the table has executed every registration (this answers the
registration-timing hazard of an empty registry in a fresh process). An UNSEALED table refuses ``admit``, ``realize``,
``require`` and ``catalogue`` with :class:`TableNotSealedError`; registering into a sealed one raises
:class:`SealedTableError`. A new response for one method touches only that method's realizer module.

**Exact-type lookup, never the MRO.** A subclass is a new declaration and owes its own entry: a method-only condition of
another method must not be realized silently as its parent.

**One door.** :meth:`RealizerTable.realize` calls :meth:`~RealizerTable.admit` (through :meth:`~RealizerTable.require`)
and raises :class:`~.errors.BoundaryRefused` iff the admission is :class:`~.admission.Refused`; no realizer body raises on
its own for an admission reason, so the query and the realization cannot disagree by construction (law L-D). The one
declared second door is POSING: the eigenvalue entries refuse a boundary source before any realization, a refusal of the
QUESTION and not of the declaration.

**The unregistered door.** A unit with no entry is ``Refused(STATE_CONE, ...)`` when its ``return_fraction < 0`` and the
table's :class:`StateCone` is ``NONNEGATIVE_FLUX`` (so SN refuses diffusion's zero flux on the state-cone axis without
importing it: SN may not import diffusion, an L3 to L3 edge), and ``Refused(PLUMBING, ..., table.unregistered_issue)``
otherwise, naming the table's admitted keys. PLUMBING owes an issue, and which issue tracks a method's missing entries is
the method's to say: the table is built with a required ``unregistered_issue`` (CP, MoC and MC cite #654, whose tables
keep the old string registries' coverage; SN and diffusion cite #659, whose tables register every shared unit and meet
this door only with another method's own condition; the orchestrator's ruling, 2026-10-23).

**Two spaces, two type parameters.** ``RealizerTable[SpaceT, OpT, SourceSpaceT]``: ``realize`` takes the METHOD space of the face
(``SpaceT``; ``SNMethodSpace`` for SN) and returns the method's operator type ``OpT``, while ``sample_inflow`` takes the trace's
inflow half-space (``SourceSpaceT``; SN passes ``trace.inflow_space(face)``), the space on which a source :math:`q` lives. They are
different spaces, so one parameter could not type both and a typed call of the source verb needed a ``type: ignore``.

**Sampling.** The ONE source verb is :meth:`RealizerTable.sample_inflow`: ``require`` the source, then the registered
:meth:`SourceRealizer.sample`, whose default body is :func:`sample_terms`: ``q = uniform * ones + sum_i w_i sampler_i.sample(space)``.

**Gotchas.**

* The sampler contract is exact: ``sample_terms`` raises ``ValueError`` when a sampler returns anything but
  ``space.shape`` (a wrong-shape array would otherwise be broadcast silently).
* Lookup is by EXACT type everywhere, the functor's branch included: a ``WallCondition`` subclass is neither the
  functor's domain nor registrable (``_register_response`` refuses it), so it is admitted as an unregistered unit.
* The source is met for EVERY condition type with a nonzero source (the base contract), not only the shared wall.
* ``admit`` takes a declaration (a deck, a ``BoundaryCondition`` or an ``InflowTerms``) and raises ``TypeError`` for
  anything else, a ``WallDisplacement`` included.
* The `weight` a realizer receives is the chain's weight in the network; it is ``1.0`` for the ``ZeroNetwork`` key (which
  carries no scalar) and for a method-only condition key.
* A pure inflow (empty network, nonzero source) consults the ``ZeroNetwork`` realizer: its linear image is the zero map.
* ``realize_direction`` is NOT minted (no caller until #529's first signed-sum consumer).

Layer: ``transport`` (L2). Imports ``numerics`` (the ``BlockRole`` stamp), ``geometry`` (the deck vocabulary) and its
siblings.
"""

from __future__ import annotations

import operator
from abc import ABC, abstractmethod
from collections.abc import Mapping
from enum import Enum
from functools import reduce
from typing import TYPE_CHECKING, ClassVar, Generic, TypeVar, cast, get_args, overload

import numpy as np

from orpheus.geometry.deck import DeckTransformation, FacePairing, pairing_of
from orpheus.numerics.operator import BlockRole
from orpheus.transport.boundary_conditions.admission import Admission, FaceSite, Realized, RefusalAxis, Refused
from orpheus.transport.boundary_conditions.condition import BoundaryCondition, WallCondition
from orpheus.transport.boundary_conditions.errors import BoundaryError, BoundaryRefused
from orpheus.transport.boundary_conditions.inflow import InflowTerms
from orpheus.transport.boundary_conditions.network import Chain, ZeroNetwork

if TYPE_CHECKING:
    from orpheus.numerics.operator import LinearOperator
    from orpheus.numerics.space import FunctionSpace

SpaceT = TypeVar("SpaceT")
OpT = TypeVar("OpT")
#: The space the source verb samples on: the INFLOW half-space :math:`\Gamma_-` of a face (SN passes ``trace.inflow_space(face)``),
#: a different space from the method space ``SpaceT`` the response verb realizes on, so it is a type parameter of its own
#: (hand-over (a) of slice S2b, closed in slice S3).
SourceSpaceT = TypeVar("SourceSpaceT", bound="FunctionSpace")
K = TypeVar("K")

__all__ = [
    "DeckRealizer",
    "RealizerTable",
    "RefusingRealizer",
    "ResponseRealizer",
    "SealedTableError",
    "SourceRealizer",
    "StateCone",
    "TableNotSealedError",
    "sample_terms",
    "stamp_boundary_role",
]

class TableNotSealedError(BoundaryError):
    """A query reached a table whose defining module has not called ``seal()`` (registrations may be missing)."""


class SealedTableError(BoundaryError):
    """A registration reached a sealed table."""


class StateCone(Enum):
    """The cone a method's trace state lives in (an Enum, not a bool that a forgetful method gets for free)."""

    NONNEGATIVE_FLUX = "psi >= 0 on the trace"  # SN, CP, MoC, MC
    SIGNED_CURRENT = "a signed partial current"  # diffusion


def sample_terms(terms: InflowTerms, space: FunctionSpace) -> np.ndarray:
    r""":math:`q` on ``space``: ``uniform * ones + sum_i w_i sampler_i.sample(space)`` (the default body of the source verb)."""
    q = terms.uniform * np.ones(space.shape)
    for weight, profile in terms.sampled:
        sample = np.asarray(profile.sampler.sample(space))
        if sample.shape != space.shape:
            raise ValueError(f"{type(profile.sampler).__name__}.sample returned shape {sample.shape}, not exactly space.shape {space.shape}")
        q = q + weight * sample
    return q


def _units(declared: BoundaryCondition) -> tuple[tuple[type, object, float], ...]:
    """The ``(key, unit, weight)`` triples a condition decomposes into: the ONE decomposition ``admit`` (the meet) and ``realize`` (the sum) share.

    The exact shared type ``WallCondition`` is the functor's domain: one triple per chain in ``CHAIN_ORDER`` with its
    weight, or the ``ZeroNetwork`` marker (weight one) when the network is empty, whatever the source. Any other
    condition is one unit under its own exact type (a subclass of ``WallCondition`` is therefore a unit of its own,
    never the functor's).
    """
    if type(declared) is WallCondition:
        if not declared.network.weights:
            return ((ZeroNetwork, ZeroNetwork(), 1.0),)
        return tuple((type(chain), chain, weight) for weight, chain in declared.network.terms)
    return ((type(declared), declared, 1.0),)


def stamp_boundary_role(op: LinearOperator) -> LinearOperator:
    r"""Stamp a realized boundary law with the :attr:`BlockRole.BOUNDARY` role.

    A COPY of ``orpheus.transport.boundary_realization.stamp_boundary_role`` (the old module is left for WF2.4, which
    deletes it; this copy then stands alone). Edit both or neither until then.

    A realized law is a boundary-block leaf (``A_ss`` only) on the
    direct-sum transport state ``V = V_bulk ⊕ V_boundary`` — it maps the
    outflow trace to the inflow trace with no bulk action (Issue #208 /
    Wave O). The role lives on the INSTANCE (the realized op is a generic
    numerics primitive — :class:`TensorProductOperator`,
    :class:`OperatorSum`, :class:`IdentityOperator`, … — that plays
    the BOUNDARY role only in this realization context), so it is stamped
    at the producer sites — the functional method realizers — through this
    ONE helper (``coding-elegance`` Pattern 7 / Pattern 2: the stamp is a
    shared realizer-layer concept, defined once).

    **Every realizable law is stamped, with no exception** (P3,
    2026-08-05). The affine split survives (:ref:`bc-affine-source-channel`);
    what is realized is the law ``γ₋ψ = L γ₊ψ + q``: a realizer realizes
    ``L``, and prescribed inflow's ``L`` is the zero morphism — which is
    linear, and is the same object vacuum realizes to. ``q`` reaches the
    solve through the boundary-source channel.

    ⚠ **The absence of a stamp never excluded a leaf from anything.** The
    boundary operator collects every face's law with no ``block_role``
    filter, so an unstamped operator entered ``B`` regardless and delivered
    its source there. Read this stamp as a *classification* that consumers
    may query, never as a gate that keeps something out.
    """
    op.block_role = BlockRole.BOUNDARY
    return op


class RealizerTable(Generic[SpaceT, OpT, SourceSpaceT]):
    """The realizers of ONE method, keyed by unit type; registered by class keyword, then sealed (see the module docstring)."""

    def __init__(self, method: str, *, state_cone: StateCone, unregistered_issue: int) -> None:
        if not isinstance(unregistered_issue, int) or isinstance(unregistered_issue, bool):
            raise TypeError(f"{method}: unregistered_issue is the int of the issue tracking the units this table does not register, got {unregistered_issue!r}")
        self.method = method
        self.state_cone = state_cone
        self.unregistered_issue = unregistered_issue
        self._responses: dict[type, ResponseRealizer] = {}
        self._source: SourceRealizer[SourceSpaceT] | None = None
        self._deck: DeckRealizer | None = None
        self._sealed = False

    # registration ------------------------------------------------------------------------------------------------

    def _open(self) -> None:
        if self._sealed:
            raise SealedTableError(f"{self.method}: the table is sealed; register before the defining module's last line")

    def _register_response(self, key: type, realizer: ResponseRealizer) -> None:
        self._open()
        legal = key in get_args(Chain) or key is ZeroNetwork or (issubclass(key, BoundaryCondition) and not issubclass(key, WallCondition))
        if not legal:
            raise TypeError(f"{self.method}: {key!r} is not a unit a response realizer registers under (a chain, ZeroNetwork or a method-only BoundaryCondition)")
        if key in self._responses:
            raise ValueError(f"{self.method}: {key.__name__} is registered twice")
        self._responses[key] = realizer

    def _register_source(self, realizer: SourceRealizer[SourceSpaceT]) -> None:
        self._open()
        if self._source is not None:
            raise ValueError(f"{self.method}: a source realizer is registered twice")
        self._source = realizer

    def _register_deck(self, realizer: DeckRealizer) -> None:
        self._open()
        if self._deck is not None:
            raise ValueError(f"{self.method}: a deck realizer is registered twice")
        self._deck = realizer

    def seal(self) -> None:
        """Close the table: called on the last line of the defining module; later registration raises."""
        self._sealed = True

    def _require_sealed(self) -> None:
        if not self._sealed:
            raise TableNotSealedError(f"{self.method}: the table is not sealed; its defining module has not finished registering")

    # queries -----------------------------------------------------------------------------------------------------

    @property
    def _admitted_names(self) -> list[str]:
        names = [key.__name__ for key in self._responses]
        if self._source is not None:
            names.append(InflowTerms.__name__)
        return sorted(names)

    def catalogue(self) -> Mapping[type, Admission]:
        """The type-level UPPER BOUND over the table's keys: the front end's drop-down, before a problem exists."""
        self._require_sealed()
        entries: dict[type, Admission] = {key: realizer.status for key, realizer in self._responses.items()}
        if self._source is not None:
            entries[InflowTerms] = self._source.status
        return entries

    def _unregistered(self, unit: object) -> Refused:
        if isinstance(unit, BoundaryCondition) and unit.return_fraction < 0.0 and self.state_cone is StateCone.NONNEGATIVE_FLUX:
            return Refused(RefusalAxis.STATE_CONE, f"a return fraction {unit.return_fraction} < 0 is outside the cone psi >= 0 of {self.method}")
        return Refused(
            RefusalAxis.PLUMBING,
            f"{type(unit).__name__} has no entry in the table of {self.method}, which admits {self._admitted_names}",
            self.unregistered_issue,
        )

    def _admit_unit(self, key: type, unit: object, weight: float, site: FaceSite, space: SpaceT | None) -> Admission:
        realizer = self._responses.get(key)
        return self._unregistered(unit) if realizer is None else realizer.admit(unit, weight, site, space)

    def _admit_source(self, terms: InflowTerms, site: FaceSite, space: SourceSpaceT | None) -> Admission:
        if self._source is None:
            return Refused(RefusalAxis.PLUMBING, f"{self.method} registers no source realizer", self.unregistered_issue)
        return self._source.admit(terms, site, space)

    def _admit_deck(self, deck: DeckTransformation, site: FaceSite, space: SpaceT | None) -> Admission:
        if self._deck is None:
            return Refused(RefusalAxis.PLUMBING, f"{self.method} registers no deck realizer", self.unregistered_issue)
        return self._deck.admit(deck, pairing_of(deck, site.face, dimension=site.dimension), site, space)

    @overload
    def admit(self, declared: DeckTransformation | BoundaryCondition, site: FaceSite, space: SpaceT | None = None) -> Admission: ...
    @overload
    def admit(self, declared: InflowTerms, site: FaceSite, space: SourceSpaceT | None = None) -> Admission: ...
    def admit(self, declared: DeckTransformation | BoundaryCondition | InflowTerms, site: FaceSite, space: SpaceT | SourceSpaceT | None = None) -> Admission:
        """ONE door, per face and per value. ``space=None`` is the type- and site-level answer (the front end).

        Any condition is the meet over the admission of its units (:func:`_units`: each nonzero chain with its weight,
        ``ZeroNetwork`` whenever a shared wall's network is EMPTY, whatever the source, or the one method-only unit) and
        of its source if that is nonzero, for EVERY condition type (the source is read through the base contract). An
        ``InflowTerms`` is the source channel's own unit (what :meth:`sample_inflow` asks about). Anything else is a
        ``TypeError``: a refusal is a statement about a declaration.

        ``space`` is the space of the declaration's own verb: the method space for a deck or a condition, the source (inflow
        half-)space for a bare ``InflowTerms``. The ``@overload`` pair types the two pairings and refuses the two crossed ones
        (an ``InflowTerms`` with a method space, a deck or condition with an inflow space); the implementation's union and its
        ``cast`` calls are the single place the pairing is trusted, because the runtime ``isinstance`` split below is what the
        overloads state. The source EMBEDDED in a condition is admitted at the site level (``space=None``) and never with the
        method space: which site accepts a source is a fact of the site, and the method space is not the space the source is
        sampled on. A bare ``InflowTerms`` is admitted with its inflow space; both paths reach the same ``SourceRealizer.admit``,
        which therefore receives ``space=None`` for an embedded source and the inflow space for a bare one (pinned by
        ``test_realizer_table_source_space.py``).
        """
        return self._admit(declared, site, space)

    def _admit(self, declared: DeckTransformation | BoundaryCondition | InflowTerms, site: FaceSite, space: SpaceT | SourceSpaceT | None) -> Admission:
        """The body of :meth:`admit` behind its overloads, shared with :meth:`require` (whose overloads are the same pair)."""
        self._require_sealed()
        if isinstance(declared, DeckTransformation):
            return self._admit_deck(declared, site, cast("SpaceT | None", space))
        if isinstance(declared, InflowTerms):
            return self._admit_source(declared, site, cast("SourceSpaceT | None", space))
        if not isinstance(declared, BoundaryCondition):
            raise TypeError(f"{self.method}: admit takes a DeckTransformation, a BoundaryCondition or an InflowTerms, got {declared!r}")
        members = [self._admit_unit(key, unit, weight, site, cast("SpaceT | None", space)) for key, unit, weight in _units(declared)]
        if not declared.source.is_zero:
            members.append(self._admit_source(declared.source, site, None))
        return Admission.meet(*members)

    def _refusal(self, declared: object, site: FaceSite, admission: Refused) -> BoundaryRefused:
        message = f"{self.method} realizes {self._admitted_names}; not {type(declared).__name__} ({admission.axis.name}: {admission.reason})"
        return BoundaryRefused(message, declared, method=self.method, face=site.face, admission=admission)

    @overload
    def require(self, declared: DeckTransformation | BoundaryCondition, site: FaceSite, space: SpaceT | None = None) -> None: ...
    @overload
    def require(self, declared: InflowTerms, site: FaceSite, space: SourceSpaceT | None = None) -> None: ...
    def require(self, declared: DeckTransformation | BoundaryCondition | InflowTerms, site: FaceSite, space: SpaceT | SourceSpaceT | None = None) -> None:
        """:meth:`admit`, then the raise of :class:`BoundaryRefused` iff refused: the ONE definition of the refusal."""
        admission = self._admit(declared, site, space)
        if isinstance(admission, Refused):
            raise self._refusal(declared, site, admission)

    def realize(self, declared: DeckTransformation | BoundaryCondition, site: FaceSite, space: SpaceT) -> OpT:
        """The method's operator for ``declared`` on ``site``: ``admit`` first, then the functor on the network."""
        self.require(declared, site, space)
        if isinstance(declared, DeckTransformation):
            assert self._deck is not None  # require() passed: a deck realizer is registered  # noqa: S101
            return self._deck.realize(declared, pairing_of(declared, site.face, dimension=site.dimension), site, space)
        assert isinstance(declared, BoundaryCondition)  # require() passed: admit raised TypeError for anything else  # noqa: S101
        return reduce(operator.add, (self._responses[key].realize(unit, weight, site, space) for key, unit, weight in _units(declared)))

    def sample_inflow(self, terms: InflowTerms, site: FaceSite, space: SourceSpaceT) -> np.ndarray:
        """The ONE source verb: ``require`` the source, then the registered :meth:`SourceRealizer.sample`."""
        self.require(terms, site, space)
        assert self._source is not None  # require() passed: a source realizer is registered  # noqa: S101
        return self._source.sample(terms, site, space)


class ResponseRealizer(ABC, Generic[K, SpaceT, OpT]):
    """A method's realizer of one unit (a chain, ``ZeroNetwork`` or a method-only condition); registers by class keyword."""

    #: The catalogue entry: the upper bound of :meth:`admit` over all values and sites.
    status: ClassVar[Admission]

    def __init_subclass__(cls, *, table: RealizerTable | None = None, realizes: type | None = None, **kwargs: object) -> None:
        super().__init_subclass__(**kwargs)
        if (table is None) != (realizes is None):
            raise TypeError(f"{cls.__name__}: pass both table= and realizes= to register, or neither for an intermediate base")
        if table is not None and realizes is not None:
            table._register_response(realizes, cls())

    def admit(self, unit: K, weight: float, site: FaceSite, space: SpaceT | None = None) -> Admission:
        """Default: :attr:`status`; value-, face- and discretization-dependent refusals override."""
        return self.status

    @abstractmethod
    def realize(self, unit: K, weight: float, site: FaceSite, space: SpaceT) -> OpT:
        """The operator of ``weight * unit`` on ``site``."""


class RefusingRealizer(ResponseRealizer[K, SpaceT, OpT]):
    """A unit a table registers only to REFUSE, with a reason sharper than the unregistered door's (a named axis, a site-dependent cell).

    Its admission is a ``Refused`` at every value and site, so the table's ``realize`` (which runs ``require`` first) never reaches the
    body: the one unreachable body is here, instead of one per method. A registered subclass must state a refused :attr:`status`
    (checked at class creation); an ``admit`` override that returned an admission would reach :meth:`realize` and raise.
    """

    def __init_subclass__(cls, *, table: RealizerTable | None = None, realizes: type | None = None, **kwargs: object) -> None:
        if table is not None and not isinstance(getattr(cls, "status", None), Refused):
            raise TypeError(f"{cls.__name__}: a RefusingRealizer states a Refused status, got {getattr(cls, 'status', None)!r}")
        super().__init_subclass__(table=table, realizes=realizes, **kwargs)

    def realize(self, unit: K, weight: float, site: FaceSite, space: SpaceT) -> OpT:
        raise AssertionError(f"{type(self).__name__}: the unit is refused at admission, so it is never realized")


class SourceRealizer(ABC, Generic[SourceSpaceT]):
    """Exactly one per table that admits a source: the realizer of the key ``InflowTerms``. A source is an array, never an operator."""

    status: ClassVar[Admission]

    def __init_subclass__(cls, *, table: RealizerTable | None = None, **kwargs: object) -> None:
        super().__init_subclass__(**kwargs)
        if table is not None:
            table._register_source(cls())

    def admit(self, terms: InflowTerms, site: FaceSite, space: SourceSpaceT | None = None) -> Admission:
        """Default: :attr:`status`."""
        return self.status

    def sample(self, terms: InflowTerms, site: FaceSite, space: SourceSpaceT) -> np.ndarray:
        """Default: :func:`sample_terms`; SN overrides to run the source-placement check on the sampled ``q``."""
        return sample_terms(terms, space)


class DeckRealizer(ABC, Generic[SpaceT, OpT]):
    """Exactly one per table: the realizer of the decks, one exhaustive match over ``DeckTransformation``."""

    def __init_subclass__(cls, *, table: RealizerTable | None = None, **kwargs: object) -> None:
        super().__init_subclass__(**kwargs)
        if table is not None:
            table._register_deck(cls())

    @abstractmethod
    def admit(self, deck: DeckTransformation, pairing: FacePairing, site: FaceSite, space: SpaceT | None = None) -> Admission:
        """How well the method couples the faces of ``pairing``."""

    @abstractmethod
    def realize(self, deck: DeckTransformation, pairing: FacePairing, site: FaceSite, space: SpaceT) -> OpT:
        """The method's gluing operator for ``pairing``."""
