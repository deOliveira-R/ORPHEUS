r"""Abstract base for boundary conditions: :class:`BoundaryTraceLaw`.

A :class:`BoundaryTraceLaw` is a **pure descriptor** for an affine
boundary law in the form

.. math::

    \gamma_- \psi \;=\; R\,G\,\gamma_+ \psi \;+\; q.

It carries the law's parameters (``axis``, ``albedo``, ``source``,
…), the five universal ``assert_*`` invariants (per Grand Report v3
§16A.12), the registry plumbing (via
:class:`~orpheus.numerics.registry.RegistryMixin`), and a minimal
algebra (``+``, ``-``, ``*``, ``/``, ``-``) that returns
:class:`~orpheus.geometry.boundary._composition.LawSum` /
:class:`~orpheus.geometry.boundary._composition.LawScaled` nodes.
It is **NOT** a callable operator — there is no ``apply`` method
and no :class:`~orpheus.numerics.operator.LinearOperator`
inheritance. The realizer is the **sole** bridge from descriptor to
operator (see
:class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer.realize`).

This is the post-Issue-#186 (Scope B3 + β2) architecture: the §16A.3
three-layer split (descriptor / realizer / operator) is enforced by
the type system rather than by convention. Direct
``law.apply(psi)`` calls fail at the type-checker / linter level —
the only way to get a callable from a :class:`BoundaryTraceLaw` is to
realize it.

Wave 7 of the original 12-wave refactor merged the legacy
``BoundaryOperator`` ABC into :class:`BoundaryTraceLaw`. During the
transition the package's ``__init__.py`` carried a deprecated
``BoundaryOperator = BoundaryTraceLaw`` alias so the production
import sites kept working unchanged; that alias was retired in Wave O
step O.4a.1 once every consumer had migrated to the canonical name.

References
----------

* Grand Report v3 §16A.1–3 (affine boundary form + three-layer
  decomposition), §16A.5 (trace-correct vacuum), §16A.12 (universal
  invariants), §15A.2 (sweep-cycle detection).
* ``.claude/plans/transient-giggling-cake.md`` -- Wave 3 / Wave 4 /
  Wave 7 briefs (the 12-wave refactor).
* ``.claude/plans/bc-trace-law-descriptor-cleanup.md`` -- this scope
  (Issue #186 B3 + β2: drop ``apply`` from descriptors, ship
  :class:`LawSum` / :class:`LawScaled`).
"""

from __future__ import annotations

from abc import ABC
from typing import TYPE_CHECKING, Any, ClassVar, Optional

import numpy as np

from orpheus.numerics.registry import RegistryMixin

from ._source import InflowSourceSpec, NoSource

if TYPE_CHECKING:
    from orpheus.numerics.operator import LinearOperator
    from orpheus.numerics.quadrature import Quadrature

    from ._composition import LawNode, LawScaled, LawSum


__all__ = ["BoundaryTraceLaw", "law_permutes_ordinates"]


def law_permutes_ordinates(law: "BoundaryTraceLaw") -> bool:
    r"""Does this law's realized composite :math:`R\,G` permute the angular
    index?

    The question four production sites ask, and the reason it has one home.

    Why a FUNCTION over both tiers, and not ``geometry_map.permutes_ordinates``
    ---------------------------------------------------------------------------

    Until campaign phase **B3.4b** every specular pairing lived in :math:`G`,
    so the four callers each spelled ``law.geometry_map.permutes_ordinates``
    inline and were right. The 2026-08-01 ruling put a pairing in :math:`R`
    too — a polished wall's specular return is *constitutive*, not a symmetry
    of the domain — so ``AlbedoBoundary(α, SpecularReturn(a))`` permutes
    ordinates with ``G = IdentityMap``. Four inline spellings then became four
    half-right answers, and since that law equals ``ReflectiveBoundary(a, α)``
    as a matrix, each was a place where two identical operators behaved
    differently:

    * :func:`~orpheus.sn.loss_representation.sweep_schedule` — the reflecting-face
      set, which decides the ``B_lower`` / ``B_upper`` schedule split. A face
      that couples ordinate :math:`n` to its mirror partner needs lagging
      whichever tier the pairing sits in.
    * :mod:`~orpheus.sn.acceleration.dsa` — the low-order admission guard.
    * ``_has_ruled_corner_action`` — whether the off-quadrature
      :math:`\mu = \pm 1` ray is expressible.
    * ``RadialCharacteristicBoundaryOperator._reflect_corner`` — the swap
      itself.

    Why not a shared Protocol member
    --------------------------------

    The tidy alternative — declare ``permutes_ordinates`` on
    :class:`~orpheus.geometry.boundary.BoundaryResponseKernel` so both tiers
    answer uniformly — is exactly wrong.
    :class:`~orpheus.geometry.boundary.SpecularReemission` already carries
    ``is_adjointable``; the extra member would complete
    :class:`~orpheus.geometry.boundary.BoundaryGeometryMap`'s member set, so a
    response would satisfy the geometry Protocol **structurally**.
    ``tests/geometry/test_boundary_factors.py`` asserts those two tiers are
    disjoint, precisely to stop a response from posing as a geometry — the
    conflation phase B3.0 corrected. A convenience member is not worth
    disarming that guard, so the pairing is asked of each tier in that tier's
    own vocabulary and joined here.

    .. note::

       No shipped mesh can carry an albedo law yet
       (``SNMesh.BOUNDARY_OPERATOR_REGISTRY`` omits it, and all four callers
       read ``sn_mesh.bc[face].law``), so this is a **latent** correction: it
       fires the day issue **#189** registers the law. It is made now because
       B3.4b is what makes the four spellings wrong, and leaving a known-false
       predicate behind a registry gate is how a landmine gets planted.
    """
    from ._factors import SpecularReemission

    return bool(
        law.geometry_map.permutes_ordinates
        or isinstance(law.response_kernel, SpecularReemission)
    )


# ═══════════════════════════════════════════════════════════════════════
# BoundaryTraceLaw — pure descriptor for an affine boundary law.
# ═══════════════════════════════════════════════════════════════════════


class BoundaryTraceLaw(RegistryMixin, ABC):
    r"""Method-agnostic boundary law in the affine form
    :math:`\gamma_- \psi = R\,G\,\gamma_+ \psi + q`.

    Three properties are declared for the affine form's factors,
    but only ONE of them is real today:

    * :attr:`source` -- the prescribed inflow :math:`q`, defaulting
      to :class:`NoSource` (the homogeneous case). **Populated**
      (:class:`~orpheus.geometry.boundary.PrescribedInflow`
      overrides it) and **read** — by
      :meth:`assert_source_lives_on_incoming_trace` here and by the
      SN realizer's prescribed-inflow arm.
    * :attr:`geometry_map` -- the geometric operator :math:`G` (a
      permutation, pushforward, angular average, spatial wrap).
      **Unpopulated**: every concrete law inherits the ``None``
      below, and no production code reads it.
    * :attr:`response_kernel` -- the scalar amplitude / kernel
      :math:`R` (albedo, white-current scaling). **Unpopulated**
      likewise — the realizers reach the same number through
      ``law.albedo`` instead.

    Campaign phase **B1** mints the typed ``G`` / ``R``
    specification objects (Grand Report v3 §16A.2's
    ``BoundaryGeometryMap`` / ``BoundaryResponseKernel``) and
    populates the two empty properties across all seven laws; the
    declarations are kept for that landing rather than retired.

    The seven concrete subclasses are ``VacuumInflow``,
    ``ReflectiveBoundary``, ``WhiteBoundary``, ``AlbedoBoundary``,
    ``PeriodicBoundary``, ``PrescribedInflow`` and
    ``ZeroFluxBoundary``. This ABC ships the universal ``assert_*``
    invariants (per §16A.12), the registry plumbing, and a minimal
    algebra that composes laws into :class:`LawSum` /
    :class:`LawScaled` nodes.

    No ``apply``
    ------------
    Descriptors are **not** callable. The §16A.3 three-layer
    architecture demands that operatorhood live in a separate layer
    (realized via
    :meth:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer.realize`),
    and Issue #186 cleanup (B3 + β2) makes this a type-level
    constraint: there is no abstract ``apply`` method to satisfy,
    nor :class:`LinearOperator` inheritance to inherit one
    from. Calling ``law.apply(psi)`` raises ``AttributeError`` at
    runtime; type checkers flag it statically.

    Realization is the only bridge:

    .. code-block:: python

        from orpheus.geometry.boundary import ReflectiveBoundary
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer
        from orpheus.sn.mesh.method_space import SNMethodSpace

        law = ReflectiveBoundary(axis="x", albedo=0.5)
        ms = SNMethodSpace.minimal(quad)
        op = SNBoundaryRealizer().realize(law, ms)
        psi_in = op.apply(psi_out)   # 1-arg LinearOperator

    Rank-N composition via :class:`LawSum` / :class:`LawScaled`:

    .. code-block:: python

        composed = 0.3 * ReflectiveBoundary(axis="x") + 0.7 * WhiteBoundary(...)
        # composed is LawSum(LawScaled(0.3, ...), LawScaled(0.7, ...))
        # Still NOT callable. Realize the tree (the walker is
        # method-blind — pass the method's own realizer):
        op_tree = realize_recursively(composed, ms, SNBoundaryRealizer())
        # op_tree is OperatorSum(ScaledOperator(0.3, ...), ...).

    Registry
    --------
    :class:`BoundaryTraceLaw` IS the single registry root for all
    boundary conditions. Concrete subclasses self-register via the
    ``key=`` class-creation kwarg
    (``class VacuumInflow(BoundaryTraceLaw, key="vacuum"): ...``).
    The :attr:`registry` ``ClassVar`` is the dict mapping kind
    strings to concrete classes.
    """

    registry: ClassVar[dict[str, type["BoundaryTraceLaw"]]] = {}

    @classmethod
    def _registry_base(cls) -> type:
        return BoundaryTraceLaw

    @property
    def kind(self) -> str:
        r"""The registry key this law self-registered under.

        ONE definition for the whole family, derived from the
        :attr:`~orpheus.numerics.registry.RegistryMixin.key` ``ClassVar`` that
        ``__init_subclass__(key=...)`` already populates — so the tag can never
        drift from the registry that indexes the class, and cannot be set to an
        arbitrary string at construction.

        Before this, ``kind`` was a mutable dataclass FIELD on
        :class:`~orpheus.geometry.boundary.VacuumInflow` (so
        ``VacuumInflow(kind="banana")`` constructed a law whose tag matched no
        registry entry), a property on two other laws, and **absent entirely on
        four** — ``White``/``Albedo``/``Periodic``/``ZeroFlux`` raised
        ``AttributeError``, which is why production reads it defensively as
        ``getattr(law, "kind", None)``.

        .. note::

           Reading this tag to decide behaviour is the stringly-typed dispatch
           that campaign phase **B2** removes in favour of structural queries on
           the law's ``geometry_map`` / ``response_kernel``. This property makes
           the tag *honest* in the meantime; it is not an endorsement of
           dispatching on it.
        """
        key = type(self).key
        if key is None:  # pragma: no cover — every concrete law passes key=
            raise TypeError(
                f"{type(self).__name__} did not self-register a boundary-law "
                f"key; concrete laws are declared as "
                f"`class X(BoundaryTraceLaw, key=\"...\")`."
            )
        return key

    @property
    def geometry_map(self) -> Any:
        r"""The deck transformation :math:`G : \Gamma_+ \to \Gamma_-`.

        **Populated on every concrete law since campaign phase B1**, and read
        by production: the sweep schedule's reflecting-face set, the DSA
        admission guard, the ray-corner predicate and the diffusion realizer's
        periodic refusal all ask it structural questions instead of comparing
        ``kind`` strings (phase B2 moved them). The ``None`` here is the ABC's
        default for a law that declares nothing; no shipped law uses it.

        Membership is decided by the two tests in
        :mod:`~orpheus.geometry.boundary._factors` — multiplicativity
        (necessary) and *is it a quotient of the domain* (sufficient) — whence
        exactly one of :math:`G`, :math:`R` is non-trivial.
        """
        return None

    @property
    def response_kernel(self) -> Any:
        r"""The constitutive response :math:`R : \Gamma_- \to \Gamma_-`.

        **Populated on every concrete law since campaign phase B1.** The
        diffusion realizer reads :attr:`~orpheus.geometry.boundary.ScalarResponse.amplitude`
        as its entire first stage, the SN realizer dispatches the albedo family
        on the kernel's TYPE (B3.4b), and the vacuum-detection sites read
        ``is_zero``. Every one of them formerly reached ``law.albedo`` as a
        bare float or compared a ``kind`` tag.
        """
        return None

    @property
    def source(self) -> InflowSourceSpec:
        r"""The prescribed inflow :math:`q`. Default:
        :class:`NoSource`."""
        return NoSource()

    # ------------------------------------------------------------------
    # Descriptor-tree algebra (Issue #186 / β2). Returns
    # LawSum / LawScaled nodes; never returns an operator.
    # ------------------------------------------------------------------

    def __add__(self, other: "LawNode") -> "LawSum":
        from ._composition import LawSum
        return LawSum(self, other)

    def __radd__(self, other: "LawNode") -> "LawSum":
        from ._composition import LawSum
        return LawSum(other, self)

    def __sub__(self, other: "LawNode") -> "LawSum":
        from ._composition import LawScaled, LawSum
        return LawSum(self, LawScaled(-1.0, other))

    def __rsub__(self, other: "LawNode") -> "LawSum":
        from ._composition import LawScaled, LawSum
        return LawSum(other, LawScaled(-1.0, self))

    def __mul__(self, scalar: float) -> "LawScaled":
        from ._composition import LawScaled
        if not isinstance(scalar, (int, float)):
            return NotImplemented
        return LawScaled(float(scalar), self)

    __rmul__ = __mul__

    def __truediv__(self, scalar: float) -> "LawScaled":
        from ._composition import LawScaled
        if not isinstance(scalar, (int, float)):
            return NotImplemented
        return LawScaled(1.0 / float(scalar), self)

    def __neg__(self) -> "LawScaled":
        from ._composition import LawScaled
        return LawScaled(-1.0, self)

    # ------------------------------------------------------------------
    # Universal invariants (§16A.12, §27.6) -- no-op defaults; concrete
    # BCs override the relevant ones.
    # ------------------------------------------------------------------

    def assert_inflow_outflow_classification(
        self, quadrature: "Quadrature"
    ) -> None:
        r"""Every ordinate is either inflow or outflow (no tangential).

        Default: no-op. Laws that REQUIRE a strict partition (e.g.
        reflective) override to raise
        :class:`~orpheus.geometry.boundary._errors.IncomingOutgoingTraceClassificationError`
        on tangential ordinates.
        """

    def assert_outgoing_leakage_unconstrained(
        self, quadrature: "Quadrature"
    ) -> None:
        r"""The outgoing-trace flux is not constrained by the BC.

        Default: no-op. Constraint laws (Dirichlet outflow,
        prescribed cell-edge interface) override.
        """

    def assert_geometry_map_measure_preserving(
        self, quadrature: "Quadrature"
    ) -> None:
        r"""The geometric map preserves the angular measure
        :math:`w(\Omega)\,|\Omega\cdot n|`.

        Default: no-op. Permutation-based laws (reflective)
        override.
        """

    def assert_response_positive_if_declared(self) -> None:
        r"""If a response kernel is declared, it produces
        non-negative output.

        Default: no-op. Albedo / white override.
        """

    def assert_source_lives_on_incoming_trace(
        self,
        quadrature: "Quadrature",
        inflow_indices: Optional[np.ndarray] = None,
    ) -> None:
        r"""The source :math:`q` is nonzero only on :math:`\Gamma_-`
        (ERR-047).

        The affine form's :math:`q` term lives on the incoming trace by
        definition. An
        :class:`~orpheus.geometry.boundary._source.InflowSourceSpec`
        fills whatever block *shape* it is handed (it carries no trace
        knowledge), so the delivered :math:`q \in \Gamma_-` guarantee
        rests on the realizer knowing :math:`\Gamma_-` at all. This check
        certifies that it DOES, whenever that matters:

        * :math:`q \equiv 0` (the :class:`NoSource` default of every
          homogeneous law, or a zero-valued constant) — trivially on
          :math:`\Gamma_-`; certifiable either way.
        * :math:`q \not\equiv 0` **and** ``inflow_indices is None`` —
          the realizer cannot name :math:`\Gamma_-`, so the delivered
          source would write into slots the sweep silently discards (the
          ERR-047 missing-source bug: total inflow SHORT of the user's
          intent). Raises
          :class:`~orpheus.geometry.boundary._errors.BoundarySourceNotOnIncomingTraceError`.
        * :math:`q \not\equiv 0` **and** the inflow set is supplied — the
          realizer sizes the operator's CODOMAIN from it, so the source is
          asked for exactly :math:`|\Gamma_-|` rows and the delivered
          :math:`q` vanishes off the incoming trace because there is
          nowhere else to write. The end-to-end postcondition is pinned by
          the ERR-047 catcher test.

        .. note::

           **This is a PRESENCE check, not an entry-wise one**, and since
           campaign phase **B3.4a** it does not certify a mask. Pre-B3.4a
           the realizer emitted a full-face block and then ZEROED every
           non-inflow ordinate; the guarantee rested on that erasure, and
           this docstring said so. The narrowing to
           :math:`\Gamma_+ \to \Gamma_-` dissolved the mask along with the
           codomain that needed it — an erasure became an absence — so what
           is certified here is that the face can NAME its inflow set, not
           that anything is subsequently masked out.

        The support probe evaluates the source on the per-ordinate
        shape ``(N,)`` — the one axis the trace partition lives on;
        the :class:`InflowSourceSpec` contract ("fill exactly the
        requested shape") makes the probe representative of any
        trailing-axis block the realizer later requests.
        """
        probe = self.source.evaluate((int(quadrature.N),))
        if not np.any(probe):
            return
        if inflow_indices is None:
            from ._errors import BoundarySourceNotOnIncomingTraceError

            raise BoundarySourceNotOnIncomingTraceError(
                f"{type(self).__name__} carries a nonzero inflow source "
                f"({type(self.source).__name__}) but no inflow-trace "
                f"indices were provided, so the realized q cannot be "
                f"certified to live on Γ_- — an unmasked source writes "
                f"into outflow slots the sweep discards (ERR-047). "
                f"Realize against a method space carrying "
                f"inflow_indices (SNMethodSpace.for_face with a trace, "
                f"or explicit inflow_indices=).",
                # The registry key, matching all six sibling raise sites
                # (``law="albedo"`` / ``"reflective"`` / ``"white"``). This was
                # ``type(self).__name__`` — the only site spelling the tag as a
                # CLASS NAME ("PrescribedInflow" vs "prescribed_inflow"), so an
                # error-tag consumer saw one law under a different vocabulary.
                # ``self.kind`` is now the one derivation (B0.1), so the drift
                # cannot recur.
                law=self.kind,
            )

    # ------------------------------------------------------------------
    # Realize-time certification (#52) -- the ONE aggregate the method
    # realizers call so the §16A.12 invariants are PRODUCTION guards,
    # not test-only helpers.
    # ------------------------------------------------------------------

    def assert_realizable(
        self,
        quadrature: "Quadrature",
        *,
        inflow_indices: Optional[np.ndarray] = None,
    ) -> None:
        r"""Fire every law-level invariant against the realization data.

        The §16A.12 universal invariants and the error catalog's
        lessons (ERR-041/042/043/045/046/047) all end the same way:
        *the contract must be CHECKED at construction/realization, not
        only verified by downstream balance* — by then the miscount has
        propagated. This template method is the production seam that
        honors those lessons: a method realizer calls it ONCE at entry
        so the law arrives at its primitive construction already
        certified against the actual quadrature and trace data it is
        about to be wired to.

        .. warning::

           **The seam has exactly one production caller today** —
           :meth:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer.realize`.
           The diffusion realizer never calls it, so a law realized
           through the diffusion arm is NOT certified: e.g.
           ``WhiteBoundary(albedo=5.0)`` reaches
           :math:`\mathcal{A} = 5.0` there with no sub-Markov error.
           Read "every law arrives already certified" as scoped to
           the SN arm until the diffusion realizer adopts the seam.

        The base body fires the five universal invariants — four of
        which are empty and two of those overridden by nobody, so in
        practice the aggregate's teeth are
        :meth:`assert_source_lives_on_incoming_trace` (all laws),
        :meth:`assert_geometry_map_measure_preserving`
        (:class:`ReflectiveBoundary`) and
        :meth:`assert_response_positive_if_declared`
        (:class:`WhiteBoundary` / :class:`AlbedoBoundary`). Concrete
        laws EXTEND via ``super().assert_realizable(...)`` plus their
        law-specific invariants — :class:`ReflectiveBoundary` adds the
        involution and inflow→outflow table checks,
        :class:`WhiteBoundary` / :class:`AlbedoBoundary` add the
        sub-Markov bound — so the realizer needs no per-law knowledge
        of WHICH invariants exist.

        Parameters
        ----------
        quadrature
            The angular quadrature the law is being realized against.
        inflow_indices
            The face's inflow-ordinate indices when the method space
            carries them (``None`` on quadrature-only spaces). Consumed
            by the ERR-047 source-trace check; a nonzero source without
            an inflow set is uncertifiable and raises there.
        """
        self.assert_inflow_outflow_classification(quadrature)
        self.assert_outgoing_leakage_unconstrained(quadrature)
        self.assert_geometry_map_measure_preserving(quadrature)
        self.assert_response_positive_if_declared()
        self.assert_source_lives_on_incoming_trace(
            quadrature, inflow_indices
        )

    # ------------------------------------------------------------------
    # Method realisation hook -- guidance raise; the realizers are the
    # bridge (each owned by its method-mesh since #290 P7b).
    # ------------------------------------------------------------------

    def realize(self, method_space: Any) -> "LinearOperator":
        r"""Realise this law against a method-specific space.

        The default body raises :class:`NotImplementedError`. The
        canonical realisation path is via a method-specific realizer:

        .. code-block:: python

            from orpheus.sn.boundary.realizer import SNBoundaryRealizer
            from orpheus.sn.mesh.method_space import SNMethodSpace
            op = SNBoundaryRealizer().realize(law, SNMethodSpace.minimal(quad))
            psi_in = op.apply(psi_out)   # 1-arg LinearOperator

        Concrete laws MAY override to delegate to a specific
        realizer when the method space's ``method_name`` is known,
        but no current production caller routes through this hook —
        each method-mesh's ``realize_boundary_law`` arm invokes its
        own realizer directly (or :func:`realize_recursively` with an
        explicit realizer for rank-N composition).
        """
        raise NotImplementedError(
            f"{type(self).__name__}.realize: route through a "
            "method-specific realizer (e.g. "
            "SNBoundaryRealizer().realize(law, method_space)) or "
            "use realize_recursively for descriptor trees."
        )
