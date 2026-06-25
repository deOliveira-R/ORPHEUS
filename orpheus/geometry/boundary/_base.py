r"""Abstract base for boundary conditions: :class:`BoundaryTraceLaw`.

A :class:`BoundaryTraceLaw` is a **pure descriptor** for an affine
boundary law in the form

.. math::

    \gamma_- \psi \;=\; R\,G\,\gamma_+ \psi \;+\; q.

It carries the law's parameters (``axis``, ``albedo``, ``source``,
…), the five universal ``assert_*`` invariants (per Grand Report v3
§16A.12), the :attr:`creates_sweep_cycle` ``ClassVar`` (per
§15A.2), the registry plumbing (via
:class:`~orpheus.numerics.registry.RegistryMixin`), and a minimal
algebra (``+``, ``-``, ``*``, ``/``, ``-``) that returns
:class:`~orpheus.geometry.boundary._composition.LawSum` /
:class:`~orpheus.geometry.boundary._composition.LawScaled` nodes.
It is **NOT** a callable operator — there is no ``apply`` method
and no :class:`~orpheus.numerics.operator.LinearOperator`
inheritance. The realizer is the **sole** bridge from descriptor to
operator (see
:class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer.realize`).

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
from typing import TYPE_CHECKING, Any, ClassVar

from orpheus.numerics.registry import RegistryMixin

from ._source import InflowSourceSpec, NoSource

if TYPE_CHECKING:
    from orpheus.numerics.operator import LinearOperator
    from orpheus.numerics.quadrature import Quadrature

    from ._composition import LawNode, LawScaled, LawSum


__all__ = ["BoundaryTraceLaw"]


# ═══════════════════════════════════════════════════════════════════════
# BoundaryTraceLaw — pure descriptor for an affine boundary law.
# ═══════════════════════════════════════════════════════════════════════


class BoundaryTraceLaw(RegistryMixin, ABC):
    r"""Method-agnostic boundary law in the affine form
    :math:`\gamma_- \psi = R\,G\,\gamma_+ \psi + q`.

    Three first-class properties:

    * :attr:`geometry_map` -- the geometric operator :math:`G` (a
      permutation, pushforward, angular average, spatial wrap).
    * :attr:`response_kernel` -- the scalar amplitude / kernel
      :math:`R` (albedo, white-current scaling).
    * :attr:`source` -- the prescribed inflow :math:`q`, defaulting
      to :class:`NoSource` (the homogeneous case).

    Concrete subclasses (``VacuumInflow``, ``ReflectiveBoundary``,
    ``WhiteBoundary``, ``AlbedoBoundary``, ``PeriodicBoundary``,
    ``PrescribedInflow``) populate these. This ABC ships the
    universal ``assert_*`` invariants (per §16A.12), the
    :attr:`creates_sweep_cycle` ``ClassVar`` signal used by the SN
    sweep planner, the registry plumbing, and a minimal algebra that
    composes laws into :class:`LawSum` / :class:`LawScaled` nodes.

    No ``apply``
    ------------
    Descriptors are **not** callable. The §16A.3 three-layer
    architecture demands that operatorhood live in a separate layer
    (realized via
    :meth:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer.realize`),
    and Issue #186 cleanup (B3 + β2) makes this a type-level
    constraint: there is no abstract ``apply`` method to satisfy,
    nor :class:`LinearOperator` inheritance to inherit one
    from. Calling ``law.apply(psi)`` raises ``AttributeError`` at
    runtime; type checkers flag it statically.

    Realization is the only bridge:

    .. code-block:: python

        from orpheus.geometry.boundary import ReflectiveBoundary
        from orpheus.sn.boundary_realizer import (
            SNBoundaryRealizer, SNMethodSpace,
        )

        law = ReflectiveBoundary(axis="x", albedo=0.5)
        ms = SNMethodSpace.minimal(quad)
        op = SNBoundaryRealizer().realize(law, ms)
        psi_in = op.apply(psi_out)   # 1-arg LinearOperator

    Rank-N composition via :class:`LawSum` / :class:`LawScaled`:

    .. code-block:: python

        composed = 0.3 * ReflectiveBoundary(axis="x") + 0.7 * WhiteBoundary(...)
        # composed is LawSum(LawScaled(0.3, ...), LawScaled(0.7, ...))
        # Still NOT callable. Realize the tree:
        op_tree = realize_recursively(composed, ms)
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

    creates_sweep_cycle: ClassVar[bool] = False
    r"""Whether the realised operator creates a cycle in the SN
    sweep DAG. ``True`` on ``ReflectiveBoundary`` and
    ``PeriodicBoundary``; ``False`` on all other laws. Used by
    §15A.2 cycle detection (Grand Report v3 lines 2155–2160)."""

    @classmethod
    def _registry_base(cls) -> type:
        return BoundaryTraceLaw

    @property
    def geometry_map(self) -> Any:
        r"""The geometric operator :math:`G` (permutation,
        pushforward, ...).

        Default: ``None``. Concrete BCs populate.
        """
        return None

    @property
    def response_kernel(self) -> Any:
        r"""The scalar amplitude / kernel :math:`R` (albedo, white
        scaling).

        Default: ``None``. Concrete BCs populate.
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
        self, quadrature: "Quadrature"
    ) -> None:
        r"""The source :math:`q` is nonzero only on
        :math:`\Gamma_-`.

        Default: no-op (concrete laws inherit unless they need a
        stricter check). The :class:`NoSource` default trivially
        satisfies the invariant because :math:`q \equiv 0`.
        """

    # ------------------------------------------------------------------
    # Method realisation hook -- delegates to BoundaryRealizerRegistry.
    # ------------------------------------------------------------------

    def realize(self, method_space: Any) -> "LinearOperator":
        r"""Realise this law against a method-specific space.

        The default body raises :class:`NotImplementedError`. The
        canonical realisation path is via a method-specific realizer:

        .. code-block:: python

            from orpheus.sn.boundary_realizer import (
                SNBoundaryRealizer, SNMethodSpace,
            )
            op = SNBoundaryRealizer().realize(law, SNMethodSpace.minimal(quad))
            psi_in = op.apply(psi_out)   # 1-arg LinearOperator

        Concrete laws MAY override to delegate to a specific
        realizer when the method space's ``method_name`` is known,
        but no current production caller routes through this hook —
        callers invoke :meth:`SNBoundaryRealizer.realize` directly
        (or :func:`realize_recursively` for rank-N composition).
        """
        raise NotImplementedError(
            f"{type(self).__name__}.realize: route through a "
            "method-specific realizer (e.g. "
            "SNBoundaryRealizer().realize(law, method_space)) or "
            "use realize_recursively for descriptor trees."
        )
