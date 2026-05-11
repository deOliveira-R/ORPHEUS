r"""The :class:`BoundaryTraceLaw` ABC -- method-agnostic boundary law.

Per Grand Report v3 §16A.3 lines 2841-2860, the boundary condition
splits into three layers:

1. **Trace structure** -- :class:`~orpheus.numerics.trace_space.InflowTraceSpace`
   / :class:`~orpheus.numerics.trace_space.OutflowTraceSpace` (Wave 2).
2. **Physical law** -- this module's :class:`BoundaryTraceLaw`: the
   method-agnostic affine map

   .. math::

       \gamma_- \psi \;=\; R\,G\,\gamma_+ \psi \;+\; q

   with three terms (geometry map :math:`G`, response kernel
   :math:`R`, source :math:`q`).
3. **Method realisation** -- ``BoundaryRealizer`` (Wave 5): one
   concrete realiser per transport method (SN, MoC, MC, CP,
   diffusion) that turns a :class:`BoundaryTraceLaw` into a
   :class:`~orpheus.numerics.operator.LinearOperator` consumable by
   that method's sweep / solver.

This wave ships the ABC and the universal ``assert_*`` invariants
(§16A.12 + §27.6). Concrete BCs ship in Wave 7 (rename / internal
delegation of the legacy operators).

The registry follows the same pattern as
:class:`~orpheus.geometry.boundary.BoundaryOperator`: each concrete
subclass declares ``key=`` in its class statement and self-registers
via :class:`~orpheus.numerics.registry.RegistryMixin`. The two
registries are kept independent by overriding :meth:`_registry_base`
to return :class:`BoundaryTraceLaw` (vs. ``BoundaryOperator``) --
verified by ``tests/geometry/test_boundary_trace_law.py``.

References
----------

* Grand Report v3 §16A.1-3 (affine boundary form + three-layer
  decomposition), §16A.12 (universal invariants), §15A.2
  (sweep-cycle detection).
* ``.claude/plans/transient-giggling-cake.md`` -- Wave 3 brief.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any, ClassVar, Optional

import numpy as np

from orpheus.numerics.operator import CAP_APPLY, LinearOperatorMixin
from orpheus.numerics.registry import RegistryMixin

from ._source import BoundarySource, NoSource

if TYPE_CHECKING:
    from orpheus.numerics.operator import LinearOperator
    from orpheus.numerics.space import FunctionSpace
    from orpheus.sn.quadrature import AngularQuadrature


__all__ = ["BoundaryTraceLaw"]


class BoundaryTraceLaw(LinearOperatorMixin, RegistryMixin, ABC):
    r"""Method-agnostic boundary law in the affine form
    :math:`\gamma_- \psi = R\,G\,\gamma_+ \psi + q`.

    Three first-class properties:

    * :attr:`geometry_map` -- the geometric operator :math:`G` (a
      permutation, pushforward, angular average, spatial wrap).
    * :attr:`response_kernel` -- the scalar amplitude / kernel
      :math:`R` (albedo, white-current scaling).
    * :attr:`source` -- the prescribed inflow :math:`q`, defaulting
      to :class:`NoSource` (the homogeneous case).

    Concrete subclasses (Wave 7) populate these. This ABC ships the
    universal ``assert_*`` invariants (per §16A.12) and a
    :attr:`creates_sweep_cycle` ClassVar signal used by the SN
    sweep planner to detect reflective / periodic faces that create
    cycles in the dependency graph.

    Notes
    -----
    The :meth:`apply` signature in this wave still carries the
    legacy multi-arg ``(psi_out, ...)`` form so concretes can
    delegate to the existing
    :class:`~orpheus.geometry.boundary.BoundaryOperator` body
    during Wave 7. Wave 10 drops the trailing ``quadrature``
    argument once all consumers go through the realiser path.

    Registry
    --------
    :class:`BoundaryTraceLaw` and :class:`BoundaryOperator`
    maintain SEPARATE registries -- each registry root re-declares
    its own ``registry: ClassVar[dict]`` and overrides
    :meth:`_registry_base` to return itself. This isolation is
    load-bearing for Wave 7's rename of concrete BCs onto the
    :class:`BoundaryTraceLaw` hierarchy: new concretes register
    under the new root without colliding with the legacy
    ``BoundaryOperator.registry`` entries.
    """

    registry: ClassVar[dict[str, type["BoundaryTraceLaw"]]] = {}
    capabilities: ClassVar[frozenset[str]] = frozenset({CAP_APPLY})

    creates_sweep_cycle: ClassVar[bool] = False
    r"""Whether the realised operator creates a cycle in the SN
    sweep DAG. ``True`` on ``ReflectiveBoundary`` and
    ``PeriodicBoundary`` (Wave 7); ``False`` on all other laws.
    Used by §15A.2 cycle detection (Grand Report v3 lines
    2155-2160)."""

    @classmethod
    def _registry_base(cls) -> type:
        return BoundaryTraceLaw

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        return None

    @property
    def range(self) -> Optional["FunctionSpace"]:
        return None

    @property
    def geometry_map(self) -> Any:
        r"""The geometric operator :math:`G` (permutation,
        pushforward, ...).

        Default: ``None``. Concrete BCs (Wave 7) populate.
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
    def source(self) -> BoundarySource:
        r"""The prescribed inflow :math:`q`. Default:
        :class:`NoSource`."""
        return NoSource()

    # ------------------------------------------------------------------
    # Universal invariants (§16A.12, §27.6) -- no-op defaults; concrete
    # BCs override the relevant ones in Wave 7.
    # ------------------------------------------------------------------

    def assert_inflow_outflow_classification(
        self, quadrature: "AngularQuadrature"
    ) -> None:
        r"""Every ordinate is either inflow or outflow (no tangential).

        Default: no-op. Laws that REQUIRE a strict partition (e.g.
        reflective) override to raise
        :class:`~orpheus.geometry.boundary._errors.IncomingOutgoingTraceClassificationError`
        on tangential ordinates.
        """

    def assert_outgoing_leakage_unconstrained(
        self, quadrature: "AngularQuadrature"
    ) -> None:
        r"""The outgoing-trace flux is not constrained by the BC.

        Default: no-op. Constraint laws (Dirichlet outflow,
        prescribed cell-edge interface) override.
        """

    def assert_geometry_map_measure_preserving(
        self, quadrature: "AngularQuadrature"
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
        self, quadrature: "AngularQuadrature"
    ) -> None:
        r"""The source :math:`q` is nonzero only on
        :math:`\Gamma_-`.

        Default: no-op (concrete laws inherit unless they need a
        stricter check). The :class:`NoSource` default trivially
        satisfies the invariant because :math:`q \equiv 0`.
        """

    # ------------------------------------------------------------------
    # Method realisation hook (Wave 5 wires this through the registry).
    # ------------------------------------------------------------------

    def realize(self, method_space: Any) -> "LinearOperator":
        r"""Realise this law against a method-specific space.

        Defers to ``BoundaryRealizerRegistry`` in Wave 5; for Wave
        3 this raises :class:`NotImplementedError` because the
        realiser layer hasn't shipped yet.
        """
        raise NotImplementedError(
            "BoundaryTraceLaw.realize: realiser registry ships in "
            "Wave 5 (feature/boundary-realizer-protocol). See plan "
            "transient-giggling-cake.md."
        )

    # ------------------------------------------------------------------
    # apply contract -- abstract in this wave (concrete BCs in Wave 7
    # delegate to the realiser path).
    # ------------------------------------------------------------------

    @abstractmethod
    def apply(
        self, psi_out: np.ndarray, *args: Any, **kwargs: Any,
    ) -> np.ndarray:
        r"""Boundary-law application.

        Concrete signature defined by subclasses. During the Wave
        7 transition, concrete BCs match the legacy 2-arg form
        ``(psi_out, quadrature)`` so they can delegate to the
        existing
        :class:`~orpheus.geometry.boundary.BoundaryOperator` body.
        Wave 10 collapses the signature to 1-arg ``(psi)`` once
        the realiser captures the quadrature at construction.
        """
        ...
