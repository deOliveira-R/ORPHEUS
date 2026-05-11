r"""Abstract base for boundary conditions: :class:`BoundaryTraceLaw`.

Wave 7 of the boundary-operator refactor merged the legacy
``BoundaryOperator`` ABC (originally defined in this file alongside
:class:`BoundaryTraceLaw`) into :class:`BoundaryTraceLaw`. The legacy
ABC's separate registry, separate ``apply(psi_out, quadrature)``
contract, and separate ``_registry_base`` override are gone; the new
ABC's ``apply(psi_out, *args, **kwargs)`` accommodates both the legacy
2-arg call from production sweeps and the Wave 10 1-arg call once
consumers route through the realizer.

The :pyattr:`__init__.py` of the package re-exports ``BoundaryOperator
= BoundaryTraceLaw`` as a deprecated alias so the 8 production import
sites and the existing tests that ``from orpheus.geometry.boundary
import BoundaryOperator`` keep working unchanged. The alias is
scheduled for removal in a future cleanup wave (see Wave 11).

The remaining content of this module is :class:`BoundaryTraceLaw`
itself — the method-agnostic affine boundary law in
:math:`\gamma_- \psi = R\,G\,\gamma_+ \psi + q` form, with its three
first-class properties (``geometry_map``, ``response_kernel``,
``source``), its five ``assert_*`` universal invariants (per Grand
Report v3 §16A.12), its :pydata:`creates_sweep_cycle` ClassVar, its
``realize`` hook into the Wave-5 :class:`BoundaryRealizerRegistry`,
and its abstract ``apply`` contract.

References
----------

* Grand Report v3 §16A.1-3 (affine boundary form + three-layer
  decomposition), §16A.12 (universal invariants), §15A.2
  (sweep-cycle detection).
* ``.claude/plans/transient-giggling-cake.md`` -- Wave 3 / Wave 4 /
  Wave 7 briefs.
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


# ═══════════════════════════════════════════════════════════════════════
# BoundaryTraceLaw — the ONE ABC (post Wave 7 ABC-merge).
# ═══════════════════════════════════════════════════════════════════════


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
    The :meth:`apply` signature carries ``(psi_out, *args, **kwargs)``
    so concretes can accept the legacy 2-arg form
    ``(psi_out, quadrature)`` (consumed by production sweeps in
    Waves 7-9) or the 1-arg form ``(psi)`` post-Wave 10 (once the
    realiser captures the quadrature at construction).

    Registry
    --------
    :class:`BoundaryTraceLaw` IS the single registry root for all
    boundary conditions. The pre-Wave-7 separate ``BoundaryOperator``
    registry has been merged — there is now ONE
    :pyattr:`BoundaryTraceLaw.registry` dict containing every
    concrete BC (``vacuum`` / ``reflective`` / ``white`` /
    ``periodic`` / ``albedo`` / ``mixed`` / ``prescribed_inflow``).
    The legacy ``BoundaryOperator`` symbol is now an alias of this
    class via the package ``__init__.py`` (deprecated; remove in a
    future cleanup wave).
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
    def apply(self, psi_out: np.ndarray) -> np.ndarray:
        r"""Apply the boundary law to the outgoing flux trace.

        Concrete subclasses MAY accept additional positional /
        keyword arguments in their signature — typically
        ``quadrature: AngularQuadrature | None = None`` for laws
        whose geometric operator is quadrature-dependent (see
        :class:`~orpheus.geometry.boundary.reflective.ReflectiveBoundary`
        and
        :class:`~orpheus.geometry.boundary.white.WhiteBoundary`).
        The abstract signature is intentionally strict 1-arg to
        signal the modernized contract: the canonical application
        path is through
        :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer.realize`,
        which captures any required method-specific metadata
        (quadrature, mesh face, inflow trace) at realization time
        and returns a 1-arg
        :class:`~orpheus.numerics.operator.LinearOperator`.

        Liskov substitutability: a concrete
        ``apply(self, psi_out, quadrature=None)`` is compatible
        with the abstract ``apply(self, psi_out)`` because the
        additional parameter has a default and is positionally /
        keyword optional. Direct callers that want the realized
        :math:`\Omega \cdot \hat n < 0`-correct semantics should
        go through :meth:`SNBoundaryRealizer.realize`; the bare
        :meth:`apply` body of each concrete BC documents whether
        the direct-call path is exact (e.g.
        :class:`AlbedoBoundary`) or a backward-compat fallback
        (e.g. :class:`VacuumInflow` returning legacy zeros-all in
        the absence of per-face inflow indices).

        Issue #176 / C176.4: previous Wave-3 signature carried
        ``*args, **kwargs`` to accommodate the legacy 2-arg form
        during the Waves 7-10 migration; that flexibility is no
        longer needed now that the realizer-or-direct contract is
        fully shipped.
        """
        ...
