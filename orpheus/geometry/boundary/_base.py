r"""Abstract bases for boundary conditions: legacy + trace-law.

This module hosts **two** abstract bases that coexist during the
Wave 3-Wave 11 transition:

1. :class:`BoundaryOperator` -- the legacy ABC, 2-arg
   ``apply(psi_out, quadrature)`` signature. The five concrete
   subclasses (:class:`~orpheus.geometry.boundary.vacuum.VacuumBoundaryOperator`,
   :class:`~orpheus.geometry.boundary.reflective.SpecularBoundaryOperator`,
   :class:`~orpheus.geometry.boundary.white.WhiteBoundaryOperator`,
   :class:`~orpheus.geometry.boundary.periodic.PeriodicBoundaryOperator`,
   :class:`~orpheus.geometry.boundary.albedo.AlbedoBoundaryOperator`)
   plus the rank-N composer
   (:class:`~orpheus.geometry.boundary.mixed.MixedBoundaryOperator`)
   all inherit from this base today.

2. :class:`BoundaryTraceLaw` -- the new method-agnostic ABC for the
   §16A.3 three-layer decomposition. Concrete BCs will be rebased
   onto this ABC in Wave 7 (rename + internal delegation).

The two ABCs maintain SEPARATE registries: each declares its own
``registry: ClassVar[dict]`` and overrides
:meth:`_registry_base` to return itself. This isolation is
load-bearing for Wave 7: new concretes can register under
:class:`BoundaryTraceLaw` without colliding with the legacy
:attr:`BoundaryOperator.registry` entries.

Wave 3 added :class:`BoundaryTraceLaw`. Wave 4 (this file's
expansion) moved :class:`BoundaryOperator` here from the package
``__init__.py`` as part of the source-layout split. The per-BC
concretes moved into per-BC submodules in the same wave.

References
----------

* Grand Report v3 §16A.1-3 (affine boundary form + three-layer
  decomposition), §16A.12 (universal invariants), §15A.2
  (sweep-cycle detection).
* ``.claude/plans/transient-giggling-cake.md`` -- Wave 3 / Wave 4
  briefs.
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


__all__ = ["BoundaryOperator", "BoundaryTraceLaw"]


# ═══════════════════════════════════════════════════════════════════════
# Legacy ABC (will be deprecated in Wave 7 when concretes rebase onto
# BoundaryTraceLaw).
# ═══════════════════════════════════════════════════════════════════════


class BoundaryOperator(LinearOperatorMixin, RegistryMixin, ABC):
    r"""Abstract base for a tensor-decomposed boundary condition.

    A :class:`BoundaryOperator` is the runtime representation of the operator

    .. math::

        R = \sum_\alpha G_\alpha \otimes A_\alpha

    that consumes the outgoing angular flux at a face and produces the
    incoming angular flux. Concrete implementations may be **rank-1**
    (a single :math:`G \otimes A` term: vacuum, specular, white,
    periodic, albedo) or **rank-N**
    (:class:`~orpheus.geometry.boundary.mixed.MixedBoundaryOperator` --
    a sum of the rank-1 primitives).

    The :class:`~orpheus.sn.quadrature.AngularQuadrature` argument
    lets each concrete BC reach into the quadrature's structural
    metadata (``reflection_index(axis)``, weights, level structure)
    without the contract exposing those as required attributes -- the
    responsibility for being able to query the quadrature lives in
    the BC, not in every consumer.

    Issue 9.6 lift
    ==============

    The earlier ``BoundaryOperator(Protocol)`` declaration was
    duck-typed; Issue 9.6 lifted it into a concrete ABC that:

    * inherits :class:`~orpheus.numerics.operator.LinearOperatorMixin`
      so concrete BCs participate in the operator algebra
      (``+`` / ``*`` / ``@``);
    * inherits :class:`~orpheus.numerics.registry.RegistryMixin` so
      concrete BCs self-register under a string ``key=`` class-creation
      kwarg, accessible via :meth:`create`;
    * canonicalises the application method as :meth:`apply` (the
      Issue 9.5 rename retained the old ``apply_to_incoming`` name;
      Issue 9.6 dropped it in favour of the operator-algebra
      vocabulary).

    Each concrete subtype declares its :pydata:`capabilities`
    frozenset. The minimum contract is ``frozenset({"apply"})``;
    :class:`~orpheus.geometry.boundary.reflective.SpecularBoundaryOperator`
    ships ``apply_transpose`` in addition (load-bearing for
    sensitivity-analysis adjoint pipelines).
    """

    registry: ClassVar[dict[str, type["BoundaryOperator"]]] = {}
    capabilities: ClassVar[frozenset[str]] = frozenset({CAP_APPLY})

    @classmethod
    def _registry_base(cls) -> type:
        return BoundaryOperator

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        """Boundary-trace space the BC consumes (outgoing flux).

        Returns ``None`` on the abstract base -- concrete subtypes
        whose trace dimensions are determined at construction time
        may override to expose them. The Issue 9.6 ship-state leaves
        this as ``None`` for backward compatibility; the function
        spaces become non-trivial when an SN solver explicitly
        constructs them with the live ordinate / group counts.
        """
        return None

    @property
    def range(self) -> Optional["FunctionSpace"]:
        """Boundary-trace space the BC produces (incoming flux).

        See :attr:`domain` for the ``None`` semantics.
        """
        return None

    @abstractmethod
    def apply(
        self,
        psi_out: np.ndarray,
        quadrature: "AngularQuadrature",
    ) -> np.ndarray:
        r"""Compute the incoming angular flux from the outgoing.

        Parameters
        ----------
        psi_out : np.ndarray
            Angular flux at the boundary face, indexed over all
            ordinates of ``quadrature``. Shape ``(N_ord, ...)`` where
            the trailing axes are typically energy groups.
        quadrature : AngularQuadrature
            The angular quadrature; lets the BC query reflection
            partners, weights, and level structure.

        Returns
        -------
        np.ndarray
            Incoming angular flux at the boundary face, same shape
            as ``psi_out``. For ordinates whose direction cosine is
            *outgoing* at this face (and thus do not have an
            incoming value), the returned entries are zero -- sweeps
            consume only the incoming entries.
        """
        ...


# ═══════════════════════════════════════════════════════════════════════
# Wave-3 ABC -- BoundaryTraceLaw (method-agnostic affine boundary law).
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
