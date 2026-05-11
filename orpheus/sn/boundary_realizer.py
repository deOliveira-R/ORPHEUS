r"""SN-specific BoundaryRealizer.

Dispatches by ``isinstance(law, ...)`` to the Wave-0 / Wave-1 primitives
that realise each legacy :class:`~orpheus.geometry.boundary.BoundaryOperator`
subclass as a single-arg :class:`~orpheus.numerics.operator.LinearOperator`
consumable by the SN sweep + Krylov path.

This realizer accepts the LEGACY :class:`BoundaryOperator` subclasses
(Wave 5 ship-state). Wave 7 will rebase concretes onto
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` and update the
isinstance checks accordingly.

The realisation map (legacy class → Wave-0 / Wave-1 primitive)
================================================================

* :class:`~orpheus.geometry.boundary.vacuum.VacuumBoundaryOperator` →
  :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`
  with the per-face inflow indices from the method space's
  :class:`~orpheus.numerics.trace_space.InflowTraceSpace`. **Semantic
  correction** (Wave 5 risk register entry, plan §16A.5): the legacy
  body zeroes ALL ordinates, but the §16A.10 trace-correct
  implementation zeroes ONLY the inflow ordinates so the outflow
  trace survives. The new behaviour is the right one; downstream
  consumers (sensitivity adjoints) need the outflow trace
  preserved.
* :class:`~orpheus.geometry.boundary.reflective.SpecularBoundaryOperator(axis, albedo)` →
  ``albedo * PermutationOperator(quadrature.reflection_index(axis))``
  (with the ``albedo=1.0`` fast path returning the bare
  :class:`PermutationOperator` to preserve bit-identity vs legacy).
* :class:`~orpheus.geometry.boundary.white.WhiteBoundaryOperator(axis, outward_sign, albedo)` →
  ``albedo * AngularAverageOperator.from_quadrature(quadrature, axis, outward_sign)``
  (with the ``albedo=1.0`` fast path).
* :class:`~orpheus.geometry.boundary.albedo.AlbedoBoundaryOperator(0.0)` →
  :class:`~orpheus.numerics.operator.ZeroOperator`.
* :class:`~orpheus.geometry.boundary.albedo.AlbedoBoundaryOperator(1.0)` →
  :class:`~orpheus.numerics.operator.IdentityOperator`.
* :class:`~orpheus.geometry.boundary.albedo.AlbedoBoundaryOperator(α)` for α∉{0,1} →
  ``α * IdentityOperator()``.
* :class:`~orpheus.geometry.boundary.periodic.PeriodicBoundaryOperator` →
  :class:`~orpheus.numerics.operator.PeriodicWrapOperator`.
* :class:`~orpheus.geometry.boundary.mixed.MixedBoundaryOperator(components)` →
  recursive realization: ``sum(coeff * realize(prim, method_space) for coeff, prim in components)``.

References
----------

* ``.claude/plans/transient-giggling-cake.md`` Wave 5 — C5.3 (functional
  realizer brief), §16A.5 (vacuum-trace semantic correction).
* Grand Report v3 §16A.3 lines 2841–2860 (realizer-as-third-layer
  motivation), §16A.10 lines 3160-3175 (vacuum-trace correct
  representation).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

import numpy as np

from orpheus.geometry.boundary import (
    AlbedoBoundaryOperator,
    BoundaryError,
    BoundaryRealizerRegistry,
    MixedBoundaryOperator,
    PeriodicBoundaryOperator,
    SpecularBoundaryOperator,
    VacuumBoundaryOperator,
    WhiteBoundaryOperator,
)
from orpheus.numerics.operator import (
    IdentityOperator,
    IncomingOrdinateMaskTensor,
    LinearOperator,
    PeriodicWrapOperator,
    PermutationOperator,
    ScaledOperator,
    ZeroOperator,
)
from orpheus.sn.angular_operator import AngularAverageOperator

if TYPE_CHECKING:
    from orpheus.geometry.boundary import BoundaryOperator
    from orpheus.sn.quadrature import AngularQuadrature


__all__ = ["SNBoundaryRealizer", "SNMethodSpace"]


@dataclass(frozen=True)
class SNMethodSpace:
    """Minimal SN method space — carries ``quadrature`` + optional ``face`` + optional ``inflow_indices``.

    Wave 8 will extend this with the full mesh + trace metadata
    (``orpheus/sn/method_space.py``) and let :class:`SNMesh`
    populate :attr:`inflow_indices` internally from its held
    :class:`~orpheus.numerics.trace_space.InflowTraceSpace`. For
    Wave 5, only the realizer's immediate needs are populated:

    * :attr:`quadrature` — used by every realizer dispatch path
      (reflection index, mu_n arrays, weights).
    * :attr:`face` — the face label (``"left"``, ``"right"``,
      ``"xmin"``, …) used for the vacuum-BC inflow-indices lookup.
      Informational at Wave 5 (the realizer doesn't consult it
      directly; ``inflow_indices`` is the load-bearing field).
    * :attr:`inflow_indices` — per-face inflow indices (1-D int
      array). Computed externally (typically from
      :meth:`~orpheus.numerics.trace_space.InflowTraceSpace.inflow_indices_for_face`)
      and passed in. Required for
      :class:`VacuumBoundaryOperator` realisation; ``None`` for
      laws that don't need it.
    """

    quadrature: "AngularQuadrature"
    face: Optional[str] = None
    inflow_indices: Optional[np.ndarray] = None

    @classmethod
    def minimal(cls, quadrature: "AngularQuadrature") -> "SNMethodSpace":
        """Quadrature-only method space (for non-vacuum laws / tests)."""
        return cls(quadrature=quadrature)


@BoundaryRealizerRegistry.register("SN")
class SNBoundaryRealizer:
    r"""Functional realizer for SN BCs.

    Dispatches by ``isinstance(law, ...)`` to map each legacy
    :class:`~orpheus.geometry.boundary.BoundaryOperator` subclass
    onto a Wave-0 / Wave-1 primitive composed with scalar
    amplitudes (per :ref:`bc-tensor-decomposition` and the
    realisation map in this module's docstring).

    The realizer instance carries no state; ``method_name`` is the
    only attribute. The :meth:`realize` method is stateless — it
    can be called freely from any context.
    """

    method_name: str = "SN"

    def realize(
        self,
        law: "BoundaryOperator",
        method_space: SNMethodSpace,
    ) -> LinearOperator:
        """Realize ``law`` for SN as a 1-arg :class:`LinearOperator`."""
        quad = method_space.quadrature

        if isinstance(law, VacuumBoundaryOperator):
            if method_space.inflow_indices is None:
                raise BoundaryError(
                    "SNBoundaryRealizer cannot realize "
                    "VacuumBoundaryOperator without inflow_indices in "
                    "method_space. Wave 8's SNMesh wiring populates "
                    "this from "
                    "InflowTraceSpace.inflow_indices_for_face. For "
                    "now, supply inflow_indices explicitly via "
                    "SNMethodSpace(quadrature=quad, face=..., "
                    "inflow_indices=...).",
                    law="vacuum",
                )
            return IncomingOrdinateMaskTensor(
                inflow_indices=method_space.inflow_indices,
                n_ordinates=quad.N,
                axis=0,
            )

        if isinstance(law, SpecularBoundaryOperator):
            perm = quad.reflection_index(law.axis)
            base = PermutationOperator(perm, axis=0)
            if law.albedo == 1.0:
                # Fast path: bare permutation preserves bit-identity
                # with legacy when albedo is exactly 1.0.
                return base
            return ScaledOperator(float(law.albedo), base)

        if isinstance(law, WhiteBoundaryOperator):
            base = AngularAverageOperator.from_quadrature(
                quad, law.axis, law.outward_sign,
            )
            if law.albedo == 1.0:
                return base
            return ScaledOperator(float(law.albedo), base)

        if isinstance(law, AlbedoBoundaryOperator):
            if law.albedo == 0.0:
                return ZeroOperator()
            if law.albedo == 1.0:
                return IdentityOperator()
            return ScaledOperator(float(law.albedo), IdentityOperator())

        if isinstance(law, PeriodicBoundaryOperator):
            return PeriodicWrapOperator()

        if isinstance(law, MixedBoundaryOperator):
            if not law.components:
                # Degenerate: empty composition → zero operator.
                return ZeroOperator()
            # Recursively realize each component, then sum via
            # the OperatorSum algebra (LinearOperatorMixin.__add__).
            result: Optional[LinearOperator] = None
            for coeff, primitive in law.components:
                realized = self.realize(primitive, method_space)
                if coeff == 1.0:
                    term: LinearOperator = realized
                else:
                    term = ScaledOperator(float(coeff), realized)
                result = term if result is None else (result + term)
            assert result is not None  # guarded by ``if not law.components``
            return result

        raise BoundaryError(
            f"SNBoundaryRealizer cannot realize "
            f"{type(law).__name__} — no isinstance dispatch case. "
            f"Available cases: VacuumBoundaryOperator, "
            f"SpecularBoundaryOperator, WhiteBoundaryOperator, "
            f"AlbedoBoundaryOperator, PeriodicBoundaryOperator, "
            f"MixedBoundaryOperator.",
            law=type(law).__name__,
        )
