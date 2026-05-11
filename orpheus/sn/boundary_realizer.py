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
* :class:`~orpheus.geometry.boundary.prescribed_inflow.PrescribedInflow(source)` →
  :class:`~orpheus.sn.angular_operator.IncomingSourceOperator(source)`
  — apply returns ``source.evaluate(probe_inflow_trace)``, ignoring
  the outgoing flux. The rank-0 affine BC (Wave 7 / C7.5).

Wave 11 removed the ``MixedBoundaryOperator`` composer from the
dispatch table. Rank-N compositions are now expressed via Wave-0
``OperatorSum``-algebra over already-realised leaves; the convenience
walker for ``BoundaryTraceLaw``-rooted expression trees lives in
:func:`orpheus.sn.boundary_realize.realize_recursively`.

References
----------

* ``.claude/plans/transient-giggling-cake.md`` Wave 5 — C5.3 (functional
  realizer brief), §16A.5 (vacuum-trace semantic correction), Wave 11
  (mixed-BC removal + tree-walking realisation).
* Grand Report v3 §16A.3 lines 2841–2860 (realizer-as-third-layer
  motivation), §16A.10 lines 3160-3175 (vacuum-trace correct
  representation).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from orpheus.geometry.boundary import (
    AlbedoBoundaryOperator,
    BoundaryError,
    BoundaryRealizerRegistry,
    PeriodicBoundaryOperator,
    PrescribedInflow,
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
from orpheus.sn.angular_operator import (
    AngularAverageOperator,
    IncomingSourceOperator,
)

# Wave 8 (C8.1): SNMethodSpace moved to its own dedicated module
# (orpheus.sn.method_space) so it can carry mesh + trace metadata
# without cluttering the realizer file. Re-export here so existing
# consumers ``from orpheus.sn.boundary_realizer import SNMethodSpace``
# keep working unchanged.
from orpheus.sn.method_space import SNMethodSpace

if TYPE_CHECKING:
    from orpheus.geometry.boundary import BoundaryOperator


__all__ = ["SNBoundaryRealizer", "SNMethodSpace"]


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

        if isinstance(law, PrescribedInflow):
            # Wave-7 addition: rank-0 affine source. The operator's
            # apply ignores the outgoing flux and returns
            # source.evaluate(probe_trace). Wave 8 will extend
            # IncomingSourceOperator with face / inflow-indices
            # plumbing once SNMethodSpace carries them.
            return IncomingSourceOperator(law.source)

        raise BoundaryError(
            f"SNBoundaryRealizer cannot realize "
            f"{type(law).__name__} — no isinstance dispatch case. "
            f"Available cases: VacuumBoundaryOperator, "
            f"SpecularBoundaryOperator, WhiteBoundaryOperator, "
            f"AlbedoBoundaryOperator, PeriodicBoundaryOperator, "
            f"PrescribedInflow. For rank-N compositions built via "
            f"Wave-0 OperatorSum/ScaledOperator algebra over "
            f"BoundaryTraceLaw leaves, use "
            f"orpheus.sn.boundary_realize.realize_recursively.",
            law=type(law).__name__,
        )
