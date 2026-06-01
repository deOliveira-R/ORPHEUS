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
  :class:`~orpheus.numerics.spaces.trace_space.TraceSpace`. **Semantic
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
  — apply returns ``source.evaluate(psi_out.shape)``, ignoring the
  outgoing flux. The rank-0 affine BC (Wave 7 / C7.5).

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
    TensorProductOperator,
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
                    "TraceSpace.inflow_indices_for_face. For "
                    "now, supply inflow_indices explicitly via "
                    "SNMethodSpace(quadrature=quad, face=..., "
                    "inflow_indices=...).",
                    law="vacuum",
                )
            # Wave T step T.1 — vacuum BC lift to 2-factor
            # TensorProductOperator, mirroring the D-B+1 specular pattern.
            # IncomingOrdinateMaskTensor acts on the ordinate axis
            # (axis=0); trailing axes (group, spatial / face) broadcast
            # through the IdentityOperator factor.  The TP fold reduces
            # to the bare mask's ``apply`` (IdentityOperator.apply
            # returns ``x`` unchanged), so bit-identity with the pre-T.1
            # single-axis form is preserved.
            return (
                IncomingOrdinateMaskTensor(
                    inflow_indices=method_space.inflow_indices,
                    n_ordinates=quad.N,
                    axis=0,
                )
                & IdentityOperator()
            )

        if isinstance(law, SpecularBoundaryOperator):
            perm = quad.reflection_index(law.axis)
            # Depth B step D-B+1 — first production tensor-network instance.
            # The grand-report §16A.10 BC decomposition is
            # ``B = G_patch ⊗ K_omega ⊗ K_g``; for specular reflection the
            # only non-trivial factor is the angular permutation. The
            # group axis (and any face axis a higher layer composes) are
            # identity. Pre-D-B+1 this was a single-axis
            # ``PermutationOperator(perm, axis=0)`` whose implicit numpy
            # broadcast played the role of ``I_group``; promoting to a
            # ``TensorProductOperator`` makes the algebra type-visible
            # so adjoint distributivity, composition distributivity,
            # and (downstream) ``assert_separable`` light up without
            # changing the matvec output.
            base = (
                PermutationOperator(perm, axis=0) & IdentityOperator()
            )  # 2-factor TensorProductOperator
            if law.albedo == 1.0:
                # Fast path: bare permutation TP preserves bit-identity
                # with legacy when albedo is exactly 1.0. The fold
                # ``IdentityOperator.apply(np.take(x, perm, axis=0))``
                # reduces to ``np.take(x, perm, axis=0)`` — same bytes.
                return base
            return ScaledOperator(float(law.albedo), base)

        if isinstance(law, WhiteBoundaryOperator):
            # Wave T step T.1 — white BC lift to 2-factor
            # TensorProductOperator, mirroring D-B+1.
            # ``AngularAverageOperator`` acts on the ordinate axis (axis 0);
            # IdentityOperator broadcasts on trailing axes (group,
            # spatial / face).  Bit-identity preserved because
            # IdentityOperator.apply returns x unchanged.
            base = (
                AngularAverageOperator.from_quadrature(
                    quad, law.axis, law.outward_sign,
                )
                & IdentityOperator()
            )
            if law.albedo == 1.0:
                return base
            return ScaledOperator(float(law.albedo), base)

        if isinstance(law, AlbedoBoundaryOperator):
            if law.albedo == 0.0:
                return ZeroOperator()
            if law.albedo == 1.0:
                return IdentityOperator()
            # Wave T step T.1 — albedo BC lift.  For α ∉ {0,1} the
            # action is uniform attenuation across all axes; lifting
            # the inner identity to a 2-factor TP (I & I) makes the
            # §16A.10 algebra type-visible while remaining a no-op at
            # the apply level (both IdentityOperator factors return
            # ``x`` unchanged; the ``ScaledOperator`` wrapper supplies
            # the α multiplication).
            return ScaledOperator(
                float(law.albedo),
                IdentityOperator() & IdentityOperator(),
            )

        if isinstance(law, PeriodicBoundaryOperator):
            # Wave T step T.1 — periodic BC lift to 2-factor
            # TensorProductOperator.  The current PeriodicWrapOperator
            # body is identity-with-copy (the SN sweep handles the
            # spatial wrap via its face-pair indexing); the TP wrap
            # makes the §16A.10 algebra type-visible without changing
            # the matvec output.  When PeriodicWrapOperator gains a
            # non-trivial spatial-pushforward (follow-up issue), the
            # second factor will carry that structure.
            return PeriodicWrapOperator() & IdentityOperator()

        if isinstance(law, PrescribedInflow):
            # Wave-7 addition: rank-0 affine source. The operator's
            # apply ignores the outgoing flux and returns
            # source.evaluate(psi_out.shape). Wave 8 will extend
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
