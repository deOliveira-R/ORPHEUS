r"""The realization seam: BoundaryRealizer Protocol + the rank-N walker.

Per Grand Report v3 §16A.3 lines 2841–2860, the third layer of the
boundary architecture (after trace structure + physical law) is the
**method realisation**: a per-transport-method strategy that turns a
BC descriptor (:class:`BoundaryTraceLaw`) into a
:class:`~orpheus.numerics.operator.LinearOperator` consumable by that
method's sweep / solver.

The three-layer split matters because the same physical law (vacuum,
specular reflection, white, …) is realized by **structurally
different** linear operators in each transport method:

* SN realizes vacuum as a sparse per-ordinate
  :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor` on
  the inflow ordinates of the affected face;
* diffusion realizes the same vacuum law as the Marshak albedo row
  :math:`J^- = 0` (the :class:`~orpheus.numerics.operator.ZeroOperator`
  collapse);
* MoC would realize it as zeroing the entering track boundary fluxes;
  MC by killing particles that exit; CP as zero rows in the
  boundary-to-region coupling matrix; …

This module carries the three method-agnostic pieces of that seam:

* the :class:`BoundaryRealizer` **Protocol** — the shape every
  method realizer satisfies (``realize(law, method_space) →
  LinearOperator``). The concrete realizers ship in their method
  packages and are owned by their method-meshes:
  :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` and
  :class:`~orpheus.diffusion.boundary_realizer.DiffusionBoundaryRealizer`.
* :func:`stamp_boundary_role` — the shared helper every functional
  realizer applies to its realized operators so they carry
  :attr:`~orpheus.numerics.operator.BlockRole.BOUNDARY`.
* :func:`realize_recursively` — the **rank-N composition walker**,
  the sole type transformer from a descriptor tree
  (:class:`BoundaryTraceLaw` leaves under :class:`LawSum` /
  :class:`LawScaled` composers) to an operator tree
  (:class:`~orpheus.numerics.operator.OperatorSum` /
  :class:`~orpheus.numerics.operator.ScaledOperator` around realized
  1-arg leaves). Method-blind: the leaf realizer is a REQUIRED
  parameter — geometry knows the walk, never the method.

The registry that used to live here (a historical note)
=======================================================

Wave 5 shipped a ``BoundaryRealizerRegistry`` (string-keyed
``method_name → realizer``) plus ``NotImplementedError`` stub
realizers auto-registered by ``import orpheus.{moc,mc,cp}``, holding
the dispatch architecture "end-to-end" for methods that had not
adopted the unified BC model. #290 P7b **dissolved it**: no production
consumer ever resolved a realizer by name — you hold the method-mesh,
therefore you hold its realizer (each mesh's ``realize_boundary_law``
arm calls its own realizer directly, and the walker takes the realizer
as an argument). The string indirection carried registration-timing
hazards (a fresh process's empty registry, masked in-suite by
process-global state) for zero payoff; dissolving it deleted the
hazard class. The walker moved here from ``orpheus.sn.boundary``
in the same carve — see :mod:`orpheus.transport.method` for the
``TransportMethod`` Protocol minted over the method-meshes.

References
----------

* Grand Report v3 §16A.4 lines 2864–2876 (``realize(law, method_space)``
  vocabulary), §16A.3 (realizer-as-third-layer motivation).
* ``.claude/plans/diffusion_integration_290.md`` §P7b (the walker
  move + registry dissolution).
* :doc:`/theory/foundations/boundary_conditions` § "Rank-N boundaries".
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol, TypeVar, runtime_checkable

from orpheus.numerics.operator import BlockRole, OperatorSum, ScaledOperator

from ._base import BoundaryTraceLaw
from ._composition import LawScaled, LawSum

if TYPE_CHECKING:
    from orpheus.numerics.operator import LinearOperator

    from ._composition import LawNode


__all__ = [
    "BoundaryRealizer",
    "realize_recursively",
    "stamp_boundary_role",
]

#: The method space a realizer consumes. CONTRAVARIANT because it appears
#: only in parameter position (``realize(law, method_space)``): a realizer
#: accepting a WIDER space stands in wherever a narrower one is required.
#: Unbound on purpose — geometry does not know the method spaces (they live
#: in ``orpheus.sn`` / ``orpheus.diffusion``), and a bound here would invert
#: the dependency the P7b move established.
MethodSpaceT_contra = TypeVar("MethodSpaceT_contra", contravariant=True)

#: The same space at function scope, where variance does not apply. Its job
#: is to TIE ``method_space`` to ``realizer`` in :func:`realize_recursively`
#: so a mismatched pair (an ``SNMethodSpace`` handed to a
#: ``DiffusionBoundaryRealizer``) is a static error rather than a runtime
#: ``AttributeError`` deep inside a realizer arm.
MethodSpaceT = TypeVar("MethodSpaceT")


def stamp_boundary_role(op: "LinearOperator") -> "LinearOperator":
    r"""Stamp a realized boundary law with the :attr:`BlockRole.BOUNDARY` role.

    A realized law is a boundary-block leaf (``A_ss`` only) on the
    direct-sum transport state ``V = V_bulk ⊕ V_boundary`` — it maps the
    outflow trace to the inflow trace with no bulk action (Issue #208 /
    Wave O). The role lives on the INSTANCE (the realized op is a generic
    numerics primitive — :class:`TensorProductOperator`,
    :class:`ScaledOperator`, :class:`IdentityOperator`, … — that plays
    the BOUNDARY role only in this realization context), so it is stamped
    at the producer sites — the functional method realizers
    (:class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer`,
    :class:`~orpheus.diffusion.boundary_realizer.DiffusionBoundaryRealizer`)
    — through this ONE helper (``coding-elegance`` Pattern 7 / Pattern 2:
    the stamp is a shared realizer-layer concept, defined once). The
    rank-0 affine ``PrescribedInflow`` source is deliberately NOT
    stamped: it is the boundary *source* ``q.boundary``, not a linear
    boundary operator ``B``.
    """
    op.block_role = BlockRole.BOUNDARY
    return op


@runtime_checkable
class BoundaryRealizer(Protocol[MethodSpaceT_contra]):
    r"""Method-specific realisation of a boundary law.

    Implementors live in per-method packages
    (:class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer`,
    :class:`~orpheus.diffusion.boundary_realizer.DiffusionBoundaryRealizer`)
    and are owned by their method-meshes: each mesh's
    ``realize_boundary_law`` arm — the per-method hook of the
    :class:`~orpheus.transport.method.TransportMethod` Protocol —
    instantiates its own realizer directly, and rank-N composition
    passes one explicitly to :func:`realize_recursively`.

    The :meth:`realize` method takes a method-agnostic boundary law
    (the BC descriptor :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`)
    and a method-specific *method space* — a lightweight container
    holding whatever discretisation metadata the realizer needs
    (quadrature, mesh, trace masks, …) — and returns a
    :class:`~orpheus.numerics.operator.LinearOperator` whose 1-arg
    :meth:`apply(psi)` maps the outflow trace to the inflow trace.

    The intent of the Protocol is **structural** (the descriptor →
    operator shape), not algorithmic — every realizer's :meth:`realize`
    body is method-specific.

    Generic in the method space
    ---------------------------
    The Protocol is parameterised by the method space it consumes, so
    :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` satisfies
    ``BoundaryRealizer[SNMethodSpace]`` and
    :class:`~orpheus.diffusion.boundary_realizer.DiffusionBoundaryRealizer`
    satisfies ``BoundaryRealizer[DiffusionMethodSpace]``. A *common-member*
    Protocol would have been the wrong shape: the two spaces share only
    ``face: str | None`` by type (``mesh`` and ``trace`` are method-specific
    types, and ``quadrature`` / ``inflow_indices`` are SN-only), so a
    structural intersection would describe neither realizer honestly.

    The parameter is what lets :func:`realize_recursively` REQUIRE that the
    space and the realizer agree. Before it, both were ``Any`` — looser than
    either implementation, which already annotate their own concrete types —
    so handing an ``SNMethodSpace`` to a ``DiffusionBoundaryRealizer``
    type-checked and failed only at runtime.

    ``isinstance`` still works against the UNSUBSCRIPTED name
    (``isinstance(r, BoundaryRealizer)``, as the realizer conformance tests
    do); subscripted runtime checks are rejected by Python itself.

    Attributes
    ----------
    method_name : str
        Stable identifier (``"SN"``, ``"diffusion"``, …) naming the
        method in diagnostics. Never used for dispatch.
    """

    method_name: str

    def realize(
        self,
        law: "BoundaryTraceLaw",
        method_space: MethodSpaceT_contra,
    ) -> "LinearOperator":
        """Return a 1-arg :class:`LinearOperator` realising ``law`` for this method."""
        ...


# ─────────────────────────────────────────────────────────────────────
# Rank-N composition walker — descriptor tree → operator tree.
#
# Production realizes single BCs directly (each method-mesh's
# ``realize_boundary_law`` → its realizer); this walker is the rank-N
# composition entry point (the Marshak ``0.3·Reflective + 0.7·White``
# partial-current BC), exercised by the law-composition wall
# ``tests/geometry/test_law_composition.py``. Moved here from
# ``orpheus.sn.boundary`` at #290 P7b, when the second functional
# realizer (diffusion, #290 P3) made the walk's realizer-genericity
# real: the walk (LawSum/LawScaled → OperatorSum/ScaledOperator) is
# method-blind, and the leaf realizer is the caller's to supply.
# ─────────────────────────────────────────────────────────────────────


def realize_recursively(
    node: "LawNode",
    method_space: MethodSpaceT,
    realizer: "BoundaryRealizer[MethodSpaceT]",
) -> "LinearOperator":
    r"""Walk a descriptor tree and realise it as a Wave-0 operator tree.

    The descriptor tree's leaves are
    :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` instances;
    internal nodes are
    :class:`~orpheus.geometry.boundary._composition.LawSum` /
    :class:`~orpheus.geometry.boundary._composition.LawScaled`. The
    output tree has the same shape with
    :class:`~orpheus.numerics.operator.OperatorSum` /
    :class:`~orpheus.numerics.operator.ScaledOperator` composers and
    1-arg :class:`~orpheus.numerics.operator.LinearOperator` leaves
    realised by ``realizer``.

    This is the ONLY function that transforms a descriptor into an
    operator. The §16A.3 three-layer split is type-checkable through
    this signature: input in the law-type family, output in the
    operator-type family.

    Parameters
    ----------
    node
        Descriptor-tree root. Must be a
        :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` leaf or a
        :class:`~orpheus.geometry.boundary._composition.LawSum` /
        :class:`~orpheus.geometry.boundary._composition.LawScaled`
        composer.
    method_space
        Method space passed verbatim to ``realizer.realize`` at every
        leaf. Its type is the realizer's:
        :class:`~orpheus.sn.mesh.method_space.SNMethodSpace` for the
        SN realizer,
        :class:`~orpheus.diffusion.method_space.DiffusionMethodSpace`
        for the diffusion realizer. Geometry stays method-blind — the
        walker never inspects it.
    realizer
        The leaf realizer (a :class:`BoundaryRealizer`). REQUIRED:
        there is no method-agnostic default — geometry does not know
        the methods. (The pre-P7b SN default died with the walker's
        move out of ``orpheus.sn``.)

    Returns
    -------
    :class:`~orpheus.numerics.operator.LinearOperator`
        The realised 1-arg operator with the composer structure of
        ``node`` preserved (every leaf replaced by its realisation).

    Raises
    ------
    TypeError
        If a node is not a
        :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`,
        :class:`~orpheus.geometry.boundary._composition.LawSum`, or
        :class:`~orpheus.geometry.boundary._composition.LawScaled`. The
        error names the offending type.

    Examples
    --------
    Realise the standard Marshak mixed boundary for SN::

        from orpheus.geometry.boundary import (
            ReflectiveBoundary, WhiteBoundary, realize_recursively,
        )
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer
        from orpheus.sn.mesh.method_space import SNMethodSpace
        from orpheus.numerics.quadrature import Quadrature

        ms = SNMethodSpace.minimal(Quadrature.gauss_legendre(4))
        law = (
            0.3 * ReflectiveBoundary(axis="x")
            + 0.7 * WhiteBoundary(axis="x", outward_sign=+1)
        )
        realised = realize_recursively(law, ms, SNBoundaryRealizer())
        # OperatorSum of ScaledOperators around realised leaves.

    The same tree under the diffusion realizer (every leaf collapses
    to its albedo-family scalar on the partial-current trace)::

        from orpheus.diffusion.boundary_realizer import DiffusionBoundaryRealizer
        from orpheus.diffusion.method_space import DiffusionMethodSpace

        op = realize_recursively(
            law, DiffusionMethodSpace.minimal(), DiffusionBoundaryRealizer(),
        )
    """
    if isinstance(node, BoundaryTraceLaw):
        # Leaf: dispatch through the caller's realizer.
        return realizer.realize(node, method_space)

    if isinstance(node, LawScaled):
        # Recurse on the inner descriptor, wrap with ScaledOperator.
        inner_op = realize_recursively(node.inner, method_space, realizer)
        return ScaledOperator(node.scalar, inner_op)

    if isinstance(node, LawSum):
        # Recurse on both operands, reassemble via OperatorSum.
        a_op = realize_recursively(node.a, method_space, realizer)
        b_op = realize_recursively(node.b, method_space, realizer)
        return OperatorSum(a_op, b_op)

    raise TypeError(
        f"realize_recursively expected BoundaryTraceLaw | LawSum | "
        f"LawScaled (the descriptor-tree type family), got "
        f"{type(node).__name__}. Operator-tree composers (e.g. "
        f"OperatorProduct, TensorProductOperator) are not valid "
        f"inputs — realize each descriptor first, then compose the "
        f"results via Wave-0 operator algebra."
    )
