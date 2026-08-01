r"""Diffusion BoundaryRealizer — the albedo-family realization.

Functional realizer for diffusion (#290 P3; closes issue #182 — the
Wave-5 stub this file carried until 2026-07-03). Owned by its
method-mesh:
:meth:`DiffusionMesh.realize_boundary_law
<orpheus.diffusion.augmented_mesh.DiffusionMesh.realize_boundary_law>`
— the diffusion arm of the
:class:`~orpheus.transport.method.TransportMethod` hook — instantiates
it directly (#290 P7b dissolved the Wave-5 realizer registry).

Every diffusion boundary condition is an albedo
================================================

Diffusion's boundary state is the scalar partial-current pair
``(J⁺, J⁻)`` on the
:class:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace`
— ONE outflow DOF and ONE inflow DOF per face per group (#290 ruling
2). A linear homogeneous boundary law therefore collapses to a single
scalar response per face per group:

.. math::

    J^- \;=\; \mathcal{A}\, J^+ ,

and :meth:`DiffusionBoundaryRealizer.realize` is the composition of
two named maps: *law → albedo scalar* 𝒜 (the physics table below),
then *albedo scalar → operator* (the structure-keyed collapse
mirroring the SN realizer's ``AlbedoBoundary`` branch). The realized
operator maps the outflow half of a face slot to the inflow half; it
is stamped :attr:`~orpheus.numerics.operator.BlockRole.BOUNDARY` via
:func:`~orpheus.geometry.boundary.stamp_boundary_role` (the ``A_ss``
leaf of the composite state, consumed by the P4
``DiffusionBoundaryOperator`` assembly).

The realization table (law → 𝒜 → operator)
===========================================

.. note::

   **The table is derived, not transcribed** (campaign phase B2.2).
   The first map is one line — ``𝒜 = law.response_kernel.amplitude`` —
   not the five-arm ``isinstance`` ladder it used to be. Every row
   below then falls out of the response factor each law already
   declares, and the reason the four value rows looked alike is now
   stated rather than coincidental: **at P1 the angular geometry map
   is integrated out of the half-range moments**, so a specular
   mirror, a Lambertian re-emission and an identity deck element all
   leave the same scalar on the partial-current trace. (This is the
   tier-2 case of :ref:`bc-method-realizability`: the realization is
   EXACT, and the projection is simply not FAITHFUL — P1 cannot tell
   those laws apart, which is a different statement from approximating
   any of them. A "null map" is not in the list because B3.0 retired
   it: the zero map is not a bijection and so was never a geometry.) The two refusals are
   exactly the laws whose remaining factor is *not* a per-face scalar
   — a spatial relabeling and an affine source. Read the table as the
   physics; read ``response_kernel`` as the single source it is
   computed from.

===============================  =======  =========================================
Law                              𝒜        Realized operator
===============================  =======  =========================================
``VacuumInflow``                 ``0``    :class:`~orpheus.numerics.operator.ZeroOperator`
``ReflectiveBoundary(α)``        ``α``    ``IdentityOperator`` (α=1) / ``α·I``
``WhiteBoundary(α)``             ``α``    ``IdentityOperator`` (α=1) / ``α·I``
``AlbedoBoundary(α)``            ``α``    ``Zero`` (α=0) / ``I`` (α=1) / ``α·I``
``ZeroFluxBoundary``             ``−1``   ``ScaledOperator(-1, IdentityOperator)``
``PeriodicBoundary``             —        REFUSED (needs the opposite-face wrap)
``PrescribedInflow``             —        REFUSED (affine source, not an operator)
===============================  =======  =========================================

* **Vacuum is Marshak** (#290 ruling 3): :math:`\mathcal{A} = 0` IS
  the Marshak zero-incoming-current condition :math:`J^- = 0`, which
  in :math:`(\phi, J)` variables reads
  :math:`\phi_\Gamma/4 - J/2 = 0`, i.e. the classic Robin form
  :math:`\phi + 2D\,\nabla\phi\cdot\hat n = 0` under Fick's law
  (Marshak 1947; Bell & Glasstone 1970 §3.4). The zero-flux Dirichlet
  condition the legacy island mis-registered as ``"vacuum"`` is a
  DIFFERENT law — :class:`~orpheus.geometry.boundary.ZeroFluxBoundary`,
  :math:`\mathcal{A} = -1`.
* **White coincides with reflective at P1** (documented per the plan):
  specular and Lambertian return differ only in the ANGULAR
  redistribution of the returned particles, which the half-range
  ℓ=0 moments integrate out — both preserve the returned current,
  :math:`J^- = \alpha J^+`. The distinction is real in transport
  (SN realizes them as a permutation vs a cosine-weighted average)
  and vanishes in any P1-closed method by construction.
* **Zero-flux is deliberately outside the sub-Markov range**
  :math:`[0, 1]`: partial-current positivity :math:`J^\pm \ge 0` is
  a property of the PHYSICAL laws, not a type invariant — see
  :class:`~orpheus.geometry.boundary.ZeroFluxBoundary`.

Refusals
========

* :class:`~orpheus.geometry.boundary.PeriodicBoundary` couples the
  face to its OPPOSITE face (``J⁻(xmin) = J⁺(xmax)`` and vice versa)
  — not a per-face rank-1 albedo but a trace-block permutation. No
  diffusion consumer exists today; the realizer refuses loudly rather
  than realizing the wrong (identity-like) thing. The wrap lands with
  a consumer, on the P4 boundary-operator assembly where the whole
  trace block is in scope.
* :class:`~orpheus.geometry.boundary.PrescribedInflow` is the rank-0
  AFFINE law :math:`J^- = q`: its realization is the boundary
  *source* ``q.boundary``, not a linear boundary operator ``B``
  (exactly the SN split — SN realizes it as an
  ``IncomingSourceOperator`` and does NOT stamp it BOUNDARY). The
  diffusion source arm lands with the P5 fixed-source wiring; until
  then the refusal keeps the operator/source split honest.

Rank-N composition (Marshak mixes, partial reflection) goes through
the descriptor-tree walker with this realizer at the leaves::

    from orpheus.geometry.boundary import realize_recursively
    op = realize_recursively(
        0.3 * ReflectiveBoundary(axis="x") + 0.7 * AlbedoBoundary(albedo=0.5),
        DiffusionMethodSpace.minimal(),
        realizer=DiffusionBoundaryRealizer(),
    )

(The walker generalized over leaf realizers at #290 P3 and moved to
its method-blind home in ``geometry/boundary/`` at the #290 P7b
``TransportMethod`` mint.)

References
----------

* Marshak, R. E. (1947). "Note on the spherical harmonic method as
  applied to the Milne problem for a sphere". *Phys. Rev.* 71, 443.
* Bell, G. I. & Glasstone, S. (1970). *Nuclear Reactor Theory*.
  Van Nostrand Reinhold. §3.4 (diffusion-theory boundary conditions).
* ``.claude/plans/diffusion_crosswalk.md`` — the ``(J⁺, J⁻)``
  convention contract (rows, signs, the BC-family table this module
  implements).
"""

from __future__ import annotations

from orpheus.geometry.boundary import (
    BoundaryError,
    BoundaryTraceLaw,
    PrescribedInflow,
    SpatialWrap,
    stamp_boundary_role,
)
from orpheus.numerics.operator import (
    IdentityOperator,
    LinearOperator,
    ZeroOperator,
)

from .method_space import DiffusionMethodSpace


__all__ = ["DiffusionBoundaryRealizer"]


def _albedo_operator(albedo: float) -> LinearOperator:
    r"""The structure-keyed collapse of :math:`J^- = \mathcal{A} J^+`.

    Mirrors the SN realizer's ``AlbedoBoundary`` branch exactly:
    :math:`\mathcal{A} = 0` collapses to the structural
    :class:`ZeroOperator`, :math:`\mathcal{A} = 1` to the structural
    :class:`IdentityOperator`, anything else to
    ``ScaledOperator(𝒜, IdentityOperator)``. The operator is
    face-slot-shape-agnostic — a scalar multiple of the identity on
    the ``(ng, *face_spatial)`` outflow data, whatever the face.
    """
    if albedo == 0.0:
        return ZeroOperator()
    if albedo == 1.0:
        return IdentityOperator()
    return float(albedo) * IdentityOperator()


class DiffusionBoundaryRealizer:
    r"""Functional realizer for diffusion BCs (#290 P3, closes #182).

    Realizes every :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`
    as the albedo-family operator :math:`J^- = \mathcal{A} J^+` on the
    scalar partial-current trace — see the module docstring for the
    law → 𝒜 table, the Marshak-vacuum and white≡reflective physics
    notes, and the two refusals.

    The realizer instance carries no state; ``method_name`` is the
    only attribute. :meth:`realize` is stateless — it can be called
    freely from any context.
    """

    method_name: str = "diffusion"

    def realize(
        self,
        law: "BoundaryTraceLaw",
        method_space: DiffusionMethodSpace,
    ) -> LinearOperator:
        r"""Realize ``law`` for diffusion as a 1-arg :class:`LinearOperator`.

        The returned operator maps the outflow half of a face slot
        (``J⁺``, shape ``(ng, *face_spatial)``) to the inflow half
        (``J⁻``), stamped
        :attr:`~orpheus.numerics.operator.BlockRole.BOUNDARY`. The
        rank-1 albedo family reads nothing from ``method_space``
        (mirroring SN's albedo branch — a
        :meth:`DiffusionMethodSpace.minimal` suffices); the parameter
        is the uniform :class:`~orpheus.geometry.boundary.BoundaryRealizer`
        Protocol surface the P4/P5 wiring threads per-face spaces
        through.
        """
        albedo = self._partial_current_albedo(law)
        return stamp_boundary_role(_albedo_operator(albedo))

    @staticmethod
    def _partial_current_albedo(law: "BoundaryTraceLaw") -> float:
        r"""The law's albedo-family response 𝒜 in :math:`J^- = \mathcal{A} J^+`.

        The physics table (module docstring) as code. Note the two
        DISTINCT albedos in play for reflective/white: the law's own
        ``albedo`` field α is the return AMPLITUDE (specular resp.
        diffuse), and at P1 the partial-current response equals it,
        :math:`\mathcal{A} = \alpha` — because the angular structure
        distinguishing the two return patterns is integrated out of
        the half-range moments.

        **𝒜 IS the law's response factor** — one read, not a table.
        Until campaign phase B2 this method was a five-arm
        ``isinstance`` ladder returning ``0.0`` / ``-1.0`` /
        ``float(law.albedo)``; every arm returned exactly
        ``law.response_kernel.amplitude``, which is why the B1 gate could
        pin the two bit-identical law by law
        (``tests/geometry/test_boundary_factors.py``). Collapsing it
        states the reason the arms looked alike: **at P1 the angular
        geometry map is integrated out of the half-range moments**, so
        a specular mirror, a Lambertian average, an identity and a null
        map all leave the same scalar on the partial-current trace, and
        only the response survives.

        Which is also exactly why the two refusals below are the
        refusals: they are the laws whose remaining factor is NOT a
        per-face scalar — a *spatial* relabeling (which couples two
        faces and so cannot be integrated away) and an affine source
        (which is not a linear response at all).
        """
        # The fallthrough guard the ``isinstance`` ladder used to provide as
        # its last arm. Collapsing the ladder onto ``response_kernel`` moved
        # the failure mode: an object that is not a law — or a
        # ``BoundaryTraceLaw`` subclass that never populated its factors, which
        # the ABC still permits (the default is ``None``) — used to fall
        # through to a named ``BoundaryError`` and would otherwise now die on a
        # bare ``AttributeError`` deep in the read. Restored here, first,
        # because everything below dereferences a factor.
        response = getattr(law, "response_kernel", None)
        if response is None:
            raise BoundaryError(
                f"DiffusionBoundaryRealizer cannot realize "
                f"{type(law).__name__} — it declares no `response_kernel`, so "
                f"there is no 𝒜 to realize. Diffusion realizes a law THROUGH "
                f"that factor: VacuumInflow (𝒜=0, Marshak), ReflectiveBoundary "
                f"(𝒜=α), WhiteBoundary (𝒜=α, P1-coincident with reflective), "
                f"AlbedoBoundary (𝒜=α), ZeroFluxBoundary (𝒜=−1). For rank-N "
                f"compositions use "
                f"orpheus.geometry.boundary.realize_recursively with "
                f"realizer=DiffusionBoundaryRealizer().",
                law=getattr(type(law), "key", None) or type(law).__name__,
            )

        # The one surviving TYPE test, and it is essential rather than a tag
        # smell: what disqualifies the prescribed family is that its inflow is
        # a free parameter q, NOT a function of the outflow — which is true of
        # the FAMILY whatever q currently holds. Testing the source VALUE
        # instead (``isinstance(law.source, NoSource)``) would silently realize
        # ``PrescribedInflow()`` at its default zero source as 𝒜 = 0, i.e.
        # turn a refusal into a value on a path (#290 P5) that does not exist
        # yet. Measured while writing B2 — it is the trap in this collapse.
        if isinstance(law, PrescribedInflow):
            # ── REFUSAL AXIS: none — this one is PLUMBING ────────────────
            # Unlike the two structural refusals (SN's zero-flux, on the
            # state-cone axis; the spatial wrap below, on the topological
            # one), nothing about P1 makes a prescribed inflow
            # unrepresentable: q is a VECTOR in Γ₋, and Γ₋ here is one
            # number per face, which the trace carries perfectly well. This
            # refusal disappears the day #290 P5 wires the fixed-source arm,
            # with no theory changing. Stated so a future reader does not
            # mistake it for a statement about diffusion's fidelity.
            raise BoundaryError(
                f"DiffusionBoundaryRealizer cannot realize "
                f"{type(law).__name__}: J⁻ = q is the rank-0 AFFINE law — "
                f"its realization is the boundary source q.boundary, "
                f"not a linear boundary operator B (the same "
                f"operator/source split SN keeps). The diffusion "
                f"fixed-source arm lands with the solver wiring "
                f"(#290 P5).",
                law=type(law).key or type(law).__name__,
            )
        if isinstance(law.geometry_map, SpatialWrap):
            # ── REFUSAL AXIS: spatial / topological ──────────────────────
            # NOT an angular-resolution refusal, which is the one a reader
            # expects from a diffusion realizer. Periodic's deck
            # transformation is a TRANSLATION: it acts on the spatial
            # coordinate and leaves angle untouched, so every angular
            # discretization is trivially equivariant under it and P1 loses
            # nothing angularly. What blocks it is that the wrap couples a
            # face to its OPPOSITE face, and this realizer's codomain is a
            # per-face scalar 𝒜 — a block-diagonal object with no slot for
            # cross-face coupling. A method could resolve angle perfectly and
            # still refuse this.
            raise BoundaryError(
                f"DiffusionBoundaryRealizer cannot realize "
                f"{type(law).__name__}: its geometry map is a SPATIAL "
                f"wrap, which couples a face to its OPPOSITE face (a "
                f"trace-block permutation) — the one geometry P1 cannot "
                f"integrate away into a per-face albedo J⁻ = 𝒜·J⁺. The "
                f"wrap lands with the boundary-operator assembly when a "
                f"diffusion consumer exists (#290 P4 seam).",
                law=type(law).key or type(law).__name__,
            )
        return response.amplitude
