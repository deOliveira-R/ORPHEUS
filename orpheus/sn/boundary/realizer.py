r"""SN-specific BoundaryRealizer.

Dispatches by ``isinstance(law, ...)`` to the Wave-0 / Wave-1 primitives
that realise each :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`
concrete as a single-arg :class:`~orpheus.numerics.operator.LinearOperator`
consumable by the SN sweep + Krylov path.

The isinstance dispatch keys on the Grand Report v3 canonical
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` concretes
(:class:`VacuumInflow`, :class:`ReflectiveBoundary`,
:class:`WhiteBoundary`, :class:`AlbedoBoundary`,
:class:`PeriodicBoundary`, :class:`PrescribedInflow`). The deprecated
``*BoundaryOperator`` aliases (Wave 5 ship-state) were retired in
Wave O step O.4a.1.

The realisation map (law concrete → Wave-0 / Wave-1 primitive)
================================================================

Since campaign phase **B3.2** a realized SN law is typed
:math:`\Gamma_+ \to \Gamma_-` — it consumes the **outflow** half-trace and
produces the **inflow** half-trace, exactly as the affine form
:math:`\gamma_-\psi = R\,G\,\gamma_+\psi + q` says and as the diffusion arm
already did. The consumer composes ``ι₋ ∘ law ∘ γ₊``; nothing is computed and
then discarded. Realizing a law therefore needs the face's
:attr:`~orpheus.sn.mesh.method_space.SNMethodSpace.outflow_indices` (its
domain) alongside ``inflow_indices`` (its codomain) — face-orientation data
that a quadrature alone cannot supply, so ``SNMethodSpace.minimal`` no longer
suffices for the two laws below.

* :class:`~orpheus.geometry.boundary.vacuum.VacuumInflow` → the **zero map**
  :math:`\Gamma_+ \to \Gamma_-` (a :class:`~orpheus.numerics.operator.ZeroOperator`
  carrying both space hooks). Vacuum's whole content is :math:`R = 0`; with
  the domain narrowed there is nothing else to represent.

  Pre-B3.2 this was an :class:`IncomingOrdinateMaskTensor` — a full-face
  projector onto the OUTFLOW subspace whose preserved rows the consumer then
  discarded. Two campaign phases documented that survival as having "no
  consumer today"; the narrowing removes the question rather than answering
  it, because those rows are no longer in the operator's domain.
* :class:`~orpheus.geometry.boundary.reflective.ReflectiveBoundary(axis, albedo)` →
  ``albedo * PermutationOperator(local_perm)`` on the REDUCED ordinate axis,
  where ``local_perm = γ₊.to_local(reflection_index(axis)[inflow])`` — row
  :math:`j` of the image reads the mirror of the :math:`j`-th inflow ordinate,
  at that ordinate's position inside :math:`\Gamma_+`. (The remap must go
  through ``to_local``: on a slab the mirror REVERSES order, so a hand-written
  ``arange`` is wrong there — see the ``TraceRestrictionOperator`` docstring.)
  The ``albedo=1.0`` fast path returns the bare
  :class:`PermutationOperator` TP.
* :class:`~orpheus.geometry.boundary.white.WhiteBoundary(axis, outward_sign, albedo)` →
  ``albedo * AngularAverageOperator.from_quadrature(quadrature, axis, outward_sign)``
  (with the ``albedo=1.0`` fast path), narrowed at **B3.4a** to contract over
  :math:`\Gamma_+` and re-emit on :math:`\Gamma_-`. Two things dissolved with
  the codomain: the operator's private ``> 0.0`` outflow test (it now reads the
  one face-name → signed-projection primitive against ``TANGENTIAL_EPS``, so it
  can no longer disagree with the trace space), and the full-face broadcast that
  wrote outflow rows nobody read. The law's declared ``axis``/``outward_sign``
  is cross-checked against the installation face's :math:`\Gamma_+` — see
  :func:`_checked_angular_average`.
* :class:`~orpheus.geometry.boundary.albedo.AlbedoBoundary(0.0)` →
  :class:`~orpheus.numerics.operator.ZeroOperator`.
* :class:`~orpheus.geometry.boundary.albedo.AlbedoBoundary(1.0)` →
  :class:`~orpheus.numerics.operator.IdentityOperator`.
* :class:`~orpheus.geometry.boundary.albedo.AlbedoBoundary(α)` for α∉{0,1} →
  ``α * IdentityOperator()``.
* :class:`~orpheus.geometry.boundary.periodic.PeriodicBoundary` →
  :class:`~orpheus.numerics.operator.PeriodicWrapOperator`.
* :class:`~orpheus.geometry.boundary.prescribed_inflow.PrescribedInflow(source)` →
  :class:`~orpheus.sn.boundary.angular.IncomingSourceOperator` — apply ignores
  the outgoing flux and asks the source spec to fill :math:`|\Gamma_-|` rows.
  The rank-0 affine BC (Wave 7 / C7.5). **B3.4a** retired the inflow MASK
  (#52 / ERR-047) along with the codomain it corrected: the rows it zeroed are
  no longer emitted on, so :math:`q \in \Gamma_-` holds by typing rather than by
  erasure. Its unmasked companion branch retired with it — post-B3.2 a method
  space without face data cannot supply :math:`\gamma_+` either, so that
  fallback was already unreachable here.
* :class:`~orpheus.geometry.boundary.zero_flux.ZeroFluxBoundary` →
  **REFUSED** (:class:`~orpheus.geometry.boundary.BoundaryError`):
  the Dirichlet idealization :math:`\phi_\Gamma = 0` is the
  albedo-family member :math:`\mathcal{A} = -1`, and a negative
  angular inflow is unrepresentable in a non-negative :math:`\psi`.
  It realizes only under the diffusion realizer (#290 P3).

Realization is the certification point (#52): every dispatch first
fires the law's :meth:`~orpheus.geometry.boundary.BoundaryTraceLaw.assert_realizable`
aggregate (the §16A.12 universal invariants + each law's specific
checks — reflective's measure/involution/inflow→outflow table trio,
white/albedo's sub-Markov bound, the universal ERR-047 source-trace
certification), and the vacuum arm cross-checks the claimed
``inflow_indices`` against the face-name geometry (ERR-041). The
error catalog's recurring lesson operationalized: contracts are
checked where the law meets the trace, not discovered by downstream
balance defects.

Wave 11 removed the ``MixedBoundaryOperator`` composer from the
dispatch table. Rank-N compositions are now expressed via Wave-0
``OperatorSum``-algebra over already-realised leaves; the walker for
``BoundaryTraceLaw``-rooted expression trees is the method-blind
:func:`~orpheus.geometry.boundary.realize_recursively` (moved to
``geometry/boundary/`` at #290 P7b — pass
``realizer=SNBoundaryRealizer()`` at the leaves).

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

import numpy as np

from orpheus.geometry.boundary import (
    AlbedoBoundary,
    BoundaryError,
    BoundaryTraceLaw,
    PeriodicBoundary,
    PrescribedInflow,
    ReflectiveBoundary,
    VacuumAppliedToOutgoingTraceError,
    VacuumInflow,
    WhiteBoundary,
    ZeroFluxBoundary,
    stamp_boundary_role,
)
from orpheus.numerics.spaces.angular_trace_space import (
    TANGENTIAL_EPS,
    build_omega_dot_n,
)
from orpheus.numerics.operator import (
    IdentityOperator,
    LinearOperator,
    PeriodicWrapOperator,
    PermutationOperator,
    TraceRestrictionOperator,
    ZeroOperator,
)
from .angular import (
    AngularAverageOperator,
    IncomingSourceOperator,
)

if TYPE_CHECKING:
    from collections.abc import Callable

    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.mesh.method_space import SNMethodSpace


__all__ = ["SNBoundaryRealizer"]


def _zero_rows(n_rows: int) -> "Callable[[object], np.ndarray]":
    r"""A :class:`ZeroOperator` space hook emitting ``n_rows`` ordinates.

    The zero map between two DIFFERENT spaces must emit the zero of the space
    it lands in, not an echo of the one it came from — see
    :class:`~orpheus.numerics.operator.ZeroOperator`'s own note. Trailing axes
    (group, spatial) are carried through unchanged.

    The parameter is typed ``object`` because a hook CONSUMED by the operator
    is contravariant in its input: it must accept whatever the operator hands
    it. This mirrors ``_ray_source_zero`` in :mod:`orpheus.sn.solver`, the
    codebase's existing space-hook spelling.
    """
    def build(x: object) -> "np.ndarray":
        arr = np.asarray(x)
        return np.zeros((n_rows,) + arr.shape[1:], dtype=arr.dtype)

    return build


def _outflow_restriction(
    method_space: "SNMethodSpace", law_key: str,
) -> "TraceRestrictionOperator":
    r"""The face's :math:`\gamma_+` — a narrowed law's DOMAIN.

    Since B3.2 an SN boundary law is typed :math:`\Gamma_+ \to \Gamma_-`, so
    realizing one requires knowing which ordinates are outflow AT THIS FACE.
    That is face-orientation data, not quadrature data: a law realized without
    a face has no :math:`\Gamma_+` to restrict to, and the demand is
    structural rather than an implementation detail. (Pre-B3.2 the specular
    branch needed only ``quadrature.reflection_index``, which is why a
    quadrature-only ``SNMethodSpace.minimal`` sufficed for it.)

    Built here rather than fetched from the trace's cache because realization
    runs **once per face at mesh construction**, not per sweep — the cached
    accessor exists for the per-sweep consumer. Reading the method space's own
    field also keeps a hand-built space workable, which is the reason
    ``inflow_indices`` is a field too.
    """
    if method_space.outflow_indices is None:
        raise BoundaryError(
            f"SNBoundaryRealizer cannot realize {law_key!r} without "
            f"outflow_indices: since B3.2 a boundary law is typed Γ₊ → Γ₋, "
            f"and its DOMAIN is the outflow half-trace of a particular face — "
            f"which a quadrature alone cannot name. Construct via "
            f"SNMethodSpace.for_face(quadrature=..., face=..., trace=...), "
            f"which populates both half-traces, or pass outflow_indices "
            f"explicitly alongside inflow_indices.",
            law=law_key,
        )
    return TraceRestrictionOperator(
        np.sort(np.asarray(method_space.outflow_indices, dtype=np.intp)),
        n_total=method_space.quadrature.N,
        axis=0,
    )


def _checked_angular_average(
    law: "WhiteBoundary",
    quadrature: "Quadrature",
    method_space: "SNMethodSpace",
    gamma_out: "TraceRestrictionOperator",
) -> "AngularAverageOperator":
    r"""White's Lambertian kernel, with the law's ORIENTATION cross-checked
    against the face it is being installed on (the ERR-041 pattern).

    :class:`~orpheus.geometry.boundary.WhiteBoundary` declares its own
    ``axis`` / ``outward_sign``; the method space independently names the
    face. Those are two encodings of ONE orientation, and until B3.4a
    nothing compared them — a white law declared for ``x`` and installed on
    ``ymin`` averaged over the wrong hemisphere and reported nothing. The
    check is the same shape as the vacuum arm's (``ERR-041``: swapped
    inflow/outflow face annotation), against index SETS rather than sizes,
    because ``|Γ₊| == |Γ₋|`` on every quadrature in the tree makes a size
    comparison Mode-12 blind.

    On the canonical ``SNMesh.realize_boundary_law`` path both encodings
    derive from the same face label, so the guard is green by construction;
    it bites on hand-built method spaces and on a mis-declared law.
    """
    law_face = f"{law.axis}{'max' if law.outward_sign == +1 else 'min'}"
    omega_dot_n = build_omega_dot_n(quadrature, (law_face,))[0]
    law_outflow = np.flatnonzero(omega_dot_n > +TANGENTIAL_EPS)
    if not np.array_equal(law_outflow, gamma_out.indices):
        raise BoundaryError(
            f"WhiteBoundary declares axis={law.axis!r}, "
            f"outward_sign={law.outward_sign:+d} — i.e. the face "
            f"{law_face!r}, whose Γ₊ is {law_outflow.tolist()} — but it is "
            f"being installed where Γ₊ is "
            f"{np.asarray(gamma_out.indices).tolist()} "
            f"(face={method_space.face!r}). The Lambertian averages over "
            f"the hemisphere its OWN orientation names, so a mismatch "
            f"averages the wrong ordinates silently.",
            law="white",
        )
    return AngularAverageOperator.from_quadrature(
        quadrature, law.axis, law.outward_sign,
    )


class SNBoundaryRealizer:
    r"""Functional realizer for SN BCs.

    Dispatches by ``isinstance(law, ...)`` to map each
    :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` concrete
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
        law: "BoundaryTraceLaw",
        method_space: SNMethodSpace,
    ) -> LinearOperator:
        """Realize ``law`` for SN as a 1-arg :class:`LinearOperator`.

        The concrete return is always a :class:`LinearOperator`
        subclass (a generic numerics primitive —
        :class:`TensorProductOperator` / :class:`ScaledOperator` /
        :class:`IncomingSourceOperator` / …); the narrower mixin type
        (vs the bare :class:`LinearOperator` Protocol) lets the
        :func:`~orpheus.geometry.boundary.stamp_boundary_role`
        ``block_role`` instance-stamp type-check and keeps the
        realized-law primitives assignable nominally rather than
        through the (invariant) Protocol.
        """
        if isinstance(law, ZeroFluxBoundary):
            # ── REFUSAL AXIS: state-cone / sign ──────────────────────────
            # NOT an angular-resolution refusal and NOT a topological one —
            # see `bc-method-realizability` in orpheus.geometry.boundary. SN resolves
            # than diffusion does and couples faces just fine; what it cannot
            # do is hold this VALUE. The diffusion trace is a signed partial
            # current, so 𝒜 = −1 lives in its state space; the SN trace is an
            # angular flux with ψ ≥ 0, and a negative inflow is outside that
            # cone. The same law is realizable there and unrealizable here for
            # a reason that has nothing to do with either method's fidelity.
            #
            # #290 P3 (user ruling 3): φ_Γ = 0 is a diffusion-level
            # Dirichlet idealization — albedo 𝒜 = −1 in the
            # partial-current basis. A transport BC prescribes the
            # inflow trace from the outflow trace, and a NEGATIVE
            # angular inflow is unrepresentable in a non-negative ψ,
            # so SN refuses rather than realizing something else.
            raise BoundaryError(
                "SNBoundaryRealizer cannot realize ZeroFluxBoundary: "
                "the zero-flux Dirichlet condition (φ_Γ = 0) is the "
                "albedo-family member 𝒜 = −1 in the partial-current "
                "basis, and a negative angular inflow has no transport "
                "realization (ψ ≥ 0 ⟹ J± ≥ 0). Use VacuumInflow for "
                "the physical zero-incoming law (J⁻ = 0), or realize "
                "zero-flux with the diffusion realizer.",
                law="zero_flux",
            )

        quad = method_space.quadrature

        # #52 — the §16A.12 invariants fire HERE, at the seam where the
        # law meets the actual quadrature + trace data (the error
        # catalog's recurring lesson: check at realization, not only by
        # downstream balance). One aggregate call; each law extends the
        # base template with its own invariants (reflective's three
        # table checks, white/albedo's sub-Markov bound, the universal
        # ERR-047 source-trace certification). Certification is a LAW
        # contract: a non-law object (and a LawSum/LawScaled tree,
        # which composes laws but is not one) falls through to the
        # loud dispatch-failure raise below instead.
        if isinstance(law, BoundaryTraceLaw):
            law.assert_realizable(
                quad, inflow_indices=method_space.inflow_indices
            )

        if isinstance(law, VacuumInflow):
            if method_space.inflow_indices is None:
                raise BoundaryError(
                    "SNBoundaryRealizer cannot realize "
                    "VacuumInflow without inflow_indices in "
                    "method_space. Wave 8's SNMesh wiring populates "
                    "this from "
                    "AngularTraceSpace.inflow_indices_for_face. For "
                    "now, supply inflow_indices explicitly via "
                    "SNMethodSpace(quadrature=quad, face=..., "
                    "inflow_indices=...).",
                    law="vacuum",
                )
            # #52 / ERR-041 — trace-orientation guard. The claimed
            # inflow_indices and the face NAME are two independent
            # encodings of the same orientation (the annotation-swap
            # bug class desynchronizes them: indices built from the
            # OUTGOING set of this face). Cross-check against the
            # signed projection the face name alone implies via the
            # single face→normal primitive. On the canonical
            # SNMesh.realize_boundary_law path both encodings derive
            # from the trace space, so the guard is green by
            # construction; it bites on hand-built method spaces —
            # exactly where the annotation swap lives. A method space
            # with no face carries no independent orientation truth,
            # so the guard cannot fire there (documented escape:
            # quadrature-only spaces with explicit indices are the
            # caller's responsibility).
            if method_space.face is not None:
                omega_dot_n = build_omega_dot_n(
                    quad, (method_space.face,)
                )[0]
                claimed = np.asarray(
                    method_space.inflow_indices, dtype=np.intp
                )
                outgoing = omega_dot_n[claimed] > +TANGENTIAL_EPS
                if np.any(outgoing):
                    offending = claimed[outgoing]
                    raise VacuumAppliedToOutgoingTraceError(
                        f"VacuumInflow at face "
                        f"{method_space.face!r}: inflow_indices "
                        f"{offending.tolist()} are OUTGOING ordinates "
                        f"(Ω·n̂ > 0) at this face — the vacuum law "
                        f"prescribes γ₋ψ = 0 on the INCOMING trace "
                        f"only; zeroing outflow slots discards the "
                        f"flux the sweep just computed (ERR-041: "
                        f"swapped inflow/outflow face annotation).",
                        law="vacuum",
                    )
            # B3.2 — the domain is Γ₊, so vacuum is the honest ZERO MAP.
            #
            # Vacuum's whole content is ``R = 0``; its ``G`` is the identity
            # deck element (it fixes no geometry). With the law typed
            # ``Γ₊ → Γ₋`` there is nothing left to represent but the zero map
            # between those two spaces — and the question the old realization
            # raised, "why does a projector preserve rows nobody reads?",
            # dissolves rather than needing an answer: those rows are no
            # longer in the operator's domain at all.
            #
            # Pre-B3.2 this returned ``IncomingOrdinateMaskTensor(...) &
            # IdentityOperator()`` — a full-face projector onto the OUTFLOW
            # subspace, whose preserved rows ``_reflect_trace`` then discarded.
            # ``ZeroOperator``'s symmetric space hooks are the designed
            # mechanism for a genuine map BETWEEN spaces: the forward emits the
            # zero of Γ₋, the transpose the zero of Γ₊. Relying on the
            # endomorphic ``0.0 * x`` echo would be wrong in principle and
            # merely lucky in practice (|Γ₊| = |Γ₋| on every reachable
            # fixture, so the shapes coincide — an accident, not a contract).
            gamma_out = _outflow_restriction(method_space, "vacuum")
            return stamp_boundary_role(
                ZeroOperator(
                    # ``reportArgumentType``: ``ZeroOperator`` is generic over
                    # ``Vector``, and MEASURED (2026-07-31, pyright + numpy
                    # stubs) ``np.ndarray`` does not satisfy that protocol
                    # statically — ``Self``-typed ``__add__`` will not bind
                    # against numpy's overloads. The gap is upstream, not at
                    # this call site: every ndarray-level operator in
                    # ``numerics.operator`` sidesteps it by being declared
                    # UNPARAMETERIZED (``PermutationOperator(LinearOperator)``),
                    # an option a generic's own hooks do not have. Runtime
                    # conformance is real; only the static bind fails.
                    codomain_zero=_zero_rows(  # type: ignore[reportArgumentType]
                        method_space.inflow_indices.size
                    ),
                    transpose_zero=_zero_rows(  # type: ignore[reportArgumentType]
                        gamma_out.n_restricted
                    ),
                )
            )

        if isinstance(law, ReflectiveBoundary):
            # B3.2 — the specular mirror, narrowed to Γ₊ → Γ₋.
            #
            # The mirror sends the inflow ordinate ``inflow[j]`` to the outflow
            # ordinate ``perm[inflow[j]]`` (the ERR-045 invariant the law layer
            # certifies at realization: a reflection maps inflow to OUTFLOW,
            # never to itself). So the narrowed operator's row ``j`` reads the
            # position of that outflow ordinate INSIDE Γ₊ — which is exactly
            # what ``γ₊.to_local`` computes, and what a hand-written
            # ``arange`` would get wrong: on a slab the mirror REVERSES order,
            # so ``perm[inflow] = [3, 2]`` maps to local ``[1, 0]``.
            gamma_out = _outflow_restriction(method_space, "reflective")
            inflow = np.asarray(method_space.inflow_indices, dtype=np.intp)
            perm = quad.reflection_index(law.axis)
            if inflow.size != gamma_out.n_restricted:
                raise BoundaryError(
                    f"SNBoundaryRealizer cannot narrow ReflectiveBoundary at "
                    f"face {method_space.face!r}: |Γ₋| = {inflow.size} but "
                    f"|Γ₊| = {gamma_out.n_restricted}. A specular mirror is a "
                    f"BIJECTION between the two half-traces, so a face where "
                    f"they differ in size has no specular realization — which "
                    f"means the quadrature's tangential band has swallowed "
                    f"ordinates asymmetrically at this face.",
                    law="reflective",
                )
            # ``to_local`` raises if the mirror sent an inflow ordinate
            # anywhere but Γ₊ — the ERR-045 violation, caught here as a
            # crossed-index-set error rather than as silent wrong physics.
            local_perm = gamma_out.to_local(perm[inflow])
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
            # Depth B step D-B+1 — the §16A.10 BC decomposition
            # ``B = G_patch ⊗ K_omega ⊗ K_g``; for specular reflection the only
            # non-trivial factor is the angular permutation, now expressed on
            # the REDUCED ordinate axis. The group axis (and any face axis a
            # higher layer composes) are identity, which keeps adjoint and
            # composition distributivity type-visible.
            base = (
                PermutationOperator(local_perm, axis=0) & IdentityOperator()
            )  # 2-factor TensorProductOperator
            if law.albedo == 1.0:
                # Fast path: the bare permutation TP, whose fold
                # ``IdentityOperator.apply(np.take(x, local_perm, axis=0))``
                # reduces to the single gather.
                return stamp_boundary_role(base)
            return stamp_boundary_role(float(law.albedo) * base)

        if isinstance(law, WhiteBoundary):
            # B3.4a — the Lambertian kernel, narrowed to Γ₊ → Γ₋.
            #
            # ``AngularAverageOperator`` now derives BOTH half-traces from
            # the one face-name → signed-projection primitive, classified
            # against ``TANGENTIAL_EPS`` — so its old private ``> 0.0``
            # outflow test is gone, and with it the disagreement with the
            # trace space on every quadrature carrying tangential ordinates.
            gamma_out = _outflow_restriction(method_space, "white")
            base = (
                _checked_angular_average(law, quad, method_space, gamma_out)
                & IdentityOperator()
            )
            if law.albedo == 1.0:
                return stamp_boundary_role(base)
            return stamp_boundary_role(float(law.albedo) * base)

        if isinstance(law, AlbedoBoundary):
            if law.albedo == 0.0:
                return stamp_boundary_role(ZeroOperator())
            if law.albedo == 1.0:
                return stamp_boundary_role(IdentityOperator())
            # Wave T step T.1 — albedo BC lift.  For α ∉ {0,1} the
            # action is uniform attenuation across all axes; lifting
            # the inner identity to a 2-factor TP (I & I) makes the
            # §16A.10 algebra type-visible while remaining a no-op at
            # the apply level (both IdentityOperator factors return
            # ``x`` unchanged; the ``ScaledOperator`` wrapper supplies
            # the α multiplication).
            return stamp_boundary_role(
                float(law.albedo) * (IdentityOperator() & IdentityOperator())
            )

        if isinstance(law, PeriodicBoundary):
            # Wave T step T.1 — periodic BC lift to 2-factor
            # TensorProductOperator.  The current PeriodicWrapOperator
            # body is identity-with-copy (the SN sweep handles the
            # spatial wrap via its face-pair indexing); the TP wrap
            # makes the §16A.10 algebra type-visible without changing
            # the matvec output.  When PeriodicWrapOperator gains a
            # non-trivial spatial-pushforward (follow-up issue), the
            # second factor will carry that structure.
            return stamp_boundary_role(PeriodicWrapOperator() & IdentityOperator())

        if isinstance(law, PrescribedInflow):
            # B3.4a — the rank-0 affine source, narrowed to Γ₊ → Γ₋.
            #
            # The operator's apply ignores the outgoing flux and asks the
            # source to fill |Γ₋| rows. The inflow MASK (#52 / ERR-047)
            # dissolved with the codomain: the rows it used to zero — every
            # outflow and tangential ordinate — are no longer in the
            # codomain to be emitted on, so "q lives on Γ₋" is now a
            # property of the TYPE rather than of an erasure. Same
            # dissolution as vacuum's projector at B3.2.
            #
            # ``_outflow_restriction`` is called for its GUARD, not its
            # value: prescribed inflow ignores its input, but a law that
            # cannot name its own domain is not realized, it is guessed.
            _outflow_restriction(method_space, "prescribed_inflow")
            if method_space.inflow_indices is None:
                raise BoundaryError(
                    "SNBoundaryRealizer cannot realize PrescribedInflow "
                    "without inflow_indices: since B3.4a the delivered q "
                    "is sized from Γ₋ — the source spec fills whatever "
                    "shape it is handed, so the face must supply that "
                    "shape. Construct via SNMethodSpace.for_face(...).",
                    law="prescribed_inflow",
                )
            return IncomingSourceOperator(
                law.source,
                n_inflow=int(np.asarray(method_space.inflow_indices).size),
            )

        raise BoundaryError(
            f"SNBoundaryRealizer cannot realize "
            f"{type(law).__name__} — no isinstance dispatch case. "
            f"Available cases: VacuumInflow, "
            f"ReflectiveBoundary, WhiteBoundary, "
            f"AlbedoBoundary, PeriodicBoundary, "
            f"PrescribedInflow. For rank-N compositions built via "
            f"Wave-0 OperatorSum/ScaledOperator algebra over "
            f"BoundaryTraceLaw leaves, use "
            f"orpheus.geometry.boundary.realize_recursively with "
            f"realizer=SNBoundaryRealizer().",
            law=type(law).__name__,
        )


# The rank-N composition walker (``realize_recursively``) lived here
# until #290 P7b. With the second functional realizer shipped
# (diffusion, #290 P3) the walk's realizer-genericity became real and
# the walker moved to its method-blind home,
# :func:`orpheus.geometry.boundary.realize_recursively` — the leaf
# realizer is now a REQUIRED argument (``SNBoundaryRealizer()`` for the
# SN path). Single BCs never routed through it: ``SNMesh``'s
# ``realize_boundary_law`` arm (the ``TransportMethod`` hook) calls
# ``SNBoundaryRealizer().realize`` directly.
