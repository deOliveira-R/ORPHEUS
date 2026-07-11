r"""The realized boundary law ``B`` as a whole-trace BOUNDARY-block operator.

Wave O (Issue #208) extracts the boundary conditions from the streaming
operator ``L`` so that ``(L_full + C − S − F − B)\psi = q`` is the canonical
transport algebra. ``B`` is the **boundary-block** operator on the direct-sum
transport state ``V = V_bulk ⊕ V_boundary``: a 2×2 block matrix with only the
``A_ss`` (boundary → boundary) block non-zero. It maps the **outflow** trace to
the **inflow** trace via the per-face realized boundary laws (reflective /
vacuum / white / albedo / periodic), with **no bulk action**.

Block structure
===============

On ``V = V_bulk ⊕ V_boundary`` the four operator families are::

    C, S, F  →  [ A_bb  0 ]   (BULK   — bulk → bulk only)
                [ 0     0 ]

    L_full   →  [ A_bb  A_bs ] (FULL   — streaming couples bulk ↔ boundary)
                [ A_sb  0    ]

    B        →  [ 0     0   ] (BOUNDARY — boundary → boundary only, ``A_ss``)
                [ 0     A_ss]

The whole boundary block is the direct sum ``B = B_a + B_b`` of the two
per-system boundaries (a block-composed system's boundary is the direct sum of
its per-system boundary blocks over the composite biproduct):
:class:`SNBoundaryOperator` (``B_a``) is System A's **trace** ``A_ss`` leaf, and
:class:`RadialCharacteristicBoundaryOperator` (``B_b``) is System B's **ψ½
ray-corner** ``A_ss`` leaf (present only on a seed-carrying mesh). As a sibling
``−B`` of ``L`` in the :class:`~orpheus.numerics.operator.OperatorSum` algebra it
supplies the reflective coupling that ``L`` previously absorbed inside its sweep
(the ``inflow = bc.apply(outflow)`` re-apply); the outer Krylov / SI loop then
drives the boundary **consistency residual** ``ψ.inflow − B·ψ.outflow − q.inflow → 0``.

Construction
============

The per-face boundary laws already live on the
:class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` in the face-name-keyed ``bc`` dict
(each entry a :class:`~orpheus.geometry.boundary._bound_compat._BoundBoundaryOperator`
wrapping a realized law that carries :attr:`BlockRole.BOUNDARY`). The whole-trace
``B`` is the block-diagonal composition over the mesh's true boundary faces: for
each face present in the trace it applies that face's law to that face's slot.
``B`` is therefore block-diagonal over faces — it never mixes faces.

See :ref:`operator-algebra` and the Wave O plan
``.claude/plans/wave_o_operator_typing.md`` (step O.4a.2).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple, Optional

import numpy as np

from orpheus.numerics.operator import (
    BlockRole,
    LinearOperator,
    MissingAdjoint,
    SystemRole,
    adjointable,
)

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping

    from orpheus.numerics.space import FunctionSpace
    from orpheus.sn.loss_representation.sweep_schedule import SweepSchedule
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.transport.fields._bases import (
        RadialCharacteristicBoundaryField,
        RadialCharacteristicField,
    )
    from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
    from orpheus.transport.full_field import FullField
    from orpheus.transport.radial_characteristic_composite import (
        RadialCharacteristicComposite,
    )
    from orpheus.transport.source_sinks import (
        AngularBoundarySourceSink,
        RadialCharacteristicBoundarySourceSink,
    )


__all__ = [
    "BoundarySplit",
    "RadialCharacteristicBoundaryOperator",
    "SNBoundaryOperator",
    "SNMaskedBoundaryOperator",
]

#: The outer-face law kinds with a RULED ψ½ corner action (RULING P1's ray
#: carrier). Single source for BOTH
#: :attr:`RadialCharacteristicBoundaryOperator.is_adjointable` (a ruled kind's
#: corner map is Euclidean-adjointable) AND
#: :meth:`RadialCharacteristicBoundaryOperator._reflect_corner` (an unruled kind
#: is loud-deferred). white / albedo / periodic are absent — their off-quadrature
#: μ = ±1 re-emission is a design decision not yet ruled (2.5d plan-of-record);
#: add a kind here AND its branch in ``_reflect_corner`` when it is ruled.
_RULED_CORNER_KINDS = frozenset({"vacuum", "reflective"})


def _zero_bulk_source(mesh: "SNMesh"):
    r"""The zero-bulk ``A_ss`` carrier every boundary block emits.

    Sized from the MESH (not ``zeros_like`` the input) so the carrier is correct
    for any bulk representation (full-angular
    :class:`~orpheus.transport.fields.angular_flux.AngularFlux` OR the Phase-5a
    windowed :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`);
    it carries the scheme's spatial-moment width (#240 D5b-S3) so it composes
    element-wise with the moment-carrying ``(L+C)ψ`` in the ``(L+C − S − B)ψ``
    matvec.  Single source (Cardinal Rule 2) for the zero bulk both ``B_a``
    (:class:`SNBoundaryOperator`) and ``B_b``
    (:class:`RadialCharacteristicBoundaryOperator`) emit.
    """
    from orpheus.transport.source_sinks import AngularSourceSink

    return AngularSourceSink.zeros_on(
        mesh, spatial_moments=mesh.scheme.spatial_basis_per_axis,
    )


class SNBoundaryOperator(LinearOperator):
    r"""``B_a`` — System A's (trace) boundary law, the SN ``A_ss`` block.

    The boundary operator of the transport system (System A of the 2×2 coupled
    block operator — bulk⊕trace): block-diagonal over the mesh's true boundary
    faces, ``B_a.apply(ψ)`` returns a
    :class:`~orpheus.transport.full_field.FullField` with **zero bulk**
    and, on each face, ``bc[<face>].apply(ψ.boundary.face_view(<face>))`` — the
    per-face realized boundary law applied to that face's trace slot. It composes
    as ``−B`` in ``(L_full + C − S − F − B)`` (acting on the same
    :class:`~orpheus.transport.full_field.FullField` carrier as ``L``/``C``/``S``/``F``).

    On a **seed-carrying** composite the whole boundary block is the direct sum
    ``B = B_a + B_b``, where ``B_b`` is the sibling
    :class:`RadialCharacteristicBoundaryOperator` — System B's (ψ½ ray-corner)
    boundary (RULING P1: a block-composed system's boundary is the direct sum of
    per-system boundary blocks over the composite biproduct; the off-diagonal
    structure is keyed to face physics — reflection is a per-system
    endomorphism ⇒ block-diagonal). ``B_a`` emits a **present-zero** ray block so
    ``B_a + B_b`` sums bit-identically; the ray corner is entirely ``B_b``'s.

    The role is :attr:`BlockRole.BOUNDARY`; the domain and codomain are the
    mesh's composite carrier
    :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`
    (``sn_mesh.full_field_space``) — the SAME space ``L``/``C``/``S``/``F``
    report, so the :class:`~orpheus.numerics.operator.OperatorSum` composition
    guard accepts ``(L + C - S - F - B)`` (Wave O / O.2b R5). ``B`` acts on the
    composite as the ``A_ss`` block (zero bulk; non-zero only on the trace
    block, where the cosine-weighted ``|Ω·n|·w`` partial-current metric lives).
    That block metric is what makes the Hilbert adjoint ``B.H`` the physically
    correct partial-current adjoint — the one channel by which the white-BC
    adjoint becomes available. (Before O.2b R5 ``B`` advertised the bare
    ``sn_mesh.angular_trace`` here, inconsistent with :meth:`apply` already consuming /
    emitting a full :class:`~orpheus.transport.full_field.FullField`.)

    Capabilities
    ------------

    ``apply`` always. ``apply_transpose`` is advertised iff EVERY per-face law
    advertises it — reflective (involution), vacuum, periodic and albedo do;
    **white does NOT** (it is self-adjoint only under the ``|Ω·n|·w`` metric, so
    its adjoint routes through ``B.H`` on the weighted trace space at O.2). The
    intersection rule keeps ``apply_transpose`` honest: it is reachable only when
    every block can honour it.

    Parameters
    ----------
    sn_mesh : SNMesh
        The augmented geometry — carries the per-face boundary laws
        (the face-name-keyed ``bc`` dict) and the unified trace space (same instance the
        composite carrier is bound to; the mesh-identity invariant of
        :class:`~orpheus.sn.operators.streaming.StreamingOperator` applies here too).
    """

    block_role = BlockRole.BOUNDARY

    def __init__(self, sn_mesh: "SNMesh") -> None:
        self.sn_mesh = sn_mesh

    @property
    def _face_laws(self) -> dict[str, LinearOperator]:
        """Map each true boundary face → its per-face realized law.

        Read from ``sn_mesh.bc`` for the faces the trace carries
        (slab ``xmin``/``xmax``; curvilinear ``xmax`` only; 2-D Cartesian
        all four) — the dict and the trace layout share their keys by
        construction (both derived from ``face_labels``, C4 / #220).
        Single source of truth — the laws are the same objects
        the sweep consumes, so ``B`` cannot drift from the realized BCs.
        """
        return {
            face: self.sn_mesh.bc[face]
            for face in self.sn_mesh.angular_trace.layout.faces
        }

    @property
    def is_adjointable(self) -> bool:
        # B = ⊕ per-face laws; the composite adjoint exists iff EVERY face law
        # is adjointable (reflective / vacuum / periodic / albedo are; white is
        # NOT — self-adjoint only under |Ω·n|·w, routed via B.H). The per-face
        # intersection rule, computed recursively like every composite
        # predicate. is_invertible inherits base False (a BC reflection map
        # is not invertible).
        laws = self._face_laws.values()
        return bool(laws) and all(law.is_adjointable for law in laws)

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        # The composite carrier (NOT the bare trace): ``B.apply`` consumes /
        # emits a full FullField (zero bulk + reflected trace), so the
        # advertised space must be the bulk ⊕ trace composite — matching the
        # ``L``/``C``/``S``/``F`` siblings for the OperatorSum composition
        # guard, and carrying the block-diagonal G-adjoint metric ``B.H``
        # reads. Wave O / O.2b R5.
        return self.sn_mesh.full_field_space

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return self.sn_mesh.full_field_space

    def _reflect_trace(
        self, boundary: "AngularBoundaryFlux", method: str,
        faces: "Iterable[str] | None" = None,
        rows: "Mapping[str, np.ndarray] | None" = None,
    ) -> "AngularBoundarySourceSink":
        r"""Core ``A_ss`` action on the trace ALONE — apply each face's law
        (``method`` ∈ {apply, apply_transpose}) to that face's slot, project onto
        the codomain row, and return a boundary-only
        :class:`~orpheus.transport.source_sinks.AngularBoundarySourceSink`.

        ``B`` is the ``A_ss`` block ``V_outflow → V_inflow``: it maps the
        **outflow** trace to the **inflow** trace, so the forward action must be
        non-zero **only on the inflow ordinate slots** of each face. The realized
        per-face law (a :class:`~orpheus.numerics.operator.PermutationOperator`
        for reflective, :class:`AngularAverageOperator` for white, …) is a
        *full-face* operator: e.g. the specular permutation also maps the input
        inflow slots onto the output **outflow** slots (``R·ψ.inflow``). The
        legacy sweep only ever read the inflow slots of ``bc.apply(face)`` so
        that spurious outflow output was harmless — but ``−B`` as a sibling of
        ``L`` in ``(L_full + C − S − F − B)`` reads the WHOLE boundary block, and
        a non-zero outflow emission would corrupt the outflow-definition residual
        ``ψ.outflow − streamed`` (it is supposed to carry no ``B`` term — see the
        block matrix in :mod:`orpheus.sn.operators.boundary`). So the forward
        action is projected onto ``inflow_indices_for_face`` (the consistency
        row): ``B_face = P_inflow ∘ law``.  The Euclidean transpose is therefore
        ``B_faceᵀ = lawᵀ ∘ P_inflow`` — the INPUT is restricted to the forward's
        codomain rows and the full ``lawᵀ`` image is written (B.2d d3).
        Output-projecting ``lawᵀ`` onto the outflow rows instead extracts a
        law's DIAGONAL block: for the vacuum mask (zero-on-inflow ⊕
        identity-on-outflow) that spelled a spurious ``+1`` outflow diagonal
        where the forward is the ZERO map — caught by the A2a grid-reciprocity
        arm on the het-VACUUM sphere.  Off-diagonal permutation laws
        (reflective / albedo) are bit-identical under either spelling, which is
        why every reflective-fixture gate stayed green over the wrong one.
        (The metric-correct Hilbert adjoint ``B.H`` under ``|Ω·n|·w`` is Wave O
        step O.2; this Euclidean ``apply_transpose`` is the un-weighted shadow.)

        This is the **single source of truth** for the boundary reflection: both
        the full-field :meth:`apply` (lifted onto a zero-bulk carrier) and the
        trace-only :meth:`reflect_into_inflow` (the direct-loop inflow seed) route
        through it, so the two cannot drift (Cardinal Rule 2).
        """
        from orpheus.transport.source_sinks import AngularBoundarySourceSink

        # Single mesh source (mesh-identity invariant — see class docstring):
        # the output buffers, the trace selectors, and ``_face_laws`` ALL read
        # ``self.sn_mesh``, so a mismatched input trace cannot desync the
        # projection from the buffer geometry.
        mesh = self.sn_mesh
        trace = mesh.angular_trace
        out_boundary = AngularBoundarySourceSink.zeros_on(mesh)
        # ``faces=None`` reflects every boundary face (the whole-trace ``B``);
        # a face subset restricts the reflection to those faces — the Phase 3
        # Gauss-Seidel octant-group schedule reflects only the just-swept
        # group's reflective OUTGOING faces, leaving the rest of the inflow
        # trace untouched (zero in this returned sink).  ``B`` is block-diagonal
        # over faces, so the subset action is the EXACT restriction (no
        # cross-face coupling is dropped).
        #
        # ``rows`` (#226 step 2) restricts WITHIN a face: per face, only the
        # given ordinate rows of the codomain projection are emitted (a subset
        # of the inflow rows — the schedule-split ``B_lower``/``B_upper``
        # halves of :meth:`split`).  A face absent from ``rows`` emits nothing.
        # Row-granular restriction is exact for the same reason the face
        # restriction is: the projected action writes each target row
        # independently.
        face_laws = self._face_laws
        if faces is not None:
            unknown = set(faces) - set(face_laws)
            if unknown:
                raise ValueError(
                    f"_reflect_trace: face(s) {sorted(unknown)} are not "
                    f"boundary faces of this mesh; available faces: "
                    f"{sorted(face_laws)}."
                )
            face_laws = {f: face_laws[f] for f in faces}
        if rows is not None:
            face_laws = {f: law for f, law in face_laws.items() if f in rows}
        for face, law in face_laws.items():
            face_in = boundary.face_view(face)
            sel = (
                rows[face] if rows is not None
                else trace.inflow_indices_for_face(face)
            )
            if method == "apply":
                full = law.apply(face_in)
                out_boundary.face_view(face)[sel] = full[sel]
            else:
                # ``(P_sel ∘ law)ᵀ = lawᵀ ∘ P_sel``: mask the INPUT to the
                # forward's codomain rows, write the FULL transpose image.
                # The checked bridge licenses the raw verb (spec §39.1) —
                # unreachable-in-practice because :attr:`is_adjointable`
                # gates the composite eagerly, but the per-face raise keeps
                # the refusal loud if a caller bypasses the predicate.
                if not adjointable(law):
                    raise MissingAdjoint(
                        f"SNBoundaryOperator.apply_transpose: face {face!r} "
                        f"law {type(law).__name__} has no Euclidean "
                        f"transpose — reachable only when every face law is "
                        f"adjointable (see is_adjointable)."
                    )
                masked = np.zeros_like(face_in)
                masked[sel] = face_in[sel]
                out_boundary.face_view(face)[...] = law.apply_transpose(masked)
        return out_boundary

    def _apply_faces(
        self, psi: "FullField", method: str,
        rows: "Mapping[str, np.ndarray] | None" = None,
    ) -> "FullField":
        r"""Lift the trace-only :meth:`_reflect_trace` onto the full
        :class:`~orpheus.transport.full_field.FullField` carrier with **zero
        bulk** — ``B_a``, the System-A (trace) boundary block on ``V = V_bulk ⊕
        V_boundary``.  #257 S8a: history-free (the matvec leaf is a base
        arrow ``FullField -> FullField``; the comonad lives on the driver).

        ``B_a`` touches only the trace; the ray-corner boundary is the sibling
        :class:`RadialCharacteristicBoundaryOperator` (``B_b``) — SYSTEM B's
        own boundary block, living at the coupled grid's (B,B) slot (RULING P1
        — a block-composed system's boundary is the direct sum of per-system
        boundary blocks; see the module docstring). B.2d: System B is its own
        composite, so ``B_a`` neither reads nor pads a ray block.
        """
        from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
        from orpheus.transport.full_field import FullField

        mesh = self.sn_mesh
        if psi.interior.mesh is not mesh:
            raise ValueError(
                "SNBoundaryOperator.apply: input field and operator must "
                "share the same SNMesh instance (mesh-identity invariant); "
                f"got field mesh {psi.interior.mesh!r} vs operator mesh {mesh!r}."
            )
        # Role parse at the composite boundary: ``B_a`` reads a FLUX trace
        # (``_reflect_trace`` applies the boundary law to outflow flux), but
        # the ``FullField.boundary`` slot erases the role (the F2-sibling
        # erasure — #289). A source-role
        # trace arriving here is a caller error worth raising loudly.
        trace = psi.boundary
        if not isinstance(trace, AngularBoundaryFlux):
            raise TypeError(
                f"SNBoundaryOperator: the input composite's boundary must "
                f"be an AngularBoundaryFlux trace; got {type(trace).__name__}."
            )
        return FullField(
            interior=_zero_bulk_source(mesh),
            boundary=self._reflect_trace(trace, method, rows=rows),
        )

    def apply(self, psi: "FullField") -> "FullField":
        r"""Forward action ``B_a·ψ`` — per-face boundary law on the trace, zero bulk."""
        return self._apply_faces(psi, "apply")

    def reflect_into_inflow(
        self, boundary: "AngularBoundaryFlux",
        faces: "Iterable[str] | None" = None,
    ) -> "AngularBoundarySourceSink":
        r"""Trace-only forward reflection ``B·ψ.outflow`` projected onto the
        inflow row — the ``A_ss`` action expressed on the boundary trace ALONE.

        Returns a boundary-only
        :class:`~orpheus.transport.source_sinks.AngularBoundarySourceSink` whose
        **inflow** ordinate slots carry the per-face reflected outflow (``R·G``
        for reflective, the angular average for white, zero for vacuum) and whose
        outflow slots are zero. It is :meth:`apply` without the zero-bulk carrier
        — the entry the direct fixed-source SI loop and the final eigenvalue
        reconstruction sweep use to seed ``ψ.inflow = B·ψ.outflow`` on a bare
        boundary buffer, without fabricating a throwaway zero-bulk field just to
        reach the ``A_ss`` block.

        ``faces`` (Phase 3 Gauss-Seidel): ``None`` (default) reflects every
        boundary face — the whole-trace Jacobi reflect used by the fixed-source
        SI loop and the final reconstruction sweep.  A face subset restricts the
        reflection to those faces: the octant-group G-S schedule reflects only
        the just-swept group's reflective OUTGOING faces between octant-group
        sweeps, so a later group reads the fresh reflected inflow (the
        ``(L+C−B_lower)⁻¹`` forward substitution).  ``B`` is block-diagonal over
        faces → the subset is the exact restriction.
        """
        return self._reflect_trace(boundary, "apply", faces=faces)

    def reflect_inflow_inplace(
        self, boundary_flux: "AngularBoundaryFlux",
        faces: "Iterable[str] | None" = None,
    ) -> None:
        r"""In place: overwrite each face's inflow rows with the reflected
        outflow — ``ψ.inflow ← (B_a·ψ)|_{\rm inflow}``, face-restrictable.

        The MUTATING façade over :meth:`reflect_into_inflow` (single source —
        both route through :meth:`_reflect_trace`), matching the sweep
        substrate's reflect signature
        (``Callable[[AngularBoundaryFlux, tuple[str, ...]], None]``): the
        :func:`~orpheus.sn.loss_representation._sweep_scheduled` inter-group
        reflect passes THIS bound method (#226 step 2 — the reified
        ``M = (L+C−B_lower)`` supplies ``boundary.reflect_inflow_inplace``),
        and the whole-trace form (``faces=None``) serves the direct
        fixed-source SI loop + the eigenvalue reconstruction sweep via
        :func:`orpheus.sn.solver._reflect_outflow_into_inflow`.

        Trace-only: the ψ½ ray-corner analogue is
        :meth:`RadialCharacteristicBoundaryOperator.reflect_corner_inplace`
        (``B_b``), the System-B boundary — the reconstruction sweep calls BOTH,
        one per system (RULING P1).
        """
        reflected = self.reflect_into_inflow(boundary_flux, faces=faces)
        trace = self.sn_mesh.angular_trace
        selected = (
            boundary_flux.layout.faces if faces is None else faces
        )
        for face in selected:
            inflow = trace.inflow_indices_for_face(face)
            boundary_flux.face_view(face)[inflow] = (
                reflected.face_view(face)[inflow]
            )

    def split(self, schedule: "SweepSchedule") -> "BoundarySplit":
        r"""Split ``B = B_lower + B_upper`` under ``schedule``'s octant order
        (#226 §17 W2 — the regular matrix splitting of the boundary G-S).

        ``B_lower`` carries exactly the (face, inflow-row) couplings the
        scheduled sweep realizes IN-sweep (rows whose octant group is swept
        strictly after the face's reflect —
        :meth:`~orpheus.sn.loss_representation.sweep_schedule.SweepSchedule.lower_inflow_rows`);
        ``B_upper`` carries the complement (the cyclic back-edges plus every
        row of a never-reflected face — vacuum, white, albedo, periodic),
        lagged by the SI driver as an external gain.  The partition is exact:
        the specular map has no octant-diagonal, and the two row sets are
        complementary within each face's inflow by construction here.

        Returns a named pair so the two construction sites cannot be swapped
        silently: ``M = (L + C) - parts.lower`` and ``gains = (S, parts.upper)``.
        The Jacobi schedule yields an empty lower support (``B_lower = 0``,
        ``B_upper = B``) — the degenerate that recovers the plain lagged-``B``
        iteration.
        """
        lower_rows = schedule.lower_inflow_rows(self.sn_mesh)
        trace = self.sn_mesh.angular_trace
        upper_rows = {
            face: np.setdiff1d(
                trace.inflow_indices_for_face(face),
                lower_rows.get(face, np.empty(0, dtype=np.intp)),
            )
            # The same face set the block-diagonal law iterates (single
            # source — the trace layout and ``bc`` share keys by construction).
            for face in self._face_laws
        }
        return BoundarySplit(
            lower=SNMaskedBoundaryOperator(self, lower_rows, schedule),
            upper=SNMaskedBoundaryOperator(self, upper_rows, schedule),
        )

    def apply_transpose(self, psi: "FullField") -> "FullField":
        r"""Euclidean transpose ``Bᵀ·ψ`` — per-face ``apply_transpose``, zero bulk.

        Reachable only when every per-face law is adjointable (see
        :attr:`is_adjointable`). The white BC has
        no Euclidean transpose; its physically-correct adjoint is ``B.H`` under
        the ``|Ω·n|·w`` trace metric (Wave O step O.2).
        """
        return self._apply_faces(psi, "apply_transpose")


class RadialCharacteristicBoundaryOperator(LinearOperator):
    r"""``B_b`` — System B's (ψ½ ray-corner) boundary law, the ray ``A_ss`` block.

    The boundary operator of the radial-characteristic system (System B of the
    2×2 coupled block operator — the ψ½ ray). A first-class sibling of
    :class:`SNBoundaryOperator` (``B_a``, System A's trace boundary), typed —
    since the B.2b re-type — on **System B's own carrier**: domain = codomain =
    ``sn_mesh.radial_characteristic_composite_space``, acting
    ``RadialCharacteristicComposite → RadialCharacteristicComposite`` (reads
    the boundary member's FLUX corners, emits boundary-member SOURCE corners;
    the interior member is a zero source — "B_b touches the trace/bulk" is now
    UNSPELLABLE, Pattern 4). The system's boundary is still the direct sum of
    per-system boundary blocks (RULING P1) — realized since B.2d as the two
    grid entries ``B_a`` at (A,A) and ``B_b`` at (B,B) of the within-group
    gain grid (:func:`orpheus.sn.coupled_system.build_within_group_system`);
    the B.2b FullField-summed adapter composition is retired. Unconstructable
    on a seedless mesh (System B does not exist there — the ctor guard
    mirrors ``A_BA``'s).

    The action is the ``(R, μ = ∓1)`` corner reflection that closes the ray's
    r = R boundary. The outer face law realizes on the trace carrier as an
    (ordinate ⊗ group) OPERATOR (``B_a``'s per-face law), but at the
    **off-quadrature** μ = ±1 ray it cannot act — so ``B_b`` applies the per-KIND
    specular FACT directly (RULING P1's "one law, two carriers": the
    carrier-indexed realizer that would produce both arms from one law is the
    named seam, built when DSA adds the 3rd carrier).

    Capabilities
    ------------

    ``apply`` always. ``apply_transpose`` iff the outer ray-face law is
    Euclidean-adjointable — reflective (involution) and vacuum are; white /
    albedo / periodic are the loud-deferred set (no ruled off-quadrature corner
    action). Per-leaf :attr:`is_adjointable` (NOT the whole-trace intersection —
    ``B_b`` has exactly one face, the outer radius).

    Adjoint metric (RULING P2)
    --------------------------

    ``B_b`` advertises the **Euclidean** :meth:`apply_transpose` only — no
    per-leaf ``.H``. The ray corner gauge is symmetric (``g₊ = g₋ = V(R)`` — both
    corners at r = R), so ``B_b.H = G_sd⁻¹ B_bᵀ G_sd = B_bᵀ``: the Euclidean
    transpose IS the ``G_sd = V_cell`` Hilbert adjoint. (``.H`` is realized ONCE
    at the composite via ``G⁺·apply_transpose·G``; a Euclidean block adjoint on
    System B is metric-correct here BECAUSE of the symmetric gauge — the
    ``G_sd``-reciprocity gate in ``test_psi_half_coupling.py`` pins the symmetry
    that keeps Mode-12 closed.)

    Parameters
    ----------
    sn_mesh : SNMesh
        The augmented geometry (seed-carrying — 1-D curvilinear). Carries the
        outer-face law ``sn_mesh.bc["xmax"]`` and the ray space
        ``sn_mesh.radial_characteristic_space`` (the mesh-identity invariant of
        :class:`SNBoundaryOperator` applies here too).
    """

    block_role = BlockRole.BOUNDARY
    # B_b is System B's boundary — it acts within the ray system alone
    # (present-zero bulk and trace); campaign step 4a.
    system_role = SystemRole.B

    def __init__(self, sn_mesh: "SNMesh") -> None:
        if sn_mesh.radial_characteristic_composite_space is None:
            raise ValueError(
                "RadialCharacteristicBoundaryOperator: the mesh carries no "
                "ψ½ ray (radial_characteristic_composite_space is None) — a "
                "seedless mesh (Cartesian / cylinder, R12a) has no System B, "
                "hence no ray-corner boundary block. B_b exists only on a "
                "seed-carrying mesh (the sphere)."
            )
        self.sn_mesh = sn_mesh

    @property
    def is_adjointable(self) -> bool:
        # Per-leaf (NOT B_a's whole-trace intersection): B_b has exactly one
        # face, the outer radius. Reflective (involution) + vacuum are
        # Euclidean-adjointable; white / albedo / periodic are the loud-deferred
        # set (:meth:`_reflect_corner` raises — no ruled off-quadrature corner
        # action → no transpose). is_invertible inherits base False. (The old
        # seedless defensive arm is dead under the ctor guard — retired.)
        kind = getattr(self.sn_mesh.bc["xmax"], "kind", None)
        return kind in _RULED_CORNER_KINDS

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        # System B's own member space (B.2b DP1; non-None by the ctor guard).
        # The B.2c CoupledOperator grid type-checks the (B, B) placement
        # against it; the FullField-summed production ``B = B_a + B_b`` rides
        # the transient adapter, which declares ``full_field_space``.
        return self.sn_mesh.radial_characteristic_composite_space

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return self.sn_mesh.radial_characteristic_composite_space

    def _reflect_corner(
        self, seed: "RadialCharacteristicBoundaryField", method: str,
    ) -> "RadialCharacteristicBoundarySourceSink":
        r"""The ``A_ss`` CORNER action on System B's boundary member (R13, 2.5d).

        The (R, μ = ∓1) corner pair closes the ray boundary on a seed-carrying
        mesh: the inward seed leg's r = R inflow is BC data, and for a
        specular-reflective outer face the reflected partner of the outward ray
        μ = +1 is EXACTLY the inward one μ = −1 (its own mirror — an
        off-quadrature ray, so the per-face law OPERATOR cannot act on it; the
        specular fact is applied directly). Forward:
        ``out.corner(level, −1) = ψ½.corner(level, +1)`` per carried level;
        Euclidean transpose: ``out.corner(level, +1) = χ̄.corner(level, −1)``.
        The opposite corners stay ZERO (``B_b`` touches only the inflow row /
        its transpose image — the exact ``_reflect_trace`` projection
        discipline); since the B.2b re-type the input IS the boundary member
        (the cells never enter — structural, not zeroed).

        Law dispatch — on the declared law KIND (the #186 shim's ``kind`` tag,
        the same registry key ``sn_mesh.bc[face] == "reflective"`` comparisons
        read), NOT on the realized operator's composition tree (which is an
        (ordinate ⊗ group) operator over the QUADRATURE rows and structurally
        cannot act on the off-quadrature μ = ±1 ray):

        * ``"vacuum"`` — no corner emission (the block stays all-zero).
        * ``"reflective"`` — the specular corner swap above.
        * anything else (white / albedo / periodic / prescribed) —
          **loud-deferred** (:class:`NotImplementedError`) per the 2.5d
          plan-of-record (e.g. white re-emission at the off-quadrature ray needs
          the ``|Ω·n|``-weighted outflow average for μ = −1, not yet ruled).
        """
        from orpheus.transport.source_sinks import (
            RadialCharacteristicBoundarySourceSink,
        )

        # A seed-carrying mesh is 1-D curvilinear: exactly ONE boundary face
        # (the outer radius renders as ``xmax``) carries the law.
        law = self.sn_mesh.bc["xmax"]
        kind = getattr(law, "kind", None)
        if kind not in _RULED_CORNER_KINDS:  # single source with is_adjointable
            raise NotImplementedError(
                f"RadialCharacteristicBoundaryOperator: the outer-face law kind "
                f"{kind!r} has no ruled corner action yet (white / albedo / "
                f"periodic at the off-quadrature μ = ±1 ray — loud-deferred, "
                f"2.5d plan-of-record)."
            )
        out = RadialCharacteristicBoundarySourceSink.zeros_on(seed.mesh)
        # vacuum ⇒ zero corner emission (the all-zero ``out`` falls through);
        # reflective ⇒ the specular swap (the mirror of μ = +1 is exactly μ = −1).
        if kind == "reflective":
            for level in seed.levels:
                if method == "apply":
                    out.corner(level, -1)[...] = seed.corner(level, +1)
                else:  # apply_transpose — the Euclidean mirror image
                    out.corner(level, +1)[...] = seed.corner(level, -1)
        return out

    def _apply_faces(
        self, ray: "RadialCharacteristicComposite", method: str,
    ) -> "RadialCharacteristicComposite":
        r"""``B_b`` on System B's own carrier: the ray-corner action
        :meth:`_reflect_corner` on the boundary member, a zero-source interior.

        Since the B.2b re-type there is NO bulk/trace padding — the composite
        has no such slots (Pattern 4: "B_b touches the trace" is unspellable).
        The interior member is a zero SOURCE (``B_b`` writes only the corner);
        the production driver consumes this block natively at the (B,B) slot
        of the within-group gain grid (B.2d).
        """
        from orpheus.transport.fields.radial_characteristic_boundary_flux import (
            RadialCharacteristicBoundaryFlux,
        )
        from orpheus.transport.radial_characteristic_composite import (
            RadialCharacteristicComposite,
        )
        from orpheus.transport.source_sinks import (
            RadialCharacteristicInteriorSourceSink,
        )

        mesh = self.sn_mesh
        # The shared System-B block-boundary parse (carrier class +
        # mesh-identity — one parse body across A_BB / A_AB / B_b, B.2c).
        RadialCharacteristicComposite.require_member(
            ray,
            mesh=mesh,
            context=f"RadialCharacteristicBoundaryOperator.{method}",
        )
        # Role parse at the block boundary (the #289-F2 discipline, relocated
        # from the erased FullField slot to the boundary MEMBER — the composite
        # slots are role-erased): ``B_b`` reflects a FLUX corner. A source- /
        # displacement-role member arriving here is a caller error worth
        # raising loudly.
        if not isinstance(ray.boundary, RadialCharacteristicBoundaryFlux):
            raise TypeError(
                f"RadialCharacteristicBoundaryOperator: the input composite's "
                f"boundary must be a RadialCharacteristicBoundaryFlux corner; "
                f"got {type(ray.boundary).__name__}."
            )
        return RadialCharacteristicComposite(
            interior=RadialCharacteristicInteriorSourceSink.zeros_on(mesh),
            boundary=self._reflect_corner(ray.boundary, method),
        )

    def apply(self, ray: "RadialCharacteristicComposite") -> "RadialCharacteristicComposite":
        r"""Forward action ``B_b·ψ½`` — the ray-corner reflection on System B."""
        return self._apply_faces(ray, "apply")

    def apply_transpose(
        self, ray: "RadialCharacteristicComposite",
    ) -> "RadialCharacteristicComposite":
        r"""Euclidean transpose ``B_bᵀ·ψ½`` — the mirror-image corner swap.

        Reachable iff the outer ray-face law is adjointable (see
        :attr:`is_adjointable`). Euclidean = the ``G_sd = V_cell`` Hilbert adjoint
        because the corner gauge is symmetric (RULING P2).
        """
        return self._apply_faces(ray, "apply_transpose")

    def reflect_corner_inplace(
        self, radial_characteristic: "RadialCharacteristicField",
    ) -> None:
        r"""In place: overwrite the ψ½ inflow-corner slots with the law's corner
        action on its OUTFLOW corners — ``ψ½.corner(p, −1) ← (B_b·ψ½).corner(p,
        −1)`` (vacuum ⇒ 0, reflective ⇒ the specular swap).

        The ray-corner analogue of :meth:`SNBoundaryOperator.reflect_inflow_inplace`
        (#282 route (a)): the final eigenvalue reconstruction sweep + the direct
        fixed-source SI loop call BOTH — ``B_a`` for the trace, ``B_b`` for the
        ray — one per system (RULING P1). The input is the UNIFIED walk buffer
        (the walk internals march it through Phase C/4e — the demote ruling);
        it is bridged onto the boundary member so the corner LAW has ONE body
        (:meth:`_reflect_corner`). Non-``None`` by signature since B.2b: a
        seedless B_b is unconstructable, so the caller guards presence
        (see ``_reflect_boundary_inplace``).
        """
        from orpheus.transport.radial_characteristic_composite import (
            RadialCharacteristicComposite,
        )

        boundary_member = RadialCharacteristicComposite.from_unified(
            radial_characteristic,
        ).boundary
        corner_reflected = self._reflect_corner(boundary_member, "apply")
        space = radial_characteristic.space
        for level in space.levels:
            space.corner_view(
                radial_characteristic.values, level, -1,
            )[...] = corner_reflected.corner(level, -1)


# The B.2b transient ``_RayBoundaryFullFieldGain`` adapter RETIRED at B.2d:
# the driver iterate is the CoupledField pair and ``B_b`` sits block-native
# in the (B, B) slot of the within-group gain grid
# (:func:`orpheus.sn.coupled_system.build_within_group_system`) — the
# FullField-summed ``B = B_a + B_b`` composition it bridged is gone.


class SNMaskedBoundaryOperator(LinearOperator["FullField", "FullField"]):
    r"""One half of the schedule split ``B = B_lower + B_upper`` — the
    whole-trace :class:`SNBoundaryOperator` restricted to a per-face set of
    inflow ordinate ROWS (#226 §17 W2).

    The restriction composes a row projection with ``B``'s codomain
    projection: per face, only the given ordinate rows of the reflected
    inflow are emitted; every other slot (and the bulk, as for ``B``) is
    zero.  Which rows belong to which half is SCHEDULE-order semantics
    (:meth:`~orpheus.sn.loss_representation.sweep_schedule.SweepSchedule.lower_inflow_rows`),
    so the instance carries its :attr:`schedule` — the reified
    ``M = (L+C−B_lower)`` reads the walk order off its lower operand rather
    than pairing a foreign schedule with a mismatched mask.  Construct via
    :meth:`SNBoundaryOperator.split` (the named pair keeps lower/upper from
    swapping silently); the exactness of the partition is pinned by the
    W2-split gate.

    A masked half is NOT invertible and does not advertise a transpose
    (``B_lowerᵀ`` masks input rows, not output rows — mint it when the
    adjoint-inverse carve #280 produces a consumer), so it is apply-only
    and the two-axis contract holds by the base defaults.
    """

    block_role = BlockRole.BOUNDARY

    def __init__(
        self,
        inner: "SNBoundaryOperator",
        rows: "Mapping[str, np.ndarray]",
        schedule: "SweepSchedule",
    ) -> None:
        #: The whole-trace boundary law this is a row restriction of.
        self.inner = inner
        #: Per-face inflow ordinate rows this half emits (global ordinate
        #: indices into each face's ``(N, …)`` trace slot).
        self.rows = rows
        #: The octant-order schedule the row split was derived from.
        self.schedule = schedule

    @property
    def sn_mesh(self) -> "SNMesh":
        return self.inner.sn_mesh

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        return self.inner.domain

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return self.inner.codomain

    def apply(self, psi: "FullField") -> "FullField":
        r"""``B_half·ψ`` — the per-face law projected onto :attr:`rows`, zero bulk."""
        return self.inner._apply_faces(psi, "apply", rows=self.rows)

    def reflect_rows_inplace(
        self, boundary_flux: "AngularBoundaryFlux", faces: "Iterable[str]",
    ) -> None:
        r"""In place, ADDITIVE, on :attr:`rows` only:
        ``bf[f][rows] += (B·bf)[f][rows]`` for each given face.

        The inter-group row update of the reified forward substitution
        (#226 §17 W2): solving :math:`M z = y` on a strictly-lower inflow
        row reads :math:`z_{\rm in} = y_{\rm row} + (B z)_{\rm row}` — the
        buffer already holds the seed :math:`y_{\rm row}` (nothing else
        writes a lower row before its face's reflect), so ACCUMULATING the
        fresh reflection completes the inhomogeneous row exactly.  This is
        what makes ``M.inverse()`` exact for arbitrary data on the INFLOW
        rows — i.e. on the source subspace ``{y : y.outflow-rows = 0}``
        (every production rhs; the sweep substrate re-derives the
        outflow-definition rows, a family-wide property shared with
        ``(L+C).solve`` — see the W2 gate module
        ``tests/sn/solve/test_gauss_seidel_reification.py`` and spec §13)
        — not merely on production's zero-lower-inflow-row subspace; and
        restricting to :attr:`rows` leaves the upper (lagged) rows carrying
        the seed the splitting
        :math:`\psi_{k+1} = M^{-1}(q + B_{\rm upper}\psi_k)` says they
        carry — the returned trace IS the splitting's honest iterate.
        (The dissolved resolvent's whole-face OVERWRITE dropped
        :math:`y_{\rm row}` — benign in production where it is zero on a
        reflective face, but O(1)-wrong as an inverse; and it stamped
        fresh values onto rows the iterate defines as lagged.)

        Contrast :meth:`SNBoundaryOperator.reflect_inflow_inplace` — the
        whole-face ASSIGNMENT ``ψ.inflow ← B·ψ.outflow`` between sweeps,
        which is the right semantics for the direct fixed-source SI loop
        and the reconstruction sweep (there the inflow is wholly recomputed
        each sweep, not a solved unknown of a linear row).
        """
        selected = {
            face: self.rows[face]
            for face in faces
            if face in self.rows and np.asarray(self.rows[face]).size
        }
        if not selected:
            return
        reflected = self.inner._reflect_trace(
            boundary_flux, "apply", faces=tuple(selected), rows=selected,
        )
        for face, rows in selected.items():
            boundary_flux.face_view(face)[rows] += (
                reflected.face_view(face)[rows]
            )

    def __repr__(self) -> str:
        n_rows = sum(int(np.asarray(r).size) for r in self.rows.values())
        return (
            f"SNMaskedBoundaryOperator({self.inner!r}, "
            f"rows={n_rows} over {len(self.rows)} faces, "
            f"schedule={self.schedule.kind!r})"
        )


class BoundarySplit(NamedTuple):
    """The named ``B = B_lower + B_upper`` pair from :meth:`SNBoundaryOperator.split`."""

    lower: SNMaskedBoundaryOperator
    upper: SNMaskedBoundaryOperator
