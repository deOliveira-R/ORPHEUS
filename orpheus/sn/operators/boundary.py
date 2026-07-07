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

:class:`SNBoundaryOperator` is the ``A_ss`` leaf. As a sibling ``−B`` of ``L``
in the :class:`~orpheus.numerics.operator.OperatorSum` algebra it supplies the
reflective coupling that ``L`` previously absorbed inside its sweep (the
``inflow = bc.apply(outflow)`` re-apply); the outer Krylov / SI loop then drives
the boundary **consistency residual** ``ψ.inflow − B·ψ.outflow − q.inflow → 0``.

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
)

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping

    from orpheus.numerics.space import FunctionSpace
    from orpheus.sn.loss_representation.sweep_schedule import SweepSchedule
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.transport.fields._bases import RadialCharacteristicField
    from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
    from orpheus.transport.full_field import FullField
    from orpheus.transport.source_sinks import (
        AngularBoundarySourceSink,
        RadialCharacteristicSourceSink,
    )


__all__ = ["BoundarySplit", "SNBoundaryOperator", "SNMaskedBoundaryOperator"]


class SNBoundaryOperator(LinearOperator):
    r"""Whole-trace boundary law ``B`` — the ``A_ss`` block of the SN algebra.

    Block-diagonal over the mesh's true boundary faces: ``B.apply(ψ)`` returns a
    :class:`~orpheus.transport.full_field.FullField` with **zero bulk**
    and, on each face, ``bc[<face>].apply(ψ.boundary.face_view(<face>))`` — the
    per-face realized boundary law applied to that face's trace slot. It composes
    as ``−B`` in ``(L_full + C − S − F − B)`` (acting on the same
    :class:`~orpheus.transport.full_field.FullField` carrier as ``L``/``C``/``S``/``F``).

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
        row); the Euclidean transpose ``Bᵀ: V_inflow → V_outflow`` is projected
        onto ``outflow_indices_for_face`` accordingly. (The metric-correct Hilbert
        adjoint ``B.H`` under ``|Ω·n|·w`` is Wave O step O.2; this Euclidean
        ``apply_transpose`` is the un-weighted shadow, not yet a live consumer.)

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
            full = getattr(law, method)(face_in)
            if rows is not None:
                target = rows[face]
            else:
                target = (
                    trace.inflow_indices_for_face(face)
                    if method == "apply"
                    else trace.outflow_indices_for_face(face)
                )
            out_boundary.face_view(face)[target] = full[target]
        return out_boundary

    def _apply_faces(
        self, psi: "FullField", method: str,
        rows: "Mapping[str, np.ndarray] | None" = None,
    ) -> "FullField":
        r"""Lift the trace-only :meth:`_reflect_trace` onto the full
        :class:`~orpheus.transport.full_field.FullField` carrier with **zero
        bulk** — the ``A_ss`` block as an operator on ``V = V_bulk ⊕
        V_boundary``.  #257 S8a: history-free (the matvec leaf is a base arrow
        ``FullField -> FullField``; the comonad lives on the driver).
        """
        from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
        from orpheus.transport.full_field import FullField
        from orpheus.transport.source_sinks import AngularSourceSink

        mesh = self.sn_mesh
        if psi.bulk.mesh is not mesh:
            raise ValueError(
                "SNBoundaryOperator.apply: input field and operator must "
                "share the same SNMesh instance (mesh-identity invariant); "
                f"got field mesh {psi.bulk.mesh!r} vs operator mesh {mesh!r}."
            )
        # Role parse at the composite boundary: ``B`` reads a FLUX trace
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
            # Zero bulk source, sized from the MESH — not ``zeros_like(psi.bulk)``
            # — so the carrier is correct whatever representation the input bulk
            # carries (full-angular :class:`AngularFlux` OR the Phase-5a windowed
            # :class:`HarmonicMomentFlux`).  ``B`` is the boundary block ``A_ss``:
            # it reads the trace and emits zero bulk regardless.  The zero bulk
            # carries the input's spatial-moment width (#240 D5b-S3) so it
            # composes element-wise with the moment-carrying ``(L+C)ψ`` in the
            # SI / Krylov ``(L+C − S − B)ψ`` matvec.  DD/Step → no factor.
            bulk=AngularSourceSink.zeros_on(
                mesh, spatial_moments=mesh.scheme.spatial_basis_per_axis,
            ),
            boundary=self._reflect_trace(trace, method, rows=rows),
            radial_characteristic=self._reflect_radial_characteristic(
                psi.radial_characteristic, method,
            ),
        )

    def _reflect_radial_characteristic(
        self, seed: "RadialCharacteristicField | None", method: str,
    ) -> "RadialCharacteristicSourceSink | None":
        r"""The ``A_ss`` CORNER arm on the starting-direction block (R13, 2.5d).

        The (R, μ = ∓1) corner pair completes the boundary block on a
        seed-carrying mesh: the inward seed leg's r = R inflow is BC data,
        and for a specular-reflective outer face the reflected partner of
        the outward starting direction μ = +1 is EXACTLY the inward one
        μ = −1 (its own mirror — an off-quadrature ray, so the per-face
        law OPERATOR cannot act on it; the specular fact is applied
        directly). Forward: ``out.corner(level, −1) = ψ½.corner(level, +1)``
        per carried level; Euclidean transpose: ``out.corner(level, +1) =
        χ̄.corner(level, −1)``. The cells legs and the opposite corners
        stay ZERO (``B`` touches only the inflow row / its transpose
        image — the exact `_reflect_trace` projection discipline).

        Law dispatch — on the declared law KIND (the #186 shim's ``kind``
        tag, the same registry key ``sn_mesh.bc[face] == "reflective"``
        comparisons read), NOT on the realized operator's composition
        tree: the realized law is an (ordinate ⊗ group) operator over the
        QUADRATURE rows and structurally cannot act on the off-quadrature
        μ = ±1 ray — the corner action is a per-KIND angular fact:

        * ``"vacuum"`` — no corner emission (the block stays all-zero).
        * ``"reflective"`` — the specular corner swap above (the mirror
          of μ = +1 is exactly μ = −1).
        * anything else (white / albedo / periodic / prescribed) —
          **loud-deferred** (`NotImplementedError`) per the 2.5d
          plan-of-record: e.g. the white re-emission at the
          off-quadrature ray needs the |Ω·n|-weighted outflow average
          evaluated for μ = −1, a design decision not yet ruled.

        ``None`` passes through (a seedless composite — Cartesian /
        cylinder, R12a). NOTE the metric-invisibility honest scope
        (§16.A A4): G3 cannot see this arm; the Euclidean Mᵀ chain at
        2.5b and the §16.C round-trip are its catchers.
        """
        if seed is None:
            return None
        from orpheus.transport.source_sinks import RadialCharacteristicSourceSink

        # A seed-carrying mesh is 1-D curvilinear: exactly ONE boundary
        # face (the outer radius renders as ``xmax``) carries the law.
        law = self._face_laws["xmax"]
        kind = getattr(law, "kind", None)
        out = np.zeros_like(seed.values)
        space = seed.space
        if kind == "vacuum":
            pass  # zero corner emission.
        elif kind == "reflective":
            for level in space.levels:
                if method == "apply":
                    space.corner_view(out, level, -1)[:] = (
                        space.corner_view(seed.values, level, +1)
                    )
                else:  # apply_transpose — the Euclidean mirror image
                    space.corner_view(out, level, +1)[:] = (
                        space.corner_view(seed.values, level, -1)
                    )
        else:
            raise NotImplementedError(
                f"SNBoundaryOperator starting-direction corner arm: the "
                f"outer-face law kind {kind!r} has no ruled corner action "
                f"yet (white / albedo / periodic at the off-quadrature "
                f"μ = ±1 ray — loud-deferred, 2.5d plan-of-record)."
            )
        return RadialCharacteristicSourceSink(
            values=out, space=space, mesh=seed.mesh,  # type: ignore[attr-defined]
        )

    def apply(self, psi: "FullField") -> "FullField":
        r"""Forward action ``B·ψ`` — per-face boundary law on the trace, zero bulk."""
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
        *,
        radial_characteristic: "RadialCharacteristicField | None" = None,
    ) -> None:
        r"""In place: overwrite each face's inflow rows with the reflected
        outflow — ``ψ.inflow ← (B·ψ)|_{\rm inflow}``, face-restrictable.

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

        ``radial_characteristic`` (#282 route (a)): a ψ½ carrier whose
        inflow-corner slots are overwritten with the law's corner action
        on its OUTFLOW corners — ``ψ½.corner(p, −1) ← (B·ψ½).corner(p,
        −1)`` through the SAME :meth:`_reflect_radial_characteristic` arm
        the matvec uses (vacuum ⇒ 0, reflective ⇒ the specular corner
        swap).  The seed analogue of the trace overwrite above, for the
        reconstruction sweep's given-data corner seeding.
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
        if radial_characteristic is not None:
            corner_reflected = self._reflect_radial_characteristic(
                radial_characteristic, "apply",
            )
            assert corner_reflected is not None  # seed given ⇒ arm emits
            space = radial_characteristic.space
            for level in space.levels:
                space.corner_view(
                    radial_characteristic.values, level, -1,
                )[...] = space.corner_view(corner_reflected.values, level, -1)

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
