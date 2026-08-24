r"""Boundary-trace source field :math:`q(\vec r_\Gamma, \hat\Omega, g)`.

The L2 typed wrapper for the **prescribed inflow source** materialised
as a stored field on the boundary trace — the :math:`q` term of the
affine boundary law

.. math::

    \gamma_-\psi \;=\; R\,G\,\gamma_+\psi \;+\; q ,

realised as concrete per-face values packed into the unified
:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` flat layout
(``(layout.total_size,)``). It is the *source* role leaf of
:class:`~orpheus.transport.fields._bases.AngularBoundaryField`, sibling to
:class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux` (flux)
and :class:`~orpheus.transport.residuals.angular_boundary_residual.AngularBoundaryResidual`
(residual).

Recipe → snapshot: relationship to :class:`InflowSourceSpec`
============================================================

There are TWO distinct objects for the inflow :math:`q`, related as
**recipe → snapshot** (NOT duplicates):

* :class:`~orpheus.geometry.boundary._source.InflowSourceSpec` (the
  geometry-layer **generator**, formerly named ``BoundarySource``):
  a lazy, per-face, ``runtime_checkable`` Protocol — ``evaluate(shape)
  -> ndarray`` — that produces inflow values on demand from a bare
  shape, with no trace/mesh handle (deliberately decoupled from the
  trace; the impls :class:`NoSource` / :class:`ConstantInflowSource`
  derive their output from construction-time data only). It is the
  *recipe* a declared affine boundary law carries; the mesh-BC bridge
  :meth:`from_mesh_laws` is what turns it into the snapshot below.
* :class:`AngularBoundarySourceSink` (THIS class, the L2 transport **field**):
  the *eager, whole-boundary, mesh-bound, role-typed snapshot* — the
  materialised :math:`q` as a stored
  :class:`~orpheus.numerics.field.Field` on ``mesh.angular_trace``.

Two renames produced the current name. First (B.3) the geometry
Protocol ``BoundarySource`` was renamed to ``InflowSourceSpec``,
freeing ``BoundarySource`` for this field leaf and completing the
``{Angular,Scalar,Boundary}×{Flux,Source,Residual}`` role grid
(issues #205 / #201). Then (B.5) the whole **source** column was
renamed ``Source`` → ``SourceSink`` so the name carries the signed
rate-density semantics — a *source* (production) and a *sink* (an
operator-loss output such as :math:`L\psi`) are the same quantity with
opposite sign, and the role-leaf type holds both. Hence
``AngularBoundarySourceSink``.

Consumer: every operator's ``.apply`` boundary + the SI/Krylov source (B.5.2)
=============================================================================

B.5.2 (#208) made this leaf **consumer-driven**. Every SN operator's
``.apply`` output boundary IS a :class:`AngularBoundarySourceSink`, because the
operator output is ``Aψ`` — a source/sink, NOT a residual (a residual
arises only from an explicit ``from_balance`` of the output against a
source; the boundary mirrors the bulk's ``AngularSourceSink``).
:class:`StreamingOperator` emits the non-zero boundary block (the bare
post-extraction emission: outflow self-consistency defect + inflow
identity); :class:`SNBoundaryOperator` (``B``) emits ``B·ψ.outflow`` on
the inflow slots; ``Collision`` / ``Scattering`` / ``Fission`` (and the
``F = 0`` ``ZeroOperator`` codomain zero) emit the auto-allocated zero.
The Krylov matvec ``(L+C).apply − (S+B).apply − F.apply`` and the SI rhs
``F.apply + (S+B).apply + q_ext`` both compose as CLOSED
``AngularBoundarySourceSink`` sums (the :class:`TimedFullField` class gate
requires every operator-output boundary to share this class), and
``q_ext.boundary`` is likewise :class:`AngularBoundarySourceSink` — the
prescribed inflow source (zero for vacuum / reflective). The completed
boundary role grid mirrors the bulk exactly::

    .apply        →  AngularBoundarySourceSink   (Aψ — a source/sink)
    .solve        →  AngularBoundaryFlux         (the swept solution trace)
    from_balance  →  AngularBoundaryResidual     (the defect — O.2 honest driver)

A prescribed ``q`` is built directly from known per-face arrays via the
ergonomic :meth:`prescribed_inflow` generator (``{face: (N, ng)}`` →
only the inflow ordinate slots written, the rest zero) — the single
source of truth the non-vacuum MMS and the splitting-invariance probe
consume (it supersedes the zero-allocate + per-face
``face_view(...)[inflow] = …`` slot-fill loop). The lower-level inherited
:meth:`~orpheus.transport.fields._bases.AngularBoundaryField.from_face_arrays`
(every face, full slot incl. outflow) remains for non-inflow uses; the
operator-output zeros allocate on ``mesh.angular_trace`` (CS4b S5).

The **recipe → snapshot bridge** is :meth:`AngularBoundarySourceSink.from_specs`
— materialise per-face :class:`InflowSourceSpec` recipes onto the trace by
evaluating each at its face's :math:`\Gamma_-` shape and packing the flat
layout. It is a DISTINCT route from :meth:`prescribed_inflow`, which serves
the known-per-face-array case; ``from_specs`` serves the lazy-recipe case, and
delegates its packing to ``prescribed_inflow`` so the "inflow slots written,
everything else zero" rule lives in exactly one place.

It **evaluates each spec at** ``(|Γ₋|,) + trailing``, so the packed result is
the spec's own output placed on that face's inflow rows and nowhere else — the
two routes to :math:`q` agree by construction rather than by transcription.

(Until **P3** this paragraph justified the shape by pointing at the retired
``IncomingSourceOperator.apply``, which asked for the same shape from inside
the realized operator. That operator is gone — it was an affine map in a
linear slot, and it delivered :math:`q` into the ``B`` block a second time —
so this channel is now the ONLY route, and the shape is its own contract
rather than a match against something else's. See
:ref:`bc-affine-source-channel`.)

.. note::

   This bridge stood DEFERRED under ``feedback_unify_after_two_instances``
   until 2026-08-05, waiting for a consumer that both declared a non-trivial
   ``InflowSourceSpec`` and drove a sweep consuming a typed boundary-source
   field. It was built ahead of that trigger by explicit user ruling, to close
   the affine boundary channel (``.claude/plans/archive/affine_boundary_source_channel.md``):
   the design of record already said prescribed inflow's realization IS
   ``q.boundary``, and the missing bridge was the reason SN had the fence
   (no ``BlockRole.BOUNDARY`` stamp) without the channel. The deferred sketch
   named it ``from_spec(spec, mesh)``, singular; the realization takes a
   per-face mapping, because a boundary's faces carry different laws.

Units (B.4 — declared as the ``UNITS`` class constant)
======================================================

The boundary trace is **all-flux**: the prescribed inflow ``q`` is a
flux added to :math:`\gamma_-\psi`, so it carries the angular-flux units
:math:`[1/(\mathrm{cm^2 \cdot s \cdot sr})]`
(:data:`~orpheus.numerics.units.ANGULAR_FLUX_UNITS`) — the SAME as
``AngularBoundaryFlux`` and ``AngularBoundaryResidual``. So on the trace, unlike the
bulk, a *source* does NOT pick up the volumetric ``cm⁻³``. eV-free per
the binned-energy convention. Same units, different role — the gate is
class identity. See :mod:`orpheus.numerics.units`.

References
==========

* ``.claude/plans/field_role_typing_view_g.md`` — Phase B step B.3.
* Grand Report v3 §16A.1 (affine boundary form), §5.5 (Field hierarchy).
* ``coding-elegance`` Pattern 2 (single source of truth — the
  ``InflowSourceSpec`` rename removes the same-name-two-concepts smell),
  Pattern 4 (illegal states unrepresentable — cross-role boundary
  arithmetic raises by class identity).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from orpheus.numerics.moment_layout import AVERAGE_MOMENT
from orpheus.numerics.units import ANGULAR_FLUX_UNITS, Unit
from orpheus.transport.fields._bases import AngularBoundaryField

if TYPE_CHECKING:
    from collections.abc import Mapping

    from numpy.typing import NDArray

    from orpheus.geometry.boundary._source import InflowSourceSpec

    from orpheus.sn.mesh.augmented_mesh import SNMesh


__all__ = ["AngularBoundarySourceSink"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class AngularBoundarySourceSink(AngularBoundaryField):
    r"""L2 boundary-trace inflow source — the *source* role leaf of
    :class:`~orpheus.transport.fields._bases.AngularBoundaryField`.

    Parameters
    ----------
    values : NDArray
        Flat backing buffer, shape ``(space.layout.total_size,)``.
    space : AngularTraceSpace
        The unified boundary
        :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`
        (canonically ``mesh.angular_trace``), carrying the per-geometry
        :class:`~orpheus.numerics.face_layout.FaceLayout`.
    mesh : SNMesh
        The SN phase-space carrier (cross-mesh-arithmetic guard).

    Notes
    -----
    A thin role leaf: all storage, validation, algebra, per-face access
    (:meth:`face_view`), and the
    :meth:`~orpheus.transport.fields._bases.AngularBoundaryField.from_face_arrays`
    packer are inherited from :class:`AngularBoundaryField`. The leaf carries
    no source-specific behaviour beyond its class identity — which is
    exactly what Field's Layer-1 gate uses to keep boundary source, flux,
    and residual arithmetic from silently mixing. Note that all three
    boundary leaves share the SAME ``AngularTraceSpace`` (``mesh.angular_trace``); the
    space comparison would pass, so it is the **class** gate (not the
    space gate) that rejects ``AngularBoundarySourceSink ± AngularBoundaryFlux`` /
    ``AngularBoundarySourceSink ± AngularBoundaryResidual``. Same-class arithmetic is
    closed.

    See the module docstring for the recipe→snapshot relationship to the
    geometry-layer :class:`InflowSourceSpec` generator (formerly named
    ``BoundarySource``) and the :meth:`from_specs` bridge that materialises it.
    """

    #: Dimensional identity (View-G, B.4): the boundary is all-flux, so
    #: ``1/(cm²·s·sr)`` — :data:`~orpheus.numerics.units.ANGULAR_FLUX_UNITS`,
    #: shared with ``AngularBoundaryFlux`` / ``AngularBoundaryResidual`` (same units,
    #: different role → class gate). Metadata, not the gate.
    UNITS: ClassVar[Unit] = ANGULAR_FLUX_UNITS

    # ── Ergonomic generator: the prescribed-inflow source ────────────

    @classmethod
    def prescribed_inflow(
        cls,
        mesh: "SNMesh",
        face_values: "Mapping[str, NDArray]",
    ) -> "AngularBoundarySourceSink":
        r"""Build a prescribed-inflow source :math:`q` from per-face values.

        The ergonomic generator for the affine-BC inhomogeneous term
        :math:`q` in :math:`\gamma_-\psi = R\,G\,\gamma_+\psi + q`: for
        each face named in ``face_values``, the **inflow** ordinate slots
        are written from the given per-face slot and **every other slot
        is left zero** (the trailing axes span the group, the transverse
        spatial axis, and — for a moment-resolved LD trace, #251 — the
        ``2^{d-1}``-transverse-moment axis). The outflow slots of a prescribed-inflow
        source are physically meaningless — the sweep determines outflow,
        not the source — so they are *unrepresentable* here by
        construction (``coding-elegance`` Pattern 4). Faces absent from
        ``face_values`` are vacuum (all-zero).

        This is the prescribed-inflow specialisation of the general
        :meth:`~orpheus.transport.fields._bases.AngularBoundaryField.from_face_arrays`
        (which requires EVERY face and writes the FULL per-face slot,
        outflow included): here only the faces that carry incoming flux
        need be named, and only their inflow ordinates are honoured. It
        supersedes the zero-allocate + per-face
        ``face_view(...)[inflow] = …`` slot-fill loop that every
        prescribed-inflow consumer (the non-vacuum MMS, the splitting-
        invariance probe) previously hand-rolled — the single source of
        truth for materialising a prescribed inflow onto the trace
        (Cardinal Rule 2).

        Note the recipe→snapshot distinction (see the module docstring):
        this builds the snapshot directly from **known per-face arrays**,
        NOT from a lazy
        :class:`~orpheus.geometry.boundary._source.InflowSourceSpec`
        recipe — :meth:`from_specs` is the latter path, and it delegates its
        packing here so this method stays the single place that knows the
        "inflow slots written, everything else zero" rule.

        Parameters
        ----------
        mesh : SNMesh
            The SN phase-space carrier; ``mesh.angular_trace`` supplies the
            :class:`~orpheus.numerics.face_layout.FaceLayout` and the
            per-face inflow ordinate index sets
            (:meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.inflow_indices_for_face`).
        face_values : Mapping[str, NDArray]
            ``{face_name: values}`` where ``values`` is the full per-face
            slot, shape ``(N, ng, *face_shape[, 2^{d-1}])`` over all
            ordinates (the trailing transverse-moment axis present only
            for a moment-resolved LD trace, #251). Only the inflow
            ordinate rows are read; the remainder is ignored.

        Returns
        -------
        AngularBoundarySourceSink
            The materialised :math:`q` on ``mesh.angular_trace`` — inflow slots
            of the named faces set, everything else zero.

        Raises
        ------
        ValueError
            If ``mesh.angular_trace is None`` (a trace-less 2-D cylindrical mesh,
            which has no SN sweep); if a key is not a face of the layout;
            or if a per-face array shape does not match the ``(N, ng)``
            layout slot.
        """
        trace = mesh.angular_trace
        if trace is None:
            raise ValueError(
                f"{cls.__name__}.prescribed_inflow: mesh has no AngularTraceSpace "
                f"(mesh.angular_trace is None — trace-less 2-D cylindrical). A "
                f"boundary source cannot be built without a trace."
            )
        bss = cls.zeros(cls._face_space_of(mesh))
        known = set(trace.layout.faces.keys())
        for face, values in face_values.items():
            if face not in known:
                raise ValueError(
                    f"{cls.__name__}.prescribed_inflow: {face!r} is not a "
                    f"face of the layout; available: {sorted(known)!r}"
                )
            view = bss.face_view(face)  # full per-face slot, mem-shared with bss
            arr = np.asarray(values, dtype=float)
            inflow = trace.inflow_indices_for_face(face)
            # A moment-resolved LD trace slot is (N, ng, *transverse, 2^{d-1})
            # (#251).  The caller may supply EITHER the full moment slot OR a
            # SCALAR (N, ng, *transverse) array (the existing scalar prescribed
            # inflow — every face-AVERAGE caller).  Discriminate by shape:
            #   * full slot  → write the whole slot (DD/Step's scalar slot is the
            #     full slot, so this is the byte-identical path);
            #   * scalar onto a moment slot → seed the AVERAGE moment (slot 0),
            #     the transverse slope moments stay zero (the scalar default is
            #     blind to the along-face variation — the Leg-B asymmetry).
            # ``[inflow]`` selects the inflow ordinate ROWS (axis 0) and spans
            # all trailing axes (a trailing ``, :`` would assume exactly one).
            if arr.shape == view.shape:
                view[inflow] = arr[inflow]
            elif view.ndim > arr.ndim and arr.shape == view.shape[:-1]:
                view[inflow, ..., AVERAGE_MOMENT] = arr[inflow]
            else:
                raise ValueError(
                    f"{cls.__name__}.prescribed_inflow: face {face!r} values "
                    f"shape {arr.shape!r} does not match the layout slot "
                    f"shape {view.shape!r} (expected the full per-face slot, or "
                    f"a scalar {view.shape[:-1]!r} onto a moment-resolved slot)."
                )
        return bss

    # ── The recipe → snapshot bridge ─────────────────────────────────

    @classmethod
    def from_specs(
        cls,
        mesh: "SNMesh",
        face_specs: "Mapping[str, InflowSourceSpec]",
    ) -> "AngularBoundarySourceSink":
        r"""Materialise lazy :class:`InflowSourceSpec` **recipes** onto the trace.

        The recipe→snapshot bridge the module docstring describes: each face's
        generator is evaluated at that face's :math:`\Gamma_-` shape and the
        result packed into the flat layout, turning a per-face lazy rule into
        one eager, mesh-bound, role-typed :math:`q`.

        ⭐ **Evaluated at** ``(|Γ₋|,) + trailing``, which is this channel's own
        contract. Until **P3** the shape was justified by matching the retired
        ``IncomingSourceOperator.apply``, and that match is exactly what let the
        carve retire it: the operator's whole body was
        ``source.evaluate((n_inflow,) + trailing)``, so nothing was lost by
        deleting the wrapper. Since P3 this is the ONLY route to :math:`q`
        (:ref:`bc-affine-source-channel`).

        Packing is delegated to :meth:`prescribed_inflow` rather than repeated
        here: that method is the single source of truth for "inflow slots
        written, everything else zero", including the moment-resolved LD case
        (#251). This method's only job is the recipe evaluation.

        .. note::

           The deferred sketch in the module docstring named this
           ``from_spec(spec, mesh)`` — ONE spec for the whole boundary. The
           realization takes a **per-face mapping**, because a mesh's faces
           carry different laws: a single spec applied boundary-wide would
           impose an inflow on vacuum and reflective faces too. Argument order
           follows :meth:`prescribed_inflow` (mesh first, mapping second) so
           the two read alike.

        Parameters
        ----------
        mesh : SNMesh
            The SN phase-space carrier; ``mesh.angular_trace`` supplies the
            layout and the per-face inflow ordinate index sets.
        face_specs : Mapping[str, InflowSourceSpec]
            ``{face_name: spec}``. Faces absent from the mapping are vacuum
            (all-zero), matching :meth:`prescribed_inflow`.

        Returns
        -------
        AngularBoundarySourceSink
            The materialised :math:`q` — inflow slots of the named faces set
            from their recipes, everything else zero.

        Raises
        ------
        ValueError
            If ``mesh.angular_trace is None``; if a key is not a face of the
            layout; or if a spec returns an array of the wrong shape (the
            :class:`InflowSourceSpec` contract's one invariant).
        """
        trace = mesh.angular_trace
        if trace is None:
            raise ValueError(
                f"{cls.__name__}.from_specs: mesh has no AngularTraceSpace "
                f"(mesh.angular_trace is None — trace-less 2-D cylindrical). A "
                f"boundary source cannot be built without a trace."
            )
        known = set(trace.layout.faces.keys())
        template = cls.zeros(cls._face_space_of(mesh))
        face_values: "dict[str, NDArray]" = {}
        for face, spec in face_specs.items():
            if face not in known:
                raise ValueError(
                    f"{cls.__name__}.from_specs: {face!r} is not a face of the "
                    f"layout; available: {sorted(known)!r}"
                )
            slot_shape = template.face_view(face).shape
            inflow = trace.inflow_indices_for_face(face)
            # ⭐ P6: hand the spec Γ₋(f) ITSELF, not a shape tuple. The space is
            # the discretization — its rows ARE ``inflow``, its ``directions``
            # are those rows' Ω (so an angular source is spellable), and its
            # ``shape`` is what a hand-rolled
            # ``(|Γ₋|,) + slot_shape[1:]`` used to recompute here. That
            # duplicate is gone: one object now answers "which rows", "which
            # directions" and "what shape", so they cannot disagree.
            space = trace.inflow_space(face)
            wanted = tuple(int(s) for s in space.shape)
            values = np.asarray(spec.evaluate(space), dtype=float)
            if values.shape != wanted:
                raise ValueError(
                    f"{cls.__name__}.from_specs: face {face!r} spec "
                    f"{spec!r} returned shape {values.shape!r}, expected "
                    f"{wanted!r} (= the shape of {space.name}). The "
                    f"InflowSourceSpec contract's ONE invariant is that the "
                    f"returned array has exactly the space's shape."
                )
            # Scatter onto the full slot so `prescribed_inflow` — which reads
            # `[inflow]` from a full-slot array — stays the only code that
            # knows the packing. The zeros in the outflow rows are never read.
            full = np.zeros(slot_shape, dtype=float)
            full[inflow] = values
            face_values[face] = full
        return cls.prescribed_inflow(mesh, face_values)

    @classmethod
    def from_mesh_laws(cls, mesh: "SNMesh") -> "AngularBoundarySourceSink":
        r"""⭐ The DECLARED boundary conditions' :math:`q` — the user's path.

        Reads each face's realized law (``mesh.bc[face].law``) and materialises
        the inflow source of every
        :class:`~orpheus.geometry.boundary.PrescribedInflow` among them. Faces
        carrying any other law contribute nothing: vacuum, reflective, white,
        albedo and periodic are all :math:`q = 0` — their entire content is the
        linear factor :math:`L`.

        **Why this method exists, and why it is the top of the ladder.** The
        four constructors below it each take the source in a form the *caller*
        already holds — zeros, full slots, known inflow arrays, lazy recipes.
        This one takes it from the form the *problem* holds: a boundary
        condition someone declared. A test or a driver that reaches for
        :meth:`prescribed_inflow` directly is supplying by hand what the
        declaration should have produced, and is therefore exercising a path no
        user travels (user ruling, 2026-08-05:

            *Tests must route through the machinery that a user would exercise
            without bypassing code functionality. Or else it's not testing the
            path the users go through.*

        ) — so a manufactured-solution driver with a non-vacuum inflow declares
        a ``PrescribedInflow`` and comes through here.

        The chain is ``from_mesh_laws → from_specs → prescribed_inflow``, each
        step delegating down, so the packing rule is stated exactly once.

        .. note::

           A non-trivial inflow is reachable ONLY by constructing the law
           object; ``BC.params`` is ``dict[str, float]``, so a ``BC(...)`` tag
           can never carry a manufactured solution restricted to a face. That
           is structural, not a gap — and it is separate from **#189**, which
           is about registering ``prescribed_inflow`` as a tag *kind* for the
           constant case.

        Returns
        -------
        AngularBoundarySourceSink
            The declared :math:`q` — all-zero when no face declares a
            prescribed inflow, which is the overwhelmingly common case
            (a zero trace field on ``mesh.angular_trace``).

        Raises
        ------
        ValueError
            If ``mesh.angular_trace is None`` (trace-less 2-D cylindrical).
        """
        from orpheus.geometry.boundary import PrescribedInflow

        trace = mesh.angular_trace
        if trace is None:
            raise ValueError(
                f"{cls.__name__}.from_mesh_laws: mesh has no AngularTraceSpace "
                f"(mesh.angular_trace is None — trace-less 2-D cylindrical). A "
                f"boundary source cannot be built without a trace."
            )
        specs = {}
        for face in trace.layout.faces:
            law = getattr(mesh.bc[face], "law", None)
            if isinstance(law, PrescribedInflow):
                specs[face] = law.source
        return cls.from_specs(mesh, specs)
