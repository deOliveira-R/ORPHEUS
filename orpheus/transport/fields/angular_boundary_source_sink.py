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

Declaration → snapshot: the inflow terms and this field
=====================================================

There are TWO distinct objects for the inflow :math:`q`, related as **declaration → snapshot** (NOT duplicates):

* :class:`~orpheus.transport.boundary_conditions.InflowTerms` (the transport-layer **declaration**, the affine part of a
  condition): :math:`q` as a weighted sum of uniform and sampled profiles, with no trace or mesh handle. The SAMPLING of it onto
  a half-trace is the method's, through ONE verb, ``RealizerTable.sample_inflow(terms, site, space)``: the table's ``require``
  (a source the method cannot carry is refused at the door), then the registered
  :class:`~orpheus.transport.boundary_conditions.SourceRealizer`'s ``sample``.
* :class:`AngularBoundarySourceSink` (THIS class, the L2 transport **field**): the *eager, whole-boundary, mesh-bound,
  role-typed snapshot*, the materialised :math:`q` as a stored :class:`~orpheus.numerics.field.Field` on
  ``mesh.angular_trace``. The bridge from the first to the second is :meth:`from_binding`: it reads the faces whose bound
  condition carries a nonzero source (``binding.sources``), samples each on that face's :math:`\Gamma_-` through the method's
  realizer table, and packs the flat layout through the same inflow-slot rule :meth:`prescribed_inflow` states.

The geometry-layer generator ``InflowSourceSpec.evaluate(shape)`` and its bridge ``from_specs`` are retired: a second sampling verb
beside ``sample_inflow`` would be two spellings of one quantity, and the old generator was a realization welded onto a
declaration.

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

The **declaration → snapshot bridge** is :meth:`AngularBoundarySourceSink.from_binding`
— materialise the declared per-face :class:`~orpheus.transport.boundary_conditions.InflowTerms` onto the trace by sampling
each at its face's :math:`\Gamma_-` (through the method's realizer table) and packing the flat
layout. It is a DISTINCT route from :meth:`prescribed_inflow`, which serves
the known-per-face-array case; ``from_binding`` serves the declared-terms case, and
delegates its packing to the same body so the "inflow slots written,
everything else zero" rule lives in exactly one place.

It **samples each term at** ``|Γ₋|`` rows (the space the sampler is handed), so the packed result is
the sampler's own output placed on that face's inflow rows and nowhere else — the
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
   source and drove a sweep consuming a typed boundary-source
   field. It was built ahead of that trigger by explicit user ruling, to close
   the affine boundary channel (``.claude/plans/archive/affine_boundary_source_channel.md``):
   the design of record already said prescribed inflow's realization IS
   ``q.boundary``, and the missing bridge was the reason SN had the fence
   (no ``BlockRole.BOUNDARY`` stamp) without the channel. It takes a
   per-face mapping (``binding.sources``), because a boundary's faces carry different
   conditions. At the wave-2 carve it moved from the geometry-layer ``from_specs`` /
   ``from_problem_laws`` pair (``InflowSourceSpec.evaluate``) onto the declared terms.

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
* ``coding-elegance`` Pattern 2 (single source of truth — one sampling verb removes the two-spellings-of-q smell),
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
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping

    from numpy.typing import NDArray

    from orpheus.numerics.spaces.angular_trace_space import AngularTraceSpace
    from orpheus.transport.boundary_conditions import Binding, FaceSite, RealizerTable
    from orpheus.transport.mesh.material_mesh import MaterialMesh


__all__ = ["AngularBoundarySourceSink"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class AngularBoundarySourceSink(AngularBoundaryField, flux=AngularBoundaryFlux):
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

    See the module docstring for the declaration→snapshot relationship to the
    transport-layer :class:`~orpheus.transport.boundary_conditions.InflowTerms` and the
    :meth:`from_binding` bridge that materialises it.
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
        problem: "MaterialMesh",
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

        Note the declaration→snapshot distinction (see the module docstring):
        this builds the snapshot directly from **known per-face arrays**,
        NOT from a declared
        :class:`~orpheus.transport.boundary_conditions.InflowTerms` —
        :meth:`from_binding` is the latter path, and it delegates its
        packing to the same body so this stays the single place that knows the
        "inflow slots written, everything else zero" rule.

        Parameters
        ----------
        problem : MaterialMesh
            The material mesh; ``problem.angular_trace`` supplies the
            trace: its
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
            If the carrier owns no angular trace (not an SN carrier); if a key is not a face of the layout;
            or if a per-face array shape does not match the ``(N, ng)``
            layout slot.
        """
        return cls._pack_inflow(cls._face_space_of(problem), face_values)

    @classmethod
    def _pack_inflow(cls, trace: "AngularTraceSpace", face_values: "Mapping[str, NDArray]") -> "AngularBoundarySourceSink":
        r"""The ONE packing rule: the inflow ordinate rows of each named face are written, everything else is zero.

        Shared by :meth:`prescribed_inflow` (known per-face arrays) and :meth:`from_binding` (declared terms sampled through the
        method's realizer table), so the rule lives in exactly one place.
        """
        bss = cls.zeros(trace)
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

    # ── The declaration → snapshot bridge ────────────────────────────

    @classmethod
    def from_binding(
        cls,
        binding: "Binding",
        trace: "AngularTraceSpace",
        realizers: "RealizerTable",
        face_site: "Callable[[str], FaceSite]",
    ) -> "AngularBoundarySourceSink":
        r"""⭐ The DECLARED boundary conditions' :math:`q` — the user's path.

        Reads ``binding.sources`` (the tagged faces whose bound condition has a nonzero source, read through
        ``BoundaryCondition.source`` of every condition type), samples each face's :class:`~orpheus.transport.boundary_conditions.InflowTerms`
        on that face's inflow half-trace :math:`\Gamma_-` through the method's own realizer table
        (``realizers.sample_inflow(terms, site, trace.inflow_space(face))``, the ONE source verb), and packs the result into the
        flat layout with the same rule :meth:`prescribed_inflow` states: inflow slots written, everything else zero. Faces whose
        condition carries no source contribute nothing: vacuum, specular, white and every deck are :math:`q = 0`, their entire
        content is the linear factor :math:`R`.

        **Why this method exists, and why it is the top of the ladder.** :meth:`prescribed_inflow` takes the source in a form the
        *caller* already holds (known per-face arrays). This one takes it from the form the *problem* holds: a boundary condition
        someone declared. A test or a driver that reaches for :meth:`prescribed_inflow` directly is supplying by hand what the
        declaration should have produced, and is therefore exercising a path no user travels (user ruling, 2026-08-05: *tests
        must route through the machinery that a user would exercise without bypassing code functionality*). The chain is
        ``from_binding -> sample_inflow -> prescribed_inflow``, so the packing rule is stated exactly once.

        Parameters
        ----------
        binding : Binding
            The method-mesh's binding of the declared conditions to its faces.
        trace : AngularTraceSpace
            The problem's unified boundary trace (canonically ``mesh.angular_trace``).
        realizers : RealizerTable
            The METHOD's realizer table (``TransportMethod.boundary_realizers``); transport names no method table, the method
            passes its own. A table that refuses a source on the face (a radial outer face of SN, #653) raises
            :class:`~orpheus.transport.boundary_conditions.BoundaryRefused` here as it did at construction.
        face_site : Callable[[str], FaceSite]
            The method-mesh's ``face_site`` hook: where a face sits, the table's site key.

        Returns
        -------
        AngularBoundarySourceSink
            The declared :math:`q` on ``trace``: all-zero when no face carries a source, which is the overwhelmingly common case.

        Raises
        ------
        ValueError
            If a sampler returns an array that is not exactly the shape of the face's :math:`\Gamma_-` (the sampler contract's one
            invariant), or a source names a face the layout does not carry.
        """
        face_values: "dict[str, NDArray]" = {}
        known = set(trace.layout.faces.keys())
        template = cls.zeros(trace)
        for face, terms in binding.sources.items():
            if face not in known:
                raise ValueError(
                    f"{cls.__name__}.from_binding: the source of {face!r} names no face of the layout; available: {sorted(known)!r}"
                )
            # ⭐ the sampler is handed Γ₋(f) ITSELF, not a shape tuple: the space is the discretization — its rows ARE the inflow
            # ordinates, its ``directions`` are those rows' Ω (so an angular source is spellable), and its ``shape`` is exactly what
            # the sample must return.
            values = np.asarray(realizers.sample_inflow(terms, face_site(face), trace.inflow_space(face)), dtype=float)
            # Scatter onto the full slot so the packing rule stays the only code that knows where inflow rows live.
            full = np.zeros(template.face_view(face).shape, dtype=float)
            full[trace.inflow_indices_for_face(face)] = values
            face_values[face] = full
        return cls._pack_inflow(trace, face_values)
