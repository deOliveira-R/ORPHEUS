r"""DiffusionMesh — the diffusion phase space (mesh + materials + realized trace).

:class:`DiffusionMesh` is the diffusion **method-mesh**: a
:class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` (the
method-agnostic mesh + materials DATA) augmented with the diffusion
method's BEHAVIOR —

* the **scalar boundary-trace space** :attr:`scalar_trace` carrying the
  per-face ``(J⁺, J⁻)`` partial-current pairs under the face-AREA
  metric (#290 P2; crosswalk contract);
* the **composite carrier** :attr:`full_field_space`
  :math:`V_{\rm bulk} \oplus V_{\rm trace}` that the operator family
  :math:`A = L + C - S - B` acts on (#290 P4);
* the per-face boundary conditions **realized at construction** —
  :attr:`bc` maps each boundary face to the realized albedo operator
  :math:`J^- = \mathcal{A}\,J^+` that :data:`~orpheus.diffusion.boundary_realizer.DIFFUSION_REALIZERS`
  produces from the :attr:`binding` (the declared ``boundary_conditions``
  joined to the carrier's tagged faces, and its decks). Vacuum is the empty
  network, :math:`\mathcal{A} = 0` (tier 1; the Marshak closure is the method's);
  the zero-flux Dirichlet idealization is its own method-only condition,
  :class:`~orpheus.diffusion.zero_flux.ZeroFluxCondition`, :math:`\mathcal{A} = -1`.

It is the structural sibling of
:class:`~orpheus.sn.problem.SNProblem` (mesh + materials +
quadrature + sweep machinery): ONE method-agnostic data carrier, one
method layer per transport method. The two method-meshes are the
witnesses of the :class:`~orpheus.transport.method.TransportMethod`
Protocol (minted #290 P7b) — conformance is structural, checked where
the shared :func:`~orpheus.transport.method.realize_boundary`
body is applied below; the Protocol class is never imported here.

The data/behavior axis (the P7a restoration)
============================================

``material_mesh.py``'s module docstring declares the axis: a
:class:`MaterialMesh` is method-agnostic transport *data*; the method
layer (quadrature, stencil, **boundary trace**, closures) is
*behavior*. #290 P4, having no diffusion method-mesh to put them on,
parked ``scalar_trace`` + the scalar composite on :class:`MaterialMesh`
itself — same file, opposite claims. P7a resolves the contradiction by
minting THIS class and reclaiming both members: **MaterialMesh does not
know what a trace is** (:class:`SNProblem` knows angular traces;
:class:`DiffusionMesh` knows scalar traces). The type-vs-property rule
is satisfied honestly: two non-isomorphic augmentations of the same
data carrier (quadrature machinery vs trace + realized BCs), real
morphisms applied to each (promotion, BC realization, operator/field
consumption), and the illegal state "a diffusion phase space with
unresolved or unrealizable boundary conditions" is unrepresentable —
every constructible :class:`DiffusionMesh` carries realized operators for
exactly its boundary faces.

Construction semantics
======================

Construction mirrors ``SNProblem._init_core`` order: the method-agnostic
data block first (:meth:`MaterialMesh._init_data` — bit-identical), then
the method layer's admission gates (1-D only today — slab / cylinder /
sphere through the mesh's own areas + volumes), then
the trace, then the binding (:func:`~orpheus.transport.boundary_conditions.bind`
of the declared conditions to this carrier's faces) and the realized ``bc`` dict.
Every gate fires AT CONSTRUCTION, not lazily inside a solver step — operators
built on a bad phase space are action-at-a-distance otherwise: an unbound tag, an
unused condition, a condition outside the sub-Markov cone and a declaration the
diffusion table does not admit (a source, a translation deck) each raise here.

The promotion classmethod :meth:`from_material_mesh` takes the
``boundary_conditions`` keyword only (the tags are part of the data carrier, in
the axes' roles; the CONDITIONS they name are the method's input, as in
``SNProblem.from_material_mesh(mm, quadrature, scheme=..., boundary_conditions=...)``), defaulting to the carrier's own
``binding.declared`` when it has one.
It has no quadrature to inject, so the promotion signatures still do not unify
and should not.

Because :class:`SNProblem` *is a* :class:`MaterialMesh`, an SN phase space
promotes directly — ``DiffusionMesh.from_material_mesh(problem)``
builds the diffusion phase space over the SAME axes / materials /
mat_map, realizing the diffusion reading of the same physical BC tags.
That is exactly the DSA construction path (#2): :math:`A_{\rm diff}`
assembles on the promoted mesh, sharing the geometry with the SN sweep
it accelerates.

Layer (``tests/gates/test_layer_imports.py``): L3 ``diffusion``. Imports
``numerics`` (spaces), ``transport`` (the MaterialMesh base + axis helpers,
the conditions package and the shared
:func:`~orpheus.transport.method.realize_boundary` body),
and its ``diffusion`` siblings (realizer table + method space).
"""

from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from orpheus.diffusion.boundary_realizer import DIFFUSION_REALIZERS
from orpheus.diffusion.method_space import DiffusionMethodSpace
from orpheus.mesh import Mesh1D, Mesh2D
from orpheus.numerics.spaces.scalar_trace_space import ScalarTraceSpace
from orpheus.transport.boundary_conditions import BoundaryConditions, bind
from orpheus.transport.method import MethodMesh, TransportMethod, realize_boundary
from orpheus.mesh import (
    axes_from_legacy_mesh,
    face_labels as _face_labels,
    face_shape as _face_shape,
)
from orpheus.transport.mesh.material_mesh import MaterialMesh

if TYPE_CHECKING:
    from collections.abc import Mapping

    from orpheus.data.materials import Materials
    from orpheus.data.macro_xs.mixture import Mixture
    from orpheus.geometry.deck import DeckTransformation
    from orpheus.numerics.operator import LinearOperator
    from orpheus.numerics.spaces.full_field_space import FullFieldSpace
    from orpheus.mesh import Axis1D
    from orpheus.transport.boundary_conditions import Binding, BoundaryCondition, RealizerTable


__all__ = ["DiffusionMesh"]


class DiffusionMesh(MethodMesh):
    r"""Diffusion method-mesh: MaterialMesh + scalar trace + realized BCs.

    See the module docstring for the data/behavior axis, the
    construction semantics, and the promotion / DSA paths.

    Parameters
    ----------
    mesh : Mesh1D
        Base geometry + zoning + per-face roles (a deck or an integer
        tag; 1-D; slab / cylinder / sphere through its own areas +
        volumes — a multi-D mesh is refused at construction). The same inbound
        surface as :class:`MaterialMesh`; a ``Mesh2D`` is accepted by
        the signature for a uniform refusal message, never admitted.
    materials : dict mapping material id to Mixture
        Macroscopic cross sections keyed by the mesh's ``mat_ids``
        values; the diffusion coefficient reads through the #290 P1
        seam (``Mixture.diffusion_coefficient``).
    boundary_conditions : BoundaryConditions or mapping of int tag to BoundaryCondition
        Keyword-only, default empty. The condition each tag of the mesh's
        faces stands for (the twin of ``materials``); an unbound tag, an
        unused condition and a condition the carrier's colouring symmetry
        does not fix are refused by :func:`~orpheus.transport.boundary_conditions.bind`.

    Attributes
    ----------
    scalar_trace : ScalarTraceSpace
        The per-face ``(J⁺, J⁻)`` boundary-trace space (face-AREA
        metric; a curvilinear pole is not a face and never appears).
    binding : Binding
        The declared conditions joined to this mesh's faces
        (``SNProblem.binding`` parity).
    boundary_realizers : RealizerTable
        Class attribute: :data:`~orpheus.diffusion.boundary_realizer.DIFFUSION_REALIZERS`,
        the one door that admits and realizes a declaration on a face.
    bc : dict[str, LinearOperator]
        Per-face REALIZED boundary operators (``SNProblem.bc`` parity), keyed
        by face name — the albedo operators :math:`J^- = \mathcal{A}
        J^+`. Built from the SAME ``face_labels`` inventory as the
        trace, so operator coverage ≡ face coverage by construction.
    """

    #: The method's realizer table (the ``SNProblem.boundary_realizers`` precedent): the one door that answers "is this declaration
    #: admitted on this face" and "realize it". It replaces the string-tag registry ``BOUNDARY_OPERATOR_REGISTRY``.
    boundary_realizers: ClassVar["RealizerTable"] = DIFFUSION_REALIZERS

    def __init__(
        self,
        mesh: Mesh1D | Mesh2D,
        materials: "Materials | Mapping[int, Mixture]",
        sigma_t_cell: np.ndarray | None = None,
        *,
        boundary_conditions: "BoundaryConditions | Mapping[int, BoundaryCondition]" = BoundaryConditions.empty(),
    ) -> None:
        # The legacy inbound surface (MaterialMesh parity): convert the
        # mesh declaration to the canonical axis tuple ONCE, extract the
        # material assignment, and run the one construction body the
        # promotion path shares.
        self._init_core(
            axes=axes_from_legacy_mesh(mesh),
            mesh=mesh,
            mat_map=mesh.mat_ids if isinstance(mesh, Mesh1D) else mesh.mat_map,
            materials=materials,
            sigma_t_cell=sigma_t_cell,
            boundary_conditions=boundary_conditions,
        )

    def _init_core(
        self,
        *,
        axes: "tuple[Axis1D, ...]",
        mesh: Mesh1D | Mesh2D | None,
        mat_map: np.ndarray | None,
        materials: "Materials | Mapping[int, Mixture]",
        sigma_t_cell: np.ndarray | None = None,
        boundary_conditions: "BoundaryConditions | Mapping[int, BoundaryCondition]",
    ) -> None:
        r"""The ONE construction body both surfaces funnel into
        (``SNProblem._init_core`` parity).

        Data block first, then the diffusion method layer: admission
        gates → trace → binding → realized boundary operators.
        """
        # ── Method-agnostic DATA block → MaterialMesh base ──
        MaterialMesh._init_data(
            self, axes=axes, mesh=mesh, mat_map=mat_map, materials=materials,
            sigma_t_cell=sigma_t_cell,
        )

        # ── Diffusion method layer (BEHAVIOR atop the data) ──
        # Admission gates — the phase space's own contract, fired at
        # construction (P7a moved the 1-D refusal here from the P5
        # LeakageOperator-constructs-first ordering hack):
        if self.ndim != 1:
            raise ValueError(
                f"DiffusionMesh supports 1-D meshes today (slab / "
                f"cylinder / sphere); got a {self.ndim}-D mesh. The N-D "
                f"diffusion stencil is a deliberate extension seam, not "
                f"a broadcast."
            )
        # ``self.mesh`` is a ``Mesh1D`` here: every d≤2 constructor carries
        # a mesh (``mesh is None`` means d≥3, refused one gate above). The
        # mesh-less infinite-medium 1-cell carrier a second gate once
        # refused here retired at the CS4c coda (2026-09-08); that arm had
        # no reachable input — the ``ndim`` gate fired first on every
        # mesh-less carrier the tree can build.

        # The scalar boundary-trace space: the face inventory from
        # ``face_labels`` (the pole of a radial axis is not a face and
        # never appears), each face's surface measure from ``areas``
        # (slab 1, cylinder 2πR, sphere 4πR²) as the trace metric.
        labels = _face_labels(self.axes)
        faces = [
            (label.face_name, _face_shape(self.axes, label))
            for label in labels
        ]
        face_areas = {
            label.face_name: float(
                self.areas[0] if label.endpoint == "min" else self.areas[-1]
            )
            for label in labels
        }
        self._scalar_trace = ScalarTraceSpace.for_faces(
            faces, self.ng, face_areas,
        )

        # Bind the declared exterior conditions to THIS carrier's faces (its ``face_roles`` and its colouring ``symmetry``, read
        # inside ``bind`` once), then realize the binding through the ONE shared TransportMethod body (#290 P7b): the face loop
        # over the binding's faces is method-generic, and :meth:`realize_boundary_role` below is the diffusion arm. The face
        # inventory with a trace IS the ``bc`` inventory by construction: operator coverage ≡ trace-face coverage, with no
        # separate validation to drift.
        self.binding: "Binding" = bind(self, boundary_conditions)
        self.bc: "dict[str, LinearOperator]" = realize_boundary(self)

    # ── Promotion (the data/behavior join) ────────────────────────────

    @classmethod
    def from_material_mesh(
        cls,
        material_mesh: MaterialMesh,
        *,
        boundary_conditions: "BoundaryConditions | Mapping[int, BoundaryCondition] | None" = None,
    ) -> "DiffusionMesh":
        r"""Promote a :class:`MaterialMesh` to a solvable diffusion phase space.

        The diffusion arm of the data/behavior join
        (``SNProblem.from_material_mesh`` parity) — and the natural
        consumer of cross-section homogenization: a homogenized
        :class:`MaterialMesh` promotes here to re-solve the coarsened
        problem in diffusion theory. The axes carry only the TAGS; the
        CONDITIONS they name are the method's input, given by
        ``boundary_conditions``. ``None`` (the default) means "the carrier's
        own": a carrier that has a binding (an
        :class:`~orpheus.sn.problem.SNProblem`, a :class:`DiffusionMesh`)
        passes its ``binding.declared``, so its conditions are re-read by the
        diffusion table; any other carrier gets the empty set (every tag then
        raises ``UnboundTagError``, a tagged carrier must say what its tags
        mean). An explicit value always wins. The axes / mesh / mat_map /
        materials pass through verbatim; ``_init_core`` re-derives the data
        block bit-identically from them.

        An :class:`~orpheus.sn.problem.SNProblem` promotes directly (it IS a
        MaterialMesh) — the DSA construction path (#2): the diffusion
        operator family assembles on the promoted mesh, sharing geometry and
        the declared conditions with the SN sweep. A condition the diffusion
        table does not admit (a source, a translation deck) raises
        ``BoundaryRefused`` at promotion.
        """
        if boundary_conditions is None:
            # A method-mesh (the TransportMethod Protocol, structural) carries a binding; a bare carrier carries tags only.
            boundary_conditions = material_mesh.binding.declared if isinstance(material_mesh, TransportMethod) else BoundaryConditions.empty()
        obj = cls.__new__(cls)
        obj._init_core(
            axes=material_mesh.axes,
            mesh=material_mesh.mesh,
            mat_map=material_mesh.mat_map,
            materials=material_mesh.materials,
            sigma_t_cell=material_mesh.sigma_t_cell,
            boundary_conditions=boundary_conditions,
        )
        return obj

    # ── Boundary condition realization ─────────────────────────────────
    #
    # The face loop (the binding's faces, in inventory order, a deck or the
    # condition a tag is bound to) lives in the ONE shared TransportMethod
    # body, :func:`~orpheus.transport.method.realize_boundary` (#290 P7b — it
    # replaced the twin ``SNProblem._resolve_bcs`` / ``DiffusionMesh._resolve_bcs``
    # loops). Only the genuinely diffusion-specific arm remains here:

    def realize_boundary_role(self, role: "DeckTransformation | BoundaryCondition", face: str) -> "LinearOperator":
        r"""Realize one declaration (a deck, or the condition the face's tag is bound to) on ``face`` — the diffusion arm
        of the :class:`~orpheus.transport.method.TransportMethod` hook.

        Called per face with a trace by the shared :func:`~orpheus.transport.method.realize_boundary` body. Builds the
        :class:`~orpheus.diffusion.method_space.DiffusionMethodSpace` for the face (validated against the trace's inventory)
        and hands the declaration to :attr:`boundary_realizers` ``.realize``, which admits first (a declaration the diffusion
        table does not admit on this face, a source or a translation deck, raises
        :class:`~orpheus.transport.boundary_conditions.BoundaryRefused` at construction) and returns the operator
        :math:`\mathcal{A}` on the outflow partial current (1-arg ``apply``, consumed by
        :class:`~orpheus.diffusion.operators.DiffusionBoundaryOperator`). A mixed wall realizes as the functor's SUM of its
        chains' scalars; each summand is a leaf-stamped scaled identity, so the sum carries
        :attr:`~orpheus.numerics.operator.BlockRole.BOUNDARY` by the numerics role join and the stamp lives at ONE site
        (:func:`~orpheus.diffusion.boundary_realizer._albedo_operator`).
        """
        method_space = DiffusionMethodSpace.for_face(trace=self.scalar_trace, face=face)
        return self.boundary_realizers.realize(role, self.face_site(face), method_space)

    def with_cross_sections(self, sigma_t_cell: np.ndarray) -> "DiffusionMesh":
        r"""A NEW diffusion phase space over the same mesh with the per-cell :math:`\sigma_t` datum replaced.

        The re-spelling carries the DECLARED conditions (``boundary_conditions=self.binding.declared``), as
        :meth:`SNProblem.with_cross_sections <orpheus.sn.problem.SNProblem.with_cross_sections>` does: the base
        spelling re-runs the constructor without them, and a carrier whose tags are unbound is refused.
        """
        assert self.mesh is not None  # the 1-D gate of ``_init_core`` refused every carrier without a mesh (type narrowing, not a guard)  # noqa: S101
        return type(self)(self.mesh, self.materials, sigma_t_cell=sigma_t_cell, **self._respelling_keywords())

    # ── identity: the declared conditions are part of the phase space ──
    # The key is :class:`~orpheus.transport.method.MethodMesh`'s (the material mesh's key plus the declared CONDITIONS in tag
    # order), the ONE spelling that SNProblem shares: two problems over the same mesh with different conditions are different
    # phase spaces, as the laws in the axes' key were before the carve.

    # ── The method layer's spaces ──────────────────────────────────────

    @property
    def scalar_trace(self) -> "ScalarTraceSpace":
        r"""The scalar boundary-trace space — per-face ``(J⁺, J⁻)`` pairs.

        Built once at construction from the mesh's own axis data
        (``SNProblem.angular_trace`` locus parity, quadrature-free):
        the face inventory from
        :func:`~orpheus.mesh.axis.face_labels`, each face's
        surface measure from :attr:`~MaterialMesh.areas` as the trace
        metric. Consumed by the scalar trace family
        (:class:`~orpheus.transport.fields.scalar_boundary_flux.ScalarBoundaryFlux`
        and siblings) and the boundary machinery (#290 P2–P4).
        """
        return self._scalar_trace

    @cached_property
    def full_field_space(self) -> "FullFieldSpace":
        r"""The scalar composite carrier :math:`V_{\rm bulk} \oplus V_{\rm trace}`.

        The function space of the scalar-composite operator family —
        the loss :math:`A = L + C - S - B` and every bulk ⊕ boundary
        composite whose bulk is a
        :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux` and
        whose boundary is the ``(J⁺, J⁻)`` scalar trace. The exact
        mirror of :attr:`SNProblem.full_field_space
        <orpheus.sn.problem.SNProblem.full_field_space>` with
        the angular measure integrated out — the block-diagonal Hilbert
        metric :math:`G` is

        * **bulk** :math:`G_{\rm bulk} = V_{\rm cell}` — the spatial
          volume measure, carried by the spatial axis of the carrier's
          cached :attr:`bulk_space` mint (the scalar flux is already
          angle-integrated: no ``w_n``);
        * **trace** — the face-AREA metric already carried by
          :attr:`scalar_trace` (the surface measure of angle-integrated
          partial currents).

        The interior member IS the carrier's cached axis-built
        :attr:`bulk_space` (CS4b S4 — the Q2 three-family unification's
        diffusion leg, forced by the F4 derived-space form: the mint's
        members must BE the mints the FIELDS ride, or the composite a
        field derives never content-equals the carrier's; `[M]` the
        retired hand-built ``"scalar_bulk"`` twin was metric-bit-identical
        to ``bulk_space``, so the re-point moved nothing but the name).

        Cached: immutable for a given mesh.
        """
        from orpheus.numerics.spaces.full_field_space import FullFieldSpace

        return FullFieldSpace.from_blocks(self.bulk_space, self.scalar_trace)
