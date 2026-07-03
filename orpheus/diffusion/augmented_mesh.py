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
* the per-face boundary laws **realized at construction** —
  :attr:`bc` maps each boundary face to the realized albedo operator
  :math:`J^- = \mathcal{A}\,J^+` produced by the
  :class:`~orpheus.diffusion.boundary_realizer.DiffusionBoundaryRealizer`
  from the mesh axes' own :class:`~orpheus.geometry.mesh.BC`
  declarations (#290 P3 semantics; ruling 3: ``"vacuum"`` MEANS Marshak
  :math:`J^- = 0`, the zero-flux Dirichlet idealization is its own
  honestly-named ``"zero_flux"`` law).

It is the structural sibling of
:class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` (mesh + materials +
quadrature + sweep machinery): ONE method-agnostic data carrier, one
method layer per transport method. The two method-meshes are the
witnesses over which the ``TransportMethod`` Protocol is minted
(#290 P7b) — conformance structural, never imported here.

The data/behavior axis (the P7a restoration)
============================================

``material_mesh.py``'s module docstring declares the axis: a
:class:`MaterialMesh` is method-agnostic transport *data*; the method
layer (quadrature, stencil, **boundary trace**, closures) is
*behavior*. #290 P4, having no diffusion method-mesh to put them on,
parked ``scalar_trace`` + the scalar composite on :class:`MaterialMesh`
itself — same file, opposite claims. P7a resolves the contradiction by
minting THIS class and reclaiming both members: **MaterialMesh does not
know what a trace is** (:class:`SNMesh` knows angular traces;
:class:`DiffusionMesh` knows scalar traces). The type-vs-property rule
is satisfied honestly: two non-isomorphic augmentations of the same
data carrier (quadrature machinery vs trace + realized BCs), real
morphisms applied to each (promotion, BC realization, operator/field
consumption), and the illegal state "a diffusion phase space with
unresolved or unrealizable boundary conditions" is unrepresentable —
every constructible :class:`DiffusionMesh` carries realized laws for
exactly its boundary faces.

Construction semantics
======================

Construction mirrors ``SNMesh._init_core`` order: the method-agnostic
data block first (:meth:`MaterialMesh._init_data` — bit-identical), then
the method layer's admission gates (1-D only today — slab / cylinder /
sphere through the mesh's own areas + volumes; a bounded geometry —
the mesh-less infinite-medium carrier has no boundary to trace), then
the trace, then the realized ``bc`` dict. Every gate fires AT
CONSTRUCTION, not lazily inside a solver step — operators built on a
bad phase space are action-at-a-distance otherwise.

The promotion classmethod :meth:`from_material_mesh` takes **no extra
parameters**: the boundary conditions come from the axes' own ``BC``
tags, already part of the data carrier. Contrast
``SNMesh.from_material_mesh(mm, quadrature, scheme=...)`` — an honest
asymmetry (diffusion has no quadrature to inject); the promotion
signatures do not unify and should not.

Because :class:`SNMesh` *is a* :class:`MaterialMesh`, an SN phase space
promotes directly — ``DiffusionMesh.from_material_mesh(sn_mesh)``
builds the diffusion phase space over the SAME axes / materials /
mat_map, realizing the diffusion reading of the same physical BC tags.
That is exactly the DSA construction path (#2): :math:`A_{\rm diff}`
assembles on the promoted mesh, sharing the geometry with the SN sweep
it accelerates.

Layer (``tests/test_layer_imports.py``): L3 ``diffusion``. Imports
``geometry`` (laws + BC tags), ``numerics`` (spaces), ``transport``
(the MaterialMesh base + axis helpers), and its ``diffusion`` siblings
(realizer + method space).
"""

from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from orpheus.diffusion.boundary_realizer import DiffusionBoundaryRealizer
from orpheus.diffusion.method_space import DiffusionMethodSpace
from orpheus.geometry.boundary import (
    AlbedoBoundary,
    ReflectiveBoundary,
    VacuumInflow,
    ZeroFluxBoundary,
)
from orpheus.geometry.mesh import BC, Mesh1D, Mesh2D
from orpheus.numerics.face_layout import AXIS_NAMES
from orpheus.numerics.spaces.scalar_trace_space import ScalarTraceSpace
from orpheus.transport.mesh.axis import (
    axes_from_legacy_mesh,
    face_labels as _face_labels,
    face_shape as _face_shape,
)
from orpheus.transport.mesh.material_mesh import MaterialMesh

if TYPE_CHECKING:
    from orpheus.data.macro_xs.mixture import Mixture
    from orpheus.geometry.boundary import BoundaryTraceLaw
    from orpheus.numerics.operator import LinearOperator
    from orpheus.numerics.spaces.full_field_space import FullFieldSpace
    from orpheus.transport.mesh.axis import Axis1D, FaceLabel


__all__ = ["DiffusionMesh"]


class DiffusionMesh(MaterialMesh):
    r"""Diffusion method-mesh: MaterialMesh + scalar trace + realized BCs.

    See the module docstring for the data/behavior axis, the
    construction semantics, and the promotion / DSA paths.

    Parameters
    ----------
    mesh : Mesh1D
        Base geometry + zoning + per-face ``BC`` declarations (1-D;
        slab / cylinder / sphere through its own areas + volumes — a
        multi-D mesh is refused at construction). The same inbound
        surface as :class:`MaterialMesh`; a ``Mesh2D`` is accepted by
        the signature for a uniform refusal message, never admitted.
    materials : dict mapping material id to Mixture
        Macroscopic cross sections keyed by the mesh's ``mat_ids``
        values; the diffusion coefficient reads through the #290 P1
        seam (``Mixture.diffusion_coefficient``).

    Attributes
    ----------
    scalar_trace : ScalarTraceSpace
        The per-face ``(J⁺, J⁻)`` boundary-trace space (face-AREA
        metric; a curvilinear pole is not a face and never appears).
    bc : dict[str, LinearOperator]
        Per-face REALIZED boundary laws (``SNMesh.bc`` parity), keyed
        by face name — the albedo operators :math:`J^- = \mathcal{A}
        J^+`. Built from the SAME ``face_labels`` inventory as the
        trace, so law coverage ≡ face coverage by construction.
    """

    #: BC-tag → boundary-law class — the diffusion method's registry
    #: (the ``SNMesh.BOUNDARY_OPERATOR_REGISTRY`` precedent). ``"white"``
    #: is deliberately absent: at the P1 level white coincides with
    #: reflective (the P3 realizer's coincidence note) — declare
    #: ``reflective`` or ``albedo``.
    BOUNDARY_OPERATOR_REGISTRY: "ClassVar[dict[str, type[BoundaryTraceLaw]]]" = {
        "vacuum": VacuumInflow,
        "reflective": ReflectiveBoundary,
        "albedo": AlbedoBoundary,
        "zero_flux": ZeroFluxBoundary,
    }

    def __init__(
        self,
        mesh: Mesh1D | Mesh2D,
        materials: "dict[int, Mixture]",
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
        )

    def _init_core(
        self,
        *,
        axes: "tuple[Axis1D, ...]",
        mesh: Mesh1D | Mesh2D | None,
        mat_map: np.ndarray | None,
        materials: "dict[int, Mixture]",
    ) -> None:
        r"""The ONE construction body both surfaces funnel into
        (``SNMesh._init_core`` parity).

        Data block first, then the diffusion method layer: admission
        gates → trace → realized boundary laws.
        """
        # ── Method-agnostic DATA block → MaterialMesh base ──
        MaterialMesh._init_data(
            self, axes=axes, mesh=mesh, mat_map=mat_map, materials=materials,
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
        if self.mesh is None:
            raise ValueError(
                "DiffusionMesh requires a bounded geometry (a Mesh1D "
                "with boundary faces); got a mesh-less MaterialMesh "
                "(the infinite-medium 1-cell carrier). An infinite "
                "homogeneous medium has no boundary trace — use the "
                "homogeneous solver."
            )

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

        # Resolve boundary conditions from the per-axis declarations.
        self._resolve_bcs()

    # ── Promotion (the data/behavior join) ────────────────────────────

    @classmethod
    def from_material_mesh(cls, material_mesh: MaterialMesh) -> "DiffusionMesh":
        r"""Promote a :class:`MaterialMesh` to a solvable diffusion phase space.

        The diffusion arm of the data/behavior join
        (``SNMesh.from_material_mesh`` parity) — and the natural
        consumer of cross-section homogenization: a homogenized
        :class:`MaterialMesh` promotes here to re-solve the coarsened
        problem in diffusion theory. NO extra parameters: the boundary
        conditions come from the axes' own ``BC`` declarations, already
        part of the data carrier (contrast SN's quadrature + scheme —
        an honest asymmetry). The axes / mesh / mat_map / materials
        pass through verbatim; ``_init_core`` re-derives the data block
        bit-identically from them.

        An :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` promotes
        directly (it IS a MaterialMesh) — the DSA construction path
        (#2): the diffusion operator family assembles on the promoted
        mesh, sharing geometry and BC declarations with the SN sweep.
        """
        obj = cls.__new__(cls)
        obj._init_core(
            axes=material_mesh.axes,
            mesh=material_mesh.mesh,
            mat_map=material_mesh.mat_map,
            materials=material_mesh.materials,
        )
        return obj

    # ── Boundary condition resolution ──────────────────────────────────

    def _resolve_bcs(self) -> None:
        r"""Resolve per-axis-declared BCs into realized albedo operators.

        ONE loop over the same ``face_labels`` inventory the trace was
        built from (``SNMesh._resolve_bcs`` parity): the BC declaration
        for each label comes from its axis
        (``axes[label.axis_index].bc[label.endpoint]``; ``None``
        defaults to ``BC("reflective")`` — the infinite-lattice /
        eigenvalue convention), is realized by :meth:`_resolve_one`,
        and lands in :attr:`bc` under the label's face name. The face
        inventory IS the BC inventory by construction: law coverage ≡
        trace-face coverage, with no separate validation to drift.
        """
        default = BC("reflective")
        self.bc: "dict[str, LinearOperator]" = {
            label.face_name: self._resolve_one(
                self.axes[label.axis_index].bc[label.endpoint] or default,
                label,
            )
            for label in _face_labels(self.axes)
        }

    def _resolve_one(self, bc: BC, label: "FaceLabel") -> "LinearOperator":
        r"""Realize one face's BC tag: tag → typed law → albedo operator.

        Builds the :class:`~orpheus.diffusion.method_space.DiffusionMethodSpace`
        for the face (validated against the trace's inventory) and
        hands the typed law to
        :meth:`DiffusionBoundaryRealizer.realize
        <orpheus.diffusion.boundary_realizer.DiffusionBoundaryRealizer.realize>`.
        The realized operator applies :math:`\mathcal{A}` to the
        outflow partial current (1-arg ``apply`` — consumed by
        :class:`~orpheus.diffusion.operators.DiffusionBoundaryOperator`).
        """
        law = self._law_from_tag(bc, label)
        method_space = DiffusionMethodSpace.for_face(
            mesh=self, face=label.face_name,
        )
        return DiffusionBoundaryRealizer().realize(law, method_space)

    def _law_from_tag(self, bc: BC, label: "FaceLabel") -> "BoundaryTraceLaw":
        r"""Construct the typed boundary law a :class:`BC` tag declares.

        Ruling-3 semantics: ``"vacuum"`` → :class:`VacuumInflow`
        (Marshak :math:`J^- = 0`, :math:`\mathcal{A} = 0`);
        ``"zero_flux"`` → the honestly-named Dirichlet idealization
        (:math:`\mathcal{A} = -1`). A reflective law reflects across
        the face's own axis (the SN convention); at the P1 level only
        the albedo = 1 scalar survives.
        """
        law_cls = self.BOUNDARY_OPERATOR_REGISTRY.get(bc.kind)
        if law_cls is None:
            supported = ", ".join(
                f"'{k}'" for k in sorted(self.BOUNDARY_OPERATOR_REGISTRY)
            )
            raise ValueError(
                f"DiffusionMesh does not support boundary condition "
                f"'{bc.kind}' on face '{label.face_name}'. "
                f"Supported: {supported}."
            )
        if law_cls is ReflectiveBoundary:
            return ReflectiveBoundary(
                axis=AXIS_NAMES[label.axis_index], albedo=1.0,
            )
        if law_cls is AlbedoBoundary:
            try:
                return AlbedoBoundary(albedo=float(bc.params["albedo"]))
            except KeyError as exc:
                raise ValueError(
                    f"BC('albedo') on face '{label.face_name}' requires an "
                    f"'albedo' parameter; got params={bc.params!r}."
                ) from exc
        return law_cls()

    # ── The method layer's spaces ──────────────────────────────────────

    @property
    def scalar_trace(self) -> "ScalarTraceSpace":
        r"""The scalar boundary-trace space — per-face ``(J⁺, J⁻)`` pairs.

        Built once at construction from the mesh's own axis data
        (``SNMesh.angular_trace`` locus parity, quadrature-free):
        the face inventory from
        :func:`~orpheus.transport.mesh.axis.face_labels`, each face's
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
        mirror of :attr:`SNMesh.full_field_space
        <orpheus.sn.mesh.augmented_mesh.SNMesh.full_field_space>` with
        the angular measure integrated out — the block-diagonal Hilbert
        metric :math:`G` is

        * **bulk** :math:`G_{\rm bulk} = V_{\rm cell}` — the spatial
          volume measure, stored ``(1, *spatial)`` so it broadcasts
          across the energy-group axis of the ``(ng, *spatial)`` bulk
          (the scalar flux is already angle-integrated: no ``w_n``);
        * **trace** — the face-AREA metric already carried by
          :attr:`scalar_trace` (the surface measure of angle-integrated
          partial currents).

        Cached: immutable for a given mesh.
        """
        from orpheus.numerics.space import FunctionSpace
        from orpheus.numerics.spaces.full_field_space import FullFieldSpace

        V = np.asarray(self.volumes, dtype=float)  # (*spatial)
        bulk_space = FunctionSpace(
            name="scalar_bulk",
            shape=(self.ng, *self.spatial_shape),
            inner_product_weights=V.reshape((1, *V.shape)),
        )
        return FullFieldSpace.from_blocks(bulk_space, self.scalar_trace)
