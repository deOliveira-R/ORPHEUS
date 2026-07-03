r"""The scalar boundary-trace function space — partial currents on a mesh boundary.

The quadrature-free sibling of
:class:`~orpheus.numerics.spaces.trace_space.TraceSpace` for methods whose
boundary state is already angle-integrated (diffusion; CP/MoC scalar traces
to follow). Where the SN trace stores the angular flux DENSITY per ordinate
under the partial-current metric :math:`G_s = |\Omega\cdot\hat n|\odot w_n`,
the scalar trace stores the two half-range ℓ=0 moments themselves — the
partial currents

.. math::

   J^+_g \;=\; \int_{\Omega\cdot\hat n > 0} (\Omega\cdot\hat n)\,
               \psi_g \, d\Omega,
   \qquad
   J^-_g \;=\; \int_{\Omega\cdot\hat n < 0} |\Omega\cdot\hat n|\,
               \psi_g \, d\Omega,

face-local against the OUTWARD normal (:math:`J^+` leaves the domain,
:math:`J^-` enters — at every face). Under the P1 closure they carry the
Cauchy data of the elliptic operator:
:math:`\phi_\Gamma = 2(J^+ + J^-)`, :math:`J = J^+ - J^-` — one output DOF
and one input DOF per face per group, which is exactly what makes an
inconsistent Cauchy pair unrepresentable and every boundary law a member of
the albedo family :math:`J^- = \mathcal{A}\,J^+` (see
``.claude/plans/diffusion_crosswalk.md``, the P2 convention contract).

**Deliberately NOT a ``TraceSpace`` subclass.** The angular trace's
``omega_dot_n`` table and quadrature coupling are SN's type-safety — making
them optional to admit a scalar realization would let an SN boundary field
be built on a table-less space (an illegal state the current type forbids).
The two traces share only the :class:`~orpheus.numerics.face_layout.FaceLayout`
flat-buffer discipline, which lives one level down and is reused unchanged.

**Metric.** The inner-product weight is the boundary face AREA (slab 1,
cylinder :math:`2\pi R`, sphere :math:`4\pi R^2`) — the surface measure, so
the trace pairing is the boundary integral :math:`\oint_\Gamma x\,y\,dA` of
already-angle-integrated quantities. The angular weights are NOT re-applied
here: they are integrated out of :math:`J^\pm` by construction; the SN→scalar
reduction (the DSA restriction, #2) is owned by the moment operator that
performs it, not by this metric.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import ClassVar, Mapping, Optional, Sequence

import numpy as np

from orpheus.numerics.face_layout import FaceLayout
from orpheus.numerics.space import FunctionSpace


__all__ = ["ScalarTraceSpace"]


@dataclass(frozen=True)
class ScalarTraceSpace(FunctionSpace):
    r"""Boundary space of per-face ``(J⁺, J⁻)`` partial-current pairs.

    One concrete space for the whole boundary :math:`\Gamma` of a material
    mesh. Each face slot has shape ``(2, ng, *face_spatial)`` — component
    axis 0 (:attr:`OUTFLOW_ROW` = :math:`J^+`, :attr:`INFLOW_ROW` =
    :math:`J^-`), group axis 1. The slot-shape and row conventions are
    owned HERE (single site — the crosswalk's storage-layout table);
    :meth:`for_faces` is the only constructor that builds them.

    Parameters
    ----------
    name, shape, inner_product_weights
        Inherited from :class:`~orpheus.numerics.space.FunctionSpace`.
        ``name`` is ``"scalar_trace"``, ``shape`` the whole-boundary flat
        shape ``(layout.total_size,)``, and the metric the per-face AREA
        broadcast over each slot (see module docstring).
    layout : FaceLayout
        The flat-buffer descriptor (which faces exist, per-face shapes,
        offsets). ``compare=False`` leaf-data — identity stays
        ``(name, shape)``, matching the angular ``TraceSpace`` convention.
    """

    #: Component-row convention of every face slot (axis 0). Row 0 is the
    #: OUTFLOW partial current J⁺ (leaves the domain through the outward
    #: normal), row 1 the INFLOW J⁻. Single source of truth — consumers
    #: (``PartialCurrent`` views, the P3 realizer, the P4 boundary
    #: operator) index through these names, never literal 0/1.
    OUTFLOW_ROW: ClassVar[int] = 0
    INFLOW_ROW: ClassVar[int] = 1

    layout: Optional[FaceLayout] = field(
        default=None, kw_only=True, repr=False, compare=False,
    )

    @property
    def face_names(self) -> tuple[str, ...]:
        """Ordered face names (flat-buffer concatenation order)."""
        if self.layout is None:
            raise RuntimeError(
                "ScalarTraceSpace has no FaceLayout; build it via "
                "ScalarTraceSpace.for_faces, not the bare constructor."
            )
        return tuple(self.layout.faces)

    @classmethod
    def for_faces(
        cls,
        faces: Sequence[tuple[str, tuple[int, ...]]],
        ng: int,
        face_areas: Mapping[str, float],
    ) -> "ScalarTraceSpace":
        r"""Build the scalar trace space from a face inventory.

        Parameters
        ----------
        faces : sequence of ``(face_name, face_spatial_shape)``
            Ordered boundary-face inventory — canonically derived from
            :func:`orpheus.transport.mesh.axis.face_labels` /
            :func:`~orpheus.transport.mesh.axis.face_shape` on the mesh's
            own axes (the pole of a radial axis is not a face and never
            appears here). ``face_spatial_shape`` is ``()`` for a 1-D
            mesh, ``(n_edge_cells,)`` for a 2-D edge.
        ng : int
            Number of energy groups (axis 1 of every slot).
        face_areas : mapping face_name → area
            The boundary surface measure per face (slab 1, cylinder
            :math:`2\pi R`, sphere :math:`4\pi R^2`) — becomes the
            diagonal inner-product weight over that face's slot.

        Raises
        ------
        ValueError
            If ``face_areas`` does not cover exactly the named faces, or
            any area is non-positive (a degenerate boundary measure).
        """
        names = [name for name, _ in faces]
        if set(face_areas) != set(names):
            raise ValueError(
                f"face_areas keys {sorted(face_areas)!r} do not match the "
                f"face inventory {sorted(names)!r}"
            )
        named_shapes = [
            (name, (2, int(ng), *face_spatial)) for name, face_spatial in faces
        ]
        layout = FaceLayout.from_named_shapes(named_shapes)
        weights = np.zeros((int(layout.total_size),), dtype=float)
        for name, slot in layout.faces.items():
            area = float(face_areas[name])
            if area <= 0.0:
                raise ValueError(
                    f"face_areas[{name!r}] = {area:g} must be positive — a "
                    f"boundary face with non-positive area is degenerate."
                )
            weights[slot.offset : slot.offset + slot.flat_size] = area
        return cls(
            name="scalar_trace",
            shape=(int(layout.total_size),),
            inner_product_weights=weights,
            layout=layout,
        )

    # ── Identity (mirror TraceSpace / FullFieldSpace) ─────────────────
    #
    # @dataclass(frozen=True) would regenerate __eq__/__hash__ over the
    # fields; explicit delegation restores the (name, shape) identity
    # convention. ``layout`` is already excluded via ``compare=False``.

    def __eq__(self, other: object) -> bool:
        return FunctionSpace.__eq__(self, other)

    def __hash__(self) -> int:
        return FunctionSpace.__hash__(self)

    def __repr__(self) -> str:
        return f"ScalarTraceSpace({self.name!r}, shape={self.shape})"

    def slot_shape(self, face: str) -> tuple[int, ...]:
        """The ``(2, ng, *face_spatial)`` shape of one face slot."""
        if self.layout is None:
            raise RuntimeError(
                "ScalarTraceSpace has no FaceLayout; build it via "
                "ScalarTraceSpace.for_faces, not the bare constructor."
            )
        if face not in self.layout.faces:
            raise ValueError(
                f"Unknown face {face!r}; available: {tuple(self.layout.faces)}"
            )
        return self.layout.faces[face].shape
