r"""Trace function spaces with per-face inflow / outflow masks.

A *trace space* is the :class:`FunctionSpace` that lives on the
boundary :math:`\partial\Omega` of the spatial domain, restricted to a
directional half: either the **inflow** half :math:`\Gamma_- = \{(x,
\Omega) \in \partial\Omega \times S^2 : \Omega \cdot \hat n(x) < 0\}`
or the **outflow** half :math:`\Gamma_+ = \{(x, \Omega) : \Omega \cdot
\hat n(x) > 0\}`. Here :math:`\hat n(x)` is the outward unit normal to
:math:`\partial\Omega` at the boundary point :math:`x`.

The native trace-space refinement of :func:`boundary_trace_space`
shipped in Issue 9.6 carries one additional piece of information: a
**per-face directional mask** marking which ordinates are inflow
(resp. outflow) at each face. The mask is the discrete realisation of
the predicate :math:`\mathrm{sign}(\Omega_n \cdot \hat n_f)` and is
the missing ingredient for downstream consumers that need to *cleanly*
distinguish inflow vs outflow per face:

* :class:`SNBoundaryRealizer` (Wave 5 of the
  ``transient-giggling-cake`` plan) calls
  :meth:`InflowTraceSpace.inflow_indices_for_face` to build
  per-face :class:`IncomingOrdinateMaskTensor` operators for the
  vacuum-BC sparse mask.
* :class:`BoundaryTraceLaw` (Wave 3) uses the trace space to assert
  that the source :math:`q` lives on :math:`\Gamma_-` only.
* Curvilinear Krylov solves that need to cleanly mask inflow
  ordinates per face (the current matrix-free path's missing
  piece).

Design — frozen-dataclass inheritance with ABC
==============================================

:class:`FunctionSpace` is ``@dataclass(frozen=True)`` with fields
``(name, shape, inner_product_weights)`` and ``__eq__`` /
``__hash__`` based only on ``(name, shape)``. The trace-space
subclasses preserve that pattern: they are themselves frozen
dataclasses, ABC inheritance is fine in Python ≥ 3.10. The
per-face mask is carried as an :class:`Optional` :class:`np.ndarray`
field with ``default=None, repr=False, compare=False`` so that:

* The mask does NOT pollute :meth:`__eq__` / :meth:`__hash__`. Two
  inflow trace spaces with the same ``(name, shape)`` compare equal
  regardless of which mesh-and-quadrature they were built from. This
  matches the abstract-vector-space framing where the identity of a
  trace space is the type tag, not the discretisation incident on
  it.
* Equality of the directional tag (``inflow`` vs ``outflow``) IS
  preserved through the ``name`` field — ``"trace_inflow"`` vs
  ``"trace_outflow"`` — so the two sub-types compare unequal even
  at matching shape.
* The frozen-dataclass machinery still raises
  :class:`dataclasses.FrozenInstanceError` on any attempted
  mutation of ``name`` / ``shape``.

The ``face_names`` field is also stored on the instance (with
``compare=False, repr=False``) so that
:meth:`inflow_indices_for_face` can look up the face index from a
name string.

Geometric convention
====================

For each face with outward unit normal :math:`\hat n_f`, an ordinate
:math:`\Omega_n` is:

* **Inflow** iff :math:`\Omega_n \cdot \hat n_f < -\epsilon`
  (direction points INTO the domain).
* **Outflow** iff :math:`\Omega_n \cdot \hat n_f > \epsilon`
  (direction points OUT of the domain).
* **Tangential** iff :math:`|\Omega_n \cdot \hat n_f| \leq \epsilon`
  — excluded from BOTH masks. ``\epsilon`` is a small tolerance
  (default ``1e-12``) accounting for floating-point round-off in
  quadrature nodes that are nominally axis-aligned.

Coord-system coverage
=====================

* **Mesh1D** — supported for Cartesian (slab), spherical, and
  cylindrical. All three share the same ``("left", "right")`` face
  structure; the outward normal is the unit vector along the radial
  / Cartesian-x axis, and :class:`Quadrature` is the shared
  quadrature with ``mu_x`` as the direction cosine along that axis.
* **Mesh2D** — Cartesian only. 2-D cylindrical (axisymmetric
  ``(r, z)``) is not implemented in ORPHEUS today and continues to
  raise :class:`NotImplementedError`; when a 2-D cylindrical sweep
  ships, the face-normal table will need a ``rmin`` / ``rmax`` /
  ``zmin`` / ``zmax`` entry set and the mask predicate will need to
  handle azimuthal averaging.

References
----------

* See ``.claude/plans/transient-giggling-cake.md`` "Wave 2 — Trace
  spaces with per-face inflow mask" for the multi-wave context
  (Cartesian-only ship-state) and
  ``.claude/plans/curvilinear-realizer-and-2arg-cleanup.md`` for
  the Issue #188 / C188.1+C188.2 curvilinear extension that
  lifted the original ``NotImplementedError`` on ``Mesh1D``
  with ``coord ∈ {SPHERICAL, CYLINDRICAL}``.
* Lewis, E.E. & Miller, W.F. (1993). *Computational Methods of
  Neutron Transport*. American Nuclear Society. §3.7 — boundary
  trace operators in the discrete-ordinates setting.
"""

from __future__ import annotations

from abc import ABC
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional, Sequence

import numpy as np
from numpy.typing import NDArray

from .space import FunctionSpace

if TYPE_CHECKING:
    from orpheus.geometry.mesh import Mesh1D, Mesh2D
    from orpheus.numerics.quadrature import Quadrature


__all__ = [
    "TraceSpace",
    "InflowTraceSpace",
    "OutflowTraceSpace",
]


# Tolerance below which |Ω · n| is considered tangential (the
# direction grazes the face and contributes to neither inflow nor
# outflow). Set at the same scale as Lebedev-quadrature node round-off
# for axis-aligned ordinates such as (0, 0, ±1).
_TANGENTIAL_EPS = 1e-12


# ─────────────────────────────────────────────────────────────────────
# Per-mesh face → outward-normal table
# ─────────────────────────────────────────────────────────────────────
#
# Each entry maps a face name to ``(axis_index, sign)`` where
# ``axis_index`` selects ``(mu_x, mu_y, mu_z)[axis_index]`` and
# ``sign`` is the outward-normal sign (±1). An ordinate is inflow at
# the face iff ``sign * mu[axis_index] < -eps`` (Ω · n < 0).


_MESH1D_FACE_NORMALS: dict[str, tuple[int, int]] = {
    # bc_left: outward normal is -x  → sign = -1, axis = 0
    "left": (0, -1),
    # bc_right: outward normal is +x → sign = +1, axis = 0
    "right": (0, +1),
}


_MESH2D_FACE_NORMALS: dict[str, tuple[int, int]] = {
    "xmin": (0, -1),
    "xmax": (0, +1),
    "ymin": (1, -1),
    "ymax": (1, +1),
}


def _quadrature_axis(quadrature: "Quadrature", axis: int) -> np.ndarray:
    """Return ``mu_x`` for axis=0, ``mu_y`` for axis=1, ``mu_z`` for axis=2.

    Falls back to a zero array when the requested axis is not present
    on the quadrature object — this lets 1-D :class:`Quadrature`
    (which has ``mu_x`` populated and ``mu_y == 0``) feed into the
    same predicate logic without special-casing.
    """
    if axis == 0:
        return np.asarray(quadrature.mu_x)
    if axis == 1:
        return np.asarray(getattr(quadrature, "mu_y", np.zeros_like(quadrature.mu_x)))
    if axis == 2:
        return np.asarray(getattr(quadrature, "mu_z", np.zeros_like(quadrature.mu_x)))
    raise ValueError(f"axis must be 0, 1, or 2; got {axis}")


def _build_per_face_mask(
    mesh: "Mesh1D | Mesh2D",
    quadrature: "Quadrature",
    faces: Sequence[str],
    *,
    direction: str,
    eps: float = _TANGENTIAL_EPS,
) -> NDArray:
    """Build the ``(n_faces, n_ordinates)`` boolean mask.

    Parameters
    ----------
    direction : ``"inflow"`` | ``"outflow"``
        Which half-space the mask selects.
    """
    # Local import to avoid cycle with geometry.mesh.
    from orpheus.geometry.coord import CoordSystem
    from orpheus.geometry.mesh import Mesh1D, Mesh2D

    if isinstance(mesh, Mesh1D):
        # All 1-D coord systems share the ("left", "right") face
        # structure with outward normals along the chosen axis
        # (Cartesian-x for slab; radial for spherical / cylindrical).
        # The quadrature's mu_x is the direction cosine along that
        # axis for both — the GaussLegendre1D adapter is shared.
        # The predicate sign(Ω·n) < -eps therefore applies identically
        # to all three coord systems. Issue #188 (curvilinear realizer
        # enablement) lifted the earlier NotImplementedError guard.
        face_table = _MESH1D_FACE_NORMALS
    elif isinstance(mesh, Mesh2D):
        if mesh.coord != CoordSystem.CARTESIAN:
            raise NotImplementedError(
                "InflowTraceSpace.from_mesh_and_quadrature for curvilinear "
                "Mesh2D (coord=CYLINDRICAL) is deferred until a "
                "curvilinear Krylov consumer arrives. "
                "See plan Wave 2 (transient-giggling-cake.md)."
            )
        face_table = _MESH2D_FACE_NORMALS
    else:
        raise TypeError(
            f"mesh must be Mesh1D or Mesh2D, got {type(mesh).__name__}"
        )

    n_ord = int(quadrature.N)
    mask = np.zeros((len(faces), n_ord), dtype=bool)

    for f_idx, face_name in enumerate(faces):
        if face_name not in face_table:
            raise ValueError(
                f"Unknown face name {face_name!r} for {type(mesh).__name__}; "
                f"valid: {sorted(face_table)}"
            )
        axis, sign = face_table[face_name]
        mu = _quadrature_axis(quadrature, axis)
        # Ω · n = sign * mu (n = sign * e_axis; ordinate component is mu)
        omega_dot_n = sign * mu
        if direction == "inflow":
            mask[f_idx] = omega_dot_n < -eps
        elif direction == "outflow":
            mask[f_idx] = omega_dot_n > +eps
        else:  # pragma: no cover — internal use only
            raise ValueError(
                f"direction must be 'inflow' or 'outflow'; got {direction!r}"
            )

    return mask


# ─────────────────────────────────────────────────────────────────────
# Trace space type hierarchy
# ─────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class TraceSpace(FunctionSpace, ABC):
    r"""Abstract base for the inflow / outflow trace spaces.

    Carries the per-face directional mask
    (``inflow_mask`` on :class:`InflowTraceSpace`, ``outflow_mask`` on
    :class:`OutflowTraceSpace`) plus the ordered tuple of face names
    that the mask rows correspond to. Both pieces are excluded from
    the dataclass-generated equality and hash to preserve
    :class:`FunctionSpace`'s ``(name, shape)`` identity convention.

    See module docstring for the geometric convention.
    """

    # The ABC tag is purely informational at runtime — Python lets us
    # instantiate this directly if we wanted to, but the factory
    # functions only emit the concrete subclasses below.


@dataclass(frozen=True)
class InflowTraceSpace(TraceSpace):
    r"""Inflow trace space :math:`\Gamma_- = \{(x, \Omega) \in
    \partial\Omega : \Omega \cdot \hat n(x) < 0\}`.

    The ``inflow_mask`` array carries the per-face directional
    predicate:

    .. math::

        \text{inflow\_mask}[f, n] = (\Omega_n \cdot \hat n_f < 0)

    Tangential ordinates (with :math:`|\Omega \cdot \hat n_f|` below
    a small tolerance) are FALSE in the inflow_mask; they live in
    neither :math:`\Gamma_-` nor :math:`\Gamma_+`.

    Parameters
    ----------
    name, shape, inner_product_weights
        Inherited from :class:`FunctionSpace`. ``name`` is
        ``"trace_inflow"`` to distinguish this from
        :class:`OutflowTraceSpace`.
    inflow_mask : NDArray of bool, shape ``(n_faces, n_ordinates)``
        Per-face inflow predicate. Built by
        :meth:`from_mesh_and_quadrature`.
    face_names : tuple of str
        Ordered tuple of face names corresponding to the rows of
        :attr:`inflow_mask`. Used by
        :meth:`inflow_indices_for_face` for name-based lookup.
    """

    inflow_mask: Optional[NDArray] = field(
        default=None, repr=False, compare=False,
    )
    face_names: tuple[str, ...] = field(
        default=(), repr=False, compare=False,
    )

    @classmethod
    def from_mesh_and_quadrature(
        cls,
        mesh: "Mesh1D | Mesh2D",
        quadrature: "Quadrature",
        faces: Sequence[str] | None = None,
        ng: int = 1,
    ) -> "InflowTraceSpace":
        r"""Build the per-face inflow mask from mesh face normals +
        quadrature direction cosines.

        Parameters
        ----------
        mesh : Mesh1D | Mesh2D
            Spatial mesh. :class:`Mesh1D` contributes the faces
            ``("left", "right")``; :class:`Mesh2D` contributes
            ``("xmin", "xmax", "ymin", "ymax")``.
        quadrature : Quadrature
            Angular quadrature exposing ``mu_x`` (always) and
            ``mu_y`` / ``mu_z`` (when applicable). Each is a
            shape-``(n_ordinates,)`` array of direction cosines.
        faces : Sequence[str], optional
            Restrict to a subset of faces. If ``None``, all faces
            available on the mesh type are used (2 for
            :class:`Mesh1D`, 4 for :class:`Mesh2D`).
        ng : int, default 1
            Number of energy groups. Sets the trailing axis of the
            trace-space shape ``(n_ordinates, ng)``.

        Raises
        ------
        NotImplementedError
            If ``mesh`` is a curvilinear :class:`Mesh2D`
            (``coord=CYLINDRICAL``); 2-D cylindrical sweeps are not
            implemented in ORPHEUS today. All :class:`Mesh1D`
            coord systems (Cartesian / spherical / cylindrical)
            are supported.
        """
        # Local import.
        from orpheus.geometry.mesh import Mesh1D, Mesh2D

        if faces is None:
            if isinstance(mesh, Mesh1D):
                faces = ("left", "right")
            elif isinstance(mesh, Mesh2D):
                faces = ("xmin", "xmax", "ymin", "ymax")
            else:
                raise TypeError(
                    f"mesh must be Mesh1D or Mesh2D, got {type(mesh).__name__}"
                )
        face_names = tuple(faces)

        mask = _build_per_face_mask(
            mesh, quadrature, face_names, direction="inflow",
        )
        return cls(
            name="trace_inflow",
            shape=(int(quadrature.N), int(ng)),
            inflow_mask=mask,
            face_names=face_names,
        )

    def inflow_indices_for_face(self, face: str) -> np.ndarray:
        """Return the ordinate indices for which ``inflow_mask[face]`` is True.

        Parameters
        ----------
        face : str
            Face name (must be one of :attr:`face_names`).

        Returns
        -------
        np.ndarray of int
            1-D integer array of ordinate indices. Empty if no
            ordinate is inflow at the requested face.
        """
        if self.inflow_mask is None:
            raise RuntimeError(
                "InflowTraceSpace was constructed without an inflow_mask; "
                "use InflowTraceSpace.from_mesh_and_quadrature() instead of "
                "the bare dataclass constructor."
            )
        try:
            f_idx = self.face_names.index(face)
        except ValueError as exc:
            raise ValueError(
                f"Unknown face {face!r}; available: {self.face_names}"
            ) from exc
        return np.flatnonzero(self.inflow_mask[f_idx])


@dataclass(frozen=True)
class OutflowTraceSpace(TraceSpace):
    r"""Outflow trace space :math:`\Gamma_+ = \{(x, \Omega) \in
    \partial\Omega : \Omega \cdot \hat n(x) > 0\}`.

    The ``outflow_mask`` array carries the per-face directional
    predicate:

    .. math::

        \text{outflow\_mask}[f, n] = (\Omega_n \cdot \hat n_f > 0)

    Tangential ordinates are FALSE in the outflow_mask (excluded
    from both halves). See :class:`InflowTraceSpace` for the
    geometric convention and the per-face normal table.

    Parameters
    ----------
    name, shape, inner_product_weights
        Inherited from :class:`FunctionSpace`. ``name`` is
        ``"trace_outflow"``.
    outflow_mask : NDArray of bool, shape ``(n_faces, n_ordinates)``
        Per-face outflow predicate.
    face_names : tuple of str
        Ordered face names corresponding to ``outflow_mask`` rows.
    """

    outflow_mask: Optional[NDArray] = field(
        default=None, repr=False, compare=False,
    )
    face_names: tuple[str, ...] = field(
        default=(), repr=False, compare=False,
    )

    @classmethod
    def from_mesh_and_quadrature(
        cls,
        mesh: "Mesh1D | Mesh2D",
        quadrature: "Quadrature",
        faces: Sequence[str] | None = None,
        ng: int = 1,
    ) -> "OutflowTraceSpace":
        """Build the per-face outflow mask. See
        :meth:`InflowTraceSpace.from_mesh_and_quadrature` for
        argument semantics."""
        from orpheus.geometry.mesh import Mesh1D, Mesh2D

        if faces is None:
            if isinstance(mesh, Mesh1D):
                faces = ("left", "right")
            elif isinstance(mesh, Mesh2D):
                faces = ("xmin", "xmax", "ymin", "ymax")
            else:
                raise TypeError(
                    f"mesh must be Mesh1D or Mesh2D, got {type(mesh).__name__}"
                )
        face_names = tuple(faces)

        mask = _build_per_face_mask(
            mesh, quadrature, face_names, direction="outflow",
        )
        return cls(
            name="trace_outflow",
            shape=(int(quadrature.N), int(ng)),
            outflow_mask=mask,
            face_names=face_names,
        )

    def outflow_indices_for_face(self, face: str) -> np.ndarray:
        """Return ordinate indices for which ``outflow_mask[face]`` is True."""
        if self.outflow_mask is None:
            raise RuntimeError(
                "OutflowTraceSpace was constructed without an outflow_mask; "
                "use OutflowTraceSpace.from_mesh_and_quadrature() instead of "
                "the bare dataclass constructor."
            )
        try:
            f_idx = self.face_names.index(face)
        except ValueError as exc:
            raise ValueError(
                f"Unknown face {face!r}; available: {self.face_names}"
            ) from exc
        return np.flatnonzero(self.outflow_mask[f_idx])
