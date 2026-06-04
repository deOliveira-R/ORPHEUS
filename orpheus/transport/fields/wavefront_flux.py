r"""Typed interior face-flux cochain :class:`WavefrontFlux` (#205 / #208).

L2 typed wrapper for the SN wavefront sweep's **interior** cell-face angular
fluxes — the per-ordinate values on every interior face the wavefront
propagates across, from the inflow trace to the outflow trace. The interior
1-cochain :math:`C^1_{\rm int}`; the interior sibling of
:class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux` (the boundary
1-cochain :math:`C^1_\partial`). Together they biproduct-decompose the full
face cochain :math:`C^1 = C^1_{\rm int} \oplus C^1_\partial`.

What this kills
===============

Pre-Wave-O the interior face fluxes were RAW ephemeral numpy arrays
(``psi_x = np.zeros((N, ng, nx+1, ny))`` / ``psi_y`` at
``orpheus/sn/sweep.py`` and ``orpheus/sn/operator.py``) — flux-bearing
tensors with NO type identity (``coding-elegance`` Pattern 3, "an unnamed
quantity is evidence the physics is wrong"), and the boundary seed/absorb
was four un-named raw edge-slice copies (``psi_x[:, :, 0, :] =
boundary.face_view("xmin")``). :class:`WavefrontFlux` names the field AND
names the trace operator the sweep applies by hand:

* :meth:`seed` is the **injection** :math:`\iota_*` (boundary inflow trace →
  interior domain-edge slots);
* :meth:`absorb` is the **pullback** :math:`\iota^*` (interior domain-edge
  slots → boundary outflow trace).

The biproduct law :math:`\iota^* \circ \iota_* = \mathrm{id}` on the boundary
chain (project-after-inject: seed then absorb) IS the "absorption = identity"
fact, now provable (``test_wavefront_flux.py``) rather than coincidental.

Field + views, NOT per-face objects
====================================

Storage is a SINGLE flat backing buffer (``space.layout.total_size``); the
per-axis face fields are zero-copy reshape views (:meth:`face`). The
cross-domain-attacker REJECTED a per-face object on three grounds
(vectorization/L16 — the ``(N_oct, ng, n_diag)`` wavefront batch is the unit
of operation; the cochain frame is storage-granularity-indifferent;
biproduct consistency). Per-axis access stays a view, exactly as
:class:`BoundaryFlux` uses ``face_view``.

Single role (flux)
==================

The interior cochain is the transient off-diagonal of a per-octant
triangular factor :math:`L_{\rm oct}` that is re-formed each sweep, so it is
**flux-only** — it has no source/residual role (the role grid is a 0-cochain
cell concept; only the boundary 1-chain inherits it, via the BC residual).
``UNITS = ANGULAR_FLUX_UNITS`` (``1/(cm²·s·sr)`` — the trace is all-flux,
shared with ``AngularFlux`` / ``BoundaryFlux``).

Honest scope (do NOT oversell)
==============================

This is a representation / elegance win: a named field, a typed
:math:`\iota_*`/:math:`\iota^*`, illegal-states-unrepresentable. It does NOT
change asymptotic cost, recover the SI rate, or enable parallelism on its own
(the seed/absorb stays an inherent cheap ``O(boundary faces)`` copy at the
persistent-boundary / ephemeral-interior lifetime split — true zero-copy is
precluded because ``BoundaryFlux`` persists across SI iterations while
``WavefrontFlux`` is rebuilt each sweep). It is the clean substrate the SI
Gauss-Seidel recovery (``si_gauss_seidel_recovery.md``) then lands on.

Unify-after-two
===============

:class:`BoundaryFlux` (boundary, instance 1) and :class:`WavefrontFlux`
(interior, instance 2) share the flat-buffer + :class:`FaceLayout` +
``layout`` read-through + factory machinery. Per
``feedback_unify_after_two_instances`` the shared ``FaceField`` ABC is lifted
only AFTER both concrete instances exist (a later elegance-enforcer pass);
this leaf carries the machinery directly for now — the deliberate
second-instance build that licenses the unification.

References
==========

* ``.claude/plans/wavefront_flux_foundation.md`` (the carve plan; §0 frame,
  §1 hierarchy, §3 the type + biproduct, §4 phasing — this is Phase 1).
* ``.claude/agent-memory/cross-domain-attacker/field_role_typing_faceflux_frames.md``
  (the cochain native-frame validation).
* :class:`~orpheus.numerics.spaces.interior_face_space.InteriorFaceSpace`
  (the layout-carrying space).
* ``coding-elegance`` Pattern 1 (read-as-the-math: ``wf.seed(bf)`` /
  ``wf.absorb(bf)`` ARE :math:`\iota_*` / :math:`\iota^*`), Pattern 3 (name
  the quantity), Pattern 4 (illegal states unrepresentable).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

from numpy.typing import NDArray

from orpheus.numerics.field import Field
from orpheus.numerics.spaces.interior_face_space import InteriorFaceSpace
from orpheus.numerics.units import ANGULAR_FLUX_UNITS, Unit

if TYPE_CHECKING:
    from orpheus.numerics.face_layout import FaceLayout
    from orpheus.sn.geometry import SNMesh
    from orpheus.transport.fields.boundary_flux import BoundaryFlux


__all__ = ["WavefrontFlux"]


#: Per-axis face-field names, indexed by axis (matches
#: :data:`InteriorFaceSpace._AXIS_NAMES`).
_AXIS_NAMES: tuple[str, ...] = ("x", "y", "z")

#: Boundary face name → (axis, edge index along that axis) for the interior
#: cochain's domain-edge slots. The min edge of axis ``a`` is index ``0`` of
#: the per-axis field; the max edge is index ``-1`` (= ``dims[a]``). Drives
#: the typed :math:`\iota_*` / :math:`\iota^*` — iterating the BOUNDARY's
#: faces, so only faces that exist (slab: xmin/xmax; curvilinear: xmax;
#: 2-D: all four) are touched.
_FACE_TO_EDGE: dict[str, tuple[int, int]] = {
    "xmin": (0, 0), "xmax": (0, -1),
    "ymin": (1, 0), "ymax": (1, -1),
    "zmin": (2, 0), "zmax": (2, -1),
}


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class WavefrontFlux(Field):
    r"""L2 interior face-flux cochain :math:`C^1_{\rm int}` — the *flux* role
    leaf on an :class:`~orpheus.numerics.spaces.interior_face_space.InteriorFaceSpace`.

    Parameters
    ----------
    values : NDArray
        Flat backing buffer, shape ``(space.layout.total_size,)``.
    space : InteriorFaceSpace
        The interior face space (Euclidean inner product) carrying the
        per-axis :class:`~orpheus.numerics.face_layout.FaceLayout`.
    mesh : SNMesh
        The SN phase-space carrier (the cross-mesh-arithmetic guard).

    Notes
    -----
    Build via :meth:`zeros_on` / :meth:`from_mesh`. Like the sweep's raw
    ``psi_x`` / ``psi_y``, the backing ``values`` buffer is mutated in place
    across the wavefront walk (the ``frozen`` dataclass freezes the
    attribute binding, not the array contents — the established pattern, cf.
    ``BoundaryFlux.face_view`` write-through).
    """

    mesh: "SNMesh"

    #: Dimensional identity (View-G): the interior trace stores flux values,
    #: so ``1/(cm²·s·sr)`` — shared with ``AngularFlux`` / ``BoundaryFlux``
    #: (the face cochain is all-flux). Metadata, not the arithmetic gate.
    UNITS: ClassVar[Unit] = ANGULAR_FLUX_UNITS

    # ── Construction validation ──────────────────────────────────────

    def __post_init__(self) -> None:
        super().__post_init__()  # Field: values.shape == space.shape.
        if (
            not isinstance(self.space, InteriorFaceSpace)
            or self.space.layout is None
        ):
            raise TypeError(
                f"{type(self).__name__} requires an InteriorFaceSpace carrying "
                f"a FaceLayout; got space={self.space!r}. Build via "
                f"{type(self).__name__}.zeros_on / from_mesh."
            )
        expected = (self.space.layout.total_size,)
        if self.values.shape != expected:
            raise ValueError(
                f"{type(self).__name__}: values.shape {self.values.shape!r} "
                f"does not match (layout.total_size,) = {expected!r}"
            )

    # ── Algebra extension (over Field) ───────────────────────────────

    def _check_partner(self, other: object) -> None:
        r"""Add the mesh-binding guard on top of Field's class/space gate."""
        super()._check_partner(other)
        if self.mesh is not other.mesh:  # type: ignore[attr-defined]
            raise ValueError(
                f"{type(self).__name__} arithmetic across distinct SNMesh "
                f"instances is forbidden — the field is mesh-bound."
            )

    # ── Per-axis access (zero-copy views into the flat buffer) ───────

    @property
    def layout(self) -> "FaceLayout":
        r"""The interior :class:`FaceLayout`, read off the space (A.5)."""
        return self.space.layout  # type: ignore[attr-defined]

    @property
    def axes(self) -> tuple[int, ...]:
        r"""The active axes (one per interior face-field): ``(0,)`` for 1-D,
        ``(0, 1)`` for 2-D, ``(0, 1, 2)`` for 3-D."""
        return tuple(range(len(self.layout.faces)))

    def face(self, axis: int) -> NDArray:
        r"""Return the per-axis face-normal field as a zero-copy reshape view.

        ``face(0)`` is the x-normal field ``(N, ng, nx+1, ...)``, ``face(1)``
        the y-normal field ``(N, ng, nx, ny+1, ...)``. The returned ndarray
        shares memory with :attr:`values` (writes propagate) — the hot-path
        wavefront walk indexes these views with byte-identical fancy-indexing
        to the legacy raw ``psi_x`` / ``psi_y``.

        Raises
        ------
        KeyError
            If ``axis`` is not an active axis of this layout.
        """
        name = _AXIS_NAMES[axis]
        if name not in self.layout.faces:
            raise KeyError(
                f"{type(self).__name__}: no face-field for axis {axis} "
                f"({name!r}) in layout; available: {list(self.layout.faces)!r}"
            )
        return self.layout.faces[name].slice_view(self.values)

    # ── Trace operators ι_* (seed) / ι* (absorb) ──────────────────────

    def _edge_slot(self, face_name: str) -> tuple[NDArray, tuple]:
        r"""Locate the interior domain-edge slot matching a boundary face.

        Returns ``(face_field, index_tuple)`` such that
        ``face_field[index_tuple]`` IS the domain-edge slice that pairs with
        ``boundary.face_view(face_name)`` (xmin → axis-0 index 0, xmax →
        axis-0 index -1, ymin → axis-1 index 0, ...). The ``2 + a`` offset
        skips the ``(N, ng)`` prefix of the per-axis face field.

        Single source of truth for the edge-locating math shared by
        :math:`\iota_*` (:meth:`seed`) and :math:`\iota^*` (:meth:`absorb`)
        — a drift here would desync the injection from the projection (the
        off-diagonal-of-:math:`L_{\rm oct}` bridge; ``coding-elegance``
        Pattern 2).
        """
        a, idx = _FACE_TO_EDGE[face_name]
        fa = self.face(a)
        sl = [slice(None)] * fa.ndim
        sl[2 + a] = idx
        return fa, tuple(sl)

    def seed(self, boundary: "BoundaryFlux") -> None:
        r"""Apply :math:`\iota_*`: inject the boundary trace into the interior
        domain-edge slots, in place.

        For each face of ``boundary``'s layout, copies the boundary face
        values into this cochain's matching domain-edge slot. Iterates the
        BOUNDARY's faces, so only faces that exist are touched (slab:
        xmin/xmax; 1-D curvilinear: xmax; 2-D: all four). The injection of
        the biproduct :math:`C^1 = C^1_{\rm int} \oplus C^1_\partial`.
        """
        for name in boundary.layout.faces:
            fa, sl = self._edge_slot(name)
            fa[sl] = boundary.face_view(name)

    def absorb(self, boundary: "BoundaryFlux") -> None:
        r"""Apply :math:`\iota^*`: pull the interior domain-edge slots back
        into the boundary trace, in place.

        The dual of :meth:`seed` — copies each interior domain-edge slot into
        the matching boundary face. After the wavefront walk the edge slots
        hold the streamed outflow (the inflow ordinate slots retain the
        seed), so this delivers the raw outflow trace to ``boundary``. The
        projection of the biproduct; ``absorb`` after ``seed`` (with no walk
        between) is the identity on the boundary chain
        (:math:`\iota^* \circ \iota_* = \mathrm{id}`).
        """
        for name in boundary.layout.faces:
            fa, sl = self._edge_slot(name)
            boundary.face_view(name)[:] = fa[sl]

    # ── Construction factories ───────────────────────────────────────

    @classmethod
    def zeros_on(cls, mesh: "SNMesh") -> "WavefrontFlux":
        r"""Construct a zero interior cochain sized to ``mesh``.

        Sources ``space = InteriorFaceSpace.from_mesh(mesh)`` and delegates to
        :meth:`~orpheus.numerics.field.Field.zeros`. The per-sweep allocator
        (the layout build is O(1) — a couple of FaceSlots — so this is cheap
        to call each sweep; cf. the raw ``np.zeros`` it replaces).
        """
        space = InteriorFaceSpace.from_mesh(mesh)
        return cls.zeros(space, mesh=mesh)

    @classmethod
    def from_mesh(cls, values: NDArray, mesh: "SNMesh") -> "WavefrontFlux":
        r"""Construct from a flat interior buffer + mesh, deriving the space."""
        space = InteriorFaceSpace.from_mesh(mesh)
        return cls(values=values, space=space, mesh=mesh)
