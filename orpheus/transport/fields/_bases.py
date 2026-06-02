r"""Storage-base ABCs for the typed transport field vocabulary (B.1).

This module is the single source of truth for the machinery that every
typed transport field used to redeclare leaf-by-leaf (Cardinal Rule 2).
Before B.1, ``AngularFlux`` and ``AngularSourceSink`` each independently
carried an identical ``mesh`` field, ``(N, ng, nx, ny)`` shape-check,
mesh-binding ``_check_partner``, ``from_mesh``/``from_ndarray``, and
``N/ng/nx/ny`` read-throughs — and ``ScalarFlux`` / ``ScalarSourceSink``
mirrored the same pattern on ``(ng, nx, ny)``. The repetition IS the
architectural smell; these bases are the consolidation.

The storage × role × locus hierarchy
=====================================

The field vocabulary (issues #205 / #201) is a grid of **locus**
(bulk vs boundary) × **storage family** (angular / scalar / moment /
trace) × **role** (flux / source / residual). This module provides the
*locus + storage-family* axis as ABCs; the *role* leaves (``AngularFlux``,
``AngularSourceSink``, ...) sit beneath them::

    Field (numerics, L1 — values + space + dunder algebra)
     ├─ BulkField (ABC)           mesh + mesh-binding + ng/nx/ny + abstract _phase_space_shape
     │   ├─ AngularField (ABC)    + N + (N,ng,nx,ny) from_mesh/_ndarray, parametrized by _SPACE_NAME
     │   │   ├─ AngularFlux           role leaf  (flux)
     │   │   └─ AngularSourceSink     role leaf  (source; renamed from PerOrdinateSource in B.2)
     │   ├─ ScalarField (ABC)     + (ng,nx,ny) from_mesh/_ndarray, parametrized by _SPACE_NAME
     │   │   ├─ ScalarFlux            role leaf  (flux)
     │   │   └─ ScalarSourceSink       role leaf  (source; renamed from IsotropicSource in B.2)
     │   └─ MomentField (ABC)     family marker; the moment shape is leaf-specific
     │       └─ HarmonicMomentField   role leaf  (flux-only for now)
     └─ BoundaryField (ABC)       mesh + mesh-binding + TraceSpace contract + layout/face_view + factories
         ├─ BoundaryFlux              role leaf  (flux)
         ├─ BoundarySourceSink            role leaf  (source; B.3 — orpheus.transport.source_sinks)
         └─ BoundaryResidual          role leaf  (residual; B.3 — orpheus.transport.residuals)

Parametrization (no twin paths)
===============================

The per-family phase-space shape is the one abstract hook
:meth:`BulkField._phase_space_shape`, used by the shared
``__post_init__`` validator. The Angular/Scalar families additionally
expose a ``from_mesh`` classmethod parametrized by the leaf's
``_SPACE_NAME`` :class:`~typing.ClassVar` (so ``AngularFlux``'s space is
named ``"angular_flux"`` and ``AngularSourceSink``'s ``"angular_source_sink"``,
preserving the pre-B.1 space identities bit-for-bit). ``MomentField`` and
``BoundaryField`` build their spaces differently (a TensorProductSpace
keyed on ``L``; the mesh's cached ``TraceSpace``) and so do not use
``_SPACE_NAME``.

References
----------

* ``.claude/plans/field_role_typing_view_g.md`` — Phase B (field
  vocabulary), step B.1 (storage-base ABCs).
* Grand Report v3 §5.5 (Field hierarchy), §32.5 (Field primitive spec).
* ``coding-elegance`` Pattern 2 (single source of truth), Pattern 4
  (illegal states unrepresentable), Pattern 5 (build the primitive).
"""

from __future__ import annotations

from abc import abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar, Mapping

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.field import Field
from orpheus.numerics.space import FunctionSpace
from orpheus.numerics.spaces.trace_space import TraceSpace

if TYPE_CHECKING:
    from orpheus.numerics.face_layout import FaceLayout
    from orpheus.sn.geometry import SNMesh


__all__ = [
    "BulkField",
    "AngularField",
    "ScalarField",
    "MomentField",
    "BoundaryField",
]


# ═══════════════════════════════════════════════════════════════════════
# Bulk locus
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class BulkField(Field):
    r"""Bulk-locus storage base — a mesh-bound :class:`Field` on the grid.

    Carries the machinery shared by every bulk transport field: the
    :class:`~orpheus.sn.geometry.SNMesh` binding, the
    cross-mesh-arithmetic guard (Layer-1.5: even same-class same-space
    fields on *distinct* meshes are non-additive), and the ``ng/nx/ny``
    read-throughs. The per-family phase-space shape is the single
    abstract hook :meth:`_phase_space_shape`, consumed by the shared
    construction validator.

    Abstract — instantiate a concrete role leaf (``AngularFlux``,
    ``ScalarFlux``, ...). The ``mesh`` field is annotated under
    ``TYPE_CHECKING`` (L2 field, L3 mesh) and duck-typed at runtime.
    """

    mesh: "SNMesh"

    # ── Construction validation ──────────────────────────────────────

    @abstractmethod
    def _phase_space_shape(self) -> tuple[int, ...]:
        r"""The expected shape of ``values`` / ``space`` for this family.

        Implemented per storage family (``AngularField`` →
        ``(N, ng, nx, ny)``; ``ScalarField`` → ``(ng, nx, ny)``;
        ``HarmonicMomentField`` → ``(L+1, 2L+1, ng, nx, ny)``). The
        single source of truth for the shape contract.
        """
        raise NotImplementedError

    def __post_init__(self) -> None:
        super().__post_init__()  # Field: values.shape == space.shape.
        expected = self._phase_space_shape()
        if self.space.shape != expected:
            raise ValueError(
                f"{type(self).__name__}: space.shape {self.space.shape!r} "
                f"does not match the phase-space shape {expected!r}"
            )

    # ── Algebra extension (over Field) ───────────────────────────────

    def _check_partner(self, other: object) -> None:
        r"""Add the mesh-binding guard on top of Field's class/space gate.

        Two bulk fields built on DISTINCT :class:`SNMesh` instances are
        non-additive even when same-class and same-shape — the mesh
        carries per-cell geometry (volumes, edges, BC tags) that two
        same-shape meshes may disagree on. Silently mixing across them
        produces a physically meaningless result.
        """
        super()._check_partner(other)
        if self.mesh is not other.mesh:  # type: ignore[attr-defined]
            raise ValueError(
                f"{type(self).__name__} arithmetic across distinct SNMesh "
                f"instances is forbidden — the field is mesh-bound."
            )

    # ── Metadata read-throughs ───────────────────────────────────────

    @property
    def ng(self) -> int:
        r"""Number of energy groups (delegated to ``mesh.ng``)."""
        return self.mesh.ng

    @property
    def nx(self) -> int:
        r"""Number of cells along x (delegated to ``mesh.nx``)."""
        return self.mesh.nx

    @property
    def ny(self) -> int:
        r"""Number of cells along y (delegated to ``mesh.ny``)."""
        return self.mesh.ny


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class AngularField(BulkField):
    r"""Per-ordinate bulk family on ``(N, ng, nx, ny)``.

    The storage base for the angular role leaves (``AngularFlux`` and
    ``AngularSourceSink`` today; ``AngularResidual`` joins after B.3).
    Subclasses declare a ``_SPACE_NAME``
    :class:`~typing.ClassVar` that names the :class:`FunctionSpace`
    built by :meth:`from_mesh` (preserving each leaf's pre-B.1 space
    identity). Abstract — instantiate a concrete leaf.
    """

    #: The :class:`FunctionSpace` ``name`` for this leaf (e.g.
    #: ``"angular_flux"``). Set on each concrete role leaf; absent on
    #: this abstract base (``from_mesh`` would raise if called on it).
    _SPACE_NAME: ClassVar[str]

    @classmethod
    def _shape_for_mesh(cls, mesh: "SNMesh") -> tuple[int, int, int, int]:
        r"""The ``(N, ng, nx, ny)`` phase-space shape for ``mesh``."""
        return (mesh.quad.N, mesh.ng, mesh.nx, mesh.ny)

    def _phase_space_shape(self) -> tuple[int, ...]:
        return type(self)._shape_for_mesh(self.mesh)

    @classmethod
    def from_mesh(cls, values: NDArray, mesh: "SNMesh"):
        r"""Construct from raw values + mesh, deriving the space.

        The space is ``FunctionSpace(name=cls._SPACE_NAME,
        shape=(N, ng, nx, ny))`` — single source of truth for both the
        leaf's space identity and the construction shape.
        """
        space = FunctionSpace(
            name=cls._SPACE_NAME, shape=cls._shape_for_mesh(mesh),
        )
        return cls(values=values, space=space, mesh=mesh)

    @classmethod
    def from_ndarray(cls, arr: NDArray, mesh: "SNMesh"):
        r"""Test-ergonomic alias for :meth:`from_mesh`."""
        return cls.from_mesh(arr, mesh)

    @property
    def N(self) -> int:  # noqa: N802 — matches AngularQuadrature.N
        r"""Number of angular ordinates (delegated to ``mesh.quad.N``)."""
        return self.mesh.quad.N


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class ScalarField(BulkField):
    r"""Scalar bulk family on ``(ng, nx, ny)``.

    The storage base for the scalar role leaves (``ScalarFlux`` and
    ``ScalarSourceSink`` today; ``ScalarResidual`` joins after B.3).
    Parametrized by the leaf's ``_SPACE_NAME``. Abstract — instantiate
    a concrete leaf.
    """

    #: The :class:`FunctionSpace` ``name`` for this leaf (e.g.
    #: ``"scalar_flux"``). Set on each concrete role leaf.
    _SPACE_NAME: ClassVar[str]

    @classmethod
    def _shape_for_mesh(cls, mesh: "SNMesh") -> tuple[int, int, int]:
        r"""The ``(ng, nx, ny)`` phase-space shape for ``mesh``."""
        return (mesh.ng, mesh.nx, mesh.ny)

    def _phase_space_shape(self) -> tuple[int, ...]:
        return type(self)._shape_for_mesh(self.mesh)

    @classmethod
    def from_mesh(cls, values: NDArray, mesh: "SNMesh"):
        r"""Construct from raw values + mesh, deriving the space."""
        space = FunctionSpace(
            name=cls._SPACE_NAME, shape=cls._shape_for_mesh(mesh),
        )
        return cls(values=values, space=space, mesh=mesh)

    @classmethod
    def from_ndarray(cls, arr: NDArray, mesh: "SNMesh"):
        r"""Test-ergonomic alias for :meth:`from_mesh`."""
        return cls.from_mesh(arr, mesh)


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class MomentField(BulkField):
    r"""Moment-space bulk family (ABC — see ``HarmonicMomentField``).

    A thin family marker under :class:`BulkField`: the moment indexing
    and shape are leaf-specific (a real-spherical-harmonic moment field
    carries the ``(L+1, 2L+1, ng, nx, ny)`` layout keyed on a truncation
    order ``L``), so :meth:`BulkField._phase_space_shape` stays abstract
    here and is implemented on the leaf. The marker exists so the
    ``bulk: BulkField`` typing and the storage×role×locus grid are
    complete; a second moment representation would lift its shared
    machinery here (``feedback_unify_after_two_instances``).
    """


# ═══════════════════════════════════════════════════════════════════════
# Boundary locus
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class BoundaryField(Field):
    r"""Boundary-locus storage base — a mesh-bound :class:`Field` on the
    unified :class:`~orpheus.numerics.spaces.trace_space.TraceSpace`.

    Carries the machinery shared by every boundary trace field
    (``BoundaryFlux``, ``BoundarySourceSink``, ``BoundaryResidual`` — the
    latter two minted in B.3): the :class:`SNMesh` binding + cross-mesh guard, the
    TraceSpace contract (the space IS the trace and carries the
    :class:`~orpheus.numerics.face_layout.FaceLayout`), the read-through
    :attr:`layout` property, per-face :meth:`face_view` access, and the
    :meth:`zeros_for_sn_mesh` / :meth:`from_face_arrays` factories (all
    via ``cls`` so they construct the concrete subclass).

    Storage is a SINGLE flat backing buffer (shape
    ``(layout.total_size,)``); per-face access is slice-view, no copies.
    Inherits Field's same-class/same-space dunder algebra. Abstract —
    instantiate a concrete role leaf.
    """

    mesh: "SNMesh"

    # ── Construction validation ──────────────────────────────────────

    def __post_init__(self) -> None:
        super().__post_init__()
        # A.5: the space IS the unified TraceSpace and carries the
        # FaceLayout. Enforce the contract (illegal-states-unrepresentable):
        # a boundary field on a layout-less space cannot do face_view.
        if not isinstance(self.space, TraceSpace) or self.space.layout is None:
            raise TypeError(
                f"{type(self).__name__} requires a TraceSpace carrying a "
                f"FaceLayout (A.5 re-home); got space={self.space!r}. Build "
                f"via {type(self).__name__}.zeros_for_sn_mesh / "
                f"from_face_arrays, or pass mesh.trace as the space."
            )
        expected = (self.space.layout.total_size,)
        if self.values.shape != expected:
            raise ValueError(
                f"{type(self).__name__}: values.shape {self.values.shape!r} "
                f"does not match (layout.total_size,) = {expected!r}"
            )
        if self.space.shape != expected:
            raise ValueError(
                f"{type(self).__name__}: space.shape {self.space.shape!r} "
                f"does not match (layout.total_size,) = {expected!r}"
            )

    # ── Algebra extension (over Field) ───────────────────────────────

    def _check_partner(self, other: object) -> None:
        super()._check_partner(other)
        # Mesh-binding override — two boundary fields can share a space
        # (``"sn_trace"``, same total_size — TraceSpace.__eq__ is on
        # ``(name, shape)``) but differ in mesh identity.
        if self.mesh is not other.mesh:  # type: ignore[attr-defined]
            raise ValueError(
                f"{type(self).__name__} arithmetic across distinct SNMesh "
                f"instances is forbidden — the field is mesh-bound."
            )
        if self.layout is not other.layout:  # type: ignore[attr-defined]
            # Belt-and-suspenders: operands sourced from the same mesh
            # share the cached ``mesh.trace`` (one layout identity), so
            # this fires only for hand-built operands. Fall back to
            # structural equality — same face names + shapes + offsets.
            if self.layout != other.layout:  # type: ignore[attr-defined]
                raise ValueError(
                    f"{type(self).__name__} layout mismatch — operands have "
                    f"structurally distinct FaceLayouts."
                )

    # ── Per-face access (slice views into the flat buffer) ───────────

    @property
    def layout(self) -> "FaceLayout":
        r"""The per-geometry :class:`FaceLayout`, read off the space.

        A.5 re-home: the layout lives on the :class:`TraceSpace`
        (:attr:`TraceSpace.layout`), not as a separate field attribute.
        This read-through property preserves the ``boundary.layout``
        access surface (``boundary.layout.faces``, ``.total_size``).
        """
        return self.space.layout  # type: ignore[attr-defined]

    def face_view(self, name: str) -> NDArray:
        r"""Return a per-face slice view into the flat backing buffer.

        The returned ndarray shares memory with :attr:`values`.

        Raises
        ------
        KeyError
            If ``name`` is not a face in this layout.
        """
        if name not in self.layout.faces:
            raise KeyError(
                f"{type(self).__name__}: no face named {name!r} in layout; "
                f"available faces: {list(self.layout.faces)!r}"
            )
        return self.layout.faces[name].slice_view(self.values)

    @property
    def face_views(self) -> Mapping[str, NDArray]:
        r"""Mapping ``{face_name: slice_view}`` for every face in the layout.

        All views memory-shared with :attr:`values`.
        """
        return {name: self.face_view(name) for name in self.layout.faces}

    # ── Construction factories (via ``cls`` — build the subclass) ────

    @classmethod
    def zeros_for_sn_mesh(cls, mesh: "SNMesh"):
        r"""Construct a zero boundary field sized to ``mesh``.

        Sources ``space = mesh.trace`` (the cached unified TraceSpace).
        """
        space = mesh.trace
        if space is None:
            raise ValueError(
                f"{cls.__name__}.zeros_for_sn_mesh: mesh has no TraceSpace "
                f"(mesh.trace is None — only trace-less 2-D cylindrical "
                f"meshes, which have no SN sweep, hit this). A boundary "
                f"field cannot be built without a boundary trace."
            )
        return cls(values=np.zeros(space.shape[0]), space=space, mesh=mesh)

    @classmethod
    def from_face_arrays(
        cls, mesh: "SNMesh", face_arrays: Mapping[str, NDArray],
    ):
        r"""Construct from per-face ndarrays, packing into the flat layout.

        ``face_arrays`` must cover EVERY face in the mesh's layout; each
        per-face ndarray's shape must match the layout slot's shape.

        Raises
        ------
        ValueError
            If ``face_arrays`` keys differ from the mesh's layout faces,
            or any per-face ndarray's shape mismatches the layout slot.
        """
        space = mesh.trace
        if space is None:
            raise ValueError(
                f"{cls.__name__}.from_face_arrays: mesh has no TraceSpace "
                f"(mesh.trace is None — trace-less 2-D cylindrical). A "
                f"boundary field cannot be built without a boundary trace."
            )
        layout = space.layout
        provided = set(face_arrays.keys())
        expected = set(layout.faces.keys())
        if provided != expected:
            raise ValueError(
                f"{cls.__name__}.from_face_arrays: face_arrays keys "
                f"{sorted(provided)!r} do not match mesh's layout faces "
                f"{sorted(expected)!r}"
            )
        flat = np.zeros(layout.total_size)
        for name, slot in layout.faces.items():
            arr = face_arrays[name]
            if arr.shape != slot.shape:
                raise ValueError(
                    f"{cls.__name__}.from_face_arrays: face {name!r} array "
                    f"shape {arr.shape!r} does not match layout slot "
                    f"shape {slot.shape!r}"
                )
            slot.slice_view(flat)[:] = arr
        return cls(values=flat, space=space, mesh=mesh)
