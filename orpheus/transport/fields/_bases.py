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
     │       └─ HarmonicMomentFlux   role leaf  (flux-only for now)
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
from typing import TYPE_CHECKING, ClassVar, Mapping, Self

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.field import Field
from orpheus.numerics.space import FunctionSpace
from orpheus.numerics.spaces.spatial_moment_space import (
    SpatialMomentSpace,
    spatial_moment_tail,
)
from orpheus.numerics.spaces.spherical_harmonic_space import (
    SphericalHarmonicSpace,
)
from orpheus.numerics.spaces.trace_space import TraceSpace

if TYPE_CHECKING:
    from orpheus.numerics.face_layout import FaceLayout
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.transport.mesh.material_mesh import MaterialMesh


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
    :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` binding, the
    cross-mesh-arithmetic guard (Layer-1.5: even same-class same-space
    fields on *distinct* meshes are non-additive), and the ``ng/nx/ny``
    read-throughs. The per-family phase-space shape is the single
    abstract hook :meth:`_phase_space_shape`, consumed by the shared
    construction validator.

    Abstract — instantiate a concrete role leaf (``AngularFlux``,
    ``ScalarFlux``, ...). The ``mesh`` field is annotated
    :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` (the
    method-agnostic mesh+materials carrier) under ``TYPE_CHECKING`` and
    duck-typed at runtime — bulk fields read only ``ng``/``spatial_shape``/
    ``ndim`` (all ``MaterialMesh`` data), so they live on any material mesh
    including a meshless single-region one (#267 / #276). The SN-only
    ``mesh.quad`` access is confined to :class:`AngularField`, which
    ``cast``\ s to :class:`SNMesh` at the two instance reads.
    """

    mesh: "MaterialMesh"

    # ── Construction validation ──────────────────────────────────────

    @abstractmethod
    def _phase_space_shape(self) -> tuple[int, ...]:
        r"""The expected shape of ``values`` / ``space`` for this family.

        Implemented per storage family (``AngularField`` →
        ``(N, ng, nx, ny)``; ``ScalarField`` → ``(ng, nx, ny)``;
        ``HarmonicMomentFlux`` → ``(L+1, 2L+1, ng, nx, ny)``). The
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

    def _check_partner(self, other: Field) -> Self:
        r"""Add the mesh-binding guard on top of Field's class/space gate.

        Two bulk fields built on DISTINCT :class:`SNMesh` instances are
        non-additive even when same-class and same-shape — the mesh
        carries per-cell geometry (volumes, edges, BC tags) that two
        same-shape meshes may disagree on. Silently mixing across them
        produces a physically meaningless result.
        """
        partner = super()._check_partner(other)
        if self.mesh is not partner.mesh:
            raise ValueError(
                f"{type(self).__name__} arithmetic across distinct SNMesh "
                f"instances is forbidden — the field is mesh-bound."
            )
        return partner

    # ── Optional spatial-moment factor (#240 D5b-S3-A0) ──────────────

    @staticmethod
    def _compose_spatial_moments(
        space: FunctionSpace, mesh: "MaterialMesh", spatial_moments_per_axis: int,
    ) -> FunctionSpace:
        r"""Append the optional within-cell spatial-moment factor to ``space``.

        Composes a
        :class:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace`
        factor (the within-cell tensor-Legendre DG basis #240 D5b-S3) onto
        ``space`` via the tensor-product ``*`` — EXACTLY as
        :meth:`HarmonicMomentFlux.from_mesh_and_L` composes the
        angular-moment :class:`SphericalHarmonicSpace`. The factor is the
        Linear-Discontinuous closure's spatial-slope carrier that travels
        between source-iteration sweeps (the diffusion-limit-consistent
        scattering source :math:`\Sigma_s \otimes I_{\rm spatial}`).

        **Gated on ``spatial_moments_per_axis > 1`` (construct-general /
        select-narrow, #240 D5b-S3-A0).** The total within-cell moment
        count is ``spatial_moments_per_axis ** mesh.ndim``. When the count
        is ``1``, the ``spatial_moment_tail`` policy returns ``()`` and
        this method returns ``space`` UNCHANGED — the field space stays
        BYTE-IDENTICAL to its pre-S3 shape (the backward-compat invariant,
        single-sourced from
        :func:`orpheus.numerics.moment_layout.face_moment_tail` via
        :func:`~orpheus.numerics.spaces.spatial_moment_space.spatial_moment_tail`).

        ``spatial_moments_per_axis`` is an EXPLICIT parameter (the
        ``spatial_moments`` factory parameter, default ``1`` everywhere),
        NOT auto-read from ``mesh.scheme.spatial_basis_per_axis``. This is
        the construct-general / select-narrow discipline: the CAPABILITY to
        carry the axis exists, but every current call site passes the
        default ``1`` so NO production field — DD, Step, OR LD — carries
        the axis yet. The iterate / cell-emit / source seams that thread
        the scheme's ``spatial_basis_per_axis`` here (selecting the axis
        for LD) are the NEXT sub-step (S3-A). Reading the scheme by default
        would silently widen LD field shapes and break LD byte-identity
        before the consumers that fill the axis exist — a Pattern-4
        violation (an axis no producer fills is an illegal state).
        """
        n_moments = spatial_moments_per_axis ** mesh.ndim
        # "append iff > 1" — single-sourced; () at n==1 → no factor, byte-id.
        if spatial_moment_tail(n_moments) == ():
            return space
        return space * SpatialMomentSpace.from_per_axis(
            spatial_moments_per_axis, mesh.ndim,
        )

    @staticmethod
    def _spatial_moment_tail_of(space: FunctionSpace) -> tuple[int, ...]:
        r"""The trailing spatial-moment shape suffix carried by ``space``, or ``()``.

        Reads the optional
        :class:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace`
        factor OFF a composed space — the space is the single source of
        truth for the moment width, so :meth:`_phase_space_shape` derives
        the expected widened shape from here rather than re-threading the
        factory's ``spatial_moments`` parameter into a stored field
        (Angular/Scalar leaves carry no such field — only the windowed
        :class:`HarmonicMomentFlux` does, because its ``L`` field already
        breaks the uniform-signature contract).

        Returns ``()`` for a non-composed / DD-default space (no factor →
        byte-identical validation prefix), and ``(per_axis ** ndim,)`` when
        a :class:`SpatialMomentSpace` factor is present.
        """
        find_factor = getattr(space, "find_factor", None)
        if find_factor is None:
            return ()  # a bare FunctionSpace (DD default) — no factor.
        try:
            factor = find_factor(SpatialMomentSpace)
        except KeyError:
            return ()
        return factor.shape

    @property
    def spatial_moments_per_axis(self) -> int:
        r"""The within-cell spatial-moment count per axis carried by this field.

        Reads the ``per_axis`` parameter OFF the optional
        :class:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace`
        factor on this field's space (the single source of truth for the moment
        width, #240 D5b-S3-A0).  Returns ``1`` for a non-composed / DD-default
        space.  Producers that derive a moment-carrying child field (e.g.
        :meth:`AngularFlux.integrate_angular`,
        :meth:`HarmonicMomentFlux.scalar_flux`) pass this as the child's
        ``spatial_moments`` so the moment axis is propagated as a TYPED factor,
        not an opaque widened ndarray."""
        from orpheus.numerics.spaces.spatial_moment_space import SpatialMomentSpace

        find_factor = getattr(self.space, "find_factor", None)
        if find_factor is None:
            return 1
        try:
            factor = find_factor(SpatialMomentSpace)
        except KeyError:
            return 1
        return factor.per_axis

    # ── Metadata read-throughs ───────────────────────────────────────

    @property
    def ng(self) -> int:
        r"""Number of energy groups (delegated to ``mesh.ng``)."""
        return self.mesh.ng

    # C5.2 (#225): the ``nx``/``ny`` read-throughs are RETIRED — a
    # field keyed on ``(nx, ny)`` silently truncates a 3-D tensor.
    # Spatial shape reads are rank-generic: ``mesh.spatial_shape``.


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

    # Narrowed to ``SNMesh`` (covariant override of ``BulkField.mesh:
    # MaterialMesh``, #267): an angular field is per-ORDINATE, so it is
    # meaningless without a quadrature and ALWAYS lives on an ``SNMesh``. The
    # narrowing keeps ``mesh.quad`` honest here (no cast) and lets the operators
    # read ``angular_field.mesh`` as an ``SNMesh`` directly.
    mesh: "SNMesh"

    #: The :class:`FunctionSpace` ``name`` for this leaf (e.g.
    #: ``"angular_flux"``). Set on each concrete role leaf; absent on
    #: this abstract base (``from_mesh`` would raise if called on it).
    _SPACE_NAME: ClassVar[str]

    @classmethod
    def _shape_for_mesh(cls, mesh: "SNMesh") -> tuple[int, ...]:
        r"""The ``(N, ng, *spatial)`` phase-space shape for ``mesh``.

        Spatial rank equals ``mesh.ndim`` — ``(N, ng, nx)`` for a 1-D
        mesh, ``(N, ng, nx, ny)`` for 2-D, ``(N, ng, nx, ny, nz)`` for
        3-D. The spatial tail is ``mesh.spatial_shape`` (no phantom
        ``ny=1`` axis on the 1-D path).
        """
        return (mesh.quad.N, mesh.ng, *mesh.spatial_shape)

    @classmethod
    def _space_for_mesh(
        cls, mesh: "SNMesh", *, spatial_moments: int = 1,
    ) -> FunctionSpace:
        r"""The leaf's :class:`FunctionSpace` for ``mesh`` (name + shape).

        Single source of truth for the leaf's space identity, shared by
        :meth:`from_mesh` and :meth:`zeros_on`.

        ``spatial_moments`` (default ``1``) is the optional within-cell
        spatial-moment basis size per axis (#240 D5b-S3-A0). At the default
        ``1`` the space is the EXACT pre-S3
        ``FunctionSpace(name=cls._SPACE_NAME, shape=(N, ng, *spatial))`` —
        byte-identical for DD/Step AND LD (no current caller passes > 1).
        At ``> 1`` a
        :class:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace`
        factor is composed on (see :meth:`BulkField._compose_spatial_moments`).
        """
        base = FunctionSpace(
            name=cls._SPACE_NAME, shape=cls._shape_for_mesh(mesh),
        )
        return cls._compose_spatial_moments(base, mesh, spatial_moments)

    def _phase_space_shape(self) -> tuple[int, ...]:
        # Base ``(N, ng, *spatial)`` prefix + the optional spatial-moment
        # tail read off the field's own space (the single source of truth
        # for the moment width). ``()`` for a DD-default space → the
        # validator cross-checks EXACTLY the pre-S3 shape (byte-identical).
        return (
            *type(self)._shape_for_mesh(self.mesh),
            *type(self)._spatial_moment_tail_of(self.space),
        )

    @classmethod
    def from_mesh(cls, values: NDArray, mesh: "SNMesh", *, spatial_moments: int = 1):
        r"""Construct from raw values + mesh, deriving the space.

        The space is ``FunctionSpace(name=cls._SPACE_NAME,
        shape=(N, ng, *spatial))`` — single source of truth for both the
        leaf's space identity and the construction shape. ``spatial_moments``
        (default ``1``, byte-identical) optionally composes the within-cell
        spatial-moment factor (#240 D5b-S3-A0).
        """
        return cls(
            values=values,
            space=cls._space_for_mesh(mesh, spatial_moments=spatial_moments),
            mesh=mesh,
        )

    @classmethod
    def zeros_on(cls, mesh: "SNMesh", *, spatial_moments: int = 1):
        r"""Construct a zero field of this leaf sized to ``mesh`` (B.5.A).

        The bulk-locus zero factory: derives the space from ``mesh`` and
        delegates to :meth:`~orpheus.numerics.field.Field.zeros`. The
        uniform leaf-side allocator that
        :meth:`~orpheus.transport.timed_full_field.TimedFullField.zeros`
        calls; replaces the retired ``SNMesh.zeros_*`` mesh-side factories.
        ``spatial_moments`` (default ``1``, byte-identical) optionally
        composes the within-cell spatial-moment factor (#240 D5b-S3-A0).
        """
        return cls.zeros(
            cls._space_for_mesh(mesh, spatial_moments=spatial_moments), mesh=mesh,
        )

    @classmethod
    def from_ndarray(cls, arr: NDArray, mesh: "SNMesh"):
        r"""Test-ergonomic alias for :meth:`from_mesh`."""
        return cls.from_mesh(arr, mesh)

    @property
    def N(self) -> int:  # noqa: N802 — matches Quadrature.N
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
    def _shape_for_mesh(cls, mesh: "MaterialMesh") -> tuple[int, ...]:
        r"""The ``(ng, *spatial)`` phase-space shape for ``mesh``.

        Spatial rank equals ``mesh.ndim`` — ``(ng, nx)`` for a 1-D mesh,
        ``(ng, nx, ny)`` for 2-D (no phantom ``ny=1`` on the 1-D path).
        """
        return (mesh.ng, *mesh.spatial_shape)

    @classmethod
    def _space_for_mesh(
        cls, mesh: "MaterialMesh", *, spatial_moments: int = 1,
    ) -> FunctionSpace:
        r"""The leaf's :class:`FunctionSpace` for ``mesh`` (name + shape).

        Single source of truth for the leaf's space identity, shared by
        :meth:`from_mesh` and :meth:`zeros_on`.

        ``spatial_moments`` (default ``1``) is the optional within-cell
        spatial-moment basis size per axis (#240 D5b-S3-A0). At the default
        ``1`` the space is the EXACT pre-S3
        ``FunctionSpace(name=cls._SPACE_NAME, shape=(ng, *spatial))`` —
        byte-identical. The :class:`ScalarSourceSink` scattering-source
        accumulator is the carrier that selects ``> 1`` at S3-A so the
        slope rows can hold :math:`\Sigma_s \cdot \hat\phi`.
        """
        base = FunctionSpace(
            name=cls._SPACE_NAME, shape=cls._shape_for_mesh(mesh),
        )
        return cls._compose_spatial_moments(base, mesh, spatial_moments)

    def _phase_space_shape(self) -> tuple[int, ...]:
        # Base ``(ng, *spatial)`` prefix + the optional spatial-moment tail
        # read off the field's own space (single source of truth). ``()``
        # for a DD-default space → byte-identical validation.
        return (
            *type(self)._shape_for_mesh(self.mesh),
            *type(self)._spatial_moment_tail_of(self.space),
        )

    @classmethod
    def from_mesh(cls, values: NDArray, mesh: "MaterialMesh", *, spatial_moments: int = 1):
        r"""Construct from raw values + mesh, deriving the space.

        ``spatial_moments`` (default ``1``, byte-identical) optionally
        composes the within-cell spatial-moment factor (#240 D5b-S3-A0).
        """
        return cls(
            values=values,
            space=cls._space_for_mesh(mesh, spatial_moments=spatial_moments),
            mesh=mesh,
        )

    @classmethod
    def zeros_on(cls, mesh: "MaterialMesh", *, spatial_moments: int = 1):
        r"""Construct a zero field of this leaf sized to ``mesh`` (B.5.A).

        The bulk-locus zero factory: derives the space from ``mesh`` and
        delegates to :meth:`~orpheus.numerics.field.Field.zeros`. The
        uniform leaf-side allocator that
        :meth:`~orpheus.transport.timed_full_field.TimedFullField.zeros`
        calls; replaces the retired ``SNMesh.zeros_*`` mesh-side factories.
        ``spatial_moments`` (default ``1``, byte-identical) optionally
        composes the within-cell spatial-moment factor (#240 D5b-S3-A0).
        """
        return cls.zeros(
            cls._space_for_mesh(mesh, spatial_moments=spatial_moments), mesh=mesh,
        )

    @classmethod
    def from_ndarray(cls, arr: NDArray, mesh: "MaterialMesh"):
        r"""Test-ergonomic alias for :meth:`from_mesh`."""
        return cls.from_mesh(arr, mesh)


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class MomentField(BulkField):
    r"""Real-spherical-harmonic moment-space bulk family (storage base).

    The storage base for the moment role leaves —
    :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
    (the ``(FluxRole, MomentField)`` flux state) and
    :class:`~orpheus.transport.source_sinks.harmonic_moment_source_sink.HarmonicMomentSourceSink`
    (the bare source/sink). It carries the construction machinery the two
    share: the ``L`` truncation-order + ``spatial_moments`` fields, the
    ``(L+1, 2L+1, ng, *spatial[, …])`` :meth:`_phase_space_shape`, the
    ``L``-match :meth:`_check_partner`, and the
    :class:`~orpheus.numerics.space.TensorProductSpace`-building
    :meth:`from_mesh_and_L` / :meth:`zeros_for_mesh_and_L` /
    :meth:`from_ndarray` factories.

    A moment field is **method-agnostic** (a moment field on the
    spherical-harmonic ⊗ cell-group phase space, keyed on the truncation
    order ``L``), so this base — unlike the Angular/Scalar families — does
    NOT use a single ``_SPACE_NAME`` (its space is a TensorProductSpace).
    Instead the per-leaf cell-group factor name is the
    :attr:`_CELL_GROUP_NAME` :class:`~typing.ClassVar`, so the two role
    leaves carry distinct space identities (matching the
    ``AngularFlux``/``AngularSourceSink`` precedent). Abstract — instantiate
    a concrete role leaf.

    This lift happened when the second moment representation arrived
    (``feedback_unify_after_two_instances``): the machinery used to live on
    the lone ``HarmonicMomentField`` leaf; the Frame-campaign P4
    source/sink sibling triggered the "clean before extending" pass.
    """

    #: Maximum harmonic order retained. Determines the leading two axes'
    #: sizes: ``values.shape[:2] == (L+1, 2L+1)``. Encoded in ``space.shape``
    #: AND kept as a top-level field for ergonomic hot-path read access
    #: (avoids a per-read composition-tree traversal of
    #: ``space.find_factor(SphericalHarmonicSpace).L``).
    L: int

    #: Optional within-cell spatial-moment basis size per axis (#240
    #: D5b-S3-A0). Default ``1`` — the cell-average closure, byte-identical
    #: to the pre-S3 ``(L+1, 2L+1, ng, *spatial)`` shape. At ``> 1`` (the
    #: LD windowed iterate) a trailing ``spatial_moments ** ndim`` axis rides
    #: on every moment so the in-sweep ``moment_buf`` can carry the
    #: within-cell slopes the diffusion-limit-consistent scattering source
    #: needs between sweeps. Single-sourced "append iff > 1" via
    #: :func:`~orpheus.numerics.spaces.spatial_moment_space.spatial_moment_tail`.
    spatial_moments: int = 1

    # Narrowed to ``SNMesh`` (covariant override of ``BulkField.mesh:
    # MaterialMesh``, #267): a moment field IS the angular flux's
    # harmonic-moment iterate (φ_ℓ^m), an SN construct, so it always lives on
    # an ``SNMesh`` — operators read ``moment_field.mesh`` as ``SNMesh`` directly.
    mesh: "SNMesh"

    #: The :class:`FunctionSpace` ``name`` of the cell-group factor in this
    #: leaf's ``SphericalHarmonicSpace(L) ⊗ CellGroup`` space — distinguishes
    #: the flux leaf (``"cell_group"``, preserving its pre-P4 identity) from
    #: the source/sink leaf (``"cell_group_source_sink"``). The SH factor and
    #: shape are identical across leaves; the class-identity gate is the
    #: arithmetic guard, this name keeps the two spaces non-``==``.
    _CELL_GROUP_NAME: ClassVar[str] = "cell_group"

    # ── Construction validation ──────────────────────────────────────

    def _phase_space_shape(self) -> tuple[int, ...]:
        r"""The ``(L+1, 2L+1, ng, *spatial[, spatial_moments**ndim])`` moment shape.

        Implements :meth:`BulkField._phase_space_shape`; the shared
        :meth:`BulkField.__post_init__` validator consumes it. The leading
        two axes encode the harmonic truncation order ``L``; the spatial
        tail is rank-``d`` via ``mesh.spatial_shape`` (no phantom ``ny``).
        At :attr:`spatial_moments` ``> 1`` a trailing within-cell
        spatial-moment axis of length ``spatial_moments ** ndim`` is appended
        (#240 D5b-S3-A0; "append iff > 1" single-sourced from
        :func:`~orpheus.numerics.spaces.spatial_moment_space.spatial_moment_tail`,
        so the default ``1`` leaves the shape byte-identical and AGREES with
        the factor :meth:`from_mesh_and_L` composes).
        """
        n_moments = self.spatial_moments ** self.mesh.ndim
        return (
            self.L + 1, 2 * self.L + 1,
            self.mesh.ng, *self.mesh.spatial_shape,
            *spatial_moment_tail(n_moments),
        )

    # ── Algebra extension (over BulkField) ───────────────────────────

    def _check_partner(self, other: Field) -> Self:
        r"""Add the ``L``-match on top of BulkField's class/space/mesh gate.

        :meth:`BulkField._check_partner` already rejects on class identity,
        space equality, AND mesh identity (mesh-bound). This override adds an
        explicit ``L`` match for a clearer error message at the
        truncation-mismatch site (the space check also catches it via shape
        mismatch, but the message is less specific).
        """
        partner = super()._check_partner(other)
        if self.L != partner.L:
            raise ValueError(
                f"{type(self).__name__} arithmetic requires matching L; "
                f"got self.L={self.L}, other.L={partner.L}."
            )
        return partner

    # ── Construction factories ───────────────────────────────────────

    @classmethod
    def from_mesh_and_L(
        cls, values: NDArray, mesh: "SNMesh", L: int, *, spatial_moments: int = 1,
    ):
        r"""Construct from raw values + mesh + L, deriving the
        :class:`TensorProductSpace`.

        Builds the space as ``SphericalHarmonicSpace.from_L(L) * CellGroup``
        where ``CellGroup`` is a plain :class:`FunctionSpace` named
        :attr:`_CELL_GROUP_NAME` with the mesh's ``(ng, *spatial)`` shape —
        the moment-axis structure is type-visible through the composition
        tree (queryable via ``space.find_factor(SphericalHarmonicSpace).L``
        per Issue #207).

        ``spatial_moments`` (default ``1``, byte-identical #240 D5b-S3-A0)
        optionally composes a within-cell
        :class:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace`
        factor on AFTER the cell-group space — EXACTLY the same ``*``
        composition that adds the angular ``SphericalHarmonicSpace`` ("append
        iff > 1", single-sourced, matching :meth:`_phase_space_shape`).
        """
        sh_space = SphericalHarmonicSpace.from_L(L)
        cell_group_space = FunctionSpace(
            name=cls._CELL_GROUP_NAME,
            shape=(mesh.ng, *mesh.spatial_shape),
        )
        space = cls._compose_spatial_moments(
            sh_space * cell_group_space, mesh, spatial_moments,
        )
        return cls(
            values=values, space=space, mesh=mesh, L=L,
            spatial_moments=spatial_moments,
        )

    @classmethod
    def zeros_for_mesh_and_L(
        cls, mesh: "SNMesh", L: int, *, spatial_moments: int = 1,
    ):
        r"""Construct a zero moment field at order ``L`` sized to ``mesh`` (B.5.A).

        The moment-space parallel of the bulk leaves' :meth:`zeros_on`,
        mirroring :meth:`from_mesh_and_L` with a zero buffer. The extra ``L``
        makes the signature non-uniform, so this is deliberately NOT named
        ``zeros_on`` — a moment field is never a
        :class:`~orpheus.transport.timed_full_field.TimedFullField` composite
        slot, so it does not need the uniform allocator interface.

        ``spatial_moments`` (default ``1``, byte-identical #240 D5b-S3-A0)
        sizes the optional within-cell spatial-moment axis to match
        :meth:`from_mesh_and_L`.
        """
        n_moments = spatial_moments ** mesh.ndim
        values = np.zeros(
            (L + 1, 2 * L + 1, mesh.ng, *mesh.spatial_shape,
             *spatial_moment_tail(n_moments)),
        )
        return cls.from_mesh_and_L(
            values, mesh, L, spatial_moments=spatial_moments,
        )

    @classmethod
    def from_ndarray(cls, arr: NDArray, mesh: "SNMesh", L: int):
        r"""Test-ergonomic alias for :meth:`from_mesh_and_L`.

        Per Depth B plan §3.7, every typed field exposes
        ``from_ndarray(arr, ...)`` as the migration path from the retired
        ``apply(np.ndarray)`` singledispatch handlers (D-I). For the moment
        family the second argument is the mesh and the third the truncation
        order ``L``.
        """
        return cls.from_mesh_and_L(arr, mesh, L)


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
    :meth:`zeros_on` / :meth:`from_face_arrays` factories (all
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
                f"via {type(self).__name__}.zeros_on / "
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

    def _check_partner(self, other: Field) -> Self:
        partner = super()._check_partner(other)
        # Mesh-binding override — two boundary fields can share a space
        # (``"sn_trace"``, same total_size — TraceSpace.__eq__ is on
        # ``(name, shape)``) but differ in mesh identity.
        if self.mesh is not partner.mesh:
            raise ValueError(
                f"{type(self).__name__} arithmetic across distinct SNMesh "
                f"instances is forbidden — the field is mesh-bound."
            )
        if self.layout is not partner.layout:
            # Belt-and-suspenders: operands sourced from the same mesh
            # share the cached ``mesh.trace`` (one layout identity), so
            # this fires only for hand-built operands. Fall back to
            # structural equality — same face names + shapes + offsets.
            if self.layout != partner.layout:
                raise ValueError(
                    f"{type(self).__name__} layout mismatch — operands have "
                    f"structurally distinct FaceLayouts."
                )
        return partner

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
    def from_mesh(cls, values: NDArray, mesh: "SNMesh"):
        r"""Construct from a flat trace buffer + mesh, sourcing ``mesh.trace``.

        The boundary-locus counterpart to the bulk families'
        :meth:`AngularField.from_mesh`, giving every locus the SAME
        "construct this leaf from raw ``values`` + ``mesh``" surface. Unlike
        :meth:`from_face_arrays` (a per-face dict, packed into the flat
        layout) it takes the already-packed ``(layout.total_size,)`` buffer
        directly. This is the reconstruction path the named composition
        :meth:`~orpheus.numerics.field.Field._from_balance` uses to re-wrap a
        differenced trace buffer in the residual role — so a boundary residual
        lands on the same ``mesh.trace`` space as every other boundary leaf.

        Raises
        ------
        ValueError
            If ``mesh.trace is None`` (a trace-less 2-D cylindrical mesh,
            which has no SN sweep) — a boundary field cannot be built
            without a boundary trace.
        """
        space = mesh.trace
        if space is None:
            raise ValueError(
                f"{cls.__name__}.from_mesh: mesh has no TraceSpace "
                f"(mesh.trace is None — only trace-less 2-D cylindrical "
                f"meshes, which have no SN sweep, hit this). A boundary "
                f"field cannot be built without a boundary trace."
            )
        return cls(values=values, space=space, mesh=mesh)

    @classmethod
    def zeros_on(cls, mesh: "SNMesh"):
        r"""Construct a zero boundary field sized to ``mesh`` (B.5.A).

        Sources ``space = mesh.trace`` (the cached unified TraceSpace) and
        delegates to :meth:`~orpheus.numerics.field.Field.zeros`. The
        uniform leaf-side allocator that
        :meth:`~orpheus.transport.timed_full_field.TimedFullField.zeros`
        calls; replaces the retired ``SNMesh.zeros_boundary_flux`` factory.
        """
        space = mesh.trace
        if space is None:
            raise ValueError(
                f"{cls.__name__}.zeros_on: mesh has no TraceSpace "
                f"(mesh.trace is None — only trace-less 2-D cylindrical "
                f"meshes, which have no SN sweep, hit this). A boundary "
                f"field cannot be built without a boundary trace."
            )
        return cls.zeros(space, mesh=mesh)

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
        if layout is None:
            raise ValueError(
                f"{cls.__name__}.from_face_arrays: mesh.trace carries no "
                f"FaceLayout (a bare-constructor TraceSpace, not "
                f"TraceSpace.from_quadrature_and_layout). A boundary field "
                f"cannot be packed without a face layout."
            )
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
