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

The field vocabulary (issues #205 / #201; the #290 P2.5 axis-coherence
ruling) is a grid of THREE orthogonal axes: **locus** {Bulk,
Boundary(field) / Trace(space)} × **family** {Angular, Scalar, Moment}
× **role** {Flux, SourceSink, Residual}. (A fourth role, Displacement,
existed until campaign 1 CS3, 2026-08-19 — flux lives in V now, so
differences are same-typed and the sibling family retired.) A bulk leaf is
named ``<Family><Role>``, a boundary leaf ``<Family>Boundary<Role>`` —
"Boundary" is the locus qualifier, never a fourth family. This module
provides the *locus + family* axes as ABCs; the *role* leaves
(``AngularFlux``, ``AngularSourceSink``, ...) sit beneath them::

    Field (numerics, L1 — values + space + dunder algebra)
     ├─ BulkField (ABC)           codim-0 (cell centres): mesh-binding + ng + _phase_space_shape
     │   ├─ AngularField (ABC)    + N + (N,ng,*spatial) from_mesh/_ndarray, parametrized by _SPACE_NAME
     │   │   ├─ AngularFlux           role leaf  (flux)
     │   │   └─ AngularSourceSink     role leaf  (source; renamed from PerOrdinateSource in B.2)
     │   ├─ ScalarField (ABC)     + (ng,*spatial) from_mesh/_ndarray, parametrized by _SPACE_NAME
     │   │   ├─ ScalarFlux            role leaf  (flux)
     │   │   └─ ScalarSourceSink       role leaf  (source; renamed from IsotropicSource in B.2)
     │   └─ MomentField (ABC)     family marker; the moment shape is leaf-specific
     │       └─ HarmonicMomentFlux   role leaf  (flux-only for now)
     └─ FaceField[K] (ABC)        codim-1 (faces/edges): flat single-buffer + FaceLayout[K]
         │                        slice-views + mesh/layout guards + from_mesh/zeros_on via
         │                        _face_space_of. STRUCTURE only — the metric descends PER LEAF
         │                        (spatial |Ω·n̂|·w; pole V_cell), never on this ABC (ERR-067).
         ├─ BoundaryField (ABC, FaceField[str])   SPATIAL faces (keyed by name) + from_face_arrays;
         │   │                    the FullField boundary-slot discriminator (the pole is NOT one)
         │   ├─ AngularBoundaryField (ABC)   mesh: SNMesh + AngularTraceSpace (mesh.angular_trace)
         │   │   ├─ AngularBoundaryFlux          role leaf  (flux)
         │   │   ├─ AngularBoundarySourceSink    role leaf  (source; B.3 — orpheus.transport.source_sinks)
         │   │   ├─ AngularBoundaryResidual      role leaf  (residual; B.3 — orpheus.transport.residuals)
         │   └─ ScalarBoundaryField (ABC)    ScalarTraceSpace (DiffusionMesh.scalar_trace; #290 P2/P7a)
         │       ├─ ScalarBoundaryFlux           role leaf  (flux — the per-face (J⁺, J⁻) pair)
         ├─ RadialCharacteristicInteriorField (ABC, FaceField[(level,sign)])  ANGULAR edge — the ψ½
         │   │                    marched cells (μ = μ_start; #282 route (a); System B's interior);
         │   │                    a FaceField SIBLING of BoundaryField, never a child
         │   ├─ RadialCharacteristicInteriorFlux          role leaf  (flux — the ψ½ state)
         │   ├─ RadialCharacteristicInteriorSourceSink    role leaf  (source — q½ cells)
         │   ├─ RadialCharacteristicInteriorResidual      role leaf  (residual)
         └─ RadialCharacteristicBoundaryField (ABC, FaceField[(level,sign)])  the r = R ψ½ corner
             │                    (System B's boundary; the unified 3-tuple base retired at 4e)
             ├─ RadialCharacteristicBoundaryFlux          role leaf  (flux — corner data/defect)
             ├─ RadialCharacteristicBoundarySourceSink    role leaf  (source — corner datum)
             ├─ RadialCharacteristicBoundaryResidual      role leaf  (residual)

Parametrization (no twin paths)
===============================

The per-family phase-space shape is the one abstract hook
:meth:`BulkField._phase_space_shape`, used by the shared
``__post_init__`` validator. The Angular/Scalar families additionally
expose a ``from_mesh`` classmethod parametrized by the leaf's
``_SPACE_NAME`` :class:`~typing.ClassVar` (so ``AngularFlux``'s space is
named ``"angular_flux"`` and ``AngularSourceSink``'s ``"angular_source_sink"``,
preserving the pre-B.1 space identities bit-for-bit). ``MomentField`` and the
``BoundaryField`` families build their spaces differently (a
TensorProductSpace keyed on ``L``; the mesh's cached trace —
``mesh.angular_trace`` / ``mesh.scalar_trace`` via
:meth:`FaceField._face_space_of`) and so do not use
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
from collections.abc import Hashable
from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar, Generic, Mapping, Self, TypeVar

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.field import Field
from orpheus.numerics.moment_layout import cell_moment_count
from orpheus.numerics.space import FunctionSpace
from orpheus.numerics.spaces.spatial_moment_space import (
    SpatialMomentSpace,
    spatial_moment_tail,
)
from orpheus.numerics.spaces.spherical_harmonic_space import (
    SphericalHarmonicSpace,
)
from orpheus.numerics.spaces.scalar_trace_space import ScalarTraceSpace
from orpheus.numerics.spaces.angular_trace_space import AngularTraceSpace
from orpheus.numerics.spaces.radial_characteristic_space import (
    RadialCharacteristicBoundarySpace,
    RadialCharacteristicInteriorSpace,
)

if TYPE_CHECKING:
    from orpheus.diffusion.augmented_mesh import DiffusionMesh
    from orpheus.numerics.face_layout import FaceLayout
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.transport.mesh.material_mesh import MaterialMesh


#: The face-key type of a :class:`FaceField` and its
#: :class:`~orpheus.numerics.face_layout.FaceLayout`: ``str`` face names for
#: the spatial :class:`BoundaryField`, the ``(level, sign)`` tuple for
#: the ψ½ split loci. Bounded by
#: :class:`~collections.abc.Hashable` — the face key is a mapping key.
K = TypeVar("K", bound=Hashable)


__all__ = [
    "BulkField",
    "AngularField",
    "ScalarField",
    "MomentField",
    "FaceField",
    "BoundaryField",
    "AngularBoundaryField",
    "ScalarBoundaryField",
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
    ``mesh.quad`` access is confined to the angular family, whose
    :class:`AngularField` / :class:`MomentField` / :class:`AngularBoundaryField`
    declarations covariantly NARROW ``mesh`` back to ``SNMesh`` (no casts
    — the narrowing lives in the field declarations, #267).
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
        the construct-general / select-narrow discipline: the CALLER
        selects. Since S3-A landed, the iterate / cell-emit / source seams
        (the SI cold starts, ``coupled_system``, ``windowing``) thread the
        scheme's ``spatial_basis_per_axis`` here explicitly, so LD
        production fields DO carry the axis while DD/Step (per_axis == 1)
        get no factor. Reading the scheme by DEFAULT here would still be
        wrong: the widening must remain the decision of the seams that
        also FILL the axis — a Pattern-4 concern (an axis no producer
        fills is an illegal state).
        """
        n_moments = cell_moment_count(spatial_moments_per_axis, mesh.ndim)
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

    The storage base for the angular role leaves (``AngularFlux``,
    ``AngularSourceSink``, ``AngularResidual``).
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

    def _integrate_angular_values(self) -> "NDArray":
        r"""The ONE moment-0 reduction body :math:`\sum_n w_n\,(\cdot)_n`.

        Contracts the leading ``N`` axis with the quadrature weights
        ``w_n`` — ``(N, ng, nx, ny[, 2^d]) → (ng, nx, ny[, 2^d])``. The
        ``ng...`` einsum is spatial-moment-axis-agnostic, so a
        φ̂-carrying field reduces to a φ̂-carrying scalar (#240 D5b-S3);
        the moment width propagates as a TYPED factor read off this
        field's space, not an opaque axis.

        Role leaves wrap this body in their own scalar type
        (:meth:`AngularFlux.integrate_angular
        <orpheus.transport.fields.angular_flux.AngularFlux.integrate_angular>`
        → ``ScalarFlux``; until campaign 1 CS3 the displacement sibling
        wrapped the same body — a linear reduction is its own tangent
        map). Values-level, role-blind, single source of truth for the
        canonical angular reduction (the DSA restriction ``R`` rides it,
        #2).
        """
        return np.einsum(
            "n,ng...->g...", self.mesh.quad.weights, self.values,
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

    The storage base for the scalar role leaves (``ScalarFlux``,
    ``ScalarSourceSink``, ``ScalarResidual``).
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
    (the moment-space flux state) and
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
        n_moments = cell_moment_count(self.spatial_moments, self.mesh.ndim)
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
        n_moments = cell_moment_count(spatial_moments, mesh.ndim)
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
# Face locus (codim-1) — the shared flat-buffer discipline
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class FaceField(Field, Generic[K]):
    r"""Codim-1 face storage base — a mesh-bound flat-buffer :class:`Field`
    on a layout-bearing face space (method- and locus-agnostic).

    The single parent of every **codim-1** transport field — the faces and
    edges bounding a phase-space cell, as opposed to :class:`BulkField`'s
    codim-0 cell centres. It owns, ONCE, the flat-buffer discipline the two
    codim-1 loci share: the single-buffer storage contract
    (``values.shape == (layout.total_size,)``), per-face slice-view access
    keyed by the layout key ``K`` (a ``str`` face name for the spatial trace,
    a ``(level, sign, part)`` tuple for the ψ½ pole edge), the cross-mesh /
    cross-layout arithmetic guards, the read-through :attr:`layout` property,
    and the :meth:`from_mesh` / :meth:`zeros_on` factories (via the single
    :meth:`_face_space_of` hook, so they build the concrete subclass on ITS
    mesh's cached face space).

    **STRUCTURE only — no metric on this ABC.** Exactly as :class:`BulkField`
    carries no metric (each bulk leaf's ``V·w`` lives on its own space), a
    face field's Hilbert metric descends **per leaf**, on the leaf's face
    space: the spatial trace carries the partial-current ``|Ω·n̂|·w`` metric
    (:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`),
    the ψ½ pole the SPD radial cell-volume STATE metric :math:`V_{\rm cell}`
    (the split ψ½ spaces in :mod:`orpheus.numerics.spaces.radial_characteristic_space`).
    The through-flux coefficient is NOT the state metric (ERR-067, vv Mode 12)
    — so this ABC deliberately unifies the *layout*, never the *measure*.

    Two codim-1 loci realize the base as **siblings** (NOT parent/child — the
    :class:`~orpheus.transport.full_field.FullField` composite discriminates
    its boundary slot by ``isinstance(·, BoundaryField)``, a test the ψ½ pole
    must FAIL):

    * :class:`BoundaryField` — the **spatial** faces (``FaceField[str]``): the
      boundary of the SPATIAL domain, keyed by face name.
    * the ψ½ split loci (:class:`RadialCharacteristicInteriorField` /
      :class:`RadialCharacteristicBoundaryField`) — the **angular** edge
      (``FaceField[tuple[int, int, str]]``): the ψ½ pole seed, the boundary of
      the ANGULAR domain (:math:`\mu = \mu_{\rm start}`), keyed by
      ``(level, sign, part)``.

    The base ``mesh`` annotation is the method-agnostic
    :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`; each family
    covariantly NARROWS it to its method-mesh (:class:`SNMesh` /
    :class:`DiffusionMesh`) — the #267 ``BulkField → AngularField``
    discipline. Storage is a SINGLE flat backing buffer; per-face access is
    slice-view, no copies. Inherits Field's same-class/same-space dunder
    algebra. Abstract — instantiate a concrete role leaf of one of the loci.
    """

    mesh: "MaterialMesh"

    # ── The per-family face-space source ─────────────────────────────

    @classmethod
    @abstractmethod
    def _face_space_of(cls, mesh: "MaterialMesh") -> FunctionSpace:
        r"""Return ``mesh``'s cached face space for this family.

        The single hook the factories construct through: the spatial
        families read ``mesh.angular_trace`` (raising on the trace-less 2-D
        cylindrical mesh) / ``mesh.scalar_trace`` (a :class:`DiffusionMesh`
        member, raising on a bare MaterialMesh — #290 P7a); the
        starting-direction family reads the split ψ½ mesh spaces
        (R12a-keyed). MUST return a layout-bearing space or raise
        :class:`ValueError` with the family's own diagnosis.
        """
        raise NotImplementedError

    # ── Construction validation ──────────────────────────────────────

    def __post_init__(self) -> None:
        super().__post_init__()  # Field: values.shape == space.shape.
        # The space IS the face space and carries the FaceLayout
        # (illegal-states-unrepresentable): a face field on a layout-less
        # space cannot do face_view. Families ADD their space-type narrowing
        # (AngularTraceSpace / ScalarTraceSpace / the split ψ½ spaces) on
        # top of this structural floor.
        layout = getattr(self.space, "layout", None)
        if layout is None:
            raise TypeError(
                f"{type(self).__name__} requires a layout-bearing face "
                f"space (a trace / starting-direction space built via its "
                f"factory); got space={self.space!r}. Build via "
                f"{type(self).__name__}.zeros_on / from_mesh."
            )
        expected = (layout.total_size,)
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
        # Mesh-binding override — two face fields can share a space (same
        # name, same total_size — space __eq__ is on ``(name, shape)``) but
        # differ in mesh identity. The message names the operand's own mesh
        # type (SNMesh for the angular families, DiffusionMesh for the scalar
        # trace).
        if self.mesh is not partner.mesh:
            raise ValueError(
                f"{type(self).__name__} arithmetic across distinct "
                f"{type(self.mesh).__name__} instances is forbidden — the "
                f"field is mesh-bound."
            )
        if self.layout is not partner.layout:
            # Belt-and-suspenders: operands sourced from the same mesh share
            # the cached space (one layout identity), so this fires only for
            # hand-built operands. Fall back to structural equality — same
            # keys + shapes + offsets.
            if self.layout != partner.layout:
                raise ValueError(
                    f"{type(self).__name__} layout mismatch — operands have "
                    f"structurally distinct FaceLayouts."
                )
        return partner

    # ── Per-face access (slice views into the flat buffer) ───────────

    @property
    def layout(self) -> "FaceLayout[K]":
        r"""The per-geometry :class:`FaceLayout`, read off the space.

        The layout lives on the face space (``space.layout``), not as a
        separate field attribute. This read-through property preserves the
        ``field.layout`` access surface (``field.layout.faces``,
        ``.total_size``).
        """
        return self.space.layout  # type: ignore[attr-defined]

    def face_view(self, key: K) -> NDArray:
        r"""Return a per-face slice view into the flat backing buffer.

        The returned ndarray shares memory with :attr:`values`.

        Raises
        ------
        KeyError
            If ``key`` is not a face in this layout.
        """
        if key not in self.layout.faces:
            raise KeyError(
                f"{type(self).__name__}: no face keyed {key!r} in layout; "
                f"available: {list(self.layout.faces)!r}"
            )
        return self.layout.faces[key].slice_view(self.values)

    @property
    def face_views(self) -> Mapping[K, NDArray]:
        r"""Mapping ``{face_key: slice_view}`` for every face in the layout.

        All views memory-shared with :attr:`values`.
        """
        return {key: self.face_view(key) for key in self.layout.faces}

    # ── Construction factories (via ``cls`` — build the subclass) ────

    @classmethod
    def from_mesh(cls, values: NDArray, mesh: "MaterialMesh"):
        r"""Construct from a flat face buffer + mesh, sourcing the family's
        cached face space (:meth:`_face_space_of`).

        The codim-1 counterpart to the bulk families' ``from_mesh``, giving
        every locus the SAME "construct this leaf from raw ``values`` +
        ``mesh``" surface. It takes the already-packed ``(layout.total_size,)``
        buffer directly. This is the reconstruction path the named
        composition :meth:`~orpheus.numerics.field.Field._from_balance` uses
        to re-wrap a differenced face buffer in the residual role — so a face
        residual lands on the same cached face space as every other leaf of
        its family.
        """
        return cls(values=values, space=cls._face_space_of(mesh), mesh=mesh)

    @classmethod
    def zeros_on(cls, mesh: "MaterialMesh"):
        r"""Construct a zero face field sized to ``mesh`` (B.5.A).

        Sources the space from :meth:`_face_space_of` (the family's cached
        face space) and delegates to
        :meth:`~orpheus.numerics.field.Field.zeros` — the uniform leaf-side
        allocator.
        """
        return cls.zeros(cls._face_space_of(mesh), mesh=mesh)


# ═══════════════════════════════════════════════════════════════════════
# Boundary locus (spatial faces)
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class BoundaryField(FaceField[str]):
    r"""Spatial-face storage base — the SPATIAL locus of :class:`FaceField`.

    The boundary of the SPATIAL domain: codim-1 faces keyed by name
    (``"xmin"`` / ``"xmax"`` / ...), under the partial-current
    ``|Ω·n̂|·w`` metric. The flat-buffer discipline (storage, slice views,
    mesh/layout guards, :meth:`~FaceField.from_mesh` /
    :meth:`~FaceField.zeros_on`) is inherited from :class:`FaceField`; this
    intermediate (a) adds the spatial-only :meth:`from_face_arrays` per-face
    packer, and (b) is the type the
    :class:`~orpheus.transport.full_field.FullField` composite discriminates
    its boundary slot by — the ψ½ pole (the RadialCharacteristic split loci, a
    :class:`FaceField` SIBLING) is deliberately NOT a ``BoundaryField``.

    Two storage families realize the SPATIAL locus (#290 P2; named per the
    P2.5 axis-coherence ruling — family-qualified, uniform role tokens):

    * :class:`AngularBoundaryField` — the ANGULAR family (``mesh`` narrowed to
      :class:`SNMesh`, space to the quadrature-coupled
      :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`).
    * :class:`ScalarBoundaryField` — the SCALAR family (``mesh`` narrowed to
      :class:`DiffusionMesh`, space to
      :class:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace`;
      per-face ``(J⁺, J⁻)`` partial-current pairs).

    A boundary trace is method BEHAVIOR (#290 P7a): the bare MaterialMesh
    data carrier owns no trace, so a boundary field on one is
    unrepresentable (:meth:`~FaceField._face_space_of` raises). Abstract —
    instantiate a concrete role leaf of one of the families.
    """

    # ── The spatial-only per-face packer (str-keyed) ─────────────────

    @classmethod
    def from_face_arrays(
        cls, mesh: "MaterialMesh", face_arrays: Mapping[str, NDArray],
    ):
        r"""Construct from per-face ndarrays, packing into the flat layout.

        ``face_arrays`` must cover EVERY face in the mesh's layout; each
        per-face ndarray's shape must match the layout slot's shape. Spatial
        faces only — the ψ½ pole builds through its role factories, not a
        per-face dict.

        Raises
        ------
        ValueError
            If ``face_arrays`` keys differ from the mesh's layout faces,
            or any per-face ndarray's shape mismatches the layout slot.
        """
        space = cls._face_space_of(mesh)
        layout = getattr(space, "layout", None)
        if layout is None:
            raise ValueError(
                f"{cls.__name__}.from_face_arrays: the trace space carries "
                f"no FaceLayout (a bare-constructor space, not one built "
                f"via its factory). A trace field cannot be packed "
                f"without a face layout."
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


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class AngularBoundaryField(BoundaryField):
    r"""Angular boundary-trace storage base — the SN family of the
    :class:`BoundaryField` locus.

    Carries what is ANGULAR about the locus: the :class:`SNMesh`
    binding (covariant narrowing of the base ``mesh: MaterialMesh``)
    and the space narrowing to the unified quadrature-coupled
    :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` (the
    ``|Ω·n̂|⊙w_n`` partial-current metric + the ``omega_dot_n``
    selector table). All storage/factory/guard machinery is inherited
    from :class:`BoundaryField`. The concrete role leaves are
    ``AngularBoundaryFlux``, ``AngularBoundarySourceSink``, ``AngularBoundaryResidual``.
    Abstract — instantiate a role leaf.
    """

    mesh: "SNMesh"
    # The static twin of the __post_init__ isinstance gate below (the
    # ``mesh: SNMesh`` covariant-narrowing idiom): an angular boundary
    # field's space IS the quadrature-coupled AngularTraceSpace, so
    # consumers of its atoms (``omega_dot_n``, the face layout, the
    # partial-current metric) type-check without re-narrowing.
    space: AngularTraceSpace

    # ── Construction validation (angular narrowing) ──────────────────

    def __post_init__(self) -> None:
        # A.5: the space IS the unified angular AngularTraceSpace. The
        # isinstance narrowing runs FIRST so the family-specific
        # message fires for the common misuse (a bare FunctionSpace or
        # a scalar trace passed where the angular trace belongs).
        if not isinstance(self.space, AngularTraceSpace):
            raise TypeError(
                f"{type(self).__name__} requires an AngularTraceSpace carrying a "
                f"FaceLayout (A.5 re-home); got space={self.space!r}. Build "
                f"via {type(self).__name__}.zeros_on / "
                f"from_face_arrays, or pass mesh.angular_trace as the space."
            )
        super().__post_init__()

    # ── The angular trace-space source (``mesh.angular_trace``) ──────────────

    @classmethod
    def _face_space_of(cls, mesh: "MaterialMesh") -> AngularTraceSpace:
        space = getattr(mesh, "angular_trace", None)
        if space is None:
            raise ValueError(
                f"{cls.__name__}: mesh has no AngularTraceSpace (mesh.angular_trace is "
                f"None — only trace-less 2-D cylindrical meshes, which "
                f"have no SN sweep, hit this). A boundary field cannot "
                f"be built without a boundary trace."
            )
        return space


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class ScalarBoundaryField(BoundaryField):
    r"""Scalar boundary-trace storage base — the ``(J⁺, J⁻)`` family of
    the :class:`BoundaryField` locus (#290 P2).

    Carries what is SCALAR about the locus: the space narrowing to
    :class:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace`
    (per-face partial-current pairs under the face-AREA metric), the
    :class:`DiffusionMesh` binding (covariant narrowing of the base
    ``mesh: MaterialMesh`` — the exact :class:`AngularBoundaryField` /
    ``SNMesh`` discipline; #290 P7a), and the trace-space source
    ``mesh.scalar_trace``. A scalar trace lives on the DIFFUSION phase
    space — when DSA (#2) restricts an SN solve, the SN mesh promotes
    (``DiffusionMesh.from_material_mesh(sn_mesh)`` — an ``SNMesh`` IS a
    ``MaterialMesh``) and :math:`A_{\rm diff}`'s fields bind to the
    promoted mesh. The concrete role leaves are
    :class:`~orpheus.transport.fields.scalar_boundary_flux.ScalarBoundaryFlux`
    (the state — and, since campaign 1 CS3, its own iterate increments:
    differences are same-typed); source/residual siblings join when their
    operator codomains demand them (#290 P4). Abstract — instantiate a
    role leaf.
    """

    mesh: "DiffusionMesh"

    # ── Construction validation (scalar narrowing) ────────────────────

    def __post_init__(self) -> None:
        # The family's space narrowing runs FIRST so the family-specific
        # message fires for the common misuse (an angular AngularTraceSpace or
        # a bare FunctionSpace passed where the scalar trace belongs).
        if not isinstance(self.space, ScalarTraceSpace):
            raise TypeError(
                f"{type(self).__name__} requires a ScalarTraceSpace (the "
                f"(J⁺, J⁻) partial-current trace); got space="
                f"{self.space!r}. Build via {type(self).__name__}.zeros_on "
                f"/ from_face_arrays, or pass mesh.scalar_trace as the "
                f"space."
            )
        super().__post_init__()

    # ── The scalar trace-space source (``mesh.scalar_trace``) ─────────

    @classmethod
    def _face_space_of(cls, mesh: "MaterialMesh") -> ScalarTraceSpace:
        space = getattr(mesh, "scalar_trace", None)
        if space is None:
            raise ValueError(
                f"{cls.__name__}: mesh carries no scalar trace "
                f"(mesh.scalar_trace) — a bare MaterialMesh is the "
                f"method-agnostic DATA carrier; the (J⁺, J⁻) boundary "
                f"trace is diffusion method BEHAVIOR (#290 P7a). Build "
                f"the field on a DiffusionMesh "
                f"(DiffusionMesh.from_material_mesh promotes)."
            )
        return space


# ═══════════════════════════════════════════════════════════════════════
# The SPLIT ψ½ loci — System B's interior ⊕ boundary (Phase B)
# ═══════════════════════════════════════════════════════════════════════
#
# The coupled-block campaign poses the ψ½ ray as System B — its OWN
# interior ⊕ boundary composite — split into two INDEPENDENT loci, exactly
# as the spatial domain splits into BulkField (interior) and BoundaryField
# (boundary). They are FaceField SIBLINGS, not parent/child of each other,
# and each keys by (level, sign) on its own split space. The concrete role
# leaves (…Flux state / …SourceSink emission) live in fields/ and
# source_sinks/ (the …Displacement increment family retired at CS3). The historical
# UNIFIED base (cells ⊕ corner interleaved on one
# FaceField[(level, sign, part)] buffer) retired at 4e, when the fused
# (L+C) walk went split-native — System B's composite, which took the
# freed ``RadialCharacteristicField`` name at 4e-e1b, is now the ONLY ψ½
# representation.


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class RadialCharacteristicInteriorField(FaceField[tuple[int, int]]):
    r"""System B's INTERIOR locus — the ψ½ cells (Phase B split).

    The ``(ng, nx)`` half-angle flux at every radial cell, per seed-carrying
    ``(level, sign)`` leg — the marched interior state that
    :class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator`
    (A_BB) evolves (μ∂_r + σ_t), that A_AB reads (inward leg) and A_BA writes.
    The interior half of the ψ½ split; a ``FaceField[tuple[int, int]]`` keyed by
    ``(level, sign)`` on the
    :class:`~orpheus.numerics.spaces.radial_characteristic_space.RadialCharacteristicInteriorSpace`
    (SPD ``G_sd = V_cell`` state metric). A :class:`FaceField` **sibling** of the
    boundary locus :class:`RadialCharacteristicBoundaryField`, NOT a child (like
    :class:`BulkField` vs :class:`BoundaryField`). ``mesh`` is :class:`SNMesh` and
    the space source is the R12a-keyed ``mesh.radial_characteristic_interior_space``
    — construction on a non-carrying mesh is unrepresentable (the factory raises;
    the composite spells absence as ``None``). The concrete role leaves are
    ``RadialCharacteristicInteriorFlux`` (state — same-class signed differences
    carry the iterate increment since campaign 1 CS3, 2026-08-19, when the
    displacement sibling retired) and
    ``RadialCharacteristicInteriorSourceSink`` (emission). Abstract
    — instantiate a role leaf.
    """

    mesh: "SNMesh"
    space: RadialCharacteristicInteriorSpace

    def __post_init__(self) -> None:
        if not isinstance(self.space, RadialCharacteristicInteriorSpace):
            raise TypeError(
                f"{type(self).__name__} requires a "
                f"RadialCharacteristicInteriorSpace (read off "
                f"mesh.radial_characteristic_interior_space); got "
                f"space={self.space!r}."
            )
        super().__post_init__()  # FaceField: Field shape + layout-bearing check.

    @classmethod
    def _face_space_of(
        cls, mesh: "MaterialMesh",
    ) -> RadialCharacteristicInteriorSpace:
        r"""Return ``mesh``'s cached ψ½ interior space, or raise (R12a)."""
        space = getattr(mesh, "radial_characteristic_interior_space", None)
        if space is None:
            raise ValueError(
                f"{cls.__name__}: mesh carries no "
                f"radial_characteristic_interior_space — no μ-level consumes "
                f"independent starting-direction state (R12a; Cartesian and the "
                f"production cylinder rules land here). System B's interior block "
                f"is spelled absent, never a zero-DOF field."
            )
        return space

    @property
    def levels(self) -> tuple[int, ...]:
        r"""The seed-carrying μ-level indices (read off the space)."""
        return self.space.levels

    def cells(self, level: int, sign: int) -> NDArray:
        r"""The ``(ng, nx)`` ψ½ cells view for ``(level, sign)`` — memory-shared."""
        return self.space.slot_view(self.values, level, sign)


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class RadialCharacteristicBoundaryField(FaceField[tuple[int, int]]):
    r"""System B's BOUNDARY locus — the ψ½ r = R corner (Phase B split).

    The ``(ng,)`` half-angle flux at the outer radius r = R, per seed-carrying
    ``(level, sign)`` leg — System B's BC locus, on which
    :class:`~orpheus.sn.operators.boundary.RadialCharacteristicBoundaryOperator`
    (B_b) acts. Inflow corner (``sign = -1``) is the given BC data; outflow corner
    (``sign = +1``) is the defect row (ruling R13). The boundary half of the split
    of the ψ½ split; a ``FaceField[tuple[int, int]]`` keyed
    by ``(level, sign)`` on the
    :class:`~orpheus.numerics.spaces.radial_characteristic_space.RadialCharacteristicBoundarySpace`
    (``G = V(r = R)`` corner gauge). A :class:`FaceField` **sibling** of the
    interior locus :class:`RadialCharacteristicInteriorField`, NOT a child.
    ``mesh`` is :class:`SNMesh` and the space source is the R12a-keyed
    ``mesh.radial_characteristic_boundary_space``. The concrete role leaves are
    ``RadialCharacteristicBoundaryFlux`` (state — same-class signed differences
    carry the iterate increment since campaign 1 CS3, 2026-08-19, when the
    displacement sibling retired) and
    ``RadialCharacteristicBoundarySourceSink`` (emission). Abstract
    — instantiate a role leaf.
    """

    mesh: "SNMesh"
    space: RadialCharacteristicBoundarySpace

    def __post_init__(self) -> None:
        if not isinstance(self.space, RadialCharacteristicBoundarySpace):
            raise TypeError(
                f"{type(self).__name__} requires a "
                f"RadialCharacteristicBoundarySpace (read off "
                f"mesh.radial_characteristic_boundary_space); got "
                f"space={self.space!r}."
            )
        super().__post_init__()  # FaceField: Field shape + layout-bearing check.

    @classmethod
    def _face_space_of(
        cls, mesh: "MaterialMesh",
    ) -> RadialCharacteristicBoundarySpace:
        r"""Return ``mesh``'s cached ψ½ boundary space, or raise (R12a)."""
        space = getattr(mesh, "radial_characteristic_boundary_space", None)
        if space is None:
            raise ValueError(
                f"{cls.__name__}: mesh carries no "
                f"radial_characteristic_boundary_space — no μ-level consumes "
                f"independent starting-direction state (R12a). System B's "
                f"boundary block is spelled absent, never a zero-DOF field."
            )
        return space

    @property
    def levels(self) -> tuple[int, ...]:
        r"""The seed-carrying μ-level indices (read off the space)."""
        return self.space.levels

    def corner(self, level: int, sign: int) -> NDArray:
        r"""The ``(ng,)`` r = R corner view for ``(level, sign)`` — memory-shared.

        Inflow corner (``sign = -1``): the given-data / identity row;
        outflow corner (``sign = +1``): the defect row (ruling R13).
        """
        return self.space.slot_view(self.values, level, sign)

