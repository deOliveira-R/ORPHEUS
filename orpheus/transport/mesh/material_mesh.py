r"""Mesh + materials — the method-agnostic transport data carrier.

:class:`MaterialMesh` is the "mesh + materials" middle type the codebase
was missing.  Between the geometry :class:`~orpheus.geometry.mesh.Mesh1D`
/ :class:`~orpheus.geometry.mesh.Mesh2D` (which carry material *ids* but
no :class:`~orpheus.data.macro_xs.mixture.Mixture` cross sections) and a
method-specific phase space such as
:class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` (mesh + materials + *quadrature* +
sweep machinery) there was no carrier for *just* mesh + materials.

The abstraction axis is **data vs behavior**:

* :class:`MaterialMesh` is the **method-agnostic transport state / data**
  — the axes/mesh, the per-material :class:`Mixture` dict, the material
  map, cell volumes, the uniform group count :attr:`ng`, the natural
  volume measure, and the macroscopic XS field built from them.
* The **method layer** (angular quadrature + sweep/streaming stencil +
  boundary trace + closures) is *behavior*.  A method-specific mesh
  **is a** :class:`MaterialMesh` that adds that behavior:
  ``SNMesh(MaterialMesh)`` for the data, conforming structurally to a
  future ``TransportMethod`` Protocol for the behavior (minted when a
  2nd method-mesh — ``MoCMesh`` / ``CPMesh`` — exists to share the
  contract; deferred per the defer-until-2-instances rule).

This is the layer where cross-section **homogenization** lands: a
fine-mesh :class:`~orpheus.sn.solution.Solution` plus a coarse
:class:`~orpheus.geometry.mesh.Mesh` produce a homogenized
:class:`MaterialMesh` (flux·volume-weighted collapse), which a transport
method can then *promote* back to a solvable phase space
(:meth:`SNMesh.from_material_mesh`).

Layer (``tests/test_layer_imports.py``): L2 ``transport``.  It imports
only ``geometry`` (legacy mesh shapes), ``numerics`` (the volume
measure), ``data`` (the :class:`Mixture` type, ``TYPE_CHECKING``), and
its sibling :mod:`~orpheus.transport.mesh.axis` /
:mod:`~orpheus.transport.mesh.material_xs_field` modules.  It imports no
L3 method package — which is exactly what let it be promoted out of
``orpheus.sn``.
"""

from __future__ import annotations

from functools import reduce
from typing import TYPE_CHECKING

import numpy as np

from orpheus.geometry import Mesh1D, Mesh2D
from orpheus.transport.mesh.axis import (
    Axis1D,
    AxisMesh,
    axes_from_legacy_mesh,
    coord_system as _axis_coord_system,
    spatial_shape as _axis_spatial_shape,
)

if TYPE_CHECKING:
    from orpheus.data.macro_xs.mixture import Mixture
    from orpheus.geometry import CoordSystem
    from orpheus.transport.mesh.material_xs_field import MaterialXSField


__all__ = ["InconsistentMaterialsError", "MaterialMesh"]


class InconsistentMaterialsError(ValueError):
    """Raised when a materials dict has inconsistent metadata.

    Currently triggered when materials disagree on ``ng`` (number of
    energy groups).  A :class:`MaterialMesh` requires a uniform group
    structure across all materials in its ``mat_map`` because every
    transport operator that consumes the mesh assumes one ``ng``.  A
    homogenization / energy-condensation step must precede method-mesh
    construction if the input materials carry different group
    structures.
    """


class MaterialMesh:
    r"""Method-agnostic mesh + materials carrier.

    Axis-primary (C5.1, #225): the canonical spatial representation is
    :attr:`axes` — a tuple of :class:`~orpheus.transport.mesh.axis.Axis1D`
    — from which all shape metadata derives.  Constructed either from a
    legacy :class:`~orpheus.geometry.mesh.Mesh1D` /
    :class:`~orpheus.geometry.mesh.Mesh2D` (converted to axes once at the
    inbound boundary; the legacy object is retained as :attr:`mesh` for
    the consumers still reading through it) or axis-natively via
    :meth:`from_axes`.

    Parameters
    ----------
    mesh : Mesh1D or Mesh2D
        Base geometry.  Its material assignment (``mat_ids`` on
        :class:`Mesh1D`, ``mat_map`` on :class:`Mesh2D`) keys into
        ``materials``.
    materials : dict mapping material id to Mixture
        Macroscopic cross sections keyed by the integer ids appearing in
        the mesh's material map.  The authoritative source of truth for
        both cross sections and the group count :attr:`ng`.  All
        materials must agree on ``ng`` — heterogeneous group structures
        are a homogenization-step concern that must precede method-mesh
        construction.

    Attributes
    ----------
    mesh : Mesh1D or Mesh2D or None
        Inbound provenance / legacy adapter (``None`` for axis-native
        d≥3 meshes).
    axes : tuple of Axis1D
        Canonical spatial representation.
    materials : dict mapping material id to Mixture
        The materials dict passed at construction (single source of
        truth).
    mat_map : np.ndarray
        Material-id assignment, shape :attr:`spatial_shape`.
    ng : int
        Number of energy groups, derived from materials and validated
        for consistency.
    """

    def __init__(
        self,
        mesh: Mesh1D | Mesh2D,
        materials: "dict[int, Mixture]",
    ) -> None:
        # Legacy inbound surface (C5.1 axis-primary inversion, #225):
        # convert the Mesh1D / Mesh2D declaration to the canonical axis
        # tuple ONCE at the boundary, extract the material assignment the
        # axes cannot carry (``mat_ids`` on Mesh1D, ``mat_map`` on
        # Mesh2D), and run the same data-construction body as
        # :meth:`from_axes`.
        self._init_data(
            axes=axes_from_legacy_mesh(mesh),
            mesh=mesh,
            mat_map=mesh.mat_ids if isinstance(mesh, Mesh1D) else mesh.mat_map,
            materials=materials,
        )

    def _init_data(
        self,
        *,
        axes: tuple[Axis1D, ...],
        mesh: Mesh1D | Mesh2D | None,
        mat_map: np.ndarray | None,
        materials: "dict[int, Mixture]",
    ) -> None:
        r"""The ONE data-construction body both surfaces funnel into.

        Subclasses (``SNMesh``) call this from their own ``_init_core``
        to populate the method-agnostic data block, then layer their
        behavior (quadrature, sweep stencil, boundary trace) on top.
        Every line here is bit-for-bit the data block formerly inlined
        in ``SNMesh._init_core`` (C5.1) — the split is a pure
        relocation, not a semantic change.
        """
        # ``materials`` is REQUIRED: a MaterialMesh without materials has
        # no ``ng`` and no XS field — an illegal state (coding-elegance
        # Pattern 4).  The single source of truth for the material dict
        # lives here; every operator reads materials + ``ng`` from the
        # mesh, not from a parallel argument.
        self.mesh = mesh
        self.materials: "dict[int, Mixture]" = materials
        # The axis tuple is the PRIMARY representation (C5.1): stored
        # verbatim — never round-tripped through a legacy mesh and
        # re-derived — it is the canonical dim-agnostic ground truth for
        # spatial_shape / coord / axis_widths.
        self.axes: tuple[Axis1D, ...] = tuple(axes)
        # ``np.diff(edges)`` is bitwise identical to the legacy spellings
        # it replaces (``Mesh1D.widths`` / ``Mesh2D.dx`` / ``Mesh2D.dy``).
        # ``axis_widths`` is THE single spelling of per-axis cell widths,
        # positional-by-axis.
        self.axis_widths: tuple[np.ndarray, ...] = tuple(
            np.diff(ax.edges) for ax in self.axes
        )
        # Material assignment: the one construction payload the axes do
        # not carry. ``None`` (axis-native default) → single material
        # with id 0; shape MUST match spatial_shape (parse, don't
        # validate downstream).
        if mat_map is None:
            mat_map = np.zeros(self.spatial_shape, dtype=int)
        else:
            mat_map = np.asarray(mat_map, dtype=int)
            if mat_map.shape != self.spatial_shape:
                raise ValueError(
                    f"MaterialMesh: mat_map shape {mat_map.shape} must "
                    f"match spatial_shape={self.spatial_shape}"
                )
        self.mat_map: np.ndarray = mat_map
        # Cell volumes / radial face areas stay dataclass-owned while the
        # adapter is present (preserves the Mesh1D curvilinear formulas +
        # the ``precomputed_volumes`` ULP escape hatch bit-identically).
        # Axis-native (mesh-less, d≥3 — all-Cartesian by construction):
        # the cell volume is the tensor-product cell measure, the
        # iterated outer product of the per-axis widths. 2-D per-face
        # areas have a different shape and feed no matvec caller — None.
        if mesh is not None:
            self._volumes: np.ndarray = mesh.volumes
            self._areas: np.ndarray | None = (
                mesh.areas if self.ndim == 1 else None
            )
        else:
            self._volumes = reduce(np.multiply.outer, self.axis_widths)
            self._areas = None
        # ``nx`` = spatial_shape[0] sugar.
        self.nx: int = self.spatial_shape[0]
        # The whole-mesh coordinate system derives from the axes
        # (multi-axis tuples are all-Cartesian by construction).
        self.coord: "CoordSystem" = _axis_coord_system(self.axes)

        # ── Materials consistency validation (Issue #197 PR-TYPED-0) ──
        # Two checks at construction time:
        #   (1) every material id used in ``mat_map`` must have an entry
        #       in ``materials`` — otherwise downstream code would key
        #       into an undefined material.
        #   (2) all materials must agree on ``ng`` — one uniform group
        #       structure; heterogeneous ``ng`` is a homogenization-step
        #       concern that must precede method-mesh construction.
        # Both surface at construction, NOT lazily — the failure mode
        # (operators built on a bad mesh) is action-at-a-distance
        # otherwise.
        self._validate_materials()
        # Trigger ``ng`` property's consistency check eagerly so
        # mismatched-ng materials raise at construction time.
        _ = self.ng

    # ── Meshless construction (infinite homogeneous medium) ───────────

    @classmethod
    def from_materials(cls, materials: "dict[int, Mixture]") -> "MaterialMesh":
        r"""Meshless single-cell carrier for an infinite homogeneous medium.

        Builds the degenerate :class:`MaterialMesh` an *infinite-medium*
        problem lives on: one Cartesian cell of unit width holding a
        single material (id ``0``), with **no** legacy mesh adapter
        (:attr:`mesh` is ``None``).  This is the phase space the
        homogeneous :math:`k_\infty` solver assembles the transport
        operators on — :math:`(C - K_\text{iso} - F/k)\,\varphi = 0` with
        streaming :math:`L` dropped (zero in an infinite medium) — so the
        whole infinite-medium spectrum is computed through the *same*
        operator algebra the meshed S\ :sub:`N` solver uses, not a bespoke
        matrix.

        The lone cell carries material id ``0``, so ``materials`` MUST
        contain key ``0`` (the canonical id for a one-region medium); any
        additional entries are retained in the dict but unused by the
        single cell.  The cell has unit volume, so :attr:`volume_measure`
        weights it ``1.0`` — :math:`k_\infty` is a production/absorption
        *ratio*, invariant to the cell volume, so the unit choice is
        immaterial to the eigenvalue and keeps reaction rates equal to the
        bare :math:`\langle\Sigma,\varphi\rangle` group contractions.

        Parameters
        ----------
        materials : dict mapping material id to Mixture
            Macroscopic cross sections; MUST contain key ``0``.  The
            single source of truth for both the cross sections and the
            group count :attr:`ng` (derived, never passed — no twin).

        Returns
        -------
        MaterialMesh
            A 1-cell, single-region, mesh-less carrier
            (``spatial_shape == (1,)``, ``mat_map == [0]``,
            ``mesh is None``).
        """
        obj = cls.__new__(cls)
        # One Cartesian cell of unit width.  BCs default to None
        # (mesh-level reflective — the physical infinite-medium closure)
        # but are never read: the homogeneous solver drops the streaming
        # operator, so no boundary trace is ever applied.  ``mat_map=None``
        # → a single material with id 0 (``np.zeros((1,), int)``).
        obj._init_data(
            axes=(AxisMesh(edges=np.array([0.0, 1.0])),),
            mesh=None,
            mat_map=None,
            materials=materials,
        )
        return obj

    # NOTE: a GENERAL axis-native ``MaterialMesh.from_axes`` (arbitrary
    # cell count / coordinate system) is still intentionally NOT provided.
    # ``SNMesh.from_axes`` already exists with a different
    # (quadrature-bearing) signature, so a base ``from_axes`` here would be
    # an incompatible override; and the only axis-native consumer so far is
    # the meshless 1-cell :meth:`from_materials` above (which manufactures
    # its own trivial axis). ``.homogenize`` still builds via the legacy
    # ``MaterialMesh(coarse_mesh, materials)`` ctor. Defer the general form
    # until a real N-cell consumer exists (defer-until-≥2-instances).

    # ── Materials validation ──────────────────────────────────────────

    def _validate_materials(self) -> None:
        """Validate the materials dict against the mat_map.

        Every material id referenced in ``self.mat_map`` MUST appear as a
        key in ``self.materials``.  Failure surfaces here at construction
        time, not lazily inside a solver step.

        Raises
        ------
        ValueError
            If any ``mat_map`` id is missing from ``materials``; the
            error message shows both sets so the user can see the gap.
        """
        if not self.materials:
            raise ValueError(
                "MaterialMesh requires a non-empty materials dict; got "
                f"materials={self.materials!r}."
            )
        used_ids = set(int(x) for x in np.unique(self.mat_map))
        available_ids = set(self.materials.keys())
        missing = used_ids - available_ids
        if missing:
            raise ValueError(
                f"MaterialMesh.mat_map references material ids "
                f"{sorted(missing)} that are NOT in materials "
                f"(available ids: {sorted(available_ids)}; used ids: "
                f"{sorted(used_ids)})."
            )

    # ── Properties ────────────────────────────────────────────────────

    @property
    def ng(self) -> int:
        """Number of energy groups; uniform across all materials.

        Derived from ``self.materials``; the single source of truth for
        the group count.  All materials must agree on ``ng`` — a
        method-mesh requires one uniform group structure across the mesh.

        Raises
        ------
        InconsistentMaterialsError
            If materials disagree on ``ng``.  A homogenization /
            condensation step must precede method-mesh construction in
            that case.
        ValueError
            If ``self.materials`` is empty (caught at construction time
            by ``_validate_materials``).
        """
        if not self.materials:
            raise ValueError(
                "MaterialMesh.ng undefined — no materials.  Construct "
                "MaterialMesh(mesh, materials) with a non-empty materials "
                "dict."
            )
        ngs = {m.ng for m in self.materials.values()}
        if len(ngs) != 1:
            raise InconsistentMaterialsError(
                f"MaterialMesh requires uniform ng across all materials; "
                f"got ng values {sorted(ngs)} in materials dict with keys "
                f"{sorted(self.materials.keys())}.  Homogenize / condense "
                f"to a common group structure before method-mesh "
                f"construction."
            )
        return ngs.pop()

    @property
    def volumes(self) -> np.ndarray:
        """Cell volumes, shape ``spatial_shape`` (rank ``ndim``)."""
        return self._volumes

    @property
    def volume_measure(self):
        r"""Cell-volume :class:`~orpheus.numerics.measure.DiscreteMeasure`.

        The natural integration measure :math:`\mu_V = \sum_i V_i\,
        \delta_{c_i}` for volume-integrated rates (keff
        production/absorption; the flux·volume weighting of
        homogenization).  Consumers read THIS property, not
        ``mesh.mesh.volume_measure`` — the legacy mesh adapter is a
        construction provenance detail, not a data path.

        Delegates to the legacy dataclass's measure while the adapter is
        present (bit-identity: same atoms, same construction — including
        the ``precomputed_volumes`` escape hatch and the curvilinear
        volume formulas the dataclass owns).  Axis-native
        (``self.mesh is None``, d≥3): the rank-d analogue — atoms are the
        cell-centre tuples ordered with ``np.meshgrid(..., indexing='ij')``
        (the same layout ``volumes.ravel()`` exposes), weights the
        flattened cell volumes.
        """
        if self.mesh is not None:
            return self.mesh.volume_measure
        from orpheus.numerics.measure import DiscreteMeasure
        centers = [0.5 * (ax.edges[:-1] + ax.edges[1:]) for ax in self.axes]
        grids = np.meshgrid(*centers, indexing="ij")
        nodes = np.stack([g.ravel() for g in grids], axis=-1)
        return DiscreteMeasure(
            nodes=nodes,
            weights=self.volumes.ravel(),
            support=f"spatial_R{self.ndim}",
        )

    @property
    def areas(self) -> np.ndarray:
        """Face areas at each radial edge, shape (nx+1,) (1-D meshes).

        Sourced from :attr:`Mesh1D.areas`.  Cartesian slab returns an
        array of ones; cylinder returns :math:`2\\pi r`; sphere returns
        :math:`4\\pi r^2`.

        Raises
        ------
        AttributeError
            If accessed on a 2-D mesh (face areas have a different shape
            and are not consumed by today's matvec callers).
        """
        if self._areas is None:
            raise AttributeError(
                "MaterialMesh.areas is not defined for 2-D meshes; "
                "face-area data lives in the underlying Mesh2D directly."
            )
        return self._areas

    @property
    def ndim(self) -> int:
        """Number of spatial dimensions; equals ``len(self.axes)``."""
        return len(self.axes)

    @property
    def spatial_shape(self) -> tuple[int, ...]:
        r"""Per-axis cell counts ``(n_0, n_1, ...)``.

        The canonical dim-agnostic shape descriptor.  Every dim-agnostic
        shape reader (typed-field factories, pack convention, sweep DAG)
        reads from here. ``self.nx`` is sugar for ``spatial_shape[0]``.
        """
        return _axis_spatial_shape(self.axes)

    # ── Macroscopic XS field ──────────────────────────────────────────

    def material_xs_field(self) -> "MaterialXSField":
        """Build the macroscopic XS field from this mesh's materials.

        Returns a
        :class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField`
        wrapping the per-material :class:`Mixture` data plus this mesh's
        ``mat_map`` — the single source of truth for both per-cell and
        per-material XS access used by every transport operator.

        Lazy import of :mod:`.material_xs_field` to avoid a circular
        dependency at module import time.
        """
        from orpheus.transport.mesh.material_xs_field import MaterialXSField
        return MaterialXSField.from_mesh(self)
