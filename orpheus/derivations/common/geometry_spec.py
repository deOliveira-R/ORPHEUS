r"""Method-agnostic geometry specification for cross-method case construction.

:class:`GeometrySpec` is a specification — a small, immutable description
of ``(geometry_kind, critical_dimension_{mfp,cm}, n_groups, mat_id,
bc_left, bc_right, regions)`` — that BOTH production solvers AND continuous
reference solvers can consume at the case-description boundary:

* **Production solvers** (``solve_cp``, ``solve_sn``, etc.) consume
  the discretised :class:`~orpheus.geometry.mesh.Mesh1D` produced by
  :meth:`GeometrySpec.build`.
* **Continuous reference solvers** (``solve_fn_*``,
  ``solve_greens_function_*``, etc.) consume the descriptor scalars
  directly via :attr:`GeometrySpec.domain_extent_cm`,
  :attr:`GeometrySpec.critical_dimension_cm`,
  :attr:`GeometrySpec.bc_left`, :attr:`GeometrySpec.bc_right` —
  with no need to discretise into cells.

The class lived inside
``orpheus.derivations.continuous.sood_registry.la13511`` until the
2026-05-03 R0.5 promotion (geometry-handling unification audit at
``.claude/scratch/geometry_handling_unification_audit.md`` §"Q3.b
Promote-and-unify"). The audit identified that:

1. The descriptor is already method-agnostic — its attributes
   describe the SHAPE of the domain, not any particular solver.
2. Both production-side (sood_registry → ``Mesh1D`` → ``solve_cp``)
   and reference-side (cross-method adapter → ``solve_fn_slab``)
   call sites benefit from a single canonical specification.
3. The promotion is independent of the trajectory_resolvent
   R1/R2/R3 chord-oracle / power-iterate refactors.

Originally named ``MeshTemplate`` (the name described what the class
*produced* — a :class:`Mesh1D` for discrete consumers). Renamed to
``GeometrySpec`` on 2026-05-03 to describe what the class *is* (a
geometry specification — geometry kind, critical dimension, materials,
BCs, group count) — the name now matches the abstraction's
method-agnostic role and is symmetric across discrete and
reference-method consumers.

Multi-region support
--------------------

Step 3 of the input-cleanup track (2026-05-04) added
:class:`Region` and the :attr:`GeometrySpec.regions` field. This
promotes multi-region geometry — reflected slabs, layered
sphere/cylinder reactors — to a first-class descriptor, retiring
the legacy ``case.notes`` ``key=value`` parsing path used by the
cross-method reflected-slab adapter.

Convention:

* ``regions`` is an ORDERED tuple of :class:`Region`. The order is
  **inside-out** for sphere / cylinder (innermost layer first; outer
  layer last) and **left-to-right** for slabs (region at ``x=0``
  first; region at ``x=L`` last).
* Each region's ``outer_thickness_{mfp,cm}`` is a literal layer
  thickness (NOT cumulative). The sum over regions equals the
  total domain extent (i.e. the slab full width or the
  sphere/cylinder outer radius).
* ``regions=None`` (the default) reproduces the single-region
  behaviour — ``build()`` yields a homogeneous mesh of material
  ``mat_id`` filling the full :attr:`domain_extent_cm`.

References
----------

* :doc:`/skills/algebra-of-record` §"Branch 1 / Branch 2 separation"
  — :class:`GeometrySpec` lives at the math layer, NOT inside any
  specific solver's algebra-of-record.
* ``.claude/scratch/geometry_handling_unification_audit.md`` —
  the audit memo that identified the promotion opportunity.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.factories import Zone, mesh1d_from_zones
from orpheus.geometry.mesh import BC, Mesh1D


_GEOMETRY_TO_COORD: dict[str, CoordSystem] = {
    "slab": CoordSystem.CARTESIAN,
    "sphere": CoordSystem.SPHERICAL,
    "cylinder": CoordSystem.CYLINDRICAL,
}


@dataclass(frozen=True)
class Region:
    r"""One layer of a multi-region :class:`GeometrySpec`.

    A :class:`Region` carries a material ID and an outer-layer
    thickness in two unit conventions (mean free paths and cm) that
    must agree under the case's primary :math:`\Sigma_t`. The dual
    representation mirrors :class:`GeometrySpec`'s
    ``critical_dimension`` pair — both forms travel together so
    consumers can pick the unit they need without re-deriving through
    ``truth_value / sigma_t``.

    Parameters
    ----------
    mat_id : int
        Material identifier matching a key in
        ``case.materials`` / ``case.registry_case.materials``. By
        convention, :data:`tests.cross_method.cases` uses ``0`` for
        the multiplying core and ``1`` for reflectors / moderators
        (other case sets have their own conventions; the geometry
        layer is agnostic).
    outer_thickness_mfp : float | None
        Layer thickness in mean free paths
        (:math:`\Delta_{\rm mfp} = \Sigma_t \cdot \Delta_{\rm cm}`).
        Either ``outer_thickness_mfp`` OR ``outer_thickness_cm``
        must be set; both may be set when the dual form has been
        pre-computed.
    outer_thickness_cm : float | None
        Layer thickness in cm. Same convention as
        :attr:`outer_thickness_mfp`.

    Notes
    -----
    Thicknesses are LITERAL (not cumulative). The
    :class:`GeometrySpec.regions` ordering — inside-out for
    sphere/cylinder, left-to-right for slab — is enforced at the
    :class:`GeometrySpec` level, not here. A :class:`Region` is
    just a layer description.

    Single-region geometries (sphere with one core, slab with one
    homogeneous medium) do NOT need :class:`Region` — they're
    described directly by ``GeometrySpec.mat_id`` +
    ``GeometrySpec.critical_dimension_cm``. :class:`Region` is for
    multi-region cases (NM 1980 reflected slab, layered fuel /
    reflector spheres, etc.).
    """

    mat_id: int
    outer_thickness_mfp: float | None = None
    outer_thickness_cm: float | None = None

    def __post_init__(self) -> None:
        if (
            self.outer_thickness_mfp is None
            and self.outer_thickness_cm is None
        ):
            raise ValueError(
                "Region: at least one of outer_thickness_mfp or "
                "outer_thickness_cm must be set."
            )
        if (
            self.outer_thickness_cm is not None
            and self.outer_thickness_cm <= 0.0
        ):
            raise ValueError(
                f"Region: outer_thickness_cm must be > 0; "
                f"got {self.outer_thickness_cm!r}"
            )
        if (
            self.outer_thickness_mfp is not None
            and self.outer_thickness_mfp <= 0.0
        ):
            raise ValueError(
                f"Region: outer_thickness_mfp must be > 0; "
                f"got {self.outer_thickness_mfp!r}"
            )


@dataclass(frozen=True)
class GeometrySpec:
    """Geometry + critical-dimension specification for a benchmark case.

    Method-agnostic: stores the shape of the domain (slab /
    sphere / cylinder / infinite / ISLC), the critical dimension as
    published, the group count, and the boundary conditions. The
    :meth:`build` method constructs a concrete :class:`Mesh1D` at the
    requested refinement using
    :func:`orpheus.geometry.factories.mesh1d_from_zones`.

    For ``geometry == "infinite"`` no mesh exists: :meth:`build`
    raises :class:`ValueError`. Infinite-medium consumers should use
    :func:`...sood_registry.builders.build_materials` and ignore the
    mesh entirely.

    Parameters
    ----------
    geometry : str
        One of ``"infinite"``, ``"slab"``, ``"sphere"``, ``"cylinder"``,
        ``"ISLC"``.
    critical_dimension_mfp : float | None
        Critical radius (sphere/cylinder) or half-thickness (slab) in
        mean free paths. ``None`` for infinite-medium cases. For
        multi-region cases this typically holds the published core
        critical scalar (e.g. NM 1980 :math:`\tau`); see notes.
    critical_dimension_cm : float | None
        Same as above but in cm.
    n_groups : int
        Number of energy groups (1, 2, 3, or 6).
    mat_id : int
        Material identifier used by the **single-region** path
        (when ``regions is None``). For multi-region cases the
        per-region :attr:`Region.mat_id` takes precedence and this
        field is ignored.
    bc_left : BC
        Inner / left boundary condition. For sphere/cylinder default
        is :attr:`BC.reflective` (centreline); for slab default is
        :attr:`BC.vacuum`.
    bc_right : BC
        Outer / right boundary condition. Default
        :attr:`BC.vacuum`.
    regions : tuple[Region, ...] | None
        Ordered layer stack for multi-region geometries. ``None``
        (default) means "single homogeneous region of ``mat_id``
        filling the full :attr:`domain_extent_cm`." Ordering is
        **inside-out** (sphere / cylinder) or **left-to-right**
        (slab). For symmetric reflected slabs (NM 1980 / Sood Table
        10) the natural ordering is
        ``(reflector, core, reflector)``. The sum of
        :attr:`Region.outer_thickness_cm` over all regions defines
        :attr:`domain_extent_cm`.

    Notes
    -----
    **Single-region default.** When ``regions is None`` the spec
    describes a single homogeneous medium of material ``mat_id``
    filling the full domain. This is the original behaviour — every
    pre-Step-3 case continues to work without modification.

    **Multi-region domain extent.** For multi-region slabs the
    natural domain is the FULL slab (left edge at ``x = 0``, right
    edge at ``x = Σ thicknesses``). This contrasts with the
    single-region slab convention where ``domain_extent_cm`` is
    derived from ``2 * critical_dimension_cm`` (core
    half-thickness). When ``regions`` is set,
    ``domain_extent_cm`` is computed from the region thicknesses,
    NOT from ``critical_dimension_cm``. The
    :attr:`critical_dimension_cm` field then stores a
    method-specific scalar (e.g. the published core half-thickness
    for reflected slab) that consumers may use as a diagnostic but
    is NOT the mesh extent.

    **Mesh density.** :meth:`build`'s ``n_cells`` argument is the
    TOTAL cell budget; for multi-region cases it is distributed
    across regions in proportion to their cm thickness with a floor
    of one cell per region. Callers needing finer per-region
    control should split the build manually using
    :func:`~orpheus.geometry.factories.mesh1d_from_zones`.
    """

    geometry: str
    critical_dimension_mfp: float | None
    critical_dimension_cm: float | None
    n_groups: int
    mat_id: int = 0
    bc_left: BC = field(default_factory=lambda: BC.vacuum)
    bc_right: BC = field(default_factory=lambda: BC.vacuum)
    regions: tuple[Region, ...] | None = None

    def __post_init__(self) -> None:
        if self.geometry not in {"infinite", "slab", "sphere", "cylinder", "ISLC"}:
            raise ValueError(f"Unknown geometry {self.geometry!r}")
        if self.geometry == "infinite":
            if self.critical_dimension_cm is not None:
                raise ValueError(
                    "infinite geometry must have critical_dimension_cm=None"
                )
            if self.regions is not None:
                raise ValueError(
                    "infinite geometry cannot carry regions (no spatial mesh)"
                )
        else:
            if self.critical_dimension_cm is None:
                raise ValueError(
                    f"geometry={self.geometry!r} requires critical_dimension_cm"
                )
        if self.regions is not None:
            self._validate_regions()

    def _validate_regions(self) -> None:
        """Check region consistency for multi-region specs."""
        regions = self.regions
        assert regions is not None  # caller-guarded
        if len(regions) == 0:
            raise ValueError(
                "GeometrySpec.regions must be non-empty when set; "
                "use regions=None for single-region."
            )
        if self.geometry in {"infinite", "ISLC"}:
            # infinite already raised above; ISLC is currently
            # unsupported for any mesh path.
            raise ValueError(
                f"GeometrySpec.regions is not supported for "
                f"geometry={self.geometry!r}."
            )
        # Every region must carry an outer_thickness_cm so we can
        # compute the geometric extent. The mfp form is optional but
        # recommended for cross-method consumers.
        cm_thicknesses: list[float] = []
        for k, region in enumerate(regions):
            if region.outer_thickness_cm is None:
                raise ValueError(
                    f"GeometrySpec.regions[{k}] is missing "
                    f"outer_thickness_cm — required to construct the "
                    f"mesh. (outer_thickness_mfp alone is not enough; "
                    f"the mesh lives in cm.)"
                )
            cm_thicknesses.append(float(region.outer_thickness_cm))
        total_cm = sum(cm_thicknesses)
        if total_cm <= 0.0:
            raise ValueError(
                f"GeometrySpec.regions sum of outer_thickness_cm must "
                f"be > 0; got {total_cm}"
            )

    @property
    def domain_extent_cm(self) -> float:
        """Total mesh extent in cm — what :meth:`build` actually constructs.

        Conventions:

        * **Slab, single-region** (``regions is None``):
          ``2 * critical_dimension_cm`` (full symmetric slab
          ``[0, 2a]`` with vacuum BCs at both ends — F_N convention).
        * **Sphere / cylinder, single-region** (``regions is None``):
          ``critical_dimension_cm`` (radius; mesh is ``[0, R]`` with
          reflective BC at the centre and vacuum at the outer surface
          unless overridden).
        * **Multi-region** (``regions is not None``): the **sum** of
          :attr:`Region.outer_thickness_cm` over all regions. For
          slabs this is the FULL slab width (left-to-right); for
          sphere/cylinder this is the outer radius.
        * **Infinite / ISLC**: undefined, raises.
        """
        if self.geometry == "infinite":
            raise ValueError("infinite geometry has no domain_extent")
        if self.geometry == "ISLC":
            raise NotImplementedError("ISLC domain_extent not implemented")
        if self.regions is not None:
            return float(
                sum(
                    r.outer_thickness_cm
                    for r in self.regions
                    if r.outer_thickness_cm is not None
                )
            )
        assert self.critical_dimension_cm is not None  # narrowed above
        if self.geometry == "slab":
            return 2.0 * float(self.critical_dimension_cm)
        return float(self.critical_dimension_cm)

    def build(self, n_cells: int = 64) -> Mesh1D:
        r"""Construct a :class:`Mesh1D` at the published critical dimension.

        Conventions (see :attr:`domain_extent_cm`):

        * **Slab** (``geometry == "slab"``, single-region): builds the
          FULL symmetric slab ``[0, 2a]`` where
          ``a = critical_dimension_cm``. This is the F_N method's
          natural domain — :math:`a` is the half-thickness, the
          published critical configuration is the full slab. Default
          BCs are vacuum at both ends.
        * **Sphere / cylinder** (single-region): builds ``[0, R]``
          where ``R = critical_dimension_cm``. Default BCs:
          reflective at ``r = 0`` (centreline / axis), vacuum at the
          outer surface.
        * **Multi-region**: builds a multi-zone mesh by passing the
          ordered :attr:`regions` to
          :func:`~orpheus.geometry.factories.mesh1d_from_zones`. The
          ``n_cells`` budget is split across regions proportionally
          to their cm thickness with a floor of 1 cell per region.

        Parameters
        ----------
        n_cells : int
            Total number of cells. For single-region the mesh has
            exactly ``n_cells`` equal-volume cells. For multi-region
            the budget is distributed across regions (see notes).

        Returns
        -------
        Mesh1D
            A 1-D mesh with the requested cells, materials populated
            from :attr:`mat_id` (single-region) or :attr:`regions`
            (multi-region), and BCs from this specification.

        Raises
        ------
        ValueError
            If ``geometry == "infinite"``.
        NotImplementedError
            If ``geometry == "ISLC"``.
        """
        if self.geometry == "infinite":
            raise ValueError(
                "infinite-medium cases have no mesh; consume materials directly"
            )
        if self.geometry == "ISLC":
            raise NotImplementedError(
                "ISLC mesh construction is not implemented in the first slice"
            )
        coord = _GEOMETRY_TO_COORD[self.geometry]

        if self.regions is None:
            zones = [Zone(
                outer_edge=self.domain_extent_cm,
                mat_id=self.mat_id,
                n_cells=n_cells,
            )]
        else:
            zones = self._build_multi_region_zones(n_cells)

        mesh = mesh1d_from_zones(zones, coord=coord)
        # Re-stamp with BC fields (mesh1d_from_zones doesn't set BCs).
        return Mesh1D(
            edges=mesh.edges,
            mat_ids=mesh.mat_ids,
            coord=mesh.coord,
            precomputed_volumes=mesh.precomputed_volumes,
            bc_left=self.bc_left,
            bc_right=self.bc_right,
        )

    def _build_multi_region_zones(self, n_cells: int) -> list[Zone]:
        """Distribute ``n_cells`` across the region stack.

        Each region gets ``round(n_cells * Δ_k / Σ Δ)`` cells with a
        floor of 1; the largest (cm) region absorbs any rounding
        residual so the total adds up to ``max(n_cells, len(regions))``
        — preserving the per-region floor without silently inflating
        n_cells when a tiny region is rounded up.
        """
        regions = self.regions
        assert regions is not None  # caller-guarded
        total = self.domain_extent_cm
        # Allocate cells proportional to cm thickness, floor 1.
        raw_alloc: list[int] = []
        for region in regions:
            assert region.outer_thickness_cm is not None  # validated
            ratio = float(region.outer_thickness_cm) / total
            raw_alloc.append(max(1, round(n_cells * ratio)))
        # Reconcile rounding so allocations sum to n_cells (or at
        # least len(regions) if n_cells < len(regions)).
        target = max(n_cells, len(regions))
        delta = target - sum(raw_alloc)
        if delta != 0:
            # Push residual onto the largest (cm) region.
            largest = max(
                range(len(regions)),
                key=lambda k: regions[k].outer_thickness_cm,
            )
            raw_alloc[largest] = max(1, raw_alloc[largest] + delta)
        # Build cumulative-edge zones.
        zones: list[Zone] = []
        cumulative = 0.0
        for region, n in zip(regions, raw_alloc):
            assert region.outer_thickness_cm is not None
            cumulative += float(region.outer_thickness_cm)
            zones.append(Zone(
                outer_edge=cumulative,
                mat_id=region.mat_id,
                n_cells=n,
            ))
        return zones


__all__ = ["GeometrySpec", "Region"]
