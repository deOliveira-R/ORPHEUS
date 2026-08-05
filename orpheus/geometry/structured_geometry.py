r"""Structured 1-D geometry — pure geometry layer, no mesh concerns.

A :class:`StructuredGeometry` describes the SHAPE of a 1-D problem domain:

* **What kind of geometry** (slab / cylinder / sphere → coordinate system).
* **How many regions** (single homogeneous medium or layered stack).
* **What boundary conditions** apply at the geometry's endpoints.

It is "structured" in the memory-layout sense: regions can be traversed
by a simple for-loop with no database lookup. This is the
counterpart to a future ``ConstructiveSolidGeometry`` (CSG) that would
describe unstructured / Boolean-composed shapes — but those need
different machinery (BVH lookups, ray–surface intersection) that is
out of scope for the first slice.

The class lives at the **geometry layer**, not the mesh layer, and
deliberately knows nothing about cell counts, discretization methods,
or numerics. The same :class:`StructuredGeometry` can be discretized
multiple ways for different studies (mesh refinement,
uniform-vs-equal-volume comparison, future temperature-aware schemes)
by passing different :class:`RegionMesh` tuples to
:meth:`Mesh1D.from_geometry`.

Architectural role
------------------

In the project's two-role solver split:

* **Reference solution generators** (``Billiard``, ``MomentSpace``,
  ``Spectrum``, ``BasisSpace``) consume
  :class:`StructuredGeometry` directly via their ``__init__`` —
  no mesh, no cell counts.
* **Discrete production solvers** (``solve_cp``, ``solve_sn``,
  ``solve_moc``, ``solve_mc``) consume a :class:`Mesh1D` built via
  :meth:`Mesh1D.from_geometry`. The discretization description (per-
  region cell counts and methods) is supplied at that build step,
  not on the geometry.

The geometry → mesh transition is the single explicit point where
discretization information enters the pipeline. Above it: pure
geometry + materials. Below it: solver-specific augmented meshes
(``CPMesh``, ``SNMesh``) built by each solver.

Geometry kinds and orbit-space classification
---------------------------------------------

The supported tags are uppercase three-letter mnemonics:

==========  ===============  ==================  ==============
Tag         Coordinate       Endpoints (BC)      Centreline
==========  ===============  ==================  ==============
``"SLB"``   Cartesian        2 (left, right)     n/a
``"CYL"``   Cylindrical      1 (outer)           implicit reflective
``"SPH"``   Spherical        1 (outer)           implicit reflective
==========  ===============  ==================  ==============

The endpoint count comes from the orbit-space classification of the
geometry's billiard table: SLB has two flat surfaces (orbit-space
rank 2); CYL and SPH have one outer surface plus an implicit
centreline reflection (orbit-space rank 1). This is the same
classification the trajectory_resolvent ``Billiard`` class uses.

When future geometries land that genuinely have two surfaces
(``HSPH`` for hollow sphere, ``ANN`` for annulus), they extend
this map with two-endpoint BC tuples.

Examples
--------

A bare-critical Sood sphere — single region, one outer BC::

    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.872),),
        bcs=(BC.vacuum,),
    )

A reflected-slab NM-1980 case — three regions (reflector | core |
reflector), two outer BCs::

    geom = StructuredGeometry(
        geometry="SLB",
        regions=(
            Region(mat_id=1, outer_thickness_cm=0.5),  # reflector left
            Region(mat_id=0, outer_thickness_cm=2.0),  # core
            Region(mat_id=1, outer_thickness_cm=0.5),  # reflector right
        ),
        bcs=(BC.vacuum, BC.vacuum),
    )

A PWR pin cell via the Wigner–Seitz factory::

    geom = StructuredGeometry.wigner_seitz_pin_cell(
        r_fuel=0.9, r_clad=1.1, pitch=3.6,
    )

Building a mesh from any of these — discretization specified at the
mesh layer::

    from orpheus.geometry.mesh import Mesh1D, RegionMesh
    mesh = Mesh1D.from_geometry(geom, region_meshes=(
        RegionMesh(n_cells=10),
        RegionMesh(n_cells=3),
        RegionMesh(n_cells=7),
    ))
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .coord import CoordSystem
from .mesh import BC

if TYPE_CHECKING:
    from .boundary import BoundaryTraceLaw


# ═══════════════════════════════════════════════════════════════════════
# Geometry-tag → coordinate-system map
# ═══════════════════════════════════════════════════════════════════════

_GEOMETRY_TO_COORD: dict[str, CoordSystem] = {
    "SLB": CoordSystem.CARTESIAN,
    "CYL": CoordSystem.CYLINDRICAL,
    "SPH": CoordSystem.SPHERICAL,
}

# Endpoint counts by geometry tag (orbit-space rank for the BC tuple).
# SLB: 2 endpoints (left, right). CYL/SPH: 1 endpoint (outer);
# the centreline is implicit reflective.
_GEOMETRY_TO_N_ENDPOINTS: dict[str, int] = {
    "SLB": 2,
    "CYL": 1,
    "SPH": 1,
}


# ═══════════════════════════════════════════════════════════════════════
# Region — one layer of a StructuredGeometry
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class Region:
    """A layer of a :class:`StructuredGeometry`. Pure geometry — no mesh
    concerns.

    A region carries a material identifier and the literal cm thickness
    of its outer-bound layer. It does NOT carry a cell count or a
    discretization method — those live on :class:`RegionMesh` at the
    mesh layer and are supplied when building a :class:`Mesh1D` from
    the geometry.

    Parameters
    ----------
    mat_id : int
        Material identifier matching a key in the ``materials:
        dict[int, Mixture]`` payload that downstream consumers receive
        alongside the geometry.
    outer_thickness_cm : float
        Literal thickness of this region in cm. NOT cumulative
        (see :class:`StructuredGeometry` for ordering convention).
        Must be positive.

    Examples
    --------

    Single-region homogeneous sphere::

        Region(mat_id=0, outer_thickness_cm=2.872)

    Three-region reflected slab (the tuple ordering is left-to-right
    on the slab; for sphere/cylinder it is inside-out)::

        (
            Region(mat_id=1, outer_thickness_cm=0.5),  # reflector left
            Region(mat_id=0, outer_thickness_cm=2.0),  # core
            Region(mat_id=1, outer_thickness_cm=0.5),  # reflector right
        )

    See Also
    --------
    StructuredGeometry : The class that holds a tuple of Regions.
    RegionMesh : The mesh-layer per-region descriptor (n_cells +
        method) that pairs with a Region at mesh-build time.
    """

    mat_id: int
    outer_thickness_cm: float

    def __post_init__(self) -> None:
        if not isinstance(self.mat_id, (int, np.integer)):
            raise TypeError(
                f"Region.mat_id must be int, got {type(self.mat_id).__name__}"
            )
        if not isinstance(self.outer_thickness_cm, (int, float)):
            raise TypeError(
                f"Region.outer_thickness_cm must be a real number, "
                f"got {type(self.outer_thickness_cm).__name__}"
            )
        if self.outer_thickness_cm <= 0.0:
            raise ValueError(
                f"Region.outer_thickness_cm must be > 0; "
                f"got {self.outer_thickness_cm!r}"
            )


# ═══════════════════════════════════════════════════════════════════════
# StructuredGeometry — 1-D structured geometry (the geometry layer)
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class StructuredGeometry:
    """1-D structured geometry — pure geometry, no mesh, no materials data.

    The geometry layer sits between the registry (problem definition)
    and the mesh layer (discretization). It carries:

    * Geometry kind (``"SLB"`` / ``"CYL"`` / ``"SPH"``) — determines
      the coordinate system and the orbit-space endpoint count.
    * An ordered tuple of :class:`Region` instances — the layered
      material structure. Single-region cases are 1-tuples.
    * A tuple of :class:`BC` instances — boundary conditions at the
      geometry's endpoints. The tuple length matches the geometry's
      endpoint count (SLB → 2, CYL/SPH → 1).

    What is NOT here
    ----------------

    * **No cell counts.** Discretization is a mesh concern. The same
      :class:`StructuredGeometry` can be meshed differently for
      different studies. Supply :class:`RegionMesh` instances to
      :meth:`Mesh1D.from_geometry` at mesh-build time.
    * **No critical-dimension or mfp values.** Those are
      registry-truth artifacts (a published critical configuration
      "happens to be" this size). They live on registry truth records
      (e.g., ``La13511Truth.critical_dimension_mfp``), not on the
      geometry.
    * **No energy-group count.** That comes from the materials:
      ``mixture.SigT.shape[0]``.
    * **No infinite-medium kind.** Solving an infinite medium means
      using :func:`solve_homogeneous_infinite` (no geometry) or
      :meth:`MomentSpace.solve_kinf` (no geometry). A "functionally
      infinite" finite domain can be expressed via reflective BCs.

    Parameters
    ----------
    geometry : str
        One of ``"SLB"`` (Cartesian slab), ``"CYL"`` (cylindrical),
        ``"SPH"`` (spherical).
    regions : tuple[Region, ...]
        Ordered layer stack. Ordering is **inside-out** for sphere /
        cylinder (innermost region first), **left-to-right** for
        slab. Always a tuple — single-region geometries use a
        1-tuple ``(Region(...),)``.
    bcs : tuple[BC, ...]
        Boundary conditions at the geometry's endpoints. Length
        matches :attr:`n_endpoints`. For ``"SLB"``: 2 BCs (left,
        right). For ``"CYL"`` / ``"SPH"``: 1 BC (outer surface; the
        centreline is implicit reflective).

    Notes
    -----

    **Slab convention.** :attr:`domain_extent_cm` is the FULL slab
    width (sum of all region thicknesses, end-to-end). This is the
    production-natural convention — a user with "a slab of thickness
    L" passes ``L`` directly. F_N's natural half-thickness ``a = L /
    2`` is computed by reference solvers internally.

    **Materials coupling.** The geometry knows ``mat_id`` integers
    but not the materials themselves. Solvers receive the
    ``materials: dict[int, Mixture]`` payload alongside the geometry
    and resolve the dict at solve time. This decoupling means the
    same materials can be reused with different geometries and vice
    versa.

    **Boundary-condition resolution.** :class:`BC` is a tag (``kind`` +
    optional ``params``). Each solver's augmented mesh resolves what
    a tag means via its own ``BC_REGISTRY``. The geometry layer never
    interprets BCs.

    Examples
    --------

    See module docstring.

    See Also
    --------
    Region : One layer of a multi-region geometry.
    Mesh1D.from_geometry : Build a discrete mesh from this geometry.
    """

    geometry: str
    regions: tuple[Region, ...]
    bcs: "tuple[BC | BoundaryTraceLaw, ...]"

    def __post_init__(self) -> None:
        if self.geometry not in _GEOMETRY_TO_COORD:
            raise ValueError(
                f"StructuredGeometry.geometry must be one of "
                f"{sorted(_GEOMETRY_TO_COORD)}; got {self.geometry!r}"
            )
        if not isinstance(self.regions, tuple):
            raise TypeError(
                f"StructuredGeometry.regions must be a tuple, "
                f"got {type(self.regions).__name__}"
            )
        if len(self.regions) == 0:
            raise ValueError(
                "StructuredGeometry.regions must be non-empty; "
                "single-region geometries use a 1-tuple."
            )
        for k, region in enumerate(self.regions):
            if not isinstance(region, Region):
                raise TypeError(
                    f"StructuredGeometry.regions[{k}] must be a Region, "
                    f"got {type(region).__name__}"
                )

        if not isinstance(self.bcs, tuple):
            raise TypeError(
                f"StructuredGeometry.bcs must be a tuple, "
                f"got {type(self.bcs).__name__}"
            )
        expected_n_bcs = _GEOMETRY_TO_N_ENDPOINTS[self.geometry]
        if len(self.bcs) != expected_n_bcs:
            raise ValueError(
                f"StructuredGeometry: geometry={self.geometry!r} "
                f"requires {expected_n_bcs} BC(s); got {len(self.bcs)}."
            )
        # A declaration is EITHER a ``BC`` tag or an already-typed
        # ``BoundaryTraceLaw``. The law arm exists because a tag cannot
        # express a law that carries a FUNCTION: ``BC.params`` is
        # ``dict[str, float]``, so a prescribed inflow whose source is a
        # manufactured solution restricted to a face has no tag spelling.
        # Before this arm such a law could only be installed by mutating a
        # constructed mesh's resolved ``bc`` dict — which the public solver
        # entry points then DISCARDED, because they rebuild the method mesh
        # from the raw geometry. Declaring on the geometry is what makes the
        # law survive that rebuild.
        #
        # Imported lazily: ``orpheus.geometry.boundary`` transitively loads
        # THIS module, so a top-level import cycles.
        from orpheus.geometry.boundary import BoundaryTraceLaw

        for k, bc in enumerate(self.bcs):
            if not isinstance(bc, (BC, BoundaryTraceLaw)):
                raise TypeError(
                    f"StructuredGeometry.bcs[{k}] must be a BC tag or a "
                    f"BoundaryTraceLaw instance, got {type(bc).__name__}"
                )

    @property
    def coord(self) -> CoordSystem:
        """Coordinate system implied by :attr:`geometry`."""
        return _GEOMETRY_TO_COORD[self.geometry]

    @property
    def n_endpoints(self) -> int:
        """Number of geometry endpoints carrying a BC.

        SLB → 2 (left, right). CYL/SPH → 1 (outer; centreline is
        implicit reflective at the coordinate origin).
        """
        return _GEOMETRY_TO_N_ENDPOINTS[self.geometry]

    @property
    def domain_extent_cm(self) -> float:
        """Total geometric extent in cm — sum of region thicknesses.

        For ``"SLB"`` this is the full slab width (production
        convention). For ``"CYL"`` / ``"SPH"`` this is the outer
        radius. Always positive.
        """
        return float(sum(r.outer_thickness_cm for r in self.regions))

    # ─────────────────────────────────────────────────────────────────
    # Construction helpers — earn their keep when the conceptual
    # transformation is non-trivial
    # ─────────────────────────────────────────────────────────────────

    @classmethod
    def wigner_seitz_pin_cell(
        cls,
        *,
        r_fuel: float = 0.9,
        r_clad: float = 1.1,
        pitch: float = 3.6,
        bcs: tuple[BC, ...] = (BC("white"),),
    ) -> "StructuredGeometry":
        r"""Wigner–Seitz equivalent pin-cell geometry.

        Replaces a square unit cell of side ``pitch`` with a cylinder
        of equal cross-sectional area:

        .. math::

            r_{\rm cell} = \frac{\rm pitch}{\sqrt{\pi}}

        The resulting cylindrical geometry has three regions:
        fuel pellet (``mat_id=2``), cladding (``mat_id=1``), and
        coolant (``mat_id=0``). Region thicknesses are the *literal
        cm extents* of each annulus:

        * fuel:    ``r_fuel - 0`` (inner-most layer)
        * clad:    ``r_clad - r_fuel``
        * coolant: ``r_cell - r_clad``

        Default outer BC is :class:`BC("white") <BC>` — the unit-cell
        symmetry assumption that maps a periodic lattice to a single
        cell with isotropic re-entry probability.

        Parameters
        ----------
        r_fuel : float
            Outer radius of the fuel pellet (cm). Default 0.9.
        r_clad : float
            Outer radius of the cladding (cm). Default 1.1.
        pitch : float
            Square unit-cell side length (cm). Default 3.6.
        bcs : tuple[BC, ...]
            Outer-surface boundary condition. CYL geometry has 1
            endpoint, so this is a 1-tuple. Default ``(BC("white"),)``.

        Returns
        -------
        StructuredGeometry
            A CYL geometry with three regions in the conventional
            (fuel, clad, coolant) ordering.

        See Also
        --------
        Mesh1D.from_geometry : Discretize this geometry into a mesh.
        """
        r_cell = pitch / np.sqrt(np.pi)
        return cls(
            geometry="CYL",
            regions=(
                Region(mat_id=2, outer_thickness_cm=float(r_fuel)),
                Region(mat_id=1, outer_thickness_cm=float(r_clad - r_fuel)),
                Region(mat_id=0, outer_thickness_cm=float(r_cell - r_clad)),
            ),
            bcs=bcs,
        )

    @classmethod
    def pwr_slab_half_cell(
        cls,
        *,
        fuel_half: float = 0.9,
        clad_thick: float = 0.2,
        cool_thick: float = 0.7,
        bcs: tuple[BC, ...] = (BC("reflective"), BC("reflective")),
    ) -> "StructuredGeometry":
        r"""Cartesian 1-D PWR half-cell geometry: fuel | clad | coolant.

        Encodes the conceptual symmetry of a square PWR unit cell about
        the fuel centreline. The geometry starts at ``x = 0`` (the
        symmetry plane through the fuel centre) and extends outward
        through ``fuel_half`` (half the fuel thickness), ``clad_thick``,
        and ``cool_thick`` to the unit-cell boundary at the coolant
        side.

        Material IDs follow the project convention:
        ``2`` = fuel, ``1`` = clad, ``0`` = coolant.

        Default BCs are reflective on both ends — the standard infinite-
        lattice eigenvalue convention. Override ``bcs`` to (reflective,
        white) for an isolated unit cell with isotropic re-entry on
        the coolant boundary, or to (reflective, vacuum) for a
        right-vacuum boundary configuration.

        Parameters
        ----------
        fuel_half : float
            Half-thickness of the fuel slab (cm). Default 0.9.
        clad_thick : float
            Cladding thickness (cm). Default 0.2.
        cool_thick : float
            Coolant thickness (cm). Default 0.7.
        bcs : tuple[BC, ...]
            ``(bc_left, bc_right)`` BCs. Default
            ``(BC("reflective"), BC("white"))``.

        Returns
        -------
        StructuredGeometry
            An SLB geometry with three regions in the conventional
            (fuel, clad, coolant) ordering.

        See Also
        --------
        Mesh1D.from_geometry : Discretize this geometry into a mesh.
        wigner_seitz_pin_cell : The cylindrical equivalent.
        """
        return cls(
            geometry="SLB",
            regions=(
                Region(mat_id=2, outer_thickness_cm=float(fuel_half)),
                Region(mat_id=1, outer_thickness_cm=float(clad_thick)),
                Region(mat_id=0, outer_thickness_cm=float(cool_thick)),
            ),
            bcs=bcs,
        )


__all__ = [
    "Region",
    "StructuredGeometry",
]
