r"""Augmented geometry for S\ :sub:`N` discrete ordinates transport.

:class:`SNMesh` wraps a :class:`~geometry.mesh.Mesh1D` or
:class:`~geometry.mesh.Mesh2D` and precomputes the coordinate-specific
streaming stencil used by the transport sweep.

Three coordinate systems are supported: Cartesian (1D/2D), spherical
(1D), and cylindrical (1D).  Curvilinear geometries precompute angular
redistribution coefficients (:math:`\alpha`), the geometry factor
:math:`\Delta A/w`, and Morel--Montry angular closure weights.
"""

from __future__ import annotations

import warnings
from typing import ClassVar, Iterator, TYPE_CHECKING

import numpy as np

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.geometry.boundary import (
    BoundaryTraceLaw,
    ReflectiveBoundary,
    VacuumInflow,
)
from orpheus.geometry.boundary._bound_compat import _BoundBoundaryOperator
from orpheus.geometry.reduced_operator import (
    ReducedStreamingOperator,
    cylindrical_streaming,
    slab_streaming,
    spherical_streaming,
)
from .axis import (
    Axis1D,
    AxisCoord,
    AxisMesh,
    FaceLabel,
    RadialAxisMesh,
    axes_from_legacy_mesh,
    face_labels as _axis_face_labels,
    face_outflow_ordinates as _axis_face_outflow_ordinates,
    face_shape as _axis_face_shape,
    legacy_mesh_from_axes,
    n_unknowns_flat as _axis_n_unknowns_flat,
    spatial_shape as _axis_spatial_shape,
)
from .boundary_realizer import SNBoundaryRealizer
from .method_space import SNMethodSpace
from orpheus.numerics.quadrature import Quadrature
from .spatial.cell_update import CellUpdate, CellVisit
from .spatial.diamond import DiamondDifference
from .spatial.pole_angular_closure import (
    IdentityAngularClosure,
    MorelMontryAngularSweep,
    PoleAngularClosure,
    default_angular_closure_class,
)
from .sweep_graph import OctantLabel, SweepDependencyGraph

if TYPE_CHECKING:
    from orpheus.data.macro_xs.mixture import Mixture
    from .angular_flux import AngularFlux
    from .boundary_flux import BoundaryFlux
    from .harmonic_moment_field import HarmonicMomentField
    from .material_xs_field import MaterialXSField
    from .scalar_flux import ScalarFlux
    from .sources import IsotropicSource, PerOrdinateSource


class InconsistentMaterialsError(ValueError):
    """Raised when a materials dict has inconsistent metadata.

    Currently triggered when materials disagree on ``ng`` (number of
    energy groups).  :class:`SNMesh` requires a uniform group structure
    across all materials in its ``mat_map`` because every operator that
    consumes ``sn_mesh`` (L, C, S, F) assumes one ``ng``.  A
    homogenization step must precede SNMesh construction if the input
    materials carry different group structures.
    """


# ═══════════════════════════════════════════════════════════════════════
# SNMesh
# ═══════════════════════════════════════════════════════════════════════

class SNMesh:
    """Augmented geometry for the discrete ordinates method.

    Wraps a :class:`~geometry.mesh.Mesh1D` or :class:`~geometry.mesh.Mesh2D`
    and precomputes the streaming stencil (diamond-difference coefficients
    that depend only on geometry and angular quadrature, not on cross
    sections).

    For Cartesian geometry the stencil stores:

    * ``streaming_x[n, i]`` = :math:`2|\\mu_{x,n}| / \\Delta x_i`
    * ``streaming_y[n, j]`` = :math:`2|\\mu_{y,n}| / \\Delta y_j`

    For future curvilinear geometries, additional curvature terms
    (:math:`\\alpha_n / r_i`) will be stored in ``self.curvature``.

    Parameters
    ----------
    mesh : Mesh1D or Mesh2D
        Base geometry.
    quadrature : AngularQuadrature
        Angular quadrature (Gauss–Legendre, Lebedev, etc.).
    materials : dict mapping material id to Mixture
        Macroscopic cross sections keyed by the integer ids appearing
        in ``mesh.mat_ids`` / ``mesh.mat_map``.  Required (Issue #197
        PR-TYPED-0).  The authoritative source of truth for both
        cross sections and the group count :attr:`ng`; every operator
        that consumes ``sn_mesh`` (L, C, S, F) reads materials from
        here, not from a parallel argument.  All materials must agree
        on ``ng`` — heterogeneous group structures are a
        homogenization-step concern that must precede SNMesh
        construction.

    Attributes
    ----------
    materials : dict mapping material id to Mixture
        The materials dict passed at construction (single source of
        truth).
    ng : int
        Number of energy groups, derived from materials and validated
        for consistency.
    BOUNDARY_OPERATOR_REGISTRY : dict[str, type[BoundaryTraceLaw]]
        Supported boundary-condition kinds (Wave 8 / C8.3). Values
        are :class:`BoundaryTraceLaw` subclasses (``VacuumInflow``,
        ``ReflectiveBoundary``) realized at ``_resolve_one`` time via
        :class:`SNBoundaryRealizer` for every supported mesh
        (1-D Cartesian, 1-D spherical, 1-D cylindrical, 2-D
        Cartesian) and wrapped in :class:`_BoundBoundaryOperator`
        for compatibility with the SN-side call surface.
    bc_left, bc_right : _BoundBoundaryOperator
        Resolved BC operator at the left/right (1-D) boundaries.
        All supported paths: a :class:`_BoundBoundaryOperator` shim
        wrapping the realized 1-arg :class:`LinearOperator` and
        carrying a free-form ``kind`` tag so ``sn_mesh.bc_left ==
        "vacuum"`` style comparisons keep working. Issue #188 /
        C188.3 removed the curvilinear bypass; the realizer is now
        applied uniformly for 1-D Cartesian, 1-D spherical, 1-D
        cylindrical, and 2-D Cartesian meshes.
    bc_xmin, bc_xmax, bc_ymin, bc_ymax : _BoundBoundaryOperator
        Same conventions for the 4 faces of a 2-D mesh. On 1-D
        meshes ``bc_xmin`` / ``bc_xmax`` alias ``bc_left`` /
        ``bc_right``, and ``bc_ymin`` / ``bc_ymax`` are realized
        :class:`ReflectiveBoundary(axis="y")` placeholders routed
        through the realizer with
        :meth:`SNMethodSpace.minimal(quad)`; for GL1D the realized
        op is an identity :class:`PermutationOperator` (no-op).
    """

    BOUNDARY_OPERATOR_REGISTRY: ClassVar[dict[str, type[BoundaryTraceLaw]]] = {
        "vacuum": VacuumInflow,
        "reflective": ReflectiveBoundary,
    }
    # Wave 8 (C8.3): values are the LAW CLASSES themselves (not factory
    # functions), looked up at ``_resolve_one`` time. The pre-Wave-8
    # factory functions (``_sn_vacuum_boundary_operator`` /
    # ``_sn_reflective_boundary_operator``) are gone; their job is now
    # done by ``_resolve_one`` directly, dispatching to
    # :class:`SNBoundaryRealizer` for every supported mesh. Issue #188
    # / C188.3 removed the curvilinear bypass: the realizer is now
    # applied uniformly for 1-D Cartesian, 1-D spherical, 1-D
    # cylindrical, and 2-D Cartesian meshes.
    #
    # The 5 other kinds the realizer handles today (``white``,
    # ``periodic``, ``albedo``, ``prescribed_inflow``, ``mixed``) are
    # NOT registered here — adding them requires SN-sweep-side wiring
    # (sweep cycles for periodic, etc.) that is out of scope for
    # Wave 8. Future expansion is mechanical: add the law class as a
    # value here, ensure the realizer dispatch handles it (it does),
    # and add an SN-side test that the sweep behaves correctly.

    def __init__(
        self,
        mesh: Mesh1D | Mesh2D,
        quadrature: AngularQuadrature,
        materials: "dict[int, Mixture]",
        cell_update: CellUpdate | None = None,
        pole_angular_closure: PoleAngularClosure | None = None,
    ) -> None:
        # Issue #197 PR-TYPED-0: ``materials`` is a REQUIRED parameter.
        # SNMesh IS the SN phase space (mesh × quadrature × material
        # group structure); constructing it without materials would
        # leave ``.ng`` undefined and downstream operators with no
        # cross-section keying.  The single source of truth for the
        # material dict moves from SNSolver to SNMesh: every operator
        # that consumes ``sn_mesh`` (L, C, S, F) reads materials and
        # ``ng`` from here, not from a parallel argument.  See
        # coding-elegance Pattern 4 (illegal states unrepresentable)
        # and Pattern 7 (normalize at definition site).
        self.mesh = mesh
        self.quad = quadrature
        self.materials: "dict[int, Mixture]" = materials
        # R-1 Phase A C1: every SNMesh exposes an axis-tuple view of its
        # dimensionality. The axis tuple is the canonical
        # dim-agnostic ground truth for spatial_shape / face_labels /
        # face_shape / face_outflow_ordinates / n_unknowns_flat (the
        # properties below); legacy ``nx``/``ny``/``dx``/``dy`` survive
        # this layer as deprecated shims and retire in a followup.
        # The adapter handles both 1-D coord systems (Cartesian /
        # spherical / cylindrical) and 2-D Cartesian — see
        # :func:`orpheus.sn.axis.axes_from_legacy_mesh`. From C6 onward
        # every dim-agnostic operator reads its shape from
        # :attr:`SNMesh.axes`, NOT from ``mesh.coord`` /
        # ``mesh.edges_x`` / ``mesh.bc_xmin``.
        self.axes: tuple[Axis1D, ...] = axes_from_legacy_mesh(mesh)
        # Cell-update strategy (Wave D Round 2 Issue #161). Defaults to
        # :class:`DiamondDifference`, which reproduces the existing
        # inlined sweep math bit-identically — every regression snapshot
        # at ``tests/sn/regression/snapshots/`` was generated with DD
        # and continues to match bit-for-bit when the unified sweep
        # dispatches via ``self.cell_update.update(...)``.  Wave C-extension
        # will ship LD / EC / Step strategies; users will then pass
        # ``cell_update=LinearDiscontinuous()`` etc. at construction time.
        self.cell_update: CellUpdate = (
            cell_update if cell_update is not None else DiamondDifference()
        )
        # Issue #168 Phase C retired the Phase A ``BoundaryFaceFlux``
        # Protocol entirely. The sweep-frame apply matvec now uses WDD
        # propagation cell-by-cell with the BC trace law applied at
        # the boundary edge; the Phase A ``DDExtrapolation`` strategy
        # is subsumed by this rewrite. The
        # ``boundary_face_flux`` field is gone from the constructor.
        # Angular-redistribution closure (Issue #168 Phase B + D).
        #
        # Phase D (Issue #168 ERR-026 closure, 2026-05-12) DEFAULT FLIP:
        # the default angular closure is now
        # :class:`MorelMontryAngularSweep` — the canonical Hébert §3.9.4
        # per-cell M-M weighted DD angular recurrence with the Carlson
        # coupled-pole seed (Eqs. 3.432-3.435) supplied by
        # :class:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonInwardSweep`.
        # The Carlson seed closes the M-M recurrence's pole-face
        # initialisation gap, making the per-ordinate flat-flux
        # residual identically zero on sphere Gate 1.1 MMS (per
        # ``tests/sn/test_phase_c_gates.py::test_apply_curvilinear_per_ordinate_flat_flux_residual``).
        # The default flip pairs the canonical pole-angular closure
        # with its canonical spatial-closure partner.
        #
        # PR-TYPED-6.5 + Step 7 ship two angular-closure strategies:
        #
        # * :class:`MorelMontryAngularSweep` (default for curvilinear)
        #   — canonical Hébert §3.9.4 form with Carlson coupled-pole
        #   seed.  Closes ERR-026 on sphere; preserves cylindrical
        #   Gate 1.1 regression-stability.
        # * :class:`IdentityAngularClosure` (default for Cartesian) —
        #   no angular redistribution (flat geometry has no Hébert
        #   §3.9.4 term); ships neutral zeros for the per-cell
        #   contribution interface.
        #
        # Used by the unified matvec
        # (:func:`transport_operator_matvec_unified`) via
        # :meth:`StreamingOperator.apply` (the new Resolution A leaf)
        # and, through Step 7's rewire, also by
        # :meth:`SNStreamingOperator.apply` (the legacy bundle, in
        # transition until ``solve_sn`` migrates to the (L+C) algebra).
        #
        # PR-TYPED-6.5 Phase 2.9: instantiation is deferred until AFTER
        # the ``match mesh.coord:`` block populates ``self.reduced`` /
        # ``self._volumes`` / ``self.dx`` (the data the strategies bind
        # to).  Cartesian gets :class:`IdentityAngularClosure`; sphere
        # and cylinder get :class:`MorelMontryAngularSweep`.  See the
        # ``self.pole_angular_closure = …`` line after ``self._resolve_bcs(...)``
        # below.
        self._user_supplied_closure = pole_angular_closure

        # Normalise to (nx, ny) shaped arrays for both 1-D and 2-D
        if isinstance(mesh, Mesh1D):
            self.nx: int = mesh.N
            self.ny: int = 1
            self.dx: np.ndarray = mesh.widths
            self.dy: np.ndarray = np.array([1.0])
            self.mat_map: np.ndarray = mesh.mat_ids.reshape(mesh.N, 1)
            self._volumes: np.ndarray = mesh.volumes.reshape(mesh.N, 1)
            self._areas: np.ndarray | None = mesh.areas
        else:
            self.nx = mesh.nx
            self.ny = mesh.ny
            self.dx = mesh.dx
            self.dy = mesh.dy
            self.mat_map = mesh.mat_map
            self._volumes = mesh.volumes
            # 2-D mesh has per-face areas of a different shape; not
            # consumed by today's matvec callers — leave None.
            self._areas = None

        # Dispatch stencil setup by coordinate system.
        #
        # Curvilinear connection-coefficient math (sphere / cylinder) lives
        # in :mod:`orpheus.geometry.reduced_operator` (Wave B Issue 6) so
        # MoC and CP can consume the same primitive — Cardinal Rule 2
        # forbids duplicating it on each solver-side mesh class.  Cartesian
        # streaming_x / streaming_y stencils are SN-specific (DD denominator
        # precomputation) and stay local to ``_setup_cartesian``.
        #
        # ``self.reduced`` is the canonical accessor every downstream
        # consumer should bind to: ``sn_mesh.reduced.streaming_terms(
        # cell_idx, dir_idx, mu_level_idx)`` returns the per-(cell,
        # direction) packet a sweep cell update needs.  Wave D Round 2
        # (Issue 12) and Wave E migrate the 6 production read sites
        # (sweep.py + solver.py) to read through it.  Until then, the
        # backward-compat ``@property`` accessors below preserve the
        # legacy attribute names with a ``DeprecationWarning``.
        match mesh.coord:
            case CoordSystem.CARTESIAN:
                self._setup_cartesian()
                # Slab also gets a ``ReducedStreamingOperator`` for
                # completeness so ``sn_mesh.reduced`` is always populated;
                # the slab variant carries empty curvature arrays and
                # ``requires_upstream_angular_state = False``.
                if isinstance(mesh, Mesh1D):
                    self.reduced: ReducedStreamingOperator = slab_streaming(
                        mesh, quadrature,
                    )
                else:
                    self.reduced = None  # type: ignore[assignment]
            case CoordSystem.CYLINDRICAL:
                assert isinstance(mesh, Mesh1D)
                self.reduced = cylindrical_streaming(mesh, quadrature)
                self.curvature: str | None = "cylindrical"
                # 2-D Cartesian-style streaming arrays not used here.
                self.streaming_x: np.ndarray | None = None
                self.streaming_y: np.ndarray | None = None
                # Wave 2 sweep-DAG precompute is 2-D Cartesian only;
                # curvilinear sweeps walk the cell graph differently
                # (see ``dag_walk``). MoC will define its own
                # primitive (fiber bundles + solution sheaves).
                self.sweep_graphs: dict[OctantLabel, SweepDependencyGraph] | None = None
            case CoordSystem.SPHERICAL:
                assert isinstance(mesh, Mesh1D)
                self.reduced = spherical_streaming(mesh, quadrature)
                self.curvature = "spherical"
                self.streaming_x = None
                self.streaming_y = None
                self.sweep_graphs = None

        # Resolve boundary conditions from mesh declarations
        self._resolve_bcs(mesh)

        # ── Materials consistency validation (Issue #197 PR-TYPED-0) ──
        # Two checks at construction time:
        #
        #   (1) every material id used in ``mat_map`` must have an
        #       entry in ``materials`` — otherwise downstream code
        #       would key into an undefined material.
        #   (2) all materials must agree on ``ng`` — SN requires one
        #       uniform group structure; heterogeneous ``ng`` is a
        #       homogenization-step concern that must precede SNMesh
        #       construction, not a silent runtime issue.
        #
        # Both checks surface at SNMesh construction, NOT lazily —
        # the failure mode (operators built on a bad SNMesh) is
        # action-at-a-distance otherwise.
        self._validate_materials()
        # Trigger ``ng`` property's consistency check eagerly so
        # mismatched-ng materials raise at construction time.
        _ = self.ng

        # ── Pole-angular closure binding (PR-TYPED-6.5 Phase 2.9) ──
        # All upstream state needed by the closure constructors is now
        # available (``self.reduced``, ``self._volumes``, ``self.dx``,
        # ``self.quad``, ``self.ng``).  If the user supplied a closure
        # at construction, use it verbatim; otherwise instantiate the
        # default-by-coord-system bound to ``self``.
        if self._user_supplied_closure is not None:
            self.pole_angular_closure: PoleAngularClosure = (
                self._user_supplied_closure
            )
        else:
            closure_cls = default_angular_closure_class(mesh.coord)
            self.pole_angular_closure = closure_cls(self)
        # Drop the temporary attribute now that the closure is bound.
        del self._user_supplied_closure

    # ── Boundary condition resolution ─────────────────────────────────

    def _resolve_bcs(self, mesh: Mesh1D | Mesh2D) -> None:
        r"""Resolve geometry-declared BCs into Wave-8 realized operators.

        ``None`` on the mesh defaults to ``BC("reflective")`` (infinite
        lattice / eigenvalue convention).

        Wave 8 (C8.3) + Issue #188 / C188.3: each face attribute
        carries a :class:`_BoundBoundaryOperator` shim wrapping the
        1-arg :class:`LinearOperator` produced by
        :class:`SNBoundaryRealizer`. The realizer is applied uniformly
        for 1-D Cartesian, 1-D spherical, 1-D cylindrical, and 2-D
        Cartesian meshes (C188.1 + C188.2 lifted the Mesh1D
        curvilinear guard on
        :meth:`InflowTraceSpace.from_mesh_and_quadrature`; C188.3
        removes the matching curvilinear bypass here).
        """
        default = BC("reflective")

        # Build the inflow / outflow trace spaces ONCE per SNMesh.
        # Issue #188 / C188.1+C188.2 extended
        # :meth:`InflowTraceSpace.from_mesh_and_quadrature` to all
        # 1-D coord systems (spherical / cylindrical share the
        # ``("left", "right")`` face structure with slab). 2-D
        # cylindrical :class:`Mesh2D` is the only mesh that still
        # raises (no 2-D cylindrical SN sweep exists in ORPHEUS
        # today); it stays at ``_inflow_trace = None``.
        self._inflow_trace = None
        self._outflow_trace = None
        build_trace = (
            isinstance(mesh, Mesh1D)
            or (isinstance(mesh, Mesh2D) and mesh.coord == CoordSystem.CARTESIAN)
        )
        if build_trace:
            from orpheus.numerics.trace_space import (
                InflowTraceSpace,
                OutflowTraceSpace,
            )
            faces = self._face_names_for_mesh(mesh)
            self._inflow_trace = InflowTraceSpace.from_mesh_and_quadrature(
                mesh, self.quad, faces=faces,
            )
            self._outflow_trace = OutflowTraceSpace.from_mesh_and_quadrature(
                mesh, self.quad, faces=faces,
            )

        if isinstance(mesh, Mesh1D):
            self.bc_left = self._resolve_one(mesh.bc_left or default, "left")
            self.bc_right = self._resolve_one(mesh.bc_right or default, "right")
            # Expose 2-D-style attributes for uniform sweep access.
            # The ``y`` faces of a 1-D mesh are degenerate (no y
            # dimension); 1-D sweeps don't consume them. The 1-D
            # trace space's ``face_names == ("left", "right")``, so
            # ``SNMethodSpace.for_face(face="ymin", inflow_trace=...)``
            # would raise — we route through the realizer with
            # :meth:`SNMethodSpace.minimal(self.quad)` instead. For
            # :class:`GaussLegendre1D` (the 1-D quadrature)
            # ``reflection_index("y")`` returns the identity
            # permutation (every ordinate is its own partner because
            # ``mu_y == 0``), so the realized op is a no-op
            # :class:`PermutationOperator`. The realizer's
            # :class:`ReflectiveBoundary` branch does NOT read
            # ``inflow_indices``, only ``law.axis`` and
            # ``quad.reflection_index``, so passing the minimal
            # method space is safe.
            #
            # Any future consumer that calls ``apply`` on these
            # placeholders sees the same uniform 1-arg contract as
            # every other ``bc_*`` attribute; the production 2-D
            # Cartesian path at :func:`solution_to_angular_flux` does
            # consume ``bc_ymin.apply(...)`` on 1-D meshes when
            # ``curvature is None`` but discards the result because
            # the inner ``if mu_y[n] > 1e-15`` filter is false for
            # every ordinate on GL1D.
            self.bc_xmin = self.bc_left
            self.bc_xmax = self.bc_right
            y_method_space = SNMethodSpace.minimal(self.quad)
            y_realized = SNBoundaryRealizer().realize(
                ReflectiveBoundary(axis="y", albedo=1.0),
                y_method_space,
            )
            self.bc_ymin = _BoundBoundaryOperator(
                y_realized, kind="reflective",
            )
            self.bc_ymax = _BoundBoundaryOperator(
                y_realized, kind="reflective",
            )
        else:
            self.bc_xmin = self._resolve_one(mesh.bc_xmin or default, "xmin")
            self.bc_xmax = self._resolve_one(mesh.bc_xmax or default, "xmax")
            self.bc_ymin = self._resolve_one(mesh.bc_ymin or default, "ymin")
            self.bc_ymax = self._resolve_one(mesh.bc_ymax or default, "ymax")
            self.bc_left = self.bc_xmin
            self.bc_right = self.bc_xmax

    @staticmethod
    def _face_names_for_mesh(mesh: Mesh1D | Mesh2D) -> tuple[str, ...]:
        """Ordered face-name tuple matching the trace-space row order."""
        if isinstance(mesh, Mesh1D):
            return ("left", "right")
        return ("xmin", "xmax", "ymin", "ymax")

    def _resolve_one(self, bc: BC, face: str):
        r"""Realize a BC on the given face.

        Build an :class:`SNMethodSpace` carrying the precomputed
        inflow / outflow traces, hand it to
        :class:`SNBoundaryRealizer.realize`, wrap the 1-arg result in
        :class:`_BoundBoundaryOperator` so the SN-side call surface
        sees a uniform 1-arg ``apply(psi)`` contract.

        Issue #188 / C188.3: every supported mesh (1-D Cartesian,
        1-D spherical, 1-D cylindrical, 2-D Cartesian) routes
        through the realizer here. The pre-C188.3 curvilinear
        bypass — which wrapped the raw 2-arg
        :class:`BoundaryTraceLaw` with a bound quadrature — is
        gone, made redundant by C188.1+C188.2's curvilinear
        :class:`InflowTraceSpace` support.
        """
        law_cls = self.BOUNDARY_OPERATOR_REGISTRY.get(bc.kind)
        if law_cls is None:
            supported = ", ".join(
                f"'{k}'" for k in sorted(self.BOUNDARY_OPERATOR_REGISTRY)
            )
            raise ValueError(
                f"SN solver does not support boundary condition '{bc.kind}' "
                f"on face '{face}'. Supported: {supported}."
            )
        # Construct the law instance with face-appropriate axis for
        # reflective; the others take no kwargs.
        if law_cls is ReflectiveBoundary:
            axis = "y" if face in ("ymin", "ymax") else "x"
            law = law_cls(axis=axis, albedo=1.0)
        else:
            law = law_cls()

        method_space = SNMethodSpace.for_face(
            mesh=self.mesh,
            quadrature=self.quad,
            face=face,
            inflow_trace=self._inflow_trace,
            outflow_trace=self._outflow_trace,
        )
        realized = SNBoundaryRealizer().realize(law, method_space)
        return _BoundBoundaryOperator(realized, kind=bc.kind)

    # ── Materials validation ──────────────────────────────────────────

    def _validate_materials(self) -> None:
        """Validate the materials dict against the mat_map.

        Every material id referenced in ``self.mat_map`` MUST appear
        as a key in ``self.materials``.  Failure surfaces here at
        construction time, not lazily inside a solver step.

        Raises
        ------
        ValueError
            If any ``mat_map`` id is missing from ``materials``; the
            error message shows both sets so the user can see the gap.
        """
        if not self.materials:
            raise ValueError(
                "SNMesh requires a non-empty materials dict; got "
                f"materials={self.materials!r}."
            )
        used_ids = set(int(x) for x in np.unique(self.mat_map))
        available_ids = set(self.materials.keys())
        missing = used_ids - available_ids
        if missing:
            raise ValueError(
                f"SNMesh.mat_map references material ids "
                f"{sorted(missing)} that are NOT in materials "
                f"(available ids: {sorted(available_ids)}; used ids: "
                f"{sorted(used_ids)})."
            )

    # ── Properties ────────────────────────────────────────────────────

    @property
    def ng(self) -> int:
        """Number of energy groups; uniform across all materials.

        Derived from ``self.materials``; the single source of truth
        for the group count.  All materials must agree on ``ng`` —
        SN requires one uniform group structure across the mesh.

        Raises
        ------
        InconsistentMaterialsError
            If materials disagree on ``ng``.  A homogenization step
            must precede SNMesh construction in that case.
        ValueError
            If ``self.materials`` is empty (caught at construction
            time by ``_validate_materials``).
        """
        if not self.materials:
            raise ValueError(
                "SNMesh.ng undefined — no materials.  Construct "
                "SNMesh(mesh, quadrature, materials) with a non-empty "
                "materials dict."
            )
        ngs = {m.ng for m in self.materials.values()}
        if len(ngs) != 1:
            raise InconsistentMaterialsError(
                f"SNMesh requires uniform ng across all materials; got "
                f"ng values {sorted(ngs)} in materials dict with keys "
                f"{sorted(self.materials.keys())}.  Homogenize to a "
                f"common group structure before SNMesh construction."
            )
        return ngs.pop()

    @property
    def volumes(self) -> np.ndarray:
        """Cell volumes, shape (nx, ny)."""
        return self._volumes

    @property
    def areas(self) -> np.ndarray:
        """Face areas at each radial edge, shape (nx+1,) (1-D meshes).

        Sourced from :attr:`Mesh1D.areas` (computed eagerly at mesh
        construction via :func:`orpheus.geometry.coord.compute_areas_1d`).
        Cartesian slab returns an array of ones; cylinder returns
        :math:`2\\pi r`; sphere returns :math:`4\\pi r^2`.

        Raises
        ------
        AttributeError
            If accessed on a 2-D mesh (face areas have a different shape
            and are not consumed by today's matvec callers).
        """
        if self._areas is None:
            raise AttributeError(
                "SNMesh.areas is not defined for 2-D meshes; "
                "face-area data lives in the underlying Mesh2D directly."
            )
        return self._areas

    @property
    def is_1d(self) -> bool:
        """True if this is a 1-D mesh (ny == 1)."""
        return self.ny == 1

    # ── Dim-agnostic geometry primitives (R-1 Phase A C1) ─────────────

    @property
    def ndim(self) -> int:
        """Number of spatial dimensions; equals ``len(self.axes)``."""
        return len(self.axes)

    @property
    def spatial_shape(self) -> tuple[int, ...]:
        r"""Per-axis cell counts ``(n_0, n_1, ...)``.

        The canonical dim-agnostic shape descriptor. For 1-D meshes
        returns ``(nx,)``; for 2-D Cartesian returns ``(nx, ny)``;
        for 3-D (followup) would return ``(nx, ny, nz)``.

        Every dim-agnostic shape reader (typed-field factories, pack
        convention, sweep DAG) reads from here, NOT from ``self.nx``
        / ``self.ny`` (the legacy shims that retire in a followup).
        """
        return _axis_spatial_shape(self.axes)

    @property
    def face_labels(self) -> tuple[FaceLabel, ...]:
        r"""Canonical boundary-face inventory.

        Each :class:`~orpheus.sn.axis.FaceLabel` carries an
        ``axis_index`` and an ``endpoint`` label, derived from the
        per-axis endpoints. Cartesian 1-D returns 2 labels; spherical
        / cylindrical 1-D returns 1 label (the pole is NOT a face —
        see :class:`~orpheus.sn.axis.Axis1D` docstring); 2-D Cartesian
        returns 4 labels; synthetic 3-D Cartesian would return 6.

        The iteration order is the canonical concatenation order for
        :meth:`AngularFlux.to_flat` (C3) and the canonical iteration
        order for :attr:`BoundaryFlux.face_buffers` (C4).
        """
        return _axis_face_labels(self.axes)

    def face_shape(self, label: FaceLabel) -> tuple[int, ...]:
        r"""Spatial shape of the boundary face identified by ``label``.

        The face lies in the codimension-1 hyperplane spanned by the
        axes other than ``label.axis_index``; its shape is the
        per-axis cell count of those axes in axis-index order.
        """
        return _axis_face_shape(self.axes, label)

    def face_outflow_ordinates(self, label: FaceLabel) -> np.ndarray:
        r"""Ordinate indices whose direction-cosine is OUTWARD at face ``label``.

        At the ``max`` / ``outer`` endpoint of an axis, outflow is
        :math:`\mu_{axis} > +10^{-15}`; at ``min``, outflow is
        :math:`\mu_{axis} < -10^{-15}`. Ordinates exactly tangent to
        the face contribute neither inflow nor outflow.

        This method is the canonical producer for the per-face
        outflow mask used by the pack convention (C3),
        :class:`BoundaryFlux.face_buffers` (C4), and the sweep DAG
        face-trace state (C5).
        """
        return _axis_face_outflow_ordinates(self.axes, label, self.quad)

    @property
    def n_unknowns_flat(self) -> int:
        r"""Total flat-vector size for typed :class:`AngularFlux`.

        The pack convention (C3) is the direct-sum decomposition
        :math:`V = V_\text{cells} \oplus \bigoplus_\ell V_{\text{face}, \ell}`;
        ``n_unknowns_flat`` is the dimension of that vector space.
        Cells contribute :math:`N \cdot n_g \cdot \prod_i n_i`; each
        face :math:`\ell` contributes
        :math:`|\text{outflow}_\ell| \cdot n_g \cdot \prod_{i \ne \text{axis}(\ell)} n_i`.
        """
        return _axis_n_unknowns_flat(self.axes, self.quad, self.ng)

    @classmethod
    def from_axes(
        cls,
        axes: tuple[Axis1D, ...],
        quadrature: "AngularQuadrature",
        materials: "dict[int, Mixture]",
        *,
        mat_map: np.ndarray | None = None,
        cell_update: CellUpdate | None = None,
        pole_angular_closure: PoleAngularClosure | None = None,
    ) -> "SNMesh":
        r"""Build an :class:`SNMesh` from an axis tuple.

        The dim-agnostic constructor surface introduced by R-1 Phase A
        C1. The axis tuple is round-tripped through a legacy
        :class:`Mesh1D` / :class:`Mesh2D` so the existing constructor
        body (BC realization, streaming stencil, materials validation,
        pole-angular closure binding) is reached uniformly regardless
        of which surface the caller used. Layer A keeps ``SNMesh.mesh``
        as a :class:`Mesh1D` / :class:`Mesh2D`; the round-trip retires
        with the legacy mesh dataclass in a later phase.

        Parameters
        ----------
        axes : tuple of :class:`~orpheus.sn.axis.Axis1D`
            Per-axis 1-D mesh descriptors. Length 1 → 1-D mesh;
            length 2 → 2-D Cartesian mesh. Length ≥3 is NOT supported
            in C1 — no ``Mesh3D`` dataclass exists today (D9 of the
            ultraplan; the 3-D admission gate exercises the pure
            shape functions in :mod:`orpheus.sn.axis` directly).
        quadrature : :class:`AngularQuadrature`
            Angular quadrature.
        materials : dict[int, Mixture]
            Materials dict keyed by material id; same contract as the
            legacy constructor.
        mat_map : ndarray or None
            Material-id assignment. Shape ``spatial_shape``. Defaults
            to all-zeros (single material with id 0).
        cell_update : CellUpdate or None
            Cell-update strategy. Defaults to :class:`DiamondDifference`.
        pole_angular_closure : PoleAngularClosure or None
            Override the default pole-angular closure
            (curvilinear → :class:`MorelMontryAngularSweep`,
            Cartesian → :class:`IdentityAngularClosure`).
        """
        mesh = legacy_mesh_from_axes(axes, mat_map=mat_map)
        return cls(
            mesh=mesh,
            quadrature=quadrature,
            materials=materials,
            cell_update=cell_update,
            pole_angular_closure=pole_angular_closure,
        )

    def material_xs_field(self) -> "MaterialXSField":
        """Build the macroscopic XS field from this mesh's materials.

        Issue #197 PR-TYPED-1.  Returns a :class:`MaterialXSField`
        wrapping the per-material :class:`Mixture` data plus this
        mesh's ``mat_map`` — the single source of truth for both
        per-cell and per-material XS access used by every SN operator
        (L, C, S, F).

        Lazy import of :mod:`.material_xs_field` to avoid a circular
        dependency (the module imports :class:`SNMesh` via
        ``TYPE_CHECKING``).
        """
        from .material_xs_field import MaterialXSField
        return MaterialXSField.from_mesh(self)

    # ── Typed-field factories (Issue #197 PR-TYPED-2) ─────────────────

    def zeros_angular_flux(self) -> "AngularFlux":
        r"""Build a zero :class:`AngularFlux` sized to this mesh.

        Returns an :class:`~orpheus.sn.angular_flux.AngularFlux` of
        shape ``(N, ng, nx, ny)`` filled with zeros.  Use as the
        ``angular`` initial guess in inner-loop iterations.
        """
        from .angular_flux import AngularFlux
        return AngularFlux(
            np.zeros((self.quad.N, self.ng, self.nx, self.ny)),
            self,
        )

    def zeros_scalar_flux(self) -> "ScalarFlux":
        r"""Build a zero :class:`ScalarFlux` sized to this mesh.

        Returns a :class:`~orpheus.sn.scalar_flux.ScalarFlux` of shape
        ``(ng, nx, ny)`` filled with zeros.  Use as the ``phi``
        initial guess in inner-loop iterations or as the zero
        isotropic source on the first sweep.
        """
        from .scalar_flux import ScalarFlux
        return ScalarFlux(
            np.zeros((self.ng, self.nx, self.ny)), self,
        )

    def zeros_boundary_flux(self) -> "BoundaryFlux":
        r"""Build a zero :class:`BoundaryFlux` sized to this mesh.

        Delegates to :meth:`BoundaryFlux.zeros`; allocates only the
        face / persistent buffers the mesh's geometry consumes (slab
        gets two 1-D faces, curvilinear gets one outer face,
        2-D Cartesian gets the two persistent ``(N, ng, nx+1, ny)``
        buffers).
        """
        from .boundary_flux import BoundaryFlux
        return BoundaryFlux.zeros(self)

    def zeros_isotropic_source(self) -> "IsotropicSource":
        r"""Build a zero :class:`IsotropicSource` sized to this mesh.

        Returns an :class:`~orpheus.sn.sources.IsotropicSource` of
        shape ``(ng, nx, ny)`` filled with zeros.  The principled
        starting point for source-iteration accumulation: P0 in-scatter,
        (n,2n), and fission contributions all add into this buffer.
        """
        from .sources import IsotropicSource
        return IsotropicSource(
            np.zeros((self.ng, self.nx, self.ny)), self,
        )

    def zeros_per_ordinate_source(self) -> "PerOrdinateSource":
        r"""Build a zero :class:`PerOrdinateSource` sized to this mesh.

        Returns an :class:`~orpheus.sn.sources.PerOrdinateSource` of
        shape ``(N, ng, nx, ny)`` filled with zeros.  The principled
        starting point for the :math:`P_\ell \ge 1` Galerkin
        reconstruction + MMS external-source accumulation buffer.
        """
        from .sources import PerOrdinateSource
        return PerOrdinateSource(
            np.zeros((self.quad.N, self.ng, self.nx, self.ny)), self,
        )

    def zeros_harmonic_moments(self, L: int) -> "HarmonicMomentField":
        r"""Build a zero :class:`HarmonicMomentField` at order ``L``.

        Issue #197 PR-TYPED-4 — companion to
        :meth:`zeros_isotropic_source` / :meth:`zeros_per_ordinate_source`
        for the moment-space carrier consumed by the
        :math:`R \cdot \Lambda \cdot M` Galerkin pipeline.

        Parameters
        ----------
        L : int
            Maximum harmonic order; result has shape
            ``(L+1, 2L+1, ng, nx, ny)``.
        """
        from .harmonic_moment_field import HarmonicMomentField
        values = np.zeros(
            (L + 1, 2 * L + 1, self.ng, self.nx, self.ny),
        )
        return HarmonicMomentField(values, self, L)

    # ── Sweep DAG traversal ───────────────────────────────────────────

    _DEGENERATE_ABS_ETA_THRESHOLD: ClassVar[float] = 1e-15

    def dag_walk(
        self,
        *,
        ordinate_idx: int | None = None,
        direction_sign: int | None = None,
        mu_level_idx: int | None = None,
    ) -> Iterator[CellVisit]:
        r"""Walk the per-ordinate cell DAG in topological order.

        Issue #196 Phase G Step 2.6 (Q3): the single canonical iteration
        primitive for 1-D sweeps.  Yields visits either for a single
        ordinate or for all ordinates of a sweep direction under one
        XOR signature.

        The SN sweep on a given ordinate is forward substitution on
        the block-triangular streaming + collision operator under the
        ordinate's DAG ordering.  This method yields the per-cell
        visit packets in that DAG order; the consumer folds over the
        packets, threading the spatial-upstream face flux through the
        accumulator and writing the per-cell angular state into a
        persistent array.

        Exactly one of ``ordinate_idx`` or ``direction_sign`` must be
        supplied (XOR):

        * ``ordinate_idx=n`` — yields visits for a single ordinate.
          For slab/sphere: ``n`` is the global ordinate index.  For
          cylindrical: ``n`` is the within-level azimuthal index
          :math:`m \in [0, M)` and ``mu_level_idx`` MUST also be
          supplied; the signed :math:`\eta` resolves through
          ``quad.level_indices[mu_level_idx][n]``.
        * ``direction_sign=±1`` — yields visits for the sweep
          direction (``+1`` outward, ``-1`` inward).  Cell ordering
          depends ONLY on the direction sign (and level for
          cylindrical), so any ordinate in the correct sign class
          yields the same cell sequence; this branch picks a
          non-degenerate representative.

        SN-specific by design.  MoC will not consume this method —
        its mathematical structure is fiber bundles + solution
        sheaves, a different DAG shape.  Premature abstraction
        avoided per Cardinal Rule 2.

        Parameters
        ----------
        ordinate_idx : int | None
            See above; mutually exclusive with ``direction_sign``.
        direction_sign : int | None
            See above; mutually exclusive with ``ordinate_idx``.
        mu_level_idx : int | None
            For cylindrical geometry: which :math:`\mu`-level the
            ordinate (subset) belongs to.  ``None`` for slab/sphere;
            required for cylindrical.

        Yields
        ------
        CellVisit
            One per cell, in topological order.  The packet's
            :attr:`face_area_downstream` is float: ``1.0`` for slab,
            ``0.0`` for cylindrical pure-azimuthal degenerate,
            physical face area for sphere / non-degenerate cylinder.

        Raises
        ------
        ValueError
            If neither or both of ``ordinate_idx`` / ``direction_sign``
            are supplied; if ``direction_sign not in (+1, -1)``; if
            called on a 2-D Cartesian mesh (no
            :class:`ReducedStreamingOperator`); if a cylindrical
            mesh is queried without ``mu_level_idx``; or if no
            non-degenerate representative ordinate exists for
            ``direction_sign``.

        Notes
        -----
        2-D Cartesian wavefront scheduling is intentionally not
        encapsulated here — its anti-diagonal vectorisation
        operates on cell slices, not per-cell visits.
        """
        if (ordinate_idx is None) == (direction_sign is None):
            raise ValueError(
                "dag_walk requires exactly one of `ordinate_idx` or "
                "`direction_sign`."
            )
        if self.reduced is None:
            raise ValueError(
                "dag_walk is only defined for meshes with a "
                "ReducedStreamingOperator (1-D Cartesian, spherical, "
                "or cylindrical).  2-D Cartesian wavefront sweeps "
                "use anti-diagonal scheduling, not per-cell visits."
            )
        coord = self.reduced.coord
        if coord is CoordSystem.CYLINDRICAL and mu_level_idx is None:
            raise ValueError(
                "cylindrical dag_walk requires mu_level_idx."
            )

        # Direction-keyed branch: resolve a non-degenerate
        # representative ordinate, then delegate to the
        # ordinate-keyed branch (single source of truth — Pattern 2).
        if direction_sign is not None:
            if direction_sign not in (+1, -1):
                raise ValueError(
                    f"direction_sign must be +1 or -1; got "
                    f"{direction_sign}"
                )
            ordinate_idx = self._representative_ordinate(
                direction_sign, mu_level_idx,
            )

        # Ordinate-keyed branch.
        if coord is CoordSystem.CARTESIAN:
            yield from self._iter_cartesian_visits(ordinate_idx)
            return
        if coord is CoordSystem.SPHERICAL:
            yield from self._iter_spherical_visits(ordinate_idx)
            return
        if coord is CoordSystem.CYLINDRICAL:
            yield from self._iter_cylindrical_visits(
                ordinate_idx, mu_level_idx,
            )
            return
        raise ValueError(  # pragma: no cover — exhaustive match above
            f"Unknown coord system: {coord!r}"
        )

    def dag_walk_cell_indices(
        self,
        *,
        direction_sign: int,
        mu_level_idx: int | None = None,
    ) -> Iterator[int]:
        r"""Lightweight twin of :meth:`dag_walk` — yields just cell indices.

        Consumers that build their own per-cell algebra from primitives
        (the unified matvec ``transport_operator_matvec_unified``) only
        need the cell traversal order, not the full
        :class:`~orpheus.sn.spatial.cell_update.CellVisit` packet.

        Eliminates per-cell-per-call ``ReducedStreamingOperator.streaming_terms()``
        construction + frozen-dataclass overhead.  PR-TYPED-6c profiling
        showed this was ~14% of matvec time on slab, ~18% on cylinder
        — all building a packet the matvec discards.

        Cell-iteration order matches :meth:`dag_walk`:

        * Slab, sphere, cylinder non-degenerate: ``range(nx)`` for
          :math:`\mu_n \ge 0`, ``range(nx-1, -1, -1)`` for :math:`\mu_n < 0`.
        * Cylindrical pure-azimuthal degenerate
          (:math:`|\eta_n| < 10^{-15}`): ``range(nx)`` regardless of
          ``direction_sign`` — same as :meth:`dag_walk`.
        """
        if self.reduced is None:
            raise ValueError(
                "dag_walk_cell_indices is only defined for meshes with a "
                "ReducedStreamingOperator (1-D Cartesian, spherical, "
                "or cylindrical)."
            )
        coord = self.reduced.coord
        if coord is CoordSystem.CYLINDRICAL and mu_level_idx is None:
            raise ValueError(
                "cylindrical dag_walk_cell_indices requires mu_level_idx."
            )
        if direction_sign not in (+1, -1):
            raise ValueError(
                f"direction_sign must be +1 or -1; got {direction_sign}"
            )

        # Resolve the representative ordinate's signed primary cosine.
        ordinate_idx = self._representative_ordinate(
            direction_sign, mu_level_idx,
        )
        if coord is CoordSystem.CYLINDRICAL:
            level_indices = self.quad.level_indices  # type: ignore[attr-defined]
            global_n = int(level_indices[mu_level_idx][ordinate_idx])
            mu_n = float(self.quad.mu_x[global_n])
        else:
            mu_n = float(self.quad.mu_x[ordinate_idx])

        # Cylindrical degenerate ordinates iterate forward regardless of sign.
        if (
            coord is CoordSystem.CYLINDRICAL
            and abs(mu_n) < self._DEGENERATE_ABS_ETA_THRESHOLD
        ):
            yield from range(self.nx)
            return

        if mu_n >= 0:
            yield from range(self.nx)
        else:
            yield from range(self.nx - 1, -1, -1)

    def _representative_ordinate(
        self,
        direction_sign: int,
        mu_level_idx: int | None,
    ) -> int:
        """Pick a non-degenerate ordinate matching the direction sign.

        Cell ordering in :meth:`dag_walk` depends only on
        ``direction_sign`` (and the level for cylindrical), so any
        non-degenerate ordinate in the correct sign class produces
        the same cell sequence.  The degenerate :math:`|\\eta| <
        10^{-15}` ordinates are excluded because they iterate forward
        regardless of sign and would not match the bulk direction's
        signed iteration.
        """
        assert self.reduced is not None
        coord = self.reduced.coord
        eps = self._DEGENERATE_ABS_ETA_THRESHOLD
        if coord is CoordSystem.CYLINDRICAL:
            level_indices = self.quad.level_indices  # type: ignore[attr-defined]
            level_ords = np.asarray(level_indices[mu_level_idx])
            eta_at_level = self.quad.eta[level_ords]
            if direction_sign == +1:
                cand = np.where(eta_at_level > +eps)[0]
            else:
                cand = np.where(eta_at_level < -eps)[0]
            if cand.size == 0:
                raise ValueError(
                    f"No non-degenerate ordinate in cylindrical level "
                    f"{mu_level_idx} satisfies "
                    f"direction_sign={direction_sign}."
                )
            return int(cand[0])
        mu_x = self.quad.mu_x
        if direction_sign == +1:
            cand = np.where(mu_x > +eps)[0]
        else:
            cand = np.where(mu_x < -eps)[0]
        if cand.size == 0:
            raise ValueError(
                f"No non-degenerate ordinate satisfies "
                f"direction_sign={direction_sign} in this quadrature."
            )
        return int(cand[0])

    def _iter_cartesian_visits(
        self,
        ordinate_idx: int,
    ) -> Iterator[CellVisit]:
        """Yield slab (1-D Cartesian) visits in sweep direction.

        Order: forward (cell 0 → nx-1) for :math:`\\mu \\ge 0`,
        backward for :math:`\\mu < 0`.  Slab carries
        ``face_area_downstream = 1.0`` (neutral curvature; Issue
        #196 Phase G Step 2.5) so the unified cell-balance helper
        consumes one geometry-blind number.
        """
        assert self.reduced is not None
        mu_n = float(self.quad.mu_x[ordinate_idx])
        cell_indices = (
            range(self.nx) if mu_n >= 0 else range(self.nx - 1, -1, -1)
        )
        for i in cell_indices:
            st = self.reduced.streaming_terms(
                cell_idx=i, direction_idx=ordinate_idx,
            )
            yield CellVisit(
                cell_idx=i,
                streaming_terms=st,
                face_area_downstream=1.0,
            )

    def _iter_spherical_visits(
        self,
        ordinate_idx: int,
    ) -> Iterator[CellVisit]:
        """Yield spherical visits in sweep direction.

        Outward (:math:`\\mu \\ge 0`): cell 0 → nx-1, downstream face
        is the outer face ``A[i+1]``.  Inward (:math:`\\mu < 0`):
        cell nx-1 → 0, downstream face is the inner face ``A[i]``.
        """
        assert self.reduced is not None
        mu_n = float(self.quad.mu_x[ordinate_idx])
        if mu_n >= 0:
            cell_indices = range(self.nx)
            select_outer = True
        else:
            cell_indices = range(self.nx - 1, -1, -1)
            select_outer = False
        for i in cell_indices:
            st = self.reduced.streaming_terms(
                cell_idx=i, direction_idx=ordinate_idx,
            )
            face_downstream = (
                st.face_area_outer if select_outer else st.face_area_inner
            )
            yield CellVisit(
                cell_idx=i,
                streaming_terms=st,
                face_area_downstream=face_downstream,
            )

    def _iter_cylindrical_visits(
        self,
        ordinate_idx: int,
        mu_level_idx: int,
    ) -> Iterator[CellVisit]:
        """Yield cylindrical visits in sweep direction for one level.

        ``ordinate_idx`` is the within-level azimuthal index
        :math:`m \\in [0, M)`.  The global ordinate is resolved via
        ``quad.level_indices[mu_level_idx][ordinate_idx]``.

        * :math:`\\eta_n \\ge 0` outward: cell 0 → nx-1, downstream
          face is the outer face.
        * :math:`\\eta_n < 0` inward: cell nx-1 → 0, downstream
          face is the inner face.
        * :math:`|\\eta_n| < 10^{-15}` pure-azimuthal degenerate:
          forward iteration (so the angular M-M closure runs in a
          natural order) but ``face_area_downstream`` is ``None`` —
          no spatial face flow.
        """
        assert self.reduced is not None
        level_indices = self.quad.level_indices  # type: ignore[attr-defined]
        global_n = int(level_indices[mu_level_idx][ordinate_idx])
        eta_n = float(self.quad.eta[global_n])
        abs_eta = abs(eta_n)

        if abs_eta < self._DEGENERATE_ABS_ETA_THRESHOLD:
            # Pure-azimuthal degenerate: no spatial flow.  Iterate
            # forward so the angular M-M closure runs in a natural
            # order; ``face_area_downstream = 0.0`` signals "no
            # spatial flow" to the strategy (geometric truth — the
            # cell has no radial face on this ordinate).  Issue #196
            # Phase G Step 2.5: replaced ``None`` with the
            # geometrically-correct float ``0.0``.
            for i in range(self.nx):
                st = self.reduced.streaming_terms(
                    cell_idx=i,
                    direction_idx=ordinate_idx,
                    mu_level_idx=mu_level_idx,
                )
                yield CellVisit(
                    cell_idx=i,
                    streaming_terms=st,
                    face_area_downstream=0.0,
                )
            return

        if eta_n >= 0:
            cell_indices = range(self.nx)
            select_outer = True
        else:
            cell_indices = range(self.nx - 1, -1, -1)
            select_outer = False
        for i in cell_indices:
            st = self.reduced.streaming_terms(
                cell_idx=i,
                direction_idx=ordinate_idx,
                mu_level_idx=mu_level_idx,
            )
            face_downstream = (
                st.face_area_outer if select_outer else st.face_area_inner
            )
            yield CellVisit(
                cell_idx=i,
                streaming_terms=st,
                face_area_downstream=face_downstream,
            )

    # ── Stencil setup ─────────────────────────────────────────────────

    def _setup_cartesian(self) -> None:
        """Precompute Cartesian diamond-difference streaming coefficients.

        These are the purely geometric parts of the DD denominator:

        .. math::

            \\text{denom} = \\Sigma_t + \\frac{2|\\mu_x|}{\\Delta x}
                            + \\frac{2|\\mu_y|}{\\Delta y}

        Precomputing avoids per-ordinate per-cell divisions in the
        inner sweep loop.
        """
        mu_x = self.quad.mu_x
        mu_y = self.quad.mu_y

        # streaming_x[n, i] = 2|μ_x[n]| / dx[i] — shape (N_ord, nx)
        self.streaming_x: np.ndarray = (
            2.0 * np.abs(mu_x)[:, None] / self.dx[None, :]
        )
        # streaming_y[n, j] = 2|μ_y[n]| / dy[j] — shape (N_ord, ny)
        self.streaming_y: np.ndarray = (
            2.0 * np.abs(mu_y)[:, None] / self.dy[None, :]
        )

        # Curvature terms (None for Cartesian — placeholder for curvilinear)
        self.curvature = None

        # Per-octant sweep dependency graphs (Wave 2 / C2.4) — one
        # ``SweepDependencyGraph`` per in-plane octant ``(sign_x,
        # sign_y) ∈ {-1, +1}²`` (the four streaming octants).
        # Pure-z ordinates with ``sign_x == 0 == sign_y`` are handled
        # by the wavefront sweep's ``Q/Σ_t`` short-circuit and have
        # no entry here. The graphs depend only on cell topology +
        # octant sign convention — independent of fluxes / sources /
        # iteration state — so they are precomputed once at mesh
        # construction and reused across every sweep call.
        self.sweep_graphs: dict[OctantLabel, SweepDependencyGraph] = {
            OctantLabel(sx, sy): SweepDependencyGraph.from_cartesian_2d(
                nx=self.nx, ny=self.ny, label=OctantLabel(sx, sy),
            )
            for sx in (-1, +1)
            for sy in (-1, +1)
        }

    # ── Backward-compat property accessors ────────────────────────────
    #
    # These properties route to ``self.reduced`` (the
    # :class:`ReducedStreamingOperator` instance built at construction).
    # The 6 production read sites in ``sweep.py`` and ``solver.py``
    # continue reading via these names; Wave D Round 2 (Issue 12) and
    # Wave E migrate them to ``self.reduced.streaming_terms(...)`` and
    # the deprecated properties are removed.
    #
    # Each property emits a ``DeprecationWarning`` on access so consumers
    # surface the migration path during normal test runs.  No
    # ``filterwarnings`` config in :file:`pyproject.toml` treats
    # ``DeprecationWarning`` as an error, so existing tests are unaffected.

    @property
    def face_areas(self) -> np.ndarray:
        """[Deprecated] Cell face areas. Use ``self.reduced.face_areas`` instead."""
        warnings.warn(
            "SNMesh.face_areas is deprecated; "
            "use SNMesh.reduced.face_areas instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        assert self.reduced.face_areas is not None
        return self.reduced.face_areas

    @property
    def delta_A(self) -> np.ndarray:
        """[Deprecated] Face-area differences. Use ``self.reduced.delta_A`` instead."""
        warnings.warn(
            "SNMesh.delta_A is deprecated; "
            "use SNMesh.reduced.delta_A instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        assert self.reduced.delta_A is not None
        return self.reduced.delta_A

    # The six curvature-specific accessors (``alpha_half``,
    # ``redist_dAw``, ``tau_mm``, and the cylindrical per-level
    # analogues) added in Wave D Round 1.1 as DeprecationWarning
    # shims were retired in Wave E Round 2 (Issue #164) along with
    # the BiCGSTAB FD-operator API surface that was their only
    # consumer.  Use ``self.reduced.<name>`` directly.
