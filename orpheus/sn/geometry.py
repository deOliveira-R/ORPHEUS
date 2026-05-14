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
from typing import ClassVar, Iterator

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
from .boundary_realizer import SNBoundaryRealizer
from .method_space import SNMethodSpace
from .quadrature import AngularQuadrature
from .spatial.cell_update import CellUpdate, CellVisit
from .spatial.diamond import DiamondDifference
from .spatial.pole_angular_closure import (
    BaileyFlatFluxRedist,
    LegacyTauSymmetricInterpolation,
    MorelMontryAngularSweep,
    PoleAngularClosure,
)
from .sweep_graph import OctantLabel, SweepDependencyGraph


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

    Attributes
    ----------
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
        cell_update: CellUpdate | None = None,
        pole_angular_closure: PoleAngularClosure | None = None,
    ) -> None:
        self.mesh = mesh
        self.quad = quadrature
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
        # Phase B ships three strategies (mirror of Phase A's two
        # boundary-face-flux strategies):
        #
        # * :class:`MorelMontryAngularSweep` (Phase D default) —
        #   canonical Hébert §3.9.4 form with Carlson coupled-pole
        #   seed.  Closes ERR-026 on sphere; preserves cylindrical
        #   Gate 1.1 regression-stability.
        # * :class:`LegacyTauSymmetricInterpolation` — pre-Phase-B
        #   inlined τ-symmetric form.  Bit-identical regression
        #   preservation against the curvilinear snapshots generated
        #   under Phase A.  Carries Defect 3 by design — the
        #   factor-of-two angular truncation gap on angularly-varying
        #   :math:`\\psi` survives.  Reachable via explicit user opt-in.
        # * :class:`BaileyFlatFluxRedist` — the algebraic flat-flux
        #   collapse equivalent (only on flat ψ).  Used by the L1
        #   flat-flux-identity test for ablation studies.
        #
        # Used by the spherical / cylindrical
        # ``transport_operator_matvec_*`` paths via
        # :meth:`SNStreamingOperator.apply`. Cartesian path is unaffected
        # — there is no angular redistribution on slab — and ignores
        # this attribute.
        self.pole_angular_closure: PoleAngularClosure = (
            pole_angular_closure if pole_angular_closure is not None
            else MorelMontryAngularSweep()
        )

        # Normalise to (nx, ny) shaped arrays for both 1-D and 2-D
        if isinstance(mesh, Mesh1D):
            self.nx: int = mesh.N
            self.ny: int = 1
            self.dx: np.ndarray = mesh.widths
            self.dy: np.ndarray = np.array([1.0])
            self.mat_map: np.ndarray = mesh.mat_ids.reshape(mesh.N, 1)
            self._volumes: np.ndarray = mesh.volumes.reshape(mesh.N, 1)
        else:
            self.nx = mesh.nx
            self.ny = mesh.ny
            self.dx = mesh.dx
            self.dy = mesh.dy
            self.mat_map = mesh.mat_map
            self._volumes = mesh.volumes

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
                # (see ``iter_cell_visits``). MoC will define its own
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

    # ── Properties ────────────────────────────────────────────────────

    @property
    def volumes(self) -> np.ndarray:
        """Cell volumes, shape (nx, ny)."""
        return self._volumes

    @property
    def is_1d(self) -> bool:
        """True if this is a 1-D mesh (ny == 1)."""
        return self.ny == 1

    # ── Sweep DAG traversal ───────────────────────────────────────────

    _DEGENERATE_ABS_ETA_THRESHOLD: ClassVar[float] = 1e-15

    def dag_walk(
        self,
        direction_sign: int,
        mu_level_idx: int | None = None,
    ) -> Iterator[CellVisit]:
        r"""Walk the per-ordinate cell DAG in topological order.

        Issue #196 Phase G Step 2.5: the unified iteration primitive
        for 1-D sweeps.  Subsumes :meth:`iter_cells_by_direction` and
        :meth:`iter_cell_visits` — both stay as thin aliases for the
        migration window.

        The SN sweep on a given ordinate is forward substitution on
        the block-triangular streaming + collision operator under the
        ordinate's DAG ordering.  This method yields the per-cell
        visit packets in that DAG order; the consumer folds over the
        packets, threading the spatial-upstream face flux through the
        accumulator and writing the per-cell angular state into a
        persistent array.

        Parameters
        ----------
        direction_sign : int
            ``+1`` for outward (:math:`\mu \ge 0`); ``-1`` for inward
            (:math:`\mu < 0`).
        mu_level_idx : int | None
            For cylindrical geometry: which :math:`\mu`-level the
            ordinate subset belongs to.  ``None`` for slab and
            sphere; required for cylindrical.

        Yields
        ------
        CellVisit
            One per cell, in topological order.  The packet's
            :attr:`face_area_downstream` is float (Issue #196 Step
            2.5 retired ``None``): ``1.0`` for slab, ``0.0`` for
            cylindrical pure-azimuthal degenerate, physical face
            area for sphere / non-degenerate cylinder.
        """
        return self.iter_cells_by_direction(
            direction_sign, mu_level_idx=mu_level_idx,
        )

    def iter_cell_visits(
        self,
        ordinate_idx: int,
        mu_level_idx: int | None = None,
    ) -> Iterator[CellVisit]:
        r"""Yield cells in DAG-topological order for one ordinate.

        The SN sweep is a topological sort of the **directed cell
        graph** where edges are oriented by
        :math:`\mathrm{sign}(\Omega \cdot \hat n_{\text{face}})`.
        This method encapsulates the sweep-direction resolution:
        inward vs outward (1D curvilinear), per-level traversal
        (cylindrical), and the cylindrical pure-azimuthal
        degenerate case.

        For each visit, the :class:`CellVisit` packet contains:

        * The cell index (so the orchestrator knows which cell it
          is working on).
        * The pure-geometric :class:`StreamingTerms` (from
          :meth:`ReducedStreamingOperator.streaming_terms`).
        * The sweep-resolved :attr:`face_area_downstream` (which
          face is outgoing for this visit) — ``None`` for slab and
          for the cylindrical pure-azimuthal degenerate case.

        SN-specific by design.  MoC will not consume this method —
        its mathematical structure is fiber bundles + solution
        sheaves, a different DAG shape.  Premature abstraction
        avoided per Cardinal Rule 2.

        Parameters
        ----------
        ordinate_idx : int
            For slab and sphere: global ordinate index.  Sign of
            ``mu_x[ordinate_idx]`` determines sweep direction
            (outward for :math:`\mu \ge 0`, inward for
            :math:`\mu < 0`).

            For cylindrical: within-level azimuthal index
            :math:`m \in [0, M)`.  The signed :math:`\eta` and the
            global ordinate are resolved through
            ``quad.level_indices[mu_level_idx][ordinate_idx]``.

        mu_level_idx : int | None
            For cylindrical geometry: which :math:`\mu`-level the
            ordinate belongs to.  ``None`` for slab and sphere; a
            ``ValueError`` is raised if missing for cylindrical.

        Yields
        ------
        CellVisit
            One per cell, in topological order for this ordinate.
            For 1-D Cartesian (slab) the order is forward
            (cell 0 → nx-1) for :math:`\mu \ge 0` and backward for
            :math:`\mu < 0`.

        Raises
        ------
        ValueError
            If called on a 2-D Cartesian mesh (no
            :class:`ReducedStreamingOperator`), or if a cylindrical
            mesh is queried without ``mu_level_idx``.

        Notes
        -----
        2-D Cartesian wavefront scheduling is intentionally not
        encapsulated here — its anti-diagonal vectorisation
        operates on cell slices, not per-cell visits.  Wave
        C-extension's LD / EC / Step rollout will revisit this
        when the per-cell ``CellUpdate`` Protocol grows a
        slice-vectorised companion.
        """
        if self.reduced is None:
            raise ValueError(
                "iter_cell_visits is only defined for meshes with a "
                "ReducedStreamingOperator (1-D Cartesian, spherical, "
                "or cylindrical).  2-D Cartesian wavefront sweeps "
                "use anti-diagonal scheduling, not per-cell visits."
            )

        coord = self.reduced.coord

        if coord is CoordSystem.CARTESIAN:
            yield from self._iter_cartesian_visits(ordinate_idx)
            return

        if coord is CoordSystem.SPHERICAL:
            yield from self._iter_spherical_visits(ordinate_idx)
            return

        if coord is CoordSystem.CYLINDRICAL:
            if mu_level_idx is None:
                raise ValueError(
                    "cylindrical iter_cell_visits requires "
                    "mu_level_idx."
                )
            yield from self._iter_cylindrical_visits(
                ordinate_idx, mu_level_idx,
            )
            return

        raise ValueError(  # pragma: no cover — exhaustive match above
            f"Unknown coord system: {coord!r}"
        )

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

    def iter_cells_by_direction(
        self,
        direction_sign: int,
        mu_level_idx: int | None = None,
    ) -> Iterator[CellVisit]:
        r"""Yield cells in DAG-topological order for a sweep direction.

        Issue #168 Phase C — companion to :meth:`iter_cell_visits`
        that surfaces what the cell-visit graph already encodes via
        direction sign **without committing to a specific ordinate**.
        The sweep-frame matvec rewrite (per :doc:`/theory/discrete_ordinates`
        ``phase-c-apply-sweep-equivalence`` subsection) operates on
        whole ordinate subsets simultaneously using
        ``outgoing_mask = quad.mu_x > +eps`` and
        ``incoming_mask = quad.mu_x < -eps``; the cell ordering only
        depends on the direction sign, not on the specific ordinate
        within a sign class.

        For Cartesian / spherical 1D the cell ordering is fully
        determined by ``direction_sign`` — every ordinate with
        :math:`\mu \ge 0` visits cells ``0 → nx-1`` and every ordinate
        with :math:`\mu < 0` visits cells ``nx-1 → 0``.

        For cylindrical, the within-level azimuthal-DAG topology means
        ``mu_level_idx`` MUST be supplied — different :math:`\mu`-levels
        have different per-level cell-graph traversal patterns at the
        pure-azimuthal degenerate :math:`|\eta_n| < 10^{-15}` boundary.

        Foundation test (``tests/sn/test_iter_cells_by_direction.py``)
        pins bit-equivalence to ``iter_cell_visits(ordinate_idx=n)``
        for every representative ordinate ``n`` of the chosen direction.

        Parameters
        ----------
        direction_sign : int
            ``+1`` for outward (:math:`\mu \ge 0`) sweep direction;
            ``-1`` for inward (:math:`\mu < 0`).
        mu_level_idx : int | None
            For cylindrical geometry: which :math:`\mu`-level the
            ordinate subset belongs to.  ``None`` for slab and sphere;
            a :class:`ValueError` is raised if missing for cylindrical.

        Yields
        ------
        CellVisit
            One per cell, in topological order for the direction. The
            packet's :attr:`cell_idx` is the only field the
            sweep-frame matvec consumes; ``streaming_terms`` and
            ``face_area_downstream`` are populated from a
            representative ordinate for backward-compatibility with
            other consumers but the sweep-frame matvec recomputes
            streaming from intrinsic cell geometry directly.

        Raises
        ------
        ValueError
            If ``direction_sign not in (+1, -1)``, or if called on a
            2-D Cartesian mesh (no :class:`ReducedStreamingOperator`),
            or if a cylindrical mesh is queried without
            ``mu_level_idx``, or if no representative ordinate exists
            for ``direction_sign`` in this quadrature.
        """
        if direction_sign not in (+1, -1):
            raise ValueError(
                f"direction_sign must be +1 or -1; got {direction_sign}"
            )
        if self.reduced is None:
            raise ValueError(
                "iter_cells_by_direction is only defined for meshes "
                "with a ReducedStreamingOperator (1-D Cartesian, "
                "spherical, or cylindrical)."
            )
        coord = self.reduced.coord
        if coord is CoordSystem.CYLINDRICAL and mu_level_idx is None:
            raise ValueError(
                "cylindrical iter_cells_by_direction requires "
                "mu_level_idx."
            )

        # Pick a representative ordinate matching the direction sign.
        # The cell-visit graph's cell ordering depends ONLY on the
        # direction sign (and the level for cylindrical), so any
        # ordinate in the correct sign class yields the same cell
        # sequence. We delegate to ``iter_cell_visits`` to keep the
        # implementation a single source of truth.
        # Pick a representative ordinate strictly in the requested
        # sign class — exclude degenerate |η| < 10⁻¹⁵ ordinates because
        # the cylindrical pure-azimuthal degenerate path iterates
        # forward independent of sign (see :meth:`_iter_cylindrical_visits`)
        # and would yield a sequence that disagrees with the bulk
        # direction's signed iteration.
        eps = self._DEGENERATE_ABS_ETA_THRESHOLD

        if coord is CoordSystem.CYLINDRICAL:
            level_indices = self.quad.level_indices  # type: ignore[attr-defined]
            level_ords = np.asarray(level_indices[mu_level_idx])
            eta_at_level = self.quad.mu_x[level_ords]
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
            # ``ordinate_idx`` for cylindrical ``iter_cell_visits`` is
            # the within-level azimuthal index, not the global index.
            yield from self.iter_cell_visits(
                ordinate_idx=int(cand[0]),
                mu_level_idx=mu_level_idx,
            )
            return

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
        yield from self.iter_cell_visits(ordinate_idx=int(cand[0]))

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
        eta_n = float(self.quad.mu_x[global_n])
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
