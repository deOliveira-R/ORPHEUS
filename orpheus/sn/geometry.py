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
from typing import Any, ClassVar, Iterator

import numpy as np

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.geometry.boundary import (
    BoundaryOperator,
    SpecularBoundaryOperator,
    VacuumBoundaryOperator,
)
from orpheus.geometry.reduced_operator import (
    ReducedStreamingOperator,
    cylindrical_streaming,
    slab_streaming,
    spherical_streaming,
)
from .quadrature import AngularQuadrature
from .spatial.boundary_face_flux import BoundaryFaceFlux, DDExtrapolation
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
# SN boundary condition factories
# ═══════════════════════════════════════════════════════════════════════
#
# Factories take a solver-agnostic ``BC(kind, params)`` declaration and
# return a concrete :class:`BoundaryOperator` instance carrying the tensor
# decomposition :math:`R = \\sum_\\alpha G_\\alpha \\otimes A_\\alpha`
# the sweep needs (see :mod:`orpheus.geometry.boundary`).
#
# The 1-D faces (``left`` / ``xmin``, ``right`` / ``xmax``) reflect
# across the radial / x-axis, so the SN factories always pin
# ``axis="x"`` for ``SpecularBoundaryOperator``. 2-D faces dispatch by the y-axis
# variants on ``ymin`` / ``ymax``.

def _sn_vacuum_boundary_operator(sn_mesh: "SNMesh", bc: BC, face: str) -> BoundaryOperator:
    """Zero incoming angular flux at this face."""
    return VacuumBoundaryOperator()


def _sn_reflective_boundary_operator(sn_mesh: "SNMesh", bc: BC, face: str) -> BoundaryOperator:
    """Specular reflection: ψ_in(Ω) = ψ_out(Ω_reflected)."""
    axis = "y" if face in ("ymin", "ymax") else "x"
    return SpecularBoundaryOperator(axis=axis, albedo=1.0)


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
    BOUNDARY_OPERATOR_REGISTRY : dict[str, Callable]
        Supported boundary condition kinds. Each value is a factory
        ``(sn_mesh, bc, face) -> resolved_kind``.  Docstrings on
        the factories serve as descriptions for programmatic query::

            >>> {k: v.__doc__ for k, v in SNMesh.BOUNDARY_OPERATOR_REGISTRY.items()}
    bc_left, bc_right : str
        Resolved BC kind for the left/right (1-D) boundaries.
    bc_xmin, bc_xmax, bc_ymin, bc_ymax : str
        Resolved BC kinds for the four faces of a 2-D mesh.
    """

    BOUNDARY_OPERATOR_REGISTRY: ClassVar[dict[str, Any]] = {
        "vacuum": _sn_vacuum_boundary_operator,
        "reflective": _sn_reflective_boundary_operator,
    }
    # The factories return :class:`BoundaryOperator` instances carrying the
    # tensor decomposition :math:`R = \sum_\alpha G_\alpha \otimes
    # A_\alpha` (see :mod:`orpheus.geometry.boundary`); the sweep calls
    # ``resolved_bc.apply(...)`` directly, so the dispatch
    # is by object identity / Protocol method, not by string tag.
    # Wave C/D will extend this registry with ``"white"`` / ``"periodic"``
    # / ``"albedo"`` entries — the primitives already exist in
    # :mod:`orpheus.geometry.boundary`; only the sweep plumbing remains.

    def __init__(
        self,
        mesh: Mesh1D | Mesh2D,
        quadrature: AngularQuadrature,
        cell_update: CellUpdate | None = None,
        boundary_face_flux: BoundaryFaceFlux | None = None,
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
        # Boundary face-flux strategy (Issue #168 Phase A). Defaults to
        # :class:`DDExtrapolation`, the one-sided second-order DD diamond
        # extrapolation that addresses Defect 1 of Issue #168 (the
        # cell-centre-as-outer-face-flux truncation in the curvilinear
        # FD operator). Used by the spherical / cylindrical
        # ``transport_operator_matvec_*`` paths via
        # :meth:`SNStreamingOperator.apply`. Cartesian path is unaffected
        # — the upwind FD stencil there has no symmetric closure to
        # break — and ignores this attribute.
        self.boundary_face_flux: BoundaryFaceFlux = (
            boundary_face_flux if boundary_face_flux is not None
            else DDExtrapolation()
        )
        # Angular-redistribution closure (Issue #168 Phase B). Defaults
        # to :class:`LegacyTauSymmetricInterpolation`, the bit-for-bit
        # reproduction of the pre-Phase-B inlined τ-symmetric form —
        # preserves the curvilinear regression-snapshot bit-identity
        # contract and the per-ordinate flat-flux invariant that the
        # ERR-026 evidence test
        # (``tests/sn/test_sweep_operator_inconsistency.py``) and the
        # Phase A flat-flux invariants depend on.
        #
        # Phase B ships three strategies (mirror of Phase A's two
        # boundary-face-flux strategies):
        #
        # * :class:`LegacyTauSymmetricInterpolation` (default) —
        #   pre-Phase-B inlined math, bit-identical regression
        #   preservation.  Carries Defect 3 by design — the
        #   factor-of-two angular truncation gap on angularly-varying
        #   :math:`\\psi` survives so future verification probes can
        #   cross-check against the documented behaviour.
        # * :class:`BaileyFlatFluxRedist` — the algebraic flat-flux
        #   collapse equivalent (only on flat ψ).  Used by the L1
        #   flat-flux-identity test.
        # * :class:`MorelMontryAngularSweep` — canonical Hébert §3.9.4
        #   per-cell M-M weighted DD angular recurrence.  Closes
        #   Defect 3 on angularly-varying :math:`\\psi` but breaks
        #   per-ordinate flat-flux balance and does NOT yet pair with
        #   the apply matvec's spatial closure to give a clean
        #   :math:`\\mathcal{O}(h^2)` MMS rate.  The full ERR-026
        #   closure requires a follow-up that aligns the apply
        #   matvec's spatial closure with the sweep's WDD form
        #   (design memo §6.4 / §11).
        #
        # Used by the spherical / cylindrical
        # ``transport_operator_matvec_*`` paths via
        # :meth:`SNStreamingOperator.apply`. Cartesian path is unaffected
        # — there is no angular redistribution on slab — and ignores
        # this attribute.
        self.pole_angular_closure: PoleAngularClosure = (
            pole_angular_closure if pole_angular_closure is not None
            else LegacyTauSymmetricInterpolation()
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
        """Resolve geometry-declared BCs into :class:`BoundaryOperator` instances.

        ``None`` on the mesh defaults to ``BC("reflective")`` (infinite
        lattice / eigenvalue convention). Each face attribute carries
        the concrete :class:`BoundaryOperator` whose ``apply``
        method the sweep invokes directly — no string-kind dispatch
        downstream.
        """
        default = BC("reflective")

        if isinstance(mesh, Mesh1D):
            self.bc_left: BoundaryOperator = self._resolve_one(
                mesh.bc_left or default, "left",
            )
            self.bc_right: BoundaryOperator = self._resolve_one(
                mesh.bc_right or default, "right",
            )
            # Expose 2-D-style attributes for uniform sweep access.
            # The ``y`` faces of a 1-D mesh are degenerate; we tag them
            # with a default reflective ``SpecularBoundaryOperator`` so 2-D-style
            # consumers still get a Protocol-conformant object.
            self.bc_xmin: BoundaryOperator = self.bc_left
            self.bc_xmax: BoundaryOperator = self.bc_right
            self.bc_ymin: BoundaryOperator = SpecularBoundaryOperator(axis="y", albedo=1.0)
            self.bc_ymax: BoundaryOperator = SpecularBoundaryOperator(axis="y", albedo=1.0)
        else:
            self.bc_xmin = self._resolve_one(
                mesh.bc_xmin or default, "xmin",
            )
            self.bc_xmax = self._resolve_one(
                mesh.bc_xmax or default, "xmax",
            )
            self.bc_ymin = self._resolve_one(
                mesh.bc_ymin or default, "ymin",
            )
            self.bc_ymax = self._resolve_one(
                mesh.bc_ymax or default, "ymax",
            )
            self.bc_left = self.bc_xmin
            self.bc_right = self.bc_xmax

    def _resolve_one(self, bc: BC, face: str) -> BoundaryOperator:
        """Look up a single BC in the registry; raise on unsupported kind."""
        factory = self.BOUNDARY_OPERATOR_REGISTRY.get(bc.kind)
        if factory is None:
            supported = ", ".join(f"'{k}'" for k in sorted(self.BOUNDARY_OPERATOR_REGISTRY))
            raise ValueError(
                f"SN solver does not support boundary condition '{bc.kind}' "
                f"on face '{face}'. Supported: {supported}."
            )
        return factory(self, bc, face)

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
        backward for :math:`\\mu < 0`.  Slab DD does not read face
        areas, so ``face_area_downstream`` is ``None``.
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
                face_area_downstream=None,
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
        eta_n = float(self.quad.mu_x[global_n])
        abs_eta = abs(eta_n)

        if abs_eta < self._DEGENERATE_ABS_ETA_THRESHOLD:
            # Pure-azimuthal degenerate: no spatial flow.  Iterate
            # forward so the angular M-M closure runs in a natural
            # order; ``face_area_downstream = None`` signals "no
            # spatial flow" to the strategy.
            for i in range(self.nx):
                st = self.reduced.streaming_terms(
                    cell_idx=i,
                    direction_idx=ordinate_idx,
                    mu_level_idx=mu_level_idx,
                )
                yield CellVisit(
                    cell_idx=i,
                    streaming_terms=st,
                    face_area_downstream=None,
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
