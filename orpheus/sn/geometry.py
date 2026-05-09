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
from typing import Any, ClassVar

import numpy as np

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.geometry.boundary import (
    ResolvedBC,
    SpecularBC,
    VacuumBC,
)
from orpheus.geometry.reduced_operator import (
    ReducedStreamingOperator,
    cylindrical_streaming,
    slab_streaming,
    spherical_streaming,
)
from .quadrature import AngularQuadrature
from .spatial.cell_update import CellUpdate
from .spatial.diamond import DiamondDifference


# ═══════════════════════════════════════════════════════════════════════
# SN boundary condition factories
# ═══════════════════════════════════════════════════════════════════════
#
# Factories take a solver-agnostic ``BC(kind, params)`` declaration and
# return a concrete :class:`ResolvedBC` instance carrying the tensor
# decomposition :math:`R = \\sum_\\alpha G_\\alpha \\otimes A_\\alpha`
# the sweep needs (see :mod:`orpheus.geometry.boundary`).
#
# The 1-D faces (``left`` / ``xmin``, ``right`` / ``xmax``) reflect
# across the radial / x-axis, so the SN factories always pin
# ``axis="x"`` for ``SpecularBC``. 2-D faces dispatch by the y-axis
# variants on ``ymin`` / ``ymax``.

def _sn_bc_vacuum(sn_mesh: "SNMesh", bc: BC, face: str) -> ResolvedBC:
    """Zero incoming angular flux at this face."""
    return VacuumBC()


def _sn_bc_reflective(sn_mesh: "SNMesh", bc: BC, face: str) -> ResolvedBC:
    """Specular reflection: ψ_in(Ω) = ψ_out(Ω_reflected)."""
    axis = "y" if face in ("ymin", "ymax") else "x"
    return SpecularBC(axis=axis, albedo=1.0)


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
    BC_REGISTRY : dict[str, Callable]
        Supported boundary condition kinds. Each value is a factory
        ``(sn_mesh, bc, face) -> resolved_kind``.  Docstrings on
        the factories serve as descriptions for programmatic query::

            >>> {k: v.__doc__ for k, v in SNMesh.BC_REGISTRY.items()}
    bc_left, bc_right : str
        Resolved BC kind for the left/right (1-D) boundaries.
    bc_xmin, bc_xmax, bc_ymin, bc_ymax : str
        Resolved BC kinds for the four faces of a 2-D mesh.
    """

    BC_REGISTRY: ClassVar[dict[str, Any]] = {
        "vacuum": _sn_bc_vacuum,
        "reflective": _sn_bc_reflective,
    }
    # The factories return :class:`ResolvedBC` instances carrying the
    # tensor decomposition :math:`R = \sum_\alpha G_\alpha \otimes
    # A_\alpha` (see :mod:`orpheus.geometry.boundary`); the sweep calls
    # ``resolved_bc.apply_to_incoming(...)`` directly, so the dispatch
    # is by object identity / Protocol method, not by string tag.
    # Wave C/D will extend this registry with ``"white"`` / ``"periodic"``
    # / ``"albedo"`` entries — the primitives already exist in
    # :mod:`orpheus.geometry.boundary`; only the sweep plumbing remains.

    def __init__(
        self,
        mesh: Mesh1D | Mesh2D,
        quadrature: AngularQuadrature,
        cell_update: CellUpdate | None = None,
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
            case CoordSystem.SPHERICAL:
                assert isinstance(mesh, Mesh1D)
                self.reduced = spherical_streaming(mesh, quadrature)
                self.curvature = "spherical"
                self.streaming_x = None
                self.streaming_y = None

        # Resolve boundary conditions from mesh declarations
        self._resolve_bcs(mesh)

    # ── Boundary condition resolution ─────────────────────────────────

    def _resolve_bcs(self, mesh: Mesh1D | Mesh2D) -> None:
        """Resolve geometry-declared BCs into :class:`ResolvedBC` instances.

        ``None`` on the mesh defaults to ``BC("reflective")`` (infinite
        lattice / eigenvalue convention). Each face attribute carries
        the concrete :class:`ResolvedBC` whose ``apply_to_incoming``
        method the sweep invokes directly — no string-kind dispatch
        downstream.
        """
        default = BC("reflective")

        if isinstance(mesh, Mesh1D):
            self.bc_left: ResolvedBC = self._resolve_one(
                mesh.bc_left or default, "left",
            )
            self.bc_right: ResolvedBC = self._resolve_one(
                mesh.bc_right or default, "right",
            )
            # Expose 2-D-style attributes for uniform sweep access.
            # The ``y`` faces of a 1-D mesh are degenerate; we tag them
            # with a default reflective ``SpecularBC`` so 2-D-style
            # consumers still get a Protocol-conformant object.
            self.bc_xmin: ResolvedBC = self.bc_left
            self.bc_xmax: ResolvedBC = self.bc_right
            self.bc_ymin: ResolvedBC = SpecularBC(axis="y", albedo=1.0)
            self.bc_ymax: ResolvedBC = SpecularBC(axis="y", albedo=1.0)
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

    def _resolve_one(self, bc: BC, face: str) -> ResolvedBC:
        """Look up a single BC in the registry; raise on unsupported kind."""
        factory = self.BC_REGISTRY.get(bc.kind)
        if factory is None:
            supported = ", ".join(f"'{k}'" for k in sorted(self.BC_REGISTRY))
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

    @property
    def alpha_half(self) -> np.ndarray:
        """[Deprecated] α dome (sphere). Use ``self.reduced.alpha_half`` instead."""
        warnings.warn(
            "SNMesh.alpha_half is deprecated; "
            "use SNMesh.reduced.alpha_half instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        assert self.reduced.alpha_half is not None
        return self.reduced.alpha_half

    @property
    def redist_dAw(self) -> np.ndarray:
        """[Deprecated] ΔA/w (sphere). Use ``self.reduced.redist_dAw`` instead."""
        warnings.warn(
            "SNMesh.redist_dAw is deprecated; "
            "use SNMesh.reduced.redist_dAw instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        assert self.reduced.redist_dAw is not None
        return self.reduced.redist_dAw

    @property
    def tau_mm(self) -> np.ndarray:
        """[Deprecated] Morel-Montry weights (sphere). Use ``self.reduced.tau_mm`` instead."""
        warnings.warn(
            "SNMesh.tau_mm is deprecated; "
            "use SNMesh.reduced.tau_mm instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        assert self.reduced.tau_mm is not None
        return self.reduced.tau_mm

    @property
    def alpha_per_level(self) -> list[np.ndarray]:
        """[Deprecated] α per μ-level (cylinder). Use ``self.reduced.alpha_per_level`` instead."""
        warnings.warn(
            "SNMesh.alpha_per_level is deprecated; "
            "use SNMesh.reduced.alpha_per_level instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        assert self.reduced.alpha_per_level is not None
        return self.reduced.alpha_per_level

    @property
    def redist_dAw_per_level(self) -> list[np.ndarray]:
        """[Deprecated] ΔA/w per μ-level (cylinder). Use ``self.reduced.redist_dAw_per_level`` instead."""
        warnings.warn(
            "SNMesh.redist_dAw_per_level is deprecated; "
            "use SNMesh.reduced.redist_dAw_per_level instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        assert self.reduced.redist_dAw_per_level is not None
        return self.reduced.redist_dAw_per_level

    @property
    def tau_mm_per_level(self) -> list[np.ndarray]:
        """[Deprecated] Morel-Montry weights per μ-level (cylinder). Use ``self.reduced.tau_mm_per_level`` instead."""
        warnings.warn(
            "SNMesh.tau_mm_per_level is deprecated; "
            "use SNMesh.reduced.tau_mm_per_level instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        assert self.reduced.tau_mm_per_level is not None
        return self.reduced.tau_mm_per_level
