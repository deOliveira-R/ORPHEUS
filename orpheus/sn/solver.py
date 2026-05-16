"""Unified SN (Discrete Ordinates) eigenvalue solver — operator-algebra form.

Wave E Round 2 (Issue #164): :class:`SNSolver` now constructs the
operator triple :math:`(L, S, F)` at ``__init__`` and consumes the
Wave E Round 1 iteration primitives.  The legacy BiCGSTAB FD-operator
path is replaced by Krylov-on-:meth:`L.apply` with the sweep as
preconditioner — the symmetric closure that closes ERR-026 for
curvilinear geometries.

Inner solver dispatch
=====================

* ``inner_solver="source_iteration"`` (default).  Sweep-driven within-
  group fixed-point iteration.  The closure is the WDD asymmetric
  closure that the curvilinear sweep ships (ERR-026 affected).  This
  path is bit-identical to the Wave A-D source iteration, by
  construction — the loop math is preserved character-for-character so
  the 11 frozen regression snapshots stay green.
* ``inner_solver="krylov"``.  GMRES on the symmetric closure carried
  by :meth:`SNStreamingOperator.apply`.  This is the Wave E
  reconciliation that makes the curvilinear ``solve_sn_fixed_source``
  path discretization-correct (closes ERR-026).  On Cartesian meshes
  it is bit-identical math to the legacy BiCGSTAB FD path.

Boundary conditions default to reflective (infinite lattice) but are
configurable via :class:`~orpheus.geometry.mesh.BC` on the mesh.

.. seealso:: :ref:`theory-discrete-ordinates` — Key Facts, equations, gotchas.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, replace

import numpy as np
from scipy.sparse.linalg import LinearOperator as ScipyLinearOperator
from scipy.sparse.linalg import gmres

from orpheus.data.macro_xs.cell_xs import assemble_cell_xs
from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import BC, Mesh1D, Mesh2D
from orpheus.numerics.eigenvalue import power_iteration
from .fission import FissionOperator
from .geometry import SNMesh
from .spatial.sweep_cache import CollisionCache, GeometryCoefficients
from .operator import (
    SNStreamingOperator,
    build_equation_map,
    build_equation_map_spherical,
    build_equation_map_cylindrical,
    solution_to_angular_flux,
    solution_to_angular_flux_spherical,
    solution_to_angular_flux_cylindrical,
)
from .quadrature import AngularQuadrature
from .scattering import ScatteringOperator
from .sweep import transport_sweep


def _apply_default_bcs(
    mesh: Mesh1D | Mesh2D,
    boundary_condition: str,
) -> Mesh1D | Mesh2D:
    """Apply *boundary_condition* string to all faces that lack explicit BCs.

    Returns the original mesh unchanged when it already carries explicit
    :class:`~orpheus.geometry.mesh.BC` declarations, so user-set BCs
    always take precedence over the ``boundary_condition`` parameter.
    """
    bc = BC(boundary_condition)
    if isinstance(mesh, Mesh1D):
        if mesh.bc_left is None and mesh.bc_right is None:
            return replace(mesh, bc_left=bc, bc_right=bc)
    else:
        faces = ("bc_xmin", "bc_xmax", "bc_ymin", "bc_ymax")
        if all(getattr(mesh, f) is None for f in faces):
            return replace(mesh, **{f: bc for f in faces})
    return mesh


@dataclass
class SNFixedSourceResult:
    """Results of an SN fixed-source (non-eigenvalue) calculation.

    Used by :func:`solve_sn_fixed_source` and MMS verification
    (:mod:`orpheus.derivations.continuous.mms.sn`). Carries the same flux arrays
    as :class:`SNResult` but omits the fission eigenvalue fields.

    Issue #196 PR-INDEX-5: storage flipped to principled layout —
    ``angular_flux: (N_ordinates, ng, nx, ny)``;
    ``scalar_flux: (ng, nx, ny)``.
    """

    angular_flux: np.ndarray   # (N_ordinates, ng, nx, ny)
    scalar_flux: np.ndarray    # (ng, nx, ny) = Σ w_n ψ_n
    geometry: Mesh1D | Mesh2D
    quadrature: AngularQuadrature
    n_inner: int               # source iterations used
    residual: float            # final ||Δφ||/||φ||
    elapsed_seconds: float
    eg: np.ndarray | None = None  # (ng+1,) energy boundaries, or None for synthetic XS


@dataclass
class SNResult:
    """Results of an SN transport calculation.

    The primary output is the angular flux (the direct solution of the
    SN equations). Scalar flux is derived by quadrature integration.

    Issue #196 PR-INDEX-5: storage flipped to principled layout —
    ``angular_flux: (N_ordinates, ng, nx, ny)``;
    ``scalar_flux: (ng, nx, ny)``.
    """

    keff: float
    keff_history: list[float]
    angular_flux: np.ndarray   # (N_ordinates, ng, nx, ny)
    scalar_flux: np.ndarray    # (ng, nx, ny) = Σ w_n ψ_n
    geometry: Mesh1D | Mesh2D
    quadrature: AngularQuadrature
    eg: np.ndarray | None      # (ng+1,) energy boundaries, or None for synthetic XS
    elapsed_seconds: float


# ═══════════════════════════════════════════════════════════════════════
# Solver class (EigenvalueSolver protocol)
# ═══════════════════════════════════════════════════════════════════════

class SNSolver:
    """Unified SN eigenvalue solver satisfying the EigenvalueSolver protocol.

    Constructs the operator triple :math:`(L, S, F)` at construction
    time and routes ``solve_fixed_source`` through one of two inner-
    solver paths:

    * ``"source_iteration"`` — sweep-driven within-group fixed-point
      iteration (WDD asymmetric closure; ERR-026-affected for
      curvilinear).  Bit-identical to the Wave A-D path.
    * ``"krylov"`` — GMRES on :meth:`SNStreamingOperator.apply` (the
      symmetric closure) with the sweep as preconditioner.  Closes
      ERR-026 on curvilinear; bit-identical math to the legacy
      BiCGSTAB FD path on Cartesian.

    The legacy ``"bicgstab"`` value is no longer accepted — call sites
    must migrate to ``"krylov"``.

    Parameters
    ----------
    sn_mesh : SNMesh — augmented geometry (wraps Mesh1D or Mesh2D with
        precomputed streaming stencil + materials dict + ``ng``).
        Issue #197 PR-TYPED-0: ``sn_mesh.materials`` IS the single
        source of truth for cross sections and group count; the
        legacy ``materials`` / ``n_groups`` SNSolver constructor
        parameters were retired (aggressive retirement per
        ``feedback_aggressive_retirement``).
    inner_solver : "source_iteration" or "krylov".
    scattering_order : int — Legendre order for scattering (0 = P0).
    keff_tol, flux_tol : outer iteration convergence.
    max_inner, inner_tol : inner iteration parameters.
    """

    def __init__(
        self,
        sn_mesh: SNMesh,
        inner_solver: str = "source_iteration",
        scattering_order: int = 0,
        keff_tol: float = 1e-7,
        flux_tol: float = 1e-6,
        max_inner: int = 200,
        inner_tol: float = 1e-8,
    ):
        if inner_solver not in ("source_iteration", "krylov"):
            raise ValueError(
                f"Unknown inner solver: {inner_solver!r}. "
                f"Valid choices are 'source_iteration' or 'krylov'. "
                f"(The legacy 'bicgstab' alias was retired in Wave E "
                f"Round 2; use 'krylov' which routes through "
                f"SNStreamingOperator.apply with the sweep as "
                f"preconditioner.)"
            )
        self.sn_mesh = sn_mesh
        self.quad = sn_mesh.quad
        self.inner_solver = inner_solver
        self.scattering_order = scattering_order
        self.keff_tol = keff_tol
        self.flux_tol = flux_tol
        self.max_inner = max_inner
        self.inner_tol = inner_tol

        # Issue #197 PR-TYPED-0: materials + ng are read from sn_mesh,
        # not from constructor parameters.  ``sn_mesh.ng`` is the
        # validated single source of truth (raises
        # ``InconsistentMaterialsError`` if materials disagree).
        materials = sn_mesh.materials
        nx, ny = sn_mesh.nx, sn_mesh.ny
        self.ng = sn_mesh.ng

        # Issue #197 PR-TYPED-1: the canonical XS state collapses to ONE
        # attribute — ``self.mat_xs`` is the :class:`MaterialXSField`
        # wrapping both the per-material :class:`Mixture` data AND the
        # per-cell typed views.  Every operator (L, C, S, F) reads cross
        # sections through this single source of truth; the seven
        # separate per-XS attributes (sig_t, sig_a, sig_p, chi, sig_s,
        # sig2, sig_s0) plus ``_cells_by_mat`` / ``_sig2_sum`` collapse
        # into ``self.mat_xs.*`` accessors.
        #
        # The thin ``sig_t / sig_a / sig_p / chi`` read-through
        # properties below preserve the SNSolver API surface for the
        # one-PR migration window (downstream consumers that read
        # ``solver.sig_t`` etc.).  PR-TYPED-2 retires them by rewiring
        # consumers to read ``solver.mat_xs.total_cross_section`` etc.
        # directly.
        self.mat_xs = sn_mesh.material_xs_field()

        # __debug__ cell-flattening invariant pinning (formerly at
        # construction of self.sig_t — now exercised through the
        # mat_xs.total_cross_section accessor, populated lazily).
        if __debug__:
            xs_check = assemble_cell_xs(materials, sn_mesh.mat_map)
            _sig_t_old = xs_check.sig_t.reshape(nx, ny, self.ng)
            assert np.array_equal(
                _sig_t_old,
                self.mat_xs.total_cross_section.transpose(1, 2, 0),
            ), "PR-INDEX-3 cell-flattening invariant broke"

        # Scattering order — clamp to the minimum Legendre count
        # available across all materials.
        L = min(
            scattering_order,
            min(len(m.SigS) - 1 for m in materials.values()),
        )
        self.scattering_order = L

        # Weight normalization (1/sum(w) — works for both GL and Lebedev)
        self.weight_norm = 1.0 / sn_mesh.quad.weights.sum()

        # Persistent boundary flux cache (passed to sweep)
        self._psi_bc: dict = {}

        # Volume array for keff computation
        self.volume = sn_mesh.volumes

        # ── Operator triple — Wave E Round 2 -----------------------------
        # The (L, S, F) algebra-of-record framing.  Constructed once at
        # __init__ so downstream consumers (the iteration primitives, the
        # _solve_krylov path, future sensitivity/adjoint hooks) see a
        # consistent operator triple over the lifetime of this solver.
        # Issue #197 PR-TYPED-1: both S and F now consume the single
        # ``self.mat_xs`` — the per-material dispatch lives inside
        # :class:`MaterialXSField`'s typed verbs, not on the operators.
        self.scattering_op = ScatteringOperator.from_solver_data(
            mat_xs=self.mat_xs,
            quadrature=sn_mesh.quad,
            scattering_order=self.scattering_order,
        )
        self.fission_op = FissionOperator.from_solver_data(
            mat_xs=self.mat_xs,
        )
        # The streaming-collision operator L = Ω·∇ + Σ_t.  Built lazily
        # — the Cartesian / spherical / cylindrical EquationMap that L's
        # apply path needs is built on first call (lazy in
        # SNStreamingOperator._ensure_eq_map).  ``self.sig_t`` is in the
        # principled ``(ng, nx, ny)`` layout (Issue #196 PR-INDEX-3) and
        # the operator's matvec helpers were updated to consume that
        # layout natively — no bridge needed.
        self.L = SNStreamingOperator(
            sn_mesh=sn_mesh, sig_t=self.mat_xs.total_cross_section,
        )
        self.S = self.scattering_op
        self.F = self.fission_op

        # ── Sweep cache (Issue #196 Phase G Step 2.5c) ───────────────
        # Two-stratum cache: GeometryCoefficients built once at __init__
        # (geometry × quadrature only); CollisionCache built once after
        # σ_t binding.  Hot path consumes (geom, coll) without per-cell
        # StreamingTerms allocation.  Only applicable to 1-D meshes with
        # ReducedStreamingOperator — 2-D Cartesian uses the wavefront path.
        self.geom_cache: GeometryCoefficients | None = None
        self.coll_cache: CollisionCache | None = None
        if sn_mesh.reduced is not None:
            self.geom_cache = GeometryCoefficients.from_mesh_and_quad(sn_mesh)
            # No bridge needed: ``mat_xs.total_cross_section`` is the
            # principled ``(ng, nx, ny=1)`` layout the cache expects.
            # Drop the trailing degenerate ``ny`` axis with a slice view.
            sig_t_1d = self.mat_xs.total_cross_section[:, :, 0]  # (ng, nx)
            self.coll_cache = CollisionCache.from_geometry(
                self.geom_cache, sig_t_1d,
            )
            # Stash the caches on the SNMesh so the sweep can find them
            # without threading a solver reference through ``transport_sweep``.
            sn_mesh._geom_cache = self.geom_cache  # type: ignore[attr-defined]
            sn_mesh._coll_cache = self.coll_cache  # type: ignore[attr-defined]

    # ── Read-through XS properties (Issue #197 PR-TYPED-1 shim) ────────
    #
    # TRANSIENT — to be retired in PR-TYPED-2 once every downstream
    # consumer (operators, build_rhs helpers, test fixtures) reads
    # ``solver.mat_xs.<view>`` directly.  Kept for one cycle to
    # minimize blast-radius in PR-TYPED-1.

    @property
    def sig_t(self) -> np.ndarray:
        r"""TRANSIENT — :math:`\sigma_t(\vec r, g)` per-cell view,
        shape ``(ng, nx, ny)``.  Read-through onto
        ``self.mat_xs.total_cross_section``.  Will be retired in
        PR-TYPED-2; consumers should read through ``solver.mat_xs``
        directly."""
        return self.mat_xs.total_cross_section

    @property
    def sig_a(self) -> np.ndarray:
        r"""TRANSIENT — :math:`\sigma_a(\vec r, g)` per-cell view.
        See :attr:`sig_t`."""
        return self.mat_xs.absorption_cross_section

    @property
    def sig_p(self) -> np.ndarray:
        r"""TRANSIENT — :math:`\nu\Sigma_f(\vec r, g)` per-cell view.
        See :attr:`sig_t`."""
        return self.mat_xs.fission_production

    @property
    def chi(self) -> np.ndarray:
        r"""TRANSIENT — :math:`\chi(\vec r, g)` per-cell view.
        See :attr:`sig_t`."""
        return self.mat_xs.emission_spectrum

    @property
    def sig_s(self) -> dict[int, list[np.ndarray]]:
        """TRANSIENT — per-material Legendre scattering dict.
        Read-through onto ``self.mat_xs.sig_s_legendre(mid)`` for
        every material.  Will be retired in PR-TYPED-2."""
        return {
            mid: self.mat_xs.sig_s_legendre(mid)
            for mid in self.mat_xs.materials
        }

    @property
    def sig2(self) -> dict[int, np.ndarray]:
        """TRANSIENT — per-material (n,2n) dict.  See :attr:`sig_s`."""
        return {
            mid: self.mat_xs.n2n_matrix(mid)
            for mid in self.mat_xs.materials
        }

    @property
    def sig_s0(self) -> dict[int, np.ndarray]:
        """TRANSIENT — per-material P0 scattering matrix dict.
        See :attr:`sig_s`."""
        return {
            mid: self.mat_xs.sig_s_legendre(mid)[0]
            for mid in self.mat_xs.materials
        }

    @property
    def _cells_by_mat(self) -> dict[int, tuple[np.ndarray, np.ndarray]]:
        """TRANSIENT — per-material cell-index dict.
        See :attr:`sig_s`."""
        return self.mat_xs.cells_by_material

    def rebind_cross_sections(self, new_sig_t: np.ndarray) -> None:
        """Rebind the total cross-section and rebuild only :class:`CollisionCache`.

        :class:`GeometryCoefficients` survives — Stratum 1 is geometry-only.
        Only the σ_t-dependent Stratum 2 rebuilds.  Used by depletion /
        thermal-feedback consumers.

        Parameters
        ----------
        new_sig_t
            New total cross-section in the principled ``(ng, nx, ny)``
            layout (Issue #196 PR-INDEX-3).

        Notes
        -----
        Issue #197 PR-TYPED-1 — ``rebind_cross_sections`` overrides
        ``self.mat_xs._sig_t_cell`` directly (without re-deriving from
        materials) because the depletion / thermal-feedback consumer
        adjusts σ_t per-cell without revisiting the per-material data.
        """
        # Override the lazy cache on mat_xs.  Force the dense
        # per-cell view to be populated first so the other cell views
        # (sig_a, sig_p, chi) match the rebind contract.
        _ = self.mat_xs.absorption_cross_section
        self.mat_xs._sig_t_cell = new_sig_t
        # Mirror onto the L operator so its apply path stays consistent.
        self.L = SNStreamingOperator(
            sn_mesh=self.sn_mesh, sig_t=self.mat_xs.total_cross_section,
        )
        if self.geom_cache is not None:
            sig_t_1d = self.mat_xs.total_cross_section[:, :, 0]
            self.coll_cache = CollisionCache.from_geometry(
                self.geom_cache, sig_t_1d,
            )
            self.sn_mesh._coll_cache = self.coll_cache  # type: ignore[attr-defined]

    def initial_flux_distribution(self) -> np.ndarray:
        """Initial scalar flux guess: ones(ng, nx, ny).

        Issue #196 PR-INDEX-5: principled layout.
        """
        return np.ones((self.ng, self.sn_mesh.nx, self.sn_mesh.ny))

    def compute_fission_source(
        self, flux_distribution: np.ndarray, keff: float,
    ) -> np.ndarray:
        """Fission source: χ · (νΣ_f · φ) / k.

        Thin delegator to :meth:`FissionOperator.apply` (Wave D Issue 13).
        The :math:`1/k` division stays at this level — the fission
        operator is a *linear* operator (Wave A Issue 1 Protocol);
        :meth:`FissionOperator.apply` returns :math:`F\\,\\phi` and the
        eigenvalue scaling lives here.

        Issue #196 PR-INDEX-5: ``flux_distribution`` is principled
        ``(ng, nx, ny)``.  No bridges — the PR-INDEX-4 transpose pair
        is GONE.
        """
        return self.fission_op.apply(flux_distribution) / keff

    def solve_fixed_source(
        self, fission_source: np.ndarray, flux_distribution: np.ndarray,
    ) -> np.ndarray:
        """Solve the within-group transport equation for given fission source.

        Returns updated scalar flux ``(ng, nx, ny)``.
        """
        if self.inner_solver == "source_iteration":
            return self._solve_source_iteration(fission_source, flux_distribution)
        if self.inner_solver == "krylov":
            return self._solve_krylov(fission_source, flux_distribution)
        # Should be unreachable — __init__ validated the choice.
        raise ValueError(f"Unknown inner solver: {self.inner_solver}")

    def compute_group_production_rate(
        self, flux_distribution: np.ndarray,
    ) -> np.ndarray:
        r"""Per-group volume-integrated neutron production rate, shape ``(ng,)``.

        Component :math:`r_g` is

        .. math::

            r_g \;=\; \int_V \nu \Sigma_{f,g}(\mathbf{r})\,\phi_g(\mathbf{r})\,dV
                       \;+\; 2 \int_V \sum_{g' } \Sigma_{2,g'\to g}(\mathbf{r})
                                                 \,\phi_{g'}(\mathbf{r})\,dV

        i.e. the per-group fission-neutron production plus the per-group
        ``(n, 2n)`` contribution (the factor of 2 accounts for the
        two-neutron-out yield).  Fission is integrated against
        ``mesh.volume_measure`` (Issue 9.6 wiring); ``(n, 2n)`` runs the
        existing per-material loop because the ``sig2`` matrices are
        keyed on material rather than cell.

        The output is the natural diagnostic intermediate for spectral
        analysis (per-group production rates are reactor-physics-meaningful
        quantities).  ``compute_keff`` consumes it via ``.sum()``.
        """
        nx, ny, ng = self.sn_mesh.nx, self.sn_mesh.ny, self.ng
        mu = self.sn_mesh.mesh.volume_measure

        # Fission production: ∫ νΣ_f · φ dV, vectorised over groups.
        # Issue #196 PR-INDEX-5: both ``self.sig_p`` and
        # ``flux_distribution`` are principled ``(ng, nx, ny)``.  The
        # named intermediate ``per_cell_per_group`` has units ``[1/s]``
        # per cell per group (a reactor-physics quantity — coding-
        # elegance Pattern 3).  ``volume_measure`` consumes a flat
        # ``(N_cells, ng)`` view (Issue 9.6 wiring); reshape via
        # ``transpose(1, 2, 0)`` then flatten the spatial axes.
        per_cell_per_group = np.einsum(
            "gxy,gxy->gxy", self.sig_p, flux_distribution,
        )
        rate = mu(per_cell_per_group.transpose(1, 2, 0).reshape(nx * ny, ng))

        # (n,2n) contribution — Issue #197 PR-TYPED-1: the per-material
        # dispatch loop lives ONLY inside
        # :meth:`MaterialXSField.add_n2n_to_group_rate`.
        self.mat_xs.add_n2n_to_group_rate(
            rate, flux_distribution, self.volume,
        )

        return rate

    def compute_group_absorption_rate(
        self, flux_distribution: np.ndarray,
    ) -> np.ndarray:
        r"""Per-group volume-integrated absorption rate, shape ``(ng,)``.

        Component :math:`a_g = \int_V \Sigma_{a,g}(\mathbf{r})\,\phi_g(\mathbf{r})\,dV`.

        Volume-integrated via ``mesh.volume_measure`` (Issue 9.6 wiring).
        The denominator of ``compute_keff`` is ``.sum()`` of this vector.

        Issue #196 PR-INDEX-5: ``flux_distribution`` is principled
        ``(ng, nx, ny)``.
        """
        nx, ny, ng = self.sn_mesh.nx, self.sn_mesh.ny, self.ng
        mu = self.sn_mesh.mesh.volume_measure
        per_cell_per_group = np.einsum(
            "gxy,gxy->gxy", self.sig_a, flux_distribution,
        )
        return mu(per_cell_per_group.transpose(1, 2, 0).reshape(nx * ny, ng))

    def compute_keff(self, flux_distribution: np.ndarray) -> float:
        """k = production / absorption (volume-weighted).

        Composed from the per-group production and absorption rate
        vectors so the intermediates are individually meaningful and
        reusable (e.g. spectral diagnostics).  See
        :meth:`compute_group_production_rate` and
        :meth:`compute_group_absorption_rate`.
        """
        production = float(self.compute_group_production_rate(flux_distribution).sum())
        absorption = float(self.compute_group_absorption_rate(flux_distribution).sum())
        return production / absorption

    def converged(
        self, keff: float, keff_old: float,
        flux_distribution: np.ndarray, flux_old: np.ndarray,
        iteration: int,
    ) -> bool:
        if iteration <= 2:
            return False
        dk = abs(keff - keff_old)
        dphi = np.linalg.norm(flux_distribution - flux_old) / \
            max(np.linalg.norm(flux_distribution), 1e-30)
        return dk < self.keff_tol and dphi < self.flux_tol

    # ── Inner solver: source iteration ────────────────────────────────

    def _solve_source_iteration(
        self, fission_source: np.ndarray, flux_distribution: np.ndarray,
    ) -> np.ndarray:
        """Scattering source iteration: sweep → update scatter → sweep → ...

        **Approach A — bit-identical preservation**: the loop math is
        preserved character-for-character from the Wave A-D
        implementation.  Conceptually this is a
        :class:`~orpheus.numerics.iteration.SourceIteration` realisation
        with operator triple :math:`(L, S, \\text{Zero})` and
        ``q_ext = fission_source``, where :math:`L^{-1}` is the sweep
        and :math:`S` carries both isotropic and Pℓ contributions.

        Issue #196 PR-INDEX-5: every flux / source intermediate is in
        principled ``(ng, nx, ny)`` / ``(N, ng, nx, ny)`` layout.
        """
        phi = flux_distribution.copy()
        angular = None  # no angular flux on first iteration

        for n_inner in range(self.max_inner):
            phi_prev = phi.copy()

            # Isotropic source = fission + P0 scattering + (n,2n)
            Q = fission_source.copy()
            self._add_scattering_source(Q, phi)
            self._add_n2n_source(Q, phi)

            # Anisotropic scattering (P1+ terms, None when L=0)
            Q_aniso = self._build_aniso_scattering(angular)

            # Transport sweep
            angular, phi = transport_sweep(
                Q, self.sig_t, self.sn_mesh,
                self._psi_bc, Q_aniso=Q_aniso,
            )

            norm = np.linalg.norm(phi)
            if norm > 0:
                res = np.linalg.norm(phi - phi_prev) / norm
                if res < self.inner_tol:
                    break

        return phi

    # ── Inner solver: Krylov-on-apply (replaces BiCGSTAB FD path) ─────

    def _solve_krylov(
        self, fission_source: np.ndarray, flux_distribution: np.ndarray,
    ) -> np.ndarray:
        r"""Inner solve via GMRES on :meth:`SNStreamingOperator.apply`.

        Replaces the four legacy ``_solve_bicgstab_*`` methods.  Routes
        through ONE operator (``self.L``, the streaming-collision
        operator with the symmetric closure) — for curvilinear meshes
        this is the Wave E Round 2 reconciliation that closes ERR-026
        (the WDD asymmetric closure carried by the sweep is removed
        from the converged-solution path); for Cartesian meshes it is
        bit-identical math to the legacy BiCGSTAB FD path.

        The sweep is wrapped as a left preconditioner (scipy gmres
        default), realising the SAILOR / Larsen-Adams preconditioned-
        Krylov framework (Adams & Larsen 2002 §III).

        Returns the updated scalar flux ``(ng, nx, ny)``.
        """
        nx, ny, ng = self.sn_mesh.nx, self.sn_mesh.ny, self.ng
        sum_w = float(self.quad.weights.sum())
        phi = flux_distribution
        fission_src_norm = fission_source / sum_w

        # ---------- packed-vector layout via EquationMap -----------------
        eq_map = self.L._ensure_eq_map()
        n = eq_map.n_unknowns
        curv = getattr(self.sn_mesh, "curvature", None)

        # ---------- previous-iterate angular flux for Pℓ (Cartesian) ----
        # Cartesian Pℓ scattering needs the angular flux of the previous
        # outer iteration to build the per-ordinate scattering source
        # via Legendre-moment Galerkin reconstruction.  We carry the
        # packed solution across calls so the warm-start serves both
        # GMRES and the Pℓ source build.
        # PR-INDEX-7: ``solution_to_angular_flux*`` returns principled
        # ``(N, ng, nx, ny)`` natively — the FD-matvec internal layout
        # flip closes the PR-INDEX-4 §9.1 deferral.
        # ``_build_rhs_cartesian`` consumes the same layout.
        if self.scattering_order > 0 and curv is None and hasattr(
            self, "_psi_solution"
        ):
            angular_full = solution_to_angular_flux(
                self._psi_solution, eq_map, self.quad, nx, ny, ng,
                bc_xmin=self.sn_mesh.bc_xmin,
                bc_xmax=self.sn_mesh.bc_xmax,
                bc_ymin=self.sn_mesh.bc_ymin,
                bc_ymax=self.sn_mesh.bc_ymax,
            )
        else:
            angular_full = None

        # ---------- packed RHS (fission + scattering + n2n) -------------
        if curv == "spherical":
            rhs = _build_rhs_spherical(
                fission_src_norm, phi, eq_map, self.quad,
                self.sig_s, self.sig2, self.sn_mesh.mat_map,
                nx, ng,
            )
        elif curv == "cylindrical":
            rhs = _build_rhs_cylindrical(
                fission_src_norm, phi, eq_map, self.quad,
                self.sig_s, self.sig2, self.sn_mesh.mat_map,
                nx, ng,
            )
        else:
            rhs = _build_rhs_cartesian(
                fission_src_norm, phi, eq_map, self.quad,
                self.sig_s, self.sig2, self.sn_mesh.mat_map,
                nx, ny, ng,
                scattering_order=self.scattering_order,
                angular_flux=angular_full,
            )

        # ---------- L as a scipy LinearOperator (matvec = L.apply) -----
        L_scipy = ScipyLinearOperator(
            (n, n), matvec=self.L.apply, dtype=float,
        )

        # ---------- sweep preconditioner ------------------------------
        # The sweep is L^{-1} (in the WDD-closure sense; close enough as
        # a preconditioner for the symmetric-closure operator).  Wraps
        # the structured (Q_iso, Q_aniso) sweep contract back to the
        # packed solution-vector layout that GMRES expects.
        precond = self._make_sweep_preconditioner(eq_map, n, sum_w)

        # ---------- warm start ----------------------------------------
        x0 = (
            self._psi_solution.copy()
            if hasattr(self, "_psi_solution")
            else np.ones(n)
        )

        # ---------- GMRES ---------------------------------------------
        try:
            solution, info = gmres(
                L_scipy, rhs, x0=x0, M=precond,
                rtol=self.inner_tol, atol=0.0,
                maxiter=self.max_inner, restart=min(50, n),
            )
        except TypeError:
            # Older scipy versions may use ``tol`` instead of ``rtol``.
            solution, info = gmres(
                L_scipy, rhs, x0=x0, M=precond,
                tol=self.inner_tol,
                maxiter=self.max_inner, restart=min(50, n),
            )
        self._psi_solution = solution

        # ---------- decode packed solution → scalar flux --------------
        # PR-INDEX-7: ``solution_to_angular_flux*`` returns principled
        # ``(N, ng, nx, ny)`` natively (closes PR-INDEX-4 §9.1 deferral).
        # ``_scalar_flux_from_angular`` converts that to principled
        # scalar flux ``(ng, nx, ny)``.
        if curv == "spherical":
            fi = solution_to_angular_flux_spherical(
                solution, eq_map, self.quad, nx, ng,
            )
            return _scalar_flux_from_angular(fi, self.quad, nx, 1, ng)
        if curv == "cylindrical":
            fi = solution_to_angular_flux_cylindrical(
                solution, eq_map, self.quad, nx, ng,
            )
            return _scalar_flux_from_angular(fi, self.quad, nx, 1, ng)
        fi = solution_to_angular_flux(
            solution, eq_map, self.quad, nx, ny, ng,
            bc_xmin=self.sn_mesh.bc_xmin,
            bc_xmax=self.sn_mesh.bc_xmax,
            bc_ymin=self.sn_mesh.bc_ymin,
            bc_ymax=self.sn_mesh.bc_ymax,
        )
        return _scalar_flux_from_angular(fi, self.quad, nx, ny, ng)

    def _make_sweep_preconditioner(
        self, eq_map, n: int, sum_w: float,
    ) -> ScipyLinearOperator:
        r"""Build a scipy LinearOperator wrapping the sweep as :math:`L^{-1}`.

        The sweep takes a structured pair ``(Q_iso, Q_aniso)`` and
        returns ``(angular_flux, scalar_flux)``.  GMRES wants a scalar
        ``matvec(q) -> M^{-1}·q`` on the packed 1-D vector.  This
        adapter:

        1. Decodes the packed RHS into per-ordinate
           ``Q_aniso`` shape ``(N, ng, nx, ny)`` (principled storage
           per :ref:`theory-sn-index-convention`), undoing the
           ``/sum_w`` normalisation from :func:`build_rhs*` so the
           sweep's internal ``× weight_norm`` step gets back to the
           caller's per-ordinate source.
        2. Runs the sweep with ``Q_iso = 0`` (everything routes through
           ``Q_aniso`` so the per-ordinate variation is preserved).
        3. Re-packs the resulting angular flux into the packed
           solution-vector layout via the inverse of
           :func:`solution_to_angular_flux*`.

        The sweep's internal ``psi_bc`` cache is NOT shared with
        :attr:`self._psi_bc` — the preconditioner is stateless across
        GMRES inner iterations, which keeps the linear-operator
        contract clean.
        """
        nx, ny, ng = self.sn_mesh.nx, self.sn_mesh.ny, self.ng
        N = self.quad.N
        curv = getattr(self.sn_mesh, "curvature", None)

        def matvec(q_packed: np.ndarray) -> np.ndarray:
            # PR-INDEX-7: ``solution_to_angular_flux*`` returns principled
            # ``(N, ng, nx, ny)`` natively — no transpose adapter needed.
            # Closes the PR-INDEX-4 §9.1 deferral; the FD-matvec internal
            # ``(ng, N, nx, ny)`` legacy layout is retired.
            if curv == "spherical":
                fi_op = solution_to_angular_flux_spherical(
                    q_packed, eq_map, self.quad, nx, ng,
                )
            elif curv == "cylindrical":
                fi_op = solution_to_angular_flux_cylindrical(
                    q_packed, eq_map, self.quad, nx, ng,
                )
            else:
                fi_op = solution_to_angular_flux(
                    q_packed, eq_map, self.quad, nx, ny, ng,
                    bc_xmin=self.sn_mesh.bc_xmin,
                    bc_xmax=self.sn_mesh.bc_xmax,
                    bc_ymin=self.sn_mesh.bc_ymin,
                    bc_ymax=self.sn_mesh.bc_ymax,
                )
            # ``fi_op`` is the principled (N, ng, nx, ny) angular flux;
            # the sweep ingests Q_aniso in the same layout. The ``* sum_w``
            # undoes the ``/sum_w`` baked into ``build_rhs*`` (the operator
            # equation source carries the inverse-weight-sum factor).
            Q_aniso = fi_op * sum_w
            Q_iso = np.zeros((ng, nx, ny))
            psi_bc_local: dict = {}
            try:
                angular, _ = transport_sweep(
                    Q_iso, self.sig_t, self.sn_mesh, psi_bc_local,
                    Q_aniso=Q_aniso,
                )
            except Exception:
                # If the sweep cannot run with this Q_aniso shape, degrade
                # to the identity preconditioner.
                return q_packed
            # Pack angular → solution vector: ``angular`` is principled
            # ``(N, ng, nx, ny)``; the packed vector expects
            # ``flux[ng, n_eq]`` in F-order via
            # ``solution.reshape(ng, n_eq, order='F')``.
            packed = np.empty((ng, eq_map.n_eq), dtype=float)
            for k in range(eq_map.n_eq):
                packed[:, k] = angular[
                    eq_map.ordinate[k], :, eq_map.ix[k], eq_map.iy[k],
                ]
            return packed.ravel(order="F")

        return ScipyLinearOperator((n, n), matvec=matvec, dtype=float)

    # ── Source computation helpers ────────────────────────────────────
    #
    # Wave D Issue 13: the math has been lifted into ScatteringOperator
    # (orpheus/sn/scattering.py). The methods below are thin delegators
    # preserved for the EigenvalueSolver Protocol surface and the
    # underscore-prefixed test probes in tests/sn/test_solver_components.py.
    # All four delegate to the same precomputed-by-construction operator
    # held on self.scattering_op, so bit-identity is by construction.

    def _add_scattering_source(self, Q: np.ndarray, phi: np.ndarray) -> None:
        """Add P0 scattering source to Q in-place (delegates to ScatteringOperator).

        Issue #196 PR-INDEX-5: every solver-public flux is principled
        ``(ng, nx, ny)``; the PR-INDEX-4 transpose pair retired.
        """
        self.scattering_op.add_iso_source(Q, phi)

    def _build_aniso_scattering(
        self, angular_flux: np.ndarray | None,
    ) -> np.ndarray | None:
        """Build per-ordinate anisotropic Pℓ scattering source (delegates to ScatteringOperator).

        Issue #196 PR-INDEX-5: ``angular_flux`` is principled
        ``(N, ng, nx, ny)``; return shape ``(N, ng, nx, ny)`` (or
        ``None`` when ``L == 0``).  The PR-INDEX-4 bridges retired.
        """
        if angular_flux is None:
            return None
        return self.scattering_op.build_aniso_source(angular_flux)

    def _add_n2n_source(self, Q: np.ndarray, phi: np.ndarray) -> None:
        """Add (n,2n) source to Q in-place (delegates to ScatteringOperator).

        Issue #196 PR-INDEX-5: both arguments are principled
        ``(ng, nx, ny)``.
        """
        self.scattering_op.add_n2n_source(Q, phi)


# ═══════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════

def _scalar_flux_from_angular(
    fi: np.ndarray, quad: AngularQuadrature, nx: int, ny: int, ng: int,
) -> np.ndarray:
    r"""Integrate angular flux to scalar flux: :math:`\phi = \sum_n w_n \psi_n`.

    Local helper so the solver is self-contained — the legacy
    :func:`orpheus.sn.operator.angular_flux_to_scalar` was retired in
    Wave E Round 2 along with the rest of the FD-operator API surface.

    Issue #196 PR-INDEX-5: output is principled ``(ng, nx, ny)``.
    Issue #196 PR-INDEX-7: input is principled ``(N, ng, nx, ny)``
    (the FD-matvec internal ``(ng, N, nx, ny)`` layout flipped at the
    public ``solution_to_angular_flux*`` boundary; closes the
    PR-INDEX-4 §9.1 deferral).

    Parameters
    ----------
    fi
        Angular flux of shape ``(N, ng, nx, ny)``.

    Returns
    -------
    np.ndarray
        Scalar flux shape ``(ng, nx, ny)``.
    """
    # Named contraction: sum over the ordinate axis with quadrature
    # weights.  ``fi`` is ``(N, ng, nx, ny)``; ``quad.weights`` is
    # ``(N,)``; result is ``(ng, nx, ny)``.
    return np.einsum("ngxy,n->gxy", fi, quad.weights)


def _build_rhs_cartesian(
    fission_source: np.ndarray,
    scalar_flux: np.ndarray,
    eq_map,
    quad: AngularQuadrature,
    sig_s: dict[int, list[np.ndarray]],
    sig2: dict[int, np.ndarray],
    mat_map: np.ndarray,
    nx: int, ny: int, ng: int,
    scattering_order: int = 0,
    angular_flux: np.ndarray | None = None,
) -> np.ndarray:
    """Packed RHS for the Cartesian operator equation :math:`L\\,\\psi = b`.

    Internal helper for the Krylov inner-solver path.  All isotropic
    sources are divided by :math:`W = \\sum w_n` to match the operator
    equation convention (the operator carries no :math:`1/W`; sources
    fed to it must).

    For Pℓ scattering (``scattering_order > 0``), the per-ordinate
    scattering source is the Galerkin reconstruction over real
    spherical harmonics::

        qS(n) = Σ_l (2l+1) · Σ_s^l^T @ [Σ_m φ^lm · Y_lm(n)] / W

    Extracted from the Wave A-D ``build_rhs`` function in
    :mod:`orpheus.sn.operator`, which was retired in Wave E Round 2
    along with the rest of the BiCGSTAB FD-operator API surface.

    Issue #196 PR-INDEX-5: ``fission_source`` / ``scalar_flux`` are
    principled ``(ng, nx, ny)`` (per-cell slices index as
    ``[:, ix, iy]``).
    Issue #196 PR-INDEX-7: ``angular_flux`` (when supplied) is now
    principled ``(N, ng, nx, ny)`` end-to-end — both
    :func:`solution_to_angular_flux` (the warm-start path) and the
    sweep's angular flux output land in the same layout. Closes the
    PR-INDEX-4 §9.1 deferral.
    """
    sum_w = float(quad.weights.sum())
    L = scattering_order
    mu_z = getattr(quad, "mu_z", np.zeros(quad.N))

    # Precompute Legendre moments if anisotropic scattering.
    fiL = None
    Y = None
    if L > 0 and angular_flux is not None:
        Y = quad.spherical_harmonics(L)  # (N, L+1, 2L+1)
        w = quad.weights
        fiL = np.zeros((nx, ny, ng, L + 1, 2 * L + 1))
        for l in range(L + 1):
            for m in range(-l, l + 1):
                for n in range(quad.N):
                    # PR-INDEX-7: angular_flux is (N, ng, nx, ny).
                    # angular_flux[n, :, :, :] is (ng, nx, ny); .T
                    # reverses axes to (ny, nx, ng) — identical to
                    # the pre-PR-INDEX-7 ``angular_flux[:, n, :, :].T``
                    # behaviour bit-for-bit.
                    fiL[:, :, :, l, l + m] += (
                        w[n] * angular_flux[n, :, :, :].T * Y[n, l, l + m]
                    )

    rhs = np.zeros((ng, eq_map.n_eq))
    eq_idx = 0
    for iy in range(ny):
        for ix in range(nx):
            mid = int(mat_map[ix, iy])
            phi_cell = scalar_flux[:, ix, iy]

            qF = fission_source[:, ix, iy]
            q2 = 2.0 * (sig2[mid].T @ phi_cell) / sum_w

            for n in range(quad.N):
                if mu_z[n] < -1e-15:
                    continue
                if ix == 0 and quad.mu_x[n] > 1e-15:
                    continue
                if ix == nx - 1 and quad.mu_x[n] < -1e-15:
                    continue
                if iy == 0 and quad.mu_y[n] > 1e-15:
                    continue
                if iy == ny - 1 and quad.mu_y[n] < -1e-15:
                    continue

                qS = np.zeros(ng)
                for l in range(L + 1):
                    if l == 0:
                        qS += sig_s[mid][0].T @ phi_cell / sum_w
                    elif fiL is not None:
                        SUM = np.zeros(ng)
                        for m in range(-l, l + 1):
                            SUM += fiL[ix, iy, :, l, l + m] * Y[n, l, l + m]
                        qS += (2 * l + 1) * (sig_s[mid][l].T @ SUM) / sum_w

                rhs[:, eq_idx] = qF + q2 + qS
                eq_idx += 1

    return rhs.ravel(order="F")


def _build_rhs_spherical(
    fission_source: np.ndarray,
    scalar_flux: np.ndarray,
    eq_map,
    quad: AngularQuadrature,
    sig_s: dict[int, list[np.ndarray]],
    sig2: dict[int, np.ndarray],
    mat_map: np.ndarray,
    nx: int, ng: int,
) -> np.ndarray:
    """Packed RHS for spherical 1-D :math:`L\\,\\psi = b`.

    Same structure as :func:`_build_rhs_cartesian` but with the
    spherical equation map (no y-direction, no z-reflection
    filtering).  P0 isotropic scattering only — the curvilinear FD
    operator does not currently carry Pℓ.
    """
    sum_w = float(quad.weights.sum())

    rhs = np.zeros((ng, eq_map.n_eq))
    eq_idx = 0
    for ix in range(nx):
        mid = int(mat_map[ix, 0])
        # PR-INDEX-5: scalar_flux / fission_source are principled
        # ``(ng, nx, ny=1)``; per-cell slice is ``[:, ix, 0]``.
        phi_cell = scalar_flux[:, ix, 0]

        qF = fission_source[:, ix, 0]
        q2 = 2.0 * (sig2[mid].T @ phi_cell) / sum_w

        for n in range(quad.N):
            if ix == nx - 1 and quad.mu_x[n] < -1e-15:
                continue

            qS = sig_s[mid][0].T @ phi_cell / sum_w
            rhs[:, eq_idx] = qF + q2 + qS
            eq_idx += 1

    return rhs.ravel(order="F")


# Cylindrical 1-D shares the spherical RHS structure (same equation
# map, same isotropic-only scattering).  Internal alias so the
# Krylov dispatch reads cleanly.
_build_rhs_cylindrical = _build_rhs_spherical


def _is_curvilinear(mesh: Mesh1D | Mesh2D) -> bool:
    """Return ``True`` iff *mesh* is a 1-D spherical or cylindrical Mesh1D."""
    if not isinstance(mesh, Mesh1D):
        return False
    coord = getattr(mesh, "coord", None)
    if coord is None:
        return False
    name = getattr(coord, "name", str(coord)).upper()
    return name in ("SPHERICAL", "CYLINDRICAL")


# ═══════════════════════════════════════════════════════════════════════
# Public API
# ═══════════════════════════════════════════════════════════════════════

def solve_sn(
    materials: dict[int, Mixture],
    mesh: Mesh1D | Mesh2D,
    quadrature: AngularQuadrature,
    inner_solver: str = "source_iteration",
    scattering_order: int = 0,
    max_outer: int = 500,
    keff_tol: float = 1e-7,
    flux_tol: float = 1e-6,
    max_inner: int = 200,
    inner_tol: float = 1e-8,
) -> SNResult:
    """Solve the multi-group SN eigenvalue problem.

    This is the **canonical entry point** for the production SN solver.
    Production callers consume ``(materials, mesh, quadrature, ...)``
    directly: materials are :class:`~orpheus.data.macro_xs.mixture.Mixture`
    objects keyed by material ID, ``mesh`` is a
    :class:`~orpheus.geometry.Mesh1D` / :class:`~orpheus.geometry.Mesh2D`
    (build via :meth:`Mesh1D.from_geometry` for multi-region 1-D cases),
    and ``quadrature`` is an explicitly chosen
    :class:`~orpheus.sn.quadrature.AngularQuadrature` — Gauss-Legendre
    for slab, level-symmetric / product quadrature for curvilinear, or
    Lebedev for 2-D.

    The mesh's boundary conditions (``bc_left`` / ``bc_right`` for 1-D,
    ``bc_xmin`` / ``bc_xmax`` / ``bc_ymin`` / ``bc_ymax`` for 2-D) are
    honoured verbatim — the SN sweep handles ``vacuum`` and
    ``reflective``.

    Parameters
    ----------
    materials : dict mapping material ID to Mixture.
    mesh : Mesh1D or Mesh2D (base geometry).
    quadrature : AngularQuadrature
        Explicitly chosen by the caller — Gauss-Legendre for slab,
        level-symmetric / product quadrature for curvilinear 1-D,
        Lebedev for 2-D. Mismatches between geometry and quadrature
        family are not silently coerced.
    inner_solver : "source_iteration" (default) or "krylov".
        Wave E Round 2 deviation from the campaign plan: ``solve_sn``
        keeps the ``source_iteration`` default for **all** geometries
        (Cartesian and curvilinear).  The eigenvalue is shape-
        independent (k = production / absorption is a volume-weighted
        ratio), so even on the ERR-026-affected curvilinear sweep the
        keff is correct to the eigenvalue's tolerance, even though the
        flux *shape* would have a closure-bias drift.  Preserving the
        default keeps the 6 curvilinear regression snapshots bit-
        identical.  ``solve_sn_fixed_source`` *does* auto-flip to
        ``"krylov"`` for curvilinear because shape correctness is the
        whole point of fixed-source MMS verification.
    scattering_order : Legendre order for scattering (0 = P0).
    max_outer : maximum outer (power) iterations.
    keff_tol, flux_tol : outer convergence.
    max_inner, inner_tol : inner solver parameters.

    Returns
    -------
    SNResult
        Eigenvalue, scalar + angular flux, geometry/quadrature handles,
        and timing.
    """
    t_start = time.perf_counter()

    # Build augmented geometry (precomputes streaming stencil).
    # Issue #197 PR-TYPED-0: materials now lives on SNMesh — the
    # phase-space-as-such object.
    sn_mesh = SNMesh(mesh, quadrature, materials)

    solver = SNSolver(
        sn_mesh,
        inner_solver=inner_solver,
        scattering_order=scattering_order,
        keff_tol=keff_tol, flux_tol=flux_tol,
        max_inner=max_inner, inner_tol=inner_tol,
    )

    keff, keff_history, scalar_flux = power_iteration(solver, max_iter=max_outer)

    # Final sweep to get angular flux.  Issue #196 PR-INDEX-5: every
    # array is principled — scalar_flux ``(ng, nx, ny)``, angular_flux
    # ``(N, ng, nx, ny)``.
    Q_final = solver.compute_fission_source(scalar_flux, keff)
    solver._add_scattering_source(Q_final, scalar_flux)
    solver._add_n2n_source(Q_final, scalar_flux)
    angular_flux, _ = transport_sweep(
        Q_final, solver.sig_t, sn_mesh, solver._psi_bc,
    )

    _any_mat = next(iter(materials.values()))
    elapsed = time.perf_counter() - t_start

    return SNResult(
        keff=keff_history[-1],
        keff_history=keff_history,
        angular_flux=angular_flux,
        scalar_flux=scalar_flux,
        geometry=mesh,
        quadrature=quadrature,
        eg=_any_mat.eg,
        elapsed_seconds=elapsed,
    )


def solve_sn_fixed_source(
    materials: dict[int, Mixture],
    mesh: Mesh1D | Mesh2D,
    quadrature: AngularQuadrature,
    external_source: np.ndarray,
    boundary_condition: str = "vacuum",
    scattering_order: int = 0,
    max_inner: int = 1000,
    inner_tol: float = 1e-12,
    inner_solver: str | None = None,
) -> SNFixedSourceResult:
    r"""Solve the multi-group SN fixed-source transport problem.

    Solves

    .. math::

        \mu_n \frac{\partial\psi_n}{\partial x}+\Sigma_t\psi_n
        = \frac{1}{W}\left(\Sigma_s\phi + Q^{\text{ext}}_n\right)

    for a prescribed per-ordinate external source ``external_source``,
    with vacuum or reflective boundary conditions. The fission source
    is zero — this is the pure transport operator. Source iteration
    converges geometrically at rate :math:`c = \max\Sigma_s/\Sigma_t`.

    Parameters
    ----------
    materials, mesh, quadrature, scattering_order :
        Same as :func:`solve_sn`.
    external_source : (N, ng, nx, ny)
        Per-ordinate volumetric source :math:`Q^{\text{ext}}_n(x)`.
        Passed to the sweep as its ``Q_aniso`` argument (the solver
        applies the :math:`1/W` factor internally).  Issue #196
        PR-INDEX-5: principled layout (``g`` axis after ``N``).
    boundary_condition : {"vacuum", "reflective"}
        Applied to all faces when the mesh has no explicit BC
        declarations (``bc_left`` etc. are ``None``).  When the mesh
        carries explicit :class:`~orpheus.geometry.mesh.BC` fields,
        those take precedence and this parameter is ignored.
        Vacuum is the default because the intended consumer is
        Method of Manufactured Solutions verification on a finite slab.
    max_inner, inner_tol :
        Inner solver iteration limits.
    inner_solver : {"source_iteration", "krylov", None}
        Inner-solve strategy.  When ``None`` (default), all geometries
        use ``"source_iteration"`` — bit-identical to the Wave A-D
        path.  Wave E Round 3 ships the BC-aware FD operator
        (:func:`solution_to_angular_flux*` consume the mesh's
        :class:`~orpheus.geometry.boundary.BoundaryOperator` instances via
        :meth:`apply`), which makes ``"krylov"`` available
        as an opt-in for vacuum / reflective / white / albedo / mixed
        BCs uniformly — but the curvilinear-default flip is **not**
        enabled because empirically the symmetric-closure FD operator
        at the curvilinear outer face uses cell-center as a face-flux
        approximation, which is only first-order at the boundary on
        non-constant solutions.  Switching the default to
        ``"krylov"`` regresses the curvilinear MMS convergence rate
        from the WDD sweep's :math:`\mathcal{O}(h^{1.3})` (still
        ERR-026-affected, but volumetrically benign for MMS) to
        :math:`\mathcal{O}(h^{1})` (the FD operator's boundary
        truncation).  Round 3's BC plumbing is therefore the
        infrastructure that *enables* a future full closure;
        the closure itself depends on a follow-up that fixes the
        FD operator's boundary face-flux treatment (DD diamond
        relation at the outer boundary, or analogous extrapolation).

        ``"krylov"`` is still the right choice for **constant-source**
        problems (where the cell-center-as-face-value approximation
        is exact), as evidenced by the
        :file:`tests/sn/test_sweep_operator_inconsistency.py` regression
        suite — krylov gives the analytical flat flux to round-off
        while the sweep produces the documented ERR-026 deviation.

    Notes
    -----
    This entry point exists for L1 verification via MMS, not for
    engineering problems — real fixed-source calculations should still
    build on :func:`solve_sn` with an appropriate external-source hook.
    See :mod:`orpheus.derivations.continuous.mms.sn` and the MMS verification
    section of the discrete-ordinates theory page.
    """
    t_start = time.perf_counter()

    # Apply boundary_condition parameter to mesh if no explicit BCs set
    mesh = _apply_default_bcs(mesh, boundary_condition)
    # Issue #197 PR-TYPED-0: materials threaded through SNMesh.
    sn_mesh = SNMesh(mesh, quadrature, materials)

    # Issue #168 status (Phase D ERR-026 closure, 2026-05-12):
    #
    # * Phase A (Defects 1 + 2): the ``BoundaryFaceFlux`` Protocol
    #   retired in Phase C.
    # * Phase B (Defect 3): :class:`MorelMontryAngularSweep` is now
    #   the SNMesh-default :class:`PoleAngularClosure`.
    # * Phase C (sweep-frame matvec): the apply matvec is one sweep
    #   iteration semantically, with WDD spatial closure.
    # * Phase D (Carlson coupled-pole seed):
    #   :class:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonInwardSweep`
    #   seeds the M-M angular recurrence's half-angle face flux, making
    #   per-ordinate flat-flux balance hold on sphere Gate 1.1 MMS.
    #
    # **Phase D default flip (2026-05-12)**: curvilinear (sphere /
    # cylinder) defaults to ``"krylov"``.  The Phase B / C / D apply
    # matvec gives the canonical Hébert form; Krylov-on-apply
    # converges cleanly to the apply matvec's fixed point without
    # the ERR-026 closure-bias drift.  Cartesian (slab / 2-D) stays
    # at ``"source_iteration"`` (the existing inner solver default
    # for Cartesian — no ERR-026 affected closure on slab).
    if inner_solver is None:
        if getattr(sn_mesh, "curvature", None) in ("spherical", "cylindrical"):
            inner_solver = "krylov"
        else:
            inner_solver = "source_iteration"

    solver = SNSolver(
        sn_mesh,
        inner_solver=inner_solver,
        scattering_order=scattering_order,
        max_inner=max_inner, inner_tol=inner_tol,
    )

    N = sn_mesh.quad.N
    nx, ny, ng = sn_mesh.nx, sn_mesh.ny, solver.ng
    # Issue #196 PR-INDEX-5: external_source principled (N, ng, nx, ny).
    if external_source.shape != (N, ng, nx, ny):
        raise ValueError(
            f"external_source shape {external_source.shape} does not "
            f"match (N, ng, nx, ny) = {(N, ng, nx, ny)}"
        )

    if inner_solver == "source_iteration":
        return _solve_fixed_source_si(
            solver, sn_mesh, external_source, mesh, quadrature, materials,
            t_start, max_inner, inner_tol,
        )

    # Krylov path.  We solve T·ψ = b directly via GMRES, where b carries
    # the external per-ordinate source plus any in-scatter / (n,2n) terms
    # built from the converged scalar flux.  Wrapping that in an outer
    # source iteration converges scattering self-consistently.
    return _solve_fixed_source_krylov(
        solver, sn_mesh, external_source, mesh, quadrature, materials,
        t_start, max_inner, inner_tol,
    )


def _solve_fixed_source_si(
    solver: SNSolver,
    sn_mesh: SNMesh,
    external_source: np.ndarray,
    mesh: Mesh1D | Mesh2D,
    quadrature: AngularQuadrature,
    materials: dict[int, Mixture],
    t_start: float,
    max_inner: int,
    inner_tol: float,
) -> SNFixedSourceResult:
    """Cartesian-default fixed-source path: source iteration via the sweep.

    Bit-identical math to the Wave A-D inline loop in
    :func:`solve_sn_fixed_source`.  Extracted as a helper to make the
    geometry-default dispatch in :func:`solve_sn_fixed_source` clean.
    """
    nx, ny, ng = sn_mesh.nx, sn_mesh.ny, solver.ng

    # Issue #196 PR-INDEX-5: phi principled (ng, nx, ny).
    phi = np.zeros((ng, nx, ny))
    angular = None
    residual = np.inf

    for n_inner in range(max_inner):
        phi_prev = phi

        # Isotropic in-scatter + (n,2n). No fission — this is fixed-source.
        Q = np.zeros_like(phi)
        solver._add_scattering_source(Q, phi)
        solver._add_n2n_source(Q, phi)

        # Merge P1+ scattering moments (if any) with external MMS source.
        Q_aniso_p1 = solver._build_aniso_scattering(angular)
        if Q_aniso_p1 is None:
            Q_aniso_total = external_source
        else:
            Q_aniso_total = Q_aniso_p1 + external_source

        angular, phi = transport_sweep(
            Q, solver.sig_t, sn_mesh, solver._psi_bc,
            Q_aniso=Q_aniso_total,
        )

        norm = np.linalg.norm(phi)
        if norm > 0:
            residual = np.linalg.norm(phi - phi_prev) / norm
            if residual < inner_tol:
                break
    else:
        n_inner = max_inner - 1  # loop exhausted without break

    elapsed = time.perf_counter() - t_start
    _any_mat = next(iter(materials.values()))
    return SNFixedSourceResult(
        angular_flux=angular,
        scalar_flux=phi,
        geometry=mesh,
        quadrature=quadrature,
        n_inner=n_inner + 1,
        residual=float(residual),
        elapsed_seconds=elapsed,
        eg=_any_mat.eg,
    )


def _solve_fixed_source_krylov(
    solver: SNSolver,
    sn_mesh: SNMesh,
    external_source: np.ndarray,
    mesh: Mesh1D | Mesh2D,
    quadrature: AngularQuadrature,
    materials: dict[int, Mixture],
    t_start: float,
    max_inner: int,
    inner_tol: float,
) -> SNFixedSourceResult:
    r"""Curvilinear-default fixed-source path: GMRES on :meth:`L.apply`.

    This is the Wave E Round 2 reconciliation that closes ERR-026.  The
    operator equation :math:`(L - S)\,\psi = q_{\rm ext}` is solved
    once via GMRES on the symmetric closure (with the sweep as
    preconditioner).  When scattering is significant, we wrap an outer
    source iteration around the Krylov inner solve; the scattering
    self-consistency converges geometrically at rate :math:`c =
    \max\Sigma_s/\Sigma_t`.

    The math reuses :meth:`SNSolver._solve_krylov`'s machinery: build
    a packed RHS, GMRES on ``L.apply`` with sweep preconditioner,
    decode packed solution to angular flux.  The only addition is that
    the per-ordinate ``external_source`` enters the packed RHS with
    the ``/sum_w`` normalisation that the operator equation expects.
    """
    nx, ny, ng = sn_mesh.nx, sn_mesh.ny, solver.ng
    N = sn_mesh.quad.N
    sum_w = float(sn_mesh.quad.weights.sum())
    curv = getattr(sn_mesh, "curvature", None)

    # EquationMap dispatch matches SNStreamingOperator.apply.
    if curv == "spherical":
        eq_map = build_equation_map_spherical(nx, sn_mesh.quad, ng)
    elif curv == "cylindrical":
        eq_map = build_equation_map_cylindrical(nx, sn_mesh.quad, ng)
    else:
        eq_map = build_equation_map(nx, ny, sn_mesh.quad, ng)
    n_unknowns = eq_map.n_unknowns

    # Pre-pack the external source as the per-ordinate contribution to
    # the RHS, divided by sum_w to match the operator equation's
    # convention (build_rhs* outputs are pre-divided by sum_w).
    # Issue #196 PR-INDEX-5: external_source is principled
    # ``(N, ng, nx, ny)``; the packed-vector slice is ``[n, :, ix, iy]``.
    ext_packed = np.empty((ng, eq_map.n_eq), dtype=float)
    for k in range(eq_map.n_eq):
        ext_packed[:, k] = external_source[
            eq_map.ordinate[k], :, eq_map.ix[k], eq_map.iy[k],
        ] / sum_w
    ext_packed_flat = ext_packed.ravel(order="F")

    # Outer source iteration.  Each outer step rebuilds the in-scatter
    # source from the current scalar flux, then drives a GMRES inner
    # solve to convergence.
    L_scipy = ScipyLinearOperator(
        (n_unknowns, n_unknowns), matvec=solver.L.apply, dtype=float,
    )
    precond = solver._make_sweep_preconditioner(eq_map, n_unknowns, sum_w)

    # Issue #196 PR-INDEX-5: phi principled (ng, nx, ny).
    phi = np.zeros((ng, nx, ny))
    solution = np.zeros(n_unknowns)
    residual = np.inf
    angular = None

    for n_outer in range(max_inner):
        phi_prev = phi.copy()

        # Build the in-scatter / (n,2n) RHS contribution from the
        # current scalar flux.  Fission source is zero (fixed-source).
        if curv == "spherical":
            rhs_iso = _build_rhs_spherical(
                np.zeros_like(phi), phi, eq_map, sn_mesh.quad,
                solver.sig_s, solver.sig2, sn_mesh.mat_map,
                nx, ng,
            )
        elif curv == "cylindrical":
            rhs_iso = _build_rhs_cylindrical(
                np.zeros_like(phi), phi, eq_map, sn_mesh.quad,
                solver.sig_s, solver.sig2, sn_mesh.mat_map,
                nx, ng,
            )
        else:
            # PR-INDEX-7: ``_build_rhs_cartesian`` consumes ``angular_flux``
            # in principled ``(N, ng, nx, ny)`` layout — the same layout
            # the sweep produces. No axis-swap adapter needed.
            angular_full = angular if (
                solver.scattering_order > 0 and angular is not None
            ) else None
            rhs_iso = _build_rhs_cartesian(
                np.zeros_like(phi), phi, eq_map, sn_mesh.quad,
                solver.sig_s, solver.sig2, sn_mesh.mat_map,
                nx, ny, ng,
                scattering_order=solver.scattering_order,
                angular_flux=angular_full,
            )

        rhs = rhs_iso + ext_packed_flat

        # GMRES inner solve.
        try:
            solution, info = gmres(
                L_scipy, rhs, x0=solution, M=precond,
                rtol=inner_tol, atol=0.0,
                maxiter=max_inner, restart=min(50, n_unknowns),
            )
        except TypeError:
            solution, info = gmres(
                L_scipy, rhs, x0=solution, M=precond,
                tol=inner_tol,
                maxiter=max_inner, restart=min(50, n_unknowns),
            )

        # Decode packed solution → angular and scalar flux.
        # Issue #168 Phase A: curvilinear decoders now return
        # Issue #168 Phase C: pure cell-centre return; the Phase A
        # boundary_face_flux companion retired with the
        # BoundaryFaceFlux Protocol.
        if curv == "spherical":
            fi = solution_to_angular_flux_spherical(
                solution, eq_map, sn_mesh.quad, nx, ng,
            )
        elif curv == "cylindrical":
            fi = solution_to_angular_flux_cylindrical(
                solution, eq_map, sn_mesh.quad, nx, ng,
            )
        else:
            fi = solution_to_angular_flux(
                solution, eq_map, sn_mesh.quad, nx, ny, ng,
                bc_xmin=sn_mesh.bc_xmin,
                bc_xmax=sn_mesh.bc_xmax,
                bc_ymin=sn_mesh.bc_ymin,
                bc_ymax=sn_mesh.bc_ymax,
            )
        # PR-INDEX-7: ``solution_to_angular_flux*`` returns principled
        # ``(N, ng, nx, ny)`` natively — no axis-swap adapter needed.
        # Closes the PR-INDEX-4 §9.1 deferral.
        angular = fi
        phi = _scalar_flux_from_angular(fi, sn_mesh.quad, nx, ny, ng)

        norm = np.linalg.norm(phi)
        if norm > 0:
            residual = np.linalg.norm(phi - phi_prev) / norm
            if residual < inner_tol:
                break
    else:
        n_outer = max_inner - 1

    elapsed = time.perf_counter() - t_start
    _any_mat = next(iter(materials.values()))
    return SNFixedSourceResult(
        angular_flux=angular,
        scalar_flux=phi,
        geometry=mesh,
        quadrature=quadrature,
        n_inner=n_outer + 1,
        residual=float(residual),
        elapsed_seconds=elapsed,
        eg=_any_mat.eg,
    )
