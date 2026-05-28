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
from dataclasses import replace

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
from orpheus.numerics.quadrature import Quadrature
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


# Issue #197 PR-TYPED-5: SNFixedSourceResult + SNResult RETIRED.
# Both solver entry points now return a typed :class:`Solution`
# (orpheus.sn.solution); the two legacy data bags collapse into one
# (coding-elegance Pattern 7 — single source of truth) with the
# discrimination living in Solution.is_eigenvalue() / is_fixed_source().
from .solution import IterationHistory, Solution


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

        # Persistent boundary flux state (passed to sweep).
        # Issue #197 PR-TYPED-2: typed :class:`BoundaryFlux` replaces
        # the stringly-typed ``psi_bc: dict``.  Per-face buffers
        # become named attributes; typos surface as AttributeError.
        self._boundary_flux = sn_mesh.zeros_boundary_flux()

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

    # Issue #197 PR-TYPED-2: the PR-TYPED-1 transient read-through
    # shims (``sig_t``, ``sig_a``, ``sig_p``, ``chi``, ``sig_s``,
    # ``sig2``, ``sig_s0``, ``_cells_by_mat``) have been RETIRED.
    # Every downstream consumer now reads ``self.mat_xs.*`` directly:
    # ``mat_xs.total_cross_section``,
    # ``mat_xs.absorption_cross_section``,
    # ``mat_xs.fission_production``,
    # ``mat_xs.emission_spectrum``,
    # ``mat_xs.sig_s_legendre(mid)``,
    # ``mat_xs.n2n_matrix(mid)``,
    # ``mat_xs.cells_by_material``.

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
        # Issue #196 PR-INDEX-5: both ``mat_xs.fission_production`` and
        # ``flux_distribution`` are principled ``(ng, nx, ny)``.  The
        # named intermediate ``per_cell_per_group`` has units ``[1/s]``
        # per cell per group (a reactor-physics quantity — coding-
        # elegance Pattern 3).  ``volume_measure`` consumes a flat
        # ``(N_cells, ng)`` view (Issue 9.6 wiring); reshape via
        # ``transpose(1, 2, 0)`` then flatten the spatial axes.
        # Issue #197 PR-TYPED-2: consumes mat_xs directly (no shim).
        per_cell_per_group = np.einsum(
            "gxy,gxy->gxy", self.mat_xs.fission_production, flux_distribution,
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
        # Issue #197 PR-TYPED-2: consumes mat_xs directly (no shim).
        per_cell_per_group = np.einsum(
            "gxy,gxy->gxy", self.mat_xs.absorption_cross_section, flux_distribution,
        )
        return mu(per_cell_per_group.transpose(1, 2, 0).reshape(nx * ny, ng))

    def compute_production_rate(self, flux_distribution: np.ndarray) -> float:
        r"""Total volume-integrated neutron production rate (scalar).

        :math:`P(\phi) = \sum_g r_g = \sum_g \int_V \nu \Sigma_{f,g} \phi_g\,dV
        + 2 \sum_g \int_V \sum_{g'} \Sigma_{2,g'\to g} \phi_{g'} \,dV`,
        i.e. ``.sum()`` of :meth:`compute_group_production_rate`.

        This is the canonical scale anchor for the SN eigenmode:
        :func:`orpheus.numerics.eigenvalue.power_iteration` renormalises
        :math:`\phi` to unit production rate at each outer step
        (ERR-052), which gives callers a stable physical handle —
        scaling to an absolute flux at a target reactor power becomes a
        single multiplication by :math:`P_{\text{target}} / \kappa`
        where :math:`\kappa` is the energy released per fission.
        """
        return float(self.compute_group_production_rate(flux_distribution).sum())

    def compute_keff(self, flux_distribution: np.ndarray) -> float:
        """k = production / absorption (volume-weighted).

        Composed from the per-group production and absorption rate
        vectors so the intermediates are individually meaningful and
        reusable (e.g. spectral diagnostics).  See
        :meth:`compute_group_production_rate` and
        :meth:`compute_group_absorption_rate`.
        """
        production = self.compute_production_rate(flux_distribution)
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
        r"""Inner within-group solve via :class:`SourceIteration` on typed AngularFlux.

        R-1 Step E (2026-05-19) — carved onto
        :class:`~orpheus.numerics.iteration.SourceIteration` consuming
        the same typed-flux operator triple as the Krylov carve:

        .. math::

            A \;=\; L + C\,, \quad
            S \;=\; \tfrac{1}{W}\,\text{full multi-group scatter}\,, \quad
            F \;=\; 0_{\rm wg}

        on :class:`~orpheus.sn.angular_flux.AngularFlux`.  The
        iteration step is

        .. math::

            \psi_{n+1} \;=\; (L + C)^{-1}\,\bigl(S\,\psi_n
                                                + F\,\psi_n
                                                + q_{\rm ext}\bigr)

        where ``(L + C).solve`` IS the WDD sweep
        (:class:`~orpheus.sn.operator.InvertibleOperator` from R-1
        Step C).  Phase 1.2 — the previous iterate :math:`\psi_n`
        travels into the sweep via the explicit ``initial_guess``
        kwarg on :meth:`InvertibleOperator.solve`; the M-M closure's
        ``psi_half_seed`` strategy reads it to derive the curvilinear
        Carlson coupled-pole seed.  ``SourceIteration._solve_with_seed``
        bridges the call.

        Scope
        =====

        1-D only (slab + sphere + cylinder).  2-D Cartesian raises
        :class:`NotImplementedError` (Phase A absorbs the 2-D B1''
        face layout).
        """
        if self.sn_mesh.reduced is None:
            raise NotImplementedError(
                "R-1 Step E — 2-D Cartesian SI carve deferred to "
                "Phase A.  The B1'' face block is 1-D-only; 2-D needs "
                "a separate 4-face layout (xmin, xmax, ymin, ymax)."
            )

        from .angular_flux import AngularFlux
        from .operator import (
            CollisionOperator, StreamingOperator,
            _copy_boundary_face_state,
        )
        from orpheus.transport.sources import PerOrdinateSource
        from orpheus.numerics.iteration import SourceIteration
        from orpheus.numerics.operator import ZeroOperator

        # ── Build the typed RHS ──────────────────────────────────────
        # R-1 Step 4 A1 — q_ext as per-ordinate density via the canonical
        # ``PerOrdinateSource.from_isotropic`` factory.  The /W projection
        # lives at the factory boundary (Pattern 7 producer-side
        # normalisation); the legacy
        # ``q_per_ord = (fission_source / sum_w)[None, :, :, :]``
        # broadcast pattern is GONE.
        q_ext_per_ord = PerOrdinateSource.from_isotropic(
            fission_source, self.sn_mesh,
        )
        q_ext_typed = AngularFlux(q_ext_per_ord.values, self.sn_mesh)
        # Pre-populate q_ext.boundary with the persistent partner-flux
        # state ``self._boundary_flux`` from the previous outer call.
        # This is the load-bearing plumbing for reflective BCs: the
        # FIRST inner SI iteration has no ``initial_guess`` (cold
        # start), so without this seeding the sweep would see a zero
        # BC trace (= vacuum), and the SI fixed point shifts away from
        # the true reflective answer.  Subsequent iterations update
        # the partner flux via the ``initial_guess.boundary`` thread
        # that :meth:`InvertibleOperator.solve` reads.
        _copy_boundary_face_state(self._boundary_flux, q_ext_typed.boundary)

        # ── Build the typed operator triple ─────────────────────────
        L_leaf = StreamingOperator(
            self.sn_mesh, self.mat_xs.total_cross_section,
        )
        C_t = CollisionOperator(
            self.sn_mesh, self.mat_xs.total_cross_section,
        )
        LC = L_leaf + C_t  # InvertibleOperator — apply + solve

        # R-1 Step 4 A1 — ``ScatteringOperator.apply`` now returns
        # per-ordinate density at the producer boundary (the
        # ``/sum_w`` projection lives inside the singledispatched
        # ``apply``).  The legacy consumer-side rescale
        # ``S_normalised = self.scattering_op / sum_w`` is GONE.

        si = SourceIteration(
            LC,
            self.scattering_op,
            ZeroOperator(),
            max_iter=self.max_inner,
            tol=self.inner_tol,
        )

        # ── Warm start (typed) ──────────────────────────────────────
        # SourceIteration._solve_with_seed forwards the previous
        # iterate to InvertibleOperator.solve via the explicit
        # ``initial_guess`` kwarg (Phase 1.2 — the M-M closure's
        # ``psi_half_seed`` strategy reads it to derive the Carlson
        # coupled-pole seed; the previous iterate's boundary trace
        # seeds the reflective-BC partner-flux state).
        initial_guess = getattr(self, "_psi_typed", None)

        psi_typed, _residuals = si.solve(
            q_ext_typed, initial_guess=initial_guess,
        )
        self._psi_typed = psi_typed

        # Reduce angular → scalar flux for the eigenvalue outer's contract.
        return psi_typed.integrate_angular().values

    # ── Inner solver: Krylov on (L+C-S)·ψ = q_ext (R-1 Step D carve) ──

    def _solve_krylov(
        self, fission_source: np.ndarray, flux_distribution: np.ndarray,
    ) -> np.ndarray:
        r"""Inner solve via GMRES on :math:`(L + C - S)\,\psi = q_{\rm ext}`.

        R-1 Step D (2026-05-19) — carved onto
        :class:`~orpheus.numerics.iteration.KrylovAcceleration` consuming
        the operator triple

        .. math::

            A \;=\; L + C\,, \quad
            S \;=\; \text{full multi-group scatter}\,, \quad
            F \;=\; 0_{\rm wg}

        on typed :class:`~orpheus.sn.angular_flux.AngularFlux`.  The
        composite ``L + C`` returns an
        :class:`~orpheus.sn.operator.InvertibleOperator` (R-1 Step C);
        its ``.solve`` IS the WDD sweep but R-1 ships GMRES
        UNPRECONDITIONED (``preconditioner=None``) per user direction
        ("consolidating the foundational architecture; the block-inverse
        face preconditioner is `issue #200
        <https://github.com/deOliveira-R/ORPHEUS/issues/200>`_").  The
        sweep-as-preconditioner reactivation lives there.

        The within-group fission ``F`` is zero — the fission source
        comes in as the external ``q_{\rm ext}`` per the eigenvalue
        outer/within-group inner decomposition (Lewis & Miller §6.4).

        Scope
        =====

        1-D meshes only (slab + sphere + cylinder).  2-D Cartesian
        raises :class:`NotImplementedError` — the typed-flux B1''
        face block was designed 1-D-only (one ``xmin_face`` + one
        ``xmax_face``); 2-D needs a 4-face layout that's deferred to
        Phase A.

        Returns the updated scalar flux ``(ng, nx, ny)``.
        """
        if self.sn_mesh.reduced is None:
            raise NotImplementedError(
                "R-1 Step D — 2-D Cartesian Krylov carve deferred to "
                "Phase A.  The B1'' face block is 1-D-only; 2-D needs a "
                "separate 4-face layout (xmin, xmax, ymin, ymax).  Use a "
                "1-D mesh (slab / sphere / cylinder) or switch "
                "inner_solver='source_iteration'."
            )

        from .angular_flux import AngularFlux
        from .operator import CollisionOperator, StreamingOperator
        from orpheus.transport.sources import PerOrdinateSource
        from orpheus.numerics.iteration import KrylovAcceleration
        from orpheus.numerics.operator import ZeroOperator

        # ── Build the typed RHS ──────────────────────────────────────
        # R-1 Step 4 A1 — q_ext as per-ordinate density via the canonical
        # ``PerOrdinateSource.from_isotropic`` factory.  /W lives at the
        # factory boundary (Pattern 7 producer-side normalisation).
        q_ext_per_ord = PerOrdinateSource.from_isotropic(
            fission_source, self.sn_mesh,
        )
        q_ext_typed = AngularFlux(q_ext_per_ord.values, self.sn_mesh)

        # ── Build the typed operator triple ─────────────────────────
        # InvertibleOperator via __add__ dispatch (R-1 Step C).
        L_leaf = StreamingOperator(
            self.sn_mesh, self.mat_xs.total_cross_section,
        )
        C_t = CollisionOperator(
            self.sn_mesh, self.mat_xs.total_cross_section,
        )
        LC = L_leaf + C_t  # InvertibleOperator: apply + solve

        # R-1 Step 4 A1 — ``ScatteringOperator.apply`` returns
        # per-ordinate density at the producer boundary.  The legacy
        # consumer-side rescale ``S_normalised = self.scattering_op /
        # sum_w`` is GONE.

        # NOTE on the preconditioner — R-1 ships GMRES UNPRECONDITIONED
        # (explicit identity) per user direction.  Issue #200 tracks the
        # block-inverse face preconditioner re-enablement.  Post-Phase-1.2
        # the silent-fallback bug that motivated the explicit identity is
        # GONE: ``InvertibleOperator.solve`` takes an explicit
        # ``initial_guess`` kwarg (the M-M Carlson seed reads it), and
        # :class:`KrylovAcceleration`'s ``preconditioner=None`` default
        # invokes ``L.solve(q)`` with no ``initial_guess`` → cold-start
        # M-M seed (deterministic, no garbage).  Sweep-as-preconditioner
        # is functionally safe to re-enable.  The explicit identity is
        # kept here to preserve the L1 anchor's bit-identity until #200
        # ships the production preconditioner choice.
        N = self.quad.N
        nx, ny, ng = self.sn_mesh.nx, self.sn_mesh.ny, self.ng
        krylov = KrylovAcceleration(
            LC,
            self.scattering_op,
            ZeroOperator(),
            preconditioner=lambda q: q,  # explicit identity — issue #200 tracks re-enable
            tol=self.inner_tol,
            max_iter=self.max_inner,
            restart=min(50, N * ng * nx * ny),
        )

        # ── Warm start (typed) ──────────────────────────────────────
        initial_guess = getattr(self, "_psi_typed", None)

        psi_typed, _residuals = krylov.solve(
            q_ext_typed, initial_guess=initial_guess,
        )
        self._psi_typed = psi_typed

        # Reduce angular → scalar flux for the eigenvalue outer's contract.
        return psi_typed.integrate_angular().values

    def _make_sweep_preconditioner(
        self, eq_map, n: int,
    ) -> ScipyLinearOperator:
        r"""Build a scipy LinearOperator wrapping the sweep as :math:`L^{-1}`.

        The sweep takes a single :class:`PerOrdinateSource` (R-1 Step 4
        A1) and returns ``(angular_flux, scalar_flux)``.  GMRES wants a
        scalar ``matvec(q) -> M^{-1}·q`` on the packed 1-D vector.  This
        adapter:

        1. Decodes the packed RHS into a per-ordinate source carrier of
           shape ``(N, ng, nx, ny)`` (principled storage per
           :ref:`theory-sn-index-convention`).  Both the packed-vector
           convention AND the sweep's expected convention are
           per-ordinate density after R-1 Step 4 A1, so the decoded
           values feed the sweep directly with no rescale.
        2. Runs the sweep with the single per-ordinate source.
        3. Re-packs the resulting angular flux into the packed
           solution-vector layout via the inverse of
           :func:`solution_to_angular_flux*`.

        The sweep's internal :class:`BoundaryFlux` is NOT shared with
        :attr:`self._boundary_flux` — the preconditioner is stateless
        across GMRES inner iterations, which keeps the linear-operator
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
            # ``fi_op`` is the principled (N, ng, nx, ny) angular flux
            # decoded from the packed RHS.  ``build_rhs*`` pre-applied
            # ``/sum_w`` to the operator-equation source, so ``fi_op``
            # IS in per-ordinate density magnitude.
            #
            # R-1 Step 4 A1 — the legacy ``* sum_w`` bridge that undid
            # the sweep's old internal ``/W`` is GONE; the sweep now
            # consumes per-ordinate density directly via a single
            # :class:`PerOrdinateSource` parameter.
            from orpheus.transport.sources import PerOrdinateSource
            source = PerOrdinateSource.from_mesh(fi_op, self.sn_mesh)
            boundary_flux_local = self.sn_mesh.zeros_boundary_flux()
            try:
                angular, _ = transport_sweep(
                    source, self.mat_xs.total_cross_section,
                    self.sn_mesh, boundary_flux_local,
                )
            except Exception:
                # If the sweep cannot run with this aniso shape, degrade
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


# ─────────────────────────────────────────────────────────────────────
# Retired in P1.7 (moment-space + layering plan):
#   `_build_rhs_cartesian`, `_build_rhs_spherical`,
#   `_build_rhs_cylindrical`  (the latter was an alias of spherical)
#
# All three were packed-1D RHS builders for the legacy BiCGSTAB
# FD-operator path.  Production call graph (verified by
# `tests/sn/test_fixed_source_g1.py::TestNoLegacyMachineryInCallPath::
# test_no_legacy_eq_map_or_decoder_in_g1_path`): ZERO callers — the
# G1 Krylov path goes through typed `KrylovAcceleration` +
# `(L + C).apply` instead of building a packed-1D RHS.
#
# The cartesian variant carried an inline `(2*l+1)` literal at its
# Pℓ-source build, duplicating the canonical addition-theorem factor
# now sourced from
# `SphericalHarmonicSpace.addition_theorem_factor` via
# `HarmonicMomentReconstruction.from_spherical_harmonic_space`.  Per
# the moment-space plan §P1.3 "exactly one place" claim, retiring
# this dead duplicate is required.
#
# The spherical/cylindrical variants had no (2*l+1) duplicate (they
# were P0-isotropic-only), but they share the same dead-code status
# and retire alongside cartesian for cleanliness — superseded code is
# noise per `feedback_aggressive_retirement`.
# ─────────────────────────────────────────────────────────────────────


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
) -> Solution:
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
    Solution
        Typed return carrying eigenvalue, typed
        :class:`~orpheus.sn.angular_flux.AngularFlux` +
        :class:`~orpheus.sn.scalar_flux.ScalarFlux` +
        :class:`~orpheus.sn.boundary_flux.BoundaryFlux` fields plus an
        :class:`~orpheus.sn.solution.IterationHistory` carrying the
        eigenvalue trajectory.  Issue #197 PR-TYPED-5 — the legacy
        :class:`SNResult` data bag was retired in favour of the
        unified :class:`Solution` type that covers both eigenvalue and
        fixed-source problems.
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
    # R-1 Step 4 A1: typed source via the canonical
    # :meth:`PerOrdinateSource.from_isotropic` factory (Pattern 7
    # producer-side normalisation — /W projection at the factory
    # boundary).
    from orpheus.transport.sources import PerOrdinateSource
    Q_final = solver.compute_fission_source(scalar_flux, keff)
    solver._add_scattering_source(Q_final, scalar_flux)
    solver._add_n2n_source(Q_final, scalar_flux)
    angular_flux, _ = transport_sweep(
        PerOrdinateSource.from_isotropic(Q_final, sn_mesh),
        solver.mat_xs.total_cross_section, sn_mesh,
        solver._boundary_flux,
    )

    elapsed = time.perf_counter() - t_start

    # Issue #197 PR-TYPED-5: build typed Solution.  The bare ndarrays
    # produced by power_iteration + the final sweep get wrapped in
    # AngularFlux + ScalarFlux at this single boundary; downstream
    # consumers read typed fields.
    # D-H.1b (2026-05-28): Solution.angular_flux is now a typed
    # TimedFullField composite (bulk + boundary + history). Build
    # via the legacy-AngularFlux adapter — solver._boundary_flux is
    # still a legacy BoundaryFlux instance; the L2 BoundaryFlux is
    # constructed inside TimedFullField.from_legacy_angular_flux via
    # the legacy → L2 adapter.
    from .angular_flux import AngularFlux
    from orpheus.transport.fields.scalar_flux import ScalarFlux
    from orpheus.transport.timed_full_field import TimedFullField
    history = IterationHistory(
        keff_history=tuple(keff_history),
        n_outer=len(keff_history),
        converged=True,
    )
    legacy_psi = AngularFlux(
        angular_flux, sn_mesh, boundary=solver._boundary_flux,
    )
    return Solution(
        angular_flux=TimedFullField.from_legacy_angular_flux(legacy_psi),
        scalar_flux=ScalarFlux.from_mesh(scalar_flux, sn_mesh),
        mesh=sn_mesh,
        keff=float(keff_history[-1]),
        history=history,
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
) -> Solution:
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
        Per-ordinate volumetric source :math:`Q^{\text{ext}}_n(x)` in
        **per-ordinate density magnitude** (R-1 Step 4 A1 convention).
        Callers with an iso scalar source :math:`Q(\vec r, g)` should
        project to per-ordinate via
        :meth:`~orpheus.sn.sources.PerOrdinateSource.from_isotropic`
        before passing (the :math:`1/W` projection lives at the
        producer boundary per Pattern 7).  The sweep does NOT apply
        ``/W`` internally.  Issue #196 PR-INDEX-5: principled layout
        (``g`` axis after ``N``).
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
) -> Solution:
    """Cartesian-default fixed-source path: source iteration via the sweep.

    Bit-identical math to the Wave A-D inline loop in
    :func:`solve_sn_fixed_source`.  Extracted as a helper to make the
    geometry-default dispatch in :func:`solve_sn_fixed_source` clean.
    """
    from orpheus.transport.sources import PerOrdinateSource
    from .angular_flux import AngularFlux
    from orpheus.transport.fields.scalar_flux import ScalarFlux
    nx, ny, ng = sn_mesh.nx, sn_mesh.ny, solver.ng

    # Issue #196 PR-INDEX-5: phi principled (ng, nx, ny).
    phi = np.zeros((ng, nx, ny))
    angular = None
    residual = np.inf
    converged_flag = False
    flux_residuals: list[float] = []

    for n_inner in range(max_inner):
        phi_prev = phi

        # Isotropic in-scatter + (n,2n). No fission — this is fixed-source.
        # R-1 Step 4 A1 — single per-ordinate source feeds the sweep.
        # Producer-side normalisation: ``Q_iso`` (iso scalar magnitude)
        # projects to per-ord via :meth:`PerOrdinateSource.from_isotropic`;
        # ``Q_aniso_p1`` is per-ord density already
        # (:meth:`ScatteringOperator.build_aniso_source` applies the
        # ``/W`` at the producer post-A1); ``external_source`` is
        # per-ord density by API contract (the caller of
        # :func:`solve_sn_fixed_source` projects iso sources via
        # :meth:`PerOrdinateSource.from_isotropic` before passing).
        Q_iso = np.zeros_like(phi)
        solver._add_scattering_source(Q_iso, phi)
        solver._add_n2n_source(Q_iso, phi)
        iso_per_ord = PerOrdinateSource.from_isotropic(Q_iso, sn_mesh)

        # Merge per-ord pieces — Pℓ scattering moments + external source.
        Q_aniso_p1 = solver._build_aniso_scattering(angular)
        total_values = iso_per_ord.values + external_source
        if Q_aniso_p1 is not None:
            total_values = total_values + Q_aniso_p1
        source = PerOrdinateSource.from_mesh(total_values, sn_mesh)

        # R-1 Step 0: thread previous-iter angular flux as initial_guess
        # for the trace-space Carlson seed.
        initial_guess = (
            AngularFlux(angular, sn_mesh) if angular is not None else None
        )
        angular, phi = transport_sweep(
            source, solver.mat_xs.total_cross_section, sn_mesh,
            solver._boundary_flux,
            initial_guess=initial_guess,
        )

        norm = np.linalg.norm(phi)
        if norm > 0:
            residual = np.linalg.norm(phi - phi_prev) / norm
            flux_residuals.append(float(residual))
            if residual < inner_tol:
                converged_flag = True
                break
    else:
        n_inner = max_inner - 1  # loop exhausted without break

    # Issue #197 PR-TYPED-5: build typed Solution at the boundary.
    # mesh / quadrature / materials remain unconsumed by Solution — the
    # typed fluxes carry the SNMesh reference, which transitively exposes
    # those handles via .mesh.{mesh, quad, materials}.
    del mesh, quadrature, materials  # retained as kwargs for API stability
    history = IterationHistory(
        flux_residuals=tuple(flux_residuals),
        n_inner=n_inner + 1,
        converged=converged_flag,
    )
    # D-H.1b (2026-05-28): Solution.angular_flux is now a TimedFullField.
    from orpheus.transport.timed_full_field import TimedFullField
    legacy_psi = AngularFlux(angular, sn_mesh, boundary=solver._boundary_flux)
    return Solution(
        angular_flux=TimedFullField.from_legacy_angular_flux(legacy_psi),
        scalar_flux=ScalarFlux.from_mesh(phi, sn_mesh),
        mesh=sn_mesh,
        keff=None,
        history=history,
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
) -> Solution:
    r"""Curvilinear-default fixed-source path: typed :class:`KrylovAcceleration`.

    R-1 Step 4 Step G1 (2026-05-22) — carved onto
    :class:`~orpheus.numerics.iteration.KrylovAcceleration` consuming
    the operator triple

    .. math::

        A \;=\; L + C\,, \quad
        S \;=\; \text{full multi-group scatter}\,, \quad
        F \;=\; 0_{\rm wg}

    on typed :class:`~orpheus.sn.angular_flux.AngularFlux`.  The
    composite ``L + C`` returns an
    :class:`~orpheus.sn.operator.InvertibleOperator`; its ``.solve``
    IS the WDD sweep, but R-1 ships GMRES UNPRECONDITIONED
    (explicit identity) — issue #200 tracks the block-inverse face
    preconditioner re-enablement.

    The outer Picard wrap on the scattering source dissolves at G1:
    typed ``KrylovAcceleration`` solves :math:`(L+C-S)\,\psi =
    q_{\rm ext}` directly via GMRES, putting scattering
    self-consistency INSIDE the Krylov polynomial.  Fission ``F`` is
    zero — this is the pure fixed-source path; the eigenvalue
    outer/within-group decomposition lives in :func:`solve_sn`.

    The dependency audit + retirement sequence eliminates the
    packed-1D vector machinery this function previously consumed:

    * ``build_equation_map_*`` (retires G3e)
    * ``solution_to_angular_flux_*`` (retires G3e)
    * ``_build_rhs_{cartesian, spherical, cylindrical}`` —
      **RETIRED** in P1.7 of the moment-space + layering plan (the
      Cartesian variant carried an inline ``(2*l+1)`` literal that
      duplicated
      :attr:`~orpheus.numerics.spaces.SphericalHarmonicSpace.addition_theorem_factor`;
      all three were already dead code in production).
    * ``_make_sweep_preconditioner`` (retires G3a)
    * ``solver.L`` / ``SNStreamingOperator`` (retires G3f)

    Scope
    =====

    1-D meshes only (slab + sphere + cylinder).  2-D Cartesian raises
    :class:`NotImplementedError` — the typed-flux B1'' face block was
    designed 1-D-only.  SURPRISE-5 of the dependency audit identifies
    :func:`_solve_fixed_source_si` as the principled landing zone for
    2-D Cartesian fixed-source: it's geometry-agnostic via
    :func:`transport_sweep` and handles 2-D natively.  Phase A absorbs
    a typed 2-D Krylov.
    """
    # 2-D Cartesian deferral — mirrors :meth:`SNSolver._solve_krylov`
    # (the eigenvalue inner): 2-D Cartesian fixed-source Krylov has no
    # typed-flux B1'' face block in 1-D-only shape.  SURPRISE-5 of the
    # dependency audit identifies SI as the principled landing zone
    # (geometry-agnostic via :func:`transport_sweep`).  Phase A absorbs
    # typed 2-D Krylov.
    if sn_mesh.reduced is None:
        raise NotImplementedError(
            "R-1 Step 4 G1 — 2-D Cartesian fixed-source Krylov carve "
            "deferred to Phase A.  The typed B1'' face block is "
            "1-D-only.  Use `inner_solver=\"source_iteration\"` for 2-D "
            "Cartesian (the SI carve in `_solve_fixed_source_si` is "
            "geometry-agnostic and handles 2-D via "
            "`transport_sweep`'s wavefront dispatch — see SURPRISE-5 "
            "of `.claude/plans/r1_step4_g_dependency_audit.md`)."
        )

    from .angular_flux import AngularFlux
    from orpheus.transport.fields.scalar_flux import ScalarFlux
    from .operator import CollisionOperator, StreamingOperator
    from orpheus.transport.sources import PerOrdinateSource
    from orpheus.numerics.iteration import KrylovAcceleration
    from orpheus.numerics.operator import ZeroOperator

    nx, ny, ng = sn_mesh.nx, sn_mesh.ny, solver.ng
    N = sn_mesh.quad.N

    # ── Build the typed RHS ──────────────────────────────────────────
    # R-1 Step 4 A1 — ``external_source`` is per-ordinate density (the
    # producer-side ``/sum_w`` projection lives at
    # :meth:`PerOrdinateSource.from_isotropic` for iso scalar sources;
    # by the time we get here the caller has projected).  Wrap into
    # a typed :class:`AngularFlux` for the path-forward Krylov.
    q_ext_per_ord = PerOrdinateSource.from_mesh(external_source, sn_mesh)
    q_ext_typed = AngularFlux(q_ext_per_ord.values, sn_mesh)

    # ── Build the typed operator triple ─────────────────────────────
    # InvertibleOperator via __add__ dispatch (R-1 Step C).  L + C
    # admits .apply (matvec) and .solve (sweep) capabilities.  Scattering
    # is the full multi-group operator; F = 0 for fixed-source.
    L_leaf = StreamingOperator(sn_mesh, solver.mat_xs.total_cross_section)
    C_t = CollisionOperator(sn_mesh, solver.mat_xs.total_cross_section)
    LC = L_leaf + C_t

    # GMRES UNPRECONDITIONED (explicit identity) — issue #200 tracks
    # the block-inverse face preconditioner re-enablement.  Post-Phase-1.2
    # the silent-fallback bug class is structurally unreachable:
    # ``InvertibleOperator.solve`` is a pure function of
    # ``(rhs, initial_guess, boundary)``; ``KrylovAcceleration``'s
    # ``preconditioner=None`` default would invoke ``L.solve(q)`` with
    # ``initial_guess=None`` → cold-start M-M seed (deterministic).
    # Explicit identity is kept until #200 lands a curvilinear-quality
    # preconditioner.
    krylov = KrylovAcceleration(
        LC,
        solver.scattering_op,
        ZeroOperator(),
        preconditioner=lambda q: q,
        tol=inner_tol,
        max_iter=max_inner,
        restart=min(50, N * ng * nx * ny),
    )

    psi_typed, residuals = krylov.solve(q_ext_typed)
    angular = psi_typed.values
    phi = psi_typed.integrate_angular().values
    converged_flag = bool(residuals) and residuals[-1] < inner_tol
    n_outer = len(residuals)
    flux_residuals = [float(r) for r in residuals]

    # Issue #197 PR-TYPED-5: build typed Solution at the boundary.
    # R-1 Step 4 G1 — ``psi_typed`` is already a typed AngularFlux with
    # the matvec's B1'' face residual on its boundary; reuse directly.
    del mesh, quadrature, materials  # retained as kwargs for API stability
    history = IterationHistory(
        flux_residuals=tuple(flux_residuals),
        n_inner=n_outer + 1,
        converged=converged_flag,
    )
    # D-H.1b (2026-05-28): Solution.angular_flux is now a TimedFullField.
    # psi_typed is the Krylov adapter output — still a legacy AngularFlux
    # (the Krylov adapter unravels back to the template type, which is the
    # legacy q_ext_typed). Wrap into TimedFullField at the Solution
    # boundary via the legacy adapter.
    from orpheus.transport.timed_full_field import TimedFullField
    return Solution(
        angular_flux=TimedFullField.from_legacy_angular_flux(psi_typed),
        scalar_flux=ScalarFlux.from_mesh(phi, sn_mesh),
        mesh=sn_mesh,
        keff=None,
        history=history,
    )
