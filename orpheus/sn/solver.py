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
  by the algebraic composition
  :class:`~orpheus.sn.operator.InvertibleOperator` (= ``L + C``).
  This is the Wave E reconciliation that makes the curvilinear
  ``solve_sn_fixed_source`` path discretization-correct (closes
  ERR-026).  On Cartesian meshes it is bit-identical math to the
  legacy BiCGSTAB FD path.

Boundary conditions default to reflective (infinite lattice) but are
configurable via :class:`~orpheus.geometry.mesh.BC` on the mesh.

.. seealso:: :ref:`theory-discrete-ordinates` — Key Facts, equations, gotchas.
"""

from __future__ import annotations

import time
from dataclasses import replace

import numpy as np
from scipy.sparse.linalg import gmres

from orpheus.data.macro_xs.cell_xs import assemble_cell_xs
from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import BC, Mesh1D, Mesh2D
from orpheus.numerics.eigenvalue import power_iteration
from .fission import FissionOperator
from .geometry import SNMesh
from .spatial.sweep_cache import CollisionCache, GeometryCoefficients
from .operator import (
    CollisionOperator,
    StreamingOperator,
)
from orpheus.numerics.quadrature import Quadrature
from .scattering import ScatteringOperator
from .sweep import transport_sweep
from orpheus.transport.fields.boundary_flux import BoundaryFlux


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


def _zero_within_group_fission(psi: "object") -> "object":
    r"""Codomain zero for the within-group zero-fission slot (``F = 0``).

    A fission operator maps FLUX → SOURCE, so its zero action emits a *zero
    source* on BOTH blocks — an ``AngularSourceSink`` bulk AND a
    ``BoundarySourceSink`` boundary — NOT a flux echo of the iterate ``ψ`` (the
    B.5.2 ``ZeroOperator`` codomain fix: a zero operator returns the zero of its
    CODOMAIN, not of its domain).  Wired into ``ZeroOperator(codomain_zero=...)``
    so the typed RHS ``S.apply(ψ) + F.apply(ψ) + q_ext`` and the Krylov matvec
    ``L.apply − S.apply − F.apply`` stay CLOSED source/sink sums on BOTH blocks
    (every operator output is a source/sink — its boundary is
    ``BoundarySourceSink``, mirroring the bulk; the residual only arises from
    ``from_balance``).  The within-group fission source proper enters as
    ``q_ext`` (the eigenvalue outer's contribution); ``F`` itself is
    structurally zero here.
    """
    from orpheus.transport.source_sinks import (
        AngularSourceSink,
        BoundarySourceSink,
    )
    from orpheus.transport.timed_full_field import TimedFullField

    return TimedFullField.zeros(
        bulk=AngularSourceSink, boundary=BoundarySourceSink, mesh=psi.bulk.mesh,
    )


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
    * ``"krylov"`` — GMRES on
      :class:`~orpheus.sn.operator.InvertibleOperator` (the algebraic
      composition ``L + C``; symmetric closure) with the sweep as
      preconditioner.  Closes ERR-026 on curvilinear; bit-identical
      math to the legacy BiCGSTAB FD path on Cartesian.

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
                f"InvertibleOperator (L + C) with the sweep as "
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
        self._boundary_flux = BoundaryFlux.zeros_on(sn_mesh)

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
        # The full transport operator :math:`L = \Omega\cdot\nabla + \Sigma_t`
        # as the algebraic sum :class:`StreamingOperator` +
        # :class:`CollisionOperator` = :class:`InvertibleOperator`.
        # :meth:`InvertibleOperator.apply` returns ``(L_stream + C)·ψ`` in
        # :class:`TimedFullField` form (the typed direct-sum carrier);
        # :meth:`InvertibleOperator.solve` consumes ``initial_guess`` for
        # the Carlson seed (R-1 Phase 1.2 unification).
        self.L = (
            StreamingOperator(sn_mesh, self.mat_xs.total_cross_section)
            + CollisionOperator(sn_mesh, self.mat_xs.total_cross_section)
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
        self.L = (
            StreamingOperator(self.sn_mesh, self.mat_xs.total_cross_section)
            + CollisionOperator(self.sn_mesh, self.mat_xs.total_cross_section)
        )
        if self.geom_cache is not None:
            sig_t_1d = self.mat_xs.total_cross_section[:, :, 0]
            self.coll_cache = CollisionCache.from_geometry(
                self.geom_cache, sig_t_1d,
            )
            self.sn_mesh._coll_cache = self.coll_cache  # type: ignore[attr-defined]

    @property
    def _scattering_with_boundary_op(self):
        r"""The scattering operator folded with the realized boundary law ``B``.

        Wave O (#208) BC-extraction: the canonical SN loss is
        ``(L + C − S − F − B)`` where ``B`` (the realized reflective / albedo /
        white law) is a first-class sibling of ``L``. ``B`` CANNOT join the
        preconditioner ``L + C`` — :class:`~orpheus.numerics.operator.OperatorSum`
        does not propagate ``CAP_SOLVE``, so ``L + C − B`` would strip the
        ``.solve`` (sweep) the Krylov preconditioner / SI splitting needs. ``B``
        therefore folds into the **subtracted** ``S`` argument of the iteration
        drivers, whose matvec is ``L.apply − S.apply − F.apply``: passing
        ``S + B`` makes that ``(L + C).apply − (S + B).apply − F.apply
        = (L + C − S − F − B)``. The reflective inflow trace is then driven by
        the boundary consistency residual ``ψ.inflow − B·ψ.outflow`` instead of
        the intra-sweep ``bc.apply`` re-application (the deleted keystone).

        The :class:`OperatorSum` domain-compatibility check skips because
        ``ScatteringOperator.domain`` is ``None`` (it predates function-space
        tagging) — the check fires only when BOTH operands declare non-``None``
        domains that differ, so it is symmetric in the operands and the sum
        order is irrelevant here (``B + S`` would skip identically).
        **O.2 forcing function / tripwire:** the moment ``ScatteringOperator``
        gains a ``domain`` (the typing wave will give it one), this fold throws
        ``IncompatibleOperatorComposition`` at construction — ``S`` (bulk space)
        and ``B`` (trace space) live on different function spaces. That failure
        IS the signal to land the honest composition (drivers consuming the
        whole loss operator ``L+C−S−F−B`` on the direct-sum carrier) — Wave O
        step O.2. Single source of truth for the ``S + B`` fold (Cardinal
        Rule 2): every Krylov site + the eigenvalue SI driver reads this.

        TWIN ROUTE (same ``−B`` coupling, different plumbing): the SI/Krylov
        DRIVER path delivers ``B·ψ.outflow`` through ``rhs.boundary`` via this
        fold; the DIRECT fixed-source loops + the final eigenvalue
        reconstruction sweep — which have no driver to fold into — deliver it
        via :func:`_reflect_outflow_into_inflow` instead. Both call the
        identical :class:`~orpheus.sn.boundary_operator.SNBoundaryOperator`
        (``B`` is single-sourced); only the seed-delivery plumbing differs. The
        two routes COLLAPSE at O.2 (drivers take ``L+C−S−F−B`` directly, the
        direct loops route through the driver, the helper retires).
        """
        from orpheus.sn.boundary_operator import SNBoundaryOperator
        return self.scattering_op + SNBoundaryOperator(self.sn_mesh)

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

        ALL geometries — slab, sphere, cylinder, AND 2-D Cartesian
        (Wave O "2-D SI Phase A", 2026-06-04).  The 2-D Cartesian
        eigenvalue SI inner is geometry-agnostic: it is the structural
        twin of :meth:`_solve_krylov` (the live 2-D eigenvalue path) —
        identical composite RHS, identical operator triple
        (``LC = StreamingOperator + CollisionOperator``,
        ``self._scattering_with_boundary_op`` for the ``S + B`` fold,
        zero within-group fission), identical
        ``psi_typed.bulk.integrate_angular()`` reduction — differing
        ONLY in the driver (:class:`SourceIteration` vs
        :class:`KrylovAcceleration`), and neither driver carries any
        geometry dependence.  The reflective coupling rides the BARE
        sweep via the ``S + B`` fold on the 4-face
        :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`
        (:class:`SNBoundaryOperator` is natively 4-face —
        xmin/xmax/ymin/ymax — and is the SAME operator the working 2-D
        Krylov path uses).

        The legacy "B1'' face block" that the 2-D path was once said to
        lack is RETIRED — it never existed as code on this branch; it
        was a 1-D boundary-closure name fully superseded by the L2
        ``BoundaryFlux`` + ``SNBoundaryOperator`` bare-boundary
        architecture (Wave O O.4a.2 / O.4b).  The historical
        ``NotImplementedError`` guard predated that architecture by two
        weeks (the body's ``S + B`` fold landed 2026-06-03) and was
        never lifted; ``solve_sn`` defaults ``inner_solver="source_-
        iteration"`` for every geometry, so the stale guard broke the
        DEFAULT 2-D Cartesian eigenvalue entry point.

        Verified: 2-D SI ≡ Krylov ≡ closed-form ``k_inf`` (1g→1.5,
        2g→1.875, 4g→1.4878) to machine precision; heterogeneous
        non-flat 2-D flux SHAPE agrees SI-vs-Krylov to ~1e-9
        (``tests/sn/eigenvalue/test_keff_2d.py::TestSIKrylov2DEquivalence``).
        """
        from .operator import CollisionOperator, StreamingOperator
        from orpheus.transport.fields.angular_flux import (
            AngularFlux,
        )
        from orpheus.transport.fields.boundary_flux import (
            BoundaryFlux,
        )
        from orpheus.transport.source_sinks import (
            AngularSourceSink,
            BoundarySourceSink,
        )
        from orpheus.transport.timed_full_field import TimedFullField
        from orpheus.numerics.iteration import SourceIteration
        from orpheus.numerics.operator import ZeroOperator

        # ── Build the composite RHS ─────────────────────────────────
        # R-1 Step 4 A1 — q_ext as per-ordinate density via the canonical
        # ``AngularSourceSink.from_isotropic`` factory.  The /W projection
        # lives at the factory boundary (Pattern 7 producer-side
        # normalisation); the legacy
        # ``q_per_ord = (fission_source / sum_w)[None, :, :, :]``
        # broadcast pattern is GONE.
        q_ext_per_ord = AngularSourceSink.from_isotropic(
            fission_source, self.sn_mesh,
        )
        # D-H.1c stage 2 — q_ext composite carries:
        #   * bulk = per-ordinate source values on the L2 AngularFlux.
        #   * boundary = ZERO (the EXTERNAL boundary source — zero for
        #     vacuum/reflective; non-zero only for prescribed inflow).
        # Wave O (#208) O.4a.2 — the partner-flux seeding
        # (``boundary = self._boundary_flux``) is RETIRED: the reflective
        # inflow is no longer pre-staged into the source.  It is driven by
        # the sibling ``−B`` folded into the subtracted ``S`` argument
        # (``self._scattering_with_boundary_op``): each SI iterate adds
        # ``B·ψ.outflow`` to ``rhs.boundary``, which ``(L+C).solve``'s bare
        # sweep reads as the inflow seed (operator.py
        # ``_solve_timed_full_field`` seeds from ``rhs.boundary``).  The
        # boundary inflow is thus a live solved unknown carried in
        # ``ψ.boundary``, not an externally-recomputed partner trace.
        q_ext_composite = TimedFullField(
            # B.5.2: q_ext IS a source — carry the AngularSourceSink bulk AND
            # the BoundarySourceSink inflow trace (zero for vacuum/reflective;
            # prescribed inflow otherwise). The SI rhs F.apply + (S+B).apply +
            # q_ext closes on BoundarySourceSink (operator outputs are sources).
            bulk=q_ext_per_ord,
            boundary=BoundarySourceSink.zeros_on(self.sn_mesh),
            _history=(),
            history_depth=2,
        )

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
            # Wave O #208 — ``S + B`` fold: the SI source becomes
            # ``(S+B)ψ + Fψ + q`` so the reflective inflow ``B·ψ.outflow``
            # rides in ``rhs.boundary`` and the bare ``(L+C).solve`` sweep
            # reads it as the inflow seed. (B cannot join the LC
            # preconditioner — OperatorSum drops CAP_SOLVE.)
            self._scattering_with_boundary_op,
            ZeroOperator(codomain_zero=_zero_within_group_fission),
            max_iter=self.max_inner,
            tol=self.inner_tol,
        )

        # ── Warm start (composite) ──────────────────────────────────
        # SourceIteration._solve_with_seed forwards the previous
        # iterate to InvertibleOperator.solve via the explicit
        # ``initial_guess`` kwarg (Phase 1.2 — the M-M closure's
        # ``psi_half_seed`` strategy reads it to derive the Carlson
        # coupled-pole seed; the previous iterate's boundary trace
        # seeds the reflective-BC partner-flux state).  Post-D-H.1c
        # stage 2 ``self._psi_typed`` is a :class:`TimedFullField`;
        # the type propagates through the iteration primitive via the
        # ravellable protocol.
        initial_guess = getattr(self, "_psi_typed", None)
        if initial_guess is None:
            # B.5.2: cold-start iterate is a FLUX composite, decoupled from
            # q_ext's now-AngularSourceSink type.  Bit-identical to the prior
            # cold start (_zeros_like(q_ext): all-zeros, BoundaryFlux boundary)
            # — only the bulk CLASS flips source→flux (the iterate's true role).
            initial_guess = TimedFullField.zeros(
                bulk=AngularFlux, boundary=BoundaryFlux, mesh=self.sn_mesh,
            )

        psi_typed, _residuals = si.solve(
            q_ext_composite, initial_guess=initial_guess,
        )
        self._psi_typed = psi_typed

        # Reduce angular → scalar flux for the eigenvalue outer's contract.
        return psi_typed.bulk.integrate_angular().values

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

        D-H.2-C4d (2026-05-29) — 2-D Cartesian unblocked.  The L2
        :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`
        face layout is the natural 4-face descriptor (xmin / xmax /
        ymin / ymax) that the legacy 1-D B1'' face block lacked; the
        L2-native ``StreamingOperator._apply_2d_cartesian`` kernel
        operates on it directly.

        Returns the updated scalar flux ``(ng, nx, ny)``.
        """

        from .operator import CollisionOperator, StreamingOperator
        from orpheus.transport.fields.angular_flux import (
            AngularFlux,
        )
        from orpheus.transport.fields.boundary_flux import (
            BoundaryFlux,
        )
        from orpheus.transport.source_sinks import (
            AngularSourceSink,
            BoundarySourceSink,
        )
        from orpheus.transport.timed_full_field import TimedFullField
        from orpheus.numerics.iteration import KrylovAcceleration
        from orpheus.numerics.operator import ZeroOperator

        # ── Build the composite RHS ─────────────────────────────────
        # R-1 Step 4 A1 — q_ext as per-ordinate density via the canonical
        # ``AngularSourceSink.from_isotropic`` factory.  /W lives at the
        # factory boundary (Pattern 7 producer-side normalisation).
        # D-H.1c stage 2 — TimedFullField bulk + zero boundary (the
        # Krylov path does NOT pre-seed q_ext.boundary; reflective-BC
        # state threads through ``initial_guess`` per the audit §5
        # contract).
        q_ext_per_ord = AngularSourceSink.from_isotropic(
            fission_source, self.sn_mesh,
        )
        q_ext_composite = TimedFullField(
            # B.5.2: q_ext IS a source — carry the AngularSourceSink bulk AND
            # the BoundarySourceSink inflow trace (zero for vacuum/reflective;
            # prescribed inflow otherwise). The SI rhs F.apply + (S+B).apply +
            # q_ext closes on BoundarySourceSink (operator outputs are sources).
            bulk=q_ext_per_ord,
            boundary=BoundarySourceSink.zeros_on(self.sn_mesh),
            _history=(),
            history_depth=2,
        )

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
            # Wave O #208 — ``S + B`` fold: the matvec ``(L+C).apply −
            # (S+B).apply − F.apply`` IS ``(L+C−S−F−B)``, with the realized
            # boundary law ``B`` as a sibling of ``L`` (BC extraction).
            self._scattering_with_boundary_op,
            ZeroOperator(codomain_zero=_zero_within_group_fission),
            preconditioner=lambda q: q,  # explicit identity — issue #200 tracks re-enable
            tol=self.inner_tol,
            max_iter=self.max_inner,
            # D-H.1e (2026-05-28): restart sized to the full problem
            # (do NOT clamp at 50).  Pre-D-H.1e clamp truncated the
            # Krylov subspace below the natural domain size, leaving
            # GMRES structurally unconverged on curvilinear meshes with
            # ``N*ng*nx*ny > 50``.  See ERR-053.
            restart=N * ng * nx * ny,
        )

        # ── Warm start (composite) ──────────────────────────────────
        # Post-D-H.1c stage 2: ``self._psi_typed`` is a TimedFullField;
        # the Krylov ravellable protocol detects the composite via
        # ``to_flat`` / ``from_flat`` (D-H.1b.1) and threads it through
        # the matvec / unravel cycle natively.
        initial_guess = getattr(self, "_psi_typed", None)
        if initial_guess is None:
            # B.5.2: cold-start iterate is a FLUX composite, decoupled from
            # q_ext's now-AngularSourceSink type.  x0 stays all-zeros
            # (bit-identical); the flux template fixes the Krylov return type.
            initial_guess = TimedFullField.zeros(
                bulk=AngularFlux, boundary=BoundaryFlux, mesh=self.sn_mesh,
            )

        psi_typed, _residuals = krylov.solve(
            q_ext_composite, initial_guess=initial_guess,
        )
        self._psi_typed = psi_typed

        # Reduce angular → scalar flux for the eigenvalue outer's contract.
        return psi_typed.bulk.integrate_angular().values

    # D-J (2026-05-29): ``_make_sweep_preconditioner`` retired —
    # the method built a scipy LinearOperator wrapping the sweep as
    # the legacy packed-vector GMRES preconditioner.  Zero callers
    # remained post-D-K.1: the production Krylov path
    # (``_solve_fixed_source_krylov``) wraps the L+C composite in a
    # KrylovAcceleration consumer that handles the to_flat/from_flat
    # ravellable protocol on :class:`TimedFullField` natively.

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
        ``None`` when ``L == 0``).

        D-I.2 (2026-05-29): the underlying
        :meth:`ScatteringOperator.build_aniso_source` retired its
        bare-ndarray arm and now accepts only typed
        :class:`AngularFlux`.  This delegator wraps the inbound
        ``angular_flux`` ndarray as :class:`AngularFlux` at the helper
        boundary, then unwraps the resulting
        :class:`AngularSourceSink` via ``.values`` so the bare-ndarray
        external contract is preserved for legacy consumers in
        :func:`_solve_fixed_source_si` and the verification probes in
        :mod:`tests.sn.test_scattering_operator`.
        """
        if angular_flux is None:
            return None
        from orpheus.transport.fields.angular_flux import AngularFlux
        typed = AngularFlux.from_mesh(angular_flux, self.sn_mesh)
        result = self.scattering_op.build_aniso_source(typed)
        if result is None:
            return None
        return result.values

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
    typed L2 ``AngularFlux`` boundary; closes the PR-INDEX-4 §9.1
    deferral).  The legacy ``solution_to_angular_flux*`` decoder
    family retired at D-J (2026-05-30).

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


def _reflect_outflow_into_inflow(boundary_flux, sn_mesh: SNMesh) -> None:
    r"""In-place: fill each face's inflow ordinate slots with the realized
    boundary law applied to that face's outflow trace — the ``−B`` reflective
    coupling, externalised for the bare ``transport_sweep`` (Wave O #208 O.4a.2).

    The bare sweep reads the inflow ordinate slots of its boundary buffer as
    the inflow seed; it no longer re-applies ``bc`` to the outflow internally.
    The SI driver path supplies ``B·ψ.outflow`` through ``rhs.boundary`` (the
    ``S + B`` fold), but the DIRECT fixed-source SI loop
    (:func:`_solve_fixed_source_si`) and the final eigenvalue reconstruction
    sweep (:func:`solve_sn`) do not route through that driver — they call this
    helper to set ``ψ.inflow = B·ψ.outflow`` on the buffer before each sweep,
    via the canonical whole-trace :class:`~orpheus.sn.boundary_operator.SNBoundaryOperator`
    (single source of truth — the same ``B`` the matvec / SI driver consume).

    For vacuum ``B = 0`` so the inflow slots stay zero (bit-identical to the
    pre-extraction ``bc.apply`` of a vacuum law); for reflective/white/albedo
    it is the same ``R·G`` reflection the pre-extraction sweep applied at entry,
    merely relocated to the caller.

    This is the TWIN of the driver route :meth:`SNSolver._scattering_with_boundary_op`
    (the ``S + B`` fold): both deliver the SAME ``−B`` coupling via the SAME
    ``SNBoundaryOperator``, differing only in plumbing (this helper writes the
    buffer's inflow slots directly; the fold rides ``B·ψ.outflow`` in
    ``rhs.boundary``). The two COLLAPSE at Wave O step O.2 — when the iteration
    drivers take the whole loss operator ``L+C−S−F−B`` directly, the direct
    loops route through the driver and THIS HELPER RETIRES (the removal trigger).
    """
    from orpheus.sn.boundary_operator import SNBoundaryOperator
    from orpheus.transport.fields.angular_flux import AngularFlux
    from orpheus.transport.timed_full_field import TimedFullField

    probe = TimedFullField(
        bulk=AngularFlux.from_mesh(
            np.zeros(
                (sn_mesh.quad.N, sn_mesh.ng, sn_mesh.nx, sn_mesh.ny),
            ),
            sn_mesh,
        ),
        boundary=boundary_flux,
        _history=(),
        history_depth=2,
    )
    reflected = SNBoundaryOperator(sn_mesh).apply(probe).boundary
    trace = sn_mesh.trace
    for face in boundary_flux.layout.faces:
        inflow = trace.inflow_indices_for_face(face)
        boundary_flux.face_view(face)[inflow] = reflected.face_view(face)[inflow]


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
    # :meth:`AngularSourceSink.from_isotropic` factory (Pattern 7
    # producer-side normalisation — /W projection at the factory
    # boundary).
    from orpheus.transport.source_sinks import AngularSourceSink
    from orpheus.transport.fields.boundary_flux import (
        BoundaryFlux as _BoundaryFlux,
    )
    Q_final = solver.compute_fission_source(scalar_flux, keff)
    solver._add_scattering_source(Q_final, scalar_flux)
    solver._add_n2n_source(Q_final, scalar_flux)
    # Wave O #208 O.4a.2 — BARE final reconstruction sweep: seed its inflow
    # from the CONVERGED boundary trace.  The inner solves drive the inflow as
    # a live unknown in ``ψ.boundary`` (carried on ``solver._psi_typed`` from
    # the last inner solve), so ``solver._boundary_flux`` is no longer the
    # partner-flux carrier (it stays all-zeros).  Reflect the converged outflow
    # into the inflow slots via −B (no-op for vacuum; idempotent here since the
    # converged inflow already equals ``B·ψ.outflow``), then sweep.
    converged = getattr(solver, "_psi_typed", None)
    final_boundary = (
        converged.boundary if converged is not None
        else _BoundaryFlux.zeros_on(sn_mesh)
    )
    # Wave O #208 O.4b Phase E2 — the 2-D wavefront is now BARE (reads the
    # given inflow, no in-sweep bc.apply), so the reflective coupling is the
    # external -B for 2-D too.  The guard is lifted: _reflect_outflow_into_inflow
    # is geometry-agnostic (iterates boundary_flux.layout.faces via the canonical
    # SNBoundaryOperator — verified 2-D-ready) and idempotent here (the converged
    # inflow already equals B·ψ.outflow); vacuum stays a no-op (B = 0).
    _reflect_outflow_into_inflow(final_boundary, sn_mesh)
    angular_flux, _ = transport_sweep(
        AngularSourceSink.from_isotropic(Q_final, sn_mesh),
        solver.mat_xs.total_cross_section, sn_mesh,
        final_boundary,
    )

    elapsed = time.perf_counter() - t_start

    # Issue #197 PR-TYPED-5: build typed Solution.  The bare ndarrays
    # produced by power_iteration + the final sweep get wrapped in
    # typed fields at this single boundary; downstream consumers read
    # typed fields.
    # D-H.1c stage 2 (2026-05-28): Solution.angular_flux is a
    # TimedFullField composite (bulk + boundary + history) constructed
    # DIRECTLY from L2 types.  The legacy-AngularFlux adapter wrap is
    # GONE — ``angular_flux`` is a bare ndarray from the final sweep;
    # we wrap it once into the composite carrier.
    from orpheus.transport.fields.angular_flux import (
        AngularFlux,
    )
    from orpheus.transport.fields.boundary_flux import (
        BoundaryFlux,
    )
    from orpheus.transport.fields.scalar_flux import ScalarFlux
    from orpheus.transport.timed_full_field import TimedFullField
    history = IterationHistory(
        keff_history=tuple(keff_history),
        n_outer=len(keff_history),
        converged=True,
    )
    return Solution(
        angular_flux=TimedFullField(
            bulk=AngularFlux.from_mesh(angular_flux, sn_mesh),
            # Wave O #208 O.4a.2: the converged boundary trace from the final
            # bare sweep (inflow = B·ψ.outflow, outflow = streamed).
            boundary=final_boundary,
            _history=(),
            history_depth=2,
        ),
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
        :meth:`~orpheus.sn.sources.AngularSourceSink.from_isotropic`
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
        path.  ``"krylov"`` is available as an opt-in for vacuum /
        reflective / white / albedo / mixed
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
    from orpheus.transport.source_sinks import AngularSourceSink
    from orpheus.transport.fields.angular_flux import (
        AngularFlux,
    )
    from orpheus.transport.fields.boundary_flux import (
        BoundaryFlux,
    )
    from orpheus.transport.fields.scalar_flux import ScalarFlux
    from orpheus.transport.timed_full_field import TimedFullField
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
        # projects to per-ord via :meth:`AngularSourceSink.from_isotropic`;
        # ``Q_aniso_p1`` is per-ord density already
        # (:meth:`ScatteringOperator.build_aniso_source` applies the
        # ``/W`` at the producer post-A1); ``external_source`` is
        # per-ord density by API contract (the caller of
        # :func:`solve_sn_fixed_source` projects iso sources via
        # :meth:`AngularSourceSink.from_isotropic` before passing).
        Q_iso = np.zeros_like(phi)
        solver._add_scattering_source(Q_iso, phi)
        solver._add_n2n_source(Q_iso, phi)
        iso_per_ord = AngularSourceSink.from_isotropic(Q_iso, sn_mesh)

        # Merge per-ord pieces — Pℓ scattering moments + external source.
        Q_aniso_p1 = solver._build_aniso_scattering(angular)
        total_values = iso_per_ord.values + external_source
        if Q_aniso_p1 is not None:
            total_values = total_values + Q_aniso_p1
        source = AngularSourceSink.from_mesh(total_values, sn_mesh)

        # R-1 Step 0: thread previous-iter angular flux as initial_guess
        # for the trace-space Carlson seed.  D-H.2-C5 phase 2: composite
        # TimedFullField directly (legacy AngularFlux retired); the
        # sweep's ``_initial_guess_values`` duck-types on ``.bulk``.
        if angular is not None:
            initial_guess = TimedFullField(
                bulk=AngularFlux.from_mesh(angular, sn_mesh),
                boundary=BoundaryFlux.zeros_on(sn_mesh),
                _history=(),
                history_depth=2,
            )
        else:
            initial_guess = None
        # Wave O #208 O.4a.2 / O.4b E2 — BARE sweep (1-D AND 2-D): reflect the
        # persisted outflow (``solver._boundary_flux`` outflow slots, from the
        # previous iteration's sweep) into the inflow slots via −B BEFORE the
        # sweep.  The bare sweep reads the inflow slots directly; it no longer
        # re-applies ``bc`` at entry.  (First iteration: outflow = 0 ⟹
        # inflow = B·0 = 0, matching the pre-extraction cold start.)
        # O.4b E1 made the 2-D Cartesian wavefront sweep bare too, so the guard
        # is lifted — ``_reflect_outflow_into_inflow`` is geometry-agnostic
        # (canonical ``SNBoundaryOperator``, verified 2-D-ready) and vacuum
        # stays a no-op (B = 0).
        _reflect_outflow_into_inflow(solver._boundary_flux, sn_mesh)
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
    # D-H.1c stage 2 (2026-05-28): Solution.angular_flux constructed
    # DIRECTLY as a TimedFullField composite.  ``angular`` is the bare
    # ndarray returned by transport_sweep; ``solver._boundary_flux`` is
    # the legacy BoundaryFlux carrying the converged partner-flux
    # trace.  No legacy AngularFlux intermediate.
    from orpheus.transport.fields.angular_flux import (
        AngularFlux,
    )
    from orpheus.transport.fields.boundary_flux import (
        BoundaryFlux,
    )
    from orpheus.transport.timed_full_field import TimedFullField
    return Solution(
        angular_flux=TimedFullField(
            bulk=AngularFlux.from_mesh(angular, sn_mesh),
            # D-H.2-C2: ``solver._boundary_flux`` is L2 directly.
            boundary=solver._boundary_flux,
            _history=(),
            history_depth=2,
        ),
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

    The dependency audit + retirement sequence eliminated the
    packed-1D vector machinery this function previously consumed:

    * ``build_equation_map_*`` — **RETIRED** in D-J (2026-05-30).
    * ``solution_to_angular_flux_*`` — **RETIRED** in D-J (2026-05-30).
    * ``pack_with_traces`` — **RETIRED** in D-J (2026-05-30).
    * ``_build_rhs_{cartesian, spherical, cylindrical}`` —
      **RETIRED** in P1.7 of the moment-space + layering plan (the
      Cartesian variant carried an inline ``(2*l+1)`` literal that
      duplicated
      :attr:`~orpheus.numerics.spaces.SphericalHarmonicSpace.addition_theorem_factor`;
      all three were already dead code in production).
    * ``_make_sweep_preconditioner`` — **RETIRED** pre-D-J.
    * ``SNStreamingOperator`` — **RETIRED** in D-K (commit
      ``dadf4e8``).  ``self.L`` now binds to the algebraic
      composition ``StreamingOperator + CollisionOperator``
      (= :class:`InvertibleOperator`).

    Scope
    =====

    1-D meshes only (slab + sphere + cylinder).  2-D Cartesian fixed-
    source Krylov is not yet wired/verified, so it raises
    :class:`NotImplementedError` — but the geometry-agnostic
    :func:`_solve_fixed_source_si` is the working 2-D fixed-source path
    (it routes through :func:`transport_sweep`'s wavefront dispatch and
    handles 2-D natively, exactly as the now-live 2-D eigenvalue SI
    inner does).  Use ``inner_solver="source_iteration"`` for 2-D
    Cartesian fixed-source.  (Un-gating 2-D fixed-source Krylov is a
    small follow-on: the eigenvalue Krylov inner :meth:`SNSolver.
    _solve_krylov` already solves the 2-D operator, so the same
    machinery applies — it needs its own de-risk + SI≡Krylov equivalence
    pin before the guard is lifted.)
    """
    # 2-D Cartesian deferral — fixed-source Krylov for 2-D is not yet
    # verified.  This is NOT the dead "B1'' face block" excuse (that
    # legacy 1-D boundary closure was retired by the bare-boundary
    # architecture — see :meth:`SNSolver._solve_source_iteration`): the
    # 4-face machinery exists and the 2-D eigenvalue Krylov path uses it.
    # The deferral is simply that this fixed-source Krylov variant lacks
    # its own verification; the geometry-agnostic SI path is the
    # principled 2-D fixed-source landing zone (SURPRISE-5 of the
    # dependency audit).
    if sn_mesh.reduced is None:
        raise NotImplementedError(
            "R-1 Step 4 G1 — 2-D Cartesian fixed-source Krylov is not "
            "yet wired/verified.  Use `inner_solver=\"source_iteration\"` "
            "for 2-D Cartesian fixed-source (`_solve_fixed_source_si` is "
            "geometry-agnostic and handles 2-D via `transport_sweep`'s "
            "wavefront dispatch)."
        )

    from orpheus.transport.fields.angular_flux import (
        AngularFlux,
    )
    from orpheus.transport.fields.boundary_flux import (
        BoundaryFlux,
    )
    from orpheus.transport.fields.scalar_flux import ScalarFlux
    from .operator import CollisionOperator, StreamingOperator
    from orpheus.transport.source_sinks import (
        AngularSourceSink,
        BoundarySourceSink,
    )
    from orpheus.transport.timed_full_field import TimedFullField
    from orpheus.numerics.iteration import KrylovAcceleration
    from orpheus.numerics.operator import ZeroOperator

    nx, ny, ng = sn_mesh.nx, sn_mesh.ny, solver.ng
    N = sn_mesh.quad.N

    # ── Build the composite RHS ──────────────────────────────────────
    # R-1 Step 4 A1 — ``external_source`` is per-ordinate density (the
    # producer-side ``/sum_w`` projection lives at
    # :meth:`AngularSourceSink.from_isotropic` for iso scalar sources;
    # by the time we get here the caller has projected).
    # D-H.1c stage 2 — TimedFullField composite for the path-forward
    # Krylov (zero boundary; reflective-BC state flows through
    # ``initial_guess`` not ``rhs.boundary``).
    q_ext_per_ord = AngularSourceSink.from_mesh(external_source, sn_mesh)
    q_ext_composite = TimedFullField(
        # B.5.2: q_ext IS a source — carry the AngularSourceSink bulk AND the
        # BoundarySourceSink inflow trace (zero for vacuum/reflective). The
        # Krylov matvec composes operator-output sources; q_ext is raveled
        # type-agnostically as the GMRES rhs `b`.
        bulk=q_ext_per_ord,
        boundary=BoundarySourceSink.zeros_on(sn_mesh),
        _history=(),
        history_depth=2,
    )

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
        # Wave O #208 — ``S + B`` fold (see SNSolver._scattering_with_boundary_op):
        # the matvec ``(L+C).apply − (S+B).apply − F.apply`` IS ``(L+C−S−F−B)``.
        solver._scattering_with_boundary_op,
        ZeroOperator(codomain_zero=_zero_within_group_fission),
        preconditioner=lambda q: q,
        tol=inner_tol,
        max_iter=max_inner,
        # D-H.1e (2026-05-28): see :meth:`SNSolver._solve_krylov`'s
        # matching note.  Restart at full problem size; the legacy
        # ``min(50, ...)`` clamp left GMRES structurally truncated on
        # any mesh with ``N*ng*nx*ny > 50``.  ERR-053.
        restart=N * ng * nx * ny,
    )

    # B.5.2: pass a FLUX initial_guess so the Krylov solution_template (the
    # return type) is a flux; x0 stays all-zeros (bit-identical to the prior
    # initial_guess=None cold start).
    psi_typed, residuals = krylov.solve(
        q_ext_composite,
        initial_guess=TimedFullField.zeros(
            bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh,
        ),
    )
    # D-H.1c stage 2 — psi_typed is a TimedFullField (the Krylov
    # ravellable protocol unravels back to the template type, which is
    # ``q_ext_composite``).  Read bulk for scalar reduction.
    phi = psi_typed.bulk.integrate_angular().values
    converged_flag = bool(residuals) and residuals[-1] < inner_tol
    n_outer = len(residuals)
    flux_residuals = [float(r) for r in residuals]

    # Issue #197 PR-TYPED-5: build typed Solution at the boundary.
    # R-1 Step 4 G1 — ``psi_typed`` carries the Krylov-converged
    # composite with the matvec's B1'' face residual on its boundary;
    # reuse directly.
    del mesh, quadrature, materials  # retained as kwargs for API stability
    history = IterationHistory(
        flux_residuals=tuple(flux_residuals),
        n_inner=n_outer + 1,
        converged=converged_flag,
    )
    # D-H.1c stage 2 (2026-05-28): psi_typed IS already a TimedFullField;
    # no adapter wrap at the Solution boundary.
    return Solution(
        angular_flux=psi_typed,
        scalar_flux=ScalarFlux.from_mesh(phi, sn_mesh),
        mesh=sn_mesh,
        keff=None,
        history=history,
    )
