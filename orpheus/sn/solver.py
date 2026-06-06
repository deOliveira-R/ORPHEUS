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
from collections.abc import Iterable
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
from orpheus.transport.timed_full_field import TimedFullField


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


def _within_group_triple(solver: "SNSolver") -> tuple:
    r"""The within-group loss decomposition ``(L + C, S, B)`` — the invertible
    resolvent plus its two lagged coupling gains.

    Single source of truth (Cardinal Rule 2 / coding-elegance Pattern 2 + 5)
    for the operators EVERY within-group solve consumes — eigenvalue SI
    (:meth:`SNSolver._solve_source_iteration`), eigenvalue Krylov
    (:meth:`SNSolver._solve_krylov`), fixed-source SI
    (:func:`_solve_fixed_source_si`), and fixed-source Krylov
    (:func:`_solve_fixed_source_krylov`).  The four solves differ ONLY in the
    driver (:class:`SourceIteration` vs :class:`KrylovAcceleration`), the
    ``q_ext`` (fission source vs external source), and the returned contract —
    NOT in this decomposition.

    The driver consumes them via the variadic ``Driver(resolvent, *gains)``
    shape (Wave O #208 O.2a): the matvec is the honest
    ``(L + C).apply − S.apply − B.apply`` ≡ ``(L + C − S − B)·ψ`` and the SI
    rhs is ``q_ext + S.apply(ψ) + B.apply(ψ)``.

    * ``L + C`` is the :class:`~orpheus.sn.operator.InvertibleOperator`
      (``.apply`` = matvec, ``.solve`` = the WDD sweep) — the resolvent.
    * ``S`` (:class:`~orpheus.sn.scattering.ScatteringOperator`) is the BULK
      scattering coupling.  Producer-side ``/W`` normalisation (R-1 Step 4 A1)
      lives inside :meth:`ScatteringOperator.apply`; no consumer-side rescale.
    * ``B`` (:class:`~orpheus.sn.boundary_operator.SNBoundaryOperator`) is the
      BOUNDARY reflective coupling, delivered as a SEPARATE first-class gain
      (Wave O #208 O.2a — the transitional ``S + B`` fold is RETIRED).  Its
      ``B·ψ.outflow`` lands on ``rhs.boundary``, which the bare
      ``(L + C).solve`` sweep reads as the inflow seed.  ``B`` stays separate
      from ``S`` because it lives on the trace — it cannot join the ``L + C``
      preconditioner (:class:`OperatorSum` drops ``CAP_SOLVE``) and the
      cosine-weighted ``|Ω·n|·w`` adjoint metric (Wave O O.2) lives on ``B``'s
      trace domain, not the bulk.

    Within-group fission is zero (it enters as ``q_ext`` per the eigenvalue
    outer / within-group decomposition), so there is no fission gain here.
    """
    from .operator import CollisionOperator, StreamingOperator
    from .boundary_operator import SNBoundaryOperator

    sn_mesh = solver.sn_mesh
    L = StreamingOperator(sn_mesh, solver.mat_xs.total_cross_section)
    C = CollisionOperator(sn_mesh, solver.mat_xs.total_cross_section)
    return (
        L + C,
        solver.scattering_op,
        SNBoundaryOperator(sn_mesh),
    )


def _within_group_krylov(
    LC: "object", *gains: "object",
    n_dof: int, max_iter: int, tol: float,
):
    r"""GMRES driver on the within-group loss operator ``(L + C − S − B)``.

    Single source of truth (Cardinal Rule 2 / Phase 1 R2) for the
    :class:`~orpheus.numerics.iteration.KrylovAcceleration` construction shared
    by the eigenvalue (:meth:`SNSolver._solve_krylov`) and fixed-source
    (:func:`_solve_fixed_source_krylov`) Krylov paths — they previously built
    byte-identical instances.  ``*gains`` are the lagged couplings (the
    scattering ``S`` and boundary reflection ``B`` from
    :func:`_within_group_triple`); the matvec is ``LC.apply − Σ gᵢ.apply``.

    GMRES is UNPRECONDITIONED (explicit identity) per `issue #200
    <https://github.com/deOliveira-R/ORPHEUS/issues/200>`_ (the block-inverse
    face preconditioner re-enablement).  ``restart`` is sized to the FULL
    problem ``n_dof = N·ng·nx·ny`` — the legacy ``min(50, …)`` clamp left GMRES
    structurally truncated on any mesh with ``n_dof > 50`` (ERR-053).
    """
    from orpheus.numerics.iteration import KrylovAcceleration

    return KrylovAcceleration(
        LC, *gains,
        preconditioner=lambda q: q,  # explicit identity — issue #200 tracks re-enable
        tol=tol, max_iter=max_iter,
        restart=n_dof,
    )


class _GaussSeidelResolvent:
    r"""SI-only resolvent folding the BOUNDARY reflection ``B`` into the 2-D
    wavefront sweep via an octant-group Gauss-Seidel ``SweepSchedule`` (Phase 3
    sub-step 3c).

    The plain SI resolvent is ``(L+C)`` with ``B`` lagged as an external gain
    (inter-sweep Jacobi — ``B·ψₙ`` rides ``rhs.boundary``).  This resolvent
    INSTEAD seeds ``B·ψₙ`` itself and re-reflects each octant group's outgoing
    reflective faces BETWEEN group sweeps, so a later group reads the fresh
    current-iterate inflow — the ``(L+C−B_lower)⁻¹`` forward substitution that
    recovers the intra-sweep reflective coupling Wave O externalised.

    Honest SCOPE (Phase 3 spike — issues #2 / #215): this folds the BOUNDARY
    coupling only — a MODEST reflective-SI rate gain (~0.86–0.92× on the B-2g
    configs; measured).  The dominant within-group SCATTERING ``c``-mode is NOT
    folded (it cannot be folded into a directional sweep — that is consistent
    DSA / Krylov territory, #2).  The converged fixed point is IDENTICAL to the
    Jacobi SI (only the spectral rate changes — ``vv-principles`` Mode 9);
    Krylov is splitting-invariant and unaffected.

    2-D Cartesian ONLY (``_sweep_2d_scheduled``).  1-D meshes route to the
    Jacobi resolvent (:func:`_select_si_resolvent` never constructs this for
    1-D — boundary G-S is a no-op there AND the scan is not a wavefront).

    Satisfies the :class:`~orpheus.numerics.iteration.SourceIteration` ``L``
    contract: ``capabilities`` advertises ``{CAP_APPLY, CAP_SOLVE}`` and
    ``solve(rhs, *, initial_guess=…)`` runs the seed-then-overwrite loop.
    ``apply`` delegates to ``(L+C)`` — it is NEVER called by SourceIteration
    (which only invokes ``.solve``); the boundary G-S is purely a solve-time
    concern.
    """

    def __init__(self, invertible, boundary_op, schedule) -> None:
        self._invertible = invertible    # (L+C) InvertibleOperator at σ_t
        self._boundary_op = boundary_op  # SNBoundaryOperator (the same B)
        self._schedule = schedule        # SweepSchedule.gauss_seidel(sn_mesh)
        self.sn_mesh = invertible.sn_mesh

    @property
    def capabilities(self) -> frozenset[str]:
        from orpheus.numerics.operator import CAP_APPLY, CAP_SOLVE

        return frozenset({CAP_APPLY, CAP_SOLVE})

    def apply(self, psi):
        # Unused by SourceIteration (it calls only .solve); delegate to (L+C)
        # so the CAP_APPLY contract holds. The boundary G-S is solve-time only.
        return self._invertible.apply(psi)

    def solve(self, rhs, *, initial_guess=None):
        r"""``(L+C−B_lower)⁻¹ rhs`` via the seed-then-overwrite G-S sweep.

        Seeds ``boundary_buf = rhs.boundary + B·ψₙ`` (the lagged whole-trace
        reflection of the previous iterate — the SAME seed the Jacobi path gets
        via the external ``B`` gain, only here ``B`` lives in the resolvent),
        then runs :func:`_sweep_2d_scheduled` with the G-S schedule and the
        face-restricted ``−B`` reflect between octant-group sweeps.
        """
        from orpheus.transport.fields.angular_flux import AngularFlux
        from orpheus.transport.fields.boundary_flux import BoundaryFlux
        from orpheus.transport.timed_full_field import TimedFullField
        from .sweep import _sweep_2d_scheduled

        sn_mesh = self.sn_mesh
        trace = sn_mesh.trace

        # boundary_buf = rhs.boundary (external inflow) + B·ψₙ (lagged reflect
        # of the previous iterate's outflow).
        boundary_buf = BoundaryFlux.zeros_on(sn_mesh)
        for face in boundary_buf.layout.faces:
            if face in rhs.boundary.layout.faces:
                boundary_buf.face_view(face)[:] = rhs.boundary.face_view(face)
        if initial_guess is not None:
            seed = self._boundary_op.reflect_into_inflow(initial_guess.boundary)
            for face in boundary_buf.layout.faces:
                inflow = trace.inflow_indices_for_face(face)
                boundary_buf.face_view(face)[inflow] += seed.face_view(face)[inflow]

        # The per-group −B reflect (face-restricted, in place) — the SAME
        # whole-trace helper the reconstruction sweep uses, single-sourced.
        def _reflect(bf, faces):
            _reflect_outflow_into_inflow(bf, sn_mesh, faces=faces)

        angular, _scalar = _sweep_2d_scheduled(
            rhs.bulk.values,
            self._invertible.sigma,
            sn_mesh,
            boundary_buf,
            schedule=self._schedule,
            reflect=_reflect,
        )
        return TimedFullField(
            bulk=AngularFlux.from_mesh(angular, sn_mesh),
            boundary=boundary_buf,
            _history=(),
            history_depth=rhs.history_depth,
        )


def _select_si_resolvent(LC, S, B, sn_mesh, inner_schedule: str):
    r"""Pick the ``(resolvent, gains)`` for the within-group SI driver per
    ``inner_schedule`` — the single source of truth for the Jacobi/G-S choice.

    * ``"jacobi"`` (or any 1-D mesh) → ``(L+C, (S, B))``: ``B`` lagged as an
      external gain (inter-sweep Jacobi — today's path, every geometry).
    * ``"gauss_seidel"`` on a 2-D Cartesian mesh → ``(GaussSeidelResolvent, (S,))``:
      ``B`` folded INTO the resolvent (the octant-group G-S forward
      substitution).  ``S`` stays a lagged gain in BOTH (only the boundary
      coupling gets G-S; the sweep never re-scatters mid-sweep).

    1-D falls back to Jacobi: boundary G-S is a no-op on the scattering-
    dominated 1-D regime AND the 1-D scan is not a wavefront
    (:func:`_GaussSeidelResolvent` is 2-D-only).  The converged fixed point is
    identical either way — this only selects the SI spectral rate.
    """
    if inner_schedule not in ("jacobi", "gauss_seidel"):
        raise ValueError(
            f"Unknown inner_schedule: {inner_schedule!r}. "
            f"Valid choices are 'gauss_seidel' (boundary G-S, 2-D Cartesian) "
            f"or 'jacobi' (the splitting-invariant control)."
        )
    if inner_schedule == "gauss_seidel" and sn_mesh.reduced is None:
        from .sweep_schedule import SweepSchedule

        schedule = SweepSchedule.gauss_seidel(sn_mesh)
        return _GaussSeidelResolvent(LC, B, schedule), (S,)
    return LC, (S, B)


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

        # Phase 3 measurement seam: total inner (within-group) SI/Krylov
        # iterations consumed across the eigenvalue outer loop — the
        # measurand for the SI spectral-rate / Gauss-Seidel-recovery
        # diagnostics.  Accumulated per outer step in
        # ``_solve_source_iteration`` / ``_solve_krylov`` and read by
        # ``solve_sn`` into ``IterationHistory.total_inner_iterations``.
        # Fresh per solve (a new SNSolver is built per ``solve_sn`` call).
        self._total_inner_iterations = 0

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
        identical composite RHS, identical loss decomposition
        (``LC = StreamingOperator + CollisionOperator`` plus the
        scattering ``S`` and boundary ``B`` coupling gains —
        :func:`_within_group_triple`, zero within-group fission), identical
        ``psi_typed.bulk.integrate_angular()`` reduction — differing
        ONLY in the driver (:class:`SourceIteration` vs
        :class:`KrylovAcceleration`), and neither driver carries any
        geometry dependence.  The reflective coupling rides the BARE
        sweep via the ``B`` coupling gain on the 4-face
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
        # the sibling ``−B`` delivered as a SEPARATE coupling gain (O.2a —
        # ``_within_group_triple`` returns ``(L+C, S, B)``): each SI iterate
        # adds ``B·ψ.outflow`` to ``rhs.boundary``, which ``(L+C).solve``'s
        # bare sweep reads as the inflow seed (operator.py
        # ``_solve_timed_full_field`` seeds from ``rhs.boundary``).  The
        # boundary inflow is thus a live solved unknown carried in
        # ``ψ.boundary``, not an externally-recomputed partner trace.
        q_ext_composite = TimedFullField(
            # B.5.2: q_ext IS a source — carry the AngularSourceSink bulk AND
            # the BoundarySourceSink inflow trace (zero for vacuum/reflective;
            # prescribed inflow otherwise). The SI rhs q_ext + S.apply + B.apply
            # closes on BoundarySourceSink (operator outputs are sources).
            bulk=q_ext_per_ord,
            boundary=BoundarySourceSink.zeros_on(self.sn_mesh),
            _history=(),
            history_depth=2,
        )

        # ── Build the within-group operator triple (single source of
        # truth — ``_within_group_triple``; shared with the Krylov and
        # fixed-source paths). ───────────────────────────────────────
        LC, S, B = _within_group_triple(self)
        si = SourceIteration(
            LC, S, B,
            max_iter=self.max_inner, tol=self.inner_tol,
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
        # Phase 3 measurement seam: accumulate this outer step's inner SI
        # iterate count (the eigenvalue path's only window onto the inner
        # spectral rate — see IterationHistory.total_inner_iterations).
        self._total_inner_iterations += len(_residuals)
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
            # prescribed inflow otherwise). The SI rhs q_ext + S.apply + B.apply
            # closes on BoundarySourceSink (operator outputs are sources).
            bulk=q_ext_per_ord,
            boundary=BoundarySourceSink.zeros_on(self.sn_mesh),
            _history=(),
            history_depth=2,
        )

        # ── Build the within-group operator triple + Krylov driver
        # (single source of truth — ``_within_group_triple`` /
        # ``_within_group_krylov``; shared with the SI and fixed-source
        # paths). ─────────────────────────────────────────────────────
        LC, S, B = _within_group_triple(self)
        krylov = _within_group_krylov(
            LC, S, B,
            n_dof=self.quad.N * self.ng * self.sn_mesh.nx * self.sn_mesh.ny,
            max_iter=self.max_inner, tol=self.inner_tol,
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
        # Phase 3 measurement seam: accumulate this outer step's inner
        # Krylov iterate count (see IterationHistory.total_inner_iterations).
        self._total_inner_iterations += len(_residuals)
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


def _reflect_outflow_into_inflow(
    boundary_flux, sn_mesh: SNMesh, faces: "Iterable[str] | None" = None,
) -> None:
    r"""In-place: fill each face's inflow ordinate slots with the realized
    boundary law applied to that face's outflow trace — the ``−B`` reflective
    coupling, externalised for the bare ``transport_sweep`` (Wave O #208 O.4a.2).

    The bare sweep reads the inflow ordinate slots of its boundary buffer as
    the inflow seed; it no longer re-applies ``bc`` to the outflow internally.
    The SI driver path supplies ``B·ψ.outflow`` through ``rhs.boundary`` (as the
    ``B`` coupling gain), but the DIRECT fixed-source SI loop
    (:func:`_solve_fixed_source_si`) and the final eigenvalue reconstruction
    sweep (:func:`solve_sn`) do not route through that driver — they call this
    helper to set ``ψ.inflow = B·ψ.outflow`` on the buffer before each sweep,
    via the canonical whole-trace :class:`~orpheus.sn.boundary_operator.SNBoundaryOperator`
    (single source of truth — the same ``B`` the matvec / SI driver consume).

    For vacuum ``B = 0`` so the inflow slots stay zero (bit-identical to the
    pre-extraction ``bc.apply`` of a vacuum law); for reflective/white/albedo
    it is the same ``R·G`` reflection the pre-extraction sweep applied at entry,
    merely relocated to the caller.

    This is the TWIN of the driver route: the within-group SI/Krylov drivers
    deliver the SAME ``−B`` coupling via the SAME
    :class:`~orpheus.sn.boundary_operator.SNBoundaryOperator`, as a first-class
    coupling gain (Wave O O.2a — :func:`_within_group_triple` returns
    ``(L+C, S, B)``).  The two differ only in plumbing: this helper writes the
    buffer's inflow slots directly; the driver's ``B`` gain rides
    ``B·ψ.outflow`` in ``rhs.boundary``.  The DRIVER route no longer needs this
    helper (O.2a collapsed it); it survives for the final eigenvalue
    reconstruction sweep (which has no driver to route through) AND Phase 3's
    octant-group Gauss-Seidel resolvent (which calls it face-restricted between
    octant-group sweeps — see ``faces``).

    ``faces`` (Phase 3 G-S): ``None`` (default) reflects EVERY boundary face —
    the whole-trace ``−B`` used by the reconstruction sweep + the SI seed. A
    face subset restricts the reflect to those faces' inflow slots, leaving the
    rest untouched: the G-S resolvent absorbs a just-swept group's outgoing
    reflective faces into ``boundary_flux``, then calls this to reflect ONLY
    those faces' outflow into inflow, so the next group reads the fresh
    current-iterate inflow. ``B`` is block-diagonal over faces ⟹ exact
    restriction.
    """
    from orpheus.sn.boundary_operator import SNBoundaryOperator

    # Trace-only ``A_ss`` action — no zero-bulk probe (the bulk was only ever a
    # carrier to reach ``B``'s boundary block). ``reflect_into_inflow`` is the
    # canonical trace-level entry; it shares ``_reflect_trace`` with ``B.apply``
    # so the helper and the matvec / SI driver cannot drift.
    reflected = SNBoundaryOperator(sn_mesh).reflect_into_inflow(
        boundary_flux, faces=faces,
    )
    trace = sn_mesh.trace
    selected = boundary_flux.layout.faces if faces is None else faces
    for face in selected:
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
        total_inner_iterations=solver._total_inner_iterations,
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


def _build_fixed_source_rhs(
    external_source: "np.ndarray | TimedFullField",
    sn_mesh: SNMesh,
) -> TimedFullField:
    r"""Normalize the external source into the composite RHS ``q = q_bulk ⊕ q_∂``.

    A fixed-source problem's RHS is the composite source
    :class:`~orpheus.transport.timed_full_field.TimedFullField` — a bulk
    :class:`~orpheus.transport.source_sinks.AngularSourceSink` paired with a
    boundary :class:`~orpheus.transport.source_sinks.BoundarySourceSink`
    (the prescribed inflow :math:`q` of the affine BC
    :math:`\gamma_-\psi = R\,G\,\gamma_+\psi + q`). This is the ONE object
    that represents a source everywhere in the solve; this helper is its
    single construction point (Cardinal Rule 2 — the SI and Krylov inner
    paths both consume what it returns, rather than each re-deriving it).

    ``external_source`` is accepted in two forms:

    * ``np.ndarray`` of shape ``(N, ng, nx, ny)`` — the per-ordinate-density
      BULK source only; the boundary is vacuum (all-zero). The original
      form; every pre-existing caller keeps working unchanged.
    * :class:`TimedFullField` — the full COMPOSITE source (bulk + a possibly
      non-vacuum prescribed-inflow boundary, e.g. from
      :meth:`BoundarySourceSink.prescribed_inflow`). Its leaf values are
      re-homed onto ``sn_mesh``: the trace/grid layout is deterministic from
      ``(mesh, quadrature, materials)``, so this is an exact values-copy onto
      the solve's own mesh instance — required because the within-group
      operators are built on ``sn_mesh`` and :class:`TimedFullField` algebra
      enforces mesh identity.
    """
    from orpheus.transport.source_sinks import (
        AngularSourceSink,
        BoundarySourceSink,
    )

    N = sn_mesh.quad.N
    nx, ny, ng = sn_mesh.nx, sn_mesh.ny, sn_mesh.ng
    expected = (N, ng, nx, ny)

    if isinstance(external_source, TimedFullField):
        bulk_values = np.asarray(external_source.bulk.values)
        trace_size = int(sn_mesh.trace.layout.total_size)
        boundary_values = external_source.boundary.values
        if boundary_values.size != trace_size:
            raise ValueError(
                f"_build_fixed_source_rhs: composite boundary source has "
                f"{boundary_values.size} values, but sn_mesh.trace expects "
                f"{trace_size} (layout mismatch — the composite must be built "
                f"on the same mesh / quadrature / materials)."
            )
        boundary = BoundarySourceSink.from_mesh(boundary_values.copy(), sn_mesh)
    else:
        bulk_values = np.asarray(external_source)
        if bulk_values.dtype == object:
            # A stray non-array, non-TimedFullField object (e.g. a bare
            # AngularSourceSink) — np.asarray wraps it as a 0-d object array.
            # Reject explicitly rather than failing the shape check obscurely.
            raise TypeError(
                f"_build_fixed_source_rhs: external_source must be an "
                f"(N, ng, nx, ny) array (bulk-only / vacuum) or a "
                f"TimedFullField composite source; got "
                f"{type(external_source).__name__}"
            )
        boundary = BoundarySourceSink.zeros_on(sn_mesh)

    # Issue #196 PR-INDEX-5: bulk source principled (N, ng, nx, ny).
    if bulk_values.shape != expected:
        raise ValueError(
            f"fixed-source bulk shape {bulk_values.shape} does not match "
            f"(N, ng, nx, ny) = {expected}"
        )
    return TimedFullField(
        bulk=AngularSourceSink.from_mesh(bulk_values, sn_mesh),
        boundary=boundary,
    )


def solve_sn_fixed_source(
    materials: dict[int, Mixture],
    mesh: Mesh1D | Mesh2D,
    quadrature: AngularQuadrature,
    external_source: "np.ndarray | TimedFullField",
    boundary_condition: str = "vacuum",
    scattering_order: int = 0,
    max_inner: int = 1000,
    inner_tol: float = 1e-12,
    inner_solver: str | None = None,
    inner_schedule: str = "gauss_seidel",
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
    external_source : (N, ng, nx, ny) ndarray OR TimedFullField
        The fixed-source RHS, in either of two forms (normalised by
        :func:`_build_fixed_source_rhs`):

        * ``np.ndarray`` of shape ``(N, ng, nx, ny)`` — the per-ordinate
          volumetric BULK source :math:`Q^{\text{ext}}_n(x)` in
          **per-ordinate density magnitude** (R-1 Step 4 A1 convention),
          with a **vacuum** boundary. Callers with an iso scalar source
          :math:`Q(\vec r, g)` should project to per-ordinate via
          :meth:`~orpheus.transport.source_sinks.AngularSourceSink.from_isotropic`
          before passing (the :math:`1/W` projection lives at the producer
          boundary per Pattern 7). The sweep does NOT apply ``/W``
          internally. Issue #196 PR-INDEX-5: principled layout (``g`` axis
          after ``N``).
        * :class:`~orpheus.transport.timed_full_field.TimedFullField` — the
          full **composite** source ``q = q_bulk ⊕ q_∂`` (a bulk
          :class:`~orpheus.transport.source_sinks.AngularSourceSink` paired
          with a boundary
          :class:`~orpheus.transport.source_sinks.BoundarySourceSink`). This
          is how a **non-vacuum prescribed inflow** is supplied — build the
          boundary via
          :meth:`~orpheus.transport.source_sinks.BoundarySourceSink.prescribed_inflow`
          (the affine-BC inhomogeneous term :math:`q`, consumed by the sweep
          as the inflow seed). The legacy array form is exactly the
          bulk-only / vacuum special case of this composite.
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
    inner_schedule : {"gauss_seidel", "jacobi"}
        Source-iteration BOUNDARY splitting (Phase 3, ``inner_solver=
        "source_iteration"`` only).  ``"gauss_seidel"`` (default) folds the
        reflective coupling ``B`` into an octant-group Gauss-Seidel resolvent
        (2-D Cartesian) — re-reflecting each octant group's outgoing reflective
        faces between group sweeps so a later group reads the fresh
        current-iterate inflow (a modest reflective-SI rate gain, ~0.86–0.92×
        on B-mixture configs).  ``"jacobi"`` lags ``B`` fully (the
        splitting-invariant control).  The converged fixed point is IDENTICAL
        for both — this selects only the SI spectral rate.  1-D meshes always
        fall back to Jacobi (boundary G-S is a no-op on the scattering-
        dominated 1-D regime, and the scan is not a wavefront).  The dominant
        within-group SCATTERING rate is unchanged either way — that is Krylov /
        consistent-DSA territory (issue #2).

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

    # Normalise the external source (raw bulk array OR composite
    # TimedFullField) into the single composite RHS ``q = q_bulk ⊕ q_∂`` the
    # inner paths consume (Cardinal Rule 2 — one construction point; shape
    # validation lives inside the helper).
    q_ext_composite = _build_fixed_source_rhs(external_source, sn_mesh)

    if inner_solver == "source_iteration":
        return _solve_fixed_source_si(
            solver, sn_mesh, q_ext_composite, mesh, quadrature, materials,
            t_start, max_inner, inner_tol, inner_schedule=inner_schedule,
        )

    # Krylov path.  We solve T·ψ = b directly via GMRES, where b carries
    # the external per-ordinate source plus any in-scatter / (n,2n) terms
    # built from the converged scalar flux.  Wrapping that in an outer
    # source iteration converges scattering self-consistently.
    return _solve_fixed_source_krylov(
        solver, sn_mesh, q_ext_composite, mesh, quadrature, materials,
        t_start, max_inner, inner_tol,
    )


def _solve_fixed_source_si(
    solver: SNSolver,
    sn_mesh: SNMesh,
    q_ext_composite: TimedFullField,
    mesh: Mesh1D | Mesh2D,
    quadrature: AngularQuadrature,
    materials: dict[int, Mixture],
    t_start: float,
    max_inner: int,
    inner_tol: float,
    inner_schedule: str = "gauss_seidel",
) -> Solution:
    r"""Fixed-source path via the :class:`SourceIteration` primitive.

    Phase 1 (R1, 2026-06-04) — carved onto the SAME
    :class:`~orpheus.numerics.iteration.SourceIteration` primitive the
    eigenvalue inner :meth:`SNSolver._solve_source_iteration` uses, on the
    identical operator triple

    .. math::

        A \;=\; L + C\,, \quad
        S \;=\; \text{full multi-group scatter} + B\,, \quad
        F \;=\; 0_{\rm wg}

    differing ONLY in ``q_ext`` (the EXTERNAL source here vs the fission
    source in the eigenvalue inner) and the returned contract (a full typed
    :class:`Solution` here vs an angular-integrated scalar flux there).

    The previous hand-rolled ``for n_inner`` loop is RETIRED.  It rebuilt the
    scattering source each iterate via ``_add_scattering_source`` +
    ``_add_n2n_source`` + ``_build_aniso_scattering`` and drove ``−B`` through
    :func:`_reflect_outflow_into_inflow` — both are now subsumed by the
    primitive's coupling gains ``S.apply(ψ_n) + B.apply(ψ_n)``:
    :meth:`ScatteringOperator.apply` (the ``TimedFullField`` branch)
    recomputes the IDENTICAL ``(P0 in-scatter + (n,2n))/W + Pℓ`` bulk source,
    and the boundary gain ``B``
    (:class:`~orpheus.sn.boundary_operator.SNBoundaryOperator`, a first-class
    coupling in :func:`_within_group_triple`) delivers the reflective
    ``B·ψ.outflow`` through ``rhs.boundary`` which the bare ``(L + C).solve``
    sweep reads as the inflow seed (single source of truth — Cardinal Rule 2).  The whole-trace
    :func:`_reflect_outflow_into_inflow` route is no longer needed on this
    path; it survives for the eigenvalue reconstruction sweep + Phase 3's
    octant-restricted Gauss-Seidel variant.

    Geometry-agnostic (slab / sphere / cylinder / 2-D Cartesian): the
    within-group solve carries no geometry dependence, exactly as the
    eigenvalue SI inner (Wave O "2-D SI Phase A").

    Equivalence note (vv-principles §bit-identity): the converged fixed point
    is identical to the retired loop's (same operators, same ``S`` and ``B``
    coupling gains, same WDD sweep), but the iteration TRAJECTORY differs — the primitive
    stops on the composite ``‖ψ_{n+1} − ψ_n‖ / ‖ψ_{n+1}‖`` residual (the full
    angular + boundary iterate, the same metric the verified eigenvalue inner
    uses) rather than the scalar-flux ``‖φ − φ_prev‖ / ‖φ‖``.  Converged ``φ``
    therefore agrees to ``~inner_tol`` (principled-equivalence), and
    ``history.n_inner`` / ``flux_residuals`` reflect the composite metric.
    """
    from orpheus.transport.fields.angular_flux import (
        AngularFlux,
    )
    from orpheus.transport.fields.boundary_flux import (
        BoundaryFlux,
    )
    from orpheus.transport.fields.scalar_flux import ScalarFlux
    from orpheus.numerics.iteration import SourceIteration

    # ``q_ext_composite`` is the normalised composite RHS ``q = q_bulk ⊕ q_∂``
    # built once by :func:`_build_fixed_source_rhs` (Cardinal Rule 2 — the SI
    # and Krylov paths share it).  The bulk is the per-ordinate-density
    # external source; the boundary is the prescribed inflow (zero for
    # vacuum / reflective — the reflective inflow rides ``rhs.boundary`` via
    # the ``B`` coupling gain, NOT ``q_ext``; a NON-vacuum prescribed inflow
    # is carried in ``q_ext_composite.boundary``).  Scattering (P0 + Pℓ +
    # (n,2n)) is NOT pre-staged — the primitive's ``S`` operator recomputes
    # it each iterate.

    # ── Build the within-group operator triple (single source of truth —
    # ``_within_group_triple``; identical to the eigenvalue inner). ────
    LC, S, B = _within_group_triple(solver)
    # Phase 3 sub-step 3c: ``inner_schedule`` selects the (resolvent, gains)
    # splitting.  "gauss_seidel" (default, 2-D Cartesian) folds B into a
    # scheduled resolvent → gains (S,); "jacobi" (or any 1-D mesh) keeps
    # resolvent (L+C) + gains (S, B) — today's inter-sweep Jacobi.  Same
    # converged fixed point; only the SI spectral rate differs.
    resolvent, gains = _select_si_resolvent(LC, S, B, sn_mesh, inner_schedule)
    si = SourceIteration(
        resolvent, *gains,
        max_iter=max_inner, tol=inner_tol,
    )

    # Cold-start FLUX iterate (x0 = zeros; the flux template fixes the return
    # type).  Fixed-source is a single solve — no eigenvalue outer to
    # warm-start from (cf. the eigenvalue inner's ``self._psi_typed``).
    psi_typed, residuals = si.solve(
        q_ext_composite,
        initial_guess=TimedFullField.zeros(
            bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh,
        ),
    )
    phi = psi_typed.bulk.integrate_angular().values
    converged_flag = bool(residuals) and residuals[-1] < inner_tol
    flux_residuals = [float(r) for r in residuals]

    # Issue #197 PR-TYPED-5: build typed Solution at the boundary.
    # mesh / quadrature / materials remain unconsumed by Solution — the
    # typed fluxes carry the SNMesh reference, which transitively exposes
    # those handles via .mesh.{mesh, quad, materials}.
    del mesh, quadrature, materials  # retained as kwargs for API stability
    history = IterationHistory(
        flux_residuals=tuple(flux_residuals),
        n_inner=len(residuals),
        total_inner_iterations=len(residuals),
        converged=converged_flag,
    )
    # ``psi_typed`` IS the converged TimedFullField composite (bulk angular +
    # boundary trace) — return it directly, exactly as the fixed-source Krylov
    # path does.  No legacy ``solver._boundary_flux`` writeback (the boundary
    # trace lives on ``psi_typed.boundary``).
    return Solution(
        angular_flux=psi_typed,
        scalar_flux=ScalarFlux.from_mesh(phi, sn_mesh),
        mesh=sn_mesh,
        keff=None,
        history=history,
    )


def _solve_fixed_source_krylov(
    solver: SNSolver,
    sn_mesh: SNMesh,
    q_ext_composite: TimedFullField,
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

    ALL geometries — slab, sphere, cylinder, AND 2-D Cartesian (Phase 1
    "gap" un-gate, 2026-06-04).  The 2-D Cartesian fixed-source Krylov is
    the structural twin of the live 2-D eigenvalue Krylov inner
    :meth:`SNSolver._solve_krylov`: identical operator triple
    (:func:`_within_group_triple`) and identical
    :class:`~orpheus.numerics.iteration.KrylovAcceleration` driver
    (:func:`_within_group_krylov`), differing ONLY in ``q_ext`` (the
    external source here vs the fission source there) — and it is the twin
    of the geometry-agnostic 2-D fixed-source SI path.  Verified: the
    converged per-ordinate flux hits the closed-form streaming equilibrium
    ``q/Σ_t`` on a homogeneous reflective box, and the SI≡Krylov flux
    shape agrees on a heterogeneous non-flat case
    (``tests/sn/solve/test_fixed_source_2d_equivalence.py``).  The legacy
    ``NotImplementedError`` 2-D guard (and the "B1'' face block" excuse it
    carried) is RETIRED — the 4-face bare-boundary architecture (Wave O
    O.4b) makes the path geometry-agnostic by construction.
    """
    from orpheus.transport.fields.angular_flux import (
        AngularFlux,
    )
    from orpheus.transport.fields.boundary_flux import (
        BoundaryFlux,
    )
    from orpheus.transport.fields.scalar_flux import ScalarFlux

    nx, ny, ng = sn_mesh.nx, sn_mesh.ny, solver.ng
    N = sn_mesh.quad.N

    # ``q_ext_composite`` is the normalised composite RHS ``q = q_bulk ⊕ q_∂``
    # built once by :func:`_build_fixed_source_rhs` (Cardinal Rule 2). B.5.2:
    # q_ext IS a source — bulk per-ordinate-density ``AngularSourceSink`` +
    # boundary ``BoundarySourceSink`` prescribed inflow (zero for vacuum /
    # reflective — the reflective inflow rides ``initial_guess`` /
    # ``rhs.boundary``, not ``q_ext``; a NON-vacuum prescribed inflow IS
    # carried in ``q_ext_composite.boundary``). The Krylov matvec composes
    # operator-output sources; q_ext is raveled type-agnostically as the
    # GMRES rhs ``b``.

    # ── Build the within-group operator triple + Krylov driver (single
    # source of truth — ``_within_group_triple`` / ``_within_group_krylov``;
    # shared with the eigenvalue Krylov + SI paths). ──────────────────
    LC, S, B = _within_group_triple(solver)
    krylov = _within_group_krylov(
        LC, S, B,
        n_dof=N * ng * nx * ny,
        max_iter=max_inner, tol=inner_tol,
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
        total_inner_iterations=n_outer + 1,
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
