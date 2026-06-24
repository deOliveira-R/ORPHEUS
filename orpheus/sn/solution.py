r"""Typed return type for the SN transport solvers.

Issue #197 PR-TYPED-5 — :class:`Solution` and :class:`IterationHistory`
replace the legacy bare-dataclass pair
:class:`~orpheus.sn.solver.SNFixedSourceResult` /
:class:`~orpheus.sn.solver.SNResult`.

The legacy types were data bags: they carried bare ``np.ndarray`` flux
fields, an opaque ``geometry`` / ``quadrature`` pair, and an ad-hoc set
of diagnostic scalars (``n_inner``, ``residual``, ``keff_history``,
``elapsed_seconds``).  Two distinct dataclasses for the two problem
kinds (fixed-source vs eigenvalue) duplicated the shape contract in two
places — a twin path waiting for drift (``coding-elegance`` anti-pattern
1: "two implementations of the same mathematical quantity").

Under PR-TYPED-5, **one** typed :class:`Solution` covers both problem
kinds.  The discrimination lives in two methods —
:meth:`Solution.is_eigenvalue` and :meth:`Solution.is_fixed_source` —
that read the optional :attr:`keff`.  The convergence-trajectory
diagnostics live on a separate :class:`IterationHistory`, exposed
through method-style accessors :meth:`Solution.dominance_ratio` and
:meth:`Solution.converged`.

Reads as the math (``coding-elegance`` Pattern 1 — match the algebra of
the domain)::

    sol = solve_sn(...)
    sol.is_eigenvalue()                    # True / False
    sol.dominance_ratio()                  # |k_n / k_{n-1}|
    sol.reaction_rate_density(sig_a)       # σ · φ
    sol.compare(other, rtol=1e-12)         # SolutionDiff

Mesh-binding consistency
========================

:class:`Solution.__post_init__` validates that every typed field
(:attr:`angular_flux`, :attr:`scalar_flux`, :attr:`boundary_flux`)
carries the SAME :class:`SNMesh` instance as the
:attr:`Solution.mesh`.  Mesh identity is checked via ``is`` — sharing a
mesh by value rather than by reference is forbidden by construction
(``coding-elegance`` Pattern 4 — illegal states unrepresentable).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .geometry import SNMesh
    from orpheus.geometry import Mesh1D
    from orpheus.transport.fields.boundary_flux import BoundaryFlux
    from orpheus.transport.fields.scalar_flux import ScalarFlux
    from orpheus.transport.mesh.material_mesh import MaterialMesh
    from orpheus.transport.timed_full_field import TimedFullField


__all__ = ["IterationHistory", "Solution", "SolutionDiff"]


@dataclass(frozen=True)
class IterationHistory:
    r"""Convergence trajectory diagnostics for an SN solve.

    Carries the per-outer / per-inner trajectory exposed by the legacy
    :class:`SNResult.keff_history` (list) plus the integer
    :attr:`SNFixedSourceResult.n_inner` + :attr:`residual` scalars — now
    bundled into one frozen container with method-style diagnostics.

    Parameters
    ----------
    keff_history : tuple of float
        Per-outer-iteration eigenvalue trajectory.  Empty tuple for
        fixed-source problems (no eigenvalue to track).
    flux_residuals : tuple of float
        Per-iteration relative-flux residual trajectory
        :math:`\lVert\phi_n - \phi_{n-1}\rVert / \lVert\phi_n\rVert`.
        Empty when the solver does not track this signal.
    n_inner : int or None
        Number of inner (within-group) iterations consumed.  ``None``
        for paths that do not surface this count (e.g. pure eigenvalue
        outer iteration).
    total_inner_iterations : int or None
        Total inner (within-group) iterations summed across ALL outer
        iterations.  Unlike :attr:`n_inner` (which is ``None`` on the
        eigenvalue outer path — there is no single inner solve to point
        at), this is populated by BOTH paths: the eigenvalue path
        accumulates it across the power-iteration outer loop, the
        fixed-source path sets it equal to :attr:`n_inner` (one inner
        solve).  It is the measurand for the SI spectral-rate /
        Gauss-Seidel-recovery diagnostics (Phase 3).  ``None`` only when
        no path surfaced a count.
    n_outer : int or None
        Number of outer (power) iterations consumed.  ``None`` for
        fixed-source problems.
    converged : bool
        ``True`` when the iteration met its convergence tolerance,
        ``False`` when ``max_iter`` was exhausted.

    Notes
    -----
    The trajectories are **tuple** rather than list — frozen dataclasses
    with mutable fields are an anti-pattern (``coding-elegance``
    Pattern 4: illegal states unrepresentable).  Callers that want a
    list pass through :func:`list` at the call site.
    """

    keff_history: tuple[float, ...] = ()
    flux_residuals: tuple[float, ...] = ()
    n_inner: int | None = None
    total_inner_iterations: int | None = None
    n_outer: int | None = None
    converged: bool = True

    def dominance_ratio(self) -> float | None:
        r"""Return :math:`|k_n - k_{n-1}| / |k_{n-1}|` in the late-iteration limit.

        The dominance ratio approximates the spectral gap
        :math:`|k_1 / k_0|` that controls power-iteration convergence.
        Returns ``None`` when fewer than 2 entries are recorded — the
        ratio is undefined for a single-point history.
        """
        if len(self.keff_history) < 2:
            return None
        prev = self.keff_history[-2]
        if prev == 0.0:
            return None
        return abs(self.keff_history[-1] - prev) / abs(prev)

    def latest_keff(self) -> float | None:
        """Return the last eigenvalue, or ``None`` for empty history."""
        return self.keff_history[-1] if self.keff_history else None

    def latest_residual(self) -> float | None:
        """Return the last recorded flux residual, or ``None`` if absent."""
        return self.flux_residuals[-1] if self.flux_residuals else None


@dataclass(frozen=True)
class Solution:
    r"""Canonical return type for :func:`solve_sn` / :func:`solve_sn_fixed_source`.

    Bundles the typed flux fields + boundary state + eigenvalue (when
    eigenvalue problem) + iteration trajectory into one frozen
    dataclass.  Replaces :class:`SNFixedSourceResult` + :class:`SNResult`
    (which were data bags without methods).  One type covers both
    fixed-source and eigenvalue problems via optional :attr:`keff` and
    :attr:`history`.

    Parameters
    ----------
    angular_flux : TimedFullField
        Composite iteration state — pure-Field
        :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
        bulk paired with pure-Field
        :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`
        boundary trace plus time-derivative history.

        Under D-H.1b (2026-05-28), this field carries a
        :class:`~orpheus.transport.timed_full_field.TimedFullField`
        composite. The bulk (``angular_flux.bulk``) is the
        per-ordinate angular flux :math:`\psi(\vec r, \hat\Omega_n, g)`
        on shape ``(N, ng, nx, ny)``; the boundary trace is exposed
        via the :attr:`Solution.boundary_flux` delegate property.
    scalar_flux : ScalarFlux
        Scalar flux field
        :math:`\phi(\vec r, g) = \sum_n w_n \psi_n`,
        shape ``(ng, nx, ny)``.
    mesh : SNMesh
        The phase-space mesh.  Both typed flux fields above MUST carry
        this same :class:`SNMesh` instance (validated at construction).
    keff : float or None
        Multiplication eigenvalue.  ``None`` for fixed-source problems.
    history : IterationHistory or None
        Convergence-trajectory diagnostics.  ``None`` for paths that
        do not surface them.

    Notes
    -----
    Pre-D-H, :attr:`angular_flux` was the legacy
    :class:`orpheus.sn.angular_flux.AngularFlux` (a bulk Field that
    ALSO owned the boundary face state and iteration history). Under
    D-H.1b that conflation dissolves: the composite-state container
    :class:`TimedFullField` holds the bulk + boundary + history trio
    as a structured composite. The :attr:`Solution.boundary_flux`
    delegate (line below) becomes a thin read-through to
    ``self.angular_flux.boundary`` (now a typed L2 BoundaryFlux).

    Examples
    --------
    Reads as the math:

    >>> sol = solve_sn(materials, mesh, quadrature)             # doctest: +SKIP
    >>> sol.is_eigenvalue()                                     # doctest: +SKIP
    True
    >>> sol.dominance_ratio()                                   # doctest: +SKIP
    1.2e-08
    >>> sol.reaction_rate_density(materials[0].sig_a)           # doctest: +SKIP
    array(...)
    """

    angular_flux: "TimedFullField"
    scalar_flux: "ScalarFlux"
    mesh: "SNMesh"
    keff: float | None = None
    history: IterationHistory | None = None

    def __post_init__(self) -> None:
        # Mesh-binding consistency: every typed field MUST share the
        # exact SNMesh instance (Pattern 4 — illegal states
        # unrepresentable; the typed field's invariants are leveraged
        # at Solution construction time so consumers downstream cannot
        # see a Solution with mismatched meshes).
        # D-H.1b: angular_flux is now a TimedFullField composite —
        # the mesh is on the bulk (and validated against boundary at
        # TimedFullField construction).
        if self.angular_flux.bulk.mesh is not self.mesh:
            raise ValueError(
                "Solution: angular_flux.bulk.mesh is not Solution.mesh "
                "(typed-field mesh-identity contract broken — every "
                "field must reference the same SNMesh instance)."
            )
        if self.scalar_flux.mesh is not self.mesh:
            raise ValueError(
                "Solution: scalar_flux.mesh is not Solution.mesh "
                "(typed-field mesh-identity contract broken — every "
                "field must reference the same SNMesh instance)."
            )

    # ── boundary_flux as a delegate property ─────────────────────────
    #
    # The composite :class:`~orpheus.transport.timed_full_field.TimedFullField`
    # owns the boundary face state via :attr:`TimedFullField.boundary`.
    # ``Solution.boundary_flux`` is a thin read-through that delegates
    # to the composite's owned boundary trace.

    @property
    def boundary_flux(self) -> "BoundaryFlux":
        r"""Boundary face state — delegate to :attr:`TimedFullField.boundary`.

        The composite :class:`~orpheus.transport.timed_full_field.TimedFullField`
        owns the boundary face state via :attr:`TimedFullField.boundary`.
        This property provides the canonical ``sol.boundary_flux``
        access path.

        See :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`
        for the per-geometry flat-layout contract.
        """
        return self.angular_flux.boundary

    # ── Discrimination ───────────────────────────────────────────────

    def is_eigenvalue(self) -> bool:
        """Return ``True`` when this is an eigenvalue-problem solution."""
        return self.keff is not None

    def is_fixed_source(self) -> bool:
        """Return ``True`` when this is a fixed-source-problem solution."""
        return self.keff is None

    # ── Convergence diagnostics (delegate to history) ────────────────

    def dominance_ratio(self) -> float | None:
        r"""Return the eigenvalue dominance ratio, or ``None`` if unavailable."""
        return self.history.dominance_ratio() if self.history else None

    def converged(self) -> bool:
        """Return whether the solver iteration met its convergence tolerance."""
        return self.history.converged if self.history else True

    def keff_history_list(self) -> list[float]:
        """Return the eigenvalue trajectory as a list (for plotting)."""
        return list(self.history.keff_history) if self.history else []

    @property
    def keff_history(self) -> list[float]:
        """Legacy alias for :meth:`keff_history_list`.

        The canonical access path is :attr:`history` + ``.keff_history``;
        this property keeps existing test fixtures (``len(result.keff_history)``,
        plotting code) working without a one-by-one migration of ~10 call
        sites.  The trajectory is a list, not a tuple, so ``len`` and
        slice operations behave as they did under the legacy
        :class:`SNResult.keff_history: list[float]`.
        """
        return self.keff_history_list()

    # ── Reaction-rate accessor (Pattern 1 — read as math) ─────────────

    def reaction_rate_density(self, xs: np.ndarray) -> np.ndarray:
        r"""Compute the per-cell reaction-rate density :math:`\sigma \cdot \phi`.

        Parameters
        ----------
        xs : np.ndarray
            Cross-section array, shape ``(ng, *spatial)`` (a
            :class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField`-shaped
            ndarray).  Per Issue #197 Wave 1 design, ``ReactionRate`` is
            a ``NewType`` over ``np.ndarray`` and is NOT promoted to a
            dataclass — staying close to the natural ``ng × *spatial``
            layout the consumers (depletion, k-eff diagnostics, response
            functionals) want directly.

        Returns
        -------
        np.ndarray
            Per-cell rate density, shape ``(ng, *spatial)``.  Units:
            :math:`\rm 1/(cm^3 \cdot s)` — the reaction-rate density
            per group.
        """
        return xs * self.scalar_flux.values

    # ── Spatial homogenization (a domain operation on the solution) ──

    def homogenize(self, coarse: "Mesh1D") -> "MaterialMesh":
        r"""Flux·volume-weighted spatial homogenization onto a coarse mesh.

        Collapse this (fine) solution's per-cell cross sections onto the
        cells of a coarser mesh, preserving every reaction rate.  Each
        coarse cell :math:`R` becomes its own effective material whose
        cross sections are the flux·volume-weighted averages of the fine
        cells it contains:

        .. math::

            \Sigma_{R,g} \;=\;
            \frac{\sum_{i \in R} V_i\,\phi_{i,g}\,\Sigma_{i,g}}
                 {\sum_{i \in R} V_i\,\phi_{i,g}}
            \qquad
            \Phi_{R,g} \;=\; \sum_{i \in R} V_i\,\phi_{i,g}

        so that the volume-integrated reaction rate is preserved exactly:
        :math:`\Sigma_{R,g}\,\Phi_{R,g} = \sum_{i \in R} V_i\,\Sigma_{i,g}\,
        \phi_{i,g}`.  This is the **space-only**, **mesh-COUPLED** half of
        the condense/homogenize asymmetry law: the coarse cells carry the
        homogenized materials, so geometry + materials are born together
        (returned as a :class:`MaterialMesh`).  Energy is **not** condensed
        — the group structure (``eg``) is carried through unchanged.

        The matrix channels (per-Legendre :math:`\Sigma_{s,\ell}[g',g]` and
        :math:`\Sigma_{2n}[g',g]`, indexed ``[g_from, g_to]``) weight by the
        **source** group :math:`g'` flux (the group whose population drives
        the out-scatter), and :math:`\chi` is the **production-weighted**
        convex average (weight :math:`p_i = \sum_g \nu\Sigma_{f,i,g}\,
        \phi_{i,g}\,V_i`) — a convex combination of simplices, hence a
        simplex (validated by :class:`Mixture.__post_init__`).  Because every
        removal channel collapses with the *same* per-(R, g) weight, the
        definitional total-XS balance :math:`\Sigma_t = \Sigma_c + \Sigma_L
        + \Sigma_f + \mathrm{rowsum}(\Sigma_{s0}) + \mathrm{rowsum}(\Sigma_{2n})`
        is preserved cell-by-cell when the fine materials balance.

        Parameters
        ----------
        coarse : Mesh1D
            The coarse target mesh.  Must share this solution's outer
            boundary; its internal cell edges must align with fine-cell
            edges (each coarse cell is a contiguous union of fine cells).
            Its own ``mat_ids`` are ignored — homogenization assigns one
            fresh effective material per coarse cell.

        Returns
        -------
        MaterialMesh
            The coarse mesh carrying the homogenized materials (one
            :class:`Mixture` per coarse cell, keyed by coarse-cell index).
            Promote to a solvable SN phase space with
            :meth:`~orpheus.sn.geometry.SNMesh.from_material_mesh`.

        Notes
        -----
        1-D only in this slice (the SN side, no vectorization loss — the
        group axis rides the integrand).  A group with zero region flux
        (:math:`\Phi_{R,g} = 0`) yields a zero effective XS for that
        (R, g) — there is no reaction rate to preserve there.
        """
        import dataclasses

        from scipy.sparse import csr_matrix

        from orpheus.geometry import Mesh1D
        from orpheus.data.macro_xs.mixture import Mixture
        from orpheus.numerics.frame import Frame
        from orpheus.transport.mesh.material_mesh import MaterialMesh

        fine = self.mesh
        if fine.ndim != 1:
            raise NotImplementedError(
                "Solution.homogenize is 1-D in this slice; got ndim="
                f"{fine.ndim}."
            )

        fine_edges = np.asarray(fine.axes[0].edges, dtype=float)
        coarse_edges = np.asarray(coarse.edges, dtype=float)
        if not (
            np.isclose(coarse_edges[0], fine_edges[0])
            and np.isclose(coarse_edges[-1], fine_edges[-1])
        ):
            raise ValueError(
                "Solution.homogenize: coarse mesh must share the fine "
                f"mesh's outer boundary [{fine_edges[0]}, {fine_edges[-1]}]; "
                f"got [{coarse_edges[0]}, {coarse_edges[-1]}]."
            )
        n_coarse = coarse_edges.size - 1
        ng = fine.ng

        # Homogenisation is the L²(φV)-orthogonal (Galerkin) projection of the
        # fine cross-section field onto the coarse cells — routed through the
        # discrete Frame, NOT hand-rolled.  K (trial) = the coarse mesh's
        # cell-indicator basis (which the mesh YIELDS); L (measure) = the fine
        # mesh's geometric volume measure.  The flux is NOT a measure weight:
        # by Radon–Nikodym (μ_φV = φ·μ_V) it is an integrand MULTIPLIER, so the
        # per-group flux rides a trailing tensor axis through ONE group-
        # independent frame.  ``analysis`` is the region integral ∫_R (·) dV;
        # the /Φ_R normalisation is the coarse Gram-inverse the coefficient
        # space carries as its metric (``apply_inverse_metric``, whose
        # Moore–Penrose pseudo-inverse zeroes empty / zero-flux regions).
        frame = Frame(coarse.indicator_basis(), fine.volume_measure)
        analysis = frame.analysis
        phi = np.asarray(self.scalar_flux.values, dtype=float).T       # (n_fine, ng)
        region_flux_integral = analysis.apply(phi)                     # Φ_{R,g}
        flux_space = dataclasses.replace(
            frame.basis_space, inner_product_weights=region_flux_integral,
        )

        mat_of_fine = np.asarray(fine.mat_map, dtype=int).ravel()      # (n_fine,)
        materials = fine.materials

        def _gather_vector(attr: str) -> np.ndarray:
            """Per-fine-cell view of a ``(ng,)`` Mixture channel — (n_fine, ng)."""
            return np.array([getattr(materials[m], attr) for m in mat_of_fine])

        def _collapse_vector(channel_fine: np.ndarray) -> np.ndarray:
            """Flux·volume-weighted collapse of a vector channel — (n_coarse, ng).

            The region reaction rate ``∫_R φΣ dV`` (``analysis`` of the
            flux-weighted integrand) over the region flux integral ``Φ_R`` (the
            coarse Gram-inverse, via the coefficient space's metric).
            """
            rate = analysis.apply(phi * channel_fine)
            return flux_space.apply_inverse_metric(rate)

        def _collapse_matrix(channel_fine: np.ndarray) -> np.ndarray:
            """Collapse a ``[g_from, g_to]`` channel, weighting by SOURCE group.

            ``channel_fine`` is ``(n_fine, ng, ng)``; the out-scatter from
            source group ``g_from`` scales with that group's flux, so the weight
            rides the ``g_from`` (first matrix) axis.  ``apply_inverse_metric``
            broadcasts the ``Φ_R[:, g_from]`` metric over the trailing ``g_to``
            axis (its trailing-axis metric-broadcast), so the source-group
            normalisation falls out for free.
            """
            weighted = phi[:, :, None] * channel_fine
            rate = analysis.apply(weighted)
            return flux_space.apply_inverse_metric(rate)

        sig_t = _collapse_vector(_gather_vector("SigT"))
        sig_c = _collapse_vector(_gather_vector("SigC"))
        sig_l = _collapse_vector(_gather_vector("SigL"))
        sig_f = _collapse_vector(_gather_vector("SigF"))
        sig_p = _collapse_vector(_gather_vector("SigP"))

        n_legendre = max(len(materials[m].SigS) for m in materials)

        def _gather_legendre(order: int) -> np.ndarray:
            """Per-fine-cell dense ``Σ_{s,ℓ}`` — (n_fine, ng, ng); zero-pad short lists."""
            return np.array([
                np.asarray(materials[m].SigS[order].todense())
                if order < len(materials[m].SigS) else np.zeros((ng, ng))
                for m in mat_of_fine
            ])

        sig_s = [_collapse_matrix(_gather_legendre(l)) for l in range(n_legendre)]
        sig2 = _collapse_matrix(
            np.array([np.asarray(materials[m].Sig2.todense()) for m in mat_of_fine])
        )

        # χ: production-weighted convex average (a convex combination of
        # simplices → a simplex).  Same projection, a different multiplier — the
        # production density p_i = Σ_g νΣ_{f,i,g} φ_{i,g}; routed through the
        # SAME ``analysis`` (which supplies the volume V_i) it gives the region
        # production Σ_{i∈R} V_i p_i, the Gram of the χ-collapse coarse space.
        chi_fine = _gather_vector("chi")
        prod_density = (_gather_vector("SigP") * phi).sum(axis=1)       # (n_fine,)
        region_production = analysis.apply(prod_density)               # (n_coarse,)
        prod_space = dataclasses.replace(
            frame.basis_space, inner_product_weights=region_production,
        )
        chi = prod_space.apply_inverse_metric(
            analysis.apply(prod_density[:, None] * chi_fine)
        )

        eg = next(iter(materials.values())).eg
        homogenized = {
            region: Mixture(
                SigC=sig_c[region], SigL=sig_l[region], SigF=sig_f[region],
                SigP=sig_p[region], SigT=sig_t[region],
                SigS=[csr_matrix(sig_s[l][region]) for l in range(n_legendre)],
                Sig2=csr_matrix(sig2[region]), chi=chi[region], eg=eg,
            )
            for region in range(n_coarse)
        }

        # One fresh effective material per coarse cell (mat_ids = 0..n-1);
        # the coarse geometry (edges / coord / BC) carries through, volumes
        # re-derived from the coarse edges.
        coarse_mesh = Mesh1D(
            edges=coarse_edges,
            mat_ids=np.arange(n_coarse, dtype=int),
            coord=coarse.coord,
            bc_left=coarse.bc_left,
            bc_right=coarse.bc_right,
        )
        return MaterialMesh(coarse_mesh, homogenized)

    # ── Comparison ──────────────────────────────────────────────────

    def compare(self, other: "Solution", *, rtol: float = 1e-12) -> "SolutionDiff":
        r"""Return a field-by-field difference summary against ``other``.

        Compares :attr:`keff` (when both have it) and the
        :math:`L^\infty`-norms of the flux deltas.  Useful for
        regression / refactor consistency checks.

        Parameters
        ----------
        other : Solution
            The reference solution to compare against.
        rtol : float
            Relative tolerance for the ``within_tolerance`` flag.

        Returns
        -------
        SolutionDiff
            Field-by-field summary.
        """
        if self.mesh is not other.mesh:
            raise ValueError(
                "Solution.compare: meshes differ — cross-mesh comparison "
                "is not defined under the typed-field contract."
            )

        if self.keff is not None and other.keff is not None:
            keff_abs = abs(self.keff - other.keff)
        else:
            keff_abs = None

        ang_diff = self.angular_flux.bulk.values - other.angular_flux.bulk.values
        sca_diff = self.scalar_flux.values - other.scalar_flux.values
        ang_linf = float(np.abs(ang_diff).max()) if ang_diff.size else 0.0
        sca_linf = float(np.abs(sca_diff).max()) if sca_diff.size else 0.0

        sca_norm = float(np.abs(other.scalar_flux.values).max())
        flux_ok = (sca_norm == 0.0) or (sca_linf <= rtol * sca_norm)
        keff_ok = (keff_abs is None) or (
            abs(other.keff or 1.0) > 0.0
            and keff_abs <= rtol * abs(other.keff or 1.0)
        )

        return SolutionDiff(
            keff_abs=keff_abs,
            angular_flux_linf=ang_linf,
            scalar_flux_linf=sca_linf,
            within_tolerance=bool(flux_ok and keff_ok),
        )


@dataclass(frozen=True)
class SolutionDiff:
    r"""Result of :meth:`Solution.compare`.

    Parameters
    ----------
    keff_abs : float or None
        Absolute eigenvalue difference :math:`|k_a - k_b|`, or
        ``None`` when at least one solution is fixed-source.
    angular_flux_linf : float
        :math:`L^\infty` norm of the angular-flux delta.
    scalar_flux_linf : float
        :math:`L^\infty` norm of the scalar-flux delta.
    within_tolerance : bool
        Aggregate verdict: True iff every available channel met the
        comparison rtol.
    """

    keff_abs: float | None
    angular_flux_linf: float
    scalar_flux_linf: float
    within_tolerance: bool
