r"""Continuous reference solution type for ORPHEUS verification.

This module defines the contract every verification reference in
ORPHEUS will commit to, starting in Phase 0 of the verification
campaign. A :class:`ContinuousReferenceSolution` is a
**mesh-independent, mathematically self-contained** solution to a
governing equation, derived by SymPy/mpmath from the transport,
diffusion, integral, or stochastic form that a target solver
discretises. It stands alone: you can evaluate it without running
any solver, without iteration, without a mesh, to arbitrary precision.

Vocabulary discipline
---------------------

Following Oberkampf & Roache, we use *verification* and *reference
solution* strictly:

- **Verification** — proving mathematical correctness of a code
  against a reference derived from the governing equation itself.
- **Reference solution** — an analytical or semi-analytical function
  of the independent variables, derived by pure mathematics, that
  any code claiming to solve the same equation must reproduce in
  the limit of zero discretisation error.
- **Benchmark** — a code-to-code comparison (L4 in the V&V ladder).
  **Never** used in this module. Collections like Sood et al. 1999
  and Ganapol 2008 use the word "benchmark" in their legacy titles;
  their *contents* are analytical reference solutions and we cite
  them as such.

See :doc:`/verification/reference_solutions` for the full treatment.

Operator-form taxonomy
----------------------

Every :class:`ContinuousReferenceSolution` commits to **one** operator
form. A reference solution of form ``"differential-sn"`` is a solution
to :math:`\mu\,\partial\psi/\partial x + \Sigma_t\psi = \mathrm{RHS}`
and has **no business** being consumed by a diffusion test, because
diffusion solves a different equation. Tests assert agreement between
the operator form of the reference and the equation form their solver
discretises.

Consumers
---------

:class:`ContinuousReferenceSolution` objects are produced by the
``derive_*`` functions inside ``orpheus.derivations.*`` modules and
registered in :mod:`orpheus.derivations.reference_values` under
``continuous_get(name)``. Tests consume them via
:meth:`ContinuousReferenceSolution.phi_on_mesh` to evaluate the
continuous reference on whatever discretisation the test chose.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Literal

import numpy as np


OperatorForm = Literal[
    "homogeneous",           # infinite medium; k from matrix eigenvalue
    "differential-sn",       # μ ∂ψ/∂x + Σ_t ψ = RHS (discrete ordinates)
    "differential-moc",      # dψ/ds + Σ_t ψ = Q along characteristics
    "diffusion",             # −∇·(D∇φ) + Σ_r φ = S
    "integral-peierls",      # φ = (c/2)∫ E₁ ⊛ φ + S (Peierls / CP)
    "stochastic-transport",  # integro-differential Boltzmann via sampling
]


BoundaryCondition = Literal[
    "vacuum",
    "reflective",
    "white",       # isotropic return — CP convention
    "marshak",     # diffusion mixed BC
    "periodic",
]


GeometryType = Literal[
    "homogeneous",   # infinite medium; no spatial variable
    "slab",          # 1D Cartesian [0, L]
    "sphere-1d",     # 1D radial on [0, R]
    "cylinder-1d",   # 1D radial on [0, R] with azimuthal symmetry
    "cartesian-2d",  # 2D Cartesian [0, Lx] × [0, Ly]
]


# ═══════════════════════════════════════════════════════════════════════
# Provenance and problem specification
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class Provenance:
    """Mathematical pedigree of a reference solution.

    Every :class:`ContinuousReferenceSolution` carries a ``Provenance``
    that tells the reader **where the mathematics comes from** —
    which textbook, which equation, which SymPy expression. This is
    the audit trail for the verification claim.

    Attributes
    ----------
    citation : str
        Primary literature reference with chapter/section/equation
        numbers, e.g. ``"Case & Zweifel 1967, Ch. 3 Eq. (3.3)"`` or
        ``"Sood, Forster & Parsons 1999 (LA-13511), Problem 1"``.
    derivation_notes : str
        Free-text exposition of the derivation steps: why this
        ansatz, what simplifying assumptions, what the residual is.
        Lives in the source module, surfaces via Sphinx.
    sympy_expression : str or None
        String form of the symbolic expression when the solution is
        a SymPy closed form. None for mpmath-quadrature-backed or
        transcendental-root-finder references.
    precision_digits : int or None
        Working precision of the mpmath evaluation. None means the
        reference is a closed-form analytical expression evaluated
        at IEEE double precision (~16 digits).
    """

    citation: str
    derivation_notes: str
    sympy_expression: str | None = None
    precision_digits: int | None = None


@dataclass(frozen=True)
class ProblemSpec:
    """Full mathematical specification of a verification problem.

    Carries everything a test needs to build its own mesh, materials,
    and solver invocation. The reference solution is evaluated against
    the same problem but without any discretisation.

    Attributes
    ----------
    materials : dict[int, Any]
        Material ID to :class:`orpheus.data.macro_xs.mixture.Mixture`
        mapping. Typed as ``Any`` to avoid a circular import; the
        contract is "anything SNSolver / CPSolver accepts."
    geometry_type : GeometryType
        One of the supported geometry tags; determines what
        ``geometry_params`` must contain.
    geometry_params : dict
        Geometry-specific dimensions. For ``"slab"``:
        ``{"length": float, "thicknesses": list, "mat_ids": list}``.
        For ``"sphere-1d"`` / ``"cylinder-1d"``:
        ``{"radius": float, "radii": list, "mat_ids": list}``.
        For ``"cartesian-2d"``: ``{"lx": float, "ly": float, ...}``.
    boundary_conditions : dict[str, BoundaryCondition]
        Face-labelled BCs. For 1D slab, keys are ``"left"``, ``"right"``.
        For 1D radial, key is ``"outer"`` (inner is always ``r=0`` symmetry).
    external_source : callable or None
        Fixed external source ``Q(x, g) -> array``. None for
        eigenvalue problems (k-eigenvalue fission source is computed
        by the solver itself, not the reference).
    is_eigenvalue : bool
        True if the problem is k-eigenvalue, False for fixed-source.
    n_groups : int
        Number of energy groups.
    """

    materials: dict[int, Any]
    geometry_type: GeometryType
    geometry_params: dict[str, Any]
    boundary_conditions: dict[str, BoundaryCondition]
    external_source: Callable[..., np.ndarray] | None = None
    is_eigenvalue: bool = True
    n_groups: int = 1


# ═══════════════════════════════════════════════════════════════════════
# Callable signatures for continuous fields
# ═══════════════════════════════════════════════════════════════════════

ScalarFluxFn = Callable[[np.ndarray, int], np.ndarray]
"""Signature ``(x: array, g: int) -> array`` returning scalar flux
at the requested points for energy group ``g``. For homogeneous
problems, ``x`` is ignored and the flat spectrum vector is returned
broadcast over the input shape."""

AngularFluxFn = Callable[[np.ndarray, np.ndarray, int], np.ndarray]
"""Signature ``(x: array, mu: array, g: int) -> array`` returning
angular flux. Only populated for operator forms where the reference
is inherently angle-resolved (differential-sn with an MMS ansatz that
keeps ψ_n angle-dependent). Most references use ``phi`` only and
leave ``psi=None``."""


# ═══════════════════════════════════════════════════════════════════════
# ContinuousReferenceSolution — the core type
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class ContinuousReferenceSolution:
    r"""A mesh-independent, mathematically self-contained reference solution.

    This is the central type of the ORPHEUS verification campaign.
    An instance is a **callable** representation of the exact solution
    to a governing equation at a specific problem configuration,
    derived by SymPy/mpmath from the equation itself. Solvers under
    test are verified by comparing their discretised output to this
    reference on a mesh of the test's choosing.

    Attributes
    ----------
    name : str
        Unique identifier, e.g. ``"sn_slab_2eg_2rg_continuous"``.
        Follows the convention ``<method>_<geometry>_<Neg>_<Nrg>[_<tag>]``.
    problem : ProblemSpec
        Full mathematical specification — the test uses this to
        build the matching mesh/solver invocation.
    operator_form : OperatorForm
        Which equation form the reference commits to. Tests assert
        the target solver discretises the same form.
    k_eff : float or None
        Eigenvalue for k-problems, None for fixed-source.
    phi : ScalarFluxFn
        Continuous scalar flux callable. Always defined.
    psi : AngularFluxFn or None
        Continuous angular flux. Only when the operator form carries
        angular resolution (differential-sn with an angle-dependent
        ansatz).
    provenance : Provenance
        Mathematical pedigree — citations, SymPy expression, precision.
    equation_labels : tuple[str, ...]
        Sphinx ``:label:`` IDs this reference exercises. Nexus builds
        test ↔ equation edges from these.
    vv_level : str or None
        ``"L1"`` for all reference-solution-driven tests, ``"L0"``
        for kernel-primitive identity checks.
    description : str
        Short human-readable summary shown in the V&V matrix.
    tolerance : str
        Expected agreement form, e.g. ``"< 1e-10"`` for eigenvalue,
        ``"O(h²)"`` for spatial convergence studies.

    Notes
    -----
    This type is **frozen** — once built, it cannot be mutated. Callable
    fields (``phi``, ``psi``) close over the SymPy/mpmath state that
    defines them, so they are self-contained after construction.

    A legacy :class:`~orpheus.derivations._types.VerificationCase`
    representation is available via :meth:`as_verification_case` for
    backward-compatible registry access during the migration.
    """

    name: str
    problem: ProblemSpec
    operator_form: OperatorForm
    phi: ScalarFluxFn
    provenance: Provenance
    k_eff: float | None = None
    psi: AngularFluxFn | None = None
    equation_labels: tuple[str, ...] = ()
    vv_level: str | None = "L1"
    description: str = ""
    tolerance: str = ""

    # ── Convenience: evaluate phi on a mesh ──────────────────────────

    def phi_on_mesh(self, mesh: Any, group: int = 0) -> np.ndarray:
        """Evaluate the reference scalar flux at every cell centre of ``mesh``.

        This is the hot-path call from tests: feed the test's own
        :class:`orpheus.geometry.Mesh1D` (or 2D variant) and get a
        cell-average (midpoint) evaluation. For operators with
        sharp boundary layers the cell-centred midpoint rule is only
        :math:`O(h^{2})`-accurate; higher-order tests should use
        :meth:`phi_cell_average` instead.

        Parameters
        ----------
        mesh : Mesh1D or Mesh2D
            Any object exposing a ``centers`` attribute or property.
        group : int
            Energy group index.

        Returns
        -------
        ndarray
            Scalar flux at each cell centre.
        """
        centers = np.asarray(mesh.centers, dtype=float)
        return np.asarray(self.phi(centers, group), dtype=float)

    def phi_cell_average(
        self, mesh: Any, group: int = 0, n_quad: int = 8,
    ) -> np.ndarray:
        r"""Evaluate the reference scalar flux as **cell averages**.

        For each cell :math:`[x_{i-1/2}, x_{i+1/2}]` returns
        :math:`\frac{1}{\Delta x_i}\int_{x_{i-1/2}}^{x_{i+1/2}}\phi(x)\,dx`
        via Gauss–Legendre quadrature of order ``n_quad``. Use this
        when the test compares against a finite-volume cell-average
        output (CP, diffusion FV) rather than a cell-centre value
        — otherwise the midpoint rule adds :math:`O(h^{2})` error
        on top of the solver's own discretisation error and can
        masquerade as a convergence failure.

        Parameters
        ----------
        mesh : Mesh1D
            1D mesh with ``edges`` attribute.
        group : int
            Energy group index.
        n_quad : int
            Gauss–Legendre quadrature order per cell. Default 8 is
            exact for polynomial integrands up to degree 15 — more
            than enough for smooth reference solutions.

        Returns
        -------
        ndarray
            Shape ``(n_cells,)``, the exact cell averages to GL
            quadrature precision.
        """
        edges = np.asarray(mesh.edges, dtype=float)
        widths = np.diff(edges)
        nodes, weights = np.polynomial.legendre.leggauss(n_quad)
        # Map GL nodes from [-1, 1] to each cell [x_{i-1/2}, x_{i+1/2}]
        half = 0.5 * widths[:, None]                    # (nc, 1)
        mids = 0.5 * (edges[:-1] + edges[1:])[:, None]  # (nc, 1)
        x_quad = mids + half * nodes[None, :]           # (nc, n_quad)
        phi_q = np.asarray(
            self.phi(x_quad.ravel(), group), dtype=float,
        ).reshape(x_quad.shape)
        # Cell average: (1/Δx) * (Δx/2) * Σ w_k phi(x_k) = 0.5 * Σ w_k phi
        return 0.5 * np.sum(weights[None, :] * phi_q, axis=1)

    # ── Backward-compat bridge ───────────────────────────────────────

    def as_verification_case(self) -> Any:
        """Return a legacy :class:`VerificationCase` pointing at this reference.

        Used during the migration so the existing registry and tests
        (which still key on ``VerificationCase``) keep working while
        downstream derivations are retrofitted one at a time.
        """
        from ._types import VerificationCase

        if self.k_eff is None:
            raise ValueError(
                f"Reference {self.name!r} has no k_eff; legacy "
                "VerificationCase bridge is eigenvalue-only."
            )

        return VerificationCase(
            name=self.name,
            k_inf=self.k_eff,
            method=_operator_form_to_legacy_method(self.operator_form),
            geometry=_geometry_to_legacy(self.problem.geometry_type),
            n_groups=self.problem.n_groups,
            n_regions=_infer_n_regions(self.problem),
            materials=self.problem.materials,
            geom_params=self.problem.geometry_params,
            latex=self.provenance.derivation_notes,
            description=self.description,
            tolerance=self.tolerance,
            vv_level=self.vv_level,
            equation_labels=self.equation_labels,
        )


# ═══════════════════════════════════════════════════════════════════════
# Internal helpers for the legacy bridge
# ═══════════════════════════════════════════════════════════════════════

def _operator_form_to_legacy_method(form: OperatorForm) -> str:
    """Map an operator form to the legacy ``method`` string used by
    :class:`VerificationCase`."""
    return {
        "homogeneous": "homo",
        "differential-sn": "sn",
        "differential-moc": "moc",
        "diffusion": "dif",
        "integral-peierls": "cp",
        "stochastic-transport": "mc",
    }[form]


def _geometry_to_legacy(geo: GeometryType) -> str:
    """Map a geometry tag to the legacy ``geometry`` string."""
    return {
        "homogeneous": "--",
        "slab": "slab",
        "sphere-1d": "sph1D",
        "cylinder-1d": "cyl1D",
        "cartesian-2d": "cart2D",
    }[geo]


def _infer_n_regions(problem: ProblemSpec) -> int:
    """Best-effort region count from the problem geometry params."""
    gp = problem.geometry_params
    if "mat_ids" in gp:
        return len(gp["mat_ids"])
    if "thicknesses" in gp:
        return len(gp["thicknesses"])
    if "radii" in gp:
        return len(gp["radii"])
    return 1
