r"""Mathematical billiards — phase-space dynamics with boundary reflection.

This module exposes :class:`Billiard`, the math-rich frame that names
the structure shared by every Variant α reference solver in
:mod:`~orpheus.derivations.continuous.trajectory_resolvent`. The
implementation is a thin facade over the geometry-specific
``solve_greens_function_*`` entry points; what it adds is a single
named *concept* (the billiard system) that covers all six geometries
(slab, sphere, cylinder, slab-asymmetric, hollow sphere, annulus) at
both orbit-space classes (one-surface compact and two-surface).

Mathematical billiards in one paragraph
---------------------------------------

A *billiard* in the Birkhoff–Sinai sense is a dynamical system on a
domain :math:`M \subset \mathbb R^d` (the "billiard table") with a
specified reflection law on :math:`\partial M`. A particle moves
freely between collisions and reflects at the boundary. The flow is
deterministic, time-reversible, volume-preserving (Liouville), and
possibly ergodic. When :math:`M` is bounded and the reflection is
specular (angle of incidence equals angle of reflection), the
billiard is *closed*; when the reflection has albedo :math:`\alpha
\in [0, 1)`, the billiard is *absorbing* and orbits decay
geometrically with rate :math:`\alpha`. Linear (neutron) transport in
a homogeneous geometry IS a mathematical billiard with the
*streaming* flow on the interior and a *partial-albedo* reflection
law on :math:`\partial M`. The Variant α resolvent
:math:`T = (I - S)^{-1}` is the Birkhoff transfer operator's
resolvent for that billiard, summed in closed form via the Neumann
series :math:`\sum_n (\alpha\, e^{-\tau})^n`.

Why "Billiard"
--------------

The class name encodes a piece of mathematics every reader benefits
from learning: that **the rank-1 / rank-2 closure formula in
:mod:`.variant_alpha_core` is not a transport-physics special trick
but the analytic resolvent of a Birkhoff transfer operator on a
billiard table**. The same operator algebra appears in classical
mechanics (gas-of-particles in a box), in semiclassical quantum
chaos (Gutzwiller trace formula on hyperbolic billiards), and in
ergodic theory (Sinai's dispersing billiards). The transport-physics
specialization we ship here is one application of a much wider
mathematical theory; naming the class :class:`Billiard` makes that
visible at every call site.

The full theory page section
:ref:`billiards-and-the-variant-alpha-resolvent` derives the
billiard ↔ Variant α correspondence formally. The class docstring
below is a compressed walk-through so the reader who lands on the
class definition first still gets the load-bearing intuition.

Architectural role
------------------

:class:`Billiard` is to :mod:`~orpheus.derivations.continuous.trajectory_resolvent`
what :class:`MomentSpace` (Galerkin half-range projection) is to
:mod:`~orpheus.derivations.continuous.fn_method`, what
:class:`Spectrum` (Case singular eigenfunction expansion) will be to
:mod:`~orpheus.derivations.continuous.singular_eigenfunction`, and
what ``CPMesh`` is to the production collision-probability solver.
Each class is a method-specific computational specialization of the
abstract :class:`~orpheus.derivations.common.geometry_spec.GeometrySpec`.
The unifying picture lives at
:ref:`billiards-and-the-variant-alpha-resolvent`.

R2 hindsight refactor
---------------------

This module lands the second hindsight refactor (R2 in
``.claude/plans/trajectory_resolvent_hindsight_refactor.md``) for the
Variant α family. R1 already collapsed the 11 byte-identical
power-iteration outer loops into
:func:`~orpheus.derivations.continuous.trajectory_resolvent.power_iteration.power_iterate_variant_alpha`.
R2 unifies the 12 per-geometry result dataclasses behind the SHARED
cross-method types
:class:`~orpheus.derivations.common.solution_types.CriticalSolution`
/
:class:`~orpheus.derivations.common.solution_types.FluxSolution`
(introduced by the parallel ``MomentSpace`` agent) and exposes the
billiard concept as a first-class object the user can construct,
introspect, and dispatch on.

Bit-equality with R1 is preserved exactly: every
``solve_critical()`` call delegates to the existing
``solve_greens_function_*`` entry point and re-packs its result into
a unified dataclass. No FP arithmetic moves; no quadrature changes;
no closures rewire. The 84-test cross-method regression net + the
205-test trajectory_resolvent suite must agree to IEEE-754 exact
between pre- and post-refactor.

References
----------

- Birkhoff, G.D. (1927). *Dynamical Systems*, Ch. 6 (billiard
  flows). Reissued by AMS Colloquium Publications, 1966.
- Sinai, Ya.G. (1970). "Dynamical systems with elastic reflections
  and their applications." *Russ. Math. Surveys* **25**, 137–189.
  DOI: 10.1070/RM1970v025n02ABEH003794.
- Bunimovich, L.A. (1979). "On the ergodic properties of nowhere
  dispersing billiards." *Commun. Math. Phys.* **65**, 295–312.
- Chernov, N. and Markarian, R. (2006). *Chaotic Billiards*, AMS
  Mathematical Surveys and Monographs **127**.
- Sanchez, R. (2002). "Treatment of Boundary Conditions in
  Trajectory-Based Deterministic Transport Methods." *Nucl. Sci.
  Eng.* **140**, 23–50. The transport-physics specialization of
  billiard theory we implement here. The Variant α resolvent
  :math:`T = (I - S)^{-1}` is Sanchez's Eq. (A4) in our notation.
- Sanchez, R. (1986). *Transp. Theor. Stat. Phys.* **14**.
  DOI: 10.1080/00411458608210456. The :math:`\omega_1 = 0`
  (isotropic-scattering) precursor.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.geometry_spec import GeometrySpec

# Use the SHARED cross-method result types — these define the contract
# every math-heart class (MomentSpace for fn_method, Billiard for
# trajectory_resolvent) populates. Method-specific rich results live
# under ``metadata["raw_result"]`` (the legacy result dataclass) and
# ``metadata["psi"]``, ``metadata["phi"]``, ``metadata["iterations"]``,
# ``metadata["mesh"]`` for direct field access.
from orpheus.derivations.common.solution_types import (
    CriticalSolution,
    FluxSolution,
)


# Re-export so consumers can do
#   from orpheus.derivations.continuous.trajectory_resolvent import (
#       Billiard, CriticalSolution, FluxSolution,
#   )
# even though the canonical home is ``derivations.common.solution_types``.
__all__ = ["Billiard", "CriticalSolution", "FluxSolution"]


# ─────────────────────────────────────────────────────────────────────
# Result-shape contract
# ─────────────────────────────────────────────────────────────────────
#
# Per the schema in ``derivations.common.solution_types``, the unified
# CriticalSolution is small (eigenvalue, eigenvalue_kind,
# parameter_value, parameter_kind, converged, metadata). Method-specific
# rich data — the angular flux ``psi``, scalar flux ``phi``, iteration
# count, and mesh metadata — lives in ``metadata`` under the canonical
# keys consumed by the cross-method protocol
# (``tests.cross_method.adapters``):
#
# - ``metadata["raw_result"]`` — the legacy per-geometry result
#   dataclass (e.g. ``GreensFunctionResult``); preserves bit-equal
#   access to every original field.
# - ``metadata["psi"]`` / ``metadata["phi"]`` — direct convenience
#   handles to the angular and scalar flux arrays.
# - ``metadata["iterations"]`` — power-iteration step count.
# - ``metadata["mesh"]`` — dict of grid arrays (``r_nodes``, ``x_nodes``,
#   ``mu_nodes``, ``mu_axial_nodes``, ``phi_az_nodes``,
#   ``region_at_node``).
# - ``metadata["n_groups"]`` — group count (1 for 1G; G for MG / MR).
# - ``metadata["geometry_kind"]`` — the Billiard's geometry tag, for
#   downstream dispatch.
#
# For Variant α power iteration, ``parameter_kind = "fixed_geometry"``:
# the geometry is fixed by the caller (we don't root-find a critical
# parameter), so the parameter_value field carries the geometry's
# characteristic length (``R`` for sphere/cylinder, ``L`` for slab,
# ``R_out`` for hollow geometries). Cross-method comparators that
# need the eigenvalue at this configuration read it from
# ``solution.eigenvalue`` directly.


# ─────────────────────────────────────────────────────────────────────
# Billiard system
# ─────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class Billiard:
    r"""A mathematical billiard — phase-space dynamics with boundary
    reflection.

    A *billiard* in the Birkhoff–Sinai sense is a dynamical system
    consisting of:

    1. **A domain :math:`M`** (the "billiard table") — here the
       spatial geometry (slab, sphere, cylinder, hollow sphere,
       annulus, asymmetric slab; classified by the orbit-space M/G
       structure — see
       :ref:`orbit-space-m-g-classification` for the precise
       taxonomy).
    2. **A boundary reflection law** at :math:`\partial M` — here
       the BC parametrized by :math:`\alpha`:

       - specular (:math:`\alpha = 1`) ↔ pure billiard, ergodic at
         high density;
       - vacuum (:math:`\alpha = 0`) ↔ open billiard, particles
         escape;
       - partial-albedo (:math:`\alpha \in (0, 1)`) ↔ scattering
         billiard (energy / particle loss at boundary).
    3. **Free flight inside the domain**, governed by the streaming
       operator :math:`\Omega \cdot \nabla + \Sigma_t` along
       characteristic curves (rays / geodesics in the appropriate
       metric).
    4. **A transfer operator** :math:`S` that maps boundary-flux
       moments to boundary-flux moments under one full bounce
       period, and its resolvent :math:`T = (I - S)^{-1}`.

    The closed-form multi-bounce summation
    :math:`\sum_n \alpha^n e^{-n\tau}` that Variant α evaluates IS
    the Birkhoff transfer-operator resolvent for the billiard with
    absorption rate :math:`\alpha` and chord-length :math:`\tau`.
    When :math:`\alpha = 1` (specular boundary) and :math:`\tau
    \to 0` (high-density medium) the billiard becomes ergodic and
    :math:`T` diverges; when :math:`\alpha < 1` (absorbing
    boundary), :math:`S` is a contraction (:math:`\|S\| < 1`) and
    the Neumann series for :math:`(I - S)^{-1}` converges
    geometrically with rate :math:`\alpha`.

    The orbit-space M/G classification (slab / sphere / cylinder /
    hollow_sphere / annulus / slab_asymmetric) determines the
    **rank** of the closure operator: rank-1 closure for
    one-endpoint orbit spaces (closed sphere, closed cylinder),
    rank-2 closure for two-endpoint orbit spaces (slab,
    slab_asymmetric, hollow_sphere, annulus). The rank is the
    bond dimension of the open MPO on the bounce-event lattice
    (see the theory page
    :ref:`billiards-and-the-variant-alpha-resolvent`).

    The billiard ↔ resolvent correspondence in two lines
    --------------------------------------------------------

    A trajectory in the billiard is an alternating sequence of
    *streaming segments* (free flight along
    :math:`\Omega \cdot \nabla \psi + \Sigma_t \psi = q`) and
    *boundary reflections* (multiplication by :math:`\alpha`). The
    *transfer operator* :math:`S` advances the trajectory by ONE
    bounce period: streaming attenuation :math:`e^{-\tau_{\rm
    period}}` followed by reflection :math:`\alpha`. The integrated
    multi-bounce contribution at the entry point is the geometric
    series

    .. math::

        \sum_{n=0}^{\infty} S^n = (I - S)^{-1} = T,

    and the closed-form sum :math:`T = 1/(1 - \alpha\,e^{-\tau})`
    is the reason Variant α can replace an iterative bounce
    discretization with a single algebraic step. The full derivation
    is at :ref:`billiards-and-the-variant-alpha-resolvent`.

    Construction
    ------------

    Use the :meth:`Billiard.from_problem` factory. Direct
    construction is supported but not recommended — the factory
    validates the materials/geometry combination and resolves the
    closure rank from the orbit-space class.

    The class is *frozen* (immutable) so that a single billiard
    instance can be safely cached and reused across multiple solves
    with different ``alpha`` overrides via
    :meth:`Billiard.with_alpha`.

    Cross-method analog
    -------------------

    :class:`Billiard` is to ``trajectory_resolvent`` what:

    - ``MomentSpace`` (Galerkin half-range projection on Legendre
      moments) is to
      :mod:`~orpheus.derivations.continuous.fn_method`. The fn_method
      *projects* the Peierls integral equation onto Legendre moments
      and solves a finite-dimensional moment system. The billiard's
      discrete transfer operator is morally the same machinery seen
      from the spatial side.
    - ``Spectrum`` (Case singular eigenfunction expansion) will be
      to :mod:`~orpheus.derivations.continuous.singular_eigenfunction`.
      The Case method *expands* the angular flux on the operator's
      eigenfunctions; the billiard's resolvent expansion
      :math:`\sum_n S^n` is the time-domain counterpart of Case's
      spectral expansion.
    - ``CPMesh`` is to the production collision-probability solver
      what :class:`Billiard` is to the reference family.

    These siblings are all method-specific computational
    specializations of the abstract
    :class:`~orpheus.derivations.common.geometry_spec.GeometrySpec`.
    See :ref:`billiards-and-the-variant-alpha-resolvent` for the
    framing.

    Attributes
    ----------
    geometry_kind : str
        One of ``"sphere"``, ``"cylinder"``, ``"slab"``,
        ``"slab_asymmetric"``. Selected internally from the
        :class:`GeometrySpec` and the alpha-dict shape; determines
        which underlying solver is dispatched. (The internal
        dispatcher additionally supports ``"sphere_mr"``,
        ``"hollow_sphere"``, ``"annulus"`` for legacy callers, but
        :meth:`from_problem` cannot construct those today — they
        await Step 3's multi-region :class:`GeometrySpec` extension.)
    materials : dict[int, Mixture]
        Production-protocol cross-section payload, keyed by
        material ID. The Protocol-conformant view of the same
        cross sections that the underlying solver consumes as a
        flat numpy-array dict.
    geometry_spec : GeometrySpec
        Method-agnostic geometry + boundary specification. The
        :class:`Billiard` derives :attr:`geometry_kind` and
        :attr:`alpha_payload` from this spec plus the factory's
        :paramref:`alpha` argument.
    xs_payload : dict
        Internal solver-facing cross-section payload (numpy arrays
        / scalars). Format depends on :attr:`geometry_kind`:

        - 1G case: ``{"sigma_t": float, "sigma_s": float,
          "nu_sigma_f": float}``.
        - MG case: ``{"sigma_t": (G,), "sigma_s": (G,G),
          "nu_sigma_f": (G,), "chi": (G,) optional}``.

        Synthesized from :attr:`materials` by
        :func:`_mixture_to_solver_xs_payload`. End users should
        prefer :attr:`materials` for the Protocol view.
    geometry_payload : dict
        Internal solver-facing geometry kwargs that pin the domain
        :math:`M`:

        - ``sphere``: ``{"R": float}``.
        - ``cylinder``: ``{"R": float}``.
        - ``slab``: ``{"L": float}``.
        - ``slab_asymmetric``: ``{"L": float}``.

        Synthesized from :attr:`geometry_spec` by
        :func:`_geometry_payload_for_solver`. End users should
        prefer :attr:`geometry_spec` for the Protocol view.
    alpha_payload : dict
        Reflection law on :math:`\partial M`:

        - One-surface compact (sphere, cylinder, slab[symmetric]):
          ``{"alpha": float}``.
        - Two-surface (slab_asymmetric): ``{"alpha_left": float,
          "alpha_right": float}``.
    quadrature : dict
        Quadrature orders. Standard keys: ``n_r`` / ``n_x``,
        ``n_mu`` / ``n_mu_axial``, ``n_phi_az``, ``n_traj_quad``.
        Defaults are inherited from each geometry's solver if the
        key is absent.
    closure_rank : int
        ``1`` for one-endpoint orbit spaces, ``2`` for two-endpoint
        orbit spaces. Auto-detected from :attr:`geometry_kind` by
        :meth:`from_problem`.

    Examples
    --------
    >>> import numpy as np
    >>> from orpheus.data.macro_xs.mixture import Mixture  # doctest: +SKIP
    >>> from orpheus.derivations.common.geometry_spec import GeometrySpec
    >>> from orpheus.derivations.continuous.trajectory_resolvent import (
    ...     Billiard,
    ... )
    >>> from orpheus.geometry.mesh import BC

    Closed homogeneous sphere with k_eff = k_inf:

    >>> spec = GeometrySpec(  # doctest: +SKIP
    ...     geometry="sphere",
    ...     critical_dimension_mfp=2.5,
    ...     critical_dimension_cm=5.0,
    ...     n_groups=1,
    ...     bc_left=BC.reflective,
    ...     bc_right=BC.reflective,
    ... )
    >>> b = Billiard.from_problem(  # doctest: +SKIP
    ...     materials={0: pu_mixture},
    ...     geometry_spec=spec,
    ...     alpha=1.0,
    ...     quadrature={"n_r": 24, "n_mu": 24, "n_traj_quad": 64},
    ... )
    >>> sol = b.solve_critical()  # doctest: +SKIP
    >>> sol.eigenvalue            # doctest: +SKIP
    1.0  # k_inf = nu_sigma_f / (sigma_t - sigma_s) = 0.1 / 0.1

    Asymmetric slab with independent wall reflectivities (the
    ``slab_asymmetric`` family is selected automatically when
    :paramref:`alpha` is a dict containing ``alpha_left`` /
    ``alpha_right``):

    >>> slab_spec = GeometrySpec(  # doctest: +SKIP
    ...     geometry="slab",
    ...     critical_dimension_mfp=1.25,
    ...     critical_dimension_cm=2.5,
    ...     n_groups=1,
    ...     bc_left=BC.vacuum,
    ...     bc_right=BC.vacuum,
    ... )
    >>> b = Billiard.from_problem(  # doctest: +SKIP
    ...     materials={0: pu_mixture},
    ...     geometry_spec=slab_spec,
    ...     alpha={"alpha_left": 0.7, "alpha_right": 0.9},
    ... )
    >>> sol = b.solve_critical()  # doctest: +SKIP
    """

    geometry_kind: str
    xs_payload: dict[str, Any]
    geometry_payload: dict[str, Any]
    alpha_payload: dict[str, float]
    quadrature: dict[str, int]
    closure_rank: int
    materials: dict[int, Mixture] = field(default_factory=dict)
    geometry_spec: GeometrySpec | None = None
    method_name: str = "trajectory_resolvent"

    @classmethod
    def from_problem(
        cls,
        *,
        materials: dict[int, Mixture],
        geometry_spec: GeometrySpec,
        alpha: float | dict[str, float] = 1.0,
        quadrature: dict[str, int] | None = None,
    ) -> "Billiard":
        r"""Construct a billiard from problem-statement inputs.

        This factory consumes the production-protocol pair
        ``(materials: dict[int, Mixture], geometry_spec: GeometrySpec)``
        — the same shape every TransportSolver-conformant class
        accepts — and:

        1. infers the internal :attr:`geometry_kind` from the spec
           (and the alpha-dict shape, for the slab-asymmetric
           branch);
        2. builds the solver-facing :attr:`geometry_payload` and
           :attr:`xs_payload` from the spec / mixture pair;
        3. normalizes :paramref:`alpha` into the per-geometry
           :attr:`alpha_payload`;
        4. resolves :attr:`closure_rank` from the orbit-space class.

        Parameters
        ----------
        materials : dict[int, Mixture]
            Production-protocol cross sections, keyed by material
            ID. Today's Billiard factory consumes only the active
            material at key ``0``; multi-region (``sphere_mr``)
            support is gated on Step 3's multi-region GeometrySpec
            extension.
        geometry_spec : GeometrySpec
            Method-agnostic geometry + boundary specification.
            Supported families: ``"slab"``, ``"sphere"``,
            ``"cylinder"``. The ``slab_asymmetric`` family is
            selected automatically when :paramref:`alpha` is a
            dict carrying ``alpha_left`` / ``alpha_right`` keys.
        alpha : float or dict, default 1.0
            Boundary reflectivity. Float: applied symmetrically
            (one-endpoint billiards) or as ``alpha_left =
            alpha_right`` (two-endpoint billiards). Dict: passed
            through verbatim. For asymmetric slabs use
            ``{"alpha_left": ..., "alpha_right": ...}``.
        quadrature : dict, optional
            Quadrature orders. Keys depend on geometry. Defaults
            are inherited from each geometry's solver.

        Returns
        -------
        Billiard
            A frozen billiard instance ready for
            :meth:`solve_critical` / :meth:`solve_fixed_source`.
        """
        if not _is_mixture_dict(materials):
            raise ValueError(
                "Billiard.from_problem: materials must be a "
                "dict[int, Mixture]; got keys "
                f"{list((materials or {}).keys())!r}."
            )

        geometry_kind = _infer_geometry_kind_from_spec(
            geometry_spec, alpha
        )
        geometry_payload = _geometry_payload_for_solver(
            geometry_kind, geometry_spec
        )
        xs_payload = _mixture_to_solver_xs_payload(
            materials, geometry_kind
        )

        # Normalize alpha into the per-geometry payload.
        if isinstance(alpha, dict):
            alpha_payload: dict[str, float] = dict(alpha)
        elif geometry_kind == "slab_asymmetric":
            alpha_payload = {
                "alpha_left": float(alpha),
                "alpha_right": float(alpha),
            }
        else:
            alpha_payload = {"alpha": float(alpha)}

        # Resolve closure rank from orbit-space class.
        closure_rank = 2 if geometry_kind == "slab_asymmetric" else 1

        return cls(
            geometry_kind=geometry_kind,
            xs_payload=xs_payload,
            geometry_payload=geometry_payload,
            alpha_payload=alpha_payload,
            quadrature=dict(quadrature or {}),
            closure_rank=closure_rank,
            materials=dict(materials),
            geometry_spec=geometry_spec,
            method_name="trajectory_resolvent",
        )

    def with_alpha(self, alpha: float | dict[str, float]) -> "Billiard":
        """Return a copy with a different boundary reflectivity."""
        return Billiard.from_problem(
            materials=self.materials,
            geometry_spec=self.geometry_spec,
            alpha=alpha,
            quadrature=self.quadrature,
        )

    # ─────────────────────────────────────────────────────────────────
    # solve methods — thin facades over the existing solvers
    # ─────────────────────────────────────────────────────────────────

    def solve_critical(
        self,
        *,
        max_iter: int | None = None,
        tol: float | None = None,
        initial_psi: np.ndarray | None = None,
        initial_k: float | None = None,
    ) -> CriticalSolution:
        r"""Solve the *k*-eigenproblem on this billiard.

        Dispatches to the geometry-specific
        ``solve_greens_function_*`` entry point and re-packs the
        result into the SHARED cross-method
        :class:`~orpheus.derivations.common.solution_types.CriticalSolution`.
        Bit-equal with the underlying entry point's return value
        (every FP operation runs in the same order).

        Schema
        ------

        The returned :class:`CriticalSolution` carries:

        - ``eigenvalue`` — the converged :math:`k_{\rm eff}`.
        - ``eigenvalue_kind`` — always ``"k_eff"``.
        - ``parameter_value`` — the geometry's characteristic length
          (:math:`R` for sphere/cylinder/sphere_mr, :math:`L` for
          slab/slab_asymmetric, :math:`R_{\rm out}` for hollow
          geometries).
        - ``parameter_kind`` — always ``"fixed_geometry"`` (Variant α
          reports :math:`k_{\rm eff}` at a fixed configuration; no
          critical root-find).
        - ``converged`` — whether the power iteration met its
          tolerance.
        - ``metadata`` — a dict with the rich method-specific data:

          * ``"raw_result"`` — the legacy per-geometry result
            dataclass (preserved bit-for-bit).
          * ``"psi"`` — angular flux array.
          * ``"phi"`` — scalar flux array.
          * ``"iterations"`` — power-iteration step count.
          * ``"mesh"`` — dict of grid arrays (``r_nodes`` /
            ``x_nodes``, ``mu_nodes``, ``mu_axial_nodes``,
            ``phi_az_nodes``, ``region_at_node`` as applicable).
          * ``"n_groups"`` — group count.
          * ``"geometry_kind"`` — the Billiard's geometry tag.

        Parameters
        ----------
        max_iter : int, optional
            Maximum power-iteration count. Defaults to the
            underlying solver's default (200 for 1G, 300 for MG,
            500 for MR).
        tol : float, optional
            Relative-:math:`k` convergence tolerance. Defaults to
            ``1e-10`` for 1G, ``1e-9`` for MG / MR.
        initial_psi : ndarray, optional
            Initial angular flux iterate. Shape must match the
            geometry's phase-space discretization.
        initial_k : float, optional
            Initial :math:`k_{\rm eff}` estimate. Defaults to
            ``k_inf`` from the homogenized infinite-medium balance.

        Returns
        -------
        CriticalSolution
            The shared cross-method solution type — read
            ``solution.eigenvalue`` for :math:`k_{\rm eff}` and
            ``solution.metadata["psi"]`` /
            ``solution.metadata["phi"]`` for the flux fields.
        """
        return _dispatch_critical(
            self,
            max_iter=max_iter,
            tol=tol,
            initial_psi=initial_psi,
            initial_k=initial_k,
        )

    def solve_fixed_source(
        self,
        external_source: np.ndarray,
        *,
        max_iter: int | None = None,
        tol: float | None = None,
    ) -> FluxSolution:
        r"""Solve the fixed-source problem on this billiard.

        Currently only supported for :attr:`geometry_kind` ==
        ``"sphere_mr"`` (the Garcia 2021 multi-region sphere
        benchmark family). Other geometries raise
        :class:`NotImplementedError`. A future extension will mirror
        the per-geometry power-iteration unification with a fixed-
        source iteration driver.

        Parameters
        ----------
        external_source : ndarray
            Per-region per-group external source
            :math:`Q^{\rm ext}` (per unit volume per second, NOT
            per steradian — the operator divides by :math:`4\pi`
            internally). Shape ``(n_regions, G)`` for
            :paramref:`geometry_kind` ``"sphere_mr"``.
        max_iter : int, optional
            Maximum fixed-source iterations. Default 500.
        tol : float, optional
            Relative scalar-flux convergence tolerance. Default
            ``1e-8``.

        Returns
        -------
        FluxSolution
            The shared cross-method solution type — read
            ``solution.metadata["psi"]`` / ``solution.metadata["phi"]``
            for the full flux fields and ``solution.scalar_flux``
            for the leading-group scalar flux profile.
        """
        return _dispatch_fixed_source(
            self,
            external_source=external_source,
            max_iter=max_iter,
            tol=tol,
        )


# ─────────────────────────────────────────────────────────────────────
# Internal dispatchers — bridge from Billiard to the per-geometry solvers
# ─────────────────────────────────────────────────────────────────────


def _filter_kwargs(d: dict[str, Any], allowed: tuple[str, ...]) -> dict[str, Any]:
    """Return only the allowed keys (drops unspecified defaults)."""
    return {k: d[k] for k in allowed if k in d}


def _dispatch_critical(
    billiard: Billiard,
    *,
    max_iter: int | None,
    tol: float | None,
    initial_psi: np.ndarray | None,
    initial_k: float | None,
) -> CriticalSolution:
    """Route to the appropriate ``solve_greens_function_*`` entry."""
    # Lazy imports — keep Billiard importable from doc-build contexts
    # that may not yet have all geometry modules wired.
    from orpheus.derivations.continuous.trajectory_resolvent import (
        greens_function as gf_sphere,
        greens_function_cylinder as gf_cyl,
        greens_function_slab as gf_slab,
        greens_function_slab_asymmetric as gf_slab_asym,
        greens_function_hollow_sphere as gf_hollow,
        greens_function_annulus as gf_annulus,
    )

    g = billiard.geometry_kind
    mats = billiard.xs_payload
    geom = billiard.geometry_payload
    a = billiard.alpha_payload
    q = billiard.quadrature

    # Build common iteration kwargs (skip None — let the underlying
    # solver use its own defaults).
    iter_kwargs: dict[str, Any] = {}
    if max_iter is not None:
        iter_kwargs["max_iter"] = max_iter
    if tol is not None:
        iter_kwargs["tol"] = tol
    if initial_psi is not None:
        iter_kwargs["initial_psi"] = initial_psi

    # Determine multi-group from materials shape.
    sigma_t = np.asarray(mats["sigma_t"])
    is_mg = sigma_t.ndim >= 1 and sigma_t.size > 1

    # Resolve the configuration parameter the eigenvalue is associated
    # with. Variant α is "fixed_geometry" — the user picks the geometry
    # and we report k_eff at that fixed config (no critical root-find).
    parameter_kind = "fixed_geometry"
    if g in ("sphere", "cylinder"):
        parameter_value = float(geom["R"])
    elif g in ("slab", "slab_asymmetric"):
        parameter_value = float(geom["L"])
    elif g in ("hollow_sphere", "annulus"):
        parameter_value = float(geom["R_out"])
    elif g == "sphere_mr":
        parameter_value = float(np.asarray(geom["radii"])[-1])
    else:
        parameter_value = 0.0
    wrap_extra = dict(
        parameter_value=parameter_value,
        parameter_kind=parameter_kind,
        geometry_kind=g,
    )

    if g == "sphere":
        if is_mg:
            res = gf_sphere.solve_greens_function_sphere_mg(
                R=geom["R"],
                sigma_t=mats["sigma_t"],
                sigma_s=mats["sigma_s"],
                nu_sigma_f=mats["nu_sigma_f"],
                chi=mats.get("chi"),
                alpha=a["alpha"],
                **_filter_kwargs(q, ("n_r", "n_mu", "n_traj_quad")),
                **iter_kwargs,
                **({"initial_k": initial_k} if initial_k is not None else {}),
            )
            return _wrap_mg(res, n_groups=int(sigma_t.size), **wrap_extra)
        else:
            res = gf_sphere.solve_greens_function_sphere(
                R=geom["R"],
                sigma_t=float(sigma_t),
                sigma_s=float(np.asarray(mats["sigma_s"])),
                nu_sigma_f=float(np.asarray(mats["nu_sigma_f"])),
                alpha=a["alpha"],
                **_filter_kwargs(q, ("n_r", "n_mu", "n_traj_quad")),
                **iter_kwargs,
            )
            return _wrap_1g(res, n_groups=1, **wrap_extra)

    if g == "sphere_mr":
        res = gf_sphere.solve_greens_function_sphere_mr(
            radii=geom["radii"],
            sigma_t=mats["sigma_t"],
            sigma_s=mats["sigma_s"],
            nu_sigma_f=mats["nu_sigma_f"],
            chi=mats.get("chi"),
            alpha=a["alpha"],
            **_filter_kwargs(q, ("n_r", "n_mu", "n_traj_quad")),
            **iter_kwargs,
            **({"initial_k": initial_k} if initial_k is not None else {}),
        )
        # Multi-region always carries (G, n_r, n_mu) shape — psi_g.
        n_groups = int(np.atleast_2d(np.asarray(mats["sigma_t"])).shape[-1]
                       if np.asarray(mats["sigma_t"]).ndim > 1
                       else 1)
        return _wrap_mr(res, n_groups=n_groups, **wrap_extra)

    if g == "cylinder":
        if is_mg:
            res = gf_cyl.solve_greens_function_cylinder_mg(
                R=geom["R"],
                sigma_t=mats["sigma_t"],
                sigma_s=mats["sigma_s"],
                nu_sigma_f=mats["nu_sigma_f"],
                chi=mats.get("chi"),
                alpha=a["alpha"],
                **_filter_kwargs(q, (
                    "n_r", "n_mu_axial", "n_phi_az", "n_traj_quad",
                )),
                **iter_kwargs,
                **({"initial_k": initial_k} if initial_k is not None else {}),
            )
            return _wrap_cyl_mg(res, n_groups=int(sigma_t.size), **wrap_extra)
        else:
            res = gf_cyl.solve_greens_function_cylinder(
                R=geom["R"],
                sigma_t=float(sigma_t),
                sigma_s=float(np.asarray(mats["sigma_s"])),
                nu_sigma_f=float(np.asarray(mats["nu_sigma_f"])),
                alpha=a["alpha"],
                **_filter_kwargs(q, (
                    "n_r", "n_mu_axial", "n_phi_az", "n_traj_quad",
                )),
                **iter_kwargs,
            )
            return _wrap_cyl_1g(res, n_groups=1, **wrap_extra)

    if g == "slab":
        if is_mg:
            res = gf_slab.solve_greens_function_slab_mg(
                L=geom["L"],
                sigma_t=mats["sigma_t"],
                sigma_s=mats["sigma_s"],
                nu_sigma_f=mats["nu_sigma_f"],
                chi=mats.get("chi"),
                alpha=a["alpha"],
                **_filter_kwargs(q, ("n_x", "n_mu", "n_traj_quad")),
                **iter_kwargs,
                **({"initial_k": initial_k} if initial_k is not None else {}),
            )
            return _wrap_slab_mg(res, n_groups=int(sigma_t.size), **wrap_extra)
        else:
            res = gf_slab.solve_greens_function_slab(
                L=geom["L"],
                sigma_t=float(sigma_t),
                sigma_s=float(np.asarray(mats["sigma_s"])),
                nu_sigma_f=float(np.asarray(mats["nu_sigma_f"])),
                alpha=a["alpha"],
                **_filter_kwargs(q, ("n_x", "n_mu", "n_traj_quad")),
                **iter_kwargs,
            )
            return _wrap_slab_1g(res, n_groups=1, **wrap_extra)

    if g == "slab_asymmetric":
        if is_mg:
            res = gf_slab_asym.solve_greens_function_slab_asymmetric_mg(
                L=geom["L"],
                sigma_t=mats["sigma_t"],
                sigma_s=mats["sigma_s"],
                nu_sigma_f=mats["nu_sigma_f"],
                chi=mats.get("chi"),
                alpha_left=a["alpha_left"],
                alpha_right=a["alpha_right"],
                **_filter_kwargs(q, ("n_x", "n_mu", "n_traj_quad")),
                **iter_kwargs,
                **({"initial_k": initial_k} if initial_k is not None else {}),
            )
            return _wrap_slab_mg(res, n_groups=int(sigma_t.size), **wrap_extra)
        else:
            res = gf_slab_asym.solve_greens_function_slab_asymmetric(
                L=geom["L"],
                sigma_t=float(sigma_t),
                sigma_s=float(np.asarray(mats["sigma_s"])),
                nu_sigma_f=float(np.asarray(mats["nu_sigma_f"])),
                alpha_left=a["alpha_left"],
                alpha_right=a["alpha_right"],
                **_filter_kwargs(q, ("n_x", "n_mu", "n_traj_quad")),
                **iter_kwargs,
            )
            return _wrap_slab_1g(res, n_groups=1, **wrap_extra)

    if g == "hollow_sphere":
        if is_mg:
            res = gf_hollow.solve_greens_function_hollow_sphere_mg(
                R_in=geom["R_in"],
                R_out=geom["R_out"],
                sigma_t=mats["sigma_t"],
                sigma_s=mats["sigma_s"],
                nu_sigma_f=mats["nu_sigma_f"],
                chi=mats.get("chi"),
                alpha_in=a["alpha_in"],
                alpha_out=a["alpha_out"],
                **_filter_kwargs(q, ("n_r", "n_mu", "n_traj_quad")),
                **iter_kwargs,
                **({"initial_k": initial_k} if initial_k is not None else {}),
            )
            return _wrap_mg(res, n_groups=int(sigma_t.size), **wrap_extra)
        else:
            res = gf_hollow.solve_greens_function_hollow_sphere(
                R_in=geom["R_in"],
                R_out=geom["R_out"],
                sigma_t=float(sigma_t),
                sigma_s=float(np.asarray(mats["sigma_s"])),
                nu_sigma_f=float(np.asarray(mats["nu_sigma_f"])),
                alpha_in=a["alpha_in"],
                alpha_out=a["alpha_out"],
                **_filter_kwargs(q, ("n_r", "n_mu", "n_traj_quad")),
                **iter_kwargs,
            )
            return _wrap_1g(res, n_groups=1, **wrap_extra)

    if g == "annulus":
        if is_mg:
            res = gf_annulus.solve_greens_function_annulus_mg(
                R_in=geom["R_in"],
                R_out=geom["R_out"],
                sigma_t=mats["sigma_t"],
                sigma_s=mats["sigma_s"],
                nu_sigma_f=mats["nu_sigma_f"],
                chi=mats.get("chi"),
                alpha_in=a["alpha_in"],
                alpha_out=a["alpha_out"],
                **_filter_kwargs(q, (
                    "n_r", "n_mu_axial", "n_phi_az", "n_traj_quad",
                )),
                **iter_kwargs,
                **({"initial_k": initial_k} if initial_k is not None else {}),
            )
            return _wrap_cyl_mg(res, n_groups=int(sigma_t.size), **wrap_extra)
        else:
            res = gf_annulus.solve_greens_function_annulus(
                R_in=geom["R_in"],
                R_out=geom["R_out"],
                sigma_t=float(sigma_t),
                sigma_s=float(np.asarray(mats["sigma_s"])),
                nu_sigma_f=float(np.asarray(mats["nu_sigma_f"])),
                alpha_in=a["alpha_in"],
                alpha_out=a["alpha_out"],
                **_filter_kwargs(q, (
                    "n_r", "n_mu_axial", "n_phi_az", "n_traj_quad",
                )),
                **iter_kwargs,
            )
            return _wrap_cyl_1g(res, n_groups=1, **wrap_extra)

    raise NotImplementedError(
        f"solve_critical not implemented for geometry_kind={g!r}"
    )


def _dispatch_fixed_source(
    billiard: Billiard,
    *,
    external_source: np.ndarray,
    max_iter: int | None,
    tol: float | None,
) -> FluxSolution:
    """Route fixed-source solves to the appropriate solver."""
    from orpheus.derivations.continuous.trajectory_resolvent import (
        greens_function as gf_sphere,
    )

    g = billiard.geometry_kind
    mats = billiard.xs_payload
    geom = billiard.geometry_payload
    a = billiard.alpha_payload
    q = billiard.quadrature

    iter_kwargs: dict[str, Any] = {}
    if max_iter is not None:
        iter_kwargs["max_iter"] = max_iter
    if tol is not None:
        iter_kwargs["tol"] = tol

    if g == "sphere_mr":
        res = gf_sphere.solve_greens_function_sphere_mr_fixed_source(
            radii=geom["radii"],
            sigma_t=mats["sigma_t"],
            sigma_s=mats["sigma_s"],
            external_source=external_source,
            alpha=a["alpha"],
            **_filter_kwargs(q, ("n_r", "n_mu", "n_traj_quad")),
            **iter_kwargs,
        )
        n_groups = int(
            np.asarray(external_source).reshape(-1, 1).shape[1]
            if np.asarray(external_source).ndim > 1
            else 1
        )
        # Collapse phi_g per group to a representative scalar flux for
        # the shared FluxSolution `scalar_flux` field. For multi-group
        # we expose phi_g[0] (the fast group, which carries the source
        # in the Garcia 2021 benchmarks); the full per-group profile
        # lives under metadata["phi_g"].
        phi_for_shared = (
            res.phi_g if res.phi_g.ndim == 1 else res.phi_g[0]
        )
        metadata = {
            "raw_result": res,
            "psi": res.psi_g,
            "phi": res.phi_g,
            "phi_g": res.phi_g,
            "iterations": int(res.iterations),
            "converged": bool(res.converged),
            "mesh": {
                "r_nodes": res.r_nodes,
                "mu_nodes": res.mu_nodes,
                "region_at_node": res.region_at_node,
            },
            "n_groups": n_groups,
            "geometry_kind": g,
        }
        return FluxSolution(
            spatial_nodes=res.r_nodes,
            scalar_flux=phi_for_shared,
            angular_flux=res.psi_g,
            angular_nodes=res.mu_nodes,
            eigenvalue=0.0,
            eigenvalue_kind="none",
            spatial_units="cm",
            metadata=metadata,
        )

    raise NotImplementedError(
        f"solve_fixed_source not yet implemented for "
        f"geometry_kind={g!r}. Currently only 'sphere_mr' is supported. "
        "(R2 deferred broader fixed-source unification — see "
        "trajectory_resolvent_hindsight_refactor.md R5.)"
    )


# ─────────────────────────────────────────────────────────────────────
# Result-shape adaptors — populate the SHARED CriticalSolution schema
# (eigenvalue, eigenvalue_kind, parameter_value, parameter_kind,
# converged, metadata) from each geometry-specific legacy dataclass.
# Bit-equality is preserved end-to-end: every numerical value lands in
# the unified result with the same FP bit pattern as in the legacy
# return.
# ─────────────────────────────────────────────────────────────────────


def _build_critical(
    res: Any,
    *,
    psi: np.ndarray,
    phi: np.ndarray,
    n_groups: int,
    parameter_value: float,
    parameter_kind: str,
    geometry_kind: str,
    mesh: dict[str, Any],
) -> CriticalSolution:
    """Construct the shared CriticalSolution from a legacy result.

    The legacy ``res`` is preserved in ``metadata["raw_result"]`` for
    callers who want bit-equal access to every original field; the
    canonical convenience handles ``psi``, ``phi``, ``iterations``,
    ``mesh``, ``n_groups``, ``geometry_kind`` are also populated.
    """
    metadata = {
        "raw_result": res,
        "psi": psi,
        "phi": phi,
        "iterations": int(res.iterations),
        "mesh": mesh,
        "n_groups": int(n_groups),
        "geometry_kind": geometry_kind,
    }
    return CriticalSolution(
        eigenvalue=float(res.k_eff),
        eigenvalue_kind="k_eff",
        parameter_value=float(parameter_value),
        parameter_kind=parameter_kind,
        converged=bool(res.converged),
        metadata=metadata,
    )


def _wrap_1g(
    res: Any,
    *,
    n_groups: int,
    parameter_value: float,
    parameter_kind: str,
    geometry_kind: str,
) -> CriticalSolution:
    """Wrap a 1G geometry-specific result with (r_nodes, mu_nodes)."""
    return _build_critical(
        res,
        psi=res.psi,
        phi=res.phi,
        n_groups=n_groups,
        parameter_value=parameter_value,
        parameter_kind=parameter_kind,
        geometry_kind=geometry_kind,
        mesh={
            "r_nodes": res.r_nodes,
            "mu_nodes": res.mu_nodes,
        },
    )


def _wrap_mg(
    res: Any,
    *,
    n_groups: int,
    parameter_value: float,
    parameter_kind: str,
    geometry_kind: str,
) -> CriticalSolution:
    """Wrap an MG geometry-specific result with (r_nodes, mu_nodes).

    Used by sphere MG and hollow_sphere MG (both share the
    psi_g/phi_g/r_nodes/mu_nodes shape).
    """
    return _build_critical(
        res,
        psi=res.psi_g,
        phi=res.phi_g,
        n_groups=n_groups,
        parameter_value=parameter_value,
        parameter_kind=parameter_kind,
        geometry_kind=geometry_kind,
        mesh={
            "r_nodes": res.r_nodes,
            "mu_nodes": res.mu_nodes,
        },
    )


def _wrap_mr(
    res: Any,
    *,
    n_groups: int,
    parameter_value: float,
    parameter_kind: str,
    geometry_kind: str,
) -> CriticalSolution:
    """Wrap a multi-region sphere result with region_at_node."""
    return _build_critical(
        res,
        psi=res.psi_g,
        phi=res.phi_g,
        n_groups=n_groups,
        parameter_value=parameter_value,
        parameter_kind=parameter_kind,
        geometry_kind=geometry_kind,
        mesh={
            "r_nodes": res.r_nodes,
            "mu_nodes": res.mu_nodes,
            "region_at_node": res.region_at_node,
        },
    )


def _wrap_cyl_1g(
    res: Any,
    *,
    n_groups: int,
    parameter_value: float,
    parameter_kind: str,
    geometry_kind: str,
) -> CriticalSolution:
    """Wrap cylinder / annulus 1G result (3D angular phase-space)."""
    return _build_critical(
        res,
        psi=res.psi,
        phi=res.phi,
        n_groups=n_groups,
        parameter_value=parameter_value,
        parameter_kind=parameter_kind,
        geometry_kind=geometry_kind,
        mesh={
            "r_nodes": res.r_nodes,
            "mu_axial_nodes": res.mu_axial_nodes,
            "phi_az_nodes": res.phi_az_nodes,
        },
    )


def _wrap_cyl_mg(
    res: Any,
    *,
    n_groups: int,
    parameter_value: float,
    parameter_kind: str,
    geometry_kind: str,
) -> CriticalSolution:
    """Wrap cylinder / annulus MG result."""
    return _build_critical(
        res,
        psi=res.psi_g,
        phi=res.phi_g,
        n_groups=n_groups,
        parameter_value=parameter_value,
        parameter_kind=parameter_kind,
        geometry_kind=geometry_kind,
        mesh={
            "r_nodes": res.r_nodes,
            "mu_axial_nodes": res.mu_axial_nodes,
            "phi_az_nodes": res.phi_az_nodes,
        },
    )


def _wrap_slab_1g(
    res: Any,
    *,
    n_groups: int,
    parameter_value: float,
    parameter_kind: str,
    geometry_kind: str,
) -> CriticalSolution:
    """Wrap slab / slab_asymmetric 1G result."""
    return _build_critical(
        res,
        psi=res.psi,
        phi=res.phi,
        n_groups=n_groups,
        parameter_value=parameter_value,
        parameter_kind=parameter_kind,
        geometry_kind=geometry_kind,
        mesh={
            "x_nodes": res.x_nodes,
            "mu_nodes": res.mu_nodes,
        },
    )


def _wrap_slab_mg(
    res: Any,
    *,
    n_groups: int,
    parameter_value: float,
    parameter_kind: str,
    geometry_kind: str,
) -> CriticalSolution:
    """Wrap slab / slab_asymmetric MG result."""
    return _build_critical(
        res,
        psi=res.psi_g,
        phi=res.phi_g,
        n_groups=n_groups,
        parameter_value=parameter_value,
        parameter_kind=parameter_kind,
        geometry_kind=geometry_kind,
        mesh={
            "x_nodes": res.x_nodes,
            "mu_nodes": res.mu_nodes,
        },
    )


# ─────────────────────────────────────────────────────────────────────
# Factory helpers — translate Protocol inputs (Mixture + GeometrySpec)
# into the solver-facing ``xs_payload`` / ``geometry_payload`` shapes
# that the per-geometry ``solve_greens_function_*`` entry points
# consume. The internal payload shapes are an implementation detail:
# end users address Billiard via :attr:`materials` / :attr:`geometry_spec`
# (the Protocol surface).
# ─────────────────────────────────────────────────────────────────────


def _is_mixture_dict(materials: Any) -> bool:
    """Return True iff *materials* looks like ``dict[int, Mixture]``."""
    if not isinstance(materials, dict) or not materials:
        return False
    first_key = next(iter(materials))
    if not isinstance(first_key, int):
        return False
    return isinstance(materials[first_key], Mixture)


def _infer_geometry_kind_from_spec(
    spec: GeometrySpec,
    alpha: float | dict[str, float],
) -> str:
    """Infer Billiard's internal geometry_kind from a GeometrySpec.

    Returned values are one of: ``"slab"``, ``"slab_asymmetric"``,
    ``"sphere"``, ``"cylinder"``. The asymmetric-slab branch is
    selected when *alpha* is a dict carrying ``alpha_left`` /
    ``alpha_right`` keys.
    """
    g = spec.geometry
    if g == "infinite" or g == "ISLC":
        raise ValueError(
            f"Billiard.from_problem cannot construct on geometry_spec "
            f"{g!r}; supported: slab, sphere, cylinder."
        )
    if g == "slab":
        if isinstance(alpha, dict):
            keys = set(alpha.keys())
            if "alpha_left" in keys or "alpha_right" in keys:
                return "slab_asymmetric"
        return "slab"
    if g == "sphere":
        return "sphere"
    if g == "cylinder":
        return "cylinder"
    raise ValueError(f"Billiard cannot infer geometry_kind from {g!r}")


def _geometry_payload_for_solver(
    geometry_kind: str,
    spec: GeometrySpec,
) -> dict[str, Any]:
    """Build the per-geometry payload the dispatchers consume.

    Maps the method-agnostic :class:`GeometrySpec` onto the
    geometry-specific kwargs of the underlying
    ``solve_greens_function_*`` entry points (``R`` for curvilinear
    one-endpoint billiards, ``L`` for slab and slab_asymmetric).
    """
    if spec.critical_dimension_cm is None:
        raise ValueError(
            f"Billiard requires geometry_spec.critical_dimension_cm "
            f"to be set; got None for geometry {geometry_kind!r}."
        )
    if geometry_kind in ("sphere", "cylinder"):
        return {"R": float(spec.critical_dimension_cm)}
    if geometry_kind in ("slab", "slab_asymmetric"):
        return {"L": float(spec.domain_extent_cm)}
    raise ValueError(
        f"_geometry_payload_for_solver: unsupported {geometry_kind!r}"
    )


def _mixture_to_solver_xs_payload(
    materials: dict[int, Mixture],
    geometry_kind: str,
) -> dict[str, Any]:
    """Translate a Mixture-keyed dict to the solver-facing XS payload.

    The ``solve_greens_function_*`` entry points consume raw numpy
    arrays (``sigma_t``, ``sigma_s``, ``nu_sigma_f``, ``chi``); this
    helper extracts those from the production-protocol
    :class:`Mixture` at ``materials[0]``. 1G payloads collapse to
    Python floats (preserving the underlying solver's scalar entry
    points); MG payloads carry the full ``(G,)`` and ``(G, G)``
    arrays.
    """
    if 0 not in materials:
        raise ValueError(
            f"Billiard.from_problem: materials must contain key 0 "
            f"(the active mat_id). Got keys {sorted(materials.keys())}."
        )
    mix = materials[0]
    sig_t = np.asarray(mix.SigT)
    sig_s_p0 = np.asarray(mix.SigS[0].todense())
    nu_sigma_f = np.asarray(mix.SigP)
    chi = np.asarray(mix.chi)
    if sig_t.size == 1:
        return {
            "sigma_t": float(sig_t[0]),
            "sigma_s": float(sig_s_p0[0, 0]),
            "nu_sigma_f": float(nu_sigma_f[0]),
        }
    return {
        "sigma_t": sig_t,
        "sigma_s": sig_s_p0,
        "nu_sigma_f": nu_sigma_f,
        "chi": chi,
    }
