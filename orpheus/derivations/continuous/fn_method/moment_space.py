r"""The F_N moment space — math-heart class for the F_N method.

This module is the **2nd concrete instance of the math-heart pattern**
across the project (the 1st is ``Billiard`` in
:mod:`orpheus.derivations.continuous.trajectory_resolvent`). The two
instances exist deliberately as the precondition for designing the
unifying Protocol over math-heart classes — per the project's "unify
after two instances" discipline, the Protocol is gated on two
working instances, and ``MomentSpace`` makes ``Billiard`` no longer
a one-off.

The single class :class:`MomentSpace` encapsulates everything the
F_N method needs to commit to before computing a solution:

* **Geometry** — slab vs sphere (cylinder is intentionally out of
  pillar; see :mod:`...singular_eigenfunction.cylinder`),
* **Materials** — the cross sections specifying :math:`c =
  (\Sigma_s + \nu\Sigma_f)/\Sigma_t`,
* **Boundary conditions** — bare-critical (vacuum) is the
  shipping configuration; reflected slabs use the parallel
  :class:`...slab.reflected.SlabReflectedFNResult` path,
* **F_N order** :math:`N` — the truncation of the moment basis,
* **Flux-reconstruction strategy** — none, Atkinson product-Nyström
  (the recommended ERR-036-fixed path), or the legacy plain
  Gauss–Legendre path (kept for diagnostic comparison).

The class is intentionally a thin computational specialisation of
:class:`~orpheus.derivations.common.geometry_spec.GeometrySpec` for
the F_N method — it is the *F_N moment space*, not the *F_N solver*.
The solve methods (:meth:`MomentSpace.solve_critical`,
:meth:`MomentSpace.solve_fixed_source`) are thin wrappers around the
existing function-level API in :mod:`...slab.one_group`,
:mod:`...sphere.one_group`, :mod:`...multi_group.k_inf`,
:mod:`...slab.flux_reconstruction`, and
:mod:`...sphere.flux_reconstruction`. The function-level API stays
as the load-bearing implementation; ``MomentSpace`` is the
class-level facade that owns the math-rich documentation and
populates the cross-method shared
:class:`~orpheus.derivations.common.solution_types.CriticalSolution`
/ :class:`~orpheus.derivations.common.solution_types.FluxSolution`
result types.

Architectural role and pattern
==============================

``MomentSpace`` is to F_N what ``Billiard`` is to
trajectory_resolvent: a **method-specific computational space**
mounted on a method-agnostic :class:`GeometrySpec`. Both classes
answer the same triplet of questions:

1. *"What is the critical configuration?"* —
   :meth:`solve_critical` returns a
   :class:`CriticalSolution` carrying the eigenvalue + the
   parameter that locates criticality.
2. *"Given a configuration, what is the flux shape?"* —
   :meth:`reconstruct_flux` returns a :class:`FluxSolution`
   carrying scalar and (optionally) angular flux.
3. *"At a given fixed configuration, what is the eigenvalue?"* —
   :meth:`solve_fixed_geometry_eigenvalue` reports
   :math:`k_\text{eff}` at a non-critical configuration.

The shared result types are the load-bearing piece of the
unification. Behavioural Protocols across the math-heart classes
remain deferred per the "unify after two instances" memo.

References
----------

* Siewert, Benoist (1979) "The F_N Method in Neutron-Transport
  Theory. Part I: Theory and Applications", *Nuclear Science and
  Engineering* **69**, 156–169.
* Grandjean, Siewert (1979) "The F_N Method in Neutron-Transport
  Theory. Part II: Applications and Numerical Results", *Nuclear
  Science and Engineering* **69**, 161–168.
* Siewert, Thomas (1986) "Neutron Transport Calculations in Bare
  and Reflected Spheres of Multiplying Material", *Nuclear Science
  and Engineering* **94**, 264–270.
* Kaper, Lindeman, Leaf (1974) "Benchmark Values for the Slab
  and Sphere Criticality Problem in One-Group Neutron Transport
  Theory", *Nuclear Science and Engineering* **54**, 94–98.
* Atkinson, K.E. (1997) *The Numerical Solution of Integral
  Equations of the Second Kind*. Cambridge University Press.
* :doc:`/theory/fn_method` — the canonical theory exposition; the
  "Mathematical structure of the F_N moment space" section is the
  rich-narrative companion to this module's docstrings.
* :doc:`/theory/reference_solvers` § "Three meanings of the
  Green's function" — locates F_N (Galerkin moment projection) and
  trajectory_resolvent (Billiard ray-tracing) within the same
  Green's-function landscape.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional, cast

import numpy as np

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.geometry_spec import GeometrySpec
from orpheus.derivations.common.solution_types import (
    CriticalSolution,
    FluxSolution,
)
from orpheus.derivations.continuous.sood_registry.extractors import (
    mixture_to_fn_arrays,
)


FluxReconstructionStrategy = Literal["none", "atkinson_nystrom", "legacy_gl"]


@dataclass(frozen=True)
class MomentSpace:
    r"""The F_N moment space — Galerkin half-range projection of the angular flux.

    The F_N method (Siewert–Benoist 1979) attacks the one-speed
    transport equation by **Galerkin projection of the angular flux
    :math:`\psi(r, \mu)` onto a finite-dimensional moment basis**, then
    collocation on a discrete set of :math:`\mu` values to close the
    system. ``MomentSpace`` is the data class that owns the basis
    (the F_N order :math:`N`), the geometry (slab vs sphere — the two
    differ by exactly one sign in the boundary attenuation term, see
    :doc:`/theory/fn_method` § "geometry_sign duality"), the materials
    (encoded as :math:`c = (\Sigma_s + \nu\Sigma_f)/\Sigma_t` for the
    1G isotropic case), and the boundary condition (bare-critical
    vacuum is the shipping configuration).

    The mathematical setup
    ----------------------

    1. **The angular flux** :math:`\psi(r, \mu)` lives in
       :math:`L^2([0, R] \times [-1, 1])` — a Hilbert space of
       square-integrable functions on the spatial domain :math:`\times`
       angular sphere. For slab and sphere problems with axial /
       spherical symmetry, the angular dependence collapses to a
       single polar cosine :math:`\mu \in [-1, 1]`. Cylinder
       geometry retains an azimuthal dependence and is **not in
       the F_N pillar** — see
       :mod:`orpheus.derivations.continuous.singular_eigenfunction.cylinder`
       for the Mitsis–Westfall–Metcalf alternative pillar.

    2. **Half-range decomposition.** The F_N method works in the
       Case singular-eigenfunction representation
       (Case–Zweifel 1967) where the angular flux on the boundary
       splits cleanly into outgoing (:math:`\mu > 0`) and incoming
       (:math:`\mu < 0`) half-range projections. The **Galerkin
       moments** are the weighted half-range integrals that the
       F_N basis spans:

       .. math::

           B_\alpha(\xi) = \int_0^1
               \frac{\mu^\alpha}{\mu - \xi}\, d\mu
           \quad,\quad
           A_\alpha(\xi) = \int_0^1
               \frac{\mu^\alpha\, e^{-2 a / \mu}}{\mu + \xi}\, d\mu

       (Siewert–Benoist Eq. 5–6 for slab; Siewert–Thomas Eq. 47–48
       for sphere). The **moment space** of order :math:`N` is the
       :math:`(N+1)`-dimensional vector space of expansion
       coefficients :math:`a_0, a_1, \ldots, a_N` such that

       .. math::

           \psi(-a, -\mu) = \sum_{\alpha=0}^{N} a_\alpha\, \mu^\alpha
           \quad\text{for}\quad \mu \in [0, 1].

       The boundary angular flux is approximated as a polynomial in
       the half-range cosine :math:`\mu`, with the polynomial
       degree :math:`N` setting the F_N truncation.

    3. **Galerkin orthogonality.** Under projection of
       Case–Zweifel's full-range completeness relation onto the
       moment basis, the residual of the truncated angular flux is
       constrained to be orthogonal to :math:`\{\mu^0, \ldots,
       \mu^N\}` against the half-range weight. The orthogonality
       condition is what enforces *consistency* of the
       finite-dimensional approximation — the truncation error is
       smaller in the moment basis than the trial-function class
       would predict, because the residual is projected away on
       every basis element.

    4. **Collocation closure.** With :math:`N+1` unknowns
       :math:`(a_0, \ldots, a_N)`, we need :math:`N+1` equations.
       The F_N method evaluates the boundary integral identity at
       :math:`N+1` collocation points :math:`\xi_\beta`. Different
       choices give different convergence rates but converge to
       the same limit:

       * **Slab (Grandjean–Siewert prescription)**: :math:`\xi_0 =
         \nu_0` (the Case discrete eigenvalue), :math:`\xi_1 = 0`
         (the half-range origin), :math:`\xi_2 = 1` (the
         half-range maximum), and :math:`N - 2` interior points
         equally spaced in :math:`(0, 1)`.
       * **Sphere (Siewert–Thomas Eq. 38a prescription)**:
         :math:`\xi_0 = \nu_0`, with the remaining :math:`N`
         points on the **Chebyshev-interior** grid
         :math:`\xi_k = \cos((2k - 1)\pi / (2N + 2))` for
         :math:`k = 1, \ldots, N`. The sphere needs the
         Chebyshev grid because the slab Grandjean–Siewert grid
         (with :math:`\xi = 0` and :math:`\xi = 1` endpoints)
         creates a rank-deficient system under the sphere's
         ``geometry_sign = -1`` matrix entries.

       The collocation matrix is :math:`(N+1) \times (N+1)`
       complex (the :math:`\xi_0 = \nu_0 = i u_0` row is genuinely
       imaginary for multiplying media :math:`c > 1`), generally
       dense, but small at the operating point :math:`N = 8`–:math:`12`.

    5. **The critical condition.** Criticality is the existence of
       a non-trivial null vector :math:`(a_0, \ldots, a_N)` of the
       collocation matrix :math:`M`, which requires
       :math:`\det M = 0`. The F_N method root-finds
       :math:`\det M(a) = 0` (slab) or :math:`\det M(R) = 0`
       (sphere) on the configuration parameter to locate the
       critical configuration.

    6. **Truncation error analysis.** The Galerkin projection
       error in :math:`L^2` decays as :math:`O(N^{-p})` where
       :math:`p` depends on the smoothness class of the true
       angular flux on the half-range
       :math:`\psi(-a, -\mu) \in C^p([0, 1])`. For the bare-critical
       slab and sphere, the angular flux at the boundary is
       analytic in :math:`\mu` away from :math:`\mu = 0` (the
       grazing-ray limit), so the convergence is super-algebraic
       — Grandjean–Siewert Table XI shows :math:`N = 5` reaches
       4–5 sig figs and :math:`N = 8` reaches :math:`10^{-5}` to
       :math:`10^{-6}` absolute accuracy on the critical
       half-thickness across the :math:`c`-sweep.

    7. **Why N = 8–12 is the typical operating point.** Empirical
       experiment + the smoothness assumption: at :math:`N = 8`
       the system is small enough to assemble in microseconds and
       large enough to give 6-digit agreement with the Sood
       LA-13511 truth values. Beyond :math:`N \approx 16` the
       collocation matrix becomes ill-conditioned (the Chebyshev
       columns become near-linearly-dependent in finite
       precision), so refinement saturates at the precision floor
       rather than continuing to gain digits.

    8. **Multi-region / multi-group extensions.** Linearity of the
       moment projection extends:

       * **Multi-group** to a tensor product of (group-coupling
         matrix) :math:`\otimes` (moment basis). Sood Eq. 76 and
         the Siewert–Thomas 1986 2G machinery realise this; the
         current ``MomentSpace`` exposes the 1G specialisation
         and delegates the multi-group :math:`k_\infty` evaluation
         to :func:`...multi_group.k_inf.compute_kinf_*`.
       * **Multi-region** to a block transfer matrix between
         region interfaces. Neshat–Maiorino 1980 implements this
         for the 1G reflected slab; access via the
         ``"reflected_slab"`` factory variant.

    9. **Connection to flux reconstruction.** The F_N coefficients
       :math:`a_\alpha` define :math:`\psi(\pm a, \mu)` on the
       boundary; reconstructing :math:`\phi(z)` or :math:`\psi(z,
       \mu)` in the interior requires either:

       * **KLL Fredholm iteration** (Path B, Kaper–Lindeman–Leaf
         1974) — structurally independent of F_N; iterates the
         scalar-flux integral equation with the F_N boundary
         angular flux as the source. Implemented in
         :func:`...slab.flux_reconstruction.slab_scalar_flux_kll`
         and the sphere analogue.
       * **F_N projection iteration** (Path A.i) — uses the F_N
         coefficients directly as the angular flux throughout the
         domain. Requires the Atkinson product-Nyström treatment
         of the log-singular Peierls kernel (ERR-036) to reach
         :math:`10^{-5}` accuracy; the legacy plain Gauss–Legendre
         path saturates at :math:`1`–:math:`7\,\%` due to silent
         diagonal truncation of :math:`E_1(0) = +\infty`. See
         :ref:`fn-method-atkinson-product-nystrom` and the
         ``flux_reconstruction`` field of this class.

    10. **The relationship to** ``Billiard`` **(trajectory_resolvent).**
        Both methods solve the same boundary-value problem on the
        same physical configuration, but attack different
        mathematical structures:

        * The F_N method works in the **Case eigenfunction
          spectrum** — the angular flux is projected onto a
          finite-dimensional moment basis and the eigenvalue
          falls out of a small dense determinant condition.
        * Trajectory_resolvent (``Billiard``) works in **phase
          space** — the angular flux is carried along bouncing
          characteristics and the eigenvalue is extracted from
          power iteration on the discretised resolvent.

        The Sanchez–Chandrasekhar three-meanings taxonomy (see
        :doc:`/theory/reference_solvers` § three-meanings)
        locates both methods within the same Green's-function
        landscape: the F_N method is the *spectral* realisation
        of the resolvent (eigenfunction expansion in the Case
        spectrum), while trajectory_resolvent is the
        *path-integral* realisation (sum over bouncing
        characteristics weighted by attenuation). They cross-check
        each other above the trusted-library line — both consume
        ``numpy``/``scipy``, neither shares any in-house primitive
        — and their cross-method gates anchor the verification
        chain for the Sood LA-13511 truth set.

    Cross-method analog
    -------------------

    ``MomentSpace`` is to fn_method what ``Billiard`` is to
    trajectory_resolvent (a billiard system with multi-bounce
    resolvent), what ``Spectrum`` will eventually be to
    singular_eigenfunction (Case singular eigenfunction expansion),
    what ``LegendreBlock`` will be to carlvik_galerkin
    (Dahl–Sjöstrand block-matrix linearisation), and what
    ``CPMesh`` is to the production collision-probability solver.
    These siblings are method-specific computational
    specialisations of the abstract ``GeometrySpec``; see
    :doc:`/theory/reference_solvers` § "method-specific spaces"
    for the unification design (gated on ≥2 instances, of which
    ``MomentSpace`` and ``Billiard`` are the first two).

    Construction
    ------------

    Use the factory :meth:`from_problem` for the canonical
    construction path. A direct constructor call is exposed for
    diagnostic / advanced use but is not the recommended public
    API.

    Parameters
    ----------
    geometry : :class:`GeometrySpec`
        Method-agnostic geometry specification. The ``geometry``
        attribute MUST be ``"slab"`` or ``"sphere"`` —
        ``"infinite"`` (k_inf only), ``"cylinder"`` (out of pillar),
        and ``"ISLC"`` (not implemented) are rejected at
        construction time.
    materials : dict[int, Mixture]
        Production-protocol materials, keyed by material ID. The
        :attr:`GeometrySpec.mat_id` field selects the active
        mixture for the F_N solve.
    fn_order : int
        F_N truncation order :math:`N`. The collocation system has
        size :math:`(N+1) \times (N+1)`. Defaults to 9 (the typical
        operating point).
    flux_reconstruction : :data:`FluxReconstructionStrategy`
        Which flux-reconstruction path to use when
        :meth:`reconstruct_flux` is called:

        * ``"atkinson_nystrom"`` (default) — Atkinson 1972/1997
          §4.2 product-Simpson rule on the log-singular kernel,
          with closed-form treatment of the diagonal. Recommended;
          fixes ERR-036.
        * ``"legacy_gl"`` — plain Gauss–Legendre quadrature on
          the (z, μ) tensor product. Diagnostic only; saturates
          at 1–7\,\% accuracy due to silent diagonal truncation.
        * ``"none"`` — :meth:`reconstruct_flux` raises; only
          :meth:`solve_critical` is callable.
    """

    geometry: GeometrySpec
    materials: dict[int, Mixture]
    fn_order: int = 9
    flux_reconstruction: FluxReconstructionStrategy = "atkinson_nystrom"

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __post_init__(self) -> None:
        """Validate the F_N method's structural preconditions.

        The F_N method as shipped applies to:

        * **Slab** (Siewert–Benoist 1979 / Grandjean–Siewert 1979).
        * **Sphere** (Siewert–Thomas 1986).
        * **Infinite medium** (Sood Eq. 19/28/29/76 — k_inf only).

        It is **explicitly out of pillar** for cylinder geometry —
        Westfall–Metcalf 1972 documents that the Mitsis-style
        Wiener-Hopf reduction is non-convergent for the bare
        cylinder. Cylinder critical dimensions ship via
        :mod:`...singular_eigenfunction.cylinder` on the
        Mitsis–Westfall–Metcalf Fredholm pillar.
        """
        if self.geometry.geometry not in {"slab", "sphere", "infinite"}:
            raise ValueError(
                f"MomentSpace supports geometry ∈ {{slab, sphere, infinite}}, "
                f"got {self.geometry.geometry!r}. Cylinder is out of pillar "
                f"(Westfall-Metcalf 1972 — see singular_eigenfunction.cylinder)."
            )
        if self.fn_order < 0:
            raise ValueError(f"fn_order must be ≥ 0, got {self.fn_order}")
        if self.geometry.mat_id not in self.materials:
            raise ValueError(
                f"materials dict missing mat_id={self.geometry.mat_id} "
                f"required by geometry_spec; got keys "
                f"{sorted(self.materials.keys())}"
            )
        if self.flux_reconstruction not in {"none", "atkinson_nystrom", "legacy_gl"}:
            raise ValueError(
                f"flux_reconstruction must be one of "
                f"{{none, atkinson_nystrom, legacy_gl}}, got "
                f"{self.flux_reconstruction!r}"
            )

    @classmethod
    def from_problem(
        cls,
        materials: dict[int, Mixture],
        geometry: GeometrySpec,
        *,
        fn_order: int = 9,
        flux_reconstruction: FluxReconstructionStrategy = "atkinson_nystrom",
    ) -> "MomentSpace":
        r"""Construct a :class:`MomentSpace` from production-protocol inputs.

        This is the recommended public construction path. It accepts
        the same ``materials: dict[int, Mixture]`` + ``GeometrySpec``
        pair that the production CP/SN/MOC solvers consume, so a
        single problem definition can be solved by both production
        machinery and the F_N reference without re-deriving any
        cross sections.

        Parameters
        ----------
        materials : dict[int, Mixture]
            Production-protocol materials. Multi-region cases use
            multiple keys (one per ``mat_id``); bare-critical
            cases typically have a single key (``mat_id == 0``).
        geometry : :class:`GeometrySpec`
            Method-agnostic geometry. ``geometry.geometry`` must be
            ``"slab"``, ``"sphere"``, or ``"infinite"``. Cylinder
            is rejected (out of pillar).
        fn_order : int, default 9
            F_N order :math:`N`.
        flux_reconstruction : :data:`FluxReconstructionStrategy`, default ``"atkinson_nystrom"``
            Which flux-reconstruction path to use.

        Returns
        -------
        :class:`MomentSpace`

        Raises
        ------
        ValueError
            If geometry is cylinder/ISLC, ``fn_order < 0``, the
            ``materials`` dict is missing ``geometry.mat_id``, or
            ``flux_reconstruction`` is not one of the supported
            strategies.
        """
        return cls(
            geometry=geometry,
            materials=materials,
            fn_order=fn_order,
            flux_reconstruction=flux_reconstruction,
        )

    # ------------------------------------------------------------------
    # Derived primary parameter
    # ------------------------------------------------------------------

    @property
    def c(self) -> float:
        r"""Mean number of secondaries per collision, :math:`c`.

        For 1G isotropic-scattering problems
        :math:`c = (\Sigma_s + \nu\Sigma_f)/\Sigma_t`. This is the
        F_N method's primary input parameter — the entire collocation
        machinery (dispersion relation, X-function, moment integrals)
        is parametrised by :math:`c` alone.

        For multi-group problems (handled internally by
        :func:`...multi_group.k_inf.compute_kinf_*`), the scalar
        :math:`c` is not the right primitive — the multi-group
        machinery consumes :math:`(\Sigma_t, \Sigma_s, \nu\Sigma_f,
        \chi)` arrays directly. Multi-group ``MomentSpace`` instances
        for which this property is accessed will raise.

        Raises
        ------
        ValueError
            If the active mixture has more than one energy group.
        """
        mixture = self.materials[self.geometry.mat_id]
        sigma_t, sigma_s_p0, nu_sigma_f, _chi = mixture_to_fn_arrays(mixture)
        if sigma_t.shape[0] != 1:
            raise ValueError(
                f"MomentSpace.c requires a 1G mixture (got "
                f"{sigma_t.shape[0]}G). Multi-group problems should call "
                f"solve_kinf() directly."
            )
        return float((sigma_s_p0[0, 0] + nu_sigma_f[0]) / sigma_t[0])

    @property
    def n_groups(self) -> int:
        """Number of energy groups in the active mixture."""
        mixture = self.materials[self.geometry.mat_id]
        return int(np.asarray(mixture.SigT).shape[0])

    # ------------------------------------------------------------------
    # solve_critical — the canonical critical-configuration call
    # ------------------------------------------------------------------

    def solve_critical(
        self,
        *,
        n_bracket: Optional[int] = None,
        bisect_tol: float = 1e-12,
        max_bisect: int = 80,
    ) -> CriticalSolution:
        r"""Solve the critical configuration for the F_N moment space.

        For a **bare-critical slab** (``geometry == "slab"``,
        ``bc_left == bc_right == BC.vacuum``) this returns the
        critical half-thickness :math:`a` in mean free paths,
        wrapped as a :class:`CriticalSolution` with
        ``parameter_kind == "half_thickness_mfp"`` and
        ``eigenvalue == 1.0`` (the eigenvalue at criticality).

        For a **bare-critical sphere** (``geometry == "sphere"``,
        ``bc_right == BC.vacuum``) this returns the critical
        radius :math:`R_c` in mean free paths,
        ``parameter_kind == "radius_mfp"``.

        For an **infinite-medium** problem (``geometry ==
        "infinite"``) the "critical configuration" is just
        :math:`k_\infty`; the result has ``parameter_kind ==
        "k_inf_only"`` and ``parameter_value == 0.0`` (no spatial
        configuration parameter).

        The underlying solve delegates to:

        * :func:`...slab.one_group.solve_fn_slab_bare_critical` for
          slab,
        * :func:`...sphere.one_group.solve_fn_sphere_bare_critical`
          for sphere,
        * :func:`...multi_group.k_inf.compute_kinf_1g` (and the 2G/mG
          variants) for infinite medium.

        The rich method-specific result (``SlabFNResult`` /
        ``SphereFNResult`` / scalar :math:`k_\infty`) is preserved
        in :attr:`CriticalSolution.metadata` under the key
        ``"raw_result"`` for callers that need access to the F_N
        expansion coefficients, the dispersion-relation root
        :math:`u_0`, or the determinant residual.

        Parameters
        ----------
        n_bracket : int | None, default None
            Initial bracket-scan resolution. ``None`` means use the
            geometry-specific default (slab: 400; sphere: 800), which
            preserves bit-equality with the function-level API's
            default behaviour. Pass an explicit value to override.
        bisect_tol : float, default 1e-12
            Bisection tolerance on the configuration parameter.
        max_bisect : int, default 80
            Maximum bisection iterations.

        Returns
        -------
        :class:`CriticalSolution`

        Notes
        -----
        Bit-equality with the function-level API: this method
        ``solve_critical()`` MUST produce the SAME float results
        as a direct call to
        :func:`solve_fn_slab_bare_critical(c=..., n_modes=fn_order, ...)`
        / :func:`solve_fn_sphere_bare_critical(...)` etc. The
        class-level call is a thin facade. Verified by
        :mod:`tests.derivations.test_moment_space` (the foundation
        gate that pins the bit-equality invariant).
        """
        geom = self.geometry.geometry
        mixture = self.materials[self.geometry.mat_id]
        sigma_t, sigma_s_p0, nu_sigma_f, chi = mixture_to_fn_arrays(mixture)
        n_groups = sigma_t.shape[0]

        if geom == "infinite":
            return self._solve_critical_infinite(
                sigma_t, sigma_s_p0, nu_sigma_f, chi, n_groups
            )

        if n_groups != 1:
            raise NotImplementedError(
                "MomentSpace.solve_critical for slab/sphere is currently "
                "1G-only — the multi-group F_N spatial extension "
                "(Siewert-Thomas 1986 2G machinery) is not yet wired "
                "through this class facade. For 1G slab/sphere, this "
                "produces SlabFNResult / SphereFNResult bit-for-bit; "
                "for multi-group infinite medium use k_inf-only solves."
            )

        c = float((sigma_s_p0[0, 0] + nu_sigma_f[0]) / sigma_t[0])
        if c <= 1.0:
            raise ValueError(
                f"F_N bare-critical {geom} requires c > 1 "
                f"(multiplying medium); got c={c} from mixture "
                f"sigma_s + nu_sigma_f = {sigma_s_p0[0, 0] + nu_sigma_f[0]}, "
                f"sigma_t = {sigma_t[0]}."
            )

        if geom == "slab":
            return self._solve_critical_slab(c, n_bracket, bisect_tol, max_bisect)
        if geom == "sphere":
            return self._solve_critical_sphere(c, n_bracket, bisect_tol, max_bisect)

        raise NotImplementedError(
            f"MomentSpace.solve_critical: unhandled geometry {geom!r}"
        )

    def _solve_critical_infinite(
        self,
        sigma_t: np.ndarray,
        sigma_s_p0: np.ndarray,
        nu_sigma_f: np.ndarray,
        chi: np.ndarray,
        n_groups: int,
    ) -> CriticalSolution:
        r"""Infinite-medium :math:`k_\infty` via the multi-group formulae.

        Delegates to :func:`...multi_group.k_inf.compute_kinf_*`
        based on group count. The group convention is the
        ORPHEUS-ordered ``[from, to]`` scattering matrix; Sood
        relabelling happens inside the multi-group functions.
        """
        from .multi_group.k_inf import (
            compute_kinf_1g,
            compute_kinf_2g_general,
            compute_kinf_mg,
        )

        if n_groups == 1:
            k_inf = compute_kinf_1g(
                float(sigma_t[0]),
                float(sigma_s_p0[0, 0]),
                float(nu_sigma_f[0]),
            )
            metadata = {
                "n_groups": 1,
                "method": "compute_kinf_1g",
                "raw_result": k_inf,
            }
        elif n_groups == 2:
            # Use the general (Eq 28) formula — it reduces to Eq 29 when
            # there's no upscatter, so it's the single source of truth.
            k_inf = compute_kinf_2g_general(sigma_t, sigma_s_p0, nu_sigma_f, chi)
            metadata = {
                "n_groups": 2,
                "method": "compute_kinf_2g_general",
                "raw_result": k_inf,
            }
        else:
            k_inf = compute_kinf_mg(sigma_t, sigma_s_p0, nu_sigma_f, chi)
            metadata = {
                "n_groups": n_groups,
                "method": "compute_kinf_mg",
                "raw_result": k_inf,
            }

        return CriticalSolution(
            eigenvalue=k_inf,
            eigenvalue_kind="k_inf",
            parameter_value=0.0,
            parameter_kind="k_inf_only",
            converged=True,
            metadata=metadata,
        )

    def _solve_critical_slab(
        self,
        c: float,
        n_bracket: Optional[int],
        bisect_tol: float,
        max_bisect: int,
    ) -> CriticalSolution:
        from .slab.one_group import solve_fn_slab_bare_critical

        # Pass-through kwargs only when explicitly set; otherwise the
        # underlying solver's defaults take effect (preserves bit-equality
        # for callers that don't pass overrides).
        kwargs: dict = {"c": c, "n_modes": self.fn_order,
                        "bisect_tol": bisect_tol, "max_bisect": max_bisect}
        if n_bracket is not None:
            kwargs["n_bracket"] = n_bracket
        res = solve_fn_slab_bare_critical(**kwargs)
        return CriticalSolution(
            eigenvalue=1.0,
            eigenvalue_kind="k_eff",
            parameter_value=float(res.a_critical_mfp),
            parameter_kind="half_thickness_mfp",
            converged=True,  # bisection always converges given bracket
            metadata={
                "n_groups": 1,
                "method": "solve_fn_slab_bare_critical",
                "n_modes": int(res.N),
                "c": float(res.c),
                "nu0": float(res.nu0),
                "determinant_residual": complex(res.determinant_residual),
                "n_bracket_points": int(res.n_bracket_points),
                "raw_result": res,
            },
        )

    def _solve_critical_sphere(
        self,
        c: float,
        n_bracket: Optional[int],
        bisect_tol: float,
        max_bisect: int,
    ) -> CriticalSolution:
        from .sphere.one_group import solve_fn_sphere_bare_critical

        kwargs: dict = {"c": c, "n_modes": self.fn_order,
                        "bisect_tol": bisect_tol, "max_bisect": max_bisect}
        if n_bracket is not None:
            kwargs["n_bracket"] = n_bracket
        res = solve_fn_sphere_bare_critical(**kwargs)
        return CriticalSolution(
            eigenvalue=1.0,
            eigenvalue_kind="k_eff",
            parameter_value=float(res.R_critical_mfp),
            parameter_kind="radius_mfp",
            converged=True,
            metadata={
                "n_groups": 1,
                "method": "solve_fn_sphere_bare_critical",
                "n_modes": int(res.N),
                "c": float(res.c),
                "nu0": float(res.nu0),
                "determinant_residual": complex(res.determinant_residual),
                "n_bracket_points": int(res.n_bracket_points),
                "raw_result": res,
            },
        )

    # ------------------------------------------------------------------
    # solve_kinf — convenience for k_inf-only access
    # ------------------------------------------------------------------

    def solve_kinf(self) -> float:
        r"""Return :math:`k_\infty` for the active mixture.

        Convenience wrapper: equivalent to
        ``self.solve_critical().eigenvalue`` when ``geometry ==
        "infinite"``, or the underlying mixture's :math:`k_\infty`
        regardless of geometry (since :math:`k_\infty` is a material
        property independent of configuration).
        """
        mixture = self.materials[self.geometry.mat_id]
        sigma_t, sigma_s_p0, nu_sigma_f, chi = mixture_to_fn_arrays(mixture)
        n_groups = sigma_t.shape[0]
        from .multi_group.k_inf import (
            compute_kinf_1g,
            compute_kinf_2g_general,
            compute_kinf_mg,
        )
        if n_groups == 1:
            return compute_kinf_1g(
                float(sigma_t[0]),
                float(sigma_s_p0[0, 0]),
                float(nu_sigma_f[0]),
            )
        if n_groups == 2:
            return compute_kinf_2g_general(
                sigma_t, sigma_s_p0, nu_sigma_f, chi
            )
        return compute_kinf_mg(sigma_t, sigma_s_p0, nu_sigma_f, chi)

    # ------------------------------------------------------------------
    # Flux reconstruction
    # ------------------------------------------------------------------

    def reconstruct_flux(
        self,
        *,
        n_panels: int = 256,
        z_eval: Optional[np.ndarray] = None,
    ) -> FluxSolution:
        r"""Reconstruct the scalar flux from F_N moments.

        For a bare-critical slab or sphere, the F_N moments
        :math:`(a_0, \ldots, a_N)` returned by
        :meth:`solve_critical` define the boundary angular flux
        :math:`\psi(\pm a, \mu)` as a polynomial in :math:`\mu`.
        Reconstructing the interior scalar flux :math:`\phi(z)`
        requires evaluating the Peierls integral

        .. math::

            \phi(z) = \frac{c}{2} \int_{-a}^{a}
                E_1(|z - z'|)\, \phi(z')\, dz' + \text{boundary terms}

        which is a Fredholm equation of the second kind with a
        log-singular kernel :math:`E_1(\tau)`. The choice of
        quadrature for the singular kernel sets the achievable
        accuracy:

        * ``"atkinson_nystrom"`` — Atkinson 1972/1997 §4.2 product-
          Simpson rule. Integrates :math:`\int \log|t - s|\, P_2(s)\,
          ds` analytically against the piecewise-quadratic Lagrange
          basis on each Simpson panel; the smooth remainder
          :math:`R(\tau) = E_1(\tau) + \log\tau` is integrated with
          standard Simpson. Achieves :math:`O(h^4 \log h)`
          convergence on the operator (de Hoog & Weiss 1973).
          **Recommended; fixes ERR-036.**

        * ``"legacy_gl"`` — plain Gauss–Legendre on the (z, μ)
          tensor product. Saturates at 1–7\,\% accuracy due to
          silent diagonal truncation of :math:`E_1(0) = +\infty`
          (the :math:`\mu`-quadrature collapses to a finite
          :math:`\sim 2\log(n_\mu)` value at :math:`\tau = 0`).
          See Numerical Bug Signatures § Signature 6 (the textbook
          fingerprint :math:`\mathrm{err} \cdot n / \log n \approx
          \mathrm{const}`). Diagnostic only — kept for comparison.

        * ``"none"`` — :meth:`reconstruct_flux` raises. Use this
          when you only need :meth:`solve_critical` to skip the
          flux-reconstruction machinery.

        Parameters
        ----------
        n_panels : int, default 256
            Number of Simpson panels for ``"atkinson_nystrom"``;
            number of Gauss–Legendre nodes for ``"legacy_gl"``. The
            two paths use different quadrature meanings; the
            parameter name is unified to keep the API simple.
        z_eval : np.ndarray | None
            Spatial nodes (in mfp) at which to evaluate the
            reconstructed flux. If ``None``, defaults to a uniform
            grid on :math:`[-a, a]` (slab) or :math:`[0, R]` (sphere)
            with ``n_panels`` nodes.

        Returns
        -------
        :class:`FluxSolution`

        Raises
        ------
        RuntimeError
            If :meth:`solve_critical` has not been called yet (the
            F_N coefficients are needed for reconstruction).
        ValueError
            If ``flux_reconstruction == "none"``.

        Notes
        -----
        Currently 1G slab and 1G sphere are implemented. Multi-group
        flux reconstruction lives behind future work; the symbolic
        framework for it is documented in
        :func:`...origins.fn_flux_reconstruction_derivations.derive_*`
        but the Branch-2 wiring is part of the multi-group F_N
        spatial-extension follow-up.
        """
        if self.flux_reconstruction == "none":
            raise ValueError(
                "MomentSpace was constructed with flux_reconstruction='none'; "
                "reconstruct_flux is not callable. Reconstruct with a different "
                "MomentSpace instance or change the flux_reconstruction strategy."
            )

        # Solve criticality first to get the F_N expansion coefficients.
        critical = self.solve_critical()
        raw = critical.metadata.get("raw_result")
        if raw is None:
            raise RuntimeError(
                "MomentSpace.reconstruct_flux: solve_critical did not "
                "return a raw_result; cannot reconstruct without F_N "
                "expansion coefficients."
            )

        geom = self.geometry.geometry
        if geom == "slab":
            return self._reconstruct_slab(raw, n_panels, z_eval, critical)
        if geom == "sphere":
            return self._reconstruct_sphere(raw, n_panels, z_eval, critical)
        raise NotImplementedError(
            f"MomentSpace.reconstruct_flux: unhandled geometry {geom!r}. "
            f"Infinite-medium has no flux to reconstruct."
        )

    def _reconstruct_slab(
        self,
        raw,
        n_panels: int,
        z_eval: Optional[np.ndarray],
        critical: CriticalSolution,
    ) -> FluxSolution:
        from .slab.flux_reconstruction import (
            slab_scalar_flux_fn_projection,
            slab_scalar_flux_fn_projection_atkinson,
        )

        a = float(raw.a_critical_mfp)
        if z_eval is None:
            z_eval = np.linspace(-a, a, n_panels)
        z_eval_arr = np.asarray(z_eval, dtype=float)

        if self.flux_reconstruction == "atkinson_nystrom":
            phi = slab_scalar_flux_fn_projection_atkinson(
                raw, z_eval_arr, n_panels=n_panels
            )
        elif self.flux_reconstruction == "legacy_gl":
            phi = slab_scalar_flux_fn_projection(
                raw, z_eval_arr, n_quad_z=n_panels
            )
        else:
            raise ValueError(
                f"unhandled flux_reconstruction strategy "
                f"{self.flux_reconstruction!r}"
            )

        return FluxSolution(
            spatial_nodes=z_eval_arr,
            scalar_flux=np.asarray(phi, dtype=float),
            angular_flux=None,
            angular_nodes=None,
            eigenvalue=critical.eigenvalue,
            eigenvalue_kind=critical.eigenvalue_kind,
            spatial_units="mfp",
            metadata={
                "geometry": "slab",
                "flux_reconstruction": self.flux_reconstruction,
                "n_panels": n_panels,
                "n_modes": self.fn_order,
                "a_critical_mfp": a,
                "c": critical.metadata.get("c"),
            },
        )

    def _reconstruct_sphere(
        self,
        raw,
        n_panels: int,
        z_eval: Optional[np.ndarray],
        critical: CriticalSolution,
    ) -> FluxSolution:
        from .sphere.flux_reconstruction import (
            solve_kll_sphere_continuum_coefficient,
            sphere_scalar_flux_kll,
        )

        R = float(raw.R_critical_mfp)
        if z_eval is None:
            z_eval = np.linspace(0.0, R, n_panels)
        z_eval_arr = np.asarray(z_eval, dtype=float)

        # Sphere flux reconstruction uses the KLL Fredholm path
        # (Path B in the fn_method theory page) — Path A.i for sphere
        # would require its own Atkinson product-Nyström treatment of
        # the spherical Peierls kernel, which is not yet shipped.
        # KLL is structurally independent of F_N; it iterates the
        # Fredholm equation directly in (R, c) without consuming the
        # F_N expansion coefficients (the F_N solve here only locates
        # criticality — the KLL path then reconstructs the eigenmode
        # at that critical radius from a different mathematical route).
        c = float(critical.metadata.get("c", self.c))
        kll = solve_kll_sphere_continuum_coefficient(
            R, c, n_nodes=n_panels
        )
        phi = sphere_scalar_flux_kll(kll, z_eval_arr)

        return FluxSolution(
            spatial_nodes=z_eval_arr,
            scalar_flux=np.asarray(phi, dtype=float),
            angular_flux=None,
            angular_nodes=None,
            eigenvalue=critical.eigenvalue,
            eigenvalue_kind=critical.eigenvalue_kind,
            spatial_units="mfp",
            metadata={
                "geometry": "sphere",
                "flux_reconstruction": "kll_fredholm",
                "n_panels": n_panels,
                "n_modes": self.fn_order,
                "R_critical_mfp": R,
                "c": critical.metadata.get("c"),
            },
        )


__all__ = [
    "FluxReconstructionStrategy",
    "MomentSpace",
]
