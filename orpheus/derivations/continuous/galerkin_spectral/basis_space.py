r"""The Galerkin spectral basis space — math-heart class for galerkin_spectral.

This module is the **4th concrete instance of the math-heart pattern**
across the project, sibling to:

* :class:`~orpheus.derivations.continuous.trajectory_resolvent.Billiard`
  — the bouncing-trajectory transfer-operator resolvent (1st instance).
* :class:`~orpheus.derivations.continuous.fn_method.MomentSpace` — the
  F_N Galerkin half-range projection (2nd instance).
* The forthcoming ``Spectrum`` (Case singular eigenfunction expansion,
  3rd instance, in
  :mod:`~orpheus.derivations.continuous.singular_eigenfunction`).
* :class:`BasisSpace` — Legendre-Galerkin spectral expansion of
  Carlvik's integral equation (this module, 4th and current
  instance).

The single class :class:`BasisSpace` encapsulates everything the
Galerkin spectral method needs to commit to before computing a
solution:

* **Geometry** — slab vs sphere (cylinder is intentionally out of
  pillar; Westfall–Metcalf 1972 documents the Mitsis-style Wiener-
  Hopf reduction is non-convergent for the bare cylinder, so the
  cylinder benchmark ships through
  :mod:`~orpheus.derivations.continuous.singular_eigenfunction.cylinder`).
* **Materials** — the production-protocol cross sections specifying
  :math:`c = (\Sigma_s + \nu\Sigma_f)/\Sigma_t` and the linearly
  anisotropic scattering mean cosine
  :math:`\bar\mu = \Sigma_{s,1} / \Sigma_{s,0}`.
* **Boundary conditions** — bare-critical (vacuum) is the shipping
  configuration; reflected variants are deferred follow-up work.
* **Basis order** :math:`N` — the truncation of the Legendre-Galerkin
  moment basis. Dahl–Sjöstrand 1979 recommend :math:`N = 9` to
  saturate against their 7-significant-figure tables.

The class is a thin computational specialisation of
:class:`~orpheus.derivations.common.geometry_spec.GeometrySpec` for
the Galerkin spectral method — it is the *Galerkin spectral basis
space*, not the *Galerkin spectral solver*. The solve methods
(:meth:`BasisSpace.solve_critical`,
:meth:`BasisSpace.solve_full_spectrum`) are thin wrappers around
the existing function-level API in :mod:`...slab.one_group_anisotropic`
and :mod:`...sphere.one_group_anisotropic`. The function-level API
stays as the load-bearing implementation; ``BasisSpace`` is the
class-level facade that owns the math-rich documentation and
populates the cross-method shared
:class:`~orpheus.derivations.common.solution_types.CriticalSolution`
result type.

Architectural role and the math-heart pattern
=============================================

``BasisSpace`` is to ``galerkin_spectral`` what ``MomentSpace`` is to
``fn_method`` and what ``Billiard`` is to ``trajectory_resolvent``: a
**method-specific computational space** mounted on a method-agnostic
:class:`GeometrySpec`. All three classes answer the same question
*"What is the critical configuration?"* by populating a
:class:`CriticalSolution` carrying the eigenvalue + the parameter
that locates criticality, but each commits to a different
**mathematical attack** on the underlying transport equation:

* The F_N method (:class:`MomentSpace`) works in the **Case
  eigenfunction spectrum** — the angular flux on the boundary is
  projected onto a half-range moment basis, and criticality falls
  out of a small dense determinant condition on the collocation
  matrix.
* trajectory_resolvent (:class:`Billiard`) works in **phase space**
  — the angular flux is carried along bouncing characteristics and
  the eigenvalue is extracted from power iteration on the
  discretised Birkhoff transfer operator's resolvent.
* Galerkin spectral (:class:`BasisSpace`) works in the **full-range
  Legendre orthogonal-polynomial basis** — the spatial scalar flux
  is expanded on Legendre polynomials, the angular variable is
  pre-integrated against the kernel, and the eigenvalue spectrum
  falls out of a *standard non-symmetric matrix eigenproblem* on
  the Galerkin matrix.

The same physical problem solved three ways. The shared result
types (:class:`CriticalSolution` /
:class:`~orpheus.derivations.common.solution_types.FluxSolution`) and
the :class:`~orpheus.derivations.common.solver_protocol.TransportSolver`
Protocol are the load-bearing piece of the unification.

Distinguishing :class:`BasisSpace` from :class:`MomentSpace`
-----------------------------------------------------------

The two classes share the *moment-basis idea* but differ in:

1. **Half-range vs full-range basis.** ``MomentSpace`` uses
   half-range Legendre moments (basis on :math:`[0, 1]`,
   :math:`\mu^\alpha`); ``BasisSpace`` uses full-range Legendre
   polynomials (basis on :math:`[-1, +1]`, :math:`P_n(\mu)` for
   slab or :math:`P_n` of the reduced flux :math:`r\phi(r)` for
   sphere).
2. **Galerkin vs collocation closure.** ``MomentSpace`` collocates
   the half-range identity at :math:`N+1` discrete :math:`\xi`-
   points (Grandjean–Siewert grid for slab, Siewert–Thomas
   Chebyshev grid for sphere) — the trial space and test space
   *differ* (collocation = Petrov–Galerkin with delta tests).
   ``BasisSpace`` projects onto the same Legendre basis on both
   sides — true Galerkin orthogonality, trial = test.
3. **Dominant eigenvalue vs full spectrum.** ``MomentSpace``
   computes only the fundamental :math:`(c, a)` root via
   :math:`\det M = 0`. ``BasisSpace`` returns the **full
   :math:`2N`-eigenvalue spectrum** of the Eq. (4) block matrix,
   including complex-conjugate pairs at high anisotropy
   (Dahl–Sjöstrand Figs. 3–6) — unique among the project's
   verification pillars.
4. **Anisotropic scattering treatment.** ``MomentSpace`` ships
   1G isotropic only at the class facade (multi-group F_N is in
   flight). ``BasisSpace`` ships **linearly anisotropic
   scattering** (Sood ``*-1-1-SL/SP``) directly via the
   :math:`\bar\mu` parameter on the Eq. (4) block matrix.

References
----------

* Galerkin (1915) "Series solution of some problems of elastic
  equilibrium of rods and plates" — the variational principle
  origin (residual orthogonal to trial space).
* Carlvik, I. (1968) "A method for calculating collision
  probabilities in general cylindrical geometry and applications
  to flux distributions and Dancoff factors." *Nuclear Science
  and Engineering* **31**, 295–300. The integral-form recurrences
  for :math:`A_{m,n}` and :math:`B_{m,n}`.
* Dahl, E.B. and Sjöstrand, N.G. (1979) "Eigenvalue spectrum of
  multiplying slabs and spheres for monoenergetic neutrons with
  anisotropic scattering." *Nuclear Science and Engineering*
  **69**, 114–125. The block-matrix linearisation Eq. (4) that
  turns the :math:`c`-search into a single standard non-symmetric
  eigenproblem.
* Strang, G. and Fix, G. (1973) *An Analysis of the Finite
  Element Method*. Prentice-Hall. Chapter 1: Galerkin
  orthogonality + best-approximation theory + truncation-error
  analysis (Céa's lemma).
* Sood/Forster/Parsons (1999) LANL report LA-13511. The
  ``*-1-1-SL/SP`` benchmark cases (linearly anisotropic
  slab/sphere with :math:`P_1` scattering) that ``BasisSpace``
  targets.
* :doc:`/theory/galerkin_spectral` — the canonical theory
  exposition. The "Galerkin spectral basis space and the
  variational principle" section is the rich-narrative companion
  to this module's docstrings.
* :doc:`/theory/transport_solver_protocol` — the shared Protocol
  ``BasisSpace`` conforms to alongside ``MomentSpace``,
  ``Billiard``, and the discrete-mesh adapters.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.geometry_spec import GeometrySpec
from orpheus.derivations.common.solution_types import CriticalSolution
from orpheus.derivations.continuous.sood_registry.extractors import (
    mixture_to_fn_arrays,
)


# ════════════════════════════════════════════════════════════════════════
# Material extraction helpers
# ════════════════════════════════════════════════════════════════════════
#
# Galerkin spectral consumes (c, mu_bar, d) — a TRIPLET of scalars,
# not the full Mixture API. The conversion is ORPHEUS-Mixture →
# (sigma_t, sigma_s_p0, sigma_s_p1, nu_sigma_f) → (c, mu_bar). This
# helper bridges the gap; the function-level API stays unchanged.


def _mixture_to_c_and_mubar(mixture: Mixture) -> tuple[float, float]:
    r"""Extract :math:`(c, \bar\mu)` from a 1G ``Mixture``.

    For 1G linearly anisotropic scattering:

    .. math::

       c &= (\Sigma_s + \nu\Sigma_f) / \Sigma_t, \\
       \bar\mu &= \Sigma_{s,1} / \Sigma_s
       \quad\text{(P}_1\text{ anisotropy mean cosine)}.

    Note that Dahl–Sjöstrand's :math:`\bar\mu` is the *all-secondaries*
    mean cosine — :math:`(\Sigma_{s,1} + \nu\Sigma_{f,1}) / (c\,\Sigma_t)`
    if fission is anisotropic — whereas the project's typical
    ``mu_bar = SigS_p1[0,0] / SigS_p0[0,0]`` is the
    *scattering-only* anisotropy. For isotropic fission (``SigP_p1
    = 0``, the canonical case) the two definitions coincide. The
    convention bridge is documented in the
    `Wave 2-C Dahl-Sjostrand carlvik_galerkin` memo.
    """
    sigma_t = float(mixture.SigT[0])
    sigma_s_p0 = float(mixture.SigS[0][0, 0])
    nu_sigma_f = float(mixture.SigP[0])
    c = (sigma_s_p0 + nu_sigma_f) / sigma_t

    # Linearly anisotropic — SigS[1] is the P_1 moment if present.
    if len(mixture.SigS) >= 2:
        sigma_s_p1 = float(mixture.SigS[1][0, 0])
        mu_bar = sigma_s_p1 / sigma_s_p0 if sigma_s_p0 != 0.0 else 0.0
    else:
        mu_bar = 0.0
    return c, mu_bar


# ════════════════════════════════════════════════════════════════════════
# BasisSpace — math-heart class for the Galerkin spectral method
# ════════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class BasisSpace:
    r"""The Galerkin spectral basis space — orthogonal-polynomial expansion.

    The *Galerkin spectral method*
    (Carlvik 1968 integral form; Dahl–Sjöstrand 1979 spectral
    construction) attacks the one-speed neutron transport equation
    with linearly anisotropic scattering by **expansion of the
    spatial scalar flux on a Legendre orthogonal-polynomial basis**,
    angular pre-integration against the kernel, and Galerkin
    projection onto the same basis to produce a Galerkin matrix
    eigensystem.

    The mathematical setup
    ----------------------

    1. **Choose an orthogonal basis for the spatial variable**
       :math:`x \in [-1, +1]` (rescaled so the slab half-thickness
       :math:`a` enters via the kernel). Dahl–Sjöstrand pick the
       Legendre family :math:`\{P_n(x)\}` because:

       * They are orthogonal under the unit weight on the
         interval: :math:`\int_{-1}^{+1} P_m P_n\, dx =
         \tfrac{2}{2n+1}\delta_{mn}`.
       * Slab scalar-flux symmetry :math:`\phi(-x) = \phi(x)`
         restricts the expansion to **even-Legendre**
         polynomials :math:`P_0, P_2, \ldots, P_{2(N-1)}`. For
         sphere, the reduced flux :math:`r\,\phi(r)` is
         antisymmetric in :math:`r \to -r` so the expansion
         restricts to **odd-Legendre** polynomials :math:`P_1,
         P_3, \ldots, P_{2N-1}`.
       * The kernel :math:`E_1(\Sigma_t |x - x'|)` (slab) and
         the spherical Bickley-Naylor analog reduce to closed-
         form / quadrature-tractable inner products against
         Legendre moments — the integrand is polynomial-in-
         :math:`x \otimes` log-singular-in-:math:`|x-x'|`,
         exactly the form Carlvik's recurrences exploit.

    2. **Expand the scalar flux**:

       .. math::

           \phi(x) = \sum_{n} F_n\, (2n+1)\, P_n(x),

       (Dahl–Sjöstrand Eq. 2). The coefficients :math:`F_n` are
       the **Galerkin moments** of the scalar flux. The factor
       :math:`(2n+1)` is the orthonormalising weight; alternative
       conventions absorb it into the moment definition.

    3. **Galerkin orthogonality principle** (Galerkin 1915):
       the residual of the truncated expansion is orthogonal to
       the trial space:

       .. math::

           \langle R[\phi_N], P_m \rangle = 0
           \quad m = 0, 2, \ldots, 2(N-1)\ \text{(slab)}.

       Equivalently — the projected operator is **best-approximated
       in the Legendre basis**: among all length-:math:`N` Legendre
       coefficient vectors, the one that solves the Galerkin
       system has minimum :math:`L^2` error against the true
       transport-operator solution. This is Céa's lemma in
       finite-element language.

    4. **The trial space and the test space coincide** (true
       Galerkin) — this is what distinguishes Galerkin spectral
       from neighbouring methods:

       * **Collocation** (Petrov–Galerkin with delta-function tests):
         test = :math:`\delta(\xi - \xi_\beta)` — used by
         ``MomentSpace`` (F_N) at the
         Grandjean–Siewert / Siewert–Thomas grids.
       * **Petrov–Galerkin** (test ≠ trial): the F_N method when
         viewed in the full-range basis, since the half-range
         moments :math:`\mu^\alpha` are *not* the Legendre
         polynomials but a different polynomial family.
       * **Least-squares** (test = residual): different normal
         equations, but the matrix is :math:`A^T A` instead of
         :math:`A`.

       True Galerkin's symmetry property — :math:`A_{mn} = A_{nm}`
       when the underlying operator is self-adjoint — is only
       partially exploited here because the Carlvik integral
       form mixes :math:`E_1, E_3` kernels and the
       boundary-chord rank-1 term, which is *not* self-adjoint
       in the :math:`L^2` inner product. The matrix is
       non-symmetric; the eigenproblem is non-Hermitian. But
       Galerkin orthogonality still gives best-approximation
       in the Galerkin (energy) norm, which is what guarantees
       the convergence rate.

    5. **Carlvik's integral form** (Carlvik 1968):
       project the Boltzmann transport equation onto :math:`P_m(x)`,
       perform the angular integration analytically against the
       scattering kernel, and obtain a closed integral equation
       in Galerkin moments:

       .. math::

           2\, F_m = c\, d \sum_n \big[ A_{m,n}(a)
              - 3\bar\mu (c - 1)\, B_{m,n}(a)\big]\, F_n,

       (Dahl–Sjöstrand Eq. 3, with prefactor convention adopted
       from Eq. 1) where :math:`d = 2a` is the slab thickness
       in mean free paths and the matrix elements are

       .. math::

           A_{m,n}(a) &= \tfrac{2n+1}{2}
              \int_{-1}^{+1}\!\!\int_{-1}^{+1}
              P_m(x) P_n(y) E_1(a |x-y|)\, dy\, dx, \\
           B_{m,n}(a) &= \tfrac{2n+1}{2}\bigg[
              \int_{-1}^{+1}\!\!\int_{-1}^{+1}
              P_m(x) P_n(y) E_3(a |x-y|)\, dy\, dx \\
              &\quad - \big(\textstyle\int_{-1}^{+1} P_m(x)\,
                K_q(x; a)\, dx\big)
              \big(\textstyle\int_{-1}^{+1} P_n(y) y^q\, dy\big)
              \bigg],

       with :math:`K_q(x; a) = \tfrac{1}{2}[E_3(a|1-x|) +
       (-1)^q E_3(a|1+x|)]`, :math:`q = 0` for slab and
       :math:`q = 1` for sphere. The matrix is structurally a
       discretisation of a **Volterra-type integral operator**:
       the kernel :math:`E_1(a|x-y|)` is the analytical Peierls
       angular average, and the rank-1 :math:`B_{m,n}` correction
       handles the boundary-chord term. The exponential integrals
       :math:`E_n` decay polynomially for large :math:`a|x-y|`
       and behave like :math:`-\gamma_E - \log\tau` near the
       diagonal; the Carlvik recurrences exploit this structure
       to evaluate the matrix elements stably (see
       :mod:`...core.carlvik_recurrences`).

    6. **Block-matrix linearisation** (Dahl–Sjöstrand Eq. 4): the
       matrix Eq. (3) is :math:`(A - 3\bar\mu (c - 1) B) F =
       (1/(c d)) F`, which is a *nonlinear-in-c* eigenproblem
       (because the LHS has a :math:`(c-1)` factor multiplying
       :math:`B`). Dahl–Sjöstrand introduce the auxiliary
       variable :math:`K = c F` and observe the system

       .. math::

           \begin{pmatrix} G & H \\ I & 0 \end{pmatrix}
           \begin{pmatrix} F \\ K \end{pmatrix}
           = \frac{1}{c}
           \begin{pmatrix} F \\ K \end{pmatrix}
           \quad\text{with}\quad
           G = d (A + 3\bar\mu B),\ \ H = -3\bar\mu d B,

       is a **standard non-symmetric eigenproblem of double
       dimension** :math:`2N \times 2N`. Solving it with
       :func:`scipy.linalg.eig` returns ALL :math:`2N`
       eigenvalues :math:`(1/c_j)` in one call.

    7. **Eigenvalue spectrum and convergence** (Dahl–Sjöstrand
       1979): the Galerkin matrix has a real fundamental
       eigenvalue :math:`c_{\rm crit}` (the smallest positive
       real :math:`c_j` — the secondaries-per-collision ratio at
       criticality) and :math:`2N - 1` higher eigenvalues that
       are typically real for :math:`\bar\mu \le 0.10` and
       develop **complex-conjugate pairs** for high modes at
       :math:`\bar\mu \ge 0.15` (Dahl–Sjöstrand Figs. 3–6). The
       full spectrum carries information unavailable from F_N
       (which gives only the fundamental by construction).

       Truncation error decays as :math:`O(N^{-p})` with
       :math:`p` set by the smoothness class of the true
       scalar flux on the interval. Bare-critical
       configurations have analytic flux away from the
       boundary endpoints, so convergence is super-algebraic
       — Dahl–Sjöstrand Tables I, II reach 7-significant-figure
       agreement at :math:`N = 9, n_{\rm quad} = 128` across
       the full :math:`(\bar\mu, d)` parameter sweep.

    8. **Geometry-specific reductions**:

       * **Slab** (:math:`q = 0`, even-Legendre basis): the
         scalar flux :math:`\phi(x)` is symmetric :math:`\phi(-x)
         = \phi(x)`, so only :math:`P_0, P_2, \ldots, P_{2(N-1)}`
         appear. The kernel :math:`E_1(a|x-y|)` couples even
         moments through the symmetric integral.
       * **Sphere** (:math:`q = 1`, odd-Legendre basis): work
         in the *reduced* radial flux :math:`r\,\phi(r)`,
         which satisfies a slab-like integral equation under
         the substitution :math:`r \in [0, R] \mapsto x \in
         [-1, +1]`. The reduced flux is antisymmetric in
         :math:`r`, so only :math:`P_1, P_3, \ldots, P_{2N-1}`
         appear. The boundary-chord parameter :math:`q = 1`
         gives :math:`B_{m,n}` the antisymmetric pairing
         :math:`E_3(a|1-x|) - E_3(a|1+x|)`.

       The same Galerkin machinery serves both geometries; they
       differ by exactly two structural choices (basis parity +
       :math:`q`) and the same Eq. (4) block matrix is solved.

    9. **Cylinder is out of pillar.** The Carlvik integral form
       extends formally to cylinder via the 3-D Bickley-Naylor
       :math:`Ki_n` family, but the Wiener–Hopf reduction is
       non-convergent for the bare cylinder (Westfall–Metcalf
       1972). The cylinder benchmark
       (Sood ``Ua-1-0-CY``) ships through the
       Mitsis–Westfall–Metcalf Fredholm pillar in
       :mod:`...singular_eigenfunction.cylinder`, NOT through
       this class.

    Where this sits in the V&V triangle
    -----------------------------------

    Per :doc:`/skills/algebra-of-record`, the Galerkin spectral
    method is a **Pillar-2 (semi-analytical)** verification
    method. The Branch-1 algebra-of-record lives in
    :mod:`...origins.derivations` (eight :math:`V_{cg}` SymPy
    derivations); the Branch-2 production code lives in
    :mod:`...slab.one_group_anisotropic` and
    :mod:`...sphere.one_group_anisotropic`. ``BasisSpace`` is
    the math-rich facade above Branch-2; it does NOT introduce
    new algebra above the trusted-library line.

    The cross-check anchor is :class:`MomentSpace` (F_N) at
    :math:`\bar\mu = 0` (isotropic limit) — the two methods
    descend from completely different mathematical structures
    (full-range Legendre Galerkin vs half-range moment
    collocation) but both reduce the one-speed transport
    equation, so agreement at :math:`10^{-6}` is strong evidence
    of correctness. The L1 cross-check tests live in
    :mod:`tests.derivations.test_carlvik_galerkin_xverif_fn`.

    Connection to the Sanchez–Chandrasekhar three-meanings taxonomy
    ---------------------------------------------------------------

    The Sanchez–Chandrasekhar three-meanings taxonomy (see
    :doc:`/theory/reference_solvers` § three-meanings) locates
    each of the three math-heart classes within the same
    Green's-function landscape:

    * The F_N method (:class:`MomentSpace`) is the **spectral
      half-range realisation** of the resolvent — eigenfunction
      expansion in the Case spectrum, projected onto half-range
      Legendre moments and collocated.
    * Trajectory_resolvent (:class:`Billiard`) is the
      **path-integral (time-domain) realisation** — sum over
      bouncing characteristics weighted by the streaming
      attenuation :math:`\sum_n \alpha^n e^{-n\tau}`.
    * Galerkin spectral (:class:`BasisSpace`) is the **full-range
      polynomial-basis realisation** — the resolvent is
      represented in a finite-dimensional Legendre basis on the
      *spatial* variable, with the angular variable pre-integrated
      analytically into the kernel before the spatial Galerkin
      projection.

    All three are different angles of attack on the same Boltzmann
    equation; their cross-checks anchor the Sood ``*-1-0-SL/SP``
    and ``*-1-1-SL/SP`` truth values from three structurally-
    independent directions.

    Construction
    ------------

    Use the factory :meth:`from_problem` for the canonical
    construction path. Direct constructor calls are exposed for
    diagnostic / advanced use but are not the recommended public
    API.

    Parameters
    ----------
    geometry_spec : :class:`GeometrySpec`
        Method-agnostic geometry specification. The ``geometry``
        attribute MUST be ``"slab"`` or ``"sphere"`` —
        ``"cylinder"`` (out of pillar; Westfall–Metcalf 1972),
        ``"infinite"`` (no spatial discretisation needed —
        :math:`k_\infty` is a transfer-matrix property; use
        :func:`...common.eigenvalue.kinf_homogeneous`), and
        ``"ISLC"`` (not implemented) are rejected at construction
        time.
    materials : dict[int, Mixture]
        Production-protocol materials, keyed by material ID. The
        :attr:`GeometrySpec.mat_id` field selects the active
        mixture. The Galerkin spectral method as shipped is 1G
        only; multi-group support is deferred follow-up.
    basis_order : int, default 9
        Galerkin truncation :math:`N`. The block matrix has size
        :math:`2N \times 2N`. Defaults to 9 (Dahl–Sjöstrand's
        recommended saturating choice for 7-significant-figure
        agreement against Tables I, II).
    n_quad : int, default 128
        Outer Gauss-Legendre quadrature order for matrix-element
        evaluation. Defaults to 128 (the value Dahl–Sjöstrand
        used to reproduce their tables).
    method_name : str
        ``"galerkin_spectral"`` — the stable tag for the
        :class:`~orpheus.derivations.common.solver_protocol.TransportSolver`
        Protocol's ``method_name`` field. NOT a constructor
        argument; auto-set on every instance.

    Examples
    --------

    Bare-critical multiplying slab at :math:`\bar\mu = 0` (isotropic):

    >>> from orpheus.derivations.common.geometry_spec import GeometrySpec
    >>> from orpheus.derivations.common.xs_library import make_mixture
    >>> from orpheus.derivations.continuous.galerkin_spectral import (
    ...     BasisSpace,
    ... )
    >>> import numpy as np
    >>> mix = make_mixture(
    ...     sig_t=np.array([1.0]),
    ...     sig_c=np.array([0.2]),
    ...     sig_f=np.array([0.4]),
    ...     nu=np.array([2.0]),
    ...     chi=np.array([1.0]),
    ...     sig_s=np.array([[0.0]]),
    ... )
    >>> geom = GeometrySpec(
    ...     geometry="slab",
    ...     critical_dimension_mfp=1.0,
    ...     critical_dimension_cm=1.0,
    ...     n_groups=1,
    ... )
    >>> bs = BasisSpace.from_problem(materials={0: mix}, geometry=geom)
    >>> sol = bs.solve_critical(d=2.0)  # doctest: +SKIP
    >>> sol.eigenvalue, sol.parameter_value  # doctest: +SKIP
    (1.0, 1.2771018)

    Full eigenvalue spectrum (the Galerkin spectral distinguishing
    feature):

    >>> result = bs.solve_full_spectrum(d=2.0, mu_bar=0.30)  # doctest: +SKIP
    >>> result.eigenvalue_spectrum.shape  # doctest: +SKIP
    (18,)  # 2 * basis_order = 2*9 for n_modes=9
    """

    geometry_spec: GeometrySpec
    materials: dict[int, Mixture]
    basis_order: int = 9
    n_quad: int = 128
    method_name: str = field(default="galerkin_spectral", init=False)

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __post_init__(self) -> None:
        """Validate the Galerkin spectral method's structural preconditions.

        The method ships for **slab** and **sphere** geometries with
        **vacuum BCs** (bare-critical) and **1G linearly anisotropic
        scattering**. Out-of-scope cases are rejected with explicit
        error messages naming the alternative pillar to use.
        """
        geom = self.geometry_spec.geometry
        if geom not in {"slab", "sphere"}:
            if geom == "cylinder":
                raise ValueError(
                    "BasisSpace does not support cylinder geometry — "
                    "Westfall-Metcalf 1972 documents the Mitsis-style "
                    "Wiener-Hopf reduction is non-convergent for the bare "
                    "cylinder. Use "
                    "orpheus.derivations.continuous.singular_eigenfunction"
                    ".cylinder for the cylinder benchmark."
                )
            if geom == "infinite":
                raise ValueError(
                    "BasisSpace does not support infinite-medium cases — "
                    "k_inf is a transfer-matrix property with no spatial "
                    "Galerkin expansion. Use "
                    "orpheus.derivations.common.eigenvalue.kinf_homogeneous "
                    "or MomentSpace.solve_kinf for k_inf."
                )
            raise ValueError(
                f"BasisSpace supports geometry ∈ {{slab, sphere}}, "
                f"got {geom!r}."
            )
        if self.basis_order < 1:
            raise ValueError(
                f"basis_order must be >= 1, got {self.basis_order}"
            )
        if self.n_quad < 1:
            raise ValueError(f"n_quad must be >= 1, got {self.n_quad}")
        if self.geometry_spec.mat_id not in self.materials:
            raise ValueError(
                f"materials dict missing mat_id={self.geometry_spec.mat_id} "
                f"required by geometry_spec; got keys "
                f"{sorted(self.materials.keys())}"
            )
        # Galerkin spectral as shipped is 1G only.
        mixture = self.materials[self.geometry_spec.mat_id]
        if int(np.asarray(mixture.SigT).shape[0]) != 1:
            raise ValueError(
                "BasisSpace as shipped supports only 1G linearly "
                "anisotropic scattering. Multi-group Galerkin spectral "
                "is deferred follow-up; for multi-group k_inf use "
                "MomentSpace.solve_kinf or the multi-group F_N machinery."
            )

    @classmethod
    def from_problem(
        cls,
        materials: dict[int, Mixture],
        geometry: GeometrySpec,
        *,
        basis_order: int = 9,
        n_quad: int = 128,
    ) -> "BasisSpace":
        r"""Construct a :class:`BasisSpace` from production-protocol inputs.

        This is the recommended public construction path. It accepts
        the same ``(materials: dict[int, Mixture], GeometrySpec)``
        pair that the production CP/SN/MOC solvers consume, so a
        single problem definition can be solved by both the
        production machinery and the Galerkin spectral reference.

        Parameters
        ----------
        materials : dict[int, Mixture]
            Production-protocol materials. Single-region cases use
            one key (typically ``mat_id == 0``); multi-region is
            deferred follow-up (the Galerkin spectral integral form
            assumes homogeneous medium).
        geometry : :class:`GeometrySpec`
            Method-agnostic geometry. ``geometry.geometry`` must be
            ``"slab"`` or ``"sphere"``.
        basis_order : int, default 9
            Galerkin truncation :math:`N`.
        n_quad : int, default 128
            Outer Gauss-Legendre quadrature order.

        Returns
        -------
        :class:`BasisSpace`

        Raises
        ------
        ValueError
            If geometry is cylinder/infinite/ISLC,
            ``basis_order < 1``, ``n_quad < 1``, the ``materials``
            dict is missing ``geometry.mat_id``, or the active
            mixture is not 1G.
        """
        return cls(
            geometry_spec=geometry,
            materials=materials,
            basis_order=basis_order,
            n_quad=n_quad,
        )

    # ------------------------------------------------------------------
    # Derived primary parameters (c, mu_bar) — material-side scalars
    # ------------------------------------------------------------------

    @property
    def c(self) -> float:
        r"""Mean number of secondaries per collision, :math:`c`.

        For 1G linearly anisotropic scattering,
        :math:`c = (\Sigma_s + \nu\Sigma_f) / \Sigma_t`. This is the
        Galerkin spectral method's primary multiplication parameter
        — the entire Eq. (4) block matrix is parametrised by
        :math:`c` (through the eigenvalue) and :math:`\bar\mu`
        (through the matrix entries).
        """
        c, _ = _mixture_to_c_and_mubar(
            self.materials[self.geometry_spec.mat_id]
        )
        return c

    @property
    def mu_bar(self) -> float:
        r"""Linearly anisotropic mean scattering cosine :math:`\bar\mu`.

        Defined as :math:`\bar\mu = \Sigma_{s,1} / \Sigma_{s,0}` (the
        :math:`P_1`-anisotropy ratio, scattering-only convention).
        Returns ``0.0`` if the mixture has no :math:`P_1` moment
        (purely isotropic scattering).

        Note that Dahl–Sjöstrand's :math:`\bar\mu` is the
        all-secondaries mean cosine, which differs from this
        scattering-only convention when fission is anisotropic. For
        isotropic fission (the canonical case) the two coincide.
        """
        _, mu_bar = _mixture_to_c_and_mubar(
            self.materials[self.geometry_spec.mat_id]
        )
        return mu_bar

    @property
    def n_groups(self) -> int:
        """Number of energy groups in the active mixture (always 1 here)."""
        mixture = self.materials[self.geometry_spec.mat_id]
        return int(np.asarray(mixture.SigT).shape[0])

    # ------------------------------------------------------------------
    # solve_critical — Protocol-conforming interface
    # ------------------------------------------------------------------

    def solve_critical(
        self,
        *,
        d: Optional[float] = None,
        mu_bar: Optional[float] = None,
    ) -> CriticalSolution:
        r"""Solve the critical configuration via the Galerkin spectral method.

        Per the Dahl–Sjöstrand 1979 formulation, the eigenvalue
        :math:`c_j(\bar\mu, d)` is a function of the **geometry
        parameters** :math:`(\bar\mu, d)` only, NOT of the cross
        sections directly. The cross sections enter via the
        comparison ``c_material vs c_critical(geometry)``: at
        criticality :math:`c_{\rm material} = c_{\rm crit}(\bar\mu,
        d)`; the implied :math:`k_{\rm eff}` of the same geometry at
        non-critical materials is :math:`k_{\rm eff} =
        c_{\rm material} / c_{\rm crit}(\bar\mu, d)`.

        The :meth:`solve_critical` interface returns the
        ``c_critical`` reported by the underlying solver and packages
        it as a :class:`CriticalSolution` with:

        * ``eigenvalue = c_critical`` — the secondaries-per-collision
          ratio at criticality (the Galerkin spectral natural
          eigenvalue; NOT :math:`k_{\rm eff}`).
        * ``eigenvalue_kind = "c_critical"`` — distinguishes from
          :math:`k_{\rm eff}` / :math:`k_\infty`. Cross-method
          comparators MUST check this tag before comparing.
        * ``parameter_value = d`` — the slab thickness or sphere
          diameter (in mfp) at which criticality is reported.
        * ``parameter_kind = "domain_extent_cm"`` — semantic tag for
          the geometry's full extent. (For Galerkin spectral the
          natural unit is mean free paths, but the
          :class:`ParameterKind` literal does not include
          ``"domain_extent_mfp"``; the :attr:`metadata` dict carries
          the explicit ``"d_mfp"`` key.)

        Note on the eigenvalue convention
        ---------------------------------

        Galerkin spectral reports :math:`c_{\rm critical}`, not
        :math:`k_{\rm eff}`. The other math-heart classes
        (:class:`MomentSpace`, :class:`Billiard`) report
        :math:`k_{\rm eff}` or :math:`k_\infty`. The
        :attr:`CriticalSolution.eigenvalue_kind` field disambiguates;
        cross-method comparators that need :math:`k_{\rm eff}` from
        a Galerkin spectral solve compute :math:`k_{\rm eff} =
        c_{\rm material} / c_{\rm crit}` themselves using
        :attr:`BasisSpace.c` for :math:`c_{\rm material}`.

        Parameters
        ----------
        d : float, optional
            Slab thickness in mean free paths (full slab :math:`d =
            2 a`) or sphere diameter :math:`d = 2 R`. If ``None``,
            inferred from
            :attr:`GeometrySpec.critical_dimension_mfp` (the
            published critical configuration). Either ``d`` or the
            :class:`GeometrySpec` must supply the dimension.
        mu_bar : float, optional
            Override the linearly-anisotropic mean cosine. If
            ``None``, uses :attr:`BasisSpace.mu_bar` derived from
            the mixture's :math:`P_1` moment.

        Returns
        -------
        :class:`CriticalSolution`
            Carrying ``eigenvalue = c_critical``,
            ``eigenvalue_kind = "c_critical"``, ``parameter_value =
            d``, ``parameter_kind = "domain_extent_cm"``, and rich
            method-specific diagnostics in ``metadata``:

            * ``"raw_result"`` — the underlying
              ``CarlvikGalerkinSlabResult`` /
              ``CarlvikGalerkinSphereResult`` (preserves bit-equal
              access to the eigenvector + spectrum).
            * ``"d_mfp"`` — the mfp dimension explicitly.
            * ``"mu_bar"`` — the :math:`\bar\mu` used.
            * ``"basis_order"`` — the Galerkin order :math:`N`.
            * ``"n_quad"`` — the quadrature order.
            * ``"eigenvalue_spectrum"`` — the full :math:`2N`-complex
              spectrum (Galerkin spectral's distinguishing feature).
            * ``"c"`` — the material's :math:`c` for downstream
              comparison.

        Notes
        -----
        Bit-equality with the function-level API: this method MUST
        produce the SAME float results as a direct call to
        :func:`solve_galerkin_spectral_slab(c=..., d=..., mu_bar=...,
        n_modes=basis_order, n_quad=n_quad)` /
        :func:`solve_galerkin_spectral_sphere(...)`. The class-level
        call is a thin facade. Verified by
        :mod:`tests.derivations.test_galerkin_spectral_basis_space`
        (the foundation gate that pins the bit-equality invariant).
        """
        if d is None:
            d = self._d_from_geometry()
        if mu_bar is None:
            mu_bar = self.mu_bar

        c_material = self.c

        if self.geometry_spec.geometry == "slab":
            return self._solve_critical_slab(c_material, d, mu_bar)
        if self.geometry_spec.geometry == "sphere":
            return self._solve_critical_sphere(c_material, d, mu_bar)
        raise NotImplementedError(  # pragma: no cover (validated above)
            f"unhandled geometry {self.geometry_spec.geometry!r}"
        )

    def _d_from_geometry(self) -> float:
        r"""Infer the spatial dimension parameter ``d`` from the GeometrySpec.

        For both slab and sphere, ``d`` in the Galerkin spectral
        formulation is the **full extent in mfp**:

        * Slab: ``d = 2 * a`` where :math:`a` is the half-thickness.
          The :class:`GeometrySpec` for slab stores
          :math:`a = \mathrm{critical\_dimension\_mfp}` (per the
          F_N convention; see
          :class:`GeometrySpec.domain_extent_cm`); we double it.
        * Sphere: ``d = 2 * R`` where :math:`R` is the radius.
          Same doubling rule.

        The unit is mean free paths since the Galerkin matrix
        eigenvalue depends only on the dimensionless :math:`d` (not
        on cm).
        """
        crit_mfp = self.geometry_spec.critical_dimension_mfp
        if crit_mfp is None:
            raise ValueError(
                "GeometrySpec.critical_dimension_mfp is None; either "
                "set it on the GeometrySpec or pass d= explicitly to "
                "solve_critical."
            )
        return 2.0 * float(crit_mfp)

    def _solve_critical_slab(
        self,
        c_material: float,
        d: float,
        mu_bar: float,
    ) -> CriticalSolution:
        from .slab.one_group_anisotropic import solve_galerkin_spectral_slab

        res = solve_galerkin_spectral_slab(
            c=c_material,
            d=d,
            mu_bar=mu_bar,
            n_modes=self.basis_order,
            n_quad=self.n_quad,
            return_full_spectrum=True,
        )
        return CriticalSolution(
            eigenvalue=float(res.c_critical),
            eigenvalue_kind="c_critical",
            parameter_value=float(d),
            parameter_kind="domain_extent_cm",
            converged=True,
            metadata={
                "raw_result": res,
                "geometry": "slab",
                "d_mfp": float(d),
                "mu_bar": float(mu_bar),
                "basis_order": int(self.basis_order),
                "n_quad": int(self.n_quad),
                "eigenvalue_spectrum": res.eigenvalue_spectrum,
                "c": float(c_material),
                "method": "solve_galerkin_spectral_slab",
            },
        )

    def _solve_critical_sphere(
        self,
        c_material: float,
        d: float,
        mu_bar: float,
    ) -> CriticalSolution:
        from .sphere.one_group_anisotropic import (
            solve_galerkin_spectral_sphere,
        )

        res = solve_galerkin_spectral_sphere(
            c=c_material,
            d=d,
            mu_bar=mu_bar,
            n_modes=self.basis_order,
            n_quad=self.n_quad,
            return_full_spectrum=True,
        )
        return CriticalSolution(
            eigenvalue=float(res.c_critical),
            eigenvalue_kind="c_critical",
            parameter_value=float(d),
            parameter_kind="domain_extent_cm",
            converged=True,
            metadata={
                "raw_result": res,
                "geometry": "sphere",
                "d_mfp": float(d),
                "mu_bar": float(mu_bar),
                "basis_order": int(self.basis_order),
                "n_quad": int(self.n_quad),
                "eigenvalue_spectrum": res.eigenvalue_spectrum,
                "c": float(c_material),
                "method": "solve_galerkin_spectral_sphere",
            },
        )

    # ------------------------------------------------------------------
    # solve_full_spectrum — distinguishing feature of Galerkin spectral
    # ------------------------------------------------------------------

    def solve_full_spectrum(
        self,
        *,
        d: Optional[float] = None,
        mu_bar: Optional[float] = None,
    ) -> dict[str, Any]:
        r"""Return the full :math:`2N`-eigenvalue spectrum of the
        Eq. (4) block matrix.

        This is the **distinguishing feature** of the Galerkin
        spectral pillar relative to the other math-heart classes:
        :class:`MomentSpace` (F_N) returns only the dominant
        :math:`(c, a)` root, and :class:`Billiard`
        (trajectory_resolvent) returns only the power-iteration
        fundamental :math:`k_{\rm eff}`. ``BasisSpace`` returns
        all :math:`2N` eigenvalues — the fundamental real
        :math:`c_{\rm crit}`, the higher real eigenvalues
        corresponding to higher Legendre modes, and the
        complex-conjugate pairs that appear at high anisotropy
        (Dahl–Sjöstrand Figs. 3–6).

        Use this entry point when you need the full eigenvalue
        spectrum for Dahl–Sjöstrand-style spectrum diagnostics —
        e.g., reproducing the ``μ̄ = 0.30`` complex-pair
        appearance at high modes or computing the spectral gap
        between fundamental and first harmonic for stability
        analysis.

        Parameters
        ----------
        d : float, optional
            See :meth:`solve_critical`.
        mu_bar : float, optional
            See :meth:`solve_critical`.

        Returns
        -------
        dict
            * ``"c_critical"``: ``float`` — the fundamental real
              eigenvalue.
            * ``"eigenvalue_spectrum"``: ``np.ndarray[complex]``,
              length ``2 * basis_order``, sorted by ascending
              :math:`\Re c_j`.
            * ``"eigenvectors"``: ``np.ndarray[float]``, shape
              ``(basis_order, 2 * basis_order)``, one upper-half
              eigenvector :math:`F_j` per column.
            * ``"raw_result"``: the underlying
              ``CarlvikGalerkinSlabResult`` /
              ``CarlvikGalerkinSphereResult``.
            * ``"d_mfp"``, ``"mu_bar"``, ``"basis_order"``,
              ``"n_quad"`` — diagnostic metadata.
        """
        critical = self.solve_critical(d=d, mu_bar=mu_bar)
        raw = critical.metadata["raw_result"]
        return {
            "c_critical": float(critical.eigenvalue),
            "eigenvalue_spectrum": raw.eigenvalue_spectrum,
            "eigenvectors": raw.eigenvectors,
            "raw_result": raw,
            "d_mfp": critical.metadata["d_mfp"],
            "mu_bar": critical.metadata["mu_bar"],
            "basis_order": critical.metadata["basis_order"],
            "n_quad": critical.metadata["n_quad"],
        }


__all__ = ["BasisSpace"]
