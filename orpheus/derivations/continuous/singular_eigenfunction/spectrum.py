r"""The transport operator's singular eigenfunction spectrum (Case 1960).

This module is the **3rd concrete instance of the math-heart pattern**
across the project. The 1st is ``Billiard`` in
:mod:`orpheus.derivations.continuous.trajectory_resolvent` (the
Birkhoff transfer-operator resolvent on a billiard); the 2nd is
``MomentSpace`` in
:mod:`orpheus.derivations.continuous.fn_method` (Galerkin half-range
projection on Legendre moments). With ``Spectrum`` landing the third
sibling, the math-heart pattern crosses the "≥3 instances" threshold:
three structurally-distinct mathematical attacks on the same
boundary-value transport problem, all sharing the same method-agnostic
:class:`~orpheus.geometry.structured_geometry.StructuredGeometry` /
:class:`~orpheus.derivations.common.solution_types.CriticalSolution`
contract.

The single class :class:`Spectrum` encapsulates everything the Case
singular-eigenfunction method needs to commit to before computing a
solution:

* **Geometry** — slab, sphere, or cylinder. The three differ in
  *which* spectral structure dominates the boundary problem:

  - **Slab** (Atalay 1997 Eq 46): linearly-anisotropic Case
    eigenvalue :math:`\nu_0 = i u_0` with reflection
    coefficient :math:`R \in [0, 1)`. Even-mode criticality.
  - **Sphere** (Atalay 1997 Eq 54): the *parity flip*
    :math:`\psi(x, \mu) = -\psi(-x, -\mu)` (Mitsis 1963)
    reduces sphere to slab odd-mode.  Same dispersion function
    :math:`\Lambda(\nu) = 1 - c\nu\,\mathrm{atanh}(1/\nu)`.
  - **Cylinder** (Westfall–Metcalf 1972, isotropic only): the
    Bessel-K kernel introduces extra spectral structure;
    needs the full Mitsis–WM Fredholm iteration with
    Mitsis–Zweifel singular subtraction. The dispersion
    function is medium-property identical to slab/sphere
    (V_se-cyl.1).

* **Materials** — the cross sections specifying :math:`c =
  (\Sigma_s + \nu\Sigma_f)/\Sigma_t` and (when the
  package supports it) the linear-anisotropy moment :math:`f_1`.

* **Boundary condition** — the reflection coefficient :math:`R`,
  read by :func:`_extract_R_refl` off the geometry's OUTER ``BC``
  (:attr:`StructuredGeometry.bcs` ``[-1]``): ``BC.vacuum`` for
  :math:`R = 0`, and ``BC("partial", params={"albedo": R})`` for any
  other :math:`R`.
  Vacuum (:math:`R = 0`) for bare-critical configurations;
  partial reflection :math:`R \in (0, 1)` for the Atalay
  reflected-slab / reflected-sphere benchmarks. The cylinder
  pillar ships ``R = 0`` only (bare); the Westfall–Metcalf
  reflected-cylinder extension is not yet shipped.

* **Number of modes / quadrature size** :math:`n` — Atalay's
  mpmath ``maxdegree`` for the half-range moments + the
  bracket-scan resolution. Westfall–Metcalf's GL grid order
  for the Fredholm iteration. The class collapses both into
  a single ``n_modes`` field (typical operating point: 8–24,
  depending on geometry and required accuracy).

The class is intentionally a thin computational specialisation of the
pure-geometry
:class:`~orpheus.geometry.structured_geometry.StructuredGeometry` for
the singular-eigenfunction family. The solve methods
(:meth:`Spectrum.solve_critical`, :meth:`Spectrum.solve_fixed_source`)
are thin wrappers around the existing function-level API in
:mod:`...slab.one_group`, :mod:`...sphere.one_group`, and
:mod:`...cylinder.one_group`. The function-level API stays as the
load-bearing implementation; ``Spectrum`` is the class-level facade
that owns the math-rich documentation and populates the cross-method
shared :class:`CriticalSolution` / :class:`FluxSolution` result types.

Architectural role
==================

``Spectrum`` is to ``singular_eigenfunction`` what:

- ``Billiard`` (Birkhoff transfer-operator resolvent on a
  bouncing-ray phase space) is to ``trajectory_resolvent``;
- ``MomentSpace`` (Galerkin half-range Legendre projection) is to
  ``fn_method``.

All three are method-specific computational specialisations of the
method-agnostic :class:`StructuredGeometry`. They answer the same
triplet of questions:

1. *"What is the critical configuration?"* — :meth:`solve_critical`
   returns a :class:`CriticalSolution`.
2. *"Given a configuration, what is the flux shape?"* —
   :meth:`solve_fixed_source` returns a :class:`FluxSolution`
   (currently delegates to the cylinder's
   :meth:`compute_scalar_flux`; slab + sphere flux reconstruction
   ships separately under :mod:`...fn_method.slab.flux_reconstruction`
   per the Path B / KLL discipline).
3. *"At a given fixed configuration, what is the eigenvalue?"* —
   inferred via :meth:`solve_critical` for criticality, or
   :meth:`solve_kinf` for the infinite-medium spectrum.

The shared :class:`CriticalSolution` / :class:`FluxSolution` types
are the load-bearing piece of the unification — they make Spectrum,
MomentSpace, Billiard substitutable at the cross-method comparison
boundary (``tests/cross_method/adapters.py``). A *behavioural*
Protocol over the math-heart classes was
tried and retired: the Phase-D ``TransportSolver`` (in
``orpheus.derivations.common.solver_protocol``) conflated continuous
reference generators with discrete production solvers, which have
functionally different roles, so the shared surface is the result
types plus direct
:class:`~orpheus.geometry.structured_geometry.StructuredGeometry`
construction — no Protocol.

References
----------

* Case, K.M. (1960) "Elementary Solutions of the Transport Equation
  and Their Applications", *Annals of Physics* **9**, 1-23. THE
  foundational paper.
* Case, K.M., Zweifel, P.F. (1967) *Linear Transport Theory*,
  Addison-Wesley. The canonical exposition; Chapter 4 covers
  half-range completeness.
* Mika, J.R. (1961) "The Initial-Value Problem in One-Velocity
  Transport Theory", *Nukleonik* **3**, 49–55. Proves Case's
  expansion-completeness theorem.
* Inönü, E. (1973) "Application of singular eigenfunctions in the
  finite slab problem", *Transp. Theor. Stat. Phys.* **3**, 137.
  Gives the X-function for finite media.
* Mitsis, G.F. (1963) "Transport Solutions to the Monoenergetic
  Critical Problems", ANL-6787. Sphere via parity flip.
* Westfall, R.M., Metcalf, D.R. (1973) "Singular Eigenfunction
  Solution of the Monoenergetic Neutron Transport Equation for
  Finite Radially Reflected Critical Cylinders", *Nucl. Sci. Eng.*
  **52**, 1-11. Cylinder via Bessel-K + addition theorem.
* Atalay, M.A. (1997) "The reflected slab and sphere criticality
  problem with anisotropic scattering in one-speed neutron
  transport theory", *Prog. Nucl. Energy* **31**(3), 229-252.
  Linear anisotropy + half-range bi-orthogonality + first-order
  Fredholm iteration.
* :doc:`/theory/references/singular_eigenfunction` — the canonical theory
  exposition; the "Case spectrum and the expansion theorem" section
  is the rich-narrative companion to this module's docstring.
* :doc:`/theory/references/index` § :ref:`reference-solvers-three-meanings`
  — locates Spectrum (γ: singular-eigenfunction angular Green's),
  Billiard (α: trajectory resolvent), and MomentSpace (γ: F_N is also
  in the singular-eigenfunction pillar but uses a *different*
  collocation projection) within the same Green's-function landscape.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import numpy as np

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.solution_types import (
    CriticalSolution,
    FluxSolution,
)
from orpheus.geometry.structured_geometry import StructuredGeometry


__all__ = ["Spectrum"]


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------


def _extract_R_refl(geom: StructuredGeometry) -> float:
    r"""Extract the reflection coefficient :math:`R` from the geometry's
    outer BC.

    Reads :attr:`StructuredGeometry.bcs` ``[-1]`` (the outer endpoint;
    for SPH / CYL it is the only entry, for SLB it is the right-hand
    surface).

    Convention:

    * ``BC.vacuum``  → :math:`R = 0` (bare).
    * ``BC.reflective`` → :math:`R = 1`, which the singular-
      eigenfunction solvers reject (the slab thickness / sphere
      radius drops out of the criticality condition under perfect
      reflection — Atalay omits :math:`R = 1` columns from Tables
      2-5 for this reason).
    * Any other ``BC`` with ``params={"albedo": R}`` → :math:`R` from
      params. The custom partial-reflection BC is the entry point
      for Atalay reflected-slab / reflected-sphere benchmarks.
    """
    bc = geom.bcs[-1]  # outer boundary; SLB → right, SPH/CYL → only entry
    if bc.kind == "vacuum":
        return 0.0
    if bc.kind == "reflective":
        return 1.0  # caller will get a ValueError from the solver
    # Custom BC: read albedo from params.
    if "albedo" in bc.params:
        return float(bc.params["albedo"])
    raise ValueError(
        f"Spectrum: cannot extract reflection coefficient from "
        f"outer BC={bc!r}. Use BC.vacuum (R=0), or "
        f'BC("partial", params={{"albedo": R}}) with R in [0, 1).'
    )


def _extract_f1(mixture: Mixture) -> float:
    r"""Extract the linear-anisotropy mean cosine :math:`f_1` from a Mixture.

    Convention:

    * If ``mixture.SigS`` has only :math:`\Sigma_{s,0}` (P_0
      isotropic), :math:`f_1 = 0`.
    * If ``mixture.SigS`` has :math:`\Sigma_{s,1}` (P_1 anisotropic),
      :math:`f_1 = \Sigma_{s,1} / \Sigma_{s,0}` for the 1G case
      (mean cosine of scattering, the Atalay 1997 :math:`f_1`).

    The Atalay validity bound :math:`c \le 1 + 1/(3 f_1)` (Eq 5)
    is checked at the solver-entry layer, not here.

    The cylinder pillar (Westfall–Metcalf) is **isotropic only**:
    if a non-zero :math:`f_1` is detected for a cylinder geometry,
    the construction fails fast.

    Returns
    -------
    f1 : float
        The linear-anisotropy mean cosine, in [0, 1).

    Raises
    ------
    NotImplementedError
        If the mixture has higher-than-P_1 moments (P_2+) — only
        isotropic and linearly-anisotropic scattering are in pillar.
    """
    if not mixture.SigS:
        raise ValueError("Mixture has no scattering matrices (SigS is empty)")
    if len(mixture.SigS) == 1:
        return 0.0
    if len(mixture.SigS) == 2:
        sig_s_p0 = mixture.SigS[0].toarray().astype(float)
        sig_s_p1 = mixture.SigS[1].toarray().astype(float)
        # 1G: f1 = SigS1 / SigS0. For multi-group this would need
        # careful handling — currently the singular-eigenfunction
        # pillar is 1G only above the function level.
        if sig_s_p0.shape != (1, 1):
            raise NotImplementedError(
                f"Spectrum: linear-anisotropy extraction "
                f"is currently 1G-only (got {sig_s_p0.shape[0]}G mixture). "
                f"Multi-group anisotropic singular-eigenfunction is not "
                f"yet shipped."
            )
        sigma_s_0 = float(sig_s_p0[0, 0])
        if sigma_s_0 <= 0.0:
            return 0.0
        return float(sig_s_p1[0, 0]) / sigma_s_0
    raise NotImplementedError(
        f"Spectrum: P_{len(mixture.SigS) - 1} scattering "
        f"is out of pillar (only P_0 isotropic and P_1 linearly-"
        f"anisotropic are in scope)."
    )


# ---------------------------------------------------------------------
# Spectrum class
# ---------------------------------------------------------------------


@dataclass(frozen=True)
class Spectrum:
    r"""The transport operator's singular eigenfunction spectrum (Case 1960).

    The *singular eigenfunction expansion* (Case 1960; Case-Zweifel
    1967) decomposes the angular flux :math:`\psi(r, \mu)` into
    eigenfunctions of the transport operator's resolvent. The
    mathematical setup, in seven parts:

    1. **The transport operator's eigenvalue problem.**  For the
       homogeneous one-speed transport equation
       :math:`\mu\,\partial_x \psi + \Sigma_t \psi = c\,\Sigma_t\,
       \bar{\psi}`, seek separated solutions of the form
       :math:`\psi(x, \mu) = \phi_\nu(\mu)\,e^{-\Sigma_t x / \nu}`.
       Substitution and division by the exponential gives the
       *Case-eigenfunction equation*

       .. math::

           (\nu - \mu)\,\phi_\nu(\mu) = c\,\nu \int_{-1}^{1}
                                        \phi_\nu(\mu')\,d\mu'.

       The function :math:`\phi_\nu(\mu)` is the **singular
       eigenfunction** at *eigenvalue* :math:`\nu` of the transport
       operator's spectrum.

    2. **The spectrum decomposes into two parts**:

       * **Discrete eigenvalues** :math:`\pm\nu_0` — the diffusion-like
         modes; :math:`\nu_0` satisfies the Case dispersion relation

         .. math::

             \Lambda(\nu_0) := 1 - c\,\nu_0\,\mathrm{atanh}(1/\nu_0) = 0.

         For :math:`c < 1` (sub-multiplying), :math:`\nu_0 > 1` is
         real (a slow exponential decay); the transport operator
         describes diffusing-and-decaying neutrons. For :math:`c > 1`
         (super-multiplying / fissioning), :math:`\nu_0 = i u_0` is
         purely imaginary (oscillatory diffusion-like mode); the
         operator has growing modes whose criticality condition is
         what reactor physics computes. The *medium-property*
         identity of :math:`\Lambda(\nu)` across slab, sphere, and
         cylinder geometries — verified at V_se-cyl.1 — is the
         load-bearing reason the dispersion-root primitive
         :func:`...fn_method.core.dispersion.case_nu0` is shared
         between :mod:`...fn_method` and
         :mod:`...singular_eigenfunction` below the trusted-library
         line. **The eigenvalue is a property of the medium, not the
         geometry.**

       * **Continuum** :math:`\nu \in [-1, 1]` — the streaming-like
         modes that handle short-range / boundary-layer behaviour.
         For :math:`\nu \in (-1, 1)` and :math:`\mu = \nu`, the
         eigenfunction equation has a *singular*
         :math:`(\nu - \mu)^{-1}` factor — hence "singular
         eigenfunctions". The correct treatment is via the
         distribution

         .. math::

             \phi_\nu(\mu) = \mathrm{P} \frac{c\,\nu}{2 (\nu - \mu)}
                             + \lambda(\nu)\,\delta(\nu - \mu),

         where :math:`\mathrm{P}` is the Cauchy principal value and
         :math:`\lambda(\nu) = 1 - c\,\nu\,\mathrm{atanh}(\nu)` is
         the dispersion function on the continuum. Both pieces — the
         principal-value singularity and the delta-function residue —
         are essential. Truncating either gives a manifestly wrong
         answer (see V_se-cyl.8 for the canonical Mitsis-Zweifel
         singular-subtraction structural identity).

       The full spectrum is :math:`\Sigma = \{\pm\nu_0\} \cup
       [-1, 1]`. This is *the spectrum*; the data class
       :class:`Spectrum` carries the medium and geometry that locate
       a particular instance of it.

    3. **The expansion theorem** (Case 1960 Theorem 1; Mika 1961
       completeness): any admissible angular flux on the
       homogeneous problem can be written as

       .. math::

           \psi(x, \mu) = a_+\,\phi_{\nu_0}(\mu)\,e^{-x/\nu_0}
                        + a_-\,\phi_{-\nu_0}(\mu)\,e^{x/\nu_0}
                        + \int_{-1}^{1} A(\nu)\,\phi_\nu(\mu)\,
                          e^{-x/\nu}\,d\nu

       with three sets of expansion coefficients:

       * :math:`a_+` (dominant outgoing diffusion mode),
       * :math:`a_-` (incoming diffusion mode — vanishes for
         physical solutions on a half-space, but is non-zero
         on a finite slab),
       * :math:`A(\nu)` (continuum density on :math:`[-1, 1]` —
         the analog of a Fourier amplitude for the continuum
         eigenfunctions).

       The expansion is **complete** (Mika 1961) and **unique** —
       a fact whose proof needs the half-range completeness theorem
       below. The expansion theorem is the singular-eigenfunction
       counterpart of Galerkin / Legendre projection in the F_N
       method; in both cases an angular flux is decomposed onto a
       spectrum of operator eigenfunctions, but in F_N the spectrum
       is *truncated* to a finite Legendre basis whereas here it is
       *complete*. This is what makes singular-eigenfunction
       criticality conditions exact (modulo the dispersion-root
       computation), where F_N introduces a truncation error
       :math:`O(N^{-p})`.

    4. **Half-range completeness** (Case-Zweifel 1967 Ch. 4;
       Inönü 1973). At a boundary surface, the angular flux's
       relevant function space splits into incoming
       (:math:`\mu > 0`) and outgoing (:math:`\mu < 0`) half-ranges.
       The X-function of Inönü 1973,

       .. math::

           X(\mu) = \exp\!\left[
               \frac{c\mu}{2}
               \int_0^1 \frac{\mathrm{atanh}(\mu')\,d\mu'}
                              {(\mu' - \mu)(c\mu' \mathrm{atanh}(\mu') - 1)}
               \right],

       is the projection operator that **completes** the half-range
       basis: any half-range function admits a unique expansion onto
       :math:`\{X(\mu)^{-1}\,\phi_{\nu_0}(\mu),\, X(\mu)^{-1}\,
       \phi_\nu(\mu) : \nu \in (0, 1)\}`. The X-function carries
       the medium dependence (the integrand depends on :math:`c`)
       and is what makes finite-geometry boundary conditions
       solvable in closed form. ORPHEUS implements :math:`X(\mu)` in
       :func:`...fn_method.core.x_function.x_function_atalay`
       (with the ERR-037 :math:`\mu = \tanh(t)` pole-cancellation
       substitution; see Numerical Bug Signatures § Signature 7 for
       the canonical instance of why this matters).

    5. **Geometry-specific reductions.** The expansion theorem +
       half-range completeness give the **machinery**; the geometry
       supplies the **boundary conditions** that pin the expansion
       coefficients :math:`(a_+, a_-, A(\nu))`. Each geometry the
       :class:`Spectrum` class supports is one such reduction:

       * **Slab** (``geometry == "slab"``) — Atalay 1997 Eq. 46.
         The X-function on :math:`(0, 1]` determines the boundary
         coefficients. The *even-mode* criticality condition for
         a symmetric slab of half-thickness :math:`d` and
         specular-style reflection coefficient :math:`R` is the
         scalar arctan-equation root in :math:`d`. The same
         spectrum :math:`\Sigma = \{\pm\nu_0\} \cup [-1, 1]`
         applies regardless of slab width; the criticality
         condition factors through the dispersion relation
         (medium) and the X-function residues (boundary).
       * **Sphere** (``geometry == "sphere"``) — Atalay 1997 Eq.
         54, derived via the *parity flip* trick (Mitsis 1963;
         Atalay Eq. 47): the antisymmetric BC :math:`\psi(x, \mu)
         = -\psi(-x, -\mu)` reduces sphere to the *odd-mode*
         counterpart of the slab problem. The structural changes
         from slab to sphere are surgical: the kernel
         :math:`T(R, \mu)` is replaced by :math:`T_1(R, \mu)`
         (a sign flip in the second exponential), the K-moments
         become L-moments (same kernel structure, different sign),
         and the LHS criticality term swaps :math:`\sin \leftrightarrow
         \cos` and :math:`R \to -R`. **The discrete eigenvalue
         :math:`\nu_0` and the continuum on :math:`[-1, 1]` are
         identical** — the same ``Spectrum`` instance, different
         boundary projection.
       * **Cylinder** (``geometry == "cylinder"``) — Westfall–Metcalf
         1972 (WM-72) for **isotropic scattering only**. The
         radial Bessel-K kernel from the cylindrical addition
         theorem introduces additional spectral structure beyond
         the basic Case eigenfunctions: the Bessel modes
         :math:`I_0(R/\nu)` couple radial and angular dimensions
         in a way the slab/sphere geometries don't. Resolution
         needs the **Mitsis–WM Fredholm iteration** (WM-72 Eqs
         30–32) with **Mitsis-Zweifel singular subtraction**
         (V_se-cyl.8) to handle the principal-value + delta
         residue cleanly. Solver achieves 6–7 digit accuracy at
         :math:`n_\text{grid} = 24` (matching WM-72 Table II's
         own quoted precision).

    6. **Linear anisotropy** (Atalay 1997, slab + sphere).
       Replacing isotropic scattering :math:`c \int \psi\,d\mu` by
       linearly-anisotropic :math:`c \int \psi (1 + 3 f_1
       \mu\mu')\,d\mu` — where :math:`f_1` is the mean cosine of
       scattering — adds a single extra term to the dispersion
       relation but preserves the spectrum's structure: still one
       discrete pair :math:`\pm \nu_0` plus a continuum
       on :math:`[-1, 1]`. The new spectrum is exactly the same
       picture; the operator is just slightly different. Atalay's
       :math:`K_j` and :math:`L_j` half-range moments
       (V_se-cyl Atalay derivations) integrate the new structure
       explicitly. The validity condition :math:`c \le 1 + 1/(3 f_1)`
       (Atalay Eq 5) bounds the regime where the transport operator
       has only one pair of discrete modes; outside that band
       complex-conjugate eigenvalue pairs appear (Dahl-Sjöstrand
       1979) and Atalay's first-order Fredholm iteration fails to
       detect them. The cylinder pillar (WM-72) is **isotropic only**;
       linearly-anisotropic cylinder is research-grade and not in
       the package.

    7. **Connection to Billiard, MomentSpace, and the
       three-meanings taxonomy.**  Three separate ways to compute
       the same Green's function on the same physical problem,
       differing in *which mathematical object* you decompose:

       * ``Spectrum`` (this class): decompose the **operator
         spectrum** directly. Locates :math:`\nu_0` from the
         dispersion relation, projects boundary data via the
         X-function, integrates over the continuum. The
         eigenfunctions :math:`\phi_\nu(\mu)` are *operator
         eigenfunctions*; the Green's function is the resolvent
         :math:`(L - \nu I)^{-1}` evaluated at boundary data.
       * ``Billiard`` (:mod:`...trajectory_resolvent`): trace **rays
         through phase space**. Each bouncing trajectory contributes
         :math:`\alpha^n e^{-n\tau}` (Birkhoff transfer-operator
         resolvent on the billiard table). The Green's function is
         the *path-integral* sum :math:`\sum_n S^n = (I - S)^{-1}`.
       * ``MomentSpace`` (:mod:`...fn_method`): project the
         **boundary angular flux** onto a finite Legendre moment
         basis, collocate. The Green's function is implicit in the
         determinant condition :math:`\det M = 0`; truncation is
         to :math:`N+1` moments.

       The Sanchez–Chandrasekhar three-meanings taxonomy
       (:doc:`/theory/references/index` § :ref:`reference-solvers-three-meanings`)
       locates ``Spectrum`` under meaning **(γ): singular-eigenfunction
       angular Green's function** — directly construct
       :math:`G(\tau, \tau'; \mu, \mu')` as a sum over ν-spectrum
       eigenfunctions weighted by X-function residues. ``Billiard``
       is meaning **(α): trajectory resolvent**. ``MomentSpace``
       is also under (γ) but uses a different (Galerkin-projection)
       *closure*. Cross-checks between the three pillars are L1
       evidence per ``vv-principles`` § "structural independence"
       — the three integrands are structurally distinct (ν-spectrum
       vs ray-traced phase space vs Legendre moments), so triple
       agreement is not coincidence.

    Why this is **not** a Protocol implementation
    ============================================

    Per the project's "unify after two instances" memory, a
    behavioural ``TransportSolver`` Protocol (in
    ``orpheus.derivations.common.solver_protocol``) was designed only
    AFTER ``MomentSpace`` and ``Billiard`` worked, with ``Spectrum``
    as the third instance meant to validate it empirically. The
    third instance instead **falsified** it: Phase D deleted the
    Protocol, because it conflated continuous reference generators
    (which consume a
    :class:`~orpheus.geometry.structured_geometry.StructuredGeometry`
    directly through a frozen ``__init__``) with discrete production
    solvers (which consume ``(materials, mesh, params)`` through the
    canonical free functions) — two functionally different roles that
    a single Protocol could only blur.

    What survives is the *structural* contract, and it is the
    stronger one: every math-heart class takes the same
    ``geometry`` / ``materials`` construction pair and returns the
    same :class:`CriticalSolution` / :class:`FluxSolution` types, so
    the cross-method adapters can hold any of them without a
    nominal type — and ``Spectrum`` inherits from no ABC or Protocol
    at all.

    Construction (Phase D)
    ----------------------

    Direct construction with a :class:`StructuredGeometry` and
    a ``materials: dict[int, Mixture]`` payload::

        from orpheus.geometry.structured_geometry import (
            Region, StructuredGeometry,
        )
        from orpheus.geometry.mesh import BC

        geom = StructuredGeometry(
            geometry="SLB",
            regions=(Region(mat_id=0, outer_thickness_cm=L_full),),
            bcs=(BC.vacuum, BC.vacuum),
        )
        spec = Spectrum(geometry=geom, materials={0: mix})
        sol = spec.solve_critical()

    Parameters
    ----------
    geometry : :class:`StructuredGeometry`
        Pure-geometry layer object. Tag MUST be ``"SLB"`` / ``"SPH"``
        / ``"CYL"``. Singular-eigenfunction criticality requires a
        finite spatial domain (no infinite-medium tag).
    materials : dict[int, Mixture]
        Production-protocol materials, keyed by material ID. The
        single-region geometry's mat_id selects the active mixture.
    n_modes : int
        Quadrature size for the half-range moments (Atalay slab /
        sphere) or the Mitsis-WM Fredholm grid (cylinder). Defaults
        to 8.
    """

    geometry: StructuredGeometry
    materials: dict[int, Mixture]
    n_modes: int = 8

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __post_init__(self) -> None:
        r"""Validate the singular-eigenfunction method's structural preconditions.

        The package as shipped applies to:

        * **Slab** (Atalay 1997 with linear anisotropy, reflected /
          partial-reflective boundary).
        * **Sphere** (Atalay 1997 via parity flip).
        * **Cylinder** (Westfall–Metcalf 1972, isotropic only,
          bare-critical only).
        """
        if self.geometry.geometry not in {"SLB", "SPH", "CYL"}:
            raise ValueError(
                f"Spectrum supports geometry ∈ {{SLB, SPH, CYL}}, "
                f"got {self.geometry.geometry!r}. Infinite-medium k_inf "
                f"is out of pillar — use MomentSpace.solve_kinf for "
                f"k_inf computations."
            )
        if self.n_modes < 2:
            raise ValueError(f"n_modes must be ≥ 2, got {self.n_modes}")
        if self._mat_id not in self.materials:
            raise ValueError(
                f"materials dict missing mat_id={self._mat_id} "
                f"required by geometry; got keys "
                f"{sorted(self.materials.keys())}"
            )
        # Eager validation: cylinder + anisotropic mixture is a
        # construction-time error (the per-method solver would raise
        # at solve time, but we want the failure surfaced at the
        # facade boundary so callers know the geometry/material
        # combination is out of pillar).
        if self.geometry.geometry == "CYL":
            mix = self.materials.get(self._mat_id)
            if mix is not None and len(mix.SigS) > 1:
                sig_s_p1 = mix.SigS[1].toarray().astype(float)
                if np.any(np.abs(sig_s_p1) > 0.0):
                    raise NotImplementedError(
                        "Spectrum: cylinder + linear "
                        "anisotropy is out of pillar (Westfall-Metcalf "
                        "1972 covers isotropic only). Use slab or sphere "
                        "for anisotropic problems."
                    )

    # ------------------------------------------------------------------
    # Material accessors
    # ------------------------------------------------------------------

    @property
    def _mat_id(self) -> int:
        """Active mat_id — the single region's material identifier."""
        return self.geometry.regions[0].mat_id

    @property
    def _mixture(self) -> Mixture:
        """The active :class:`Mixture` for this geometry's mat_id."""
        return self.materials[self._mat_id]

    # ------------------------------------------------------------------
    # Derived primary parameters
    # ------------------------------------------------------------------

    @property
    def c(self) -> float:
        r"""Mean number of secondaries per collision, :math:`c`.

        For 1G isotropic-scattering problems
        :math:`c = (\Sigma_s + \nu\Sigma_f)/\Sigma_t`.

        Raises
        ------
        ValueError
            If the active mixture has more than one energy group.
        """
        mixture = self._mixture
        sig_t = np.asarray(mixture.SigT, dtype=float)
        if sig_t.shape[0] != 1:
            raise ValueError(
                f"Spectrum.c requires a 1G mixture (got "
                f"{sig_t.shape[0]}G). Multi-group singular-eigenfunction "
                f"is out of pillar."
            )
        sig_s_p0 = mixture.SigS[0].toarray().astype(float)
        nu_sig_f = np.asarray(mixture.SigP, dtype=float)
        return float((sig_s_p0[0, 0] + nu_sig_f[0]) / sig_t[0])

    @property
    def f1(self) -> float:
        r"""Linear-anisotropy mean cosine :math:`f_1`.

        Extracted from :math:`\Sigma_{s,1} / \Sigma_{s,0}` of the active
        mixture.
        """
        return _extract_f1(self._mixture)

    @property
    def n_groups(self) -> int:
        """Number of energy groups in the active mixture."""
        return int(np.asarray(self._mixture.SigT).shape[0])

    # ------------------------------------------------------------------
    # solve_critical
    # ------------------------------------------------------------------

    def solve_critical(
        self,
        *,
        n_bracket: Optional[int] = None,
        bisect_tol: float = 1e-10,
        max_bisect: int = 80,
        mode: int = 1,
        radius_min: Optional[float] = None,
        radius_max: Optional[float] = None,
    ) -> CriticalSolution:
        r"""Solve the critical configuration via singular-eigenfunction expansion.

        For a **bare or partially-reflective slab** (``geometry ==
        "slab"``) this returns the critical half-thickness :math:`d`
        in mean free paths via Atalay 1997 Eq 46.
        ``parameter_kind == "half_thickness_mfp"``.

        For a **bare or partially-reflective sphere** (``geometry ==
        "sphere"``) this returns the critical radius :math:`R_c` in
        mean free paths via Atalay 1997 Eq 54.
        ``parameter_kind == "radius_mfp"``.

        For a **bare cylinder** (``geometry == "cylinder"``,
        ``BC.vacuum`` outer) this returns the critical radius
        :math:`R_c` via Westfall-Metcalf 1973 Eq 32 — the full
        Mitsis-WM Fredholm iteration with Mitsis-Zweifel singular
        subtraction. Cylinder ships isotropic-only and bare-only.

        The underlying solve delegates to:

        * :func:`...slab.one_group.solve_case_method_slab_critical`
          for slab,
        * :func:`...sphere.one_group.solve_case_method_sphere_critical`
          for sphere,
        * :func:`...cylinder.one_group.solve_singular_eigenfunction_cylinder_bare_critical`
          for cylinder.

        The rich method-specific result
        (:class:`CaseMethodSlabResult` /
        :class:`CaseMethodSphereResult` /
        :class:`CylinderSingularEigenfunctionResult`) is preserved in
        :attr:`CriticalSolution.metadata` under the key
        ``"raw_result"`` for callers that need access to the
        :math:`\nu_0` dispersion root, the K_j / L_j moments, or the
        Fredholm-iteration continuum amplitudes :math:`A'(\nu)`.

        Parameters
        ----------
        n_bracket : int | None, default None
            Bracket-scan resolution for slab / sphere root-finding.
            ``None`` means use the underlying solver's default (200
            for slab, 200 for sphere).  Cylinder uses a fixed
            64-point sweep regardless of this value.
        bisect_tol : float, default 1e-10
            Bisection tolerance on the configuration parameter.
        max_bisect : int, default 80
            Maximum bisection iterations.
        mode : int, default 1
            Eigenvalue mode for slab / sphere (1 = fundamental, 2 =
            first overtone, etc).  Ignored for cylinder (the WM-72
            sweep uses the first sign change).
        radius_min, radius_max : float | None, default None
            Bracket bounds for sphere / cylinder radius search. ``None``
            uses the underlying solver's default (slab: ``[0.05, 60]``;
            sphere: ``[0.10, 80]``; cylinder: ``[0.1, 30]``).

        Returns
        -------
        :class:`CriticalSolution`

        Notes
        -----
        Bit-equality with the function-level API: this method MUST
        produce the SAME float results as a direct call to the
        underlying ``solve_*`` function with matched kwargs. The
        class-level call is a thin facade. Verified by
        :mod:`tests.derivations.test_singular_eigenfunction_spectrum`
        (the foundation gate that pins the bit-equality invariant).
        """
        tag = self.geometry.geometry
        mixture = self._mixture
        sig_t = np.asarray(mixture.SigT, dtype=float)
        n_groups = sig_t.shape[0]

        if n_groups != 1:
            raise NotImplementedError(
                "Spectrum.solve_critical is currently 1G-only. "
                "The singular-eigenfunction multi-group spatial "
                "extension is not yet shipped."
            )

        sig_s_p0 = mixture.SigS[0].toarray().astype(float)
        nu_sig_f = np.asarray(mixture.SigP, dtype=float)
        c = float((sig_s_p0[0, 0] + nu_sig_f[0]) / sig_t[0])
        if c <= 1.0:
            raise ValueError(
                f"Singular-eigenfunction bare-critical {tag} requires c > 1 "
                f"(multiplying medium); got c={c}."
            )

        if tag == "SLB":
            return self._solve_critical_slab(
                c, n_bracket, bisect_tol, max_bisect, mode
            )
        if tag == "SPH":
            return self._solve_critical_sphere(
                c, n_bracket, bisect_tol, max_bisect, mode, radius_min, radius_max
            )
        if tag == "CYL":
            return self._solve_critical_cylinder(
                c, bisect_tol, radius_min, radius_max
            )
        raise NotImplementedError(  # pragma: no cover
            f"Spectrum.solve_critical: unhandled geometry {tag!r}"
        )

    def _solve_critical_slab(
        self,
        c: float,
        n_bracket: Optional[int],
        bisect_tol: float,
        max_bisect: int,
        mode: int,
    ) -> CriticalSolution:
        from .slab.one_group import solve_case_method_slab_critical

        R_refl = _extract_R_refl(self.geometry)
        f1 = self.f1

        kwargs: dict[str, Any] = {
            "c": c,
            "R": R_refl,
            "f1": f1,
            "bisect_tol": bisect_tol,
            "max_bisect": max_bisect,
            "mode": mode,
            "maxdegree": self.n_modes,
        }
        if n_bracket is not None:
            kwargs["n_bracket"] = n_bracket
        res = solve_case_method_slab_critical(**kwargs)
        return CriticalSolution(
            eigenvalue=1.0,
            eigenvalue_kind="k_eff",
            parameter_value=float(res.d_critical_mfp),
            parameter_kind="half_thickness_mfp",
            converged=True,
            metadata={
                "n_groups": 1,
                "method": "solve_case_method_slab_critical",
                "c": float(res.c),
                "R": float(res.R),
                "f1": float(res.f1),
                "u0": float(res.u0),
                "nu_bar": float(res.nu_bar),
                "z0": float(res.z0),
                "K_moments": dict(res.K_moments),
                "eq46_residual": float(res.eq46_residual),
                "n_bracket_points": int(res.n_bracket_points),
                "n_bisect_iters": int(res.n_bisect_iters),
                "raw_result": res,
            },
        )

    def _solve_critical_sphere(
        self,
        c: float,
        n_bracket: Optional[int],
        bisect_tol: float,
        max_bisect: int,
        mode: int,
        radius_min: Optional[float],
        radius_max: Optional[float],
    ) -> CriticalSolution:
        from .sphere.one_group import solve_case_method_sphere_critical

        R_refl = _extract_R_refl(self.geometry)
        f1 = self.f1

        kwargs: dict[str, Any] = {
            "c": c,
            "R_refl": R_refl,
            "f1": f1,
            "bisect_tol": bisect_tol,
            "max_bisect": max_bisect,
            "mode": mode,
        }
        if n_bracket is not None:
            kwargs["n_bracket"] = n_bracket
        if radius_min is not None:
            kwargs["radius_min"] = radius_min
        if radius_max is not None:
            kwargs["radius_max"] = radius_max
        res = solve_case_method_sphere_critical(**kwargs)
        return CriticalSolution(
            eigenvalue=1.0,
            eigenvalue_kind="k_eff",
            parameter_value=float(res.R_critical_mfp),
            parameter_kind="radius_mfp",
            converged=True,
            metadata={
                "n_groups": 1,
                "method": "solve_case_method_sphere_critical",
                "c": float(res.c),
                "R_refl": float(res.R_refl),
                "f1": float(res.f1),
                "u0": float(res.u0),
                "nu_bar": float(res.nu_bar),
                "z0": float(res.z0),
                "L_moments": dict(res.L_moments),
                "eq54_residual": float(res.eq54_residual),
                "n_bracket_points": int(res.n_bracket_points),
                "n_bisect_iters": int(res.n_bisect_iters),
                "raw_result": res,
            },
        )

    def _solve_critical_cylinder(
        self,
        c: float,
        bisect_tol: float,
        radius_min: Optional[float],
        radius_max: Optional[float],
    ) -> CriticalSolution:
        from .cylinder.one_group import (
            solve_singular_eigenfunction_cylinder_bare_critical,
        )

        R_refl = _extract_R_refl(self.geometry)
        if R_refl != 0.0:
            raise NotImplementedError(
                f"Spectrum cylinder solve_critical: only bare cylinder "
                f"(R_refl=0) is in pillar (Westfall-Metcalf 1972 limit); "
                f"got R_refl={R_refl}."
            )

        kwargs: dict[str, Any] = {
            "c": c,
            "n_grid": self.n_modes,
            "bisect_tol": bisect_tol,
        }
        if radius_min is not None:
            kwargs["R_min"] = radius_min
        if radius_max is not None:
            kwargs["R_max"] = radius_max
        res = solve_singular_eigenfunction_cylinder_bare_critical(**kwargs)
        return CriticalSolution(
            eigenvalue=1.0,
            eigenvalue_kind="k_eff",
            parameter_value=float(res.r_c_mfp),
            parameter_kind="radius_mfp",
            converged=bool(res.converged),
            metadata={
                "n_groups": 1,
                "method": (
                    "solve_singular_eigenfunction_cylinder_bare_critical"
                ),
                "c": float(res.c),
                "u_0": float(res.u_0),
                "nu_0": complex(res.nu_0),
                "criticality_residual": float(res.criticality_residual),
                "iterations": int(res.iterations),
                "n_grid": int(res.n_grid),
                "A_prime": np.asarray(res.A_prime).copy(),
                "raw_result": res,
            },
        )

    # ------------------------------------------------------------------
    # solve_fixed_source
    # ------------------------------------------------------------------

    def solve_fixed_source(
        self,
        q: Optional[np.ndarray] = None,
        *,
        n_eval: int = 64,
    ) -> FluxSolution:
        r"""Reconstruct the flux at the critical configuration.

        For the singular-eigenfunction pillar, "fixed source" is a
        slight misnomer: the canonical use is
        **flux reconstruction at the critical configuration**, where
        the source :math:`q` is the implicit fission source that
        sustains criticality. The semi-physical interpretation is
        that the converged flux is the dominant Case eigenfunction
        at the critical eigenvalue :math:`k_\text{eff} = 1`, and the
        ``q`` argument is preserved in the API for future extensions
        (e.g., subcritical fixed-source problems via the Case
        spectrum's continuum integral).

        Currently implemented for **cylinder only** — Westfall-Metcalf
        Eq 22 / V_se-cyl.7 reconstruction:

        .. math::

            \rho(r) = J_0(r/u_0)

        — the dominant Case eigenfunction with imaginary :math:`\nu_0
        = i u_0` for :math:`c > 1`. Slab + sphere flux reconstruction
        is **out of pillar** for the singular-eigenfunction package
        (the F_N / KLL Fredholm iteration in
        :mod:`...fn_method.slab.flux_reconstruction` /
        :mod:`...fn_method.sphere.flux_reconstruction` is the
        structurally-independent flux reconstruction; using it here
        would violate the structural-independence rule).

        Parameters
        ----------
        q : np.ndarray, optional
            Source distribution (currently ignored; reserved for
            future fixed-source extensions). At the critical
            configuration the flux is the dominant Case eigenfunction
            and the source is implicit (the fission source).
        n_eval : int, default 64
            Number of spatial evaluation points.

        Returns
        -------
        :class:`FluxSolution`

        Raises
        ------
        NotImplementedError
            For slab / sphere geometries — flux reconstruction in
            those is owned by the F_N pillar
            (:mod:`...fn_method.slab.flux_reconstruction` /
            :mod:`...fn_method.sphere.flux_reconstruction`).
        """
        del q  # reserved for future use
        tag = self.geometry.geometry
        if tag != "CYL":
            raise NotImplementedError(
                f"Spectrum.solve_fixed_source: flux reconstruction for "
                f"{tag!r} is owned by the F_N pillar "
                f"(orpheus.derivations.continuous.fn_method) — using it "
                f"here would violate the structural-independence rule. "
                f"Use MomentSpace.reconstruct_flux for slab / sphere."
            )

        # Cylinder: solve criticality first, then reconstruct.
        critical = self.solve_critical()
        raw = critical.metadata["raw_result"]
        r_c = float(raw.r_c_mfp)
        r_eval = np.linspace(0.0, r_c, n_eval)
        rho = np.asarray(raw.compute_scalar_flux(r_eval), dtype=float).ravel()
        return FluxSolution(
            spatial_nodes=r_eval,
            scalar_flux=rho,
            angular_flux=None,
            angular_nodes=None,
            eigenvalue=critical.eigenvalue,
            eigenvalue_kind=critical.eigenvalue_kind,
            spatial_units="mfp",
            metadata={
                "geometry": "cylinder",
                "method": "wm72_bare_cylinder_J0",
                "u_0": float(raw.u_0),
                "r_c_mfp": r_c,
                "c": critical.metadata.get("c"),
            },
        )
