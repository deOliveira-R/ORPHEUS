.. _theory-homogeneous:

=============================================
Homogeneous Infinite-Medium Reactor
=============================================

.. contents:: Contents
   :local:
   :depth: 3


Key Facts
=========

**Read this before modifying the homogeneous solver.**

- **The infinite medium is the point in phase space.** ORPHEUS poses it
  on energy alone, the point in position and in direction: it holds one
  material and no geometry, since a geometry admits spatial dimension,
  more than one material and a direction chart. Position drops out by the
  translation symmetry of the medium and its data, direction by their
  rotation symmetry, and the reduction is exact for everything the
  definition admits (every
  :class:`~orpheus.data.macro_xs.mixture.Mixture` is isotropic, and every
  admitted question prefers no direction). A uniform problem whose source
  prefers a direction is not this object: it is the spatial marginal of a
  body problem (:ref:`infinite-medium-point-in-phase-space`)
- Balance: :math:`\mathbf{A}\phi = \frac{1}{k}\mathbf{F}\phi` where the loss
  matrix is :math:`\mathbf{A} = \text{diag}(\Sigma_t) - \Sigma_{s0}^T - 2\Sigma_2^T`
  and the production dyad is :math:`\mathbf{F} = \chi \otimes (\nu\Sigma_f)`
- **(n,2n) convention**: the :math:`(n,2n)` reaction is a **loss-side
  multiplicity-2 transfer** — it lives ONLY in :math:`\mathbf{A}` (as
  :math:`-2\Sigma_2^T`), NEVER in the production :math:`\mathbf{F}`. The two
  emitted neutrons are redistributed by :math:`2\Sigma_2`; they are not
  produced with the fission spectrum :math:`\chi`. Production is
  :math:`\nu\Sigma_f` only. (Double-counting :math:`2\,\text{colsum}(\Sigma_2)`
  into production — the retired bespoke bug — moves :math:`\kinf` by
  :math:`\sim 0.43` on the asymmetric-:math:`\Sigma_2` ``homo_2eg_n2n`` case.)
- 1-group: :math:`k = \nu\Sigma_f / \Sigma_a` (exact, no iteration)
- Multi-group: :math:`\kinf = \langle \nu\Sigma_f,\, \mathbf{A}^{-1}\chi \rangle`
  with flux :math:`\boldsymbol{\phi} \propto \mathbf{A}^{-1}\chi`. The
  multiplication operator :math:`\mathbf{K} = \mathbf{A}^{-1}\mathbf{F} =
  (\mathbf{A}^{-1}\chi)\,(\nu\Sigma_f)^{\mathsf T}` has **rank one**, so this
  is its only non-zero eigenvalue and, for a physical medium, its dominant
  one (proof: :ref:`homogeneous-rank-one-route`). The solver performs **one**
  LU solve, ``MatrixInverseOperator(problem.loss).apply(production.emission_spectrum)``,
  and one :math:`G`-term contraction with the production-rate co-vector. It
  calls **no eigen-solver**: the last bit of a dense eigen-solve belongs to
  the platform's LAPACK, not to the problem (`[M]` 2026-09-30, a macOS
  update moved ``geev``'s :math:`\kinf` by 1 ULP on an unchanged tree). There
  is **no power iteration** either: the 0-D problem is solved exactly, not by
  the iterative :func:`~orpheus.numerics.eigenvalue.power_iteration` the
  spatially-coupled solvers use (see :ref:`direct-eigensolve-solve`,
  :ref:`three-eigenvalue-engines`)
- **The rank one is the model's, not the solver's.** A mixture carries ONE
  fission spectrum :math:`\chi`, the production-weighted average of its
  isotopes' spectra taken at flat flux, so :math:`\mathbf{F} = \chi \otimes
  \nu\Sigma_f` has rank one by that approximation; the exact operator over
  :math:`K` fissile isotopes has rank up to :math:`K`
  (`#549 <https://github.com/deOliveira-R/ORPHEUS/issues/549>`_). A
  non-rank-one :math:`\mathbf{F}` is not spellable on this path: the fission
  kernel requires one-dimensional factors
- **Homogeneous is the FIRST production consumer of**
  :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
  (taxonomy step 5b). Constructing the matrix inverse *explicitly* — rather
  than calling the structure-keyed ``loss.inverse()``, which would return the
  **iterative** :class:`~orpheus.numerics.green_operator.GreenOperator`
  splitting — **is** the direct-realization strategy choice, encoded as a type
  rather than a flag. :func:`~orpheus.numerics.eigenvalue.direct_eigenvalue`
  (the ``(A, F)``-posed sibling engine) is **no longer on the homogeneous call
  path**
- **The problem is a HUB** (:ref:`the-problem-hub`).
  :class:`~orpheus.homogeneous.solver.HomogeneousProblem` — a frozen
  dataclass over one
  :class:`~orpheus.data.macro_xs.mixture.Mixture` — is the place the
  consumed objects live: the pose, the kernel-tier material fields, the
  cross-section fields, the bound operators and the reaction-rate
  co-vectors, each a per-instance ``cached_property`` minted from the
  mixture **and from nothing else**.
  :func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite` reads
  the hub and computes; it constructs no data of its own. (CS4c coda,
  ruling R-c1, 2026-09-08 — the hub is ``SNProblem``'s analogue on the
  infinite path; see :ref:`homogeneous-development-history`.)
- **Nothing is fabricated on the path.** There is no carrier, no
  ``[0, 1]`` edges, no invented node and no coordinate system: the
  infinite medium has no mesh, so the problem states none. Until the
  coda the cross sections came through a *meshless* single-cell
  :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` built by a
  ``from_materials`` factory whose edges, node and Cartesian chart were
  consumed by nothing; that factory retired with its last consumer.
- **A is assembled from the transport operators**, not a bespoke matrix:
  :math:`\mathbf{A} = C - K_\mathrm{iso}`, posed on
  :math:`V_E \otimes V_{\rm pt}` (the energy axis tensored with the
  quotient point) and reading its cross sections off fields the hub
  mints **on that same pose**. The space was the problem's own from
  campaign 1 CS4a (K2) — minted from the MIXTURE by
  :func:`~orpheus.homogeneous.solver._pose_space` and threaded into every
  arm — and the coda made the DATA follow it, so a field posed on a space
  nothing checks is no longer spellable here. With
  :math:`C = \text{diag}(\Sigma_t)` and
  :math:`K_\mathrm{iso} = \Sigma_{s0}^T + 2\Sigma_2^T` supplied by
  :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicScattering`
  and :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicN2N`.
  Streaming :math:`L` is identically zero in an infinite medium and is dropped,
  so the whole spectrum runs through the SAME operator algebra the meshed SN
  solver uses (cross-model single source, Cardinal Rule 2; campaign #276)
- This is the reference eigenvalue for ALL solvers on homogeneous problems
- Tolerance: :math:`\kinf` matches the registry's analytical values to
  :math:`< 10^{-12}`, and the **exact** answer for its float inputs
  (rational arithmetic, no rounding anywhere) within a forward-error bound
  derived per case, 3.75 to 43.2 ULP; `[M]` 2026-10-01 the measured
  deviation is at most 0.86 ULP on the eight shipped producing mixtures
  (:ref:`homogeneous-exact-reference`). There is no byte pin: the one LU
  solve on the path is the platform library's, and its last bits are not an
  ORPHEUS property
- **Gotcha**: this eigenvalue is flux-shape independent — it tests nothing
  about spatial or angular discretization


Overview
========

The infinite homogeneous medium is the simplest model in reactor physics.
Position drops out because the medium and its data are the same at every
point (translation symmetry), direction drops out because they prefer no
direction (rotation symmetry), and the neutron transport equation reduces
to a pure **energy balance** (:ref:`infinite-medium-point-in-phase-space`).
The only unknowns are the
**neutron energy spectrum** :math:`\phi(E)` and the **infinite
multiplication factor** :math:`\kinf`.

Despite its simplicity, the homogeneous model is the foundation on which
all other solvers build:

- It is the **first module** students encounter in the ORPHEUS
  curriculum, introducing the multi-group eigenvalue problem and its
  direct dense solution.
- The **cross-section preparation pipeline** — isotope loading,
  sigma-zero self-shielding, interpolation, macroscopic summation — is
  exercised here and reused unchanged by every subsequent solver (SN,
  MoC, CP, Monte Carlo, diffusion).
- Analytical eigenvalues for 1-, 2-, and 4-group homogeneous media
  serve as **verification benchmarks** for all deterministic solvers.

This chapter derives the infinite-medium eigenvalue problem from first
principles, describes the cross-section preparation pipeline, and
presents the direct rank-one solve used to compute :math:`\kinf`
and :math:`\phi(E)`.

The solver is the single function
:func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite`.  It reads
the problem's hub, :class:`~orpheus.homogeneous.solver.HomogeneousProblem`
— which assembles the loss operator from the model-shared transport
operators (see :ref:`direct-eigensolve`) — reads the dominant eigenpair of
:math:`\mathbf{A}^{-1}\mathbf{F}` off its rank-one structure (one solve and
one contraction, :ref:`homogeneous-rank-one-route`), and returns a
:class:`~orpheus.homogeneous.solver.HomogeneousResult`.


From the Boltzmann Equation to the Infinite Medium
====================================================

The Boltzmann Transport Equation
---------------------------------

The starting point is the steady-state neutron transport equation in its
integro-differential form :cite:`Duderstadt1976`:

.. math::
   :label: boltzmann

   \hat{\Omega} \cdot \nabla \psi(\mathbf{r}, \hat{\Omega}, E)
   + \Sigma_\mathrm{t}(\mathbf{r}, E) \, \psi(\mathbf{r}, \hat{\Omega}, E)
   = \int_0^\infty \!\!\int_{4\pi}
     \Sigma_\mathrm{s}(\mathbf{r}, E' \!\to\! E, \hat{\Omega}' \!\to\! \hat{\Omega})
     \, \psi(\mathbf{r}, \hat{\Omega}', E') \, d\Omega' \, dE'
   + \frac{\chi(E)}{4\pi \, k}
     \int_0^\infty \nu\Sigma_\mathrm{f}(\mathbf{r}, E')
     \, \phi(\mathbf{r}, E') \, dE'

.. vv-status: boltzmann documented

Here :math:`\psi(\mathbf{r}, \hat{\Omega}, E)` is the :term:`angular flux`,
:math:`\phi(\mathbf{r}, E) = \int_{4\pi} \psi \, d\Omega` is the scalar
flux, :math:`\chi(E)` is the fission spectrum, and :math:`k` is the
multiplication factor eigenvalue.


.. _infinite-medium-point-in-phase-space:

Simplification for the Infinite Homogeneous Medium
----------------------------------------------------

Eq. :eq:`boltzmann` is posed on the neutron's **phase space**, the
product of three factors: a position :math:`\mathbf{r} \in
\mathbb{R}^3` (a point of space), a direction of motion
:math:`\hat{\Omega} \in S^2` (a point of the unit sphere of directions)
and an energy :math:`E \in (0, \infty)`. The infinite homogeneous medium
removes the first two factors and keeps the third. Each removal has its
own condition, and neither condition is "the medium has no boundary":
both are **symmetries of the data**. The subsections below state which
symmetry removes which factor, and the result is the definition this
chapter rests on (:ref:`infinite-medium-definition`): in ORPHEUS the
infinite medium is the point in position **and** in direction, posed on
energy alone.


Two spaces, and the group that acts on them
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Position and direction are two different spaces. A position is a point
of :math:`\mathbb{R}^3`. A direction is a unit vector, a point of
:math:`S^2`, and it carries no position. A set of directions is
described on :math:`S^2` alone: the hemisphere
:math:`\{\hat{\Omega} : \hat{\Omega}\cdot\hat{n} > 0\}` is defined by a
reference **vector** :math:`\hat{n}`, not by a surface that splits
space, so describing it needs no boundary and no geometry.

The rigid motions of space, the Euclidean group
:math:`E(3) = \mathbb{R}^3 \rtimes O(3)`, act on phase space. An element
:math:`g = (\mathbf{a}, R)`, a translation by :math:`\mathbf{a}` and a
rotation or reflection :math:`R`, moves a neutron at
:math:`(\mathbf{r}, \hat{\Omega}, E)` to
:math:`(R\mathbf{r} + \mathbf{a},\, R\hat{\Omega},\, E)`, and acts on a
function of phase space by

.. math::

   (g\psi)(\mathbf{r}, \hat{\Omega}, E)
   = \psi\bigl(R^{\mathsf T}(\mathbf{r} - \mathbf{a}),\,
               R^{\mathsf T}\hat{\Omega},\, E\bigr).

A translation moves the position and leaves the direction alone; a
rotation moves both together. Write :math:`T` for the operator of
Eq. :eq:`boltzmann` (streaming, collision, scattering, and the fission
term), so that a source problem reads :math:`T\psi = q`. Each term
commutes with part of the group, under a condition on the medium:

- **Streaming** :math:`\hat{\Omega}\cdot\nabla` commutes with every
  element of :math:`E(3)`, unconditionally:
  :math:`\hat{\Omega}\cdot\nabla(g\psi) = g\,(\hat{\Omega}\cdot\nabla\psi)`,
  because rotating the direction and the gradient together leaves their
  dot product unchanged. It does **not** commute with a rotation of the
  direction alone; acting on direction it couples spherical-harmonic
  degree :math:`\ell` to :math:`\ell \pm 1`
  (:ref:`spherical-harmonics-eigenbasis`).
- **Collision, scattering and fission** commute with every translation
  if and only if the cross sections do not depend on position: the
  **homogeneous** medium. A homogeneous medium must fill all of space,
  because no translation maps a bounded body onto itself. This is the
  only role infiniteness plays.
- **Collision and scattering** commute with every rotation if and only
  if :math:`\Sigma_t` does not depend on direction and the scattering
  kernel depends on the two directions only through the cosine
  :math:`\mu_0 = \hat{\Omega}'\cdot\hat{\Omega}`: the **isotropic**
  medium. Bell and Glasstone name the exceptions, a moving medium and a
  single crystal (:cite:`BellGlasstone1970` §2.6a, p. 102). The fission
  term commutes with every rotation because its emission
  :math:`\chi/4\pi` is isotropic.

The consequence is one lemma. *If* :math:`T` *commutes with every
element of a group* :math:`G`, *and the source problem has one solution
(a subcritical medium), then the solution has exactly the symmetry of
the source:* :math:`\mathrm{Stab}(\psi) = \mathrm{Stab}(q)`, where
:math:`\mathrm{Stab}(f) = \{g \in G : gf = f\}`. The proof is three
lines. From :math:`Tg = gT` follows :math:`gT^{-1} = T^{-1}g`. If
:math:`gq = q`, then :math:`g\psi = gT^{-1}q = T^{-1}gq = T^{-1}q = \psi`.
Conversely, if :math:`g\psi = \psi`, then
:math:`gq = gT\psi = Tg\psi = T\psi = q`. A symmetric source therefore
forces a symmetric solution, an asymmetric source forces an asymmetric
one, and the solution never gains or loses a symmetry that the data do
not have.


The first collapse: position, by the translations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In a homogeneous medium with a source that does not depend on position,
every translation leaves the data unchanged. By the lemma the solution
does not depend on position either: :math:`\nabla\psi = 0`, and the
streaming term vanishes identically. This is exact; it is not an
approximation that neglects a small leakage.

Infiniteness alone does not give it. An isotropic point source in an
infinite homogeneous medium produces a flux that falls with the distance
from the source (:cite:`BellGlasstone1970` §2.2f), because the
translations map the medium onto itself but move the source.
Uniformity is the translation symmetry of the **data**.

Every position is then equivalent to every other, and the position
factor collapses to a point: the quotient of :math:`\mathbb{R}^3` by its
translations. What remains is posed on direction and energy:

.. math::
   :label: inf-med-direction-energy

   \Sigma_t(E)\,\psi(\hat{\Omega}, E)
   = \int_0^\infty \!\! \int_{4\pi}
       \Sigma_s(E' \!\to\! E,\, \hat{\Omega}'\cdot\hat{\Omega})\,
       \psi(\hat{\Omega}', E')\, d\Omega'\, dE'
     + \frac{\chi(E)}{4\pi\,k} \int_0^\infty
       \nu\Sigma_f(E')\,\phi(E')\,dE'
     + q(\hat{\Omega}, E) .

.. (vv-status rationale) Literature-transcribed reduction: Eq. boltzmann with
   the streaming term removed by translation invariance (the position quotient),
   with an external source added. A statement of what the infinite medium's
   transport equation is before the direction collapse, not a solver claim; no
   ORPHEUS module poses this direction-and-energy equation.
.. vv-status: inf-med-direction-energy documented

The eigenvalue question takes :math:`q = 0`; a source question in a
non-multiplying medium drops the fission term. The direction sphere
**survives** this collapse.

For the eigenvalue question the data are the cross sections alone, and
the fission emission is built from them, so the translations act on it
as they act on the medium. The question asks for the
translation-invariant mode. In the plane-geometry Fourier family
:math:`e^{-iBx}\,\psi(B, \mu, E)` of the :math:`B_N` method that mode is
the member :math:`B = 0` (:cite:`BellGlasstone1970` §4.5c,
Eqs. (4.65)–(4.66)). The spaces chapter records the same family from the
side of the spaces: the spatial axis is quotiented to a one-point axis
whose unit weight is the per-unit-volume convention, and the buckling
:math:`B` parameterises the intermediate members
(:ref:`spaces-quotient-family`; clause 1 of
:ref:`spaces-collapse-doctrine-standing`).


The second collapse: direction, by the rotations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On a function that no longer depends on position, a rotation acts on
the direction alone, and the scattering operator of an isotropic medium
commutes with it. Its kernel is **zonal**, a function of
:math:`\hat{\Omega}'\cdot\hat{\Omega}` only, so it expands in Legendre
polynomials (:cite:`BellGlasstone1970` §2.6a, Eqs. (2.77)–(2.78)):

.. math::

   \Sigma_s(E' \!\to\! E, \mu_0)
   = \sum_{\ell=0}^{\infty} \frac{2\ell + 1}{4\pi}\,
     \Sigma_{s,\ell}(E' \!\to\! E)\, P_\ell(\mu_0),
   \qquad
   \Sigma_{s,\ell}(E' \!\to\! E) = 2\pi \int_{-1}^{1}
     \Sigma_s(E' \!\to\! E, \mu_0)\, P_\ell(\mu_0)\, d\mu_0 .

By the Funk–Hecke theorem every spherical harmonic of degree
:math:`\ell` is an eigenfunction of the scattering operator, with
eigenvalue :math:`\Sigma_{s,\ell}` (:eq:`sh-funk-hecke-eigenvalue`; the
eigenvalue does not depend on the order :math:`m` by Schur's lemma,
:ref:`spherical-harmonics-eigenbasis`). Expand the angular flux and the
source in the real harmonics of the project's convention
(:ref:`spherical-harmonics`),

.. math::

   \psi(\hat{\Omega}, E) = \sum_{\ell=0}^{\infty} \frac{2\ell + 1}{4\pi}
     \sum_{m=-\ell}^{\ell} \psi_\ell^m(E)\, Y_\ell^m(\hat{\Omega}),
   \qquad
   \psi_\ell^m(E) = \int_{4\pi} \psi(\hat{\Omega}, E)\,
     Y_\ell^m(\hat{\Omega})\, d\Omega ,

so that :math:`\psi_0^0 = \phi` and :math:`(\psi_1^{-1}, \psi_1^0,
\psi_1^1)` are the components of the current :math:`\mathbf{J}`, and
project Eq. :eq:`inf-med-direction-energy` onto each
:math:`Y_\ell^m`. The equation splits into one equation per degree and
order, with no coupling between them:

.. math::
   :label: inf-med-moment-decoupling

   \Sigma_t(E)\,\psi_\ell^m(E)
   - \int_0^\infty \Sigma_{s,\ell}(E' \!\to\! E)\,\psi_\ell^m(E')\,dE'
   = q_\ell^m(E)
     + \delta_{\ell 0}\,\frac{\chi(E)}{k} \int_0^\infty
       \nu\Sigma_f(E')\,\phi(E')\,dE' .

.. (vv-status rationale) Literature-transcribed decoupling: the projection of
   inf-med-direction-energy onto the real spherical harmonics, using the
   Legendre expansion of a zonal kernel (Bell & Glasstone 1970 Eqs. 2.77-2.78)
   and the Funk-Hecke eigenvalue (sh-funk-hecke-eigenvalue). A derivation step,
   not a solver claim; the homogeneous solver poses only its l = 0 member.
.. vv-status: inf-med-moment-decoupling documented

The fission emission appears only at :math:`\ell = 0` because it is
isotropic. (An :math:`(n,2n)` emission is a second zonal kernel and
enters every degree in the same way, :math:`\Sigma_{s,\ell} \to
\Sigma_{s,\ell} + 2\,\Sigma_{2,\ell}`.)

**When the source prefers no direction**, :math:`q_\ell^m = 0` for every
:math:`\ell \ge 1`. The eigenvalue question qualifies, since its only
emission is fission. Then every equation with :math:`\ell \ge 1` is
homogeneous, and its operator :math:`\Sigma_t - \Sigma_{s,\ell}` is
invertible whenever the :math:`\ell = 0` removal operator
:math:`\Sigma_t - \Sigma_{s,0}` is, and that one is the loss operator
whose inverse this chapter's solve applies
(:ref:`homogeneous-rank-one-route`). The reason is a bound: the differential cross section is non-negative and
:math:`|P_\ell| \le 1`, so :math:`|\Sigma_{s,\ell}| \le \Sigma_{s,0}`
entry by entry. In one group this is plain arithmetic,
:math:`\Sigma_t - \Sigma_{s,\ell} \ge \Sigma_t - \Sigma_s = \Sigma_a > 0`.
In the multigroup form, with :math:`D = \mathrm{diag}(\Sigma_t)`, the
comparison theorem for non-negative matrices gives
:math:`\rho\bigl(D^{-1}|\boldsymbol{\Sigma}_{s,\ell}^{\mathsf T}|\bigr)
\le \rho\bigl(D^{-1}\boldsymbol{\Sigma}_{s,0}^{\mathsf T}\bigr) < 1`, and
the right-hand inequality is the statement that the :math:`\ell = 0`
removal matrix is a non-singular M-matrix. The left-hand side bounds
:math:`\rho\bigl(D^{-1}\boldsymbol{\Sigma}_{s,\ell}^{\mathsf T}\bigr)`,
so :math:`D - \boldsymbol{\Sigma}_{s,\ell}^{\mathsf T}
= D\,\bigl(I - D^{-1}\boldsymbol{\Sigma}_{s,\ell}^{\mathsf T}\bigr)` is
invertible by its Neumann series. So :math:`\psi_\ell^m = 0` for
every :math:`\ell \ge 1`, the angular flux is isotropic,
:math:`\psi = \phi/4\pi`, and the :math:`\ell = 0` equation is the
**energy balance** on which the rest of this chapter is built:

.. math::
   :label: inf-hom-balance

   \Sigt{} \phi(E)
   = \int_0^\infty \Sigma_{\mathrm{s},0}(E' \!\to\! E) \, \phi(E') \, dE'
     + \frac{\chi(E)}{k}
       \int_0^\infty \nu\Sigma_\mathrm{f}(E') \, \phi(E') \, dE'


.. implements:: inf-hom-balance
   :by: orpheus.homogeneous.solver.HomogeneousProblem.loss

   **Implemented by** 7 sites. Every symbol that executes this
   equation's arithmetic is declared, not only the canonical one: a
   test is adjudicated against the transcription it actually ran, so
   declaring a single site would refute the tests that exercise the
   others.

.. implements:: inf-hom-balance
   :by: orpheus.homogeneous.solver.solve_homogeneous_infinite

.. implements:: inf-hom-balance
   :by: orpheus.derivations.common.eigenvalue._infinite_medium_matrices

.. implements:: inf-hom-balance
   :by: orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous

.. implements:: inf-hom-balance
   :by: orpheus.derivations.common.eigenvalue.kinf_homogeneous

.. implements:: inf-hom-balance
   :by: orpheus.derivations.continuous.analytical.homogeneous.derive_1g

.. implements:: inf-hom-balance
   :by: orpheus.derivations.continuous.analytical.homogeneous.derive_1g_continuous

where :math:`\Sigma_{\mathrm{s},0}` is the isotropic scattering kernel.

.. note::

   In the infinite homogeneous medium, scattering merely redistributes
   neutrons in energy.  It does not change the total production-to-loss
   ratio, so :math:`\kinf` depends only on the fission and absorption
   cross sections.  Scattering does, however, determine the **shape** of
   the neutron spectrum :math:`\phi(E)` — specifically the 1/E
   slowing-down region and the thermal Maxwellian peak.

The anisotropic moments :math:`\Sigma_{s,\ell \ge 1}` are present in a
material's data and inert here: they act only on modes the data do not
excite. The isotropy of the flux is therefore the **data's**, not the
medium's alone. The same isotropic medium driven by a source that prefers
a direction keeps that direction
(:ref:`infinite-medium-anisotropic-counterexample`).


.. _infinite-medium-definition:

ORPHEUS's infinite medium: the point in position and in direction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"A point" has two meanings in this reduction. The quotient by the
translations collapses **position** to a point, and the direction sphere
at that point survives (Eq. :eq:`inf-med-direction-energy`). The
direction sphere collapses as well only when the whole problem, data
included, is invariant under rotations (Eq. :eq:`inf-med-moment-decoupling`
with a source that prefers no direction).

**Definition.** In ORPHEUS the :term:`infinite medium` is the point in
both. It is posed on **energy alone**, and Eq. :eq:`inf-hom-balance` is
its whole transport equation. It holds a material and nothing else. It
has **no geometry**, because a geometry admits spatial dimension, and
with spatial dimension come more than one material (a map from regions to
materials) and a direction chart (a frame in which directions are
written as coordinates). The infinite medium has none of the three.

The definition is exact for everything it admits, for two reasons:

- **Every material is isotropic.** A
  :class:`~orpheus.data.macro_xs.mixture.Mixture` stores its scattering as
  the Legendre moments :math:`\Sigma_{s,\ell}` of a kernel in
  :math:`\mu_0` alone (``SigS[l]``, one transfer matrix per order;
  :ref:`scattering-matrix-convention`), so the rotation condition on the
  medium holds by construction. A medium with a preferred axis, such as a
  single crystal, cannot be represented.
- **Every question it admits is invariant under rotations**: an
  eigenvalue of the medium, such as :math:`\kinf`, or a source that
  depends on no coordinate. With no direction chart there is no way to
  state a source that prefers a direction.

The shipped solver poses exactly this object.
:class:`~orpheus.homogeneous.solver.HomogeneousProblem` holds one
:class:`~orpheus.data.macro_xs.mixture.Mixture` and poses its operators
on the energy axis tensored with the quotient point,
:math:`V_E \otimes V_{\rm pt}` (:ref:`the-problem-hub`). It has no
angular axis, and its spatial slot is the one-point axis of the
translation quotient.

.. dropdown:: First got wrong: the infinite medium as the point in position only
   :color: muted

   **What was tried.** Two readings preceded the definition above.

   #. This section's earlier text listed three conditions. The first read,
      verbatim, "Infinite geometry — no boundaries, so the flux is
      spatially uniform". The third derived an isotropic angular flux from
      "isotropic sources" without saying where the isotropy comes from.
   #. In the design of the reference specification, part of the
      reference-solution architecture that is step 2 of
      `#405 <https://github.com/deOliveira-R/ORPHEUS/issues/405>`_, the
      infinite medium was spelled as a geometry value with no extent (a
      review finding of 2026-09-25).
      The discussion that refuted it took eight exchanges, and it turned
      on a second reading: the infinite medium as the point in
      **position** only, with the direction sphere surviving at that
      point, so that it would admit data with a preferred direction (a
      hemisphere of directions, a source that prefers an axis, a beam).

   **Why each failed.**

   - *No boundary does not mean uniform.* Uniformity is the translation
     symmetry of the data, not a consequence of infiniteness; an
     isotropic point source in an infinite medium is the counterexample
     (the first collapse above).
   - *An isotropic flux needs an isotropic medium and isotropic data.* The
     direction sphere collapses only when the scattering kernel and the
     source are both invariant under rotations (the second collapse
     above).
   - *A hemisphere needs no geometry.* It is defined by a reference vector
     on :math:`S^2`, not by a surface in space.
   - *A geometry with no placement and no finiteness is not a geometry.*
     It defines nothing a geometry defines, so it carries no mathematical
     information, and it misdirects: it admits spatial dimension, and
     with it more than one material and a direction chart.
   - *A beam is not a direction-dependent source.* A beam is localised in
     space, so it needs spatial dimension (:ref:`infinite-medium-beam`).
   - *The direction-dependent problem that does exist is another object.*
     It is the spatial marginal of a body problem, the retraction of a
     body solution along the spatial axis
     (:ref:`infinite-medium-spatial-marginal`), so it needs no
     infinite-medium type of its own.

   **What replaced it.** The definition above, the user's ruling of
   2026-10-02 (the History table at the end of this page).


.. _infinite-medium-anisotropic-counterexample:

What the definition excludes: a uniform source that prefers a direction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Take one speed, an unbounded homogeneous medium, and a volumetric source
that is the same at every position but emits preferentially along a unit
vector :math:`\hat{n}`:

.. math::

   q(\hat{\Omega}) = \frac{Q}{4\pi}\,
     \bigl(1 + 3a\,\hat{\Omega}\cdot\hat{n}\bigr),
   \qquad |a| \le \tfrac{1}{3} .

Its total emission is :math:`Q = \int_{4\pi} q\,d\Omega`, its first moment
is :math:`\mathbf{q}_1 = \int_{4\pi} \hat{\Omega}\,q\,d\Omega = Q a\,\hat{n}`,
and the bound on :math:`a` keeps :math:`q \ge 0`. Every translation leaves
this source unchanged, so the first collapse applies and the problem is
Eq. :eq:`inf-med-direction-energy`. The rotations do not: only the
rotations about :math:`\hat{n}` fix the source, so by the lemma the
solution has exactly that axial symmetry, and the direction sphere
survives.

*A pure absorber.* With no scattering, Eq. :eq:`inf-med-direction-energy`
is algebraic:

.. math::

   \psi(\hat{\Omega}) = \frac{q(\hat{\Omega})}{\Sigma_t},
   \qquad
   \mathbf{J} = \int_{4\pi} \hat{\Omega}\,\psi\,d\Omega
              = \frac{Q a}{\Sigma_t}\,\hat{n} .

The same result is the integral of the source along the backward
characteristic, :math:`\int_0^\infty e^{-\Sigma_t s}\, q(\hat{\Omega})\,ds`.
A uniform current flows everywhere, in the direction the source prefers,
with no gradient of anything to drive it.

*Isotropic scattering.* With :math:`\Sigma_s(\mu_0) = \Sigma_s/4\pi`,

.. math::
   :label: inf-med-anisotropic-source

   \psi(\hat{\Omega})
   = \frac{1}{\Sigma_t}\left[\frac{\Sigma_s\,\phi}{4\pi}
     + q(\hat{\Omega})\right],
   \qquad
   \phi = \frac{Q}{\Sigma_a} .

.. (vv-status rationale) Closed-form worked example (one speed, isotropic
   scattering, a uniform source with a dipole term): a teaching counterexample
   to the definition of the infinite medium, solved from
   inf-med-direction-energy. No ORPHEUS solver poses it, since the definition
   excludes it; not a solver claim.
.. vv-status: inf-med-anisotropic-source documented

Integrating over directions gives :math:`\Sigma_t\phi = \Sigma_s\phi + Q`,
hence :math:`\phi = Q/\Sigma_a`. The **collided** part
:math:`\Sigma_s\phi/(4\pi\Sigma_t)` is isotropic: one isotropic scatter
erases the direction. The **uncollided** part
:math:`q(\hat{\Omega})/\Sigma_t` keeps it, so the current is still
:math:`\mathbf{J} = Q a\,\hat{n}/\Sigma_t`.

*Any isotropic medium.* Eq. :eq:`inf-med-moment-decoupling` solves the
general case degree by degree,
:math:`\psi_\ell^m = q_\ell^m/(\Sigma_t - \Sigma_{s,\ell})`. At
:math:`\ell = 1` it gives, exactly and with no diffusion approximation,

.. math::

   \mathbf{J} = \frac{\mathbf{q}_1}{\Sigma_t - \Sigma_{s,1}}
              = \frac{\mathbf{q}_1}{\Sigma_{\rm tr}} .

The current that a direction-preferring source sustains is its first
moment times the :term:`transport mean free path`
:math:`1/\Sigma_{\rm tr}`, the relaxation length of the :math:`\ell = 1`
mode (:ref:`diffusion-transport-xs-relaxation`). Isotropic scattering is
the case :math:`\Sigma_{s,1} = 0`, where the length is :math:`1/\Sigma_t`.

The problem is well posed and its answer is exact, but it is not an
infinite medium in ORPHEUS's sense. It is posed on direction and energy,
not on energy alone. Its data name a vector :math:`\hat{n}`, so stating
them needs a direction chart. Its answer carries a vector field,
:math:`\mathbf{J}`, which an energy-only problem cannot hold.


.. _infinite-medium-spatial-marginal:

The resolution: a spatial marginal, not a specification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The direction-dependent problem has a home, and it is not the infinite
medium. Eq. :eq:`boltzmann` is **Eulerian**: it balances the neutrons at
a fixed point of space. The direction-dependent problem is
**Lagrangian**: follow each neutron and forget where it is. Its unknown
is the **spatial marginal** of a body problem,

.. math::

   \bar{\psi}(\hat{\Omega}, E, t)
   = \int_{\mathbb{R}^3} \psi(\mathbf{r}, \hat{\Omega}, E, t)\, dV ,

which a Monte Carlo code computes as a tally over all positions.

To find its equation, integrate the time-dependent transport equation
over all space. The streaming term is a divergence,
:math:`\hat{\Omega}\cdot\nabla\psi = \nabla\cdot(\hat{\Omega}\,\psi)`
(:cite:`BellGlasstone1970` §1.1e, Eq. (1.16)), so by the divergence
theorem its integral over a region is the angular outflow through the
region's surface, :math:`\oint (\hat{\Omega}\cdot\hat{\mathbf{n}}_A)\,
\psi\, dA`. This is the step that turns streaming into leakage in the
conservation relation, Bell and Glasstone's Eq. (1.19). For a pulse of
neutrons released at :math:`t = 0` from a bounded region into an
unbounded medium, every neutron is within :math:`vt` of that region, so a
surface enclosing the whole population carries no neutrons and the
leakage is zero. (For a steady source in a subcritical unbounded medium
the flux decays exponentially with distance and the leakage through a
receding surface vanishes in the limit.) What remains is closed in
direction and energy. In one speed, with the path length :math:`s = vt`
travelled since the release, it reads

.. math::
   :label: inf-med-spatial-marginal

   \frac{\partial \bar{\psi}}{\partial s}(\hat{\Omega}, s)
   + \Sigma_t\,\bar{\psi}(\hat{\Omega}, s)
   = \int_{4\pi} \Sigma_s(\hat{\Omega}'\cdot\hat{\Omega})\,
     \bar{\psi}(\hat{\Omega}', s)\, d\Omega'
   + \bar{q}(\hat{\Omega}, s) .

.. (vv-status rationale) Literature-transcribed: the spatial integral of the
   time-dependent one-speed transport equation (Bell & Glasstone 1970
   Eq. (1.16) and the divergence-theorem step of Eq. (1.19)); the
   space-independent equation of their Ch. 1, Exercise 16. No ORPHEUS module
   evolves it; not a solver claim.
.. vv-status: inf-med-spatial-marginal documented

This is the **spatially homogeneous** transport equation. Bell and
Glasstone pose it, for a source-free infinite medium with direction and
energy retained, as the "space-independent neutron transport equation"
(:cite:`BellGlasstone1970` Ch. 1, Exercise 16, p. 61). It is the
zero-wavenumber member of the Fourier family of the first collapse,
because :math:`\int \psi\, dV` is the spatial Fourier transform at
:math:`B = 0`. Its steady form is Eq. :eq:`inf-med-direction-energy`,
and the uniform direction-preferring problem above is the same equation
read per unit volume.

In ORPHEUS's vocabulary the marginal is the **retraction of a body
solution along the spatial axis**, the arrow
:class:`~orpheus.numerics.operator.AxisRetractionOperator` that
integrates a field over one axis with that axis's measure
(:ref:`spaces-collapse-pair`). Two different collapses of the spatial
axis therefore meet at the same equation, and they are different
operations:

- the **quotient** by the translations, which acts on a
  translation-invariant field and normalises it, per unit volume, since a
  non-zero constant on :math:`\mathbb{R}^3` cannot be integrated (clause
  1 of :ref:`spaces-collapse-doctrine-standing`); this is the infinite
  medium's position collapse;
- the **retraction**, which acts on a body solution and integrates it,
  which is possible because the body, or the pulse, is bounded.

The direction-dependent problem is the second. It needs no specification
type of its own: it is posed as a body problem, solved, and retracted
along the spatial axis.


.. _infinite-medium-beam:

A beam is not a direction-dependent source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A beam is collimated: its source is concentrated near a line in space,
for example
:math:`q(\mathbf{r}, \hat{\Omega}) = \delta^2(\mathbf{r}_\perp)\,
\delta(\hat{\Omega} - \hat{\Omega}_0)` per unit length along the line,
with :math:`\mathbf{r}_\perp` the position transverse to it. It is
localised in space, so it needs spatial dimension and cannot be posed on
a point. The uniform direction-preferring source of the counterexample
has no spatial structure at all: its anisotropy lives entirely on
:math:`S^2`.

How far from its source a beam's direction is remembered is set by the
:term:`transport mean free path`
:math:`\lambda_{\rm tr} = 1/\Sigma_{\rm tr}`, the relaxation length of
the :math:`\ell = 1` mode of Eq. :eq:`inf-med-spatial-marginal`
(:ref:`diffusion-transport-xs-relaxation`). A detector a few transport
mean free paths from the beam axis sees no trace of the beam's
direction. The angular distribution there is the asymptotic one of
diffusion theory, whose remaining anisotropy is the current that the
local flux gradient drives (Fick's law), not the beam's;
Duderstadt and Hamilton state the validity condition of diffusion theory
as being "over several mean free paths from any sources or boundaries in
a weakly absorbing medium" (:cite:`Duderstadt1976` §4-IV, p. 138). With
isotropic scattering only the uncollided neutrons remember the
direction, and the length is :math:`1/\Sigma_t`.


.. _infinite-medium-reflective-images:

The reflective slab is the infinite medium for mirror-symmetric data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A slab :math:`0 \le x \le a` with both faces reflective represents an
infinite medium by images. Write :math:`\mu = \hat{\Omega}\cdot\hat{e}_x`.
Reflection across the face :math:`x = 0` maps a neutron at
:math:`(x, \mu)` to :math:`(-x, -\mu)`, and a reflective face is exactly
the statement that the solution is unchanged by that map. Unfold the
slab across its faces again and again: each image of the slab is a copy
of the same medium, and the image of a uniform source :math:`q(\mu)` is
:math:`q(-\mu)`. The unfolded problem is the infinite medium whose source
alternates between :math:`q(\mu)` and :math:`q(-\mu)` from one slab
width to the next.

- **If the data are mirror-symmetric**, :math:`q(\mu) = q(-\mu)`, the
  unfolded source is uniform, the unfolded problem is the infinite
  homogeneous medium, and the reflective slab reproduces it exactly: the
  flux is flat in :math:`x` and equals the infinite-medium flux. A flux
  that is not flat in a reflective slab with uniform, mirror-symmetric
  data is therefore numerical error (an unconverged iteration or a
  defect), never physics.
- **If they are not**, the unfolded source alternates with period
  :math:`2a` and the solution is not uniform. The reflective slab is then
  a different problem from the uniform direction-preferring medium of
  the counterexample: a mirror forces zero net current through its plane,
  while that medium carries the uniform current
  :math:`\mathbf{J} = \mathbf{q}_1/\Sigma_{\rm tr}`.

The same argument holds in a box with every face reflective, with the
group generated by the face mirrors in place of the single reflection.
Data invariant under every rotation and reflection are invariant under
that group, so every question ORPHEUS's infinite medium admits is
reproduced exactly by an all-reflective slab or box. Two gates assert it
for the eigenvalue question:
``tests/gates/sn/solve/test_d3_admission.py::test_kinf_3d_equals_2d_equals_1d_homogeneous_reflective``
(S\ :sub:`N` at :math:`d = 1, 2, 3`, two and four groups,
:math:`|k - \kinf| < 10^{-8}` against the registry's infinite-medium
eigenvalue) and
``tests/gates/diffusion/test_solver.py::TestInfiniteMedium::test_reflective_slab_reproduces_k_infinity``
(diffusion on a non-uniform mesh, two and three groups,
:math:`|k - \kinf| < 10^{-11}` against
:func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite`, and a flux
flat in space to :math:`10^{-10}` of its maximum).


Multi-Group Energy Discretisation
==================================

Group-Averaged Cross Sections
------------------------------

The continuous energy variable is discretised into :math:`G` groups.
The fastest group carries the highest energies, the last group the
lowest (thermal neutrons): the descending boundary array
:math:`E_0 > E_1 > \cdots > E_G` defines the grid stored in
:attr:`Mixture.eg` for production cases (XS computed from ENDF
:class:`~orpheus.data.micro_xs.isotope.Isotope` data via
:func:`compute_macro_xs`). This is the :ref:`canonical fast-first
energy-group convention <canonical-group-convention>`; the code index
runs :math:`g = 0` (fastest) to :math:`g = G-1` (thermal), while the
1-based labels used in the equations below (group 1 = fastest) are a
presentation choice for the slowing-down algebra. For synthetic verification cases (Sood-style
abstract XS, MMS test mixtures), :attr:`Mixture.eg` is ``None`` —
there is no real grid, only a discrete set of group cross-sections.
Per-energy diagnostics (:term:`lethargy` widths, flux-per-energy plots,
spectrum-weighted condensation) require the grid to be populated and
gracefully skip the synthetic-XS path.

The group flux is the integral over the group's energy interval:

.. math::
   :label: group-flux

   \phi_g = \int_{E_g}^{E_{g-1}} \phi(E) \, dE

.. vv-status: group-flux documented

Group-averaged cross sections are **flux-weighted** averages:

.. math::
   :label: group-xs

   \Sigt{g} = \frac{1}{\phi_g} \int_{E_g}^{E_{g-1}} \Sigt{}(E) \, \phi(E) \, dE

.. vv-status: group-xs documented

In practice, these averages are pre-computed and stored in the
421-group HELIOS library that ships with ORPHEUS.  The library provides
cross sections tabulated at several background cross section
(:math:`\sigma_0`) values; the sigma-zero iteration (see
:ref:`sigma-zero-iteration`) selects the appropriate value for each
isotope and group.


.. _mg-eigenvalue-problem:

The Multi-Group Neutron Balance
--------------------------------

Substituting group-averaged quantities into Eq. :eq:`inf-hom-balance`
gives the **multi-group neutron balance** for group :math:`g`:

.. math::
   :label: mg-balance

   \Sigt{g} \, \phi_g
   = \sum_{g'=1}^{G} \Sigs{g' \to g} \, \phi_{g'}
     + \frac{\chi_g}{k} \sum_{g'=1}^{G} \nSigf{g'} \, \phi_{g'}


.. implements:: mg-balance
   :by: orpheus.cp.solver.CPSolver._compute_balance_residual

   **Implemented by** 12 sites. Every symbol that executes this
   equation's arithmetic is declared, not only the canonical one: a
   test is adjudicated against the transcription it actually ran, so
   declaring a single site would refute the tests that exercise the
   others.

   Four of the twelve were re-pointed on 2026-09-04 (#426 step 2), when
   the two collision-gain channels collapsed onto one family: the
   per-material P\ :sub:`0` verb and the scalar energy binding's
   ``apply`` moved to the shared
   :class:`~orpheus.transport.material_field.TransferMaterialField` /
   :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicTransfer`
   cores, ``LegendreMomentScattering`` became
   :class:`~orpheus.transport.operators.transfer.LegendreMomentTransfer`
   — the ONE :math:`\Lambda` both channels use — and the retired
   ``N2NMomentOperator`` row is now the :math:`(n,2n)` term itself,
   :class:`~orpheus.transport.operators.n2n.N2NOperator`, since its
   moment factor is no longer a class of its own.

.. implements:: mg-balance
   :by: orpheus.homogeneous.solver.HomogeneousProblem.loss

.. implements:: mg-balance
   :by: orpheus.homogeneous.solver.solve_homogeneous_infinite

.. implements:: mg-balance
   :by: orpheus.moc.core.MOCSolver.solve_fixed_source

.. implements:: mg-balance
   :by: orpheus.transport.material_field.TransferMaterialField.add_p0_source

.. implements:: mg-balance
   :by: orpheus.transport.operators.isotropic_transfer.IsotropicFission.apply

.. implements:: mg-balance
   :by: orpheus.transport.operators.isotropic_transfer.IsotropicTransfer.apply

.. implements:: mg-balance
   :by: orpheus.transport.operators.transfer.LegendreMomentTransfer

.. implements:: mg-balance
   :by: orpheus.transport.operators.n2n.N2NOperator

.. implements:: mg-balance
   :by: orpheus.transport.operators.scattering.ScatteringOperator

.. implements:: mg-balance
   :by: orpheus.derivations.common.eigenvalue._infinite_medium_matrices

.. implements:: mg-balance
   :by: orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous

The first term on the right is in-scattering from all groups
(including self-scattering :math:`g' = g`), and the second is the
fission source weighted by the fission spectrum :math:`\chi_g`.


Matrix Form
------------

Collecting all :math:`G` group equations into vectors and matrices:

.. math::
   :label: matrix-eigenvalue

   \mathbf{A} \, \boldsymbol{\phi}
   = \frac{1}{k} \, \mathbf{F} \, \boldsymbol{\phi}

where the **removal matrix** and **fission matrix** are:

.. math::
   :label: removal-matrix

   \mathbf{A} = \mathrm{diag}(\Sigt{g})
                - \boldsymbol{\Sigma}_{\mathrm{s}}^T
                - 2 \, \boldsymbol{\Sigma}_2^T

.. math::
   :label: fission-matrix

   \mathbf{F} = \boldsymbol{\chi} \otimes \nu\boldsymbol{\Sigma}_\mathrm{f}

Here :math:`\boldsymbol{\Sigma}_{\mathrm{s}}` is the :math:`G \times G`
scattering transfer matrix (:math:`P_0` component) and
:math:`\boldsymbol{\Sigma}_2` is the :math:`(n,2n)` transfer matrix. The
production matrix :math:`\mathbf{F}` is the **rank-1 dyad**
:math:`\boldsymbol{\chi} \otimes \nu\boldsymbol{\Sigma}_\mathrm{f}`
embodied by
:class:`~orpheus.transport.operators.isotropic_transfer.IsotropicFission`
— a group contraction onto the production rate followed by a broadcast
across the emission spectrum :math:`\boldsymbol{\chi}`.

.. note::

   **The** :math:`(n,2n)` **reaction appears ONLY in the loss matrix**
   :math:`\mathbf{A}` (as :math:`-2\boldsymbol{\Sigma}_2^T`), never in
   the production :math:`\mathbf{F}`.  The :math:`(n,2n)` event is a
   **loss-side multiplicity-2 transfer**: one neutron of group
   :math:`g'` is removed and **two** neutrons are deposited into the
   scattering system with the :math:`(n,2n)` energy-transfer kernel
   :math:`2\boldsymbol{\Sigma}_2(g' \!\to\! g)`.  The factor of two is
   the emission multiplicity; the transpose puts the *source* group on
   the row exactly as for :math:`\boldsymbol{\Sigma}_{\mathrm{s}}^T`
   (see :ref:`scattering-matrix-convention`).

   The two emitted neutrons are **not** produced with the fission
   spectrum :math:`\boldsymbol{\chi}` — they carry the :math:`(n,2n)`
   transfer kernel, not :math:`\boldsymbol{\chi}`.  Production is
   :math:`\nu\boldsymbol{\Sigma}_\mathrm{f}` only.  This matches the
   analytical oracle
   :func:`~orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous`
   (:math:`\mathbf{A} = \text{diag}(\Sigma_t) - (\Sigma_s + 2\Sigma_2)^T`,
   :math:`\mathbf{F} = \chi \otimes \nu\Sigma_f`) and the collision-probability
   oracle :func:`~orpheus.derivations.common.eigenvalue.kinf_from_cp`.

   .. warning::

      A retired bespoke formulation put :math:`(n,2n)` in **both**
      matrices — :math:`+2\,\text{colsum}(\boldsymbol{\Sigma}_2)` in the
      production numerator as well as :math:`-2\boldsymbol{\Sigma}_2^T`
      in the loss.  That double-counts the :math:`(n,2n)` neutrons.  On
      the asymmetric-:math:`\boldsymbol{\Sigma}_2` ``homo_2eg_n2n`` case
      it moves :math:`\kinf` from the correct ``1.6532`` to ``2.08`` — a
      :math:`\sim 0.43` error, far above the FP floor.  See the
      :ref:`direct-eigensolve` section for the live assembly and the
      ``homo_2eg_n2n`` de-vacuum case.

   The loss matrix :math:`\mathbf{A} = C - K_\mathrm{iso}` is assembled
   from the transport operators
   :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicScattering`
   (:math:`\Sigma_{s0}^T`) and
   :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicN2N`
   (:math:`2\Sigma_2^T`); the production dyad :math:`\mathbf{F}` is the
   :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicFission`
   rank-1 form (materialised densely as
   :math:`\chi \otimes \nu\Sigma_f`).  See
   :func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite`.

The eigenvalue :math:`k = \kinf` is the largest eigenvalue of the
generalised problem :eq:`matrix-eigenvalue`.  By the Perron–Frobenius
theorem :cite:`Hebert2009`, the dominant eigenvector :math:`\boldsymbol{\phi}`
is the unique non-negative solution — the **fundamental mode** — which
is the physically meaningful neutron spectrum.


.. _scattering-matrix-convention:

Scattering Matrix Convention
-----------------------------

The scattering transfer matrix :math:`\boldsymbol{\Sigma}_{\mathrm{s}}`
is stored in the **from-row, to-column** convention:

.. math::
   :label: sigs-convention

   (\boldsymbol{\Sigma}_{\mathrm{s}})_{g',g}
   = \Sigs{g' \to g}

.. vv-status: sigs-convention documented

That is, row :math:`g'` gives the **source group** and column :math:`g`
gives the **destination group**.  A downscatter-only matrix is therefore
**upper-triangular**: non-zero entries only on or above the diagonal,
because :math:`\Sigs{g' \to g} = 0` when :math:`g < g'` (no neutrons
scatter from thermal to fast).  Its transpose — the form that acts on
:math:`\boldsymbol{\phi}` in the in-scattering sum below — is
correspondingly lower-triangular, which is why the two-group operator
:eq:`two-group-A` carries its off-diagonal entry
:math:`-\Sigs{1 \to 2}` below the diagonal.

The neutron balance :eq:`mg-balance` requires the **in-scattering** into
group :math:`g` from all groups :math:`g'`:

.. math::
   :label: sigs-in-scatter-transpose

   \sum_{g'} \Sigs{g' \to g} \phi_{g'}
   = \bigl(\boldsymbol{\Sigma}_{\mathrm{s}}^T \cdot \boldsymbol{\phi}\bigr)_g

.. vv-status: sigs-in-scatter-transpose documented
.. Representational convention identity: the in-scatter sum equals the transpose
.. matvec Sig_s^T phi (the from-row / to-column convention). Its terminal use is
.. the removal matrix (removal-matrix), verified end-to-end by the multi-group
.. homogeneous chain (tests/gates/homogeneous/test_homogeneous.py verifies
.. "removal-matrix", >=2 groups per the ERR-002 warning). A convention identity,
.. not a separate solver claim.

This is why the removal matrix :eq:`removal-matrix` uses the
**transpose** :math:`\boldsymbol{\Sigma}_{\mathrm{s}}^T`: the transpose
converts column :math:`g` (destination) to row :math:`g`, making the
matrix-vector product give the in-scattering rate per group.

.. warning::

   Getting the transpose wrong is a common source of bugs (see
   ERR-002 in the error catalog).  For symmetric scattering matrices
   (e.g., 1-group self-scatter), the transpose is invisible, and the
   bug only manifests in multi-group problems with asymmetric
   down-scatter.  This is why verification must always include
   :math:`\geq 2` groups.

The :class:`~data.macro_xs.mixture.Mixture` stores ``SigS`` as a
list of :math:`G \times G` sparse matrices, one per Legendre order.
The :math:`P_0` component ``SigS[0]`` is used by the homogeneous
solver; the higher orders are used by transport solvers with
anisotropic scattering (SN, MoC).


Analytical Solutions
=====================

One-Group Theory
-----------------

For a single energy group, the matrices reduce to scalars.  The
scattering terms cancel (a neutron scattered in group 1 remains in
group 1), and the eigenvalue problem gives immediately:

.. math::
   :label: one-group-kinf

   \kinf = \frac{\nu \Sigf{}}{\Siga{}}

.. verifies:: one-group-kinf
   :by: orpheus.derivations.continuous.analytical.homogeneous.derive_1g

   Verified analytically (exact closed-form ratio) against the
   ``homo_1eg`` :class:`~orpheus.derivations.common.verification_case.VerificationCase`.

This is the most fundamental result in reactor physics.  It states that
:math:`\kinf` is the ratio of neutron production to neutron absorption,
which is the definition of the infinite multiplication factor :cite:`Stacey2007`.

For the connection to the **four-factor formula**: in a single-material
homogeneous medium the thermal utilisation :math:`f = 1`, the resonance
escape probability :math:`p = 1` (no spatial heterogeneity), and the
fast fission factor :math:`\varepsilon = 1`, so :math:`\kinf = \eta \cdot f
\cdot p \cdot \varepsilon = \eta`.

**Numerical example** (from :func:`orpheus.derivations.continuous.analytical.homogeneous.derive_1g`):
:math:`\Sigt{} = 1.0`, :math:`\Sigma_\mathrm{c} = 0.2`,
:math:`\Sigf{} = 0.3`, :math:`\nu = 2.5`,
:math:`\Sigs{} = 0.5` cm\ :sup:`-1`:

.. math::

   \kinf = \frac{2.5 \times 0.3}{0.2 + 0.3} = 1.500000


Two-Group Theory
-----------------

For two energy groups (fast and thermal) with downscatter only
(:math:`\chi = [1, 0]`, no upscatter from thermal to fast), the
matrices are:

.. math::
   :label: two-group-A

   \mathbf{A} = \begin{pmatrix}
     \Sigt{1} - \Sigs{1 \to 1} & 0 \\
     -\Sigs{1 \to 2} & \Sigt{2} - \Sigs{2 \to 2}
   \end{pmatrix}


.. implements:: two-group-A
   :by: orpheus.homogeneous.solver.HomogeneousProblem.loss

   **Implemented by** 3 sites. Every symbol that executes this
   equation's arithmetic is declared, not only the canonical one: a
   test is adjudicated against the transcription it actually ran, so
   declaring a single site would refute the tests that exercise the
   others.

.. implements:: two-group-A
   :by: orpheus.derivations.common.eigenvalue._infinite_medium_matrices

.. implements:: two-group-A
   :by: orpheus.derivations.continuous.analytical.homogeneous.derive_2g

.. math::
   :label: two-group-F

   \mathbf{F} = \begin{pmatrix}
     \nu_1 \Sigf{1} & \nu_2 \Sigf{2} \\
     0 & 0
   \end{pmatrix}


.. implements:: two-group-F
   :by: orpheus.transport.operators.isotropic_transfer.IsotropicFission

   **Implemented by** 3 sites. Every symbol that executes this
   equation's arithmetic is declared, not only the canonical one: a
   test is adjudicated against the transcription it actually ran, so
   declaring a single site would refute the tests that exercise the
   others. (The transport-side site was ``FissionOperator`` until CS4c
   step 4 rebound the channel — the infinite-medium problem carries no
   angular axis, so the *energy* binding is the one that executes this
   dyad; :ref:`fission-as-dyad`.)

.. implements:: two-group-F
   :by: orpheus.derivations.common.eigenvalue._infinite_medium_matrices

.. implements:: two-group-F
   :by: orpheus.derivations.continuous.analytical.homogeneous.derive_2g

Note that :math:`\mathbf{A}` is lower-triangular because there is no
upscatter (:math:`\Sigs{2 \to 1} = 0`).  This makes the inverse
analytical:

.. math::
   :label: two-group-Ainv

   \mathbf{A}^{-1} = \begin{pmatrix}
     \dfrac{1}{\Sigma_{\mathrm{r},1}} & 0 \\[8pt]
     \dfrac{\Sigs{1 \to 2}}{\Sigma_{\mathrm{r},1} \, \Sigma_{\mathrm{r},2}}
     & \dfrac{1}{\Sigma_{\mathrm{r},2}}
   \end{pmatrix}


.. implements:: two-group-Ainv
   :by: orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator

   **Implemented by** 3 sites. Every symbol that executes this
   equation's arithmetic is declared, not only the canonical one: a
   test is adjudicated against the transcription it actually ran, so
   declaring a single site would refute the tests that exercise the
   others.

.. implements:: two-group-Ainv
   :by: orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator.as_matrix

.. implements:: two-group-Ainv
   :by: orpheus.derivations.continuous.fn_method.origins.k_inf_derivations.derive_kinf_mg_matrix_form

where :math:`\Sigma_{\mathrm{r},g} = \Sigt{g} - \Sigs{g \to g}` is
the **removal cross section** for group :math:`g` (total minus
in-group scattering = absorption + out-scattering).

The eigenvalue matrix :math:`\mathbf{M} = \mathbf{A}^{-1}\mathbf{F}`
is:

.. math::
   :label: two-group-M

   \mathbf{M} = \begin{pmatrix}
     \dfrac{\nu_1 \Sigf{1}}{\Sigma_{\mathrm{r},1}}
     & \dfrac{\nu_2 \Sigf{2}}{\Sigma_{\mathrm{r},1}} \\[8pt]
     \dfrac{\Sigs{1 \to 2}\, \nu_1 \Sigf{1}}
           {\Sigma_{\mathrm{r},1}\,\Sigma_{\mathrm{r},2}}
     & \dfrac{\Sigs{1 \to 2}\, \nu_2 \Sigf{2}}
             {\Sigma_{\mathrm{r},1}\,\Sigma_{\mathrm{r},2}}
   \end{pmatrix}


.. implements:: two-group-M
   :by: orpheus.homogeneous.solver.solve_homogeneous_infinite

   **Implemented by** 3 sites. Every symbol that executes this
   equation's arithmetic is declared, not only the canonical one: a
   test is adjudicated against the transcription it actually ran, so
   declaring a single site would refute the tests that exercise the
   others.

.. implements:: two-group-M
   :by: orpheus.numerics.eigenvalue.direct_eigenvalue

.. implements:: two-group-M
   :by: orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous

Because the fission source enters only group 1 (:math:`\chi = [1, 0]`,
so the second row of :math:`\mathbf{F}` is zero), the term
:math:`\nu_2\Sigf{2}/\Sigma_{\mathrm{r},2}` does **not** appear in
:math:`M_{22}`: group 2 absorbs and down-scatters but produces no
fission emission of its own.  Consequently the two rows of
:math:`\mathbf{M}` are proportional — :math:`\mathbf{M}` is **rank 1**
(as it must be, since :math:`\mathbf{F} = \boldsymbol{\chi}\otimes\nu\Sigf{}`
is a rank-1 dyad) — and its only non-zero eigenvalue is its trace.

The characteristic equation :math:`\det(\mathbf{M} - \lambda\mathbf{I}) = 0`
gives a quadratic in :math:`\lambda`:

.. math::
   :label: two-group-charpoly

   \lambda^2 - \bigl(M_{11} + M_{22}\bigr)\lambda
   + \bigl(M_{11}M_{22} - M_{12}M_{21}\bigr) = 0


.. implements:: two-group-charpoly
   :by: orpheus.derivations.continuous.fn_method.origins.k_inf_derivations.derive_kinf_mg_matrix_form

   **Implemented by** the one site in the tree that executes this
   equation's arithmetic.

whose roots are:

.. math::
   :label: two-group-roots

   \lambda_{\pm} = \frac{(M_{11} + M_{22})
                   \pm \sqrt{(M_{11} - M_{22})^2 + 4 M_{12} M_{21}}}{2}


.. implements:: two-group-roots
   :by: orpheus.homogeneous.solver.solve_homogeneous_infinite

   **Implemented by** 2 sites. Every symbol that executes this
   equation's arithmetic is declared, not only the canonical one: a
   test is adjudicated against the transcription it actually ran, so
   declaring a single site would refute the tests that exercise the
   others. The solver executes the **rank-one specialisation**:
   :math:`\det\mathbf{M} = 0`, so the root
   :math:`\lambda_+ = \operatorname{tr}\mathbf{M} = \langle\nu\Sigma_f,
   \mathbf{A}^{-1}\chi\rangle` (derived below the worked example); the
   reference oracle runs a dense eigen-solve on the materialized
   :math:`\mathbf{M}`.

.. implements:: two-group-roots
   :by: orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous

The dominant root :math:`\lambda_+` is :math:`\kinf`.

**Worked numerical example** (from :func:`orpheus.derivations.continuous.analytical.homogeneous.derive_2g`):

.. list-table::
   :header-rows: 1
   :widths: 15 15 15 15 15 15 15

   * - :math:`g`
     - :math:`\Sigt{}`
     - :math:`\Sigma_\mathrm{c}`
     - :math:`\Sigf{}`
     - :math:`\nu`
     - :math:`\Sigs{g \to g}`
     - :math:`\Sigs{1 \to 2}`
   * - 1
     - 0.50
     - 0.01
     - 0.01
     - 2.50
     - 0.38
     - 0.10
   * - 2
     - 1.00
     - 0.02
     - 0.08
     - 2.50
     - 0.90
     - ---

The removal cross sections are :math:`\Sigma_{\mathrm{r},1} = 0.50 - 0.38 = 0.12`
and :math:`\Sigma_{\mathrm{r},2} = 1.00 - 0.90 = 0.10`.

The eigenvalue matrix entries are:

.. math::

   M_{11} &= \frac{2.50 \times 0.01}{0.12} = 0.208\overline{3} \\[4pt]
   M_{12} &= \frac{2.50 \times 0.08}{0.12} = 1.6\overline{6} \\[4pt]
   M_{21} &= \frac{0.10 \times 2.50 \times 0.01}{0.12 \times 0.10} = 0.208\overline{3} \\[4pt]
   M_{22} &= \frac{0.10 \times 2.50 \times 0.08}{0.12 \times 0.10} = 1.6\overline{6}

The two rows are identical (rank-1 :math:`\mathbf{M}`), so the
characteristic polynomial :eq:`two-group-charpoly` factors as
:math:`\lambda\,(\lambda - \operatorname{tr}\mathbf{M}) = 0` and the
dominant root :eq:`two-group-roots` is simply the trace:

.. math::

   \kinf = \lambda_+ = M_{11} + M_{22}
         = 0.208\overline{3} + 1.6\overline{6} = 1.8750000000

The second eigenvalue is :math:`\lambda_- = 0`.  This is exact and
structural, not a coincidence of these cross sections: the production
matrix :math:`\mathbf{F} = \boldsymbol{\chi} \otimes \nu\Sigma_f` is a
**rank-1 dyad** (fission emits with the single spectrum
:math:`\boldsymbol{\chi}`), so :math:`\mathbf{M} = \mathbf{A}^{-1}\mathbf{F}`
is also rank 1 — it has exactly one non-zero eigenvalue,
:math:`\kinf`, regardless of the group count.

The general root formula :eq:`two-group-roots` shows why the trace is the
answer and not an accident of these numbers. A rank-one :math:`2\times 2`
matrix has :math:`\det\mathbf{M} = M_{11}M_{22} - M_{12}M_{21} = 0`, so the
discriminant is

.. math::

   (M_{11} - M_{22})^2 + 4M_{12}M_{21}
   \;=\; (M_{11} + M_{22})^2 - 4\det\mathbf{M}
   \;=\; (\operatorname{tr}\mathbf{M})^2 ,

and the two roots are
:math:`\lambda_\pm = (\operatorname{tr}\mathbf{M} \pm
\lvert\operatorname{tr}\mathbf{M}\rvert)/2`, that is
:math:`\{\operatorname{tr}\mathbf{M},\, 0\}`. The solver reads that root
without forming :math:`\mathbf{M}` at all: with
:math:`\mathbf{M} = (\mathbf{A}^{-1}\boldsymbol{\chi})(\nu\Sigma_f)^{\mathsf T}`,
:math:`\operatorname{tr}\mathbf{M} = \langle\nu\Sigma_f,
\mathbf{A}^{-1}\boldsymbol{\chi}\rangle`, one solve and one dot product at
any group count (:ref:`homogeneous-rank-one-route`). There is no
iteration whose convergence rate would depend on a dominance ratio, and no
eigen-solver whose last bit depends on the platform.

.. note::

   The large :math:`\kinf` in these analytical benchmarks reflects the
   synthetic cross sections chosen for verification, not a physical
   reactor.  The cross sections are deliberately simple to enable exact
   symbolic solutions.


Four-Group Theory
------------------

For four groups (fast, epithermal, thermal-1, thermal-2) with a full
downscatter cascade and fission in all groups, the characteristic
polynomial is degree 4.  It needs no closed form, because
:math:`\mathbf{M} = \mathbf{A}^{-1}\mathbf{F}` is rank one exactly as in
two groups: its only non-zero root is
:math:`\operatorname{tr}\mathbf{M} = \langle\nu\Sigma_f,
\mathbf{A}^{-1}\boldsymbol{\chi}\rangle`, and for decimal cross sections that
number is **rational**.

**How the registry computes it.**
:func:`orpheus.derivations.continuous.analytical.homogeneous.derive_4g`
calls :func:`~orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous`,
which forms :math:`\mathbf{M}` with :func:`numpy.linalg.solve` and takes its
dominant eigenvalue with :func:`numpy.linalg.eig`: a float64 LAPACK
computation, not a SymPy one (`[M]` 2026-10-01, reading the function; SymPy
in that module only typesets the 1- and 2-group matrices for the generated
LaTeX).

**Result**, three readings of one quantity (`[M]` 2026-10-01):

.. list-table::
   :header-rows: 1
   :widths: 38 34 28

   * - what is computed
     - value
     - how
   * - exact, for the **decimal** data of ``derive_4g``
     - :math:`31243/21000 = 1.48776190476190476\ldots`
     - rational elimination
   * - exact, for the **float64** inputs the solvers receive
     - :math:`1.48776190476190405\ldots`
     - rational elimination
       (:mod:`orpheus.derivations.common.exact_homogeneous`)
   * - the registry's ``case.k_inf``
     - ``1.4877619047619044``
     - ``numpy.linalg.eig`` of the float :math:`\mathbf{M}`
   * - the production solver
     - ``1.487761904761904``
     - one LU solve and one dot product

The first two rows differ by about 3.2 ULP because ``0.01`` and the other
decimal cross sections are not binary numbers: the float inputs pose a
slightly different problem, and the exact answer to THAT problem is the
only thing a floating-point solver can be held to
(:ref:`homogeneous-exact-reference`). To ten digits all four read

.. math::

   \kinf = 1.4877619048

.. warning::

   The analytical eigenvalues are computed from the **same matrix
   structure** as the numerical solver.  This is a code-verification
   test (does the code correctly implement the matrix algebra?), not a
   physics-validation test.  Independent validation requires comparison
   to a different code or to experimental data (see the MATLAB
   reference values in the demo scripts).


.. _xs-preparation:

Cross-Section Preparation
==========================

Before any solver can run, the macroscopic cross sections must be
assembled from isotopic data.  This section describes the pipeline
implemented in :mod:`data.macro_xs`, which is exercised for the first
time in the homogeneous module and reused by all subsequent solvers.


.. _xs-pipeline-overview:

Pipeline Overview
------------------

The cross-section preparation follows five steps:

1. **Load isotopes** — read the 421-group microscopic cross-section
   library for each nuclide at the desired temperature.
2. **Compute number densities** — convert mass densities and
   compositions to number densities in the library's unit system.
3. **Sigma-zero iteration** — find the self-consistent background cross
   section for each isotope and group (self-shielding).
4. **Interpolate** — evaluate microscopic cross sections at the
   converged sigma-zero values.
5. **Sum to macroscopic** — weight by number densities and sum over
   isotopes to obtain the :class:`~data.macro_xs.mixture.Mixture`.

This pipeline is encapsulated in
:func:`~data.macro_xs.mixture.compute_macro_xs`.


Number Densities
-----------------

The atomic number density of species :math:`i` (in :math:`1/(\text{barn}
\cdot \text{cm})`) is:

.. math::
   :label: number-density

   N_i = \frac{\rho_i}{m_u \, A_i}


.. implements:: number-density
   :by: orpheus.data.macro_xs.recipes._number_density

   **Implemented by** the one site in the tree that executes this
   equation's arithmetic.

where :math:`\rho_i` is the partial mass density in
:math:`\text{g}/\text{cm}^3`, :math:`m_u = 1.660538 \times 10^{-24}` g
is the atomic mass unit, and :math:`A_i` is the atomic weight.  The
factor :math:`10^{-24}` converts the natural units
(:math:`\text{cm}^{-3}`) to the library units
(:math:`1/(\text{barn} \cdot \text{cm})`).

For aqueous solutions, the water density is obtained from the IAPWS-IF97
steam tables via ``pyXSteam``.  See
:func:`~data.macro_xs.recipes.aqueous_uranium` and
:func:`~data.macro_xs.recipes.pwr_like_mix`.


.. _sigma-zero-iteration:

Sigma-Zero Self-Shielding
---------------------------

Cross sections in the resonance region depend strongly on the
**background cross section** :math:`\sigma_{0,i,g}` — a measure of how
"dilute" isotope :math:`i` is relative to its neighbours.  The
background cross section is defined as :cite:`Bondarenko1964`:

.. math::
   :label: sigma-zero

   \sigma_{0,i,g}
   = \frac{\Sigma_\mathrm{escape} + \displaystyle\sum_{j \ne i}
           N_j \, \sigma_{\mathrm{t},j,g}}{N_i}

where :math:`\Sigma_\mathrm{escape}` is the escape cross section
(zero for an infinite homogeneous medium) and the sum runs over all
other isotopes in the mixture.

**Physical meaning**: when :math:`\sigma_0` is large (dilute limit or
strong moderator), the resonance peaks are fully resolved and the
effective cross section is close to the infinite-dilution value.  When
:math:`\sigma_0` is small (concentrated heavy absorber), the neutron
flux is depressed at resonance energies — **self-shielding** — and the
effective cross section is reduced.

The definition :eq:`sigma-zero` is implicit: :math:`\sigma_{\mathrm{t},j,g}`
itself depends on :math:`\sigma_{0,j,g}` through the library
interpolation tables.  The solution is obtained by **fixed-point
iteration**:

1. Initialise :math:`\sigma_0` to a large value (:math:`10^{10}` barns,
   the infinite-dilution limit).
2. Interpolate :math:`\sigma_{\mathrm{t},j,g}` from the library at the
   current :math:`\sigma_0`.
3. Recompute :math:`\sigma_0` from Eq. :eq:`sigma-zero`.
4. Repeat until :math:`\|\sigma_0^{(n)} - \sigma_0^{(n-1)}\| < 10^{-6}`.

Convergence is fast (typically 3--5 iterations) because the dependence
of :math:`\sigma_\mathrm{t}` on :math:`\sigma_0` is weak and monotonic.
This is implemented in :func:`~data.macro_xs.sigma_zeros.solve_sigma_zeros`.

.. note::

   For an **infinite homogeneous** medium, :math:`\Sigma_\mathrm{escape}
   = 0`.  The sigma-zero depends only on the other isotopes in the
   mixture.  For **heterogeneous** cells (fuel pins), the escape cross
   section :math:`\Sigma_e = \Sigma_\mathrm{pot} / \bar{\ell}
   \approx S/(4V)` accounts for spatial self-shielding via the
   equivalence theory of Bondarenko :cite:`Bondarenko1964`.


Cross-Section Interpolation
-----------------------------

The 421-group library tabulates microscopic cross sections at discrete
:math:`\sigma_0` base points (e.g., :math:`10^0, 10^1, \ldots, 10^{10}`
barns).  Once the sigma-zero iteration converges, the cross section at
the converged :math:`\sigma_0` is obtained by **log-linear interpolation**
in :math:`\log_{10}(\sigma_0)` space:

.. math::
   :label: xs-interp

   \sigma_{x,g}(\sigma_0) \approx \sigma_{x,g}(\sigma_0^{(a)})
   + \frac{\log_{10} \sigma_0 - \log_{10} \sigma_0^{(a)}}
          {\log_{10} \sigma_0^{(b)} - \log_{10} \sigma_0^{(a)}}
     \bigl[\sigma_{x,g}(\sigma_0^{(b)}) - \sigma_{x,g}(\sigma_0^{(a)})\bigr]

where :math:`\sigma_0^{(a)}` and :math:`\sigma_0^{(b)}` are the
bracketing base points.  This is performed by
:func:`~data.macro_xs.interpolation.interp_xs_field` for scalar
cross sections and
:func:`~data.macro_xs.interpolation.interp_sig_s` for scattering
matrices.


Macroscopic Summation
----------------------

The macroscopic cross section for reaction :math:`x` in group :math:`g`
is the density-weighted sum over all isotopes:

.. math::
   :label: macro-sum

   \Sigma_{x,g} = \sum_{i=1}^{I} N_i \, \sigma_{x,i,g}

The following reaction types are assembled:

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - Attribute
     - Reaction
     - Notes
   * - ``SigC``
     - :math:`(n,\gamma)` capture
     - Radiative capture
   * - ``SigL``
     - :math:`(n,\alpha)` loss
     - Charged-particle emission
   * - ``SigF``
     - :math:`(n,f)` fission
     - Fission cross section
   * - ``SigP``
     - Production
     - :math:`\nu\Sigf{}`, summed over **fissile** isotopes only
   * - ``SigS``
     - Scattering matrices
     - One :math:`G \times G` sparse matrix per Legendre order
   * - ``Sig2``
     - :math:`(n,2n)` matrix
     - :math:`G \times G` sparse transfer matrix
   * - ``SigT``
     - Total
     - :math:`\Sigma_\mathrm{c} + \Sigma_\mathrm{L} + \Sigma_\mathrm{f}
       + \text{rowsum}(\Sigma_\mathrm{s}^{P_0})
       + \text{rowsum}(\Sigma_2)`
   * - ``chi``
     - Fission spectrum
     - Taken from first fissile isotope (simplification)

The **absorption cross section** — the diagnostic one-group balance
ratio reported alongside :math:`\kinf` (and the denominator of the
classical production/absorption form of the eigenvalue, Eq.
:eq:`keff-update`) — is not stored directly but computed as a derived
property (:attr:`~data.macro_xs.mixture.Mixture.absorption_xs`):

.. math::
   :label: absorption-xs

   \Siga{g} = \Sigf{g} + \Sigma_{\mathrm{c},g}
            + \Sigma_{\mathrm{L},g}
            + \text{rowsum}(\boldsymbol{\Sigma}_{2,g})

This includes fission (neutron is absorbed to produce fission
fragments), radiative capture :math:`(n,\gamma)`, charged-particle
emission :math:`(n,\alpha)`, and the :math:`(n,2n)` reaction (where
one neutron is "absorbed" and two are emitted, for a net gain of one).

The result is stored in a :class:`~data.macro_xs.mixture.Mixture`
dataclass, which is the universal input to all ORPHEUS solvers.


Neutron Spectrum Physics
=========================

The shape of the neutron energy spectrum :math:`\phi(E)` in a
homogeneous medium is controlled by the competition between
moderation (slowing-down) and absorption.  Three distinct energy
regions are visible in the spectrum plots:

Fast Region (:math:`E > 0.1` MeV)
-----------------------------------

Neutrons are born in fission with a spectrum peaked around 2 MeV.
At these energies, scattering is nearly isotropic in the
centre-of-mass frame and the mean logarithmic energy loss per
collision with hydrogen is :math:`\xi = 1`.  The fission source
produces the characteristic fast peak.

For heavy nuclei like :sup:`238`\ U, the elastic energy loss
per collision is very small (:math:`\xi \approx 2/A`), so
the spectrum in the fast range is close to the fission spectrum
:math:`\chi(E)`.

Slowing-Down Region (:math:`1\;\text{eV} < E < 0.1\;\text{MeV}`)
-------------------------------------------------------------------

In this intermediate range, neutrons are slowed by elastic
scattering (primarily with hydrogen).  In the absence of absorption,
the slowing-down equation yields the well-known **1/E flux** law:

.. math::
   :label: one-over-E

   \phi(E) = \frac{S}{\xi \Sigt{}} \cdot \frac{1}{E}

.. vv-status: one-over-E documented

where :math:`S` is the slowing-down source (neutrons entering from
above) and :math:`\xi` is the mean logarithmic energy decrement.
On a **flux-per-lethargy** plot, the 1/E region appears as a
horizontal plateau:

.. math::
   :label: flux-per-lethargy-plateau

   \frac{\phi}{du} = \frac{\phi(E) \cdot E}{\Delta u}
   \propto \frac{1}{E} \cdot E = \text{const}

.. vv-status: flux-per-lethargy-plateau documented
.. Definitional physics identity: the 1/E slowing-down flux appears flat on a
.. per-lethargy plot (phi/du ~ (1/E)*E = const), the plotting-convention sibling
.. of the 1/E law (one-over-E). A spectral-physics teaching identity, not a
.. solver claim.

This is why flux-per-lethargy is the standard representation: it
makes the slowing-down region flat, and deviations (resonance dips
from :sup:`238`\ U, thermal peak) are immediately visible.

Resonance absorption (:sup:`238`\ U capture resonances) creates
**dips** in the spectrum throughout this range.  The sigma-zero
self-shielding (see :ref:`sigma-zero-iteration`) accounts for the
flux depression in the resonance peaks.

Thermal Region (:math:`E < 1` eV)
------------------------------------

Below about 1 eV, neutrons reach thermal equilibrium with the
moderator atoms.  The thermal flux approaches a **Maxwell–Boltzmann
distribution** at the moderator temperature :math:`T`:

.. math::
   :label: maxwellian

   \phi_\mathrm{th}(E) \propto E \, \exp\!\left(-\frac{E}{k_B T}\right)

.. vv-status: maxwellian documented

which peaks at :math:`E_\mathrm{peak} = k_B T`.  At room temperature
(294 K), :math:`k_B T = 0.0253` eV, producing the characteristic
thermal peak near 0.025 eV.

At higher moderator temperatures (e.g., 600 K in a PWR), the peak
shifts to higher energies and broadens — this is Doppler broadening
of the moderator distribution, which affects the thermal spectrum
shape and hence the thermal cross sections.

Absorber poisons (e.g., boron) selectively remove thermal neutrons,
depressing the thermal peak.  This is clearly visible comparing the
aqueous spectrum (no boron, strong thermal peak) with the PWR-like
spectrum (4000 ppm B, suppressed thermal peak) in the
:ref:`example-problems` section below.


.. _power-iteration-algorithm:
.. _direct-eigensolve:

The Eigenvalue Solution
=======================

The eigenvalue problem :eq:`matrix-eigenvalue`,
:math:`\mathbf{A}\boldsymbol{\phi} = \tfrac{1}{k}\mathbf{F}\boldsymbol{\phi}`,
asks for the dominant eigenpair :math:`(\kinf, \boldsymbol{\phi})`.  How
it is solved depends on whether the problem couples space:

- **Spatially-coupled solvers** (SN, CP, MoC, diffusion) cannot afford a
  dense inverse of the full loss operator, so they sweep/solve the loss
  once per outer step and drive :math:`k` up the dominant mode by
  **power iteration** on the fission source :cite:`Hebert2009`.  Those
  realisations live in the spatial theory pages (e.g.
  :eq:`cp-keff-update`, :eq:`moc-keff-update`); this section is the
  shared conceptual hub they cross-reference.
- **The infinite homogeneous medium** has no spatial coupling — the loss
  matrix :math:`\mathbf{A}` is a single :math:`G \times G` dense block —
  so the eigenpair is taken **directly**, with no iteration, by
  :func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite`: because
  the fission operator has rank one, the eigenpair is one linear solve and
  one inner product (:ref:`homogeneous-rank-one-route`).

The remainder of this section describes that direct solve.


.. _direct-eigensolve-assembly:

Assembling the loss matrix from the transport operators
-------------------------------------------------------

The defining design decision of campaign **#276** is that the
infinite-medium loss matrix is **not** a bespoke energy matrix — it is
the meshed SN solver's own loss operator
:math:`\mathbf{A} = C - K_\mathrm{iso}` evaluated on the degenerate
phase space :math:`V_E \otimes V_{\rm pt}`.  This is Cardinal Rule 2
(cross-model single source) applied to the simplest model in the
curriculum: there is exactly one place in ORPHEUS where the isotropic
in-scatter source :math:`\Sigma_{s0}^T\phi + 2\Sigma_2^T\phi` is
assembled, and the homogeneous problem reuses it rather than
re-implementing the same algebra.

The construction proceeds in five steps, and since the CS4c coda
(2026-09-08) all five live on the problem's own hub,
:class:`~orpheus.homogeneous.solver.HomogeneousProblem`, as cached
per-instance state; the solver only reads them.  Steps 1 and 2 are
listed separately because they answer two different questions —
*which space?* and *which numbers?* — but they now have **one** answer
between them: **the mixture determines both**, and each consumed object
is minted from it and from nothing else.

1. **Pose the space.**
   :func:`~orpheus.homogeneous.solver._pose_space` mints
   :math:`V_E \otimes V_{\rm pt}` from the **mixture** — the energy axis
   through the one energy-arm rule
   (:meth:`EnergyAxis.from_materials
   <orpheus.numerics.axis.EnergyAxis.from_materials>`, the same rule
   :attr:`MaterialMesh.bulk_space
   <orpheus.transport.mesh.material_mesh.MaterialMesh.bulk_space>`
   routes through, so the two spellings cannot diverge) tensored with
   the explicit **quotient point**, a one-element spatial axis carrying
   the COUNTING weight. That weight *is* the normalized "per unit
   volume" density convention of clause 1
   (:ref:`spaces-quotient-family`), and it is what the post-processing
   reaction-rate pairings consume
   (:ref:`homogeneous-rates-and-normalisation`).

   ⚠ **This mint is the whole story, and it has been since two
   separate corrections.** Campaign 1 CS4a (**K2**) stopped the SPACE
   being read off a carrier; the CS4c coda (**R-c1**) stopped the DATA
   being read off one. A page that says the space comes off a carrier is
   describing the pre-K2 tree; a page that says a carrier supplies the
   cross sections is describing the pre-coda tree. What survives of the
   old reading is a *reference*, not a source: a genuine unit-width
   one-cell ``Mesh1D`` carrier's
   :attr:`~orpheus.transport.mesh.material_mesh.MaterialMesh.bulk_space`
   still mints a space that is ``==``
   :attr:`HomogeneousProblem.space
   <orpheus.homogeneous.solver.HomogeneousProblem.space>`, and the
   identity-bridge gate (``tests/gates/homogeneous/test_operator_spaces.py``
   G2.1) keeps that equality honest — the two spellings route through
   the same energy-arm rule and the same one-cell volume, so they cannot
   silently diverge.

2. **Supply the cross sections — from the mixture, onto the pose.**
   The hub mints two tiers directly from the
   :class:`~orpheus.data.macro_xs.mixture.Mixture`.  The **kernel** tier
   is the representation-free material data — the scattering and
   :math:`(n,2n)` Legendre transfer stacks
   (:attr:`~orpheus.homogeneous.solver.HomogeneousProblem.scattering`,
   :attr:`~orpheus.homogeneous.solver.HomogeneousProblem.n2n`) and the
   fission datum :math:`\chi \otimes \nu\Sigma_f`
   (:attr:`~orpheus.homogeneous.solver.HomogeneousProblem.fission`), each
   carried over the one-cell
   :attr:`~orpheus.homogeneous.solver.HomogeneousProblem.layout`
   ``{0: (arange(1),)}`` — material ``0`` on cell ``0``.  The **field**
   tier is the three
   :class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`
   objects :math:`\Sigma_t`, :math:`\Sigma_a` and :math:`\nu\Sigma_f`
   (:attr:`~orpheus.homogeneous.solver.HomogeneousProblem.total_cross_section_field`,
   :attr:`~orpheus.homogeneous.solver.HomogeneousProblem.absorption_cross_section_field`,
   :attr:`~orpheus.homogeneous.solver.HomogeneousProblem.fission_production_field`),
   each reshaped to :math:`(n_g, 1)` and **born on the pose of step 1**.

   *Born*, not re-posed, is the load-bearing word.  A field carries its
   own space, and that space is the measure authority every downstream
   pairing reads (:ref:`homogeneous-rates-and-normalisation`).  When the
   fields were minted elsewhere and rebound afterwards, "is this field on
   the right space?" was a question someone had to keep answering
   correctly; minting them on the pose makes the wrong answer
   *unspellable* — there is no other space in scope to mint them on.
   That matters more than it looks, because
   :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
   does **not** validate its coefficient's space against its own ends
   (`[M]` zero ``coefficient.space`` reads in its body), which is exactly
   why a carrier-minted coefficient could ride the collision diagonal for
   years without any gate objecting.

   .. note::

      **What this replaced, and why it went.** Until the coda a
      single-cell, single-region **mesh-less** ``MaterialMesh`` was
      fabricated by a ``from_materials`` factory, and its
      ``material_xs_field()`` facade supplied the operators.  That
      carrier existed only to be a shape the data could hang on: its
      ``[0, 1]`` edges, its node at ``0.5`` and its Cartesian chart were
      read by nothing on the path — the falsifiable tell of objective O1
      ("no fabricated data reaches the operators"), and the reason the
      retirement is a *re-source* rather than a re-baseline (`[M]` the
      byte gate of the day read 8 of 8 across the change).  The factory retired with
      its last consumer, and three arms that existed only to serve it
      went with it: ``areas``' "no faces at all" case, the
      :math:`S_N`-promotion refusal and the diffusion bounded-geometry
      refusal.  Each had become **input-less**: `[M]` the carrier hierarchy
      has exactly one producer of ``mesh = None``
      (``SNProblem.from_axes``, above :math:`d = 2`) and every other
      constructor requires a mesh, so with the homogeneous path no
      longer building one there is no mesh-less :math:`d \leq 2`
      carrier left for those arms to receive.  ``mesh is None``
      therefore carries ONE meaning again — the :math:`d \geq 3`
      axis-native carrier — and CS4b S7's ``ndim`` discrimination has
      nothing left to discriminate
      (:ref:`homogeneous-development-history`).

3. **Collision diagonal** :math:`C = \mathrm{diag}(\Sigma_t)`, the
   :attr:`~orpheus.homogeneous.solver.HomogeneousProblem.collision`
   binding of
   :attr:`~orpheus.homogeneous.solver.HomogeneousProblem.total_cross_section_field`
   with the pose as both domain and codomain.

4. **Isotropic energy transfer**
   :math:`K_\mathrm{iso} = \Sigma_{s0}^T + 2\Sigma_2^T`, the action of the
   two model-shared operators

   .. math::
      :label: fission-source

      K_\mathrm{iso} \;=\;
      \underbrace{\Sigma_{s0}^{T}}_{\text{\scriptsize :class:`IsotropicScattering`}}
      \;+\;
      \underbrace{2\,\Sigma_2^{T}}_{\text{\scriptsize :class:`IsotropicN2N`}}

   where :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicScattering`
   realises :math:`\Sigma_{s0}^T` (the in-scatter source matrix, the stored
   ``[g_from, g_to]`` transfer transposed) and
   :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicN2N`
   realises :math:`2\Sigma_2^T` (the loss-side multiplicity-2 transfer).
   The composed loss operator :math:`\mathbf A = C - K_\mathrm{iso}` is
   the hub's
   :attr:`~orpheus.homogeneous.solver.HomogeneousProblem.loss`, and it is
   held **un-materialized** (an
   :class:`~orpheus.numerics.operator.OperatorSum`) — the consumer chooses
   the realization (taxonomy step 5b).  Because every arm poses on the
   same space, the ``OperatorSum`` guard *validates* the sum rather than
   skipping it.  Its dense :math:`(n_g, n_g)` form is **not** assembled
   term-by-term from per-material blocks; it is produced
   one layer later, by the operator's own
   :meth:`~orpheus.numerics.operator.LinearOperator.as_matrix` apply-to-basis
   (:ref:`matrix-inverse-operator`) on the one-cell pose — the
   ``(ng, 1)`` basis shape **derived from the operators' threaded domain**
   (the mixture-minted pose of step 1; campaign 1 CS1 gave the operators a
   real space to derive it from, CS4a K2 made that space the problem's
   rather than a carrier's, and the CS4c coda made the hub its owner —
   before CS1, every consumer passed ``basis_shape=(ng, 1)`` by hand) —
   **inside the**
   :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
   **constructor** (one eager materialization + LU factorization, then one
   backsolve against :math:`\chi`; see :ref:`direct-eigensolve-solve`). (The operators'
   :meth:`~orpheus.transport.operators.isotropic_transfer.IsotropicScattering.dense_per_material`
   accessor — the transpose read straight off the stored cross sections — is
   a storage-side *oracle* used by the verification gates as a
   structurally-independent cross-check, **not** a production assembly path.)

5. **Drop streaming.**  In an infinite medium the streaming operator
   :math:`L` is identically zero (:math:`\nabla\psi = 0`), so it is
   omitted from the sum.  What remains,
   :math:`\mathbf{A} = C - K_\mathrm{iso} = \mathrm{diag}(\Sigma_t)
   - \Sigma_{s0}^T - 2\Sigma_2^T`, is exactly the removal matrix
   :eq:`removal-matrix`.

.. note::

   The label :eq:`fission-source` historically named the per-iteration
   fission source of the retired power iteration,
   :math:`\mathbf{Q}_f = (\boldsymbol{\chi}/k)\,\nu\Sigma_f\cdot\boldsymbol{\phi}`.
   Under the direct method there is no iterate :math:`k^{(n)}` and no
   reassembled source; the production source is the single application of
   the dyad :math:`\mathbf{F}\boldsymbol{\phi}` (see
   :ref:`direct-eigensolve-solve`).  The label is retained on the
   isotropic-transfer assembly :math:`K_\mathrm{iso}` — the energy
   redistribution that *was* the in-scatter half of the old source — so
   the verification edge from the homogeneous test suite continues to pin
   the operator algebra that produces the source.


.. _direct-eigensolve-solve:

The fission dyad and the rank-one solve
---------------------------------------

The production matrix is the rank-1 dyad
:math:`\mathbf{F} = \boldsymbol{\chi} \otimes \nu\Sigma_f`
:eq:`fission-matrix`, so solving the loss out of it is solving the loss out
of ONE column:

.. math::
   :label: fixed-source-solve

   \mathbf{M} \;=\; \mathbf{A}^{-1}\mathbf{F}
   \;=\; \mathbf{A}^{-1}\,\bigl(\boldsymbol{\chi}\otimes\nu\Sigma_f\bigr)
   \;=\; \bigl(\mathbf{A}^{-1}\boldsymbol{\chi}\bigr)\otimes\nu\Sigma_f

i.e. the loss matrix is **solved out** of the production once, against the
emission spectrum :math:`\boldsymbol{\chi}` that spans the range of
:math:`\mathbf{F}` (one solve, not :math:`G`, and never an explicit
:math:`[\mathbf{A}^{-1}]`), giving the column factor
:math:`\mathbf{u} = \mathbf{A}^{-1}\boldsymbol{\chi}` of the
:math:`G \times G` eigenvalue matrix :math:`\mathbf{M}`.  The eigenpair
follows directly:

.. math::
   :label: keff-update

   \kinf \;=\; \lambda_{\max}(\mathbf{M})
   \;=\; \bigl\langle \nu\Sigma_f,\, \mathbf{A}^{-1}\boldsymbol{\chi} \bigr\rangle,
   \qquad
   \boldsymbol{\phi} \;\propto\; \mathbf{A}^{-1}\boldsymbol{\chi},


.. implements:: keff-update
   :by: orpheus.homogeneous.solver.solve_homogeneous_infinite

   **Implemented by** 4 sites. Every symbol that executes this
   equation's arithmetic is declared, not only the canonical one: a
   test is adjudicated against the transcription it actually ran, so
   declaring a single site would refute the tests that exercise the
   others. The solver and the exact rational reference execute the
   right-hand form (one solve, one inner product); the dense engine and
   the registry's oracle execute the middle form (an eigen-solve of the
   materialized :math:`\mathbf{M}`; ``direct_eigenvalue`` delegates the
   extraction to
   :func:`~orpheus.numerics.eigenvalue.dominant_eigenpair`).

.. implements:: keff-update
   :by: orpheus.numerics.eigenvalue.direct_eigenvalue

.. implements:: keff-update
   :by: orpheus.derivations.common.exact_homogeneous.exact_pencil_eigenpair

.. implements:: keff-update
   :by: orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous

the dominant eigenpair, which for a physical medium is the unique
non-negative solution — the **fundamental mode** of the Perron–Frobenius
theorem :cite:`Hebert2009`.  Here the theorem is not invoked to *select*
the mode: the rank-one structure *constructs* it, and the next subsection
proves that it is the dominant one and that it lies in the positive cone.

.. _homogeneous-rank-one-route:

Why one solve is the whole eigenproblem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The claim is that :eq:`keff-update` is exact — not an approximation of the
eigenvalue, not the first step of an iteration — whenever
:math:`\mathbf{F}` is the dyad :eq:`fission-matrix`.  Four steps.

**1. The operator has rank one.**  For every group vector
:math:`\mathbf{x}`, :math:`\mathbf{F}\mathbf{x} = \boldsymbol{\chi}\,
\langle\nu\Sigma_f, \mathbf{x}\rangle`: fission collects ONE number, the
production rate, and emits it with ONE spectrum.  So

.. math::

   \mathbf{K}\mathbf{x} \;=\; \mathbf{A}^{-1}\mathbf{F}\mathbf{x}
   \;=\; \mathbf{u}\,\langle\nu\Sigma_f, \mathbf{x}\rangle,
   \qquad \mathbf{u} = \mathbf{A}^{-1}\boldsymbol{\chi},

that is :math:`\mathbf{K} = \mathbf{u}\,(\nu\Sigma_f)^{\mathsf T}`, whose
range is the line spanned by :math:`\mathbf{u}`.

**2. Its one non-zero eigenpair.**  Apply :math:`\mathbf{K}` to
:math:`\mathbf{u}` itself:
:math:`\mathbf{K}\mathbf{u} = \mathbf{u}\,\langle\nu\Sigma_f,\mathbf{u}\rangle`,
so :math:`\mathbf{u}` is an eigenvector with eigenvalue
:math:`k = \langle\nu\Sigma_f, \mathbf{u}\rangle`.  Conversely, if
:math:`\mathbf{K}\mathbf{x} = \lambda\mathbf{x}` with :math:`\lambda \neq 0`
then :math:`\mathbf{x} = \mathbf{K}\mathbf{x}/\lambda` lies in the range of
:math:`\mathbf{K}`, the line through :math:`\mathbf{u}`; so every eigenvector
with a non-zero eigenvalue is a multiple of :math:`\mathbf{u}`, and its
eigenvalue is :math:`k`.  The characteristic polynomial is
:math:`\lambda^{G-1}(\lambda - k)`, which is the general-:math:`G` form of the
two-group factorisation :math:`\lambda(\lambda - \operatorname{tr}\mathbf{M})`
above, and :math:`k = \operatorname{tr}\mathbf{K}`.

**3. It is the dominant one, and the mode is in the cone.**  The loss
matrix :math:`\mathbf{A} = \operatorname{diag}(\Sigma_t) - \Sigma_{s0}^{\mathsf T}
- 2\Sigma_2^{\mathsf T}` has non-positive off-diagonal entries
:math:`-\Sigma_{s0}(g'\!\to g) - 2\Sigma_2(g'\!\to g)`, so it is a
Z-matrix.  Column :math:`g'` of :math:`\mathbf{A}` sums to
:math:`\Sigma_{t,g'} - \sum_g \bigl[\Sigma_{s0}(g'\!\to g) +
2\Sigma_2(g'\!\to g)\bigr]`: one collision in group :math:`g'` removes one
neutron and returns that many to the scattering system.  When every column
sum is positive (each collision returns fewer neutrons than it removes —
the sufficient condition of strict column diagonal dominance),
:math:`\mathbf{A}` is a non-singular M-matrix, and a non-singular M-matrix
is inverse-positive, :math:`\mathbf{A}^{-1} \geq 0` entrywise (Berman and
Plemmons, *Nonnegative Matrices in the Mathematical Sciences*, SIAM 1994,
Chapter 6).  With :math:`\boldsymbol{\chi} \geq 0` and
:math:`\nu\Sigma_f \geq 0` this gives :math:`\mathbf{u} \geq 0` and
:math:`k \geq 0`, strictly positive as soon as a fissile group receives
flux.  Every other eigenvalue is :math:`0`, so :math:`k` is the spectral
radius: the dominant eigenvalue, and :math:`\mathbf{u}` is the
non-negative fundamental mode.  The flux is the gauged representative of
that ray, :math:`\boldsymbol{\phi} = 100\,\mathbf{u}/k`
(:eq:`normalisation`).

**4. What is and is not checked at run time.**  The dense route's guards
have nothing left to act on: an eigen-solver can return a complex dominant
eigenvalue or an eigenvector of either sign, which is why
``dominant_eigenpair`` rejected the first and sign-normalised the second;
here :math:`k` is a real inner product and the sign of :math:`\mathbf{u}`
is fixed by :math:`\boldsymbol{\chi}`, not by a solver's arbitrary
orientation.  The cone membership of step 3 is a theorem for a physical
:math:`\mathbf{A}`, and `[M]` 2026-10-01 nothing in
:func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite` asserts it:
an unphysical loss matrix would return whatever
:math:`\langle\nu\Sigma_f, \mathbf{A}^{-1}\boldsymbol{\chi}\rangle` is.  A
non-fissile mixture (:math:`\nu\Sigma_f = 0`) reads :math:`k = 0`, which
the gauge's zero-reading refusal stops before a division by zero
(:class:`~orpheus.numerics.gauge.ScaleGauge`).

**The rank one is the model's.**  Every step above rests on
:math:`\mathbf{F}` being one column times one row, and that is a property
of the cross-section data, not of the transport physics.  A mixture
carries a single fission spectrum: the production-weighted average of its
isotopes' spectra, taken at flat flux.  The exact operator over :math:`K`
fissile isotopes is :math:`\mathbf{X}\mathbf{N}^{\mathsf T}` (spectra as the
columns of :math:`\mathbf{X}`, production cross sections as the columns of
:math:`\mathbf{N}`), of rank up to :math:`K`, and its :math:`k` is the
dominant eigenvalue of the :math:`K \times K` matrix
:math:`\mathbf{N}^{\mathsf T}\mathbf{A}^{-1}\mathbf{X}` — the same
construction one dimension up
(`#549 <https://github.com/deOliveira-R/ORPHEUS/issues/549>`_).  On this
path a non-rank-one :math:`\mathbf{F}` cannot be spelled: the fission kernel
the hub builds requires one-dimensional factors, so the theorem's
hypothesis is enforced by a type rather than assumed by the solver.

**Why no eigen-solver.**  The last bit of a dense eigen-solve belongs to
the platform's linear-algebra library, not to the problem.  `[M]`
2026-09-30: after macOS 27.0.1 replaced Accelerate, :func:`numpy.linalg.eig`
(LAPACK ``geev``) returned a :math:`\kinf` 1 ULP away from the previous
release's on two of the eight shipped producing mixtures (``homo_4eg``,
``mixture_A_4g``) with no ORPHEUS change (the evidence entry is "AP38
platform bit pin" on :doc:`/development/evidence/vv-anti-patterns`).
The rank-one route has no eigen-solver, so the only platform primitive on
the path is one LU solve (``getrf`` then ``getrs``), whose bytes are
deliberately not pinned (:ref:`homogeneous-exact-reference`).  It is also
the more accurate route: `[M]` 2026-10-01, against the exact rational
:math:`\kinf` of the float inputs over the eight mixtures, a ``geev``
eigen-solve of the materialized :math:`[\mathbf{K}]` reads at most
:math:`+1.55` ULP and the rank-one route at most :math:`0.86` ULP.

The solver's own body is now four lines — build the problem, take the
fission operator, solve the loss out of its emission spectrum, contract with
its production rate:

.. code-block:: python

   problem = HomogeneousProblem(mix)
   fission = problem.production                                             # F = |χ⟩⟨νΣf|
   phi = MatrixInverseOperator(problem.loss).apply(fission.emission_spectrum)  # u = A⁻¹χ, the Strategy
   k_inf = float(fission.production_rate.evaluate(phi).item())             # k∞ = ⟨νΣf, u⟩

:attr:`IsotropicFission.emission_spectrum
<orpheus.transport.operators.isotropic_transfer.IsotropicFission.emission_spectrum>`
is the dyad's column factor and
:attr:`~orpheus.transport.operators.isotropic_transfer.IsotropicFission.production_rate`
its row factor; the operator's
:attr:`~orpheus.transport.operators.isotropic_transfer.IsotropicFission.kernel`
is built from the same two objects, so the solver reads the very arrays the
dyad applies.  The production rate is evaluated on the posed
:math:`(n_g, 1)` column that ``apply`` returns, which matters: the
reaction-rate functional silently broadcasts a flat vector
(:ref:`homogeneous-rates-and-normalisation`).  The two factors are spelled
from the operator's own factors, never from ``mix.chi`` and
``mix.nu_sig_f``, so k is the contraction of the dyad the hub posed — and it
is measure-free: the pose's point weight enters neither factor
(`[M]` the ``test_operator_spaces`` G2.5 leg, which doubles the point weight,
keeps :math:`\kinf` unchanged).

The pencil — the problem's terminal object, the pair
:math:`(\mathbf{A}, \mathbf{F})` — is posed on the one space
:attr:`~orpheus.homogeneous.solver.HomogeneousProblem.space`:

.. code-block:: python

   # orpheus/homogeneous/solver.py — HomogeneousProblem
   loss           = collision - (isotropic_scattering + isotropic_n2n)   # A = C − K_iso
   production     = IsotropicFission(fission, domain=space, codomain=space)   # F = χ ⊗ νΣ_f
   pencil         = OperatorPencil(loss, production)                     # the terminal object
   eigen_posing   = EigenPosing(pencil, K_MAP)                           # the QUESTION

⭐ **Since 2026-09-13 the pencil is a TYPE, not a manner of speaking.**
The last two lines are new: :attr:`HomogeneousProblem.pencil
<orpheus.homogeneous.solver.HomogeneousProblem.pencil>` is an
:class:`~orpheus.numerics.pencil.OperatorPencil` and
:attr:`~orpheus.homogeneous.solver.HomogeneousProblem.eigen_posing` an
:class:`~orpheus.numerics.posing.EigenPosing` — and they are **the same
two types** the S\ :sub:`N` hub poses over its own operators
(:ref:`sn-the-problem-poses-its-pencil`), which is why they live in
:mod:`orpheus.numerics` and not in either method's package.  The
addition is purely additive: the solver reads ``problem.loss`` and
``problem.production``, the pencil's two members, and no arithmetic of
:math:`k_\infty` passes through the pencil type.  The full articulation — what a pencil IS,
its degree contract, and why it carries no inverse — is
:ref:`the-operator-pencil`.

This hub is also where the type is cheapest to *interrogate*, because
its operators are dense enough to materialise.  ``[M]`` over the four
shipped mixture families at 1, 2 and 4 groups,
:attr:`~orpheus.numerics.pencil.OperatorPencil.rhs_rank` reads **1** for
the fissile ``A`` family and **0** for ``B``/``C``/``D`` — the
Weierstrass–Kronecker count of *finite* eigenvalues, which for a rank-1
:math:`\mathbf{F}` — the dyad :eq:`fission-matrix` — says the whole
k-problem is **one-dimensional**, on
:math:`\operatorname{range}\mathbf{F}`.  The closed form that follows,
:math:`k_\infty = \operatorname{tr}(\mathbf{A}^{-1}\mathbf{F}) = \langle
\nu\Sigma_f,\, \mathbf{A}^{-1}\chi\rangle`, is what the pencil's own gate
uses as its structurally-independent reference against
:func:`~orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous`
(``[M]`` bit-exact at 1g and 2g; **one ulp** apart at 4g — absolute
:math:`2.2\times10^{-16}` on :math:`k_\infty = 1.4878`, relative
:math:`1.5\times10^{-16}` — so it is pinned at ``rtol=1e-13``, never
``array_equal``).  The same closed form is the production route
(:ref:`homogeneous-rank-one-route`).

The resolvent :math:`\mathbf{K} = \mathbf{A}^{-1}\mathbf{F}` is NOT a
property of the problem: *how* the pencil is inverted is a Strategy
choice (the consumers campaign's ruling R-cc2/R-cc5, 2026-09-12 — a
Problem builds its pencil as its last step; a resolvent picks an
inversion of it), so the inversion lives in the solver, which inverts
:math:`\mathbf{A}` against the one column :math:`\boldsymbol{\chi}` and
composes no :math:`\mathbf{K}` at all. (Until
2026-09-12 the hub carried it as a ``multiplication`` property — its own
docstring already called the explicit inverse "the strategy choice".)
⚠ And the pencil landing did **not** change that: it deliberately
carries no ``inverse`` and no ``resolvent`` method, precisely so that
this separation cannot be undone by a convenience property.

(The operators pose on the MIXTURE-MINTED Energy ⊗ point space — the
problem's own physics names its space, :doc:`spaces` — and since the coda
the *fields they carry* are born on it too, so no rebinding step stands
between the mint and the operator.  ``MatrixInverseOperator``
and ``as_matrix`` **derive** the ``(ng, 1)`` basis shape from the
threaded domain; the pre-CS1 idiom passed ``basis_shape=(ng, 1)``
explicitly at both sites because the operators of that era carried no
space to derive it from.)

.. note::

   **Which fission binding the production factor is (CS4c step 4,
   2026-08-30).**  The line above named
   :class:`~orpheus.transport.operators.fission.FissionOperator` until
   step 4 rebound the channel as **two bindings of one datum**
   (:ref:`fission-as-dyad`): the *energy* binding
   :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicFission`
   — the rank-1 dyad on the scalar flux, and the one this solver, the
   S\ :sub:`N` k-outer and the 1-D diffusion solver all consume — and
   the *angular* binding ``FissionOperator``, the frame's
   :math:`\ell = 0` conjugation of the same dyad on a posed angular
   composite, which only a solver with an angular axis can pose.
   (⚠ Since the consumers campaign's step 2, 2026-09-13, the instance
   S\ :sub:`N`'s k-outer consumes is the ``isotropic_energy`` face
   *derived from* its hub's one angular :math:`F`, not a separate mint —
   :ref:`sn-one-fission-per-problem`.  This solver and diffusion still
   bind the energy operator directly, because neither has an angular
   composite to derive it from.)
   The infinite-medium problem has no angular axis, so the energy
   binding is not merely sufficient here — it is the honest one, and
   the angular binding now **refuses** a scalar carrier at construction
   rather than silently accepting it.  The arithmetic is unchanged: the
   dyad, its ``outer`` reduction order, and therefore :math:`k_\infty`
   are the same object under both names.

:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
materializes and LU-factors the loss operator **once** at construction
(one :func:`scipy.linalg.lu_factor` of the :math:`G \times G` block,
produced by the operator's own
:meth:`~orpheus.numerics.operator.LinearOperator.as_matrix`), and its
:meth:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator.apply`
is one LU backsolve against the held factors.  The solver applies it once,
to the emission spectrum: :math:`\mathbf{u} = \mathbf{A}^{-1}\boldsymbol{\chi}`,
the loss **solved out** of the production's one column, never inverted into
an explicit :math:`[\mathbf{A}^{-1}]` and never multiplied into a
:math:`[\mathbf{K}]`.  The production-rate co-vector then reads
:math:`\kinf = \langle\nu\Sigma_f, \mathbf{u}\rangle`.  The route is
gated as a route, not only by its output:
``tests/gates/homogeneous/test_homogeneous.py::test_the_rank_one_route_is_the_homogeneous_call_path``
replaces every dense eigen driver (``dominant_eigenpair`` and the
``numpy``/``scipy`` ``eig``/``eigvals``) by a decoy that raises and requires
the solve to succeed with the same answer, and replaces
``MatrixInverseOperator.apply`` by a decoy returning :math:`2\mathbf{A}^{-1}x`
and requires :math:`\kinf` to double exactly and the gauged flux to stay
bit-identical (a power-of-two scale commutes with every rounding).

**Explicit direct realization — the strategy choice as a type.**  The
homogeneous solver is the **first production consumer** of
:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
(taxonomy step 5b).  It constructs the matrix inverse *explicitly* rather
than calling the structure-keyed ``loss.inverse()`` — which, reading the
operand tree :math:`C - K_\mathrm{iso}` (a sum with an invertible leading
collision diagonal :math:`C`), would return the **iterative**
:class:`~orpheus.numerics.green_operator.GreenOperator` preconditioned
splitting.  For a 0-D loss operator that is a single small dense block the
iterative splitting is the wrong realization and an exact direct inverse is
right; encoding that decision as the *type* ``MatrixInverseOperator`` — not a
``strategy=`` flag on ``.inverse()`` — is the taxonomy §3 strategy-override
seam realized honestly (the type **is** the choice).

The ``(A, F)``-posed engine
:func:`~orpheus.numerics.eigenvalue.direct_eigenvalue` is the **dense
sibling** of this route: it forms the resolvent from a dense ``(A, F)`` pair
via :func:`numpy.linalg.solve` and extracts the dominant eigenpair with
:func:`~orpheus.numerics.eigenvalue.dominant_eigenpair`
(:func:`numpy.linalg.eig`, the largest-real eigenpair, a sign convention, and
a refusal of a complex dominant eigenvalue).  The homogeneous route shares
neither step with it, and `[M]` 2026-10-06 neither engine has a caller in
``orpheus/`` outside :mod:`orpheus.numerics` itself: the derivation oracle
``kinf_and_adjoint_spectrum_homogeneous``, formerly the only caller of
``dominant_eigenpair``, moved to the reference kernel
(:class:`~orpheus.derivations.common.dense_pencil.DensePencil`) when the
references were insulated from production numerics (2026-10-06).
Both remain the dense members of the three-engine family
(:ref:`three-eigenvalue-engines`), and ``direct_eigenvalue`` is the
cross-engine oracle of
``tests/gates/homogeneous/test_homogeneous.py::test_kinf_matches_direct_eigenvalue_engine_of_the_assembled_pair``:
the solver's rank-one :math:`\kinf` against the dense engine's eigenvalue of
the SAME assembled :math:`(\mathbf{A}, \mathbf{F})`, at ``rtol=1e-12``.  That
gate is stronger than its own docstring says.  The two sides still share
their input, the assembled pair, but no longer an eigen-solver: one reads a
trace off a solve, the other runs ``geev``, so the comparison is
independent on the derivation axis (``instrument-doctrine`` X4) and a rewire
regression in either reds it.

.. note::

   The labels :eq:`fixed-source-solve` and :eq:`keff-update` historically
   named the per-iteration fixed-source solve
   (:math:`\mathbf{A}\boldsymbol{\phi}^{(n)} = \mathbf{Q}_f^{(n)}`) and
   the production/absorption eigenvalue ratio of the retired power
   iteration.  They are retained on the **direct** analogues: the
   loss solve against the emission spectrum,
   :math:`\mathbf{u} = \mathbf{A}^{-1}\boldsymbol{\chi}` (the single solve
   that replaces the per-iteration sequence), and the eigenvalue
   :math:`\kinf = \langle\nu\Sigma_f, \mathbf{u}\rangle` (the converged limit
   the iteration ratio approached).  The classical
   production/absorption form
   :math:`k = (\nu\Sigma_f\cdot\phi)/(\Sigma_a\cdot\phi)` remains a valid
   one-group balance identity — it is the per-group balance the
   :class:`~data.macro_xs.mixture.Mixture.absorption_xs` property reports
   alongside :math:`\kinf` — but it is not the computational path.

.. note::

   Because :math:`\mathbf{A}` is a single small dense block, the whole
   solve is one :func:`scipy.linalg.lu_factor` of :math:`\mathbf{A}`, one
   ``lu_solve`` backsolve against :math:`\boldsymbol{\chi}`, and one
   :math:`G`-term inner product.  There is **no inner iteration, no outer
   iteration and no eigen-solver** — the homogeneous solver is the one
   deterministic solver in ORPHEUS with no iteration at all.  This is what
   makes it the instantaneous reference eigenvalue for every other solver on a
   homogeneous problem.

.. _homogeneous-exact-reference:

The exact reference, and the bound the solve is held to
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**What the solve is compared with.**  Every float64 is a dyadic rational, so
the loss matrix and the dyad built from the FLOAT cross sections are exact
rational matrices, and their eigenproblem has an exact rational answer.
:mod:`orpheus.derivations.common.exact_homogeneous` computes it with
:class:`fractions.Fraction`: :math:`\mathbf{A}` assembled exactly as
:math:`\operatorname{diag}(\Sigma_t) - \Sigma_{s0}^{\mathsf T} -
2\Sigma_{2,0}^{\mathsf T}`, :math:`\mathbf{A}^{-1}` by Gauss–Jordan elimination,
:math:`\mathbf{u} = \mathbf{A}^{-1}\boldsymbol{\chi}`,
:math:`\kinf = \langle\nu\Sigma_f, \mathbf{u}\rangle`, the flux in the
production gauge :math:`\boldsymbol{\phi} = 100\,\mathbf{u}/\kinf`, and the
one-group condensed cross sections
:math:`\bar\sigma_x = \langle\Sigma_x, \boldsymbol{\phi}\rangle /
\langle 1, \boldsymbol{\phi}\rangle`.  No rounding occurs anywhere.  This is
"the exact answer for the float inputs", which differs from the exact
answer for the decimal cross sections a library tabulates by a few ULP (the
four-group table under "Analytical Solutions" above shows both); the former
is the only thing a floating-point solver can be held to.

**Why it is independent of the solver it checks.**  The reference uses the
same rank-one identity as the solver, so agreement alone could be two
spellings of one formula (``instrument-doctrine`` X4).  It therefore
certifies itself against the DEFINING equations, in exact arithmetic, before
any comparison with production:
:math:`\mathbf{A}\mathbf{A}^{-1} = \mathbf{I}`,
:math:`\mathbf{F}\boldsymbol{\phi} = \kinf\,\mathbf{A}\boldsymbol{\phi}` with
zero residual, :math:`\operatorname{tr}(\mathbf{A}^{-1}\mathbf{F}) = \kinf`,
the gauge, and the cone.  It reaches no production solver, operator or
LAPACK routine; it shares only the float inputs with production, which is
the point.  ``test_the_exact_reference_certifies_itself`` runs the
certificate on every case.

**The bound, derived rather than chosen.**  The gate
``tests/gates/homogeneous/test_kinf_exact_reference.py`` asserts that
production's :math:`\kinf`, flux and condensed cross sections lie within a
forward-error bound of the exact values, evaluated per case in exact
rationals from the solve's own LU factors.  With :math:`u = 2^{-53}` the unit
roundoff and :math:`\gamma_m = mu/(1 - mu)` (Higham, *Accuracy and Stability
of Numerical Algorithms*, 2nd ed., SIAM 2002, §3.1):

1. **Assembly.**  :math:`\mathbf{E} = \hat{\mathbf{A}} - \mathbf{A}`, the
   rounding of the float assembly, is computed exactly.
2. **Solve.**  The computed :math:`\hat{\mathbf{u}}` satisfies
   :math:`(\hat{\mathbf{A}} + \boldsymbol{\Delta})\hat{\mathbf{u}} =
   \boldsymbol{\chi}` with :math:`\lvert\boldsymbol{\Delta}\rvert \le
   \gamma_{3G}\,\mathbf{P}^{\mathsf T}\lvert\hat{\mathbf{L}}\rvert
   \lvert\hat{\mathbf{U}}\rvert` (Higham Theorem 9.4), read from the computed
   factors of the same :math:`\hat{\mathbf{A}}` in the same process, so the
   bound holds for whatever LAPACK the host has.
3. **Forward error.**  With :math:`\mathbf{W} = \lvert\mathbf{A}^{-1}\rvert
   (\lvert\mathbf{E}\rvert + \gamma_{3G}\mathbf{P}^{\mathsf T}
   \lvert\hat{\mathbf{L}}\rvert\lvert\hat{\mathbf{U}}\rvert)` and
   :math:`r = \lVert\mathbf{W}\rVert_\infty < 1` (asserted),

   .. math::

      \lvert\hat{\mathbf{u}} - \mathbf{u}\rvert
      \;\le\; (\mathbf{I} - \mathbf{W})^{-1}\mathbf{W}\lvert\mathbf{u}\rvert
      \;\le\; \mathbf{W}\lvert\mathbf{u}\rvert
      + \tfrac{r}{1-r}\,\bigl\lVert\mathbf{W}\lvert\mathbf{u}\rvert\bigr\rVert_\infty\,\mathbf{1},

   the componentwise form of Higham Theorem 7.4 with no first-order
   truncation (the Neumann series of a non-negative matrix of norm
   :math:`r` majorises the inverse).
4. **Inner product.**  A :math:`G`-term dot product in any summation order is
   within :math:`\gamma_G\lvert\mathbf{a}\rvert^{\mathsf T}\lvert\mathbf{b}\rvert`
   of the exact one (Higham Eq. 3.5); the gate uses :math:`\gamma_{G+1}`,
   admitting one extra rounding per term for the pairing's weight multiply.
5. **Gauge and condensation.**  Two roundings in
   :math:`\hat{\mathbf{u}}\cdot\mathrm{fl}(100/\hat n)`, and, because every
   term is non-negative (the cone), relative errors that add in the ratio
   :math:`\bar\sigma_x`, plus one rounding for its division.

The derivation is written out in full in the gate's docstring.  Evaluated
per case (`[M]` 2026-10-01; ULP of the exact value):

.. list-table:: The derived bound, per case
   :header-rows: 1
   :widths: 30 12 30 14 14

   * - case
     - :math:`\kinf`
     - flux, per group
     - :math:`\bar\sigma_{\rm prod}`
     - :math:`\bar\sigma_{\rm abs}`
   * - ``homo_1eg``, ``mixture_A_1g``
     - 3.75
     - 5.2
     - 18.8
     - 12.5
   * - ``homo_2eg``, ``homo_2eg_with_eg``, ``mixture_A_2g``
     - 18.4
     - 24.0, 34.4
     - 71.1
     - 75.2
   * - ``homo_2eg_n2n``
     - 17.7
     - 23.7, 30.5
     - 45.1
     - 76.9
   * - ``homo_4eg``, ``mixture_A_4g``
     - 43.2
     - 47.3, 79.5, 52.4, 69.8
     - 173.4
     - 115.1

`[M]` 2026-10-01 the rank-one route's :math:`\kinf` sits at most
:math:`0.86` ULP from the exact value over the eight cases (``homo_1eg`` and
``mixture_A_1g`` exactly, the 2-group cases :math:`-0.86`, ``homo_2eg_n2n``
:math:`+0.73`, the 4-group cases :math:`-0.45`), and every flux component
within 7 % of its bound; the gate's 40 rows pass.

**What the bound can and cannot see.**  A worst-case bound holds for every
rounding pattern of an LU solve, so it is loose by construction; its value
is that no platform can move it.  The structural mutations a rewire would
introduce red it by orders of magnitude (`[M]` 2026-10-01, the
test-architect's battery over the 40 rows: a dropped dot-product term reds
6, the :math:`(n,2n)` term dropped from :math:`\mathbf{A}` reds 3 — the one
case that carries it — a non-transposed scattering matrix reds 18, a gauge
target off by :math:`2^{-40}` reds 8).  A 1-ULP change of one :math:`\chi`
entry does **not** red it, and cannot: it moves the exact :math:`\kinf` by
about 1 ULP, inside every case's bound, and the rounding of the LU solve
alone is entitled to more.  The resolution in :math:`\chi` is 8 ULP of the
entry at one group, 16 at two, between 32 and 64 at four.  Nothing on the
homogeneous route is pinned at the bit, deliberately: its one platform
primitive is the LU solve, and bit-level platform independence is a claim
this route does not make.


.. _spectral-invisibility:

Spectral invisibility: what the eigenvalue gate cannot see
----------------------------------------------------------

Two natural mistakes in the operator spelling — swapping the factor **order**
of the resolvent (:math:`\mathbf{F}\mathbf{A}^{-1}` instead of
:math:`\mathbf{A}^{-1}\mathbf{F}`), and **transposing** the materialized
resolvent — are *spectrally invisible*: they move :math:`\kinf` by **exactly
zero**.  Every :math:`k`-level gate (the cross-engine equivalence, the
registry's analytical anchor) is therefore structurally **blind** to them.  The
reason is a pair of standard linear-algebra identities, and understanding them
dictates *which* gate must catch the bug.

.. note::

   **Scope.**  The production solve forms no resolvent
   (:ref:`homogeneous-rank-one-route`); the identities below are about the
   composition ``MatrixInverseOperator(loss) @ production`` as an operator
   of the algebra, which the object gate named below still constructs and
   which any consumer materializing :math:`\mathbf{A}^{-1}\mathbf{F}` meets.
   On the rank-one route neither mistake has a single-step analogue that
   leaves :math:`\kinf` unmoved `[R]`: transposing :math:`\mathbf{A}` alone
   gives :math:`\langle \mathbf{A}^{-1}\nu\Sigma_f, \boldsymbol{\chi}\rangle`,
   which differs from :math:`\kinf` unless :math:`\mathbf{A}` is symmetric,
   and so does exchanging the roles of :math:`\boldsymbol{\chi}` and
   :math:`\nu\Sigma_f`.  The two together are the **adjoint** problem,
   whose eigenvalue is :math:`\kinf` exactly and whose eigenvector,
   :math:`\mathbf{A}^{-\mathsf T}\nu\Sigma_f`, is the importance rather than
   the flux; the flux rows of the exact-reference gate
   (:ref:`homogeneous-exact-reference`) are what see that pair, and a
   non-transposed scattering matrix reds 18 of its 40 rows `[M]`.

**Factor-order swap is a similarity transform.**  For any invertible
:math:`\mathbf{A}`,

.. math::
   :label: resolvent-similarity

   \mathbf{A}\,\bigl(\mathbf{A}^{-1}\mathbf{F}\bigr)\,\mathbf{A}^{-1}
   \;=\; \mathbf{F}\mathbf{A}^{-1},

.. vv-status: resolvent-similarity documented
.. Structural linear-algebra identity (similarity of the swapped resolvent),
.. NOT a solver claim; explains why the factor-order/transpose mutations are
.. spectrally invisible. Verifiable content is the object-level matrix gate
.. ``test_K_operator_as_matrix_is_the_resolvent`` (rtol=1e-12) named below.

so :math:`\mathbf{F}\mathbf{A}^{-1} = \mathbf{A}\,\mathbf{M}\,\mathbf{A}^{-1}`
is **similar** to :math:`\mathbf{M} = \mathbf{A}^{-1}\mathbf{F}` with
similarity matrix :math:`\mathbf{A}`.  Similar matrices share their entire
spectrum, so :math:`\lambda_{\max}(\mathbf{F}\mathbf{A}^{-1}) =
\lambda_{\max}(\mathbf{A}^{-1}\mathbf{F}) = \kinf` **identically**.  The
eigenvector is *not* invariant — if
:math:`\mathbf{M}\boldsymbol{\phi} = \kinf\boldsymbol{\phi}` then
:math:`(\mathbf{F}\mathbf{A}^{-1})(\mathbf{A}\boldsymbol{\phi}) =
\kinf(\mathbf{A}\boldsymbol{\phi})`, i.e. the mode maps
:math:`\boldsymbol{\phi} \mapsto \mathbf{A}\boldsymbol{\phi}` — but the
*eigenvalue* is untouched.

**Transpose preserves the spectrum.**  A matrix and its transpose share a
characteristic polynomial,
:math:`\det(\mathbf{M}^{\mathsf T} - \lambda\mathbf{I}) =
\det\bigl((\mathbf{M} - \lambda\mathbf{I})^{\mathsf T}\bigr) =
\det(\mathbf{M} - \lambda\mathbf{I})`, so
:math:`\lambda_{\max}(\mathbf{M}^{\mathsf T}) = \lambda_{\max}(\mathbf{M}) =
\kinf` (the eigenvectors of :math:`\mathbf{M}^{\mathsf T}` being the *left*
eigenvectors of :math:`\mathbf{M}`).

Both mutations were verified to give :math:`|\Delta\kinf| = 0.0` exactly.  A
value gate that reads only the eigenvalue — or even the sign-normalised
eigenvector, whose *shape* is fixed only up to the
:math:`\boldsymbol{\phi} \mapsto \mathbf{A}\boldsymbol{\phi}` remapping above —
cannot see either bug.  The committed catcher is therefore an **object-level**
gate, not a spectral one: ``test_K_operator_as_matrix_is_the_resolvent``
asserts the *materialized matrix itself* equals the reference resolvent,

.. math::
   :label: resolvent-object-gate

   [\mathbf{K}] \;=\; \texttt{np.linalg.solve}(\mathbf{A},\,\mathbf{F})
   \qquad (\text{rtol} = 10^{-12}),


.. implements:: resolvent-object-gate
   :by: orpheus.numerics.eigenvalue.direct_eigenvalue

   **Implemented by** 4 sites. Every symbol that executes this
   equation's arithmetic is declared, not only the canonical one: a
   test is adjudicated against the transcription it actually ran, so
   declaring a single site would refute the tests that exercise the
   others.

.. implements:: resolvent-object-gate
   :by: orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator

.. implements:: resolvent-object-gate
   :by: orpheus.numerics.operator.LinearOperator.as_matrix

.. implements:: resolvent-object-gate
   :by: orpheus.numerics.operator.OperatorProduct

and both mutations move :math:`[\mathbf{K}]` by :math:`O(1)`
(:math:`\mathbf{F}\mathbf{A}^{-1} \neq \mathbf{A}^{-1}\mathbf{F}` unless
:math:`\mathbf{A}` and :math:`\mathbf{F}` commute;
:math:`\mathbf{M}^{\mathsf T} \neq \mathbf{M}` unless :math:`\mathbf{M}` is
symmetric).  The general lesson is **pin the object, not just its spectrum**: a
value gate can be blind to an entire mutation class for structural (here
spectral-similarity) reasons, so a resolvent-forming operator earns an
intrinsic gate on the matrix it produces, above and beyond the eigenvalue it
feeds.


.. _three-eigenvalue-engines:

Why a direct engine: the three eigenvalue realisations
------------------------------------------------------

The dominant eigenpair of :math:`\mathbf{A}^{-1}\mathbf{F}` is what *every*
deterministic ORPHEUS solver ultimately wants; what differs is the
**realisation**.  :mod:`orpheus.numerics.eigenvalue` ships three siblings of
the same generalised eigenproblem
:math:`\mathbf{A}\boldsymbol{\phi} = \tfrac{1}{k}\mathbf{F}\boldsymbol{\phi}`:

.. list-table:: The three eigenvalue engines (:mod:`orpheus.numerics.eigenvalue`)
   :header-rows: 1
   :widths: 26 22 52

   * - Engine
     - Convergence
     - When it is the right realisation
   * - :func:`~orpheus.numerics.eigenvalue.power_iteration`
     - iterative, **linear** (rate :math:`|k_1/k_0|`)
     - Large, **sweep-only** loss operators that are never densely formed
       (SN, CP, MoC, diffusion).  These only *apply* :math:`\mathbf{A}^{-1}`
       — a :term:`sweep` or Krylov inner solve — and drive :math:`k` up the dominant
       mode through the
       :class:`~orpheus.numerics.eigenvalue.EigenvalueSolver` Protocol, which
       sees only a normalised-source fixed point, never a dense matrix.
   * - :func:`~orpheus.numerics.eigenvalue.direct_eigenvalue`
     - **exact** (one LAPACK shot)
     - **Small, densifiable** operators posed as a dense
       :math:`(\mathbf{A}, \mathbf{F})` pair — few-group / few-region
       problems.  Forms the dense resolvent
       :math:`\mathbf{A}^{-1}\mathbf{F}` via :func:`numpy.linalg.solve` and
       delegates the extraction to
       :func:`~orpheus.numerics.eigenvalue.dominant_eigenpair`.  The direct
       (non-iterative) sibling of ``power_iteration``.  The 0-D homogeneous
       medium does not use it: its rank-one fission makes the eigenpair one
       solve and one inner product (:ref:`direct-eigensolve-solve`).
   * - :func:`~orpheus.numerics.eigenvalue.rayleigh_quotient_iteration`
     - iterative, **superlinear** (locally quadratic)
     - Polishing an eigenpair *estimate* to the eigenpair NEAREST its
       Rayleigh quotient — **not** necessarily the dominant one (warm-start
       near the mode you want).  The bordered / augmented-Newton form, in
       which the previous iterate enters as the normalisation **row**.  Not
       yet wired into a meshed solver — that integration, and its use as the
       adjoint-:math:`\phi^*` vehicle, is
       `#277 <https://github.com/deOliveira-R/ORPHEUS/issues/277>`_.

The infinite-medium :math:`\kinf` takes **none of the three**: it is the
closed form :math:`\langle\nu\Sigma_f, \mathbf{A}^{-1}\boldsymbol{\chi}\rangle`
(:ref:`homogeneous-rank-one-route`), one LU solve and one inner product.  An
iterative engine would only approximate, at a convergence tolerance, an answer
the closed form gives to the rounding of one solve; and a dense eigen-solve
would add an eigen-solver whose last bit is the platform library's.  The
closed form exists because the fission operator has rank one: the resolvent
:math:`\mathbf{A}^{-1}\mathbf{F} = (\mathbf{A}^{-1}\boldsymbol{\chi})\,
(\nu\Sigma_f)^{\mathsf T}` has a single non-zero eigenvalue and :math:`G-1`
exact zeros, so its dominance ratio is :math:`0` and even
:func:`~orpheus.numerics.eigenvalue.power_iteration` would converge in one
step from any start with a component along :math:`\mathbf{A}^{-1}\boldsymbol{\chi}`.

.. note::

   **The rank-one structure is a hypothesis, and the type enforces it.**
   The route is exact only while :math:`\mathbf{F}` is one column times one
   row.  That holds for every mixture the tree can pose: the mixture carries
   one fission spectrum, and the fission kernel the hub builds refuses
   factors that are not one-dimensional, so a multi-spectrum production term
   is not spellable on this path rather than silently mishandled.  The
   physics does not promise it: the flat-flux production-weighted
   :math:`\chi` is a data-layer approximation, and the exact operator over
   :math:`K` fissile isotopes has rank up to :math:`K`
   (`#549 <https://github.com/deOliveira-R/ORPHEUS/issues/549>`_).

**The pure-math verification of the engines.**  All three engines — and the
shared :func:`~orpheus.numerics.eigenvalue.dominant_eigenpair` extraction the
dense one delegates to — are verified against a **transport-unrelated,
hand-derived closed-form eigenproblem** — :math:`\mathbf{M} =
V\operatorname{diag}(\lambda)V^{-1}` with chosen eigenpairs, and the rank-1
closed form :math:`k = v^{\mathsf T} A^{-1} u` — in the pure-math gate
``tests/gates/numerics/test_eigenvalue.py`` (the closed-form eigenproblem, the direct
``dominant_eigenpair`` surface with its one-home relocation proofs, and the RQI
gates).  This is a **closed-form** reference: V&V pillar 1, the *only* pillar
that proves an eigenvalue (MMS is source-driven and cannot).  In that file no
reference VALUE is produced by :func:`numpy.linalg.eig`: `[M]` 2026-10-01 its
four calls to it are a precondition that the host's LAPACK still returns a
negative-sum eigenvector on the sign-convention fixture, the bodies of two
stand-ins that neuter ``dominant_eigenpair`` to prove its guards have one
home, and an executable copy of the engine's body that exists only to show
each named mutation reddens its gate.

Once an engine is pinned against domain-independent ground truth it is
**trusted machinery**, and two routes may both call it without contamination,
provided they assemble its input differently.  That is the relationship
between :func:`~orpheus.numerics.eigenvalue.direct_eigenvalue` and the
registry's oracle
:func:`~orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous`:
both end in :func:`numpy.linalg.eig`, one on the pair the hub assembles
through the transport operators (:ref:`direct-eigensolve-assembly`), the other
on :math:`\mathbf{A} = \operatorname{diag}(\Sigma_t) - (\Sigma_{s} +
2\Sigma_2)^{\mathsf T}` assembled by the fused route.  The production solver
calls no eigen-solver, so its comparisons with both are independent on the
derivation axis as well.  The registry's analytical values are therefore
float64 LAPACK eigenvalues, accurate to a few ULP and platform-dependent in
their last bits; the structurally independent value the solver is held to at
the bit scale is the exact rational reference
(:ref:`homogeneous-exact-reference`).


.. _homogeneous-rates-and-normalisation:

Reaction rates, flux normalisation, and the one-group condensation
-------------------------------------------------------------------

The eigenvector :math:`\boldsymbol{\phi}` is determined only up to a
scalar multiple.  After the solve,
:func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite` normalises
the flux so that the **fission** production rate is 100 n/cm\ :sup:`3`/s:

.. math::
   :label: normalisation

   \boldsymbol{\phi} \leftarrow \boldsymbol{\phi} \times
   \frac{100}{\nu\boldsymbol{\Sigma}_\mathrm{f} \cdot \boldsymbol{\phi}}


.. implements:: normalisation
   :by: orpheus.homogeneous.solver.solve_homogeneous_infinite

   **Implemented by** 3 sites, each executing one part of the arithmetic:
   the solver, which builds the section and applies it; the hub's
   production-rate co-vector below, which computes the denominator
   :math:`\nu\Sigma_f\cdot\boldsymbol{\phi}` as a typed pairing on the pose;
   and — since step 3 U1, 2026-09-17 — the generic section
   :meth:`ScaleGauge.apply <orpheus.numerics.gauge.ScaleGauge.apply>`,
   which is where the division :math:`100 / (\nu\Sigma_f\cdot\phi)` now
   lives.  Every symbol that executes this equation's arithmetic is
   declared, not only the canonical one: a test is adjudicated against the
   transcription it actually ran, and
   ``tests/gates/homogeneous/test_homogeneous.py::test_post_solve_production_rate_is_100``
   runs all three.  (The solver's own body carried the ``phi * (100 / …)``
   update by hand until U1; the arithmetic is unchanged, and the
   exact-reference gate pins the gauged flux
   (:ref:`homogeneous-exact-reference`).  And until the CS4c coda, 2026-09-08, both of
   the then-two directives named the solver — a duplicate from the
   94-equation declaration pass whose "2 sites" body counted one symbol
   twice.)

.. implements:: normalisation
   :by: orpheus.homogeneous.solver.HomogeneousProblem.production_rate

.. implements:: normalisation
   :by: orpheus.numerics.gauge.ScaleGauge.apply

The normalisation denominator is the **fission** production rate
:math:`\nu\Sigma_f\cdot\boldsymbol{\phi}` only — consistent with the
production matrix :math:`\mathbf{F} = \boldsymbol{\chi}\otimes\nu\Sigma_f`.
The :math:`(n,2n)` neutrons are **not** in this denominator: they are a
loss-side transfer folded into :math:`\mathbf{A}` as
:math:`2\Sigma_2^T`, not a production channel (see
:ref:`scattering-matrix-convention` and the note under the production
matrix :eq:`fission-matrix` above).

The rescale is a RECORDED section, and the answer carries it
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Equation :eq:`normalisation` is a **choice**, and until step 3 U1
(2026-09-17) nothing recorded that it had been made.  The eigenvector is a
ray; picking the member with :math:`\nu\Sigma_f\cdot\varphi = 100` is
choosing a *section* of the :math:`(\mathbb{R}_+,\times)`-quotient, and a
consumer holding only :math:`\varphi` could not tell which section had
been applied — a live hazard, because ``[M]`` the tree applies **four
different functionals** under the one name "production rate" and this
page's target (:math:`100`) is the only one of the four that is not
:math:`1` (:ref:`the-gauge-section`).

So the line now names its own section:

.. code-block:: python

   # orpheus/homogeneous/solver.py — solve_homogeneous_infinite
   gauge = ScaleGauge(production_rate.evaluate, 100.0)   # the SECTION, as a value
   phi_column = gauge.apply(phi.reshape(ng, 1))          # φ · 100 / ⟨νΣf, φ⟩

and the answer carries it.
:class:`~orpheus.homogeneous.solver.HomogeneousResult` gained an
:class:`~orpheus.numerics.outcome.EigenOutcome` as its **first field** —
the k-eigen question over the 0-D pencil
(:attr:`~orpheus.homogeneous.solver.HomogeneousProblem.eigen_posing`), the
gauged representative, :math:`\lambda = k_\infty` with its one-point
trajectory, and the :class:`~orpheus.numerics.gauge.ScaleGauge` above.
Two consequences for readers of the result type:

* :attr:`~orpheus.homogeneous.solver.HomogeneousResult.k_inf` and
  :attr:`~orpheus.homogeneous.solver.HomogeneousResult.flux` are now
  **properties**, not stored fields: ``k_inf`` is ``outcome.keff`` and
  ``flux`` is a *view* of the outcome's state
  (``outcome.state[:, 0]`` — ``[M]`` ``np.shares_memory`` is ``True``, so
  there is one storage and the historical names survive unchanged).  No
  reader needed editing: ``[M]`` the type has exactly ONE construction
  site tree-wide (the solver itself) and nothing anywhere assigns to
  either name, so turning two fields into two properties is invisible
  outside the mint.
* the section is interrogable after the fact: ``result.outcome.gauge``
  reports the functional **as the object that ran** and the target it
  landed on, so ``gauge.functional(outcome.state)`` re-reads
  :math:`100` — ``[M]`` 2026-10-01 ``100.00000000000001`` at 2 groups and
  exactly ``100.0`` at 4 (one rounding in the division and :math:`G`
  in the re-read, so a last-place difference either way is the arithmetic,
  not a defect).

⚠ **The state is the posed** :math:`(n_g, 1)` **column, and the shape is
load-bearing.**  This is the gotcha the verification design surfaced, and
it is silent in the dangerous direction:
:class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
**accepts a flat** :math:`(n_g,)` **vector without complaint and returns a
different number.**  ``[M]`` 2026-10-01 on the fissile ``A`` family, with
the state already gauged so the column reads its target to the last place:

.. list-table:: The same functional on the same flux, two shapes
   :header-rows: 1
   :widths: 14 28 28 30

   * - groups
     - posed :math:`(n_g, 1)` column
     - flat :math:`(n_g,)` vector
     - the solve's :math:`k_\infty` (to reproduce the row)
   * - 1
     - ``100.0``
     - ``100.0``
     - ``1.5``
   * - 2
     - ``100.00000000000001``
     - ``200.0``
     - ``1.8750000000000004``
   * - 4
     - ``100.0``
     - ``411.2729251352303``
     - ``1.487761904761904``

**The mechanism is a silent broadcast, and it is worth spelling out
because the wrong number is not obviously wrong.**  The cross-section
field's values carry the pose's spatial axis, so they are shaped
:math:`(n_g, 1)`.  Pairing them with a :math:`(n_g, 1)` column is the
elementwise product the equation means.  Pairing them with a flat
:math:`(n_g,)` vector broadcasts instead — to the :math:`(n_g \times n_g)`
**outer product** — whose sum is
:math:`\bigl(\sum_g \nu\Sigma_{f,g}\bigr)\bigl(\sum_g \varphi_g\bigr)`
rather than :math:`\sum_g \nu\Sigma_{f,g}\,\varphi_g`.  ``[M]``
reproducing that product by hand recovers the flat column of the table to
the last digits at every group count.

Three readings follow, and each is a trap in its own right:

* **at 1 group the two coincide** — the outer product is :math:`1\times 1`.
  Every law written on this gauge is therefore stated at
  :math:`\geq 2` groups: a 1-group fixture is green under **either**
  shape and pins nothing.
* **the error is not a clean factor of** :math:`n_g`.  It is
  :math:`2.000` at two groups and :math:`4.113` at four, because the
  wrong expression is a product of two sums and the right one is a sum of
  products — so "the answer is off by the group count" is a rule that
  would hold on the 2-group fixture and break on the 4-group one.
* **the sibling primitive is strict.**  :meth:`EigenPosing.rayleigh
  <orpheus.numerics.posing.EigenPosing.rayleigh>` refuses the flat vector
  *loudly* on the same input.  The asymmetry — one numerics primitive
  raises, the reaction-rate functional broadcasts — is the thing to carry
  away; it is why the gauge's state is the posed column and why the law
  that checks it re-reads the functional rather than trusting the shape.

**Nothing about the arithmetic moved**, and that is checkable rather than
asserted.  The retired line read
``phi * (100.0 / production_rate.evaluate(phi.reshape(ng, 1)))``;
:meth:`ScaleGauge.apply <orpheus.numerics.gauge.ScaleGauge.apply>` is
``state * (target / functional(state))`` — the *same two operations in the
same order*, with the divisor computed on the same posed column, so
bit-identity is structural and not a coincidence of the fixture.  Only the
SHAPE the multiplication is written on changed (the flat vector became the
:math:`(n_g, 1)` view), and scaling by a scalar is elementwise.  ``[M]``
re-running the retired expression beside the shipped one on the fissile
``A`` family at 1, 2 and 4 groups: the flux is ``np.array_equal`` on
**3 of 3**, :math:`k_\infty` compares ``==``, and both condensed rates
compare ``==``.  The standing adjudicator of the gauged flux and both
rates is ``tests/gates/homogeneous/test_kinf_exact_reference.py``, which
holds each to the exact rational answer for its float inputs within a
derived bound (:ref:`homogeneous-exact-reference`); a gauge target off by
:math:`2^{-40}` reds its flux rows in all eight cases `[M]`.

⚠ One aliasing consequence, deliberate and gated: ``result.flux`` is a
**view** of ``result.outcome.state``, where before it was an independent
array.  Writing into it in place now writes into the outcome.  The gate
asserts the sharing (``np.shares_memory``) rather than tolerating it,
because the alternative — a defensive copy — would reintroduce the second
storage the outcome exists to remove.

**And one property is now a theorem the answer can restate.**  Because the
outcome carries the posing, the result can be asked for the *posing's own*
eigenvalue estimate — the Rayleigh quotient :math:`\langle
w, F\varphi\rangle / \langle w, A\varphi\rangle` of
:eq:`posing-balance-functional` — beside the eigenvalue the solve
returned.  ``[M]`` 2026-10-01 on this direct solve the two agree to a few
ulp of :math:`k`: :math:`4.4\times 10^{-16}` (2 ulp) against
:math:`k_\infty = 1.8750000000000004` and :math:`8.9\times 10^{-16}`
(4 ulp) against :math:`1.487761904761904` — there is no iteration here for
a convergence residual to hide in, so the gap is the rounding of the
quotient's two matrix-vector products and one solve.
This gap is the quantity
:class:`~orpheus.numerics.outcome.ExitReport`'s ``rayleigh_gap``
member exists to record — a diagnostic, deliberately never asserted at
construction (:ref:`the-solution-outcome` gives the two independent
reasons).  ``[M]`` no homogeneous result carries an exit report at U1: the
TYPE ships, the evaluator that mints one lands with the S\ :sub:`N`
entries at U2, so the number above is obtained by asking the outcome
directly (``result.outcome.rayleigh()``).

The rate is a typed integrated co-vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since campaign 1 CS4b (step S7, EE-1) that denominator — and the two
post-processing rates beside it — are not a hand-written contraction.
All three evaluations go through
:class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`,
the volume-integrated reaction rate

.. math::

   R_x(\varphi)
   \;=\; \int_V \sum_g \Sigma_{x,g}(\vec r)\,\varphi_g(\vec r)\,\mathrm{d}V
   \;=\; \sum_{\text{cells}} V_{\rm cell}\,
          \langle \Sigma_x, \varphi\rangle(\text{cell}),

which is the :math:`\varphi^\dagger = 1` **degenerate** of the
adjoint-weighted homogenization bilinear
:math:`\langle \varphi^\dagger, M[\Sigma_x]\varphi\rangle` (theorem
T1 of the algebra of record,
:mod:`orpheus.derivations.common.homogenization`; the weighted form is
live — pass ``adjoint=``). On the infinite-medium pose the spatial sum
has one term and :math:`V_{\rm cell}` is the quotient point's unit
weight, so the object degenerates to the bare group contraction — but
it is the *same* object the meshed solvers use, which is the point:
one functional, one contraction, no 0-D special case.

.. important::

   **The rate co-vectors read the POSE's measure — and since the CS4c
   coda that is true by construction, not by a rebinding step.** The
   functional's measure authority is its cross section's space (the
   :math:`\sigma`\ ↔geometry pairing tier), while the total-flux leg
   below is the *pose's* pairing. Bind the rates to a foreign space and
   the two legs would read **different** measures — so re-weighting the
   pose would move one and not the other, and the condensed cross
   sections would stop being ratios of commensurable quantities.

   The hub's cross-section fields are **born on the pose**
   (:ref:`direct-eigensolve-assembly`, step 2), so the two legs share a
   measure with nothing in between: the rate co-vectors
   :attr:`~orpheus.homogeneous.solver.HomogeneousProblem.production_rate`
   and
   :attr:`~orpheus.homogeneous.solver.HomogeneousProblem.absorption_rate`
   simply wrap fields that were never posed anywhere else.  The gate
   that adjudicates this is unchanged, and it is worth reading exactly,
   because its assertion is *not* "the rates double": the CS4a-R
   **G2.5** leg scales the pose's point weight by two and pins,
   bit-exactly (``×2`` and ``÷2`` are exact in binary floating point),
   that the normalised flux **halves** — proof that the pairing consults
   the space's measure at all — while both condensed cross sections stay
   **unchanged**, because :math:`\bar\sigma_x = \langle \Sigma_x,
   \varphi\rangle / \langle 1, \varphi\rangle` is the same pairing
   top and bottom and the weight cancels (the CS4a-R **XD-6**
   intensivity ruling).  Until CS4a-R that leg asserted the rates
   DOUBLE — the covariant reading the pre-review spelling shipped — and
   the ruling replaced it; a page that still says G2.5 "requires the
   rates to move" is quoting the pre-XD-6 gate.

   .. note::

      **What changed here at the coda (2026-09-08), and why the
      distinction matters.** Until then the fields were minted on a
      carrier's space and the solver *rebound* them —
      ``replace(field, space=space)`` — immediately before wrapping.
      That was correct, and this page recorded it as load-bearing:
      the rebinding existed precisely to satisfy G2.5. What replaced it
      is not a shorter spelling of the same step but the **removal of
      the step**:
      with no second space in scope there is no rebinding to forget, and
      the class of bug G2.5 catches becomes unspellable on this path
      rather than merely caught. `[M]` the change is bit-identical
      end-to-end — the byte gate of the day read 8 of 8 on
      :math:`k_\infty`, the flux bytes and both rates — because the
      pose and the retired carrier's mint were content-equal all along;
      that equality is what made the rebinding safe, and its removal is
      why nothing had to be re-baselined.

**The total flux is NOT a reaction rate.** The denominator of the
one-group condensation is

.. math::

   \langle 1, \varphi\rangle \;=\; \int_V \sum_g \varphi_g\,\mathrm{d}V ,

which carries no cross section: it is the pose's own **integration
co-vector**, so it stays
:meth:`FunctionSpace.inner_product
<orpheus.numerics.space.FunctionSpace.inner_product>` against a field
of ones rather than being dressed as an
:class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
with a unit cross section. `[M]` on the counting-weighted quotient
point it is bit-identical to ``float(phi.sum())``, which is what it
replaced.

The two condensed one-group cross sections
:attr:`~orpheus.homogeneous.solver.HomogeneousResult.sig_prod` and
:attr:`~orpheus.homogeneous.solver.HomogeneousResult.sig_abs` are then
**same-pairing ratios**

.. math::

   \bar\sigma_x \;=\;
   \frac{\langle \Sigma_x, \varphi\rangle}{\langle 1, \varphi\rangle},

and being ratios of two integrals against the *same* measure they are
**measure-invariant**: a future pose carrying a different point weight
moves numerator and denominator together and leaves
:math:`\bar\sigma_x` fixed. That is the CS4a-R **XD-6** ruling — a
quantity documented as a cross section must not scale with the point
weight — and its gate is the same ×2 pose mutation (G2.5), asserting
*the normalised flux halves, the ratios stay*: the halving is what
proves the pairing reads the measure at all (without it the gate would
be compatible with "the rate reads nothing"), and the two fixed ratios
are XD-6 itself. The counting weight of step 1 is a
convention the posing function states; it is not a contract these
ratios depend on.

Post-processing reads three energy-grid diagnostics off the mixture's
:class:`~orpheus.data.energy_grid.EnergyGrid` value object (campaign #276
P4-F — the group geometry lives on the grid, not re-derived in the solver)
and stores them on
:class:`~orpheus.homogeneous.solver.HomogeneousResult`:

- **Representative energy** — the plot abscissa:
  :attr:`~orpheus.homogeneous.solver.HomogeneousResult.representative_energy`
  :math:`= \bar E_g = \sqrt{E_g^{\mathrm{up}}\,E_g^{\mathrm{lo}}}`, the
  **geometric** group centre.
- **Flux per unit energy**: :math:`\phi_g / \Delta E_g`
  (:attr:`~orpheus.homogeneous.solver.HomogeneousResult.flux_per_energy`),
  with :math:`\Delta E_g = E_g^{\mathrm{up}} - E_g^{\mathrm{lo}}`.
- **Flux per unit lethargy**: :math:`\phi_g / \Delta u_g`
  (:attr:`~orpheus.homogeneous.solver.HomogeneousResult.flux_per_lethargy`),
  with :math:`\Delta u_g = \ln\!\bigl(E_g^{\mathrm{up}} / E_g^{\mathrm{lo}}\bigr)`.

Here :math:`E_g^{\mathrm{up}} =` ``edges[g]`` and :math:`E_g^{\mathrm{lo}} =`
``edges[g+1]`` are the upper / lower bounds of group :math:`g` under the
**fast-first descending** convention (group :math:`0` is the highest-energy
group, boundaries strictly decreasing; see
:ref:`canonical-group-convention`).

.. note::

   **Why the geometric centre (the P4-F correction).**  The spectrum is
   plotted as flux-per-lethargy on a **logarithmic** energy abscissa
   (``semilogx``).  The natural centre of a group on a log axis is the
   **geometric** mean :math:`\sqrt{E^{\mathrm{up}}E^{\mathrm{lo}}}`, which
   sits at the midpoint of the group's *lethargy* interval — exactly where a
   flux-per-lethargy value belongs.  The arithmetic midpoint
   :math:`\tfrac{1}{2}(E^{\mathrm{up}} + E^{\mathrm{lo}})` is biased
   **high** by the AM–GM inequality
   (:math:`\tfrac{1}{2}(a+b) \ge \sqrt{ab}`, the gap widening as a group
   spans more decades), so it plots each point to the right of the lethargy
   centre.  Before P4-F the result carried the (wrong) arithmetic midpoint;
   P4-F renamed the field ``eg_mid`` → ``representative_energy`` and moved it
   to the geometric centre.  The change is **purely the abscissa** — the
   flux *values* are unchanged, only the energy each group is plotted *at*
   moved.  At a thermal floor (:math:`E_g^{\mathrm{lo}} = 0` eV) the
   geometric mean degenerates;
   :attr:`~orpheus.data.energy_grid.EnergyGrid.representative_energy` falls
   back to half the upper edge there (still strictly inside the group), while
   the lethargy width :math:`\Delta u_g \to +\infty` is genuinely unbounded.

For synthetic verification mixtures with no physical energy grid
(:attr:`Mixture.eg` is ``None``) all three diagnostics are ``None`` and the
:attr:`~orpheus.homogeneous.solver.HomogeneousResult.flux_per_energy` /
:attr:`~orpheus.homogeneous.solver.HomogeneousResult.flux_per_lethargy`
properties raise — :math:`\kinf` and the flux spectrum are still
well-defined, only the per-energy plotting path is unavailable.


.. _example-problems:

Example Problems
=================

Aqueous Uranium Solution Reactor
----------------------------------

The simplest physical problem: water with dissolved uranium-235
(1000 ppm) at room temperature (294 K) and atmospheric pressure.  This
models a bare, infinite, aqueous homogeneous reactor — a configuration
historically important for early criticality experiments.

The mixture contains only three isotopes: :sup:`1`\ H, :sup:`16`\ O,
and :sup:`235`\ U.  Water provides the moderation (hydrogen
down-scatter) and :sup:`235`\ U the fission source.  The water density
is obtained from the IAPWS-IF97 steam tables.

See :func:`~data.macro_xs.recipes.aqueous_uranium`.

.. plot::
   :caption: Neutron spectrum for an aqueous uranium solution reactor
             (:math:`k_\infty \approx 1.036`).  The thermal Maxwellian
             peak near 0.025 eV, the 1/E slowing-down region, and the
             fast fission peak above 1 MeV are clearly visible.

   import numpy as np
   import matplotlib.pyplot as plt
   import warnings
   warnings.filterwarnings('ignore')

   from orpheus.data.macro_xs.recipes import aqueous_uranium
   from orpheus.homogeneous import solve_homogeneous_infinite

   mix = aqueous_uranium(temp_K=294, pressure_MPa=0.1, u_conc_ppm=1000.0)
   result = solve_homogeneous_infinite(mix)

   fig, ax = plt.subplots()
   ax.semilogx(result.representative_energy, result.flux_per_lethargy, 'b-', linewidth=1.2)
   ax.set_xlabel('Energy (eV)')
   ax.set_ylabel(r'Flux per unit lethargy $\phi / \Delta u$')
   ax.set_title(
       rf'Aqueous U Solution — $k_\infty$ = {result.k_inf:.5f}'
   )
   ax.set_xlim(1e-3, 1e7)
   ax.grid(True, alpha=0.3)
   plt.tight_layout()


PWR-Like Homogenised Cell
---------------------------

A more realistic problem: a PWR unit cell (UO\ :sub:`2` fuel, Zircaloy
cladding, borated water) **volume-homogenised** into a single mixture.
This is not a physically realisable configuration, but it exercises the
full cross-section pipeline with 12 isotopes, self-shielding of
:sup:`238`\ U resonances, and boron absorption.

The geometric homogenisation uses volume fractions from the pin-cell
geometry:

.. math::
   :label: pin-cell-volume-fractions

   f_\mathrm{fuel} = \frac{r_\mathrm{fuel}^2}{r_\mathrm{cell}^2}, \quad
   f_\mathrm{clad} = \frac{r_\mathrm{clad,out}^2 - r_\mathrm{clad,in}^2}
                           {r_\mathrm{cell}^2}, \quad
   f_\mathrm{cool} = \frac{r_\mathrm{cell}^2 - r_\mathrm{clad,out}^2}
                           {r_\mathrm{cell}^2}

.. vv-status: pin-cell-volume-fractions documented
.. Definitional geometric formula: the Wigner-Seitz pin-cell volume fractions
.. consumed by data.macro_xs.recipes.pwr_like_mix. A textbook geometry
.. definition, not a solver claim.

where :math:`r_\mathrm{cell} = p / \sqrt{\pi}` is the Wigner–Seitz
equivalent radius for a square lattice of pitch :math:`p`.

The mixture includes: :sup:`235`\ U, :sup:`238`\ U, :sup:`16`\ O (fuel),
five Zr isotopes (:sup:`90,91,92,94,96`\ Zr), :sup:`1`\ H,
:sup:`16`\ O (coolant), :sup:`10`\ B, :sup:`11`\ B.

See :func:`~data.macro_xs.recipes.pwr_like_mix`.

.. plot::
   :caption: Neutron spectrum for the PWR-like homogenised mixture
             (:math:`k_\infty \approx 1.014`).  Compared to the
             aqueous solution, the thermal peak is suppressed by boron
             absorption and the :sup:`238`\ U resonance self-shielding
             is visible in the epithermal range.

   import numpy as np
   import matplotlib.pyplot as plt
   import warnings
   warnings.filterwarnings('ignore')

   from orpheus.data.macro_xs.recipes import pwr_like_mix
   from orpheus.homogeneous import solve_homogeneous_infinite

   mix = pwr_like_mix()
   result = solve_homogeneous_infinite(mix)

   fig, ax = plt.subplots()
   ax.semilogx(result.representative_energy, result.flux_per_lethargy, 'r-', linewidth=1.2)
   ax.set_xlabel('Energy (eV)')
   ax.set_ylabel(r'Flux per unit lethargy $\phi / \Delta u$')
   ax.set_title(
       rf'PWR-Like Homogenised Cell — $k_\infty$ = {result.k_inf:.5f}'
   )
   ax.set_xlim(1e-3, 1e7)
   ax.grid(True, alpha=0.3)
   plt.tight_layout()


Comparison
-----------

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Property
     - Aqueous U Solution
     - PWR-Like Mixture
   * - :math:`\kinf`
     - 1.03596
     - 1.01357
   * - Fuel
     - Dissolved :sup:`235`\ U (1000 ppm)
     - UO\ :sub:`2` (3% enrichment)
   * - Moderator
     - Light water (294 K)
     - Borated water (600 K, 4000 ppm B)
   * - Isotopes
     - 3
     - 12
   * - Self-shielding
     - Negligible (:sup:`235`\ U dilute)
     - Significant (:sup:`238`\ U resonances)
   * - MATLAB reference
     - 1.03596
     - 1.01357


.. _infinite-medium-verification-pins:

Verification — what pins this chapter
=====================================

The homogeneous solver's verification evidence — the registry's analytical
:math:`\kinf` eigenvalues
(:mod:`orpheus.derivations.continuous.analytical.homogeneous`; every value,
the 1-group one included, is a float64 :func:`numpy.linalg.eig` eigenvalue of
the oracle's own fused assembly, which at one group is the closed form
:math:`\nu\Sigma_f/\Sigma_a` to one rounding), the
multi-group matrix-eigenvalue chain, the exact rational reference with its
derived bound (:ref:`homogeneous-exact-reference`), and the two 421-group
industrial cross-checks against the legacy MATLAB implementation — lives in
the verification part: :doc:`/theory/verification/homogeneous`. The same
:class:`~orpheus.derivations.common.verification_case.VerificationCase` objects
serve both that chapter's LaTeX equations and the test suite. The
auto-generated :doc:`/theory/verification/matrix` reports per-equation test
coverage; :ref:`theory-verification` carries the part-wide principles and
harness contracts.


Comparison with Spatially-Dependent Solvers
============================================

The homogeneous infinite-medium solver sits at the simplest end of the
solver hierarchy.  The following table compares it with the
spatially-dependent solvers available in ORPHEUS:

.. list-table::
   :header-rows: 1
   :widths: 20 20 20 20 20

   * - Aspect
     - Homogeneous
     - Collision Probability
     - Discrete Ordinates
     - Diffusion
   * - Spatial dependence
     - None
     - Region-averaged
     - Mesh-resolved
     - Mesh-resolved
   * - Angular dependence
     - None (isotropic)
     - Integrated out
     - Discrete ordinates
     - Fick's law
   * - Transport operator
     - :math:`\mathbf{A}^{-1}` (direct)
     - :math:`P_\infty` matrix
     - Diamond-difference sweep
     - Implicit solve
   * - Inner iterations
     - None
     - None
     - Scattering source
     - None
   * - Typical convergence
     - Direct (no iteration)
     - 10--20 outer
     - 20--50 outer
     - 100+ outer
   * - Eigenvalue computed
     - :math:`\kinf`
     - :math:`\kinf` (lattice)
     - :math:`\kinf` (lattice)
     - :math:`\keff` (core)
   * - Implementation
     - :func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite`
     - :class:`CPSolver`
     - :class:`SNSolver`
     - :class:`DiffusionSolver`

.. _homogeneous-development-history:

Development history
===================

Reverse-chronological changelog of the **architectural** milestones of
*this page's subject* — how the infinite-medium problem is posed, where
its consumed objects live, and what supplies its data.  It is
deliberately short: this chapter's machinery is mostly inherited, so the
space layer's own changelog is on :doc:`spaces`, the operator algebra's
on :doc:`operator_algebra`, and the S\ :sub:`N` solver's (which
introduced most of the shared operators this page reuses) is
:ref:`sn-development-history`.  Only milestones that changed *this
problem's* posing are recorded here; iteration-rate work, gate counts and
intermediate replans are deliberately omitted.  **Merge status is a git
question, never a frozen note** — every entry lands with its commit
hash, and ``git`` outranks this column.

.. list-table::
   :header-rows: 1
   :widths: 10 50 12 28

   * - When
     - Architectural milestone
     - Issue
     - Where
   * - 2026-10-02
     - **The infinite medium is defined as the point in phase space.**
       The user's ruling after the review of the reference
       specification's step 8: the infinite medium is posed on energy
       alone, the point in position and in direction, and holds a
       material and no geometry.  The reduction from Eq.
       :eq:`boltzmann` was rewritten to say which symmetry removes which
       factor of phase space (:ref:`infinite-medium-point-in-phase-space`);
       its earlier three conditions, the first of which read "Infinite
       geometry — no boundaries, so the flux is spatially uniform", are
       recorded in the "first got wrong" box there.  The
       direction-dependent problem is documented as the spatial marginal
       of a body problem, and the transport mean free path as the
       relaxation length of anisotropy
       (:ref:`diffusion-transport-xs-relaxation`).  No code changed.
     - `#405 <https://github.com/deOliveira-R/ORPHEUS/issues/405>`_
     - the ruling: ``a076f2c7``; this text: ``65b533fe``
   * - 2026-10-01
     - **The eigenpair is read off the rank-one fission; no eigen-solver
       remains on the path, and the byte pin retires for an exact
       reference.**  The trigger was a platform drift, not a defect:
       `[M]` after macOS 27.0.1 was installed on 2026-09-30 (Accelerate),
       ``geev`` returned :math:`\kinf` 1 ULP from its previous value on
       ``homo_4eg`` and ``mixture_A_4g`` with no ORPHEUS change, and the
       homogeneous byte pin went red.  The user's ruling the same day:
       compute :math:`\kinf = \langle\nu\Sigma_f, \mathbf{A}^{-1}\chi\rangle`
       from the rank-one dyad, and re-pose the pin as an exact rational
       reference with a ULP bound derived from the LU error bound, never a
       hand-picked number.

       **The route that retired.**  From taxonomy step 5b until this row
       the solver composed ``K = MatrixInverseOperator(loss) @ production``,
       materialized :math:`[\mathbf{K}]` column by column (:math:`G` LU
       backsolves) and took its eigenpair with
       :func:`~orpheus.numerics.eigenvalue.dominant_eigenpair`
       (:func:`numpy.linalg.eig`, largest real part, a
       :math:`\sum\phi \ge 0` sign convention, a refusal of a complex
       dominant).  Step 5b's re-baseline note justified the operator
       spelling by "the closed-form SymPy :math:`\kinf` of
       ``test_kinf_exact``"; that anchor was never SymPy — the registry
       values come from ``kinf_and_spectrum_homogeneous``
       (``numpy.linalg.solve`` then ``numpy.linalg.eig``), so it shared the
       eigen-solver's LAPACK with the route it anchored.  The page's
       four-group paragraph made the same false claim ("SymPy's symbolic
       eigenvalue solver").  Both are corrected in the body above.

       **The pin that retired.**  ``test_byte_stability.py`` and its
       fixture ``cs1_prewiring.json`` (captured at ``24a991ba``) compared
       :math:`\kinf` as hex, the flux as raw bytes and both condensed cross
       sections as hex for the eight producing mixtures.  A byte of a
       LAPACK output pins the platform; its successor,
       ``test_kinf_exact_reference.py``, holds the same quantities to the
       exact answer of their float inputs within a per-case derived bound
       (:ref:`homogeneous-exact-reference`).  The shared population list
       moved to ``tests/gates/homogeneous/_homogeneous_population.py``.
       `[M]` the rank-one route reads at most 0.86 ULP from the exact
       :math:`\kinf` over the eight cases, where ``geev`` read up to 1.55.

       **What rode along.**  :attr:`IsotropicFission.emission_spectrum
       <orpheus.transport.operators.isotropic_transfer.IsotropicFission.emission_spectrum>`
       names the dyad's column factor, so the solver and the operator's
       kernel read one array.  The Gauss-rule half of the same platform
       drift is recorded on :ref:`gauss-rules-correctly-rounded`.
     - `#549 <https://github.com/deOliveira-R/ORPHEUS/issues/549>`_
       (the rank one is the model's)
     - ``473a68de``, merged at ``a34f8c2c``
   * - 2026-09-08
     - **The infinite-medium problem gets its HUB, and the fabricated
       carrier retires** (campaign 1 residue, CS4c *coda*; rulings R-c1
       and R-c2, the user, 2026-09-08).  Two commits, in that order.

       **(1) The hub.**  Ruled verbatim (R-c1): *"The homogeneous
       problem needs a hub, just like the function SNMesh (future
       SNProblem) currently fulfills, to act as the place the consumed
       objects live (and a save state)."*
       :class:`~orpheus.homogeneous.solver.HomogeneousProblem` is that
       hub — a frozen dataclass over one
       :class:`~orpheus.data.macro_xs.mixture.Mixture`, every consumed
       object a per-instance ``cached_property`` minted from the mixture
       and from nothing else: the pose
       :attr:`~orpheus.homogeneous.solver.HomogeneousProblem.space` and
       the one-cell
       :attr:`~orpheus.homogeneous.solver.HomogeneousProblem.layout`; the
       kernel-tier material fields
       (:attr:`~orpheus.homogeneous.solver.HomogeneousProblem.scattering`,
       :attr:`~orpheus.homogeneous.solver.HomogeneousProblem.n2n`,
       :attr:`~orpheus.homogeneous.solver.HomogeneousProblem.fission`);
       the three cross-section fields **born on the pose**; the bound
       operators
       (:attr:`~orpheus.homogeneous.solver.HomogeneousProblem.collision`,
       the two isotropic transfers and their sum,
       :attr:`~orpheus.homogeneous.solver.HomogeneousProblem.loss`,
       :attr:`~orpheus.homogeneous.solver.HomogeneousProblem.production`,
       and — until 2026-09-12 — a ``multiplication`` resolvent, moved to
       the solver by the consumers campaign's R-cc5: a resolvent is a
       Strategy object, the hub owns the PENCIL);
       and the two typed rate co-vectors.
       :func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite`
       became a thin reader.  Three things retired in the same commit:
       the private ``_assemble_loss_operator`` helper (into
       :attr:`~orpheus.homogeneous.solver.HomogeneousProblem.loss`), the
       two ``replace(..., space=space)`` field re-poses (the fields are
       born on the pose, so there is nothing to re-pose), and the last
       production call to the fabricated carrier.  ⭐ *Per-instance* is
       the design decision, not an implementation detail: a module-scope
       memo would mask every decoy a gate installs on the pose, so the
       hub's identity claim is ``is``-stability **within** one instance,
       and ``==``-equality (with equal hashes) **across** hubs built over
       equal mixtures — ``is`` across hubs is deliberately not claimed.

       **(2) The retirement.**  With its last consumer gone,
       ``MaterialMesh.from_materials`` — the factory that fabricated a
       mesh-less one-cell carrier with ``[0, 1]`` edges, a node at
       ``0.5`` and a Cartesian chart that `[M]` nothing consumed — was
       deleted, and with it the three arms that existed only to serve it:
       ``areas``' "no faces at all" case, the :math:`S_N`-promotion
       refusal and the diffusion bounded-geometry refusal.  Each was
       **input-less**, and the argument is a closure rather than a
       census, which is what makes it durable: `[M]` the whole carrier
       hierarchy has exactly **one** producer of ``mesh = None`` —
       :meth:`SNProblem.from_axes
       <orpheus.sn.problem.SNProblem.from_axes>`, which
       synthesises a legacy adapter when ``len(axes) <= 2`` and returns
       ``None`` only above that — while every other constructor
       (``MaterialMesh.__init__``, ``SNProblem.__init__``,
       ``DiffusionMesh.__init__``) takes a mesh as a *required*
       argument.  So once the homogeneous path stopped building one, no
       producer of a mesh-less :math:`d \leq 2` carrier existed at all,
       and three refusals were left refusing nothing — which is
       ``plan-authoring`` §6c's mirror (a gate that ships with no case
       to catch, arrived at from the other direction: a case that
       retired out from under its gate).  ``mesh is None`` therefore
       carries ONE meaning again — the :math:`d \geq 3` axis-native
       carrier — and CS4b S7's ``ndim`` discrimination (G7.3) became the
       **singleton law** ``mesh is None ⟹ ndim ≥ 3``, gated with a
       positive control and every :math:`d \leq 2` construction path.
       A companion gate makes the retirement *unspellable* rather than
       merely done (``not hasattr(cls, "from_materials")`` over
       ``MaterialMesh`` / ``SNProblem`` / ``DiffusionMesh``); the homonym
       :meth:`EnergyAxis.from_materials
       <orpheus.numerics.axis.EnergyAxis.from_materials>` — a different
       object, and the one energy-arm rule both spellings of the pose
       route through — survives untouched.  Test fixtures that had used
       the fabricated carrier for convenience moved to ONE shared home,
       a **genuine** unit-width one-cell ``Mesh1D`` carrier, `[M]`
       like-for-like on every read surface the migrated sites touch.

       **(3) Why this is a re-source, not a re-baseline.**  `[M]`
       bit-identical across both commits: the byte gate 8 of 8 on
       :math:`k_\infty`, the flux bytes and both condensed cross
       sections, against a fixture captured before the campaign and never
       regenerated; and the operator tier (:math:`A` and :math:`F`) 8 of
       8 ``array_equal`` against a capture frozen on the pre-carve tree.
       The hub re-sources the *same* kernels, einsums, dyad and pairing —
       which is exactly what objective **O1** predicted, since the
       fabricated data was consumed by nothing.

       ⚠ **Interim home, stated.** The hub lives in the solver module.
       Carving it into a standalone ``Problem`` module with a thin
       ``Problem → Solution`` solver is the consumers campaign's work,
       alongside the same split for ``SNMesh`` → ``SNProblem``; the
       ruling names that as the long-term shape.
     - —
     - ``5caad3d6`` (the hub), ``39e7f32f`` (the retirement)
