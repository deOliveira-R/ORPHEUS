.. _theory-spectral-resolvent:

==========================================================================
Spectral resolvent — closed-form scalar Green's kernel (reserved)
==========================================================================

.. note::

   **Stub-grade theory page — implementation deferred.** This page
   reserves ``:label:`` anchors and the canonical-references
   bibliography for the spectral-resolvent reference solver. The
   Python folder exists at
   :mod:`orpheus.derivations.continuous.spectral_resolvent` with a
   README pinning the literature; the implementation is queued.

   When implementation lands, the archivist will expand this page with
   full derivations, V&V evidence tables, and structural-independence
   cross-checks against the sister methods.


Why this folder exists
======================

This is **Meaning (β)** in the
:ref:`reference-solvers-three-meanings` taxonomy: it constructs the
**same scalar Green's kernel** as :ref:`theory-trajectory-resolvent`
(Variant α), but via a **structurally independent integrand**:

- Variant α (``trajectory_resolvent/``): traces characteristic rays,
  integrates :math:`\tau` along each ray, sums multi-bounce via
  :math:`T = (I - S)^{-1}`. Integrand is exponential attenuation
  along characteristics.
- Spectral resolvent (this folder): closed-form spectral
  μ-integration of the within-medium angular Green's function.
  Integrand is :math:`T(\mu)\cdot\cosh(\rho\mu) \cdot \exp(-2R\mu)`
  closed-form per Sanchez 1986 Eq. (A6) / PS-1982 Eq. (21).

Both produce the same physical scalar kernel
:math:`G(\rho \to \rho')` under specialisation to the homogeneous
medium. The two-path agreement is treated as an **L1 cross-check
anchor** — agreement is meaningful evidence because the integrands
are structurally distinct (ray-traced vs spectral-μ).

Today ORPHEUS realises (β) only **indirectly** via the trajectory
construction. A direct closed-form evaluator of PS-1982 Eq. (21) is
the headline gap for hardening the (α)-vs-(β) cross-check.


Canonical references
====================

.. list-table::
   :header-rows: 1
   :widths: 32 18 50

   * - Reference
     - Local PDF
     - Provides
   * - **Sanchez, R. (1986).** *Transp. Theor. Stat. Phys.*
       vol. 15, no. 3, pp. 333–343.
       DOI: 10.1080/00411458608210456.
     - ✓
     - Eq. (A1)–(A7): spectral closed-form scalar kernel for sphere
       with linearly anisotropic scattering, specular + diffuse BC.
   * - **Pomraning, G.C., Siewert, C.E. (1982).**
       *J. Quant. Spectrosc. Radiat. Transfer* vol. 28, no. 6,
       pp. 503–506.
     - ✓
     - Eq. (21): explicit Fredholm integral form for sphere with
       isotropic scattering. The direct implementation roadmap.
   * - **Mitsis, G.J. (1963)** ANL-6787 PhD thesis.
     - ✓
     - Sphere = antisymmetric slab via parity reduction (Sec. IV.5).


.. _spectral-resolvent-implementation-roadmap:

Implementation roadmap
======================

Detailed planning lives at
``.claude/scratch/sanchez_chandrasekhar_gap.md`` (literature-researcher
memo). The headline tasks:

1. Implement PS-1982 Eq. (21) directly via
   :func:`mpmath.quad` over :math:`T(\mu)` and :math:`\cosh(\rho\mu)`
   integrands.
2. Extend to Sanchez 1986 linearly anisotropic via Eq. (A6) and the
   :math:`J_R(\rho')` reduced-current solver (Eq. 9–10).
3. Wire into the verification matrix with the (α)-vs-(β) cross-check
   harness.

Cross-references
================

- :ref:`reference-solvers-three-meanings` — taxonomy positioning.
- :ref:`theory-trajectory-resolvent` — sister method (α).
- :ref:`theory-singular-eigenfunction` — independent anchor (γ).
- ``orpheus.derivations.continuous.spectral_resolvent`` — Python
  package (currently empty; implementation pending).
