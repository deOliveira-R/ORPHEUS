.. _path-integral:

=========================================================
The Transport Path Integral: One Object, Five Methods
=========================================================

.. contents:: Contents
   :local:
   :depth: 1


.. Machine header — the ``nexus-meta`` schema for this page.  Ingestion is
.. PENDING nexus#1 Phase 2 (the directive is not yet registered), so the
.. schema is rendered here as a collapsed dropdown and machine-consumed
.. later.  The header is PROVISIONAL: it declares the page's structural
.. role now; the content-heavy fields are completed at Phase H.

.. dropdown:: Machine header — ``nexus-meta`` schema (PROVISIONAL · role · thesis · axes)
   :color: muted

   .. code-block:: yaml

      module: transport
      concept: path_integral
      status: "SCAFFOLD — labelled anchors + synopses only; full authoring at Phase H (corpus plan §3.6)"
      role: "root of the transport corpus; parent of methods/index — the frame every method derives FROM"
      thesis: >
        the five transport methods are five discretizations of ONE object
        (the sum over neutron histories); the reaction operators C, S, F are
        method-invariant AT FIXED MULTIGROUP DATA; a method chooses how to
        realize the propagator (L+C)^-1 — with diffusion the one exception
        (a limit of the object, not a quadrature of it)
      taxonomy_axes:
        A1: "how (L+C)^-1 is realized   — sweep-DAG / track / region-pair kernel / sampled / not-at-all (a limit)"
        A2: "where S is resummed        — outer Neumann (SI) / direct inverse / exact spectral"
        A3: "angular representation     — ordinates / harmonics / continuous / Case ν-spectrum"
      shared_operators: [MultiplicationOperator, IsotropicScattering, IsotropicN2N, FissionOperator]   # orpheus.transport.operators — the SAME objects across methods
      eigenvalue_posing: "k and α are properties of the OPERATOR, posed before any discretization; every method inherits the posing"
      gated_on: ["#298", "#299", "Phase-I literature survey"]
      authored_at_phase_H: [name_chain, generator_splitting_table, three_axis_method_table, method_placement_map, eigenvalue_spectral_yield]
      depends_on: [operator_algebra, discretization]
      parent_of: [methods/index]


.. note:: **This page is a labelled scaffold, not the authored chapter.**

   Its section anchors exist so that every method chapter and the
   :doc:`transport-methods entry </theory/methods/index>` can cross-link
   the root **now**; each section below carries a one-paragraph synopsis of
   *what it will establish* and a marker for the deferred authoring. The
   full derivation is scheduled for **Phase H** of the documentation-corpus
   restructure — gated on two upstream algebra corrections (tracked as
   ``#298`` and ``#299``) and a literature-survey pass, per corpus plan
   ``§3.6``.

This is the **root of the transport corpus** — the page every method
derives *from*, and the parent of the :doc:`transport-methods entry
</theory/methods/index>`. Where that entry states the differential
transport equation the deterministic methods discretize, this page answers
the prior question: *what is the one object that all of the methods
approximate, which part of it is shared, which part varies, on which axes
do methods differ, and where does each one land?*

The thesis, in one line: **the transport methods are not five different
subjects but five discretizations of one object — the sum over neutron
histories.** The operator algebra is powerful precisely because the
*reaction* operators are the same objects in every method — a shared-code
fact, not an analogy:
:class:`~orpheus.transport.operators.MultiplicationOperator`,
:class:`~orpheus.transport.operators.IsotropicScattering` and the fission
operator are the *same Python objects* imported by S\ :sub:`N`, diffusion
and the infinite-medium solver from
:mod:`orpheus.transport.operators`. A method is a choice of how to realize
what remains.


.. _path-integral-one-object:

1. The one object
=================

Every ORPHEUS method computes the **first moment of one branching
stochastic process** — a neutron released into the medium, streaming
between collisions and branching at them, with the :term:`scalar flux` the
expected track-length density summed over all histories. Linearity is what
lets the branching walk collapse to a single weighted path (the
many-to-one / spine decomposition), so the *same* object admits three
readings: a **collision-order series** (the Peierls collision-number
expansion — the name ORPHEUS already uses), a **resolvent** of the
transport generator, and a **Monte Carlo expectation** over sampled
histories. Fission does not break the path reading; it makes each history's
multiplicative weight exceed one.

*Authored at Phase H.*


.. _path-integral-invariant:

2. What is invariant — the reaction operators
=============================================

Collision (:math:`C`), scattering (:math:`S`) and fission (:math:`F`) are
**method-invariant reaction operators at fixed multigroup data**. This is
demonstrated by shared code, not asserted by analogy — the same operator
objects are imported across S\ :sub:`N`, diffusion and the infinite-medium
solver. The scope condition is load-bearing and will be stated precisely:
once cross sections are condensed with a method's *own* flux, its :math:`S`
is no longer literally the same operator — the invariance is a statement
*at fixed data*, and condensation is a solution-weighted projection owned by
no single operator.

*Authored at Phase H.*


.. _path-integral-streaming:

3. What varies — realizing the propagator
=========================================

What a method actually chooses is **how to realize the propagator**
:math:`(L+C)^{-1}` — the inverse that carries a neutron from one collision
to the next. Here streaming carries its intuitive **Lagrangian** meaning:
leakage along the flight path between interactions. (The complementary
**Eulerian** reading — that :math:`\hat\Omega\cdot\nabla\psi` is the
divergence of the angular current — is taken up on the :doc:`diffusion page
</theory/methods/diffusion_1d>`, where it becomes the continuity law that
Fick's law closes.) The realization is *not* the only thing that changes
across methods: the split of the generator into :math:`L`, :math:`C` and
:math:`S` is itself method-dependent, and **diffusion is the one genuine
exception** — a limit of the object rather than a :term:`quadrature` of it.

*Authored at Phase H.*


.. _path-integral-generator-splitting:

4. The branch point — how the generator is split
================================================

The deepest branch between methods is a choice of **generator splitting**
:math:`\mathcal{A} = \mathcal{A}_0 + P`: which physics rides in the
propagated part and which becomes the source or jump kernel. Three
splittings organize the whole family — a **killing** split (the total cross
section rides in the attenuation functional; scattering becomes the source;
this is the collision-order series the deterministic methods sum), a
**jump** split (the stochastic process *is* the answer — analog Monte
Carlo), and a **majorized-jump** split (a virtual self-scatter flattens the
collision rate — Woodcock / delta-tracking). These are **different
splittings**, not equivalent readings of one process; the
deterministic-versus-stochastic distinction falls along this line.

*Authored at Phase H.*


.. _path-integral-axes:

5. The three orthogonal axes
============================

Once the object is fixed, a method is located by **three orthogonal
choices**, not a single dichotomy:

- **A1 — how** :math:`(L+C)^{-1}` **is realized:** a :term:`sweep` over a cell
  dependency DAG (S\ :sub:`N`), exact exponential attenuation along tracks
  (method of characteristics), a region-to-region kernel (collision
  probability), sampled histories (Monte Carlo), or **not realized at all**
  (diffusion — a limit).
- **A2 — where** :math:`S` **is resummed:** an outer Neumann iteration
  (source iteration), a direct inverse, or an exact spectral resummation.
- **A3 — the angular representation:** :term:`discrete ordinates <ordinate>`, spherical
  harmonics, a continuous direction, or a Case discrete-plus-continuum
  spectrum. The angular representation is also what fixes the **trace** of
  the boundary law :math:`B`.

*Authored at Phase H.*


.. _path-integral-method-map:

6. Where each method lands
==========================

With the axes fixed, each method is a **point in their product space** —
and neighbours that a textbook keeps in separate chapters turn out
adjacent. S\ :sub:`N` and the method of characteristics sit on the **same
side**: a :term:`diamond-difference <diamond difference>` sweep is a rational (Padé) approximant of the
exponential that characteristics integrate exactly, so "negative flux in an
optically thick cell" is the pole of that approximant, not a scheme
pathology. The operator algebra :math:`A = L + C - S - B` is the form this
frame takes **on a deterministic angular–spatial grid**; it is *not* the
universal form — collision probability folds the boundary condition into
its kernel, and the Case / F\ :sub:`N` spectral methods carry no separate
:math:`C`, :math:`S`, :math:`F` at all.

*Authored at Phase H.*


.. _path-integral-eigenvalue:

7. Posing the eigenvalue problem
================================

The multiplication eigenvalue (:math:`k`) and the time eigenvalue
(:math:`\alpha`) are properties of the **operator**, posed *before any
discretization* — so every method inherits the same posing. This is also
where the path reading **yields to a spectral statement**: the
generation-by-generation history sum is a Perron–Frobenius / Krein–Rutman
statement about the mean-offspring operator :math:`A^{-1}F` (with
:math:`\alpha` the Malthusian growth rate), and in a multiplying medium the
naive sum over histories does **not** simply converge — the :math:`1/k`
rescaling that the eigenvalue posing supplies is exactly what makes the
generation series summable.

*Authored at Phase H.*
