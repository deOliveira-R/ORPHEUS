r"""Carlvik-Galerkin spectral solver for monoenergetic anisotropic transport.

Mounted on Dahl & Sjostrand (1979) NSE 69, 114-125 [DahlSjostrand1979]_:
**Legendre-Galerkin spectral expansion of Carlvik's integral equation**
for one-speed transport in a homogeneous slab or sphere with linearly
anisotropic scattering and vacuum boundary conditions.

This is the third structurally-independent verification pillar in this
project (Pillar 2 — semi-analytical):

* :mod:`orpheus.derivations.continuous.fn_method` is also Pillar 2
  but uses **boundary collocation in Case singular eigenfunctions**.
* :mod:`orpheus.derivations.continuous.singular_eigenfunction` is
  Pillar 2 via **Mitsis / Westfall-Metcalf integral-equation
  Fredholm iteration**.
* This package is Pillar 2 via **Galerkin spectral expansion in the
  spatial variable** (even-Legendre for slab, odd-Legendre for the
  reduced flux :math:`r\phi(r)` in sphere) of Carlvik's integral
  form.

Three completely different mathematical pillars solving the same
physical problems → high-confidence redundant cross-check.

What the method delivers
------------------------

The defining feature: the linearization Eq. (4) of Dahl-Sjostrand
gives the **full eigenvalue spectrum** (typically 6-11 eigenvalues
per case, including complex-conjugate pairs at high anisotropy and
high modes), not just the fundamental. This is unique among the
verification pillars in the project — F_N gives only the dominant
mode by construction.

Structurally independent of every F_N result in
:mod:`...fn_method` — the Galerkin spectral form of Carlvik's
integral equation is a different mathematical object than Case's
singular-eigenfunction expansion projected onto F_N's collocation
points.

Literature
----------

* Carlvik 1968, *Nucl. Sci. Eng.* **31**, 295. The original derivation
  of the recurrences for :math:`A_{m,n}` and :math:`B_{m,n}`.
  Dahl-Sjostrand explicitly flag a **typographical sign error in
  Carlvik's Eq. (4b)** — Dahl-Sjostrand's recurrences are the
  corrected master.
* Dahl-Sjostrand 1979, *Nucl. Sci. Eng.* **69**, 114-125. The matrix
  eigenvalue formulation Eq. (3) and the block-matrix linearization
  Eq. (4) that turns the c-search into a single standard eigenproblem.
* Sood/Forster/Parsons 1999, LANL LA-13511. The benchmark test set
  whose ``*-1-1-SL/SP`` cases (linearly anisotropic slab and sphere,
  P_1) this package targets.

Layout
------

* :mod:`.origins` — Branch-1 SymPy: closed-form derivations of the
  Galerkin algebra, low-order :math:`A_{m,n}` and :math:`B_{m,n}`
  matrix elements, and the block-matrix linearization. Foundation
  tests at ``tests/derivations/test_galerkin_spectral_symbolic.py``.
* :mod:`.core` — shared production primitives: numerical evaluation
  of the Carlvik recurrences (:mod:`.core.carlvik_recurrences`) and
  assembly of the Eq.(3) / Eq.(4) matrices
  (:mod:`.core.galerkin_matrix`).
* :mod:`.slab` — Branch-2 production: :func:`.slab.solve_galerkin_spectral_slab`
  for bare-critical multiplying slabs with linearly anisotropic
  scattering. Reproduces Dahl-Sjostrand Table II.
* :mod:`.sphere` — Branch-2 production: :func:`.sphere.solve_galerkin_spectral_sphere`
  for bare-critical spheres. Reproduces Dahl-Sjostrand Table I.

References
----------

.. [DahlSjostrand1979] Dahl, E. B. & Sjostrand, N. G. (1979).
   "Eigenvalue spectrum of multiplying slabs and spheres for
   monoenergetic neutrons with anisotropic scattering."
   *Nuclear Science and Engineering* **69**, 114-125.
   DOI: 10.13182/NSE69-114.
"""
from __future__ import annotations
