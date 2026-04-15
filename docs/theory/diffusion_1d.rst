.. _theory-diffusion-1d:

==================
1D Diffusion (P1)
==================

The diffusion equation is the lowest-order angular approximation of
the neutron transport equation: the P1 (first spherical-harmonic)
expansion truncated to the scalar flux. It is not a transport
solver in the strict sense — it discards angular information
beyond the current — but it is the workhorse of reactor design and
its verification is a mathematical problem in its own right.

This page carries the continuous reference solutions for the
ORPHEUS diffusion module and the equation labels the verification
tests point at. It is deliberately scoped to **1D plane geometry**:
multi-dimensional diffusion and cylindrical/spherical variants are
tracked as follow-ups.

.. contents::
   :local:
   :depth: 2


Key Facts
=========

- The 1D multigroup diffusion equation in plane geometry is a
  second-order elliptic boundary-value problem in :math:`\phi_g(x)`
  with vacuum Dirichlet conditions at the slab edges:

  .. math::
     :label: diffusion-operator

     -\frac{d}{dx}\!\left(D_g(x)\,\frac{d\phi_g}{dx}\right)
       + \Sigma_{r,g}(x)\,\phi_g
       = \sum_{g' \ne g} \Sigma_{s,g' \to g}\,\phi_{g'}
       + \frac{1}{k}\,\chi_g\sum_{g'}\nu\Sigma_{f,g'}\,\phi_{g'}.

- The diffusion coefficient in each region is computed from the
  transport cross section:

  .. math::
     :label: diffusion-coefficient

     D_g(x) \;=\; \frac{1}{3\,\Sigma_{\text{tr},g}(x)}.

- The removal cross section is the total minus in-group
  scattering:

  .. math::

     \Sigma_{r,g} \;=\; \Sigma_{t,g} - \Sigma_{s,g\to g}.

- Standard discretisation is central finite difference on a
  uniform mesh; the design spatial order is
  :math:`\mathcal O(h^{2})` for smooth cross sections and
  :math:`\mathcal O(h)` at material interfaces that do not lie on
  cell faces.

- Boundary conditions: ORPHEUS uses hard Dirichlet
  :math:`\phi_g(0) = \phi_g(L) = 0` for all groups — the zero-flux
  variant of the vacuum condition (contrast with the
  extrapolation-distance Marshak form). This choice is
  intentional: it lets the analytical reference solutions below be
  pure sinusoids without an extrapolation-length fudge, so the
  spatial convergence test isolates finite-difference error.


.. _diffusion-1g-bare-slab:

1-group bare slab
=================

The simplest verification configuration: a homogeneous slab of
thickness :math:`L` with zero-flux vacuum boundaries and a
single energy group. The diffusion equation collapses to

.. math::

   -D\,\phi''(x) + \Sigma_r\,\phi(x)
      = \frac{1}{k}\,\nu\Sigma_f\,\phi(x)

with :math:`\phi(0) = \phi(L) = 0`. Separation of variables
gives the eigenfunction

.. math::
   :label: bare-slab-eigenfunction

   \phi(x) \;=\; \sin\!\left(\frac{\pi x}{L}\right)

and the **geometric buckling**

.. math::
   :label: bare-slab-buckling

   B^{2} \;=\; \left(\frac{\pi}{L}\right)^{2}.

Substituting :eq:`bare-slab-eigenfunction` into the diffusion
equation yields the eigenvalue condition

.. math::
   :label: bare-slab-critical-equation

   D\,B^{2} + \Sigma_r \;=\; \frac{1}{k}\,\nu\Sigma_f,

which solves to

.. math::

   k \;=\; \frac{\nu\Sigma_f}{D\,B^{2} + \Sigma_r}.

Because the eigenfunction is independent of group in the
multigroup generalisation (all groups share the same spatial
:math:`\sin(\pi x/L)` shape), multigroup reduces to a
:math:`ng \times ng` matrix eigenvalue problem in the spectrum
vector — exactly what
:func:`orpheus.derivations.homogeneous.kinf_and_spectrum_homogeneous`
solves, plus an extra ``D B²`` removal term on the diagonal of
:math:`\mathbf{A}`.

This is a **T1 analytical reference**: no integration, no
quadrature, no iteration. See
:func:`orpheus.derivations.diffusion.derive_1rg_continuous` for the
Phase-0 :class:`~orpheus.derivations.ContinuousReferenceSolution`
that carries :math:`k_{\text{eff}}` and the continuous
multigroup eigenfunction callable.


.. _diffusion-2rg-fuel-reflector:

2-group fuel + reflector slab
=============================

A more demanding verification problem: fuel surrounded by a
reflector, both treated with 2-group diffusion, with vacuum
boundaries on the outer faces. The eigenfunction is no longer
a single sine — it is a linear combination of exponential
(or trigonometric) modes determined by the material properties
in each region, matched across the fuel/reflector interface.

Region system matrix
--------------------

Define :math:`\mathbf y(x) = [\phi_1, \phi_2, J_1, J_2]^{T}`
with currents :math:`J_g = -D_g\,\phi_g'`. In each homogeneous
region the multigroup diffusion equation is a first-order linear
ODE system

.. math::
   :label: diffusion-transfer-matrix

   \frac{d\mathbf y}{dx} \;=\; \mathbf S\,\mathbf y,
   \qquad
   \mathbf S \;=\; \begin{pmatrix} \mathbf 0 & -\mathbf D^{-1} \\
                                    -\mathbf M & \mathbf 0 \end{pmatrix},

where

.. math::

   \mathbf M \;=\; \text{diag}(\Sigma_{t,g})
                   - \Sigma_{s}^{T}
                   - \frac{1}{k}\,\chi\otimes(\nu\Sigma_f).

Propagation across a region of thickness :math:`t` is the matrix
exponential :math:`\mathbf T(t) = \exp(\mathbf S\,t)`: given the
state :math:`\mathbf y(x_0)` at the left edge, the state at the
right edge is :math:`\mathbf y(x_0+t) = \mathbf T(t)\,\mathbf y(x_0)`.
Continuity of :math:`\phi_g` and :math:`J_g` across an interface
is automatic when we carry :math:`\mathbf y = [\boldsymbol\phi;
\mathbf J]` and do not separately enforce the matching
conditions.

.. math::
   :label: diffusion-interface-matching

   \phi_g^{\text{fuel}}(x_i) \;=\; \phi_g^{\text{refl}}(x_i),
   \qquad
   -D_g^{\text{fuel}}\,\phi_g'^{\,\text{fuel}}(x_i)
     \;=\; -D_g^{\text{refl}}\,\phi_g'^{\,\text{refl}}(x_i).

The composed transfer matrix

.. math::

   \mathbf T_{\text{total}}(k) \;=\;
     \mathbf T_{\text{refl}}(H_{\text{refl}}; k)\,\cdot\,
     \mathbf T_{\text{fuel}}(H_{\text{fuel}}; k)

maps :math:`\mathbf y(0)` to :math:`\mathbf y(L)`, where
:math:`L = H_{\text{fuel}} + H_{\text{refl}}`.

Vacuum boundary conditions — transcendental eigenvalue
------------------------------------------------------

At :math:`x=0` the vacuum condition fixes
:math:`\phi_g(0) = 0` but leaves :math:`J_g(0)` free, so the
initial state is :math:`\mathbf y(0) = [\mathbf 0_{ng};\,\mathbf J_0]`
for unknown :math:`\mathbf J_0 \in \mathbb R^{ng}`. Propagating
through the full slab gives

.. math::

   \mathbf y(L) \;=\; \mathbf T_{\text{total}}(k)\,
                      \begin{pmatrix}\mathbf 0\\ \mathbf J_0\end{pmatrix}.

Imposing :math:`\boldsymbol\phi(L) = \mathbf 0` picks out the
upper-right :math:`ng\times ng` block of
:math:`\mathbf T_{\text{total}}`:

.. math::
   :label: diffusion-transcendental

   \mathbf T_{\text{total}}(k)_{[0{:}ng,\ ng{:}2ng]}\,\mathbf J_0
     \;=\; \mathbf 0.

A non-trivial :math:`\mathbf J_0` exists iff the block is
singular. The eigenvalue condition is therefore the
transcendental equation

.. math::

   \det\!\left(
      \mathbf T_{\text{total}}(k)_{[0{:}ng,\ ng{:}2ng]}
   \right) \;=\; 0,

solved by :func:`scipy.optimize.brentq` to machine precision
after bracketing by a coarse scan over :math:`k`. The analogous
setup for reflective boundaries (with :math:`\mathbf J(0) =
\mathbf J(L) = \mathbf 0`) is already used in
:mod:`orpheus.derivations.sn_heterogeneous` and gives the
lower-left block :math:`[ng{:}2ng,\ 0{:}ng]` instead.

This is a **T2 semi-analytical reference**: the transfer-matrix
composition uses :func:`scipy.linalg.expm` which is accurate to
double precision, and the eigenvalue is found to ~\ :math:`10^{-12}`
via brentq. No spatial mesh, no iteration-count knob.

Back-substitution for continuous :math:`\phi(x)`
-------------------------------------------------

Once :math:`k` is found, the null vector :math:`\mathbf J_0` is
extracted from the singular block :eq:`diffusion-transcendental`
via SVD. The continuous flux at any :math:`x` in the slab is then

.. math::
   :label: diffusion-back-substitution

   \mathbf y(x) \;=\; \mathbf T_{\text{left}}(x; k)\,
                      \begin{pmatrix}\mathbf 0\\ \mathbf J_0\end{pmatrix},

where :math:`\mathbf T_{\text{left}}(x; k)` is the composition of
region transfer matrices from :math:`0` up to :math:`x` (crossing
any interfaces along the way, and using a partial-region transfer
matrix for the final stretch). Extracting the first :math:`ng`
components gives :math:`\boldsymbol\phi(x)`.

The back-substituted :math:`\phi_g(x)` is **mesh-independent**:
the test chooses its own cell centres, calls
:meth:`~orpheus.derivations.ContinuousReferenceSolution.phi_on_mesh`,
and compares the diffusion solver's output to the continuous
reference at exactly those points. See
:func:`orpheus.derivations.diffusion.derive_2rg_continuous`.


Verification
============

- **Bare slab L1 eigenvalue** — the Phase-0 continuous reference
  ``dif_slab_2eg_1rg`` pulls :math:`k` from the analytical matrix
  eigenvalue and :math:`\phi_g(x)` from :eq:`bare-slab-eigenfunction`;
  the diffusion solver must reproduce both to better than
  :math:`10^{-10}` for the eigenvalue and :math:`\mathcal O(h^{2})`
  for the flux shape.
- **Fuel+reflector L1/L2** — the continuous reference
  ``dif_slab_2eg_2rg`` produces :math:`k` from
  :eq:`diffusion-transcendental` and :math:`\phi_g(x)` from
  :eq:`diffusion-back-substitution`. Comparison against the
  solver at successive mesh refinements gives the measured order
  of the finite-difference discretisation. This replaces the
  Richardson-extrapolated reference that previously served this
  role (see the verification-campaign audit in
  :doc:`/verification/reference_solutions`).

Both cases live under ``operator_form="diffusion"`` in the
Phase-0 registry.


References
==========

- Bell, G. I. and Glasstone, S., *Nuclear Reactor Theory*,
  Van Nostrand Reinhold, 1970. Chapter 7 covers multigroup
  diffusion theory; §7.4 specifically treats the slab
  eigenvalue problem and the two-region interface matching.
- Stacey, W. M., *Nuclear Reactor Physics*, 3rd ed., Wiley,
  2018. Ch. 3 on diffusion theory and Ch. 8 on multigroup
  formulation.
- Duderstadt, J. J. and Hamilton, L. J., *Nuclear Reactor
  Analysis*, Wiley, 1976. Ch. 5 (one-group) and Ch. 7
  (multigroup) — the transfer matrix formulation is
  spelled out explicitly.
