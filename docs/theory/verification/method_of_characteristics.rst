.. _theory-verification-method-of-characteristics:

=========================
Method of Characteristics
=========================

.. Seeded at task #10 stage V3: the per-solver sections of the
   dissolved ``docs/theory/verification.rst`` plus the Verification
   and MMS Verification chapters of
   :doc:`/theory/methods/method_of_characteristics` (case
   descriptions moved here; the method chapter keeps a thin pointer
   and its ERR-019 investigation history).  Prose carried verbatim;
   V5–V7 author the part-level framing.

Flat-source analytical eigenvalue
=================================

The MOC method solves the transport equation along characteristic rays.  The
characteristic ODE for :term:`angular flux` along a ray with direction :math:`\hat\Omega`:

.. math::

   \frac{d\psi}{ds} + \Sigma_t \psi = \frac{Q}{4\pi}

The segment-average solution:

.. math::

   \bar\psi = \psi_{\rm in} \frac{1 - e^{-\Sigma_t \ell}}{\Sigma_t \ell}
   + \frac{Q}{4\pi\Sigma_t} \left(1 - \frac{1 - e^{-\Sigma_t \ell}}{\Sigma_t \ell}\right)

For a homogeneous medium with isotropic incoming flux
(:math:`\psi_{\rm in} = Q/(4\pi\Sigma_t)`), this simplifies to
:math:`\bar\psi = Q/(4\pi\Sigma_t)` regardless of segment length — the
flat-source solution is exact.  The eigenvalue is then :math:`k = \nu\Sigma_f/\Sigma_a`.

.. include:: ../../_generated/moc_derivation.rst


Test suite
==========

102 Tests Across Four Levels
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 10 20 15 35

   * - Level
     - File(s)
     - Tests
     - Coverage
   * - L0
     - ``tests/moc/test_quadrature.py``
     - 24
     - Weight sums, TY values, shapes, validation
   * - L0
     - ``tests/moc/test_ray_tracing.py``
     - 20
     - Ray-circle intersection, region ID, segments, volume, links
   * - L0
     - ``tests/moc/test_verification.py``
     - 27
     - Single-track attenuation, equilibrium flux, fission-only, (n,2n),
       scatter isolation, geometric invariants, protocol compliance, volume
       tracking, boundary conditions
   * - L1
     - ``tests/moc/test_moc.py``
     - 6
     - Homogeneous eigenvalue (1G/2G/4G), heterogeneous (slow)
   * - L1
     - ``tests/moc/test_properties.py``
     - 4
     - Particle balance, positivity, flux consistency, thermal depression
   * - L1
     - ``tests/moc/test_verification.py``
     - 13
     - Eigenvalue + flux ratio, heterogeneous + monotonicity, particle
       balance (2G), flux positivity (all groups), material sensitivity
   * - L2
     - ``tests/moc/test_verification.py``
     - 4
     - Ray spacing convergence, azimuthal convergence, polar convergence
   * - XV
     - ``tests/moc/test_verification.py``
     - 2
     - MOC vs CP cross-verification (slow)

Run the full suite (excluding slow tests)::

   pytest tests/moc/ -v -k "not slow"


Homogeneous Infinite Medium
---------------------------

For homogeneous geometry with reflective BCs, the flux is spatially flat
and :math:`k_{\text{eff}} = \lambda_{\max}(A^{-1}F)`.

.. list-table::
   :header-rows: 1
   :widths: 10 14 20 20

   * - Groups
     - :math:`k_\infty`
     - Error (8 azi, TY-3)
     - Tolerance
   * - 1
     - 1.5000
     - :math:`< 10^{-15}` (exact)
     - :math:`< 10^{-4}`
   * - 2
     - 1.8750
     - :math:`\sim 10^{-6}`
     - :math:`< 10^{-4}`
   * - 4
     - 1.4878
     - :math:`\sim 10^{-6}`
     - :math:`< 10^{-4}`

The 1-group case is exact because :math:`k = \nu\Sigma_f / \Sigma_a`
is independent of the angular flux distribution.  The multi-group cases
converge through the power iteration.


Why Homogeneous Verification Is Necessary But Insufficient
----------------------------------------------------------

For a homogeneous medium, the flat-source approximation is **exact**:
the source is genuinely spatially flat.  Consequences:

1. :math:`\Delta\psi = 0` everywhere (all boundary fluxes equal
   :math:`Q / \Sigt{}`).
2. Angular integration weights cancel in the eigenvalue ratio.
3. The transport sweep contributes nothing --- :math:`\phi = 4\pi Q / \Sigt{}`.

This means homogeneous tests **cannot detect**:

- Wrong angular integration weights (ERR-019)
- Wrong scattering convention (SigS vs SigS.T)
- Wrong boundary condition linking
- Wrong ray spacing / volume conservation

**Rule:** Every transport solver must be verified on heterogeneous
multi-group problems.  Homogeneous success gives false confidence.


Heterogeneous Cross-Verification
--------------------------------

For a 2-region cylindrical pin cell (fuel + coolant), the MOC
eigenvalue is compared against the CP solver on the same geometry:

.. list-table::
   :header-rows: 1
   :widths: 15 20 20 20

   * - Groups
     - MOC :math:`k_{\text{eff}}`
     - CP :math:`k_{\text{eff}}`
     - Gap
   * - 2
     - 0.6078
     - 0.6072
     - 0.001
   * - 4
     - 0.5399
     - 0.5394
     - 0.001

The ~0.1% gap is consistent with the white-BC (CP) vs reflective-BC
(MOC) approximation difference.  This is verified by
``tests/moc/test_verification.py::TestXVCrossVerification``.


Convergence Properties
----------------------

All convergence studies use a 2-region 2G pin cell (fuel + coolant,
:math:`r_{\text{fuel}} = 0.5` cm, pitch = 2.0 cm) with materials A
(fissile) and B (scatterer) from the verification library.

**Ray spacing convergence** (:math:`N_\varphi = 32`, TY-3):

.. list-table::
   :header-rows: 1
   :widths: 15 15 15 15 15

   * - :math:`t_s` (cm)
     - Tracks
     - :math:`k_{\text{eff}}`
     - :math:`|\Delta k|`
     - Ratio
   * - 0.16
     - 528
     - 0.608256
     -
     -
   * - 0.08
     - 1036
     - 0.608021
     - :math:`2.36 \times 10^{-4}`
     -
   * - 0.04
     - 2048
     - 0.608129
     - :math:`1.09 \times 10^{-4}`
     - 2.2
   * - 0.02
     - 4092
     - 0.608038
     - :math:`9.10 \times 10^{-5}`
     - 1.2
   * - 0.01
     - 8168
     - 0.608029
     - :math:`9.23 \times 10^{-6}`
     - 9.9

The eigenvalue converges with decreasing ray spacing.  The convergence
is not perfectly monotonic (non-monotone behaviour between 0.08 and 0.04
reflects the discrete nature of track placement on the annular boundary),
but the refinement from 0.02 to 0.01 shows a clear order-of-magnitude
improvement.

**Azimuthal convergence** (:math:`t_s = 0.02` cm, TY-3):

.. list-table::
   :header-rows: 1
   :widths: 15 15 15 15

   * - :math:`N_\varphi`
     - Tracks
     - :math:`k_{\text{eff}}`
     - :math:`|\Delta k|`
   * - 4
     - 524
     - 0.605797
     -
   * - 8
     - 1028
     - 0.607311
     - :math:`1.51 \times 10^{-3}`
   * - 16
     - 2048
     - 0.607959
     - :math:`6.47 \times 10^{-4}`
   * - 32
     - 4092
     - 0.608038
     - :math:`7.96 \times 10^{-5}`
   * - 64
     - 8188
     - 0.608057
     - :math:`1.89 \times 10^{-5}`

The azimuthal convergence shows a clean reduction factor of ~2--8×
per doubling of :math:`N_\varphi`.  The error at :math:`N_\varphi = 64`
is dominated by the ray spacing, not the angular discretisation.

**Polar convergence** (:math:`N_\varphi = 16`, :math:`t_s = 0.03` cm):

.. list-table::
   :header-rows: 1
   :widths: 15 15 15 20

   * - :math:`N_p`
     - :math:`k_{\text{eff}}`
     - :math:`|\Delta k|` from TY-1
     - Quadrature points
   * - TY-1
     - 0.607538
     -
     - 1 per half-space
   * - TY-2
     - 0.607632
     - :math:`9.4 \times 10^{-5}`
     - 2 per half-space
   * - TY-3
     - 0.607790
     - :math:`2.5 \times 10^{-4}`
     - 3 per half-space

The polar convergence is weak: TY-1 to TY-3 gives only a
:math:`\sim 2.5 \times 10^{-4}` shift.  This is because the TY
quadrature is specifically optimised for the Bickley-function integrals
that arise in MOC, so even TY-1 (a single polar angle) gives a
surprisingly good approximation.  **Recommendation:** TY-3 is the
default; increasing to TY-2 or TY-1 is acceptable for quick
exploratory runs.


.. _moc-mms-verification:

MMS Verification
================

The Method of Manufactured Solutions (MMS) provides an L1-level
verification of the flat-source MOC spatial operator.  Unlike the
SN MMS (which injects a per-ordinate source into
:func:`~orpheus.sn.solver.solve_sn_fixed_source`), the MOC MMS requires
**per-segment, per-angle manufactured sources** along each
characteristic because the streaming residual is angle-dependent.

Ansatz
------

A smooth radial scalar flux centred at :math:`(P/2, P/2)`:

.. math::
   :label: moc-mms-psi-ref

   \phi_{\text{ref}}(r)
   = 1 + A\,\cos\!\bigl(\pi\,r^{2}/R^{2}\bigr),
   \qquad r = \sqrt{(x - P/2)^{2} + (y - P/2)^{2}},
   \quad R = P/2

The :math:`r^{2}` argument ensures :math:`C^{\infty}` smoothness at
:math:`r = 0`.  With :math:`A = 0.3` the ansatz is positive
everywhere (minimum :math:`1 - A = 0.7`).  The radial form is
deliberate: any product-of-cosines Cartesian ansatz centred at the
cell centre integrates to **exactly zero** over any concentric
annulus, making the FSR averages degenerate.

Manufactured source along each characteristic
---------------------------------------------

With isotropic angular flux
:math:`\psi_{\text{ref}} = \phi_{\text{ref}}/(4\pi)`, substituting
into the characteristic ODE
:math:`\mathrm{d}\psi/\mathrm{d}s + \Sigma_t\psi = Q` gives:

.. math::
   :label: moc-mms-qext

   Q_{\text{mms}}(x,y,\varphi_a,\theta_p) =
     \frac{1}{4\pi}\bigl[
       \sin\theta_p\,(\cos\varphi_a\,\partial_x\phi_{\text{ref}}
                     + \sin\varphi_a\,\partial_y\phi_{\text{ref}})
       + \Sigma_t\,\phi_{\text{ref}}
     \bigr]

The streaming term depends on track direction :math:`(\varphi_a, \theta_p)`.
When averaged over :math:`4\pi` it vanishes (both
:math:`\sum\omega_a\cos\varphi_a = 0` by quadrature symmetry, and
forward + backward sweeps cancel), so **an isotropic per-FSR source
cannot carry this information**.  The MMS sweep injects the source
per segment at the segment midpoint.

Why isotropic source fails
--------------------------

The isotropic average of :math:`Q_{\text{mms}}` over all angles is
:math:`\Sigma_t\,\phi_{\text{ref}}/(4\pi)`.  If fed as a per-FSR
source into a standard MOC fixed-source solve, the resulting scalar
flux differs from :math:`\phi_{\text{ref}}` by a transport correction
:math:`\nabla\cdot\mathbf{J} / \Sigma_t` that does not converge to
zero with FSR refinement.

Scalar flux reconstruction
--------------------------

Each segment has its own :math:`\tilde q_{\text{seg}}`, so the
standard Boyd Eq. 45 equilibrium term
:math:`4\pi Q_i / \Sigma_t` cannot be used directly (it assumes
:math:`Q` is constant per FSR).  The equilibrium value is instead
computed analytically:

.. math::
   :label: moc-mms-reference-equilibrium

   \langle\phi_{\text{ref}}\rangle_i
   = 1 + \frac{A\,R^{2}}{\pi\,(r_2^{2} - r_1^{2})}
     \bigl[\sin(\pi\,r_2^{2}/R^{2}) - \sin(\pi\,r_1^{2}/R^{2})\bigr]

and the transport correction :math:`\delta\phi_i` from the sweep
is added:

.. math::
   :label: moc-scalar-flux-reconstruction

   \phi_i = \langle\phi_{\text{ref}}\rangle_i
          + \frac{\delta\phi_i}{A_i\,\Sigma_t}

The :math:`\delta\phi` is the only term carrying flat-source spatial
discretisation error; it converges at :math:`\mathcal{O}(h^{2})` or
better as the FSR mesh is refined with fixed track spacing and angular
quadrature.

Convergence evidence
--------------------

.. list-table::
   :header-rows: 1

   * - :math:`N_{\text{ann}}`
     - :math:`\|e\|_{L^{2}}` (inner FSRs)
     - Measured order
   * - 4
     - :math:`5.8 \times 10^{-2}`
     - —
   * - 16
     - :math:`7.4 \times 10^{-3}`
     - 3.0
   * - 64
     - :math:`1.1 \times 10^{-3}`
     - 2.7

The measured order exceeds :math:`\mathcal{O}(h^{2})` because the
midpoint-rule source integral benefits from the symmetry of the
exponential attenuation kernel (even-order error terms cancel).

Implementation
--------------

- :mod:`orpheus.derivations.continuous.mms.moc` — MMS case, standalone sweep,
  continuous-reference registration
- ``tests/moc/test_mms.py`` — L1 convergence consumer test
- The outermost FSR (square-border region) is excluded from the
  convergence measurement because its complex geometry has a fixed
  track-sampling error that does not converge with FSR refinement.
