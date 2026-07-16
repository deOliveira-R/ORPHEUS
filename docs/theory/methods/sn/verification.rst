Verification
============

.. note:: **Verification slice — automation pending.**  The per-page V&V
   table (equation label × test × level × ERR coverage, auto-filtered
   from the ``tests/_harness`` registry to this page's equation labels) is
   blocked on Nexus equation-label ↔ test linking; until it ships, the
   verification below is **hand-authored**.  The project-wide matrix lives
   at :doc:`/verification/matrix`.

.. _sn-mms-verification:

Method of Manufactured Solutions (1D slab)
-------------------------------------------

Homogeneous and heterogeneous eigenvalue tests verify :math:`\keff`
--- a scalar. They do not tell us whether the **spatial operator**
itself converges at the design order :math:`\mathcal O(h^{2})` of
diamond difference.  The Method of Manufactured Solutions closes
that gap by constructing a fixed-source problem whose exact angular
flux is known in closed form, so the error against the prescribed
flux is pure spatial-discretisation error.

**Ansatz.**  For a vacuum-BC slab of length :math:`L` in one energy
group, pick an isotropic angular flux

.. math::
   :label: sn-mms-psi

   \psi_n(x) = \frac{1}{W}\,A(x),
   \qquad A(x) = \sin\!\left(\frac{\pi x}{L}\right),

where :math:`W = \sum_n w_n = 2` for Gauss--Legendre.  Because
:math:`A(0) = A(L) = 0`, every ordinate vanishes at both faces ---
the vacuum boundary conditions are satisfied automatically, with no
inflow bookkeeping required on the caller side.  Since :math:`\psi_n`
is independent of ordinate, the scalar flux recovered by any
quadrature order is *exactly* :math:`\phi(x) = A(x)` --- the test
isolates spatial error from angular quadrature error.

**Manufactured source.**  Substituting :eq:`sn-mms-psi` into the
discrete ordinates transport equation :eq:`transport-cartesian`
(with the :math:`1/W` convention ORPHEUS uses),

.. math::

   \mu_n\,\frac{\partial\psi_n}{\partial x} + \Sigma_t\,\psi_n
   = \frac{1}{W}\!\left(\Sigma_s\,\phi + Q^{\text{ext}}_n\right),

and solving algebraically for :math:`Q^{\text{ext}}_n` gives

.. math::
   :label: sn-mms-qext

   Q^{\text{ext}}_n(x)
   = \mu_n\,A'(x) + \bigl(\Sigma_t - \Sigma_s\bigr)\,A(x)
   = \mu_n\,\frac{\pi}{L}\cos\!\left(\frac{\pi x}{L}\right)
     + \bigl(\Sigma_t - \Sigma_s\bigr)\sin\!\left(\frac{\pi x}{L}\right).

The :math:`W` factor cancels cleanly because the ansatz was already
divided by :math:`W`, so what we hand the solver is the full residual
without any additional rescaling.  The expression is per-ordinate and
linear in :math:`\mu_n`: a constant isotropic external source *cannot*
drive a non-trivial manufactured flux because the streaming term
:math:`\mu_n\,\psi'_n` is odd in :math:`\mu`.  That is the fundamental
reason MMS for SN requires the :math:`Q_{\rm aniso}` plumbing path ---
no "cheat" with a cell-by-cell isotropic source exists.

**Why :math:`\sin(\pi x/L)`?**  The ansatz is smooth
(:math:`C^{\infty}`) so all derivatives of the exact solution exist
and DD's :math:`\mathcal O(h^{2})` truncation error dominates.  It
vanishes at both boundaries for free.  Its derivatives do not collapse
to a polynomial --- a cubic ansatz, for instance, has a constant
second derivative so the DD truncation term :math:`\psi'''` would be
zero and the error could disappear for a non-physical reason,
hiding bugs.  Trigonometric or exponential ansätze have bounded
but non-zero derivatives of every order and therefore expose the
leading truncation term cleanly.

**Implementation.**  The case is built by
:func:`orpheus.derivations.continuous.mms.sn.build_1d_slab_mms_case` and
consumed by :func:`orpheus.sn.solve_sn_fixed_source`.  The latter
accepts a per-ordinate external source of shape
:math:`(N, n_g, n_x, n_y)` (Issue #196 PR-INDEX-5 principled layout,
the ``g`` axis directly after ``N``) and threads it through the sweep's
:math:`Q_{\rm aniso}` slot --- merging additively with any P1+
scattering contribution the solver itself builds.  This bare-array form
is the **bulk-only / vacuum** special case of the composite source
``q = q_bulk ⊕ q_∂`` the solver also accepts (see
:ref:`sn-composite-fixed-source`); this isotropic slab MMS is
vacuum-automatic, so the boundary leaf is identically zero.  Vacuum
boundary conditions are applied via the mesh-level BC infrastructure
described in :ref:`boundary-conditions`:
:func:`solve_sn_fixed_source` defaults its ``boundary_condition``
parameter to ``"vacuum"`` and the internal helper
``_apply_default_bcs`` stamps :attr:`BC.vacuum <geometry.mesh.BC.vacuum>`
onto every face of the mesh that lacks an explicit BC declaration.
:class:`SNMesh` then resolves these to the ``"vacuum"`` kind string,
which the sweep reads directly.  In the 1-D cumprod path, the
recurrence starts from zero; in the 2-D wavefront path, the
reflective-partner copy is skipped, leaving incoming-face angular
fluxes at their zero initialisation (which is correct because no
code path writes the incoming-face slot of any ordinate except the
reflection step itself).

.. note::

   Before the BC infrastructure was introduced, the then-production
   ``transport_sweep`` entry accepted a ``boundary_condition: str``
   parameter directly.  That parameter has been removed --- BCs now flow
   through the mesh → SNMesh resolution path.  The description above
   reflects the current implementation.

**Measured convergence.**  With
:math:`\Sigma_t = 1\ \mathrm{cm^{-1}}`,
:math:`\Sigma_s = 0.5\ \mathrm{cm^{-1}}`,
:math:`L = 5\ \mathrm{cm}`, Gauss--Legendre :math:`S_{16}`:

.. list-table::
   :header-rows: 1
   :widths: 10 20 20

   * - :math:`n_{\rm cells}`
     - :math:`\|\phi_h - \phi_{\rm ex}\|_{L^{2}}`
     - measured order
   * - 10
     - :math:`2.17\!\times\!10^{-3}`
     - ---
   * - 20
     - :math:`5.40\!\times\!10^{-4}`
     - 2.01
   * - 40
     - :math:`1.35\!\times\!10^{-4}`
     - 2.00
   * - 80
     - :math:`3.37\!\times\!10^{-5}`
     - 2.00
   * - 160
     - :math:`8.42\!\times\!10^{-6}`
     - 2.00

Successive ratios hit :math:`4.00\pm0.02`, i.e. the measured order
is exactly the design order of diamond difference.  The L1 test
:func:`tests.sn.verification.mms.test_mms.test_sn_1d_slab_mms_converges_second_order`
asserts a slightly loose ``order > 1.9`` bracket to leave room for
round-off at the finest mesh.

**Risk points / things that can go wrong.**

- *Vacuum BC not honoured.*  If the reflective-partner copy is not
  skipped, incoming-face angular flux at the boundary is non-zero
  (the reflected outgoing from the opposite sweep) and the
  manufactured solution no longer satisfies the discrete problem.
  Symptom: :math:`\mathcal O(1)` error at the coarsest mesh; no
  convergence regardless of refinement.
- *Wrong normalisation for* :math:`Q_{\rm ext}`.  The solver's
  :math:`Q_{\rm aniso}` slot is divided by :math:`W` internally;
  the ansatz has a :math:`1/W` prefactor; the two must cancel.
  If the derivation forgets the :math:`W` cancellation, the
  measured flux is a factor of :math:`W` off but still converges at
  order 2 --- sneaky.  Guard: the second test in ``test_mms.py``
  cross-checks the algebraic symmetry of :eq:`sn-mms-qext`.
- *Non-smooth ansatz.*  A discontinuous material or a piecewise
  linear ansatz degrades the observed order to :math:`\mathcal O(h)`.
  The homogeneous sinusoid avoids both.
- *1-group vs multigroup.*  Because the manufactured flux is isotropic
  and there is no fission in the fixed-source problem, 1 group is
  sufficient --- the degeneracy warning about 1-group eigenvalue
  tests does not apply, since no :math:`\keff` enters.  Multigroup
  and heterogeneous MMS extensions are tracked as follow-ups for
  richer operator coverage.

**Follow-ups.**  MMS for :doc:`/theory/methods/method_of_characteristics`, diffusion,
and spherical / cylindrical curvilinear SN is tracked in GitHub
Issues (see ``type:feature level:L1``).  The curvilinear sweeps
need their own ansatz because their vacuum BC plumbing is not
yet wired up. **Heterogeneous and multigroup SN MMS is covered
by the next subsection.**


.. _sn-mms-heterogeneous-verification:

Heterogeneous MMS — 2-group continuous-:math:`\Sigma` slab
-----------------------------------------------------------

The homogeneous MMS case above verifies the Cartesian 1D SN
sweep for a *single-material* slab. To verify the multigroup
operator on a **heterogeneous** problem --- where each cell can
have different cross sections and the scatter matrix couples
groups across positions --- the Method of Manufactured Solutions
is extended in Phase 2.1a of the verification campaign with two
deliberate choices:

1. **Continuous (smooth)** :math:`\Sigma_{t,g}(x)` and
   :math:`\Sigma_{s,g\to g'}(x)` instead of piecewise-constant
   material regions. Discontinuous :math:`\Sigma` at interfaces
   that do not coincide with cell faces degrades diamond
   difference from :math:`\mathcal O(h^{2})` to
   :math:`\mathcal O(h)`, which would contaminate the
   spatial-convergence measurement with interface-treatment
   artefacts. With smooth :math:`\Sigma(x)` the diamond-
   difference operator hits its design :math:`\mathcal O(h^{2})`
   order exactly --- the convergence study becomes a clean test
   of the operator itself. This follows Salari & Knupp
   SAND2000-1444 §6, the canonical MMS reference for
   heterogeneous verification.
2. **Per-group amplitudes** :math:`\mathbf c = (c_1, c_2)` in
   the ansatz, so the scalar flux has a non-trivial group
   spectrum and the downscatter source term in the manufactured
   :math:`Q^{\text{ext}}` is non-zero. A bug that transposes
   the scatter matrix or drops a cross-group source term
   produces an incorrect :math:`\phi_2` that the convergence
   test catches immediately.

**Ansatz.**  The homogeneous ansatz carries over, now with a
per-group amplitude:

.. math::
   :label: sn-mms-hetero-psi

   \psi_{n,g}(x) \;=\; \frac{c_g}{W}\,A(x),
   \qquad A(x) \;=\; \sin\!\left(\frac{\pi x}{L}\right),

where :math:`W = \sum_n w_n` is the quadrature weight sum. The
scalar flux in each group is
:math:`\phi_g(x) = c_g\,A(x)`, so the amplitudes
:math:`\mathbf c` literally are the group fluxes at the slab
midpoint (where :math:`A` peaks). With
:math:`\mathbf c = (1.0, 0.3)` the two groups are linearly
independent and the downscatter coupling is visible.

Both groups share the same *spatial* mode :math:`\sin(\pi x/L)`
--- this is the fundamental mode of the bare slab and is exactly
the shape that emerges from separation of variables in the
diffusion limit. The heterogeneous SN problem would in general
have each group living in its own spatial harmonic, but we
*choose* the shared-mode ansatz as the manufactured target and
derive the non-trivial :math:`Q^{\text{ext}}` that makes it
satisfy the transport equation. The test then measures how well
the numerical SN sweep reproduces this prescribed shape.

**Manufactured source.**  Substituting :eq:`sn-mms-hetero-psi`
into the multigroup discrete-ordinates transport equation

.. math::

    \mu_n\,\frac{\partial\psi_{n,g}}{\partial x}
        + \Sigma_{t,g}(x)\,\psi_{n,g}
    \;=\; \frac{1}{W}\!\left(
        \sum_{g'}\Sigma_{s,g'\to g}(x)\,\phi_{g'}(x)
      + Q^{\text{ext}}_{n,g}(x)
    \right)

and solving algebraically for :math:`Q^{\text{ext}}`:

.. math::
   :label: sn-mms-hetero-qext

   Q^{\text{ext}}_{n,g}(x) \;=\;
       \mu_n\,c_g\,A'(x)
     + c_g\,\Sigma_{t,g}(x)\,A(x)
     \;-\; \sum_{g'}\Sigma_{s,g'\to g}(x)\,c_{g'}\,A(x).

The :math:`W` factor cancels between the ansatz's :math:`1/W`
prefactor and the solver's own :math:`1/W` convention on the
isotropic and anisotropic source slots, so :eq:`sn-mms-hetero-qext`
is the residual hand-delivered to the sweep without any
additional rescaling.

**Structure of the source.**  The streaming term
:math:`\mu_n\,c_g\,A'(x)` is odd in :math:`\mu` and carries the
only angular dependence, which is why SN MMS fundamentally
needs the per-ordinate ``Q_aniso`` plumbing path. The removal
term :math:`c_g\,\Sigma_{t,g}(x)\,A(x)` is diagonal in group
index. The **in-scatter** sum
:math:`\sum_{g'}\Sigma_{s,g'\to g}\,c_{g'}\,A(x)` is the only
term that couples groups, and for :math:`g=2` in the default
2-group setup it contributes
:math:`-\Sigma_{s,1\to 2}(x)\,c_1\,A(x)` --- the thermal source
depends on the fast amplitude through the downscatter cross
section, exactly as the multigroup scatter assembly in the
sweep does.

**Canonical cross sections.**  The reference uses smooth
profiles on :math:`[0, L]`:

.. math::

    \Sigma_{t,1}(x) &= 1.0 + 0.2\sin(\pi x/L), \\
    \Sigma_{t,2}(x) &= 2.0 + 0.3\cos(\pi x/L), \\
    \Sigma_{s,1\to 1}(x) &= 0.3 + 0.1\sin(\pi x/L), \\
    \Sigma_{s,1\to 2}(x) &= 0.2 + 0.05\sin(\pi x/L), \\
    \Sigma_{s,2\to 2}(x) &= 1.5 + 0.15\sin(\pi x/L), \\
    \Sigma_{s,2\to 1}(x) &= 0.

These give :math:`\Sigma_{a,1}(x) = 0.5 + 0.05\sin(\pi x/L) > 0`
trivially and
:math:`\Sigma_{a,2}(x) = 0.5 + 0.3\cos(\pi x/L) - 0.15\sin(\pi x/L)`,
bounded below by :math:`0.5 - \sqrt{0.3^{2} + 0.15^{2}} \approx
0.165 > 0`, so the cross sections are physical everywhere. The
scattering ratios :math:`c_g = \Sigma_{s,\text{tot},g}/\Sigma_{t,g}`
stay around :math:`0.5` for both groups, which means source
iteration converges geometrically at rate :math:`\sim 0.5^n`
per sweep.

**Per-cell material construction.**  The solver consumes the
continuous :math:`\Sigma(x)` by creating **one material per cell**
with cross sections evaluated at the cell centre
:math:`x_i = (x_{i-1/2} + x_{i+1/2})/2`. The midpoint rule for
the cell-average cross section is :math:`\mathcal O(h^{2})`-
accurate on smooth :math:`\Sigma`, matching the diamond-
difference design order and not degrading the measured
convergence rate. The number of materials scales with mesh
refinement, so each mesh in the convergence study builds a
fresh materials dictionary via
:meth:`orpheus.derivations.continuous.mms.sn.SNSlab2GHeterogeneousMMSCase.build_materials`.

**Measured convergence.**  With default parameters
(:math:`L = 5\,\text{cm}`, :math:`\mathbf c = (1.0, 0.3)`,
Gauss--Legendre :math:`S_{16}`):

.. list-table::
   :header-rows: 1
   :widths: 10 20 20 20

   * - :math:`n_{\rm cells}`
     - :math:`\|\phi_1 - \phi_{1,\rm ex}\|_{L^{2}}`
     - :math:`\|\phi_2 - \phi_{2,\rm ex}\|_{L^{2}}`
     - measured order
   * - 20
     - :math:`3.71\!\times\!10^{-4}`
     - :math:`3.38\!\times\!10^{-4}`
     - ---
   * - 40
     - :math:`9.25\!\times\!10^{-5}`
     - :math:`8.45\!\times\!10^{-5}`
     - 2.00
   * - 80
     - :math:`2.31\!\times\!10^{-5}`
     - :math:`2.11\!\times\!10^{-5}`
     - 2.00
   * - 160
     - :math:`5.78\!\times\!10^{-6}`
     - :math:`5.28\!\times\!10^{-6}`
     - 2.00

Both groups hit the design order independently, confirming
that the multigroup scatter coupling is correctly exercised.
The L1 test
:func:`tests.sn.verification.mms.test_mms_heterogeneous.test_sn_heterogeneous_mms_converges_second_order`
asserts ``> 1.9`` to leave round-off headroom at the finest
mesh.

**What this replaces.** Before Phase 2.1a, the heterogeneous
SN verification was
``orpheus.derivations.continuous.cases.sn._derive_sn_heterogeneous``, which
computed the reference :math:`k_{\text{eff}}` by running the
SN solver itself at four mesh refinements and Richardson-
extrapolating the eigenvalue sequence. That is a **T3 circular
self-test** in the verification-campaign taxonomy: the solver
verifies against its own extrapolated output, so any consistent
bug in the SN sweep that affects all mesh refinements the same
way is invisible to the test. The heterogeneous MMS reference
above breaks the circularity: the reference comes from the
manufactured-solution algebra, not from the solver.

**Complementary eigenvalue verification.** The MMS test
verifies the **spatial operator** on a heterogeneous problem
but does not exercise the eigenvalue iteration. Phase 2.1b
lands a Case singular-eigenfunction eigenvalue reference --- see
:ref:`sn-case-heterogeneous-verification` --- that restores
eigenvalue-heterogeneous coverage for the SN solver (T2
semi-analytical, from the first-order Boltzmann equation
itself, no diffusion approximation).


.. _sn-mms-2d-verification:

2D Cartesian MMS — separable sinusoidal ansatz
-----------------------------------------------

Phase 3.1 of the verification campaign extends the MMS spatial-operator
verification to **two Cartesian dimensions**.  The 1D slab MMS tests
verify the :math:`\mu\,\partial\psi/\partial x` streaming term in
isolation; this section adds :math:`\mu_y\,\partial\psi/\partial y`
and confirms that the 2D wavefront sweep
(:func:`orpheus.sn.loss_representation._sweep_jacobi`) with diamond-difference
closure achieves its design :math:`\mathcal O(h^{2})` convergence rate.

**Ansatz.**  On a rectangle :math:`[0, L_x] \times [0, L_y]` with
vacuum boundary conditions:

.. math::
   :label: sn-mms-2d-psi

   \psi_n(x, y) \;=\; \frac{1}{W}\,A(x, y),
   \qquad A(x, y) \;=\; \sin\!\left(\frac{\pi x}{L_x}\right)
                         \sin\!\left(\frac{\pi y}{L_y}\right).

The ansatz is **isotropic in angle** --- every ordinate carries the
same angular flux amplitude --- so the scalar flux recovered by any
quadrature set equals :math:`\phi(x, y) = A(x, y)` exactly.  This
design is deliberate: it isolates **spatial** discretisation error from
angular quadrature error, exactly as in the 1D case
(:eq:`sn-mms-psi`).

The separable sinusoidal ansatz vanishes on all four domain edges
(:math:`x = 0`, :math:`x = L_x`, :math:`y = 0`, :math:`y = L_y`),
so vacuum BCs are satisfied automatically for every ordinate.

**Manufactured source.**  Substituting :eq:`sn-mms-2d-psi` into the
2D Cartesian transport equation :eq:`transport-cartesian-2d` and
solving for the residual:

.. math::
   :label: sn-mms-2d-qext

   Q^{\text{ext}}_n(x, y) \;=\;
       \mu_{x,n}\,\frac{\partial A}{\partial x}
     + \mu_{y,n}\,\frac{\partial A}{\partial y}
     + (\Sigma_t - \Sigma_s)\,A(x, y)

where the partial derivatives of the separable ansatz are:

.. math::

   \frac{\partial A}{\partial x} =
       \frac{\pi}{L_x}\cos\!\left(\frac{\pi x}{L_x}\right)
       \sin\!\left(\frac{\pi y}{L_y}\right), \qquad
   \frac{\partial A}{\partial y} =
       \sin\!\left(\frac{\pi x}{L_x}\right)
       \frac{\pi}{L_y}\cos\!\left(\frac{\pi y}{L_y}\right).

The manufactured source :eq:`sn-mms-2d-qext` is angle-dependent through
:math:`\mu_{x,n}` and :math:`\mu_{y,n}` (streaming terms) and
angle-independent in the removal term :math:`(\Sigma_t - \Sigma_s) A`.
It enters the solver through the ``Q_aniso`` external source slot in
:func:`orpheus.sn.solve_sn_fixed_source`.

**Quadrature.**  2D problems use Lebedev spherical quadrature
(:meth:`Quadrature.lebedev(17)
<orpheus.numerics.quadrature.Quadrature.lebedev>`, order 17 = 110 ordinates).
Because the ansatz is isotropic in angle, the quadrature-level angular
integration is exact for *any* quadrature set --- the spatial
convergence study isolates spatial error exclusively.

**Measured convergence.**  Four mesh refinements on a
:math:`5 \times 5\,\text{cm}` square domain with
:math:`\Sigma_t = 1.0`, :math:`\Sigma_s = 0.5`:

.. list-table::
   :header-rows: 1

   * - :math:`n_x = n_y`
     - L2 error
     - Order
   * - 10
     - :math:`5.50 \times 10^{-3}`
     -
   * - 20
     - :math:`1.37 \times 10^{-3}`
     - 2.01
   * - 40
     - :math:`3.41 \times 10^{-4}`
     - 2.00
   * - 80
     - :math:`8.53 \times 10^{-5}`
     - 2.00

The measured order is indistinguishable from 2.00 across all
refinements, confirming that the 2D wavefront sweep preserves the
diamond-difference design order.

**Code pointers.**

- Derivation:
  :class:`orpheus.derivations.continuous.mms.sn.SN2DCartesianMMSCase` and
  :func:`orpheus.derivations.continuous.mms.sn.build_2d_cartesian_mms_case`.
- Test:
  :func:`tests.sn.verification.mms.test_mms_2d.test_sn_2d_cartesian_mms_converges_second_order`.
- Sweep:
  :func:`orpheus.sn.loss_representation._sweep_jacobi` (the 2D diamond-difference
  kernel verified by this test).

**Why this test matters.**  The existing 2D SN tests
(:mod:`tests.sn.sweep.cartesian_2d.test_discrete_ordinates_2d`) are L2 self-convergence
tests with real cross sections that verify the solver as a black box.
This MMS test is more incisive: it provides a **closed-form reference
flux** and asserts the **design convergence order** of the spatial
discretisation.  A bug that corrupts the 2D DD cell-average formula
(e.g. swapping :math:`\Delta x` and :math:`\Delta y`, mis-indexing the
wavefront anti-diagonal, or computing face fluxes with the wrong
sign) would break the :math:`\mathcal O(h^{2})` rate while possibly
still converging at some reduced order — the MMS test catches this
immediately, while a self-convergence test might not.

**Gotchas.**

- *Ordinates with* :math:`\mu_x = \mu_y = 0`.  The Lebedev set
  includes purely :math:`z`-directed ordinates.  For these, the
  streaming terms vanish, and the sweep reduces to
  :math:`\psi = Q/\Sigma_t`.  The manufactured source formula
  handles this correctly because both :math:`\mu_{x,n}` and
  :math:`\mu_{y,n}` multiply the gradient terms.
- *Aspect ratio.*  The test uses :math:`L_x = L_y` (square domain).
  A non-square domain would work identically — the separable ansatz
  is parameterised by :math:`L_x` and :math:`L_y` independently.
  Phase 3.2 extends to 2-group with heterogeneous materials (below).


.. _sn-mms-2d-2g-verification:

2D Cartesian 2-group heterogeneous MMS
----------------------------------------

Phase 3.2 combines the 2D geometry from Phase 3.1 with the
smooth-:math:`\Sigma` heterogeneous approach from Phase 2.1a.  The
cross sections are smooth 2D functions :math:`\Sigma(x, y)` so the
diamond-difference design order :math:`\mathcal O(h^{2})` is preserved
(no interface degradation).

**Ansatz.**  Per-group amplitudes :math:`c_g` with the same 2D shape:

.. math::
   :label: sn-mms-2d-2g-psi

   \psi_{n,g}(x, y) = \frac{c_g}{W}\,A(x, y), \qquad
   A(x, y) = \sin(\pi x/L_x)\,\sin(\pi y/L_y),

giving :math:`\phi_g(x, y) = c_g\,A(x, y)` with
:math:`\mathbf c = (1.0, 0.3)`.

**Manufactured source.**  From the 2D multigroup transport equation:

.. math::
   :label: sn-mms-2d-2g-qext

   Q^{\text{ext}}_{n,g}(x, y) =
       \mu_{x,n}\,c_g\,\partial_x A
     + \mu_{y,n}\,c_g\,\partial_y A
     + \Sigma_{t,g}(x, y)\,c_g\,A
     - \sum_{g'}\Sigma_{s,g'\to g}(x, y)\,c_{g'}\,A.

The thermal (:math:`g = 2`) source couples to :math:`c_1` through
the downscatter term :math:`\Sigma_{s,1\to 2}(x, y)\,c_1\,A`, which
exercises the multigroup scatter assembly in the 2D sweep.

**Cross-section profiles.**  The 2D functions extend the 1D
Phase-2.1a profiles (see :ref:`sn-mms-heterogeneous-verification`)
with a mild :math:`y`-dependent modulation:

- :math:`\Sigma_{t,1}(x,y) = 1.0 + 0.2\sin(\pi x/L_x) + 0.1\cos(\pi y/L_y)`
- :math:`\Sigma_{t,2}(x,y) = 2.0 + 0.3\cos(\pi x/L_x) + 0.1\sin(\pi y/L_y)`

Scattering cross sections carry a :math:`0.05\cos(\pi y/L_y)` modulation.
All :math:`\Sigma_a > 0` bounds from the 1D case are preserved because
the :math:`y`-modulation amplitudes (0.1, 0.05) are smaller than the
1D absorption margin (:math:`\sim 0.165`).

**Measured convergence.**  Four refinements on a :math:`5 \times 5` cm
square:

.. list-table::
   :header-rows: 1

   * - :math:`n_x = n_y`
     - L2 error (g=1)
     - Order (g=1)
     - L2 error (g=2)
     - Order (g=2)
   * - 10
     - :math:`3.79 \times 10^{-3}`
     -
     - :math:`2.85 \times 10^{-3}`
     -
   * - 20
     - :math:`9.41 \times 10^{-4}`
     - 2.01
     - :math:`7.09 \times 10^{-4}`
     - 2.01
   * - 40
     - :math:`2.35 \times 10^{-4}`
     - 2.00
     - :math:`1.77 \times 10^{-4}`
     - 2.00
   * - 80
     - :math:`5.87 \times 10^{-5}`
     - 2.00
     - :math:`4.42 \times 10^{-5}`
     - 2.00

Both groups achieve the design :math:`\mathcal O(h^{2})` rate.

**Code pointers.**

- Derivation:
  :class:`orpheus.derivations.continuous.mms.sn.SN2DCartesian2GHeterogeneousMMSCase`
  and :func:`orpheus.derivations.continuous.mms.sn.build_2d_cartesian_heterogeneous_mms_case`.
- Test:
  :func:`tests.sn.verification.mms.test_mms_2d.test_sn_2d_cartesian_2g_heterogeneous_mms_converges_second_order`.


.. _sn-mms-p1-verification:

P1 anisotropic scattering MMS
-------------------------------

Phase 3.5 verifies that the P\ :sub:`N` anisotropic scattering
source assembly (:ref:`pn-scattering`) preserves
:math:`\mathcal O(h^{2})` convergence. All previous MMS tests use
isotropic (P0) scattering; this test exercises the P1 slot
:math:`\Sigma_s^{(1)}` through a weakly angle-dependent ansatz.

**Ansatz.** On a 1D vacuum-BC slab :math:`[0, L]`:

.. math::
   :label: sn-mms-p1-psi

   \psi_n(x) = \frac{1}{W}\bigl(A(x) + \alpha\,\mu_n\,B(x)\bigr)

with :math:`A(x) = B(x) = \sin(\pi x/L)` and small
:math:`\alpha = 0.1`. The scalar flux is :math:`\phi(x) = A(x)`
(the :math:`\mu`-odd term integrates to zero), and the P1 current
is :math:`J(x) = \alpha\,B(x)/3` (using
:math:`\sum w_n\mu_n^2 = 2/3` for Gauss–Legendre on
:math:`[-1, 1]`).

**Manufactured source.** Substituting :eq:`sn-mms-p1-psi` into
the 1D transport equation with P1 scattering and solving for
the residual:

.. math::
   :label: sn-mms-p1-qext

   Q^{\text{ext}}_n(x) =
       \mu_n\,A'(x)
     + (\Sigma_t - \Sigma_s^{(0)})\,A(x)
     + \alpha\,\mu_n\,(\Sigma_t - \Sigma_s^{(1)})\,B(x)
     + \alpha\,\mu_n^2\,B'(x).

The first two terms are the isotropic MMS source from
:eq:`sn-mms-qext`. The third term comes from the P1 scattering
slot :math:`3\,\Sigma_s^{(1)}\,\mu_n\,J(x)` in the transport
equation, and the fourth from the :math:`\mu_n`-weighted
streaming of :math:`B(x)`.

**Measured convergence.** Four refinements with
:math:`\Sigma_t = 1.0`, :math:`\Sigma_s^{(0)} = 0.5`,
:math:`\Sigma_s^{(1)} = 0.2`, :math:`\alpha = 0.1`:

.. list-table::
   :header-rows: 1

   * - :math:`n_{\text{cells}}`
     - L2 error
     - Order
   * - 20
     - :math:`6.15 \times 10^{-4}`
     -
   * - 40
     - :math:`1.53 \times 10^{-4}`
     - 2.00
   * - 80
     - :math:`3.84 \times 10^{-5}`
     - 2.00
   * - 160
     - :math:`9.59 \times 10^{-6}`
     - 2.00

**Code pointers.**

- Derivation:
  :class:`orpheus.derivations.continuous.mms.sn.SNP1AnisoMMSCase` and
  :func:`orpheus.derivations.continuous.mms.sn.build_p1_aniso_mms_case`.
- Test:
  :func:`tests.sn.verification.mms.test_mms_aniso.test_sn_p1_aniso_mms_converges_second_order`.
- P1 assembly:
  :meth:`orpheus.sn.solver.SNSolver._build_aniso_scattering`.


.. _sn-mms-curvilinear-isotropic-verification:

Curvilinear isotropic MMS — radial DD-closure probe
----------------------------------------------------

Phase 3.4 of the verification campaign extends the slab MMS
(:eq:`sn-mms-psi` / :eq:`sn-mms-qext`) to 1-D **spherical** and
1-D **cylindrical** geometries with the simplest non-trivial trial
solution that respects the vacuum-at-outer and symmetry-at-origin
boundary conditions: an **isotropic** ansatz
:math:`\psi_n(r) = A(r)/W`.  By construction the angular
redistribution operator vanishes on this ansatz
(:math:`(1-\mu^2)/r \cdot \partial\psi/\partial\mu = 0` for the
sphere; :math:`-(1/r)\,\partial(\xi\psi)/\partial\varphi = 0` for
the cylinder), so the only spatial-discretisation error that drives
the measured convergence rate is the **radial DD closure**.  The
isotropic case is therefore the focused L1 probe for the
streaming + removal path; the angular redistribution path is
covered by the companion anisotropic case
(:ref:`sn-mms-curvilinear-aniso-verification` below — a deliberate
pairing that defeats the ``vv-principles`` Mode 7 "MMS
simplification bias" failure mode).

**Spherical isotropic ansatz.**  For a vacuum-BC sphere of radius
:math:`R` with reflective inner BC at :math:`r=0` in one energy
group, pick

.. math::
   :label: sn-mms-spherical-psi

   \psi_n(r) = \frac{1}{W}\,A(r),
   \qquad A(r) = \sin\!\left(\frac{\pi r}{R}\right),

with :math:`W = \sum_n w_n = 2` for symmetric Gauss--Legendre.
Because :math:`A(0) = A(R) = 0`, every ordinate vanishes at both
the symmetry centre and the vacuum outer face — both BC kinds are
satisfied automatically.  Since :math:`\psi_n` is independent of
ordinate, the scalar flux recovered by any quadrature order is
exactly :math:`\phi(r) = A(r)`.

**Spherical manufactured source.**  Substituting
:eq:`sn-mms-spherical-psi` into :eq:`transport-spherical` and
using that :math:`(1-\mu^2)\,\partial_\mu\psi/r \equiv 0` for an
isotropic flux gives

.. math::
   :label: sn-mms-spherical-qext

   Q^{\text{ext}}_n(r)
        = \mu_n\,A'(r)
        + (\Sigma_t - \Sigma_s)\,A(r)
        = \mu_n\,\frac{\pi}{R}\cos\!\left(\frac{\pi r}{R}\right)
          + (\Sigma_t - \Sigma_s)\sin\!\left(\frac{\pi r}{R}\right).

This is structurally identical to the slab source
:eq:`sn-mms-qext` with :math:`x \to r` — the spherical
:math:`(2/r)\partial_r` curvature term and the angular
redistribution term both vanish on the isotropic ansatz, leaving
the per-ordinate streaming + removal balance as the residual.

**Cylindrical isotropic ansatz.**  The radial direction cosine for
1-D cylindrical is :math:`\eta_n = \sin\theta_n \cos\varphi_n`.
Use

.. math::
   :label: sn-mms-cylindrical-psi

   \psi_n(r) = \frac{1}{W}\,A(r),
   \qquad A(r) = \sin\!\left(\frac{\pi r}{R}\right),

with the same :math:`W = \sum_n w_n` for the cylindrical Product or
LS quadrature.  Symmetric Product quadrature gives
:math:`\sum_n w_n \eta_n = 0`, so :math:`\phi(r) = A(r)` exactly.

**Cylindrical manufactured source.**

.. math::
   :label: sn-mms-cylindrical-qext

   Q^{\text{ext}}_n(r)
        = \eta_n\,A'(r) + (\Sigma_t - \Sigma_s)\,A(r).

The cylindrical curvature term :math:`-(1/r)\,\partial(\xi\psi)/\partial\varphi`
vanishes by isotropy of :math:`A(r)`, the same way the spherical
:math:`(1-\mu^2)/r \cdot \partial_\mu\psi` vanishes; the radial
streaming :math:`\eta_n A'(r)` and the removal
:math:`(\Sigma_t - \Sigma_s)A(r)` carry the residual.

**Risk point — Mode 7 ansatz bias.**  Per ``vv-principles`` failure
Mode 7 ("MMS simplification bias"), the isotropic ansatz is
deliberately structured to NULL the angular redistribution path.
A passing :math:`\mathcal{O}(h^2)` convergence here is necessary
evidence for the radial DD closure but it is *not* sufficient for
the full curvilinear sweep — ERR-026 (the curvilinear sweep WDD
flux-shape bug) is mathematically invisible to this MMS because the
redistribution term that ERR-026 lives on cancels by ansatz
construction.  The companion anisotropic case
(:ref:`sn-mms-curvilinear-aniso-verification`) is the load-bearing
sufficient evidence for the full sweep; both are required.

**Code pointers.**

- Derivation:
  :class:`orpheus.derivations.continuous.mms.sn.SNSphericalMMSCase`,
  :class:`orpheus.derivations.continuous.mms.sn.SNCylindricalMMSCase`,
  :func:`orpheus.derivations.continuous.mms.sn.build_spherical_mms_case`,
  :func:`orpheus.derivations.continuous.mms.sn.build_cylindrical_mms_case`.
- Tests:
  :func:`tests.sn.verification.mms.test_mms_curvilinear.test_sn_spherical_mms_converges_second_order`
  (sphere) and
  :func:`tests.sn.verification.mms.test_mms_curvilinear.test_sn_cylindrical_mms_converges_second_order`
  (cylinder), both ``catches("ERR-058")``.  **Their ``xfail`` markers
  came off 2026-06-12 with the ERR-058 closure-seed fix** (Issue #195).
  Post-fix the ladders are clean second-order with SI :math:`\equiv`
  Krylov bit-identical — sphere ``[1.49e-2, 3.73e-3, 9.28e-4, 2.31e-4,
  5.74e-5]`` (orders 2.00–2.01), cylinder ``[2.16e-3, 5.39e-4, 1.35e-4,
  3.37e-5]`` (orders 2.00); the magnitude band
  :math:`[10^{-8}, 10^{-3}]` is met (sphere :math:`n_x\ge 80`,
  cylinder :math:`n_x\ge 40`).  Through the bug era (Wave E Round 3
  2026-05 → ERR-058) they were ``xfail(strict=True)`` — the
  now-superseded "pre-asymptotic transient" reading; see
  :ref:`sn-err-058-closure-seed-closeout`.


.. _sn-mms-curvilinear-aniso-verification:

Curvilinear anisotropic MMS — angular redistribution probe
-----------------------------------------------------------

Phase 3.6 closes the **angular-redistribution coverage gap** in the
curvilinear MMS verification chain. The existing isotropic
1D-spherical (:class:`SNSphericalMMSCase`) and 1D-cylindrical
(:class:`SNCylindricalMMSCase`) MMS cases use the ansatz
:math:`\psi_n(r) = A(r)/W` (no :math:`\mu`-dependence). For that
ansatz, **the angular-redistribution operator is identically zero**:
:math:`(1-\mu^2)/r \cdot \partial\psi/\partial\mu = 0` for the sphere,
:math:`-(1/r)\,\partial(\xi\psi)/\partial\varphi = 0` for the
cylinder. The hardest math the curvilinear sweep performs — where
ERR-026 (curvilinear sweep WDD wrong fixed point) lives — is
mathematically invisible to the isotropic MMS because it cancels by
ansatz construction.

This is the ``vv-principles`` failure mode #7 ("MMS simplification
bias") — the MMS test cannot catch a bug class because the ansatz
nulls it. The defence is **not**
to replace the isotropic case (it remains the right probe for the
non-redistribution paths) but to **pair** it with a companion case
whose ansatz activates redistribution. The two cases together let a
narrow-down diagnosis route a failing convergence rate to either
the streaming/removal path (only isotropic fails) or the
redistribution path (only anisotropic fails).

**Spherical anisotropic ansatz**

.. math::
   :label: sn-mms-spherical-aniso-psi

   \psi_n(r) = \frac{1}{W}\bigl(A(r) + B(r)\,\mu_n\bigr),
   \qquad
   A(r) = \sin\!\left(\frac{\pi r}{R}\right),
   \qquad
   B(r) = \frac{r}{R}\Bigl(1 - \frac{r}{R}\Bigr)
            \cos\!\left(\frac{\pi r}{R}\right).

Both :math:`A` and :math:`B` vanish at :math:`r \in \{0, R\}`, so
**every** ordinate satisfies the symmetry BC at :math:`r=0` and the
vacuum BC at :math:`r=R`, regardless of the sign of :math:`\mu_n`.
The :math:`B(r)\,\mu_n` coefficient is non-trivial: ordinates with
opposite sign of :math:`\mu_n` differ in sign of the
angular-flux contribution, but both still vanish at the boundaries.

The choice :math:`B(r) = (r/R)(1-r/R)\cos(\pi r/R)` is **not**
algebraically reducible to a multiple of :math:`A(r)` — the
:math:`(r/R)(1-r/R)` envelope and the :math:`\cos(\pi r/R)` factor
produce a derivative :math:`B'(r)` whose extrema do not co-locate
with :math:`A'(r)`'s extrema, so the redistribution term
:math:`(1-\mu_n^2)\,B/r` cannot be absorbed into a renormalisation
of the streaming term.

The discrete scalar flux is :math:`\phi(r) = A(r)` because
:math:`\sum_n w_n \mu_n = 0` for any symmetric Gauss-Legendre
quadrature — the :math:`B \mu` term integrates to zero in the
scalar moment.

**Spherical manufactured source**

Substituting :eq:`sn-mms-spherical-aniso-psi` into
:eq:`transport-spherical` and solving for the residual:

.. math::
   :label: sn-mms-spherical-aniso-qext

   Q^{\text{ext}}_n(r) =
        \mu_n\,A'(r)
      + \mu_n^2\,B'(r)
      + (1 - \mu_n^2)\,\frac{B(r)}{r}
      + (\Sigma_t - \Sigma_s)\,A(r)
      + \Sigma_t\,\mu_n\,B(r).

The first and fourth terms are the isotropic-MMS source from
:eq:`sn-mms-qext` adapted to spherical. The **second term**
(:math:`\mu_n^2 B'(r)` — :math:`\mu`-weighted streaming of the
anisotropic profile) and the **third term**
(:math:`(1-\mu_n^2)\,B/r` — angular redistribution) are
load-bearing: they are precisely what the isotropic case lacks.
The fifth (:math:`\Sigma_t\,\mu_n B`) comes from the removal
operator acting on the :math:`B \mu` part of :math:`\psi_n`.

**Cylindrical anisotropic ansatz**

The radial direction cosine for cylindrical 1D is :math:`\eta_n =
\sin\theta_n \cos\varphi_n`; the azimuthal partner that drives
the redistribution is :math:`\xi_n = \sin\theta_n \sin\varphi_n`.
Use:

.. math::
   :label: sn-mms-cylindrical-aniso-psi

   \psi_n(r) = \frac{1}{W}\bigl(A(r) + B(r)\,\eta_n\bigr),

with the same :math:`A(r),\,B(r)` shapes. Symmetric ProductQuadrature
gives :math:`\sum_n w_n \eta_n = 0`, so :math:`\phi(r) = A(r)`.

**Cylindrical manufactured source**

Substituting :eq:`sn-mms-cylindrical-aniso-psi` into
:eq:`transport-cylindrical` (treating :math:`\eta_n` and
:math:`\xi_n` as the :math:`\varphi`-dependent functions
:math:`\sin\theta\cos\varphi` and :math:`\sin\theta\sin\varphi`)
and solving for the residual:

.. math::
   :label: sn-mms-cylindrical-aniso-qext

   Q^{\text{ext}}_n(r) =
        \eta_n\,A'(r)
      + \eta_n^2\,B'(r)
      + \xi_n^2\,\frac{B(r)}{r}
      + (\Sigma_t - \Sigma_s)\,A(r)
      + \Sigma_t\,\eta_n\,B(r).

The :math:`\xi_n^2\,B/r` redistribution term is the cylindrical analog
of the sphere's :math:`(1-\mu_n^2)\,B/r`. Both come from the
same operator — angular redistribution of the linearly-:math:`\mu`
(or linearly-:math:`\eta`) ansatz — and both vanish for any
isotropic ansatz.

**Spatial-convergence claims.**  Diamond-Difference is design-order
:math:`\mathcal{O}(h^2)` in the cell width (:eq:`dd-cartesian-1d` /
:eq:`dd-curvilinear-scalar`); the curvilinear anisotropic L1 claim
asserts that the **measured** scalar-flux error against the
manufactured solution :eq:`sn-mms-spherical-aniso-psi` /
:eq:`sn-mms-cylindrical-aniso-psi` falls at the same rate.  For the
sphere,

.. math::
   :label: sn-mms-spherical-aniso-spatial-convergence

   \bigl\|\phi_h(r) - A(r)\bigr\|_{L^2(\Omega)}
        \;=\; \mathcal{O}(h^2)
        \qquad \text{as } h = R/n_x \to 0\,,

with the convergence ORDER (slope of :math:`\log\|\phi_h - A\|`
versus :math:`\log h` over the last two mesh halvings) the L1
acceptance criterion ``min(orders[-2:]) > 1.9``.  The cylindrical
analogue,

.. math::
   :label: sn-mms-cylindrical-aniso-spatial-convergence

   \bigl\|\phi_h(r) - A(r)\bigr\|_{L^2(\Omega)}
        \;=\; \mathcal{O}(h^2)
        \qquad \text{as } h = R/n_x \to 0\,,

uses the same acceptance criterion on the cylindrical-aniso
ansatz.  Both labels are consumed by the
:file:`tests/sn/verification/mms/test_curvilinear_aniso_convergence.py`
gate-3 tests, which **stay ``xfail``** — but, post-ERR-058 (Issue
#195), no longer for the wrong-fixed-point reason.  The ERR-058
closure-seed fix recovered :math:`\mathcal{O}(h^2)` *spatial*
convergence (the isotropic ladders are clean second-order; see
:ref:`sn-err-058-closure-seed-closeout`).  These **anisotropic** rows
remain xfail because the angle-varying ansatz hits the
**fixed-quadrature angular floor** of the half-angle thread
interpolation: under spatial refinement at fixed quadrature the error
converges to a floor (sphere S16 :math:`\approx 7\mathrm{e}{-4}`,
cylinder :math:`n_\mu{=}4` :math:`\approx 1.9\mathrm{e}{-2}`) that
drops only under *quadrature* refinement.  Their markers are re-pinned
to the `Issue #229
<https://github.com/deOliveira-R/ORPHEUS/issues/229>`_
quadrature-aware retune (the regression gate that flips them to
unexpected-pass when the retune lands).  See
:ref:`sn-err-058-aniso-floor` for the floor-vs-quadrature evidence.

**Verification chain (Branch 1 / Branch 2)**

Per the ``algebra-of-record`` discipline (Branch-1 SymPy reference,
Branch-2 numpy production, structurally-independent L1 cross-check):

- **Branch 1 (SymPy)**:
  :func:`orpheus.derivations.continuous.mms.sn.derive_spherical_anisotropic_mms`
  and
  :func:`orpheus.derivations.continuous.mms.sn.derive_cylindrical_anisotropic_mms`
  substitute the ansatz into the transport operator symbolically and
  prove ``simplify(LHS - RHS) == 0``. Foundation tests:
  :file:`tests/derivations/test_sn_mms_anisotropic_symbolic.py` (one
  ``@pytest.mark.foundation`` test per ``derive_*`` function).
- **Branch 2 (vectorised numpy)**:
  :class:`orpheus.derivations.continuous.mms.sn.SNSphericalAnisotropicMMSCase`
  and
  :class:`orpheus.derivations.continuous.mms.sn.SNCylindricalAnisotropicMMSCase`
  evaluate :eq:`sn-mms-spherical-aniso-qext` and
  :eq:`sn-mms-cylindrical-aniso-qext` per ordinate using vectorised
  numpy.
- **L1 cross-check (the gate)**: the Branch-2 numerical
  :math:`Q^{\text{ext}}_n(r_i)` agrees with Branch-1 SymPy-evaluated
  :math:`Q^{\text{ext}}_n(r_i)` (via :func:`sympy.lambdify`) to
  :math:`\sim 10^{-16}` (max absolute) on a sample mesh in both
  geometries. Tested in
  :func:`tests.derivations.test_sn_mms_anisotropic_symbolic.test_spherical_aniso_numerical_qext_matches_sympy`
  and the cylindrical sibling.

**Code pointers**

- Derivations: :class:`orpheus.derivations.continuous.mms.sn.SNSphericalAnisotropicMMSCase`,
  :class:`orpheus.derivations.continuous.mms.sn.SNCylindricalAnisotropicMMSCase`,
  :func:`orpheus.derivations.continuous.mms.sn.build_spherical_anisotropic_mms_case`,
  :func:`orpheus.derivations.continuous.mms.sn.build_cylindrical_anisotropic_mms_case`.
- Symbolic factory:
  :func:`orpheus.derivations.continuous.mms.sn._spherical_anisotropic_symbolic`,
  :func:`orpheus.derivations.continuous.mms.sn._cylindrical_anisotropic_symbolic`.
- Foundation tests:
  :file:`tests/derivations/test_sn_mms_anisotropic_symbolic.py`.
- Consumer L1 convergence test (Phase-0 work, separate branch):
  ``tests/sn/_l1/test_mms_spherical_anisotropic_dd_convergence_O_h2.py``
  (planned).


.. _sn-composite-fixed-source:

The composite fixed-source API — :math:`q = q_{\text{bulk}} \oplus q_\partial`
------------------------------------------------------------------------------

.. admonition:: Key Facts
   :class: important

   - **A fixed source is a source EVERYWHERE.** The right-hand side of
     a fixed-source transport problem is not just a volumetric bulk
     source — it is a source on the whole phase space: a bulk
     :math:`q_{\text{bulk}}` *and* a boundary (prescribed-inflow)
     :math:`q_\partial`. ORPHEUS represents this as the composite
     ``q = q_bulk ⊕ q_∂``, the direct sum of the two role-typed leaves.
   - **The carrier is the object we already have.**
     :func:`~orpheus.sn.solver.solve_sn_fixed_source` accepts the
     composite as a
     :class:`~orpheus.transport.timed_full_field.TimedFullField` — the
     **same** typed direct-sum carrier the SI / Krylov inner already
     flows through internally. This is *not* a new type; it is
     ergonomics to **generate** the right object (Cardinal Rule 2 — we
     already have the right concepts).
   - **A source, role-distinguished from a flux by its leaf types.**
     The composite's bulk leaf is an
     :class:`~orpheus.transport.source_sinks.AngularSourceSink` and its
     boundary leaf a
     :class:`~orpheus.transport.source_sinks.AngularBoundarySourceSink` — the
     *source* column of the role grid (see :ref:`bc-extraction-operator-output-typing`).
     The iterate / solution it produces is a *flux* (``AngularFlux`` ⊕
     ``AngularBoundaryFlux``). Same carrier shape, different role; the class
     gate keeps source and flux arithmetic from silently mixing.
   - **The legacy array is the bulk-only / vacuum special case.** Passing
     the historical ``(N, ng, nx, ny)`` ndarray is *exactly* the
     composite with an all-zero (vacuum) boundary. All 37 pre-existing
     callers keep working bit-unchanged.
   - **One construction point.** The private helper
     :func:`~orpheus.sn.solver._build_fixed_source_rhs` is the single
     place the RHS composite is built (Cardinal Rule 2 — it collapsed a
     ``q_ext_composite`` build that previously lived in *both* the SI
     and Krylov inner paths).
   - **The ergonomic boundary generator.**
     :meth:`~orpheus.transport.source_sinks.AngularBoundarySourceSink.prescribed_inflow`
     writes ONLY the inflow ordinate slots of the named faces (outflow
     slots of a prescribed inflow are physically meaningless →
     unrepresentable by construction, ``coding-elegance`` Pattern 4),
     leaving everything else zero. It is the known-per-face-array route;
     the lazy ``InflowSourceSpec``-recipe route (``from_spec``) is a
     distinct, still-deferred bridge.

The fixed-source right-hand side is a source on the whole phase space
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A fixed-source SN problem solves the affine within-group system
:eq:`si-within-group-operator-eq`

.. math::

   (L + C - S - B)\,\psi = q,

and the right-hand side :math:`q` is **not** a bulk volumetric source
alone. It has two pieces, one per phase-space locus:

* the **bulk** source :math:`q_{\text{bulk}}(\vec r, \hat\Omega, g)` —
  the per-ordinate volumetric external source :math:`Q^{\text{ext}}_n`
  on every cell;
* the **boundary** source :math:`q_\partial` — the prescribed inflow,
  the inhomogeneous term :math:`q` of the affine boundary law
  :eq:`affine-bc-form` :math:`\gamma_-\psi = R\,G\,\gamma_+\psi + q`,
  living on the inflow ordinate slots of the boundary trace.

A vacuum boundary is simply :math:`q_\partial \equiv 0`; a non-vacuum
prescribed inflow is a non-zero :math:`q_\partial`. The natural object
is therefore the **direct sum** of the two:

.. math::

   q \;=\; q_{\text{bulk}} \,\oplus\, q_\partial,

an object that "represents the source everywhere". ORPHEUS already has
exactly this carrier: the
:class:`~orpheus.transport.timed_full_field.TimedFullField`, the typed
bulk⊕boundary(⊕history) direct sum that the within-group SI and Krylov
inner paths *already* pass around (the matvec
:math:`(L+C)\psi - (S+B)\psi - F\psi` and the SI rhs
:math:`F\psi + (S+B)\psi + q_{\text{ext}}` are CLOSED ``TimedFullField``
sums). The field-role-typing work did **not** introduce a new source
type — it surfaced the carrier we already had and added the ergonomics
to *generate* it (Cardinal Rule 2: we have the right object, we just
need a better way to build it).

Source vs flux — same carrier, different role
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The composite source and the angular-flux solution share the
``TimedFullField`` carrier *shape* but differ in their leaf **types**,
which encode the role:

.. list-table:: The composite carrier's two roles
   :header-rows: 1
   :widths: 22 39 39

   * - Locus
     - Source role (the RHS ``q``)
     - Flux role (the iterate / solution ``ψ``)
   * - bulk
     - :class:`~orpheus.transport.source_sinks.AngularSourceSink`
     - :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
   * - boundary
     - :class:`~orpheus.transport.source_sinks.AngularBoundarySourceSink`
     - :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`

The role-leaf types are the gate. A *source* and a *flux* are never the
same field even when they ride the same carrier — the
:class:`~orpheus.transport.timed_full_field.TimedFullField` class gate
(via :class:`~orpheus.numerics.field.Field`) rejects
``AngularSourceSink ± AngularFlux`` and the boundary analogue, so the
"RHS is a source, the iterate is a flux" distinction is illegal to mix
by construction. The completed boundary role grid mirrors the bulk
exactly (:ref:`bc-extraction-operator-output-typing`): an operator's
``.apply`` output is a *source/sink* (:math:`A\psi`), its ``.solve``
output is a *flux* (the swept solution trace), and a ``from_balance``
defect is a *residual*. ``q_\partial`` is a ``AngularBoundarySourceSink``
because a prescribed inflow IS a source added to :math:`\gamma_-\psi`,
not the swept solution.

The two accepted forms of ``external_source``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:func:`~orpheus.sn.solver.solve_sn_fixed_source` accepts the
``external_source`` argument in either of two forms, normalised by the
single helper :func:`~orpheus.sn.solver._build_fixed_source_rhs`:

#. **A bare ``np.ndarray`` of shape** :math:`(N, n_g, n_x, n_y)` — the
   per-ordinate-density **bulk** source only, with a **vacuum**
   boundary. This is the original form, and it is *exactly* the
   composite with an all-zero boundary leaf
   (``AngularBoundarySourceSink.zeros_on(sn_mesh)``). Every one of the 37
   pre-existing callers passes this form and keeps working bit-for-bit
   unchanged (the vacuum path is verified bit-identical).
#. **A full** :class:`~orpheus.transport.timed_full_field.TimedFullField`
   **composite** ``q = q_bulk ⊕ q_∂`` — the route for a **non-vacuum
   prescribed inflow**. Its leaf values are re-homed onto the solve's
   own ``sn_mesh``: the trace / grid layout is deterministic from
   ``(mesh, quadrature, materials)``, so this is an exact values-copy
   onto the solve's mesh instance, required because the within-group
   operators are built on ``sn_mesh`` and ``TimedFullField`` algebra
   enforces mesh identity.

.. code-block:: python

   from orpheus.sn import solve_sn_fixed_source
   from orpheus.sn.mesh.augmented_mesh import SNMesh
   from orpheus.transport.source_sinks import (
       AngularSourceSink, AngularBoundarySourceSink,
   )
   from orpheus.transport.timed_full_field import TimedFullField

   sn = SNMesh(mesh, quadrature, materials)

   # Bulk volumetric source, per-ordinate density (N, ng, nx, ny).
   q_bulk = AngularSourceSink.from_mesh(Q_ext, sn)
   # Prescribed inflow: only the named faces' inflow ordinate slots.
   q_bndry = AngularBoundarySourceSink.prescribed_inflow(
       sn, {"xmin": gamma_minus_xmin, "xmax": gamma_minus_xmax},
   )
   q = TimedFullField(bulk=q_bulk, boundary=q_bndry)

   result = solve_sn_fixed_source(materials, mesh, quadrature, q)

The legacy ``solve_sn_fixed_source(materials, mesh, quadrature, Q_ext)``
with a bare ``Q_ext`` array is identical to the above with a vacuum
``q_bndry`` (``AngularBoundarySourceSink.zeros_on(sn)``).

The single construction point — Cardinal Rule 2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before the field-role-typing work, the SI inner and the Krylov inner
each built their own ``q_ext_composite`` from the bulk array (the same
``AngularSourceSink.from_isotropic`` / ``from_mesh`` projection paired
with a zero boundary). That was a shared concept living in two places —
precisely the smell Cardinal Rule 2 flags.
:func:`~orpheus.sn.solver._build_fixed_source_rhs` collapses both into
one construction point: ``solve_sn_fixed_source`` calls it once, and
**both** inner paths consume what it returns. The helper:

* validates the bulk shape against :math:`(N, n_g, n_x, n_y)` (Issue
  #196 PR-INDEX-5 principled layout — the ``g`` axis directly after
  ``N``);
* for a bare array, pairs the bulk
  :class:`~orpheus.transport.source_sinks.AngularSourceSink` with a
  vacuum ``AngularBoundarySourceSink``;
* for a composite, re-homes the leaf values onto the solve's
  ``sn_mesh`` (with a layout-size guard on the boundary trace), and
  raises a descriptive ``ValueError`` if the composite was built on an
  incompatible mesh / quadrature / materials.

The validation, the projection, and the vacuum-default boundary now
live in exactly one function. The SI and Krylov paths differ only in
the inner solve they run, not in how the RHS is assembled.

The ergonomic boundary generator — ``prescribed_inflow``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generating the boundary leaf :math:`q_\partial` is the part the
ergonomics target. The classmethod
:meth:`AngularBoundarySourceSink.prescribed_inflow(mesh, {face: (N, ng) values}) <orpheus.transport.source_sinks.AngularBoundarySourceSink.prescribed_inflow>`
builds the prescribed inflow from known per-face arrays:

* for each face named in the mapping, it writes ONLY the **inflow**
  ordinate slots from the given :math:`(N, n_g)` array;
* **every other slot is left zero** — the outflow ordinate slots of a
  named face, and every slot of an unnamed (vacuum) face.

**Why outflow is unrepresentable (Pattern 4).** The outflow ordinate
slots of a *prescribed-inflow source* are physically meaningless: the
sweep determines the outflow trace, the source does not. Writing them
would be an illegal state. Rather than accept-then-ignore them,
:meth:`~orpheus.transport.source_sinks.AngularBoundarySourceSink.prescribed_inflow`
makes the illegal state **unrepresentable by construction** — it reads
the inflow ordinate index set
(:meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.inflow_indices_for_face`)
and copies *only* those rows, so an accidentally-populated outflow row
in the caller's array simply cannot reach the field. This is
``coding-elegance`` Pattern 4 (illegal states unrepresentable). It
supersedes the ``zeros_on`` + nested
``face_view(face)[inflow] = …`` slot-fill loop that every
prescribed-inflow consumer (the non-vacuum MMS, the splitting-invariance
probe) previously hand-rolled — the single source of truth for
materialising a prescribed inflow onto the trace (Cardinal Rule 2).

**The recipe → snapshot distinction (vs ``from_spec``).** There are two
distinct routes to a boundary source, related as *recipe → snapshot*,
not as duplicates:

.. list-table:: Two routes to ``q_∂`` — known arrays vs lazy recipe
   :header-rows: 1
   :widths: 28 36 36

   * -
     - :meth:`~orpheus.transport.source_sinks.AngularBoundarySourceSink.prescribed_inflow`
     - ``AngularBoundarySourceSink.from_spec`` (deferred)
   * - Input
     - known per-face ``(N, ng)`` arrays
     - a lazy
       :class:`~orpheus.geometry.boundary._source.InflowSourceSpec`
       recipe (``evaluate(shape) -> ndarray``)
   * - When
     - the inflow values are already computed (the MMS case)
     - the inflow is described by a per-face recipe evaluated on demand
   * - Status
     - **shipped** — the route the 4.6 MMS and the T4 probe use
     - **deferred** (``unify-after-two`` — no recipe-driven consumer
       that drives a typed boundary-source sweep yet)

The 4.6 MMS uses
:meth:`~orpheus.transport.source_sinks.AngularBoundarySourceSink.prescribed_inflow`
because it has explicit per-face arrays
(:math:`\gamma_-\psi = (A + \mu_n B)/W`); it does not need the lazy
recipe bridge, which waits for its first genuine consumer per the
``unify-after-two`` discipline.

Why this is ergonomics, not new types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The entire change is the **generation** of objects that already
existed in the codebase, not the introduction of new ones:

* The carrier :class:`~orpheus.transport.timed_full_field.TimedFullField`
  pre-dates this work — it is what the inner solve already flows.
* The leaf types
  :class:`~orpheus.transport.source_sinks.AngularSourceSink` and
  :class:`~orpheus.transport.source_sinks.AngularBoundarySourceSink` pre-date
  this work — they are the *source* column of the role grid.
* The affine boundary law :eq:`affine-bc-form` and the ``q.boundary``
  slot pre-date this work.

What was missing was the *ergonomic generator* for a non-vacuum
boundary leaf and a *public entry point* that accepts the composite.
:meth:`~orpheus.transport.source_sinks.AngularBoundarySourceSink.prescribed_inflow`
and the second accepted form of
:func:`~orpheus.sn.solver.solve_sn_fixed_source` supply exactly those —
no more. This is the operational meaning of "we already have the right
objects/concepts — we just need better ergonomics to generate them":
the abstraction is unchanged; only the surface that builds it is
better. The first consumer is the non-vacuum prescribed-inflow MMS of
:ref:`sn-mms-nonvacuum`.


.. _sn-mms-nonvacuum:

Non-vacuum prescribed-inflow MMS (Phase 4 / O.2b 4.6)
------------------------------------------------------

.. admonition:: Key Facts
   :class: important

   - **What this section adds.** The entire pre-existing MMS catalog
     (:ref:`sn-mms-curvilinear-isotropic-verification`,
     :ref:`sn-mms-curvilinear-aniso-verification`) is
     *vacuum-automatic*: every manufactured ansatz vanishes at both
     boundaries, so the inflow trace :math:`\gamma_-\psi \equiv 0` on
     every ordinate and the prescribed-inflow source slot
     ``q.boundary`` is **identically zero** in all of them. Phase 4 /
     O.2b sub-step 4.6 fills that gap with a manufactured solution that
     is **non-zero at the outer face**, lighting the
     :math:`q.\text{boundary} \neq 0` path for the first time.
   - **The ansatz is the proven P1 element.** :math:`\psi_n = (A + \mu_n
     B)/W` is the same truncated-Legendre :math:`P_0 \oplus P_1` form
     used by the Phase 3.6 anisotropic cases
     (:eq:`sn-mms-spherical-aniso-psi`). 4.6 changes **only the
     boundary trace** — :math:`A,B` are chosen non-vanishing at the
     outer face — and reuses the verified angular structure. Linear in
     :math:`\mu` *fully* (not partially) activates the curvilinear
     redistribution; the question "do we need :math:`\mu^2`?" is
     answered **no** below.
   - **Two manufactured sources, derived from the continuous
     operator.** The slab source :eq:`sn-mms-nonvacuum-qext` has **no**
     redistribution term (the Cartesian operator lacks the
     :math:`\partial_\mu` coupling); the sphere source
     :eq:`sn-mms-nonvacuum-sph-qext` is the **same closed form** as the
     Phase 3.6 vacuum case — only :math:`A,B` differ. The spherical
     residual therefore lives in **one place**
     (:func:`~orpheus.derivations.continuous.mms.sn._spherical_anisotropic_symbolic`,
     Cardinal Rule 2).
   - **HAZARD H1 (sphere pole regularity).** :math:`B(0)=0` is a HARD
     constraint — without it the redistribution :math:`(1-\mu^2)B/r \to
     \infty` at :math:`r=0`. The :math:`(r/R)` prefactor on the sphere
     :math:`B` enforces it; the slab has no pole, so a slab-style
     :math:`B(0)\neq 0` is fine there but **wrong** on the sphere.
   - **The affine-BC-to-RHS framing.** Prescribed inflow IS the
     inhomogeneous term :math:`q` of the affine boundary law
     :eq:`affine-bc-form` :math:`\gamma_-\psi = R\,G\,\gamma_+\psi + q`,
     carried in the ``q.boundary`` slot
     (:class:`~orpheus.transport.source_sinks.AngularBoundarySourceSink`) and
     consumed directly by :math:`(L+C)\text{.solve}` as the sweep
     inflow seed. **No** :class:`~orpheus.geometry.boundary.PrescribedInflow`
     mesh-BC bridge is touched, and **no** ``from_spec`` recipe bridge
     is needed — the inflow is supplied as the boundary leaf of the
     **composite source** ``q = q_bulk ⊕ q_∂`` that
     :func:`~orpheus.sn.solver.solve_sn_fixed_source` now accepts
     directly (see :ref:`sn-composite-fixed-source`).
   - **The convergence rows drive the public composite-source API.**
     The 4.6 MMS no longer assembles the within-group operator triple
     by hand: :func:`~orpheus.sn.solver.solve_sn_fixed_source` accepts
     a :class:`~orpheus.transport.timed_full_field.TimedFullField`
     composite source, so each case bundles its manufactured bulk and
     prescribed-inflow boundary into one ``case.fixed_source(sn)``
     call. The migration off the operator-triple bypass *is* the
     retirement (retirement = test migration).
   - **The load-bearing assertion is the converged VALUE, not the
     rate.** Per the ``vv-principles`` skill anti-pattern #5 (rate is
     necessary, not sufficient), a silently-dropped ``q.boundary``
     converges cleanly at :math:`\mathcal{O}(h^2)` to the **wrong**,
     boundary-zero limit. Only the flux-value-vs-:math:`A(x)` check —
     with :math:`A` non-zero at the boundary (:math:`a_0>0`) — sees it.
   - **T3 (sphere) ships ``xfail(strict)``, re-scoped to Issue #229.**
     The slab rows are clean :math:`\mathcal{O}(h^2)` with value match.
     The sphere row's anisotropic :math:`(A+B\mu)/W` ansatz is
     angle-varying, so after ERR-058 (#195) recovered the spatial
     :math:`\mathcal{O}(h^2)` convergence it now hits the
     fixed-quadrature **angular floor** of the half-angle thread
     interpolation (sphere S16 floor :math:`\approx 7\mathrm{e}{-4}`,
     above the band) — NOT the old #195 plateau, and NOT a non-vacuum
     machinery failure (the boundary value *is* honoured).  The marker
     moved from #195 to the
     `Issue #229 <https://github.com/deOliveira-R/ORPHEUS/issues/229>`_
     quadrature-aware retune.  Both the sphere and slab rows now run the
     curvilinear/Cartesian **source-iteration** default of
     :func:`~orpheus.sn.solver.solve_sn_fixed_source` (SI :math:`\equiv`
     Krylov bit-identical post-ERR-058); the composite-source API
     delivers the prescribed inflow identically to every splitting (T4).
     The green companion T3g provides live structural coverage of the
     inflow + redistribution paths now.

This section narrates the Branch-1 SymPy algebra-of-record
(:mod:`orpheus.derivations.continuous.mms.sn`), the Branch-2 numpy
factories, and the L1 / foundation gates that verify the prescribed-
inflow discretisation. The verification chain follows the
``algebra-of-record`` discipline (Branch-1 SymPy reference, Branch-2
numpy production, structurally-independent L1 cross-check).

.. list-table:: 4.6 verification gates (measured, ``-O`` mode)
   :header-rows: 1
   :widths: 6 30 10 12 32

   * - Gate
     - Description
     - Level
     - Pillar
     - Status / evidence
   * - V_nonvac-slab
     - Slab substitution identity ``simplify(W·LHS − Σ_s φ − Q) == 0``
     - foundation
     - MMS (1C)
     - PASS — :func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_v_nonvac_slab_substitution_identity`
   * - V_nonvac-sph
     - Sphere substitution identity (reuses the 3.6 spherical residual)
     - foundation
     - MMS (1C)
     - PASS — :func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_v_nonvac_sph_substitution_identity`
   * - Decision-A pin
     - Parameterised :math:`A=B=` ``None`` reproduces 3.6 vacuum shapes byte-for-byte
     - foundation
     - regression
     - PASS — :func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_existing_spherical_aniso_still_passes_after_parameterization`
   * - L1 xcheck (slab)
     - Branch-2 numpy :math:`Q^{\text{ext}}` == lambdified SymPy (≤1e-13)
     - foundation
     - MMS (1C)
     - PASS — :func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_slab_nonvacuum_numerical_qext_matches_sympy`
   * - L1 xcheck (sphere)
     - Branch-2 numpy :math:`Q^{\text{ext}}` == lambdified SymPy (≤1e-13)
     - foundation
     - MMS (1C)
     - PASS — :func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_sphere_nonvacuum_numerical_qext_matches_sympy`
   * - T1 (slab 1g)
     - DD :math:`\mathcal{O}(h^2)` + converged value + inflow honoured
     - L1
     - MMS (1C)
     - PASS — orders ``[2.04, 2.01]``, finest L2 ~1.2e-3, max\|φ−A\| ~8e-5
   * - T2 (slab 2g asym)
     - As T1, asymmetric downscatter :math:`\Sigma_s` (ERR-002 hazard)
     - L1
     - MMS (1C)
     - PASS — g0 ``[2.04, 2.01]``, g1 ``[2.05, 2.01]``
   * - T3 (sphere)
     - Curvilinear redistribution under non-vacuum inflow
     - L1
     - MMS (1C)
     - **xfail(strict)** on #229 — aniso angular floor ≈7e-4 (spatial
       O(h²) recovered by ERR-058; boundary value honoured)
   * - T3g (sphere)
     - Inflow honoured at :math:`r=R` + redistribution source live (green now)
     - foundation
     - structural
     - PASS — :func:`tests.sn.verification.analytical.test_mms_prescribed_inflow.test_sphere_nonvacuum_inflow_honoured_and_redistribution_live`
   * - T4 (Mode 9)
     - SI-Jacobi ≡ SI-Gauss-Seidel ≡ Krylov honour ``q.boundary``
     - foundation
     - splitting-invariance
     - PASS — pairwise reldiffs 1.3e-13 … 5.6e-13

Why the existing MMS catalog is vacuum-automatic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every manufactured solution already in the SN verification ladder —
the isotropic :math:`\psi_n = A(r)/W` cases
(:ref:`sn-mms-curvilinear-isotropic-verification`) and the
anisotropic :math:`\psi_n = (A + \mu_n B)/W` cases
(:ref:`sn-mms-curvilinear-aniso-verification`) — was built with
:math:`A` and :math:`B` chosen to **vanish at both boundaries**. For
the canonical 3.6 sphere, :math:`A(r) = \sin(\pi r/R)` gives
:math:`A(0) = A(R) = 0`, and :math:`B(r) = (r/R)(1-r/R)\cos(\pi r/R)`
gives :math:`B(0) = B(R) = 0`. The slab isotropic case likewise uses
:math:`A(x) = \sin(\pi x/L)`.

The consequence is structural and total. The inflow trace of the
manufactured solution on any face is

.. math::

   \gamma_- \psi_n \big|_{\text{face}}
       = \frac{1}{W}\bigl(A(x_{\text{face}})
                         + \mu_n B(x_{\text{face}})\bigr)
       = \frac{1}{W}\bigl(0 + \mu_n \cdot 0\bigr) = 0
       \qquad \text{for every ordinate } n.

So the affine boundary law :eq:`affine-bc-form`
:math:`\gamma_-\psi = R\,G\,\gamma_+\psi + q` collapses to its
homogeneous (vacuum) form for these cases — and the inhomogeneous
inflow term :math:`q` is identically zero. The existing cases verify
the **interior** spatial / angular operator and the **homogeneous**
vacuum BC, but they say *nothing* about the prescribed-inflow path,
where a non-zero :math:`q` is pushed into the right-hand side. That
path is the one O.2b's field-role-typing work makes a first-class
boundary trace: an inhomogeneous inflow injected as the boundary
*source* slot. Until 4.6, no MMS row exercised it.

The fix is the smallest possible structural delta: keep the proven P1
angular form, keep the proven interior operator, and change **only**
:math:`A,B` so that they are non-zero at the outer face. Then
:math:`\gamma_-\psi \neq 0`, and the converged scalar flux
:math:`\phi(x) = A(x)` is non-zero at the boundary — exactly the
property the verification needs (see "The converged-value assertion"
below).

The ansatz — the P1 element and why linear-in-:math:`\mu` is enough
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The slab ansatz and its manufactured source:

.. math::
   :label: sn-mms-nonvacuum-psi

   \psi_n(x) = \frac{1}{W}\bigl(A(x) + \mu_n B(x)\bigr),
   \qquad A(x) = a_0 + a_1\sin(kx),\quad B(x) = b_0\cos(kx),
   \quad a_0 > 0.

.. (vv-status rationale) definition: Definitional ansatz — the
.. manufactured angular flux is *imposed*, not solved for. Its
.. correctness as a reference is established by the source identity
.. :eq:`sn-mms-nonvacuum-qext` (SymPy ``simplify == 0``), not by a
.. property of this expression alone.
.. vv-status: sn-mms-nonvacuum-psi documented

The form :math:`\psi_n = (A + \mu_n B)/W` is the truncated Legendre
:math:`P_0 \oplus P_1` element: :math:`P_0(\mu) = 1` carries the
isotropic amplitude :math:`A`, and :math:`P_1(\mu) = \mu` carries the
first-moment amplitude :math:`B`. This is the **native angular basis**
of the SN closure — the Carlson half-angle pole closure folds the
moment source through :math:`P_\ell(-1) = (-1)^\ell`, so a linear-in-
:math:`\mu` input is exactly the lowest non-trivial moment with a
non-zero :math:`\partial_\mu`.

**Why linear-in-:math:`\mu` fully activates the redistribution.** The
curvilinear angular-redistribution operator is
:math:`\tfrac{1-\mu^2}{r}\,\partial_\mu\psi`. With :math:`\psi` linear
in :math:`\mu`, the angular derivative is a non-zero **constant** in
:math:`\mu`,

.. math::

   \frac{\partial \psi_n}{\partial \mu} = \frac{B(r)}{W} \neq 0,

and multiplying by the redistribution dome :math:`(1-\mu^2)` produces
a genuinely :math:`\mu^2`-structured term :math:`(1-\mu^2)B/r`. The
discrete closure that realises this operator — the Morel–Montry
half-angle recurrence with the Carlson :math:`\mu=-1` seed (see
:ref:`sn-pole-angular-closure-protocol`) — is **linear** in
:math:`\psi`. A linear operator is fully probed by any input that is
non-constant in its argument; the linear-in-:math:`\mu` ansatz is
non-constant in :math:`\mu`, so it exercises the entire linear
redistribution map (including the half-angle recurrence and the
second-moment coupling).

A quadratic-in-:math:`\mu` (P2) ansatz term would add **no** new
structural coverage of the redistribution. Because the closure is
linear, a quadratic input only changes *which point* in the operator's
already-fully-probed range you land on — it does not reach any term
the linear input misses. (A P2 term *would* additionally exercise the
:math:`\sum_n w_n \mu_n^2` quadrature-exactness, but that is a
property of the quadrature, not of the redistribution operator, and is
already covered elsewhere.) This settles the "do we need
:math:`\mu^2`?" question definitively: **no.** The verdict is recorded
in the cross-domain-attacker frame analysis (memo
``phase4_o2b_4_6_mms_ansatz_frame.md``, Q1/Q2) and is empirically
consistent with Phase 3.6, which uses exactly this linear-in-:math:`\mu`
ansatz and whose gate tests carry ``catches("ERR-026")`` — the
redistribution-bug catcher.

**The scalar flux is :math:`A`.** Because Gauss–Legendre (and every
symmetric quadrature on :math:`\mu \in [-1,1]`) satisfies
:math:`\sum_n w_n \mu_n = 0`, the first-moment term integrates out of
the scalar moment:

.. math::

   \phi(x) = \frac{1}{W}\sum_n w_n\bigl(A(x) + \mu_n B(x)\bigr)
           = \frac{1}{W}\Bigl(A(x)\sum_n w_n
                            + B(x)\underbrace{\sum_n w_n \mu_n}_{=\,0}\Bigr)
           = A(x).

This discrete identity is verified directly in
:func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_slab_nonvacuum_phi_equals_A_under_quadrature`
(≤1e-14 on a sample mesh). The reference scalar flux for the
convergence rows is therefore :math:`\phi_{\text{chosen}}(x) = A(x)`,
which — because :math:`a_0>0` — is **non-zero at the boundary**.

The manufactured slab source, derived from the continuous operator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The slab SN transport operator (per ordinate, 1-group) is the
first-order streaming-plus-collision form
:math:`\mu\,\partial_x\psi + \Sigma_t\psi
= \tfrac{1}{W}(\Sigma_s\phi + Q^{\text{ext}})`. The Cartesian operator
has **no** angular-derivative term — the slab geometry produces no
angular redistribution (:eq:`transport-cartesian`). Substituting the
ansatz :eq:`sn-mms-nonvacuum-psi` and solving for the residual source:

.. math::

   \mu\,\frac{\partial}{\partial x}
       \frac{A + \mu B}{W}
   + \Sigma_t\,\frac{A + \mu B}{W}
   &= \frac{1}{W}\bigl(\Sigma_s\,A + Q^{\text{ext}}_n\bigr) \\[1mm]
   \frac{1}{W}\bigl(\mu A' + \mu^2 B'\bigr)
   + \frac{\Sigma_t}{W}\bigl(A + \mu B\bigr)
   &= \frac{1}{W}\bigl(\Sigma_s A + Q^{\text{ext}}_n\bigr),

where :math:`\phi = A` was used on the right. Multiplying through by
:math:`W` and isolating :math:`Q^{\text{ext}}_n` gives the closed form

.. math::
   :label: sn-mms-nonvacuum-qext

   Q^{\text{ext}}_n(x) = \mu_n A'(x) + \mu_n^2 B'(x)
                       + (\Sigma_t - \Sigma_s) A(x)
                       + \Sigma_t\,\mu_n B(x).

.. (vv-status rationale) derivation: A closed form obtained by
.. symbolic substitution into the continuous slab operator; verified
.. by SymPy ``simplify(W·LHS − Σ_s φ − Q) == 0`` and cross-checked
.. against the Branch-2 numpy evaluation to ≤1e-13.
.. vv-status: sn-mms-nonvacuum-qext verified

Note that there is **no** :math:`(1-\mu^2)B/r` term — the slab operator
simply does not generate it. The :math:`\mu^2 B'` term *is* present
(streaming the first-moment piece :math:`\mu B` gives
:math:`\mu \cdot \mu B' = \mu^2 B'`), so the slab still exercises the
second-moment streaming closure, just not the angular redistribution.
The :math:`(\Sigma_t - \Sigma_s)A` term is the within-group removal
net of isotropic self-scatter (:math:`c = \Sigma_s/\Sigma_t < 1`), and
:math:`\Sigma_t\,\mu_n B` is the collision of the first-moment piece.

The Branch-1 algebra-of-record is
:func:`~orpheus.derivations.continuous.mms.sn.derive_nonvacuum_slab_mms`
(building on
:func:`~orpheus.derivations.continuous.mms.sn._nonvacuum_slab_symbolic`),
which performs the substitution symbolically and proves
``simplify(W·LHS − Σ_s·φ − Q_closed) == 0``. Because the slab operator
lacks redistribution, it is a *genuinely different* operator from the
sphere and gets its own fresh symbolic pair — it cannot reuse the
spherical residual (which carries the :math:`\partial_\mu` term the
slab does not have).

**Multi-group generalisation (T2).** The slab case is multi-group-
capable. Each group carries a per-group amplitude :math:`c_g` scaling
the shared shape, :math:`A_g(x) = c_g(a_0 + a_1\sin kx)` and
:math:`B_g(x) = c_g\,b_0\cos kx`, and the source picks up the
in-scatter term

.. math::

   Q^{\text{ext}}_{n,g}(x) = \mu_n A_g'(x) + \mu_n^2 B_g'(x)
       + \Sigma_{t,g}\,A_g(x) + \Sigma_{t,g}\,\mu_n B_g(x)
       - \sum_{g'} \Sigma_s[g', g]\,A_{g'}(x).

The in-scatter sum uses the ORPHEUS scattering convention
``SigS[g_from, g_to]``, so the in-scatter source is
:math:`(\Sigma_s^\top\phi)_g = \sum_{g'}\Sigma_s[g', g]\,A_{g'}` — the
**transpose-active** term where the ERR-002 group-swap hazard lives.
T2 uses a 2-group **asymmetric downscatter-only** :math:`\Sigma_s`
(:math:`\Sigma_s[0,1]\neq 0`, :math:`\Sigma_s[1,0]=0`) so a transposed
scattering matrix would produce a detectably wrong group ratio
(the 1-group-degeneracy rule — multi-group with asymmetric :math:`\Sigma_s` is
mandatory, ``vv-principles`` anti-pattern #3 and failure-mode #6). The
1-group T1 path is the degenerate :math:`c_{\text{groups}} = (1.0,)`
reduction of the same dataclass.

The manufactured spherical source — the Cardinal-Rule-2 reuse
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The spherical ansatz (pole-regular, non-vacuum at :math:`r=R`):

.. math::
   :label: sn-mms-nonvacuum-sph-psi

   \psi_n(r) = \frac{1}{W}\bigl(A(r) + \mu_n B(r)\bigr),
   \quad A(r) = a_0 + a_1\sin(kr),\quad
   B(r) = \frac{r}{R}\bigl[b_0 + b_1\cos(kr)\bigr],\quad B(0)=0.

.. (vv-status rationale) definition: Definitional ansatz — imposed,
.. not solved. Correctness as a reference rests on the source
.. identity :eq:`sn-mms-nonvacuum-sph-qext` (SymPy ``simplify == 0``)
.. and HAZARD H1 (:math:`B(0)=0`), verified in the foundation gate.
.. vv-status: sn-mms-nonvacuum-sph-psi documented

The spherical SN operator carries the angular-redistribution term
(:eq:`transport-spherical`):
:math:`\mu\,\partial_r\psi + \tfrac{1-\mu^2}{r}\,\partial_\mu\psi
+ \Sigma_t\psi = \tfrac{1}{W}(\Sigma_s\phi + Q^{\text{ext}})`.
Substituting :eq:`sn-mms-nonvacuum-sph-psi` (with
:math:`\partial_\mu\psi = B/W`, so the redistribution term is
:math:`\tfrac{1-\mu^2}{r}\cdot\tfrac{B}{W}`) and isolating the source
gives the **same closed form** as the Phase 3.6 vacuum case
(:eq:`sn-mms-spherical-aniso-qext`):

.. math::
   :label: sn-mms-nonvacuum-sph-qext

   Q^{\text{ext}}_n(r) = \mu_n A'(r) + \mu_n^2 B'(r)
                       + (1-\mu_n^2)\,\frac{B(r)}{r}
                       + (\Sigma_t-\Sigma_s) A(r)
                       + \Sigma_t\,\mu_n B(r).

.. (vv-status rationale) derivation: A closed form obtained by
.. symbolic substitution into the continuous spherical operator;
.. verified by SymPy ``simplify == 0`` (reusing the 3.6 spherical
.. residual machinery) and cross-checked against the Branch-2 numpy
.. evaluation to ≤1e-13.
.. vv-status: sn-mms-nonvacuum-sph-qext verified

The structural point is that :eq:`sn-mms-nonvacuum-sph-qext` and the
Phase 3.6 :eq:`sn-mms-spherical-aniso-qext` are *byte-identical* closed
forms — only the radial profiles :math:`A,B` plugged into them differ.
The spherical-operator residual is therefore derived in **exactly one
place**:
:func:`~orpheus.derivations.continuous.mms.sn._spherical_anisotropic_symbolic`,
which now takes optional ``A=None, B=None`` arguments. With no
arguments it reproduces the Phase 3.6 vacuum shapes byte-for-byte (the
decision-A regression pin verifies this in
:func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_existing_spherical_aniso_still_passes_after_parameterization`);
with the 4.6 non-vacuum shapes it re-proves the residual for free
(:func:`~orpheus.derivations.continuous.mms.sn.derive_nonvacuum_spherical_mms`).
This is Cardinal Rule 2 in action — one source of truth for the
spherical transport-operator residual, shared between the vacuum and
non-vacuum cases.

HAZARD H1 — sphere pole regularity demands :math:`B(0)=0`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The redistribution term :math:`(1-\mu^2)B(r)/r` has an explicit
:math:`1/r` factor. As :math:`r\to 0`, this diverges as
:math:`(1-\mu^2)\,B(0)/r \to \infty` **unless** :math:`B(0)=0`. So on
the sphere, :math:`B(0)=0` is a **hard regularity constraint** — not a
stylistic preference. A naive slab-style choice :math:`B(x) = b_0\cos
kx` gives :math:`B(0) = b_0 \neq 0`, which is **fine on the slab** (no
pole, no :math:`1/r`) but **wrong on the sphere** (it manufactures a
non-integrable :math:`1/r` singularity at the centre that the
continuous solution does not actually have).

The 4.6 sphere therefore uses :math:`B(r) = (r/R)[b_0 + b_1\cos kr]`.
The :math:`(r/R)` prefactor forces :math:`B(0) = 0` (pole-regular: the
redistribution :math:`(1-\mu^2)B/r = (1-\mu^2)[b_0+b_1\cos kr]/R` is
*finite* at :math:`r=0`), while leaving :math:`B(R) = b_0 + b_1\cos kR
\neq 0` (the non-vacuum first-moment structure at the outer inflow
face). The amplitude :math:`A(r) = a_0 + a_1\sin kr` needs **no** such
prefactor: :math:`A` has no :math:`1/r` companion in the operator, so
:math:`A(0) = a_0` finite is perfectly regular at the pole, and
:math:`a_0>0` makes :math:`A(R)\neq 0` (non-vacuum). HAZARD H1 is
verified in
:func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_v_nonvac_sph_pole_regularity_and_nonvacuum`
(:math:`B(0)=0`, :math:`A(0)=\tfrac12`, :math:`B(R)\neq 0`; concretely
at :math:`kR=\pi/2`: :math:`A(R)=\tfrac34`, :math:`B(R)=\tfrac{3}{10}`).

The non-vacuum lever — :math:`a_0>0` is the entire novelty
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every other choice in 4.6 (the P1 angular form, the interior operator,
the redistribution term, the quadrature, the BC machinery) is shared
with Phase 3.6. The *single* new ingredient is :math:`a_0 > 0`, which
makes :math:`A` — and hence the inflow trace
:math:`\gamma_-\psi_n = (A + \mu_n B)/W` — **non-zero at the outer
face**. That non-zero trace is what lights up the prescribed-inflow
``q.boundary`` path. Strip :math:`a_0` back to zero and 4.6 degenerates
to Phase 3.6 (vacuum-automatic). The non-vacuum-ness is pinned by the
foundation test
:func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_v_nonvac_slab_ansatz_nonvanishing_at_faces`
(:math:`A(0)=a_0>0`) so the verification cannot silently drift back to
the vacuum regime.

The affine-BC-to-RHS framing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Prescribed inflow is **not** a special solver mode — it is the
inhomogeneous term of the universal affine boundary law
(:ref:`affine-bc-form`). The general boundary trace law is

.. math::

   \gamma_-\psi = R\,G\,\gamma_+\psi + q,

where :math:`q \in \Gamma_-` is the **prescribed inflow source**. For a
vacuum boundary :math:`R=0` and :math:`q\equiv 0`. For a manufactured
non-vacuum inflow, :math:`q = \gamma_-\psi_{\text{chosen}}` — the
imposed inflow trace — pushed to the right-hand side of the
discretised within-group system.

In ORPHEUS this :math:`q` is carried by the ``q.boundary`` slot, a
:class:`~orpheus.transport.source_sinks.AngularBoundarySourceSink` field whose
inflow-ordinate entries hold :math:`\gamma_-\psi = (A + \mu_n B)/W` per
face per group. The within-group fixed point is the **affine** system

.. math::

   (L + C - S - B)\,\psi = q,
   \qquad q = q_{\text{ext}}
            + (\text{prescribed inflow in } q.\text{boundary}),

and the inflow term is consumed directly by :math:`(L+C)\text{.solve}`
as the sweep inflow seed. This is the cleanest possible realisation:
the inhomogeneous BC term is *just another source* on the RHS.

**No ``from_spec`` / ``PrescribedInflow``-BC bridge is touched.** A
:class:`~orpheus.geometry.boundary.PrescribedInflow` mesh-BC descriptor
*does* exist (the rank-0 affine BC), but it is a *different surface* —
it declares a prescribed inflow at mesh-construction time as a
first-class boundary condition. The 4.6 MMS deliberately does **not**
use it: the manufactured inflow is injected as the ``q.boundary``
source slot, which is exactly the affine-:math:`q`-to-RHS path, and is
the surface O.2b's field-role-typing work targets. The mesh BCs for the
4.6 cases are plain **vacuum** — the inflow lives entirely in
``q.boundary``.

The public composite-source API drives the convergence rows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The 4.6 convergence rows drive the **public** fixed-source entry point
:func:`~orpheus.sn.solver.solve_sn_fixed_source` directly. Earlier in
this work that entry point hardcoded a vacuum ``q.boundary``
(``zeros_on``) — it had no way to carry a prescribed inflow — and the
rows therefore took an operator-triple *bypass*: assembling
:math:`(L+C)`, :math:`S`, :math:`B` by hand via
:func:`~orpheus.sn.coupled_system.build_within_group_system` and driving
them with :func:`~orpheus.numerics.iteration.SourceIteration`. That bypass
is **retired**. The field-role-typing work gave
:func:`~orpheus.sn.solver.solve_sn_fixed_source` a second accepted
source form — the full **composite source** ``q = q_bulk ⊕ q_∂``
represented by a
:class:`~orpheus.transport.timed_full_field.TimedFullField` (see
:ref:`sn-composite-fixed-source` for the API in full). Each case now
bundles its manufactured bulk (:meth:`external_source`) and its
prescribed-inflow boundary (:meth:`prescribed_inflow`) into one
``case.fixed_source(sn)`` and passes it straight to the public solver::

    result = solve_sn_fixed_source(
        materials, mesh, case.quadrature, case.fixed_source(sn),
        max_inner=1000, inner_tol=1e-13,
    )

Migrating the rows off the bypass onto the public API *is* the
retirement (retirement = test migration — the new code is what gets
tested). The slab rows take the SI (1-D Jacobi) inner; the sphere row
takes the curvilinear **source-iteration** default (post-ERR-058, #195;
SI :math:`\equiv` Krylov bit-identical — the curvilinear ``"krylov"``
default was reverted, see :ref:`sn-err-058-closure-seed-closeout`); both
honour the prescribed inflow identically (verified by T4 below).

**The flux/source space bridge — now INTERNAL to the solve (B.5.2).**
The composite RHS lives in **source** space (an
:class:`~orpheus.transport.source_sinks.AngularSourceSink` bulk plus a
:class:`~orpheus.transport.source_sinks.AngularBoundarySourceSink` boundary),
while the iterate :math:`\psi` and the returned solution live in
**flux** space (an :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
bulk plus a :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
boundary). The source-iteration / Krylov inner therefore needs a
flux-typed ``initial_guess`` to template the solution space — without
it, ``S.apply`` would hit an ``AngularSourceSink`` that has no
``integrate_angular`` method. That seed
(``TimedFullField.zeros(bulk=AngularFlux, boundary=AngularBoundaryFlux,
mesh=sn)``) is now built **inside**
:func:`~orpheus.sn.solver.solve_sn_fixed_source`, not hand-passed by
the test. The field-role-typing distinction — the iterate is a *flux*,
the RHS is a *source* — survives intact; it has simply moved behind the
public API where it belongs.

The converged-value assertion — rate is necessary, not sufficient
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the load-bearing verification design choice, and it is a
direct application of ``vv-principles`` anti-pattern #5 ("NEVER read
'convergence rate is correct' as 'result is correct' — verify the
converged-to value; :math:`\mathcal{O}(h^2)` to the wrong limit is
still :math:`\mathcal{O}(h^2)`") and the necessity hierarchy H4.

Consider the failure mode the test must catch: a bug (or a refactor
regression) that silently **drops the prescribed inflow** — a solve
that runs with ``q.boundary = 0`` despite the manufactured non-vacuum
inflow. That degenerate solve is **still a perfectly consistent
fixed-source problem** — it just solves the *vacuum-BC* version of the
same interior source. It converges cleanly at :math:`\mathcal{O}(h^2)`
to a *different*, boundary-zero scalar flux. A rate-only test passes
it. The only assertion that sees the dropped inflow is the one that
checks the **converged value against** :math:`A(x)`, because
:math:`A(x)` is non-zero at the boundary (:math:`a_0>0`) while the
dropped-inflow limit is zero there — a discrepancy of order
:math:`a_0 \approx 0.5` at the faces, dwarfing the pointwise
convergence error (~8e-5).

The slab T1/T2 rows therefore make **three** assertions per group:
(1) the rate ``orders > 1.9`` (DD design order on a smooth ansatz);
(2) the finest-mesh :math:`\phi_{\text{num}}` matches :math:`A(x)` to
``rtol=atol=5e-3`` — with a guard asserting the reference is genuinely
non-vacuum (:math:`|A(0)|, |A(L)| > 0.1`) so the value check is
discriminating; and (3) an inflow-honoured spot-check that the solved
trace slot equals the imposed :math:`\gamma_-\psi = (A + \mu_n B)/W` to
``rtol=1e-9``. Only the combination is a meaningful test of the
prescribed-inflow path.

The Mode-7 activates/nulls map
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``vv-principles`` failure-mode #7 ("MMS simplification bias") requires
every multi-dimensional MMS test to **declare** which operator terms
its ansatz activates and which it nulls — and to ship an
angularly-non-trivial companion whenever the nulled set includes a term
covered by an active ERR-NNN. The 4.6 declaration:

.. list-table:: Mode-7 term map — slab vs sphere under the (A+μB)/W ansatz
   :header-rows: 1
   :widths: 40 30 30

   * - Operator term
     - Slab (Cartesian)
     - Sphere (spherical)
   * - streaming :math:`\mu A'` (isotropic)
     - **activates**
     - **activates**
   * - streaming :math:`\mu^2 B'` (second moment)
     - **activates**
     - **activates**
   * - angular redistribution :math:`(1-\mu^2)B/r`
     - **nulls** (no :math:`\partial_\mu` term)
     - **activates** (the ERR-026 path)
   * - within-group scatter :math:`\Sigma_s\phi` (:math:`c<1`)
     - **activates**
     - **activates**
   * - collision :math:`\Sigma_t\,\mu B` (first moment)
     - **activates**
     - **activates**
   * - 2G group transfer :math:`\Sigma_s^\top` (asymmetric)
     - **activates** (T2)
     - n/a (1g)
   * - prescribed non-vacuum inflow :math:`\gamma_-\psi \neq 0`
     - **activates** (both faces, :math:`a_0>0`)
     - **activates** (:math:`r=R` face)
   * - fission
     - **nulls** (non-fissile; MMS proves no eigenvalue)
     - **nulls**

The slab **nulls the angular redistribution** — the Cartesian operator
has no :math:`\partial_\mu` coupling. Redistribution is exactly where
ERR-026 (the curvilinear sweep WDD wrong-fixed-point bug) lives, so a
slab-only 4.6 would be a textbook Mode-7 trap: it would verify the
prescribed-inflow path while being structurally blind to the hardest
math the curvilinear sweep performs. The **sphere companion is
therefore mandatory** (NEVER ship slab-only — ERR-026 territory). The
sphere activates the redistribution term under non-vacuum inflow,
closing the Mode-7 declaration.

T3 (sphere) — why ``xfail(strict)``, now on Issue #229
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The sphere row
(:func:`tests.sn.verification.analytical.test_mms_prescribed_inflow.test_mms_prescribed_inflow_sphere_activates_redistribution`)
ships ``@pytest.mark.xfail(strict=True)`` with ``catches("ERR-026")``.
The reason is **not** that the non-vacuum machinery fails.

.. note:: **Re-scoped (2026-06-12, Issue #195 CLOSED → Issue #229).**

   This row's xfail was originally attributed to the #195
   "pre-asymptotic transient" plateau (the stagnation table below).
   The ERR-058 closure-seed fix **closed** the curvilinear
   wrong-fixed-point family, so the stagnation is gone — the isotropic
   curvilinear DD interior is now :math:`\mathcal{O}(h^2)`-consistent
   (see :ref:`sn-err-058-closure-seed-closeout`).  T3's ansatz, however,
   is the *anisotropic* :math:`(A(r)+B(r)\mu)/W`, which is angle-varying
   and therefore hits the **fixed-quadrature angular floor** of the
   half-angle thread interpolation (sphere S16 floor
   :math:`\approx 7\mathrm{e}{-4}`, above the band).  The marker
   **stays ``xfail(strict)``** but is now pinned to the
   `Issue #229 <https://github.com/deOliveira-R/ORPHEUS/issues/229>`_
   quadrature-aware retune, NOT the #195 plateau; it flips to
   unexpected-pass when #229 lands (the regression gate for the retune).
   The stagnation table below is preserved as **bug-era evidence**; its
   "pre-asymptotic" interpretation is superseded.

**Bug-era stagnation (pre-ERR-058).**  The slab was pole-free and
converged perfectly: orders ``[2.04, 2.01]``, finest L2 ~1.2e-3,
pointwise ``max|φ−A|`` ~8e-5, boundary value matched.  The sphere L2
(volume-weighted), by contrast, **stagnated** mesh-independently — the
plateau that refuted the "pre-asymptotic transient" premise:

.. list-table:: T3 sphere volume-weighted L2 error (bug-era plateau, pre-ERR-058)
   :header-rows: 1
   :widths: 25 25 25 25

   * - :math:`n_c`
     - 20
     - 40
     - 80
   * - :math:`\|\phi_h - A\|_{L^2(V)}`
     - 2.37e-2
     - 2.42e-2
     - 2.43e-2

The observed "orders" were ≈ :math:`-0.02` to :math:`-0.006` — the
error was *not* decreasing under refinement, the plateau signature
ERR-058 diagnosed.  Post-ERR-058 the spatial convergence is recovered;
the residual gap on this *anisotropic* row is the #229 angular floor,
which DOES drop under quadrature refinement (sphere S16
:math:`\to` S32 halves the floor) — the structural test that the
remaining gap is angular, not a wrong fixed point.

Crucially, the **boundary value is honoured** (always was): the
finest-mesh :math:`\phi[-1] \approx 0.7499 \approx A(R) = 0.75`, and
the inflow-trace spot check passes.  The non-vacuum prescribed-inflow
machinery *works*; the remaining xfail is purely the angular-floor
budget of the anisotropic ansatz under fixed quadrature.

**T3g — the green structural companion.** Because T3 is xfail (now on
#229), it provides *no live* convergence coverage of the 4.6 machinery.
The green companion
:func:`tests.sn.verification.analytical.test_mms_prescribed_inflow.test_sphere_nonvacuum_inflow_honoured_and_redistribution_live`
fills that gap with two non-convergence-dependent claims that pass
*now*: (1) the prescribed inflow at :math:`r=R` is honoured per inflow
ordinate (:math:`\gamma_-\psi = (A(R) + \mu_n B(R))/W` with :math:`A(R)
> 0` non-vacuum); and (2) the redistribution source
:math:`(1-\mu^2)B(r)/r` is non-zero on the mesh interior (the ERR-026
term is live under the 4.6 ansatz — :math:`B(r)\neq 0` on the open
interval, with :math:`B(0)=0` pole-regular). T3g is the live structural
guarantee that the Mode-7 sphere companion exercises the redistribution
path even while the convergence row is parked on #229.

T4 (vv Mode 9) — splitting invariance of the prescribed inflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The consistency floor the convergence rows trust is that a non-zero
prescribed inflow is honoured **identically** by the three operator
splittings of the affine within-group system: SI-Jacobi (the resolvent
:math:`L+C` with lagged gains :math:`S, B`), SI-Gauss–Seidel (the
:math:`B_{\rm lower}`-folding reified :math:`M`,
:class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`,
with lagged gains :math:`S, B_{\rm upper}`), and Krylov (the matvec
:math:`L+C-S-B`). All three are different reduction trees of the *same*
affine fixed point :math:`(L+C-S-B)\psi = q`, so they MUST reach the
same :math:`\psi` (``vv-principles`` Mode 9 — verify splittings reach
the same fixed point under anisotropic / :math:`B\neq 0` stressing).
This is a **foundation** test, not an L1 claim: no theory-page
:math:`:label:` is being verified — it pins that three reduction trees
of one affine operator agree on one RHS, which is a software invariant
(``foundation`` NEVER carries ``verifies()``).

T4
(:func:`tests.sn.verification.analytical.test_prescribed_inflow_consistency.test_prescribed_inflow_consistency_si_jacobi_gs_krylov`)
runs two configs. The ``slab_1d`` config (SI is always Jacobi in 1-D)
makes **SI ≡ Krylov** the discriminating pair. The
``cart2d_reflective_y`` config adds reflective-:math:`y` faces so
:math:`B \neq 0` — which is what makes **SI-Jacobi vs SI-Gauss–Seidel**
distinct (G-S folds :math:`B` into the resolvent; Jacobi lags it). The
:math:`B\neq 0`-plus-prescribed-inflow combination is the only config
where the :math:`B`-folding path runs *with* a non-zero boundary source
(the ERR-056 neighbourhood). Measured pairwise reldiffs: 1.3e-13 …
5.6e-13 — comfortably under the 1e-11 ceiling, which itself leaves
headroom for the FP-non-associativity of three reduction trees
(bounded by :math:`\text{iter} \times \text{ULP}` per the
``vv-principles`` bit-identity criteria).

The test carries explicit anti-latent-dud preconditions (the
splitting-invariance check is vacuous if all three trivially agree on
:math:`\psi \equiv 0`): the inflow slot must actually be written
(:math:`>0`), the inflow must non-trivially drive the interior
(:math:`\max|\psi| > 10^{-3}`), and the 2-D row must actually select the
:math:`B_{\rm lower}`-folding reified :math:`M`
(:class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`,
not silently fall back to Jacobi) with an explicit reflective-:math:`y`
``Mesh2D`` BC.

Verification chain — Branch 1 / Branch 2 / L1 cross-check
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Following the ``algebra-of-record`` discipline:

- **Branch 1 (SymPy, State 1C — MMS).** The manufactured sources
  :eq:`sn-mms-nonvacuum-qext` (slab) and
  :eq:`sn-mms-nonvacuum-sph-qext` (sphere) are derived by substituting
  the imposed ansatz into the continuous operator and solving for the
  residual, symbolically. The slab pair is
  :func:`~orpheus.derivations.continuous.mms.sn._nonvacuum_slab_symbolic`
  /
  :func:`~orpheus.derivations.continuous.mms.sn.derive_nonvacuum_slab_mms`;
  the sphere reuses
  :func:`~orpheus.derivations.continuous.mms.sn._spherical_anisotropic_symbolic`
  with the 4.6 shapes
  (:func:`~orpheus.derivations.continuous.mms.sn._nonvacuum_spherical_AB`)
  via
  :func:`~orpheus.derivations.continuous.mms.sn.derive_nonvacuum_spherical_mms`.
  Each ``derive_*`` proves ``simplify(W·LHS − Σ_s·φ − Q_closed) == 0``.
  Foundation gate:
  :file:`tests/derivations/test_sn_mms_nonvacuum_symbolic.py`.
- **Branch 2 (vectorised numpy).** The factories
  :class:`~orpheus.derivations.continuous.mms.sn.SNSlabNonVacuumMMSCase`
  and
  :class:`~orpheus.derivations.continuous.mms.sn.SNSphericalNonVacuumMMSCase`
  (built by
  :func:`~orpheus.derivations.continuous.mms.sn.build_slab_nonvacuum_mms_case`,
  :func:`~orpheus.derivations.continuous.mms.sn.build_slab_2g_nonvacuum_mms_case`,
  and
  :func:`~orpheus.derivations.continuous.mms.sn.build_sphere_nonvacuum_mms_case`)
  evaluate the closed-form source per ordinate using vectorised numpy.
  Each carries a ``prescribed_inflow(sn)`` method returning the
  ``q.boundary`` :class:`~orpheus.transport.source_sinks.AngularBoundarySourceSink`,
  and a ``fixed_source(sn)`` bundler returning the composite
  ``q = q_bulk ⊕ q_∂``
  :class:`~orpheus.transport.timed_full_field.TimedFullField` the public
  solver consumes (see :ref:`sn-composite-fixed-source`).
- **L1 cross-check (the gate).** The Branch-2 numpy
  :math:`Q^{\text{ext}}_n` is bit-equal (≤1e-13 max absolute) to the
  Branch-1 SymPy closed form evaluated via :func:`sympy.lambdify` on a
  sample mesh, for both geometries. The two branches are *structurally
  independent above the trusted-library line* — Branch 1 is
  ``lambdify``-d SymPy, Branch 2 is hand-written numpy — so agreement
  catches a copy error between the symbolic derivation and the
  numerical implementation. Tested in
  :func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_slab_nonvacuum_numerical_qext_matches_sympy`
  and the spherical sibling.

**Structural independence (L11).** The chosen scalar flux
:math:`\phi = A` is *imposed* analytically; the source :math:`Q^{\text{ext}}`
is SymPy-derived (not generated by the solver's own primitives); the
numpy ``external_source`` is then cross-checked bit-equal to the
lambdified SymPy. The reference is structurally independent of the code
under test — the manufactured source does not pass through any of the
solver's discretisation primitives, so the L1 convergence rows are a
genuine test, not a tautology.

**What this section does NOT verify.** Per the three pillars
(``vv-principles``), MMS is a *source-driven* problem: it verifies the
convergence order (a math claim) and the flux shape (a model claim,
because the source is structurally independent), but it **cannot**
verify an eigenvalue. The 4.6 mixtures are non-fissile by construction
and there is no eigenvalue claim anywhere in this section. The
prescribed-inflow verification is a forward-only, fixed-source result.


.. _sn-case-heterogeneous-verification:

Heterogeneous eigenvalue — Case singular-eigenfunction method
--------------------------------------------------------------

Phase 2.1b of the verification campaign closes the last
heterogeneous gap in the SN verification ladder: the
**eigenvalue iteration** on a 1-group two-region reflective
slab, verified against a semi-analytical reference derived
from the discrete-:math:`S_N` slope matrix itself --- no
diffusion approximation, no cross-code comparison, no
Richardson self-test.

The reference is produced by
:func:`orpheus.derivations.continuous.cases.sn.derive_sn_heterogeneous_continuous`
and consumed by
:func:`tests.sn.eigenvalue.test_heterogeneous_transport.test_sn_2region_reflective_case_eigenvalue`
(eigenvalue) and
:func:`tests.sn.eigenvalue.test_heterogeneous_transport.test_sn_2region_reflective_flux_shape`
(scalar flux shape). The Phase 2.1a smooth-:math:`\Sigma` MMS
test verifies the **spatial operator** at :math:`\mathcal O(h^{2})`
design order; this section's Case method verifies the
**eigenvalue** iteration at the material-interface-degraded
:math:`\mathcal O(h)` rate expected for diamond-difference on
piecewise-constant :math:`\Sigma`.

**Motivation: why a second verification path.** The Phase 2.1a
MMS test deliberately uses smooth :math:`\Sigma(x)` to avoid
interface degradation and hit the :math:`\mathcal O(h^{2})`
design order of diamond difference. That is the right choice
for verifying the spatial operator, but it **cannot** exercise
the heterogeneous-interface regime where material
discontinuities force the operator into its interface-layer
behaviour --- the regime where a significant fraction of
production solver bugs live (including ERR-025 — see the
:ref:`homogeneous / uniform-rescale gotcha <sn-homogeneous-degeneracy-gotcha>`
for the mechanism by which it hid). The Case singular-eigenfunction
method provides the complementary reference: an eigenvalue
solution with genuine material-interface discontinuities, built
from the transport equation without running the solver.

**Operator.** The 1-group 1D slab SN transport equation in a
single region with cross sections
:math:`(\Sigma_t, \Sigma_s, \nu\Sigma_f)` and reflective BCs
is, per ordinate,

.. math::
   :label: sn-case-per-ordinate

   \mu_n\,\frac{d\psi_n}{dx} + \Sigma_t\,\psi_n
     \;=\; \frac{c_\text{eff}(k)}{W}\,\phi,
   \qquad
   \phi = \sum_m w_m\,\psi_m,
   \qquad
   c_\text{eff}(k) = \Sigma_s + \frac{\nu\Sigma_f}{k},

where :math:`W = \sum_m w_m`. Substituting the scalar-flux
definition and stacking the angular flux into
:math:`\mathbf y \in \mathbb R^N` (for Gauss--Legendre order
:math:`N`), the system becomes a first-order constant-coefficient
ODE

.. math::
   :label: sn-case-slope-matrix

   \frac{d\mathbf y}{dx} \;=\; \mathbf S(k)\,\mathbf y,
   \qquad
   \mathbf S(k)[n, m] \;=\; \frac{1}{\mu_n}
       \left(-\Sigma_t\,\delta_{nm}
             + \frac{c_\text{eff}(k)}{W}\,w_m\right).

Note the **row-scaling** :math:`1/\mu_n`: the slope matrix is
generally non-symmetric even for symmetric GL quadrature,
because the angular ODE has different "speeds" for different
ordinates.

**Per-region spatial modes.** For each region (fuel at
:math:`x \in [0, H_A]` and moderator at :math:`x \in [H_A, L]`),
diagonalise :math:`\mathbf S(k)`:

.. math::
   :label: sn-case-spatial-modes

   \mathbf S(k)\,\mathbf v_i \;=\; \lambda_i\,\mathbf v_i,
   \qquad i = 1,\ldots,N,

via :func:`numpy.linalg.eig`. For subcritical regions
(:math:`c_\text{eff}(k) < 1`, typical moderator) the eigenvalues
come in :math:`\pm` real pairs. For supercritical regions
(:math:`c_\text{eff}(k) > 1`, fuel at :math:`k` below
:math:`k_{\infty,\text{fuel}}`) some pairs are
complex-conjugate. Each real eigenvalue gives one exponential
mode :math:`\exp(\lambda\,x)\,\mathbf v`; each complex-conjugate
pair gives two real modes built from the canonical
:math:`\cos/\sin/\Re/\Im` combination.

**Real bounded basis.** The naive unbounded basis
:math:`\exp(\lambda\,x)\,\mathbf v` is catastrophically
ill-conditioned for optically thick slabs --- the Phase 1.2
diffusion investigation history records the ``expm``-based
transfer-matrix composition dying from :math:`\text{cond}
\sim 10^{17}` on an 80-cm slab, finding spurious roots with
:math:`\mathcal O(10^{-3})` null-vector residuals rather than
machine-precision zeros. The fix, ported verbatim to Phase 2.1b,
is to **anchor each mode at the nearer region edge**:

.. math::
   :label: sn-case-real-basis

   m^{\text{real}}_j(x) &\;=\; \exp(\lambda_j\,\xi_j)\,\mathbf v_j,
       \qquad
       \xi_j = \begin{cases}
         x - L_\text{reg} & \lambda_j \ge 0 \;\;\text{(anchor right)} \\
         x                & \lambda_j < 0 \;\;\text{(anchor left)}
       \end{cases} \\[1mm]
   m^{\text{c}}_j(x) &\;=\; e^{\Re\lambda_j\,\xi_j}\,
       \bigl(\cos(\Im\lambda_j\,\xi_j)\,\mathbf v_{R,j}
          - \sin(\Im\lambda_j\,\xi_j)\,\mathbf v_{I,j}\bigr), \\
   m^{\text{s}}_j(x) &\;=\; e^{\Re\lambda_j\,\xi_j}\,
       \bigl(\sin(\Im\lambda_j\,\xi_j)\,\mathbf v_{R,j}
          + \cos(\Im\lambda_j\,\xi_j)\,\mathbf v_{I,j}\bigr),

where :math:`\mathbf v_j = \mathbf v_{R,j} + i\,\mathbf v_{I,j}`
is the complex eigenvector. Every mode is bounded by
:math:`|\mathbf v_j|` on its region, so the assembled matching
matrix has :math:`\mathcal O(1)` entries.

**Matching matrix.** For the 2-region reflective slab the
coefficient vector has dimension :math:`2N` (one real mode per
eigenvalue per region). The linear constraints are:

.. math::
   :label: sn-case-matching-matrix

   &\text{Reflective at } x = 0:\quad
      \psi^A_n(0) - \psi^A_{N-1-n}(0) = 0,
      \qquad n \in [0, N/2) \\[1mm]
   &\text{Interface at } x = H_A:\quad
      \psi^A_n(H_A) - \psi^B_n(H_A) = 0,
      \qquad n \in [0, N) \\[1mm]
   &\text{Reflective at } x = L:\quad
      \psi^B_n(L) - \psi^B_{N-1-n}(L) = 0,
      \qquad n \in [0, N/2)

:math:`N/2 + N + N/2 = 2N` equations in :math:`2N` unknowns.
The partner index :math:`N-1-n` is the Gauss--Legendre
reflection pairing (ordinates sorted by ascending :math:`\mu`).
The eigenvalue condition is
:math:`\det\mathbf C(k) = 0`.

**Root finding.** :func:`scipy.optimize.brentq` on
:math:`\det\mathbf C(k)` over a coarse :math:`k`-scan, with
sign-change bracketing, refines every candidate to
``xtol=1e-14``. But :func:`numpy.linalg.eig`'s eigenvalue
ordering is not a continuous function of :math:`k` --- at
parameter values where two per-region eigenvalues cross, the
eigenvalue labels permute discontinuously, and
:math:`\det\mathbf C(k)` flips sign by permutation rather than
by passing through zero. brentq will "converge" to such
spurious points.

**Physical validation.** Every candidate root is rebuilt via
SVD of :math:`\mathbf C(k)`, and the null vector's reflective-BC
residuals at :math:`x = 0` and :math:`x = L`, and the interface
continuity residual at :math:`x = H_A`, are explicitly
reconstructed and checked against a dimensionless tolerance
relative to the peak angular flux:

.. math::
   :label: sn-case-physical-validation

   \|\psi(0, +\mu_n) - \psi(0, -\mu_n)\| / \|\psi\|_\text{peak}
     &< \text{tol} \\
   \|\psi^A(H_A) - \psi^B(H_A)\| / \|\psi\|_\text{peak}
     &< \text{tol} \\
   \|\psi(L, +\mu_n) - \psi(L, -\mu_n)\| / \|\psi\|_\text{peak}
     &< \text{tol}

Only candidates passing all three are accepted; the fundamental
is the largest validated root. This is the SN analogue of the
Phase 1.2 diffusion physical validation (same pattern, different
operator).

**Back-substitution.** Once :math:`k_\text{fund}` is found,
the null vector at that :math:`k` is the coefficient vector in
the :math:`2N`-dimensional real basis. Evaluation of
:math:`\phi(x) = \sum_n w_n\,\psi_n(x)` at any point reduces to
a linear combination of a handful of bounded exponential or
trigonometric modes:

.. math::
   :label: sn-case-back-substitution

   \psi(x) = \begin{cases}
     \sum_j c^A_j\,m^A_j(x) & x \le H_A \\[1mm]
     \sum_j c^B_j\,m^B_j(x - H_A) & x > H_A
   \end{cases},
   \qquad
   \phi(x) = \sum_n w_n\,\psi_n(x).

All modes are bounded by :math:`\mathcal O(1)`, so
:math:`\phi(x)` is stable to machine precision.

**The Phase 2.1b diagnostic configuration.** The canonical
test problem is the ``A`` + ``B`` 1-group mixture pair from
:mod:`orpheus.derivations.common.xs_library`:

.. list-table::
   :header-rows: 1
   :widths: 15 15 15 15 15

   * - Region
     - :math:`\Sigma_t`
     - :math:`\Sigma_s`
     - :math:`\nu\Sigma_f`
     - :math:`k_\infty`
   * - A (fuel)
     - 1.0
     - 0.5
     - 0.75
     - 1.5
   * - B (moderator)
     - 2.0
     - 1.9
     - 0
     - ---

with :math:`H_A = H_B = 0.5\,\text{cm}`, reflective BCs on both
outer edges, :math:`S_8` Gauss--Legendre quadrature. The
resulting Case reference is

.. math::

   k_\text{eff}^{\text{Case}}(S_8) = 1.2746160417

--- the exact discrete-:math:`S_8` eigenvalue. For
cross-validation, the same configuration run through ORPHEUS's
:func:`~orpheus.cp.solver.solve_cp` (1D slab E\ :sub:`3` kernel,
completely independent numerical path) gives
:math:`k^{\text{CP}} = 1.2744284665` --- agreement to
:math:`\sim 2\times 10^{-4}`, well below the :math:`\mathcal O(1\%)`
difference that typically exists between discrete-SN and
continuous-angle formulations. This cross-check is used only as
a sanity input, not as a verification crutch.

**Measured convergence.** With :math:`S_8`, refining
:math:`n_\text{per}` per region:

.. list-table::
   :header-rows: 1
   :widths: 15 25 15

   * - :math:`n_\text{per}`
     - :math:`k_\text{solve}`
     - :math:`|k_\text{solve} - k_\text{Case}|`
   * - 20
     - 1.2746074093
     - :math:`\sim 8.6\!\times\!10^{-6}`
   * - 40
     - 1.2746138837
     - :math:`\sim 2.2\!\times\!10^{-6}`
   * - 80
     - 1.2746155022
     - :math:`\sim 5.4\!\times\!10^{-7}`
   * - 160
     - 1.2746159068
     - :math:`\sim 1.3\!\times\!10^{-7}`
   * - 320
     - 1.2746160080
     - :math:`\sim 3.4\!\times\!10^{-8}`

Each refinement roughly halves the error, confirming the
:math:`\mathcal O(h)` rate expected at a material interface with
piecewise-constant :math:`\Sigma`. The finest-mesh residual of
:math:`3.4 \times 10^{-8}` is **machine-precision agreement**
between two independent mathematical constructions (the Case
matching-matrix + back-substitution reference and the
diamond-difference sweep-based power iteration); both
implementations solve the same discrete-:math:`S_N` spectral
problem and agree to within the BiCGSTAB-compatible
truncation.

**Contrast with Phase 2.1a.** The Phase 2.1a MMS section hits
:math:`\mathcal O(h^{2})` because it uses smooth
:math:`\Sigma(x)`; the Phase 2.1b Case section hits
:math:`\mathcal O(h)` because it uses piecewise-constant
:math:`\Sigma(x)` with a genuine material interface. Both are
correct for their respective regimes. The degradation from
:math:`h^{2}` to :math:`h` at the interface is the standard
Salari--Knupp result for DD on discontinuous coefficients, and
is the **reason** Phase 2.1a deliberately chose smooth
:math:`\Sigma` to isolate the spatial operator.


Homogeneous Infinite Medium
----------------------------

For homogeneous geometry with reflective BCs, the flux is spatially flat
and :math:`\keff = \lambda_{\max}(A^{-1}F)`.  This is geometry-independent
--- Cartesian, spherical, and cylindrical must all give the same
:math:`\keff`.

.. list-table::
   :header-rows: 1
   :widths: 10 14 19 19 19 19

   * - Groups
     - :math:`\kinf`
     - Cartesian (GL S8)
     - Spherical (GL S8)
     - Cylindrical (Prod 4x8)
     - Cylindrical (LS S4)
   * - 1
     - 1.5000
     - exact
     - exact
     - exact
     - exact
   * - 2
     - 1.8750
     - exact
     - exact
     - exact
     - exact
   * - 4
     - 1.4878
     - exact
     - exact
     - exact
     - exact

All entries are exact to machine precision.  Spherical 2G/4G results
(previously showing ~1% error) are now exact thanks to the M-M angular
closure weights.

Heterogeneous Convergence
--------------------------

For a cylindrical fuel (r < 0.5) + moderator (r < 1.0) geometry with
Product(4x8) quadrature:

.. list-table::
   :header-rows: 1
   :widths: 20 25 25

   * - Cells/region
     - :math:`\keff` (1G)
     - :math:`\Delta k` from previous
   * - 5
     - 0.9769
     -
   * - 10
     - 0.9842
     - +0.0073
   * - 20
     - 0.9874
     - +0.0032

:math:`\keff` converges monotonically toward the CP reference
(0.9955).  The ~1% residual gap is the white-BC (CP) vs reflective-BC
(SN) approximation difference, consistent with the slab geometry
findings.

For the 2G heterogeneous resolution test, Product(4x8) and Product(8x8)
agree to :math:`< 0.01\%` (keff = 0.7227 for both), confirming
angular convergence.

Why 1-Group Verification Is Degenerate
----------------------------------------

For 1 energy group, the eigenvalue is:

.. math::

   k = \frac{\nSigf{}}{\Sigma_a}

This is a scalar ratio independent of the spatial or angular flux
distribution.  Consequences:

- Angular weight errors scale all flux equally --- cancels in :math:`k`.
- Wrong scattering convention --- no inter-group coupling to distort.
- Wrong flux shape --- does not matter; :math:`k` is a material property.

Only multi-group problems have a flux-shape-dependent eigenvalue:
:math:`k = (\nSigf{} \cdot \phi) / (\Sigma_a \cdot \phi)` where the
dot product weights each group differently.  A wrong group ratio (from
angular errors, normalization errors, or convergence failures) directly
shifts :math:`\keff`.

**Rule:** Every transport solver must be verified on at least 2-group
problems.  1-group success gives false confidence.

Spatial and Angular Convergence
--------------------------------

The diamond-difference scheme converges at :math:`O(h^2)` with mesh
refinement.  Gauss--Legendre quadrature shows spectral convergence in
angle.  Both are verified in ``test_sn_1d.py``.

Property Tests
---------------

For all geometries:

- **Particle balance**: production / absorption :math:`= \keff`
- **Flux non-negativity**: :math:`\phi \geq 0` everywhere
- **Angular flux at** :math:`r = 0` **all positive** (curvilinear only)
- **Multi-group eigenvector not flat**: flux spectrum differs between
  fuel and moderator (catches 1G-degenerate bugs)

Run the full suite::

   pytest tests/sn/ -v -m "not slow"


Numerical Sensitivities
========================

:math:`\keff` Sensitivity Table (421-Group Heterogeneous PWR Slab)
-------------------------------------------------------------------

All cases: 10 cells, :math:`\delta = 0.2` cm, material layout
``[fuel x 5, clad x 1, cool x 4]``, P0 scattering, 421 energy groups.

.. list-table::
   :header-rows: 1
   :widths: 50 15 35

   * - Configuration
     - :math:`\keff`
     - Notes
   * - 1D GL S16, BiCGSTAB (FD operator)
     - 1.03882
     - True 1D, 16 ordinates
   * - 1D Lebedev 110, source iteration (DD sweep)
     - 1.04294
     - 1D mesh, 2D quadrature
   * - 2D (10x2) Lebedev 110, source iter (DD sweep)
     - 1.04294
     - Pseudo-2D, full volumes
   * - 2D (10x2) Lebedev 110, BiCGSTAB (FD)
     - 1.04007
     - Pseudo-2D, full volumes
   * - 2D (10x2) Lebedev 110, BiCGSTAB, half-volumes
     - 1.04192
     - MATLAB convention
   * - **MATLAB reference**
     - **1.04188**
     - 2D Lebedev, FD, half-volumes

Sources of Variation
---------------------

1. **Angular quadrature** (GL vs Lebedev): ~0.004 difference.
   GL S16 integrates 1D angular flux with 16 points on :math:`[-1,1]`.
   Lebedev 110 integrates over the unit sphere --- more angular
   resolution but different effective weights per :math:`\mu_x`
   direction.  On a coarse heterogeneous mesh, these give different
   eigenvalues.

2. **Spatial discretisation** (DD sweep vs FD gradient): ~0.003
   difference.  Source iteration uses the DD wavefront sweep
   (:math:`T^{-1}`).  BiCGSTAB uses the explicit FD transport operator
   (:math:`T`).  Both are :math:`O(h)` on this mesh but with different
   truncation error constants.

3. **Boundary volume weighting**: ~0.002 difference (full vs half).
   The MATLAB code halves boundary cell volumes.  With ``ny=2`` and
   materials uniform in *y*, only the *x*-direction halving (fuel edge,
   coolant edge) affects :math:`\keff`.  This is an artifact of the
   pseudo-2D implementation: a true 1D calculation has no *y*-volumes.

4. **Inner convergence**: source iteration with ``max_inner=200``,
   ``inner_tol=1e-8`` does not fully converge for 421 groups (spectral
   radius ~0.97).  BiCGSTAB fully converges the inner solve in ~100
   Krylov iterations.

Matching the MATLAB Result
---------------------------

The MATLAB code uses: 2D Lebedev 110 on a 10x2 mesh, explicit FD
operator with BiCGSTAB, boundary half-volumes, P0 scattering.

The BiCGSTAB path with half-volumes reproduces 1.04192 vs MATLAB's
1.04188 (:math:`4 \times 10^{-5}` agreement).  The residual difference
is from floating-point details in cross-section processing.

The cleanest reference is the **1D GL BiCGSTAB** result (1.03882): no
pseudo-2D artifacts, well-conditioned angular quadrature, fully
converged inner solve.


