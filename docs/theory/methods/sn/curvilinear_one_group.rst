.. _sn-curvilinear-one-group:

Curvilinear one-group: angle enters the walk
============================================

This chapter opens Part B of the S\ :sub:`N` book: **curvilinear
geometry** — the 1-D sphere and the axially-infinite cylinder — at one
energy group.  Exactly one thing is new.  In Cartesian geometry a
neutron's direction is constant in flight, so ordinates couple only
through the sources; on a curved coordinate frame the *same straight
flight line* changes its **local** direction coordinate as it advances
— the radial cosine (:math:`\mu` on the sphere, :math:`\eta` on the
cylinder) drifts toward its most-outward value as the ray passes its
point of closest approach.  The transport operator therefore gains an
**angular redistribution** term — a conservative derivative in the
direction coordinate — and ordinates become **sequentially coupled**,
along the full :math:`\mu` sequence (sphere) or within each
:math:`\mu`-level (cylinder).  The walk acquires an angular thread.

The chain of the book repeats on the new axis:

1. **the invariant** — sinks = sources on a cell that is now a shell
   :math:`\times` an angular cell: curvature moves neutrons between
   neighbouring ordinates without absorbing them, so on a flat
   isotropic flux the redistribution must cancel the streaming *per
   ordinate* → pose the conservative balance :eq:`conservative-form`
   and force the geometry factor :math:`\Delta A/w` from that
   invariant;
2. **the operator** — :math:`A = L + C - S - B` keeps its honest shape
   and :math:`L+C` stays **lower-triangular**: sweep order now threads
   the angular cascade (seeded at each level's starting direction,
   with :math:`\alpha_{1/2} = 0` closing the recursion) through the
   radial march;
3. **the matrix picture** — the per-cell system gains one more
   upstream state (the angular half-flux :math:`\psi_{n-1/2}`) and
   one more closure weight (the Morel--Montry :math:`\tau`);
4. **the strategy encodings** — the cell-update contract
   (:ref:`cell-update-strategies`) is unchanged: the same Step/DD/LD
   algebra serves slab and curvilinear alike, only the populated
   fields of the streaming packet change.

What does *not* change: space stays 1-D radial (the multi-D walk of
:doc:`cartesian_multid` is not needed), the energy axis stays a single
group (multigroup returns in the next Part-B chapter), and the closure
machinery is still the generic
:doc:`/theory/foundations/discretization`, cross-linked never
re-derived.

.. admonition:: Key Facts
   :class: tip

   * The curvilinear balance :eq:`balance-general` adds two structures
     to the slab balance: the **redistribution cascade**
     (:eq:`alpha-recursion`, :math:`\alpha_{1/2} = 0`, a non-negative
     dome peaking near :math:`\mu = 0`) and the **geometry factor**
     :math:`\Delta A_i/w_n` — forced by per-ordinate flat-flux
     consistency, not optional (without it the solver manufactures
     angular anisotropy that *worsens* under refinement near
     :math:`r = 0`; failure mode #3).
   * The **streaming-equilibrium identity** :eq:`streaming-equilibrium`
     (:math:`\phi = Q/(\Sigma_t(1-c))`, :math:`\psi_n = \phi/W`) is
     the canonical L0 gate — asserted **per ordinate**, never via
     particle balance (telescoping sums hide per-ordinate errors; vv
     anti-pattern #8).
   * The angular closure is **weighted diamond** :eq:`wdd-closure`
     with **Morel--Montry weights** :eq:`mm-weights`: :math:`\tau`
     drives the diffusion-limit contamination :math:`\beta` to
     machine zero, eliminating the M-M flux dip.  **Sphere
     unclamped**; **cylinder clipped** to :math:`[\tfrac12, 1]`
     (product quadratures put :math:`\tau = 0` on the most-inward
     azimuthal ordinate).
   * :math:`\tau` is an **angular-scheme property the closure owns**
     (#236 Phase 2): produced solely by the pole-angular closure and
     delivered to the stateless spatial scheme as
     :class:`~orpheus.transport.spatial.scheme.CellVisit` **data**
     (:math:`c_{\rm in}`, :math:`c_{\rm out}`, :math:`\tau`),
     stamped at one production site — never as a closure dependency
     (the spatial :math:`\otimes` angular separation).
   * The error axes obey the **geometry-split law**
     :eq:`sn-space-angle-separability`: Cartesian **separates**
     (:math:`E \approx E_{\rm space} + E_{\rm angle}`), curvilinear
     **gates** (:math:`E \approx \max(E_{\rm space}, E_{\rm
     angle})`) — fine-:math:`h` accuracy cannot be harvested at
     coarse :math:`N`; pinned by the ST5 characterisation gate.
   * The solved per-cell form is :eq:`dd-solve` — the WDD closure
     substituted into the balance: one more upstream state
     (:math:`\psi_{n-\frac12}`) in the numerator and the
     :math:`(\Delta A_i/w_n)\,c_{\rm out}` shift in the
     denominator, relative to the slab update.

.. _balance-curvilinear:

Curvilinear Balance Equation (Spherical and Cylindrical)
=========================================================

Derivation from the Continuous PDE
------------------------------------

Start with the general 1D curvilinear transport equation.  In
conservative form for a coordinate :math:`r` with face area
:math:`A(r)` and volume element :math:`V`:

.. math::
   :label: conservative-form

   \frac{\mu_n}{V_i}
   \bigl[A_{i+\frac12}\psi_{i+\frac12} - A_{i-\frac12}\psi_{i-\frac12}\bigr]
   + \frac{1}{V_i}
   \bigl[\alpha_{n+\frac12}\psi_{n+\frac12} - \alpha_{n-\frac12}\psi_{n-\frac12}\bigr]
   + \Sigt{} \psi_{n,i} = S_i

.. vv-status: conservative-form documented

where the streaming cosine is :math:`\mu_n` for spherical and
:math:`\eta_m` for cylindrical, and :math:`S_i = Q_i / W` is the
isotropic source density divided by the quadrature weight sum.

**Step 1: Integrate the PDE over a spatial cell.**

For spherical geometry, integrating :eq:`transport-spherical` over the
shell :math:`[r_{i-1/2}, r_{i+1/2}]` and using the divergence theorem
on the radial streaming gives:

.. math::

   \mu_n \bigl[A_{i+\frac12}\psi_{i+\frac12} - A_{i-\frac12}\psi_{i-\frac12}\bigr]
   + \int_{V_i} \frac{1-\mu^2}{r} \frac{\partial\psi}{\partial\mu}\, dV
   + \Sigt{} V_i \psi_{n,i} = S_i V_i

For cylindrical geometry, integrating :eq:`transport-cylindrical` over
the annular shell gives:

.. math::

   \eta_m \bigl[A_{i+\frac12}\psi_{i+\frac12} - A_{i-\frac12}\psi_{i-\frac12}\bigr]
   - \int_{V_i} \frac{1}{r} \frac{\partial(\xi\psi)}{\partial\varphi}\, dV
   + \Sigt{} V_i \psi_{m,i} = S_i V_i

**Step 2: Discretise the angular redistribution.**

The angular integral is discretised as a finite difference in the
ordinate index.  For spherical:

.. math::

   \int_{V_i} \frac{1-\mu^2}{r}\frac{\partial\psi}{\partial\mu}\, dV
   \;\approx\;
   \alpha_{n+\frac12}\psi_{n+\frac12} - \alpha_{n-\frac12}\psi_{n-\frac12}

For cylindrical (per :math:`\mu`-level):

.. math::

   -\int_{V_i} \frac{1}{r}\frac{\partial(\xi\psi)}{\partial\varphi}\, dV
   \;\approx\;
   \alpha_{m+\frac12}\psi_{m+\frac12} - \alpha_{m-\frac12}\psi_{m-\frac12}

**Step 3: Apply the geometry factor** :math:`\Delta A / w`.

The raw discretisation above does NOT preserve per-ordinate flat-flux
consistency.  The correct form from [BaileyMorelChang2010]_ includes the
geometry factor :math:`\Delta A_i / w_n`:

.. math::
   :label: balance-general

   \mu_n
   \bigl[A_{i+\frac12}\psi_{i+\frac12} - A_{i-\frac12}\psi_{i-\frac12}\bigr]
   + \frac{\Delta A_i}{w_n}
   \bigl[\alpha_{n+\frac12}\psi_{n+\frac12} - \alpha_{n-\frac12}\psi_{n-\frac12}\bigr]
   + \Sigt{} V_i \psi_{n,i} = S_i V_i

where :math:`\Delta A_i = A_{i+1/2} - A_{i-1/2}`.  This is the curvilinear
balance form of [BaileyMorelChang2010]_ for both spherical and cylindrical
geometry.

Note why :eq:`dd-cartesian-1d` has no :math:`\alpha` or :math:`\Delta A`
terms: in Cartesian geometry the face area is unity (:math:`A = 1`), so
:math:`\Delta A = 0`, and there is no curvature to redistribute angular
flux.

The Alpha Redistribution Coefficients
======================================

The :math:`\alpha` coefficients encode how the angular flux redistributes
between neighbouring ordinates due to the geometry curvature.  They are
defined recursively:

.. math::
   :label: alpha-recursion

   \alpha_{n+\frac12} = \alpha_{n-\frac12} - w_n \mu_n

with the boundary condition :math:`\alpha_{1/2} = 0`.

For **spherical** geometry, all :math:`N` ordinates form a single
sequence sorted by :math:`\mu` (most negative to most positive).
The :math:`\alpha` values form a **non-negative dome**: they rise while
:math:`\mu < 0`, peak near :math:`\mu = 0`, and fall back to zero at
:math:`\mu = 1`.  The endpoint condition
:math:`\alpha_{N+1/2} = 0` is guaranteed by Gauss--Legendre
antisymmetry: :math:`\sum_n w_n \mu_n = 0`.

For **cylindrical** geometry, each :math:`\mu`-level has its own
independent :math:`\alpha` sequence.  On level :math:`p`, the ordinates
are sorted by increasing :math:`\eta` (radial direction cosine), and the
recursion uses :math:`\eta` instead of :math:`\mu`:

.. math::
   :label: alpha-cylindrical

   \alpha_{p,m+\frac12} = \alpha_{p,m-\frac12} - w_m \eta_m

This is the [BaileyMorelChang2010]_ curvilinear :math:`\alpha`-recursion.
Each level's :math:`\alpha` values form an independent dome from
:math:`\eta = -\sin\theta` to :math:`\eta = +\sin\theta`.

**Dome shape properties:**

- :math:`\alpha_{n+1/2} \geq 0` for all :math:`n` (non-negative dome).
- The peak occurs near the ordinate where :math:`\mu_n` (or
  :math:`\eta_m`) crosses zero.
- The dome height scales with the quadrature weight sum: higher-order
  quadratures have narrower but taller domes.
- Non-negativity ensures the denominator of the DD equation is
  unconditionally positive, guaranteeing numerical stability.

The code stores these on the reduced streaming operator —
``mesh.reduced.alpha_half`` (spherical, shape ``(N+1,)``) and
``mesh.reduced.alpha_per_level`` (cylindrical, list of ``(M+1,)``
arrays); they are genuinely geometric and stay on the geometry side
(:ref:`sn-tau-step-c-closeout`).

The Geometry Factor and Why It Is Needed
=========================================

The geometry factor :math:`\Delta A_i / w_n` in :eq:`balance-general`
is the key to correct curvilinear transport.  Without it, the balance
equation violates **per-ordinate flat-flux consistency**: for a spatially
uniform, isotropic flux :math:`\psi = \text{const}`, the streaming and
redistribution terms should cancel exactly for EACH ordinate
individually.

**Proof of consistency.**

Set :math:`\psi_{n,i} = \psi_{n+1/2} = \psi_{n-1/2} = \psi_0` (flat
in both space and angle) and :math:`\psi_{i+1/2} = \psi_{i-1/2} = \psi_0`
(flat in space).  The streaming term becomes:

.. math::

   \mu_n \bigl[A_{i+\frac12} - A_{i-\frac12}\bigr] \psi_0
   = \mu_n \,\Delta A_i\, \psi_0

The redistribution term with the :math:`\Delta A/w` factor becomes:

.. math::

   \frac{\Delta A_i}{w_n}
   \bigl[\alpha_{n+\frac12} - \alpha_{n-\frac12}\bigr] \psi_0
   = \frac{\Delta A_i}{w_n} (-w_n \mu_n) \psi_0
   = -\mu_n \,\Delta A_i\, \psi_0

where we used the recursion :eq:`alpha-recursion`:
:math:`\alpha_{n+1/2} - \alpha_{n-1/2} = -w_n \mu_n`.  The two terms
cancel exactly, giving :math:`\Sigt{} \psi_0 = S_0`, which is the
correct homogeneous solution.

**Without** the :math:`\Delta A/w` factor (i.e., using
:math:`[\alpha_{n+1/2}\psi_{n+1/2} - \alpha_{n-1/2}\psi_{n-1/2}]`
directly), the redistribution term for flat flux is
:math:`(-w_n \mu_n)\psi_0`, but the streaming term is
:math:`\mu_n \Delta A_i \psi_0`.  These differ by a factor of
:math:`\Delta A_i`, so consistency only holds in the limit
:math:`\Delta A_i \to 0` (i.e., at the origin or on an infinitely fine
mesh).

**Consequence of the missing factor:**  The solver creates artificial
angular anisotropy that *worsens* with mesh refinement near :math:`r = 0`
(where :math:`\Delta A_i` is smallest but non-zero).  This manifests as
a flux spike at the origin in fixed-source problems and as divergent
eigenvalues in heterogeneous eigenvalue problems.

The code precomputes this factor as ``mesh.reduced.redist_dAw``
(spherical, shape ``(nx, N)``) and ``mesh.reduced.redist_dAw_per_level``
(cylindrical, list of ``(nx, M)`` arrays).

The Streaming-Equilibrium Identity (canonical L0 gate)
=======================================================

The flat-flux consistency proof above is the *per-cell* statement of a
*global* exact solution that the verification suite leans on harder
than any other: for a *homogeneous* medium with a *uniform isotropic*
source :math:`Q` per group and boundaries that sustain a flat flux
(reflective faces, or an infinite/periodic medium), the discrete
transport equation has the exact fixed point

.. math::
   :label: streaming-equilibrium

   \phi \;=\; \frac{Q}{\Sigma_t\,\bigl(1 - c\bigr)},
   \qquad
   \psi_n \;=\; \frac{\phi}{W}
   \quad \forall n,
   \qquad
   W \equiv \sum_n w_n,
   \quad
   c \equiv \frac{\Sigma_{s0}}{\Sigma_t},

per group. For a pure-attenuation configuration (no scattering in the
residual, :math:`c = 0`) this reduces to :math:`\phi = Q/\Sigma_t` with
the per-ordinate angular flux :math:`\psi_n = Q/(W\,\Sigma_t)`.

**Why this is the canonical L0 gate.** Substituting the flat
:math:`\psi_n` into the discrete balance :eq:`balance-general` nulls
the streaming and redistribution terms *per ordinate* (the proof of
consistency above), leaving the pure collision balance
:math:`\Sigma_t\,\psi_n = Q/W + \Sigma_{s0}\,\phi/W`. Every term that
a discretisation can get wrong — a missing :math:`\Delta A/w` factor
(failure mode #3, the flux spike at :math:`r=0`), a wrong
:math:`\alpha` recursion (mode #4), a face-index slip (mode #5), a
weight-normalisation drift (:math:`1/W` vs :math:`1/4\pi`) — breaks
the identity at machine precision, with no discretisation error to
hide behind. The assertion is **per-ordinate**, never
particle-balance: telescoping global balance holds by construction
even when per-ordinate balance is wrong (the canonical ERR-006 hide;
vv-principles anti-pattern #8).

The identity holds in every geometry ORPHEUS supports (slab, sphere,
cylinder, 2-D/3-D Cartesian) and at both algebraic access points —
the sweep (``solve``: given :math:`Q`, recover the flat
:math:`\psi`) and the matvec (``apply``: given the flat :math:`\psi`,
recover the residual :math:`Q` with no spurious boundary or pole
contribution). Tests declare it via
``@pytest.mark.verifies("streaming-equilibrium")``.

The Morel--Montry Flux Dip
============================

Even with the correct :math:`\Delta A/w` factor, the standard
diamond-difference closure (equal weight :math:`\tau = 0.5`) introduces
a flux error near :math:`r = 0` known as the **Morel--Montry flux dip**
[MorelMontry1984]_.

The standard DD angular closure is:

.. math::

   \psi_{n,i} = \frac{1}{2}(\psi_{n-\frac12} + \psi_{n+\frac12})

This can be rewritten as:

.. math::

   \psi_{n+\frac12} = 2\psi_{n,i} - \psi_{n-\frac12}

The contamination factor :math:`\beta` ([BaileyMorelChang2010]_) quantifies
the coupling between the leading-order scalar flux and the first-order
current in the asymptotic diffusion limit.  For spherical geometry:

.. math::

   \beta = \frac{1}{2} \sum_{n=1}^{N} \mu_n
   \bigl[\alpha_{n+\frac12}\, \mu_{n+\frac12}
        - \alpha_{n-\frac12}\, \mu_{n-\frac12}\bigr]

where :math:`\mu_{n\pm 1/2}` are the angular cell-edge cosines.  For
cylindrical, the equivalent is a per-level sum using :math:`\eta` and
:math:`\eta_{m\pm 1/2}`.  When :math:`\beta \neq 0`, the discrete
S\ :sub:`N` equations satisfy a **contaminated** diffusion equation near
:math:`r = 0`, producing the artificial flux dip (or spike).

The module :mod:`orpheus.derivations.discrete.sn.contamination`
computes :math:`\beta`
for any quadrature and geometry.  With the correct :math:`\Delta A/w`
factor AND Morel--Montry weights, :math:`\beta \sim 10^{-16}`
(machine zero) for both spherical and cylindrical.

Weighted Diamond Difference (WDD) and Morel--Montry Weights
=============================================================

The Morel--Montry (M-M) angular closure replaces the equal-weight DD
with position-dependent weights :math:`\tau_n` [BaileyMorelChang2010]_ Eq. 43:

.. math::
   :label: wdd-closure

   \psi_{n,i} = (1 - \tau_n)\,\psi_{n-\frac12} + \tau_n\,\psi_{n+\frac12}

Solving for the angular face flux:

.. math::
   :label: wdd-face

   \psi_{n+\frac12}
   = \frac{\psi_{n,i} - (1 - \tau_n)\,\psi_{n-\frac12}}{\tau_n}

The M-M weights are defined as:

.. math::
   :label: mm-weights

   \tau_n = \frac{\mu_n - \mu_{n-\frac12}}{\mu_{n+\frac12} - \mu_{n-\frac12}}

where :math:`\mu_{n\pm 1/2}` are the angular cell edges.

**Spherical cell edges:**  :math:`\mu_{1/2} = -1`,
:math:`\mu_{N+1/2} = +1`, and interior edges by weight-sum:
:math:`\mu_{n+1/2} = \mu_{n-1/2} + w_n`.  This is exact for
Gauss--Legendre quadrature because the weights correspond to the
:math:`\mu`-space widths of the angular cells.

**Cylindrical cell edges:**  :math:`\eta_{1/2} = -\sin\theta`,
:math:`\eta_{M+1/2} = +\sin\theta`, and interior edges at
**midpoints** of consecutive :math:`\eta` values:
:math:`\eta_{m+1/2} = (\eta_m + \eta_{m+1})/2`.
The weight-sum approach is NOT used for cylindrical because the
quadrature weights are uniform in :math:`\varphi`-space (not
:math:`\eta`-space): the Product quadrature spaces :math:`\varphi`
equally, but :math:`\eta = \sin\theta\cos\varphi` is
cosine-distributed, so equal :math:`\varphi`-widths map to unequal
:math:`\eta`-widths.  The midpoint approach gives a proper partition
of :math:`[-\sin\theta, +\sin\theta]`.

For the Product quadrature with equally-spaced :math:`\varphi`,
ordinates come in **pairs** with the same :math:`|\eta|` but opposite
:math:`\xi` (e.g., :math:`\varphi = \pi/4` and :math:`\varphi = 7\pi/4`
both give :math:`\eta = \sin\theta/\sqrt{2}`).  The midpoint between
paired ordinates equals their shared :math:`\eta`, creating zero-width
angular cells.  The resulting :math:`\tau` alternates between 0.5
(DD, for the left member of each pair) and 1.0 (step, for the right
member).  This alternating pattern is correct but could be smoothed by
using a Gauss-type azimuthal quadrature with distinct :math:`\eta`
values (see `GitHub Issue #1 <https://github.com/deOliveira-R/ORPHEUS/issues/1>`_).

The M-M weights force the contamination factor :math:`\beta` to **machine
zero** (verified: :math:`\beta \sim 10^{-16}`), completely eliminating
the Morel--Montry flux dip.

**Clipping (geometry-dependent since W1, 2026-06-13).**  The raw
weight :eq:`mm-weights` is the unique weight exact for a flux linear in
:math:`\mu` (Bailey-Morel-Chang 2010 Eq. 43), admissible range
:math:`\tau \in [0, 1]`.

* **Sphere** uses :math:`\tau_n` **unclamped**.  On Gauss--Legendre
  quadrature :math:`\tau_n \in [0.39, 0.61]` (never 0), so the closure
  is positive without a clamp.  The former :math:`[0.5, 1.0]` clamp was
  an over-conservative, mis-cited positivity floor that re-floored the
  anisotropic solution; the full vindication is at
  :ref:`sn-curvilinear-aniso-norm-reconciliation`.
* **Cylinder** retains the clip :math:`\tau_n \to
  \mathrm{clip}(\tau_n, 0.5, 1.0)`: product / level-symmetric
  quadratures give :math:`\tau_n = 0` exactly on the most-inward
  azimuthal ordinate (a structural :math:`\div 0` block; the alternating
  0.5 / 1.0 pattern described above is the clamp acting on these
  zero-width angular cells).

The clamp lives **inside the angular closure**, in
:func:`~orpheus.sn.sweep.pole_angular_closure.morel_montry_tau_per_level`
(sphere unclamped, cylinder clipped to :math:`[\tfrac12, 1]`).  The
closure exposes the resulting weight per global ordinate through
:attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.tau_per_ordinate`
(spherical: a single :math:`(N,)` array; cylindrical: the per-level
weights gathered to the global ordinate order).  Issue #236 Phase 2
retired the parallel geometry-side :math:`\tau` producer that formerly
baked these weights onto :class:`StreamingTerms`; :math:`\tau` is now
produced **solely** by the closure (see :ref:`sn-tau-c-on-cellvisit-live`).

.. _sn-tau-closure-owned:

τ is an angular-scheme property — the closure owns it (Issue #236 Phase 2)
------------------------------------------------------------------------------

.. todo:: Archivist expansion needed.

   The Morel--Montry weight :math:`\tau` :eq:`mm-weights` is a function of
   the quadrature :math:`(\mu, w, \text{levels})` ALONE — an
   ANGULAR-scheme property, not a geometry one.  Issue #236 Phase 2
   (Step A) relocates :math:`\tau` PRODUCTION onto the angular closure:
   :func:`~orpheus.sn.sweep.pole_angular_closure.morel_montry_tau_per_level`
   produces :math:`\tau` from the quadrature the closure already binds,
   and
   :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
   consumes its OWN :math:`\tau` for the matvec contribution (P0) instead
   of reading it back from the streaming-geometry factory
   (:func:`~orpheus.geometry.reduced_operator.spherical_streaming` /
   :func:`~orpheus.geometry.reduced_operator.cylindrical_streaming`).

   Step A is BIT-IDENTICAL: the producer is a 0-ULP line-for-line replica
   of the factory arithmetic (sphere unclamped, cylinder clamped to
   :math:`[\tfrac12, 1]`), so the geometry factory still bakes an
   IDENTICAL :math:`\tau` for the sweep path while the carve de-risks by
   parallel-run-and-compare.  The producer-equivalence gate (Leg 1)
   ``tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py`` pins
   the closure-produced :math:`\tau` to (a) the geometry-factory value
   (0-ULP) AND (b) the structurally-independent reference
   :func:`~orpheus.derivations.discrete.sn.contamination.morel_montry_weights`
   (a different code path to the same BMC-2010-Eq.43 weight).  The
   Cartesian :class:`~orpheus.sn.sweep.pole_angular_closure.IdentityAngularClosure`
   supplies the neutral :math:`\tau = 1` WITHOUT a geometry branch — the
   closure TYPE is the dispatch.

   Step B (a later dispatch) retires the geometry-side
   :math:`\tau` producer and consolidates the four-site
   ``c_in``/``c_out`` duplication onto the closure.

.. _sn-closure-c-constants-owned:

c_in / c_out are angular-closure constants — Step B1 (one site folded)
-------------------------------------------------------------------------

.. todo:: Archivist expansion needed.

   The weighted-diamond constants

   .. math::

      c_{\rm out}[m] &= \frac{\alpha_{m+1/2}}{\tau_m}, \\
      c_{\rm in}[m]  &= \frac{1-\tau_m}{\tau_m}\,\alpha_{m+1/2}
                        + \alpha_{m-1/2}

   are an ANGULAR-closure property: a function of the closure's own
   :math:`\alpha`-dome and :math:`\tau` weight :eq:`mm-weights`
   :eq:`dd-mm-closure-constants`.  Issue #236 Phase 2 Step B consolidates
   the FOUR independent inline rebuilds of this pair onto the closure,
   which already computes it once at construction (per :math:`\mu`-level,
   :math:`(M_p,)` arrays in ``_c_in_per_level`` / ``_c_out_per_level``).

   Step B1 (this dispatch) folds the ONE free seam — the
   :class:`~orpheus.sn.sweep.cache.GeometryCoefficients` populator
   (:meth:`~orpheus.sn.sweep.cache.GeometryCoefficients.from_mesh_and_quad`),
   which already holds ``sn_mesh`` and so reads
   :attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.c_out_per_ordinate`
   /
   :attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.c_in_per_ordinate`
   with zero plumbing.  The accessor pair is PUBLIC and polymorphic on the
   base
   :class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`:
   :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
   returns its precomputed per-level :math:`c` gathered to the
   :math:`(N,)` global-ordinate order; the Cartesian
   :class:`~orpheus.sn.sweep.pole_angular_closure.IdentityAngularClosure`
   returns the NEUTRAL zeros (:math:`\alpha=0,\ \tau=1 \Rightarrow c=0`).
   The dispatch is by closure TYPE, not by a ``coord ==`` branch in the
   cache.

   Step B1 is BIT-IDENTICAL: the closure computes :math:`c` from
   closure-:math:`\tau` (0-ULP equal to geometry-:math:`\tau`, pinned by
   the Step-A Leg-1 gate) and the SAME :math:`\alpha` the populator read,
   so the closure's per-level :math:`c` equals the inline :math:`c`
   bit-for-bit; the per-level :math:`\to (N,)` gather is a pure
   permutation (no arithmetic).  The anchor gate
   ``tests/sn/sweep/core/test_cache.py::test_cache_populator_matches_cell_balance_terms``
   pins the cache ``denom`` (which carries :math:`(\Delta A/w)\,c_{\rm out}`)
   to :func:`~orpheus.transport.spatial.cell_balance.cell_balance_terms` at
   ``rtol=1e-14``, and the curvilinear regression snapshots stay unmoved.

   The remaining THREE inline ``c`` rebuild sites (they need CellVisit
   threading) are later dispatches (Step B2 / B3 / C).  See
   :mod:`orpheus.sn.sweep.pole_angular_closure` for the canonical
   accessor and :mod:`orpheus.sn.sweep.cache` for the folded
   consumer.

.. _sn-closure-c-on-cellvisit:

c_in / c_out reach the stateless DD scheme as CellVisit data — Step B2
-------------------------------------------------------------------------

.. todo:: Archivist expansion needed.

   Step B2 folds the SECOND of the four inline ``c_in`` / ``c_out``
   rebuild sites — the matvec-twin residual
   :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.residual`
   (formerly rebuilding :math:`c_{\rm out} = \alpha_{\rm out}/\tau`,
   :math:`c_{\rm in} = (1-\tau)/\tau\,\alpha_{\rm out} + \alpha_{\rm in}`
   inline from the geometry-owned :class:`StreamingTerms`).

   The architectural crux: :class:`~orpheus.transport.spatial.diamond.DiamondDifference`
   is deliberately STATELESS — it reads only the
   :class:`~orpheus.transport.spatial.scheme.CellVisit` packet + the
   :class:`~orpheus.transport.spatial.scheme.UpstreamState`, never the mesh or
   the angular closure.  So the closure-owned :math:`c` cannot reach
   ``DD.residual`` by coupling DD to the closure object (that would break
   the spatial :math:`\otimes` angular separation — the SPATIAL scheme
   must not see the ANGULAR closure's type).  Instead the constants travel
   as DATA: the :class:`~orpheus.transport.spatial.scheme.CellVisit` gains two
   angular-closure-owned fields
   (:attr:`~orpheus.transport.spatial.scheme.CellVisit.c_in` /
   :attr:`~orpheus.transport.spatial.scheme.CellVisit.c_out`, distinct in
   provenance from the geometry-owned
   :attr:`~orpheus.transport.spatial.scheme.CellVisit.streaming_terms`), and the
   single production site
   :meth:`~orpheus.sn.mesh.augmented_mesh.SNMesh._make_cell_visit` — through which
   ALL four ``dag_walk`` yield paths funnel (Pattern 2, no per-site
   divergence) — stamps them from
   :attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.c_in_per_ordinate`
   /
   :attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.c_out_per_ordinate`
   indexed by the GLOBAL ordinate (``direction_idx`` for slab / sphere,
   ``level_indices[p][m]`` for cylinder — mirroring
   :meth:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator.streaming_terms`).
   ``DD.residual`` then reads ``visit.c_in`` / ``visit.c_out``; the
   :math:`(\Delta A/w)`-scaled assembly that follows is byte-unchanged —
   only the SOURCE of :math:`c` moved.

   Step B2 also completes the matvec's typed-consumer binding (Issue
   #226): the unified SN matvec reads ``sn_mesh.pole_angular_closure``
   typed against the
   :class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`
   ABC and drives the angular path through
   :meth:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.precompute_psi_state`,
   :meth:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.cell_contribution`,
   and
   :meth:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.angular_adjoint`.
   These were declared ``@abstractmethod`` on the ABC (matching the
   precedent where
   :class:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase` declares
   ``update`` / ``residual`` abstract so ``mesh.scheme`` consumers see the
   full contract) — making the ABC the COMPLETE strategy contract instead
   of declaring only ``__call__``.

   Step B2 is BIT-IDENTICAL: ``visit.c_in`` / ``visit.c_out``
   (closure-sourced) equal the former inline values 0-ULP — the closure
   computes :math:`c` from closure-:math:`\tau` (0-ULP equal to
   geometry-:math:`\tau`, Step-A Leg-1) and the SAME geometry-:math:`\alpha`,
   and the per-level :math:`\to (N,)` gather is a pure permutation.  The
   matvec residual path (fed by ``DD.residual``) stays bit-for-bit on the
   ``tests/sn/sweep/curvilinear/test_unified_matvec_{sphere,cylinder}.py``
   twin and on the DriftWarning-escalating
   ``tests/sn/sweep/core`` + ``tests/sn/solve`` snapshots.  The remaining
   TWO ``c`` rebuild sites
   (:func:`~orpheus.transport.spatial.cell_balance.cell_balance_terms` for the
   ``DD.update`` solve path; the geometry-side :math:`\tau` producer) are
   Step B3 / C.  See :mod:`orpheus.transport.spatial.scheme` for the CellVisit
   fields and :mod:`orpheus.sn.mesh.augmented_mesh` for the production stamp.

   B2 review fixes (finishing pass).  THREE follow-ups landed after the
   carve, all bit-identical (0-ULP):

   * **Per-ordinate gather cached (L16).** The public accessors
     :attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.c_in_per_ordinate`
     /
     :attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.c_out_per_ordinate`
     re-ran the full :math:`(N,)` per-level :math:`\to` global gather on
     EVERY access, so the per-visit stamp made the visit-producing loop
     :math:`O(N^2\,n_x)`.  The gather is a pure permutation of immutable
     per-level data, so it is now computed ONCE in each mesh-bound
     ``__init__`` (shared
     :meth:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase._build_per_ordinate_cache`,
     called by both ``MorelMontryAngularSweep`` and
     ``IdentityAngularClosure``) and the accessors return the read-only
     cache (``setflags(write=False)`` guards the shared :math:`(N,)` view
     B1's ``GeometryCoefficients`` populator holds).  Measured on a
     ``sphere N=32 nx=200`` walk: :math:`\sim 32\,\text{ms} \to \sim
     22\,\text{ms}` per sweep (:math:`\sim 1.46\times`), value-identical.

   * **Committed production-stamp catcher (vv L11 Mode 11).** The original
     carve had NO committed test exercising ``_make_cell_visit``'s c-stamp
     — the matvec twin reads the closure's ``cell_contribution`` directly
     (never ``DD.residual``), and the diamond / cell-balance fixtures stamp
     visits with a SURROGATE.  A wrong global-ordinate map (a
     ``c_in``:math:`\leftrightarrow`\ ``c_out`` swap, a mis-scattered
     cylinder level block) would ship silently.
     ``tests/sn/sweep/core/test_cell_visit_c_stamp.py`` walks a REAL
     production ``dag_walk`` (sphere + multi-level cylinder + slab) and
     asserts every ``visit.c_in`` / ``visit.c_out`` equals the constants
     recomputed INLINE from that visit's OWN
     :class:`~orpheus.geometry.reduced_operator.StreamingTerms` at 0-ULP
     (the hand-transcribed independent reference, not the closure's own
     ``c`` — vv L11).  Mutation-verified: the ``c_in``\ /\ ``c_out`` swap
     reddens the sphere + cylinder cases.

   * **Test-surrogate dedup (Pattern 2).** The byte-identical
     ``_c_from_streaming_terms`` (``test_diamond.py``) and ``_visit_c``
     (``test_cell_balance_for_streaming.py``) hand-recomputes were unified
     into one shared ``tests/sn/sweep/core/_c_surrogate.py`` consumed by
     both files and the new catcher.

.. _sn-tau-c-on-cellvisit-live:

The live sweep + scan consume closure-owned τ / c — Step B3
---------------------------------------------------------------

Step B1 folded the one free seam (the cache populator); Step B2 carried
the redistribution constants :math:`c_{\rm in}` / :math:`c_{\rm out}`
onto the :class:`~orpheus.transport.spatial.scheme.CellVisit` so the
*apply-direction* residual could read them as data.  Step B3 is the
step that makes the **live** paths consume the closure-owned weight:
the per-cell sweep solve, the matvec solve, and the CumprodScan
fast path now all read the angular weight :math:`\tau` :eq:`mm-weights`
(and the derived constants :eq:`dd-mm-closure-constants`) off the
closure rather than off the geometry-owned ``StreamingTerms.tau_mm``.
After B3 there is **no live reader** of ``StreamingTerms.tau_mm``
anywhere in the sweep, scan, or matvec — precisely the precondition
that let Step C delete the geometry-side :math:`\tau` producer (the two
parallel producers could be reduced to one only once nothing live
depended on the soon-to-be-deleted one; that retirement has now landed
— see the close-out at the end of this section).

This is the third of the four c-fold sites (after the cache populator
in B1 and the residual twin in B2) and the fifth :math:`\tau` consumer.
Like its predecessors it is **bit-identical (0-ULP)**: the carve moves
the *source* of an already-correct number, it does not change the
number.  The sections below derive why :math:`\tau` belongs to the
angular closure, why the constants must travel as visit *data* rather
than as a closure *reference*, why the field default is :math:`1.0`
(and not the more obvious :math:`0.0`), how the three live consumers
share one operator, and what makes the fold provably bit-identical and
therefore a safe regression floor for Step C.

τ is an angular-scheme property the closure owns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Morel--Montry weight

.. math::

   \tau_n = \frac{\mu_n - \mu_{n-\frac12}}{\mu_{n+\frac12} - \mu_{n-\frac12}}

(:eq:`mm-weights`, Bailey--Morel--Chang 2010 Eq. 43) is built **entirely**
from the angular quadrature: the ordinate cosine :math:`\mu_n`, the
neighbouring angular-cell edges :math:`\mu_{n\pm 1/2}`, and — for the
cylinder — the :math:`\mu`-level partition that groups ordinates.  Not
one of those inputs is a property of the *spatial* streaming geometry:
:math:`\tau` does not depend on the cell volume :math:`V_i`, the face
areas :math:`A_{i\pm 1/2}`, the surface-curvature redistribution area
:math:`\Delta A_i`, or the radial mesh at all.  It is a number attached
to an **ordinate**, not to a **cell**.

That :math:`\tau` had historically lived on
:class:`~orpheus.geometry.reduced_operator.StreamingTerms` (as
``tau_mm``) was an accident of where the curvilinear sweep was first
assembled — the streaming-geometry factory happened to be the object
in scope when the weighted-diamond closure was wired in, so it baked
the angular weight in alongside the genuinely geometric
:math:`\alpha`-dome and face areas.  The architectural correction
(Issue #236 Phase 2) is to give the weight back to its owner: the
**pole-angular closure**
(:class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`),
which already binds the quadrature and already computes :math:`\tau`
from it in :func:`~orpheus.sn.sweep.pole_angular_closure.morel_montry_tau_per_level`.
The closure exposes it through one public, polymorphic accessor,
:attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.tau_per_ordinate`,
a read-only :math:`(N,)` array indexed by the global ordinate.  The
two concrete strategies answer it differently *by type*, with no
``coord ==`` branch anywhere:

* :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
  returns its own per-level :math:`\tau` gathered to the global order;
* :class:`~orpheus.sn.sweep.pole_angular_closure.IdentityAngularClosure`
  (Cartesian) returns the neutral :math:`\tau \equiv 1` — there is no
  angular redistribution in slab geometry, so the M-M weight reduces to
  its identity element (see below).

This is the same both-sites mint B1 applied to the :math:`c`-accessors:
:attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.tau_per_ordinate`
is declared on the
:class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`
ABC, so the contract is complete for every consumer.  (Phase B
originally declared these accessors on **both** the
``@runtime_checkable`` ``PoleAngularClosure`` Protocol and the ABC, to
serve structural-typing and nominal-inheritance consumers alike;
Issue #236 Phase 2 B2 retyped every consumer onto the ABC and Issue
#248 deleted the now-orphaned Protocol, so the ABC is the single
declaration site.)  The gather itself is a pure permutation
of the immutable per-level data, hoisted once into each mesh-bound
``__init__`` via the shared
:meth:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase._build_per_ordinate_cache`
(renamed from ``_build_c_per_ordinate_cache`` now that it gathers three
constants — :math:`c_{\rm in}`, :math:`c_{\rm out}`, and :math:`\tau`
— rather than two); the accessor returns the cached read-only view, so
the per-visit lookup is :math:`O(1)`.

The spatial ⊗ angular separation forbids coupling DD to the closure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If :math:`\tau` is owned by the angular closure but consumed in the
spatial cell update, the obvious move would be to hand the closure
object to the cell-update scheme so it can ask for
:math:`\tau`.  That move is **forbidden by design**, and the reason is
the load-bearing architectural fact of the SN sweep.

:class:`~orpheus.transport.spatial.diamond.DiamondDifference` is a **stateless
spatial discretization scheme**.  It reads only the per-cell
:class:`~orpheus.transport.spatial.scheme.CellVisit` packet and the
sweep-resolved :class:`~orpheus.transport.spatial.scheme.UpstreamState`; it
never sees the mesh, the quadrature, or the angular closure.  The whole
point of the spatial :math:`\otimes` angular product is that the
spatial scheme is interchangeable (diamond difference, linear
discontinuous, ...) without knowing *which* angular treatment sits on
the other axis of the tensor product, and the angular closure is
interchangeable (Morel--Montry, identity, a future Carlson variant)
without knowing the spatial scheme.  Coupling
:class:`~orpheus.transport.spatial.diamond.DiamondDifference` to
:class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`
would collapse that product into a Cartesian-vs-curvilinear conditional
inside the spatial scheme — exactly the geometry dispatch the unified
body was built to delete.

So the constants travel as **data**, not as a **dependency**.  The
:class:`~orpheus.transport.spatial.scheme.CellVisit` packet — which the
orchestrator already populates per cell and per ordinate — carries the
angular-closure-owned numbers as plain ``float`` fields:
:attr:`~orpheus.transport.spatial.scheme.CellVisit.c_in` and
:attr:`~orpheus.transport.spatial.scheme.CellVisit.c_out` (added in B2) and now
:attr:`~orpheus.transport.spatial.scheme.CellVisit.tau` (B3).  They are stamped
at exactly **one** production site,
:meth:`~orpheus.sn.mesh.augmented_mesh.SNMesh._make_cell_visit`, through which all
four ``dag_walk`` yield paths (slab, sphere, cylinder, cylindrical
pure-azimuthal degenerate) funnel — Pattern 2, no per-site divergence.
That site reads the closure's per-global-ordinate accessors and stamps:

.. code-block:: python

   closure = self.pole_angular_closure
   return CellVisit(
       cell_idx=cell_idx,
       streaming_terms=st,
       face_area_downstream=face_area_downstream,
       c_in=float(closure.c_in_per_ordinate[global_ordinate]),
       c_out=float(closure.c_out_per_ordinate[global_ordinate]),
       tau=float(closure.tau_per_ordinate[global_ordinate]),
   )

where ``global_ordinate`` is the global ordinate index resolved the
same way :meth:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator.streaming_terms`
resolves it (``direction_idx`` for slab / sphere,
``level_indices[mu_level_idx][m]`` for cylinder).  The spatial scheme
downstream sees only ``visit.tau`` / ``visit.c_in`` / ``visit.c_out``;
it has no idea a closure produced them.  The provenance is recorded in
the field docstrings (the constants are distinct in origin from the
geometry-owned :attr:`~orpheus.transport.spatial.scheme.CellVisit.streaming_terms`),
but the *type system* never lets the spatial axis reach across to the
angular axis.

Why the CellVisit.tau default is 1.0, not 0.0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:attr:`~orpheus.transport.spatial.scheme.CellVisit.tau` defaults to
:math:`1.0`, not the zero a numeric field usually defaults to.  The
default is the value the slab / Cartesian path leaves in place — the
identity closure supplies it through ``tau_per_ordinate`` — and the
choice is forced by the angular recurrence the value feeds.

:math:`\tau = 1` is the **neutral element** of the Morel--Montry weight.
With :math:`\tau_n = 1` the WDD angular closure
:eq:`wdd-closure` becomes
:math:`\psi_{n,i} = \tau_n\,\psi_{n+\frac12} + (1-\tau_n)\,\psi_{n-\frac12}
= \psi_{n+\frac12}`, i.e. the step (fully-outgoing) closure, and the
outgoing-face recurrence

.. math::

   \psi^a_{\rm out} = \frac{\bar\psi - (1-\tau)\,\psi^a_{\rm in}}{\tau}
   \;\xrightarrow{\;\tau = 1\;}\;
   \frac{\bar\psi - 0}{1} = \bar\psi

reduces to the **identity** in :math:`\bar\psi` — exactly what slab
geometry needs, where there is no angular redistribution and the
"angular-out" state is just the cell average.  Likewise the
denominator constant :math:`c_{\rm out} = \alpha_{n+1/2}/\tau` and the
scan split :math:`1/\tau`, :math:`(1-\tau)/\tau` are all well-defined
and reduce to their slab values (:math:`0`, :math:`1`, :math:`0`
respectively, since :math:`\alpha = 0` on the slab).

A :math:`0.0` default, by contrast, is a **divide-by-zero landmine**.
Every consumer of :math:`\tau` divides by it:
:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.update`'s angular
thread divides :math:`(\bar\psi - (1-\tau)\psi^a_{\rm in})` by
:math:`\tau`; the cache derives ``tau_inv = 1/tau``.  A visit that
reached the curvilinear branch with an un-stamped :math:`\tau = 0`
would produce a silent ``inf``/``nan`` rather than a loud error.  The
:math:`1.0` default makes the *un-stamped* visit compute the
**identity** transformation — the safe no-op — which is both correct
for the slab (where the stamp is genuinely neutral) and a benign
fallback that surfaces a missing stamp as a *wrong-but-finite* answer a
regression snapshot catches, rather than a ``nan`` that propagates.

The three live consumers and the L21 framing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Three live paths read the closure-owned :math:`\tau` (or the constants
derived from it) after B3.  The first two are the **solve** and
**apply** directions of the same per-cell linear system; the third is
the vectorized scan form of the same recurrence — the apply / sweep
duality this page calls the **L21 twin-path** (two applications of the
*same* operator; cf. :ref:`phase-c-sweep-frame-matvec` for the
apply-direction matvec that is the twin of the curvilinear sweep).

**(1) The scalar solve helper.**
:func:`~orpheus.transport.spatial.cell_balance.cell_balance_terms` — the
:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.update` solve
direction — no longer rebuilds :math:`c_{\rm in}` / :math:`c_{\rm out}`
from ``st.alpha_* / st.tau_mm``.  Its signature now takes them as
keyword inputs, and
:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.update` supplies
them straight off the visit::

   terms = cell_balance_terms(
       visit.streaming_terms, visit.face_area_downstream, total_xs,
       upstream_state, c_in=visit.c_in, c_out=visit.c_out,
   )

The helper now reads :math:`\tau` **not at all** — it consumes only the
already-derived :math:`c` constants, which is the right factoring: the
cell-balance denominator :eq:`dd-solve` needs
:math:`(\Delta A_i / w_n)\,c_{\rm out}`, and the upstream numerator
needs :math:`(\Delta A_i / w_n)\,c_{\rm in}\,\psi_{n-\frac12}`, neither
of which references :math:`\tau` once :math:`c` is in hand.

**(2) The angular recurrence.** The other half of
:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.update` — the
Morel--Montry outgoing-angular-face thread — *does* need the raw
:math:`\tau`:

.. math::
   :label: dd-mm-angular-recurrence

   \psi^a_{\rm out}
   = \frac{\bar\psi - (1 - \tau)\,\psi^a_{\rm in}}{\tau}

and reads it from :attr:`~orpheus.transport.spatial.scheme.CellVisit.tau`
(stamped by
:meth:`~orpheus.sn.mesh.augmented_mesh.SNMesh._make_cell_visit` from
:attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.tau_per_ordinate`)
rather than from ``visit.streaming_terms.tau_mm``.  This is the line
the :math:`1.0` default protects.

**(3) The CumprodScan split.** The vectorized fast path — the cumulative
product that replaces the per-cell Python loop for the curvilinear
sweep — needs the same recurrence in a form amenable to a forward scan.
:class:`~orpheus.sn.sweep.cache.GeometryCoefficients` sources
:math:`\tau` from
:attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.tau_per_ordinate`
and precomputes the split

.. math::
   :label: dd-mm-scan-split

   \texttt{tau\_inv} = \frac{1}{\tau},
   \qquad
   \texttt{mm\_a\_in\_coeff} = \frac{1 - \tau}{\tau},

consumed at the loss-representation scan recurrence (in
:mod:`orpheus.sn.loss_representation`) as

.. math::

   \psi^a_{\rm out}
   = \texttt{tau\_inv}\cdot\bar\psi
     - \texttt{mm\_a\_in\_coeff}\cdot\psi^a_{\rm in}
   = \frac{\bar\psi}{\tau} - \frac{1-\tau}{\tau}\,\psi^a_{\rm in},

which is algebraically identical to :eq:`dd-mm-angular-recurrence` — the
same operator, applied in the vectorized frame.  Precomputing
:math:`1/\tau` and :math:`(1-\tau)/\tau` is a legitimate
perform-once-at-construction hoist (L16): the closure exposes only the
**primitive** :math:`\tau` (Pattern 5 — build the primitive, not the
product), and each consumer derives the trivial :math:`1/\tau`,
:math:`(1-\tau)/\tau`, :math:`\alpha_{\rm out}/\tau` algebra it needs at
its own definition site.  The scan derivation lives in the cache; the
recurrence consumes it.  That keeps the closure's surface minimal (one
weight) while letting two structurally-different consumers (a scalar
Python recurrence and a vectorized NumPy scan) each shape it to their
reduction tree.

Why the fold is bit-identical, and the regression floor for Step C
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Step B3 is **bit-identical (0-ULP)**.  The argument has two legs, both
already established by the earlier carve steps:

#. **Closure-:math:`\tau` is 0-ULP equal to geometry-:math:`\tau`.** The
   closure's
   :func:`~orpheus.sn.sweep.pole_angular_closure.morel_montry_tau_per_level`
   is a line-for-line replica of the geometry factory's :math:`\tau`
   arithmetic (Step A), pinned by the **Leg-1 producer-equivalence
   gate** ``tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py``
   to the geometry-factory value (0-ULP) *and* to the
   structurally-independent reference
   :func:`~orpheus.derivations.discrete.sn.contamination.morel_montry_weights`
   (a different code path to the same BMC-2010-Eq. 43 weight).  Reading
   :math:`\tau` from the closure therefore yields exactly the same
   ``float64`` bits the former ``st.tau_mm`` read carried.

#. **The per-level :math:`\to (N,)` gather is a pure permutation.** No
   arithmetic happens between the closure's per-level :math:`\tau` and
   the global-ordinate :math:`(N,)` view the accessor returns — only a
   reindex.  So every derived quantity (:math:`c_{\rm in}`,
   :math:`c_{\rm out}`, ``tau_inv``, ``mm_a_in_coeff``) is bit-for-bit
   what the inline rebuilds produced.

This was confirmed not only by the regression snapshots but **in
process**: at sphere (single :math:`\mu`-level) and cylinder
(multi-level) configurations,
:math:`\lvert\texttt{visit.tau} - \texttt{st.tau\_mm}\rvert`,
:math:`\lvert\texttt{closure.tau\_per\_ordinate} - \texttt{st.tau\_mm}\rvert`,
and ``np.array_equal`` on the cache's ``tau_inv`` / ``mm_a_in_coeff``
against the geometry-:math:`\tau`-derived split were **all exactly
zero**.  The DriftWarning-escalating ``tests/sn/sweep/core``,
``tests/sn/sweep``, and ``tests/sn/solve`` snapshots stayed unmoved
(588 + 60 green), with zero drift escalation.

The bit-identity guarantee is what makes B3 the **regression floor for
Step C**.  Two committed catchers pin the new live paths so that the
retirement of the parallel geometry-:math:`\tau` producer (now landed,
Step C) cannot silently break them:

* The **Leg-1 producer-equivalence gate** pins
  ``closure.tau_per_ordinate`` to the BMC-2010-Eq. 43 reference, so the
  closure remains the correct sole producer after the geometry one is
  deleted.

* The **production-stamp catcher**
  ``tests/sn/sweep/core/test_cell_visit_c_stamp.py`` walks a real
  production ``dag_walk`` (sphere, multi-level cylinder, slab) and — in
  its dedicated :math:`\tau` arm — asserts every ``visit.tau`` equals
  the **independently recomputed** Morel--Montry weight for that visit's
  ordinate at 0-ULP.  The reference is
  :func:`~orpheus.derivations.discrete.sn.contamination.morel_montry_weights`
  (a structurally-independent code path to the same BMC-2010-Eq. 43
  weight, not the closure comparing against itself — vv L11): the test
  pins the stamp's *ordinate map*, the complement of what Leg-1 pins
  (the producers' *values*).  Before Step C this catcher used the
  *geometry-produced* ``st.tau_mm`` as its independent reference; when
  Step C deleted that field the oracle was re-pointed onto
  ``morel_montry_weights`` (with the cylinder clamp replicated), keeping
  it geometry-:math:`\tau`-free and still structurally independent of
  the closure under test.  This arm was added specifically because B3
  made the
  :attr:`~orpheus.transport.spatial.scheme.CellVisit.tau` stamp **live** while
  the existing named twins never call the rewired reader (vv L11
  Mode 11): a mutation stamping ``tau = ... * 1.1`` drifts the converged
  cylinder scalar flux by :math:`\sim 0.2\,\%` with **no** other test
  red, so the dedicated arm is the only committed catcher of a
  :math:`\tau`-stamp ordinate-map error.  Mutation-verified RED on the
  :math:`\times 1.1` stamp across sphere + cylinder + slab; GREEN clean.

* The **seam-6 scan catcher**
  ``tests/sn/sweep/core/test_affine_carve_baseline.py`` reddens on a
  corruption of the CumprodScan :math:`\tau` split, pinning the third
  live consumer.

With these in place, **Step C has now deleted** the geometry-side
:math:`\tau` producer.  The :math:`\tau` blocks inside
:func:`~orpheus.geometry.reduced_operator.spherical_streaming`
(the ``mu_edge`` weight-sum loop)
:math:`\,/\,`
:func:`~orpheus.geometry.reduced_operator.cylindrical_streaming` (the
per-level ``eta_edge`` loop) and the slab synthetic were excised, and
the now-orphaned ``StreamingTerms.tau_mm``, ``StreamingTerms.alpha_in``
/ ``alpha_out`` (whose sole readers were the c-rebuild sites B1--B3 just
retired), ``ReducedStreamingOperator.tau_mm``, and
``ReducedStreamingOperator.tau_mm_per_level`` dataclass fields were
dropped — confident that nothing live depended on them.
See :mod:`orpheus.sn.sweep.pole_angular_closure` for the
``tau_per_ordinate`` accessor and the three-constant cache,
:mod:`orpheus.transport.spatial.scheme` for the
:attr:`~orpheus.transport.spatial.scheme.CellVisit.tau` field,
:mod:`orpheus.sn.mesh.augmented_mesh` for the single production stamp, and
:mod:`orpheus.sn.sweep.cache` for the scan split.

.. _sn-tau-step-c-closeout:

Step C close-out — the geometry-side τ producer is retired
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The retirement is a **surgical field excision plus a test-oracle
migration**, not a blind delete: the τ producer was already dead in
*production* (B3 left zero live readers), but it remained load-bearing
as a *test oracle* — the regression-floor catchers used the
geometry-produced ``st.tau_mm`` as their structurally-independent
reference.  Migrate-then-delete preserved the floor:

#. **The oracles were re-pointed first.** The production-stamp catcher
   and the surviving producer-equivalence legs now read
   :func:`~orpheus.derivations.discrete.sn.contamination.morel_montry_weights`
   — a different code path to the same BMC-2010-Eq. 43 weight, unclamped
   on both geometries, with the cylinder clamp :math:`\mathrm{clip}(\cdot,
   \tfrac12, 1)` replicated in the test surrogate.  This kept the floor
   green *while geometry-:math:`\tau` was still present*, proving the
   migration faithful before any deletion.  The two
   ``*_equals_geometry_factory_0ulp`` legs of
   ``test_tau_producer_equivalence.py`` (which compared closure-:math:`\tau`
   against the soon-to-be-deleted factory output) became vacuous and were
   retired; the independent-reference and clamp-difference legs survive.

#. **The producers were excised surgically.** The τ blocks are
   *interleaved* with still-live outputs:
   :func:`~orpheus.geometry.reduced_operator.spherical_streaming` shares
   its ``mu_edge`` array with the live ``mu_start`` (the Hébert §3.9.4
   starting direction :math:`\mu_{1/2} = -1.0`), and
   :func:`~orpheus.geometry.reduced_operator.cylindrical_streaming`
   shares its per-level loop with the live ``mu_start_per_level``.  A
   whole-function deletion would have been wrong; only the τ statements
   were removed.  The :math:`\alpha`-dome (``alpha_half`` /
   ``alpha_per_level``), the redistribution factor ``redist_dAw``, the
   face areas, and the starting-direction edges all **stay on the
   geometry operator** — they are genuinely geometric.

#. **The deletion was proven inert.** The bit-identity regression gates
   (run under an escalated ``DriftWarning``) showed **zero** failures
   across the sweep / scan / matvec suites, and the test-count delta
   reconciled exactly to the four retired ``*_equals_geometry_factory``
   legs and the two retired ``test_reduced_operator.py``
   τ-bit-identical tests — no silent test loss.  After deletion the
   re-pointed catcher was **mutation-verified RED** (a :math:`\times 1.1`
   stamp and a ``c_in`` :math:`\leftrightarrow` ``c_out`` swap both
   reddened it against the independent oracle), confirming the migrated
   catcher is a real catcher reading the independent reference, not a
   tautology against the closure.

.. note::

   The legacy ``__call__``-argument ``tau_mm`` on the unbound
   :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
   path (``MorelMontryAngularSweep(sn_mesh=None)``, where :math:`\tau` was
   passed as a runtime argument because the closure is not mesh-bound) was
   a **separate surface** that **survived Step C** unchanged — it was the
   closure's own runtime parameter, not the geometry-side field the carve
   retired.  It was subsequently retired under
   `Issue #248 <https://github.com/deOliveira-R/ORPHEUS/issues/248>`_
   (landed in this same re-staging, which deleted the strategy
   ``__call__`` bundle entirely; see the contract-evolution note at
   :ref:`sn-pole-angular-closure-protocol`).

.. _sn-space-angle-separability-section:

Space ⊗ angle separability — the (spatial ⊗ angular) product capstone (Issue #236 Phase 3)
-----------------------------------------------------------------------------------------------

This section closes the Issue #236 *(spatial* :math:`\otimes` *angular)
product* narrative on the theory page.  The campaign had three phases:

* **Phase 1 — pairing validity.**  The spatial closure (the diamond /
  weighted-diamond cell update of :eq:`dd-curvilinear-scalar`) and the
  angular closure (the Morel--Montry weight :math:`\tau` of
  :eq:`mm-weights`, the redistribution dome
  :eq:`alpha-recursion`) are two distinct, independently-selectable
  axes — a genuine tensor product, with separate injection points on
  :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` (``cell_update=`` for the spatial
  scheme, ``pole_angular_closure=`` for the angular scheme).
* **Phase 2 — :math:`\tau`-ownership carve.**  The angular weight
  :math:`\tau` was moved off the geometry operator and onto the angular
  closure, so the angular axis literally *owns* its own discretisation
  knob (:ref:`sn-tau-closure-owned` through
  :ref:`sn-tau-step-c-closeout`).  That carve made the product
  *structural in the type system*, not merely conceptual.
* **Phase 3 — separability characterisation (this section).**  Having
  established the product and given each axis its own knob, the final
  question is: *how do the two error contributions combine?*  The answer
  is geometry-dependent, and it is the campaign's headline claim.

The decomposition is pinned permanently by the L1 MMS characterisation
gate :mod:`tests.sn.verification.mms.test_space_angle_separability`,
which carries ``@pytest.mark.verifies("sn-space-angle-separability")``
against :eq:`sn-space-angle-separability` below.

The space–angle error decomposition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Write the total SN discretisation error of the scalar flux as a function
of the two refinement parameters — the spatial mesh size
:math:`h \sim 1/n_{\rm cells}` and the angular (quadrature) order
:math:`N` (the ordinate count, or the azimuthal order :math:`n_\varphi`
in the cylinder).  Let :math:`E_{\rm space}(h)` be the error of the
spatial closure at infinite quadrature order and :math:`E_{\rm angle}(N)`
the error of the angular closure on an exactly-resolved spatial mesh.
The campaign's headline result is the geometry-split law

.. math::
   :label: sn-space-angle-separability

   E(h, N) \;\approx\;
   \begin{cases}
     E_{\rm space}(h) \,+\, E_{\rm angle}(N),
       & \text{Cartesian (slab / 2-D / 3-D): \textbf{separates}},\\[1.2ex]
     \max\!\bigl(E_{\rm space}(h),\, E_{\rm angle}(N)\bigr),
       & \text{curvilinear (sphere / cylinder): \textbf{gates}}.
   \end{cases}

The two regimes are distinguished operationally by the **mixed second
difference** — the discrete :math:`\partial^2 E / \partial h\,\partial N`
evaluated on a two-quadrature error table over a coarse/fine pair of
mesh sizes:

.. math::
   :label: sn-space-angle-cross-term

   M \;=\; E[h_1, N_1] - E[h_1, N_2] - E[h_2, N_1] + E[h_2, N_2],
   \qquad
   \frac{|M|}{\max(\Delta E_h,\, \Delta E_N)} \;
   \begin{cases}
     \ll 1, & \text{separable (additive),}\\
     = \mathcal{O}(1), & \text{gated (coupled).}
   \end{cases}

For an additively separable error, :math:`E[h,N] = f(h) + g(N)` exactly,
so the cross-term telescopes to zero: :math:`M = f(h_1) + g(N_1) -
f(h_1) - g(N_2) - f(h_2) - g(N_1) + f(h_2) + g(N_2) = 0`.  A non-zero
:math:`M` is therefore a *direct, mechanism-anchored* measurement that
the two axes interact — that the second mixed partial of the error
surface does not vanish.  This is the quantity the ST5 gate measures.

.. (vv-status rationale) Characterisation claim, now tested: both the
   law :eq:`sn-space-angle-separability` and its discriminator
   :eq:`sn-space-angle-cross-term` describe the STRUCTURE of the
   discretisation-error surface (the regime discrimination), not a
   solver eigenvalue or flux VALUE.  This is an L1 MMS-convergence-
   structure (math) claim per vv-principles — MMS does not reach the
   eigenvalue layer, so neither label is, or ever becomes, an
   eigenvalue / flux-value claim.  The ST5 characterisation gate
   ``test_space_angle_separability.py`` now carries the
   ``verifies`` markers for both labels, so each is ``documented`` AND
   ``tested``: the verifying ``tests`` edge is the characterisation
   gate (an L1 MMS gate), not a closed-form / semi-analytical value
   reference.  What each gate leg verifies:
     * :eq:`sn-space-angle-separability` (the geometry-split decomposition
       law) is pinned by all six legs as a *positive* signature — the
       Cartesian legs assert separability (mixed-second-difference
       :math:`\to 0`, N-independent spatial rate); the curvilinear legs
       assert gating (N-gated spatial rate).  The marker is FILE-level.
     * :eq:`sn-space-angle-cross-term` (the mixed-second-difference
       discriminator :math:`M`) is the gate's measured instrument: the
       three legs that assert directly on :math:`|M|/\max`
       (``test_cartesian_slab_iso_space_angle_separable``,
       ``test_cartesian_slab_p1_aniso_floor_n_independent``,
       ``test_sphere_cross_term_large_discriminates_from_cartesian``)
       carry the per-test ``verifies`` marker.  The discriminator is a
       quantity the gate *measures against a declared threshold*, not a
       passive derivation step, so the ``tested`` edge is a real
       coverage claim.
   The posture mirrors the pole-cell characterisation gate this gate is
   modelled on (#233).
.. vv-status: sn-space-angle-separability documented
.. vv-status: sn-space-angle-separability tested
.. vv-status: sn-space-angle-cross-term documented
.. vv-status: sn-space-angle-cross-term tested

Why the two axes factorize: LMM-1987 (spatial) × BMC-2010 (angular)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The decomposition is not an empirical accident — it is forced by the
structure of the asymptotic *diffusion-limit-consistency* literature,
which is **literally split into a spatial paper and an angular paper**.
This split is the strongest possible evidence that the consistency
conditions live on two separate axes.

* **The spatial condition** — Larsen, Morel & Miller, "Asymptotic
  solutions of numerical transport problems in optically thick,
  diffusive regimes," *Journal of Computational Physics*
  **69(2):283--324 (1987)**, DOI
  `10.1016/0021-9991(87)90170-7 <https://doi.org/10.1016/0021-9991(87)90170-7>`_
  (and Part II, Larsen & Morel, JCP **83(2):212--236 (1989)**, DOI
  10.1016/0021-9991(89)90229-5).  LMM analyse the *spatial* differencing
  scheme's diffusion limit (cells scaled so they are not optically
  thin): a scheme whose discrete spatial limit is itself a valid
  diffusion discretisation (linear-discontinuous, weighted-diamond with
  the right closure) is "substantially more accurate" than one without
  (bare diamond difference).  This is a condition **on the spatial axis
  alone** — the angular order does not enter.

* **The angular condition** — Bailey, Morel & Chang,
  [BaileyMorelChang2010]_, "The asymptotic diffusion-limit accuracy of
  S\ :sub:`N` angular differencing schemes," *Nuclear Science and
  Engineering* **165(2):149--169 (2010)**, DOI
  `10.13182/NSE08-66 <https://doi.org/10.13182/NSE08-66>`_.  BMC analyse
  the SN equations **discretised only in angle, with space kept
  continuous** (their analysis deliberately removes spatial differencing
  to isolate the angular error).  They prove that the angular axis
  carries its *own* diffusion-limit condition, independent of the
  spatial one.  Their p. 151 statement is the separability fact in the
  authors' own words: the spatial half "has been shown by Larsen, Morel,
  and Miller," while "retaining full first-order consistency can be
  important for **angular** discretisations" — the angular contribution
  they introduce.

The two conditions factorise.  In the leading-order (:math:`\varepsilon^0`)
diffusion limit, *any* weighted-diamond angular weight (step, diamond,
Morel--Montry) preserves consistency — BMC Eqs. (23)--(25).  The
**first-order** (:math:`\varepsilon^1`) limit carries a contamination
term :math:`\beta` (BMC Eq. 40), a **purely angular** functional of the
redistribution coefficients and quadrature,
:math:`\beta = \sum_m \mu_m\bigl[\alpha_{m+1/2}\mu_{m+1/2} -
\alpha_{m-1/2}\mu_{m-1/2}\bigr]`, which vanishes *only* for the
Morel--Montry weights (BMC Eq. 43, the weight of :eq:`mm-weights`).
Because :math:`\beta` depends on no spatial quantity, the angular
condition is provable on its own axis — exactly as the spatial
condition is provable on its own axis.  The diffusion limit needs
**both**:

.. math::

   \text{accurate diffusion limit}
   \;\;\Longleftrightarrow\;\;
   \underbrace{(\text{LMM spatial condition})}_{\text{depends on the spatial scheme only}}
   \;\;\wedge\;\;
   \underbrace{(\text{BMC angular condition},\ \beta = 0)}_{\text{depends on the angular weights only}}.

This conjunction of two single-axis conditions is *why* the Cartesian
error separates additively (each axis contributes its own consistency
defect, and the two defects add) and *why* a bad pairing can still break
the limit (independence of *selection* is not independence of
*consequence* — both conditions must hold simultaneously).

.. note::

   The literature's double-use of the name "linear-discontinuous" (LD)
   is itself evidence of the two-axis structure: LMM and every
   spatial-scheme paper list LD as a *spatial* scheme, while Lathrop
   (2000) lists "linear-discontinuous" among his *angular* differencing
   schemes.  The same trial-space name applies on either axis; the
   ORPHEUS registries disambiguate by axis (a spatial cell-update vs an
   angular closure), never a single ``LD`` enum.  This is the #158
   (spatial scheme) vs #6 (LD *angular* finite elements) distinction.

Cartesian separates, curvilinear gates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The geometry split of :eq:`sn-space-angle-separability` follows
mechanically from *where the angular redistribution term lives*.

**Cartesian — additive separation.**  In slab / 2-D / 3-D Cartesian
geometry the curvilinear angular-redistribution term
:math:`\frac{1-\mu^2}{r}\,\partial_\mu\psi` is **absent**
(:math:`r \to \infty`; there is no :math:`1/r`).  The angular closure is
the :class:`~orpheus.sn.sweep.pole_angular_closure.IdentityAngularClosure`,
which contributes no redistribution: each ordinate's spatial sweep is
fully independent of every other ordinate.  The Cartesian cell update
(:eq:`dd-2d-balance-form` / the slab balance) consumes only the per-axis
streaming ratios and :math:`\Sigma_t` — **no** :math:`\tau`, **no**
angular state.  The spatial error and the angular (quadrature) error are
generated by disjoint mechanisms, so they add:
:math:`E(h,N) \approx E_{\rm space}(h) + E_{\rm angle}(N)`, and the
mixed partial vanishes.  The operational signature is that the spatial
convergence **rate** is the same at every quadrature order
(N-independent O(h\ :sup:`2`)).

**Curvilinear — multiplicative gating.**  In the sphere / cylinder the
Morel--Montry angular thread

.. math::

   \psi_{n+\frac12} \;=\;
       \frac{\overline{\psi}_n - (1-\tau_n)\,\psi_{n-\frac12}}{\tau_n}

couples the ordinates *sequentially within a* :math:`\mu`-*level*, and
the coupling enters the **shared cell-balance denominator** of
:eq:`dd-curvilinear-scalar`: the redistribution divisor
:math:`(\Delta A_i / w_n)\,c_{\rm out}` (with
:math:`c_{\rm out} = \alpha_{\rm out}/\tau_n`) sits in the *same*
denominator that produces the spatial cell average
:math:`\overline{\psi}_{n,i}` the spatial closure then uses.  The
angular interpolation error of the :math:`\tau`-thread therefore
**caps** the accuracy the spatial closure can deliver: at a coarse
quadrature, refining :math:`h` cannot drive the cell average below the
angular floor, because the angular term contaminates the denominator
the spatial refinement acts through.  Hence the error *gates*:
:math:`E(h,N) \approx \max(E_{\rm space}(h), E_{\rm angle}(N))`.  You
cannot harvest fine-:math:`h` accuracy at coarse :math:`N`; both axes
must advance together.  The mechanism is documented in detail at
:ref:`sn-tau-c-on-cellvisit-live` (why :math:`\tau` is an angular
property that nonetheless flows through the spatial denominator) and the
shared-denominator algebra is :eq:`dd-curvilinear-scalar`.

.. warning::

   The gating is a property of *today's* curvilinear closure (the 1-D
   :math:`\eta`-march Morel--Montry thread), not a law of nature.  A
   future 2-D angular closure (#229) that resolves the
   :math:`(\eta,\varphi)` azimuthal variation the 1-D march cannot
   thread, or a higher-order spatial scheme (#158 / #6), would *lift*
   the gating — at which point the curvilinear error would begin to
   separate.  The ST5 gate is designed so that lifting the floor reddens
   the gating assertions (the coarse-N saturated h-ratio rises toward
   the O(h\ :sup:`2`) value), signalling that the regime changed; that
   redding is the intended signal to *re-tune* the gate to the new,
   better regime, **not** a regression.

Measured cross-term evidence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The decomposition was established empirically by the four
``diag_sep_*`` probes (reproduced bit-for-bit after the Phase-2
:math:`\tau` carve) and is now pinned by the ST5 gate.  The measured
mixed-second-difference :math:`|M|/\max` spans **three orders of
magnitude** between the separable Cartesian and the gated sphere — a
clean discrimination band, not a brittle exact number.

.. list-table:: Measured space–angle error structure (2026-06-18,
   ``nx ∈ {20, 40, 80}``)
   :header-rows: 1
   :widths: 18 14 30 24 14

   * - Geometry
     - Regime
     - Scalar-flux L2 ladder (coarse → fine quadrature)
     - Spatial h-ratios (coarse-N / fine-N)
     - :math:`|M|/\max`
   * - Slab (isotropic)
     - **separates**
     - N=4 ``[5.42e-4, 1.35e-4, 3.38e-5]`` · N=16 ``[5.40e-4, 1.35e-4, 3.37e-5]``
     - ``[4.01, 4.00]`` / ``[4.01, 4.00]`` (N-independent O(h²))
     - **0.0047**
   * - Slab (P1 aniso)
     - **separates**
     - N=4 floor ``6.80e-3`` · N=16 floor ``6.79e-3`` (flat, angular floor)
     - flat at both N (floor N-independent to <0.3 %)
     - **0.0038**
   * - Cylinder
     - **gates**
     - :math:`n_\varphi`\ =8 ``[1.95e-2, 1.91e-2, 1.90e-2]`` · :math:`n_\varphi`\ =16 ``[8.05e-3, 7.47e-3, 7.37e-3]``
     - ``[1.02, 1.00]`` (saturated); azimuthal floor drops 2.58× at :math:`n_\varphi` 8→16
     - **0.019** (small only because :math:`E \approx E_{\rm angle}` swamps)
   * - Sphere
     - **gates**
     - N=8 ``[1.47e-2, 5.40e-3, 4.69e-3]`` · N=32 ``[1.50e-2, 3.71e-3, 9.29e-4]``
     - ``[2.71, 1.15]`` (saturates) / ``[4.04, 4.00]`` (O(h²) recovers)
     - **0.411**

The reading of the table:

* **Slab (both rows): separable.**  The spatial h-ratio is :math:`\approx
  4` (O(h\ :sup:`2`)) at *every* quadrature order — the spatial rate is
  blind to :math:`N`.  The isotropic row has a genuine O(h\ :sup:`2`)
  window; the P1-anisotropic row sits at a flat MMS/angular floor that
  is the *same* at every :math:`N`.  Both have :math:`|M|/\max \le
  0.005` — the cross-term vanishes whether or not the angular axis is
  active.  The P1 row is the load-bearing control: separability survives
  an *active* angular term, so it is not an artefact of the isotropic
  degeneracy.

* **Sphere: gating, the discriminator.**  At coarse N=8 the finest
  spatial h-ratio collapses to **1.15** — refinement saturates at the
  angular floor.  At fine N=32 the *same* spatial ladder recovers
  :math:`\approx 4.00` (O(h\ :sup:`2`)).  The spatial rate *depends on*
  :math:`N` — the defining gating fact — and the cross-term
  :math:`|M|/\max = 0.411` sits three orders above the Cartesian
  ceiling.

* **Cylinder: the extreme of gating.**  There is no pre-floor
  O(h\ :sup:`2`) window at any practical azimuthal order — the
  :math:`(\eta,\varphi)` variation a 1-D :math:`\eta`-march cannot
  thread exactly (#229) — so :math:`E \approx E_{\rm angle}(n_\varphi)`
  and the spatial h-ratio is :math:`\approx 1` at fixed :math:`n_\varphi`.
  The positive signature is the floor's azimuthal scaling: it drops
  :math:`2.58\times` when :math:`n_\varphi` doubles.  (The cylinder's
  small :math:`|M|/\max` is *not* evidence of separability — it is small
  because the angular floor so dominates that the spatial delta
  :math:`\Delta E_h` in the denominator is itself near zero; the gating
  is read from the *saturation* and the *azimuthal scaling*, not the
  cross-term magnitude.)

.. note::

   The scalar (weight-summed) L2 of the table is, by construction, blind
   to a *wrong angular closure* — the Morel--Montry :math:`\alpha`-dome
   telescopes under :math:`\sum_n w_n \psi_n` (vv-principles L27 / the
   per-ordinate-flat-flux discipline).  Because the curvilinear gating is
   *itself* an angular-closure phenomenon, the ST5 gate adds a
   **per-ordinate** leg
   (``test_curvilinear_gating_per_ordinate_not_blind``) that reproduces
   the sphere gating signature from the max-over-ordinates per-ordinate
   L2 (N=8 finest h-ratio 1.16 saturates; N=32 recovers ≈3) — so the
   gate cannot be telescoped blind to a future angular-closure
   regression.  That leg corrects a measured 1/W normalisation trap:
   ``case.psi_exact(r, μ_n)`` returns :math:`A(r) + B(r)\mu_n` *without*
   the :math:`1/W` factor by its own contract, while the solver stores
   the per-ordinate flux *with* it — the reference must be divided by
   :math:`W = \sum_n w_n` before comparison, else a 2× mismatch swamps
   the metric.

The #233 pole-cell × #229 azimuthal-floor interference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The gating law has a concrete consequence for the two open curvilinear
defects.  The spatial pole-cell :math:`\mathcal{O}(h)` order (#233,
documented at :ref:`sn-pole-cell-spatial-closure` and minted as ERR-059
in :ref:`sn-curvilinear-aniso-norm-reconciliation`) and the azimuthal
angular floor (#229) are **not independent contributors** to the
curvilinear error — they *interfere through the gating*.

The mechanism: the angular thread's interpolation error sets the floor
that spatial refinement saturates at.  So the pole-cell spatial defect
(#233) is only *visible* — only the dominant error — once the angular
floor (#229) has been pushed below it by refining :math:`N`.  At a
coarse quadrature the #229 angular floor *masks* the #233 pole-cell
order entirely (the spatial ladder saturates before the
:math:`\mathcal{O}(h)` pole-cell term emerges); only at a fine
quadrature does the spatial ladder run long enough for the pole-cell
order to surface.  This is precisely why the sphere N=8 ladder saturates
at 1.15 while N=32 recovers O(h\ :sup:`2`): the same spatial closure,
read through two different angular floors.

This interference is the reason the two issues must be characterised
*together* rather than as separate spatial and angular bugs, and the
reason a fix to one cannot be validated in isolation: lifting the #229
angular floor (a 2-D angular closure) would *expose* the #233 pole-cell
order that the floor currently masks.  The gating law
:eq:`sn-space-angle-separability` makes this dependency explicit — the
curvilinear error is :math:`\max(E_{\rm space}, E_{\rm angle})`, so
whichever defect is larger *hides* the other.

The permanent pin: the ST5 characterization gate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The decomposition is pinned permanently by
:mod:`tests.sn.verification.mms.test_space_angle_separability` (Issue
#236 Phase 3, sub-task ST5), an **L1** MMS characterisation gate modelled
on the pole-cell characterisation gate
``test_curvilinear_pole_cell_characterization.py`` (#233).  It carries
``@pytest.mark.verifies("sn-space-angle-separability")`` against
:eq:`sn-space-angle-separability`.  Its six legs pin both regimes as
*positive* signatures (never as an xfail-pending-fix):

* **Cartesian, separable** — ``test_cartesian_slab_iso_space_angle_separable``
  (N-independent O(h\ :sup:`2`) spatial rate, :math:`|M|/\max < 0.05`)
  and ``test_cartesian_slab_p1_aniso_floor_n_independent`` (the active-
  angular-axis control: the P1 floor is N-independent and the cross-term
  stays :math:`\approx 0`).
* **Curvilinear, gating** —
  ``test_sphere_spatial_rate_is_quadrature_gated`` (the discriminator:
  coarse-N saturates, fine-N recovers O(h\ :sup:`2`); the proven
  ``@catches("ERR-026")`` catcher via the fine-N O(h\ :sup:`2`)-recovery
  assertion), ``test_sphere_cross_term_large_discriminates_from_cartesian``
  (:math:`|M|/\max > 0.15`),
  ``test_cylinder_spatial_saturates_at_azimuthal_floor`` (spatial
  saturation + azimuthal floor scaling), and
  ``test_curvilinear_gating_per_ordinate_not_blind`` (the L27
  angular-aware per-ordinate leg, also ``@catches("ERR-026")``).

The gate is *characterisation*, not calcification: if a future 2-D
angular closure (#229) or higher-order spatial scheme (#158 / #6) lifts
the curvilinear gating, the gating assertions are designed to redden so
the regime change is *signalled* and the gate is re-tuned to the new
(better) regime — they are not xfails awaiting a fix.  The ``@slow``
mark reflects that the curvilinear solves dominate the ~2 s wall-clock,
not that the gate is optional.

Substituting the WDD Closure into the Balance Equation
=======================================================

Combining the balance equation :eq:`balance-general` with the WDD
angular closure :eq:`wdd-closure` and the standard spatial DD
(:math:`\psi_{n,i} = \frac{1}{2}(\psi_{\rm in}^s + \psi_{\rm out}^s)`,
:math:`\psi_{\rm out}^s = 2\psi_{n,i} - \psi_{\rm in}^s`), define:

.. math::

   c_{\rm out} &= \frac{\alpha_{n+\frac12}}{\tau_n} \\[6pt]
   c_{\rm in}  &= \frac{(1-\tau_n)}{\tau_n}\,\alpha_{n+\frac12}
                 + \alpha_{n-\frac12}

The cell-average angular flux is then:

.. math::
   :label: dd-solve

   \psi_{n,i} = \frac{
       S_i V_i
       + |\mu_n|(A_{\rm in} + A_{\rm out})\,\psi_{\rm in}^s
       + \dfrac{\Delta A_i}{w_n}\, c_{\rm in}\, \psi_{n-\frac12}
   }{
       2|\mu_n|\, A_{\rm out}^s
       + \dfrac{\Delta A_i}{w_n}\, c_{\rm out}
       + \Sigt{} V_i
   }

where the superscript :math:`s` denotes spatial face fluxes, and
:math:`A_{\rm in}`, :math:`A_{\rm out}` are the cell face areas in the
direction of neutron travel (see :ref:`sweep-algorithm` for their
definition).  This is the equation the per-cell
update solves —
:func:`~orpheus.transport.spatial.cell_balance.cell_balance_terms`,
consumed by
:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.update` —
and, in vectorized scan form, the CumprodScan fast path
(:ref:`sn-tau-c-on-cellvisit-live`).

Geometry Comparison
====================

.. list-table::
   :header-rows: 1
   :widths: 15 28 28 29

   * - Aspect
     - Cartesian
     - Spherical
     - Cylindrical
   * - Streaming cosine
     - :math:`\mu`
     - :math:`\mu`
     - :math:`\eta` (radial)
   * - Face area :math:`A`
     - 1 (per unit area)
     - :math:`4\pi r^2`
     - :math:`2\pi r`
   * - Volume :math:`V`
     - :math:`\Delta x`
     - :math:`\tfrac{4}{3}\pi(r_{\rm out}^3 - r_{\rm in}^3)`
     - :math:`\pi(r_{\rm out}^2 - r_{\rm in}^2)`
   * - :math:`\Delta A`
     - 0
     - :math:`4\pi(r_{\rm out}^2 - r_{\rm in}^2)`
     - :math:`2\pi(r_{\rm out} - r_{\rm in})`
   * - Redistribution
     - None
     - :math:`+(\Delta A/w)\,[\alpha\psi]`
     - :math:`+(\Delta A/w)\,[\alpha\psi]`
   * - :math:`\alpha` scope
     - N/A
     - Global (all :math:`N` ordinates)
     - Per :math:`\mu`-level
   * - :math:`\alpha` recursion variable
     - N/A
     - :math:`\mu`
     - :math:`\eta`
   * - Quadrature required
     - GL or Lebedev
     - GL
     - Product or Level-Sym
