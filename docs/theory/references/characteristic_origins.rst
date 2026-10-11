.. _theory-characteristic-origins:

==============================================================================
The characteristic reference's algebra of record — the laws it rests on
==============================================================================

.. contents:: Contents
   :local:
   :depth: 2


.. Machine header — the ``nexus-meta`` schema for this page (PROVISIONAL).

.. dropdown:: Machine header — ``nexus-meta`` schema (PROVISIONAL)
   :color: muted

   .. code-block:: yaml

      module: derivations
      concept: algebra of record of the characteristic reference, specular line closure, rank-1 resolvent, rank-2 resolvent, closed-body identity, white-wall equivalence, vacuum reduction, impact parameter, impact-parameter partition, obliquity, multi-region segments, piecewise optical depth, method of images, retired trajectory-resolvent family
      role: "the Branch-1 algebra of record of the characteristic reference: the SymPy identities of orpheus.derivations.continuous.characteristic.origins (18 derive functions in 6 modules) and the line-geometry, closure and eigenvalue laws the characteristic reference's tests verify numerically, 37 labelled equations, each with the function or the test that verifies it and the place on the characteristic reference page that uses it; and the history of the trajectory-resolvent family (Variant alpha) that the characteristic reference replaced"
      code: [orpheus.derivations.continuous.characteristic.origins, orpheus.derivations.continuous.characteristic.origins.specular, orpheus.derivations.continuous.characteristic.closure, orpheus.derivations.continuous.characteristic.transport, orpheus.geometry.chord]
      depends_on: [characteristic, chart_and_chord]
      related: [characteristic, peierls_nystrom, error_catalog]


Key facts
=========

- **What this page is.** The algebra of record of the characteristic
  reference (:ref:`theory-characteristic-reference`): every equation the
  reference rests on that carries a ``@pytest.mark.verifies`` marker
  today, with the SymPy function or the test that verifies it and the
  section of the characteristic reference page that uses it. There are
  37 such equations `[M]` (2026-10-10, an AST pass over every
  ``verifies(...)`` call in ``tests/`` and ``orpheus/``, 1 105 files
  parsed; the census is the table of
  :ref:`characteristic-origins-census`).
- **Two kinds of evidence, stated per section.** 20 of the 37 are
  *proved symbolically*: an identity of SymPy algebra, a function
  ``derive_*`` in
  :mod:`orpheus.derivations.continuous.characteristic.origins.specular`
  returning a ``dict`` whose ``"pass"`` entry a foundation test asserts.
  20 are *verified numerically*: a test of the characteristic reference,
  of the geometric kernel's chord (``tests/gates/geometry/test_chord.py``)
  or of the Peierls Nyström family's production primitives asserts the
  law on computed output against an independently written value (a hand
  table, an mpmath march, a published number, a closed form). Three
  labels carry both (``V-alpha-2``, ``slab-V-alpha-2``,
  ``slab-architecture``), so 17 are proved only, 17 verified only.
- **The SymPy modules.** Six modules, 18 functions:
  ``greens_function`` (the sphere: three),
  ``greens_function_cylinder`` (six),
  ``greens_function_slab`` (two),
  ``greens_function_slab_asymmetric`` (three),
  ``greens_function_hollow_sphere`` (two) and
  ``greens_function_annulus`` (two). They import SymPy and nothing from
  the characteristic reference, so a proof here shares no code with the
  reference it grounds (the ``algebra-of-record`` skill, "Structural
  independence applies above the trusted-library line"). They are not
  rendered in the API reference; their docstrings are read in the source.
- **One definition per law (X4).** Where a law here is a special case of
  a law the geometric kernel or the characteristic reference defines in
  general, the general statement is the definition and the equation here
  is its specialisation to one body, kept because a test marker names
  it. The line through a concentric body is defined by
  :eq:`geometry-line-crossing-law`, its segment lengths by
  :eq:`geometry-chord-segment-lengths` and its obliquity by
  :eq:`geometry-cylinder-axial-factor` (:ref:`theory-chart-and-chord`);
  the period of a line by :eq:`characteristic-transit-rank`, its closure
  by :eq:`characteristic-closure` and the integrals along it by
  :eq:`characteristic-traversal-integrals`
  (:ref:`theory-characteristic-reference`). Each section below names the
  general law its equations specialise.
- **The labels are addresses.** Every label here begins
  ``peierls-greens-``, the name of the angle-resolved Green's-function
  campaign that wrote them for the trajectory-resolvent family. The
  prefix names no object of the characteristic reference; the labels
  keep it because test markers name them, and a label is renamed only
  together with its markers.
- **The family these equations came from is retired.** The
  trajectory-resolvent family, which solved the same problems by power
  iteration on an angle-resolved collocation grid, was deleted on
  2026-10-10 (P1 step (e2) of the characteristic reference campaign,
  ``5aa8cb88``, #405). What it was, what failed in it and why the
  characteristic reference replaced it is the history at the end of this
  page (:ref:`characteristic-origins-history`).


.. _characteristic-origins-notation:

Notation: the family's symbols and the reference's
==================================================

The equations below keep the symbols they were derived in. Each names an
object of the characteristic reference:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Symbol here
     - The characteristic reference's object
   * - :math:`\alpha`, :math:`\alpha_L`, :math:`\alpha_R`,
       :math:`\alpha_{\rm in}`, :math:`\alpha_{\rm out}`
     - the specular amplitude :math:`a_k` of the wall a traversal exits at
       (:ref:`characteristic-walls`): 1 for a mirror, 0 for vacuum, the
       albedo of a partial specular wall.
   * - :math:`b`
     - the impact parameter, the least orbit coordinate along a line
       (:eq:`geometry-line-crossing-law`).
   * - :math:`\mu_{\rm surf}` (sphere)
     - :math:`\sqrt{1 - (b/R)^2}`, the cosine to the normal at which the
       line meets the outer wall.
   * - :math:`s_{\rm in\!-\!plane} = \sqrt{1 - \mu_{\rm axial}^2}`
       (cylinder)
     - the projected speed :math:`|P\Omega| = \sin\theta`, whose reciprocal
       is the obliquity (:eq:`geometry-cylinder-axial-factor`).
   * - :math:`L_0`, :math:`L_{\rm first}`
     - the backward distance from a point to the wall its line entered
       at, along :math:`-\Omega`.
   * - :math:`L_p`, :math:`L_{\rm period}`
     - the length of one traversal of a solid body's line, wall to wall.
   * - :math:`\tau_{\rm step}`, :math:`\tau_{\rm period}`
     - the optical depth :math:`\tau_k` of a traversal
       (:eq:`characteristic-traversal-integrals`).
   * - :math:`F(r, \mu)`
     - the source integral along the first leg, from the entry wall to the
       point: the angular flux of the vacuum Volterra block
       (:ref:`characteristic-volterra`).
   * - :math:`B`, :math:`B_{LR}`, :math:`B_{RL}`
     - the outflow :math:`B_k` of a traversal, its source integral
       attenuated to its exit (:eq:`characteristic-traversal-integrals`).
   * - :math:`\psi_{\rm surf}`, :math:`\psi_L^+`, :math:`\psi_R^-`,
       :math:`\psi_{\rm in}^{\rm out}`, :math:`\psi_{\rm out}^{\rm in}`
     - the inflow :math:`\psi^{\rm in}_k` to a traversal
       (:eq:`characteristic-closure`).
   * - :math:`T(\mu_{\rm surf})`, :math:`T = (I - S)^{-1}`
     - the resolvent of the line's cycle: :math:`1/(1 - \Pi)` with the
       cycle product :math:`\Pi` (:ref:`characteristic-closure-section`).
   * - :math:`q`, :math:`q_g`
     - the isotropic emission density per steradian, the unknown of the
       multigroup system (:eq:`characteristic-pencil`).

A *rank-1* line is a line whose period is one traversal: every line of a
solid sphere or cylinder, and a shell's line that misses the cavity. A
*rank-2* line has two: a slab's line between two returning walls and a
shell's line through the cavity (:ref:`characteristic-period`).


.. _characteristic-origins-sphere:

The solid sphere: the closed body, the white wall and vacuum
============================================================

*Kind: proved symbolically* (``greens_function``, three functions; the
foundation tests of
``tests/gates/derivations/test_peierls_greens_function_symbolic.py``).
*General laws specialised:* :eq:`characteristic-closure` at rank 1 and
the angular flux on a line (:ref:`characteristic-angular-flux`).

The angular flux on a line of a solid sphere
--------------------------------------------

Take a point at radius :math:`r` and a direction of cosine :math:`\mu`
to the outward radius. Its line has the impact parameter
:math:`b = r\sqrt{1 - \mu^2}` (a special case of
:eq:`geometry-line-crossing-law`, where :math:`b` is the distance from
the centre to the line), meets the outer wall :math:`R` at the cosine
:math:`\mu_{\rm surf} = \sqrt{R^2 - b^2}/R`, and its traversal is the
chord :math:`L_p = 2R\mu_{\rm surf} = 2\sqrt{R^2 - b^2}`. The backward
distance from the point to the wall the line entered at is

.. math::

   L_0(r, \mu) = r\mu + \sqrt{R^2 - r^2(1 - \mu^2)},

the larger root of :math:`|x - s\,\Omega|^2 = R^2` in the backward
arc length :math:`s`. Write :math:`F` for the source integral over the
first leg and :math:`B` for the outflow of one traversal,

.. math::

   F(r, \mu) = \int_0^{L_0} q\bigl(|x - s\,\Omega|\bigr)\,
       e^{-\Sigma_t s}\,\mathrm d s, \qquad
   B(\mu_{\rm surf}) = \int_0^{L_p} q\bigl(c(s)\bigr)\,
       e^{-\Sigma_t (L_p - s)}\,\mathrm d s,

the special cases of :eq:`characteristic-traversal-integrals` on a
homogeneous sphere. The intensity entering the line at the wall is
returned by a mirror (amplitude 1) after one traversal, attenuated by
:math:`e^{-\Sigma_t L_p}`, plus that traversal's outflow:

.. math::
   :label: peierls-greens-surface-fixed-point

   \psi_{\rm surf} = B(\mu_{\rm surf}) +
       e^{-\Sigt{}\,L_p}\,\psi_{\rm surf}.

This is :eq:`characteristic-closure` with :math:`m = 1` and
:math:`a_0 = 1`. Its least solution is the geometric series over the
returns, :math:`\psi_{\rm surf} = B/(1 - e^{-\Sigma_t L_p})`; with a
specular amplitude :math:`\alpha` the cycle product is
:math:`\Pi = \alpha e^{-\Sigma_t L_p}` and the inflow is
:math:`\alpha B/(1 - \alpha e^{-\Sigma_t L_p})`, the rank-1 resolvent
:math:`(1 - \alpha e^{-\tau})^{-1}` that Sanchez 1986 Appendix Eq. (A4)
:cite:`SanchezTTSP1986` and Pomraning–Siewert 1982 Eq. (14)
:cite:`PomraningSiewert1982` write by direct integration of the transfer
equation. The angular flux at the point is the first leg plus the inflow
carried down it:

.. math::
   :label: peierls-greens-function-architecture

   \boxed{\;
   \psi(r_i, \mu_q) = F(r_i, \mu_q) +
       e^{-\Sigt{}\,L_0(r_i,\mu_q)}\,
       T(\mu_{\rm surf})\,B(\mu_{\rm surf})
   \;}

with :math:`T(\mu_{\rm surf}) = 1/(1 - e^{-\Sigma_t L_p})`. On the
characteristic reference this is :meth:`TraversalRule.angular_flux
<orpheus.derivations.continuous.characteristic.transport.TraversalRule.angular_flux>`
on a rank-1 line (:ref:`characteristic-angular-flux`).

V_α1: a closed homogeneous sphere has the infinite medium's eigenvalue
-----------------------------------------------------------------------

.. math::
   :label: peierls-greens-V-alpha-1

   (K \cdot 1)(r, \mu) = \omega_0,
   \qquad \omega_0 = \frac{\Sigs{}}{\Sigt{}},

where :math:`K` is the scattering operator of the homogeneous sphere
behind a mirror: the angular flux produced by the emission
:math:`q = \Sigma_s\,\psi_{\rm trial}`. The constant is an eigenfunction
of :math:`K` with eigenvalue :math:`\omega_0`: the transport of a
constant emission :math:`q` on the closed sphere is the constant
:math:`q/\Sigma_t`. A constant flux is therefore the fundamental mode of
the fission problem exactly when its collisions balance its emissions,
:math:`\Sigma_t\phi = \Sigma_s\phi + \nu\Sigma_f\phi/k`, which gives
:math:`k = \nu\Sigma_f/(\Sigma_t - \Sigma_s) = \nu\Sigma_f/\Sigma_a`,
the infinite medium's. The proof is
:func:`~orpheus.derivations.continuous.characteristic.origins.specular.greens_function.derive_operator_constant_trial_closed_sphere`,
in three steps.

*Step 1, the surface fixed point.* With :math:`q = \Sigma_s` constant,
:math:`B = (\Sigma_s/\Sigma_t)(1 - e^{-\Sigma_t L_p})`, and
:eq:`peierls-greens-surface-fixed-point` becomes
:math:`(1 - e^{-\Sigma_t L_p})\,\psi_{\rm surf} =
(\Sigma_s/\Sigma_t)(1 - e^{-\Sigma_t L_p})`, so
:math:`\psi_{\rm surf} = \omega_0` whatever :math:`L_p`: the chord
cancels. (``test_v_alpha1_surface_fixed_point_solves_to_q_over_sigma_t``.)

*Step 2, the first leg cancels.* With the same source,
:math:`F = (\Sigma_s/\Sigma_t)(1 - e^{-\Sigma_t L_0})`, and
:eq:`peierls-greens-function-architecture` gives

.. math::

   \psi(r, \mu) = \frac{\Sigs{}}{\Sigt{}}\,(1 - e^{-\Sigt{}\,L_0})
       + e^{-\Sigt{}\,L_0}\,\frac{\Sigs{}}{\Sigt{}} = \omega_0 ,

independent of :math:`L_0`, hence of :math:`(r, \mu)`.
(``test_v_alpha1_total_psi_is_independent_of_first_leg``.)

*Step 3, the eigenvalue.* The kernel of :math:`K` is positive, so the
constant, a positive eigenfunction, belongs to the spectral radius
(Perron–Frobenius): :math:`\omega_0` is the dominant eigenvalue of
:math:`K` and :math:`k = \nu\Sigma_f/\Sigma_a`.
(``test_v_alpha1_operator_on_constant_gives_omega_0``,
``test_v_alpha1_overall_pass``.)

**Where the reference uses it.** The closed homogeneous body is the
:math:`k_\infty` floor of the characteristic reference's multigroup
system: a closed body of one mixture reads :math:`k = \rho(A^{-1}F)` with
a flat flux in the infinite medium's group ratio
(:eq:`peierls-greens-cylinder-mr-kinf` below;
:ref:`characteristic-galerkin-system`). V_α1 is the one-group statement
proved on the continuous closure, before any discretisation.

**What V_α1 cannot see.** The identity holds for any :math:`L_0` and any
:math:`L_p`, so a wrong first-leg length or a wrong chord passes it: a
constant emission evaluated anywhere on a line is the same constant. The
family's first prototype integrated the first leg over the forward
distance :math:`\sqrt{R^2 - r^2(1-\mu^2)} - r\mu` instead of the
backward one, and V_α1 stayed green; a vacuum sphere against
Pomraning–Siewert's reference found it, 6 % in :math:`k`
(:ref:`characteristic-origins-history-failures`). A closed-form identity
on a constant source verifies the closure's algebra, never the geometry
it is fed; the geometry is the kernel's, verified against closed forms
(:ref:`theory-chart-and-chord`).

V_α2: at rank 1 the specular closure is the white wall's
---------------------------------------------------------

The rank-1 term of the specular transmission matrix of the Peierls
Nyström family (:ref:`theory-peierls-nystrom`), in the isotropic mode
:math:`\tilde P_0 = 1`, is

.. math::
   :label: peierls-greens-T00-integrand

   T_{00}^{\rm sphere} = 2 \int_0^1 \mu\,\tilde P_0(\mu)^2\,
       e^{-2\Sigt{} R \mu}\,\mathrm d\mu
   = 2 \int_0^1 \mu\,e^{-2\Sigt{} R \mu}\,\mathrm d\mu,

and the self-collision probability of the sphere's surface, Hébert's
:math:`P_{ss}` (Hébert 2009 §3.8.5 :cite:`Hebert2020`), is the polar
integral :math:`2\int_0^{\pi/2}\sin\theta'\cos\theta'
e^{-2\tau_R\cos\theta'}\,\mathrm d\theta'`. SymPy integrates each in its
own variable and both reduce to one closed form:

.. math::
   :label: peierls-greens-V-alpha-2

   T_{00}^{\rm sphere} = P_{ss}^{\rm sphere} =
       \frac{1 - (1 + 2\tau_R)\,e^{-2\tau_R}}{2\,\tau_R^{2}},
   \qquad \tau_R = \Sigt{}\,R.

The proof is
:func:`~orpheus.derivations.continuous.characteristic.origins.specular.greens_function.derive_T00_equals_P_ss_sphere`:
path A integrates the transfer-matrix definition in :math:`\mu`, path B
the escape integral in :math:`\theta'` with its :math:`\sin\theta'`
Jacobian, and ``simplify`` of their difference is 0. The two paths are
different integrals reaching one value, which is what makes the
equality evidence: an earlier version built both sides from one SymPy
literal, and its ``simplify(LHS - RHS) == 0`` was a tautology
(:ref:`characteristic-origins-history-failures`). The production
primitives of the Peierls Nyström family reproduce it numerically,
``compute_T_specular_sphere`` at rank 1 against ``compute_P_ss_sphere``
over :math:`\tau_R \in \{0.5, 1, 2.5, 5, 10\}`
(``tests/gates/derivations/test_peierls_specular_primitives.py::test_v_alpha2_sphere_T00_equals_Pss_via_production_primitives``,
which verifies :eq:`peierls-greens-V-alpha-2` beside the four SymPy rows).

**What it means.** At rank 1 the specular closure and the white wall's
:math:`(1 - P_{ss})^{-1}` act on the isotropic mode identically: a mirror
and a white wall are one closure on a closed homogeneous sphere's
fundamental mode, and differ on every other.

**Where the reference uses it.** :math:`P_{ss}` is the wall-to-wall
transmission of a white sphere: the characteristic reference's white
walls' gates assert its line part's transmission against :math:`P_{ss}`
on the sphere and the cylinder and against :math:`2E_3(\tau)` on the slab
(:eq:`peierls-greens-slab-V-alpha-2`), the escape probability two ways
beside it (:ref:`characteristic-wall-coupling`, "The white walls'
gates"); and the gate
``test_the_white_and_specular_laws_are_their_closed_forms_and_differ``
asserts that at :math:`\alpha = 0.5` a white and a specular sphere
differ, which is where V_α2's rank-1 coincidence ends.

V_α3: an absorbing wall returns nothing
---------------------------------------

Sanchez 1986 Eq. (A6) :cite:`SanchezTTSP1986` writes the boundary part
of the angle-integrated kernel of the homogeneous sphere,

.. math::

   g_h(\rho'\to\rho) = 2\alpha \int_{\mu_0}^{1} T(\mu_-)\,\mu_*^{-1}\,
       \cosh(\rho\mu)\,\cosh(\rho'\mu_*)\,e^{-2 a \mu_-}\,\mathrm d\mu ,

with :math:`a` the optical radius and :math:`T(\mu_-)` the rank-1
resolvent. The leading :math:`2\alpha` makes the whole integrand
proportional to the amplitude:

.. math::
   :label: peierls-greens-V-alpha-3

   g_h(\rho' \to \rho)\bigr|_{\alpha = 0} = 0,

so the full kernel :math:`\bar g_2 + g_h` reduces to the vacuum kernel
:math:`\bar g_2` with no branch. The proof is
:func:`~orpheus.derivations.continuous.characteristic.origins.specular.greens_function.derive_alpha_zero_kernel_reduction`.
On the characteristic reference vacuum is a wall of amplitude 0, not a
missing wall: the period runs through it and its inflow is exactly 0,
bitwise (``test_every_amplitude_zero_adds_nothing``,
:ref:`characteristic-closure-section`). The angle-integrated kernel
:math:`g_h` itself is never assembled there: its diagonal is
hypersingular (:ref:`characteristic-origins-history-family`).


.. _characteristic-origins-cylinder:

The solid cylinder: the obliquity and the in-plane chord
========================================================

*Kind:* the in-plane speed, the impact parameter and the period chord are
*proved symbolically* (``greens_function_cylinder``,
``derive_bounce_period_chord_cylinder``); the first-leg chord, the
closure, the angular flux and the V_α2 identity are *verified numerically* on the
characteristic reference and the geometric kernel. *General laws
specialised:* :eq:`geometry-cylinder-axial-factor`,
:eq:`geometry-line-crossing-law`, :eq:`characteristic-closure`.

The obliquity, the impact parameter and the period
--------------------------------------------------

A cylinder's line has an axial cosine :math:`\mu_{\rm axial}` and, in the
plane normal to the axis, an azimuth :math:`\varphi_{\rm az}` measured
from the outward radius at the point. Its image in that plane is traversed
at the speed

.. math::
   :label: peierls-greens-cylinder-in-plane-speed

   s_{\rm in\!-\!plane}(\mu_{\rm axial}) =
       \sin\theta_{\rm axis} = \sqrt{1 - \mu_{\rm axial}^2},

the cylinder's row of :eq:`geometry-cylinder-axial-factor`: an in-plane
length :math:`\ell` is the 3-D length :math:`\ell/\sin\theta`. The
kernel stores the speed and sums it from the direction's components,
never as :math:`1 - \mu_{\rm axial}^2`, which loses digits near the axis
(:ref:`chart-and-chord-obliquity`). A specular reflection on the
cylinder's wall keeps :math:`\mu_{\rm axial}` and the in-plane impact
parameter

.. math::
   :label: peierls-greens-cylinder-impact-parameter

   b(r, \varphi_{\rm az}) = r\,|\sin\varphi_{\rm az}|,

the distance from the axis to the line's image, so every traversal of a
line is the same in-plane chord at :math:`b`, lifted by the same
obliquity:

.. math::
   :label: peierls-greens-cylinder-bounce-period

   L_{\rm period}(b, \mu_{\rm axial}) =
       \frac{2\sqrt{R^2 - b^2}}{\sqrt{1 - \mu_{\rm axial}^2}}.

:func:`~orpheus.derivations.continuous.characteristic.origins.specular.greens_function_cylinder.derive_bounce_period_chord_cylinder`
derives :math:`L_{\rm period}` two ways (from the impact parameter, and
from the surface tangent at the bounce point) and proves the two agree;
its foundation test,
``test_v_alpha1_cyl_bounce_period_chord_two_derivations_agree`` in
``test_peierls_greens_function_cylinder_symbolic.py``, verifies all three
labels. This is the cylinder's row of :eq:`geometry-chord-segment-lengths`
with one region, both sides of the closest approach summed.

The first leg
-------------

The backward in-plane distance from a point to the wall its line entered
at, and its 3-D length, are

.. math::
   :label: peierls-greens-cylinder-trajectory

   L_{\rm 2D, first}(r, \varphi_{\rm az}) =
       r\,\cos\varphi_{\rm az} +
       \sqrt{R^2 - r^2 \sin^2\varphi_{\rm az}}, \qquad
   L_{0, \rm 3D} =
       \frac{L_{\rm 2D, first}}{s_{\rm in\!-\!plane}}.

*Verified numerically* by
``tests/gates/geometry/test_chord.py::test_the_multi_region_segments_are_the_hand_written_table``:
from a start at :math:`r\cos\varphi = 0.8`, :math:`r\sin\varphi = b`, the
kernel's half-line slots beyond the start sum to
:math:`(0.8 + \sqrt{R^2 - b^2})/|P\Omega|`, each slot asserted on its own
against a hand-written closed form, the obliquity on the rows with
:math:`\Omega_z = 0.8`. On the characteristic reference the first leg is
not a separate object: it is the part of the line's transit behind the
point (:ref:`characteristic-angular-flux`).

The closure and the angular flux
--------------------------------

A line of a solid cylinder has rank 1, and its inflow is the rank-1 case
of :eq:`characteristic-closure` with the cylinder's traversal:

.. math::
   :label: peierls-greens-cylinder-T

   \psi_{\rm surf}(b, \mu_{\rm axial}) =
       \frac{\alpha\,B(b, \mu_{\rm axial})}
            {1 - \alpha\,e^{-\Sigma_t L_{\rm period}}},

with :math:`B` the traversal's outflow. The leading :math:`\alpha` is the
amplitude of the wall the traversal exits at; the denominator's is the
cycle product of the one-traversal period. The angular flux at a point is

.. math::
   :label: peierls-greens-cylinder-architecture

   \psi(r_i, \mu_{\rm axial}, \varphi_{\rm az}) =
       F(r_i, \mu_{\rm axial}, \varphi_{\rm az})
       + e^{-\Sigma_t L_{0, \rm 3D}}\,\psi_{\rm surf}(b, \mu_{\rm axial}),

and the scalar flux its integral over the directions,
:math:`\phi(r) = \int_{-1}^{1}\mathrm d\mu_{\rm axial}
\int_0^{2\pi}\mathrm d\varphi_{\rm az}\,\psi`, which on the reference is a
point reading over the lines through the point
(:eq:`characteristic-quadrature`).

*Verified numerically.*
:eq:`peierls-greens-cylinder-T` is asserted by
``test_characteristic_transport.py::test_the_closed_angular_flux_is_the_unfolded_backward_path``,
case ``cylinder_solid_partial`` (a partial mirror on an oblique line,
against an mpmath march wall by wall, to :math:`10^{-13}`); its first red
is the obliquity dropped.
:eq:`peierls-greens-cylinder-architecture` is asserted by
``test_characteristic_reading.py::test_the_reading_at_a_general_point_of_a_cylinder_is_the_mpmath_route``
(``slow``): off the axis of a three-region cylinder, under a partial
mirror, the reading against a two-dimensional mpmath route over the
in-plane disc, to :math:`10^{-12}`; its vacuum rows null
:math:`\psi_{\rm surf}`, so the partial-mirror row is the witness.

V_α2 and V_α3 on the cylinder
-----------------------------

The cylinder's rank-1 specular term and the self-collision probability
of its surface are one integral:

.. math::
   :label: peierls-greens-cylinder-V-alpha-2

   T_{00}^{\rm cyl} = P_{ss}^{\rm cyl}
   = \frac{4}{\pi}\int_0^{\pi/2}\cos\alpha\,
     \mathrm{Ki}_3\bigl(2\Sigma_t R\cos\alpha\bigr)\,\mathrm d\alpha ,

with :math:`\alpha` here the in-plane angle of the chord to the inward
normal at the wall, not an amplitude.
:func:`~orpheus.derivations.continuous.characteristic.origins.specular.greens_function_cylinder.derive_T00_equals_P_ss_cylinder`
proves it at the level of the integrand. The cylinder has no elementary
closed form: both sides reach
:math:`(4/\pi)\cos\alpha\,\mathrm{Ki}_3(2\Sigma_t R\cos\alpha)` through
the Bickley–Naylor identity
:math:`\mathrm{Ki}_3(x) = \int_0^{\pi/2}\sin^2\beta\,e^{-x/\sin\beta}\,\mathrm d\beta`,
applied as a substitution. Because both sides end at the same special
function, the SymPy proof shares an upstream identity (the ERR-032 risk
of the ``vv-principles`` skill, anti-pattern #7), and the evidence is the
numerical row: ``compute_T_specular_cylinder_3d`` at rank 1 against
``compute_P_ss_cylinder``, two production primitives with different
intermediate quantities, at :math:`\tau_R \in \{0.5, 1, 2.5, 5, 10\}`
(``test_v_alpha2_cyl_T00_equals_Pss_via_production_primitives``, the
label's one verifier; the SymPy rows carry no marker).
:func:`~orpheus.derivations.continuous.characteristic.origins.specular.greens_function_cylinder.derive_alpha_zero_kernel_reduction_cylinder`
proves that the closure's leading :math:`\alpha` removes the wall's
contribution at :math:`\alpha = 0`; it carries no label.


.. _characteristic-origins-multiregion:

Lines through several regions
=============================

*Kind:* the piecewise optical depth, the multi-region closure and the
homogeneous reduction are *proved symbolically* (``greens_function_cylinder``,
three functions); the segment law, the emission's region and the four
eigenvalue laws are *verified numerically*. *General laws specialised:*
:eq:`geometry-chord-segment-lengths`, :eq:`geometry-crossing-order`,
:eq:`characteristic-traversal-integrals`, :eq:`characteristic-pencil`.

The segments
------------

Along a cylinder line's image, at in-plane arc length
:math:`s_{\rm 2D}` from a start at radius :math:`r_{\rm start}`, the
radius is

.. math::
   :label: peierls-greens-cylinder-mr-trajectory-segments

   r(s_{\rm 2D})^{2} = r_{\rm start}^{2}
                       - 2\,r_{\rm start}\cos\varphi_{\rm az}\,s_{\rm 2D}
                       + s_{\rm 2D}^{2},
   \qquad s_{\rm 2D} \in [0, L_{\rm 2D, first}],

the conic of :eq:`geometry-line-crossing-law` measured from the start
instead of the closest approach. It meets an interior circle
:math:`R_k` at :math:`s_{\rm 2D} = r_{\rm start}\cos\varphi_{\rm az} \pm
\sqrt{R_k^2 - r_{\rm start}^2\sin^2\varphi_{\rm az}}`, real iff
:math:`R_k \ge b`. The kernel never forms a segment length as a
difference of these roots: it writes
:math:`(r_{k+1}^2 - r_k^2)/(h_{k+1} + h_k)`, with no cancellation however
close the line passes to a surface, and it reads each segment's region
from the crossing that opens it, never by locating a point
(:eq:`geometry-chord-segment-lengths`, :eq:`geometry-crossing-order`).

*Verified numerically* by
``test_chord.py::test_the_multi_region_segments_are_the_hand_written_table``:
each slot end of the kernel's chord through a three-region body is where
this conic meets a breakpoint, on the half-lines from two starting radii,
against a reference written in the test as differences of float
half-chords, a spelling the kernel does not use.

The emission reads its own region
---------------------------------

The isotropic emission density of a piecewise-homogeneous body jumps at
every interface, because the cross sections do, while the scalar flux is
continuous. Each segment of a line therefore reads the piece of the
emission that belongs to its own region:

.. math::
   :label: peierls-greens-mr-regionwise-source

   q(r(s)) \;=\; Q_k\bigl(r(s)\bigr),
   \qquad s \in (s_a, s_b) \text{ of a segment in region } k,

with :math:`Q_k` the emission's representation on region :math:`k`
alone. On the characteristic reference :math:`Q_k` is a polynomial per
panel, and the panel ends refine the body's breakpoints, so no panel
straddles an interface: ``PanelBasis`` refuses a partition whose panel
ends do not contain the breakpoints (:ref:`characteristic-panel-basis`).

*Verified numerically* by
``test_characteristic_transport.py::test_the_outflow_of_a_per_region_polynomial_is_its_line_integral_attenuated_to_the_exit``:
a line's outflow against an mpmath line integral of an emission that is a
different cubic in each region, jumping at both interfaces, with
:math:`\Sigma_t` jumping too, to :math:`10^{-13}` relative, on lines of
both ranks through every chart. It also catches ERR-090, the family's one
cubic spline across the interfaces (:ref:`characteristic-origins-history-failures`).

The piecewise optical depth
---------------------------

The obliquity does not depend on the segment, so it factors out of the
optical depth of a cylinder traversal:

.. math::
   :label: peierls-greens-cylinder-mr-piecewise-tau

   \tau_{\rm period}(b, \mu_{\rm axial}) =
       \frac{1}{\sqrt{1 - \mu_{\rm axial}^{2}}}\,
       \sum_k \Sigma_{t,k}\,\Delta s_{\rm 2D,\,k},

with :math:`\Delta s_{\rm 2D,k}` the in-plane length of segment :math:`k`.
:func:`~orpheus.derivations.continuous.characteristic.origins.specular.greens_function_cylinder.derive_piecewise_3d_optical_depth_cylinder_mr`
proves the factoring (the sum of :math:`\Sigma_{t,k}\Delta s_{\rm 3D,k}`
with each segment lifted equals the lifted sum), and
:func:`~orpheus.derivations.continuous.characteristic.origins.specular.greens_function_cylinder.derive_homogeneous_limit_reducibility_cylinder_mr`
proves that with every :math:`\Sigma_{t,k}` equal the sum collapses to
:math:`\Sigma_t L_{\rm period}`. Both carry this label. The factoring is
what lets the kernel solve the in-plane chord once and rescale it for
every axial cosine (:ref:`chart-and-chord-obliquity`). Inside a segment,
the attenuation at a point carries the optical depth of every earlier
segment: :math:`\tau(s) = \tau_{\rm back} + \Sigma_{t,k}(s - s_a)/\sin\theta`;
restarting it at each region drops the earlier segments' attenuation
(the battery arm ``attenuation-restarted-at-region`` reddens the
equal-material interface rows).

The multi-region closure and its homogeneous limit
--------------------------------------------------

The rank-1 closure is unchanged when the traversal crosses several
regions; only its two inputs become sums over the segments:

.. math::
   :label: peierls-greens-cylinder-mr-bounce-sum-piecewise

   \psi_{\rm surf}(b, \mu_{\rm axial}) =
       \frac{\alpha\,B}{1 - \alpha\,e^{-\tau_{\rm period}}}, \qquad
   B = \sum_k \int_{s_{\rm 2D,\,a}^{(k)}}^{s_{\rm 2D,\,b}^{(k)}}
         \frac{q(r(s_{\rm 2D}))\,e^{-\tau(s_{\rm 2D})}}{
               \sqrt{1 - \mu_{\rm axial}^{2}}}\,
         \mathrm d s_{\rm 2D}.

With every region of one material the piecewise sums are the homogeneous
ones:

.. math::
   :label: peierls-greens-cylinder-mr-homogeneous-reduction

   \forall k:\;\;\Sigma_{t,k} = \Sigma_t \;\;\Longrightarrow\;\;
   \begin{cases}
       \tau_{\rm period}^{\rm MR}
           = \displaystyle\sum_k \Sigma_t\,\Delta s_{\rm 3D,\,k}
           = \Sigma_t\,L_{\rm period}, \\[1ex]
       B^{\rm MR}
           = \displaystyle\sum_k \int_{s_a^{(k)}}^{s_b^{(k)}}
                 q(r(s))\,e^{-\Sigma_t s_{\rm 3D}}\,
                 \mathrm d s_{\rm 3D}
           = \displaystyle\int_0^{L_{\rm period}}
                 q(r(s))\,e^{-\Sigma_t s_{\rm 3D}}\,
                 \mathrm d s_{\rm 3D}.
   \end{cases}

:func:`~orpheus.derivations.continuous.characteristic.origins.specular.greens_function_cylinder.derive_two_region_constant_source_consistency_cylinder_mr`
proves both on two regions with a constant source: the exact piecewise
outflow

.. math::

   B = \frac{q}{\Sigma_{t,1}}\!\left(1 - e^{-\Sigma_{t,1}\ell_1}\right)
     + \frac{q}{\Sigma_{t,2}}\,e^{-\Sigma_{t,1}\ell_1}\!
       \left(1 - e^{-\Sigma_{t,2}\ell_2}\right)

reduces under :math:`\Sigma_{t,1} = \Sigma_{t,2}`, :math:`\alpha = 1` to
:math:`\psi_{\rm surf} = q/\Sigma_t`, V_α1's value. Its test,
``test_v_alpha1_cyl_mr_two_region_constant_source_homogeneous_limit``,
verifies both labels. On the characteristic reference the same reduction
is a gate on the assembled block: an interface between two equal
materials is invisible
(``test_characteristic_assembly.py::test_an_interface_between_equal_materials_is_invisible_on_a_cylinder``,
``slow``, :math:`3.1 \times 10^{-15}`).

The four eigenvalue laws of a multi-region body
-----------------------------------------------

*Verified numerically*, each by a test of the characteristic reference.

**A closed body reads the infinite medium.** A closed body of one mixture
loses no neutron, so its eigenvalue is the infinite medium's, whatever its
shape:

.. math::
   :label: peierls-greens-cylinder-mr-kinf

   {\bf A} = \mathrm{diag}(\Sigma_t) - \Sigma_s^{\rm T},
   \qquad
   {\bf F} = \chi\otimes\nu\Sigma_f,

with :math:`k_\infty = \rho({\bf A}^{-1}{\bf F})` and :math:`\Sigma_s`
stored :math:`[g_{\rm from}, g_{\rm to}]`, so the in-scatter source is
:math:`\Sigma_s^{\rm T}\phi`.
``test_characteristic_system.py::test_a_closed_body_reads_k_inf_with_a_flat_flux_in_the_infinite_mediums_group_ratio``
asserts :math:`k = \rho({\bf A}^{-1}{\bf F})` to :math:`10^{-12}` and a
flux flat to :math:`10^{-10}` in the infinite medium's group ratio, on
spheres, slabs and shells (the cylinder rows ``slow``). Its two-group
mixtures have asymmetric scattering, so the scattering matrix transposed
moves :math:`k` by order 1; the one-group rows cannot see that (the
``vv-principles`` skill, one-group degeneracy).

**A published critical cylinder.** The one-group bare cylinder Ua-1-0-CY
of Sood, Forster and Parsons 2003 (Table 10, :math:`c = 1.30`) is
critical at the printed radius:

.. math::
   :label: peierls-greens-cylinder-mr-wm72-vacuum

   k\bigl(K=1,\,\alpha=0,\,c=1.30,\,r = r_c^{\rm Sood}\bigr) = 1
   \quad \text{within} \quad
   \Bigl|\frac{\mathrm dk}{\mathrm dr}\Bigr|\,\frac{\delta r_{\rm print}}{2}
   + \epsilon_{\rm ladder},
   \qquad
   r_c^{\rm WM\text{-}72} = r_c^{\rm Sood}\ \text{to}\ 2 \times 10^{-6},

with :math:`\delta r_{\rm print}` a unit of the printed radius's last
digit and :math:`\epsilon_{\rm ladder}` the reference's error estimated
from its own ladder.

``test_characteristic_independent_references.py::test_k_is_one_at_the_one_group_cylinders_published_critical_radius``
(``slow``) asserts that band at the printed radius 1.72500292 mfp
(no fixed :math:`10^{-5}`: the family's row applied one), and that the Westfall–Metcalf singular-eigenfunction critical
radius (:ref:`theory-singular-eigenfunction`) equals Sood's to
:math:`2 \times 10^{-6}`: two truths of one quantity, neither through a
Bickley–Naylor function. `[M]` 2026-10-10 (the step (e1b) record,
``probes/probe_ua_cy.log``): :math:`k - 1 = -8.4 \times 10^{-6}` and
:math:`-2.7 \times 10^{-7}` at rungs 2 and 3. Its fast-tier twin,
``test_k_is_one_at_the_one_group_cylinders_published_critical_radius_in_the_fast_tier``,
reads rung 2 to :math:`3 \times 10^{-5}`.

**The scalar flux is continuous at an interface.** The angular flux has
a kink at an interface (the emission jumps with the cross sections), but
its integral over directions does not:

.. math::
   :label: peierls-greens-cylinder-mr-interface-continuity

   \phi(r) = \int_{-1}^{1}\!\mathrm d\mu_{\rm axial}\,
       \int_0^{2\pi}\!\mathrm d\varphi_{\rm az}\;
       \psi(r, \mu_{\rm axial}, \varphi_{\rm az})

is continuous at every interior radius :math:`r_k`.
``test_characteristic_albedos.py::test_the_modes_scalar_flux_is_continuous_across_a_material_interface``
reads the fundamental mode at :math:`r_k(1 \mp \epsilon)`: the jump
falls by more than 30 times as :math:`\epsilon` falls from
:math:`10^{-3}` to :math:`10^{-5}` and again to :math:`10^{-7}`, and is
below :math:`10^{-5}` relative at :math:`10^{-7}` (`[M]` 2026-10-10, the
sphere at rung 3: :math:`6.2 \times 10^{-4}`, :math:`1.1 \times 10^{-5}`,
:math:`1.6 \times 10^{-7}` at :math:`r = 1`). The point reading
transports the emission to each side, so the two readings share nothing
but the emission; a point value read from the Galerkin coefficients
instead, discontinuous at panel ends, reddens it. It is blind to any
error in the transport itself, which keeps the reading continuous.

**The reference contracts on every resolution axis.**

.. math::
   :label: peierls-greens-cylinder-mr-quadrature-convergence

   \forall\,\text{axis} \in \{\text{degree},\,\text{layers},\,
       \text{line},\,\text{along}\}:
   \qquad
   |\Delta k^{(i+1)}| < |\Delta k^{(i)}|
   \quad \text{until } |\Delta k| < 10^{-12}\,|k|,

with each axis varied alone from the door's working point: the panels'
polynomial degree, the grading layers toward walls and interfaces, the
points per piece of the line measure (where grazing lines live) and the
points per arc-length piece along a line.
``test_characteristic_convergence.py::test_k_contracts_on_every_resolution_axis_and_the_working_point_is_within_the_old_floor``
asserts it on two-group heterogeneous bodies at intermediate albedos,
each angular step at most :math:`10^{-2}` of the one before it, and the
working point's error estimated from the steps below
:math:`3 \times 10^{-4}`; ``test_k_contracts_in_the_panel_degree_on_a_two_region_cylinder``
(``slow``) asserts the degree axis on a two-region cylinder. A ladder
estimates an error and never certifies one; the row's teeth are the
contraction (:ref:`characteristic-evidence`).


.. _characteristic-origins-slab:

The slab
========

*Kind:* the V_α2 identity and the vacuum reduction are *proved
symbolically* (``greens_function_slab``, two functions); the first leg,
the closure and the V_α2 identity again are *verified numerically*, on
the kernel's chord, the characteristic reference and a production
primitive. *General laws specialised:*
:eq:`geometry-chord-segment-lengths` (the slab's row),
:eq:`characteristic-closure` at rank 2.

The first leg
-------------

A slab's line at cosine :math:`\mu` to the normal crosses the slab
:math:`[0, L]` monotonically. The backward distance from :math:`x` to the
wall behind it is

.. math::
   :label: peierls-greens-slab-trajectory

   L_{\rm first}(x, \mu) = \begin{cases}
       x / \mu       & \mu > 0 \text{ (came from } x = 0 \text{)} \\
       (L-x) / |\mu| & \mu < 0 \text{ (came from } x = L \text{)}
   \end{cases}.

The position at backward arc length :math:`s` is :math:`x - \mu s`: the
arc length times the cosine, never the arc length alone (ERR-034,
:ref:`characteristic-origins-history-failures`). *Verified numerically*
by ``test_chord.py::test_the_slab_chord_in_both_orientations``: the
crossing parameter of the wall behind the point on each orientation of an
asymmetric slab, :math:`-r_0/0.6` and :math:`r_3/0.6`.

The closure on a symmetric slab
-------------------------------

A slab's line between two returning walls has rank 2: the transit and its
reverse. At equal amplitudes, **and** with equal outflows of the two
transits, :math:`B_{LR} = B_{RL} = B` (a source symmetric about the
mid-plane :math:`x = L/2`), the cycle's least solution is

.. math::
   :label: peierls-greens-slab-T

   \psi_{\rm surf}(\mu) = \frac{\alpha\,B(\mu)}
                                {1 - \alpha\,e^{-\Sigma_t L/|\mu|}},
   \qquad \alpha_L = \alpha_R = \alpha,\quad B_{LR} = B_{RL} = B,

with :math:`B(\mu)` the outflow of **one** transit, the source integral
over :math:`L/|\mu|`, not over the out-and-back path. It is the
symmetric reduction of the rank-2 resolvent
:eq:`peierls-greens-slab-asym-resolvent` below: the determinant
:math:`1 - \alpha^2 e^{-2\tau}` factors as
:math:`(1 - \alpha e^{-\tau})(1 + \alpha e^{-\tau})`, and the second
factor cancels against the numerator's off-diagonal entries only when
the two outflows are equal. Without that hypothesis the inflows differ,
:math:`\psi_L^+ - \alpha B_{LR}/(1 - \alpha e^{-\tau}) =
\alpha^2 e^{-\tau}(B_{LR} - B_{RL})/(1 - \alpha^2 e^{-2\tau})`, and the
symmetric slab's inflow is the general rank-2 cycle
:eq:`peierls-greens-slab-asym-closure`. The angular-flux row's
``slab_symmetric_half`` case runs a source that is not symmetric, so it
verifies that general cycle (``characteristic-closure``), of which this
equation is the symmetric reduction. The
heuristic closure :math:`\alpha B_{\rm period}/(1 - \alpha^2 e^{-2\tau})`,
built by analogy with the rank-1 sphere, agrees with this form only at
:math:`\alpha \in \{0, 1\}` (ERR-035,
:ref:`characteristic-origins-history-failures`). The angular flux at a
point is

.. math::
   :label: peierls-greens-slab-architecture

   \psi(x_i, \mu) =
       F(x_i, \mu)
       + e^{-\Sigma_t L_{\rm first}}\,\psi_{\rm surf}(\mu),

with :math:`\phi(x) = 2\pi\int_{-1}^{1}\psi(x, \mu)\,\mathrm d\mu`.

:func:`~orpheus.derivations.continuous.characteristic.origins.specular.greens_function_slab.derive_alpha_zero_kernel_reduction_slab`
proves that the wall's term vanishes at :math:`\alpha = 0`
(``test_v_alpha3_slab_psi_surf_vanishes_at_alpha_zero`` verifies
:eq:`peierls-greens-slab-architecture`), on the closure
:eq:`peierls-greens-slab-T` with the one-transit optical depth. The
limit is decided by the leading :math:`\alpha` alone, so it holds for
ERR-035's heuristic closure too and cannot tell the two apart; the
closure itself is asserted below. The row withdrawn under #506,
``test_characteristic_nystrom_withdrawn.py::test_the_vacuum_slabs_k_is_the_nystrom_slabs``,
also names it: the vacuum slab :math:`[0, 10]` reads the Peierls Nyström
slab's :math:`k` (an E\ :sub:`1` kernel with no lines) to
:math:`5 \times 10^{-5}`. The closure :eq:`peierls-greens-slab-T` is
*verified numerically*, under its hypotheses, by
``test_characteristic_transport.py::test_a_symmetric_slabs_surface_inflow_is_the_one_transit_closed_form``:
a homogeneous slab with one albedo (0.3, 0.5, 0.9) on both faces and the
source :math:`q(x) = 1 + 0.8(x - L/2)^2`, symmetric about the mid-plane;
the inflow at each traversal's entry against the closed form written in
mpmath along the transit (40 digits), on a slab of width 2.3 with
:math:`\Sigma_t = 0.7` at :math:`\mu = 0.6`, both traversals' inflows to
:math:`10^{-13}` relative. It reddens when the source is made
non-symmetric (3 of 3 rows: the hypothesis is load-bearing) and under
ERR-035's heuristic (3 of 3; `[M]` 2026-10-10, at :math:`\alpha = 0.5`
the heuristic reads 1.0143 against 0.9819). The general cycle for any source, ERR-035's
regime included, is the angular-flux row's ``slab_symmetric_half`` case
(``verifies("characteristic-closure")``).

V_α2 on the slab
----------------

.. math::
   :label: peierls-greens-slab-V-alpha-2

   T_{00}^{\rm slab} = P_{\rm ss}^{\rm slab} = 2\,E_3(\Sigma_t L)

:func:`~orpheus.derivations.continuous.characteristic.origins.specular.greens_function_slab.derive_T00_equals_P_ss_slab`
proves it in a hybrid of closed form and arbitrary precision (states 1A
and 1B of the ``algebra-of-record`` skill). SymPy fails to integrate
:math:`\int_0^1 2\mu\,e^{-\tau/\mu}\,\mathrm d\mu` at the endpoint
:math:`\mu \to 0^+` (``Add object cannot be interpreted as an integer``,
the skill's choke mode 1), so the substitution :math:`u = 1/\mu` is
verified symbolically to map the integrand into the definition
:math:`E_3(x) = \int_1^\infty t^{-3}e^{-xt}\,\mathrm dt`, and mpmath
integrates both original integrands (in :math:`\mu` and in
:math:`\theta`) against ``mpmath.expint(3, τ)`` at six optical
thicknesses from 0.1 to 10, to :math:`10^{-12}`. The production primitive
row, ``test_v_alpha2_slab_T00_equals_2E3_via_production_primitive``,
compares ``compute_T_specular_slab`` at rank 1 with ``mpmath.expint``;
it and the four rows of ``test_peierls_greens_function_slab_symbolic.py``
(``test_v_alpha2_slab_substitution_algebra_holds`` among them) verify
this label.
:math:`2E_3(\tau)` is the face-to-face transmission of a white slab, which
the characteristic reference's white walls' gates assert
(:ref:`characteristic-wall-coupling`).


.. _characteristic-origins-rank-two:

Two returning walls: the rank-2 cycle
=====================================

*Kind:* the resolvent and the closures on constant sources are *proved
symbolically* (``greens_function_slab_asymmetric``,
``greens_function_hollow_sphere``, ``greens_function_annulus``); the
method of images, the impact-parameter partition and the through-cavity
resolvents are *verified numerically* on the characteristic reference.
*General law specialised:* :eq:`characteristic-closure` at :math:`m = 2`,
and the rank of :eq:`characteristic-transit-rank`.

The resolvent of a slab between two walls
-----------------------------------------

The two inflows of a rank-2 line, :math:`\psi_L^+` entering at the left
wall and :math:`\psi_R^-` at the right, are coupled by one transit of
optical depth :math:`\tau = \Sigma_t L/|\mu|` and the amplitude of the
wall at its end. The step from one inflow to the next is the
anti-diagonal matrix
:math:`S = \bigl(\begin{smallmatrix}0 & \alpha_L e^{-\tau}\\
\alpha_R e^{-\tau} & 0\end{smallmatrix}\bigr)` (a transit never returns
to the wall it left), and the sum over returns is

.. math::
   :label: peierls-greens-slab-asym-resolvent

   T(\alpha_L, \alpha_R, \tau) = (I - S)^{-1}
       = \frac{1}{1 - \alpha_L\,\alpha_R\,e^{-2\tau}}
         \begin{pmatrix}
             1                       & \alpha_L\,e^{-\tau} \\
             \alpha_R\,e^{-\tau}     & 1
         \end{pmatrix}.

:func:`~orpheus.derivations.continuous.characteristic.origins.specular.greens_function_slab_asymmetric.derive_rank2_resolvent_slab_asymmetric`
inverts :math:`I - S` with ``Matrix.inv()``, checks the closed form, and
reduces it to the symmetric pair, to one vacuum wall and to two
(``test_v_alpha2_slab_asym_determinant_canonical_form`` and its
siblings). The determinant is the cycle product of
:eq:`characteristic-closure`, :math:`1 - \gamma_0\gamma_1`; the second
form of that equation is this matrix applied to the two returns
:math:`r_k = a_k B_k`. The inflows satisfy

.. math::
   :label: peierls-greens-slab-asym-closure

   \psi_L^+(\mu) &= \alpha_L\,B_{LR}(\mu) + \alpha_L\,e^{-\tau}\,
                       \psi_R^-(\mu), \\
   \psi_R^-(\mu) &= \alpha_R\,B_{RL}(\mu) + \alpha_R\,e^{-\tau}\,
                       \psi_L^+(\mu),

with the one-transit outflows
:math:`B_{LR} = \int_0^{L/|\mu|} q(|\mu|s)\,e^{-\Sigma_t s}\,\mathrm ds`
and :math:`B_{RL} = \int_0^{L/|\mu|} q(L - |\mu|s)\,e^{-\Sigma_t s}\,\mathrm ds`.
:func:`~orpheus.derivations.continuous.characteristic.origins.specular.greens_function_slab_asymmetric.derive_operator_constant_trial_closed_slab_asymmetric`
proves that on a constant source with both walls mirrors both inflows are
:math:`q/\Sigma_t` (the factor :math:`(1 - e^{-\tau})(1 + e^{-\tau})/(1 -
e^{-2\tau}) = 1` collapses them). The angular flux at a point is

.. math::
   :label: peierls-greens-slab-asym-architecture

   \psi(x_i, \mu) = F(x_i, \mu)
                   + e^{-\Sigma_t L_{\rm first}}\,
                     \psi_{\rm surface}(\mu),

with :math:`\psi_{\rm surface} = \psi_L^+` for :math:`\mu > 0` and
:math:`\psi_R^-` for :math:`\mu < 0`;
:func:`~orpheus.derivations.continuous.characteristic.origins.specular.greens_function_slab_asymmetric.derive_alpha_zero_kernel_reduction_slab_asymmetric`
proves that at :math:`\alpha_L = \alpha_R = 0`, :math:`S = 0`,
:math:`T = I`, and the flux is the first leg alone.

**Which wall's amplitude multiplies which inflow.** The inflow to a
traversal carries the amplitude of the wall the *previous* traversal
exits at. Swapping the pairing is invisible wherever the two walls have
one amplitude; the characteristic reference's gates draw distinct
amplitudes (0.3 and 0.6) for that reason (:ref:`characteristic-closure-section`,
"The albedo pairing").

The method of images
--------------------

.. math::
   :label: peierls-greens-slab-asym-method-of-images

   \boxed{\;
   k_{\rm eff}\!\big[\text{slab } [0, L],\, \alpha_L = 1,\, \alpha_R = 0\big]
   \;=\;
   k_{\rm eff}\!\big[\text{slab } [0, 2L],\, \alpha_L = \alpha_R = 0\big]
   \;}

The fundamental mode of the symmetric vacuum slab :math:`[0, 2L]` is even
about :math:`x = L`, so its centre plane carries exactly a mirror's
condition: the half slab with a mirror at the cut and vacuum at the other
wall has the same eigenvalue, and its flux is the right half of the
doubled slab's. *Verified numerically* by two rows of
``test_characteristic_albedos.py``:
``test_a_mirror_is_the_symmetry_plane_of_the_doubled_vacuum_slab``
(:math:`k` to :math:`10^{-8}`) and
``test_the_mirrored_slabs_flux_is_twice_the_doubled_slabs_right_half``
(:math:`\phi_{\rm half}(x) = 2\,\phi_{\rm doubled}(1 + x)` at four points,
every group, to :math:`10^{-4}`, both gauged to a production of 1 over
their body). Both carry ``catches("ERR-034")``.

**What the identity cannot see.** It compares two posings of one
reference, so an error the two posings share passes it. The retired
family satisfied it to :math:`10^{-7}` while both its sides were
:math:`5.3 \times 10^{-4}` off the value the characteristic reference and
the Peierls Nyström slab converge to
(:ref:`characteristic-origins-history-failures`;
:ref:`characteristic-successors`). The identity catches what a defect does
to the relation between the two posings, as ERR-034 does; it is not a
value against an independent reference.

The hollow sphere and the annulus: the impact-parameter partition
-----------------------------------------------------------------

A line through a shell :math:`[R_{\rm in}, R_{\rm out}]` meets the cavity
iff its impact parameter is below the inner radius. On the hollow sphere
and on the annulus that parameter is

.. math::
   :label: peierls-greens-hollow-sph-impact-parameter-partition

   b(r, \mu) = r\sqrt{1 - \mu^2}

.. math::
   :label: peierls-greens-annulus-impact-parameter-partition

   b(r, \varphi_{\rm az}) = r\,|\sin\varphi_{\rm az}|

and the partition at :math:`b = R_{\rm in}` decides the rank: a line with
:math:`b \ge R_{\rm in}` misses the cavity, has one transit from the
outer wall to the outer wall and rank 1, and is a solid body's line at the
same impact parameter; a line with :math:`b < R_{\rm in}` has two
transits, outer wall to inner and inner to outer, and rank 2. At
:math:`b = R_{\rm in}` exactly the line touches the inner surface, and a
tangency is not a crossing (:ref:`chart-and-chord-tangency`): rank 1.
The characteristic reference derives the rank from the period of the
line and never tests :math:`b` against :math:`R_{\rm in}`
(:eq:`characteristic-transit-rank`; :ref:`characteristic-period`, "The
ranks, derived"). *Verified numerically* by
``test_characteristic_closure.py::test_the_period_matches_the_hand_counted_table``,
whose hollow rows sit below, at, one ulp below and above the inner radius,
the walls counted by hand.

The shell's traversal through the cavity has the optical depth
:math:`\tau_{\rm step}(b) = \Sigma_t\,(\sqrt{R_{\rm out}^2 - b^2} -
\sqrt{R_{\rm in}^2 - b^2})` on the sphere; on the annulus the in-plane
length is the same and the obliquity lifts it:

.. math::
   :label: peierls-greens-annulus-3d-chord-scaling

   \tau_{\rm step}^{\rm annulus}(b, \mu_{\rm axial})
        = \frac{\tau_{\rm step}^{\rm hollow\,sph}(b)}
               {\sqrt{1 - \mu_{\rm axial}^2}}
        = \Sigma_t \cdot
          \frac{\sqrt{R_{\rm out}^2 - b^2} - \sqrt{R_{\rm in}^2 - b^2}}
               {\sqrt{1 - \mu_{\rm axial}^2}}.

:func:`~orpheus.derivations.continuous.characteristic.origins.specular.greens_function_annulus.derive_3d_chord_scaling_annulus`
proves it and that the outer-only period scales the same way
(``test_v_alpha2_annulus_aux_3d_chord_scaling_step``). It is the shell's
row of :eq:`geometry-chord-segment-lengths` (the cavity slot) divided by
the obliquity, and it fails if any code path integrates the axial
direction out before the closure.

The through-cavity rank-2 resolvent of the two shells is the slab's with
:math:`\tau_{\rm step}` in place of :math:`\Sigma_t L/|\mu|`:

.. math::
   :label: peierls-greens-hollow-sph-through-rank2

   T(\alpha_{\rm in}, \alpha_{\rm out}, \tau_{\rm step})
       = \frac{1}{1 - \alpha_{\rm in}\,\alpha_{\rm out}\,
                     e^{-2\tau_{\rm step}}}
         \begin{pmatrix}
             1                                 & \alpha_{\rm in}\,
                                                  e^{-\tau_{\rm step}} \\
             \alpha_{\rm out}\,e^{-\tau_{\rm step}}    & 1
         \end{pmatrix}

.. math::
   :label: peierls-greens-annulus-through-rank2

   T(\alpha_{\rm in}, \alpha_{\rm out}, \tau_{\rm step})
       = \frac{1}{1 - \alpha_{\rm in}\,\alpha_{\rm out}\,
                     e^{-2\tau_{\rm step}}}
         \begin{pmatrix}
             1                                   & \alpha_{\rm in}\,
                                                    e^{-\tau_{\rm step}} \\
             \alpha_{\rm out}\,e^{-\tau_{\rm step}}      & 1
         \end{pmatrix},

the second with :math:`\tau_{\rm step}` from
:eq:`peierls-greens-annulus-3d-chord-scaling`. *Verified numerically* by
``test_characteristic_transport.py::test_the_closed_angular_flux_is_the_unfolded_backward_path``,
cases ``sphere_hollow_cavity`` and ``cylinder_hollow_cavity``: the angular
flux on a line through the cavity, with amplitudes 0.3 inside and 0.6
outside, against an mpmath march wall by wall to :math:`10^{-13}`; the
2 × 2 matrix is the least solution of that line's cycle written in closed
form.

At :math:`b = R_{\rm in}` the inner chord :math:`\sqrt{R_{\rm in}^2 - b^2}`
vanishes and the through-cavity closure meets the outer-only one: the two
branches are continuous across the partition. The angular flux at a point
of either shell is

.. math::
   :label: peierls-greens-hollow-sph-architecture

   \psi(r, \mu) = F(r, \mu) + e^{-\Sigma_t\,L_{\rm first}(r, \mu)}
                                \cdot \psi_{\rm surface}(r, \mu),

.. math::
   :label: peierls-greens-annulus-architecture

   \psi(r, \mu_{\rm axial}, \varphi_{\rm az}) =
       F(r, \mu_{\rm axial}, \varphi_{\rm az}) +
       e^{-\Sigma_t\,L_{\rm first}^{\rm 3D}(r, \mu_{\rm axial},
                                                  \varphi_{\rm az})}
       \cdot \psi_{\rm surface}(r, \mu_{\rm axial}, \varphi_{\rm az}),

with :math:`\psi_{\rm surface}` the rank-1 inflow for :math:`b \ge R_{\rm in}`
and, for :math:`b < R_{\rm in}`, the inflow at the wall the point's
transit entered at (the inner wall for an outward point, the outer wall
for an inward one).
:func:`~orpheus.derivations.continuous.characteristic.origins.specular.greens_function_hollow_sphere.derive_operator_constant_trial_closed_hollow_sphere`
and
:func:`~orpheus.derivations.continuous.characteristic.origins.specular.greens_function_annulus.derive_operator_constant_trial_closed_annulus`
prove that both branches give :math:`\psi = q/\Sigma_t` on a constant
source with both walls mirrors, by two cancellations: the rank-1
:math:`(1 - x)/(1 - x)` and the rank-2
:math:`(1 - e^{-\tau})(1 + e^{-\tau})/(1 - e^{-2\tau})`, and that a wall
that leaks gives strictly less.
:func:`~orpheus.derivations.continuous.characteristic.origins.specular.greens_function_hollow_sphere.derive_alpha_zero_kernel_reduction_hollow_sphere`
proves that at :math:`\alpha_{\rm in} = \alpha_{\rm out} = 0` both
branches reduce to the first leg.

**A vacuum cavity is second order in its radius on a sphere, first order
on a cylinder.** A vacuum cavity absorbs the lines through it, the
fraction :math:`(R_{\rm in}/R)^2` of the lines on the sphere. `[M]`
2026-10-10 (the step (e1b) record, ``probes/probe_rin.py``): the distance
of :math:`k` to the solid sphere's is :math:`1.9 \times 10^{-2}`,
:math:`2.2 \times 10^{-4}`, :math:`1.8 \times 10^{-6}` at
:math:`R_{\rm in}/R = 10^{-1}, 10^{-2}, 10^{-3}`, and the cylinder's
:math:`9.4 \times 10^{-2}`, :math:`1.1 \times 10^{-2}`,
:math:`1.2 \times 10^{-3}`.
``test_characteristic_albedos.py::test_a_vacuum_shells_k_tends_to_the_solid_bodys_as_its_cavity_closes``
asserts the order, two-sided, and :math:`k_{\rm hollow} < k_{\rm solid}`.


.. _characteristic-origins-census:

The census: every live label and its verifier
=============================================

`[M]` 2026-10-10: an AST pass over the decorators of every function in
``tests/`` (1 105 Python files with ``orpheus/``), collecting the string
arguments of ``verifies(...)``; positive control
``peierls-greens-slab-T``, found. 36 of the 67 labels the
trajectory-resolvent page carried are named by a marker, and one more,
``peierls-greens-cylinder-V-alpha-2``, was minted here for a row whose
marker named the closure while it asserts the :math:`T_{00}` identity;
this page carries those 37 and no other ``:label:``. Kind P is proved symbolically, N
verified numerically. Test paths are under ``tests/gates/derivations/``
unless they begin ``geometry/``.

.. list-table::
   :header-rows: 1
   :widths: 34 6 60

   * - Label (``peierls-greens-`` plus)
     - Kind
     - Verifier
   * - ``surface-fixed-point``
     - P
     - ``test_peierls_greens_function_symbolic.py::test_v_alpha1_surface_fixed_point_solves_to_q_over_sigma_t``
   * - ``function-architecture``
     - P
     - ``test_peierls_greens_function_symbolic.py::test_v_alpha1_total_psi_is_independent_of_first_leg``
   * - ``V-alpha-1``
     - P
     - ``test_peierls_greens_function_symbolic.py``: the two rows above,
       ``test_v_alpha1_operator_on_constant_gives_omega_0``,
       ``test_v_alpha1_overall_pass``
   * - ``T00-integrand``
     - P
     - ``test_peierls_greens_function_symbolic.py::test_v_alpha2_T00_matches_hebert_via_matrix_path``
   * - ``V-alpha-2``
     - P, N
     - ``test_peierls_greens_function_symbolic.py``: the row above,
       ``test_v_alpha2_Pss_matches_hebert_via_polar_path``,
       ``test_v_alpha2_T00_equals_Pss_closed_forms``,
       ``test_v_alpha2_overall_pass``;
       ``test_peierls_specular_primitives.py::test_v_alpha2_sphere_T00_equals_Pss_via_production_primitives``
   * - ``V-alpha-3``
     - P
     - ``test_peierls_greens_function_symbolic.py::test_v_alpha3_g_h_vanishes_at_alpha_zero``
   * - ``cylinder-in-plane-speed``, ``cylinder-impact-parameter``,
       ``cylinder-bounce-period``
     - P
     - ``test_peierls_greens_function_cylinder_symbolic.py::test_v_alpha1_cyl_bounce_period_chord_two_derivations_agree``
   * - ``cylinder-trajectory``
     - N
     - ``geometry/test_chord.py::test_the_multi_region_segments_are_the_hand_written_table``
   * - ``cylinder-T``
     - N
     - ``test_characteristic_transport.py::test_the_closed_angular_flux_is_the_unfolded_backward_path``
   * - ``cylinder-V-alpha-2``
     - N
     - ``test_peierls_specular_primitives.py::test_v_alpha2_cyl_T00_equals_Pss_via_production_primitives``
   * - ``cylinder-architecture``
     - N
     - ``test_characteristic_reading.py::test_the_reading_at_a_general_point_of_a_cylinder_is_the_mpmath_route``
       (``slow``)
   * - ``cylinder-mr-trajectory-segments``
     - N
     - ``geometry/test_chord.py::test_the_multi_region_segments_are_the_hand_written_table``
   * - ``mr-regionwise-source``
     - N
     - ``test_characteristic_transport.py::test_the_outflow_of_a_per_region_polynomial_is_its_line_integral_attenuated_to_the_exit``
   * - ``cylinder-mr-piecewise-tau``
     - P
     - ``test_peierls_greens_function_cylinder_symbolic.py::test_v_alpha1_cyl_mr_homogeneous_reducibility``,
       ``::test_v_alpha1_cyl_mr_piecewise_3d_optical_depth``
   * - ``cylinder-mr-bounce-sum-piecewise``,
       ``cylinder-mr-homogeneous-reduction``
     - P
     - ``test_peierls_greens_function_cylinder_symbolic.py::test_v_alpha1_cyl_mr_two_region_constant_source_homogeneous_limit``
   * - ``cylinder-mr-kinf``
     - N
     - ``test_characteristic_system.py::test_a_closed_body_reads_k_inf_with_a_flat_flux_in_the_infinite_mediums_group_ratio``
   * - ``cylinder-mr-wm72-vacuum``
     - N
     - ``test_characteristic_independent_references.py::test_k_is_one_at_the_one_group_cylinders_published_critical_radius``
       (``slow``) and ``::..._in_the_fast_tier``
   * - ``cylinder-mr-interface-continuity``
     - N
     - ``test_characteristic_albedos.py::test_the_modes_scalar_flux_is_continuous_across_a_material_interface``
   * - ``cylinder-mr-quadrature-convergence``
     - N
     - ``test_characteristic_convergence.py::test_k_contracts_on_every_resolution_axis_and_the_working_point_is_within_the_old_floor``,
       ``::test_k_contracts_in_the_panel_degree_on_a_two_region_cylinder``
       (``slow``)
   * - ``slab-trajectory``
     - N
     - ``geometry/test_chord.py::test_the_slab_chord_in_both_orientations``
   * - ``slab-T``
     - N
     - ``test_characteristic_transport.py::test_a_symmetric_slabs_surface_inflow_is_the_one_transit_closed_form``
   * - ``slab-architecture``
     - P, N
     - ``test_peierls_greens_function_slab_symbolic.py::test_v_alpha3_slab_psi_surf_vanishes_at_alpha_zero``;
       ``test_characteristic_nystrom_withdrawn.py::test_the_vacuum_slabs_k_is_the_nystrom_slabs``
       (withdrawn under #506)
   * - ``slab-V-alpha-2``
     - P, N
     - ``test_peierls_greens_function_slab_symbolic.py``:
       ``test_v_alpha2_slab_substitution_algebra_holds``,
       ``test_v_alpha2_slab_closed_form_equals_canonical``,
       ``test_v_alpha2_slab_numerical_path_independence``,
       ``test_v_alpha2_slab_overall_pass``;
       ``test_peierls_specular_primitives.py::test_v_alpha2_slab_T00_equals_2E3_via_production_primitive``
   * - ``slab-asym-resolvent``
     - P
     - ``test_peierls_greens_function_slab_asymmetric_symbolic.py::test_v_alpha2_slab_asym_determinant_canonical_form``
   * - ``slab-asym-closure``
     - P
     - ``test_peierls_greens_function_slab_asymmetric_symbolic.py::test_v_alpha1_slab_asym_psi_L_plus_at_closed_BC_equals_q_over_sigma_t``
   * - ``slab-asym-architecture``
     - P
     - ``test_peierls_greens_function_slab_asymmetric_symbolic.py::test_v_alpha3_slab_asym_psi_surf_vanishes_at_vacuum_vacuum``
   * - ``slab-asym-method-of-images``
     - N
     - ``test_characteristic_albedos.py::test_a_mirror_is_the_symmetry_plane_of_the_doubled_vacuum_slab``,
       ``::test_the_mirrored_slabs_flux_is_twice_the_doubled_slabs_right_half``
   * - ``hollow-sph-impact-parameter-partition``,
       ``annulus-impact-parameter-partition``
     - N
     - ``test_characteristic_closure.py::test_the_period_matches_the_hand_counted_table``
   * - ``hollow-sph-through-rank2``, ``annulus-through-rank2``
     - N
     - ``test_characteristic_transport.py::test_the_closed_angular_flux_is_the_unfolded_backward_path``
   * - ``hollow-sph-architecture``
     - P
     - ``test_peierls_greens_function_hollow_sphere_symbolic.py``:
       ``test_v_alpha1_hollow_sph_outer_only_closed_BC_equals_q_over_sigma_t``,
       ``test_v_alpha3_hollow_sph_psi_surf_vanishes_at_vacuum_vacuum``
   * - ``annulus-3d-chord-scaling``
     - P
     - ``test_peierls_greens_function_annulus_symbolic.py::test_v_alpha2_annulus_aux_3d_chord_scaling_step``
   * - ``annulus-architecture``
     - P
     - ``test_peierls_greens_function_annulus_symbolic.py::test_v_alpha1_annulus_outer_only_closed_BC_equals_q_over_sigma_t``

The 18 ``derive_*`` functions by module:

- ``greens_function`` (sphere): ``derive_operator_constant_trial_closed_sphere``,
  ``derive_T00_equals_P_ss_sphere``, ``derive_alpha_zero_kernel_reduction``;
- ``greens_function_cylinder``: ``derive_bounce_period_chord_cylinder``,
  ``derive_T00_equals_P_ss_cylinder``,
  ``derive_alpha_zero_kernel_reduction_cylinder``,
  ``derive_homogeneous_limit_reducibility_cylinder_mr``,
  ``derive_piecewise_3d_optical_depth_cylinder_mr``,
  ``derive_two_region_constant_source_consistency_cylinder_mr``;
- ``greens_function_slab``: ``derive_T00_equals_P_ss_slab``,
  ``derive_alpha_zero_kernel_reduction_slab``;
- ``greens_function_slab_asymmetric``:
  ``derive_operator_constant_trial_closed_slab_asymmetric``,
  ``derive_rank2_resolvent_slab_asymmetric``,
  ``derive_alpha_zero_kernel_reduction_slab_asymmetric``;
- ``greens_function_hollow_sphere``:
  ``derive_operator_constant_trial_closed_hollow_sphere``,
  ``derive_alpha_zero_kernel_reduction_hollow_sphere``;
- ``greens_function_annulus``:
  ``derive_operator_constant_trial_closed_annulus``,
  ``derive_3d_chord_scaling_annulus``.

Five more were deleted in P1 step (e2) because they proved the same
identity as one of these on another geometry (the closed cylinder's and
the closed slab's constant-source identity, the hollow sphere's and the
annulus's rank-2 inversion, the annulus's vacuum reduction).


.. _characteristic-origins-gotchas:

Gotchas
=======

- **A SymPy identity on a constant source cannot see the geometry it is
  fed.** On a constant emission every chord length cancels from the
  closed-body identities, so a wrong first-leg length, a wrong period, a
  position advanced by the wrong factor all pass V_α1 and its siblings.
  The geometry has its own gates on the kernel, against closed forms; a
  defect of the geometry inside a closure is caught by a non-uniform
  source at an intermediate amplitude, never by a closed body.
- **Two derivations ending at one special function share it.** The
  cylinder's V_α2 reaches :math:`\mathrm{Ki}_3` on both sides through
  the same Bickley–Naylor identity, so its SymPy proof is not
  structurally independent; the production-primitive row is the
  evidence.
- **A label's name is not its verifier's claim.** ``peierls-greens-slab-T``
  and ``peierls-greens-cylinder-T`` state the specular closure, and until
  2026-10-10 the V_α2 rows of the slab and the cylinder also named them
  while asserting the :math:`T_{00}` identity; those rows now verify
  ``slab-V-alpha-2`` and ``cylinder-V-alpha-2``. The cylinder's closure
  is verified by the angular-flux row; the slab's has a row of its own,
  because :eq:`peierls-greens-slab-T` holds only for a source symmetric
  about the mid-plane while the angular-flux row's slab case is not. Read the test body, not the
  label, before crediting a marker (the ``vv-principles`` skill, "Log
  every caught bug").
- **An identity between two posings of one reference shares their
  error.** The method of images, the closing cavity's order and the
  equal-material interface compare a reference with itself; they catch a
  defect that breaks the relation and are blind to one that moves both
  sides.
- **The origins modules are not in the API reference.** A
  ``:func:`` role to a ``derive_*`` function renders as plain text; the
  docstrings are read in the source. Their docstrings still describe the
  retired family's solvers in places; the code is the SymPy body, not the
  prose around it.


.. _characteristic-origins-history:

History
=======

This section records the family these equations were derived for, what
it computed, what failed in it and why the characteristic reference
replaced it. The migration of its tests is
:ref:`characteristic-successors`; each catalogued defect's full entry is
in :doc:`/theory/verification/error_catalog`.

.. _characteristic-origins-history-family:

The trajectory-resolvent family (Variant α), 2026-04 to 2026-10-10
------------------------------------------------------------------

**Why it was built.** The Peierls Nyström family's specular sphere
(Phase 4, ``boundary="specular_multibounce"``) closes the boundary with a
rank-:math:`N` projection, and that closure carries a truncation bias that
does not vanish at the ranks production reaches: on the fuel-A-like
closed sphere :math:`(R, \Sigma_t, \Sigma_s, \nu\Sigma_f) = (5, 0.5,
0.38, 0.025)`, :math:`k_\infty = 0.208\overline{3}`, rank 1 (identically
the white wall's :math:`(1 - P_{ss})^{-1}`, which is V_α2) read
0.20777642780, rank 2 0.20781489733, rank 3 0.20808891764, errors of
0.27 %, 0.25 % and 0.12 %. On the two-region Class B sphere of Issue #132
(fuel inner, moderator outer at radii 0.5 and 1.0) its rank-2 closure
read 1.015 against a homogenised :math:`k_\infty` of 0.648, a 57 %
catastrophe (:ref:`peierls-rank-n-class-b-mr-mg-falsification`). Phase 5
(Issue #133, closed 2026-04-28) tried to discretise Sanchez 1986 Eq.
(A6) :cite:`SanchezTTSP1986` instead, the angle-integrated kernel
:math:`g_h(\rho'\to\rho)`, by Nyström sampling, and found that kernel is
not of Fredholm second kind: at the discrete diagonal the angle-resolved
factors carry a :math:`1/(\cos\omega_i\cos\omega_j)` Jacobian whose
poles meet on the visibility cone, a non-integrable
:math:`1/(\mu^2 - \mu_0^2)` Hadamard finite-part singularity
(:ref:`peierls-continuous-mu-retreat`). Four paths were closed for good:
Nyström sampling of :math:`g_h`, Galerkin double integration of it over a
Lagrange basis, a bounce-resolved expansion as a Nyström target (the
singularity persists with no bounce at all, so it is in the first-flight
factors, not in the bounce sum) and Hadamard regularisation (gauge
ambiguous, unbounded). The lesson the family was built on: *do not
integrate over angle before the boundary is closed*. Work with the
angle-resolved kernel of Sanchez Eqs. (A1) and (A5), :math:`t = \bar t +
t_h`, along single characteristics, where every piece is finite, and
integrate over angle only at the end, where the bounce sum's
:math:`1/\mu` is integrable. The angle-integrated kernel and the
angle-resolved one are not different discretisations of one operator:
they are different operators with the same physical content, which is
why Phase 4's basis projection produced sane numbers from a kernel
Phase 5 could not sample. The characteristic reference keeps the lesson:
it closes each line before it integrates over lines.

**What it computed.** The angular flux on a collocation grid of nodes
:math:`(r_i, \mu_q)` (the sphere), :math:`(r_i, \mu_{\rm axial},
\varphi_{\rm az})` (the cylinder) or :math:`(x_i, \mu)` (the slab), by
:eq:`peierls-greens-function-architecture`: a first-leg integral, a
period integral and the closed-form bounce sum, each evaluated by
composite Gauss–Legendre along the chord (``n_traj``, 32 to 128 nodes),
the emission between radial nodes read from a cubic spline through the
nodal values. The eigenvalue came from power iteration with a Rayleigh
quotient, started from a guess :math:`k` passed in (``initial_k``,
assigned as ``k_eff = float(initial_k)`` at seven sites). It covered six
geometries on two orbit-space classes (:ref:`orbit-space-m-g-classification`):
the sphere, the cylinder and the symmetric slab, the asymmetric slab, the
hollow sphere and the annulus, one material each, with multi-group and,
on the sphere and the cylinder, multi-region arms. The code was twelve
numeric modules in
``orpheus/derivations/continuous/trajectory_resolvent/``: fifteen
``solve_greens_function_*`` entry points; a shared closure module,
``variant_alpha_core``, with a rank-1 closure carrying an
``alpha_per_period`` keyword and a rank-2 closure; seven chord-oracle
classes behind a ``ChordOracle`` protocol, one per geometry and region
count; a family ``power_iteration``; a ``Billiard`` facade dispatching on
an eight-case string tag, ``geometry_kind``, with a ``closure_rank``
field; and a reference reading, ``TrajectoryResolventDerivation``
(:ref:`characteristic-origins-history-reading`).

**The billiard frame.** The family named its closure the resolvent of a
Birkhoff transfer operator on a mathematical billiard (Birkhoff 1927
[#birkhoff_1927]_; Sinai 1970 [#sinai_1970]_; Chernov and Markarian 2006
[#chernov_markarian_2006]_): free flight between collisions with the
boundary, a reflection law
:math:`\psi^+(\Omega') = \alpha\,\psi^-(\Omega)` with
:math:`\Omega' = \Omega - 2(\Omega\cdot\hat n)\hat n`, a transfer operator
:math:`(Sf)(q, v) = \alpha e^{-\tau(q,v)} f(\Phi(q, v))` on the Poincaré
section of the boundary and its Neumann series
:math:`T = (I - S)^{-1} = \sum_n S^n`, rank 1 on a one-wall body and the
anti-diagonal 2 × 2 on a two-wall one. The frame predicted, before any
two-wall code existed, that rank 2 would cover every two-wall body
whatever its curvature and that the impact parameter would separate the
two ranks on a shell; both held on the slab, the hollow sphere and the
annulus. The characteristic reference states the same object without the
tag: :eq:`characteristic-closure` is the least solution of a line's cycle
for every rank, and the rank is derived from the walls the line meets.

**What it established.** `[M]` at the time, each in the deleted gates:
the closed homogeneous sphere at :math:`k_\infty` to :math:`10^{-15}` in
one power iteration (V_α1 reproduced numerically; the cylinder, slab and
shells the same, 4e-16 to 9.3e-16); the vacuum sphere against the
Pomraning–Siewert 1982 Eq. (21) Nyström reference
:cite:`PomraningSiewert1982` to :math:`10^{-4}` at optical radii 2.5 and 5
and secondaries per collision
:math:`c = (\Sigma_s + \nu\Sigma_f)/\Sigma_t` of 0.45 and 0.65 (0.4 and
0.6 counting scattering alone, as the probe log labels them); the
two-group closed bodies at
:math:`k_\infty` to :math:`10^{-9}`; the Issue #132 sphere at
:math:`k \approx 0.735`, above the homogenised :math:`k_\infty` of 0.648
(the fuel is concentrated inside) and below 1, where Phase 4's rank-2
closure read 1.015 and its rank-1 closure 0.551; with no
rank-:math:`N` closure to go wrong (a regression gate against a known
pathology, not a verification); the
three-region Williams 1991 Case 1 fixed-source sphere of Garcia 2021
Table 5 :cite:`Garcia2021` within :math:`1.05 \times 10^{-3}` inside and
:math:`2.3 \times 10^{-2}` at the vacuum surface (after ERR-090); and the
one-group bare cylinder Ua-1-0-CY at the Westfall–Metcalf critical radius
to :math:`4.8 \times 10^{-6}`.

.. _characteristic-origins-history-failures:

What failed, and why
--------------------

Each item is a defect or a misreading, with the measurement that found
it and the mechanism that hid it. The catalogue entries carry the full
record.

- **The first leg was integrated forward.** The B-phase prototype took
  :math:`L_0 = \sqrt{R^2 - r^2(1-\mu^2)} - r\mu`, the forward distance to
  the wall, where the integral form needs the backward one,
  :math:`r\mu + \sqrt{R^2 - r^2(1-\mu^2)}`. V_α1 stayed green, because
  the closure cancels :math:`L_0` on a constant source; the vacuum sphere
  against Pomraning–Siewert found it, 6 % in :math:`k`. A closed-form
  identity on a constant source cannot see the geometry it is fed.
- **V_α2 was a tautology.** The first SymPy proofs of
  :math:`T_{00} = P_{ss}` (sphere and cylinder) built both sides from one
  expression literal, so ``simplify(LHS - RHS) == 0`` held by
  construction; qa caught it in the cylinder's Phase-1 review. Commit
  ``4d87840`` rebuilt each side from its own definition (the
  transfer-matrix integral in :math:`\mu`, the escape integral in
  :math:`\theta'`), added the production-primitive rows, and on the slab
  took the hybrid closed-form and mpmath route. The lesson: two
  derivations are evidence only when they travel different mathematical
  paths, not different code paths around one identity.
- **ERR-034: the slab's position was advanced by the arc length.**
  ``x_traj = x - s`` instead of :math:`x - \mu s`, with the attenuation
  correct. A constant source does not care where on the chord it is read,
  so V_α1 on the closed slab was blind; the vacuum slab's
  :math:`0.45 < k/k_\infty < 0.85` band absorbed a 21 % offset (0.130
  against 0.157 at :math:`\tau_L = 5`); and the slab's vacuum
  self-convergence floor, about :math:`5 \times 10^{-4}`, was written up
  as a slow angular quadrature. The method of images found it, 5 % in
  :math:`k`, because the half slab's mode peaks at the mirror and the
  defect moved the peak to :math:`x \approx 0.16`. After this fix and
  ERR-035's the floor fell 56 times, to :math:`8.85 \times 10^{-6}`, and
  the vacuum slab's :math:`k` rose from about 0.130 to 0.157. The lesson: a
  convergence rate below what the method's theory predicts is a defect's
  fingerprint until shown otherwise.
- **ERR-035: the symmetric slab's closure was built by analogy.** The
  slab was the first body with two walls per period, and its closure was
  taken from the rank-1 sphere by substitution: the per-period
  reflection :math:`\alpha^2` (the keyword ``alpha_per_period``) and the
  out-and-back outflow, :math:`\alpha B_{\rm period}/(1 - \alpha^2
  e^{-2\tau})`. The first-principles rank-2 closure is
  :eq:`peierls-greens-slab-T`, and the identity
  :math:`(1 - \alpha e^{-\tau})(1 + \alpha e^{-\tau}) = 1 - \alpha^2
  e^{-2\tau}` shows where the heuristic's denominator came from; that
  form, like :eq:`peierls-greens-slab-T`, holds only for a source
  symmetric about the mid-plane, and otherwise the general rank-2 cycle
  applies. The two agree at :math:`\alpha \in \{0, 1\}`, which is every
  regime the slab's 22 tests exercised (V_α1 on a constant source,
  vacuum), and differ by :math:`1.3 \times 10^{-4}` relative in
  :math:`k` at :math:`\alpha = 0.5`. The
  asymmetric slab's reduce-to-symmetric row found it; the fix delegated
  the symmetric slab to the rank-2 solver and deleted the heuristic. The
  lesson, now in the ``algebra-of-record`` skill: a closure carried to a
  new geometry by parameter substitution is an algebraic claim that
  needs a derivation on a non-uniform source at an intermediate
  amplitude.
- **ERR-090: one cubic spline across the interfaces.** The multi-region
  oracles read the emission density along each chord from one cubic
  spline through every radial node. The emission jumps at an interface,
  and a smooth cubic across a jump overshoots on both sides with an
  amplitude refinement does not shrink, so the eigenvalue converged in
  :math:`n_r` at about first order and not monotonically: the
  heterogeneous closed sphere (fuel A | moderator B | fuel A, two groups)
  read 1.358083, 1.361371, 1.379031 at :math:`n_r` = 24, 36, 48 with one
  spline, and 1.383737, 1.381917, 1.381293 with one spline per region,
  against the discrete-ordinates solve's recorded 1.38108
  (1.381079639…, ``.claude/plans/reference_p2_spec.md``). It hid behind no radial
  ladder, bands of 2 % to 15 % set from the readings, an interface gate
  whose ceiling was relaxed from :math:`10^{-2}` to
  :math:`5 \times 10^{-2}` to absorb a jump it called a floor, and two
  test defects that agreed with it (a hand-typed stale eigenvalue, and
  cell centres placed uniformly on an equal-volume mesh). The fix was one
  spline per region; on the characteristic reference the class is
  unspellable, because no panel crosses a breakpoint.
- **The angular rules ignored the tangencies.** The chord integrals have
  a square-root kink in angle where a line grazes an interior interface
  (:math:`b = R_k`), at an angle that moves with the node radius, and
  Gauss–Legendre in :math:`\mu` or :math:`\varphi` does not place it. The
  convergence was algebraic and non-monotone: the cylinder's eigenvalue
  read 1.17058, 1.22131, 1.23093, 1.23326, 1.23158 at 8 to 128 azimuthal
  nodes, so no finite ladder bounded its error (#516), and the interface
  continuity gate read jumps from :math:`1.9 \times 10^{-3}` to
  :math:`2.0 \times 10^{-2}` across :math:`n_r` = 24 to 72, catching only
  gross region-indexing errors. The characteristic reference's line rule
  is graded toward every tangency in the impact parameter, with the
  substitution :math:`y = \sqrt{b^2 - r_k^2}` (:ref:`characteristic-line-rule`).
- **ERR-091: one group reported for every group.** ``Billiard``'s
  multi-region sphere fixed-source arm read the group count as
  ``reshape(-1, 1).shape[1]``, which is 1 for every array, and returned
  group 0's flux as the scalar flux. The arm was unreachable until P1
  step 2b of the reference-solution campaign made the multi-region sphere
  constructible, and a one-group source cannot see the defect. The
  lesson: an arm no constructor reaches is untested, and making it
  reachable is when its first success path is reviewed.
- **The method of images held while both sides were off.** `[M]`
  2026-10-10 (the step (e1b) record, ``probes/probe_pred_succ.log``): the
  family's one-group mirror-vacuum slab :math:`[0, 1]` read
  :math:`k = 0.0584228` at its row's resolution; the characteristic
  reference reads 0.0584536527 at rungs 5 and 6, and the Peierls Nyström
  slab :math:`[0, 2]`, an E\ :sub:`1` kernel with no lines, converges onto
  that value (0.05845381, 0.05845368). The family's row asserted the
  identity between its two posings to :math:`10^{-7}` while both were
  :math:`5.3 \times 10^{-4}` low: an agreement of two readings that share
  the reference's error (the ``vv-principles`` skill, anti-pattern #7).
- **The closing-cavity band was below the physics.** The family's
  hollow sphere at :math:`R_{\rm in} = 10^{-3}R` matched the solid sphere
  to :math:`10^{-9}`, and its row asserted :math:`10^{-7}`; a vacuum
  cavity removes a fraction :math:`(R_{\rm in}/R)^2` of the lines, and the
  characteristic reference measures the distance as
  :math:`1.8 \times 10^{-6}` there (:ref:`characteristic-origins-rank-two`).
  `[HYPOTHESIS]` The family under-resolved the few lines through the
  tiny cavity.
- **Pomraning–Siewert does not converge at an optical radius of 25.**
  `[M]` 2026-10-10 (``probes/probe_ps1982.log``): its power iteration
  reaches the 200-iteration cap at 30 and 40 nodes, so the family's
  thick-sphere row, 2e-3 against it, measured that iteration and not the
  family; the row's docstring attributed the gap to the family's cubic
  spline. The successor bands the thick sphere against :math:`k_\infty`
  only and compares with Pomraning–Siewert at optical radii 2.5 and 5.
- **The structure the inventory found.** `[M]` 2026-10-05 (the plan's
  inventory, at ``a336bde4``): the region index, the chord at an impact
  parameter and the sphere's first leg in 3, 3 and 2 copies identical up
  to renaming; the one-region multi-region sphere oracle equal to the
  sphere oracle to :math:`2.3 \times 10^{-16}`; ``compute_resolvent_T_rank2``
  with 0 callers while the rank-2 closure re-spelled its determinant
  inline; :math:`\Sigma_t` passed twice to the cylinder oracle and read
  once; ``closure_rank`` 2 for the asymmetric slab only, although the
  shells close at rank 2; the ``ChordOracle`` protocol used only by
  ``isinstance`` checks. Of 298 tests about the family, 17 were
  delegation tautologies, 6 compared one oracle under two drivers and 27
  asserted that a closed body reads :math:`k_\infty`, which holds for any
  geometry; about 11 compared with an independent value. Every solver
  gate used a symmetric Gauss–Legendre axial rule, blind to a mutation
  that reverses the cylinder's obliquity across the cosines. The optical
  depth of a line through concentric shells was written six times across
  the tree, the cylinder's obliquity three ways and the billiard closure
  five.

.. _characteristic-origins-history-reading:

The family's reference reading
------------------------------

From #405 P2 step 7b.2.3 (``a21b6f8e``) to P1 step (d) of the
characteristic reference campaign (``d9425977``), the multi-region sphere
and cylinder were ``ReferenceSolution`` values that the A|B|A
S\ :sub:`N` rows read, with no certificate. The reading's state was the
emission density at the radial nodes that the final power iterate was
transported from,
:math:`q_g(r_i)/4\pi = [\sum_{g'}\Sigma_{s,g'\to g}\phi_{g'} + \chi_g
\sum_{g'}\nu\Sigma_{f,g'}\phi_{g'}/k]/4\pi`, and a reading at any radius
was its transport by the solver's own operator, divided by the iterate's
fission rate to read it in the solve's gauge,
:math:`\phi_g^{\rm ext}(r) = F^{-1}\int_{4\pi}(K_\alpha q_g/4\pi)\,\mathrm
d\Omega`. That is Atkinson's Nyström interpolation formula
:cite:`Atkinson1997` (eq. 4.1.6), the extension of a nodal solution by
the integral equation itself, chosen over interpolating the nodal flux by
the user's ruling of 2026-10-03, "one reference, one answer"; the
characteristic reference's point reading is the same principle
(:ref:`characteristic-reading`). Three decisions there carry over:

- **One transport.** The chord oracle took the spline's knots and the
  evaluation radii as one argument; a carve (``at=``) separated them, so
  one body served the power iteration and the reading. At the knots the
  extension equalled the solve's final angular flux bit for bit.
- **Split at every tangency.** The angular integral was split where a ray
  grazes a knot or an interface sphere (the sphere at
  :math:`\mu = \pm\sqrt{1 - (\rho/r)^2}`; the cylinder in azimuth, with
  the polar angle in :math:`\theta`). `[M]` (the qa review of step
  7b.2.2): 16, 32 and 64 points per piece moved a sphere point value by
  :math:`2.9 \times 10^{-7}`, :math:`3.4 \times 10^{-8}`,
  :math:`9 \times 10^{-11}`; a missing split no ray crosses gives a
  silently wrong value.
- **Three readings compared.** On the A|B|A shape metric the nodal spline
  (E0), the extension under the solver's rule (E1) and under the split
  rule (E2) differed by :math:`3.2 \times 10^{-4}` to
  :math:`6.4 \times 10^{-3}`; E2 was chosen, and the S\ :sub:`N` shape gap
  moved from :math:`4.361 \times 10^{-3}` to :math:`4.325 \times 10^{-3}`
  on the sphere.

Every reading was ``Uncertified``: the solver derived no bound on its own
error. `[M]` on the Garcia Table 5 sphere its error was about
:math:`2 \times 10^{-4}`, about :math:`10^4` times a rigorous reading's
bound, and every knot of the spline was a singular sphere, so no panel
bound existed for the field as built (#566; the cylinder also #516).

.. _characteristic-origins-history-why:

Why the characteristic reference replaced it
--------------------------------------------

On 2026-10-05, answering two elegance findings on a hoist of the cylinder
oracle's chord, the user ruled: "this machinery is old. check if there is
a better way to architect the general idea before patching old code"
(the plan ``.claude/plans/characteristic_reference_architecture.md``,
"The ruling that opened this plan"). The audit that followed found the
family's objects to be instances of a few that were each written several
times, and the characteristic reference is those objects, once:

- **The geometry** is the geometric kernel's (:ref:`theory-chart-and-chord`):
  one crossing law per chart, segment lengths without cancellation,
  regions read from the crossing order, the obliquity as one factor for
  three charts. The family's seven oracle classes and their copied
  helpers have no successor of their own.
- **The rank is derived**, from the period of the unfolded line through
  the walls' partners (:eq:`characteristic-transit-rank`), where the
  family tagged it (``geometry_kind``, ``closure_rank``) and the shells
  branched on :math:`b` inside their oracles.
- **One closure for every rank** (:eq:`characteristic-closure`), its
  rank-1 and rank-2 cases the equations of this page, with white walls
  added as a finite-rank update the family did not have
  (:ref:`characteristic-resolvent`).
- **The emission lives on a panel basis per region**, so ERR-090's class
  cannot be written, and the transport is a Galerkin assembly over lines
  whose eigenvalue is the dense pencil's (:eq:`characteristic-pencil`):
  no power iteration, no initial guess of :math:`k`, no convergence flag.
- **The line rule is graded from the group's optical scale**, toward
  walls, interfaces, tangencies and grazing, where the family used fixed
  Gauss–Legendre counts that its ladders could not bound.
- **One reference posed from the specification**,
  ``CharacteristicDerivation(specification, resolution)``, answers the
  eigen, source and response questions with flux integrals and point
  values (:ref:`characteristic-door`), where the family had fifteen entry
  points and a facade.

The family was built rung by rung beside its successor, compared with it
on the S\ :sub:`N` fixtures (P1 step (c)), replaced as the S\ :sub:`N`
rows' reference (step (d)), its tests re-posed (step (e1b)) and deleted
(step (e2)).

.. rubric:: References for the billiard frame

.. [#birkhoff_1927] Birkhoff, G.D. (1927). *Dynamical Systems*,
   Ch. 6 (billiard flows). Reissued by AMS Colloquium Publications,
   1966.
.. [#sinai_1970] Sinai, Ya.G. (1970). "Dynamical systems with elastic
   reflections and their applications." *Russ. Math. Surveys* **25**,
   137–189. DOI: 10.1070/RM1970v025n02ABEH003794.
.. [#chernov_markarian_2006] Chernov, N. and Markarian, R. (2006).
   *Chaotic Billiards*. AMS Mathematical Surveys and Monographs
   **127**.

Changelog
---------

.. list-table::
   :header-rows: 1
   :widths: 14 58 14 14

   * - Date
     - Decision
     - Commit
     - Issue
   * - 2026-04-28
     - Phase 5 closed: the angle-integrated Sanchez kernel is
       hypersingular at the diagonal; the angle-resolved Green's function
       (Variant α) is adopted as a parallel reference for the closed
       sphere.
     - ``4dc03cf2``
     - #133
   * - 2026-05-02
     - The family grows to six geometries on two orbit-space classes: the
       vacuum and multi-group sphere, the multi-region sphere (k and fixed
       source), the cylinder, the shared closure module, the slab, the
       asymmetric slab with the rank-2 closure, the hollow sphere and the
       annulus. ERR-034 and ERR-035 found and fixed; V_α2 rebuilt on
       independent paths.
     - ``efbae9c``, ``92d4f10``, ``166b9ae``, ``7cb5bc6``, ``4d87840``
     - #129, #132
   * - 2026-05-12
     - The multi-region cylinder (Phase 1b); the composite per-region
       radial quadrature.
     - ``37e3e299``, ``2d3e7f2d``
     - #168
   * - 2026-09-26
     - ERR-090: one spline per region for the emission density; the
       angular rules' tangency kinks left open.
     - ``6bcea45c``
     - #516
   * - 2026-10-03
     - The multi-region sphere and cylinder become uncertified reference
       solutions read by the natural extension of their emission density.
     - ``a21b6f8e``
     - #405, #566
   * - 2026-10-10
     - P1 step (e) of the characteristic reference campaign: the SymPy
       derivations move to ``characteristic/origins/`` (e1a); the family's
       tests are re-posed on the characteristic reference (e1b); the
       family is deleted (e2). This page replaces the family's page and
       keeps the 36 labels test markers name; the V_α2 rows' markers move
       off the closure labels, and ``peierls-greens-cylinder-V-alpha-2`` is
       added for the cylinder's row.
     - ``f66c1c45``, ``44303919``, ``5aa8cb88``
     - #405, #572
