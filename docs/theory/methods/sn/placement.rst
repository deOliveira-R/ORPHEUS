.. _sn-placement:

Placement: why discrete ordinates
=================================

Every method book in the canon opens by deriving its method; none opens
by justifying it. The five method families sit in five consecutive
silos — Hébert's chapters :cite:`Hebert2009` derive P\ :sub:`N`, CP,
S\ :sub:`N` and MoC each from scratch, and the only comparative language
in the entire span compares Monte Carlo *estimators to each other* — so
the question an engineer actually faces first, *"when do I reach for*
S\ :sub:`N` *instead of CP, or MoC, or diffusion?"*, has no chapter
anywhere. The review that comes closest, Sanchez & McCormick
:cite:`Sanchez1982`, is organized by *formulation* — the approximation
axis, not the named-method axis — and its introduction is the canon's one
genuinely comparative page: it states the selection criteria (the
spatial and angular detail required, the scattering anisotropy, the
heterogeneity, the optical size of the regions — and, "perhaps, most
importantly," it concedes, "which computer programs are readily
available!") and the formulation-level rule: integro-differential
forms for optically large media, integral forms for optically thin
ones. But it never descends from formulations to a method-level
decision, and the review ends with an appendix — no concluding
assessment exists. This chapter is the missing chapter, for this book:
the engineering trade-space of discrete ordinates, built on that
review's criteria and the corpus root's axes.

The structural half of the answer is already settled one level up and
is **not repeated here**: the corpus root
(:doc:`/theory/foundations/path_integral`) locates every method on
three independent axes — S\ :sub:`N` realizes the propagator
:math:`(L+C)^{-1}` as a cell-DAG :term:`sweep` with
rational-approximant transmission (A1), resums scattering by outer
Neumann iteration or Krylov (A2), and represents angle by
:term:`ordinates <ordinate>` (A3) — the axes are defined at
:ref:`path-integral-axes`, and the landing of every method on them at
:ref:`path-integral-method-map`. This chapter answers the question
the axes cannot: **when is that point in the product space the right
place to stand** — what the position costs, what it buys, and which
problem regimes reward it.

.. admonition:: Key Facts
   :class: tip

   * **What S**\ :sub:`N` **buys**: the resolved :term:`angular flux`
     itself (not just its moments), on a mesh, at deterministic cost —
     a **matrix-free direct solve** per direction (the sweep is forward
     substitution, never an assembled matrix), anisotropic scattering
     to arbitrary Legendre order, and the canon's most mature
     acceleration theory (DSA, Krylov) :cite:`AdamsLarsen2002`.
   * **What it pays**: **ray effects** — the canon's own verdict is
     "the most significant of all S\ :sub:`N` deficiencies," with
     "very little … accomplished during the past 40 years"
     :cite:`LarsenMorel2010`; angular storage
     :math:`O(N_{\mathrm{cells}} N_{\mathrm{ord}} N_g)`; a sweep
     schedule that wants structured meshes; and the multigroup data
     dependence every deterministic method shares.
   * **The decision in one line**: reach for S\ :sub:`N` when you need
     the angular flux resolved *everywhere* on a structured or
     curvilinear mesh at deterministic cost — exact-geometry lattice
     work goes to MoC/CP, whole-core few-group work goes to
     diffusion/nodal methods, and data-exact or geometry-exact
     reference work goes to Monte Carlo.
   * **S**\ :sub:`N` :math:`\equiv` **P**\ :sub:`N−1` in slab geometry
     with :math:`N`-point Gauss quadrature — ordinates versus
     harmonics is a *representation* choice on the root's A3 axis, not
     an accuracy chasm; the real differences are operational
     (sweepability, positivity, ray effects versus Gibbs
     oscillations).

What S\ :sub:`N` buys, and what it pays
---------------------------------------

The identity that drives every entry in the trade-space: S\ :sub:`N`
**retains the angular flux** :math:`\psi(\mathbf r, \hat\Omega, E)` as
its unknown (:ref:`sn-synopsis`), and it pays for that resolution with
a *collocation* in angle. Everything S\ :sub:`N` is good at, and
everything it suffers from, follows from those two facts.

**The buy side.** Because angle is collocated rather than coupled, the
within-group operator decouples direction by direction in Cartesian
geometry, and each direction's spatial system is triangular under the
sweep order. (Curvilinear geometry couples ordinates through the
angular-redistribution term — yet the *augmented* space-and-angle walk
remains triangular, so the economics survive the coupling;
:doc:`curvilinear_one_group`.) Either way the
propagator is applied by **forward substitution, with no assembled
matrix**, at :math:`O(N_{\mathrm{cells}} N_{\mathrm{ord}} N_g)` per
source iteration (the root's A1 cell; certified per (mesh, closure,
boundary) triple in this book,
:doc:`loss_representation`). That is the cheapest exact
propagator-application in the deterministic family, and the structural
reason is one Sanchez & McCormick state as a matter of *balance
locality* :cite:`Sanchez1982`: the integro-differential form rests on a
**local** neutron balance, so its matrices are sparse, cheaply
assembled, and solvable with only a fragment in memory at a time — in
S\ :sub:`N`'s case, never assembled at all. It is also the reason the
iteration/acceleration literature — spectral radius
:math:`\rho = c`, DSA's :math:`0.2247c`, the Krylov rescue — is
*built on* S\ :sub:`N` :cite:`AdamsLarsen2002`: the sweep is the
preconditioner every scheme wraps. Anisotropic scattering enters
through Legendre moments at any order the data carries, with the
Galerkin frame making the moment map exactly invertible on Gauss
quadratures :cite:`LarsenMorel2010`. And the answer is a *field*:
every cell, every direction, every group, in one solve — the
substrate that perturbation theory, adjoint weighting
(:doc:`adjoint`), and reaction-rate functionals consume directly.

**The pay side.** Ray effects are the structural price of collocation
— in the 1982 review's words, the angular discretization couples the
spatial and angular approximations strongly enough "to produce
space-angular nuisances such as the ray effect" :cite:`Sanchez1982`:
a localized source in a weakly collided medium emits into a continuum
of directions, but the discrete set can only transport along its
ordinates, so the computed flux develops unphysical ridges aligned
with the quadrature directions. The failure regime is precise —
optically thin, absorbing media with localized sources (ducts, voids,
detector ports) — and the canon's forty-year verdict stands: quadrature
refinement converges slowly, P\ :sub:`N` is *not* a cure ("they rarely
provide the correct solutions to problems for which S\ :sub:`N`
solutions exhibit ray effects"), and the probable future is angular
adaptation :cite:`LarsenMorel2010`. The second price is memory: keeping
:math:`\psi` costs :math:`O(N_{\mathrm{cells}} N_{\mathrm{ord}} N_g)`
storage, which is exactly why this book's sweep machinery ships a
rolling-window sweep mode that avoids storing the full angular field
unless a consumer asks for it (:doc:`loss_representation`). The third is the sweep
schedule itself: the triangularity that makes the solve cheap is a
*theorem about the mesh* (:ref:`path-integral-method-map`), and
unstructured or cyclically-coupled configurations must first be
decomposed before they sweep.

The trade-space, method by method
---------------------------------

Each comparison below states what the *other* method's position on the
root's axes buys it, where it beats S\ :sub:`N`, where it loses, and
the regime where the choice flips.

Versus collision probability — the integral sibling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CP stands at the opposite pole of the angular axis: it **integrates
angle out** before discretizing anything, collapsing the transport
problem to a scalar-flux integral equation whose kernel — the
collision probabilities — carries the exact streaming attenuation
between every region pair (:doc:`/theory/methods/collision_probability`).
What that buys: *no ray effects at all* (angle was never collocated),
*no sweep ordering* (the kernel is assembled, not marched), exactness
in angle within its assumptions, and its **white (isotropic-return)
boundary** folded into the kernel in closed form — the rank-1
Sherman–Morrison update the root page names as the move S\ :sub:`N`
has yet to borrow (:ref:`path-integral-method-map`, Issue #300). What
it pays: that closed-form boundary is itself an approximation — an
isotropic re-entry standing in for a true reflective edge, a gap the
CP chapter quantifies against reflective S\ :sub:`N`, whose angular
trace reflects *exactly*; the
kernel is **dense** — the 1982 review derives this from balance
locality, the mirror of S\ :sub:`N`'s sparsity: the integral equation
rests on a *global* balance along each direction, so it is "strongly
coupled," its full matrices are assembled by "expensive evaluation of
transcendental functions," and the complete matrix must be held and
solved globally :cite:`Sanchez1982`. The coupling grows as
:math:`O(N_{\mathrm{reg}}^2)` per group — a hard domain-size ceiling.
The angular exactness is likewise conditional: it holds "provided the
scattering anisotropy is low (isotropic or linearly anisotropic)," and
in the integral formulation the equation count "dramatically
increases with the degree of anisotropy" — whereas the
integro-differential family absorbs arbitrary anisotropy by modifying
the collision term alone :cite:`Sanchez1982`. And the exact geometry
is bought per configuration: each new geometry class needs its own
specialized integration kernel, where a mesh-based method approximates
*any* configuration with enough zones.

The flip is domain size and angular detail — the review's own
formulation-level rule, integro-differential for optically large
media, integral for optically thin ones :cite:`Sanchez1982`, is this
flip stated at the family level. Small, geometrically intricate,
scattering-dominated domains — the classical lattice cell — reward
CP's exact kernel, forgive its dense coupling, and rarely need more
than the scalar flux CP produces directly ("which usually is all that
is needed," the review notes); within the integral family the same
review crowns CP "by far the most used" of the variants
:cite:`Sanchez1982`. This is the WIMS tradition :cite:`WIMSD` and the
workhorse role Hébert's lattice codes give it :cite:`Hebert2009`.
Large domains, strong anisotropy, or any question whose *answer* is
angular (streaming through a duct, interface currents, ray-level
detail) flip the choice to S\ :sub:`N`.

Versus the method of characteristics — the same-side sibling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The root page proves S\ :sub:`N` and MoC are **one method family**:
both collocate angle, both realize the propagator exactly on their
grid, and the spatial closures of the sweeping methods are rational
approximants of the very exponential :math:`e^{-\tau}` that MoC
integrates exactly — one Padé table, different entries
(:ref:`path-integral-method-map`). The placement question between them
is therefore *not* accuracy doctrine; it is **geometry versus cost**.
MoC buys arbitrary geometry: its rays are traced through the true cell
shapes :cite:`Askew1972`, which is why modern lattice physics chose it
— exact pin, clad, and moderator boundaries with no stair-stepping
:cite:`KnottYamamoto2010`. The 1982 review already saw the shape of
this settlement: it names the method of characteristics as *the*
example of the "hybrid" trend — integro-differential transport inside
each region, integral linkage of subregions through their interface
angular fluxes — whose twofold advantage it states precisely: streaming
is well approximated (so subregions can be larger than an
integro-differential mesh would allow), and interface-only coupling
keeps the matrices sparse and iteration-friendly :cite:`Sanchez1982`. It pays in tracking: segment storage and
per-sweep work scale with the track discretization, and the cost of
carrying that machinery to full 3-D is the standing reason industrial
MoC lives in 2-D planes fused with axial low-order methods. S\
:sub:`N` buys the opposite: on structured Cartesian and 1-D
curvilinear meshes the sweep needs *no stored geometry beyond the mesh
itself*, extends to 3-D at linear cost in cells, and — through the
:math:`\alpha`-redistribution machinery this book derives
(:doc:`curvilinear_one_group`) — handles curved coordinate systems
that ray tracing handles only by brute force.

The flip: exact 2-D geometry → MoC; structured meshes, 3-D
deterministic work, or curvilinear 1-D → S\ :sub:`N`. When the mesh
*is* structured, the two methods' answers converge to each other by
construction — they share the propagator, differing only in the Padé
entry — which is why this corpus treats MoC as S\ :sub:`N`'s
verification partner on shared idealizations rather than a rival
(:doc:`/theory/methods/method_of_characteristics`).

Versus P\ :sub:`N` — the other angular representation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

P\ :sub:`N` differs from S\ :sub:`N` on exactly one axis — A3: angle by
global spherical-harmonics expansion instead of collocation — and the
slab-geometry equivalence **S**\ :sub:`N` :math:`\equiv`
**P**\ :sub:`N−1` (with :math:`N`-point Gauss quadrature) shows how
thin that difference is at the level of *content*: the S\ :sub:`N`
scattering matrix is a similarity transform of the Legendre one, and
the two methods carry the same angular information
:cite:`LarsenMorel2010`. The differences are operational, and they
all trace to *support*: harmonics are global in angle, ordinates are
local.

- A global basis has **no sweep**: the P\ :sub:`N` moment equations
  couple all moments through the streaming term, producing coupled
  elliptic-flavored systems rather than per-direction triangular ones
  — the propagator is not realized (the root's A1 column marks the
  family "—"), and the solve machinery is a different world.
- A global basis **cannot ray-effect** — its artifacts are the dual
  failure: **Gibbs oscillations** at angular discontinuities (beam
  boundaries, material interfaces, voids), where truncated harmonic
  expansions ring and can drive the reconstructed flux negative.
- A local basis is **positivity-friendlier**: each ordinate's flux is
  a physical flux, and closures control its sign cell by cell
  (:ref:`path-integral-method-map`); a truncated moment set has no
  per-direction sign to control.

The flip is the angular character of the solution: smooth, diffusive,
near-isotropic transport rewards P\ :sub:`N` (and its industrial
simplification SP\ :sub:`N`, the whole-core compromise —
:doc:`/theory/references/spn_method`); transport with fronts, beams,
and streaming paths rewards ordinates. In slab geometry the choice is
provably cosmetic; in multi-D the operational differences decide.

Versus diffusion — the limit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Diffusion is not a competitor on the same footing — the root page is
precise about this: it is the **one method that is a limit of the
object rather than a quadrature of it**
(:ref:`path-integral-streaming`). The placement question is therefore
*when the limit is enough*: optically thick, scattering-dominated
regions, gradients mild on the scale of a mean free path, and a few
mean free paths away from boundaries, sources, and streaming channels
— the classical validity conditions
:cite:`BellGlasstone1970,Duderstadt1976`. Where they hold, diffusion
is orders of magnitude cheaper than any transport method and just as
accurate; whole-core steady-state analysis lives there, in few-group
nodal diffusion fed by lattice-homogenized data — the two-step
paradigm :cite:`Stacey2007,KnottYamamoto2010`.

Two structural bridges keep S\ :sub:`N` and diffusion from being mere
neighbors:

1. **The thick diffusion limit** is a property a *spatial closure* can
   possess: an S\ :sub:`N` scheme that has it produces accurate
   diffusive answers even on cells many mean free paths thick —
   provided the cells resolve the *diffusion* length — so a
   limit-possessing S\ :sub:`N` degrades gracefully *into* diffusion
   exactly where diffusion is valid :cite:`LarsenMorel2010`. The
   scheme verdicts (which closures pass, which fail) are part of this
   book's closure story (:doc:`/theory/foundations/discretization`).
2. **DSA** makes diffusion S\ :sub:`N`'s *inner partner*: the
   diffusion operator is the natural preconditioner for the scattering
   iteration precisely because it captures the modes the sweep
   attenuates worst :cite:`AdamsLarsen2002` — so in production
   S\ :sub:`N`, diffusion is not the alternative but a component
   (Issue #2 tracks ORPHEUS's DSA).

The flip: if the whole problem satisfies the validity conditions, use
diffusion and be done; S\ :sub:`N` earns its cost where the conditions
*break* — near strong absorbers, voids, control elements, boundaries —
and in the reference role, checking the diffusion answer's error
budget from above (:doc:`/theory/methods/diffusion_1d`).

Versus Monte Carlo — the other splitting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Monte Carlo differs at the *deepest* branch the root page names — the
generator splitting itself (:ref:`path-integral-generator-splitting`):
the collision is read as a jump in a sampled process, not a source in
an iterated series, and no quantity is discretized at all. The error
anatomies are therefore disjoint, and the comparison is genuinely
complementary rather than competitive:

- **MC's buy**: continuous-energy data (no multigroup approximation —
  the one bias *every* deterministic method shares and MC alone
  escapes), exact geometry, exact angle (no ray effects), and a
  statistical error bar that comes with the answer
  :cite:`Lux1991`. ORPHEUS's MC rides the majorized-jump splitting —
  Woodcock delta tracking :cite:`Woodcock1965` — whose exactness the
  root page states as a change of measure, not an approximation.
- **MC's price**: the answer is an *estimate at tallied locations*,
  not a field — global detail costs variance everywhere, rare-event
  responses (deep penetration, small detectors) need variance
  reduction to converge at all, and derivative/adjoint information
  requires dedicated machinery rather than falling out of the solve.
- **S**\ :sub:`N`'s **counter**: the full deterministic field —
  every cell, direction, and group at once, with adjoints and
  perturbation chains as linear-algebra corollaries — at the price of
  discretization bias in angle, space, and energy.

The flip follows the question's shape: *"what is the response, with
what confidence?"* → MC; *"what does the field look like, and how
does it move when the data moves?"* → S\ :sub:`N`. In this corpus the
two are verification partners: a converged S\ :sub:`N` sweep and an
MC track-length tally estimate **the same first moment** — the root
page's opening identity — so their agreement at fixed data is a
cross-splitting check that neither method can run alone.

The regime map
--------------

The trade-space, folded to a table. "First reach" is the default the
regime rewards; S\ :sub:`N`'s role is stated honestly, including where
that role is *secondary*.

.. list-table::
   :header-rows: 1
   :widths: 24 18 32 26

   * - Problem regime
     - First reach
     - Why
     - S\ :sub:`N`'s role
   * - Lattice / assembly spectral calculation (2-D, exact geometry,
       fine groups)
     - MoC (modern), CP (classical)
     - exact pin geometry (MoC); cheap exact kernel on a small domain
       (CP)
     - reference partner on structured idealizations of the same
       problem
   * - Whole-core steady state (3-D, few groups)
     - nodal diffusion / SP\ :sub:`N`
     - the diffusion limit is valid over most of the core, at orders
       less cost
     - the transport reference that prices the low-order model's
       error; the DSA partner inside production S\ :sub:`N`
   * - Deep-penetration shielding
     - S\ :sub:`N`, or MC with variance reduction
     - a deterministic field beats vanishing tally statistics — but
       ray effects are S\ :sub:`N`'s own duct/void hazard
     - the primary deterministic tool
   * - Streaming-dominated media (ducts, voids, beams)
     - MC, or S\ :sub:`N` with high order / adaptation
     - the exact-angle representation has no ray effects
     - usable with care; this is S\ :sub:`N`'s characteristic failure
       regime
   * - Curvilinear 1-D benchmarks and method study
     - S\ :sub:`N`
     - the angular-redistribution machinery makes spheres and
       cylinders native; exhaustive parameter scans are cheap
     - this book's home ground (:doc:`curvilinear_one_group`)
   * - Data-exact / geometry-exact reference
     - MC
     - continuous energy and exact geometry — the biases the
       deterministic family cannot shed
     - the cross-splitting verification partner
   * - Exact-in-angle reference values
     - Case / F\ :sub:`N` reference solvers
     - closed-form resummation of the collision series
     - consumer of their benchmark values
       (:doc:`/theory/references/index`)

Where ORPHEUS stands
--------------------

This corpus does not merely *describe* the trade-space — it ships
every side of it over one operator algebra, which changes what the
comparisons mean. Because the reaction operators are shared code
(:ref:`path-integral-invariant`) and one cross-section pipeline feeds
every solver, a cross-method comparison in ORPHEUS is *at fixed
multigroup data by construction* — the scope condition that makes
"the methods agree" a statement about methods rather than about data.
S\ :sub:`N` is the flagship realization (this book, taught through
the broadening progression); CP is the production integral sibling;
diffusion is the landed limit sharing the same operator leaves; MC is
the Woodcock sampler of the same first moment; the Case-family
reference solvers pin the exact-in-angle corner. The placement
decisions this chapter describes are therefore *runnable experiments*
here, not literature claims.

One routing note closes the chapter: this page decides *whether*
S\ :sub:`N` is the right tool **before** you commit; once you are
inside the method and something misbehaves, the router's
symptom-to-chapter table (:ref:`sn-symptom-table`) is the map that
takes over.
