.. _verification-reference-architecture:

=============================================================================
Readings, claims and certificates: the reference-solution architecture
=============================================================================

.. Machine header. Owner of: the ontology of a reading, a claim and a
   certificate (#405 P2); the readings by role; the observables as read
   off an answer; how a reference's enclosure is established and why no
   ladder certifies; the reference certificate's derived state; the
   reference solution read on demand; the verification certificate, the
   order verification and their floors; the uncertified reading and its
   explicit comparison; the chain-of-problems error ontology (the user's
   G6 draft and its literature). Not owned here: the binding
   vocabulary and operator-form contract (reference_solutions.rst), the
   doctrine and the ladder (principles.rst), the API reference
   (docs/api/reference.rst), the trajectory resolvent's reading
   (references/characteristic_origins.rst,
   characteristic-origins-history-reading). Plan of record:
   .claude/plans/reference_cache.md; specification of the gates:
   .claude/plans/reference_p2_spec.md. Written at P2's close,
   2026-10-03, against branch feature/reference-uncertified-reading at
   a21b6f8e.

A verification compares two numbers, and the comparison is only as good as
what each number claims. Most of the suite compares them inline,
``abs(k - k_ref) < tol`` (`[M]` 38 lines in 22 test files, the P2 census
of 2026-10-02), against a reference value that carries no statement of
its own accuracy. A float cannot say whether it is exact, a
printed table entry, a converged iterate or an estimate, and the test that
reads it cannot either. This chapter documents the objects that now carry
that statement: what one reads off an answer (the **observable**), what the
answer says about it (the **reading**), the claim a reference makes about
the exact answer of its own equation (the **claim**, established exactly or
by a derived bound), the evidence that can refute such a claim (the
**reference certificate**), and the comparison of a production reading
with a reference that may stand as one (the **verification certificate**).
Its last part is the error ontology the design rests on: an answer is a
node on a chain of problems that starts at the specification, and every
error is one link of that chain.

.. contents::
   :local:
   :depth: 2


Key facts
=========

* **A reading is typed by whose claim it is** (the user's ruling "G3",
  2026-10-02). A reference's reading is a guarantee about the exact answer
  of its own equation, an :class:`~orpheus.numerics.enclosure.Enclosure`, a
  cited :class:`~orpheus.reference.reading.Printed` value, or, where its
  family cannot yet derive a bound, an
  :class:`~orpheus.reference.reading.Uncertified` value. A production
  method's reading is a self-report that verification puts on trial,
  :class:`~orpheus.numerics.outcome.Measured`. The two sums share no member.
* **One closed set of observables**, read by every answer through one verb,
  ``answer.read(observable)``: a flux integral (the primitive), a ratio of
  two linear observables, the eigenvalue, a point value. The answer is the
  receiver because the answers grow and the observables do not (ruling
  "G1").
* **A reference's enclosure is established in two ways only**:
  :class:`~orpheus.reference.certificate.Exact` (a symbolic value and the
  self-check that proves it) and
  :class:`~orpheus.reference.certificate.DerivedBound` (an a posteriori
  bound with computable constants). **No refinement ladder certifies**
  (the user's ruling of 2026-10-03): a ladder can refute a reference,
  never certify one, because a systematic bias is invisible to it.
* **The reference certificate's state is derived, never stored**:
  ``Valid``, ``Invalid`` with every failed check's reason, or the
  ``Withdrawal`` itself. Only a reference whose certificate is ``Valid`` may
  anchor a verification.
* **A reference is read on demand through its natural extension**, one
  answer per reference (for an integral-equation reference, the transport
  integral of its emission density), never a second nodal answer, and P2
  stores no field: certification is decoupled from caching.
* **Verification has two floors**: the reference's bound is at most a tenth
  of the tolerance, and production's algebraic error is small against the
  tolerance (recorded under an agreement verdict, required under an order
  verdict).
* **Honest scope.** P2 ships exactly one certified reference family, the
  exact infinite medium; the trajectory-resolvent sphere and cylinder
  were reference solutions with no certificate, so every row that read
  them was an explicit uncertified comparison or a strict xfail on the
  verbs' own refusal, and none was a verification claim. Since P1 step (e2)
  of #405 (2026-10-10) the family is deleted and its successor, the
  characteristic reference, is uncertified in the same way. Deriving each
  family's bound is phase P4's work (#566). Every other
  reference in the tree is still a
  :class:`~orpheus.derivations.ContinuousReferenceSolution` (`[M]` 19
  construction sites under ``orpheus/``, 2026-10-03, against 2 of
  ``ReferenceSolution``), governed by the contract of
  :doc:`reference_solutions` until P4 migrates it.


The ontology: an answer, an observable, a reading, a claim, a certificate
=========================================================================

The objects, from the bottom up, with the package each lives in (the
layer rule is :ref:`architecture-layering`: physics-free values in
:mod:`orpheus.numerics`, anything naming a specification or a citation in
the input-tier :mod:`orpheus.reference`, which the derivations import and
never the converse).

.. list-table:: The objects of the reference architecture
   :header-rows: 1
   :widths: 24 46 30

   * - Object
     - What it is
     - Home
   * - observable
     - one number to read off an answer, the same number whichever answer
       is asked; a closed sum of four members
     - :mod:`orpheus.numerics.observable`
   * - answer
     - anything with a ``read(observable)`` verb: a production solution, a
       reference solution, a published table
     - each method's own package
   * - reading
     - what one answer says about one observable, typed by whose claim it is
     - :mod:`orpheus.numerics.enclosure`, :mod:`orpheus.reference.reading`,
       :mod:`orpheus.numerics.outcome`
   * - establishment
     - how a reference's enclosure of one observable is known: exactly, or by
       a derived bound
     - :mod:`orpheus.reference.certificate`
   * - claim
     - one observable's declared target and its establishment
     - :mod:`orpheus.reference.certificate`
   * - reference certificate
     - a reference's claims, the evidence that can refute them, its standing,
       and its derived state
     - :mod:`orpheus.reference.certificate`
   * - published solution
     - the values a cited work prints, per observable, each with its own
       citation
     - :mod:`orpheus.reference.published`
   * - reference solution
     - a reference method's answer: a specification, a derivation evaluated
       on demand, and a certificate where the family has one
     - :mod:`orpheus.reference.solution`
   * - verification certificate
     - one comparison of a production reading with a valid reference's
       enclosure, and its verdict; returned, never stored
     - :mod:`orpheus.reference.verification`

**Why "certificate" names two objects and no more** (the user's ruling
G4). A certificate is a guarantee someone other than the answer's author
can check: the reference certificate (what a reference claims and the
evidence that can refute it) and the verification certificate (one
comparison and its verdict). Production's exit record is not one: it is
the solver's report on its own discrete solve,
:class:`~orpheus.numerics.outcome.ExitReport`, and the verification
certificate reads it as an input. For the same reason the member of the
evidence sum that records a bound the convergence-claim check asserted is
:class:`~orpheus.numerics.outcome.Asserted`. A group action's
``OrbitCertificate`` keeps its name: it certifies a symmetry, a different
subject.


Readings by role
================

A reading is a claim about one observable, and whose claim it is decides
its type. The user's ruling G3 (2026-10-02): a **reference's** claim is a
guarantee about the exact answer of its own equation; a **production**
method's claim is a self-report that verification puts on trial. The two
are separate sums, so a production reading cannot be passed where a
reference's claim is required.

.. list-table:: The readings
   :header-rows: 1
   :widths: 18 18 34 30

   * - Reading
     - Whose claim
     - What it guarantees
     - Can it anchor a verification?
   * - ``Enclosure(value, bound)``
     - a reference
     - the exact answer lies in :math:`[v - b, v + b]`
     - yes, when its reference's certificate is ``Valid``
   * - ``Printed(text, citation)``
     - a cited author
     - nothing ORPHEUS computed; the author's printed precision, with the
       place it is printed
     - never directly: it is an anchor of a reference certificate
   * - ``Uncertified(value)``
     - a reference whose family derives no bound
     - nothing
     - no; it is compared explicitly, never verified against
   * - ``Measured(value)``
     - a production method
     - nothing; it is the thing on trial
     - it is the other side

The reference sum is
:data:`~orpheus.reference.reading.ReferenceReading`
``= Enclosure | Printed | Uncertified``; the production sum is
:data:`~orpheus.numerics.outcome.ProductionReading` ``= Measured``, the
same class the iteration record uses, defined once. A statistical
interval (Monte Carlo's ``Estimated``, not built) would join the
production sum only, so an ``Estimated`` reference is unspellable twice
over: Monte Carlo is never a reference (#505), and the type that would
carry its interval lives on the other side. A stochastic reference would
need a finite-sample bound (a Hoeffding or Bernstein inequality, with the
score bounded a priori and the failure probability declared); that is
admissible in principle and not built.

The enclosure
-------------

``Enclosure(value, bound)`` claims that an exact real number :math:`x`
satisfies :math:`|x - v| \le b`, in exact real arithmetic on the two
doubles. An exact value is the enclosure with bound 0. A value enters as
an enclosure once, at the place it is known exactly:
:meth:`~orpheus.numerics.enclosure.Enclosure.about` takes a rational
:math:`x` and a rational radius :math:`\rho` and returns

.. math::

   v = \mathrm{fl}(x), \qquad
   b = \mathrm{up}\!\left(\rho + \lvert \mathrm{fl}(x) - x\rvert\right),

where :math:`\mathrm{fl}` rounds to the nearest double and
:math:`\mathrm{up}` rounds up to a double: the sum is computed exactly in
rationals, so the enclosure contains the whole closed ball. A bare float
never divides an enclosure; the only way a number becomes a claim is this
one door.

**The quotient.** A ratio observable reads as the quotient of its operands'
readings, so the quotient of two enclosures must enclose every quotient of
their members. On a box that excludes :math:`y = 0` the map
:math:`(x, y) \mapsto x / y` is monotone in each argument, so its extremes
over the box are among the four corner quotients. With
:math:`[a_-, a_+]` and :math:`[b_-, b_+]` the outward-rounded ends of the
two operands,

.. math::
   :label: reference-enclosure-quotient

   \frac{[a]}{[b]} = \mathrm{Enclosure}\Bigl(\frac{v_a}{v_b},\;
     \max\bigl(\mathrm{up}(c_+ - v_a/v_b),\ \mathrm{up}(v_a/v_b - c_-)\bigr)\Bigr),
   \quad
   c_- = \mathrm{down}\Bigl(\min_{\alpha, \beta \in \{-, +\}} \frac{a_\alpha}{b_\beta}\Bigr),
   \quad
   c_+ = \mathrm{up}\Bigl(\max_{\alpha, \beta \in \{-, +\}} \frac{a_\alpha}{b_\beta}\Bigr),

.. vv-status: reference-enclosure-quotient documented

.. (vv-status rationale) definitional: the enclosure quotient as built by
   Enclosure.__truediv__; its laws (containment, the centre, tightness,
   the refusals) are pinned by the foundation rows R1.2-R1.5 of
   tests/gates/numerics/test_enclosure.py, which by the harness contract
   carry no verifies marker.

where every rounded operation is followed by one outward step
(``math.nextafter``), which covers its half-unit error. The centre is the
quotient of the centres, **bit for bit**, so a ratio's reading is the ratio
of its readings (a midpoint centre would move with the operands' bounds).
A denominator whose rounded enclosure holds zero has no quotient
(``ZeroDivisionError``), and an overflow is refused rather than returned as
an infinite bound.

**Why tightness is a law, not a courtesy.** A bound of :math:`10^{300}`
contains everything and is valid; it would also fail every verification
floor. So the gates hold two sides: containment (R1.2: `[M]` 3010 of 3010
population rows, including 3000 seeded draws with magnitudes from
:math:`10^{-30}` to :math:`10^{31}`, contain the exact rational quotient
at all four corners) and tightness (R1.4). The tightness bound is derived:
each operand end carries one rounding and one outward step (3 u relative,
u the unit roundoff), a corner quotient is perturbed by 6 u, rounded
(u) and stepped (2 u), 9 ulp in all, and the half-width's subtraction is
rounded and stepped, 1.5 ulp more: **10.5 ulp** of the largest corner.
`[M]` the worst over the population is 7.84 ulp (2026-10-03,
``scratch/reference_architecture/p2/ta/probe_tight.py``). The linearised
bound :math:`(b_a\lvert v_b\rvert + b_b\lvert v_a\rvert)/v_b^2` is not an
enclosure: it under-encloses whenever the denominator's relative width is
not small (`[M]` the battery arm that installs it reddens the containment
law).

The printed value
-----------------

A value printed in a cited work is its **decimal text**:
``Printed(text, citation)``. The trailing zeros a float drops are the
precision claim (``1.0`` and ``1.00`` are two claims), so the text is the
one source; :class:`~decimal.Decimal` canonicalises it (``1.0E0``, ``+1.0``
and ``10E-1`` are one claim), and the value (the correctly rounded double)
and the half unit in the last printed digit are derived from it, never
stored beside it. A stored value with a digit count beside it would be
two fields that can disagree, and "digits" is ambiguous between
significant figures and decimal places.
:meth:`~orpheus.reference.reading.Printed.enclosure` turns the printed
claim into an outward-rounded enclosure,
:math:`[t - h, t + h] \subseteq` the enclosure with :math:`h` the half
unit, so a published value can be compared with an enclosure by one
interval law. The citation must name the place the value is printed (a
table, a page): 8 of the 47 Sood cases print their values in different
places, so a citation is per value, never pooled.

The uncertified value
---------------------

``Uncertified(value)`` is a reference's value where its family cannot
derive a bound on its distance to the exact answer (the user's ruling of
2026-10-03, step 7b: "I think it's a good decision until we face something
that might contradict it (if we face it)"). It has no ``enclosure()``, so
it cannot pass where a guarantee is required. A quotient with an
uncertified operand is uncertified, ``Uncertified(a.value / b.value)``,
and drops the other operand's bound. The ruling chose a third member over
two alternatives: reading such a reference as an enclosure with an
estimated bound (a ladder estimate passing for a guarantee, which the
step-5 ruling forbids) and refusing to read it at all (every existing test
against the trajectory resolvent would have lost its comparison, and the
weaker claim would have left the code instead of being spelled in it).


The observables
===============

An observable names one number to read off the answer of a question. It
holds no mesh and no physics, so the same observable is read by a
production solution on its mesh, by a reference solution in its own
representation, and by a published table. The closed sum (the user's
ruling of 2026-10-03, Q2) is
:data:`~orpheus.numerics.observable.Observable`
``= FluxIntegral | Ratio | Eigenvalue | PointValue``, each member
``@final`` (a subclass would pass every ``isinstance`` door while a
``match`` over the sum sent it to its parent's arm), and each a content
value, so two spellings of one observable are one value and an observable
can key a cache.

The flux integral, the primitive
--------------------------------

:class:`~orpheus.numerics.observable.FluxIntegral` ``(weight)`` is the
linear functional of the scalar flux

.. math::
   :label: reference-flux-integral

   \langle w, \phi\rangle = \sum_g \int_V w_g(\vec r)\,\phi_g(\vec r)\,dV,

.. vv-status: reference-flux-integral documented

.. (vv-status rationale) definitional: the flux integral is the definition of
   an observable, implemented by each answer's read; its readers are pinned by
   foundation rows (R7.8, R7b2.8, R7b2.9), which carry no verifies marker.

with the weight a mesh-free function over position and group (a
:class:`~orpheus.numerics.mesh_free_function.RegionwiseConstant` table or
a :class:`~orpheus.numerics.mesh_free_function.Symbolic` expression). It is
the primitive because every linear observable of the scalar flux is one:
a reaction rate :math:`\langle w, T\phi\rangle` is the flux integral of the
weight :math:`T^{\top} w` (for a removal cell, the weight times its cross
section; for an emission cell, the weight contracted with the emission
spectrum), and a region-averaged group flux is the flux integral of the
region's and the group's indicator divided by the region's volume. The
spelling ``Rate(cells, weight)`` is therefore not a member: it will be a
constructor the **specification** resolves into a flux integral from its
materials, as it resolves a question's keys, so one functional has one
canonical spelling after resolution (step 8, which waits for the removal
cells of posing unit 3, #526). On the infinite medium, which has no
extent, the measure is per unit volume and the flux integral reads
:math:`\sum_g w_g\varphi_g`; its total over the medium is not finite.

**Why a weight cannot depend on the direction.** A flux integral pairs
with the scalar flux, which has no direction. A weight that depends on
:math:`\mu` or on the azimuth has no meaning there, and each reader would
have to decide one (integrate the weight over the sphere of directions?
read it at :math:`\mu = 1`?). The observable refuses such a weight at
construction, before any specification exists, through
:meth:`~orpheus.numerics.mesh_free_function.Symbolic.depends_on`, the one
definition of independence: ``sin(φ)² + cos(φ)²`` spells :math:`\varphi`
and is admitted, because it does not depend on it. A functional of the
angular flux is a different observable, not a flux integral.

The ratio of two linear observables
-----------------------------------

:class:`~orpheus.numerics.observable.Ratio` ``(numerator, denominator)``
has operands of the sum ``Linear = FluxIntegral | PointValue``. An eigen
answer's representative is defined up to a scale :math:`s` (its gauge),
and a reading of a linear observable scales with it,
:math:`J(s\phi) = s\,J(\phi)`. A quotient of two linear observables is
homogeneous of degree 0 in :math:`s`, so it is the same number whatever
gauge the answer was picked at, which is what lets a reference and a
production solution, gauged differently, be compared on it. **A ratio does
not nest**: ``Ratio(Ratio(a, b), c)`` would be of degree :math:`-1`, a
gauge-dependent number dressed as a gauge-free one. Restricting the operands
to the linear members makes every ratio of degree 0 by construction; a
ratio of two ratios would be of degree 0 too, and no consumer needs it. The rule "a ratio reads as the
quotient of its operands' readings" is spelled once,
:meth:`~orpheus.numerics.observable.Ratio.quotient`, which every answer's
``read`` calls (gate R7.19, a route gate: rebinding it to a reversed decoy
moves the ratio readings of every answer at once).

The eigenvalue and the point value
----------------------------------

:class:`~orpheus.numerics.observable.Eigenvalue` ``()`` is a datum of an
eigen question's answer, not a functional of its flux, and it is read in
the chart of the question's parameter (the k chart for the question posed
along the fission emission), so k and 1/k cannot be confused.
:class:`~orpheus.numerics.observable.PointValue` ``(position, group)`` is
the scalar flux of one group at one position on a one-dimensional axis;
its weight would be a Dirac delta, which is no function, so it is its own
member, and an answer with no point evaluation (production's cell
averages) refuses it.

**Admission is decided once.** Whether an observable fits a problem (a
weight's group and region counts, the coordinates a symbolic weight may
depend on, a point's group and position, an eigenvalue asked of an eigen
question) is decided by the specification,
:func:`~orpheus.specification.specification.admit_observable`, never inside
each answer. A published solution calls it for every observable it
prints; a reference solution calls it on every read. A production answer
holds no specification until phase P4's projection, so its ``read`` checks
what it can see (the weight's group count against its flux).

**The observable and its mesh realisation are two objects** (ruling G1).
The mesh-free observable is the value; production realises it on its own
space one level down, as a co-vector paired with its cell values. The
posing plan's earlier rename "Observable → Functional" predates G2 and is
superseded: the mesh-bound
:class:`~orpheus.numerics.functional.Functional` stays as production's
realisation.


How a reference's enclosure is established, and why no ladder certifies
=======================================================================

A reference's enclosure of one observable is established in exactly two
ways, the sum ``Establishment = Exact | DerivedBound``.

**Exact.** :class:`~orpheus.reference.certificate.Exact` ``(expression,
by)`` holds the exact value as a SymPy expression (its ``srepr`` text, a
finite real constant with no free symbol) and names the self-check that
proves it. A rational is exact: its enclosure is
``Enclosure.about(x)``, bound :math:`\lvert\mathrm{fl}(x) - x\rvert`
rounded up. Any other constant (a quadratic irrational qualifies) is
evaluated by SymPy's ``evalf`` in **strict** mode at the working precision
of 60 digits, which raises rather than return a value it cannot certify
to the requested digits; the enclosure's radius is one unit in the
sixtieth digit. A constant whose value SymPy cannot certify (cancellation:
the Machin zero :math:`\arctan\tfrac12 + \arctan\tfrac13 - \pi/4`, exactly 0
without SymPy proving it) raises
:class:`~orpheus.reference.certificate.Uncertifiable`, the one refusal of
``Exact`` that means "no bound can be derived"; it carries the
non-strict value at the working precision, which a derivation may read
as uncertified. Every other refusal (not a constant, not real, not
provably finite) is a defect of the expression and propagates. The
working precision is not a guarantee: at SymPy's default 15 digits a
cancelling sum reads its rounding noise, sign included (`[M]` gate R6.8:
the weight :math:`(z + 10^{-130}, 0)` reads :math:`-1.21\times10^{-122}`
against the true :math:`3.34\times10^{-128}` at 15 digits, and within
:math:`7.9\times10^{-17}` relative at 60), and any fixed precision fails
for a deep enough offset (the strict xfail at :math:`10^{-250}`; #568
certifies such a value with an absolute radius from SymPy's accuracy
record).

**DerivedBound.** :class:`~orpheus.reference.certificate.DerivedBound`
``(enclosure, method)`` is an a posteriori bound with computable constants.
The canonical example for an integral-equation reference is Atkinson's
Nyström residual bound :cite:`Atkinson1997` (Theorem 4.1.2, eq. 4.1.33,
pp. 107–108):

.. math::

   \lVert x - x_n\rVert_\infty \le \lVert(\lambda - \mathcal K_n)^{-1}\rVert\,
   \lVert(\mathcal K - \mathcal K_n)x\rVert_\infty ,

whose second factor is the quadrature error of the integrand
:math:`K(t, \cdot)\,x(\cdot)` (eq. 4.1.35). A series tail with a proved
contraction, interval arithmetic and an integrator's own error control
qualify too. An estimate does not: the tail of a Chebyshev expansion, a
Richardson extrapolant and a residual without its inverse's norm are
estimates.

The user's ruling: no ladder certifies
--------------------------------------

A refinement ladder (a sequence showing the theoretical rate, with or
without a Richardson extrapolant) is not an establishment. The user's
ruling of 2026-10-03: Richardson extrapolation "is a good indicator for
some things, but a bad indicator of correctness. A systematic bias for
example will be invisible to it IF richardson is the standard and not part
of a battery of checks." The consequences, as ruled:

* **No** ``ConvergedLadder`` **establishment and no** ``Extrapolated``
  **reading.** Gate R5c.9 asserts, by an AST census over ``orpheus/``,
  that no class of either name, nor ``RichardsonEstimate`` or
  ``LadderBound``, exists, with ``DerivedBound`` found as its positive
  control.
* **A refinement sequence is falsifying evidence only.** Correct
  enclosures at every resolution all contain the exact answer, so they
  must share a point; a disjoint pair proves one of them wrong. The
  observed orders are reported and decide nothing. ERR-006 is the
  standing reason: a curvilinear sweep converged at exactly
  :math:`O(h^2)` to a wrong limit, which a self-referential gate cannot
  see (:ref:`richardson-extrapolation` retired Richardson as a source of
  reference values for the same reason; this ruling extends it to a
  reference's own refinement).
* **A family that cannot derive its bound has no** ``Valid``
  **certificate**, so it cannot anchor a verification, and deriving each
  family's bound is phase P4's work.

The literature agrees on why an estimate is not a bound. Oberkampf and
Trucano, on series-type analytical references: "One cannot estimate the
accuracy of these type analytical solutions by simply comparing how much
the solution changes by adding one more term" (the 2007 NEA invited
paper, printed p. 83). A Chebyshev tail estimate needs only computed
coefficients and can be fooled by aliasing: :math:`T_{128}` takes the value
1 at every point of the Chebyshev grids with 17, 33 and 65 points
:cite:`AurentzTrefethen2017` (p. 5; `[M]` reproduced on those grids, not on
129, in the literature note). A rigorous bound needs a property of the
function that must be proved: analyticity in a Bernstein ellipse
:math:`E_\rho` with :math:`\lvert f\rvert \le M` there gives
:math:`\lvert a_n\rvert \le 2M\rho^{-n}` and a truncation error at most
:math:`2M/((\rho - 1)\rho^n)` :cite:`Trefethen2008` (Theorem 4.2, eq. 4.7;
Theorem 4.3, eq. 4.11).

.. list-table:: A bound and an estimate on one function, :math:`f(x) = 1/(2 - x)` on :math:`[-1, 1]`, :math:`\rho = 3`, :math:`M = 3` (`[M]` 2026-10-03, mpmath at 40 digits, 401 points)
   :header-rows: 1
   :widths: 16 28 28 28

   * - :math:`n`
     - true sup error of the truncated series
     - tail sum :math:`\sum_{j > n}\lvert a_j\rvert` (an estimate)
     - Bernstein bound (a guarantee)
   * - 5
     - 5.838e-4
     - 5.838e-4
     - 1.235e-2
   * - 10
     - 8.063e-7
     - 8.063e-7
     - 5.081e-5
   * - 20
     - 1.538e-12
     - 1.538e-12
     - 8.604e-10

The tail sum equals the error here because every coefficient is positive
and the function is analytic; nothing in the computed coefficients says
so. The bound is a property of :math:`(\rho, M)`, the continuous object,
and holds whatever the coefficients look like.

.. dropdown:: First got wrong: a ladder as a certifier, one expansion per region, a stored panel bound
   :color: muted

   Three candidates for how a reference states its accuracy were built or specified in P2 and
   withdrawn; each failed for a structural reason that still holds.

   * **The ladder as a certifier** (the first text of the specification's
     R5c.5). It needed two declared inputs, an acceptance band on the
     observed order and a safety factor on the Richardson estimate; neither
     had a derivation (Roache's grid-convergence factor of 1.25 is from a
     paywalled source not read), and the ruling above struck the method.
     The replacement law is band-free by construction: the family of
     enclosures must share a point.
   * **One Chebyshev expansion per region as the answer's representation**
     (the representation prototype, 2026-10-02,
     ``scratch/reference_architecture/p2/representation/report.md``). `[M]`
     Transport scalar fluxes carry a logarithmically singular derivative at
     every scattering free surface and on both sides of every material
     interface, and a :math:`d^{3/2}` term outside a sphere's interface
     (tangent rays); one expansion per region converges algebraically and
     every tail estimate under-predicts the error, by up to
     :math:`5.4\times10^{4}`. The plan's "geometric tail" refusal could not
     tell :math:`k^{-3}` from :math:`\rho^k` over the visible coefficients
     (`[REFUTED 2026-10-02]` as a certifying rule).
   * **Graded piecewise-Chebyshev panels with a stored per-panel bound**
     (step 4, specified and gated, then withdrawn). The panels, graded toward
     the singular set the specification determines, were reliable on six
     fixtures (`[M]` true error over estimate 0.011 to 0.19) and refused the
     panel holding an undeclared singular point (6 of 6). Their bound, "twice
     the sum of the resolved tail", is an estimate under the step-5 ruling, and
     a rigorous one was measured unachievable as built (the reading-bound
     prototype, in the reference-solution section below). The
     specification, drafts and battery are kept in
     ``scratch/reference_architecture/p2/ta/step4/`` as phase P4's starting
     material, where panels return as each reference method's discretisation
     of its emission density.



The reference certificate and its derived state
===============================================

:class:`~orpheus.reference.certificate.ReferenceCertificate`
``(claims, corroborations, refinements, standing)`` holds, per observable,
a :class:`~orpheus.reference.certificate.Claim` ``(target, established)``:
the declared target its bound must meet and its establishment, so a target
without an establishment is unspellable. A ratio is never claimed (it
reads as the quotient of its operands). The target :math:`\tau` is declared
per reference: it cannot be borrowed from a solver bound that does not
exist, and the representation's own error, when P4 stores one, must sit at
:math:`\tau/10` or below.

**The evidence that can refute a claim.**

* A :class:`~orpheus.reference.certificate.Corroboration`
  ``(anchor, independence)``: a current
  :class:`~orpheus.reference.published.PublishedSolution` of the same
  problem, with a note saying why it is independent of the reference. An
  anchor that shares no observable with the claims is refused: a vacuous
  corroboration would read as evidence. A withdrawn publication, an
  enclosure, a float, a bare printed value and a production ``Measured``
  cannot be spelled as anchors; Monte Carlo and a fine-mesh production run
  never are (#505).
* A :class:`~orpheus.reference.certificate.Refinement` ``(observable,
  members)``: the observable's enclosures at a strictly decreasing sequence
  of resolution parameters, falsifying evidence only; its
  ``observed_orders()`` are reported and never decide.

**The state is derived on every read**
(:attr:`~orpheus.reference.certificate.ReferenceCertificate.state`):

1. the :class:`~orpheus.reference.withdrawal.Withdrawal` itself, when the
   standing is one (a maintainer ruling recorded in an issue; the same
   class the test harness and the derivations' run-time lock read,
   :ref:`vv-withdrawn-generators`);
2. otherwise :class:`~orpheus.reference.certificate.Invalid` with every
   failed check's reason: a claim's bound above its target; a refinement
   whose enclosures share no point; a claim that misses the refinement's
   common part; a claim disjoint from an anchor's printed interval; or the
   whole family of one observable sharing no point;
3. otherwise :class:`~orpheus.reference.certificate.Valid`.

The agreement checks are one family law. Every enclosure the certificate
holds for one observable (the claim, each anchor's printed interval, each
refinement member) encloses the one exact value, so together they must
share a point. On a line a family of intervals shares a point exactly when
every pair meets (Helly's theorem in one dimension), so one test over the
family, :func:`~orpheus.numerics.enclosure.common_part`, decides every
pairwise agreement at once. It is decided on the **outward-rounded ends**
(:meth:`~orpheus.numerics.enclosure.Enclosure.ends`), so two enclosures
whose exact intervals meet are never called disjoint; it can miss a
disjointness of about one ulp, which the qa review of step 5 (finding F4)
recorded as the price of never refuting a correct claim. A printed half
unit :math:`5\cdot10^{-k}` is never a double, so an exactly touching
anchor cannot be built.

**The published solution.**
:class:`~orpheus.reference.published.PublishedSolution` ``(specification,
printed, standing)`` maps observables to the
:class:`~orpheus.reference.reading.Printed` values a work prints for a
specification. It reads what it prints and raises
:class:`~orpheus.reference.published.NotPrinted` for anything else; a
withdrawn publication refuses every read, naming its issue. An erratum is
not a withdrawal: it is a cited correction to one value, so it is data, a
``Printed`` citing the erratum in place of the value it corrects. A
published solution is an anchor of a reference certificate and never the
reference a production solution is verified against (the user's ruling of
2026-09-25).


The reference solution, read on demand through its natural extension
====================================================================

:class:`~orpheus.reference.solution.ReferenceSolution`
``(specification, derivation, certificate)`` is a reference method's
answer (ruling G5: a separate type from production's
:class:`~orpheus.sn.solution.Solution`, whose form follows its purpose; the
two share only the reading protocol). All three fields are required, and
``certificate`` is ``None`` for a family with no bound. Its
:class:`~orpheus.reference.solution.Derivation` is a protocol with one
verb, ``evaluate(observable) -> Evaluation``, where ``Evaluation =
Establishment | Uncertified``.

**The user's rulings of 2026-10-03** ("they generally match. You can go
ahead with these rulings"), each binding here:

1. **A reference's reading is its natural extension** (for an
   integral-equation reference, the transport integral of its emission
   density), never a second, nodal answer: one reference, one answer.
2. **P2 builds the architecture without a rigorous reading layer.** Readings
   are evaluated on demand by the derivation with ordinary float quadrature.
   The stored-panel field (step 4) is withdrawn from P2; graded panels
   return in P4 as each reference method's discretisation of its emission
   density, which is what makes a derived solver bound and cheap bounded
   readings possible (#566 for the trajectory resolvent).
3. **Certification is decoupled from caching.** A family without a derived
   bound has no ``Valid`` certificate and cannot anchor a verification
   certificate; phase P3's cache stores its readings regardless, so the
   development cadence improves now (landed:
   :ref:`verification-reference-cache`, where the memo sits below the
   certificate and a cached evaluation is re-checked against every claim
   exactly as a computed one). Existing tests comparing against
   uncertified references keep their current tolerances until each family
   is certified in P4.
4. **No python-flint dependency** until the first family needs rigorous
   arithmetic.

**Construction.** Every claimed observable passes ``admit_observable``; a
claimed ``Ratio`` is refused; every anchor answers the holder's
specification (a certificate holds no specification, so "the same problem"
is enforced by the solution that holds it: without the check a
certificate claiming :math:`k = 1/3` on a slab, anchored by an
infinite-medium publication printing 0.33333, reads ``Valid``, `[M]` the qa
review of step 5, finding F2). Each claimed observable is read once, through the same route
``read`` uses: it must not read ``Uncertified`` (a claim needs a bound), and
the claim's enclosure and the derivation's own must share a point, or the
reference is refused ("disagrees"). A claim and the derivation's
establishment of the same observable are one quantity in two places
(``instrument-doctrine`` X4), so they are checked once, at construction,
and ``read`` is a pure evaluation.

**Reading.** ``read(observable)`` admits the observable, reads a ratio as
the quotient of its operands' readings, and otherwise evaluates it: an
``Exact`` or a ``DerivedBound`` reads as its enclosure, an
``Uncertified`` as itself. Nothing about the answer is stored: an
observable a test poses later (the cell averages of its own mesh) is read
the same way as one the certificate declared. ``read`` answers on an
``Invalid`` or ``Withdrawn`` certificate; only the verification verbs
refuse those. The reference solution is not a content value yet: phase P3
keys it by its derivation's execution trace.

**Why on demand, measured.** The reading-bound prototype (2026-10-03,
``scratch/reference_architecture/p2/reading_bounds/report.md``) built
both readings the user asked for on the trajectory-resolvent sphere of
Garcia 2021, Table 5 :cite:`Garcia2021`. `[M]` An on-demand rigorous
reading (ball-arithmetic Gauss–Legendre along each ray, a Bernstein-ellipse
remainder in :math:`\mu`) read the point value at :math:`r = 3.5` as
:math:`6.7259244498717 \pm 3.0\times10^{-8}` against a true error of
:math:`5.9\times10^{-10}` (bound over true error 51), refused on a node
and on a missing split that a ray crosses, and was **silently wrong** on a
missing split that no ray crosses (:math:`0.11708 \pm 2.8\times10^{-14}`
against 0.12092): the split set is a derived input that must be complete.
It cost 12 minutes per point value; the 40 cell averages of one mesh were
estimated at weeks in Python (`[R]`). Rigorous stored panels were not
achievable as built: the family's per-region cubic spline of the emission
density makes every knot a singular sphere. And the binding error was the
solver's (about :math:`2\times10^{-4}` against Garcia, :math:`10^{4}` times
the reading bound), whose derivable bound (at most 100 times the residual)
was useless. A rigorous reading of an unbounded solver is not an
enclosure end to end, which is why P2 reads with ordinary quadrature and
leaves the bound to P4.


Verification: the certificate, the order, and the two floors
============================================================

:func:`~orpheus.reference.verification.verify_agreement` ``(answer,
observable, reference, tolerance, algebraic_error)`` returns a
:class:`~orpheus.reference.verification.VerificationCertificate`. It holds
the reference ``Valid`` **before reading anything**
(:class:`~orpheus.reference.verification.ReferenceNotValid` names a
missing certificate, an ``Invalid`` state with its first reason, or a
withdrawal with its issue), reads the reference's enclosure (refusing an
uncertified one, :class:`~orpheus.reference.verification.ReadingUncertified`),
and calls ``answer.read(observable)`` itself, so the reading verb is the
only producer of a production reading and a diagnostic that happens to be a
``Measured`` cannot be passed in its place (the elegance review of step 1,
note S4). With :math:`m` production's reading, :math:`v_{\rm ref}` and
:math:`b_{\rm ref}` the reference's enclosure and :math:`\tau` the
tolerance, both decided exactly in rationals:

.. math::
   :label: reference-verification-floor

   b_{\rm ref} \le \tau / 10,

.. math::
   :label: reference-verification-agreement

   \lvert m - v_{\rm ref}\rvert + b_{\rm ref} \le \tau .

The first is the **reference floor**: a reference too loose for the
tolerance verifies nothing, and ``require()`` raises
:class:`~orpheus.reference.verification.ReferenceTooLoose` before it
judges agreement. The second bounds the total error of production's
reading, whatever its split between discretisation and algebra, because
the exact answer lies within :math:`b_{\rm ref}` of :math:`v_{\rm ref}`.
The algebraic error is therefore **recorded** under an agreement verdict
and never required: the verb takes it as an explicit
:data:`~orpheus.numerics.outcome.Evidence` argument with no default, and
every fixture today passes ``NotYet(564, …)``, because no member of the
exit report is a stage-A error estimate of a functional (#564).

:func:`~orpheus.reference.verification.verify_order` returns an
:class:`~orpheus.reference.verification.OrderVerification`: production's
observed order of accuracy over a refinement, against a valid reference.
It attributes error to the discretisation, so it needs the **second
floor**, production's algebraic error small against the tolerance, as a
requirement: every answer's algebraic error must be established
(``Measured`` or ``Asserted``, at most :math:`\tau/10`; ``NotYet`` and
``NotApplicable`` are refused as
:class:`~orpheus.reference.verification.Unestablished`), the reference's
bound must be at most :math:`\tau/10`, and every error must be resolved,
at least :math:`\tau`. Each error :math:`e_i = m_i - v_{\rm ref}` is then
known only within :math:`u_i`, the reference bound plus that answer's
algebraic error, so each observed order is an **interval**:

.. math::
   :label: reference-order-interval

   p_i \in \left[\frac{\log\bigl((\lvert e_i\rvert - u_i) / (\lvert e_{i+1}\rvert + u_{i+1})\bigr)}{\log(h_i / h_{i+1})},\;
                \frac{\log\bigl((\lvert e_i\rvert + u_i) / (\lvert e_{i+1}\rvert - u_{i+1})\bigr)}{\log(h_i / h_{i+1})}\right],

.. vv-status: reference-order-interval documented

.. (vv-status rationale) definitional: the observed-order interval as
   OrderVerification computes it; pinned by the foundation row R7.6.

and the verdict holds only when every interval lies inside the declared
band about the declared order and the errors keep one sign. Point
estimates are not enough: with both floors at their limit, a true order of
1.7 can read 1.963, inside :math:`2 \pm 0.1` (`[M]` the qa review of step
7a). The order verdict is a claim about production against a
reference; it is not a reference certifying itself, so the ruling that no
ladder certifies does not reach it.

Both results are returned and never stored: they are a test's evidence,
not data. **Declared limit:** nothing yet checks that the answer answers
the reference's specification; production results hold no specification
until phase P4's projection of a specification onto a mesh, so the caller
pairs them.

.. vv-status: reference-verification-floor documented
.. vv-status: reference-verification-agreement documented

.. (vv-status rationale) definitional: the floor and the agreement law are
   the definitions of the verification certificate's two verdicts. Their
   gates (R7.3, R7.4 in tests/gates/reference/test_verification.py) are
   foundation rows on the law's arithmetic, which by the harness contract
   carry no verifies marker.


The uncertified reading and the explicit comparison
===================================================

A reference with no derived bound is still a reference: the trajectory
resolvent answers the A|B|A sphere to about :math:`10^{-4}` in k
(self-convergence, not independent; ERR-090), and the tests that compared
with it before P2 should not lose their comparison. The user's ruling of
2026-10-03 keeps the comparison and spells it as the weaker claim it is:
:func:`~orpheus.reference.verification.compare_uncertified`
``(answer, observable, reference, tolerance)``, with no defaults, returns
an :class:`~orpheus.reference.verification.UncertifiedComparison` whose
``agrees`` is :math:`\lvert m - v\rvert \le \tau`, decided exactly, with no
bound on the reference's own error used. ``require()`` raises
:class:`~orpheus.reference.verification.Disagreement` naming the
"unverified reference value".

It compares exactly the readings the verification verbs cannot anchor on,
and refuses the rest, so the two verbs cover every state of a reference in
good standing:

.. list-table:: Which verb a reading admits
   :header-rows: 1
   :widths: 34 33 33

   * - The reference and its reading
     - ``verify_agreement`` / ``verify_order``
     - ``compare_uncertified``
   * - no certificate; any reading (an ``Enclosure`` included, its bound
       unused)
     - refused, ``ReferenceNotValid``
     - compared
   * - ``Valid`` certificate; the observable reads ``Uncertified``
     - refused, ``ReadingUncertified``
     - compared
   * - ``Valid`` certificate; the observable reads an ``Enclosure``
     - verified
     - refused, ``ReadingCertified``, naming ``verify_agreement``
   * - ``Invalid`` or ``Withdrawn`` certificate
     - refused, ``ReferenceNotValid``
     - refused before anything is read

**Why the comparison refuses a certified reading.** If it accepted one, a
test written against an uncertified family would keep making the weaker
claim after the family's certificate landed, silently and forever: the
weaker verb would still be green, and nothing would ask for the
verification the reference can now anchor. Refusing it turns the
certification into a red in every such test, which forces each up to
``verify_agreement`` the day its family is certified. The refusal is
keyed on the certificate, not on the reading's type: an enclosure read
from a reference with no certificate (a family that derives a bound but
declares no targets) is compared, because no verb could otherwise read it.

.. warning::

   **The forbidden sentence**: "the sphere rows verify the S\ :sub:`N`
   solve against the trajectory resolvent". Since #405 P2 step 7b.2.3
   (``a21b6f8e``) the A|B|A sphere rows are ``compare_uncertified``
   (since P1 step (d) of the characteristic reference campaign,
   ``d9425977``, against the characteristic reference, at tolerances
   computed from its ladders), and an uncertified comparison is never a
   verification claim: two methods agreeing to a tolerance, with no bound
   on one of them, is consistency evidence, not correctness evidence
   (``vv-principles`` anti-pattern #1 is the limiting case). The correct
   framing: "the sphere rows compare the S\ :sub:`N` solve with an
   uncertified semi-analytical reference at a tolerance set from both
   methods' measured ladders; they would catch a gross defect (a vacuum
   law realised for the reflective one reddens both rows), and they cannot
   verify either side". The sphere and cylinder rows return as verifications when
   phase P4 gives the reference a derived bound (#566).


An ontology of approximations and errors for verification
=========================================================

The user's draft "G6" (2026-10-02) is the working ontology, to be refined
as the families are certified. It answers what an "approximation" is in a
verification context, where an answer sits, and what each certificate
bounds. The user seeded it with three levels (a truncation at the
continuous level, the transition from continuous to discrete, and the
discrete level: discrete steps and iterative schemes) and ruled that
verification neither cares about nor can assert physical approximations:
validation does. The literature was read for it on 2026-10-02
(``scratch/reference_architecture/p2/error_ontology_literature.md``);
page locators below are the ones checked on the rendered page there.

1. The specification defines exactness
--------------------------------------

Whatever the specification poses is exact for verification. A "physical
approximation" is a choice of specification, and a comparison across that
choice is validation. Larsen and Morel state the instance for energy
:cite:`LarsenMorel2010` (printed p. 6): "if the physical cross sections are
histograms in :math:`E`, then the multigroup transport equations are
exact"; multigroup is exact for a problem whose data are group histograms
(a verification problem) and an approximation only against a
continuous-energy specification, whose error depends on an unknown
weighting flux and has no order of accuracy. Oberkampf and Trucano argue
the same boundary from the other side: "Mixing physics modeling and
numerical solution approximations is, in our view, as bad a[s] mixing
different dimensional units" (the 2007 NEA invited paper, printed
pp. 82–83). A group structure fixed in the problem statement is the model;
a group structure refined toward a continuous-energy limit is a
discretisation; the specification must say which. The governing equation
belongs to the specification in principle (the user's ruling, deferred as
#563 until the methods besides S\ :sub:`N` need it).

2. An approximation is a consistent family
------------------------------------------

An approximation is a family of problems :math:`P_h` with a parameter
:math:`h` whose limit recovers the problem it approximates; its error,
per observable :math:`J`, is :math:`J(u_h) - J(u)`. A method whose limit is
a different problem is not approximating but **defective**: a consistency
defect (a scheme converging to the wrong limit, the first term of
Oberkampf and Trucano's eq. 2, :math:`\lVert u_{\rm exact} -
u_{h,\tau\to0}\rVert`, "primarily the concern of the numerical analyst",
:cite:`OberkampfTrucano2002`, SAND2002-0529 report p. 29), or an
inconsistently discretised acceleration that moves the fixed point:
Adams and Larsen record that such schemes "are not true acceleration
methods; upon iterative convergence, they do not yield the unaccelerated
solution of the spatially discrete S\ :sub:`N` equations"
:cite:`AdamsLarsen2002` (printed p. 70). Verification's object is exactly
this distinction: show that each approximation's error vanishes at its
claimed rate, and to the right limit.

3. The answers live on a chain of problems
------------------------------------------

From the specification to the printed number runs a chain of problems
:math:`P_0 \to P_1 \to \cdots \to P_n`, each link one approximation, and
the total error of an observable telescopes over the links:

.. math::
   :label: reference-error-chain

   J(u_n) - J(u_0) = \sum_{k=1}^{n} \bigl(J(u_k) - J(u_{k-1})\bigr).

.. vv-status: reference-error-chain documented

.. (vv-status rationale) definitional: the telescoping identity that
   defines the chain-of-problems ontology; no function implements it, and
   no test can verify an identity of bookkeeping.

The two-link case is Arioli, Liesen, Międlar and Strakoš's eq. 12
:cite:`ArioliLiesenMiedlarStrakos2013` (preprint p. 6): the total error is
the discretisation error plus the algebraic error,
:math:`u - u_h^{(n)} = (u - u_h) + (u_h - u_h^{(n)})`, and their goal-oriented
form (eqs. 31–35, p. 17) splits a functional's error the same way. Their
energy-norm Pythagoras (eq. 13) needs a symmetric coercive form and
Galerkin orthogonality; a transport S\ :sub:`N` discretisation is neither,
so only the triangle-inequality form applies in ORPHEUS. Oberkampf and
Trucano's validation decomposition (eqs. 13–14, report p. 59) is the same
chain extended upward: :math:`E_1` measurement and :math:`E_2` modelling
belong to validation, :math:`E_3` (the consistency limit) and :math:`E_4`
(finite :math:`h`, iteration and arithmetic) to verification.

A reference that is exact for an intermediate problem with its parameter
held fixed (an S\ :sub:`8`-exact eigenvalue, a flat-source collision
probability) sits at a **node** of the chain. "The equation an answer
declares" is that node. Two answers are compared at a common node, and
verification of production measures the links between the node and
production's answer. P2 admits the continuous node only (the empty set of
fixed approximations); discrete-exact nodes wait for content identity on
the quadrature and the spatial scheme (posing unit 6, #529).

4. Two axes classify a link: the stage and the mechanism
--------------------------------------------------------

The **stages** are the user's three levels: (C) continuous, before any
discretisation; (D) continuous to discrete; (A) the discrete problem
solved inexactly. The **mechanisms** are projection or truncation (onto a
finite subspace or a finite series), quadrature (an integral replaced by a
weighted sum: Nyström, the S\ :sub:`N` angular rule, the multigroup
collapse of a continuous-energy specification; Monte Carlo is quadrature
with random nodes), iteration (finitely many steps of a convergent
sequence) and arithmetic (floating point).

.. list-table:: Stage × mechanism, with an ORPHEUS instance per cell that has one
   :header-rows: 1
   :widths: 16 28 28 28

   * - Mechanism
     - (C) continuous
     - (D) continuous to discrete
     - (A) discrete, solved inexactly
   * - projection / truncation
     - a spectral or series reference truncated at :math:`N` terms
     - a Galerkin method (the two coincide here, below)
     - —
   * - quadrature
     - a reference's own ray or angular integral (the trajectory resolvent's
       reading)
     - the S\ :sub:`N` angular rule; Nyström; the spatial scheme's cell
       integrals
     - —
   * - iteration
     - the Neumann (collision-order) series truncated: a series reference
     - —
     - source iteration, Krylov, the power iteration
   * - arithmetic
     - —
     - —
     - floating point, a floor refinement does not lower

The same mechanism recurs at several stages, and that is not a defect of
the table: the Neumann series truncated after :math:`\ell` collisions is a
stage-C object when it represents the exact solution and a stage-A object
when it is source iteration on the discrete system (the source-iteration
iterate :math:`\ell` is the flux of particles scattered at most
:math:`\ell - 1` times, :cite:`AdamsLarsen2002` eq. 1.10, printed p. 8).
Projection is stage C for a spectral reference and stage D for a Galerkin
method, and for a projection method the two coincide up to stability
constants: :math:`x_n \to x` if and only if :math:`\mathcal P_n x \to x`,
"with exactly the same speed" :cite:`Atkinson1997` (Theorem 3.1.1,
eq. 3.1.35, pp. 55–57). For Nyström they do not: the discretisation error
is the quadrature residual of :math:`K x` propagated by a uniformly bounded
inverse (Theorem 4.1.2), a property of the integrand's integrability by the
rule, not of the solution's approximability.

5. The mechanism decides how a link's error is bounded or estimated
-------------------------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 18 42 40

   * - Mechanism
     - Bound (a guarantee)
     - Estimate, and how it fails
   * - projection / truncation
     - coefficient decay from proved analyticity, :math:`2M/((\rho - 1)\rho^n)`
       :cite:`Trefethen2008` (Thm 4.3)
     - the tail of computed coefficients; fooled by aliasing
       :cite:`AurentzTrefethen2017` and by algebraic decay (the
       representation prototype, above)
   * - quadrature
     - a rule's error term with a bounded derivative or an analyticity
       region; Gauss in :math:`E_\rho`
     - node refinement at the rule's order; unreliable where the integrand
       has an unresolved singularity (ray effects: "the cause ... is not
       numerical but is due to the discrete ordinates formulation itself",
       :cite:`Lathrop1968`, p. 357, so convergence in :math:`N` is
       non-smooth and an asymptotic-order estimator misleads)
   * - iteration
     - a proved contraction :math:`\lVert\mathcal A\rVert < 1`:
       :math:`\lVert x - x^{(k)}\rVert \le \frac{\lVert\mathcal A\rVert}{1 - \lVert\mathcal A\rVert}\lVert x^{(k)} - x^{(k-1)}\rVert`
     - the increment alone: **false convergence**,
       :math:`\lVert\psi - \psi^{(\ell)}\rVert \approx \frac{\sigma}{1 - \sigma}\lVert\psi^{(\ell)} - \psi^{(\ell-1)}\rVert`
       :cite:`AdamsLarsen2002` (eqs. 1.18–1.19, printed p. 10): at
       :math:`\sigma = 0.999` the error is three orders of magnitude above
       the criterion
   * - arithmetic
     - condition number times the unit roundoff
     - none needed; it is a floor refinement does not lower
   * - random quadrature (Monte Carlo)
     - a finite-sample inequality with a bounded score (not built)
     - the central-limit standard error, a random error and a production
       claim only (G3)

The contraction bound in the iteration row is a one-line corollary of
Atkinson's geometric series theorem (Theorem A.1, p. 516:
:math:`\lVert\mathcal A\rVert < 1 \Rightarrow \lVert(I - \mathcal
A)^{-1}\rVert \le 1/(1 - \lVert\mathcal A\rVert)`), derived for this
chapter; Atkinson does not print it. The Adams–Larsen form replaces an
operator norm by the spectral radius :math:`\sigma`, so it is an
asymptotic estimate, not a guaranteed bound, for a non-normal iteration
operator.

6. Inputs are approximated too
------------------------------

Representing the data (a source or a cross-section field projected onto a
mesh, a curved boundary polygonised) is a link of its own, neither a
truncation of the solution nor algebra: the data-oscillation term
:math:`\eta_O` of the estimator in Arioli et al.'s eq. 37 (preprint
p. 19). The specification holds each datum in its most lossless
expression and emits it in the forms a method accepts (ruling G5), so this
link is named where it happens.

7. What is not an approximation
-------------------------------

* The **mode** an eigen question selects and the **gauge** of its
  representative are posed by the question, not errors: an observable
  that depends on the gauge is excluded at the observable level (a ratio
  of linears, above).
* **Programming errors** are what verification exists to find (Oberkampf
  and Trucano's "unacknowledged" errors, report pp. 19 and 28), not a link
  of the chain.
* The **choice of observable** decides the rate (functionals can
  super-converge), so every error and every bound is per observable, as
  the certificates are.

8. Where the certificates sit on the chain
------------------------------------------

* The **reference certificate** bounds the links inside the reference's
  own computation, up to its node.
* Production's self-report, the
  :class:`~orpheus.numerics.outcome.ExitReport`, estimates its stage-A
  links (iteration, arithmetic); it is an input to verification, not a
  certificate.
* The **verification certificate** measures the remaining links, the
  discretisation, against the reference within its bound; the reference
  floor keeps the reference's own links below a tenth of the tolerance,
  and the algebraic floor keeps production's stage-A links apart from the
  discretisation the order verdict attributes error to.

**Where the literature disagrees, declared and not resolved by citation.**
Oberkampf and Trucano's eq. 2 cuts the error at the consistency limit
(:math:`u_{\rm exact} - u_{h\to0}` against everything at finite
:math:`h`, iteration and arithmetic, round-off included); Arioli et al.
cut at the exact discrete solution (:math:`u - u_h` against
:math:`u_h - u_h^{(n)}`). The user's stages D and A follow Arioli. The
AIAA guide that Oberkampf and Trucano quote (report p. 29) puts modelling
simplification and discretisation in one class ("acknowledged error");
the 2007 paper, and this ontology, separate them by the specification.
Lathrop's "not numerical" (1968) means "not spatial, iterative or
round-off"; in verification vocabulary ray effects are an angular
discretisation error.

**Not read** (paywalled, so nothing here rests on them): Oberkampf and
Roy (2010), Roy (2005), Roache (1998), Babuška and Oden (2004), Hoeffding
(1963), and chapters 7–8 of Trefethen's *Approximation Theory and
Approximation Practice*; Trefethen (2008) states the needed theorems. Their
references are in the literature note's ``NEEDS``.


The first instances
===================

The exact infinite medium, the one certified family
---------------------------------------------------

:func:`~orpheus.derivations.common.exact_homogeneous.exact_infinite_medium_reference`
answers one question, the fundamental k-eigenvalue at the physical point
on an infinite-medium specification, through
:class:`~orpheus.derivations.common.exact_homogeneous.ExactInfiniteMediumDerivation`.
The fission operator :math:`F = \chi(\nu\Sigma_f)^{\top}` has rank one, so
:math:`k_\infty = \langle\nu\Sigma_f, A^{-1}\chi\rangle` with eigenvector
:math:`A^{-1}\chi` (:ref:`homogeneous-rank-one-route`), solved in exact
rational arithmetic on the float inputs and certified against the defining
equations (:math:`AA^{-1} = I`, :math:`F\varphi = kA\varphi` with zero
residual). Its flux is read at a production density of 100 per unit
volume, the production being the one the question declares
(``Eigen.gauge``: by default fission and (n,2n) emission, since
2026-10-08; :ref:`structured-geometry-question-values-gauge`). Every
rational reading (k, a flux integral of a tabulated weight) is ``Exact``; a symbolic weight whose value SymPy cannot certify
reads ``Uncertified``. Its certificate claims k and each group's flux at a
target of :math:`10^{-12}` relative (the rationals' only error is their
rounding to a double) and is ``Valid``.

**Honest scope.** Gate R6.6 compares the reference against
``exact_infinite_medium_of``, the function the factory itself calls, so it
is a threading gate (gauge, group order, weighted sum, chart), not a value
gate; the k∞ value rests on ``tests/gates/homogeneous/test_kinf_exact_reference.py``,
certified against the defining equations. The claim-versus-establishment
check is vacuous for this family by construction (both come from one
derivation, X4); gate R6.3 exercises it with a test double.

**The first verification.**
:meth:`~orpheus.homogeneous.solver.HomogeneousResult.read`, the first
production reading (step 7a), reads k∞, flux integrals per unit volume in
the gauge :math:`\langle\nu\Sigma_f, \varphi\rangle = 100`, ratios through
``Ratio.quotient``, and refuses a point value. Gate R7.9 verifies
production against the exact medium with ``verify_agreement`` on three
mixtures at :math:`10^{-13}` relative (`[M]` 2026-10-03, production within
:math:`1.65\times10^{-16}` relative, the specification's measurement; qa
measured a worst :math:`2.56\times10^{-16}` over 5 fissile mixtures), with
the X1 negative: :math:`\nu\Sigma_f` scaled by :math:`1 + 10^{-11}`
disagrees. Declared blindness: the 4-group fixture has equal group fluxes,
so it is blind to a group reversal, which the 2-group rows catch.

The trajectory resolvent, uncertified
-------------------------------------

The multi-region trajectory-resolvent sphere and cylinder were reference
solutions with no certificate from step 7b.2.2 (``311c6148``) until the
family's deletion in P1 step (e2) of #405 (2026-10-10), every reading
``Uncertified``, read as the transport integral of their emission
density. Their reading, its one transport, its splits, its measured
quadrature errors and its costs are recorded with the family's history:
:ref:`characteristic-origins-history-reading`.

The migration of step 7b.2.3 (``a21b6f8e``) re-posed every row that called
the retired test-side helper ``certify_agreement``:

.. list-table:: The A|B|A rows after step 7b.2.3 (`[M]` 2026-10-03, ``python -O -m pytest``)
   :header-rows: 1
   :widths: 30 40 30

   * - Row
     - Now
     - Reading
   * - phase C sphere, k
     - ``compare_uncertified`` at :math:`4\times10^{-3} \times 1.38`
       (the recorded k truncated to three figures)
     - S\ :sub:`N` 1.381079639349644, reference 1.381169542313358
   * - phase C sphere, shape (80 observables)
     - ``compare_uncertified`` per S\ :sub:`N` cell and group, the cell
       average gauged to unit total fission production, at
       :math:`2\times10^{-2} \times 0.249`
     - worst :math:`\lvert m - v\rvert` 1.078e-3 against
       :math:`\tau` = 4.98e-3
   * - phase C cylinder, k and shape; standoff sweep and refinement;
       unified k
     - ``verify_agreement``, strict xfail on ``ReferenceNotValid`` (no
       certificate); an XPASS needs the family certified, not a test edit
     - not compared
   * - phase C cylinder, standoff and unified RECORDs
     - k keys only, read through ``read`` (the cylinder's 80-ratio shape
       set costs about 38 minutes, above the user's 10-minute cap)
     - cylinder k 1.231036749830859, bit-identical to the record

**Declared limit.** A 1 % perturbation
of the sphere reference's moderator emission density at reading time only
moves the shape gap from 1.08e-3 to 1.92e-3, under the row's absolute
tolerance of 4.98e-3: below the row's resolution. Its catcher was the old
reading's fixed-point identity (R7b2.2.2, deleted with the trajectory
resolvent in P1 step (e2) of #405); its successor is the characteristic
reference's fundamental mode satisfying its pencil
(``tests/gates/derivations/test_characteristic_system.py``). Tightening the
row is phase P4's.

No migrated row carries a ``verifies`` marker: an uncertified comparison
is never a verification claim, and a strict xfail verifies nothing. The
label ``sn-curvilinear-characteristic-reference-crosscheck`` is verified by 3
rows (`[M]` the regenerated matrix, 2026-10-03), which verify its
homogeneous reduction only
(:ref:`sn-curvilinear-characteristic-reference-crosscheck-section`).

**Since P1 step (d)** of the characteristic reference campaign
(``d9425977``, 2026-10-10) none of the rows in the table above reads the
trajectory resolvent: the sphere's k and shape rows compare against the
characteristic reference at :math:`p = 5`, at :math:`4\times10^{-5}`
(tightened from :math:`4\times10^{-3}`) and :math:`2\times10^{-2}`; the
cylinder rows' strict xfails refuse the characteristic reference at the
door default (#566); and the RECORD's ``k_ref`` was re-baselined to
1.2317294528902538 while its three S\ :sub:`N` keys did not move. The
table records the rows as step 7b.2.3 left them. Outside the family's own
gates (the Garcia 2021 rows among them), the trajectory resolvent's last
readers are the corroboration rows of step (c), until step (e) deletes
the family. The estimates,
the tolerances and the re-baseline: :ref:`characteristic-sn-rows`.

The production readings
-----------------------

Two production answers read observables, each returning ``Measured``:
:meth:`~orpheus.homogeneous.solver.HomogeneousResult.read` (above) and
:meth:`~orpheus.sn.solution.Solution.read` on the forward
S\ :sub:`N` solution (step 7b.2.1, ``65b0de93``). The S\ :sub:`N` reading
pairs the cell-average flux with the weight's cell integrals,

.. math::

   \langle w, \phi\rangle_h = \sum_{g, i} \phi_{g, i} \int_{V_i} w_g\,dV,

which is exact for a piecewise-constant answer. The co-vector is a mesh
fact, :meth:`~orpheus.mesh.structured.Mesh1D.cell_integrals`: a
regionwise-constant weight is read through the mesh's region labels
(:meth:`~orpheus.mesh.structured.Mesh1D.per_cell`) times the stored cell
volumes; a symbolic weight is integrated exactly in :math:`r` by SymPy (one
antiderivative per group, evaluated at the edges as the binary rationals
they are, then rounded once and scaled by the measure constant), split at
the weight's steps (:meth:`~orpheus.numerics.mesh_free_function.Symbolic.steps`).
The region labels exist for this reading: a region is not a material
(the two outer regions of A|B|A share one material), so a per-region
weight read through the cells' materials would merge regions
(:ref:`the structured-geometry chapter <structured-geometry-mesh>`). Gate
R7b2.9.5 reads a reflective homogeneous slab's k against the exact
infinite medium's reading at :math:`10\times` the eigenvalue tolerance
(`[M]` :math:`\lvert\Delta k\rvert = 6.8\times10^{-14}` at
``keff_tol`` :math:`10^{-12}`): a THEOREM row on raw readings, not
``verify_agreement``, because the slab's specification and the infinite
medium's are different specifications that nothing can pair until P4.
Declared limits: a linear-discontinuous answer's slopes are not read
(#571), and a 2-D mesh refuses every flux integral (#569).


What no gate in P2 can see
==========================

* Whether a production author's printed text was transcribed correctly
  from the paper: that is the review of a ``PublishedSolution``, and the
  per-value citation is what makes it checkable.
* Whether a production answer answers the reference's specification: the
  caller pairs them until P4.
* Whether the singular set a quadrature splits at is complete: a missing
  split no ray crosses is silently wrong (the reading-bound prototype,
  above). The retired trajectory resolvent derived its splits from its
  knots and interfaces, checked by the brute-force rows R7b2.2.1 against
  an unsplit fine rule; those rows were deleted with it in P1 step (e2),
  and the characteristic reference's reading is checked by the flux
  integral of a step weight instead (the sphere supports of
  ``tests/gates/sn/verification/analytical/test_phase_c_crosscheck.py``).


Declared limits and their issues
================================

.. list-table::
   :header-rows: 1
   :widths: 12 58 30

   * - Issue
     - The limit
     - Phase
   * - #566
     - a derived error bound for the characteristic reference (the
       trajectory resolvent's, until its deletion on 2026-10-10); until
       then its readings are uncertified
     - P4
   * - #516
     - the cylinder trajectory-resolvent reference was not converged in
       :math:`n_r` at the cross-check's resolution, its angular rules
       meeting the interfaces' tangency kinks; the family was deleted on
       2026-10-10, and the characteristic reference's line rule is graded
       toward every tangency
     - P4
   * - #563
     - the governing equation belongs to the specification
     - after the methods besides S\ :sub:`N`
   * - #564
     - no production estimator of the algebraic error of a functional; every
       fixture passes ``NotYet``
     - open
   * - #565
     - the admissibility evidence stores the measured k in ``Asserted``'s
       bound slot (a weld)
     - open
   * - #567
     - ``Symbolic`` admits group expressions SymPy proves non-real or cannot
       prove finite
     - open
   * - #568
     - a cancelling constant is read uncertified instead of certified by an
       absolute radius from SymPy's accuracy record
     - P4
   * - #569
     - ``Mesh2D`` stores the composite material map; its cell integrals are
       not built
     - open
   * - #570
     - ``InnerProductFunctional.evaluate`` broadcasts a weight over unmatched
       axes silently
     - open
   * - #571
     - the S\ :sub:`N` reading of a linear-discontinuous answer drops the
       slopes
     - open
   * - #572
     - ``Billiard``'s geometry-kind tag is branched on at five sites and its
       solve result is an untyped metadata entry (moot since P1 step (e2)
       of the characteristic reference campaign deleted ``Billiard``)
     - open
   * - #573
     - an anisotropic reference as isotropic plus an anisotropic correction
       (the user's idea; the trajectory resolvent is verified for what it
       gives, isotropic transport)
     - later
   * - #505
     - Monte Carlo is never a reference
     - standing
   * - #526
     - ``Rate(cells, weight)`` resolved by the specification (step 8) waits
       for the removal cells of posing unit 3
     - with #526
   * - #529
     - discrete-exact nodes (an S\ :sub:`8`-exact eigenvalue) need content
       identity on the quadrature and the scheme
     - posing unit 6

Out of P2 by design: ``Estimated`` (no Monte Carlo consumer), a typed
citation locator (no consumer compares locators), the reference cache
(phase P3, :ref:`verification-reference-cache`) and the certified families (phase P4).


Development history
===================

Reverse-chronological changelog of phase P2 of #405, the plan
``.claude/plans/reference_cache.md`` and its gate specification
``.claude/plans/reference_p2_spec.md``. Every step followed one protocol:
the test-architect's gates and their first reds on a detached worktree,
the main agent's code, a battery of mutations on the real code, and a qa
and elegance review round. Merge status is a git question: the rows dated
2026-10-03 after step 7a sit on the branch
``feature/reference-uncertified-reading`` at the time of writing, and
``git merge-base --is-ancestor <hash> main`` outranks this column.

.. list-table::
   :header-rows: 1
   :widths: 11 55 10 24

   * - When
     - Milestone
     - Issue
     - Where
   * - 2026-10-03
     - **Step 7b.2.3: the migration; the test-side** ``certify_agreement``
       **retires.** The A|B|A sphere rows become ``compare_uncertified``
       at their 2026-09-26 tolerances, the cylinder rows strict xfails on
       ``ReferenceNotValid``, the RECORDs read through ``read``. The sphere's
       shape reading moved from the nodal spline (4.361e-3 relative) to the
       natural extension (4.325e-3) with its tolerance unchanged; the
       cylinder RECORD dropped its two shape keys (the user's 10-minute cost
       ruling); the ``verifies`` marker left the 7 migrated functions, so
       the cross-check label's verifying rows went from 9 to 3 (`[M]` the
       regenerated matrix). ``Billiard`` refuses an anisotropic or (n,2n)
       body material in its payload readers (the user: "Methods should be
       verified for what they give").
     - #405, #566, #516, #573
     - ``a21b6f8e``
   * - 2026-10-03
     - **Step 7b.2.2: the trajectory-resolvent sphere and cylinder as
       reference solutions**, lazy, with no certificate, read as the natural
       extension through one transport (the chord oracle's ``at=`` carve);
       one ``emission_density``, the solves keep their last source.
     - #405, #566
     - ``311c6148``
   * - 2026-10-03
     - **Step 7b.2.1: the S**\ :sub:`N` **production reading**,
       ``Solution.read`` through ``Mesh1D.cell_integrals``.
     - #405, #569, #571
     - ``65b0de93``
   * - 2026-10-03
     - ``Symbolic.without`` **substitutes and never simplifies** (a periodic
       step lost every period after the first; ERR-096).
     - #405
     - ``c861970c``
   * - 2026-10-03
     - **Step 7b.2.0: the mesh carries region labels**; ``mat_ids`` is
       derived (the user's ruling: a mesh holds no geometry, because an
       external mesh has none).
     - #405, #522, #569
     - ``fa38de31``
   * - 2026-10-03
     - **A flux integral refuses a direction-dependent weight; one ratio
       rule.** Before it, the specification's admission admitted 4 of 4
       direction-dependent weight shapes on the sphere and the cylinder
       (`[M]` the test-architect's ``probe4``).
     - #405
     - ``0cb78c84``
   * - 2026-10-03
     - **Step 7b.1: the uncertified reading** (the user's ruling of the
       same day): ``Uncertified``, ``Derivation.evaluate``,
       ``Uncertifiable``, ``compare_uncertified``. The review round keyed
       the comparison's refusal on the certificate (its first version
       refused every enclosure, which left a reference with a derived bound
       and no certificate with no verb).
     - #405, #567, #568
     - ``380c0109``
   * - 2026-10-03
     - **Step 7a: verification.** ``VerificationCertificate`` and
       ``OrderVerification`` (observed orders as intervals, after qa found a
       true order of 1.7 reading 1.963 under point estimates);
       ``HomogeneousResult.read``, the first production reading, verified
       against the exact infinite medium.
     - #405, #564
     - ``b445b68a``, ``9db19fd8``, ``05f27209``; merged at ``17782c7e``
   * - 2026-10-03
     - **Step 6: the reference solution, read on demand; the exact infinite
       medium**, the first certified family.
     - #405
     - ``5af745a8``, ``4b724f04``, ``61c03a32``; merged at ``e06eb95c``
   * - 2026-10-03
     - **The reading-bound rulings**: the reading is the natural extension;
       no rigorous reading layer in P2; certification decoupled from
       caching; step 4 (the stored panel field) withdrawn to P4.
     - #405, #566
     - plan ``99862745``
   * - 2026-10-03
     - **Step 5: published solutions and the reference certificate**, with
       ``Exact`` and ``DerivedBound`` only, after the user's ruling that no
       ladder certifies struck ``ConvergedLadder``; one family agreement law
       (Helly in one dimension); ``Withdrawal`` moves into
       ``orpheus.reference``.
     - #405, #506
     - ``1db08615``, ``f649e6c2``; merged at ``a3ff64d0``
   * - 2026-10-03
     - **Step 3: the observables**, a closed sum; the review round made a
       ratio's operands linear, so ratios do not nest.
     - #405
     - ``9f8d1f7f``, ``cf9893e8``; merged at ``8f6549e3``
   * - 2026-10-03
     - **Steps 1 and 2: the readings by role; the exit record is a report.**
       ``Enclosure`` (outward-rounded quotient), ``Printed`` (its decimal text,
       after the test-architect's question Q1 refuted a stored value beside
       a digit count), the two disjoint sums; ``ExitCertificate`` renamed
       ``ExitReport`` and ``Certified`` renamed ``Asserted``, so that
       "certificate" names two objects where it had named four (`[M]` the P2
       census: the exit record, the orbit certificate, the test-side helper,
       the planned reference certificate).
     - #405
     - ``86748f35``, ``6e04f1f0``, ``ddb5e198``; merged at ``55870ddc``
   * - 2026-10-03
     - **P2 opened**: the census, and the user's rulings G1 to G6 of
       2026-10-02 (the design as ruled).
     - #405, #563
     - plan ``317418da``
