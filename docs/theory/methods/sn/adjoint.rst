.. _sn-adjoint:

Adjoint transport: the dual operators
=====================================

This chapter is the S\ :sub:`N` book's adjoint rung — the dual
operators that the perturbation-theory, detector-sensitivity, and
adjoint-weighted-homogenisation chains consume.  The adjoint chain has
three layers, and **all three are now landed**:

* The **walk adjoint** — the loss composite's transpose
  :math:`(L+C)^{\mathsf T} = L^{\mathsf T} + C` — is machinery of the
  loss *representations* and lives with them: the orientation axis,
  its swap law, and the deferral ledger
  (:ref:`loss-rep-orientation-two-frames`), realised as the Wave-O
  analytic reverse-direction matvec
  (`#280 <https://github.com/deOliveira-R/ORPHEUS/issues/280>`_ /
  `#310 <https://github.com/deOliveira-R/ORPHEUS/issues/310>`_).  This
  chapter points, it does not re-derive.
* The **scattering adjoint** :math:`S^{\mathsf T}`
  (:ref:`sn-scattering-adjoint`, the #276 P3 record): free by frame
  conjugation, with no per-geometry derivation to verify (closes
  `#118 <https://github.com/deOliveira-R/ORPHEUS/issues/118>`_).
* The **daggered posing and the adjoint flux** :math:`\psi^*`
  solving :math:`A_{\rm loss}^{\dagger}\,\psi^* =
  \tfrac1k\,F^{\dagger}\,\psi^*` (with
  :math:`A_{\rm loss} = L+C-S-B`) landed at **#276 A4/A5**
  (:ref:`sn-adjoint-daggered-posing`): the whole eigenproblem is posed
  by DAGGER-ing the forward operator triple through
  :func:`~orpheus.numerics.iteration.KEigenvalue`, and the importance
  map :math:`\varphi^*` rides a role-typed
  :class:`~orpheus.sn.solution.AdjointSolution` carrier
  (:ref:`sn-adjoint-carrier`).  The :math:`\varphi^*` consumers —
  adjoint-weighted homogenisation (frame-machinery P6,
  `#51 <https://github.com/deOliveira-R/ORPHEUS/issues/51>`_ /
  `#281 <https://github.com/deOliveira-R/ORPHEUS/issues/281>`_) and
  perturbation theory / response estimation — are now **unblocked**
  and grow the chapter as they land (:ref:`sn-adjoint-consumers`).

The **spine of the chapter is the route decision**
(:ref:`sn-adjoint-route`): ORPHEUS poses the adjoint by transposing
the *discrete* forward operator, **not** by discretising the
*continuous* adjoint.  Every property below — the exact
:math:`k^{\dagger} = k` identity, the reciprocity that holds at finite
:math:`N` and :math:`h`, the absence of any adjoint-specific loop or
sweep — is a consequence of that one choice.

.. admonition:: Key Facts
   :class: tip

   * **The route (the spine).**  The adjoint is the exact **discrete
     transpose** of the forward operator triple —
     ``KEigenvalue((L+C).H, (S+B).H, F.H)`` (the daggered resolvent,
     gain, and fission; the loss :math:`A_{\rm loss}^{\dagger} =
     (L{+}C).\mathtt{H} - (S{+}B).\mathtt{H}` is formed inside) fed to
     the UNCHANGED
     :func:`~orpheus.numerics.eigenvalue.power_iteration`.  There is
     **no** discretise-then-adjoint step, so duality holds EXACTLY at
     finite :math:`N` and :math:`h` and :math:`k^{\dagger} = k` is an
     exact algebraic identity, not a converged agreement
     (:ref:`sn-adjoint-route`).  The textbook :math:`\mu`-reversal
     (continuous route) survives ONLY as a slab oracle, never in
     production.
   * **Three transposes, one landmine** (:ref:`sn-adjoint-three-transposes`).
     (1) the **Euclidean matrix transpose** :math:`A^{\mathsf T}` (the
     scattering group-transpose :math:`S^{\mathsf T}`, #118, and the
     reverse-scan walk transpose, #280); (2) the **Hilbert / G-metric
     adjoint** :math:`A^{\dagger} = A.\mathtt{H} = G^{-1}A^{\mathsf T}G`
     (:ref:`operator-adjoint`); (3) the **continuous adjoint operator**
     (whose spatial signature is :math:`\mu`-reversal).  Conflating
     them is the #1 way to a plausible-but-wrong adjoint.
   * **The operator algebra IS the implementation.**  No
     adjoint-specific solver code exists anywhere: ``.H`` is the exact
     discrete Hilbert adjoint of every leaf, and the swap law
     :math:`(A^{\dagger})^{-1} = (A^{-1})^{\dagger}`
     (:eq:`loss-rep-adjoint-inverse-swap`) makes :math:`(L+C).\mathtt{H}`
     invertible for free — the daggered inner solve rides the
     reverse-scan transpose sweep behind ``A.H.inverse()``.
   * **The G-metric is a free parameter that no eigenvalue gate can
     see.**  :math:`A.\mathtt{H} = G'^{-1}A^{\mathsf T}G'` is
     metric-similar to :math:`A^{\mathsf T}` for ANY invertible
     :math:`G'`, so :math:`k^{\dagger}` is EXACTLY invariant even under
     a wrong metric.  A metric bug is **k-blind but vector-visible**:
     the coupled-sphere defining-law residual row is its sole catcher
     (:ref:`sn-adjoint-verification`; the ERR-067 family).
   * **The importance carrier is a TYPE, not a flag.**
     :class:`~orpheus.sn.solution.AdjointSolution` is a sibling of
     :class:`~orpheus.sn.solution.Solution` under a role-agnostic base;
     the forward-physics operations (``homogenize`` / ``condense`` /
     ``reaction_rate_density``) are **structurally absent** on it,
     because there is no reaction rate to preserve on an importance map
     (:ref:`sn-adjoint-carrier`).  ``.importance`` is the domain-named
     alias for :math:`\varphi^*`.
   * **The scattering adjoint** :math:`S^{\mathsf T}` is assembled from
     **leaf transposes** of the frame conjugation,
     :math:`(R \circ (\Lambda + N_{2n}) \circ M)^{\mathsf T}
     = M^{\mathsf T} \circ (\Lambda + N_{2n})^{\mathsf T} \circ
     R^{\mathsf T}` — no per-geometry derivation (#276 P3, closes
     #118), reciprocity-pinned against the structurally *independent*
     forward fast-path.  The forward-adjoint **asymmetry is
     principled**: the forward source keeps the scalar fast-path for
     SI-sweep performance; the adjoint — not the hot path — rides the
     validated frame form, which is what makes the reciprocity gate a
     genuine cross-check rather than a tautology.

The continuous adjoint problem and importance
=============================================

Before the discrete machinery, the physics.  The adjoint transport
equation, its importance interpretation, and the reciprocity duality
are the *targets* the discrete construction must reproduce — and the
route decision (:ref:`sn-adjoint-route`) is precisely a claim about how
faithfully it reproduces them.

The adjoint transport equation
------------------------------

The forward within-group transport equation streams a particle along
:math:`\hat\Omega`, removes it at rate :math:`\Sigma_t`, and gains it
back through in-scatter and fission.  The **continuous adjoint** is the
formal transpose of that operator under the phase-space inner product
:math:`\langle a,b\rangle = \int_{\mathcal D}\!\mathrm dV
\int_{4\pi}\!\mathrm d\Omega\,\sum_g a\,b`:

.. math::
   :label: sn-adjoint-continuous

   -\,\hat\Omega\cdot\nabla\psi^*_g
   + \Sigma_{t,g}\,\psi^*_g
   \;=\; \sum_{g'} \Sigma_{s,\,g\to g'}\,\psi^*_{g'}
   \;+\; \frac1k\,\nu\Sigma_{f,g}\sum_{g'}\chi_{g'}\,\psi^*_{g'} .

.. (vv-status rationale) Literature-transcribed definitional identity: the
   continuous adjoint transport equation (Bell & Glasstone §6, Lewis & Miller
   §6).  It states the CONTINUOUS target the discrete construction reproduces;
   it is NOT a per-term solver claim about ORPHEUS code.  Its verifiable
   discrete counterpart is the daggered eigenproblem :eq:`sn-adjoint-eigenproblem`
   and the reciprocity duality :eq:`sn-adjoint-duality`, which the daggered
   posing satisfies EXACTLY (P1.3/P1.2 certification rows).
.. vv-status: sn-adjoint-continuous documented

This is the **eigenvalue** (criticality) adjoint, with the fission gain
scaled by :math:`1/k` and :math:`k^{\dagger} = k` (below).  The
**fixed-source (importance)** adjoint drops the fission term and adds a
prescribed adjoint source — the detector response —
:math:`-\hat\Omega\cdot\nabla\psi^*_g + \Sigma_{t,g}\psi^*_g -
\sum_{g'}\Sigma_{s,\,g\to g'}\psi^*_{g'} = q^*_g = \Sigma_{d,g}`.  Three
sign/role changes distinguish either form from the forward equation, and
each is the continuous face of a discrete transpose the code must
realise:

* **The streaming term flips sign**,
  :math:`\hat\Omega\cdot\nabla \to -\hat\Omega\cdot\nabla` (equivalently
  :math:`\mu\to-\mu`): the adjoint particle propagates *against* the
  physical flow.  This is the **:math:`\mu`-reversal**.
* **The scattering kernel transposes** its energy transfer,
  :math:`\Sigma_{s,\,g'\to g}\to\Sigma_{s,\,g\to g'}`: the adjoint
  particle downscatters where the forward particle upscattered.  With
  the ORPHEUS convention ``SigS[g_from, g_to]``, the forward in-scatter
  source is :math:`\Sigma_s^{\mathsf T}\varphi` and the adjoint source
  is :math:`\Sigma_s\varphi^*` — the transpose is dropped.
* **The fission term swaps emission and production**,
  :math:`\chi_g\,\nu\Sigma_{f,g'} \to \nu\Sigma_{f,g}\,\chi_{g'}`: the
  adjoint particle is "born" with the production spectrum
  :math:`\nu\Sigma_f` and "weighted" by the emission spectrum
  :math:`\chi`.  This :math:`\chi\leftrightarrow\nu\Sigma_f` swap is the
  canonical adjoint-fission trap.

Importance — the interpretation of the adjoint flux
---------------------------------------------------

The adjoint flux is not a flux.  :math:`\psi^*(\vec r,\hat\Omega,g)` is
the **importance**: the expected contribution to a detector response
:math:`\langle\Sigma_d,\psi\rangle` from a single particle introduced
at the phase-space point :math:`(\vec r,\hat\Omega,g)`
(:cite:`BellGlasstone1970` §6; :cite:`LewisMiller1984` §6;
:cite:`Lux1991` for the Monte-Carlo importance).  A neutron born deep
in the fuel with an inward-pointing direction is *more important* to
the chain reaction than one born at the vacuum boundary heading out —
importance is direction- and energy-resolved, and on a finite system
it is genuinely :math:`\mu`-asymmetric, which is exactly the content
the :math:`\mu`-reversal carries.

For the eigenvalue problem the adjoint source vanishes
(:math:`q^* = 0`) and the importance is the fundamental **adjoint
eigenmode** :math:`\varphi^*` — the reactor's "worth function": the
first-order sensitivity of :math:`k` to a perturbation introduced at
each phase-space point.  This is why perturbation theory and
generalised perturbation theory (:ref:`sn-adjoint-consumers`) are the
adjoint flux's native consumers.

Reciprocity and the fundamental duality
----------------------------------------

The single identity from which every adjoint application descends: for
a forward fixed-source solve :math:`A_{\rm loss}\,\psi = q` and an
adjoint solve :math:`A_{\rm loss}^{\dagger}\,\psi^* = \Sigma_d` driven
by a detector response :math:`\Sigma_d`,

.. math::
   :label: sn-adjoint-duality

   \langle \Sigma_d,\,\psi\rangle
   \;=\;
   \langle \psi^*,\, q\rangle ,

.. The label sn-adjoint-duality is a verifies-target (a solver claim:
.. the discrete reciprocity the daggered posing satisfies EXACTLY).
.. It carries NO vv-status sentinel by the wired⟹no-sentinel convention —
.. the pinning L1 gate (test_sn_adjoint_entries.TestSolveSnAdjointFixedSource.
.. test_duality_cross_group_source_detector) carries the @verifies marker
.. (wired at the A6 marker commit).

where the discrete pairing is the solution G-inner-product
:math:`\langle\cdot,\cdot\rangle_G` (the discrete realisation of the
phase-space integral): the detector reading computed from the forward
flux equals the source
weighted by the importance.  The proof is one line —
:math:`\langle\psi^*, q\rangle = \langle\psi^*, A_{\rm loss}\psi\rangle
= \langle A_{\rm loss}^{\dagger}\psi^*, \psi\rangle = \langle\Sigma_d,
\psi\rangle` — using only the defining adjoint relation
:math:`\langle\psi^*, A\psi\rangle = \langle A^{\dagger}\psi^*,
\psi\rangle`.  It is exact for the *continuous* operators, and — this
is the whole point of the route decision — it is exact for the
*discrete* operators too, at finite :math:`N` and :math:`h`, because
:math:`A_{\rm loss}^{\dagger}` is built as the exact transpose of the
discrete :math:`A_{\rm loss}` (:ref:`sn-adjoint-route`).  The
verification rows (P1.2, :ref:`sn-adjoint-verification`) exercise the
un-hideable case — source and detector in *different* energy groups
*and* regions — so the reciprocity forces the downscatter chain
through :math:`S^{\mathsf T}` and a dropped transpose reddens O(1).

The adjoint eigenproblem and the exact eigenvalue identity
----------------------------------------------------------

The forward criticality problem is :math:`A_{\rm loss}\,\psi =
\tfrac1k\,F\,\psi` with :math:`M = A_{\rm loss}^{-1}F` and
:math:`k = \lambda_{\max}(M)`.  The adjoint criticality problem is its
transpose,

.. math::

   A_{\rm loss}^{\dagger}\,\psi^* \;=\; \frac{1}{k^{\dagger}}\,
   F^{\dagger}\,\psi^* ,
   \qquad
   M^{\dagger} = \bigl(A_{\rm loss}^{\dagger}\bigr)^{-1}F^{\dagger},
   \qquad
   k^{\dagger} = \lambda_{\max}(M^{\dagger}) .

**The eigenvalue is unchanged:** :math:`k^{\dagger} = k`, an *exact
algebraic identity*.  With the Hilbert adjoint
:math:`A^{\dagger} = G^{-1}A^{\mathsf T}G` (and likewise
:math:`F^{\dagger} = G^{-1}F^{\mathsf T}G`),

.. math::

   M^{\dagger}
   = \bigl(G^{-1}A^{\mathsf T}G\bigr)^{-1}\!\bigl(G^{-1}F^{\mathsf T}G\bigr)
   = G^{-1}\,A^{-\mathsf T}F^{\mathsf T}\,G
   = G^{-1}\,\bigl(F A^{-1}\bigr)^{\mathsf T}\,G ,

so :math:`M^{\dagger}` is similar to :math:`(FA^{-1})^{\mathsf T}`,
whose spectrum equals that of :math:`FA^{-1}`, which equals that of
:math:`A^{-1}F = M` (the shared non-zero spectrum of :math:`XY` and
:math:`YX`).  Hence :math:`\operatorname{eig}(M^{\dagger}) =
\operatorname{eig}(M)` and :math:`k^{\dagger} = k` *for any invertible
metric* :math:`G`.

The **eigenVECTOR is not** :math:`\varphi`, however.  :math:`\varphi^*`
is the dominant *left* eigenvector of :math:`M`, equivalently the
dominant right eigenvector of :math:`(A^{\mathsf T})^{-1}F^{\mathsf T}`.
That the value is shared but the vector is not is the load-bearing
Mode-12 fact of this whole chapter: **a** :math:`k^{\dagger} = k`
**gate carries zero information about** :math:`\varphi^*`, so every
eigenvalue check is *designed-green* on the entire adjoint mutation
class.  The vector must be pinned separately, and by a reference of the
right structural form — the :math:`(A^{\mathsf T})^{-1}F^{\mathsf T}`
spectrum, **not** :math:`\operatorname{eig}(M^{\mathsf T})`
(:ref:`sn-adjoint-verification` records the factor-order trap in full).

.. _sn-adjoint-route:

The route decision — the exact discrete transpose
=================================================

This is the spine of the chapter.  Every property above — the exact
:math:`k^{\dagger} = k`, the reciprocity that holds at finite mesh, the
absence of any adjoint-specific solver code — is a consequence of a
single design choice about *when* to transpose.

Two routes to a discrete adjoint
--------------------------------

There are two ways to obtain a discrete adjoint solver, and they do
**not** commute:

#. **Discretise-then-adjoint** (the textbook continuous route).  Take
   the continuous adjoint equation :eq:`sn-adjoint-continuous` —
   :math:`\mu`-reversal, kernel transpose, :math:`\chi\leftrightarrow
   \nu\Sigma_f` swap, all written in continuous form — and discretise
   it from scratch with its own sweep, its own upwinding, its own
   angular closure.
#. **Adjoint-of-the-discrete** (the ORPHEUS route).  Take the *already
   discretised* forward operator :math:`A_{\rm loss}` and form its
   exact matrix transpose, wrapped in the phase-space metric:
   :math:`A_{\rm loss}^{\dagger} = G^{-1}A_{\rm loss}^{\mathsf T}G`.

The two agree only in the limit :math:`h\to 0`, :math:`N\to\infty`.  At
finite resolution they differ by a discretisation error, and route (1)
carries a subtle, expensive failure: the discretised-continuous-adjoint
operator is **not** the transpose of the discrete forward operator, so
the discrete reciprocity :eq:`sn-adjoint-duality` holds only to
:math:`\mathcal O(h^p)`, and :math:`k^{\dagger} \ne k` except in the
limit.  A perturbation-theory or GPT chain built on route (1) inherits
that inconsistency as a spurious first-order error term — the adjoint
is "correct" only asymptotically, exactly where it is least useful for
a sensitivity estimate on a real, coarse mesh.

Why ORPHEUS transposes the discrete operator
--------------------------------------------

ORPHEUS takes route (2): it poses the adjoint eigenproblem by
**DAGGER-ing the forward operator triple**, feeding
:func:`~orpheus.numerics.iteration.KEigenvalue` the daggered triple
``((L+C).H, (S+B).H, F.H)`` — the daggered RESOLVENT, gain, and
fission; the loss dagger is formed inside the posing — and running
the **unchanged**
:func:`~orpheus.numerics.eigenvalue.power_iteration`
(:func:`~orpheus.sn.solver.solve_sn_adjoint`).  The consequences are
exactly the properties route (1) cannot deliver at finite resolution:

* **Reciprocity holds EXACTLY at finite** :math:`N` **and** :math:`h`.
  :math:`A_{\rm loss}^{\dagger}` *is* the discrete transpose, so
  :math:`\langle\Sigma_d,\psi\rangle = \langle\psi^*,q\rangle` is an
  algebraic identity of the discrete system, not an
  :math:`\mathcal O(h^p)` approximation.
* :math:`k^{\dagger} = k` **is an exact algebraic identity**
  (derived above), not a converged agreement of two independent
  discretisations.
* **There is no adjoint-specific loop or sweep code anywhere.**  The
  operator algebra IS the implementation: ``.H`` is the exact discrete
  Hilbert adjoint of every leaf, composed through
  :meth:`OperatorSum.apply_transpose
  <orpheus.numerics.operator.OperatorSum>` /
  :meth:`OperatorProduct.apply_transpose
  <orpheus.numerics.operator.OperatorProduct.apply_transpose>`; the
  **swap law** :math:`(A^{\dagger})^{-1} = (A^{-1})^{\dagger}`
  (:eq:`loss-rep-adjoint-inverse-swap`) makes :math:`(L+C).\mathtt{H}`
  invertible by routing to :math:`(L+C).\mathtt{inverse()}.\mathtt{H}`,
  so the daggered inner solve rides the reverse-scan transpose sweep
  (:ref:`loss-rep-orientation-two-frames`) behind ``(L+C).H.inverse()``
  with no new machinery.

This is Cardinal Rule 2 (architecture) paying a physics dividend: the
adjoint is not a parallel implementation to keep in sync with the
forward one — it is the *same* operators, transposed.  A bug in the
forward streaming operator and its adjoint cannot silently diverge,
because there is only one operator.

The continuous route, kept only as a slab oracle
------------------------------------------------

Route (1) is not discarded — it is kept as a **verification oracle**.
On the 1-D slab, where the geometry is simple enough that the
:math:`\mu`-reversed continuous adjoint discretises cleanly, the
continuous route provides a *structurally independent* reference for
the discrete-transpose adjoint (the fuller-view-oracle exception of the
project's retirement discipline).  It is never a production path: a
user never solves the adjoint by flipping :math:`\mu` and re-sweeping.
The distinction matters because :math:`\mu`-reversal *looks* like the
adjoint and is the single most common way to build a
plausible-but-wrong one — which is the subject of the next section.

.. _sn-adjoint-three-transposes:

The three transposes — the recurring landmine
=============================================

The word "transpose" names three different objects in the adjoint
chain, and conflating any two of them produces a plausible-but-wrong
adjoint that passes every eigenvalue gate.  This is the **#1 landmine**
of the whole carve.  The symptom table (:ref:`sn-symptom-table`) routes
"the adjoint reciprocity gate reds" straight here.

(1) The Euclidean matrix transpose
----------------------------------

:math:`A^{\mathsf T}` — the plain linear-algebra transpose of the
discrete operator's matrix, carrying **no metric**.  It appears in two
places in ORPHEUS, both realised without ever forming the matrix:

* the **scattering group-transpose** :math:`S^{\mathsf T}`
  (:ref:`sn-scattering-adjoint`), deliberately Euclidean — the frame
  conjugation transposes leaf-by-leaf, closes #118;
* the **reverse-scan walk transpose** :math:`(L+C)^{\mathsf T}` — the
  loss representation's orientation axis
  (:ref:`loss-rep-orientation-two-frames`): reversed DAG order + face
  in↔out swap over the *same* per-ordinate cell graph.

The walk transpose is the object the loss-representation warning
contrasts against :math:`\mu`-reversal: reversing the *within-octant
cell order* is what gives :math:`(L+C)^{\mathsf T}`, **not** reflecting
the angular quadrature.  Both :math:`S^{\mathsf T}` and
:math:`(L+C)^{\mathsf T}` are category (1).

(2) The Hilbert (G-metric) adjoint
----------------------------------

:math:`A^{\dagger} = A.\mathtt{H} = G^{-1}A^{\mathsf T}G` — the
metric-weighted adjoint, the *physical* one, and the object the
daggered posing actually uses.  It rides *on top of* the Euclidean
transpose: the function space owns the metric :math:`G` and the ``.H``
wrapper (:class:`~orpheus.numerics.operator.LinearOperator`,
:ref:`operator-adjoint`) applies :math:`G` before and :math:`G^{-1}`
after the leaf's ``apply_transpose``.  The metric is the phase-space
measure (:eq:`g-adjoint-block-metric`):

.. math::

   G \;=\; \operatorname{diag}\bigl(G_{\rm bulk},\,G_{\rm trace},\,
   G_{\rm sd}\bigr),
   \qquad
   G_{\rm bulk} = V_{\rm cell}\,w_n,
   \quad
   G_{\rm trace} = |\Omega\cdot\hat n_f|\,w_n,
   \quad
   G_{\rm sd} = V_{\rm cell}
   \qquad(\text{sd = the starting-direction ray — the spelling the}
   \text{tests and the error catalog use}),

the bulk volume·weight block :math:`V_{\rm cell}\,w_n`, the
partial-current trace block :math:`|\Omega\cdot\hat n_f|\,w_n` (with a
pseudo-inverse on the singular grazing-ordinate trace), and — on a
carrying (sphere) mesh — the System-B ray block :math:`G_{\rm sd} =
V_{\rm cell}`.  The sweep itself carries **no metric code**: the metric
enters only at the space boundary, so the same ``.H`` wrapper serves a
flat spherical-harmonic metric and a composite ``FullField`` metric
alike.

(3) The continuous adjoint operator
-----------------------------------

The adjoint of the *continuous* transport operator
(:eq:`sn-adjoint-continuous`), whose discrete signature is
**:math:`\mu`-reversal** (reflecting the angular quadrature) plus the
continuous kernel transpose.  This is what the discretise-then-adjoint
route (1) would build.  ORPHEUS does **not** use it in production — it
survives only as the slab oracle.  The loss-representation warning's
"NOT :math:`\mu`-reversal and NOT the continuous transport adjoint" is
precisely the statement that the walk transpose (category 1) must not
be confused with category (3).

The taxonomy reconciles the two framings a reader will meet elsewhere:
the loss-representation carve's "three transposes" are {the walk's
Euclidean transpose, :math:`\mu`-reversal (which it identifies with
the continuous adjoint — one object in its framing), and the Hilbert
G-adjoint riding on top}, and the thin pre-A6 Key Facts named
the trio *Euclidean / Hilbert / walk-orientation* — the walk-orientation
transpose being simply the streaming realisation of category (1).  All
three framings are the same taxonomy: **Euclidean** (bare
:math:`A^{\mathsf T}`, including both :math:`S^{\mathsf T}` and the
walk), **Hilbert** (:math:`G^{-1}A^{\mathsf T}G`, riding on top), and
**continuous** (whose signature is :math:`\mu`-reversal).

The G-metric free-parameter lesson
----------------------------------

The most dangerous property of the Hilbert adjoint, and the one every
future adjoint session must internalise: **the metric** :math:`G` **is a
free parameter that no eigenvalue gate can ever see.**  Because
:math:`A.\mathtt{H} = G'^{-1}A^{\mathsf T}G'` is *metric-similar* to
:math:`A^{\mathsf T}` for **any** invertible :math:`G'`, the daggered
spectrum — and therefore :math:`k^{\dagger}` — is EXACTLY invariant
even under a **wrong** metric.  A metric bug is **k-blind but
vector-visible**.

This is not a hypothetical.  The ghost :math:`G_{\rm sd} \equiv 0`
defect (ERR-067) put the System-B seed rows in :math:`\ker G` and made
``A.H`` a *wrong* adjoint for any non-zero seed, while every eigenvalue
gate stayed green.  The sole catcher is the coupled-sphere
defining-law residual row (:ref:`sn-adjoint-verification`): dropping
:math:`G_{\rm sd} = V_{\rm cell} \to 1` leaves
:math:`|k^{\dagger}_{\rm mut} - k_{\rm fwd}| = 2.6\times10^{-11}`
(EXACTLY k-blind) while the vector residual reds O(1) at
:math:`2.35`.  **NEVER** certify a metric by an eigenvalue — pin the
adjoint *vector* against a metric built independently from raw mesh
data.

.. _sn-scattering-adjoint:

The scattering adjoint, free from the harmonic frame
====================================================

The loss composite's analytic adjoint is hard — sign-flipping the upwind
direction, transposing the M–M closure, re-deriving the per-level azimuthal
redistribution, each an AI-failure-mode trap.  It is carried as the
**orientation axis** of the loss-representation machinery
(:ref:`loss-rep-orientation-two-frames`): the Wave-O analytic
reverse-direction matvec, its swap law, and the deferral ledger — landed
after a dense-transpose interim.  The **scattering**
operator :math:`S` is the counterexample: campaign **#276 P3** (commit
``15185e5``, closes
`#118 <https://github.com/deOliveira-R/ORPHEUS/issues/118>`_) made
:math:`S^{T}` fall out **for free**, because :math:`S` is already written as
a harmonic-frame conjugation.

The modernised in-scatter source is ONE frame-conjugated operator
(:attr:`~orpheus.transport.operators.scattering.ScatteringOperator.full_scatter_kernel`):

.. math::
   :label: sn-scattering-adjoint-kernel

   \mathrm{full\_scatter\_kernel}
   \;=\; R \circ (\Lambda_{\ell\ge 0} + N_{2n}) \circ M ,

.. (vv-status rationale) Representational identity: the frame-conjugation
   definition of the full P0+ℓ≥1+(n,2n) in-scatter kernel (analysis M,
   moment-space transfer Λ+N₂ₙ, reconstruction R).  Its verifiable content —
   the frame form reproduces the independent scalar fast-path forward source —
   is the ``@pytest.mark.foundation`` gate
   ``tests/sn/operators/test_scattering_adjoint.py::TestFullScatterKernel::test_reproduces_forward_scattering_source``
   (rtol 1e-12); the gate is unwired, so the label stays ``documented``
   with the gate named here (wiring backlog: #309).
.. vv-status: sn-scattering-adjoint-kernel documented

where :math:`M` / :math:`R` are the angular frame's analysis /
reconstruction faces, :math:`\Lambda_{\ell\ge 0}` is the per-:math:`\ell`
moment-space group transfer
(:class:`~orpheus.transport.operators.scattering.LegendreMomentScattering`),
and :math:`N_{2n}` is the distinct :math:`(n,2n)` multiplication channel
(:class:`~orpheus.transport.operators.scattering.N2NMomentOperator`) —
summed with :math:`\Lambda` in moment space and conjugated by the frame
*together* (one analysis, one reconstruction) for the WHOLE
P0 + :math:`\ell\ge1` + :math:`(n,2n)` source.  Its transpose is therefore
the product transpose

.. math::
   :label: sn-scattering-adjoint-kernel-transpose

   \mathrm{full\_scatter\_kernel}^{T}
   \;=\; M^{T} \circ (\Lambda + N_{2n})^{T} \circ R^{T},

.. (vv-status rationale) Structural / representational identity: the product
   transpose assembled from the leaf transposes (no per-geometry derivation).
   Its verifiable content is the Euclidean reciprocity ⟨kernel ψ, c⟩ =
   ⟨ψ, kernelᵀ c⟩, pinned by the ``@pytest.mark.foundation`` gate
   ``tests/sn/operators/test_scattering_adjoint.py::TestFullScatterKernel::test_full_kernel_euclidean_reciprocity``
   (scalar + LD trailing spectator) — foundation gates carry no
   ``verifies(...)`` by design.
.. vv-status: sn-scattering-adjoint-kernel-transpose documented

which :meth:`OperatorProduct.apply_transpose
<orpheus.numerics.operator.OperatorProduct.apply_transpose>` assembles from
the leaf transposes — the frame's :math:`M^{T}` / :math:`R^{T}` faces (landed
in the Frame/Basis carve), the per-:math:`\ell` group transpose
:math:`\Lambda^{T}`, and :math:`N_{2n}^{T}` — with **no per-geometry
derivation to verify** (the trap the streaming adjoint above could not
avoid).  The :term:`per-ordinate <ordinate>` adjoint scattering source is then

.. math::
   :label: sn-scattering-adjoint-source

   S^{T}\chi \;=\; \tfrac{1}{W}\,\mathrm{full\_scatter\_kernel}^{T}\,\chi ,

.. (vv-status rationale) Definitional identity: the per-ordinate adjoint
   scattering source (the producer-side 1/W transposing as the scalar it is).
   Its verifiable content is the LOAD-BEARING per-group Euclidean reciprocity
   ⟨Sψ,χ⟩=⟨ψ,Sᵀχ⟩ — the frame-form Sᵀ cross-checked against the structurally
   INDEPENDENT scalar fast-path S — plus the S.apply_transpose == (1/W)·kernelᵀ
   wiring gate, both ``@pytest.mark.foundation`` in
   ``tests/sn/operators/test_scattering_adjoint.py::TestFullScatterKernel``;
   both gates are unwired, so the label stays ``documented`` with the
   gates named here (wiring backlog: #309).
.. vv-status: sn-scattering-adjoint-source documented

the producer-side :math:`1/W` transposing as the scalar it is
(:math:`(A/W)^{T} = A^{T}/W`).
:class:`~orpheus.transport.operators.scattering.ScatteringOperator` now
reports ``is_adjointable = True`` (it has a working ``apply_transpose``),
and the old "no ``apply_transpose``" class-docstring confession is
retired.

**Forward fast-path, adjoint frame-path — and why the asymmetry is
principled.**  The production FORWARD source keeps the scalar fast-path
(:attr:`~orpheus.transport.operators.scattering.ScatteringOperator.isotropic_kernel`
for P0 + :math:`(n,2n)`, and the per-:math:`\ell` ``build_aniso_source``)
for SI-sweep performance; the **adjoint** — not the hot path — rides the
validated frame form instead.  The two are thus structurally *different*
representations of the same operator, which is exactly what makes the
verification a genuine cross-check rather than a tautology: the per-group
Euclidean reciprocity
:math:`\langle S\psi, \chi\rangle = \langle\psi, S^{T}\chi\rangle`
(``tests/sn/operators/test_scattering_adjoint.py``,
``TestFullScatterKernel::test_S_euclidean_reciprocity``) pins the frame-form
:math:`S^{T}` against the *independent* scalar fast-path :math:`S`, and the
forward equivalence
:math:`(1/W)\,\mathrm{full\_scatter\_kernel}.\mathrm{apply} \equiv
S.\mathrm{apply}` holds to :math:`\sim 10^{-12}`.

.. note::

   This :math:`S^{T}` is the **Euclidean** transpose (the plain
   group-and-angle matvec adjoint), NOT the metric Hilbert adjoint
   :math:`S^{\dagger} = G^{-1}S^{T}G` — that angular-Gram weighting is the
   :attr:`~orpheus.numerics.operator.LinearOperator.H` wrapper's job.  The
   campaign and commit name it "S†" colloquially; the precise object the
   operator computes is the transpose.

This is the discrete scattering adjoint the SN adjoint chain builds on: the
adjoint flux :math:`\psi^{*}` solving :math:`(L+C-S)^{T}\psi^{*} = q^{*}`,
adjoint-weighted homogenisation, perturbation theory, and detector
sensitivity all need :math:`S^{T}`.  Its companion forward step (campaign
**#276 P2**, commit ``dcea43a``) routes the SN forward *isotropic* source
through the same model-shared
:class:`~orpheus.transport.operators.isotropic_scattering.IsotropicScattering`
(:math:`\Sigma_{s0}`) and
:class:`~orpheus.transport.operators.isotropic_scattering.IsotropicN2N`
(:math:`2\Sigma_{2n}`) operators (0-ULP bit-identical), so the
:math:`K_\mathrm{iso}` energy operator — which also assembles the
infinite-medium loss matrix (:ref:`direct-eigensolve-assembly`) — is one
cross-model source.  These model-shared operators live in
:mod:`orpheus.transport.operators`.

.. _sn-adjoint-daggered-posing:

The daggered posing
===================

The adjoint entries are :func:`~orpheus.sn.solver.solve_sn_adjoint`
(eigenvalue) and
:func:`~orpheus.sn.solver.solve_sn_adjoint_fixed_source`
(importance / detector), module-level siblings of the forward family.
Neither spells adjoint physics: both hand
:func:`~orpheus.numerics.iteration.KEigenvalue` /
:class:`~orpheus.numerics.iteration.SourceIteration` a **daggered**
operator triple and let the operator algebra do the rest.

The daggered eigenproblem
-------------------------

The shared build ``_adjoint_posing_parts`` takes the forward
within-group system off the single construction site
:func:`~orpheus.sn.coupled_system.build_within_group_system` and returns
its daggerable parts — the invertible resolvent :math:`(L+C)`, the
summed coupling gain :math:`(S+B)`, and the fission operator :math:`F`.
The CALLER daggers each with ``.H`` and poses

.. math::
   :label: sn-adjoint-eigenproblem

   A_{\rm loss}^{\dagger}\,\psi^* \;=\; \frac1k\,F^{\dagger}\,\psi^*
   \qquad\Longleftrightarrow\qquad
   \mathtt{KEigenvalue}\bigl((L{+}C).\mathtt{H},\;(S{+}B).\mathtt{H},\;
   F.\mathtt{H}\bigr),

with :math:`A_{\rm loss}^{\dagger} = (L+C).\mathtt{H} - (S+B).\mathtt{H}`
fed to the **unchanged** :func:`~orpheus.numerics.eigenvalue.power_iteration`
(the adjoint row of the eigenvalue-posing table,
:mod:`orpheus.numerics.eigenvalue`).  The within-group loss splits and
each term daggers independently:

* :math:`(L+C).\mathtt{H}` is invertible **for free** by the swap law
  (:eq:`loss-rep-adjoint-inverse-swap`): ``(L+C).H.inverse()`` routes to
  ``(L+C).inverse().H``, i.e. the reverse-scan transpose sweep
  (:ref:`loss-rep-orientation-two-frames`), so the daggered inner solve
  is the reversed walk with no new solver code.
* :math:`S.\mathtt{H}` is the frame-conjugated scattering transpose
  (:ref:`sn-scattering-adjoint`) wrapped in the angular metric.
* :math:`B.\mathtt{H}` is the boundary-law transpose; the reflective
  and vacuum traces transpose structurally — an adjoint vacuum is the
  transpose of the forward vacuum, never a user-facing BC flip.
* :math:`F.\mathtt{H}` is the metric-wrapped fission transpose whose
  Euclidean core is the :math:`\chi\leftrightarrow\nu\Sigma_f` role swap
  :math:`F^{\mathsf T}\psi^* = \nu\Sigma_f\,(\chi\cdot\psi^*)`
  (:class:`~orpheus.transport.operators.fission.FissionOperator`).

Because :math:`k^{\dagger} = k` is exact, the entry returns
``keff`` equal to the forward eigenvalue to iteration tolerance, and
its ``angular_flux`` is the true discrete adjoint (importance) flux
:math:`\psi^*` — verified against the closed-form
:math:`(A^{\mathsf T})^{-1}F^{\mathsf T}` spectrum
(:ref:`sn-adjoint-verification`).

The coupled (sphere) posing
---------------------------

On a carrying mesh — the sphere, whose half-angle starting-direction
seed is first-class System-B state — the posing is a 2×2 block operator
over System A (the transport bulk ⊕ trace) and System B (the
radial-characteristic ray).  The gain is the builder's own coupled gain
grid :math:`N`; the fission operator is lifted to the coupled grid:

.. math::

   F_{\rm posed} \;=\;
   \begin{pmatrix}
     F & 0 \\[2pt]
     A_{BA}^{\rm fis} & 0_{BB}
   \end{pmatrix},

where the :math:`(B,A)` block :math:`A_{BA}^{\rm fis}` is the **fission
ray fold** — the kernel-generic
:class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicEmission`
carrying the fission kernel (the operator spelling of the coupled
fission seed's :math:`q_{1/2}` assembly).  On the eigen-:math:`M`
operator that fold **belongs** in the posing; the within-group gain
keeps it out (HAZARD 5), not the eigenproblem.  The :math:`(B,B)` block
is the genuine **zero map** ray-flux → ray-source: the :math:`w = 0`
closed rays carry no quadrature weight, so they never source fission —
spelled with the space-typed
:class:`~orpheus.numerics.operator.ZeroOperator` (``codomain_zero`` and
its dual ``transpose_zero`` hook) so both the forward grid and its
dagger emit the source-classed ray zero.  The ray-leg's
``solve_transpose`` output is duality-typed to the adjoint FLUX (the
dual of a source under the G-pairing is the adjoint flux), the exact
sibling of the within-group ``StreamingCollisionOperator.solve_transpose``
fix.

The dual lift asymmetry
-----------------------

The fixed-source entry
(:func:`~orpheus.sn.solver.solve_sn_adjoint_fixed_source`) exposes the
sharpest test of "adjoint ≠ forward-with-a-sign-flipped-source": the
source lift.  The forward isotropic-source lift divides by the
quadrature weight sum, :math:`q \mapsto \tfrac1W\,\mathbf 1_\Omega\,q`;
the **adjoint** detector-response lift does **not** —
:math:`\Sigma_d \mapsto \mathbf 1_\Omega\,\Sigma_d`, a plain angle-flat
broadcast with **no** :math:`w_n`, **no** :math:`1/W`.  The two lifts
are duals of *different* maps.  Under the bulk metric :math:`G = V\,w_n`
the plain broadcast is exactly the dual of the scalar-flux extraction,

.. math::

   \langle \mathbf 1_\Omega\,\Sigma_d,\,\psi\rangle_G
   \;=\; \sum_{\rm cells} V\,\Sigma_d\,\varphi
   \;=\; \langle \Sigma_d,\,\varphi\rangle_V ,

the detector-response functional — whereas the forward :math:`1/W` lift
is the dual of *source injection*.  This asymmetry IS the content of the
P1.2 reciprocity gate: the entries duality row cross-checks the
detector side against the hand volume sum :math:`\sum V\,\Sigma_d\,
\varphi`, pinning the angle-flat lift as exactly the adjoint of the
extraction.  (The daggered **coupled** fixed-source arm — a carrying
mesh with System B — is a typed, loud refusal at #276 A4: it has no
consumer or gate yet and lands with its first consumer rather than
shipping unexercised.  The eigenvalue entry covers carrying meshes.)

.. _sn-adjoint-carrier:

The adjoint flux carrier — ``AdjointSolution``
==============================================

Where does :math:`\varphi^*` live?  Not on the forward
:class:`~orpheus.sn.solution.Solution`.  Campaign #276 A5 split the
solution carrier along a **role axis** into sibling leaves under a
role-agnostic base
(:class:`~orpheus.sn.solution.SolutionBase` →
{:class:`~orpheus.sn.solution.Solution`,
:class:`~orpheus.sn.solution.AdjointSolution`}).

Role is a type; problem kind is a property
------------------------------------------

The solution family discriminates along **two independent axes that use
deliberately different mechanisms**:

* **Problem kind** (fixed-source vs eigenvalue) is a **property** — one
  carrier covers both via the optional ``keff``, because the two kinds
  share every realisation *and* every operation (homogenising a
  fixed-source flux is as meaningful as homogenising an eigenmode).  A
  type here would be ceremony.
* **Solution role** (forward vs adjoint) is a **type**.  The roles
  share the carrier — same fields, same packaging convention (both
  route through the one scalar- and role-agnostic ``_package_solution``
  tail) — but **not the operation set**.

The base is deliberately non-instantiable (a role-less solution is not
a value that exists); a capability-*removing* subclass
(``AdjointSolution`` inheriting ``Solution`` and hiding
``homogenize``) would violate Liskov substitutability, so the two roles
are **siblings** under the base, not parent/child.

The forward physics is structurally absent
------------------------------------------

The forward-physics operations —
:meth:`~orpheus.sn.solution.Solution.homogenize`,
:meth:`~orpheus.sn.solution.Solution.condense`,
:meth:`~orpheus.sn.solution.Solution.reaction_rate_density` — live on
:class:`~orpheus.sn.solution.Solution` **and only there**.  They
interpret ``scalar_flux`` as the flux :math:`\phi` and collapse cross
sections *preserving reaction rates* — an operation whose subject is
the forward flux.  **An importance map has no reaction rate to
preserve**, so those methods do not exist on
:class:`~orpheus.sn.solution.AdjointSolution` at all.  The absence is
**structural**, not a runtime refusal: ``adj.homogenize`` is an
``AttributeError``, and the wrong physics is *unspellable*
(``coding-elegance`` Pattern 4 — illegal states unrepresentable).  The
adjoint's ``scalar_flux`` is the importance :math:`\varphi^* = \sum_n
w_n\,\psi^*_n`, exposed under the domain-named alias
:attr:`~orpheus.sn.solution.AdjointSolution.importance` — one storage,
two vocabularies.

The design tension, ruled forward-looking
-----------------------------------------

The type-minting criterion (``coding-standards``: mint a type iff
:math:`\ge 2` non-isomorphic realisations **and** a non-identity
morphism is applied to it) technically **fails** for :math:`\varphi^*`
in isolation: it is byte-for-byte the same storage as :math:`\phi`, and
the adjoint-weighting :math:`\langle\varphi^*, \Sigma\varphi\rangle` is
a bilinear applied to the *pair* at a call site, not a morphism on
:math:`\varphi^*` alone.  The testability axis therefore favoured
leaving :math:`\varphi^*` an unmarked ``Solution``.  The USER **ruled
otherwise (#276 A5, Option 3 — mint the type)**, on the *trajectory*:
the forthcoming adjoint-machinery family — perturbation theory
:math:`\langle\varphi^*, \delta A\,\varphi\rangle`, generalised
perturbation / response estimation, and #281 adjoint-weighted
homogenisation — earns :math:`\varphi^*` a signature-level carrier that
makes its role legible at every boundary.  This is a *correctness /
forward-design* judgement that overrode the local testability
recommendation, and it is recorded here so a future reader does not
"simplify" the type back into a property and lose the intent.

The #281 adjoint-weighted homogenisation API
--------------------------------------------

The homogenize / condense asymmetry is resolved cleanly by the ratified
#281 (P6-B2) API: the forward machinery stays forward-only, and the
adjoint enters as an **optional test weight**,

.. code-block:: python

   # LANDED (P6 #281) — with its Petrov-Galerkin implementation and gates
   Solution.homogenize(coarse, *, adjoint: AdjointSolution | None = None)

— ``None`` keeps today's flux-weighted (Galerkin-degenerate) collapse
bit-identically; a real :math:`\varphi^*` makes the collapse
eigenvalue-consistent per the worth-zeroing taxonomy of the algebra of
record (:mod:`orpheus.derivations.common.homogenization` — the test
weight is the bilinear PRODUCT :math:`\varphi^*\!\odot\varphi` for the
vector channels, the exact angular pairing for :math:`\Sigma_t`, the
per-pair sink×source rule for the matrices, and the mixed-fold factored
fission rule; the same parameter on :meth:`~orpheus.sn.solution.Solution.condense`
runs the B&G-convention bilinear condensation).  The adjoint is the
*test weight* of the forward collapse, never its subject — which is
exactly why the forward trio is structurally absent on the adjoint
leaf.  The full taxonomy narrative and gates live in the frame chapter
(:ref:`sn-homogenization-why-petrov-galerkin`) and the SN verification
slice.

.. _sn-adjoint-verification:

Verification — how :math:`\psi^*` is certified
==============================================

The adjoint flux is **not** verified by MMS.  MMS is a source-driven
pillar that reaches flux-shape and convergence-order but **cannot**
verify an eigenvalue (``vv-principles`` — the pillars); the daggered
:math:`k` and :math:`\varphi^*` need **closed-form** references.  Every
Phase-1 value gate is L1, anchored to a reference that terminates in
``np.linalg.eig`` or the reciprocity identity — never in another
ORPHEUS solver.  φ\* is "correct" only when the whole chain below is
green and every named mutation reddens its named gate under
``python -O``.

The gate map
------------

.. list-table:: The adjoint certification battery (measured on Mixture A)
   :header-rows: 1
   :widths: 14 20 30 36

   * - Gate
     - Claim layer / pillar
     - Reference (structurally independent)
     - What it pins (measured)
   * - **P1.2** duality
     - model (fixed-source) / closed-form
     - the reciprocity identity :eq:`sn-adjoint-duality`; two
       independent solves
     - :math:`\langle\Sigma_d,\psi\rangle = \langle\psi^*,q\rangle` on a
       2G asymmetric-SigS vacuum slab, source and detector in
       DIFFERENT groups AND regions; detector side hand-checked against
       :math:`\sum V\Sigma_d\varphi` (pins the angle-flat lift)
   * - **P1.3** :math:`k^{\dagger}=k`
     - eigenvalue / closed-form
     - ``kinf_homogeneous`` (triple equality, terminates in
       ``np.linalg.eig``)
     - :math:`k^{\dagger}=k_{\rm fwd}=k_\infty` on ∞ 2G+4G, a
       heterogeneous reflective slab, AND the coupled sphere; teeth:
       :math:`F^{\dagger}\!\to\!F`, :math:`S^{\dagger}\!\to\!S`,
       :math:`L^{\dagger}\!\to\!L` each shift :math:`k`
   * - **P1.4** spectrum
     - eigenvalue + flux-shape / closed-form
     - the dominant right eigenvector of
       :math:`(A^{\mathsf T})^{-1}F^{\mathsf T}`
       (``kinf_and_adjoint_spectrum_homogeneous``)
     - the 4G adjoint energy spectrum :math:`\varphi^*_{\rm cf} =
       [0.470, 0.486, 0.518, 0.524]` (:math:`\ne\varphi`, asserted);
       :math:`F^{\dagger}\!\to\!F` reds the spectrum O(1)
   * - **P1.5** bi-orthogonality
     - flux-shape (intrinsic law) / closed-form
     - the spectral decomposition of :math:`M` and :math:`M^{\dagger}`
       (both from ``np.linalg.eig``)
     - the cross-Gram :math:`\langle\psi^*_i, F\varphi_j\rangle` is
       diagonal; for the rank-1 :math:`F=\chi\otimes\nu\Sigma_f` this is
       the degenerate one-nonzero-entry form (both zero mechanisms
       :math:`F\varphi_j=0` and :math:`\chi\cdot\psi^*_i=0` asserted)
   * - **sphere** :math:`\varphi^*`-shape
     - flux-shape / dense forward-probe
     - a dense FORWARD-probed :math:`(A_{\rm loss}, F)` + a raw-data
       coupled :math:`G` (both structurally independent of the ``.H``
       reverse-scan under test)
     - the coupled defining-law residual
       :math:`\|A_{\rm loss}^{\mathsf T}(G\psi^*) -
       F^{\mathsf T}(G\psi^*)/k\|` at rel floor :math:`1.2\times10^{-10}`
       vs gate :math:`10^{-7}` (:math:`n=140`); anti-vacuity
       :math:`|\Delta k| = 3.3\times10^{-11}`

The k rows verify the daggered **eigenproblem**
:eq:`sn-adjoint-eigenproblem`; the P1.2 row verifies the reciprocity
**duality** :eq:`sn-adjoint-duality`.  The full narrative and mutation
teeth live in the V&V slice (:ref:`sn-adjoint-verification-slice`); the
gate code is
``tests/sn/solve/test_sn_adjoint_certification.py`` and
``tests/sn/solve/test_sn_adjoint_entries.py``.

The Mode-12 accounting — what :math:`k` can and cannot see
----------------------------------------------------------

Because :math:`\operatorname{eig}(M^{\dagger}) =
\operatorname{eig}(M)` by construction — the identity lives on the
ITERATION operator :math:`M = A_{\rm loss}^{-1}F` (every factor
daggered; the derivation above) — a :math:`k^{\dagger}=k` gate is
**designed-green** (``vv-principles`` Mode 12) on entire classes of
error.  Getting the boundary exactly right is load-bearing — this
campaign twice caught a wrong "why" here:

* :math:`k` **is EXACTLY blind** to (i) the **factor-order / similarity
  family** (:math:`\operatorname{eig}(M^{\mathsf T}) =
  \operatorname{eig}(M)` — transposing *all* factors is a similarity),
  (ii) **all vector content**, and (iii) **the G-metric itself**
  (:math:`G'^{-1}A^{\mathsf T}G'` is metric-similar to
  :math:`A^{\mathsf T}` for any invertible :math:`G'`).  No tolerance,
  mesh refinement, or regime change can expose these through
  :math:`k`; the committed catchers are the spectrum row, the
  bi-orthogonality row, the duality pairing, and the sphere vector row
  — functionals **outside** the eigenvalue stabiliser.
* :math:`k` **is NOT blind** to a single **leaf-transpose drop**
  (:math:`F^{\dagger}\!\to\!F`, :math:`S^{\dagger}\!\to\!S`,
  :math:`L^{\dagger}\!\to\!L`): transposing *one* factor is **not** a
  similarity of the pencil, and :math:`k` measurably moves —
  :math:`F^{\dagger}=F` shifts :math:`k` from :math:`1.488` to
  :math:`0.171` on the 4G ∞ fixture (the FULL SN-solve measurement; the angular-collapsed 0-D closed-form proxy of the same mutation gives 0.153 — cite the solve, not the proxy).  So the k-equality rows *are*
  legitimate teeth for drops (with the visibility preconditions:
  asymmetric SigS, :math:`\chi\not\parallel\nu\Sigma_f`, spatial
  structure), while the factor-order and metric classes stay the vector
  rows' exclusive province.

**The factor-order trap.**  The P1.4 reference must be the dominant
right eigenvector of :math:`(A^{\mathsf T})^{-1}F^{\mathsf T}`, **not**
:math:`\operatorname{eig}(M^{\mathsf T})`.  The two are similar
(conjugation by :math:`A^{\mathsf T}`) so every :math:`k` check passes
on both — but for the rank-1 :math:`F`, the dominant eigenvector of
:math:`M^{\mathsf T} = F^{\mathsf T}A^{-\mathsf T}` degenerates to
**exactly** :math:`\widehat{\nu\Sigma_f}` (:math:`F^{\mathsf T}x \propto
\nu\Sigma_f` for all :math:`x`), a reference carrying **zero A-physics**.
The wrong reference was caught by the SN daggered solve on first contact
(the corrected law is recorded in
:func:`~orpheus.derivations.common.eigenvalue.kinf_and_adjoint_spectrum_homogeneous`).
And **the metric is caught by nothing but the sphere vector row**: the
:math:`G_{\rm sd} = V_{\rm cell} \to 1` drop leaves
:math:`|k^{\dagger}_{\rm mut}-k_{\rm fwd}| = 2.6\times10^{-11}`
(EXACTLY k-blind) while the residual reds to :math:`2.35` (the ERR-067
family — a metric bug in production that no eigenvalue gate could see).

.. _sn-adjoint-consumers:

Consumers and horizon
=====================

The adjoint flux is a means, not an end.  With :math:`\varphi^*` landed,
its consumers are unblocked:

* **Adjoint-weighted (eigenvalue-consistent) homogenisation** (#281 P6,
  frame-machinery :ref:`galerkin-projection` / #51).  Flux-weighted
  homogenisation is the Galerkin degenerate :math:`\varphi^*=\varphi`;
  the adjoint weight makes it genuinely Petrov-Galerkin, so the coarse
  :math:`k` becomes **first-order stationary** — the homogenisation
  error is :math:`\mathcal O(\delta\varphi^2)` rather than
  :math:`\mathcal O(\delta\varphi)`.  This lands with its C1/C2/C3 gates
  in P6 via the ``Solution.homogenize(..., adjoint=...)`` parameter
  (:ref:`sn-adjoint-carrier`).
* **Perturbation theory and GPT.**  The first-order worth
  :math:`\delta k \propto \langle\varphi^*, \delta A\,\varphi\rangle /
  \langle\varphi^*, F\varphi\rangle` and generalised perturbation /
  response-sensitivity estimation are the adjoint eigenmode's native
  applications — the trajectory the A5 type ruling was made for.

The **honest deferral ledger.**  Two arms are callable-but-deferred by
design, each a typed refusal rather than an unexercised path:

* the daggered **coupled fixed-source** arm (a carrying mesh with
  System B) refuses loud in
  :func:`~orpheus.sn.solver.solve_sn_adjoint_fixed_source` — no consumer
  or gate yet; the eigenvalue entry covers carrying meshes;
* the Gauss–Seidel **schedule-reverse** transpose (#310 R7) has no
  consumer, so a
  :class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`
  over it stays non-adjointable (:ref:`loss-rep-orientation-two-frames`,
  the deferral ledger).

Development history
===================

* **#276 P3** (commit ``15185e5``, closes #118) — the **scattering
  adjoint** :math:`S^{\mathsf T}`, free from the harmonic frame.  The
  modernised in-scatter source became one frame-conjugated operator, so
  its transpose falls out as the product transpose with no per-geometry
  derivation.  :class:`~orpheus.transport.operators.scattering.ScatteringOperator`
  gained ``apply_transpose`` and now reports ``is_adjointable = True``.
* **#280 / #310** — the **walk adjoint** / orientation axis.  #280 added
  the reverse-scan inner solve :math:`(L+C)^{-\mathsf T}b` and the swap
  law :math:`(A^{\dagger})^{-1} = (A^{-1})^{\dagger}`; #310 completed the
  ``loss_action_transpose`` grid over every registered
  scheme × representation (DD 1-D/2-D/3-D, LD 1-D/2-D), retiring the
  transpose residue (:ref:`loss-rep-orientation-two-frames`).
* **#276 A4** (merged @ ``065a0e5d``) — the **daggered posing
  activation**.  ``KEigenvalue((L+C).H, (S+B).H, F.H)`` runs through the
  unchanged ``power_iteration``; the entries
  :func:`~orpheus.sn.solver.solve_sn_adjoint` /
  :func:`~orpheus.sn.solver.solve_sn_adjoint_fixed_source` land; the
  coupled sphere posing (fission ray fold + space-typed
  ``ZeroOperator`` + duality-typed ``solve_transpose``) lands; the P1.4
  reference is corrected from the factor-order-degenerate
  :math:`\operatorname{eig}(M^{\mathsf T})` to
  :math:`(A^{\mathsf T})^{-1}F^{\mathsf T}`.
* **#276 A5** (merged @ ``a24380ca``) — the **role-typed carrier**.
  :class:`~orpheus.sn.solution.SolutionBase` →
  {:class:`~orpheus.sn.solution.Solution`,
  :class:`~orpheus.sn.solution.AdjointSolution`}, the forward physics
  structurally absent on the adjoint, ``importance`` alias, and the
  coupled-sphere :math:`\varphi^*`-shape row closing the sphere
  :math:`k` row's honest-scope gap.
* **#276 A6** — this chapter: the route decision, the three-transposes
  taxonomy, the daggered-posing mechanics, the carrier ruling, and the
  verification narrative, closing #276.

References
==========

The adjoint transport equation, importance interpretation, and
reciprocity duality follow :cite:`BellGlasstone1970` (§6, importance
and the adjoint) and :cite:`LewisMiller1984` (§6, adjoint transport,
reciprocity, and perturbation theory); the Monte-Carlo importance
reading is :cite:`Lux1991`.  The discrete Hilbert (G-metric) adjoint is
derived in :ref:`operator-adjoint`; the reverse-scan walk transpose and
the swap law in :ref:`loss-rep-orientation-two-frames`.
