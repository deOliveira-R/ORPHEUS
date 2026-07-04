.. _loss-representations:

==================================================================
Selectable Representations of the S\ :sub:`N` Loss Operator
==================================================================

This is the **capstone architecture page** for the within-group
S\ :sub:`N` loss operator :math:`(L+C)` — the streaming-plus-collision
operator whose inverse *is* the transport sweep. It documents the
final state of issue #222 (the *sweep-strategy carve* + the S6
operator/representation re-layering, complete 2026-06-11 at the
Fork-B2 default flip): one lower-triangular operator, several
*algorithms* that realise it, a single source of truth for which
algorithm a mesh gets, and the L21 theorem that makes "matvec is the
same operator as the sweep" a type fact rather than a coincidence.

It deliberately does **not** re-derive the cell-level mathematics.
Those live in :doc:`discrete_ordinates`:

* the WDD cell balance and closure — :eq:`dd-cartesian-1d`,
  :eq:`dd-cartesian-2d`, :eq:`dd-2d-balance-form`;
* the 1-D cumprod recurrence — :ref:`sweep-cumprod`;
* the 2-D anti-diagonal wavefront — :ref:`sweep-wavefront`;
* the three-layer **walk / level-op / kernel** stack — the
  :ref:`sweep-dispatch-relayering` section (S6.4(e)).

…and the operator algebra lives in :doc:`operator_algebra`:

* the typed forward/solve/adjoint actions —
  :eq:`operator-apply`, :eq:`operator-solve`,
  :eq:`operator-apply-transpose`;
* the interior face-flux cochain :math:`C^1_{\rm int}` and its
  succession after the typed ``WavefrontFlux`` retirement —
  :ref:`wavefront-flux-cochain`.

This page is the layer *above* both: the **representation** layer that
selects and unifies them.

.. contents:: On this page
   :local:
   :depth: 2


Key Facts
=========

.. admonition:: Key Facts — the loss-representation architecture
   :class: important

   * **One operator, several schedules.** :math:`(L+C)` is
     **lower-triangular** under the upwind (sweep) cell ordering. Its
     two actions are one object: :math:`SOLVE = (L+C)^{-1}q` is forward
     substitution (the *transport sweep*); :math:`APPLY = (L+C)\psi`
     is the row action (the *Krylov matvec*). A
     :class:`~orpheus.sn.loss_representation.LossRepresentation` is a
     **schedule** for traversing that triangular structure, not a
     different operator.

   * **The four representations**
     (:mod:`orpheus.sn.loss_representation`), each a stateless frozen
     ``@dataclass`` carrying only the mesh:

     - :class:`~orpheus.sn.loss_representation.CumprodScan` — the 1-D
       Blelloch parallel-prefix scan (slab + sphere + cylinder via one
       body); the **1-D production default**.
     - :class:`~orpheus.sn.loss_representation.ScanMarch` —
       :math:`\mathrm{scan}(x)\circ\mathrm{march}(y)`; the **multi-D
       Cartesian production default since the S6.9 Fork-B2 flip**.
     - :class:`~orpheus.sn.loss_representation.MovingFrontierWindow` —
       the anti-diagonal wavefront carrying a rolling
       :math:`(d{-}1)`-frontier; a **selectable peer** (the d=2 default
       through S6.9).
     - :class:`~orpheus.sn.loss_representation.FullFieldWavefront` —
       the same DAG schedule retaining the **whole** interior cochain;
       the **verification oracle**, explicit-select only.

   * **Selection is a single source of truth.**
     :meth:`~orpheus.sn.loss_representation.LossRepresentation.supports`
     returns :class:`~orpheus.sn.loss_representation.Compatibility`
     ``(ok, reason)``; :func:`~orpheus.sn.loss_representation.default_for`
     returns the first match in the ordered registry
     :data:`~orpheus.sn.loss_representation.LOSS_REPRESENTATIONS`.
     Illegal ``(representation, mesh)`` pairings are unrepresentable —
     ``__post_init__`` re-checks ``supports`` and raises
     :class:`~orpheus.sn.loss_representation.IncompatibleRepresentation`.

   * **L21 — one walk, one instance.** A single d-generic
     ``_OctantWalk._interior_walk`` frame serves sweep AND matvec for
     every multi-D representation, forked only by a **kernel object** +
     **emit policy** (never a boolean). The operator holds ONE
     representation instance consumed by ``apply``, ``solve``, and the
     Gauss–Seidel resolvent. Both are *type facts*, pinned by spy tests.

   * **The Fork-B2 decision (the WHY of the default).** ScanMarch is
     **0.55–0.84×** the window's per-sweep time at **identical** peak
     memory; the window's memory claim only ever held against the
     full-field oracle (1.3–1.4× both). The window was **kept**, not
     retired — a genuinely different schedule over the same operator is
     the point of selectability.

   * **Bit-identity vs principled-equivalence.** Different *schedules*
     are NOT bit-comparable (FP-association differs by construction) →
     ``nulp`` / solver-tol Mode-9 gates. Different *storage policies* of
     the same schedule ARE bit-identical → ``array_equal`` oracles. The
     cell kernel's explicit left fold ``((Σ_t + s_0) + s_1)`` is the FP
     reduction tree of record (sha256-source-pinned).

   * **Gotcha — the operator subtracts** :math:`C` **once.**
     ``loss_action`` MUST return the **full** :math:`(L+C)\psi`, NOT
     :math:`L\psi`. The operator applies the only glue
     :math:`L = (L+C) - C` (Resolution A). A leaf returning
     :math:`L\psi` would double-count the collision diagonal.


.. _loss-rep-native-frame:

The native mathematical frame: a lower-triangular operator
==========================================================

The within-group discrete-ordinates balance, for one ordinate
:math:`\Omega_n` and one cell, is the WDD relation derived in
:doc:`discrete_ordinates` (:eq:`dd-cartesian-2d`). Collect every cell
and every ordinate into the within-group operator

.. math::
   :label: loss-rep-LpC

   (L+C)\,\psi \;=\; q,
   \qquad
   L \;=\; \Omega\cdot\nabla\big|_{\rm WDD},
   \quad
   C \;=\; \sigma_t\,\odot,

where :math:`L` is the discretised streaming operator (the upwind
WDD difference relations) and :math:`C` is the collision diagonal
(multiply-by-:math:`\sigma_t`). The decisive structural fact:

.. admonition:: :math:`(L+C)` is lower-triangular under the sweep ordering
   :class: note

   Order the unknowns cell-by-cell in **upwind (inflow-to-outflow)
   order** for a given ordinate. Then cell :math:`i`'s balance depends
   only on faces that are *already known* — the domain inflow, or the
   outflow of strictly-upstream cells. In that ordering the matrix of
   :math:`(L+C)` is **lower-triangular** (block-lower-triangular over
   groups). A lower-triangular system has two canonical operations,
   and they are the same operator viewed two ways:

   .. list-table::
      :header-rows: 1
      :widths: 14 30 30 26

      * - Operation
        - Linear-algebra name
        - S\ :sub:`N` name
        - Code entry
      * - :math:`SOLVE`
        - forward substitution :math:`(L+C)^{-1}q`
        - the **transport sweep**
        - ``LossRepresentation.sweep``
      * - :math:`APPLY`
        - the row action :math:`(L+C)\psi`
        - the **Krylov matvec**
        - ``LossRepresentation.loss_action``

This is the **L21 invariant** referenced throughout the SN code: *the
sweep and the matvec are two actions of the SAME operator* — they are
not independent code paths that happen to agree, and they must never be
allowed to drift. Forward substitution and the row action share the
triangular factor; they differ only in which is the unknown.

Resolution A — the operator's only glue
---------------------------------------

The representation returns the **full** within-group loss action
:math:`(L+C)\psi`. The operator
(:meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply`,
:eq:`operator-apply`) then applies the *only* remaining algebra glue,

.. math::
   :label: loss-rep-resolution-a

   L\,\psi \;=\; (L+C)\,\psi \;-\; \sigma_t \odot \psi.\mathrm{bulk}
   \qquad\text{(Resolution A: } L = (L+C) - C\text{)},

subtracting the collision diagonal :math:`C = \sigma_t\,\odot`
**exactly once**. Before the S6.3 re-layering this :math:`-C`
subtraction was duplicated five times across five private ``_apply_*``
bodies; collapsing it to one site is a single-source-of-truth win, but
it imposes the **load-bearing contract** on every representation:

.. warning::

   ``loss_action`` MUST return :math:`(L+C)\psi`, **not** :math:`L\psi`.
   A leaf that returned the bare streaming action :math:`L\psi` would
   make the operator subtract :math:`C` a *second* time — a
   double-counted collision diagonal, a silent sign-and-magnitude error
   (vv-principles failure Mode 3, *missing/duplicated factor*). The
   convention is pinned by ``test_loss_action_convention.py``: the
   non-tautological anchor checks that for a flat reflective field
   :math:`L\psi_{\rm flat} = 0`, so :math:`(L+C)\psi = \sigma_t\psi` —
   proving the action is the *full* :math:`(L+C)` loss, not bare
   :math:`L`, and cross-checks the :math:`-C` glue against an
   independent collision multiplier :math:`C = M[\sigma_t]`
   (a :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`).

.. _loss-rep-removal-form-matvec:

The composite owns its matvec — the removal form :math:`C(\sigma_r)`
------------------------------------------------------------------------

So far the matvec :math:`(L+C)\psi` has been described as the matvec twin of
the sweep, applied by :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply`
(:eq:`loss-rep-resolution-a`). But the *composite* operator
:class:`~orpheus.sn.operators.streaming.InvertibleOperator` — the :math:`A = (L+C)`
object whose :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.solve` *is* the
sweep — has its own matvec too, and getting that matvec to single-source its
diagonal :math:`\sigma` from the *same* place :meth:`solve` does is the
substance of issue #240 Phase 2 Step B. This subsection records that carve: a
**principled-equivalence** refactor (not a bug fix — no wrong value ever
shipped) that turns a value-correct-by-coincidence leaf sum into a
correct-by-construction action, and in doing so foreclosed a latent trap that
the **removal form** :math:`\mathrm{InvertibleOperator}(L(\sigma_t),
C(\sigma_r))` would have made observable.

The :math:`\sigma` parameter is now explicit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The protocol signature of the two matvec actions changed: both
:meth:`~orpheus.sn.loss_representation.LossRepresentation.loss_action` and
:meth:`~orpheus.sn.loss_representation.LossRepresentation.loss_action_transpose`
now take the collision diagonal :math:`\sigma` **explicitly** as their first
argument (``loss_action(sigma, psi)``), exactly symmetric with the sweep door
``sweep(Q, sig_t, ...)``. Before the carve the representation read
``operator.sigma_t`` *off an operator handle* it was passed — so the leaf, not
its caller, decided which :math:`\sigma` the matvec realised. After the carve
the **caller single-sources** :math:`\sigma`:

* :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply` passes its own
  ``sigma_t`` (and subtracts it back via Resolution A, :eq:`loss-rep-resolution-a`);
* the **new** :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.apply` /
  :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.apply_transpose` overrides pass
  the composite's *own* diagonal ``self.sigma`` — the SAME array
  :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.solve` threads into the WDD
  sweep — and realise :math:`M(\sigma)\psi` *directly*, via
  ``self.loss_representation.loss_action(self.sigma, psi)``, instead of inheriting
  the :meth:`~orpheus.numerics.operator.OperatorSum.apply` leaf sum.

This is the **one-instance theorem** (:ref:`loss-rep-one-walk-one-instance`)
extended to the diagonal: the operator already held ONE representation instance
shared by :meth:`apply` and :meth:`solve`; now it also holds ONE
:math:`\sigma`, supplied by the same operand whichever door is exercised. The
removal form is what makes that single-sourcing *matter*.

The affine-in-:math:`\sigma` structure of the forward matvec
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The inherited :meth:`~orpheus.numerics.operator.OperatorSum.apply` is the
literal leaf sum :math:`L.\mathrm{apply}(\psi) + C.\mathrm{apply}(\psi)`. The
crux of the carve is that this leaf sum is value-correct **only by
coincidence** — a coincidence that follows from a structural property of the
*forward* (apply) direction that does **not** hold in the sweep direction.

Start from the WDD cell balance derived in :doc:`discrete_ordinates`. In the
forward (apply) direction the cell average :math:`\bar\psi` is **known**: it is
the input. The matvec is the *row action* — it reconstructs the residual that
the balance would have to absorb,

.. math::
   :label: loss-rep-affine-cell

   M(\sigma)\,\psi\big|_{\rm cell}
   \;=\;
   \mathrm{denom}\cdot\bar\psi \;-\; \mathrm{numer}_{\rm upstream},
   \qquad
   \mathrm{denom} \;=\; \underbrace{\textstyle\sum_a c_a}_{\text{streaming}}
                        \;+\; \sigma,

with :math:`c_a = 2|\mu_a|/\Delta a` the per-axis scheme-scaled streaming
coupling (:eq:`loss-rep-scanmarch-solve`; the diamond :math:`2` is the
scheme's). Because :math:`\bar\psi` is known,
:math:`\sigma` enters :eq:`loss-rep-affine-cell` only through the *additive*
term :math:`\mathrm{denom}\cdot\bar\psi = (\sum_a c_a)\bar\psi +
\sigma\bar\psi`. Collecting the :math:`\sigma`-independent part into a
``streaming_action``,

.. math::
   :label: loss-rep-affine

   M(\sigma)\,\psi
   \;=\;
   \mathrm{streaming\_action}(\psi) \;+\; \sigma\cdot\psi
   \qquad\text{(forward matvec is AFFINE in }\sigma\text{).}

This is the decisive fact. **The matvec is affine in** :math:`\sigma`: a clean
additive :math:`+\,\sigma\cdot\psi`, never a :math:`1/\mathrm{denom}`.

Contrast the **sweep** direction, where :math:`\bar\psi` is the *unknown*. The
cell-average solve is :math:`\bar\psi = \mathrm{numer}/\mathrm{denom}` — and
now :math:`\sigma` sits in the **denominator** (:math:`\mathrm{denom} = \sum_a
c_a + \sigma`), so :math:`\bar\psi` is a *rational* function of :math:`\sigma`,
**non-affine**. The whole non-affinity of the loss operator in :math:`\sigma`
lives in the sweep direction; the apply direction is affine because it never
inverts the denominator. (This asymmetry is why the round-trip
:math:`\mathrm{apply}\circ\mathrm{solve}` cannot detect a :math:`\sigma`-routing
error in ``apply`` alone — see the verification subsection below.)

The affine structure :eq:`loss-rep-affine` is precisely what makes the
inherited leaf sum value-correct. Each leaf reads its **own** diagonal:
:math:`L.\mathrm{apply}` returns :math:`M(\sigma_t)\psi - \sigma_t\cdot\psi =
\mathrm{streaming\_action}(\psi)` (Resolution A subtracts :math:`L`'s own
:math:`\sigma_t`), and :math:`C.\mathrm{apply}` returns :math:`\sigma_r\cdot\psi`
(the collision leaf's own :math:`\sigma_r`). Summing,

.. math::
   :label: loss-rep-leaf-sum

   \underbrace{L.\mathrm{apply}(\psi)}_{M(\sigma_t)\psi - \sigma_t\psi
                                       \;=\; \mathrm{streaming\_action}(\psi)}
   \;+\;
   \underbrace{C.\mathrm{apply}(\psi)}_{\sigma_r\cdot\psi}
   \;=\;
   \mathrm{streaming\_action}(\psi) + \sigma_r\cdot\psi
   \;=\;
   M(\sigma_r)\,\psi.

So the leaf sum **does** compute :math:`M(\sigma_r)\psi`, the right value for
the removal form — but it gets there by sourcing :math:`\sigma_t` from
``L.sigma_t`` (the streaming leaf), cancelling it, and then *re-adding*
:math:`\sigma_r` from ``C``. The value is right; the **source is wrong** —
:math:`\sigma` is sourced from a *different operand* than the one
:meth:`solve` uses, and only the affine structure :eq:`loss-rep-affine` keeps
that from showing up as a wrong number.

The two-source hazard
~~~~~~~~~~~~~~~~~~~~~~~

The latent defect is a **twin-source**. The composite's two doors source the
diagonal :math:`\sigma` from two *different* operands:

.. list-table:: Where the two doors source :math:`\sigma` (pre-carve)
   :header-rows: 1
   :widths: 24 38 38

   * - Door
     - Sources :math:`\sigma` from
     - Realises
   * - :meth:`InvertibleOperator.solve`
     - the collision leaf ``C`` (``self.sigma``)
     - :math:`(L+C(\sigma_C))^{-1}q` (the sweep)
   * - inherited ``OperatorSum.apply``
     - the streaming leaf ``L`` (``L.sigma_t``)
     - :math:`\mathrm{streaming\_action} + \sigma_C\cdot\psi`
       via :eq:`loss-rep-leaf-sum`

They **agree** today only because production always builds :math:`L` and
:math:`C` from the *same* :math:`\sigma_t` (``L.sigma_t == C.sigma``). The
moment they differ — the removal form — the two doors realise actions sourced
from different diagonals, held in agreement purely by the affine coincidence
:eq:`loss-rep-affine`. This is the ``coding-elegance`` Smell of a *value
correct by coincidence rather than by construction*: the inherited
``OperatorSum.apply`` is a twin-path realisation of the composite's own
action, re-deriving the streaming part through :math:`\sigma_t` only to cancel
it.

The override (#240 Step B) collapses the twin-source. Both
:meth:`~orpheus.sn.operators.streaming.InvertibleOperator.apply` and
:meth:`~orpheus.sn.operators.streaming.InvertibleOperator.apply_transpose` now read the
composite's **own** ``self.sigma`` — the SAME array
:meth:`~orpheus.sn.operators.streaming.InvertibleOperator.solve` threads into the sweep —
and call ``loss_action(self.sigma, psi)`` directly:

.. code-block:: python

   class InvertibleOperator(OperatorSum):
       def apply(self, psi):                         # (L+C)ψ = M(σ)ψ
           _require_typed_composite("InvertibleOperator.apply", self.sn_mesh, psi)
           return self.loss_representation.loss_action(self.sigma, psi)

       def apply_transpose(self, phi):               # (L+C)ᵀφ = M(σ)ᵀφ
           _require_typed_composite(..., self.sn_mesh, phi)
           return self.loss_representation.loss_action_transpose(self.sigma, phi)

This is ``coding-elegance`` Pattern 2 (single source of truth): **one**
``loss_action``, **one** source of :math:`\sigma`, consumed identically by
both directions. The composite never asks a leaf for a :math:`\sigma`-bearing
action it must then undo. The input contract (typed composite + the
mesh-identity invariant) is itself single-sourced through the module-level
``_require_typed_composite`` helper that :meth:`StreamingOperator.apply` now
shares. The multi-D Cartesian *adjoint* still raises the representation's
deferral (``ScanMarch.loss_action_transpose`` is a ``NotImplementedError`` —
the reverse-mode multi-D sweep is the O.2b follow-on) rather than silently
routing around it, so :meth:`apply_transpose` is correct-by-construction or a
loud refusal, never a silent wrong answer.

The removal form makes the distinction observable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The **removal form** is the within-group regime that makes the two-source
distinction visible. It folds the within-group self-scatter into the
collision diagonal,

.. math::
   :label: loss-rep-removal-sigma

   \sigma_r \;=\; \sigma_t \;-\; \Sigma_{s,0}^{\,g\to g},

so the composite becomes :math:`\mathrm{InvertibleOperator}(L(\sigma_t),
C(\sigma_r))` with :math:`\sigma_C = \sigma_r \neq \sigma_t` (Adams & Larsen
2002 §III — the within-group iteration framing in which the sweep
:math:`(L+C)^{-1}` *is* the transport operator; #200). On this form the
sweep solves :math:`(L+C(\sigma_r))^{-1}`, and the matvec MUST realise the
matching :math:`M(\sigma_r)\psi` — its inverse twin (L21). With
:math:`\sigma_t \neq \sigma_r` the question "**which** :math:`\sigma` does the
matvec realise?" finally has two different answers, and the override is the one
that single-sources :math:`\sigma_r` from ``C`` exactly as :meth:`solve` does.

There is **no production caller of the removal form yet**: the consumer that
would build :math:`\sigma_r` is the within-group self-scatter fold of issue
#200, which is not wired. The collision multiplier :math:`C = M[\sigma]`
(a :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`)
already accepts either :math:`\sigma_t` or :math:`\sigma_r` (it carries no
interpretation flag — both are :math:`(\mathrm{ng}, \ldots)` arrays applied as
:math:`\sigma\cdot\psi`), and a :math:`\sigma_r`-*sweep* is **not** a correct
within-group accelerator for anisotropic flux (it inverts the
diagonal-in-angle removal :math:`\sigma_r\,I`, not the isotropic-projection
self-scatter :math:`\Sigma_{s,0}P_{\rm iso}` — issue #215, 46–56 % errors on
vacuum/heterogeneous). The carve is therefore *prophylactic*: the override
forecloses the latent twin-source **before** any consumer can trip it. It
would have become a genuine numerical trap only if a future refactor made
``L.apply`` *non-affine* in :math:`\sigma` (e.g. an :math:`L`-leaf that itself
inverts a denominator) — at which point :eq:`loss-rep-leaf-sum` would no
longer collapse to :math:`M(\sigma_r)\psi`. The override removes that coupling
by construction.

Verification — value-preserving re-association, structurally cross-checked
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because the carve preserves the matvec *value* and only re-associates the
floating-point reduction, it is verified as a #240-Phase-2-Step-A-shape
change (a value-preserving re-association), against the
:ref:`vv-principles three criteria <loss-rep-bit-vs-principled>`. The gate is
``tests/sn/operators/test_removal_form_matvec_sweep.py``, ``foundation``-tagged
and ``verifies("loss-rep-resolution-a")`` (the carve *sharpens* the existing
Resolution-A equation rather than minting a new label — the composite-owns-its-matvec
convention is the same :eq:`loss-rep-resolution-a` glue, now single-sourcing
:math:`\sigma`). It has four parts:

(a) **The structural teeth** —
    ``test_invertible_apply_is_M_of_C_sigma_bit_identical`` (and its adjoint
    sibling, slab/sphere/cyl) demand
    :math:`(L+C).\mathrm{apply}(\psi) = M(\sigma_r)\psi`
    **bit-identically** (``array_equal``, 0 ULP) against a structurally
    *independent* reference: a SEPARATE :class:`~orpheus.sn.operators.streaming.StreamingOperator`
    whose **own** :math:`\sigma_t` *is* :math:`\sigma_r`, so its
    ``loss_action(σ_r, ψ)`` is unambiguously :math:`M(\sigma_r)\psi` (no
    removal/leak ambiguity — the leaf reads its own diagonal). Under the
    override, ``op.apply`` IS that same ``loss_action`` call on the same walk →
    0 ULP. Under the pre-fix leak it is the leaf sum :eq:`loss-rep-leaf-sum`:
    value-equal but a *different* reduction tree → :math:`\le 2` ULP off
    ``array_equal`` → it fails. **These are the only gates that discriminate
    the carve.** They were ``xfail(strict=True)`` until the override landed;
    the marker is removed and they now pass plainly. The 2-D adjoint is
    *excluded* because ``ScanMarch.loss_action_transpose`` is the deferral
    raise.

(b) **The value ground** — the matvec≡sweep round-trip
    (``test_removal_form_matvec_sweep_roundtrip``: slab/cyl/2-D, sphere
    excluded) and the affine-in-:math:`\sigma` value characterisation
    (``test_removal_form_apply_value_equals_M_of_sigma_r``). These PASS under
    BOTH the leak and the override — the leak is value-correct
    (:eq:`loss-rep-affine`), so the round-trip
    :math:`\mathrm{apply}\circ\mathrm{solve}\approx\mathrm{id}` cannot tell
    leaky from override (a vv Mode-9 degeneracy *in the round-trip gate's
    design* — the teeth therefore live in (a), not here). Sphere is excluded
    from the bare round-trip because the curvilinear Morel–Montry sweep reads
    the previous iterate for the Carlson coupled-pole closure (ERR-058
    family), so a single :meth:`solve` is not the one-shot inverse of
    :meth:`apply`; the sphere matvec≡sweep claim rides the teeth gate (a) plus
    the production-:math:`\sigma` fixed-point bridge.

(c) **The production-:math:`\sigma` invariant** —
    ``test_production_sigma_apply_value_preserved`` builds the *production*
    composite (:math:`\sigma_C = \sigma_t`) and an explicit
    ``L.apply + C.apply`` leaf-sum reference, and asserts they agree to
    ``nulp=6``. The override **drops** the redundant
    :math:`(x - \sigma\psi) + \sigma\psi` round-trip the leaf sum carries
    (:math:`L.\mathrm{apply}` returns :math:`\mathrm{loss\_action}(\sigma_t) -
    \sigma_t\psi`, then :math:`C.\mathrm{apply}` adds :math:`\sigma_t\psi`
    back), so the override is the **more accurate** path. The drift is
    FP-non-associativity (reduction-depth :math:`\times` ULP), measured
    :math:`\le 2` ULP on slab/sphere (the 1-D ``_OneDimScanWalk``),
    :math:`\le 4` on cylinder, :math:`\le 5` on 2-D Cartesian (the
    ``_OctantWalk``). The gate's tolerance was itself re-baselined off a
    ``bitexact=True`` spec: that spec came from a probe comparing two
    *override-style* paths (0 ULP), but the real gate compares
    override-vs-leaf-sum, where bit-identity is *not* a property of the pair
    (the override removes the round-trip). All three vv criteria hold — named
    intermediate (``loss_action``), verified against :math:`M(\sigma_r)` by the
    teeth (a), drift bounded by reduction depth.

(d) **The negative control** —
    ``test_removal_form_nonpositive_sigma_r_rejected`` confirms
    :math:`\mathrm{InvertibleOperator}(L, C(\sigma_r\le 0))` raises at
    construction (vv L11: a contract-validation path needs BOTH a positive
    must-not-raise case — the removal cases — AND a negative must-raise case).

The eigenvalue cross-check (the closed-form :math:`\kinf = \nu\Sigma_f /
\Sigma_a` for a homogeneous reflective medium with the within-group
self-scatter folded into :math:`\sigma_r`,
``test_removal_form_kinf_independent_reference_2g``) is the *structurally
independent* reference required by vv criterion 2 for the removal regime; it is
``xfail`` until the #200 solver entry exists. The blast-radius split is itself
the structural-independence evidence that the apply re-baseline is principled:
APPLY-path snapshots (the cyl/2-D matvec golden) re-associate by :math:`\le 5`
ULP, but SWEEP/SOLVE snapshots (slab/sphere) stay **bit-identical** — they ride
:meth:`solve`, which never touched the override. And the 2-D Krylov converged
flux (which *does* drive the override through GMRES) agrees with the 2-D SI
converged flux (which does **not** use the override) to :math:`3.9\times10^{-12}`
relative — two structurally distinct iteration paths landing on the same fixed
point, the cross-check that the apply-direction re-baseline did not move the
physics.

The two actions bottom out, for every Cartesian representation, in **one
pure kernel pair** on
:class:`~orpheus.transport.spatial.diamond.DiamondDifference`:

* :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.cell_kernel_batch`
  (solve) — :math:`\bar\psi = \mathrm{numer}/\mathrm{denom}`;
* :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.residual_kernel_batch`
  (apply) — :math:`r = \mathrm{denom}\cdot\bar\psi - \mathrm{numer}`,

with the explicit **left fold**
:math:`((\Sigma_t + s_0) + s_1) + \cdots` as the IEEE-754 reduction
tree of record (see :ref:`loss-rep-bit-vs-principled`).


.. _loss-rep-four:

The four representations — four schedules of one operator
=========================================================

A :class:`~orpheus.sn.loss_representation.LossRepresentation` is a
stateless frozen dataclass (its only field is the
:class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` it was selected for). Each is a
distinct **schedule** over the same lower-triangular
:math:`(L+C)` — a different topological linearisation of the identical
cell dependencies. They are *algorithms*, not operators.

.. list-table:: The four loss representations
   :header-rows: 1
   :widths: 22 16 22 20 20

   * - Representation
     - ``supports``
     - Schedule
     - Storage
     - Role
   * - :class:`~orpheus.sn.loss_representation.CumprodScan`
     - 1-D, any geometry
     - Blelloch prefix scan
     - chain recurrence
     - 1-D **production default**
   * - :class:`~orpheus.sn.loss_representation.ScanMarch`
     - 1-D OR 2-D Cartesian
     - ``scan(x) ∘ march(y)``
     - :math:`(d{-}1)`-slab per line
     - 2-D **production default**
   * - :class:`~orpheus.sn.loss_representation.MovingFrontierWindow`
     - Cartesian, d = 2
     - anti-diagonal wavefront
     - rolling :math:`(d{-}1)`-frontier
     - selectable **peer**
   * - :class:`~orpheus.sn.loss_representation.FullFieldWavefront`
     - Cartesian, any d
     - anti-diagonal wavefront
     - full interior cochain
     - verification **oracle**


CumprodScan — the 1-D chain
---------------------------

A 1-D mesh is a total order, so the sweep is a **chain prefix scan**.
:class:`~orpheus.sn.loss_representation.CumprodScan` is intrinsically
1-D — a prefix scan needs a total order, and there is no "2-D chain
scan" (this is *legitimate* d-specificity by the algorithm's nature,
not a narrow crutch). The geometry difference (slab vs sphere vs
cylinder) is absorbed by the two-stratum sweep cache, so all three
geometries share **one body**
(:meth:`~orpheus.sn.loss_representation._OneDimScanWalk.sweep` →
:func:`~orpheus.sn.spatial.scan.ordinate_scan`); the curvilinear
Morel–Montry angular redistribution folds into the scan's affine
source. The recurrence and its closed-form cumprod/cumsum solution are
derived at :ref:`sweep-cumprod`; the pair-monoid algebra that justifies
the closed form is documented on
:func:`~orpheus.sn.spatial.scan.ordinate_scan`.

CumprodScan is conditioning-robust by construction: the closed-form
backend handles the pole reset (ERR-054) and denormal underflow
(ERR-057) that a naive cumprod would lose to gradual underflow.


ScanMarch — the row-march that reuses the scan
----------------------------------------------

:class:`~orpheus.sn.loss_representation.ScanMarch` reframes the d-D DD
sweep as forward substitution along the sweep axis — *the same
first-order linear scan* — **marched over the transverse axes**:

.. math::
   :label: loss-rep-scanmarch

   \mathrm{ScanMarch}
   \;=\;
   \begin{cases}
     \mathrm{scan}(x) & d = 1 \quad(\text{the } s_y = 0 \text{ degeneration})\\[2pt]
     \mathrm{scan}(x)\circ\mathrm{march}(y) & d = 2 \\[2pt]
     \mathrm{scan}(x)\circ\mathrm{march}(y,z) & d = 3.
   \end{cases}

Within each transverse row the diamond-difference x-face recurrence is
the **same Blelloch scan** that
:class:`~orpheus.sn.loss_representation.CumprodScan` uses
(:func:`~orpheus.sn.spatial.scan._scanmarch_row`); the transverse
coupling rides the affine source. The per-axis streaming the row-march
feeds is the **raw down-face coefficient** :math:`g_a = |\mu_a|/\Delta a`;
since #240 Phase 2 Step D5a the diamond factor :math:`2 = 1/w_{\rm DD}`
that scales it lives in the *scheme*, not inline in the march
(:ref:`loss-rep-scanmarch-coefficient-model`). The scheme returns the
per-axis cell-balance coupling

.. math::
   :label: loss-rep-scanmarch-solve

   c_a \;=\; 2\,g_a \;=\; \frac{2|\mu_a|}{\Delta a},
   \qquad
   S \;=\; \Sigma_t + c_x + c_y,

and from it the affine x-scan coefficients
:math:`(\alpha, \beta)` that close the row,

.. math::
   :label: loss-rep-scanmarch-solve-affine

   \alpha \;=\; 2\,c_x\cdot\mathrm{inverse\_denom} - 1,
   \qquad
   \beta \;=\; \big(Q + c_y\,\psi_{y,\rm in}\big)\,
               \frac{\mathrm{inverse\_denom}}{w},
   \qquad
   \mathrm{inverse\_denom} \;=\; \frac1S,
   \quad w = \tfrac12,

where the diagonal :math:`S`, its reciprocal
:math:`\mathrm{inverse\_denom}`, the diamond transmission :math:`\alpha`,
the blend weight :math:`w`, and the transverse coupling :math:`c_y` all
come from a single scheme call,
:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.cartesian_scan_coefficients`;
the affine source :math:`\beta` is the generic
:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.source_emission`
(:math:`\beta = QV_{\rm eff}\cdot\mathrm{inverse\_denom}/w` with
:math:`QV_{\rm eff} = Q + c_y\,\psi_{y,\rm in}`). The :math:`\times
\mathrm{inverse\_denom}` reciprocal form — never a :math:`\div S`
division — is the same coefficient model the 1-D
:class:`~orpheus.sn.loss_representation.CumprodScan` rides
(:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.affine_scan_coefficients`),
so the row-march **unifies** the 1-D scan (its degenerate
:math:`c_y = 0` case) and the 2-D row-march in one primitive, is DAG-free
(no graph is built), and is now scheme-generic over every
``transverse_coupling_is_facewise`` closure rather than hard-wired to
Diamond Difference. At :math:`w = \tfrac12` these reduce to the legacy
inline DD values byte-for-byte (:math:`\alpha = 2c_x/S - 1`,
:math:`\beta = 2(Q + c_y\psi_{y,\rm in})/S`), :math:`\div\tfrac12` being an
exact power-of-2 :math:`\times 2`; the cell then closes generically in
the blend weight via
:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.cell_average`
(:math:`\bar\psi = (1-w)\psi^{\rm in}_x + w\,\psi^{\rm out}_x`) and
:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.outgoing_face_from_average`
(:math:`\psi^{\rm out}_y = (\bar\psi - (1-w)\psi^{\rm in}_y)/w`).

The matvec twin reconstructs the interior x-faces from the *known*
probe :math:`\bar\psi` through the scheme's apply-direction **reflection
scan**
(:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.reflect_scan_coefficients`):
because :math:`\bar\psi` is known, the WDD closure
:math:`\psi^{\rm out}_x = 2\bar\psi - \psi^{\rm in}_x` is itself a
first-order recurrence — a **pure-reflection scan** with the recurrence
form of
:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.outgoing_face_from_average`
at :math:`w = \tfrac12`,

.. math::
   :label: loss-rep-scanmarch-apply

   \alpha \;=\; -\frac{1-w}{w} \;=\; -1,
   \qquad
   \beta \;=\; \frac{\bar\psi}{w} \;=\; 2\bar\psi
   \qquad(\text{DD},\ w = \tfrac12),

evaluated by :func:`~orpheus.sn.spatial.scan._x_scan_faces`. The per-cell
residual then rides the uniform :math:`\div V` matvec kernel
:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.residual_kernel_batch`
that every facewise closure shares,

.. math::
   :label: loss-rep-scanmarch-apply-residual

   r_{i,j} \;=\;
   (\Sigma_t + c_x + c_y)\,\bar\psi_{i,j}
   \;-\; c_x\,\psi^{\rm in}_{x,i,j}
   \;-\; c_y\,\psi^{\rm in}_{y,i,j}
   \;\equiv\; (L+C)\bar\psi
   \quad(\text{at zero source}),

from which :eq:`loss-rep-resolution-a` subtracts :math:`\Sigma_t\bar\psi`
to give :math:`L\bar\psi`. ScanMarch additionally inherits the
conditioning robustness of ``ordinate_scan`` (ERR-054 / ERR-057 handled
per line for free) and is the natural home for the flux-independent
``a_attenuation`` two-stratum cache the wavefront lacks (the #206
follow-on).


.. _loss-rep-scanmarch-coefficient-model:

ScanMarch on the coefficient model (the #239 lift)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before #240 Phase 2 Step D5a (committed ``66dbd9a``, 2026-06-16) the 2-D
row-march **hard-coded** Diamond Difference's
:math:`\alpha`/:math:`\beta`/:math:`D` arithmetic inline:
``ScanMarch._sweep_interior`` and ``ScanMarch._loss_action_interior`` carried
the diamond factor :math:`2 = 1/w_{\rm DD}` and the cell-average blend
:math:`w = \tfrac12` literally in their bodies, so the row-march worked for
**DD only** — no other closure could ride it polymorphically. D5a folds those
two interior kernels onto the **scheme coefficient model** (#239) that the 1-D
:class:`~orpheus.sn.loss_representation.CumprodScan` has ridden since #158
Inc B. The scheme now *owns* the discretization constants; the row-march
becomes a pure schedule that asks the scheme for its coefficients and consumes
generic reconstruction ops. This is the carve that makes the scan-march
**scheme-generic over every** ``transverse_coupling_is_facewise`` closure
(:ref:`loss-rep-scanmarch-facewise`) — the equations :eq:`loss-rep-scanmarch-solve`,
:eq:`loss-rep-scanmarch-solve-affine` and :eq:`loss-rep-scanmarch-apply` above
are the post-fold (coefficient-model) form.

**The litmus that drove the fold** (``coding-elegance`` Pattern 2 / the
coefficient-model litmus): *if an explicit-matrix representation would have to
re-derive a calculation the sweep does, that calculation belongs on the
scheme.* The diamond :math:`2` and the blend :math:`w` are precisely such
shared calculations — the DD cell balance and the DD cell-average closure —
yet they sat duplicated in the row-march body, the 1-D scan cache, and the
batch kernels. D5a moves them to the one place they are defined.

**Raw streaming in, scheme-scaled coupling out.** The per-axis streaming the
row-march feeds the scheme is the **raw down-face coefficient**
:math:`g_a = |\mu_a|/\Delta a` — *not* the pre-scaled :math:`2|\mu_a|/\Delta a`
the inline code used to thread. The scheme applies its own diamond factor and
returns the scaled coupling :math:`c_a = 2 g_a` (:eq:`loss-rep-scanmarch-solve`).
Both interior kernels go through one scheme door per direction:

.. list-table:: The two scheme doors D5a routes the row-march through
   :header-rows: 1
   :widths: 18 30 52

   * - Direction
     - Scheme producer
     - Returns / consumes
   * - SOLVE
     - :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.cartesian_scan_coefficients`
       ``(s_scan, s_transverse, sig_t)``
     - the affine x-scan tuple
       :math:`(a,\ \mathrm{inverse\_denom},\ w,\ \text{transverse\_couplings})`
       with :math:`\mathrm{scan\_diag} = 2 g_x`, :math:`c_\perp = 2 g_\perp`,
       :math:`S = \Sigma_t + \mathrm{scan\_diag} + \sum_\perp c_\perp`,
       :math:`\mathrm{inverse\_denom} = 1/S`,
       :math:`a = 2\,\mathrm{scan\_diag}\cdot\mathrm{inverse\_denom} - 1`,
       :math:`w = \tfrac12`
   * - APPLY
     - :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.reflect_scan_coefficients`
       ``(psi_bar)``
     - the pure-reflection :math:`(\alpha = -1,\ \beta = 2\bar\psi)`
       recurrence (:eq:`loss-rep-scanmarch-apply`) — the :math:`(\alpha,\beta)`
       form of
       :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.outgoing_face_from_average`
       at :math:`w = \tfrac12`, :math:`w`-generic via
       :meth:`~orpheus.transport.spatial.diamond.DiamondDifference._reflection_coeffs`

**SOLVE — the row asks, the scheme answers, the cell closes generically.**
Per y-row, ``ScanMarch._sweep_interior`` calls
:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.cartesian_scan_coefficients`
to obtain :math:`(a, \mathrm{inverse\_denom}, w, (c_y,))`. The transverse-y
coupling :math:`c_y` does **not** stay separate — it folds into the affine
source: the effective volumetric source is
:math:`QV_{\rm eff} = Q + c_y\,\psi_{y,\rm in}`, and the additive scan
coefficient is the generic
:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.source_emission`,
:math:`\beta = \mathrm{source\_emission}(QV_{\rm eff}, \mathrm{inverse\_denom},
w) = QV_{\rm eff}\cdot\mathrm{inverse\_denom}/w`, while
:math:`\alpha = a`. The in-row x-face recurrence
:math:`\psi^{\rm out}_x = \alpha\,\psi^{\rm in}_x + \beta` is the **same
Blelloch scan** (:func:`~orpheus.sn.spatial.scan._x_scan_faces`) the 1-D
``CumprodScan`` runs; the cell then closes **generically in the blend weight**
:math:`w` via the base reconstruction staticmethods —
:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.cell_average`
recovers :math:`\bar\psi = (1-w)\psi^{\rm in}_x + w\,\psi^{\rm out}_x` and
:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.outgoing_face_from_average`
sheds the downstream y-face :math:`\psi^{\rm out}_y = (\bar\psi -
(1-w)\psi^{\rm in}_y)/w` (the next row's :math:`\psi_{y,\rm in}`). This entire
closure runs through
:func:`~orpheus.sn.spatial.scan._scanmarch_row`, which gained a ``w`` parameter
in D5a (it was a hard-coded ``0.5``) and now consumes the scheme-supplied
blend.

**APPLY — reconstruct off the known probe, residual through the uniform
matvec kernel.** Per y-row, ``ScanMarch._loss_action_interior`` reconstructs
the interior x-faces from the *known* probe :math:`\bar\psi` via the scheme's
reflection scan
:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.reflect_scan_coefficients`
(:eq:`loss-rep-scanmarch-apply`), then evaluates the per-cell residual **and**
the transverse-y outflow through the uniform :math:`\div V` matvec kernel
:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.residual_kernel_batch`
(:eq:`loss-rep-scanmarch-apply-residual`) — the same batched kernel the
anti-diagonal wavefront family uses for its apply direction. This kernel is the
natural home for the matvec because the apply direction has a *concrete*
:math:`\bar\psi`, so the per-cell residual is independent of its neighbours and
needs no closed-form scan; only the SOLVE direction is x-coupled (cell
:math:`i`'s inflow is cell :math:`i-1`'s outflow) and therefore needs the
closed-form affine scan coefficients.

**Why APPLY reuses a kernel but SOLVE needs a new producer.** The asymmetry is
structural, not incidental:

.. list-table::
   :header-rows: 1
   :widths: 14 43 43

   * - Direction
     - Coupling structure
     - Scheme surface
   * - APPLY
     - :math:`\bar\psi` is the **input** → per-cell residual is independent →
       a batched ÷V kernel applies directly
     - reuse the existing
       :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.residual_kernel_batch`
   * - SOLVE
     - :math:`\bar\psi` is the **unknown**, x-coupled cell-to-cell → cannot
       use an independent-batch kernel
     - new closed-form scan producer
       :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.cartesian_scan_coefficients`

**Single-sourcing the DD constants.** The fold rests on three pieces of
single-source-of-truth plumbing on
:class:`~orpheus.transport.spatial.diamond.DiamondDifference`:

* the cell-balance diagonal :math:`S = \Sigma_t + \sum_a 2 g_a` is built once
  by the shared private
  :meth:`~orpheus.transport.spatial.diamond.DiamondDifference._cartesian_streaming_diagonal`,
  consumed by **all three** Cartesian producers —
  :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.cell_kernel_batch`,
  :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.residual_kernel_batch`,
  and the new
  :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.cartesian_scan_coefficients`
  — as an **explicit left fold** :math:`((\Sigma_t + 2 g_0) + 2 g_1) + \cdots`
  (NOT ``sum()``), the IEEE-754 reduction tree of record
  (:ref:`loss-rep-bit-vs-principled`);
* the blend weight :math:`w = \tfrac12` is the named module constant ``_DD_W``
  — the literal ``0.5`` has exactly one definition site, so no body can drift
  it;
* the apply-direction reflection arithmetic
  :math:`\alpha = -(1-w)/w,\ \beta = \bar\psi/w` is the :math:`w`-generic
  :meth:`~orpheus.transport.spatial.diamond.DiamondDifference._reflection_coeffs`,
  which DD's
  :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.reflect_scan_coefficients`
  calls at ``_DD_W``, so the future Step closure inherits it for free.

The curvilinear ``affine_scan_coefficients`` merge — unifying the 1-D-chain
producer and the Cartesian row-march producer into one surface — is deliberately
**deferred to #242**: ``affine_scan_coefficients`` is 1-D-chain-shaped (it
carries curvilinear geometry arguments ``abs_mu``/``A_down``/``A_total``/``dA_w``/``c_out``/``V``
with no transverse slot), whereas ``cartesian_scan_coefficients`` needs the
transverse :math:`c_\perp` in the diagonal and the raw scan-axis :math:`g`, so
they are kept as separate methods for now (a ``d``-generic ``s_transverse``
tuple keeps the latter ready for d=3, where the row would sum multiple
transverse couplings).

**The bit-identity / principled-equivalence posture.** The carve is a
**principled ~1-ULP re-baseline** (``vv-principles``
§":ref:`Bit-identity vs principled-equivalence <loss-rep-bit-vs-principled>`",
precedent #158-B1), not a bit-identical move on the 2-D path. The reason is
specific: the pre-D5a ``ScanMarch._sweep_interior`` was the **one remaining**
place in the SN stack still using a :math:`\div D` **division**
(:math:`\alpha = 2 s_x / D`, :math:`\beta = 2(\cdots)/D`), *not* the
:math:`\times\,\mathrm{inverse\_denom}` **reciprocal** the 1-D
``CumprodScan`` already rode. Folding the 2-D solve onto the coefficient model
switches :math:`\div D \to \times(1/D)`, a genuine FP-association change — the
*same* sanctioned re-baseline the 1-D path took at #158 Inc B — and the matvec
likewise re-associates as it joins the ÷V ``residual_kernel_batch`` reduction.

.. admonition:: What re-baselined and what stayed byte-identical
   :class: note

   The fold re-associates :math:`\sim 1` ULP on the **2-D solve AND the 2-D
   matvec**; the converged *values* are unchanged. The discriminating
   evidence (verified against the three ``vv-principles`` criteria, not
   old-vs-new proximity):

   * **principled** — named coefficient-model ops, zero inline DD in
     ``_sweep_interior`` / ``_loss_action_interior`` / ``_scanmarch_row``;
   * **structurally-independent reference, two angles** — (i) the
     ``ScanMarch ≡ FullFieldWavefront`` oracle
     (``test_scan_march_equivalence.py``, :eq:`loss-rep-scanmarch`) pins the
     value to the analytical :math:`k_\infty = 1.875` /
     :math:`\varphi = Q/\Sigma_t` grounds, and (ii) ``test_keff_2d`` pins the
     homogeneous :math:`k_\infty = \nu\Sigma_f/\Sigma_a` (a closed-form,
     :math:`\ge 2`-group eigenvalue ground — never a 1-group degeneracy);
   * **two structurally-independent code paths agree** — the post-fold 2-D
     SI ``.solve`` :math:`\varphi` equals the Krylov ``.apply``
     :math:`\varphi` to :math:`\approx 2.8\times10^{-12}` relative (the SI
     iterate never touches the matvec override; the Krylov matvec does), and
     the drift is FP-non-associativity (iteration-count :math:`\times` ULP for
     SI, GMRES-tol :math:`\times` ULP for Krylov);
   * **the 1-D scan path is byte-identical** — it already rode
     :math:`\times\,\mathrm{inverse\_denom}`, so the fold did not touch it.
     This is the **negative control**: the slab SHA in
     ``test_affine_carve_bit_identity.py`` is byte-unchanged, proving the
     re-baseline is confined to the 2-D path the division-to-reciprocal switch
     actually moved.

This matches the affine-carve file's own history: the 2-D Krylov golden hashes
re-baselined at #240 Phase 2 Step B for the apply re-association; D5a extends
that to the **SI-2D solve** arm too, because the solve fold is what switched the
division to the reciprocal.

**Verification gates.** The carve is pinned by, in order:

.. list-table::
   :header-rows: 1
   :widths: 28 18 54

   * - Gate
     - Tag
     - What it pins
   * - ``test_scan_march_equivalence.py``
     - ``foundation``,
       ``verifies("loss-rep-scanmarch", -solve, -apply)``
     - ScanMarch :math:`\equiv` FullFieldWavefront oracle at ``nulp`` /
       absolute tolerance, on :math:`\ge 2`-group heterogeneous anisotropic
       **non-square** meshes (the x↔y swap moat, failure Mode 2). The existing
       parametrize already covers the Mode-9 degeneracy break (a non-square,
       multi-group, reflective tuple) — no new case was needed.
   * - ``test_streaming_operator.py::TestT4bPreT4RegressionSnapshot``
     - regression
     - the cart2d matvec frozen reference plus three slab arms; the cart2d
       ``*_apply_bulk`` snapshots were **regenerated in-commit** to the
       post-D5a value (the regenerate-in-commit discipline,
       :ref:`loss-rep-bit-vs-principled`), boundary byte-identical.
   * - ``test_affine_carve_bit_identity.py``
     - regression (strict)
     - the negative control — the 1-D slab SHA is byte-unchanged; the two
       2-D SHAs (``si_2d_p1_aniso_het``, ``krylov_2d_p1_aniso_het``)
       re-baselined off the division-to-reciprocal switch.

The two equation labels :eq:`loss-rep-scanmarch-solve` and
:eq:`loss-rep-scanmarch-apply` are the gate's ``verifies(...)`` targets, so the
rewrite above keeps both labels (the bodies changed to the coefficient-model
form; the *claims* — the scan-march coefficients and their solve/apply
realizations — are unchanged and the oracle still pins them).

**The zero-inline-DD consequence.** With the fold complete,
``ScanMarch._sweep_interior``, ``ScanMarch._loss_action_interior`` and
``_scanmarch_row`` carry **no** inline ``2.0*``, **no** ``Diamond`` reference,
and **no** hard-coded ``0.5`` in code (only in docstrings as :math:`w = \tfrac12`
notes). Polymorphic readiness is therefore satisfied *by construction*: the
moment a ``Step.cartesian_scan_coefficients`` /
``Step.reflect_scan_coefficients`` pair lands, Step rides the row-march with
**zero** further ``ScanMarch`` change. This also gives the routing gate
(:ref:`loss-rep-scanmarch-facewise`) its second, independent enforcement leg:
Linear-Discontinuous's exclusion from the row-march is now enforced both by
``ScanMarch.supports`` (the trait gate) **and** by the structural *absence* of
inline DD — the two routes agree, so a non-facewise scheme can neither be
selected for the row-march nor silently compute DD inside it.


MovingFrontierWindow — the rolling-frontier wavefront
-----------------------------------------------------

:class:`~orpheus.sn.loss_representation.MovingFrontierWindow` is the
anti-diagonal (level-scheduled) wavefront over the per-octant
:class:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph` derived at
:ref:`sweep-wavefront`. It carries only a **rolling**
:math:`(d{-}1)`-frontier of interior face fluxes (a 2-diagonal at d=2),
advanced anti-hyperplane by anti-hyperplane via
:meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed`. The
frontier is the moving realisation of the interior face cochain
:math:`C^1_{\rm int}` — its theory and the post-``WavefrontFlux``
succession are at :ref:`wavefront-flux-cochain`. Its historical claim
to fame was a **~30 % peak-memory win over the full-field oracle**; the
Fork-B2 measurement (:ref:`loss-rep-fork-b2`) showed that the *same*
memory profile is shared by the row-march, which is why the window is
now a selectable peer rather than the default.


FullFieldWavefront — the verification oracle
--------------------------------------------

:class:`~orpheus.sn.loss_representation.FullFieldWavefront` walks the
**same** anti-diagonal DAG schedule as the window, but retains the
**whole** interior face cochain — the *fuller view*. It seeds only the
octant's domain in-edge slot (``_octant_face_cochain``); by the upwind
invariant every other slot is written before any read, so the
zero-initialised buffer is byte-identical to the historical whole-trace
:math:`\iota_*` seed. It is slower and more memory-hungry by design —
its purpose is to be the reference the d-specific production
optimisations are cross-checked against:

* the 1-D :class:`~orpheus.sn.loss_representation.CumprodScan`
  (principled-equivalence at nulp);
* the 2-D :class:`~orpheus.sn.loss_representation.MovingFrontierWindow`
  (``window ≡ full`` bit-identity);
* the multi-D :class:`~orpheus.sn.loss_representation.ScanMarch`
  (principled-equivalence at nulp).

Keeping a *fuller view of the concept* as a pinned verification oracle
is the deliberate **aggressive-retirement exception**: the production
paths could not be cross-checked structurally without it. It is the one
genuinely d-generic body (``supports`` is any-d Cartesian), so it is
also the admission spine for synthetic d=3 correctness tests before any
3-axis mesh exists (the angular quadrature is already 3-cosine with all
8 sign-octants — what is missing at d=3 is ``Mesh3D``, not the
quadrature), and the representation a d=3 Cartesian mesh falls through
``default_for`` to (C3.6: :class:`ScanMarch` narrowed its ``supports``
to the d=2 truth of its kernels).


.. _loss-rep-selection:

Selection: one predicate, three consumers
==========================================

Applicability is a **declared, queryable capability** — "make illegal
states unrepresentable" applied to method selection. Each representation
answers one classmethod:

.. code-block:: python

   class CumprodScan(_LossRepresentation):
       @classmethod
       def supports(cls, mesh):
           return Compatibility(
               mesh.is_1d and mesh.scheme.is_affine_scannable,
               "requires a 1-D mesh with an affine-scannable cell-update scheme",
           )

   class ScanMarch(_LossRepresentation):
       @classmethod
       def supports(cls, mesh):
           # 1-D arm: SINGLE-axis prefix-scannability (LD's 1-D scan IS valid).
           if mesh.is_1d:
               return Compatibility(
                   mesh.scheme.is_affine_scannable,
                   "requires an affine-scannable cell-update scheme on a "
                   "1-D mesh (any geometry)",
               )
           # d≥2 arm: the DISTINCT cross-axis separability claim (DD/Step,
           # NOT LD — see the scheme-gate subsection below).
           return Compatibility(
               mesh.is_cartesian
               and mesh.ndim == 2
               and mesh.scheme.transverse_coupling_is_facewise,
               "2-D scan-march requires a scheme whose transverse coupling is "
               "facewise (separable into independent per-axis 1-D scans) — the "
               "slopeless cell-average closures (Diamond Difference, Step); "
               "Linear-Discontinuous's bilinear slope coupling needs the "
               "wavefront (the d≥3 row-march kernels are deferred — the "
               "full-field spine serves d≥3)",
           )

   class _DAGWavefront(_LossRepresentation):       # MovingFrontierWindow's base
       @classmethod
       def supports(cls, mesh):
           return Compatibility(
               mesh.is_cartesian and mesh.ndim == 2,
               "requires Cartesian geometry, d = 2",
           )

   class FullFieldWavefront(_DAGWavefront):
       @classmethod
       def supports(cls, mesh):
           return Compatibility(mesh.is_cartesian, "requires Cartesian geometry")

The compatibility signal is the *genuine* criterion — the coordinate
system (:attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.is_cartesian`, i.e.
``curvature is None``), the dimensionality
(:attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.ndim`), **and the cell-update
scheme's capability traits** — **not** the ``sweep_graphs is None``
substrate proxy that the pre-carve code keyed on. The two scan
representations read a *scheme* trait, not just geometry:
``CumprodScan`` and the 1-D arm of ``ScanMarch`` require
:attr:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.is_affine_scannable`,
and the :math:`d \ge 2` arm of ``ScanMarch`` requires the **distinct**
:attr:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.transverse_coupling_is_facewise`
— the split that closes the silent 2-D Linear-Discontinuous misroute
(:ref:`loss-rep-scanmarch-facewise`).
:class:`~orpheus.sn.loss_representation.Compatibility` is an
``(ok, reason)`` pair; the ``reason`` lets a teaching frontend gray-out
a method *and explain why* ("Moving-frontier window — requires Cartesian
geometry, d = 2"), which is pedagogically load-bearing — ORPHEUS teaches
reactor physics.

**One predicate, three consumers (single source of truth):**

#. **Frontend** —
   ``[R for R in LOSS_REPRESENTATIONS if R.supports(mesh).ok]`` lists
   the applicable methods. A cylinder (non-Cartesian) → only
   ``CumprodScan`` and ``ScanMarch``; the dropdown shows exactly those.

#. **Factory default** —
   :func:`~orpheus.sn.loss_representation.default_for` returns the
   **first** entry in the ordered registry
   :data:`~orpheus.sn.loss_representation.LOSS_REPRESENTATIONS` whose
   ``supports`` admits the mesh, falling back to the oracle so it is
   never stuck:

   .. list-table:: ``default_for`` outcomes (first-supports-match) — **facewise scheme (DD / Step)**
      :header-rows: 1
      :widths: 22 40 22

      * - mesh
        - applicable (registry order)
        - ``default_for``
      * - Cart-1D
        - ``{CumprodScan, ScanMarch, FullField}``
        - ``CumprodScan``
      * - Cart-2D
        - ``{ScanMarch, MovingFrontierWindow, FullField}``
        - ``ScanMarch``
      * - Cart-3D
        - ``{ScanMarch, FullField}`` (window is d=2 only)
        - ``ScanMarch``
      * - Cyl/Sph-1D
        - ``{CumprodScan, ScanMarch}``
        - ``CumprodScan``

   .. admonition:: The outcome depends on the **scheme**, not only the mesh
      :class: note

      The table above is for a mesh carrying a **facewise** cell-update
      scheme (Diamond Difference — the production default — or Step). Since
      #240 D5-0 the ``supports`` predicates also read the scheme's capability
      traits, so a non-facewise scheme changes which representations apply:

      .. list-table::
         :header-rows: 1
         :widths: 26 40 22

         * - mesh + scheme
           - applicable (registry order)
           - ``default_for``
         * - Cart-2D + **DD/Step** (facewise)
           - ``{ScanMarch, MovingFrontierWindow, FullField}``
           - ``ScanMarch``
         * - Cart-2D + **LD** (slope-wise)
           - ``{MovingFrontierWindow, FullField}`` — ``ScanMarch`` **refuses**
           - ``MovingFrontierWindow``

      A 2-D **Linear-Discontinuous** mesh is refused by ``ScanMarch.supports``
      (its row-march interior runs inline Diamond Difference, which would
      silently drop LD's slope) and falls through to a **wavefront**
      (``MovingFrontierWindow`` / ``FullFieldWavefront``). The wavefront's LD
      kernel is itself d=1-only today, so the 2-D LD sweep raises an honest
      ``NotImplementedError`` rather than returning DD values — see
      :ref:`loss-rep-scanmarch-facewise` for why DD/Step are facewise and LD
      is not.

   The **registry order is the policy**:

   .. code-block:: python

      LOSS_REPRESENTATIONS = (
          CumprodScan,            # 1-D production default
          ScanMarch,              # 2-D Cartesian production default
          MovingFrontierWindow,   # selectable peer (d = 2)
          FullFieldWavefront,     # never-stuck any-d oracle fallback
      )

   At d=1 ``CumprodScan`` wins (registered first; ``ScanMarch`` would
   also apply but degenerates to the same scan with a march shell). At
   d=2 ``ScanMarch`` wins — the **Fork-B2 flip** (2026-06-11) that moved
   ``ScanMarch`` ahead of ``MovingFrontierWindow`` in the tuple. At d=3
   (when ``Mesh3D`` lands) the first two refuse and ``default_for``
   falls through to the d-generic ``FullFieldWavefront`` spine — pinned
   today at the ``supports`` level
   (``test_unified_sweep_dispatch.py::TestD3SupportsMatrix``, C3.6).
   The day a d=3 row-march or window lands, widening its ``supports``
   is the *only* change needed for Cart-3D's default to upgrade — one
   predicate, no caller touched.

#. **Construction guard** — ``_LossRepresentation.__post_init__``
   re-runs ``supports(mesh)`` and raises
   :class:`~orpheus.sn.loss_representation.IncompatibleRepresentation`
   on a false verdict, so even a bypassed UI cannot build an illegal
   pairing. Combined with the frozen-dataclass immutability, the
   ``(representation, mesh)`` pairing is *correct by construction*.

That ``supports`` predicate **is** the ``is_1d`` / ``curvature``
dispatch that the pre-carve code scattered across ``transport_sweep``
plus five operator gates — now declared **once** per representation.
"Add 3-D window support" becomes "extend one representation, widen one
predicate," not a hunt through every call site.


.. _loss-rep-scanmarch-facewise:

The scan-march scheme gate: facewise vs slope-wise transverse coupling
----------------------------------------------------------------------

The :math:`d \ge 2` arm of
:meth:`~orpheus.sn.loss_representation.ScanMarch.supports` reads a **scheme**
trait,
:attr:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.transverse_coupling_is_facewise`,
not just geometry. This subsection records *why* the row-march needs a
genuinely different qualification than the 1-D scan does, and how a single
conflated predicate silently misrouted a 2-D Linear-Discontinuous mesh into
computing Diamond Difference — issue #240 Phase 2 Step D5-0, a
**routing-honesty fix** that closed a live correctness hole by making the
illegal pairing unrepresentable.

.. admonition:: The bug, in one line
   :class: important

   Before D5-0 the :math:`d \ge 2` scan-march was gated on the SINGLE-axis
   trait ``is_affine_scannable``, which Linear-Discontinuous **satisfies**.
   So ``default_for(2-D LD mesh)`` selected ``ScanMarch``, whose row-march
   interior runs *inline Diamond Difference with no scheme dispatch* — a 2-D
   LD problem **silently computed DD**, dropping LD's bilinear slope (a
   vv-principles failure Mode 2, *variable / formulation swap*, masquerading
   as a correct answer). D5-0 mints the **distinct** cross-axis trait
   ``transverse_coupling_is_facewise`` and narrows the :math:`d \ge 2` arm to
   read it, so the LD/scan-march pairing can no longer be formed.

Two different questions: 1-D prefix-scannability vs cross-axis separability
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The conflation at the root of the bug is that **one** trait was made to
answer **two** structurally distinct questions. They look similar — both are
"can a scan consume this scheme?" — but they live on different axes:

.. list-table:: The two scan-family qualifications a scheme can carry
   :header-rows: 1
   :widths: 30 34 36

   * - Trait
     - The question it answers
     - What it licenses
   * - :attr:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.is_affine_scannable`
     - **Single-axis**: is the cell-average an affine function of a *single*
       upstream face flux, :math:`\psi_{\rm out} = a\,\psi_{\rm in} + b`
       (Blelloch §1.5)?
     - The 1-D prefix scan (``CumprodScan``; the 1-D arm of ``ScanMarch``).
   * - :attr:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.transverse_coupling_is_facewise`
     - **Cross-axis**: in :math:`d \ge 2`, does a NON-swept axis couple
       through a *0th-order face value* (so the :math:`d`-D closure is
       tensor-product separable into independent per-axis scans)?
     - The :math:`d \ge 2` row-march (``scan(x) ∘ march(y)``).

The decisive fact is that **these answers can differ for the same scheme**.
A scheme can be perfectly prefix-scannable along one axis (the slope it
carries is eliminated *locally*, leaving a clean single-upstream recurrence)
yet have a multi-dimensional closure that does **not** separate across axes.
Linear-Discontinuous is exactly that scheme:
``is_affine_scannable = True`` (a 1-D fact) but
``transverse_coupling_is_facewise = False`` (a :math:`d \ge 2` fact). Diamond
Difference, by contrast, satisfies *both* — which is precisely why the bug
hid: on the production-default scheme the two traits coincide, so the
conflated predicate gave the right answer on every shipped path and only the
LD path exposed the gap.

Why the row-march needs cross-axis separability
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The row-march :math:`\mathrm{scan}(x) \circ \mathrm{march}(y)`
(:eq:`loss-rep-scanmarch`) sweeps one transverse row at a time. Within a row
it runs the *same* first-order linear scan
(:func:`~orpheus.sn.spatial.scan._scanmarch_row`) the 1-D
:class:`~orpheus.sn.loss_representation.CumprodScan` uses, and the coupling to
the previously-marched row enters **the scan's affine source** as a single
term :math:`c_y\,\psi_{y,\rm in}` (:eq:`loss-rep-scanmarch-solve-affine`). Read
off the solve coefficient
:math:`\beta = (Q + c_y\,\psi_{y,\rm in})\,\mathrm{inverse\_denom}/w`: the
*only* way axis :math:`y`
reaches the x-scan is through that one number :math:`\psi_{y,\rm in}` — the
0th-order **face value** entering the row from below. The march absorbs every
non-swept axis into the scan source as one scalar trace per cell.

This is correct **if and only if** the cell's transverse coupling really is a
single face value. Formally, the row-march is exact when the :math:`d`-D cell
closure is **tensor-product separable**:

.. math::
   :label: loss-rep-facewise-separable

   M_{d}\;=\;\bigoplus_{a}\, M^{(1)}_{a}
   \quad\Longleftrightarrow\quad
   \text{the }d\text{-D update}\ =\ \text{independent per-axis 1-D scans,
   chained by scalar face traces,}

so that the per-axis updates commute and a scan along :math:`x` marched over
:math:`y` reconstructs the same solution the joint :math:`d`-D system would.
``transverse_coupling_is_facewise`` is the name for exactly this property:
*the non-swept axis enters as a 0th-order Dirichlet face trace, not as a
higher-order moment the swept-axis row must consume.*

Diamond Difference / Step: facewise (separable)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Diamond Difference (and Step, when it is built) is a **slopeless
cell-average** closure: each cell carries a single unknown, the cell-average
flux :math:`\bar\psi`, with the downstream face reconstructed from it by the
diamond mean :math:`\psi^{\rm out} = 2\bar\psi - \psi^{\rm in}`. Its
:math:`d`-D balance (:eq:`dd-cartesian-2d`) is the standard tensor-product
central-difference structure (Lewis & Miller §§4.5, 8): the coupling from a
transverse axis :math:`y` is the **single 0th-order face value**
:math:`c_y\,\psi_{y,\rm in}` (:math:`c_y = 2 g_y`) that folds straight into the
scan source. In the
batched DD kernel
(:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.cell_kernel_batch`) this
is the explicit per-axis left fold
:math:`\mathrm{numer} = Q + \sum_a 2 g_a\,\psi^{\rm in}_a` — every axis
contributes one additive face term, exactly the separable form
:eq:`loss-rep-facewise-separable`. So a row-march along :math:`x` that carries
each :math:`y`-face value into the source is **exact** for DD/Step, and the
scheme declares ``transverse_coupling_is_facewise = True``
(:class:`~orpheus.transport.spatial.diamond.DiamondDifference`).

Linear Discontinuous: slope-wise (non-separable)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Linear Discontinuous (DG-P1-upwind) represents the in-cell angular flux as a
**linear function** rather than a constant, so each cell carries a *second*
spatial moment, the slope :math:`\hat\psi`, alongside the average
:math:`\bar\psi` (the 1-D two-moment system and its Schur-complement
reduction are derived in the
:class:`~orpheus.transport.spatial.linear_discontinuous.LinearDiscontinuous`
docstring). In **1-D** that slope is a *local* quantity:
the Schur complement of the per-cell :math:`2\times2` eliminates
:math:`\hat\psi` analytically, leaving the clean single-upstream affine
recurrence :math:`\psi_{\rm out} = a\,\psi_{\rm in} + b` — which is exactly
why LD is ``is_affine_scannable`` along one axis.

In :math:`d \ge 2` that elimination no longer decouples the axes. The cell's
transverse face flux **varies linearly along the in-cell swept coordinate**
(that is the defining feature of a P1 representation), so the coupling from a
non-swept axis is a **1st-order slope moment** :math:`\hat\psi_y`, *not* a
single face value. The slope row of the swept axis must consume it; the
:math:`d`-D closure is an irreducible, axis-coupled per-cell block whose Schur
complement does **not** diagonalize across axes — it fails the separability
:eq:`loss-rep-facewise-separable`. There is no single scalar trace the march
can fold into the scan source that would carry the transverse slope. LD
therefore declares ``transverse_coupling_is_facewise = False`` (it inherits
the conservative default;
:class:`~orpheus.transport.spatial.linear_discontinuous.LinearDiscontinuous`), and a
multi-dimensional LD problem must ride the **DAG wavefront**, which resolves
the genuine joint cell dependencies, not the scan-march.

.. note::

   **Scope.** The multi-D LD cell closure itself — the per-cell
   :math:`d`-dimensional bilinear system — is **not shipped** (it is the #240
   D5b / D6 deliverable, in design). The Cartesian multi-D LD is the
   tensor-product **bilinear (Upstream LD / UBLD)** object,
   :math:`2^d` moments on the basis :math:`\{1, x, y, xy\}` — not the
   simplex :math:`1+d` form (which fails the thick-diffusion limit on
   quadrilaterals). This page records only *that* LD's transverse coupling
   is a slope moment (hence non-separable), which is all the routing gate
   needs; the discretization of the multi-D LD block is deferred to its own
   theory section when D5b lands.

Making the illegal pairing unrepresentable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fix is a textbook application of ``coding-elegance`` Pattern 4
(*illegal states unrepresentable*): the bug was a **dimensional predicate
standing in for a scheme-capability question**. Pre-D5-0 the :math:`d \ge 2`
arm read ``is_cartesian and ndim == 2`` — a pure geometry test that admitted
*any* 2-D Cartesian scheme, including the one whose interior the row-march
cannot honour. Minting the distinct ``transverse_coupling_is_facewise`` and
narrowing the arm to read it (:meth:`~orpheus.sn.loss_representation.ScanMarch.supports`,
the refreshed code block above) means a 2-D LD mesh now answers
``supports().ok == False``, so ``default_for`` never returns ``ScanMarch`` for
it and the construction guard (``__post_init__``) rejects an explicit
``ScanMarch(2-D LD mesh)``. The pairing is gone by construction.

The trait is declared on **both** the concrete base
:class:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase` (the typed-and-
defaulted single source of truth) and the
``@runtime_checkable`` :class:`~orpheus.transport.spatial.scheme.DiscretizationScheme`
Protocol (kept symmetric with the other three capability traits). The
``@runtime_checkable`` Protocol carries a known footgun: on Python 3.12+,
``isinstance`` validates member *presence*, not type — so a scheme declaring
``transverse_coupling_is_facewise = "yes"`` would pass the ``isinstance``
check and then read *truthy* in ``supports``, re-opening a **narrower**
instance of the very silent-misroute the trait was minted to close. That
footgun is shut by ``TestCapabilityTraitsAreGenuineBools`` in
``tests/sn/sweep/core/test_discretization_scheme_protocol.py``: every
*registered* production scheme is asserted to declare all four capability
traits as a **genuine** ``bool`` (``isinstance(value, bool)`` — rejecting both
truthy ``int`` and ``np.bool_``), so a non-bool trait fails the foundation
gate rather than mis-routing in production.

.. admonition:: The default is conservative opt-in
   :class: note

   Both scan-family capability traits default to ``False`` on the base
   (``is_affine_scannable`` and ``transverse_coupling_is_facewise`` alike). A
   scheme must **declare** the capability to claim it; a newly-added scheme
   that forgets to set ``transverse_coupling_is_facewise`` is therefore
   *safely excluded* from the :math:`d \ge 2` scan-march and falls back to
   the wavefront — "slow but correct," never "fast but wrong." This is the
   exclusionary default that makes the opt-in safe.

Why scheme-named, not strategy-named
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The trait is named for the **scheme property** (separable transverse
coupling), *not* for the consuming strategy (a rejected alternative was
``is_scan_march_compatible``). Naming a scheme property after one of its
consumers is a frame-leak: it bakes the current caller into the abstraction
and forces a rename the moment a second consumer appears. And a second
consumer is already confirmed — the **diffusion ADI / line-SOR
preconditioner** (#240's next consumer) decides whether it can sweep one axis
at a time by reading *exactly this same predicate*: a line-relaxation that
marches one axis is correct only when the transverse coupling is facewise. By
naming the trait for the property, that future consumer reuses it with no
rename. The strategy-independence of the property is itself a tested fact:
``TestSchemeTraitProbe`` reads the trait directly off the scheme class with
**no** ``ScanMarch``, ``supports``, or mesh in scope, proving the answer is a
genuine scheme property a second consumer can ask of a scheme in isolation.

Verification — the routing-honesty gates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because D5-0 is a *routing* change (it alters which representation a mesh
selects, not any computed value), its gates are **selection** and
**refusal** assertions, all ``foundation``-tagged (software-structure
invariants, no theory ``:label:``) and ``-O``-safe (``pytest.fail``, never
bare ``assert`` — vv-principles failure Mode 8). Two test files carry them:

* In ``tests/sn/sweep/core/test_unified_sweep_dispatch.py``,
  ``TestD3SupportsMatrix`` pins the routing honesty four ways:

  - ``test_scan_march_refuses_2d_non_facewise_scheme_fake`` — on a synthetic
    scheme, ``supports`` **refuses** ``facewise=False`` and **admits**
    ``facewise=True`` (the anti-pattern-#11 *both-directions* pair: a
    refusal gate that only ever refuses validates the raising, not the
    invariant);
  - ``test_scan_march_refuses_2d_ld_real_mesh`` — on a real
    ``SNMesh(Mesh2D, LS-S4, scheme=LinearDiscontinuous())``,
    ``ScanMarch.supports(...).ok is False`` (the CONFIRMED-LIVE misroute,
    now closed);
  - ``test_2d_ld_default_for_routes_to_wavefront`` —
    ``default_for(2-D LD)`` lands on a wavefront
    (``MovingFrontierWindow`` / ``FullFieldWavefront``), **never**
    ``ScanMarch``;
  - ``test_2d_ld_sweep_raises_not_silently_dd`` — the **headline honesty
    claim**: ``solve_sn_fixed_source`` on a 2-D LD mesh now **raises**
    ``NotImplementedError`` (the loud d=1-only signal) rather than returning
    DD numbers via the inline-DD row-march. *Silent-wrong became
    loud-not-yet-implemented.*

* ``TestSchemeTraitProbe`` (same file) is the strategy-free probe described
  above: DD reports ``True`` standalone, LD reports ``False`` standalone,
  and ``test_facewise_distinct_from_affine_scannable`` asserts the two traits
  **diverge on LD** (``is_affine_scannable=True``,
  ``transverse_coupling_is_facewise=False``) while coinciding on DD — the
  split is *observable* only on a scheme where the two answers differ, which
  is LD. Were they to coincide on LD, the conflation that drove the misroute
  would still be latent.

* In ``tests/sn/sweep/core/test_discretization_scheme_protocol.py``,
  ``TestCapabilityTraitsAreGenuineBools`` is the genuine-``bool`` teeth that
  shut the ``@runtime_checkable`` presence-only footgun (above).

The change is **bit-identical for every exercised path**: it is a routing
predicate that touches no computed flux on any path that ran before (a 2-D LD
mesh previously *ran* — wrongly — and now *raises*; every other path is
unchanged). The strict bit-identity gate's pre-existing set was unmoved; the
only delta is the seven added D5-0 tests (closeout memo
``.claude/agent-memory/method-implementer/issue_240_phase2_step_d5_0_closeout.md``).

The honest interim state and what closes it
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

D5-0 deliberately stops at *routing*. It does **not** supply the multi-D LD
kernel — so a 2-D LD mesh, having been correctly diverted to the wavefront,
hits the wavefront's LD ``cell_kernel_batch``, which is itself d=1-only and
raises the existing ``NotImplementedError("…supports d=1 (slab/1-D) only;
got d=2…")`` from
:class:`~orpheus.transport.spatial.linear_discontinuous.LinearDiscontinuous`. This is
the **correct interim state**: a loud "not yet implemented" is strictly better
than a silent wrong answer (a different scheme's values returned under LD's
name). The follow-on:

* **D5b** supplies the multi-D bilinear LD kernel (the UBLD per-cell block),
  closing the wavefront raise — at which point 2-D LD *runs* on the wavefront,
  correctly, and ``ScanMarch`` still (correctly) refuses it.
* **D5a** folds the DD scan-march onto the coefficient model so that
  ``ScanMarch`` becomes scheme-generic — at which point LD's exclusion from the
  row-march is enforced by the *absence* of inline DD, not only by
  ``supports``; the trait gate and the structural absence then agree by two
  independent routes.


.. _loss-rep-one-walk-one-instance:

The one-walk and one-instance theorems
======================================

The structural payoff of the S6 re-layering is two theorems, each a
*type fact* enforced by construction and pinned by a discriminating spy
test (not a tautology — each test FAILED at the pre-carve HEAD and
flipped to PASS in its landing commit).

One walk (S6.4)
---------------

ONE d-generic frame —
:meth:`orpheus.sn.loss_representation._OctantWalk._interior_walk` —
serves **both** the sweep and the matvec for every multi-D
representation. For each octant it projects the quadrature octant to its
in-plane signs, branches the pure-z degenerate octants, derives the
per-axis in/out domain faces, reads the octant's inflow, runs the
interior traversal, and sheds the outflow. The two directions fork
**only** at:

* the **cell kernel** — the per-octant interior traversal the calling
  representation supplies (the window's frontier walk, the scan-march's
  row-march, the oracle's full cochain), in its solve
  (:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.cell_kernel_batch`)
  or apply
  (:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.residual_kernel_batch`)
  direction; and
* the **emit policy** — what the direction accumulates (the sweep's
  angular/moment output via ``_SweepEmit``; the matvec's
  :math:`(L+C)\psi` bulk plus the O.4b boundary defect).

.. admonition:: The direction is an OBJECT, never a boolean
   :class: warning

   The solve/apply fork is carried by **kernel and emit objects**
   (``_CellSolve`` / ``_CellResidual``,
   :class:`~orpheus.sn.loss_representation._SweepEmit`), never by an
   ``is_solve`` flag threaded down the walk. A boolean direction flag is
   the coding-elegance Smell #3 anti-pattern (special-case dispatch
   masquerading as abstraction). ``test_one_octant_walk.py`` carries an
   **AST tripwire** that parses ``_OctantWalk``'s source and fails if any
   ``is_solve`` / ``is_apply`` / ``is_matvec`` name appears — so the
   carve cannot silently degrade back into the flag pattern.

The discriminating test
``test_sweep_and_loss_action_hit_one_octant_walk`` spies the call-time
``self`` of ``_interior_walk`` and asserts that both ``L.apply`` and
``A.solve`` on a 2-D Cartesian mesh exercise it. The three-layer stack
beneath the walk (storage walk / level operation / pure kernel pair) is
documented at :ref:`sweep-dispatch-relayering`; the graph layer
(:class:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.for_shape`,
per-shape ``lru_cache`` of immutable ``MappingProxyType`` octant→DAG
maps) is family-owned, so the mesh stays pure geometry.

One instance (S6.5)
-------------------

The operator holds **one** representation instance —
:attr:`StreamingOperator.loss_representation
<orpheus.sn.operators.streaming.StreamingOperator.loss_representation>` (a
``cached_property`` = ``default_for(mesh)``) — consumed by:

* :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply` (the matvec
  :math:`(L+C)\psi`);
* :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.solve` (the forward
  substitution :math:`(L+C)^{-1}q`), via the delegating
  :attr:`InvertibleOperator.loss_representation
  <orpheus.sn.operators.streaming.InvertibleOperator.loss_representation>`
  property; and
* the boundary Gauss–Seidel resolvent.

Because representations are stateless frozen dataclasses, "one instance"
is a **structural type-fact goal** — the L21 invariant becomes
construction-enforced — not a performance fix (the per-shape DAG cache
already de-duplicated the heavy state). The discriminating test
``test_apply_and_solve_share_one_representation_instance`` spies the
call-time ``self`` of *both* doors and asserts object identity; it
**failed pre-S6.5** (the doors each called ``default_for`` independently
— two distinct frozen-dataclass ids per solve), and flips PASS once both
consume the operator's one instance.

.. note::

   **Deliberate scope boundary.** The module-level
   :func:`~orpheus.sn.loss_representation.transport_sweep` REMAINS the
   operator-free functional entry — it still selects via ``default_for``
   because its one production caller (the ``solve_sn``
   post-convergence reconstruction) has no operator in scope. The
   one-instance theorem is about the *operator's* doors; an
   operator-free functional caller legitimately mints its own.


.. _loss-rep-fork-b2:

The S6.9 measurement and the Fork-B2 decision
=============================================

The carve built ``ScanMarch`` as an **opt-in** representation (Fork B1)
with the default unchanged, precisely so the default could be decided on
*measured* numbers rather than a plausibility argument — the governing
principle: *construct each strategy as general as its algorithm
naturally allows; select narrow; specialize only on measured internal
cost.* The S6.9 benchmark
(``derivations/diagnostics/diag_s69_scanmarch_vs_window_bench.py``,
median over repeats, ``python -O``; full table in #222 comment
4683241855) measured ``ScanMarch`` / ``MovingFrontierWindow`` ratios:

.. list-table:: ScanMarch ÷ MovingFrontierWindow (lower is faster / leaner; median, ``-O``)
   :header-rows: 1
   :widths: 30 16 16 16

   * - config
     - sweep
     - matvec
     - peak mem
   * - 24×24 LS4 2g
     - 0.61
     - 0.55
     - 0.98
   * - 48×48 LS4 2g
     - 0.57
     - 0.59
     - 0.99
   * - 96×96 LS4 2g
     - 0.58
     - 0.60
     - 0.99
   * - 48×48 LS8 2g
     - 0.69
     - 0.67
     - 0.99
   * - 48×48 LS16 2g
     - 0.84
     - 0.78
     - 0.99
   * - 48×48 LS8 4g
     - 0.71
     - 0.72
     - 0.99

End-to-end fixed-source (48×48 LS8 2G heterogeneous): **10.5 s vs
12.8 s = 0.82×**.

Three findings shaped the decision:

#. **No memory edge at d=2.** The rolling frontier and the row-march are
   *both* :math:`O(\text{row})` working set at d=2 (``peak ≈ 0.98–0.99``).
   The window's memory advantage only ever held against the
   **full-field oracle** (which is ~1.3–1.4× both). The window's reason
   for being the default — peak-memory — does **not** distinguish it from
   the row-march.

#. **The ScanMarch advantage narrows with angular order** (LS4 ~0.58 →
   LS16 ~0.84). The scan's per-ordinate closed-form work scales with the
   ordinate count :math:`N`, while the level-batched wavefront amortises
   across ordinates per anti-hyperplane — so the gap closes as :math:`N`
   grows. The win is real across the measured range but is not
   asymptotically unbounded.

#. **End-to-end dilution to 0.82×.** The per-sweep kernel win is
   amortised by the scattering / moment-projection / reflect overhead
   that surrounds the within-group solve, so the wall-clock win on a full
   fixed-source solve is smaller than the bare-sweep ratio.

The decision (Fork-B2, USER, 2026-06-11):

.. admonition:: Fork-B2 decision (verbatim)
   :class: important

   *"There is no need to retire the moving frontier window method. We can
   have multiple methods (this is the whole point of having them
   selectable), as long as they are different methods, with slightly
   different principles, and the code is proper. But you can flip the
   default to ScanMarch."*

So the flip is the **one-line registry reorder** (move ``ScanMarch``
ahead of ``MovingFrontierWindow`` in
:data:`~orpheus.sn.loss_representation.LOSS_REPRESENTATIONS`);
``default_for`` is unchanged (still first-supports-match); the 1-D
default is untouched (``CumprodScan`` stays first). The window is
**kept** as a genuinely-different schedule (anti-diagonal wavefront vs
row-march) over the same lower-triangular operator — its end-to-end
coverage rides the forced-window leg of
``test_scan_march_end_to_end.py`` plus the explicit ``window ≡ full``
oracles. This is the **architecture-over-implementation** discipline
(Cardinal Rule 2): selectability is the whole point, and two proper
methods with different principles are an asset, not redundancy to be
pruned.


.. _loss-rep-bit-vs-principled:

Verification architecture: bit-identity vs principled-equivalence
=================================================================

The carve's verification turns on the vv-principles distinction between
**bit-identity** (an implementation property) and
**principled-equivalence** (a math property), applied at two different
granularities.

Across *schedules*: principled-equivalence
------------------------------------------

Different schedules (CumprodScan vs ScanMarch vs the wavefront) are
**NOT bit-comparable**: the row-march and the anti-diagonal reconstruct
the same cell dependencies *in a different order*, and IEEE-754 addition
is non-associative, so the converged values differ at FP-association.
Demanding ``array_equal`` across schedules would be the *wrong* gate. The
cross-schedule gates are therefore:

* **nulp / absolute-tolerance oracle** —
  ``test_scan_march_equivalence.py`` pins ``ScanMarch`` against the
  unconditionally-stable
  :class:`~orpheus.sn.loss_representation.FullFieldWavefront` oracle
  (G2.c), on ≥2G heterogeneous **anisotropic** configs with **non-square**
  meshes (the x↔y swap moat — failure Mode 2). The tolerance is absolute
  (``rtol=1e-11``, ``atol=1e-12``), not nulp, because near-zero boundary
  shed elements amplify a ~1e-15 absolute difference into a spurious
  ~16000-ULP reading.

* **solver-tol Mode-9 FP-invariance** —
  ``test_scan_march_end_to_end.py`` drives the full production solvers
  with the schedules swapped and asserts the **converged fixed point** is
  schedule-invariant to solver tolerance (a row-march MUST NOT move
  :math:`\psi^*` or :math:`k` — only the per-sweep FP-association, which
  the outer iteration washes out). Per vv-principles failure Mode 9, the
  FP-invariance is verified on configs that *break* the degenerate
  coincidence, never the isotropic-reflective box:

  - **G4.a** — P1-anisotropic + heterogeneous (fuel|moderator) + vacuum
    streaming, with the non-flat degenerate-gate guard;
  - **G4.b** — all-reflective + a level-symmetric cubature, so the
    outflow shed order is load-bearing (the ERR-056 shared-face failure
    class); the honest **d=2 limitation** is stated in the test — the
    full diagonal-cubature shared-face stressor is a d=3 case, deferred
    until ``Mesh3D`` exists (the quadrature is already 3-cosine; the
    schedule-level d=3 shared-face assignment is pinned synthetically by
    ``test_sweep_schedule_nd.py``, C3.6);
  - **G6** — the closed-form
    :math:`\kinf = \lambda_{\max}(A^{-1}F)` (homogeneous, ≥2G — no
    1-group eigenvalue evidence, per the cardinal-rule degeneracy) plus
    SI ≡ Krylov flux-shape agreement on the heterogeneous non-flat 2G
    config.

  The forcing is a **context manager** (a fixture would force the
  *reference* leg too and compare the window to itself), with explicit
  non-vacuity counters asserting the forced path actually ran. Since the
  Fork-B2 flip the polarity is inverted: the *default* leg runs
  ``ScanMarch`` and the *forced* leg runs ``MovingFrontierWindow`` — so
  the window peer gets its end-to-end coverage from the same gate, and
  the FP-invariance claim (being symmetric) is unchanged in meaning.

.. note::

   **MMS reaches the flux-shape layer, not the eigenvalue.** The
   ScanMarch verification matches the vv-principles claim taxonomy: the
   oracle and FP-invariance gates establish *convergence-order* and
   *flux-shape* invariance under the schedule swap; the **eigenvalue**
   anchor is the *closed-form* :math:`\kinf` leg (a structural-independence
   ground), not the SI≡Krylov twin agreement (which is necessary but not
   sufficient — two ORPHEUS paths agreeing is cross-implementation
   agreement, not correctness).

Within a *schedule*: bit-identity
---------------------------------

Two storage *policies* of the **same** schedule must agree to the byte —
they run the identical cell math in the identical level order, differing
only in how much of the cochain they retain. These ARE ``array_equal``
oracles:

* **window ≡ full** — ``test_2d_full_field_oracle.py`` (end-to-end
  sweep + matvec) and ``test_sweep_graph_window_equivalence.py`` (the
  graph-level walk, d=1 / d=2 / synthetic d=3) pin
  :class:`~orpheus.sn.loss_representation.MovingFrontierWindow` exactly
  against :class:`~orpheus.sn.loss_representation.FullFieldWavefront`.

* **the kernel pair is the FP reduction tree of record** —
  ``test_cell_kernel_batch.py::TestKernelSourceOfRecord`` freezes a
  ``sha256`` of the source of
  :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.cell_kernel_batch`
  and
  :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.residual_kernel_batch`.
  Their explicit left fold ``((Σ_t + s_0) + s_1) + …`` is
  bit-identity-load-bearing: an algebraically-equivalent rearrangement
  (a ``sum()`` instead of the fold) passes every value-tolerance test yet
  silently invalidates the 1-ULP regression contract. This is the **one
  legitimate source-hash pin** in the SN stack — the kernels *are* the
  reduction tree every byte-identity anchor inherits from.

* **the affine converged-bytes golden** —
  ``test_affine_carve_bit_identity.py`` freezes a ``sha256`` of the
  production default's converged ``angular_flux`` / ``scalar_flux``
  bytes. At the Fork-B2 flip the four 2-D hashes were **regenerated in
  the flip commit** (a schedule change shifts the converged bytes at
  FP-association level — principled-equivalent, not a numerics change),
  with a history block naming the output-identity evidence (the G4
  Mode-9 gates + the G2.c nulp oracle); the 1-D slab hashes were
  byte-unchanged (the flip's blast-radius pin). The discipline is
  **regenerate-in-commit, never pin stale**.

Structural spies
----------------

The L21 theorems are pinned by the one-walk
(``test_one_octant_walk.py``) and one-instance
(``test_one_representation_instance.py``) call-time-``self`` spies, both
``foundation``-tagged software-structure invariants (no theory
``:label:``), both ``-O``-safe (``pytest.fail``, never bare ``assert``,
so they fire under the canonical ``python -O`` invocation — vv-principles
failure Mode 8).


.. _loss-rep-history:

History and rationale: what was tried, measured, decided
========================================================

The carve arc (the WHY of the final shape)
------------------------------------------

The carve replaced a *scattered, procedural* dispatch — the same 1-D vs
multi-D decision spelled three different ways (``transport_sweep``
branching on ``reduced is not None``; the matvec branching in five
operator gates on ``not is_1d``; the oracle reachable only through
hand-built test adapters) — with one polymorphic abstraction. The
phases, each independently bit-identical-gated:

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Phase
     - What it did
   * - **S0**
     - ``test-architect`` verification plan (the proactive dispatch for
       an operator-algebra carve crossing a subsystem boundary).
   * - **S1**
     - The skeleton: protocol + ``_DAGWavefront`` base + the four
       leaves as **thin wrappers** over the existing sweep code +
       ``is_cartesian`` + ``supports`` / ``default_for`` / registry;
       ``transport_sweep`` rewired. Bit-identical.
   * - **S2**
     - The matvec side: ``strategy.residual`` collapses the **5 matvec
       gates** to one delegating call. Bit-identical.
   * - **S3**
     - Solve-vs-solve equivalence **retires the hand adapters**; the
       full-field spine becomes the genuine d-generic oracle.
   * - **S4**
     - The window generalised to ``frontier_dim = d−1`` (a point at d=1,
       a line at d=2, a surface at d=3); d=2 stays bit-identical.
   * - **S5.1**
     - ``ScanMarch`` built oracle-pinned as **opt-in** (Fork B1) — the
       default unchanged until measured.
   * - **S6**
     - The re-layering: S6.2 rename ``SweepStrategy → LossRepresentation``
       / ``residual → loss_action``; S6.3 the walk moved **off** the
       operator (``operator.py`` became pure algebra, the ``−C`` glue
       collapsed 5×→1×); S6.4 **one walk** + family-owned DAG +
       ``diamond.py`` pure kernel pair (``sweep.py`` dissolved,
       ``WavefrontFlux`` retired with the cochain succession); S6.5
       **one instance**.
   * - **S6.9**
     - Measure + the Fork-B2 default flip (this page's
       :ref:`loss-rep-fork-b2`).

Rejected alternatives (and why)
-------------------------------

* **An enum threaded into** ``transport_sweep``. Rejected: it adds a
  *second* branch axis to a function that already branches on
  dimensionality — cyclomatic complexity, not abstraction. The
  polymorphic-dispatch representation declares the choice once per type.

* **A boolean** ``is_solve`` **flag in the shared walk.** Rejected
  (coding-elegance Smell #3): the direction is carried by kernel and
  emit **objects**; an AST tripwire enforces the absence of the flag.

* **Retiring the window.** Rejected by the Fork-B2 decision: a
  genuinely-different schedule over the same operator is the *point* of
  selectability, and the window is correct and proper on its own.

The ``WavefrontFlux`` succession
--------------------------------

The S6.4 carve **retired** the typed ``WavefrontFlux`` field and its
``InteriorFaceSpace``, but the *concept* — the interior 1-cochain
:math:`C^1_{\rm int}` (with :math:`C^1 = C^1_{\rm int}\oplus C^1_\partial`
remaining valid theory) — survives in two realisations: the rolling
``_MovingFrontier`` front (the window) and the
``_octant_face_cochain`` history (the oracle). The cochain theory anchor
:ref:`wavefront-flux-cochain` is **kept**; its derivation is preserved as
the theory of both realisations.

The deferred extension point
----------------------------

The designed-but-deferred fifth representation is an **ExplicitMatrix**:
the sparse-assembled :math:`(L+C)`, whose ``sweep`` is
``scipy.sparse.linalg.spsolve_triangular``. It is the proof that the
representation abstraction is genuinely a *set of schedules over one
operator* — an assembled lower-triangular matrix is the most literal
realisation of the native frame. It slots into the registry with one
``supports`` predicate when a use case (e.g. a teaching demonstration of
the triangular structure, or a direct-solve cross-check) motivates it.


Literature
==========

The cell-level mathematics these representations *schedule* is sourced
from the standard discrete-ordinates references, anchored in the
:class:`~orpheus.transport.spatial.diamond.DiamondDifference` docstring and
:doc:`discrete_ordinates`:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Reference
     - Used for
   * - Lewis & Miller (1984), *Computational Methods of Neutron
       Transport*, §4.5 / §5.3 / §6.4
     - The Morel–Montry angular closure (§4.5), the Diamond / weighted-DD
       / Step / LD closures and the negative-flux failure mode (§5.3),
       and the canonical MMS ansatz set (§6.4).
   * - Hébert (2009), *Applied Reactor Physics*, Ch. 3 §3.9.4
     - The curvilinear S\ :sub:`N` cell-balance and DD difference
       relations.
   * - Blelloch (1990), *Prefix Sums and Their Applications*,
       CMU-CS-90-190 §1.5
     - The first-order linear-recurrence scan
       (:func:`~orpheus.sn.spatial.scan.ordinate_scan`) that
       ``CumprodScan`` and ``ScanMarch`` both reuse — the pair-monoid
       closed form.
   * - Adams & Larsen (2002), *Fast iterative methods for
       discrete-ordinates particle transport calculations*, §III
     - The within-group iteration framing in which the sweep
       :math:`(L+C)^{-1}` is the transport operator and the matvec
       :math:`(L+C)` its twin.
   * - Maginot, Ragusa & Morel (2016), *A non-negative moment-preserving
       spatial discretization scheme for solving the
       discrete-ordinates equations* (Upstream / UBLD)
     - The irreducible multi-dimensional coupling of the
       Linear-Discontinuous slope moments — the tensor-product bilinear
       (UBLD) :math:`d`-D LD closure that the routing gate
       (:ref:`loss-rep-scanmarch-facewise`) excludes from the row-march.
   * - Adams (2001), *Discontinuous finite element transport solutions in
       thick diffusive problems*, NSE 137(3)
     - The thick-diffusion-limit behaviour establishing why the multi-D LD
       closure is the bilinear :math:`\{1, x, y, xy\}` form (not the simplex
       :math:`1+d`), i.e. why the transverse coupling is a slope moment, not
       a face value.


See also
========

* :doc:`discrete_ordinates` — the cell-level WDD algebra, the 1-D
  cumprod recurrence (:ref:`sweep-cumprod`), the 2-D anti-diagonal
  wavefront (:ref:`sweep-wavefront`), and the three-layer walk /
  level-op / kernel stack (:ref:`sweep-dispatch-relayering`).
* :doc:`operator_algebra` — the Wave-O typed operator algebra
  (:eq:`operator-apply`, :eq:`operator-solve`,
  :eq:`operator-apply-transpose`) and the interior face-flux cochain
  (:ref:`wavefront-flux-cochain`).
* :mod:`orpheus.sn.loss_representation` — the module that implements
  every representation, the selection layer, and the orchestration
  entry.
* ``.claude/plans/sn_sweep_strategy.md`` — the authoritative locked
  design (decisions, verification strategy, phases S0–S6).
