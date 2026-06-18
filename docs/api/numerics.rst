Numerical Methods (``numerics``)
================================

The :mod:`orpheus.numerics` package holds the algorithm-agnostic
numerical primitives that every deterministic solver shares. Its
job is to keep one copy of "how to converge an eigenvalue
problem" in the codebase rather than replicating the loop in each
of the SN, CP, MOC, diffusion, and homogeneous drivers.

.. contents::
   :local:
   :depth: 2


Power Iteration
---------------

The criticality eigenvalue problem

.. math::

   A\,\phi \;=\; \frac{1}{k}\,F\,\phi

has a spectrum of eigenvalues
:math:`k_0 > k_1 > k_2 > \dots`. Only the **dominant eigenvalue**
:math:`k_0 = k_{\rm eff}` and its eigenvector :math:`\phi_0` are
physically meaningful: by the Perron–Frobenius theorem
:math:`\phi_0` is the unique non-negative eigenvector, while all
higher harmonics change sign in space.

:func:`~orpheus.numerics.eigenvalue.power_iteration` converges to
:math:`(k_0, \phi_0)` by repeatedly applying the transport
operator to an estimate of :math:`\phi`:

.. math::

   \phi^{(n+1)} \;=\; A^{-1}\,\frac{1}{k^{(n)}}\,F\,\phi^{(n)},
   \qquad
   k^{(n+1)} \;=\; \frac{\lVert F\,\phi^{(n+1)}\rVert}
                         {\lVert L\,\phi^{(n+1)}\rVert}.

The convergence rate is governed by the **dominance ratio**
:math:`|k_1/k_0|`; problems with a narrow spectral gap (large
lattices, near-critical systems with weakly coupled regions)
converge slowly and may benefit from Chebyshev or Wielandt
acceleration — not currently implemented in ORPHEUS.

**Normalisation.**
The returned eigenvector has *arbitrary* absolute scale. Power
iteration preserves shape but not magnitude — callers that need
absolute flux (e.g. for power calibration, dose calculations) must
post-normalise, typically by fixing the total integral fission
source or the total power deposition.


The EigenvalueSolver Protocol
-----------------------------

Every deterministic solver plugs into the same power iteration
loop by implementing the
:class:`~orpheus.numerics.eigenvalue.EigenvalueSolver` protocol.
The protocol has five methods and one structural contract:

* ``initial_flux_distribution`` — return a flux guess. Most
  solvers use a flat unit array; MOC uses a cell-averaged flat
  angular flux.
* ``compute_fission_source`` — build
  :math:`Q_f = \chi\,(\nu\Sigma_f\,\phi)/k`. Pure function of the
  current flux and eigenvalue.
* ``solve_fixed_source`` — apply :math:`A^{-1}` to the fission
  source. **Scattering and (n,2n) sources are assembled *inside*
  this method** because they need to be updated between inner
  iterations (source iteration in SN, Gauss–Seidel in CP, etc.).
  This is the single most important structural decision in the
  protocol: it lets each solver manage its own inner iteration
  strategy without leaking through to the outer loop.
* ``compute_keff`` — update the eigenvalue from the current
  :math:`\phi`. For reflective lattices the leakage term is zero;
  for whole-core diffusion it is not.
* ``converged`` — stopping test. Typical tolerance
  :math:`10^{-6}` on :math:`|\Delta k|`; richer tests on flux
  L2 norm are also used.

**Reference implementations** (each satisfies the protocol and is
tested against the power-iteration loop without any solver-specific
glue):

* :class:`orpheus.cp.solver.CPSolver` — collision probability.
* :class:`orpheus.sn.solver.SNSolver` — discrete ordinates.
* :class:`orpheus.moc.solver.MOCSolver` — method of characteristics.
* :class:`orpheus.diffusion.solver.DiffusionSolver` — 1-D two-group
  diffusion.
* :class:`orpheus.homogeneous.solver.HomogeneousSolver` — infinite
  homogeneous medium.


Operator Algebra (Wave A)
-------------------------

The :mod:`orpheus.numerics.operator` module installs the matrix-free
operator-algebra primitives consumed by every solver. See
:ref:`operator-algebra` for the design rationale, capability-set
semantics, and tensor-product algebra.

Tensor-product primitives (Wave 0 of SN performance plan, Wave T
consumers landed May 2026):

* :class:`~orpheus.numerics.operator.DiagonalOperator` — diagonal
  multiplication on a tagged tensor axis. The ``AngularWeightMatrix``
  :math:`W` of Grand Report v3 §9 is
  ``DiagonalOperator.from_measure(quad.measure, axis=0)``.
* :class:`~orpheus.numerics.operator.TensorProductOperator` —
  per-axis tensor product :math:`A \otimes B \otimes \cdots`. Built
  via the ``&`` dunder; carries axis tags and the closure laws
  :math:`(A \otimes B)^* = A^* \otimes B^*`,
  :math:`(A \otimes B) \circ (C \otimes D) = (A \circ C) \otimes
  (B \circ D)`. See :ref:`tensorial-framing`.
* :class:`~orpheus.numerics.operator.SumOfTensorProductsOperator`
  — :math:`\sum_k A_k \otimes B_k \otimes \cdots`; the §15.2
  *aspirational* canonical scattering / streaming form. **Zero
  production consumers** today; see :ref:`wave-t-tensor-network`
  for why (the MA-Q1 master condition: coupled physics falls back to
  :class:`OperatorSum` over bespoke leaves, not SOTP).
* :class:`~orpheus.numerics.operator.OperatorSum` — the additive
  composer :math:`A + B`. Wave T promoted this to load-bearing
  status (T.3 scattering kernel; the T.4 per-direction streaming split
  was retired in #238 — the fused matvec walks both directions in one
  pass); see :ref:`wave-t-shape-table`.
* :class:`~orpheus.numerics.operator.RankOneOperator` — rank-1
  outer product :math:`|\ell\rangle\langle r|` on a tagged axis;
  native to the multigroup fission emission
  :math:`F = \chi \otimes \nu\Sigma_f` (Wave T T.2,
  :attr:`FissionOperator.kernel`).


Consumer matrix for tensor-product primitives (post Wave T)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The table below lists the production consumers as of Wave T close-out
(May 2026). The full architectural narrative — including the MA-Q1
master condition that decides between :class:`TensorProductOperator`,
:class:`SumOfTensorProductsOperator`, and
:class:`OperatorSum`-of-bespoke-leaves — is at
:ref:`wave-t-tensor-network`.

.. list-table:: Production consumers of the tensor-product primitives
   :header-rows: 1
   :widths: 28 14 28 30

   * - Primitive
     - Consumers (count)
     - Examples
     - Notes
   * - :class:`~orpheus.numerics.operator.TensorProductOperator`
     - 6
     - 5 BC realizers (vacuum / specular / white / albedo / periodic
       via Wave T T.1 ``& IdentityOperator()`` wrap); fission kernel
       (Wave T T.2,
       :attr:`FissionOperator.kernel = RankOneOperator(χ, νΣ_f,
       axis=0) & IdentityOperator()`)
     - Six clean-TP production instances. The MA-Q1 master condition
       is satisfied: each consumer factors as disjoint per-axis
       operations.
   * - :class:`~orpheus.numerics.operator.SumOfTensorProductsOperator`
     - 0
     - (no production consumers)
     - The §15.2 SOTP form is aspirational; Wave T T.3 (scattering)
       and T.4 (streaming) both fell back to :class:`OperatorSum`
       over bespoke leaves per the MA-Q1 master condition. See
       :ref:`wave-t-tensor-network` for the per-substep rationale.
   * - :class:`~orpheus.numerics.operator.OperatorSum`
     - 1 (load-bearing, post Wave T)
     - Wave T T.3 scattering kernel
       (:attr:`ScatteringOperator.kernel` =
       ``reduce(add, kernel_summands)`` over per-ℓ
       :class:`_PerLegendreOrderScattering`)
     - The T.4 streaming per-direction split (``M_spatial`` as an
       :class:`OperatorSum` of two per-direction summands) was retired
       in #238 — it had no production consumer; the fused matvec
       (:meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk`)
       walks both directions in ONE bidirectional pass.
       See :ref:`wave-t-orchestrated-apply`.
   * - :class:`~orpheus.numerics.operator.RankOneOperator`
     - 1
     - Wave T T.2 fission kernel
       (:attr:`FissionOperator.kernel` first factor of the
       :class:`TensorProductOperator`)
     - Encodes the group-axis contraction-then-broadcast
       :math:`(F\,\phi)_g = \chi_g\,\sum_{g'}\nu\Sigma_{f,g'}\,\phi_{g'}`
       as a typed primitive.

**New SN-side bespoke leaves shipped in Wave T** (private to
:mod:`orpheus.sn`, accessed via public properties; documented for
forward-reference by Wave O typing work):

* :class:`orpheus.sn.scattering._PerLegendreOrderScattering` — per-ℓ
  summand of :attr:`ScatteringOperator.kernel`. Public access via
  :attr:`ScatteringOperator.kernel_summands`. Implements
  :math:`R_\ell \circ \Lambda_\ell \circ M_\ell` for one Legendre
  order.

.. note::

   Wave T also shipped a per-direction streaming split (``M_spatial`` /
   ``M_angular_redist`` over ``_SpatialSweepDirection`` /
   ``_MSpatialOperatorSum`` / ``AngularRedistributionOperator``) as
   separately-applicable typed leaves to anticipate Wave-O adjoint /
   DSA consumers. **#238 retired that split**: no production code ever
   applied the leaves separately (the #240 adjoint uses the fused
   ``loss_action_transpose``), so the streaming + curvilinear
   Morel–Montry angular redistribution is computed IN-SWEEP inside the
   fused matvec
   (:meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk`),
   verified end-to-end by the anisotropic curvilinear MMS
   (:ref:`sn-mms-curvilinear-aniso-verification`).

The leading-underscore primitives are intentionally private (the
public surface is via the :attr:`~orpheus.sn.scattering.ScatteringOperator.kernel`
/ :attr:`~orpheus.sn.scattering.ScatteringOperator.kernel_summands`
properties on the operator classes). Wave O (`Issue #208
<https://github.com/deOliveira-R/ORPHEUS/issues/208>`_) will introduce
``BulkOperator`` / ``FullOperator`` / ``BoundaryOperator`` Protocols
that may promote some of these to public status if a downstream
consumer surfaces.

Discrete measures and partition (Wave 0 of SN performance plan):

* :class:`~orpheus.numerics.measure.DiscreteMeasure` — atomic
  measure :math:`\mu = \sum_i w_i\,\delta_{x_i}` on a measurable
  space; carries integration, tensor product, direct sum,
  pushforward, restriction, and the new
  :meth:`~orpheus.numerics.measure.DiscreteMeasure.partition_by`
  primitive (the inverse of direct sum). See
  :ref:`discrete-measure-partition`.
* :class:`~orpheus.numerics.measure.DiscreteMeasurePartition` —
  a partition entry returned by
  :meth:`partition_by`. Carries label, indices into the parent,
  and the restricted measure.

Galerkin / Petrov-Galerkin projection (Wave 0 of SN performance plan):

* :class:`~orpheus.numerics.projection.ProjectionOperator` —
  most-general projection ABC.
* :class:`~orpheus.numerics.projection.GalerkinProjection` —
  Galerkin discipline (test space = trial space) ABC.
* :class:`~orpheus.numerics.projection.PetrovGalerkinProjection`
  — Petrov-Galerkin discipline (test space ≠ trial space)
  ABC. Sibling of Galerkin; concrete subclasses land with energy
  condensation (§17).
* :class:`~orpheus.numerics.projection.HarmonicMomentProjection`
  — concrete Galerkin projection on real spherical harmonics.
* :class:`~orpheus.numerics.projection.HarmonicMomentReconstruction`
  — paired reconstruction with the addition-theorem
  :math:`(2\ell+1)` factor.

See :ref:`galerkin-projection` for the discipline narrative,
cross-method consumer table, and the naming-hierarchy rationale.

Real spherical harmonics:

* :func:`~orpheus.numerics.spherical_harmonics.evaluate_real_sh`
  — :math:`Y_\ell^m(\hat\Omega_n)` evaluator; uses the
  no-:math:`4\pi/(2\ell+1)`-prefactor convention so the addition
  theorem reads :math:`\sum_m Y_\ell^m Y_\ell^m = P_\ell(\Omega
  \cdot \Omega')`. See :ref:`spherical-harmonics`.


Field algebra (Depth B, step D-A)
---------------------------------

.. _field-algebra:

The :class:`~orpheus.numerics.field.Field` ABC is the L1 algebraic
base of every typed transport field — angular flux, scalar flux,
spherical-harmonic moments, boundary traces, sources, residuals. It
codifies the Grand Report v3 §5.5 / §32.5 prescription:

   *Every typed transport field is the pair ``(values, space)`` with
   closed same-CLASS, same-SPACE arithmetic.*

The ABC carries the dunder algebra (``+``, ``-``, ``-`` unary, scalar
``*``, scalar ``/``) and the diagnostics (``linf``, ``l2``,
``inner_product``) inherited unchanged by every concrete subclass.
Subclasses add domain-specific fields (``mesh``, ``boundary``,
``history``) on top of the ``(values, space)`` base; the algebra is
inherited verbatim via :func:`dataclasses.replace`. The same
hand-coded dunder skeleton that previously lived in six separate
classes (``AngularFlux``, ``ScalarFlux``, ``HarmonicMomentField``,
``BoundaryFlux``, ``IsotropicSource``, ``PerOrdinateSource``) is
consolidated here — Cardinal Rule 2 (single source of truth).

**Three-layer dimensional enforcement.** Dimensional consistency is
gated at three layers, each with a different cost / coverage trade-off:

* **Layer 1 — class identity.**
  :meth:`~orpheus.numerics.field.Field._check_partner` rejects
  ``type(self) is not type(other)`` before any value comparison. This
  is the *primary* gate: even when units match (an ``AngularSourceSink``
  and an ``AngularResidual`` may both carry
  :math:`1/(\mathrm{cm^2 \cdot s \cdot sr \cdot eV})`),
  cross-class arithmetic raises by construction. Same units gives
  PERMISSION to add in linear algebra; it does not give MEANING.
* **Layer 2 — construction-time dimensional check.** Solvers like
  :class:`~orpheus.numerics.iteration.SourceIteration` do a single
  :math:`O(1)` ``pint.Unit.dimensionality`` comparison per operator at
  ``__init__`` to verify the operator algebra is dimensionally sound
  before any iteration runs. Cost: microseconds per build. ALWAYS
  runs (both default ``python -O -m pytest`` and ``pytest``).
* **Layer 3 — defensive assert in dunders.** Inside
  ``_check_partner``, ``assert self.space.units == other.space.units``
  catches the rare class/units misdesign (two instances of the same
  class whose spaces nonetheless carry inconsistent unit STRINGS — e.g.
  one in :math:`1/\mathrm{cm}^2`, one in :math:`1/\mathrm{m}^2` — same
  dimensionality, different scaling). Stripped in ``-O`` mode; defense
  in depth during development.

Together these layers make dimensional-mismatch bugs unrepresentable
without paying the cost of full ``pint.Quantity`` arithmetic on every
ndarray operation.

**Class identity for cross-class same-units operations.** Two mechanisms
work together when distinct Field subclasses share a dimensional
signature (e.g. ``AngularResidual`` and ``AngularSourceSink`` both carry
``1/(cm³·s·sr)``):

* **The class gate forbids cross-class arithmetic** even when units
  match — ``angular_residual - angular_source`` RAISES. Same units grant
  permission to add in linear algebra, not meaning.
* **A role transition goes through a named composition.** The transport
  balance :math:`r = A\psi - q` differences two *same-class*
  ``AngularSourceSink`` operands (the operator output :math:`A\psi` and
  the external source :math:`q`), but the *result* is a residual (a
  defect), not a source. Bare same-class subtraction typechecks yet
  MIS-TYPES the defect; each residual leaf's ``from_balance`` factory
  makes the transition explicit and lands the result in the correct
  residual class:

  .. code-block:: python

      # same-class subtraction typechecks but MIS-TYPES the defect:
      wrong = operator_output - q_ext           # AngularSourceSink (!)

      # named composition — correctly typed as the residual role:
      residual = AngularResidual.from_balance(lhs=operator_output, rhs=q_ext)

The named-composition discipline IS what makes the Field algebra
sound under physical interpretation — ``coding-elegance`` Pattern 4
(illegal states unrepresentable) takes its strictest form here.

See the Depth B plan
(``.claude/plans/depth_b_field_on_function_space.md``) §3.2 for the
ABC spec, §3.7 for the singledispatch policy that consumes ``Field``
in operator apply, §5 for the Layer 2 construction-time check, and
§7.5 for the full CI matrix.


API Reference
-------------

.. automodule:: orpheus.numerics.eigenvalue
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: orpheus.numerics.field
   :members:
   :undoc-members:
   :show-inheritance:

The :class:`~orpheus.numerics.operator.LinearOperator` Protocol,
its capability-set semantics, and the composition / tensor-product
primitives are documented at :ref:`operator-algebra` (theory page).
The Galerkin / Petrov-Galerkin projection ABCs and the concrete
:class:`HarmonicMomentProjection` / :class:`HarmonicMomentReconstruction`
pair are documented at :ref:`galerkin-projection`. The
:math:`Y_\ell^m` evaluator and the no-:math:`4\pi/(2\ell+1)`-prefactor
convention are documented at :ref:`spherical-harmonics`. The
:meth:`partition_by` primitive on
:class:`~orpheus.numerics.measure.DiscreteMeasure` is documented at
:ref:`discrete-measure-partition`. The theory pages contain the full
mathematical narrative; per-symbol API docstrings live in the
modules themselves and are accessible via the standard
``orpheus.numerics`` import path.

Function spaces — trace and interior face cochains
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :class:`~orpheus.numerics.space.FunctionSpace` ABC carries the
``(name, shape, inner_product_weights)`` identity for every typed
field. Two **face** spaces specialise it for the SN sweep's
codim-1 (face) quantities, both re-homing their
:class:`~orpheus.numerics.face_layout.FaceLayout` onto the space as
``compare=False`` leaf-data (the A.5 "layout-on-space" pattern):

* :class:`~orpheus.numerics.spaces.trace_space.TraceSpace` — the
  **boundary** trace cochain :math:`C^1_\partial`, carrying
  ``omega_dot_n`` for the inflow / outflow directional partition. The
  space of :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`.
* the **interior** face cochain :math:`C^1_{\rm int}` — through S6.4
  the typed ``InteriorFaceSpace`` (the trace space *minus*
  ``omega_dot_n``; flux-only, axis-parametric) carrying
  ``WavefrontFlux``; RETIRED at S6.4(f) when the walk re-layering
  dissolved the type's boundary algebra into the shared octant frame.
  The cochain concept's live realizations are the rolling front
  (``_MovingFrontier``) and the per-octant full-cochain buffers
  (``FullFieldWavefront._octant_face_cochain``).

Together the boundary and interior cochains biproduct-decompose the
full face cochain
:math:`C^1 = C^1_{\rm int} \oplus C^1_\partial` (Issue #205 Phase 5).
The full cochain frame, the :math:`\iota_*` / :math:`\iota^*` trace
operators, the flux-only-single-role rationale, and the succession
history are documented at :ref:`wavefront-flux-cochain` (theory page).
