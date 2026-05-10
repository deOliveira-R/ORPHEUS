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

Tensor-product primitives (Wave 0 of SN performance plan):

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
  canonical scattering / streaming form.

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


API Reference
-------------

.. automodule:: orpheus.numerics.eigenvalue
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
