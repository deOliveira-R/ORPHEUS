Transport Operator Algebra (``transport``)
==========================================

Reference for the method-agnostic **transport operator-algebra** layer —
the L2 operators and functionals shared across the SN / CP / MoC solvers
(:mod:`orpheus.transport.operators` and the reaction-rate functionals).
These are the concrete :class:`~orpheus.numerics.operator.LinearOperator`
and :class:`~orpheus.numerics.functional.Functional` leaves of the §5.6
suffix law; the full mathematical narrative — the integral-kernel /
multiplication / functional partition, the fission rank-1 dyad, and the
scattering spectral sum — lives at :ref:`operator-algebra` (theory page).

The TransportMethod Protocol
----------------------------

The structural Protocol over the **method-mesh layer** (#290 P7b):
``SNMesh`` and ``DiffusionMesh`` conform without importing it, and the
ONE shared ``resolve_boundary_conditions`` body turns each mesh's
per-axis :class:`~orpheus.geometry.mesh.BC` declarations into realized
boundary operators through the per-method ``realize_boundary_law``
hook. The module docstring carries the full design record (the
two-witness genesis, the instance-surface-only ruling, and the
realizer-registry dissolution rationale); the boundary-architecture
narrative lives at :doc:`/theory/foundations/boundary_conditions`.

.. automodule:: orpheus.transport.method
   :members:
   :undoc-members:
   :show-inheritance:

The reaction-rate functionals
-----------------------------

The §5.6 reaction-rate co-vector :math:`\langle\Sigma_x,\cdot\rangle` and
its volume-integrated companion. ``ReactionRateFunctional`` is the
per-cell density (fission's rank-1 row-factor; see :ref:`fission-as-dyad`);
``IntegratedReactionRate`` is the scalar total behind the SN
:math:`k`-eigenvalue (:ref:`integrated-reaction-rate-keff`). Both
specialise the generic numerics
:class:`~orpheus.numerics.functional.InnerProductFunctional`.

.. automodule:: orpheus.transport.reaction_rate_functional
   :members:
   :undoc-members:
   :show-inheritance:

The integral-kernel operators
-----------------------------

The §5.6 **Kernel** refinement and the two emission kernels carry rich
``.. math::`` docstrings written for in-prose cross-reference; their
mathematical narrative lives in the theory page rather than being
``automodule``-rendered here (matching the
:mod:`orpheus.numerics.operator` convention). Per-symbol API docstrings
are accessible via the standard import path.

* :class:`~orpheus.transport.operators.integral_kernel_operator.IntegralKernelOperator`
  — the §5.6 Kernel Protocol (a refinement of
  :class:`~orpheus.numerics.operator.LinearOperator` adding a ``kernel``
  member). See :ref:`integral-kernel-category`.
* :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicFission`
  — the multigroup fission source :math:`F` in its **energy** binding,
  realized as the rank-1 dyad
  :math:`|\chi\rangle\langle\nu\Sigma_f| = \texttt{outer}(\chi,\,
  \mathrm{ReactionRateFunctional}(\nu\Sigma_f))` whose matvec routes
  through the production-rate functional (:ref:`fission-as-dyad`).  This
  is the class the S\ :sub:`N` k-outer, the 1-D diffusion solver and the
  infinite-medium solver all consume.
* :class:`~orpheus.transport.operators.fission.FissionOperator` — the
  same datum's **angular** binding: the harmonic frame's
  :math:`\ell = 0` conjugation of that dyad on a posed angular
  composite, retaining the energy binding as its middle factor exactly
  as ``N2NOperator`` retains ``IsotropicN2N``.  Since CS4c step 4 its
  Euclidean transpose is factor reversal of the cached conjugated
  product, with no fission-side :math:`w`-arithmetic
  (:ref:`sn-fission-binding-adjoint`).
* :class:`~orpheus.transport.operators.transfer.TransferOperator` — the
  **collision-gain core**, and since #426 step 2 (2026-09-04) the one
  home of both gains' arithmetic: the angular binding
  :math:`T_c = \tfrac1W R\,\Lambda_c\,M` of a channel's Legendre
  transfer stack with its yield, exposing the
  :math:`R\circ\Lambda\circ M` spectral kernel (the spherical-harmonic
  eigenbasis of the zonal kernel, Funk–Hecke) plus the local :math:`P_0`
  fast path.  Since CS4c step 3 its exact constructor retains a
  per-material kernel FIELD, the two minted frame faces, and its two
  mandatory ends — no cross-section facade and no quadrature
  (:ref:`scattering-binding-cs4c`).
* :class:`~orpheus.transport.operators.scattering.ScatteringOperator` and
  :class:`~orpheus.transport.operators.n2n.N2NOperator` — the two
  **terms** of the within-group algebra
  :math:`A = L + C - S - N_{2n} - B`, as thin ROLE subclasses of that
  core.  ⭐ **A role is two class constants and no code**: ``channel``
  (a ``ClassVar`` holding the
  :class:`~orpheus.transport.material_field.TransferMaterialField`
  constructor for its channel) and ``isotropic_binding`` (the P0 energy
  binding it lifts).  An AST gate
  (``tests/transport/test_transfer_roles.py``) refuses **any** method on
  a role — no ``apply``, no ``apply_transpose``, no classmethod, no
  field — so the twin path the collapse removed cannot regrow one
  override at a time.  The tier-2 mint
  (:meth:`~orpheus.transport.operators.transfer.TransferOperator.from_solver_data`)
  and the exact constructor are both the CORE's, and read the role's
  two constants.  They differ in the **yield** :math:`y` alone —
  :math:`y_S = 1`, :math:`y_{2n} = 2`
  (:data:`~orpheus.transport.kernels.N2N_MULTIPLICITY`) — and
  ``N2N.apply(ψ) == 2·S'.apply(ψ)`` over the same stack is bit-exact.
  :math:`N_{2n}` was extracted from :math:`S` at CS4c step 3 because
  the channel's bundling (scattering-like vs production-like) is
  context-dependent and must not be decided at the operator level; the
  kernel tier names the mathematical OBJECT and the operator tier names
  the TERM (:ref:`sn-n2n-adjoint`).

  ⛔ Until 2026-09-04 :math:`N_{2n}` was described here as *"the
  isotropic lift of its own energy binding"*.  That was a truthful
  description of a :math:`P_0` MODEL imposed at this tier, and the
  model was a defect (ERR-082): the tape stores seven Legendre moments
  for MT=16, the same order as elastic.
* :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
  — the §5.7 collision operator :math:`C = M[\sigma_t]`, diagonal
  (pointwise) in every axis (:eq:`multiplication-operator-action`).
