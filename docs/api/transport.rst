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
* :class:`~orpheus.transport.operators.fission.FissionOperator` — the
  multigroup fission source :math:`F`, realized as the rank-1 dyad
  :math:`|\chi\rangle\langle\nu\Sigma_f| = \texttt{outer}(\chi,\,
  \mathrm{ReactionRateFunctional}(\nu\Sigma_f))` whose matvec routes
  through the production-rate functional (:ref:`fission-as-dyad`).
* :class:`~orpheus.transport.operators.scattering.ScatteringOperator` —
  the scattering source :math:`S`, exposing the
  :math:`R\circ\Lambda\circ M` spectral kernel (the spherical-harmonic
  eigenbasis of the zonal scattering kernel, Funk–Hecke) plus the local
  :math:`P_0` component.  Since CS4c step 3 its exact constructor
  retains a per-material kernel FIELD, the two minted frame faces, and
  its two mandatory ends — no cross-section facade and no quadrature
  (:ref:`scattering-binding-cs4c`).
* :class:`~orpheus.transport.operators.n2n.N2NOperator` — the
  :math:`(n,2n)` source :math:`N_{2n}` on the angular composite: the
  isotropic lift of its own energy binding, extracted from :math:`S`
  at CS4c step 3 because the channel's bundling (scattering-like vs
  production-like) is context-dependent and must not be decided at the
  operator level (:ref:`sn-n2n-adjoint`).
* :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
  — the §5.7 collision operator :math:`C = M[\sigma_t]`, diagonal
  (pointwise) in every axis (:eq:`multiplication-operator-action`).
