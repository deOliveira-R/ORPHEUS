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
  composite, built on the shared ``AngularLift`` base below.  Since CS4c
  step 4 its Euclidean transpose is factor reversal of the cached
  conjugated product, with no fission-side :math:`w`-arithmetic
  (:ref:`sn-fission-binding-adjoint`).  Since CS4c step 5 (#450) its
  exact constructor holds the
  :class:`~orpheus.transport.material_field.FissionMaterialField`, the
  two minted :math:`L = 0` faces and its two ends — and **derives** the
  energy binding (``F.isotropic_energy``) from that one datum rather
  than being handed it, so the two bindings cannot be wired
  inconsistently.
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
  (:ref:`scattering-binding-cs4c`).  Since CS4c step 5 the
  :math:`\ell \ge 1` body is selected ONCE, at construction, from which
  end of the retained analysis face the domain's interior is: none when
  the bound datum carries no moment above :math:`\ell = 0`, the fused
  :math:`R\,\Lambda\,M` kernel on the angular end, the explicit typed
  grid path on the moment end (:ref:`cs4c-ends-select-the-body`).
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

The lift — how a binding meets its carrier
------------------------------------------

Two modules landed at CS4c step 5 (2026-09-04) carrying the machinery
that makes *"each binding acts through the body its ends select"* a
structural fact rather than a convention.  Neither is
``automodule``-rendered here (the operator convention above); the theory
is at :ref:`cs4c-ends-select-the-body`.

* ``orpheus.transport.operators.lift`` — the family's ONLY two carrier
  parses and the composite-entry verb.  ``admit_composite(op, x,
  end=...)`` admits the :class:`~orpheus.transport.full_field.FullField`
  riding the bound end's interior space; ``admit_array`` admits the bare
  :class:`numpy.ndarray` of a plain end's shape; every other operand is a
  typed refusal naming the operator, both spaces, and the way out
  (``on_moment_domain()`` for a moment iterate, ``BulkLift`` for a
  composite fed to a plain binding).  ``lift_bulk_action`` is *a bulk
  action enters the composite by extension by zero on the trace*, spelled
  once for the whole family; ``embed_bulk_assembly`` is its assembly twin
  (index-identity on the bulk block of the flat ``[bulk C-ravel | trace]``
  layout, zero trace rows and columns).  ``BulkLift`` packages both as an
  operator, and is what the 1-D diffusion solver builds over its
  plain-bound energy bindings.
* ``orpheus.transport.operators.angular_lift`` — ``AngularLift``, the
  :math:`\ell = 0` lift :math:`R_0\,E\,M_0/W` of an ENERGY binding onto
  the angular composite, and the shared base of ``TransferOperator``
  (→ :math:`S`, :math:`N_{2n}`) and ``FissionOperator`` (ruling R-1:
  ``{S, N₂ₙ} | {F}`` share the base and differ by the datum, the energy
  binding they derive, and whether an :math:`\ell \ge 1` part exists).
  Its subclass contract is three members and no verbs — ``data_ng``,
  ``_bind_energy``, ``_frame_form`` — and its public surface is
  ``apply(FullField) -> FullField``, ``apply_transpose`` likewise, plus
  ``isotropic_energy`` (the derived middle factor every scalar consumer
  reads) and ``on_moment_domain()`` (the same datum and faces re-bound to
  consume the moment composite).

The role pair — flux ↔ source/sink, declared once
-------------------------------------------------

``orpheus.transport.fields._bases`` gained two names at the same step,
and they are what let the lift verb above name a ROLE instead of a leaf:

* ``FieldRole`` — the two roles a transport field can play on one carrier
  space (``FLUX``, ``SOURCE_SINK``).  A flux is a state; a source/sink is
  a signed rate density.
* ``RolePair`` — the mixin that carries the partnership.  It is declared
  once per pair, on the source/sink leaf's class statement
  (``class AngularSourceSink(AngularField, flux=AngularFlux)``), and
  ``__init_subclass__`` registers BOTH directions into one shared
  mapping — so the map is a bijection by construction and a second
  source/sink naming an already-partnered flux is a :class:`TypeError` at
  import time.  Consumers: ``role_partner(role)`` (the leaf class playing
  a role on this carrier; asking for a leaf's own role returns that leaf),
  ``role()``, and ``into_role(role, values, space=...)`` — *"same space,
  same family fields, the other role's class"*, with the family fields
  carried across via :func:`dataclasses.fields`.

`[M]` 2026-09-04: **7** declared pairs; **16** carriers with no partner
(the 10 family/locus ABCs, the 5 residual leaves, and
``CrossSectionField``), for which asking raises a :class:`TypeError`
naming them.  The full roster and the reason the residual column is
outside the partnership are at :ref:`role-partner-declaration`.
