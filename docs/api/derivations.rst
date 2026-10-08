Analytical Derivations (``derivations``)
========================================

Reference for the :mod:`orpheus.derivations` package — the single
source of truth for analytical reference eigenvalues used across the
V&V ladder. Each derivation module is a SymPy- or closed-form
calculation that emits :class:`~orpheus.derivations.common.verification_case.VerificationCase`
objects carrying the analytical :math:`k_\infty` (or :math:`k_\text{eff}`),
the material definitions, the geometry parameters, a LaTeX trace of the
derivation, and the V&V level.

Tests pull these cases through
:func:`orpheus.derivations.reference_values.get` (exposed at the
package top level as ``orpheus.derivations.get``), and the conftest hook
propagates ``vv_level`` / ``equation_labels`` from the case onto the
parametrized test node so that a single ``@pytest.mark.verifies(...)``
covers every consumer automatically.

Submodules
----------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Module
     - Purpose
   * - :mod:`~orpheus.derivations.continuous.analytical.homogeneous`
     - Infinite-medium :math:`k_\infty` for 1/2/4-group synthetic XS.
   * - :mod:`~orpheus.derivations.continuous.flat_source_cp.slab`
     - Closed-form slab collision-probability eigenvalues.
   * - :mod:`~orpheus.derivations.continuous.flat_source_cp.cylinder`
     - Closed-form cylindrical-annulus collision-probability eigenvalues.
   * - :mod:`~orpheus.derivations.continuous.cases.diffusion`
     - 1-D / 2-region two-group diffusion eigenvalues matching the
       ``diffusion`` solver's core geometry.
   * - :mod:`~orpheus.derivations.common.kernels`
     - Exponential-integral (``E_3``) and Bickley–Naylor
       (``Ki_3``/``Ki_4``) kernels shared by the slab and cylinder
       derivations, plus the :func:`chord_half_lengths` chord-segment
       primitive consumed by chord-impact-parameter integrals.
   * - :mod:`~orpheus.derivations.common.quadrature`
     - Unified 1-D quadrature contract
       (:class:`~orpheus.derivations.common.quadrature.Quadrature1D` value
       object) and its primitive constructors:
       :func:`gauss_legendre`,
       :func:`gauss_legendre_visibility_cone`,
       :func:`composite_gauss_legendre`, :func:`gauss_laguerre`.
       Plus the sibling
       :class:`~orpheus.derivations.common.quadrature.AdaptiveQuadrature1D`
       (no-fixed-nodes adaptive rule built via
       :func:`~orpheus.derivations.common.quadrature.adaptive_mpmath`).
   * - :mod:`~orpheus.derivations.common.dense_pencil`
     - The reference kernel's dense linear algebra:
       :class:`~orpheus.derivations.common.dense_pencil.DensePencil`, the
       pencil :math:`F c = k\,L c` in weak form, with its fundamental mode
       (refused, :class:`~orpheus.derivations.common.dense_pencil.NoFundamentalMode`,
       unless real, positive, strictly dominant and single-signed), its
       full spectrum, its adjoint (the transposed pair) and the least
       solution of :math:`L x = F x + s` on the unknowns the source
       reaches (refused,
       :class:`~orpheus.derivations.common.dense_pencil.NoLeastSolution`,
       unless that reached part is subcritical; zero elsewhere). The references' own copy of what production
       spells in :mod:`orpheus.numerics.eigenvalue`, kept apart on purpose
       (:ref:`architecture-reference-insulation`); the theory is
       :ref:`verification-reference-kernel`.
   * - :mod:`~orpheus.derivations.continuous.characteristic`
     - The characteristic reference, built rung by rung: the walls of a
       concentric body read from their laws' factors, with the reference's
       own tag registry
       (:mod:`~orpheus.derivations.continuous.characteristic.walls`), and
       the line part of the boundary resolvent, the period of each line's
       unfolded path and the least solution of its cycle
       (:mod:`~orpheus.derivations.continuous.characteristic.closure`), the
       panel basis the emission density and the flux are represented in
       (:mod:`~orpheus.derivations.continuous.characteristic.basis`), and
       the transport along each line on that basis: the traversal
       integrals, the vacuum Volterra block and the angular flux
       (:mod:`~orpheus.derivations.continuous.characteristic.transport`),
       one group's transport block, assembled by Galerkin over the lines of
       the chart's line domain with the white walls' coupling
       (:mod:`~orpheus.derivations.continuous.characteristic.assembly`),
       each region's total cross section and emission matrices with the
       regions each group emits in
       (:mod:`~orpheus.derivations.continuous.characteristic.cross_sections`),
       the multigroup Galerkin system on the emission space, with its k
       pencil, its source pencil and their adjoints
       (:mod:`~orpheus.derivations.continuous.characteristic.system`), and
       the geometric, exponential and hp gradings every rule of the package
       places its piece ends by
       (:mod:`~orpheus.derivations.continuous.characteristic.grading`);
       the theory is :ref:`theory-characteristic-reference`.
   * - :mod:`~orpheus.derivations.common.quadrature_recipes`
     - Geometry-aware quadrature recipes:
       :func:`chord_quadrature` (impact-parameter integrals on
       concentric annular geometries) and
       :func:`observer_angular_quadrature` (kink-aware ω-sweeps from
       an internal observer).
   * - :mod:`~orpheus.derivations.common.xs_library`
     - Synthetic cross-section library (``_FUEL_XS``, ``_REFL_XS``, …)
       that guarantees derivation cases and solver tests use the exact
       same numbers.
   * - :mod:`~orpheus.derivations.common.verification_case`
     - :class:`VerificationCase` dataclass + ``VVLevel`` literal.
   * - :mod:`~orpheus.derivations.reference_values`
     - Lazy registry and lookup helpers (``get``, ``all_names``,
       ``by_geometry``, ``by_groups``, ``by_method``).
   * - :mod:`~orpheus.derivations.common.withdrawal`
     - The lock on a withdrawn reference generator:
       :func:`~orpheus.derivations.common.withdrawal.withdrawn_generator`
       and its refusal
       :class:`~orpheus.derivations.common.withdrawal.GeneratorWithdrawn`,
       lifted by ``ORPHEUS_RUN_WITHDRAWN``, keyed on the value
       :class:`~orpheus.reference.withdrawal.Withdrawal`
       ``(reason, issue)``, which lives in the reference package
       (see :ref:`vv-withdrawn-generators`).
   * - :mod:`~orpheus.derivations.common.reference_body`
     - The one reading of a
       :class:`~orpheus.geometry.structured_geometry.StructuredGeometry`
       as the body a continuous reference generator (``Spectrum``,
       ``MomentSpace``, ``BasisSpace``, ``Billiard``) solves on.
       :func:`~orpheus.derivations.common.reference_body.reference_body`
       is total and knows no solver: after merging adjacent intervals
       of one material into runs, it returns exactly one of
       :class:`~orpheus.derivations.common.reference_body.HomogeneousBody`
       ``(coord, extent_cm, mat_id)``,
       :class:`~orpheus.derivations.common.reference_body.HollowBody`
       ``(coord, inner_radius_cm, outer_radius_cm, mat_id)``,
       :class:`~orpheus.derivations.common.reference_body.ReflectedSlab`
       ``(core_width_cm, reflector_width_cm, core_mat_id,
       reflector_mat_id)`` or
       :class:`~orpheus.derivations.common.reference_body.LayeredBody`
       ``(coord, breakpoints, mat_ids)``; ``ReferenceBody`` is their
       union. The boundary laws are read by
       :func:`~orpheus.derivations.common.reference_body.specular_albedo`
       ``(law, owner=)``, the one reader of a law as a specular albedo
       (vacuum 0, a mirror 1, a partial specular law its albedo; any
       other law refused),
       :func:`~orpheus.derivations.common.reference_body.specular_albedos`
       (every boundary point of a geometry, inner first) and
       :func:`~orpheus.derivations.common.reference_body.require_vacuum`.
       Each generator refuses a shape or a law it does not serve
       through
       :func:`~orpheus.derivations.common.reference_body.refuse_unserved`
       ``(what, owner=, missing=)``, which raises
       ``NotImplementedError`` naming the owner, the refused
       configuration (a body is worded by
       :func:`~orpheus.derivations.common.reference_body.describe`) and
       the solver it lacks (the table of which generator serves which
       shape under which laws: :ref:`structured-geometry-reference-body`).

The angular measure and its two lifts (Branch 1)
------------------------------------------------

The continuous measure on the direction sphere, its retraction, section
and pullback, with the mass :math:`4\pi` derived by integration and never
typed; the theory is :ref:`structured-geometry-two-lifts-branch-1`.

.. automodule:: orpheus.derivations.common.angular_measure
   :members:

Reference-value registry
------------------------

.. automodule:: orpheus.derivations.reference_values
   :members:

Verification case type
----------------------

.. automodule:: orpheus.derivations.common.verification_case
   :members:

Withdrawn reference generators
------------------------------

.. automodule:: orpheus.derivations.common.withdrawal
   :members:

The body a reference generator solves on
----------------------------------------

.. automodule:: orpheus.derivations.common.reference_body
   :members:

Homogeneous
-----------

.. automodule:: orpheus.derivations.continuous.analytical.homogeneous
   :members:

Slab Collision Probability
--------------------------

.. automodule:: orpheus.derivations.continuous.flat_source_cp.slab
   :members:
   :exclude-members: _XS_A, _XS_B

Cylindrical Collision Probability
---------------------------------

.. automodule:: orpheus.derivations.continuous.flat_source_cp.cylinder
   :members:
   :exclude-members: _XS_A, _XS_B

Diffusion
---------

.. automodule:: orpheus.derivations.continuous.cases.diffusion
   :members:

The reference kernel's dense pencil
-----------------------------------

.. automodule:: orpheus.derivations.common.dense_pencil
   :members:

The characteristic reference
----------------------------

The package, built rung by rung beside the trajectory-resolvent family;
the theory is :ref:`theory-characteristic-reference`. Its eight modules
are documented below; the package re-exports the public names of every
module but the gradings, which are imported by module.

.. automodule:: orpheus.derivations.continuous.characteristic

.. automodule:: orpheus.derivations.continuous.characteristic.walls
   :members:

.. automodule:: orpheus.derivations.continuous.characteristic.closure
   :members:

.. automodule:: orpheus.derivations.continuous.characteristic.basis
   :members:

.. automodule:: orpheus.derivations.continuous.characteristic.transport
   :members:

.. automodule:: orpheus.derivations.continuous.characteristic.assembly
   :members:

.. automodule:: orpheus.derivations.continuous.characteristic.cross_sections
   :members:

.. automodule:: orpheus.derivations.continuous.characteristic.system
   :members:

.. automodule:: orpheus.derivations.continuous.characteristic.grading
   :members:

Kernels
-------

.. automodule:: orpheus.derivations.common.kernels
   :members:

Quadrature contract
-------------------

.. automodule:: orpheus.derivations.common.quadrature
   :members:

Quadrature recipes (geometry-aware)
-----------------------------------

.. automodule:: orpheus.derivations.common.quadrature_recipes
   :members:

Cross-section library
---------------------

.. automodule:: orpheus.derivations.common.xs_library
   :members:
