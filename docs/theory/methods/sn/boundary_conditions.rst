.. _boundary-conditions:

Boundary Conditions
===================

Infrastructure
--------------

Boundary conditions are declared on the **geometry mesh** and resolved
by the **solver's augmented mesh** at construction time.  This two-stage
design separates physics intent (what condition to apply) from solver
mechanics (how to enforce it in the :term:`sweep`).

**Stage 1 --- Geometry declaration.**
:class:`~geometry.mesh.Mesh1D` carries ``bc_left: BC | None`` and
``bc_right: BC | None``; :class:`~geometry.mesh.Mesh2D` carries
``bc_xmin``/``bc_xmax``/``bc_ymin``/``bc_ymax: BC | None``.
:class:`~geometry.mesh.BC` is a frozen dataclass with two fields:

- ``kind: str`` --- an identifier such as ``"vacuum"``, ``"reflective"``,
  or ``"white"``.
- ``params: dict[str, float]`` --- optional numeric parameters
  (e.g. ``{"albedo": 0.7}``).

Convenience instances are available for the common cases:
:attr:`BC.vacuum <geometry.mesh.BC.vacuum>`,
:attr:`BC.reflective <geometry.mesh.BC.reflective>`, and
:attr:`BC.white <geometry.mesh.BC.white>`.
When a face is left as ``None``, the solver applies its own default
(reflective for the SN solver, matching the infinite-lattice /
eigenvalue convention).

**Stage 2 --- Solver resolution via the BC realizer.**
:class:`SNMesh` owns a class-level
:attr:`~SNMesh.BOUNDARY_OPERATOR_REGISTRY` mapping kind strings to
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` **subclasses**
(post Wave 8 of the trace-law refactor in
``.claude/plans/transient-giggling-cake.md``)::

    BOUNDARY_OPERATOR_REGISTRY = {
        "vacuum":     VacuumInflow,
        "reflective": ReflectiveBoundary,
    }

The registry values are the law classes themselves, not factory
functions. The pre-refactor ``_sn_vacuum_boundary_operator`` /
``_sn_reflective_boundary_operator`` factories were retired; their
job is now done by :meth:`SNMesh._resolve_one`, which dispatches
through :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer`
**uniformly** for every supported mesh (1-D Cartesian, 1-D
spherical, 1-D cylindrical, 2-D Cartesian) — see
:ref:`bc-sn-resolution-table` below. Issue #188 (curvilinear support
in the trace space — then named ``InflowTraceSpace``, now the unified
:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`) and Issue
#176 (drop 2-arg ``apply`` + simplify shim) collapsed the pre-cleanup
Cartesian-vs-curvilinear bypass into a single realizer-routed
path; details at :ref:`bc-curvilinear-realizer-unification`.

During ``SNMesh.__init__``, each face's :class:`~geometry.mesh.BC`
is looked up in the registry.  If the kind is not found, a
``ValueError`` lists the supported kinds.  For curvilinear
geometries (spherical, cylindrical), only ``"reflective"`` and
``"vacuum"`` are currently supported --- requesting any other kind
on a curvilinear mesh raises ``ValueError``.

The two-arg legacy interface ``bc.apply(angular_flux_outgoing,
quadrature)`` is **retired entirely** post Issue #186 / B3 + β2
(2026-05-11). Concrete BC :meth:`apply` methods no longer exist;
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` is a **pure
descriptor** with no callable interface. Rank-N composition uses
the descriptor-tree algebra
(:class:`~orpheus.geometry.boundary.LawSum` /
:class:`~orpheus.geometry.boundary.LawScaled`) with
:func:`~orpheus.geometry.boundary.realize_recursively` as the
sole descriptor→operator type transformer. See
:ref:`bc-trace-law-descriptor-model` for the design rationale.

The resolved BCs at ``sn_mesh.bc["xmin"]`` etc. expose the uniform
1-arg contract through the
:class:`~orpheus.geometry.boundary._bound_compat._BoundBoundaryOperator`
shim — internal to the package, not in
:attr:`orpheus.geometry.boundary.__all__`. Post Issue #186 the
shim is a **strict 1-arg passthrough** (extra args raise
:class:`TypeError`); it carries the originating ``BC.kind`` string
so the ``sn_mesh.bc["xmin"] == "vacuum"`` diagnostic
comparison continues to evaluate True iff the underlying law is
:class:`VacuumInflow`. (C4 / #220 re-keyed this surface from the
per-attribute ``sn_mesh.bc_left`` to the face-name-keyed
:attr:`SNMesh.bc` dict — see :ref:`bc-face-name-carve`.)
See :ref:`bc-tensor-decompositions` below
for the operator-algebra view and
:ref:`theory-boundary-conditions` for the full trace-law /
realizer architecture.

**Backward compatibility.**
:func:`solve_sn_fixed_source` still accepts a ``boundary_condition: str``
parameter (default ``"vacuum"``).  Internally it calls
``_apply_default_bcs(mesh, boundary_condition)``, which applies the
string to **all faces** that lack explicit :class:`~geometry.mesh.BC`
declarations.  When the mesh already carries explicit BCs, the parameter
is silently ignored --- mesh-level declarations always take precedence.
:func:`solve_sn` (the eigenvalue entry point) does not expose a
``boundary_condition`` parameter; eigenvalue problems use whatever the
mesh declares (defaulting to reflective on all faces).

.. note::

   Before this infrastructure existed, the SN solver hardcoded
   reflective BCs on all faces and the then-production ``transport_sweep``
   entry accepted a ``boundary_condition: str`` parameter.  That parameter
   has been removed --- BCs now flow exclusively through the mesh → SNMesh
   resolution path described above.

Supported Types
---------------

**Reflective** (specular reflection).
At the outer boundary :math:`r = R` (or :math:`x = L`), the incoming
flux for :term:`ordinate` :math:`n` is set to the outgoing flux of its reflected
partner:

.. math::
   :label: reflective-bc

   \psi_n^{\rm in} = \psi_{n'}^{\rm out}

where :math:`n'` is the reflected partner ordinate (negating the
appropriate direction cosine).  Reflective partner indices are precomputed
by each :term:`quadrature`'s :meth:`reflection_index` method.  This is the
default for eigenvalue problems (infinite lattice / infinite medium).
The CP solver uses white (isotropic) BCs instead; see
:ref:`white-bc-quality` for a comparison showing the ~1% gap between
the two approaches.

**Vacuum** (zero incoming flux).
All incoming angular fluxes at the face are set to zero:

.. math::
   :label: vacuum-bc

   \psi_n^{\rm in} = 0

.. (vv-status rationale) definition: Definitional / notation introduction. Operational rule ψ_n^in = 0 for vacuum boundary; semantics exercised by every vacuum-BC test (test_boundary_conditions, MMS suite); no isolated identity to verify.
.. vv-status: vacuum-bc documented


In the 1-D cumprod path, this means the recurrence starts from zero
instead of the reflected outgoing flux.  In the 2-D wavefront sweep,
the reflective-partner copy is skipped, leaving incoming-face angular
fluxes at their zero initialisation.  Vacuum BCs are the natural
choice for fixed-source MMS verification on finite slabs (see
:ref:`sn-mms-verification`).

.. _bc-tensor-decompositions:

Boundary conditions as tensor decompositions
---------------------------------------------

The boundary conditions used by :class:`SNMesh` are concrete instances
of a more general tensor-decomposed framing, defined in
:mod:`orpheus.geometry.boundary`. A boundary condition is a linear
operator :math:`R` mapping the outgoing angular flux at a face to the
incoming angular flux:

.. math::
   :label: bc-tensor-decomposition

   \psi_{\rm in}(\Omega)
   = (R\,\psi_{\rm out})(\Omega)
   = \sum_\alpha \bigl(G_\alpha\,\psi_{\rm out}\bigr)(\Omega) \cdot A_\alpha,

.. (vv-status rationale) Definitional/literature-transcribed framing (a BC
   as R = sum_alpha G_alpha A_alpha, Lewis & Miller 1984 §3.4). The concrete
   rank-N primitives (vacuum/reflective/white) carry their own verification.
.. vv-status: bc-tensor-decomposition documented

where :math:`G_\alpha` is a **geometric operator** (permutation,
pushforward, angular average, spatial wrap) and :math:`A_\alpha` is a
**scalar amplitude** (typically an :term:`albedo` :math:`\in [0, 1]`).

This is the same algebra Lewis & Miller (1984) §3.4 use to introduce
boundary conditions in transport: every BC of practical interest is
either rank-1 (one :math:`G \otimes A` term) or a finite linear
combination of rank-1 primitives (rank-N). The implemented primitives
are:

.. list-table:: Implemented :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` primitives (Wave 7 vocabulary)
   :widths: 20 30 25 25
   :header-rows: 1

   * - Class
     - :math:`G_\alpha`
     - :math:`A_\alpha`
     - Rank / wired into ``solve_sn``
   * - :class:`~orpheus.geometry.boundary.VacuumInflow`
     - :math:`0` (no operator)
     - 0
     - 0 / yes
   * - :class:`~orpheus.geometry.boundary.ReflectiveBoundary`
     - permutation under reflection axis
     - albedo (1 = perfect)
     - 1 / yes
   * - :class:`~orpheus.geometry.boundary.WhiteBoundary`
     - cosine-weighted hemispheric average
     - albedo
     - 1 / no (Wave C)
   * - :class:`~orpheus.geometry.boundary.PeriodicBoundary`
     - spatial pushforward (caller-supplied)
     - 1
     - 1 / no (Wave C/D)
   * - :class:`~orpheus.geometry.boundary.AlbedoBoundary`
     - identity in angle
     - albedo
     - 1 / no (building block)
   * - :class:`~orpheus.geometry.boundary.PrescribedInflow`
     - 0
     - 0
     - 0 with :math:`q \neq 0` / no (rank-0 source-only affine BC)

The pre-Wave-7 names ``VacuumBoundaryOperator`` /
``SpecularBoundaryOperator`` / ``WhiteBoundaryOperator`` /
``PeriodicBoundaryOperator`` / ``AlbedoBoundaryOperator`` were
**retired in Wave O step O.4a.1**; their canonical successors
(:class:`~orpheus.geometry.boundary.VacuumInflow` /
:class:`~orpheus.geometry.boundary.ReflectiveBoundary` /
:class:`~orpheus.geometry.boundary.WhiteBoundary` /
:class:`~orpheus.geometry.boundary.PeriodicBoundary` /
:class:`~orpheus.geometry.boundary.AlbedoBoundary`) are the sole
live names in :mod:`orpheus.geometry.boundary`.
``MixedBoundaryOperator`` was **retired in Wave 11**; rank-N
(Marshak, partial-current) boundaries are now expressed via the
**descriptor-tree algebra**
(:class:`~orpheus.geometry.boundary.LawSum` /
:class:`~orpheus.geometry.boundary.LawScaled`) on the unrealised
laws:

.. code-block:: python

   tree = 0.3 * spec + 0.7 * white            # LawSum of LawScaled
   op = realize_recursively(tree, method_space)  # OperatorSum of ScaledOperator
   psi_in = op.apply(psi_out)

See :ref:`bc-rank-n-algebra` for the closed algebra and the
:ref:`bc-realize-recursively` walker that lowers a descriptor
tree to a Wave-0 operator tree.

The abstract base :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`
is a **pure descriptor** post Issue #186 / B3 + β2 (2026-05-11)
— it has **no** :meth:`apply` method. The :class:`LinearOperator`
inheritance that historically supplied ``apply`` was removed; the
concrete laws likewise carry no ``apply`` / ``apply_transpose``
methods. The §16A.3 three-layer architecture (descriptor /
realizer / operator) is now enforced by the type system: a static
type checker rejects ``law.apply(...)`` at the linter level
without running the program. The full retrospective on the
predecessor Option A and β1 forms (and why each was rejected) is
at :ref:`bc-trace-law-descriptor-model`.

The tensor framing pays off architecturally because partial-current
Marshak boundaries (:math:`R = c_1 \, G_{\rm refl} + c_2 \, G_{\rm
diff}`, Bell & Glasstone 1970 §1.5) and multi-region interface
couplings are all instances of the same algebra: pick the geometric
operators, pick the amplitudes, sum. New BCs are one
:class:`BoundaryTraceLaw` subclass + one
``BOUNDARY_OPERATOR_REGISTRY`` entry away — no sweep edits per BC.

.. _bc-sn-resolution-table:

SN BC resolution table
----------------------

The :meth:`SNMesh._resolve_one` dispatch is summarized below.
Each row maps the user-facing :class:`~orpheus.geometry.mesh.BC`
kind string to (a) the resolved :class:`BoundaryTraceLaw`
subclass, (b) the :class:`SNBoundaryRealizer.realize` output
operator, and (c) the ``creates_sweep_cycle`` flag used by §15A.2
sweep-cycle detection (see :ref:`bc-sweep-cycle`). The realizer
dispatch is **uniform** across every supported mesh — 1-D
Cartesian / spherical / cylindrical and 2-D Cartesian — post
Issue #188 + #176.

.. list-table:: BC.kind → law class → realized SN operator
   :header-rows: 1
   :widths: 16 22 36 12 14

   * - ``BC.kind``
     - Law class
     - Realized SN operator
     - α
     - ``creates_sweep_cycle``
   * - ``"vacuum"``
     - :class:`~orpheus.geometry.boundary.VacuumInflow`
     - :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`
       (per-face inflow indices)
     - —
     - ``False``
   * - ``"reflective"``
     - :class:`~orpheus.geometry.boundary.ReflectiveBoundary`
     - :class:`~orpheus.numerics.operator.PermutationOperator`
       (``quadrature.reflection_index(axis)``)
     - 1 (fast path)
     - **``True``**
   * - ``"reflective"``
     -
     - ``α * PermutationOperator``
       (:class:`~orpheus.numerics.operator.ScaledOperator`)
     - α ≠ 1
     - **``True``**
   * - ``"white"``
     - :class:`~orpheus.geometry.boundary.WhiteBoundary`
     - :class:`~orpheus.sn.boundary.angular.AngularAverageOperator`
     - 1 (fast path)
     - ``False``
   * - ``"white"``
     -
     - ``α * AngularAverageOperator``
     - α ≠ 1
     - ``False``
   * - ``"periodic"``
     - :class:`~orpheus.geometry.boundary.PeriodicBoundary`
     - :class:`~orpheus.numerics.operator.PeriodicWrapOperator`
     - 1
     - **``True``**
   * - ``"albedo"``
     - :class:`~orpheus.geometry.boundary.AlbedoBoundary`
     - :class:`~orpheus.numerics.operator.ZeroOperator`
     - 0
     - ``False``
   * - ``"albedo"``
     -
     - :class:`~orpheus.numerics.operator.IdentityOperator`
     - 1
     - ``False``
   * - ``"albedo"``
     -
     - ``α * IdentityOperator``
     - α ∉ {0, 1}
     - ``False``
   * - ``"prescribed_inflow"``
     - :class:`~orpheus.geometry.boundary.PrescribedInflow`
     - :class:`~orpheus.sn.boundary.angular.IncomingSourceOperator`
       (source.evaluate; ignores outgoing flux)
     - —
     - ``False``

The :meth:`SNMesh._resolve_one` dispatch constructs the resolved
operator via :meth:`SNBoundaryRealizer.realize(law, method_space)`
where the ``method_space`` is built by
:meth:`SNMethodSpace.for_face` carrying the precomputed unified
:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` (built
once at :class:`SNMesh` construction for every supported mesh).
The reflective branch derives its reflection axis from the face's
own :class:`~orpheus.transport.mesh.axis.FaceLabel` —
``AXIS_NAMES[label.axis_index]`` — so the partner is correct at any
dimension by construction (C4 / #220; see
:ref:`bc-face-name-latent-d3-bug`). The
:class:`~orpheus.geometry.boundary._bound_compat._BoundBoundaryOperator`
shim wraps the result with a ``kind`` tag for the
``sn_mesh.bc["xmin"] == "vacuum"`` string-equality surface.

.. note::

   **The 1-D y-face placeholders were retired in C4 / #220.** Pre-C4,
   a slab :class:`SNMesh` carried a pair of realized no-op
   ``ReflectiveBoundary(axis="y")`` operators at ``bc_ymin`` /
   ``bc_ymax`` so cross-dimensional code could read them without
   coord-system gating — but **no production code ever read them**
   (a 1-D mesh's ``trace.layout.faces`` is ``("xmin", "xmax")``).
   C4 makes them unrepresentable: a slab has no y-axis in its
   :attr:`~SNMesh.axes` tuple, so
   :func:`~orpheus.transport.mesh.axis.face_labels` emits no y-label and
   :attr:`SNMesh.bc` has no y-entry — ``slab.bc["ymin"]`` is a
   :class:`KeyError`, not a no-op. See
   :ref:`bc-face-name-carve-what-retired` for the full retirement
   record (the pre-C4 "why the placeholders were once safe"
   rationale is preserved there).

**Pre-cleanup history.** Before Issue #188 + #176 (closed
2026-05-11), curvilinear ``Mesh1D`` bypassed the realizer because the
trace factory (then ``InflowTraceSpace.from_mesh_and_quadrature``, since
C5.3 the geometry-blind
:meth:`AngularTraceSpace.from_quadrature_and_layout
<orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.from_quadrature_and_layout>`)
raised :class:`NotImplementedError` on those coord systems; the
``_BoundBoundaryOperator`` shim carried a dual mode where the
``quadrature=`` kwarg, when non-``None``, bound an
``AngularQuadrature`` and forwarded ``inner.apply(psi,
bound_quad)`` to the legacy 2-arg :class:`BoundaryTraceLaw` body.
The bypass and dual-mode are gone; details and the algebraic
sequence ("Issue #188 unblocks Issue #176") at
:ref:`bc-curvilinear-realizer-unification`.

Inner Boundary (Curvilinear)
-----------------------------

At :math:`r = 0`:

- The face area :math:`A(0) = 0`, so **no spatial flux crosses the
  origin**.  The spatial incoming flux for outward-sweeping ordinates is
  zero.
- The **angular redistribution provides the inward-to-outward
  transition**: flux entering as an inward-directed ordinate
  (:math:`\mu < 0` or :math:`\eta < 0`) is redistributed to outward
  ordinates (:math:`\mu > 0` or :math:`\eta > 0`) through the
  :math:`\alpha` coupling.

This means the curvilinear sweep does not need an explicit boundary
condition at :math:`r = 0` --- the geometry handles it naturally.
Curvilinear sweeps currently only support reflective BCs on the outer
face; this is enforced by the validation in
:meth:`SNMesh._resolve_one`.


