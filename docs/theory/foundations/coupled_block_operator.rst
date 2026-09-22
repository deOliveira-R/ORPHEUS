.. _coupled-block:

==========================
The Coupled Block Operator
==========================

.. contents:: Contents
   :local:
   :depth: 2


.. Machine header — the ``nexus-meta`` schema for this page (PROVISIONAL).
.. This page was extracted verbatim from ``operator_algebra.rst`` (#231
.. Phase 3); content is unchanged. The schema is provisional pending a
.. full re-audit of the split corpus.

.. dropdown:: Machine header — ``nexus-meta`` schema (PROVISIONAL)
   :color: muted

   .. code-block:: yaml

      module: transport
      concept: coupled_block_operator
      role: "the 2x2 coupled block operator for the curvilinear within-group problem — the psi-half ray as a first-class System B; the four named blocks (A_AA the loss composite, seed / emission / march couplings), the N-general block machinery, and the structure-keyed block solve"
      depends_on: [operator_algebra]
      related: [operator_tensor_network, boundary_conditions]
      status: "extracted from operator_algebra.rst; content verbatim, provisional header"


This page develops the **coupled block operator** — the 2×2 block system
that gives the curvilinear-S\ :sub:`N` within-group problem its capstone
form, promoting the :term:`starting-direction <starting direction>` flux :math:`\psi_{1/2}` (the ψ½
ray) to a first-class **System B** alongside the angular-flux **System
A**. It is the block-level generalisation of the operator algebra
developed in :doc:`/theory/foundations/operator_algebra`: where that
page types the operators of a single within-group system (the loss
composite :math:`A = L + C - S - N_{2n} - B`), this page types the four blocks
:math:`A_{AA}`, :math:`A_{AB}`, :math:`A_{BA}`, :math:`A_{BB}` of the
coupled 2×2 operator — whose diagonal :math:`A_{AA}` block IS that loss
composite, and whose off-diagonal blocks are the named seed and emission
couplings, so that a wrong pairing is **unconstructable by type**.

.. _coupled-block-operator:

The coupled block operator — the ψ½ ray as System B
===================================================

The curvilinear-S\ :sub:`N` within-group problem is the operator-surface
taxonomy's capstone: a **2×2 coupled block operator** over two systems.
It is the block-level generalisation of the three-layer surface
(:ref:`capability-set-semantics`) — ``apply`` becomes the block
matvec, ``assemble`` scatters each block at its DOF offset, and ``solve``
runs structure-keyed block substitution — and it is where the
curvilinear starting-direction flux :math:`\psi_{1/2}` finally earns a
first-class home. The campaign that built it (``coupled_block_operator``,
GH #280 the walk unification / #282 the direct ψ½ seed) replaced a
*fused* implementation — the ψ½ seed hand-rolled inside the :term:`sweep` and
inside the model-generic scattering gain — with an explicit block system
in which every coupling is a named operator and a wrong pairing is
**unconstructable**. This section documents the posing, the four blocks,
the N-general machinery, the one production spelling, the block solve,
the convergence certificate, and the swap law. The starting-direction
physics (the pole as a straight characteristic, the M1/M2/M3 metric
distinction, R12a presence) lives in
:ref:`sn-direct-seed-solve` in
:doc:`/theory/methods/sn/curvilinear_one_group`;
here we document the **algebra** the physics is posed in.

Why ψ½ is a system, not a kwarg — the two-point radial BVP
----------------------------------------------------------

At the two closed :math:`\mu = \pm 1` rays of a curvilinear μ-level the
angular-redistribution coefficient :math:`1-\mu^2` vanishes
(:math:`\alpha_{1/2}=0`, Hébert 3.423), so the streaming–collision
balance for :math:`\psi_{1/2}` **decouples** from the α-cascade and
reduces to a plain :term:`diamond-difference <diamond difference>` recurrence in radius —
:math:`\mu\,\partial_r + \sigma_t`, a *straight characteristic* (Hébert
§3.9.4). This ODE is not a scalar boundary datum the bulk sweep can
consume as a kwarg: it is a **two-point boundary-value problem** in its
own right, carrying **two** boundary conditions —

* **r = R Dirichlet** — the outer-face inflow corner (:math:`\mu = -1`)
  is *given* data: 0 for vacuum, the reflected outflow corner for
  reflective;
* **r = 0 pole continuation** — :math:`\psi_{1/2}^{+}(0) =
  \psi_{1/2}^{-}(0)`, the inward leg's exit face is the outward leg's
  entry face.

A quantity closed by its own two-point BVP is a **system**, not a
parameter. Posing it as one lets the algebra answer "what invariant
tests it?" (the march is the exact inverse of its own recurrence),
"how is it coupled?" (two named off-diagonal operators), and "how is it
solved?" (block substitution) — the four operator-algebra questions the
fused kwarg spelling could not even ask.

The within-group augmented system is therefore

.. math::
   :label: coupled-block-2x2

   \begin{bmatrix} A_{AA} & A_{AB} \\ A_{BA} & A_{BB} \end{bmatrix}
   \begin{bmatrix} \psi_A \\ \psi_B \end{bmatrix}
   \;=\;
   \begin{bmatrix} q_A \\ q_B \end{bmatrix},
   \qquad A_{ij} : V_j \to W_i ,

.. vv-status: coupled-block-2x2 documented

with **System A** the transport bulk ⊕ trace (the angular-flux
:class:`~orpheus.transport.full_field.FullField`, ``BulkField`` bulk
angular flux ⊕ its spatial ``BoundaryField`` trace) and **System B** the
ψ½ radial-characteristic ray (the ``RadialCharacteristicField`` cells at
each radial cell, carrying its two BC data). This is not a new tensor
type — it is the composite biproduct re-partitioned,
:math:`\mathrm{Mat}_2(\mathrm{Mat}_2(\mathcal C)) \cong
\mathrm{Mat}_4(\mathcal C)`. The two-system membership is carried by the
:class:`~orpheus.numerics.operator.SystemRole` axis (below), orthogonal to
the bulk↔boundary :class:`~orpheus.numerics.operator.BlockRole` that
refines System A.

The four blocks
---------------

Each block is a named production operator with its own code home,
invariants, and adjoint. The diagonal blocks are the systems'
self-operators; the off-diagonals are the couplings. All four live in
:mod:`orpheus.sn.operators.radial_characteristic` (System A's self-block
alone has no explicit object — it is the driver-level composite).

.. list-table:: The 2×2 blocks and their realisations
   :header-rows: 1
   :widths: 8 30 34 28

   * - Block
     - Object
     - What it is
     - Invariants
   * - :math:`A_{AA}`
     - the driver composite :math:`L + C - S - N_{2n} - B_a`
     - System A's self-block: streaming + collision − bulk scattering
       gain − bulk :math:`(n,2n)` emission gain − trace boundary gain
     - carries :class:`BlockRole` ``FULL`` by the join lattice;
       ``solve`` = the WDD sweep
   * - :math:`A_{AB}`
     - :class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicSeeding`
     - ray → bulk seed injection (the Morel–Montry angular seed) —
       **cell-local angular**, no spatial cell-cell coupling
     - :math:`\sigma`-**independent** (a type fact); ``is_adjointable``,
       ``is_invertible`` ``False`` (rectangular)
   * - :math:`A_{BA}`
     - :class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicEmission`
     - bulk → ray coupling, :math:`\mathrm{Fold}\circ K\circ\int d\mu`
       (the emission the within-group gain lags)
     - ``is_adjointable`` (its transpose IS the S-adjoint pullback);
       ``is_invertible`` ``False``
   * - :math:`A_{BB}`
     - :class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator`
     - the banded radial DD march :math:`\mu\,\partial_r + \sigma_t`
     - ``is_invertible`` = ``is_adjointable`` = ``True``; the direct
       Carlson march IS :math:`A_{BB}^{-1}`

**A naming gotcha to spell out.** ``A_BB`` (the class
:class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator`)
is the **bare** radial march :math:`\mu\,\partial_r + \sigma_t`. System
B's self-block in the *loss* grid is :math:`A_{BB} - B_b`, the march minus
the ray-corner boundary gain :math:`B_b`
(:class:`~orpheus.sn.operators.boundary.RadialCharacteristicBoundaryOperator`)
— **exactly parallel to System A**, whose self-block is
:math:`A_{AA} = L + C - S - N_{2n} - B_a`. Both systems factor as
``(transport operator) − (bulk/reaction gains) − (boundary gain)``, and
both boundary gains :math:`B_a`, :math:`B_b` are lagged (they live in
:math:`N`, below), each reflecting through its **own** trace so the
:math:`-B` term arrives as given data (:ref:`bc-extraction`).

The self-block :math:`A_{BB}` — the march is its own inverse
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because the recurrence is triangular in radius (:math:`\rho = 0` — the
#284 forward-substitution certificate, measured :math:`A_{ss} = 5.0`
self-coupling, :math:`A_{sb} = 0` exact), the two-leg Carlson march
(:func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`,
Hébert 3.434–3.435) is the **exact** direct inverse :math:`A_{BB}^{-1}` —
no iteration, no previous-iterate seed. :meth:`RadialCharacteristicOperator.solve
<orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator.solve>`
marches the inward :math:`\mu=-1` leg from the r=R inflow corner to the
pole, then rides the *same* engine on reversed cell data (orientation
carried by the data, never a flag) from the pole-continued face out to
the r=R outflow corner. All four action surfaces are realised —
:meth:`apply` / :meth:`apply_transpose` (the forward
:math:`A_{BB}^{\mathsf T}`), :meth:`solve` / :meth:`solve_transpose`
(:math:`(A_{BB}^{-1})^{\mathsf T}`) — so ``is_invertible`` and
``is_adjointable`` both read ``True`` and the involution web
``inverse().solve == apply`` closes. The ``apply ∘ solve`` outflow-corner
defect closes to ``0.0`` bit-exactly; the cell round-trip is
principled-equivalent at ~FP ULP (the forward's :math:`2/\Delta r` and
the march's :math:`\Delta r\,\sigma + 2` reassociate). Since campaign
step 4e the production ``(L+C)`` walk routes its ray solve **through**
this operator — the former in-sweep inline march is retired, so the
two-leg orchestration lives in exactly one place (Cardinal Rule 2).

The seed block :math:`A_{AB}` — a σ-independent cell-local angular coupling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At each radial cell :math:`i`, the ray value :math:`\psi_{1/2}(i)` seeds
the Morel–Montry **weighted** angular recurrence
(:cite:`MorelMontry1984`; implemented form :cite:`BaileyMorelChang2010`
Eqs. (42)/(43) — Hébert defines no :math:`\tau`, see
:ref:`sn-tau-source-of-record`)

.. math::
   :label: coupled-ab-seed

   \psi_{m+1/2,\,i}
     = \frac{\psi_{m,\,i} - (1-\tau_m)\,\psi_{m-1/2,\,i}}{\tau_m},
   \qquad \psi_{-1/2,\,i} \equiv \psi_{1/2}(i),

run over :term:`ordinates <ordinate>` :math:`m` at a **fixed** cell :math:`i`. The upstream
half-flux enters that cell's balance as the angular numerator
:math:`(\Delta A/w)\,c_{\rm in}\,\psi_{m-1/2,\,i}`. So the seed at cell
:math:`i` couples ONLY to the bulk ordinates at the *same* cell — there
is **no spatial cell-cell coupling** (that is :math:`A_{BB}`'s job). This
is what separates :math:`A_{AB}` from :math:`A_{BB}`: :math:`A_{AB}` is
cell-local angular and realises both directions HERE as thin wraps of the
single-sourced Morel–Montry closure (:meth:`apply` zeroes the bulk to
isolate its block by linearity; :meth:`apply_transpose` runs the local
gather :math:`\bar n = -\bar o/V` then reverses the recurrence to the seed
cotangent ``seed_cells_bar`` — the Euclidean transpose
:math:`A_{AB}^{\mathsf T}`). It needs **no** :math:`\sigma_t`: with the
bulk zeroed the collision/streaming terms drop out and only the
σ-independent angular numerator survives, so the coupling is a pure
function of the mesh geometry and :term:`quadrature`. That σ-independence is a
**type fact** — the constructor takes only the mesh.

.. vv-status: coupled-ab-seed documented

The emission block :math:`A_{BA}` — the Schur fold ``Fold ∘ K ∘ integrate``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The bulk-within-group flux emits an isotropic source that seeds the ray.
:math:`A_{BA}`
(:class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicEmission`)
composes three factors,

.. math::
   :label: coupled-ba-emission

   A_{BA} \;=\; \underbrace{\mathrm{Fold}}_{\text{Reconstruction}}
            \;\circ\; \underbrace{K}_{\text{emission kernel}}
            \;\circ\; \underbrace{\textstyle\int\! d\mu}_{\text{integrate}} ,

the angular integral of the bulk flux to :math:`\phi_0`, the operator's
isotropic emission kernel :math:`K`, and the reconstruction of that
:math:`\ell = 0` moment at the closed rays. The kernel is a **dependency
injection**: :math:`A_{BA}` is generic over any ``ndarray → ndarray``
emitter carrying ``apply``/``apply_transpose``. The production
instantiation passes the solver-composed :math:`K_{\rm iso} =
\Sigma_{s0} + \nu_{2n}\Sigma_{2n}` — assembled at the ONE within-group
construction site as ``S.isotropic_energy + N2N.isotropic_energy``
(:attr:`TransferOperator.isotropic_energy
<orpheus.transport.operators.transfer.TransferOperator.isotropic_energy>`
on each — CS4c §14.1: the two cached energy bindings the bulk gains
also consume), so
the emission is single-sourced per LEAF (one shared kernel object per
channel, never a twin re-implementation of
:math:`\Sigma_{s0}^{\mathsf T}\phi`).

.. note::

   **The composition site moved, and the reason is the sharper record
   than the move.**  Until 2026-08-30 this read *"passes the scattering
   operator's* ``isotropic_kernel`` *— the same shared object the bulk
   scattering gain uses"*: :math:`S` owned an accessor that summed its
   own :math:`P_0` energy binding with the :math:`(n,2n)` one and
   handed the sum out.  That accessor was the **operator-level
   bundling** the CS4c §14.1 ruling forbids — a scattering operator
   deciding, on every consumer's behalf, that :math:`(n,2n)` groups
   with scattering.  It retired with the extraction, and the sum is now
   written where the grouping is a legitimate local choice: the
   builder, at the one construction site, composes
   ``S.isotropic_energy + N2N.isotropic_energy`` and hands the
   resulting :class:`~orpheus.numerics.operator.OperatorSum` to the
   emission.  Single-sourcing is preserved per LEAF — the two energy
   bindings are the same objects the bulk gains hold — and what is
   gained is that a consumer wanting :math:`\Sigma_{s0}` *without* the
   multiplicity can now spell it.

   ⚠ The :math:`(n,2n)` half was spelled ``N2N.energy`` until #426 step
   2 (2026-09-04); it is ``N2N.isotropic_energy`` now, because both
   gains are roles of one
   :class:`~orpheus.transport.operators.transfer.TransferOperator` and
   that accessor is the core's.

   ⭐ **And the same day (CS4c step 5) the SIBLING slot came up to this
   level.**  This within-group site had always passed the *operators*;
   the daggered eigen-posing's own :math:`(B,A)` fission fold passed
   ``F.kernel`` — the rank-1 tensor product *inside* the energy binding —
   so one dependency-injection slot was fed two different levels of one
   abstraction.  It now takes ``F.isotropic_energy``, the operator, and
   :meth:`IsotropicFission.apply_transpose
   <orpheus.transport.operators.isotropic_transfer.IsotropicFission.apply_transpose>`
   is reached instead of stepped over on the adjoint path
   (:ref:`sn-adjoint-coupled-posing` in
   :doc:`/theory/methods/sn/adjoint`).  Both spellings computed identical
   numbers, which is why no value record could have flagged it — a claim
   about the ROUTE needs an instrument on the route
   (``vv-principles`` #26).  ⭐ **And this** :math:`K_{\rm iso}`
   **is** :math:`\ell = 0` **by PHYSICS, so step 2 did not touch it.**
   The ray seed is driven by the scalar flux; what the emission needs
   is each gain's :math:`P_0` energy binding, and it would need exactly
   that even if the solve ran at :math:`L = 6`.  The step-2 change —
   the :math:`(n,2n)` gain reading its Legendre stack at the solve's
   order — lives on the BULK arm, where the angular shape of the source
   is a real degree of freedom.  A reader who reads
   :ref:`the (n,2n) P0-truncation record <sn-n2n-p0-truncation>` and
   comes here looking for a second truncation to fix will not find
   one.

.. vv-status: coupled-ba-emission documented

The **fold factor** :math:`\mathrm{Fold}`
(:class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicReconstruction`)
is the 1-D angular reconstruction (Hébert Eq. 3.432) *sampled* at the
closed rays :math:`\mu = \pm 1`:

.. math::
   :label: coupled-ba-fold

   \bar q_{1/2}(\mu = \pm 1)
     \;=\; \sum_\ell \frac{2\ell + 1}{2}\,q_\ell\,(\pm 1)^\ell ,

the same :math:`\tfrac{2\ell+1}{2}\,(\pm 1)^\ell` weights the angular
frame reconstructs with (:math:`P_\ell(\pm 1) = (\pm 1)^\ell`), evaluated
at the rays rather than the quadrature nodes. Both the forward fold and
its Euclidean transpose route through the ONE math kernel
:func:`~orpheus.numerics.spaces.radial_characteristic_space.fold_moments_to_radial_characteristic`
(Cardinal Rule 2 — the :math:`P_1(-1) = -1` sign is spelled once). The
transpose is the exact **broadcast/sum adjoint pair**: the forward folds
one moment source onto every carried level and both signs (a broadcast);
the transpose *sums* the per-level, per-sign ray cotangents back into
moment space. Corners stay zero — the fold is volumetric; the
inflow-corner datum is :math:`B_b`'s job.

.. vv-status: coupled-ba-fold documented

**The LIFT — S/F pure-bulk, the driver composes the emission.** Before
campaign step 4c the model-generic scattering/fission composites
hand-rolled the ψ½ seed inside their own ``apply`` — a curvilinear-SN
augmentation welded into a model-generic gain. The **lift** (the Wave-O
#208 pattern that separated :math:`B` from :math:`S`) reversed that: it
made :math:`S`/:math:`F` **pure bulk** and posed the coupling as a
first-class operator the driver lags as a separate gain. The
within-group gain rides :math:`(S,\ N_{2n},\ A_{BA},\ B_a)`, and
:math:`A_{BA}`'s transpose is exactly the ``w · K_isoᵀ(Reconstructionᵀ
χ_seed)`` bulk pullback the S-adjoint used to carry inline — moved HERE,
single-sourced, so ``S.apply_transpose`` is pure bulk and
:math:`(L+C-S-N_{2n}-A_{BA}-B).H` reconstructs the monolithic adjoint.
Because
the forward keeps :math:`K_{\rm iso}` while the adjoint rides the fold
transpose, the reciprocity gate is a genuine cross-check of two
structurally-different representations of one operator, not a tautology.

.. note::

   **Fission does NOT flow through** :math:`A_{BA}` **(HAZARD 5).** The
   emission is kernel-generic and accepts the fission dyad
   :math:`\chi \otimes \nu\Sigma_f` as a smoke-verified second kernel, but
   fission's *production* ray seed rides the outer :math:`q_{\rm ext}`
   seam as a **direct** :math:`\mathrm{Fold}`: fission is the eigenvalue
   outer source, so its :math:`K\circ\int d\mu` is already pre-computed as
   the fission source, and routing it through the full
   :math:`\mathrm{Fold}\circ K\circ\int d\mu` would **double-apply**
   :math:`K\circ\int`. Within-group fission is zero (it enters as
   :math:`q_{\rm ext}` per the eigenvalue-outer / within-group split), so
   :math:`A_{BA}`'s production shape is the *scatter* coupling; the
   genericity keeps the emission a clean dependency injection, not a claim
   that fission wires through it.

   ⚠ **Two grids, one builder** (since 2026-09-13).
   :func:`~orpheus.sn.coupled_system.build_within_group_system` now poses
   the **production** :math:`F` on the same carrier as the loss and
   stores it as
   :attr:`WithinGroupSystem.production
   <orpheus.sn.coupled_system.WithinGroupSystem.production>`, and *that*
   grid's :math:`(B,A)` row **is** the fission ray fold — the
   kernel-generic emission carrying ``F.isotropic_energy``.  HAZARD 5 is
   unchanged and is about which grid the fold belongs in: out of the
   within-group **gain** :math:`N` (where it would double-apply
   :math:`K\circ\int`), into the eigen-:math:`M` **posing** (where it
   belongs).  The forward path is likewise unchanged — its ray seed is
   still the direct :math:`\mathrm{Fold}` of the already-assembled
   fission source (``_radial_characteristic_fission_seed``).  What moved
   is only that the posing is built at this one site instead of inside
   the adjoint entry (:ref:`sn-one-fission-per-problem`).

.. _coupled-block-n-general-machinery:

The N-general block machinery
-----------------------------

The 2×2 ψ½ system is **instance #1** of a semantics-agnostic N-system
machinery in :mod:`orpheus.numerics.coupled_system` — nothing there knows
transport, rays, or meshes. A coupled system is N sub-systems solved
together: the state is the direct sum of the systems' fields
(:class:`~orpheus.numerics.coupled_system.CoupledField`), the space is
their product (:class:`~orpheus.numerics.coupled_system.CoupledSpace`),
and the operator is the N×N grid of blocks
(:class:`~orpheus.numerics.coupled_system.CoupledOperator`) with the block
matvec

.. math::
   :label: coupled-block-matvec

   y_i \;=\; \sum_j A_{ij}\, x_j

as the *only* spelling. A missing block (``None``) **is** the zero map,
structurally — block existence doubles as coupling sparsity, and no
zero-padding arithmetic ever runs.

.. vv-status: coupled-block-matvec documented

**Why a typed grid over a padded** ``OperatorSum`` **(the rejected
route).** The alternative — a flat sum of same-space operators, each
padding the blocks it does not touch with present-zeros — was **rejected**
(campaign re-scope, 2026-07-10): padding keeps *wrong multiplications
representable*. A padded block accepts any composite; nothing at the type
level stops a ray operator from receiving a bulk field or an emission
from landing in the trace slot, and ``system_role`` tags are runtime
metadata, not type constraints. The honest object is a typed block
operator over a typed block vector, where every term of
:eq:`coupled-block-matvec` type-checks per block and a wrong pairing is
**unconstructable** (``coding-elegance`` Pattern 1 ∘ Pattern 4 — the
algebra is the syntax).

**The three block-level modes** mirror the three-layer operator surface,
one level up: ``apply`` (the block matvec, Krylov action), ``assemble``
(each block's sparse emission scattered at its ``(row_i, col_j)`` DOF
offset into one flat matrix via :func:`scipy.sparse.block_array`, closing
over :class:`~orpheus.numerics.assembled_operator.SparseAssembledOperator`
— the same axis as :ref:`operator-algebra-assembly-axis`), and ``solve``
(the structure-keyed direct solve, below). The offsets ARE the block
structure — :attr:`CoupledSpace.system_slices
<orpheus.numerics.coupled_system.CoupledSpace.system_slices>` is the
scoped local→global DOF map ``assemble`` needs, the same layout
:meth:`CoupledField.to_flat` packs. Because :class:`CoupledField` itself
satisfies the ravellable ``to_flat``/``from_flat`` protocol, every
``restart = n_dof = template.to_flat().size`` sizing site tracks the
coupled dimension automatically (the ERR-053 GMRES-truncation family is
closed by conformance, not per-site edits).

**The Hilbert adjoint comes free — Mode-12 closure by construction.**
:class:`CoupledOperator` implements only the Euclidean
:meth:`apply_transpose` — the transposed grid :math:`(A^{\mathsf T})_{ji}
= (A_{ij})^{\mathsf T}` — and carries a :class:`CoupledSpace` whose metric
methods dispatch member-wise. The metric adjoint :math:`A.H = G^{+}
A^{\mathsf T} G` is then realised ONCE by the existing
:class:`~orpheus.numerics.operator.AdjointOperator` wrapper. No adjoint
code lives in the block machinery — which is exactly what keeps the block
adjoint Mode-12-closed: a hand-rolled "Euclidean block ``.H``" that skips
the metric conjugation is the ERR-067 reopening (``vv-principles``
Mode 12). The block ``.H`` inheriting the member metrics is the campaign's
highest-value verification row.

.. _coupled-block-system-restriction:

The system restriction and extension by zero — one member's split pair
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A coupled field is a tuple of member fields, one per system. Two maps move
one member into and out of it: the **system restriction** selects member
:math:`i`, and its transpose, **extension by zero**, writes a member-:math:`i`
field into a coupled field whose other members are zero. This subsection is
the corpus's one definition of the pair; the conceptual view
(:ref:`architecture-conceptual-view`) and the S\ :sub:`N` pages that pose the
fission operator through it point here. It is defined beside the coupled
space it acts on, and not beside the boundary's trace restriction
(:ref:`bc-domain-narrowing`), because nothing about a system member is a
boundary: the two pairs share a structure (below), not a subject.

The subscripted :math:`r_i` is the code's spelling and the S\ :sub:`N`
pages' (:math:`r_{\rm bulk}` is :math:`r_1`, System A); a bare :math:`r` on
this page is the radial coordinate of the ψ½ two-point problem, and the
subscript is what tells them apart.

Every measurement in this subsection is taken on one fixture, the carrying
arm of ``tests/sn/operators/test_step2_posed_fission_anchors.py``: a
two-region sphere (edges ``np.linspace(0.01, 2.0, 9)``, mixtures A and B at
2 groups, ``Quadrature.gauss_legendre(8)``, isotropic scattering), whose
``problem.system.space`` couples System A, a ``FullField`` of 144 degrees
of freedom, with System B, the ψ½ ray, of 36. A draw is a coupled field
of standard-normal entries from ``np.random.default_rng(s)``, over 200
consecutive seeds per reading; the probes are
``scratch/_followups/sysres_probe*.py`` (untracked), run 2026-09-22.

**The definition.** Let :math:`V = V_1 \oplus \cdots \oplus V_N` be the
coupled space, carrying the block-diagonal metric
:math:`G = \operatorname{diag}(G_1, \dots, G_N)` that
:class:`~orpheus.numerics.coupled_system.CoupledSpace` applies member by
member, each :math:`G_j` being the metric of the member space :math:`V_j`.
For each member :math:`i`:

.. (vv-status rationale) Definitional: it states what the one class
   SystemRestrictionOperator applies (apply selects member i; apply_transpose
   is the scatter, extension by zero). Not a solver claim: no flux, no
   eigenvalue, no discretization error. Its defining behaviour is pinned by
   the foundation tests of TestSystemRestriction in
   tests/numerics/test_coupled_operator.py (member selection into the
   1-member coupling, extension with the dual-role zeros, the Euclidean
   pairing <r x, y> = <x, iota y>); foundation tests carry no verifies(...)
   marker, so the label is sentineled.
.. vv-status: coupled-block-system-restriction-pair documented

.. math::
   :label: coupled-block-system-restriction-pair

   r_i : V \to V_i,\quad r_i\,(x_1, \dots, x_N) = x_i
   \qquad\text{and}\qquad
   \iota_i = r_i^{\mathsf T} : V_i \to V,\quad
   \iota_i\,v = (0, \dots, 0, \underbrace{v}_{\text{slot } i}, 0, \dots, 0) .

.. implements:: coupled-block-system-restriction-pair
   :by: orpheus.numerics.coupled_system.SystemRestrictionOperator.apply

   **Implemented by** the restriction's ``apply``, which returns member
   :math:`i` of the coupled field as a one-member coupled field.

.. implements:: coupled-block-system-restriction-pair
   :by: orpheus.numerics.coupled_system.SystemRestrictionOperator.apply_transpose

   **Implemented by** the restriction's ``apply_transpose``, the extension
   by zero: the input written into slot :math:`i` of a fresh cotangent-side
   zero of the whole coupled space.

On the flat layout a coupled field packs (the member offsets
:attr:`CoupledSpace.system_slices
<orpheus.numerics.coupled_system.CoupledSpace.system_slices>`), :math:`r_i`
is the gather over the contiguous index block :math:`S_i` of member
:math:`i`, and :math:`\iota_i` is its scatter: the trace restriction pair
:eq:`bc-trace-restriction-pair` with :math:`S = S_i`. `[M]` on the fixture,
the production restriction and a
:class:`~orpheus.numerics.operator.TraceRestrictionOperator` built on
``system_slices[i]`` agree under ``np.array_equal`` in 200 of 200 draws, for
both members and in both directions, while the gather over the block
shifted by one row disagrees (``sysres_probe4.py``). What the coupled class
adds over the flat gather is typing: its ends are coupled spaces and its
values are member fields.

**The laws.** The pair obeys the laws of a *biproduct*, the direct sum of an
additive category, whose structure maps the category-theory literature
calls projections (here :math:`r_i`) and injections (here
:math:`\iota_i`); the corpus states the same laws for the face cochain at
:eq:`wavefront-cochain-biproduct-laws`. The fourth law is what the metric
adds:

.. (vv-status rationale) Identities of the pair, not a solver claim, and
   computed by nothing: no production code forms a sum of iota_i r_i or
   compares r_i.H with iota_i, so the label carries a declared
   no-implementation (below). Gated: the split law r_1 iota_1 = id and the
   idempotence of iota_1 r_1 on the FIRST member (the code's system=0), by
   the foundation tests of TestSystemRestriction
   (tests/numerics/test_coupled_operator.py), on toy members. Since
   2026-09-22 every law is gated on BOTH members (system=0 and system=1)
   of toy members with non-unit weights: the four laws by
   tests/numerics/test_coupled_operator.py::TestSystemRestrictionLaws
   (test_split_law, test_cross_law, test_completeness,
   test_metric_adjoint_is_extension_on_the_metric_range, each carrying
   verifies(coupled-block-system-restriction-laws)), the last including a
   member whose metric has a kernel, where r_i* = iota_i P_range(G_i)
   differs from iota_i; the two consequences by
   test_projector_is_idempotent and test_projector_is_G_orthogonal. Each
   gate reddens under an in-process mutation of the law it names (a
   wrong-member restriction, an extension writing into another member, a
   scaled extension, a dropped Riesz leg), and its docstring names the
   mutation class it is blind to. The 200-draw measurements on the
   production sphere carrier ([M] 2026-09-22, the fixture above) agree.
.. vv-status: coupled-block-system-restriction-laws documented

.. math::
   :label: coupled-block-system-restriction-laws

   \begin{aligned}
   r_i \circ \iota_i &= \mathrm{id}_{V_i},
   &\qquad
   r_j \circ \iota_i &= 0 \quad (j \neq i), \\
   \sum_{i=1}^{N} \iota_i \circ r_i &= \mathrm{id}_{V},
   &\qquad
   r_i^{*} &= \iota_i \circ P_{\operatorname{range}(G_i)} .
   \end{aligned}

.. no-implementation:: coupled-block-system-restriction-laws
   :kind: identity

   **Nothing implements this**: the laws are consequences of the definition
   :eq:`coupled-block-system-restriction-pair` and of the block-diagonal
   metric, and no production code evaluates any of them as an identity.
   Nothing forms the sum of the :math:`\iota_i r_i`, and the one production
   use daggers a composition that contains the restriction as a whole
   (`[M]` 2026-09-22, the hub's posed production's ``.H`` is an
   :class:`~orpheus.numerics.operator.AdjointOperator` whose inner operator
   is the composition, ``sysres_probe5.py``), so it applies
   :math:`\iota_1 = r_1^{\mathsf T}` between the Riesz legs that sandwich
   the composition's transpose and never forms :math:`r_1^{*}` on its own.

The first three are read off the definition. The restriction :math:`r_j`
reads slot :math:`j` of :math:`\iota_i v`, which holds :math:`v` when
:math:`j = i` and zero otherwise. The composite :math:`\iota_i r_i x` keeps
slot :math:`i` of :math:`x` and zeroes every other slot, so the sum over
:math:`i` recovers each slot exactly once. The fourth follows from the
adjoint's Riesz composition :eq:`spaces-adjoint-riesz-composition`,
:math:`r_i^{*} = \sharp_V \circ r_i^{\mathsf T} \circ \flat_{V_i}`,
applied to :math:`v \in V_i`:

.. math::

   r_i^{*}\,v
   \;=\; G^{+}\,\bigl(0, \dots, 0, G_i\, v, 0, \dots, 0\bigr)
   \;=\; \bigl(0, \dots, 0, G_i^{+} G_i\, v, 0, \dots, 0\bigr)
   \;=\; \iota_i\bigl(P_{\operatorname{range}(G_i)}\, v\bigr) .

The second step uses that the pseudo-inverse of a block-diagonal matrix is
the block-diagonal matrix of the blocks' pseudo-inverses,
:math:`G^{+} = \operatorname{diag}(G_1^{+}, \dots, G_N^{+})` (the four
Penrose conditions hold block by block, and products and transposes of
block-diagonal matrices are block-diagonal); the third is the Riesz round
trip :eq:`spaces-riesz-round-trip`. On a member whose metric is nonsingular
the projector is the identity, :math:`r_i^{*} = \iota_i`, and
:math:`r_i\, r_i^{*} = \mathrm{id}_{V_i}`: the restriction is a
coisometry, and **its adjoint is its extension by zero**. On a member
whose metric has a kernel (a trace block with tangential ordinates, for
instance) the two agree on the range of that metric and differ on its
kernel, exactly as the Riesz round trip does.

Two consequences follow. The composite :math:`P_i = \iota_i \circ r_i` is
idempotent, :math:`P_i^2 = \iota_i (r_i \iota_i) r_i = P_i`, and
self-adjoint in :math:`G`, because the metric is block-diagonal:

.. math::

   \langle P_i a, b\rangle_G
   \;=\; \langle a_i, b_i\rangle_{G_i}
   \;=\; \langle a, P_i b\rangle_G .

So :math:`P_i` is the :math:`G`-orthogonal projector onto the copy of
:math:`V_i` inside :math:`V`, the copies of different members are mutually
orthogonal, :math:`\langle P_i a, (I - P_i) b\rangle_G = 0`, and
:math:`\sum_i P_i = I` is the third law restated.

`[M]` on the fixture, 200 draws per row, each law read on the flat values
of the coupled fields. The first four rows were read on a source-role field
and again on a flux-role field, with the same result; each control is a
deliberately wrong reading that must fail
(``sysres_probe.py``, ``sysres_probe6.py``, ``sysres_probe_controls.py``):

.. list-table:: The laws of the system restriction, measured on the fixture
   :header-rows: 1
   :widths: 36 64

   * - law
     - reading
   * - :math:`r_i \circ \iota_i = \mathrm{id}`
     - bitwise equal in 200 of 200 draws, for :math:`i = 1` and
       :math:`i = 2`, on a source-role member and on a flux-role member;
       control: against a doubled member, 0 of 200
   * - :math:`r_j \circ \iota_i = 0`, :math:`j \neq i`
     - exact zeros in 200 of 200 draws, both off-diagonal pairs; control:
       :math:`r_1 \circ \iota_1` reads all-zero in 0 of 200
   * - :math:`\sum_i \iota_i \circ r_i = \mathrm{id}`
     - bitwise equal in 200 of 200 draws; control: the sum with one term
       only, 0 of 200
   * - :math:`P_i^2 = P_i`
     - bitwise equal in 200 of 200 draws, both members
   * - :math:`\langle P_i a, (I - P_i) b\rangle_G = 0`
     - largest absolute value ``0.0`` over 200 draws and both members, with
       :math:`G` materialised column by column from
       :meth:`CoupledSpace.apply_metric
       <orpheus.numerics.coupled_system.CoupledSpace.apply_metric>` and its
       off-member blocks exactly zero
   * - :math:`r_i^{*} = \iota_i`
     - both member metrics nonsingular (diagonal; smallest eigenvalue
       :math:`7.3\times10^{-3}` for System A and :math:`7.3\times10^{-2}`
       for System B); ``r.H`` against ``r.apply_transpose`` differs by at
       most :math:`2.6\times10^{-16}` relative, and is bitwise equal in 0
       of 200 draws for System A and 7 of 200 for System B: the round-off
       of :math:`G_i^{+} G_i v`, not a different map

**Why extension by zero, and why one class.** The first law makes
:math:`\iota_i` a *section* of :math:`r_i`, a right inverse, in the sense
:ref:`spaces-collapse-pair-naming` fixes. It is not the only one: every map
:math:`s(v) = (A_1 v, \dots, v, \dots, A_N v)`, with :math:`v` in slot
:math:`i` and arbitrary linear :math:`A_j : V_i \to V_j` in the other
slots, satisfies :math:`r_i \circ s = \mathrm{id}`. The metric picks one out.
Take :math:`a` with only member :math:`i` nonzero in
:math:`\langle s\,r_i a, b\rangle_G = \langle a, s\,r_i b\rangle_G`: the left
side is :math:`\langle a_i, b_i\rangle_{G_i} + \sum_{j \neq i} \langle A_j
a_i, b_j\rangle_{G_j}` and the right side is :math:`\langle a_i,
b_i\rangle_{G_i}`, so :math:`s \circ r_i` is :math:`G`-self-adjoint exactly
when :math:`G_j A_j = 0` for every :math:`j \neq i` (the converse is the same
computation read backwards). With nonsingular member
metrics that forces :math:`A_j = 0`: extension by zero is the one section
whose projector is orthogonal, and it is the restriction's adjoint.

That is why the family carries one class and not two. Because the section
is the transpose, :class:`~orpheus.numerics.coupled_system.SystemRestrictionOperator`
serves both arrows, ``apply`` the restriction and ``apply_transpose`` the
extension by zero. The axis collapse pair (:ref:`spaces-collapse-pair`) is
the contrasting kind: its retraction is not a coisometry,
:math:`R\,R^{\dagger} = (\Sigma w)\,\mathrm{id}`, so its section
:math:`E = R^{\dagger}/\Sigma w` differs from the adjoint by that scalar and
is a class of its own,
:class:`~orpheus.numerics.operator.AxisSectionOperator`. The trace
restriction is the system restriction's kind, a gather whose scatter is its
transpose.

"Embedding" is not an operator name in this corpus
(:ref:`spaces-collapse-pair-naming`), and the argument recorded there
applies unchanged: every section :math:`s` of :math:`r_i` is injective, so
"embedding" names a property all of them share and cannot pick out the one
that is the restriction's adjoint. The pair is named by what each arrow
does, restriction and extension by zero.

**The class.** ``SystemRestrictionOperator(space, system=i)`` realizes the
pair on a :class:`~orpheus.numerics.coupled_system.CoupledSpace` (the code
counts members from 0, so the code's ``system=0`` is :math:`r_1`):

- **Born bound, onto a one-member coupling.** Its domain is the whole
  coupled space; its codomain is the one-member coupled space
  ``CoupledSpace.from_systems([V_i])`` built at construction, not the
  member space :math:`V_i` itself, so :math:`r_i x` is a one-member
  :class:`~orpheus.numerics.coupled_system.CoupledField` and
  :math:`\iota_i` consumes one. `[M]` on the fixture, for both members, the
  codomain's single member is the parent's member-space object (an
  ``is`` comparison) while the codomain compares unequal to that member
  space, its name being ``coupled(<member name>)``. The codomain's zero
  factories are derived from the parent's.
- **Refusals.** ``apply`` refuses an argument that is not a coupled field
  (``TypeError``) and a coupled field of the wrong arity (``ValueError``);
  a member index outside :math:`0, \dots, N-1` is refused at construction
  (``ValueError``, "out of range").
- **The extension is the cotangent map.** ``apply_transpose`` mints the
  zeros it extends with through the space's cotangent seam,
  :meth:`CoupledSpace.dual_zeros
  <orpheus.numerics.coupled_system.CoupledSpace.dual_zeros>`, so they carry
  the source role, the dual of the flux role: a transpose consumes the
  codomain's cotangent and emits the domain's. `[M]` on the fixture,
  extending System A leaves System B as a
  ``RadialCharacteristicInteriorSourceSink`` and
  ``RadialCharacteristicBoundarySourceSink`` pair. On a coupled space built
  without a dual factory, ``apply_transpose`` raises the space's
  ``RuntimeError`` naming the missing "DUAL zero-element factory".
- **The adjoint is not a method of its own.** ``r.H`` is the generic
  :class:`~orpheus.numerics.operator.AdjointOperator`,
  :math:`\sharp \circ \text{dual} \circ \flat` (:ref:`spaces-riesz-legs`),
  with ``apply_transpose`` as its middle factor; the fourth law is its
  theorem, and the table above its measurement.

.. warning::

   **The extension by zero is role-consistent only on a source-role
   input.** Its other members are always source-role zeros, and the class
   of those zeros is load-bearing, not incidental: the class docstring
   records (`[M]` 2026-08-22) that extending with the flux-role zeros put a
   flux-classed ray zero into the daggered fission chain, where the
   cross-class arithmetic gate refused the sum. `[M]` on the fixture, on a
   source-role coupled field :math:`y` everything composes, and the typed
   sum :math:`\iota_1 r_1 y + \iota_2 r_2 y` equals :math:`y`. On a
   flux-role coupled field :math:`x`, :math:`\iota_i r_i x` is a mixed-role
   composite for both members (member :math:`i` flux-role, the other a
   source-role zero), and so is ``r.H`` applied to :math:`r_1 x`; adding
   :math:`x` to any of the three, or forming the typed sum
   :math:`\iota_1 r_1 x + \iota_2 r_2 x`, is refused with a ``TypeError``
   by the same gate. The laws hold on the flat values of a flux-role field
   (the table above), but no extension by zero whose zeros carry the flux
   role is built. Whether a Riesz leg should change the role class, so
   that ``r.H`` on a flux-role member returns an all-flux composite, is
   open as #496. The one production consumer never meets the flux-role
   case: the daggered posed production applies the composition's transpose,
   :math:`\iota_1 \circ \text{stack}^{\mathsf T}`, between its Riesz
   legs, the stack's transpose already emits a source-role System-A member,
   and `[M]` the hub's posed production's ``.H`` returns an all-source-role
   coupled field on either role of input.

**Where it is used.** `[M]` 2026-09-22, an AST census of call sites spelled
``SystemRestrictionOperator(...)`` finds 1 in ``orpheus/`` (352 files
parsed) and 3 in ``tests/`` (612 files), all 3 in
``tests/numerics/test_coupled_operator.py``; the same census finds 6
``CoupledOperator(...)`` call sites in ``orpheus/``. The production site is
:func:`~orpheus.sn.coupled_system.build_within_group_system`, which poses
the fission operator on a carrying mesh as the composition
:eq:`sn-posed-production-carrying`, :math:`[\,[F], [E_F]\,] \circ r_1`. The
reason is a physical fact recorded where the composition is read: fission
annihilates the ray system, and an annihilated input is a restriction, not a
zero block (:ref:`sn-adjoint-coupled-posing`). The restriction and that
posing landed together (S4-amendment A2, 2026-08-22, ``6fc247fb``),
replacing a :math:`2\times2` grid whose zero right-hand column carried two
hand-written zero closures; the transpose of the composition,
:math:`\iota_1 \circ [\,F^{\mathsf T}, E_F^{\mathsf T}\,]`, is where the
extension by zero does its work.

**What the gates pin, and what is only measured.** The 8 foundation tests
of ``TestSystemRestriction`` in ``tests/numerics/test_coupled_operator.py``
pin, on toy members: member selection into the one-member codomain; the
split law :math:`r \circ \iota = \mathrm{id}`; extension with the dual-role
zeros; the Euclidean pairing :math:`\langle r x, y\rangle = \langle x,
\iota y\rangle`; the idempotence of :math:`\iota \circ r`; the two refusals
(an unwired dual seam, a member index out of range); and
``is_adjointable``. Every construction site among them restricts onto
the first member, or refuses a third. The rest was measured and not gated
until 2026-09-22; since then
``tests/numerics/test_coupled_operator.py::TestSystemRestrictionLaws``
gates, on both members, the cross law,
the completeness sum, the metric adjoint in its full form
:math:`r_i^{*} = \iota_i \circ P_{\operatorname{range}(G_i)}` (on a member
with a singular metric it is not :math:`\iota_i`), the idempotence and the
:math:`G`-orthogonality of :math:`P_i`, each red under a mutation of the law
it names.

The two-system role lattice — :class:`SystemRole`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`~orpheus.numerics.operator.SystemRole` tags which of the two
systems an operator's action maps between —
:attr:`A` (within System A), :attr:`B` (within System B: the self-block
:math:`A_{BB}` and the ray boundary :math:`B_b`), or :attr:`COUPLED` (an
off-diagonal, or the assembled 2×2: :math:`A_{AB}`, :math:`A_{BA}`, and
the grid). Reading each role as the *set* of systems its action touches
(:math:`A = \{A\}`, :math:`B = \{B\}`, :math:`\text{COUPLED} = \{A,B\}`),
a sum touches the union: "same role stays, anything different becomes
COUPLED". This is the two-system analogue of the bulk↔boundary
:class:`~orpheus.numerics.operator.BlockRole` join, with ``COUPLED``
playing the top-of-lattice role ``FULL`` plays there — a deliberate twin
kept while only two role axes exist (the collapse trigger is a third
parallel axis, e.g. a DSA/multiphysics role). Model-generic leaves
(``C``/``S``/``F``, every diffusion/CP/MoC operator) leave ``system_role``
at ``None`` — the conservative "not part of the ψ½ augmentation" reading —
so a System-A composite must be **stamped** ``SystemRole.A`` explicitly at
its composition site (the model-generic members' honest ``None`` would
otherwise poison the join to ``None``).

The one production spelling — ``build_within_group_system``
-----------------------------------------------------------

The joint system has exactly ONE construction site:
:func:`~orpheus.sn.coupled_system.build_within_group_system`, which
returns the frozen :class:`~orpheus.sn.coupled_system.WithinGroupSystem`
record — **four** members: the loss grid, the
:class:`~orpheus.numerics.coupled_system.CoupledSpace` it is typed
against, the bound LEAVES the grid is the signed sum of
(:class:`~orpheus.sn.coupled_system.SNLossFactors`, by role), and — since
2026-09-13 — the posed production
:attr:`~orpheus.sn.coupled_system.WithinGroupSystem.production`, which is
:math:`F` on that same carrier so that the pair :math:`(A, F)` is a
**pencil** on one space (:ref:`the-operator-pencil`).  All four are built
from the SAME piece objects (one ``L+C``, one ``S``, one ``B_a``, one
``B_b``, …).

⭐ **And that builder now has exactly one CALLER**, which is a stronger
statement than "one construction site".  Since 2026-09-13 it is invoked
only by the :func:`~functools.cached_property`
:attr:`SNProblem.system <orpheus.sn.problem.SNProblem.system>`, so
a Problem's joint system is built once and every consumer reads THAT
object: ``[M]`` a counting spy over one forward eigenvalue solve read the
outer-iteration count before the change and reads **1** after it
(:ref:`sn-the-problem-poses-its-pencil`).  The two ``scattering_op=`` /
``n2n_op=`` keyword arguments that let a caller inject its own leaves —
the seam through which the Problem's :math:`A` was a function of
something other than the Problem's data — are **deleted** in the same
change.

The **splitting** :math:`A = M - N` (Hackbusch 2016 §11 —
block partitionings; a *splitting*, **not** a *regular* splitting in
Varga's sense, see :ref:`sn-boundary-gs-not-regular`) is *not* one of
those members: a posing does not choose which leaf is inverted, so since
2026-09-13 it is a Strategy VALUE labelled from the leaves
(:class:`~orpheus.sn.splitting.Splitting`,
:ref:`sn-splitting-is-a-strategy-value`), and the grids below are what
that value derives on a carrying mesh.
This builder is what retired the former ``_within_group_triple`` /
``_lagged_gains`` construction pair. The grid and its
:class:`CoupledSpace` are emitted **together**, aligned by construction
(RULING P1: a mismatched operator/space pairing is unconstructable).

The **loss grid** carries the loss-sign convention, and the two
off-diagonals differ in sign — a trap worth spelling out:

.. math::
   :label: coupled-loss-grid

   A \;=\;
   \begin{bmatrix}
     L + C - S - N_{2n} - B_a & +\,\text{Seeding} \\[2pt]
     -\,\text{Emission} & A_{BB} - B_b
   \end{bmatrix}

The :math:`(A,B)` seed is **positive** (its ``apply`` already emits the
seed's term of the fused bulk residual row — the loss sign is *internal*
to the operator, matching the in-sweep placement it wraps); the
:math:`(B,A)` emission is **negated** (a
:class:`~orpheus.numerics.operator.ScaledOperator` — the emission is a
gain, so the ray equation :math:`(A_{BB} - B_b)\psi_B -
\text{Emission}\,\psi_A = q_B` carries the minus).

.. vv-status: coupled-loss-grid documented

**The splitting.** The SI/Krylov drivers do not consume the loss grid
raw — they consume its splitting :math:`A = M - N`,

.. math::
   :label: coupled-mn-splitting

   M = \begin{bmatrix} L+C & +\,\text{Seeding} \\ \mathbf 0 & A_{BB}
       \end{bmatrix},
   \qquad
   N = \begin{bmatrix} S + N_{2n} + B_a & \mathbf 0 \\
       +\,\text{Emission} & B_b
       \end{bmatrix},

with :math:`M` the **sweepable part** inverted every step and :math:`N`
the lagged coupling gains (all signs positive — gains on the rhs
:math:`\text{rhs} = q + N\psi`; the loss grid's :math:`-\text{Emission}`,
:math:`-B_b` minus signs are the :math:`M - N` complement, not a
contradiction). Note :math:`M`'s :math:`(B,B)` entry is the **bare** march
:math:`A_{BB}` while :math:`N`'s is :math:`B_b`, so
:math:`A(B,B) = A_{BB} - B_b` recovers the loss self-block. The drivers
iterate :math:`\psi \leftarrow M^{-1}(q + N\psi)` (SI) or GMRES on
:math:`(M-N)\psi = q`.

.. vv-status: coupled-mn-splitting documented

**Presence is structural (R12a).** System B exists **only** on a
seed-carrying mesh — the 1-D sphere and, since the Q5.6.3 admission flip,
*every admitted cylinder*, whose folded rule carries on every
:math:`\mu`-level: the mesh's ray spaces are ``None``
:math:`\iff` non-carrying, and the ``RadialCharacteristic*`` block
constructors *refuse* seedless meshes. So a seed-carrying mesh builds the
2×2; a seedless mesh — a Cartesian chart, or a curvilinear rule with no
carrying level (a class cylindrical admission refuses since Q5.6.3; a
:math:`\mu = -1`-noded sphere rule reaches it, the sphere arm having no
admission gate) — builds the 1×1
:math:`[[A_{AA}]]` and the splitting degrades to the bare :math:`(L+C,\
(S, N_{2n}, B_a))` the seedless driver paths consume zero-touch. "Applying a
System-B block on a non-carrying mesh" is not a runtime branch — it is an
object that **does not exist**. Since the B.2d eviction ``FullField`` is a
pure 2-block composite (:math:`\psi_A` cannot carry ray state at all), the
old dead-slot double-count hazard dissolved structurally: the coupled flat
dimension is the honest two-system sum, no dead padding.

The block-triangular solve and the materialised-LU EXTRACT
----------------------------------------------------------

On a carrying mesh :math:`M` is an **honest upper-triangular**
:class:`CoupledOperator` grid, so its ``solve`` is matrix-free block
back-substitution — and the substitution order is exactly the curvilinear
sweep order:

.. math::
   :label: coupled-block-substitution

   \psi_B \;=\; A_{BB}^{-1}\,q_B ,
   \qquad
   \psi_A \;=\; (L+C)^{-1}\bigl(q_A - \text{Seeding}\,\psi_B\bigr).

System B's march runs first (its exact direct inverse — one pass, no
seed), then System A's bulk sweep runs on its **ray-decoupled** channel
:math:`q_A - \text{Seeding}\,\psi_B`. This is one body for all four
orientation × transpose combinations: the transpose of a triangular grid
is triangular the other way (:math:`(A^{\mathsf T})_{ij} =
(A_{ji})^{\mathsf T}`), so the same substitution runs with the visit
order flipped, guarded per block by
:class:`~orpheus.numerics.operator.MissingAdjoint`. The inverse **object**
is :class:`~orpheus.numerics.coupled_system.CoupledSubstitutionOperator`
(the taxonomy §12 wrap-delegate — its ``apply`` delegates to the grid's
``solve``, ``initial_guess`` accepted and dropped since a direct
substitution has nothing to seed).

.. vv-status: coupled-block-substitution documented

**The materialise/LU EXTRACT (step 5a).** The *full* 2×2 loss grid is not
triangular (:math:`A_{BA} \neq 0`), so it takes the second route:
:meth:`CoupledOperator.inverse` materialises the assembled matrix and
LU-factors it via
:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
(one eager LU, then back-substitution per solve). This EXTRACT is
**principled-equivalent, not bit-identical** to the production
splitting-iteration — the matrix reduction tree differs — and a naive
extraction returns O(1) garbage (the sweep treats inflow/seed rows as
*given data*, so the row-contract must be preserved). It is the oracle the
swap-law gates ride; production solves stay the splitting iteration on a
:class:`~orpheus.sn.splitting.Splitting` value minted from the record's
factors. The **iterative** splitting solve
(block-Jacobi / block-Gauss-Seidel over :math:`A = M - N`) deliberately
stays with the drivers — convergence is spectral
(:math:`\rho(M^{-1}N) < 1`), never a structural capability.

Why the :math:`(B,A)` block is ``None`` in :math:`M` — the Schur/lag argument
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ray↔bulk coupling splits cleanly by *which direction iterates*. The
seed :math:`A_{AB}` sits **inside** :math:`M` (the :math:`(A,B)` slot),
because within :math:`(L+C)` the ray→bulk coupling is block-triangular and
solved **directly** in one pass — measured :math:`A_{sb} = 0` exact (no
bulk→seed feedback inside :math:`L+C`), :math:`A_{bs} = 7.5` (the M-M
seed feed). The emission :math:`A_{BA}` sits **outside** the resolvent, in
:math:`N` — it is the bulk *scattering* gain's ray arm, and ONLY the
scattering coupling iterates: measured :math:`S_{sb} = 0.183` (bulk→ray
source), :math:`S_{bs} = 0` exact (ψ½ carries zero moment weight), with
outer spectral radius :math:`\rho(M^{-1}N) = 0.371 < c` (Adams–Larsen).
Placing :math:`A_{BA}` in :math:`M` would fold a lagged gain into the
one-pass resolvent — it belongs on the rhs, lagged. So the :math:`(B,A)`
slot of :math:`M` is structurally zero and the resolvent stays
block-**triangular** (direct); the coupling that genuinely iterates rides
:math:`N`. "Deciding where the block goes IS operator algebra": the same
four posed operators would let :math:`A_{BA}` sit folded, explicit, or in
a DSA preconditioner — a composition choice the machinery supports.

The ρ-honest stop and the lag-death certificate
-----------------------------------------------

Splitting the emission into a lagged gain creates a verification hazard:
an iteration whose fixed point does **not** solve the equation (a stale or
lagged block inside :math:`M` — the #282 class) can still report
"converged" to any :math:`\lVert\Delta\psi\rVert`-family stop. Step 5c
closes this in two pieces.

**The running stop is the ρ-honest equation residual** (a deliberate
re-interpretation of ``tol`` — :class:`~orpheus.numerics.iteration.SourceIteration`):

.. math::
   :label: coupled-free-identity-residual

   r_{n} \;=\; A\,\psi_{n} - q_{\rm ext}
         \;=\; \mathrm{rhs}_{n-1} - \mathrm{rhs}_{n}
         \;=\; \textstyle\sum_i g_i\,(\psi_{n-1} - \psi_{n}),
   \qquad
   {\rm res}_n = \frac{\lVert r_n\rVert_2}
                      {\max(\lVert q_{\rm ext}\rVert_2,\ 10^{-30})} .

The **free-identity** spelling :math:`\mathrm{rhs}_{n-1} -
\mathrm{rhs}_{n}` (retained from the loop's own bookkeeping — zero
marginal cost, checked before the next apply) is exact when the step
operator is an exact inverse of :math:`M`: then :math:`M\psi_n =
\mathrm{rhs}_{n-1}` and :math:`A\psi_n - q = N(\psi_{n-1} - \psi_n)`. The
residual is preferred over :math:`\lVert\Delta\psi\rVert` because the
increment understates the true error by :math:`1/(1-\rho)` (Banach) AND is
blind to the lag-death class; the residual claim is
contraction-rate-independent.

.. vv-status: coupled-free-identity-residual documented

**The certificate closes the exact-**:math:`M` **assumption.** The
free identity itself *assumes* exact :math:`M` — the hole the driver-level
**convergence certificate**
:class:`~orpheus.sn.solver.ConvergenceCertificateError` closes. Once, at a
*claimed* exit, :func:`~orpheus.sn.solver.evaluate_residual` measures the
true :math:`r = A\psi - q` through a real forward apply — the only
measurement an in-:math:`M` lag cannot fool — and raises if the defect
exceeds :math:`10\times\text{tol}`. It is wired on every full-angular arm
(the coupled sphere, the seedless un-windowed SI, both Krylov paths). The
headline gate is a stale-zero-:math:`\psi_B` in-:math:`M` lag mutation: the
running stop reports convergent while the certificate raises
``match="lag-death"`` — with control legs on both sides, the classifier's
asymmetry proof (#282 measured the cold-residual defect at ~5e5). The
typed residual carries System B as its **own** member
(:class:`~orpheus.transport.residuals.radial_characteristic_interior_residual.RadialCharacteristicInteriorResidual`
⊕ its boundary sibling), so a System-A-only residual cannot silently drop
a wrong seed row — the Mode-12 (b) closure. This is an additive diagnostic
+ certificate, NOT in the convergence path (see :ref:`affine-typed-residual`).

The swap law on the grid
------------------------

The inverse-as-operator taxonomy's swap law :math:`A.H.\text{inverse}()
\equiv A.\text{inverse}().H` extends to the grid: the joint adjoint IS the
grid's transposed substitution. :class:`CoupledSubstitutionOperator`
advertises ``is_adjointable`` iff every coupling block transposes AND every
diagonal carries the direct ``solve_transpose`` verb (the #280 two-factor
discipline per block), and its :meth:`apply_transpose` delegates to
:meth:`CoupledOperator.solve_transpose` (the transposed block
substitution). The step-6 collapse sharpened the single-system predicate
to match: :attr:`SweepOperator.is_adjointable
<orpheus.sn.operators.sweep_operator.SweepOperator.is_adjointable>` is now
the two-factor ``isinstance(inner, StreamingCollisionOperator) and
inner.is_adjointable`` — the B.2d carrying-mesh *third* factor retired,
because with no ψ½ legs anywhere a bare :math:`(L+C)` on a carrying mesh IS
unambiguously the ray-decoupled :math:`(A,A)` block (the type says what you
hold). The joint adjoint's one home is the grid's transposed substitution;
there is no adjoint SI driver in production, so the :math:`A_{BA}`/:math:`A_{AB}`
pullbacks are reached only by the ``.H`` reciprocity gates — which is why
they are realised now and pinned by **nonzero-seed-cotangent** gates (a
present-zero seed would hide a lost pullback).

The one-channel collapse (step 6) — the walk IS the :math:`(A,A)` block
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With the four blocks posed and the solve routed through the grid, the
walk's fused ψ½ joint-leg channel became dead production code and was
**retired** (step 6, net −803 lines). The two forward/transpose kwarg
pairs, the three-function presence-guard family, the walk's in-solve
System-B engines, and the operator-free ``transport_sweep`` wrapper are
all gone. Every walk surface is now the ray-decoupled :math:`(A,A)`
diagonal block, on **every** mesh, with no seed parameters: the matvec
substitutes a zero seed into the Morel–Montry thread (bit-identical to the
retired dead-slot arithmetic — the closure reads zeros), and the transpose
**discards** the thread cotangent (a fixed zero input's cotangent
propagates nowhere; the :math:`\text{Seeding}^{\mathsf T}` pullback is the
explicit grid block). The eigenvalue finalize re-routes through the SAME
:func:`~orpheus.sn.coupled_system.build_within_group_system` posing every
driver consumes — and since #448 it does not even rebuild it: it reads the
:class:`~orpheus.sn.solver.InnerSolve` record the last inner solve left
behind, so it holds the driver's own operator instance and the driver's own
gains (:ref:`sn-finalize-one-step`).  On the SI arm that operator is the
**un-windowed** forward :math:`M`, and it now *is* a named attribute:
``inner.splitting.implicit``.

⭐ That sentence read *"which is not always* ``.implicit_operator`` *— the
boundary-Gauss-Seidel schedule splits it — so the precise statement is the
splitting the inner solve drove, not a named attribute"* until 2026-09-13,
and the caveat was the whole diagnosis: the record's named attribute could
disagree with what the driver ran, because the Gauss-Seidel arm re-derived
a second splitting behind it.  The consumers campaign's step 2 removed the
disagreement at the source — the splitting the inner solve drove is a
value the inner solve *carries*
(:ref:`sn-splitting-is-a-strategy-value`) — so the precise statement and
the named attribute are now the same thing.

The mesh remains the single authority on presence;
what changed is that nothing *checks* against it anymore — the type system
carries the biconditional. The narrative of the walk's own view of this
collapse lives in :doc:`/theory/methods/sn/loss_representation`.

Numerical evidence
------------------

The block architecture is a **principled-equivalent** re-association of
the fused implementation — same operator algebra, different reduction
tree — so the sphere re-baselines are FP-grain, sphere-only, and the
seedless geometries stay bit-identical (the strongest leak tripwire).

.. list-table:: Measured equivalence (the fused → block-sum re-association class)
   :header-rows: 1
   :widths: 40 30 30

   * - Check
     - Result
     - What it pins
   * - :math:`A_{BB}` forward-substitution
       :math:`\lVert\text{solve}\circ\text{apply}(\psi)-\psi\rVert`
     - :math:`3.5\text{e-}16`
     - the march is the exact direct inverse (#284, :math:`\rho=0`)
   * - block-triangular certificate :math:`A_{sb}`
     - :math:`0` exact
     - no bulk→seed feedback inside :math:`L+C` (direct one-pass)
   * - EXTRACT vs dense-LU of assembled :math:`(L+C)`
     - :math:`5.5\text{e-}16`
     - the materialise/LU route is the principled oracle
   * - scatter LIFT vs the old monolithic ``S.apply``
     - :math:`1.2\text{e-}15`
     - :math:`S_{\rm bulk} + A_{BA}` ≡ the pre-lift composite
   * - finalize equivalence — slab / cylinder, every member
     - **bitwise**
     - seedless geometries are untouched (leak tripwire)
   * - finalize equivalence — sphere :math:`k`, scalar, both ψ½ members
     - **bitwise**
     - :math:`A_{BB}`'s march ≡ the retired in-solve march exactly
   * - finalize equivalence — sphere bulk + trace movers
     - rel-max :math:`8.4\text{e-}16`
     - the M-substitution vs fused-interleave re-association only

.. note::

   ⚠ **The three "finalize equivalence" rows are step-6-relative and are
   not a live property of the finalize.**  They measure the fused →
   block-sum re-association *of a reconstruction that no longer exists*:
   #448 (2026-09-06) replaced the finalize's body outright — it is now
   one step of the driven map rather than a fresh solve of a hand-built
   source (:ref:`sn-finalize-one-step`,
   :doc:`ERR-083 </theory/verification/error_catalog>`).  Re-running these
   checks today would compare a different object, and at
   ``scattering_order ≥ 1`` the pre-#448 finalize was ``[M]`` 8.776e-02
   from the fixed point, so "bitwise" here means *bitwise across the
   step-6 carve*, never *correct*.  The numbers stay because they are the
   evidence step 6 was leak-free; the caption is what needed the date.

The outer spectral radius :math:`\rho(M^{-1}N) = 0.371` and the measured
block algebra (:math:`A_{bs} = 7.5`, :math:`A_{ss} = 5.0`, :math:`S_{sb} =
0.183`, :math:`S_{bs} = 0`) confirm the direct/iterative split the API
reflects. The full within-group wall (tests/sn + tests/numerics, not-slow
serial) is 3080/0 through the step-6 collapse.

The DSA seam — issue #2 as a future consumer
--------------------------------------------

The coupled-block machinery is not a primitive-without-a-product beyond
the ψ½ instance: **consistent DSA (Issue #2)** is its next consumer. A
diffusion-synthetic-acceleration step poses the SN transport system and a
low-order diffusion correction as a **coupled system** — an N-system
:class:`CoupledOperator` whose diagonal blocks are the transport loss and
the diffusion loss operator :math:`A_{\rm diff} = L + C - S - B` (which
already exists, :mod:`orpheus.diffusion.operators`, #290 P4 — four
terms, carrying no separate :math:`N_{2n}`, because diffusion's
:math:`S` is itself the ``OperatorSum``
``IsotropicScattering + IsotropicN2N``, i.e. the SAME two energy leaves
the S\ :sub:`N` grid keeps apart, bundled at the *composition* site
where a bundling is a legitimate local choice; see the note at
:eq:`coupled-ba-emission`) and whose
off-diagonals are the restriction :math:`R` (transport residual → coarse)
and prolongation :math:`P` (coarse correction → fine) operators. Those
:math:`R`/:math:`P` operators **do not exist yet** — they are the future
``SystemRole`` members Issue #2 will add; the seam is *named*, not built.
When they land, DSA consumes the same three modes (the block matvec, the
``assemble`` route the diffusion resolvent already rides via
:ref:`operator-algebra-assembly-axis`, and the block ``solve``), and it
consumes the transport residual :func:`~orpheus.sn.solver.evaluate_residual`
already types — DSA computes that residual, the diffusion solve corrects
it, and the correction :math:`\to 0` at convergence (so the accelerator is
correctness-safe by construction). The Wave-T DSA row
(:ref:`tensor-network-decomposition`) tracks the same seam from the
loss-leaf side.
