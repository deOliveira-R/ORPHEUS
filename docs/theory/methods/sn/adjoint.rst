.. _sn-adjoint:

Adjoint transport: the dual operators
=====================================

This chapter is the S\ :sub:`N` book's adjoint rung — the transposes
that the perturbation-theory, detector-sensitivity, and
adjoint-weighted-homogenisation chains consume.  It is deliberately
thin today, and honest about why: of the adjoint chain's three
layers, two are landed and documented (one of them elsewhere), and
the third is in flight.

* The **walk adjoint** — the loss composite's transpose
  :math:`(L+C)^{\mathsf T} = L^{\mathsf T} + C` — is machinery of the
  loss *representations* and lives with them: the orientation axis,
  its swap law, and the deferral ledger
  (:ref:`loss-rep-orientation-two-frames`), realized as the Wave-O
  analytic reverse-direction matvec.  This chapter points, it does
  not re-derive.
* The **scattering adjoint** :math:`S^{\mathsf T}` — the section
  below, the #276 P3 record: free by frame conjugation, with no
  per-geometry derivation to verify.
* The **daggered posing and the adjoint flux** :math:`\psi^*`
  solving :math:`(L+C-S)^{\mathsf T}\psi^* = q^*` (the boundary gain
  :math:`B` transposes alongside and is elided in this shorthand —
  the full reflective-BC statement is
  :math:`(L+C-S-B)^{\mathsf T}`) are **not yet
  landed** (`#276 <https://github.com/deOliveira-R/ORPHEUS/issues/276>`_
  A4/A5); the :math:`\varphi^*` consumers — adjoint-weighted
  homogenisation (frame-machinery P6,
  `#51 <https://github.com/deOliveira-R/ORPHEUS/issues/51>`_) and the
  response posing
  (`#281 <https://github.com/deOliveira-R/ORPHEUS/issues/281>`_) —
  unblock when they land, and this chapter grows with them.

.. admonition:: Key Facts
   :class: tip

   * :math:`S^{\mathsf T}` is assembled from **leaf transposes** of
     the frame conjugation,
     :math:`(R \circ (\Lambda + N_{2n}) \circ M)^{\mathsf T}
     = M^{\mathsf T} \circ (\Lambda + N_{2n})^{\mathsf T} \circ
     R^{\mathsf T}` — no per-geometry derivation (#276 P3, closes
     #118), reciprocity-pinned against the structurally *independent*
     forward fast-path.
   * It is the **Euclidean transpose**, not the Hilbert adjoint: the
     angular-Gram weighting :math:`S^{\dagger} = G^{-1}S^{\mathsf T}G`
     is the ``.H`` wrapper's job.  The campaign name "S†" is
     colloquial.
   * The **walk adjoint is not documented here** — it is the
     orientation axis of the loss representations
     (:ref:`loss-rep-orientation-two-frames`); conflating the three
     transposes (Euclidean, Hilbert, walk-orientation) is the #1
     landmine of that carve.
   * The forward-adjoint **asymmetry is principled**: the forward
     source keeps the scalar fast-path for SI-sweep performance; the
     adjoint — not the hot path — rides the validated frame form.
     That structural difference is what makes the reciprocity gate a
     genuine cross-check.

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

References
==========

* Lewis, E. E., & Miller, W. F. (1984).  *Computational Methods of
  Neutron Transport.*  §10 (adjoint transport, reciprocity, and
  perturbation theory).
