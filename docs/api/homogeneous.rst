Homogeneous Infinite-Medium Solver
====================================

The :mod:`orpheus.homogeneous` package solves the multi-group
eigenvalue problem in an infinite homogeneous medium — the
simplest reactor physics configuration and the foundation for
every spatial solver in ORPHEUS. Because the flux is uniform and
all streaming terms vanish, the transport equation collapses to a
single :math:`G \times G` dense eigenvalue problem solved directly
(no iteration), with an exact analytical structure that makes it the
go-to harness for L0 / L1 verification of cross-section libraries and
scattering-matrix conventions.

.. contents::
   :local:
   :depth: 2

.. seealso::

   :ref:`theory-homogeneous` — full derivation, scattering
   convention, and worked examples.


Eigenvalue Problem
------------------

With no spatial dependence, the multi-group transport equation
reduces to

.. math::

   \underbrace{\bigl(\operatorname{diag}(\Sigma_t)
   - \Sigma_{s0}^{\mathsf T}
   - 2\,\Sigma_{2}^{\mathsf T}\bigr)}_{\mathbf{A}}\,\phi
   \;=\; \frac{1}{k_\infty}\,
   \underbrace{\bigl(\chi \otimes \nu\Sigma_f\bigr)}_{\mathbf{F}}\,\phi,

where :math:`\Sigma_{s0}` is the :math:`P_0` (isotropic) scattering
matrix, :math:`\Sigma_2` is the (n,2n) cross-section matrix (stored
separately because each collision produces two neutrons), and
:math:`\chi` is the prompt fission spectrum. The :math:`(n,2n)`
reaction is a **loss-side multiplicity-2 transfer**: the
:math:`-2\,\Sigma_2^{\mathsf T}` term in the loss matrix
:math:`\mathbf{A}` removes the incident neutron and redistributes the
two emitted neutrons by the :math:`(n,2n)` transfer kernel. It does
**not** enter the production matrix :math:`\mathbf{F}` — the two
neutrons are not produced with the fission spectrum :math:`\chi`, so
production is :math:`\nu\Sigma_f` only. (A retired bespoke formulation
added :math:`2\,\mathrm{colsum}(\Sigma_2)` to the production numerator,
double-counting the :math:`(n,2n)` neutrons; see
:ref:`theory-homogeneous` for the de-vacuum case.)

**Scattering convention.**
:attr:`~orpheus.data.macro_xs.mixture.Mixture.SigS` stores matrices
in ``SigS[g_from, g_to]`` order — **the source uses the transpose**,
:math:`Q_{\rm scatter} = \Sigma_{s}^{\mathsf T}\phi`. The same
transpose appears in the removal matrix
:math:`\Sigma_{s0}^{\mathsf T}` above. This is the single convention
every ORPHEUS solver follows; mis-transposing is the most common
bug when porting from other codes and is caught by L0 spectrum
tests on asymmetric scattering matrices.


Implementation
--------------

The solver is the single function
:func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite`. It
assembles the loss matrix from the **model-shared transport
operators** over a meshless single-cell
:class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` (campaign
#276, Cardinal Rule 2 — the infinite-medium spectrum runs through the
same operator algebra as the meshed SN solver, not a bespoke matrix),
then takes the dominant eigenpair directly:

* **Loss matrix** :math:`\mathbf{A} = C - K_\mathrm{iso} =
  \operatorname{diag}(\Sigma_t) - \Sigma_{s0}^{\mathsf T} -
  2\Sigma_2^{\mathsf T}`, with :math:`C = \operatorname{diag}(\Sigma_t)`
  the collision diagonal and :math:`K_\mathrm{iso}` the sum of the
  model-shared
  :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicScattering`
  (:math:`\Sigma_{s0}^{\mathsf T}`) and
  :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicN2N`
  (:math:`2\Sigma_2^{\mathsf T}`) operators. The composed operator
  :math:`C - K_\mathrm{iso}` is materialised densely via its own
  :meth:`~orpheus.numerics.operator.LinearOperator.as_matrix`
  apply-to-basis (``basis_shape=(ng, 1)``). Streaming :math:`L \equiv 0` in
  an infinite medium and is dropped.
* **Production dyad** :math:`\mathbf{F} = \chi \otimes \nu\Sigma_f`,
  the rank-1 form of
  :class:`~orpheus.transport.operators.fission.FissionOperator`, likewise
  materialised densely via its own
  :meth:`~orpheus.numerics.operator.LinearOperator.as_matrix`.
* **Eigenpair** :math:`k_\infty = \lambda_{\max}(\mathbf{A}^{-1}\mathbf{F})`
  and the dominant right eigenvector, computed by one
  :func:`numpy.linalg.solve` (to apply :math:`\mathbf{A}^{-1}`) plus
  one :func:`numpy.linalg.eig`. There is **no inner or outer
  iteration**.

The function normalises the flux so that the **fission** production
rate :math:`\nu\Sigma_f\cdot\phi` equals :math:`100\ {\rm n/cm^3/s}`
(production is :math:`\nu\Sigma_f` only — :math:`(n,2n)` lives in
:math:`\mathbf{A}`, not :math:`\mathbf{F}`), computes the one-group
collapsed production and absorption cross sections, and packages
everything into a
:class:`~orpheus.homogeneous.solver.HomogeneousResult`.

**Why no iteration is needed.**
In an infinite homogeneous medium there is no spatial eigenmode
spectrum to filter and no spatial coupling to invert iteratively —
the loss operator is a single :math:`G \times G` dense block, so the
fundamental eigenpair is taken in closed form by the dense
eigensolver. This makes the homogeneous solver the one deterministic
solver in ORPHEUS with no iteration at all, and the instantaneous
reference eigenvalue for every other solver on a homogeneous problem.


API Reference
-------------

.. automodule:: orpheus.homogeneous.solver
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:
