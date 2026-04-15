1-D Two-Group Diffusion Solver
================================

Reference for the :mod:`orpheus.diffusion` package — a 1-D two-group
finite-difference diffusion solver with a fuel/reflector core geometry.
The underlying theory (two-group eigenvalue problem, power iteration,
boundary treatment) is covered in the diffusion theory chapter and
verified against analytical reference solutions built by
:mod:`orpheus.derivations.diffusion`.

See :ref:`theory-verification` for the L0/L1 verification cases that
exercise this solver, and :func:`orpheus.derivations.diffusion.solver_cases`
for the set of :class:`~orpheus.derivations._types.VerificationCase`
instances a test consumes via ``ref(name)``.

Solver
------

.. automodule:: orpheus.diffusion.solver
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:
