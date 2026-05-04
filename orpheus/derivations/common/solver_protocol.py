r"""The :class:`TransportSolver` Protocol — the unifying contract every
math-heart class and discrete-mesh adapter conforms to.

This module is the **structural unification** of the math-heart pattern
realised across the project:

* :class:`~orpheus.derivations.continuous.trajectory_resolvent.Billiard`
  — the bouncing-trajectory operator on phase space.
* :class:`~orpheus.derivations.continuous.fn_method.MomentSpace` — the
  Galerkin half-range projection on the F_N moment basis.
* The forthcoming ``Spectrum`` (Case singular eigenfunction expansion in
  :mod:`~orpheus.derivations.continuous.singular_eigenfunction`) and
  ``BasisSpace`` (Dahl–Sjöstrand Legendre-Galerkin block-matrix
  linearisation in
  :mod:`~orpheus.derivations.continuous.galerkin_spectral`).
* The production discrete-mesh solvers
  :class:`~orpheus.cp.solver.CPSolver` and
  :class:`~orpheus.sn.solver.SNSolver`, which conform natively
  via their ``from_problem`` factories (Step 4 of the input-cleanup
  track, 2026-05-04). The earlier test-only adapter scaffold
  (``CPMeshAdapter`` / ``SNMeshAdapter``) was retired when the
  production classes took on the Protocol surface directly.

Per the project's *unify-after-two-instances* discipline (see
``feedback_unify_after_two_instances.md``), the Protocol was *gated* on
two working math-heart instances: ``MomentSpace`` and ``Billiard``. With
both shipped, the variation patterns are observable and the Protocol
codifies what they share.

Three axes of carving
---------------------

The project's cross-method architecture carves the transport problem on
three orthogonal axes (see :doc:`/theory/transport_solver_protocol`):

1. **What** — the materials (``dict[int, Mixture]``) and the domain
   (``GeometrySpec``). These are *method-agnostic* — the same problem
   description served to a CP solver, an SN solver, and an F_N moment
   space. The Protocol exposes both via :attr:`materials` and
   :attr:`geometry_spec`.

2. **How** — the method-specific computational specialisation. F_N
   commits to a moment basis order :math:`N`; trajectory_resolvent
   commits to a billiard table classified by orbit-space rank;
   discrete CP/SN commit to a discretisation. Each Protocol-conforming
   class owns its own method kwargs (``fn_order`` for ``MomentSpace``,
   ``quadrature`` for ``Billiard``, ``DiscretizationSpec`` for
   :class:`~orpheus.cp.solver.CPSolver` /
   :class:`~orpheus.sn.solver.SNSolver`); the Protocol does NOT
   dictate them — that would couple unrelated methods.

3. **What is asked** — :meth:`solve_critical` returns a
   :class:`~orpheus.derivations.common.solution_types.CriticalSolution`;
   future :meth:`solve_fixed_source` returns a
   :class:`~orpheus.derivations.common.solution_types.FluxSolution`.
   These shared result types are the cross-method comparator's
   substrate; the Protocol's behavioural contract is "anything that
   returns a ``CriticalSolution`` from ``solve_critical(materials,
   geometry_spec)`` is a transport solver, regardless of pillar."

Why a Protocol, not an ABC
--------------------------

A :class:`typing.Protocol` (PEP-544) is **structural typing** — any
class that has the listed attributes / methods automatically conforms,
without inheriting from a base. This matches our reality:

* ``Billiard`` and ``MomentSpace`` were designed independently before
  the Protocol existed; forcing them under a common ABC would have
  required retroactive subclassing.
* Production CP / SN classes (``CPSolver``, ``SNSolver``) take on
  Protocol conformance through their ``from_problem`` factories
  without inheriting from any base — the legacy
  ``__init__`` + function-level ``solve_cp`` / ``solve_sn`` API
  remains intact for callers that don't need the cross-method
  surface.
* The Protocol is :func:`runtime_checkable`, so ``isinstance(x,
  TransportSolver)`` works for the schema-gate tests in
  :mod:`tests.derivations.test_transport_solver_protocol` without a
  registration step.

Foundation-tier, not L1
-----------------------

Conformance to the Protocol is a software invariant — it does NOT
verify a physics claim. The L1 chains stay where they were
(``test_fn_la13511_*.py``, ``test_peierls_greens_function_*.py``,
``tests.cross_method.test_eigenvalue``); the Protocol unifies the
*shape* those L1-verified components must expose. Tests in
:mod:`tests.derivations.test_transport_solver_protocol` are
foundation-tagged accordingly.

References
----------

* :doc:`/theory/transport_solver_protocol` — the rich-narrative theory
  page (three-axis carve, ABNF surface, cross-method positioning).
* :doc:`/skills/algebra-of-record` § "Branch 1 / Branch 2 separation"
  — the Protocol lives at the math layer, not inside any algebra-of-
  record. Bit-equality preservation through the Protocol facade is
  what guarantees Branch-1 ⊥ Branch-2 cross-checks remain valid.
* :doc:`/skills/vv-principles` § "Hierarchical claim taxonomy" — the
  software-invariant tier this Protocol's tests live at.
"""
from __future__ import annotations

from typing import Protocol, runtime_checkable

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.geometry_spec import GeometrySpec
from orpheus.derivations.common.solution_types import CriticalSolution


@runtime_checkable
class TransportSolver(Protocol):
    r"""Structural contract for a one-speed-or-multi-group transport solver.

    A class conforms to :class:`TransportSolver` iff it exposes:

    * :attr:`materials` — the production-protocol cross-section payload
      (``dict[int, Mixture]``), keyed by material ID. Reflects the
      *what* axis of the three-axis carve.
    * :attr:`geometry_spec` — the method-agnostic
      :class:`~orpheus.derivations.common.geometry_spec.GeometrySpec`
      describing the spatial domain + boundary conditions. Reflects
      the *what* axis.
    * :attr:`method_name` — a stable string tag identifying the
      computational pillar (``"trajectory_resolvent"``, ``"fn_method"``,
      ``"singular_eigenfunction"``, ``"galerkin_spectral"``, ``"cp"``,
      ``"sn"``, ...). Used by cross-method comparators to group
      results without case-by-case isinstance checks.
    * :meth:`solve_critical` — returns a
      :class:`CriticalSolution`. The eigenvalue
      semantics (``k_eff`` vs ``k_inf``) is recorded in the result's
      ``eigenvalue_kind`` field — the Protocol does NOT prescribe
      which.

    Non-Protocol surface — kept on each conformant class
    -----------------------------------------------------

    Method kwargs (``fn_order``, ``n_quad``, ``alpha``, ``DiscretizationSpec``,
    flux-reconstruction strategy, …) are NOT part of the Protocol.
    Each class owns its own construction surface; the Protocol is
    deliberately small. The discipline "match the Protocol to what
    Billiard + MomentSpace already expose; if a Protocol field has no
    instance using it, drop it" was applied during design.

    Production discrete-mesh solvers
    --------------------------------

    The Protocol covers BOTH continuous reference solvers (where the
    L1 truth lives) AND production discrete solvers — the latter
    via the ``CPSolver.from_problem`` / ``SNSolver.from_problem``
    factories that mint the production classes onto the Protocol
    natively (Step 4 of the input-cleanup track, 2026-05-04). The
    discrete side does NOT replace continuous-reference verification
    — it enables L4 *agreement* gates between production and
    reference after each side has its own L1 chain. Per
    :doc:`/skills/algebra-of-record` § "Structural independence
    applies above the trusted-library line", the production CP / SN
    solvers consume their own scipy / numpy primitives via separate
    code paths and share NO in-house primitive with the
    continuous-reference family above the trusted-library line.

    Examples
    --------

    Continuous-reference math-heart classes::

        from orpheus.derivations.continuous.trajectory_resolvent import (
            Billiard,
        )
        from orpheus.derivations.continuous.fn_method import MomentSpace

        b: TransportSolver = Billiard.from_problem(
            materials={0: pu_mixture},
            geometry_spec=GeometrySpec(geometry="sphere", ...),
            alpha=0.5,
        )
        m: TransportSolver = MomentSpace.from_problem(
            materials={0: pu_mixture},
            geometry=GeometrySpec(geometry="slab", ...),
            fn_order=10,
        )
        assert isinstance(b, TransportSolver)
        assert isinstance(m, TransportSolver)
        assert b.method_name == "trajectory_resolvent"
        assert m.method_name == "fn_method"

    Notes on the iter_pillars tag
    ------------------------------

    ``method_name`` is a *string*, not an enum, on purpose: new
    pillars (e.g. spectral_resolvent, mc_reference) can register
    themselves without modifying this module. Cross-method
    comparators that need to group by pillar should use string
    equality and accept "unknown" for tags they don't recognize.

    Attributes
    ----------
    materials : dict[int, Mixture]
        Production-protocol cross sections, keyed by material ID. The
        same payload shape consumed by ``solve_cp`` / ``solve_sn`` /
        ``solve_moc``. Multi-region cases use multiple keys; bare-
        critical / single-zone cases typically have one key
        (``mat_id == 0``).
    geometry_spec : GeometrySpec
        The domain + boundary conditions. Method-agnostic; both the
        continuous-reference and the discrete-mesh paths consume it.
    method_name : str
        Stable tag identifying the computational pillar.
    """

    materials: dict[int, Mixture]
    geometry_spec: GeometrySpec
    method_name: str

    def solve_critical(self) -> CriticalSolution:
        r"""Solve the critical configuration on this transport problem.

        Returns the
        :class:`~orpheus.derivations.common.solution_types.CriticalSolution`
        carrying the eigenvalue (``k_eff`` or ``k_inf``), the
        configuration parameter (critical half-thickness, critical
        radius, fixed-geometry configuration scalar), and the
        method's diagnostic metadata.

        The eigenvalue / parameter conventions (``parameter_kind``,
        ``eigenvalue_kind``) are NOT prescribed by the Protocol —
        each method may report its natural quantity. Cross-method
        comparators read these fields off the result before
        comparing.

        Note that :meth:`solve_critical` MAY accept method-specific
        keyword overrides on each conforming class (e.g. ``fn_order``
        override on ``MomentSpace``, quadrature overrides on
        ``Billiard``). The Protocol's contract is only that the
        *bare* call ``self.solve_critical()`` returns a
        :class:`CriticalSolution` — kwargs are off-Protocol but
        permitted.

        Returns
        -------
        CriticalSolution
            See
            :class:`~orpheus.derivations.common.solution_types.CriticalSolution`.
        """
        ...


# A registry of stable method-name tags that ship today. New pillars
# extend this list when they land; the registry exists for the
# foundation-tagged "registry is complete" gate that catches drop-outs
# (e.g. someone removing ``method_name`` from a math-heart class).
KNOWN_TRANSPORT_SOLVERS: tuple[str, ...] = (
    "trajectory_resolvent",
    "fn_method",
    "singular_eigenfunction",
    "galerkin_spectral",
    "cp",
    "sn",
)
"""Stable tag set covering today's TransportSolver-conformant classes.

Used by ``test_solver_registry_is_complete`` to detect silent
removal of ``method_name`` from a math-heart or discrete-mesh
adapter. New pillars (singular_eigenfunction, galerkin_spectral,
moc, mc) join this tuple as they land. The string tag is the same
one each class returns from :attr:`TransportSolver.method_name`.
"""


__all__ = [
    "TransportSolver",
    "KNOWN_TRANSPORT_SOLVERS",
]
