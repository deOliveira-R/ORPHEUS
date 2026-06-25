r"""Boundary conditions for transport solvers — trace-law architecture.

This package ships the post-12-wave **trace-law / realizer / tensor-
algebra** architecture for transport-equation boundary conditions
described by Grand Report v3 §16, §16A, and §16A.10. The narrative
theory page is :doc:`/theory/boundary_conditions`; the design plan
is ``.claude/plans/transient-giggling-cake.md`` (the 12-wave
"transient-giggling-cake" effort).


Architecture: three layers, two registries
==========================================

A boundary condition on a transport interface is, in the
method-agnostic sense, a single affine map on the boundary trace

.. math::

    \gamma_- \psi \;=\; R\,G\,\gamma_+ \psi \;+\; q,
    :label: bc-affine-form-init

where :math:`\gamma_\pm` are the **inflow / outflow trace operators**
that restrict :math:`\psi(\mathbf{r}, \Omega)` to the directional
halves of :math:`\partial\Omega`; :math:`G : \Gamma_+ \to \Gamma_+`
is the **geometric map** (a permutation, pushforward, angular
average, or spatial wrap); :math:`R : \Gamma_+ \to \Gamma_-` is the
**response kernel** (a scalar amplitude in :math:`[0, 1]` for
sub-Markov BCs); and :math:`q \in \Gamma_-` is the optional
**prescribed inflow source** (zero for the homogeneous case).

The §16A.3 decomposition splits this map into three concrete layers:

1. **Trace structure** — the unified typed
   :class:`~orpheus.numerics.spaces.trace_space.TraceSpace`
   (every supported mesh — 1-D Cartesian / spherical / cylindrical
   + 2-D Cartesian — post Issue #188; 2-D cylindrical
   :class:`Mesh2D` is the only mesh that still raises, deferred
   until a 2-D cylindrical SN sweep ships). It carries the signed
   :math:`\Omega_n \cdot \hat n_f` per face; inflow and outflow are
   *selectors* over the sign predicate
   :math:`\mathrm{sign}(\Omega_n \cdot \hat n_f)`.
2. **Boundary law** — :class:`BoundaryTraceLaw` ABC + six concrete
   subclasses, each carrying the three first-class properties
   (:attr:`~BoundaryTraceLaw.geometry_map`,
   :attr:`~BoundaryTraceLaw.response_kernel`,
   :attr:`~BoundaryTraceLaw.source`) and the :attr:`creates_sweep_cycle`
   ``ClassVar`` signal. Method-agnostic.
3. **Method realisation** —
   :class:`~orpheus.geometry.boundary.BoundaryRealizer` Protocol +
   one functional implementation
   (:class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer`) plus
   four stubs (``MoCBoundaryRealizer``, ``MCBoundaryRealizer``,
   ``CPBoundaryRealizer``, ``DiffusionBoundaryRealizer``) holding
   the dispatch architecture in place for future modernisation of
   each method.

Two registries connect the layers:

* **Law registry** — :attr:`BoundaryTraceLaw.registry` (a class-level
  dict maintained by
  :class:`~orpheus.numerics.registry.RegistryMixin`). Keyed by the
  ``BC.kind`` string (``"vacuum"``, ``"reflective"``, ``"white"``,
  ``"periodic"``, ``"albedo"``, ``"prescribed_inflow"``). Concrete
  laws self-register at import time via the ``key=`` class-creation
  kwarg.
* **Realizer registry** — :class:`BoundaryRealizerRegistry` (a
  stand-alone class with a class-level ``_registry`` dict). Keyed
  by method name (``"SN"``, ``"MoC"``, ``"MC"``, ``"CP"``,
  ``"diffusion"``). Each realizer self-registers via the
  :meth:`BoundaryRealizerRegistry.register` decorator at module
  import time.

The two are **disjoint by design** (§16A.11). They describe two
orthogonal extension axes: adding a new BC type is one law class +
one registry entry (every existing realizer adds a dispatch branch);
adding a new transport method is one realizer class + one registry
entry (every existing law gets one new realisation branch).


Concrete laws (Wave 7 vocabulary, post-#186 descriptor model)
=============================================================

Each of the six concrete :class:`BoundaryTraceLaw` subclasses is a
**frozen dataclass descriptor** — instantiate it to declare the
law's parameters; realise it (via
:class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer`) to
obtain a 1-arg callable :class:`LinearOperator`. None of the
concrete laws has an :meth:`apply` method; the realiser is the
sole bridge. The canonical SN-realised representation per law:

* :class:`VacuumInflow` (registry key ``"vacuum"``, deprecated alias
  ``VacuumBoundaryOperator``) — the rank-0 case :math:`R = G = q
  = 0`. SN realises to
  :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`
  with the per-face inflow indices; this zeroes **only the
  inflow ordinates** and preserves the outflow trace (the §16A.5
  trace-correct representation). The pre-refactor zeros-all
  body that lived on the Option-A standalone ``apply`` path was
  **deleted in Issue #186 / B3 + β2** — there is no longer an
  alternate path; the inflow-only mask is the unique vacuum
  semantics. :attr:`creates_sweep_cycle` = ``False``.
* :class:`ReflectiveBoundary(axis, albedo)` (registry key
  ``"reflective"``, deprecated alias ``SpecularBoundaryOperator``)
  — :math:`R = G_{\text{refl}} \cdot \alpha` with the geometric
  operator the ordinate permutation
  ``quadrature.reflection_index(axis)``. SN realises to
  :class:`~orpheus.numerics.operator.PermutationOperator` (α=1
  fast path) or
  ``ScaledOperator(α, PermutationOperator)`` (α ≠ 1).
  :attr:`creates_sweep_cycle` = **``True``** — signals §15A.2
  sweep-cycle detection.
* :class:`WhiteBoundary(axis, outward_sign, albedo)` (registry key
  ``"white"``, deprecated alias ``WhiteBoundaryOperator``) —
  :math:`R = G_{\text{diff}} \cdot \alpha` with the geometric
  operator the cosine-weighted Lambertian average over the outgoing
  hemisphere. SN realises to
  :class:`~orpheus.sn.angular_operator.AngularAverageOperator` (α=1
  fast path) or scaled.
  :attr:`creates_sweep_cycle` = ``False``.
* :class:`PeriodicBoundary` (registry key ``"periodic"``, deprecated
  alias ``PeriodicBoundaryOperator``) — :math:`R` is a spatial
  pushforward (wrap-around to the opposite face) with
  :math:`\alpha = 1`. SN realises to
  :class:`~orpheus.numerics.operator.PeriodicWrapOperator`
  (currently an angular identity; the spatial-pushforward extension
  is tracked as a follow-up under ``module:sn``).
  :attr:`creates_sweep_cycle` = **``True``**.
* :class:`AlbedoBoundary(albedo)` (registry key ``"albedo"``,
  deprecated alias ``AlbedoBoundaryOperator``) — :math:`R = I \cdot
  \alpha` with the geometric operator the angular identity. SN
  realises to :class:`~orpheus.numerics.operator.ZeroOperator` (α=0),
  :class:`~orpheus.numerics.operator.IdentityOperator` (α=1), or
  ``ScaledOperator(α, IdentityOperator)`` (α ∉ {0, 1}).
  :attr:`creates_sweep_cycle` = ``False``.
* :class:`PrescribedInflow(source)` (registry key
  ``"prescribed_inflow"``, Wave 7 addition) — the rank-0 affine BC
  :math:`R = G = 0`, :math:`q \neq 0`. SN realises to
  :class:`~orpheus.sn.angular_operator.IncomingSourceOperator(source)`
  whose :meth:`apply` ignores the outgoing flux and returns
  ``source.evaluate(psi_out.shape)``. Consumes a
  :class:`InflowSourceSpec` (typically :class:`NoSource` for the
  homogeneous case or :class:`ConstantInflowSource(value)` for a
  uniform inflow).
  :attr:`creates_sweep_cycle` = ``False``.


Rank-N boundaries via the descriptor-tree algebra
=================================================

Rank-N (Marshak, partial-current) boundary conditions are
**not** a dedicated class. They are built directly through the
**descriptor-tree algebra** that
:class:`BoundaryTraceLaw` exposes via its ``+`` / ``-`` / ``*`` /
``/`` / unary ``-`` dunders. These return
:class:`~orpheus.geometry.boundary._composition.LawSum` /
:class:`~orpheus.geometry.boundary._composition.LawScaled` nodes —
pure descriptor structures with **no** ``apply`` method. The
:func:`~orpheus.sn.boundary_realize.realize_recursively` walker is
the **sole** type transformer from descriptor tree to operator
tree (Issue #186 / B3 + β2, 2026-05-11):

.. code-block:: python

    from orpheus.geometry.boundary import (
        ReflectiveBoundary, WhiteBoundary,
    )
    from orpheus.sn.boundary_realize import realize_recursively
    from orpheus.sn.boundary_realizer import SNMethodSpace

    # Build the descriptor tree (no realisation yet).
    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    white = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
    marshak_law = 0.3 * spec + 0.7 * white
    # marshak_law is LawSum(LawScaled(0.3, ReflectiveBoundary(...)),
    #                       LawScaled(0.7, WhiteBoundary(...))).
    # Not callable: marshak_law.apply does NOT exist.

    # Realise the tree at one face.
    ms = SNMethodSpace.for_face(...)
    marshak_op = realize_recursively(marshak_law, ms)
    # marshak_op is:
    #   OperatorSum(ScaledOperator(0.3, PermutationOperator(...)),
    #               ScaledOperator(0.7, AngularAverageOperator(...)))
    psi_in = marshak_op.apply(psi_out)   # 1-arg LinearOperator

The output is a Wave-0
:class:`~orpheus.numerics.operator.OperatorSum` of
:class:`~orpheus.numerics.operator.ScaledOperator` wrappers around
realised Wave-0 primitives, consumable by the SN sweep / Krylov
path via the uniform 1-arg :meth:`apply`.

The descriptor-tree and operator-tree algebras are **two separate
type families**: the descriptor tree
(:class:`~orpheus.geometry.boundary.LawSum` /
:class:`~orpheus.geometry.boundary.LawScaled` over
:class:`BoundaryTraceLaw` leaves) carries no ``apply``; the
operator tree (:class:`OperatorSum` / :class:`ScaledOperator` over
:class:`LinearOperator` leaves) carries :meth:`apply`. The two
algebras never inter-compose — mixing a :class:`LawNode` with an
already-realised :class:`LinearOperator` via ``+`` is a static-type
error; the caller MUST :func:`realize_recursively` the descriptor
tree first. The static type checker enforces this distinction —
no runtime contract is needed.

Two retired predecessors converged on this design:

* **Wave 11**'s ``MixedBoundaryOperator(components: list[tuple[
  float, BoundaryOperator]])`` class delayed realisation until
  apply-time; deleted because the delayed-realisation pattern
  broke down once vacuum needed per-face inflow indices that the
  bare-law container had no access to.
* **β1** (the Issue #186 / B3 interim that kept
  :class:`LinearOperator` on :class:`BoundaryTraceLaw`):
  ``0.3 * spec + 0.7 * white`` produced an :class:`OperatorSum`
  with raw-law leaves whose realisation was deferred to
  apply-time. β1 was algebraically equivalent to β2 but conflated
  the two type families — the resulting :class:`OperatorSum`
  instance looked like an operator tree to the type checker
  while behaving as a descriptor tree at runtime.

See Bell & Glasstone 1970 §1.5 for the physics; see
:doc:`/theory/boundary_conditions` § "Descriptor-tree algebra for
rank-N boundaries" and § "The trace-law descriptor model
(Issue #186 / B3 + β2)" for the architectural retrospective.


Descriptor-tree composition module
==================================

The :mod:`_composition` submodule houses the descriptor-tree
composition dataclasses introduced in Issue #186:

* :class:`LawSum(a, b)` — frozen dataclass representing
  :math:`a + b` over two :data:`LawNode` operands. Closed under
  ``+`` / ``-`` / ``*`` / ``/`` / unary ``-``.
* :class:`LawScaled(scalar, inner)` — frozen dataclass
  representing :math:`\alpha \cdot \mathrm{inner}`. Closed under
  the same dunders; ``LawScaled.__mul__`` performs **constant
  folding** (``α * (β * x) → (α β) * x`` at construction).
* :data:`LawNode` — type alias
  ``Union[BoundaryTraceLaw, LawSum, LawScaled]``, the input type
  of :func:`realize_recursively`.

All three are re-exported below for downstream consumers. None
has an :meth:`apply` method — the descriptor tree becomes
callable only after :func:`realize_recursively`. See
:mod:`._composition` for the algebra closure details.


Universal invariants + named-error catalog
==========================================

Each :class:`BoundaryTraceLaw` subclass overrides the five universal
``assert_*`` invariants on the ABC (per Grand Report v3 §16A.12 +
§27.6) where applicable. Eight typed errors in :mod:`._errors`
replace the pre-refactor generic :class:`ValueError` raises:

* :class:`IncomingOutgoingTraceClassificationError` — ERR-040, fires
  on tangential ordinates where strict partition was required.
* :class:`VacuumAppliedToOutgoingTraceError` — ERR-041, fires when
  vacuum is mistakenly applied on :math:`\Gamma_+`.
* :class:`BoundaryGeometryMapNotMeasurePreservingError` — ERR-042,
  fires when the geometric map :math:`G` does not preserve
  :math:`w(\Omega)|\Omega \cdot \hat n|`.
* :class:`BoundaryResponseNotPositiveError` — ERR-043, fires on
  negative response kernel output.
* :class:`ReflectionNotInvolutiveError` — ERR-044, fires when
  ``perm[perm] != arange``.
* :class:`ReflectionDidNotMapInflowToOutflowError` — ERR-045, fires
  when a reflection maps an inflow ordinate to itself.
* :class:`SubmarkovViolationError` — ERR-046, fires when
  :math:`\alpha > 1` on a sub-Markov BC (albedo, white).
* :class:`BoundarySourceNotOnIncomingTraceError` — ERR-047, fires on
  a :class:`InflowSourceSpec` with nonzero entries on :math:`\Gamma_+`.

All eight extend :class:`ValueError` (via the
:class:`BoundaryError` base) for backward compatibility with
existing ``except ValueError`` consumers. The
``@pytest.mark.catches("ERR-NNN")`` decorators on the relevant
fault-injection tests in
:mod:`tests.geometry.test_bc_universal_invariants` pin the error
firing under the right conditions.


Cross-method realizer stubs
===========================

The :class:`BoundaryRealizer` Protocol has one functional
implementation today —
:class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer`. Four
stub realizers in :mod:`orpheus.moc.boundary_realizer`,
:mod:`orpheus.mc.boundary_realizer`,
:mod:`orpheus.cp.boundary_realizer`, and
:mod:`orpheus.diffusion.boundary_realizer` self-register at import
time but raise :class:`NotImplementedError` from :meth:`realize`.
The stubs hold the dispatch architecture in place so a future
modernisation session for each method can ship a functional body
without touching the SN side. Per the "Unify after two instances"
rule, a shared ``BoundaryRealizerBase`` ABC is deferred until the
second functional realizer ships — see
:doc:`/theory/boundary_conditions` § "Cross-method realizer stubs".


Package layout (Wave 4 source-layout split, post-#186 descriptor cleanup)
=========================================================================

* :mod:`_base` -- :class:`BoundaryTraceLaw` ABC (Wave 3 / Wave 7
  ABC merge; Issue #186 / B3 + β2 dropped
  :class:`~orpheus.numerics.operator.LinearOperator`
  inheritance and the abstract ``apply``).
* :mod:`_composition` -- :class:`LawSum` /
  :class:`LawScaled` / :data:`LawNode` descriptor-tree composition
  dataclasses (Issue #186 / B3 + β2).
* :mod:`_errors` -- :class:`BoundaryError` + 8 typed subclasses
  (Wave 3 / ERR-040..ERR-047).
* :mod:`_source` -- :class:`InflowSourceSpec` Protocol +
  :class:`NoSource` + :class:`ConstantInflowSource` (Wave 3).
* :mod:`_realizer` -- :class:`BoundaryRealizer` Protocol +
  :class:`BoundaryRealizerRegistry` (Wave 5).
* :mod:`_bound_compat` -- ``_BoundBoundaryOperator`` shim wrapping
  the realised operator with a ``kind`` string tag (Wave 8
  introduced; Wave 9 added curvilinear bound-quadrature path,
  retired Issue #176 / C176.1; Issue #186 / C-B3.4 dropped the
  ``*_extra, **_kw`` swallow). Now a **strict 1-arg** passthrough.
  **Internal** — not in this module's :attr:`__all__`.
* :mod:`vacuum` -- :class:`VacuumInflow`.
* :mod:`reflective` -- :class:`ReflectiveBoundary`.
* :mod:`white` -- :class:`WhiteBoundary`.
* :mod:`periodic` -- :class:`PeriodicBoundary`.
* :mod:`albedo` -- :class:`AlbedoBoundary`.
* :mod:`prescribed_inflow` -- :class:`PrescribedInflow` (Wave 7).

The pre-Wave-7 ``mixed.py`` submodule (carrying the now-retired
``MixedBoundaryOperator``) was **deleted in Wave 11**; the registry
no longer contains a ``"mixed"`` key and a test pins that absence
(:func:`tests.geometry.test_boundary.test_registry_contains_all_primitives`).

The Wave-7 deprecated aliases (``VacuumBoundaryOperator`` =
:class:`VacuumInflow`, etc.) are re-exported below for backward
compatibility with consumers that import the old names. They will
be removed in a future cleanup wave once every consumer migrates.


Cross-references
================

* :doc:`/theory/boundary_conditions` — the master theory page
  carrying the full §16A.3 three-layer decomposition, the affine
  form derivation, the dual-registry design, the universal
  invariants, the named-error catalog, the vacuum semantic
  correction (§16A.5), the Wave-0 rank-N algebra, the
  worked end-to-end example, and the Cartesian / curvilinear
  split.
* :doc:`/theory/discrete_ordinates` § "Boundary Conditions" —
  the SN-side consumption: how :class:`SNMesh._resolve_bcs` reads
  the law registry, builds the method space, and dispatches
  through :class:`SNBoundaryRealizer` to produce the resolved
  1-arg :class:`~orpheus.numerics.operator.LinearOperator`.
* :doc:`/theory/operator_algebra` § "Boundary conditions as
  Wave-0 / Wave-1 primitives" — the operator-algebra view of the
  realisation: every realised BC IS a Wave-0
  :class:`LinearOperator` composable with the rest of the algebra.
* :mod:`orpheus.numerics.operator` — Wave-0 primitives
  (:class:`PermutationOperator`, :class:`IncomingOrdinateMaskTensor`,
  :class:`PeriodicWrapOperator`, plus the
  :class:`OperatorSum` / :class:`ScaledOperator` composers that
  appear in the *operator* tree).
* :mod:`._composition` — :class:`LawSum` / :class:`LawScaled` /
  :data:`LawNode` composers that appear in the *descriptor* tree
  (Issue #186 / B3 + β2).
* :mod:`._base` — :class:`BoundaryTraceLaw` ABC: pure descriptor
  with no ``apply``; carries the minimal algebra dunders that
  build descriptor-tree nodes.
* :mod:`orpheus.sn.boundary_realizer` —
  :class:`SNBoundaryRealizer` (functional realizer dispatching by
  ``isinstance``; the leaf descriptor → operator transformer).
* :mod:`orpheus.sn.boundary_realize` —
  :func:`realize_recursively` walker, the **type transformer**
  from a descriptor tree
  (``BoundaryTraceLaw | LawSum | LawScaled``) to an operator tree
  (Wave-0 :class:`OperatorSum` / :class:`ScaledOperator` around
  realised 1-arg :class:`LinearOperator` leaves).


References
----------

* Grand Report v3 §15.2, §16, §16A.1–5, §16A.10–12, §26A.4, §27.6,
  §28. Source: ``.claude/plans/neutron_transport_grand_report_v3.md``.
* The 12-wave refactor plan:
  ``.claude/plans/transient-giggling-cake.md`` (Waves 0-12).
* The post-Wave-12 cleanup plans:

  - ``.claude/plans/curvilinear-realizer-and-2arg-cleanup.md`` —
    Issue #188 + #176 (curvilinear realiser unification + Option-A
    interim).
  - ``.claude/plans/bc-trace-law-descriptor-cleanup.md`` —
    Issue #186 / B3 + β2 (descriptor-model cleanup).
* Lewis, E. E. & Miller, W. F. (1984). *Computational Methods of
  Neutron Transport*. American Nuclear Society. §3.4 (boundary
  conditions in transport).
* Bell, G. I. & Glasstone, S. (1970). *Nuclear Reactor Theory*.
  Van Nostrand Reinhold. §1.5 (albedo, white, and Marshak boundary
  conditions).
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Abstract base — Wave 7 merged the legacy ``BoundaryOperator`` ABC into
# :class:`BoundaryTraceLaw`. Wave O step O.4a.1 retired the deprecated
# ``BoundaryOperator`` alias (all consumers moved to ``BoundaryTraceLaw``;
# the name is now free for the ``BlockRole.BOUNDARY`` marker in
# :mod:`orpheus.numerics.operator`).
# ---------------------------------------------------------------------------

from ._base import BoundaryTraceLaw

# ---------------------------------------------------------------------------
# Issue #186 (B3 + β2) -- descriptor-tree composition. LawSum / LawScaled
# form a closed algebra over BoundaryTraceLaw | LawSum | LawScaled, used
# for rank-N boundary composition (e.g. ``0.3 * spec + 0.7 * white``).
# Realised to a Wave-0 operator tree via
# :func:`orpheus.sn.boundary_realize.realize_recursively`.
# ---------------------------------------------------------------------------

from ._composition import LawNode, LawScaled, LawSum

# ---------------------------------------------------------------------------
# Wave 3 -- typed error catalog (ERR-040..ERR-047).
# ---------------------------------------------------------------------------

from ._errors import (
    BoundaryError,
    BoundaryGeometryMapNotMeasurePreservingError,
    BoundaryResponseNotPositiveError,
    BoundarySourceNotOnIncomingTraceError,
    IncomingOutgoingTraceClassificationError,
    ReflectionDidNotMapInflowToOutflowError,
    ReflectionNotInvolutiveError,
    SubmarkovViolationError,
    VacuumAppliedToOutgoingTraceError,
)

# ---------------------------------------------------------------------------
# Wave 3 -- prescribed-inflow source.
# ---------------------------------------------------------------------------

from ._source import InflowSourceSpec, ConstantInflowSource, NoSource

# ---------------------------------------------------------------------------
# Wave 5 -- BoundaryRealizer Protocol + per-method registry.
# ---------------------------------------------------------------------------

from ._realizer import (
    BoundaryRealizer,
    BoundaryRealizerRegistry,
    BoundaryRealizerRegistryError,
)

# ---------------------------------------------------------------------------
# Concrete BCs -- split into per-BC submodules in Wave 4, renamed to
# Grand Report v3 vocabulary in Wave 7. Wave O step O.4a.1 retired the
# deprecated ``*BoundaryOperator`` aliases (all consumers moved to the
# canonical names below).
# ---------------------------------------------------------------------------

from .albedo import AlbedoBoundary
from .periodic import PeriodicBoundary
from .prescribed_inflow import PrescribedInflow
from .reflective import ReflectiveBoundary
from .vacuum import VacuumInflow
from .white import WhiteBoundary


__all__ = [
    # Abstract base
    "BoundaryTraceLaw",
    # Descriptor-tree composition (Issue #186)
    "LawNode",
    "LawScaled",
    "LawSum",
    # Errors
    "BoundaryError",
    "BoundaryGeometryMapNotMeasurePreservingError",
    "BoundaryResponseNotPositiveError",
    "BoundarySourceNotOnIncomingTraceError",
    "IncomingOutgoingTraceClassificationError",
    "ReflectionDidNotMapInflowToOutflowError",
    "ReflectionNotInvolutiveError",
    "SubmarkovViolationError",
    "VacuumAppliedToOutgoingTraceError",
    # Source
    "InflowSourceSpec",
    "ConstantInflowSource",
    "NoSource",
    # Wave 5 -- realizer Protocol + registry
    "BoundaryRealizer",
    "BoundaryRealizerRegistry",
    "BoundaryRealizerRegistryError",
    # Concrete BCs (Wave 7 canonical names)
    "AlbedoBoundary",
    "PeriodicBoundary",
    "PrescribedInflow",
    "ReflectiveBoundary",
    "VacuumInflow",
    "WhiteBoundary",
]
