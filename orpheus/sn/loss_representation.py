r"""Selectable representations of the S\ :sub:`N` loss operator :math:`(L+C)`.

The within-group transport solve :math:`\psi = (L+C)^{-1} q` (the *sweep*)
and its operator twin :math:`(L+C)\,\psi` (the *matvec*) admit several
distinct *algorithms*, each natural for a different mesh:

* a **1-D parallel-prefix scan** (Blelloch 1990 §1.5) — the geometry-blind
  chain recurrence (slab + sphere + cylinder), :func:`._sweep_1d_unified`;
* a **multi-D wavefront walk** over the per-octant anti-hyperplane DAG
  (:meth:`SweepDependencyGraph.apply`), in two buffer policies — a
  full-field buffer (the slow, readable verification oracle) and a rolling
  :math:`(d{-}1)`-frontier window (the fast production path).

Historically the choice between them was a *scattered, procedural* branch
spelled three different ways: ``transport_sweep`` branched on
``sn_mesh.reduced is not None``; the matvec branched in five operator gates
on ``not sn_mesh.is_1d``; and the full-field oracle was reachable only
through hand-built test adapters.  Adding a method, a dimensionality, or a
frontend meant touching all three — an enum-style branch repeated at every
call site (cyclomatic complexity, not abstraction).

This module replaces that with a first-class :class:`LossRepresentation`: each
algorithm is an object that carries **both** the forward ``sweep`` and (from
Phase S2) the ``loss_action`` matvec twin, plus a **declared, queryable**
:meth:`~LossRepresentation.supports` predicate.  The operator and the
``transport_sweep`` dispatcher select one via :func:`default_for` and then
call it branchlessly.

The hierarchy
=============

.. code-block:: text

    LossRepresentation (Protocol: sweep, loss_action, supports)
    ├── _DAGWavefront            ── Cartesian anti-hyperplane DAG family
    │   ├── FullFieldWavefront     buffer = full field     · the ORACLE
    │   └── MovingFrontierWindow   buffer = rolling frontier · production opt
    └── CumprodScan             ── 1-D chain prefix scan, any geometry

``FullFieldWavefront`` and ``MovingFrontierWindow`` consume the **same**
per-octant DAG (the family-owned ``sweep_graphs`` accessor, cached per
mesh shape since S6.4(c)) — they are two *buffer policies* over one
anti-hyperplane walk, already pinned bit-identical by the C3.2b
``window ≡ full`` oracle.  ``CumprodScan`` builds no DAG: a 1-D chain is a
total order, the Blelloch closed form needs no graph.

The governing principle
========================

    *Construct each strategy as general as its algorithm naturally allows.
    Select narrow.  Specialize the implementation only on a measured
    internal performance cost.*

Three separable layers, never conflated:

* **Construct general (capability).**  ``CumprodScan`` is intrinsically 1-D
  (a prefix scan needs a total order → a chain → 1-D; there is no "2-D chain
  scan") — legitimately d-specific *by the algorithm's nature*.
  ``FullFieldWavefront`` and ``MovingFrontierWindow`` are naturally
  d-general (the moving frontier is the :math:`(d{-}1)`-dim rolling slab: a
  point at d=1, a line at d=2, a surface at d=3).
* **Select narrow (policy).**  Whether we *offer / recommend / default* a
  strategy at a given ``(geometry, ndim)`` lives in
  :meth:`~LossRepresentation.supports` / :func:`default_for`, independent of the
  code's capability.  "Don't pick the window at d=1, pick the scan" is a
  recommendation, *not* a reason to leave the window unable to express d=1.
* **Specialize on measured cost only.**  The sole justification to restrict
  an implementation's d-range is a *measured* hot-path regression.

Selection is a single source of truth
======================================

:meth:`~LossRepresentation.supports` returns :class:`Compatibility` — an
``(ok, reason)`` pair.  The same predicate serves three consumers:

#. a (future) teaching frontend — ``[S for S in LOSS_REPRESENTATIONS if
   S.supports(mesh).ok]`` grays-out an inapplicable method *and explains
   why* (pedagogically load-bearing — ORPHEUS teaches reactor physics);
#. the factory :func:`default_for` — picks the best *available* production
   optimization, falling back to the full-field spine when no optimization
   exists, so it is never stuck;
#. the construction guard — :meth:`_LossRepresentation.__post_init__` raises
   :class:`IncompatibleRepresentation` on an illegal pairing, so even a bypassed
   UI cannot build one.

The compatibility signal is the genuine criterion — the coordinate system
(:attr:`SNMesh.is_cartesian`) and the dimensionality (:attr:`SNMesh.ndim`)
— NOT the ``sweep_graphs is None`` substrate proxy.

Phase status (the sweep-strategy carve, plan ``sn_sweep_strategy.md``)
====================================================================

* **S1 (this commit): SWEEP side.**  The protocol, the hierarchy, the three
  strategies as *thin wrappers* over the existing sweep functions, the
  selection layer, and the ``transport_sweep`` rewire.  Bit-identical: each
  strategy wraps exactly the path the old branch chose, so ``default_for``
  reproduces the legacy dispatch byte-for-byte.
* **S2: MATVEC side.**  Each strategy gains ``loss_action`` (the
  :math:`(L+C)\psi` twin); the five operator gates collapse to
  ``strategy.loss_action(...)``.
* **S3: generalize the oracle.**  ``FullFieldWavefront`` is the genuine
  d-generic spine (sweep ``_sweep_full_field`` + matvec
  ``FullFieldWavefront.loss_action`` — since S6.3 both walk
  :meth:`SweepDependencyGraph.residual`; ``supports`` is any-d Cartesian,
  wires the d=1 cumprod-vs-spine equivalence).
* **S4: generalize the window** to ``frontier_dim = d-1``.
* **S6.2 (this commit): RENAME.**  ``SweepStrategy → LossRepresentation``,
  ``residual → loss_action`` — the abstraction is the :math:`(L+C)` loss
  operator's *representation* (a matrix-free traversal; the
  cross-domain-attacker named it "representation" across four expert
  frames).  Pure rename, bit-identical (the bodies still delegate to the
  operator's ``_apply_*``).  S6.3 moves the walk *into* ``loss_action``
  (returns :math:`(L+C)\psi`; the operator's ``−C`` becomes the only glue).

See also
========

* ``.claude/plans/sn_sweep_strategy.md`` — the authoritative design (the
  locked decisions, the verification strategy, phases S0–S6).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Protocol, runtime_checkable

import numpy as np

# This module wraps the sweep bodies, so it imports ``sweep`` at load time.
# The back-edge — ``transport_sweep`` reaching here for ``default_for`` — is a
# function-local (lazy) import on the ``sweep`` side, which breaks the cycle.
# S6.3: the matvec ``(L+C)ψ`` walks moved OFF the operator INTO the
# representations' ``loss_action`` — so ``_x_scan_faces`` (the shared row-march
# face primitive) is imported here too.
from .sweep import (
    _scanmarch_row,
    _sweep_1d_unified,
    _sweep_2d_scheduled,
    _sweep_2d_wavefront,
    _sweep_full_field,
    _x_scan_faces,
)
from .sweep_graph import OctantLabel, SweepDependencyGraph
from .sweep_schedule import SweepSchedule

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping

    from orpheus.numerics.projection import MomentProjection
    from orpheus.transport.fields.angular_flux import AngularFlux
    from orpheus.transport.fields.boundary_flux import BoundaryFlux
    from orpheus.transport.timed_full_field import TimedFullField

    from .geometry import SNMesh
    from .operator import StreamingOperator
    from .sweep_schedule import OctantSweep, OctantSweepGroup


# ═══════════════════════════════════════════════════════════════════════
# Selection vocabulary
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class Compatibility:
    """Whether a strategy applies to a mesh, with a human-readable reason.

    ``ok`` is the machine-checkable verdict; ``reason`` is the explanation a
    frontend shows when graying-out an inapplicable method ("Moving-frontier
    window — requires Cartesian geometry, d = 2") and the message the
    construction guard raises.  ``reason`` is the empty string when
    ``ok is True`` (no explanation needed).
    """

    ok: bool
    reason: str


class IncompatibleRepresentation(ValueError):
    """A :class:`LossRepresentation` was constructed for a mesh it cannot sweep.

    Raised by the construction guard (:meth:`_LossRepresentation.__post_init__`)
    so that an illegal ``(strategy, mesh)`` pairing is unrepresentable —
    even if a caller bypasses :func:`default_for`.
    """


@runtime_checkable
class LossRepresentation(Protocol):
    r"""One algorithm for the within-group transport solve and its twin.

    A strategy is constructed *for a mesh* (``Strategy(mesh)``); the
    construction guard rejects an incompatible pairing.  It then exposes:

    * :meth:`sweep` — one forward substitution :math:`\psi = (L+C)^{-1} q`;
    * :meth:`loss_action` — the matvec twin :math:`(L+C)\,\psi` *(added in
      Phase S2)*;
    * :meth:`supports` — the (classmethod) selection predicate.
    """

    def sweep(
        self,
        Q: "np.ndarray",
        sig_t: "np.ndarray",
        boundary_flux: "BoundaryFlux",
        *,
        initial_guess: "AngularFlux | TimedFullField | None" = None,
        moment_projection: "MomentProjection | None" = None,
    ) -> "tuple[np.ndarray, np.ndarray]":
        """Perform one within-group transport sweep on this strategy's mesh."""
        ...

    def loss_action(
        self, operator: "StreamingOperator", psi: "TimedFullField",
    ) -> "TimedFullField":
        r"""The forward within-group loss action :math:`(L+C)\,\psi` for this geometry.

        The sweep's operator-twin (L21 — sweep and matvec are different
        applications of the SAME operator): the sweep solves
        :math:`(L+C)^{-1} q`, this APPLIES :math:`(L+C)`.  **Return the FULL loss
        :math:`(L+C)\psi`, NOT :math:`L\psi`** — the operator
        (:meth:`~orpheus.sn.operator.StreamingOperator.apply`) subtracts the
        collision diagonal :math:`C = \sigma_t\odot` exactly ONCE (Resolution A
        :math:`L = (L+C) - C`).  A leaf that returned :math:`L\psi` would make the
        operator subtract :math:`C` a SECOND time (a double-counted collision
        diagonal).  ``operator`` supplies :math:`\sigma_t` and the per-geometry
        walk machinery.
        """
        ...

    def loss_action_transpose(
        self, operator: "StreamingOperator", phi: "TimedFullField",
    ) -> "TimedFullField":
        r"""The adjoint loss action :math:`(L+C)^{\mathsf T}\,\phi` for this geometry.

        Return the FULL adjoint loss :math:`(L+C)^{\mathsf T}\phi` (the operator
        subtracts the self-adjoint diagonal :math:`C` in
        :meth:`~orpheus.sn.operator.StreamingOperator.apply_transpose`).  Raises
        :class:`NotImplementedError` for representations whose adjoint is deferred
        (the multi-D Cartesian reverse sweep — O.2b lands the 1-D reverse sweep
        first).  Never a silent wrong answer.
        """
        ...

    @classmethod
    def supports(cls, mesh: "SNMesh") -> Compatibility:
        """Whether this strategy can sweep ``mesh`` (the selection layer)."""
        ...


# ═══════════════════════════════════════════════════════════════════════
# Common base — the construction guard (illegal pairings unrepresentable)
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class _LossRepresentation:
    """Base for every concrete strategy: holds the mesh + the guard.

    A frozen dataclass carrying the one piece of state every strategy needs
    (the :class:`SNMesh` it was selected for) and the construction guard
    that makes an incompatible pairing unrepresentable.
    """

    mesh: "SNMesh"

    @classmethod
    def supports(cls, mesh: "SNMesh") -> Compatibility:
        """The selection predicate — every concrete strategy implements it."""
        raise NotImplementedError(
            f"{cls.__name__} must implement supports()"
        )

    def sweep(
        self,
        Q: "np.ndarray",
        sig_t: "np.ndarray",
        boundary_flux: "BoundaryFlux",
        *,
        initial_guess: "AngularFlux | TimedFullField | None" = None,
        moment_projection: "MomentProjection | None" = None,
    ) -> "tuple[np.ndarray, np.ndarray]":
        """One within-group sweep — every concrete strategy implements it."""
        raise NotImplementedError(
            f"{type(self).__name__} must implement sweep()"
        )

    def loss_action(
        self, operator: "StreamingOperator", psi: "TimedFullField",
    ) -> "TimedFullField":
        """The forward loss action ``(L+C)ψ`` — every concrete leaf implements it.

        Returns the FULL within-group loss ``(L+C)ψ`` (NOT ``Lψ``); the operator
        subtracts ``C = σ_t⊙`` in ``apply`` (Resolution A ``L = (L+C) − C``).
        """
        raise NotImplementedError(
            f"{type(self).__name__} must implement loss_action()"
        )

    def loss_action_transpose(
        self, operator: "StreamingOperator", phi: "TimedFullField",
    ) -> "TimedFullField":
        """The adjoint loss action ``(L+C)ᵀφ`` — 1-D implemented or a deferral raise.

        Returns the FULL adjoint loss ``(L+C)ᵀφ`` (the operator subtracts ``C`` in
        ``apply_transpose``).
        """
        raise NotImplementedError(
            f"{type(self).__name__} must implement loss_action_transpose()"
        )

    def __post_init__(self) -> None:
        compat = type(self).supports(self.mesh)
        if not compat.ok:
            raise IncompatibleRepresentation(
                f"{type(self).__name__} cannot sweep this mesh "
                f"(ndim={self.mesh.ndim}, curvature={self.mesh.curvature!r}): "
                f"{compat.reason}."
            )


# ═══════════════════════════════════════════════════════════════════════
# _OctantWalk — THE in-plane octant traversal (S6.4, #222)
# ═══════════════════════════════════════════════════════════════════════

#: Spatial axis names, positional-by-axis — the same axis order as
#: :attr:`OctantLabel.signs` and every per-axis kernel tuple.
_AXIS_NAMES = ("x", "y", "z")


def _inflow_faces(signs_eff: tuple[int, ...]) -> tuple[str, ...]:
    """Per-axis domain faces an octant's streaming ENTERS through.

    An octant streaming in the ``+a`` direction enters at the ``a``-min face
    (``("xmin", "ymin")`` for the ``(+1, +1)`` octant); a ``−a`` octant at
    the ``a``-max face.  ``signs_eff`` carries the EFFECTIVE signs (grazing
    ``0`` already mapped to ``+1`` — the streaming coefficient is zero, so
    the WDD result is sign-independent).
    """
    return tuple(
        f"{_AXIS_NAMES[a]}min" if s >= 0 else f"{_AXIS_NAMES[a]}max"
        for a, s in enumerate(signs_eff)
    )


def _outflow_faces(signs_eff: tuple[int, ...]) -> tuple[str, ...]:
    """Per-axis domain faces an octant's streaming EXITS through (the
    opposite of :func:`_inflow_faces`, axis by axis)."""
    return tuple(
        f"{_AXIS_NAMES[a]}max" if s >= 0 else f"{_AXIS_NAMES[a]}min"
        for a, s in enumerate(signs_eff)
    )


@dataclass(frozen=True)
class _ApplyOperands:
    r"""Problem data of the APPLY direction :math:`(L+C)\,\psi`.

    What every apply-direction interior kernel consumes, bundled once per
    :meth:`_OctantWalk.loss_action` call.  ``probe`` is the apply target
    :math:`\bar\psi` (the matvec input).  ``Q_zero`` is the zero volumetric
    source — the matvec evaluates the loss *action*, not a balance; kernels
    whose walk signature requires a source slot (the windowed graph walk)
    pass it through.  ``str_axes`` is the per-axis streaming tuple
    :math:`2|\mu_a|/\Delta a` over ``range(ndim)`` — positional-by-axis like
    every kernel tuple, so axis ``a``'s coefficients pair with axis ``a``'s
    faces by construction.
    """

    probe: "np.ndarray"                  # (N, ng, *spatial) — ψ̄, the apply target
    sig_t: "np.ndarray"                  # (ng, *spatial)
    str_axes: tuple["np.ndarray", ...]   # d arrays, each (N, n_a)
    Q_zero: "np.ndarray"                 # (1, ng, *spatial)


@dataclass(frozen=True)
class _SolveOperands:
    r"""Problem data of the SOLVE direction :math:`(L+C)^{-1} q`.

    The sweep's mirror of :class:`_ApplyOperands`: the solve direction is
    driven by the GIVEN per-ordinate volumetric source ``Q`` (the unknown is
    :math:`\bar\psi`), where the apply direction is driven by the given probe
    :math:`\bar\psi` (no source).  Same positional-by-axis ``str_axes``
    convention.
    """

    Q: "np.ndarray"                      # (N, ng, *spatial) — per-ordinate source
    sig_t: "np.ndarray"                  # (ng, *spatial)
    str_axes: tuple["np.ndarray", ...]   # d arrays, each (N, n_a)


@dataclass(frozen=True)
class _SweepEmit:
    r"""Solve-direction OUTPUT mode — angular field XOR harmonic moments.

    The Phase 5c output DI (which buffers are given selects the mode —
    mirroring :meth:`SweepDependencyGraph.apply_windowed`'s contract), made
    a TYPE: the construction guard rejects a mixed or empty mode, so an
    illegal half-wired output is unrepresentable.

    * **angular** — ``angular_flux`` ``(N, ng, *spatial)`` written per
      octant + ``scalar_flux`` ``(ng, *spatial)`` accumulated
      :math:`\sum_n w_n \psi_n`.
    * **moment** — ``moment_buf`` ``(L+1, 2L+1, ng, *spatial)`` accumulated
      :math:`\phi_\ell^m \mathrel{+}= \sum_n w_n Y_\ell^m \psi_n` with the
      octant harmonics ``Y`` ``(N, L+1, 2L+1)``; the full angular field is
      never materialized (the ~3× peak-memory win; the scalar is subsumed,
      ``moment_buf[0, 0]``).

    The pure-z volumetric balance emits through :meth:`pure_z`; the
    interior kernels accumulate at their own granularity (per
    anti-hyperplane for the window, per row for the scan-march) reading the
    mode off these buffers.
    """

    weights: "np.ndarray"                       # (N,)
    angular_flux: "np.ndarray | None" = None    # (N, ng, *spatial)
    scalar_flux: "np.ndarray | None" = None     # (ng, *spatial)
    moment_buf: "np.ndarray | None" = None      # (L+1, 2L+1, ng, *spatial)
    Y: "np.ndarray | None" = None               # (N, L+1, 2L+1)

    def __post_init__(self) -> None:
        angular = (self.angular_flux is not None) and (self.scalar_flux is not None)
        moment = (self.moment_buf is not None) and (self.Y is not None)
        if angular == moment:
            raise ValueError(
                "_SweepEmit: exactly ONE output mode must be wired — either "
                "(angular_flux AND scalar_flux) or (moment_buf AND Y)."
            )

    def pure_z(self, oct_idx: "np.ndarray", psi_avg: "np.ndarray") -> None:
        """Emit the pure-z volumetric balance ``ψ = Q/Σ_t`` (no faces).

        The accumulations use ``buf[...] +=`` (item-level in-place add, the
        same ufunc as a bare ``+=``) — a bare ``self.buf +=`` would rebind
        the attribute and trip the frozen dataclass.
        """
        if self.moment_buf is None:
            self.angular_flux[oct_idx] = psi_avg
            self.scalar_flux[...] += np.einsum(
                "ng...,n->g...", psi_avg, self.weights[oct_idx],
            )
        else:
            self.moment_buf[...] += np.einsum(
                "nlm,ng...,n->lmg...", self.Y[oct_idx], psi_avg,
                self.weights[oct_idx],
            )


@dataclass(frozen=True)
class _OctantWalk:
    r"""THE in-plane octant traversal of the Cartesian loss operator.

    The sweep (forward substitution :math:`(L+C)^{-1} q`) and the matvec
    (:math:`(L+C)\,\psi`) traverse the SAME octant decomposition: project the
    quadrature octant to its in-plane signs, branch the pure-z degenerate
    octants, derive the per-axis in/out domain faces, read the octant's
    inflow, run the interior traversal, shed the outflow.  The two
    directions fork ONLY at

    * the **cell kernel** — the per-octant interior traversal the calling
      representation supplies: the window's frontier walk, the scan-march's
      row-march, the oracle's full cochain — each in its solve
      (:meth:`~orpheus.sn.spatial.diamond.DiamondDifference.cell_kernel_batch`)
      or apply
      (:meth:`~orpheus.sn.spatial.diamond.DiamondDifference.residual_kernel_batch`)
      direction; and
    * the **emit policy** — what the direction accumulates: the sweep's
      angular/moment output; the matvec's :math:`(L+C)\psi` bulk + the O.4b
      boundary defect.

    NEVER a boolean ``is_solve`` flag — the direction is carried by the
    kernel/emit OBJECTS (the anti-degradation tripwire in
    ``tests/sn/operators/test_one_octant_walk.py`` enforces this shape).

    Dimension-generic from birth: signs / faces / inflow / captures are
    per-axis tuples over ``mesh.ndim`` — at d = 2 byte-identical to the
    legacy x/y spelling (pinned by the ``window ≡ full`` oracles).

    S6.4 staging: sub-step (a) routes the window + scan-march MATVEC frames
    through this walk; (b) brings the sweep frames in (the one-walk spy
    test flips xfail → xpass); (d) folds the full-field oracle.
    """

    mesh: "SNMesh"

    def _interior_walk(
        self,
        sweeps: "tuple[OctantSweep, ...]",
        *,
        inflow_of: "Callable[[str], np.ndarray]",
        shed: "Callable[[str, np.ndarray, np.ndarray], None]",
        pure_z: "Callable[[np.ndarray], None]",
        interior: "Callable[[np.ndarray, tuple[int, ...], tuple[np.ndarray, ...]], tuple[np.ndarray, ...]]",
    ) -> None:
        r"""THE shared octant frame (the one-walk seam, S6.4).

        For each octant sweep unit: project the label to the in-plane signs,
        dispatch pure-z degenerates to ``pure_z``, derive the effective signs
        and per-axis domain faces, read the octant's inflow via
        ``inflow_of(face)[oct_idx]``, run ``interior`` (returning the
        per-axis outflow captures), and ``shed`` each capture into its
        outflow face.  Both public directions route through here — the
        matvec since sub-step (a), the sweep from sub-step (b) — so
        "matvec ≡ sweep is one walk" is a code fact, not a test-maintained
        coincidence.
        """
        ndim = self.mesh.ndim
        for sweep in sweeps:
            oct_idx = np.asarray(sweep.indices)
            signs = sweep.label.signs[:ndim]
            if not any(signs):
                # Pure-z degenerate octant: no in-plane streaming — no
                # faces, no boundary interaction. The direction's policy
                # handles the volumetric balance.
                pure_z(oct_idx)
                continue
            # Grazing (sign 0) ordinates ride the +1 sweep direction: the
            # streaming coefficient is zero, the WDD result sign-independent
            # (matches the legacy sx_eff/sy_eff mapping).
            signs_eff = tuple(+1 if s == 0 else s for s in signs)
            inflow = tuple(
                inflow_of(face)[oct_idx] for face in _inflow_faces(signs_eff)
            )
            capture = interior(oct_idx, signs_eff, inflow)
            for face, capture_a in zip(_outflow_faces(signs_eff), capture):
                shed(face, oct_idx, capture_a)

    def sweep_group(
        self,
        group: "OctantSweepGroup",
        *,
        operands: _SolveOperands,
        emit: _SweepEmit,
        boundary_flux: "BoundaryFlux",
        interior: "Callable[[_SolveOperands, _SweepEmit, np.ndarray, tuple[int, ...], tuple[np.ndarray, ...]], tuple[np.ndarray, ...]]",
    ) -> None:
        r"""The SOLVE-direction frame for ONE octant group (S6.4 sub-step (b)).

        One forward-substitution pass over the group's octants on the SAME
        :meth:`_interior_walk` frame the matvec uses — the L21 unification.
        The Jacobi / Gauss-Seidel splitting lives one level up (the schedule
        loop's inter-group reflect, :func:`~orpheus.sn.sweep._sweep_2d_scheduled`);
        this frame is the bare per-group sweep, blind to the boundary
        coupling.  The calling representation supplies ONLY its interior
        kernel::

            interior(operands, emit, oct_idx, signs_eff, inflow) -> capture

        Boundary coupling via the LIVE ``boundary_flux``: each octant reads
        its inflow off the trace and sheds its outflow back into it as the
        walk advances.  Distinct octants own DISJOINT ordinate slices of a
        face, so an octant's outflow write never clobbers another octant's
        inflow — the Jacobi single-group call is bit-identical to the legacy
        per-octant loop, and the Gauss-Seidel schedule reflects the
        just-shed outflow between groups so a later group reads the fresh
        current-iterate inflow off the SAME trace (the
        :math:`(L+C-B_{\rm lower})^{-1}` forward substitution).

        The pure-z degenerate octants take the volumetric balance
        :math:`\psi = Q_n / \Sigma_t` straight into the emit policy — no
        faces, no boundary interaction.
        """
        def pure_z(oct_idx: "np.ndarray") -> None:
            emit.pure_z(oct_idx, operands.Q[oct_idx] / operands.sig_t)

        def run_interior(
            oct_idx: "np.ndarray",
            signs_eff: tuple[int, ...],
            inflow: tuple["np.ndarray", ...],
        ) -> tuple["np.ndarray", ...]:
            return interior(operands, emit, oct_idx, signs_eff, inflow)

        def shed(face: str, oct_idx: "np.ndarray", capture_a: "np.ndarray") -> None:
            boundary_flux.face_view(face)[oct_idx] = capture_a

        self._interior_walk(
            group.sweeps,
            inflow_of=boundary_flux.face_view,
            shed=shed,
            pure_z=pure_z,
            interior=run_interior,
        )

    def loss_action(
        self,
        operator: "StreamingOperator",
        psi: "TimedFullField",
        interior: "Callable[[_ApplyOperands, np.ndarray, tuple[int, ...], tuple[np.ndarray, ...]], tuple[np.ndarray, tuple[np.ndarray, ...]]]",
    ) -> "TimedFullField":
        r"""The APPLY-direction frame :math:`(L+C)\,\psi` (S6.4 sub-step (a)).

        Owns everything the matvec frames previously duplicated in lockstep:
        the probe / accumulator setup, the pure-z branch
        (:math:`(L+C)\bar\psi = \Sigma_t\,\bar\psi` — no in-plane streaming,
        so :math:`L\bar\psi = 0` after the operator's :math:`-C`), the
        per-octant inflow read from the GIVEN trace, the outflow capture,
        the O.4b active-trace boundary residual, and the typed assembly.
        The calling representation supplies ONLY its interior kernel::

            interior(operands, oct_idx, signs_eff, inflow)
                -> (LpC_octant, capture)

        Boundary semantics (BARE — O.4b Phase E): each octant reads its
        inflow from the GIVEN trace ``psi.boundary`` (NO ``bc.apply`` — the
        reflective coupling is the sibling ``−B``); the domain-edge outflow
        is captured into ``streamed`` (OUTFLOW slots only).  The output
        boundary is the O.4b active-trace residual: OUTFLOW slots → defect
        ``streamed − given``; INFLOW slots → identity ``given``.

        Returns the FULL loss :math:`(L+C)\bar\psi` (NOT :math:`L\bar\psi`);
        :meth:`~orpheus.sn.operator.StreamingOperator.apply` subtracts
        :math:`\Sigma_t\,\bar\psi` exactly once (Resolution A).
        """
        from orpheus.transport.source_sinks import (
            AngularSourceSink, BoundarySourceSink,
        )
        from orpheus.transport.timed_full_field import TimedFullField

        sn_mesh = self.mesh
        ndim = sn_mesh.ndim
        sig_t = operator.sigma_t
        ng = sig_t.shape[0]
        spatial = sig_t.shape[1:]
        probe = psi.bulk.values
        operands = _ApplyOperands(
            probe=probe,
            sig_t=sig_t,
            str_axes=tuple(sn_mesh.streaming(a) for a in range(ndim)),
            Q_zero=np.zeros((1, ng, *spatial)),
        )

        # (L+C)·ψ̄ accumulator; ``apply`` subtracts Σ_t·ψ̄ → bare-streaming Lψ̄.
        LpC = np.zeros((sn_mesh.quad.N, ng, *spatial))
        trace = sn_mesh.trace
        boundary = psi.boundary
        streamed = {
            face: np.zeros_like(boundary.face_view(face))
            for face in trace.face_names
        }

        def pure_z(oct_idx: "np.ndarray") -> None:
            # (L+C)·ψ̄ = Σ_t·ψ̄ for the in-plane-degenerate ordinates.
            LpC[oct_idx] = sig_t * probe[oct_idx]

        def run_interior(
            oct_idx: "np.ndarray",
            signs_eff: tuple[int, ...],
            inflow: tuple["np.ndarray", ...],
        ) -> tuple["np.ndarray", ...]:
            LpC_oct, capture = interior(operands, oct_idx, signs_eff, inflow)
            LpC[oct_idx] = LpC_oct
            return capture

        def shed(face: str, oct_idx: "np.ndarray", capture_a: "np.ndarray") -> None:
            streamed[face][oct_idx] = capture_a

        (jacobi_group,) = SweepSchedule.jacobi(sn_mesh).groups
        self._interior_walk(
            jacobi_group.sweeps,
            inflow_of=boundary.face_view,
            shed=shed,
            pure_z=pure_z,
            interior=run_interior,
        )

        # Boundary-block residual (O.4b — the active trace).
        out_boundary = BoundarySourceSink.zeros_on(sn_mesh)
        for face in trace.face_names:
            given = boundary.face_view(face)
            out_idx = trace.outflow_indices_for_face(face)
            in_idx = trace.inflow_indices_for_face(face)
            if out_idx.size:
                out_boundary.face_view(face)[out_idx] = (
                    streamed[face][out_idx] - given[out_idx]
                )
            if in_idx.size:
                out_boundary.face_view(face)[in_idx] = given[in_idx]

        return TimedFullField(
            bulk=AngularSourceSink.from_mesh(LpC, sn_mesh),
            boundary=out_boundary,
            _history=(),
            history_depth=psi.history_depth,
        )


# ═══════════════════════════════════════════════════════════════════════
# CumprodScan — the 1-D chain prefix scan (any geometry)
# ═══════════════════════════════════════════════════════════════════════


class CumprodScan(_LossRepresentation):
    r"""1-D parallel-prefix scan — slab, sphere, cylinder via one body.

    Intrinsically 1-D: a prefix scan needs a total order (a chain).  The
    geometry difference is absorbed by the two-stratum cache, so slab +
    sphere + cylinder share THE SAME scan expression
    (:func:`._sweep_1d_unified` → :func:`~orpheus.sn.spatial.scan.ordinate_scan`).
    The default production path for every 1-D mesh.
    """

    @classmethod
    def supports(cls, mesh: "SNMesh") -> Compatibility:
        return Compatibility(mesh.is_1d, "requires a 1-D mesh")

    def sweep(
        self,
        Q: "np.ndarray",
        sig_t: "np.ndarray",
        boundary_flux: "BoundaryFlux",
        *,
        initial_guess: "AngularFlux | TimedFullField | None" = None,
        moment_projection: "MomentProjection | None" = None,
    ) -> "tuple[np.ndarray, np.ndarray]":
        if moment_projection is not None:
            # Moment output is the 2-D windowed-SI peak-memory optimization;
            # 1-D / curvilinear meshes stay full-angular (the Morel–Montry
            # Carlson seed reads the per-ordinate iterate; lesson L21).
            raise ValueError(
                "CumprodScan.sweep: moment output (moment_projection given) "
                "is 2-D Cartesian only — 1-D/curvilinear meshes stay "
                "full-angular (the Morel–Montry Carlson seed reads the "
                "per-ordinate iterate; lesson L21)."
            )
        return _sweep_1d_unified(
            Q, sig_t, self.mesh, boundary_flux, initial_guess=initial_guess,
        )

    def loss_action(
        self, operator: "StreamingOperator", psi: "TimedFullField",
    ) -> "TimedFullField":
        r"""1-D forward loss action ``(L+C)ψ`` — the geometry-blind spatial sum.

        S6.3: returns the FULL within-group loss ``(L+C)ψ``; the operator
        subtracts the collision diagonal ``C = σ_t⊙`` in :meth:`apply`.  The
        1-D ``(L+C)`` action IS the operator's spatial operator-sum
        ``operator.M_spatial`` (per-direction streaming + collision-share); the
        angular Morel–Montry redistribution rides inside ``_compute_LpC``.
        """
        return operator.M_spatial._compute_LpC(psi)

    def loss_action_transpose(
        self, operator: "StreamingOperator", phi: "TimedFullField",
    ) -> "TimedFullField":
        r"""1-D adjoint loss action ``(L+C)ᵀφ`` — the reverse spatial sum.

        S6.3: returns ``(L+C)ᵀφ`` (the operator subtracts ``C`` in
        :meth:`apply_transpose`).  ``operator.M_spatial._compute_LpC_transpose``
        carries the curvilinear angular SECOND triangular factor
        (``closure.angular_adjoint``) — so the spatial reverse NEVER silently
        drops the angular adjoint (pinned by ``test_g_adjoint_reciprocity``
        sphere/cyl, -O-firing since S6.3a).
        """
        return operator.M_spatial._compute_LpC_transpose(phi)


# ═══════════════════════════════════════════════════════════════════════
# _DAGWavefront — the Cartesian anti-hyperplane DAG family
# ═══════════════════════════════════════════════════════════════════════


class _DAGWavefront(_LossRepresentation):
    r"""Base for the two buffer policies over the per-octant DAG walk.

    ``FullFieldWavefront`` (full-field buffer; the oracle) and
    ``MovingFrontierWindow`` (rolling :math:`(d{-}1)`-frontier; the
    production optimization) both walk the **same** per-octant
    anti-hyperplane DAG (:attr:`sweep_graphs`) with the same
    diamond-difference cell kernel.  They differ only in *how much* of the
    interior face cochain they retain — a storage policy, pinned
    bit-identical by the ``window ≡ full`` oracle.

    S6.4(c) — the family OWNS the DAG: ``sweep_graphs`` is THIS base's
    accessor over the per-shape cache
    :meth:`SweepDependencyGraph.for_shape`, NOT a mesh attribute.  The mesh
    is pure geometry; DAG-free representations (:class:`CumprodScan`,
    :class:`ScanMarch`) never mention the substrate — the historical
    curvilinear ``mesh.sweep_graphs = None`` slot (an illegal state) is
    unrepresentable.

    The DAG walk is naturally d-general; in S1 both wrappers cover exactly
    the existing 2-D Cartesian sweep (``supports`` below).  S3 widens the
    oracle to any-d Cartesian; S4 widens the window to ``d ≥ 2``.
    """

    @classmethod
    def supports(cls, mesh: "SNMesh") -> Compatibility:
        return Compatibility(
            mesh.is_cartesian and mesh.ndim == 2,
            "requires Cartesian geometry, d = 2",
        )

    @property
    def sweep_graphs(self) -> "Mapping[OctantLabel, SweepDependencyGraph]":
        r"""The per-octant DAG family for this representation's mesh shape.

        Routes to the per-shape cache (the graphs depend only on cell
        topology + octant signs, so same-shape meshes share byte-identical
        graphs); treat the mapping as immutable.
        """
        return SweepDependencyGraph.for_shape(self.mesh.spatial_shape)

    def loss_action_transpose(
        self, operator: "StreamingOperator", phi: "TimedFullField",
    ) -> "TimedFullField":
        r"""The multi-D Cartesian adjoint is DEFERRED (shared by both policies).

        The forward matvec works for both buffer policies, but the
        reverse-direction adjoint sweep is not yet wired (O.2b landed the 1-D
        reverse sweep first).  Raises :class:`NotImplementedError` — the mesh
        is compatible, only the adjoint *feature* is deferred (so this is NOT
        an :class:`IncompatibleRepresentation`).  Never a silent wrong answer.
        """
        raise NotImplementedError(
            "StreamingOperator.apply_transpose: the multi-D Cartesian adjoint "
            "is deferred (O.2b lands the 1-D reverse sweep first)."
        )


class MovingFrontierWindow(_DAGWavefront):
    r"""Production wavefront sweep — rolling :math:`(d{-}1)`-frontier buffer.

    The default production path for 2-D Cartesian meshes
    (:func:`._sweep_2d_wavefront`).  Carries only the rolling frontier of
    interior face fluxes (a 2-diagonal at d=2), the ~3× peak-memory win over
    the full-field oracle.  Generalized to ``frontier_dim = d-1`` in S4.
    """

    def sweep(
        self,
        Q: "np.ndarray",
        sig_t: "np.ndarray",
        boundary_flux: "BoundaryFlux",
        *,
        initial_guess: "AngularFlux | TimedFullField | None" = None,
        moment_projection: "MomentProjection | None" = None,
    ) -> "tuple[np.ndarray, np.ndarray]":
        return _sweep_2d_wavefront(
            Q, sig_t, self.mesh, boundary_flux,
            moment_projection=moment_projection,
            interior=self._sweep_interior,
        )

    def _sweep_interior(
        self,
        operands: _SolveOperands,
        emit: _SweepEmit,
        oct_idx: "np.ndarray",
        signs_eff: tuple[int, ...],
        inflow: tuple["np.ndarray", ...],
    ) -> tuple["np.ndarray", ...]:
        r"""Rolling-frontier interior kernel, SOLVE direction, one octant.

        Drives
        :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.apply_windowed`
        (the windowed walk of the solve cell kernel
        :meth:`~orpheus.sn.spatial.diamond.DiamondDifference.cell_kernel_batch`)
        over this octant's DAG, emitting per anti-hyperplane into the
        :class:`_SweepEmit` mode buffers.  Returns the per-axis domain-edge
        outflow ``capture``.
        """
        graph = self.sweep_graphs[OctantLabel(signs_eff)]
        ng = operands.sig_t.shape[0]
        spatial = operands.sig_t.shape[1:]
        capture = tuple(np.empty_like(face) for face in inflow)
        # Angular mode allocates a per-octant angular buffer (scattered into
        # the global field below); moment mode accumulates directly into the
        # shared moment tensor per anti-hyperplane, so NO per-octant angular
        # field is materialized (the Phase 5c peak-memory win).
        angular_flux_oct = (
            np.zeros((oct_idx.size, ng, *spatial))
            if emit.moment_buf is None else None
        )
        graph.apply_windowed(
            cell_update=self.mesh.cell_update,
            inflow=inflow,
            Q_octant=operands.Q[oct_idx],
            sig_t=operands.sig_t,
            str_axes_octant=tuple(s[oct_idx] for s in operands.str_axes),
            weights_octant=emit.weights[oct_idx],
            capture=capture,
            angular_flux_octant=angular_flux_oct,
            scalar_flux_buf=emit.scalar_flux,
            moment_buf=emit.moment_buf,
            Y_octant=None if emit.Y is None else emit.Y[oct_idx],
        )
        if emit.moment_buf is None:
            emit.angular_flux[oct_idx] = angular_flux_oct
        return capture

    def loss_action(
        self, operator: "StreamingOperator", psi: "TimedFullField",
    ) -> "TimedFullField":
        r"""2-D Cartesian forward loss action ``(L+C)ψ`` via the rolling-frontier window.

        S6.4 sub-step (a): routes through the shared :class:`_OctantWalk`
        apply frame (the ONE octant traversal — octant projection, pure-z
        branch, boundary I/O, the O.4b boundary residual), supplying only the
        rolling-frontier interior kernel :meth:`_loss_action_interior`
        (:meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.residual_windowed`
        — the apply-direction walk of the SAME per-octant wavefront DAG and
        the SAME diamond-difference closure the 2-D sweep uses; matvec ≡
        sweep, ONE discretization, L21).  Returns ``(L+C)ψ̄``;
        :meth:`~orpheus.sn.operator.StreamingOperator.apply` subtracts
        ``Σ_t·ψ̄`` (the collision diagonal ``C``) to recover the
        bare-streaming ``Lψ̄``.
        """
        return _OctantWalk(self.mesh).loss_action(
            operator, psi, self._loss_action_interior,
        )

    def _loss_action_interior(
        self,
        operands: _ApplyOperands,
        oct_idx: "np.ndarray",
        signs_eff: tuple[int, ...],
        inflow: tuple["np.ndarray", ...],
    ) -> tuple["np.ndarray", tuple["np.ndarray", ...]]:
        r"""Rolling-frontier interior kernel, APPLY direction, one octant.

        Drives
        :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.residual_windowed`
        (the windowed walk of the apply cell kernel
        :meth:`~orpheus.sn.spatial.diamond.DiamondDifference.residual_kernel_batch`)
        over this octant's DAG.  Returns ``(LpC_octant, capture)`` — the
        octant's ``(L+C)ψ̄`` block and the per-axis domain-edge outflow.
        """
        graph = self.sweep_graphs[OctantLabel(signs_eff)]
        ng = operands.sig_t.shape[0]
        spatial = operands.sig_t.shape[1:]
        LpC_oct = np.zeros((oct_idx.size, ng, *spatial))
        capture = tuple(np.empty_like(face) for face in inflow)
        graph.residual_windowed(
            cell_update=self.mesh.cell_update,
            inflow=inflow,
            psi_avg_probe_octant=operands.probe[oct_idx],
            Q_octant=operands.Q_zero,
            sig_t=operands.sig_t,
            str_axes_octant=tuple(s[oct_idx] for s in operands.str_axes),
            residual_octant=LpC_oct,
            capture=capture,
        )
        return LpC_oct, capture


class FullFieldWavefront(_DAGWavefront):
    r"""Verification-oracle wavefront sweep — the dimension-generic SPINE.

    Walks the same per-octant DAG as :class:`MovingFrontierWindow` but
    retains the FULL interior face cochain (the fuller view).  Slower and
    more memory-hungry — its purpose is verification: ONE body for d=1 (slab)
    and d=2 (Cartesian), and the reference the d-specific production
    optimizations are cross-checked against — the 1-D :class:`CumprodScan`
    (principled-equivalence at nulp) and the 2-D :class:`MovingFrontierWindow`
    (``window ≡ full`` bit-identity).  Never the production default (the
    window wins at d=2, the scan at d=1); selected explicitly by oracle tests.

    Walks the d-generic :func:`._sweep_full_field` (sweep); since S6.3 it holds
    the matvec ``(L+C)ψ`` walk in its own :meth:`loss_action` (both walk
    :meth:`SweepDependencyGraph.residual` directly).  ``supports`` is any-d
    Cartesian (S3) — the spine is genuinely dimension-generic, unlike the
    d=2-only window.
    """

    @classmethod
    def supports(cls, mesh: "SNMesh") -> Compatibility:
        # Override the _DAGWavefront family's d=2-only predicate: the spine is
        # the genuine d-generic oracle (it walks the per-octant DAG for any
        # Cartesian d via the d-generic ``graph.residual``).
        return Compatibility(mesh.is_cartesian, "requires Cartesian geometry")

    def sweep(
        self,
        Q: "np.ndarray",
        sig_t: "np.ndarray",
        boundary_flux: "BoundaryFlux",
        *,
        initial_guess: "AngularFlux | TimedFullField | None" = None,
        moment_projection: "MomentProjection | None" = None,
    ) -> "tuple[np.ndarray, np.ndarray]":
        if moment_projection is not None:
            raise ValueError(
                "FullFieldWavefront.sweep: the full-field oracle does not "
                "implement moment output — use MovingFrontierWindow for the "
                "windowed-SI moment path."
            )
        return _sweep_full_field(Q, sig_t, self.mesh, boundary_flux)

    def loss_action(
        self, operator: "StreamingOperator", psi: "TimedFullField",
    ) -> "TimedFullField":
        r"""Forward loss action ORACLE ``(L+C)ψ`` — the full-field DAG walk (d-generic).

        S6.3: the verification-oracle matvec ``(L+C)ψ`` walk LIVES here (moved off
        the operator).  The matvec counterpart of
        :func:`~orpheus.sn.sweep._sweep_full_field`: it carries the FULL interior
        cochain as a typed
        :class:`~orpheus.transport.fields.wavefront_flux.WavefrontFlux` and walks
        via :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.residual` (the
        full-field walk sharing the SAME cell kernel as ``residual_windowed``), so
        the MATH cannot drift from :meth:`MovingFrontierWindow.loss_action` — only
        storage.  Returns ``(L+C)ψ̄``;
        :meth:`~orpheus.sn.operator.StreamingOperator.apply` subtracts ``Σ_t·ψ̄``.
        Sole purpose: verification (production is the window / the 1-D scan).

        Dimension-generic: the per-axis tuples are built over ``range(ndim)``
        (``streaming(a)`` / ``WavefrontFlux.face(a)``), and the quadrature octant
        label (which lives in the full angular space) is projected to the in-plane
        ``label[:ndim]``.
        """
        from orpheus.sn.sweep_graph import OctantLabel
        from orpheus.transport.fields.wavefront_flux import WavefrontFlux
        from orpheus.transport.source_sinks import (
            AngularSourceSink, BoundarySourceSink,
        )
        from orpheus.transport.timed_full_field import TimedFullField

        sn_mesh = operator.sn_mesh
        quad = sn_mesh.quad
        N = quad.N
        ndim = sn_mesh.ndim
        sig_t = operator.sigma_t
        ng = sig_t.shape[0]
        spatial = sig_t.shape[1:]                        # (nx,) d=1; (nx, ny) d=2
        probe = psi.bulk.values
        cell_update = sn_mesh.cell_update
        Q_zero = np.zeros((1, ng, *spatial))

        # (L+C)·ψ̄ accumulator; ``apply`` subtracts Σ_t·ψ̄ → Lψ̄.
        LpC = np.zeros((N, ng, *spatial))
        trace = sn_mesh.trace
        boundary = psi.boundary
        # The FULL interior cochain (typed) — ι_*-seeded whole-trace.  The
        # per-axis face tuple is born from the cochain's OWN axis map
        # (``WavefrontFlux.face(a)`` over ``WavefrontFlux.axes``); the streaming
        # tuple is the axis-keyed ``sn_mesh.streaming(a)`` over the SAME
        # ``range(ndim)`` — so ``str_axes[a]`` cannot drift from ``psi_faces[a]``.
        wavefront = WavefrontFlux.zeros_on(sn_mesh)
        psi_faces = tuple(wavefront.face(a) for a in wavefront.axes)
        str_axes = tuple(sn_mesh.streaming(a) for a in range(ndim))
        wavefront.seed(boundary)

        for octant in quad.octants:
            oct_idx = octant.indices
            signs = octant.label[:ndim]                  # in-plane projection
            if not any(signs):                           # pure-z degenerate (d≥2 only)
                LpC[oct_idx] = sig_t * probe[oct_idx]
                continue
            signs_eff = tuple(+1 if s == 0 else s for s in signs)
            graph = self.sweep_graphs[OctantLabel(signs_eff)]
            psi_faces_oct = tuple(pf[oct_idx].copy() for pf in psi_faces)
            LpC_oct = np.zeros((oct_idx.size, ng, *spatial))
            graph.residual(
                cell_update=cell_update,
                psi_faces_octant=psi_faces_oct,
                psi_avg_probe_octant=probe[oct_idx],
                Q_octant=Q_zero, sig_t=sig_t,
                str_axes_octant=tuple(s[oct_idx] for s in str_axes),
                residual_octant=LpC_oct,
            )
            for a in wavefront.axes:
                psi_faces[a][oct_idx] = psi_faces_oct[a]
            LpC[oct_idx] = LpC_oct

        # The post-walk domain-edge outflow (full interior cochain edge slots).
        streamed = {face: wavefront.edge_view(face) for face in trace.face_names}
        out_boundary = BoundarySourceSink.zeros_on(sn_mesh)
        for face in trace.face_names:
            given = boundary.face_view(face)
            out_idx = trace.outflow_indices_for_face(face)
            in_idx = trace.inflow_indices_for_face(face)
            if out_idx.size:
                out_boundary.face_view(face)[out_idx] = (
                    streamed[face][out_idx] - given[out_idx]
                )
            if in_idx.size:
                out_boundary.face_view(face)[in_idx] = given[in_idx]

        return TimedFullField(
            bulk=AngularSourceSink.from_mesh(LpC, sn_mesh),
            boundary=out_boundary,
            _history=(),
            history_depth=psi.history_depth,
        )


# ═══════════════════════════════════════════════════════════════════════
# ScanMarch — the row-march + x-scan schedule (1-D scan ∘ transverse march)
# ═══════════════════════════════════════════════════════════════════════


class ScanMarch(_LossRepresentation):
    r"""Scan-march sweep — ``scan(x)`` marched over the transverse axes (#222).

    Reframes the d-D diamond-difference sweep as forward substitution along the
    sweep axis — the first-order linear scan
    :func:`~orpheus.sn.spatial.scan.ordinate_scan` — marched over the transverse
    axes: ``scan(x)`` at d=1, ``scan(x) ∘ march(y)`` at d=2.  ONE primitive that
    **unifies** the 1-D :class:`CumprodScan` (its degenerate ``s_y = 0`` case)
    and the 2-D row-march: the within-row x-face recurrence is the SAME Blelloch
    scan, the transverse coupling rides the affine source (the row-march
    interior kernel :meth:`_sweep_interior`, S6.4(b) — the former private
    ``_sweep_2d_scanmarch`` frame dissolved into the shared
    :class:`_OctantWalk` + the Jacobi schedule).

    A different *schedule* from the :class:`_DAGWavefront` family (row-march vs
    anti-diagonal) over the SAME lower-triangular solve — principled-equivalent
    at nulp, pinned against the :class:`FullFieldWavefront` oracle (issue #222).
    Its production value: it reuses the conditioning-robust ``ordinate_scan``
    per line (the ERR-054 pole reset + the ERR-057 denormal underflow handled
    for free) and is the natural home for the flux-independent ``a_attenuation``
    cache the wavefront lacks (#206).

    Selection — ``is_1d OR is_cartesian``: 1-D any geometry (the chain scan; the
    curvilinear Morel–Montry angular thread folds into the source) AND Cartesian
    any d (the row-march).  Currently **OPT-IN**: registered after the production
    window, so :func:`default_for` keeps the legacy choice (1-D → ``CumprodScan``,
    2-D → ``MovingFrontierWindow``) until a measured speedup justifies the flip
    (``scan_march_verification.md`` Fork B1 → B2).  The 2-D matvec twin
    (``loss_action``) + the d≥3 recursive transverse march land in S5.1b/S5.2.
    """

    @classmethod
    def supports(cls, mesh: "SNMesh") -> Compatibility:
        return Compatibility(
            mesh.is_1d or mesh.is_cartesian,
            "requires a 1-D mesh (any geometry) or Cartesian geometry",
        )

    def sweep(
        self,
        Q: "np.ndarray",
        sig_t: "np.ndarray",
        boundary_flux: "BoundaryFlux",
        *,
        initial_guess: "AngularFlux | TimedFullField | None" = None,
        moment_projection: "MomentProjection | None" = None,
    ) -> "tuple[np.ndarray, np.ndarray]":
        if self.mesh.is_1d:
            # d=1 ⇒ ``scan(x)`` with no transverse march: the unified 1-D body
            # (slab + curvilinear via the two-stratum cache; the Morel–Montry
            # Carlson angular thread folds into the scan's affine source).  This
            # is the ``s_y = 0`` degeneration of the 2-D scan-march.
            if moment_projection is not None:
                raise ValueError(
                    "ScanMarch.sweep: moment output (moment_projection given) "
                    "is 2-D Cartesian only — 1-D/curvilinear meshes stay "
                    "full-angular (the Morel–Montry Carlson seed reads the "
                    "per-ordinate iterate; lesson L21)."
                )
            return _sweep_1d_unified(
                Q, sig_t, self.mesh, boundary_flux, initial_guess=initial_guess,
            )
        # 2-D ⇒ the row-march sweep = the Jacobi schedule × the scan-march
        # interior kernel on the SAME schedule loop the window uses (S6.4(b):
        # the former private ``_sweep_2d_scanmarch`` frame dissolved into the
        # shared walk — and the Gauss-Seidel schedule composes for free, the
        # inter-group reflect being kernel-agnostic).
        return _sweep_2d_scheduled(
            Q, sig_t, self.mesh, boundary_flux,
            schedule=SweepSchedule.jacobi(self.mesh),
            reflect=None,
            moment_projection=moment_projection,
            interior=self._sweep_interior,
        )

    def _sweep_interior(
        self,
        operands: _SolveOperands,
        emit: _SweepEmit,
        oct_idx: "np.ndarray",
        signs_eff: tuple[int, ...],
        inflow: tuple["np.ndarray", ...],
    ) -> tuple["np.ndarray", ...]:
        r"""Row-march interior kernel, SOLVE direction, one octant.

        Marches the y-rows in the octant's y-sweep order: within each row the
        diamond-difference x-face recurrence is the first-order linear scan
        (:func:`~orpheus.sn.sweep._scanmarch_row` with the *solve*
        coefficients ``α = 2s_x/D − 1``, ``β = 2(Q + s_y ψ_{y,in})/D``), the
        transverse-y coupling riding the affine source.  Emits per row into
        the :class:`_SweepEmit` mode buffers.  Returns the per-axis
        domain-edge outflow ``(capture_x, out_y)``.

        The flux-independent ``α``/``D`` are computed PER LINE (the
        :math:`(d{-}1)`-slab working set; mesh-memoising them — the #206
        single-source cache — is the measured follow-on).
        """
        sig_t = operands.sig_t                      # (ng, nx, ny)
        ng, nx, ny = sig_t.shape
        sx_eff, sy_eff = signs_eff
        inflow_x, inflow_y = inflow                 # (N_oct, ng, ny) / (N_oct, ng, nx)
        s_x, s_y = (s[oct_idx] for s in operands.str_axes)
        Q_oct = operands.Q[oct_idx]                 # (N_oct, ng, nx, ny)
        w_oct = emit.weights[oct_idx]               # (N_oct,)
        Y_oct = None if emit.Y is None else emit.Y[oct_idx]
        N_oct = oct_idx.size

        x_reverse = sx_eff < 0
        capture_x = np.empty((N_oct, ng, ny))       # domain x-outflow, per y-row
        angular_oct = (
            np.zeros((N_oct, ng, nx, ny)) if emit.moment_buf is None else None
        )

        # March the y-rows in the octant's y-sweep order, threading ψ_y.
        psi_y_in = inflow_y                          # (N_oct, ng, nx) — row-0 inflow
        out_y = psi_y_in                             # last-row out_y (ny ≥ 1 → set below)
        y_rows = range(ny) if sy_eff >= 0 else range(ny - 1, -1, -1)
        for j in y_rows:
            # D = σ_t + s_x + s_y on this row; the cell-kernel left-fold order
            # ((σ_t + s_x) + s_y) is bit-id-load-bearing WITHIN a schedule only.
            D_row = (
                sig_t[None, :, :, j]                 # (1, ng, nx)
                + s_x[:, None, :]                    # (N_oct, 1, nx)
                + s_y[:, j][:, None, None]           # (N_oct, 1, 1)
            )                                         # (N_oct, ng, nx)
            alpha = 2.0 * s_x[:, None, :] / D_row - 1.0
            beta = (
                2.0 * (Q_oct[:, :, :, j] + s_y[:, j][:, None, None] * psi_y_in)
                / D_row
            )
            psi_avg_row, out_y, x_outflow = _scanmarch_row(
                alpha, beta, inflow_x[:, :, j], psi_y_in, x_reverse,
            )
            psi_y_in = out_y
            capture_x[:, :, j] = x_outflow
            if emit.moment_buf is None:
                angular_oct[:, :, :, j] = psi_avg_row
                emit.scalar_flux[:, :, j] += np.einsum(
                    "ngi,n->gi", psi_avg_row, w_oct,
                )
            else:
                emit.moment_buf[:, :, :, :, j] += np.einsum(
                    "nlm,ngi,n->lmgi", Y_oct, psi_avg_row, w_oct,
                )

        if emit.moment_buf is None:
            emit.angular_flux[oct_idx] = angular_oct
        # x-outflow is each row's last x-scan value (captured above); the
        # y-outflow is the LAST-marched row's out_y.
        return (capture_x, out_y)

    def loss_action(
        self, operator: "StreamingOperator", psi: "TimedFullField",
    ) -> "TimedFullField":
        r"""Forward loss action ``(L+C)ψ`` — the row-march apply (L21: sweep & matvec, ONE operator).

        S6.3: the matvec ``(L+C)ψ`` walk LIVES here.  1-D → the geometry-blind
        spatial operator-sum ``operator.M_spatial`` (the ``s_y = 0`` degeneration
        of the 2-D scan-march).  2-D Cartesian → the row-march reconstruction of
        the interior faces from the KNOWN probe via
        :func:`~orpheus.sn.sweep._x_scan_faces` with the apply coefficients
        ``α = −1``, ``β = 2 ψ̄`` (a pure-reflection scan: since ψ̄ is known the WDD
        closure ``out_x = 2ψ̄ − in_x`` IS a first-order recurrence).  The per-cell
        residual is ``(Σ_t + s_x + s_y)·ψ̄ − s_x·in_x − s_y·in_y`` (``= (L+C)ψ̄`` at
        zero source); :meth:`~orpheus.sn.operator.StreamingOperator.apply`
        subtracts ``Σ_t·ψ̄`` → ``Lψ̄``.

        Principled-equivalent (NOT bit-identical) to
        :meth:`MovingFrontierWindow.loss_action`: the row-march and the
        anti-diagonal reconstruct the SAME faces in a different order, so the
        residual agrees to FP-association.  The :class:`FullFieldWavefront`
        oracle pins it (G2.c).  S6.4 sub-step (a): the octant frame + the
        O.4b boundary-residual block are the shared :class:`_OctantWalk`
        apply frame (the former Fork-B1 lockstep duplication of
        :meth:`MovingFrontierWindow.loss_action` is GONE) — this class
        supplies only the row-march interior kernel
        :meth:`_loss_action_interior`.
        """
        if self.mesh.is_1d:
            # d=1 ⇒ scan(x) with no transverse march: the spatial operator-sum
            # (the s_y = 0 degeneration of the 2-D scan-march).
            return operator.M_spatial._compute_LpC(psi)
        return _OctantWalk(self.mesh).loss_action(
            operator, psi, self._loss_action_interior,
        )

    def _loss_action_interior(
        self,
        operands: _ApplyOperands,
        oct_idx: "np.ndarray",
        signs_eff: tuple[int, ...],
        inflow: tuple["np.ndarray", ...],
    ) -> tuple["np.ndarray", tuple["np.ndarray", ...]]:
        r"""Row-march interior kernel, APPLY direction, one octant.

        Marches the y-rows in the octant's y-sweep order, reconstructing the
        x-faces from the KNOWN probe via :func:`~orpheus.sn.sweep._x_scan_faces`
        with the apply coefficients ``α = −1``, ``β = 2ψ̄`` (a pure-reflection
        scan) and threading the transverse-y face ``out_y = 2ψ̄ − ψy_in`` row
        to row.  Returns ``(LpC_octant, capture)`` with the per-axis
        domain-edge outflow ``(cap_x, out_y)`` — the x-outflow captured per
        row, the y-outflow the LAST-marched row's ``out_y``.
        """
        sig_t = operands.sig_t                      # (ng, nx, ny)
        ng, nx, ny = sig_t.shape
        sx_eff, sy_eff = signs_eff
        inflow_x, inflow_y = inflow                 # (N_oct, ng, ny) / (N_oct, ng, nx)
        s_x, s_y = (s[oct_idx] for s in operands.str_axes)  # (N_oct, nx) / (N_oct, ny)
        probe_oct = operands.probe[oct_idx]         # (N_oct, ng, nx, ny)
        N_oct = oct_idx.size

        x_reverse = sx_eff < 0
        alpha_reflect = np.full((N_oct, ng, nx), -1.0)   # α = −1, reused per row
        LpC_oct = np.empty((N_oct, ng, nx, ny))
        cap_x = np.empty((N_oct, ng, ny))            # domain x-outflow, per y-row

        # March the y-rows in the octant's y-sweep order, reconstructing the
        # transverse-y face ψ_y from the KNOWN probe (out_y = 2ψ̄ − ψy_in).
        psi_y_in = inflow_y                          # (N_oct, ng, nx) — row-0 inflow
        out_y = psi_y_in                             # last-row out_y (ny ≥ 1 → set below)
        y_rows = range(ny) if sy_eff >= 0 else range(ny - 1, -1, -1)
        for j in y_rows:
            psi_bar_row = probe_oct[:, :, :, j]      # (N_oct, ng, nx)
            # Reconstruct the x-faces from the probe: out_x = 2ψ̄ − in_x.
            in_x_row, _out_x_row, x_outflow = _x_scan_faces(
                alpha_reflect, 2.0 * psi_bar_row, inflow_x[:, :, j], x_reverse,
            )
            D_row = (
                sig_t[None, :, :, j]                 # (1, ng, nx)
                + s_x[:, None, :]                    # (N_oct, 1, nx)
                + s_y[:, j][:, None, None]           # (N_oct, 1, 1)
            )
            LpC_oct[:, :, :, j] = (
                D_row * psi_bar_row
                - s_x[:, None, :] * in_x_row
                - s_y[:, j][:, None, None] * psi_y_in
            )
            out_y = 2.0 * psi_bar_row - psi_y_in
            psi_y_in = out_y
            cap_x[:, :, j] = x_outflow
        return LpC_oct, (cap_x, out_y)

    def loss_action_transpose(
        self, operator: "StreamingOperator", phi: "TimedFullField",
    ) -> "TimedFullField":
        """Adjoint loss action — 1-D wired ``(L+C)ᵀφ``; multi-D Cartesian deferred."""
        if self.mesh.is_1d:
            return operator.M_spatial._compute_LpC_transpose(phi)
        raise NotImplementedError(
            "ScanMarch.loss_action_transpose: the multi-D Cartesian adjoint is "
            "deferred (O.2b lands the 1-D reverse sweep first; the multi-D "
            "adjoint follows the forward scan-march matvec, S5.1b+)."
        )


# ═══════════════════════════════════════════════════════════════════════
# Registry + factory — the single selection source of truth
# ═══════════════════════════════════════════════════════════════════════

#: The CONCRETE strategy leaves (never the abstract ``_LossRepresentation`` /
#: ``_DAGWavefront`` bases) — :func:`default_for` constructs whichever it
#: picks, so only buildable strategies belong here.
#:
#: Selection priority order: the 1-D scan, the 2-D production window, the
#: d-general scan-march, then the full-field oracle.  :func:`default_for`
#: returns the FIRST that applies, so the legacy production optimizations win
#: at d=1/d=2 (the scan-march is OPT-IN — issue #222 Fork B1, default
#: unchanged), the scan-march is the d≥3 Cartesian production primitive, and
#: the oracle is the never-stuck final fallback.
LOSS_REPRESENTATIONS: tuple[type[_LossRepresentation], ...] = (
    CumprodScan,
    MovingFrontierWindow,
    ScanMarch,
    FullFieldWavefront,
)


def default_for(mesh: "SNMesh") -> LossRepresentation:
    """Select the default sweep strategy for ``mesh``.

    Returns the first strategy in :data:`LOSS_REPRESENTATIONS` whose
    :meth:`~LossRepresentation.supports` admits ``mesh`` — the best *available*
    production optimization, falling back to the spine.  This reproduces the
    legacy dispatch exactly: 1-D → :class:`CumprodScan`; 2-D Cartesian →
    :class:`MovingFrontierWindow`.

    Raises
    ------
    IncompatibleRepresentation
        If no strategy applies.  Unreachable for any constructible mesh
        (every 1-D mesh → ``CumprodScan``; every 2-D Cartesian →
        ``MovingFrontierWindow``).  A hypothetical 3-D Cartesian mesh selects
        the d-general ``ScanMarch`` (production), with ``FullFieldWavefront``
        (the oracle) as the never-stuck fallback.
    """
    for cls in LOSS_REPRESENTATIONS:
        if cls.supports(mesh).ok:
            return cls(mesh)
    raise IncompatibleRepresentation(
        f"no sweep strategy supports this mesh "
        f"(ndim={mesh.ndim}, curvature={mesh.curvature!r}, "
        f"is_cartesian={mesh.is_cartesian})."
    )
