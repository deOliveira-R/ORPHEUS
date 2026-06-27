r"""Selectable representations of the S\ :sub:`N` loss operator :math:`(L+C)`.

The within-group transport solve :math:`\psi = (L+C)^{-1} q` (the *sweep*)
and its operator twin :math:`(L+C)\,\psi` (the *matvec*) admit several
distinct *algorithms*, each natural for a different mesh:

* a **1-D parallel-prefix scan** (Blelloch 1990 §1.5) — the geometry-blind
  chain recurrence (slab + sphere + cylinder), :meth:`._OneDimScanWalk.sweep`;
* a **multi-D wavefront walk** over the per-octant anti-hyperplane DAG
  (:meth:`SweepDependencyGraph.walk_full` /
  :meth:`~SweepDependencyGraph.walk_windowed`), in two buffer policies — a
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

Carve history (the sweep-strategy carve, COMPLETE 2026-06-11)
==============================================================

The S0–S6.9 arc that produced this module — full narrative + rationale on
the theory page :doc:`/theory/loss_representations` (§History):

* **S1–S4**: the protocol + thin-wrapper strategies (bit-identical rewire),
  the ``loss_action`` matvec twin (the five operator gates collapsed), the
  d-generic ``FullFieldWavefront`` oracle, the ``frontier_dim = d-1`` window.
* **S6.2–S6.5**: ``SweepStrategy → LossRepresentation`` (the abstraction IS
  the :math:`(L+C)` operator's *representation* — a matrix-free traversal);
  the walk moved INTO ``loss_action`` (returns :math:`(L+C)\psi`; the
  operator's ``−C`` is the only glue); ONE ``_OctantWalk`` frame for sweep
  AND matvec; the family-owned per-shape DAG cache; ONE representation
  instance per operator (L21 as a type fact).
* **S6.9 Fork B2 (2026-06-11)**: the multi-D Cartesian default flipped
  window → ``ScanMarch`` on measured numbers (sweep 0.57–0.84× at identical
  peak memory); the window KEPT as a selectable peer (user decision —
  multiple genuinely-different schedules are the point of selectability).

See also
========

* :doc:`/theory/loss_representations` — the capstone architecture page
  (the native lower-triangular frame, the four schedules, the selection
  SSOT, the one-walk/one-instance theorems, the Fork-B2 evidence).
* ``.claude/plans/sn_sweep_strategy.md`` — the authoritative design (the
  locked decisions, the verification strategy, phases S0–S6.9).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Protocol, runtime_checkable

import numpy as np

# S6.4(f): the orchestration (transport_sweep + the schedule loop + the 1-D
# unified body) lives IN this module — ``sweep.py`` dissolved; selector and
# bodies share one home, so the historical load-time import cycle is gone.
from orpheus.geometry import CoordSystem
from orpheus.numerics.moment_layout import (
    face_moment_count,
    is_moment_valued_by_flat_rank,
    is_moment_valued_by_rank,
)

from ..spatial._ubld import (
    AVERAGE_MOMENT,
    face_moment_tail,
    octant_moment_frame_signs,
)
from ..spatial.scheme import UpstreamState
from ..spatial.psi_half_angle_seed import CarlsonSweepContext
from ..spatial.scan import _scanmarch_row, _x_scan_faces, ordinate_scan
from ..spatial.sweep_cache import CollisionCache, GeometryCoefficients
from orpheus.transport.mesh.axis import AXIS_NAMES
from .sweep_graph import (
    OctantLabel,
    SweepDependencyGraph,
    _CellResidual,
    _CellSolve,
    _reframe,
)
from .sweep_schedule import SweepSchedule

#: Source-free ``Q_cells`` for the matvec apply (the operator action ``(L+C)ψ̄``
#: carries no volumetric source).  Shared read-only ``(1,1,1)`` zero broadcast
#: into ``residual_kernel_batch`` — the kernel only READS ``Q_cells`` (it never
#: mutates it), so one shared instance is safe and avoids a per-cell allocation.
_MATVEC_ZERO_SOURCE = np.zeros((1, 1, 1))


def frame_signs_for(
    scheme: "DiscretizationSchemeBase", signs: tuple[int, ...],
) -> "np.ndarray | None":
    r"""Sweep⇄global moment-frame sign vector bound to a scheme — or ``None``.

    The ONE binding of an octant's ``signs`` and the scheme's
    ``spatial_basis_per_axis`` to the single-source involution
    :func:`~orpheus.sn.spatial._ubld.octant_moment_frame_signs` (#240 D5b-S3 —
    the diffusion-limit root cause; the primitive owns the ``None``-for-DD/Step
    no-op convention).  Hoisted to module level so BOTH the
    :class:`_LossRepresentation` cell-op frames (via
    :meth:`_LossRepresentation._moment_frame_signs`) AND the
    :class:`_OneDimScanWalk` 1-D scan/matvec sites call ONE binding —
    ``_OneDimScanWalk`` does not inherit ``_LossRepresentation`` and so cannot
    reach the method (Pattern 2 — the last binding-dup the per-axis read was
    repeated across).
    """
    return octant_moment_frame_signs(signs, scheme.spatial_basis_per_axis)


if TYPE_CHECKING:
    from collections.abc import Callable, Mapping

    from orpheus.numerics.frame import FrameBase
    from orpheus.transport.fields.angular_flux import AngularFlux
    from orpheus.transport.fields.boundary_flux import BoundaryFlux
    from orpheus.transport.full_field import FullField
    from orpheus.transport.source_sinks import AngularSourceSink, BoundarySourceSink
    from orpheus.transport.timed_full_field import TimedFullField

    from ..mesh.augmented_mesh import SNMesh
    from ..operators.streaming import StreamingOperator
    from ..spatial.scheme import DiscretizationSchemeBase
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
        initial_guess: "FullField | None" = None,
        moment_frame: "FrameBase | None" = None,
    ) -> "tuple[np.ndarray, np.ndarray]":
        """Perform one within-group transport sweep on this strategy's mesh."""
        ...

    def loss_action(
        self, sigma: "np.ndarray", psi: "FullField",
    ) -> "FullField":
        r"""The forward within-group loss action :math:`(L+C)\,\psi` for this geometry.

        The sweep's operator-twin (L21 — sweep and matvec are different
        applications of the SAME operator): the sweep solves
        :math:`(L+C)^{-1} q`, this APPLIES :math:`(L+C)`.  **Return the FULL loss
        :math:`(L+C)\psi`, NOT :math:`L\psi`** — the operator
        (:meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply`) subtracts the
        collision diagonal :math:`C = \sigma\odot` exactly ONCE (Resolution A
        :math:`L = (L+C) - C`).  A leaf that returned :math:`L\psi` would make the
        operator subtract :math:`C` a SECOND time (a double-counted collision
        diagonal).  ``sigma`` is the ``(ng, ...)`` diagonal coefficient the matvec
        realises (#240 Phase 2 Step B — passed EXPLICITLY, symmetric with
        :meth:`sweep`'s ``sig_t``, so the composite — not the leaf — single-sources
        :math:`\sigma`); the per-geometry walk machinery is on ``self.mesh``.
        """
        ...

    def streaming_action(self, psi: "FullField") -> "FullField":
        r"""The pure σ-free streaming action :math:`L\,\psi = \Omega\cdot\nabla\psi`.

        The genuine pure-:math:`L` leaf — the spatial streaming + curvilinear
        angular redistribution discretization with NO collision diagonal.  The
        within-group WDD matvec is AFFINE in :math:`\sigma` in the forward
        direction (:math:`(L+C)\psi = \text{streaming\_action}(\psi) +
        \sigma\cdot\psi`; the curvilinear Carlson coupled-pole seed's
        :math:`\sigma`-dependence is exactly the collision diagonal it injects,
        so it cancels into :math:`\sigma\cdot\psi`).  This method names that
        σ-free primitive (``coding-elegance`` Pattern 3) and single-sources the
        ONE streaming discretization through :meth:`loss_action` at
        :math:`\sigma = 0` (Pattern 2 — the streaming walk lives ONCE in
        ``loss_action``; there is no twin σ-free discretization).
        :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply` calls this directly
        (#257 S8b) so :math:`L` reads no :math:`\sigma`: the collision diagonal
        :math:`C = M[\sigma_t]` is the separate shared multiplier leaf, and the
        composition :math:`L + C` recovers the full loss.
        """
        ...

    def loss_action_transpose(
        self, sigma: "np.ndarray", phi: "FullField",
    ) -> "FullField":
        r"""The adjoint loss action :math:`(L+C)^{\mathsf T}\,\phi` for this geometry.

        Return the FULL adjoint loss :math:`(L+C)^{\mathsf T}\phi` (the operator
        subtracts the self-adjoint diagonal :math:`C` in
        :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply_transpose`).  Raises
        :class:`NotImplementedError` for representations whose adjoint is deferred
        (the multi-D Cartesian reverse sweep — O.2b lands the 1-D reverse sweep
        first).  Never a silent wrong answer.  ``sigma`` is the ``(ng, ...)``
        diagonal coefficient, passed EXPLICITLY (#240 Phase 2 Step B).
        """
        ...

    def streaming_action_transpose(self, phi: "FullField") -> "FullField":
        r"""The pure σ-free adjoint streaming action :math:`L^{\mathsf T}\,\phi`.

        The transpose sibling of :meth:`streaming_action` (#257 S8b): the σ-free
        :math:`L^{\mathsf T}` leaf, single-sourced through
        :meth:`loss_action_transpose` at :math:`\sigma = 0`.  Used by
        :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply_transpose`.  Inherits
        the deferral contract (multi-D Cartesian raises, never a silent wrong
        answer).
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

    # ── pure-L streaming primitive (#257 S8b) ────────────────────────────
    # The σ-free streaming leaf, single-sourced through loss_action at σ = 0.
    # The within-group WDD matvec is affine in σ (the curvilinear Carlson
    # coupled-pole seed's σ-dependence is exactly the collision diagonal it
    # injects, so it cancels into σ·ψ), hence loss_action(0, ψ) == Ω·∇ψ.  The
    # streaming discretization lives ONCE in loss_action; there is no twin
    # σ-free walk (coding-elegance Pattern 2).  These concrete defaults work
    # for every strategy because each subclass overrides loss_action /
    # loss_action_transpose, and Python method resolution dispatches the
    # self.loss_action call to the right walk.

    if TYPE_CHECKING:
        # Abstract signatures for the type checker only — every concrete
        # strategy implements these (the base is never instantiated directly).
        # Declared under TYPE_CHECKING so they do NOT create a runtime method
        # that the subclass overrides would "obscure" (reportRedeclaration);
        # at runtime ``streaming_action``'s ``self.loss_action`` resolves to the
        # concrete subclass via normal MRO.
        def loss_action(
            self, sigma: "np.ndarray", psi: "FullField",
        ) -> "FullField": ...

        def loss_action_transpose(
            self, sigma: "np.ndarray", phi: "FullField",
        ) -> "FullField": ...

    def streaming_action(self, psi: "FullField") -> "FullField":
        r"""Pure σ-free forward streaming :math:`L\,\psi = \Omega\cdot\nabla\psi`.

        See the :meth:`LossRepresentation.streaming_action` protocol docstring.
        """
        return self.loss_action(self._zero_sigma_for(psi), psi)

    def streaming_action_transpose(self, phi: "FullField") -> "FullField":
        r"""Pure σ-free adjoint streaming :math:`L^{\mathsf T}\,\phi`.

        See :meth:`LossRepresentation.streaming_action_transpose`.  Inherits the
        deferral contract of :meth:`loss_action_transpose`.
        """
        return self.loss_action_transpose(self._zero_sigma_for(phi), phi)

    def _zero_sigma_for(self, field: "FullField") -> "np.ndarray":
        r"""The zero diagonal coefficient :math:`\sigma = 0` matching ``field``'s
        ``(ng, *spatial)`` shape — the σ-free streaming probe.

        The group count is read off the carrier (``field.bulk.values`` is the
        ``(N, ng, *spatial)`` angular flux); the spatial shape is the mesh's.
        """
        ng = int(field.bulk.values.shape[1])
        return np.zeros((ng, *self.mesh.spatial_shape))

    @property
    def _n_face_moments(self) -> int:
        r"""Per-face transverse moment count :math:`(\text{per\_axis})^{d-1}`.

        The interior face cochain (full or windowed) carries this many moments
        per face for the selected scheme: ``1`` for the slopeless cell-average
        closures (DD/Step → byte-identical scalar faces) and ``2^{d-1}`` for the
        bilinear UBLD Linear-Discontinuous closure (#240 D5b — d=2: 2).  Reads
        the multi-moment face-cochain width from the single-source
        :func:`~orpheus.numerics.moment_layout.face_moment_count` (shared with the
        trace producer :meth:`~orpheus.sn.mesh.augmented_mesh.SNMesh.boundary_face_layout`)."""
        return face_moment_count(self.mesh.scheme.spatial_basis_per_axis, self.mesh.ndim)

    def _moment_frame_signs(
        self, signs_eff: tuple[int, ...],
    ) -> "np.ndarray | None":
        r"""Octant sweep⇄global moment-frame sign vector — or ``None`` for DD/Step.

        The cell kernel works in the per-ordinate SWEEP frame; the iterate
        ``φ̂`` + its scattering source ``Σ_s·φ̂`` live in the GLOBAL frame, so the
        ``_CellSolve`` / ``_CellResidual`` level operations re-sign the slope
        moments at the octant boundary (#240 D5b-S3 — the diffusion-limit root
        cause).  Thin delegate to the module-level :func:`frame_signs_for`
        binding (this octant's ``signs_eff`` + the scheme's ``per_axis``), the
        ONE site that binds the involution to a scheme — shared with the
        :class:`_OneDimScanWalk` 1-D scan/matvec sites (which cannot reach this
        method).
        """
        return frame_signs_for(self.mesh.scheme, signs_eff)

    @property
    def _spatial_moment_tail(self) -> tuple[int, ...]:
        r"""Trailing per-CELL spatial-moment axis shape :math:`(\text{per\_axis})^d`.

        The bulk-field analogue of :attr:`_n_face_moments` (the FACE tail is
        ``per_axis^{d-1}``; the CELL tail is ``per_axis^d``).  A multi-moment
        closure (LD, ``per_axis > 1``) carries this trailing axis on the
        iterate / source / probe / residual / angular-octant buffers so the
        spatial-moment iterate ``φ̂`` travels between sweeps (#240 D5b-S3 — the
        unified moment matvec).  DD/Step (``per_axis == 1``) → ``()`` (no axis;
        every buffer stays byte-identical — the negative control).  Single
        source via :func:`~orpheus.numerics.moment_layout.face_moment_tail` (the same
        "append iff > 1" policy ``spatial_moment_tail`` delegates to), fed the
        per-CELL count ``per_axis^d``."""
        per_axis = self.mesh.scheme.spatial_basis_per_axis
        return face_moment_tail(per_axis ** self.mesh.ndim)

    def _inflow_to_moments(
        self, inflow: tuple["np.ndarray", ...],
    ) -> tuple["np.ndarray", ...]:
        r"""Carry the per-face domain inflow as a ``2^{d-1}``-moment object.

        A multi-moment closure's cochain consumes a ``2^{d-1}``-transverse-moment
        domain inflow per face (per-axis Legendre order ``[bar, slope, …]``); the
        boundary trace supplies it.  This is the boundary twin of Leg A's bulk
        :func:`~orpheus.sn.solver._lift_external_source_to_moments` — it
        rank-DISCRIMINATES the incoming face against the flat (moment-free) face
        rank ``d + 1`` (a scalar face ``(N_oct, ng, *transverse)`` has rank
        ``2 + (d−1) = d + 1``; transverse carries ``d−1`` axes):

        * **single-moment closure** (DD/Step, ``n == 1``): identity — the trailing
          axis is absent, every buffer byte-identical.
        * **scalar inflow** (a vacuum face / the existing scalar prescribed
          inflow): widen — zeros buffer, the AVERAGE moment (slot 0) seeded by
          the scalar, the transverse slopes ZERO (a scalar trace carries no
          along-face variation; the Leg-B asymmetry — the scalar default is
          correctly blind to the transverse slope).
        * **moment-resolved inflow** (the #251 moment-resolved boundary trace,
          ``2^{d-1}``-valued): PASS THROUGH — the producer (the widened trace,
          ``geometry.boundary_face_layout``) already carries the projected
          transverse face-slope; this method threads it unchanged.  Validates the
          trailing transverse-moment width == ``2^{d-1}`` (a clear ValueError
          otherwise — the moment-resolved relaxation must not swallow a real shape
          bug; ``coding-elegance`` Pattern 4).

        The projection that fills slot 1 lives at the call site (the MMS / the
        trace producer), NEVER here — production accepts the moment-resolved face,
        does not compute it (Pattern 6, L11 structural independence, exactly Leg A).
        """
        n = self._n_face_moments
        if n == 1:
            return inflow
        # A scalar face is (N_oct, ng, *transverse); transverse carries d−1 axes,
        # so its flat rank is 2 + (d − 1) = ndim + 1.  A moment-resolved face
        # carries one MORE (the trailing 2^{d-1}-moment) axis (#251).
        flat_face_ndim = self.mesh.ndim + 1
        widened = []
        for face in inflow:
            if is_moment_valued_by_flat_rank(face, flat_face_ndim):
                _assert_face_moment_width(
                    face, n, where=f"{type(self).__name__}._inflow_to_moments",
                )
                widened.append(face)            # thread the projected slope through
            else:
                buf = np.zeros((*face.shape, n))
                buf[..., AVERAGE_MOMENT] = face  # average moment ← scalar inflow
                widened.append(buf)
        return tuple(widened)

    def sweep(
        self,
        Q: "np.ndarray",
        sig_t: "np.ndarray",
        boundary_flux: "BoundaryFlux",
        *,
        initial_guess: "FullField | None" = None,
        moment_frame: "FrameBase | None" = None,
    ) -> "tuple[np.ndarray, np.ndarray]":
        """One within-group sweep — every concrete strategy implements it."""
        raise NotImplementedError(
            f"{type(self).__name__} must implement sweep()"
        )

    def loss_action(
        self, sigma: "np.ndarray", psi: "FullField",
    ) -> "FullField":
        """The forward loss action ``(L+C)ψ`` — every concrete leaf implements it.

        Returns the FULL within-group loss ``(L+C)ψ`` (NOT ``Lψ``); the operator
        subtracts ``C = σ⊙`` in ``apply`` (Resolution A ``L = (L+C) − C``).
        ``sigma`` is the diagonal coefficient, passed explicitly (#240 Step B).
        """
        raise NotImplementedError(
            f"{type(self).__name__} must implement loss_action()"
        )

    def loss_action_transpose(
        self, sigma: "np.ndarray", phi: "FullField",
    ) -> "FullField":
        """The adjoint loss action ``(L+C)ᵀφ`` — 1-D implemented or a deferral raise.

        Returns the FULL adjoint loss ``(L+C)ᵀφ`` (the operator subtracts ``C`` in
        ``apply_transpose``).  ``sigma`` is the diagonal coefficient (#240 Step B).
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

def _assert_face_moment_width(face: "np.ndarray", n: int, *, where: str) -> None:
    r"""Guard a moment-resolved face inflow's trailing transverse-moment width.

    A moment-resolved boundary inflow MUST carry exactly ``n = 2^{d-1}``
    transverse moments, or the relaxation silently MIS-BROADCASTS: a width-1
    trailing axis fans the single moment across all ``n`` slots — a wrong-physics
    seed numpy does NOT reject.  Single source of the face-moment-width contract,
    shared by the windowed :meth:`_LossRepresentation._inflow_to_moments` and the
    FFW-oracle :meth:`FullFieldWavefront._octant_face_cochain` seed — the twin
    face-cochain entry points must validate identically (#251).
    """
    if face.shape[-1] != n:
        raise ValueError(
            f"{where}: moment-resolved face inflow trailing width "
            f"{face.shape[-1]} != the per-face transverse-moment count "
            f"2^(d-1) = {n}.  A moment-resolved boundary inflow must carry "
            f"exactly 2^(d-1) transverse moments."
        )


def _inflow_faces(signs_eff: tuple[int, ...]) -> tuple[str, ...]:
    """Per-axis domain faces an octant's streaming ENTERS through.

    An octant streaming in the ``+a`` direction enters at the ``a``-min face
    (``("xmin", "ymin")`` for the ``(+1, +1)`` octant); a ``−a`` octant at
    the ``a``-max face.  ``signs_eff`` carries the EFFECTIVE signs (grazing
    ``0`` already mapped to ``+1`` — the streaming coefficient is zero, so
    the WDD result is sign-independent).
    """
    return tuple(
        f"{AXIS_NAMES[a]}min" if s >= 0 else f"{AXIS_NAMES[a]}max"
        for a, s in enumerate(signs_eff)
    )


def _outflow_faces(signs_eff: tuple[int, ...]) -> tuple[str, ...]:
    """Per-axis domain faces an octant's streaming EXITS through (the
    opposite of :func:`_inflow_faces`, axis by axis)."""
    return tuple(
        f"{AXIS_NAMES[a]}max" if s >= 0 else f"{AXIS_NAMES[a]}min"
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


def _moment_broadcast_sigma(
    sig: "np.ndarray", moment_valued: "np.ndarray",
) -> "np.ndarray":
    r"""``σ_t`` reshaped to broadcast over a trailing spatial-moment axis.

    The pure-z degenerate ordinates have no in-plane streaming, so the cell is
    collision-only and the loss couples to ``σ_t`` alone — the SOLVE balance is
    ``ψ̄ = Q / σ_t`` and its matvec twin is ``(L+C)ψ̄ = σ_t · ψ̄``.  At a
    multi-moment closure (LD, #240 D5b-S3) ``Q`` / ``ψ̄`` carry a trailing
    ``2^d`` spatial-moment axis that ``σ_t`` ``(ng, *spatial)`` lacks; each
    moment is scaled by the SAME scalar ``1/σ_t`` (resp. ``σ_t``), so ``σ_t``
    gains a length-1 trailing axis to broadcast.  DD/Step (no moment axis) →
    ``sig`` unchanged, byte-identical.

    The SINGLE source of the pure-z moment-broadcast convention: both the
    sweep arm (:math:`Q/σ_t`) and the matvec arm (:math:`σ_t·ψ̄`) call this, so
    the L21 twin (sweep ≡ matvec are two applications of the same collision-only
    operator) CANNOT diverge on the moment-axis reshape.  The moment-axis
    discriminator is :func:`~orpheus.numerics.moment_layout.is_moment_valued_by_rank`
    (single-sourced with the ``_CellSolve`` cell-solve source-reframe gate).
    """
    return sig[..., None] if is_moment_valued_by_rank(moment_valued, sig) else sig


@dataclass(frozen=True)
class _SweepEmit:
    r"""Solve-direction OUTPUT mode — angular field XOR harmonic moments.

    The Phase 5c output DI (which buffers are given selects the mode —
    mirroring the windowed walk's historical output contract), made
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
        for sweep in sweeps:
            oct_idx = np.asarray(sweep.indices)
            # The schedule's ``_octant_sweep`` is the SOLE in-plane projection
            # site, so the label carries exactly ``mesh.ndim`` signs — no
            # re-truncation here (a second silent projection could mask a
            # mis-sized label; a wrong length now fails loud at the face zips).
            signs = sweep.label.signs
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
        loop's inter-group reflect, :func:`~orpheus.sn.loss_representation._sweep_scheduled`);
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
            # ψ = Q/Σ_t for the in-plane-degenerate ordinates.  The trailing 2^d
            # spatial-moment axis (#240 D5b-S3) is broadcast on Σ_t by the SHARED
            # helper (the SAME convention the matvec twin's σ_t·ψ̄ rides — L21).
            q = operands.Q[oct_idx]
            emit.pure_z(oct_idx, q / _moment_broadcast_sigma(operands.sig_t, q))

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
        sigma: "np.ndarray",
        psi: "FullField",
        interior: "Callable[[_ApplyOperands, np.ndarray, tuple[int, ...], tuple[np.ndarray, ...]], tuple[np.ndarray, tuple[np.ndarray, ...]]]",
    ) -> "FullField":
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
        :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply` subtracts
        :math:`\Sigma_t\,\bar\psi` exactly once (Resolution A).
        """
        from orpheus.transport.full_field import FullField
        from orpheus.transport.source_sinks import (
            AngularSourceSink, BoundarySourceSink,
        )

        sn_mesh = self.mesh
        ndim = sn_mesh.ndim
        ng = sigma.shape[0]
        spatial = sigma.shape[1:]
        probe = psi.bulk.values
        # The matvec is intrinsically moment-valued (#240 D5b-S3): a multi-moment
        # closure (LD) carries a trailing 2^d spatial-moment axis on the probe
        # / accumulators so the apply returns the full moment residual.  DD/Step
        # (per_axis == 1) → ``()`` tail, every buffer byte-identical.  The tail
        # is read OFF the probe (its space already carries the SpatialMomentSpace
        # factor — the iterate is the single source of truth for the width).
        per_axis = sn_mesh.scheme.spatial_basis_per_axis
        moment_tail = face_moment_tail(per_axis ** ndim)
        operands = _ApplyOperands(
            probe=probe,
            sig_t=sigma,
            str_axes=tuple(sn_mesh.streaming(a) for a in range(ndim)),
            Q_zero=np.zeros((1, ng, *spatial, *moment_tail)),
        )

        # (L+C)·ψ̄ accumulator; ``apply`` subtracts Σ_t·ψ̄ → bare-streaming Lψ̄.
        LpC = np.zeros((sn_mesh.quad.N, ng, *spatial, *moment_tail))
        trace = sn_mesh.trace
        boundary = psi.boundary
        streamed = {
            face: np.zeros_like(boundary.face_view(face))
            for face in trace.face_names
        }

        def pure_z(oct_idx: "np.ndarray") -> None:
            # (L+C)·ψ̄ = σ·ψ̄ for the in-plane-degenerate ordinates.  The trailing
            # 2^d spatial-moment axis (#240 D5b-S3) is broadcast on σ by the SAME
            # shared helper the sweep twin's Q/σ_t rides (L21 — the qa Concern A
            # blocker: a quadrature WITH pure-z ordinates + a moment-valued probe
            # broadcast-crashed here without the guard the sweep already had).
            probe_oct = probe[oct_idx]
            LpC[oct_idx] = _moment_broadcast_sigma(sigma, probe_oct) * probe_oct

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

        return FullField(
            bulk=AngularSourceSink.from_mesh(LpC, sn_mesh, spatial_moments=per_axis),
            boundary=out_boundary,
        )


# ═══════════════════════════════════════════════════════════════════════
# CumprodScan — the 1-D chain prefix scan (any geometry)
# ═══════════════════════════════════════════════════════════════════════


class CumprodScan(_LossRepresentation):
    r"""1-D parallel-prefix scan — slab, sphere, cylinder via one body.

    Intrinsically 1-D: a prefix scan needs a total order (a chain).  The
    geometry difference is absorbed by the two-stratum cache, so slab +
    sphere + cylinder share THE SAME scan expression
    (:meth:`._OneDimScanWalk.sweep` → :func:`~orpheus.sn.spatial.scan.ordinate_scan`).
    The default production path for every 1-D mesh.
    """

    @classmethod
    def supports(cls, mesh: "SNMesh") -> Compatibility:
        return Compatibility(
            mesh.is_1d and mesh.scheme.is_affine_scannable,
            "requires a 1-D mesh with an affine-scannable cell-update scheme",
        )

    def sweep(
        self,
        Q: "np.ndarray",
        sig_t: "np.ndarray",
        boundary_flux: "BoundaryFlux",
        *,
        initial_guess: "FullField | None" = None,
        moment_frame: "FrameBase | None" = None,
    ) -> "tuple[np.ndarray, np.ndarray]":
        if moment_frame is not None:
            # Moment output is the 2-D windowed-SI peak-memory optimization;
            # 1-D / curvilinear meshes stay full-angular (the Morel–Montry
            # Carlson seed reads the per-ordinate iterate; lesson L21).
            raise ValueError(
                "CumprodScan.sweep: moment output (moment_frame given) "
                "is 2-D Cartesian only — 1-D/curvilinear meshes stay "
                "full-angular (the Morel–Montry Carlson seed reads the "
                "per-ordinate iterate; lesson L21)."
            )
        return _OneDimScanWalk(self.mesh).sweep(
            Q, sig_t, boundary_flux, initial_guess=initial_guess,
        )

    def loss_action(
        self, sigma: "np.ndarray", psi: "FullField",
    ) -> "FullField":
        r"""1-D forward loss action ``(L+C)ψ`` — the geometry-blind spatial sum.

        S6.3 / #206 Phase C: returns the FULL within-group loss ``(L+C)ψ``; the
        operator subtracts the collision diagonal ``C = σ⊙`` in :meth:`apply`.
        The matvec walk LIVES in :meth:`._OneDimScanWalk.loss_action` (the
        apply-direction twin of the sweep — L21 "matvec ≡ sweep"); the angular
        Morel–Montry redistribution + Carlson pole seed ride through
        ``pole_angular_closure`` there (NOT re-inlined).
        """
        return _OneDimScanWalk(self.mesh).loss_action(sigma, psi)

    def loss_action_transpose(
        self, sigma: "np.ndarray", phi: "FullField",
    ) -> "FullField":
        r"""1-D adjoint loss action ``(L+C)ᵀφ`` — the reverse spatial sum.

        S6.3 / #206 Phase C: returns ``(L+C)ᵀφ`` (the operator subtracts ``C`` in
        :meth:`apply_transpose`).  The transpose walk LIVES in
        :meth:`._OneDimScanWalk.loss_action_transpose`, which carries the
        curvilinear angular SECOND triangular factor (``closure.angular_adjoint``)
        — so the spatial reverse NEVER silently drops the angular adjoint
        (pinned by ``test_g_adjoint_reciprocity`` sphere/cyl, -O-firing).
        """
        return _OneDimScanWalk(self.mesh).loss_action_transpose(sigma, phi)


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
        self, sigma: "np.ndarray", phi: "FullField",
    ) -> "FullField":
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
    r"""Wavefront sweep — rolling :math:`(d{-}1)`-frontier buffer.

    The anti-diagonal (level-scheduled) sweep over the per-octant DAG,
    carrying only the rolling frontier of interior face fluxes (a 2-diagonal
    at d=2) — the ~30 % peak-memory win over the full-field oracle.
    Generalized to ``frontier_dim = d-1`` in S4 — the windowed WALK is
    d=3-pinned at the graph layer (``walk_windowed ≡ walk_full`` bit-id,
    ``test_sweep_graph_window_equivalence``), while ``supports`` stays
    conservatively d=2 (select narrow) until a d≥3 compute path + mesh
    exist; widen it WITH a measured d=3 profile, not before.

    A SELECTABLE PEER since the S6.9 Fork-B2 flip (2026-06-11, #222): the
    multi-D Cartesian production default is now :class:`ScanMarch` (measured
    1.2–1.8× faster at identical peak memory), and this representation is
    kept as a genuinely different schedule over the same lower-triangular
    operator (user decision: multiple proper methods ARE the point of
    selectability).  Its end-to-end coverage rides the forced-window gates in
    ``tests/sn/solve/test_scan_march_end_to_end.py`` + the explicit
    window≡full oracles.
    """

    def sweep(
        self,
        Q: "np.ndarray",
        sig_t: "np.ndarray",
        boundary_flux: "BoundaryFlux",
        *,
        initial_guess: "FullField | None" = None,
        moment_frame: "FrameBase | None" = None,
    ) -> "tuple[np.ndarray, np.ndarray]":
        return _sweep_jacobi(
            Q, sig_t, self.mesh, boundary_flux,
            moment_frame=moment_frame,
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
        :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed`
        with the ``_CellSolve`` level operation (the windowed walk of the
        solve cell kernel
        :meth:`~orpheus.sn.spatial.diamond.DiamondDifference.cell_kernel_batch`)
        over this octant's DAG, emitting per anti-hyperplane into the
        :class:`_SweepEmit` mode buffers.  Returns the per-axis domain-edge
        outflow ``capture``.
        """
        graph = self.sweep_graphs[OctantLabel(signs_eff)]
        ng = operands.sig_t.shape[0]
        spatial = operands.sig_t.shape[1:]
        # Multi-moment closures (LD's bilinear UBLD) carry a trailing
        # 2^{d-1}-moment face axis on the cochain; the domain inflow + capture
        # widen to match (#240 D5b).  Single-moment (DD/Step) → identity.
        n_face_moments = self._n_face_moments
        inflow = self._inflow_to_moments(inflow)
        capture = tuple(np.empty_like(face) for face in inflow)
        # Angular mode allocates a per-octant angular buffer (scattered into
        # the global field below); moment mode accumulates directly into the
        # shared moment tensor per anti-hyperplane, so NO per-octant angular
        # field is materialized (the Phase 5c peak-memory win).  At a
        # multi-moment closure (LD) the angular buffer carries the trailing 2^d
        # spatial-moment axis (the φ̂ iterate, #240 D5b-S3); DD/Step → ``()``.
        angular_flux_oct = (
            np.zeros((oct_idx.size, ng, *spatial, *self._spatial_moment_tail))
            if emit.moment_buf is None else None
        )
        graph.walk_windowed(
            level_op=_CellSolve(
                scheme=self.mesh.scheme,
                weights_octant=emit.weights[oct_idx],
                angular_flux_octant=angular_flux_oct,
                scalar_flux_buf=emit.scalar_flux,
                moment_buf=emit.moment_buf,
                Y_octant=None if emit.Y is None else emit.Y[oct_idx],
                moment_frame_signs=self._moment_frame_signs(signs_eff),
            ),
            inflow=inflow,
            Q_octant=operands.Q[oct_idx],
            sig_t=operands.sig_t,
            str_axes_octant=tuple(s[oct_idx] for s in operands.str_axes),
            capture=capture,
            n_face_moments=n_face_moments,
        )
        if emit.moment_buf is None:
            emit.angular_flux[oct_idx] = angular_flux_oct
        # The domain-edge capture carries the trailing 2^{d-1}-transverse-moment
        # axis at a multi-moment closure; the moment-resolved boundary trace
        # (#251 — geometry.boundary_face_layout) STORES it whole (the outflow
        # moments land in the now-moment-shaped slot, no longer collapsed to the
        # average).  DD/Step (n_face_moments == 1) → no moment axis, byte-identical.
        return capture

    def loss_action(
        self, sigma: "np.ndarray", psi: "FullField",
    ) -> "FullField":
        r"""2-D Cartesian forward loss action ``(L+C)ψ`` via the rolling-frontier window.

        S6.4 sub-step (a): routes through the shared :class:`_OctantWalk`
        apply frame (the ONE octant traversal — octant projection, pure-z
        branch, boundary I/O, the O.4b boundary residual), supplying only the
        rolling-frontier interior kernel :meth:`_loss_action_interior`
        (:meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed`
        × ``_CellResidual`` — the apply-direction walk of the SAME per-octant
        wavefront DAG and the SAME diamond-difference closure the 2-D sweep
        uses; matvec ≡ sweep, ONE discretization, L21).  Returns ``(L+C)ψ̄``;
        :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply` subtracts
        ``σ·ψ̄`` (the collision diagonal ``C``) to recover the
        bare-streaming ``Lψ̄``.
        """
        return _OctantWalk(self.mesh).loss_action(
            sigma, psi, self._loss_action_interior,
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
        :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed`
        with the ``_CellResidual`` level operation (the windowed walk of the
        apply cell kernel
        :meth:`~orpheus.sn.spatial.diamond.DiamondDifference.residual_kernel_batch`)
        over this octant's DAG.  Returns ``(LpC_octant, capture)`` — the
        octant's ``(L+C)ψ̄`` block and the per-axis domain-edge outflow.
        """
        graph = self.sweep_graphs[OctantLabel(signs_eff)]
        ng = operands.sig_t.shape[0]
        spatial = operands.sig_t.shape[1:]
        # The unified moment matvec (#240 D5b-S3): a multi-moment closure carries
        # the trailing 2^d spatial-moment axis on the residual + the 2^{d-1}-moment
        # face cochain, exactly as the SOLVE direction (_sweep_interior).  The
        # probe already carries the moment axis (operands.probe).  DD/Step
        # (per_axis == 1) → ``()`` tail + n_face_moments == 1, byte-identical.
        n_face_moments = self._n_face_moments
        LpC_oct = np.zeros((oct_idx.size, ng, *spatial, *self._spatial_moment_tail))
        inflow = self._inflow_to_moments(inflow)
        capture = tuple(np.empty_like(face) for face in inflow)
        graph.walk_windowed(
            level_op=_CellResidual(
                scheme=self.mesh.scheme,
                psi_avg_probe_octant=operands.probe[oct_idx],
                residual_octant=LpC_oct,
                moment_frame_signs=self._moment_frame_signs(signs_eff),
            ),
            inflow=inflow,
            Q_octant=operands.Q_zero,
            sig_t=operands.sig_t,
            str_axes_octant=tuple(s[oct_idx] for s in operands.str_axes),
            capture=capture,
            n_face_moments=n_face_moments,
        )
        # Domain-edge capture carries the trailing 2^{d-1}-transverse-moment axis
        # at a multi-moment closure; the moment-resolved boundary trace (#251)
        # STORES it whole into the now-moment-shaped ``streamed`` slot (which
        # inherits the widened face_view shape) — the B-residual emit below then
        # carries the outflow moments.  DD/Step → no moment axis, byte-identical.
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

    Since S6.4(d) BOTH directions route through the shared
    :class:`_OctantWalk` frame (sweep via the kernel-parameterized schedule
    loop, matvec via the apply frame) — this class supplies only the
    full-cochain interior kernels (:meth:`_sweep_interior` /
    :meth:`_loss_action_interior`, walking the d-generic
    :meth:`SweepDependencyGraph.walk_full` × the ``_CellSolve`` /
    ``_CellResidual`` level operations).
    ``supports`` is any-d Cartesian (S3) — the spine is genuinely
    dimension-generic, unlike the d=2-only window.
    """

    @classmethod
    def supports(cls, mesh: "SNMesh") -> Compatibility:
        # Override the _DAGWavefront family's d=2-only predicate: the spine is
        # the genuine d-generic oracle (it walks the per-octant DAG for any
        # Cartesian d via the d-generic ``graph.residual``).
        return Compatibility(mesh.is_cartesian, "requires Cartesian geometry")

    # ── The full-cochain boundary embedding (shared by both kernels) ──

    @staticmethod
    def _octant_face_cochain(
        spatial: tuple[int, ...],
        signs_eff: tuple[int, ...],
        inflow: tuple["np.ndarray", ...],
        n_face_moments: int = 1,
    ) -> tuple["np.ndarray", ...]:
        r"""Allocate one octant's FULL per-axis interior face cochain, with the
        domain in-edges seeded from the octant's inflow.

        Axis ``a``'s buffer carries ``n_a + 1`` face slots along its own axis
        (the fuller view — every interior face is retained, vs the window's
        rolling frontier).  Only the IN-edge slot is seeded: by the upwind
        invariant every other slot (interior + out-edge) is written before any
        read, so a zero initialization is byte-identical to the historical
        whole-trace ι_* seed.

        ``n_face_moments`` is the per-face transverse moment count
        :math:`(\text{per\_axis})^{d-1}` (DD/Step: 1; LD-2D: 2 — #240 D5b).  At
        ``> 1`` each face buffer carries a trailing moment axis; the IN-edge is
        seeded from the domain inflow.  The inflow is rank-DISCRIMINATED against
        the flat face rank ``d + 1`` (a scalar face ``(N_oct, ng, *transverse)``
        has rank ``2 + (d−1) = d + 1``), exactly as the windowed twin's
        :meth:`_inflow_to_moments`:

        * a SCALAR inflow (a vacuum face) seeds the AVERAGE moment (slot 0); the
          transverse slope moments stay zero (the scalar default is blind to the
          along-face variation — the Leg-B asymmetry);
        * a MOMENT-RESOLVED inflow (the #251 moment-resolved boundary trace,
          ``2^{d-1}``-valued) seeds ALL ``2^{d-1}`` moments — the projected
          transverse face-slope is threaded into the cochain unchanged.

        For the S2 foundation gates the domain inflow is VACUUM / zero, so the
        whole moment object is zero.  At ``n_face_moments == 1`` the trailing axis
        is ABSENT — DD's rank-r face buffers are byte-identical.
        """
        N_oct, ng = inflow[0].shape[0], inflow[0].shape[1]
        ndim = len(spatial)
        tail = face_moment_tail(n_face_moments)
        # A scalar face is (N_oct, ng, *transverse) — rank 2 + (d−1) = ndim + 1;
        # a moment-resolved face carries one MORE trailing axis (#251).
        flat_face_ndim = ndim + 1
        faces = []
        for a in range(ndim):
            face_shape = list(spatial)
            face_shape[a] += 1
            buf = np.zeros((N_oct, ng, *face_shape, *tail))
            in_edge = [slice(None)] * (2 + ndim)
            in_edge[2 + a] = 0 if signs_eff[a] >= 0 else spatial[a]
            if n_face_moments == 1:
                buf[tuple(in_edge)] = inflow[a]
            elif is_moment_valued_by_flat_rank(inflow[a], flat_face_ndim):
                # Moment-resolved inflow (#251): seed ALL 2^{d-1} moments — the
                # projected transverse face-slope threads into the cochain.  Guard
                # the width as the windowed twin does (else a width-1 trailing axis
                # silently broadcasts the single moment across all n slots).
                _assert_face_moment_width(
                    inflow[a], n_face_moments,
                    where="FullFieldWavefront._octant_face_cochain",
                )
                buf[tuple(in_edge)] = inflow[a]
            else:
                # Scalar inflow: seed the average moment (slot 0); the transverse
                # slope moments stay zero (vacuum → all zero; the scalar default
                # is blind to the along-face variation — the Leg-B asymmetry).
                buf[(*in_edge, AVERAGE_MOMENT)] = inflow[a]
            faces.append(buf)
        return tuple(faces)

    @staticmethod
    def _edge_outflow(
        psi_faces_oct: tuple["np.ndarray", ...],
        spatial: tuple[int, ...],
        signs_eff: tuple[int, ...],
    ) -> tuple["np.ndarray", ...]:
        """Extract the per-axis domain OUT-edge slots (the octant's shed
        outflow) from the walked cochain.

        At a multi-moment closure the captured slot carries the trailing
        ``2^{d-1}`` moment axis; the domain-edge consumer (the boundary trace)
        takes the average moment.  For S2 the captured outflow feeds only the
        per-axis ``capture`` arrays the walk discards on a vacuum domain edge.
        """
        ndim = len(spatial)
        capture = []
        for a in range(ndim):
            out_edge = [slice(None)] * (2 + ndim)
            out_edge[2 + a] = spatial[a] if signs_eff[a] >= 0 else 0
            capture.append(psi_faces_oct[a][tuple(out_edge)])
        return tuple(capture)

    # ── The two directions' interior kernels ──────────────────────────

    def sweep(
        self,
        Q: "np.ndarray",
        sig_t: "np.ndarray",
        boundary_flux: "BoundaryFlux",
        *,
        initial_guess: "FullField | None" = None,
        moment_frame: "FrameBase | None" = None,
    ) -> "tuple[np.ndarray, np.ndarray]":
        if moment_frame is not None:
            raise ValueError(
                "FullFieldWavefront.sweep: the full-field oracle does not "
                "implement moment output — use MovingFrontierWindow for the "
                "windowed-SI moment path."
            )
        # S6.4(d): the oracle sweep = the Jacobi schedule × the full-cochain
        # kernel on the SAME schedule loop + walk frame as production (the
        # former private ``_sweep_full_field`` frame dissolved).
        return _sweep_jacobi(
            Q, sig_t, self.mesh, boundary_flux,
            moment_frame=None,
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
        r"""Full-cochain interior kernel, SOLVE direction, one octant.

        Drives :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_full`
        with the ``_CellSolve`` level operation over the octant's
        complete per-axis face cochain — the fuller view the window replaces
        with a rolling frontier (the ``window ≡ full`` bit-identity anchor).
        Angular emit only (the oracle has no moment mode — guarded in
        :meth:`sweep`).
        """
        sig_t = operands.sig_t
        ng = sig_t.shape[0]
        spatial = sig_t.shape[1:]
        n_face_moments = self._n_face_moments
        graph = self.sweep_graphs[OctantLabel(signs_eff)]
        psi_faces_oct = self._octant_face_cochain(
            spatial, signs_eff, inflow, n_face_moments,
        )
        # The angular octant buffer carries the trailing 2^d spatial-moment axis
        # at a multi-moment closure (the φ̂ iterate, #240 D5b-S3); DD/Step →
        # ``()`` tail, byte-identical.
        angular_oct = np.zeros((oct_idx.size, ng, *spatial, *self._spatial_moment_tail))
        graph.walk_full(
            level_op=_CellSolve(
                scheme=self.mesh.scheme,
                weights_octant=emit.weights[oct_idx],
                angular_flux_octant=angular_oct,
                scalar_flux_buf=emit.scalar_flux,
                moment_frame_signs=self._moment_frame_signs(signs_eff),
            ),
            psi_faces_octant=psi_faces_oct,
            Q_octant=operands.Q[oct_idx],
            sig_t=sig_t,
            str_axes_octant=tuple(s[oct_idx] for s in operands.str_axes),
        )
        emit.angular_flux[oct_idx] = angular_oct
        capture = self._edge_outflow(psi_faces_oct, spatial, signs_eff)
        # The domain-edge capture carries the trailing 2^{d-1}-transverse-moment
        # axis at a multi-moment closure; the moment-resolved boundary trace
        # (#251) STORES it whole (oracle twin of the windowed _sweep_interior).
        # DD/Step → no moment axis, byte-identical.
        return capture

    def loss_action(
        self, sigma: "np.ndarray", psi: "FullField",
    ) -> "FullField":
        r"""Forward loss action ORACLE ``(L+C)ψ`` — the full-field DAG walk (d-generic).

        S6.4(d): routes through the shared :class:`_OctantWalk` apply frame,
        supplying the full-cochain interior kernel
        :meth:`_loss_action_interior`
        (:meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_full` ×
        ``_CellResidual`` — the full-field walk sharing the SAME cell kernel
        as the windowed walk,
        so the MATH cannot drift from
        :meth:`MovingFrontierWindow.loss_action` — only storage).  Returns
        ``(L+C)ψ̄``; :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply`
        subtracts ``σ·ψ̄``.  Sole purpose: verification (production is the
        window / the 1-D scan).
        """
        return _OctantWalk(self.mesh).loss_action(
            sigma, psi, self._loss_action_interior,
        )

    def _loss_action_interior(
        self,
        operands: _ApplyOperands,
        oct_idx: "np.ndarray",
        signs_eff: tuple[int, ...],
        inflow: tuple["np.ndarray", ...],
    ) -> tuple["np.ndarray", tuple["np.ndarray", ...]]:
        r"""Full-cochain interior kernel, APPLY direction, one octant.

        Drives :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_full`
        with the ``_CellResidual`` level operation (the full-field walk of
        the apply cell kernel
        :meth:`~orpheus.sn.spatial.diamond.DiamondDifference.residual_kernel_batch`)
        over the octant's complete face cochain.  Returns
        ``(LpC_octant, capture)``.
        """
        sig_t = operands.sig_t
        ng = sig_t.shape[0]
        spatial = sig_t.shape[1:]
        graph = self.sweep_graphs[OctantLabel(signs_eff)]
        # The unified moment matvec oracle (#240 D5b-S3): the full cochain carries
        # the 2^{d-1}-moment faces + the residual the 2^d spatial-moment axis.
        n_face_moments = self._n_face_moments
        psi_faces_oct = self._octant_face_cochain(
            spatial, signs_eff, inflow, n_face_moments,
        )
        LpC_oct = np.zeros((oct_idx.size, ng, *spatial, *self._spatial_moment_tail))
        graph.walk_full(
            level_op=_CellResidual(
                scheme=self.mesh.scheme,
                psi_avg_probe_octant=operands.probe[oct_idx],
                residual_octant=LpC_oct,
                moment_frame_signs=self._moment_frame_signs(signs_eff),
            ),
            psi_faces_octant=psi_faces_oct,
            Q_octant=operands.Q_zero,
            sig_t=sig_t,
            str_axes_octant=tuple(s[oct_idx] for s in operands.str_axes),
        )
        capture = self._edge_outflow(psi_faces_oct, spatial, signs_eff)
        # The domain-edge capture carries the trailing 2^{d-1}-transverse-moment
        # axis at a multi-moment closure; the moment-resolved boundary trace
        # (#251) STORES it whole (oracle twin of the windowed _loss_action_interior).
        # DD/Step → no moment axis, byte-identical.
        return LpC_oct, capture


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

    Selection — ``is_1d OR (is_cartesian AND ndim == 2)``: 1-D any geometry
    (the chain scan; the curvilinear Morel–Montry angular thread folds into
    the source) AND 2-D Cartesian (the row-march).  The d≥3 row-march
    (``scan(x)∘march(y, z)`` — a raster march over the transverse
    hyperplane) is the algorithm's natural generalization but the interior
    kernels unpack d=2 today, so ``supports`` tells the truth (C3.6:
    construct general, SELECT NARROW) and a d≥3 Cartesian mesh
    (constructible since C5.5/#225 via the mesh-less
    ``SNMesh.from_axes``) falls through ``default_for`` to the genuinely
    d-generic :class:`FullFieldWavefront` spine instead of misrouting
    here (pinned LIVE by ``TestD3SupportsMatrix``; the d≥3 kernel
    generalization is #227).  Widen this predicate WITH the kernel generalization,
    never before it.

    **The 2-D Cartesian PRODUCTION DEFAULT since the S6.9 Fork-B2 flip
    (2026-06-11, #222)** — the measured basis: sweep 0.57–0.84× / matvec
    0.55–0.78× the window's time at IDENTICAL peak memory (the rolling
    frontier has no memory edge over the row-march at d=2; both are ~0.7×
    the full-field oracle's peak), end-to-end fixed-source 0.82×.
    1-D still selects ``CumprodScan`` (registered first; same scan primitive,
    no march shell).  Mode-9 FP-invariance vs the window is pinned end-to-end
    by ``tests/sn/solve/test_scan_march_end_to_end.py``.
    """

    @classmethod
    def supports(cls, mesh: "SNMesh") -> Compatibility:
        # The 1-D arm reads ``is_affine_scannable`` (single-axis prefix
        # scannability — LD's 1-D scan IS valid here).  The d≥2 arm reads the
        # DISTINCT ``transverse_coupling_is_facewise`` (cross-axis
        # separability): the row-march ``scan(x) ∘ march(y)`` is exact ONLY
        # when the transverse coupling folds into the scan source as a
        # 0th-order face value (DD/Step), NOT when it is a 1st-order slope
        # moment (LD's bilinear multi-D closure).  Conflating the two — a 1-D
        # trait licensing a multi-D schedule — silently misroutes a 2-D LD
        # mesh into the inline-DD row-march (#240 D5-0); the split keeps the
        # selection honest.
        if mesh.is_1d:
            return Compatibility(
                mesh.scheme.is_affine_scannable,
                "requires an affine-scannable cell-update scheme on a "
                "1-D mesh (any geometry)",
            )
        return Compatibility(
            mesh.is_cartesian
            and mesh.ndim == 2
            and mesh.scheme.transverse_coupling_is_facewise,
            "2-D scan-march requires a scheme whose transverse coupling is "
            "facewise (separable into independent per-axis 1-D scans) — the "
            "slopeless cell-average closures (Diamond Difference, Step); "
            "Linear-Discontinuous's bilinear slope coupling needs the "
            "wavefront (the d≥3 row-march kernels are deferred — the "
            "full-field spine serves d≥3)",
        )

    def sweep(
        self,
        Q: "np.ndarray",
        sig_t: "np.ndarray",
        boundary_flux: "BoundaryFlux",
        *,
        initial_guess: "FullField | None" = None,
        moment_frame: "FrameBase | None" = None,
    ) -> "tuple[np.ndarray, np.ndarray]":
        if self.mesh.is_1d:
            # d=1 ⇒ ``scan(x)`` with no transverse march: the unified 1-D body
            # (slab + curvilinear via the two-stratum cache; the Morel–Montry
            # Carlson angular thread folds into the scan's affine source).  This
            # is the ``s_y = 0`` degeneration of the 2-D scan-march.
            if moment_frame is not None:
                raise ValueError(
                    "ScanMarch.sweep: moment output (moment_frame given) "
                    "is 2-D Cartesian only — 1-D/curvilinear meshes stay "
                    "full-angular (the Morel–Montry Carlson seed reads the "
                    "per-ordinate iterate; lesson L21)."
                )
            return _OneDimScanWalk(self.mesh).sweep(
                Q, sig_t, boundary_flux, initial_guess=initial_guess,
            )
        # 2-D ⇒ the row-march sweep = the Jacobi schedule × the scan-march
        # interior kernel on the SAME schedule loop the window uses (S6.4(b):
        # the former private ``_sweep_2d_scanmarch`` frame dissolved into the
        # shared walk — and the Gauss-Seidel schedule composes for free, the
        # inter-group reflect being kernel-agnostic).
        return _sweep_jacobi(
            Q, sig_t, self.mesh, boundary_flux,
            moment_frame=moment_frame,
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
        x-face recurrence is the first-order linear scan
        (:func:`~orpheus.sn.spatial.scan._scanmarch_row`) whose coefficients
        come from the scheme's
        :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.cartesian_scan_coefficients`
        — the row body carries NO inline diamond ``2`` or blend ``w`` (the
        scheme owns them; #240 D5a / #239 the coefficient-model lift).  The
        transverse-y coupling ``c_y · ψ_{y,in}`` rides the affine source via
        :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.source_emission`.
        Emits per row into the :class:`_SweepEmit` mode buffers.  Returns the
        per-axis domain-edge outflow ``(capture_x, out_y)``.

        The flux-independent scan coefficients are computed PER LINE (the
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
        scheme = self.mesh.scheme
        psi_y_in = inflow_y                          # (N_oct, ng, nx) — row-0 inflow
        out_y = psi_y_in                             # last-row out_y (ny ≥ 1 → set below)
        y_rows = range(ny) if sy_eff >= 0 else range(ny - 1, -1, -1)
        for j in y_rows:
            # #239 coefficient-model lift: the row-march asks the SCHEME for its
            # affine x-scan coefficients (the diamond ``2 = 1/w_DD`` + the
            # transverse coupling ``c_y = 2 g_y`` live in the scheme, NOT inline
            # here), so the scan-march is generic over every facewise closure
            # (DD today; Step rides it free once Step exists).  The transverse-y
            # ``g_y`` is the KNOWN-upstream 0th-order face value the row absorbs
            # into the affine source (``transverse_coupling_is_facewise``).
            a_scan, inverse_denom, w_row, (c_y,) = scheme.cartesian_scan_coefficients(
                s_scan=s_x[:, None, :],              # (N_oct, 1, nx) RAW g_x
                s_transverse=(s_y[:, j][:, None, None],),  # (N_oct, 1, 1) RAW g_y
                reaction_xs=sig_t[None, :, :, j],    # (1, ng, nx) Σ_t on this row
            )
            # Affine source b = source_emission(Q + c_y·ψ_y, inverse_denom, w):
            # the transverse-y direct term folds into the effective source.
            beta = scheme.source_emission(
                Q_oct[:, :, :, j] + c_y * psi_y_in, inverse_denom, w_row,
            )
            psi_avg_row, out_y, x_outflow = _scanmarch_row(
                a_scan, beta, inflow_x[:, :, j], psi_y_in, w_row, x_reverse,
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
        self, sigma: "np.ndarray", psi: "FullField",
    ) -> "FullField":
        r"""Forward loss action ``(L+C)ψ`` — the row-march apply (L21: sweep & matvec, ONE operator).

        S6.3: the matvec ``(L+C)ψ`` walk LIVES here.  1-D → the geometry-blind
        :meth:`_OneDimScanWalk._apply_walk` (the ``s_y = 0`` degeneration of the
        2-D scan-march).  2-D Cartesian → the row-march reconstruction of
        the interior faces from the KNOWN probe via
        :func:`~orpheus.sn.spatial.scan._x_scan_faces` with the apply coefficients
        ``α = −1``, ``β = 2 ψ̄`` (a pure-reflection scan: since ψ̄ is known the WDD
        closure ``out_x = 2ψ̄ − in_x`` IS a first-order recurrence).  The per-cell
        residual is ``(Σ_t + s_x + s_y)·ψ̄ − s_x·in_x − s_y·in_y`` (``= (L+C)ψ̄`` at
        zero source); :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply`
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
            # d=1 ⇒ scan(x) with no transverse march: the 1-D apply-direction
            # walk (#206 Phase C — the s_y = 0 degeneration of the 2-D
            # scan-march; the matvec walk lives in _OneDimScanWalk.loss_action).
            return _OneDimScanWalk(self.mesh).loss_action(sigma, psi)
        return _OctantWalk(self.mesh).loss_action(
            sigma, psi, self._loss_action_interior,
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
        x-faces from the KNOWN probe via :func:`~orpheus.sn.spatial.scan._x_scan_faces`
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
        scheme = self.mesh.scheme
        LpC_oct = np.empty((N_oct, ng, nx, ny))
        cap_x = np.empty((N_oct, ng, ny))            # domain x-outflow, per y-row
        s_x_row = s_x[:, None, :]                    # (N_oct, 1, nx) RAW g_x — row-invariant

        # March the y-rows in the octant's y-sweep order.  #239 coefficient-model
        # lift: the scheme owns the diamond ``2`` and the blend ``w`` — the
        # x-faces reconstruct off the probe via the scheme's reflection scan
        # (``α = −1``, ``β = 2ψ̄`` for DD), and the per-cell residual + the
        # transverse-y outflow ride the scheme's ``residual_kernel_batch`` (the
        # ÷V matvec kernel every facewise closure shares).  No inline DD.
        psi_y_in = inflow_y                          # (N_oct, ng, nx) — row-0 inflow
        out_y = psi_y_in                             # last-row out_y (ny ≥ 1 → set below)
        y_rows = range(ny) if sy_eff >= 0 else range(ny - 1, -1, -1)
        for j in y_rows:
            psi_bar_row = probe_oct[:, :, :, j]      # (N_oct, ng, nx)
            # Reconstruct the x-faces from the KNOWN probe (scheme reflection
            # scan ``ψ_out = α·ψ_in + β``; for DD ``out_x = 2ψ̄ − in_x``).
            alpha_reflect, beta_reflect = scheme.reflect_scan_coefficients(
                psi_bar_row,
            )
            in_x_row, _out_x_row, x_outflow = _x_scan_faces(
                alpha_reflect, beta_reflect, inflow_x[:, :, j], x_reverse,
            )
            # Per-cell residual (L+C)ψ̄ + the transverse-y outflow via the
            # scheme's ÷V matvec kernel: residual = (Σ_t + Σ 2g_a)ψ̄ − Σ 2g_a·in_a
            # at zero source; psi_out is the (out_x, out_y) face pair.
            residual_row, (_out_x_face, out_y) = scheme.residual_kernel_batch(
                psi_bar=psi_bar_row,
                psi_in=(in_x_row, psi_y_in),
                s_axes=(s_x_row, s_y[:, j][:, None, None]),
                reaction_xs=sig_t[:, :, j],
                Q_cells=np.zeros((1, ng, nx)),
            )
            LpC_oct[:, :, :, j] = residual_row
            psi_y_in = out_y
            cap_x[:, :, j] = x_outflow
        return LpC_oct, (cap_x, out_y)

    def loss_action_transpose(
        self, sigma: "np.ndarray", phi: "FullField",
    ) -> "FullField":
        """Adjoint loss action — 1-D wired ``(L+C)ᵀφ``; multi-D Cartesian deferred."""
        if self.mesh.is_1d:
            # #206 Phase C: the 1-D transpose walk lives in _OneDimScanWalk.
            return _OneDimScanWalk(self.mesh).loss_action_transpose(sigma, phi)
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
#: Selection priority order: the 1-D scan, the 2-D Cartesian scan-march
#: (the production default since the S6.9 Fork-B2 flip, 2026-06-11), the
#: wavefront window (a selectable peer, d=2), then the full-field oracle
#: (any-d Cartesian — the spine a d≥3 mesh falls through to, C3.6).
#: :func:`default_for` returns the FIRST that applies — the registry ORDER
#: is the default-selection policy, single-sourced here.
LOSS_REPRESENTATIONS: tuple[type[_LossRepresentation], ...] = (
    CumprodScan,
    ScanMarch,
    MovingFrontierWindow,
    FullFieldWavefront,
)


def default_for(mesh: "SNMesh") -> LossRepresentation:
    """Select the default sweep strategy for ``mesh``.

    Returns the first strategy in :data:`LOSS_REPRESENTATIONS` whose
    :meth:`~LossRepresentation.supports` admits ``mesh`` — the best *available*
    production optimization, falling back to the spine: 1-D →
    :class:`CumprodScan`; multi-D Cartesian → :class:`ScanMarch` (the S6.9
    Fork-B2 flip, 2026-06-11 — measured 0.57–0.84× the window's sweep time at
    identical peak memory, issue #222).  :class:`MovingFrontierWindow` stays a
    selectable peer (a genuinely different schedule — anti-diagonal wavefront
    vs row-march — kept by user decision: multiple proper methods ARE the
    point of selectability), pinned end-to-end by the forced-window gates in
    ``test_scan_march_end_to_end.py``.

    Raises
    ------
    IncompatibleRepresentation
        If no strategy applies.  Unreachable for any constructible mesh
        (every 1-D mesh → ``CumprodScan``; every 2-D Cartesian mesh →
        ``ScanMarch``; a d≥3 Cartesian mesh — axis-native via
        ``SNMesh.from_axes`` since C5.5 (#225) —
        → ``FullFieldWavefront``, the never-stuck any-d spine).
    """
    for cls in LOSS_REPRESENTATIONS:
        if cls.supports(mesh).ok:
            return cls(mesh)
    raise IncompatibleRepresentation(
        f"no sweep strategy supports this mesh "
        f"(ndim={mesh.ndim}, curvature={mesh.curvature!r}, "
        f"is_cartesian={mesh.is_cartesian})."
    )


# ═══════════════════════════════════════════════════════════════════════
# ORCHESTRATION — the transport-sweep entry + the schedule loop + the 1-D
# unified scan body (relocated from the dissolved ``sweep.py`` at S6.4(f):
# the walks live with their owner, and the historical lazy-import cycle is
# gone).
#
# References (carried with the 1-D body):
#
# * Hébert, A. (2009). *Applied Reactor Physics*. Ch. 3 §3.9.4
#   (pp. 141-144) — curvilinear SN cell-balance + DD difference relations.
# * Bailey, Morel & Chang (2010). NSE 165(2):149-169 — M-M clamp.
# * Lewis & Miller (1984). *Computational Methods of Neutron Transport.*
#   §4.5 (curvilinear DD); §5.3 (DD/WDD/Step/LD); §6.4 (sweep ordering).
# * Blelloch (1990). CMU-CS-90-190 §1.5 — first-order linear recurrence
#   closed form (the cumprod/cumsum scan the 1-D body evaluates).
# ═══════════════════════════════════════════════════════════════════════


def transport_sweep(
    source: "AngularSourceSink",
    sig_t: np.ndarray,
    sn_mesh: "SNMesh",
    boundary_flux: "BoundaryFlux",
    *,
    initial_guess: "FullField | None" = None,
    moment_frame: "FrameBase | None" = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Perform one full transport sweep.

    Boundary conditions are read from ``sn_mesh`` (resolved at
    construction time from the geometry mesh's :class:`BC` declarations).
    The cell-update strategy is read from ``sn_mesh.scheme``
    (defaults to :class:`DiamondDifference`).

    Single-source contract (R-1 Step 4 A1)
    --------------------------------------

    The sweep consumes ONE :class:`AngularSourceSink` carrying the
    combined per-ordinate source magnitude
    :math:`Q_n(\\vec r, g)` — whatever combination of iso (P0 + n2n +
    fission) and aniso (P_ℓ ≥ 1) the caller produced.  The producer
    side has already applied the :math:`1/W` projection (Pattern 7
    per ``coding-elegance`` SKILL.md §"Convention crosswalk template",
    lesson L18); the sweep does NOT apply ``/W`` internally to ANY
    part of the source.

    External iso scalar sources :math:`Q(\\vec r, g)` (e.g. user-
    supplied fixed-source problems) project to per-ordinate at
    construction time via
    :meth:`~orpheus.transport.source_sinks.AngularSourceSink.from_isotropic`.
    Scattering-generated sources project at the producer boundary
    via the singledispatched
    :meth:`~orpheus.transport.operators.scattering.ScatteringOperator.apply`.  Fission-
    generated sources project at the producer boundary via
    :meth:`~orpheus.transport.operators.fission.FissionOperator.apply`.

    The legacy two-parameter convention (``iso_source: ScalarSourceSink``
    + ``aniso_source: AngularSourceSink | None`` with sweep-internal
    ``/W``) is GONE.  See `#205
    <https://github.com/deOliveira-R/ORPHEUS/issues/205>`_ for the
    cross-method field architecture that will further refine the
    typed contract.

    History (chronological)
    -----------------------

    * Issue #196 PR-INDEX-5: PUBLIC contract is principled
      ``(ng, nx, ny)`` / ``(N, ng, nx, ny)``.
    * Issue #197 PR-TYPED-2: ``psi_bc: dict`` retired in favour of
      typed :class:`BoundaryFlux`.
    * Issue #197 PR-TYPED-3: typed :class:`ScalarSourceSink` /
      :class:`AngularSourceSink` inputs.
    * Issue #197 PR-TYPED-4: bare-``np.ndarray`` overload retired.
    * R-1 Step 0 (2026-05-19): curvilinear Carlson seed derives from
      ``initial_guess`` (= previous iterate; ``None`` fallback uses
      the in-iteration source angular average).
    * R-1 Step 4 A1 (2026-05-21): ``iso_source`` parameter retired;
      sweep takes one :class:`AngularSourceSink` carrying the combined
      per-ordinate source magnitude.  Producer-side ``/W`` projection
      everywhere.

    Parameters
    ----------
    source : AngularSourceSink
        Per-ordinate volumetric source, shape ``(N, ng, nx, ny)``.
        Convention: **per-ordinate magnitude** (the producer has
        already applied any required ``/W`` projection).
    sig_t : np.ndarray
        Total cross-section, shape ``(ng, nx, ny)``.
    sn_mesh : SNMesh
        :class:`SNMesh` carrying geometry, BCs, quadrature, cell-update
        strategy.
    boundary_flux : BoundaryFlux
        Persistent :class:`BoundaryFlux` (mutated in place).  Build a
        zero-initialised instance via ``BoundaryFlux.zeros_on(sn_mesh)``.
    initial_guess : FullField or None, optional
        Previous-iteration angular flux estimate, used for the
        curvilinear Carlson coupled-pole seed and the per-ordinate
        spatial-upstream inflow at the pole cell.  Ignored on slab.
        ``None`` (default) selects the in-iteration fallback seed.

        The timeless
        :class:`~orpheus.transport.full_field.FullField` composite (the
        history-bearing
        :class:`~orpheus.transport.timed_full_field.TimedFullField`
        iterate passes via inheritance).  The kernel reads only the
        per-ordinate bulk ndarray ``.bulk.values`` via the
        container-agnostic :func:`_initial_guess_values` extractor — it
        is history-blind.
    moment_frame : FrameBase or None, optional
        Phase 5c moment OUTPUT mode (2-D Cartesian ONLY — raises on a 1-D
        mesh).  ``None`` (default) returns the full per-ordinate angular flux
        (every full-angular consumer).  When given (the windowed-SI path), the
        2-D sweep accumulates the harmonic moment tensor per anti-diagonal and
        returns it INSTEAD of the angular field — the full ``(N, ng, nx, ny)``
        field is never materialized (the ~3× peak-memory win).  See
        :func:`_sweep_scheduled`.

    Returns
    -------
    bulk
        ``moment_frame is None`` → ``angular_flux`` ``(N, ng, nx, ny)``.
        Given → the harmonic moment tensor ``(L+1, 2L+1, ng, nx, ny)``.
    scalar_flux
        ``(ng, nx, ny)`` :math:`\\sum_n w_n \\psi_n` in angular mode; ``None`` in
        moment mode (the scalar IS :math:`\\phi_0^0` = ``moments[0, 0]``).

    Dispatch:

    The sweep algorithm is a first-class, selectable sweep strategy
    (``orpheus.sn.loss_representation.LossRepresentation``; ``default_for`` picks the
    default for the mesh).  This replaces the historical scattered branch —
    the ``reduced is not None`` test here and the five ``not is_1d`` gates in
    the operator algebra — with one polymorphic selection:

    * 1-D meshes → ``CumprodScan`` (:meth:`._OneDimScanWalk.sweep`; slab, sphere,
      cylinder — one body via the two-stratum cache).
    * multi-D Cartesian → ``ScanMarch`` (row-march ``scan(x)∘march(y)``; the
      S6.9 Fork-B2 default since 2026-06-11 — measured 0.57–0.84× the
      window's sweep time at identical peak memory).
      ``MovingFrontierWindow`` (anti-diagonal scheduling) stays selectable.

    The ``moment_frame`` guard (moment output is 2-D Cartesian only)
    now lives in ``CumprodScan.sweep`` — the strategy that cannot produce it
    carries its own refusal.

    (The strategy architecture's rendered API + theory page land in
    Sweep-Strategy carve phase S5; this docstring names the symbols as
    literals until then.)
    """
    Q = _unwrap_source(source)
    # S6.4(f): selector and orchestration share this module — the historical
    # sweep ↔ loss_representation lazy-import cycle is GONE.
    return default_for(sn_mesh).sweep(
        Q, sig_t, boundary_flux,
        initial_guess=initial_guess,
        moment_frame=moment_frame,
    )


def _unwrap_source(source: "AngularSourceSink") -> np.ndarray:
    """Unwrap typed :class:`AngularSourceSink` to bare ndarray.

    Issue #197 PR-TYPED-4 — strict typed input.  R-1 Step 4 A1 collapsed
    the iso / aniso parameter pair into a single per-ordinate source.
    The internal hot path consumes bare ndarray; this helper performs
    the unwrap once at the public boundary.
    """
    from orpheus.transport.source_sinks import AngularSourceSink
    if not isinstance(source, AngularSourceSink):
        raise TypeError(
            f"transport_sweep: source must be "
            f"AngularSourceSink (R-1 Step 4 A1); got "
            f"{type(source).__name__}"
        )
    return source.values


def _initial_guess_values(
    initial_guess: "FullField | None",
) -> "np.ndarray | None":
    """Extract the per-ordinate bulk ndarray from either container type.

    D-H.1c stage 4 — the sweep kernel reads only the per-ordinate bulk
    values (shape ``(N, ng, nx, ny)``) at two sites: the M-M Carlson
    coupled-pole seed (transposed slice at level p) and the pole-cell
    spatial-upstream inflow (single-cell slice at ordinate global_n).

    The :class:`~orpheus.transport.full_field.FullField` composite (and
    its :class:`~orpheus.transport.timed_full_field.TimedFullField`
    subclass) exposes the ndarray under ``.bulk.values``.  This helper
    centralises the access so the kernel stays container-agnostic and
    history-blind.

    Parameters
    ----------
    initial_guess : FullField or None
        Container carrying the previous iterate, OR ``None`` for
        cold-start.

    Returns
    -------
    np.ndarray or None
        The per-ordinate ``(N, ng, nx, ny)`` ndarray, or ``None``
        when ``initial_guess`` is ``None``.

    Notes
    -----
    Uses duck-typing on ``.bulk`` to detect the composite — avoids a
    runtime import of :class:`TimedFullField` (which would create a
    circular-dependency risk through transport↔sn).
    """
    if initial_guess is None:
        return None
    # Composite container exposes ``.bulk`` (an L2 AngularFlux); the
    # legacy bundle does not.
    bulk = getattr(initial_guess, "bulk", None)
    if bulk is not None:
        return bulk.values
    return initial_guess.values  # type: ignore[union-attr]


@dataclass(frozen=True)
class _OneDimScanWalk:
    r"""The shared 1-D-scan frame — the 1-D analogue of :class:`_OctantWalk`.

    Owns the geometry-blind 1-D SN sweep (the SOLVE direction), shared by
    :meth:`CumprodScan.sweep` and the :class:`ScanMarch` 1-D branch (#206
    Phase B — pure relocation of the former free helpers
    ``_sweep_1d_unified`` / ``_ensure_geom_cache`` / ``_ensure_coll_cache``
    / ``_run_1d_sweep`` into this frame, bit-identical). Like ``_OctantWalk``
    it is a frozen ``mesh`` holder; the per-ordinate cache stash, the slab
    joint-batch + curvilinear per-ordinate bodies, and the two-stratum cache
    ensure/stash live here. The MATVEC (apply direction) attaches in Phase C
    as the per-ordinate apply-kernel (the α=−1, β=2ψ̄ scan), mirroring
    ``_OctantWalk``'s cell-kernel injection.
    """

    mesh: "SNMesh"

    def sweep(
        self,
        Q: np.ndarray,
        sig_t: np.ndarray,
        boundary_flux: "BoundaryFlux",
        *,
        initial_guess: "FullField | None" = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""Geometry-blind 1-D SN sweep — three numpy tensor ops per ordinate.

        Replaces ``_sweep_1d_cartesian`` and ``_sweep_1d_curvilinear`` with
        one body driven by the two-stratum precomputed cache.  Slab,
        sphere, and cylinder share THE SAME scan expression; the
        cache abstracts the geometry difference (slab carries neutral
        curvature values; the M-M angular thread and Carlson seed run only
        when ``geom.level_ordinates is not None``).

        Per-ordinate hot path
        ---------------------

        1. ``b = 2 · (QV_chain + ang_contrib) · coll.inverse_denom[n]``
           — per-cell (in chain order) affine additive coefficient.
        2. ``psi_face = ordinate_scan(coll.a_attenuation[n], b, psi_in)``
           — the Blelloch closed form, three numpy ops internally.
        3. ``psi_avg = 0.5 · (psi_in_chain + psi_face)``
           — DD spatial closure.

        For the rare degenerate cylindrical pure-azimuthal ordinate
        (``geom.is_degenerate[n] == True``, ``|η| < 10^{-15}``), the scan
        is meaningless and the slow per-cell ``scheme.update`` path
        runs instead.

        Cache provenance
        ----------------

        The cache is stashed on ``self.mesh`` by :class:`SNSolver.__init__`.
        If the sweep is invoked outside the solver (e.g. ad-hoc tests),
        the cache is built lazily on first call and held on the mesh.

        Bit-identity contract
        ---------------------

        The cache-driven path produces algebraically the SAME values as the
        per-cell ``scheme.update`` reference iteration (the Pattern 2
        dual-view contract).  The cache's ``a_attenuation`` field IS the
        per-ordinate sequence of transmission coefficients that
        Step 2.5b's ``affine_coefficients`` builder produced — but
        precomputed once at solver construction rather than rebuilt every
        sweep.
        The Pattern 2 dual-view test
        (``tests/sn/spatial/test_sweep_cache.py``) pins this at
        ``rtol=1e-13`` across the parametrised geometry × ng × source
        grid.  Slab regression snapshots stay bit-identical at
        ``rtol=1e-12``.
        """
        geom = self._ensure_geom_cache()
        coll = self._ensure_coll_cache(sig_t, geom)
        return self._run(
            Q, sig_t, boundary_flux, geom, coll,
            initial_guess=initial_guess,
        )

    def loss_action(
        self, sigma: "np.ndarray", psi: "FullField",
    ) -> "FullField":
        r"""1-D forward loss action ``(L+C)ψ`` — the matvec, apply direction.

        #206 Phase C: ``(L+C)ψ`` via the shared :meth:`_apply_walk` (the
        apply-direction twin of :meth:`sweep` — L21 "matvec ≡ sweep"). Returns
        the FULL ``(L+C)ψ``;
        :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply` subtracts ``σ·ψ``
        ONCE to recover ``Lψ`` (Resolution A).  ``sigma`` is the diagonal
        coefficient, passed explicitly (#240 Step B); :meth:`_apply_walk` already
        took a plain ``sigma_t`` array, so the change is signature-only here.
        """
        from orpheus.transport.full_field import FullField
        from orpheus.transport.source_sinks import AngularSourceSink

        m_cell, m_boundary = self._apply_walk(sigma, psi)
        # The matvec output carries the trailing 2^d spatial-moment axis at a
        # multi-moment closure (the φ̂ iterate, #240 D5b-S3); the typed wrap
        # selects the SpatialMomentSpace factor.  DD/Step → no factor, byte-id.
        return FullField(
            bulk=AngularSourceSink.from_mesh(
                m_cell, self.mesh,
                spatial_moments=self.mesh.scheme.spatial_basis_per_axis,
            ),
            boundary=m_boundary,
        )

    def _apply_walk(
        self, sigma: "np.ndarray", psi: "FullField",
    ) -> "tuple[np.ndarray, BoundarySourceSink]":
        r"""The 1-D apply-direction walk — the fused ``(L+C)ψ`` single emission.

        #206 Phase C: relocated verbatim off
        ``_MSpatialOperatorSum._compute_LpC`` (single-emission matvec).
        The apply direction is the structural twin of :meth:`sweep` (L21 "matvec
        ≡ sweep"). The sweep SOLVES ``(L+C)⁻¹q``; this APPLIES ``(L+C)ψ`` to a
        KNOWN probe ψ̄.  Since the apply has a concrete ψ̄ it rides a
        scheme-specific density residual (#158 the coefficient model):

        * **Cartesian (every affine scheme)** —
          ``scheme.residual_kernel_batch`` (the ÷V ``g=|μ|/Δ`` kernel)
          returns the density residual AND the outgoing face in one call.  DD and
          LD route through it UNIFORMLY (#158/#240) — DD reproduces its diamond
          march, LD its Schur residual, with no scheme branch.  ``s = 2|μ|/Δ`` is
          scheme-agnostic (LD's kernel halves it to ``g = |μ|/Δ`` internally).
          (DD's Cartesian matvec re-associates ~1 ULP vs the pre-#240
          ``cell_balance`` path on non-power-of-2 widths — a deliberate
          principled-equivalence re-baseline, not a special case.)
        * **Curvilinear (DD-only)** — the ``cell_balance_for_streaming`` density
          path (``m_full = (denom·ψ̄ − numer_upstream)/V``) carrying the
          Morel–Montry angular redistribution thread IN-SWEEP (NOT a pure
          ``(a, inverse_denom, w)`` coefficient), with DD's diamond march
          ``out = 2ψ̄ − in`` inlined.  Curvilinear SN is DD-only (the LD
          curvilinear closure is unpublished), so this is a single-occupant
          geometry, not a polymorphism gap.  The angular-redistribution term is
          verified end-to-end by the anisotropic curvilinear MMS
          (``tests/sn/verification/mms/test_curvilinear_aniso_convergence.py``,
          ``catches("ERR-026")``); #238 retired the separately-applicable
          ``M_angular_redist`` leaf that re-walked here only to isolate this
          term (it had no production consumer).

        Returns ``(m_full_cell, m_boundary)`` — the bulk buffer in cell-first
        ``(N, ng, nx)`` layout + the boundary trace block (OUTFLOW =
        self-consistency defect, INFLOW = identity; NO BC reflection — the
        sibling ``−B`` carries it). The Morel–Montry angular redistribution +
        the Carlson coupled-pole seed (curvilinear) ride through
        ``pole_angular_closure`` (ERR-058 / #195 — NEVER re-inlined). ``sigma``
        is the ``(ng, nx)`` group diagonal cross-section, passed directly — the
        frame needs no operator handle (it may be ``σ_t`` OR the removal ``σ_r``
        — the caller single-sources it; the frame never assumes ``σ_t``).  Since #240 Phase 2 Step B the protocol
        :meth:`loss_action` / :meth:`loss_action_transpose` ALSO take ``sigma``
        directly (symmetric with :meth:`sweep`'s ``sig_t``), so the diagonal is
        single-sourced by the CALLER (``StreamingOperator.apply`` passes its
        ``sigma_t``; ``InvertibleOperator.apply`` passes the composite's
        diagonal ``self.sigma``) — the frame never reads it off an operator
        handle.
        """
        from orpheus.transport.source_sinks import BoundarySourceSink
        from ..spatial.cell_balance import cell_balance_for_streaming
        from ..spatial.pole_angular_closure import MorelMontryAngularSweep

        sn_mesh = self.mesh
        psi_view = psi.bulk.values
        quad = sn_mesh.quad
        N = quad.N
        ng = psi_view.shape[1]
        nx = sn_mesh.nx
        eps = 1e-15
        curvature_raw = getattr(sn_mesh, "curvature", None)
        curvature = curvature_raw if curvature_raw is not None else "cartesian"

        if curvature not in ("spherical", "cylindrical", "cartesian"):
            raise ValueError(f"Unknown curvature: {curvature!r}")
        if curvature == "cartesian" and not sn_mesh.is_1d:
            raise NotImplementedError(
                "_OneDimScanWalk._apply_walk: multi-D Cartesian is not "
                "handled by the 1-D scan walk; multi-D Cartesian routes "
                "through ScanMarch / _OctantWalk (the production default), "
                "not this frame."
            )

        pole_angular_closure = sn_mesh.pole_angular_closure
        if pole_angular_closure is None and curvature != "cartesian":
            pole_angular_closure = MorelMontryAngularSweep()

        mu_x = quad.mu_x
        level_indices: tuple[np.ndarray, ...] = pole_angular_closure.level_indices
        A = sn_mesh.areas

        psi_g_first = psi_view.swapaxes(0, 1)
        # The unified moment matvec (#240 D5b-S3): a multi-moment closure (LD)
        # carries a trailing 2^d spatial-moment axis on the iterate.  In 1-D
        # d=1, 2^d = per_axis (2 for LD) and the FACE cochain is 2^{d-1} = 1
        # (scalar), so only the cell probe / residual carries the axis.  DD/Step
        # (per_axis == 1) → ``()`` tail, every buffer byte-identical.  The width
        # is read OFF the iterate's space (the single source of truth).
        per_axis = sn_mesh.scheme.spatial_basis_per_axis
        moment_tail = face_moment_tail(per_axis ** sn_mesh.ndim)
        out_g_first = np.zeros((ng, N, nx, *moment_tail))

        V = sn_mesh.volumes
        sigma_gx = sigma

        boundary = psi.boundary
        trace = sn_mesh.trace
        has_inner_face = "xmin" in boundary.layout.faces
        face_outer = boundary.face_view("xmax")
        face_inner = boundary.face_view("xmin") if has_inner_face else None

        if curvature == "cartesian" and face_inner is None:
            raise ValueError(
                "Slab geometry requires psi.boundary.xmin_face to be "
                "populated."
            )

        # Wave O O.4a.2: the pole-closure inflow estimate is the GIVEN outer
        # inflow trace (``precompute_psi_state`` reads its μ<0 / inward
        # ordinates), not the reflected forward outflow.
        outer_inflow_estimate = face_outer
        psi_state = pole_angular_closure.precompute_psi_state(
            psi_view, sigma_t=sigma_gx,
            bc_outer_inflow_estimate=outer_inflow_estimate,
        )

        # The d=1 matvec probe (the iterate) is GLOBAL-frame; the residual
        # kernel works in the SWEEP frame, so for a BACKWARD sweep
        # (``direction_sign < 0``) the LD slope moment must be re-signed on the
        # probe IN and the residual OUT — the d=1 counterpart of the d≥2
        # ``_CellResidual`` frame map (#240 D5b-S3 root cause).  ``direction_sign``
        # IS the d=1 ``signs_eff``, so this rides the SAME single-source
        # binding (:func:`frame_signs_for`) as the d≥2 sites (``None`` for
        # DD/Step → byte-identical).

        def _sweep_direction(
            direction_sign: int,
            psi_face_in_init: np.ndarray,
        ) -> np.ndarray:
            frame_signs = frame_signs_for(sn_mesh.scheme, (direction_sign,))
            # The d=1 scan probe + residual are genuine moment buffers at a
            # multi-moment closure (LD); scalar at DD/Step (frame_signs None →
            # the reframe is a short-circuit no-op regardless).
            is_moment_valued = sn_mesh.scheme.is_multi_moment
            outflow_at_end = np.zeros((ng, N))
            for p, level_idx in enumerate(level_indices):
                level_idx_arr = np.asarray(level_idx)
                mu_level = mu_x[level_idx_arr]
                within_mask = (
                    mu_level > +eps if direction_sign > 0
                    else mu_level < -eps
                )
                if not np.any(within_mask):
                    continue
                global_dir = level_idx_arr[within_mask]
                abs_mu = np.abs(mu_x[global_dir])
                within_positions = np.where(within_mask)[0]

                cell_indices = list(
                    sn_mesh.dag_walk_cell_indices(
                        direction_sign=direction_sign, mu_level_idx=p,
                    )
                )
                if not cell_indices:
                    continue
                psi_face_in = psi_face_in_init[global_dir, :].T

                for i in cell_indices:
                    psi_cell = psi_g_first[:, global_dir, i]
                    if curvature == "cartesian":
                        # Coefficient model (#158/#240): the Cartesian matvec rides
                        # the scheme's group-2 ÷V kernel ``residual_kernel_batch``
                        # UNIFORMLY — DD reproduces its diamond march, LD its
                        # bilinear UBLD residual, with NO scheme branch (the kernel
                        # returns BOTH the moment residual and the outgoing face,
                        # the apply twin of the scan solve).  ``s_axes`` is the RAW
                        # down-face streaming ``g = |μ|·A_down/V = |μ|/Δ`` (slab
                        # A_down=1, V=Δ).  Source-free apply (``Q_cells = 0``).
                        # Cartesian has no Morel–Montry angular redistribution
                        # thread (the curvilinear arm below carries it).
                        # NOTE(#240): ``abs_mu / V[i]`` re-derives the raw ``g`` that
                        # ``SNMesh.streaming(0)`` already produces — a Pattern-2 dup;
                        # single-sourcing it (``streaming(0)[global_dir, i]``) is a
                        # deferred follow-up pending a widths-vs-volumes bit-id check.
                        #
                        # The unified moment matvec (#240 D5b-S3): the probe carries
                        # the trailing 2^d moment axis (``moment_tail``); the d=1
                        # FACE is scalar (2^{d-1} = 1).  ``swapaxes(0, 1)[:, :, None]``
                        # maps ``(ng, K[, 2^d]) → (K, ng, 1[, 2^d])`` (the ``[None]``
                        # inserts the n_diag axis BEFORE the moment axis) — agnostic
                        # over the trailing pack (DD scalar → byte-identical).
                        probe_cell = _reframe(
                            np.swapaxes(psi_cell, 0, 1)[:, :, None], frame_signs,
                            is_moment_valued=is_moment_valued,
                        )
                        resid, (psi_out_cell,) = (
                            sn_mesh.scheme.residual_kernel_batch(
                                psi_bar=probe_cell,
                                psi_in=(psi_face_in.T[:, :, None],),
                                s_axes=((abs_mu / V[i])[:, None, None],),
                                reaction_xs=sigma_gx[:, i][None, :, None],
                                Q_cells=_MATVEC_ZERO_SOURCE,   # source-free apply
                            )
                        )
                        resid = _reframe(
                            resid, frame_signs, is_moment_valued=is_moment_valued,
                        )
                        # resid (K, ng, 1[, 2^d]) → (ng, K[, 2^d]); the outgoing
                        # face is scalar (K, ng, 1) → (ng, K).
                        out_g_first[:, global_dir, i] = np.swapaxes(
                            resid[:, :, 0], 0, 1,
                        )
                        psi_face_in = psi_out_cell[:, :, 0].T             # (ng, K)
                        continue
                    # Curvilinear matvec — the ``cell_balance`` density path
                    # carrying the Morel–Montry angular thread (NOT a pure
                    # (a, inverse_denom, w) coefficient, so it cannot ride the
                    # coefficient-model kernel above).  Curvilinear SN is DD-only
                    # (the LD curvilinear closure is unpublished — guarded in
                    # ``LinearDiscontinuous.affine_scan_coefficients``), so DD's
                    # diamond march ``out = 2ψ̄ − in`` inlined here is a
                    # single-occupant geometry, NOT a polymorphism gap.
                    angular_denom_term, angular_numer_upstream = (
                        pole_angular_closure.cell_contribution(
                            psi_state, i, p, within_positions,
                        )
                    )
                    A_downstream = A[i + 1] if direction_sign > 0 else A[i]
                    denom, numer_upstream = cell_balance_for_streaming(
                        abs_mu=abs_mu,
                        A_downstream=A_downstream,
                        A_total=A[i] + A[i + 1],
                        total_xs=sigma_gx[:, i],
                        volume=V[i],
                        psi_face_in=psi_face_in,
                        angular_denom_term=angular_denom_term,
                        angular_numer_upstream=angular_numer_upstream,
                    )
                    m_full = (denom * psi_cell - numer_upstream) / V[i]
                    out_g_first[:, global_dir, i] = m_full
                    psi_face_in = 2.0 * psi_cell - psi_face_in   # DD diamond march
                outflow_at_end[:, global_dir] = psi_face_in
            return outflow_at_end

        # Wave O O.4a.2 — KEYSTONE DELETED. The backward sweep seeds from the
        # GIVEN outer inflow trace (``face_outer``'s μ<0 / inward ordinates),
        # NOT from the forward sweep's own reflected outflow
        # (``inflow_full = bc_outer.apply(outflow_at_boundary.T)``). This
        # decouples bulk ↔ boundary inside one matvec call: the reflective
        # coupling moves to the sibling −B, and the outer Krylov/SI loop drives
        # the inflow consistency ``ψ.inflow − B·ψ.outflow → 0``.
        outflow_at_inner = _sweep_direction(-1, face_outer)

        if curvature != "cartesian":
            # Carlson coupled-pole spatial seed (ERR-058 a, Issue #195):
            # at r = 0 the outward characteristic is the CONTINUATION of
            # the inward one — ψ(0, +μ) = ψ(0, −μ) — so the +1 sweep's
            # pole-face seed is the −1 sweep's pole-face outflow at the
            # mirror ordinate (already computed above: data, propagated
            # from the outer boundary, lower-triangular).  The historical
            # innermost-CELL-CENTRE read ψ(Δr/2) was O(h)-wrong on
            # non-flat profiles (exact on flat ψ — which is why every
            # flat-flux gate stayed green) and is retired; this is the
            # #192-deferred "inward-determines-outward" pole condition.
            pole_face_seed = outflow_at_inner.T[quad.reflection_index("x")]
        else:
            # Slab: read the GIVEN inner inflow trace (the forward sweep's
            # μ>0 seed at xmin) directly. Wave O O.4a.2 — the BC reflection
            # is NOT re-derived here; it moves to the sibling −B.
            assert face_inner is not None  # guaranteed by the cartesian guard
            pole_face_seed = face_inner
        outflow_at_boundary = _sweep_direction(+1, pole_face_seed)

        # Degenerate-ordinate branch (cylinder).
        degenerate_mask = np.abs(mu_x) < eps
        if np.any(degenerate_mask):
            global_deg = np.where(degenerate_mask)[0]
            deg_level: list[int] = []
            deg_within: list[int] = []
            for n_global in global_deg:
                for p, lvl in enumerate(level_indices):
                    lvl_arr = np.asarray(lvl)
                    pos = np.where(lvl_arr == n_global)[0]
                    if pos.size > 0:
                        deg_level.append(p)
                        deg_within.append(int(pos[0]))
                        break
            n_deg = global_deg.size
            abs_mu_deg = np.abs(mu_x[global_deg])
            zero_face = np.zeros((ng, n_deg))
            for i in range(nx):
                angular_denom_term = np.empty(n_deg)
                angular_numer_upstream = np.empty((ng, n_deg))
                for col_idx in range(n_deg):
                    denom_one, numer_one = pole_angular_closure.cell_contribution(
                        psi_state, i, deg_level[col_idx],
                        np.array([deg_within[col_idx]]),
                    )
                    angular_denom_term[col_idx] = denom_one[0]
                    angular_numer_upstream[:, col_idx] = numer_one[:, 0]

                psi_cell = psi_g_first[:, global_deg, i]
                denom, numer_upstream = cell_balance_for_streaming(
                    abs_mu=abs_mu_deg,
                    A_downstream=0.0,
                    A_total=A[i] + A[i + 1],
                    total_xs=sigma_gx[:, i],
                    volume=V[i],
                    psi_face_in=zero_face,
                    angular_denom_term=angular_denom_term,
                    angular_numer_upstream=angular_numer_upstream,
                )
                m_full = (denom * psi_cell - numer_upstream) / V[i]
                out_g_first[:, global_deg, i] = m_full

        m_cell = out_g_first.swapaxes(0, 1)

        # Wave O O.4a.2 — the boundary block of (L+C) carries the two trace
        # DIAGONALS of the block matrix; the off-diagonal −B is a sibling
        # operator (so this matvec contains NO BC reflection):
        #   * OUTFLOW slots — the self-consistency defect
        #     ``ψ.outflow − streamed`` (the r_outflow row's I·ψ.outflow
        #     diagonal minus L_out,b·ψ.bulk). UNCHANGED from pre-extraction;
        #     kept as ``computed − stored`` so the vacuum path is bit-identical
        #     (the per-row sign is free — q.outflow ≡ 0, the outflow trace is a
        #     pure definition with no source).
        #   * INFLOW slots — the identity ``ψ.inflow`` (the r_inflow row's
        #     I·ψ.inflow diagonal). NEW at O.4a.2. The sibling −B adds
        #     −B·ψ.outflow, so the full (L+C−S−F−B) inflow residual is
        #     ``ψ.inflow − B·ψ.outflow`` (the consistency the outer loop drives
        #     to q.inflow, the prescribed inflow / zero for vacuum+reflective).
        # The outflow / inflow ordinate sets are the disjoint sign(Ω·n)
        # partitions read from the unified TraceSpace selector (single source
        # of truth) — A.4 retired the inline ``mu_x > ±eps`` masks.
        m_boundary = BoundarySourceSink.zeros_on(sn_mesh)
        outer_outflow = trace.outflow_indices_for_face("xmax")
        if outer_outflow.size:
            m_boundary.face_view("xmax")[outer_outflow, :] = (
                outflow_at_boundary[:, outer_outflow].T
                - face_outer[outer_outflow, :]
            )
        outer_inflow = trace.inflow_indices_for_face("xmax")
        if outer_inflow.size:
            m_boundary.face_view("xmax")[outer_inflow, :] = (
                face_outer[outer_inflow, :]
            )
        if face_inner is not None:
            inner_outflow = trace.outflow_indices_for_face("xmin")
            if inner_outflow.size:
                m_boundary.face_view("xmin")[inner_outflow, :] = (
                    outflow_at_inner[:, inner_outflow].T
                    - face_inner[inner_outflow, :]
                )
            inner_inflow = trace.inflow_indices_for_face("xmin")
            if inner_inflow.size:
                m_boundary.face_view("xmin")[inner_inflow, :] = (
                    face_inner[inner_inflow, :]
                )

        return m_cell, m_boundary

    def loss_action_transpose(
        self, sigma: "np.ndarray", phi: "FullField",
    ) -> "FullField":
        r"""1-D adjoint loss action ``(L+C)ᵀφ`` — the matvec transpose.

        #206 Phase C: the reverse-mode adjoint of :meth:`loss_action`,
        relocated verbatim off ``_MSpatialOperatorSum._compute_LpC_transpose``
        (Wave O / O.2b, #208). The forward matvec is a forward-substitution
        sweep (lower-triangular in cell-visit order, with the Morel–Montry
        angular recurrence + Carlson pole seed forming a SECOND triangular
        factor in the ordinate index); its Euclidean transpose is the
        reverse-substitution sweep:

        * reversed cell traversal (the DD face-flux chain
          ``ψ_face_in ← 2·ψ_cell − ψ_face_in`` transposed);
        * the boundary block SWAPPED — the forward FULL operator reads the
          inflow trace and writes the outflow trace, so the transpose reads
          OUTFLOW cotangents and writes INFLOW cotangents;
        * the angular factor reversed, delegated to ``closure.angular_adjoint``
          (zero for the slab identity closure) — NEVER re-inlined; the
          Carlson coupled-pole seed adjoint routes through the mirror
          permutation.

        Every coefficient is ψ-independent (geometry + σ_t): ``denom`` is
        recomputed through the SAME ``cell_balance_for_streaming`` /
        ``cell_contribution`` the forward uses (Pattern 2 — no twin algebra).
        Returns ``(L+C)ᵀφ``;
        :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply_transpose`
        subtracts ``σ_t·φ`` ONCE (Resolution A, ``C`` a self-adjoint diagonal).
        Pinned by the G-adjoint reciprocity gate ``test_g_adjoint_reciprocity``
        (slab / sphere / cylinder, -O-firing) + its L11 wrong-trace-metric
        negative control.
        """
        from orpheus.transport.full_field import FullField
        from orpheus.transport.source_sinks import AngularSourceSink, BoundarySourceSink
        from ..spatial.cell_balance import cell_balance_for_streaming

        sn_mesh = self.mesh
        quad = sn_mesh.quad
        N = quad.N
        ng = phi.bulk.values.shape[1]
        nx = sn_mesh.nx
        eps = 1e-15
        curvature_raw = getattr(sn_mesh, "curvature", None)
        curvature = curvature_raw if curvature_raw is not None else "cartesian"
        if curvature == "cartesian" and not sn_mesh.is_1d:
            raise NotImplementedError(
                "_OneDimScanWalk.loss_action_transpose: the multi-D Cartesian "
                "adjoint is deferred (O.2b lands the 1-D reverse sweep first; "
                "the multi-D reverse sweep is a later Wave-O sub-step)."
            )

        closure = sn_mesh.pole_angular_closure
        mu_x = quad.mu_x
        # Mirror-ordinate permutation for the coupled-pole seed adjoint
        # (curvilinear only; cheap to build unconditionally).
        mirror = quad.reflection_index("x")
        A = sn_mesh.areas
        V = sn_mesh.volumes
        sgx = sigma                                  # (ng, nx)
        trace = sn_mesh.trace
        has_inner_face = "xmin" in phi.boundary.layout.faces

        out_bar = phi.bulk.values.swapaxes(0, 1)   # (ng, N, nx)
        fo = phi.boundary.face_view("xmax")                       # (N, ng)
        fi = phi.boundary.face_view("xmin") if has_inner_face else None

        psi_bar = np.zeros((ng, N, nx))
        fo_bar = np.zeros((N, ng))
        fi_bar = np.zeros((N, ng)) if has_inner_face else None
        numer_bar = [
            np.zeros((ng, np.asarray(li).size, nx))
            for li in closure.level_indices
        ]

        # ── reverse the boundary writeback (mirror _compute_LpC m_boundary) ──
        # m.outflow = (swept outflow) − ψ.outflow;  m.inflow = ψ.inflow.
        outflow_boundary_bar = np.zeros((ng, N))    # +1 sweep outflow → xmax
        outflow_inner_bar = np.zeros((ng, N))       # −1 sweep outflow → xmin (slab); pole-discarded (curv)
        oo = trace.outflow_indices_for_face("xmax")
        oi = trace.inflow_indices_for_face("xmax")
        if oo.size:
            outflow_boundary_bar[:, oo] += fo[oo].T
            fo_bar[oo] += -fo[oo]
        if oi.size:
            fo_bar[oi] += fo[oi]
        if has_inner_face:
            io = trace.outflow_indices_for_face("xmin")
            ii = trace.inflow_indices_for_face("xmin")
            if io.size:
                outflow_inner_bar[:, io] += fi[io].T
                fi_bar[io] += -fi[io]
            if ii.size:
                fi_bar[ii] += fi[ii]

        # ── ψ-independent angular_denom_term source (dummy state) ──
        psi_state_coef = closure.precompute_psi_state(
            np.zeros((N, ng, nx)),
            sigma_t=sgx,
            bc_outer_inflow_estimate=np.zeros((N, ng)),
        )

        # ── reverse the spatial DD sweeps (both directions, per level) ──
        for p, level_idx in enumerate(closure.level_indices):
            level_idx = np.asarray(level_idx)
            mu_level = mu_x[level_idx]
            for s in (+1, -1):
                within = np.where(
                    mu_level > +eps if s > 0 else mu_level < -eps
                )[0]
                if within.size == 0:
                    continue
                gd = level_idx[within]
                abs_mu = np.abs(mu_x[gd])
                cells = list(sn_mesh.dag_walk_cell_indices(
                    direction_sign=s, mu_level_idx=p,
                ))
                if not cells:
                    continue
                f_bar = (
                    outflow_boundary_bar[:, gd] if s > 0
                    else outflow_inner_bar[:, gd]
                ).copy()
                for i in reversed(cells):
                    A_downstream = A[i + 1] if s > 0 else A[i]
                    A_total = A[i] + A[i + 1]
                    angular_denom_term, _ = closure.cell_contribution(
                        psi_state_coef, i, p, within,
                    )
                    denom, _ = cell_balance_for_streaming(
                        abs_mu=abs_mu,
                        A_downstream=A_downstream,
                        A_total=A_total,
                        total_xs=sgx[:, i],
                        volume=V[i],
                        psi_face_in=np.zeros((ng, within.size)),
                        angular_denom_term=angular_denom_term,
                        angular_numer_upstream=np.zeros((ng, within.size)),
                    )                                       # (ng, n_mask)
                    ob = out_bar[:, gd, i]
                    # reverse psi_face_in = 2·psi_cell − psi_face_in_old
                    psi_bar[:, gd, i] += 2.0 * f_bar
                    f_bar = -f_bar
                    # reverse m = (denom·ψ − |μ|A_total·psi_face_in − angular_numer)/V
                    psi_bar[:, gd, i] += denom * ob / V[i]
                    f_bar += -(abs_mu * A_total)[None, :] * ob / V[i]
                    numer_bar[p][:, within, i] += -ob / V[i]
                # reverse the sweep seed
                if s > 0:
                    if curvature != "cartesian":
                        # adjoint of the Carlson coupled-pole seed: the
                        # forward +1 seed reads the −1 sweep's pole-face
                        # outflow at the mirror ordinate, so the seed
                        # cotangent routes into the −1 reversal's INITIAL
                        # outflow cotangent (mirror partners live in the
                        # same level; the s=−1 pass below reads it).
                        outflow_inner_bar[:, mirror[gd]] += f_bar
                    else:
                        fi_bar[gd] += f_bar.T               # slab +1 seed = ψ.inflow[xmin]
                else:
                    fo_bar[gd] += f_bar.T                   # −1 seed = ψ.inflow[xmax]

        # ── reverse the angular factor (delegated; zero for the slab closure) ──
        psi_ang_bar, bc_ang_bar = closure.angular_adjoint(
            tuple(numer_bar), sigma_t=sgx,
        )
        psi_bar += psi_ang_bar
        fo_bar += bc_ang_bar

        # ── assemble the typed composite ──
        m_boundary = BoundarySourceSink.zeros_on(sn_mesh)
        m_boundary.face_view("xmax")[...] = fo_bar
        if has_inner_face:
            m_boundary.face_view("xmin")[...] = fi_bar
        return FullField(
            bulk=AngularSourceSink.from_mesh(
                psi_bar.swapaxes(0, 1), sn_mesh,
            ),
            boundary=m_boundary,
        )

    def _ensure_geom_cache(self) -> GeometryCoefficients:
        """Return the geometry cache, building it on first use if absent."""
        cache = getattr(self.mesh, "_geom_cache", None)
        if cache is None:
            cache = GeometryCoefficients.from_mesh_and_quad(self.mesh)
            self.mesh._geom_cache = cache  # type: ignore[attr-defined]
        return cache

    def _ensure_coll_cache(
        self,
        sig_t: np.ndarray,
        geom: GeometryCoefficients,
    ) -> CollisionCache:
        """Return the collision cache, building it on first use if absent.

        The expected invariant (per cache-invariance test #4) is that the
        cache is constructed by :class:`SNSolver.__init__` and consumed by
        every sweep without rebuild.  Ad-hoc test callers may bypass the
        solver — in that case the cache is built lazily here.

        No bridge needed under PR-INDEX-3: ``sig_t`` arrives as principled
        ``(ng, nx, ny=1)`` and the cache consumes ``(ng, nx)`` — a single
        slice on the degenerate ``ny`` axis suffices.
        """
        cache = getattr(self.mesh, "_coll_cache", None)
        if cache is None:
            # 1-D meshes: sig_t is the principled (ng, nx) layout the cache
            # expects natively (rank-d (N, ng, *spatial); no phantom ny axis).
            sig_t_1d = sig_t  # (ng, nx)
            cache = CollisionCache.from_geometry(geom, sig_t_1d, self.mesh.scheme)
            self.mesh._coll_cache = cache  # type: ignore[attr-defined]
        return cache

    def _run(
        self,
        Q: np.ndarray,
        sig_t: np.ndarray,
        boundary_flux: "BoundaryFlux",
        geom: GeometryCoefficients,
        coll: CollisionCache,
        *,
        initial_guess: "FullField | None" = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Inner body of the unified 1-D sweep.

        Issue #196 PR-INDEX-1 through PR-INDEX-5: internal arrays AND the
        public ``transport_sweep`` signature both carry the principled
        ``(N, ng, nx, ny=1)`` layout (energy ``g`` is the *second* axis,
        NOT trailing; see :ref:`theory-sn-index-convention`).  No
        entry/exit transposes are required at the public boundary —
        caller-side principled-layout inputs flow directly through the
        sweep body.

        Issue #196 PR-INDEX-2: :class:`CollisionCache` fields carry the
        principled ``(N, ng, nx)`` layout natively; the bridge transposes at
        the cache-access sites are gone.  :class:`GeometryCoefficients` stays
        on ``(N, nx)`` / ``(N,)`` shapes — no group axis, no flip needed.

        Splits cleanly into setup (BC inflow, source pre-scale, Carlson
        seed when curvilinear) and a per-direction or per-ordinate scan:

        * **SLAB** (joint-batch): ordinates within a chain direction are
          independent (no M-M angular thread), so one
          :func:`ordinate_scan` call per chain handles the entire chain's
          ordinates at once with shape ``(nx, K, ng)`` where ``K`` is the
          number of ordinates in that direction (``N/2`` for symmetric GL).
          Exactly 2 scan calls per sweep regardless of ``N`` or ``ng``.

        * **CURVILINEAR** (sphere/cylinder, per-ordinate): the M-M angular
          thread couples ordinates sequentially within a μ-level (the
          Hébert §3.9.4 Eqs. 3.437/3.439 recurrence reads
          ``psi_angle[chain]`` updated by the *previous* ordinate in the
          level).  One ``ordinate_scan`` per ordinate per level — unchanged
          from PR-INDEX-1's pre-state.  A future parallel-prefix
          reformulation of the M-M recurrence could unlock joint-batch for
          curvilinear too (research-level; deferred per plan §7).
        """
        quad = self.mesh.quad
        N = quad.N
        nx = self.mesh.nx
        ng = Q.shape[1]                                          # (N, ng, nx, ny=1)
        weights = quad.weights
        mu = quad.mu_x

        # ── Entry layout — the public contract is the principled
        # (N, ng, *spatial) = (N, ng, nx) for 1-D (no phantom ny axis).
        Q_per_ord = Q                                            # (N, ng, nx)
        sig_t_p = sig_t                                          # (ng, nx)
        V = self.mesh.volumes                                      # (nx,) — no group axis
        scheme = self.mesh.scheme

        coord = self.mesh.reduced.coord  # type: ignore[union-attr]
        is_slab = coord is CoordSystem.CARTESIAN
        is_sphere = coord is CoordSystem.SPHERICAL

        # ── Spatial-moment width (the unified moment matvec, #240 D5b-S3 OWED-2)
        # A multi-moment closure (LD, ``per_axis > 1``) carries a trailing 2^d
        # spatial-moment axis on the iterate / source / output so the within-cell
        # slope iterate ``φ̂`` travels between sweeps (the SCAN analogue of the
        # DAG's ``_CellSolve`` φ̂ accumulation).  DD/Step (``per_axis == 1``) →
        # ``()`` tail, every buffer + every recurrence byte-identical (the
        # negative control).  Single source via ``face_moment_tail`` (the same
        # "append iff > 1" policy the rest of the carve keys on).
        per_axis = scheme.spatial_basis_per_axis
        moment_tail = face_moment_tail(per_axis ** self.mesh.ndim)
        is_moment = moment_tail != ()

        # A multi-moment closure lifts a FLAT scalar source onto the average
        # moment (slope 0), exactly as the DAG's ``_ubld_system`` does (#240
        # D5b-S3): the scattering-slope source ``Σ_s·φ̂`` arrives as a genuine
        # ``(N, ng, nx, 2^d)`` moment source, but an external / flat source (the
        # two-paths oracle, a manufactured Q̄) is rank-3 and lifts here so the
        # scan and the DAG consume the SAME source convention.  Discriminated by
        # RANK (a genuine moment source carries exactly ONE extra axis), the same
        # contract the DAG kernel uses.
        if is_moment and Q_per_ord.ndim == 3:
            lifted = np.zeros((N, ng, nx, *moment_tail))
            lifted[..., AVERAGE_MOMENT] = Q_per_ord
            Q_per_ord = lifted

        # ── Common pre-scale ──────────────────────────────────────────────
        # R-1 Step 4 A1 — single per-ordinate source.  The producer applied
        # ``1/W`` already; the sweep multiplies by cell volume V only.
        # No iso/aniso distinction internally — every WDD recurrence
        # consumes the same ``QV_per_ord``.  When the source carries the
        # spatial-moment axis (LD's ``[Q̄, Q̂]``), V broadcasts over it (each
        # moment scaled by the cell volume — the ×V source-moment convention
        # ``s_bar = Q̄·V`` / ``s_hat = Q̂·V`` the d=1 closed form consumes).
        V_b = V[None, None, :, None] if is_moment else V[None, None, :]
        QV_per_ord = Q_per_ord * V_b                             # (N, ng, nx[, 2^d])

        # Internal principled layout — angular flux (N, ng, nx[, 2^d]),
        # scalar flux (ng, nx[, 2^d]) working buffer (ny added at return).
        angular_flux = np.zeros((N, ng, nx, *moment_tail))
        scalar_flux = np.zeros((ng, nx, *moment_tail))

        # ── BC inflow + per-level Carlson seed (curvilinear only) ─────────
        #
        # Wave O (#208) O.4a.2 — BARE SWEEP: the entry ``bc_*.apply`` is GONE.
        # The reflective coupling ``ψ.inflow = B·ψ.outflow`` is no longer
        # re-applied inside the sweep; it is supplied by the CALLER as the
        # ``−B`` source term (the SI driver folds ``S + B`` into the source;
        # the direct fixed-source loops + the final reconstruction reflect the
        # persisted outflow into the inflow slots via ``SNBoundaryOperator``
        # before each sweep — see ``solver.py``). The sweep now reads the
        # SEEDED inflow trace DIRECTLY: the incoming-ordinate slots of the
        # face view ARE the inflow seed, and the outgoing-ordinate slots are
        # persisted in place after the sweep. Reading the inflow ords (before)
        # and writing the outflow ords (after) touch DISJOINT ordinate sets,
        # so aliasing the face view is safe.
        if is_slab:
            # D-H.2-C2: L2 :class:`BoundaryFlux` provides writable per-face
            # views via :meth:`face_view`.  Slab layout has both ``xmin``
            # and ``xmax`` slots (shape ``(N, ng)`` each); writes through
            # the view propagate to the flat backing buffer.  Per-cell-call
            # outflow persistence below (``xmin_face[ords] = ...``)
            # mutates these views in place.
            xmin_face = boundary_flux.face_view("xmin")   # (N, ng)
            xmax_face = boundary_flux.face_view("xmax")  # (N, ng)
            inflow_left = xmin_face    # incoming-ord slots = seeded inflow
            inflow_right = xmax_face  # incoming-ord slots = seeded inflow
            levels = [None]
            level_ordinates_list = [list(range(N))]
            bc_outer = None
            sigma_t_gx = None
            dr = None
            inflow_full = None
        else:
            # D-H.2-C2: 1-D curvilinear layout has only the outer radial
            # ``xmax`` face (the geometric pole at r=0 is a regularity
            # condition, not a BC face).  Writable view into the L2 flat
            # backing buffer.
            bc_outer = boundary_flux.face_view("xmax")  # (N, ng)
            inflow_full = bc_outer  # incoming-ord slots = seeded inflow (bare sweep)

            # Per-level Carlson coupled-pole seed delegates to the M-M
            # closure's ``psi_half_seed`` strategy — the SAME strategy the
            # matvec uses (Pattern 2 single source of truth).  Pre-amble
            # only stashes the (σ_t, Δr) bundle the per-level context
            # needs; the level loop builds the per-level CarlsonSweepContext
            # and calls ``closure.psi_half_seed(psi_level, context)``.
            # See ``MorelMontryAngularSweep.precompute_psi_state`` (the
            # matvec entry point) for the symmetric routing.
            sigma_t_gx = sig_t_p                                  # (ng, nx)
            dr = self.mesh.axis_widths[0]

            if is_sphere:
                levels = [None]
                level_ordinates_list = [list(range(N))]
            else:
                level_indices = quad.level_indices  # type: ignore[attr-defined]
                levels = list(range(len(level_indices)))
                level_ordinates_list = [list(level_indices[p]) for p in levels]

            inflow_left = inflow_right = None
            xmin_face = xmax_face = None

        # ── SLAB joint-batch fast path ────────────────────────────────────
        #
        # Slab has no M-M angular thread, no degenerate ordinates, and one
        # chain per direction.  Group ordinates by chain direction and run
        # ONE ordinate_scan per chain.  Exactly 2 scan calls per sweep.
        if is_slab:
            # Partition ordinates by direction sign (μ ≥ 0 → forward chain).
            forward_mask = mu >= 0
            forward_ords = np.where(forward_mask)[0]
            backward_ords = np.where(~forward_mask)[0]

            for direction_sign, ords in ((+1, forward_ords), (-1, backward_ords)):
                if ords.size == 0:
                    continue
                K = ords.size

                # Chain order is identical across ordinates in one direction
                # for slab — pick from the first ordinate.
                chain = geom.chain_idx[ords[0]]                   # (nx,)
                inv = geom.chain_idx_inv[ords[0]]                 # (nx,)

                # Per-ordinate inflow (cells degenerate, group axis full).
                psi_in_chain = (
                    inflow_left[ords] if direction_sign == +1
                    else inflow_right[ords]
                )                                                  # (K, ng)

                # Per-ordinate source in chain order — R-1 Step 4 A1's
                # single-source convention: ``QV_per_ord`` already encodes
                # per-ordinate magnitude × cell volume.  Slice the K
                # ordinates and reorder along the chain axis.
                QV_full_chain = QV_per_ord[ords][:, :, chain]      # (K, ng, nx[, 2^d])

                # Cache fields are (N, ng, nx) natively under PR-INDEX-2.
                # Indexed slice [ords] yields (K, ng, nx) — no transpose.
                inv_denom_chain = coll.inverse_denom[ords]         # (K, ng, nx)
                a_atten_chain = coll.a_attenuation[ords]           # (K, ng, nx)
                w_chain = coll.face_blend_weight[ords]             # (K, ng, nx)

                if not is_moment:
                    # ── Slopeless (DD/Step) flat-source scan — byte-identical ──
                    # Affine source emission b = QV·inverse_denom/w (#158
                    # coefficient model — DD's 2·QV·inv is the w=½ case).
                    # (K, ng, nx); ordinate_scan wants the cell axis leading.
                    b_chain = self.mesh.scheme.source_emission(
                        QV_full_chain, inv_denom_chain, w_chain,
                    )
                    a_scan = np.transpose(a_atten_chain, (2, 0, 1))   # (nx, K, ng)
                    b_scan = np.transpose(b_chain, (2, 0, 1))         # (nx, K, ng)
                    w_scan = np.transpose(w_chain, (2, 0, 1))         # (nx, K, ng)

                    # ONE scan call per chain — joint-batched over (K, ng).
                    psi_face_chain_scan = ordinate_scan(
                        a_scan, b_scan, psi_in_chain,
                    )                                                  # (nx, K, ng)

                    # Spatial closure ψ̄ = (1−w)ψ_in + w·ψ_out (DD's ½-mean is
                    # w=½) — face-in shifts upstream by 1.
                    psi_face_in_chain = np.empty_like(psi_face_chain_scan)
                    psi_face_in_chain[0] = psi_in_chain
                    psi_face_in_chain[1:] = psi_face_chain_scan[:-1]
                    psi_avg_scan = self.mesh.scheme.cell_average(
                        psi_face_in_chain, psi_face_chain_scan, w_scan,
                    )
                    # (nx, K, ng) → per-ordinate (ng, nx) via reorder.
                    psi_avg_per_ord = np.transpose(psi_avg_scan, (1, 2, 0))  # (K,ng,nx)
                    # Scatter back to cell-index order + write angular_flux,
                    # accumulate scalar_flux.
                    psi_avg_cell_order = psi_avg_per_ord[:, :, inv]   # (K, ng, nx)
                    angular_flux[ords, :, :] = psi_avg_cell_order
                    w_ords = weights[ords]                            # (K,)
                    scalar_flux += np.einsum(
                        "k,kgx->gx", w_ords, psi_avg_cell_order,
                    )
                    psi_face_out = psi_face_chain_scan[-1]            # (K, ng)
                else:
                    # ── Multi-moment (LD) slope-source scan ─────────────────
                    # The full LD with the threaded scattering-slope source
                    # ``Σ_s·φ̂`` (#240 D5b-S3 OWED-2): the scan propagates the
                    # downstream FACE (scalar — the d=1 face is 2^{d-1}=1) along
                    # the chain with a slope-augmented affine source, then
                    # reconstructs the (ψ̄, ψ̂) cell moments per cell — both
                    # single-sourced through the d=1 closed form (Pattern 2).
                    #
                    # FRAME (the diffusion-limit root cause, ERR-061): the cell
                    # kernel works in the per-ordinate SWEEP frame, but the
                    # iterate ``φ̂`` + the scattering source ``Σ_s·φ̂`` live in the
                    # GLOBAL frame, so for a BACKWARD sweep the slope moment must
                    # be re-signed — global→sweep IN on the source moments,
                    # sweep→global OUT on the (ψ̄, ψ̂) result — exactly the SAME
                    # single-source binding ``_CellSolve``/``_apply_walk`` ride
                    # at d≥1.  The scalar OUTGOING FACE stays sweep-frame (it
                    # propagates along the chain, never crossing into the iterate).
                    frame_signs = frame_signs_for(scheme, (direction_sign,))  # (2^d,)
                    # Inside the multi-moment (LD) branch: the source ``QV`` and
                    # the reconstructed ``(ψ̄, ψ̂)`` are genuine moment buffers
                    # (``is_moment`` ≡ ``scheme.is_multi_moment``).
                    # Source moments in chain order; global→sweep IN.
                    QV_chain_sweep = _reframe(
                        QV_full_chain, frame_signs, is_moment_valued=is_moment,
                    )
                    s_bar = QV_chain_sweep[..., AVERAGE_MOMENT]        # (K, ng, nx)
                    s_hat = QV_chain_sweep[..., 1]                     # (K, ng, nx)

                    # The per-cell d=1 closed form — the ONE LD algebra handle
                    # (slope fold shared with the matvec/per-cell Schur).
                    abs_mu_c = geom.abs_mu[ords][:, None, None]        # (K, 1, 1)
                    A_down_c = geom.A_down[ords][:, None, :]           # (K, 1, nx)
                    V_c = geom.V[ords][:, None, :]                     # (K, 1, nx)
                    # Σ_t chain-ordered (ng, nx) → broadcast (1, ng, nx) over K.
                    sig_t_chain = sig_t_p[:, chain][None, :, :]        # (1, ng, nx)
                    cf = scheme.moment_scan_closure(
                        abs_mu=abs_mu_c, A_down=A_down_c, V=V_c,
                        reaction_xs=sig_t_chain,
                    )
                    # Face-chain affine source b = flat emission + slope term.
                    b_chain = (
                        scheme.source_emission(s_bar, inv_denom_chain, w_chain)
                        + cf.scan_slope_face_source(
                            V_c, s_hat, inv_denom_chain, w_chain,
                        )
                    )                                                  # (K, ng, nx)
                    a_scan = np.transpose(a_atten_chain, (2, 0, 1))   # (nx, K, ng)
                    b_scan = np.transpose(b_chain, (2, 0, 1))         # (nx, K, ng)

                    # ONE scan call per chain — the scalar downstream FACE.
                    psi_face_chain_scan = ordinate_scan(
                        a_scan, b_scan, psi_in_chain,
                    )                                                  # (nx, K, ng)
                    psi_face_in_chain = np.empty_like(psi_face_chain_scan)
                    psi_face_in_chain[0] = psi_in_chain
                    psi_face_in_chain[1:] = psi_face_chain_scan[:-1]
                    # (nx, K, ng) → (K, ng, nx): the per-cell upstream face.
                    psi_in_cell = np.transpose(psi_face_in_chain, (1, 2, 0))

                    # Reconstruct the (ψ̄, ψ̂) cell moments — the scan twin of the
                    # per-cell Schur (sweep frame).
                    psi_bar, psi_hat = cf.scan_reconstruct(
                        V_c, s_bar, s_hat, psi_in_cell,
                    )                                                  # (K, ng, nx)
                    mom_sweep = np.stack([psi_bar, psi_hat], axis=-1)  # (K,ng,nx,2)
                    # sweep→global OUT on the moment (the involution's inverse =
                    # itself).
                    mom_global = _reframe(
                        mom_sweep, frame_signs, is_moment_valued=is_moment,
                    )
                    # Scatter back to cell-index order; write moment angular flux,
                    # accumulate the moment scalar flux φ̂ = Σ_n w_n ψ̂_n.
                    mom_cell_order = mom_global[:, :, inv, :]          # (K,ng,nx,2)
                    angular_flux[ords] = mom_cell_order
                    w_ords = weights[ords]                            # (K,)
                    scalar_flux += np.einsum(
                        "k,kgxp->gxp", w_ords, mom_cell_order,
                    )
                    psi_face_out = psi_face_chain_scan[-1]            # (K, ng) scalar

                # Persist outflow at the appropriate boundary face — the
                # last chain output is the (scalar) outgoing-face flux on that
                # side (the d=1 face cochain is moment-free for both closures).
                if direction_sign == +1:
                    xmax_face[ords] = psi_face_out             # (K, ng)
                else:
                    xmin_face[ords] = psi_face_out             # (K, ng)

        # ── CURVILINEAR per-ordinate path ─────────────────────────────────
        #
        # M-M angular thread couples ordinates sequentially within a level
        # (psi_angle[chain] is updated by ordinate m → consumed by m+1).
        # Joint-batch over ordinates is blocked; loop stays per-ordinate.
        else:
            # M-M closure owns the per-level Carlson coupled-pole seed
            # (Pattern 2 single source of truth — the matvec consumes the
            # SAME ``psi_half_seed`` strategy via
            # :meth:`MorelMontryAngularSweep.precompute_psi_state`).  The
            # strategy reads ψ at the level's ordinates (the previous
            # iterate when ``initial_guess`` is supplied; zeros on cold
            # start) and emits the cell-centred half-angle face flux
            # ``φ̄_{1/2,i,g}`` per Hébert §3.9.4 Eqs. (3.432)-(3.435).
            closure = self.mesh.pole_angular_closure
            # D-H.1c stage 4 — container-agnostic bulk extraction: the
            # FullField composite (and its TimedFullField subclass) exposes
            # the per-ordinate iterate under ``.bulk.values``; history-blind.
            ig_values = _initial_guess_values(initial_guess)
            if ig_values is not None:
                # (N, ng, nx) → (ng, N, nx)
                psi_g_first = ig_values.swapaxes(0, 1)
            else:
                psi_g_first = None

            # Carlson coupled-pole spatial seed (ERR-058, Issue #195): each
            # inward (μ<0) ordinate's pole-face outflow is captured here and
            # consumed as the spatial seed of its MIRROR outward (μ>0)
            # ordinate — the r=0 continuity ψ(0, +μ) = ψ(0, −μ).  Mirror
            # partners share a level, and the M-M thread sweeps inward
            # ordinates first, so the captured value is always data.
            mirror = self.mesh.quad.reflection_index("x")
            pole_outflow = np.zeros((mu.size, ng))

            for p_idx, level in enumerate(levels):
                ordinates_in_level = level_ordinates_list[p_idx]
                ords_arr = np.asarray(ordinates_in_level)
                mu_in_level = mu[ords_arr]
                most_inward_global = int(ords_arr[int(np.argmin(mu_in_level))])
                bc_outer_value = inflow_full[most_inward_global, :]  # type: ignore[index]
                level_weights = weights[ords_arr]
                if psi_g_first is not None:
                    psi_level = psi_g_first[:, ords_arr, :]          # (ng, M_p, nx)
                else:
                    psi_level = np.zeros((ng, ords_arr.size, nx))
                carlson_ctx = CarlsonSweepContext(
                    sigma_t=sigma_t_gx,
                    dr=dr,
                    mu_quad=mu_in_level.copy(),
                    weights=level_weights.copy(),
                    bc_outer_value=bc_outer_value,
                    mu_start=float(geom.mu_start[most_inward_global]),
                )
                phi_aux = closure.psi_half_seed(psi_level, carlson_ctx)  # (ng, nx)
                psi_angle = phi_aux.copy()                            # (ng, nx) — principled

                for m_local, global_n in enumerate(ordinates_in_level):
                    mu_n = mu[global_n]
                    w_n = weights[global_n]
                    chain = geom.chain_idx[global_n]

                    # Per-ordinate source assembly (R-1 Step 4 A1):
                    # ``QV_per_ord[global_n]`` is the per-ordinate source ×
                    # cell volume for ordinate ``global_n``, shape (ng, nx).
                    QV_full = QV_per_ord[global_n]                  # (ng, nx)
                    QV_chain = QV_full[:, chain]                    # (ng, nx)

                    # Per-ordinate spatial-upstream inflow (ng,).
                    if mu_n < 0:
                        psi_in = inflow_full[global_n]               # type: ignore[index]
                    elif geom.is_degenerate[global_n]:
                        # Degenerate (μ_r = 0) ordinate: no radial streaming —
                        # the spatial-upstream slot is inert.
                        psi_in = np.zeros(ng)
                    else:
                        # Coupled-pole seed (ERR-058 a): the mirror inward
                        # ordinate's pole-face outflow (captured below in
                        # this level's earlier M-M steps) — the r = 0
                        # continuity ψ(0, +μ) = ψ(0, −μ).  Mirrors the
                        # matvec's seed (Pattern 2 — the sweep/matvec pair
                        # stays ONE discrete system).  The historical
                        # pole-CELL-centre read of the previous iterate was
                        # O(h)-wrong on non-flat profiles (exact on flat ψ —
                        # invisible to every flat-flux gate).
                        psi_in = pole_outflow[mirror[global_n]]

                    # Degenerate cyl-axis ordinate: slow per-cell path.
                    if geom.is_degenerate[global_n]:
                        ordinate_idx = global_n if is_sphere else m_local
                        visits = list(self.mesh.dag_walk(
                            ordinate_idx=ordinate_idx,
                            mu_level_idx=level,
                        ))
                        for visit in visits:
                            i = visit.cell_idx
                            upstream = UpstreamState(
                                spatial_upstream=psi_in,
                                angular_upstream=psi_angle[:, i],
                            )
                            # scheme.update expects per-cell (ng,)
                            # arrays — sig_t / source slice on the cell axis.
                            result = scheme.update(
                                visit=visit,
                                total_xs=sig_t_p[:, i],
                                source=QV_full[:, i],
                                upstream_state=upstream,
                            )
                            psi = result.cell_average_flux           # (ng,)
                            psi_angle[:, i] = result.outgoing_angular_state
                            angular_flux[global_n, :, i] = psi
                            scalar_flux[:, i] += w_n * psi
                        continue

                    # Non-degenerate fast path: per-ordinate scan (ng, nx).
                    # psi_angle on (ng, nx); chain reorders the nx axis.
                    psi_a_in_chain = psi_angle[:, chain].copy()      # (ng, nx)
                    ang_contrib = (
                        geom.dA_w[global_n] * geom.c_in[global_n]
                    )[None, :] * psi_a_in_chain                       # (ng, nx)

                    # Cache fields are (N, ng, nx) natively under PR-INDEX-2.
                    # Indexed slice [global_n] yields (ng, nx) — no transpose.
                    inv_denom_p = coll.inverse_denom[global_n]       # (ng, nx)
                    a_atten_p = coll.a_attenuation[global_n]         # (ng, nx)
                    w_p = coll.face_blend_weight[global_n]           # (ng, nx)
                    # Affine source emission b = (QV + angular)·inverse_denom/w
                    # (#158 coefficient model — the Morel–Montry angular
                    # redistribution rides the volumetric source; DD's
                    # 2·(QV+ang)·inv is the w=½ case).
                    b = self.mesh.scheme.source_emission(
                        QV_chain + ang_contrib, inv_denom_p, w_p,
                    )  # (ng, nx)

                    # ordinate_scan: leading axis is the scan/cell axis.
                    # Pass (nx, ng) — transpose from (ng, nx).
                    psi_face_chain = ordinate_scan(
                        a_atten_p.T, b.T, psi_in,
                    )                                                 # (nx, ng)
                    if mu_n < 0:
                        # Coupled-pole capture: an inward chain ends at the
                        # pole; its final face value seeds the mirror outward
                        # ordinate's chain (consumed above).
                        pole_outflow[global_n] = psi_face_chain[-1]

                    # Spatial closure ψ̄ = (1−w)ψ_in + w·ψ_out (w_p.T → (nx, ng)
                    # to match the scan layout) — vectorised cell-average.
                    psi_face_in_chain = np.empty_like(psi_face_chain)
                    psi_face_in_chain[0] = psi_in
                    psi_face_in_chain[1:] = psi_face_chain[:-1]
                    psi_avg_chain = self.mesh.scheme.cell_average(
                        psi_face_in_chain, psi_face_chain, w_p.T,
                    )
                    # Principled view: (ng, nx).
                    psi_avg_chain_p = psi_avg_chain.T                # (ng, nx)

                    # M-M angular thread (curvilinear-only).
                    psi_angle_out_chain_p = (
                        geom.tau_inv[global_n] * psi_avg_chain_p
                        - geom.mm_a_in_coeff[global_n] * psi_a_in_chain
                    )                                                 # (ng, nx)
                    psi_angle[:, chain] = psi_angle_out_chain_p

                    # Scatter back to cell-index order + writes.
                    inv = geom.chain_idx_inv[global_n]
                    psi_avg_p = psi_avg_chain_p[:, inv]              # (ng, nx)
                    angular_flux[global_n, :, :] = psi_avg_p
                    scalar_flux += w_n * psi_avg_p

                    # Persist outflow at the outer face for outward ordinates.
                    if mu_n >= 0 and abs(mu_n) >= self.mesh._DEGENERATE_ABS_ETA_THRESHOLD:
                        bc_outer[global_n] = psi_face_chain[-1]      # (ng,)

        # ── Exit — PR-INDEX-5: caller consumes principled layout ──────────
        # R-1 Step 0: NO iteration-cache write-back.  The caller threads the
        # returned ``angular_flux`` as ``initial_guess`` to the NEXT sweep —
        # that's all the "history" needed.  The matvec already operates this
        # way; the sweep now mirrors it.
        # angular_flux is (N, ng, nx); scalar_flux is (ng, nx) — the principled
        # (N, ng, *spatial) / (ng, *spatial) public contract (no phantom ny).
        return angular_flux, scalar_flux


def _sweep_scheduled(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: "SNMesh",
    boundary_flux: "BoundaryFlux",
    *,
    schedule: "SweepSchedule",
    reflect: "Callable[[BoundaryFlux, tuple[str, ...]], None] | None" = None,
    moment_frame: "FrameBase | None" = None,
    interior: "Callable" ,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Polymorphic schedule-driven 2-D sweep (Phase 3 sub-step 3c; S6.4(b)
    kernel-parameterized).

    ONE uniform sweep-and-reflect loop parameterized by ``schedule`` (the
    Jacobi / Gauss-Seidel splitting — there is NO ``if jacobi/gs`` branch;
    the splitting IS the schedule) **and by the representation's solve
    ``interior`` kernel** (S6.4(b) — the per-group octant frame is the shared
    :meth:`~orpheus.sn.loss_representation._OctantWalk.sweep_group`; this
    loop does not know HOW an octant's interior is traversed).  Because the
    inter-group reflect is kernel-agnostic, any representation's kernel
    composes with any schedule — e.g. the scan-march gains Gauss-Seidel for
    free.

    1. The GIVEN inflow ``boundary_flux`` (carrying ``B·ψₙ`` — the lagged
       reflection of the previous iterate, prepared by the caller) is read
       per-octant by the walk; there is no separate whole-trace seed.
    2. ``for group in schedule.groups``: walk the group's octants
       (:meth:`_OctantWalk.sweep_group` × ``interior``), sheding each octant's
       outflow into ``boundary_flux`` (the ι* absorb, per-octant during the
       walk). Then if ``reflect`` is given AND the group has reflective
       outgoing faces, apply ``reflect`` (the ``−B`` outflow→inflow
       reflection, in place, face-restricted) so a LATER group reads the
       fresh current-iterate inflow directly off ``boundary_flux`` (the
       ``(L+C−B_lower)⁻¹`` forward substitution) — no re-seed needed, the
       next walk reads the trace fresh.

    * **Jacobi** (``reflect=None``, one all-octants group) — every octant reads
      the frozen seed; the inter-group reflect never fires. This is exactly the
      bare sweep :func:`_sweep_jacobi` passes.
    * **Gauss-Seidel** (``reflect`` = the face-restricted ``−B``, one group per
      in-plane octant) — later groups see earlier groups' fresh reflected
      outflow. The SI scheduled resolvent supplies both: its ``.solve`` seeds
      ``B·ψₙ`` onto ``boundary_flux`` then calls this with the G-S schedule +
      the reflect closure. The walk's per-octant shed populates the fresh
      outflow that ``reflect`` then maps to inflow (replacing the storage-A
      ``absorb``-before-``reflect`` step).

    The converged fixed point is INVARIANT under ``schedule`` (any consistent
    splitting of ``(L+C−S−B)ψ=q`` shares ψ\*); only the SI spectral rate
    changes. NOTE (Phase 3 spike, issue #2/#215): this folds the BOUNDARY
    coupling ``B`` only — a modest reflective-SI rate gain. The dominant
    within-group SCATTERING ``c``-mode is NOT folded here (it cannot be folded
    into a directional sweep); that is consistent DSA / Krylov territory.

    Phase 5b storage-B: the interior face cochain is a rolling 2-diagonal
    moving-frontier window (O(N·ng·nx)), carried inside
    :meth:`SweepDependencyGraph.walk_windowed` per octant — NOT the full
    O(N·ng·nx·ny) per-axis field. The full-field walk
    (:meth:`SweepDependencyGraph.walk_full` with the solve level operation,
    on the full per-axis face cochain —
    the ``FullFieldWavefront`` kernel since S6.4(d)) is retained as the
    bit-identity verification oracle (see the ``window ≡ full-field``
    test); the converged solution is unchanged.

    Phase 5c moment output: when ``moment_frame`` is given (the 2-D
    Cartesian windowed-SI path), the walk accumulates the harmonic moment tensor
    ``(L+1, 2L+1, ng, nx, ny)`` per anti-diagonal directly — the full
    per-ordinate angular OUTPUT ``(N, ng, nx, ny)`` is never materialized (the
    ~3× linear peak-memory win; the persistent SI iterate is already moments,
    5a).  Returns ``(moment_buf, None)`` — the scalar flux is :math:`\phi_0^0`
    = ``moment_buf[0, 0]`` (``Y_0^0 = 1``), read off the tensor, NOT returned
    separately (the angular-mode scalar is an independent array; ``None`` keeps
    the modes' second slot from being mistaken).  Principled-equivalence, NOT
    bit-identity: the cross-octant ``+=`` reorders the ordinate sum vs the
    post-sweep flat :attr:`~orpheus.numerics.frame.FrameBase.analysis` projection
    reduce (≤ 4 ULP de-risk).  ``moment_frame is None`` (every full-angular
    consumer — reconstruction, Krylov, 1-D) returns ``(angular_flux,
    scalar_flux)`` exactly as before.
    """
    # d-generic buffer setup (S6.4(d): the full-field SPINE routes through
    # this loop too, so the orchestrator must admit d = 1 as well as d = 2;
    # the historical ``ng, nx, ny`` unpack was a 2-D hardcode).
    ng = sig_t.shape[0]
    spatial = sig_t.shape[1:]
    N = sn_mesh.quad.N
    weights = sn_mesh.quad.weights
    # The output buffers carry the trailing 2^d spatial-moment axis at a
    # multi-moment closure (the φ̂ iterate accumulated by ``_CellSolve``; #240
    # D5b-S3) — both the angular field (FFW oracle / non-windowed) and the
    # harmonic-moment tensor (windowed production).  DD/Step (per_axis == 1) →
    # ``()`` tail, every buffer byte-identical (the negative control).
    moment_tail = face_moment_tail(
        sn_mesh.scheme.spatial_basis_per_axis ** sn_mesh.ndim
    )
    if moment_frame is None:
        angular_flux = np.zeros((N, ng, *spatial, *moment_tail))
        scalar_flux = np.zeros((ng, *spatial, *moment_tail))
        moment_buf = None
        emit = _SweepEmit(
            weights=weights, angular_flux=angular_flux, scalar_flux=scalar_flux,
        )
    else:
        # The moment buffer carries the frame table's (L+1, 2L+1) harmonic
        # block (table shape (N, L+1, 2L+1)); size it from there so the buffer
        # and the projection table agree by construction — and the sweep stays
        # basis-agnostic (it reads the order off the table it already consumes,
        # not a basis-specific attribute).
        n_l, n_m = moment_frame.table.shape[1:3]
        moment_buf = np.zeros((n_l, n_m, ng, *spatial, *moment_tail))
        emit = _SweepEmit(
            weights=weights, moment_buf=moment_buf, Y=moment_frame.table,
        )

    operands = _SolveOperands(
        Q=Q, sig_t=sig_t,
        str_axes=tuple(sn_mesh.streaming(a) for a in range(sn_mesh.ndim)),
    )
    walk = _OctantWalk(sn_mesh)
    for group in schedule.groups:
        walk.sweep_group(
            group,
            operands=operands,
            emit=emit,
            boundary_flux=boundary_flux,
            interior=interior,
        )
        if reflect is not None and group.reflect_faces:
            # G-S inter-group reflect (a no-op for the Jacobi schedule, whose
            # sole group carries no reflect_faces): the walk already shed this
            # group's fresh outflow into ``boundary_flux``; reflect ``−B``
            # (outflow→inflow, in place, face-restricted) so the NEXT group
            # reads the fresh current-iterate reflected inflow off the trace.
            reflect(boundary_flux, group.reflect_faces)

    if moment_frame is None:
        return angular_flux, scalar_flux
    # Moment mode: (moments, None).  The scalar IS φ_0^0 = ``moments[0, 0]``
    # (Y_0^0 = 1), read off the tensor by the caller — NOT returned separately
    # (returning the live ``moment_buf[0, 0]`` view invites aliasing; the
    # angular-mode scalar is an independent array, so a None here keeps the two
    # modes' second slot from being mistaken for the same kind of value).
    return moment_buf, None


def _sweep_jacobi(
    Q: np.ndarray,
    sig_t: np.ndarray,
    sn_mesh: "SNMesh",
    boundary_flux: "BoundaryFlux",
    *,
    moment_frame: "FrameBase | None" = None,
    interior: "Callable",
) -> tuple[np.ndarray, np.ndarray]:
    r"""The bare multi-D sweep = the **Jacobi** octant schedule × one
    interior kernel (renamed from ``_sweep_2d_wavefront`` at C3.6 — the
    body has been d-generic since S6.4(d), and it is the JACOBI spelling,
    not a wavefront-specific one: all three multi-D representations'
    ``sweep`` doors route through here, each supplying its own interior).

    ONE group (all octants), NO inter-group reflect — delegates to the
    polymorphic :func:`_sweep_scheduled` with ``reflect=None``.  All octants
    read the same frozen inflow seed (**BARE**, Wave O #208 O.4b E1: the
    octant-incoming face slots come from the GIVEN ``boundary_flux``
    trace, no ``bc.apply`` — the reflective coupling ``B`` is delivered
    externally between sweeps, so this is the pure bulk solve
    ``ψ = (L+C)^{-1} q``).  The Gauss-Seidel SI resolvent calls the SAME
    ``_sweep_scheduled`` orchestrator with the per-in-plane-octant
    schedule + a ``−B`` reflect closure; Jacobi and G-S differ ONLY in
    the schedule object (the splitting is selected once, never by an
    ``if`` in the loop).

    S6.5 (#222): ``interior`` is REQUIRED — every caller names the
    representation instance whose kernel runs (production threads the
    operator's ONE instance; a defaulted kernel here would be a
    construction door outside ``default_for``).  Direct test consumers
    use the first-class ``MovingFrontierWindow(mesh).sweep(...)`` /
    ``ScanMarch(mesh).sweep(...)`` instead of this bare entry.

    Layout / history: ``Q`` is principled ``(N, ng, *spatial)``, ``sig_t``
    ``(ng, *spatial)`` (R-1 A1: producer-side-projected magnitude, no
    internal ``/W``; #196 PR-INDEX-5: principled face buffers,
    principled-equivalent to the legacy layout per vv-principles
    § bit-identity-vs-principled).
    """
    return _sweep_scheduled(
        Q, sig_t, sn_mesh, boundary_flux,
        schedule=SweepSchedule.jacobi(sn_mesh),
        reflect=None,
        moment_frame=moment_frame,
        interior=interior,
    )


__all__ = [
    # the selection layer
    "Compatibility",
    "IncompatibleRepresentation",
    "LossRepresentation",
    "LOSS_REPRESENTATIONS",
    "default_for",
    # the concrete representations
    "CumprodScan",
    "FullFieldWavefront",
    "MovingFrontierWindow",
    "ScanMarch",
    # the public sweep entry (relocated from the dissolved ``sweep.py``)
    "transport_sweep",
]
