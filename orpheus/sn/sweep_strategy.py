r"""Selectable S\ :sub:`N` sweep strategies — polymorphic sweep/matvec dispatch.

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

This module replaces that with a first-class :class:`SweepStrategy`: each
algorithm is an object that carries **both** the forward ``sweep`` and (from
Phase S2) the ``residual`` matvec twin, plus a **declared, queryable**
:meth:`~SweepStrategy.supports` predicate.  The operator and the
``transport_sweep`` dispatcher select one via :func:`default_for` and then
call it branchlessly.

The hierarchy
=============

.. code-block:: text

    SweepStrategy (Protocol: sweep, residual [S2], supports)
    ├── _DAGWavefront            ── Cartesian anti-hyperplane DAG family
    │   ├── FullFieldWavefront     buffer = full field     · the ORACLE
    │   └── MovingFrontierWindow   buffer = rolling frontier · production opt
    └── CumprodScan             ── 1-D chain prefix scan, any geometry

``FullFieldWavefront`` and ``MovingFrontierWindow`` consume the **same**
per-octant ``sweep_graphs`` DAG — they are two *buffer policies* over one
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
  :meth:`~SweepStrategy.supports` / :func:`default_for`, independent of the
  code's capability.  "Don't pick the window at d=1, pick the scan" is a
  recommendation, *not* a reason to leave the window unable to express d=1.
* **Specialize on measured cost only.**  The sole justification to restrict
  an implementation's d-range is a *measured* hot-path regression.

Selection is a single source of truth
======================================

:meth:`~SweepStrategy.supports` returns :class:`Compatibility` — an
``(ok, reason)`` pair.  The same predicate serves three consumers:

#. a (future) teaching frontend — ``[S for S in SWEEP_STRATEGIES if
   S.supports(mesh).ok]`` grays-out an inapplicable method *and explains
   why* (pedagogically load-bearing — ORPHEUS teaches reactor physics);
#. the factory :func:`default_for` — picks the best *available* production
   optimization, falling back to the full-field spine when no optimization
   exists, so it is never stuck;
#. the construction guard — :meth:`_SweepStrategy.__post_init__` raises
   :class:`IncompatibleStrategy` on an illegal pairing, so even a bypassed
   UI cannot build one.

The compatibility signal is the genuine criterion — the coordinate system
(:attr:`SNMesh.is_cartesian`) and the dimensionality (:attr:`SNMesh.ndim`)
— NOT the ``sweep_graphs is None`` substrate proxy.

Phase status (the SweepStrategy carve, plan ``sn_sweep_strategy.md``)
====================================================================

* **S1 (this commit): SWEEP side.**  The protocol, the hierarchy, the three
  strategies as *thin wrappers* over the existing sweep functions, the
  selection layer, and the ``transport_sweep`` rewire.  Bit-identical: each
  strategy wraps exactly the path the old branch chose, so ``default_for``
  reproduces the legacy dispatch byte-for-byte.
* **S2: MATVEC side.**  Each strategy gains ``residual`` (the
  :math:`(L+C)\psi` twin); the five operator gates collapse to
  ``strategy.residual(...)``.
* **S3: generalize the oracle.**  ``FullFieldWavefront`` is the genuine
  d-generic spine (wraps :meth:`SweepDependencyGraph.apply` via the
  d-generic ``_sweep_full_field`` / ``_apply_full_field``, ``supports`` is
  any-d Cartesian, wires the d=1 cumprod-vs-spine equivalence).
* **S4: generalize the window** to ``frontier_dim = d-1``.

See also
========

* ``.claude/plans/sn_sweep_strategy.md`` — the authoritative design (the
  locked decisions, the verification strategy, phases S0–S5).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Protocol, runtime_checkable

# This module wraps the sweep bodies, so it imports ``sweep`` at load time.
# The back-edge — ``transport_sweep`` reaching here for ``default_for`` — is a
# function-local (lazy) import on the ``sweep`` side, which breaks the cycle.
from .sweep import (
    _sweep_1d_unified,
    _sweep_2d_scanmarch,
    _sweep_2d_wavefront,
    _sweep_full_field,
)

if TYPE_CHECKING:
    import numpy as np

    from orpheus.numerics.projection import MomentProjection
    from orpheus.transport.fields.angular_flux import AngularFlux
    from orpheus.transport.fields.boundary_flux import BoundaryFlux
    from orpheus.transport.timed_full_field import TimedFullField

    from .geometry import SNMesh
    from .operator import StreamingOperator


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


class IncompatibleStrategy(ValueError):
    """A :class:`SweepStrategy` was constructed for a mesh it cannot sweep.

    Raised by the construction guard (:meth:`_SweepStrategy.__post_init__`)
    so that an illegal ``(strategy, mesh)`` pairing is unrepresentable —
    even if a caller bypasses :func:`default_for`.
    """


@runtime_checkable
class SweepStrategy(Protocol):
    r"""One algorithm for the within-group transport solve and its twin.

    A strategy is constructed *for a mesh* (``Strategy(mesh)``); the
    construction guard rejects an incompatible pairing.  It then exposes:

    * :meth:`sweep` — one forward substitution :math:`\psi = (L+C)^{-1} q`;
    * :meth:`residual` — the matvec twin :math:`(L+C)\,\psi` *(added in
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

    def residual(
        self, operator: "StreamingOperator", psi: "TimedFullField",
    ) -> "TimedFullField":
        r"""The forward matvec twin :math:`L\,\psi` for this geometry.

        The sweep's operator-twin (L21 — sweep and matvec are different
        applications of the SAME operator): the sweep solves
        :math:`(L+C)^{-1} q`, this applies :math:`L`.  ``operator`` is the
        :class:`~orpheus.sn.operator.StreamingOperator` (it supplies
        :math:`\sigma_t` and the per-geometry matvec body the strategy
        selects).
        """
        ...

    def residual_transpose(
        self, operator: "StreamingOperator", phi: "TimedFullField",
    ) -> "TimedFullField":
        r"""The adjoint matvec twin :math:`L^{\mathsf T}\,\phi` for this geometry.

        Raises :class:`NotImplementedError` for strategies whose adjoint is
        deferred (the multi-D Cartesian reverse sweep — O.2b lands the 1-D
        reverse sweep first).  Never a silent wrong answer.
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
class _SweepStrategy:
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

    def residual(
        self, operator: "StreamingOperator", psi: "TimedFullField",
    ) -> "TimedFullField":
        """The forward matvec twin — every concrete strategy implements it."""
        raise NotImplementedError(
            f"{type(self).__name__} must implement residual()"
        )

    def residual_transpose(
        self, operator: "StreamingOperator", phi: "TimedFullField",
    ) -> "TimedFullField":
        """The adjoint matvec twin — implemented (1-D) or a deferral raise."""
        raise NotImplementedError(
            f"{type(self).__name__} must implement residual_transpose()"
        )

    def __post_init__(self) -> None:
        compat = type(self).supports(self.mesh)
        if not compat.ok:
            raise IncompatibleStrategy(
                f"{type(self).__name__} cannot sweep this mesh "
                f"(ndim={self.mesh.ndim}, curvature={self.mesh.curvature!r}): "
                f"{compat.reason}."
            )


# ═══════════════════════════════════════════════════════════════════════
# CumprodScan — the 1-D chain prefix scan (any geometry)
# ═══════════════════════════════════════════════════════════════════════


class CumprodScan(_SweepStrategy):
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

    def residual(
        self, operator: "StreamingOperator", psi: "TimedFullField",
    ) -> "TimedFullField":
        """1-D forward matvec — the operator's geometry-blind ``L·ψ``."""
        return operator._apply_1d(psi)

    def residual_transpose(
        self, operator: "StreamingOperator", phi: "TimedFullField",
    ) -> "TimedFullField":
        """1-D adjoint matvec — the operator's reverse-sweep ``Lᵀ·φ``."""
        return operator._apply_1d_transpose(phi)


# ═══════════════════════════════════════════════════════════════════════
# _DAGWavefront — the Cartesian anti-hyperplane DAG family
# ═══════════════════════════════════════════════════════════════════════


class _DAGWavefront(_SweepStrategy):
    r"""Base for the two buffer policies over the per-octant DAG walk.

    ``FullFieldWavefront`` (full-field buffer; the oracle) and
    ``MovingFrontierWindow`` (rolling :math:`(d{-}1)`-frontier; the
    production optimization) both walk the **same** ``mesh.sweep_graphs``
    anti-hyperplane DAG with the same diamond-difference cell kernel.  They
    differ only in *how much* of the interior face cochain they retain — a
    storage policy, pinned bit-identical by the ``window ≡ full`` oracle.

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

    def residual_transpose(
        self, operator: "StreamingOperator", phi: "TimedFullField",
    ) -> "TimedFullField":
        r"""The multi-D Cartesian adjoint is DEFERRED (shared by both policies).

        The forward matvec works for both buffer policies, but the
        reverse-direction adjoint sweep is not yet wired (O.2b landed the 1-D
        reverse sweep first).  Raises :class:`NotImplementedError` — the mesh
        is compatible, only the adjoint *feature* is deferred (so this is NOT
        an :class:`IncompatibleStrategy`).  Never a silent wrong answer.
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
        )

    def residual(
        self, operator: "StreamingOperator", psi: "TimedFullField",
    ) -> "TimedFullField":
        """2-D forward matvec — the operator's windowed ``L·ψ`` (storage-B).

        Wraps :meth:`~orpheus.sn.operator.StreamingOperator._apply_2d_cartesian`
        (the rolling-frontier ``residual_windowed`` walk — the matvec twin of
        :func:`._sweep_2d_wavefront`; ONE discretization, L21).  WRAP not
        relocate, so the A2D-1 source-hash pin stays a free tripwire.
        """
        return operator._apply_2d_cartesian(psi)


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

    Wraps the d-generic :func:`._sweep_full_field` / the operator's
    :meth:`~orpheus.sn.operator.StreamingOperator._apply_full_field` (both walk
    :meth:`SweepDependencyGraph.apply` / ``.residual`` directly).  ``supports``
    is any-d Cartesian (S3) — the spine is genuinely dimension-generic, unlike
    the d=2-only window.
    """

    @classmethod
    def supports(cls, mesh: "SNMesh") -> Compatibility:
        # Override the _DAGWavefront family's d=2-only predicate: the spine is
        # the genuine d-generic oracle (it walks the per-octant DAG for any
        # Cartesian d via the d-generic ``graph.apply``/``.residual``).
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

    def residual(
        self, operator: "StreamingOperator", psi: "TimedFullField",
    ) -> "TimedFullField":
        """Forward matvec ORACLE — the full-field ``L·ψ`` (d-generic).

        Wraps
        :meth:`~orpheus.sn.operator.StreamingOperator._apply_full_field`
        (the full interior-cochain walk, matvec twin of
        :func:`._sweep_full_field`).  The fuller-view reference the windowed
        matvec is cross-checked against; never the production default.
        """
        return operator._apply_full_field(psi)


# ═══════════════════════════════════════════════════════════════════════
# ScanMarch — the row-march + x-scan schedule (1-D scan ∘ transverse march)
# ═══════════════════════════════════════════════════════════════════════


class ScanMarch(_SweepStrategy):
    r"""Scan-march sweep — ``scan(x)`` marched over the transverse axes (#222).

    Reframes the d-D diamond-difference sweep as forward substitution along the
    sweep axis — the first-order linear scan
    :func:`~orpheus.sn.spatial.scan.ordinate_scan` — marched over the transverse
    axes: ``scan(x)`` at d=1, ``scan(x) ∘ march(y)`` at d=2.  ONE primitive that
    **unifies** the 1-D :class:`CumprodScan` (its degenerate ``s_y = 0`` case)
    and the 2-D row-march: the within-row x-face recurrence is the SAME Blelloch
    scan, the transverse coupling rides the affine source
    (:func:`~orpheus.sn.sweep._sweep_2d_scanmarch`).

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
    (``residual``) + the d≥3 recursive transverse march land in S5.1b/S5.2.
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
        return _sweep_2d_scanmarch(
            Q, sig_t, self.mesh, boundary_flux,
            moment_projection=moment_projection,
        )

    def residual(
        self, operator: "StreamingOperator", psi: "TimedFullField",
    ) -> "TimedFullField":
        """Forward matvec twin — 1-D wired; the 2-D scan-march matvec is S5.1b."""
        if self.mesh.is_1d:
            return operator._apply_1d(psi)
        raise NotImplementedError(
            "ScanMarch.residual: the 2-D scan-march matvec twin lands in S5.1b "
            "(a new StreamingOperator._apply_2d_cartesian_scanmarch, keeping the "
            "A2D-1 source-hash free-green); the forward sweep is verified against "
            "the FullFieldWavefront oracle first."
        )

    def residual_transpose(
        self, operator: "StreamingOperator", phi: "TimedFullField",
    ) -> "TimedFullField":
        """Adjoint matvec twin — 1-D wired; the multi-D Cartesian adjoint is deferred."""
        if self.mesh.is_1d:
            return operator._apply_1d_transpose(phi)
        raise NotImplementedError(
            "ScanMarch.residual_transpose: the multi-D Cartesian adjoint is "
            "deferred (O.2b lands the 1-D reverse sweep first; the multi-D "
            "adjoint follows the forward scan-march matvec, S5.1b+)."
        )


# ═══════════════════════════════════════════════════════════════════════
# Registry + factory — the single selection source of truth
# ═══════════════════════════════════════════════════════════════════════

#: The CONCRETE strategy leaves (never the abstract ``_SweepStrategy`` /
#: ``_DAGWavefront`` bases) — :func:`default_for` constructs whichever it
#: picks, so only buildable strategies belong here.
#:
#: Selection priority order: the 1-D scan, the 2-D production window, the
#: d-general scan-march, then the full-field oracle.  :func:`default_for`
#: returns the FIRST that applies, so the legacy production optimizations win
#: at d=1/d=2 (the scan-march is OPT-IN — issue #222 Fork B1, default
#: unchanged), the scan-march is the d≥3 Cartesian production primitive, and
#: the oracle is the never-stuck final fallback.
SWEEP_STRATEGIES: tuple[type[_SweepStrategy], ...] = (
    CumprodScan,
    MovingFrontierWindow,
    ScanMarch,
    FullFieldWavefront,
)


def default_for(mesh: "SNMesh") -> SweepStrategy:
    """Select the default sweep strategy for ``mesh``.

    Returns the first strategy in :data:`SWEEP_STRATEGIES` whose
    :meth:`~SweepStrategy.supports` admits ``mesh`` — the best *available*
    production optimization, falling back to the spine.  This reproduces the
    legacy dispatch exactly: 1-D → :class:`CumprodScan`; 2-D Cartesian →
    :class:`MovingFrontierWindow`.

    Raises
    ------
    IncompatibleStrategy
        If no strategy applies.  Unreachable for any constructible mesh
        (every 1-D mesh → ``CumprodScan``; every 2-D Cartesian →
        ``MovingFrontierWindow``).  A hypothetical 3-D Cartesian mesh selects
        the d-general ``ScanMarch`` (production), with ``FullFieldWavefront``
        (the oracle) as the never-stuck fallback.
    """
    for cls in SWEEP_STRATEGIES:
        if cls.supports(mesh).ok:
            return cls(mesh)
    raise IncompatibleStrategy(
        f"no sweep strategy supports this mesh "
        f"(ndim={mesh.ndim}, curvature={mesh.curvature!r}, "
        f"is_cartesian={mesh.is_cartesian})."
    )
