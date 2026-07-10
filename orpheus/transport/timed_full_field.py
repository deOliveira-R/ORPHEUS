r"""Cross-method generic container: the timed (history-bearing) full field.

L2 (transport, method-agnostic). :class:`TimedFullField` is the
**cofree comonad** ``Cofree(FullField, depth=d)`` over the timeless
:class:`~orpheus.transport.full_field.FullField` carrier — it pairs the
*current* timeless ``bulk ⊕ boundary`` frame with a rotating history
buffer of prior timeless frames. This is the composite transport state
an SN / CP / MoC / diffusion solver **iterates** over.

The split (the #217 deliverable)
================================

The timeless algebra — the ``bulk ⊕ boundary`` pair, the six
vector-space dunders, the flat-vector protocol (``to_flat`` /
``from_flat``), ``copy``, ``zeros``, and the bulk/boundary/mesh
construction validation — lives ONCE in the base
:class:`~orpheus.transport.full_field.FullField` (DRY). This class ADDS
only the time-aware part:

* ``_history`` / ``history_depth`` — the rotating shift-register.
* :meth:`advance` / :meth:`at_lag` / :meth:`history_length` — the
  shift-register verbs (the only time-aware behaviour; only the drivers
  use them).
* The ``history_depth >= 0`` validation (extends the base's
  bulk/boundary/mesh checks).
* The :meth:`_recombine` override that carries ``history_depth`` and an
  EMPTY history into algebra results, plus the ``zeros`` / ``from_flat``
  classmethod specialisations that thread ``history_depth``.

Only the iteration / time-stepping drivers see the comonad (the
history); the operator algebra reads a timeless
:class:`~orpheus.transport.full_field.FullField` frame in and writes one
out (it is blind to the history). That is the structural reason the
split is forced — see :mod:`orpheus.transport.full_field` for the full
rationale.

Why the bulk + boundary + history split is the right L2 abstraction
====================================================================

The pre-D-H :class:`~orpheus.sn.angular_flux.AngularFlux` conflated
THREE concepts on a single class:

1. **Volumetric (bulk) flux values** — :math:`\psi(\vec r, \hat\Omega, g)`
   on the interior cell-centred grid.
2. **Boundary trace values** — :math:`\psi|_\Gamma` on the inflow /
   outflow faces (required across iterations for reflective BCs).
3. **Iteration / time history** — the rotating :math:`(\psi^{n-1},
   \psi^{n-2}, \dots)` buffer used by BDF time-derivative stencils
   and by Anderson-acceleration / Krylov inner-product reductions.

Per Cardinal Rule 2 (shared concepts → shared abstraction), this
split is NOT SN-specific. Concepts (1)+(2) are the timeless
:class:`~orpheus.transport.full_field.FullField`; concept (3) is the
history buffer this class adds. Every transport method has the same
trio:

* **SN**: bulk =
  :class:`~orpheus.transport.fields.angular_flux.AngularFlux` on
  ``(N, ng, nx, ny)``; boundary =
  :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux` on the
  flat face-trace layout.
* **CP** (future): bulk =
  :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux` on
  ``(ng, nx, ny)``; boundary = region-interface current (TBD).
* **MoC** (future): bulk = ``RayField`` on ``(track, segment, ng)``;
  boundary = ``RayEndpointField`` per-track-family phase angles. See
  ``[[project_moc_structure]]`` for the fiber-bundle framing.
* **Diffusion** (future): bulk =
  :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`; boundary
  = surface-current ``J = -D ∇φ`` on :math:`\partial\Omega`.

Algebra contract
================

Same-class arithmetic propagates to both bulk and boundary members and
carries an empty history (history is iteration metadata, NOT algebraic
state):

.. code-block:: python

    a + b = TimedFullField(interior=a.interior + b.interior,
                           boundary=a.boundary + b.boundary,
                           _history=(),  # algebra drops iteration metadata
                           history_depth=a.history_depth)

This is realised via the inherited dunders routed through the
overridden :meth:`_recombine`, so the algebra is written once on
:class:`~orpheus.transport.full_field.FullField`. The
``advance(new_bulk, new_boundary)`` method is the canonical state-step
verb that rotates the history shift-register; :meth:`at_lag` reads
historical states for time-derivative stencils.

Cross-class arithmetic is rejected by an isinstance-and-bulk-type-
identity check in :meth:`~orpheus.transport.full_field.FullField._check_partner`:
an SN ``TimedFullField`` cannot add to a CP ``TimedFullField`` even
though both are ``TimedFullField`` instances, because their bulk types
differ (``AngularFlux`` vs ``ScalarFlux``).

Grep signal
===========

The name ``TimedFullField`` is a deliberately strong grep signal —
three load-bearing tokens with no false-positive matches in the
codebase:

* **Timed** — has history-aware time-derivative capability via
  :meth:`at_lag` / :meth:`advance`. Distinguishes from the timeless
  :class:`~orpheus.transport.full_field.FullField` base.
* **Full** — covers the full domain (bulk + boundary). Distinguishes
  from a bulk-only :class:`~orpheus.transport.fields.angular_flux.AngularFlux`.
* **Field** — algebraic composition object with dunder arithmetic.
  Distinguishes from non-algebraic state containers.

References
==========

* GH **issue #217** — the timeless-``FullField`` extraction (the base
  this class now inherits).
* :mod:`orpheus.transport.full_field` — the timeless base (carries the
  algebra + the cofree-comonad rationale).
* Depth B plan §3.8 (Container architecture). Located at
  ``.claude/plans/depth_b_field_on_function_space.md``.
* ``[[project_transport_state_container]]`` — the original
  architectural decision memo (named ``TransportState`` pre-final-
  naming; renamed to ``TimedFullField`` 2026-05-28 per the user's
  grep-signal preference). A SEPARATE, later ``TransportState`` —
  the #257 S4 structural Protocol — was also retired (by #257 S4.5,
  superseded by the concrete :class:`~orpheus.transport.full_field.FullField`
  base); the name ``TransportState`` is now unused in the codebase.
"""

from __future__ import annotations

from dataclasses import dataclass

from orpheus.transport.fields._bases import (
    BulkField,
    BoundaryField,
    RadialCharacteristicField,
)
from orpheus.transport.full_field import FullField

__all__ = ["TimedFullField"]


@dataclass(frozen=True, kw_only=True)
class TimedFullField(FullField):
    r"""Cofree comonad over :class:`~orpheus.transport.full_field.FullField`.

    The history-bearing composite — a current timeless ``bulk ⊕
    boundary`` frame plus a rotating history of prior frames. The
    iterating flux state of every transport solve.

    Parameters
    ----------
    interior : BulkField
        The volumetric / bulk field (inherited from
        :class:`~orpheus.transport.full_field.FullField`). Any
        :class:`~orpheus.transport.fields._bases.BulkField` subclass.
    boundary : AngularBoundaryField
        The boundary partner field on the trace of ``interior``'s domain
        (inherited).
    _history : tuple[FullField, ...], optional
        Rotating history of prior frames, most recent first. Default
        ``()`` (no history). Updated via :meth:`advance`. The elements
        are timeless :class:`~orpheus.transport.full_field.FullField`
        snapshots (a historical frame is itself timeless).
    history_depth : int, optional
        Maximum number of historical frames to retain. Default ``2``
        (current + two lags). Bump for higher-order BDF stencils.

    Notes
    -----
    The algebra (the six vector-space dunders, ``to_flat`` /
    ``from_flat``, ``copy``, ``zeros``, the bulk/boundary/mesh
    validation) lives on the
    :class:`~orpheus.transport.full_field.FullField` base. This class
    overrides :meth:`_recombine` so algebra results are
    ``TimedFullField`` with empty history + preserved ``history_depth``
    (#217), and specialises :meth:`zeros` / :meth:`from_flat` to thread
    ``history_depth``.

    History
    -------

    * :meth:`advance(new_bulk, new_boundary)`: rotates the shift-
      register (current → lag-1, new → current).
    * :meth:`at_lag(k)`: reads the state at lag ``k`` (0 = current,
      1 = previous, ...). Used by BDF stencils:
      ``dpsi_dt ≈ (state.at_lag(0).interior - state.at_lag(1).interior) / dt``.
    """

    _history: tuple["FullField", ...] = ()
    history_depth: int = 2

    # ── Construction ─────────────────────────────────────────────────

    @classmethod
    def zeros(
        cls,
        *,
        interior: type[BulkField],
        boundary: type[BoundaryField],
        mesh: "object",
        history_depth: int = 2,
        radial_characteristic: type[RadialCharacteristicField] | None = None,
    ) -> "TimedFullField":
        r"""Allocate a zero composite from the bulk + boundary leaf TYPES (B.5.A).

        Specialises
        :meth:`~orpheus.transport.full_field.FullField.zeros` by threading
        ``history_depth`` (the timeless base has no such concept). The
        bulk / boundary leaf CLASSES are zero-allocated on ``mesh`` via
        their own :meth:`zeros_on` — the SN-specific ``(AngularFlux,
        AngularBoundaryFlux)`` composition lives at the SN call site, not here.
        Replaces the retired ``SNMesh.zeros_timed_full_field``.

        Parameters
        ----------
        interior : type[BulkField]
            The bulk leaf CLASS to instantiate (must expose
            ``zeros_on(mesh)``).
        boundary : type[BoundaryField]
            The boundary leaf CLASS to instantiate (must expose
            ``zeros_on(mesh)``).
        mesh : object
            The phase-space carrier passed through to each leaf's
            ``zeros_on`` (duck-typed — no transport→mesh hard
            dependency).
        history_depth : int, optional
            History buffer depth (default 2; see the class docstring).
        radial_characteristic : type[RadialCharacteristicField], optional
            The starting-direction leaf CLASS; presence is MESH-keyed
            (allocated iff ``mesh.radial_characteristic_space`` is
            non-``None`` — the R12a predicate). See
            :meth:`FullField.zeros`.
        """
        # Delegate leaf zero-allocation to the timeless base (the SINGLE
        # source of the duck-typed ``zeros_on`` calls AND of the
        # mesh-keyed seed-presence resolution), then re-wrap with the
        # history metadata — so neither lives duplicated here.
        base = FullField.zeros(
            interior=interior,
            boundary=boundary,
            mesh=mesh,
            radial_characteristic=radial_characteristic,
        )
        return cls(
            interior=base.interior,
            boundary=base.boundary,
            radial_characteristic=base.radial_characteristic,
            _history=(),
            history_depth=history_depth,
        )

    # ── Construction validation ──────────────────────────────────────

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.history_depth < 0:
            raise ValueError(
                f"TimedFullField: history_depth must be non-negative; "
                f"got {self.history_depth}"
            )

    # ── Polymorphic recombine hook (carries history metadata) ────────

    def _recombine(
        self,
        *,
        interior: BulkField,
        boundary: BoundaryField,
        radial_characteristic: RadialCharacteristicField | None,
    ) -> "TimedFullField":
        r"""Rebuild a ``TimedFullField`` with empty history + preserved depth.

        Overrides
        :meth:`~orpheus.transport.full_field.FullField._recombine` so the
        inherited vector-space dunders (defined ONCE on the base) return
        a :class:`TimedFullField` for a ``TimedFullField`` operand —
        carrying ``history_depth`` and an EMPTY history (#217: algebra
        results carry empty history; history is iteration metadata, not
        algebraic state). ``radial_characteristic`` is the REQUIRED keyword
        of the base hook — see the base docstring's silent-drop note.
        """
        return TimedFullField(
            interior=interior,
            boundary=boundary,
            radial_characteristic=radial_characteristic,
            _history=(),
            history_depth=self.history_depth,
        )

    # ── History shift-register (iteration metadata) ──────────────────

    def advance(
        self,
        new_bulk: BulkField,
        new_boundary: BoundaryField,
        new_radial_characteristic: RadialCharacteristicField | None = None,
    ) -> "TimedFullField":
        r"""Push the current frame into history; install the new one as current.

        Returns a fresh :class:`TimedFullField` with:

        * ``bulk = new_bulk``, ``boundary = new_boundary`` (and, on a
          seed-carrying composite, ``radial_characteristic =
          new_radial_characteristic``) — the current frame;
        * ``_history = (current_snapshot, *self._history)[: history_depth]``
          where ``current_snapshot`` is the pre-advance timeless frame.

        The history-of-history is trimmed to ``history_depth`` (current
        frame is NOT counted in the depth — depth=2 means current +
        2 lags). A historical frame is a timeless
        :class:`~orpheus.transport.full_field.FullField` snapshot of the
        blocks at that step.

        Parameters
        ----------
        new_bulk : BulkField
            The new current bulk field. Must match ``type(self.interior)``.
        new_boundary : AngularBoundaryField
            The new current boundary field. Must match
            ``type(self.boundary)``.
        new_radial_characteristic : RadialCharacteristicField, optional
            The new current ψ½ block. Presence must MATCH the current
            frame's (a seed-carrying iterate cannot silently drop its
            block, an unseeded one cannot grow it — R12a presence is a
            structural fact of the mesh, not of the step), and the type
            must match when present.

        Returns
        -------
        TimedFullField
            The advanced state with rotated history.
        """
        if type(new_bulk) is not type(self.interior):
            raise TypeError(
                f"TimedFullField.advance: new_bulk type "
                f"{type(new_bulk).__name__} does not match current "
                f"bulk type {type(self.interior).__name__}."
            )
        if type(new_boundary) is not type(self.boundary):
            raise TypeError(
                f"TimedFullField.advance: new_boundary type "
                f"{type(new_boundary).__name__} does not match current "
                f"boundary type {type(self.boundary).__name__}."
            )
        if (self.radial_characteristic is None) != (new_radial_characteristic is None):
            raise TypeError(
                f"TimedFullField.advance: starting-direction presence must "
                f"match the current frame (current: "
                f"{self.radial_characteristic is not None}, new: "
                f"{new_radial_characteristic is not None}) — presence is a "
                f"structural fact of the mesh (R12a), not of the step."
            )
        if new_radial_characteristic is not None and type(
            new_radial_characteristic
        ) is not type(self.radial_characteristic):
            raise TypeError(
                f"TimedFullField.advance: new_radial_characteristic type "
                f"{type(new_radial_characteristic).__name__} does not match "
                f"current type {type(self.radial_characteristic).__name__}."
            )
        current_snapshot = FullField(
            interior=self.interior,
            boundary=self.boundary,
            radial_characteristic=self.radial_characteristic,
        )
        new_history = (current_snapshot, *self._history)[: self.history_depth]
        return TimedFullField(
            interior=new_bulk,
            boundary=new_boundary,
            radial_characteristic=new_radial_characteristic,
            _history=new_history,
            history_depth=self.history_depth,
        )

    def at_lag(self, lag: int) -> "FullField":
        r"""Read the frame at iteration lag.

        ``lag=0`` returns ``self`` (the current frame); ``lag=1``
        returns the most recently advanced-from frame (a timeless
        :class:`~orpheus.transport.full_field.FullField` snapshot); etc.

        Parameters
        ----------
        lag : int
            Non-negative lag depth. ``lag <= len(self._history)``.

        Returns
        -------
        FullField
            The frame at the requested lag (``self`` at lag 0, a timeless
            ``FullField`` snapshot otherwise).

        Raises
        ------
        IndexError
            If ``lag`` exceeds the available history depth.
        ValueError
            If ``lag`` is negative.
        """
        if lag < 0:
            raise ValueError(
                f"TimedFullField.at_lag: lag must be non-negative; "
                f"got {lag}."
            )
        if lag == 0:
            return self
        if lag - 1 >= len(self._history):
            raise IndexError(
                f"TimedFullField.at_lag({lag}): only "
                f"{len(self._history)} historical states retained "
                f"(history_depth={self.history_depth})."
            )
        return self._history[lag - 1]

    @property
    def history_length(self) -> int:
        r"""Number of historical states currently retained."""
        return len(self._history)

    # ``from_flat`` is INHERITED from
    # :meth:`~orpheus.transport.full_field.FullField.from_flat`: the base
    # is generic over the template type and routes reconstruction through
    # :meth:`_recombine`, so ``TimedFullField.from_flat`` already returns a
    # ``TimedFullField`` with the template's ``history_depth`` + empty
    # history. No override needed (and no Liskov param-narrowing ignore).
