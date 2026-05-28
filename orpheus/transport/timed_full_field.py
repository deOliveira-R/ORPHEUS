r"""Cross-method generic container: bulk field + boundary field + time history.

L2 (transport, method-agnostic). Holds the composite transport state
that an SN / CP / MoC / diffusion solver iterates over — a **bulk**
field paired with its **boundary** field partner, plus a rotating
history buffer for time-derivative stencils.

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
split is NOT SN-specific. Every transport method has the same
trio:

* **SN**: bulk = :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
  on ``(N, ng, nx, ny)``; boundary =
  :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux` on the
  flat face-trace layout.
* **CP** (future): bulk =
  :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux` on
  ``(ng, nx, ny)``; boundary = region-interface current (TBD when CP
  field migration lands in Phase 3.4+).
* **MoC** (future): bulk = ``RayField`` on ``(track, segment, ng)``;
  boundary = ``RayEndpointField`` per-track-family phase angles. See
  ``[[project_moc_structure]]`` for the fiber-bundle framing.
* **Diffusion** (future): bulk =
  :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`; boundary
  = surface-current ``J = -D ∇φ`` on :math:`\partial\Omega`.

:class:`TimedFullField` is the cross-method generic container that
holds this trio. The container is **field-like** (composes via Field-
inherited dunders propagated to both members) but is NOT a
:class:`~orpheus.numerics.field.Field` subclass at the typed-class
level — its natural backing space is a
:class:`~orpheus.numerics.space.DirectSumSpace` (the direct sum of
bulk and boundary spaces), and :class:`DirectSumSpace` is DEFERRED
to Phase 3.6 per the ``[[feedback_unify_after_two_instances]]`` rule
(unify once kinetics' ``flux ⊕ precursors`` ships the second use
case). When 3.6 lands, :class:`TimedFullField` may promote to a
``Field(values=DirectSumStorage, space=DirectSumSpace(...))`` form —
non-breaking, the public API stays the same.

Algebra contract
================

Same-class arithmetic propagates to both bulk and boundary members:

.. code-block:: python

    a + b = TimedFullField(bulk=a.bulk + b.bulk,
                           boundary=a.boundary + b.boundary,
                           _history=(),  # algebra drops iteration metadata
                           history_depth=a.history_depth)

History is iteration metadata, NOT algebraic state — algebra results
carry empty history. The ``advance(new_bulk, new_boundary)`` method is
the canonical state-step verb that rotates the history shift-register;
:meth:`at_lag` reads historical states for time-derivative stencils.

Cross-class arithmetic is rejected by an isinstance-and-bulk-type-
identity check in :meth:`_check_partner`: an SN ``TimedFullField``
cannot add to a CP ``TimedFullField`` even though both are
``TimedFullField`` instances, because their bulk types differ
(``AngularFlux`` vs ``ScalarFlux``).

Grep signal
===========

The name ``TimedFullField`` is a deliberately strong grep signal —
three load-bearing tokens with no false-positive matches in the
codebase:

* **Timed** — has history-aware time-derivative capability via
  :meth:`at_lag` / :meth:`advance`. Distinguishes from a stateless
  flux.
* **Full** — covers the full domain (bulk + boundary). Distinguishes
  from a bulk-only :class:`~orpheus.transport.fields.angular_flux.AngularFlux`.
* **Field** — algebraic composition object with dunder arithmetic.
  Distinguishes from non-algebraic state containers.

References
==========

* Depth B plan §3.8 (Container architecture). Located at
  ``.claude/plans/depth_b_field_on_function_space.md``.
* Grand Report v3 §5.5 (Field hierarchy), §5.3 (Space hierarchy
  including ``DirectSumSpace``), §6.1 (Space dunders).
* ``coding-elegance`` Pattern 5 (build the right primitive) — the
  cross-method container IS the right primitive; SN/CP/MoC/diffusion
  consumers all share it.
* ``[[project_transport_state_container]]`` — the original
  architectural decision memo (named ``TransportState`` pre-final-
  naming; renamed to ``TimedFullField`` 2026-05-28 per the user's
  grep-signal preference).
"""

from __future__ import annotations

from dataclasses import dataclass, field as dc_field, replace
from typing import TYPE_CHECKING

from orpheus.numerics.field import Field

if TYPE_CHECKING:
    from numpy.typing import NDArray


__all__ = ["TimedFullField"]


@dataclass(frozen=True, kw_only=True)
class TimedFullField:
    r"""Cross-method generic container: bulk field + boundary field + history.

    Parameters
    ----------
    bulk : Field
        The volumetric / bulk field. Any
        :class:`~orpheus.numerics.field.Field` subclass — typically
        :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
        for SN, :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
        for CP / diffusion, ``RayField`` for MoC.
    boundary : Field
        The boundary partner field on the trace of ``bulk``'s domain.
        Typically :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`
        for SN; method-specific for other methods.
    _history : tuple[TimedFullField, ...], optional
        Rotating history of prior states, most recent first.
        Default ``()`` (no history). Updated via :meth:`advance`.
    history_depth : int, optional
        Maximum number of historical states to retain. Default ``2``
        (current + one lag = enough for first-order time derivative).
        Bump for higher-order BDF stencils.

    Notes
    -----
    NOT a :class:`Field` subclass at the typed-class level — its
    natural backing is a :class:`DirectSumSpace` Field, which is
    deferred to Phase 3.6. Ships as a structured composite with
    delegate-style dunders that propagate to bulk + boundary; future
    P3.6 promotion is non-breaking.

    Algebra
    -------

    * Same-class ``+``, ``-``: propagate to bulk + boundary members.
      Algebra results have empty history (history is iteration
      metadata, not algebraic state).
    * Scalar ``*``, ``/``, unary ``-``: propagate to bulk + boundary.
    * Cross-class arithmetic is REJECTED: ``_check_partner`` enforces
      identical TimedFullField wrapping AND identical bulk/boundary
      class identity (catches SN-vs-CP cross-method composition even
      though both wrap as ``TimedFullField``).

    History
    -------

    * :meth:`advance(new_bulk, new_boundary)`: rotates the shift-
      register (current → lag-1, new → current).
    * :meth:`at_lag(k)`: reads the state at lag ``k`` (0 = current,
      1 = previous, ...). Used by BDF stencils:
      ``dpsi_dt ≈ (state.at_lag(0).bulk - state.at_lag(1).bulk) / dt``.
    """

    bulk: Field
    boundary: Field
    _history: tuple["TimedFullField", ...] = ()
    history_depth: int = 2

    # ── Construction validation ──────────────────────────────────────

    def __post_init__(self) -> None:
        if not isinstance(self.bulk, Field):
            raise TypeError(
                f"TimedFullField: bulk must be a Field; got "
                f"{type(self.bulk).__name__}"
            )
        if not isinstance(self.boundary, Field):
            raise TypeError(
                f"TimedFullField: boundary must be a Field; got "
                f"{type(self.boundary).__name__}"
            )
        # Mesh-identity check (where both members carry a ``mesh``
        # attribute — the cross-method generic contract). For SN both
        # AngularFlux and BoundaryFlux carry ``mesh: SNMesh``; other
        # methods follow the same convention.
        bulk_mesh = getattr(self.bulk, "mesh", None)
        boundary_mesh = getattr(self.boundary, "mesh", None)
        if bulk_mesh is not None and boundary_mesh is not None:
            if bulk_mesh is not boundary_mesh:
                raise ValueError(
                    "TimedFullField: bulk and boundary must share mesh "
                    "identity (both bound to the same mesh instance); "
                    f"got bulk.mesh={bulk_mesh!r}, "
                    f"boundary.mesh={boundary_mesh!r}"
                )
        if self.history_depth < 0:
            raise ValueError(
                f"TimedFullField: history_depth must be non-negative; "
                f"got {self.history_depth}"
            )

    # ── Algebra (propagates to bulk + boundary; history dropped) ─────

    def _check_partner(self, other: object) -> None:
        r"""Reject ill-formed binary operations.

        Layer 1 — class identity at the container level AND at the
        bulk / boundary member level. This catches cross-method
        composition (SN ``TimedFullField`` + CP ``TimedFullField``)
        where both operands wrap as ``TimedFullField`` but their bulk
        member types differ.
        """
        if type(other) is not TimedFullField:
            raise TypeError(
                f"TimedFullField arithmetic requires a same-class "
                f"partner; got {type(other).__name__}."
            )
        if type(self.bulk) is not type(other.bulk):  # type: ignore[attr-defined]
            raise TypeError(
                f"TimedFullField bulk type mismatch: "
                f"{type(self.bulk).__name__} vs "
                f"{type(other.bulk).__name__}."  # type: ignore[attr-defined]
            )
        if type(self.boundary) is not type(other.boundary):  # type: ignore[attr-defined]
            raise TypeError(
                f"TimedFullField boundary type mismatch: "
                f"{type(self.boundary).__name__} vs "
                f"{type(other.boundary).__name__}."  # type: ignore[attr-defined]
            )

    def __add__(self, other: "TimedFullField") -> "TimedFullField":
        self._check_partner(other)
        return TimedFullField(
            bulk=self.bulk + other.bulk,
            boundary=self.boundary + other.boundary,
            _history=(),
            history_depth=self.history_depth,
        )

    def __sub__(self, other: "TimedFullField") -> "TimedFullField":
        self._check_partner(other)
        return TimedFullField(
            bulk=self.bulk - other.bulk,
            boundary=self.boundary - other.boundary,
            _history=(),
            history_depth=self.history_depth,
        )

    def __neg__(self) -> "TimedFullField":
        return TimedFullField(
            bulk=-self.bulk,
            boundary=-self.boundary,
            _history=(),
            history_depth=self.history_depth,
        )

    def __mul__(self, scalar: float) -> "TimedFullField":
        return TimedFullField(
            bulk=self.bulk * float(scalar),
            boundary=self.boundary * float(scalar),
            _history=(),
            history_depth=self.history_depth,
        )

    def __rmul__(self, scalar: float) -> "TimedFullField":
        return self.__mul__(scalar)

    def __truediv__(self, scalar: float) -> "TimedFullField":
        return TimedFullField(
            bulk=self.bulk / float(scalar),
            boundary=self.boundary / float(scalar),
            _history=(),
            history_depth=self.history_depth,
        )

    # ── History shift-register (iteration metadata) ──────────────────

    def advance(
        self,
        new_bulk: Field,
        new_boundary: Field,
    ) -> "TimedFullField":
        r"""Push current ``(bulk, boundary)`` into history; install new as current.

        Returns a fresh :class:`TimedFullField` with:

        * ``bulk = new_bulk``, ``boundary = new_boundary`` (current frame)
        * ``_history = (current_snapshot, *self._history)[: history_depth]``
          where ``current_snapshot`` is the pre-advance state with
          empty inner history.

        The history-of-history is trimmed to ``history_depth`` (current
        frame is NOT counted in the depth — depth=2 means current +
        2 lags).

        Parameters
        ----------
        new_bulk : Field
            The new current bulk field. Must match ``type(self.bulk)``.
        new_boundary : Field
            The new current boundary field. Must match
            ``type(self.boundary)``.

        Returns
        -------
        TimedFullField
            The advanced state with rotated history.
        """
        if type(new_bulk) is not type(self.bulk):
            raise TypeError(
                f"TimedFullField.advance: new_bulk type "
                f"{type(new_bulk).__name__} does not match current "
                f"bulk type {type(self.bulk).__name__}."
            )
        if type(new_boundary) is not type(self.boundary):
            raise TypeError(
                f"TimedFullField.advance: new_boundary type "
                f"{type(new_boundary).__name__} does not match current "
                f"boundary type {type(self.boundary).__name__}."
            )
        current_snapshot = TimedFullField(
            bulk=self.bulk,
            boundary=self.boundary,
            _history=(),
            history_depth=self.history_depth,
        )
        new_history = (current_snapshot, *self._history)[: self.history_depth]
        return TimedFullField(
            bulk=new_bulk,
            boundary=new_boundary,
            _history=new_history,
            history_depth=self.history_depth,
        )

    def at_lag(self, lag: int) -> "TimedFullField":
        r"""Read the state at iteration lag.

        ``lag=0`` returns ``self`` (the current frame); ``lag=1``
        returns the most recently advanced-from frame; etc.

        Parameters
        ----------
        lag : int
            Non-negative lag depth. ``lag <= len(self._history)``.

        Returns
        -------
        TimedFullField
            The state at the requested lag.

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

    # ── Flat-vector protocol (Krylov / scipy.gmres adapter) ──────────
    #
    # Direct-sum flat representation: ``concat(bulk.values.ravel(),
    # boundary.values)``. The boundary values are already a flat 1-D
    # ndarray (per :class:`BoundaryFlux`'s flat-backing storage); the
    # bulk values are reshaped via :meth:`ndarray.ravel`. The Krylov
    # adapter at :mod:`orpheus.numerics.iteration` consumes this flat
    # representation as the GMRES iterate vector; round-trip exactness
    # is the load-bearing invariant.

    def to_flat(self) -> "NDArray":
        r"""Pack the current frame ``(bulk, boundary)`` into a flat 1-D vector.

        The packed layout is ``[bulk.values.ravel(), boundary.values]``
        — the direct-sum representation of the composite. The history
        is NOT packed (algebra results carry empty history; Krylov
        operates on the current frame only).

        Returns
        -------
        np.ndarray
            1-D ``float64`` ndarray of size
            ``bulk.values.size + boundary.values.size``.
        """
        import numpy as np

        return np.concatenate([
            self.bulk.values.ravel(),
            self.boundary.values,  # already 1-D (BoundaryFlux flat storage)
        ])

    @classmethod
    def from_flat(
        cls, flat: "NDArray", template: "TimedFullField",
    ) -> "TimedFullField":
        r"""Reconstruct a :class:`TimedFullField` from a flat 1-D vector + template.

        The ``template`` provides the shapes and types: ``bulk`` is
        reshaped to ``template.bulk.values.shape`` and constructed with
        the same ``space`` / ``mesh`` as the template; ``boundary``
        likewise. History is empty (Krylov iterates carry no
        iteration history).

        Parameters
        ----------
        flat : np.ndarray
            1-D vector matching ``template.to_flat()`` in size.
        template : TimedFullField
            Source of structural metadata (shapes, spaces, meshes,
            history_depth).

        Returns
        -------
        TimedFullField
            Reconstructed composite, empty history.
        """
        n_bulk = template.bulk.values.size
        n_boundary = template.boundary.values.size
        expected_total = n_bulk + n_boundary
        if flat.size != expected_total:
            raise ValueError(
                f"TimedFullField.from_flat: flat.size = {flat.size} "
                f"does not match template total size "
                f"{n_bulk} + {n_boundary} = {expected_total}"
            )
        bulk_values = flat[:n_bulk].reshape(template.bulk.values.shape)
        boundary_values = flat[n_bulk:]
        new_bulk = replace(template.bulk, values=bulk_values)
        new_boundary = replace(template.boundary, values=boundary_values)
        return cls(
            bulk=new_bulk,
            boundary=new_boundary,
            _history=(),
            history_depth=template.history_depth,
        )

    # ── Legacy adapter (D-H.1b transitional; retires in D-H.2) ───────

    @classmethod
    def from_legacy_angular_flux(
        cls, legacy_psi: object,
    ) -> "TimedFullField":
        r"""Construct from a legacy :class:`orpheus.sn.angular_flux.AngularFlux`.

        One-way adapter for D-H.1b consumer migration: extracts the
        legacy ``psi.values`` (bulk), ``psi.boundary`` (legacy
        BoundaryFlux), and ``psi.history_depth`` and constructs an L2
        composite. Internal calls use the legacy → L2 adapters on
        :class:`~orpheus.transport.fields.angular_flux.AngularFlux` and
        :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`.

        **Transitional**: this method retires when D-H.2 deletes the
        legacy :class:`orpheus.sn.angular_flux.AngularFlux` class. Do
        not use in new code — construct directly via
        :meth:`SNMesh.zeros_timed_full_field` or by passing typed L2
        bulk + boundary explicitly.

        Parameters
        ----------
        legacy_psi : object
            Legacy ``orpheus.sn.angular_flux.AngularFlux`` instance.
            Duck-typed: must expose ``values``, ``mesh``, ``boundary``,
            and ``history_depth`` attributes.

        Returns
        -------
        TimedFullField
            L2 composite with empty history (the legacy bundle's
            ``_history`` is NOT propagated — Krylov / SI iterates carry
            no history, and the legacy history was rarely populated in
            practice; future BDF kinetics work will use
            :meth:`advance` to build history).
        """
        from orpheus.transport.fields.angular_flux import AngularFlux as L2AngularFlux
        from orpheus.transport.fields.boundary_flux import BoundaryFlux as L2BoundaryFlux

        mesh = legacy_psi.mesh  # type: ignore[attr-defined]
        history_depth = getattr(legacy_psi, "history_depth", 2)
        return cls(
            bulk=L2AngularFlux.from_mesh(legacy_psi.values, mesh),  # type: ignore[attr-defined]
            boundary=L2BoundaryFlux.from_legacy_sn(
                legacy_psi.boundary, mesh,  # type: ignore[attr-defined]
            ),
            _history=(),
            history_depth=history_depth,
        )

    # ── Diagnostics ──────────────────────────────────────────────────

    def copy(self) -> "TimedFullField":
        r"""Return a deep copy of the current frame (history dropped).

        Snapshots the current ``(bulk, boundary)`` with owned ndarrays
        and empty history. Used by callers that need a stable iterate
        without the history-rotation aliasing.
        """
        return TimedFullField(
            bulk=self.bulk.copy(),
            boundary=self.boundary.copy(),
            _history=(),
            history_depth=self.history_depth,
        )
