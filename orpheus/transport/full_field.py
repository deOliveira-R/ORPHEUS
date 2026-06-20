r"""The timeless composite carrier: bulk field ⊕ boundary field.

L2 (transport, method-agnostic). :class:`FullField` is the **timeless**
``bulk ⊕ boundary`` carrier that the transport operator algebra acts on:
a flux flows in, a source / residual flows out, and *no time* enters the
picture. It is the object-image of the operator algebra
:math:`(L + C - S - F - B)\,\psi = q` — every operator leaf maps a
:class:`FullField` to a :class:`FullField` (the codomain inner role may
change, flux → source, but the composite carrier type does not).

The cofree-comonad framing (the #217 split)
============================================

:class:`~orpheus.transport.timed_full_field.TimedFullField` is the
**cofree comonad** ``Cofree(FullField, depth=d)`` over this base: it
pairs the *current* timeless frame with a rotating history buffer of
prior timeless frames. Only the iteration / time-stepping drivers see
the comonad (the history); the operator algebra is blind to it — it
reads a timeless frame in and writes a timeless frame out.

This is why the split is structurally FORCED rather than aesthetic:

* A **static source** ``q = q_bulk ⊕ q_∂`` has no time — it never
  ``advance``\ s. Typing it as a history-bearing ``TimedFullField``
  hands it verbs (``advance`` / ``at_lag``) it must never use ("a type
  error of altitude" — the value carries an unused history tail).
* An **iterating flux state** ``ψ^n`` DOES advance through iterations /
  time steps — it is the comonad, ``TimedFullField``.

So the operator-algebra carrier is the timeless :class:`FullField`; the
driver-level iterate is the timed :class:`TimedFullField`. This module
holds the algebra ONCE (DRY) — :class:`TimedFullField` inherits it.

Why the bulk + boundary split is the right L2 abstraction
=========================================================

The pre-D-H :class:`~orpheus.sn.angular_flux.AngularFlux` conflated
volumetric flux values with boundary trace values. Per Cardinal Rule 2
(shared concepts → shared abstraction), the bulk + boundary split is
NOT SN-specific — every transport method has the same pair:

* **SN**: bulk =
  :class:`~orpheus.transport.fields.angular_flux.AngularFlux` on
  ``(N, ng, *spatial)``; boundary =
  :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux` on the
  flat face-trace layout.
* **CP / MoC / diffusion** (future): their own bulk + boundary leaves.

Algebra contract
================

Same-class arithmetic propagates to both bulk and boundary members:

.. code-block:: python

    a + b = FullField(bulk=a.bulk + b.bulk,
                      boundary=a.boundary + b.boundary)

The six vector-space dunders (``+``, ``-``, unary ``-``, scalar ``*``,
``scalar *``, ``/ scalar``) live here ONCE and route the recombined
``(bulk, boundary)`` pair through :meth:`_recombine` — a single
polymorphic hook. :class:`FullField`'s hook rebuilds a bare
:class:`FullField`; :class:`TimedFullField` OVERRIDES it to rebuild a
``TimedFullField`` with empty history and preserved ``history_depth``
(#217: "algebra results carry empty history"). One definition of the
algebra, the correct concrete return type for each subclass.

Cross-class arithmetic is rejected at two layers: :meth:`_check_partner`
rejects a partner whose type is not ``type(self)`` (an SN ``FullField``
cannot add to a CP ``FullField`` whose bulk type differs — caught at the
member level too); the member-level leaf dunders enforce role / units /
mesh / space (the #208 affine torsor gate ``flux + flux → TypeError``,
the cross-mesh guard) by delegation, so the composite does NOT
pre-check member types (that would block the legitimate composite
torsor ``flux + displacement → flux``).

Grep signal
===========

``FullField`` — the *full* domain (bulk + boundary), *timeless*. Its
history-bearing subclass keeps the strong three-token grep signal
``TimedFullField`` (Timed + Full + Field).

References
==========

* GH **issue #217** — the timeless-``FullField`` extraction this
  delivers (the composite source is the first timeless consumer).
* ``.claude/plans/issue_257_coefficient_field_promotion.md`` §S4.5 —
  the locked cofree-comonad design.
* Grand Report v3 §5.5 (Field hierarchy), §5.3 (``DirectSumSpace`` /
  :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`).
* ``coding-elegance`` Pattern 2 (the algebra lives ONCE in the base via
  the recombine hook) + Pattern 5 (``FullField`` is the right
  primitive; ``TimedFullField`` composes it with a history buffer).
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, TypeVar

import numpy as np

from orpheus.transport.fields._bases import BoundaryField, BulkField

if TYPE_CHECKING:
    from numpy.typing import NDArray


__all__ = ["FullField"]

#: Self-type for the composite dunders. Bound to :class:`FullField`, so a
#: subclass (:class:`~orpheus.transport.timed_full_field.TimedFullField`)
#: inherits the algebra typed as ``(other: Subclass) -> Subclass`` — the
#: self-type idiom mirrors :class:`~orpheus.numerics.field.Field`'s ``T``.
#: This is what lets ``TimedFullField`` satisfy the
#: :class:`~orpheus.numerics.vector.Vector` ``Self`` contract while the
#: algebra is defined ONCE on the base.
T = TypeVar("T", bound="FullField")


@dataclass(frozen=True, kw_only=True)
class FullField:
    r"""Timeless composite carrier: a bulk field paired with its boundary partner.

    Parameters
    ----------
    bulk : BulkField
        The volumetric / bulk field. Any
        :class:`~orpheus.transport.fields._bases.BulkField` subclass —
        typically :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
        for SN, :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
        for CP / diffusion, ``RayField`` for MoC.
    boundary : BoundaryField
        The boundary partner field on the trace of ``bulk``'s domain.
        Typically :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`
        for SN; method-specific for other methods.

    Notes
    -----
    NOT a :class:`~orpheus.numerics.field.Field` subclass at the
    typed-class level — its natural backing is a
    :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`
    (direct-sum) Field, deferred to Phase 3.6. Ships as a structured
    composite with delegate-style dunders that propagate to bulk +
    boundary; the future P3.6 promotion is non-breaking.

    Algebra
    -------

    * Same-class ``+``, ``-``: propagate to bulk + boundary members.
    * Scalar ``*``, ``/``, unary ``-``: propagate to bulk + boundary.
    * Cross-class arithmetic is REJECTED: :meth:`_check_partner` enforces
      identical ``type(self)`` wrapping; the member leaves enforce
      role / units / mesh / space by delegation.

    The six dunders route their recombined pair through the
    :meth:`_recombine` hook so a :class:`TimedFullField` operand yields a
    :class:`TimedFullField` (empty history), not a bare
    :class:`FullField`.
    """

    bulk: BulkField
    boundary: BoundaryField

    # ── Construction ─────────────────────────────────────────────────

    @classmethod
    def zeros(
        cls,
        *,
        bulk: type[BulkField],
        boundary: type[BoundaryField],
        mesh: "object",
    ) -> "FullField":
        r"""Allocate a zero timeless composite from the bulk + boundary leaf TYPES.

        Generic over the method's leaf types: the caller passes the bulk
        and boundary :class:`~orpheus.numerics.field.Field` *subclasses*
        (SN passes
        :class:`~orpheus.transport.fields.angular_flux.AngularFlux` /
        :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`; CP
        / MoC will pass their own), and each is zero-allocated on
        ``mesh`` via its own :meth:`zeros_on`. This keeps the
        cross-method-generic container free of any hard-wired leaf type.

        Parameters
        ----------
        bulk : type[BulkField]
            The bulk leaf CLASS to instantiate (must expose
            ``zeros_on(mesh)``).
        boundary : type[BoundaryField]
            The boundary leaf CLASS to instantiate (must expose
            ``zeros_on(mesh)``).
        mesh : object
            The phase-space carrier passed through to each leaf's
            ``zeros_on`` (duck-typed — no transport→mesh hard
            dependency).
        """
        return cls(
            bulk=bulk.zeros_on(mesh),  # type: ignore[attr-defined]
            boundary=boundary.zeros_on(mesh),  # type: ignore[attr-defined]
        )

    # ── Construction validation ──────────────────────────────────────

    def __post_init__(self) -> None:
        if not isinstance(self.bulk, BulkField):
            raise TypeError(
                f"{type(self).__name__}: bulk must be a BulkField; got "
                f"{type(self.bulk).__name__}"
            )
        if not isinstance(self.boundary, BoundaryField):
            raise TypeError(
                f"{type(self).__name__}: boundary must be a BoundaryField; "
                f"got {type(self.boundary).__name__}"
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
                    f"{type(self).__name__}: bulk and boundary must share "
                    "mesh identity (both bound to the same mesh instance); "
                    f"got bulk.mesh={bulk_mesh!r}, "
                    f"boundary.mesh={boundary_mesh!r}"
                )

    # ── Polymorphic recombine hook (Pattern 2 — the algebra lives once) ─

    def _recombine(self: T, *, bulk: BulkField, boundary: BoundaryField) -> T:
        r"""Rebuild a composite of the SAME concrete type from a recombined pair.

        The single polymorphic hook the six vector-space dunders route
        their result through. The base spelling is ``replace(self, ...)``
        — provably ``T`` (a :class:`FullField` recombines to a
        :class:`FullField`, no history to drop).
        :class:`TimedFullField` OVERRIDES this to rebuild a
        ``TimedFullField`` carrying its ``history_depth`` and an EMPTY
        history (#217: algebra results carry empty history — ``replace``
        would copy the history, which is wrong for the timed subclass).
        This keeps the dunders defined ONCE while preserving the correct
        concrete return type per subclass (DRY + Liskov-correct).
        """
        return replace(self, bulk=bulk, boundary=boundary)

    # ── Algebra (propagates to bulk + boundary via _recombine) ───────

    def _check_partner(self, other: object) -> None:
        r"""Reject a partner that is not a :class:`FullField` (timed or timeless).

        Layer 1 at the CONTAINER level only. The bulk / boundary
        member-level algebra (the affine gate ``flux + flux → TypeError``,
        the torsor ``flux + displacement → flux``, the displacement mint
        ``flux − flux → displacement``, cross-class rejection, the
        cross-mesh guard) is the SINGLE SOURCE OF TRUTH on the leaves —
        ``__add__`` / ``__sub__`` delegate to ``self.bulk ± other.bulk``
        and ``self.boundary ± other.boundary``, where the leaf dunders
        enforce role-correctness (#208). Pre-checking member types here
        would BLOCK the legitimate composite torsor (``flux`` bulk +
        ``displacement`` bulk), so the redundant member pre-checks are
        intentionally absent.

        The accepted partner is ANY :class:`FullField` — a timeless base
        or a timed subclass. This is load-bearing for the
        time-derivative stencil ``state.at_lag(0) - state.at_lag(1)``: the
        current frame is a :class:`TimedFullField` while a historical
        frame from :meth:`~orpheus.transport.timed_full_field.TimedFullField.at_lag`
        is a timeless :class:`FullField` snapshot, and the two must
        subtract. The CONCRETE result type is governed by
        :meth:`_recombine` (``timed − timeless → timed`` via the timed
        hook; ``timeless − timed → timeless`` via the base hook), so
        accepting the base partner does not weaken the "algebra results
        carry empty history" guarantee. ``state + 42`` (a non-``FullField``
        partner) is still rejected here; the leaf gate then rejects any
        bulk/boundary/mesh/units mismatch by delegation.
        """
        if not isinstance(other, FullField):
            raise TypeError(
                f"{type(self).__name__} arithmetic requires a same-class "
                f"partner; got {type(other).__name__}."
            )

    def __add__(self: T, other: T) -> T:
        self._check_partner(other)
        return self._recombine(
            bulk=self.bulk + other.bulk,
            boundary=self.boundary + other.boundary,
        )

    def __sub__(self: T, other: T) -> T:
        self._check_partner(other)
        return self._recombine(
            bulk=self.bulk - other.bulk,
            boundary=self.boundary - other.boundary,
        )

    def __neg__(self: T) -> T:
        return self._recombine(bulk=-self.bulk, boundary=-self.boundary)

    def __mul__(self: T, scalar: float) -> T:
        return self._recombine(
            bulk=self.bulk * float(scalar),
            boundary=self.boundary * float(scalar),
        )

    def __rmul__(self: T, scalar: float) -> T:
        return self.__mul__(scalar)

    def __truediv__(self: T, scalar: float) -> T:
        return self._recombine(
            bulk=self.bulk / float(scalar),
            boundary=self.boundary / float(scalar),
        )

    # ── Flat-vector protocol (Krylov / scipy.gmres adapter) ──────────
    #
    # Direct-sum flat representation: ``concat(bulk.values.ravel(),
    # boundary.values)``. The boundary values are already a flat 1-D
    # ndarray (per :class:`BoundaryFlux`'s flat-backing storage); the
    # bulk values are reshaped via :meth:`ndarray.ravel`. The Krylov
    # adapter at :mod:`orpheus.numerics.iteration` consumes this flat
    # representation as the GMRES iterate vector; round-trip exactness
    # is the load-bearing invariant.
    #
    # NB taking ``boundary.values`` WITHOUT a ravel encodes the SN
    # convention that ``BoundaryFlux`` stores an already-flat trace; a
    # future method whose boundary leaf is not flat-backed would
    # generalise this to ``boundary.values.ravel()`` (deferred per
    # ``coding-elegance`` Pattern 6 — one boundary leaf today).

    def to_flat(self) -> "NDArray":
        r"""Pack ``(bulk, boundary)`` into a flat 1-D vector.

        The packed layout is ``[bulk.values.ravel(), boundary.values]``
        — the direct-sum representation of the composite.

        Returns
        -------
        np.ndarray
            1-D ``float64`` ndarray of size
            ``bulk.values.size + boundary.values.size``.
        """
        return np.concatenate([
            self.bulk.values.ravel(),
            self.boundary.values,  # already 1-D (BoundaryFlux flat storage)
        ])

    @classmethod
    def from_flat(
        cls, flat: "NDArray", template: T,
    ) -> T:
        r"""Reconstruct a composite from a flat 1-D vector + template.

        The ``template`` provides the shapes, types, AND the concrete
        composite class: ``bulk`` is reshaped to
        ``template.bulk.values.shape`` and rebuilt with the same ``space``
        / ``mesh`` as the template; ``boundary`` likewise; the pair is
        then routed through ``template``'s :meth:`_recombine` hook, so the
        result is the SAME concrete type as ``template`` — a
        :class:`TimedFullField` template yields a ``TimedFullField`` with
        preserved ``history_depth`` and an EMPTY history (Krylov iterates
        carry no history), a timeless :class:`FullField` template yields a
        :class:`FullField`. One definition, correct concrete type per
        subclass (``coding-elegance`` Pattern 2 — the same hook the
        algebra dunders use).

        Parameters
        ----------
        flat : np.ndarray
            1-D vector matching ``template.to_flat()`` in size.
        template : FullField
            Source of structural metadata (shapes, spaces, meshes) AND the
            concrete return type.

        Returns
        -------
        FullField
            Reconstructed composite of the SAME concrete type as
            ``template``.
        """
        n_bulk = template.bulk.values.size
        n_boundary = template.boundary.values.size
        expected_total = n_bulk + n_boundary
        if flat.size != expected_total:
            raise ValueError(
                f"{cls.__name__}.from_flat: flat.size = {flat.size} "
                f"does not match template total size "
                f"{n_bulk} + {n_boundary} = {expected_total}"
            )
        bulk_values = flat[:n_bulk].reshape(template.bulk.values.shape)
        boundary_values = flat[n_bulk:]
        new_bulk = replace(template.bulk, values=bulk_values)
        new_boundary = replace(template.boundary, values=boundary_values)
        return template._recombine(bulk=new_bulk, boundary=new_boundary)

    # ── Diagnostics ──────────────────────────────────────────────────

    def copy(self: T) -> T:
        r"""Return a deep copy with owned ndarrays.

        Snapshots ``(bulk, boundary)`` with owned ndarrays. Used by
        callers that need a stable iterate without aliasing. Routes
        through :meth:`_recombine`, so a :class:`TimedFullField` copy is a
        :class:`TimedFullField` with EMPTY history (the existing
        ``copy`` drops history — bit-identical behaviour).
        """
        return self._recombine(
            bulk=self.bulk.copy(),
            boundary=self.boundary.copy(),
        )
