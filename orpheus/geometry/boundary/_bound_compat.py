r"""The realized boundary law — a realized operator PAIRED with its descriptor.

After Issue #188 (curvilinear trace support) every
:class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` BC routes through
:class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` to a 1-arg
:class:`~orpheus.numerics.operator.LinearOperator`. After Issue #186
(B3 + β2) the descriptor model removed the last 2-arg call sites
in tests, so this shim is now strict 1-arg. It adds three surfaces to
a realized operator:

* **The law it was realized from** (:attr:`law`, campaign phase B2.0).
  The three-layer architecture is descriptor → realizer → operator, and
  before B2.0 the descriptor was *dropped* on the way out: the shim kept
  ``kind=law.key``, a **string**, and nothing else. So ``sn_mesh.bc[face]``
  could answer *"what were you declared as?"* but not *"what do you
  DO?"* — and every consumer needing the latter had to re-derive it
  from the tag. Carrying the law makes the structural questions
  answerable at the object (``law.geometry_map.permutes_ordinates``,
  ``law.response_kernel.is_zero``) instead of by string comparison.

* **String-equality compatibility** (``bc == "vacuum"``) via the
  :attr:`kind` property. Several SN-side tests and five production
  sites rely on this comparison; phase B2.2 repoints them at the law's
  factors and retires the surface. Until then :attr:`kind` is a
  *read-through* of the law's registry key — no longer a stored copy
  that could drift from the law it claims to describe.

* **Structural-predicate delegation** to the wrapped operator so
  consumers composing the shim with other Wave-0 primitives inherit
  the right two-axis truth (e.g. a :class:`PermutationOperator` brings
  invertibility + a working ``apply_transpose``, an
  :class:`IncomingOrdinateMaskTensor` only adjointability).

History
=======

The original Wave-8/9 implementation carried an optional
``quadrature=`` kwarg that, when non-``None``, bound an
AngularQuadrature and forwarded ``inner.apply(psi, bound_quad)`` to
a legacy 2-arg :class:`BoundaryTraceLaw` body. That mode existed
ONLY because Wave 2 deferred curvilinear trace support, which forced
:class:`SNMesh._resolve_bcs` to bypass the realizer for spherical /
cylindrical meshes. Issue #188 lifted that deferral;
Issue #176 dropped the bound-quadrature mode; Issue #186 then
dropped the ``*_extra, **_kw`` swallow on :meth:`apply` since the
underlying need (legacy 2-arg test call sites) was eliminated by
the descriptor-model cleanup.

Internal — not re-exported
==========================

This module is internal to the SN-side wiring at
``SNMesh.realize_boundary_law`` (the SN arm of the #290 P7b
``TransportMethod`` hook). It is NOT re-exported from
:mod:`orpheus.geometry.boundary` because no consumer outside that
construction site should wrap operators in this shim.

References
----------

* ``.claude/plans/transient-giggling-cake.md`` Wave 8 / C8.3
  (original brief — dual-mode shim).
* ``.claude/plans/curvilinear-realizer-and-2arg-cleanup.md`` Phase 2
  / C176.1 (dual-mode → 1-arg-with-extras).
* ``.claude/plans/bc-trace-law-descriptor-cleanup.md`` C-B3.4
  (this commit — strict 1-arg).
* Issues #188 (curvilinear support unblock), #176 (legacy 2-arg cleanup),
  and #186 (descriptor model + final tightening).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

import numpy as np

from orpheus.numerics.operator import (
    LinearOperator,
    MissingAdjoint,
    NotInvertible,
    adjointable,
    invertible,
)


if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace

    from ._base import BoundaryTraceLaw


__all__ = ["_BoundBoundaryOperator"]


class _BoundBoundaryOperator(LinearOperator):
    r"""A realized BC operator paired with the law it was realized from.

    The wrapped :attr:`inner` is a 1-arg realized
    :class:`~orpheus.numerics.operator.LinearOperator` produced by
    :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer`. The
    shim's :meth:`apply` and :meth:`apply_transpose` forward ``psi``
    to the inner operator. Issue #186 (B3 + β2) cleanup: extra
    positional / keyword args raise :class:`TypeError` rather than
    being silently swallowed.

    The structural predicates, ``apply_transpose``, ``domain`` /
    ``codomain``, and the :class:`LinearOperator` dunders delegate to
    :attr:`inner`; :attr:`law` is the descriptor :attr:`inner` was
    realized from, and :attr:`kind` reads that law's registry key.

    Parameters
    ----------
    inner
        Realized 1-arg :class:`LinearOperator` — what the law *does*
        on this face's trace space.
    law
        The :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`
        ``inner`` was realized from — what the law *means*, and the
        only thing that can answer a structural question about it.
        **Required**: a realized boundary law that cannot say which law
        it realizes is exactly the state phase B2 exists to delete, so
        it is not constructible.
    """

    def __init__(
        self,
        inner: "LinearOperator",
        law: "BoundaryTraceLaw",
    ) -> None:
        self.inner = inner
        self.law = law

    @property
    def kind(self) -> Optional[str]:
        r"""The realized law's REGISTRY KEY — a read-through, not a copy.

        Reads ``type(self.law).key``, deliberately **not**
        :attr:`~orpheus.geometry.boundary.BoundaryTraceLaw.kind`. The two
        diverge for exactly one law: a partially-reflecting
        :class:`~orpheus.geometry.boundary.reflective.ReflectiveBoundary`
        (``albedo != 1``) reports ``kind == "partial"`` — mirroring the
        ``BC("partial", albedo=…)`` *declaration* vocabulary that
        ``BC.to_alpha`` accepts (B0.1 ruling) — while its ``key`` stays
        ``"reflective"`` for every albedo.

        The key is what the pre-B2.0 shim stored, so this property is
        behaviour-identical to it; reading ``law.kind`` here would
        silently drop partially-reflecting faces out of the sweep
        schedule's reflective set. **Phase B2.2 retires this surface
        entirely** — the five production sites that compare it against
        string literals are asking structural questions that
        :attr:`law`'s two factors answer directly.

        ``None`` for a law with no registry slot (an ad-hoc subclass
        declared without ``key=``); the string comparison in
        :meth:`__eq__` is then False against every string.
        """
        return type(self.law).key

    # ── Function-space forwarding ────────────────────────────────────
    #
    # A transparent handle must forward the space tags too: an operator
    # composed with this shim type-checks against ``inner``'s spaces, and
    # the base ``LinearOperator`` default of ``None`` would silently skip
    # the composability check (see ``LinearOperator.domain``) rather than
    # perform it. Every other wrapper in the codebase forwards these;
    # this one did not until B2.0.

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        return getattr(self.inner, "domain", None)

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return getattr(self.inner, "codomain", None)

    @property
    def is_invertible(self) -> bool:
        # Delegate the structural predicates to the realized inner law —
        # the wrapper is a transparent typed handle, so ``B``'s per-face
        # truth rides on the inner op's own predicates, NOT the base
        # ``LinearOperator`` default ``False`` (which would silently break
        # the boundary aggregator's ``all(law.is_adjointable …)`` rule).
        return self.inner.is_invertible

    @property
    def is_adjointable(self) -> bool:
        return self.inner.is_adjointable

    def inverse(self) -> "LinearOperator":
        r"""Forward :meth:`inverse` to the realized inner law (#226 step 1).

        The shim is a transparent typed handle: whenever the inner law
        advertises ``is_invertible`` (a reflective
        :class:`~orpheus.numerics.operator.PermutationOperator` does), its
        inverse OPERATOR must be reachable through the shim exactly as the
        predicate is — otherwise the shim would advertise a promise
        (``is_invertible=True``) it cannot deliver. Raises the inner law's
        own guard on a non-invertible law (a vacuum mask projects).
        """
        if not invertible(self.inner):
            raise NotInvertible(
                f"_BoundBoundaryOperator.inverse: the realized "
                f"{type(self.inner).__name__} law is not invertible."
            )
        return self.inner.inverse()

    @property
    def block_role(self):  # type: ignore[override]
        # Forward the realized law's block-role classification (Issue #208
        # / Wave O) so ``isinstance(sn_mesh.bc["xmin"], BoundaryOperator)``
        # reads the inner op's role. The realized boundary laws carry
        # ``BlockRole.BOUNDARY``; the rank-0 affine PrescribedInflow source
        # carries ``None`` (it is ``q.boundary``, not a linear ``B``).
        # ``getattr``-safe: an inner op without the attribute reports None.
        return getattr(self.inner, "block_role", None)

    def apply(self, psi: np.ndarray) -> np.ndarray:
        return self.inner.apply(psi)

    def apply_transpose(self, psi: np.ndarray) -> np.ndarray:
        # Guard-narrow licenses the delegated call (Design C); the shim
        # stays transparent — the refusal is the inner law's truth.
        if not adjointable(self.inner):
            raise MissingAdjoint(
                f"_BoundBoundaryOperator.apply_transpose: the realized "
                f"{type(self.inner).__name__} law is not adjointable."
            )
        return self.inner.apply_transpose(psi)

    # NO ``solve`` (carve P4 — the promised deletion): solving with the
    # shim is ``.inverse().apply(b)``, the forwarded inner law's inverse.

    def __eq__(self, other: object) -> bool:
        if isinstance(other, str):
            return self.kind == other
        return NotImplemented

    def __hash__(self) -> int:
        return id(self)
