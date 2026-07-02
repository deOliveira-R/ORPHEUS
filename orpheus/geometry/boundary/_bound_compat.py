r"""Strict 1-arg passthrough wrapping a realized BC operator with a kind tag.

After Issue #188 (curvilinear trace support) every
:class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` BC routes through
:class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` to a 1-arg
:class:`~orpheus.numerics.operator.LinearOperator`. After Issue #186
(B3 + β2) the descriptor model removed the last 2-arg call sites
in tests, so this shim is now strict 1-arg. Its purpose is to add
two thin surfaces to a realized operator:

* **String-equality compatibility** (``bc == "vacuum"``) via the
  :attr:`kind` attribute. Several SN-side tests + the BC-resolution
  diagnostic in :mod:`orpheus.sn.solver` rely on this comparison.

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
:meth:`SNMesh._resolve_bcs`. It is NOT re-exported from
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

import numpy as np

from orpheus.numerics.operator import (
    LinearOperator,
    MissingAdjoint,
    NotInvertible,
    adjointable,
    invertible,
)


__all__ = ["_BoundBoundaryOperator"]


class _BoundBoundaryOperator(LinearOperator):
    r"""Strict 1-arg passthrough wrapping a realized BC operator with a kind tag.

    The wrapped :attr:`inner` is a 1-arg realized
    :class:`~orpheus.numerics.operator.LinearOperator` produced by
    :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer`. The
    shim's :meth:`apply` and :meth:`apply_transpose` forward ``psi``
    to the inner operator. Issue #186 (B3 + β2) cleanup: extra
    positional / keyword args raise :class:`TypeError` rather than
    being silently swallowed.

    The structural predicates, ``apply_transpose``, and the
    :class:`LinearOperator` dunders delegate to :attr:`inner`.
    The optional :attr:`kind` tag carries the originating
    :class:`~orpheus.geometry.mesh.BC` kind string and is the basis
    for the ``sn_mesh.bc["xmin"] == "reflective"`` style comparisons that
    several SN tests and the BC-resolution diagnostic rely on.

    Parameters
    ----------
    inner
        Realized 1-arg :class:`LinearOperator`.
    kind
        Optional free-form string tag (``"vacuum"`` / ``"reflective"``
        / etc.). When non-``None``, the shim compares equal to that
        string via :meth:`__eq__`.
    """

    def __init__(
        self,
        inner: "LinearOperator",
        kind: str | None = None,
    ) -> None:
        self.inner = inner
        self.kind = kind

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
