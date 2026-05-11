r"""1-arg passthrough wrapping a realized BC operator with a kind tag.

After Issue #188 (curvilinear `InflowTraceSpace` support) every
:class:`~orpheus.sn.geometry.SNMesh` BC routes through
:class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer` to a 1-arg
:class:`~orpheus.numerics.operator.LinearOperator`. The dual-mode
bound-quadrature variant the Wave-8/9 shim used to bridge the
curvilinear bypass is no longer needed; this shim is now a thin
passthrough that adds three things to a realized operator:

* **String-equality compatibility** (``bc == "vacuum"``) via the
  :attr:`kind` attribute. Several SN-side tests + the BC-resolution
  diagnostic in :mod:`orpheus.sn.solver` rely on this comparison.

* **Pass-through of extra positional / keyword args** to
  :meth:`apply` / :meth:`apply_transpose`. Any test that still calls
  ``bc.apply(psi, quad)`` (the legacy 2-arg form before Issue #176
  finished cleaning everything up) keeps working; the extra args are
  silently ignored because the realized inner op only needs ``psi``.

* **Capabilities delegation** to the wrapped operator so consumers
  composing the shim with other Wave-0 primitives inherit the right
  feature set (e.g. a :class:`PermutationOperator` brings
  ``apply_transpose``, an :class:`IncomingOrdinateMaskTensor` brings
  only ``apply``).

History
=======

The original Wave-8/9 implementation carried an optional
``quadrature=`` kwarg that, when non-``None``, bound an
AngularQuadrature and forwarded ``inner.apply(psi, bound_quad)`` to
a legacy 2-arg :class:`BoundaryTraceLaw` body. That mode existed
ONLY because Wave 2 deferred curvilinear :class:`InflowTraceSpace`,
which forced :class:`SNMesh._resolve_bcs` to bypass the realizer for
spherical / cylindrical meshes. Issue #188 lifted that deferral;
Issue #176 then dropped the bound-quadrature mode from this shim
because no production-issued shim carries
``_quadrature is not None`` after C188.3.

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
  / C176.1 (this commit — dual-mode → 1-arg-only).
* Issues #188 (curvilinear support unblock) and #176 (this cleanup).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

from orpheus.numerics.operator import LinearOperatorMixin

if TYPE_CHECKING:
    from orpheus.numerics.operator import LinearOperator


__all__ = ["_BoundBoundaryOperator"]


class _BoundBoundaryOperator(LinearOperatorMixin):
    r"""1-arg passthrough wrapping a realized BC operator with a kind tag.

    The wrapped :attr:`inner` is a 1-arg realized
    :class:`~orpheus.numerics.operator.LinearOperator` produced by
    :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer`. The
    shim's :meth:`apply` forwards ``psi`` to ``inner.apply(psi)`` and
    swallows any extra positional / keyword args — the latter
    accommodates lingering ``bc.apply(psi, quad)`` test sites during
    the Issue #176 cleanup.

    Capabilities, ``apply_transpose``, and the ``LinearOperatorMixin``
    dunders delegate to :attr:`inner`. The optional :attr:`kind` tag
    carries the originating :class:`~orpheus.geometry.mesh.BC`
    kind string and is the basis for the
    ``sn_mesh.bc_xmin == "reflective"`` style comparisons that
    several legacy SN tests and the BC-resolution diagnostic rely
    on.

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
    def capabilities(self) -> frozenset[str]:  # type: ignore[override]
        return self.inner.capabilities

    def apply(self, psi: np.ndarray, *_extra: Any, **_kw: Any) -> np.ndarray:
        return self.inner.apply(psi)

    def apply_transpose(
        self, psi: np.ndarray, *_extra: Any, **_kw: Any,
    ) -> np.ndarray:
        return self.inner.apply_transpose(psi)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, str):
            return self.kind == other
        return NotImplemented

    def __hash__(self) -> int:
        return id(self)
