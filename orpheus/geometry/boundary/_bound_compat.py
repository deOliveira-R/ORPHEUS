r"""Backward-compat wrapper around a realized 1-arg LinearOperator.

During the Waves 8-10 transition, SN-side ``bc_*`` attributes on
:class:`~orpheus.sn.geometry.SNMesh` are realized via
:class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer` into 1-arg
:class:`~orpheus.numerics.operator.LinearOperator` instances. The 13
production call sites at :mod:`orpheus.sn.sweep` +
:mod:`orpheus.sn.operator` still pass the 2-arg ``(psi, quad)``
signature through ``bc.apply(...)``. Wave 9 migrates those to 1-arg;
Wave 10 drops this shim. Until then, every realized BC at SNMesh
construction time is wrapped in :class:`_BoundBoundaryOperator` whose
``apply(psi, *_extra, **_kw)`` swallows the extra quadrature arg and
delegates to the realized op's 1-arg ``apply(psi)``.

The shim is intentionally minimal: it forwards :meth:`apply`,
:meth:`apply_transpose`, :attr:`capabilities`, and the standard
:class:`~orpheus.numerics.operator.LinearOperatorMixin` dunders.

Internal — not re-exported
==========================

This module is internal to the SN-side wiring at
:meth:`SNMesh._resolve_bcs`. It is NOT re-exported from
:mod:`orpheus.geometry.boundary` because no consumer outside that
construction site has any business wrapping operators in this shim
— if you find yourself reaching for it from production code, that
production code is the call site Wave 9-10 will migrate, and the
shim should NOT cement that legacy signature into the public API.

References
----------

* ``.claude/plans/transient-giggling-cake.md`` Wave 8 -- C8.3
  (this shim's brief in the SNMesh wiring), Wave 11 (deletion).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

from orpheus.numerics.operator import LinearOperatorMixin

if TYPE_CHECKING:
    from orpheus.numerics.operator import LinearOperator


__all__ = ["_BoundBoundaryOperator"]


class _BoundBoundaryOperator(LinearOperatorMixin):
    r"""Transitional 2-arg shim around a 1-arg realized operator.

    Forwards every ``apply(psi, *extra, **kw)`` call to
    ``inner.apply(psi)``, ignoring extra positional/keyword args
    (typically the legacy ``quadrature`` argument that Waves 9-10
    will remove). Same shape for :meth:`apply_transpose`.

    The :attr:`capabilities` property delegates to the wrapped
    operator, so a consumer composing this with other linear
    operators inherits whatever the realized op advertises (e.g. a
    :class:`PermutationOperator` brings ``apply_transpose``, a
    :class:`IncomingOrdinateMaskTensor` brings only ``apply``).

    The shim ALSO carries an optional ``kind`` tag (a free-form
    string like ``"vacuum"`` / ``"reflective"``) and implements
    ``__eq__`` against strings — preserves the legacy SN-side
    ``sn_mesh.bc_xmin == "reflective"`` comparison that several
    tests + the BC-resolution diagnostic rely on. Wave 9-10 may
    drop this string-comparison surface once those tests migrate
    to direct kind queries.

    Parameters
    ----------
    inner
        The realized :class:`LinearOperator` to wrap. Must support
        at least the :attr:`CAP_APPLY` capability.
    kind
        Optional free-form string tag carrying the originating
        :class:`BC` kind (``"vacuum"`` / ``"reflective"`` /
        ``"partial"`` / ...). When non-``None`` the shim compares
        equal to that string via :meth:`__eq__`.
    """

    def __init__(
        self, inner: "LinearOperator", kind: str | None = None,
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
        # Legacy string-kind comparison surface — preserves
        # ``sn_mesh.bc_xmin == "vacuum"`` style assertions.
        if isinstance(other, str):
            return self.kind == other
        return NotImplemented

    def __hash__(self) -> int:
        # Hashable by identity so it remains usable as a dict key /
        # set element. Using ``id(self)`` keeps each instance
        # distinct.
        return id(self)
