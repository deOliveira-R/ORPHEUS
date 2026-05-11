r"""Backward-compat wrapper around realized 1-arg / legacy 2-arg ops.

During the Waves 8-10 transition, SN-side ``bc_*`` attributes on
:class:`~orpheus.sn.geometry.SNMesh` present a UNIFORM 1-arg
``bc.apply(psi)`` contract to the production call sites at
:mod:`orpheus.sn.sweep` + :mod:`orpheus.sn.operator`. Two concrete
backings hide behind that contract, both wrapped in this shim:

* **Cartesian path (Wave 8 / C8.3):** ``inner`` is a 1-arg realized
  :class:`~orpheus.numerics.operator.LinearOperator` produced by
  :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer`. The shim
  forwards ``apply(psi, *_extra, **_kw) → inner.apply(psi)`` and
  swallows any extra args.

* **Curvilinear path (Wave 9 / C9.0):** :class:`InflowTraceSpace` does
  not yet support spherical / cylindrical meshes (Wave 2 deferral), so
  the SN realizer skips and ``_resolve_one`` returns the raw legacy
  :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` subclass
  instance whose ``apply`` is 2-arg. To keep the 1-arg contract uniform
  across ALL ``bc_*`` attributes, the shim is built with a
  ``quadrature=`` kwarg that binds the AngularQuadrature once; the
  shim then forwards ``apply(psi) → inner.apply(psi, bound_quad)``.
  The bound quadrature is bit-identical to the legacy direct call
  because it is the SAME object the call sites used to pass.

Wave 9 migrates the 13 call sites from ``bc.apply(psi, quad)`` to
``bc.apply(psi)``; Wave 10 drops this shim entirely once
:class:`InflowTraceSpace` covers curvilinear geometry and the
:class:`BoundaryTraceLaw` ``apply`` bodies are dropped to 1-arg.

The shim forwards :meth:`apply`, :meth:`apply_transpose`,
:attr:`capabilities`, and the standard
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
    from orpheus.sn.quadrature import AngularQuadrature


__all__ = ["_BoundBoundaryOperator"]


class _BoundBoundaryOperator(LinearOperatorMixin):
    r"""Transitional shim presenting a uniform 1-arg ``apply`` contract.

    Two backings (see module docstring):

    * ``quadrature=None`` (Cartesian path): ``inner`` is a 1-arg
      realized :class:`LinearOperator`; ``apply(psi) → inner.apply(psi)``.
    * ``quadrature=<AngularQuadrature>`` (curvilinear path): ``inner``
      is a legacy 2-arg :class:`BoundaryTraceLaw` instance; the bound
      quadrature is forwarded as the 2nd positional arg so the call
      site sees a uniform 1-arg signature ``apply(psi)``.

    In both modes the shim's ``apply(psi, *_extra, **_kw)`` accepts
    extra positional / keyword args and ignores them — Wave 9
    intermediate states where some call sites still pass ``quad`` keep
    working while the migration is mid-flight. Same shape for
    :meth:`apply_transpose`.

    The :attr:`capabilities` property delegates to the wrapped
    operator, so a consumer composing this with other linear
    operators inherits whatever the realized op advertises (e.g. a
    :class:`PermutationOperator` brings ``apply_transpose``, an
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
        The wrapped operator. For the Cartesian path: a 1-arg
        realized :class:`LinearOperator`. For the curvilinear path:
        a 2-arg :class:`BoundaryTraceLaw` instance.
    quadrature
        Optional :class:`AngularQuadrature` bound to the shim. When
        non-``None`` the shim's :meth:`apply` forwards
        ``inner.apply(psi, quadrature)``; when ``None`` it forwards
        ``inner.apply(psi)``.
    kind
        Optional free-form string tag carrying the originating
        :class:`BC` kind (``"vacuum"`` / ``"reflective"`` /
        ``"partial"`` / ...). When non-``None`` the shim compares
        equal to that string via :meth:`__eq__`.
    """

    def __init__(
        self,
        inner: "LinearOperator",
        quadrature: "AngularQuadrature | None" = None,
        kind: str | None = None,
    ) -> None:
        self.inner = inner
        self._quadrature = quadrature
        self.kind = kind

    @property
    def capabilities(self) -> frozenset[str]:  # type: ignore[override]
        return self.inner.capabilities

    def apply(self, psi: np.ndarray, *_extra: Any, **_kw: Any) -> np.ndarray:
        if self._quadrature is None:
            return self.inner.apply(psi)
        return self.inner.apply(psi, self._quadrature)

    def apply_transpose(
        self, psi: np.ndarray, *_extra: Any, **_kw: Any,
    ) -> np.ndarray:
        if self._quadrature is None:
            return self.inner.apply_transpose(psi)
        return self.inner.apply_transpose(psi, self._quadrature)

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
