r"""Per-ordinate residual field :math:`r(\vec r, \hat\Omega_n, g)`.

The L2 typed wrapper for the **per-ordinate transport-balance
residual** — the quantity produced when the discrete transport
operator is applied to a flux and the right-hand source is subtracted:

.. math::

    r_n \;=\; \bigl[(L + C - S - F)\,\psi\bigr]_n \;-\; q_n .

Same storage shape as
:class:`~orpheus.transport.fields.angular_flux.AngularFlux` and
:class:`~orpheus.transport.sources.angular_source.AngularSource`
(``(N, ng, nx, ny)``), but a structurally distinct physical quantity:
it is the *defect* of the balance, not a flux and not an external
source.

Why a distinct class (the "dimensional sin")
============================================

Pre-refactor, ``AngularFlux`` was the carrier for FOUR physically
distinct quantities at once: the flux :math:`\psi`
(:math:`1/(\mathrm{cm^2 \cdot s \cdot sr \cdot eV})`), the operator
outputs :math:`L\psi / S\psi / F\psi`, the external source
:math:`q_{\rm ext}`, and this iteration residual (the latter three all
:math:`1/(\mathrm{cm^3 \cdot s \cdot sr \cdot eV})`). The epicentre is
``numerics/iteration.py`` (``rhs = F.apply(ψ) + S.apply(ψ) + q_ext``)
and ``sn/operator.py``. Collapsing four quantities onto one type
erases the distinctions physics makes and is the bug surface this role
grid removes.

:class:`AngularResidual` and
:class:`~orpheus.transport.sources.angular_source.AngularSource` share
a dimensional signature (both rate densities,
:math:`1/(\mathrm{cm^3 \cdot s \cdot sr \cdot eV})`) — yet they are
DIFFERENT classes. Field's Layer-1 class-identity gate
(:meth:`~orpheus.numerics.field.Field._check_partner`) therefore
rejects ``angular_residual - angular_source`` even though the units
align: *same units give permission to add in linear algebra; they do
not give meaning*. The balance that combines a residual and a source
goes through an explicit **named composition** (planned
``IterationResidual.from_balance(lhs, rhs)``, B.5 / Issue #201), never
a bare cross-class ``-``. See :class:`~orpheus.numerics.field.Field`
"Class identity for cross-class same-units operations".

Units (B.4 — declared as the ``UNITS`` class constant)
======================================================

Per-ordinate reaction-rate-density defect:
:math:`[1/(\mathrm{cm^3 \cdot s \cdot sr})]`
(:data:`~orpheus.numerics.units.ANGULAR_RATE_UNITS`) — the
``AngularFlux`` units times an inverse length (the streaming
:math:`\hat\Omega\cdot\nabla` and collision :math:`\sigma_t` both
carry :math:`1/\mathrm{cm}`); eV-free per the binned-energy convention.
Under View-G the units are NOT a space property; they ARE the
role-leaf's ``UNITS`` class constant (the design commitment that
*class identity is units identity*). See :mod:`orpheus.numerics.units`.

References
==========

* ``.claude/plans/field_role_typing_view_g.md`` — Phase B (field
  vocabulary), step B.3 (new role leaves); B.5 (named-composition
  ``IterationResidual.from_balance``, dimensional-sin rewire).
* Grand Report v3 §5.5 (Field hierarchy), §32.5 (Field primitive spec).
* ``coding-elegance`` Pattern 4 (illegal states unrepresentable) — a
  residual and a source can no longer be silently mixed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import ANGULAR_RATE_UNITS, Unit
from orpheus.transport.fields._bases import AngularField


__all__ = ["AngularResidual"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class AngularResidual(AngularField):
    r"""L2 per-ordinate transport-balance residual on ``(N, ng, nx, ny)``.

    Parameters
    ----------
    values : NDArray
        Field values of shape ``(N, ng, nx, ny)`` in the principled
        layout (Issue #196 PR-INDEX-5/7).
    space : FunctionSpace
        The function space. Must have ``shape == (mesh.quad.N, mesh.ng,
        mesh.nx, mesh.ny)``. Use :meth:`from_mesh` to derive
        automatically.
    mesh : SNMesh
        The SN phase-space carrier.

    Notes
    -----
    A thin role leaf: all storage, validation, algebra, and the
    ``from_mesh`` / ``from_ndarray`` factories are inherited from
    :class:`~orpheus.transport.fields._bases.AngularField` /
    :class:`~orpheus.transport.fields._bases.BulkField`. This leaf
    carries no residual-specific behaviour beyond its class identity —
    which is exactly what Field's Layer-1 gate uses to keep residual
    arithmetic from silently mixing with
    :class:`~orpheus.transport.sources.angular_source.AngularSource`
    (same shape AND same units) or
    :class:`~orpheus.transport.fields.angular_flux.AngularFlux` (same
    shape, different units). Same-class arithmetic is closed.

    Unlike :class:`AngularFlux` (``zeros_for_sn_mesh``) and
    :class:`AngularSource` (``from_isotropic``), :class:`AngularResidual`
    intentionally exposes NO bespoke factory: a residual is *produced*
    by an operator balance (B.5 ``IterationResidual.from_balance``), not
    constructed from thin air. Tests and producers build it via the
    inherited :meth:`from_mesh`.
    """

    #: The :class:`FunctionSpace` name for this leaf — distinct from
    #: ``"angular_flux"`` / ``"angular_source"`` so the space identity
    #: tracks the role. (Class identity is the primary gate; the space
    #: name keeps the identity legible in diagnostics.) All
    #: mesh/shape/algebra/factory machinery is inherited from
    #: :class:`AngularField` / :class:`BulkField`.
    _SPACE_NAME: ClassVar[str] = "angular_residual"

    #: Dimensional identity (View-G, B.4): per-ordinate rate density
    #: ``1/(cm³·s·sr)`` — :data:`~orpheus.numerics.units.ANGULAR_RATE_UNITS`,
    #: shared with ``AngularSource`` (same units, different role → the gate
    #: is class identity, NOT units). Metadata, not the gate.
    UNITS: ClassVar[Unit] = ANGULAR_RATE_UNITS
