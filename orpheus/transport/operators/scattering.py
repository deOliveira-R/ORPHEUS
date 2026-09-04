r"""The scattering gain :math:`S` — the scattering TERM of the within-group algebra.

This module names the **scattering gain** :math:`S` of the honest
within-group operator algebra :math:`A = L + C - S - N_{2n} - B` (streaming
:math:`L`, collision :math:`C`, scattering gain :math:`S`, the
first-class :math:`(n,2n)` gain :math:`N_{2n}` —
:mod:`orpheus.transport.operators.n2n` — boundary :math:`B`; the
eigenvalue problem is :math:`K = A^{-1}F`). The operator letters are
pinned in ``docs/theory/conventions/notation.rst §notation-symbol-table``;
the multigroup posing in ``docs/theory/methods/sn/slab_multigroup.rst
§sn-scattering-fission-operators``.

:math:`S` is the SCATTERING channel alone — the angular binding of the
mixture's ``SigS`` Legendre stack (yield 1) at the solve's order:

* **P0 isotropic in-scatter** :math:`\Sigma_s^0(g'\to g)\,\phi_{g'}`;
* **Pℓ Galerkin reconstruction** on the frame's basis (real spherical
  harmonics on a sphere rule, Legendre on a 1-D rule) from the
  angular-flux moments, :math:`\ell\ge 1`.

**Every verb lives on the core** (#426 step 2, 2026-09-04 — the F3
ruling: *the kernel tier names the mathematical object, the operator
tier names the TERM*). :class:`ScatteringOperator` is a thin role
subclass of
:class:`~orpheus.transport.operators.transfer.TransferOperator` whose
only content is two class constants — which channel the ONE tier-2 mint
reads (the facade's scattering channel) and the P0 energy binding it
lifts (:class:`~orpheus.transport.operators.isotropic_scattering.IsotropicScattering`)
— and no code.
The body this module carried until 2026-09-04 — faces, arms, the
:math:`R\Lambda M` kernel, the transposes, the fold family — is the
core's, once, shared with the :math:`(n,2n)` term; an AST gate keeps the
role thin. (The **(n,2n) doubling** itself lived here until CS4c step 3
— §14.1 extracted it as the first-class
:class:`~orpheus.transport.operators.n2n.N2NOperator`, because its
bundling — with :math:`S` for anisotropy studies, with :math:`F` for
production accounting — is context-dependent and must not be an
operator-level decision.)

The theory lives in the book — one concept, one home:

* P0 in-scatter, the (n,2n) fold, the Pℓ reconstruction, and the
  :math:`1/W` normalisation chain —
  ``docs/theory/methods/sn/slab_multigroup.rst §mg-inscatter-source``,
  §pn-scatter, §n2n-source, §pn-scatter-rlm.
* The no-prefactor :math:`Y_\ell^m` convention and the Funk–Hecke
  eigenbasis (why :math:`S = R\circ\Lambda\circ M`) —
  ``docs/theory/foundations/spherical_harmonics.rst
  §spherical-harmonics-eigenbasis``.
* The §5.6 integral-kernel reading and the apply-only capability
  surface — ``docs/theory/foundations/operator_algebra.rst
  §integral-kernel-category``, §capability-set-semantics.
* The Euclidean adjoint :math:`S^{T}` (forward fast-path vs adjoint
  frame-form) — ``docs/theory/methods/sn/adjoint.rst
  §sn-scattering-adjoint-source``.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING, ClassVar

from orpheus.transport.material_field import TransferMaterialField
from orpheus.transport.operators.isotropic_scattering import (
    IsotropicScattering,
    IsotropicTransfer,
)
from orpheus.transport.operators.transfer import TransferOperator

if TYPE_CHECKING:
    from orpheus.transport.mesh.material_xs_field import MaterialXSField


__all__ = ["ScatteringOperator"]


class ScatteringOperator(TransferOperator):
    r"""Scattering source operator :math:`S` (P0 + Pℓ) — the scattering role of the transfer binding.

    The :math:`S` instance of
    :class:`~orpheus.transport.operators.transfer.TransferOperator`: the
    angular binding of the scattering channel's field (yield 1) on the
    posed composite at the solve's ``scattering_order``. Build instances
    with the core's :meth:`~TransferOperator.from_solver_data` from a
    :class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField` +
    the posed composite space (the quadrature is reached through the
    space's angular axis — the CS5 generator channel; no ``quadrature``
    field survives); the exact ctor is the core's too. This class is two
    constants: the channel the mint reads and the P0 binding it lifts.

    The **(n,2n) channel is NOT here** (§14.1, ruled 2026-08-30): it is
    the sibling role
    :class:`~orpheus.transport.operators.n2n.N2NOperator`, and the
    within-group algebra spells ``− S − N₂ₙ`` explicitly. Its P0 energy
    binding (:attr:`isotropic_energy`) is
    :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicScattering`
    — the leaf the solver's K_iso sums with the :math:`(n,2n)` term's.
    """

    channel: ClassVar[Callable[["MaterialXSField"], "TransferMaterialField"]] = (
        TransferMaterialField.scattering
    )
    isotropic_binding: ClassVar[type[IsotropicTransfer]] = IsotropicScattering
