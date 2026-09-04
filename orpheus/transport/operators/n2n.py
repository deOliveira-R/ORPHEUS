r"""The :math:`(n,2n)` source operator :math:`N_{2n}` on the angular composite — first-class, anisotropic.

**Why its own operator (CS4c §14.1, ruled 2026-08-30).** The
:math:`(n,2n)` channel is scattering-like (a group transfer carrying its
own anisotropy) AND production-like (it carries the multiplicity
:math:`\nu_{2n}`): its bundling is CONTEXT-dependent — with :math:`S` when
scattering anisotropy is the interesting axis, with :math:`F` when
production accounting is — and a context-dependent bundling must not be
decided at the operator level. So the within-group algebra spells it
explicitly,

.. math::

    A \;=\; L + C - S - N_{2n} - B,

and any bundling is a solver-side :class:`~orpheus.numerics.operator.OperatorSum`
grouping. (Before the extraction the term hid inside
:class:`~orpheus.transport.operators.scattering.ScatteringOperator`'s iso
accumulator — the operator-level commitment this module retired.)

**What it is, algebraically (#426 step 2, 2026-09-04).** The
:math:`N_{2n}` role of
:class:`~orpheus.transport.operators.transfer.TransferOperator`: the
angular binding :math:`N_{2n} = R\,\Lambda_{2n}\,M/W` of the mixture's
``Sig2`` Legendre stack with yield :math:`\nu_{2n} = 2`
(:data:`~orpheus.transport.kernels.N2N_MULTIPLICITY`), built on the SAME
posed space, at the SAME ``scattering_order`` and through the SAME
interned frame as :math:`S` — the :math:`(n,2n)` stack is brought to the
scattering stack's order (a shorter stack is exactly zero above its own
``NL``, the evaluation's statement; ruling O-1/§4.3). The :math:`\ell = 0`
half rides the reaction-rate fast path through the P0 energy binding
:class:`~orpheus.transport.operators.isotropic_transfer.IsotropicN2N`
(:attr:`~orpheus.transport.operators.transfer.TransferOperator.isotropic_energy`
— the K_iso leaf the ray seed's emission sums with :math:`S`'s, ℓ = 0 by
physics); the :math:`\ell \ge 1` half is the frame-conjugated
redistribution, and its transpose is the product chain's reversal
(:attr:`~orpheus.transport.operators.transfer.TransferOperator.full_transfer_kernel`;
``docs/theory/methods/sn/adjoint.rst`` §sn-n2n-adjoint-source). The
production accounting the channel adds — the k-balance's net removal
and the ERR-052 scale anchor — reads the P0 binding's field verb
(:meth:`~orpheus.transport.material_field.TransferMaterialField.add_to_group_rate`),
which carries the yield.

**What it was until 2026-09-04.** ORPHEUS MODELLED the :math:`(n,2n)`
emission as isotropic: this module re-spelled the scattering operator's
arms with ``aniso = None``, minted its frame at :math:`L = 0`, and held a
single-matrix kernel — while the tapes store NL = 7 Legendre moments for
MT=16, the same order as elastic (`[M]` 2026-08-31), and since #426 step
1 ``Mixture.Sig2`` carried all of them. `[M]` 2026-09-03 the discarded
:math:`\ell = 1` moment is worth **−414 Δk·10⁵** on a Be-reflected fast
slab (forward-peaked emission in the reflector sends the emitted pair
outward; 99.9 % of the effect is Be-9's; the ladder converges by
:math:`\ell = 1`) — issue #426, ``docs/theory/methods/sn/adjoint.rst``
§sn-n2n-p0-truncation. That measurement is what retired the twin: the
two terms are now two instances of one binding, and this module is the
role — its channel constant, its P0 binding and its name, no code (an
AST gate, ``tests/transport/test_transfer_roles.py``, keeps it so).

Carrier arms are the core's: composite ``FullField`` (bulk-only; zero
trace), per-ordinate
:class:`~orpheus.transport.fields.angular_flux.AngularFlux`, the windowed
:class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
(whose :math:`\ell=0` moment IS the scalar flux), and the P0-only
:class:`~orpheus.transport.fields.scalar_flux.ScalarFlux` arm in iso
scalar magnitude. A scalar consumer that wants the ENERGY binding alone
builds :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicN2N`
directly.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING, ClassVar

from orpheus.transport.material_field import TransferMaterialField
from orpheus.transport.operators.isotropic_transfer import (
    IsotropicN2N,
    IsotropicTransfer,
)
from orpheus.transport.operators.transfer import TransferOperator

if TYPE_CHECKING:
    from orpheus.transport.mesh.material_xs_field import MaterialXSField

__all__ = ["N2NOperator"]


class N2NOperator(TransferOperator):
    r"""The :math:`(n,2n)` source operator :math:`N_{2n}` — the :math:`(n,2n)` role of the transfer binding.

    The :math:`N_{2n}` instance of
    :class:`~orpheus.transport.operators.transfer.TransferOperator`: the
    angular binding of the :math:`(n,2n)` channel's field (yield 2) on
    the posed composite at the solve's ``scattering_order`` — the same
    faces, arms, kernel and transposes as :math:`S`, over a different
    datum (module docstring). Build instances with the core's
    :meth:`~TransferOperator.from_solver_data` — the SAME interned frame
    :math:`S` is built on, so the two gains share one metric and the
    (n,2n) stack meets the solve at the scattering order (ruling O-1: a
    channel that stores fewer orders is zero above them); the exact ctor
    is the core's too, and the tier-2 equivalence family pins the two
    spellings equal. Its P0
    energy binding (:attr:`isotropic_energy`) is
    :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicN2N`
    — the leaf the solver's K_iso and k-balance read.
    """

    channel: ClassVar[Callable[["MaterialXSField"], "TransferMaterialField"]] = (
        TransferMaterialField.n2n
    )
    isotropic_binding: ClassVar[type[IsotropicTransfer]] = IsotropicN2N
