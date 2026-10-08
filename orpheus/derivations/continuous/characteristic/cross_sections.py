r"""The cross sections of each region, as the Galerkin system reads them.

The transport block of a group (:class:`~.assembly.GroupTransport`) carries
the group's total cross section; the system multiplies it by the emission,
which is a matrix per region, indexed ``[to, from]``: the scattering with
the (n,2n) emission, and the fission production. The scattering comes from
the one assembly site of the emission,
:func:`~orpheus.derivations.common.eigenvalue.group_emission`, which the
infinite-medium reference reads too; the fission is kept as its two factors,
:math:`\chi` and :math:`\nu\Sigma_f`, and formed as that function forms it,
so that the adjoint cross sections exchange them.

**The emission support.** Group :math:`g` emits in region :math:`r` iff
some group scatters or fissions into it there: a non-zero row
``scattering[r, g, :]`` or ``fission[r, g, :]``. It is the set of columns
the group's transport block needs, exactly (the user's ruling of
2026-10-07, P1 step (b) fourth rung): a group's total cross section does
not decide it, since a region transparent in :math:`g` can still emit into
:math:`g`, and a region that absorbs in :math:`g` with nothing emitted into
it needs no column.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from functools import cached_property

import numpy as np

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.eigenvalue import group_emission


def _refuse_anisotropy(mixture: Mixture, region: int) -> None:
    r"""Refuse a mixture whose scattering or (n,2n) emission has a non-zero higher Legendre order.

    SCOPE-BOUNDARY[guard] machinery: the anisotropic emission (an angular basis on each line, a moment per order).
    ruling: the user, 2026-10-07, P1 step (b) fourth rung, Q2 (`.claude/plans/characteristic_reference_architecture.md`): the refusal is read here, at the mixtures, the one door.
    revisit: when a consumer poses an anisotropic reference problem.
    """
    for name, stack in (("SigS", mixture.SigS), ("Sig2", mixture.Sig2)):
        for order, block in enumerate(stack[1:], start=1):
            if block.count_nonzero():
                raise NotImplementedError(
                    f"the characteristic reference serves isotropic emission only: region {region}'s "
                    f"{name} has a non-zero Legendre order {order}"
                )


@dataclass(frozen=True, eq=False)
class RegionCrossSections:
    r"""The total cross section and the emission of each region, for :math:`G` groups.

    Built from one mixture per region by :meth:`of`, the way in: the
    :class:`~orpheus.data.macro_xs.mixture.Mixture` is the input boundary
    that checks the data's physics (a balance, a fission spectrum). The
    constructor checks shapes, finiteness and signs only. The arrays are
    read-only.

    The fission emission is stored as its two factors, the spectrum it
    emits into and the production it is driven by, so that the adjoint
    cross sections (:meth:`transposed`) swap them and each keeps its
    meaning: the adjoint's production is the forward spectrum.

    Attributes
    ----------
    total:
        :math:`\Sigma_t`, ``(n, G)``.
    scattering:
        :math:`(\Sigma_s + 2\Sigma_2)^T`, ``(n, G, G)``, indexed ``[region, to, from]``.
    spectrum:
        :math:`\chi`, ``(n, G)``: the groups fission emits into, zero in a region that does not produce.
    production:
        :math:`\nu\Sigma_f`, ``(n, G)``: the neutrons fission emits per unit flux of each group.
    """

    total: np.ndarray
    scattering: np.ndarray
    spectrum: np.ndarray
    production: np.ndarray

    def __post_init__(self) -> None:
        n, groups = np.shape(self.total)
        shapes = {"scattering": (n, groups, groups), "spectrum": (n, groups), "production": (n, groups)}
        for name, shape in shapes.items():
            if np.shape(getattr(self, name)) != shape:
                raise ValueError(
                    f"the {name} of {n} regions in {groups} groups has shape {shape}; got {np.shape(getattr(self, name))}"
                )
        for name in ("total", *shapes):
            array = np.array(getattr(self, name), dtype=float)
            if not np.all(np.isfinite(array)) or array.min(initial=0.0) < 0.0:
                raise ValueError(f"the {name} is finite and non-negative")
            array.flags.writeable = False
            object.__setattr__(self, name, array)

    @classmethod
    def of(cls, mixtures: Sequence[Mixture]) -> "RegionCrossSections":
        """The cross sections of the regions whose mixtures are ``mixtures``, in region order."""
        for region, mixture in enumerate(mixtures):
            _refuse_anisotropy(mixture, region)
        emission = [group_emission(m.SigS[0], m.SigP, m.chi, m.Sig2[0]) for m in mixtures]
        return cls(
            np.stack([m.SigT for m in mixtures]),
            np.stack([e.scattering for e in emission]),
            np.stack([m.chi for m in mixtures]),
            np.stack([m.SigP for m in mixtures]),
        )

    @property
    def n_groups(self) -> int:
        """The number :math:`G` of groups."""
        return self.total.shape[1]

    @cached_property
    def fission(self) -> np.ndarray:
        r""":math:`\chi \otimes \nu\Sigma_f`, ``(n, G, G)``, indexed ``[region, to, from]``: the product :func:`~orpheus.derivations.common.eigenvalue.group_emission` forms, per region."""
        fission = self.spectrum[:, :, None] * self.production[:, None, :]
        fission.flags.writeable = False
        return fission

    def transposed(self) -> "RegionCrossSections":
        r"""The adjoint cross sections: the same total, the scattering with ``[to, from]`` swapped, the spectrum and the production exchanged.

        For isotropic emission the transport is self-adjoint (reciprocity),
        so the adjoint flux of a problem is the flux of the problem posed on
        these cross sections, and its emission support is where something
        is emitted OUT of a group.
        """
        return RegionCrossSections(self.total, self.scattering.swapaxes(1, 2), self.production, self.spectrum)

    def emission_support(self) -> np.ndarray:
        r"""Whether group :math:`g` emits in region :math:`r`, ``(n, G)``: a non-zero transfer into :math:`g` there."""
        return ((self.scattering != 0.0) | (self.fission != 0.0)).any(axis=-1)


__all__ = ["RegionCrossSections"]
