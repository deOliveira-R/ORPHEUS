r"""The boundary closure along each line: the period of the unfolded path and its least solution.

The boundary resolvent :math:`P = P_0 + E\,(I - T)^{-1} X` is realised in two
parts. This module holds the line part, :class:`LinePeriod`: on every wall
whose return is specular (a mirror, a partial mirror, vacuum as amplitude 0)
or a periodic wrap, the returned path is a line again, congruent to the one
that left, so :math:`T` is diagonal over lines and its block on one line is
a cycle of at most two traversals. The diffuse walls couple every line to
every other and are the second part, a finite-rank update over the walls.

**The period.** Unfold the line through its walls. A *traversal* is one of
the line's transits (:attr:`~orpheus.geometry.chord.Chord.transits`) read
forward or reversed. The traversal after one that exits at wall :math:`w` is
the one entering at the wall's partner (:meth:`~.walls.Walls.partner_at`):
the wall itself under a mirror, the opposite wall under a wrap. Of the two
candidates that enter there, the forward transit is taken; both trace the
same path in the orbit space, because the chart's motion keeps the impact
parameter and the projected speed. The sequence closes on its first
traversal after :math:`m \le 2` steps, and :math:`m` is the line's rank: 1 on
a solid body, on a shell's ray that misses the cavity and on the periodic
slab; 2 on a shell's ray through the cavity and on a slab between mirrors;
0 on a line parallel to the level sets. No case is tagged.

**The closure.** With :math:`B_k` the source integral of traversal
:math:`k` attenuated to its exit, :math:`\tau_k` its optical depth and
:math:`a_k` the specular amplitude of the wall it exits at, the inflow at
each traversal's entry solves the cycle

.. math::

   \psi^{\rm in}_{k+1} = a_k \bigl(e^{-\tau_k}\,\psi^{\rm in}_k + B_k\bigr),

whose least non-negative solution is the Neumann series. It equals the
closed form where the cycle product :math:`\Pi = \prod_k a_k e^{-\tau_k} < 1`
and is exactly 0 on a lossless trapped line (:math:`\Pi = 1`) that carries no
source; a source there has no finite solution and is refused. The amplitude
of the wall traversal :math:`k` exits at multiplies the inflow to
traversal :math:`k + 1`: it is the wall at which the backward path from
:math:`k + 1` reflects, the albedo pairing.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from orpheus.geometry.chord import Chord

from .walls import Walls

#: A slab line has one transit, read in both directions; a radial line's
#: successor is always one of its own (at most two) transits read forward. So a
#: period closes within two traversals.
_MAX_PERIOD = 2


class TrappedSource(ValueError):
    """A source on a lossless trapped line: the least solution is infinite."""


def _directed(chord: Chord) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r"""The four candidate traversals of each line: entry wall, exit wall, presence, ``(..., 4)``.

    Candidates 0 and 1 are the transits 0 and 1 read forward, 2 and 3 the same
    transits read reversed; a traversal is named by its candidate index.
    """
    transits = chord.transits
    entry = np.concatenate([transits.entry_wall, transits.exit_wall], axis=-1)
    exit_ = np.concatenate([transits.exit_wall, transits.entry_wall], axis=-1)
    present = np.concatenate([transits.present, transits.present], axis=-1)
    return entry, exit_, present


@dataclass(frozen=True, eq=False)
class LinePeriod:
    r"""The period of the unfolded path of each line of a chord, in order.

    Every array is ``(..., 2)``, one column per traversal of the period;
    ``present`` masks the line's :attr:`rank` traversals.

    Attributes
    ----------
    chord:
        The chord whose lines are unfolded.
    candidate:
        Each traversal as a candidate of :func:`_directed`: transit ``candidate % 2``,
        read reversed when ``candidate >= 2``.
    amplitude:
        The specular amplitude of the wall the traversal exits at; 0 where absent.
    present:
        Whether the period has the traversal.
    """

    chord: Chord
    candidate: np.ndarray
    amplitude: np.ndarray
    present: np.ndarray

    @classmethod
    def of(cls, chord: Chord, walls: Walls) -> "LinePeriod":
        r"""The period of each line of ``chord`` through ``walls``."""
        entry, exit_, available = _directed(chord)
        has = chord.transits.present[..., 0]

        def successor(current: np.ndarray) -> np.ndarray:
            leaving = np.take_along_axis(exit_, current[..., None], axis=-1)[..., 0]
            target = walls.partner_at(leaving, where=has)
            match = available & (entry == target[..., None])
            if np.any(has & ~match.any(axis=-1)):
                raise RuntimeError("no traversal of a line enters at the partner of the wall it left: the chord's transits and the walls disagree")
            # argmax takes the first match: a forward candidate before a reversed one
            return np.argmax(match, axis=-1)

        first = np.zeros(has.shape, dtype=int)
        second = successor(first)
        rank = np.where(has, np.where(second == first, 1, 2), 0)
        if np.any((rank == 2) & (successor(second) != first)):
            raise RuntimeError("a line's unfolded path did not close within two traversals: the chord's transits and the walls disagree")

        candidate = np.stack([first, second], axis=-1)
        present = np.arange(_MAX_PERIOD) < rank[..., None]
        leaving = np.take_along_axis(exit_, candidate, axis=-1)
        return cls(
            chord=chord,
            candidate=np.where(present, candidate, 0),
            amplitude=walls.specular_at(leaving, where=present).astype(float),
            present=present,
        )

    @property
    def rank(self) -> np.ndarray:
        """The number of traversals in each line's period, ``(...,)``: 0, 1 or 2."""
        return self.present.sum(axis=-1)

    @property
    def transit(self) -> np.ndarray:
        """The index of the transit (into ``chord.transits``) each traversal reads."""
        return self.candidate % 2

    @property
    def reversed(self) -> np.ndarray:
        """Whether each traversal reads its transit backward."""
        return self.present & (self.candidate >= 2)

    def _wall(self, table: np.ndarray) -> np.ndarray:
        walls = np.take_along_axis(table, self.candidate, axis=-1)
        return np.where(self.present, walls, self.chord.partition.no_interface)

    @property
    def entry_wall(self) -> np.ndarray:
        """The breakpoint each traversal enters at; the transits' absent code where absent."""
        return self._wall(_directed(self.chord)[0])

    @property
    def exit_wall(self) -> np.ndarray:
        """The breakpoint each traversal exits at; the transits' absent code where absent."""
        return self._wall(_directed(self.chord)[1])

    def optical_depth(self, sigma_t: np.ndarray) -> np.ndarray:
        r"""The optical depth :math:`\tau_k` of each traversal, ``(..., 2)``; 0 where absent.

        ``sigma_t`` is the total cross section of each region, ``(n,)``; the
        exteriors (codes :math:`n` and :math:`n + 1`) are void. A transit and
        its reverse cross the same slots.
        """
        chord = self.chord
        n = chord.partition.n_regions
        sigma_t = np.asarray(sigma_t, dtype=float)
        if sigma_t.shape != (n,):
            raise ValueError(f"sigma_t holds one total cross section per region, shape ({n},); got {sigma_t.shape}")
        transits = chord.transits
        first = np.take_along_axis(transits.first_slot, self.transit, axis=-1)
        stop = np.take_along_axis(transits.stop_slot, self.transit, axis=-1)
        slot = np.arange(chord.slot_region.shape[-1])
        crossed = self.present[..., None] & (slot >= first[..., None]) & (slot < stop[..., None])
        # an exterior slot inside a transit is untraversed; a void slot adds 0 even at infinite length
        region = chord.slot_region
        sigma = np.where(region < n, sigma_t[np.minimum(region, n - 1)], 0.0)[..., None, :]
        attenuating = crossed & (sigma > 0.0)
        length = np.broadcast_to(chord.slot_length[..., None, :], attenuating.shape)
        depth = np.multiply(length, sigma, out=np.zeros(attenuating.shape), where=attenuating)
        return depth.sum(axis=-1)

    def inflow(self, optical_depth: np.ndarray, outflow: np.ndarray) -> np.ndarray:
        r"""The inflow at each traversal's entry, the least solution of the period's cycle.

        ``optical_depth`` is :math:`\tau_k`, ``(..., 2)``; ``outflow`` is
        :math:`B_k`, the traversal's source integral attenuated to its exit,
        ``(..., 2, *rest)``. Returns ``(..., 2, *rest)``, 0 where absent.

        An absent traversal is the cycle's unit (gain 1, nothing returned), so
        one expression serves every rank: the inflow to :math:`k` is what
        traversal :math:`k - 1` returns, plus what traversal :math:`k - 2`
        returns carried once through :math:`k - 1`, over :math:`1 - \Pi`
        (indices modulo :data:`_MAX_PERIOD`). :math:`1 - \Pi` is formed as
        ``-expm1(log Pi)``, so a nearly lossless line keeps its digits.
        Raises :class:`TrappedSource` where :math:`\Pi = 1` and the outflow
        is not zero.
        """
        outflow = np.asarray(outflow, dtype=float)
        rest = (1,) * (outflow.ndim - self.present.ndim)
        axis = -1 - len(rest)                                    # the period axis

        def lifted(per_traversal: np.ndarray) -> np.ndarray:     # (..., 2) -> (..., 2, *rest)
            return per_traversal.reshape(per_traversal.shape + rest)

        with np.errstate(divide="ignore"):
            log_gain = np.where(self.present, np.log(self.amplitude) - optical_depth, 0.0)
        one_minus_product = -np.expm1(log_gain.sum(axis=-1))
        one_minus_product = one_minus_product.reshape(one_minus_product.shape + (1,) + rest)
        returned = np.where(lifted(self.present), lifted(self.amplitude) * outflow, 0.0)
        gain = lifted(np.exp(log_gain))
        around = np.roll(returned, 1, axis=axis) + np.roll(gain, 1, axis=axis) * np.roll(returned, 2, axis=axis)
        numerator = np.where(lifted(self.present), around, 0.0)
        trapped = np.broadcast_to(one_minus_product == 0.0, numerator.shape)
        if np.any(trapped & (numerator != 0.0)):
            raise TrappedSource(
                "a source on a lossless trapped line (every amplitude 1, every optical depth 0) "
                "has no finite inflow"
            )
        return np.divide(numerator, one_minus_product, out=np.zeros(numerator.shape), where=~trapped)

__all__ = ["LinePeriod", "TrappedSource"]
