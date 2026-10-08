r"""The boundary closure along each line: the period of the unfolded path and its least solution.

The boundary resolvent :math:`P = P_0 + E\,(I - T)^{-1} X` is realised in two
parts, and this module holds both. The line part is :class:`LinePeriod`: on
every wall whose return is specular (a mirror, a partial mirror, vacuum as
amplitude 0) or a periodic wrap, the returned path is a line again,
congruent to the one that left, so :math:`T` is diagonal over lines and its
block on one line is a cycle of at most two traversals. The diffuse walls
couple every line to every other and are the second part, a finite-rank
update over the walls, :class:`WallCoupling`.

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
from math import pi

import numpy as np

from orpheus.geometry.chord import Chord, ConcentricPartition

from .walls import Walls

#: The full solid angle. A line's quadrature weight carries its inverse, so the lines carry every flux as
#: 4 pi times a flux per steradian: one convention, read by the block (:mod:`.assembly`), the walls' injection
#: (:class:`DiffuseWalls`) and the reading (:mod:`~orpheus.derivations.continuous.characteristic.reading`).
FULL_SOLID_ANGLE = 4.0 * pi

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
        Each traversal as a candidate of ``_directed``: transit ``candidate % 2``,
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

    def forward_traversal(self, transit: np.ndarray) -> np.ndarray:
        r"""The traversal of the period that reads transit ``transit`` forward, ``(..., q)`` for ``transit`` ``(..., q)``.

        Every transit a line makes is read forward by its period: the first
        traversal is transit 0 forward, and the successor rule prefers a
        forward candidate, so a line's second transit (through a cavity) is
        reached forward. A point's flux lives on the forward reading; the
        reversed traversals carry the cycle for the opposite line.
        """
        transit = np.asarray(transit)
        forward = (self.candidate[..., None, :] == transit[..., None]) & (self.present & ~self.reversed)[..., None, :]
        if not np.all(forward.any(axis=-1)):
            raise RuntimeError("a transit of a line is not read forward by its period: the period and the chord disagree")
        return np.argmax(forward, axis=-1)

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

    def inflow(self, optical_depth: np.ndarray, outflow: np.ndarray, arriving: np.ndarray | None = None) -> np.ndarray:
        r"""The inflow at each traversal's entry, the least solution of the period's cycle.

        ``optical_depth`` is :math:`\tau_k`, ``(..., 2)``; ``outflow`` is
        :math:`B_k`, the traversal's source integral attenuated to its exit,
        ``(..., 2, *rest)``. Returns ``(..., 2, *rest)``, 0 where absent.

        ``arriving`` is :math:`s_k`, a flux injected at each traversal's
        entry from outside the line part (a diffuse wall's re-entry), shape
        as ``outflow``; it is not multiplied by the wall's amplitude.

        An absent traversal is the cycle's unit (gain 1, nothing returned), so
        one expression serves every rank: the flux entering :math:`k` is
        :math:`e_k = a_{k-1} B_{k-1} + s_k`, and the inflow is
        :math:`(e_k + g_{k-1} e_{k-1}) / (1 - \Pi)` with
        :math:`g = a\,e^{-\tau}` (indices modulo ``_MAX_PERIOD``).
        :math:`1 - \Pi` is formed as ``-expm1(log Pi)``, so a nearly lossless
        line keeps its digits. Raises :class:`TrappedSource` where
        :math:`\Pi = 1` and the entering flux is not zero.
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
        entering = np.roll(returned, 1, axis=axis)
        if arriving is not None:
            entering = entering + np.where(lifted(self.present), np.asarray(arriving, dtype=float), 0.0)
        gain = lifted(np.exp(log_gain))
        around = entering + np.roll(gain, 1, axis=axis) * np.roll(entering, 1, axis=axis)
        numerator = np.where(lifted(self.present), around, 0.0)
        trapped = np.broadcast_to(one_minus_product == 0.0, numerator.shape)
        if np.any(trapped & (numerator != 0.0)):
            raise TrappedSource(
                "a source on a lossless trapped line (every amplitude 1, every optical depth 0) "
                "has no finite inflow"
            )
        return np.divide(numerator, one_minus_product, out=np.zeros(numerator.shape), where=~trapped)


@dataclass(frozen=True, eq=False)
class DiffuseWalls:
    r"""The walls that re-emit isotropically, keyed on a panel partition, and the flux a unit current entering each injects.

    Attributes
    ----------
    breakpoint:
        The breakpoint of each diffuse wall, ``(W,)``.
    amplitude:
        :math:`\alpha`, the diffuse amplitude of each wall, ``(W,)``.
    quarter_area:
        :math:`D_w = A_w/4`, each wall's area over four, ``(W,)``.
    """

    breakpoint: np.ndarray
    amplitude: np.ndarray
    quarter_area: np.ndarray

    def __post_init__(self) -> None:
        w = self.breakpoint.shape
        if len(w) != 1 or self.amplitude.shape != w or self.quarter_area.shape != w:
            raise ValueError(
                f"one breakpoint, amplitude and quarter area per wall; got {self.breakpoint.shape}, "
                f"{self.amplitude.shape}, {self.quarter_area.shape}"
            )

    @classmethod
    def of(cls, walls: Walls, partition: ConcentricPartition) -> "DiffuseWalls":
        r"""The diffuse walls of ``walls``, keyed on ``partition`` (a refinement of the body's), with :math:`D = \pi A / 4\pi`.

        A unit isotropic current entering a wall of area :math:`A` is
        :math:`1/(\pi A)` per steradian, which the lines carry as
        :data:`FULL_SOLID_ANGLE` times that: :math:`1/D`.
        """
        diffuse = [w for w in walls.on(partition).walls if w.diffuse > 0.0]
        at = np.array([w.breakpoint for w in diffuse], dtype=int)
        area = partition.chart.measure_density(np.asarray(partition.breakpoints)[at])
        return cls(at, np.array([w.diffuse for w in diffuse], dtype=float), pi * area / FULL_SOLID_ANGLE)

    def injected(self, period: LinePeriod) -> np.ndarray:
        r"""A unit current entering each wall, as the flux it injects at each traversal's entry, ``(..., 2, W)``.

        The injection :math:`1/D` of :meth:`of`.
        """
        entering = period.present[..., None] & (period.entry_wall[..., None] == self.breakpoint)
        return entering / self.quarter_area


@dataclass(frozen=True, eq=False)
class WallCoupling:
    r"""The diffuse part of the boundary resolvent: the walls that re-emit isotropically, coupled through the line part.

    With :math:`W` diffuse walls, of amplitudes :math:`\alpha`, the emission
    leaves through them as the partial currents :math:`U^{\mathsf T} q`, each
    wall returns :math:`\alpha` of what reaches it, the line part carries a
    unit isotropic current entering at :math:`w` to the fraction
    :math:`T_{w'w}` leaving at :math:`w'`, and a unit current entering at
    :math:`w` produces the flux moments :math:`R_{:,w}`. The returned
    currents solve :math:`j = \alpha(U^{\mathsf T} q + T j)`, so the block's
    update is

    .. math::

        R\,\alpha\,(I - T\alpha)^{-1}\,U^{\mathsf T}.

    Reciprocity makes :math:`R = U D^{-1}` on the emission support, with
    :math:`D = \mathrm{diag}(A_w/4)` and :math:`A_w` the wall's area, so the
    update is the symmetric :math:`U D^{-1}\alpha(I - T\alpha)^{-1}U^{\mathsf T}`
    there; :math:`R` is computed directly because its rows cover every panel
    while :math:`U`'s cover only the emission support.

    **The loss, not the difference.** A current entering at :math:`w` is
    absorbed, leaks through a wall that does not return it, or leaves at a
    diffuse wall: :math:`\sum_{w'} T_{w'w} + \ell_w = 1`, with every term a
    sum of non-negative parts. :math:`I - T\alpha` is formed from the loss
    :math:`\ell` and the off-diagonal transmissions, so its diagonal
    :math:`(1 - \alpha_w) + \alpha_w(\ell_w + \sum_{w' \ne w} T_{w'w})` is
    never a difference of near-equal numbers, and the solve replaces one
    row by the sum of all rows, the balance :math:`(1 - \alpha) + \alpha\ell`, so
    the total current, the mode nearly singular when little is lost, is
    read from the balance (measured 2026-10-06: with only the loss-formed
    diagonal, two white walls still missed conservation by 7.7e-6 at
    :math:`\Sigma_t = 10^{-12}`, through the rounded off-diagonal entries). The subtraction
    :math:`1 - T_{ww}` amplified rounding by the inverse of the absorption
    (measured 2026-10-06 by the elegance review on a closed white sphere:
    conservation off by 1.7e-7 at :math:`\Sigma_t = 10^{-9}`, 1.3e-4 at
    :math:`10^{-12}`), and on a body that absorbs nothing it missed
    singularity by 1.1e-16, returning a flux of 2e16 for a source that has
    none.

    Attributes
    ----------
    response:
        :math:`R`, the flux moments of a unit isotropic current entering at each wall, ``(N, W)``.
    escape:
        :math:`U`, the current leaving at each wall from each emission function, ``(M, W)``.
    transmission:
        :math:`T`, ``(W, W)``, the fraction of a current entering at the column's wall that leaves at the row's.
        Its diagonal is not read by :attr:`returning`: the balance makes it :math:`1 - \ell_w - \sum_{w' \ne w}
        T_{w'w}`, and the diagonal is formed from that, without the subtraction; it is kept for the gates.
    loss:
        :math:`\ell`, ``(W,)``, the fraction of a current entering at each wall that is absorbed or leaks.
    walls:
        The diffuse walls: their breakpoints, their amplitudes :math:`\alpha` and their injection.
    """

    response: np.ndarray
    escape: np.ndarray
    transmission: np.ndarray
    loss: np.ndarray
    walls: DiffuseWalls

    def __post_init__(self) -> None:
        w = self.walls.amplitude.shape[0]
        if (
            self.response.shape[1:] != (w,)
            or self.escape.shape[1:] != (w,)
            or self.transmission.shape != (w, w)
            or self.loss.shape != (w,)
        ):
            raise ValueError(
                f"a coupling of {w} walls has response (N, {w}), escape (M, {w}), transmission ({w}, {w}) and "
                f"loss ({w},); got {self.response.shape}, {self.escape.shape}, {self.transmission.shape}, "
                f"{self.loss.shape}"
            )

    @property
    def returning(self) -> np.ndarray:
        r""":math:`I - T\alpha`, ``(W, W)``, its diagonal formed from the loss."""
        alpha = self.walls.amplitude
        passed_on = self.transmission - np.diag(np.diag(self.transmission))
        kept = (1.0 - alpha) + alpha * (self.loss + passed_on.sum(axis=0))
        return np.diag(kept) - passed_on * alpha[None, :]

    @property
    def currents(self) -> np.ndarray:
        r"""The currents the diffuse walls return per emission function, :math:`\alpha(I - T\alpha)^{-1}U^{\mathsf T}`, ``(W, M)``.

        Where every :math:`\alpha = 1` and every :math:`\ell = 0` the body
        loses nothing and :math:`I - T\alpha` is exactly singular: a source
        reaching the walls raises :class:`TrappedSource` (no finite flux
        exists), and no source (an escape that is zero or empty) returns no
        current.
        """
        alpha = self.walls.amplitude
        if alpha.size and np.all(alpha == 1.0) and np.all(self.loss == 0.0):
            # I - T alpha is exactly singular: a source reaching the walls has no finite flux, and no source has none
            if np.any(self.escape != 0.0):
                raise TrappedSource(
                    "a source in a body that absorbs nothing, behind walls that return everything, has no finite flux"
                )
            return np.zeros((alpha.size, self.escape.shape[0]))
        # The balance row: the column sums of I - T alpha are (1 - alpha) + alpha loss, known without
        # cancellation; the sum of the rows replaces the last, so the total current, the mode that is
        # nearly singular when little is lost, is solved from the balance and not from rounded entries.
        balanced, sources = self.returning, self.escape.T.copy()
        if alpha.size:
            balanced[-1] = (1.0 - alpha) + alpha * self.loss
            sources[-1] = self.escape.T.sum(axis=0)
        return alpha[:, None] * np.linalg.solve(balanced, sources)

    def on_emission(self, emission: np.ndarray, walls: np.ndarray) -> np.ndarray:
        r"""A functional of the stacked sources, read on the emission alone: ``emission`` ``(..., M)`` plus ``walls`` ``(..., W)`` through :attr:`currents`, ``(..., M)``.

        A unit current entering wall :math:`w` contributes ``walls[..., w]``,
        and the emission returns :attr:`currents` of it: the one fold of the
        walls onto the emission, for the block (its diffuse update
        :math:`R\,\alpha\,(I - T\alpha)^{-1}U^{\mathsf T}` is
        ``on_emission(line, response) - line``) and the reading alike.
        """
        return emission + walls @ self.currents


__all__ = ["FULL_SOLID_ANGLE", "DiffuseWalls", "LinePeriod", "TrappedSource", "WallCoupling"]
