r"""The transport along each line on the panel basis: the traversals' source integrals and the Volterra block.

Along one line the angular flux from an emission density :math:`q` is the
Volterra integral :math:`\psi(s) = \int_{-\infty}^{s} e^{-\tau(s', s)}\,q(s')\,\mathrm{d}s'`
plus the inflow the boundary closure returns at each traversal's entry
(:class:`~.closure.LinePeriod`). On the panel basis every quantity the
closure and the Galerkin assembly need is a linear functional of the basis
coefficients, and :class:`TraversalRule` computes them:

* the **outflow** :math:`B_k`, the source integral of each basis function
  over traversal :math:`k` attenuated to its exit, which feeds
  :meth:`~.closure.LinePeriod.inflow`;
* the **entry response** :math:`A_k`, each basis function integrated over
  the traversal attenuated from its entry, the test that an inflow at the
  entry is paired with; it is the outflow of the reversed traversal;
* the **Volterra block** :math:`\sum_L w_L \sum_j \int\!\!\int_{s' < s}
  u_i(s)\,e^{-\tau(s', s)}\,u_j(s')`, the vacuum part of the Galerkin
  matrix, accumulated over a batch of lines with the caller's weights;
* the **angular flux** at parameters on the line, the per-point transport
  the reading shares with the assembly.

**The pieces.** The lines are chorded through the basis's panel partition
(:attr:`~.basis.PanelBasis.partition`), so each slot lies in one panel, and
every slot is cut into pieces by two gradings:

* toward both ends, exponentially, where the attenuation concentrates every
  integrand: piece ends at 1, 2, 4, ..., 64 mean free paths from each end,
  the rest one middle piece. A slot of optical width at most 2 is not cut.
* on a cylinder or a sphere, toward the end nearer the line's closest
  approach, halving, down to that end's distance from the complex branch
  points of the orbit coordinate :math:`c(t) = \sqrt{b^2 + |P\Omega|^2 (t - t^*)^2}`.
  They lie at :math:`t^* \pm i\,b/|P\Omega|`, a distance :math:`c/|P\Omega|`
  from a point of orbit coordinate :math:`c`, so each piece stays a fixed
  number of its widths from them and Gauss-Legendre in arc length converges
  geometrically on every piece (hp grading). On a line through the centre
  (:math:`b = 0`) the coordinate is linear at the closest approach, so the
  two slots meeting there are not graded; a slot starting at an interface
  :math:`r_k > 0` still is.

Each line keeps only its live pieces, in chord order.

**One attenuated integral.** Every quantity is built from one body,
:meth:`TraversalRule._attenuated`: the panel's functions integrated between
two distances along a slot, attenuated to the second, on a rule graded
exponentially toward it. A piece's integral attenuated to its end (or to its
start), carried along the transit, and the integral from a piece's start to a
point inside it, are three uses of it. The integral out of a wide middle
piece lives in its last few mean free paths, which an ungraded rule on the
piece's own nodes misses at O(1) (`[M]` qa's ``probe_thick_psi.py``,
2026-10-06: 0.54 in the flux just past the middle piece of a 1000-mfp slot).

**The Volterra block** pairs the flux at a piece's own nodes with the test
functions there: the flux at a node is what the transit's earlier pieces
carry in, attenuated from the piece's start, plus the integral from the
piece's start to the node. The outer rule integrates the smooth product
:math:`u_i \psi_j`; the inner integral is the graded one.

The design and its rulings: ``.claude/plans/characteristic_reference_architecture.md``
(P1 step (b), second rung).
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np

from orpheus.derivations.common.quadrature import gauss_legendre
from orpheus.geometry.chord import RadialImage
from orpheus.geometry.line import Line

from .basis import PanelBasis
from .closure import LinePeriod
from .walls import Walls

#: Piece ends at 2^k mean free paths, k = 0..6, from a graded end: beyond 64, e^-64 is below double precision.
_DOUBLINGS = 7

#: An interval of optical width at most this is one Gauss panel: e^-2 is integrated to machine precision.
_THIN = 2.0

#: Halvings toward a branch point: depths 2^-k of the slot, k = 1..52 (np.finfo(float).nmant).
_BRANCH_LAYERS = np.finfo(float).nmant


def _toward(stop: np.ndarray, start: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    r"""The ends of the intervals of ``[start, stop]`` graded exponentially toward ``stop``, sorted, ``(..., K + 2)``.

    Depths :math:`2^k/\Sigma` from ``stop``, clipped to the interval; none where
    the interval's optical width is at most :data:`_THIN`.
    """
    width = np.abs(stop - start)
    thick = sigma * width > _THIN
    mean_free_path = 1.0 / np.where(thick, sigma, 1.0)
    depth = np.where(thick[..., None], np.minimum(mean_free_path[..., None] * 2.0 ** np.arange(_DOUBLINGS), width[..., None]), 0.0)
    inward = np.sign(start - stop)[..., None]
    return np.sort(np.concatenate([start[..., None], stop[..., None] + inward * depth, stop[..., None]], axis=-1), axis=-1)


@dataclass(frozen=True, eq=False)
class _Slots:
    r"""Each slot of each line, flattened over the line batch, ``(L, S)``; 0 on a slot no transit traverses.

    Attributes
    ----------
    start:
        The line parameter at the slot's start.
    end:
        The line parameter at the slot's end, the closing crossing's.
    length:
        The slot's 3-D length, cancellation-free.
    panel:
        The slot's panel.
    transit:
        The transit the slot is traversed in, 0 or 1.
    member:
        Whether a transit traverses the slot.
    near_end:
        +1 where the slot's end is the one nearer the line's closest
        approach, -1 where its start is, 0 on a line with no closest
        approach (the slab).
    branch_distance:
        The distance from that end to the orbit coordinate's complex branch
        points, :math:`c_{\rm near}/|P\Omega|`.
    """

    start: np.ndarray
    end: np.ndarray
    length: np.ndarray
    panel: np.ndarray
    transit: np.ndarray
    member: np.ndarray
    near_end: np.ndarray
    branch_distance: np.ndarray


@dataclass(frozen=True, eq=False)
class _Pieces:
    r"""The live pieces of each line in chord order, ``(L, J)``; a dead tail pads the shorter lines.

    Attributes
    ----------
    slot:
        The slot holding the piece (0 on the tail).
    lower, upper:
        The piece's ends, distances from its slot's start (equal on the tail).
    panel, transit:
        The piece's panel and transit (0 on the tail).
    live:
        Whether the piece is one.
    """

    slot: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    panel: np.ndarray
    transit: np.ndarray
    live: np.ndarray


@dataclass(frozen=True, eq=False)
class TraversalRule:
    r"""The quadrature along each line's traversals, on a :class:`~.basis.PanelBasis`.

    Built by :meth:`of`, which chords the lines through the basis's panel
    partition and reads their period through the walls re-keyed onto it, so
    the period, the basis and the cross sections cannot disagree.

    Attributes
    ----------
    period:
        The period of each line, on the panel chord.
    basis:
        The basis.
    sigma:
        The total cross section of each panel, ``(P,)``.
    points:
        Gauss-Legendre points per piece.
    inner_points:
        Gauss-Legendre points per interval of the inner rule (the attenuated integral).
    """

    period: LinePeriod
    basis: PanelBasis
    sigma: np.ndarray
    points: int
    inner_points: int

    def __post_init__(self) -> None:
        chord = self.period.chord
        if chord.partition is not self.basis.partition:
            raise ValueError("the period's lines are chorded through the basis's own panel partition")
        sigma = np.asarray(self.sigma, dtype=float)
        if sigma.shape != (self.basis.n_panels,) or not np.all(np.isfinite(sigma)) or np.any(sigma < 0.0):
            raise ValueError(
                f"sigma holds one finite, non-negative total cross section per panel, shape ({self.basis.n_panels},); "
                f"got {sigma.shape}"
            )
        object.__setattr__(self, "sigma", sigma)
        if self.points < 1 or self.inner_points < 1:
            raise ValueError(f"a rule has at least one point; got {self.points}, {self.inner_points}")
        self._refuse_overflowed_lengths()

    def _refuse_overflowed_lengths(self) -> None:
        r"""Refuse a traversed slot whose 3-D length overflowed.

        **ELEGANCE-DEBT[guard]** #582: retires when the kernel's chord refuses
        (or resolves) a line whose projected speed is subnormal, so that every
        traversed slot it returns has a finite length.
        """
        if not np.all(np.isfinite(np.where(self._slots.member, self._slots.length, 0.0))):
            raise ValueError(
                "a traversed slot has an infinite 3-D length (a line within an underflow of parallel to the "
                "level sets) and no quadrature (#582)"
            )

    @classmethod
    def of(
        cls, lines: Line, basis: PanelBasis, walls: Walls, sigma_t: np.ndarray, points: int, inner_points: int
    ) -> "TraversalRule":
        r"""The rule along ``lines`` through ``basis`` and ``walls``, with ``sigma_t`` the total cross section of each REGION, ``(n,)``."""
        period = LinePeriod.of(basis.partition.chord(lines), walls.on(basis.partition))
        return cls(period, basis, basis.on_panels(sigma_t), points, inner_points)

    # ── the line batch, flattened ────────────────────────────────────────

    @property
    def _batch(self) -> tuple[int, ...]:
        return self.period.present.shape[:-1]

    def _flat(self, per_line: np.ndarray) -> np.ndarray:
        """``(*batch, ...)`` to ``(L, ...)``."""
        return per_line.reshape(-1, *per_line.shape[len(self._batch):])

    def _unflat(self, per_line: np.ndarray) -> np.ndarray:
        """``(L, ...)`` to ``(*batch, ...)``."""
        return per_line.reshape(*self._batch, *per_line.shape[1:])

    # ── the slots ────────────────────────────────────────────────────────

    @cached_property
    def _slots(self) -> _Slots:
        chord = self.period.chord
        transits = chord.transits
        slot = np.arange(chord.slot_region.shape[-1])
        in_range = (slot >= transits.first_slot[..., None]) & (slot < transits.stop_slot[..., None])
        traversed = (chord.slot_length > 0.0) & (chord.slot_region < self.basis.n_panels)
        in_transit = transits.present[..., None] & in_range & traversed[..., None, :]      # (..., 2, S)
        member = in_transit.any(axis=-2)
        crossings = chord.crossings
        start = np.where(member, crossings.parameter[..., :-1], 0.0)
        end = np.where(member, crossings.parameter[..., 1:], 0.0)
        length = np.where(member, chord.slot_length, 0.0)
        panel = np.where(member, chord.slot_region, 0)
        transit = np.argmax(in_transit, axis=-2)
        image = chord.image
        if isinstance(image, RadialImage):
            # the orbit coordinate at each end: the breakpoint's where the crossing is made, b at the closest approach
            r = np.asarray(chord.partition.breakpoints)
            b = image.impact_parameter[..., None]
            at_open = np.where(crossings.present[..., :-1], r[crossings.breakpoint[..., :-1]], b)
            at_close = np.where(crossings.present[..., 1:], r[crossings.breakpoint[..., 1:]], b)
            near_end = np.where(member, np.where(at_close < at_open, 1, -1), 0)
            speed = np.where(image.parallel, 1.0, image.speed)[..., None]
            branch_distance = np.where(member, np.minimum(at_open, at_close) / speed, 0.0)
        else:
            near_end, branch_distance = np.zeros_like(panel), np.zeros_like(length)
        return _Slots(*(self._flat(a) for a in (start, end, length, panel, transit, member, near_end, branch_distance)))

    def _branch_edges(self) -> np.ndarray:
        r"""Piece ends graded toward each slot's end nearer the closest approach, ``(L, S, K)``; 0 where ungraded.

        On a slot ending at the closest approach the branch distance is
        :math:`b/|P\Omega|`; on one starting at the crossing of a small radius
        :math:`r_k` (a small cavity or a small inner region) it is
        :math:`r_k/|P\Omega|`. The pieces halve toward that end until they
        reach it. `[M]` 2026-10-06: one piece at 16 points missed by 2.5e-12
        at :math:`b = 10^{-4}` (the test-architect's ``probe_turn2.py``) and by
        8.5e-7 on a hollow sphere of cavity radius 0.01 (qa's ``probe_branch.py``).

        A change of variable :math:`c = b + (c_{\rm far} - b)u^2` on the slots
        ending at the closest approach was tried first and retired: it leaves
        a branch point in the Jacobian at :math:`u = \pm i\sqrt{2b/(c_{\rm far} - b)}`,
        and it covered only those slots.
        """
        s = self._slots
        depth = s.length[..., None] * 0.5 ** np.arange(1, _BRANCH_LAYERS + 1)
        scale = s.branch_distance[..., None]
        graded = (s.near_end[..., None] != 0) & (scale > 0.0) & (depth > scale)
        at = np.where(s.near_end[..., None] > 0, s.length[..., None] - depth, depth)
        return np.where(graded, at, 0.0)

    # ── the pieces ───────────────────────────────────────────────────────

    @cached_property
    def _pieces(self) -> _Pieces:
        s = self._slots
        sigma = self.sigma[s.panel]
        half = s.length / 2.0
        ends = np.sort(
            np.concatenate(
                [_toward(np.zeros_like(half), half, sigma), _toward(s.length, half, sigma), self._branch_edges()], axis=-1
            ),
            axis=-1,
        )
        lower, upper = ends[..., :-1], ends[..., 1:]                                       # (L, S, M)
        live = (upper > lower) & s.member[..., None]
        slot = np.broadcast_to(np.arange(lower.shape[-2])[:, None], lower.shape)

        def flat(a: np.ndarray) -> np.ndarray:
            return a.reshape(a.shape[0], -1)

        # the live pieces first, each line's in chord order (slot-major, then along the slot)
        count = max(int(flat(live).sum(axis=-1).max(initial=0)), 1)
        order = np.argsort(~flat(live), axis=-1, kind="stable")[:, :count]

        def take(a: np.ndarray) -> np.ndarray:
            return np.take_along_axis(flat(a), order, axis=-1)

        live_j, slot_j = take(live), take(slot)
        return _Pieces(
            slot=np.where(live_j, slot_j, 0),
            lower=np.where(live_j, take(lower), 0.0),
            upper=np.where(live_j, take(upper), 0.0),
            panel=np.where(live_j, np.take_along_axis(s.panel, slot_j, axis=-1), 0),
            transit=np.where(live_j, np.take_along_axis(s.transit, slot_j, axis=-1), 0),
            live=live_j,
        )

    @cached_property
    def _piece_sigma(self) -> np.ndarray:
        return self.sigma[self._pieces.panel]

    @cached_property
    def _piece_depth(self) -> np.ndarray:
        r"""The optical depth of each piece, ``(L, J)``."""
        p = self._pieces
        return self._piece_sigma * (p.upper - p.lower)

    # ── placing nodes and the one attenuated integral ────────────────────

    def _orbit_coordinate(self, t: np.ndarray) -> np.ndarray:
        """The orbit coordinate at parameters ``t`` ``(L, ...)`` of each line."""
        flat = t.reshape(*self._batch, -1)
        return self.period.chord.orbit_coordinate_at(flat).reshape(t.shape)

    def _slot_start(self, slot: np.ndarray) -> np.ndarray:
        """The start parameter of slot ``slot`` ``(L, ...)`` of each line."""
        return self._of_slot(self._slots.start, slot)

    @staticmethod
    def _of_slot(per_slot: np.ndarray, slot: np.ndarray) -> np.ndarray:
        """A per-slot field ``(L, S)`` read at slot ``slot`` ``(L, ...)``."""
        return np.take_along_axis(per_slot, slot.reshape(per_slot.shape[0], -1), axis=-1).reshape(slot.shape)

    def _nodes(self, slot: np.ndarray, lower: np.ndarray, upper: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        r"""``n`` Gauss-Legendre nodes on each ``[lower, upper]`` along slot ``slot``: distance, orbit coordinate, weight.

        ``slot`` broadcasts against ``lower`` and ``upper``, ``(L, ...)``.
        Returns arrays ``(L, ..., n)``.
        """
        rule = gauss_legendre(-1.0, 1.0, n)
        lo, hi = lower[..., None], upper[..., None]
        distance = lo + (hi - lo) * (rule.pts + 1.0) / 2.0
        weight = rule.wts * (hi - lo) / 2.0
        start = self._slot_start(np.broadcast_to(slot, lower.shape))[..., None]
        return distance, self._orbit_coordinate(start + distance), weight

    def _attenuated(self, slot: np.ndarray, start: np.ndarray, stop: np.ndarray) -> np.ndarray:
        r"""The panel's functions integrated between distances ``start`` and ``stop`` along slot ``slot``, attenuated to ``stop``.

        :math:`\int u_m(c(s))\,e^{-\Sigma |{\rm stop} - s|}\,\mathrm{d}s` over the
        interval between the two, on a rule graded exponentially toward
        ``stop``. ``slot`` broadcasts against ``start`` and ``stop``,
        ``(L, ...)``. Returns ``(L, ..., p + 1)``.
        """
        slot = np.broadcast_to(slot, start.shape)
        panel = np.take_along_axis(self._slots.panel, slot.reshape(slot.shape[0], -1), axis=-1).reshape(slot.shape)
        sigma = self.sigma[panel]
        ends = _toward(stop, start, sigma)
        lower, upper = ends[..., :-1], ends[..., 1:]
        used = (upper > lower).reshape(-1, upper.shape[-1]).any(axis=0)
        lower, upper = lower[..., used], upper[..., used]
        distance, orbit, weight = self._nodes(slot[..., None], lower, upper, self.inner_points)
        values = self.basis.values(orbit, np.broadcast_to(panel[..., None, None], orbit.shape))
        factor = weight * np.exp(-sigma[..., None, None] * np.abs(stop[..., None, None] - distance))
        return np.einsum("...kq,...kqi->...i", factor, values)

    @cached_property
    def _to_end(self) -> np.ndarray:
        """Each piece's integral attenuated to its end, ``(L, J, p + 1)``."""
        p = self._pieces
        return self._attenuated(p.slot, p.lower, p.upper)

    @cached_property
    def _to_start(self) -> np.ndarray:
        """Each piece's integral attenuated to its start, ``(L, J, p + 1)``."""
        p = self._pieces
        return self._attenuated(p.slot, p.upper, p.lower)

    # ── along each transit ───────────────────────────────────────────────

    def _transit_depth(self, upstream: bool) -> np.ndarray:
        r"""The optical depth between each piece and its transit's entry (``upstream``) or exit, ``(L, J)``.

        A cumulative sum of the depths of the transit's pieces before (or
        after) it, shifted by one piece: every term is non-negative and no
        difference of depths is formed.
        """
        p = self._pieces
        depth = np.where(p.live, self._piece_depth, 0.0)
        own = np.where(p.transit[..., None] == np.arange(2), depth[..., None], 0.0)    # (L, J, 2)
        none = np.zeros_like(own[:, :1])
        if upstream:
            summed = np.concatenate([none, np.cumsum(own, axis=-2)[:, :-1]], axis=-2)
        else:
            summed = np.concatenate([np.cumsum(own[:, ::-1], axis=-2)[:, ::-1][:, 1:], none], axis=-2)
        return np.take_along_axis(summed, p.transit[..., None], axis=-1)[..., 0]

    def _per_transit(self, transmission: np.ndarray, local: np.ndarray) -> np.ndarray:
        """Each transit's pieces' ``local`` integrals weighted by ``transmission`` and summed, ``(L, 2, N)``."""
        p = self._pieces
        n_lines = p.live.shape[0]
        out = np.zeros((n_lines, 2, self.basis.size))
        line = np.broadcast_to(np.arange(n_lines)[:, None, None], local.shape)
        transit = np.broadcast_to(p.transit[..., None], local.shape)
        np.add.at(out, (line, transit, self.basis.columns(p.panel)), np.where(p.live, transmission, 0.0)[..., None] * local)
        return out

    @cached_property
    def _carried(self) -> np.ndarray:
        r"""What each transit's earlier pieces carry into each piece's start, ``(L, J, N)``.

        A scan: carrying by the transmission of each piece in turn never
        forms a difference of optical depths.
        """
        p = self._pieces
        n_lines, n_pieces = p.live.shape
        carried = np.zeros((n_lines, 2, self.basis.size))
        before = np.zeros((n_lines, n_pieces, self.basis.size))
        lines = np.arange(n_lines)
        columns = self.basis.columns(p.panel)
        transmission = np.exp(-self._piece_depth)
        for j in range(n_pieces):
            k = p.transit[:, j]
            before[:, j] = carried[lines, k]
            step = transmission[:, j, None] * carried[lines, k]
            np.add.at(step, (lines[:, None], columns[:, j]), self._to_end[:, j])
            carried[lines, k] = np.where(p.live[:, j, None], step, carried[lines, k])
        return before

    # ── the functionals ──────────────────────────────────────────────────

    @property
    def optical_depth(self) -> np.ndarray:
        r"""The optical depth :math:`\tau_k` of each traversal of the period, ``(..., 2)``."""
        return self.period.optical_depth(self.sigma)

    @cached_property
    def _exit_integral(self) -> np.ndarray:
        """Each transit read forward: its integral attenuated to its exit, ``(L, 2, N)``."""
        return self._per_transit(np.exp(-self._transit_depth(upstream=False)), self._to_end)

    @cached_property
    def _entry_integral(self) -> np.ndarray:
        """Each transit read forward: its integral attenuated from its entry, ``(L, 2, N)``."""
        return self._per_transit(np.exp(-self._transit_depth(upstream=True)), self._to_start)

    def _per_traversal(self, forward: np.ndarray, backward: np.ndarray) -> np.ndarray:
        """A per-transit quantity read along each traversal of the period: ``forward`` or, reversed, ``backward``."""
        period = self.period
        transit = self._flat(period.transit)[..., None]
        read = np.where(
            self._flat(period.reversed)[..., None],
            np.take_along_axis(backward, transit, axis=-2),
            np.take_along_axis(forward, transit, axis=-2),
        )
        return self._unflat(np.where(self._flat(period.present)[..., None], read, 0.0))

    def outflow(self) -> np.ndarray:
        r"""The source integral :math:`B_k` of each basis function, attenuated to traversal :math:`k`'s exit, ``(..., 2, N)``."""
        return self._per_traversal(self._exit_integral, self._entry_integral)

    def entry_response(self) -> np.ndarray:
        r"""Each basis function integrated over traversal :math:`k`, attenuated from its entry, ``(..., 2, N)``.

        The outflow of the reversed traversal: reading a transit backward
        swaps its exit and its entry.
        """
        return self._per_traversal(self._entry_integral, self._exit_integral)

    def volterra(self, line_weight: np.ndarray) -> np.ndarray:
        r"""The vacuum Galerkin block :math:`\sum_L w_L \int u_i\,\psi^{\rm vac}_j`, ``(N, N)``, over each line's transits read forward.

        ``line_weight`` is :math:`w_L`, broadcast to the line batch. A line's
        flux lives on its own transits in its own direction; the period's
        reversed traversals carry its cycle and belong to the opposite line.
        """
        p = self._pieces
        weight = self._flat(np.broadcast_to(np.asarray(line_weight, dtype=float), self._batch))
        weight = np.where(p.live, weight[:, None], 0.0)                                      # (L, J)
        distance, orbit, node_weight = self._nodes(p.slot, p.lower, p.upper, self.points)   # (L, J, n)
        tests = self.basis.values(orbit, np.broadcast_to(p.panel[..., None], orbit.shape))  # (L, J, n, p+1)
        sigma = self._piece_sigma[..., None]
        # the flux at each node: carried in from the piece's start, plus the piece's own integral up to the node
        entering = node_weight * np.exp(-sigma * (distance - p.lower[..., None]))
        carried = np.einsum("ljq,ljqi,ljn->ljin", entering, tests, self._carried)
        own = self._attenuated(p.slot[..., None], np.broadcast_to(p.lower[..., None], distance.shape), distance)
        triangle = np.einsum("ljq,ljqi,ljqm->ljim", node_weight, tests, own)
        columns = self.basis.columns(p.panel)                                                 # (L, J, p+1)
        block = np.zeros((self.basis.size, self.basis.size))
        np.add.at(block, columns, weight[..., None, None] * carried)
        np.add.at(block, (columns[..., :, None], columns[..., None, :]), weight[..., None, None] * triangle)
        return block

    def angular_flux(self, t: np.ndarray, inflow: np.ndarray) -> np.ndarray:
        r"""The angular flux from each basis function at parameters ``t`` on the line, read forward, ``(..., q, N)``.

        ``t`` is ``(..., q)`` on the caller's lines; ``inflow`` is the period's
        inflow at each traversal's entry per basis function, ``(..., 2, N)``
        (:meth:`~.closure.LinePeriod.inflow`). Each point lies on a transit of
        its line: inside the body, never in a cavity or beyond a wall. The
        flux is the inflow attenuated from the transit's entry, plus what the
        transit's earlier pieces carry in, plus the piece's own integral up to
        the point.
        """
        p = self._pieces
        t = self._flat(np.asarray(t, dtype=float))                                          # (L, q)

        def at(per_piece: np.ndarray, piece: np.ndarray) -> np.ndarray:
            return np.take_along_axis(per_piece, piece, axis=1)

        begins = np.where(p.live, self._slot_start(p.slot) + p.lower, np.inf)               # (L, J)
        piece = np.sum(begins[:, None, :] <= t[..., None], axis=-1) - 1                     # (L, q)
        found = piece >= 0
        piece = np.maximum(piece, 0)
        slot, low = at(p.slot, piece), at(p.lower, piece)
        # a slot ends at its closing crossing: the parameter, not start + length, which round apart by an ulp
        if not np.all(found & at(p.live, piece) & (t <= self._of_slot(self._slots.end, slot))):
            raise ValueError("a point at which the angular flux is read lies on a transit of its line")
        distance = np.minimum(t - self._slot_start(slot), self._of_slot(self._slots.length, slot))
        traversal = self._flat(self.period.forward_traversal(self._unflat(at(p.transit, piece))))
        entering = np.take_along_axis(self._flat(np.asarray(inflow, dtype=float)), traversal[..., None], axis=-2)
        within = at(self._piece_sigma, piece) * (distance - low)
        upstream = at(self._transit_depth(upstream=True), piece)
        carried = np.take_along_axis(self._carried, piece[..., None], axis=1)
        own = np.zeros(t.shape + (self.basis.size,))
        np.put_along_axis(own, self.basis.columns(at(p.panel, piece)), self._attenuated(slot, low, distance), axis=-1)
        psi = np.exp(-(upstream + within))[..., None] * entering + np.exp(-within)[..., None] * carried + own
        return self._unflat(psi)


__all__ = ["TraversalRule"]
