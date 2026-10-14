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

**The speed along a line.** A line's parameter :math:`t` need not be its
arc length :math:`s`: the rule carries the speed :math:`v = \mathrm ds/\mathrm dt`
(:attr:`TraversalRule.speed`), so every optical depth is
:math:`\Sigma v\,\Delta t` and every length element is :math:`v\,\mathrm dt`.
A cylinder's 3-D line at polar angle :math:`\theta` is its in-plane line read
at :math:`v = 1/\sin\theta`, so one in-plane line carries every polar angle
as a family of speeds: its pieces and their node counts are formed once,
and the rules of the attenuation, their nodes and the basis values there
are formed per speed, since the attenuation is each rule's measure. A line parametrised by its own arc length, any line the kernel
chords in three dimensions, has speed 1.

**The pieces.** The lines are chorded through the basis's panel partition
(:attr:`~.basis.PanelBasis.partition`), or a refinement of it at the levels
a flux is read at, so each slot lies in one panel. On a cylinder or a
sphere each slot is cut by the hp walk toward its end nearer the line's
closest approach (:func:`~.grading.branch_walk`): the orbit coordinate
:math:`c(t) = \sqrt{b^2 + |P\Omega|^2 (t - t^*)^2}` has complex branch
points at :math:`t^* \pm i\,b/|P\Omega|`, and every piece is no wider than
its distance to them over the tightness :math:`c = 2`
(:data:`~.grading.TIGHTNESS`), so its integrand is analytic on a fixed
Bernstein ellipse. A slot of an :attr:`~.basis.PanelBasis.even` panel, of
a line through the centre (:math:`b = 0`, where :math:`c` is linear on each
slot) or of a slab is entire and is not cut. Each line keeps only its live
pieces, in chord order. The attenuation cuts nothing: it is the measure of
each piece's rule.

**One integrand, one rule per piece and speed.** Along a piece every
quantity is a smooth factor against the attenuation to one of its ends, so
each is integrated by the Gauss rule of that attenuation
(:class:`~.attenuation_rule.AttenuationRule`) at the piece's optical
half-width, one rule per member of the speed family: positive weights,
exact for the factor's polynomials of degree :math:`2m - 1` at every optical
width, the factor evaluated at the rule's own nodes. The node count
:math:`m` is derived per piece (:func:`~.grading.gauss_nodes`), from the
Bernstein ellipse through the piece's branch point and the basis's bound on
it (:meth:`~.basis.PanelBasis.bound`), and does not depend on the width;
the block's triangle, a product of two panel functions, takes its own
count. An entire piece takes the count its degree in the parameter asks:
:math:`p + 1` (:math:`2p + 1` for the triangle) on an even panel, whose
functions are of degree :math:`2p` in :math:`t`, and
:math:`\lceil (p + 1)/2 \rceil` (:math:`p + 1`) where :math:`c` is linear.

**The nodes are placed from the piece's end.** A node is an offset from
the attenuated end (:class:`~.attenuation_rule.AttenuationRule`), and the
orbit coordinate there is read from the lines anchored at that end: exactly
at a crossing (:meth:`~orpheus.geometry.chord.Chord.rays_at`), and inside a
slot advanced (:meth:`~orpheus.geometry.chord.Rays.advanced`) from the
slot's near end, a crossing or the line's closest approach
(:meth:`~orpheus.geometry.chord.Chord.rays_at_closest_approach`), never
through the line's parameter. Only the live piece ends are read
(:meth:`~orpheus.geometry.chord.Chord.take`).

**The functionals.** The outflow, the entry response and the flux at a
crossing are one-level integrals of each piece (:meth:`_PieceSet.to_end`,
:meth:`_PieceSet.to_start`) carried along the transit. **The Volterra
block** pairs the flux with the test functions: what the transit's earlier
pieces carry into a piece (:attr:`TraversalRule._carried`, a scan of the
pieces' integrals to their ends) is integrated against the tests toward the
piece's start, and the piece's own triangle
:math:`\int\!\!\int_{t' < t} u_i(t)\,e^{-\Sigma v (t - t')}\,u_j(t')` is
the rule of the triangle (:meth:`~.attenuation_rule.AttenuationRule.triangle`).

The design and its rulings: ``.claude/plans/characteristic_reference_architecture.md``
(P1 step (b), second rung; P3 step 3f, product integration along lines).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum
from functools import cached_property

import numpy as np

from orpheus.geometry.chord import Chord, RadialImage, Rays
from orpheus.geometry.line import Line

from .attenuation_rule import AttenuationRule
from .basis import PanelBasis
from .grading import bernstein_parameter, branch_walk, exact_nodes, gauss_nodes
from .closure import LinePeriod
from .walls import Walls


@dataclass(frozen=True, eq=False)
class _Packing:
    r"""Each line's live entries packed first, in order: the restriction of a per-entry field to them, and back.

    The entries are ``(L, O, K)``: per line, :math:`O` owners of :math:`K`
    entries each (a slot's pieces, an integral's intervals). The packing
    holds :math:`M` per line, the most live entries on one line (at least
    1), so an array gathered through it is padded only to the costliest
    line of the batch, not to the batch's costliest owner. On the sphere and
    the slab the lines of a batch are ordered by projected speed, a proxy for
    their cost (:meth:`~.lines.Lines.ordered`), which keeps the remaining
    padding small; the cylinder's in-plane lines all have projected speed 1,
    so the order leaves them as the impact rule made them.

    Attributes
    ----------
    order:
        The flat index of each packed entry among its line's :math:`O K`, ``(L, M)``.
    owner:
        The owner of each packed entry, ``(L, M)``; 0 on the padding.
    live:
        Whether the packed entry is live, ``(L, M)``.
    owners:
        :math:`O`.
    """

    order: np.ndarray
    owner: np.ndarray
    live: np.ndarray
    owners: int

    @classmethod
    def of(cls, live: np.ndarray) -> "_Packing":
        """The packing of ``live`` ``(L, O, K)``, kept in each line's flat (owner-major) order."""
        flat = live.reshape(live.shape[0], -1)
        count = max(int(flat.sum(axis=-1).max(initial=0)), 1)
        order = np.argsort(~flat, axis=-1, kind="stable")[:, :count]
        packed = np.take_along_axis(flat, order, axis=-1)
        return cls(order, np.where(packed, order // live.shape[-1], 0), packed, live.shape[1])

    def gather(self, per_entry: np.ndarray, fill: float | int) -> np.ndarray:
        """A per-entry field ``(L, O, K)`` on the packed entries, ``fill`` on the padding, ``(L, M)``."""
        taken = np.take_along_axis(per_entry.reshape(per_entry.shape[0], -1), self.order, axis=-1)
        return np.where(self.live, taken, fill)

    def scatter(self, per_packed: np.ndarray) -> np.ndarray:
        """Packed values ``(L, M, ...)`` summed onto their owners, ``(L, O, ...)``: the transpose of :meth:`gather`.

        A line's live entries are packed in owner order, so each owner's are
        one run and the sum is a segmented reduction, each run summed in
        order.
        """
        n_lines = per_packed.shape[0]
        out = np.zeros((n_lines * self.owners, *per_packed.shape[2:]))
        target = (np.arange(n_lines)[:, None] * self.owners + self.owner)[self.live]           # non-decreasing
        if target.size:
            runs = np.flatnonzero(np.diff(target, prepend=-1))
            out[target[runs]] = np.add.reduceat(per_packed[self.live], runs, axis=0)
        return out.reshape(n_lines, self.owners, *per_packed.shape[2:])


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
    height:
        :math:`b/|P\Omega|`, the distance in the parameter from the line to
        the orbit coordinate's complex branch points
        :math:`t^* \pm i\,b/|P\Omega|`; 0 on the slab.
    approach:
        The parameter distance from the closest approach :math:`t^*` to the
        slot's near end, :math:`y_k/|P\Omega|` from the half-chord
        :math:`y_k` at the crossing there (0 where the near end is
        :math:`t^*` itself); 0 on the slab.
    anchor:
        The column of :attr:`~orpheus.geometry.chord.Chord.crossings` at the
        slot's anchor, the end every piece end is measured from: its near end,
        its start on the slab.
    anchored_on_crossing:
        Whether the anchor is a crossing the line makes; where not, it is the
        closest approach.
    """

    start: np.ndarray
    end: np.ndarray
    length: np.ndarray
    panel: np.ndarray
    transit: np.ndarray
    member: np.ndarray
    near_end: np.ndarray
    height: np.ndarray
    approach: np.ndarray
    anchor: np.ndarray
    anchored_on_crossing: np.ndarray

    @property
    def direction(self) -> np.ndarray:
        r"""The slot's direction along the line from its anchor, ``(L, S)``: -1 where the anchor is its end, else +1."""
        return np.where(self.anchor > np.arange(self.start.shape[-1]), -1.0, 1.0)

    @property
    def origin(self) -> np.ndarray:
        r"""The parameter of each slot's anchor."""
        return np.where(self.direction < 0, self.end, self.start)

    @property
    def first(self) -> np.ndarray:
        r"""The slot's start, as a signed distance from its anchor: :math:`-\ell` where the anchor is its end, else 0."""
        return np.where(self.direction < 0, -self.length, 0.0)

    @property
    def last(self) -> np.ndarray:
        r"""The slot's end, as a signed distance from its anchor: 0 where the anchor is its end, else :math:`\ell`."""
        return np.where(self.direction < 0, 0.0, self.length)


@dataclass(frozen=True, eq=False)
class _Pieces:
    r"""The live pieces of each line in chord order, ``(L, J)``; a dead tail pads the shorter lines.

    Attributes
    ----------
    slot:
        The slot holding the piece (0 on the tail).
    lower, upper:
        The piece's ends, signed distances along the line from its slot's anchor (:attr:`_Slots.origin`):
        the hp walk's own distances, never a difference with the slot's length (equal on the tail).
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


class _End(Enum):
    r"""A piece's end a one-level integral is attenuated to: the :class:`_Pieces` field holding it, and the direction its nodes lie in along the line."""

    UPPER = ("upper", -1.0)
    LOWER = ("lower", 1.0)

    @property
    def toward(self) -> float:
        r"""The direction along the line from the end into its piece: behind the upper end, ahead of the lower."""
        return self.value[1]

    def of(self, pieces: _Pieces) -> np.ndarray:
        r"""The end's signed distance from each piece's slot anchor, ``(L, J)``."""
        return getattr(pieces, self.value[0])


@dataclass(frozen=True, eq=False)
class _EndPoints:
    r"""The live pieces' ends on their lines, each read from where it is exact; a padded entry is no point.

    Attributes
    ----------
    groups:
        Each group of ends read one way (from a crossing, or from the closest approach): its points
        (:class:`~orpheus.geometry.chord.Rays`, ``(P,)``) and their flat places among the pieces ``(L, J)``.
    """

    groups: tuple[tuple[Rays, np.ndarray], ...]

    @cached_property
    def at(self) -> np.ndarray:
        r"""The flat places of every live end among the pieces, ``(P,)``, in the order :meth:`orbit_coordinate_at` returns them."""
        return np.concatenate([at for _, at in self.groups]) if self.groups else np.zeros(0, dtype=int)

    def orbit_coordinate_at(self, distance: np.ndarray) -> np.ndarray:
        r"""The orbit coordinate at the distances ``distance`` ``(P, q)`` ahead of each live end (:attr:`at`), ``(P, q)``."""
        if not self.groups:
            return np.zeros_like(distance)
        sizes = np.cumsum([at.size for _, at in self.groups])[:-1]
        return np.concatenate([points.orbit_coordinate_at(part) for (points, _), part in zip(self.groups, np.split(distance, sizes))])


@dataclass(frozen=True, eq=False)
class _PieceSet:
    r"""A rule's slots cut into pieces by the hp walk, and the quadrature on them.

    Each slot is cut toward its end nearer the line's closest approach
    (:func:`~.grading.branch_walk`); an entire slot is one piece. Every
    integral over a piece is a smooth factor against the attenuation to one
    of its ends, integrated by the Gauss rule of that attenuation
    (:class:`~.attenuation_rule.AttenuationRule`), per member of the speed
    family, with the node counts of :attr:`node_counts`.
    """

    rule: "TraversalRule"
    pieces: _Pieces

    @classmethod
    def of(cls, rule: "TraversalRule") -> "_PieceSet":
        r"""The pieces of ``rule``'s slots."""
        s = rule._slots
        cut = branch_walk(np.where(rule._walked, s.length, 0.0), s.height, s.approach)          # from the near end, the anchor
        ends = np.sort(np.concatenate([s.first[..., None], s.last[..., None], s.direction[..., None] * cut], axis=-1), axis=-1)
        lower, upper = ends[..., :-1], ends[..., 1:]                                       # (L, S, M)
        live = (upper > lower) & s.member[..., None]
        # the live pieces first, each line's in chord order (slot-major, then along the slot)
        packing = _Packing.of(live)
        slot = packing.owner
        return cls(
            rule,
            _Pieces(
                slot=slot,
                lower=packing.gather(lower, 0.0),
                upper=packing.gather(upper, 0.0),
                panel=np.where(packing.live, rule._of_slot(s.panel, slot), 0),
                transit=np.where(packing.live, rule._of_slot(s.transit, slot), 0),
                live=packing.live,
            ),
        )

    @cached_property
    def sigma(self) -> np.ndarray:
        r"""The total cross section of each piece, per unit arc length, ``(L, J)``."""
        return self.rule.sigma[self.pieces.panel]

    @cached_property
    def depth(self) -> np.ndarray:
        r"""The optical depth of each piece, per family member, ``(F, L, J)``."""
        p = self.pieces
        return self.rule._speed[..., None] * self.sigma * (p.upper - p.lower)

    @cached_property
    def begins(self) -> np.ndarray:
        r"""The line parameter at each piece's beginning, ``(L, J)``; :math:`\infty` on the padding.

        Its slot's start crossing where it begins there, exactly; else its
        slot's anchor plus its signed distance.
        """
        rule, p = self.rule, self.pieces
        s = rule._slots
        at_start = p.lower == rule._of_slot(s.first, p.slot)
        begins = np.where(at_start, rule._slot_start(p.slot), rule._of_slot(s.origin, p.slot) + p.lower)
        return np.where(p.live, begins, np.inf)

    @cached_property
    def half_length(self) -> np.ndarray:
        r""":math:`v h/2`, each piece's arc length over its reference coordinate's :math:`[-1, 1]`, per family member, ``(F, L, J)``."""
        p = self.pieces
        return self.rule._speed[..., None] * (p.upper - p.lower) / 2.0

    @cached_property
    def node_counts(self) -> tuple[np.ndarray, np.ndarray]:
        r"""The nodes of each piece's rules, ``(L, J)`` each: for one panel function, and for the triangle's product of two.

        A piece cut by the hp walk takes the count of :func:`~.grading.gauss_nodes`:
        its branch points :math:`t^* \pm i\,b/|P\Omega|` lie at
        :math:`z_0 = 2(d + i\,b/|P\Omega|)/h` in its coordinate, :math:`d` the
        distance from :math:`t^*` to its middle and :math:`h` its width, and on
        the Bernstein ellipse :math:`\rho` inside them the orbit coordinate is
        at most :math:`|P\Omega|\,|(b/|P\Omega|,\; d + (h/2)(\rho + 1/\rho)/2)|`,
        whose panel functions the basis bounds (:meth:`~.basis.PanelBasis.bound`;
        squared for the triangle). An entire piece takes its degree's count.
        """
        rule, p = self.rule, self.pieces
        s, degree = rule._slots, rule.basis.degree
        even = rule.basis.even[p.panel]
        walked = p.live & rule._of_slot(rule._walked, p.slot)
        # an even panel's functions are of degree 2p in the parameter; elsewhere the coordinate is linear, degree p
        one = np.where(even, exact_nodes(2 * degree), exact_nodes(degree))
        two = np.where(even, exact_nodes(4 * degree), exact_nodes(2 * degree))
        if not walked.any():
            return one, two
        height, approach, speed = (rule._of_slot(a, p.slot)[walked] for a in (
            s.height, s.approach, rule._projected_speed[:, None] * np.ones_like(s.length)))
        lower, upper, panel = p.lower[walked], p.upper[walked], p.panel[walked]
        width = upper - lower
        middle = (lower + upper) / 2.0
        distance = approach + np.abs(middle)                                                  # from t* to the piece's middle
        singular = bernstein_parameter(2.0 * (distance + 1j * height) / width)

        def majorant(rho: np.ndarray) -> np.ndarray:
            reach = distance + width / 2.0 * (rho + 1.0 / rho) / 2.0
            return rule.basis.bound(speed * np.hypot(height, reach), panel)

        one, two = one.copy(), two.copy()
        one[walked] = gauss_nodes(singular, majorant)
        two[walked] = gauss_nodes(singular, lambda rho: majorant(rho) ** 2)
        return one, two

    @cached_property
    def _upper(self) -> _EndPoints:
        r"""Each live piece's upper end on its line (:meth:`_end_points`)."""
        return self._end_points(_End.UPPER)

    @cached_property
    def _lower(self) -> _EndPoints:
        r"""Each live piece's lower end on its line (:meth:`_end_points`)."""
        return self._end_points(_End.LOWER)

    def _end_points(self, end: _End) -> _EndPoints:
        r"""Each live piece's ``end`` on its line, read from where it is exact (:class:`_EndPoints`).

        A piece end on one of its slot's crossings is read from that crossing
        (:meth:`~orpheus.geometry.chord.Chord.rays_at`). Every other end is
        advanced (:meth:`~orpheus.geometry.chord.Rays.advanced`) by its own
        signed distance from its slot's anchor, the end the hp walk measured
        from: the crossing there, or the closest approach itself
        (:meth:`~orpheus.geometry.chord.Chord.rays_at_closest_approach`) where
        the slot ends without one, so a point near :math:`t^*` is never
        advanced from a far crossing, which forms :math:`r^2 - b^2` by
        cancellation. Only the live ends are read, each on the chord of its
        own line (:meth:`~orpheus.geometry.chord.Chord.take`): a line with no
        piece places no point.
        """
        rule, p = self.rule, self.pieces
        s, chord = rule._slots, rule.period.chord
        present = rule._flat(chord.crossings.present)                                      # (L, S + 1)
        at = end.of(p)
        first, last, anchor, anchored = (rule._of_slot(a, p.slot) for a in (s.first, s.last, s.anchor, s.anchored_on_crossing))
        at_start = (at == first) & np.take_along_axis(present, p.slot, axis=-1)
        at_end = (at == last) & np.take_along_axis(present, p.slot + 1, axis=-1)
        column = np.where(at_start, p.slot, np.where(at_end, p.slot + 1, anchor))
        advance = np.where(at_start | at_end, 0.0, at)
        on_crossing = p.live & (at_start | at_end | anchored)
        line = np.broadcast_to(np.arange(p.live.shape[0])[:, None], p.live.shape)

        def from_crossing(chosen: np.ndarray) -> Rays:
            return chord.take(line[chosen]).rays_at(column[chosen])

        def from_closest_approach(chosen: np.ndarray) -> Rays:
            return chord.take(line[chosen]).rays_at_closest_approach()

        groups = tuple(
            (read(chosen).advanced(advance[chosen]), np.flatnonzero(chosen))
            for chosen, read in ((on_crossing, from_crossing), (p.live & ~on_crossing, from_closest_approach))
            if chosen.any()
        )
        return _EndPoints(groups)

    def _values(self, end: _End, offset: np.ndarray) -> np.ndarray:
        r"""The panel functions at reference offsets ``offset`` ``(F, L, J, *q)`` from each piece's ``end``, ``(F, L, J, *q, p + 1)``; 0 on the padding.

        An offset :math:`\delta` is the distance :math:`\delta h/2` from the
        end toward the piece's other end (:attr:`_End.toward`). The basis is
        evaluated at the live ends' nodes only.
        """
        p = self.pieces
        members, nodes = offset.shape[0], offset.shape[3:]
        ell = end.toward * offset * ((p.upper - p.lower) / 2.0).reshape(1, *p.live.shape, *(1,) * len(nodes))
        ahead = np.moveaxis(ell, 0, 2).reshape(p.live.size, -1)                             # (L J, F q)
        points = self._upper if end is _End.UPPER else self._lower
        at = points.at
        orbit = points.orbit_coordinate_at(ahead[at])                                         # (P, F q)
        live = self.rule.basis.values(orbit, np.broadcast_to(p.panel.reshape(-1)[at][:, None], orbit.shape))
        values = np.zeros((p.live.size, *orbit.shape[1:], self.rule.basis.per_panel))
        values[at] = live
        values = values.reshape(*p.live.shape, members, *nodes, self.rule.basis.per_panel)
        return np.moveaxis(values, 2, 0)

    def _one_level(self, end: _End) -> np.ndarray:
        r"""Each piece's integral attenuated to its ``end``, ``(F, L, J, p + 1)``."""
        rule = AttenuationRule.of(self.depth / 2.0, self.node_counts[0])
        values = self._values(end, rule.offset)
        return self.half_length[..., None] * np.einsum("fljq,fljqi->flji", rule.weight, values)

    @cached_property
    def to_end(self) -> np.ndarray:
        r"""Each piece's integral attenuated to its end, ``(F, L, J, p + 1)``."""
        return self._one_level(_End.UPPER)

    @cached_property
    def to_start(self) -> np.ndarray:
        r"""Each piece's integral attenuated to its start, ``(F, L, J, p + 1)``."""
        return self._one_level(_End.LOWER)

    @cached_property
    def triangle(self) -> np.ndarray:
        r"""Each piece's Volterra triangle :math:`\int\!\!\int_{t' < t} u_a(t)\,e^{-\Sigma v (t - t')}\,u_b(t')\,v^2\,\mathrm dt'\,\mathrm dt`, ``(F, L, J, p + 1, p + 1)``.

        The test :math:`u_a` at each node's lead point, the source :math:`u_b`
        at its lag point upstream (:class:`~.attenuation_rule.TriangleRule`),
        both read from the piece's end.
        """
        rule = AttenuationRule.of(self.depth / 2.0, self.node_counts[1]).triangle()
        test, source = self._values(_End.UPPER, rule.lead), self._values(_End.UPPER, rule.lag)
        paired = np.einsum("fljqk,fljqka,fljqkb->fljab", rule.weight, test, source)
        return self.half_length[..., None, None] ** 2 * paired

    def transit_depth(self, upstream: bool) -> np.ndarray:
        r"""The optical depth between each piece and its transit's entry (``upstream``) or exit, ``(F, L, J)``.

        A cumulative sum of the depths of the transit's pieces before (or
        after) it, shifted by one piece: every term is non-negative and no
        difference of depths is formed.
        """
        p = self.pieces
        depth = np.where(p.live, self.depth, 0.0)
        own = np.where(p.transit[..., None] == np.arange(2), depth[..., None], 0.0)    # (F, L, J, 2)
        none = np.zeros_like(own[..., :1, :])
        if upstream:
            summed = np.concatenate([none, np.cumsum(own, axis=-2)[..., :-1, :]], axis=-2)
        else:
            summed = np.concatenate([np.cumsum(own[..., ::-1, :], axis=-2)[..., ::-1, :][..., 1:, :], none], axis=-2)
        return np.take_along_axis(summed, p.transit[None, ..., None], axis=-1)[..., 0]

    def per_transit(self, transmission: np.ndarray, local: np.ndarray) -> np.ndarray:
        """Each transit's pieces' ``local`` integrals weighted by ``transmission`` and summed, ``(F, L, 2, N)``."""
        p = self.pieces
        in_transit = (p.transit[:, None, :] == np.arange(2)[:, None]) & p.live[:, None, :]      # (L, 2, J)
        return self.on_columns(np.where(in_transit, transmission[:, :, None, :], 0.0)[..., None] * local[:, :, None])

    @cached_property
    def _on_panel(self) -> np.ndarray:
        r"""Whether each live piece lies in each basis panel, ``(L, J, P)``, as 0 or 1."""
        p = self.pieces
        return ((p.panel[..., None] == np.arange(self.rule.basis.n_panels)) & p.live[..., None]).astype(float)

    def on_columns(self, per_piece: np.ndarray) -> np.ndarray:
        r"""Integrals per piece and panel function, ``(F, L, ..., J, p + 1)``, summed onto the basis columns of each piece's panel, ``(F, L, ..., N)``."""
        indicator = self._on_panel.reshape(self._on_panel.shape[:1] + (1,) * (per_piece.ndim - 4) + self._on_panel.shape[1:])
        per_panel = np.swapaxes(per_piece, -1, -2) @ indicator                                  # (F, L, ..., p + 1, P)
        return np.swapaxes(per_panel, -1, -2).reshape(*per_panel.shape[:-2], self.rule.basis.size)


@dataclass(frozen=True, eq=False)
class TraversalRule:
    r"""The quadrature along each line's traversals, on a :class:`~.basis.PanelBasis`.

    Built by :meth:`of`, which chords the lines through the basis's panel
    partition and reads their period through the walls re-keyed onto it, so
    the period, the basis and the cross sections cannot disagree; or by
    :meth:`on_chord`, from a chord through a refinement of that partition,
    whose extra levels are where the flux is read (:meth:`at`).

    Attributes
    ----------
    period:
        The period of each line, on a chord through the basis's panel partition or a refinement of it.
    basis:
        The basis.
    sigma:
        The total cross section of each panel, ``(P,)``.
    speed:
        :math:`v = \mathrm ds/\mathrm dt`, the arc length per unit line parameter, positive: a scalar, every line
        read at it, or a family ``(F,)``, every line read at each member, which every functional returns first,
        ``(F, *batch, ...)``. Never one per line: a family and a per-line field of equal length would be told apart
        by nothing. 1 on lines parametrised by arc length. The reciprocal of the kernel's projected speed
        :math:`|P\Omega| = \mathrm dt/\mathrm ds` on a lifted line (:func:`~.lines.axial_lift`).
    """

    period: LinePeriod
    basis: PanelBasis
    sigma: np.ndarray
    speed: np.ndarray | float = 1.0

    def __post_init__(self) -> None:
        self._panel_of_region                                                              # a refinement of the panels, or refused
        sigma = np.asarray(self.sigma, dtype=float)
        if sigma.shape != (self.basis.n_panels,) or not np.all(np.isfinite(sigma)) or np.any(sigma < 0.0):
            raise ValueError(
                f"sigma holds one finite, non-negative total cross section per panel, shape ({self.basis.n_panels},); "
                f"got {sigma.shape}"
            )
        object.__setattr__(self, "sigma", sigma)
        speed = np.asarray(self.speed, dtype=float)
        if speed.ndim > 1 or not np.all(np.isfinite(speed)) or np.any(speed <= 0.0):
            raise ValueError(f"a speed is finite and positive, a scalar or a family (F,) shared by every line; got shape {speed.shape}")
        object.__setattr__(self, "speed", speed)
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
        cls,
        lines: Line,
        basis: PanelBasis,
        walls: Walls,
        sigma_t: np.ndarray,
        level: tuple[np.ndarray, np.ndarray] | None = None,
        speed: np.ndarray | float = 1.0,
    ) -> "TraversalRule":
        r"""The rule along ``lines`` through ``basis`` and ``walls``, with ``sigma_t`` the total cross section of each REGION, ``(n,)``.

        ``level`` is each line's exact level, a radius and its half-chord
        there (:meth:`~orpheus.geometry.chord.ConcentricPartition.chord`);
        ``speed`` is :attr:`speed`.
        """
        return cls.on_chord(basis.partition.chord(lines, level), basis, walls, sigma_t, speed)

    @classmethod
    def on_chord(
        cls,
        chord: Chord,
        basis: PanelBasis,
        walls: Walls,
        sigma_t: np.ndarray,
        speed: np.ndarray | float = 1.0,
    ) -> "TraversalRule":
        r"""The rule along the lines of ``chord``, a chord of ``basis``'s panel partition (:meth:`of`, for a chord already formed, such as :meth:`~orpheus.geometry.chord.ConcentricPartition.rays_from`'s)."""
        period = LinePeriod.of(chord, walls.on(chord.partition))
        return cls(period, basis, basis.on_panels(sigma_t), speed)

    # ── the line batch, flattened, and the family of speeds over it ──────

    @property
    def _lines(self) -> tuple[int, ...]:
        """The line batch's shape."""
        return self.period.present.shape[:-1]

    @property
    def _family(self) -> tuple[int, ...]:
        """The family's shape, leading every functional: :attr:`speed`'s, ``()`` or ``(F,)``."""
        return np.shape(self.speed)

    @cached_property
    def _speed(self) -> np.ndarray:
        """The speed of each member of the family on each line, ``(F, L)``; one member for a scalar speed."""
        speed = np.reshape(self.speed, (-1, 1))
        return np.broadcast_to(speed, (speed.shape[0], int(np.prod(self._lines))))

    def _flat(self, per_line: np.ndarray) -> np.ndarray:
        """``(*batch, ...)`` to ``(L, ...)``; the line count is spelled, so a field with no entry per line flattens too."""
        return per_line.reshape(int(np.prod(self._lines)), *per_line.shape[len(self._lines):])

    def _unflat(self, per_line: np.ndarray) -> np.ndarray:
        """``(L, ...)`` to ``(*batch, ...)``."""
        return per_line.reshape(self._lines + per_line.shape[1:])

    def _unflat_family(self, per_member: np.ndarray) -> np.ndarray:
        """``(F, L, ...)`` to ``(*family, *batch, ...)``."""
        return per_member.reshape(self._family + self._lines + per_member.shape[2:])

    # ── the slots ────────────────────────────────────────────────────────

    @cached_property
    def _panel_of_region(self) -> np.ndarray:
        r"""The basis panel each region of the chord's partition lies in, ``(m,)``: the identity when the two are one partition.

        The chord's partition coarsens onto the panels, or the rule is refused
        (:meth:`~orpheus.geometry.chord.ConcentricPartition.coarsening`).
        """
        return self.period.chord.partition.coarsening(self.basis.partition)

    @cached_property
    def _slots(self) -> _Slots:
        chord = self.period.chord
        transits = chord.transits
        slot = np.arange(chord.slot_region.shape[-1])
        in_range = (slot >= transits.first_slot[..., None]) & (slot < transits.stop_slot[..., None])
        traversed = (chord.slot_length > 0.0) & (chord.slot_region < chord.partition.n_regions)
        in_transit = transits.present[..., None] & in_range & traversed[..., None, :]      # (..., 2, S)
        member = in_transit.any(axis=-2)
        crossings = chord.crossings
        start = np.where(member, crossings.parameter[..., :-1], 0.0)
        end = np.where(member, crossings.parameter[..., 1:], 0.0)
        length = np.where(member, chord.slot_length, 0.0)
        panel = np.where(member, self._panel_of_region[np.where(member, chord.slot_region, 0)], 0)
        transit = np.argmax(in_transit, axis=-2)
        image = chord.image
        if isinstance(image, RadialImage):
            # the orbit coordinate at each end: the breakpoint's where the crossing is made, b at the closest approach
            r = np.asarray(chord.partition.breakpoints)
            b = image.impact_parameter[..., None]
            at_open = np.where(crossings.present[..., :-1], r[crossings.breakpoint[..., :-1]], b)
            at_close = np.where(crossings.present[..., 1:], r[crossings.breakpoint[..., 1:]], b)
            near_end = np.where(member, np.where(at_close < at_open, 1, -1), 0)
            # the near end's distance from t*: its crossing's half-chord from the exact level over |P Omega|, 0 at t* itself
            anchor = np.where(near_end > 0, 1, 0) + np.arange(start.shape[-1])
            near_present = np.take_along_axis(crossings.present, anchor, axis=-1)
            near_radius = np.where(near_present, r[np.take_along_axis(crossings.breakpoint, anchor, axis=-1)], 0.0)
            speed = self._unflat(self._projected_speed)[..., None]
            approach = np.where(member & near_present, image.half_chord_at(near_radius) / speed, 0.0)
            height = np.where(member, b / speed, 0.0)
        else:
            near_end, height, approach = np.zeros_like(panel), np.zeros_like(length), np.zeros_like(length)
            anchor = np.broadcast_to(np.arange(start.shape[-1]), start.shape)
            near_present = np.take_along_axis(crossings.present, anchor, axis=-1)
        fields = (start, end, length, panel, transit, member, near_end, height, approach, anchor, near_present)
        return _Slots(*(self._flat(np.asarray(a)) for a in fields))

    @cached_property
    def _projected_speed(self) -> np.ndarray:
        r""":math:`|P\Omega|` of each line, ``(L,)``: 1 on a line parallel to the level sets, and on the slab."""
        image = self.period.chord.image
        if isinstance(image, RadialImage):
            return self._flat(np.where(image.parallel, 1.0, image.speed))
        return np.ones(int(np.prod(self._lines)))

    @cached_property
    def _walked(self) -> np.ndarray:
        r"""Whether each slot is cut by the hp walk, ``(L, S)``: a traversed slot of a line with a branch point off it (:math:`b > 0`), on a panel that is not :attr:`~.basis.PanelBasis.even`."""
        s = self._slots
        return s.member & (s.height > 0.0) & (s.near_end != 0) & ~self.basis.even[s.panel]

    # ── the pieces ───────────────────────────────────────────────────────

    @cached_property
    def _pieces(self) -> _PieceSet:
        r"""The pieces every functional reads (:class:`_PieceSet`)."""
        return _PieceSet.of(self)

    @property
    def piece_counts(self) -> np.ndarray:
        r"""The number of live pieces on each line, ``(*batch,)``: a line's share of :attr:`extent`, before the padding."""
        return self._unflat(self._pieces.pieces.live.sum(axis=-1))

    @property
    def extent(self) -> int:
        r"""The number of piece slots the rule's arrays hold at most: family members times lines times the most pieces on one line (:attr:`piece_counts`).

        Every node array scales with it, so a caller bounds a batch's memory by
        it; a rule of one line holds its whole family, whatever its extent.
        """
        return int(self._pieces.pieces.live.size) * self._speed.shape[0]

    # ── reading per-slot fields at pieces ────────────────────────────────

    def _slot_start(self, slot: np.ndarray) -> np.ndarray:
        """The start parameter of slot ``slot`` ``(L, ...)`` of each line."""
        return self._of_slot(self._slots.start, slot)

    @staticmethod
    def _of_slot(per_slot: np.ndarray, slot: np.ndarray) -> np.ndarray:
        """A per-slot field ``(L, S)`` read at slot ``slot`` ``(L, ...)``."""
        return np.take_along_axis(per_slot, slot.reshape(per_slot.shape[0], -1), axis=-1).reshape(slot.shape)

    # ── the block's scan ─────────────────────────────────────────────────

    @cached_property
    def _carried(self) -> np.ndarray:
        r"""What each transit's earlier pieces carry into each piece's start, ``(F, L, J, N)``.

        A scan: carrying by the transmission of each piece in turn never
        forms a difference of optical depths.
        """
        pieces = self._pieces
        p = pieces.pieces
        n_members = self._speed.shape[0]
        n_lines, n_pieces = p.live.shape
        carried = np.zeros((n_members, n_lines, 2, self.basis.size))
        before = np.zeros((n_members, n_lines, n_pieces, self.basis.size))
        members, lines = np.arange(n_members)[:, None, None], np.arange(n_lines)
        columns = self.basis.columns(p.panel)
        transmission = np.exp(-pieces.depth)
        for j in range(n_pieces):
            k = p.transit[:, j]
            entering = carried[:, lines, k]                                                 # (F, L, N)
            before[:, :, j] = entering
            step = transmission[:, :, j, None] * entering
            step[members, lines[:, None], columns[:, j]] += pieces.to_end[:, :, j]           # a piece's columns are distinct
            carried[:, lines, k] = np.where(p.live[:, j, None], step, entering)
        return before

    # ── the functionals ──────────────────────────────────────────────────

    @property
    def optical_depth(self) -> np.ndarray:
        r"""The optical depth :math:`\tau_k` of each traversal of the period, ``(*family, *batch, 2)``."""
        per_member = np.reshape(self.speed, self._family + (1,) * (len(self._lines) + 1))
        return per_member * self.period.optical_depth(self.sigma[self._panel_of_region])

    def _exit_integral(self) -> np.ndarray:
        """Each transit read forward: its integral attenuated to its exit, ``(F, L, 2, N)``."""
        pieces = self._pieces
        return pieces.per_transit(np.exp(-pieces.transit_depth(upstream=False)), pieces.to_end)

    def _entry_integral(self) -> np.ndarray:
        """Each transit read forward: its integral attenuated from its entry, ``(F, L, 2, N)``."""
        pieces = self._pieces
        return pieces.per_transit(np.exp(-pieces.transit_depth(upstream=True)), pieces.to_start)

    def _per_traversal(self, forward: Callable[[], np.ndarray], backward: Callable[[], np.ndarray]) -> np.ndarray:
        """A per-transit quantity read along each traversal of the period: ``forward()`` or, reversed, ``backward()``.

        Each side is formed only where the period reads a transit that way: a
        solid body's line is read forward alone.
        """
        period = self.period
        transit = self._flat(period.transit)[None, ..., None]
        reversed_, read_forward = self._flat(period.reversed), self._flat(period.present & ~period.reversed)
        shape = self._speed.shape + (2, self.basis.size)
        read = np.zeros(shape)
        if read_forward.any():
            read = np.where(read_forward[..., None], np.take_along_axis(forward(), transit, axis=-2), read)
        if reversed_.any():
            read = np.where(reversed_[..., None], np.take_along_axis(backward(), transit, axis=-2), read)
        return self._unflat_family(read)

    def outflow(self) -> np.ndarray:
        r"""The source integral :math:`B_k` of each basis function, attenuated to traversal :math:`k`'s exit, ``(*family, *batch, 2, N)``."""
        return self._per_traversal(self._exit_integral, self._entry_integral)

    def entry_response(self) -> np.ndarray:
        r"""Each basis function integrated over traversal :math:`k`, attenuated from its entry, ``(*family, *batch, 2, N)``.

        The outflow of the reversed traversal: reading a transit backward
        swaps its exit and its entry.
        """
        return self._per_traversal(self._entry_integral, self._exit_integral)

    def volterra(self, line_weight: np.ndarray) -> np.ndarray:
        r"""The vacuum Galerkin block :math:`\sum_L w_L \int u_i\,\psi^{\rm vac}_j`, ``(N, N)``, over each line's transits read forward.

        ``line_weight`` is :math:`w_L`, broadcast to the family over the line
        batch; the family's members are summed. A line's flux lives on its
        own transits in its own direction; the period's reversed traversals
        carry its cycle and belong to the opposite line. Both integrals are
        over arc length, :math:`v\,\mathrm dt`.
        """
        pieces = self._pieces
        p = pieces.pieces
        n_members = self._speed.shape[0]
        weight = np.broadcast_to(np.asarray(line_weight, dtype=float), self._family + self._lines).reshape(n_members, -1)
        weight = np.where(p.live, weight[..., None], 0.0)                                    # (F, L, J)
        # the flux the transit's earlier pieces carry into each piece, attenuated from its start, against the tests ...
        carried = np.einsum("flj,flji,fljn->ljin", weight, pieces.to_start, self._carried)    # (L, J, p+1, N)
        # ... and the piece's own triangle
        own = np.einsum("flj,fljab->ljab", weight, pieces.triangle)                          # (L, J, p+1, p+1)
        columns = self.basis.columns(p.panel)                                                 # (L, J, p+1)
        block = np.zeros((self.basis.size, self.basis.size))
        np.add.at(block, columns, carried)
        np.add.at(block, (columns[..., :, None], columns[..., None, :]), own)
        return block

    def angular_flux(self, t: np.ndarray, inflow: np.ndarray) -> np.ndarray:
        r"""The angular flux from each basis function at parameters ``t`` on the line, read forward, ``(*family, ..., q, N)``.

        ``t`` is ``(..., q)`` on the caller's lines; ``inflow`` is the period's
        inflow at each traversal's entry per basis function, ``(*family, ..., 2, N)``
        (:meth:`~.closure.LinePeriod.inflow`). Each point lies on a transit of
        its line: inside the body, never in a cavity or beyond a wall. The
        flux is the inflow attenuated from the transit's entry, plus each
        upstream piece's integral to its end, carried to the point (:meth:`at`).
        """
        return self.at(t).flux(inflow)

    def at(self, t: np.ndarray) -> "LineReading":
        r"""The transport to parameters ``t`` ``(..., q)`` on the lines, read forward (:class:`LineReading`).

        Each point is a crossing on a transit of its line: a slot's start, or a
        transit's exit, a piece end there, so the flux is
        each upstream piece's one-level integral (:meth:`_PieceSet.to_end`)
        carried to the point, with no scan over the basis. The depth from each
        piece to the point is summed from the point upstream, so no difference
        of depths is formed. A point inside a slot is refused: its lines are
        chorded through the panels refined at its level (:meth:`on_chord`).
        The place of a point on its line is the same for every member of the
        family.
        """
        exits, slots = self._pieces, self._slots
        p = exits.pieces
        t = self._flat(np.asarray(t, dtype=float))                                          # (L, q)
        begins = exits.begins                                                                 # (L, J)
        held = np.sum(begins[:, None, :] <= t[..., None], axis=-1) - 1                       # the piece holding each point
        piece = np.maximum(held, 0)
        slot = np.take_along_axis(p.slot, piece, axis=-1)
        # a slot ends at its closing crossing: the parameter, not start + length, which round apart by an ulp
        first = self._of_slot(slots.first, p.slot)
        at_start = (np.take_along_axis(p.lower == first, piece, axis=-1)) & (t == self._of_slot(slots.start, slot))
        at_exit = t == self._of_slot(slots.end, slot)
        if not np.all((held >= 0) & np.take_along_axis(p.live, piece, axis=-1) & (at_start | at_exit)):
            raise ValueError(
                "a point at which the flux is read is a crossing on a transit of its line: chord the lines through the "
                "panels refined at its level"
            )
        transit = np.take_along_axis(p.transit, piece, axis=-1)                              # (L, q)
        upstream = p.live[:, None, :] & (p.transit[:, None, :] == transit[..., None]) & (begins[:, None, :] < t[..., None])
        depth = np.where(upstream, exits.depth[:, :, None, :], 0.0)                           # (F, L, q, J)
        beyond = np.cumsum(depth[..., ::-1], axis=-1)[..., ::-1]                             # from each piece's start to the point
        after = np.concatenate([beyond[..., 1:], np.zeros_like(beyond[..., :1])], axis=-1)   # from each piece's end to the point
        carried = np.where(upstream, np.exp(-after), 0.0)                                    # (F, L, q, J)
        # each upstream piece's integral carried to the point, gathered onto its panel's columns
        vacuum = exits.on_columns(carried[..., None] * exits.to_end[:, :, None])           # (F, L, q, N)
        traversal = self.period.forward_traversal(self._unflat(transit))
        return LineReading(traversal, self._unflat_family(np.exp(-beyond[..., 0])), self._unflat_family(vacuum))


@dataclass(frozen=True, eq=False)
class LineReading:
    r"""The transport to points on lines, read forward: :math:`\psi = e^{-\tau_{\rm in}}\,\psi^{\rm in}_k + \psi^{\rm vac}`.

    The inflow term is a functional of the period's inflow, the vacuum term
    of the basis functions, so a caller can pair the inflow with any set of
    sources (the emission's basis functions, a current injected at a wall)
    while the vacuum term stays on the basis (:meth:`TraversalRule.at`).

    Attributes
    ----------
    traversal:
        The traversal of the period each point lies on, read forward, ``(*batch, q)``: the same for every member of
        the family.
    transmission:
        :math:`e^{-\tau_{\rm in}}`, the attenuation from that traversal's entry to the point, ``(*family, *batch, q)``.
    vacuum:
        :math:`\psi^{\rm vac}`, the flux each basis function produces at the point along its own transit,
        ``(*family, *batch, q, N)``.
    """

    traversal: np.ndarray
    transmission: np.ndarray
    vacuum: np.ndarray

    def inflow_term(self, inflow: np.ndarray) -> np.ndarray:
        r""":math:`e^{-\tau_{\rm in}}\,\psi^{\rm in}_k` per source, ``(*family, *batch, q, S)``, for an inflow per traversal and source ``(*family, *batch, 2, S)``."""
        inflow = np.asarray(inflow, dtype=float)
        traversal = self.traversal[..., None]
        traversal = traversal.reshape((1,) * (inflow.ndim - traversal.ndim) + traversal.shape)
        return self.transmission[..., None] * np.take_along_axis(inflow, traversal, axis=-2)

    def flux(self, inflow: np.ndarray) -> np.ndarray:
        r"""The angular flux from each basis function, ``(..., q, N)``, for its inflow ``(..., 2, N)``."""
        return self.inflow_term(inflow) + self.vacuum


__all__ = ["LineReading", "TraversalRule"]
