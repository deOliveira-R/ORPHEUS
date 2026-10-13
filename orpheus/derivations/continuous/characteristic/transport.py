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
as a family of speeds: its pieces, its nodes and its basis values are formed
once, and only the attenuations and the length weights depend on the
speed. A line parametrised by its own arc length, any line the kernel
chords in three dimensions, has speed 1.

**The pieces.** The lines are chorded through the basis's panel partition
(:attr:`~.basis.PanelBasis.partition`), or a refinement of it at the levels
a flux is read at, so each slot lies in one panel, and every slot is cut at
its middle and into pieces by two gradings:

* exponentially toward the ends it is read at, where the attenuation
  concentrates every integrand: piece ends at :math:`h\rho^k` mean free
  paths of the family's steepest attenuation from each such end, until 64
  mean free paths of its gentlest (:func:`~.grading.exponential_ends`), with
  :math:`h` and :math:`\rho` derived from the node count
  (:class:`~.grading.ExponentialGrading`). A slot of optical width at most
  :math:`h` at the steepest is not graded. Three sets of
  pieces are graded so (:class:`_Toward`): toward each slot's end, for an
  outflow and the flux at a crossing; toward its start, for the entry
  response; toward both ends, for the block's Volterra integral.
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

**One integrand, two levels.** Every attenuated integral along a line has
the one integrand :math:`v\,w\,e^{-\Sigma v d}` (:func:`_attenuated_length`).
A quantity read at a slot end the pieces are graded toward (the outflow, the
entry response, the flux at a crossing) is integrated on each piece's own
nodes, one level (:meth:`_PieceSet.to_end`, :meth:`_PieceSet.to_start`). A
quantity read where no piece end is graded toward, the flux at a node inside
a piece, is integrated on an inner rule graded toward that node, two levels
(:meth:`TraversalRule._attenuated`): the integral out of a wide middle piece
lives in its last few mean free paths, which an ungraded rule on the piece's
own nodes misses at O(1) (`[M]` qa's ``probe_thick_psi.py``, 2026-10-06:
0.54 in the flux just past the middle piece of a 1000-mfp slot).

**The Volterra block** pairs the flux at a piece's own nodes with the test
functions there: the flux at a node is what the transit's earlier pieces
carry in, attenuated from the piece's start, plus the integral from the
piece's start to the node. The outer rule integrates the smooth product
:math:`u_i \psi_j`; the inner integral is the graded one, and what the
earlier pieces carry is each piece's graded integral to its end, scanned
along the transit (:attr:`TraversalRule._carried`).

The design and its rulings: ``.claude/plans/characteristic_reference_architecture.md``
(P1 step (b), second rung).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum
from functools import cached_property

import numpy as np

from orpheus.derivations.common.quadrature import gauss_legendre
from orpheus.geometry.chord import Chord, RadialImage
from orpheus.geometry.line import Line

from .basis import PanelBasis
from .grading import ExponentialGrading, exponential_ends, halvings
from .closure import LinePeriod
from .walls import Walls


def _attenuated_length(speed: np.ndarray, weight: np.ndarray, sigma: np.ndarray, distance: np.ndarray) -> np.ndarray:
    r""":math:`v\,w\,e^{-\Sigma v d}`: a node's arc-length element :math:`v\,w` (its weight :math:`w` in the line parameter at speed
    :math:`v`), attenuated over the parameter distance :math:`d` at cross section :math:`\Sigma`; the arguments broadcast.

    The one integrand of every attenuated integral along a line: the one-level
    integrals of a piece set (:meth:`_PieceSet._integral`), the graded inner rule
    (:meth:`TraversalRule._attenuated`) and the flux carried into a node of the
    block (:meth:`TraversalRule.volterra`).
    """
    return speed * weight * np.exp(-(speed * sigma) * distance)


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


class _Toward(Enum):
    r"""The ends of every slot a set of pieces is graded toward: ``(start, end)``.

    ``EXITS``: each slot's end, where a transit read forward exits and where
    a point on a level is read (:meth:`TraversalRule.at`). ``ENTRIES``: each
    slot's start, where a transit read backward exits. ``BOTH``: the block's
    Volterra integral, which reads the flux at every node and so meets the
    layers at both ends. Every set is graded for the rule's
    :attr:`~TraversalRule.points` (:class:`~.grading.ExponentialGrading`).
    The block's inner rule is graded for its own node count,
    :attr:`~TraversalRule.inner_points`: at the fixed ratio 2 it replaced, a
    piece's integral to its end at :math:`\Sigma = 30`, speeds 1 and 100,
    against 24 inner points, missed by 1.0e-6 at 4 inner points, 5.8e-13 at 8
    and 3.3e-15 at 12 (`[M]` 2026-10-12, the elegance review,
    ``scratch/characteristic_architecture/p3/step3/elegance/probe_shape.py``).
    """

    EXITS = (False, True)
    ENTRIES = (True, False)
    BOTH = (True, True)


@dataclass(frozen=True, eq=False)
class _PieceSet:
    r"""A rule's slots cut into pieces graded toward some of their ends (:class:`_Toward`), and the quadrature on them.

    Each slot is graded exponentially toward the ends it is read at
    (:func:`~.grading.exponential_ends`, over the family's range of
    attenuations) and, halving, toward the end nearer the line's closest
    approach (:meth:`TraversalRule._branch_edges`). The integrals over it are
    one level, on each piece's :attr:`TraversalRule.points` nodes, and exact
    to the rule's order where they are read at a graded end.
    """

    rule: "TraversalRule"
    pieces: _Pieces

    @classmethod
    def of(cls, rule: "TraversalRule", toward: _Toward) -> "_PieceSet":
        r"""The pieces of ``rule``'s slots graded toward the ends ``toward``."""
        s = rule._slots
        steepest, gentlest = rule._attenuations(rule.sigma[s.panel])
        grading = ExponentialGrading.of(rule.points)
        to_start, to_end = toward.value
        # graded toward both ends, each half of a slot is graded toward its own; toward one, the whole slot is. A slot
        # is always cut at its middle: `[M]` 2026-10-11 (``scratch/characteristic_architecture/p3/step3/impl/
        # exits_accuracy.py``, a three-region cylinder's outflow at 8 points against 24), the whole slot graded toward
        # its end alone left a thin slot one piece, 7.9e-11 off on a 3-D line; the halves graded toward the end
        # alone left the far half one piece, 1.2e-13 off on the in-plane line at speeds 1 to 100
        half = s.length / 2.0
        middle = half if to_start and to_end else (np.zeros_like(s.length) if to_end else s.length)
        graded = [rule._branch_edges(), half[..., None]]
        if to_start:
            graded.append(exponential_ends(np.zeros_like(s.length), middle, steepest, gentlest, grading))
        if to_end:
            graded.append(exponential_ends(s.length, middle, steepest, gentlest, grading))
        ends = np.sort(np.concatenate([np.zeros_like(s.length)[..., None], s.length[..., None], *graded], axis=-1), axis=-1)
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
    def nodes(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        r"""The :attr:`TraversalRule.points` Gauss-Legendre nodes of every piece: distance from the slot's start, weight, and the basis values there, ``(L, J, n)`` and ``(L, J, n, p + 1)``."""
        p, rule = self.pieces, self.rule
        distance, orbit, weight = rule._nodes(p.slot, p.lower, p.upper, rule.points)
        return distance, weight, rule.basis.values(orbit, np.broadcast_to(p.panel[..., None], orbit.shape))

    def _integral(self, distance_to: np.ndarray) -> np.ndarray:
        r"""Each piece's integral attenuated over ``distance_to`` from each node ``(L, J, n)``, on its own nodes, ``(F, L, J, p + 1)``."""
        _, weight, values = self.nodes
        attenuated = _attenuated_length(self.rule._speed[..., None, None], weight, self.sigma[..., None], distance_to)
        return (attenuated[..., None, :] @ values)[..., 0, :]

    @cached_property
    def to_end(self) -> np.ndarray:
        r"""Each piece's integral attenuated to its end, ``(F, L, J, p + 1)``: read at the ends the pieces are graded toward.

        `[M]` 2026-10-11 (P3 step 3, ``scratch/characteristic_architecture/p3/step3/impl/one_level_probe.py``):
        the outflow at 8 points on pieces graded toward both ends, against the
        two-level rule at 24, equals the two-level outflow on the
        three-region, thin and thick cylinders at one speed (5.7e-15,
        1.5e-14, 2.4e-15), and reads 7.6e-13 against 9.6e-14 on a thick
        cylinder at speeds 1 to 1000 (1.7e-14 both at 12 points). The flux
        at a piece end inside a slot is not such a read: the block's scan
        carries the graded integral (:attr:`TraversalRule._to_end`).
        """
        distance, _, _ = self.nodes
        return self._integral(self.pieces.upper[..., None] - distance)

    @cached_property
    def to_start(self) -> np.ndarray:
        r"""Each piece's integral attenuated to its start, ``(F, L, J, p + 1)``: read at the ends the pieces are graded toward."""
        distance, _, _ = self.nodes
        return self._integral(distance - self.pieces.lower[..., None])

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
    points:
        Gauss-Legendre points per piece.
    inner_points:
        Gauss-Legendre points per interval of the inner rule (the attenuated integral).
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
    points: int
    inner_points: int
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
        if self.points < 1 or self.inner_points < 1:
            raise ValueError(f"a rule has at least one point; got {self.points}, {self.inner_points}")
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
        points: int,
        inner_points: int,
        level: tuple[np.ndarray, np.ndarray] | None = None,
        speed: np.ndarray | float = 1.0,
    ) -> "TraversalRule":
        r"""The rule along ``lines`` through ``basis`` and ``walls``, with ``sigma_t`` the total cross section of each REGION, ``(n,)``.

        ``level`` is each line's exact level, a radius and its half-chord
        there (:meth:`~orpheus.geometry.chord.ConcentricPartition.chord`);
        ``speed`` is :attr:`speed`.
        """
        return cls.on_chord(basis.partition.chord(lines, level), basis, walls, sigma_t, points, inner_points, speed)

    @classmethod
    def on_chord(
        cls,
        chord: Chord,
        basis: PanelBasis,
        walls: Walls,
        sigma_t: np.ndarray,
        points: int,
        inner_points: int,
        speed: np.ndarray | float = 1.0,
    ) -> "TraversalRule":
        r"""The rule along the lines of ``chord``, a chord of ``basis``'s panel partition (:meth:`of`, for a chord already formed, such as :meth:`~orpheus.geometry.chord.ConcentricPartition.rays_from`'s)."""
        period = LinePeriod.of(chord, walls.on(chord.partition))
        return cls(period, basis, basis.on_panels(sigma_t), points, inner_points, speed)

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
        return per_line.reshape(*self._lines, *per_line.shape[1:])

    def _unflat_family(self, per_member: np.ndarray) -> np.ndarray:
        """``(F, L, ...)`` to ``(*family, *batch, ...)``."""
        return per_member.reshape(*self._family, *self._lines, *per_member.shape[2:])

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
        :math:`r_k/|P\Omega|`. The pieces halve toward that end until each is
        no wider than its distance to the branch point (:func:`~.grading.halvings`). `[M]` 2026-10-06: one piece at 16 points missed by 2.5e-12
        at :math:`b = 10^{-4}` (the test-architect's ``probe_turn2.py``) and by
        8.5e-7 on a hollow sphere of cavity radius 0.01 (qa's ``probe_branch.py``).

        A change of variable :math:`c = b + (c_{\rm far} - b)u^2` on the slots
        ending at the closest approach was tried first and retired: it leaves
        a branch point in the Jacobian at :math:`u = \pm i\sqrt{2b/(c_{\rm far} - b)}`,
        and it covered only those slots.
        """
        s = self._slots
        reach = (s.near_end != 0) & (s.branch_distance > 0.0)
        layers = halvings(s.length, np.where(reach, s.branch_distance, s.length))   # 0 where ungraded
        level = np.arange(1, np.finfo(float).nmant + 1)
        depth = s.length[..., None] * 0.5**level
        at = np.where(s.near_end[..., None] > 0, s.length[..., None] - depth, depth)
        return np.where(level <= layers[..., None], at, 0.0)

    # ── the pieces ───────────────────────────────────────────────────────

    @cached_property
    def _toward_exits(self) -> _PieceSet:
        r"""The pieces graded toward each slot's end: a forward transit's outflow and the flux at a point (:meth:`at`)."""
        return _PieceSet.of(self, _Toward.EXITS)

    @cached_property
    def _toward_entries(self) -> _PieceSet:
        r"""The pieces graded toward each slot's start: a transit's integral attenuated from its entry."""
        return _PieceSet.of(self, _Toward.ENTRIES)

    @cached_property
    def _toward_both(self) -> _PieceSet:
        r"""The pieces graded toward both ends of each slot: the block's Volterra integral (:meth:`volterra`)."""
        return _PieceSet.of(self, _Toward.BOTH)

    @property
    def piece_counts(self) -> np.ndarray:
        r"""The number of live pieces on each line graded toward the exits, ``(*batch,)``: a line's share of :attr:`extent`, before the padding.

        Counted on the set every consumer reads (the outflow, and the flux at a
        point), so sizing a batch builds no set its consumer does not; the
        block's set graded toward both ends is of the same order.
        """
        return self._unflat(self._toward_exits.pieces.live.sum(axis=-1))

    def _attenuations(self, sigma: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        r"""The steepest and the gentlest attenuation per unit parameter over the family, for cross sections ``sigma`` ``(L, ...)``.

        :math:`\Sigma v` at each line's fastest and slowest speed: the range
        an exponential grading must resolve (:func:`~.grading.exponential_ends`).
        """
        speed = self._speed
        fastest, slowest = (v.reshape(v.shape + (1,) * (sigma.ndim - 1)) for v in (speed.max(axis=0), speed.min(axis=0)))
        return sigma * fastest, sigma * slowest

    @property
    def extent(self) -> int:
        r"""The number of piece slots the rule's arrays hold at most: family members times lines times the most pieces on one line (:attr:`piece_counts`).

        Every node array scales with it, so a caller bounds a batch's memory by
        it; a rule of one line holds its whole family, whatever its extent.
        """
        return int(self._toward_exits.pieces.live.size) * self._speed.shape[0]

    # ── placing nodes and the one attenuated integral ────────────────────

    def _orbit_coordinate(self, t: np.ndarray) -> np.ndarray:
        """The orbit coordinate at parameters ``t`` ``(L, ...)`` of each line."""
        flat = t.reshape(*self._lines, -1)
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

        :math:`\int u_m(c(t))\,e^{-\Sigma v |{\rm stop} - t|}\,v\,\mathrm{d}t` over the
        interval between the two, per member of the family of speeds
        :math:`v`, on one rule graded exponentially toward ``stop`` over the
        family's range of attenuations. ``slot`` broadcasts against ``start``
        and ``stop``, ``(L, ...)``. Returns ``(F, L, ..., p + 1)``.

        Only the live intervals of the graded rule are evaluated, each line's
        packed first (:class:`_Packing`): an integral over a thin stretch has
        one interval and one over a thick stretch up to
        :math:`\lceil\log_\rho(64\,v_{\max}/(h\,v_{\min}))\rceil + 2`
        (:class:`~.grading.ExponentialGrading`'s :math:`h` and :math:`\rho`
        at :attr:`inner_points`), and padding every integral to the batch's thickest evaluated 4 to 6 times the
        live nodes (`[M]` 2026-10-07, #586: three-region and tau = 30 white
        cylinders, 8 points, a sixteenth of the lines). The nodes and the
        basis values there are formed once for the whole family.
        """
        shape = np.broadcast_shapes(np.shape(slot), start.shape, stop.shape)
        n_lines = shape[0]
        slot, start, stop = (np.broadcast_to(a, shape).reshape(n_lines, -1) for a in (slot, start, stop))
        panel = self._of_slot(self._slots.panel, slot)                                    # (L, O)
        sigma = self.sigma[panel]
        ends = exponential_ends(stop, start, *self._attenuations(sigma), ExponentialGrading.of(self.inner_points))   # (L, O, K + 2)
        intervals = _Packing.of(ends[..., 1:] > ends[..., :-1])
        owner = intervals.owner
        lower, upper = intervals.gather(ends[..., :-1], 0.0), intervals.gather(ends[..., 1:], 0.0)   # (L, M)
        distance, orbit, weight = self._nodes(self._of_slot(slot, owner), lower, upper, self.inner_points)
        values = self.basis.values(orbit, np.broadcast_to(self._of_slot(panel, owner)[..., None], orbit.shape))
        to_stop = np.abs(self._of_slot(stop, owner)[..., None] - distance)
        attenuated = _attenuated_length(self._speed[..., None, None], weight, self._of_slot(sigma, owner)[..., None], to_stop)
        # per interval, every member's attenuated length weights against the basis values there: (L, M, F, q) @ (L, M, q, p + 1)
        per_interval = np.moveaxis(attenuated, 0, -2) @ values                               # (L, M, F, p + 1)
        per_owner = intervals.scatter(per_interval)                                          # (L, O, F, p + 1)
        return np.moveaxis(per_owner, -2, 0).reshape(-1, *shape, self.basis.per_panel)

    @cached_property
    def _to_end(self) -> np.ndarray:
        r"""Each piece graded toward both ends, its integral attenuated to its end on the graded inner rule, ``(F, L, J, p + 1)``: exact at every piece end, for the block's scan (:attr:`_carried`)."""
        p = self._toward_both.pieces
        return self._attenuated(p.slot, p.lower, p.upper)

    # ── the block's scan ─────────────────────────────────────────────────

    @cached_property
    def _carried(self) -> np.ndarray:
        r"""What each transit's earlier pieces carry into each piece's start, ``(F, L, J, N)``.

        A scan: carrying by the transmission of each piece in turn never
        forms a difference of optical depths.
        """
        p = self._toward_both.pieces
        n_members = self._speed.shape[0]
        n_lines, n_pieces = p.live.shape
        carried = np.zeros((n_members, n_lines, 2, self.basis.size))
        before = np.zeros((n_members, n_lines, n_pieces, self.basis.size))
        members, lines = np.arange(n_members)[:, None, None], np.arange(n_lines)
        columns = self.basis.columns(p.panel)
        transmission = np.exp(-self._toward_both.depth)
        for j in range(n_pieces):
            k = p.transit[:, j]
            entering = carried[:, lines, k]                                                 # (F, L, N)
            before[:, :, j] = entering
            step = transmission[:, :, j, None] * entering
            step[members, lines[:, None], columns[:, j]] += self._to_end[:, :, j]            # a piece's columns are distinct
            carried[:, lines, k] = np.where(p.live[:, j, None], step, entering)
        return before

    # ── the functionals ──────────────────────────────────────────────────

    @property
    def optical_depth(self) -> np.ndarray:
        r"""The optical depth :math:`\tau_k` of each traversal of the period, ``(*family, *batch, 2)``."""
        per_member = np.reshape(self.speed, self._family + (1,) * (len(self._lines) + 1))
        return per_member * self.period.optical_depth(self.sigma[self._panel_of_region])

    def _exit_integral(self) -> np.ndarray:
        """Each transit read forward: its integral attenuated to its exit, on the pieces graded toward the exits, ``(F, L, 2, N)``."""
        exits = self._toward_exits
        return exits.per_transit(np.exp(-exits.transit_depth(upstream=False)), exits.to_end)

    def _entry_integral(self) -> np.ndarray:
        """Each transit read forward: its integral attenuated from its entry, on the pieces graded toward the entries, ``(F, L, 2, N)``."""
        entries = self._toward_entries
        return entries.per_transit(np.exp(-entries.transit_depth(upstream=True)), entries.to_start)

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
        both = self._toward_both
        p = both.pieces
        n_members = self._speed.shape[0]
        weight = np.broadcast_to(np.asarray(line_weight, dtype=float), self._family + self._lines).reshape(n_members, -1)
        weight = np.where(p.live, weight[..., None], 0.0)                                    # (F, L, J)
        distance, node_weight, tests = both.nodes                                             # (L, J, n), (L, J, n, p+1)
        speed = self._speed[..., None, None]                                                 # (F, L, 1, 1)
        weighted_length = weight[..., None] * speed * node_weight                            # (F, L, J, n): w v dt
        # the flux at each node, summed over the family with the weights: carried in from the piece's start ...
        entering = weight[..., None] * _attenuated_length(speed, node_weight, both.sigma[..., None], distance - p.lower[..., None])
        carried = np.moveaxis(entering, 0, -1) @ np.moveaxis(self._carried, 0, -2)          # (L, J, n, N)
        # ... plus the piece's own integral up to the node
        own = self._attenuated(p.slot[..., None], np.broadcast_to(p.lower[..., None], distance.shape), distance)
        triangle = (weighted_length[..., None] * own).sum(axis=0)                           # (L, J, n, p+1)
        # the test functions paired with that flux at the nodes
        paired = np.swapaxes(tests, -1, -2)                                                  # (L, J, p+1, n)
        columns = self.basis.columns(p.panel)                                                 # (L, J, p+1)
        block = np.zeros((self.basis.size, self.basis.size))
        np.add.at(block, columns, paired @ carried)
        np.add.at(block, (columns[..., :, None], columns[..., None, :]), paired @ triangle)
        return block

    def angular_flux(self, t: np.ndarray, inflow: np.ndarray) -> np.ndarray:
        r"""The angular flux from each basis function at parameters ``t`` on the line, read forward, ``(*family, ..., q, N)``.

        ``t`` is ``(..., q)`` on the caller's lines; ``inflow`` is the period's
        inflow at each traversal's entry per basis function, ``(*family, ..., 2, N)``
        (:meth:`~.closure.LinePeriod.inflow`). Each point lies on a transit of
        its line: inside the body, never in a cavity or beyond a wall. The
        flux is the inflow attenuated from the transit's entry, plus what the
        transit's earlier pieces carry in, plus the piece's own integral up to
        the point (:meth:`at`).
        """
        return self.at(t).flux(inflow)

    def at(self, t: np.ndarray) -> "LineReading":
        r"""The transport to parameters ``t`` ``(..., q)`` on the lines, read forward (:class:`LineReading`).

        Each point is a crossing on a transit of its line: a slot's start, or a
        transit's exit. The pieces are graded toward it there, so the flux is
        each upstream piece's one-level integral (:meth:`_PieceSet.to_end`, on the pieces graded toward the exits)
        carried to the point, with no scan over the basis. The depth from each
        piece to the point is summed from the point upstream, so no difference
        of depths is formed. A point inside a slot is refused: its lines are
        chorded through the panels refined at its level (:meth:`on_chord`).
        The place of a point on its line is the same for every member of the
        family.
        """
        exits, slots = self._toward_exits, self._slots
        p = exits.pieces
        t = self._flat(np.asarray(t, dtype=float))                                          # (L, q)
        begins = np.where(p.live, self._slot_start(p.slot) + p.lower, np.inf)               # (L, J)
        held = np.sum(begins[:, None, :] <= t[..., None], axis=-1) - 1                       # the piece holding each point
        piece = np.maximum(held, 0)
        slot = np.take_along_axis(p.slot, piece, axis=-1)
        # a slot ends at its closing crossing: the parameter, not start + length, which round apart by an ulp
        at_start = (np.take_along_axis(p.lower, piece, axis=-1) == 0.0) & (t == self._of_slot(slots.start, slot))
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
