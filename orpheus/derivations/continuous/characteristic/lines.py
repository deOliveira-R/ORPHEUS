r"""The rules over lines a body's quadratures are built from, and the sources a line carries.

Two measures integrate the same transport over lines: the lines of space,
for the Galerkin block (:class:`~.assembly.LineRule`), and the directions
at one point, for the reading (:class:`~.reading.PointRule`). Both are
built from the one-dimensional rules of this module, so the gradings that
resolve the transport are written once:

* the **impact parameter** :math:`b` on the radial charts, panel by panel
  under the visibility-cone substitution :math:`y = \sqrt{r_{k+1}^2 - b^2}`
  (:func:`impact_rule`), graded toward every singularity of the line's
  transport in :math:`y` (:class:`ImpactPanels`);
* the cylinder's **polar angle** and the slab's **cosine**, graded toward
  grazing and toward the normal (:func:`polar_rule`, :func:`cosine_rule`).

Each impact node carries its panel's top and its exact half-chord there
(:class:`ImpactRule`), which its line is chorded from
(:meth:`~orpheus.geometry.chord.ConcentricPartition.chord`'s ``level``):
the line itself stores :math:`b` to an ulp, which near a tangency loses
the half-chord (#590).

The design and its rulings: ``.claude/plans/characteristic_reference_architecture.md``
(P1 step (b), third rung, and rung 5b).
"""

from __future__ import annotations

from collections.abc import Callable
from collections.abc import Iterator
from dataclasses import dataclass, replace

import numpy as np

from orpheus.derivations.common.quadrature import Quadrature1D, composite_gauss_legendre, gauss_legendre
from orpheus.geometry.chart import Chart
from orpheus.geometry.chord import ConcentricPartition

from .basis import PanelBasis

from .closure import DiffuseWalls
from .grading import VANISHING_DEPTH, exponential_ends, graded_ends, halvings
from .transport import TraversalRule
from .walls import Walls


def _grazing_ends(speed_to_coordinate: Callable[[np.ndarray], np.ndarray], thinnest: float, thickest: float) -> list[float]:
    r"""The ends of a direction coordinate's interval, graded toward its two features.

    The coordinate runs from the grazing direction (projected speed
    :math:`v = |P\Omega| = 0`) to the normal one (:math:`v = 1`), and
    ``speed_to_coordinate`` maps :math:`v` to it (the slab's cosine is
    :math:`v` itself, the cylinder's polar angle :math:`\arcsin v`). A line
    crossing an optical width :math:`\tau` normally crosses
    :math:`\tau / v` along it, so:

    * toward :math:`v = 0`, every panel's attenuation :math:`e^{-\tau_P s}`
      changes until :math:`\tau_P s` reaches the depth where it vanishes
      (:data:`~.grading.VANISHING_DEPTH`, 64): halved down to
      :math:`v = \tau_{\min}/64`, with ``thinnest`` the thinnest absorbing
      panel's width. Stopping at :math:`\tau_{\min}` left the totals exact
      (measured 2026-10-06,
      ``scratch/characteristic_architecture/p1_step_b3/margin/study.py``, to
      5e-16) but a slab's block entry off by 2.0e-6 at 8 points (qa's
      ``q3``: basis layers 6, :math:`\Sigma = 0.01`; 2 more halvings gave
      2.6e-10, 6 gave 9.1e-13);
    * toward :math:`v = 1`, :math:`e^{-\tau s}` with :math:`s = 1/v` is an
      exponential layer of width :math:`1/\tau` in :math:`s`: the
      exponential ends of the traversal rule at :math:`2^k` over the
      thickest normal optical depth ``thickest``, on :math:`s \in [1, \infty)`.
      Without them a thick slab's transmission missed :math:`2E_3(\tau)` by
      4.1e-11 at :math:`\tau = 8` and 1.8e-4 at :math:`\tau = 30`, 8 points
      (measured 2026-10-06 by the test-architect, ``gates/slab_tw.py``); with
      them on :math:`s \in [1, 2]` only, 3.8e-11 remained at :math:`\tau = 8`,
      from the layer's tail below :math:`\mu = 1/2`.

    With no absorbing panel the coordinate is not graded. At most
    ``np.finfo(float).nmant`` halvings.
    """
    if thinnest <= 0.0 or not np.isfinite(thinnest):
        return [0.0, float(speed_to_coordinate(np.array(1.0)))]
    layers = int(halvings(1.0, thinnest / VANISHING_DEPTH))
    toward_grazing = graded_ends(0.0, 1.0, True, False, layers, 0.5)
    toward_normal = 1.0 / exponential_ends(np.array(1.0), np.array(np.inf), np.array(thickest))
    speeds = np.unique(np.concatenate([[0.0, 1.0], toward_grazing, toward_normal]))
    return list(np.unique(speed_to_coordinate(speeds)))


@dataclass(frozen=True, eq=False)
class ImpactRule:
    r"""An impact-parameter rule: its nodes :math:`b`, its weights for :math:`\mathrm db`, and each node's exact level.

    ``top`` is the upper end :math:`r_{k+1}` of the node's panel and
    ``half_chord`` is :math:`y = \sqrt{r_{k+1}^2 - b^2}` there, the
    substitution's own variable where the panel is substituted: the pair is
    the line's exact level (:meth:`~orpheus.geometry.chart.RadialImage.at_level`),
    and a consumer reads :math:`\sqrt{c^2 - b^2} = \sqrt{(c - r_{k+1})(c + r_{k+1}) + y^2}`
    for :math:`c \ge r_{k+1}` without the cancellation of :math:`c - b`.
    """

    b: np.ndarray
    weights: np.ndarray
    top: np.ndarray
    half_chord: np.ndarray


def impact_rule(ends: np.ndarray, sigma: np.ndarray, points: int, toward_tangency: np.ndarray, panels: int | None = None) -> ImpactRule:
    r"""The impact-parameter rule on :math:`[0, r_n]`, panel by panel of the panel partition ``ends``.

    On each panel :math:`[r_k, r_{k+1}]` the variable is the chord
    half-length :math:`y = \sqrt{r_{k+1}^2 - b^2}` (the visibility-cone
    substitution, which absorbs the square-root end at :math:`r_{k+1}`),
    graded toward each of its singularities until every piece is no wider
    than its distance to them (hp), and at :math:`2^j` mean free paths of
    the thickest panel from :math:`[r_k, r_{k+1}]` outward:

    * toward :math:`y = 0`, the nearest feature of the transport at the
      tangency (:class:`ImpactPanels`): the next radius out, just past
      the panel in :math:`b` but :math:`\sqrt{2 r\,\delta}` away in :math:`y`
      (measured 2026-10-06: on a sphere with a 1e-3 first region, Gauss in
      :math:`b` on the wide middle panel [0.52, 0.88] left 1e-6 at 8 points,
      the neighbours at 0.952 and 1.0), the turning slot's layer, and the
      closure's pole;
    * toward :math:`b = r_k`, the point :math:`b = 0` at :math:`y = r_{k+1}`: the
      Abel transform of the turning panel's odd modes and the
      substitution's Jacobian :math:`y/b` are singular there.

    The panel touching :math:`b = 0` is smooth in :math:`b^2` (its turning
    panel is even, or a cavity), so its lower half is plain Gauss in
    :math:`b` and its upper half is graded as above. ``sigma`` is each
    panel's total cross section, ``(P,)``. With ``panels``, only the first
    that many panels carry nodes; the gradings still read the panels above
    them (a point's rule, :mod:`~orpheus.derivations.continuous.characteristic.reading`).

    ``toward_tangency`` is that distance per panel.
    """
    pts, wts, tops, chords = [], [], [], []
    for k, (lo, hi) in enumerate(zip(ends[:-1], ends[1:])):
        if panels is not None and k >= panels:
            break
        if lo == 0.0:
            half = gauss_legendre(0.0, hi / 2.0, points)
            pts.append(half.pts)
            wts.append(half.wts)
            tops.append(np.full(points, hi))
            chords.append(np.sqrt((hi - half.pts) * (hi + half.pts)))
            lo = hi / 2.0
        span = float(np.sqrt((hi - lo) * (hi + lo)))
        to_centre = lo * lo / (hi + span)                         # r_{k+1} - span, without cancellation
        y_ends = [
            *graded_ends(0.0, span, True, False, halvings(span, float(toward_tangency[k])), 0.5),
            *graded_ends(0.0, span, False, True, halvings(span, to_centre), 0.5),
            *exponential_ends(np.array(0.0), np.array(span), np.array(float(np.max(sigma[k:])))),
        ]
        y = composite_gauss_legendre(np.unique(np.clip(y_ends, 0.0, span)), points)
        # b^2 = r_k^2 + (span^2 - y^2): near the lower end b -> r_k without the cancellation of hi^2 - y^2
        b = np.sqrt(lo * lo + (span - y.pts) * (span + y.pts))
        pts.append(b)
        wts.append(y.wts * y.pts / b)
        tops.append(np.full(b.shape, hi))
        chords.append(y.pts)
    return ImpactRule(*(np.concatenate(x) for x in (pts, wts, tops, chords)))


@dataclass(frozen=True, eq=False)
class ImpactPanels:
    r"""The impact parameter's panels, and per panel the features of a line's transport near its tangency :math:`b = r_{k+1}`.

    A line of impact parameter just below :math:`r_{k+1}`, at
    :math:`y = \sqrt{r_{k+1}^2 - b^2} \to 0`, turns in panel :math:`k` across
    a slot of in-plane length :math:`2y`, then crosses the shells above, of
    in-plane optical depth :math:`\tau_{\rm out}`, and reaches the outer
    wall. Along it, at projected speed :math:`s`, three features sit near
    :math:`y = 0`:

    * the **next radius out**, :math:`b = r_{k+2}`, at the imaginary
      :math:`y = \pm i\sqrt{r_{k+2}^2 - r_{k+1}^2}` (its chord's square root);
    * the **layer**: the turning slot transmits :math:`e^{-2\Sigma_k y/s}`,
      which changes over :math:`s/(2\Sigma_k)`;
    * the **pole**: the period's cycle product
      :math:`\Pi = a\,e^{-(\tau_{\rm out} + 2\Sigma_k y)/s}` reaches 1 at
      :math:`y = -(s\,(-\ln a) + \tau_{\rm out})/(2\Sigma_k)`, a pole of the
      closure :math:`1/(1 - \Pi)` (:class:`~.closure.LinePeriod`) that nears
      the tangency where the shells above are thin or void and :math:`a \to 1`.

    Their distances from the tangency are, in order, the next radius's
    half-chord, :math:`s/(2\Sigma_k)` and
    :math:`(s\,(-\ln a) + \tau_{\rm out})/(2\Sigma_k)`: the first is the
    body's, the last two depend on the line's own speed :math:`s`. So the
    speed-free data are held here per panel, read once, and graded per
    speed by :meth:`distances`; :meth:`rule` builds the impact rule at a
    speed. On the cylinder each polar angle's impact rule is graded at its
    own :math:`\sin\theta` (#587, :func:`impact_per_polar`), and on the
    sphere at 1.
    The pole is absent where the outer wall returns nothing specularly and
    where it returns everything (amplitude 1: the closure is regular, the
    source integral vanishing with :math:`1 - \Pi`); a void panel's
    transport does not change with :math:`y`.

    `[M]` 2026-10-08, rung 5b (qa and the test-architect), each at 8 line
    points before the layer and the pole were graded: a sphere's reading on
    its wall under :math:`a = 0.99` off by 9.3e-5 and the block's
    :math:`1^{\mathsf T}K1` by 5.3e-10; a void outer region under
    :math:`a = 0.99` off by 5.1e-4 at the interface and 2.6e-6 in the block;
    a three-region cylinder off by 7.9e-8 on its vacuum wall and 1.0e-6 at
    an interface. Grading every polar angle at the slowest one instead of
    its own cost the ``ABA`` cylinder 308 480 lines per block against
    102 320, for the same k to 8e-15 (2026-10-08,
    ``scratch/characteristic_architecture/p1_step_c/``).

    Attributes
    ----------
    ends:
        The impact panels' ends (:func:`impact_panels`), ``(P + 1,)``.
    sigma:
        Each panel's total cross section, ``(P,)``.
    to_next_radius:
        The half-chord of the next radius out at each panel top, ``(P,)``; infinite at the outermost.
    outside:
        :math:`\tau_{\rm out}`, the in-plane optical depth above each panel top along its tangent line, ``(P,)``.
    pole:
        :math:`-\ln a` of the outer wall, infinite where the closure has no pole.
    """

    ends: np.ndarray
    sigma: np.ndarray
    to_next_radius: np.ndarray
    outside: np.ndarray
    pole: float

    @classmethod
    def of(cls, chart: Chart, ends: np.ndarray, sigma: np.ndarray, amplitude: float) -> "ImpactPanels":
        r"""The impact panels of the panel partition ``ends`` of cross sections ``sigma`` (:func:`impact_panels`), under an outer
        wall of specular amplitude ``amplitude``.

        :math:`\tau_{\rm out}` and the next radius's half-chord are read from
        the kernel's chord of the tangent lines :math:`b = r_{k+1}` through the
        panels, at unit speed, so no chord length is spelled twice.
        """
        ends, sigma = impact_panels(np.asarray(ends, dtype=float), np.asarray(sigma, dtype=float))
        tops = ends[1:]
        unit_speed = np.zeros((tops.size, len(chart.line_domain().axes)))
        unit_speed[:, 0] = tops
        if unit_speed.shape[1] == 2:
            unit_speed[:, 1] = np.pi / 2.0                                       # the cylinder's line normal to the axis
        chord = ConcentricPartition(chart, tuple(ends)).chord(chart.line_domain().lines(unit_speed))
        region = np.minimum(chord.slot_region, sigma.size - 1)
        outside = np.sum(np.where(chord.slot_region < sigma.size, sigma[region] * chord.traversed_length, 0.0), axis=-1)
        next_radius = np.append(chord.image.half_chord_at(ends)[np.arange(tops.size - 1), np.arange(2, ends.size)], np.inf)
        pole = -np.log(amplitude) if 0.0 < amplitude < 1.0 else np.inf
        return cls(ends, sigma, next_radius, outside, float(pole))

    def distances(self, speed: float) -> np.ndarray:
        r"""Per panel, the distance in :math:`y` from the tangency to the nearest feature of a line of projected speed ``speed``, ``(P,)``.

        The layer's reach is the speed itself; the pole's, where the closure
        has one, is :math:`s(-\ln a) + \tau_{\rm out}`.
        """
        reach = speed if np.isinf(self.pole) else np.minimum(speed, speed * self.pole + self.outside)
        with np.errstate(divide="ignore"):
            transport = np.where(self.sigma > 0.0, reach / (2.0 * self.sigma), np.inf)
        return np.minimum(self.to_next_radius, transport)

    def rule(self, points: int, speed: float, below: float | None = None) -> ImpactRule:
        r"""The impact rule (:func:`impact_rule`) graded at the projected speed ``speed``, ``points`` per piece.

        With ``below``, only the panels below that orbit coordinate carry
        nodes (a point's rule, which has the point among its ends).
        """
        panels = None if below is None else int(np.searchsorted(self.ends, below))
        return impact_rule(self.ends, self.sigma, points, self.distances(speed), panels)


def outer_amplitude(walls: Walls, partition: ConcentricPartition) -> float:
    r"""The specular amplitude of the body's outer wall, :math:`a` of :class:`ImpactPanels`."""
    return float(walls.on(partition).specular_at(np.array([partition.n_regions]))[0])


@dataclass(frozen=True, eq=False)
class OpticalScale:
    r"""The optical scale a rule over lines is graded from: the thinnest absorbing panel and the normal depth across."""

    thinnest: float
    across: float

    @classmethod
    def of(cls, ends: np.ndarray, sigma: np.ndarray) -> "OpticalScale":
        r"""The scale of the panels ``ends`` of total cross sections ``sigma`` ``(P,)``."""
        optical_width = sigma * np.diff(ends)
        return cls(float(optical_width[optical_width > 0.0].min(initial=np.inf)), float(optical_width.sum()))


def impact_panels(ends: np.ndarray, sigma: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    r"""The impact parameter's panels and their cross sections: a cavity's lines, :math:`b < r_0`, cross a void panel :math:`[0, r_0]`."""
    if ends[0] > 0.0:
        return np.concatenate([[0.0], ends]), np.concatenate([[0.0], sigma])
    return ends, sigma


def polar_rule(scale: OpticalScale, points: int) -> Quadrature1D:
    r"""The cylinder's polar angle :math:`\theta \in [0, \pi/2]`, graded toward grazing and toward the normal (:func:`_grazing_ends`)."""
    return composite_gauss_legendre(_grazing_ends(np.arcsin, scale.thinnest, 2.0 * scale.across), points)


def cosine_rule(scale: OpticalScale, points: int) -> tuple[np.ndarray, np.ndarray]:
    r"""The slab's cosine :math:`\mu \in [-1, 1]`: nodes and weights, each half graded toward grazing and toward the normal (:func:`_grazing_ends`)."""
    half = composite_gauss_legendre(_grazing_ends(np.asarray, scale.thinnest, scale.across), points)
    return np.concatenate([-half.pts[::-1], half.pts]), np.concatenate([half.wts[::-1], half.wts])


def impact_per_polar(impact_at: Callable[[float], ImpactRule], polar: Quadrature1D):
    r"""The cylinder's iterated rule over :math:`(b, \theta)` for :math:`\mathrm db\,\mathrm d\theta`: coordinates ``(L, 2)``,
    weights ``(L,)``, levels.

    For each polar node :math:`\theta_j` the impact rule ``impact_at`` builds
    at that line's own projected speed :math:`\sin\theta_j`, so each polar
    angle is graded toward its own features (:meth:`ImpactPanels.distances`)
    and never at a slower one's (#587). A role's own density over the lines
    multiplies the weights after (:class:`~.assembly.LineRule`,
    :class:`~.reading.PointRule`).
    """
    coordinates, weights, tops, half_chords = [], [], [], []
    for theta, polar_weight in zip(polar.pts, polar.wts, strict=True):
        impact = impact_at(float(np.sin(theta)))
        coordinates.append(np.stack([impact.b, np.full(impact.b.shape, theta)], axis=-1))
        weights.append(impact.weights * polar_weight)
        tops.append(impact.top)
        half_chords.append(impact.half_chord)
    return np.concatenate(coordinates), np.concatenate(weights), (np.concatenate(tops), np.concatenate(half_chords))


@dataclass(frozen=True, eq=False)
class Lines:
    r"""A weighted set of the oriented lines through a body, for one group, on a panel basis.

    Its lines are the chart's :class:`~orpheus.geometry.chart.LineDomain` at
    :attr:`coordinates`, each of weight :attr:`weights`, chorded at their
    exact :attr:`levels` on a radial chart. Two roles hold one, each with its
    own measure: the lines of space (:class:`~.assembly.LineRule`), whose
    weight is the quadrature weight times the invariant density over
    :math:`4\pi`, for the Galerkin block; and the directions at one point
    (:class:`~.reading.PointRule`), whose weight is
    :math:`\mathrm d\Omega/4\pi`, for the reading. Both are graded from
    the group's optical scale by the same one-dimensional rules
    (:mod:`.lines`).

    Attributes
    ----------
    basis:
        The basis.
    walls:
        The body's walls, keyed on the body's partition.
    sigma_t:
        The group's total cross section of each region, ``(n,)``.
    coordinates:
        The lines' coordinates on the chart's line domain, ``(L, k)``.
    weights:
        :math:`w_L`, ``(L,)``.
    levels:
        Each line's exact level, a radius and its half-chord there, each ``(L,)``
        (:meth:`~orpheus.geometry.chart.RadialImage.at_level`): required on a radial chart, ``None`` on the slab,
        whose lines have none.
    chunk:
        The most lines one traversal rule takes.
    budget:
        The most piece slots (:attr:`TraversalRule.extent`) one traversal rule holds; a chunk over it is halved.
        Smaller is faster, the arrays staying in cache, until the chunks shrink to a few lines and the per-rule
        overhead dominates (measured 2026-10-06, white cylinders at 8 points: two regions, 8512 lines, budget
        4096 took 126 s and 7.9 GB, 256 took 32 s and 0.6 GB; one region at tau = 0.01, 15 360 lines graded to
        theta = 3e-7, budget 256 did not finish in 25 minutes, 1024 took 57 s and 1.3 GB, 4096 took 72 s and 3.5 GB,
        the lines ordered by projected speed). Re-measured 2026-10-07 after the attenuated integral stopped
        padding to the batch's thickest stretch (#586), one-region white cylinders at 8 points, 12 along each
        line: at tau = 30, 18 432 lines, budget 512 took 29.4 s and 0.7 GB, 1024 took 27.8 s and 0.8 GB, 4096
        took 31.4 s and 2.0 GB, 16 384 took 37.2 s and 6.0 GB; at tau = 0.01, 512 took 17.5 s and 0.5 GB, 1024
        took 15.6 s and 0.7 GB.
    """

    basis: PanelBasis
    walls: Walls
    sigma_t: np.ndarray
    coordinates: np.ndarray
    weights: np.ndarray
    levels: tuple[np.ndarray, np.ndarray] | None
    chunk: int
    budget: int

    def __post_init__(self) -> None:
        if self.coordinates.shape[:1] != self.weights.shape:
            raise ValueError(f"one weight per line; got {self.coordinates.shape} coordinates and {self.weights.shape} weights")
        if self.chunk < 1 or self.budget < 1:
            raise ValueError(f"a chunk holds at least one line and the budget one piece; got {self.chunk}, {self.budget}")
        radial = self.basis.regions.chart.acts_on_kept_space
        if radial != (self.levels is not None):
            raise ValueError("a radial chart's lines carry their exact levels and a slab's carry none (#590)")
        if self.levels is not None and any(level.shape != self.weights.shape for level in self.levels):
            raise ValueError(f"one level per line; got {[level.shape for level in self.levels]} for {self.weights.shape} weights")

    def ordered(self) -> "Lines":
        r"""The same lines ordered by projected speed: a line's pieces grow as it nears grazing, and a chunk pads every
        line to its longest, so lines of like cost share a chunk."""
        domain = self.basis.regions.chart.line_domain()
        order = np.argsort(self.basis.regions.chart.projected_speed(domain.lines(self.coordinates).direction), kind="stable")
        levels = None if self.levels is None else (self.levels[0][order], self.levels[1][order])
        return replace(self, coordinates=self.coordinates[order], weights=self.weights[order], levels=levels)

    def chunks(self, points: int, inner_points: int) -> Iterator[tuple[TraversalRule, slice]]:
        r"""The traversal rule of each chunk of lines, with its slice, ``points`` and ``inner_points`` along each line.

        A chunk holds at most :attr:`chunk` lines, and one whose rule exceeds
        :attr:`budget` piece slots (:attr:`TraversalRule.extent`) is halved
        until it fits.
        """
        domain = self.basis.regions.chart.line_domain()
        pending = [(k, min(k + self.chunk, len(self.coordinates))) for k in range(0, len(self.coordinates), self.chunk)]
        while pending:
            start, stop = pending.pop(0)
            level = None if self.levels is None else (self.levels[0][start:stop], self.levels[1][start:stop])
            traversals = TraversalRule.of(
                domain.lines(self.coordinates[start:stop]), self.basis, self.walls, self.sigma_t, points, inner_points, level
            )
            if traversals.extent > self.budget and stop - start > 1:
                middle = (start + stop) // 2
                pending[:0] = [(start, middle), (middle, stop)]
                continue
            yield traversals, slice(start, stop)


@dataclass(frozen=True, eq=False)
class StackedSources:
    r"""The sources a traversal rule carries, stacked: the emission functions, then a unit current entering each diffuse wall.

    Each array is ``(..., 2, M + W)``: per traversal of the period and per
    source. The layout is spelled here and read back by :meth:`split`.
    """

    outflow: np.ndarray
    arriving: np.ndarray
    inflow: np.ndarray

    @classmethod
    def of(cls, rule: TraversalRule, columns: np.ndarray, walls: DiffuseWalls) -> "StackedSources":
        r"""The stacked sources of ``rule``'s lines: the basis functions ``columns`` and the walls ``walls``.

        The emission functions leave each traversal as their outflow
        :math:`B_k` and inject nothing; the walls' currents inject
        :math:`1/D` at the traversals entering them and leave nothing. The
        inflow is the period's least solution of both
        (:meth:`~.closure.LinePeriod.inflow`).
        """
        entering = walls.injected(rule.period)                                                 # (..., 2, W)
        depth = rule.optical_depth
        outflow = np.concatenate([rule.outflow()[..., columns], np.zeros(entering.shape)], axis=-1)
        arriving = np.concatenate([np.zeros(depth.shape + (columns.size,)), entering], axis=-1)
        return cls(outflow, arriving, rule.period.inflow(depth, outflow, arriving))

    @staticmethod
    def split(stacked: np.ndarray, emission: int) -> tuple[np.ndarray, np.ndarray]:
        r"""A quantity per stacked source ``(..., M + W)`` as its emission part ``(..., M)`` and its walls' part ``(..., W)``."""
        return stacked[..., :emission], stacked[..., emission:]


__all__ = [
    "ImpactRule", "Lines", "OpticalScale", "StackedSources", "ImpactPanels", "cosine_rule", "impact_panels", "impact_per_polar", "impact_rule",
    "outer_amplitude", "polar_rule",
]
