r"""The Galerkin assembly over lines: one group's transport block on the panel basis.

The transport block of one group is the bilinear form

.. math::

    K_{ij} = \int u_i(x)\,\phi[u_j](x)\,\mathrm{d}V,

with :math:`\phi[q]` the scalar flux an isotropic emission density :math:`q`
produces through the body and its walls. Every flux path is a line, so the
volume integral is an integral over the oriented lines of space modulo the
chart's group (:meth:`~orpheus.geometry.chart.Chart.line_domain`), each line
carrying its own transport (:class:`~.transport.TraversalRule`):

.. math::

    K = \sum_L w_L \Bigl(V_L + \sum_{k\ \text{forward}} A_k \otimes \mathrm{in}_k\Bigr)
        + R\,\alpha\,(I - T\alpha)^{-1}U^{\mathsf T},

where :math:`w_L` is the line's quadrature weight times the invariant
density over :math:`4\pi` (so :math:`K` is the scalar-flux operator of an
isotropic emission: on a closed homogeneous body :math:`K\mathbf 1 =
W\mathbf 1/\Sigma_t`), :math:`V_L` the line's Volterra block, :math:`A_k`
the entry response of its forward traversals and :math:`\mathrm{in}_k` the
inflow the specular part of the walls returns
(:meth:`~.closure.LinePeriod.inflow`). The last term is the diffuse walls'
(:class:`~.closure.WallCoupling`).

**The emission support.** The columns of :math:`K` span the panels of the
regions that emit; a void region never emits, and a basis function there
would be a source on the lossless lines a mirror traps in the void, which
has no finite flux. The rows span every panel: the flux is read everywhere.

**The rule over lines** (:class:`LineRule`) is graded from the group's
optical scale, by the law the traversal rule follows along a line: a piece
is no wider than the feature it must resolve. Each rule's grading was found
missing by qa on inputs the first fixed counts did not reach
(``.claude/plans/characteristic_reference_architecture.md``, P1 step (b),
third rung; measured 2026-10-06):

* **the grazing direction** (the slab's cosine, the cylinder's polar angle),
  where a line through a panel of optical width :math:`\tau_P` changes on
  the scale :math:`|P\Omega| \sim \tau_P`: halved toward 0 down to the
  thinnest absorbing panel's :math:`\tau_{\min}`. A fixed 12 halvings left the near-void slab's
  collision probability off by 4e-2 to 5e-1, and plain Gauss in the polar
  angle left the cylinder's off by 3.2e-4 at :math:`\tau = 0.01`;
* **the rim of a radius** (sphere, cylinder): each impact panel
  :math:`[r_k, r_{k+1}]` under the visibility-cone substitution, whose
  variable is the chord half-length :math:`y = \sqrt{r_{k+1}^2 - b^2}`,
  graded in :math:`y` at :math:`2^j` mean free paths of the thickest panel
  the line crosses there (the exponential ends of the traversal rule). An
  ungraded rim left the sphere's surface-to-surface probability off by
  4.6e-3 at :math:`\tau = 100`;
* **a small inner radius**: a line turning in the panel above :math:`r_k`
  integrates the panel's odd modes from :math:`c = b`, an Abel transform
  with a :math:`b^{2m+2}\log b` singularity at :math:`b = 0`, a distance
  :math:`r_k` from the panel; the impact panel is halved toward
  :math:`r_k` until each piece is no wider than that distance (the hp law
  of ERR-099). On the panel touching :math:`b = 0` the turning panel is the
  even one, or a cavity, and the integrand is smooth in :math:`b^2`: one
  split at the middle keeps the substitution's Jacobian off :math:`b = 0`.

The lines are taken in chunks, so a batch's node arrays fit in memory; a
chunk whose rule exceeds the piece budget is halved until it fits.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator
from dataclasses import dataclass
from math import pi

import numpy as np

from orpheus.derivations.common.quadrature import composite_gauss_legendre, gauss_legendre
from orpheus.geometry.chart import LineShape

from .basis import PanelBasis
from .grading import VANISHING_DEPTH, exponential_ends, graded_ends, halvings
from .closure import WallCoupling
from .transport import TraversalRule
from .walls import Walls


#: The full solid angle: a line's weight carries its inverse, so the lines carry fluxes as 4pi times a flux per steradian.
_FULL_SOLID_ANGLE = 4.0 * pi


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


def _impact_rule(ends: np.ndarray, sigma: np.ndarray, points: int) -> tuple[np.ndarray, np.ndarray]:
    r"""The impact-parameter rule on :math:`[0, r_n]`, panel by panel of the panel partition ``ends``.

    On each panel :math:`[r_k, r_{k+1}]` the variable is the chord
    half-length :math:`y = \sqrt{r_{k+1}^2 - b^2}` (the visibility-cone
    substitution, which absorbs the square-root end at :math:`r_{k+1}`),
    graded toward each of its singularities until every piece is no wider
    than its distance to them (hp), and at :math:`2^j` mean free paths of
    the thickest panel from :math:`[r_k, r_{k+1}]` outward:

    * toward :math:`y = 0`, the next radius out, :math:`b = r_{k+2}`, at the
      imaginary :math:`y = \pm i\sqrt{r_{k+2}^2 - r_{k+1}^2}` (its chord's square
      root, just past the panel in :math:`b` but :math:`\sqrt{2 r\,\delta}`
      away in :math:`y`; measured 2026-10-06: on a sphere with a 1e-3 first
      region, Gauss in :math:`b` on the wide middle panel [0.52, 0.88] left
      1e-6 at 8 points, the neighbours at 0.952 and 1.0);
    * toward :math:`b = r_k`, the point :math:`b = 0` at :math:`y = r_{k+1}`: the
      Abel transform of the turning panel's odd modes and the
      substitution's Jacobian :math:`y/b` are singular there.

    The panel touching :math:`b = 0` is smooth in :math:`b^2` (its turning
    panel is even, or a cavity), so its lower half is plain Gauss in
    :math:`b` and its upper half is graded as above. ``sigma`` is each
    panel's total cross section, ``(P,)``.
    """
    pts, wts = [], []
    for k, (lo, hi) in enumerate(zip(ends[:-1], ends[1:])):
        if lo == 0.0:
            half = gauss_legendre(0.0, hi / 2.0, points)
            pts.append(half.pts)
            wts.append(half.wts)
            lo = hi / 2.0
        span = float(np.sqrt((hi - lo) * (hi + lo)))
        beyond = float(np.sqrt((ends[k + 2] - hi) * (ends[k + 2] + hi))) if k + 2 < len(ends) else np.inf
        to_centre = lo * lo / (hi + span)                         # r_{k+1} - span, without cancellation
        y_ends = [
            *graded_ends(0.0, span, True, False, halvings(span, beyond), 0.5),
            *graded_ends(0.0, span, False, True, halvings(span, to_centre), 0.5),
            *exponential_ends(np.array(0.0), np.array(span), np.array(float(np.max(sigma[k:])))),
        ]
        y = composite_gauss_legendre(np.unique(np.clip(y_ends, 0.0, span)), points)
        # b^2 = r_k^2 + (span^2 - y^2): near the lower end b -> r_k without the cancellation of hi^2 - y^2
        b = np.sqrt(lo * lo + (span - y.pts) * (span + y.pts))
        pts.append(b)
        wts.append(y.wts * y.pts / b)
    return np.concatenate(pts), np.concatenate(wts)


@dataclass(frozen=True, eq=False)
class GroupTransport:
    r"""One group's transport block: its line part and its diffuse walls' coupling.

    Attributes
    ----------
    line:
        :math:`\sum_L w_L (V_L + \sum_k A_k \otimes \mathrm{in}_k)`, ``(N, M)``: rows on every panel, columns on the emission support.
    coupling:
        The diffuse walls' coupling, its update ``(N, M)``.
    support:
        The basis indices of the emission support, ``(M,)``: the block's columns.
    """

    line: np.ndarray
    coupling: WallCoupling
    support: np.ndarray

    def __post_init__(self) -> None:
        n, m = self.coupling.response.shape[0], self.coupling.escape.shape[0]
        if self.line.shape != (n, m) or self.support.shape != (m,):
            raise ValueError(
                f"a block of {n} rows over {m} emission columns has line ({n}, {m}) and support ({m},); "
                f"got {self.line.shape} and {self.support.shape}"
            )

    @property
    def block(self) -> np.ndarray:
        r"""The transport block :math:`K`, ``(N, M)``."""
        return self.line + self.coupling.update


@dataclass(frozen=True, eq=False)
class LineRule:
    r"""A quadrature over the oriented lines through a body, for one group, on a panel basis.

    Built by :meth:`of`, graded from the group's optical scale. Its lines are
    the chart's :class:`~orpheus.geometry.chart.LineDomain` at
    :attr:`coordinates`, each of weight :attr:`weights`: the quadrature
    weight times the invariant density over :math:`4\pi`.

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
    chunk: int
    budget: int

    def __post_init__(self) -> None:
        if self.coordinates.shape[:1] != self.weights.shape:
            raise ValueError(f"one weight per line; got {self.coordinates.shape} coordinates and {self.weights.shape} weights")
        if self.chunk < 1 or self.budget < 1:
            raise ValueError(f"a chunk holds at least one line and the budget one piece; got {self.chunk}, {self.budget}")

    @classmethod
    def of(
        cls,
        basis: PanelBasis,
        walls: Walls,
        sigma_t: np.ndarray,
        points: int,
        chunk: int = 512,
        budget: int = 1024,
    ) -> "LineRule":
        r"""The rule for the group of total cross sections ``sigma_t`` (per region), ``points`` per piece of every coordinate."""
        sigma_t = np.asarray(sigma_t, dtype=float)
        domain = basis.regions.chart.line_domain()
        ends = np.asarray(basis.partition.breakpoints)
        sigma = basis.on_panels(sigma_t)
        optical_width = sigma * np.diff(ends)
        thinnest = float(optical_width[optical_width > 0.0].min(initial=np.inf))
        across = float(optical_width.sum())                    # the normal optical depth across the body, once
        # the lines through a cavity, b < r_0, cross no panel there: a void impact panel [0, r_0]
        cavity = ends[0] > 0.0
        b_ends = np.concatenate([[0.0], ends]) if cavity else ends
        b_sigma = np.concatenate([[0.0], sigma]) if cavity else sigma
        match domain.shape:
            case LineShape.IMPACT:
                b, weights = _impact_rule(b_ends, b_sigma, points)
                coordinates = b[:, None]
            case LineShape.IMPACT_POLAR:
                b, b_weights = _impact_rule(b_ends, b_sigma, points)
                polar = composite_gauss_legendre(_grazing_ends(np.arcsin, thinnest, 2.0 * across), points)
                b_grid, theta = np.meshgrid(b, polar.pts, indexing="ij")
                coordinates = np.stack([b_grid.ravel(), theta.ravel()], axis=-1)
                weights = np.outer(b_weights, polar.wts).ravel()
            case LineShape.COSINE:
                half = composite_gauss_legendre(_grazing_ends(np.asarray, thinnest, across), points)
                coordinates = np.concatenate([-half.pts[::-1], half.pts])[:, None]
                weights = np.concatenate([half.wts[::-1], half.wts])
        # ordered by projected speed: a line's pieces grow as it nears grazing, and a chunk pads every line to its
        # longest, so lines of like cost share a chunk
        order = np.argsort(domain.chart.projected_speed(domain.lines(coordinates).direction), kind="stable")
        coordinates, weights = coordinates[order], weights[order]
        return cls(basis, walls, sigma_t, coordinates, weights * domain.density(coordinates) / _FULL_SOLID_ANGLE, chunk, budget)

    def _rules(self, points: int, inner_points: int) -> Iterator[tuple[TraversalRule, np.ndarray]]:
        r"""Each chunk's traversal rule with its line weights, a chunk over the budget halved until it fits."""
        domain = self.basis.regions.chart.line_domain()
        pending = [(k, min(k + self.chunk, len(self.weights))) for k in range(0, len(self.weights), self.chunk)]
        while pending:
            start, stop = pending.pop(0)
            rule = TraversalRule.of(
                domain.lines(self.coordinates[start:stop]), self.basis, self.walls, self.sigma_t, points, inner_points
            )
            if rule.extent > self.budget and stop - start > 1:
                middle = (start + stop) // 2
                pending[:0] = [(start, middle), (middle, stop)]
                continue
            yield rule, self.weights[start:stop]

    def transport(self, points: int, inner_points: int, support: np.ndarray | None = None) -> GroupTransport:
        r"""The group's transport block, ``points`` and ``inner_points`` along each line (:class:`TraversalRule`).

        ``support`` masks the emitting regions, ``(n,)``; by default the regions with :math:`\Sigma_t > 0`.

        **Two tallies.** The line part, :math:`R` and :math:`U` are read on
        each line's forward traversals: a point of the line domain is an
        orbit of oriented lines, and its forward reading is one orientation
        of it (the reversed traversals carry the opposite line's cycle). The
        diffuse walls' transmission and loss are a balance, so they are
        tallied over every traversal of the period, where it closes per
        line: what enters a traversal is absorbed in it
        (:math:`\mathrm{in}\,(1 - e^{-\tau})`, by ``expm1``) or leaves it, and
        what leaves is returned along the line, leaks through a wall that
        does not return it, or reaches a diffuse wall (which returns nothing
        specularly: a wall is never both, :class:`~.walls.Wall`). Every
        traversal counts each orbit twice on the slab, whose cosine rule is
        symmetric in sign, and once on the sphere and the cylinder; the
        fractions are ratios of one tally, so the count cancels. They are
        per current the rule actually injects, so a line rule that misses
        the unit current by its quadrature error moves neither. :math:`R` is
        per NOMINAL unit current, the injection 1/D with D = A_w/4 from the
        chart's one density, so the wall's area is read there and nowhere
        else; the two normalisations differ by that quadrature error
        (measured 2026-10-06 by the elegance review: the forward injected
        tally misses 1 by 2.7e-14 on the cylinder at 8 points, 1.1e-7 at 4).
        Dividing :math:`R` by its own tally instead cancels the area out of
        the block, which the gate on the one density reddens.
        """
        basis, sigma_t = self.basis, self.sigma_t
        emitting = sigma_t > 0.0 if support is None else np.asarray(support, dtype=bool)
        columns = np.flatnonzero(np.repeat(basis.on_panels(emitting), basis.per_panel))
        keyed = self.walls.on(basis.partition)
        diffuse = [w for w in keyed.walls if w.diffuse > 0.0]
        walls_at = np.array([w.breakpoint for w in diffuse], dtype=int)
        # D = diag(A_w / 4): a unit isotropic current entering a wall of area A is 1/(pi A) per steradian,
        # which the lines carry as _FULL_SOLID_ANGLE times that (their weights hold its inverse): 1/D
        area = basis.regions.chart.measure_density(np.asarray(basis.partition.breakpoints)[walls_at])
        quarter_area = pi * area / _FULL_SOLID_ANGLE
        n, m, d = basis.size, columns.size, walls_at.size
        # the sources stacked: the emission functions, then a unit current entering at each diffuse wall
        response, escape = np.zeros((n, m + d)), np.zeros((m, d))
        line_volterra = np.zeros((n, m))
        injected, lost, reached = np.zeros(d), np.zeros(d), np.zeros((d, d))
        for rule, weight in self._rules(points, inner_points):
            period = rule.period
            depth = rule.optical_depth
            forward = weight[:, None] * (period.present & ~period.reversed)                   # (L, 2)
            entering_at = period.present[..., None] & (period.entry_wall[..., None] == walls_at)  # (L, 2, d)
            outflow = np.concatenate([rule.outflow()[..., columns], np.zeros(entering_at.shape)], axis=-1)
            arriving = np.concatenate([np.zeros(depth.shape + (m,)), entering_at / quarter_area], axis=-1)
            inflow = period.inflow(depth, outflow, arriving)                                    # (L, 2, m + d)
            line_volterra += rule.volterra(weight)[:, columns]
            response += np.einsum("lk,lki,lks->is", forward, rule.entry_response(), inflow)
            leaving = np.exp(-depth)[..., None] * inflow[..., :m] + outflow[..., :m]
            exits = forward[..., None] * (period.exit_wall[..., None] == walls_at)              # (L, 2, d)
            escape += np.einsum("lkw,lks->sw", exits, leaving)
            # the injections' balance over every traversal: entering = absorbed + leaked + reaching a diffuse wall
            every = weight[:, None] * period.present                                           # (L, 2)
            carried = inflow[..., m:]                                                           # (L, 2, d)
            at_diffuse = period.exit_wall[..., None] == walls_at                                # (L, 2, d)
            leaks = (1.0 - period.amplitude) * ~at_diffuse.any(axis=-1)
            injected += np.einsum("lk,lkd->d", every, arriving[..., m:])
            lost += np.einsum("lk,lkd->d", every * -np.expm1(-depth), carried)
            lost += np.einsum("lk,lkd->d", every * leaks * np.exp(-depth), carried)
            reached += np.einsum("lkw,lkd->wd", every[..., None] * at_diffuse, np.exp(-depth)[..., None] * carried)
        line = line_volterra + response[:, :m]
        coupling = WallCoupling(
            response[:, m:], escape, reached / injected, lost / injected, np.array([w.diffuse for w in diffuse])
        )
        return GroupTransport(line, coupling, columns)


__all__ = ["GroupTransport", "LineRule"]
