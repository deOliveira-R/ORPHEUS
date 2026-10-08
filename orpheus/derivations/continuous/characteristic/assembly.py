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
regions that emit in the group. A region that emits nothing in it needs no
column, and must not get one when it is void in the group: a basis function
there would be a source on the lossless lines a mirror traps in the void,
which has no finite flux. A region void in the group that does emit in it
(another group scatters or fissions into it there) does get a column, and
under a mirror the trapped source is refused, because its flux is infinite.
Without the emission, :meth:`LineRule.transport` takes the regions with
:math:`\Sigma_t > 0`; the multigroup system passes each group's emission
support (:class:`~.system.GalerkinSystem`). The rows span every panel: the
flux is read everywhere.

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

from collections.abc import Iterator
from dataclasses import dataclass, replace

import numpy as np

from orpheus.geometry.chart import LineShape
from orpheus.numerics.content import ContentIdentity

from .basis import PanelBasis
from .closure import FULL_SOLID_ANGLE, DiffuseWalls, WallCoupling
from .lines import (
    Lines, OpticalScale, StackedSources, cosine_rule, impact_by_polar, impact_panels, impact_rule, outer_amplitude, polar_rule,
    tangency_distances,
)
from .transport import TraversalRule
from .walls import Walls


@dataclass(frozen=True, eq=False)
class TransportResolution(ContentIdentity):
    r"""The resolution of a transport block: the line rule's points per piece, and the traversal rule's along each line.

    Attributes
    ----------
    line_points:
        The points per piece of every coordinate of the line rule (:meth:`LineRule.of`).
    points:
        The points per arc-length piece of the traversals' integrals (:meth:`LineRule.transport`).
    inner_points:
        The points per piece of the Volterra triangle's inner rule (:meth:`LineRule.transport`).
    """

    line_points: int
    points: int
    inner_points: int

    def __post_init__(self) -> None:
        if min(self.line_points, self.points, self.inner_points) < 1:
            raise ValueError(f"every rule takes at least one point per piece; got {self}")


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
        return self.coupling.on_emission(self.line, self.coupling.response)


@dataclass(frozen=True, eq=False)
class LineRule:
    r"""The lines of space through a body, for one group, weighted for the Galerkin block (:meth:`transport`).

    The role of a :class:`~.lines.Lines` whose weight is the quadrature
    weight times the invariant density over :math:`4\pi`, built by
    :meth:`of`; the directions at a point are the other role
    (:class:`~.reading.PointRule`).

    Attributes
    ----------
    lines:
        The weighted lines.
    """

    lines: Lines

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
        scale = OpticalScale.of(ends, sigma)
        b_ends, b_sigma = impact_panels(ends, sigma)
        amplitude = outer_amplitude(walls, basis.partition)
        levels = None
        match domain.shape:
            case LineShape.IMPACT:
                impact = impact_rule(b_ends, b_sigma, points, tangency_distances(basis.regions.chart, b_ends, b_sigma, amplitude, 1.0))
                coordinates, weights, levels = impact.b[:, None], impact.weights, (impact.top, impact.half_chord)
            case LineShape.IMPACT_POLAR:
                polar = polar_rule(scale, points)
                slowest = float(np.sin(polar.pts.min()))
                impact = impact_rule(b_ends, b_sigma, points, tangency_distances(basis.regions.chart, b_ends, b_sigma, amplitude, slowest))
                coordinates, weights, levels = impact_by_polar(impact, polar, impact.weights, polar.wts)
            case LineShape.COSINE:
                mu, weights = cosine_rule(scale, points)
                coordinates = mu[:, None]
        weights = weights * domain.density(coordinates) / FULL_SOLID_ANGLE
        return cls(Lines(basis, walls, sigma_t, coordinates, weights, levels, chunk, budget).ordered())

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
        lines = self.lines
        basis, sigma_t = lines.basis, lines.sigma_t
        emitting = sigma_t > 0.0 if support is None else np.asarray(support, dtype=bool)
        columns = np.flatnonzero(np.repeat(basis.on_panels(emitting), basis.per_panel))
        diffuse = DiffuseWalls.of(lines.walls, basis.partition)
        walls_at = diffuse.breakpoint
        n, m, d = basis.size, columns.size, walls_at.size
        # the sources stacked: the emission functions, then a unit current entering at each diffuse wall
        response, escape = np.zeros((n, m + d)), np.zeros((m, d))
        line_volterra = np.zeros((n, m))
        injected, lost, reached = np.zeros(d), np.zeros(d), np.zeros((d, d))
        for rule, chunk in lines.chunks(points, inner_points):
            weight = lines.weights[chunk]
            period = rule.period
            depth = rule.optical_depth
            forward = weight[:, None] * (period.present & ~period.reversed)                   # (L, 2)
            sources = StackedSources.of(rule, columns, diffuse)                                 # (L, 2, m + d)
            outflow, arriving, inflow = sources.outflow, sources.arriving, sources.inflow
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
        coupling = WallCoupling(response[:, m:], escape, reached / injected, lost / injected, diffuse)
        return GroupTransport(line, coupling, columns)


__all__ = ["GroupTransport", "LineRule", "TransportResolution"]
