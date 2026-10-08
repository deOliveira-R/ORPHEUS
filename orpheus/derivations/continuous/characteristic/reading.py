r"""The reading at a point: the transported emission's scalar flux there, and its angular flux.

The Galerkin flux :math:`\phi_h` is the L2 projection of the transported
flux :math:`\mathcal K q` (:math:`W\phi_h = Kq` says :math:`(W\phi_h)_i =
\int u_i\,\mathcal K q`), so its value at a point carries the basis's
projection error. The reading transports the converged emission once more,
to the point (the iterated Galerkin solution, the user's ruling of
2026-10-08):

.. math::

    \phi(x) = (\mathcal K q)(x) = \int_{S^2} \psi(x, \Omega)\,\mathrm d\Omega,

with :math:`\psi` the transport of :math:`q` along the line through
:math:`x` in direction :math:`\Omega`, its walls' returns included. It is
the same transport the Galerkin block is assembled from, under a second
test measure: the block integrates :math:`u_i\,\psi` over the lines of
space, the reading integrates :math:`\psi` over the directions at one
point, so :math:`\int u_i\,\phi\,\mathrm dV = (Kq)_i` (the verification
spec's C6).

**The lines through a point are the line rule's.** A line through the point
:math:`x` at orbit coordinate :math:`c` with impact parameter :math:`b \le c`
is congruent, under the chart's group, to the line domain's canonical line
of the same coordinates (:meth:`~orpheus.geometry.chart.LineDomain.lines`),
and :math:`x` lies on it at the two signed positions :math:`\pm\sqrt{c^2 -
b^2}` from its closest approach: one inward, one outward. So the directions
at :math:`x` are integrated over the line domain's own coordinates, by the
line rule's own constructors (:func:`~.lines.impact_rule`,
:func:`~.lines.polar_rule`, :func:`~.lines.cosine_rule`) on the panel
partition with :math:`c` inserted as one more end, keeping the impact
panels below :math:`c` (the user's ruling of 2026-10-08), and the lines are
one :class:`~.lines.Lines`, weighted by the point's measure and held by
the point's role (:class:`PointRule`), as the block's role holds the lines
of space (:class:`~.assembly.LineRule`). The gradings the point needs are
the ones the line rule already makes:

* a **tangency** :math:`b = r_k` is a panel end, under the visibility
  substitution;
* the **grazing direction** at the point is the top panel's end
  :math:`b = c`, where the substitution's variable
  :math:`y = \sqrt{c^2 - b^2}` is :math:`c|\mu|` exactly;
* a point **near a wall** or an interface is the impact rule's grading
  toward :math:`y = 0`, at the next radius out and at the tangency's own
  scales (:func:`~.lines.tangency_distances`; the boundary layer, the
  spec's regime 3); on the slab the point splits its panel, so the thinnest
  optical width that grades the cosine toward 0 includes the point's
  distance to each face.

**The weights are the point's measure,** :math:`\mathrm d\Omega/4\pi` in the
line coordinates: the direction box's density
(:attr:`~orpheus.geometry.chart.DirectionDomain.density`) over :math:`4\pi`,
times the Jacobian from the line coordinates to the box's:

* the sphere: :math:`\mu = \pm\sqrt{c^2 - b^2}/c`, so
  :math:`|\mathrm d\mu| = b\,\mathrm db / (c\sqrt{c^2 - b^2})`, one weight for
  each of the two signs;
* the cylinder: :math:`\sin\alpha = b/c` and :math:`w = \cos\theta`, so
  :math:`\mathrm d\alpha\,\mathrm dw = \mathrm db\,\sin\theta\,\mathrm d\theta
  / \sqrt{c^2 - b^2}`, one weight for each branch :math:`\alpha` and
  :math:`\pi - \alpha`;
* the slab: the cosine itself, read once on each line.

On the singular stratum (the sphere's centre, the cylinder's axis) every
direction has :math:`b = 0`: the rule is the line :math:`b = 0` read at its
closest approach, the direction box a point or the axial cosine.

**Every half-chord is exact.** Each line carries its node's level, the
panel's top and the substitution's own variable :math:`y` there
(:class:`~.lines.ImpactRule`), and is chorded from it (#590): the weight's
:math:`\sqrt{c^2 - b^2}` and the point's parameters
(:meth:`~orpheus.geometry.chart.RadialImage.parameters_at`, the method the
chord forms its crossings by) are
:math:`\sqrt{(c - r_{k+1})(c + r_{k+1}) + y^2}`, so no difference
:math:`c - b` is formed near the grazing direction, and a point on a wall
or an interface is read at that crossing to the last bit. A line built
through the point in a given direction (:func:`angular_flux`) carries the
level :math:`(c, c|\Omega_x|/|P\Omega|)`, exact from the direction.

**The diffuse walls.** The walls that re-emit isotropically return the
currents :math:`j = \alpha(I - T\alpha)^{-1}U^{\mathsf T}q`
(:attr:`~.closure.WallCoupling.currents`), the currents the Galerkin block's
update carries. The reading transports a unit current entering each such
wall to the point, as the block's response does
(:class:`~.lines.StackedSources`), and folds it onto the emission through
:math:`j` (:meth:`~.closure.WallCoupling.on_emission`).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from orpheus.geometry.chart import Chart, DirectionShape, LineShape, half_chord
from orpheus.geometry.line import Line

from .assembly import GroupTransport, TransportResolution
from .basis import PanelBasis
from .closure import FULL_SOLID_ANGLE
from .lines import (
    ImpactRule, Lines, OpticalScale, StackedSources, cosine_rule, impact_by_polar, impact_panels, impact_rule, outer_amplitude,
    polar_rule, tangency_distances,
)
from .transport import TraversalRule
from .walls import Walls

#: Below this orbit coordinate a point's half-chords and weights underflow (their squares are subnormal): read the centre.
_SMALLEST_POINT = float(np.sqrt(np.finfo(float).tiny))

#: The slab's grazing floor: below it a direction's crossing parameters are too large to locate the point on its line.
_SLAB_GRAZING = float(np.sqrt(np.finfo(float).eps))


def refuse_outside(basis: PanelBasis, orbit_coordinate: float) -> float:
    r"""The point's orbit coordinate :math:`c`, refused outside the body and where it is too small to be resolved.

    One door for both readings (:class:`PointRule`, :func:`angular_flux`).
    A positive :math:`c` below :data:`_SMALLEST_POINT` is refused: its
    squares underflow, and the centre itself is the reading to ask for.

    ELEGANCE-DEBT[guard] #582: retires when the kernel's half-chords are scaled so that their squares neither
    underflow nor overflow at extreme radii, so a point near the centre has a resolved measure.
    """
    c = float(orbit_coordinate)
    ends = np.asarray(basis.partition.breakpoints)
    if not ends[0] <= c <= ends[-1]:
        raise ValueError(f"a point is read inside the body, [{ends[0]!r}, {ends[-1]!r}]; got the orbit coordinate {c!r}")
    if 0.0 < c < _SMALLEST_POINT:
        raise ValueError(
            f"a point closer to the centre than {_SMALLEST_POINT:.3g} has half-chords that underflow; got {c!r}: read the centre"
        )
    return c


def _stacked_flux(rule: TraversalRule, t: np.ndarray, transport: GroupTransport) -> np.ndarray:
    r"""The flux at parameters ``t`` ``(L, q)`` on ``rule``'s lines, ``(L, q, M + W)``, :data:`~.closure.FULL_SOLID_ANGLE` times per steradian.

    Per stacked source (:class:`~.lines.StackedSources`): the inflow the
    period returns, attenuated to the point, plus the emission's own
    transport along the point's transit.
    """
    columns = transport.support
    sources = StackedSources.of(rule, columns, transport.coupling.walls)
    reading = rule.at(t)
    walls = np.zeros(t.shape + (transport.coupling.walls.breakpoint.size,))
    return reading.inflow_term(sources.inflow) + np.concatenate([reading.vacuum[..., columns], walls], axis=-1)


def _on_emission(stacked: np.ndarray, transport: GroupTransport) -> np.ndarray:
    """A functional of the stacked sources read on the group's emission (:meth:`~.closure.WallCoupling.on_emission`)."""
    return transport.coupling.on_emission(*StackedSources.split(stacked, transport.support.size))


@dataclass(frozen=True, eq=False)
class PointRule:
    r"""The directions at a point as the line domain's lines through it, weighted by :math:`\mathrm d\Omega/4\pi`.

    Built by :meth:`of`. Each line passes through the point once on each
    :attr:`side` of its closest approach, both readings with the line's
    weight.

    Attributes
    ----------
    lines:
        The lines, weighted by the point's measure (:class:`~.lines.Lines`).
    orbit_coordinate:
        The point's orbit coordinate :math:`c`.
    """

    lines: Lines
    orbit_coordinate: float

    @property
    def side(self) -> np.ndarray:
        r"""The sides of the closest approach each line is read at, ``(q,)``: both, off a radial chart's singular stratum; the approach itself otherwise."""
        chart = self.lines.basis.regions.chart
        on_axis = chart.directions_at(self.orbit_coordinate).on_stratum
        return np.array([-1.0, 1.0]) if chart.acts_on_kept_space and not on_axis else np.array([0.0])

    @classmethod
    def of(cls, basis: PanelBasis, walls: Walls, sigma_t: np.ndarray, orbit_coordinate: float, points: int) -> "PointRule":
        r"""The rule at the point of orbit coordinate ``orbit_coordinate``, ``points`` per piece of every coordinate.

        The lines are chunked more coarsely than the block's
        (:class:`~.assembly.LineRule`'s ``chunk`` and ``budget``): the reading
        holds no Volterra arrays, so its pieces are cheaper and the per-chunk
        overhead dominates sooner (measured 2026-10-08, a point at c = 0.37
        of white cylinders, 8 line points, 12 along each line: one region,
        5440 lines, chunk 512 and budget 1024 took 2.37 s, 4096 and 32 768
        took 1.49 s at 1.25 GB, 8192 and 131 072 took 1.76 s at 2.8 GB; three
        regions, 8512 lines, 15.0 s, 8.1 s and 8.3 s). Those line counts predate
        the grading law at every panel top (:func:`~.lines.tangency_distances`),
        which raised them to 19 584 and 44 992 lines and the readings to about
        5 s and 43 s on a loaded host (the archivist, 2026-10-08): the cost is
        #591's.
        """
        sigma_t = np.asarray(sigma_t, dtype=float)
        c = refuse_outside(basis, orbit_coordinate)
        chart = basis.regions.chart
        directions = chart.directions_at(c)
        ends = np.asarray(basis.partition.breakpoints)
        # the point is one more end: each new panel keeps the cross section of the panel it splits
        split = np.unique(np.append(ends, c))
        sigma = basis.on_panels(sigma_t)[np.searchsorted(ends, split[:-1], side="right") - 1]
        scale = OpticalScale.of(split, sigma)
        per_steradian = directions.density / FULL_SOLID_ANGLE
        amplitude = outer_amplitude(walls, basis.partition)
        levels: tuple[np.ndarray, np.ndarray] | None = None
        match chart.line_domain().shape, directions.shape:
            case LineShape.IMPACT, DirectionShape.WHOLE:
                coordinates, weights = np.zeros((1, 1)), np.array([per_steradian])
                levels = (np.zeros(1), np.zeros(1))
            case LineShape.IMPACT_POLAR, DirectionShape.AXIAL_COSINE:
                polar = polar_rule(scale, points)
                coordinates = np.stack([np.zeros_like(polar.pts), polar.pts], axis=-1)
                weights = per_steradian * np.sin(polar.pts) * polar.wts
                levels = (np.zeros_like(polar.pts), np.zeros_like(polar.pts))
            case LineShape.IMPACT, DirectionShape.COSINE:
                impact = _impact_below(chart, split, sigma, c, amplitude, 1.0, points)
                coordinates, levels = impact.b[:, None], (impact.top, impact.half_chord)
                weights = per_steradian * impact.b * impact.weights / (c * _distance_to_point(impact, c))
            case LineShape.IMPACT_POLAR, DirectionShape.ANGLE_AXIAL:
                polar = polar_rule(scale, points)
                impact = _impact_below(chart, split, sigma, c, amplitude, float(np.sin(polar.pts.min())), points)
                coordinates, weights, levels = impact_by_polar(
                    impact, polar, per_steradian * impact.weights / _distance_to_point(impact, c), np.sin(polar.pts) * polar.wts
                )
            case LineShape.COSINE, DirectionShape.COSINE:
                mu, mu_weights = cosine_rule(scale, points)
                coordinates, weights = mu[:, None], per_steradian * mu_weights
            case pair:
                raise AssertionError(f"unreachable: a chart with lines {pair[0]} and directions {pair[1]}")
        return cls(Lines(basis, walls, sigma_t, coordinates, weights, levels, chunk=4096, budget=32768).ordered(), c)

    def row(self, transport: GroupTransport, resolution: TransportResolution) -> np.ndarray:
        r"""The point's functional on the group's emission coefficients, ``(M,)``: :math:`\phi(x) = \mathrm{row}\cdot q`.

        Along each line the traversal rule is the block's (``resolution``'s
        ``points`` and ``inner_points``, :meth:`~.assembly.LineRule.transport`).
        """
        lines, side, c = self.lines, self.side, np.array([self.orbit_coordinate])
        stacked = np.zeros(transport.support.size + transport.coupling.walls.breakpoint.size)
        for rule, chunk in lines.chunks(resolution.points, resolution.inner_points):
            t = rule.period.chord.image.parameters_at(c, side)
            stacked += np.einsum("l,lqs->s", lines.weights[chunk], _stacked_flux(rule, t, transport))
        return _on_emission(stacked, transport)


def _impact_below(
    chart: Chart, ends: np.ndarray, sigma: np.ndarray, c: float, amplitude: float, slowest: float, points: int
) -> ImpactRule:
    r"""The impact rule over the panels below :math:`c` of the partition ``ends`` with :math:`c` among its ends."""
    b_ends, b_sigma = impact_panels(ends, sigma)
    toward = tangency_distances(chart, b_ends, b_sigma, amplitude, slowest)
    return impact_rule(b_ends, b_sigma, points, toward, panels=int(np.searchsorted(b_ends, c)))


def _distance_to_point(impact: ImpactRule, c: float) -> np.ndarray:
    r""":math:`\sqrt{c^2 - b^2}` at each node, from its level (:func:`~orpheus.geometry.chart.half_chord`, :class:`~.lines.ImpactRule`)."""
    return half_chord(np.asarray(c, dtype=float), impact.top, impact.half_chord)


def angular_flux(
    basis: PanelBasis,
    walls: Walls,
    sigma_t: np.ndarray,
    transport: GroupTransport,
    orbit_coordinate: float,
    coordinates: np.ndarray,
    resolution: TransportResolution,
) -> np.ndarray:
    r"""The angular flux per steradian at the point, in the directions at ``coordinates`` of its direction box, per emission function, ``(..., M)``.

    The directions are :meth:`~orpheus.geometry.chart.Chart.directions_at`'s;
    each is read on the line through the point
    (:meth:`~orpheus.geometry.line.Line.through`), chorded at the level
    :math:`(c, c|\Omega_x|/|P\Omega|)`, exact from the direction, on the
    side of its closest approach the direction faces.
    """
    c = refuse_outside(basis, orbit_coordinate)
    chart = basis.regions.chart
    directions = chart.directions_at(c)
    omega = directions.direction(coordinates).reshape(-1, 3)
    batch = np.shape(coordinates)[:-1]
    lines = Line.through(np.broadcast_to(directions.point, omega.shape), omega)
    if chart.acts_on_kept_space:
        speed = chart.projected_speed(omega)
        level = (np.full(speed.shape, c), c * np.abs(omega[:, 0]) / np.where(speed > 0.0, speed, 1.0))
    else:
        _refuse_slab_grazing(omega[:, 0])
        level = None
    traversals = TraversalRule.of(lines, basis, walls, sigma_t, resolution.points, resolution.inner_points, level)
    # past the closest approach where the direction leaves the axis (Omega_x > 0), before it otherwise
    side = np.sign(omega[:, :1])
    t = traversals.period.chord.image.parameters_at(np.array([c]), side)
    psi = _on_emission(_stacked_flux(traversals, t, transport)[:, 0], transport) / FULL_SOLID_ANGLE
    return psi.reshape(*batch, -1)


def _refuse_slab_grazing(cosine: np.ndarray) -> None:
    r"""Refuse a slab direction within :data:`_SLAB_GRAZING` of grazing.

    The slab's line through the point is parametrised from its foot, so a
    grazing direction's crossings sit at parameters of order the width over
    :math:`|\Omega_x|`, whose spacing outgrows a mean free path: `[M]`
    2026-10-08 (qa, F5), at :math:`|\Omega_x| = 10^{-15}` the angular flux
    was off by 1.6e-3, and at :math:`10^{-100}` the line was parallel and
    read 0 for :math:`q/(4\pi\Sigma)`.

    The mechanism is #585's: a slab line works in absolute positions along
    it, so its crossings are differences of large numbers.

    ELEGANCE-DEBT[guard] #585: retires when a slab line's crossings are formed relative to a point inside the body
    (the chord taking the line's foot there), so a grazing direction is located.
    """
    if np.any(np.abs(cosine) < _SLAB_GRAZING):
        raise ValueError(
            f"a slab direction within {_SLAB_GRAZING:.3g} of grazing has crossing parameters too large to locate the "
            "point on its line (#585)"
        )


__all__ = ["PointRule", "angular_flux", "refuse_outside"]
