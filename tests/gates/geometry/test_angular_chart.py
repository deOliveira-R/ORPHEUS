r"""Every coordinate system declares its angular chart, once (#405 P1 step 6, S6.18-S6.19).

A function on phase space, :math:`q(r, \mu, \varphi)`, reads its direction
in a chart: :math:`\mu` is the cosine to a polar axis and :math:`\varphi` the
azimuth about it, from a reference direction. The user's ruling of
2026-10-02: the chart is the coordinate system's local frame, declared once
by the coordinate system (:attr:`CoordSystem.angular_chart`). The polar axis
is the slab's normal or the radial direction (frame column 0); the azimuth is
measured from :math:`\hat e_z` on the slab and the cylinder (column 2); the
sphere declares no reference, because none exists (the hairy-ball theorem).

S6.18 gates the declaration on 3 of 3 members. S6.19 ties it to the one
other place the same frame is written, the quadrature's direction columns,
so the chart and the sweeps cannot disagree (X4, one definition): the
column the sweeps read as μ is the declared polar axis; on the cylinder the
column the chart leaves perpendicular to both the polar axis and the
reference is the one the quadrature names ξ = Ω·ê_φ; where the azimuth is
declared unobservable, the rule carries no azimuthal information.
"""

from __future__ import annotations

import numpy as np
import numpy.testing as npt
import pytest

from orpheus.geometry.coord import AngularChart, CoordSystem
from orpheus.numerics.quadrature import Quadrature

pytestmark = pytest.mark.foundation

_DECLARED = {
    CoordSystem.CARTESIAN: AngularChart(polar_axis=0, azimuth_reference=2, azimuth_observable=False),
    CoordSystem.CYLINDRICAL: AngularChart(polar_axis=0, azimuth_reference=2, azimuth_observable=True),
    CoordSystem.SPHERICAL: AngularChart(polar_axis=0, azimuth_reference=None, azimuth_observable=False),
}

# One SN-admitted rule per coordinate system.
_RULES = {
    CoordSystem.CARTESIAN: lambda: Quadrature.gauss_legendre(n_ordinates=8),
    CoordSystem.CYLINDRICAL: lambda: Quadrature.folded_product(4, 8),
    CoordSystem.SPHERICAL: lambda: Quadrature.gauss_legendre(n_ordinates=8),
}


@pytest.mark.parametrize("coord", list(CoordSystem), ids=lambda c: c.name)
def test_s6_18_every_coordinate_system_declares_its_chart(coord: CoordSystem) -> None:
    """3 of 3 members declare the chart; the population is the enum itself,
    so a new member without a declaration reds here (the property's ``match``
    has no fallback arm)."""
    if coord not in _DECLARED:
        pytest.fail(f"{coord.name} has no row in this gate's declared table")
    if coord.angular_chart != _DECLARED[coord]:
        pytest.fail(f"{coord.name} declares {coord.angular_chart}, expected {_DECLARED[coord]}")


def test_s6_18_the_chart_refuses_a_degenerate_frame() -> None:
    with pytest.raises(ValueError, match="polar axis is a frame column"):
        AngularChart(polar_axis=3, azimuth_reference=None, azimuth_observable=False)
    with pytest.raises(ValueError, match="other than the polar axis"):
        AngularChart(polar_axis=0, azimuth_reference=0, azimuth_observable=True)


@pytest.mark.parametrize("coord", list(CoordSystem), ids=lambda c: c.name)
def test_s6_19_a_the_sweeps_read_the_declared_polar_axis(coord: CoordSystem) -> None:
    """The legacy names the sweeps read (``mu_x`` on the slab and sphere,
    ``eta`` on the cylinder) are the declared polar column."""
    quad = _RULES[coord]()
    polar = quad.axis_cosines(coord.angular_chart.polar_axis)
    npt.assert_array_equal(quad.mu_x, polar)
    npt.assert_array_equal(quad.eta, polar)


def test_s6_19_b_the_cylinder_rule_names_the_declared_perpendicular_column() -> None:
    """With ê_polar = ê_r and ê_ref = ê_z, the column left over is ê_⊥ = ê_φ,
    and the quadrature names that column ξ = Ω·ê_φ (``Quadrature.xi``). The
    declared reference is therefore the column the rule does NOT call ξ.
    Swapping the reference to column 1 reds the row.

    Refuted 2026-10-02, FOR the question "does a re-synthesis of the columns
    from (μ, φ) pin the reference?": φ = atan2(Ω·ê_⊥, Ω·ê_ref) rebuilds the two
    columns it was computed from for ANY choice of reference, so that row was
    green under the swapped reference (the step-6 mutation battery)."""
    chart = CoordSystem.CYLINDRICAL.angular_chart
    quad = _RULES[CoordSystem.CYLINDRICAL]()
    assert chart.azimuth_reference is not None  # narrowing only: the cylinder declares one
    (perp,) = {0, 1, 2} - {chart.polar_axis, chart.azimuth_reference}
    npt.assert_array_equal(quad.xi, quad.axis_cosines(perp))


@pytest.mark.parametrize("coord", list(CoordSystem), ids=lambda c: c.name)
def test_s6_19_c_an_unobservable_azimuth_leaves_no_trace_in_the_rule(coord: CoordSystem) -> None:
    """Where the chart declares the azimuth unobservable, the rule's
    orbit-mean cosines off the polar axis are identically 0 (ERR-080: they
    are the orbit mean, not a coordinate); where it is observable, they are
    not."""
    chart = coord.angular_chart
    quad = _RULES[coord]()
    off_polar = [np.asarray(quad.mean_axis_cosine(k)) for k in (0, 1, 2) if k != chart.polar_axis]
    carries_azimuth = any(np.any(col != 0.0) for col in off_polar)
    if carries_azimuth != chart.azimuth_observable:
        pytest.fail(
            f"{coord.name}: the rule carries azimuthal information = {carries_azimuth}, "
            f"the chart declares the azimuth observable = {chart.azimuth_observable}"
        )
