"""L4 corroboration: the geometric kernel against today's independent spellings, during the migration.

**What this is.** Code-to-code agreement (``vv-principles``: L4, no correctness
content). The kernel's correctness is ``test_chart.py``, ``test_line.py``,
``test_chord.py`` and ``test_line_measure.py``. This file's value is migration
safety: when a consumer moves onto the kernel, the pair (its old spelling, the
kernel) at its own inputs shows whether the move changed an answer. Ruled
2026-10-05 (``.claude/plans/characteristic_reference_architecture.md``,
NEEDS 10): a tracked L4 harness, deleted row by row in the migration commits.
It carries no ``l``-level marker (``pyproject.toml`` defines l0-l3 and
foundation only).

**Its retirement is built in.** Before comparing, each row asserts by AST
that the old spelling's module does not import the kernel, and, at run time,
that evaluating the old spelling makes 0 calls into the kernel's entry points
(``ConcentricPartition.chord``, ``ConcentricPartition.region_containing``,
``Chart.image``, ``Chart.orbit_coordinate``). Both legs live in
``tests/gates/_corroboration.py``, shared, until step (e2), with the characteristic reference's
corroboration file (deleted with the old family). When this file landed the AST leg read direct imports in 6
of 6 shapes and missed indirect ones in 4 of 4 (qa ``p5.py``: a helper module,
attribute access through ``orpheus.geometry``, ``from orpheus import geometry``,
``importlib`` by string); the shared leg now also reads attribute chains on an
imported name (and its literal ``getattr`` spelling) and module names in strings; the helper module's shape
belonged to a transitive leg this file never called (deleted in step (e2)). The runtime leg counts every code object of the kernel modules that runs
(``sys.monitoring``), so it sees every in-process route into the kernel. When a spelling
migrates, its row turns RED (the comparison would be the kernel compared with
itself through a facade: ``retirement-audit`` D.14, D.16), and the migration
commit deletes that row. The file is gone when its last row is.

The draws are probe F1's (``scratch/characteristic_architecture/geometry_census/
_probe_geometry_twins.py``: seed 20261005, 1000 draws, R in [0.5, 3], r in [0, R],
mu in [-1, 1]), where the three old spellings agreed to 7e-14 relative.
"""
from __future__ import annotations

from collections.abc import Iterator

import numpy as np
import pytest

from orpheus.geometry.chart import Chart, RadialImage
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.line import Line
from tests.gates import _corroboration as corroboration

_EPS = np.finfo(float).eps
#: The kernel, and the entry points the runtime leg counts (``tests/gates/_corroboration.py`` holds both legs).
KERNEL = corroboration.NewSide(
    modules=("orpheus.geometry.chord", "orpheus.geometry.line", "orpheus.geometry.chart"),
    entry_points=(
        "ConcentricPartition.chord",
        "ConcentricPartition.region_containing",
        "Chart.image",
        "Chart.orbit_coordinate",
    ),
    label="the kernel",
)
# [M] 2026-10-05, seed 20261005: the largest ratio of |kernel - old| to
# eps R (r + R)/h over the 1000 draws was 1.59 for each old spelling, 1.80 on the
# restructured kernel (scaled |P Omega|). The constant was set at 16 (10x the
# first measurement) and is kept: the row stayed green at 16, and a tolerance is
# never widened without a red (8.9x margin over the larger measurement).
_JOIN_C = 16.0


def _assert_independent(module_name: str) -> None:
    """The precondition: ``module_name`` imports no kernel module (else this row is a tautology; delete it)."""
    corroboration.assert_independent([module_name], KERNEL)


def test_the_independence_precondition_sees_a_kernel_import() -> None:
    """Positive control for the AST precondition (X1), one per import shape: absolute, and relative.

    ``orpheus.geometry.chord`` imports the line and the chart absolutely; a
    module written with ``from .chord import ...`` is the relative shape,
    checked on a synthetic source through the same leg.
    """
    with pytest.raises(AssertionError, match="this corroboration row compares the kernel"):
        _assert_independent("orpheus.geometry.chord")
    with pytest.raises(AssertionError, match=r"now imports the kernel \(\{'import': \['orpheus\.geometry\.chord'"):
        corroboration.assert_source_independent(
            "a module of orpheus.geometry", "from .chord import ConcentricPartition", "orpheus.geometry", KERNEL,
        )


@pytest.fixture
def kernel_calls(monkeypatch: pytest.MonkeyPatch) -> Iterator[dict[str, int]]:
    """A counting spy on the kernel's entry points and code; the dict maps ``Class.method`` (or a function) to its call count."""
    with corroboration.spy(monkeypatch, KERNEL) as counts:
        yield counts


def _without_kernel(counts: dict[str, int], old, *args):
    """Evaluate the old spelling and require that it reached no kernel entry point (the runtime independence leg)."""
    return corroboration.without(counts, KERNEL, old, *args)


def test_the_runtime_independence_leg_sees_a_kernel_call(kernel_calls: dict[str, int]) -> None:
    """Positive control for the runtime leg (X1): an "old spelling" that routes through the kernel is refused.

    The stand-in reaches the kernel the way qa's indirect shapes do (no import
    statement naming the kernel module in the caller), through an attribute of
    the ``orpheus.geometry`` package. A kernel call made OUTSIDE the leg is not
    counted against it.
    """
    import orpheus.geometry.chord as chord_module

    def routed_old_spelling(R: float, r: float, mu: float) -> float:
        part = chord_module.ConcentricPartition(Chart(CoordSystem.SPHERICAL), (0.0, R))
        return float(np.sum(part.chord(Line.through(np.array([r, 0.0, 0.0]), np.array([-mu, np.sqrt(1 - mu * mu), 0.0]))).slot_length))

    with pytest.raises(AssertionError, match="the old spelling called the kernel"):
        _without_kernel(kernel_calls, routed_old_spelling, 2.0, 0.5, 0.3)
    _kernel_exit_distance(2.0, 0.5, 0.3)
    assert kernel_calls["ConcentricPartition.chord"] >= 2
    assert _without_kernel(kernel_calls, lambda x: x + 1.0, 1.0) == 2.0


def _draws():
    rng = np.random.default_rng(20261005)
    for _ in range(1000):
        R = rng.uniform(0.5, 3.0)
        r = rng.uniform(0, R)
        mu = rng.uniform(-1, 1)
        yield R, r, mu


def _kernel_exit_distance(R: float, r: float, mu: float) -> tuple[float, float]:
    """The backward first leg from ``(r, 0, 0)`` with ``r(s)^2 = r^2 - 2 r s mu + s^2``, and its half-chord ``h``."""
    p = np.array([r, 0.0, 0.0])
    line = Line.through(p, np.array([-mu, np.sqrt(1.0 - mu * mu), 0.0]))
    ch = ConcentricPartition(Chart(CoordSystem.SPHERICAL), (0.0, R)).chord(line)
    if not isinstance(ch.image, RadialImage):
        pytest.fail("a sphere's chord has a radial image")
    b = float(ch.image.impact_parameter)
    return float(np.sum(ch.lengths_beyond(line.parameter_of(p)))), float(np.sqrt((R - b) * (R + b)))


@pytest.mark.parametrize(
    "spelling",
    ["peierls_rho_max", "moc_root"],
)
def test_the_exit_distance_agrees_with_each_old_spelling(spelling: str, kernel_calls: dict[str, int]) -> None:
    """F1 joined: the kernel's exit distance against Peierls ``rho_max`` and the MoC root.

    The Variant-alpha ``L_back`` arm was deleted with the trajectory-resolvent
    family in step (e2) of the characteristic-reference campaign.

    Per-draw band ``16 eps R (r + R)/h`` (the problem's conditioning at
    tangency, spec §8).
    """
    if spelling == "peierls_rho_max":
        _assert_independent("orpheus.derivations.continuous.peierls_nystrom.geometry")
        from orpheus.derivations.continuous.peierls_nystrom.geometry import SPHERE_1D

        def old(R, r, mu):
            return SPHERE_1D.rho_max(r, -mu, R)
    else:
        _assert_independent("orpheus.moc.geometry")
        from orpheus.moc.geometry import _ray_circle_intersections

        def old(R, r, mu):
            return max(_ray_circle_intersections(r, 0.0, -mu, np.sqrt(1 - mu * mu), 0.0, 0.0, R))

    worst = 0.0
    for R, r, mu in _draws():
        got, h = _kernel_exit_distance(R, r, mu)
        worst = max(worst, abs(got - _without_kernel(kernel_calls, old, R, r, mu)) / (_EPS * R * (r + R) / h))
    assert worst <= _JOIN_C, worst


def test_the_bare_locator_agrees_with_the_inner_owns_spellings_and_not_with_peierls(kernel_calls: dict[str, int]) -> None:
    """F3 joined: at every breakpoint the kernel answers as the three inner-owns spellings, and unlike Peierls ``which_annulus``.

    Until step (e2) two more inner-owns spellings, the trajectory resolvent's
    ``_region_at_radius_oracle`` and ``_region_at_radius_cyl``, were compared
    here; they were deleted with their family.

    Peierls is outer-biased by its own docstring; the disagreement is asserted
    so that its migration onto the kernel shows up here as a changed answer.
    """
    for name in ("orpheus.moc.geometry", "orpheus.mc.solver", "orpheus.derivations.continuous.peierls_nystrom.geometry"):
        _assert_independent(name)
    from orpheus.derivations.continuous.peierls_nystrom.geometry import SPHERE_1D
    from orpheus.mc.solver import ConcentricPinCell
    from orpheus.moc.geometry import _identify_region

    radii = np.array([0.3, 1.1, 2.0])
    part = ConcentricPartition(Chart(CoordSystem.SPHERICAL), (0.0, *radii))
    pin = ConcentricPinCell(radii=list(radii), mat_ids=[0, 1, 2], pitch=10.0)
    for r in (0.3, 1.1, 0.7, 0.05, 1.9):
        def old_spellings(r=r):
            return (_identify_region(r, 0.0, 0.0, 0.0, radii, 3), pin.material_id_at(5.0 + r, 5.0),
                    SPHERE_1D.which_annulus(r, radii))
        moc, mc, peierls = _without_kernel(kernel_calls, old_spellings)
        kernel = int(part.region_containing(np.array(r)))
        assert kernel == moc == mc
        on_breakpoint = r in (0.3, 1.1)
        assert peierls == (kernel + 1 if on_breakpoint else kernel)
