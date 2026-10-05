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
``Chart.image``, ``Chart.orbit_coordinate``). The AST leg reads direct imports in 6 of 6 shapes and misses
indirect ones in 4 of 4 (qa ``p5.py``: a helper module, attribute access through
``orpheus.geometry``, ``from orpheus import geometry``, ``importlib`` by string);
the runtime leg sees every route that reaches the kernel. When a spelling
migrates, its row turns RED (the comparison would be the kernel compared with
itself through a facade: ``retirement-audit`` D.14, D.16), and the migration
commit deletes that row. The file is gone when its last row is.

The draws are probe F1's (``scratch/characteristic_architecture/geometry_census/
_probe_geometry_twins.py``: seed 20261005, 1000 draws, R in [0.5, 3], r in [0, R],
mu in [-1, 1]), where the three old spellings agreed to 7e-14 relative.
"""
from __future__ import annotations

import ast
import importlib
import importlib.util
import inspect

import numpy as np
import pytest

from orpheus.geometry.chart import Chart, RadialImage
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.line import Line

_EPS = np.finfo(float).eps
_KERNEL = ("orpheus.geometry.chord", "orpheus.geometry.line", "orpheus.geometry.chart")
# [M] 2026-10-05, seed 20261005: the largest ratio of |kernel - old| to
# eps R (r + R)/h over the 1000 draws was 1.59 for each old spelling, 1.80 on the
# restructured kernel (scaled |P Omega|). The constant was set at 16 (10x the
# first measurement) and is kept: the row stayed green at 16, and a tolerance is
# never widened without a red (8.9x margin over the larger measurement).
_JOIN_C = 16.0


def _assert_independent(module_name: str) -> None:
    """The precondition: ``module_name`` imports no kernel module (else this row is a tautology; delete it)."""
    module = importlib.import_module(module_name)
    tree = ast.parse(inspect.getsource(module))
    found = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.level > 0:          # a relative import, resolved
            base = importlib.util.resolve_name("." * node.level + (node.module or ""), module.__package__)
            node = ast.ImportFrom(module=base, names=node.names, level=0)
        if isinstance(node, ast.ImportFrom) and node.module:
            names = [node.module] + [f"{node.module}.{a.name}" for a in node.names]
        elif isinstance(node, ast.Import):
            names = [a.name for a in node.names]
        else:
            continue
        found += [n for n in names if any(n == k or n.startswith(k + ".") for k in _KERNEL)]
        if isinstance(node, ast.ImportFrom) and node.module == "orpheus.geometry":
            found += [a.name for a in node.names if a.name in {"Chart", "Line", "ConcentricPartition", "Chord", "Crossings"}]
    if found:
        raise AssertionError(
            f"{module_name} now imports the kernel ({found}): this corroboration row compares the kernel "
            f"with itself; delete it in the migration commit (retirement-audit D.14)"
        )


def test_the_independence_precondition_sees_a_kernel_import() -> None:
    """Positive control for the AST precondition (X1), one per import shape: absolute, and relative.

    ``orpheus.geometry.chord`` imports the line and the chart absolutely; a
    module written with ``from .chord import ...`` is the relative shape,
    checked on a synthetic source.
    """
    with pytest.raises(AssertionError, match="this corroboration row compares the kernel"):
        _assert_independent("orpheus.geometry.chord")
    node = ast.parse("from .chord import ConcentricPartition").body[0]
    assert isinstance(node, ast.ImportFrom)
    assert importlib.util.resolve_name("." * node.level + (node.module or ""), "orpheus.geometry") == "orpheus.geometry.chord"


_ENTRY_POINTS = (
    ("ConcentricPartition", "chord"),
    ("ConcentricPartition", "region_containing"),
    ("Chart", "image"),
    ("Chart", "orbit_coordinate"),
)


@pytest.fixture
def kernel_calls(monkeypatch: pytest.MonkeyPatch) -> dict[str, int]:
    """A counting spy on the kernel's entry points; the dict maps ``Class.method`` to its call count."""
    import orpheus.geometry.chart as chart_module
    import orpheus.geometry.chord as chord_module

    owners = {"ConcentricPartition": chord_module.ConcentricPartition, "Chart": chart_module.Chart}
    counts: dict[str, int] = {}
    for owner, name in _ENTRY_POINTS:
        cls = owners[owner]
        original = getattr(cls, name)
        key = f"{owner}.{name}"
        counts[key] = 0

        def counted(*args, _original=original, _key=key, **kwargs):
            counts[_key] += 1
            return _original(*args, **kwargs)

        monkeypatch.setattr(cls, name, counted)
    return counts


def _without_kernel(counts: dict[str, int], old, *args):
    """Evaluate the old spelling and require that it reached no kernel entry point (the runtime independence leg)."""
    before = dict(counts)
    value = old(*args)
    reached = {k: counts[k] - before[k] for k in counts if counts[k] != before[k]}
    if reached:
        raise AssertionError(
            f"the old spelling called the kernel ({reached}): this corroboration row compares the kernel "
            f"with itself; delete it in the migration commit (retirement-audit D.14)"
        )
    return value


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
    ["variant_alpha", "peierls_rho_max", "moc_root"],
)
def test_the_exit_distance_agrees_with_each_old_spelling(spelling: str, kernel_calls: dict[str, int]) -> None:
    """F1 joined: the kernel's exit distance against Variant-alpha ``L_back``, Peierls ``rho_max``, the MoC root.

    Per-draw band ``16 eps R (r + R)/h`` (the problem's conditioning at
    tangency, spec §8).
    """
    if spelling == "variant_alpha":
        _assert_independent("orpheus.derivations.continuous.trajectory_resolvent.chord_oracle")
        from orpheus.derivations.continuous.trajectory_resolvent.chord_oracle import _trajectory_segments_oracle

        def old(R, r, mu):
            return _trajectory_segments_oracle(r, mu, R, np.array([R]))[1]
    elif spelling == "peierls_rho_max":
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


def test_the_backward_segments_agree_with_variant_alpha(kernel_calls: dict[str, int]) -> None:
    """The multi-region backward first leg: region sequence equal, lengths within 1e-13, over 1000 draws.

    ``_trajectory_segments_oracle`` locates each segment's midpoint; the kernel
    takes regions from the crossing order and splits the region of closest
    approach into two slots, merged here. ``[M]`` 2026-10-05: 0 of 1000
    sequences differ, largest length difference 5.1e-15.
    """
    _assert_independent("orpheus.derivations.continuous.trajectory_resolvent.chord_oracle")
    from orpheus.derivations.continuous.trajectory_resolvent.chord_oracle import _trajectory_segments_oracle

    radii = np.array([0.3, 1.1, 2.0])
    part = ConcentricPartition(Chart(CoordSystem.SPHERICAL), (0.0, *radii))
    rng = np.random.default_rng(20261005)
    for _ in range(1000):
        r, mu = rng.uniform(0, 2.0), rng.uniform(-1, 1)
        segments, _ = _without_kernel(kernel_calls, _trajectory_segments_oracle, r, mu, 2.0, radii)
        p = np.array([r, 0.0, 0.0])
        line = Line.through(p, np.array([-mu, np.sqrt(1.0 - mu * mu), 0.0]))
        ch = part.chord(line)
        beyond = ch.lengths_beyond(line.parameter_of(p))
        regions, lengths = [], []
        for region, length in zip(ch.slot_region, beyond):
            if length <= 1e-12:
                continue
            if regions and regions[-1] == region:
                lengths[-1] += length
            else:
                regions.append(int(region))
                lengths.append(float(length))
        assert regions == [s[2] for s in segments]
        np.testing.assert_allclose(lengths, [s[1] - s[0] for s in segments], rtol=0, atol=1e-13)


def test_the_bare_locator_agrees_with_the_inner_owns_spellings_and_not_with_peierls(kernel_calls: dict[str, int]) -> None:
    """F3 joined: at every breakpoint the kernel answers as the five inner-owns spellings, and unlike Peierls ``which_annulus``.

    Peierls is outer-biased by its own docstring; the disagreement is asserted
    so that its migration onto the kernel shows up here as a changed answer.
    """
    for name in ("orpheus.derivations.continuous.trajectory_resolvent.chord_oracle", "orpheus.moc.geometry",
                 "orpheus.mc.solver", "orpheus.derivations.continuous.peierls_nystrom.geometry"):
        _assert_independent(name)
    from orpheus.derivations.continuous.peierls_nystrom.geometry import SPHERE_1D
    from orpheus.derivations.continuous.trajectory_resolvent.chord_oracle import (
        _region_at_radius_cyl,
        _region_at_radius_oracle,
    )
    from orpheus.mc.solver import ConcentricPinCell
    from orpheus.moc.geometry import _identify_region

    radii = np.array([0.3, 1.1, 2.0])
    part = ConcentricPartition(Chart(CoordSystem.SPHERICAL), (0.0, *radii))
    pin = ConcentricPinCell(radii=list(radii), mat_ids=[0, 1, 2], pitch=10.0)
    for r in (0.3, 1.1, 0.7, 0.05, 1.9):
        def old_spellings(r=r):
            return (_region_at_radius_oracle(r, radii), _region_at_radius_cyl(r, radii),
                    _identify_region(r, 0.0, 0.0, 0.0, radii, 3), pin.material_id_at(5.0 + r, 5.0),
                    SPHERE_1D.which_annulus(r, radii))
        va, va_cyl, moc, mc, peierls = _without_kernel(kernel_calls, old_spellings)
        kernel = int(part.region_containing(np.array(r)))
        assert kernel == va == va_cyl == moc == mc
        on_breakpoint = r in (0.3, 1.1)
        assert peierls == (kernel + 1 if on_breakpoint else kernel)
