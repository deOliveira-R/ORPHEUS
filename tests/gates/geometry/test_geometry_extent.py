r"""The geometric coordinate ``GeometryExtent(interval)`` (#405 P1 step 8, S8.8).

Specified by the test-architect (2026-10-02, ``.claude/plans/reference_p1_spec.md``
§1.8). The module is ``orpheus/geometry/extent.py``, exported from
``orpheus.geometry``.

``GeometryExtent(interval=i)`` is the width of interval ``i``, in cm, with
every interval outside it translated outward (ruling 3 of 2026-10-02): on a
bare body interval 0 is the outer size, on a reflected slab it is the core
under a fixed reflector. Step 8 resolves and validates the key only; the
chart zero, the physical value and the admissible range wait for #529
(ruling 6), so no row here moves a breakpoint.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import GeometryExtent
from orpheus.numerics.content import ContentIdentity
from tests.gates._content_identity_helpers import require
from tests.gates.specification._fixtures import slab2, slab3_repeated

pytestmark = pytest.mark.foundation


def test_s8_8_an_extent_is_a_content_value_keyed_by_its_index() -> None:
    a = GeometryExtent(1)
    require(isinstance(a, ContentIdentity), "not a ContentIdentity")
    require(a == GeometryExtent(interval=1) == GeometryExtent(np.int64(1)), "a spelling of index 1 is another key")  # pyright: ignore[reportArgumentType]  # the coercion is the subject
    require(a != GeometryExtent(0) and len({a, GeometryExtent(0)}) == 2, "the index is not content")
    require(type(GeometryExtent(np.int64(1)).interval) is int, "the index is not stored as an int")  # pyright: ignore[reportArgumentType]  # the coercion is the subject


@pytest.mark.parametrize(
    "interval,error,fragment",
    [
        pytest.param(-1, ValueError, r"non-negative, got -1", id="negative"),
        pytest.param(True, TypeError, r"the interval is an int, got bool", id="bool"),  # the shared parse_integer's wording (#559)
        pytest.param(1.0, TypeError, r"the interval is an int, got float", id="float"),
        pytest.param("outer", TypeError, r"the interval is an int, got str", id="str"),
        pytest.param(None, TypeError, r"the interval is an int, got NoneType", id="none"),
    ],
)
def test_s8_8_construction_refusals(interval, error, fragment: str) -> None:
    """``-1`` is refused rather than read as "the last interval": that would be a
    second spelling of one coordinate whose meaning moves with the interval count."""
    with pytest.raises(error, match=fragment):
        GeometryExtent(interval)


def test_s8_8_resolution_against_a_geometry() -> None:
    geometry = slab2()
    for i in (0, 1):
        require(GeometryExtent(i).resolve(geometry) == GeometryExtent(i), f"interval {i} did not resolve to itself")
    with pytest.raises(ValueError, match=r"interval 2 does not exist; the geometry has 2 interval"):
        GeometryExtent(2).resolve(geometry)


def test_s8_8_the_range_is_the_intervals_not_the_materials() -> None:
    """qa F1: ``mat_ids (1, 0, 1)`` is three intervals over two materials; index 2
    resolves and index 3 is refused naming 3 (a ``len(set(mat_ids))`` reads 2)."""
    geometry = slab3_repeated()
    require(GeometryExtent(2).resolve(geometry) == GeometryExtent(2), "the last interval did not resolve")
    with pytest.raises(ValueError, match=r"interval 3 does not exist; the geometry has 3 interval"):
        GeometryExtent(3).resolve(geometry)
