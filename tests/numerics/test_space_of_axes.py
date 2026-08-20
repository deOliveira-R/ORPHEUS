r"""``FunctionSpace.of_axes`` — composition, naming, per-axis metric, cone
(campaign 1, CS1 step 2).

Gate ids B1–B8 refer to the CS1 battery of record
(``scratch/cs1_verification_plan.md`` §2); B9–B11 land with step 3b in this
same module.

⚠ Deliberately NOT hosted in ``test_space.py`` / ``test_space_algebra.py``:
those pin the LEGACY densifying ``__mul__`` path
(``test_space_algebra.py`` asserts ``inner_product_weights`` IS the dense
``np.outer``), which CS1 keeps and CS2 retires. Two doctrines, two files,
so the CS2 retirement is a file-level move rather than a diff nobody can
read.
"""

from __future__ import annotations

import subprocess
import sys

import numpy as np
import pytest

from orpheus.numerics.axis import Axis, BasisKind, EnergyAxis
from orpheus.numerics.space import FunctionSpace

pytestmark = pytest.mark.foundation

_EDGES_2G = np.array([1.0e7, 1.0e3, 1.0e-3])


def _require(condition: bool, message: str) -> None:
    """A ``-O``-firing assertion (NOT a bare ``assert``)."""
    if not condition:
        pytest.fail(message)


def _point() -> Axis:
    """The quotient spatial point (counting weight — the density convention)."""
    return Axis("spatial", (1,), kind=BasisKind.NODAL)


def _reachable_arrays(space: FunctionSpace) -> list[np.ndarray]:
    """Every ndarray reachable from the space's own state (B4's walker)."""
    found: list[np.ndarray] = []
    stack: list[object] = list(vars(space).values())
    while stack:
        obj = stack.pop()
        if isinstance(obj, np.ndarray):
            found.append(obj)
        elif isinstance(obj, (tuple, list)):
            stack.extend(obj)
        elif isinstance(obj, Axis):
            stack.extend(vars(obj).values())
    return found


def test_of_axes_shape_is_the_concatenation() -> None:
    """B1 — shape by concatenation, over a RANK-2 axis too.

    The rank-2 member (the harmonic axis's shape, in-fence as a generic
    ``Axis``) is what makes a ``shape + (n,)`` vs ``shape + axis.shape``
    slip visible; a battery of rank-1 axes only cannot see it.
    """
    space = FunctionSpace.of_axes(
        EnergyAxis.synthetic(2),
        Axis("harmonic", (3, 5), kind=BasisKind.MODAL),
        _point(),
    )
    _require(space.shape == (2, 3, 5, 1), f"shape {space.shape} != (2, 3, 5, 1)")
    axes = space.axes
    assert axes is not None
    _require(len(axes) == 3, "the axes record must carry all three factors")


def test_of_axes_name_is_INJECTIVE_on_structural_content() -> None:
    r"""B2 ⭐ — distinct axis tuples ⟹ distinct names.

    Load-bearing because space identity is ``(name, shape)`` until S3 (Q2):
    a NAME COLLISION between different axis tuples makes two different
    spaces compare EQUAL, so the composition guard passes an ill-posed sum
    — the exact defect CS1 exists to make unspellable, reintroduced one
    layer down.

    The population deliberately contains two SAME-SHAPE pairs, where shape
    carries no information and only the name can discriminate:

    * ``synthetic(2)`` vs ``from_grid(<2-group edges>)``  (A5)
    * spatial point with weight ``1.0`` vs weight ``2.0``  (A12 / B9)
    """
    from orpheus.data.energy_grid import EnergyGrid

    tuples: list[tuple[Axis, ...]] = [
        (EnergyAxis.synthetic(2), _point()),
        (EnergyAxis.from_grid(EnergyGrid(_EDGES_2G)), _point()),  # same shapes
        (EnergyAxis.synthetic(2), Axis("spatial", (1,), weights=np.array([2.0]), kind=BasisKind.NODAL)),  # same shapes
        (EnergyAxis.synthetic(3), _point()),
        (EnergyAxis.synthetic(2),),
        (EnergyAxis.synthetic(2), Axis("spatial", (1,), kind=BasisKind.MODAL)),  # kind differs
    ]
    spaces = [FunctionSpace.of_axes(*t) for t in tuples]
    names = {s.name for s in spaces}
    _require(
        len(names) == len(tuples),
        f"name collision: {len(names)} names for {len(tuples)} distinct axis tuples",
    )
    _require(
        len(set(spaces)) == len(tuples),
        "space identity collision — two different axis tuples compare equal",
    )


def test_of_axes_name_is_deterministic_across_processes() -> None:
    """B3 — the derived name does not depend on ``PYTHONHASHSEED``.

    A name built from ``hash(...)`` of a str/tuple is per-process random.
    The in-process leg cannot see that; the subprocess leg can, and it is
    the only leg that can.
    """
    space = FunctionSpace.of_axes(EnergyAxis.synthetic(2), _point())
    twin = FunctionSpace.of_axes(EnergyAxis.synthetic(2), _point())
    _require(space.name == twin.name, "same content, two constructions, one name")

    snippet = (
        "from orpheus.numerics.axis import Axis, BasisKind, EnergyAxis\n"
        "from orpheus.numerics.space import FunctionSpace\n"
        "s = FunctionSpace.of_axes(EnergyAxis.synthetic(2), "
        "Axis('spatial', (1,), kind=BasisKind.NODAL))\n"
        "print(s.name)\n"
    )
    import os

    env = dict(os.environ, PYTHONHASHSEED="271828")
    child = subprocess.run(
        [sys.executable, "-c", snippet], capture_output=True, text=True, env=env,
    )
    _require(child.returncode == 0, f"subprocess failed: {child.stderr}")
    _require(
        child.stdout.strip() == space.name,
        f"name differs across processes: {child.stdout.strip()!r} != {space.name!r} "
        f"— a per-process hash leaked into the identity bridge",
    )


def test_of_axes_never_densifies_the_metric() -> None:
    r"""B4 ⭐ — the factor measures stay per-axis; no outer product is
    materialized.

    Three legs, because each is blind alone:

    1. **exact/structural** — ``space.inner_product_weights is None`` (the
       dense slot is not populated) and no ndarray reachable from the
       space has ``size == prod(space.shape)``;
    2. **memory-shaped** — total reachable ndarray bytes <= the per-axis
       bytes + slack, at a shape with a ``[M]`` 1000x separation
       (``(2000,) x (2000,)``: dense 32 000 000 B vs per-axis 32 000 B);
    3. **behavioural** — the metric still APPLIES correctly at that shape
       (a "never densify" implemented by dropping the metric would pass
       legs 1-2; this is that leg's control).

    ⛔ NOT done by asking a densifier to ``MemoryError``: ``[M]`` a 550 GB
    ``np.multiply.outer`` did not raise, it got the process OOM-KILLED
    (exit 137), which fails the RUN, not the TEST.
    """
    n = 2000
    w_a = np.linspace(0.5, 4.0, n)
    w_b = np.linspace(0.25, 8.0, n)
    space = FunctionSpace.of_axes(
        Axis("a", (n,), weights=w_a, kind=BasisKind.NODAL),
        Axis("b", (n,), weights=w_b, kind=BasisKind.NODAL),
    )
    # Leg 1: structural.
    _require(space.inner_product_weights is None, "dense slot populated")
    full = int(np.prod(space.shape))
    _require(
        all(arr.size != full for arr in _reachable_arrays(space)),
        "an ndarray of the FULL product size is reachable from the space",
    )
    # Leg 2: memory-shaped (per-axis 2 x 16 000 B; dense would be 32 MB).
    total_bytes = sum(arr.nbytes for arr in _reachable_arrays(space))
    _require(
        total_bytes <= 10 * (2 * n * 8),
        f"reachable ndarray bytes {total_bytes} exceed per-axis budget "
        f"— something densified",
    )
    # Leg 3: behavioural control — the metric ACTS at this shape.
    x = np.ones((n, n))
    gx = space.apply_metric(x)
    sample = np.random.default_rng(7).integers(0, n, size=(20, 2))
    for i, j in sample:
        _require(
            bool(gx[i, j] == w_a[i] * w_b[j]),
            f"metric wrong at ({i},{j}): {gx[i, j]} != {w_a[i] * w_b[j]}",
        )


def test_per_axis_metric_equals_an_INDEPENDENT_broadcast_reference() -> None:
    r"""B5 — on a weighted toy, the per-axis metric equals a reference
    built in this test.

    ⚠ The reference is written HERE from an explicit reshape, NOT from
    ``FunctionSpace._broadcast_metric`` — routing both sides through the
    production helper would make this a tautology on the very convention
    (LEADING vs trailing padding) that has already shipped a bug once.
    Non-square ``(3, 4)``, because a square toy cannot see an axis swap.
    Power-of-two weights, so every product is IEEE-exact and the
    comparison is bit-level.
    """
    w_a = np.array([2.0, 4.0, 8.0])
    w_b = np.array([0.5, 1.0, 2.0, 4.0])
    space = FunctionSpace.of_axes(
        Axis("a", (3,), weights=w_a, kind=BasisKind.NODAL),
        Axis("b", (4,), weights=w_b, kind=BasisKind.NODAL),
    )
    rng = np.random.default_rng(42)
    x = rng.standard_normal((3, 4))
    # Independent reference: explicit reshape, each weight on ITS axis.
    reference = w_a.reshape(3, 1) * (w_b.reshape(1, 4) * x)
    _require(
        bool(np.array_equal(space.apply_metric(x), reference)),
        "per-axis metric disagrees with the independent placement reference",
    )
    # Pseudo-inverse roundtrip (exact for power-of-two weights).
    _require(
        bool(np.array_equal(space.apply_inverse_metric(space.apply_metric(x)), x)),
        "inverse metric does not invert the metric",
    )


def test_metric_and_inner_product_agree_on_a_NONSQUARE_axis_space() -> None:
    r"""B6 ⭐ — the inner product equals the independently-built weighted
    pairing on a non-square space.

    The ERR-067 family, one layer up. ``FunctionSpace._diagonal_inner_product``'s
    own Notes record that ``inner_product`` and ``apply_metric`` DIVERGED
    in the tree until 2026-08-04 (trailing vs leading broadcast; ``[M]``
    456 vs 552 on a ``(3,3)`` probe, invisible whenever
    ``w.ndim >= x.ndim``). ``A.H = G^-1 A^T G`` is built from
    ``apply_metric`` while the pairing that judges it comes from
    ``inner_product``: if they disagree, the adjoint identity is false by
    construction and every reciprocity gate downstream is meaningless.

    On the per-axis path the production pairing is SINGLE-SOURCED through
    ``apply_metric`` (one spelling — the divergence is unspellable by
    construction), so the literal pair-agreement check would be a
    tautology; the gate therefore compares against a reference built HERE
    from an explicit reshape, which keeps teeth on the whole path (an
    ``apply_metric`` mutation reddens this gate through the inherited
    pairing).
    """
    w_a = np.array([2.0, 4.0, 8.0])
    w_b = np.array([0.5, 1.0, 2.0, 4.0])
    space = FunctionSpace.of_axes(
        Axis("a", (3,), weights=w_a, kind=BasisKind.NODAL),
        Axis("b", (4,), weights=w_b, kind=BasisKind.NODAL),
    )
    rng = np.random.default_rng(3)
    x = rng.standard_normal((3, 4))
    y = rng.standard_normal((3, 4))
    reference = float(np.sum(w_a.reshape(3, 1) * w_b.reshape(1, 4) * x * y))
    got = space.inner_product(x, y)
    _require(
        bool(np.isclose(got, reference, rtol=1e-13, atol=0.0)),
        f"inner product {got} != independent weighted pairing {reference}",
    )


def test_mul_threads_axes_and_does_not_fabricate_them() -> None:
    """B7 — ``(A * B).axes == A.axes + B.axes``; a legacy space on either
    side leaves the product's ``axes`` ``None``.

    The negative half is the point: inventing an axis for a space that
    never declared one would make ``has_coordinate_cone`` answer for
    spaces that have not been migrated (a false True, in the direction
    that silently ENABLES the step-4 cone consult).

    Third leg — the mixed-product BRIDGE: an axis-built factor's measure
    must ride the legacy product's dense slot, never be silently dropped
    (treating a weighted axis-built factor as Euclidean would be a value
    bug wearing a representation label).
    """
    a = FunctionSpace.of_axes(EnergyAxis.synthetic(2), _point())
    b = FunctionSpace.of_axes(Axis("x", (3,), kind=BasisKind.NODAL))
    product = a * b
    a_axes, b_axes = a.axes, b.axes
    assert a_axes is not None and b_axes is not None
    _require(product.axes == a_axes + b_axes, "axes must thread through *")
    _require(
        product.inner_product_weights is None,
        "an all-axes product must not densify",
    )

    legacy = FunctionSpace("legacy", (5,))
    _require((a * legacy).axes is None, "axes fabricated for a legacy right factor")
    _require((legacy * a).axes is None, "axes fabricated for a legacy left factor")

    # The bridge: weighted axis-built x legacy — the measure survives densely.
    weighted = FunctionSpace.of_axes(
        Axis("w", (2,), weights=np.array([2.0, 4.0]), kind=BasisKind.NODAL)
    )
    mixed = weighted * legacy
    _require(mixed.axes is None, "mixed product must not carry a partial axes record")
    dense = mixed.inner_product_weights
    assert dense is not None
    expected = np.multiply.outer(np.array([2.0, 4.0]), np.ones(5))
    _require(
        bool(np.array_equal(dense, expected)),
        "the axis-borne measure was dropped (or distorted) by the mixed product",
    )


@pytest.mark.parametrize(
    "kinds,expected",
    [
        ((BasisKind.NODAL, BasisKind.NODAL), True),
        ((BasisKind.NODAL, BasisKind.MODAL), False),
        ((BasisKind.MODAL, BasisKind.MODAL), False),
    ],
)
def test_has_coordinate_cone_follows_the_basis_kinds(
    kinds: tuple[BasisKind, BasisKind], expected: bool
) -> None:
    """B8a — all-nodal ⟹ True, any-modal ⟹ False."""
    space = FunctionSpace.of_axes(
        Axis("p", (2,), kind=kinds[0]), Axis("q", (3,), kind=kinds[1])
    )
    _require(
        space.has_coordinate_cone is expected,
        f"kinds {kinds} => {space.has_coordinate_cone}, expected {expected}",
    )


def test_has_coordinate_cone_is_None_on_a_legacy_space() -> None:
    """B8b — ``axes is None`` ⟹ ``None``, the third state.

    Three-valued deliberately: ``False`` means "provably no coordinate
    cone", ``None`` means "not migrated". Collapsing them would make the
    step-4 refusal fire on every legacy space in the tree.
    """
    _require(
        FunctionSpace("legacy", (5,)).has_coordinate_cone is None,
        "a legacy space must answer None (not migrated), never True/False",
    )


def test_dual_of_an_axis_built_space_keeps_the_measure() -> None:
    """The dual carries the SAME metric as the primal (L²-Riesz) — for an
    axis-built primal that means the axes record threads through
    ``dual()``; dropping it would silently strip the measure from every
    adjoint built against the dual."""
    w = np.array([2.0, 4.0])
    space = FunctionSpace.of_axes(
        Axis("a", (2,), weights=w, kind=BasisKind.NODAL)
    )
    dual = space.dual()
    _require(dual.axes == space.axes, "dual() must thread the axes record")
    x = np.array([1.0, 1.0])
    _require(
        bool(np.array_equal(dual.apply_metric(x), w)),
        "the dual's metric must equal the primal's",
    )


def test_axis_built_construction_guards() -> None:
    """The two illegal states of an axis-built space are refused, and
    ``of_axes`` refuses an empty factor list.

    One metric source only (per-axis measures XOR dense weights), and the
    shape must BE the axes' concatenation — both are construction bugs
    when violated, so both raise typed errors.
    """
    ax = Axis("a", (2,), weights=np.array([2.0, 4.0]), kind=BasisKind.NODAL)
    with pytest.raises(ValueError, match="one metric source"):
        FunctionSpace(
            "bad", (2,), inner_product_weights=np.array([1.0, 2.0]), axes=(ax,)
        )
    with pytest.raises(ValueError, match="concatenation"):
        FunctionSpace("bad", (3,), axes=(ax,))
    with pytest.raises(ValueError, match="at least one axis"):
        FunctionSpace.of_axes()
