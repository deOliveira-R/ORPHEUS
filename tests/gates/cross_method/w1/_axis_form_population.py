r"""The axes whose form bytes G03g compares with a record taken before the change (N03, spec 7.4).

Not a test module (no ``test_`` prefix). Every function here RETURNS data and none asserts. Each axis is built through a constructor
that exists on the tree before the unit and after it (``Axis.of_quadrature``, ``Axis.of_basis``, ``Axis.of_measure``), so the one
population is built on both trees: the record is read from it on the unchanged tree, the gate reads the same population on the built one.

The population (``[ON03-12]``): the ``Axis.of_quadrature`` axis of each of the 61 directional rules the shipped factories give (the
17 rules of the A population and the 5-factory sweep of the N02 opening, distinct by name), the ``FromBasis`` axes of the spherical
harmonic (L = 1 to 4), Legendre, interval Legendre and tensor product families, one indicator basis (a family with no Gram), and one
synthetic measure with a ``-0.0`` weight.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np

from orpheus.numerics.axis import Axis

RECORD_DIR = Path(__file__).parent / "axis_forms"
RECORD_FILE = RECORD_DIR / "axis_weights.npz"
STAMP_FILE = RECORD_DIR / "STAMP.txt"


def _rules() -> list[tuple[str, Callable[[], object]]]:
    from orpheus.numerics.quadrature.directional import Quadrature as Q

    cases: list[tuple[str, Callable[[], object]]] = []
    for n in (2, 4, 6, 8, 12, 16, 32, 64):
        cases.append((f"gauss_legendre({n})", lambda n=n: Q.gauss_legendre(n)))
    for o in (3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31):
        cases.append((f"lebedev({o})", lambda o=o: Q.lebedev(o)))
    for s in (2, 4, 6, 8, 10, 12, 14, 16):
        cases.append((f"level_symmetric({s})", lambda s=s: Q.level_symmetric(s)))
    for a in (2, 4, 6, 8):
        for b in (2, 3, 4, 5, 6, 8):
            cases.append((f"product({a},{b})", lambda a=a, b=b: Q.product(a, b)))
    for a in (2, 4, 8):
        for b in (4, 8):
            cases.append((f"folded_product({a},{b})", lambda a=a, b=b: Q.folded_product(a, b)))
    return cases


def _bases() -> list[tuple[str, Callable[[], object]]]:
    from orpheus.geometry.manifold import RealSpace
    from orpheus.numerics.basis.indicator_basis import IndicatorBasis
    from orpheus.numerics.basis.interval_legendre_basis import IntervalLegendreBasis, TensorProductBasis
    from orpheus.numerics.basis.legendre_basis import LegendreBasis
    from orpheus.numerics.basis.spherical_harmonic_basis import SphericalHarmonicBasis

    cases: list[tuple[str, Callable[[], object]]] = []
    for L in (1, 2, 3, 4):
        cases.append((f"basis/spherical_harmonic(L={L})", lambda L=L: SphericalHarmonicBasis(L=L)))
    cases.append(("basis/legendre(L=2)", lambda: LegendreBasis(L=2, axis="z")))
    cases.append(("basis/interval_legendre(n=3)", lambda: IntervalLegendreBasis(3)))
    cases.append(("basis/tensor_product(4)", lambda: TensorProductBasis((IntervalLegendreBasis(4),), "spatial_moment")))
    cases.append(("basis/tensor_product(2,2)", lambda: TensorProductBasis((IntervalLegendreBasis(2),) * 2, "spatial_moment")))
    cases.append(("basis/indicator(2 cells)", lambda: IndicatorBasis((np.array([0.0, 1.0, 2.0]),), RealSpace(1))))
    return cases


def population() -> list[tuple[str, Axis]]:
    """``(key, axis)`` for every axis of the population, in a fixed order."""
    from orpheus.geometry.manifold import RealSpace
    from orpheus.numerics.measure import DiscreteMeasure

    out: list[tuple[str, Axis]] = []
    for name, build in _rules():
        out.append((f"rule/{name}", Axis.of_quadrature(build())))  # type: ignore[arg-type]
    for name, build in _bases():
        out.append((name, Axis.of_basis(build())))  # type: ignore[arg-type]
    negzero = DiscreteMeasure(nodes=np.array([0.0, 0.5, 1.0]), weights=np.array([1.0, -0.0, 2.0]), support=RealSpace(1))
    out.append(("measure/negative-zero weight", Axis.of_measure(negzero, "x")))
    return out
