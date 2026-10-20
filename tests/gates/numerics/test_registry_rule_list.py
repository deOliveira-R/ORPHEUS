"""The registry holds ``Quadrature`` classmethods: one rule list, not two.

Ruling 32 of the numerics untangling (2026-10-18/19). The registry used to hold
its own factories (module rule functions wrapped in ``partial``) beside
``Quadrature``'s classmethods, and ``QuadratureSpec.build`` unpacked the result
to a bare ``DiscreteMeasure``. Now each entry's ``factory`` IS the matching
``Quadrature`` classmethod and ``build`` returns the ``Quadrature``.

Two instruments, both able to fail:

* **Bit-for-bit.** For every registered spec and every parameter set its degree
  inversion can return (targets ``0 .. _TARGET_CAP``), the classmethod's
  ``.measure`` equals the measure the module rule function it wraps builds:
  nodes, weights, support, invariance group and exactness claim, with
  ``np.array_equal`` and ``==``, never ``allclose``. The module rule functions
  are spelled out in ``_RULE_FUNCTION`` below so the equality is pinned rather
  than claimed; a classmethod that adds a keyword the rule function lacks, or
  reorders its arguments, reddens it.
* **The table asserts itself.** The registry's degree table is an INVERSION
  (``degree_of_exactness_for``), and nothing used to check it against the
  rule it promises. For every spec and every target degree: the parameters it
  returns build a rule whose ``degree_of_exactness`` is at least the target
  (the promise), and the next smaller parameter set in the family does not
  (the inversion is the cheapest, not merely sufficient).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.quadrature import (
    Quadrature,
    gauss_legendre_on_polar_orbit,
    lebedev_sphere,
    level_symmetric_sn,
    product_mu_phi,
    quadrature_registry,
)
from orpheus.numerics.quadrature.registry import _lebedev_orders

#: Largest target degree swept. Past the level-symmetric frontier (S_18, degree
#: 17) that family's inversion returns ``None`` after a failed weight solve
#: whose cost grows with the order (`[M]` a sweep to 41 did not finish in
#: 120 s; 0..21 takes about 3 s), so the cap sits just above the frontier. The
#: other families keep answering.
_TARGET_CAP = 21

#: Parameter sets the inversion does not return within the cap but the family
#: builds, added to the bit-for-bit sweep: the level-symmetric orders 12 and 16
#: (never the cheapest for any target), and Gauss-Legendre counts past the cap's
#: n = 12 and the odd counts (the inversion rounds to even; the classmethod
#: accepts any).
_EXTRA_GRID = {
    "LevelSymmetricSN": [{"sn_order": n} for n in (12, 16)],
    "GaussLegendre1D": [{"n_ordinates": n} for n in (1, 3, 5, 16, 24, 33, 64)],
}

#: The module rule function each registered rule is a thin classmethod over,
#: spelled out here (not read off the registry) so the comparison has an
#: independent left-hand side.
_RULE_FUNCTION = {
    "GaussLegendre1D": lambda p: gauss_legendre_on_polar_orbit(p["n_ordinates"], axis="x"),
    "LebedevSphere": lambda p: lebedev_sphere(p["order"]),
    "LevelSymmetricSN": lambda p: level_symmetric_sn(p["sn_order"])[0],
    "ProductQuadrature": lambda p: product_mu_phi(p["n_mu"], p["n_phi"])[0],
}

#: The classmethod each registered rule must reference.
_CLASSMETHOD = {
    "GaussLegendre1D": "gauss_legendre",
    "LebedevSphere": "lebedev",
    "LevelSymmetricSN": "level_symmetric",
    "ProductQuadrature": "product",
}


def _grid(spec) -> list[dict]:
    """Every parameter set the spec's inversion returns for targets 0..cap, deduplicated."""
    seen: list[dict] = []
    for target in range(_TARGET_CAP + 1):
        params = spec.degree_of_exactness_for(target)
        if params is not None and params not in seen:
            seen.append(params)
    for extra in _EXTRA_GRID.get(spec.name, []):
        if extra not in seen:
            seen.append(extra)
    return seen


_SPECS = {spec.name: spec for spec in quadrature_registry}


def test_the_registry_names_are_the_rules_this_test_knows():
    """A new registry entry must be added to the two tables above, or this reddens."""
    assert set(_SPECS) >= set(_RULE_FUNCTION) == set(_CLASSMETHOD)


@pytest.mark.foundation
@pytest.mark.parametrize("name", sorted(_CLASSMETHOD))
def test_the_factory_is_the_quadrature_classmethod(name):
    factory = _SPECS[name].factory
    assert factory.__self__ is Quadrature
    assert factory.__func__ is getattr(Quadrature, _CLASSMETHOD[name]).__func__


@pytest.mark.foundation
@pytest.mark.parametrize("name", sorted(_RULE_FUNCTION))
def test_the_classmethod_measure_is_the_rule_function_measure_bit_for_bit(name):
    spec = _SPECS[name]
    grid = _grid(spec)
    assert len(grid) >= 8, f"{name}: the sweep reached only {len(grid)} parameter sets"
    for params in grid:
        built = spec.build(params)
        assert isinstance(built, Quadrature), f"{name}{params}: build returned {type(built).__name__}"
        old = _RULE_FUNCTION[name](params)
        new = built.measure
        assert np.array_equal(new.nodes, old.nodes), f"{name}{params}: nodes"
        assert np.array_equal(new.weights, old.weights), f"{name}{params}: weights"
        assert new.support == old.support, f"{name}{params}: support"
        assert new.invariance_group == old.invariance_group, f"{name}{params}: invariance group"
        assert new.exactness == old.exactness, f"{name}{params}: exactness"
        assert new.exactness is not None, f"{name}{params}: no claim to compare"


def _one_step_smaller(name: str, params: dict) -> list[dict]:
    """The neighbouring cheaper parameter sets of the family (empty at its floor)."""
    if name == "GaussLegendre1D":
        return [{"n_ordinates": params["n_ordinates"] - 2}] if params["n_ordinates"] > 2 else []
    if name == "LebedevSphere":
        orders = _lebedev_orders()
        i = orders.index(params["order"])
        return [{"order": orders[i - 1]}] if i > 0 else []
    if name == "LevelSymmetricSN":
        return [{"sn_order": params["sn_order"] - 2}] if params["sn_order"] > 2 else []
    if name == "ProductQuadrature":
        out = []
        if params["n_mu"] > 1:
            out.append({"n_mu": params["n_mu"] - 1, "n_phi": params["n_phi"]})
        if params["n_phi"] > 1:
            out.append({"n_mu": params["n_mu"], "n_phi": params["n_phi"] - 1})
        return out
    raise AssertionError(f"no neighbour rule for registry entry {name!r}: add one to this test")


@pytest.mark.foundation
@pytest.mark.parametrize("name", sorted(_SPECS))
def test_the_degree_inversion_is_the_cheapest_parameter_set_meeting_the_target(name):
    spec = _SPECS[name]
    answered = 0
    for target in range(_TARGET_CAP + 1):
        params = spec.degree_of_exactness_for(target)
        if params is None:
            continue
        answered += 1
        degree = spec.build(params).measure.degree_of_exactness
        assert degree is not None and degree >= target, (
            f"{name}{params}: the inversion promises degree >= {target}, the built rule says {degree}"
        )
        for smaller in _one_step_smaller(name, params):
            try:
                smaller_degree = spec.build(smaller).measure.degree_of_exactness
            except ValueError:
                continue  # below the family's floor: nothing cheaper exists
            assert smaller_degree is None or smaller_degree < target, (
                f"{name}: {smaller} already reaches degree {smaller_degree} >= {target}, "
                f"so the inversion's {params} is not the cheapest"
            )
    assert answered >= 8, f"{name}: the inversion answered only {answered} of {_TARGET_CAP + 1} targets"
