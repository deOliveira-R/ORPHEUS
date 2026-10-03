r"""Rows to append to ``tests/gates/numerics/test_mesh_free_function.py`` (#405 P1 step 8, S8.11).

Specified by the test-architect (2026-10-02). ``Symbolic.depends_on(*coordinates)``
is the one simplify-decided dependence predicate; ``is_isotropic`` is
``not depends_on(mu, phi)``, the specification's sphere rule reads
``depends_on(phi)`` and its infinite-medium rule ``depends_on(r, mu, phi)``.
(The method's name is the main agent's to choose; the rows name it once, here.)
"""

from __future__ import annotations

import pytest
import sympy as sp

from orpheus.numerics.mesh_free_function import Symbolic
from tests.gates._content_identity_helpers import require

pytestmark = pytest.mark.foundation

r, mu, phi = Symbolic.r, Symbolic.mu, Symbolic.phi

_ROWS = [
    # (expression, depends on r, on mu, on phi)
    (sp.Integer(2), False, False, False),
    (r**2 + 1, True, False, False),
    (1 + mu, False, True, False),
    (sp.cos(phi), False, False, True),
    (sp.sin(phi) ** 2 + sp.cos(phi) ** 2, False, False, False),        # a free-symbols test says phi
    (mu - mu + r - r, False, False, False),
    (sp.Piecewise((1, phi < sp.pi), (0, True)), False, False, True),   # a derivative test says no
    (sp.Piecewise((1, r < 1), (2, True)), True, False, False),
    (r * mu * sp.sin(phi), True, True, True),
]
_IDS = ["2", "r^2+1", "1+mu", "cos phi", "sin^2+cos^2", "cancelled", "step in phi", "step in r", "r mu sin phi"]


@pytest.mark.rests_on("tests/gates/numerics/test_mesh_free_function.py::test_s6_15_the_anisotropy_predicate")
@pytest.mark.parametrize("expression,on_r,on_mu,on_phi", _ROWS, ids=_IDS)
def test_s8_11_the_azimuth_predicate(expression, on_r: bool, on_mu: bool, on_phi: bool) -> None:
    f = Symbolic.of(expression)
    got = (f.depends_on(r), f.depends_on(mu), f.depends_on(phi))
    require(got == (on_r, on_mu, on_phi), f"{expression}: depends on (r, mu, phi) = {got}")
    require(f.depends_on(r, mu, phi) is (on_r or on_mu or on_phi), f"{expression}: the joint predicate")
    require(f.is_isotropic is not (on_mu or on_phi), f"{expression}: is_isotropic disagrees with depends_on(mu, phi)")


def test_s8_11_one_dependent_group_makes_the_function_dependent() -> None:
    require(Symbolic.of(1, sp.cos(phi)).depends_on(phi), "a phi-dependent second group was missed")


@pytest.mark.parametrize(
    "coordinates",
    [pytest.param((), id="none"), pytest.param((sp.Symbol("r"),), id="r-without-real"),
     pytest.param((sp.Symbol("x", real=True),), id="a-stray-symbol"), pytest.param((r, sp.Symbol("t")), id="owned-and-stray")],
)
def test_s8_11_only_owned_coordinates_are_asked(coordinates) -> None:
    """An empty call would read "depends on nothing" (always False), and a symbol the
    function does not own is not one of its coordinates (``Symbol("r")`` without
    ``real=True`` is another symbol to SymPy): both are refused, keyed."""
    with pytest.raises(ValueError, match="name at least one owned coordinate"):
        Symbolic.of(r * mu).depends_on(*coordinates)
