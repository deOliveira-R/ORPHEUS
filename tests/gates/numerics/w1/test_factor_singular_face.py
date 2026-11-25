r"""G03b: the tangential face slot is quotiented by its kernel, and the round trip is the identity on ``S minus K`` only (N03, spec 7.4).

``Axis.of_quadrature(Quadrature.product(2, 4), "angular", face_axis=0)`` has the form ``Diag(abs(cos) w)``, whose kernel is the tangential
ordinates. A factor with ``sk`` whole is refused with the message naming the kernel to quotient by. On ``sk = V/tang`` the ``riesz`` image
returns the factor, and with a ``psi`` NONZERO on the tangential rows ``f.riesz().form.apply(f.form.apply(psi))`` equals ``psi`` on ``S minus K``
to ``1e-12`` and differs from ``psi`` on the tangential rows (by more than ``1e-2``: the draw ``default_rng(0)`` reads 0.947).

Exact on any host: a ``1e-12`` comparison of a computation with a computation and a lower bound of ``1e-2`` on a miss. First red at
``ff305e0f``: import-only (``orpheus.numerics.factor``); ``Axis.of_quadrature(face_axis=)`` is a ``TypeError`` there too. Mutations: compare
on ALL coordinates (reddens, the miss is 0.947), or allow a Riesz leg on the whole slot (the flat-then-sharp round trip misses ``psi`` by
more than ``1e-2``).
"""

from __future__ import annotations

import numpy as np
import pytest

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def test_singular_face_slot() -> None:
    from orpheus.numerics import factor as F
    from orpheus.numerics.axis import Axis
    from orpheus.numerics.declaration import Coord
    from orpheus.numerics.quadrature.directional import Quadrature

    ax = Axis.of_quadrature(Quadrature.product(2, 4), "angular", face_axis=0)
    n = ax.shape[0]
    w = np.asarray(ax.form.weights)
    tang = tuple(int(i) for i in np.flatnonzero(w == 0.0))
    live = [i for i in range(n) if i not in tang]
    _require(len(tang) == 4, f"the tangential set has {len(tang)} ordinates, expected 4 of product(2, 4)")
    v = F.Variance(F.Kind.MEASURE, F.VarianceChart.DENSITY)
    try:
        F.Factor(ax, v, F.Subquotient.whole(n))
    except ValueError as e:
        _require("quotient" in str(e), f"the whole-sk refusal does not name the kernel to quotient by: {e}")
    else:
        pytest.fail("a factor with a whole sk on the singular face slot was constructed")
    f = F.Factor(ax, v, F.Subquotient(Coord(n, tuple(range(n))), Coord(n, tang)))
    _require(f.riesz().riesz() == f, "riesz.riesz != id on V/tang")
    psi = np.random.default_rng(0).standard_normal(n)
    back = np.asarray(f.riesz().form.apply(f.form.apply(psi)))
    on_live = float(np.max(np.abs(back[live] - psi[live])))
    on_tang = float(np.max(np.abs(back[list(tang)] - psi[list(tang)])))
    _require(on_live <= 1e-12, f"the round trip misses psi on S minus K by {on_live:.3e}")
    _require(on_tang > 1e-2, f"the round trip equals psi on the tangential rows (miss {on_tang:.3e}): the quotient is not read, or psi is zero there")
    _require(float(np.max(np.abs(back - psi))) > 1e-2, "the round trip on ALL coordinates is the identity: the slot is not singular")
