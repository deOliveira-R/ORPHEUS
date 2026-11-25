r"""The admission table of a dense form's induced subquotient: a read refuses what the declared kernel does not decide (N03, spec 3.1).

On a FULL-RANK dense axis a ``Spanned`` kernel ``K`` over the whole space is not the form's declared kernel (``Coord(n, ())``), so the density
form of ``whole/Spanned`` is refused with "not the form's declared kernel" and is never computed as the unquotiented axis form. On a DEFICIENT
dense axis (kernel ``Spanned(n, k, 1)``): ``whole/Spanned`` with the form's own kernel reads the axis form; ``whole/Spanned`` with another
key is refused at construction (the inclusion is undecided); a density factor on ``Ann(kernel)`` with ``K = Coord(n, ())`` reads the axis form;
on another ``Ann`` it is refused (at construction by the undecided meet, or at the read with "not the annihilator of the form's declared kernel"). Exact on any host (message fragments and identity of the returned form). First red at the commit before the
admission was tested: the first leg read the axis form silently. Mutation: drop the comparison of ``K`` with the form's kernel (reddens the first leg).
"""

from __future__ import annotations

import numpy as np
import pytest

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _refusal(read) -> str:
    try:
        read()
    except ValueError as e:
        return str(e)
    return ""


def test_the_declared_kernel_decides_the_dense_table() -> None:
    from orpheus.numerics import factor as F
    from orpheus.numerics.axis import Axis, BasisKind, Declared
    from orpheus.numerics.declaration import Ann, Coord, Spanned
    from orpheus.numerics.metric import DenseMetric, Face

    n = 4
    rng = np.random.default_rng(5)
    b = rng.standard_normal((n, n))
    full = Axis("h", (n,), origin=Declared(BasisKind.NODAL, DenseMetric(b @ b.T + np.eye(n), kernel=Coord(n, ()), face=Face.METRIC)))
    density = F.Variance(F.Kind.MEASURE, F.VarianceChart.DENSITY)
    whole = Coord(n, tuple(range(n)))
    f = F.Factor(full, density, F.Subquotient(whole, Spanned(n, ("x",), 1)))
    msg = _refusal(lambda: f.form)
    _require("not the form's declared kernel" in msg, f"whole/Spanned on a full-rank dense axis read without the refusal: {msg!r}")

    c = b[:, : n - 1]
    kernel = Spanned(n, ("declared", 1), 1)
    ax = Axis("h", (n,), origin=Declared(BasisKind.NODAL, DenseMetric(c @ c.T, kernel=kernel, face=Face.METRIC)))
    own = F.Factor(ax, density, F.Subquotient(whole, kernel))
    _require(own.form is ax.form, "whole/Spanned with the form's own kernel does not read the axis form")
    # another key: the inclusion of the form's kernel in K is undecided (a key is a generator's address, not a canonical form), so it is
    # refused at construction; the read is not reached
    msg = _refusal(lambda: F.Factor(ax, density, F.Subquotient(whole, Spanned(n, ("other", 2), 1))))
    _require("subspace inclusion is undecided" in msg, f"whole/Spanned with another key: {msg!r}")
    on_ann = F.Factor(ax, density, F.Subquotient(kernel.ann(), Coord(n, ())))
    _require(on_ann.form is ax.form, "a density factor on Ann(kernel) with K = Coord(n, ()) does not read the axis form")
    msg = _refusal(lambda: F.Factor(ax, density, F.Subquotient(Ann(n, Spanned(n, ("other", 2), 1)), Coord(n, ()))).form)
    _require(
        "not the annihilator of the form's declared kernel" in msg or "undecided" in msg,
        f"another Ann as S is refused at construction (undecided meet) or at the read: {msg!r}",
    )
