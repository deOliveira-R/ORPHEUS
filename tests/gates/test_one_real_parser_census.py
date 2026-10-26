r"""The one real parser's consumers refuse what it refuses, with its message (#559; layer DAG ruling 28).

A census of a values contract over its consumers, so it lives at the root of ``tests/gates/``, beside
``test_layer_imports.py``: every site that admits a real number calls
:func:`orpheus.values.scalars.canonical_real` / :func:`~orpheus.values.scalars.exact_double`, so a NaN is refused
with the one message and ``-0.0`` is stored as ``+0.0`` everywhere. Moved on 2026-10-26 from
``tests/gates/values/test_one_real_parser.py``, which keeps the values sites of each refusal and the two AST
censuses showing the message is written once; the cases here are the consumer sites, ``RegionwiseConstant``
(numerics) and ``Mixture`` (data).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from orpheus.numerics.mesh_free_function import RegionwiseConstant
from orpheus.values.content import encode
from orpheus.values.scalars import parse_real
from tests.gates.values.test_one_real_parser import _FRAGMENT

pytestmark = pytest.mark.foundation


def _mixture_with_nan():
    from dataclasses import replace

    from orpheus.data.macro_xs.synthetic import get_mixture

    mixture = get_mixture("A", "2g")
    sig_t = np.array(mixture.SigT)
    sig_t[1] = math.nan
    return replace(mixture, SigT=sig_t)


@pytest.mark.parametrize(
    "admit",
    [
        lambda: RegionwiseConstant(np.array([[math.nan]])),
        lambda: _mixture_with_nan(),
    ],
    ids=["RegionwiseConstant", "Mixture"],
)
def test_every_site_refuses_nan_with_the_one_message(admit) -> None:
    with pytest.raises(ValueError, match=_FRAGMENT):
        admit()


def test_the_sites_agree_on_signed_zero() -> None:
    """``-0.0`` is stored as ``+0.0`` by the parsers and encodes as ``+0.0``."""
    if math.copysign(1.0, parse_real(-0.0, "x")) != 1.0:
        pytest.fail("parse_real keeps the sign of -0.0")
    table = RegionwiseConstant(np.array([[-0.0]])).values
    if math.copysign(1.0, float(table[0, 0])) != 1.0:
        pytest.fail("RegionwiseConstant keeps the sign of -0.0")
    if encode(-0.0) != encode(0.0):
        pytest.fail("the encoder separates -0.0 from +0.0")


@pytest.mark.parametrize(
    "admit",
    [
        lambda: RegionwiseConstant([[1.0, True]]),  # type: ignore[arg-type]
    ],
    ids=["RegionwiseConstant_list"],
)
def test_a_bool_entry_is_refused_by_every_parse(admit) -> None:
    """One definition of a real ENTRY: a ``bool`` hidden in a list is refused
    as the scalar parse refuses it. First red (``[M]`` 2026-10-02, the step-7
    elegance review): ``np.asarray([1.0, True])`` cast the bool to 1.0, so the
    list spelling admitted what the scalar spelling refused."""
    with pytest.raises(TypeError, match="must be a real number, got bool"):
        admit()
