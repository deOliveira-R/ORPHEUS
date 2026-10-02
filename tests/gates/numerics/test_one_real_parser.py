r"""One definition of "a real number" at L1 (#559).

The rule, NaN refused and ``-0.0`` made ``+0.0``, is
:func:`orpheus.numerics.scalars.canonical_real`, and the conversion every
admitting site makes, the rule plus the 2**53 and overflow refusals, is
:func:`~orpheus.numerics.scalars.exact_double` (``canonical_reals``
entrywise, which the encoder's array and sparse paths call). Before #559 it was spelled
three times: the content encoder's private ``_real``, the geometry layer's
``parse_real``, and ``RegionwiseConstant``'s inline parse. This file is the
witness that keeps the copies from coming back (X4):

- the NaN refusal's shortest distinctive fragment, ``which is not a number``,
  is written in exactly 1 of the production files under ``orpheus/`` (an AST
  pass over string constants, with the positive control that it finds the
  one in ``scalars.py``), and every admitting site raises it;
- the three sites agree on the canonical bits of ``-0.0``.

First red (``[M]`` 2026-10-02, the tree before #559): the fragment was
written at 2 sites (``geometry/scalars.py`` and ``mesh_free_function.py``),
and the encoder refused NaN with a message of its own.
"""

from __future__ import annotations

import ast
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

from orpheus.numerics.content import encode
from orpheus.numerics.mesh_free_function import RegionwiseConstant
from orpheus.numerics.scalars import canonical_real, parse_finite_real, parse_finite_reals, parse_real

pytestmark = pytest.mark.foundation

_ROOT = Path(__file__).resolve().parents[3]
_FRAGMENT = "which is not a number"


def _sites_writing_the_fragment() -> list[str]:
    sites = []
    for path in sorted((_ROOT / "orpheus").rglob("*.py")):
        for node in ast.walk(ast.parse(path.read_text(), filename=str(path))):
            if isinstance(node, ast.Constant) and isinstance(node.value, str) and _FRAGMENT in node.value:
                sites.append(f"{path.relative_to(_ROOT)}:{node.lineno}")
    return sites


def test_the_nan_refusal_is_written_once() -> None:
    sites = _sites_writing_the_fragment()
    if not any(site.startswith("orpheus/numerics/scalars.py:") for site in sites):
        pytest.fail(f"activation: the AST pass did not find the fragment in scalars.py; it found {sites}")
    if len(sites) != 1:
        pytest.fail(f"the NaN refusal is spelled at {len(sites)} sites, not 1: {sites}")


@pytest.mark.parametrize(
    "admit",
    [
        lambda: canonical_real(math.nan, "x"),
        lambda: parse_real(math.nan, "x"),
        lambda: parse_finite_real(math.nan, "x"),
        lambda: parse_finite_reals(np.array([1.0, math.nan]), "x"),
        lambda: encode(math.nan),
        lambda: RegionwiseConstant(np.array([[math.nan]])),
    ],
    ids=["canonical_real", "parse_real", "parse_finite_real", "parse_finite_reals", "encode", "RegionwiseConstant"],
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


def test_a_finite_parse_refuses_infinity() -> None:
    with pytest.raises(ValueError, match="is infinite"):
        parse_finite_real(math.inf, "x")
    with pytest.raises(ValueError, match=r"the entry \(1,\) is infinite"):
        parse_finite_reals(np.array([0.0, -math.inf]), "x")
    with pytest.raises(TypeError, match="must be real numbers"):
        parse_finite_reals(np.array([True]), "x")


@pytest.mark.parametrize(
    "admit",
    [
        lambda: parse_real(True, "x"),
        lambda: parse_finite_reals([1.0, True], "x"),
        lambda: RegionwiseConstant([[1.0, True]]),  # type: ignore[arg-type]
    ],
    ids=["parse_real", "parse_finite_reals_list", "RegionwiseConstant_list"],
)
def test_a_bool_entry_is_refused_by_every_parse(admit) -> None:
    """One definition of a real ENTRY: a ``bool`` hidden in a list is refused
    as the scalar parse refuses it. First red (``[M]`` 2026-10-02, the step-7
    elegance review): ``np.asarray([1.0, True])`` cast the bool to 1.0, so the
    list spelling admitted what the scalar spelling refused."""
    with pytest.raises(TypeError, match="must be a real number, got bool"):
        admit()


@pytest.mark.parametrize(
    "admit,error,fragment",
    [
        (lambda: parse_real(2**53 + 1, "x"), ValueError, "lies beyond 2\\*\\*53"),
        (lambda: encode(2**53 + 1), ValueError, "lies beyond 2\\*\\*53"),
        (lambda: encode(np.array([2**53 + 1])), ValueError, r"the entry \(0,\): .*lies beyond 2\*\*53"),
        (lambda: parse_real(10**400, "x"), ValueError, "lies beyond 2\\*\\*53"),
        (lambda: parse_real(Fraction(10**400, 3), "x"), ValueError, "beyond the range of a double"),
    ],
    ids=["parse_real_2**53+1", "encode_2**53+1", "encode_array_2**53+1", "parse_real_10**400", "parse_real_overflow"],
)
def test_the_parsers_and_the_encoder_agree_on_what_a_double_carries(admit, error, fragment: str) -> None:
    """One conversion, :func:`~orpheus.numerics.scalars.exact_double`: an
    integer a double would round is refused by the parsers as by the encoder,
    and an overflow is a keyed refusal. First red (``[M]`` 2026-10-02, qa of
    step 7): ``parse_real(2**53 + 1)`` stored ``2**53`` silently while the
    encoder refused it, and ``10**400`` raised an unkeyed ``OverflowError``."""
    with pytest.raises(error, match=fragment):
        admit()


def test_a_rank_zero_array_parses() -> None:
    parsed = parse_finite_reals(np.array(-0.0), "x")
    if parsed.shape != () or math.copysign(1.0, float(parsed)) != 1.0:
        pytest.fail(f"rank 0: got {parsed!r}")
