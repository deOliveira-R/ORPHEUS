r"""One definition of "a real number" at L1 (#559).

The rule, NaN refused and ``-0.0`` made ``+0.0``, is
:func:`orpheus.values.scalars.canonical_real`, and the conversion every
admitting site makes, the rule plus the 2**53 and overflow refusals, is
:func:`~orpheus.values.scalars.exact_double` (``canonical_reals``
entrywise, which the encoder's array and sparse paths call). Before #559 it was spelled
three times: the content encoder's private ``_real``, the geometry layer's
``parse_real``, and ``RegionwiseConstant``'s inline parse. This file is the
witness that keeps the copies from coming back (X4):

- the NaN refusal's shortest distinctive fragment, ``which is not a number``,
  is written in exactly 1 of the production files under ``orpheus/`` (an AST
  pass over string constants, with the positive control that it finds the
  one in ``scalars.py``), and every admitting site raises it;
- the three sites agree on the canonical bits of ``-0.0``.

The rows here are values' own: the AST censuses over ``orpheus/`` (they read source, they import no
consumer) and the values sites of each refusal. The consumer sites of the same refusals
(``RegionwiseConstant`` in numerics, ``Mixture`` in data) and the signed-zero agreement across the three
sites are a census of this contract over its consumers, and live at the root,
``tests/gates/test_one_real_parser_census.py`` (layer DAG ruling 28, 2026-10-26).

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

from orpheus.values.content import encode
from orpheus.values.scalars import canonical_real, parse_finite_real, parse_finite_reals, parse_real

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


def test_no_other_raise_names_nan() -> None:
    """The second filter, in another vocabulary (X2 step 5): by AST, every
    ``raise`` under ``orpheus/`` whose message text says ``NaN`` sits in
    ``scalars.py``. First red (``[M]`` 2026-10-02, the archivist's finding):
    ``Mixture`` kept its own loop, "holds NaN, which is not a cross section",
    which the fragment census above could not see."""
    sites = []
    for path in sorted((_ROOT / "orpheus").rglob("*.py")):
        for node in ast.walk(ast.parse(path.read_text(), filename=str(path))):
            if isinstance(node, ast.Raise) and node.exc is not None:
                text = " ".join(
                    c.value for c in ast.walk(node.exc) if isinstance(c, ast.Constant) and isinstance(c.value, str)
                )
                if "NaN" in text:
                    sites.append(f"{path.relative_to(_ROOT)}:{node.lineno}")
    if not any(site.startswith("orpheus/values/scalars.py:") for site in sites):
        pytest.fail(f"activation: the AST pass did not find scalars.py's raise; it found {sites}")
    outside = [site for site in sites if not site.startswith("orpheus/values/scalars.py:")]
    if outside:
        pytest.fail(f"a NaN refusal outside the one parser: {outside}")


def test_the_nan_refusal_is_written_once() -> None:
    sites = _sites_writing_the_fragment()
    if not any(site.startswith("orpheus/values/scalars.py:") for site in sites):
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
    ],
    ids=["canonical_real", "parse_real", "parse_finite_real", "parse_finite_reals", "encode"],
)
def test_every_site_refuses_nan_with_the_one_message(admit) -> None:
    with pytest.raises(ValueError, match=_FRAGMENT):
        admit()


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
    ],
    ids=["parse_real", "parse_finite_reals_list"],
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
    """One conversion, :func:`~orpheus.values.scalars.exact_double`: an
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
