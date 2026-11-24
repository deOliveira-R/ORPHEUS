r"""G02h: a trace space whose half-ranges disagree with its own ``omega_dot_n`` is refused at construction (N02 review findings).

``AngularTraceSpace.half_ranges`` repeats a fact that ``omega_dot_n`` carries (the sign of each row), is excluded from equality and from
the name digest, so a ``dataclasses.replace`` could pair them wrongly and the result would compare equal to the honest space. The check
sits in ``__post_init__`` under ``if __debug__:`` (a misuse check of the raw constructor, not a guarantee), so the refusal is read in a
subprocess WITHOUT ``-O``; the canonical ``python -O`` strips the block, and the same code run with ``-O`` accepts all three bad spaces.

Four constructions on ``Quadrature.product(4, 4)`` with faces ``xmin`` and ``xmax``: the honest ``replace`` (accepted), the half-ranges
reversed (the outflow of ``xmax`` would be the inflow), one half-range dropped (count differs from the rows), and a half-range on the
wrong number of ordinates. The three bad ones are refused with a message naming the disagreement.

Exact on any host: the table is a sign times an exact cosine and the check compares integer signs; no LAPACK or libm result is read.

First red at the base of this commit: the three bad spaces are accepted (the post-init does not exist), read by ``test_the_bad_spaces_are_refused``.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.foundation

_ROOT = Path(__file__).resolve().parents[4]

CODE = r"""
import dataclasses
from orpheus.numerics.face_layout import FaceLayout
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces import AngularTraceSpace
from orpheus.numerics.tangential import HalfRange
from orpheus.numerics.declaration import Coord

quad = Quadrature.product(4, 4)
n = int(quad.n_ordinates)
layout = FaceLayout.from_named_shapes([("xmin", (n, 2)), ("xmax", (n, 2))])
t = AngularTraceSpace.from_quadrature_and_layout(quad, layout)
short = HalfRange(Coord(3, (0,)), Coord(3, (1,)), Coord(3, (2,)), 3)
cases = (
    ("honest", dict(half_ranges=tuple(t.half_ranges))),
    ("reversed", dict(half_ranges=tuple(reversed(t.half_ranges)))),
    ("dropped", dict(half_ranges=t.half_ranges[:1])),
    ("wrong_n", dict(half_ranges=(short, short))),
)
for name, kw in cases:
    try:
        s = dataclasses.replace(t, **kw)
    except ValueError as e:
        print(name + " refused: " + str(e))
    else:
        print(name + " accepted equal=" + str(s == t))
"""


def _run(*flags: str) -> list[str]:
    env = {k: v for k, v in os.environ.items() if k != "PYTHONOPTIMIZE"}
    env["PYTHONPATH"] = os.pathsep.join([str(_ROOT), *([env["PYTHONPATH"]] if env.get("PYTHONPATH") else [])])
    r = subprocess.run([sys.executable, *flags, "-c", CODE], env=env, cwd=_ROOT, capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        pytest.fail(f"the subprocess ({flags}) exited {r.returncode}: {r.stderr[-400:]}")
    return r.stdout.strip().splitlines()


def refusal_violations(lines: list[str]) -> list[str]:
    r"""What is wrong with the printed outcome: the honest space must be accepted and each bad one refused with its reason."""
    want = {"honest": None, "reversed": "disagrees with the sign of omega_dot_n", "dropped": "half-ranges for 2 faces", "wrong_n": "ordinates"}
    got = {line.split(" ", 1)[0]: line.split(" ", 1)[1] for line in lines}
    bad: list[str] = []
    for name, reason in want.items():
        text = got.get(name, "<missing>")
        if reason is None and text != "accepted equal=True":
            bad.append(f"{name}: {text}")
        if reason is not None and not (text.startswith("refused: ") and reason in text):
            bad.append(f"{name}: {text}")
    return bad


def test_the_bad_spaces_are_refused() -> None:
    assert refusal_violations(_run()) == []  # no -O: the check stands


def test_the_refusal_reading_reds_on_the_optimized_run() -> None:
    r"""Control (X1): the same reading applied to the ``-O`` run, where the guard is stripped, reports all three bad spaces accepted."""
    lines = _run("-O")
    assert len(refusal_violations(lines)) == 3, lines
