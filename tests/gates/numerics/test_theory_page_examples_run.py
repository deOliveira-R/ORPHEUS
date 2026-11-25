r"""The code example of the axis-marginal section of ``docs/theory/foundations/spaces.rst`` runs, and its assertions hold.

The example is the page's published use of ``TensorSpace.of_axes`` with the three marginal arrows. A page that states a constructor that does not exist
(the example once called ``FunctionSpace.of_axes(Axis(...))``, retired and refusing an ``Axis``) fails here: the test reads the text of the page, compiles it
with assertions kept (``optimize=0``, so ``python -O`` does not strip them) and executes it.

First red: the page's example before the repair raised ``AttributeError: type object 'FunctionSpace' has no attribute 'of_axes'``. Mutation: replacing
``TensorSpace.of_axes`` by ``TensorSpace.of_axis`` in the page text (planted on the in-memory copy) reddens it. Exact on any host: the example's own
assertions are ``nulp`` and identity checks.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

pytestmark = pytest.mark.foundation

PAGE = Path(__file__).resolve().parents[3] / "docs" / "theory" / "foundations" / "spaces.rst"
ANCHOR = "axis label.\n\n.. code-block:: python\n\n"


def _example(text: str) -> str:
    start = text.index(ANCHOR) + len(ANCHOR)
    lines: list[str] = []
    for line in text[start:].splitlines():
        if line.strip() and not line.startswith("   "):
            break
        lines.append(line)
    return textwrap.dedent("\n".join(lines))


def _run(source: str) -> None:
    exec(compile(source, "spaces.rst:axis-marginal-example", "exec", optimize=0), {"__name__": "docs_example"})


def test_the_axis_marginal_example_runs() -> None:
    _run(_example(PAGE.read_text()))


def test_a_stale_constructor_in_the_example_is_caught() -> None:
    source = _example(PAGE.read_text())
    assert "TensorSpace.of_axes(" in source, "the example no longer calls the constructor this gate plants against"
    with pytest.raises(AttributeError):
        _run(source.replace("TensorSpace.of_axes(", "TensorSpace.of_axis("))
