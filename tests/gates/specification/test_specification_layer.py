r"""The specification's layer contract (#405 P1 step 8, S8.3).

Specified by the test-architect (2026-10-02, ``.claude/plans/reference_p1_spec.md``
§1.8). ``orpheus/specification`` is an input-tier package above ``data``,
``geometry`` and ``numerics`` and below everything else: it may import those
three and nothing above, and nothing at or below its tier imports it. The
table rows live in ``tests/gates/test_layer_imports.py`` (the edit is in this
step's spec); these two rows are the run-time footprint and the package's own
import census, the plan's ban made concrete: the tempting spelling of the
group-count rule (S8.1 (c)) imports ``MaterialMesh``, which is transport.
"""

from __future__ import annotations

import ast
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from tests.gates._content_identity_helpers import require

pytestmark = pytest.mark.foundation

_ROOT = Path(__file__).resolve().parents[3]
_ABOVE = ("orpheus.transport", "orpheus.mesh", "orpheus.derivations", "orpheus.plotting") + tuple(
    f"orpheus.{p}" for p in ("sn", "pn", "moc", "cp", "mc", "diffusion", "kinetics", "fuel",
                             "thermal_hydraulics", "homogeneous")
)

_SCRIPT = textwrap.dedent(
    """
    import sys
    import numpy as np
    import orpheus
    from orpheus.data.cells import CellCoefficient, Channel
    from orpheus.data.materials import Materials
    from orpheus.numerics.mesh_free_function import Symbolic
    from orpheus.numerics.question import Eigen, FixedSource
    from orpheus.specification import Specification
    from tests.gates.specification._fixtures import fuel, moderator, slab2

    print("FILE", orpheus.__file__)
    Specification(materials=Materials({0: fuel(), 1: moderator()}), geometry=slab2(),
                  question=FixedSource(Symbolic.of(1 + Symbolic.mu, Symbolic.r), {CellCoefficient.every(Channel.SCATTERING_EMISSION): 0.1}))
    Specification(materials=Materials({0: fuel()}), geometry=None, question=Eigen(CellCoefficient.every(Channel.FISSION_EMISSION)))
    for name in sorted(sys.modules):
        if name.startswith("orpheus."):
            print("MODULE", name)
    """
)


def test_s8_3_building_a_specification_loads_nothing_above_the_input_tier() -> None:
    """In a fresh interpreter, import the package and build two specifications
    (one with a ``Symbolic`` source, so SymPy and every resolution path run):
    no module of ``transport``, ``mesh``, ``derivations`` or an L3 package is
    loaded. The fixtures module the script imports loads only data and geometry."""
    out = subprocess.run([sys.executable, "-O", "-c", _SCRIPT], cwd=_ROOT, capture_output=True, text=True,
                         timeout=300, env={**os.environ, "PYTHONPATH": str(_ROOT)})
    require(out.returncode == 0, out.stderr[-3000:])
    lines = out.stdout.splitlines()
    files = [line.split(" ", 1)[1] for line in lines if line.startswith("FILE ")]
    require(files and files[0].startswith(str(_ROOT)), f"the subprocess imported {files} (L22)")
    loaded = [line.split(" ", 1)[1] for line in lines if line.startswith("MODULE ")]
    require(any(m.startswith("orpheus.specification") for m in loaded), "activation: the package did not load")
    above = [m for m in loaded if m.startswith(_ABOVE)]
    require(not above, f"{len(above)} modules above the input tier loaded: {above[:10]}")


def test_s8_3_the_package_imports_only_data_geometry_numerics() -> None:
    """An AST census of every module of ``orpheus/specification`` (TYPE_CHECKING
    blocks included: the tier admits no annotation edge upward either)."""
    package = _ROOT / "orpheus" / "specification"
    modules = sorted(package.rglob("*.py"))
    require(modules, "activation: the package has no module")
    allowed = ("orpheus.data", "orpheus.geometry", "orpheus.numerics", "orpheus.specification")
    stray: list[str] = []
    for path in modules:
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            names = [node.module or ""] if isinstance(node, ast.ImportFrom) and node.level == 0 else (
                [a.name for a in node.names] if isinstance(node, ast.Import) else [])
            stray += [f"{path.name}:{getattr(node, 'lineno', 0)} {n}" for n in names if n.startswith("orpheus") and not n.startswith(allowed)]
    require(not stray, f"imports above the tier: {stray}")
