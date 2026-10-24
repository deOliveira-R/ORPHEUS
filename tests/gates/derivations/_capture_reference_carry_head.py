r"""Capture the reference carry fixture from ``main`` ``5d04c296`` (wave 2, S5; plan 6.5 gate 2).

NOT a test and NOT runnable on the carved tree: it spells each wall of ``_reference_carry_cases.CASES`` as the retired
string tag ``BC(...)`` (or the typed ``WhiteBoundary``), which only the HEAD tree has. It was run ONCE, from a detached
extraction of ``5d04c296`` (``git archive 5d04c296 orpheus``), by a script outside the repository, and wrote
``reference_carry_head.json`` next to it::

    cd <outside the repo> && PYTHONPATH=<extraction> .venv/bin/python -O <this file> <output.json>

and printed ``orpheus.__file__`` first, so the capture is known to have read HEAD and not the carved tree. The fixture holds,
per case: the walls ``Walls.of(geometry)`` reads, as ``(breakpoint, specular, diffuse, partner)``; the characteristic
reference's fundamental eigenvalue at a fixed resolution, read in process (``traced_memo.bypass()``) so a stale memo entry
cannot stand in for a solve, as the hex of the float (``float.hex``), so the comparison is bit-for-bit; and, for the
specular-only bodies, ``specular_albedos(geometry)``. A solve that the HEAD reference refuses is recorded as
``{"refused": "<exception type>: <first line>"}``, so a refusal that moves is also a finding.
"""

# pyright: reportAttributeAccessIssue=false, reportArgumentType=false, reportCallIssue=false
# (the module's subject is the retired API of ``main`` ``5d04c296``: against the carved tree every import below is unresolved by design)

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _cases():
    spec = importlib.util.spec_from_file_location("_reference_carry_cases", HERE / "_reference_carry_cases.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["_reference_carry_cases"] = module
    spec.loader.exec_module(module)
    return module


def main(out: str) -> None:
    import orpheus
    from orpheus.data.cells import CellCoefficient, Channel
    from orpheus.data.materials import Materials
    from orpheus.derivations.common.reference_body import specular_albedos
    from orpheus.derivations.continuous.characteristic import (
        CharacteristicDerivation, Resolution, TransportResolution, Walls,
    )
    from orpheus.derivations.continuous.sood_registry import sood2003 as sood
    from orpheus.geometry import BC, StructuredGeometry
    from orpheus.geometry.boundary import WhiteBoundary
    from orpheus.numerics.observable import Eigenvalue
    from orpheus.numerics.question import Eigen
    from orpheus.numerics.traced_memo import bypass
    from orpheus.specification.specification import GeometrySpecification

    print("orpheus read from", orpheus.__file__, flush=True)
    cases = _cases()
    mixtures = {0: sood.PU_2_0_SL_STUB.materials[0], 1: sood.URRB_2_0_IN.materials[0]}

    def law(decl, sign):
        match decl:
            case ("vacuum",):
                return BC.vacuum
            case ("mirror",):
                return BC.reflective
            case ("partial", a):
                return BC("partial", {"albedo": a})
            case ("white", a):
                return WhiteBoundary(axis="x", outward_sign=sign, albedo=a)
            case ("periodic",):
                return BC("periodic")
        raise ValueError(decl)

    def geometry(chart, bp, faces):
        mats = cases.mat_ids(bp)
        if chart == "slab":
            return StructuredGeometry.slab(bp, mats, left=law(faces[0], -1), right=law(faces[1], +1))
        build = StructuredGeometry.sphere if chart == "sphere" else StructuredGeometry.cylinder
        if len(faces) == 1:
            return build(bp, mats, outer=law(faces[0], +1))
        return build(bp, mats, inner=law(faces[0], -1), outer=law(faces[1], +1))

    resolution = Resolution(1, 0, 0.4, TransportResolution(8), 8)
    result: dict[str, dict] = {"_resolution": {"degree": 1, "layers": 0, "ratio": 0.4, "line_points": 8, "projection_points": 8}}
    for name, chart, bp, faces in cases.CASES:
        g = geometry(chart, bp, faces)
        walls = Walls.of(g)
        row: dict = {"walls": [[w.breakpoint, w.specular, w.diffuse, w.partner] for w in walls.walls]}
        if name in cases.SPECULAR_ONLY:
            row["albedos"] = list(specular_albedos(g, owner="carry"))
        t0 = time.time()
        try:
            spec = GeometrySpecification(Materials(mixtures), g, Eigen(CellCoefficient.every(Channel.FISSION_EMISSION)))
            with bypass():
                k = CharacteristicDerivation(spec, resolution).evaluate(Eigenvalue()).value
            row["k_hex"] = float(k).hex()
        except Exception as exc:  # noqa: BLE001  a refusal is data
            row["k"] = {"refused": f"{type(exc).__name__}: {str(exc).splitlines()[0][:160]}"}
        row["seconds"] = round(time.time() - t0, 1)
        print(name, row, flush=True)
        result[name] = row
    Path(out).write_text(json.dumps(result, indent=1, sort_keys=True) + "\n")


if __name__ == "__main__":
    main(sys.argv[1])
