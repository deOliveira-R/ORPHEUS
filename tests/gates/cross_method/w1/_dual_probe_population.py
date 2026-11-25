"""The population builder of the dual-probe record (tracked in the build as tests/gates/cross_method/w1/_dual_probe_population.py).

Not a test module. Every function RETURNS data and none asserts (``python -O`` strips an ``assert`` outside a collected module). It reaches the
spaces through the harness's own tables of named getters (``_snapshot.SN_ROWS``, ``DIFFUSION_ROWS``), keeps the axis-built ones (those that carry
``axes`` or, after N04, ``factors``: a space whose ``.axes`` is not None), and applies the three byte verbs to ``space.dual()``:
``apply_metric`` and ``apply_inverse_metric`` of the harness's seeded probe vector and ``metric_diagonal()``.
"""
from __future__ import annotations

import numpy as np

from tests.gates.cross_method.w1 import _fixtures as fx, _snapshot as snap

Key = tuple[str, str]


def _arrays_of_dual(space) -> dict[str, np.ndarray]:
    dual = space.dual()
    v = snap._probe_vector(int(dual.n_coordinates))
    out = {"apply_metric": snap._apply(dual, "apply_metric", v), "apply_inverse_metric": snap._apply(dual, "apply_inverse_metric", v)}
    diagonal = dual.metric_diagonal()
    if diagonal is not None:
        out["metric_diagonal"] = snap._flat(diagonal)
    return out


def capture_dual_probes() -> dict[Key, np.ndarray]:
    """``{(fixture, "d.dual.<space>.<verb>"): flat array}`` over every axis-built harness-reached space."""
    arrays: dict[Key, np.ndarray] = {}
    def take(fixture: str, name: str, space) -> None:
        if space is None or getattr(space, "axes", None) is None:
            return
        for verb, array in _arrays_of_dual(space).items():
            arrays[(fixture, f"d.dual.{name}.{verb}")] = array
    for kind in fx.SN_KINDS:
        for ng in fx.GROUPS:
            problem, _ = fx.sn_problem(kind, ng)
            for name, getter, _groups in snap.SN_ROWS:
                take(snap.sn_fixture_name(kind, ng), name, getter(problem))
    for ng in fx.GROUPS:
        mesh, _solver, _ = fx.diffusion_problem(ng)
        for name, getter, _groups in snap.DIFFUSION_ROWS:
            take(snap.diffusion_fixture_name(ng), name, getter(mesh))
    return arrays
