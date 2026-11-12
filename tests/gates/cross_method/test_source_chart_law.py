r"""The chart law of #707 Step A: every data-side value is the spatial measure times its pre-carve density value.

The generic gate of the re-baseline protocol (spec ``scratch/true_ends/stepA/gate_spec.md`` §3). It replaces "re-read
each frozen source number" by ONE law checked on fresh anchors: for every flux → source atom, every input door and
the detector door,

.. math::

   y_{\rm now} \;=\; m \odot y_{\rm anchor}, \qquad
   m_{\rm SN\,bulk}[n, g, i, o] = V_i\,\theta^{|o|}, \quad m_{\rm trace} = 1 \ \text{(Step A)}, \quad
   m_{\rm diffusion\,bulk}[g, i] = V_i ,

with ``y_anchor`` frozen on the pre-carve tree (``source_chart_anchor.npz``, ``meta/commit``) and ``m`` built in
``_source_chart_cases`` from raw mesh volumes and a typed θ = 1/3 (no production metric object, no ``w_n``: R8).

Kind: THEOREM about the carve (the chart map is ♭ restricted to the spatial factor, a diagonal). RED today on every
non-vacuous row, by ``max|1 − m|`` (the store is the density chart); GREEN after Step A iff the integrated chart is
exactly R8's (a stray ``w_n``, a missing θ on a slope row, or a second ×V each red it by O(1)).

Tolerance, per block, ``|y_now − m⊙y_anchor| ≤ TOL_ULPS · eps · max|m⊙y_anchor|``, TOL_ULPS = 256. Derivation `[R]`:
after the carve each output element is computed integrated directly instead of as ``V·(density)``, which re-associates
one reduction per element; the deepest reduction in the population is the anisotropic scattering on 2-D LS4 P1,
``N · ng · (L+1)² = 12 · 2 · 4 = 96`` terms, plus the one product by ``m``; 256 is that depth with a ×2.5 margin. It
is normalised by the block maximum because a cancelling sum (the streaming-collision balance) has no meaningful
per-element relative error. The smallest chart error the gate must see is a missing θ (relative 2/3) or a stray
``w_n`` (≥ 0.17 on GL8): ≥ 1e12 above the band.

A row whose anchor is identically zero (the library's ``Sig2 = 0`` makes every ``N2NOperator`` block zero; the
scattering and fission trace blocks are zero by construction) is VACUOUS: it is asserted to stay zero and is not
counted as coverage (``test_population``).
"""

from __future__ import annotations

import warnings
from functools import cache

import numpy as np
import pytest

from tests.gates.cross_method import _source_chart_cases as cc

pytestmark = [pytest.mark.foundation, pytest.mark.filterwarnings("ignore::Warning")]

TOL_ULPS = 256
_EPS = np.finfo(float).eps


@cache
def _anchor() -> dict[str, np.ndarray]:
    with np.load(cc.ANCHOR) as data:
        return {k: data[k] for k in data.files}


@cache
def _problem(name: str):
    return cc.problems()[name]()


def _atom_keys() -> list[str]:
    a = _anchor()
    return sorted(k for k in a if k.split("/")[1] in cc.ATOMS)


def _nonvacuous(keys: list[str]) -> list[str]:
    a = _anchor()
    return [k for k in keys if float(np.max(np.abs(a[k]), initial=0.0)) > 0.0]


def _check(key: str, now: np.ndarray, measure: np.ndarray) -> None:
    want = measure * _anchor()[key]
    assert now.shape == want.shape, f"{key}: shape {now.shape} != {want.shape}"
    scale = float(np.max(np.abs(want), initial=0.0))
    err = float(np.max(np.abs(now - want), initial=0.0))
    if scale == 0.0:
        assert err == 0.0, f"{key}: a vacuous (zero) block became nonzero: {err:.3e}"
        return
    rel = err / scale
    today = float(np.max(np.abs(now - _anchor()[key]))) / max(float(np.max(np.abs(now))), 1e-300)
    assert rel <= TOL_ULPS * _EPS, (
        f"{key}: the stored chart is not m ⊙ density: max|now − m⊙anchor|/max|m⊙anchor| = {rel:.3e} "
        f"(band {TOL_ULPS * _EPS:.1e}); distance to the unscaled anchor {today:.3e} "
        f"(0 means the store is still the density chart)"
    )


def test_population() -> None:
    """X2: the anchor's counts, the vacuous rows named, the commit an ancestor of HEAD."""
    import subprocess

    keys = _atom_keys()
    live = _nonvacuous(keys)
    assert len(keys) == 44 and len(live) == 21, (len(keys), len(live))
    assert not any("N2NOperator" in k for k in live), "N2N became non-vacuous: re-read the VACUOUS note"
    commit = str(_anchor()["meta/commit"])
    ok = subprocess.run(["git", "merge-base", "--is-ancestor", commit, "HEAD"], capture_output=True).returncode == 0
    assert ok, f"the anchor commit {commit} is not an ancestor of HEAD"


@pytest.mark.parametrize("key", _atom_keys())
def test_atom_output_is_integrated(key: str) -> None:
    """Every flux → source atom's output: bulk ``V θ^{|o|} ⊙`` density, trace unchanged (Step A)."""
    name, atom, block = key.split("/")
    a = _anchor()
    problem = _problem(name)
    psi = cc.state_from(problem, a[f"{name}/state/bulk"], a[f"{name}/state/trace"])
    out = cc.atoms(problem)[atom].apply(psi)
    now = np.asarray(out.interior.values if block == "bulk" else out.boundary.values, dtype=float)
    measure = cc.sn_bulk_measure(problem, now.shape) if block == "bulk" else np.ones(now.shape)
    _check(key, now, measure)


def _door_keys() -> list[str]:
    return sorted(k for k in _anchor() if "/door/" in k and k.endswith("_stored"))


@pytest.mark.parametrize("key", _door_keys())
def test_input_door_integrates_once(key: str) -> None:
    """R10: a density handed to the external-source door is stored as ``V θ^{|o|} ⊙`` its lifted density, once."""
    from orpheus.sn.solver import _build_fixed_source_rhs

    name = key.split("/")[0]
    problem = _problem(name)
    given = _anchor()[key.replace("_stored", "_input")]
    now = np.asarray(cc._first_system(_build_fixed_source_rhs(given, problem)).interior.values, dtype=float)
    _check(key, now, cc.sn_bulk_measure(problem, now.shape))


def test_detector_door_integrates_once() -> None:
    """R10: the detector response Σ_d (a density) reaches the adjoint posing as ``r = V ⊙ Σ_d`` per ordinate, no w_n."""
    problem = _problem("slab_dd_p1")
    a = _anchor()
    sigma_d = a["slab_dd_p1/detector/input"]
    import orpheus.sn.solver as solver_mod

    seen = []
    original = solver_mod.SourcePosing.of_pencil

    def spy(pencil, sigma, source, *args, **kwargs):
        seen.append(source)
        return original(pencil, sigma, source, *args, **kwargs)

    solver_mod.SourcePosing.of_pencil = spy
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                solver_mod.solve_sn_adjoint_fixed_source(
                    cc._mixtures(), cc._slab_mesh(), problem.quad, sigma_d, max_inner=3, inner_tol=1e-3,
                    boundary_conditions=cc._conditions(cc._slab_mesh()))
            except Exception:
                pass
    finally:
        solver_mod.SourcePosing.of_pencil = original
    assert seen, "the detector door was not reached (SourcePosing.of_pencil never called): re-point the spy"
    now = np.asarray(cc._first_system(seen[0]).interior.values, dtype=float)
    _check("slab_dd_p1/detector/stored", now, cc.sn_bulk_measure(problem, now.shape))


@pytest.mark.parametrize("label", ["loss", "fission"])
@pytest.mark.parametrize("block", ["bulk", "trace"])
def test_diffusion_output_is_integrated(label: str, block: str) -> None:
    """Diffusion's flux → source operators: bulk ``V ⊙`` density (census P18, ``operators.py:418``), trace unchanged."""
    solver = cc.diffusion_problem()
    phi = solver.template
    member = type(phi).from_flat(_anchor()["diffusion/state"], phi)
    out = getattr(solver, label).apply(member)
    now = np.asarray(out.interior.values if block == "bulk" else out.boundary.values, dtype=float)
    if block == "bulk":
        V = cc.diffusion_volumes(solver)
        measure = np.broadcast_to(V.reshape((1,) + V.shape), now.shape)
    else:
        measure = np.ones(now.shape)
    _check(f"diffusion/{label}/{block}", now, measure)
