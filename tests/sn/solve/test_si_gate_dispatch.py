r"""C5.4 (#225) — the SI windowing / G-S gates key on GENUINE dimensionality.

vv-principles Mode 9: the pre-C5.4 gates keyed on ``sn_mesh.reduced is
None`` — a coincidence proxy that is true for EVERY multi-D Cartesian
mesh, not just 2-D. At d=3 that proxy would have silently
moment-windowed the SI iterate (the in-sweep moment emission is a 2-D
kernel; ``FullFieldWavefront`` refuses moment mode) — a corrupted-
iterate failure with no loud signal on the G-S path. These pins hold
the gate functions to the genuine conditions:

* ``_maybe_window`` → wrap ONLY at ``is_cartesian and ndim == 2``;
* ``_select_si_splitting`` G-S → ``is_cartesian and not is_1d``
  (multi-D; the schedule + scheduled sweep are d-generic since C3 —
  d=3 G-S FP-invariance is value-gated by the C5.5 Mode-9 box).

The gate functions read only ``is_cartesian`` / ``ndim`` / ``is_1d``
(plus, for a constructed G-S resolvent, the schedule inputs), so the
rows here pin the gate FUNCTIONS cheaply with a duck-typed mesh; the
constructible real-mesh twin lives in
``tests/sn/solve/test_d3_admission.py::test_d3_real_mesh_window_passthrough_and_gs_admissible``
(landed with C5.5) and the value-level FP-invariance in the Mode-9 box
there.

Assertions are ``np.testing`` / ``pytest.fail`` only (Mode-8 safe
under the canonical ``-O``).
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from orpheus.numerics.space import FunctionSpace
from orpheus.numerics.spaces.full_field_space import FullFieldSpace
from orpheus.sn.operators.sweep_operator import SweepOperator
from orpheus.sn.operators.windowing import WindowedSweep
from orpheus.sn.solver import _maybe_window
from orpheus.numerics.quadrature import Quadrature

pytestmark = [pytest.mark.foundation]


class _BaseResolvent:
    """Minimal forward stand-in whose ``.inverse()`` is what the caller
    hands to ``_maybe_window``.  #226 step 3: ``_maybe_window`` receives
    the INVERSE (a real :class:`SweepOperator`) and either composes the
    typed product ``P @ A.inverse()`` (2-D Cartesian) or passes it
    through — so the stub carries an ``inverse()`` returning a real
    :class:`SweepOperator` and the shared ``sn_mesh`` handle (the
    ``WindowedSweep`` mesh-identity guard compares it to ``P``'s)."""

    def __init__(self, sn_mesh) -> None:
        self.sn_mesh = sn_mesh

    def inverse(self) -> SweepOperator:
        return SweepOperator(self)


def _fake_mesh(*, is_cartesian: bool, ndim: int) -> SimpleNamespace:
    # F-1: BulkAnalysisOperator POSES its codomain at construction
    # (``from_blocks(face moment bulk, trace)`` — the pre-F-1 ``None``
    # codomain debt died), so the fake carries a minimal posed composite
    # rather than ``None``. Content is irrelevant to the gate — only the
    # trace block's existence is consumed.
    return SimpleNamespace(
        is_cartesian=is_cartesian, ndim=ndim,
        full_field_space=FullFieldSpace.from_blocks(
            FunctionSpace("fake_bulk", (2,)), FunctionSpace("fake_trace", (2,)),
        ),
    )


def _scattering_stub():
    # F-1: _maybe_window reads ONLY scattering_op.flux_analysis (the MINTED
    # face the windowed product's P factor is built from); the surrogate
    # honours exactly the surface the factory + BulkAnalysisOperator
    # constructor consume — ``.codomain`` (the posed moment bulk) and
    # ``.frame`` (what the fused apply would read; unused by the gate).
    # The real-face twin lives in test_2d_anisotropic_windowing (real S,
    # real mesh) per this module's cheap-gate-pin philosophy.
    return SimpleNamespace(
        flux_analysis=SimpleNamespace(
            codomain=FunctionSpace("fake_moment_bulk", (3,)),
            frame=Quadrature.level_symmetric(sn_order=4).angular_frame(0),
        ),
    )


def test_2d_cartesian_windows() -> None:
    mesh = _fake_mesh(is_cartesian=True, ndim=2)
    sweep = _BaseResolvent(mesh).inverse()
    step, windowed = _maybe_window(sweep, _scattering_stub(), mesh)
    np.testing.assert_equal(windowed, True)
    if not isinstance(step, WindowedSweep):
        pytest.fail(
            "2-D Cartesian must return the FUSED typed product "
            f"P @ A.inverse() directly (#226 step 3 — no adapter); "
            f"got {type(step).__name__}"
        )
    if step.sweep is not sweep:
        pytest.fail(
            "the product's inverse factor must BE the sweep the caller "
            "built (the composition, not a re-derivation)"
        )


@pytest.mark.parametrize("is_cartesian,ndim,label", [
    (True, 3, "d3-cartesian"),
    (True, 1, "slab"),
    (False, 1, "curvilinear"),
])
def test_non_2d_passthrough(is_cartesian: bool, ndim: int, label: str) -> None:
    """The Mode-9 guard: ONLY genuine 2-D Cartesian windows.

    The d=3 row is the load-bearing one — under the pre-C5.4
    ``reduced is None`` proxy it would have windowed.
    """
    mesh = _fake_mesh(is_cartesian=is_cartesian, ndim=ndim)
    sweep = _BaseResolvent(mesh).inverse()
    step, windowed = _maybe_window(sweep, _scattering_stub(), mesh)
    np.testing.assert_equal(windowed, False)
    if step is not sweep:
        pytest.fail(
            f"{label}: _maybe_window must pass the inverse operator "
            f"through unwrapped (got {type(step).__name__})"
        )
