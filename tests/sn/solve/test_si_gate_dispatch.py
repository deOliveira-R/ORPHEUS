r"""C5.4 (#225) — the SI windowing / G-S gates key on GENUINE dimensionality.

vv-principles Mode 9: the pre-C5.4 gates keyed on ``sn_mesh.reduced is
None`` — a coincidence proxy that is true for EVERY multi-D Cartesian
mesh, not just 2-D. At d=3 that proxy would have silently
moment-windowed the SI iterate (the in-sweep moment emission is a 2-D
kernel; ``FullFieldWavefront`` refuses moment mode) — a corrupted-
iterate failure with no loud signal on the G-S path. These pins hold
the gate functions to the genuine conditions:

* ``_maybe_window`` → wrap ONLY at ``is_cartesian and ndim == 2``;
* ``_select_si_resolvent`` G-S → ``is_cartesian and not is_1d``
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

from orpheus.sn.solver import _maybe_window, _MomentWindowedResolvent
from orpheus.numerics.quadrature import Quadrature

pytestmark = [pytest.mark.foundation]


class _BaseResolvent:
    """Minimal stand-in — _maybe_window only wraps or passes it through."""


def _fake_mesh(*, is_cartesian: bool, ndim: int) -> SimpleNamespace:
    return SimpleNamespace(is_cartesian=is_cartesian, ndim=ndim)


def _scattering_stub():
    # _maybe_window reads ONLY scattering_op.frame (the angular Frame that
    # _MomentWindowedResolvent injects); stub exactly that.
    return SimpleNamespace(
        frame=Quadrature.level_symmetric(sn_order=4).angular_frame(0),
    )


def test_2d_cartesian_windows() -> None:
    base = _BaseResolvent()
    wrapped, windowed = _maybe_window(
        base, _scattering_stub(), _fake_mesh(is_cartesian=True, ndim=2),
    )
    np.testing.assert_equal(windowed, True)
    if not isinstance(wrapped, _MomentWindowedResolvent):
        pytest.fail("2-D Cartesian must moment-window the SI iterate")


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
    base = _BaseResolvent()
    wrapped, windowed = _maybe_window(
        base, _scattering_stub(), _fake_mesh(is_cartesian=is_cartesian, ndim=ndim),
    )
    np.testing.assert_equal(windowed, False)
    if wrapped is not base:
        pytest.fail(
            f"{label}: _maybe_window must pass the resolvent through "
            f"unwrapped (got {type(wrapped).__name__})"
        )
