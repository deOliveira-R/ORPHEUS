r"""``chord_half_lengths`` against arbitrary precision, after its move onto the geometric kernel (#578, P2 group C).

:func:`~orpheus.derivations.common.kernels.chord_half_lengths` returns, for
each impact parameter :math:`y`, the half-chord a line spends in each
concentric region,
:math:`\ell_k(y) = \sqrt{r_k^2 - y^2} - \sqrt{r_{k-1}^2 - y^2}` (each root 0
past its radius). Since #578 it reads the outbound slots of the in-plane
chord; since the P2 fix round it reads the kernel's orbit-space shell function
:func:`~orpheus.geometry.chord.region_half_lengths` directly, with no line, from
the half-chords :math:`\sqrt{(r_k - y)(r_k + y)}`, the shells formed as
:math:`(r_k^2 - r_{k-1}^2)/(h_k + h_{k-1})`, without the cancellation of the
difference of two roots. Row C2 (the in-plane line's impact parameter was
``y`` exactly) was deleted with the line it pinned.

Before this file the only pins of the function were single-region (the P_ss
closed forms in ``tests/gates/cp/test_sphere_pss.py``), which never form a
shell difference, and the CP production-against-reference rows, which share
this function on both sides (X4). Verification spec:
``scratch/characteristic_architecture/p2/gate_spec.md`` §C (row C1).
"""
from __future__ import annotations

import mpmath as mp
import numpy as np
import pytest

from orpheus.derivations.common.kernels import chord_half_lengths

_EPS = np.finfo(float).eps

_PARTITIONS = {
    "three_regions": np.array([0.3, 1.1, 2.0]),
    "thin_shell": np.array([1.0, 1.0 + 1e-6, 2.0]),
    "two_thin_shells": np.array([0.5, 0.5000001, 0.5000002]),
}


def _draws(radii: np.ndarray) -> np.ndarray:
    """2000 uniform y on [0, R], near-tangent y = r_k (1 - 10^-u) with u in [1, 15], every r_k, and 0 (seed 578)."""
    rng = np.random.default_rng(578)
    near = (radii * (1.0 - 10.0 ** -rng.uniform(1.0, 15.0, (60, len(radii))))).ravel()
    return np.sort(np.concatenate([rng.uniform(0.0, radii[-1], 2000), near, radii, [0.0]]))


def _exact(radii: np.ndarray, y: np.ndarray) -> list[list]:
    mp.mp.dps = 50
    out = []
    inner = [0.0, *radii[:-1]]
    for r_out, r_in in zip(radii, inner):
        row = []
        for yy in y:
            Y, RO, RI = mp.mpf(float(yy)), mp.mpf(float(r_out)), mp.mpf(float(r_in))
            o = mp.sqrt(RO * RO - Y * Y) if Y < RO else mp.mpf(0)
            i = mp.sqrt(RI * RI - Y * Y) if Y < RI else mp.mpf(0)
            row.append(o - i)
        out.append(row)
    return out


@pytest.mark.l0
@pytest.mark.verifies("chord-length", "geometry-chord-segment-lengths")
@pytest.mark.rests_on("tests/gates/geometry/test_chord.py::test_the_half_chord_is_well_conditioned_at_near_tangency")
@pytest.mark.parametrize("name", list(_PARTITIONS))
def test_each_shell_of_the_half_chord_is_exact_to_a_few_ulp(name: str) -> None:
    """C1: every nonzero ``ell_k(y)`` within 8 eps relative of the 50-digit value; exactly 0 where ``y >= r_k``.

    `[M]` 2026-10-10 (``ta/c/probe_c.py``): the kernel spelling is within 4.6e-16 relative (4 ulp of the value) on
    all three partitions, 6552 entries each. First red: the spelling it replaced (HEAD ``a7aeebb5``, the difference
    of two square roots of ``r^2 - y^2``) is off by 7.6e-3 relative near tangency on the three-region partition,
    2.2e-6 on the thin shell, 7.3e-7 on the two thin shells.
    """
    radii = _PARTITIONS[name]
    y = _draws(radii)
    got = chord_half_lengths(radii, y)
    exact = _exact(radii, y)
    worst = 0.0
    for k, row in enumerate(exact):
        for m, e in enumerate(row):
            g = float(got[k, m])
            if e == 0:
                assert g == 0.0, f"region {k}, y = {y[m]!r}: a shell the line misses has length {g!r}"
            else:
                worst = max(worst, float(abs(mp.mpf(g) - e) / e))
    assert worst <= 8 * _EPS, f"worst relative error {worst:.2e} > 8 eps"
