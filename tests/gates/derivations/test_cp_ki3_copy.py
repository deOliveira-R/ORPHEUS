r"""G5 of the unit ``insulation`` (#641 item 1): production's Bickley-Naylor kernel ``orpheus.numerics.bickley.ki3`` and the reference's.

``Ki_3(tau) = integral_0^{pi/2} cos^2(theta) exp(-tau / cos(theta)) d(theta)``, evaluated as ``exp(-tau)`` times the degree-63 Chebyshev
interpolant of ``exp(tau) Ki_3(tau)`` on ``[0, 50]``, and clamped to 0 above 50. Production holds a COPY of the interpolant that
``orpheus.derivations.continuous.flat_source_cp.geometry`` holds (``_ki3_mp`` over ``_ki3_scaled_cheb``), rebuilt from the integral at 30
digits, never stored (a stored table is an untraceable constant). This file lives under ``tests/gates/derivations/`` because
``FORBIDDEN_EDGES["numerics"]`` forbids ``derivations``: a numerics-directory test could not import the reference.

What it proves: (a) and (b), the two copies agree bit for bit (the same 64 Chebyshev points of the first kind, the same coefficients);
(c) the copy is accurate to 3e-6 against the DEFINING integral written here (independent of the reference copy, which this file does not use
for the value); (d) the clamp. What it does NOT prove: independence of the approximation. Two copies of one 2.6e-6 approximation drift
together, so legs (a) and (b) say nothing about accuracy (X4); the independent verification chain for the kernel is #641 item 2, filed
separately.

First red at the base tip: ``ModuleNotFoundError`` for ``orpheus.numerics.bickley`` (imported at call time, so every row reds on its own).
Mutations: M-degree-62 (``_KI3_DEG`` patched to 62 and the ``lru_cache`` of ``ki3_scaled_chebyshev`` cleared, AND cleared again on teardown
after ``monkeypatch.undo()``: otherwise every later cylinder solve of the session reads the degree-62 kernel) reddens (a) and (b); leg (c)
reddens when the mutant is accurate to less than 3e-6, which the row reports. The accuracy checker's own control is ``test_the_integral_check_sees_a_wrong_kernel``,
green today: the reference passes it and ``exp(-tau)`` does not.
"""

from __future__ import annotations

import mpmath
import numpy as np
import pytest

from orpheus.derivations.continuous.flat_source_cp import geometry as reference

_GRID = np.concatenate([np.linspace(0.0, 50.0, 2000), [50.0, 50.0000001, 1e-300, 60.0]])
_ACCURACY_POINTS = np.concatenate([[0.0], np.logspace(-8, np.log10(49.9), 80)])
_RELATIVE_TOLERANCE = 3e-6


def _bickley():
    from orpheus.numerics import bickley

    return bickley


def _defining_integral(t: float) -> float:
    r""":math:`\mathrm{Ki}_3(t) = \int_0^\infty (1+u^2)^{-2} e^{-t\sqrt{1+u^2}}\,du` (the tan substitution of the defining integral), 30 digits."""
    with mpmath.workdps(30):
        return float(mpmath.quad(lambda u: (1 + u * u) ** -2 * mpmath.exp(-t * mpmath.sqrt(1 + u * u)), [0, 1, mpmath.inf]))


def _worst_relative_error(kernel) -> float:
    exact = np.array([_defining_integral(float(t)) for t in _ACCURACY_POINTS])
    return float(np.max(np.abs(kernel(_ACCURACY_POINTS) - exact) / exact))


@pytest.mark.l1
def test_the_production_kernel_equals_the_reference_kernel_bit_for_bit() -> None:
    """(a) ``ki3`` equals ``_ki3_mp`` on a 2004-point grid (``linspace(0, 50, 2000)`` and the ends ``50``, ``50.0000001``, ``1e-300``, ``60``).
    The two are one interpolant, so this says the copy is the copy; it says nothing about accuracy (leg (c) does)."""
    got = _bickley().ki3(_GRID)
    assert got.shape == _GRID.shape and np.array_equal(got, reference._ki3_mp(_GRID))


@pytest.mark.l1
def test_the_production_interpolant_has_the_reference_coefficients() -> None:
    """(b) the 64 Chebyshev coefficients and the domain are equal bit for bit, and the cached object's coefficients are read-only."""
    poly, theirs = _bickley().ki3_scaled_chebyshev(), reference._ki3_scaled_cheb()
    assert poly.coef.shape == (64,) and np.array_equal(poly.coef, theirs.coef) and np.array_equal(poly.domain, theirs.domain)
    assert not poly.coef.flags.writeable, "a mutable cached interpolant would let one assignment change every later kernel call"


@pytest.mark.l1
def test_the_production_kernel_is_accurate_against_the_defining_integral() -> None:
    """(c) the relative error against the integral written in this file is at most 3e-6 on the 81 points
    ``[0] + logspace(-8, log10(49.9), 80)`` (measured maximum 2.6e-6, at tau = 0). Independent of the reference copy."""
    assert _worst_relative_error(_bickley().ki3) <= _RELATIVE_TOLERANCE


@pytest.mark.l1
def test_the_kernel_is_zero_above_the_clamp_and_empty_in_empty_out() -> None:
    """(d) ``ki3(60) == 0.0`` exactly (``Ki_3(50)`` is about 3e-23) and the output of an empty input is empty."""
    ki3 = _bickley().ki3
    assert ki3(60.0) == 0.0 and ki3(np.array([50.0000001, 1e3])).tolist() == [0.0, 0.0]
    assert ki3(np.array([])).shape == (0,)


@pytest.mark.foundation
def test_the_integral_check_sees_a_wrong_kernel() -> None:
    """Control of leg (c), green at the base tip: the reference interpolant passes the 3e-6 check, and ``exp(-tau)`` (the sphere's kernel,
    which is not ``Ki_3``) fails it, so the check separates a right kernel from a wrong one and is not satisfied by any input."""
    assert _worst_relative_error(reference._ki3_mp) <= _RELATIVE_TOLERANCE
    assert _worst_relative_error(lambda t: np.exp(-t)) > 1e-2


@pytest.mark.l1
def test_a_degree_62_interpolant_is_seen_by_legs_a_and_b(monkeypatch: pytest.MonkeyPatch) -> None:
    r"""M-degree-62: the degree is patched and the cache cleared; the grid reads differ at (nearly) every point and the coefficient
    arrays differ. The cache is cleared AGAIN after ``monkeypatch.undo()``, in a ``finally``."""
    module = _bickley()
    honest = module.ki3(_GRID)
    try:
        monkeypatch.setattr(module, "_KI3_DEG", 62)
        module.ki3_scaled_chebyshev.cache_clear()
        mutant = module.ki3(_GRID)
        assert np.count_nonzero(mutant != reference._ki3_mp(_GRID)) > 1000
        assert not np.array_equal(module.ki3_scaled_chebyshev().coef, reference._ki3_scaled_cheb().coef)
    finally:
        monkeypatch.undo()
        module.ki3_scaled_chebyshev.cache_clear()
    assert np.array_equal(module.ki3(_GRID), honest)
