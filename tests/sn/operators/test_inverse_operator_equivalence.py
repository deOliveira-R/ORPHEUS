r"""Phase 2 keystone — ``(L+C).inverse().apply(b) == (L+C).solve(b)`` (#226).

THE proof that ``.inverse()`` is the SAME operator as today's ``.solve``, in
operator form. :meth:`SweepOperator.apply` wraps
:meth:`InvertibleOperator.solve`, so the equality is bit-exact BY CONSTRUCTION —
this gate proves the wrapping introduced ZERO change. It is necessary-not-
sufficient on its own (it proves "unchanged", not "was right"); it rides the
existing closed-form anchors (``test_removal_form_kinf_independent_reference_2g``
+ the curvilinear keff gates) for the "was right" half (lessons L2 two-anchor).

Verification spec §2 (Gate 2.1). The curvilinear Carlson coupled-pole seed is
exercised on sphere + cyl (slab is BLIND to the seed-threading recursion — only
the curvilinear angular redistribution reads the previous-iterate ``ψ`` at
``μ=−1``), PAIRED with a Mode-11 in-process seed-capture so that a seed-
insensitive converged fixed point cannot mask a dropped ``initial_guess``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.sn.operators.streaming import InvertibleOperator, StreamingOperator
from orpheus.sn.operators.sweep_operator import SweepOperator
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from tests.sn.operators.test_removal_form_matvec_sweep import (
    _cyl,
    _random_state,
    _removal_sigmas,
    _slab,
    _sphere,
)

pytestmark = [pytest.mark.foundation, pytest.mark.verifies("inverse-as-operator")]


# slab + the two curvilinear geometries, all ≥2G heterogeneous, vacuum BC (a
# single direct solve — no reflective iteration). The curvilinear pair activates
# the Carlson seed recursion the slab cannot see (§0.6).
_CASES = {"slab_2g": _slab, "sphere_2g": _sphere, "cyl_2g": _cyl}
_SEEDED_CASES = {"sphere_2g": _sphere, "cyl_2g": _cyl}


def _build(case_builder):
    """Return ``(sn, A=L+C InvertibleOperator, b=rhs source)`` for a removal case."""
    sn = case_builder(ng=2, bc="vacuum")
    _, sig_r = _removal_sigmas(sn, seed=2026)
    A = StreamingOperator(sn) + MultiplicationOperator.from_mesh(sig_r, sn)
    assert isinstance(A, InvertibleOperator)  # L + C dispatches to the resolvent
    b = _random_state(sn, seed=11)
    return sn, A, b


def _assert_fields_equal(got, ref, msg):
    """Bit-identical on BOTH the bulk flux AND every boundary face trace."""
    np.testing.assert_array_equal(got.interior.values, ref.interior.values, err_msg=msg)
    for face in ref.boundary.layout.faces:
        np.testing.assert_array_equal(
            got.boundary.face_view(face), ref.boundary.face_view(face),
            err_msg=f"{msg} (boundary face {face})",
        )


@pytest.mark.parametrize("case", _CASES)
def test_inverse_apply_equals_solve(case):
    """``A.inverse().apply(b)`` is BIT-IDENTICAL to ``A.solve(b)``, every geometry."""
    _, A, b = _build(_CASES[case])
    _assert_fields_equal(
        A.inverse().apply(b), A.solve(b),
        f"{case}: inverse().apply diverged from solve",
    )


@pytest.mark.parametrize("case", _SEEDED_CASES)
def test_inverse_apply_accepts_and_drops_seed(case):
    """The WDD sweep is an EXACT direct inverse — nothing to seed. #282 route
    (a) computes the curvilinear ψ½ starting direction DIRECTLY from the source
    (2.5d), so #280 2.5c retired the vestigial threading: ``inverse().apply``
    ACCEPTS ``initial_guess`` (the uniform ``SupportsSeededApply`` contract,
    #285) but DROPS it. The result is seed-INDEPENDENT — bit-identical to the
    un-seeded apply. The warm start lives at the iteration layer
    (``SourceIteration`` / ``GreenOperator``), never on a direct sweep."""
    sn, A, b = _build(_SEEDED_CASES[case])
    seed = _random_state(sn, seed=99)  # an arbitrary state that MUST NOT matter
    _assert_fields_equal(
        A.inverse().apply(b, initial_guess=seed),
        A.inverse().apply(b),
        f"{case}: inverse().apply CONSUMED the initial_guess — a direct "
        f"exact inverse must drop it",
    )


def test_inverse_apply_drops_seed_before_inner_solve(monkeypatch):
    """Mode-11 path pin: ``inverse().apply`` drops the seed and enters
    ``inner.solve`` with the rhs ALONE (the exact-inverse accept-and-drop,
    #280 2.5c). The STRICT spy signature ``spy(rhs)`` would ``TypeError`` if
    ``apply`` re-forwarded ``initial_guess`` — so a resurrected thread reddens
    here, it does not pass silently."""
    sn, A, b = _build(_sphere)
    seed = _random_state(sn, seed=99)
    captured = {}
    real_solve = A.solve

    def spy(rhs):  # strict: NO initial_guess kwarg — apply must not forward one
        captured["reached"] = True
        return real_solve(rhs)

    monkeypatch.setattr(A, "solve", spy)
    A.inverse().apply(b, initial_guess=seed)
    assert captured.get("reached")  # inner.solve entered, seed dropped, no crash


def test_inverse_returns_sweep_operator_surface():
    """``(L+C).inverse()`` is a :class:`SweepOperator` wrapping the forward op.

    Pins the taxonomy step-1 surface (supersedes the Phase-2 "apply-only"
    deferral pin): the invertibility axis is now DELIVERED — ``is_invertible``
    True with a faithful ``solve`` verb (``solve`` on the inverse is the
    forward matvec), and the involution ``inverse().inverse() is A`` holds by
    OBJECT IDENTITY (§13 I2). The ADJOINT axis is now DELIVERED (#280 2.5c):
    ``is_adjointable`` True — the reverse-scan ``apply_transpose``
    ``(A⁻¹)ᵀ = (Aᵀ)⁻¹``. Domain/codomain swap (equal here — ``L+C`` is
    endomorphic)."""
    _, A, _ = _build(_slab)
    inv = A.inverse()
    assert isinstance(inv, SweepOperator)
    assert inv.inner is A
    assert callable(getattr(inv, "solve", None))  # inverse family keeps solve
    assert inv.is_invertible is True and inv.is_adjointable is True
    assert inv.inverse() is A  # (A^{-1})^{-1} = A — by object identity
    assert inv.domain is A.codomain and inv.codomain is A.domain
