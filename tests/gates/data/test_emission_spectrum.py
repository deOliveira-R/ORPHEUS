"""S10a — intrinsic-property gate for ``EmissionSpectrum``.

The fission emission spectrum :math:`\\chi` is a probability simplex:
:math:`\\sum_g \\chi_g = 1` and :math:`\\chi_g \\ge 0\\ \\forall g`.
S10a introduces ``orpheus.data.emission_spectrum.EmissionSpectrum`` — a
``numpy.ndarray`` subclass that carries this invariant as methods
(``assert_normalized`` / ``is_emitting`` / ``assert_null``).

This file is the vv-anti-pattern-#11 gate: EACH validator gets BOTH a
positive leg (a correct instance MUST NOT raise) AND a negative leg (a
broken instance MUST raise), with the simplex / non-simplex arrays
constructed BY HAND (L11 structural independence — never via the
production ``make_mixture`` / ``load_isotope`` builder, which would make
the test self-referential).

Mode-8 (vv §"Compiled-out assertion"): ORPHEUS' canonical invocation is
``python -O``, which strips bare ``assert``. Every leg below asserts via
``pytest.raises`` / ``np.testing.*`` — NEVER a bare ``assert`` — so the
gate fires under ``-O``.

These are software-invariant tests (the simplex law is a probability-
distribution invariant, not a theory-page ``:label:`` equation), so the
module is ``@pytest.mark.foundation`` and carries NO ``verifies(...)``.

Tolerance design (probe done in the S10a brief): real GENDF
``load_isotope('U_235', 294)`` gives ``chi.sum() == 1.0000000000000002``
(residual 2.22e-16). ``assert_normalized`` uses
``np.isclose(self.sum(), 1.0, atol=1e-12, rtol=0)`` — ~4 orders of
headroom over the observed FP residual, while still RED on a 1e-6
physical error. The edge legs below pin both ends of that band.
"""

from __future__ import annotations

import numpy as np
import pytest

# NOTE (pre-impl): this import does not resolve until the method-implementer
# lands ``orpheus/data/emission_spectrum.py``. The whole module is guarded
# so the suite collects (skipped) pre-impl and goes green post-impl with
# zero edits here.
emission_spectrum = pytest.importorskip("orpheus.data.emission_spectrum")
EmissionSpectrum = emission_spectrum.EmissionSpectrum

pytestmark = pytest.mark.foundation


def _require(condition: bool, message: str) -> None:
    """``-O``-safe assertion (``pytest.fail`` is a call, not a stripped assert).

    Mode-8: a bare ``assert`` is a no-op under ``python -O``. Boolean /
    isinstance / shape legs route through this helper so they FIRE under
    the canonical ``-O`` invocation. (``np.testing.*`` legs already do.)
    """
    if not condition:
        pytest.fail(message)


# ══════════════════════════════════════════════════════════════════════
# ndarray-subclass behaviour (gate 3a — zero-ripple at the type level)
# ══════════════════════════════════════════════════════════════════════


class TestNdarraySubclassBehaviour:
    """EmissionSpectrum must behave as an ndarray in every existing site.

    The headline S10a claim is "the type wraps the values; the values do
    not move". These legs pin the ndarray-subclass contract that lets the
    existing ``chi[None, :, :]`` / ``chi.sum()`` / ``chi.copy()`` /
    einsum-feeding call sites keep working with zero ripple. They also
    pin that ``__array_finalize__`` was implemented (slicing / copy /
    arithmetic must return a usable array, not crash).
    """

    def test_is_ndarray_subclass(self):
        chi = EmissionSpectrum(np.array([0.6, 0.35, 0.05]))
        _require(isinstance(chi, np.ndarray), "EmissionSpectrum must subclass ndarray")

    def test_asarray_roundtrips_values(self):
        """np.asarray(EmissionSpectrum(a)) == a — coercion preserves values."""
        a = np.array([0.6, 0.35, 0.05])
        chi = EmissionSpectrum(a)
        np.testing.assert_array_equal(np.asarray(chi), a)

    def test_sum_returns_correct_value(self):
        chi = EmissionSpectrum(np.array([0.6, 0.35, 0.05]))
        # 0.6+0.35+0.05 == 1.0 up to FP; assert_allclose tolerant of ULP.
        np.testing.assert_allclose(chi.sum(), 1.0, atol=1e-12, rtol=0)

    def test_newaxis_broadcast_shape(self):
        """chi[None, :, :] is the per-cell broadcast feed — must not crash."""
        chi = EmissionSpectrum(np.array([0.6, 0.4]))
        broadcast = chi[None, :]
        _require(broadcast.shape == (1, 2), f"unexpected broadcast shape {broadcast.shape}")

    def test_copy_preserves_values(self):
        a = np.array([0.6, 0.35, 0.05])
        chi = EmissionSpectrum(a)
        np.testing.assert_array_equal(chi.copy(), a)

    def test_einsum_feed_returns_correct_value(self):
        """einsum over (chi, nu_sigma_f) — the FissionOperator feed path."""
        chi = EmissionSpectrum(np.array([1.0, 0.0]))
        nu_sigma_f = np.array([0.01, 0.08])
        phi = np.array([2.0, 3.0])
        # production rate p = sum_g nu_sigma_f_g * phi_g
        rate = np.einsum("g,g->", nu_sigma_f, phi)
        emission = chi * rate  # chi broadcast * scalar rate
        np.testing.assert_allclose(emission, np.array([0.26, 0.0]), rtol=0, atol=1e-15)


# ══════════════════════════════════════════════════════════════════════
# assert_normalized — simplex/probability law (gate 1, both legs)
# ══════════════════════════════════════════════════════════════════════


class TestAssertNormalized:
    """``assert_normalized`` validates Σχ ≈ 1 AND χ_g ≥ 0 ∀g.

    Both clauses are tested SEPARATELY: the sum clause (negative leg
    [0.6,0.35,0.10] → Σ=1.05) and the negativity clause
    ([1.2,-0.2] → Σ=1 but a negative entry) catch independent
    failure modes. Testing only one would let the other regress.
    """

    # ── positive leg (MUST NOT raise) ──────────────────────────────────

    def test_simplex_passes(self):
        chi = EmissionSpectrum(np.array([0.6, 0.35, 0.05]))
        chi.assert_normalized()  # MUST NOT raise

    # ── negative leg: sum clause ───────────────────────────────────────

    def test_sum_over_one_raises(self):
        """Σ = 1.05 > 1 — violates the simplex sum clause."""
        chi = EmissionSpectrum(np.array([0.6, 0.35, 0.10]))
        with pytest.raises((ValueError, AssertionError)):
            chi.assert_normalized()

    # ── negative leg: negativity clause (INDEPENDENT of sum) ───────────

    def test_negative_entry_raises_even_when_sum_is_one(self):
        """[1.2, -0.2] sums to 1 but has a negative entry — must raise.

        This is the clause-separation leg: a test that only checked the
        sum would PASS this array (Σ=1.0) while the spectrum is not a
        valid probability distribution.
        """
        chi = EmissionSpectrum(np.array([1.2, -0.2]))
        with pytest.raises((ValueError, AssertionError)):
            chi.assert_normalized()

    # ── tolerance-edge legs (pin the FP band) ──────────────────────────

    def test_one_ulp_residual_passes(self):
        """[0.5, 0.5+2e-16] sums to ~1+2e-16 — inside the 1e-12 band.

        Mirrors the real-GENDF U_235 residual (2.22e-16). MUST NOT raise.
        """
        chi = EmissionSpectrum(np.array([0.5, 0.5 + 2e-16]))
        chi.assert_normalized()  # MUST NOT raise

    def test_off_by_1e6_raises(self):
        """A 1e-6 physical error is OUTSIDE the band — MUST raise.

        Pins the upper end: the tolerance is loose enough for FP residual
        but tight enough to catch a real normalization bug.
        """
        chi = EmissionSpectrum(np.array([0.5, 0.5 + 1e-6]))
        with pytest.raises((ValueError, AssertionError)):
            chi.assert_normalized()


# ══════════════════════════════════════════════════════════════════════
# is_emitting — boolean property (gate 1, both legs)
# ══════════════════════════════════════════════════════════════════════


class TestIsEmitting:
    """``is_emitting`` is True iff any χ_g > 0."""

    def test_nonzero_is_emitting(self):
        # ``is_emitting`` is specified as ``bool(np.any(...))`` — a Python
        # ``bool``. Pin that exactly (an ``np.bool_`` leak would fail the
        # ``type`` check, catching a missing ``bool(...)`` coercion).
        chi = EmissionSpectrum(np.array([1.0, 0.0]))
        _require(chi.is_emitting is True, "nonzero spectrum must be emitting")
        _require(type(chi.is_emitting) is bool, "is_emitting must be a Python bool")

    def test_all_zero_is_not_emitting(self):
        chi = EmissionSpectrum(np.zeros(3))
        _require(chi.is_emitting is False, "all-zero spectrum must not be emitting")
        _require(type(chi.is_emitting) is bool, "is_emitting must be a Python bool")


# ══════════════════════════════════════════════════════════════════════
# assert_null — strict χ ≡ 0 (gate 1, both legs)
# ══════════════════════════════════════════════════════════════════════


class TestAssertNull:
    """``assert_null`` validates χ ≡ 0 (every entry exactly zero).

    This is the non-fissile contract: a non-fissile material emits no
    fission neutrons, so its spectrum is identically zero (STRICT — not
    "approximately zero"; there is no physical residual to tolerate).
    """

    def test_zeros_pass(self):
        chi = EmissionSpectrum(np.zeros(4))
        chi.assert_null()  # MUST NOT raise

    def test_any_nonzero_raises(self):
        chi = EmissionSpectrum(np.array([0.0, 0.0, 1e-9, 0.0]))
        with pytest.raises((ValueError, AssertionError)):
            chi.assert_null()
