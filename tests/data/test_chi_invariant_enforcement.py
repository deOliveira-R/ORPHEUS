"""S10a — container χ-invariant enforcement gate.

Both ``Mixture`` (``orpheus.data.macro_xs.mixture``) and ``Isotope``
(``orpheus.data.micro_xs.isotope``) gain a ``__post_init__`` that routes
``self.chi`` through ``enforce_emission_spectrum`` — coercing it to an
``EmissionSpectrum`` and running the STRICT conditional cross-check keyed
on PRODUCTION::

    enforce_emission_spectrum(chi, is_producing=self.is_producing)
    #   is_producing → chi.assert_normalized()   (simplex)
    #   not producing → chi.assert_null()         (χ ≡ 0, strict)

``is_producing`` is the emission predicate (``νΣ_f > 0``): ``Mixture``
keys it on ``SigP`` (the production XS), ``Isotope`` on ``nubar * sigF``.
χ is consumed only as a fission SOURCE (``χ·νΣ_f·φ``), so a valid simplex
is required exactly where production is nonzero. (Distinct from
``Isotope.is_fissile = np.any(sigF > 0)``, the "can it fission?" predicate
consumed by ``compute_macro_xs`` — left unchanged.)

This file is the vv-anti-pattern-#11 gate at the CONTAINER level: for
EACH container, BOTH a positive leg (correct construction → ok) AND a
negative leg (broken construction → raises). Minimal fixtures are
built BY HAND (L11 — never via the production ``make_mixture`` /
``compute_macro_xs`` / ``load_isotope`` builder for the synthetic legs;
those builders ARE the system the guard protects and reusing them would
make the test self-referential).

The real-GENDF production-path proof (``TestRealGendfConstructs``) is the
flip side: it uses the REAL builder ``load_isotope`` precisely to prove
the invariant is a GUARD (passes on valid real data), not a BREAKER.

Mode-8: every assertion uses ``pytest.raises`` / ``np.testing.*`` (the
suite runs ``-O``; bare ``assert`` is stripped).

Foundation tests (software invariant, no theory ``:label:``).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.data.micro_xs.isotope import NG, Isotope

pytestmark = pytest.mark.foundation


# ══════════════════════════════════════════════════════════════════════
# Hand-built minimal fixtures (L11 structural independence)
# ══════════════════════════════════════════════════════════════════════


def _mixture(*, sig_p: np.ndarray, chi: np.ndarray) -> Mixture:
    """Minimal 2-group ``Mixture`` with explicit SigP (production) and chi.

    The χ law keys on PRODUCTION: ``is_producing = np.any(self.SigP > 0)``,
    so the two arrays that matter for the cross-check are ``SigP`` (the
    gate) and ``chi`` (the validated value). A producing leg passes
    ``sig_p=[...>0]``; a non-producing leg passes ``sig_p=zeros``. Built
    directly (NOT via ``make_mixture``) so the test does not depend on the
    builder it guards.
    """
    ng = len(sig_p)
    zeros = np.zeros(ng)
    return Mixture(
        SigC=zeros.copy(),
        SigL=zeros.copy(),
        SigF=zeros.copy(),
        SigP=np.asarray(sig_p, dtype=float),
        SigT=np.full(ng, 1.0),
        SigS=[csr_matrix((ng, ng))],
        Sig2=csr_matrix((ng, ng)),
        chi=np.asarray(chi, dtype=float),
        eg=None,
    )


def _isotope(*, sig_f_row: np.ndarray, nubar: np.ndarray, chi: np.ndarray) -> Isotope:
    """Minimal ``Isotope`` with explicit sigF, nubar and chi.

    The χ law keys on PRODUCTION: ``is_producing = np.any(nubar * sigF > 0)``.
    A producing leg sets ``nubar > 0`` where ``sig_f_row > 0``; a
    non-producing leg passes ``nubar=zeros`` (or ``sig_f_row=zeros``). Built
    directly (NOT via ``load_isotope``). ``NG`` is the hard-coded library
    size; ``chi`` is full-length ``(NG,)`` so the simplex sum is over the
    real group count.
    """
    n_sig0 = 1
    return Isotope(
        name="TEST",
        aw=1.0,
        temp=294.0,
        eg=np.linspace(2e7, 1e-5, NG + 1),
        sig0=np.array([1e10]),
        sigC=np.zeros((n_sig0, NG)),
        sigL=np.zeros((n_sig0, NG)),
        sigF=np.asarray(sig_f_row, dtype=float).reshape(n_sig0, NG),
        sigT=np.full((n_sig0, NG), 1.0),
        nubar=np.asarray(nubar, dtype=float),
        chi=np.asarray(chi, dtype=float),
        sigS=[[csr_matrix((NG, NG))]],
        sig2=csr_matrix((NG, NG)),
    )


def _simplex_ng() -> np.ndarray:
    """A valid (NG,)-length simplex: all mass in group 0."""
    chi = np.zeros(NG)
    chi[0] = 1.0
    return chi


# ── ``-O``-safe assertion helpers (Mode-8) ─────────────────────────────
#
# The suite runs ``python -O``, which strips bare ``assert``. ``pytest.fail``
# is a function call, so it FIRES under ``-O`` — use it for the positive
# legs' substance (the fissile-state + coercion checks) so they are not
# reduced to "did the constructor throw?".


def _require(condition: bool, message: str) -> None:
    if not condition:
        pytest.fail(message)


# ══════════════════════════════════════════════════════════════════════
# Mixture cross-check (gate 2 — both legs, both branches)
# ══════════════════════════════════════════════════════════════════════


class TestMixturePostInit:
    """Mixture.__post_init__ enforces the conditional χ cross-check.

    The law keys on PRODUCTION (``is_producing = np.any(self.SigP > 0)``):
    a producing mixture carries a simplex, a non-producing one carries the
    null spectrum.
    """

    # ── producing branch ────────────────────────────────────────────────

    def test_producing_simplex_constructs(self):
        """Producing (SigP > 0) + simplex χ → constructs ok (MUST NOT raise)."""
        m = _mixture(sig_p=np.array([0.05, 0.12]), chi=np.array([1.0, 0.0]))
        _require(bool(m.is_producing), "producing mixture must report is_producing")

    def test_producing_non_simplex_raises(self):
        """Producing + non-simplex χ (Σ=1.5) → raises."""
        with pytest.raises((ValueError, AssertionError)):
            _mixture(sig_p=np.array([0.05, 0.12]), chi=np.array([1.0, 0.5]))

    # ── non-producing branch ────────────────────────────────────────────

    def test_non_producing_zeros_constructs(self):
        """Non-producing (SigP ≡ 0) + χ ≡ 0 → constructs ok (MUST NOT raise)."""
        m = _mixture(sig_p=np.zeros(2), chi=np.zeros(2))
        _require(not bool(m.is_producing), "non-producing mixture must report not is_producing")

    def test_non_producing_nonzero_raises(self):
        """Non-producing + nonzero χ → raises (the strict assert_null clause).

        This is the leg that the 77 ``get_mixture("B"/"C"/"D")`` call
        sites would trip TODAY (region B/C/D carry chi=[1.0] / [1,0] /
        [.6,.35,.05,0]). See the gate spec §Precursor — the xs_library
        regions MUST be zeroed before this branch is enabled in
        production, or the SN/CP/MoC suites red en masse.
        """
        with pytest.raises((ValueError, AssertionError)):
            _mixture(sig_p=np.zeros(2), chi=np.array([1.0, 0.0]))

    # ── coercion leg ───────────────────────────────────────────────────

    def test_chi_coerced_to_emission_spectrum(self):
        """After construction, chi is an EmissionSpectrum (type wrap)."""
        es = pytest.importorskip("orpheus.data.emission_spectrum")
        m = _mixture(sig_p=np.array([0.05, 0.12]), chi=np.array([1.0, 0.0]))
        _require(isinstance(m.chi, es.EmissionSpectrum), "Mixture.chi must coerce to EmissionSpectrum")


# ══════════════════════════════════════════════════════════════════════
# Isotope cross-check (gate 2 — both legs, both branches)
# ══════════════════════════════════════════════════════════════════════


class TestIsotopePostInit:
    """Isotope.__post_init__ enforces the conditional χ cross-check.

    The law keys on PRODUCTION (``is_producing = np.any(nubar * sigF > 0)``):
    a producing isotope (``nubar > 0`` where ``sigF > 0``) carries a
    simplex, a non-producing one carries the null spectrum.
    """

    @staticmethod
    def _producing_nubar(row: np.ndarray) -> np.ndarray:
        """nubar > 0 exactly where the fission row is nonzero (so νΣ_f > 0)."""
        nubar = np.zeros(NG)
        nubar[row > 0] = 2.5
        return nubar

    def test_producing_simplex_constructs(self):
        row = np.zeros(NG)
        row[0] = 0.3
        iso = _isotope(sig_f_row=row, nubar=self._producing_nubar(row), chi=_simplex_ng())
        _require(bool(iso.is_producing), "producing isotope must report is_producing")

    def test_producing_non_simplex_raises(self):
        row = np.zeros(NG)
        row[0] = 0.3
        bad = _simplex_ng()
        bad[1] = 0.5  # now Σ = 1.5
        with pytest.raises((ValueError, AssertionError)):
            _isotope(sig_f_row=row, nubar=self._producing_nubar(row), chi=bad)

    def test_non_producing_zeros_constructs(self):
        iso = _isotope(sig_f_row=np.zeros(NG), nubar=np.zeros(NG), chi=np.zeros(NG))
        _require(not bool(iso.is_producing), "non-producing isotope must report not is_producing")

    def test_non_producing_nonzero_raises(self):
        """Non-producing (νΣ_f ≡ 0 via nubar=0) + nonzero χ → raises.

        sigF may be nonzero (it can fission) but with ``nubar=0`` it emits
        nothing — production is zero, so χ must be null. This isolates the
        production predicate from ``is_fissile`` (which would call this
        material fissile on the sigF row).
        """
        row = np.zeros(NG)
        row[0] = 0.3
        with pytest.raises((ValueError, AssertionError)):
            _isotope(sig_f_row=row, nubar=np.zeros(NG), chi=_simplex_ng())

    def test_chi_coerced_to_emission_spectrum(self):
        es = pytest.importorskip("orpheus.data.emission_spectrum")
        row = np.zeros(NG)
        row[0] = 0.3
        iso = _isotope(sig_f_row=row, nubar=self._producing_nubar(row), chi=_simplex_ng())
        _require(isinstance(iso.chi, es.EmissionSpectrum), "Isotope.chi must coerce to EmissionSpectrum")


# ══════════════════════════════════════════════════════════════════════
# Real-GENDF production-path proof (gate 5)
# ══════════════════════════════════════════════════════════════════════

_MICRO_DIR = Path(__file__).resolve().parents[2] / "orpheus" / "data" / "micro_xs"


@pytest.mark.skipif(
    not (_MICRO_DIR / "U_235.h5").exists(),
    reason="GENDF .h5 library not present",
)
class TestRealGendfConstructs:
    """The invariant is a GUARD on real data, not a BREAKER.

    Uses the REAL ``load_isotope`` builder (deliberately — this is the
    production-path no-red proof; the synthetic legs above are the L11
    hand-built leg). The 2.22e-16 residual probed for U_235 must sit
    inside the 1e-12 band, so the fissile branch must NOT raise.
    """

    def test_u235_producing_constructs_without_raising(self):
        from orpheus.data.micro_xs import load_isotope

        iso = load_isotope("U_235", 294)  # MUST NOT raise post-impl
        _require(bool(iso.is_producing), "U_235 must report is_producing (νΣ_f > 0)")
        # Belt-and-braces: the real residual is inside the band.
        np.testing.assert_allclose(iso.chi.sum(), 1.0, atol=1e-12, rtol=0)

    @pytest.mark.parametrize("name", ["O_016", "H_001"])
    def test_non_producing_isotope_constructs_without_raising(self, name):
        """A non-producing real isotope (νΣ_f ≡ 0, χ ≡ 0) constructs ok.

        Proves the non-producing assert_null branch does not red on real
        non-producing GENDF data (whose nubar·sigF ≡ 0 and chi is
        genuinely zero).
        """
        from orpheus.data.micro_xs import load_isotope

        iso = load_isotope(name, 294)  # MUST NOT raise post-impl
        _require(not bool(iso.is_producing), f"{name} must report not is_producing")
