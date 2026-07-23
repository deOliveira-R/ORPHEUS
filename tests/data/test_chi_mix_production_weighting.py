"""S10b — production-weighted χ_mix gate.

``compute_macro_xs`` (``orpheus.data.macro_xs.mixture``) builds the
mixture fission spectrum :math:`\\chi_{\\mathrm{mix}}` as the
**production-weighted convex average** over ALL fissile isotopes,
replacing the legacy first-fissile-isotope shortcut::

    prod_i  = aDen[i] * (isotopes[i].nubar * sigF[i]).sum()   # flat-flux production
    w_i     = prod_i / sum(prod_j over fissile_indices)        # Σ w_i = 1
    chi_mix = sum(w_i * isotopes[i].chi for i in fissile_indices)

────────────────────────────────────────────────────────────────────────
HONEST SCOPE (vv Mode-7 — flat-flux representative weighting)
────────────────────────────────────────────────────────────────────────
This is NOT the truly-correct material spectrum. The exact spectral
weight is FLUX-weighted (``w_i ∝ Σ_g νΣf_{i,g} φ_g``); a static
data-prep value cannot capture φ. These gates pin the **flat-flux
production-weighted convex average** (``w_i ∝ aDen_i · Σ_g νΣf_{i,g}``,
``φ_g ≡ 1``), NOT flux-exactness. This is a DOCUMENTED approximation,
not a separate issue. The verification claim is therefore: "χ_mix is
the correct flat-flux production-weighted convex average of the
per-isotope simplices", and nothing about the eigenvalue-correctness
of the resulting spectrum.

────────────────────────────────────────────────────────────────────────
WHY THE WEIGHT MATTERS (the discriminator, gate 2)
────────────────────────────────────────────────────────────────────────
For a 2-fissile mixture with χ_0 ≠ χ_1 and non-trivial weights, the
production-weighted answer differs from BOTH:
  • the legacy first-isotope shortcut (χ_0), AND
  • a naive unweighted mean ((χ_0 + χ_1)/2).
Gate 2 picks prod_0 ≠ prod_1 and χ_0 ≠ χ_1 so all three candidate
answers are pairwise distinct → the test mutation-pins the WEIGHTING.

────────────────────────────────────────────────────────────────────────
STRUCTURAL INDEPENDENCE (L11)
────────────────────────────────────────────────────────────────────────
The reference χ_mix is hand-laid term-by-term (explicit scalar
productions p0/p1, NOT a loop re-spelling the code's ``sum(w_i*χ_i)``).
The SUT is driven by the REAL ``compute_macro_xs`` on hand-built
SYNTHETIC single-σ0 ``Isotope`` objects (single-σ0 takes the cheap
identity paths in interpolation/sigma_zeros, so sigF/nubar/chi pass
through verbatim). Reference and SUT share only numpy primitives
(below the trusted-library line), never the production χ-mix expression.

The simplex-preservation interlock (gate 1) REUSES S10a's
``EmissionSpectrum.assert_normalized`` verbatim — a convex combination
of simplices is a simplex, so the S10a source-guard fires for free on
the S10b output.

Mode-8: every assertion uses ``np.testing.*`` / ``pytest.raises`` /
``_require`` (the suite runs ``-O``; bare ``assert`` is stripped).

Foundation tests (software invariants). The χ-mix convex-average law
carries a theory label (``emission-spectrum-chi-mix``, backfilled by
#231 Phase G) and the hand-reference class is wired to it via
``verifies`` — a foundation gate that genuinely pins an equation's
content carries the mark (see the harness chapter's audit-reporting
rules).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from orpheus.data.emission_spectrum import EmissionSpectrum
from orpheus.data.macro_xs.mixture import compute_macro_xs
from orpheus.data.micro_xs.isotope import NG, Isotope
from scipy.sparse import csr_matrix

pytestmark = pytest.mark.foundation


# ══════════════════════════════════════════════════════════════════════
# ``-O``-safe assertion helper (Mode-8)
# ══════════════════════════════════════════════════════════════════════


def _require(condition: bool, message: str) -> None:
    if not condition:
        pytest.fail(message)


# ══════════════════════════════════════════════════════════════════════
# Hand-built minimal synthetic isotopes (L11 structural independence)
# ══════════════════════════════════════════════════════════════════════
#
# Single-σ0 NG-group isotopes. ``compute_macro_xs`` drives them through
# the REAL pipeline but the single-σ0 branches are identity copies
# (interpolation.py:31/64 → ``field[0].copy()``; sigma_zeros.py:78 →
# ``sigT[0, ig]``), so the per-isotope sigF/nubar/chi the test controls
# survive verbatim into the χ_mix computation. ``Isotope.__post_init__``
# (S10a) enforces the simplex/null law on construction, so a fissile
# synthetic isotope MUST carry a valid (NG,) simplex χ.
#
# ⚠ LOAD-BEARING: these synthetic isotopes carry a single Legendre
# scattering order (``sigS=[[csr_matrix]]``), so EVERY ``compute_macro_xs``
# drive MUST pass ``n_legendre=1`` (the default 3 would index
# ``iso.sigS[1]`` → IndexError). PROVEN: with ``n_legendre=1`` the drive
# succeeds and the pre-S10b HEAD returns χ_mix == χ_0 (the legacy
# shortcut), distinct from both the S10b weighted answer and the
# unweighted mean — so gate 2 is RED pre-S10b and goes GREEN exactly when
# the production-weighting lands.


def _chi_simplex(*groups_and_mass: tuple[int, float]) -> np.ndarray:
    """Build a valid (NG,) simplex with mass placed in named groups.

    e.g. ``_chi_simplex((0, 0.7), (1, 0.3))`` → mass 0.7 in g0, 0.3 in g1.
    The masses MUST sum to 1.0 (the S10a constructor guard enforces it).
    """
    chi = np.zeros(NG)
    for g, m in groups_and_mass:
        chi[g] = m
    return chi


def _fissile_isotope(
    *,
    nubar_val: float,
    sigf_val: float,
    sigf_group: int,
    chi: np.ndarray,
) -> Isotope:
    """A producing synthetic isotope: nubar*sigF > 0 in one group, simplex χ.

    ``(nubar * sigF).sum() == nubar_val * sigf_val`` (single nonzero group),
    which is the per-isotope flat-flux production density the weight uses.
    """
    sigF = np.zeros((1, NG))
    sigF[0, sigf_group] = sigf_val
    nubar = np.zeros(NG)
    nubar[sigf_group] = nubar_val
    return Isotope(
        name="FISSILE",
        aw=1.0,
        temp=294.0,
        eg=np.linspace(2e7, 1e-5, NG + 1),
        sig0=np.array([1e10]),
        sigC=np.zeros((1, NG)),
        sigL=np.zeros((1, NG)),
        sigF=sigF,
        sigT=np.full((1, NG), 1.0),
        nubar=nubar,
        chi=np.asarray(chi, dtype=float),
        sigS=[[csr_matrix((NG, NG))]],
        sig2=csr_matrix((NG, NG)),
    )


def _inert_isotope() -> Isotope:
    """A non-fissile synthetic isotope: sigF ≡ 0, null χ (S10a contract)."""
    return Isotope(
        name="INERT",
        aw=1.0,
        temp=294.0,
        eg=np.linspace(2e7, 1e-5, NG + 1),
        sig0=np.array([1e10]),
        sigC=np.zeros((1, NG)),
        sigL=np.zeros((1, NG)),
        sigF=np.zeros((1, NG)),
        sigT=np.full((1, NG), 1.0),
        nubar=np.zeros(NG),
        chi=np.zeros(NG),
        sigS=[[csr_matrix((NG, NG))]],
        sig2=csr_matrix((NG, NG)),
    )


# ══════════════════════════════════════════════════════════════════════
# Gate 1 — convex-simplex intrinsic property (L11, vv #11 interlock)
# ══════════════════════════════════════════════════════════════════════


class TestChiMixIsSimplex:
    """χ_mix is a probability simplex for ANY multi-fissile mixture.

    A convex combination of simplices is a simplex (Σ=1, χ≥0). Verified
    for several weight vectors (equal, skewed, near-degenerate) so the
    claim is not an accident of one weight. The S10a interlock leg
    confirms the S10b output is a legal ``Mixture.chi`` — the same
    ``assert_normalized`` call ``Mixture.__post_init__`` already runs.
    """

    @pytest.mark.parametrize(
        "aden, sigf0, sigf1",
        [
            (np.array([1.0, 1.0, 1.0]), 0.10, 0.10),   # equal production → w=[.5,.5]
            (np.array([3.0, 1.0, 1.0]), 0.10, 0.10),   # skewed by density → w=[.75,.25]
            (np.array([1.0, 1.0, 1.0]), 0.20, 0.001),  # near-degenerate → w≈[1,0]
        ],
    )
    def test_chi_mix_sums_to_one_and_nonnegative(self, aden, sigf0, sigf1):
        iso0 = _fissile_isotope(
            nubar_val=2.4, sigf_val=sigf0, sigf_group=0,
            chi=_chi_simplex((0, 0.7), (1, 0.3)),
        )
        iso1 = _fissile_isotope(
            nubar_val=2.8, sigf_val=sigf1, sigf_group=2,
            chi=_chi_simplex((2, 0.4), (3, 0.6)),
        )
        mix = compute_macro_xs(
            [iso0, iso1, _inert_isotope()], aden,
            n_legendre=1, fissile_indices=[0, 1],
        )
        chi_mix = np.asarray(mix.chi)
        np.testing.assert_allclose(chi_mix.sum(), 1.0, atol=1e-12, rtol=0)
        _require(bool(np.all(chi_mix >= 0)), f"χ_mix has negative entries: {chi_mix!r}")

    def test_chi_mix_passes_s10a_normalization_interlock(self):
        iso0 = _fissile_isotope(
            nubar_val=2.4, sigf_val=0.10, sigf_group=0,
            chi=_chi_simplex((0, 0.7), (1, 0.3)),
        )
        iso1 = _fissile_isotope(
            nubar_val=2.8, sigf_val=0.20, sigf_group=2,
            chi=_chi_simplex((2, 0.4), (3, 0.6)),
        )
        mix = compute_macro_xs(
            [iso0, iso1, _inert_isotope()], np.array([3.0, 1.0, 1.0]),
            n_legendre=1, fissile_indices=[0, 1],
        )
        # The same guard Mixture.__post_init__ runs (enforce_emission_spectrum):
        # asserted explicitly here so the gate-1↔S10a interlock is a named,
        # inspectable claim. MUST NOT raise.
        EmissionSpectrum(np.asarray(mix.chi)).assert_normalized()


# ══════════════════════════════════════════════════════════════════════
# Gate 2 — 2-fissile hand-reference (THE TEETH, L11)
# ══════════════════════════════════════════════════════════════════════


@pytest.mark.verifies("emission-spectrum-chi-mix")
class TestChiMixHandReference:
    """The production weighting matches an INDEPENDENT hand calc to FP.

    Discriminator: χ_0 ≠ χ_1 and prod_0 ≠ prod_1 so the answer differs
    from BOTH the legacy shortcut (χ_0) and the unweighted mean.
    """

    def test_two_fissile_matches_hand_weighted_average(self):
        # Controlled inputs (single nonzero fission group each).
        aden = np.array([2.0, 1.0])  # iso0 denser
        nubar0, sigf0, g0 = 2.5, 0.10, 0  # prod0 = 2.0 * 2.5 * 0.10 = 0.50
        nubar1, sigf1, g1 = 3.0, 0.20, 5  # prod1 = 1.0 * 3.0 * 0.20 = 0.60
        chi0 = _chi_simplex((0, 0.8), (1, 0.2))
        chi1 = _chi_simplex((5, 0.3), (6, 0.7))
        iso0 = _fissile_isotope(nubar_val=nubar0, sigf_val=sigf0, sigf_group=g0, chi=chi0)
        iso1 = _fissile_isotope(nubar_val=nubar1, sigf_val=sigf1, sigf_group=g1, chi=chi1)
        mix = compute_macro_xs([iso0, iso1], aden, n_legendre=1, fissile_indices=[0, 1])

        # INDEPENDENT reference — explicit scalars laid out term-by-term, NOT a
        # loop re-spelling the code's sum(w_i * χ_i). Shares only numpy
        # primitives with the SUT (below the trusted-library line).
        p0 = aden[0] * nubar0 * sigf0  # 0.50  (single nonzero group)
        p1 = aden[1] * nubar1 * sigf1  # 0.60
        expected = (p0 * chi0 + p1 * chi1) / (p0 + p1)

        np.testing.assert_allclose(np.asarray(mix.chi), expected, atol=1e-12, rtol=0)

        # Discriminator guards (prove the test catches the WRONG formula).
        _require(
            not np.allclose(np.asarray(mix.chi), chi0),
            "χ_mix == χ_0 → the legacy first-isotope shortcut survived",
        )
        _require(
            not np.allclose(np.asarray(mix.chi), (chi0 + chi1) / 2),
            "χ_mix == unweighted mean → production weighting not applied",
        )


# ══════════════════════════════════════════════════════════════════════
# Gate 3 — single-fissile byte-identity
# ══════════════════════════════════════════════════════════════════════


class TestSingleFissileByteIdentity:
    """χ_mix collapses to identity (w=[1]) for a single fissile isotope.

    Proves S10b is scoped to multi-fissile and leaves every single-fissile
    case bit-identical to the legacy shortcut.
    """

    def test_single_fissile_chi_equals_isotope_chi(self):
        # 3a — single fissile isotope: the weighted average collapses to
        # identity (w=[1]) → χ_mix == that isotope.chi, byte-identical.
        chi0 = _chi_simplex((0, 0.6), (1, 0.4))
        iso0 = _fissile_isotope(nubar_val=2.4, sigf_val=0.10, sigf_group=0, chi=chi0)
        mix = compute_macro_xs(
            [iso0, _inert_isotope()], np.array([1.0, 1.0]),
            n_legendre=1, fissile_indices=[0],
        )
        # np.asarray on RHS: iso0.chi is an EmissionSpectrum (S10a); compare values.
        np.testing.assert_array_equal(np.asarray(mix.chi), np.asarray(iso0.chi))

    def test_no_fissile_chi_is_zeros(self):
        # 3b — the fissile_indices=[] / no-fissile default branch
        # (borated_water / zircaloy_clad shape) is unchanged by S10b.
        mix = compute_macro_xs(
            [_inert_isotope(), _inert_isotope()], np.array([1.0, 1.0]),
            n_legendre=1, fissile_indices=[],
        )
        np.testing.assert_array_equal(np.asarray(mix.chi), np.zeros(NG))


# ══════════════════════════════════════════════════════════════════════
# Gate 5 / optional — real-GENDF smoke check (U_235.h5 + U_238.h5 present)
# ══════════════════════════════════════════════════════════════════════

_MICRO_XS = Path(__file__).resolve().parents[2] / "orpheus" / "data" / "micro_xs"
_HAS_U235_U238 = (_MICRO_XS / "U_235.h5").exists() and (_MICRO_XS / "U_238.h5").exists()


@pytest.mark.skipif(
    not _HAS_U235_U238,
    reason="real U_235.h5 + U_238.h5 macro-XS library absent",
)
class TestRealUO2Smoke:
    """Smoke: real 2-fissile UO2 χ_mix is a simplex AND moved off both inputs.

    NOT an L11 reference (gate 2's synthetic hand-calc is). Pins the
    STRUCTURAL facts on real production data: (a) the S10a interlock holds
    on the real weighted average; (b) the multi-fissile weighting actually
    FIRED (χ_mix ≠ χ_U235 AND ≠ χ_U238) — the U238 contribution is real,
    NOT the dropped-isotope shortcut. Does NOT pin the exact χ_mix value
    (that would be a brittle data-version pin).
    """

    @pytest.mark.slow
    def test_uo2_chi_mix_simplex_and_moved(self):
        from orpheus.data.micro_xs import load_isotope

        U235 = load_isotope("U_235", 900)
        U238 = load_isotope("U_238", 900)
        O16 = load_isotope("O_016", 900)
        # Direct compute_macro_xs (skip the recipe's pyXSteam/density ceremony):
        densities = np.array([0.03, 0.97, 2.0])  # representative UO2 ratios
        mix = compute_macro_xs([U235, U238, O16], densities, fissile_indices=[0, 1])
        # (a) the S10a interlock holds on the real weighted average:
        EmissionSpectrum(np.asarray(mix.chi)).assert_normalized()  # MUST NOT raise
        # (b) the multi-fissile weighting actually FIRED — χ_mix moved off
        # BOTH inputs (U238 contribution is real, NOT the dropped-isotope
        # shortcut). Structural fact only — NOT a brittle exact-value pin.
        _require(
            not np.allclose(np.asarray(mix.chi), np.asarray(U235.chi)),
            "χ_mix == χ_U235 → multi-fissile weighting did not fire",
        )
        _require(
            not np.allclose(np.asarray(mix.chi), np.asarray(U238.chi)),
            "χ_mix == χ_U238 → weighting collapsed to a single isotope",
        )
