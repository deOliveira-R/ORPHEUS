r"""The :math:`(n,2n)` yield convention at the GENDF ingest boundary.

GENDF's MF=6 records hold :math:`\sigma(E)\,y(E)\,f(E\to E')` — the transfer
matrix carries the reaction's **yield** — while MF=3 holds the un-multiplied
group cross section. For MT=16 that yield is 2, so the raw MF=6 row sum is
:math:`2\Sigma_{2n}`.

Everything downstream of :class:`~orpheus.data.micro_xs.isotope.Isotope`
reads ``sig2`` as the **reaction** matrix with no multiplicity folded in:
:attr:`~orpheus.data.macro_xs.mixture.Mixture.SigT` and
:attr:`~orpheus.data.macro_xs.mixture.Mixture.absorption_xs` add its row sum
ONCE (one neutron absorbed per event), and
:meth:`~orpheus.transport.kernels.N2NKernel.emission_matrix` applies the
factor itself (two emitted). Until #427 the ingest handed those consumers
``2*Sigma_2n``, so removal was doubled and the emission quadrupled — a
defect masked because every synthetic library mixture has ``Sig2 = 0`` and
because ``compute_macro_xs`` DERIVES ``SigT`` from the same balance identity
``assert_balanced`` then checks.

The three legs below are deliberately different in kind:

* :class:`TestTheTapeCarriesTheYield` reads the **files**, not our code. It
  is the reason the ingest divides anything out, and it is the leg that
  would change if the library were ever re-processed by GROUPR.
* :class:`TestIngestStripsTheYield` pins the **contract** the consumer set
  relies on. It is near-tautological by construction after the fix — that
  is stated rather than hidden — and it exists so a future edit to the
  ingest cannot silently re-introduce the multiplicity.
* :class:`TestEmissionIsTwiceTheReaction` is the **end-to-end** leg, and the
  one with teeth: it reads 4x under the pre-#427 ingest.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.data.macro_xs.mixture import compute_macro_xs
from scipy.sparse import coo_matrix, csr_matrix

from orpheus.data.micro_xs.gendf import (
    _extract_mf3,
    _extract_mf6,
    _GXS_DIR,
    _parse_gendf,
    _reverse_groups_2d,
    _strip_transfer_yield,
    convert_gxs,
)
from orpheus.data.micro_xs.isotope import NG
from orpheus.transport.kernels import N2NKernel

# The tape-fact leg is an L0 term verification against the primary data
# (it reads the GENDF records, not our code) and carries the equation label;
# everything else is a software invariant. Marked per class rather than at
# module scope, because `foundation` may never carry `verifies(...)`.

#: Be-9 is the discriminating nuclide: one of the most pronounced (n,2n)
#: channels in the library, 50 open incident groups, and the material the
#: #426/#427 investigation was run against.
_ISOTOPE = "BE009"


@pytest.fixture(scope="module")
def be9_tape() -> np.ndarray:
    """The parsed BE009 GENDF card array (~2 s; shared across the module)."""
    return _parse_gendf(_GXS_DIR / f"{_ISOTOPE}.GXS")


@pytest.fixture(scope="module")
def be9_isotope():
    """The ingested BE009 at 294 K, straight from the tape (not the h5 store).

    `[M]` ``convert_gxs("BE009")`` builds all four temperatures in 17.8 s; the
    module pays it ONCE (it was paid twice before the #426 rows joined).
    """
    return convert_gxs(_ISOTOPE)[0]


def _mf6_l0_rowsum(m: np.ndarray) -> np.ndarray:
    """Raw MF=6/MT=16 ell=0 row sums, straight off the tape."""
    from scipy.sparse import csr_matrix

    section = _extract_mf6(16, 1, m)
    assert section is not None, "BE009 carries a MF=6/MT=16 section"
    return np.asarray(
        csr_matrix(
            (section.moments[(0, 0)], (section.ifrom - 1, section.ito - 1)),
            shape=(NG, NG),
        ).sum(axis=1)
    ).ravel()


@pytest.mark.l0
class TestTheTapeCarriesTheYield:
    """A fact about the DATA FILES, independent of any ORPHEUS code path."""

    @pytest.mark.verifies("gendf-mf6-yield")
    def test_mf6_rowsum_is_twice_the_mf3_cross_section(self, be9_tape) -> None:
        r"""``rowsum(MF=6/MT=16) == 2 * sigma(MF=3/MT=16)``.

        This is why the ingest renormalises. If a future re-processing of
        the library stopped folding the yield into MF=6, THIS leg reddens
        first and the renormalisation becomes a no-op that must be removed.
        """
        rows = _mf6_l0_rowsum(be9_tape)
        sigma = np.asarray(_extract_mf3(16, 1, be9_tape))[0]
        live = (rows > 0.0) & (sigma > 0.0)
        assert live.sum() == 50, "Be-9's (n,2n) opens in 50 of 421 groups"
        aggregate = rows[live].sum() / sigma[live].sum()
        np.testing.assert_allclose(
            aggregate, 2.0, rtol=1e-6, atol=0.0,
            err_msg="GENDF MF=6/MT=16 must carry the (n,2n) yield of 2.",
        )


@pytest.mark.foundation
class TestIngestStripsTheYield:
    """The CONTRACT every ``sig2`` consumer relies on.

    Near-tautological by construction — ``_strip_transfer_yield`` scales the
    rows onto exactly this reference — and kept for that reason: it names
    the invariant so an ingest edit cannot quietly break it. The leg with
    real teeth is :class:`TestEmissionIsTwiceTheReaction`.
    """

    def test_isotope_sig2_rowsum_is_the_reaction_xs(self, be9_tape, be9_isotope) -> None:
        iso = be9_isotope
        # convert_gxs flips to the canonical fast-first order; so must the
        # reference (_to_canonical_group_order reverses both axes).
        sigma = np.asarray(_extract_mf3(16, 1, be9_tape))[0][::-1]
        rows = np.asarray(iso.sig2[0].sum(axis=1)).ravel()
        live = rows > 0.0
        np.testing.assert_allclose(
            rows[live], sigma[live], rtol=1e-12, atol=0.0,
            err_msg="Isotope.sig2 must be the REACTION matrix (yield removed).",
        )


@pytest.mark.foundation
class TestEmissionIsTwiceTheReaction:
    """End-to-end through the macroscopic path — the leg with teeth.

    Under the pre-#427 ingest the emission column sum reads ``4 * Sigma_2n``
    and the absorption contribution ``2 * Sigma_2n``.
    """

    @pytest.fixture(scope="class")
    def be9_mixture(self, be9_isotope):
        return compute_macro_xs([be9_isotope], np.array([0.1235]))

    def test_absorption_counts_the_collision_once(
        self, be9_mixture, be9_tape,
    ) -> None:
        sigma = np.asarray(_extract_mf3(16, 1, be9_tape))[0][::-1] * 0.1235
        contribution = np.asarray(be9_mixture.Sig2[0].sum(axis=1)).ravel()
        live = contribution > 0.0
        np.testing.assert_allclose(
            contribution[live], sigma[live], rtol=1e-12, atol=0.0,
            err_msg="absorption_xs/SigT must count the (n,2n) collision ONCE.",
        )

    def test_emission_carries_exactly_the_multiplicity(
        self, be9_mixture, be9_tape,
    ) -> None:
        sigma = np.asarray(_extract_mf3(16, 1, be9_tape))[0][::-1] * 0.1235
        emission = N2NKernel.from_mixture(be9_mixture).emission_matrix()
        live = sigma > 0.0
        np.testing.assert_allclose(
            emission[:, live].sum(axis=0),
            N2NKernel.multiplicity * sigma[live],
            rtol=1e-12, atol=0.0,
            err_msg="the (n,2n) emission must be exactly 2*Sigma_2n.",
        )

    def test_net_production_is_one_neutron_per_reaction(
        self, be9_mixture, be9_tape,
    ) -> None:
        r"""emission - removal == :math:`\Sigma_{2n}`: two out, one in."""
        sigma = np.asarray(_extract_mf3(16, 1, be9_tape))[0][::-1] * 0.1235
        emission = N2NKernel.from_mixture(be9_mixture).emission_matrix()
        removal = np.asarray(be9_mixture.Sig2[0].sum(axis=1)).ravel()
        live = sigma > 0.0
        np.testing.assert_allclose(
            emission[:, live].sum(axis=0) - removal[live], sigma[live],
            rtol=1e-12, atol=0.0,
            err_msg="net (n,2n) production must be +1 neutron per reaction.",
        )


@pytest.mark.foundation
class TestStripTransferYieldGuards:
    """The ingest guard's own witnesses — one positive, three negative.

    (vv-principles #11: a validation helper tested only against broken input
    validates the raising, never the invariant. Each negative leg is mutated
    separately — #17's granularity clause — so an arm with no witness is
    visible rather than certified by its neighbour.)
    """

    @staticmethod
    def _transfer(rows: np.ndarray):
        from scipy.sparse import csr_matrix

        m = np.zeros((NG, NG))
        m[: len(rows), 0] = rows
        return csr_matrix(m)

    def test_positive_yield_two_is_divided_out(self) -> None:
        sigma = np.zeros((1, NG)); sigma[0, :3] = [1.0, 2.0, 4.0]
        out = _strip_transfer_yield(
            [self._transfer(np.array([2.0, 4.0, 8.0]))], sigma, mt=16,
        )
        np.testing.assert_allclose(
            np.asarray(out[0].sum(axis=1)).ravel()[:3], [1.0, 2.0, 4.0],
            rtol=1e-14, atol=0.0,
        )

    def test_the_whole_stack_is_scaled_by_the_p0_diagonal(self) -> None:
        r"""One yield per incident group multiplies EVERY order (#426).

        ``σ_ℓ = σ·y·f_ℓ`` on the tape (ENDF-102 Eq. (6.1)/(6.3), NJOY Eq.
        (242)); the per-row scale is derived from ℓ = 0 and applied to ℓ ≥ 1,
        whose own row sums (``y·σ·⟨P_ℓ⟩``) could not be normalised. So a
        stack ``[2σ f, 2σ f·μ̄]`` must come back ``[σ f, σ f·μ̄]`` — the ℓ=1 /
        ℓ=0 RATIO is untouched (the two-sided law the entrywise bound
        ``|Σ_ℓ| ≤ Σ_0`` cannot check).
        """
        sigma = np.zeros((1, NG)); sigma[0, :3] = [1.0, 2.0, 4.0]
        p0 = self._transfer(np.array([2.0, 4.0, 8.0]))
        p1 = self._transfer(np.array([0.6, -1.2, 2.4]))    # 2σ f · μ̄, μ̄ = ±0.3
        out = _strip_transfer_yield([p0, p1], sigma, mt=16)
        assert len(out) == 2
        np.testing.assert_allclose(
            np.asarray(out[1].sum(axis=1)).ravel()[:3], [0.3, -0.6, 1.2],
            rtol=1e-14, atol=0.0,
        )
        dense0, dense1 = (np.asarray(m.todense()) for m in out)
        np.testing.assert_array_equal(
            dense1[:3, 0] / dense0[:3, 0], [0.3, -0.3, 0.3],
            err_msg="the ℓ=1/ℓ=0 ratio must survive the strip exactly",
        )

    def test_an_all_zero_transfer_is_returned_untouched(self) -> None:
        empty = self._transfer(np.zeros(3))
        assert _strip_transfer_yield([empty], None, mt=16)[0] is empty

    def test_a_missing_mf3_partner_is_refused(self) -> None:
        with pytest.raises(ValueError, match="MF=3 cross section is absent"):
            _strip_transfer_yield([self._transfer(np.array([2.0]))], None, mt=16)

    def test_a_vanishing_cross_section_is_refused(self) -> None:
        with pytest.raises(ValueError, match="MF=3 cross section is 0"):
            _strip_transfer_yield(
                [self._transfer(np.array([2.0]))], np.zeros((1, NG)), mt=16,
            )

    def test_a_non_integral_yield_is_refused(self) -> None:
        sigma = np.zeros((1, NG)); sigma[0, 0] = 2.0
        with pytest.raises(ValueError, match="not a positive integer yield"):
            _strip_transfer_yield(
                [self._transfer(np.array([3.0]))], sigma, mt=16,
            )


# ── #426 step 1: the ingest is lossless in ℓ, and the stack obeys the format's physics ──

_TEMP_IDX = 1  # 293.6 K on the tape ↔ load_isotope(name, 294)
_L_TAPE = 6    # NL = 7 on MT=16, `[M]` 10 of the 11 tapes carrying it
_BOUND_ROWS = [("BE009", 16), ("BE009", 2), ("U_235", 16), ("U_235", 2)]


@pytest.fixture(scope="module")
def tapes(be9_tape) -> dict[str, np.ndarray]:
    """BE009 (shared with the module) + U_235 (`[M]` 5.4 s), parsed ONCE."""
    return {"BE009": be9_tape, "U_235": _parse_gendf(_GXS_DIR / "U_235.GXS")}


@pytest.mark.l0
class TestTapePhysics:
    """REFERENCE-class: theorems of the ENDF-102 / GROUPR format, independent of every
    line the ingest carve touches. `[M]` 2026-09-03 max|Σ_ℓ|/Σ_0 over live ℓ=0 entries,
    ℓ = 1…6: BE009/MT16 0.9603→0.3190, BE009/MT2 0.9997→0.9980, U_235/MT16
    0.8977→0.5552, U_235/MT2 0.9773→0.6009 — 0 violations of 8195 / 5080 / 6067 / 895.

    ⚠ The bound is ONE-SIDED: a stray ``(2ℓ+1)`` on ℓ=1 reads ≈2.9 and a stray ×2 ≈1.9
    (both caught); a stray ``1/(2ℓ+1)`` reads 0.32 and passes. The two-sided catcher is
    :meth:`TestThreading.test_the_yield_strip_is_one_diagonal`.
    ⛔ No "moments decay monotonically in ℓ" leg: true of MT=16, FALSE of MT=221
    (`[M]` BE009 thermal ℓ=6 > ℓ=5) — a property of the channel, not of the ingest.
    """

    @pytest.mark.parametrize("name,mt", _BOUND_ROWS, ids=lambda v: str(v))
    def test_every_legendre_coefficient_is_bounded_by_its_p0(self, tapes, name, mt):
        r"""``|Σ_ℓ(g'→g)| ≤ Σ_0(g'→g)`` entrywise: ``Σ_ℓ/Σ_0 = ⟨P_ℓ(μ)⟩`` and ``|P_ℓ| ≤ 1``.

        Refuses an addition-theorem factor entering the STORED moment — the
        ``(2ℓ+1)`` belongs on ``LegendreBasis.reconstruct``, once. Threshold
        ``1 + 1e-9``, not tighter: `[M]` BE009's elastic ℓ=1 margin is 3e-4.
        """
        got = _extract_mf6(mt, _TEMP_IDX, tapes[name])
        assert got is not None, f"{name} carries no MF=6/MT={mt} section"
        sig = got.moments
        v0 = np.asarray(sig[(0, 0)])
        live = np.abs(v0) > 0.0
        assert live.sum() > 0, "no live ℓ=0 entries — the bound is vacuous"
        for l in range(1, _L_TAPE + 1):
            if (l, 0) not in sig:
                continue
            ratio = np.abs(np.asarray(sig[(l, 0)])[live]) / np.abs(v0[live])
            n_over = int((ratio > 1.0 + 1e-9).sum())
            if n_over:
                pytest.fail(
                    f"{name} MT={mt}: {n_over} of {int(live.sum())} entries have "
                    f"|Σ_{l}|/Σ_0 > 1 (max {ratio.max():.4f}) — a stored Legendre "
                    f"coefficient cannot exceed its P0; a stray (2ℓ+1) would read "
                    f"≈{2 * l + 1:.1f}"
                )

    @pytest.mark.parametrize("name", ["BE009", "U_235"])
    def test_the_n2n_transfer_only_loses_energy_at_every_order(self, tapes, name):
        r"""Strictly upper-triangular in canonical fast-first order, ℓ = 0…6.

        A threshold reaction: neither emitted neutron can exceed the incident
        energy. The ``(ifrom, ito)`` arrays are SHARED across ℓ, so every ℓ
        being upper-triangular pins that the per-ℓ reversal used the same
        arrays — a row/column swap puts the whole mass below the diagonal.
        """
        got = _extract_mf6(16, _TEMP_IDX, tapes[name])
        assert got is not None
        ifrom, ito, sig = got.ifrom, got.ito, got.moments
        for l in range(_L_TAPE + 1):
            if (l, 0) not in sig:
                continue
            rev = _reverse_groups_2d(
                csr_matrix(coo_matrix(
                    (np.asarray(sig[(l, 0)]), (ifrom - 1, ito - 1)), shape=(NG, NG),
                ).tocsr())
            ).tocoo()
            lower = int((rev.col < rev.row).sum())
            assert lower == 0, (
                f"{name} ℓ={l}: {lower} entries below the diagonal in canonical "
                f"fast-first order — the (n,2n) channel would be up-scattering"
            )


def _reference_stack(m: np.ndarray) -> list[csr_matrix]:
    """The ℓ = 0…6 (n,2n) stack built the way the ingest must build it.

    ⚠ SHARED with production (``_extract_mf6``, ``_strip_transfer_yield``,
    ``_reverse_groups_2d``): this pins the THREADING — which ℓ lands in which
    slot, both axes reversed per ℓ, the strip applied to the whole stack — not
    the parse. The parse is pinned by :class:`TestTapePhysics` (physics) and
    by :class:`TestTheTapeCarriesTheYield` (the yield identity).
    """
    section = _extract_mf6(16, _TEMP_IDX, m)
    assert section is not None, "BE009 carries no MF=6/MT=16 section"
    ifrom, ito, sig = section.ifrom, section.ito, section.moments
    raw = [
        csr_matrix(coo_matrix(
            (np.asarray(sig[(l, 0)]), (ifrom - 1, ito - 1)), shape=(NG, NG),
        ).tocsr())
        for l in range(_L_TAPE + 1)
    ]
    stripped = _strip_transfer_yield(raw, _extract_mf3(16, _TEMP_IDX, m), mt=16)
    return [_reverse_groups_2d(csr_matrix(s)) for s in stripped]


@pytest.fixture(scope="module")
def be9_reference(be9_tape) -> list[csr_matrix]:
    return _reference_stack(be9_tape)


@pytest.mark.foundation
class TestThreading:
    """THEOREM-class: the code contract — ``Isotope.sig2`` IS the tape's stack.

    §6c red-before, `[M]` 2026-09-03: pre-carve ``Isotope.sig2`` was one
    ``csr_matrix`` holding ℓ = 0, so every row here reddened.
    """

    @pytest.mark.parametrize("l", range(_L_TAPE + 1))
    def test_isotope_sig2_carries_the_tape_moment_at_every_order(
        self, be9_isotope, be9_reference, l,
    ):
        """``array_equal`` and nothing looser: a gather plus ONE diagonal multiply."""
        stack = be9_isotope.sig2
        assert len(stack) == _L_TAPE + 1, (
            f"the tape carries NL = {_L_TAPE + 1} for MT=16 on BE009; the Isotope kept {len(stack)}"
        )
        np.testing.assert_array_equal(
            np.asarray(stack[l].todense()), np.asarray(be9_reference[l].todense()),
            err_msg=f"Isotope.sig2[{l}] is not the tape's ℓ={l} (n,2n) moment",
        )

    @pytest.mark.parametrize("l", range(1, _L_TAPE + 1))
    def test_the_yield_strip_is_one_diagonal(self, be9_isotope, be9_tape, l):
        r"""The stored ℓ / ℓ=0 ratio equals the RAW tape ratio (to an ulp).

        A row-diagonal cancels in ``Σ_ℓ/Σ_0``; a per-ℓ scale (``scale**ℓ``,
        ``scale/(2ℓ+1)``, a wrong diagonal) moves it and reds HERE while the
        entrywise bound stays satisfied — the two-sided law.
        """
        stack = be9_isotope.sig2
        section = _extract_mf6(16, _TEMP_IDX, be9_tape)
        assert section is not None
        sig = section.moments
        raw0 = np.asarray(sig[(0, 0)])
        rawl = np.asarray(sig[(l, 0)])
        live = np.abs(raw0) > 0.0
        want = rawl[live] / raw0[live]
        # _reverse_groups_2d flips BOTH axes: 1-based (ifrom, ito) → (NG - ifrom, NG - ito)
        rev_rows = NG - section.ifrom[live]
        rev_cols = NG - section.ito[live]
        s0 = np.asarray(stack[0].todense())[rev_rows, rev_cols]
        sl = np.asarray(stack[l].todense())[rev_rows, rev_cols]
        # rtol 1e-14, not array_equal: (s·a)/(s·b) and a/b differ by an ulp
        # (`[M]` max 2.6e-16 relative); a per-ℓ scale moves the ratio by O(1).
        np.testing.assert_allclose(
            sl / s0, want, rtol=1e-14, atol=0.0,
            err_msg=(
                f"the ℓ={l} moment was scaled by a different diagonal than ℓ=0 — the yield "
                f"is one number per incident group and multiplies the whole emission distribution"
            ),
        )
