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
from orpheus.data.micro_xs.gendf import (
    _extract_mf3,
    _extract_mf6,
    _GXS_DIR,
    _parse_gendf,
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


def _mf6_l0_rowsum(m: np.ndarray) -> np.ndarray:
    """Raw MF=6/MT=16 ell=0 row sums, straight off the tape."""
    from scipy.sparse import csr_matrix

    ifrom, ito, sig = _extract_mf6(16, 1, m)
    return np.asarray(
        csr_matrix(
            (sig[(0, 0)], (ifrom - 1, ito - 1)), shape=(NG, NG)
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

    def test_isotope_sig2_rowsum_is_the_reaction_xs(self, be9_tape) -> None:
        iso = convert_gxs(_ISOTOPE)[0]
        # convert_gxs flips to the canonical fast-first order; so must the
        # reference (_to_canonical_group_order reverses both axes).
        sigma = np.asarray(_extract_mf3(16, 1, be9_tape))[0][::-1]
        rows = np.asarray(iso.sig2.sum(axis=1)).ravel()
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
    def be9_mixture(self):
        return compute_macro_xs(convert_gxs(_ISOTOPE)[0:1], np.array([0.1235]))

    def test_absorption_counts_the_collision_once(
        self, be9_mixture, be9_tape,
    ) -> None:
        sigma = np.asarray(_extract_mf3(16, 1, be9_tape))[0][::-1] * 0.1235
        contribution = np.asarray(be9_mixture.Sig2.sum(axis=1)).ravel()
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
        removal = np.asarray(be9_mixture.Sig2.sum(axis=1)).ravel()
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
            self._transfer(np.array([2.0, 4.0, 8.0])), sigma, mt=16,
        )
        np.testing.assert_allclose(
            np.asarray(out.sum(axis=1)).ravel()[:3], [1.0, 2.0, 4.0],
            rtol=1e-14, atol=0.0,
        )

    def test_an_all_zero_transfer_is_returned_untouched(self) -> None:
        empty = self._transfer(np.zeros(3))
        assert _strip_transfer_yield(empty, None, mt=16) is empty

    def test_a_missing_mf3_partner_is_refused(self) -> None:
        with pytest.raises(ValueError, match="MF=3 cross section is absent"):
            _strip_transfer_yield(self._transfer(np.array([2.0])), None, mt=16)

    def test_a_vanishing_cross_section_is_refused(self) -> None:
        with pytest.raises(ValueError, match="MF=3 cross section is 0"):
            _strip_transfer_yield(
                self._transfer(np.array([2.0])), np.zeros((1, NG)), mt=16,
            )

    def test_a_non_integral_yield_is_refused(self) -> None:
        sigma = np.zeros((1, NG)); sigma[0, 0] = 2.0
        with pytest.raises(ValueError, match="not a positive integer yield"):
            _strip_transfer_yield(
                self._transfer(np.array([3.0])), sigma, mt=16,
            )
