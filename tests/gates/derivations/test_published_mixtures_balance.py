"""The published cross-section tables balance; the Atalay criticality mixtures are the declared exemption.

Moved from ``tests/gates/data/test_mixture_xs_balance.py`` on 2026-10-26 (layer DAG ruling 29): the sweep
over the Sood 2003 registry and the homogeneous derive builders, and the rows pinning the Atalay 1997
criticality-parameter mixtures as intentionally imbalanced, are claims about the derivations' data, not
about ``Mixture``. The balance law itself (``Mixture.assert_balanced``, its positive and negative legs) and
the synthetic-catalogue sweep stay in data, and so does the hand-spelled residual these rows import.

Foundation tests (software invariant — no theory ``:label:``). Every assertion is ``np.testing.*`` /
``pytest.raises`` (the suite runs ``-O``).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.data.macro_xs.mixture import Mixture
from tests.gates.data.test_mixture_xs_balance import _ATOL, _balance_residual, _needs_api

pytestmark = pytest.mark.foundation


def _sood2003_mixtures():
    from orpheus.derivations.continuous.sood_registry import SOOD2003_CASES

    out = []
    for case_id, case in SOOD2003_CASES.items():
        for mat_id, mix in case.materials.items():
            out.append(pytest.param(mix, id=f"sood2003-{case_id}-m{mat_id}"))
    return out


def _homogeneous_mixtures():
    from orpheus.derivations.continuous.analytical.homogeneous import (
        derive_1g,
        derive_2g,
        derive_4g,
    )

    out = []
    for fn in (derive_1g, derive_2g, derive_4g):
        case = fn()
        for mat_id, mix in case.materials.items():
            out.append(pytest.param(mix, id=f"homog-{fn.__name__}-m{mat_id}"))
    return out


@_needs_api
class TestPhysicalTableSweep:
    """Every PUBLISHED Mixture must balance. The typo-catcher."""

    @pytest.mark.parametrize("mix", _sood2003_mixtures())
    def test_sood2003_balances(self, mix: Mixture) -> None:
        mix.assert_balanced(atol=_ATOL)

    @pytest.mark.parametrize("mix", _homogeneous_mixtures())
    def test_homogeneous_balances(self, mix: Mixture) -> None:
        mix.assert_balanced(atol=_ATOL)


class TestExemptionIntegrity:
    """Intentionally-imbalanced published Mixtures still CONSTRUCT (no __post_init__ raise)."""

    def test_atalay_criticality_mixture_constructs(self) -> None:
        """Atalay νΣf=(c-1)Σt c>1 mixture (residual≈SigF) still builds.

        Atalay is in ATALAY1997_CASES, NOT in SOOD2003_CASES, so gate 2 never
        touches it. Confirm it constructs and IS imbalanced (the exemption is
        real, not an accident of a balanced fixture).
        """
        from orpheus.derivations.continuous.sood_registry.atalay1997 import (
            _mix_iso_at_c,
        )

        mix = _mix_iso_at_c(1.30)  # MUST NOT raise on construction
        rowsum_s0 = np.array(mix.SigS[0].sum(axis=1)).ravel()
        rowsum_2 = np.array(mix.Sig2[0].sum(axis=1)).ravel()
        derived = mix.SigC + mix.SigL + mix.SigF + rowsum_s0 + rowsum_2
        residual = float(np.max(np.abs(mix.SigT - derived)))
        # It is INTENTIONALLY imbalanced (the c>1 criticality encoding).
        np.testing.assert_array_less(
            1e-3, residual,
            err_msg=f"Atalay c=1.30 must be (intentionally) imbalanced; residual={residual}",
        )

    def test_atalay_not_in_sood2003_registry(self) -> None:
        """Structural exemption: no Atalay case_id appears in SOOD2003_CASES.

        This is WHY gate 2 needs no defensive filter — the swept registry
        simply does not contain the exempt set. If a future edit adds Atalay
        to SOOD2003_CASES, this reddens BEFORE gate 2 starts failing.
        """
        from orpheus.derivations.continuous.sood_registry import SOOD2003_CASES

        atalay_ids = [cid for cid in SOOD2003_CASES if "atalay" in cid.lower()]
        np.testing.assert_array_equal(
            len(atalay_ids), 0,
            err_msg=f"Atalay must stay out of SOOD2003_CASES; found {atalay_ids}",
        )


class TestAtolBand:
    """The violator side of the atol band: the Atalay encoding stays far above atol=1e-9."""

    def test_atalay_violator_above_1e_3(self) -> None:
        from orpheus.derivations.continuous.sood_registry.atalay1997 import (
            _mix_iso_at_c,
        )

        r = _balance_residual(_mix_iso_at_c(1.30))
        np.testing.assert_array_less(
            1e-3, r,
            err_msg=f"Atalay c=1.30 residual={r} should be ≫ atol (intentional)",
        )
