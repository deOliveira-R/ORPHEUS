r"""The transfer family's LAWS at the datum tier (#426 step 2, the F2 ruling).

One kernel type carries both collision-gain channels —
:class:`~orpheus.transport.kernels.TransferKernel` ``(moments, multiplicity)``
— and the two differ by the yield :math:`y` alone. These rows pin what makes
that a law rather than a restatement of ``2``:

* **G2.2 — admission.** The yield is a positive integer (1 for scattering,
  :data:`~orpheus.transport.kernels.N2N_MULTIPLICITY` for (n,2n)); ``0``, a
  negative, a float and a bool are refused; the default is the identity.
* **G2.2 — the emission law** ``emission_matrix() == y · p0.T`` for BOTH
  ``y ∈ {1, 2}`` — the ``y = 1`` leg is what makes it a law.
* **§4.3 — ``at_order``** is the identity at the stored order (the SAME
  object), truncates below it, and PADS exact zeros above it (ruled: a
  shorter stack is complete — the evaluation's zeros, not an invention).
* **``is_isotropic``** is a structural fact of the datum: exact zeros above
  ℓ = 0 (an absent section, ``NL = 1``, a padded stack) — the predicate the
  angular binding uses to skip an :math:`R\Lambda M` product that would
  return exact zeros.
* **G2.3 — the field is ONE channel**: mixed yields are refused; the
  group-rate verb (the k-balance's accounting, ``y ∫ Σ_{c,0}ᵀ φ dV``)
  vanishes for NO channel and scales by ``y`` against a hand loop — the
  ``y = 1`` row is the in-scatter rate the code never spelled because it
  never needed to.

⚠ These are DATUM/VERB claims. The angular-binding claims (the conjugation
at the solve's order, the Be-reflected eigenvalue) live in
``tests/sn/operators/test_n2n_operator.py`` and
``tests/sn/verification/analytical/test_be_reflected_n2n_anisotropy.py``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.transport.kernels import N2N_MULTIPLICITY, TransferKernel
from orpheus.transport.material_field import TransferMaterialField

pytestmark = pytest.mark.foundation

_P0 = np.array([[0.38, 0.10], [0.05, 0.90]])
_P1 = np.array([[0.20, -0.04], [0.02, 0.30]])
_P2 = np.array([[0.06, 0.01], [-0.01, 0.09]])


def _stack(multiplicity: int = 1, *orders: np.ndarray) -> TransferKernel:
    return TransferKernel(moments=tuple(orders) or (_P0, _P1, _P2), multiplicity=multiplicity)


class TestKernelAdmission:
    def test_the_default_yield_is_the_identity(self):
        assert TransferKernel(moments=(_P0,)).multiplicity == 1

    def test_the_n2n_constant_has_one_home(self):
        assert N2N_MULTIPLICITY == 2

    @pytest.mark.parametrize("bad", [0, -1, 2.0, True], ids=["zero", "negative", "float", "bool"])
    def test_a_non_positive_or_non_integer_yield_is_refused(self, bad):
        with pytest.raises(ValueError, match="positive integer"):
            TransferKernel(moments=(_P0,), multiplicity=bad)

    def test_the_yield_survives_re_ordering(self):
        k = _stack(N2N_MULTIPLICITY)
        assert k.at_order(0).multiplicity == N2N_MULTIPLICITY
        assert k.at_order(5).multiplicity == N2N_MULTIPLICITY


class TestTheEmissionLaw:
    @pytest.mark.parametrize("y", [1, N2N_MULTIPLICITY], ids=["scattering", "n2n"])
    def test_emission_matrix_is_the_yield_times_the_p0_transpose(self, y):
        k = _stack(y)
        np.testing.assert_array_equal(k.emission_matrix(), y * _P0.T)

    def test_emission_matrix_is_a_fresh_writable_copy(self):
        k = _stack(N2N_MULTIPLICITY)
        m = k.emission_matrix()
        assert m.flags.writeable and m is not k.p0
        m[0, 0] = 99.0
        np.testing.assert_array_equal(k.p0, _P0)  # the stored datum is untouched


class TestAtOrder:
    def test_identity_at_the_stored_order_is_the_same_object(self):
        k = _stack()
        assert k.at_order(k.order) is k

    def test_below_the_stored_order_truncates(self):
        k = _stack()
        p1 = k.at_order(1)
        assert p1.order == 1
        np.testing.assert_array_equal(p1.moments[0], _P0)
        np.testing.assert_array_equal(p1.moments[1], _P1)

    def test_above_the_stored_order_pads_exact_zeros(self):
        k = _stack(N2N_MULTIPLICITY, _P0)
        wide = k.at_order(3)
        assert wide.order == 3 and wide.multiplicity == N2N_MULTIPLICITY
        np.testing.assert_array_equal(wide.moments[0], _P0)
        for l in (1, 2, 3):
            np.testing.assert_array_equal(wide.moments[l], 0.0)
            assert not wide.moments[l].flags.writeable
        # the pads are distinct buffers, not one shared array
        assert wide.moments[1] is not wide.moments[2]

    def test_a_negative_order_is_refused(self):
        with pytest.raises(ValueError, match="order >= 0"):
            _stack().at_order(-1)


class TestIsIsotropic:
    def test_a_length_one_stack_is_isotropic(self):
        assert TransferKernel(moments=(_P0,)).is_isotropic

    def test_exact_zeros_above_l0_are_isotropic(self):
        assert _stack(1, _P0, np.zeros((2, 2)), np.zeros((2, 2))).is_isotropic
        assert TransferKernel(moments=(_P0,)).at_order(6).is_isotropic

    def test_one_non_zero_entry_above_l0_is_anisotropic(self):
        l1 = np.zeros((2, 2))
        l1[1, 0] = 1e-300
        assert not _stack(1, _P0, l1).is_isotropic

    def test_the_field_is_isotropic_iff_every_material_is(self):
        layout = {0: (np.array([0]),), 1: (np.array([1]),)}
        iso = TransferMaterialField(
            per_material={0: TransferKernel(moments=(_P0,)), 1: TransferKernel(moments=(_P0,))},
            cells_by_material=layout,
        )
        assert iso.is_isotropic
        mixed = TransferMaterialField(
            per_material={0: TransferKernel(moments=(_P0,)), 1: _stack(1, _P0, _P1)},
            cells_by_material=layout,
        )
        assert not mixed.is_isotropic


class TestTheFieldIsOneChannel:
    _LAYOUT = {0: (np.array([0, 1]),), 1: (np.array([2]),)}

    def test_mixed_yields_are_refused(self):
        with pytest.raises(ValueError, match="ONE channel"):
            TransferMaterialField(
                per_material={0: _stack(1), 1: _stack(N2N_MULTIPLICITY)},
                cells_by_material=self._LAYOUT,
            )

    @pytest.mark.parametrize("y", [1, N2N_MULTIPLICITY], ids=["scattering", "n2n"])
    def test_the_group_rate_verb_scales_by_the_yield(self, y):
        r"""``rate_g += y ∫_V Σ_{c,0}(g'→g) φ_{g'} dV`` per material, against a
        hand loop (plain ``@``; no einsum, no shared dispatch — vv L11).
        For y = 1 this is the group in-scatter rate; for y = 2 the (n,2n)
        emission the k-balance subtracts from removal."""
        field = TransferMaterialField(
            per_material={0: _stack(y), 1: _stack(y, _P1, _P2)},  # material 1: a different P0
            cells_by_material=self._LAYOUT,
        )
        rng = np.random.default_rng(4)
        phi = rng.uniform(0.1, 1.0, size=(2, 3))
        volume = np.array([1.0, 1.5, 2.0])
        rate = np.zeros(2)
        field.add_to_group_rate(rate, phi, volume)
        ref = np.zeros(2)
        for ix, mid in ((0, 0), (1, 0), (2, 1)):
            p0 = field.per_material[mid].p0
            ref += y * volume[ix] * (p0.T @ phi[:, ix])
        np.testing.assert_allclose(rate, ref, rtol=1e-15, atol=0)
        assert field.multiplicity == y
