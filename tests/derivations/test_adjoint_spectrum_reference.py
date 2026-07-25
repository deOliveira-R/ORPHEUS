r"""Pin the adjoint (left-eigenvector) k∞ reference — campaign #276 A4 / P1.4.

:func:`kinf_and_adjoint_spectrum_homogeneous` is the Branch-1 closed-form
reference the SN adjoint flux-shape gate (P1.4) compares against.  This file
pins the reference itself — its DEFINING law and its comparability contract
with the forward sibling — before any SN machinery consumes it:

* the left eigen-law :math:`\varphi^{*T}\mathbf{M} = k\,\varphi^{*T}`
  against an INDEPENDENTLY assembled resolvent (explicit
  ``diag``/``outer``/``solve`` in-test, no call into the module helper);
* ``k_adjoint == k_forward`` exactly (:math:`\text{eig}(\mathbf{M}^T) =
  \text{eig}(\mathbf{M})`) at 2G and 4G;
* the 4G non-degeneracy dud-guard: :math:`\varphi^*` is genuinely ≠
  :math:`\varphi` and non-flat (Mixture A **2G** is coincidentally flat —
  the documented reason P1.4 mandates 4G), with the :math:`\chi \nparallel
  \nu\Sigma_f` precondition asserted so the downstream F†-swap teeth cannot
  silently go blind.

vv Mode-8: ``np.testing.*`` / :func:`require` only (fire under ``python -O``).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.eigenvalue import (
    kinf_and_adjoint_spectrum_homogeneous,
    kinf_and_spectrum_homogeneous,
)
from orpheus.derivations.common.xs_library import get_mixture

pytestmark = pytest.mark.foundation


def require(condition: bool, message: str) -> None:
    if not condition:
        pytest.fail(message)


def _mixture_arrays(ng_key: str):
    """Mixture A → the raw ``(Σ_t, Σ_s0, νΣ_f, χ)`` arrays the references take."""
    mix = get_mixture("A", ng_key)
    return (
        np.asarray(mix.SigT, dtype=float),
        np.asarray(mix.SigS[0].todense(), dtype=float),
        np.asarray(mix.SigP, dtype=float),
        np.asarray(mix.chi, dtype=float),
    )


class TestAdjointSpectrumReference:
    def test_left_eigen_law_4g(self):
        r"""DEFINING law: :math:`\varphi^{*T}\mathbf{M} = k\,\varphi^{*T}`.

        The resolvent is assembled INLINE (explicit ``diag``/``outer``/
        ``solve``) — the test never calls the module's assembly helper, so a
        wrong operator pair inside the reference disagrees with this spelling.
        Also pins the output convention: unit :math:`\ell^2` norm,
        non-negative sum (comparability with the forward spectrum).
        """
        sig_t, sig_s, nu_sig_f, chi = _mixture_arrays("4g")
        k, phi_star = kinf_and_adjoint_spectrum_homogeneous(
            sig_t, sig_s, nu_sig_f, chi,
        )

        M = np.linalg.solve(np.diag(sig_t) - sig_s.T, np.outer(chi, nu_sig_f))
        np.testing.assert_allclose(
            phi_star @ M, k * phi_star, rtol=1e-12, atol=1e-14,
            err_msg="φ*ᵀM ≠ kφ*ᵀ — the returned vector is not the left "
            "eigenvector of the independently assembled resolvent.",
        )
        np.testing.assert_allclose(
            float(np.linalg.norm(phi_star)), 1.0, rtol=1e-12,
            err_msg="adjoint spectrum is not ℓ²-normalised — the "
            "comparability convention with the forward spectrum is broken.",
        )
        require(
            float(phi_star.sum()) >= 0.0,
            "adjoint spectrum violates the non-negative-sum sign convention.",
        )

    @pytest.mark.parametrize("ng_key", ["2g", "4g"])
    def test_k_matches_forward_exactly(self, ng_key):
        r"""``k_adjoint == k_forward`` — :math:`\text{eig}(\mathbf{M}^T) =
        \text{eig}(\mathbf{M})` is an algebraic identity, so the two
        references must agree to LAPACK round-off, not a physics tolerance."""
        arrays = _mixture_arrays(ng_key)
        k_fwd, _ = kinf_and_spectrum_homogeneous(*arrays)
        k_adj, _ = kinf_and_adjoint_spectrum_homogeneous(*arrays)
        np.testing.assert_allclose(
            k_adj, k_fwd, rtol=1e-12,
            err_msg=f"k from eig(Mᵀ) differs from k from eig(M) at {ng_key} — "
            "the two references do not share an operator pair.",
        )

    def test_4g_adjoint_spectrum_differs_from_forward(self):
        r"""Non-degeneracy dud-guard: at 4G, :math:`\varphi^* \ne \varphi`
        and :math:`\varphi^*` is non-flat.

        P1.4 mandates 4G because Mixture A's 2G spectrum is coincidentally
        flat (forward and adjoint coincide — flux-shape-blind).  If THIS
        assert ever fires, the 4G library entry has degenerated and every
        downstream F†-swap tooth that leans on φ*≠φ is silently blind —
        fail here, at the reference, not there.
        """
        sig_t, sig_s, nu_sig_f, chi = _mixture_arrays("4g")

        # Precondition (spec §2 anti-blindness): χ ∦ νΣf, else F†=F up to
        # scale and the role swap is invisible by construction.
        chi_hat = chi / np.linalg.norm(chi)
        nsf_hat = nu_sig_f / np.linalg.norm(nu_sig_f)
        require(
            not np.allclose(chi_hat, nsf_hat, rtol=1e-6),
            "Mixture A 4g has χ ∝ νΣf — the F† role swap is invisible on "
            "this fixture; pick a mixture with non-proportional χ/νΣf.",
        )

        _, phi = kinf_and_spectrum_homogeneous(sig_t, sig_s, nu_sig_f, chi)
        _, phi_star = kinf_and_adjoint_spectrum_homogeneous(
            sig_t, sig_s, nu_sig_f, chi,
        )
        require(
            not np.allclose(phi_star, phi, rtol=1e-3),
            "4G adjoint spectrum ≈ forward spectrum — the fixture is "
            "effectively self-adjoint and P1.4's flux-shape teeth are blind.",
        )
        flat = np.full_like(phi_star, 1.0 / np.sqrt(phi_star.size))
        require(
            not np.allclose(phi_star, flat, rtol=1e-3),
            "4G adjoint spectrum is flat — the energy-adjoint content "
            "P1.4 verifies is absent from this fixture.",
        )
