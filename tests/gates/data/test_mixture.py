"""L0 verification of the macro-sum equation.

Sphinx label: ``math:equation:macro-sum`` (docs/theory/foundations/infinite_medium.rst:622)

The macro-sum equation states

.. math::

    \\Sigma_x(r) = \\sum_i N_i(r) \\, \\sigma_i(x)

i.e. every macroscopic cross section is a linear combination of the
microscopic cross sections of the constituent isotopes weighted by
their number densities.

Implementation: ``orpheus.data.macro_xs.mixture.compute_macro_xs``
lines 121-135, which computes ``SigC = sigC.T @ aDen`` (and similar
for ``SigL``, ``SigF``, ``SigS``, ``Sig2``) after resolving sigma-zero
dependencies through ``interp_xs_field``.

Why this file exists: the full ``compute_macro_xs`` pipeline is
exercised indirectly by every solver eigenvalue test in the suite,
but none of them verifies the underlying formula with a hand
calculation. The QA audit (2026-04-13) flagged ``macro-sum`` as the
only theory equation with zero test coverage even after PR-2 populated
``VerificationCase.equation_labels`` across the derivation package.
This file closes that gap.

The three hand-calculation tests run the production sum: each builds one
synthetic ``Isotope`` per row of its hand table (one sigma-zero base point,
so the interpolation returns the table and the sigma-zero solve is inert)
and calls ``compute_macro_xs``; the hand-computed values stay the
independent oracle. The micro table is written into capture, (n,alpha) and
fission at the scalings 1, 2 and 4 (exact in binary floating point), so a
dropped weight, a transposed sum or a swap of two channels each reads
differently. Until 2026-10-26 they computed ``micro_xs.T @ n_densities``
inside the test and executed no ``orpheus`` line, so their ``verifies``
marker was prose (layer DAG Phase 1b, candidate 9).

Closes issue #84.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture, compute_macro_xs
from orpheus.data.micro_xs.isotope import NG, Isotope

pytestmark = [
    pytest.mark.l0,
    pytest.mark.verifies("macro-sum"),
]


#: The scaling of the hand table per summed channel: distinct powers of two, exact in binary arithmetic.
_CHANNEL_SCALE = {"SigC": 1.0, "SigL": 2.0, "SigF": 4.0}


def _production_macro_sum(n_densities: np.ndarray, micro_xs: np.ndarray) -> dict[str, np.ndarray]:
    """The macroscopic channels ``compute_macro_xs`` sums from the hand table ``micro_xs`` (n_iso, k).

    One ``Isotope`` per row: its capture, (n,alpha) and fission tables are the row scaled by
    ``_CHANNEL_SCALE``, padded with zeros to the library's ``NG`` groups, at a single sigma-zero
    base point (``interp_xs_field`` returns the table itself). No isotope produces (``nubar`` is
    zero, ``fissile_indices`` empty), so the fission spectrum stays null. Returns the first ``k``
    groups of each channel, after checking every padded group sums to exactly zero.
    """
    n_iso, k = micro_xs.shape
    padded = np.zeros((n_iso, NG))
    padded[:, :k] = micro_xs
    grid = np.geomspace(2.0e7, 1.0e-5, NG + 1)
    isotopes = []
    for i in range(n_iso):
        row = {name: (scale * padded[i])[None, :] for name, scale in _CHANNEL_SCALE.items()}
        isotopes.append(Isotope(
            name=f"hand-{i}", aw=1.0, temp=300.0, eg=grid, sig0=np.array([1.0e10]),
            sigC=row["SigC"], sigL=row["SigL"], sigF=row["SigF"],
            sigT=row["SigC"] + row["SigL"] + row["SigF"],
            nubar=np.zeros(NG), chi=np.zeros(NG),
            sigS=[[csr_matrix((NG, NG))]],
        ))
    mix = compute_macro_xs(isotopes, n_densities, escape_xs=1.0, fissile_indices=[])
    channels = {name: np.asarray(getattr(mix, name)) for name in _CHANNEL_SCALE}
    for name, values in channels.items():
        np.testing.assert_array_equal(values[k:], 0.0, err_msg=f"{name}: a group beyond the hand table is not zero")
    return {name: values[:k] for name, values in channels.items()}


def _assert_every_channel(macro: dict[str, np.ndarray], expected: np.ndarray) -> None:
    """Each summed channel equals the hand value at its own scaling."""
    for name, scale in _CHANNEL_SCALE.items():
        np.testing.assert_allclose(macro[name], scale * expected, rtol=1e-14, err_msg=name)


def test_macro_sum_two_isotope_two_group_hand_calc():
    """Σᵢ Nᵢ·σᵢ for an asymmetric 2-isotope, 2-group mixture.

    Hand calc (all in barn, cm⁻¹)::

        N       = [N₁, N₂]       = [5.0,  3.0]
        σ₁      = [σ₁(g₀), σ₁(g₁)] = [0.1, 0.2]   # fast, thermal
        σ₂      = [σ₂(g₀), σ₂(g₁)] = [0.4, 0.6]

        Σ(g₀)   = N₁·σ₁(g₀) + N₂·σ₂(g₀)
                = 5.0·0.1 + 3.0·0.4
                = 0.5 + 1.2
                = 1.7
        Σ(g₁)   = N₁·σ₁(g₁) + N₂·σ₂(g₁)
                = 5.0·0.2 + 3.0·0.6
                = 1.0 + 1.8
                = 2.8

    Catches a transpose flip: ``micro @ N`` instead of ``micro.T @ N``
    cannot even be formed on the (n_iso, NG) table, and a dropped density
    weight reads [0.5, 0.8].
    """
    n_densities = np.array([5.0, 3.0])
    micro_xs = np.array(
        [
            [0.1, 0.2],  # isotope 1
            [0.4, 0.6],  # isotope 2
        ]
    )

    macro = _production_macro_sum(n_densities, micro_xs)

    expected = np.array([1.7, 2.8])
    _assert_every_channel(macro, expected)


def test_macro_sum_three_isotope_three_group_asymmetric():
    """Σᵢ Nᵢ·σᵢ for a 3-isotope, 3-group asymmetric mixture.

    Hand calc per group::

        g₀: 2.0·0.10 + 7.0·0.50 + 1.5·1.00 = 0.20 + 3.50 + 1.50 = 5.20
        g₁: 2.0·0.25 + 7.0·0.10 + 1.5·2.00 = 0.50 + 0.70 + 3.00 = 4.20
        g₂: 2.0·0.05 + 7.0·0.80 + 1.5·0.00 = 0.10 + 5.60 + 0.00 = 5.70

    Two isotopes have zero XS in a particular group — catches bugs
    where a non-zero isotope is dropped by a misplaced if-guard.
    """
    n_densities = np.array([2.0, 7.0, 1.5])
    micro_xs = np.array(
        [
            [0.10, 0.25, 0.05],
            [0.50, 0.10, 0.80],
            [1.00, 2.00, 0.00],
        ]
    )

    macro = _production_macro_sum(n_densities, micro_xs)

    expected = np.array([5.20, 4.20, 5.70])
    _assert_every_channel(macro, expected)


def test_macro_sum_single_isotope_reduces_to_scaling():
    """Σ = N₁·σ₁ for n_iso=1: every macro is just N scaled.

    Degenerate limit. Catches off-by-one in the reduction direction.
    """
    n_densities = np.array([4.0])
    micro_xs = np.array([[0.1, 0.25, 0.5, 1.0]])  # 1 isotope, 4 groups

    macro = _production_macro_sum(n_densities, micro_xs)

    expected = 4.0 * np.array([0.1, 0.25, 0.5, 1.0])
    _assert_every_channel(macro, expected)


def test_macro_sum_mixture_roundtrip_preserves_values():
    """Building a Mixture from hand-specified arrays preserves them bit-exact.

    Not a test of the macro-sum math per se — it's a structural check
    that the dataclass doesn't mangle arrays between construction and
    read. Part of the macro-sum verification because a silent copy or
    type-coercion bug here would still pass the math test above but
    corrupt downstream solver eigenvalues.
    """
    ng = 2
    SigC = np.array([1.7, 2.8])  # carries the result of the first test
    SigL = np.zeros(ng)
    SigF = np.array([0.1, 0.5])
    SigP = 2.5 * SigF  # synthetic nubar = 2.5
    # Synthetic downscatter scattering matrix
    SigS = csr_matrix(np.array([[0.30, 0.10], [0.00, 0.90]]))
    SigT = (
        SigC
        + SigL
        + SigF
        + np.array(SigS.sum(axis=1)).ravel()
    )

    # Synthetic abstract XS — no physical energy grid (Phase E).
    mix = Mixture(
        SigC=SigC,
        SigL=SigL,
        SigF=SigF,
        SigP=SigP,
        SigT=SigT,
        SigS=[SigS],
        Sig2=[csr_matrix((ng, ng))],
        chi=np.array([1.0, 0.0]),
    )

    np.testing.assert_allclose(mix.SigC, SigC, rtol=1e-14)
    np.testing.assert_allclose(mix.SigT, SigT, rtol=1e-14)
    np.testing.assert_allclose(mix.SigF, SigF, rtol=1e-14)
    assert mix.ng == ng
