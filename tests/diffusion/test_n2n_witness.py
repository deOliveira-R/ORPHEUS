r"""The diffusion family's (n,2n) witness (#428 F-2, landed with #426 step 2).

**The gap this closes.** `[M]` 2026-09-03 (`qa`, ``scratch/_428_four_solver_check.md``
§F-2): under ν₂ₙ ``2 → 1`` **and** ``2 → 0``, ``tests/diffusion`` read
**113 passed / 0 red**.  Cause, measured with a ``Mixture.__post_init__``
census over the run: **625** mixtures constructed, **1** with a non-zero
``Sig2``, and that one built as a side effect of the derivations registry and
never handed to a diffusion solve.  Every diffusion fixture comes from
``_fixture_materials()`` (``tests/diffusion/test_operators.py``), which calls
``make_mixture(...)`` with no ``sig_2=`` — and `[M]` ``make_mixture`` defaults
``sig_2`` to all-zero (lessons §3: the builder that silently nulls a channel).
``IsotropicN2N`` appears in 4 diffusion tests and is a ZERO kernel in all
four: Mode-10 *exercised-but-unconstrained* at family scale.

**Claim.** Layer = EIGENVALUE.  Pillar = CLOSED-FORM (the 2-group infinite
medium ``k_inf = λ_max(A⁻¹F)`` with ``A = diag(Σ_t) − (Σ_s + 2Σ_2)ᵀ``).  Claim
kind = REFERENCE.

**Structural independence — and where it stops.** The reference is
``orpheus/derivations/continuous/analytical/homogeneous.py::derive_2g_n2n``,
whose eigenvalue comes from ``derivations/common/eigenvalue.py``, which
spells its OWN literal ``2.0 *`` (#428 F-5 — a NAMED exclusion of the
multiplicity census, on purpose).  That independent literal is exactly why
the production mutation reddens this row instead of moving both sides
together.  What is NOT independent: both sides read the same ``Mixture``, so
this is not a check of the cross sections.

**Measured, `[M]` 2026-09-03 on branch ``fix/n2n-anisotropy``:**

======================  ==========================  =============
ν₂ₙ                     k (4-cell reflective slab)  rel vs k_inf
======================  ==========================  =============
2  (production)         1.6532258064516114          **2.686e-16**
1  (the ERR-015 class)  1.3194444444444438          2.019e-01
0  (channel dropped)    1.0377026074700486          3.723e-01
3  (over-count)         2.051643192488262           2.410e-01
======================  ==========================  =============

⟹ 0.005 s, and the gate's ``rtol = 1e-10`` sits **9 orders** below the weakest
mutation.  It is NOT ``slow``, and it must never become so (#428 F-4's lesson:
a channel whose only catcher is ``slow``-marked has zero coverage under the
canonical ``-m "not slow"`` gate).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations import get
from orpheus.diffusion import solve_diffusion_1d
from orpheus.geometry import BC, CoordSystem, Mesh1D

pytestmark = pytest.mark.l1

_CASE = "homo_2eg_n2n"


@pytest.fixture(scope="module")
def case():
    return get(_CASE)


@pytest.fixture(scope="module")
def reflective_mesh():
    """Zero-leakage 1-D slab: the discrete eigenvalue IS ``k_inf`` (no
    buckling, no discretisation error to account for).  4 cells is enough —
    the answer is spatially flat; more cells only cost time."""
    return Mesh1D(
        edges=np.linspace(0.0, 10.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )


@pytest.mark.verifies("matrix-eigenvalue", "removal-matrix", "fission-matrix")
def test_diffusion_reproduces_the_n2n_infinite_medium_eigenvalue(
    case, reflective_mesh,
):
    r"""``solve_diffusion_1d`` on a Σ₂ₙ ≠ 0 medium reaches the closed form.

    ACTIVATION leg first (`vv` #19 / lessons §3): assert the fixture's channel
    is actually live, or the row is the ``Sig2 = 0`` tautology every other
    diffusion fixture already is.
    """
    mix = case.materials[0]
    if mix.Sig2[0].nnz == 0:
        pytest.fail(
            f"{_CASE} lost its (n,2n) channel — this row would then be a "
            f"Sig2 = 0 tautology, which is exactly the #428 F-2 gap it exists "
            f"to close"
        )
    result = solve_diffusion_1d(
        {0: mix}, reflective_mesh,
        keff_tol=1e-13, flux_tol=1e-12, max_outer=2000,
    )
    np.testing.assert_allclose(
        result.keff, case.k_inf, rtol=1e-10,
        err_msg=(
            f"diffusion k disagrees with the closed-form k_inf on the "
            f"Σ₂ₙ ≠ 0 medium: got {result.keff!r}, expected {case.k_inf!r}. "
            f"[M] the honest value reproduces to 2.7e-16; ν₂ₙ = 1 reads "
            f"1.3194444444444438 (20.2 % low)."
        ),
    )


def test_the_diffusion_loss_operator_carries_the_channel(case, reflective_mesh):
    r"""The un-wired-arm control, and the reason F-2 was silent.

    `[M]` ``tests/diffusion/test_operators.py``'s ``_loss()`` helper
    assembles ``A = L + C − S − B`` **without** the N2N arm, while production's
    ``DiffusionSolver`` assembles ``L + C − (S + N2N) − B``
    (``orpheus/diffusion/solver.py``).  The two are bit-identical on
    ``Sig2 = 0`` — a Pattern-2 twin that only wakes when a Σ₂ₙ ≠ 0 fixture
    lands, i.e. now.  This row pins the DIFFERENCE: dropping the channel from
    the medium must move k by an amount far outside the solver tolerance, so
    the twin cannot be "corrected" into agreement.
    """
    from dataclasses import replace
    from scipy.sparse import csr_matrix

    mix = case.materials[0]
    ng = len(mix.SigT)
    # Same Σ_t (removal unchanged), channel emission removed: this is NOT the
    # physical "no (n,2n)" medium — it is the mutation the twin cannot see.
    without = replace(mix, Sig2=[csr_matrix((ng, ng))])
    k_with = solve_diffusion_1d(
        {0: mix}, reflective_mesh, keff_tol=1e-13, flux_tol=1e-12, max_outer=2000,
    ).keff
    k_without = solve_diffusion_1d(
        {0: without}, reflective_mesh, keff_tol=1e-13, flux_tol=1e-12,
        max_outer=2000,
    ).keff
    rel = abs(k_with - k_without) / k_with
    assert rel > 1e-3, (
        f"removing the (n,2n) emission moved k by only {rel:.3e} — the "
        f"diffusion loss operator is not consuming the channel, and every "
        f"row above it is a Sig2 = 0 tautology"
    )
