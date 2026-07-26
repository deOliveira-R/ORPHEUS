r"""Pin the consistent-DSA four-step derivation — Phase 3a (#2).

:mod:`orpheus.derivations.discrete.sn.dsa` is the algebra of record for
the slab weighted-diamond consistent low-order (synthetic diffusion)
operator: it EXECUTES Larsen's four-step symbolically on a general
symmetric quadrature and proves the printed operator — the interior
tridiagonal (27) with coefficients (23a–f), the updates (28), the
Marshak/reflecting boundary rows — rather than transcribing it. The
production low-order build (Phase 3b) implements exactly what these
theorems prove, and the D-gates cite them — so the theorems themselves
must stay green under the tree: this file runs each derivation as a test.

Every theorem verifies its identity in EXACT symbolic arithmetic inside
the module and raises
:class:`~orpheus.derivations.discrete.sn.dsa.DerivationError` on any
failure — there is no tolerance anywhere in this file. One test per
theorem so a regression localizes to the broken identity.

vv Mode-8 note: the module uses explicit raises (never bare ``assert``),
so the theorems fire under the canonical ``python -O`` invocation; the
tests below simply propagate.
"""
from __future__ import annotations

import pytest

from orpheus.derivations.discrete.sn import dsa

pytestmark = pytest.mark.foundation


class TestFourStepDerivation:
    """The step-by-step battery (each raises ``DerivationError`` on failure)."""

    def test_closure_weight_annihilation_identities(self):
        r"""L0γ = L1γ = L0β = L1β = 0 for every symmetric quadrature —
        with Larsen's (14b) "3" exposed as 1/W2."""
        dsa.derive_closure_weight_identities()

    def test_p_recursion_lemma(self):
        r"""L1[μv] = (2/3)L2[v] + (1/3)L0[v] identically — the
        quadrature-INDEPENDENT ⅓ mechanism, kept distinct from W2."""
        dsa.prove_p_recursion_lemma()

    def test_moment_equations_split_exact(self):
        r"""(16a–d): the explicit-moments + opaque-lagged-functionals form
        equals the raw angular reductions for the symbolic quadrature."""
        dsa.derive_moment_equations()

    def test_main_theorem_interior_row_is_larsen_27(self):
        r"""THE MAIN THEOREM: shared-edge f1 continuity ≡ Larsen (27) with
        (23a–f) under the printed convention (W0=1, W2=1/3), proven by
        exact coefficient proportionality on the constraint variety."""
        dsa.derive_interior_row()

    def test_update_relations_are_larsen_28(self):
        r"""The cell-average updates reproduce (28a/28b)."""
        dsa.derive_update_relations()

    def test_boundary_rows(self):
        r"""Marshak (38a) coefficients DERIVED as (γ_N, W2⁺/W2) — printed
        (γ_N, 1/2) — and the vacuum row closes one-sidedly onto f0;
        reflecting is f1 = 0."""
        rows = dsa.derive_boundary_rows()
        if set(rows) != {"marshak_open", "vacuum_closed", "reflecting"}:
            raise dsa.DerivationError("boundary-row set changed unexpectedly")

    def test_dd_instance_coefficients(self):
        r"""The diamond member (α = 0): σ̂_R = σ_T − σ_S0,
        D = 1/[3(σ_T − σ_S1)] — the operator Phase 3b wires."""
        dsa.derive_dd_instance()

    def test_one_sided_f1_forms(self):
        r"""The compact (25)/(26) closures the numeric builder realizes
        equal the solved edge f1 values — the builder transcribes
        nothing."""
        dsa.derive_one_sided_f1_forms()
