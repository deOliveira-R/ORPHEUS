"""Conservation probe for the rank-N hollow-sphere W matrix at σ_t = 0.

Addresses the 1.42 % plateau of Sanchez-McCormick-style rank-N closure
on hollow sphere (Issue #119).

Literature synthesis (2026-04-21): Ligou 1982, Stamm'ler 1983, Stacey 2007
all use scalar/DP-0 closures equivalent to ORPHEUS Phase F.4. None of
them presents a rank-N Legendre ladder. Sanchez 2002 uses piecewise-
constant angular sectors, not Legendre. The Sanchez-McCormick 1982
§III.F.1 Legendre-moment rank-N formulation is NOT cross-validated by
any successor textbook. See:
`.claude/agent-memory/literature-researcher/rank_n_closure_four_references_synthesis.md`.

**Recipe 1 test**: verify `W_oo[m,n] + W_io[m,n] = δ_{mn}` for the
shipped P̃_n-basis `compute_hollow_sph_transmission_rank_n` at σ_t = 0.

**Result**:
- At mode 0 (n=0): identity HOLDS exactly. W_oo[0,0] + W_io[0,0] = 1.
  (This is F.4's conservation that drives the 0.077 % scalar residual.)
- At modes n ≥ 1: identity FAILS. Diagonal sums are 0.28, 0.13, 0.09
  for n=1, 2, 3 — not 1, not even close.

**Why it fails at n ≥ 1**:

The failure is NOT a bug — it's a structural fact of the formulation.
At mode 0, P̃_0 = 1 is angle-independent, so P̃_0(c_in) = P̃_0(µ_emit) = 1
and the inner-surface angle remapping is invisible. At higher modes,
P̃_n(c_in) ≠ P̃_n(µ_emit) because c_in = √(1 − (R/r_0)²(1 − µ²_emit))
— the inner surface "sees" emission from a different angle than what
was emitted from outer. The W_io integrand therefore involves a basis
mismatch at arrival, and the naive conservation identity does not hold.

The TRUE conservation identity at σ_t = 0 involves the angle-remapping
Jacobian. A rank-N closure that correctly accounts for this mapping
would satisfy a harder identity — neither the shipped P̃_n W nor the
Sanchez µ-ortho W do so at n ≥ 1.

**Conclusion**: the rank-N hollow-sphere closure plateau is a
consequence of the angle-remapping at the inner surface not being
captured by the Legendre-moment basis. F.4 sidesteps this by living
at mode 0 where the remap is trivial. A correct rank-N closure would
need either (a) a basis that diagonalizes the c_in → µ_emit map (not
derived in the four-reference synthesis), or (b) a completely
different formulation like Sanchez 2002's piecewise-constant angular
sectors.

**Production**: F.4 scalar is the production closure for hollow
curvilinear cells. Its 0.077 % residual at σ_t·R = 5 is the
"quadrature floor" — see `diag_rank_n_closure_characterization.py`.

If promoted, this diagnostic becomes a pytest gate documenting the
conservation-identity structure for future reference.
"""
from __future__ import annotations
import sys
import numpy as np
import pytest

sys.path.insert(0, "/workspaces/ORPHEUS")

from orpheus.derivations.peierls_geometry import (
    compute_hollow_sph_transmission,
    compute_hollow_sph_transmission_rank_n,
)


R, r_0 = 5.0, 1.5


def test_f4_mode_0_conservation_holds_exactly_at_sigt_zero():
    """W_oo + W_io = 1 for mode 0 of the shipped W at σ_t = 0."""
    sig_t_zero = np.array([0.0])
    radii = np.array([R])
    W = compute_hollow_sph_transmission(r_0, R, radii, sig_t_zero, dps=20)
    # W = [[W_oo_scalar, W_oi_scalar], [W_io_scalar, W_ii]]
    assert abs(W[0, 0] + W[1, 0] - 1.0) < 1e-12, (
        f"F.4 mode-0 conservation failed: W_oo + W_io = "
        f"{W[0, 0] + W[1, 0]:.15f} ≠ 1"
    )


def test_higher_modes_fail_naive_conservation_at_sigt_zero():
    """At n ≥ 1, W_oo[n,n] + W_io[n,n] ≠ 1 at σ_t = 0.

    This is NOT a bug — it's a structural consequence of the P̃_n(c_in)
    ≠ P̃_n(µ_emit) angle remapping at the inner surface. Documented as
    a ``plateau characterization`` gate: if a future fix alters this,
    retire this test.
    """
    sig_t_zero = np.array([0.0])
    radii = np.array([R])
    N = 3
    W = compute_hollow_sph_transmission_rank_n(
        r_0, R, radii, sig_t_zero, N, dps=20,
    )
    W_oo = W[:N, :N]
    W_io = W[N:, :N]
    # Mode 0 satisfies conservation:
    assert abs(W_oo[0, 0] + W_io[0, 0] - 1.0) < 1e-12
    # Modes 1, 2 do NOT satisfy naive conservation (diag sum ≠ 1):
    for n in (1, 2):
        diag_sum = W_oo[n, n] + W_io[n, n]
        assert abs(diag_sum - 1.0) > 0.5, (
            f"Mode n={n}: diag sum = {diag_sum:.6f} — unexpectedly close "
            f"to 1. The angle-remapping structure may have changed; "
            f"investigate."
        )


def main():
    print(f"Rank-N W conservation probe at σ_t = 0, R={R}, r_0/R={r_0/R}")
    print()
    W_sca = compute_hollow_sph_transmission(
        r_0, R, np.array([R]), np.array([0.0]), dps=20,
    )
    print(f"Scalar F.4 W:\n{W_sca}")
    print(f"  W_oo + W_io = {W_sca[0,0] + W_sca[1,0]:.10f} (F.4 identity)")
    print()
    for N in (1, 2, 3, 4):
        W = compute_hollow_sph_transmission_rank_n(
            r_0, R, np.array([R]), np.array([0.0]), N, dps=20,
        )
        W_oo = W[:N, :N]
        W_io = W[N:, :N]
        print(f"N={N}: diagonal sums W_oo[n,n] + W_io[n,n]:")
        for n in range(N):
            print(f"  n={n}: {W_oo[n,n]:.6f} + {W_io[n,n]:.6f} = "
                  f"{W_oo[n,n] + W_io[n,n]:.6f}")


if __name__ == "__main__":
    main()
