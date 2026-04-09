# ORPHEUS Bug Patterns

Recurrent patterns that shape what QA checks. Each maps to a
concrete diagnostic.

## Cylindrical DD (ERR-014/015)

Wrong alpha recursion + missing delta-A/w. Passed homogeneous exact,
particle balance, conservation, non-negativity. Failed ONLY on
heterogeneous convergence (keff diverged: 1.15 -> 0.90 -> 0.52).

**Diagnostic**: Curvilinear transport needs fixed-source flat-flux test
(uniform Q, Sigma_t=1, 40+ cells). Check phi ~ Q/Sigma_t AND bounded
flux range.

## MOC Weight Bug (ERR-019)

Weight formula `4pi * w_a * w_p * t_s * sin(theta_p)` -- the 4pi and
sin(theta_p) cancel in keff for homogeneous. Only heterogeneous
multi-region exposes the error.

## Cross-Section Convention

`Mixture.SigS[l][g_from, g_to]` -- source uses transpose:
`Q = SigS^T @ phi`. Swap gives wrong multi-group, correct 1-group.

## Quadrature Weight Sums

GL: `sum(w) = 2`. Lebedev/LS/Product: `sum(w) = 4pi`. Wrong sum
cancels in keff -- catch with `phi = Q/Sigma_t` normalization test.

## The Alpha Dome

Curvilinear alpha coefficients: non-negative dome (0 -> peak -> 0).
Negative alpha -> NaN. Test: `assert np.all(alpha >= -1e-14)`.
