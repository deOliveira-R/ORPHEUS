r"""Garcia 2021 Case 1: the published scalar flux of a three-region sphere, the truth the characteristic reference reads.

R. D. M. Garcia, J. Comput. Phys. 433 (2021), Table 5, Case 1 (Williams 1991 Example 5), the converged ppP_N
column (the table's rightmost): radii (3.0, 5.0, 7.0) cm; core Sigma_t = 1.0, Sigma_s = 0.99; middle Sigma_t = 0.5,
Sigma_s = 0.30; outer Sigma_t = 2.0, Sigma_s = 1.90; internal sources (0.5, 1.0, 1.5) per cm^3 per steradian;
vacuum at r = 7. Garcia's "scalar flux" is :math:`\int_{-1}^{1}\Psi\,\mathrm d\mu` (no :math:`2\pi`) with a
per-steradian source, so a reference posed with the source as a rate reads half of it. Garcia verified this case
against Williams 1991 (integral-equation MoC) and Picca, Furfaro and Ganapol 2012 (S_N) to 3-4 significant figures
at every point.

Moved here in P1 step (e1b) from ``test_peierls_greens_function_garcia2021.py``, which step (e2) deleted with the
trajectory-resolvent family; the consumer is
``test_characteristic_reading.py::test_garcias_case_1_per_point``.
"""
from __future__ import annotations

import numpy as np

# Garcia 2021 Table 5 (Case 1 converged ppP_N; rightmost column).
GARCIA_2021_CASE1_R = np.array([
    0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0,
    3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0,
])
#: The table's values as printed (the strings carry the precision).
GARCIA_2021_CASE1_PHI_PRINTED = (
    "18.860", "18.756", "18.442", "17.911", "17.145", "16.095", "14.381",  # core
    "13.455", "13.337", "13.590", "14.361",                                # mid
    "15.532", "14.198", "10.807", "4.0763",                                # outer
)
GARCIA_2021_CASE1_PHI = np.array([float(v) for v in GARCIA_2021_CASE1_PHI_PRINTED])
#: Half a unit in each value's last printed digit, relative: Garcia's own error.
GARCIA_2021_CASE1_ROUNDING = np.array(
    [0.5 * 10.0 ** -len(v.split(".")[1]) for v in GARCIA_2021_CASE1_PHI_PRINTED]
) / GARCIA_2021_CASE1_PHI
