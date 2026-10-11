r"""The Peierls Nystrom specular primitives at rank 1: T_00 = P_ss (sphere, cylinder) and T_oi = 2 E_3 (slab).

Re-homed in step (e2) of the characteristic-reference campaign
(``.claude/plans/characteristic_reference_architecture.md``) from three files the
step deleted with the trajectory-resolvent family
(``test_peierls_greens_function_cylinder_solver.py``, ``..._slab_solver.py``,
``..._xverif.py``). The rows never read that family: they gate the SURVIVING
production primitives of :mod:`orpheus.derivations.continuous.peierls_nystrom.geometry`
(``compute_T_specular_{sphere,slab,cylinder_3d}`` against
``compute_P_ss_{sphere,cylinder}`` and mpmath's ``E_3``), so deleting them would
have left those primitives without a live gate (``retirement-audit`` C.12; the
other cylinder witness, ``test_peierls_specular_bc.py::test_specular_multibounce_cyl_rank1_equals_hebert``,
is withdrawn under #506). The bodies are unchanged; at ``196f5215`` none of the
three carried a withdrawal mark, and none carries one here.

The ``verifies`` labels (``peierls-greens-cylinder-T``, ``peierls-greens-slab-T``)
live on ``docs/theory/references/trajectory_resolvent.rst``, which step (f)
rewrites: the archivist keeps or re-homes them.
"""
from __future__ import annotations

import mpmath
import numpy as np
import pytest

from orpheus.derivations.continuous.peierls_nystrom.geometry import (
    compute_P_ss_cylinder,
    compute_P_ss_sphere,
    compute_T_specular_cylinder_3d,
    compute_T_specular_slab,
    compute_T_specular_sphere,
)


@pytest.mark.l1
@pytest.mark.parametrize(
    "tau_R",
    [0.5, 1.0, 2.5, 5.0, 10.0],
    ids=["tauR_0.5", "tauR_1.0", "tauR_2.5", "tauR_5.0", "tauR_10.0"],
)
def test_v_alpha2_sphere_T00_equals_Pss_via_production_primitives(tau_R):
    r"""V_α2 — **structurally-independent L1 cross-check**:
    :math:`T_{00}^{\rm sphere} = P_{ss}^{\rm sphere}` at the
    production-primitive level over a :math:`\tau_R` sweep.

    The SymPy V_α2 sphere proof
    (:func:`derive_T00_equals_P_ss_sphere`) reduces both sides to the
    same Hébert closed form via independent integration paths
    (transfer-matrix vs polar escape). This test pins the analogous
    identity at the **production-code level**: the matrix function
    :func:`compute_T_specular_sphere` and the scalar function
    :func:`compute_P_ss_sphere` are independent code paths (different
    functions, different intermediate quantities — matrix construction
    vs scalar P_ss). At rank-1 (``n_modes = 1``, isotropic mode
    :math:`\tilde P_0 = 1`) the matrix [0, 0] element MUST equal the
    scalar :math:`P_{ss}`.

    This is the V&V evidence the original integrand-tautology proof
    in SymPy could not provide. Each side comes from its own
    derivation in the production code; equality is a structurally-
    independent cross-check.
    """
    R = 5.0
    sig_t = np.array([tau_R / R])
    radii = np.array([R])

    P_ss_scalar = compute_P_ss_sphere(radii, sig_t, n_quad=64, dps=25)
    T_matrix = compute_T_specular_sphere(radii, sig_t, n_modes=1, n_quad=64)
    T_00 = float(T_matrix[0, 0])

    np.testing.assert_allclose(
        T_00, P_ss_scalar, rtol=1e-10, atol=1e-12,
        err_msg=(
            f"V_α2 sphere production-primitive cross-check failed at "
            f"τ_R = {tau_R}: T_00 = {T_00:.16e}, P_ss = "
            f"{P_ss_scalar:.16e}, diff = {abs(T_00 - P_ss_scalar):.3e}"
        ),
    )


@pytest.mark.l1
@pytest.mark.parametrize(
    "tau_L",
    [0.5, 1.0, 2.5, 5.0, 10.0],
    ids=["tauL_0.5", "tauL_1.0", "tauL_2.5", "tauL_5.0", "tauL_10.0"],
)
@pytest.mark.verifies("peierls-greens-slab-T")
def test_v_alpha2_slab_T00_equals_2E3_via_production_primitive(tau_L):
    r"""V_α2_slab — **structurally-independent L1 cross-check**:
    :math:`T_{\rm slab}[0, n_modes] = 2 E_3(\Sigma_t L)` via the
    production primitive :func:`compute_T_specular_slab` over a
    :math:`\tau_L` sweep.

    The slab T-matrix is block-off-diagonal:

    .. math::

       T_{\rm slab} = \begin{pmatrix} 0 & T_{oi} \\ T_{io} & 0
                      \end{pmatrix}

    with :math:`T_{io} = T_{oi}` by face symmetry on a homogeneous slab
    and :math:`T_{oi}^{(0,0)} = 2 E_3(\Sigma_t L)` at rank-1
    (``n_modes = 1``).

    The production primitive computes :math:`T_{oi}^{(0,0)}` via
    explicit half-range Gauss-Legendre on :math:`\mu \in [0, 1]`
    integrating :math:`2 \mu e^{-\Sigma_t L /\mu}`. The mpmath E_3
    is computed via ``mpmath.expint(3, τ)`` — a fundamentally
    different code path (analytical continuation of the gamma
    function). Equality of the two is the V_α2_slab numerical
    structural-independence cross-check.
    """
    L = 1.0
    sig_t = np.array([tau_L / L])  # ensures τ_L total = sig_t · L
    radii = np.array([L])

    T_matrix = compute_T_specular_slab(radii, sig_t, n_modes=1, n_quad=64)
    # Off-diagonal block T_oi (the slab face-to-face transit):
    T_oi_00 = float(T_matrix[0, 1])  # block [outer, inner] entry
    # Also check T_io = T_oi (face symmetry):
    T_io_00 = float(T_matrix[1, 0])

    canonical = float(2.0 * mpmath.expint(3, tau_L))

    np.testing.assert_allclose(
        T_oi_00, canonical, rtol=1e-10, atol=1e-12,
        err_msg=(
            f"V_α2_slab T_oi production-primitive cross-check failed at "
            f"τ_L = {tau_L}: T_oi[0,0] = {T_oi_00:.16e}, 2·E_3 = "
            f"{canonical:.16e}, diff = {abs(T_oi_00 - canonical):.3e}"
        ),
    )
    np.testing.assert_allclose(
        T_io_00, canonical, rtol=1e-10, atol=1e-12,
        err_msg=(
            f"V_α2_slab T_io production-primitive cross-check failed at "
            f"τ_L = {tau_L}: T_io[0,0] = {T_io_00:.16e}, 2·E_3 = "
            f"{canonical:.16e}"
        ),
    )


@pytest.mark.l1
@pytest.mark.parametrize(
    "tau_R",
    [0.5, 1.0, 2.5, 5.0, 10.0],
    ids=["tauR_0.5", "tauR_1.0", "tauR_2.5", "tauR_5.0", "tauR_10.0"],
)
@pytest.mark.verifies("peierls-greens-cylinder-T")
def test_v_alpha2_cyl_T00_equals_Pss_via_production_primitives(tau_R):
    r"""V_α2_cyl — **structurally-independent L1 cross-check**:
    :math:`T_{00}^{\rm cyl} = P_{ss}^{\rm cyl}` at the production-
    primitive level over a :math:`\tau_R` sweep.

    Why this test is load-bearing for cylinder V_α2: unlike sphere,
    cylinder has no elementary closed form (Ki_3 obstruction), so the
    SymPy V_α2_cyl proof
    (:func:`derive_T00_equals_P_ss_cylinder`) can only reach the
    integrand-identity level via the Bickley-Naylor bridge
    :math:`\mathrm{Ki}_3(x) = \int_0^{\pi/2}\sin^2\beta e^{-x/\sin\beta}\mathrm d\beta`
    — a known mathematical fact applied as a substitution, not a
    SymPy-derived equality. The rigorous V&V evidence for the rank-1
    Knyazev ≡ Hébert white-BC theorem must therefore come from the
    **numerical level**.

    This test pins the identity at the production-code level: the
    matrix function :func:`compute_T_specular_cylinder_3d` (Knyazev
    expansion + scattering kernel + matrix construction) and the
    scalar function :func:`compute_P_ss_cylinder` (slanted-chord polar
    integration → Ki_3 → in-plane :math:`\alpha` integration) are
    independent code paths — different functions, different
    intermediate quantities. At rank-1 (``n_modes = 1``,
    :math:`m = n = 0`, all Knyazev shifted-Legendre coefficients
    :math:`k_m = k_n = 0`) the matrix [0, 0] element MUST equal the
    scalar :math:`P_{ss}^{\rm cyl}`.

    The L1 evidence pillar for cylinder V_α2 is therefore this
    numerical cross-check, not the SymPy integrand identity (which
    must rely on the Bickley-Naylor bridge as a known identity).
    """
    R = 5.0
    sig_t = np.array([tau_R / R])
    radii = np.array([R])

    P_ss_scalar = compute_P_ss_cylinder(radii, sig_t, n_quad=64, dps=25)
    T_matrix = compute_T_specular_cylinder_3d(
        radii, sig_t, n_modes=1, n_quad=64,
    )
    T_00 = float(T_matrix[0, 0])

    np.testing.assert_allclose(
        T_00, P_ss_scalar, rtol=1e-10, atol=1e-12,
        err_msg=(
            f"V_α2_cyl production-primitive cross-check failed at "
            f"τ_R = {tau_R}: T_00 = {T_00:.16e}, P_ss = "
            f"{P_ss_scalar:.16e}, diff = {abs(T_00 - P_ss_scalar):.3e}"
        ),
    )
