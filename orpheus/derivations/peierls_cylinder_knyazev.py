r"""SymPy derivation of the 3-D-corrected rank-:math:`N` cylinder
:math:`P_{\rm esc}^{(n,3d)}` and :math:`G_{\rm bc}^{(n,3d)}` primitives
(Knyazev :math:`\mathrm{Ki}_{2+k}` expansion).

The same Knyazev expansion is the **shared kernel** behind two BC
families: the curvilinear specular closure
(:func:`compute_P_esc_cylinder_3d_mode` /
:func:`compute_G_bc_cylinder_3d_mode` for ``boundary="specular"``) and
the corrected 3-D Hébert white-BC kernel (Issue #112 Phase C — see
:mod:`orpheus.derivations.peierls_cylinder_g_bc_3d`). For that reason
this module sits **outside** the ``peierls_specular`` sub-package and
serves both clients.

Background — why mode-N cylinder needs special handling
-------------------------------------------------------

The legacy rank-N primitives evaluate the shifted Legendre
:math:`\tilde P_n(\mu_{\rm exit})` weight at the **in-plane cosine**
:math:`\mu_{\rm 2D} = (r_i\cos\omega + \rho_{\max})/R` and use the
:func:`escape_kernel_mp`-style :math:`\mathrm{Ki}_2(\tau_{\rm 2D})`
absorption of the polar integral. This is consistent only at
:math:`n = 0` (where :math:`\tilde P_0 \equiv 1` and the polar integral
trivially gives :math:`\mathrm{Ki}_2`).

For :math:`n \ge 1`, the **3-D direction cosine** at the lateral
surface is :math:`\mu_{\rm 3D} = \sin\theta_p\,\mu_{\rm 2D}` (the
in-plane component scaled by the polar projection). Inserting this
into the rank-:math:`N` partial-current moment integral and expanding
:math:`\tilde P_n` in the monomial basis :math:`\tilde P_n(y) =
\sum_k c_n^k\,y^k`, the polar integral reduces to a Bickley-Naylor
function via the substitution :math:`u = \pi/2 - \theta_p`,
:math:`\sin\theta_p = \cos u`,

.. math::
   :label: peierls-cyl-knyazev-polar-id

   \int_0^{\pi/2} \sin^{k+1}\theta_p\,
       e^{-\tau_{\rm 2D}/\sin\theta_p}\,\mathrm d\theta_p
   \;=\; \int_0^{\pi/2} \cos^{k+1}u\,e^{-\tau_{\rm 2D}/\cos u}\,\mathrm du
   \;=\; \mathrm{Ki}_{k+2}(\tau_{\rm 2D})

(standard Bickley definition :math:`\mathrm{Ki}_n(x) = \int_0^{\pi/2}
\cos^{n-1}u\,e^{-x/\cos u}\,\mathrm du`). This is the Knyazev
:math:`\mathrm{Ki}_{2+k}` expansion at orders :math:`k = 0, 1, \ldots, n`.

The result is the **rank-:math:`n` Knyazev expansion** of the cylinder
primitives:

.. math::

   P_{\rm esc}^{(n,3d)}(r_i) \;=\; \frac{1}{\pi}\!\int_0^\pi
       \sum_{k=0}^n c_n^k\,\mu_{\rm 2D}(\omega)^k\,
       \mathrm{Ki}_{k+2}\!\bigl(\tau_{\rm 2D}(\omega)\bigr)\,\mathrm d\omega

   G_{\rm bc}^{(n,3d)}(r_i) \;=\; \frac{4}{\pi}\!\int_0^\pi
       \sum_{k=0}^n c_n^k\,\mu_{\rm 2D}(\omega)^k\,
       \mathrm{Ki}_{k+2}\!\bigl(\tau_{\rm 2D}(\omega)\bigr)\,\mathrm d\omega

with the same :math:`\mu_{\rm 2D}, \tau_{\rm 2D}` integrand kernel —
:math:`P` and :math:`G` differ only by the :math:`\frac{1}{\pi}` vs
:math:`\frac{4}{\pi}` prefactor (the latter from the
:math:`(b_n / \pi) \cdot 4` factor for the inward-distribution
:math:`b_n` convention; see :func:`derive_g_prefactor`).

For :math:`n = 0`, both reduce to the existing
:func:`compute_P_esc` / :func:`compute_G_bc_cylinder_3d` primitives
because :math:`\tilde P_0 = 1`, :math:`c_0^0 = 1`, and only the
:math:`k = 0` term survives.

Production breadcrumb
---------------------

The shipped primitives live in :mod:`peierls_geometry` as
:func:`compute_P_esc_cylinder_3d_mode` and
:func:`compute_G_bc_cylinder_3d_mode`. They are pinned to this
SymPy origin by
``tests/derivations/test_peierls_cylinder_knyazev_symbolic.py``
which numerically verifies the polar integral identity
(:eq:`peierls-cyl-knyazev-polar-id`) and bit-equates the production
output with the lambdified symbolic expansion.
"""

from __future__ import annotations

import sympy as sp

from ._shifted_legendre import shifted_legendre_monomial_coefs


def derive_polar_integral_identity(k: int) -> dict:
    r"""Return SymPy expressions for the polar integral identity at
    order :math:`k`.

    Computes (symbolically) the substitution chain that lets the
    polar integral reduce to a Bickley function:

    .. math::

       \int_0^{\pi/2} \sin^{k+1}\theta\,e^{-x/\sin\theta}\,
           \mathrm d\theta
       \;=\; \int_0^{\pi/2} \cos^{k+1}u\,e^{-x/\cos u}\,\mathrm du
       \;=\; \mathrm{Ki}_{k+2}(x).

    Returns a dict with ``theta_form`` (the :math:`\sin` integrand
    in :math:`\theta_p`), ``u_form`` (after the
    :math:`u = \pi/2 - \theta_p` substitution), and the symbol
    ``x`` for downstream lambdification. Numerical verification
    against :func:`ki_n_mp` lives in the test suite — see
    ``test_polar_integral_identity``.
    """
    theta_p, u, x = sp.symbols("theta_p u x", positive=True)
    theta_form = sp.sin(theta_p) ** (k + 1) * sp.exp(-x / sp.sin(theta_p))
    u_form = sp.cos(u) ** (k + 1) * sp.exp(-x / sp.cos(u))
    return {
        "k": k,
        "theta_form": theta_form,
        "u_form": u_form,
        "theta_p": theta_p,
        "u": u,
        "x": x,
    }


def derive_p_prefactor() -> sp.Expr:
    r"""Return the symbolic :math:`1/\pi` P-prefactor for the cylinder
    mode-N :math:`P_{\rm esc}^{(n,3d)}` primitive.

    The angle-integrated partial-current moment from a unit isotropic
    source at :math:`r_i`:

    .. math::

       \bar J^{+}_n(r_i) \;=\; \int_{4\pi} \frac{e^{-\tau_{\rm 3D}}}{4\pi}\,
           \tilde P_n(\mu_{\rm 3D})\,\mathrm d\Omega_{\rm 3D}.

    With :math:`\mathrm d\Omega_{\rm 3D} = \sin\theta_p\,\mathrm d\theta_p\,
    \mathrm d\phi`, :math:`\mu_{\rm 3D} = \sin\theta_p\,\mu_{\rm 2D}(\alpha)`,
    :math:`\tau_{\rm 3D} = \tau_{\rm 2D}(\alpha)/\sin\theta_p`. The two
    discrete symmetries (:math:`\alpha \to 2\pi - \alpha` in-plane,
    :math:`\theta_p \to \pi - \theta_p` axial) together produce a factor
    of 4. Combined with the :math:`(4\pi)^{-1}` isotropic normalisation
    and the :math:`\sum_k c_n^k\,\mu_{\rm 2D}^k\,\mathrm{Ki}_{k+2}`
    polar reduction (eq. :eq:`peierls-cyl-knyazev-polar-id`):

    .. math::

       \bar J^{+}_n \;=\; \frac{1}{\pi}\!\int_0^\pi\!\sum_k c_n^k\,
                   \mu_{\rm 2D}^k\,\mathrm{Ki}_{k+2}(\tau_{\rm 2D})\,
                   \mathrm d\alpha.

    Returns the SymPy expression :math:`1/\pi`.
    """
    return sp.Rational(1) / sp.pi


def derive_g_prefactor() -> sp.Expr:
    r"""Return the symbolic :math:`4/\pi` G-prefactor for the cylinder
    mode-N :math:`G_{\rm bc}^{(n,3d)}` primitive.

    With the inward distribution
    :math:`\psi^{-}(\mu) = (b_n/\pi)\,\tilde P_n(\mu)` integrated over
    the full inward 4π solid angle:

    .. math::

       \varphi(r_i) \;=\; \frac{b_n}{\pi}\!\int_{4\pi}\!
           \tilde P_n(\mu_{\rm 3D})\,e^{-\tau_{\rm 3D}}\,
           \mathrm d\Omega_{\rm 3D}.

    Both :math:`\alpha` and :math:`\theta_p` symmetries each give a
    factor of 2 (so :math:`2 \cdot 2 = 4` overall) over the
    half-quadrant :math:`\alpha \in [0,\pi]`,
    :math:`\theta_p \in [0,\pi/2]`. Combined with the
    :math:`(b_n/\pi)` source amplitude:

    .. math::

       \varphi(r_i) \;=\; \frac{4\,b_n}{\pi}\!\int_0^\pi\!
           \sum_k c_n^k\,\mu_{\rm 2D}^k\,\mathrm{Ki}_{k+2}(\tau_{\rm 2D})\,
           \mathrm d\alpha.

    The factor :math:`G^{(n,3d)}_{\rm cyl} = \varphi/b_n` therefore has
    leading prefactor :math:`4/\pi`.

    Returns the SymPy expression :math:`4/\pi`.
    """
    return sp.Rational(4) / sp.pi


def shifted_legendre_monomial_table(
    N_max: int,
) -> dict[int, tuple[float, ...]]:
    r"""Return :math:`\{n: (c_n^0, \ldots, c_n^n)\}` for ``n = 0..N_max``.

    Thin wrapper over the canonical
    :func:`orpheus.derivations._shifted_legendre.shifted_legendre_monomial_coefs`
    that surfaces the coefficient ladder used to expand
    :math:`\tilde P_n(\mu) = \sum_k c_n^k\,\mu^k`. The Knyazev
    expansion uses these coefficients directly.
    """
    return {n: shifted_legendre_monomial_coefs(n) for n in range(N_max + 1)}
