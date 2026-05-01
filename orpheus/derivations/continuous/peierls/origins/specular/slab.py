r"""SymPy derivation of slab per-face mode-:math:`N` primitives for
the specular boundary condition.

Companion to :mod:`.r_matrix` (the curvilinear rank-:math:`N` reflection
operator). The slab geometry has TWO planar boundary surfaces (at
:math:`x = 0` and :math:`x = L`); per-face mode decomposition gives the
mode space :math:`A_{\rm slab} = A_{\rm outer} \oplus A_{\rm inner}
= \mathbb R^{2N}` and the reflection operator becomes block-diagonal:

.. math::

    R_{\rm slab} \;=\; \begin{pmatrix}
        R_{\rm spec}(N) & 0 \\
        0 & R_{\rm spec}(N)
    \end{pmatrix}_{2N \times 2N},
    \qquad R_{\rm spec}(N) \;=\; \tfrac{1}{2}\,M^{-1}.

The block-diagonal structure follows from specular being **local at
each face** — a particle leaving the outer face cannot magically appear
at the inner face except through the volume kernel, which is already
captured in :math:`K_{\rm vol}`.

Per-face primitives in closed form
----------------------------------

For slab the µ-integrals reduce to **exponential integrals**
:math:`E_n(\tau)` because the path optical depth factors as
:math:`\tau_{\rm path} = \tau_{\rm perp}/\mu`. Substituting
:math:`t = 1/\mu` shows
:math:`\int_0^1 \mu^k\,e^{-\tau_{\rm perp}/\mu}\,\mathrm d\mu
= E_{k+2}(\tau_{\rm perp})`. Combined with the shifted-Legendre
monomial expansion :math:`\tilde P_n(\mu) = \sum_k c_n^k\,\mu^k`:

.. math::
   :label: peierls-slab-spec-Pesc-no-mu-weight

   P_{\rm face}^{(n)}(x_i)
       \;=\; \tfrac{1}{2}\!\int_0^1 \tilde P_n(\mu)\,
                e^{-\tau_{\rm face}(x_i)/\mu}\,\mathrm d\mu
       \;=\; \tfrac{1}{2}\sum_{k=0}^{n} c_n^k\,
                E_{k+2}\!\bigl(\tau_{\rm face}(x_i)\bigr)

.. math::
   :label: peierls-slab-spec-Gbc-no-mu-weight

   G_{\rm face}^{(n)}(x_i)
       \;=\; 2\!\int_0^1 \tilde P_n(\mu)\,e^{-\tau_{\rm face}(x_i)/\mu}\,
                \mathrm d\mu
       \;=\; 2\sum_{k=0}^{n} c_n^k\,
                E_{k+2}\!\bigl(\tau_{\rm face}(x_i)\bigr).

The inner face uses :math:`\tau_{\rm in}(x_i) = \int_0^{x_i}
\Sigma_t(x')\,\mathrm dx'` (perpendicular optical depth from
:math:`x_i` to the face at :math:`x = 0`); same formula.

Basis-choice resolution — the no-:math:`\mu`-weight basis is canonical
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are two reasonable bases for the per-face slab specular
expansion:

- **No-:math:`\mu`-weight basis** (canonical here): both P and G use
  :math:`E_{k+2}(\tau)`. Mode 0 reduces to :math:`(1/2)\,E_2(\tau)`
  for P and :math:`2\,E_2(\tau)` for G — the **same form** as
  :func:`compute_P_esc_outer` and :func:`compute_G_bc_outer` (slab
  branch) at :math:`n = 0`. Used by the shipped
  ``boundary="specular"`` slab path.

- **µ-weighted basis** (rejected): P uses :math:`E_{k+3}(\tau)` (with
  one extra :math:`\mu` factor for the partial-current-moment cosine
  weight), G stays :math:`E_{k+2}(\tau)`. This basis has its own
  internal consistency, but mode 0 gives :math:`(1/2)\,E_3(\tau)` for
  P which does NOT match :func:`compute_P_esc_outer` — it would
  require a corresponding R rescaling and a different convention for
  the cosine weight inside :math:`R_{\rm spec}`.

The two bases are mathematically equivalent (same operator under a
basis change), but only the no-:math:`\mu`-weight one is locked-in by
this module because the production code already uses it at :math:`n =
0` (matching :func:`compute_P_esc_outer`). The :math:`n \ge 1`
production primitives are not yet implemented for slab; when a future
implementer ships them, the no-:math:`\mu`-weight basis is the
canonical contract — they should NOT switch to the µ-weighted form
without also reworking the cosine weight inside :math:`R_{\rm spec}`
and the slab-branch reduction at :math:`n = 0`.

Critical implementation gotcha — the divisor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The per-face G primitive must use the **single-face surface area
divisor** (= 1 for slab; each face has unit transverse area), NOT the
legacy combined-face divisor (= 2 from
:meth:`CurvilinearGeometry.rank1_surface_divisor`). The combined-face
divisor=2 is correct for legacy Mark, which uses ONE primitive
representing the SUM over both faces and then normalises by the TOTAL
surface area = 2.

The per-face decomposition keeps the two faces algebraically separate.
Each face block is a self-contained "outgoing × R_face × inward"
factorisation with the face's individual area = 1.

This was the entire root cause of the 2026-04-27 "rank-N specular slab
plateaus at -5%/+7%" investigation. With divisor=2 (incorrectly
inherited from the legacy Mark code path), the per-face block-diagonal
specular reflects only HALF the outgoing current at each face, leading
to the observed 5-8 % undershoot. With the correct divisor=1,
:math:`R_{\rm face} = (1/2)\,M^{-1}` recovers the rank-N specular and
converges monotonically to :math:`k_\infty` as :math:`N \to \infty`
AND as the spatial mesh refines.

Verification
------------

Pinned numerically against the production
:func:`compute_P_esc_outer`/:func:`compute_G_bc_outer` (slab branches)
at mode 0, plus the algebraic
:math:`2\,M_{\rm blockdiag}\,R_{\rm slab} - I = 0` per-face contract,
by ``tests/derivations/test_peierls_specular_slab_symbolic.py``.

Future work — :math:`n \ge 1` slab specular
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The shipped slab specular path uses :math:`n = 0` exclusively
(:func:`compute_P_esc_outer` / :func:`compute_G_bc_outer`). Higher
modes for slab are not implemented in production. When a future
implementer wires :math:`n \ge 1`, this module's
:eq:`peierls-slab-spec-Pesc-no-mu-weight` /
:eq:`peierls-slab-spec-Gbc-no-mu-weight` are the canonical contract —
the no-:math:`\mu`-weight basis is locked in.
"""

from __future__ import annotations

import sympy as sp

from .....common.shifted_legendre import shifted_legendre_monomial_coefs
from .r_matrix import build_R_specular_symbolic


def derive_slab_p_outer_n0() -> sp.Expr:
    r"""Mode-0 slab outer P primitive in closed form.

    .. math::

       P_{\rm esc, out}^{(0)}(x_i) \;=\; \tfrac{1}{2}\,E_2(\tau_{\rm out}(x_i)).

    Returned as the SymPy expression :math:`(1/2)\,E_2(\tau)` — pure
    symbol, no numerical specialisation. The substitution
    :math:`t = 1/\mu` reduces :math:`\int_0^1 e^{-\tau/\mu}\,\mathrm d\mu
    = E_2(\tau)` from the standard exponential-integral identity.
    """
    tau = sp.Symbol("tau", positive=True)
    return sp.Rational(1, 2) * sp.expint(2, tau)


def derive_slab_g_outer_n0() -> sp.Expr:
    r"""Mode-0 slab outer G primitive in closed form.

    .. math::

       G_{\rm bc, out}^{(0)}(x_i) \;=\; 2\,E_2(\tau_{\rm out}(x_i)).

    Returned as the SymPy expression :math:`2\,E_2(\tau)`; same
    :math:`t = 1/\mu` substitution as :func:`derive_slab_p_outer_n0`,
    but with the no-:math:`\mu`-weight integrand (factor of 2 surface
    normalisation absorbed in the 2 prefactor).
    """
    tau = sp.Symbol("tau", positive=True)
    return 2 * sp.expint(2, tau)


def build_R_slab_blockdiag(n_modes: int) -> sp.Matrix:
    r"""Symbolic block-diagonal slab specular reflection operator.

    Returns

    .. math::

       R_{\rm slab}(N) \;=\; \begin{pmatrix}
           R_{\rm spec}(N) & 0 \\
           0 & R_{\rm spec}(N)
       \end{pmatrix}_{2N \times 2N},

    a :math:`2N \times 2N` block-diagonal matrix where each
    :math:`N \times N` block is the curvilinear specular operator
    :math:`R_{\rm spec}(N) = (1/2)\,M^{-1}` from
    :func:`.r_matrix.build_R_specular_symbolic`.

    The block-diagonal structure follows from the locality of specular
    reflection at each face — a particle leaving the outer face cannot
    appear at the inner face (except through :math:`K_{\rm vol}`).
    """
    R_spec = build_R_specular_symbolic(n_modes)
    R_slab = sp.zeros(2 * n_modes, 2 * n_modes)
    R_slab[:n_modes, :n_modes] = R_spec
    R_slab[n_modes:, n_modes:] = R_spec
    return R_slab


def slab_p_outer_mode_n_no_mu(n_modes: int) -> sp.Expr:
    r"""Symbolic :eq:`peierls-slab-spec-Pesc-no-mu-weight` —
    no-:math:`\mu`-weight basis P primitive at mode :math:`n`:

    .. math::

       P_{\rm esc, face}^{(n)}(x_i)
           \;=\; \tfrac{1}{2}\sum_{k=0}^{n} c_n^k\,E_{k+2}(\tau).

    Used internally by the symbolic test suite and as the canonical
    formula for any future :math:`n \ge 1` implementer.
    """
    tau = sp.Symbol("tau", positive=True)
    coefs = shifted_legendre_monomial_coefs(n_modes)
    return sp.Rational(1, 2) * sum(
        sp.Rational(0).__class__(c)  # keep float fidelity
        if False else c * sp.expint(k + 2, tau)
        for k, c in enumerate(coefs)
    )


def slab_g_outer_mode_n_no_mu(n_modes: int) -> sp.Expr:
    r"""Symbolic :eq:`peierls-slab-spec-Gbc-no-mu-weight` —
    no-:math:`\mu`-weight basis G primitive at mode :math:`n`:

    .. math::

       G_{\rm bc, face}^{(n)}(x_i)
           \;=\; 2\sum_{k=0}^{n} c_n^k\,E_{k+2}(\tau).
    """
    tau = sp.Symbol("tau", positive=True)
    coefs = shifted_legendre_monomial_coefs(n_modes)
    return 2 * sum(c * sp.expint(k + 2, tau) for k, c in enumerate(coefs))
