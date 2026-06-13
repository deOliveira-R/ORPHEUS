r"""Foundation: the curvilinear coupled-pole seed's mirror is the
intra-level sign-flip partner (Issue #193, re-scoped).

Why this test exists — the re-scope
-----------------------------------

#193 originally asked for a "BC-realizer level-locality" invariant: the
cylindrical matvec was said to build its per-level ``CarlsonSweepContext``
by calling ``bc_outer.apply`` once on the full ψ then extracting the
most-inward ordinate per level — a pattern that would silently break if a
future BC coupled ordinates across μ-levels.

**That premise dissolved.** Wave O O.4a.2 deleted the ``bc_outer.apply``
keystone from both matvec twins (the reflective coupling ``B`` moved OUTSIDE
the bare sweep as a first-class sibling operator), and the surviving SI-sweep
seed builder reads the *raw* given inflow trace (``inflow_full = bc_outer``,
no apply) — so the apply/restrict commutativity #193 targeted is now vacuous
(a pure gather always commutes with restriction).

But a **narrower, genuinely load-bearing** μ-level invariant moved into the
ERR-058 coupled-pole spatial seed (Issue #195, this codebase's curvilinear
fix). Both matvec twins (``orpheus/sn/operator.py`` ``_compute_LpC`` /
``_compute_decomposition``) and the SI sweep
(``orpheus/sn/loss_representation.py``) seed the +μ (outward) pole face from
the −μ (inward) sweep's pole outflow at the **mirror ordinate**::

    pole_face_seed = outflow_at_inner.T[quad.reflection_index("x")]

This realises the Carlson coupled-pole continuity
:math:`\psi(0,+\mu_r) = \psi(0,-\mu_r)` (:ref:`the documented continuity
<sn-err-058-coupled-pole-continuity>`).  It is correct **iff**
``reflection_index("x")[n]`` is the **intra-level sign-flip partner** of
ordinate ``n``: the same μ-level (same axial cosine :math:`\mu = \mu_z`), with
:math:`\mu_x` negated and the other cosines held.  The physics demands it —
the pole continuity holds at *fixed axial direction*, so the reflected partner
must stay in the same level.

If ``reflection_index("x")`` ever returned a **cross-level** partner (a future
cubature, or a refactor of ``_compute_sphere_reflection_partners``), the
``pole_outflow[mirror[n]]`` read would grab a *different axial direction's*
value — a **silent** correctness break: on a flat ψ field the mirror value
equals the cell value, so every flat-flux gate stays green (the exact ERR-058
blindness class — vv-principles Mode 7, at the operator-internals level).  This
gate pins the invariant so a future change cannot break it silently.

The invariant holds by construction today: ``reflection_index("x")`` is
``_find_reflections(-mu_x, mu_y, mu_z, ...)`` (``directional.py``), which flips
only :math:`\mu_x` and holds :math:`\mu_y, \mu_z`.  The cylinder μ-level is keyed
on the **axial** cosine :math:`\mu = \mu_z` (``rules_sphere.py`` groups by
``|mu_z|``; ``rules_product.py`` levels are fixed-``mu_z`` GL nodes), and the
x-mirror holds :math:`\mu_z` — so the x-mirror provably preserves the level.

Why this is ``foundation``, *not* a ``verifies(...)`` of the continuity label.
The physics continuity :math:`\psi(0,+\mu)=\psi(0,-\mu)` is a representational
identity (``:vv-status: documented``), whose *verifiable* content lives in the
curvilinear-MMS operator-admission gate (``test_curvilinear_operator_admits_mms``,
``catches("ERR-058")``) and the seed-adjoint bit-identity gate — NOT in this
reflection-index data-structure contract.  This test pins the **software
invariant** the seed relies on; attaching the physics ``:label:`` here would
write a misleading ``tests`` edge promoting the equation to "covered" on the
strength of a structure check that never touches a non-flat ψ profile.  Per the
V&V harness contract (``tests/_harness/registry.py``): foundation tests never
carry ``verifies(...)``.  The equation is referenced above by ``:ref:`` for the
human reader only.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.quadrature import Quadrature

# Curvilinear-relevant cubatures.  Sphere/slab GL1D is single-level (the
# intra-level leg is trivial; it still exercises the μ_x sign-flip leg).
# Level-symmetric and product cubatures are the genuine multi-level cases the
# cylinder uses — LS4 has 2 levels (16+8), product(4,8) has 4 levels of 8.
_CUBATURES = [
    pytest.param(lambda: Quadrature.gauss_legendre(8), id="gl8"),
    pytest.param(lambda: Quadrature.gauss_legendre(16), id="gl16"),
    pytest.param(lambda: Quadrature.level_symmetric(4), id="ls4"),
    pytest.param(lambda: Quadrature.level_symmetric(6), id="ls6"),
    pytest.param(lambda: Quadrature.product(n_mu=4, n_phi=8), id="product_4x8"),
    pytest.param(lambda: Quadrature.product(n_mu=8, n_phi=8), id="product_8x8"),
]


def _level_of(quad: Quadrature) -> np.ndarray:
    """``(N,)`` array mapping each ordinate index to its μ-level index."""
    level_of = np.empty(quad.N, dtype=int)
    for p, ids in enumerate(quad.level_indices):
        level_of[np.asarray(ids, dtype=int)] = p
    return level_of


@pytest.mark.foundation
@pytest.mark.parametrize("quad_factory", _CUBATURES)
def test_x_reflection_is_intra_level_signflip_partner(quad_factory):
    r"""``reflection_index("x")`` is the intra-level :math:`\mu_x`-sign-flip
    partner — the contract the ERR-058 coupled-pole seed relies on.

    Three facets, asserted together for every ordinate ``n`` with
    ``m = reflection_index("x")[n]``:

    1. **Intra-level** — ``level(m) == level(n)``.  The coupled-pole partner
       stays in ``n``'s μ-level, so ``pole_outflow[m]`` is the *same axial
       direction's* already-computed inward outflow, never a cross-level /
       uninitialised read.
    2. **μ_x sign-flip** — ``mu_x[m] == -mu_x[n]``.
    3. **Other cosines held** — ``mu_y[m] == mu_y[n]``, ``mu_z[m] == mu_z[n]``
       — so it is specifically the *x*-reflection, not some other same-level
       direction.

    Assertions go through ``np.testing`` (function calls fire under
    ``python -O``; bare ``assert`` would be stripped — vv-principles Mode 8).
    """
    quad = quad_factory()
    mirror = np.asarray(quad.reflection_index("x"), dtype=int)

    # (2)+(3): the full node maps to itself with column 0 (μ_x) negated.
    mu = np.column_stack(
        [quad.axis_cosines(0), quad.axis_cosines(1), quad.axis_cosines(2)]
    )
    expected = mu.copy()
    expected[:, 0] *= -1.0  # -0.0 for axis-tangent μ_x=0 ordinates → still ==0
    np.testing.assert_allclose(
        mu[mirror], expected, atol=1e-12,
        err_msg=(
            "reflection_index('x') is not the μ_x sign-flip partner "
            "(μ_y/μ_z must be held) — the coupled-pole seed would read the "
            "wrong direction's pole outflow (ERR-058 class, Issue #193)."
        ),
    )

    # (1): the partner stays in the same μ-level.
    level_of = _level_of(quad)
    np.testing.assert_array_equal(
        level_of[mirror], level_of,
        err_msg=(
            "reflection_index('x') maps an ordinate to a DIFFERENT μ-level — "
            "the coupled-pole seed's pole_outflow[mirror[n]] read would grab a "
            "cross-level (wrong axial direction) value, a SILENT break "
            "invisible to flat-flux gates (ERR-058 class, Issue #193)."
        ),
    )


@pytest.mark.foundation
@pytest.mark.parametrize("quad_factory", _CUBATURES)
def test_x_reflection_is_an_involution(quad_factory):
    r"""``reflection_index("x")`` is its own inverse — applying it twice
    returns the identity.

    The coupled-pole seed pairs ``(+μ_x, -μ_x)`` ordinates; a non-involutive
    map would pair an ordinate with a partner that does not pair back, leaving
    one direction of each pole couple un-seeded.  Cheap structural companion to
    the sign-flip/level invariant above.
    """
    quad = quad_factory()
    mirror = np.asarray(quad.reflection_index("x"), dtype=int)
    np.testing.assert_array_equal(
        mirror[mirror], np.arange(quad.N),
        err_msg="reflection_index('x') is not an involution (mirror∘mirror ≠ id).",
    )
