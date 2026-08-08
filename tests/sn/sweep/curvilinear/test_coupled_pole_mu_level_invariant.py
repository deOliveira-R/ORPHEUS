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

    pole_face_seed = outflow_at_inner.T[self._ensure_pole_mirror()]

This realises the Carlson coupled-pole continuity
:math:`\psi(0,+\mu_r) = \psi(0,-\mu_r)` (:ref:`the documented continuity
<sn-err-058-coupled-pole-continuity>`).  It is correct **iff** the derived
mirror pairing sends ``n`` to its **intra-level sign-flip partner**:
the same μ-level (same axial cosine :math:`\mu = \mu_z`), with
:math:`\mu_x` negated and the other cosines held.  The physics demands it —
the pole continuity holds at *fixed axial direction*, so the reflected partner
must stay in the same level.

If the derived pairing ever returned a **cross-level** partner (a future
cubature, or a refactor of ``ordinate_permutation``'s match machinery), the
``pole_outflow[mirror[n]]`` read would grab a *different axial direction's*
value — a **silent** correctness break: on a flat ψ field the mirror value
equals the cell value, so every flat-flux gate stays green (the exact ERR-058
blindness class — vv-principles Mode 7, at the operator-internals level).  This
gate pins the invariant so a future change cannot break it silently.

The invariant holds by construction today: the seed derives its pairing at
first use (``_OneDimScanWalk._ensure_pole_mirror``) from the σ_x mirror
MOTION via ``Quadrature.ordinate_permutation`` — certified bijective and
weight-preserving, with the same match discipline the retired
``reflection_index`` table carried until G6.3 §7d — and the x-mirror's
action flips only :math:`\mu_x` and holds
:math:`\mu_y, \mu_z`.  The cylinder μ-level is keyed
on the **axial** cosine :math:`\mu = \mu_z` (``rules_sphere.py`` groups by
``|mu_z|``; ``rules_product.py`` levels are fixed-``mu_z`` GL nodes), and the
x-mirror holds :math:`\mu_z` — so the x-mirror provably preserves the level.

Why this is ``foundation``, *not* a ``verifies(...)`` of the continuity label.
The physics continuity :math:`\psi(0,+\mu)=\psi(0,-\mu)` is a representational
identity (``:vv-status: documented``), whose *verifiable* content lives in the
curvilinear-MMS operator-admission gate (``test_curvilinear_operator_admits_mms``,
``catches("ERR-058")``) and the seed-adjoint bit-identity gate — NOT in this
mirror-pairing data-structure contract.  This test pins the **software
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

from orpheus.derivations.continuous.mms.sn import _make_1g_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.geometry.boundary import SelfPairedDeck
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import solve_sn_fixed_source


def _seed_mirror(quad: Quadrature) -> np.ndarray:
    """The seed's x-mirror pairing, derived exactly as production derives it.

    This is deliberately the SAME expression
    ``_OneDimScanWalk._ensure_pole_mirror`` evaluates — the gates below pin
    properties of the datum production actually consumes, not of a
    reference lookalike (a lookalike could hold the invariant while the
    production derivation drifts).
    """
    pi = quad.ordinate_permutation(SelfPairedDeck.mirror(axis="x").motion)
    if pi is None:
        raise AssertionError(
            "no x-mirror ordinate pairing on this rule — every _CUBATURES "
            "row is expected to be σ_x-closed"
        )
    return np.asarray(pi.indices, dtype=int)

# Curvilinear-relevant cubatures.  Sphere/slab GL1D is single-level (the
# intra-level leg is trivial; it still exercises the μ_x sign-flip leg).
# The FOLDED rows are the family SNMesh(CYLINDRICAL) admits (Q5.6.3 —
# the production pole seed derives its pairing on these); the
# level-symmetric and full-product rows stay as QUADRATURE-level
# multi-level instances of the σ_x-closure contract (these foundation
# tests never build an SNMesh, and the rules remain constructible —
# only the cylindrical admission refuses them).  The fold quotients by
# σ_y, so σ_x-closure survives it by construction; these gates measure
# that it does.
_CUBATURES = [
    pytest.param(lambda: Quadrature.gauss_legendre(8), id="gl8"),
    pytest.param(lambda: Quadrature.gauss_legendre(16), id="gl16"),
    pytest.param(lambda: Quadrature.level_symmetric(4), id="ls4"),
    pytest.param(lambda: Quadrature.level_symmetric(6), id="ls6"),
    pytest.param(lambda: Quadrature.product(n_mu=4, n_phi=8), id="product_4x8"),
    pytest.param(lambda: Quadrature.product(n_mu=8, n_phi=8), id="product_8x8"),
    pytest.param(
        lambda: Quadrature.folded_product(n_mu=4, n_phi=8), id="folded_4x8",
    ),
    pytest.param(
        lambda: Quadrature.folded_product(n_mu=8, n_phi=8), id="folded_8x8",
    ),
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
    r"""The seed's x-mirror pairing is the intra-level :math:`\mu_x`-sign-flip
    partner — the contract the ERR-058 coupled-pole seed relies on.

    Three facets, asserted together for every ordinate ``n`` with
    ``m = mirror[n]`` (the pairing derived exactly as production derives it):

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
    mirror = _seed_mirror(quad)

    # (2)+(3): the full node maps to itself with column 0 (μ_x) negated.
    mu = np.column_stack(
        [quad.axis_cosines(0), quad.axis_cosines(1), quad.axis_cosines(2)]
    )
    expected = mu.copy()
    expected[:, 0] *= -1.0  # -0.0 for axis-tangent μ_x=0 ordinates → still ==0
    np.testing.assert_allclose(
        mu[mirror], expected, atol=1e-12,
        err_msg=(
            "the seed's x-mirror pairing is not the μ_x sign-flip partner "
            "(μ_y/μ_z must be held) — the coupled-pole seed would read the "
            "wrong direction's pole outflow (ERR-058 class, Issue #193)."
        ),
    )

    # (1): the partner stays in the same μ-level.
    level_of = _level_of(quad)
    np.testing.assert_array_equal(
        level_of[mirror], level_of,
        err_msg=(
            "the seed's x-mirror pairing maps an ordinate to a DIFFERENT "
            "μ-level — the coupled-pole seed's pole_outflow[mirror[n]] read "
            "would grab a cross-level (wrong axial direction) value, a SILENT "
            "break invisible to flat-flux gates (ERR-058 class, Issue #193)."
        ),
    )


@pytest.mark.foundation
@pytest.mark.parametrize("quad_factory", _CUBATURES)
def test_x_reflection_is_an_involution(quad_factory):
    r"""The seed's x-mirror pairing is its own inverse — applying it twice
    returns the identity.

    The coupled-pole seed pairs ``(+μ_x, -μ_x)`` ordinates; a non-involutive
    map would pair an ordinate with a partner that does not pair back, leaving
    one direction of each pole couple un-seeded.  Cheap structural companion to
    the sign-flip/level invariant above.
    """
    quad = quad_factory()
    mirror = _seed_mirror(quad)
    np.testing.assert_array_equal(
        mirror[mirror], np.arange(quad.N),
        err_msg="the seed's x-mirror pairing is not an involution (mirror∘mirror ≠ id).",
    )


@pytest.mark.foundation
def test_a_rule_without_the_x_mirror_refuses_before_any_march():
    r"""A σ_x-unclosed rule refuses the curvilinear sweep LOUDLY, with the
    pairing's own diagnosis — the pole seed never runs on a wrong partner.

    ``product(4, 5)`` (odd :math:`n_\varphi`) has no x-mirror closure — the
    ERR-074 fixture: its pre-certification partner table was wrong by 0.58
    in the direction cosines *and still involutive*.  Until G6.3 §7d the
    only barrier here was the retired ``reflection_index`` table's generic
    lookup ``ValueError``, raised MID-ITERATION; since 7d.2 the seed derives
    its pairing at first use (``_ensure_pole_mirror``) and refuses with the
    coupled-pole continuation's own message, before any march.

    Pinned on the message FRAGMENT, not the exception type: this is the
    re-point of the retired ``"no precomputed reflection partner"`` pin, and
    ``BoundaryError ⊂ ValueError`` taught the campaign that a type-only
    refusal pin stays green across the wrong refusal.
    """
    nx = 4
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),          # r = 0 symmetry axis
        bc_right=BC("vacuum"),             # no mirror demand from the BC tier
    )
    quad = Quadrature.product(n_mu=4, n_phi=5)
    Q = np.full((quad.N, 1, nx), 1.0)
    with pytest.raises(ValueError, match="cannot seed the r = 0 pole"):
        solve_sn_fixed_source(
            {0: _make_1g_mixture(1.0, 0.5)}, mesh, quad, Q,
            max_inner=2, inner_tol=1e-13,
        )
