r"""An infinite medium's flux is isotropic, so :math:`P_L` order cannot move it.

⛔ **This module is RED by design at** :math:`L \ge 2` **(ERR-080).** The
three strict-xfail rows below (`[M]` by AST 2026-09-02 — this line read
"two" while the third row existed) assert the physics; the tree does not
yet deliver it. They retire themselves — an ``XPASS`` under
``xfail(strict=True)`` is a FAILURE, so the marker cannot outlive the
repair.

The claim, in the domain's own terms
------------------------------------

In an infinite homogeneous medium with an isotropic source the angular
flux is isotropic. Every Legendre moment :math:`\phi_\ell` with
:math:`\ell \ge 1` is therefore **exactly zero**, and the anisotropic
scattering source

.. math::

   Q_s(\Omega) \;=\; \sum_{\ell=0}^{L} \frac{2\ell+1}{4\pi}\,
                     \Sigma_{s,\ell}\, \phi_\ell\, P_\ell(\Omega\cdot\Omega')

collapses to its :math:`\ell = 0` term whatever :math:`L` and whatever
the higher :math:`\Sigma_{s,\ell}` are. So the converged scalar flux is

.. math::

   \phi \;=\; \frac{W\,Q}{\Sigma_t - \Sigma_{s,0}}

at **every** scattering order. ``scattering_order`` is a knob on the
*representation*, not on the *answer*: that is the whole content of this
module, and it is why the higher moments here are given deliberately
large, ugly values (:math:`\Sigma_{s,1..3}` comparable to
:math:`\Sigma_{s,0}`). A row that only passes for small higher moments
would be measuring smallness, not isotropy.

What actually happens (ERR-080)
-------------------------------

``[M]`` 2026-08-31, found while building #426's missing (n,2n)
measurement:

======  ===================  ==========
``L``   converged ``phi``    rel. error
======  ===================  ==========
0       ``+4.000000000000``  1.2e-15
1       ``+4.000000000000``  4.4e-16
2       ``-3.764705882353``  1.94
3       ``+0.454164575725``  0.886
======  ===================  ==========

The :math:`L = 2` answer has the **wrong sign** and a 194 % error, from a
public, unguarded entry point.

⚠ Read the :math:`L = 3` row with care: it is *positive*, which makes it
look like a partial recovery. It is not — it is 89 % low, and the sign is
an accident of which degenerate harmonics happen to cancel at that order.
There is no monotone trend here to extrapolate, and no order at which the
answer becomes trustworthy; that is what makes a single "the higher orders
drift" summary the wrong reading of this table.

Mechanism, in one sentence: a 1-D quadrature supplies
:math:`\mu_y = \mu_z = 0` meaning *"there is no azimuthal information"*,
and ``_evaluate_real_sh`` reads it as *"the azimuth is 0"*
(``arctan2(0, 0) = 0``), so every :math:`m > 0` harmonic becomes a
non-zero constant across the ordinate set instead of being absent. The
degenerate harmonics are linearly dependent on the slab nodes, the
discrete Gram is rank-deficient, and the per-mode scattering multiplier
is no longer a function of the flux.

:math:`\ell = 1` is clean **only because it never takes that path** —
``_evaluate_real_sh`` hard-codes :math:`\ell = 1` in Cartesian form
(:math:`Y_1 \propto (\mu_z, \mu_x, \mu_y)`), where the transverse zeros
are genuinely zero. That accident is what makes the defect read as an
":math:`L \ge 2` problem" rather than as what it is, and it is why the
:math:`L \le 1` rows below are **controls, not decoration**: they prove
this fixture reaches the analytic answer, so an xfail two rows down
cannot be passing because the fixture is broken
(``vv-principles`` #17 — a battery without a positive control measures
the harness).

Why the suite never caught it
-----------------------------

``[M]`` the corpus uses ``scattering_order`` 0 (74 sites), 1 (46 sites)
and 3 (2 sites) — and **both** :math:`L = 3` sites are
``Quadrature.lebedev(17)`` on a 2-D mesh, a 3-D rule where the
fabrication does not occur. :math:`P_{\ge 2}` on a 1-D chart had no gate
anywhere in the tree. This module is that gate.

The repair
----------

Not a special case: a 1-D angular quadrature is a quadrature on the
orbit space :math:`S^2/SO(2)`, and the harmonics that survive are that
quotient's **trivial isotypic component** :math:`\{Y_\ell^0\} \cong
\{P_\ell\}`. Tracked by **#429**, planned in
``.claude/plans/angular_spaces_derived_from_symmetry.md`` (Phase 3.4 is
the fix; Phase 0.7 is this file).

Structure
---------

The four solves run in a **module-scoped fixture** on purpose. A solve
that raises is then a setup ERROR rather than an xfail, so a strict-xfail
row can only ever fail for the one documented reason — its single
assertion on the converged value.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import BC
from orpheus.geometry.mesh import Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn_fixed_source

pytestmark = pytest.mark.l1

#: Total cross-section, and the Legendre scattering moments Sigma_{s,0..3}.
#: The higher moments are deliberately LARGE (comparable to Sigma_{s,0}):
#: isotropy of the converged flux must hold identically, not asymptotically.
_SIG_T: float = 1.0
_MOMENTS: tuple[float, ...] = (0.5, 0.40, 0.45, 0.20)

_N_ORDINATES: int = 8
_N_CELLS: int = 4
_SLAB_WIDTH: float = 4.0
_SOURCE_PER_ORDINATE: float = 1.0

#: Relative tolerance on the converged flux. The inner solve is driven to
#: 1e-13, so a correct answer lands at ~1e-15; 1e-10 is ~5 orders of
#: headroom and still 9 orders tighter than the 1.9e+00 defect.
_RTOL: float = 1e-10


def _one_group_mixture(n_moments: int) -> Mixture:
    """A 1-group medium carrying exactly ``n_moments`` Legendre moments."""
    zero = np.zeros(1)
    return Mixture(
        SigC=np.array([_SIG_T - _MOMENTS[0]]),
        SigL=zero.copy(),
        SigF=zero.copy(),
        SigP=zero.copy(),
        SigT=np.array([_SIG_T]),
        SigS=[csr_matrix(np.array([[v]])) for v in _MOMENTS[:n_moments]],
        Sig2=csr_matrix((1, 1)),
        chi=zero.copy(),
        eg=None,
    )


@pytest.fixture(scope="module")
def infinite_medium_flux() -> dict[int, float]:
    r"""Converged scalar flux at each scattering order, and the reference.

    Reflective on both faces with a uniform per-ordinate source: the
    problem has no spatial and no angular structure, so the mesh and the
    quadrature only supply the total weight :math:`W`.

    Module-scoped and eager (all four orders solved once) so that any
    solve failure surfaces as a setup ERROR for every row, never as a
    strict-xfail passing for an undocumented reason.
    """
    quadrature = Quadrature.gauss_legendre(n_ordinates=_N_ORDINATES)
    mesh = Mesh1D(
        edges=np.linspace(0.0, _SLAB_WIDTH, _N_CELLS + 1),
        mat_ids=np.zeros(_N_CELLS, dtype=int),
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    source = np.full((_N_ORDINATES, 1, _N_CELLS), _SOURCE_PER_ORDINATE)

    fluxes: dict[int, float] = {}
    for order in range(len(_MOMENTS)):
        solution = solve_sn_fixed_source(
            materials={0: _one_group_mixture(order + 1)},
            mesh=mesh,
            quadrature=quadrature,
            external_source=source,
            scattering_order=order,
            inner_tol=1e-13,
            max_inner=5000,
            inner_solver="krylov",
        )
        fluxes[order] = float(np.asarray(solution.scalar_flux.values)[0, 0])

    total_weight = float(np.asarray(quadrature.weights).sum())
    fluxes["reference"] = (  # type: ignore[index]
        total_weight * _SOURCE_PER_ORDINATE / (_SIG_T - _MOMENTS[0])
    )
    return fluxes


# ─── controls: these PASS today, and they are what licenses the xfails ───


@pytest.mark.parametrize("order", [0, 1])
def test_isotropic_flux_is_analytic_at_low_order(
    infinite_medium_flux: dict[int, float], order: int
) -> None:
    r"""``L = 0, 1`` reach :math:`WQ/(\Sigma_t - \Sigma_{s,0})`.

    The positive control for the two rows below. :math:`\ell = 1` is
    clean only because ``_evaluate_real_sh`` evaluates it in Cartesian
    form and never fabricates an azimuth — so this row passing is
    evidence the FIXTURE is sound, not evidence the machinery is.
    """
    np.testing.assert_allclose(
        infinite_medium_flux[order],
        infinite_medium_flux["reference"],  # type: ignore[index]
        rtol=_RTOL,
        err_msg=(
            f"scattering_order={order} must not move the infinite-medium "
            f"flux; this row is the control for the ERR-080 xfails"
        ),
    )


# ─── the defect, one assertion per row, strict so it self-retires ───


@pytest.mark.catches("ERR-080")
@pytest.mark.xfail(
    strict=True,
    reason=(
        "ERR-080 / #429: a 1-D quadrature's mu_y = mu_z = 0 means 'no "
        "azimuthal information', but _evaluate_real_sh reads it as "
        "'azimuth 0' (arctan2(0,0)), so every m > 0 harmonic degenerates "
        "to a non-zero constant and the l >= 2 scattering source is "
        "computed from moments that are not moments. Measured: "
        "phi = -3.764705882353 against an analytic +4.000000000000. "
        "Repaired by Phase 3.4 of "
        ".claude/plans/angular_spaces_derived_from_symmetry.md; this "
        "marker MUST be deleted then (XPASS is a failure under strict)."
    ),
)
def test_isotropic_flux_is_analytic_at_p2(
    infinite_medium_flux: dict[int, float]
) -> None:
    r"""``L = 2`` must reach the same flux. ``[M]`` it returns ``-3.7647``."""
    np.testing.assert_allclose(
        infinite_medium_flux[2],
        infinite_medium_flux["reference"],  # type: ignore[index]
        rtol=_RTOL,
        err_msg="scattering_order=2 moved the infinite-medium flux (ERR-080)",
    )


@pytest.mark.catches("ERR-080")
@pytest.mark.xfail(
    strict=True,
    reason=(
        "ERR-080 / #429, same mechanism as the P2 row one order up — "
        "recorded separately because the exit gate for #429 requires the "
        "analytic value at L = 0, 1, 2 AND 3, and a single row cannot "
        "witness an order it does not run."
    ),
)
def test_isotropic_flux_is_analytic_at_p3(
    infinite_medium_flux: dict[int, float]
) -> None:
    r"""``L = 3`` must reach the same flux."""
    np.testing.assert_allclose(
        infinite_medium_flux[3],
        infinite_medium_flux["reference"],  # type: ignore[index]
        rtol=_RTOL,
        err_msg="scattering_order=3 moved the infinite-medium flux (ERR-080)",
    )


@pytest.mark.catches("ERR-080")
@pytest.mark.xfail(
    strict=True,
    reason=(
        "ERR-080 / #429, the SECOND public symptom: at higher (order, L) on "
        "a 1-D chart the fabricated harmonics do not merely corrupt the "
        "answer, they make the discrete Gram inconsistent enough that "
        "DenseMetric's Penrose guard REFUSES it -- so the public entry "
        "point RAISES. Measured: solve_sn_fixed_source(gauss_legendre(16), "
        "scattering_order=4) -> ValueError 'DenseMetric was handed an "
        "inconsistent inverse face: max|G G+ G - G| = 7.117e-04'. Note the "
        "message blames the metric, three hops downstream of the "
        "fabrication. Retires with the P2/P3 rows."
    ),
)
def test_high_order_anisotropic_slab_does_not_raise() -> None:
    r"""``GL16`` at :math:`L = 4` must solve, not raise.

    A SEPARATE failure mode from the two rows above, and it needs its own
    witness: those return a wrong number, this one refuses to return at
    all. ``vv-principles`` #17's granularity trap — one gate covering both
    would report ``xfail`` from whichever arm it happened to reach, and
    certify the other.

    The solve runs INSIDE the test (not in the module fixture) precisely
    because it currently raises: in the fixture that would be a setup
    ERROR for every row in the module rather than an ``xfail`` here.

    ``[M]`` 2026-08-31 the 3-D rules are unaffected at every ``L``, so this
    is a property of the 1-D chart, not of the order.
    """
    quadrature = Quadrature.gauss_legendre(n_ordinates=16)
    mesh = Mesh1D(
        edges=np.linspace(0.0, _SLAB_WIDTH, _N_CELLS + 1),
        mat_ids=np.zeros(_N_CELLS, dtype=int),
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    solution = solve_sn_fixed_source(
        materials={0: _one_group_mixture(len(_MOMENTS))},
        mesh=mesh,
        quadrature=quadrature,
        external_source=np.full((16, 1, _N_CELLS), _SOURCE_PER_ORDINATE),
        scattering_order=len(_MOMENTS) - 1,
        inner_tol=1e-13,
        max_inner=5000,
        inner_solver="krylov",
    )
    total_weight = float(np.asarray(quadrature.weights).sum())
    np.testing.assert_allclose(
        float(np.asarray(solution.scalar_flux.values)[0, 0]),
        total_weight * _SOURCE_PER_ORDINATE / (_SIG_T - _MOMENTS[0]),
        rtol=_RTOL,
        err_msg="GL16 at L=4 moved the infinite-medium flux (ERR-080)",
    )
