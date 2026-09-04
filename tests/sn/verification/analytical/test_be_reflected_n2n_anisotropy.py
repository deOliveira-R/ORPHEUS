r"""The #426 flagship: a Be-reflected fast slab reads the (n,2n) anisotropy
the tape carries.

The EIGENVALUE claim layer on the production library: it reads the
**tracked** ``.GXS``-derived library through ``load_isotope`` and builds its
own mixtures; it opens no path under ``scratch/`` (issue #444 / the
2026-09-03 plan-authoring surprise: an untracked artefact a gate reads is
green on exactly one checkout).

────────────────────────────────────────────────────────────────────────────
CLAIM LEDGER
────────────────────────────────────────────────────────────────────────────
========================  ==============  ===========  =========================
row                       claim layer     kind         truth source
========================  ==============  ===========  =========================
``…_p0_control_…``        eigenvalue      RECORD       frozen pre-carve k
``…_reads_the_measured…`` eigenvalue      RECORD+band  §0 measurement, qa-reproduced
``…_the_delta_is_…``      eigenvalue      REFERENCE*   forward-peaked emission (physics)
``…_ladder_has_converged`` eigenvalue     REFERENCE*   the ℓ-ladder's measured decay
``…_p0_solve_is_…``       eigenvalue      RECORD       frozen pre-carve k at L=0
========================  ==============  ===========  =========================

\* the Δ row is the only one whose EXPECTED SIGN comes from physics rather than
from a recording: (n,2n) emission in the reflector is forward-peaked
(`[M]` μ̄ = rowsum(Σ₁)/rowsum(Σ₀) = **+0.2524 … +0.4348**, mean +0.3031, positive
on **50 of 50** live Be-9 groups), so the emitted pair continues OUTWARD and
less returns to the fuel.  The P0 model emitted isotropically and therefore
OVER-estimated the reflector's return current ⟹ restoring the peak must LOWER
k.  A magnitude band is attached because there is no structurally-independent
eigenvalue reference for a 421-group heterogeneous reflected slab (lessons §4:
ORPHEUS has none) — the band is calibrated from the §0 measurement,
independently reproduced by `qa` (``scratch/_426_repro.md``: three k values to
every published digit, C0 exactly 0.0, six controls).

────────────────────────────────────────────────────────────────────────────
§6c — the RED-BEFORE, measured 2026-09-03 on the P0-model tree
────────────────────────────────────────────────────────────────────────────
=========================  =====================  ==========================
arm                        k under the P0 model   the carve's target
=========================  =====================  ==========================
``scattering_order=2``     **1.0953221881419453** 1.091199657  (Δ −412.25·1e-5)
``scattering_order=0``     **1.1587120371368607** unchanged (bit-identity leg)
Δ (full − ℓ≥1 zeroed)      **0.0** (no ℓ≥1 to     −4.12e-3
                           zero: the arms coincide)
=========================  =====================  ==========================

⟹ the value row was red by 4.14e-3 = **41× its own band**, and the Δ row red
at exactly 0.0, until #426 step 2 (2026-09-04) made the (n,2n) binding read
the tape's stack at the solve's order.

────────────────────────────────────────────────────────────────────────────
COST, and the ``slow`` ruling
────────────────────────────────────────────────────────────────────────────
`[M]` one k-solve on this fixture = **6.84 s** at ``scattering_order=2`` and
4.42 s at 0; the mixture build is 0.29 s (module-scoped, once).  Four solves
≈ **25 s**.

⛔ **NOT ``slow``, and it must not become so.**  #428 F-4 is the measured
precedent: MC's only ERR-023 catcher carried ``@pytest.mark.slow``, so under
the project's canonical ``-m "not slow"`` invocation the MC tree read *39
passed / 0 red* with the multiplicity mutated — a real catcher the gate that
matters could not see.  25 s inside a ≥90-min gate is 0.5 %.
The FAST companions, so a future session cannot demote this file and lose the
whole claim: the ℓ ≥ 1 (n,2n) term is ALSO pinned at the term tier by
``tests/transport/test_material_field.py`` / ``test_transfer_kernel.py`` (ms)
and ``tests/sn/operators/test_n2n_operator.py`` (the activation and the
yield law), and the ν₂ₙ VALUE end-to-end by the ``2eg_n2n`` rows of
``test_kinf_homogeneous.py`` (`[M]` +12 rows, +2.4 s, both functionals redden
20.2 % / 10.5 % under ν₂ₙ 2→1).

────────────────────────────────────────────────────────────────────────────
CONFIG-BLINDNESS DECLARATION (`AGENT.md` §0.6, run row by row)
────────────────────────────────────────────────────────────────────────────
* NOT flat flux — 3-region heterogeneous, vacuum both faces.
* NOT 1-group — **421** groups.
* NOT homogeneous — Be | U-235 metal | Be.
* NOT slab-degenerate for the term under test — the (n,2n) ℓ=1 source is a
  *streaming*-coupled anisotropy, fully active on a slab (unlike the
  curvilinear redistribution terms, which a slab nulls).
* The regime ACTIVATES the term: `[M]` 99.93 % of the effect is the Be-9
  reflector (Be-only −377.30 vs U-only −0.26 Δk/k₀·1e5, sum = A1 exactly),
  and Be-9's MT=16 ℓ=1 is 13× more forward-peaked than U-235's.
* ⛔ WHAT THIS FIXTURE IS BLIND TO, stated so no one over-reads it: every arm
  runs at ``scattering_order = 2``, so the ELASTIC channel is P2 in every arm
  and `[M]` the elastic P3…P6 the ingest now carries move this k by a further
  −229 / −163 / −175 / −173 Δk·1e5 (step 1's exit measurement,
  ``test_scattering_order_is_the_only_truncation.py``).  So the number here
  is the (n,2n) anisotropy ALONE at P2, never "the converged anisotropy
  answer" (`qa` finding D5).
────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations

import contextlib
import io
import warnings
from dataclasses import replace

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import compute_macro_xs
from orpheus.data.micro_xs import load_isotope
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn

pytestmark = [pytest.mark.l2, pytest.mark.catches("ERR-082")]

# ── the §0 fixture, verbatim ──────────────────────────────────────────
_N_U_METAL = 0.04894  # /b·cm, 19.1 g/cc
_N_BE = 0.1236  # /b·cm, 1.85 g/cc
_TEMP_K = 294
#: The §0 solve configuration, spelled at the ONE call site below rather than
#: splatted from a dict.  `[M]` a ``**dict`` splat costs 4 pyright errors here
#: (the checker cannot see the keys and binds them to ``mat_map``) — the
#: ``# type: ignore`` reflex would hide a real one (lessons `L44k`).
_KEFF_TOL = 1e-9
_FLUX_TOL = 1e-8
_INNER_TOL = 1e-10
_MAX_OUTER = 3000

#: `[M]` 2026-09-03 pre-carve, this file's own configuration.  These are
#: RECORDS (`numerical-bug-signatures` Sig 10): a red says *something moved*,
#: not which side is right.  Their independent corroboration is
#: ``scratch/_426_repro.md`` (a separately written instrument reproducing all
#: three k values to 9 decimals).
_K_PRE_CARVE_L2 = 1.0953221881419453
_K_PRE_CARVE_L0 = 1.1587120371368607

#: `[M]` the §0 ladder at ``scattering_order=2`` with both (n,2n) ℓ=1 and ℓ=2
#: restored (arm ``A2_lmax2``).  The plan's done-when band is ±1e-4.
_K_POST_CARVE_L2 = 1.091199657

#: Regression band for an ITERATIVE solve: SAFETY(10) × keff_tol
#: (``feedback_regression_tolerance_design.md``).  `[M]` re-running the same
#: arm at keff_tol=1e-7 moves k by 1.3e-8, so 1e-8 is the honest floor and it
#: sits 5 orders below the Δ this file measures.
_ITER_RTOL = 1e-8


def _mixture(names, densities):
    """Pure-isotope macroscopic mixture from the shipped library.

    ⚠ The σ₀ solve divides by zero on a pure isotope (``sigma_zeros.py:63``,
    ``interpolation.py:36``) and clips to σ₀ = 1 b — identically in every arm,
    so it cancels from the Δ.  The filter is NARROW and the fired warnings are
    surfaced, never blanket-suppressed (`qa` finding D4).
    """
    isos = [load_isotope(n, _TEMP_K) for n in names]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with contextlib.redirect_stdout(io.StringIO()):
            mix = compute_macro_xs(isos, np.asarray(densities), escape_xs=0.0)
    unexpected = [
        w for w in caught
        if not issubclass(w.category, RuntimeWarning)
        or ("divide by zero" not in str(w.message)
            and "invalid value" not in str(w.message))
    ]
    if unexpected:
        pytest.fail(
            "unexpected warnings from the σ₀ path — the pure-isotope clip is "
            f"the only one this fixture licenses: {[str(w.message) for w in unexpected]}"
        )
    return mix


@pytest.fixture(scope="module")
def library():
    """`[M]` 0.29 s, once for the whole module."""
    return {0: _mixture(["U_235"], [_N_U_METAL]), 1: _mixture(["BE009"], [_N_BE])}


@pytest.fixture(scope="module")
def mesh():
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(
            Region(mat_id=1, outer_thickness_cm=3.0),
            Region(mat_id=0, outer_thickness_cm=4.0),
            Region(mat_id=1, outer_thickness_cm=3.0),
        ),
        bcs=(BC.vacuum, BC.vacuum),
    )
    return Mesh1D.from_geometry(
        geom,
        region_meshes=(
            RegionMesh(n_cells=12), RegionMesh(n_cells=16), RegionMesh(n_cells=12),
        ),
    )


@pytest.fixture(scope="module")
def quadrature():
    return Quadrature.gauss_legendre(n_ordinates=8)


def _truncate_n2n(mix, lmax: int):
    r"""Zero every (n,2n) moment above ``lmax`` — **keeping the stack LENGTH**.

    ⛔ NEVER build the control by SHORTENING the list.  `[M]` the solver clamps
    ``L = min(scattering_order, min(len(m.SigS) - 1))`` — the SCATTERING stack
    alone (ruling O-1), and the (n,2n) stack is brought to that order by
    ``TransferKernel.at_order`` — so shortening would not confound today; but
    a control that varies ONE thing (the VALUES) stays honest under any future
    clamp.  Zero the values; vary one thing.
    """
    stack = mix.Sig2
    shape = stack[0].shape
    kept = [
        s if l <= lmax else csr_matrix(shape) for l, s in enumerate(stack)
    ]
    return replace(mix, Sig2=kept)


def _solve(materials, mesh, quadrature, order: int) -> float:
    with contextlib.redirect_stdout(io.StringIO()):
        sol = solve_sn(
            materials, mesh, quadrature,
            scattering_order=order,
            keff_tol=_KEFF_TOL, flux_tol=_FLUX_TOL,
            inner_tol=_INNER_TOL, max_outer=_MAX_OUTER,
        )
    history = sol.history
    if history is None or not history.fully_converged:
        pytest.fail(
            f"the arm did not fully converge (scattering_order={order}); a "
            f"starved solve degrades the RATE, not the limit, so no budget "
            f"certifies this tolerance (#340 N5)"
        )
    keff = sol.keff
    if keff is None:
        pytest.fail(f"no eigenvalue returned (scattering_order={order})")
    return float(keff)


@pytest.fixture(scope="module")
def ladder(library, mesh, quadrature):
    """Three arms at ``scattering_order=2`` — elastic P2 FIXED in all three,
    only the (n,2n) truncation varies.  `[M]` 3 × 6.84 s."""
    def arm(lmax):
        mats = {mid: _truncate_n2n(m, lmax) for mid, m in library.items()}
        return _solve(mats, mesh, quadrature, 2)

    return {"p0": arm(0), "p1": arm(1), "p2": arm(2)}


class TestTheBeReflectedAnisotropy:
    def test_the_p0_control_reproduces_the_pre_carve_eigenvalue(self, ladder):
        r"""Zeroing the (n,2n) ℓ ≥ 1 moments must return the P0 model's answer.

        This is the arm that proves the Δ below is the VALUES and not the
        machinery (`qa` control (e): the rebuild is bit-inert).  ``rtol`` is
        the iterative regression band, NOT ``array_equal``: F2/F3 route the P0
        (n,2n) emission through the shared transfer machinery, which may
        re-associate the reduction — a bit-identity claim here would be
        arithmetically impossible to honour (lessons §6, the
        which-reductions-are-reordered check).
        """
        np.testing.assert_allclose(
            ladder["p0"], _K_PRE_CARVE_L2, rtol=_ITER_RTOL,
            err_msg=(
                f"the P0-only control moved: got {ladder['p0']!r}, frozen "
                f"pre-carve {_K_PRE_CARVE_L2!r}"
            ),
        )

    def test_it_reads_the_measured_anisotropic_eigenvalue(self, ladder):
        r"""The plan's done-when 2: ``k = 1.0912 ± 1e-4`` with the SHIPPED
        library and no probe.  Under the P0 model this read 1.0953221881419453."""
        assert abs(ladder["p2"] - _K_POST_CARVE_L2) < 1e-4, (
            f"k = {ladder['p2']!r}; the §0 measurement (independently "
            f"reproduced, scratch/_426_repro.md) is {_K_POST_CARVE_L2!r} "
            f"± 1e-4.  Under the P0 model this reads {_K_PRE_CARVE_L2!r}."
        )

    def test_the_delta_is_negative_and_of_the_measured_size(self, ladder):
        r"""The PHYSICS row: forward-peaked (n,2n) emission in the reflector
        LOWERS k, by 300–500 ·1e-5 in Δk on this fixture.

        ⚠ CONVENTION (`vv` #35, and the `qa` D1 finding that produced it):
        ``pcm`` is overloaded.  This row is stated in **absolute Δk·1e5**;
        `[M]` the same pair reads −413.55 absolute / −377.56 relative /
        −346.01 in reactivity, and on the k₀ = 1.53 sibling fixture the
        absolute and reactivity readings *invert* the thin-vs-thick
        comparison.  The band below is absolute and the docstring says so.
        """
        dk = (ladder["p2"] - ladder["p0"]) * 1e5
        assert dk < 0.0, (
            f"Δk·1e5 = {dk:+.2f} — restoring a FORWARD-peaked emission "
            f"(μ̄ = +0.30 on 50 of 50 live Be-9 groups) must LOWER k; a "
            f"positive reading means the moments entered with the wrong sign"
        )
        assert 300.0 < abs(dk) < 500.0, (
            f"Δk·1e5 = {dk:+.2f}, outside the measured 300–500 band "
            f"(§0: −412.25 absolute on this fixture).  A factor-2 slip in the "
            f"multiplicity would read ≈2× (the qa linearity control measured "
            f"the response linear to 0.16 %)."
        )

    def test_the_ladder_has_converged_by_the_first_moment(self, ladder):
        r"""ℓ = 1 carries the effect; ℓ = 2 adds < 5 % of it.

        `[M]` §0: |k(ℓ≤2) − k(ℓ≤1)| = 1.30e-5 against |k(ℓ≤1) − k(P0)| =
        4.14e-3 — a ratio of 3.1e-3.  A gate at 0.05 is 16× loose and would
        red on a moments-off-by-one (which would put the ℓ=2 content into the
        ℓ=1 slot and move the increment by O(1)).
        """
        first = abs(ladder["p1"] - ladder["p0"])
        second = abs(ladder["p2"] - ladder["p1"])
        assert first > 0.0, (
            "the ℓ=1 arm did not move — the (n,2n) stack carries no ℓ=1 "
            "moment, or the binding does not read it (the P0 model)"
        )
        assert second < 0.05 * first, (
            f"the ℓ-ladder is not converging: |Δ(ℓ=2)| = {second:.3e} against "
            f"|Δ(ℓ=1)| = {first:.3e}"
        )

    def test_the_p0_solve_is_untouched_by_the_carve(
        self, library, mesh, quadrature,
    ):
        r"""``scattering_order = 0`` must read the frozen pre-carve value.

        The carve's inertness leg: at ``L = 0`` the transfer family brings
        BOTH channels to ℓ = 0, so the whole ℓ ≥ 1 machinery is out of the
        picture and the answer is the P0 model's.  A red here says the carve
        perturbed the P0 path — the one thing it must not do.

        `[M]` pre-carve 1.1587120371368607 in 4.42 s.
        """
        k = _solve(library, mesh, quadrature, 0)
        np.testing.assert_allclose(
            k, _K_PRE_CARVE_L0, rtol=_ITER_RTOL,
            err_msg=(
                f"the P0 solve moved: got {k!r}, frozen pre-carve "
                f"{_K_PRE_CARVE_L0!r}"
            ),
        )
