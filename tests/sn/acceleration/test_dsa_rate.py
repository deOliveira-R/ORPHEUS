r"""The 3c rate/stability tier — D11–D13, S2 exactness, the negative
controls, and the #215 catchers (D9, D10, M-D3-SIGMAR → ERR-070).

The correction→0 partition (verification spec §0) makes this tier
load-bearing: FP-invariance (3b's D3/D4) is structurally blind to every
accelerator-quality error, so the RATE gates here are the only catchers
for seven of the eight canonical implementation errors. Every numeric
bar below is pinned to a measured value (the 3c design scan,
2026-07-26; regenerate via ``rate_design_scan.py`` in the job tmp — the
durable copies live in ``.claude/plans/dsa_d2_characterization.md`` and
the Phase-3 roadmap).

Measured design facts (migrated to the theory page)
===================================================

The measured design facts these bars encode — the D11 one-sided
Fourier bound and the plain-SI honesty control, the reflective
Jacobi-wall lag (production rate > matrix rho, not a consistency
failure), the S2 K2=0 one-iteration exactness, the c -> 1 FP floor
and scale-free metric, the weighted-diamond partial-consistency
negative control, and the P1-DSA (d1) anisotropy-ladder payoff — are
documented as teaching tables in
``docs/theory/methods/sn/acceleration.rst`` (the
``sn-dsa-rate-and-stability`` section). This docstring records only
the test CONTRACT; the durable numeric copies live in
``.claude/plans/dsa_rate_characterization.md`` and the Phase-3
roadmap.

The #215 catchers (mutation matrix row 7)
=========================================

The σ_r-fold — realizing the within-group scattering as a diagonal
σ_r-sweep, exact only for isotropic flux — is the ONE canonical error
that changes the fixed point. Its catchers here:

* **M-D3-SIGMAR** (``TestSigmaRFoldCaught``): the folded operator's
  fixed point measurably departs the true one on an anisotropic-flux
  config (vacuum + heterogeneous — the Mode-9-honest choice; the
  isotropic reflective box is the designed-green degeneracy), far
  outside D3's 1e-7 identity band ⟹ the FP-invariance gate has the
  fold in its teeth. Files **ERR-070**.
* **D10** (``TestD10RoutingSentinel``): the structural fence — the
  foldable accessors' production consumers are exactly the definition
  site, the split layer, and the DSA low-order build; a new consumer
  (someone wiring the fold into a sweep) reds the sentinel at design
  time, before any numerics run.
"""
from __future__ import annotations

import ast
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import lu_factor
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.xs_library import get_mixture, make_mixture
from orpheus.geometry import BC
from orpheus.geometry.mesh import Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.acceleration.dsa import DSACorrection, DSALowOrderSystem
from orpheus.sn.solver import Solution, solve_sn_fixed_source

#: The A&L (3.65) continuum Fourier bound coefficient: ρ ≤ 0.2247 c.
RHO_BOUND = 0.2247

_TOL = 1e-11
#: FP-identity band (same rationale as the 3b file: both runs stop on
#: the equation residual < _TOL; ×10⁴ headroom absorbs 1/(1−ρ)).
_FP_RTOL = 1e-7


# ── Fixtures ──────────────────────────────────────────────────────────


def _mix_1g(c: float, sigma_t: float = 1.0) -> Mixture:
    return make_mixture(
        sig_t=np.array([sigma_t]),
        sig_c=np.array([sigma_t * (1.0 - c)]),
        sig_f=np.array([0.0]),
        nu=np.array([0.0]),
        chi=np.array([0.0]),
        sig_s=np.array([[sigma_t * c]]),
    )


def _uniform_solve_raw(
    c: float,
    sth: float,
    bc: tuple[str, str],
    *,
    k: int = 40,
    n_ord: int = 8,
    acceleration: str | None = None,
    inner_solver: str | None = None,
    max_inner: int = 4000,
) -> Solution:
    """1G homogeneous slab, K cells of optical thickness ``sth``."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, k * sth, k + 1),
        mat_ids=np.zeros(k, dtype=int),
        bc_left=BC(bc[0]),
        bc_right=BC(bc[1]),
    )
    return solve_sn_fixed_source(
        materials={0: _mix_1g(c)},
        mesh=mesh,
        quadrature=Quadrature.gauss_legendre(n_ordinates=n_ord),
        external_source=np.ones((n_ord, 1, k)),
        scattering_order=0,
        inner_tol=_TOL,
        max_inner=max_inner,
        inner_solver=inner_solver,
        acceleration=acceleration,
    )


@lru_cache(maxsize=None)
def _uniform_solve(
    c: float,
    sth: float,
    bc: tuple[str, str],
    *,
    k: int = 40,
    n_ord: int = 8,
    acceleration: str | None = None,
    inner_solver: str | None = None,
    max_inner: int = 4000,
) -> Solution:
    """Cached healthy-path runs (results shared across gates).

    Mutation teeth MUST call ``_uniform_solve_raw`` — a monkeypatched
    run served from this cache would be vacuously green.
    """
    return _uniform_solve_raw(
        c, sth, bc, k=k, n_ord=n_ord, acceleration=acceleration,
        inner_solver=inner_solver, max_inner=max_inner,
    )


@pytest.fixture(autouse=True)
def _no_cache_across_mutation(request):
    """The Pattern-4 cache guard: ``lru_cache`` is invisible to
    ``monkeypatch``, so a mutated run served from (or poisoning) the
    healthy cache would be a vacuous green — the one failure mode that
    silently DISARMS a gate. Any test declaring ``monkeypatch`` gets a
    cleared cache on both sides. Scope: ``_uniform_solve`` is this
    module's ONE cached wrapper (``_p1_solve`` and the teeth helpers
    are uncached BY DESIGN) — cache any new wrapper and this guard
    must clear it too, or the footgun reopens."""
    is_mutation = "monkeypatch" in request.fixturenames
    if is_mutation:
        _uniform_solve.cache_clear()
    yield
    if is_mutation:
        _uniform_solve.cache_clear()


def _residuals(solution: Solution) -> tuple[float, ...]:
    history = solution.history
    if history is None or not history.flux_residuals:
        pytest.fail("the fixed-source Solution must carry flux_residuals")
    return tuple(history.flux_residuals)


def _inners(solution: Solution) -> int:
    history = solution.history
    if history is None or history.n_inner is None:
        pytest.fail("the fixed-source Solution must carry n_inner history")
    return int(history.n_inner)


def _phi(solution: Solution) -> np.ndarray:
    v = solution.scalar_flux
    return np.asarray(getattr(v, "values", v), dtype=float)


def _rho_est(residuals: tuple[float, ...], k_tail: int = 8) -> tuple[float, float]:
    """The ρ-honest estimator: geometric mean of the tail of successive
    EQUATION-residual ratios (never the increment — Signature 9), with
    the tail's coefficient of variation as the asymptotia diagnostic."""
    r = np.asarray(residuals, dtype=float)
    r = r[r > 0.0]
    if r.size < 3:
        pytest.fail(f"too few residuals for a rate estimate ({r.size})")
    ratios = r[1:] / r[:-1]
    tail = ratios[-k_tail:] if ratios.size >= k_tail else ratios
    geo_mean = float(np.exp(np.mean(np.log(tail))))
    cv = float(np.std(tail) / np.mean(tail))
    return geo_mean, cv


# ── D11 — measured ρ vs the A&L (3.65) Fourier bound ─────────────────


class TestD11SpectralRadius:
    """D11 — the production SI+DSA spectral radius respects the
    continuum Fourier bound ρ ≤ 0.2247c (1-group is LEGITIMATE: rate is
    a flux-shape-independent claim; the Fourier analysis is itself 1G).

    Gate design (deviation from the spec's pre-measurement ±0.03 band,
    measured 2026-07-26): the discrete S8 ρ sits BELOW the continuum
    sup (0.176–0.180 vs 0.2022 at c = 0.9) — as it must; the bound is
    one-sided.  Three legs close the estimator's blind spots: the bound
    (theory), the floor (a dead/collapsed measurement cannot pass), and
    the plain-SI control (the estimator tracks the operator: ρ_est ≈ c).
    """

    pytestmark = pytest.mark.l1

    @pytest.mark.verifies("sn-dsa-consistent-fourier")
    @pytest.mark.parametrize("sth", [0.5, 1.0])
    def test_rho_below_bound_c09(self, sth):
        c = 0.9  # the scattering ratio — the bound is ρ ≤ 0.2247·c
        sol = _uniform_solve(c, sth, ("vacuum", "vacuum"),
                             acceleration="dsa")
        rho, cv = _rho_est(_residuals(sol))
        if cv >= 0.05:
            pytest.fail(f"ρ tail not asymptotic (cv = {cv:.3f})")
        if not (0.14 <= rho <= RHO_BOUND * c):
            pytest.fail(
                f"ρ_est = {rho:.4f} outside [floor 0.14, Fourier bound "
                f"0.2247·c = {RHO_BOUND * c:.4f}] (c = {c}, σ_t·h = {sth})"
            )

    def test_rho_below_bound_c05(self):
        # The c = 0.5 run converges in ~10 iterations, so the geo-mean
        # is transient-inclusive and under-reports (measured 0.074 vs a
        # plateau ~0.085) — bound + floor only, no cv guard.
        c = 0.5
        sol = _uniform_solve(c, 0.5, ("vacuum", "vacuum"),
                             acceleration="dsa")
        rho, _ = _rho_est(_residuals(sol))
        if not (0.04 <= rho <= RHO_BOUND * c):
            pytest.fail(
                f"ρ_est = {rho:.4f} outside [floor 0.04, Fourier bound "
                f"0.2247·c = {RHO_BOUND * c:.4f}] (c = {c})"
            )

    def test_plain_si_rho_matches_c(self):
        """The estimator-honesty control (doubles as M-D11-NO-ACCEL:
        without the accelerator the SAME estimator must report ρ ≈ c,
        four-and-a-half times the DSA value)."""
        sol = _uniform_solve(0.9, 0.5, ("vacuum", "vacuum"))
        rho, cv = _rho_est(_residuals(sol))
        if cv >= 0.05:
            pytest.fail(f"plain-SI ρ tail not asymptotic (cv = {cv:.3f})")
        if abs(rho - 0.9) >= 0.02:
            pytest.fail(
                f"plain-SI ρ_est = {rho:.4f} must track c = 0.9 ± 0.02 — "
                f"the estimator does not measure the iteration operator"
            )


# ── D12 — reflective stability (the historical failure mode) ─────────


class TestD12ReflectiveStability:
    """D12 — the fully-reflective box, c → 1, thick cells: exactly the
    regime where inconsistent DSA historically diverged (no leakage
    damps the transport↔diffusion mismatch). Consistent DSA must
    converge with a bounded, thickness-independent rate — AND to the
    right fixed point (a boolean convergence gate's stabiliser contains
    "converges to a wrong stable FP"; the value pairing closes it)."""

    pytestmark = pytest.mark.l1

    @pytest.mark.parametrize("c", [0.9, 0.99])
    @pytest.mark.parametrize("sth", [1.0, 5.0, 20.0])
    def test_fully_reflective_converges_bounded(self, c, sth):
        sol = _uniform_solve(c, sth, ("reflective", "reflective"),
                             acceleration="dsa")
        res = _residuals(sol)
        n = _inners(sol)
        if not (res[-1] < _TOL and n < 4000):
            pytest.fail(
                f"SI+DSA must converge on the reflective box "
                f"(c={c}, σ_t·h={sth}): n={n}, final={res[-1]:.2e}"
            )
        rho, _ = _rho_est(res)
        # Measured 0.28 (c=0.9) / 0.31 (c=0.99), flat in thickness —
        # the production reflective rate is lag-limited (the Jacobi
        # wall gain), NOT the vacuum Fourier value.
        if not (rho <= 0.35 and n <= 30):
            pytest.fail(
                f"reflective rate out of band (c={c}, σ_t·h={sth}): "
                f"ρ_est={rho:.4f} (≤0.35), n={n} (≤30)"
            )

    def test_reflective_value_pairing_c09(self):
        plain = _uniform_solve(0.9, 1.0, ("reflective", "reflective"))
        dsa = _uniform_solve(0.9, 1.0, ("reflective", "reflective"),
                             acceleration="dsa")
        np.testing.assert_allclose(
            _phi(dsa), _phi(plain), rtol=_FP_RTOL, atol=0,
            err_msg="reflective SI+DSA must share the plain-SI FP",
        )

    @pytest.mark.slow
    def test_reflective_value_pairing_c099(self):
        """The c = 0.99 pairing — the plain-SI reference alone needs
        ~3000 iterations (ρ ≈ 0.99); slow-marked, the stability rows
        above stay in the inner loop."""
        plain = _uniform_solve(0.99, 5.0, ("reflective", "reflective"),
                               max_inner=8000)
        dsa = _uniform_solve(0.99, 5.0, ("reflective", "reflective"),
                             acceleration="dsa")
        np.testing.assert_allclose(
            _phi(dsa), _phi(plain), rtol=_FP_RTOL, atol=0,
            err_msg="reflective SI+DSA must share the plain-SI FP "
                    "(c = 0.99)",
        )

    def test_lumped_removal_mass_blows_up_the_count(self, monkeypatch):
        """M-D12-INCONSISTENT-LO, in its MEASURED form: lumping the
        consistent ¼(1,2,1) removal mass onto the diagonal (the classic
        'discretize-the-reduced' finite-difference choice) degrades the
        thick-cell rate: the S4 matrix scan measured ρ 0.157 → 0.807
        (refl/vac, c = 0.99, σ_t·h = 30; at this test's S8 the
        consistent value is 0.189, the lumped ≈ 0.81 — order-
        insensitive) — a > 3× iteration-count blow-up. (Full
        DIVERGENCE, ρ up to 54.7, is the landed cell-centered operator's
        measured behaviour — the D2 report rows — and the WD pairing
        below; lumping alone is the degraded-not-divergent member.)"""
        c, sth, bc = 0.99, 30.0, ("reflective", "vacuum")
        healthy = _uniform_solve(c, sth, bc, acceleration="dsa")
        n_healthy = _inners(healthy)

        original = DSALowOrderSystem.from_sn_mesh.__func__

        def lumped_from_sn_mesh(cls, sn_mesh, *, scattering_order=0):
            system = original(cls, sn_mesh, scattering_order=scattering_order)
            # Rebuild the (1, K+1, K+1) operator with the removal mass
            # lumped onto the wall/diagonal nodes; SAME leakage, SAME G,
            # SAME discrete Marshak — the mutation isolates the
            # consistent-vs-lumped mass treatment (1G uniform config,
            # data closed over from the fixture).
            #
            # TWIN CROSS-REF (deliberate, enforcer-ruled): the non-mass
            # structure below MUST track `DSALowOrderSystem._build`
            # (orpheus/sn/acceleration/dsa.py) line for line — only the
            # mass DISTRIBUTION differs (diagonal lump vs the consistent
            # ¼(1,2,1) spread). If `_build`'s leakage, Marshak, or γ_N
            # spelling changes, port it here or this tooth compares
            # apples to oranges. Collapse trigger: if a second lumped
            # variant is ever needed, refactor to a mass-delta on
            # `system.a_low` instead of a third transcription.
            h = np.diff(np.asarray(sn_mesh.mesh.edges, dtype=float))
            k = h.shape[0]
            st = np.full(k, 1.0)
            ss = np.full(k, c)
            mu = np.asarray(sn_mesh.quad.mu_x, dtype=float)
            w = np.asarray(sn_mesh.quad.weights, dtype=float)
            omega = w / 2.0
            up = mu > 0
            gamma_n = float(omega[up] @ mu[up])
            c1 = float(omega[up] @ mu[up] ** 2) / float(omega @ mu**2)
            dh = (1.0 / (3.0 * st)) / h
            quarter = 0.25 * (st - ss) * h
            a = np.zeros((k + 1, k + 1))
            j = np.arange(1, k)
            a[j, j - 1] = -dh[j - 1]
            a[j, j] = dh[j - 1] + dh[j] + 2 * (quarter[j - 1] + quarter[j])
            a[j, j + 1] = -dh[j]
            for row, ci, side, kind in (
                (0, 0, "left", bc[0]), (k, k - 1, "right", bc[1]),
            ):
                sgn = +1.0 if side == "left" else -1.0
                f0_cols = {ci: dh[ci] + 2.0 * sgn * quarter[ci],
                           ci + 1: -dh[ci]}
                if kind == "vacuum":
                    orient = +1.0 if side == "left" else -1.0
                    a[row, row] += gamma_n
                    for col, val in f0_cols.items():
                        a[row, col] += orient * c1 * val
                else:  # reflective — the (39) row, lumped one-sided f1
                    for col, val in f0_cols.items():
                        a[row, col] += val
            return DSALowOrderSystem(
                a_low=a[None], g_map=system.g_map, _lu=(lu_factor(a),),
                _dh=system._dh, _a_coef=system._a_coef,
            )

        monkeypatch.setattr(
            DSALowOrderSystem, "from_sn_mesh",
            classmethod(lumped_from_sn_mesh),
        )
        mutated = _uniform_solve_raw(c, sth, bc, acceleration="dsa")
        n_mutated = _inners(mutated)
        if not n_mutated > 3 * n_healthy:
            pytest.fail(
                f"the lumped-removal low-order must blow up the count "
                f"(measured design ρ 0.157 → 0.807): mutated n = "
                f"{n_mutated} !> 3 × {n_healthy} — the consistent mass "
                f"is unconstrained"
            )


# ── D13 — iteration-count CI gates ───────────────────────────────────


class TestD13IterationCounts:
    """D13's HARD CI subset — the rate-regression tripwires. The full
    c × thickness × BC × solver characterization grid is the slow
    table (``TestD13CharacterizationGrid``); these four bars are the
    inner-loop floor under it."""

    pytestmark = pytest.mark.l1

    @pytest.mark.parametrize(
        "bc", [("vacuum", "vacuum"), ("reflective", "reflective")]
    )
    def test_si_dsa_speedup_floor(self, bc):
        """Measured 225→15 (vacuum) and 249→20 (reflective) at c = 0.9;
        the 0.5× floor is a loose-but-real regression tripwire."""
        plain = _uniform_solve(0.9, 0.5, bc)
        dsa = _uniform_solve(0.9, 0.5, bc, acceleration="dsa")
        if not _inners(dsa) < 0.5 * _inners(plain):
            pytest.fail(
                f"SI+DSA speedup floor lost ({bc}): "
                f"{_inners(dsa)} !< 0.5 × {_inners(plain)}"
            )

    def test_c_independence(self):
        """The whole point of DSA: the count stays ~flat as c → 1
        (measured 14 → 17; plain SI goes 225 → ~2300)."""
        n_09 = _inners(_uniform_solve(0.9, 1.0, ("vacuum", "vacuum"),
                                      acceleration="dsa"))
        n_099 = _inners(_uniform_solve(0.99, 1.0, ("vacuum", "vacuum"),
                                       acceleration="dsa"))
        if not n_099 < 1.5 * n_09:
            pytest.fail(
                f"c-independence lost: n(c=0.99) = {n_099} !< "
                f"1.5 × n(c=0.9) = 1.5 × {n_09}"
            )

    @pytest.mark.parametrize("sth", [0.1, 1.0, 10.0, 100.0])
    def test_thickness_independence(self, sth):
        """Unconditional stability, counted: the accelerated count is
        bounded across four decades of cell optical thickness (the
        consistency stress axis where partial schemes fail)."""
        n = _inners(_uniform_solve(0.9, sth, ("vacuum", "vacuum"),
                                   acceleration="dsa"))
        if not n <= 40:
            pytest.fail(
                f"thickness-independence lost: n = {n} > 40 at "
                f"σ_t·h = {sth}"
            )

    def test_krylov_preconditioner_helps(self):
        """The D4 teeth (spec: D4 proves safety, D13 proves the
        preconditioner works)."""
        plain = _uniform_solve(0.9, 0.5, ("vacuum", "vacuum"),
                               inner_solver="krylov")
        dsa = _uniform_solve(0.9, 0.5, ("vacuum", "vacuum"),
                             inner_solver="krylov", acceleration="dsa")
        if not _inners(dsa) < _inners(plain):
            pytest.fail(
                f"the DSA preconditioner must reduce GMRES iterations: "
                f"{_inners(dsa)} !< {_inners(plain)}"
            )


@pytest.mark.slow
class TestD13CharacterizationGrid:
    """The full accelerated-rate grid (characterization — its failure
    signals a rate regression, not a wrong answer; the plain-SI columns
    live in the committed characterization report, not here — a 25 000-
    iteration plain-SI column proves nothing a 0.5× floor doesn't).

    The cell metric is the SCALE-FREE decade count — iterations to a
    10-decade residual reduction — because the grid's c → 1 corner
    measured two facts an absolute-residual bar conflates (probe,
    2026-07-26):

    * **The FP floor**: at c = 0.999 the fixed-source flux scale is
      O(1/(1−c)) ≈ 10³, so the achievable ABSOLUTE residual floor
      (≈ scale × 1e-14) sits at the 1e-11 stop itself — the run
      contracts beautifully (matrix ρ ≤ 0.205 at every grid corner:
      the operator is healthy), then stalls on double precision, not
      on the rate. Callers at c → 1 must scale ``inner_tol`` with the
      flux scale; the decade count is the honest rate measurand.
    * **The double-lag reflective mode**: at σ_t·h = 100 fully
      reflective, the production splitting's rate is the LAGGED-WALL
      mode (ρ_prod ≈ 0.745, vs 0.065 for the matrix operator) —
      convergent and c-INDEPENDENT (75 iterations at c = 0.9 and
      0.99 alike), a property of the Jacobi wall lag at extreme wall
      decoupling, not of the low-order consistency (which the matrix
      ρ certifies). Characterized here; improvement filed at
      Phase-3 close.
    """

    pytestmark = pytest.mark.l1

    @pytest.mark.parametrize("bc", [("vacuum", "vacuum"),
                                    ("reflective", "reflective")])
    @pytest.mark.parametrize("c", [0.9, 0.95, 0.99, 0.999])
    @pytest.mark.parametrize("sth", [0.1, 1.0, 10.0, 100.0])
    def test_accelerated_decade_count_bounded_on_grid(self, bc, c, sth):
        sol = _uniform_solve_raw(c, sth, bc, acceleration="dsa",
                                 max_inner=200)
        res = np.asarray(_residuals(sol), dtype=float)
        below = np.nonzero(res < res[0] * 1e-10)[0]
        if below.size == 0:
            pytest.fail(
                f"grid cell (bc={bc}, c={c}, σ_t·h={sth}): no 10-decade "
                f"reduction within {len(res)} iterations (last "
                f"{res[-1]:.2e} from {res[0]:.2e}) — the reflective "
                f"column must stay finite (a stall here is a D12 "
                f"stability failure surfacing in the table)"
            )
        n_dec = int(below[0])
        # Measured extremes: ≤ 20 everywhere except the double-lag
        # reflective σ_t·h = 100 column (≈ 63–84); 100 bounds both.
        if not n_dec <= 100:
            pytest.fail(
                f"grid cell (bc={bc}, c={c}, σ_t·h={sth}): 10-decade "
                f"count {n_dec} > 100 — rate regression"
            )


# ── S2 exactness (K₂ = 0) ────────────────────────────────────────────


class TestS2Exactness:
    """The sharpest whole-chain unit test: for S2-GL the two-moment
    reduction closes the transport system EXACTLY (K₂ = 0), so one DSA
    correction annihilates the error — and each LAGGED reflective wall
    costs exactly one extra iteration before the same machine-zero
    landing (the trace arm feeds the corrected outflow; the second
    pass closes). A wrong sign anywhere in T, G, A_edge, the Marshak
    row, the trace arm, or the (28a) update breaks the exact landing.
    """

    pytestmark = pytest.mark.l1

    @pytest.mark.parametrize(
        "bc,n_expected",
        [
            (("vacuum", "vacuum"), 2),
            (("reflective", "vacuum"), 3),
            (("reflective", "reflective"), 3),
        ],
    )
    @pytest.mark.verifies("sn-dsa-s2-exactness")
    def test_one_correction_exactness(self, bc, n_expected):
        sol = _uniform_solve(0.9, 1.0, bc, n_ord=2, acceleration="dsa",
                             max_inner=50)
        res = _residuals(sol)
        n = _inners(sol)
        if not (n == n_expected and res[-1] < 1e-13):
            pytest.fail(
                f"S2 exactness broken ({bc}): n = {n} (expect exactly "
                f"{n_expected} — one correction + one lagged pass per "
                f"reflective wall), final residual {res[-1]:.2e} "
                f"(expect machine zero < 1e-13)"
            )

    def test_heterogeneous_exactness(self):
        """K₂ = 0 is a property of the ANGULAR closure, not of material
        uniformity — the exact landing must survive heterogeneity."""
        mesh = Mesh1D(
            edges=np.linspace(0.0, 20.0, 21),
            mat_ids=np.array([0] * 10 + [1] * 10),
            bc_left=BC("reflective"),
            bc_right=BC("vacuum"),
        )
        sol = solve_sn_fixed_source(
            materials={0: _mix_1g(0.9), 1: _mix_1g(0.5, sigma_t=2.0)},
            mesh=mesh,
            quadrature=Quadrature.gauss_legendre(n_ordinates=2),
            external_source=np.ones((2, 1, 20)),
            scattering_order=0,
            inner_tol=_TOL,
            max_inner=50,
            acceleration="dsa",
        )
        res = _residuals(sol)
        if not (_inners(sol) == 3 and res[-1] < 1e-13):
            pytest.fail(
                f"heterogeneous S2 exactness broken: n = {_inners(sol)}, "
                f"final = {res[-1]:.2e}"
            )


# ── The P1-DSA arm (the d₁ moment pair — R5 ruling 2026-07-26) ───────


def _mix_1g_p1(c0: float, eta: float) -> Mixture:
    """1G mixture with σ_s0 = c0·σ_t and σ_s1 = η·σ_s0, built on the
    ``Mixture`` constructor directly — ``make_mixture`` is P0-only
    (hardcodes no ℓ≥1 rows), and the ℓ=1 row IS the physics under
    test here (the lessons-L1 config trap: a make_mixture call would
    silently null it)."""
    s0 = c0
    return Mixture(
        SigC=np.array([1.0 - s0]),
        SigL=np.array([0.0]),
        SigF=np.array([0.0]),
        SigP=np.array([0.0]),
        SigT=np.array([1.0]),
        SigS=[csr_matrix(np.array([[s0]])),
              csr_matrix(np.array([[eta * s0]]))],
        Sig2=csr_matrix((1, 1)),
        chi=np.array([0.0]),
    )


def _p1_solve(
    eta: float,
    *,
    n_ord: int,
    bc: tuple[str, str] = ("vacuum", "vacuum"),
    acceleration: str | None = None,
    k: int = 40,
    max_inner: int = 3000,
) -> Solution:
    """1G homogeneous ℓ≥1 slab (c₀ = 0.9, σ_t·h = 1), sweep retaining
    ℓ = 1 — the anisotropy-ladder configuration."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, float(k), k + 1),
        mat_ids=np.zeros(k, dtype=int),
        bc_left=BC(bc[0]),
        bc_right=BC(bc[1]),
    )
    return solve_sn_fixed_source(
        materials={0: _mix_1g_p1(0.9, eta)},
        mesh=mesh,
        quadrature=Quadrature.gauss_legendre(n_ordinates=n_ord),
        external_source=np.ones((n_ord, 1, k)),
        scattering_order=1,
        inner_tol=_TOL,
        max_inner=max_inner,
        acceleration=acceleration,
    )


class TestP1DSAArm:
    """The d₁ moment-pair arm (R5 ruling, 2026-07-26: wire d₁ now).

    Consistency with the ITERATED operator applies to the restriction
    exactly as to the data row: an ℓ≥1 sweep iterates a moment PAIR,
    so the consistent correction restricts, solves, and injects the
    pair — d₁ through the (23f) g₁ columns, the (28b) cell update, and
    Larsen's (33) synthesis. Measured on the anisotropy ladder
    (σ_s1/σ_s0 ∈ {0.3, 0.6, 0.9}): the P0-only arm degrades 24/39/86
    iterations (ρ up to 0.751); the pair arm restores the FLAT
    A&L rate — 14/15/15 at ρ_est = 0.173–0.175, the same Fourier band
    as the isotropic D11 row.

    The decisive convention anchor: S2's angular space IS span{1, μ},
    so with the pair arm live one correction annihilates the FULL
    error INCLUDING the ℓ=1 gain — machine-zero landing (measured
    5.4e-15 vacuum, 4.4e-15 after the one lagged reflective pass). A
    wrong normalization anywhere in the w·μ restriction, the g₁
    columns, (28b), or the (33) synthesis breaks the exact landing.
    """

    pytestmark = pytest.mark.l1

    @pytest.mark.parametrize(
        "bc,n_expected",
        [(("vacuum", "vacuum"), 2), (("reflective", "vacuum"), 3)],
    )
    @pytest.mark.verifies("sn-dsa-synthesis")
    def test_s2_exactness_with_l1_scattering(self, bc, n_expected):
        sol = _p1_solve(0.5, n_ord=2, bc=bc, acceleration="dsa",
                        max_inner=50)
        res = _residuals(sol)
        if not (_inners(sol) == n_expected and res[-1] < 1e-13):
            pytest.fail(
                f"S2 ℓ=1 exactness broken ({bc}): n = {_inners(sol)} "
                f"(expect exactly {n_expected}), final {res[-1]:.2e} "
                f"(expect < 1e-13) — the moment-pair chain (w·μ "
                f"restriction → g₁ → (28b) → (33)) is off somewhere"
            )

    def test_p0_armed_corrector_cannot_close_s2(self, monkeypatch):
        """The arm's tooth: force the corrector back to the P0 arm on
        the SAME ℓ=1 sweep — the ℓ=1 error survives the correction and
        the exact landing breaks (this is the pre-R5 state, in
        miniature: the ladder's degradation mechanism isolated on the
        S2 exactness anchor)."""
        original = DSACorrection.from_sn_mesh.__func__

        def p0_forced(cls, sn_mesh, *, scattering_order=0):
            del scattering_order  # the mutation: the P0 arm regardless
            return original(cls, sn_mesh, scattering_order=0)

        monkeypatch.setattr(
            DSACorrection, "from_sn_mesh", classmethod(p0_forced)
        )
        sol = _p1_solve(0.5, n_ord=2, acceleration="dsa", max_inner=50)
        if not _inners(sol) > 2:
            pytest.fail(
                "a P0-armed corrector must NOT close the ℓ=1 S2 system "
                f"in one correction (n = {_inners(sol)}) — the d₁ arm "
                "is unconstrained"
            )

    def test_ladder_flattens(self):
        """The worst pre-arm rung (η = 0.9: 86 iterations, ρ 0.751)
        restored to the flat A&L band — the same bound + floor as the
        isotropic D11 row."""
        c0 = 0.9
        sol = _p1_solve(0.9, n_ord=8, acceleration="dsa")
        rho, cv = _rho_est(_residuals(sol))
        n = _inners(sol)
        if cv >= 0.05:
            pytest.fail(f"ladder ρ tail not asymptotic (cv = {cv:.3f})")
        if not (n <= 20 and 0.14 <= rho <= RHO_BOUND * c0):
            pytest.fail(
                f"the d₁ arm must flatten the anisotropy ladder: "
                f"n = {n} (≤ 20; pre-arm 86), ρ_est = {rho:.4f} "
                f"(band [0.14, {RHO_BOUND * c0:.4f}]; pre-arm 0.751)"
            )

    def test_moment_pair_injection_object_pins(self):
        """The object-level pins (D8's ℓ=1 siblings): the frame's
        slab ℓ=1 row IS w·μ (the SH table's slab component is μ
        bit-exactly), and the (33) synthesis is the exact right
        inverse of the pair restriction (R∘P = I on (φ₀, φ₁) — the
        μ-arm is moment-0-free and the iso arm moment-1-free by
        quadrature symmetry and W₂ exactness)."""
        quad = Quadrature.gauss_legendre(n_ordinates=8)
        mu = np.asarray(quad.mu_x, dtype=float)
        w = np.asarray(quad.weights, dtype=float)
        table = np.asarray(quad.angular_frame(1).table)
        np.testing.assert_array_equal(table[:, 0, 0], np.ones(mu.size))
        np.testing.assert_array_equal(table[:, 1, 1], mu)

        rng = np.random.default_rng(31)
        phi0 = rng.normal(size=(2, 5))
        phi1 = rng.normal(size=(2, 5))
        sum_w = float(w.sum())
        inv_w2 = float(1.0 / ((w / 2.0) @ mu**2))
        synth = (
            phi0[None] + inv_w2 * mu[:, None, None] * phi1[None]
        ) / sum_w
        m0 = np.einsum("n,ngk->gk", w, synth)
        m1 = np.einsum("n,ngk->gk", w * mu, synth)
        np.testing.assert_allclose(m0, phi0, rtol=1e-13, atol=1e-13)
        np.testing.assert_allclose(m1, phi1, rtol=1e-13, atol=1e-13)


# ── The partial-consistency negative control (matrix tier) ────────────


def _wd_sweep_matrix(h, st, ss, mu, w, a):
    """The weighted-diamond SI sweep matrix (ψ̄ = a·ψ_out + (1−a)·ψ_in;
    a = ½ is diamond), VACUUM walls only — the negative control's
    entire config surface (a reflective arm would be unexercised,
    unverified instrument code). Anchored below by the a = ½ composite
    S2 machine-zero, so it cannot drift from the production convention
    unnoticed."""
    k = h.shape[0]

    def sweep(phi):
        q = 0.5 * ss * phi
        phi_new = np.zeros(k)
        for m in np.argsort(mu):
            if mu[m] >= 0:
                continue
            psi_in = 0.0
            for i in range(k - 1, -1, -1):
                g1 = abs(mu[m]) / h[i]
                psi_out = (q[i] + psi_in * (g1 - (1 - a) * st[i])) / (
                    g1 + a * st[i]
                )
                phi_new[i] += w[m] * (a * psi_out + (1 - a) * psi_in)
                psi_in = psi_out
        for m in np.argsort(mu):
            if mu[m] <= 0:
                continue
            psi_in = 0.0
            for i in range(k):
                g1 = abs(mu[m]) / h[i]
                psi_out = (q[i] + psi_in * (g1 - (1 - a) * st[i])) / (
                    g1 + a * st[i]
                )
                phi_new[i] += w[m] * (a * psi_out + (1 - a) * psi_in)
                psi_in = psi_out
        return phi_new

    return np.stack([sweep(np.eye(k)[:, j]) for j in range(k)], axis=1)


def _t_dsa(a_edge, g_map, t_si):
    """T + C2E · A⁻¹ · G · [d0; 0] — the corrected error-iteration
    operator for an injected low-order (A, G) pair."""
    k = t_si.shape[0]
    c2e = np.zeros((k, k + 1))
    c2e[np.arange(k), np.arange(k)] = 0.5
    c2e[np.arange(k), np.arange(k) + 1] = 0.5
    d_stack = np.vstack([t_si - np.eye(k), np.zeros((k, k))])
    return t_si + c2e @ np.linalg.solve(a_edge, g_map @ d_stack)


def _consistent_low_order(h, st, ss, mu, w, bc):
    """The PRODUCTION low-order build (1-group arrays), so the control
    exercises the shipped operator, not a re-transcription."""
    system = DSALowOrderSystem._build(
        h, st[None], ss[None], np.zeros_like(ss)[None], mu, w, bc,
    )
    return system.a_low[0], system.g_map[0]


def _rho(m: np.ndarray) -> float:
    return float(np.abs(np.linalg.eigvals(m)).max())


class TestNegativeControlPartialConsistency:
    """The battery must detect INCONSISTENCY, not just breakage: a
    weighted-diamond sweep paired with the DD-consistent low-order is
    the Adams-Larsen partially-consistent class, and its measured ρ
    must reproduce the Table-II shape — climbing with cell thickness
    into divergence as c → 1 — while the consistent (a = ½) pairing
    stays flat on the SAME problems."""

    pytestmark = pytest.mark.foundation

    K = 40

    def _grid(self, sth, c):
        h = np.full(self.K, sth)
        st = np.full(self.K, 1.0)
        ss = np.full(self.K, c)
        return h, st, ss

    def test_dd_member_anchor(self):
        """a = ½ IS diamond: the WD instrument + production low-order
        must close the S2 system at machine zero (the K₂ = 0 anchor —
        one number ties the test-local sweep, the production build,
        and the update composition; a convention drift in ANY of them
        breaks the exact landing)."""
        q2 = Quadrature.gauss_legendre(n_ordinates=2)
        mu, w = np.asarray(q2.mu_x, float), np.asarray(q2.weights, float)
        h, st, ss = self._grid(1.0, 0.9)
        bc = ("vacuum", "vacuum")
        t_si = _wd_sweep_matrix(h, st, ss, mu, w, a=0.5)
        a_e, g_m = _consistent_low_order(h, st, ss, mu, w, bc)
        rho_s2 = _rho(_t_dsa(a_e, g_m, t_si))
        if not rho_s2 < 1e-12:
            pytest.fail(
                f"the S2 anchor must close at machine zero for the "
                f"a = ½ (diamond) member: ρ = {rho_s2:.3e}"
            )

    @pytest.mark.parametrize("sth,rho_measured", [(10.0, 1.44), (30.0, 1.78)])
    def test_wd_pairing_diverges_thick(self, sth, rho_measured):
        """a = 0.75, c = 0.99: the partial pairing DIVERGES at thick
        cells (measured ρ = 1.44 / 1.78) while the consistent pairing
        on the SAME problem stays inside the Fourier band."""
        q4 = Quadrature.gauss_legendre(n_ordinates=4)
        mu, w = np.asarray(q4.mu_x, float), np.asarray(q4.weights, float)
        h, st, ss = self._grid(sth, 0.99)
        bc = ("vacuum", "vacuum")
        a_e, g_m = _consistent_low_order(h, st, ss, mu, w, bc)

        t_wd = _wd_sweep_matrix(h, st, ss, mu, w, a=0.75)
        rho_partial = _rho(_t_dsa(a_e, g_m, t_wd))
        if not rho_partial > 1.0:
            pytest.fail(
                f"the WD(a=0.75) + DD-low-order pairing must diverge at "
                f"σ_t·h = {sth}, c = 0.99 (measured {rho_measured}); "
                f"got ρ = {rho_partial:.4f} — the battery cannot see "
                f"partial consistency"
            )

        t_dd = _wd_sweep_matrix(h, st, ss, mu, w, a=0.5)
        rho_consistent = _rho(_t_dsa(a_e, g_m, t_dd))
        if not rho_consistent < RHO_BOUND * 0.99 + 0.01:
            pytest.fail(
                f"the consistent pairing must stay in band on the same "
                f"problem: ρ = {rho_consistent:.4f}"
            )

    def test_partial_consistency_climbs_with_thickness(self):
        """The Table-II SHAPE: the partial pairing's ρ climbs an order
        of magnitude across σ_t·h ∈ [1, 30] (measured 0.20 → 1.54 at
        a = 0.6, c = 0.99); the consistent pairing's does not climb."""
        q4 = Quadrature.gauss_legendre(n_ordinates=4)
        mu, w = np.asarray(q4.mu_x, float), np.asarray(q4.weights, float)
        bc = ("vacuum", "vacuum")

        def rho_at(sth, a):
            h, st, ss = self._grid(sth, 0.99)
            a_e, g_m = _consistent_low_order(h, st, ss, mu, w, bc)
            t = _wd_sweep_matrix(h, st, ss, mu, w, a=a)
            return _rho(_t_dsa(a_e, g_m, t))

        climb_partial = rho_at(30.0, 0.6) / rho_at(1.0, 0.6)
        climb_consistent = rho_at(30.0, 0.5) / rho_at(1.0, 0.5)
        if not (climb_partial > 5.0 and climb_consistent < 1.2):
            pytest.fail(
                f"Table-II shape lost: partial climb = "
                f"{climb_partial:.2f} (must be > 5), consistent climb = "
                f"{climb_consistent:.2f} (must be < 1.2)"
            )


# ── D9 — the no-masking control ──────────────────────────────────────


class TestD9NoMasking:
    """DSA changes the RATE, never the VALUE — in BOTH directions.
    3b's D3 proved it does not move the right fixed point; D9 proves it
    faithfully reproduces a WRONG one: on a deliberately-mutated
    transport operator, the accelerated and plain solvers must agree
    (the accelerator must not launder a solver bug into a
    plausible-but-different answer)."""

    pytestmark = pytest.mark.l2

    def _solve_two_zone(self, mats, acceleration=None):
        mesh = Mesh1D(
            edges=np.linspace(0.0, 16.0, 33),
            mat_ids=np.array([0] * 16 + [1] * 16),
            bc_left=BC("vacuum"),
            bc_right=BC("vacuum"),
        )
        return solve_sn_fixed_source(
            materials=mats,
            mesh=mesh,
            quadrature=Quadrature.gauss_legendre(n_ordinates=8),
            external_source=np.ones((8, 1, 32)),
            scattering_order=0,
            inner_tol=_TOL,
            max_inner=3000,
            acceleration=acceleration,
        )

    def test_seeded_bug_reproduced_not_laundered(self):
        true_mats = {0: _mix_1g(0.9), 1: _mix_1g(0.5, sigma_t=2.0)}
        # The seeded transport bug: zone-B scattering scaled ×1.05
        # (constructed balanced, so it is a legal-but-different
        # operator — a bug in A, not in the DSA machinery).
        buggy_mats = {
            0: _mix_1g(0.9),
            1: make_mixture(
                sig_t=np.array([2.0]),
                sig_c=np.array([0.95]),
                sig_f=np.array([0.0]),
                nu=np.array([0.0]),
                chi=np.array([0.0]),
                sig_s=np.array([[1.05]]),
            ),
        }
        plain_buggy = self._solve_two_zone(buggy_mats)
        dsa_buggy = self._solve_two_zone(buggy_mats, acceleration="dsa")
        np.testing.assert_allclose(
            _phi(dsa_buggy), _phi(plain_buggy), rtol=_FP_RTOL, atol=0,
            err_msg="DSA must converge to the SAME wrong answer as "
                    "plain SI on the mutated operator",
        )
        # Activation guard (anti-Mode-10): the seeded bug must actually
        # move the FP — without this leg a null mutation passes vacuously.
        plain_true = self._solve_two_zone(true_mats)
        shift = np.max(
            np.abs(_phi(plain_buggy) - _phi(plain_true))
            / np.max(np.abs(_phi(plain_true)))
        )
        if not shift > 1e-3:
            pytest.fail(
                f"the seeded mutation is inactive (FP shift {shift:.2e}) "
                f"— the no-masking claim was not exercised"
            )


# ── M-D3-SIGMAR — the σ_r-fold in D3's teeth (ERR-070) ───────────────


def _folded(mix: Mixture) -> Mixture:
    """The #215 σ_r-fold realized as material data: the within-group
    isotropic self-scatter moved into the removal (σ_t' = σ_t −
    σ_s0^{g→g}, diagonal of Σ_s0 zeroed; cross-group transfer kept).
    Balance is preserved — the folded operator is a legal Mixture that
    is EXACT for isotropic flux and wrong for anisotropic flux."""
    s0 = np.asarray(mix.SigS[0].toarray(), dtype=float)
    diag = np.diag(s0).copy()
    return Mixture(
        SigC=np.asarray(mix.SigC, dtype=float),
        SigL=np.asarray(mix.SigL, dtype=float),
        SigF=np.asarray(mix.SigF, dtype=float),
        SigP=np.asarray(mix.SigP, dtype=float),
        SigT=np.asarray(mix.SigT, dtype=float) - diag,
        SigS=[csr_matrix(s0 - np.diag(diag))],
        Sig2=mix.Sig2,
        chi=np.asarray(mix.chi, dtype=float),
        eg=mix.eg,
    )


class TestSigmaRFoldCaught:
    """M-D3-SIGMAR — the ONE mutation-matrix row that changes the fixed
    point (spec §4 row 7): wiring the within-group solve as the
    σ_r-fold (`A_wg.solve` with the diagonal σ_r-sweep realizing
    Σ_s0·P_iso as Σ_s0·I) is exact ONLY for isotropic flux. On the
    Mode-9-honest config (vacuum + heterogeneous ⟹ anisotropic flux,
    2G) the folded operator's fixed point measurably departs the true
    one — far outside D3's 1e-7 identity band, so the FP-invariance
    gate has this class in its teeth. ERR-070 records the measured
    magnitude (the catalog's first entry for the #215 class)."""

    pytestmark = [pytest.mark.l2, pytest.mark.catches("ERR-070")]

    def _solve(self, mats, acceleration=None):
        mesh = Mesh1D(
            edges=np.linspace(0.0, 20.0, 41),
            mat_ids=np.array([0] * 20 + [1] * 20),
            bc_left=BC("vacuum"),
            bc_right=BC("vacuum"),
        )
        return solve_sn_fixed_source(
            materials=mats,
            mesh=mesh,
            quadrature=Quadrature.gauss_legendre(n_ordinates=4),
            external_source=np.ones((4, 2, 40)),
            scattering_order=0,
            inner_tol=_TOL,
            max_inner=3000,
            acceleration=acceleration,
        )

    def test_fold_moves_the_fixed_point_outside_d3s_band(self):
        true_mats = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
        fold_mats = {k: _folded(v) for k, v in true_mats.items()}

        phi_true = _phi(self._solve(true_mats))
        phi_fold = _phi(self._solve(fold_mats))
        shift = float(
            np.max(np.abs(phi_fold - phi_true) / np.max(np.abs(phi_true)))
        )
        # ERR-070: measured 43.2% peak shift on this config (11.4% /
        # 43.2% per group — #215's own configs measured 46–56%),
        # 4.3e+06× above D3's 1e-7 band.
        if not shift > 1e-2:
            pytest.fail(
                f"the σ_r-fold must move the FP on anisotropic flux "
                f"(measured shift {shift:.2e} ≤ 1e-2) — D3 would be "
                f"blind to the #215 class"
            )

        # The fold is self-consistently ACCELERATED (its own low-order
        # sees σ̂_S = 0): DSA on the folded operator reproduces the
        # folded FP — the wrongness comes from the fold, is faithfully
        # preserved by the accelerator, and only the comparison against
        # the TRUE plain FP (D3's actual assertion) exposes it.
        phi_fold_dsa = _phi(self._solve(fold_mats, acceleration="dsa"))
        np.testing.assert_allclose(
            phi_fold_dsa, phi_fold, rtol=_FP_RTOL, atol=0,
            err_msg="DSA must faithfully accelerate even the folded "
                    "operator (no laundering)",
        )


# ── D10 — the #215 routing sentinel ──────────────────────────────────


_FOLD_ACCESSORS = frozenset({"foldable_sigma", "residual_sig_s"})
#: The complete legitimate production surface for the σ_r-fold data:
#: the definition site, the split layer (whose residual_part IS the
#: safe spelling), and the DSA low-order build (the fenced consumer —
#: correction→0 makes the fold legitimate THERE and only there).
_FOLD_ALLOWLIST = frozenset({
    "orpheus/transport/mesh/material_xs_field.py",
    "orpheus/transport/operators/scattering.py",
    "orpheus/sn/acceleration/dsa.py",
})


def _fold_accessor_callers(root: Path) -> set[str]:
    """Every module under ``root`` whose AST touches a foldable
    accessor — attribute reads (``x.foldable_sigma``) OR method
    definitions (``def foldable_sigma``). Docstring/comment mentions
    are string constants and do not register (no false positives)."""
    callers: set[str] = set()
    for path in sorted(root.rglob("*.py")):
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            touches = (
                isinstance(node, ast.Attribute)
                and node.attr in _FOLD_ACCESSORS
            ) or (
                isinstance(node, ast.FunctionDef)
                and node.name in _FOLD_ACCESSORS
            )
            if touches:
                callers.add(str(path.relative_to(root.parent)))
                break
    return callers


class TestD10RoutingSentinel:
    """D10 — the structural #215 fence: the σ_r-fold data (the foldable
    accessors) must never grow a production consumer outside the three
    audited files. A sweep builder reaching for ``foldable_sigma`` reds
    this sentinel at design time — before any numerics can be wrong."""

    pytestmark = [pytest.mark.foundation, pytest.mark.catches("ERR-070")]

    def test_fold_accessors_have_no_new_production_consumer(self):
        root = Path(__file__).resolve().parents[3] / "orpheus"
        callers = _fold_accessor_callers(root)
        unexpected = callers - _FOLD_ALLOWLIST
        if unexpected:
            pytest.fail(
                f"NEW production consumer(s) of the foldable accessors: "
                f"{sorted(unexpected)}. The σ_r-fold is exact only for "
                f"isotropic flux (#215 / ERR-070) — a consumer outside "
                f"the audited set {sorted(_FOLD_ALLOWLIST)} must be "
                f"reviewed against the Mode-9 trap before extending the "
                f"allowlist."
            )
        # The wrap-fires leg (Mode 11): the fenced consumer itself IS
        # on the scanned surface — an empty caller set would mean the
        # scanner is broken, not that the tree is clean.
        if "orpheus/sn/acceleration/dsa.py" not in callers:
            pytest.fail(
                "the sentinel's scanner no longer sees the audited DSA "
                "consumer — the tripwire is scanning the wrong surface"
            )

    def test_sentinel_flags_a_planted_consumer(self, tmp_path):
        """M-D10-WIRE-FOLD: a planted module that reaches for the fold
        accessor must be flagged — the tripwire has teeth."""
        pkg = tmp_path / "orpheus"
        pkg.mkdir()
        (pkg / "planted.py").write_text(
            "def sweep(mat_xs):\n"
            "    return mat_xs.foldable_sigma()\n"
        )
        callers = _fold_accessor_callers(pkg)
        if "orpheus/planted.py" not in callers:
            pytest.fail(
                "the routing sentinel failed to flag a planted "
                "σ_r-fold consumer — the tripwire is decorative"
            )
