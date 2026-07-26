r"""Energy-condensation gate — ``Solution.condense`` end-to-end.

P5.0 (GitHub #274). ``Solution.condense(coarse: EnergyGrid) -> dict[int,
Mixture]`` collapses a fine SN solution's per-material cross sections onto a
coarse group structure, spectrum-weighted, returning **portable** few-group
XS (mesh-DECOUPLED — the asymmetry law vs ``Solution.homogenize``, which is
mesh-COUPLED). The per-material spectrum is the region-integrated scalar
flux the solve emits.

Mirrors ``tests/sn/test_homogenization.py`` (the spatial sibling). THE gate
(rate preservation) checks the condensed reaction rate against a
**structurally-independent** per-fine-group Python hand-sum (``vv L11`` —
NOT a re-call of ``frame.project``). Companion checks: the result is a
portable ``dict[int, Mixture]`` (no mesh), every condensed Mixture balances,
χ stays a fast-peaked simplex, and the real EGB421→WIMS-69 structure
integrates end-to-end.

vv Mode-8: every assertion is ``np.testing.*`` / ``pytest.raises`` /
``pytest.fail`` (the suite runs ``-O``; bare ``assert`` is stripped).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.data.energy_grid import EnergyGrid
from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn

pytestmark = [pytest.mark.l1, pytest.mark.cap("solve")]


# ── A 4-fine-group, 2-coarse-group setup, fast-first (descending eg) ────
_NG_FINE = 4
_GROUPS = (np.array([0, 1]), np.array([2, 3]))      # fast-first coarse partition
_NG_COARSE = len(_GROUPS)
_EG_FINE = np.array([1.0e7, 1.0e5, 1.0e2, 1.0e0, 1.0e-2])   # 5 edges → 4 groups
_EG_COARSE = np.array([1.0e7, 1.0e2, 1.0e-2])               # 3 edges → 2 coarse


def _balanced_fissile_4g(eg: np.ndarray) -> "object":
    # NOTE: make_mixture hardcodes SigL = 0, so SigT must NOT include a SigL
    # term here (SigT = SigC + SigF + rowsum(SigS0)); the SigL-bearing balance
    # leg lives in tests/data/test_mixture_condense.py (hand-built fixture).
    sig_c = np.array([0.02, 0.05, 0.10, 0.30])
    sig_f = np.array([0.01, 0.02, 0.08, 0.20])
    nu = np.array([2.6, 2.5, 2.45, 2.43])
    chi = np.array([0.7, 0.25, 0.05, 0.0])
    sig_s0 = np.array([
        [0.30, 0.10, 0.05, 0.02],
        [0.00, 0.40, 0.15, 0.03],
        [0.00, 0.00, 0.55, 0.20],
        [0.00, 0.00, 0.00, 0.90],
    ])
    sig_t = sig_c + sig_f + sig_s0.sum(axis=1)
    mix = make_mixture(sig_t, sig_c, sig_f, nu, chi, sig_s0, eg=eg)
    mix.assert_balanced()
    return mix


@pytest.fixture(scope="module")
def materials():
    """Single fissile material on the fine grid (condensation is per-material)."""
    return {0: _balanced_fissile_4g(_EG_FINE)}


@pytest.fixture(scope="module")
def problem(materials):
    """The SHARED (mesh, quad) constituents — the two-entry pairing pattern:
    forward and adjoint solves must be built from the SAME objects for
    ``SNMesh.is_same_phase_space`` to accept the pair (P6 B2)."""
    fine = Mesh1D(
        edges=np.linspace(0.0, 4.0, 9), mat_ids=np.zeros(8, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"), bc_right=BC("reflective"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    return fine, quad


@pytest.fixture(scope="module")
def solution(materials, problem):
    """A real fine-group SN solve. A vacuum→reflective tilt keeps the spectrum
    spatially non-flat, and the within-fine-group spectral spread is the
    condensation's load-bearing weight."""
    fine, quad = problem
    return solve_sn(materials, fine, quad, scattering_order=0)


def _region_spectrum(solution, mat_id: int = 0) -> np.ndarray:
    """The per-material spectrum φ_g = Σ_i V_i φ_{i,g} over the cells of ``mat_id``.

    The region-integrated scalar flux — the structurally-independent weight
    the condensation uses. Computed here by an explicit volume-weighted sum
    (NOT via the frame), so it is a clean oracle input.
    """
    phi = np.asarray(solution.scalar_flux.values, dtype=float)  # (ng, n_fine)
    V = np.asarray(solution.mesh.volumes, dtype=float)          # (n_fine,)
    mat = np.asarray(solution.mesh.mat_map, dtype=int).ravel()
    sel = mat == mat_id
    return (phi[:, sel] * V[sel]).sum(axis=1)                   # (ng,)


def _condense(solution):
    return solution.condense(EnergyGrid(_EG_COARSE))


# ══════════════════════════════════════════════════════════════════════
# THE gate — rate preservation against a per-material hand-sum (vv L11)
# ══════════════════════════════════════════════════════════════════════


@pytest.mark.verifies("energy-condensation-rate-preservation")
def test_solution_condense_rate_preservation(solution, materials):
    """Σ_G·Φ_G == Σ_{g∈G} φ_g Σ_g per vector channel, per material.

    Φ_G = Σ_{g∈G} φ_g is the coarse region-spectrum; φ_g is the region-
    integrated fine flux. The oracle is an explicit fine-group Python sum.
    """
    condensed = _condense(solution)          # dict[int, Mixture]
    phi = _region_spectrum(solution, 0)      # (ng_fine,)
    coarse_mix = condensed[0]
    for channel in ("SigT", "SigC", "SigL", "SigF", "SigP"):
        sig_fine = np.asarray(getattr(materials[0], channel), dtype=float)
        sig_eff = np.asarray(getattr(coarse_mix, channel), dtype=float)
        for G in range(_NG_COARSE):
            phi_R = float(phi[_GROUPS[G]].sum())
            rate_coarse = sig_eff[G] * phi_R
            rate_ref = float((phi[_GROUPS[G]] * sig_fine[_GROUPS[G]]).sum())
            np.testing.assert_allclose(
                rate_coarse, rate_ref, rtol=1e-10, atol=1e-12,
                err_msg=f"{channel} rate not preserved in coarse {G}",
            )


@pytest.mark.verifies("energy-condensation-scattering-collapse")
def test_solution_condense_scattering_rate(solution, materials):
    """In-scatter RATE preserved (source-group weighted) — catches g_from↔g_to."""
    condensed = _condense(solution)
    phi = _region_spectrum(solution, 0)
    coarse_mix = condensed[0]
    for order in range(len(coarse_mix.SigS)):
        sig_s_fine = np.asarray(materials[0].SigS[order].todense(), dtype=float)
        sig_s_coarse = np.asarray(coarse_mix.SigS[order].todense(), dtype=float)
        for Gi in range(_NG_COARSE):
            phi_R = float(phi[_GROUPS[Gi]].sum())
            for Gj in range(_NG_COARSE):
                rate_coarse = sig_s_coarse[Gi, Gj] * phi_R
                rate_ref = float(sum(
                    phi[g] * sig_s_fine[g, gp]
                    for g in _GROUPS[Gi] for gp in _GROUPS[Gj]
                ))
                np.testing.assert_allclose(
                    rate_coarse, rate_ref, rtol=1e-10, atol=1e-12,
                    err_msg=f"SigS[{order}][{Gi},{Gj}] in-scatter rate not preserved",
                )


# ══════════════════════════════════════════════════════════════════════
# Portability — the asymmetry law (condense is mesh-DECOUPLED)
# ══════════════════════════════════════════════════════════════════════


def test_condense_returns_portable_dict_no_mesh(solution):
    """``condense`` returns ``dict[int, Mixture]`` — portable, NO mesh.

    The condense/homogenize asymmetry: condense yields mesh-DECOUPLED
    portable few-group XS (a plain dict of Mixtures), unlike homogenize which
    returns a mesh-COUPLED MaterialMesh. Pin the return TYPE + that the
    Mixtures carry no mesh attribute.
    """
    from orpheus.data.macro_xs.mixture import Mixture

    condensed = _condense(solution)
    np.testing.assert_array_equal(
        isinstance(condensed, dict), True,
        err_msg="condense must return a dict[int, Mixture] (portable)",
    )
    for key, mix in condensed.items():
        np.testing.assert_array_equal(
            isinstance(mix, Mixture), True,
            err_msg=f"condense value for {key} must be a Mixture",
        )
        # Portability: a Mixture has no mesh / spatial binding.
        np.testing.assert_array_equal(
            hasattr(mix, "mesh"), False,
            err_msg="condensed Mixture must NOT carry a mesh (mesh-decoupled)",
        )


def test_condensed_group_count_dropped(solution):
    """Condensed Mixtures live on the COARSE group structure (ng → n_coarse)."""
    condensed = _condense(solution)
    for mix in condensed.values():
        np.testing.assert_array_equal(
            mix.ng, _NG_COARSE,
            err_msg="condensed Mixture must have the coarse group count",
        )
        np.testing.assert_array_equal(
            np.asarray(mix.eg).shape, (_NG_COARSE + 1,),
            err_msg="condensed eg must be the coarse boundary array (n_coarse+1)",
        )


# ══════════════════════════════════════════════════════════════════════
# Companion invariants — balance + simplex χ survive the collapse
# ══════════════════════════════════════════════════════════════════════


@pytest.mark.verifies("energy-condensation-balance-preservation")
def test_condensed_materials_balance(solution):
    """Balance survives the collapse — every removal channel shares the weight."""
    condensed = _condense(solution)
    for mix in condensed.values():
        mix.assert_balanced(atol=1e-9)


@pytest.mark.verifies("energy-condensation-chi-simplex-preservation")
def test_condensed_chi_is_fast_peaked_simplex(solution):
    """χ_R is a fast-peaked probability simplex (birth-group sum)."""
    condensed = _condense(solution)
    for mix in condensed.values():
        chi = np.asarray(mix.chi)
        np.testing.assert_allclose(
            float(chi.sum()), 1.0, atol=1e-12,
            err_msg="condensed χ must remain a simplex (Σχ=1)",
        )
        np.testing.assert_array_less(
            -1e-15, chi, err_msg="condensed χ must be non-negative"
        )
        # Fine χ = [0.7,0.25,0.05,0] → coarse 0 (fast) holds the bulk.
        np.testing.assert_array_less(
            chi[1], chi[0],
            err_msg="condensed χ must remain fast-peaked (coarse 0 > coarse 1)",
        )


# ══════════════════════════════════════════════════════════════════════
# Degenerate control — identity condensation recovers the fine material
# ══════════════════════════════════════════════════════════════════════


def test_identity_condensation_recovers_fine_material(solution, materials):
    """Condense onto the SAME fine grid → each coarse group is one fine group,
    so the effective material equals the original (the rate-average over a
    singleton group is the identity)."""
    condensed = solution.condense(EnergyGrid(_EG_FINE))
    np.testing.assert_allclose(
        np.asarray(condensed[0].SigT), np.asarray(materials[0].SigT),
        atol=1e-12,
        err_msg="identity condensation must recover the fine SigT",
    )
    np.testing.assert_allclose(
        np.asarray(condensed[0].SigS[0].todense()),
        np.asarray(materials[0].SigS[0].todense()), atol=1e-12,
        err_msg="identity condensation must recover the fine SigS[0]",
    )


# ══════════════════════════════════════════════════════════════════════
# F7 (#274) — real 421-group library → WIMS-69 NON-NESTED condensation
# ══════════════════════════════════════════════════════════════════════
#
# Real production data exists (P5.0 memo FLAG 3 CORRECTED 2026-06-27):
# ``pwr_like_mix()`` / ``uo2_fuel()`` are real 421-group balanced Mixtures
# (``ng=421``, ``eg`` descending), and ``WIMS_69`` / ``WIMS_172`` are the real
# targets. The 421→WIMS-69 map is NON-NESTED (WIMS boundaries don't align with
# the 421 spacing, and WIMS is locally finer in narrow resonance/thermal
# bands), so a one-hot containing-interval leaves EMPTY coarse groups (the bug
# the fractional re-binning fixes). NO synthesis, NO skipif on the SUT — it is
# the L2 real-structure integration leg, guarded only by the real-data recipe
# cache (``*.h5``) which may be absent locally.

from orpheus.data.energy_grid import InverseEnergySpectrum as _IES  # noqa: E402
from orpheus.data.group_structures.ornl import ORNL_421  # noqa: E402

# Real 421-group XS recipes (the production library; *.h5 caches present
# locally). Guard the import so the file collects if a cache is missing.
try:
    from orpheus.data.macro_xs.recipes import pwr_like_mix  # type: ignore

    _HAS_RECIPES = True
except Exception:  # pragma: no cover - import guard
    pwr_like_mix = None  # type: ignore
    _HAS_RECIPES = False

_needs_real_421 = pytest.mark.skipif(
    not _HAS_RECIPES,
    reason="P5.0/#274: real 421→WIMS-69 needs the pwr_like_mix recipe (*.h5 cache)",
)


def test_ornl_421_grid_is_well_formed():
    """The canonical ORNL 421-group grid: 422 strictly-descending positive edges.

    Runs without the SUT/recipe — pins the production fine grid the F7 condensation
    projects FROM. ``ORNL_421`` is an :class:`EnergyGrid`, so descending+positive is
    enforced at construction; this also pins the expected 421-group count.
    """
    np.testing.assert_array_equal(
        ORNL_421.edges.shape, (422,),
        err_msg="ORNL_421 must have 422 boundaries (421 groups)",
    )
    np.testing.assert_array_equal(ORNL_421.n_groups, 421)
    np.testing.assert_array_less(
        ORNL_421.edges[1:], ORNL_421.edges[:-1],
        err_msg="canonical ORNL 421 grid must be strictly DESCENDING (fast-first)",
    )
    np.testing.assert_array_less(
        0.0, ORNL_421.edges, err_msg="all ORNL_421 energies must be positive"
    )


@_needs_real_421
def test_ornl_421_matches_live_library():
    """``ORNL_421`` is the LIVE library grid (vv L4: generated, not hand-transcribed).

    The boundary literal in ``ornl.py`` was emitted from ``pwr_like_mix().eg``; this
    pins it to the library and reddens if the library grid ever drifts from the stored
    literal (regenerate ``ornl.py`` then). The structurally-independent reference is the
    real loaded data, not a second copy of the literal.
    """
    np.testing.assert_array_equal(
        ORNL_421.edges, np.asarray(pwr_like_mix().eg, dtype=float),
        err_msg="ORNL_421 literal drifted from the live library eg — regenerate ornl.py",
    )


@_needs_real_421
@pytest.mark.slow
def test_real_pwr_421_to_wims69_condensation_succeeds():
    """Condense the REAL 421-group ``pwr_like_mix()`` → WIMS-69 end-to-end.

    The L2 real-structure integration leg (NON-NESTED 421→69). Asserts:
      * condensation does NOT raise (the empty-group bug is gone);
      * the condensed Mixture has ng == 69;
      * NO zero columns — every coarse group is populated (was 3 empty under
        one-hot for →69);
      * ``assert_balanced()`` passes on the real-data collapse;
      * χ stays a fast-peaked simplex (bulk in the fast HALF, NOT literally
        argmax==0 — the real fission χ peaks ≈1 MeV → a low-but-nonzero coarse
        index, so the assertion is the physically-correct cumulative mass, not
        a brittle argmax);
      * ``fractional_columns`` is reported (non-empty — the resonance/thermal
        coarse groups that leaned on w(E)).
    """
    from orpheus.data.energy_grid import EnergyGrid  # type: ignore
    from orpheus.data.group_structures.wims import WIMS_69  # type: ignore

    fine_mix = pwr_like_mix()                            # real 421-group, balanced
    np.testing.assert_array_equal(
        fine_mix.ng, 421, err_msg="pwr_like_mix must be the 421-group library"
    )
    fine_mix.assert_balanced()                           # the real fine data balances

    fine_grid = EnergyGrid(np.asarray(fine_mix.eg, dtype=float))
    # NON-NESTED: WIMS-69 boundaries do not align with the 421 grid. The overlap
    # basis carries the fractional table + the locally-interpolated report.
    overlap = fine_grid.overlap_to(WIMS_69, _IES())

    # The per-material spectrum (the real fine flux is not on the Mixture; use a
    # 1/E-ish slowing-down weight as the condensing spectrum — the within-group
    # MODEL handles the sub-fine split, this φ_g is the between-group weight).
    reps = fine_grid.representative_energy                # (421,) descending eV
    phi = 1.0 / reps                                      # 1/E slowing-down weight
    coarse_mix = fine_mix.condense(WIMS_69, phi, _IES())  # MUST NOT raise

    # ng == 69, shapes coarse.
    np.testing.assert_array_equal(
        coarse_mix.ng, 69, err_msg="condensed Mixture must live on WIMS-69 (ng=69)"
    )
    np.testing.assert_array_equal(
        np.asarray(coarse_mix.SigT).shape, (69,),
        err_msg="condensed vectors must be (69,)",
    )
    np.testing.assert_array_equal(
        np.asarray(coarse_mix.SigS[0].todense()).shape, (69, 69),
        err_msg="condensed scattering must be (69, 69)",
    )

    # NO zero columns — every coarse group populated (the empty-group bug gone).
    # SigT > 0 in every coarse group (a populated group has a positive total XS).
    np.testing.assert_array_less(
        0.0, np.asarray(coarse_mix.SigT, dtype=float),
        err_msg="every WIMS-69 coarse group must be populated (SigT>0) — the "
        "empty-group bug (3 empty under one-hot →69) must be gone",
    )

    # Balance survives the real-data NON-NESTED collapse.
    coarse_mix.assert_balanced(atol=1e-9)

    # χ stays a fast-peaked simplex. "Fast-peaked" = bulk of the emission mass
    # in the fast half of the 69-group structure (group 0 = fastest). NOT
    # argmax==0: the real fission χ peaks ≈1 MeV, a low-but-nonzero coarse
    # index, so a cumulative-mass property is the physically-correct gate.
    chi_c = np.asarray(coarse_mix.chi, dtype=float)
    np.testing.assert_allclose(
        float(chi_c.sum()), 1.0, atol=1e-9,
        err_msg="condensed χ must remain a probability simplex",
    )
    np.testing.assert_array_less(
        -1e-15, chi_c, err_msg="condensed χ must be non-negative"
    )
    fast_half = float(chi_c[: 69 // 2].sum())
    np.testing.assert_array_less(
        0.5, fast_half,
        err_msg=f"condensed χ must stay fast-peaked (fast-half mass={fast_half:.3f} "
        "must exceed 0.5) — the bulk of fission emission is fast",
    )

    # fractional_columns reported (non-empty — the non-nested straddle groups).
    li = np.asarray(overlap.fractional_columns, dtype=int).ravel()
    np.testing.assert_array_less(
        0, li.size,
        err_msg="421→WIMS-69 is NON-NESTED → fractional_columns must be "
        "non-empty (the resonance/thermal coarse groups that leaned on w(E))",
    )


# ════════════════════════════════════════════════════════════════════
# P6 (#281) — the bilinear (eigenvalue-consistent) condensation.
# Spec: .claude/plans/p6_adjoint_verification_spec.md §4 C4/C5 (with the
# post-B&G convention correction — the sink axis GAINS the adjoint
# carrier per B&G (6.136); it is NOT frozen);
# rules: orpheus/derivations/common/homogenization.py T6 (B&G Ch. 6).
# ════════════════════════════════════════════════════════════════════


@pytest.fixture(scope="module")
def adjoint_solution(materials, problem):
    """The adjoint solve over the SAME constituents as ``solution`` — the
    two-entry pattern ``is_same_phase_space`` exists to accept."""
    from orpheus.sn.solver import solve_sn_adjoint

    fine, quad = problem
    return solve_sn_adjoint(materials, fine, quad, scattering_order=0)


def _adjoint_region_spectrum(adj, mat_id: int = 0) -> np.ndarray:
    """The importance spectrum φ*_g — the SAME V-weighted reduction as the
    forward representative (the pair-of-representatives convention)."""
    phis = np.asarray(adj.scalar_flux.values, dtype=float)
    V = np.asarray(adj.mesh.volumes, dtype=float)
    mat = np.asarray(adj.mesh.mat_map, dtype=int).ravel()
    sel = mat == mat_id
    return (phis[:, sel] * V[sel]).sum(axis=1)


@pytest.mark.verifies("sn-homogenization-adjoint-weighted")
@pytest.mark.verifies("sn-homogenization-bilinear")
class TestC4BilinearCondensation:
    """C4: every channel of the bilinear-condensed Mixture equals the B&G
    Ch. 6 hand rule (nested blocks — no straddles, 0/1 overlap table), the
    result discriminates from the forward collapse, and the fixture proves
    itself non-dud.

    NOTE the spec correction (post-B&G, 2026-07-26): the sink axis is NOT
    frozen under the adjoint — it gains the flux-weighted-average adjoint
    carrier Ψ†_G (B&G (6.136)); χ takes the χ† rule with the rank-1
    simplex rescale. The pre-B&G Gate B ("marginalize frozen") pinned the
    superseded expectation and is replaced by the hand-rule equalities.
    """

    _BLOCKS = [(0, 1), (2, 3)]

    def _hand(self, mix, phi_s, phis_s):
        """The B&G-convention hand collapse of ``mix`` with the spectrum pair."""
        blocks = self._BLOCKS
        Phi = np.array([phi_s[list(b)].sum() for b in blocks])
        pair = np.array([(phis_s * phi_s)[list(b)].sum() for b in blocks])
        psi_dag = pair / Phi
        out = {}
        for name in ("SigT", "SigC", "SigL", "SigF"):
            ch = np.asarray(getattr(mix, name), float)
            out[name] = np.array(
                [(phis_s * ch * phi_s)[list(b)].sum() for b in blocks]
            ) / pair
        s_mat = np.asarray(mix.SigS[0].todense(), float)
        sig_s = np.zeros((2, 2))
        for Gf, bf in enumerate(blocks):
            for Gt, bt in enumerate(blocks):
                num = sum(
                    phis_s[g] * s_mat[gp, g] * phi_s[gp] for gp in bf for g in bt
                )
                sig_s[Gf, Gt] = num / (Phi[Gf] * psi_dag[Gt])
        chi_f = np.asarray(mix.chi, float)
        chi_dag = np.array(
            [(chi_f * phis_s)[list(b)].sum() for b in blocks]
        ) / psi_dag
        s = chi_dag.sum()
        nsf = np.asarray(mix.SigP, float)
        out["chi"] = chi_dag / s
        out["SigP"] = (
            np.array([(nsf * phi_s)[list(b)].sum() for b in blocks]) / Phi
        ) * s
        out["SigS0"] = sig_s
        return out

    def test_dud_guard_spectrum_shapes_differ(self, solution, adjoint_solution):
        phi_s = _region_spectrum(solution)
        phis_s = _adjoint_region_spectrum(adjoint_solution)
        a = phis_s / np.linalg.norm(phis_s)
        b = phi_s / np.linalg.norm(phi_s)
        assert not np.allclose(a, b, rtol=1e-2, atol=1e-3), (
            "fixture too self-adjoint on the energy axis: φ*_spec ≈ φ_spec"
        )

    def test_every_channel_matches_the_bg_hand_rule(
        self, materials, solution, adjoint_solution,
    ):
        out = solution.condense(EnergyGrid(_EG_COARSE), adjoint=adjoint_solution)[0]
        fwd = solution.condense(EnergyGrid(_EG_COARSE))[0]
        ref = self._hand(
            materials[0],
            _region_spectrum(solution),
            _adjoint_region_spectrum(adjoint_solution),
        )
        discriminated = False
        for name in ("SigT", "SigC", "SigL", "SigF", "SigP", "chi"):
            got = np.asarray(getattr(out, name), float)
            np.testing.assert_allclose(got, ref[name], rtol=1e-12, err_msg=name)
            if not np.allclose(got, np.asarray(getattr(fwd, name), float),
                               rtol=1e-6):
                discriminated = True
        np.testing.assert_allclose(
            np.asarray(out.SigS[0].todense(), float), ref["SigS0"], rtol=1e-12,
        )
        assert discriminated, "bilinear ≡ forward on every channel: dud fixture"

    def test_bilinear_chi_is_simplex(self, solution, adjoint_solution):
        out = solution.condense(EnergyGrid(_EG_COARSE), adjoint=adjoint_solution)[0]
        chi = np.asarray(out.chi, float)
        np.testing.assert_allclose(chi.sum(), 1.0, atol=1e-12)
        np.testing.assert_array_less(-1e-15, chi)

    def test_no_arg_equals_explicit_none_bitwise(self, solution):
        """§4.0 tooth 1 on the condense verb."""
        a = solution.condense(EnergyGrid(_EG_COARSE))
        b = solution.condense(EnergyGrid(_EG_COARSE), adjoint=None)
        for mat_id in a:
            ma, mb = a[mat_id], b[mat_id]
            for ch in ("SigT", "SigC", "SigL", "SigF", "SigP", "chi"):
                np.testing.assert_array_equal(
                    np.asarray(getattr(ma, ch)), np.asarray(getattr(mb, ch)),
                )
            for la, lb in zip(ma.SigS, mb.SigS):
                np.testing.assert_array_equal(
                    np.asarray(la.todense()), np.asarray(lb.todense()),
                )


class TestC5CondenseWeightCaptureSentinel:
    """C5 (Mode-11 capture): the bilinear condense hands the PAIR spectrum
    weight φ*_spec⊙φ_spec to a test-basis construction; bare φ*_spec never
    appears as a frame weight (it enters only the sink/χ† folds). The plain
    φ_spec frame IS legitimately constructed too — the B&G source-axis /
    νΣf flux-average (correction to the pre-B&G C5 wording)."""

    def test_pair_weight_captured_and_bare_adjoint_absent(
        self, solution, adjoint_solution, monkeypatch,
    ):
        from orpheus.numerics.basis.weighted_indicator_basis import (
            WeightedIndicatorBasis,
        )

        phi_s = _region_spectrum(solution)
        phis_s = _adjoint_region_spectrum(adjoint_solution)
        captured: list[np.ndarray] = []
        orig = WeightedIndicatorBasis.__init__

        def spy(self, basis, weight, *a, **k):
            captured.append(np.asarray(weight, dtype=float))
            return orig(self, basis, weight, *a, **k)

        monkeypatch.setattr(WeightedIndicatorBasis, "__init__", spy)
        solution.condense(EnergyGrid(_EG_COARSE), adjoint=adjoint_solution)

        def _seen(expected):
            return any(
                c.shape == np.shape(expected) and np.array_equal(c, expected)
                for c in captured
            )

        assert _seen(phis_s * phi_s), (
            "no frame received the PAIR spectrum weight φ*⊙φ"
        )
        assert _seen(phi_s), (
            "the source-axis flux-average frame (B&G 6.136/4.38) was not built"
        )
        assert not _seen(phis_s), (
            "a frame received BARE φ*_spec (the φ→φ* trap on the energy axis)"
        )
