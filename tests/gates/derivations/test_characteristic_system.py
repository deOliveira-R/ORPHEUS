"""Gates for the characteristic reference's multigroup Galerkin system (:class:`~orpheus.derivations.continuous.characteristic.system.GalerkinSystem`).

P1 step (b), fourth rung, of the characteristic-reference campaign
(``.claude/plans/characteristic_reference_architecture.md``, "P1 step (b), fourth
rung: API sketch", ruled 2026-10-07, and the user's ruling of the same day that
the unknown is the emission density q on the supports). Verification spec
``scratch/characteristic_architecture/p1_step_b4/spec.md``: rows D0a, D0b, SU1, SU2, SU3', SU4, SU5,
XS1, D7, D5, D5b, D5d, E2, E3, E6, D11a, D11b, D2', D13, D14.

The object: the Galerkin form W_s q = (S + F/k) K q + W_s q_ext of the emission
on each group's support, with the flux read from it, W_G phi = K q. The k pencil
is (W_s - S K, F K); the source questions read (W_s, (S + F) K).

The ladder, bottom up (``rests_on`` on each row):

1. rung 3, landed: each group's block K_g (``test_characteristic_assembly.py``);
2. term level: the cross sections and their per-node layout [D0a, D0b], the
   emission support [SU1], the anisotropy refusal [XS1];
3. structure: the refusals off the support and of a mis-oriented mask [SU2, SU3'], the pencil residual [D7];
4. closed bodies against the 0-D pair written by hand [D5, D5b, D5d, E2], the
   support's two edges [SU4, SU5];
5. open bodies: the 1-group Rayleigh-Ritz rows [D11a, D11b], Sood's 2-group
   critical sizes [D2'], the source questions [E3, E6];
6. the adjoint [D13, D14].

Every closed form is written here from the ``Mixture`` arrays (``toarray``,
``np.outer``, ``np.linalg``), never through ``group_emission`` or the 0-D
helpers of ``orpheus.derivations.common.eigenvalue`` (``instrument-doctrine`` X4).
Bands are measured (`[M]` 2026-10-07, ``.venv/bin/python -O``, the spec's probes
``scratch/characteristic_architecture/p1_step_b4/ta/``); each row's first red is a
battery arm (``scratch/characteristic_architecture/p1_step_b4/battery/``).

Labels (minted 2026-10-07, ``docs/theory/references/characteristic.rst``, section
``characteristic-galerkin-system``): ``characteristic-pencil`` (D5/D5b, D5d, D2'),
``characteristic-adjoint`` (D13, D14 (i)), ``characteristic-fixed-source`` (E2, E3,
D14 (ii)), ``characteristic-one-group-bound`` (D11a, D11b). The rows on the
cylinder are ``slow`` (one-region cylinder, 26 s for two groups at 8 points).
"""
from __future__ import annotations

from decimal import Decimal
from functools import lru_cache

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.dense_pencil import NoLeastSolution
from orpheus.derivations.continuous.characteristic import (
    GalerkinSystem, PanelBasis, RegionCrossSections, TransportResolution, Wall, Walls,
)
from orpheus.derivations.continuous.sood_registry import sood2003 as sood
from orpheus.geometry.chart import Chart
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.coord import CoordSystem

pytestmark = pytest.mark.filterwarnings("error::RuntimeWarning")

_HERE = "tests/gates/derivations/test_characteristic_system.py::"
_ASSEMBLY = "tests/gates/derivations/test_characteristic_assembly.py::"
_CONSERVES = _ASSEMBLY + "test_a_closed_body_conserves_its_emission"
_OWN_GROUP = _ASSEMBLY + "test_each_group_is_coupled_through_its_own_cross_section"
_ESCAPE = _ASSEMBLY + "test_the_escape_and_transmission_probabilities_are_the_closed_forms"

_COORD = {"sphere": CoordSystem.SPHERICAL, "cylinder": CoordSystem.CYLINDRICAL, "slab": CoordSystem.CARTESIAN}
_MR3 = (0.0, 0.5, 1.5, 2.0)
_MR3H = (0.4, 0.5, 1.5, 2.0)
_SLB3 = (0.0, 0.4, 1.5, 2.3)
_MIRROR, _WHITE, _VACUUM = (1.0, 0.0), (0.0, 1.0), (0.0, 0.0)
_WORK = (3, 2, 0.4, 8, 12)          # (degree, layers, ratio, line points, points along each line): the working point


# ── the mixtures, built directly (``make_mixture`` nulls Sig2 unless given) ──


def _mixture(sig_t, sig_s, sig_2=None, nu_sig_f=None, chi=None) -> Mixture:
    """A Mixture from (from, to) tables; capture is the balance remainder, nu = 2.5."""
    sig_t = np.asarray(sig_t, dtype=float)
    groups = sig_t.size
    s = np.asarray(sig_s, dtype=float)
    s2 = np.zeros((groups, groups)) if sig_2 is None else np.asarray(sig_2, dtype=float)
    nsf = np.zeros(groups) if nu_sig_f is None else np.asarray(nu_sig_f, dtype=float)
    sig_f = nsf / 2.5
    capture = sig_t - s.sum(axis=1) - s2.sum(axis=1) - sig_f        # the balance Sigma_t = c + f + out-scatter + n2n
    if capture.min() < 0.0:
        raise ValueError(f"an inconsistent fixture: negative capture {capture}")
    return Mixture(
        SigC=capture, SigL=np.zeros(groups), SigF=sig_f, SigP=nsf, SigT=sig_t,
        SigS=[csr_matrix(s)], Sig2=[csr_matrix(s2)], chi=np.zeros(groups) if chi is None else np.asarray(chi, dtype=float),
    )


#: Upscatter only, an (n,2n) transfer in both directions, subcritical fission with chi in both groups: k_inf 0.344.
_UP2N = _mixture([1.0, 1.5], [[0.5, 0.0], [0.1, 1.2]], [[0.03, 0.02], [0.0, 0.0]], [0.05, 0.2], [0.7, 0.3])
#: A non-fissile absorber with downscatter.
_ABS = _mixture([0.8, 1.2], [[0.4, 0.2], [0.0, 0.7]])
#: Supercritical by its (n,2n) emission alone, no fission: the closed body's collision gain has radius 1.155.
_N2N = _mixture([1.0, 1.5], [[0.3, 0.1], [0.05, 0.9]], [[0.35, 0.0], [0.0, 0.4]])
#: Transparent in group 1 (Sigma_t = 0), scattering into it from group 0.
_TRANSPARENT = _mixture([1.0, 0.0], [[0.5, 0.3], [0.0, 0.0]])
#: Absorbs in both groups, scatters into both.
_ABSORBER = _mixture([0.8, 1.2], [[0.3, 0.2], [0.0, 0.6]])
#: Void in group 1, emitting in group 0 only.
_SELF_ONLY = _mixture([0.9, 0.0], [[0.6, 0.0], [0.0, 0.0]])
#: Absorbs in group 1 and scatters OUT of it (upscatter), nothing INTO it.
_OUT_ONLY = _mixture([0.5, 0.6], [[0.2, 0.0], [0.3, 0.0]])
#: Fissile with chi = (1, 0) and downscatter, no slow self-scatter.
_FAST_CHI = _mixture([0.7, 1.1], [[0.3, 0.1], [0.0, 0.0]], None, [0.1, 0.5], [1.0, 0.0])

_PU2 = sood.PU_2_0_SL_STUB.materials[0]       # downscatter, chi in both groups
_URRB = sood.URRB_2_0_IN.materials[0]         # upscatter, chi = (1, 0)
_URRC = sood.URRC_2_0_IN.materials[0]
_URR3 = sood.URR_3_0_IN.materials[0]          # three groups


def _scattering_by_hand(m: Mixture) -> np.ndarray:
    """(Sigma_s + 2 Sigma_2)^T, [to, from]."""
    return (m.SigS[0].toarray() + 2.0 * m.Sig2[0].toarray()).T


def _fission_by_hand(m: Mixture) -> np.ndarray:
    """chi (x) nu Sigma_f, [to, from]."""
    return np.outer(m.chi, m.SigP)


def _zero_d(m: Mixture) -> tuple[np.ndarray, np.ndarray]:
    """The infinite medium's loss diag(Sigma_t) - S and production F, by hand."""
    return np.diag(np.asarray(m.SigT, dtype=float)) - _scattering_by_hand(m), _fission_by_hand(m)


def _dominant(matrix: np.ndarray) -> tuple[float, np.ndarray]:
    """The dominant eigenpair of a small positive matrix, the vector normalised to sum 1."""
    values, vectors = np.linalg.eig(matrix)
    i = int(np.argmax(np.abs(values)))
    v = np.abs(np.real(vectors[:, i]))
    return float(values[i].real), v / v.sum()


# ── the bodies ───────────────────────────────────────────────────────────


def _walls(chart: str, breakpoints, laws) -> Walls:
    """One (specular, diffuse) per boundary point, inner first, or ``"P"`` for a periodic face."""
    n = len(breakpoints) - 1
    ends = (n,) if chart != "slab" and breakpoints[0] == 0.0 else (0, n)
    walls = tuple(Wall(k, 1.0, 0.0, n - k) if law == "P" else Wall(k, law[0], law[1], k)
                  for k, law in zip(ends, laws, strict=True))
    return Walls(walls, n, Chart(_COORD[chart]))


@lru_cache(maxsize=None)
def _basis(chart: str, breakpoints, degree: int, layers: int, ratio: float) -> PanelBasis:
    return PanelBasis.of(ConcentricPartition(Chart(_COORD[chart]), tuple(breakpoints)), degree, layers, ratio)


@lru_cache(maxsize=None)
def _system(chart: str, breakpoints, laws, mixtures, resolution=_WORK, source_regions=None) -> GalerkinSystem:
    degree, layers, ratio, line_points, points = resolution
    return GalerkinSystem(
        _basis(chart, breakpoints, degree, layers, ratio), _walls(chart, breakpoints, laws),
        RegionCrossSections.of(list(mixtures)), TransportResolution(line_points, points, points),
        source_regions=None if source_regions is None else np.array(source_regions, dtype=bool),
    )


#: Closed bodies: every wall returns everything (mirror, white at 1, the periodic wrap).
_CLOSED = [
    ("sphere", (0.0, 1.0), (_MIRROR,)),
    ("sphere", _MR3, (_MIRROR,)),
    ("sphere", (0.0, 1.0), (_WHITE,)),
    ("slab", _SLB3, (_MIRROR, _MIRROR)),
    ("slab", (0.0, 1.0), ("P", "P")),
    ("sphere", _MR3H, (_MIRROR, _MIRROR)),
    ("sphere", (0.4, 2.0), (_WHITE, _MIRROR)),
]
_CLOSED_IDS = ["sphere1-mirror", "sphere3-mirror", "sphere1-white", "slab3-mirrors", "slab1-periodic",
               "hollow3-mirrors", "hollow1-white-mirror"]
_CLOSED_MIXTURES = {"PU2": _PU2, "URRb": _URRB, "URR3": _URR3}


# ── 4.0 term level ───────────────────────────────────────────────────────


@pytest.mark.l0
def test_the_cross_sections_are_the_mixtures_tables() -> None:
    """[D0a; l0] ``total = SigT``, ``scattering = (SigS0 + 2 Sig2_0)^T``,
    ``fission = chi (x) nu Sigma_f``, region by region, ``array_equal``.

    The same three float operations as the assembly site (one add, one
    transpose, one outer product), so the comparison is bitwise. The (n,2n)
    leg is live only because ``_UP2N`` carries a non-zero Sig2: every Sood and
    ``xs_library`` mixture ships Sig2 = 0. First reds: scattering not
    transposed; the (n,2n) multiplicity 1; ``outer(nu Sigma_f, chi)``.
    """
    mixtures = [_UP2N, _PU2, _ABS, _FAST_CHI]
    xs = RegionCrossSections.of(mixtures)
    assert xs.n_groups == 2
    for r, m in enumerate(mixtures):
        np.testing.assert_array_equal(xs.total[r], m.SigT)
        np.testing.assert_array_equal(xs.scattering[r], _scattering_by_hand(m))
        np.testing.assert_array_equal(xs.fission[r], _fission_by_hand(m))
    # the premises the legs need: asymmetric transfer, a non-zero (n,2n), asymmetric chi (x) nu Sigma_f
    assert not np.array_equal(xs.scattering[0], xs.scattering[0].T)
    assert _UP2N.Sig2[0].count_nonzero() > 0
    assert not np.array_equal(xs.fission[0], xs.fission[0].T)


@pytest.mark.l0
@pytest.mark.rests_on(_HERE + "test_the_cross_sections_are_the_mixtures_tables",
                      _HERE + "test_the_emission_support_is_where_something_is_emitted_into_the_group")
def test_the_emission_acts_node_by_node() -> None:
    """[D0b; l0] Row (g, a) of ``scattering`` and ``fission`` (node i = support_g[a]) is zero except at the
    columns (g', i), where it is the region's [g, g'] entry; ``array_equal`` against a table built here node by node.

    The supports differ by group (``_SELF_ONLY``'s region emits in group 0
    only), so a row offset error shows. First reds: every group's column read
    at group 0 (the group offset dropped); the per-node matrix transposed
    ([r, g', g]).
    """
    system = _system("sphere", _MR3, (_VACUUM,), (_UP2N, _SELF_ONLY, _FAST_CHI), (1, 0, 0.5, 4, 4))
    xs, basis = system.cross_sections, system.basis
    n = basis.size
    supports = system.emission.supports
    assert [s.size for s in supports] == [3 * 2, 2 * 2]                     # the supports differ by group
    for name in ("scattering", "fission"):
        table = getattr(xs, name)
        expected = np.zeros((system.emission.size, 2 * n))
        row = 0
        for g, support in enumerate(supports):
            for i in support:
                for g_from in range(2):
                    expected[row, g_from * n + i] = table[basis.region[i], g, g_from]
                row += 1
        np.testing.assert_array_equal(getattr(system, name), expected)


@pytest.mark.foundation
def test_the_emission_support_is_where_something_is_emitted_into_the_group() -> None:
    """[SU1] ``emission_support()`` is the hand table, region-major ``(n, G)``: True where a scattering, (n,2n) or fission
    transfer enters group g in that region.

    Four regions, one per edge: transparent in group 1 but scattering into it
    (True, True); emitting in group 0 only (True, False); scattering OUT of
    group 1 with nothing into it (True, False); fissile with chi = (1, 0) and a
    downscatter (True, True). First reds: the support read from Sigma_t > 0
    (the first two regions flip); the union over groups (the second and third
    flip); the "from" axis read for the "to" axis (the third flips).
    """
    xs = RegionCrossSections.of([_TRANSPARENT, _SELF_ONLY, _OUT_ONLY, _FAST_CHI])
    expected = np.array([[True, True], [True, False], [True, False], [True, True]])
    np.testing.assert_array_equal(xs.emission_support(), expected)


@pytest.mark.foundation
@pytest.mark.parametrize(("stack", "fragment"), [
    ("SigS", "region 1's SigS has a non-zero Legendre order 1"),
    ("Sig2", "region 1's Sig2 has a non-zero Legendre order 1"),
])
def test_an_anisotropic_emission_is_refused_naming_the_region_and_the_order(stack: str, fragment: str) -> None:
    """[XS1] SCOPE-BOUNDARY refusal: a non-zero Legendre order >= 1 of SigS or Sig2, one leg per stack; a stack
    carrying an all-zero higher block is admitted (positive leg).

    First reds: the guard reads SigS only (the Sig2 leg); the guard reads the
    stack's length (the positive leg).
    """
    p1 = csr_matrix(np.array([[0.05, 0.0], [0.0, 0.02]]))
    zero = csr_matrix((2, 2))
    stacks = {"SigS": [_ABS.SigS[0], zero], "Sig2": [_ABS.Sig2[0], zero]}
    admitted = Mixture(SigC=_ABS.SigC, SigL=_ABS.SigL, SigF=_ABS.SigF, SigP=_ABS.SigP, SigT=_ABS.SigT,
                       SigS=stacks["SigS"], Sig2=stacks["Sig2"], chi=_ABS.chi)
    RegionCrossSections.of([_ABS, admitted])
    stacks[stack] = [stacks[stack][0], p1]
    refused = Mixture(SigC=_ABS.SigC, SigL=_ABS.SigL, SigF=_ABS.SigF, SigP=_ABS.SigP, SigT=_ABS.SigT,
                      SigS=stacks["SigS"], Sig2=stacks["Sig2"], chi=_ABS.chi)
    with pytest.raises(NotImplementedError, match=fragment):
        RegionCrossSections.of([_ABS, refused])


# ── 4.1 structure ────────────────────────────────────────────────────────


_SUPPORT_BODY = (_TRANSPARENT, _ABSORBER, _SELF_ONLY)


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_emission_support_is_where_something_is_emitted_into_the_group")
def test_a_field_off_the_support_is_refused_and_answered_once_the_support_is_widened() -> None:
    """[SU2] ``emission.restrict`` (and so ``fixed_source``) refuses a group-1 source in the layer where group 1 emits
    nothing, naming the group; posed with ``source_regions`` the same source is answered, finite and
    non-negative.

    The slab between mirrors: under the mirror sphere the widened column is a
    source on lossless trapped lines (``TrappedSource``, SU5). First red:
    ``restrict`` drops the off-support coefficients silently.
    """
    body = ("slab", (0.0, 0.5, 1.5, 2.1), (_MIRROR, _MIRROR), _SUPPORT_BODY)
    system = _system(*body)
    source = np.zeros((2, system.basis.size))
    source[1, system.basis.region == 2] = 1.0
    with pytest.raises(ValueError, match=r"off the emission support, in the regions \{1: \[2\]\} \(by group\)"):
        system.fixed_source(source)
    widened = _system(*body, source_regions=((False, False), (False, False), (False, True)))
    flux = widened.fixed_source(source)
    assert np.all(np.isfinite(flux)) and flux.min() >= 0.0 and flux[1].max() > 0.0


@pytest.mark.foundation
def test_a_source_region_mask_in_group_major_orientation_is_refused() -> None:
    """[SU3'] ``source_regions`` is region-major, ``(n, G)`` as the cross sections: a ``(G, n)`` mask is refused,
    naming the orientation, when n != G (three regions, two groups).

    The refusal is the orientation's only witness: the rung's first spelling
    was group-major, and with n == G a transposed mask has the right shape
    and is silently read as another set of regions (no shape can see it).
    First red: the shape check removed (the mask then fails to broadcast,
    with numpy's message, or is read transposed when n == G).
    """
    basis = _basis("sphere", _MR3, 1, 0, 0.5)
    walls = _walls("sphere", _MR3, (_VACUUM,))
    xs = RegionCrossSections.of([_UP2N, _ABS, _PU2])
    with pytest.raises(ValueError, match=r"region-major \(3, 2\) mask"):
        GalerkinSystem(basis, walls, xs, TransportResolution(4, 4, 4), source_regions=np.zeros((2, 3), dtype=bool))


_RESIDUAL_CASES = [
    ("sphere", _MR3, (_MIRROR,), (_URRB,) * 3),
    ("slab", _SLB3, (_MIRROR, _VACUUM), (_UP2N, _ABS, _PU2)),
    ("sphere", _MR3, (_VACUUM,), (_UP2N, _ABS, _PU2)),
]


@pytest.mark.foundation
@pytest.mark.parametrize("case", _RESIDUAL_CASES, ids=["closed-sphere", "slab", "sphere"])
@pytest.mark.rests_on(_CONSERVES)
def test_the_fundamental_mode_satisfies_its_pencil(case) -> None:
    """[D7] ``||L q - P q / k|| / ||P q / k||`` below 1e-12 for the k pencil's fundamental on the emission space.

    `[M]` <= 1.4e-14 over 28 closed bodies. First red: the mode's vector
    returned in another group order.
    """
    pencil = _system(*case).pencil
    mode = pencil.fundamental()
    produced = pencil.production @ mode.vector / mode.k
    residual = pencil.loss @ mode.vector - produced
    assert np.linalg.norm(residual) / np.linalg.norm(produced) < 1e-12


# ── 4.2 closed bodies ────────────────────────────────────────────────────


def _closed_params(bodies, ids, mixtures, marks=()):
    return [pytest.param(*body, name, id=f"{bid}-{name}", marks=marks)
            for body, bid in zip(bodies, ids, strict=True) for name in mixtures]


_CLOSED_SLOW = [("cylinder", (0.0, 1.0), (_MIRROR,)), ("cylinder", (0.0, 1.0), (_WHITE,))]


@pytest.mark.l1
@pytest.mark.verifies("characteristic-pencil")
@pytest.mark.parametrize(("chart", "breakpoints", "laws", "name"),
                         _closed_params(_CLOSED, _CLOSED_IDS, _CLOSED_MIXTURES)
                         + _closed_params(_CLOSED_SLOW, ["cylinder-mirror", "cylinder-white"], ["URRb"],
                                          marks=pytest.mark.slow))
@pytest.mark.rests_on(_CONSERVES, _OWN_GROUP, _HERE + "test_the_emission_acts_node_by_node",
                      _HERE + "test_the_fundamental_mode_satisfies_its_pencil")
def test_a_closed_body_reads_k_inf_with_a_flat_flux_in_the_infinite_mediums_group_ratio(chart, breakpoints, laws, name) -> None:
    """[D5, D5b; l1, ``characteristic-pencil``] A closed body of one mixture: k = rho(A^-1 F) to 1e-12
    relative; the flux (``flux`` of the fundamental emission) is flat to 1e-10 in each group with the infinite
    medium's group ratio to 1e-12. The flatness band holds an 11x margin over the
    largest reading on these fixtures (8.8e-12, ``ta/m10_bands.log``); qa
    measured 1.14e-10 on a closer fixture (two regions, R = 1.5) at this
    resolution, so the band is a property of THESE bodies, not of the method.

    The 0-D pair is written here from the mixture's arrays. `[M]` 2026-10-07 at
    the working point: k <= 4.4e-14, flatness <= 8.8e-12 (the white/mirror
    hollow sphere), ratio <= 1.9e-14; the cylinder 2.3e-13 and 9.0e-13. PU2
    has downscatter and chi in both groups, URRb upscatter and chi = (1, 0),
    URR3 three groups. First reds: the scattering transposed (k and the ratio
    move O(1)); the group blocks paired with another group's emission.
    """
    mixture = _CLOSED_MIXTURES[name] if name in _CLOSED_MIXTURES else _URRB
    system = _system(chart, breakpoints, laws, (mixture,) * (len(breakpoints) - 1))
    loss, production = _zero_d(mixture)
    k_inf, ratio = _dominant(np.linalg.solve(loss, production))
    mode = system.pencil.fundamental()
    assert abs(mode.k - k_inf) / k_inf < 1e-12, (mode.k, k_inf)
    flux = system.flux(mode.vector)
    mean = flux.mean(axis=1)
    assert np.max(np.abs(flux / mean[:, None] - 1.0)) < 1e-10
    np.testing.assert_allclose(mean / mean.sum(), ratio, rtol=0.0, atol=1e-12)


@pytest.mark.l1
@pytest.mark.verifies("characteristic-pencil")
@pytest.mark.parametrize(("case", "published"), [(sood.URRB_2_0_IN, "1.365821"), (sood.URRC_2_0_IN, "1.633380")],
                         ids=["URRb-2-0-IN", "URRc-2-0-IN"])
@pytest.mark.rests_on(_HERE + "test_a_closed_body_reads_k_inf_with_a_flat_flux_in_the_infinite_mediums_group_ratio")
def test_a_closed_sphere_reads_soods_k_inf_with_upscatter(case, published: str) -> None:
    """[D5d; l1, ``characteristic-pencil``] Sood 2003 pp. 85-86: k_inf with thermal upscatter, within half a
    unit of the printed last digit (5e-7). `[M]` 4.1e-7 (URRb), 6.6e-8 (URRc).

    The truth is typed here from the report, not read from the registry's
    ``truth`` field (whose value would move with an edit of the registry).
    First red: the scattering transposed (the upscatter moves to the
    downscatter slot).
    """
    assert float(case.truth.k_eff_or_kinf) == float(published)
    system = _system("sphere", (0.0, 1.0), (_MIRROR,), (case.materials[0],))
    assert abs(system.pencil.fundamental().k - float(published)) < 5e-7


_E2_BODIES = [("sphere", _MR3, (_MIRROR,)), ("slab", _SLB3, (_MIRROR, _MIRROR)), ("sphere", (0.4, 2.0), (_WHITE, _MIRROR))]


@pytest.mark.l1
@pytest.mark.verifies("characteristic-fixed-source")
@pytest.mark.parametrize("body", _E2_BODIES, ids=["sphere3-mirror", "slab3-mirrors", "hollow-white-mirror"])
@pytest.mark.rests_on(_HERE + "test_a_closed_body_reads_k_inf_with_a_flat_flux_in_the_infinite_mediums_group_ratio")
def test_a_closed_body_with_a_uniform_source_reads_the_infinite_medium_flux(body) -> None:
    """[E2; l1, ``characteristic-fixed-source``] phi = (diag Sigma_t - S - F)^-1 q, flat in each group, with
    upscatter, (n,2n) and subcritical fission all active (``_UP2N``: collision gain radius 0.872), to 1e-10.

    `[M]` 3.5e-12, 1.7e-14, 8.4e-12 (``ta/m10_bands.log``): the band holds a 12x
    margin on these bodies only (see D5b). First reds: the scattering transposed;
    the (n,2n) multiplicity 1; the fission dropped from the source pencil's
    gain.
    """
    chart, breakpoints, laws = body
    system = _system(chart, breakpoints, laws, (_UP2N,) * (len(breakpoints) - 1))
    q = np.array([1.0, 0.4])
    loss, production = _zero_d(_UP2N)
    expected = np.linalg.solve(loss - production, q)
    flux = system.fixed_source(np.outer(q, np.ones(system.basis.size)))
    assert np.max(np.abs(flux / expected[:, None] - 1.0)) < 1e-10


def _balance(system: GalerkinSystem, source: np.ndarray, flux: np.ndarray, mixtures) -> float:
    """sum_g <1, Sigma_a,g phi_g>_W / sum_g <1, q_g>_W - 1, Sigma_a written by hand from the mixtures."""
    absorption = np.array([np.asarray(m.SigT) - m.SigS[0].toarray().sum(axis=1) for m in mixtures])   # (n, G)
    mass, one = system.basis.mass, np.ones(system.basis.size)
    absorbed = sum(one @ mass @ (absorption[system.basis.region, g] * flux[g]) for g in range(flux.shape[0]))
    emitted = sum(one @ mass @ source[g] for g in range(flux.shape[0]))
    return float(absorbed / emitted - 1.0)


@pytest.mark.l1
@pytest.mark.parametrize("body", [("sphere", (_MIRROR,)), ("slab", (_MIRROR, _MIRROR)), ("sphere", (_WHITE,))],
                         ids=["sphere-mirror", "slab-mirrors", "sphere-white"])
@pytest.mark.rests_on(_HERE + "test_the_emission_support_is_where_something_is_emitted_into_the_group",
                      _HERE + "test_a_closed_body_with_a_uniform_source_reads_the_infinite_medium_flux")
def test_a_region_transparent_in_a_group_keeps_its_emission_into_it(body) -> None:
    """[SU4; l1] A closed non-multiplying body whose inner region is transparent in group 1 and scatters into
    it: the absorption equals the source, to 1e-12, and group 1's flux is positive in the transparent region.

    `[M]` 4.9e-15, 2.2e-16, 5.1e-15. First red: the support read from
    Sigma_t > 0 (-2.2e-2 on the sphere, -2.0e-1 on the slab). A balance is
    blind to a redistribution that conserves (``vv`` anti-#8); D0b and D5b
    carry that.
    """
    chart, laws = body
    mixtures = _SUPPORT_BODY[:2]
    system = _system(chart, (0.0, 0.5, 1.5), laws, mixtures)
    source = np.zeros((2, system.basis.size))
    source[0] = 1.0
    flux = system.fixed_source(source)
    assert abs(_balance(system, source, flux, mixtures)) < 1e-12
    assert flux[1][system.basis.region == 0].min() > 0.1


@pytest.mark.l1
@pytest.mark.rests_on(_HERE + "test_a_region_transparent_in_a_group_keeps_its_emission_into_it")
def test_a_layer_void_in_a_group_under_a_mirror_is_answered() -> None:
    """[SU5; l1] A sphere whose outer layer is void in group 1 and emits in group 0 only, under a mirror: the
    system is built (no ``TrappedSource``), group 1's support excludes the layer, and the body balances to 1e-12.

    Only the mirror sphere discriminates: lines with an impact parameter
    beyond the material are trapped in the layer, lossless in group 1. On the
    slab and the white sphere a column there is not refused (`[M]` 2026-10-07,
    ``m4_support``). First red: the support taken as the union over groups,
    which puts a source on those lines (``TrappedSource``).
    """
    system = _system("sphere", (0.0, 0.5, 1.5, 2.1), (_MIRROR,), _SUPPORT_BODY)
    np.testing.assert_array_equal(system.cross_sections.emission_support(),
                                  [[True, True], [True, True], [True, False]])
    source = np.zeros((2, system.basis.size))
    source[0] = 1.0
    flux = system.fixed_source(source)
    assert abs(_balance(system, source, flux, _SUPPORT_BODY)) < 1e-12


# ── 4.3 open bodies ──────────────────────────────────────────────────────


def _one_group(sig_t: float, sig_s: float, nu_sig_f: float) -> Mixture:
    return _mixture([sig_t], [[sig_s]], None, [nu_sig_f], [1.0] if nu_sig_f else [0.0])


#: 1G heterogeneous: two fissile regions round a scatterer, Sigma_t distinct.
_ONE_GROUP = (_one_group(0.6, 0.3, 0.45), _one_group(1.3, 1.1, 0.0), _one_group(0.45, 0.2, 0.35))
_RR_BODIES = [("sphere", _MR3, (_VACUUM,)), ("slab", _SLB3, (_MIRROR, _VACUUM)), ("sphere", _MR3H, (_VACUUM, _VACUUM))]


@pytest.mark.l1
@pytest.mark.verifies("characteristic-one-group-bound")
@pytest.mark.parametrize("body", _RR_BODIES, ids=["sphere", "slab", "hollow-sphere"])
@pytest.mark.rests_on(_HERE + "test_a_closed_body_reads_k_inf_with_a_flat_flux_in_the_infinite_mediums_group_ratio")
def test_one_group_k_increases_on_nested_spaces(body) -> None:
    """[D11a] The 1-group Galerkin k increases strictly on nested spaces: degree p = 1, 2, 3, 4, 6 on fixed panels
    (the even panel nests too: polynomials in c^2), each step above 1e-10.

    The theorem: on the emission support the 1G form is the Rayleigh-Ritz form
    of the self-adjoint positive operator in the weight M^-1 (the plan's
    derivation), so k is a lower bound monotone under nesting. `[M]`
    2026-10-07: smallest step 2.8e-8 (hollow sphere, 3 to 4); the line rule's
    floor k(8,12) - k(16,16) <= 1.7e-13. First red: a non-reciprocal transport
    matrix (K perturbed by 1 + 0.3 U, the battery's A20), which leaves the form
    non-symmetric and the ladder non-monotone. Declared blind (battery A15,
    `[M]` 2026-10-07): an under-integrated line rule (4 points per piece, 8
    along each line) moves k_4 by +1.2e-5 but upward with the degree, so the
    order holds; the under-integration is caught by D5 and E2 instead.
    """
    chart, breakpoints, laws = body
    ks = [_system(chart, breakpoints, laws, _ONE_GROUP, (p, 2, 0.4, 8, 12)).pencil.fundamental().k
          for p in (1, 2, 3, 4, 6)]
    assert np.all(np.diff(ks) > 1e-10), ks


@pytest.mark.l1
@pytest.mark.verifies("characteristic-one-group-bound")
@pytest.mark.parametrize(("case", "chart", "laws"), [
    (sood.UA_1_0_SP_STUB, "sphere", (_VACUUM,)),
    (sood.UA_1_0_SL_STUB, "slab", (_MIRROR, _VACUUM)),
], ids=["Ua-1-0-SP", "Ua-1-0-SL"])
@pytest.mark.rests_on(_HERE + "test_one_group_k_increases_on_nested_spaces")
def test_one_group_k_is_below_one_at_soods_critical_size(case, chart: str, laws) -> None:
    """[D11b] The Rayleigh-Ritz bound against an independent truth: k_p < 1 for p = 1..4 on one ungraded panel at
    Sood's 1G critical size, with the margin 1 - k_p above the truth's resolution.

    Sood 2003 Table 10: Ua-1-0-SP 2.4248249802 mfp, Ua-1-0-SL 0.93772556 (the
    half thickness, mirror at the centre). The truth resolves k to
    |dk/dmfp| x half a unit of its last digit: 1.7e-11 (sphere), 3.6e-9
    (slab). `[M]` the smallest margin 2.4e-6 (sphere p = 4), 3.8e-6 (slab).
    On the graded working basis the slab's p = 4 reads +4.3e-10, inside the
    truth's resolution: no margin, so the row is posed on the ungraded panel.
    First red: a non-reciprocal transport matrix (battery A20). Declared blind:
    the under-integrated line rule (A15) leaves every k_p below 1 here.
    """
    mixture = case.materials[0]
    a = case.truth.critical_dimension_mfp / float(mixture.SigT[0])
    resolution = _last_digit(case.truth.critical_dimension_mfp) * 0.5 * 1.0       # |dk/dmfp| < 1 for both rows
    for p in (1, 2, 3, 4):
        k = _system(chart, (0.0, a), laws, (mixture,), (p, 0, 0.5, 8, 12)).pencil.fundamental().k
        assert 1.0 - k > 10.0 * resolution, (p, k)


def _last_digit(value: float) -> float:
    """The unit of the last printed digit of a registry value."""
    return 10.0 ** int(Decimal(repr(value)).as_tuple().exponent)


def _critical(case, resolution):
    """The case's body at its published critical size: full slab between vacuum walls, or the vacuum sphere."""
    mixture = case.materials[0]
    a = case.truth.critical_dimension_mfp / float(mixture.SigT[0])
    if case.geometry_kind == "slab":
        return ("slab", (0.0, a), (_MIRROR, _VACUUM), (mixture,), resolution), mixture
    return ("sphere", (0.0, a), (_VACUUM,), (mixture,), resolution), mixture


_SOOD_2G = [
    (sood.PU_2_0_SL_STUB, _WORK), (sood.PU_2_0_SP_STUB, _WORK), (sood.U_2_0_SL_STUB, _WORK),
    (sood.U_2_0_SP_STUB, _WORK), (sood.URRA_2_0_SL_STUB, (4, 3, 0.4, 12, 16)),
    (sood.URRA_2_0_SP_STUB, (4, 3, 0.4, 12, 16)),
]


@pytest.mark.l1
@pytest.mark.verifies("characteristic-pencil")
@pytest.mark.parametrize(("case", "resolution"), _SOOD_2G, ids=[c.case_id for c, _ in _SOOD_2G])
@pytest.mark.rests_on(_HERE + "test_a_closed_body_reads_k_inf_with_a_flat_flux_in_the_infinite_mediums_group_ratio",
                      _ESCAPE)
def test_k_is_one_at_soods_two_group_critical_sizes(case, resolution) -> None:
    """[D2'; l1, ``characteristic-pencil``] Sood 2003 Tables 29-38 (Siewert-Thomas 2G F_N): k = 1 within
    |dk/dmfp| x half a unit of the printed critical dimension's last digit, the slope measured here.

    `[M]` (k - 1 of band): PU-SL 3.4e-7 of 8.4e-7; PU-SP 7.0e-7 of 3.6e-6;
    U-SL 1.2e-7 of 5.1e-7; U-SP 1.9e-8 of 2.4e-6; URRa-SL 6.4e-7 of 7.2e-7;
    URRa-SP 2.2e-6 of 3.6e-6. URRa needs (4, 3, 0.4 / 12, 16): its slow group
    is 40 mfp thick. UAL-2-0-SL/SP are excluded: the reference converges to
    7.4e-6 and 9.4e-6 off Sood's digits, 2.9 and 7.5 half-units, and a
    production S_N ladder lands on the reference for the slab (the spec's
    "Refuted premise (D2' for UAL)"; an issue on the registry truth). First
    reds: the scattering transposed; chi applied to the wrong group.
    """
    body, mixture = _critical(case, resolution)
    k = _system(*body).pencil.fundamental().k
    chart, (lo, hi), laws, mixtures, res = body
    stretched = _system(chart, (lo, hi * (1.0 + 1e-4)), laws, mixtures, res).pencil.fundamental().k
    slope = (stretched - k) / (hi * 1e-4 * float(mixture.SigT[0]))           # per mfp of the published dimension
    band = abs(slope) * 0.5 * _last_digit(case.truth.critical_dimension_mfp)
    assert abs(k - 1.0) < band, (k - 1.0, band)


_SUBCRITICAL = [("sphere", _MR3, (_VACUUM,)), ("slab", _SLB3, (_VACUUM, _MIRROR))]


@pytest.mark.l1
@pytest.mark.verifies("characteristic-fixed-source")
@pytest.mark.parametrize("body", _SUBCRITICAL, ids=["sphere", "slab"])
@pytest.mark.rests_on(_HERE + "test_the_fundamental_mode_satisfies_its_pencil",
                      _HERE + "test_a_closed_body_with_a_uniform_source_reads_the_infinite_medium_flux")
def test_the_source_question_returns_the_mode_from_its_own_fission_deficit(body) -> None:
    """[E3; l1, ``characteristic-fixed-source``] On a subcritical body the source (1/k - 1) F phi_k returns phi_k: the k pencil and the
    source pencil read one system. To 1e-12 of max phi_k.

    With the fission in the source pencil's gain, the emission of that source
    is (S + F/k) phi_k, the mode's own; it needs k < 1, asserted. `[M]`
    3.4e-14 (sphere, k 0.075), 3.0e-15 (slab, k 0.187). First red: the source
    pencil's scattering assembled with another pairing than the k pencil's
    (a twin assembly).
    """
    chart, breakpoints, laws = body
    system = _system(chart, breakpoints, laws, (_UP2N, _ABS, _UP2N))
    mode = system.pencil.fundamental()
    assert mode.k < 1.0
    flux = system.flux(mode.vector)
    fission = np.einsum("igh,hi->gi", system.cross_sections.fission[system.basis.region], flux)
    returned = system.fixed_source((1.0 / mode.k - 1.0) * fission)
    assert np.max(np.abs(returned - flux)) / flux.max() < 1e-12


@pytest.mark.foundation
@pytest.mark.parametrize("body", _E2_BODIES, ids=["sphere3-mirror", "slab3-mirrors", "hollow-white-mirror"])
@pytest.mark.rests_on(_HERE + "test_a_closed_body_with_a_uniform_source_reads_the_infinite_medium_flux")
def test_a_source_in_a_body_supercritical_by_its_n2n_emission_is_refused(body) -> None:
    """[E6] A fixed source in a closed body made supercritical by (n,2n) alone (no fission) is refused through
    ``NoLeastSolution``; control: the k pencil's direct solve of the same load returns a finite flux with a negative
    entry, so the refusal is not vacuous. `[M]` 2026-10-07 (``ta/m11_e6.log``): on each of the three bodies the
    direct solve's emission coefficients reach -30 and the flux read from them -20; the -20 of the spec was the flux
    form's unknown, the -30 the archivist measured is this row's (the emission).

    First red: the source question answered on the k pencil, whose gain is the
    fission alone.
    """
    chart, breakpoints, laws = body
    system = _system(chart, breakpoints, laws, (_N2N,) * (len(breakpoints) - 1))
    source = np.ones((2, system.basis.size))
    with pytest.raises(NoLeastSolution, match="spectral radius of the gain"):
        system.fixed_source(source)
    direct = np.linalg.solve(system.pencil.loss, system.emission_mass @ system.emission.restrict(source))
    assert np.all(np.isfinite(direct)) and direct.min() < 0.0


# ── 4.4 the adjoint ──────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.verifies("characteristic-adjoint")
@pytest.mark.parametrize(("chart", "breakpoints", "laws", "name"),
                         _closed_params(_CLOSED[:4], _CLOSED_IDS[:4], _CLOSED_MIXTURES))
@pytest.mark.rests_on(_HERE + "test_a_closed_body_reads_k_inf_with_a_flat_flux_in_the_infinite_mediums_group_ratio")
def test_the_closed_body_adjoint_is_the_infinite_medium_importance(chart, breakpoints, laws, name) -> None:
    """[D13 (i), (ii); l1, ``characteristic-adjoint``] The k-adjoint (``pencil.adjoint()``) of a closed body:
    its fundamental eigenvalue is k; its vector, the adjoint flux on the emission space, is flat in each group with
    the group ratio of A^-T nu Sigma_f, by hand, to 1e-12.

    First red: the flux form's transpose, whose vector is Sigma_t (.)
    phi^dagger, the adjoint collision density (`[M]` 2026-10-07: off the
    ratio by 0.10 to 0.38 on these mixtures, `ta/m1_closed.log`).
    """
    mixture = _CLOSED_MIXTURES[name]
    system = _system(chart, breakpoints, laws, (mixture,) * (len(breakpoints) - 1))
    loss, production = _zero_d(mixture)
    _, importance = _dominant(np.linalg.solve(loss.T, production.T))
    forward, adjoint = system.pencil.fundamental(), system.pencil.adjoint().fundamental()
    assert abs(adjoint.k - forward.k) / forward.k < 1e-12
    groups = system.emission.split(adjoint.vector)
    means = np.array([g.mean() for g in groups])
    assert max(np.max(np.abs(g / g.mean() - 1.0)) for g in groups) < 1e-10
    np.testing.assert_allclose(means / means.sum(), importance, rtol=0.0, atol=1e-12)


@pytest.mark.l1
@pytest.mark.verifies("characteristic-adjoint")
@pytest.mark.parametrize("body", _RR_BODIES, ids=["sphere", "slab", "hollow-sphere"])
@pytest.mark.rests_on(_HERE + "test_one_group_k_increases_on_nested_spaces")
def test_the_one_group_adjoint_is_the_forward_flux(body) -> None:
    """[D13 (iii); l1, ``characteristic-adjoint``] In one group the transport is self-adjoint: the k-adjoint's vector on the emission
    space is the forward flux on the same nodes, heterogeneous body, to 1e-10 (both unit-normalised).

    On the support W_s phi_s = K_s q and q = M phi_s, so the forward flux
    solves W_s phi = K_s M phi, the adjoint's own equation (K_s symmetric).
    First red: the flux form's transpose, whose vector is the emission M phi.
    """
    chart, breakpoints, laws = body
    system = _system(chart, breakpoints, laws, _ONE_GROUP)
    forward = system.flux(system.pencil.fundamental().vector)[0, system.emission.supports[0]]
    adjoint = system.pencil.adjoint().fundamental().vector
    np.testing.assert_allclose(adjoint, forward / np.linalg.norm(forward), rtol=0.0, atol=1e-10)


#: qa's adjoint fixture: upscatter and downscatter, an absorber that emits into group 1 only (the supports differ by
#: group), (n,2n) and fission with chi in both groups.
_ADJ_FUEL = _mixture([1.0, 1.5], [[0.5, 0.3], [0.08, 1.0]], [[0.03, 0.02], [0.0, 0.0]], [0.05 * 2.5, 0.1 * 2.4],
                     [0.7, 0.3])
_ADJ_ABSORBER = _mixture([0.8, 2.0], [[0.0, 0.3], [0.0, 0.0]])
_ADJ_REFLECTOR = _mixture([0.9, 1.4], [[0.6, 0.28], [0.0, 1.35]])
_ADJ_BODY = (_ADJ_FUEL, _ADJ_ABSORBER, _ADJ_REFLECTOR)
_ADJ_BREAKPOINTS = (0.0, 0.7, 1.1, 1.8)


def _group_transposed(mixture: Mixture) -> Mixture:
    """The mixture whose transfers run the other way: SigS, Sig2 transposed; fission chi (x) nu Sigma_f -> its
    transpose, spelled as nu Sigma_f' = chi * sum(nu Sigma_f), chi' = nu Sigma_f / sum(nu Sigma_f) (rank one)."""
    total_production = float(np.sum(mixture.SigP))
    chi = np.asarray(mixture.SigP) / total_production if total_production else np.zeros_like(mixture.SigP)
    nu_sig_f = np.asarray(mixture.chi) * total_production
    return Mixture(SigC=mixture.SigC, SigL=mixture.SigL, SigF=nu_sig_f / 2.5, SigP=nu_sig_f, SigT=mixture.SigT,
                   SigS=[csr_matrix(mixture.SigS[0].toarray().T)], Sig2=[csr_matrix(mixture.Sig2[0].toarray().T)],
                   chi=chi)


_ADJ_CASES = [
    ("sphere", (_VACUUM,)), ("sphere", ((0.0, 0.6),)),
    ("slab", (_MIRROR, _VACUUM)), ("slab", (_WHITE, _VACUUM)), ("slab", ("P", "P")),
]


@pytest.mark.l1
@pytest.mark.verifies("characteristic-adjoint")
@pytest.mark.parametrize(("chart", "laws"), _ADJ_CASES,
                         ids=["sphere-vacuum", "sphere-white0.6", "slab-mirror-vacuum", "slab-white-vacuum", "slab-periodic"])
@pytest.mark.rests_on(_HERE + "test_the_closed_body_adjoint_is_the_infinite_medium_importance")
def test_the_adjoint_is_the_forward_solve_of_the_group_transposed_problem(chart: str, laws) -> None:
    """[D13 (iv); l1, ``characteristic-adjoint``] An independent route to the adjoint: the response to a
    detector r, and the k-adjoint's fundamental, equal the FORWARD solves (``fixed_source`` of r, the fundamental flux)
    of the problem with every transfer transposed, posed with sources everywhere, read on the original supports;
    to 1e-12 (response) and 1e-10 (mode, unit-normalised).

    The transposed problem is written from the mixtures here (``_group_transposed``),
    and its system assembles its OWN blocks over the full support, so the
    agreement rests on the reciprocity of the assembled transport (K_g's
    support columns are the transposes of its support rows), which no
    algebraic identity supplies. The supports differ by group (the absorber
    emits into group 1 only). The periodic slab is closed and supercritical
    under the source question, so it carries the mode leg alone. `[M]` qa
    2026-10-07: 3.4e-15 to 8.3e-15 (``qa/p1b_adjoint_disjoint_support.py``,
    ``p4_slab_adjoint.py``). First red: a non-reciprocal perturbation of K.
    """
    transposed = tuple(_group_transposed(m) for m in _ADJ_BODY)
    system = _system(chart, _ADJ_BREAKPOINTS, laws, _ADJ_BODY)
    mirror_problem = _system(chart, _ADJ_BREAKPOINTS, laws, transposed, source_regions=((True, True),) * 3)
    supports = system.emission.supports
    assert supports[0].size != supports[1].size                          # the supports differ by group
    np.testing.assert_array_equal(system.cross_sections.scattering, np.swapaxes(mirror_problem.cross_sections.scattering, 1, 2))
    np.testing.assert_allclose(system.cross_sections.fission, np.swapaxes(mirror_problem.cross_sections.fission, 1, 2),
                               rtol=1e-15, atol=0.0)

    def on_supports(field: np.ndarray) -> np.ndarray:
        return np.concatenate([field[g, s] for g, s in enumerate(supports)])

    if laws != ("P", "P"):
        detector = np.random.default_rng(5).random((2, system.basis.size))
        response = system.response(detector)
        forward = on_supports(mirror_problem.fixed_source(detector))
        assert np.max(np.abs(response - forward)) / np.max(np.abs(forward)) < 1e-12
    adjoint = system.pencil.adjoint().fundamental()
    mode = mirror_problem.pencil.fundamental()
    assert abs(adjoint.k - mode.k) / mode.k < 1e-12
    reference = on_supports(mirror_problem.flux(mode.vector))
    np.testing.assert_allclose(adjoint.vector, reference / np.linalg.norm(reference), rtol=0.0, atol=1e-10)


_HETEROGENEOUS = ("sphere", _MR3, (_VACUUM,), (_UP2N, _ABS, _PU2))


@pytest.mark.l1
@pytest.mark.verifies("characteristic-adjoint")
@pytest.mark.rests_on(_HERE + "test_the_closed_body_adjoint_is_the_infinite_medium_importance")
def test_the_forward_and_adjoint_modes_are_biorthogonal() -> None:
    """[D14 (i); l1, ``characteristic-adjoint``] The first four forward modes x_j and adjoint modes y_i of a
    heterogeneous 2-group vacuum sphere: the adjoint eigenvalues are the forward ones, and
    |y_i^T B x_j| / (|y_i| |B x_j|) is below 1e-10 for i != j and above 1e-3 for i = j (B the production form).

    First red: the adjoint pencil left untransposed (the forward modes used as
    their own adjoints, which a non-normal pencil does not admit).
    """
    pencil = _system(*_HETEROGENEOUS).pencil
    forward, adjoint = pencil.spectrum()[:4], pencil.adjoint().spectrum()
    production = pencil.production
    for i, f in enumerate(forward):
        match = min(adjoint, key=lambda a: abs(a.eigenvalue - f.eigenvalue))
        assert abs(match.eigenvalue - f.eigenvalue) < 1e-10 * abs(f.eigenvalue)
        for j, g in enumerate(forward):
            bx = production @ g.vector
            pairing = abs(match.vector @ bx) / (np.linalg.norm(match.vector) * np.linalg.norm(bx))
            assert (pairing > 1e-3) if i == j else (pairing < 1e-10), (i, j, pairing)


@pytest.mark.l1
@pytest.mark.verifies("characteristic-fixed-source")
@pytest.mark.rests_on(_HERE + "test_the_source_question_returns_the_mode_from_its_own_fission_deficit")
def test_a_detectors_reading_of_a_source_is_its_response_paired_with_the_source() -> None:
    """[D14 (ii); l1, ``characteristic-fixed-source``] ``response`` transposes ``fixed_source``:
    <r, phi(q)>_W = phi^dagger(r)^T W_s q on the emission space, 20 seeded non-negative (q, r) pairs on a
    subcritical heterogeneous body, to 1e-12 relative.

    It checks only that: the identity holds algebraically for ANY transport
    matrix K (qa, 2026-10-07: K perturbed by a non-reciprocal factor
    1 + 0.3 U read 0.0). The physics of the adjoint (that the transpose of
    the assembled K is the adjoint transport) is D13 (iv)'s, against the
    forward solve of the group-transposed problem. First red: the
    response's load spelled K^T W_G r (the metric applied twice).
    """
    system = _system("sphere", _MR3, (_VACUUM,), (_UP2N, _ABS, _UP2N))
    rng = np.random.default_rng(20261007)
    shape = (2, system.basis.size)
    for _ in range(20):
        q, r = rng.random(shape), rng.random(shape)
        reading = r.ravel() @ system.mass @ system.fixed_source(q).ravel()
        paired = system.response(r) @ system.emission_mass @ system.emission.restrict(q)
        assert abs(reading - paired) < 1e-12 * abs(reading)
