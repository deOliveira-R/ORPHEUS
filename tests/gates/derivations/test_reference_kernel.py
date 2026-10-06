r"""The reference kernel's gates (P0.5): the dense pencil, the least solution, the composite Gauss rule.

Every matrix here is DOMAIN-FREE: a pencil built from chosen eigenpairs, a gain built
from a chosen spectrum (``instrument-doctrine`` X4: a primitive two branches rely on is
pinned against a closed form, never against an in-domain solver). The weak-form pencil
:math:`F c = k L c` is manufactured as :math:`F = L V \Lambda V^{-1}` with a non-symmetric,
well-conditioned :math:`L`, so :math:`F v_i = \lambda_i L v_i` holds by construction.

The rows read the kernel through its MODULE (``dp.DensePencil`` and its ``least_solution``),
so a mutation battery that rebinds a module attribute reaches every row.

The theory and the claims these rows verify: :ref:`verification-reference-kernel`
(``docs/theory/verification/reference_solutions.rst``), with each row's first red in its gate
table. K4, the layer allowlist, lives in ``tests/gates/test_layer_imports.py``.
"""
from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest

import orpheus.derivations.common.dense_pencil as dp
import orpheus.derivations.common.quadrature as quad

pytestmark = pytest.mark.foundation

_SEED = 20261006


# ── manufactured pencils ────────────────────────────────────────────────────


def _loss(n: int, rng: np.random.Generator) -> np.ndarray:
    """A non-symmetric, well-conditioned loss form (cond about 2-4)."""
    return np.eye(n) + 0.3 * rng.standard_normal((n, n)) / np.sqrt(n)


def _basis(n: int, rng: np.random.Generator, *, first_column: np.ndarray | None = None) -> np.ndarray:
    v = np.eye(n) + 0.25 * rng.standard_normal((n, n)) / np.sqrt(n)
    v[:, 0] = 1.0 + rng.random(n) if first_column is None else first_column
    return v


def _pencil(blocks: np.ndarray, v: np.ndarray, loss: np.ndarray) -> "dp.DensePencil":
    """F = L V Lambda V^{-1}; ``blocks`` is the real (block-diagonal) Lambda."""
    production = loss @ v @ blocks @ np.linalg.inv(v)
    return dp.DensePencil(loss, production)


def _real_blocks(*entries) -> np.ndarray:
    """Block-diagonal Lambda from reals and (re, im) pairs, the pair as [[re, im], [-im, re]]."""
    size = sum(1 if np.isscalar(e) else 2 for e in entries)
    lam = np.zeros((size, size))
    i = 0
    for e in entries:
        if np.isscalar(e):
            lam[i, i] = e
            i += 1
        else:
            re, im = e
            lam[i:i + 2, i:i + 2] = [[re, im], [-im, re]]
            i += 2
    return lam


def _parallel(a: np.ndarray, b: np.ndarray) -> float:
    """1 - |cos angle| between two real vectors: 0 when parallel."""
    return float(1.0 - abs(float(a @ b)) / (np.linalg.norm(a) * np.linalg.norm(b)))


# ── K1: the fundamental mode ────────────────────────────────────────────────


@pytest.mark.parametrize("n", [5, 40])
@pytest.mark.verifies("reference-kernel-fundamental-contract")
def test_k1_the_fundamental_is_the_chosen_dominant_eigenpair(n: int) -> None:
    """K1. The dominant eigenvalue 2.5 and its positive eigenvector V[:, 0] are returned.

    The rest of the spectrum holds a negative eigenvalue (-0.7), a complex pair (0.3 +- 0.2i)
    and, at n = 40, a geometric tail; none dominates. The band on k is 1e-13 relative (a
    QZ reduction of a cond-4 pencil; LAPACK's bits are never pinned, ``vv`` #38).
    """
    rng = np.random.default_rng(_SEED + n)
    tail = tuple(0.2 * 0.9 ** j for j in range(n - 5))
    lam = _real_blocks(2.5, 1.1, -0.7, (0.3, 0.2), *tail)
    v = _basis(n, rng)
    mode = _pencil(lam, v, _loss(n, rng)).fundamental()
    assert abs(mode.k - 2.5) <= 1e-13 * 2.5, f"k = {mode.k!r}, expected 2.5"
    assert np.all(mode.vector >= 0.0), "the fundamental's vector must be non-negative"
    assert _parallel(mode.vector, v[:, 0]) < 1e-13, "the fundamental's vector is not V[:, 0]"


def test_k1_a_within_band_negative_coefficient_is_clipped_to_exactly_zero() -> None:
    """K1, the clip: a coefficient that is 0 in exact arithmetic comes out of QZ at -1.2e-16 (seed 5, measured);
    within the pencil's sign band it is admitted and set to exactly 0.0, every other coefficient untouched."""
    rng = np.random.default_rng(5)
    first = np.array([1.0, 0.7, 0.0, 0.4, 0.9])          # an exact zero at index 2
    pencil = _pencil(_real_blocks(2.5, 1.1, -0.7, 0.4, 0.15), _basis(5, rng, first_column=first), _loss(5, rng))
    raw = np.real(pencil.spectrum()[0].vector)
    raw = raw if raw.sum() >= 0.0 else -raw
    assert -1e-14 < raw[2] < 0.0, f"premise: the exact zero comes back slightly negative, got {raw[2]!r}"
    mode = pencil.fundamental()
    assert mode.vector[2] == 0.0, f"the rounding coefficient is not clipped: {mode.vector[2]!r}"
    keep = np.arange(5) != 2
    assert np.array_equal(mode.vector[keep], raw[keep]), "the clip touched a coefficient outside the band"


# ── K1r: the refusals, one input per guard, each keyed to its own fragment ─────


_FRAGMENTS = {
    "complex": "is complex",
    "negative": "is not positive",
    "degenerate": "is not strictly dominant",
    "mixed_sign": "is not single-signed",
}


def test_k1r_the_refusal_fragments_are_disjoint() -> None:
    """Each refusal's fragment names only its own condition (so a leg cannot pass on a neighbour's message)."""
    values = list(_FRAGMENTS.values())
    for a in values:
        assert sum(a in b for b in values) == 1, f"fragment {a!r} is contained in another"


def _refused(lam: np.ndarray, first: np.ndarray | None = None) -> "dp.DensePencil":
    rng = np.random.default_rng(_SEED + 7)
    n = lam.shape[0]
    return _pencil(lam, _basis(n, rng, first_column=first), _loss(n, rng))


@pytest.mark.verifies("reference-kernel-fundamental-contract")
def test_k1r_a_negative_dominant_eigenvalue_is_refused() -> None:
    """(-3, 2, 1): the largest MODULUS is -3, so there is no fundamental; picking the largest real part (2) would answer."""
    pencil = _refused(_real_blocks(-3.0, 2.0, 1.0, 0.5))
    assert abs(pencil.spectrum()[0].eigenvalue - (-3.0)) < 1e-12, "premise: -3 leads by modulus"
    with pytest.raises(dp.NoFundamentalMode, match=_FRAGMENTS["negative"]):
        pencil.fundamental()


@pytest.mark.verifies("reference-kernel-fundamental-contract")
def test_k1r_a_degenerate_dominant_eigenvalue_is_refused() -> None:
    """(2, 2, 1, 0.5): a two-dimensional dominant eigenspace; and (2, -2, 1, 0.5): equal modulus, opposite sign."""
    for lam in (_real_blocks(2.0, 2.0, 1.0, 0.5), _real_blocks(2.0, -2.0, 1.0, 0.5)):
        with pytest.raises(dp.NoFundamentalMode, match=_FRAGMENTS["degenerate"]):
            _refused(lam).fundamental()


@pytest.mark.verifies("reference-kernel-fundamental-contract")
def test_k1r_a_mixed_sign_dominant_eigenvector_is_refused() -> None:
    """V[:, 0] with one coefficient at -1e-3 of the largest: far outside the pencil's sign band (``DensePencil.sign_band``)."""
    first = np.array([1.0, 0.8, -1e-3, 0.6, 0.9])
    with pytest.raises(dp.NoFundamentalMode, match=_FRAGMENTS["mixed_sign"]):
        _refused(_real_blocks(2.5, 1.1, -0.7, 0.4, 0.15), first).fundamental()


@pytest.mark.verifies("reference-kernel-fundamental-contract")
def test_k1r_a_complex_dominant_pair_is_refused() -> None:
    """2 +- 1e-6 i dominant. For a REAL pencil a complex dominant eigenvalue always has its conjugate at equal
    modulus, so the dominance guard would refuse it too; the complex guard runs first and owns this message."""
    with pytest.raises(dp.NoFundamentalMode, match=_FRAGMENTS["complex"]):
        _refused(_real_blocks((2.0, 1e-6), 1.0, 0.5)).fundamental()


def test_k1r_a_singular_loss_is_refused_by_the_spectrum() -> None:
    """A rank-deficient loss form: the QZ reduction returns an infinite eigenvalue."""
    loss = np.diag([1.0, 1.0, 0.0])
    with pytest.raises(ValueError, match="the loss form is singular"):
        dp.DensePencil(loss, np.eye(3)).spectrum()


# ── K5: the full spectrum ──────────────────────────────────────────────────


def test_k5_the_spectrum_holds_every_mode_by_decreasing_modulus_with_no_refusal() -> None:
    """Every chosen eigenvalue is returned, complex and negative ones included, ordered by |k|.

    The order is the discriminator: by real part the pair 0.3 +- 0.6i (modulus 0.67) would sort
    below -0.7 and 0.5; by modulus it sorts between them.
    """
    rng = np.random.default_rng(_SEED + 5)
    lam = _real_blocks(2.5, -0.7, (0.3, 0.6), 0.5)
    v = _basis(5, rng)
    modes = _pencil(lam, v, _loss(5, rng)).spectrum()
    got = np.array([m.eigenvalue for m in modes])
    expected = np.array([2.5, -0.7, 0.3 + 0.6j, 0.3 - 0.6j, 0.5])
    moduli = np.abs(got)
    assert np.all(np.diff(moduli) <= 1e-12), f"not ordered by decreasing modulus: {got}"
    assert np.allclose(sorted(got, key=lambda z: (round(abs(z), 9), z.imag)),
                       sorted(expected, key=lambda z: (round(abs(z), 9), z.imag)), rtol=0, atol=1e-12), got
    for m in modes:
        assert abs(np.linalg.norm(m.vector) - 1.0) < 1e-14
    real_modes = [m for m in modes if m.eigenvalue.imag == 0.0]
    for target, column in ((2.5, 0), (-0.7, 1), (0.5, 4)):
        m = min(real_modes, key=lambda m: abs(m.eigenvalue - target))
        assert _parallel(np.real(m.vector), v[:, column]) < 1e-12, f"mode {target}: wrong eigenvector"


# ── K6: the adjoint, in weak form ──────────────────────────────────────────


def _nonsymmetric_pencil() -> "dp.DensePencil":
    rng = np.random.default_rng(_SEED + 6)
    return _pencil(_real_blocks(2.5, 1.1, -0.7, 0.4, 0.15), _basis(5, rng), _loss(5, rng))


@pytest.mark.verifies("reference-kernel-adjoint-pencil")
def test_k6i_the_adjoint_has_the_same_spectrum() -> None:
    pencil = _nonsymmetric_pencil()
    forward = np.array([m.eigenvalue for m in pencil.spectrum()])
    backward = np.array([m.eigenvalue for m in pencil.adjoint().spectrum()])
    assert np.allclose(forward, backward, rtol=1e-12, atol=0.0), (forward, backward)


@pytest.mark.verifies("reference-kernel-biorthogonality")
def test_k6ii_adjoint_and_forward_modes_are_biorthogonal_in_the_loss_form() -> None:
    """u_i^T L v_j = 0 for i != j (relative to the diagonal), the modes paired by their (distinct) eigenvalues.

    The fixture's L is non-symmetric and V is not L-orthogonal, so an adjoint that returned the
    untransposed pencil (u_i = v_i) leaves off-diagonal entries of order 1e-1.
    """
    pencil = _nonsymmetric_pencil()
    v = np.array([np.real(m.vector) for m in pencil.spectrum()]).T
    u = np.array([np.real(m.vector) for m in pencil.adjoint().spectrum()]).T
    gram = u.T @ pencil.loss @ v
    diagonal = np.abs(np.diag(gram))
    assert np.all(diagonal > 1e-3), f"a mode is L-orthogonal to its own adjoint: {diagonal}"
    off = np.abs(gram - np.diag(np.diag(gram))) / np.sqrt(np.outer(diagonal, diagonal))
    assert off.max() < 1e-12, f"biorthogonality defect {off.max():.2e}"


@pytest.mark.verifies("reference-kernel-reciprocity")
def test_k6iii_the_source_solve_is_reciprocal_through_the_adjoint() -> None:
    """<d, x> = <x_adj, s>: x solves the pencil's forms, x_adj the adjoint's (mass = loss, gain = production)."""
    rng = np.random.default_rng(_SEED + 8)
    n = 6
    loss = _loss(n, rng)
    gain = 0.1 * rng.random((n, n))
    pencil = dp.DensePencil(loss, gain)
    adjoint = pencil.adjoint()
    s, d = rng.random(n), rng.random(n)
    x = pencil.least_solution(s)
    x_adj = adjoint.least_solution(d)
    lhs, rhs = float(d @ x), float(x_adj @ s)
    assert abs(lhs - rhs) <= 1e-13 * abs(lhs), f"<d, x> = {lhs!r}, <x_adj, s> = {rhs!r}"


# ── K2: the least solution ─────────────────────────────────────────────────


def _gain_with_spectrum(mu: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    v = _basis(len(mu), rng)
    return v @ np.diag(mu) @ np.linalg.inv(v)


@pytest.mark.parametrize("rho", [0.5, 0.9, 0.999])
@pytest.mark.verifies("reference-kernel-least-solution")
def test_k2_the_least_solution_with_unit_mass(rho: float) -> None:
    """mass = I: x = V diag(1/(1 - mu)) V^{-1} b, the Neumann sum in closed form."""
    rng = np.random.default_rng(_SEED + int(rho * 1000))
    n = 6
    mu = rho * np.array([1.0, 0.6, -0.4, 0.3, 0.1, -0.05])
    v = _basis(n, rng)
    gain = v @ np.diag(mu) @ np.linalg.inv(v)
    b = rng.random(n)
    expected = v @ np.diag(1.0 / (1.0 - mu)) @ np.linalg.solve(v, b)
    x = dp.DensePencil(np.eye(n), gain).least_solution(b)
    cond = np.linalg.cond(np.eye(n) - gain)
    assert np.max(np.abs(x - expected)) <= 50.0 * cond * np.finfo(float).eps * np.max(np.abs(expected)), (x, expected)


@pytest.mark.verifies("reference-kernel-least-solution")
def test_k2_the_least_solution_with_a_non_identity_mass_is_the_neumann_sum() -> None:
    """M x = G x + s with M != I, against the explicit Neumann sum sum_n (M^{-1}G)^n M^{-1}s, truncated at 1e-18."""
    rng = np.random.default_rng(_SEED + 11)
    n = 6
    mass = np.diag(1.0 + rng.random(n)) + 0.1 * rng.random((n, n))
    gain = 0.15 * rng.random((n, n))
    s = rng.random(n)
    step = np.linalg.solve(mass, gain)
    term = np.linalg.solve(mass, s)
    total = np.zeros(n)
    while np.max(np.abs(term)) > 1e-18 * max(1.0, np.max(np.abs(total))):
        total += term
        term = step @ term
    x = dp.DensePencil(mass, gain).least_solution(s)
    assert np.max(np.abs(x - total)) <= 1e-13 * np.max(np.abs(total)), (x, total)
    assert np.all(x >= 0.0), "a positive gain and source give a non-negative least solution"


# ── K2r: the refusals and the zero source ──────────────────────────────────


def _positive_gain(rho: float, n: int, rng: np.random.Generator) -> np.ndarray:
    a = 0.1 + rng.random((n, n))
    return rho * a / np.max(np.abs(np.linalg.eigvals(a)))


def test_k2r_a_supercritical_gain_is_refused_where_a_direct_solve_answers_negatively() -> None:
    """rho = 1.5, I - G non-singular. The control leg shows the refusal is not vacuous: the bare solve
    returns a finite flux with a negative entry."""
    rng = np.random.default_rng(_SEED + 12)
    gain = _positive_gain(1.5, 5, rng)
    s = rng.random(5)
    bare = np.linalg.solve(np.eye(5) - gain, s)
    assert np.all(np.isfinite(bare)) and bare.min() < 0.0, "premise: the unguarded solve answers, negatively"
    with pytest.raises(dp.NoLeastSolution, match="the Neumann series diverges"):
        dp.DensePencil(np.eye(5), gain).least_solution(s)


def test_k2r_the_radius_is_measured_against_the_mass() -> None:
    """rho(G) = 0.6 but rho(M^{-1} G) = 1.2 (M = I/2): refused; and rho(G) = 1.2 with M = 2I, rho(M^{-1}G) = 0.6: answered."""
    rng = np.random.default_rng(_SEED + 13)
    s = rng.random(5)
    g = _positive_gain(0.6, 5, rng)
    with pytest.raises(dp.NoLeastSolution, match="the Neumann series diverges"):
        dp.DensePencil(0.5 * np.eye(5), g).least_solution(s)
    x = dp.DensePencil(2.0 * np.eye(5), 2.0 * g).least_solution(s)
    assert np.allclose(2.0 * x - 2.0 * g @ x, s, rtol=1e-13, atol=0.0)


def test_k2r_the_margin_decides_every_draw_on_both_sides() -> None:
    """A radius inside the margin (1 - 1e-11, 1 + 1e-12) is refused; one 10 margins below 1 is answered.

    The computed radius of a gain whose radius is exactly 1 scatters (QZ, ``[M]`` 2026-10-06, the
    archivist: 20 seeds of 200 draws of this construction span 1 - 7.2e-13 to 1 + 1.8e-13), so a bare
    comparison with 1 is not decidable there; the kernel refuses every radius not below
    ``1 - SUBCRITICAL_MARGIN`` (1e-10), about 140 times beyond the worst draw. First red: with the
    margin set to 0, the radius-exactly-1 row admits draws (23 of this row's 200).
    """
    rng = np.random.default_rng(_SEED + 14)
    margin = dp.SUBCRITICAL_MARGIN
    for _ in range(200):
        s = rng.random(4)
        for inside in (1.0 + 1e-12, 1.0 - 0.1 * margin):
            gain = _gain_with_spectrum(np.array([inside, 0.5, 0.2, -0.3]), rng)
            with pytest.raises(dp.NoLeastSolution, match="the Neumann series diverges"):
                dp.DensePencil(np.eye(4), gain).least_solution(s)
        below = _gain_with_spectrum(np.array([1.0 - 10.0 * margin, 0.5, 0.2, -0.3]), rng)
        assert np.all(np.isfinite(dp.DensePencil(np.eye(4), below).least_solution(s)))


def test_k2r_a_unit_spectral_radius_is_refused_in_every_draw() -> None:
    """mu_max = 1 exactly in the construction (I - G singular): refused in every one of 200 draws."""
    rng = np.random.default_rng(_SEED + 14)
    admitted = 0
    for _ in range(200):
        gain = _gain_with_spectrum(np.array([1.0, 0.5, 0.2, -0.3]), rng)
        try:
            dp.DensePencil(np.eye(4), gain).least_solution(rng.random(4))
        except dp.NoLeastSolution:
            continue
        except np.linalg.LinAlgError:           # admitted, then the solve of the singular I - G fails untyped
            pass
        admitted += 1
    assert admitted == 0, f"{admitted} of 200 unit-radius gains admitted"


def test_k2r_a_zero_source_returns_exact_zero_at_any_radius() -> None:
    """The least solution of a source-free problem is 0 even when rho >= 1 (a lossless trapped line)."""
    rng = np.random.default_rng(_SEED + 15)
    for rho in (0.5, 1.0, 1.5):
        x = dp.DensePencil(np.eye(4), _positive_gain(rho, 4, rng)).least_solution(np.zeros(4))
        assert np.array_equal(x, np.zeros(4)), f"rho = {rho}: {x}"


# ── K3: the composite Gauss rule ───────────────────────────────────────────


_BREAKS = (Fraction(0), Fraction(3, 10), Fraction(11, 10), Fraction(2))
_N = 4                                                    # exact through degree 2n - 1 = 7 per panel
_COEFFS = (                                               # a different degree-7 polynomial per panel: a jump at each break
    (1, -2, 3, 0, 1, 0, -1, 2),
    (-3, 1, 0, 2, -1, 1, 0, -1),
    (2, 0, -1, 1, 0, -2, 1, 1),
)


def _exact_integral() -> Fraction:
    total = Fraction(0)
    for (a, b), coeffs in zip(zip(_BREAKS, _BREAKS[1:]), _COEFFS):
        total += sum(Fraction(c) * (b ** (k + 1) - a ** (k + 1)) / (k + 1) for k, c in enumerate(coeffs))
    return total


def _piecewise(x: np.ndarray) -> np.ndarray:
    edges = np.array([float(b) for b in _BREAKS])
    panel = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, len(_COEFFS) - 1)
    return np.array([np.polyval(_COEFFS[p][::-1], xi) for p, xi in zip(panel, x)])


def test_k3_the_composite_rule_integrates_a_jumping_piecewise_polynomial_exactly() -> None:
    """Degree 2n - 1 on each panel with jumps at the breakpoints; band 8 ulp of the sum of |panel integrals|."""
    rule = quad.composite_gauss_legendre([float(b) for b in _BREAKS], _N)
    assert len(rule.pts) == _N * (len(_BREAKS) - 1)
    assert np.all(np.diff(rule.pts) > 0.0), "the nodes are increasing"
    edges = [float(b) for b in _BREAKS]
    for a, b in zip(edges, edges[1:]):
        inside = (rule.pts > a) & (rule.pts < b)
        assert inside.sum() == _N, f"panel ({a}, {b}) holds {inside.sum()} nodes, not {_N}"
    got = float(rule.wts @ _piecewise(rule.pts))
    exact = float(_exact_integral())
    scale = sum(abs(float(sum(Fraction(c) * (b ** (k + 1) - a ** (k + 1)) / (k + 1) for k, c in enumerate(co))))
                for (a, b), co in zip(zip(_BREAKS, _BREAKS[1:]), _COEFFS))
    assert abs(got - exact) <= 8 * np.finfo(float).eps * scale, f"{got!r} vs exact {exact!r}"


def test_k3_the_weights_of_each_panel_sum_to_its_width() -> None:
    edges = [float(b) for b in _BREAKS]
    rule = quad.composite_gauss_legendre(edges, _N)
    for a, b in zip(edges, edges[1:]):
        inside = (rule.pts > a) & (rule.pts < b)
        assert abs(rule.wts[inside].sum() - (b - a)) <= 4 * np.finfo(float).eps * (b - a)


def test_a_fundamental_mode_holds_its_invariant_whoever_builds_it() -> None:
    """The type, not only ``fundamental()``, carries the Krein--Rutman contract; its vector is read-only.

    First red: a ``FundamentalMode`` without its ``__post_init__`` constructs ``(-1, [-1, 2])``.
    """
    for k, vector in ((-1.0, [0.6, 0.8]), (np.nan, [0.6, 0.8]), (1.0, [-0.6, 0.8]), (1.0, [0.6, 0.6])):
        with pytest.raises(dp.NoFundamentalMode, match="a fundamental"):
            dp.FundamentalMode(k, np.array(vector))
    mode = dp.FundamentalMode(1.0, np.array([0.6, 0.8]))
    with pytest.raises(ValueError, match="read-only"):
        mode.vector[0] = 0.0


@pytest.mark.parametrize(
    ("gain", "source", "expected"),
    [
        # Two decoupled classes: the supercritical one (radius 1) is not reached.
        ([[1.0, 0.0], [0.0, 0.5]], [0.0, 1.0], [0.0, 2.0]),
        # Class 0 feeds class 1 (G[1, 0] = 0.3), not back: a source on 1 never reaches 0.
        ([[1.0, 0.0], [0.3, 0.5]], [0.0, 1.0], [0.0, 2.0]),
    ],
)
@pytest.mark.verifies("reference-kernel-reach")
def test_k2_an_unreached_supercritical_class_carries_zero(gain, source, expected) -> None:
    """The least solution exists when the source reaches no class of radius >= 1 (the user's ruling, 2026-10-06).

    First red: a radius check over the whole gain refuses both rows.
    """
    x = dp.DensePencil(np.eye(2), np.array(gain)).least_solution(np.array(source))
    np.testing.assert_array_equal(x, np.array(expected))


@pytest.mark.verifies("reference-kernel-reach")
def test_k2r_a_reached_supercritical_class_is_refused() -> None:
    """A source on class 0 reaches the radius-1 class (itself), whether or not it also reaches class 1."""
    for gain in ([[1.0, 0.0], [0.0, 0.5]], [[1.0, 0.0], [0.3, 0.5]]):
        with pytest.raises(dp.NoLeastSolution, match="unknowns the source reaches"):
            dp.DensePencil(np.eye(2), np.array(gain)).least_solution(np.array([1.0, 0.0]))


@pytest.mark.verifies("reference-kernel-reach")
def test_k2_reach_is_the_downstream_closure_of_the_source() -> None:
    """Unknown j enters equation i when L[i, j] or F[i, j] is non-zero; a chain 0 -> 1 -> 2 and an island 3."""
    loss = np.eye(4)
    loss[1, 0] = -0.1
    production = np.zeros((4, 4))
    production[2, 1] = 0.2
    pencil = dp.DensePencil(loss, production)
    np.testing.assert_array_equal(pencil.reach(np.array([1.0, 0, 0, 0])), [True, True, True, False])
    np.testing.assert_array_equal(pencil.reach(np.array([0, 0, 1.0, 0])), [False, False, True, False])
    np.testing.assert_array_equal(pencil.reach(np.zeros(4)), [False] * 4)


def _conditioned_pencil(cond: float, fundamental: np.ndarray, rng: np.random.Generator) -> dp.DensePencil:
    """A pencil whose loss form has condition number ``cond`` and whose fundamental vector is chosen."""
    n = fundamental.size
    u, _ = np.linalg.qr(rng.standard_normal((n, n)))
    w, _ = np.linalg.qr(rng.standard_normal((n, n)))
    loss = u @ np.diag(np.logspace(0, -np.log10(cond), n)) @ w.T
    v = rng.random((n, n)) + 0.1
    v[:, 0] = fundamental
    spectrum = np.diag([2.0, 1.1, 0.9, 0.7, 0.5, 0.4, 0.3, 0.2])
    return dp.DensePencil(loss, loss @ v @ spectrum @ np.linalg.inv(v))


@pytest.mark.verifies("reference-kernel-fundamental-contract")
@pytest.mark.parametrize("cond", [1e2, 1e6, 1e8])
def test_k1_the_sign_band_scales_with_conditioning(cond: float) -> None:
    """Exact zeros are admitted at every conditioning; a sign change of 1e-3 is refused while it is resolvable.

    ``[M]`` qa 2026-10-06: a fixed band of 1e-10 refused 9 of 20 pencils at cond(L) near 1e8,
    whose exact zeros return at up to 5.6e-9. First red: ``sign_band`` returning 1e-10.
    """
    rng = np.random.default_rng(_SEED + 30)
    exact_zeros = np.array([0.0, 1.0, 0.5, 0.0, 0.7, 0.3, 0.0, 0.9])
    for _ in range(20):
        mode = _conditioned_pencil(cond, exact_zeros, rng).fundamental()  # admitted, not refused
        assert mode.vector.min() >= 0.0
    sign_change = exact_zeros.copy()
    sign_change[3] = -1e-3
    for _ in range(20):
        with pytest.raises(dp.NoFundamentalMode, match=_FRAGMENTS["mixed_sign"]):
            _conditioned_pencil(cond, sign_change, rng).fundamental()


def test_complex_input_is_refused_not_cast() -> None:
    """A complex form or source would lose its imaginary part to a cast; the kernel refuses it.

    First red: ``np.array(matrix, dtype=float)`` alone admits it with a ``ComplexWarning``.
    """
    with pytest.raises(ValueError, match="a real array is required"):
        dp.DensePencil(np.eye(2) * (1 + 1j), np.eye(2))
    with pytest.raises(ValueError, match="a real array is required"):
        dp.DensePencil(np.eye(2), 0.5 * np.eye(2)).least_solution(np.array([1.0, 1j]))
