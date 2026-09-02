r"""The Legendre basis on :math:`S^2/SO(2)_a` — the object that repairs ERR-080.

``solve_sn(scattering_order >= 2)`` on any 1-D chart returned a wrong answer
because a 1-D rule's frame bound the FULL real spherical harmonics to a
measure whose nodes had been forged onto :math:`S^2` as :math:`(\mu, 0, 0)`.
A point of :math:`[-1,1]` is an ORBIT of the :math:`SO(2)_x` action, not a
point of the sphere, and the harmonics read the fabricated azimuth as a real
one. The functions on that orbit space are the :math:`SO(2)`-invariant
functions on the base — the trivial isotypic component, one-dimensional in
every degree, spanned downstairs by :math:`P_\ell(\mu)`. This module is that
basis's own gate.

Claim layer and pillar (``vv-principles`` §"three pillars", ``AGENT.md`` §1.5)
------------------------------------------------------------------------------

Every row here is **term-level (L0)** against a **closed form** — the
Legendre polynomials' own defining laws: :math:`P_\ell(\pm 1)`, the
orthogonality :math:`\int_{-1}^{1} P_j P_k = \tfrac{2}{2j+1}\delta_{jk}`, the
three-term recurrence. None of it is an eigenvalue or a flux-shape claim, and
none of it needs the solver.

⚠ **The bit-identity rows are a REGRESSION contract, not a pillar** (``vv``
§bit-identity): they pin that the descended table is spelled exactly as the
harmonic table spells its :math:`m = 0` column, which is what lets the slab's
:math:`L \le 1` flux be ``array_equal`` across the repair rather than merely
close. They are paired with :func:`test_legendre_values_against_an_independent_three_term_recurrence`,
which carries the *correctness* claim from a source that shares no ``scipy``
call with either side — without it the bit rows would be a convention echo
(``vv-principles`` #22 / ``algebra-of-record`` §"structural independence
applies ABOVE the trusted-library line").

⭐ **The spelling is MEASURED, not chosen** (the verification memo's H-1).
No single ``scipy`` routine reproduces the harmonic :math:`m = 0` column
bit-for-bit: ``lpmv(0, 1, μ)`` misses at :math:`\ell = 1`, ``eval_legendre``
misses at :math:`\ell \ge 2`. :func:`test_the_descended_column_is_bit_identical_to_the_sh_m0_slot`
carries a **pure-**``lpmv`` **mutation arm** so the gate cannot pass by
echoing whatever the production spelling happens to be.
"""

from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal_nulp
from scipy.special import eval_legendre, lpmv

from orpheus.numerics.basis.legendre_basis import LegendreBasis, legendre_table
from orpheus.numerics.basis.spherical_harmonic_basis import SphericalHarmonicBasis
from orpheus.numerics.manifold import SPHERE
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces.legendre_space import LegendreSpace
from orpheus.numerics.symmetry import SubgroupOfO3

pytestmark = [pytest.mark.foundation]


#: Sphere rules whose nodes are honest unit directions — the (N, 3) coordinate
#: system of the orbit space. Built in the BODY of every test that uses them
#: (never in a ``parametrize`` argument list: a production call there runs at
#: collection, so a mutation battery reads ``rc=2 / 0 collected`` and a
#: ``^FAILED`` scanner scores it as "0 caught" — ``vv-principles`` #17).
_SPHERE_RULES: tuple[str, ...] = (
    "level_symmetric(4)",
    "level_symmetric(8)",
    "lebedev(11)",
    "lebedev(17)",
    "product(4,6)",
    "product(8,8)",
    "folded_product(4,8)",
)

_ONE_D_RULES: tuple[str, ...] = ("gauss_legendre(8)", "gauss_legendre(16)")


def _rule(label: str) -> Quadrature:
    """Build a shipped rule from its label — a LABEL is what parametrize carries."""
    match label:
        case "level_symmetric(4)":
            return Quadrature.level_symmetric(4)
        case "level_symmetric(8)":
            return Quadrature.level_symmetric(8)
        case "lebedev(11)":
            return Quadrature.lebedev(11)
        case "lebedev(17)":
            return Quadrature.lebedev(17)
        case "product(4,6)":
            return Quadrature.product(n_mu=4, n_phi=6)
        case "product(8,8)":
            return Quadrature.product(n_mu=8, n_phi=8)
        case "folded_product(4,8)":
            return Quadrature.folded_product(n_mu=4, n_phi=8)
    if label.startswith("gauss_legendre("):
        return Quadrature.gauss_legendre(int(label[15:-1]))
    raise ValueError(f"unknown rule label {label!r}")


# ══════════════════════════════════════════════════════════════════════
# §3 intrinsic laws — what a Legendre polynomial IS
# ══════════════════════════════════════════════════════════════════════


def test_endpoint_normalisation_is_exact() -> None:
    r""":math:`P_\ell(1) = 1` and :math:`P_\ell(-1) = (-1)^\ell`, at the BIT.

    The normalisation that fixes the family up to nothing — every other law
    here is invariant under :math:`P_\ell \mapsto c_\ell P_\ell`, so without
    this row a uniformly-rescaled table would satisfy the whole module.
    """
    basis = LegendreBasis(L=6)
    at_plus = basis.evaluate(np.array([1.0]))[0]
    at_minus = basis.evaluate(np.array([-1.0]))[0]

    assert np.array_equal(at_plus, np.ones(7)), f"P_l(1) = {at_plus}"
    assert np.array_equal(at_minus, (-1.0) ** np.arange(7)), f"P_l(-1) = {at_minus}"

    # NEGATIVE leg: a rescaled family satisfies orthogonality and the Gram's
    # structure and fails HERE. (vv #11 — a positive-only law is unfalsifiable.)
    rescaled = at_plus * (1.0 + 0.5 * np.arange(7))
    assert not np.array_equal(rescaled, np.ones(7))


@pytest.mark.parametrize("n_nodes,L", [(8, 2), (16, 4), (4, 3)])
def test_orthogonality_against_the_entrys_own_reference_measure(
    n_nodes: int, L: int
) -> None:
    r""":math:`\langle P_j, P_k\rangle = \tfrac{2}{2j+1}\delta_{jk}` on ``LEGENDRE``'s Gauss rule.

    The reference measure is the ENTRY's own
    (:attr:`Quotient.reference` is ``GeneratingMeasure(name='legendre',
    support=Interval(-1,1))``), so this is orthogonality against the mass-2
    normalisation the family is DEFINED by — a closed form, not another
    ORPHEUS object.
    """
    basis = LegendreBasis(L=L)
    gram = basis.mass_matrix(_rule(f"gauss_legendre({n_nodes})").measure)
    expected = np.diag(2.0 / (2.0 * np.arange(L + 1) + 1.0))

    np.testing.assert_allclose(
        gram,
        expected,
        atol=1e-14,
        err_msg=f"GL{n_nodes} L={L}: Gram is not diag(2/(2l+1))",
    )


def test_the_continuum_gram_is_the_pushforward_one_and_NOT_the_bare_mass_two() -> None:
    r"""The continuum Gram is :math:`4\pi/(2\ell+1)`, EQUAL to the harmonics' — not :math:`2/(2\ell+1)`.

    ⭐ This is the row that makes Landing A's binding metric-identical: the
    descent :math:`\mathrm{Funcs}(S^2)^{SO(2)} \to \mathrm{Funcs}(S^2/SO(2))`
    is an ISOMETRY only against the pushforward :math:`\pi_*d\Omega = 2\pi d\mu`,
    and the harmonics' own per-degree metric already IS that number. Declaring
    the bare Legendre mass-2 value instead would be a factor :math:`2\pi` off
    and would move every moment operator's end.

    Positive leg: ``array_equal`` to :attr:`SphericalHarmonicBasis.metric_per_ell`.
    Negative leg: NOT the bare value, and the ratio is exactly :math:`2\pi`.
    """
    L = 4
    legendre = LegendreBasis(L=L)
    harmonics = SphericalHarmonicBasis(L=L)

    assert np.array_equal(legendre.metric_per_ell, harmonics.metric_per_ell), (
        f"Legendre {legendre.metric_per_ell} vs SH {harmonics.metric_per_ell}"
    )

    bare = 2.0 / (2.0 * np.arange(L + 1) + 1.0)
    assert not np.allclose(legendre.metric_per_ell, bare), (
        "the continuum Gram must NOT be the bare mass-2 Legendre value"
    )
    np.testing.assert_allclose(legendre.metric_per_ell / bare, 2.0 * np.pi, rtol=0.0, atol=1e-14)


def test_the_dual_factor_is_two_ell_plus_one() -> None:
    r"""The canonical-dual (addition-theorem) factor is :math:`2\ell+1`, exactly as for the harmonics."""
    basis = LegendreBasis(L=5)
    assert np.array_equal(
        basis.addition_theorem_factor, 2.0 * np.arange(6) + 1.0
    ), f"{basis.addition_theorem_factor}"
    assert np.array_equal(
        basis.addition_theorem_factor, SphericalHarmonicBasis(L=5).addition_theorem_factor
    )


# ══════════════════════════════════════════════════════════════════════
# B1 — the discrete Gram, and the GL admissibility THEOREM
# ══════════════════════════════════════════════════════════════════════


def test_gram_is_diagonal_and_exact_against_the_legendre_gauss_rule() -> None:
    r"""B1 — ``GL_n``'s Legendre Gram is diagonal and exact for :math:`L \le n-1`.

    ``[M]`` 2026-09-02 at this tree: GL8/L=2 offdiag ``8.808e-17``, diag
    ``[2, 2/3, 0.4]`` to ``2.2e-16``; GL16/L=4 offdiag ``7.702e-16``, diag
    exact to ``4.4e-16``; GL4/L=3 offdiag ``1.994e-16``.
    """
    for n_nodes, L in ((8, 2), (16, 4), (4, 3)):
        gram = LegendreBasis(L=L).mass_matrix(_rule(f"gauss_legendre({n_nodes})").measure)
        diag = np.diag(gram)
        offdiag = float(np.max(np.abs(gram - np.diag(diag))))
        assert offdiag < 1e-14, f"GL{n_nodes} L={L}: offdiag {offdiag:.3e}"
        np.testing.assert_allclose(
            diag, 2.0 / (2.0 * np.arange(L + 1) + 1.0), atol=1e-14,
            err_msg=f"GL{n_nodes} L={L} diagonal",
        )
        assert LegendreBasis(L=L).gram_structure.value == "diagonal"


def test_a_gauss_rule_acquires_a_structurally_dead_slot_at_ell_equals_n() -> None:
    r"""B1's theorem row — the ADMISSIBILITY law of the whole 1-D family.

    :math:`P_n` vanishes identically at ``GL_n``'s nodes (they ARE its roots),
    so the Gram's :math:`(n,n)` entry is **structurally zero** — no tolerance
    is involved and no refinement helps. Consequences, all reachable from the
    public entry: **no 1-D Gauss frame can be dense AND full rank**, and
    ``gauss_legendre(2)`` at ``scattering_order >= 2`` is a rank-deficient
    frame that the DENSE arm's pseudo-inverse metric carries (user-ruled
    2026-09-02).

    ``[M]`` 2026-09-02: GL2/L=2 diag ``[2, 2/3, 0.0]`` rank 2 (offdiag
    ``5.6e-17``); GL8/L=8 diag ends ``0.0``, rank 8; GL2/L=4 offdiag
    ``7.778e-01`` and rank **2**, i.e. beyond the dead slot the Gram is no
    longer even diagonal.

    Both legs are here on purpose: the ALIVE row (``L = n-1``) and the DEAD
    row (``L = n``) differ only in the truncation order, so a table that
    silently regularised :math:`P_n` would pass the first and fail this.
    """
    for n_nodes in (2, 4, 8):
        alive = LegendreBasis(L=n_nodes - 1).mass_matrix(
            _rule(f"gauss_legendre({n_nodes})").measure
        )
        assert np.linalg.matrix_rank(alive) == n_nodes, (
            f"GL{n_nodes} at L={n_nodes - 1} must be full rank"
        )

        dead = LegendreBasis(L=n_nodes).mass_matrix(
            _rule(f"gauss_legendre({n_nodes})").measure
        )
        assert abs(dead[n_nodes, n_nodes]) < 1e-28, (
            f"GL{n_nodes}: G[{n_nodes},{n_nodes}] = {dead[n_nodes, n_nodes]:.3e} "
            f"must be the structurally dead slot"
        )
        assert np.linalg.matrix_rank(dead) == n_nodes


# ══════════════════════════════════════════════════════════════════════
# B2 — evaluate accepts BOTH honest coordinate systems of the orbit space
# ══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("label", _SPHERE_RULES[:4])
def test_evaluate_accepts_both_honest_coordinate_systems(label: str) -> None:
    r"""B2 — ``(N,)``/``(N,1)`` cosines and ``(N,3)`` directions agree, through :math:`\pi`.

    An orbit space has two charts, and a basis on it must eat points in
    either: the REALIZATION's coordinate (a 1-D rule's own nodes, values of
    :math:`\mu`) and the BASE's (unit directions of a full-sphere rule,
    pulled back along the entry's quotient map). ``[M]`` on the sphere rules
    ``π(nodes)`` is ``array_equal`` to ``axis_cosines(0)``, so the two charts
    are the same numbers and the two calls are ``array_equal``.
    """
    quadrature = _rule(label)
    basis = LegendreBasis(L=3)
    entry = SPHERE.quotient(SubgroupOfO3.SO2("x"))

    directions = np.asarray(quadrature.measure.nodes, dtype=float)
    cosines = entry.quotient_map(directions)

    assert np.array_equal(cosines, np.asarray(quadrature.axis_cosines(0), dtype=float)), (
        f"{label}: pi(nodes) is not the rule's own x-cosine"
    )

    from_base = basis.evaluate(directions)
    from_chart = basis.evaluate(cosines)
    from_column = basis.evaluate(cosines.reshape(-1, 1))

    assert from_base.shape == (directions.shape[0], 4)
    assert np.array_equal(from_base, from_chart), f"{label}: charts disagree"
    assert np.array_equal(from_chart, from_column), f"{label}: (N,) vs (N,1) disagree"


def test_a_wrong_width_is_refused_naming_BOTH_admissible_widths() -> None:
    r"""B2's negative leg — width 2 is neither chart, and the message says which two are.

    The refusal is a real ``raise``: under the canonical ``python -O`` a bare
    ``assert`` in production is stripped (``vv-principles`` Mode 8), so a
    contract spelled as one would not run in the suite that matters.
    ``pytest.raises`` here therefore pins a statement that genuinely executes.
    """
    basis = LegendreBasis(L=2)
    with pytest.raises(ValueError, match=r"\(N,\) / \(N, 1\) values of mu"):
        basis.evaluate(np.zeros((5, 2)))
    with pytest.raises(ValueError, match=r"\(N, 3\) unit directions"):
        basis.evaluate(np.zeros((5, 2)))

    # and the two charts it DOES name are admitted, so the message is not
    # merely well-worded (vv #11's positive leg on a refusal).
    basis.evaluate(np.zeros(5))
    basis.evaluate(np.zeros((5, 1)))
    basis.evaluate(np.column_stack([np.zeros(5), np.zeros(5), np.ones(5)]))


def test_points_off_the_orbit_space_are_refused_in_either_chart() -> None:
    r"""A cosine outside :math:`[-1,1]` and a non-unit direction are both refused."""
    basis = LegendreBasis(L=2)
    with pytest.raises(ValueError, match="not on the orbit space"):
        basis.evaluate(np.array([0.0, 1.5]))
    with pytest.raises(ValueError, match="not all lie on it"):
        basis.evaluate(np.array([[0.3, 0.4, 0.5]]))


def test_the_basis_states_its_orbit_space_and_its_spent_group() -> None:
    r"""The domain IS the catalogue entry, so ``invariance_group`` is DERIVED, not declared.

    ⚠ The derived answer is the lower bound :math:`SO(2)_a`; the functions are
    in fact :math:`O(2)_a`-invariant (a mirror across the polar axis does not
    move :math:`\mu`). Declaring that needs an axis-parameterised :math:`O(2)`
    member — GitHub #432. The consequence is recorded, not asserted away:
    a Legendre basis on a :math:`\sigma_b`-fold is refused by G0 today
    although it is mathematically admissible.
    """
    for axis in ("x", "y", "z"):
        basis = LegendreBasis(L=2, axis=axis)
        assert basis.domain == SPHERE.quotient(SubgroupOfO3.SO2(axis))
        assert basis.invariance_group == SubgroupOfO3.SO2(axis)
        assert basis.space == LegendreSpace.from_L(2, axis)
        assert basis.space.shape == (3,)

    with pytest.raises(ValueError, match="axis must be x/y/z"):
        LegendreBasis(L=1, axis="w")
    with pytest.raises(ValueError, match="L must be non-negative"):
        LegendreBasis(L=-1)


# ══════════════════════════════════════════════════════════════════════
# B3 / B3b — the descended column, and its INDEPENDENT reference
# ══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("label", _SPHERE_RULES)
def test_the_descended_column_is_bit_identical_to_the_sh_m0_slot(label: str) -> None:
    r"""B3 — :math:`P_\ell(\Omega\cdot\hat e_x)` IS the harmonic table's :math:`m = 0` column, at the BIT.

    ``[M]`` 2026-09-02 ``max|Δ| = 0.0`` on all 7 shipped sphere rules at
    :math:`L = 4`. That is what makes the repair a bit-identity claim on the
    slab at :math:`L \le 1` rather than a tolerance one.

    ⭐ **The mutation arm is the point.** The two sides could agree because the
    production spelling was echoed rather than because it is right, so this
    row re-spells the table with pure ``lpmv`` — the obvious "clean"
    implementation — and requires that it *breaks* the identity. ``[M]``
    ``lpmv(0, 1, μ)`` differs from :math:`\mu` by ``8.3e-17`` and
    ``eval_legendre`` differs at :math:`\ell \ge 2` by up to ``8.0e-16``: **no
    single scipy routine reproduces the column**, which is why the matched
    branching (``1.0`` / the input / ``lpmv``) is a measured constraint and
    not a taste.
    """
    quadrature = _rule(label)
    L = 4
    directions = np.asarray(quadrature.measure.nodes, dtype=float)

    descended = LegendreBasis(L=L).evaluate(directions)
    m0_column = SphericalHarmonicBasis(L=L).evaluate(directions)[
        :, np.arange(L + 1), np.arange(L + 1)
    ]

    assert np.array_equal(descended, m0_column), (
        f"{label}: max|Δ| = {np.max(np.abs(descended - m0_column)):.3e}"
    )

    # ── the mutation arms: two "obvious" re-spellings, both must MOVE ──
    mu = SPHERE.quotient(SubgroupOfO3.SO2("x")).quotient_map(directions)
    pure_lpmv = np.column_stack([lpmv(0, ell, mu) for ell in range(L + 1)])
    pure_eval = np.column_stack([eval_legendre(ell, mu) for ell in range(L + 1)])

    assert not np.array_equal(pure_lpmv, m0_column), (
        f"{label}: pure lpmv must NOT reproduce the column bit-for-bit — "
        f"if it does, this gate is a convention echo"
    )
    assert not np.array_equal(pure_eval, m0_column), (
        f"{label}: pure eval_legendre must NOT reproduce the column"
    )
    # …and both are numerically fine, which is exactly why a tolerance gate
    # here would be blind: the claim is about BITS.
    np.testing.assert_allclose(pure_lpmv, m0_column, atol=1e-14)
    np.testing.assert_allclose(pure_eval, m0_column, atol=1e-14)


@pytest.mark.parametrize("label", _ONE_D_RULES)
def test_the_slab_table_is_bit_identical_to_the_sh_m0_slot_at_low_order(label: str) -> None:
    r"""B3 on the 1-D rules at :math:`L \le 2` — the population the repair's bit-identity claim rides on.

    The slab's own nodes are cosines, so the comparison is against the
    harmonics evaluated at the ORBIT BARYCENTRE :math:`(\mu, 0, 0)` — the
    forged directions ERR-080 handed them. That construction is what the fix
    retires, and it is built HERE, by hand, precisely because production no
    longer builds it: this row measures that the numbers the slab used to get
    at :math:`\ell \le 1` are the numbers it gets now.

    ⚠ The forged points are OFF :math:`S^2`, so
    :meth:`SphericalHarmonicBasis.evaluate` refuses them (0.6) — the reference
    side here calls the private kernel directly, which is legitimate for a
    HISTORICAL comparison and is the only thing in the tree that still may.
    """
    from orpheus.numerics.basis.spherical_harmonic_basis import _evaluate_real_sh

    quadrature = _rule(label)
    L = 2
    mu = np.asarray(quadrature.measure.nodes, dtype=float).reshape(-1)

    table = LegendreBasis(L=L).evaluate(mu)
    zeros = np.zeros_like(mu)
    forged = _evaluate_real_sh(L, mu, zeros, zeros)[:, np.arange(L + 1), np.arange(L + 1)]

    assert np.array_equal(table, forged), (
        f"{label}: the m=0 column moved — max|Δ| = {np.max(np.abs(table - forged)):.3e}"
    )


def test_legendre_values_against_an_independent_three_term_recurrence() -> None:
    r"""B3b — the CORRECTNESS claim, from a source that shares no ``scipy`` call.

    B3's two sides both route through ``scipy.special.lpmv``. That is legal
    (below the trusted-library line — ``algebra-of-record``), but it means B3
    tests *agreement of conventions*, not *correctness of values*. This row
    supplies the independent leg: Bonnet's recurrence
    :math:`\ell P_\ell = (2\ell-1)\mu P_{\ell-1} - (\ell-1)P_{\ell-2}`,
    evaluated in plain numpy.

    ``[M]`` 2026-09-02 on GL16 at :math:`L = 6`: ``max|Δ| = 4.44e-16``,
    ``nulp <= 32``. Not ``array_equal`` — and it must not be: a recurrence
    accumulates its own rounding, so demanding bits here would pin an
    arithmetic order rather than a value (``vv`` §bit-identity criterion 3,
    drift bounded by reduction depth × ULP).
    """
    mu = np.asarray(_rule("gauss_legendre(16)").measure.nodes, dtype=float).reshape(-1)
    L = 6

    table = legendre_table(L, mu)

    recurrence = np.zeros_like(table)
    recurrence[:, 0] = 1.0
    recurrence[:, 1] = mu
    for ell in range(2, L + 1):
        recurrence[:, ell] = (
            (2 * ell - 1) * mu * recurrence[:, ell - 1] - (ell - 1) * recurrence[:, ell - 2]
        ) / ell

    assert_array_almost_equal_nulp(table, recurrence, nulp=32)

    # NEGATIVE leg: the sign-flipped recurrence is a different family, and
    # this row must be able to tell (vv #11 — otherwise "agrees to nulp 32"
    # is compatible with agreeing with anything).
    wrong = np.zeros_like(table)
    wrong[:, 0] = 1.0
    wrong[:, 1] = mu
    for ell in range(2, L + 1):
        wrong[:, ell] = (
            (2 * ell - 1) * mu * wrong[:, ell - 1] + (ell - 1) * wrong[:, ell - 2]
        ) / ell
    assert not np.allclose(table, wrong, atol=1e-8), (
        "the sign-flipped recurrence must NOT reproduce the table"
    )


# ══════════════════════════════════════════════════════════════════════
# B4 — the contractions, and the dual factor they carry
# ══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("label", ("gauss_legendre(8)", "gauss_legendre(16)"))
@pytest.mark.parametrize("L", (0, 1, 2))
def test_the_contractions_reproduce_the_sh_restriction_on_the_m0_slot(
    label: str, L: int
) -> None:
    r"""B4 — ``analyze`` / ``analyze_transpose`` / ``reconstruct`` restricted to :math:`m=0` are the Legendre ones.

    ``[M]`` 2026-09-02 ``array_equal``, ``max|Δ| = 0.0``, on 6 of 6 (rule,
    :math:`L`) rows. This is the contract that makes the slab's converged
    flux bit-identical across the repair at :math:`L \le 1`.

    ⭐ **The dual-factor arm.** ``reconstruct`` carries :math:`2\ell+1`; drop
    it and the reconstruction is a different operator. The arm below rebuilds
    the contraction with :math:`d_\ell = 1` and requires it to MOVE, so the
    factor is pinned by a red rather than by a restatement of its definition.
    """
    from orpheus.numerics.basis.spherical_harmonic_basis import _evaluate_real_sh

    quadrature = _rule(label)
    mu = np.asarray(quadrature.measure.nodes, dtype=float).reshape(-1)
    weights = np.asarray(quadrature.weights, dtype=float)
    basis = LegendreBasis(L=L)
    table = basis.evaluate(mu)

    zeros = np.zeros_like(mu)
    sh_table = _evaluate_real_sh(L, mu, zeros, zeros)
    diag = (np.arange(L + 1), np.arange(L + 1))

    values = np.sin(3.0 * mu) + 0.25
    coefficients = np.cos(np.arange(L + 1) + 1.0)

    # analyze
    got = basis.analyze(values, table, weights)
    want = np.einsum("n,nlm,n->lm", weights, sh_table, values)[diag]
    assert np.array_equal(got, want), f"analyze max|Δ| = {np.max(np.abs(got - want)):.3e}"

    # analyze_transpose
    got_t = basis.analyze_transpose(coefficients, table, weights)
    padded = np.zeros((L + 1, 2 * L + 1))
    padded[diag] = coefficients
    want_t = np.einsum("n,nlm,lm->n", weights, sh_table, padded)
    assert np.array_equal(got_t, want_t)

    # reconstruct, WITH the 2l+1 dual factor
    got_r = basis.reconstruct(coefficients, table)
    want_r = np.einsum(
        "nlm,l,lm->n", sh_table, SphericalHarmonicBasis(L=L).addition_theorem_factor, padded
    )
    assert np.array_equal(got_r, want_r)

    # NEGATIVE leg — drop the dual factor; the answer must move (except at
    # L = 0, where 2l+1 = 1 and no input can tell: stated, not hidden).
    dropped = np.einsum("nl,l...->n...", table, coefficients)
    if L == 0:
        assert np.array_equal(dropped, got_r), (
            "at L=0 the dual factor is 1 — this arm is STRUCTURALLY blind here"
        )
    else:
        assert not np.allclose(dropped, got_r), (
            "dropping the (2l+1) dual factor must move reconstruct"
        )

    # synthesize is the naked (factor-free) partner — pin the pair's
    # difference so neither can drift into the other.
    assert np.array_equal(basis.synthesize(coefficients, table), dropped)


def test_reconstruct_transpose_is_the_transpose_of_reconstruct() -> None:
    r"""The representation transpose is the matrix transpose — an algebraic law, not a spelling.

    :math:`\langle R\phi, v\rangle = \langle \phi, R^\top v\rangle` for every
    pair, which pins the dual factor's PLACEMENT (both faces carry it) rather
    than merely its value.
    """
    quadrature = _rule("gauss_legendre(8)")
    basis = LegendreBasis(L=3)
    table = basis.evaluate(np.asarray(quadrature.measure.nodes, dtype=float).reshape(-1))

    rng = np.random.default_rng(20260902)
    phi = rng.normal(size=4)
    v = rng.normal(size=8)

    lhs = float(basis.reconstruct(phi, table) @ v)
    rhs = float(phi @ basis.reconstruct_transpose(v, table))
    np.testing.assert_allclose(lhs, rhs, rtol=0.0, atol=1e-13)

    # NEGATIVE: a factor dropped on one side only breaks it.
    broken = np.einsum("nl,n...->l...", table, v)
    assert not np.isclose(lhs, float(phi @ broken))
