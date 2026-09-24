r"""Branch-1 gates: the octant face transmission of step, diamond and LD.

System under test: :mod:`orpheus.derivations.discrete.sn.face_transmission`,
the algebra of record in which each closure is a Petrov-Galerkin choice
(trial, test and face spaces, weak or strong inflow) integrated on the
reference cell, with lengths in mean free paths (:math:`\Sigma_t = 1`), so the
cell enters only through :math:`g_a = |\mu_a|/(\Sigma_t\Delta_a) = 1/\tau_a`.

One or more tests per claim of the module (its docstring numbers them 1 to 8),
plus the two certifiers every spectral claim rests on.  Each test runs the
module's own function, which raises ``DerivationFailed`` when a proof step does
not hold, and then checks the returned object against a reference the module
does not compute:

* the certifiers: known-true and known-false polynomials, both legs;
* the realization (1) and the module's ``transmission()``: the Kronecker UBLD of
  :mod:`~orpheus.derivations.discrete.sn.ld_ubld` composed here with plain
  ``LUsolve`` at rational points, and the literature 1-D forms;
* the spectrum (2, 3): a direct exact characteristic polynomial of :math:`T`
  at a rational point against the factors, and a numeric root sweep over six
  decades of :math:`g` against the Jury certificates;
* the conserved mode (4) and the stability function (6): the literature 1-D
  transmissions, :func:`mpmath.pade`, and a numeric :math:`|R| < 1` sweep;
* the page's closed forms (8): the page's formulas in physical variables,
  evaluated in :class:`fractions.Fraction`;
* step at :math:`d = 1`: a hand balance in :class:`fractions.Fraction`, which
  the retired ``derivations/sn_dd_face_transmission.py`` fails (it built step
  with diamond's weights, :math:`2/(2+\tau)`).

Markers: every test is ``foundation`` (a SymPy identity of the algebra of
record).  The diamond rows that state what the label
``dd-face-transmission-spectrum`` states (:math:`(2/D)\mathbf 1 w^{\mathsf T}
- I`, spectrum :math:`\{1 - 2\Sigma_t V/D\} \cup \{-1\}^{d-1}`) also carry
``verifies``: the algebra-of-record Branch-1 exception of ``vv-principles``.

Runtime: about 22 s per process, of which LD's :math:`d = 3` factorisation is
about 20 s (cached by ``functools.cache``).  Not marked ``slow``: the LD
:math:`d = 3` rows are the only evidence that LD damps every face mode in three
dimensions, where production cannot yet answer (#503), and a ``slow`` mark
would deselect them from every ``-m "not slow"`` run.
"""
from __future__ import annotations

import itertools
from fractions import Fraction

import mpmath
import numpy as np
import pytest
import sympy as sp

from orpheus.derivations.discrete.sn.face_transmission import (
    LAMBDA,
    Closure,
    SCHEMES,
    DerivationFailed,
    Imposition,
    PetrovGalerkinScheme,
    _weak_form_response,
    _is_a_stable,
    _is_positive_on_orthant,
    _jury_strictly_inside,
    cell_response,
    derive_comparison,
    derive_conserved_mode,
    derive_flat_flux_preservation,
    derive_ld_factorisation,
    derive_page_closed_form,
    derive_particle_balance,
    derive_realization,
    derive_stability_function,
    streaming_coefficients,
    transmission_spectrum,
)
from orpheus.derivations.discrete.sn.ld_ubld import (
    assemble_ubld,
    inflow_scatter_axis,
    outflow_trace_axis,
)

LABEL = "dd-face-transmission-spectrum"
DIMS = (1, 2, 3)
UPWIND = (Closure.STEP, Closure.LINEAR_DISCONTINUOUS)
ALL = tuple(Closure)
LD = Closure.LINEAR_DISCONTINUOUS

#: Generic rational points: distinct, non-integer, no two coordinates equal (Mode 2).
_POINTS = (
    (sp.Rational(7, 10), sp.Rational(13, 11), sp.Rational(2, 5)),
    (sp.Rational(1, 50), sp.Rational(9, 2), sp.Rational(3, 7)),
    (sp.Rational(40, 3), sp.Rational(1, 9), sp.Rational(5, 4)),
)
_POINT = _POINTS[0]

#: Six decades of streaming per axis: optically thick (g = 1e-3) to thin (g = 1e3).
_GRID = (sp.Rational(1, 1000), sp.Rational(1, 10), sp.Integer(1), sp.Integer(10), sp.Integer(1000))


def _fail_unless(condition: bool, message: str) -> None:
    """``pytest.fail`` unless ``condition`` holds; survives ``python -O``."""
    if not condition:
        pytest.fail(message)


def _at(d: int, point=_POINT) -> dict:
    return dict(zip(streaming_coefficients(d), point[:d]))


def _is_zero_matrix(M: sp.Matrix) -> bool:
    return M.applyfunc(sp.cancel) == sp.zeros(*M.shape)


def _expected_undamped(closure: Closure, d: int) -> tuple:
    return ((sp.Integer(-1), d - 1),) if closure is Closure.DIAMOND and d > 1 else ()


def _ubld_transmission(d: int, point) -> sp.Matrix:
    r"""LD's :math:`T` from :mod:`ld_ubld`'s Kronecker UBLD, composed here with plain ``LUsolve``.

    The module no longer builds LD from ``ld_ubld``; this is the independent
    route (``ld_ubld`` carries its own d=1 reduction and exact-on-bilinear oracles).
    """
    ones: list[sp.Expr] = [sp.Integer(1)] * d
    mus: list[sp.Expr] = list(point[:d])
    theta = sp.Rational(1, 3)
    A = assemble_ubld(ones, mus, sp.Integer(1), theta)["A"]
    E = sp.Matrix.hstack(*(inflow_scatter_axis(ones, mus, a, theta) for a in range(d)))
    R = sp.Matrix.vstack(*(outflow_trace_axis(d, a) for a in range(d)))
    return R * A.LUsolve(E)


# ── the certifiers: each must be able to say False ───────────────────────

_x, _y = sp.symbols("x y", positive=True)


@pytest.mark.foundation
@pytest.mark.parametrize(
    ("p", "certified"),
    [
        pytest.param((_x - _y) ** 2, False, id="square-vanishes-on-diagonal"),
        pytest.param(_x**2 - 3 * _x * _y + _y**2 + _x + 1, False, id="indefinite-quadratic"),
        pytest.param(-_x, False, id="negative"),
        pytest.param(sp.Integer(0), False, id="zero"),
        pytest.param(_x * _y - _x, False, id="negative-linear-term"),
        pytest.param((_x - _y) ** 2 + _x, True, id="semidefinite-plus-positive"),
        pytest.param(12 * _x**2 - 20 * _x * _y + 12 * _y**2 + 2 * _x + 2 * _y + 1, True, id="definite"),
        pytest.param(_x * _y + _x + 3, True, id="non-negative-coefficients"),
    ],
)
def test_positivity_certifier_says_false_and_true(p, certified) -> None:
    """``_is_positive_on_orthant`` refuses a certificate it cannot give, and gives the ones it can.

    ``False`` means "no certificate", so every False row is a polynomial that
    is not positive on the open orthant or that no certificate covers; the
    True rows are positive and inside the certificate's scope.
    """
    got = _is_positive_on_orthant(p, (_x, _y))
    _fail_unless(got is certified, f"certificate({p}) = {got}, expected {certified}")


@pytest.mark.foundation
@pytest.mark.parametrize(
    ("factor", "inside"),
    [
        pytest.param(LAMBDA + 1, False, id="root-at-minus-1"),
        pytest.param(LAMBDA**2 + 1, False, id="roots-at-plus-minus-i"),
        pytest.param(LAMBDA**2 - 1, False, id="roots-at-plus-minus-1"),
        pytest.param((_x + 1) * LAMBDA - 2 * _x, False, id="root-2x-over-x-plus-1-exceeds-1"),
        pytest.param((_x + 1) * LAMBDA - _x, True, id="root-x-over-x-plus-1"),
        pytest.param((2 * _x + 1) * LAMBDA - 2 * _x + 1, True, id="diamond-1d-conserved"),
        pytest.param(4 * LAMBDA**2 + LAMBDA + 1, True, id="quadratic-roots-modulus-one-half"),
        pytest.param(LAMBDA**3, None, id="cubic-refused"),
    ],
)
def test_jury_certifier_says_false_and_true(factor, inside) -> None:
    """``_jury_strictly_inside`` refuses roots on or outside the unit disk."""
    if inside is None:
        with pytest.raises(NotImplementedError, match="degree <= 2"):
            _jury_strictly_inside(sp.Poly(factor, LAMBDA, _x, _y), (_x, _y))
        return
    got = _jury_strictly_inside(sp.Poly(factor, LAMBDA, _x, _y), (_x, _y))
    _fail_unless(got is inside, f"jury({factor}) = {got}, expected {inside}")


_z = sp.Symbol("z")


@pytest.mark.foundation
@pytest.mark.parametrize(
    ("numerator", "denominator", "a_stable"),
    [
        pytest.param(1 + _z / 2, 1 - _z / 2, False, id="pole-at-z-2-not-hurwitz"),
        pytest.param(sp.Integer(1), 1 + _z**2, False, id="poles-on-the-imaginary-axis"),
        pytest.param(1 - _z + _z**2 / 2, sp.Integer(1), False, id="hurwitz-but-E-is-minus-y4-over-4"),
        pytest.param(1 - 2 * _z / 3 + _z**2 / 6, 1 + _z / 3, False, id="pade-2-1-E-negative"),
        pytest.param(sp.Integer(-1), -1 - _z, True, id="hurwitz-with-all-negative-coefficients"),
        pytest.param(sp.Integer(1), 1 + _z + _z**2 / 2, True, id="pade-0-2"),
        pytest.param(sp.Integer(1), 1 + _z + _z**2 + _z**3, None, id="cubic-denominator-refused"),
    ],
)
def test_a_stability_certifier_says_false_and_true(numerator, denominator, a_stable) -> None:
    """``_is_a_stable`` refuses a pole in the closed right half-plane and a negative E(y), and certifies the rest.

    ``pole-at-z-2`` fails the Hurwitz leg alone (its E(y) is 0); the two
    ``E`` rows fail the E-polynomial leg with a Hurwitz denominator; so each
    leg has a row only it refuses.  ``poles-on-the-imaginary-axis`` fails both
    (a zero coefficient, and E = y^4 - 2y^2).
    """
    if a_stable is None:
        with pytest.raises(NotImplementedError, match="degree <= 2"):
            _is_a_stable(numerator, denominator, _z)
        return
    got = _is_a_stable(numerator, denominator, _z)
    _fail_unless(got is a_stable, f"A-stable({numerator}/{denominator}) = {got}, expected {a_stable}")


@pytest.mark.foundation
@pytest.mark.parametrize("closure", ALL, ids=lambda c: c.value)
def test_a_stability_certifier_certifies_each_closures_own_transmission(closure: Closure) -> None:
    """True on each closure's own 1-D P and Q, written from the literature (E = y^2, 0, y^4)."""
    tau = sp.Symbol("tau", positive=True)
    numerator, denominator = sp.fraction(sp.cancel(_one_dimensional_reference(closure, tau)))
    _fail_unless(_is_a_stable(numerator, denominator, tau) is True, f"{closure.value}: {numerator}/{denominator}")


# ── claim 1: the realizations ────────────────────────────────────────────

_SHAPE = {
    Closure.STEP: (lambda d: d, lambda d: 1, Imposition.WEAK, 0),
    Closure.DIAMOND: (lambda d: d, lambda d: 1, Imposition.STRONG, -1),
    LD: (lambda d: d * 2 ** (d - 1), lambda d: 2**d, Imposition.WEAK, 0),
}


@pytest.mark.foundation
@pytest.mark.parametrize("d", DIMS)
@pytest.mark.parametrize("closure", ALL, ids=lambda c: c.value)
def test_weak_form_realizes_the_known_closure(closure: Closure, d: int) -> None:
    """The weak form's (A, E, S, R, D) is the known closure: LD is ld_ubld's UBLD, step and diamond their balances.

    Also pins the sizes and the feedthrough: :math:`D = 0` for the weak
    closures, :math:`D = -I` for diamond (strong imposition creates it).
    """
    response = derive_realization(closure, d)
    n_face, n_cell, imposition, f = _SHAPE[closure]
    _fail_unless(closure.scheme(d).inflow is imposition, f"{closure.value}: imposition {closure.scheme(d).inflow}")
    _fail_unless(
        (response.n_face_moments, response.n_cell_moments) == (n_face(d), n_cell(d)),
        f"{closure.value}, d={d}: {response.n_face_moments} face and {response.n_cell_moments} cell moments",
    )
    _fail_unless(
        response.feedthrough_scalar == f and response.feedthrough == f * sp.eye(n_face(d)),
        f"{closure.value}, d={d}: feedthrough {response.feedthrough}",
    )


@pytest.mark.foundation
def test_feedthrough_scalar_refuses_a_non_scalar_feedthrough() -> None:
    """``feedthrough_scalar`` is read only when D = f I; a non-scalar D is refused as a scope edge, not truncated."""
    import dataclasses

    response = cell_response(Closure.DIAMOND, 2)
    skewed = dataclasses.replace(response, feedthrough=sp.Matrix([[-1, 0], [0, sp.Rational(-1, 2)]]))
    with pytest.raises(NotImplementedError, match="not a rational multiple of the identity"):
        _ = skewed.feedthrough_scalar


# ── the family's admission: the scheme type, the table, the condensation ──


@pytest.mark.foundation
@pytest.mark.parametrize(
    ("trial", "test", "face", "fragment"),
    [
        pytest.param(((0, 0), (1,)), ((0, 0),), ((0,),), "distinct multi-indices of length 2", id="trial-wrong-length"),
        pytest.param(((0, 0), (0, 0)), ((0, 0),), ((0,),), "distinct multi-indices of length 2", id="trial-duplicate"),
        pytest.param(((0, 0),), ((0, 0),), ((0, 0),), "the face space is distinct multi-indices of length 1", id="face-wrong-length"),
        pytest.param(
            ((0, 0), (1, 0), (0, 1), (1, 1)), ((0, 0), (1, 0), (0, 1), (1, 1)), ((1,), (0,)),
            "face average first", id="face-not-average-first",
        ),
    ],
)
def test_scheme_refuses_malformed_spaces(trial, test, face, fragment) -> None:
    """``PetrovGalerkinScheme`` refuses a space of wrong-length or repeated multi-indices, or a face space not led by its average.

    The face-first row is LD's own spaces with the face space reordered, so
    only the average-first rule can refuse it.
    """
    with pytest.raises(DerivationFailed, match=fragment):
        PetrovGalerkinScheme(2, trial, test, face, Imposition.WEAK)


@pytest.mark.foundation
def test_scheme_admits_every_tabulated_closure() -> None:
    """The positive leg: every row of ``SCHEMES`` constructs at d = 1, 2, 3, and the table covers ``Closure``."""
    _fail_unless(set(SCHEMES) == set(Closure), f"SCHEMES rows {set(SCHEMES)}")
    for closure in Closure:
        for d in DIMS:
            scheme = closure.scheme(d)
            _fail_unless(scheme.d == d and scheme.face[0] == (0,) * (d - 1), f"{closure.value}, d={d}: {scheme}")


@pytest.mark.foundation
def test_a_closure_with_no_table_row_raises_key_error(monkeypatch) -> None:
    """``Closure.scheme`` reads ``SCHEMES``; a closure the table does not list has no scheme."""
    monkeypatch.delitem(SCHEMES, Closure.STEP)
    with pytest.raises(KeyError):
        Closure.STEP.scheme(1)


#: Passes __post_init__ (3 trial - 1 test = 2 = d x 1 constrained coordinates),
#: but both constraints read (1,0)'s trace and (1,1) is invisible on both
#: inflow faces' averages, so C_c is singular.  Before the Schur-complement
#: condensation it was accepted with D = diag(-1, 0) (the elegance-enforcer's witness).
_UNDERCONSTRAINED = PetrovGalerkinScheme(2, ((0, 0), (1, 0), (1, 1)), ((0, 0),), ((0,),), Imposition.STRONG)


@pytest.mark.foundation
def test_condensation_refuses_constraints_that_miss_a_constrained_coordinate(monkeypatch) -> None:
    """``_static_condensation`` refuses a singular C_c instead of returning a wrong feedthrough."""
    monkeypatch.setattr(Closure, "scheme", lambda self, d: _UNDERCONSTRAINED)
    with pytest.raises(DerivationFailed, match="C_c square, invertible"):
        _weak_form_response(Closure.DIAMOND, 2)


@pytest.mark.foundation
def test_condensation_admits_diamonds_constraints() -> None:
    """The positive leg: diamond's strong constraints pass the same check, and yield D = -I."""
    response = _weak_form_response(Closure.DIAMOND, 2)
    _fail_unless(response.feedthrough == -sp.eye(2), f"D = {response.feedthrough}")


# ── the module's own transmission(), pinned ──────────────────────────────


@pytest.mark.foundation
def test_ld_1d_transmission_is_the_literature_form() -> None:
    r"""The symbolic ``transmission()`` of LD at d=1 is :math:`(6g^2 - 2g)/(6g^2 + 4g + 1)`.

    That is :math:`(6 - 2\tau)/(6 + 4\tau + \tau^2)` at :math:`g = 1/\tau`
    (Morel, Wareing & Smith 1996, JCP 128:445, Eq. (74), their Padé (1,2)).
    """
    (g,) = streaming_coefficients(1)
    T = cell_response(LD, 1).transmission()
    _fail_unless(T.shape == (1, 1), f"shape {T.shape}")
    _fail_unless(
        sp.cancel(T[0, 0] - (6 * g**2 - 2 * g) / (6 * g**2 + 4 * g + 1)) == 0,
        f"LD d=1 transmission {T[0, 0]}",
    )


@pytest.mark.foundation
@pytest.mark.parametrize("point", _POINTS, ids=("p0", "p1", "p2"))
def test_ld_2d_symbolic_transmission_equals_the_ubld_composition(point) -> None:
    """The symbolic d=2 ``transmission()``, evaluated at a point, and ``transmission(at=)`` equal ld_ubld's T there."""
    response = cell_response(LD, 2)
    at = _at(2, point)
    reference = _ubld_transmission(2, point)
    _fail_unless(
        _is_zero_matrix(response.transmission().subs(at) - reference),
        f"LD d=2 symbolic transmission at {point[:2]} differs from ld_ubld's",
    )
    _fail_unless(
        _is_zero_matrix(response.transmission(at) - reference),
        f"LD d=2 transmission(at) at {point[:2]} differs from ld_ubld's",
    )


@pytest.mark.foundation
@pytest.mark.parametrize("point", _POINTS[:2], ids=("p0", "p1"))
def test_ld_3d_transmission_at_a_point_equals_the_ubld_composition(point) -> None:
    """At d=3 at a point only: the symbolic block costs about 50 s."""
    _fail_unless(
        _is_zero_matrix(cell_response(LD, 3).transmission(_at(3, point)) - _ubld_transmission(3, point)),
        f"LD d=3 transmission(at) at {point} differs from ld_ubld's",
    )


# ── claims 2 and 3: the spectrum ─────────────────────────────────────────


def _check_spectrum(closure: Closure, d: int) -> None:
    """Claims 2 and 3 at one (closure, d), with two independent witnesses."""
    spectrum = transmission_spectrum(closure, d)
    _fail_unless(
        spectrum.undamped == _expected_undamped(closure, d),
        f"{closure.value}, d={d}: undamped {spectrum.undamped}, expected {_expected_undamped(closure, d)}",
    )
    response = cell_response(closure, d)
    _fail_unless(
        spectrum.feedthrough_multiplicity >= response.n_face_moments - response.n_cell_moments,
        f"{closure.value}, d={d}: Sylvester's bound on the feedthrough multiplicity",
    )
    g = response.g
    at = _at(d)

    # Witness 1: the exact characteristic polynomial of T at a rational point,
    # by Berkowitz on the n x n matrix (no Sylvester shortcut), equals the
    # product of the factors up to a lambda-free constant.
    direct = response.transmission(at).charpoly(LAMBDA.name).as_expr()
    product = sp.Mul(*(fac.as_expr().subs(list(at.items())) ** mult for fac, mult in spectrum.factors))
    ratio = sp.cancel(direct / product)
    _fail_unless(
        ratio != 0 and not ratio.has(LAMBDA),
        f"{closure.value}, d={d}: det(lambda I - T) / prod(factors) = {ratio} depends on lambda",
    )

    # Witness 2: every non-feedthrough root over a six-decade grid is strictly
    # inside the unit disk.  Independent of the positivity certificates.
    f = response.feedthrough_scalar
    others = [fac for fac, _ in spectrum.factors if sp.expand(fac.as_expr() - (LAMBDA - f)) != 0
              and sp.expand(fac.as_expr() + (LAMBDA - f)) != 0]
    coefficient_fns = [sp.lambdify(g, sp.Poly(fac.as_expr(), LAMBDA).all_coeffs(), "mpmath") for fac in others]
    worst = 0.0
    for point in itertools.product(_GRID, repeat=d):
        for coefficients in coefficient_fns:
            roots = np.roots([complex(c) for c in coefficients(*point)])
            worst = max(worst, float(np.abs(roots).max()) if roots.size else 0.0)
    _fail_unless(worst < 1.0, f"{closure.value}, d={d}: a non-feedthrough root reaches |lambda| = {worst!r}")


@pytest.mark.foundation
@pytest.mark.verifies(LABEL)
@pytest.mark.parametrize("d", DIMS)
def test_diamond_carries_exactly_d_minus_1_undamped_face_modes(d: int) -> None:
    r"""Diamond: :math:`-1` with multiplicity :math:`d-1`, every other root strictly inside."""
    _check_spectrum(Closure.DIAMOND, d)


@pytest.mark.foundation
@pytest.mark.parametrize("d", DIMS)
@pytest.mark.parametrize("closure", UPWIND, ids=lambda c: c.value)
def test_step_and_ld_damp_every_face_mode(closure: Closure, d: int) -> None:
    r"""Step and LD (weak, no feedthrough): every root strictly inside, for every :math:`g > 0`."""
    _check_spectrum(closure, d)


@pytest.mark.foundation
def test_derive_comparison_is_the_table_of_the_nine_rows() -> None:
    """The module's own comparison, as published: diamond alone leaves face modes undamped."""
    table = derive_comparison(DIMS)
    _fail_unless(set(table) == {(c, d) for c in Closure for d in DIMS}, f"covers {sorted((c.value, d) for c, d in table)}")
    for (closure, d), spectrum in table.items():
        _fail_unless(spectrum.undamped == _expected_undamped(closure, d), f"{closure.value}, d={d}: {spectrum.undamped}")
        _fail_unless(
            spectrum.damps_every_face_mode is (closure is not Closure.DIAMOND or d == 1),
            f"{closure.value}, d={d}: damps_every_face_mode = {spectrum.damps_every_face_mode}",
        )


# ── claim 4: the conserved mode ──────────────────────────────────────────


def _one_dimensional_reference(closure: Closure, tau: sp.Expr) -> sp.Expr:
    r"""The 1-D transmissions typed from the literature, not from the module.

    Step :math:`1/(1+\tau)`; diamond :math:`(2-\tau)/(2+\tau)`; LD (Lewis and
    Miller 1989) :math:`(6-2\tau)/(6+4\tau+\tau^2)`.
    """
    return {
        Closure.STEP: 1 / (1 + tau),
        Closure.DIAMOND: (2 - tau) / (2 + tau),
        LD: (6 - 2 * tau) / (6 + 4 * tau + tau**2),
    }[closure]


def _check_conserved(closure: Closure, d: int) -> None:
    g = streaming_coefficients(d)
    conserved = derive_conserved_mode(closure, d)
    expected = _one_dimensional_reference(closure, sp.Integer(1) / sp.Add(*g))
    _fail_unless(
        sp.cancel(conserved - expected) == 0,
        f"{closure.value}, d={d}: conserved eigenvalue {conserved}, expected T_1D(1/G) = {expected}",
    )


@pytest.mark.foundation
@pytest.mark.verifies(LABEL)
@pytest.mark.parametrize("d", DIMS)
def test_diamond_conserved_mode_is_one_minus_two_sigma_v_over_d(d: int) -> None:
    r"""Diamond's absorption-damped eigenvalue is :math:`1 - 2\Sigma_t V/D = 1 - 2/(1 + 2G)`."""
    _check_conserved(Closure.DIAMOND, d)
    total = sp.Add(*streaming_coefficients(d))
    _fail_unless(
        sp.cancel(derive_conserved_mode(Closure.DIAMOND, d) - (1 - 2 / (1 + 2 * total))) == 0,
        f"diamond, d={d}: the conserved eigenvalue is not the label's 1 - 2 Sigma_t V / D",
    )


@pytest.mark.foundation
@pytest.mark.parametrize("d", DIMS)
@pytest.mark.parametrize("closure", UPWIND, ids=lambda c: c.value)
def test_step_and_ld_conserved_mode_is_their_1d_transmission(closure: Closure, d: int) -> None:
    """The factor lambda - T_1D(G) of the characteristic polynomial, with T_1D the closure's own."""
    _check_conserved(closure, d)


# ── claim 5: LD's spectrum is generated by d = 1 and d = 2 ───────────────


def _lambda_degree(fac: sp.Poly) -> int:
    return int(sp.Poly(fac.as_expr(), LAMBDA).degree())


@pytest.mark.foundation
@pytest.mark.parametrize(("d", "zero_multiplicity", "n_quadratics"), [(2, 1, 1), (3, 5, 3)])
def test_ld_factorisation_from_its_low_dimensions(d: int, zero_multiplicity: int, n_quadratics: int) -> None:
    r""":math:`\lambda^z\,\ell(G)\prod q(g_a, G - g_a)`, with z, the quadratic count and the degree pinned."""
    factors = derive_ld_factorisation(d)
    zero = sp.Poly(LAMBDA, LAMBDA, *streaming_coefficients(d))
    _fail_unless(factors.get(zero) == zero_multiplicity, f"d={d}: lambda^{factors.get(zero)}")
    _fail_unless(
        sum(m for fac, m in factors.items() if _lambda_degree(fac) == 2) == n_quadratics, f"d={d}: {factors}",
    )
    _fail_unless(
        sum(_lambda_degree(fac) * m for fac, m in factors.items()) == d * 2 ** (d - 1),
        f"d={d}: total degree is not the {d * 2 ** (d - 1)} face moments",
    )


@pytest.mark.foundation
def test_ld_factorisation_refuses_an_unclaimed_dimension() -> None:
    """The claim is stated for d in {2, 3}; d = 4 is refused, not guessed."""
    with pytest.raises(ValueError, match="claimed for d in"):
        derive_ld_factorisation(4)


# ── claim 6: the stability function ──────────────────────────────────────

_STABILITY = {
    Closure.STEP: ((0, 1), 0, True),
    Closure.DIAMOND: ((1, 1), -1, False),
    LD: ((1, 2), 0, True),
}


@pytest.mark.foundation
@pytest.mark.parametrize("closure", ALL, ids=lambda c: c.value)
def test_stability_function_is_an_a_stable_pade_approximant(closure: Closure) -> None:
    r"""Padé type, :math:`|R| < 1`, :math:`R(\infty) = f`, L-stability, and :math:`T \to fI` as :math:`g \to 0`.

    References: the literature form, :func:`mpmath.pade` of :math:`e^{-\tau}`
    (a different algorithm from the module's series check), and a numeric
    :math:`|R(\tau)| < 1` sweep over :math:`\tau \in [10^{-6}, 10^{6}]`.
    """
    sf = derive_stability_function(closure)
    (p, q), at_infinity, l_stable = _STABILITY[closure]
    tau = sp.Symbol("tau", positive=True)
    _fail_unless(sf.closure is closure and sf.pade_order == (p, q), f"{closure.value}: order {sf.pade_order}")
    _fail_unless(sf.value_at_infinity == at_infinity, f"{closure.value}: R(oo) = {sf.value_at_infinity}")
    _fail_unless(sf.l_stable is l_stable, f"{closure.value}: l_stable = {sf.l_stable}")
    symbols = tuple(sf.transmission.free_symbols)
    _fail_unless(len(symbols) == 1, f"{closure.value}: R = {sf.transmission} is not a function of tau alone")
    _fail_unless(
        sp.cancel(sf.transmission.subs(symbols[0], tau) - _one_dimensional_reference(closure, tau)) == 0,
        f"{closure.value}: R = {sf.transmission}, the literature form is {_one_dimensional_reference(closure, tau)}",
    )
    mpmath.mp.dps = 30
    taylor = [mpmath.mpf(-1) ** k / mpmath.factorial(k) for k in range(p + q + 1)]
    num_ref, den_ref = mpmath.pade(taylor, p, q)
    num, den = (sp.Poly(part, tau) for part in sp.fraction(sp.cancel(_one_dimensional_reference(closure, tau))))
    scale = den.all_coeffs()[-1]
    got = [num.coeff_monomial(tau**k) / scale for k in range(p + 1)] + [den.coeff_monomial(tau**k) / scale for k in range(q + 1)]
    for coefficient, reference in zip(got, list(num_ref) + list(den_ref)):
        _fail_unless(abs(float(coefficient) - float(reference)) < 1e-25, f"{closure.value}: Pade coefficient {coefficient} != {reference}")
    R = sp.lambdify(symbols, sf.transmission, "mpmath")
    worst = max(abs(R(mpmath.mpf(10) ** k)) for k in np.linspace(-6, 6, 49))
    _fail_unless(worst < 1, f"{closure.value}: |R(tau)| reaches {worst}")


@pytest.mark.foundation
@pytest.mark.parametrize(
    ("order", "fake", "fragment"),
    [
        pytest.param((2, 0), lambda t: 1 - t + t**2 / 2, "E\\(y\\)", id="pade-2-0-E-is-minus-y4-over-4"),
        pytest.param((2, 1), lambda t: (1 - 2 * t / 3 + t**2 / 6) / (1 + t / 3), "E\\(y\\)", id="pade-2-1-outside-the-a-stable-band"),
        pytest.param((1, 1), lambda t: (1 + t / 2) / (1 - t / 2), "agrees with exp", id="pole-at-tau-2-refused-by-the-pade-leg"),
    ],
)
def test_stability_function_refuses_a_non_a_stable_transmission(monkeypatch, order, fake, fragment) -> None:
    r"""A 1-D transmission that is not A-stable is refused, and the row names the leg that refuses it.

    The first two are genuine Padé approximants of :math:`e^{-\tau}` that pass
    the Padé legs, so only :math:`E(y) \ge 0` can refuse them (:math:`[2/0]`:
    :math:`E = -y^4/4`; :math:`[2/1]` lies outside the band
    :math:`q-2 \le p \le q`).  The third, with its pole at :math:`\tau = 2`, is
    refused one leg earlier, by the Padé agreement.  No input that passes the
    Padé legs with :math:`q \le 2` can fail the Hurwitz leg, since every Padé
    denominator of :math:`e^{-\tau}` has positive coefficients; that leg is gated
    directly on ``_is_a_stable`` (``test_a_stability_certifier_says_false_and_true``).
    """
    import orpheus.derivations.discrete.sn.face_transmission as ft

    monkeypatch.setattr(ft, "_one_dimensional_transmission", lambda closure, g: sp.cancel(fake(1 / g)))
    monkeypatch.setitem(ft.PADE_ORDER, Closure.STEP, order)
    with pytest.raises(DerivationFailed, match=fragment):
        derive_stability_function(Closure.STEP)


# ── claim 7: conservation and flat-flux preservation ─────────────────────


def _balance_point(closure: Closure, d: int):
    return _at(d) if (closure is LD and d == 3) else None


@pytest.mark.foundation
@pytest.mark.parametrize("d", DIMS)
@pytest.mark.parametrize("closure", ALL, ids=lambda c: c.value)
def test_particle_balance_from_the_four_blocks(closure: Closure, d: int) -> None:
    """Outflow minus inflow plus absorption equals the source, for every inflow and source moment.

    Symbolic in g, except LD at d=3, proved at a rational point.  The four
    blocks enter: transmission, escape, inflow to cell, source to cell.
    """
    derive_particle_balance(closure, d, _balance_point(closure, d))


@pytest.mark.foundation
@pytest.mark.parametrize("d", DIMS)
@pytest.mark.parametrize("closure", ALL, ids=lambda c: c.value)
def test_flat_flux_passes_through_the_cell(closure: Closure, d: int) -> None:
    """A uniform inflow equal to a uniform source leaves the outflow and the cell state uniform."""
    derive_flat_flux_preservation(closure, d, _balance_point(closure, d))


# ── claim 8: the page's closed forms ─────────────────────────────────────

_WIDTHS = (Fraction(7, 10), Fraction(13, 10), Fraction(2, 5))
_COSINES = (Fraction(1, 3), Fraction(2, 3), Fraction(2, 3))
_SIGMA_T = Fraction(9, 10)


def _page_matrix(closure: Closure, d: int) -> list[list[Fraction]]:
    r"""The page's :math:`\Sigma` in physical variables, in :class:`Fraction`.

    Diamond :math:`(2/D)\mathbf 1 w^{\mathsf T} - I`, :math:`w_a = 2|\mu_a|A_a`,
    :math:`D = \Sigma_t V + \sum_b w_b`; step :math:`(1/D')\mathbf 1 w'^{\mathsf T}`,
    :math:`w'_a = |\mu_a|A_a`, :math:`D' = \Sigma_t V + \sum_b w'_b`;
    :math:`A_a = V/\Delta_a`.
    """
    volume = Fraction(1)
    for width in _WIDTHS[:d]:
        volume *= width
    area = [volume / _WIDTHS[a] for a in range(d)]
    factor, minus_identity = (2, 1) if closure is Closure.DIAMOND else (1, 0)
    w = [factor * _COSINES[a] * area[a] for a in range(d)]
    D = _SIGMA_T * volume + sum(w)
    return [[factor * w[b] / D - (minus_identity if a == b else 0) for b in range(d)] for a in range(d)]


def _check_closed_form(closure: Closure, d: int) -> None:
    """One closure against its page form, through the symbolic and the point path."""
    derive_page_closed_form(closure, d)
    response = cell_response(closure, d)
    at = dict(zip(response.g, [sp.Rational(_COSINES[a] / (_SIGMA_T * _WIDTHS[a])) for a in range(d)]))
    page = sp.Matrix([[sp.Rational(x) for x in row] for row in _page_matrix(closure, d)])
    for path, T in (("symbolic", response.transmission().subs(at)), ("at", response.transmission(at))):
        _fail_unless(_is_zero_matrix(T - page), f"{closure.value}, d={d}, {path}: T = {T}, page = {page}")


@pytest.mark.foundation
@pytest.mark.verifies(LABEL)
@pytest.mark.parametrize("d", DIMS)
def test_diamond_transmission_is_the_page_closed_form(d: int) -> None:
    r""":math:`T_{\rm DD} = (2/D)\mathbf 1 w^{\mathsf T} - I`, in the page's physical variables."""
    _check_closed_form(Closure.DIAMOND, d)


@pytest.mark.foundation
@pytest.mark.parametrize("d", DIMS)
def test_step_transmission_is_the_page_closed_form(d: int) -> None:
    r""":math:`T_{\rm step} = (1/D')\mathbf 1 w'^{\mathsf T}`, with step's own weights."""
    _check_closed_form(Closure.STEP, d)


@pytest.mark.foundation
def test_page_closed_form_refuses_ld() -> None:
    """The page states no closed form for LD; the claim refuses it rather than inventing one."""
    with pytest.raises(ValueError, match="step and diamond only"):
        derive_page_closed_form(LD, 2)


# ── step at d = 1: the hand balance the predecessor failed ───────────────


@pytest.mark.foundation
@pytest.mark.catches("ERR-088")
@pytest.mark.parametrize("tau", [Fraction(1, 7), Fraction(1), Fraction(9, 4), Fraction(40)])
def test_step_1d_is_one_over_one_plus_tau_not_the_predecessors_two_over_two_plus_tau(tau: Fraction) -> None:
    r"""Step at :math:`d = 1` is :math:`1/(1+\tau)`, from the balance solved by hand in :class:`Fraction`.

    :math:`g(\psi^{\rm out} - \psi^{\rm in}) + \psi_c = 0` with :math:`\psi^{\rm out} = \psi_c`
    gives :math:`\psi^{\rm out} = g/(g+1)\,\psi^{\rm in}`.  The retired
    predecessor read :math:`2/(2+\tau)`; the first assertion keeps the row
    able to tell the two apart.
    """
    g = 1 / tau
    by_hand = g / (g + 1)
    predecessor = Fraction(2) / (2 + tau)
    _fail_unless(by_hand != predecessor, f"tau={tau}: the row cannot tell step from the predecessor")
    response = cell_response(Closure.STEP, 1)
    derived = response.transmission({response.g[0]: sp.Rational(g)})[0, 0]
    _fail_unless(
        derived == sp.Rational(by_hand),
        f"tau={tau}: step T_1D = {derived}, the hand balance gives {by_hand} (the predecessor's {predecessor})",
    )
