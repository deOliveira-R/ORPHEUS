r"""The EXACT infinite-medium eigenpair of the float cross sections, in rational arithmetic.

What it computes
----------------
For an infinite homogeneous medium the loss and production operators are the
``G × G`` matrices

.. math::

   A = \operatorname{diag}(\Sigma_t) - \Sigma_{s0}^{T} - 2\Sigma_{2,0}^{T},
   \qquad
   F = \chi\,(\nu\Sigma_f)^{T},

and the eigenvalue question is :math:`A\varphi = F\varphi / k`. Every float64
is a dyadic rational, so the matrices built from the FLOAT inputs are exact
rational matrices, and this module solves the question for them exactly with
:class:`fractions.Fraction`: no rounding occurs anywhere. The answer is the
exact value production's floating-point answer approximates — "the exact
answer for the float inputs", which differs from the exact answer for the
DECIMAL cross sections a library tabulates (the latter is not what any float
computation can be held to).

The rank-one theorem
--------------------
:math:`F` has rank one, so :math:`K = A^{-1}F = (A^{-1}\chi)(\nu\Sigma_f)^{T}`
has rank one: its only non-zero eigenvalue is

.. math::

   k_\infty = \operatorname{tr} K = \langle \nu\Sigma_f,\, A^{-1}\chi \rangle,

with eigenvector :math:`u = A^{-1}\chi` (indeed :math:`Ku = u\,\langle
\nu\Sigma_f, u\rangle`). It is the dominant eigenvalue with no condition
on :math:`A` at all: a rank-one operator's other eigenvalues are all zero.
That the eigenvector lies in the positive cone does need one: :math:`u \ge 0`
when :math:`A^{-1} \ge 0`, which holds when :math:`A` is a non-singular
M-matrix. Its off-diagonal entries are :math:`\le 0` always; it is column
diagonally dominant when every group's removal exceeds its (n,2n) gain,
:math:`\Sigma_{c,g} + \Sigma_{f,g} > \Sigma_{2,g}` (the column sum is
:math:`\Sigma_t - \Sigma_s^{\rm out} - 2\Sigma_2^{\rm out}`), which a
medium whose fast group is dominated by (n,2n) (beryllium) can violate.
:meth:`ExactInfiniteMedium.certify` checks the cone on the result instead of
assuming it.

The module does NOT rest on the theorem for its correctness: every result
certifies itself against the DEFINING equations in exact arithmetic
(:meth:`ExactInfiniteMedium.certify`): :math:`A A^{-1} = I`,
:math:`F\varphi = k\,A\varphi` with zero residual, and
:math:`\operatorname{tr}(A^{-1}F) = k`. A wrong elimination, a dropped
transpose or a wrong gauge reddens there, before any comparison with
production.

The gauge is production's: :math:`\langle\nu\Sigma_f,\varphi\rangle = 100`
n/cm³/s, so :math:`\varphi = 100\,u / k_\infty`. The one-group condensed
cross sections are the flux-weighted means
:math:`\bar\sigma_x = \langle\Sigma_x,\varphi\rangle / \langle 1,\varphi\rangle`.

Dependencies: :mod:`fractions` for the arithmetic and numpy only to read the
float arrays (and scipy-sparse ``.toarray()`` duck-typed). No production
solver, operator or LAPACK routine is reached, which is what makes this a
structurally independent reference for
:func:`orpheus.homogeneous.solver.solve_homogeneous_infinite`
(``instrument-doctrine`` X4: it shares the float INPUTS with production —
the point — and nothing above them).
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import TYPE_CHECKING, Any, assert_never

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.numerics.mesh_free_function import RegionwiseConstant, values_without_position
from orpheus.numerics.observable import Eigenvalue, FluxIntegral, Linear, PointValue
from orpheus.numerics.question import Eigen
from orpheus.reference.certificate import Claim, Exact, ReferenceCertificate, Uncertifiable
from orpheus.reference.published import Current
from orpheus.reference.reading import Uncertified
from orpheus.reference.solution import ReferenceSolution
from orpheus.specification.specification import InfiniteMediumSpecification

if TYPE_CHECKING:
    import sympy

import numpy as np

Matrix = list[list[Fraction]]
Vector = list[Fraction]


def _as_exact_vector(values: Any) -> Vector:
    return [Fraction(float(v)) for v in np.asarray(values, dtype=float).ravel()]


def _as_exact_matrix(values: Any) -> Matrix:
    dense = values.toarray() if hasattr(values, "toarray") else values
    return [[Fraction(float(v)) for v in row] for row in np.asarray(dense, dtype=float)]


def _identity(n: int) -> Matrix:
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


def _matmul(a: Matrix, b: Matrix) -> Matrix:
    return [[sum((a[i][k] * b[k][j] for k in range(len(b))), Fraction(0)) for j in range(len(b[0]))] for i in range(len(a))]


def _matvec(a: Matrix, x: Vector) -> Vector:
    return [sum((a_ij * x_j for a_ij, x_j in zip(row, x)), Fraction(0)) for row in a]


def _dot(a: Vector, b: Vector) -> Fraction:
    return sum((x * y for x, y in zip(a, b)), Fraction(0))


def exact_inverse(a: Matrix) -> Matrix:
    """Gauss–Jordan elimination in exact arithmetic (any non-zero pivot is exact)."""
    n = len(a)
    work = [row[:] + ident for row, ident in zip(a, _identity(n))]
    for col in range(n):
        pivot = next((r for r in range(col, n) if work[r][col] != 0), None)
        if pivot is None:
            raise ZeroDivisionError("the loss matrix A is singular")
        work[col], work[pivot] = work[pivot], work[col]
        lead = work[col][col]
        work[col] = [v / lead for v in work[col]]
        for r in range(n):
            if r != col and work[r][col] != 0:
                factor = work[r][col]
                work[r] = [v - factor * p for v, p in zip(work[r], work[col])]
    return [row[n:] for row in work]


@dataclass(frozen=True)
class ExactInfiniteMedium:
    """The exact eigenpair of one float mixture, and the operators it came from.

    Attributes
    ----------
    loss : the exact :math:`A` of the float inputs (``G × G``).
    loss_inverse : :math:`A^{-1}`, exact.
    chi, nu_sigma_f, sigma_a : the float inputs as exact rationals.
    u : :math:`A^{-1}\\chi`, the un-gauged eigenvector.
    k_inf : :math:`\\langle\\nu\\Sigma_f, u\\rangle`.
    flux : :math:`100\\,u/k_\\infty` (gauge :math:`\\langle\\nu\\Sigma_f,\\varphi\\rangle = 100`).
    sig_prod, sig_abs : the flux-weighted one-group :math:`\\nu\\Sigma_f` and :math:`\\Sigma_a`.
    """

    loss: Matrix
    loss_inverse: Matrix
    chi: Vector
    nu_sigma_f: Vector
    sigma_a: Vector
    u: Vector
    k_inf: Fraction
    flux: Vector
    sig_prod: Fraction
    sig_abs: Fraction

    def certify(self) -> None:
        """Raise unless every result satisfies its defining equation EXACTLY."""
        g = len(self.chi)
        if _matmul(self.loss, self.loss_inverse) != _identity(g):
            raise ArithmeticError("A A^{-1} != I")
        production = [[c * f for f in self.nu_sigma_f] for c in self.chi]
        lhs = _matvec(production, self.flux)
        rhs = [self.k_inf * v for v in _matvec(self.loss, self.flux)]
        if lhs != rhs:
            raise ArithmeticError("F phi != k A phi: the eigen-equation residual is not zero")
        resolvent = _matmul(self.loss_inverse, production)
        if sum((resolvent[i][i] for i in range(g)), Fraction(0)) != self.k_inf:
            raise ArithmeticError("tr(A^{-1} F) != k: k is not the rank-one eigenvalue")
        if _dot(self.nu_sigma_f, self.flux) != 100:
            raise ArithmeticError("the gauge <nu Sigma_f, phi> = 100 does not hold")
        if self.k_inf <= 0 or any(v < 0 for v in self.u):
            raise ArithmeticError("the eigenvector is not in the positive cone")


def exact_pencil_eigenpair(
    *,
    loss: Any,
    chi: Any,
    nu_sigma_f: Any,
    sigma_a: Any,
) -> ExactInfiniteMedium:
    r"""The exact eigenpair of the pencil :math:`(A, \chi\,(\nu\Sigma_f)^T)` for a float ``G × G`` loss matrix.

    The pencil form: ``loss`` is taken as given (e.g. a materialised
    ``loss.as_matrix()``), so the answer is the exact eigenpair of THAT
    matrix and the two production factors — "the answer from the pair alone".
    """
    a = loss if isinstance(loss, list) else _as_exact_matrix(loss)
    loss_inverse = exact_inverse(a)
    chi_x, nsf, sa = _as_exact_vector(chi), _as_exact_vector(nu_sigma_f), _as_exact_vector(sigma_a)
    u = _matvec(loss_inverse, chi_x)
    k_inf = _dot(nsf, u)
    flux = [100 * v / k_inf for v in u]
    total = sum(flux, Fraction(0))
    result = ExactInfiniteMedium(
        loss=a,
        loss_inverse=loss_inverse,
        chi=chi_x,
        nu_sigma_f=nsf,
        sigma_a=sa,
        u=u,
        k_inf=k_inf,
        flux=flux,
        sig_prod=_dot(nsf, flux) / total,
        sig_abs=_dot(sa, flux) / total,
    )
    result.certify()
    return result


def exact_infinite_medium(
    *,
    sigma_t: Any,
    scattering_p0: Any,
    n2n_p0: Any,
    chi: Any,
    nu_sigma_f: Any,
    sigma_a: Any,
) -> ExactInfiniteMedium:
    r"""The exact eigenpair for float cross sections.

    ``scattering_p0`` and ``n2n_p0`` are the P0 transfer matrices in the
    mixture's FROM-row convention (``Σ_s0[g_from, g_to]``); the loss operator
    takes their transposes, assembled here EXACTLY (no rounding, unlike the
    float assembly production performs). ``sigma_a`` is the absorption cross
    section the caller's condensation uses, taken as a float input.
    """
    t = _as_exact_vector(sigma_t)
    s0 = _as_exact_matrix(scattering_p0)
    n0 = _as_exact_matrix(n2n_p0)
    g = len(t)
    loss = [
        [(t[i] if i == j else Fraction(0)) - s0[j][i] - 2 * n0[j][i] for j in range(g)]
        for i in range(g)
    ]
    return exact_pencil_eigenpair(loss=loss, chi=chi, nu_sigma_f=nu_sigma_f, sigma_a=sigma_a)


def exact_infinite_medium_of(mixture: Any) -> ExactInfiniteMedium:
    """Read a :class:`~orpheus.data.macro_xs.mixture.Mixture`'s float data (duck-typed)."""
    return exact_infinite_medium(
        sigma_t=mixture.SigT,
        scattering_p0=mixture.SigS[0],
        n2n_p0=mixture.Sig2[0],
        chi=mixture.chi,
        nu_sigma_f=mixture.SigP,
        sigma_a=mixture.absorption_xs,
    )


__all__ = [
    "ExactInfiniteMedium",
    "exact_infinite_medium",
    "exact_pencil_eigenpair",
    "exact_infinite_medium_of",
    "exact_inverse",
]


# ── the reference solution (#405 P2 step 6) ────────────────────────────────


def _rational(value: Fraction) -> sympy.Rational:
    import sympy

    return sympy.Rational(value.numerator, value.denominator)


@dataclass(frozen=True, eq=False)
class ExactInfiniteMediumDerivation:
    r"""The exact infinite medium as a reference derivation: every reading exact.

    The answer is :class:`ExactInfiniteMedium`'s exact rationals, certified
    against their defining equations at construction. The eigenvalue reads
    :math:`k_\infty` in the k chart (the chart of the question's parameter,
    the fission emission; other charts are #529's). A flux integral
    :math:`\langle w, \varphi\rangle = \sum_g w_g \varphi_g` is read under the
    infinite medium's measure, per unit volume (stated once, on
    :mod:`orpheus.numerics.observable`), in the production gauge
    :math:`\langle\nu\Sigma_f,\varphi\rangle = 100`. Its weight is the one
    region's group values: a regionwise-constant table read exactly (every
    double is a dyadic rational), or a symbolic weight read on the medium,
    which has no coordinate (:meth:`~orpheus.numerics.mesh_free_function.Symbolic.without`,
    the same independence admission decided). Every rational reading (the
    eigenvalue, a flux integral of a tabulated weight) is
    :class:`~orpheus.reference.certificate.Exact`; a symbolic weight whose
    value SymPy cannot certify to the working precision (cancellation) reads
    :class:`~orpheus.reference.reading.Uncertified`.
    """

    medium: ExactInfiniteMedium

    def evaluate(self, observable: Eigenvalue | Linear) -> Exact | Uncertified:
        """The observable's exact value, or its float value uncertified where the exact one cannot be certified."""
        import sympy

        value = self._exact_value(observable)
        try:
            return Exact(sympy.srepr(value), "ExactInfiniteMedium.certify")
        except Uncertifiable as refusal:
            return Uncertified(refusal.approximation)

    def _exact_value(self, observable: Eigenvalue | Linear) -> sympy.Expr:
        import sympy

        flux = [_rational(f) for f in self.medium.flux]
        match observable:
            case Eigenvalue():
                value = _rational(self.medium.k_inf)
            case FluxIntegral(weight=weight):
                constants = values_without_position(weight, len(flux))
                # a SymPy Float multiplies at its own precision: replace each by the exact binary value it holds
                exact_constants = [c.xreplace({f: sympy.Rational(f) for f in c.atoms(sympy.Float)}) for c in constants]
                value = sum((c * f for c, f in zip(exact_constants, flux, strict=True)), sympy.Integer(0))
            case PointValue():
                raise ValueError("the infinite medium has no position, so a point value has no reading")
            case _:
                assert_never(observable)
        return value


def exact_infinite_medium_reference(specification: InfiniteMediumSpecification) -> ReferenceSolution:
    r"""The exact infinite medium's :class:`~orpheus.reference.solution.ReferenceSolution`.

    It answers exactly one question: the fundamental k-eigenvalue at the
    physical point, :math:`\mathrm{Eigen}` along every fission emission with
    the default mode and no offset (any other mode or point is a different
    question, refused). Its certificate claims the eigenvalue and each group's
    flux (the indicator weight of that group), each
    :class:`~orpheus.reference.certificate.Exact`, with a target of
    :math:`10^{-12}` relative (the exact rationals' only error is their
    rounding to a double).
    """
    if not isinstance(specification, InfiniteMediumSpecification):
        raise TypeError(f"the exact infinite medium answers an InfiniteMediumSpecification, got a {type(specification).__name__}")
    k_question = Eigen(CellCoefficient.every(Channel.FISSION_EMISSION).resolve(specification.materials))
    if specification.question != k_question:
        raise ValueError(
            f"the exact infinite medium answers the fundamental k-eigenvalue question at the physical point, "
            f"got {specification.question!r}"
        )
    derivation = ExactInfiniteMediumDerivation(exact_infinite_medium_of(specification.mixture))
    groups = specification.n_groups
    observables: list[Eigenvalue | Linear] = [Eigenvalue()]
    observables += [FluxIntegral(RegionwiseConstant(np.eye(groups)[g][None, :])) for g in range(groups)]
    claims = {}
    for observable in observables:
        match derivation.evaluate(observable):
            case Exact() as established:
                claims[observable] = Claim(max(1e-12 * abs(established.enclosure().value), 1e-300), established)
            case Uncertified():
                raise ValueError(f"the exact infinite medium reads {observable!r} uncertified, though its k and fluxes are rational")
    return ReferenceSolution(specification, derivation, ReferenceCertificate(claims, (), (), Current()))
