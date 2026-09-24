r"""The algebra of record for the octant face transmission of step, diamond and LD.

This module is the **source of truth** for the face-transmission claims in
``docs/theory/methods/sn/cartesian_multid.rst`` (the label
``dd-face-transmission-spectrum`` and the section
``sn-boundary-gs-not-regular``).  If a transmission claim on that page cannot
be produced by this module, it is added here first.

One family: a spatial closure is a Petrov-Galerkin choice
===========================================================

Fix one ordinate and one Cartesian cell in ``d`` dimensions.  Measure lengths
in mean free paths (:math:`\Sigma_t = 1`), so that the cell and the ordinate
enter only through the per-axis **streaming coefficients**

.. math::

    g_a \;=\; \frac{|\mu_a|}{\Delta_a} \;=\; \frac{1}{\tau_a},
    \qquad G \;=\; \sum_a g_a ,

where :math:`\tau_a = \Sigma_t \Delta_a / |\mu_a|` is the optical thickness
the ordinate sees along axis :math:`a`.  These are the coefficients the
production kernel consumes (``s_axes`` of ``cell_kernel_batch``).  On the
reference cell :math:`\xi \in [-1, 1]^d` the streaming term per unit volume
is :math:`2 g_a \partial_{\xi_a}`, and a face integral per unit volume is
:math:`g_a` times the face average.

A closure is four choices (:class:`PetrovGalerkinScheme`):

* a **trial space** :math:`V`, spanned by tensor Legendre monomials
  :math:`\phi_o(\xi) = \prod_a P_{o_a}(\xi_a)`, :math:`o \in \{0, 1\}^d`;
* a **test space** :math:`W \subseteq V`;
* a **face space** :math:`\Phi`, the moments a face carries between cells;
* how the inflow enters: **weakly**, through the upwind face term of the
  weak form, or **strongly**, as a constraint on the trial function's
  inflow trace.

=========  ===================================  ===============  =====================  ======
closure    trial :math:`V`                      test :math:`W`   face :math:`\Phi`      inflow
=========  ===================================  ===============  =====================  ======
step       constants (DG0)                      constants        constants              weak
LD         bilinear :math:`Q_1` (DG1)           :math:`Q_1`      :math:`Q_1` in d - 1   weak
diamond    :math:`1` and each :math:`\xi_a`     constants        constants              strong
=========  ===================================  ===============  =====================  ======

The weak form, for every test function :math:`v \in W`, is

.. math::

    \sum_a g_a\Big(\langle v\,\psi\rangle_{a+} - \langle v\,\psi^{-}_a\rangle_{a-}
    - 2\langle \psi\,\partial_a v\rangle\Big) + \langle v\,\psi\rangle
    \;=\; \langle v\, q\rangle ,

with :math:`\langle\cdot\rangle` the cell average,
:math:`\langle\cdot\rangle_{a\pm}` the average over the face
:math:`\xi_a = \pm 1`, and :math:`\psi^-_a` the upwind data: the inflow
:math:`\psi^{\rm in}_a` for a weak closure, the trial function's own trace
for a strong one, whose inflow instead enters through the constraints
:math:`\Pi_\Phi(\psi|_{\xi_a = -1}) = \psi^{\rm in}_a`.  Every matrix below
is an integral of Legendre monomials over the reference cell, computed here;
none is typed.

The cell response, a state-space realization
============================================

Eliminating the strongly constrained trial coordinates (those of
:math:`V \setminus W`) leaves the minimal realization of the cell, its
**cell response** (:class:`CellResponse`):

.. math::

    A\,c \;=\; E\,\psi^{\rm in} + S\,q, \qquad
    \psi^{\rm out} \;=\; R\,c + D\,\psi^{\rm in},

with :math:`c` the free trial coordinates (the cell state), :math:`A` the
cell operator, :math:`E` the inflow scatter (for a weak closure, the transpose of the upstream
trace, weighted by the face mass), :math:`S` the source mass, :math:`R` the
outflow trace, and :math:`D` the **feedthrough**: the part of the outflow the
inflow reaches without passing through the cell.  It is the transfer function
of a linear system, and its four blocks are the cell's response matrix:

* **transmission** :math:`T = R A^{-1} E + D`, inflow to outflow;
* **escape** :math:`R A^{-1} S`, source to outflow;
* **inflow to cell** :math:`A^{-1} E`; **source to cell** :math:`A^{-1} S`.

These are the objects of the response-matrix and interface-current methods
(transmission and escape probabilities), here for one cell and one ordinate.

What is proved
==============

1. **The realizations** (:func:`derive_realization`): LD's weak form equals
   the Kronecker-factor UBLD of
   :mod:`orpheus.derivations.discrete.sn.ld_ubld`, which was built by a
   different route; step's is the balance with :math:`\psi^{\rm out} =
   \psi_c`; diamond's reproduces :math:`\psi^{\rm out}_a = 2\psi_c -
   \psi^{\rm in}_a`.  The feedthrough is :math:`D = 0` for a weak closure
   and :math:`D = -I` for diamond: strong imposition is what creates it.
2. **The feedthrough eigenvalue** (:func:`transmission_spectrum`).  When
   :math:`D = fI`, Sylvester's determinant identity gives
   :math:`\det(\lambda I_n - T) = (\lambda - f)^{n-m}\det((\lambda - f)A -
   ER)/\det A`, with :math:`n` face moments and :math:`m` cell moments, so
   :math:`f` is an eigenvalue of multiplicity at least :math:`n - m`.  For
   diamond that is :math:`-1`, :math:`d - 1` times: the undamped face
   sawtooth is the feedthrough and nothing else.
3. **Every other eigenvalue lies strictly inside the unit disk**, for every
   :math:`g > 0`, all three closures, :math:`d \in \{1, 2, 3\}` (the Jury
   conditions on each factor, each with a positivity certificate).
4. **The conserved mode** (:func:`derive_conserved_mode`):
   :math:`\lambda - T_{1\text{-}D}(G)` is a factor of every closure's
   characteristic polynomial, where :math:`T_{1\text{-}D}` is the closure's
   own 1-D transmission; along :math:`g_a = t \to \infty` its root tends to 1.
5. **LD's spectrum is generated by its low dimensions**
   (:func:`derive_ld_factorisation`): at :math:`d = 3` the characteristic
   polynomial is :math:`\lambda^5\,\ell(G)\prod_a q(g_a, G - g_a)`.
6. **The stability function** (:func:`derive_stability_function`): in 1-D
   the transmission is a Padé approximant of :math:`e^{-\tau}` (step
   :math:`[0/1]`, diamond :math:`[1/1]`, LD :math:`[1/2]`); it is A-stable
   (bounded by 1 on the right half-plane, by the E-polynomial criterion) and
   strictly below 1 in modulus for every real :math:`\tau > 0`; and its value at
   :math:`\tau \to \infty` is the feedthrough :math:`f`: 0 for the weak
   closures (L-stable), :math:`-1` for diamond (A-stable only).  In every
   dimension :math:`T \to fI` as :math:`g \to 0`.  Diamond's undamped face
   mode is the trapezoidal (Crank-Nicolson) rule's undamped stiff mode, seen
   across a cell.
7. **Conservation** (:func:`derive_particle_balance`): for every inflow and
   source, outflow minus inflow plus absorption equals the source, read from
   the four blocks.  **Flat-flux preservation**
   (:func:`derive_flat_flux_preservation`): a uniform inflow equal to a
   uniform source passes through unchanged, in the outflow and in the cell.
8. **The page's closed forms** (:func:`derive_page_closed_form`):
   :math:`(2/D)\mathbf 1 w^{\mathsf T} - I` for diamond,
   :math:`(1/D')\mathbf 1 w'^{\mathsf T}` for step.

Branch 1 only.  The production counterparts are each scheme's
``cell_kernel_batch`` and ``DiscretizationSchemeBase.face_transmission_spectrum``
(``orpheus/transport/spatial/scheme.py``); the L1 cross-check lives in the
tests.

History: first written 2026-08-09 as ``derivations/sn_dd_face_transmission.py``
for #341.  That version built step's matrix with diamond's weights
(:math:`w_a = 2|\mu_a|A_a`), so its "step" was :math:`2/(2+\tau)` at
:math:`d = 1` instead of :math:`1/(1+\tau)`; its qualitative claim survived,
its eigenvalue formula did not.  Here every closure is derived from its weak
form, never typed.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from functools import cache, cached_property
from itertools import combinations, product

import sympy as sp
from sympy.polys.matrices import DomainMatrix

LAMBDA = sp.Symbol("lambda")
r"""The eigenvalue variable of every characteristic polynomial here."""

MultiIndex = tuple[int, ...]
Point = Mapping[sp.Symbol, sp.Expr]


class DerivationFailed(AssertionError):
    """A step of a proof in this module did not hold.

    Raised explicitly, so ``python -O`` cannot strip it.
    """


def _require(holds: bool, claim: str) -> None:
    if not holds:
        raise DerivationFailed(claim)


def streaming_coefficients(d: int) -> tuple[sp.Symbol, ...]:
    r"""The positive symbols :math:`g_1, \dots, g_d`, with :math:`g_a = |\mu_a|/\Delta_a`."""
    return sp.symbols(f"g1:{d + 1}", positive=True)


def reference_coordinates(d: int) -> tuple[sp.Symbol, ...]:
    r"""The reference-cell coordinates :math:`\xi_1, \dots, \xi_d \in [-1, 1]`."""
    return sp.symbols(f"xi1:{d + 1}", real=True)


# ── the reference cell: Legendre monomials, averages, face projections ──


def _monomial(o: MultiIndex, xi: tuple[sp.Symbol, ...]) -> sp.Expr:
    r""":math:`\phi_o(\xi) = \prod_a P_{o_a}(\xi_a)`, the tensor Legendre monomial."""
    return sp.Mul(*(sp.legendre(k, x) for k, x in zip(o, xi, strict=True)))


def _average(expr: sp.Expr, xi: tuple[sp.Symbol, ...]) -> sp.Expr:
    """The average over the reference cube spanned by ``xi`` (each axis's measure is dx/2)."""
    for x in xi:
        expr = sp.integrate(expr, (x, -1, 1)) / 2
    return sp.expand(expr)


def _transverse(xi: tuple[sp.Symbol, ...], axis: int) -> tuple[sp.Symbol, ...]:
    return tuple(x for a, x in enumerate(xi) if a != axis)


def _face_average(expr: sp.Expr, xi: tuple[sp.Symbol, ...], axis: int, node: int) -> sp.Expr:
    r"""The average over the face :math:`\xi_{\rm axis} = {\rm node}`."""
    return _average(expr.subs(xi[axis], node), _transverse(xi, axis))


def _face_moments(
    expr: sp.Expr, xi: tuple[sp.Symbol, ...], axis: int, node: int, face: tuple[MultiIndex, ...],
) -> list[sp.Expr]:
    r"""The coefficients of the trace on face ``(axis, node)`` in the face space: an :math:`L^2` projection."""
    transverse = _transverse(xi, axis)
    return [
        _face_average(expr * _monomial(f, transverse), xi, axis, node)
        / _average(_monomial(f, transverse) ** 2, transverse)
        for f in face
    ]


# ── the Petrov-Galerkin family ──


class Imposition(Enum):
    """How the inflow data enters the cell."""

    WEAK = "weak"
    r"""Through the upwind face term, :math:`\langle v\,\psi^{\rm in}\rangle_{a-}`."""

    STRONG = "strong"
    r"""As a constraint on the trial function, :math:`\Pi_\Phi(\psi|_{\xi_a=-1}) = \psi^{\rm in}_a`."""


@dataclass(frozen=True)
class PetrovGalerkinScheme:
    """A spatial closure on the reference cell: trial, test and face spaces, and the inflow imposition.

    The spaces are tuples of Legendre multi-indices; the face space is
    indexed over the ``d - 1`` transverse axes, in increasing axis order.
    Strong imposition fixes the trial coordinates outside the test space,
    one per inflow constraint; weak imposition has none to fix.
    """

    d: int
    trial: tuple[MultiIndex, ...]
    test: tuple[MultiIndex, ...]
    face: tuple[MultiIndex, ...]
    inflow: Imposition

    def __post_init__(self) -> None:
        for name, space, arity in (
            ("trial", self.trial, self.d), ("test", self.test, self.d), ("face", self.face, self.d - 1),
        ):
            _require(
                len(set(space)) == len(space) and all(len(o) == arity for o in space),
                f"the {name} space is distinct multi-indices of length {arity}",
            )
        _require(self.face[0] == (0,) * (self.d - 1), "the face space lists the face average first")
        _require(set(self.test) <= set(self.trial), "the test space lies in the trial space")
        constrained = len(self.trial) - len(self.test)
        expected = self.d * len(self.face) if self.inflow is Imposition.STRONG else 0
        _require(
            constrained == expected,
            f"{self.inflow.value} imposition: {constrained} trial coordinates beyond the test space, "
            f"{expected} inflow constraints",
        )


def _constants(d: int) -> tuple[MultiIndex, ...]:
    return ((0,) * d,)


def _bilinear(d: int) -> tuple[MultiIndex, ...]:
    return tuple(product((0, 1), repeat=d))


def _constant_and_axis_slopes(d: int) -> tuple[MultiIndex, ...]:
    """:math:`P_1`: the constant and one linear term per axis, no cross terms."""
    return _constants(d) + tuple(tuple(int(b == a) for b in range(d)) for a in range(d))


class Closure(Enum):
    """The three spatial closures compared here, each a Petrov-Galerkin choice (:data:`SCHEMES`)."""

    STEP = "step"
    DIAMOND = "diamond"
    LINEAR_DISCONTINUOUS = "linear_discontinuous"

    def scheme(self, d: int) -> PetrovGalerkinScheme:
        """This closure's trial, test and face spaces and inflow imposition in ``d`` dimensions."""
        trial, test, face, inflow = SCHEMES[self]
        return PetrovGalerkinScheme(d, trial(d), test(d), face(d - 1), inflow)


SCHEMES = {
    Closure.STEP: (_constants, _constants, _constants, Imposition.WEAK),
    Closure.LINEAR_DISCONTINUOUS: (_bilinear, _bilinear, _bilinear, Imposition.WEAK),
    Closure.DIAMOND: (_constant_and_axis_slopes, _constants, _constants, Imposition.STRONG),
}
r"""Each closure as (trial, test, face space, inflow imposition): step is
:math:`(Q_0, Q_0, Q_0, \text{weak})`, LD :math:`(Q_1, Q_1, Q_1, \text{weak})`,
diamond :math:`(P_1, Q_0, Q_0, \text{strong})`.  A closure with no row raises
``KeyError``."""


# ── the cell response ──


def _over_field(M: sp.Matrix, field) -> DomainMatrix:
    return DomainMatrix.from_Matrix(M).convert_to(field)


@dataclass(frozen=True)
class CellResponse:
    r"""The minimal realization of one closure on one cell, for one octant.

    ``A c = E ψ_in + S q`` and ``ψ_out = R c + D ψ_in``; see the module
    docstring.  Lengths are in mean free paths (:math:`\Sigma_t = 1`); the
    source ``q`` is given by its coefficients in the test space; face
    moments are axis-major, the average first on each face.  Every block
    accepts a point ``at`` (values of :math:`g`) and is then exact over
    :math:`\mathbb{Q}`; without one it is symbolic over :math:`\mathbb{Q}(g)`.
    """

    closure: Closure
    g: tuple[sp.Symbol, ...]
    cell_operator: sp.Matrix
    inflow_scatter: sp.Matrix
    source_mass: sp.Matrix
    outflow_trace: sp.Matrix
    feedthrough: sp.Matrix

    @property
    def d(self) -> int:
        return len(self.g)

    @property
    def n_face_moments(self) -> int:
        """Face moments over the octant's ``d`` outflow faces (equally, its inflow faces)."""
        return self.outflow_trace.rows

    @property
    def n_cell_moments(self) -> int:
        return self.cell_operator.rows

    @cached_property
    def feedthrough_scalar(self) -> sp.Rational:
        r"""The :math:`f` of :math:`D = fI`.

        SCOPE-BOUNDARY: the Sylvester route of :attr:`characteristic_factors`
        needs :math:`D = fI`, which all three closures here satisfy.  A
        closure with a general :math:`D` would need
        :math:`\det(\lambda I - T) = \det(\lambda I - D)\,
        \det(A - E(\lambda I - D)^{-1}R)/\det A`, not built; it is refused
        with ``NotImplementedError``, which is a scope edge, not a failed proof.
        """
        f = self.feedthrough[0, 0]
        if not (f.is_Rational and self.feedthrough == f * sp.eye(self.n_face_moments)):
            raise NotImplementedError(
                f"{self.closure.value}: the feedthrough is not a rational multiple of the identity; "
                "the general form det(lambda I - D) det(A - E (lambda I - D)^-1 R) / det A is not built"
            )
        return sp.Rational(f)

    def _at(self, M: sp.Matrix, at: Point | None) -> sp.Matrix:
        return M if at is None else M.subs(at)

    def _solve(self, block: sp.Matrix, at: Point | None) -> sp.Matrix:
        r""":math:`A^{-1}` applied to ``block``, exactly.

        Solved in SymPy's polynomial domains rather than on expression
        trees: LD's :math:`d = 2` transmission takes 0.06 s this way and did
        not finish in 30 s through ``LUsolve`` plus ``cancel`` `[M]` 2026-09-23.
        """
        field = sp.QQ if at is not None else sp.QQ.frac_field(*self.g)
        A = _over_field(self._at(self.cell_operator, at), field)
        return (A.inv() * _over_field(self._at(block, at), field)).to_Matrix()

    def inflow_to_cell(self, at: Point | None = None) -> sp.Matrix:
        r""":math:`A^{-1} E`: the cell state an inflow produces."""
        return self._solve(self.inflow_scatter, at)

    def source_to_cell(self, at: Point | None = None) -> sp.Matrix:
        r""":math:`A^{-1} S`: the cell state a source produces."""
        return self._solve(self.source_mass, at)

    def transmission(self, at: Point | None = None) -> sp.Matrix:
        r""":math:`T = R A^{-1} E + D`: inflow to outflow."""
        return self._at(self.outflow_trace, at) * self.inflow_to_cell(at) + self.feedthrough

    def escape(self, at: Point | None = None) -> sp.Matrix:
        r""":math:`R A^{-1} S`: source to outflow."""
        return self._at(self.outflow_trace, at) * self.source_to_cell(at)

    @cached_property
    def characteristic_factors(self) -> tuple[tuple[sp.Poly, int], ...]:
        r"""The irreducible factors of :math:`\det(\lambda I - T)`, with multiplicity.

        By Sylvester's identity :math:`\det(\lambda I_n - T) = \mu^{n-m}
        \det(\mu A - E R) / \det A` with :math:`\mu = \lambda - f` (it needs
        :math:`D = fI` and :math:`\det A \ne 0`, both asserted), so only the
        :math:`m \times m` polynomial determinant is expanded and factored,
        and the power of :math:`\lambda - f` is corrected by :math:`n - m`.
        ``det A`` is free of :math:`\lambda` and drops out.  Each factor is a
        primitive ``Poly`` in :math:`(\lambda, g)` over :math:`\mathbb{Z}`,
        oriented by :func:`_canonical`.  At :math:`d = 3` LD's direct
        :math:`12 \times 12` rational determinant did not finish in ten
        minutes `[M]` 2026-09-23; this route takes seconds.
        """
        f = self.feedthrough_scalar
        _require(
            _is_positive_on_orthant(self.cell_operator.det(method="berkowitz"), self.g),
            f"{self.closure.value}, d={self.d}: det A > 0 for every g > 0",
        )
        shifted = LAMBDA - f
        cell_det = (shifted * self.cell_operator - self.inflow_scatter * self.outflow_trace).det(
            method="berkowitz"
        )
        _, raw = sp.Poly(cell_det, LAMBDA, *self.g).factor_list()
        counts: dict[sp.Poly, int] = {}
        for poly, mult in raw:
            if poly.degree(LAMBDA) == 0:
                continue
            key = _canonical(poly.as_expr(), self.g)
            counts[key] = counts.get(key, 0) + mult
        linear = _canonical(shifted, self.g)
        counts[linear] = counts.get(linear, 0) + self.n_face_moments - self.n_cell_moments
        _require(all(mult >= 0 for mult in counts.values()), "Sylvester: the shifted power is non-negative")
        factors = tuple((fac, mult) for fac, mult in counts.items() if mult > 0)
        degree = sum(fac.degree(LAMBDA) * mult for fac, mult in factors)
        _require(degree == self.n_face_moments, f"degree {degree} == {self.n_face_moments} face moments")
        return factors


def _canonical(expr: sp.Expr, g: tuple[sp.Symbol, ...]) -> sp.Poly:
    """The primitive integer ``Poly`` of ``expr``'s numerator in (lambda, g), leading coefficient positive at g = 1.

    One canonical form, so factors compare by ``==`` and key a ``dict``.
    """
    numerator = sp.fraction(sp.together(expr))[0]
    _, primitive = sp.Poly(numerator, LAMBDA, *g).primitive()
    lead = sp.Poly(primitive.as_expr(), LAMBDA).LC().subs({ga: 1 for ga in g})
    return primitive if lead > 0 else -primitive


def _weak_form_response(closure: Closure, d: int) -> CellResponse:
    """Assemble the weak form of ``closure.scheme(d)`` and eliminate the constrained coordinates."""
    scheme = closure.scheme(d)
    g = streaming_coefficients(d)
    xi = reference_coordinates(d)
    n_face = len(scheme.face)

    u = sp.symbols(f"u0:{len(scheme.trial)}")
    psi = sum(ui * _monomial(o, xi) for ui, o in zip(u, scheme.trial))
    p = sp.symbols(f"p0:{d * n_face}")  # inflow face moments, axis-major
    inflow = [
        sum(p[a * n_face + k] * _monomial(f, _transverse(xi, a)) for k, f in enumerate(scheme.face))
        for a in range(d)
    ]
    q = sp.symbols(f"q0:{len(scheme.test)}")
    source = sum(qi * _monomial(o, xi) for qi, o in zip(q, scheme.test))

    upwind = inflow if scheme.inflow is Imposition.WEAK else [psi] * d
    residuals = []
    for o in scheme.test:
        v = _monomial(o, xi)
        streaming = sum(
            g[a] * (
                _face_average(v * psi, xi, a, 1)
                - _face_average(v * upwind[a], xi, a, -1)
                - 2 * _average(psi * sp.diff(v, xi[a]), xi)
            )
            for a in range(d)
        )
        residuals.append(sp.expand(streaming + _average(v * psi, xi) - _average(v * source, xi)))
    outflow = [moment for a in range(d) for moment in _face_moments(psi, xi, a, 1, scheme.face)]
    constraints = (
        [moment for a in range(d) for moment in _face_moments(psi, xi, a, -1, scheme.face)]
        if scheme.inflow is Imposition.STRONG else []
    )
    free = [ui for ui, o in zip(u, scheme.trial) if o in scheme.test]
    constrained = [ui for ui, o in zip(u, scheme.trial) if o not in scheme.test]
    return _static_condensation(closure, g, residuals, outflow, constraints, free, constrained, p, q)


def _linear_blocks(rows: list[sp.Expr], groups: list[list[sp.Symbol]], where: str) -> list[sp.Matrix]:
    """The Jacobian of ``rows`` with respect to each variable group, after proving ``rows`` linear in them."""
    blocks = [sp.Matrix(rows).jacobian(group) if group else sp.zeros(len(rows), 0) for group in groups]
    rebuilt = sum((blk * sp.Matrix(group) for blk, group in zip(blocks, groups) if group), sp.zeros(len(rows), 1))
    _require(
        (sp.Matrix(rows) - rebuilt).applyfunc(sp.expand) == sp.zeros(len(rows), 1),
        f"{where} is linear and homogeneous in the state, the inflow and the source",
    )
    return blocks


def _static_condensation(
    closure: Closure,
    g: tuple[sp.Symbol, ...],
    residuals: list[sp.Expr],
    outflow: list[sp.Expr],
    constraints: list[sp.Expr],
    free: list[sp.Symbol],
    constrained: list[sp.Symbol],
    p: tuple[sp.Symbol, ...],
    q: tuple[sp.Symbol, ...],
) -> CellResponse:
    r"""Eliminate the strongly constrained coordinates: the minimal realization, by a Schur complement.

    With the tested residuals :math:`K_f u_f + K_c u_c + P\psi^{\rm in} + Qq`,
    the outflow :math:`R_f u_f + R_c u_c + O\psi^{\rm in}`, and the inflow
    constraints :math:`C_f u_f + C_c u_c = \psi^{\rm in}` (none for a weak
    closure), the constrained coordinates are
    :math:`u_c = C_c^{-1}(\psi^{\rm in} - C_f u_f)`, so

    .. math::

        A = K_f - K_c C_c^{-1} C_f,\quad E = -(P + K_c C_c^{-1}),\quad S = -Q,

        R = R_f - R_c C_c^{-1} C_f,\quad D = O + R_c C_c^{-1}.

    :math:`C_c` must be square and invertible: one constraint per
    constrained coordinate, each actually imposed (a constraint that misses
    every constrained coordinate makes :math:`C_c` singular and is refused).
    The feedthrough :math:`D = R_c C_c^{-1}` of a strong closure is this
    identity: strong imposition is what creates it.
    """
    K_f, K_c, P, Q = _linear_blocks(residuals, [free, constrained, list(p), list(q)], f"{closure.value}: the residual")
    R_f, R_c, O, _no_source = _linear_blocks(outflow, [free, constrained, list(p), list(q)], f"{closure.value}: the outflow")
    _require(_no_source == sp.zeros(*_no_source.shape), f"{closure.value}: the source does not reach the outflow directly")
    if constrained:
        C_f, C_c, C_p, _ = _linear_blocks(
            [c for c in constraints], [free, constrained, list(p), list(q)], f"{closure.value}: the constraints",
        )
        _require(C_p == sp.zeros(*C_p.shape), f"{closure.value}: the constraint rows are the trace moments alone")
        _require(
            C_c.is_square and C_c.det() != 0 and not C_c.has(*g),
            f"{closure.value}: the inflow constraints fix every constrained coordinate (C_c square, invertible)",
        )
        C_c_inv = C_c.inv()
        A = K_f - K_c * C_c_inv * C_f
        E = -(P + K_c * C_c_inv)
        R = R_f - R_c * C_c_inv * C_f
        D = O + R_c * C_c_inv
    else:
        _require(not constraints, f"{closure.value}: a weak closure imposes no constraint")
        A, E, R, D = K_f, -P, R_f, O
    return CellResponse(
        closure, g, A.applyfunc(sp.expand), E.applyfunc(sp.expand), (-Q).applyfunc(sp.expand),
        R.applyfunc(sp.expand), D.applyfunc(sp.expand),
    )


@cache
def cell_response(closure: Closure, d: int) -> CellResponse:
    """The cell response of ``closure`` in ``d`` dimensions, from its weak form (cached)."""
    return _weak_form_response(closure, d)


# ── claim 1: the realizations ──


def derive_realization(closure: Closure, d: int) -> CellResponse:
    r"""The weak form reproduces each closure as it is known, by an independent route.

    * **LD**: :math:`A`, :math:`E`, :math:`S` and :math:`R` equal
      :mod:`~orpheus.derivations.discrete.sn.ld_ubld`'s Kronecker-factor
      UBLD (``assemble_ubld`` at unit widths, :math:`\mu_a \to g_a`,
      :math:`\Sigma_t = 1`, :math:`\theta = 1/3`; its per-axis inflow lift
      and outflow trace), and :math:`D = 0`.
    * **Step**: :math:`(1 + G)\psi_c = \sum_a g_a\psi^{\rm in}_a + q`,
      :math:`\psi^{\rm out}_a = \psi_c`, :math:`D = 0`.
    * **Diamond**: :math:`(1 + 2G)\psi_c = \sum_a 2g_a\psi^{\rm in}_a + q`,
      :math:`\psi^{\rm out}_a = 2\psi_c - \psi^{\rm in}_a`, :math:`D = -I`.

    Returns the response.
    """
    response = cell_response(closure, d)
    g = response.g
    if closure is Closure.LINEAR_DISCONTINUOUS:
        from orpheus.derivations.discrete.sn.ld_ubld import (
            assemble_ubld,
            inflow_scatter_axis,
            outflow_trace_axis,
        )

        ones = [sp.Integer(1)] * d
        theta = sp.Rational(1, 3)
        kron = assemble_ubld(ones, list(g), sp.Integer(1), theta)
        known = {
            "A": kron["A"],
            "E": sp.Matrix.hstack(*(inflow_scatter_axis(ones, list(g), a, theta) for a in range(d))),
            "S": kron["M"],
            "R": sp.Matrix.vstack(*(outflow_trace_axis(d, a) for a in range(d))),
            "D": sp.zeros(response.n_face_moments),
        }
    else:
        weight = 1 if closure is Closure.STEP else 2
        f = 0 if closure is Closure.STEP else -1
        known = {
            "A": sp.Matrix([[1 + weight * sum(g)]]),
            "E": sp.Matrix([[weight * ga for ga in g]]),
            "S": sp.Matrix([[1]]),
            "R": weight * sp.ones(d, 1),
            "D": f * sp.eye(d),
        }
    derived = {
        "A": response.cell_operator, "E": response.inflow_scatter, "S": response.source_mass,
        "R": response.outflow_trace, "D": response.feedthrough,
    }
    for name, block in derived.items():
        _require(
            (block - known[name]).applyfunc(sp.expand) == sp.zeros(*block.shape),
            f"{closure.value}, d={d}: {name} is the known realization's",
        )
    return response


# ── positivity on the open orthant, the tool every Jury condition needs ──


def _is_positive_on_orthant(p: sp.Expr, variables: tuple[sp.Symbol, ...]) -> bool:
    r"""A certificate that :math:`p > 0` for every positive argument, or ``False``.

    Two certificates, tried in order:

    * :math:`p \ne 0` and every coefficient is non-negative (each monomial
      is positive on the open orthant);
    * otherwise, the negative coefficients all sit in the homogeneous
      quadratic part :math:`x^{\mathsf T} Q x`, every other coefficient is
      non-negative, and either :math:`Q` is positive definite (every leading
      principal minor positive) or :math:`Q` is positive semidefinite (every
      principal minor non-negative) while the other terms are not all zero.
      The semidefinite case is the one LD needs at :math:`d = 3`, where a
      Jury condition depends on :math:`(g_a, G - g_a)` only.

    ``False`` means *no certificate found*, never *not positive*.
    """
    poly = sp.Poly(sp.expand(p), *variables)
    if poly.is_zero:
        return False
    if all(c >= 0 for c in poly.coeffs()):
        return True
    rest = [c for mon, c in poly.terms() if sum(mon) != 2]
    if not all(c >= 0 for c in rest):
        return False
    k = len(variables)
    Q = sp.zeros(k, k)
    for mon, c in poly.terms():
        if sum(mon) != 2:
            continue
        i, j = [axis for axis, e in enumerate(mon) for _ in range(e)]
        if i == j:
            Q[i, i] = c
        else:
            Q[i, j] = Q[j, i] = c / 2
    if all(Q[:r, :r].det() > 0 for r in range(1, k + 1)):
        return True
    principal = (
        Q.extract(list(rows), list(rows)).det()
        for r in range(1, k + 1)
        for rows in combinations(range(k), r)
    )
    return any(c > 0 for c in rest) and all(minor >= 0 for minor in principal)


def _jury_strictly_inside(factor: sp.Poly, variables: tuple[sp.Symbol, ...]) -> bool:
    r"""Both roots (or the root) of a degree-1 or -2 factor lie strictly inside the unit disk.

    Degree 1, :math:`a\lambda + b`: :math:`a > 0`, :math:`a + b > 0`,
    :math:`a - b > 0`.  Degree 2, :math:`a\lambda^2 + b\lambda + c`:
    :math:`a > 0`, :math:`a - c > 0`, :math:`a + c > 0`, :math:`Q(1) > 0`,
    :math:`Q(-1) > 0` (Jury 1964; necessary and sufficient for real
    coefficients).  Every condition is certified for every positive
    argument.  Higher degrees are refused, not guessed.
    """
    coeffs = sp.Poly(factor.as_expr(), LAMBDA).all_coeffs()
    if len(coeffs) == 2:
        a, b = coeffs
        conditions = (a, a + b, a - b)
    elif len(coeffs) == 3:
        a, b, c = coeffs
        conditions = (a, a - c, a + c, a + b + c, a - b + c)
    else:
        raise NotImplementedError(f"Jury conditions are written for degree <= 2, got {len(coeffs) - 1}")
    return all(_is_positive_on_orthant(condition, variables) for condition in conditions)


def _is_a_stable(numerator: sp.Expr, denominator: sp.Expr, variable: sp.Symbol) -> bool:
    r"""A certificate that :math:`R = P/Q` satisfies :math:`|R| \le 1` on :math:`\operatorname{Re} z \ge 0`, or ``False``.

    The E-polynomial criterion (Hairer & Wanner, Vol. II, §IV.3, reasoned
    citation): every pole lies in :math:`\operatorname{Re} z < 0` (:math:`Q`
    Hurwitz; for degree :math:`\le 2`, all coefficients non-zero and of one
    sign) and :math:`E(y) = |Q(iy)|^2 - |P(iy)|^2 \ge 0` for real :math:`y`
    (a polynomial in :math:`y^2` with non-negative coefficients); the maximum
    modulus principle then bounds :math:`|R|` by 1 on the right half-plane.
    ``False`` means *no certificate found*; a denominator of degree above 2
    is refused, not guessed.
    """
    z, y = sp.Symbol("z"), sp.Symbol("y", real=True)
    P, Q = numerator.subs(variable, z), denominator.subs(variable, z)
    q_coeffs = sp.Poly(Q, z).all_coeffs()
    if len(q_coeffs) > 3:
        raise NotImplementedError(f"the Hurwitz test is written for degree <= 2, got {len(q_coeffs) - 1}")
    if not (all(c > 0 for c in q_coeffs) or all(c < 0 for c in q_coeffs)):
        return False
    E = sp.expand(sp.Abs(Q.subs(z, sp.I * y)) ** 2 - sp.Abs(P.subs(z, sp.I * y)) ** 2)
    if sp.expand(sp.im(E)) != 0:
        return False
    return all(e % 2 == 0 and c >= 0 for (e,), c in sp.Poly(sp.expand(sp.re(E)), y).terms())


# ── claims 2 and 3: the spectrum, by factor ──


@dataclass(frozen=True)
class TransmissionSpectrum:
    r"""The proven spectrum of one closure's octant transmission.

    ``factors`` are :attr:`CellResponse.characteristic_factors`.
    ``undamped`` lists the eigenvalues of modulus exactly 1 for every
    :math:`g > 0`, with multiplicity; every other root is proved strictly
    inside the unit disk.
    """

    closure: Closure
    d: int
    factors: tuple[tuple[sp.Poly, int], ...]
    feedthrough_multiplicity: int
    undamped: tuple[tuple[sp.Rational, int], ...]

    @property
    def damps_every_face_mode(self) -> bool:
        return not self.undamped


@cache
def transmission_spectrum(closure: Closure, d: int) -> TransmissionSpectrum:
    """Factor the characteristic polynomial and certify where every root lies."""
    response = cell_response(closure, d)
    f = response.feedthrough_scalar
    factors = response.characteristic_factors
    feedthrough = _canonical(LAMBDA - f, response.g)
    feed_mult = dict(factors).get(feedthrough, 0)
    _require(
        feed_mult >= response.n_face_moments - response.n_cell_moments,
        f"{closure.value}, d={d}: Sylvester forces the feedthrough {f} at least "
        f"{response.n_face_moments - response.n_cell_moments} times; found {feed_mult}",
    )
    for fac, _mult in factors:
        if fac == feedthrough:
            continue
        _require(
            _jury_strictly_inside(fac, response.g),
            f"{closure.value}, d={d}: the factor {fac.as_expr()} has a root on or outside the unit disk",
        )
    undamped = ((f, feed_mult),) if abs(f) == 1 and feed_mult > 0 else ()
    return TransmissionSpectrum(closure, d, factors, feed_mult, undamped)


# ── claim 4: the conserved mode is the 1-D transmission at the total streaming ──


def _one_dimensional_transmission(closure: Closure, g: sp.Expr) -> sp.Expr:
    (g1,) = streaming_coefficients(1)
    return sp.cancel(cell_response(closure, 1).transmission()[0, 0].subs(g1, g))


def derive_conserved_mode(closure: Closure, d: int) -> sp.Expr:
    r""":math:`\lambda - T_{1\text{-}D}(G)` is a factor of the characteristic polynomial.

    :math:`T_{1\text{-}D}` is the same closure's 1-D transmission, evaluated
    at the total streaming :math:`G = \sum_a g_a`.  Read off
    :func:`transmission_spectrum`'s factors, so no degree-:math:`n`
    substitution is expanded.  Along the diagonal :math:`g_a = t \to \infty`
    (the optically thin cell) this root tends to 1.  It is not the only root
    that can approach 1 off the diagonal: LD at :math:`d = 3`,
    :math:`g = (10^6, 1, 1)`, has four roots above 0.99998 (`[M]` qa,
    2026-09-23).  Returns :math:`T_{1\text{-}D}(G)`.
    """
    g = streaming_coefficients(d)
    conserved = _one_dimensional_transmission(closure, sum(g))
    _require(
        _canonical(LAMBDA - conserved, g) in dict(transmission_spectrum(closure, d).factors),
        f"{closure.value}, d={d}: lambda - T_1D(G) is a factor of the characteristic polynomial",
    )
    t = sp.Symbol("t", positive=True)
    _require(
        sp.limit(conserved.subs({ga: t for ga in g}), t, sp.oo) == 1,
        f"{closure.value}, d={d}: along g_a = t -> oo the conserved eigenvalue tends to 1",
    )
    return conserved


# ── claim 5: LD's spectrum at d = 3 is generated by d = 1 and d = 2 ──


def _ld_quadratic() -> tuple[sp.Expr, tuple[sp.Symbol, sp.Symbol]]:
    """The one quadratic factor of LD's d=2 characteristic polynomial, in (g1, g2)."""
    spectrum = transmission_spectrum(Closure.LINEAR_DISCONTINUOUS, 2)
    quadratics = [fac for fac, _ in spectrum.factors if fac.degree(LAMBDA) == 2]
    _require(len(quadratics) == 1, "LD d=2 has exactly one quadratic factor")
    return quadratics[0].as_expr(), streaming_coefficients(2)


def _one_versus_rest_splits(d: int) -> tuple[int, ...]:
    """The axes a whose split {a} | rest is distinct: every axis at d = 3, one at d = 2."""
    return (0,) if d == 2 else tuple(range(d))


def derive_ld_factorisation(d: int = 3) -> dict[sp.Poly, int]:
    r"""LD's characteristic polynomial from its 1-D and 2-D factors.

    Proves that the irreducible factorisation of
    :math:`\det(\lambda I - T_{\rm LD})` over :math:`\mathbb{Z}[g]` is, with
    multiplicity,
    :math:`\lambda^{z}\,\ell(G)\prod_{\{a\}\mid\text{rest}} q(g_a, G - g_a)`,
    with :math:`\ell` the conserved factor (:func:`derive_conserved_mode`),
    :math:`q` the quadratic of :math:`d = 2`, the product over the distinct
    splits of the axes into one axis and the rest (one at :math:`d = 2`,
    since :math:`q` is symmetric; three at :math:`d = 3`), and :math:`z`
    filling the degree.  Returns the factors with multiplicity.  Stated for
    :math:`d \in \{2, 3\}`; :math:`d \ge 4` is not claimed.
    """
    if d not in (2, 3):
        raise ValueError(f"the factorisation is claimed for d in {{2, 3}}, got {d}")
    g = streaming_coefficients(d)
    total = sum(g)
    q, (x, y) = _ld_quadratic()
    _require(
        sp.expand(q - q.subs({x: y, y: x}, simultaneous=True)) == 0,
        "the d=2 quadratic is symmetric in its two arguments",
    )
    expected: dict[sp.Poly, int] = {}
    for a in _one_versus_rest_splits(d):
        key = _canonical(q.subs({x: g[a], y: total - g[a]}, simultaneous=True), g)
        expected[key] = expected.get(key, 0) + 1
    conserved = _canonical(LAMBDA - derive_conserved_mode(Closure.LINEAR_DISCONTINUOUS, d), g)
    expected[conserved] = expected.get(conserved, 0) + 1
    n = cell_response(Closure.LINEAR_DISCONTINUOUS, d).n_face_moments
    zero = _canonical(LAMBDA, g)
    expected[zero] = n - sum(f.degree(LAMBDA) * m for f, m in expected.items())

    _require(
        dict(transmission_spectrum(Closure.LINEAR_DISCONTINUOUS, d).factors) == expected,
        f"LD d={d}: the factorisation is lambda^z l(G) prod q(g_a, G - g_a)",
    )
    return expected


# ── claim 6: the stability function ──


PADE_ORDER: dict[Closure, tuple[int, int]] = {
    Closure.STEP: (0, 1),
    Closure.DIAMOND: (1, 1),
    Closure.LINEAR_DISCONTINUOUS: (1, 2),
}
r"""The :math:`[p/q]` Padé type of each closure's 1-D transmission in :math:`\tau`."""


@dataclass(frozen=True)
class StabilityFunction:
    r"""A closure's 1-D transmission :math:`R(\tau)`, read as a one-step method's stability function."""

    closure: Closure
    transmission: sp.Expr
    pade_order: tuple[int, int]
    value_at_infinity: sp.Rational

    @property
    def l_stable(self) -> bool:
        r"""A-stable (proved by :func:`derive_stability_function`) and :math:`R(\infty) = 0`."""
        return self.value_at_infinity == 0


def derive_stability_function(closure: Closure) -> StabilityFunction:
    r"""The 1-D transmission is an A-stable Padé approximant of :math:`e^{-\tau}` with :math:`R(\infty) = f`.

    Proves: (i) :math:`R` is of type :math:`[p/q]` and agrees with
    :math:`e^{-\tau}` through :math:`\tau^{p+q}` and not beyond, so it is
    the Padé approximant (uniqueness of the Padé table); (ii) **A-stability**,
    :math:`|R(\tau)| \le 1` for every complex :math:`\tau` with
    :math:`\operatorname{Re}\tau \ge 0`, by the E-polynomial criterion
    (Hairer & Wanner, Vol. II, §IV.3, reasoned citation): with
    :math:`R = P/Q`, every pole lies in :math:`\operatorname{Re}\tau < 0`
    (:math:`Q` Hurwitz; for degree :math:`\le 2`, all coefficients of one
    sign) and :math:`E(y) = |Q(iy)|^2 - |P(iy)|^2 \ge 0` for real :math:`y`
    (a polynomial in :math:`y^2` with non-negative coefficients), so the
    maximum modulus principle bounds :math:`|R|` by 1 on the right
    half-plane; (iii) :math:`|R(\tau)| < 1` strictly for every real
    :math:`\tau > 0` (both :math:`1 - R` and :math:`1 + R` certified
    positive); (iv) :math:`R(\infty)` equals the feedthrough :math:`f`;
    (v) in every dimension :math:`d \in \{1,2,3\}`, :math:`T \to fI` in
    the optically thick cell (:math:`g \to 0`).
    """
    tau = sp.Symbol("tau", positive=True)
    R = _one_dimensional_transmission(closure, 1 / tau)
    num, den = sp.fraction(sp.cancel(R))
    p, q = PADE_ORDER[closure]
    _require((sp.degree(num, tau), sp.degree(den, tau)) == (p, q), f"{closure.value}: R is of type [{p}/{q}]")
    error = sp.series(R - sp.exp(-tau), tau, 0, p + q + 2).removeO()
    _require(
        all(error.coeff(tau, k) == 0 for k in range(p + q + 1)),
        f"{closure.value}: R agrees with exp(-tau) through tau^{p + q}",
    )
    _require(error.coeff(tau, p + q + 1) != 0, f"{closure.value}: and not beyond")
    _require(_is_a_stable(num, den, tau), f"{closure.value}: R is A-stable (Hurwitz denominator, E(y) >= 0)")
    for bound in (1 - R, 1 + R):
        top, bottom = sp.fraction(sp.together(bound))
        _require(
            _is_positive_on_orthant(sp.expand(top * bottom), (tau,)),
            f"{closure.value}: |R(tau)| < 1 for every tau > 0",
        )
    at_infinity = sp.limit(R, tau, sp.oo)
    f = cell_response(closure, 1).feedthrough_scalar
    _require(at_infinity == f, f"{closure.value}: R(oo) = {at_infinity} is the feedthrough {f}")
    for d in (1, 2, 3):
        response = cell_response(closure, d)
        thick = response.transmission(at={ga: 0 for ga in response.g})
        _require(
            thick == f * sp.eye(response.n_face_moments),
            f"{closure.value}, d={d}: T -> f I in the optically thick cell",
        )
    return StabilityFunction(closure, R, (p, q), sp.Rational(at_infinity))


# ── claim 7: conservation and flat-flux preservation, from the four blocks ──


def _face_average_rows(response: CellResponse) -> list[int]:
    """The index of each face's average moment within the axis-major face-moment vector."""
    per_face = response.n_face_moments // response.d
    return [a * per_face for a in range(response.d)]


def derive_particle_balance(closure: Closure, d: int, at: Point | None = None) -> None:
    r"""Outflow minus inflow plus absorption equals the source, for every inflow and source.

    With :math:`\bar\cdot` the average moment, for every
    :math:`\psi^{\rm in}` and :math:`q`,

    .. math::

        \sum_a g_a\big(\bar\psi^{\rm out}_a - \bar\psi^{\rm in}_a\big)
        + \bar c \;=\; \bar q ,

    where :math:`\psi^{\rm out} = T\psi^{\rm in} + (R A^{-1} S)\,q` and
    :math:`\bar c` is the cell average of
    :math:`A^{-1} E\,\psi^{\rm in} + A^{-1} S\,q`: all four blocks enter.
    Proved identically in :math:`g`, or at the point ``at`` (LD at
    :math:`d = 3`, whose symbolic blocks take a minute).
    """
    response = cell_response(closure, d)
    g = [ga if at is None else at[ga] for ga in response.g]
    rows = _face_average_rows(response)
    T, esc = response.transmission(at), response.escape(at)
    to_cell, from_source = response.inflow_to_cell(at), response.source_to_cell(at)
    for k in range(response.n_face_moments):
        leak = sum(g[a] * (T[rows[a], k] - int(k == rows[a])) for a in range(d))
        _require(sp.cancel(leak + to_cell[0, k]) == 0, f"{closure.value}, d={d}: balance for inflow moment {k}")
    for k in range(response.source_mass.cols):
        leak = sum(g[a] * esc[rows[a], k] for a in range(d))
        _require(
            sp.cancel(leak + from_source[0, k] - int(k == 0)) == 0,
            f"{closure.value}, d={d}: balance for source moment {k}",
        )


def derive_flat_flux_preservation(closure: Closure, d: int, at: Point | None = None) -> None:
    r"""A uniform inflow equal to a uniform source passes through unchanged.

    The infinite-medium solution :math:`\psi \equiv q/\Sigma_t` is a fixed
    point of the cell: with every face average of :math:`\psi^{\rm in}`
    equal to 1 (every higher moment 0) and :math:`q = 1`, the outflow is the
    same uniform vector and the cell state is the uniform function.  Uses
    all four blocks.
    """
    response = cell_response(closure, d)
    uniform_face = sp.zeros(response.n_face_moments, 1)
    for r in _face_average_rows(response):
        uniform_face[r] = 1
    unit_source = sp.zeros(response.source_mass.cols, 1)
    unit_source[0] = 1
    outflow = response.transmission(at) * uniform_face + response.escape(at) * unit_source
    _require(
        outflow.applyfunc(sp.cancel) == uniform_face,
        f"{closure.value}, d={d}: the uniform flux passes through the faces",
    )
    cell = response.inflow_to_cell(at) * uniform_face + response.source_to_cell(at) * unit_source
    uniform_cell = sp.zeros(response.n_cell_moments, 1)
    uniform_cell[0] = 1
    _require(cell.applyfunc(sp.cancel) == uniform_cell, f"{closure.value}, d={d}: the cell state is uniform")


# ── claim 8: the page's closed forms, in the page's variables ──


def derive_page_closed_form(closure: Closure, d: int) -> sp.Matrix:
    r"""The page's closed form of ``closure``'s transmission equals the derived one.

    In the page's variables, diamond is :math:`(2/D)\mathbf 1 w^{\mathsf T} - I`
    with :math:`w_a = 2|\mu_a|A_a` and :math:`D = \Sigma_t V + \sum_b w_b`,
    and step is :math:`(1/D')\mathbf 1 w'^{\mathsf T}` with
    :math:`w'_a = |\mu_a| A_a` and :math:`D' = \Sigma_t V + \sum_b w'_b`.
    Per unit volume with :math:`\Sigma_t = 1`, :math:`|\mu_a| A_a / V = g_a`,
    so :math:`w_a \to 2g_a`, :math:`D \to 1 + 2G`, :math:`w'_a \to g_a` and
    :math:`D' \to 1 + G`.  Returns the page's matrix.  LD has no closed form
    on the page, so it is refused.
    """
    if closure is Closure.LINEAR_DISCONTINUOUS:
        raise ValueError(f"the page states a closed form for step and diamond only, not {closure.value}")
    g = sp.Matrix(streaming_coefficients(d))
    total = sum(g)
    ones = sp.ones(d, 1)
    if closure is Closure.DIAMOND:
        closed_form = (2 / (1 + 2 * total)) * ones * (2 * g).T - sp.eye(d)
    else:
        closed_form = (1 / (1 + total)) * ones * g.T
    _require(
        (cell_response(closure, d).transmission() - closed_form).applyfunc(sp.cancel) == sp.zeros(d, d),
        f"{closure.value}, d={d}: the derived transmission is the page's closed form",
    )
    return closed_form


# ── the whole comparison ──


def derive_comparison(dims: tuple[int, ...] = (1, 2, 3)) -> dict[tuple[Closure, int], TransmissionSpectrum]:
    r"""Every closure at every dimension: diamond alone leaves a face mode undamped.

    Diamond's undamped set is :math:`\{-1\}^{d-1}` (empty at :math:`d = 1`);
    step's and LD's are empty at every :math:`d`.
    """
    table = {}
    for d in dims:
        for closure in Closure:
            spectrum = transmission_spectrum(closure, d)
            expected = ((sp.Integer(-1), d - 1),) if closure is Closure.DIAMOND and d > 1 else ()
            _require(
                spectrum.undamped == expected,
                f"{closure.value}, d={d}: undamped set {spectrum.undamped}, expected {expected}",
            )
            table[(closure, d)] = spectrum
    return table


if __name__ == "__main__":
    import time

    t0 = time.time()
    for d in (1, 2, 3):
        for closure in Closure:
            derive_realization(closure, d)
    print(f"realizations proved, d = 1, 2, 3  ({time.time() - t0:.1f}s)")
    for (closure, d), spec in derive_comparison().items():
        print(f"{closure.value:22s} d={d}  undamped={spec.undamped}  feedthrough x{spec.feedthrough_multiplicity}")
    for d in (2, 3):
        derive_ld_factorisation(d)
    print(f"LD factorisation proved, d = 2, 3  ({time.time() - t0:.1f}s)")
    point = {1: (sp.Rational(7, 10),), 2: (sp.Rational(7, 10), sp.Rational(13, 11)),
             3: (sp.Rational(7, 10), sp.Rational(13, 11), sp.Rational(2, 5))}
    for d in (1, 2, 3):
        for closure in Closure:
            derive_conserved_mode(closure, d)
            symbolic = d <= 2 or closure is not Closure.LINEAR_DISCONTINUOUS
            at = None if symbolic else dict(zip(streaming_coefficients(d), point[d]))
            derive_particle_balance(closure, d, at)
            derive_flat_flux_preservation(closure, d, at)
        for closure in (Closure.DIAMOND, Closure.STEP):
            derive_page_closed_form(closure, d)
    print(f"conserved modes, balance, flat flux, closed forms proved  ({time.time() - t0:.1f}s)")
    for closure in Closure:
        sf = derive_stability_function(closure)
        print(f"{closure.value:22s} R(tau) = {sf.transmission}  Pade {sf.pade_order}  R(oo) = {sf.value_at_infinity}")
    print(f"total {time.time() - t0:.1f}s")
