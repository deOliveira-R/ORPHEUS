r"""The Lambertian, FACTORED — ``R = B ∘ C`` through the scalar-current space.

G6.3 step 3 (#330). The white/albedo response was one welded operator
:math:`\Gamma_+ \to \Gamma_-` that **withheld its transpose**
(``is_adjointable = False``, deferred to boundary phase B5) because the
Euclidean and cosine-weighted readings differ and it declined to advertise the
ambiguous one. Factoring it removes the ambiguity rather than resolving it: each
factor has ONE honest transpose, and the composite's Hilbert adjoint then falls
out of the spaces.

.. math::

    \Gamma_+(f) \;\xrightarrow{\;C\;}\; S(f)
                \;\xrightarrow{\;B\;}\; \Gamma_-(f),
    \qquad
    R^{\mathsf T} = C^{\mathsf T}B^{\mathsf T},
    \qquad
    R^{*} = G_+^{-1}R^{\mathsf T}G_- .

⭐ **This module carries the gate that justified G6** — the weighted adjoint law
on a boundary operator that does NOT commute with the trace metric. It was
unwritable before this step for two independent reasons: the operator had no
transpose, and its spaces were unbound so ``.H`` degraded to the bare Euclidean
transpose. Its control leg is the specular case, where the commutator vanishes
and weighted ≡ Euclidean; that leg lives with the specular operator, and the
contrast is stated in :func:`test_the_metric_is_NOT_invisible_here`.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.face_layout import FaceLayout
from orpheus.numerics.operator import IncompatibleOperatorComposition
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces import AngularTraceSpace
from orpheus.sn.boundary.angular import (
    AngularAverageOperator,
    IsotropicEmissionOperator,
    PartialCurrentOperator,
)

pytestmark = pytest.mark.l0

_QUADRATURES = {
    "gauss_legendre(4)": lambda: Quadrature.gauss_legendre(n_ordinates=4),
    "product(4,4)": lambda: Quadrature.product(n_mu=4, n_phi=4),
    "lebedev(17)": lambda: Quadrature.lebedev(order=17),
}
_ALL = list(_QUADRATURES)
_SEED = 20260804


class _Lambertian:
    """The factored chain and its welded predecessor, on one face."""

    def __init__(self, quad_name: str, ng: int = 3, face: str = "xmin") -> None:
        quad = _QUADRATURES[quad_name]()
        layout = FaceLayout.from_named_shapes(
            [("xmin", (int(quad.N), ng)), ("xmax", (int(quad.N), ng))]
        )
        self.trace = AngularTraceSpace.from_quadrature_and_layout(quad, layout)
        self.face = face
        self.gamma_plus = self.trace.outflow_space(face)
        self.gamma_minus = self.trace.inflow_space(face)
        self.current = self.trace.current_space(face)

        full_metric = np.asarray(self.trace.face_space(face).inner_product_weights)
        self.cos_w = full_metric[self.trace.outflow_indices_for_face(face)]
        self.norm = float(self.cos_w.sum())

        self.C = PartialCurrentOperator(
            self.cos_w, domain=self.gamma_plus, codomain=self.current,
        )
        self.B = IsotropicEmissionOperator(
            self.norm, self.gamma_minus.shape[0],
            domain=self.current, codomain=self.gamma_minus,
        )
        self.R = self.B @ self.C
        self.welded = AngularAverageOperator(
            cos_w=self.cos_w, norm=self.norm,
            n_inflow=self.gamma_minus.shape[0],
        )

    def draw(self, seed_offset: int = 0):
        rng = np.random.default_rng(_SEED + seed_offset)
        return (
            rng.standard_normal(self.gamma_plus.shape),
            rng.standard_normal(self.gamma_minus.shape),
        )


def _dense(op, in_shape) -> np.ndarray:
    """Materialise ``op`` by applying it to the standard basis."""
    n_in = int(np.prod(in_shape))
    cols = []
    for i in range(n_in):
        e = np.zeros(n_in)
        e[i] = 1.0
        cols.append(np.asarray(op.apply(e.reshape(in_shape))).reshape(-1))
    return np.column_stack(cols)


# ─────────────────────────────────────────────────────────────────────
# The factoring changed nothing
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_factored_chain_is_BIT_IDENTICAL_to_the_welded_operator(quad_name):
    r"""``(B ∘ C)ψ == R_diff ψ``, exactly — factoring is a re-spelling.

    Bit-identity is achievable here (not merely principled equivalence) because
    the operation order is unchanged: the welded body computes
    ``(cos_w·ψ).sum() / norm`` and the chain computes the same sum in ``C`` and
    the same division in ``B``. If a future refactor moves the normalisation
    into ``C``, this gate SHOULD break — the reduction would then happen after a
    per-element divide, which is a different floating-point computation.
    """
    lam = _Lambertian(quad_name)
    psi, _ = lam.draw()
    np.testing.assert_array_equal(lam.R.apply(psi), lam.welded.apply(psi))


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_composite_is_typed_end_to_end(quad_name):
    r"""``R : Γ₊(f) → Γ₋(f)``, with ``S(f)`` internal to the chain."""
    lam = _Lambertian(quad_name)
    assert lam.R.domain is lam.gamma_plus
    assert lam.R.codomain is lam.gamma_minus
    assert lam.C.codomain is lam.B.domain is lam.current


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_chain_REFUSES_to_compose_the_wrong_way_round(quad_name):
    r"""``C @ B`` does not connect — ``S(f) != Γ₋(f)``."""
    lam = _Lambertian(quad_name)
    with pytest.raises(IncompatibleOperatorComposition, match="angular_trace"):
        _ = lam.C @ lam.B


# ─────────────────────────────────────────────────────────────────────
# The transpose that did not exist
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_welded_operator_withheld_the_transpose_the_chain_supplies(quad_name):
    """The step's whole point, as a before/after on one line each.

    Not decoration: it pins that the capability genuinely changed, so a future
    reader can see the deferral was closed rather than forgotten.
    """
    lam = _Lambertian(quad_name)
    assert lam.welded.is_adjointable is False
    assert lam.R.is_adjointable is True
    assert lam.C.is_adjointable is True
    assert lam.B.is_adjointable is True


@pytest.mark.verifies("bc-response-factored-adjoint")
@pytest.mark.parametrize("quad_name", _ALL)
def test_the_transpose_matches_the_dense_transpose_of_apply(quad_name):
    r"""``Rᵀ`` verified against ``dense(R).T`` — an independent construction.

    ``apply_transpose`` is derived analytically (``Rᵀφ = cos_w·Σφ / norm``);
    the dense transpose is built by applying ``R`` to the standard basis and
    transposing the matrix. Neither uses the other, so agreement is evidence
    rather than tautology.
    """
    lam = _Lambertian(quad_name)
    _, phi = lam.draw()

    dense_t = _dense(lam.R, lam.gamma_plus.shape).T
    got = np.asarray(lam.R.apply_transpose(phi)).reshape(-1)
    expected = dense_t @ np.asarray(phi).reshape(-1)
    np.testing.assert_allclose(got, expected, rtol=0.0, atol=1e-15)

    # and the closed form the factorisation predicts
    closed = (
        lam.cos_w.reshape((-1,) + (1,) * (phi.ndim - 1))
        * np.asarray(phi).sum(axis=0) / lam.norm
    )
    np.testing.assert_allclose(
        lam.R.apply_transpose(phi), closed, rtol=0.0, atol=1e-15,
    )


# ─────────────────────────────────────────────────────────────────────
# ⭐ THE FLAGSHIP — the weighted adjoint law on a metric-sensitive operator
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.verifies("bc-response-factored-adjoint")
@pytest.mark.parametrize("quad_name", _ALL)
def test_the_WEIGHTED_adjoint_law_holds(quad_name):
    r"""⭐ ``⟨Rx, y⟩_{Γ₋} == ⟨x, R*y⟩_{Γ₊}`` — the gate that justified G6.

    Unwritable before this step for TWO independent reasons: the operator had no
    transpose at all, and its spaces were unbound so ``.H`` silently degraded to
    the Euclidean transpose — which does **not** satisfy this identity here (see
    the companion below). This is the shape of gate that would have caught
    ERR-067, where a wrong/absent metric put the error class inside the measured
    functional's own stabiliser.
    """
    lam = _Lambertian(quad_name)
    x, y = lam.draw()

    lhs = lam.gamma_minus.inner_product(lam.R.apply(x), y)
    rhs = lam.gamma_plus.inner_product(x, lam.R.H.apply(y))
    assert lhs == pytest.approx(rhs, rel=1e-13)


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_metric_is_NOT_invisible_here(quad_name):
    r"""⭐ The Lambertian does NOT commute with the trace metric — measured.

    The negative leg the flagship needs, and the reason this operator was chosen
    for it. Every OTHER shipped boundary law is Mode-12 blind: vacuum is a zero
    map, reflective and albedo-specular are weight-preserving permutations,
    periodic has opposite normals, prescribed-inflow is rank-0. For all of them
    ``[G, Aᵀ] = 0`` and dropping the metric is invisible, so a gate built on one
    of those could never distinguish "metric applied" from "metric dropped".

    `[M]` here ``‖A† − Aᵀ‖ / ‖Aᵀ‖`` is 0.209 / 0.684 / 0.612 on
    ``gauss_legendre(4)`` / ``product(4,4)`` / ``lebedev(17)`` — an O(1)
    separation, not a tolerance question.
    """
    lam = _Lambertian(quad_name)
    _, phi = lam.draw()

    hilbert = np.asarray(lam.R.H.apply(phi))
    euclidean = np.asarray(lam.R.apply_transpose(phi))
    relative = np.abs(hilbert - euclidean).max() / np.abs(euclidean).max()

    assert relative > 0.1, (
        f"the Hilbert adjoint and the Euclidean transpose differ by only "
        f"{relative:.2e} relative — this fixture can no longer tell an applied "
        f"metric from a dropped one, and the flagship gate above is blind"
    )


@pytest.mark.verifies("bc-response-factored-adjoint")
@pytest.mark.parametrize("quad_name", _ALL)
def test_dropping_the_binding_BREAKS_the_weighted_law(quad_name):
    r"""⭐ The mutation, run as a gate: unbound spaces fail the adjoint law.

    Re-creating the SAME chain with no spaces reproduces exactly the pre-G6.3
    state — ``.H`` degrades to the Euclidean transpose with no error and no
    warning — and the weighted law then FAILS. This is the standing proof that
    the binding is load-bearing rather than decorative, and it is the mutation
    ``domain := None`` written as a permanent test instead of a one-off probe.
    """
    lam = _Lambertian(quad_name)
    x, y = lam.draw()

    unbound = (
        IsotropicEmissionOperator(lam.norm, lam.gamma_minus.shape[0])
        @ PartialCurrentOperator(lam.cos_w)
    )
    assert unbound.domain is None and unbound.codomain is None

    # the forward map is unchanged — only the adjoint degrades
    np.testing.assert_array_equal(unbound.apply(x), lam.R.apply(x))

    lhs = lam.gamma_minus.inner_product(unbound.apply(x), y)
    rhs = lam.gamma_plus.inner_product(x, unbound.H.apply(y))
    assert lhs != pytest.approx(rhs, rel=1e-6), (
        "the unbound chain satisfied the weighted adjoint law, so this fixture "
        "cannot demonstrate that binding matters"
    )


# ─────────────────────────────────────────────────────────────────────
# Each factor on its own
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("quad_name", _ALL)
@pytest.mark.verifies("bc-response-factored-adjoint")
@pytest.mark.parametrize("which", ["C", "B"])
def test_each_factor_satisfies_its_own_weighted_adjoint_law(quad_name, which):
    r"""``⟨Fx, y⟩_cod == ⟨x, F*y⟩_dom`` for each link separately.

    The composite's law could hold by a cancellation of two factor-level errors;
    gating the links independently removes that escape.
    """
    lam = _Lambertian(quad_name)
    op = lam.C if which == "C" else lam.B
    domain, codomain = op.domain, op.codomain
    assert domain is not None and codomain is not None, "step 3 binds both"

    rng = np.random.default_rng(_SEED + 11)
    x = rng.standard_normal(domain.shape)
    y = rng.standard_normal(codomain.shape)

    assert codomain.inner_product(op.apply(x), y) == pytest.approx(
        domain.inner_product(x, op.H.apply(y)), rel=1e-13,
    )


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_contraction_produces_the_PARTIAL_CURRENT(quad_name):
    r"""``C`` is not an average — it is :math:`J^+`, and the units say so.

    The welded operator divided by ``norm`` inside the contraction, so its
    intermediate was a mean intensity with no name in the literature. Splitting
    the division into ``B`` (where it performs the current→intensity unit
    change) leaves ``C`` producing the outgoing partial current itself, which is
    what an albedo law scales: :math:`J^- = \alpha J^+`.
    """
    lam = _Lambertian(quad_name)
    psi = np.ones(lam.gamma_plus.shape)

    # J+ of a unit isotropic outflow is exactly the cosine-weighted measure.
    np.testing.assert_allclose(
        lam.C.apply(psi), np.full(lam.current.shape, lam.norm),
        rtol=0.0, atol=1e-13,
    )
    # and it is NOT the mean (which would be 1.0 for this input)
    assert not np.allclose(lam.C.apply(psi), 1.0)


@pytest.mark.parametrize("quad_name", _ALL)
def test_the_chain_is_CONSERVATIVE_on_a_symmetric_quadrature(quad_name):
    r"""``J⁻ == J⁺`` — the normalisation's reason for existing.

    Bell & Glasstone 1970 §1.5: the cosine-weighted outgoing current equals the
    cosine-weighted incoming current. Measured through the chain, this is
    ``Σ_{Γ₋} cos_w⁻ · (Rψ) == Σ_{Γ₊} cos_w⁺ · ψ``, which holds when the two
    hemispheres carry equal cosine-weighted measure — a property of the
    quadrature, not of the operator, and worth pinning because the
    normalisation silently assumes it.
    """
    lam = _Lambertian(quad_name)
    face_metric = np.asarray(lam.trace.face_space(lam.face).inner_product_weights)
    cos_w_in = face_metric[lam.trace.inflow_indices_for_face(lam.face)]

    psi = np.ones(lam.gamma_plus.shape)
    j_out = lam.C.apply(psi)
    j_in = (cos_w_in.reshape((-1,) + (1,) * (psi.ndim - 1)) * lam.R.apply(psi)).sum(axis=0)

    assert float(cos_w_in.sum()) == pytest.approx(lam.norm, rel=1e-13), (
        "the hemispheres carry unequal cosine-weighted measure on this "
        "quadrature, so the normalisation is not conservative and this gate "
        "is asserting something the operator does not promise"
    )
    np.testing.assert_allclose(j_in, j_out, rtol=1e-13, atol=0.0)


# ─────────────────────────────────────────────────────────────────────
# Construction guards
# ─────────────────────────────────────────────────────────────────────


def test_the_contraction_refuses_a_full_face_weight_vector():
    """A zero entry means Γ₊ was not restricted — the B3.4a twin-classifier trap.

    Pre-B3.4a the welded operator carried ``N`` weights zeroed off the outgoing
    hemisphere by its OWN test, a second outflow classifier that disagreed with
    the trace space's. Refusing non-positive entries keeps that unspellable here
    rather than merely fixed.
    """
    with pytest.raises(ValueError, match="STRICTLY positive"):
        PartialCurrentOperator(np.array([0.5, 0.0, 0.3]))


def test_the_emission_refuses_a_degenerate_normalisation():
    with pytest.raises(ValueError, match="strictly positive"):
        IsotropicEmissionOperator(0.0, 4)


def test_the_emission_refuses_an_empty_codomain():
    with pytest.raises(ValueError, match="n_inflow must be positive"):
        IsotropicEmissionOperator(1.0, 0)
