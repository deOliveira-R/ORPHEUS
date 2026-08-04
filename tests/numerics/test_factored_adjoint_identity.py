r"""The factored-adjoint identity — a composed operator's Hilbert adjoint.

Permanent home for the theorem the boundary algebra of record states as
:ref:`bc-response-factored-adjoint`. It was derived in a design session for
campaign step **G6.3** (`.claude/plans/geometric_transformation_consolidation.md`
§7h) and lived only as a scratch probe; that made the corpus claim
reproducible-in-principle but not *gated*, which is the weaker thing.

**The claim.** For a response factored as :math:`R = B \circ C` with
:math:`C : \Gamma_+ \to S` and :math:`B : S \to \Gamma_-`,

.. math::

    R^* = C^* B^*
        = \bigl(G_+^{-1} C^{\mathsf T} G_S\bigr)
          \bigl(G_S^{-1} B^{\mathsf T} G_-\bigr)
        = G_+^{-1} C^{\mathsf T} B^{\mathsf T} G_-
        = G_+^{-1} R^{\mathsf T} G_-,

so **the intermediate space's metric cancels** — provided it is invertible.

⭐ **And it telescopes at ANY chain length**, which is the general statement:

.. math::

    (A_1 \cdots A_n)^* = A_n^* \cdots A_1^*
                       = G_0^{-1} (A_1 \cdots A_n)^{\mathsf T} G_n .

**Every** interior metric cancels, not just one. That is what licenses modelling
a boundary law as a **sequence from outflow to inflow** — an expression in the
operator algebra (``@`` for the chain, ``+`` for parallel channels such as a
partly-polished wall) rather than a bespoke two-factor struct. A deck
transformation is then not an "atomic" special arm but a chain of **length 1**,
gated as such below; and no interior space needs a declared type, only a
non-degenerate metric.

**Why it matters.** It is the difference between "a response's adjoint exists"
and "a response's adjoint exists if we choose the intermediate metric
correctly". The answer is the former, and the *only* requirement on the
intermediate space is non-degeneracy. That is what makes
:class:`~orpheus.sn.boundary.angular.AngularAverageOperator`'s missing transpose
a consequence of factoring rather than a hand-rolled addition.

**What is the SUT.** Not hand-rolled matrix algebra: every adjoint below is
computed through :class:`~orpheus.numerics.space.FunctionSpace`'s **production**
:meth:`~orpheus.numerics.space.FunctionSpace.apply_metric` /
:meth:`~orpheus.numerics.space.FunctionSpace.apply_inverse_metric` — the same
two calls :class:`~orpheus.numerics.operator._AdjointOperator` makes. The dense
``C`` / ``B`` are the *inputs*, the metric machinery is what is under test.

Every fixture uses :math:`|\Gamma_+| \neq |\Gamma_-|` deliberately. The shipped
half-traces happen to have equal size on every quadrature — an accident, not a
contract — and a square fixture would let a genuinely wrong adjoint typecheck.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.space import FunctionSpace

pytestmark = pytest.mark.l0

#: Deliberately unequal — see the module docstring.
N_OUT, N_IN = 7, 5

_SEED = 20260804

#: The trace metric on each half. ⭐ ``_METRIC_PLUS`` is ALSO the contraction
#: weight, and that coupling is physics, not convenience: the Lambertian
#: contracts the outflow against :math:`|\Omega\cdot\hat n|\,w`, which IS the
#: partial-current metric the trace space installs. The structural-symmetry
#: gate below holds *because* of this identification and would fail without
#: it — see its docstring.
_METRIC_PLUS = np.random.default_rng(_SEED).uniform(0.2, 1.0, N_OUT)
_METRIC_MINUS = np.random.default_rng(_SEED + 1).uniform(0.2, 1.0, N_IN)


def _spaces(g_intermediate: float) -> tuple[FunctionSpace, FunctionSpace, FunctionSpace]:
    r"""``(Γ₊, S, Γ₋)`` with distinct, strictly-positive diagonal metrics.

    ``g_intermediate`` is the scalar metric installed on :math:`S`; the whole
    point of the theorem is that the answer does not depend on it.
    """
    return (
        FunctionSpace(
            name="probe:gamma_plus", shape=(N_OUT,),
            inner_product_weights=_METRIC_PLUS,
        ),
        FunctionSpace(
            name="probe:S", shape=(1,),
            inner_product_weights=np.array([g_intermediate]),
        ),
        FunctionSpace(
            name="probe:gamma_minus", shape=(N_IN,),
            inner_product_weights=_METRIC_MINUS,
        ),
    )


def _contraction_and_broadcast() -> tuple[np.ndarray, np.ndarray]:
    r"""``(C, B)`` — the Lambertian's shape: cosine contraction, then broadcast.

    ``C`` collapses :math:`\Gamma_+` to one number against the outflow trace
    metric, normalised; ``B`` writes that number to every :math:`\Gamma_-` row.
    Their product is rank one, which is the structure the boundary algebra's
    rank-one theorem is about.
    """
    contraction = (_METRIC_PLUS / _METRIC_PLUS.sum())[None, :]   # (1, N_OUT)
    broadcast = np.ones((N_IN, 1))                               # (N_IN, 1)
    return contraction, broadcast


def _hilbert_adjoint(
    matrix: np.ndarray, domain: FunctionSpace, codomain: FunctionSpace,
) -> np.ndarray:
    r"""``A* = G_dom⁻¹ Aᵀ G_cod``, built column-by-column THROUGH THE SPACES.

    Deliberately routed via ``codomain.apply_metric`` and
    ``domain.apply_inverse_metric`` rather than ``np.diag`` algebra: those are
    the two production calls ``_AdjointOperator.apply`` makes, so this gate
    fails if the metric machinery regresses — which is the point of putting it
    here rather than leaving it in a probe.
    """
    columns = []
    for i in range(matrix.shape[0]):
        basis = np.zeros(matrix.shape[0])
        basis[i] = 1.0
        weighted = codomain.apply_metric(basis)
        columns.append(domain.apply_inverse_metric(matrix.T @ weighted))
    return np.column_stack(columns)


# ─────────────────────────────────────────────────────────────────────
# The identity, and its independence of the intermediate metric
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.verifies("bc-response-factored-adjoint")
@pytest.mark.parametrize(
    "g_intermediate",
    [1.0, 1e-6, 3.7e5, 2.0, 1e3],
    ids=["euclidean", "tiny", "huge", "modest", "large"],
)
def test_the_factored_adjoint_equals_the_direct_one(g_intermediate):
    r"""``C* B* == G₊⁻¹RᵀG₋`` for every intermediate metric — it CANCELS.

    The parametrisation spans eleven orders of magnitude in :math:`G_S`. If the
    cancellation were not exact this gate would fail at the extremes first,
    which is why ``tiny`` and ``huge`` are present rather than a single
    representative value.
    """
    gamma_plus, scalar, gamma_minus = _spaces(g_intermediate)
    contraction, broadcast = _contraction_and_broadcast()
    response = broadcast @ contraction

    direct = _hilbert_adjoint(response, gamma_plus, gamma_minus)
    factored = (
        _hilbert_adjoint(contraction, gamma_plus, scalar)
        @ _hilbert_adjoint(broadcast, scalar, gamma_minus)
    )

    assert factored.shape == (N_OUT, N_IN)
    np.testing.assert_allclose(factored, direct, rtol=0.0, atol=1e-15)


@pytest.mark.verifies("bc-response-factored-adjoint")
def test_the_adjoint_law_holds_in_the_WEIGHTED_pairing():
    r"""``⟨Rx, y⟩_{Γ₋} == ⟨x, R*y⟩_{Γ₊}`` — the definition being claimed.

    The gate above shows two *constructions* agree; this shows the construction
    is the adjoint at all. Both legs are needed: a systematically wrong
    :math:`G_+^{-1}\cdot G_-` would satisfy the first and fail this one.
    """
    gamma_plus, _, gamma_minus = _spaces(1.0)
    contraction, broadcast = _contraction_and_broadcast()
    response = broadcast @ contraction
    adjoint = _hilbert_adjoint(response, gamma_plus, gamma_minus)

    rng = np.random.default_rng(_SEED + 2)
    x = rng.standard_normal(N_OUT)
    y = rng.standard_normal(N_IN)

    assert gamma_minus.inner_product(response @ x, y) == pytest.approx(
        gamma_plus.inner_product(x, adjoint @ y), rel=1e-14,
    )


@pytest.mark.verifies("bc-response-factored-adjoint")
def test_a_DEGENERATE_intermediate_metric_BREAKS_the_cancellation():
    r"""⭐ The negative leg — the precondition is load-bearing, not decorative.

    ``G_S = 0`` is the one case where the cancellation fails, because
    ``apply_inverse_metric`` is the Moore–Penrose pseudo-inverse and returns
    zero on the metric's null space rather than ``1/0``. Without this gate the
    claim "any non-degenerate metric works" has no evidence that
    *non-degenerate* is doing any work — and the whole reason the boundary
    ladder's HALF-traces (strictly positive) can serve as intermediates while
    the FULL per-face tier (zero on tangential ordinates) cannot.
    """
    gamma_plus, degenerate, gamma_minus = _spaces(0.0)
    contraction, broadcast = _contraction_and_broadcast()
    response = broadcast @ contraction

    direct = _hilbert_adjoint(response, gamma_plus, gamma_minus)
    factored = (
        _hilbert_adjoint(contraction, gamma_plus, degenerate)
        @ _hilbert_adjoint(broadcast, degenerate, gamma_minus)
    )

    assert np.allclose(factored, 0.0), (
        "a zero metric should annihilate the factored adjoint via the "
        "pseudo-inverse; if it does not, this gate no longer probes degeneracy"
    )
    assert not np.allclose(factored, direct)


# ─────────────────────────────────────────────────────────────────────
# The general statement: it telescopes at ANY chain length
# ─────────────────────────────────────────────────────────────────────


def _chain(dims: list[int]) -> tuple[list[np.ndarray], list[FunctionSpace]]:
    r"""A chain ``V0 → V1 → … → Vn`` with all-different dims and metrics."""
    rng = np.random.default_rng(_SEED + 7)
    spaces = [
        FunctionSpace(
            name=f"probe:chain:V{i}", shape=(d,),
            inner_product_weights=rng.uniform(0.3, 2.0, d),
        )
        for i, d in enumerate(dims)
    ]
    links = [
        rng.standard_normal((dims[i + 1], dims[i])) for i in range(len(dims) - 1)
    ]
    return links, spaces


@pytest.mark.verifies("bc-response-factored-adjoint")
@pytest.mark.parametrize(
    "dims",
    [[7, 5], [7, 1, 5], [7, 3, 5, 4], [7, 2, 6, 3, 5]],
    ids=["len1-deck", "len2-lambertian", "len3", "len4"],
)
def test_the_chain_adjoint_telescopes_at_any_length(dims):
    r"""⭐ ``(A₁…Aₙ)* = Aₙ*…A₁*`` — EVERY interior metric cancels.

    The general form of the identity above, and the reason a boundary law is an
    **expression in the operator algebra** rather than a bespoke two-factor
    struct: it is a sequence from :math:`\Gamma_+` to :math:`\Gamma_-`, and how
    many links it has is not a type distinction.

    ⭐ **``len1-deck`` is the load-bearing case.** A deck transformation was
    designed as an "atomic" arm *separate* from the composed response; this gate
    pins that it is nothing of the kind — it is the same formula at length 1. If
    a future refactor reintroduces a special path for it, this parametrisation
    is what says the special path is unnecessary.

    Interior spaces therefore need **no declared type** — only a non-degenerate
    metric, which the negative-leg gate above shows is load-bearing.
    """
    links, spaces = _chain(dims)

    composite = links[0]
    for link in links[1:]:
        composite = link @ composite
    direct = _hilbert_adjoint(composite, spaces[0], spaces[-1])

    chained = _hilbert_adjoint(links[0], spaces[0], spaces[1])
    for i, link in enumerate(links[1:], start=1):
        chained = chained @ _hilbert_adjoint(link, spaces[i], spaces[i + 1])

    assert chained.shape == (dims[0], dims[-1])
    np.testing.assert_allclose(chained, direct, rtol=0.0, atol=1e-13)

    rng = np.random.default_rng(_SEED + 8)
    x = rng.standard_normal(dims[0])
    y = rng.standard_normal(dims[-1])
    assert spaces[-1].inner_product(composite @ x, y) == pytest.approx(
        spaces[0].inner_product(x, direct @ y), rel=1e-12,
    )


@pytest.mark.verifies("bc-response-factored-adjoint")
def test_the_telescoping_survives_wildly_scaled_interior_metrics():
    r"""The cancellation is exact, not merely well-conditioned.

    An interior metric scaled by :math:`10^{-8}` against another scaled by
    :math:`4\times10^{7}` — 15 orders of magnitude apart — must give the same
    adjoint as the unscaled chain. This is what licenses the design statement
    *"the interior metric is free"*: if the cancellation only held for
    comparably-scaled metrics, the choice would be a numerical liability rather
    than a genuine freedom.
    """
    dims = [7, 3, 5, 4]
    links, spaces = _chain(dims)

    composite = links[0]
    for link in links[1:]:
        composite = link @ composite
    direct = _hilbert_adjoint(composite, spaces[0], spaces[-1])

    scaled = [
        spaces[0],
        FunctionSpace(
            name="probe:chain:V1", shape=(dims[1],),
            inner_product_weights=np.asarray(spaces[1].inner_product_weights) * 1e-8,
        ),
        FunctionSpace(
            name="probe:chain:V2", shape=(dims[2],),
            inner_product_weights=np.asarray(spaces[2].inner_product_weights) * 4e7,
        ),
        spaces[-1],
    ]
    chained = _hilbert_adjoint(links[0], scaled[0], scaled[1])
    for i, link in enumerate(links[1:], start=1):
        chained = chained @ _hilbert_adjoint(link, scaled[i], scaled[i + 1])

    np.testing.assert_allclose(chained, direct, rtol=0.0, atol=1e-13)


# ─────────────────────────────────────────────────────────────────────
# The structural symmetry that replaced a false self-adjointness claim
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.verifies("bc-response-factored-adjoint")
def test_R_and_its_adjoint_share_one_form_with_the_halves_exchanged():
    r"""⭐ What actually survives of the Lambertian's "self-adjointness".

    :class:`~orpheus.sn.boundary.angular.AngularAverageOperator` claimed to be
    *self-adjoint under the cosine-weighted inner product* until 2026-08-04.
    That was true of the pre-B3.4a full-face endomorphism and became
    **type-incoherent** when B3.4a narrowed it to :math:`\Gamma_+ \to \Gamma_-`:
    self-adjointness requires domain :math:`=` codomain, and these are spaces
    over disjoint index sets of different size.

    What survives — and is strictly more useful, because it says how to WRITE
    the transpose — is that :math:`R` and :math:`R^*` have the **same form**
    with the two half-traces exchanged: each broadcasts over its codomain a
    contraction of its domain against *that domain's own metric*, normalised.
    Hence ``apply_transpose`` is this operator's own body with the halves
    swapped.

    ⭐ **The symmetry is NOT generic — it holds because the contraction weight
    IS the domain metric**, which is physics: the Lambertian contracts against
    :math:`|\Omega\cdot\hat n|\,w`, the partial-current metric. Algebraically
    :math:`R^*[i,j] = (g_+[i]/\|g_+\|_1)\,g_-[j]/g_+[i]`, and the
    :math:`g_+[i]` cancels **only** under that identification; with an
    independent contraction weight the rows are merely proportional, not equal.
    A first draft of this gate used independent random weights and failed for
    exactly that reason — the failure was the finding.
    """
    gamma_plus, _, gamma_minus = _spaces(1.0)
    contraction, broadcast = _contraction_and_broadcast()
    response = broadcast @ contraction
    adjoint = _hilbert_adjoint(response, gamma_plus, gamma_minus)

    # Not self-adjoint — it cannot even be spelled: the shapes are transposes.
    assert response.shape == (N_IN, N_OUT)
    assert adjoint.shape == (N_OUT, N_IN)

    # Both rank one, and both CONSTANT along the codomain axis: every row of R
    # is the same functional on Γ₊, every row of R* the same one on Γ₋.
    assert np.linalg.matrix_rank(response) == 1
    assert np.linalg.matrix_rank(adjoint) == 1
    # ``assert_allclose`` does NOT broadcast, so tile explicitly rather than
    # comparing against a ``(1, n)`` row (that fails on shape alone, with
    # bit-identical values — a misleading red).
    np.testing.assert_allclose(
        response, np.broadcast_to(response[0], response.shape), atol=1e-15,
    )
    np.testing.assert_allclose(
        adjoint, np.broadcast_to(adjoint[0], adjoint.shape), atol=1e-15,
    )

    # ⭐ The mirror image, stated exactly: R's row is Γ₊'s metric normalised,
    # R*'s row is Γ₋'s metric divided by the SAME normalisation ‖g₊‖₁.
    np.testing.assert_allclose(
        response[0], _METRIC_PLUS / _METRIC_PLUS.sum(), atol=1e-15,
    )
    np.testing.assert_allclose(
        adjoint[0], _METRIC_MINUS / _METRIC_PLUS.sum(), atol=1e-15,
    )
