r"""Tests for :class:`orpheus.numerics.operator.TraceRestrictionOperator`.

Campaign phase **B3.1** (``.claude/plans/b3_domain_narrowing_crosswalk.md``).
This operator is the trace map :math:`\gamma_\pm` of the affine boundary form
:math:`\gamma_-\psi = R\,G\,\gamma_+\psi + q` — the map the theory page has
always named and the codebase had spelled three ways and typed as none.

Tagged ``@pytest.mark.foundation``: these assert an algebraic contract on a
numerics primitive, not a discretization claim, so they carry no
``verifies(...)`` equation label (the verifies⊥level doctrine).

**The first half tests the DEFINING LAWS, not the usage.** A type that carries
a mathematical concept ships a test of the identities that make it that
concept — here the restriction/scatter adjoint pair and the projector it
generates — so a future edit that keeps every call site working while breaking
:math:`\gamma\iota = I` still reds.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.operator import TraceRestrictionOperator


pytestmark = [pytest.mark.foundation]


def _as_matrix(op: TraceRestrictionOperator) -> np.ndarray:
    """Materialise the restriction by applying it to the full identity."""
    return op.apply(np.eye(op.n_total))


# ═══════════════════════════════════════════════════════════════════════
# The defining laws
# ═══════════════════════════════════════════════════════════════════════


def test_gather_selects_exactly_the_indexed_rows() -> None:
    """:math:`(\\gamma_S x)_j = x_{S(j)}` — the definition itself."""
    x = np.arange(24.0).reshape(6, 4)
    gamma = TraceRestrictionOperator(np.array([1, 3, 4]), n_total=6)
    assert np.array_equal(gamma.apply(x), x[[1, 3, 4], :])


def test_scatter_is_the_exact_matrix_transpose_of_the_gather() -> None:
    r""":math:`\iota_S = \gamma_S^{\mathsf T}` — materialised, not asserted.

    The adjoint claim is what licenses using ``apply_transpose`` wherever the
    algebra calls for :math:`\gamma^{\mathsf T}`, so it is checked against the
    dense transpose rather than against another call of the same code.
    """
    gamma = TraceRestrictionOperator(np.array([0, 2, 5]), n_total=6)
    G = _as_matrix(gamma)                       # (3, 6)
    I_restricted = np.eye(gamma.n_restricted)   # (3, 3)
    scattered = gamma.apply_transpose(I_restricted)   # (6, 3)
    assert np.array_equal(scattered, G.T)


def test_gather_after_scatter_is_the_identity_on_the_restricted_space() -> None:
    r""":math:`\gamma_S \circ \iota_S = I_m` — the pair is a section/retraction."""
    gamma = TraceRestrictionOperator(np.array([1, 2, 5]), n_total=7)
    y = np.array([[3.0, -1.0], [0.5, 2.0], [7.0, 7.0]])
    assert np.array_equal(gamma.apply(gamma.apply_transpose(y)), y)


def test_scatter_after_gather_is_an_idempotent_symmetric_projector() -> None:
    r""":math:`P_S = \iota_S \circ \gamma_S` satisfies :math:`P^2 = P = P^{\mathsf T}`.

    This is the map the boundary subsystem had written twice by hand (a
    slice-write and a dense diagonal multiply); pinning its two defining
    properties here is what lets both spellings retire onto one composition.
    """
    gamma = TraceRestrictionOperator(np.array([0, 3]), n_total=5)
    P = gamma.apply_transpose(_as_matrix(gamma))      # (5, 5)
    assert np.array_equal(P @ P, P), "projector is not idempotent"
    assert np.array_equal(P, P.T), "projector is not symmetric"
    # ...and it is the diagonal 0/1 mask on the selected rows.
    expected = np.diag([1.0, 0.0, 0.0, 1.0, 0.0])
    assert np.array_equal(P, expected)


def test_disjoint_restrictions_annihilate_each_other() -> None:
    r""":math:`\gamma_T \circ \iota_S = 0` for :math:`S \cap T = \emptyset`.

    The composition the review found stated as a curiosity —
    ``P_in ∘ P_out = 0`` — is this, and it holds for the structural reason
    that two disjoint index sets have nothing to hand each other.
    """
    inflow = TraceRestrictionOperator(np.array([0, 1]), n_total=4)
    outflow = TraceRestrictionOperator(np.array([2, 3]), n_total=4)
    y = np.array([[1.0], [2.0]])
    assert np.array_equal(
        outflow.apply(inflow.apply_transpose(y)), np.zeros((2, 1)),
    )


def test_a_three_way_partition_resolves_the_identity() -> None:
    r""":math:`\sum_S \iota_S\gamma_S = I` over a partition — and the partition
    at an SN face is **three**-way, not two.

    Inflow, outflow AND tangential: on a real cylinder under
    ``product(2, 4)``, four of eight ordinates are tangential. So
    "not inflow" is not "outflow", and a sum over only the first two halves
    falls short of the identity by exactly the tangential rows — which is why
    the re-posed boundary gates assert on ``complement(inflow)`` rather than
    on ``outflow``.
    """
    n = 8
    inflow = TraceRestrictionOperator(np.array([2, 6]), n_total=n)
    outflow = TraceRestrictionOperator(np.array([0, 4]), n_total=n)
    tangent = TraceRestrictionOperator(np.array([1, 3, 5, 7]), n_total=n)

    def projector(op: TraceRestrictionOperator) -> np.ndarray:
        return op.apply_transpose(_as_matrix(op))

    two_way = projector(inflow) + projector(outflow)
    assert not np.array_equal(two_way, np.eye(n)), (
        "inflow + outflow reproduced the identity — the fixture lost its "
        "tangential rows, and with them the whole point of this test"
    )
    three_way = two_way + projector(tangent)
    assert np.array_equal(three_way, np.eye(n))


# ``to_local`` and its gates (the round-trip, the N8 non-prefix trap, the
# crossed-set refusal) moved to the half-trace SPACE at G6.5 — the remap is
# a fact about the subspace's embedding, which the space owns. See
# ``tests/numerics/test_angular_face_trace_space.py``, section "The
# local↔global map".


def test_composes_along_a_trailing_axis_leaving_others_untouched() -> None:
    """``axis`` selects the gathered axis; the rest broadcast."""
    x = np.arange(2 * 5 * 3.0).reshape(2, 5, 3)
    gamma = TraceRestrictionOperator(np.array([1, 4]), n_total=5, axis=1)
    assert np.array_equal(gamma.apply(x), x[:, [1, 4], :])


def test_declares_adjointable_and_structurally_non_invertible() -> None:
    """A restriction discards rows: adjointable, never invertible."""
    gamma = TraceRestrictionOperator(np.array([0, 2]), n_total=4)
    assert gamma.is_adjointable is True
    assert gamma.is_invertible is False
    assert not hasattr(gamma, "inverse")


def test_neither_direction_aliases_its_input() -> None:
    """Both directions return FRESH storage — the caller's array is safe.

    Migrated at B3.3 from ``IncomingOrdinateMaskTensor``'s suite, which is
    the only one of its thirteen tests asserting a claim this battery did
    not already make. It is not idle defensiveness: ``§17.6`` recorded
    :class:`~orpheus.numerics.operator.IdentityOperator` returning its
    input **by reference**, so "an operator hands back fresh storage" is a
    real per-operator distinction in this package rather than a universal.

    Asserted both ways round — a mutation of the *output* must not reach
    back into the input, which `assert_array_equal` on the input alone
    would miss if the two shared a buffer.
    """
    gamma = TraceRestrictionOperator(np.array([0, 2]), n_total=4)

    x = np.array([1.0, 2.0, 3.0, 4.0])
    gathered = gamma.apply(x)
    assert not np.shares_memory(gathered, x)
    gathered[0] = -99.0
    assert np.array_equal(x, np.array([1.0, 2.0, 3.0, 4.0]))

    y = np.array([7.0, 8.0])
    scattered = gamma.apply_transpose(y)
    assert not np.shares_memory(scattered, y)
    scattered[0] = -99.0
    assert np.array_equal(y, np.array([7.0, 8.0]))


# ═══════════════════════════════════════════════════════════════════════
# The guards — each one prevents a measured failure mode
# ═══════════════════════════════════════════════════════════════════════


def test_unsorted_indices_are_refused() -> None:
    """Sortedness is what makes ``to_local``'s ``searchsorted`` correct."""
    with pytest.raises(ValueError, match="SORTED"):
        TraceRestrictionOperator(np.array([3, 1, 2]), n_total=5)


def test_duplicate_indices_are_refused() -> None:
    """A repeated row breaks ``γ ∘ ι = I``: the scatter drops all but one."""
    with pytest.raises(ValueError, match="unique"):
        TraceRestrictionOperator(np.array([1, 2, 2]), n_total=5)


@pytest.mark.parametrize("bad", [np.array([-1, 2]), np.array([0, 5])])
def test_out_of_range_indices_are_refused(bad) -> None:
    with pytest.raises(ValueError, match=r"lie in"):
        TraceRestrictionOperator(bad, n_total=5)


def test_non_1d_indices_are_refused() -> None:
    with pytest.raises(ValueError, match="1-D"):
        TraceRestrictionOperator(np.zeros((2, 2), dtype=np.intp), n_total=5)


def test_apply_refuses_an_input_of_the_wrong_length() -> None:
    """The guard a hand-rolled reduced permutation does not have.

    Measured: a reduced ``PermutationOperator`` fed a full-length input
    returns a same-shaped array of WRONG VALUES with no raise. That silent
    mode is what this refusal closes — and it is the negative test the
    domain-narrowing gates (RG-3b) are built on.
    """
    gamma = TraceRestrictionOperator(np.array([0, 2]), n_total=4)
    with pytest.raises(ValueError, match="DOMAIN is the full space"):
        gamma.apply(np.ones((2, 3)))      # the RESTRICTED length, wrongly fed in


def test_apply_transpose_refuses_an_input_of_the_wrong_length() -> None:
    gamma = TraceRestrictionOperator(np.array([0, 2]), n_total=4)
    with pytest.raises(ValueError, match="DOMAIN is the restricted space"):
        gamma.apply_transpose(np.ones((4, 3)))   # the FULL length, wrongly fed in


# The crossed-index-set refusal moved with ``to_local`` to the half-trace
# space's battery (G6.5).
