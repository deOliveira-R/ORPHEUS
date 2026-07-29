r"""The trace-digraph acyclicity certificate — when does the WHOLE within-group
problem sweep in one pass?

System under test: :mod:`orpheus.derivations.discrete.sn.sweep_acyclicity`, the
algebra of record for the boundary half of the sweep-triangularity question.

Why this file exists
====================

The production certificate ``triu(PᵀMP, 1) == 0``
(``test_assembly_mode.py::test_g2_walk_order_triangularity_is_exact``) pins the
**bulk** claim: for a fixed ordinate, ``L + C`` is triangular in walk order.
It cannot see the **trace** claim, because ``B`` is a sibling operator —
structurally outside ``L + C`` — so no boundary edge ever enters that matrix.

That left a load-bearing proposition asserted in prose and tested nowhere:
whether a reflecting face closes a cycle.  Two theory pages disagreed about it
in opposite directions.  These gates settle it by construction.

Marks: ``foundation`` — this is a software/structural invariant of a discrete
operator's dependency graph, not a physics claim against a reference, so it
carries no ``verifies(...)`` (the verifies⊥level doctrine).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.discrete.sn.sweep_acyclicity import (
    build_slab_trace_digraph,
    derive_slab_trace_acyclicity,
)
from orpheus.numerics.quadrature import Quadrature


def _slab(n_ordinates: int = 4):
    q = Quadrature.gauss_legendre(n_ordinates=n_ordinates)
    return np.asarray(q.mu_x, dtype=float), q.reflection_index(0)


# ── the headline: acyclicity is per-CONFIGURATION, not per-boundary-KIND ──


@pytest.mark.foundation
@pytest.mark.parametrize(
    ("bc_left", "bc_right", "expect_acyclic"),
    [
        ("vacuum", "vacuum", True),
        ("reflective", "vacuum", True),      # ← the folklore refutation
        ("vacuum", "reflective", True),      # ← and its mirror image
        ("reflective", "reflective", False),  # ← the closed LOOP
        ("periodic", "vacuum", False),       # wraps to the opposite face
    ],
)
def test_single_reflection_is_acyclic_a_closed_loop_is_not(
    bc_left, bc_right, expect_acyclic,
):
    r"""[L11 paired] A *single* reflecting face does NOT close a cycle; a
    closed loop of reflections does.

    This is the discriminating gate for the proposition "a reflective boundary
    creates a sweep cycle".  It is FALSE as stated — POSITIVE rows (one
    reflecting face, either side) are acyclic; the NEGATIVE row (both faces)
    is not.  The criterion is therefore a strongly-connected-component
    decomposition of the trace digraph, never a boolean on the boundary kind.
    """
    mu, refl = _slab()
    G = build_slab_trace_digraph(mu, refl, bc_left=bc_left, bc_right=bc_right)
    if G.is_acyclic is not expect_acyclic:
        pytest.fail(
            f"{bc_left}|{bc_right}: expected acyclic={expect_acyclic}, got "
            f"{G.is_acyclic}. Cycles: {G.describe_cycles()}"
        )


@pytest.mark.foundation
def test_acyclic_configurations_certify_a_lower_triangular_matrix():
    r"""Acyclic ⟹ a topological order exists ⟹ the dependency matrix permuted
    to it is EXACTLY lower-triangular — the same ``triu(M, 1) == 0`` test the
    production bulk certificate uses, on the trace digraph.

    Keeping ONE convention (row = equation, col = unknown it reads) across the
    bulk gate and this one is deliberate: the two certificates then compose
    into a single statement about the whole within-group operator.
    """
    mu, refl = _slab()
    for (bl, br), G in derive_slab_trace_acyclicity(mu, refl).items():
        if not G.is_acyclic:
            continue
        order = G.topological_order()
        if order is None:
            pytest.fail(f"{bl}|{br}: acyclic but no topological order")
        M = G.permuted_dependency_matrix(order)
        np.testing.assert_array_equal(
            np.triu(M, k=1), 0,
            err_msg=(
                f"{bl}|{br}: acyclic yet NOT lower-triangular in its own "
                "topological order — the SCC test and the matrix test "
                "disagree, so one of them is wrong."
            ),
        )


@pytest.mark.foundation
def test_cyclic_configuration_names_the_edges_that_close_the_loop():
    r"""The certificate must be ACTIONABLE, not just a boolean: on the
    reflective|reflective slab it reports the SCCs and attributes each to the
    specific boundary edges that closed it.

    For S4 the mirror pairing is ``n ↔ N−1−n``, so the 4 streaming ordinates
    close into 2 disjoint 4-slot loops (one per mirror pair) — each an
    ``inflow → outflow → mirror-inflow → mirror-outflow → back`` cycle.
    """
    mu, refl = _slab()
    G = build_slab_trace_digraph(
        mu, refl, bc_left="reflective", bc_right="reflective",
    )
    comps = G.cyclic_components()
    if len(comps) != 2:
        pytest.fail(f"expected 2 mirror-pair SCCs for S4, got {len(comps)}")
    for comp in comps:
        if len(comp) != 4:
            pytest.fail(f"expected 4-slot loops, got {len(comp)}: {comp}")
    described = G.describe_cycles()
    if not all("closed by boundary edges" in d for d in described):
        pytest.fail(f"cycles not attributed to their closing edges: {described}")


# ── the second finding: acyclic ≠ triangular in an ARBITRARY order ──


@pytest.mark.foundation
def test_a_single_reflection_is_one_pass_only_in_the_matching_direction():
    r"""Acyclicity says SOME order works; it does not say a GIVEN order does.

    A left-reflecting slab is one-pass in the μ<0-first order (the μ<0 sweep
    reaches the left face, and the mirror feeds strictly downstream).  A
    RIGHT-reflecting slab is equally acyclic but NOT triangular in that same
    order — it needs μ>0 first.  Both are one-pass; the sweep direction has to
    match the reflecting face.

    This is the distinction that makes "is it triangular?" an ill-posed
    question about an operator: it is a question about (operator, ORDER).
    """
    mu, refl = _slab()
    fixed = lambda G: int(  # noqa: E731 - local alias for readability
        np.triu(G.permuted_dependency_matrix(G.mu_negative_first_order()), k=1).sum()
    )
    left = build_slab_trace_digraph(mu, refl, bc_left="reflective", bc_right="vacuum")
    right = build_slab_trace_digraph(mu, refl, bc_left="vacuum", bc_right="reflective")

    # BOTH are acyclic — a one-pass order exists for each.
    if not (left.is_acyclic and right.is_acyclic):
        pytest.fail("both single-reflection slabs must be acyclic")
    # But the FIXED μ<0-first order is triangular for only one of them.
    if fixed(left) != 0:
        pytest.fail("left-reflecting slab must be triangular in μ<0-first order")
    if fixed(right) == 0:
        pytest.fail(
            "right-reflecting slab is triangular in μ<0-first order — the "
            "order-sensitivity control has no teeth"
        )


# ── teeth: the certificate must RED on a broken digraph ──


@pytest.mark.foundation
def test_teeth_dropping_the_boundary_edge_falsely_certifies_acyclic():
    r"""Mutation: build the reflective|reflective slab but suppress the
    boundary edges (the ``vacuum|vacuum`` digraph stands in).  The SCC test
    MUST flip to acyclic — proving the cyclic verdict is driven by the
    boundary edges and not by an artefact of the sweep edges.
    """
    mu, refl = _slab()
    cyclic = build_slab_trace_digraph(
        mu, refl, bc_left="reflective", bc_right="reflective",
    )
    no_bnd = build_slab_trace_digraph(
        mu, refl, bc_left="vacuum", bc_right="vacuum",
    )
    if cyclic.is_acyclic:
        pytest.fail("reflective|reflective must be cyclic")
    if not no_bnd.is_acyclic:
        pytest.fail("with no boundary edges the digraph must be acyclic")
    if len(no_bnd.boundary_edges) != 0:
        pytest.fail("vacuum|vacuum must contribute zero boundary edges")
    if len(cyclic.boundary_edges) == 0:
        pytest.fail("reflective|reflective must contribute boundary edges")


@pytest.mark.foundation
def test_white_boundary_couples_every_ordinate_on_the_face():
    r"""``white`` is not ``reflective`` with a different partner map: it
    couples EVERY outflow on the face to EVERY inflow on it (rank-1 in angle).

    Consequence, and the reason the production Gauss-Seidel schedule excludes
    white from its octant split: a single white face already couples all
    ordinates, so an octant-ordered G-S degenerates to Jacobi.  A white face
    opposite a vacuum is still acyclic (the coupling is all-to-all but still
    strictly downstream); white on BOTH faces is not.
    """
    mu, refl = _slab()
    n_stream = int(np.count_nonzero(mu))
    one = build_slab_trace_digraph(mu, refl, bc_left="white", bc_right="vacuum")
    # every outflow on the left couples to every inflow on the left
    n_out_left = sum(
        1 for s in one.slots if s.face == "left" and s.role == "outflow"
    )
    n_in_left = sum(
        1 for s in one.slots if s.face == "left" and s.role == "inflow"
    )
    if len(one.boundary_edges) != n_out_left * n_in_left:
        pytest.fail(
            f"white must be all-to-all on the face: expected "
            f"{n_out_left * n_in_left} edges, got {len(one.boundary_edges)}"
        )
    if not one.is_acyclic:
        pytest.fail("a single white face is still acyclic")
    both = build_slab_trace_digraph(mu, refl, bc_left="white", bc_right="white")
    if both.is_acyclic:
        pytest.fail("white on both faces must close a cycle")
    if n_stream != 4:
        pytest.fail(f"fixture drift: expected 4 streaming ordinates, got {n_stream}")
