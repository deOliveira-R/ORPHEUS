r"""Phase 2.5a (#280) — the "one 1-D loop walk" discriminating tests.

The structural payoff of the orientation-parametrized 1-D apply-loop frame:
the forward matvec (:math:`(L+C)\psi`, ``loss_action``) and the adjoint
matvec (:math:`(L+C)^{\mathsf T}\varphi`, ``loss_action_transpose``) on the
SAME 1-D mesh MUST route their per-cell marches through the SAME
:meth:`~orpheus.sn.loss_representation._OneDimScanWalk._loop_walk` frame
over :meth:`~orpheus.sn.loss_representation._OneDimScanWalk._dag_legs`-built
legs — they fork only at the per-cell kernel + endpoint closures, NEVER in
a duplicated leg frame.  The 1-D sibling of ``test_one_octant_walk.py``
(S6.4), minted by the 2.5a carve (gate spec
``a3_solve_transpose_verification.md`` §10.3).

Value-level proof of the relocation lives elsewhere (the FROZEN pre-carve
``walk_matvec_*`` baselines, ``tests/sn/regression/``); these tests pin the
STRUCTURE: one walk observed at runtime for both orientations, orientation
carried by objects/data (never a boolean flag), the reversal law, the
dependency order, and the leg/degenerate trichotomy.

``-O``-safe (vv Mode 8): ``pytest.fail`` only — no bare asserts.
``foundation`` — software-structure invariants (no theory ``:label:``).
"""
from __future__ import annotations

import ast
import inspect
import textwrap

import numpy as np
import pytest

from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import StreamingOperator
from tests.sn.operators.test_g_adjoint_reciprocity import (
    _make_cyl_product,
    _make_slab,
    _make_sphere,
    _random_composite,
)

pytestmark = pytest.mark.foundation


def _product_cylinder(ng: int = 1) -> SNMesh:
    """The degenerate-class cylinder (φ = π/2, 3π/2 ⇒ |μ_x| ≈ 6e-17) —
    the reciprocity file's :func:`_make_cyl_product` builder, mesh only.
    The partition pin below needs a quadrature where all three direction
    classes are non-empty (the ``level_symmetric`` rule has no degenerate
    ordinates)."""
    return _make_cyl_product(ng=ng)[0]


def test_apply_and_apply_transpose_hit_one_loop_walk(monkeypatch):
    """[L0 structural] the forward AND the adjoint matvec exercise the SAME
    ``_OneDimScanWalk._loop_walk`` (the one-loop claim, observed at runtime).

    NEGATIVE PRE-CONDITION (why this is a tripwire, not a tautology): before
    2.5a the leg orchestration was duplicated in lockstep — ``_apply_walk``'s
    per-level mask/DAG-order loop vs ``loss_action_transpose``'s — so a
    drift in one body's leg decomposition was structurally possible and
    invisible to the self-referential removal-form canaries.  Since 2.5a
    both orientations march through the ONE frame: "the adjoint walks the
    SAME DAG, reversed" is a CODE FACT.

    Sphere included deliberately: its legs exercise the mirror-routed pole
    handoff (slab legs are independent), so the spy also witnesses the
    orientation object carrying the Carlson seed routing.
    """
    from orpheus.sn.loss_representation import _OneDimScanWalk

    hits: list[str] = []
    real_walk = _OneDimScanWalk._loop_walk

    def spy(self, *args, **kwargs):
        hits.append("walked")
        return real_walk(self, *args, **kwargs)

    monkeypatch.setattr(_OneDimScanWalk, "_loop_walk", spy)

    for name, build in (("slab", _make_slab), ("sphere", _make_sphere)):
        sn, _sig_t = build(ng=2)
        rng = np.random.default_rng(20260704)
        psi = _random_composite(sn, rng)
        L = StreamingOperator.pose(sn)

        # (1) the FORWARD matvec routes through the shared loop frame:
        hits.clear()
        _ = L.apply(psi)
        if not hits:
            pytest.fail(
                f"[{name}] loss_action did NOT route through "
                "_OneDimScanWalk._loop_walk — the forward matvec has a "
                "private duplicated leg frame again (the 2.5a carve "
                "regressed)."
            )

        # (2) the ADJOINT matvec routes through the SAME loop frame:
        hits.clear()
        _ = L.apply_transpose(psi)
        if not hits:
            pytest.fail(
                f"[{name}] loss_action_transpose did NOT route through "
                "_OneDimScanWalk._loop_walk — the adjoint matvec has a "
                "private duplicated leg frame again (the 2.5a carve "
                "regressed)."
            )


def test_one_dim_walk_is_orientation_object_not_boolean():
    """[L0 structural] the 1-D frame forks on OBJECTS (leg schedule +
    kernel/endpoint closures), NOT a boolean orientation flag
    (coding-elegance Smell #3 — the ``test_one_octant_walk`` rule's
    ORIENTATION sibling, gate spec §10.3).

    A boolean dispatch (``if is_adjoint: reversed(...)``) inside the shared
    frame would re-introduce the twin-path branch the carve eliminates.
    Orientation must be carried by the DATA (``_reverse_traversal`` of the
    leg schedule) and the injected callables — the ``_ApplyOperands`` /
    ``_SolveOperands`` / ``_SweepEmit`` discipline.  Source-inspection via
    the AST's identifiers, so docstrings/comments that NAME the
    anti-pattern don't trip it.  ``-O``-safe.
    """
    from orpheus.sn.loss_representation import (
        _OneDimScanWalk,
        _WalkLeg,
        _reverse_traversal,
    )

    # The direction smells of the octant rule PLUS the orientation smells
    # of the 2.5a fold — this frame hosts solve, apply, AND adjoint doors.
    smells = {
        "is_solve", "is_apply", "is_matvec",
        "is_adjoint", "is_forward", "is_transpose", "is_reverse",
    }
    offenders: set[str] = set()
    for obj in (_OneDimScanWalk, _WalkLeg, _reverse_traversal):
        tree = ast.parse(textwrap.dedent(inspect.getsource(obj)))
        identifiers = {
            node.id for node in ast.walk(tree) if isinstance(node, ast.Name)
        } | {
            node.attr for node in ast.walk(tree)
            if isinstance(node, ast.Attribute)
        } | {
            node.arg for node in ast.walk(tree)
            if isinstance(node, (ast.arg, ast.keyword))
            and node.arg is not None
        }
        offenders |= identifiers & smells
    if offenders:
        pytest.fail(
            f"the 1-D loop-walk frame carries boolean orientation flag(s) "
            f"{sorted(offenders)} — the carve degraded into the "
            "boolean-flag anti-pattern (coding-elegance Smell #3).  "
            "Orientation MUST be carried by the leg schedule + injected "
            "kernel/endpoint objects."
        )


def test_reverse_traversal_is_the_exact_reverse():
    """[L0 structural] ``_reverse_traversal`` is reverse-mode order EXACTLY:
    legs in reverse schedule order, each leg's cell chain reversed, every
    other leg field untouched."""
    from orpheus.sn.loss_representation import (
        _OneDimScanWalk,
        _reverse_traversal,
    )

    for name, sn in (
        ("sphere", _make_sphere(ng=1)[0]),
        ("cyl-product", _product_cylinder()),
    ):
        legs = _OneDimScanWalk(sn)._dag_legs()
        rev = _reverse_traversal(legs)
        if len(rev) != len(legs):
            pytest.fail(f"[{name}] reversal changed the leg count")
        for fwd_leg, rev_leg in zip(reversed(legs), rev):
            if (
                rev_leg.mu_level_idx != fwd_leg.mu_level_idx
                or rev_leg.direction_sign != fwd_leg.direction_sign
                or not np.array_equal(rev_leg.ordinates, fwd_leg.ordinates)
                or not np.array_equal(rev_leg.within, fwd_leg.within)
                or not np.array_equal(rev_leg.abs_mu, fwd_leg.abs_mu)
            ):
                pytest.fail(
                    f"[{name}] reversal permuted leg identity/bundles "
                    f"(level {fwd_leg.mu_level_idx}, sign "
                    f"{fwd_leg.direction_sign:+d})"
                )
            if rev_leg.cells != fwd_leg.cells[::-1]:
                pytest.fail(
                    f"[{name}] leg (level {fwd_leg.mu_level_idx}, sign "
                    f"{fwd_leg.direction_sign:+d}): reversed cells "
                    f"{rev_leg.cells} != {fwd_leg.cells[::-1]}"
                )


def test_dag_legs_dependency_order_inward_before_outward():
    """[L0 structural] ``_dag_legs`` is DEPENDENCY-ordered: every −1
    (inward) leg precedes every +1 (outward) leg — the pole continuation
    ψ(0,+μ) = ψ(0,−μ) (ERR-058a) makes inward legs the predecessors — and
    μ-levels ascend within each sign class."""
    from orpheus.sn.loss_representation import _OneDimScanWalk

    for name, sn in (
        ("slab", _make_slab(ng=1)[0]),
        ("sphere", _make_sphere(ng=1)[0]),
        ("cyl-product", _product_cylinder()),
    ):
        legs = _OneDimScanWalk(sn)._dag_legs()
        signs = [leg.direction_sign for leg in legs]
        n_in = signs.count(-1)
        if signs != [-1] * n_in + [+1] * (len(signs) - n_in):
            pytest.fail(
                f"[{name}] leg schedule {signs} is not "
                "[all −1 legs, then all +1 legs]"
            )
        for sign in (-1, +1):
            levels = [
                leg.mu_level_idx for leg in legs
                if leg.direction_sign == sign
            ]
            if levels != sorted(levels):
                pytest.fail(
                    f"[{name}] sign {sign:+d} μ-levels {levels} not "
                    "ascending"
                )


def test_legs_and_degenerates_partition_the_quadrature():
    """[L0 structural] the μ-direction trichotomy PARTITIONS the ordinate
    set: leg ordinates (both signs) and the degenerate pure-azimuthal set
    are disjoint and jointly exhaustive.

    Run on the product-quadrature cylinder where the degenerate class is
    NON-EMPTY (φ = π/2, 3π/2 → |μ_x| ≈ 6e-17) — a two-threshold drift
    between the leg masks and the degenerate mask would silently drop or
    double-count exactly these ordinates (the single-source
    ``_MU_DIRECTION_EPS`` contract).
    """
    from orpheus.sn.loss_representation import _OneDimScanWalk

    sn = _product_cylinder()
    walk = _OneDimScanWalk(sn)
    leg_ordinates: list[int] = []
    for leg in walk._dag_legs():
        leg_ordinates.extend(int(n) for n in leg.ordinates)
    global_deg, deg_level, deg_within = walk._degenerate_positions()
    degenerates = [int(n) for n in global_deg]

    if not degenerates:
        pytest.fail(
            "product(n_mu=2, n_phi=4) produced NO degenerate ordinates — "
            "the partition pin needs a non-empty degenerate class "
            "(quadrature layout changed?)"
        )
    if len(deg_level) != len(degenerates) or len(deg_within) != len(degenerates):
        pytest.fail(
            "_degenerate_positions returned ragged (level, within) "
            "resolutions — every degenerate ordinate must resolve to "
            "exactly one level position"
        )
    overlap = set(leg_ordinates) & set(degenerates)
    if overlap:
        pytest.fail(
            f"ordinates {sorted(overlap)} appear on a leg AND in the "
            "degenerate set — the trichotomy double-counts"
        )
    covered = set(leg_ordinates) | set(degenerates)
    if covered != set(range(sn.quad.N)):
        pytest.fail(
            f"legs ∪ degenerates covers {sorted(covered)} != all "
            f"{sn.quad.N} ordinates — the trichotomy drops ordinates"
        )
