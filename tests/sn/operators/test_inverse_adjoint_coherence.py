r"""The #280 Phase 2.5c inverse-adjoint wiring — the swap law ``A.H.inverse() ≡
A.inverse().H`` on the SN loss ``(L+C)``.

2.5c closes the adjoint-inverse of the operator algebra: the inverse of the
Hilbert adjoint IS the adjoint of the inverse,

.. math::   (A^{*})^{-1} = (A^{-1})^{*},

an OBJECT IDENTITY (not a computed numerical equivalence).  The landed wiring
(Part A, ``orpheus/numerics/operator.py`` + ``orpheus/sn/operators/sweep_operator.py``):

* ``A = _loss(sn)`` is the ``(L+C)`` :class:`~orpheus.sn.operators.streaming.InvertibleOperator`.
* ``A.inverse()`` → :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`;
  NEW ``SweepOperator.apply_transpose(b) = inner.solve_transpose(b)`` (the 2.5b
  reverse-scan ``(A⁻¹)ᵀ = (Aᵀ)⁻¹``, the plain EUCLIDEAN transpose-solve) and
  ``SweepOperator.is_adjointable`` flips ``True`` over the ``InvertibleOperator`` arm.
* ``A.H`` → :class:`~orpheus.numerics.operator._AdjointOperator`; NEW
  ``_AdjointOperator.inverse() = inner.inverse().H`` (the swap law) and
  ``_AdjointOperator.is_invertible = invertible(inner) and adjointable(inner.inverse())``.
* ``_AdjointOperator.apply`` (UNCHANGED) already does ``G⁺·inner.apply_transpose(G·y)``,
  so the metric adjoint-solve ``A.H.inverse().apply(b) = G⁺·A.solve_transpose(G·b)``
  falls out FOR FREE — no new metric code enters the sweep.

The gates (spec ``a3_solve_transpose_verification.md`` §13 / Deliverable 3 §5):

* **G1 — forward-matvec G-reciprocity pin.** ``⟨A.apply(ψ), x⟩_G = ⟨ψ, b⟩_G``
  for ``x = A.H.inverse().apply(b)``.  The verification arithmetic exercises
  ONLY the FORWARD ``A.apply`` + the metric inner product — it NEVER calls
  ``apply_transpose`` / ``solve_transpose`` — so it is STRUCTURALLY INDEPENDENT
  of the reverse-scan (the SUT ``x`` is the only reverse-scan invocation).
  ``⟨Aψ, (A*)⁻¹b⟩_G = ⟨ψ, A*(A*)⁻¹b⟩_G = ⟨ψ, b⟩_G`` holds iff ``x`` genuinely
  solves ``A* x = b``.
* **G2 — swap-law value (object identity).** ``A.H.inverse().apply(b) ≡
  A.inverse().H.apply(b)`` — BIT-IDENTICAL (both sides construct and run the
  SAME ``A.inverse().H`` object graph on the same input).
* **G3 — adjoint-inverse round-trip.** ``A.H.apply(A.H.inverse().apply(b)) ≈ b``
  = ``A*(A*)⁻¹ b = b`` on the source-carried subspace.
* **G4 — Mode-11 wrap sentinel.** The reverse scan EXECUTES (``InvertibleOperator.
  solve_transpose`` is entered ``> 0`` times) and the FORWARD ``solve`` is NOT
  touched (counter ``== 0``) when ``A.H.inverse().apply(b)`` runs — the gold
  standard in-process wrap counter (a green gate that routed around the new line
  would leave the counter at 0 and red).
* **mutations** M-ADJ-swap (drop the ``.H`` in ``_AdjointOperator.inverse``) reds
  G1 + G2; M-ADJ-metric (skip the ``G⁺``/``G`` wrap in ``_AdjointOperator.apply``)
  reds G1.
* **predicate flips** — ``SweepOperator.is_adjointable`` True (over
  ``InvertibleOperator``); the sibling wrap-delegates ``InverseOperator`` /
  ``GreenOperator`` STAY False (the flip is SURGICAL).
* **static pin** — ``assert_type(A.H.inverse(), LinearOperator[FullField, FullField])``.

Why ``b`` is **bulk-only** (source-carried) in G1/G3.  ``solve``/``apply`` (and
their transposes) are genuine inverses only on the source-carried subspace: the
solve COMPUTES the boundary/seed outflow slots that apply treats as FREE DOFs
(the #284 subspace).  A fully-random ``b`` carries components in the metric null
space (tangential trace ``|Ω·n| = 0`` ⊕ the all-zero-ghost-metric seed) and
outside the range of ``A*``, so the reciprocity/round-trip identities are
well-posed ONLY for a bulk-only ``b`` (empirically clean to ~1e-16 on
slab/sphere/cyl; a random full ``b`` sits at O(1e-1)).  ``ψ`` stays FULL
(bulk+trace random) so the operator's trace coupling and the trace metric ARE
exercised in the pairing.  The seed block stays present-but-ZERO on carrying
meshes (sphere) — the Mode-12 G-reciprocity regime.

CRITICAL DESIGN DEVIATION (M-ADJ-metric — flagged to the main agent).  The spec
predicted the M-ADJ-metric mutation stays GREEN on the slab ("slab's metric is
trivial/uniform") and reds only sphere/cyl — a ``.H``≠Euclidean *discriminator*.
This is EMPIRICALLY FALSE: the slab metric ``G = V·w_n`` (bulk) ⊕ ``|Ω·n|·w_n``
(trace) is non-trivial (non-uniform mesh + the bulk-vs-trace weight mismatch that
streaming couples), so ``A* = G⁻¹AᵀG ≠ Aᵀ`` on the slab TOO and the mutation reds
G1 on ALL geometries (slab 0.33 / sphere 0.19 / cyl 0.13).  The ``.H``≠Euclidean
claim is TRUE everywhere (the metric wrap is load-bearing on every geometry —
STRONGER coverage than the spec's discriminator); it is simply NOT a
slab-vs-curvilinear split.  ``test_mutation_adj_metric_reds_gate1`` therefore
asserts RED on all three geometries.

``@pytest.mark.foundation`` — algebraic/object identities (no theory ``:label:``).
``-O``-safe (vv Mode 8): every value gate is ``np.testing.assert_*`` / ``pytest.fail``
(a bare ``assert`` is stripped to a NO-OP under the canonical ``python -O``).
Mutations are IN-PROCESS ``monkeypatch`` (auto-reverted) — NEVER ``git checkout``
an uncommitted file (lessons L28 / ``.claude/rules/process-discipline.md``).
"""
from __future__ import annotations

from typing import TYPE_CHECKING, assert_type

import numpy as np
import pytest

from orpheus.numerics.green_operator import GreenOperator
from orpheus.numerics.operator import (
    InverseOperator,
    LinearOperator,
    _AdjointOperator,
    invertible,
)
from orpheus.sn.operators.boundary import SNBoundaryOperator
from orpheus.sn.operators.streaming import InvertibleOperator
from orpheus.sn.operators.sweep_operator import SweepOperator
from orpheus.transport.full_field import FullField

# Reuse the landed mesh/operator/field helpers — do NOT re-derive (2.5b keystone
# file; ``_slab`` vacuum het 2G, ``_sphere`` GL seed-carrying, ``_cyl_product``
# degenerate-ordinate + live seed-fold).
from tests.sn.operators.test_loss_transpose_solve import (
    _cyl_product,
    _fresh,
    _loss,
    _sphere,
    _slab,
)

pytestmark = pytest.mark.foundation

#: Reciprocity / round-trip PASS tolerance (clean values ~1e-16; mutations red
#: at O(1e-1) — a decade of margin both ways).  NOT ``inner_tol``-scaled: the
#: reverse-scan is DIRECT (single-pass slab/cyl; the sphere seed is single-call
#: exact post-2.5d — no iterate-threading needed).
_RTOL = 1e-10
#: A mutation is CAUGHT when the residual exceeds this (all mutations red at
#: ≥ 0.13; the clean gates sit at ≤ ~1e-15).
_CAUGHT = 1e-6

# Step 6 (R-6.2): the block swap-law rows run on the NON-carrying meshes;
# on a carrying mesh the bare ``(L+C)`` is the honest ray-decoupled (A,A)
# block and IS adjointable (the predicate row below), while the JOINT swap
# law lives on the coupled sibling (the ``TestCoupledSwapLaw`` arm — the
# same G1–G4 claims on ``M``).
_MESHES = {"slab": _slab, "cyl_product": _cyl_product}


def _rand_bulk(sn, seed: int) -> FullField:
    r"""A bulk-only source-carried RHS ``b``: random bulk, ZERO boundary/seed.

    The reverse-solve is a genuine inverse of the reverse-matvec only on the
    source-carried subspace (#284: ``solve`` COMPUTES the boundary/seed outflow
    slots that ``apply`` treats as free DOFs), so the adjoint-inverse
    reciprocity/round-trip are well-posed for a bulk-only ``b``.
    """
    f = _fresh(sn)
    f.interior.values[:] = np.random.default_rng(seed).standard_normal(f.interior.values.shape)
    return f


def _rand_full(sn, seed: int) -> FullField:
    r"""A full random test vector ``ψ``: random bulk AND boundary, ZERO seed.

    Full so ``A.apply(ψ)`` exercises the streaming trace coupling and the
    ``|Ω·n|·w_n`` trace metric enters the pairing.  The seed block stays
    present-but-zero (the Mode-12 G-reciprocity regime — the all-zero-ghost-metric
    seed coupling is un-balanceable, so the law lives on the bulk⊕trace subspace).
    """
    f = _fresh(sn)
    rng = np.random.default_rng(seed)
    f.interior.values[:] = rng.standard_normal(f.interior.values.shape)
    for face in sn.angular_trace.layout.faces:
        f.boundary.face_view(face)[:] = rng.standard_normal(
            f.boundary.face_view(face).shape
        )
    return f


def _reciprocity_rel(A, psi, x, b, space) -> float:
    r"""``|⟨A.apply(ψ), x⟩_G − ⟨ψ, b⟩_G| / (…)`` — the G1 residual.

    Forward ``A.apply`` + the production ``G`` inner product ONLY; never the
    transpose path (structural independence of the reverse-scan).
    """
    lhs = space.inner_product(A.apply(psi), x)
    rhs = space.inner_product(psi, b)
    return abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)


def _adjoint_inverse(
    A: LinearOperator[FullField, FullField],
) -> LinearOperator[FullField, FullField]:
    r"""``A.H.inverse()`` through the Design-C ``invertible()`` narrowing bridge.

    ``.inverse()`` lives on :class:`~orpheus.numerics.operator.SupportsInverse`,
    not on the base ``LinearOperator`` that ``.H`` returns — so the runtime
    predicate ``invertible(A.H)`` (reading the NEW ``_AdjointOperator.is_invertible``,
    2.5c) IS the static permission to call it (``TypeGuard``, positive branch),
    and doubles as a live pin that ``A.H`` advertises invertibility.
    """
    ah = A.H
    if invertible(ah):
        return ah.inverse()
    pytest.fail(
        f"A.H must advertise is_invertible for the 2.5c swap law; got "
        f"is_invertible={getattr(ah, 'is_invertible', None)}"
    )


# ── G1 — forward-matvec G-reciprocity pin ──────────────────────────────


@pytest.mark.parametrize("geom", list(_MESHES))
def test_gate1_forward_matvec_g_reciprocity(geom):
    r"""``⟨A.apply(ψ), x⟩_G = ⟨ψ, b⟩_G`` for ``x = A.H.inverse().apply(b)``.

    The forward-matvec pin: ``x`` (the adjoint-inverse-solve output) paired
    against the FORWARD operator under ``G`` reproduces the trivial pairing
    ``⟨ψ, b⟩_G`` — provable ONLY if ``x`` genuinely solves ``A* x = b``.  Never
    calls the transpose path in the verification arithmetic, so it is
    structurally independent of the reverse-scan.
    """
    sn = _MESHES[geom]()
    A = _loss(sn)
    b = _rand_bulk(sn, 5)
    psi = _rand_full(sn, 4)
    x = _adjoint_inverse(A).apply(b)
    rel = _reciprocity_rel(A, psi, x, b, sn.full_field_space)
    if not rel < _RTOL:
        pytest.fail(
            f"[{geom}] forward-matvec G-reciprocity broken: rel={rel:.2e} "
            f"(≥ {_RTOL:.0e}) — A.H.inverse().apply(b) does not solve A* x = b"
        )


# ── G2 — swap-law value (object identity, bit-identical) ───────────────


@pytest.mark.parametrize("geom", list(_MESHES))
def test_gate2_swap_law_bit_identical(geom):
    r"""``A.H.inverse().apply(b) ≡ A.inverse().H.apply(b)`` — BIT-IDENTICAL.

    The swap law ``(A^*)^{-1} = (A^{-1})^*`` is an OBJECT IDENTITY:
    ``A.H.inverse()`` returns literally ``A.inverse().H``, so both spellings
    build and run the SAME object graph on the same ``b``.  Asserted at 0 ULP
    (``assert_array_equal``) — stronger than the spec's rtol=1e-12 allowance,
    because it is an identity of the algebra, not a numerical coincidence.
    """
    sn = _MESHES[geom]()
    A = _loss(sn)
    b = _rand_bulk(sn, 5)
    lhs = _adjoint_inverse(A).apply(b)
    rhs = A.inverse().H.apply(b)
    np.testing.assert_array_equal(
        np.asarray(lhs.interior.values), np.asarray(rhs.interior.values),
        err_msg=f"[{geom}] swap-law bulk not bit-identical",
    )
    np.testing.assert_array_equal(
        np.asarray(lhs.boundary.values), np.asarray(rhs.boundary.values),
        err_msg=f"[{geom}] swap-law boundary not bit-identical",
    )


# ── G3 — adjoint-inverse round-trip ────────────────────────────────────


@pytest.mark.parametrize("geom", list(_MESHES))
def test_gate3_adjoint_inverse_round_trip(geom):
    r"""``A.H.apply(A.H.inverse().apply(b)) ≈ b`` — ``A*(A*)⁻¹ = I``.

    The metric-conjugated inverse round-trip on the source-carried (bulk)
    subspace: routing ``.H.inverse()`` to the forward ``solve`` (no reverse)
    breaks this O(1) for the non-symmetric het+2G operator (see the M-ADJ-swap
    mutation).
    """
    sn = _MESHES[geom]()
    A = _loss(sn)
    b = _rand_bulk(sn, 5)
    rt = A.H.apply(_adjoint_inverse(A).apply(b))
    np.testing.assert_allclose(
        np.asarray(rt.interior.values), np.asarray(b.interior.values),
        rtol=_RTOL, atol=1e-11,
        err_msg=f"[{geom}] A*(A*)⁻¹ b ≠ b on the source-carried bulk",
    )


# ── G4 — Mode-11 wrap sentinel (the reverse scan must EXECUTE) ──────────


@pytest.mark.parametrize("geom", list(_MESHES))
def test_gate4_reverse_scan_executes_forward_solve_untouched(geom, monkeypatch):
    r"""``A.H.inverse().apply(b)`` ENTERS ``InvertibleOperator.solve_transpose``
    (``> 0``) and NEVER the forward ``solve`` (``== 0``).

    The Mode-11 gold standard: an in-process wrap counter on the internal
    readers.  A regression that wired ``.H.inverse().apply`` to the forward
    ``solve`` (or routed around the reverse-scan entirely) would leave the
    ``solve_transpose`` counter at 0 and red — the routed-around path cannot
    fake the wrap.
    """
    sn = _MESHES[geom]()
    A = _loss(sn)
    b = _rand_bulk(sn, 5)
    counts = {"solve": 0, "solve_transpose": 0}
    real_solve = InvertibleOperator.solve
    real_solve_transpose = InvertibleOperator.solve_transpose

    def _wrap_solve(self, rhs, *args, **kwargs):
        counts["solve"] += 1
        return real_solve(self, rhs, *args, **kwargs)

    def _wrap_solve_transpose(self, rhs, *args, **kwargs):
        counts["solve_transpose"] += 1
        return real_solve_transpose(self, rhs, *args, **kwargs)

    monkeypatch.setattr(InvertibleOperator, "solve", _wrap_solve)
    monkeypatch.setattr(InvertibleOperator, "solve_transpose", _wrap_solve_transpose)

    _adjoint_inverse(A).apply(b)

    if not counts["solve_transpose"] > 0:
        pytest.fail(
            f"[{geom}] Mode-11: A.H.inverse().apply(b) did NOT enter "
            f"InvertibleOperator.solve_transpose (count=0) — the reverse scan "
            f"was routed around; the gate is vacuous for the new line."
        )
    if not counts["solve"] == 0:
        pytest.fail(
            f"[{geom}] Mode-11: A.H.inverse().apply(b) touched the FORWARD "
            f"solve (count={counts['solve']}) — the adjoint-inverse must route "
            f"the reverse-scan, not the forward sweep."
        )


# ── M-ADJ-swap — drop the ``.H`` in ``_AdjointOperator.inverse`` ───────


@pytest.mark.parametrize("geom", list(_MESHES))
def test_mutation_adj_swap_reds_gates_1_and_2(geom, monkeypatch):
    r"""M-ADJ-swap: ``_AdjointOperator.inverse`` → ``self.inner.inverse()`` (no
    ``.H``) reds BOTH G1 and G2.

    Dropping the swap law's ``.H`` routes ``A.H.inverse()`` to the FORWARD
    ``SweepOperator`` (``A.solve``, not ``A.solve_transpose``), so
    ``A.H.inverse().apply(b) = A⁻¹b`` ≠ ``(A*)⁻¹b``.  G2 (vs the true
    ``A.inverse().H``) reds at O(1); G1's reciprocity reds because ``A A⁻¹ ≠
    A*`` for the non-symmetric het+2G operator.  In-process ``monkeypatch``,
    ``-O``-safe.
    """
    sn = _MESHES[geom]()
    A = _loss(sn)
    b = _rand_bulk(sn, 5)
    psi = _rand_full(sn, 4)
    ref = A.inverse().H.apply(b)  # the true adjoint-inverse (G2 RHS, swap-invariant)

    monkeypatch.setattr(
        _AdjointOperator, "inverse", lambda self: self.inner.inverse()
    )
    x = _adjoint_inverse(A).apply(b)

    g2 = float(np.max(np.abs(np.asarray(x.interior.values) - np.asarray(ref.interior.values))))
    g1 = _reciprocity_rel(A, psi, x, b, sn.full_field_space)
    if not (g2 > _CAUGHT and g1 > _CAUGHT):
        pytest.fail(
            f"[{geom}] M-ADJ-swap NOT caught: G2 bulk-diff={g2:.2e}, "
            f"G1 rel={g1:.2e} (need both > {_CAUGHT:.0e}) — dropping the swap "
            f"law's .H is invisible; check the het+2G non-symmetry."
        )


# ── M-ADJ-metric — skip the metric wrap in ``_AdjointOperator.apply`` ──


@pytest.mark.parametrize("geom", list(_MESHES))
def test_mutation_adj_metric_reds_gate1(geom, monkeypatch):
    r"""M-ADJ-metric: ``_AdjointOperator.apply`` → ``inner.apply_transpose(y)``
    (skip the ``G⁺``/``G`` wrap) reds G1 on ALL geometries.

    SPEC DEVIATION (flagged): the spec predicted this stays GREEN on the slab
    ("slab metric trivial") — a slab-vs-curvilinear discriminator.  EMPIRICALLY
    the slab metric ``G = V·w_n ⊕ |Ω·n|·w_n`` is non-trivial, so ``A* = G⁻¹AᵀG
    ≠ Aᵀ`` on the slab too and G1 reds on EVERY geometry (slab ~0.33, sphere
    ~0.19, cyl ~0.13).  The ``.H``≠Euclidean claim holds everywhere (the metric
    wrap is load-bearing on every geometry — stronger coverage than the intended
    discriminator).  The verification's ``⟨·,·⟩_G`` uses the production metric
    directly (NOT ``_AdjointOperator.apply``), so the monkeypatch reds ``x``
    alone → asymmetry → caught.
    """
    sn = _MESHES[geom]()
    A = _loss(sn)
    b = _rand_bulk(sn, 5)
    psi = _rand_full(sn, 4)
    monkeypatch.setattr(
        _AdjointOperator, "apply", lambda self, y: self.inner.apply_transpose(y)
    )
    x = _adjoint_inverse(A).apply(b)
    rel = _reciprocity_rel(A, psi, x, b, sn.full_field_space)
    if not rel > _CAUGHT:
        pytest.fail(
            f"[{geom}] M-ADJ-metric NOT caught (G1 rel={rel:.2e} ≤ "
            f"{_CAUGHT:.0e}) — the G⁺/G metric wrap in _AdjointOperator.apply is "
            f"not load-bearing here; A* would equal the Euclidean Aᵀ."
        )


# ── predicate flips — the ``is_adjointable`` flip is SURGICAL ───────────


def test_predicate_flip_is_surgical():
    r"""``SweepOperator.is_adjointable`` flips True (over ``InvertibleOperator``);
    the sibling wrap-delegates ``InverseOperator`` / ``GreenOperator`` STAY False.

    Forward-safety: 2.5c adds ``SweepOperator.apply_transpose`` and flips ITS
    ``is_adjointable``, keyed on ``isinstance(inner, InvertibleOperator) and
    inner.is_adjointable`` — the two OTHER ``InverseWrapMixin`` siblings must NOT
    inherit the flip (they have no reverse-scan).
    """
    sn = _slab()
    A = _loss(sn)

    sweep = A.inverse()
    if not (isinstance(sweep, SweepOperator) and sweep.is_adjointable is True):
        pytest.fail(
            f"SweepOperator over InvertibleOperator must be adjointable: "
            f"type={type(sweep).__name__}, is_adjointable={sweep.is_adjointable}"
        )

    green = (A - SNBoundaryOperator(sn)).inverse()
    if not (isinstance(green, GreenOperator) and green.is_adjointable is False):
        pytest.fail(
            f"GreenOperator must STAY non-adjointable (no reverse-scan): "
            f"type={type(green).__name__}, is_adjointable={green.is_adjointable}"
        )

    inverse_op = InverseOperator(A)
    if inverse_op.is_adjointable is not False:
        pytest.fail(
            f"InverseOperator must STAY non-adjointable (the flip is surgical to "
            f"SweepOperator): is_adjointable={inverse_op.is_adjointable}"
        )


# ── static pin (pyright-only, never run — no ``test_`` prefix) ──────────


def _static_typing_pins(endo: LinearOperator[FullField, FullField]) -> None:
    r"""STATIC pin: the adjoint-inverse ``A.H.inverse()`` composes to a
    ``LinearOperator`` with the carriers preserved.

    ``.H`` lives on the base (eager :class:`MissingAdjoint` gate, no narrowing);
    ``.inverse()`` lives on :class:`SupportsInverse`, so ``invertible(...)``
    narrows the adjoint before the call — the Design-C checked bridge (mirrors
    ``test_operators_apply_typed._typeguard_bridge_narrowing_static_pins``).
    """
    if TYPE_CHECKING:
        ah = endo.H
        if invertible(ah):
            assert_type(ah.inverse(), LinearOperator[FullField, FullField])


# ═══════════════════════════════════════════════════════════════════════
# B.2d → step 5 — the CARRYING arm: the fused wrap refuses eagerly; the
# coupled M (the honest triangular grid since 5b; the fused bridge
# deleted at 5d) carries the joint swap law (the same G1–G4 claims)
# ═══════════════════════════════════════════════════════════════════════


def _coupled_M(sn):
    from tests.sn._test_helpers import joint_m_grid

    return joint_m_grid(sn, _loss(sn))


def _coupled_bulk_b(sn, seed: int):
    # b is a SOURCE-side cotangent — System B's member rides the source
    # composite (role-honest: the transposed substitution's rhs update
    # ``b_B − Seedingᵀ·x̄_A`` is a source-role ``−``; the flux-role zeros
    # the fused bridge tolerated were a role leak its role-blind buffers
    # never noticed — step 5's member algebra does).
    from orpheus.numerics.coupled_system import CoupledField
    from orpheus.transport.radial_characteristic_field import (
        RadialCharacteristicField,
    )

    b = _rand_bulk(sn, seed)
    return CoupledField(
        systems=(b, RadialCharacteristicField.source_zeros_on(sn)),
    )


def _coupled_full_psi(sn, seed: int):
    from orpheus.numerics.coupled_system import CoupledField
    from orpheus.transport.radial_characteristic_field import (
        RadialCharacteristicField,
    )

    psi = _rand_full(sn, seed)
    return CoupledField(systems=(psi, RadialCharacteristicField.from_mesh(sn)))


def test_carrying_block_and_coupled_adjointability_coexist():
    r"""The step-6 predicate row (R-6.2 — the B.2d third factor retired
    with the two-channel collapse): on the carrying sphere the bare
    ``A.inverse()`` IS adjointable — the leg-less surface is unambiguously
    the ray-decoupled ``(A, A)`` diagonal block, whose reverse-scan is
    well-defined (the type says what you hold) — AND the COUPLED sibling
    ``M.inverse()`` is adjointable too: the JOINT adjoint's one home is
    the grid's transposed substitution, the block adjoint's the block."""
    sn = _sphere()
    A = _loss(sn)
    if not A.inverse().is_adjointable:
        pytest.fail("the (A,A)-block SweepOperator refuses adjointability on "
                    "a carrying mesh — the R-6.2 two-factor predicate "
                    "regressed to the retired B.2d third factor")
    if not getattr(A.H, "is_invertible", False):
        pytest.fail("the (A,A)-block A.H does not advertise invertibility on "
                    "a carrying mesh (the swap law's block arm)")
    M_op, _space = _coupled_M(sn)
    if not M_op.inverse().is_adjointable:
        pytest.fail("the coupled M.inverse() is not adjointable — the joint "
                    "swap law lost its home")
    if not invertible(M_op.H):
        pytest.fail("the coupled M.H does not advertise invertibility")
    # Step 5: the joint home is the SUBSTITUTION wrap over the triangular
    # grid — pin the concrete successor types (the retired fused bridge
    # must not silently come back as an untyped shim).
    from orpheus.numerics.coupled_system import (
        CoupledOperator,
        CoupledSubstitutionOperator,
    )
    if not isinstance(M_op, CoupledOperator):
        pytest.fail(f"the joint M is {type(M_op).__name__}, not the grid")
    if not isinstance(M_op.inverse(), CoupledSubstitutionOperator):
        pytest.fail("M.inverse() is not the block substitution wrap")


def test_coupled_gate1_forward_matvec_g_reciprocity():
    r"""G1 on the pair: ``⟨M.apply(ψ), x⟩_G = ⟨ψ, b⟩_G`` for
    ``x = M.H.inverse().apply(b)`` — the joint forward + the member-wise
    coupled metric only (never the transpose path)."""
    sn = _sphere()
    M_op, space = _coupled_M(sn)
    b = _coupled_bulk_b(sn, 5)
    psi = _coupled_full_psi(sn, 4)
    ah = M_op.H
    if not invertible(ah):
        pytest.fail("M.H must advertise is_invertible")
    x = ah.inverse().apply(b)
    lhs = space.inner_product(M_op.apply(psi), x)
    rhs = space.inner_product(psi, b)
    rel = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)
    if not rel < _RTOL:
        pytest.fail(f"coupled forward-matvec G-reciprocity broken: rel={rel:.2e}")


def test_coupled_gate2_swap_law_bit_identical():
    r"""G2 on the pair: ``M.H.inverse().apply(b) ≡ M.inverse().H.apply(b)``
    — bit-identical on BOTH systems (the same object graph)."""
    sn = _sphere()
    M_op, _space = _coupled_M(sn)
    b = _coupled_bulk_b(sn, 5)
    ah = M_op.H
    if not invertible(ah):
        pytest.fail("M.H must advertise is_invertible")
    lhs = ah.inverse().apply(b)
    rhs = M_op.inverse().H.apply(b)
    np.testing.assert_array_equal(
        np.asarray(lhs.systems[0].interior.values),
        np.asarray(rhs.systems[0].interior.values),
        err_msg="coupled swap-law bulk not bit-identical",
    )
    np.testing.assert_array_equal(
        np.asarray(lhs.systems[1].to_flat()),
        np.asarray(rhs.systems[1].to_flat()),
        err_msg="coupled swap-law System-B member not bit-identical",
    )


def test_coupled_gate4_reverse_scan_executes(monkeypatch):
    r"""G4 on the pair (Mode-11): ``M.H.inverse().apply(b)`` ENTERS the
    grid's transposed substitution (``CoupledOperator.solve_transpose``
    > 0) and never the forward ``solve`` (== 0)."""
    from orpheus.numerics.coupled_system import CoupledOperator

    sn = _sphere()
    M_op, _space = _coupled_M(sn)
    b = _coupled_bulk_b(sn, 5)
    counts = {"solve": 0, "solve_transpose": 0}
    real_solve = CoupledOperator.solve
    real_st = CoupledOperator.solve_transpose

    def _wrap_solve(self, rhs, *a, **k):
        counts["solve"] += 1
        return real_solve(self, rhs, *a, **k)

    def _wrap_st(self, rhs, *a, **k):
        counts["solve_transpose"] += 1
        return real_st(self, rhs, *a, **k)

    monkeypatch.setattr(CoupledOperator, "solve", _wrap_solve)
    monkeypatch.setattr(CoupledOperator, "solve_transpose", _wrap_st)
    ah = M_op.H
    if not invertible(ah):
        pytest.fail("M.H must advertise is_invertible")
    ah.inverse().apply(b)
    if not counts["solve_transpose"] > 0:
        pytest.fail("coupled M.H.inverse().apply never entered the grid's "
                    "solve_transpose — the transposed substitution was "
                    "routed around")
    if not counts["solve"] == 0:
        pytest.fail(f"coupled M.H.inverse().apply touched the FORWARD solve "
                    f"({counts['solve']}×)")


# ═══════════════════════════════════════════════════════════════════════
# B.2d d3 — A2a: the LOSS-GRID rows in the swap-law home (G-d3.4, R5/R11)
# ═══════════════════════════════════════════════════════════════════════
#
# The rows above live on M (the resolvent factor).  A2a adds the full
# within-group LOSS GRID ``A = M − N`` (``WithinGroupSystem.loss``) to the
# SAME file so the swap-law family is complete in its home: the forward
# ``.H`` reciprocity arm landed live at d3; the inverse/swap-law arm sat
# ``xfail(strict=False)`` until step 5 flipped ``grid.is_invertible``
# (the materialize/LU route — the space's zero exemplar wired at 5b) and
# was converted to the LIVE gate below at 5e (vv xfail discipline), with
# the E2 positive control ahead of the object-identity row.


def _within_group_grid(sn):
    from orpheus.sn.coupled_system import build_within_group_system

    system = build_within_group_system(sn, sn.material_xs_field())
    return system.loss, system.space


def _coupled_rand_pair(sn, seed: int):
    r"""A coupled pair with EVERY block non-zero (random ψ_A bulk+trace AND a
    random System-B member) — the grid's coupling columns must enter the
    pairing, unlike the zero-member vectors the M rows use."""
    from orpheus.numerics.coupled_system import CoupledField
    from orpheus.transport.radial_characteristic_field import (
        RadialCharacteristicField,
    )

    rng = np.random.default_rng(seed)
    ns = sn.radial_characteristic_field_space.shape[0]
    member = RadialCharacteristicField.from_flat(
        rng.standard_normal(ns),
        RadialCharacteristicField.from_mesh(sn),
    )
    return CoupledField(systems=(_rand_full(sn, seed + 1000), member))


def test_a2a_grid_forward_h_reciprocity(monkeypatch):
    r"""A2a forward arm (G-d3.4): ``⟨A·ψ, x⟩_G = ⟨ψ, A.H·x⟩_G < 1e-12`` on
    the carrying ``WithinGroupSystem.loss`` — the grid row in the swap-law
    home (all four coupling blocks enter through the random members).

    Tooth M-ADJ-metric (mirroring psi_half G-c2.5): stripping the
    CoupledSpace metric conjugation (identity ``apply_metric`` /
    ``apply_inverse_metric``) reds the defect O(1) — a hand-rolled
    EUCLIDEAN block-.H would pass the stripped form, which is exactly the
    ERR-067 reopening this tooth exists to catch."""
    from orpheus.numerics.coupled_system import CoupledSpace

    sn = _sphere()
    grid, space = _within_group_grid(sn)
    psi = _coupled_rand_pair(sn, 25)
    x = _coupled_rand_pair(sn, 26)
    lhs = space.inner_product(grid.apply(psi), x)
    rhs = space.inner_product(psi, grid.H.apply(x))
    defect = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)
    if not defect < 1e-12:
        pytest.fail(f"grid-.H G-reciprocity defect {defect:.3e} ≥ 1e-12 on "
                    f"the carrying loss grid — the metric conjugation or a "
                    f"coupling-block transpose is broken (A2a forward arm)")
    with monkeypatch.context() as m:
        m.setattr(CoupledSpace, "apply_metric", lambda self, f: f)
        m.setattr(CoupledSpace, "apply_inverse_metric", lambda self, f: f)
        lhs2 = space.inner_product(grid.apply(psi), x)
        rhs2 = space.inner_product(psi, grid.H.apply(x))
        defect2 = abs(lhs2 - rhs2) / (abs(lhs2) + abs(rhs2) + 1e-300)
    if not defect2 > 1e-3:
        pytest.fail(f"the metric-strip tooth left the defect at "
                    f"{defect2:.3e} — the grid reciprocity gate has no "
                    f"teeth (a Euclidean block-.H would pass)")


def test_a2a_grid_swap_law_inverse_arm():
    r"""A2a inverse arm (G-d3.4, LIVE since 5e): ``A.H.inverse() ≡
    A.inverse().H`` on the FULL within-group loss grid — the R5/R11 swap
    law through the materialize/LU route (``grid.inverse()`` is the
    ``MatrixInverseOperator`` EXTRACT; its ``trans=1`` backsolve is the
    adjointable arm ``_AdjointOperator.is_invertible``'s second clause
    demands). E2 positive control FIRST: the adjoint-inverse SOLVES the
    adjoint equation against the structurally-independent dense-``Aᵀ`` LU
    (a bare bit-identity of two broken constructions would be a false
    green); THEN the object-identity row (both spellings run the SAME
    ``A.inverse().H`` graph — bit-identical BY CONSTRUCTION).
    """
    sn = _sphere()
    grid, _space = _within_group_grid(sn)
    if not grid.is_invertible:
        pytest.fail("grid.is_invertible is False — the step-5 materialize/"
                    "LU route regressed (the space's zero exemplar?)")
    b = _coupled_bulk_b(sn, 5)
    # E2 — the positive control: the adjoint-inverse SOLVES the adjoint
    # equation, checked through the equation itself (metric-closed —
    # ``A.H.apply(A.H.inverse().apply(b)) ≈ b`` — no hand G-algebra to
    # drift; the LU conditioning prices the tolerance).
    x = grid.H.inverse().apply(b)
    round_trip = grid.H.apply(x)
    b_scale = float(np.abs(np.asarray(b.to_flat())).max())
    np.testing.assert_allclose(
        np.asarray(round_trip.to_flat()),
        np.asarray(b.to_flat()),
        rtol=1e-9, atol=1e-10 * b_scale,
        err_msg="E2 positive control: A.H.apply(A.H.inverse().apply(b)) "
                "does not reproduce b — the adjoint-inverse does not solve "
                "the adjoint equation",
    )
    lhs = x
    rhs = grid.inverse().H.apply(b)
    np.testing.assert_array_equal(
        np.asarray(lhs.systems[0].interior.values),
        np.asarray(rhs.systems[0].interior.values),
        err_msg="grid swap-law bulk not bit-identical",
    )
    np.testing.assert_array_equal(
        np.asarray(lhs.systems[1].to_flat()),
        np.asarray(rhs.systems[1].to_flat()),
        err_msg="grid swap-law System-B member not bit-identical",
    )
