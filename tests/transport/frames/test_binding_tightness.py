r"""XD-1 — the scattering binding's tightness gate (CS4c step 3; design
record §6 as amended, verification plan §5).

The ERR-039 family (a claimed :math:`\Pi^* = R` that wasn't) is caught by
**leg (i)**: swept as matrices under the frame's OWN space metrics
(``measure_space`` / F-0 Parseval ``basis_space`` — no hand-derived
Grams), the shipped analysis face's Hilbert adjoint equals the
reconstruction over the measure's total weight,
:math:`M^{\dagger} = R/W` — and a face re-minted with a WRONG embedding
(constant :math:`\Sigma w/N`, or unweighted) REDS by the measured G6.3
87 %-class margin.

**Leg (ii)**: multiplicativity of the binding on the Funk–Hecke
eigenstructure — ``bind(K₁∘K₂) = bind(K₁)·bind(K₂)`` — holds iff the
frame is tight on the spanned harmonics. The product operand is built
GATE-SIDE from the ℓ-eigenvalue stacks through the frame's REAL faces
(F7's defer-until-2: no kernel product API); in the stored
``[g_from, g_to]`` convention a group-bearing product would compose as
``moments₂ @ moments₁`` — here the operands are zonal eigenvalue stacks,
where composition is the elementwise product (Funk–Hecke).

⛔ NOT spelled here (§6's refutation, 2026-08-30): the route-(b)
comparison ``bind(K).H == M†K†R†`` — ``(RKM)† = M†K†R†`` is an algebraic
identity of the metric adjoint, measured ≤ 2.24e-16 under correct AND
wrong embeddings alike; it cannot red and is recorded as a structural
THEOREM only (vv #19).
⚠ Leg (ii)'s negative control is the equispaced-equal-weight rule by
MEASURED bite. The pre-carve plan (§1.2) additionally warned that
``gauss_legendre(L)`` is blind on zonal multiplicativity (≤ 5.9e-16
over 200 draws) — [M] 2026-08-30, re-measured THROUGH THE SHIPPED FACES
(F-0 Parseval metrics + the producer-side /W): it is NOT blind in this
spelling (rel. 3.0/4.3 at L = 2/3), so that hazard was a property of
the probe's raw-table binding, not of this gate's construction. The
equispaced control stays as the shipped control for its measured band;
no negative-control-of-the-control row is asserted (its premise does
not hold here — the disagreement is recorded per §4's verify rule).

**Leg (iii)** records the ℓ=0 blindness as an ASSERTED row: at L = 0
both rules read clean on both legs, so an ℓ=0 gate discriminates
nothing — the ℓ≥1 rows are the coverage, and the recorded blindness
cannot silently stop being true.

Draw discipline (vv #31): seeds pinned; the BANDS are the claim, not the
point values.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.manifold import COSINE_INTERVAL, SPHERE
from orpheus.numerics.symmetry import SubgroupOfO3
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.frames.harmonic_frame import HarmonicFrame

pytestmark = pytest.mark.foundation


def _sphere_frame(L: int) -> HarmonicFrame:
    """Lebedev order-17 — tight through degree 2L+1 for every L ≤ 3,
    with genuinely NON-uniform weights (the wrong-embedding controls
    would degenerate on an equal-weight rule)."""
    return HarmonicFrame.from_galerkin(
        Quadrature.lebedev(order=17).angular_frame(L),
    )


def _equispaced_frame(L: int, n: int) -> HarmonicFrame:
    r"""Equispaced-equal-weight on :math:`\mu` — the measured red-capable NON-tight control (total weight 2).

    ⛔ **RE-DECLARED 2026-09-02 (#429).** This fixture built its measure with
    ``support=COSINE_INTERVAL`` — the bare interval :math:`[-1,1]`. That is
    the *realization* of an angular rule's point set, not the point set: a
    1-D angular rule integrates over :math:`S^2/O(2)_x`, whose realization
    happens to be that interval. Post-repair the distinction is enforced at
    two doors, and this fixture hit **both**:

    1. ``Quadrature._harmonic_basis``'s ``match support:`` — an ``Interval``
       matches no arm, so the dispatch refuses ("harmonic machinery on a rule
       whose measure lives on [-1,1]");
    2. the frame's G0 arrow — `[M]` ``domain.realization == COSINE_INTERVAL``
       is **True** while ``domain == COSINE_INTERVAL`` is **False**, so a
       G0 comparing objects refuses where one comparing realizations would
       admit.

    Declaring the orbit space is the honest fix, not a workaround: it is
    exactly what every shipped 1-D rule now declares. ⭐ And it PAYS — with
    the Legendre basis the entry carries, this measure is the tree's only
    dense-AND-full-rank 1-D frame (`[M]` ``n = 8, L = 3``: offdiag
    ``2.222e-01``, rank 4/4, **0 dead slots**), because a theorem forbids a
    Gauss rule from being one (``GL_n``'s Legendre Gram is diagonal-and-exact
    for ``L <= n-1`` and structurally dead at :math:`\ell = n`).

    ⚠ The ``COSINE_INTERVAL``-declaring twin is kept as G0's constructed
    refusal witness — see
    :func:`test_a_bare_interval_measure_is_refused_by_the_dispatch_before_g0`.
    """
    quad = Quadrature(measure=_equispaced_measure(n, support="orbit"))
    return HarmonicFrame.from_galerkin(quad.angular_frame(L))


def _equispaced_measure(n: int, *, support: str) -> DiscreteMeasure:
    """The equispaced-equal-weight measure, declared on the orbit space or on its bare realization."""
    mu = np.linspace(-1.0, 1.0, n + 2)[1:-1]
    w = np.full(mu.size, 2.0 / mu.size)
    declared = (
        SPHERE.quotient(SubgroupOfO3.O2("x")) if support == "orbit"
        else COSINE_INTERVAL
    )
    return DiscreteMeasure(nodes=mu, weights=w, support=declared)


def _sweep(op, shape_in, n_in):
    cols = []
    for i in range(n_in):
        e = np.zeros(n_in)
        e[i] = 1.0
        cols.append(np.asarray(op.apply(e.reshape(shape_in))).ravel())
    return np.array(cols).T


def _faces(frame):
    tab = frame.table
    # ⛔ head-GENERIC since 2026-09-02 (#429): this read ``shape[1] * shape[2]``,
    # i.e. it assumed the RECTANGULAR spherical-harmonic head. A 1-D rule's
    # frame carries the FLAT Legendre head, table ``(N, L+1)``, and the old
    # spelling raised ``IndexError`` on it.
    n_modes = int(np.prod(tab.shape[1:]))
    shape = tab.shape[1:]
    M = _sweep(frame.analysis, (tab.shape[0],), tab.shape[0]) if False else None
    # analysis: ordinate -> coeff; sweep over ordinates
    cols = []
    for i in range(tab.shape[0]):
        e = np.zeros(tab.shape[0])
        e[i] = 1.0
        cols.append(np.asarray(frame.analysis.apply(e)).ravel())
    M = np.array(cols).T                       # (n_modes, N)
    R = _sweep(frame.reconstruction, shape, n_modes)   # (N, n_modes)
    return M, R


def _metric_adjoint(frame, M):
    r"""``M† = G_V⁻¹ Mᵀ G_C`` via the frame's OWN space verbs — the same
    metric route the shipped ``.H`` takes, applied to an arbitrary
    analysis MATRIX (which is what lets a doctored face be adjointed
    under the TRUE metric rather than self-consistently under its own)."""
    tab = frame.table
    n_modes = M.shape[0]
    shape = tab.shape[1:]
    cols = []
    for i in range(n_modes):
        c = np.zeros(n_modes)
        c[i] = 1.0
        gc = np.asarray(frame.basis_space.apply_metric(c.reshape(shape))).ravel()
        cols.append(
            np.asarray(frame.measure_space.apply_inverse_metric(M.T @ gc)).ravel()
        )
    return np.array(cols).T                    # (N, n_modes)


class TestLegIGalerkinOnTheFaces:
    """G-D1 — ``M† == R/W`` under the shipped embedding; REDS under a
    wrong embedding (the ERR-039 catcher). [M] at L = 1/2/3:
    correct ≤ 1.1e-15; constant 3.3e-1; unweighted 8.7."""

    @pytest.mark.parametrize("L", [1, 2, 3])
    def test_shipped_embedding_closes(self, L):
        frame = _sphere_frame(L)
        M, R = _faces(frame)
        W = float(frame.measure.weights.sum())
        target = R / W
        cols = []
        n_modes = M.shape[0]
        shape = frame.table.shape[1:]
        for i in range(n_modes):
            c = np.zeros(n_modes)
            c[i] = 1.0
            cols.append(
                np.asarray(frame.analysis.H.apply(c.reshape(shape))).ravel()
            )
        MH = np.array(cols).T
        defect = np.linalg.norm(MH - target) / np.linalg.norm(target)
        if defect > 1e-13:
            pytest.fail(f"shipped-face Galerkin defect {defect:.3e} > 1e-13")
        # And the matrix-level metric route agrees with the shipped .H
        # (the formula below is what the wrong-embedding controls ride).
        via_formula = _metric_adjoint(frame, M)
        agree = np.linalg.norm(via_formula - MH) / np.linalg.norm(MH)
        if agree > 1e-13:
            pytest.fail(
                f"the control route's metric formula drifted from the "
                f"shipped .H ({agree:.3e}) — the controls would measure "
                f"a different question",
            )

    @pytest.mark.parametrize("L", [1, 2, 3])
    @pytest.mark.parametrize("embedding", ["constant", "unweighted"])
    def test_wrong_embedding_reds(self, L, embedding):
        """Negative controls (vv #11 + #19): the doctored faces are
        hand-built HERE, so the gate is red-capable the day it lands
        (§6c) with no production change."""
        frame = _sphere_frame(L)
        _, R = _faces(frame)
        W = float(frame.measure.weights.sum())
        tab = frame.table
        N = tab.shape[0]
        n_modes = tab.shape[1] * tab.shape[2]
        Y = tab.reshape(N, n_modes)
        w = np.asarray(frame.measure.weights, dtype=float)
        emb = np.full(N, w.sum() / N) if embedding == "constant" else np.ones(N)
        M_wrong = (Y * emb[:, None]).T
        defect = (
            np.linalg.norm(_metric_adjoint(frame, M_wrong) - R / W)
            / np.linalg.norm(R / W)
        )
        if defect < 1e-1:
            pytest.fail(
                f"{embedding} embedding read {defect:.3e} < 1e-1 — the "
                f"leg cannot see the G6.3 embedding defect class",
            )


def _bind(frame, per_l):
    r"""The zonal binding ``(1/W)·R · Λ · M`` through the frame's REAL
    faces, with :math:`\Lambda` the per-ℓ eigenvalue diagonal
    (Funk–Hecke) and the producer-side :math:`1/W` exactly as the
    production kernel applies it (``S_aniso = (1/W)·RΛM``) — on a tight
    rule ``M R = W·I`` on the live harmonics, so the normalized binding
    composes: ``bind(K₁)·bind(K₂) = (1/W²)·RΛ(MR)Λ'M = bind(K₁K₂)``."""
    M, R = _faces(frame)
    # the per-degree eigenvalue is broadcast over the HEAD's own degree
    # blocks (rectangular: 2l+1 slots per degree; flat: one), read off the
    # basis's coefficient space rather than from a layout literal.
    head = frame.basis.space
    lam_table = np.zeros(head.shape)
    for l in range(head.L + 1):
        lam_table[head.degree_block(l)] = per_l[l]
    lam = lam_table.ravel()
    W = float(np.asarray(frame.measure.weights).sum())
    return (R @ (lam[:, None] * M)) / W


class TestLegIIMultiplicativity:
    r"""G-D2 — ``bind(K₁∘K₂) == bind(K₁)·bind(K₂)`` at ℓ ≥ 1: tight rule ≤ 1e-13; equispaced control ≥ 1e-3.

    Bands, not points (``vv-principles`` #31). ⭐ **RE-MEASURED 2026-09-02
    (#429)**, because the control's frame changed underneath it: the
    equispaced rule's re-declared support makes it bind the Legendre basis
    (table ``(N, L+1)``) instead of the full harmonics, so the numbers had to
    be re-taken rather than inherited. `[M]` equispaced control now reads
    ``0.4130 / 0.3481 / 0.2502`` at ``L = 1/2/3`` (was ``≥ 1e-3``); the tight
    rule reads ``3.5e-16 / 5.0e-16 / 6.9e-16``. The band is unchanged and the
    separation is ~15 orders, so the leg is if anything stronger.
    """

    @pytest.mark.parametrize("L", [1, 2, 3])
    def test_tight_rule_is_multiplicative(self, L):
        rng = np.random.default_rng(20260830)
        frame = _sphere_frame(L)
        k1 = rng.uniform(0.2, 1.0, size=L + 1)
        k2 = rng.uniform(0.2, 1.0, size=L + 1)
        left = _bind(frame, k1 * k2)
        right = _bind(frame, k1) @ _bind(frame, k2)
        rel = np.linalg.norm(left - right) / np.linalg.norm(left)
        if rel > 1e-13:
            pytest.fail(f"tight rule multiplicativity broke: {rel:.3e}")

    @pytest.mark.parametrize("L", [1, 2, 3])
    def test_non_tight_control_reds(self, L):
        rng = np.random.default_rng(20260830)
        frame = _equispaced_frame(L, n=L + 2)
        k1 = rng.uniform(0.2, 1.0, size=L + 1)
        k2 = rng.uniform(0.2, 1.0, size=L + 1)
        left = _bind(frame, k1 * k2)
        right = _bind(frame, k1) @ _bind(frame, k2)
        rel = np.linalg.norm(left - right) / np.linalg.norm(left)
        if rel < 1e-3:
            pytest.fail(
                f"equispaced control read {rel:.3e} < 1e-3 — the leg "
                f"lost its red capability",
            )


class TestLegIIIRecordedL0Blindness:
    """G-D3 — at L = 0 BOTH rules read clean: an ℓ=0 gate discriminates
    nothing; the ℓ≥1 rows above are the coverage. Asserted, not prose."""

    @pytest.mark.parametrize(
        "make",
        [_sphere_frame, lambda L: _equispaced_frame(L, 4)],
        ids=["tight", "equispaced"],
    )
    def test_l0_is_blind(self, make):
        frame = make(0)
        rng = np.random.default_rng(7)
        k1 = rng.uniform(0.2, 1.0, size=1)
        k2 = rng.uniform(0.2, 1.0, size=1)
        left = _bind(frame, k1 * k2)
        right = _bind(frame, k1) @ _bind(frame, k2)
        rel = np.linalg.norm(left - right) / max(np.linalg.norm(left), 1e-300)
        np.testing.assert_array_less(rel, 1e-13)



# ── G0's constructed refusal witness (#429 tracker 2.2, plan-authoring §6c) ──


class TestABareIntervalMeasureIsRefused:
    r"""The twin of :func:`_equispaced_frame` that declares the bare interval.

    ⛔ **§6c — G0 lands with no shipped production refusal.** `[M]` 2026-09-02:
    after the fused commit every support the dispatch selects picks a basis G0
    admits, so nothing production BUILDS is rejected. A gate whose rejected
    input does not exist ships green and structurally unable to fail. This
    class constructs the input.

    ⚠ **The fixture hits THREE guards in sequence, and the ORDER is the
    claim.** Asserting only "it raises" would be satisfied by any of them, so
    the row below pins the FIRST message PRESENT and the other two ABSENT —
    which is what makes it a wiring-order gate rather than a smoke test.
    """

    def test_a_bare_interval_measure_is_refused_by_the_dispatch_before_g0(self) -> None:
        r"""Guard 1 fires: the basis DISPATCH, not the frame's G0 arrow.

        ``Quadrature._harmonic_basis`` reads the measure's support and matches
        ``Sphere`` / a sphere ``Quotient``; a bare ``Interval`` matches
        neither, so the refusal happens before a basis exists to bind. Its
        message names the point set; G0's names a missing quotient MAP; the
        door's names a truncation order. The three are disjoint by design.
        """
        quad = Quadrature(measure=_equispaced_measure(8, support="interval"))
        with pytest.raises(NotImplementedError) as excinfo:
            quad.angular_frame(3)
        message = str(excinfo.value)

        assert "harmonic machinery on a rule whose measure lives on" in message
        assert "[-1,1]" in message
        # …and the two downstream guards did NOT fire — if either fragment
        # appears, the wiring order changed and this gate is measuring a
        # different guard than the one it names.
        assert "no quotient map" not in message
        assert "truncation order" not in message

    def test_the_realization_and_the_orbit_space_are_DIFFERENT_objects(self) -> None:
        r"""The measurement that makes the re-declaration a repair rather than a workaround.

        `[M]` 2026-09-02: the orbit space's ``realization`` **==**
        ``COSINE_INTERVAL`` (``True``) while the orbit space itself does not
        (``False``). So a G0 comparing REALIZATIONS and a G0 comparing OBJECTS
        give OPPOSITE verdicts on this one fixture — which is why the
        predicate had to be ruled rather than inferred, and why declaring the
        realization was never harmless.
        """
        entry = SPHERE.quotient(SubgroupOfO3.O2("x"))
        assert entry.realization == COSINE_INTERVAL
        assert entry != COSINE_INTERVAL

    def test_the_re_declared_twin_is_ADMITTED_and_binds_the_legendre_basis(self) -> None:
        r"""The positive leg (``vv-principles`` #11) — same nodes, same weights, honest support.

        Nothing about the RULE changed: the numbers are bit-identical. Only
        the declaration moved, and that is enough to make the frame
        constructible — which is the whole content of "angular spaces are
        derived from symmetry".
        """
        from orpheus.numerics.basis.legendre_basis import LegendreBasis

        interval = _equispaced_measure(8, support="interval")
        orbit = _equispaced_measure(8, support="orbit")
        np.testing.assert_array_equal(interval.nodes, orbit.nodes)
        np.testing.assert_array_equal(interval.weights, orbit.weights)

        frame = _equispaced_frame(3, n=8)
        assert isinstance(frame.basis, LegendreBasis)
        assert frame.table.shape == (8, 4)
        # ⭐ and it is the dense-AND-full-rank 1-D witness a Gauss rule cannot be
        gram = np.asarray(frame.discrete_gram)
        assert np.linalg.matrix_rank(gram) == 4
        assert float(np.max(np.abs(gram - np.diag(np.diagonal(gram))))) > 0.1
