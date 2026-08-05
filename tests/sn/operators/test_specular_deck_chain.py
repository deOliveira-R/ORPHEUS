r"""The specular deck transformation as a **length-1 chain**, and its binding.

Campaign step **G6.3 step 5** (#330), sibling of ``test_lambertian_chain.py``.
The Lambertian is a two-link chain :math:`\Gamma_+(f) \to S(f) \to
\Gamma_-(f)`; a measure-preserving bijection has nothing to factor, so the
specular mirror is the **same structure with one link** — no separate "atomic"
arm, just a shorter chain. ``_specular_kernel`` now returns its
:class:`~orpheus.numerics.operator.PermutationOperator` bound
:math:`\Gamma_+(f) \to \Gamma_-(f)`.

What these gates are FOR — the mutation measurement that shaped them
--------------------------------------------------------------------

Before this module existed, the binding was gated by **nothing**. Measured over
``tests/geometry tests/sn/operators -m "not slow"`` (baseline ``3 failed /
1668 passed``, ``_specular_kernel`` entered 1252×), with an in-process plugin
rebinding the kernel:

.. list-table::
   :header-rows: 1
   :widths: 58 20

   * - mutation
     - new reds
   * - **positive control** — identity permutation, binding intact
     - **+23**
   * - drop the binding entirely
     - **0**
   * - ⭐ **swap** ``domain`` ↔ ``codomain``
     - **0**
   * - bind both ends to :math:`\Gamma_+`
     - **0**

⭐ **The swap is the mutation to design for.** It is a one-word slip; it
survives :func:`~orpheus.numerics.operator.checked_space_extent` because
:math:`|\Gamma_+| = |\Gamma_-|` on every shipped quadrature; it changes no
arithmetic, so every value gate stays green. And once step 8 routes
``_reflect_trace`` through ``@`` it makes the *legal* composition raise and the
*illegal* one pass — a fault whose symptom is the inverse of its cause.
:class:`TestTheBindingIsTheRightWayRound` is its **only** possible catcher: an
assertion naming which space is which end.

**Re-run with this module present** (baseline ``3 failed / 1747 passed``), and
note *where* the reds land — a red in a gate with no plausible view of the
mutated property is not coverage (``vv-principles`` anti-pattern #18):

.. list-table::
   :header-rows: 1
   :widths: 24 12 64

   * - mutation
     - new reds
     - which rows, and why they can see it
   * - identity permutation
     - **+28**
     - the positive control — it changes the arithmetic, so value gates red
   * - ⭐ **swap**
     - **+10**
     - :class:`TestTheBindingIsTheRightWayRound` **only**. Nothing else, in
       the whole slice, can see it.
   * - both ends → :math:`\Gamma_+`
     - **+18**
     - the binding rows, plus :class:`TestTheSquareDoesNotExist` (with equal
       ends ``P @ P`` becomes *legal*) and the codomain-metric leg
   * - drop the binding
     - **+31**
     - the above plus every :class:`TestTheHilbertAdjointCarriesTheMetric`
       row — unbound, ``.H`` is the bare transpose and cannot move

⚠ The swap does **not** red the ``.H`` rows, and that is correct rather than a
gap: the swapped adjoint is :math:`G_-^{-1}P^{\mathsf T}G_+`, which the
metric-preservation theorem below makes numerically equal to the unswapped
one. The metric is blind to the swap **for the same reason** ``.H`` is blind to
the metric here.

The retired flag, and why the algebra replaced it
-------------------------------------------------

``PermutationOperator.is_involution`` was retired in this step rather than
corrected. `[M]` the narrowed local permutation satisfies
``perm[perm] == arange`` on ``gauss_legendre(4/8)``, ``product(4,4)`` and
``level_symmetric(6)`` and NOT on ``lebedev(17)`` — one physical law, an answer
tracking the quadrature's local index ordering rather than the mirror. Bound,
the question stops being answerable: :math:`P \circ P` is not an expression
between two different spaces. So the claim is now made by the algebra
(:class:`TestTheSquareDoesNotExist`), which cannot hold a wrong value because
it holds no value. See :ref:`bc-narrowed-involution`.

⚠ **Known blindness, measured not assumed.** Every *value* gate here
(bit-identity, ``.H``, the metric criterion) is blind to the swap, and
:class:`TestTheSquareDoesNotExist` is blind to it too — the spaces still differ,
so ``P @ P`` still raises. Credit the swap to
:class:`TestTheBindingIsTheRightWayRound` and to nothing else.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from orpheus.geometry.boundary import ReflectiveBoundary
from orpheus.numerics.face_layout import FaceLayout
from orpheus.numerics.operator import (
    IncompatibleOperatorComposition,
    PermutationOperator,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces import AngularTraceSpace
from orpheus.sn.boundary.realizer import (
    SNBoundaryRealizer,
    _outflow_restriction,
    _specular_kernel,
)
from orpheus.sn.mesh.method_space import SNMethodSpace

pytestmark = pytest.mark.l1

_NG = 2

# ``gauss_legendre(4)`` is deliberately absent from the default set: it has no
# tangential ordinates at all, so it nulls the very tier the half-traces exist
# to separate. It appears only where a case explicitly wants that degeneracy.
_CASES: dict[str, tuple] = {
    "gauss_legendre(8) xmin": (
        lambda: Quadrature.gauss_legendre(n_ordinates=8), "xmin", "x",
    ),
    "gauss_legendre(8) xmax": (
        lambda: Quadrature.gauss_legendre(n_ordinates=8), "xmax", "x",
    ),
    "product(4,4) xmin": (
        lambda: Quadrature.product(n_mu=4, n_phi=4), "xmin", "x",
    ),
    "level_symmetric(6) xmin": (
        lambda: Quadrature.level_symmetric(sn_order=6), "xmin", "x",
    ),
    "lebedev(17) xmin": (lambda: Quadrature.lebedev(order=17), "xmin", "x"),
}


def _fixture(case: str):
    """``(quad, trace, method_space, gamma_out, P, axis, face)`` for a case."""
    factory, face, axis = _CASES[case]
    quad = factory()
    layout = FaceLayout.from_named_shapes(
        [("xmin", (int(quad.N), _NG)), ("xmax", (int(quad.N), _NG))]
    )
    trace = AngularTraceSpace.from_quadrature_and_layout(quad, layout)
    method_space = SNMethodSpace.for_face(
        quadrature=quad, face=face, trace=trace
    )
    gamma_out = _outflow_restriction(method_space, "reflective")
    kernel = _specular_kernel(
        quad, method_space, gamma_out, axis=axis, law_key="reflective"
    )
    return quad, trace, method_space, gamma_out, kernel, axis, face


def _unbound_twin(kernel: PermutationOperator) -> PermutationOperator:
    """The same permutation with NO spaces — the control for every gate here.

    Built from a COPY of the array, so the twin shares no state with
    production: any difference in behaviour is attributable to the binding
    alone and not to aliasing.
    """
    return PermutationOperator(np.asarray(kernel.perm).copy(), axis=kernel.axis)


# ─────────────────────────────────────────────────────────────────────
# The square does not exist — the retired flag's claim, asked of the algebra
# ─────────────────────────────────────────────────────────────────────


class TestTheSquareDoesNotExist:
    r"""``P @ P`` is refused; the legal compositions are not.

    This is the whole semantic content of the ``is_involution`` retirement.
    The flag had to *answer*; the algebra gets to *refuse*.
    """

    @pytest.mark.parametrize("case", list(_CASES))
    def test_the_mirror_squared_is_not_an_expression(self, case: str) -> None:
        r"""⭐ ``P @ P`` raises, naming BOTH half-traces.

        :math:`P : \Gamma_+ \to \Gamma_-`, so its square would have to feed
        :math:`\Gamma_-` into a domain of :math:`\Gamma_+`. The error names
        both, which is what makes the diagnosis readable — shape cannot,
        since the two halves have identical shape.
        """
        _, _, _, _, kernel, _, _ = _fixture(case)
        with pytest.raises(IncompatibleOperatorComposition) as excinfo:
            _ = kernel @ kernel
        message = str(excinfo.value)
        assert "outflow" in message and "inflow" in message, (
            f"the refusal must name both half-traces so the reader can see "
            f"WHICH end mismatched; got: {message}"
        )

    @pytest.mark.parametrize("case", list(_CASES))
    def test_the_legal_round_trips_are_accepted(self, case: str) -> None:
        r"""Positive control: :math:`P P^{-1}` and :math:`P^{-1} P` compose.

        Without this leg, an operator that refused **every** composition
        would pass the row above. These two are the round trips
        :math:`\Gamma_- \to \Gamma_+ \to \Gamma_-` and its mirror, and both
        are legal by construction.
        """
        _, _, _, _, kernel, _, _ = _fixture(case)
        forward_then_back = kernel @ kernel.inverse()
        back_then_forward = kernel.inverse() @ kernel
        assert forward_then_back.domain is kernel.codomain
        assert forward_then_back.codomain is kernel.codomain
        assert back_then_forward.domain is kernel.domain
        assert back_then_forward.codomain is kernel.domain

    @pytest.mark.parametrize("case", list(_CASES))
    def test_the_unbound_twin_composes_with_itself(self, case: str) -> None:
        """Control: the refusal comes from the BINDING, not from the class.

        The same permutation array, unbound, squares without complaint. So
        ``P @ P`` raising is a statement about the spaces this operator was
        given, not about permutations in general — which is the distinction
        the retired flag could not draw.
        """
        _, _, _, _, kernel, _, _ = _fixture(case)
        twin = _unbound_twin(kernel)
        squared = twin @ twin
        assert squared.domain is None and squared.codomain is None

    def test_the_raw_index_test_says_involution_on_four_of_five(self) -> None:
        r"""⭐ The activation guard for this whole class — and it is reddenable.

        If the narrowed permutation were NOT an index-involution, the refusal
        above would be unremarkable — the square would be meaningless anyway.
        The gates matter precisely because ``perm[perm] == arange`` is
        **True** here, so a stored flag would have said "involution" and the
        algebra says "no such composition".

        ``lebedev(17)`` is excluded and gets its own control row below: it
        already answers ``False`` for the *index* reason, so including it here
        would let this guard pass on a case that proves nothing.
        """
        for case in (
            "gauss_legendre(8) xmin",
            "product(4,4) xmin",
            "level_symmetric(6) xmin",
        ):
            _, _, _, _, kernel, _, _ = _fixture(case)
            perm = np.asarray(kernel.perm)
            assert np.array_equal(perm[perm], np.arange(perm.size)), (
                f"{case}: the narrowed permutation is not an index-involution, "
                f"so the refusal of `P @ P` would be uninformative here — this "
                f"guard exists to keep the class from going quietly vacuous"
            )

    def test_lebedev_is_already_false_for_the_INDEX_reason(self) -> None:
        r"""A **control, not a catcher** — and the measurement behind the retire.

        `[M]` on ``lebedev(17)`` the narrowed permutation is *not* an index
        involution, while on the other three it is. Same physical law, same
        face: the raw index answer tracks how ``to_local``'s ``searchsorted``
        orders the locals, not the mirror. That variation is why a stored
        ``bool`` could not be right — and why the refusal below must be
        credited to the spaces rather than to this row.
        """
        _, _, _, _, kernel, _, _ = _fixture("lebedev(17) xmin")
        perm = np.asarray(kernel.perm)
        assert not np.array_equal(perm[perm], np.arange(perm.size))
        # And the refusal still holds, for the space reason, not this one.
        with pytest.raises(IncompatibleOperatorComposition):
            _ = kernel @ kernel


# ─────────────────────────────────────────────────────────────────────
# The binding is the RIGHT one — the only catcher for the swap
# ─────────────────────────────────────────────────────────────────────


class TestTheBindingIsTheRightWayRound:
    r"""⭐ Which space is which END. The module docstring's swap mutation.

    Asserted with ``is``, not ``==``:
    :meth:`~orpheus.numerics.space.FunctionSpace.__eq__` compares
    ``(name, shape)``, and the trace caches its per-face accessors, so object
    identity is available and is strictly stronger — an ``==`` gate would
    accept a freshly-derived look-alike.
    """

    @pytest.mark.parametrize("case", list(_CASES))
    def test_the_kernel_maps_outflow_to_inflow(self, case: str) -> None:
        """The ends are the face's own Γ₊ and Γ₋, by object identity."""
        _, trace, _, _, kernel, _, face = _fixture(case)
        assert kernel.domain is trace.outflow_space(face)
        assert kernel.codomain is trace.inflow_space(face)

    def test_the_binding_is_to_THIS_face_not_the_opposite_one(self) -> None:
        r"""Cross-face: :math:`\Gamma_+(\text{xmin})` is not :math:`\Gamma_+(\text{xmax})`.

        `[M]` the two have **identical shape** ``(4, 2)`` on
        ``gauss_legendre(8)``, so no shape or extent guard can distinguish
        them. This row is the load-bearing justification for binding spaces at
        all rather than checking lengths.
        """
        _, trace, _, _, kernel, _, _ = _fixture("gauss_legendre(8) xmin")
        assert kernel.domain is not trace.outflow_space("xmax")
        assert kernel.codomain is not trace.inflow_space("xmax")
        assert (
            trace.outflow_space("xmin").shape
            == trace.outflow_space("xmax").shape
        ), (
            "if these shapes ever differ, the cross-face claim above stops "
            "being the space check's exclusive territory — re-read the gate"
        )

    @pytest.mark.parametrize("case", list(_CASES))
    def test_the_gamma_plus_it_consumes_is_bound_too(self, case: str) -> None:
        r"""The realizer's own :math:`\gamma_+` carries its spaces.

        ``_outflow_restriction`` builds its restriction locally rather than
        fetching the trace's cached one (realization runs once per face; the
        cache serves the per-sweep consumer). That is a construction-cost
        decision, so it owes the same typing — otherwise the realizer hands
        its own consumer an untyped domain while typing that consumer's
        output.
        """
        _, trace, _, gamma_out, _, _, face = _fixture(case)
        assert gamma_out.domain is trace.face_space(face)
        assert gamma_out.codomain is trace.outflow_space(face)

    @pytest.mark.parametrize("case", list(_CASES))
    def test_the_whole_face_action_types_end_to_end(self, case: str) -> None:
        r"""⭐ :math:`P \circ \gamma_+ : \Gamma(f) \to \Gamma_-(f)`.

        The payoff of binding both: the full boundary action — restrict to
        the outflow half, then mirror — composes into one typed operator
        whose ends are the face's full trace and its inflow half. Neither
        binding alone gives this.
        """
        _, trace, _, gamma_out, kernel, _, face = _fixture(case)
        face_action = kernel @ gamma_out
        assert face_action.domain is trace.face_space(face)
        assert face_action.codomain is trace.inflow_space(face)


# ─────────────────────────────────────────────────────────────────────
# inverse() inverts the binding
# ─────────────────────────────────────────────────────────────────────


class TestTheInverseSwapsItsEnds:
    r""":math:`P^{-1} : \Gamma_- \to \Gamma_+` — inverted, not carried, not dropped.

    Dropping would have been the quieter bug: the return leg composes with
    anything while the forward leg is fully typed, an asymmetry no gate on
    the forward leg could see.
    """

    @pytest.mark.parametrize("case", list(_CASES))
    def test_the_ends_are_exchanged(self, case: str) -> None:
        """``inverse().domain is codomain`` and vice versa, by identity."""
        _, _, _, _, kernel, _, _ = _fixture(case)
        inverse = kernel.inverse()
        assert inverse.domain is kernel.codomain
        assert inverse.codomain is kernel.domain

    @pytest.mark.parametrize("case", list(_CASES))
    def test_inverting_twice_restores_the_binding_and_the_perm(
        self, case: str
    ) -> None:
        r""":math:`(P^{-1})^{-1} = P`, ends included.

        ``argsort`` of a permutation is exactly involutive in integer math, so
        the array claim is exact; the space claim is what step 5 added.
        """
        _, _, _, _, kernel, _, _ = _fixture(case)
        round_trip = kernel.inverse().inverse()
        assert round_trip.domain is kernel.domain
        assert round_trip.codomain is kernel.codomain
        np.testing.assert_array_equal(round_trip.perm, kernel.perm)

    @pytest.mark.parametrize("case", list(_CASES))
    def test_the_unbound_twins_inverse_stays_unbound(self, case: str) -> None:
        """Negative leg: the swap is produced by the binding, not the method.

        An ``inverse()`` that fabricated spaces from nowhere would pass the
        two rows above and fail this one.
        """
        _, _, _, _, kernel, _, _ = _fixture(case)
        inverse = _unbound_twin(kernel).inverse()
        assert inverse.domain is None and inverse.codomain is None


# ─────────────────────────────────────────────────────────────────────
# The binding costs no arithmetic
# ─────────────────────────────────────────────────────────────────────


class TestTheBindingIsFreeOfArithmetic:
    """Bound production against a genuinely unbound twin, bit-for-bit.

    The comparison is against an operator built from a **copy** of the
    permutation array, not against production re-read — the campaign has
    twice been bitten by a "two implementations" pin that had quietly become
    a value compared with itself.
    """

    @pytest.mark.parametrize("case", list(_CASES))
    def test_apply_and_transpose_are_bit_identical_to_the_unbound_twin(
        self, case: str
    ) -> None:
        """Binding adds metadata; every gathered row is the same row."""
        _, _, _, _, kernel, _, _ = _fixture(case)
        twin = _unbound_twin(kernel)
        assert twin is not kernel and twin.perm is not kernel.perm
        rng = np.random.default_rng(20260805)
        x = rng.standard_normal((kernel.n, _NG))
        np.testing.assert_array_equal(kernel.apply(x), twin.apply(x))
        np.testing.assert_array_equal(
            kernel.apply_transpose(x), twin.apply_transpose(x)
        )
        np.testing.assert_array_equal(
            kernel.inverse().apply(x), twin.inverse().apply(x)
        )


# ─────────────────────────────────────────────────────────────────────
# The metric — why .H is blind, stated as its own claim
# ─────────────────────────────────────────────────────────────────────


class TestTheMirrorPreservesTheTraceMetric:
    r"""⭐ :math:`G_{\Gamma_-} = G_{\Gamma_+} \circ \pi`, exactly.

    A specular mirror sends an outflow ordinate to the inflow ordinate with
    the *same* :math:`|\Omega\cdot\hat n|` and the *same* :math:`w_n`, so the
    trace metric is permuted rather than changed. This is
    ``vv-principles`` Mode 12's own named example, and it is why ``.H``
    coincides with the bare transpose here.

    ⚠ **Pinned as its own row deliberately.** A gate asserting only
    ``.H ≈ transpose`` would be pinning a *coincidence*: if a future
    quadrature or tangential-band change broke the criterion, that gate would
    red with no indication of why. This row reds FIRST and says which
    assumption moved.
    """

    @pytest.mark.parametrize("case", list(_CASES))
    def test_the_inflow_metric_is_the_permuted_outflow_metric(
        self, case: str
    ) -> None:
        """Bit-exact — no tolerance, because the weights are the same floats."""
        _, trace, _, _, kernel, _, face = _fixture(case)
        outflow_weights = np.asarray(
            trace.outflow_space(face).inner_product_weights, dtype=float
        )
        inflow_weights = np.asarray(
            trace.inflow_space(face).inner_product_weights, dtype=float
        )
        np.testing.assert_array_equal(
            inflow_weights.reshape(-1)[: kernel.n],
            outflow_weights.reshape(-1)[: kernel.n][np.asarray(kernel.perm)],
        )


class TestTheHilbertAdjointCarriesTheMetric:
    r"""``.H`` is :math:`G_+^{-1} P^{\mathsf T} G_-`, and the metric is load-bearing.

    ⛔ **The value leg alone would pin "the metric is ignored."** Because the
    mirror preserves the trace metric (the class above), ``.H`` numerically
    equals the bare transpose on every shipped fixture — so an implementation
    that dropped the metric entirely would pass it. The two perturbation legs
    are what make the gate about the Hilbert adjoint rather than about
    ``np.take``.
    """

    @pytest.mark.parametrize("case", list(_CASES))
    def test_the_adjoint_matches_the_transpose_on_a_metric_preserving_map(
        self, case: str
    ) -> None:
        r"""`[M]` :math:`\le` 1 nulp; tolerance 2, from the round-trip DEPTH.

        The residual is the :math:`(g \cdot y)/g` round trip — one multiply
        and one divide — not physics. ``assert_array_equal`` would be a FALSE
        RED on three of the five fixtures, and which three is an accident of
        the weights: read the 0-vs-1 split as luck, not as meaning.
        """
        _, _, _, _, kernel, _, _ = _fixture(case)
        rng = np.random.default_rng(20260805)
        y = rng.standard_normal((kernel.n, _NG))
        np.testing.assert_array_almost_equal_nulp(
            kernel.H.apply(y), kernel.apply_transpose(y), nulp=2
        )

    @pytest.mark.parametrize("case", list(_CASES))
    @pytest.mark.parametrize(
        "end, expected_ratio", [("domain", 2.0 / 3.0), ("codomain", 2.0)]
    )
    def test_scaling_one_end_metric_moves_the_adjoint_by_a_known_factor(
        self, case: str, end: str, expected_ratio: float
    ) -> None:
        r"""⭐ The negative leg — and it asserts a NUMBER, not "it moved".

        Scaling :math:`G_+` by 3 divides :math:`G_+^{-1}` by 3, so ``.H``
        moves by :math:`1 - 1/3 = 2/3`; scaling :math:`G_-` by 3 multiplies
        by 3, so it moves by :math:`3 - 1 = 2`. Both are exact and
        **quadrature-independent** — a fixture-dependent result is itself the
        signal that something other than the metric changed.
        """
        _, _, _, _, kernel, _, _ = _fixture(case)
        end_space = getattr(kernel, end)
        perturbed = replace(
            end_space,
            inner_product_weights=(
                np.asarray(end_space.inner_product_weights, dtype=float) * 3.0
            ),
        )
        rebuilt = PermutationOperator(
            np.asarray(kernel.perm).copy(),
            axis=kernel.axis,
            **{
                "domain": kernel.domain,
                "codomain": kernel.codomain,
                end: perturbed,
            },
        )
        rng = np.random.default_rng(20260805)
        y = rng.standard_normal((kernel.n, _NG))
        baseline = kernel.H.apply(y)
        moved = rebuilt.H.apply(y)
        ratio = np.max(np.abs(moved - baseline)) / np.max(np.abs(baseline))
        assert ratio == pytest.approx(expected_ratio, rel=1e-12), (
            f"scaling the {end} metric by 3 must move .H by exactly "
            f"{expected_ratio}; got {ratio}. A ratio of 0 means the metric is "
            f"being ignored — which the value leg above cannot detect, "
            f"because this mirror preserves the metric."
        )

    @pytest.mark.parametrize("case", list(_CASES))
    def test_the_weighted_adjoint_law_holds(self, case: str) -> None:
        r""":math:`\langle Px, y\rangle_{\Gamma_-} = \langle x, P^{*}y\rangle_{\Gamma_+}`.

        The defining pairing, asked of the two bound metrics directly — an
        independent route to the same claim as the transpose row, which
        compares arrays rather than inner products.
        """
        _, _, _, _, kernel, _, _ = _fixture(case)
        rng = np.random.default_rng(20260805)
        x = rng.standard_normal((kernel.n, _NG))
        y = rng.standard_normal((kernel.n, _NG))
        assert kernel.domain is not None and kernel.codomain is not None
        left = kernel.codomain.inner_product(kernel.apply(x), y)
        right = kernel.domain.inner_product(x, kernel.H.apply(y))
        assert left == pytest.approx(right, rel=1e-14, abs=1e-15)


# ─────────────────────────────────────────────────────────────────────
# The optional-binding branch, and the hole step 8 must close
# ─────────────────────────────────────────────────────────────────────


def test_a_traceless_method_space_yields_an_unbound_kernel() -> None:
    r"""The optional-binding branch is LIVE, and it is the cleanest control.

    `[M]` ``SNMethodSpace.for_face`` *without* ``trace=`` cannot reach the
    kernel — ``_outflow_restriction`` raises first — but direct dataclass
    construction can, and that path is exactly what the ``if trace is not
    None`` comment claims to serve.

    Two jobs in one row. It pins the documented contract that binding stays
    optional while the tree migrates (#330); and it is the positive control
    for every refusal in this module — same permutation, same face, same code
    path, differing **only** in whether a trace was supplied, and the square
    goes from refused to composable.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    layout = FaceLayout.from_named_shapes(
        [("xmin", (int(quad.N), _NG)), ("xmax", (int(quad.N), _NG))]
    )
    bound_trace = AngularTraceSpace.from_quadrature_and_layout(quad, layout)
    with_trace = SNMethodSpace.for_face(
        quadrature=quad, face="xmin", trace=bound_trace
    )
    traceless = replace(with_trace, trace=None)

    gamma_out = _outflow_restriction(traceless, "reflective")
    kernel = _specular_kernel(
        quad, traceless, gamma_out, axis="x", law_key="reflective"
    )
    assert kernel.domain is None and kernel.codomain is None
    assert gamma_out.domain is None and gamma_out.codomain is None
    # The control: unbound, the square composes. Bound, it raises.
    squared = kernel @ kernel
    assert squared.domain is None


@pytest.mark.xfail(
    strict=True,
    reason=(
        "G6.3 step 8 (#330): TensorProductOperator does not derive its spaces "
        "from its factors, so the binding is enforced at the inner "
        "permutation and invisible at the object the realizer returns. The "
        "committed step-3 Lambertian chain has the identical hole. Until this "
        "flips, routing _reflect_trace through `@` gates nothing, because one "
        "None short-circuits the composability check."
    ),
)
def test_the_realized_operator_carries_the_binding() -> None:
    r"""⛔ The hole step 8 must close, committed as a strict xfail.

    Body kept to the two assertions so **exactly one statement can fail and
    it is the documented one** — the operator is built above the xfail's
    reach would be wrong here, so it is built inside, but every construction
    step is independently exercised by the passing rows in this module.

    Flip-proof `[M]`: teaching ``TensorProductOperator.domain`` to forward
    ``ops[0].domain`` turns this row ``XPASS(strict)``, which is what proves
    it measures the thing step 8 will change rather than being ceremony.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    layout = FaceLayout.from_named_shapes(
        [("xmin", (int(quad.N), _NG)), ("xmax", (int(quad.N), _NG))]
    )
    trace = AngularTraceSpace.from_quadrature_and_layout(quad, layout)
    method_space = SNMethodSpace.for_face(
        quadrature=quad, face="xmin", trace=trace
    )
    realized = SNBoundaryRealizer().realize(
        ReflectiveBoundary(axis="x"), method_space
    )
    assert realized.domain is trace.outflow_space("xmin")
    assert realized.codomain is trace.inflow_space("xmin")
