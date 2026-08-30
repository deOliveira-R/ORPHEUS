r"""The Lambertian chain — the gates migrated off the retired welded operator.

Campaign step **G6.3 step 3b** (#330). ``AngularAverageOperator`` is retired;
its successor is the two-link chain committed in ``b4290873``:

.. math::

    \Gamma_+(f) \;\xrightarrow{\;C\;}\; S(f)
                \;\xrightarrow{\;B\;}\; \Gamma_-(f),
    \qquad R = B \circ C .

This module carries the claims of ``tests/sn/operators/test_angular_average_operator.py``
that survive the retirement, re-posed on the chain. The classification (16
REWIRE / 5 DELETE, one verdict per test with its reason) is
``scratch/g6_3b_test_migration.md``; what is NOT here is there too, with the
gate that subsumes it.

It deliberately does **not** restate what
``tests/sn/operators/test_lambertian_factored.py`` already gates — the
factoring's bit-identity, the composite typing, the wrong-way composition
refusal, the weighted adjoint law and its unbound-mutation leg, and the two
emission guards all live there. The complements this module adds are:

* **faces and axes.** The factored module tests ``xmin`` on three quadratures;
  every claim below runs on ``xmax`` / ``ymin`` / ``zmax`` as well, so the
  ``"{axis}{min|max}"`` outward-normal convention (the ERR-041 class) has a
  unit-level gate on both signs and all three axes.
* **a hand-spelled reference.** Every reference here is written from
  ``quad.weights`` / ``quad.mu_*`` and the face convention spelled out in
  :func:`_axis_and_sign` — never from
  :func:`~orpheus.numerics.spaces.angular_trace_space.build_omega_dot_n`, the
  primitive production reads. The arithmetic (:math:`\mathrm{sign}\cdot\mu`) is
  necessarily the same on both sides; what the hand spelling keeps independent
  is *which face carries which outward normal* and *where the tangential band
  is*.
* **the B3.4a twin classifier.** :meth:`TestTheDomainIsGammaPlus.\
  test_the_domain_is_the_tangential_band_not_the_strict_compare` is the only
  gate in the tree that discriminates ``> 0.0`` from ``> TANGENTIAL_EPS``, and
  the successor's own ``cos_w > 0`` constructor guard is provably **blind** to
  it — the mis-admitted weights are ``6.98e-16``, strictly positive. That
  measurement is pinned inside the gate.
* **the mirror-face adjoint.** :class:`TestReciprocityAgainstTheMirrorFace`
  says *what* :math:`R^{*}` is (the Lambertian kernel of the reversed face, a
  second production object) rather than that ``.H`` is self-consistent with
  itself.

⚠ **Known blindness, MEASURED rather than assumed.** Flipping the outward-normal
convention (``face_normal`` returns ``-sign``) reddens **7** of the 37 rows:
:meth:`TestTheDomainIsGammaPlus.test_the_half_traces_are_the_hand_spelled_ones`
on all five cases, the classifier row, and the ``gl8`` conservation row —
*nothing else*. Every other claim here is **symmetric in the two hemispheres**
(the average of a constant is the constant whichever half you average;
:math:`J^- = J^+` and the mirror-face reciprocity both swap into themselves), so
they cannot see the flip in principle. Credit the convention to
``test_the_half_traces_are_the_hand_spelled_ones`` and to nothing else.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from orpheus.numerics.face_layout import FaceLayout, face_opposite
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces import AngularTraceSpace
from orpheus.numerics.spaces.angular_trace_space import (
    TANGENTIAL_EPS,
    AngularFaceTraceSpace,
)
from orpheus.sn.boundary.angular import (
    IsotropicEmissionOperator,
    PartialCurrentOperator,
)

#: Group extent of every face slot. >1 so the trailing axis is never degenerate.
_NG = 3
_SEED = 20260804

#: ``(id, quadrature factory, face)``. Covers both outward signs, all three
#: axes, the 1-D slab adapter and two 3-D cubature families — the axis/face
#: dimension the factored module (``xmin`` only) does not carry.
_CASES = (
    ("lebedev17/xmax", lambda: Quadrature.lebedev(17), "xmax"),
    ("ls4/ymin", lambda: Quadrature.level_symmetric(4), "ymin"),
    ("gl8/xmax", lambda: Quadrature.gauss_legendre(8), "xmax"),
    ("ls6/zmax", lambda: Quadrature.level_symmetric(6), "zmax"),
    ("lebedev11/zmax", lambda: Quadrature.lebedev(11), "zmax"),
)
_IDS = [case[0] for case in _CASES]

#: Reduction-depth bound, not a fitted number: both currents are
#: ``|Γ₊|``-term positive-summand contractions, so the error is
#: ``|Γ₊| · eps`` — ``49 × 2.2e-16 ≈ 1.1e-14`` on the deepest fixture
#: (``lebedev(17)``). `[M]` worst measured across the five cases: 4.8e-16.
_CURRENT_RTOL = 1e-13

#: ``R.H`` divides by ``G₊`` and multiplies by ``G₋`` where the mirror chain
#: multiplies directly, so the two differ by a ``(g·x)/g`` round trip — a few
#: ULP, never bit-identical. `[M]` worst measured: 1.4e-16 relative.
_ADJOINT_RTOL = 1e-14


# ─────────────────────────────────────────────────────────────────────
# Hand-spelled reference — the face convention written out, never fetched
# ─────────────────────────────────────────────────────────────────────


def _axis_and_sign(face: str) -> tuple[str, int]:
    r"""``"{axis}{min|max}"`` → ``(axis, outward sign)``, spelled by hand.

    ``…max`` is the upper face (outward normal :math:`+\hat e_{\rm axis}`),
    ``…min`` the lower. This is the convention the production renderer
    :func:`~orpheus.numerics.face_layout.face_name` implements; writing it
    out here is what keeps the references below independent of it, so a flip
    reddens rather than cancels (the ERR-041 threat class).
    """
    axis, tag = face[:1], face[1:]
    if tag == "max":
        return axis, +1
    if tag == "min":
        return axis, -1
    raise ValueError(f"not a '{{axis}}{{min|max}}' face name: {face!r}")


def _signed_projection(quad: Quadrature, face: str) -> np.ndarray:
    r""":math:`\Omega\cdot\hat n` at ``face``, from the raw cosines."""
    axis, sign = _axis_and_sign(face)
    mu = {"x": quad.mu_x, "y": quad.mu_y, "z": quad.mu_z}[axis]
    return float(sign) * np.asarray(mu, dtype=float)


def _half_traces(quad: Quadrature, face: str) -> tuple[np.ndarray, np.ndarray]:
    r"""``(Γ₊, Γ₋)`` index arrays for ``face``, by the tangential-band rule.

    The face partition is THREE-way — :math:`\Gamma_+ \sqcup \Gamma_- \sqcup`
    tangential — so "not inflow" is never spelled as "outflow" here.
    """
    odn = _signed_projection(quad, face)
    return (
        np.flatnonzero(odn > +TANGENTIAL_EPS),
        np.flatnonzero(odn < -TANGENTIAL_EPS),
    )


def _cosine_weights(quad: Quadrature, rows: np.ndarray, face: str) -> np.ndarray:
    r"""``w_n |Ω·n̂|`` restricted to ``rows`` — the partial-current weight."""
    odn = _signed_projection(quad, face)
    return np.asarray(quad.weights, dtype=float)[rows] * np.abs(odn[rows])


# ─────────────────────────────────────────────────────────────────────
# The system under test — production's own construction
# ─────────────────────────────────────────────────────────────────────


class _Lambertian:
    r"""The bound chain :math:`R = B \circ C` on ONE face, plus its pieces.

    ⚠ The construction mirrors
    :func:`~orpheus.sn.boundary.realizer._checked_angular_average` (minus its
    ERR-041 orientation cross-check, which is the realizer's own claim and is
    gated in ``test_sn_boundary_realizer.py``). Calling the realizer here
    instead would make a unit module depend on the dispatch layer above it. If
    this file is folded into ``test_lambertian_factored.py``, this class and
    that module's ``_Lambertian`` collapse into one — see
    ``scratch/g6_3b_test_migration.md`` §6.
    """

    def __init__(self, trace: AngularTraceSpace, face: str) -> None:
        self.trace = trace
        self.face = face
        self.gamma_plus = trace.outflow_space(face)
        self.gamma_minus = trace.inflow_space(face)
        self.current = trace.current_space(face)

        face_metric = np.asarray(trace.face_space(face).inner_product_weights)
        self.cos_w = face_metric[trace.outflow_indices_for_face(face)]
        self.norm = float(self.cos_w.sum())

        self.C = PartialCurrentOperator(
            self.cos_w, domain=self.gamma_plus, codomain=self.current,
        )
        self.B = IsotropicEmissionOperator(
            self.norm, self.gamma_minus.shape[0],
            domain=self.current, codomain=self.gamma_minus,
        )
        self.R = self.B @ self.C


def _trace(quad: Quadrature, faces) -> AngularTraceSpace:
    """The production trace space over ``faces``, one slot of ``_NG`` groups."""
    layout = FaceLayout.from_named_shapes(
        [(f, (int(quad.N), _NG)) for f in faces]
    )
    return AngularTraceSpace.from_quadrature_and_layout(quad, layout)


def _one_face(quad: Quadrature, face: str) -> _Lambertian:
    return _Lambertian(_trace(quad, (face,)), face)


def _mirror_pair(quad: Quadrature, face: str) -> tuple[_Lambertian, _Lambertian]:
    """This face's chain and the reversed face's, from ONE trace space."""
    trace = _trace(quad, sorted({face, face_opposite(face)}))
    return _Lambertian(trace, face), _Lambertian(trace, face_opposite(face))


# ─────────────────────────────────────────────────────────────────────
# A — the domain IS Γ₊ (L0)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestTheDomainIsGammaPlus:
    r"""Which index set the chain reads, and against which threshold.

    MIGRATED from ``TestOutflowHalfTraceIsTheDomain``. Pre-B3.4a those rows
    asserted that a masked ``cos_w`` *erased* the incoming contribution; the
    narrowing made the erasure unspellable (the incoming ordinates are not in
    the domain) and the retirement leaves the claim exactly where it belongs —
    on the trace space that does the classifying and on the chain bound to it.
    """

    @pytest.mark.parametrize("quad_factory,face", [c[1:] for c in _CASES], ids=_IDS)
    def test_a_unit_outflow_re_emits_unit_intensity(self, quad_factory, face):
        r""":math:`R\,\mathbf 1 = \mathbf 1` — the normalisation is the
        :math:`\Gamma_+` sum.

        The surviving, honest form of the old "only the outgoing hemisphere
        contributes" row. It is not decoration: the normalisation lives in
        :class:`IsotropicEmissionOperator` and is the ONLY thing making the
        identity hold, so a norm taken over the full face slot (the pre-B3.4a
        shape) gives ``norm₊ / norm_full ≠ 1`` — an O(1) miss, on the cheapest
        row in the file.

        The probe is the old one restricted to :math:`\Gamma_+`: ``+1`` on the
        outgoing hemisphere, ``-1`` elsewhere, sliced by the HAND-SPELLED
        :math:`\Gamma_+`, so a production classification of a DIFFERENT SIZE
        fails on shape before it can fail on value.

        ⚠ It is blind to a same-size mis-classification — notably a flipped
        outward normal, `[M]` which leaves this row green because the average of
        a constant is that constant whichever hemisphere you average over. That
        is :meth:`test_the_half_traces_are_the_hand_spelled_ones`'s job, and
        only its.
        """
        quad = quad_factory()
        lam = _one_face(quad, face)
        outflow, inflow = _half_traces(quad, face)
        sign_field = np.where(_signed_projection(quad, face) > 0.0, +1.0, -1.0)
        psi = np.repeat(sign_field[:, None], _NG, axis=1)[outflow]

        result = np.asarray(lam.R.apply(psi))
        assert result.shape == (inflow.size, _NG), (
            f"emitted {result.shape}; Γ₋ on {face} is {inflow.size} ordinates "
            f"over {_NG} groups."
        )
        np.testing.assert_allclose(result, +1.0, rtol=1e-14)

    @pytest.mark.parametrize("quad_factory,face", [c[1:] for c in _CASES], ids=_IDS)
    def test_the_half_traces_are_the_hand_spelled_ones(self, quad_factory, face):
        r"""The production half-traces ARE the hand-spelled ones, and both sit
        strictly inside the face slot.

        The convention gate (ERR-041 class): ``…max`` must carry outward normal
        :math:`+\hat e_{\rm axis}`. A flipped renderer swaps the two index sets
        and this row is the first to say so, on all three axes and both signs.

        It also SUBSUMES the retired ``test_z_on_lebedev_succeeds``: the two
        ``zmax`` cases assert that a ``z`` face is realizable on a cubature with
        genuine :math:`\mu_z`, and the strict ``< N`` legs are what a revert to
        the full-face domain would break.
        """
        quad = quad_factory()
        lam = _one_face(quad, face)
        outflow, inflow = _half_traces(quad, face)

        np.testing.assert_array_equal(
            lam.trace.outflow_indices_for_face(face), outflow,
            err_msg=f"Γ₊({face}) is not the hand-spelled sign·μ > +ε set.",
        )
        np.testing.assert_array_equal(
            lam.trace.inflow_indices_for_face(face), inflow,
            err_msg=f"Γ₋({face}) is not the hand-spelled sign·μ < −ε set.",
        )
        assert lam.gamma_plus.shape[0] == outflow.size
        assert lam.gamma_minus.shape[0] == inflow.size
        assert 0 < outflow.size < quad.N and 0 < inflow.size < quad.N, (
            f"half-traces {outflow.size}/{inflow.size} are not strictly inside "
            f"the {quad.N}-ordinate face slot — a full-N domain would be the "
            f"retired welded shape."
        )
        # The endpoints know WHICH face they are: `FunctionSpace.__eq__` is
        # (name, shape), and every quadrature gives |Γ₊(xmin)| == |Γ₊(xmax)|,
        # so the face in the name is the only thing refusing a cross-face
        # composition.
        domain, codomain = lam.R.domain, lam.R.codomain
        assert isinstance(domain, AngularFaceTraceSpace)
        assert isinstance(codomain, AngularFaceTraceSpace)
        assert (domain.face, domain.role) == (face, "outflow")
        assert (codomain.face, codomain.role) == (face, "inflow")

    def test_the_domain_is_the_tangential_band_not_the_strict_compare(self):
        r"""The domain is :math:`\Gamma_+` at ``TANGENTIAL_EPS`` — NOT the
        retired private ``(sign·μ) > 0.0`` set.

        ⭐ **The load-bearing row of the migration.** It is the only gate in the
        tree that discriminates the B3.4a twin classifier, and the successor's
        own construction guard cannot replace it: on this fixture the
        mis-admitted rows carry ``cos_w ≈ 7e-16`` — **strictly positive** — so
        :class:`PartialCurrentOperator`'s ``cos_w > 0`` refusal (gated next door
        by ``test_the_contraction_refuses_a_full_face_weight_vector``) is blind
        to them. That measurement is asserted below rather than trusted, so if
        the guard ever does become sufficient this row says so.

        What catches the twin after the retirement is the **domain size**: the
        trace space classifies at the band, so the chain reads 2 rows where the
        strict compare would read 6, and a caller still classifying the old way
        is refused at apply.

        ⚠ **The fixture is CONSTRUCTED, and why is the point.** It used to be
        the shipped ``product(2, 4)``, which discriminated because its
        :math:`\varphi = \pm\pi/2` azimuths came out of
        ``np.cos(np.linspace(...))`` as :math:`\pm 6.1\text{e-}17` — i.e. the
        gate was activated by round-off. Once the azimuths became roots of unity
        (E3/E4, issue #325) those ordinates are *exactly* tangential and no
        shipped rule discriminates any more. The nudge goes into the MEASURE, so
        the trace space and the chain are both built from it and the claim is
        about production, not about a hand-made comparison array.

        The first assertion is the ACTIVATION guard — without it this row
        silently decays into a restatement of its siblings (``vv`` Mode-8
        class 7). It is what caught that decay on 2026-08-02.
        """
        shipped = Quadrature.product(n_mu=2, n_phi=4)
        nodes = np.array(shipped.measure.nodes, copy=True)
        nodes[nodes[:, 0] == 0.0, 0] = TANGENTIAL_EPS / 2.0
        quad = replace(shipped, measure=replace(shipped.measure, nodes=nodes))

        odn = _signed_projection(quad, "xmax")
        gamma_plus = np.flatnonzero(odn > +TANGENTIAL_EPS)
        retired_strict = np.flatnonzero(odn > 0.0)
        assert not np.array_equal(gamma_plus, retired_strict), (
            f"fixture no longer discriminates: Γ₊ {gamma_plus.tolist()} and the "
            f"retired '> 0.0' set {retired_strict.tolist()} agree, so this gate "
            f"cannot see the classifier the narrowing retired."
        )

        trace = _trace(quad, ("xmin", "xmax"))
        np.testing.assert_array_equal(
            trace.outflow_indices_for_face("xmax"), gamma_plus,
            err_msg="the trace space is not classifying at TANGENTIAL_EPS.",
        )
        lam = _Lambertian(trace, "xmax")
        assert lam.gamma_plus.shape[0] == gamma_plus.size, (
            f"the chain reads {lam.gamma_plus.shape[0]} ordinates; Γ₊ "
            f"(|Ω·n̂| > TANGENTIAL_EPS) has {gamma_plus.size}. A domain of "
            f"{retired_strict.size} would be the retired '> 0.0' classifier."
        )

        # `[M]` the constructor guard is BLIND here — pin it, so a future
        # strengthening of the guard is visible rather than assumed.
        mis_admitted = np.setdiff1d(retired_strict, gamma_plus)
        weights = np.asarray(
            trace.face_space("xmax").inner_product_weights
        )[mis_admitted]
        assert mis_admitted.size and np.all(weights > 0.0), (
            f"the mis-admitted rows now carry non-positive weights "
            f"({weights!r}), so PartialCurrentOperator's `cos_w > 0` guard "
            f"would refuse them — this gate's premise has changed."
        )

        # …and a caller that still classifies the old way is refused, rather
        # than silently averaging six rows against two weights.
        with pytest.raises(ValueError, match=r"expected \|Γ₊\|"):
            lam.R.apply(np.ones((retired_strict.size, _NG)))


# ─────────────────────────────────────────────────────────────────────
# B — cosine-weighted current conservation (L1)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestCosineWeightedCurrentConservation:
    r"""The current re-emitted on :math:`\Gamma_-` equals the one received on
    :math:`\Gamma_+`.

    .. math::

        J_- \;=\; \sum_{\Omega \in \Gamma_-} w\,|\Omega\cdot\hat n|\,
                  (R\psi)(\Omega)
        \;=\; \sum_{\Omega \in \Gamma_+} w\,|\Omega\cdot\hat n|\,\psi(\Omega)
        \;=\; J_+ ,

    i.e. a white face at :math:`\alpha = 1` returns exactly the partial current
    it received (Bell & Glasstone 1970 §1.5).

    ⭐ **Stronger than its counterpart next door.**
    ``test_the_chain_is_CONSERVATIVE_on_a_symmetric_quadrature`` drives
    :math:`\psi = \mathbf 1`, at which point the law collapses into its own
    activation guard (``cos_w_in.sum() == norm``) and the numerator carries no
    information. These rows drive a random positive :math:`\psi` and a strongly
    anisotropic one, so the cosine WEIGHTING is exercised and not merely the
    weight sum.

    Both reference sums are hand-spelled, so a face-convention flip reddens
    them; the quadrature's face symmetry is asserted inline as an activation
    guard rather than assumed (`[M]` it holds bit-identically on every
    quadrature × face pair in the tree, and a future asymmetric rule would make
    this operator genuinely non-conservative — which must red HERE, loudly,
    rather than be absorbed into a tolerance).
    """

    @staticmethod
    def _weights(quad, face):
        outflow, inflow = _half_traces(quad, face)
        w_plus = _cosine_weights(quad, outflow, face)
        w_minus = _cosine_weights(quad, inflow, face)
        np.testing.assert_allclose(
            w_plus.sum(), w_minus.sum(), rtol=1e-15,
            err_msg=(
                f"quadrature is not face-symmetric at {face}: the Γ₊ and Γ₋ "
                f"cosine-weighted weight sums differ, so a Lambertian face "
                f"cannot return the current it received."
            ),
        )
        return w_plus, w_minus

    @pytest.mark.parametrize("quad_factory,face", [c[1:] for c in _CASES], ids=_IDS)
    def test_the_reemitted_current_equals_the_received_current(
        self, quad_factory, face,
    ):
        """Random positive outflow — the numerator carries information."""
        quad = quad_factory()
        lam = _one_face(quad, face)
        w_plus, w_minus = self._weights(quad, face)

        rng = np.random.default_rng(_SEED)
        psi = rng.uniform(0.5, 1.5, size=lam.gamma_plus.shape)
        psi_in = np.asarray(lam.R.apply(psi))

        j_plus = (w_plus[:, None] * psi).sum(axis=0)
        j_minus = (w_minus[:, None] * psi_in).sum(axis=0)
        np.testing.assert_allclose(
            j_minus, j_plus, rtol=_CURRENT_RTOL,
            err_msg=f"J⁻ = {j_minus!r} vs J⁺ = {j_plus!r} on {face}",
        )

    def test_the_reemitted_current_with_a_strongly_anisotropic_input(self):
        r"""``exp(μ_z)·(μ_z + 2)``, per-group scaled — the one row where the
        cosine weighting is the dominant term rather than the weight sum.

        Also the only ``z``-axis conservation row, and the only one whose
        trailing (group) axis is non-degenerate, so a reduction that collapsed
        the group axis would be caught here rather than pass by symmetry.
        """
        quad = Quadrature.level_symmetric(6)
        face = "zmax"
        lam = _one_face(quad, face)
        outflow, _ = _half_traces(quad, face)
        w_plus, w_minus = self._weights(quad, face)

        shape_field = np.exp(quad.mu_z) * (quad.mu_z + 2.0)
        psi = shape_field[outflow][:, None] * np.array([1.0, 3.0, 7.0])
        psi_in = np.asarray(lam.R.apply(psi))

        j_plus = (w_plus[:, None] * psi).sum(axis=0)
        j_minus = (w_minus[:, None] * psi_in).sum(axis=0)
        np.testing.assert_allclose(
            j_minus, j_plus, rtol=_CURRENT_RTOL,
            err_msg=f"J⁻ = {j_minus!r} vs J⁺ = {j_plus!r}",
        )


# ─────────────────────────────────────────────────────────────────────
# C — reciprocity, against the mirror face's own chain (L1)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestReciprocityAgainstTheMirrorFace:
    r"""The weighted adjoint of the Lambertian IS the Lambertian of the
    REVERSED face.

    With :math:`W_\pm = \operatorname{diag}(w\,|\Omega\cdot\hat n|)` on each
    half-trace,

    .. math::

        R^{*} \;=\; W_+^{-1} R^{\mathsf T} W_-
        \;=\; \mathbf 1_{\Gamma_+} \otimes
              \frac{w\,|\Omega\cdot\hat n|}{\operatorname{norm}_+}
              \bigg|_{\Gamma_-} ,

    which is the same form with the two half-traces exchanged — i.e. exactly
    the chain built on ``face_opposite(f)``, whose :math:`\Gamma_+` IS this
    face's :math:`\Gamma_-`.

    ⭐ **Why this is not a restatement of the flagship gate next door.**
    ``test_the_WEIGHTED_adjoint_law_holds`` verifies that ``.H`` is
    *self-consistent* with the pairing; the first row here never touches ``.H``
    at all — its reference is a second production chain built from the OTHER
    face's data — so it is the one adjoint gate in the tree that is independent
    of the ``AdjointOperator`` machinery. The second row then closes the loop
    by measuring ``.H`` against it pointwise.
    """

    @staticmethod
    def _activate(lam, mirror):
        """The two preconditions, asserted rather than assumed."""
        face, opp = lam.face, mirror.face
        np.testing.assert_array_equal(
            lam.trace.outflow_indices_for_face(opp),
            lam.trace.inflow_indices_for_face(face),
            err_msg=f"Γ₊({opp}) is not Γ₋({face}); the mirror is not the adjoint's shape.",
        )
        np.testing.assert_array_equal(
            lam.trace.inflow_indices_for_face(opp),
            lam.trace.outflow_indices_for_face(face),
        )
        np.testing.assert_allclose(
            mirror.norm, lam.norm, rtol=1e-15,
            err_msg="the two faces' normalisations differ; reciprocity is not posed.",
        )
        # The two spaces carry the SAME ordinates and are still not equal —
        # `FunctionSpace.__eq__` is (name, shape) and the name carries the
        # face. This is the reciprocity's cross-face content in one line.
        assert lam.gamma_minus != mirror.gamma_plus, (
            "Γ₋(f) compares EQUAL to Γ₊(f_opposite); a cross-face composition "
            "would then type-check while being wrong."
        )

    @pytest.mark.parametrize("quad_factory,face", [c[1:] for c in _CASES], ids=_IDS)
    def test_the_weighted_adjoint_is_the_mirror_face_kernel(
        self, quad_factory, face,
    ):
        r""":math:`\langle R\psi,\varphi\rangle_{\Gamma_-} =
        \langle \psi, R'\varphi\rangle_{\Gamma_+}` with :math:`R'` the reversed
        face's chain — no ``.H`` anywhere.
        """
        lam, mirror = _mirror_pair(quad_factory(), face)
        self._activate(lam, mirror)

        rng = np.random.default_rng(_SEED + 1)
        psi = rng.uniform(0.0, 1.0, size=lam.gamma_plus.shape)
        phi = rng.uniform(0.0, 1.0, size=lam.gamma_minus.shape)

        lhs = lam.gamma_minus.inner_product(lam.R.apply(psi), phi)
        rhs = lam.gamma_plus.inner_product(psi, mirror.R.apply(phi))
        assert lhs == pytest.approx(rhs, rel=_CURRENT_RTOL), (
            f"⟨Rψ,φ⟩_W₋ = {lhs!r} vs ⟨ψ,R'φ⟩_W₊ = {rhs!r}"
        )

    @pytest.mark.verifies("bc-response-factored-adjoint")
    @pytest.mark.parametrize("quad_factory,face", [c[1:] for c in _CASES], ids=_IDS)
    def test_H_is_pointwise_the_mirror_face_kernel(self, quad_factory, face):
        r"""``R.H`` and the reversed-face chain agree ELEMENTWISE.

        This is :eq:`bc-response-factored-adjoint`
        (:math:`R^{*} = G_+^{-1}R^{\mathsf T}G_-`) measured against an
        independently constructed production object, which is a different
        structural angle from both of its siblings: the dense-basis gate
        verifies :math:`R^{\mathsf T}` (Euclidean, not :math:`R^{*}`), and the
        weighted-law gate verifies :math:`R^{*}` through a pairing that already
        uses ``.H``. This row says what :math:`R^{*}` IS.

        `[M]` agreement is 0.0 – 1.4e-16 relative over the five cases. It is not
        bit-identical by construction: ``.H`` divides by :math:`G_+` where the
        mirror multiplies, a ``(g·x)/g`` round trip worth a couple of ULP.
        Dropping the metric from ``.H`` leaves the Euclidean
        :math:`(\cos w_+/\mathrm{norm}) \otimes \mathbf 1`, which is not
        proportional to this at all — an O(1) miss.
        """
        lam, mirror = _mirror_pair(quad_factory(), face)
        self._activate(lam, mirror)

        rng = np.random.default_rng(_SEED + 2)
        phi = rng.standard_normal(lam.gamma_minus.shape)

        expected = np.asarray(mirror.R.apply(phi))
        assert np.abs(expected).max() > 1e-3, (
            "the mirror kernel returned ~0 on this draw; the comparison below "
            "would be vacuous."
        )
        np.testing.assert_allclose(
            np.asarray(lam.R.H.apply(phi)), expected, rtol=_ADJOINT_RTOL,
        )


# ─────────────────────────────────────────────────────────────────────
# D — capability predicates (L0)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestCapabilityPredicates:
    def test_the_chain_is_adjointable_but_NOT_invertible(self):
        r"""INVERTED at G6.3 step 3, and half of it is gated nowhere else.

        The retired ``TestPredicates::test_apply_only`` asserted
        ``not is_adjointable and not is_invertible``. The first half now
        asserts the OPPOSITE of the contract — factoring supplied the transpose
        the welded operator withheld — and its replacement
        (``test_the_welded_operator_withheld_the_transpose_the_chain_supplies``)
        dies with the welded operator it names.

        The second half was never restated anywhere: :math:`R` is **rank one**,
        so it is not invertible *despite being square-shaped* — and it is
        square-shaped on every reachable fixture, which is exactly why the
        claim needs an assertion rather than an inference.
        """
        lam = _one_face(Quadrature.lebedev(11), "zmax")
        assert lam.gamma_plus.shape == lam.gamma_minus.shape, (
            "activation: |Γ₊| != |Γ₋| here, so 'not invertible' would follow "
            "from the shape alone and the rank-one content would be untested."
        )
        assert lam.R.is_adjointable is True
        assert lam.R.is_invertible is False
        assert lam.C.is_adjointable is True and lam.B.is_adjointable is True


# ─────────────────────────────────────────────────────────────────────
# E — face realizability (L0)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestFaceRealizability:
    def test_a_z_face_on_a_slab_quadrature_is_a_rank_mismatch(self):
        r"""1-D Gauss-Legendre carries ``mu_z == zeros(N)``, so a ``z`` face is
        a RANK MISMATCH between the face layout and the cubature.

        RE-POSED onto its honest home. The retired row put this on
        ``AngularAverageOperator.from_quadrature``; the raise never came from
        the operator — it comes from
        :func:`~orpheus.numerics.spaces.angular_trace_space.build_omega_dot_n`,
        the single face-name → signed-projection primitive, reached here
        through the trace-space constructor that production uses. `[M]` message
        unchanged by the retirement.

        The positive control is the point of the second half: the SAME
        quadrature realizes its ``x`` faces, so the refusal is about the axis
        and not about the quadrature being unusable.
        """
        quad = Quadrature.gauss_legendre(8)
        with pytest.raises(ValueError, match="requires genuine mu_z"):
            _trace(quad, ("zmax",))
        # Positive control — an x face on the same slab rule is fine.
        lam = _one_face(quad, "xmax")
        assert lam.gamma_plus.shape[0] > 0


# ─────────────────────────────────────────────────────────────────────
# F — construction guards (L0)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestContractionGuards:
    r"""``cos_w`` must be the restriction to :math:`\Gamma_+`.

    MIGRATED from ``TestConstructionGuards``. The emission's two guards
    (non-positive ``norm``, empty codomain) are gated verbatim next door
    (``test_the_emission_refuses_a_degenerate_normalisation`` /
    ``..._an_empty_codomain``) and are NOT restated here. What is restated is
    the contraction's, because the version next door feeds it a synthetic
    ``[0.5, 0.0, 0.3]`` with no positive control — and `vv-principles`
    anti-#11: a guard tested only negatively validates the raising, not the
    invariant.
    """

    def test_a_real_masked_full_face_weight_vector_is_refused(self):
        """The pre-B3.4a shape, from a real quadrature, plus its control.

        A masked full-face ``cos_w`` — ``N`` entries zeroed off the outgoing
        hemisphere by a second classifier — is the defect the guard exists to
        name, so that is the input the negative feeds it, and the RESTRICTION
        of the same array is the positive control that the guard is not simply
        refusing everything.
        """
        quad = Quadrature.lebedev(17)
        face = "xmax"
        outflow, _ = _half_traces(quad, face)
        odn = _signed_projection(quad, face)
        masked_full = np.where(odn > 0.0, np.asarray(quad.weights) * odn, 0.0)
        assert (masked_full == 0.0).any(), (
            "activation: the probe carries no zeros, so it is not the "
            "masked full-face shape the guard is aimed at."
        )

        with pytest.raises(ValueError, match="STRICTLY positive"):
            PartialCurrentOperator(masked_full)

        # Positive control: the restriction constructs cleanly and computes.
        restricted = PartialCurrentOperator(masked_full[outflow])
        assert restricted.apply(np.ones(outflow.size)) > 0.0

    @pytest.mark.parametrize(
        "bad,match",
        [
            (np.zeros((3, 2)), "must be 1-D"),
            (np.array([]), "empty"),
        ],
        ids=["two-dimensional", "empty"],
    )
    def test_a_malformed_weight_vector_is_refused(self, bad, match):
        """NET-NEW at the migration — the welded operator's ``ndim`` guard had
        no negative test and the successor's ``size == 0`` guard is new.

        Both are the same claim as the row above (``cos_w`` is the
        :math:`\\Gamma_+` restriction, so it is a non-empty 1-D vector); they
        are here because a retirement is the moment the successor's guards get
        the coverage the predecessor's never had.
        """
        with pytest.raises(ValueError, match=match):
            PartialCurrentOperator(bad)


# ─────────────────────────────────────────────────────────────────────
# G — the operator is decoupled from its inputs and outputs (L0)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestDecoupling:
    def test_mutating_the_source_weights_does_not_move_the_operator(self):
        """MIGRATED and strengthened from ``test_quadrature_reference_not_held``.

        The retired row asserted ``op._cos_w is not quad.weights`` — an
        identity check on a private attribute, which a shallow copy would also
        satisfy. The behavioural statement is that the operator's OUTPUT does
        not move when the caller's array does, and the chain adds a guarantee
        the welded operator did not have: the exposed view is read-only, so a
        consumer cannot reach in either.

        (The retired row's two ``hasattr(op, "quadrature")`` legs are dropped —
        API-smoke on attributes that never existed.)
        """
        source = np.array([1.0, 2.0, 3.0])
        contraction = PartialCurrentOperator(source)
        before = np.asarray(contraction.apply(np.ones(3)))

        source[0] = 1000.0
        np.testing.assert_array_equal(
            np.asarray(contraction.apply(np.ones(3))), before,
            err_msg="the operator aliased the caller's weight array.",
        )
        assert contraction.cos_w is not source
        with pytest.raises(ValueError, match="read-only"):
            contraction.cos_w[0] = 5.0

    def test_the_output_is_a_fresh_writeable_array(self):
        """Calling code may mutate the output without affecting later calls.

        Not cosmetic: :meth:`IsotropicEmissionOperator.apply` ends in
        ``np.broadcast_to(...).copy()``. Drop the ``copy()`` and every caller
        gets a NON-WRITEABLE stride-0 view of one scalar — this row reds on the
        write, before it can red on the value.
        """
        lam = _one_face(Quadrature.lebedev(5), "xmax")
        psi = np.ones(lam.gamma_plus.shape)

        first = np.asarray(lam.R.apply(psi))
        assert first.flags.writeable, "the emitted array is a read-only view."
        first[0] = 999.0
        second = np.asarray(lam.R.apply(psi))
        assert not np.any(second == 999.0), (
            "a later apply saw the caller's mutation — the emission is "
            "handing out its own state."
        )


# ─────────────────────────────────────────────────────────────────────
# H — the shape contract (L0)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestShapeContract:
    def test_a_full_face_input_is_refused(self):
        r"""The whole face slot is NOT the domain.

        The un-migrated-caller threat: the shapes are plausible, the values
        would be silently wrong, and the guard is the only thing between them.
        Both wrong sizes land on the same message, so they are two inputs to
        one gate rather than two gates.
        """
        quad = Quadrature.lebedev(5)  # N = 14
        lam = _one_face(quad, "xmax")
        n_plus = lam.gamma_plus.shape[0]
        assert n_plus < quad.N  # activation: the two shapes differ

        for n_rows in (int(quad.N), n_plus + 1):
            with pytest.raises(ValueError, match=r"expected \|Γ₊\|"):
                lam.R.apply(np.ones((n_rows, _NG)))

    def test_trailing_axes_pass_through_and_the_emission_is_isotropic(self):
        r"""``(|Γ₊|, 5, 3) → (|Γ₋|, 5, 3)``, every leading-axis slice equal.

        The row count is asserted against the CODOMAIN, never against
        ``psi.shape``: `[M]` ``|Γ₊| == |Γ₋|`` on every quadrature × face pair in
        the tree, so ``result.shape == psi.shape`` sits inside the shape
        functional's invariance group (``vv`` Mode 12) and would stay green
        under a codomain that echoed the domain.
        :meth:`test_the_codomain_size_follows_n_inflow_not_the_input` breaks
        that degeneracy directly.
        """
        lam = _one_face(Quadrature.lebedev(11), "xmax")
        rng = np.random.default_rng(_SEED + 3)
        psi = rng.uniform(0.0, 1.0, size=(lam.gamma_plus.shape[0], 5, 3))

        result = np.asarray(lam.R.apply(psi))
        assert result.shape == (lam.gamma_minus.shape[0], 5, 3)
        for row in range(result.shape[0]):
            np.testing.assert_allclose(result[row], result[0], rtol=1e-14)

    def test_the_codomain_size_follows_n_inflow_not_the_input(self):
        r"""The number of emitted rows comes from the CODOMAIN datum.

        The anti-Mode-12 leg of the shape contract, and the only row in either
        module that breaks the ``|Γ₊| == |Γ₋|`` coincidence: no
        quadrature-built chain can distinguish "emits :math:`|\Gamma_-|` rows"
        from "echoes its input's leading axis", so the halves have to be made
        unequal by hand.
        """
        chain = IsotropicEmissionOperator(6.0, 5) @ PartialCurrentOperator(
            np.array([1.0, 2.0, 3.0])
        )
        out = np.asarray(chain.apply(np.array([1.0, 1.0, 1.0])))
        assert out.shape == (5,), (
            f"emitted {out.shape[0]} rows for a 3-row input with |Γ₋| = 5 — "
            f"the codomain is echoing the domain."
        )
        np.testing.assert_allclose(out, 1.0, rtol=1e-14)
