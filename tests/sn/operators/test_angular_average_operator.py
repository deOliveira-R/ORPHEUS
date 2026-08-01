r"""Tests for AngularAverageOperator (Wave 1 / C1.1).

The operator implements the cosine-weighted Lambertian average over the
outgoing hemisphere; the incoming flux at a white surface equals
α × (this average).

**It is the RESPONSE, not a geometry.** Campaign phase B3.0 corrected
the older reading — "the §15.2 geometric projection G_diff in the
decomposition R_white = G_diff ⊗ α" — on a decidable criterion: a
geometry map is the composition operator of a measure-preserving
bijection, hence **multiplicative** (``G(ψφ) = (Gψ)(Gφ)``), and an
average is never that. So white's factors are ``G = IdentityMap``
(a white face fixes no geometry) and
``R = LambertianReemission(α)``, which is what this operator realizes.

The misreading was invisible for a structural reason worth keeping in
mind here: a rank-one response **annihilates** ``G`` entirely
(``R∘G = u ⊗ (Gᵀv)`` and the Lambertian's ``v = |Ω·n|`` is invariant
under both the mirror and the periodic translation), so white's
geometry slot had no observable consequence and the physics drifted
into it. Nothing this file measures changed when the factors moved.

Migrated at **B3.4a** — the operator is now :math:`\Gamma_+ \to \Gamma_-`
====================================================================

Pre-B3.4a this was a full-face endomorphism: it took all ``N`` ordinates,
zeroed the non-outgoing ones through a masked ``cos_w``, and broadcast the
average back over all ``N``. Every gate below was written against that
shape. The narrowing changes what each one can say:

* the **mask** gates ("non-outgoing ordinates contribute zero") no longer
  describe an erasure — the non-outgoing ordinates are not in the domain, so
  the claim is now a statement about WHICH index set the domain is. That is
  strictly stronger and it is what
  :class:`TestOutflowHalfTraceIsTheDomain` asserts, with the
  ``product(2, 4)`` row pinning that the domain is the trace space's
  :math:`\Gamma_+` (classified against ``TANGENTIAL_EPS``) and NOT the
  retired private ``> 0.0`` set;
* the **conservation** gate was ``J_out(input) == J_out(broadcast output)``
  using the OUTGOING weights on both sides — self-consistent by
  construction, with teeth only against a wrong ``norm``. Narrowed, the
  honest statement is the physical one: the cosine-weighted current the
  operator RE-EMITS on :math:`\Gamma_-` equals the one it RECEIVED on
  :math:`\Gamma_+`. It subsumes the old claim;
* the **self-adjointness** gate cannot be posed verbatim — the operator is
  rectangular, so ``⟨Ax, y⟩_w == ⟨x, Ay⟩_w`` has no reading with both
  arguments in the right spaces. Its mathematical content survives as
  RECIPROCITY, and the adjoint turns out to be an object the tree already
  builds: see :class:`TestReciprocityCosineWeighted`.

Structural independence of the references
-----------------------------------------

Every reference below is spelled from ``quad.weights`` / ``quad.mu_*`` and
the ``"{axis}{min|max}"`` face convention written out by hand — never from
:func:`~orpheus.numerics.spaces.angular_trace_space.build_omega_dot_n`,
which is the primitive production reads. The arithmetic
(:math:`\mathrm{sign}\cdot\mu`) is necessarily the same on both sides; what
the hand spelling keeps independent is the CONVENTION (which face has which
outward normal) — the ERR-041 threat class — and the CLASSIFICATION
threshold, which is the B3.4a bug.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.sn.boundary.angular import AngularAverageOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces.angular_trace_space import TANGENTIAL_EPS


# ─────────────────────────────────────────────────────────────────────
# Hand-spelled half-trace reference (see module docstring)
# ─────────────────────────────────────────────────────────────────────


def _signed_projection(quad, axis: str, outward_sign: int) -> np.ndarray:
    r""":math:`\Omega\cdot\hat n` for the face ``(axis, outward_sign)`` names.

    The face convention is written out here rather than fetched: ``+1`` is the
    upper face (outward normal :math:`+\hat e_{\rm axis}`), ``-1`` the lower.
    """
    mu = {"x": quad.mu_x, "y": quad.mu_y, "z": quad.mu_z}[axis]
    return float(outward_sign) * np.asarray(mu, dtype=float)


def _half_traces(quad, axis: str, outward_sign: int):
    r"""``(Γ₊, Γ₋)`` index arrays for the face, by the tangential-band rule.

    The face partition is THREE-way — :math:`\Gamma_+ \sqcup \Gamma_- \sqcup`
    tangential — so "not inflow" is never spelled as "outflow" here.
    """
    odn = _signed_projection(quad, axis, outward_sign)
    return (
        np.flatnonzero(odn > +TANGENTIAL_EPS),
        np.flatnonzero(odn < -TANGENTIAL_EPS),
    )


def _cosine_weights(quad, rows: np.ndarray, axis: str, outward_sign: int):
    r"""``w_n |Ω·n̂|`` restricted to ``rows`` — the partial-current weight."""
    odn = _signed_projection(quad, axis, outward_sign)
    return np.asarray(quad.weights, dtype=float)[rows] * np.abs(odn[rows])


# ─────────────────────────────────────────────────────────────────────
# The domain IS Γ₊ (L0, parametrized over 3 quadratures)
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l0
class TestOutflowHalfTraceIsTheDomain:
    r"""Only the outgoing hemisphere contributes to the average.

    RE-POSED at **B3.4a**. Pre-B3.4a these fed the WHOLE face slot with
    ``+1`` on the outgoing ordinates and ``-1`` elsewhere, and asserted the
    result was ``+1`` — i.e. that the masked ``cos_w`` erased the incoming
    contribution. The narrowing makes that erasure unspellable: the incoming
    ordinates are not in the domain, so there is nothing to erase.

    What remains falsifiable, and is now the point of the class, is WHICH
    index set the domain is — see
    :meth:`test_domain_is_the_tangential_band_gamma_plus`, which pins it
    against the retired ``> 0.0`` classifier on the one fixture where the two
    disagree.
    """

    def test_lebedev_xmax_averages_the_outflow_values(self):
        quad = Quadrature.lebedev(17)
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        outflow, inflow = _half_traces(quad, "x", +1)
        # The pre-B3.4a probe, restricted: +1 on Γ₊, -1 everywhere else.
        # Feeding only γ₊ψ, the average must be +1 on every inflow row.
        psi_face = np.where(_signed_projection(quad, "x", +1) > 0, +1.0, -1.0)
        result = op.apply(psi_face[outflow])
        assert result.shape == (inflow.size,)
        np.testing.assert_allclose(result, +1.0, rtol=1e-14)

    def test_level_symmetric_ymin_averages_the_outflow_values(self):
        quad = Quadrature.level_symmetric(4)
        op = AngularAverageOperator.from_quadrature(quad, axis="y", outward_sign=-1)
        outflow, inflow = _half_traces(quad, "y", -1)
        # ymin face: outward normal −ŷ, so outgoing ⟺ mu_y < 0.
        psi_face = np.where(quad.mu_y < 0, +1.0, -1.0)
        result = op.apply(psi_face[outflow])
        assert result.shape == (inflow.size,)
        np.testing.assert_allclose(result, +1.0, rtol=1e-14)

    def test_gauss_legendre_1d_xmax_averages_the_outflow_values(self):
        quad = Quadrature.gauss_legendre(8)
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        outflow, inflow = _half_traces(quad, "x", +1)
        psi_face = np.where(quad.mu_x > 0, +1.0, -1.0)
        result = op.apply(psi_face[outflow])
        assert result.shape == (inflow.size,)
        np.testing.assert_allclose(result, +1.0, rtol=1e-14)

    def test_domain_is_the_tangential_band_gamma_plus(self):
        r"""The domain is the trace space's :math:`\Gamma_+` — NOT the retired
        private ``(outward_sign * mu) > 0.0`` set.

        This is the gate the B3.4a narrowing exists for. Pre-B3.4a
        :meth:`from_quadrature` carried its OWN outflow classifier, a strict
        ``> 0.0`` compare, which disagrees with the trace space's
        ``> TANGENTIAL_EPS`` wherever a quadrature carries tangential
        ordinates whose cosine is small-but-nonzero. `[M]` ``product(2, 4)``
        on ``xmax`` is exactly such a fixture: 4 of its 8 ordinates are
        tangential, the strict compare admits 4 rows and the tangential-band
        rule admits 2.

        The first assertion is the ACTIVATION guard — without it this test
        silently decays into a restatement of the sibling rows the moment the
        fixture stops discriminating (``vv`` Mode-8 class 7).
        """
        quad = Quadrature.product(n_mu=2, n_phi=4)
        odn = _signed_projection(quad, "x", +1)
        gamma_plus = np.flatnonzero(odn > +TANGENTIAL_EPS)
        retired_strict = np.flatnonzero(odn > 0.0)
        # ACTIVATION: the two classifiers must genuinely disagree here.
        assert not np.array_equal(gamma_plus, retired_strict), (
            f"fixture no longer discriminates: Γ₊ {gamma_plus.tolist()} and "
            f"the retired '> 0.0' set {retired_strict.tolist()} agree, so "
            f"this gate cannot see the classifier the narrowing retired."
        )
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        assert op.n_outflow == gamma_plus.size, (
            f"domain has {op.n_outflow} rows; Γ₊ (|Ω·n̂| > TANGENTIAL_EPS) "
            f"has {gamma_plus.size}. A domain of "
            f"{retired_strict.size} would be the retired '> 0.0' classifier."
        )
        # …and the strict-compare-sized input is REFUSED, so a caller that
        # still classifies the old way cannot feed this operator silently.
        with pytest.raises(ValueError, match=r"expected \|Γ₊\|"):
            op.apply(np.ones(retired_strict.size))


# ─────────────────────────────────────────────────────────────────────
# Conservation (L1)
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l1
class TestCosineWeightedCurrentConservation:
    r"""The current re-emitted on :math:`\Gamma_-` equals the one received on
    :math:`\Gamma_+`.

    RE-POSED and STRENGTHENED at **B3.4a**. The pre-B3.4a claim was
    ``Σ_{Γ₊} w|Ω·n| ψ == Σ_{Γ₊} w|Ω·n| (Aψ)`` — the OUTGOING weights applied
    to the broadcast OUTPUT. With ``A`` a full-face endomorphism that reading
    was available, but it is self-consistent by construction (both sides
    reduce to ``c · norm``), so its only teeth were against a wrong ``norm``.

    Narrowed, the physical statement is available directly and subsumes it:

    .. math::

        J_- \;=\; \sum_{\Omega \in \Gamma_-} w\,|\Omega\cdot\hat n|\,
                  (A\psi)(\Omega)
        \;=\; \sum_{\Omega \in \Gamma_+} w\,|\Omega\cdot\hat n|\,\psi(\Omega)
        \;=\; J_+ ,

    i.e. a white face at :math:`\alpha = 1` returns exactly the partial
    current it received (Bell & Glasstone 1970 §1.5). It still reds on a
    wrong ``norm`` (which rescales ``c``), and additionally on a wrong
    codomain size.

    It rides on the quadrature's face symmetry
    :math:`\sum_{\Gamma_+} w|\Omega\cdot\hat n| =
    \sum_{\Gamma_-} w|\Omega\cdot\hat n|`, asserted inline as an activation
    guard rather than assumed: `[M]` it holds BIT-IDENTICALLY on every
    quadrature × face pair in the tree, and a future asymmetric rule would
    make this operator genuinely non-conservative — which must red HERE,
    loudly, rather than be absorbed into a tolerance.
    """

    def _assert_face_symmetric(self, quad, axis, outward_sign):
        outflow, inflow = _half_traces(quad, axis, outward_sign)
        n_plus = _cosine_weights(quad, outflow, axis, outward_sign).sum()
        n_minus = _cosine_weights(quad, inflow, axis, outward_sign).sum()
        np.testing.assert_allclose(
            n_plus, n_minus, rtol=1e-15,
            err_msg=(
                "quadrature is not face-symmetric: the Γ₊ and Γ₋ "
                "cosine-weighted weight sums differ, so a Lambertian face "
                "cannot return the current it received."
            ),
        )
        return outflow, inflow

    def test_reemitted_current_equals_received_current_lebedev(self):
        quad = Quadrature.lebedev(17)
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        outflow, inflow = self._assert_face_symmetric(quad, "x", +1)
        rng = np.random.default_rng(42)
        psi_out = rng.uniform(0.5, 1.5, size=outflow.size)  # positive flux
        psi_in = op.apply(psi_out)
        j_plus = (_cosine_weights(quad, outflow, "x", +1) * psi_out).sum()
        j_minus = (_cosine_weights(quad, inflow, "x", +1) * psi_in).sum()
        np.testing.assert_allclose(j_minus, j_plus, rtol=1e-13)

    def test_reemitted_current_with_anisotropic_input(self):
        quad = Quadrature.level_symmetric(6)
        op = AngularAverageOperator.from_quadrature(quad, axis="z", outward_sign=+1)
        outflow, inflow = self._assert_face_symmetric(quad, "z", +1)
        # Strongly anisotropic input — exercises the cosine weighting.
        psi_face = np.exp(quad.mu_z) * (quad.mu_z + 2.0)
        psi_out = psi_face[outflow]
        psi_in = op.apply(psi_out)
        j_plus = (_cosine_weights(quad, outflow, "z", +1) * psi_out).sum()
        j_minus = (_cosine_weights(quad, inflow, "z", +1) * psi_in).sum()
        np.testing.assert_allclose(j_minus, j_plus, rtol=1e-13)


# ─────────────────────────────────────────────────────────────────────
# Reciprocity under the cosine-weighted inner product (L1)
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l1
class TestReciprocityCosineWeighted:
    r"""The Lambertian kernel is RECIPROCAL under the partial-current metric,
    and its weighted adjoint is the Lambertian kernel of the OPPOSITE face.

    RE-POSED at **B3.4a** from ``TestSelfAdjointnessCosineWeighted``. The old
    claim was ``⟨Ax, y⟩_w == ⟨x, Ay⟩_w`` — a self-adjointness statement that
    needs ``A`` to be an endomorphism. It held pre-B3.4a because the full-face
    operator broadcast over ALL ``N`` ordinates against the SAME masked
    ``cos_w`` it contracted with: as a matrix, ``A = 1 ⊗ (cos_w/norm)``, so
    ``W A = cos_w ⊗ (cos_w/norm) = (W A)ᵀ``. Narrowed, ``A : \Gamma_+ \to
    \Gamma_-`` is rectangular and that reading is unavailable.

    The mathematical content is reciprocity, which IS posable — and the
    adjoint is an object the tree already builds. With
    :math:`W_\pm = \operatorname{diag}(w\,|\Omega\cdot\hat n|)` on each
    half-trace, the weighted adjoint is

    .. math::

        A^{*} \;=\; W_+^{-1} A^{\mathsf T} W_-
        \;=\; \mathbf 1_{\Gamma_+} \otimes
              \frac{w\,|\Omega\cdot\hat n|}{\operatorname{norm}_+}
              \bigg|_{\Gamma_-} ,

    i.e. a Lambertian average that contracts over :math:`\Gamma_-` and
    re-emits on :math:`\Gamma_+` — which is precisely the Lambertian kernel
    of the REVERSED face, ``from_quadrature(quad, axis, -outward_sign)``,
    whenever :math:`\operatorname{norm}_+ = \operatorname{norm}_-` (the
    quadrature's face symmetry, asserted as an activation guard below and
    `[M]` bit-identical on every quadrature in the tree). Both are production
    objects; nothing is hand-built, so the gate pins the production kernel
    against the production kernel of the mirror face — which is the whole
    physical content of "white is self-adjoint in the partial-current inner
    product".

    Phase **B5** types this kernel as the rank-one :math:`u \otimes v` it is,
    at which point the adjoint becomes structurally available rather than
    coincidental; this gate is what B5 must keep green.
    """

    def _reciprocity(self, quad, axis, outward_sign, seed):
        outflow, inflow = _half_traces(quad, axis, outward_sign)
        w_plus = _cosine_weights(quad, outflow, axis, outward_sign)
        w_minus = _cosine_weights(quad, inflow, axis, outward_sign)
        # ACTIVATION: the reciprocity holds iff the two normalizations agree.
        np.testing.assert_allclose(
            w_plus.sum(), w_minus.sum(), rtol=1e-15,
            err_msg="quadrature not face-symmetric; reciprocity is not posed.",
        )
        A = AngularAverageOperator.from_quadrature(quad, axis, outward_sign)
        # The reversed face: its Γ₊ IS this face's Γ₋ and vice versa, so
        # A_star maps Γ₋ → Γ₊ exactly as the weighted adjoint must.
        A_star = AngularAverageOperator.from_quadrature(quad, axis, -outward_sign)
        assert (A_star.n_outflow, A_star.n_inflow) == (inflow.size, outflow.size), (
            "the reversed-face kernel does not map Γ₋ → Γ₊; the adjoint "
            "candidate is not even shape-compatible."
        )
        rng = np.random.default_rng(seed)
        psi = rng.uniform(0.0, 1.0, size=outflow.size)   # on Γ₊
        phi = rng.uniform(0.0, 1.0, size=inflow.size)    # on Γ₋
        lhs = (w_minus * A.apply(psi) * phi).sum()       # ⟨Aψ, φ⟩_{W₋}
        rhs = (w_plus * psi * A_star.apply(phi)).sum()   # ⟨ψ, A*φ⟩_{W₊}
        np.testing.assert_allclose(
            lhs, rhs, rtol=1e-13,
            err_msg=f"⟨Aψ,φ⟩_W₋ = {lhs!r} vs ⟨ψ,A*φ⟩_W₊ = {rhs!r}",
        )

    def test_reciprocity_lebedev_xmax(self):
        self._reciprocity(Quadrature.lebedev(17), "x", +1, seed=0)

    def test_reciprocity_level_symmetric_ymin(self):
        self._reciprocity(Quadrature.level_symmetric(4), "y", -1, seed=1)


# ─────────────────────────────────────────────────────────────────────
# Structural predicates (L0)
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l0
class TestPredicates:
    def test_apply_only(self):
        # Rank-deficient projection: neither invertible nor adjointable —
        # apply is the only verb (apply itself is universal, guarded at
        # composition time).
        quad = Quadrature.lebedev(5)
        op = AngularAverageOperator.from_quadrature(quad, axis="z", outward_sign=+1)
        assert not op.is_invertible and not op.is_adjointable


# ─────────────────────────────────────────────────────────────────────
# Axis selection — z on Lebedev OK, z on GL raises (L0)
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l0
class TestAxisDispatch:
    def test_z_on_lebedev_succeeds(self):
        r"""``axis="z"`` is realizable on a cubature with genuine ``mu_z``.

        RE-POSED at **B3.4a**: the pre-B3.4a assertion was
        ``op.n_ordinates == quad.N``, which named an attribute the narrowing
        retired AND asserted the full-face domain the narrowing removed. The
        surviving CLAIM — "z succeeds where the cubature carries a third
        cosine" — is re-stated on the two half-trace sizes, each pinned
        against the hand-spelled reference and each strictly below ``N``
        (a full-``N`` domain would be the retired shape).
        """
        quad = Quadrature.lebedev(11)
        op = AngularAverageOperator.from_quadrature(quad, axis="z", outward_sign=+1)
        outflow, inflow = _half_traces(quad, "z", +1)
        assert (op.n_outflow, op.n_inflow) == (outflow.size, inflow.size)
        assert op.n_outflow < quad.N and op.n_inflow < quad.N, (
            f"half-traces {op.n_outflow}/{op.n_inflow} are not strictly "
            f"inside the {quad.N}-ordinate face slot."
        )

    def test_z_on_gauss_legendre_raises(self):
        r"""1-D Gauss-Legendre carries ``mu_z == zeros(N)``, so a ``z`` face is
        a RANK MISMATCH between the face layout and the cubature.

        RE-POSED at **B3.4a**. Pre-B3.4a the operator classified its own
        hemisphere and reached its "no outgoing ordinates" guard, so the test
        pinned ``match="no outgoing"``. It now routes through
        :func:`~orpheus.numerics.spaces.angular_trace_space.build_omega_dot_n`
        — the single face-name → signed-projection primitive — which refuses
        the degenerate axis FIRST and by its real name. The message is
        strictly more specific, so the pattern tightens rather than widens.
        """
        quad = Quadrature.gauss_legendre(8)
        with pytest.raises(ValueError, match="requires genuine mu_z"):
            AngularAverageOperator.from_quadrature(quad, axis="z", outward_sign=+1)

    def test_invalid_axis_raises(self):
        quad = Quadrature.lebedev(5)
        with pytest.raises(ValueError, match="Unknown axis"):
            AngularAverageOperator.from_quadrature(quad, axis="w", outward_sign=+1)

    def test_invalid_outward_sign_raises(self):
        quad = Quadrature.lebedev(5)
        with pytest.raises(ValueError, match="outward_sign"):
            AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=0)


# ─────────────────────────────────────────────────────────────────────
# Construction guards (L0) — B3.4a: cos_w is the RESTRICTION to Γ₊
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l0
class TestConstructionGuards:
    r"""``cos_w`` must be the restriction to :math:`\Gamma_+`, not a masked
    full-face array.

    NET-NEW at **B3.4a**: the constructor's positivity test changed meaning
    (``>= 0`` → strictly ``> 0``) because its argument changed meaning, and
    the strengthened guard shipped with no negative test. A masked full-face
    ``cos_w`` — the pre-B3.4a shape — is now the defect the guard exists to
    name, so that is the input the negative feeds it.
    """

    def test_masked_full_face_cos_w_is_refused(self):
        quad = Quadrature.lebedev(17)
        outflow, inflow = _half_traces(quad, "x", +1)
        odn = _signed_projection(quad, "x", +1)
        masked_full = np.where(odn > 0, np.asarray(quad.weights) * odn, 0.0)
        with pytest.raises(ValueError, match="strictly positive"):
            AngularAverageOperator(
                cos_w=masked_full, norm=float(masked_full.sum()),
                n_inflow=int(inflow.size),
            )
        # Positive control: the RESTRICTION constructs cleanly.
        restricted = masked_full[outflow]
        AngularAverageOperator(
            cos_w=restricted, norm=float(restricted.sum()),
            n_inflow=int(inflow.size),
        )

    def test_non_positive_norm_is_refused(self):
        with pytest.raises(ValueError, match="norm must be positive"):
            AngularAverageOperator(cos_w=np.array([1.0, 2.0]), norm=0.0, n_inflow=2)

    def test_empty_codomain_is_refused(self):
        with pytest.raises(ValueError, match="n_inflow must be positive"):
            AngularAverageOperator(cos_w=np.array([1.0]), norm=1.0, n_inflow=0)


# ─────────────────────────────────────────────────────────────────────
# Defensive copy — operator is decoupled from the source quadrature (L0)
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l0
class TestDefensiveCopy:
    def test_quadrature_reference_not_held(self):
        quad = Quadrature.lebedev(5)
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        # No quadrature attribute on the operator (no held reference).
        assert not hasattr(op, "quadrature")
        assert not hasattr(op, "_quadrature")
        # The internal cos_w is its own copy.
        assert op._cos_w is not quad.weights
        assert op._cos_w is not quad.mu_x

    def test_output_is_fresh_array(self):
        """Calling code may mutate the output without affecting the
        operator's internal state."""
        quad = Quadrature.lebedev(5)
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        psi = np.ones(op.n_outflow)
        result1 = op.apply(psi)
        result1[0] = 999.0
        result2 = op.apply(psi)
        # result2 must NOT see the mutation (operator state unaffected).
        assert result2[0] != 999.0


# ─────────────────────────────────────────────────────────────────────
# Apply input-shape validation (L0)
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l0
class TestInputShape:
    def test_full_face_input_is_refused(self):
        r"""The whole face slot is NOT the domain.

        RE-POSED at **B3.4a** from a wrong-``N`` probe (13 rows against a
        14-ordinate quadrature). The narrowing's actual threat is a caller
        that has not migrated and still hands the operator the FULL face
        slot: the shapes are plausible, the values would be silently wrong,
        and the guard is the only thing between them.
        """
        quad = Quadrature.lebedev(5)  # N=14
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        assert op.n_outflow < quad.N  # activation: the two shapes differ
        with pytest.raises(ValueError, match=r"expected \|Γ₊\|"):
            op.apply(np.ones(quad.N))
        with pytest.raises(ValueError, match="psi.shape"):
            op.apply(np.ones(op.n_outflow + 1))

    def test_multi_dim_input_broadcasts_over_the_codomain(self):
        r"""Shape ``(|Γ₊|, 5, 3)`` input — the average is taken over the
        leading axis only and broadcast over the ``|Γ₋|`` inflow rows.

        RE-POSED at **B3.4a** from ``result.shape == psi.shape``. `[M]`
        ``|Γ₊| == |Γ₋|`` on EVERY quadrature × face pair in the tree, so that
        assertion is inside the shape functional's invariance group (``vv``
        Mode 12) — it would stay green under a codomain that echoed the
        domain. The row count is therefore asserted against ``n_inflow``, and
        :meth:`test_codomain_size_follows_n_inflow_not_the_input` breaks the
        degeneracy directly.
        """
        quad = Quadrature.lebedev(11)
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        rng = np.random.default_rng(7)
        psi = rng.uniform(0, 1, size=(op.n_outflow, 5, 3))
        result = op.apply(psi)
        assert result.shape == (op.n_inflow, 5, 3)
        # All leading-axis slices identical (Lambertian: angle-independent).
        for n in range(op.n_inflow):
            np.testing.assert_allclose(result[n], result[0], rtol=1e-14)

    def test_codomain_size_follows_n_inflow_not_the_input(self):
        r"""The number of emitted rows comes from the CODOMAIN datum, not
        from the domain.

        NET-NEW at **B3.4a**, and the anti-Mode-12 leg of the shape contract:
        every reachable fixture has ``|Γ₊| == |Γ₋|``, so no
        ``from_quadrature``-built operator can distinguish "emits ``|Γ₋|``
        rows" from "echoes its input's leading axis". Constructing with
        deliberately unequal half-trace sizes is the only way to see it.
        """
        op = AngularAverageOperator(
            cos_w=np.array([1.0, 2.0, 3.0]), norm=6.0, n_inflow=5,
        )
        assert (op.n_outflow, op.n_inflow) == (3, 5)
        out = op.apply(np.array([1.0, 1.0, 1.0]))
        assert out.shape == (5,), (
            f"emitted {out.shape[0]} rows for a 3-row input with |Γ₋| = 5 — "
            f"the codomain is echoing the domain."
        )
        np.testing.assert_allclose(out, 1.0, rtol=1e-14)


# ─────────────────────────────────────────────────────────────────────
# Legacy-vs-realiser bit-equivalence test removed in Issue #186 (B3 + β2)
# ─────────────────────────────────────────────────────────────────────
#
# The pre-#186 TestLegacyBitEquivalence class compared
# ``AngularAverageOperator.apply(psi)`` to ``WhiteBoundary.apply(psi, quad)``.
# That legacy method no longer exists; the comparison would now be
# circular (the realiser-path test pins the same equivalence in
# tests/sn/operators/test_sn_boundary_realizer.py). The Wave-1
# current-conservation + reciprocity tests above are the
# structurally-independent references for AngularAverageOperator
# correctness.
