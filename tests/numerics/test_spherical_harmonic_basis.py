r"""Tests for :class:`orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis`.

The class's :meth:`evaluate` / :meth:`evaluate_from_components` methods
return the :math:`(N, L+1, 2L+1)` table of real spherical harmonics
under the no-:math:`4\pi/(2\ell+1)`-prefactor convention used by the
SN Galerkin reconstruction. The invariants below are the test hooks
the architecture catalog (§B-table) names for the harmonic basis:

* ``Y_0^0 = 1`` (P0 normalisation).
* Hard-coded P1 values: :math:`Y_1^{-1} = \mu_z`, :math:`Y_1^0 = \mu_x`,
  :math:`Y_1^{+1} = \mu_y`.
* Addition theorem :math:`\sum_m Y_\ell^m(\Omega) Y_\ell^m(\Omega') =
  P_\ell(\Omega \cdot \Omega')` — the identity consumed by the Pℓ
  source build.
* Discrete orthogonality on a Lebedev grid:
  :math:`\sum_n w_n Y_\ell^m(\Omega_n) Y_{\ell'}^{m'}(\Omega_n) =
  (4\pi / (2\ell+1)) \delta_{\ell\ell'} \delta_{mm'}` for
  :math:`\ell + \ell' \le \deg`.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.special import eval_legendre

from orpheus.numerics.basis.spherical_harmonic_basis import (
    MirrorEvenSphericalHarmonicBasis,
    SphericalHarmonicBasis,
)
from orpheus.numerics.manifold import SPHERE
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.quadrature import Quadrature, lebedev_sphere
from orpheus.numerics.symmetry import SubgroupOfO3


def _unit(directions: np.ndarray) -> np.ndarray:
    """Row-normalise to genuine points of :math:`S^2` (0.6 refuses anything else)."""
    directions = np.asarray(directions, dtype=float)
    return directions / np.linalg.norm(directions, axis=1)[:, None]


@pytest.mark.l0
class TestEvaluateRealSHBasic:
    """L0: hard-coded P0/P1 invariants."""

    def test_l_negative_rejected_at_construction(self):
        """Negative L is semantically meaningless for spherical harmonics;
        :class:`SphericalHarmonicBasis` rejects it at ``__post_init__``.

        Tightens the contract relative to the pre-P3.2 ``evaluate_real_sh``
        free function, which silently returned an empty array on
        ``L < 0`` — a permissive behaviour with no consumer."""
        with pytest.raises(ValueError, match="L must be non-negative"):
            SphericalHarmonicBasis(L=-1)

    def test_l_zero_p0_unity(self):
        r""":math:`Y_0^0 \equiv 1` at three genuine points of :math:`S^2`.

        ⚠ The directions are NORMALISED. Until 2026-09-02 this row typed
        ``(0.5, 0, 0)`` and friends — points of the ball, not of the sphere —
        and got away with it because ``evaluate`` did not check. #429's
        tracker 0.6 wired the membership refusal in, and the claim ("the
        :math:`\ell = 0` slot is unity") is carried perfectly well by unit
        directions: the fixture was never load-bearing, the sloppiness was.
        """
        directions = _unit(np.array([[0.5, 0.3, -0.2], [-0.5, 0.1, 0.7], [1.0, 0.4, 0.4]]))
        Y = SphericalHarmonicBasis(L=0).evaluate(directions)
        assert Y.shape == (3, 1, 1)
        np.testing.assert_array_equal(Y[:, 0, 0], 1.0)

    def test_l_one_p1_hardcoded(self):
        r"""The :math:`\ell = 1` slots are the Cartesian cosines, in the chart's own order.

        ⭐ This row is why ERR-080 read as an ":math:`L \ge 2` problem": the
        :math:`\ell = 1` branch is hard-coded in Cartesian form and never
        computes an azimuth, so a forged :math:`(\mu, 0, 0)` came out clean
        here and dirty two degrees up.
        """
        directions = _unit(
            np.array([[0.6, 0.0, 0.5], [-0.3, 0.4, 0.7], [0.8, -0.5, 0.1]])
        )
        mu_x, mu_y, mu_z = directions[:, 0], directions[:, 1], directions[:, 2]
        Y = SphericalHarmonicBasis(L=1).evaluate(directions)
        assert Y.shape == (3, 2, 3)
        np.testing.assert_array_equal(Y[:, 0, 0], 1.0)
        np.testing.assert_array_equal(Y[:, 1, 0], mu_z)   # m = -1
        np.testing.assert_array_equal(Y[:, 1, 1], mu_x)   # m =  0
        np.testing.assert_array_equal(Y[:, 1, 2], mu_y)   # m = +1

        # and the un-normalised fixture this row used to carry is now the
        # REFUSAL leg — the same numbers, off the sphere (0.6).
        with pytest.raises(ValueError, match="not points of S\\^2"):
            SphericalHarmonicBasis(L=1).evaluate(
                np.array([[0.6, 0.0, 0.5], [-0.3, 0.4, 0.7], [0.8, -0.5, 0.1]])
            )


@pytest.mark.l1
class TestAdditionTheorem:
    r"""L1: :math:`\sum_m Y_\ell^m(\Omega) Y_\ell^m(\Omega') =
    P_\ell(\Omega \cdot \Omega')`.

    This is the structural identity the no-:math:`4\pi/(2\ell+1)`
    convention is *defined by*. Verifying it on a few random ordinate
    pairs gates the entire scaling factor.
    """

    @pytest.mark.parametrize("L", [2, 3, 4, 5])
    def test_addition_theorem_pairs(self, L):
        rng = np.random.default_rng(seed=2026 + L)
        basis = SphericalHarmonicBasis(L=L)
        # Sample two random unit vectors per trial, 8 trials.
        for trial in range(8):
            v = rng.standard_normal((2, 3))
            v /= np.linalg.norm(v, axis=1, keepdims=True)
            mu_x = v[:, 0]
            mu_y = v[:, 1]
            mu_z = v[:, 2]
            Y = basis.evaluate_from_components(mu_x, mu_y, mu_z)
            # cos(angle) = Ω · Ω'
            cos_alpha = float(np.dot(v[0], v[1]))
            for l in range(L + 1):
                lhs = float(np.sum(Y[0, l, : 2 * l + 1] * Y[1, l, : 2 * l + 1]))
                rhs = float(eval_legendre(l, cos_alpha))
                assert abs(lhs - rhs) < 1e-12, (
                    f"addition theorem mismatch at L={L}, l={l}, "
                    f"trial={trial}: lhs={lhs}, rhs={rhs}"
                )


@pytest.mark.l1
class TestDiscreteOrthogonalityOnLebedev:
    r"""L1: discrete orthogonality on a Lebedev quadrature.

    For the no-:math:`4\pi/(2\ell+1)`-prefactor convention,

    .. math::

        \sum_n w_n Y_\ell^m(\Omega_n) Y_{\ell'}^{m'}(\Omega_n)
        = \frac{4\pi}{2\ell+1} \delta_{\ell\ell'} \delta_{mm'}

    holds for :math:`\ell + \ell' \le \mathrm{deg}` where ``deg`` is
    the spherical-harmonic exactness degree of the Lebedev rule.
    """

    @pytest.mark.parametrize("order", [3, 5, 7])
    def test_lebedev_orthogonality(self, order):
        measure = lebedev_sphere(order)
        nodes = measure.nodes  # (N, 3)
        mu_x, mu_y, mu_z = nodes[:, 0], nodes[:, 1], nodes[:, 2]
        w = measure.weights
        deg = measure.degree_of_exactness
        L = deg // 2  # safe order so ℓ+ℓ' ≤ deg
        Y = SphericalHarmonicBasis(L=L).evaluate_from_components(mu_x, mu_y, mu_z)
        # Build flat (Y_lm) basis indexed by (l, m_offset).
        flat = []
        labels = []
        for l in range(L + 1):
            for m in range(-l, l + 1):
                flat.append(Y[:, l, l + m])
                labels.append((l, m))
        flat = np.array(flat)  # (n_basis, N)
        gram = np.einsum("an,bn,n->ab", flat, flat, w)
        for a, (la, ma) in enumerate(labels):
            for b, (lb, mb) in enumerate(labels):
                expected = (4.0 * np.pi / (2 * la + 1)) if (la == lb and ma == mb) else 0.0
                assert abs(gram[a, b] - expected) < 1e-10, (
                    f"orthogonality fails at "
                    f"(l,m)=({la},{ma}) vs ({lb},{mb}): "
                    f"got {gram[a, b]}, expected {expected}"
                )


@pytest.mark.l0
def test_on_axis_safe_no_division_by_zero():
    """When sin(theta) ≈ 0 (mu_x = ±1), the cos_phi/sin_phi expressions
    must NOT divide by zero. The implementation uses ``np.where`` to
    sentinel; verify L=3 evaluates without warnings on a pure-x input."""
    mu_x = np.array([1.0, -1.0, 0.0])
    mu_y = np.array([0.0, 0.0, 1.0])
    mu_z = np.array([0.0, 0.0, 0.0])
    with np.errstate(divide="raise", invalid="raise"):
        Y = SphericalHarmonicBasis(L=3).evaluate_from_components(mu_x, mu_y, mu_z)
    assert np.all(np.isfinite(Y))


# ══════════════════════════════════════════════════════════════════════
# 0.6 (#429) — a real spherical harmonic eats a POINT of S^2
# ══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.catches("ERR-080")
class TestEvaluateRefusesDirectionsOffTheSphere:
    r"""B11 — the membership refusal, in kind, at the basis's own door.

    ERR-080's construction padded a 1-D rule's :math:`\mu` to
    :math:`(\mu, 0, 0)` and declared the result a measure on :math:`S^2`.
    Those points have :math:`|\Omega| = |\mu| < 1` off the poles, so they are
    not points of this basis's domain, and ``arctan2(0, 0) = 0`` read the
    fabricated azimuth as a real one. The repair refuses them here rather
    than tabulating them wrong.

    ⚠ **The witness has to be BUILT** — production no longer constructs the
    forged measure (``_harmonic_frame_measure`` was deleted in the same
    commit), so the rows below assemble it from the raw
    :class:`~orpheus.numerics.measure.DiscreteMeasure` constructor. That is
    ``plan-authoring`` §6c in its exact shape: a gate whose rejected input
    the same commit removes owes a constructed witness, or it lands green and
    unable to fail.

    ⛔ **The refusal is a real** ``raise``, **not a bare** ``assert``. The
    canonical invocation is ``python -O``, which strips ``assert`` in
    ``orpheus/`` at compile time (``vv-principles`` Mode 8), so a contract
    spelled that way would be a no-op in the suite that matters.
    :func:`test_the_refusal_is_a_raise_and_survives_dash_O` verifies that
    directly, by running the guard's own arithmetic under both interpreter
    modes in a subprocess.
    """

    def test_the_forged_one_d_ordinates_are_refused(self) -> None:
        r"""``[M]`` GL8's forged nodes have norms ``0.183435 … 0.960290`` — every one off :math:`S^2` bar the poles."""
        mu = np.asarray(Quadrature.gauss_legendre(8).measure.nodes, dtype=float).reshape(-1)
        zeros = np.zeros_like(mu)
        forged = np.column_stack([mu, zeros, zeros])

        norms = np.linalg.norm(forged, axis=1)
        assert norms.min() < 0.2 and norms.max() < 0.97, f"norms {norms}"
        assert not SPHERE.contains(forged)

        with pytest.raises(ValueError, match=r"8 of 8 directions are not points of S\^2"):
            SphericalHarmonicBasis(L=2).evaluate(forged)

        # the same refusal one layer in: evaluate_from_components routes
        # through evaluate, so the component door cannot bypass it.
        with pytest.raises(ValueError, match="not points of S"):
            SphericalHarmonicBasis(L=2).evaluate_from_components(mu, zeros, zeros)

        # …and the forged MEASURE, assembled the way production used to —
        # this is the object the fix retired, kept alive here as the witness.
        forged_measure = DiscreteMeasure(
            nodes=forged,
            weights=np.asarray(Quadrature.gauss_legendre(8).weights, dtype=float),
            support=SPHERE,
        )
        with pytest.raises(ValueError, match="not points of S"):
            SphericalHarmonicBasis(L=2).mass_matrix(forged_measure)

    def test_honest_three_d_nodes_are_admitted(self) -> None:
        r"""The positive leg — every shipped sphere rule tabulates, at every :math:`L` used here.

        Without it "the basis refuses" is compatible with "the basis refuses
        everything" (``vv-principles`` #11).
        """
        for build in (
            lambda: Quadrature.level_symmetric(4),
            lambda: Quadrature.lebedev(11),
            lambda: Quadrature.product(n_mu=4, n_phi=6),
            lambda: Quadrature.folded_product(n_mu=4, n_phi=8),
        ):
            nodes = np.asarray(build().measure.nodes, dtype=float)
            assert SPHERE.contains(nodes)
            table = SphericalHarmonicBasis(L=3).evaluate(nodes)
            assert table.shape == (nodes.shape[0], 4, 7)
            assert np.all(np.isfinite(table))

    def test_a_single_off_sphere_direction_among_many_is_still_refused(self) -> None:
        r"""The count in the message is the DIAGNOSTIC — one bad row out of many must not be averaged away."""
        nodes = np.asarray(Quadrature.level_symmetric(4).measure.nodes, dtype=float).copy()
        nodes[3] *= 0.5
        with pytest.raises(ValueError, match=r"1 of 24 directions are not points"):
            SphericalHarmonicBasis(L=2).evaluate(nodes)

    def test_the_refusal_is_a_raise_and_survives_dash_O(self) -> None:
        r"""⛔ Mode 8 — the guard must fire under ``python -O``, which is the canonical runner.

        A bare ``assert`` in ``orpheus/`` is stripped at compile time by
        ``-O``, so a contract written that way ships accepting exactly the
        input it was written to refuse. This row runs the real guard in a
        FRESH ``-O`` interpreter (a subprocess: this test module is
        AST-rewritten by pytest, so an in-process check would measure
        pytest's rewriter, not the interpreter) and requires the
        ``ValueError``. A companion arm proves the same call in a plain
        interpreter behaves identically — if the two modes ever disagree,
        the guard has decayed into an ``assert``.
        """
        import subprocess
        import sys

        program = (
            "import numpy as np\n"
            "from orpheus.numerics.basis.spherical_harmonic_basis import "
            "SphericalHarmonicBasis\n"
            "forged = np.array([[0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]])\n"
            "try:\n"
            "    SphericalHarmonicBasis(L=2).evaluate(forged)\n"
            "except ValueError as exc:\n"
            "    print('RAISED' if 'not points of S^2' in str(exc) else 'WRONG')\n"
            "else:\n"
            "    print('ACCEPTED')\n"
        )
        for flags in ([sys.executable, "-O", "-c"], [sys.executable, "-c"]):
            result = subprocess.run(
                [*flags, program], capture_output=True, text=True, check=True
            )
            assert result.stdout.strip() == "RAISED", (
                f"{' '.join(flags[:-1])}: the off-sphere guard did not fire "
                f"({result.stdout.strip()!r}) — it is an inert `assert`, not a "
                f"contract (vv-principles Mode 8)"
            )


# ══════════════════════════════════════════════════════════════════════
# B6 — the σ-even mask is READ off the fold's catalogue entry
# ══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestEvenSlotMaskThroughTheEntry:
    r"""``even_slot_mask`` is the entry's descending-slot probe, not a private copy.

    Before #429 this class carried its own five-direction parity probe. The
    predicate — *"constant on the fibres of the quotient map"* — is the same
    one :class:`~orpheus.numerics.basis.descent.Descent` needs for every
    entry, so spelling it twice was a Cardinal-Rule-2 twin. The mask now
    READS :meth:`Quotient.descending_slots`, and these rows pin that the read
    is faithful and that it still depends on the group.
    """

    def test_the_mask_is_bit_identical_through_the_entrys_probe(self) -> None:
        r"""``[M]`` 2026-09-02: ``array_equal`` on 15 of 15 (3 axes × :math:`L \in \{1..5\}`) rows."""
        for axis_index, axis_name in enumerate("xyz"):
            for L in (1, 2, 3, 4, 5):
                basis = MirrorEvenSphericalHarmonicBasis(L=L, mirror_axis=axis_index)
                through_entry = SPHERE.quotient(
                    SubgroupOfO3.Mirror(axis_name)
                ).descending_slots(SphericalHarmonicBasis(L=L)).astype(float)
                assert np.array_equal(basis.even_slot_mask, through_entry), (
                    f"axis {axis_name}, L={L}"
                )

    def test_the_mask_moves_when_the_mirror_axis_moves(self) -> None:
        r"""NEGATIVE leg — the mask is a function of the GROUP, not a constant wearing one.

        ``[M]`` the mask WEIGHTS agree across axes (5 / 12 / 22 / 35 / 51 at
        :math:`L = 1..5`), so a cardinality check would be blind here; the
        PATTERNS differ, which is what this row asserts.
        """
        for L in (1, 2, 3, 4, 5):
            masks = [
                MirrorEvenSphericalHarmonicBasis(L=L, mirror_axis=a).even_slot_mask
                for a in range(3)
            ]
            assert len({int(m.sum()) for m in masks}) == 1, (
                "the cardinalities agree — a count-only gate would be blind"
            )
            for i in range(3):
                for j in range(i + 1, 3):
                    assert not np.array_equal(masks[i], masks[j]), (
                        f"L={L}: axes {i} and {j} give the SAME mask"
                    )

    def test_the_masked_basis_zeroes_exactly_the_sigma_odd_slots(self) -> None:
        r"""The mask's CONSUMER: the σ-odd columns come out as exact ``0.0``, the even ones untouched."""
        nodes = np.asarray(Quadrature.folded_product(n_mu=4, n_phi=8).measure.nodes, dtype=float)
        L = 3
        parent = SphericalHarmonicBasis(L=L).evaluate(nodes)
        masked = MirrorEvenSphericalHarmonicBasis(L=L, mirror_axis=1).evaluate(nodes)
        mask = MirrorEvenSphericalHarmonicBasis(L=L, mirror_axis=1).even_slot_mask

        assert np.array_equal(masked[:, mask == 0.0], np.zeros_like(masked[:, mask == 0.0]))
        assert np.array_equal(masked[:, mask == 1.0], parent[:, mask == 1.0])
