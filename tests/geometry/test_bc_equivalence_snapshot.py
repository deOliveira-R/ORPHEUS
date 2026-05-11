r"""BC equivalence snapshot harness (Wave 6 / C6.1).

This file is the safety net for Waves 7-11 of the boundary-operator
refactor (see ``.claude/plans/transient-giggling-cake.md``). For each of
the 8 representative BC × quadrature × input cases captured by
:mod:`tests.geometry._generate_bc_equivalence_snapshots`, two
assertions are pinned on every pytest invocation:

1. **Legacy 2-arg path.** ``bc.apply(psi_out, quadrature)`` reproduces
   the snapshot at the per-case tolerance (bit-exact for pure
   permutation / multiplication / identity bodies; ``nulp=4`` for the
   cosine-weighted averaging body; ``nulp=64`` for the
   sum-of-products mixed-BC reduction tree).
2. **Realized 1-arg path.** ``SNBoundaryRealizer().realize(bc,
   SNMethodSpace.minimal(quadrature)).apply(psi_out)`` reproduces the
   snapshot at the same per-case tolerance — proving the Wave-5
   functional realizer did not drift relative to the legacy code.

The exception is the vacuum case, whose two assertions are
SEMANTICALLY DIFFERENT (the legacy zeros every ordinate; the realizer
zeros only the inflow ordinates per :ref:`vv-principles` and Grand
Report §16A.5). The vacuum class therefore pins two distinct
expectations: the legacy ``zeros_like`` output AND the §16A.5
inflow-only mask. Wave 8 will update the *legacy* test to follow the
§16A.5 semantics (the legacy code itself disappears in Wave 11); the
*realized* test already follows §16A.5 and stays unchanged.

V&V tags
--------

* ``@pytest.mark.l1`` — cross-implementation L1 (the realizer path is a
  structurally-independent re-implementation of the legacy two-arg
  ``apply``; the snapshot pins the legacy reference exactly).
* ``@pytest.mark.regression`` — gates the Wave 7-11 production
  call-site migration. Failure means the refactor drifted numerically;
  see the Wave 6 closeout commit message for the bit-identity vs
  principled-equivalence three-criteria rubric (``vv-principles``
  skill).
* No ``@pytest.mark.verifies`` decorator — these tests verify a
  *software invariant* (refactor stability), not an equation. The
  snapshot is the contract, not a Sphinx label.

Generator
---------

The snapshots live under ``tests/geometry/snapshots/`` and are
committed to the repository. To regenerate (e.g. when a documented
semantic transition lands)::

    python -m tests.geometry._generate_bc_equivalence_snapshots

Per ``.claude/plans/transient-giggling-cake.md`` (Wave 6 risk
register), regeneration must be accompanied by an inline justification
citing one of the three principled-equivalence criteria.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from orpheus.geometry.boundary import (
    AlbedoBoundaryOperator,
    PeriodicBoundaryOperator,
    SpecularBoundaryOperator,
    VacuumBoundaryOperator,
    WhiteBoundaryOperator,
)
from orpheus.sn.boundary_realizer import SNBoundaryRealizer, SNMethodSpace
from orpheus.sn.quadrature import (
    GaussLegendre1D,
    LebedevSphere,
    LevelSymmetricSN,
)


pytestmark = [pytest.mark.l1, pytest.mark.regression]


_SNAPSHOT_DIR = Path(__file__).parent / "snapshots"
_GENERATOR_HINT = (
    "Run `python -m tests.geometry._generate_bc_equivalence_snapshots` "
    "to create."
)


def _load_or_skip(case_id: str) -> np.lib.npyio.NpzFile:
    """Load the snapshot for ``case_id`` or skip the test.

    Tests are skipped (NOT failed) if the snapshot is missing — this
    lets the test infrastructure land before the snapshots themselves
    in a decoupled-commit roll-out, and protects against running a
    stale-snapshot CI gate.
    """
    path = _SNAPSHOT_DIR / f"bc_equivalence_{case_id}.npz"
    if not path.exists():
        pytest.skip(
            f"Snapshot {path.name} not generated. {_GENERATOR_HINT}",
        )
    return np.load(path)


# ─────────────────────────────────────────────────────────────────────
# Case 1 — Vacuum × Lebedev(17)
# ─────────────────────────────────────────────────────────────────────


class TestVacuumLebedev17Snapshot:
    """Vacuum BC on Lebedev-17 quadrature.

    SEMANTIC NOTE: the legacy ``VacuumBoundaryOperator.apply`` returns
    ``np.zeros_like(psi_out)`` (zeroes EVERY ordinate). The Wave-5/8
    realizer returns an :class:`IncomingOrdinateMaskTensor` that zeros
    ONLY the inflow ordinates (outflow rows pass through unchanged) —
    the intentional §16A.5 semantic correction. This class pins BOTH
    expectations: legacy bit-exact zeros, realizer inflow-only-zeroing.

    Wave 8 will update the legacy-path assertion when the legacy code
    is deleted; this is principled per the three criteria in the
    ``vv-principles`` skill (the new behaviour is the structurally-
    independent §16A.5 reference, the FP drift is zero because the
    operator is a mask, and each intermediate IS a named, inspectable
    quantity — the inflow-index set).
    """

    case_id = "vacuum_lebedev17"

    @pytest.fixture(scope="class")
    def snapshot(self) -> np.lib.npyio.NpzFile:
        return _load_or_skip(self.case_id)

    def test_legacy_apply_returns_zeros(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Legacy ``bc.apply(psi, quad)`` returns the snapshotted zeros.

        Wave-6 ship-state assertion: bit-exact zeros everywhere. The
        snapshot's ``psi_in`` is ``np.zeros_like(psi_out)`` by
        construction.

        TODO(Wave 8): when ``SNMesh._resolve_bcs`` switches to the
        :class:`SNBoundaryRealizer`, the production legacy path
        disappears and this assertion either (a) is dropped along
        with the legacy code in Wave 11, or (b) is retained as a
        compatibility-shim regression while the shim lives. The
        semantic-transition criterion is documented in the
        ``vv-principles`` bit-identity vs principled-equivalence
        section.
        """
        quad = LebedevSphere.create(17)
        bc = VacuumBoundaryOperator()
        psi_out = snapshot["psi_out"]
        actual = bc.apply(psi_out, quad)
        # Bit-exact: legacy body is ``np.zeros_like``, single operation.
        np.testing.assert_array_equal(
            actual,
            snapshot["psi_in"],
            err_msg=(
                f"Legacy VacuumBoundaryOperator.apply drifted for "
                f"case {self.case_id!r}"
            ),
        )

    def test_realizer_zeroes_only_inflow_per_section_16A5(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Realizer (§16A.5 semantics) zeroes only the inflow ordinates.

        This is the FORWARD-LOOKING assertion: the realized
        :class:`IncomingOrdinateMaskTensor` masks only inflow rows of
        the xmin face. Outflow rows pass through unchanged. This is
        the §16A.5 intentional correction; the snapshot records the
        inflow index set so this assertion has a deterministic target
        without re-deriving the index set at test time.
        """
        quad = LebedevSphere.create(17)
        inflow_indices = snapshot["inflow_indices_xmin"]
        space = SNMethodSpace(
            quadrature=quad, face="xmin", inflow_indices=inflow_indices,
        )
        op = SNBoundaryRealizer().realize(VacuumBoundaryOperator(), space)
        psi_out = snapshot["psi_out"]
        actual = op.apply(psi_out)
        # Inflow rows are zero (§16A.5 mask action).
        np.testing.assert_array_equal(
            actual[inflow_indices], 0.0,
            err_msg=(
                "§16A.5 realizer should zero ALL inflow ordinates"
            ),
        )
        # Outflow / tangential rows pass through unchanged — the §16A.5
        # correction. Snapshot's ``psi_out`` is the source of truth for
        # the original input.
        non_inflow = np.setdiff1d(
            np.arange(quad.N), inflow_indices,
        )
        np.testing.assert_array_equal(
            actual[non_inflow], psi_out[non_inflow],
            err_msg=(
                "§16A.5 realizer should pass OUTFLOW/TANGENTIAL "
                "ordinates through unchanged"
            ),
        )


# ─────────────────────────────────────────────────────────────────────
# Case 2 — Albedo(0.5) × Lebedev(17)
# ─────────────────────────────────────────────────────────────────────


class TestAlbedo05Lebedev17Snapshot:
    """Pure albedo: ``psi_in = 0.5 * psi_out``.

    Bit-exact for both paths: the legacy body is a single multiplication;
    the realized op is :class:`ScaledOperator(0.5, IdentityOperator())`,
    same single multiplication.
    """

    case_id = "albedo_05_lebedev17"

    @pytest.fixture(scope="class")
    def snapshot(self) -> np.lib.npyio.NpzFile:
        return _load_or_skip(self.case_id)

    def test_legacy_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """``AlbedoBoundaryOperator(0.5).apply`` is bit-exact."""
        quad = LebedevSphere.create(17)
        bc = AlbedoBoundaryOperator(albedo=0.5)
        actual = bc.apply(snapshot["psi_out"], quad)
        np.testing.assert_array_equal(actual, snapshot["psi_in"])

    def test_realizer_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Realized ``ScaledOperator(0.5, IdentityOperator)`` is bit-exact.

        Multiplication is the only operation; no reduction tree change
        between legacy ``alpha * psi_out`` and ``ScaledOperator(alpha,
        IdentityOperator()).apply(psi_out)``. Bit-exact is the right
        gate.
        """
        quad = LebedevSphere.create(17)
        bc = AlbedoBoundaryOperator(albedo=0.5)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        actual = op.apply(snapshot["psi_out"])
        np.testing.assert_array_equal(actual, snapshot["psi_in"])


# ─────────────────────────────────────────────────────────────────────
# Case 3 — Specular(x, α=1) × Lebedev(17)
# ─────────────────────────────────────────────────────────────────────


class TestSpecularXLebedev17Snapshot:
    """Specular x-axis at α=1: pure index permutation, bit-exact."""

    case_id = "specular_x_lebedev17"

    @pytest.fixture(scope="class")
    def snapshot(self) -> np.lib.npyio.NpzFile:
        return _load_or_skip(self.case_id)

    def test_legacy_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Legacy ``psi_out[ref]`` permutation is bit-exact."""
        quad = LebedevSphere.create(17)
        bc = SpecularBoundaryOperator(axis="x", albedo=1.0)
        actual = bc.apply(snapshot["psi_out"], quad)
        np.testing.assert_array_equal(actual, snapshot["psi_in"])

    def test_realizer_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Realized bare :class:`PermutationOperator` is bit-exact.

        At α=1 the realizer's fast path returns the bare
        :class:`PermutationOperator` (no :class:`ScaledOperator` wrap),
        preserving bit-identity with legacy.
        """
        quad = LebedevSphere.create(17)
        bc = SpecularBoundaryOperator(axis="x", albedo=1.0)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        actual = op.apply(snapshot["psi_out"])
        np.testing.assert_array_equal(actual, snapshot["psi_in"])


# ─────────────────────────────────────────────────────────────────────
# Case 4 — Specular(y, α=0.7) × LevelSymmetric(6)
# ─────────────────────────────────────────────────────────────────────


class TestSpecularYPartial07LS6Snapshot:
    """Specular y-axis with partial albedo α=0.7.

    Multiplication-then-permutation is bit-exact: legacy does ``albedo *
    psi_out[ref]``; the realized :class:`ScaledOperator(0.7,
    PermutationOperator)` does the same (the permutation accesses
    elements, then ``ScaledOperator.apply`` multiplies). The reduction
    tree is identical so ``assert_array_equal`` is the right gate.
    """

    case_id = "specular_y_partial_07_LS6"

    @pytest.fixture(scope="class")
    def snapshot(self) -> np.lib.npyio.NpzFile:
        return _load_or_skip(self.case_id)

    def test_legacy_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Legacy ``0.7 * psi_out[ref]`` is bit-exact."""
        quad = LevelSymmetricSN.create(sn_order=6)
        bc = SpecularBoundaryOperator(axis="y", albedo=0.7)
        actual = bc.apply(snapshot["psi_out"], quad)
        np.testing.assert_array_equal(actual, snapshot["psi_in"])

    def test_realizer_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Realized ``ScaledOperator(0.7, PermutationOperator)`` is bit-exact."""
        quad = LevelSymmetricSN.create(sn_order=6)
        bc = SpecularBoundaryOperator(axis="y", albedo=0.7)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        actual = op.apply(snapshot["psi_out"])
        np.testing.assert_array_equal(actual, snapshot["psi_in"])


# ─────────────────────────────────────────────────────────────────────
# Case 5 — White(xmax, α=1) × LevelSymmetric(4)
# ─────────────────────────────────────────────────────────────────────


class TestWhiteXmaxLS4Snapshot:
    """White (Lambertian) xmax at α=1.

    The cosine-weighted hemisphere average is an N-term reduction. Wave
    1's ``TestLegacyBitEquivalence`` already pinned
    :class:`AngularAverageOperator.from_quadrature` body bit-exact at
    α=1, so the realizer's fast path (bare
    :class:`AngularAverageOperator`) is bit-equivalent today. ``nulp=4``
    is the safety margin against post-Wave-7 compositional reordering
    (per the plan's table).
    """

    case_id = "white_xmax_LS4"

    @pytest.fixture(scope="class")
    def snapshot(self) -> np.lib.npyio.NpzFile:
        return _load_or_skip(self.case_id)

    def test_legacy_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Legacy ``bc.apply`` matches the snapshot bit-for-bit.

        The snapshot WAS generated FROM this code path so the legacy
        side is necessarily bit-exact — this is the pinning gate that
        future Waves' refactors of the legacy ``WhiteBoundaryOperator``
        body must not violate.
        """
        quad = LevelSymmetricSN.create(sn_order=4)
        bc = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0)
        actual = bc.apply(snapshot["psi_out"], quad)
        np.testing.assert_array_equal(actual, snapshot["psi_in"])

    def test_realizer_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Realized bare :class:`AngularAverageOperator` matches at nulp=4.

        Wave 1 bit-equivalence is the underlying guarantee; nulp=4 is
        the documented safety margin in case a future ``ScaledOperator``
        wrap or compositional reorder shifts a single ULP. The current
        ship-state passes ``assert_array_equal``; this test would only
        relax if Wave 7+ introduces such a reordering.
        """
        quad = LevelSymmetricSN.create(sn_order=4)
        bc = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        actual = op.apply(snapshot["psi_out"])
        np.testing.assert_array_almost_equal_nulp(
            actual, snapshot["psi_in"], nulp=4,
        )


# ─────────────────────────────────────────────────────────────────────
# Case 6 — White(xmin, α=0.3) × GaussLegendre1D(8)
# ─────────────────────────────────────────────────────────────────────


class TestWhiteXminPartial03GLSnapshot:
    """White (Lambertian) xmin at α=0.3 on a 1-D GL quadrature.

    Body identical to Case 5 (cosine-weighted average), but with a
    ``ScaledOperator(0.3, ...)`` wrap on the realizer side. Legacy
    multiplies inside the body (``psi_avg[None, ...] * self.albedo``);
    the realized op multiplies on the outside. Reduction tree is the
    same; bit-exact in practice on this case, ``nulp=4`` documents
    the allowance for one-ULP drift if downstream waves reorder.

    GL is selected here because it exercises a 1-D adapter
    (different ``Sigma w`` than spherical Lebedev / level-symmetric)
    — see numerical-bug-signatures Signature 4 (quadrature-dependent
    constant hardcoded). This case is the canary against a future
    refactor introducing a ``4*pi``-hardcode regression.
    """

    case_id = "white_xmin_partial_03_GL"

    @pytest.fixture(scope="class")
    def snapshot(self) -> np.lib.npyio.NpzFile:
        return _load_or_skip(self.case_id)

    def test_legacy_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Legacy ``bc.apply`` matches snapshot bit-for-bit (the pin)."""
        quad = GaussLegendre1D.create(n_ordinates=8)
        bc = WhiteBoundaryOperator(axis="x", outward_sign=-1, albedo=0.3)
        actual = bc.apply(snapshot["psi_out"], quad)
        np.testing.assert_array_equal(actual, snapshot["psi_in"])

    def test_realizer_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Realized ``ScaledOperator(0.3, AngularAverageOperator)`` at nulp=4.

        The plan allows ``nulp=4`` for this case because the legacy
        does ``psi_avg * albedo`` (multiplication-then-broadcast) and
        the realizer does ``ScaledOperator(0.3, base).apply`` which is
        ``base.apply(psi) * 0.3`` (broadcast-then-multiplication). The
        floating-point result is identical for these N values; nulp=4
        is the safety margin against future reordering.
        """
        quad = GaussLegendre1D.create(n_ordinates=8)
        bc = WhiteBoundaryOperator(axis="x", outward_sign=-1, albedo=0.3)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        actual = op.apply(snapshot["psi_out"])
        np.testing.assert_array_almost_equal_nulp(
            actual, snapshot["psi_in"], nulp=4,
        )


# ─────────────────────────────────────────────────────────────────────
# Case 7 — Periodic × Lebedev(17)
# ─────────────────────────────────────────────────────────────────────


class TestPeriodicLebedev17Snapshot:
    """Periodic BC: identity body (``psi_out.copy()``). Bit-exact."""

    case_id = "periodic_lebedev17"

    @pytest.fixture(scope="class")
    def snapshot(self) -> np.lib.npyio.NpzFile:
        return _load_or_skip(self.case_id)

    def test_legacy_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Legacy ``psi_out.copy()`` is bit-exact."""
        quad = LebedevSphere.create(17)
        bc = PeriodicBoundaryOperator()
        actual = bc.apply(snapshot["psi_out"], quad)
        np.testing.assert_array_equal(actual, snapshot["psi_in"])

    def test_realizer_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Realized :class:`PeriodicWrapOperator` is bit-exact.

        Body is a no-op copy; no FP operations, no reduction tree.
        ``assert_array_equal`` is the right gate.
        """
        quad = LebedevSphere.create(17)
        bc = PeriodicBoundaryOperator()
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        actual = op.apply(snapshot["psi_out"])
        np.testing.assert_array_equal(actual, snapshot["psi_in"])


# ─────────────────────────────────────────────────────────────────────
# Case 8 — Mixed (0.3 Specular + 0.7 White) × LevelSymmetric(4)
# ─────────────────────────────────────────────────────────────────────


class TestMixed30Spec70WhiteLS4Snapshot:
    r"""0.3-Specular + 0.7-White mixed boundary on LS_4 — Wave-11 form.

    Realises the rank-N tensor decomposition

    .. math::

        R = 0.3 \, G_{\text{refl}} + 0.7 \, G_{\text{diff}}

    via the Wave-0 algebra: each leaf
    :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` is realised
    independently against the SN method space, and the rank-N
    composition is expressed as a Wave-0
    :class:`~orpheus.numerics.operator.OperatorSum` of
    :class:`~orpheus.numerics.operator.ScaledOperator` wrappers.

    Wave 11 removed the ``MixedBoundaryOperator`` composer, so the
    pre-Wave-11 legacy-path assertion (which exercised
    ``MixedBoundaryOperator.apply``) is gone. Only the realiser path
    survives — the snapshot ``psi_in`` is the realised
    ``apply(psi)`` of the Wave-0 composition; see
    :func:`tests.geometry._generate_bc_equivalence_snapshots._build_mixed_realized_apply`.

    The composed-vs-snapshot assertion stays at ``nulp=64`` (the FP
    safety margin for the OperatorSum reduction-tree change vs the
    pre-Wave-11 snapshot if it shifts by a ULP across the renamed
    field-access path; on the current ship-state both producer and
    consumer use bit-identical Wave-0 reduction so the actual
    distance is 0). The ``vv-principles`` bit-identity vs
    principled-equivalence three criteria are satisfied: each
    intermediate is a named Wave-0 type, the realised value is
    verified against the explicit pointwise weighted sum in
    :func:`tests.geometry.test_boundary.test_operator_sum_of_bcs_acts_as_weighted_sum`,
    and any drift is FP non-associativity over a 2-summand reduction.
    """

    case_id = "mixed_30spec_70white_LS4"

    @pytest.fixture(scope="class")
    def snapshot(self) -> np.lib.npyio.NpzFile:
        return _load_or_skip(self.case_id)

    def test_realizer_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Wave-0 ``0.3 * spec_realised + 0.7 * white_realised`` matches
        snapshot at ``nulp=64``."""
        quad = LevelSymmetricSN.create(sn_order=4)
        spec_realized = SNBoundaryRealizer().realize(
            SpecularBoundaryOperator(axis="x", albedo=1.0),
            SNMethodSpace.minimal(quad),
        )
        white_realized = SNBoundaryRealizer().realize(
            WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0),
            SNMethodSpace.minimal(quad),
        )
        composed = 0.3 * spec_realized + 0.7 * white_realized
        actual = composed.apply(snapshot["psi_out"])
        np.testing.assert_array_almost_equal_nulp(
            actual, snapshot["psi_in"], nulp=64,
        )
