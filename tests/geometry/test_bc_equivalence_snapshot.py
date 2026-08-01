r"""BC equivalence snapshot harness (Wave 6 / Issue #186 cleanup).

This file is the safety net for the boundary-operator refactor (see
``.claude/plans/transient-giggling-cake.md`` for the 12-wave plan and
``.claude/plans/bc-trace-law-descriptor-cleanup.md`` for the Issue #186
descriptor-model cleanup that removed the legacy 2-arg call path).

Wave-6 ship-state had two assertions per case: a legacy
``bc.apply(psi_out, quadrature)`` half pinning the pre-refactor body,
and a realised ``SNBoundaryRealizer().realize(...).apply(psi_out)``
half pinning the modern path. Issue #186 (B3 + β2) removed the
``apply`` method from every concrete BC descriptor, eliminating the
legacy half. Only the realiser-path assertion remains; the snapshots
themselves now record realiser-path outputs.

For each of the 8 representative BC × quadrature × input cases captured
by :mod:`tests.geometry._generate_bc_equivalence_snapshots`, one
assertion is pinned per pytest invocation:

* **Realized 1-arg path.** ``SNBoundaryRealizer().realize(bc,
  SNMethodSpace.minimal(quadrature)).apply(psi_out)`` reproduces the
  snapshot at the per-case tolerance — protecting against post-Issue-#186
  numerical drift in the realiser bodies.

The vacuum case asserts the narrowed law's SHAPE and its VALUE: since
campaign phase B3.2 the realized law is typed ``Γ₊ → Γ₋``, so vacuum is
the honest zero map and there is no "outflow pass-through" half left to
check — those rows are outside the operator's domain. (Pre-B3.2 this
case carried a 2-assertion structure, a zero-mask check AND a
pass-through check, because the §16A.5 realization was a full-face
projector that preserved the outflow rows. The shape assertion replaced
the second half, and it is load-bearing: |Γ₊| == |Γ₋| on every
quadrature in the tree, so shape alone cannot catch a law that stayed
an endomorphism — only paired with the value can it.)

V&V tags
--------

* ``@pytest.mark.foundation`` + ``@pytest.mark.regression`` — a
  frozen-reference drift gate on a software invariant (realiser-body
  stability), NOT a verification reference.
* No ``@pytest.mark.verifies`` decorator — these tests verify a
  *software invariant* (refactor stability), not an equation.

**B0.3 RELABEL (2026-07-30).** This file was marked ``@pytest.mark.l1``
and its docstring claimed *"cross-implementation L1 (the realizer path
is a structurally-independent re-implementation; the snapshot pins the
numerical reference)"*. That claim died with the legacy half: Issue
#186 removed the descriptor ``apply`` bodies, and — as the paragraph
above already recorded — **the snapshots now record realiser-path
outputs**. A baseline regenerated from the very code it gates is not a
cross-implementation reference at any level; under ``vv-principles``
(L4 is parallel to the ladder, and a self-generated baseline is not
even that) an ``l1`` label on it is level conflation: it advertises
correctness evidence while producing none.

The level marker is therefore ``foundation`` — the honest home for
"software invariant, no theory-page ``:label:``" — and ``regression``
stays, which is what the pytest marker table already calls this shape
(*"Frozen-reference snapshot tests gating numerical drift across
refactors (NOT a verification reference)"*).

**Nothing about the file's VALUE changes:** it remains the widest
mutation net in the boundary subsystem (measured: it reddened under
**9 of the 12** leaf-action mutations in the boundary review). What
changes is only the claim it makes about itself.

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
    AlbedoBoundary,
    PeriodicBoundary,
    ReflectiveBoundary,
    VacuumInflow,
    WhiteBoundary,
)
from orpheus.sn.boundary.realizer import SNBoundaryRealizer
from orpheus.sn.mesh.method_space import SNMethodSpace
from tests.sn._test_helpers import face_method_space

#: The mixed-law snapshot below sums reflective with **white**. Both are
#: narrowed (B3.2, B3.4a), so the rank-N composition claim is stateable again
#: and the ``xfail(strict=True)`` that stood here was deleted at B3.4a — as
#: its own reason text instructed. It must NOT be re-posed over albedo or
#: periodic: those are still full-face (B3.4b / B3.4c) and, because
#: ``|Γ₊| == |Γ₋|`` everywhere, such a sum does not RAISE — it runs silently
#: wrong (vv Mode 12), so the gate would be green and worthless.
from orpheus.numerics.quadrature import Quadrature


pytestmark = [pytest.mark.foundation, pytest.mark.regression]


_SNAPSHOT_DIR = Path(__file__).parent / "snapshots"
_GENERATOR_HINT = (
    "Run `python -m tests.geometry._generate_bc_equivalence_snapshots` "
    "to create."
)


def _load_snapshot(case_id: str) -> np.lib.npyio.NpzFile:
    """Load the snapshot for ``case_id``; a missing file is a FAILURE.

    B0.3 REPAIR — this helper used to ``pytest.skip`` on a missing
    snapshot, on the rationale that "the test infrastructure lands
    before the snapshots themselves in a decoupled-commit roll-out".
    That roll-out finished long ago (all 8 ``.npz`` are git-tracked),
    and the skip is now purely a way for the widest net in the
    subsystem to be **silently disabled** by a deleted or renamed
    baseline — the same shape as the skip-swallowed sentinel this
    review found next door (``vv-principles`` Mode-8, the
    SKIP-SWALLOWED class: a skip is invisible in a summary line).

    A missing baseline is a broken gate, so it reds.
    """
    path = _SNAPSHOT_DIR / f"bc_equivalence_{case_id}.npz"
    if not path.exists():
        pytest.fail(
            f"Snapshot {path.name} is MISSING — this gate is the widest "
            f"mutation net in the boundary subsystem (9 of 12 leaf "
            f"mutations) and a missing baseline silently disables it. "
            f"{_GENERATOR_HINT} Regeneration must carry an inline "
            f"justification citing one of the three "
            f"principled-equivalence criteria (vv-principles).",
        )
    return np.load(path)


# ─────────────────────────────────────────────────────────────────────
# Case 1 — Vacuum × Lebedev(17)
# ─────────────────────────────────────────────────────────────────────


class TestVacuumLebedev17Snapshot:
    """Vacuum BC on Lebedev-17 quadrature.

    Since B3.2 the realized law is typed ``Γ₊ → Γ₋``, so vacuum is the
    ZERO MAP: it consumes the outflow half-trace and emits |Γ₋| rows of
    zeros. The snapshot's ``inflow_indices_xmin`` index set is the
    independent statement of which ordinates those are, cross-checked
    against the live derivation before either is used.
    """

    case_id = "vacuum_lebedev17"

    @pytest.fixture(scope="class")
    def snapshot(self) -> np.lib.npyio.NpzFile:
        return _load_snapshot(self.case_id)

    def test_realizer_zeroes_only_inflow_per_section_16A5(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """The narrowed vacuum law is the zero map ``Γ₊ → Γ₋``.

        Pre-B3.2 this asserted the §16A.5 semantics — an
        ``IncomingOrdinateMaskTensor`` zeroing the inflow rows while the
        outflow rows passed through. The narrowing removed the second
        half rather than changing it: the preserved rows left the
        operator's domain, so vacuum's whole content is ``R = 0``.

        The snapshot records the inflow index set so this assertion has
        a deterministic target without re-deriving it at test time.
        """
        quad = Quadrature.lebedev(17)
        space = face_method_space(quad, face="xmin")
        # The frozen snapshot's own inflow index set is the independent
        # statement of WHICH ordinates are inflow at xmin; cross-check the
        # live derivation against it before using either.
        np.testing.assert_array_equal(
            space.inflow_indices, snapshot["inflow_indices_xmin"],
        )
        op = SNBoundaryRealizer().realize(VacuumInflow(), space)
        psi_out = snapshot["psi_out"]
        actual = op.apply(psi_out[space.outflow_indices])
        assert actual.shape == (space.inflow_indices.size,) + psi_out.shape[1:], (
            f"vacuum emitted {actual.shape}; the narrowed codomain is Γ₋."
        )
        np.testing.assert_array_equal(
            actual, 0.0,
            err_msg="the narrowed vacuum law is the ZERO map Γ₊ → Γ₋",
        )


# ─────────────────────────────────────────────────────────────────────
# Case 2 — Albedo(0.5) × Lebedev(17)
# ─────────────────────────────────────────────────────────────────────


class TestAlbedo05Lebedev17Snapshot:
    """Pure albedo: ``psi_in = 0.5 * psi_out``.

    Realized as :class:`ScaledOperator(0.5, IdentityOperator())`. The
    snapshot ``psi_in`` is the realiser-path output; the test pins
    bit-equality against it.
    """

    case_id = "albedo_05_lebedev17"

    @pytest.fixture(scope="class")
    def snapshot(self) -> np.lib.npyio.NpzFile:
        return _load_snapshot(self.case_id)

    def test_realizer_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Realized ``ScaledOperator(0.5, IdentityOperator)`` is bit-exact.

        Multiplication is the only operation; no reduction tree change.
        ``assert_array_equal`` is the right gate.
        """
        quad = Quadrature.lebedev(17)
        bc = AlbedoBoundary(albedo=0.5)
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
        return _load_snapshot(self.case_id)

    def test_realizer_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Realized bare :class:`PermutationOperator` is bit-exact.

        At α=1 the realizer's fast path returns the bare
        :class:`PermutationOperator` (no :class:`ScaledOperator` wrap).
        """
        quad = Quadrature.lebedev(17)
        bc = ReflectiveBoundary(axis="x", albedo=1.0)
        space = face_method_space(quad, face="xmax")
        op = SNBoundaryRealizer().realize(bc, space)
        actual = op.apply(snapshot["psi_out"][space.outflow_indices])
        # ⭐ THE BIT-IDENTITY CLAIM, at the snapshot: the frozen ``psi_in`` IS
        # the pre-B3.2 full-face image, so the narrowed law's image must be
        # that same array RESTRICTED to Γ₋ — bit-for-bit, no tolerance. The
        # reference therefore comes from the FROZEN artefact, never from
        # re-running the new code.
        np.testing.assert_array_equal(
            actual, snapshot["psi_in"][space.inflow_indices],
        )


# ─────────────────────────────────────────────────────────────────────
# Case 4 — Specular(y, α=0.7) × LevelSymmetric(6)
# ─────────────────────────────────────────────────────────────────────


class TestSpecularYPartial07LS6Snapshot:
    """Specular y-axis with partial albedo α=0.7.

    Multiplication-then-permutation realises as
    :class:`ScaledOperator(0.7, PermutationOperator)`. Reduction tree
    is multiplication-then-gather; ``assert_array_equal`` is the right
    gate.
    """

    case_id = "specular_y_partial_07_LS6"

    @pytest.fixture(scope="class")
    def snapshot(self) -> np.lib.npyio.NpzFile:
        return _load_snapshot(self.case_id)

    def test_realizer_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Realized ``ScaledOperator(0.7, PermutationOperator)`` is bit-exact."""
        quad = Quadrature.level_symmetric(sn_order=6)
        bc = ReflectiveBoundary(axis="y", albedo=0.7)
        space = face_method_space(
            quad, face="ymax", faces=("xmin", "xmax", "ymin", "ymax"),
        )
        op = SNBoundaryRealizer().realize(bc, space)
        actual = op.apply(snapshot["psi_out"][space.outflow_indices])
        # Bit-identity against the FROZEN pre-B3.2 image, restricted to Γ₋
        # (see the α=1 sibling). This row additionally carries the SCALED
        # path, so an α-fold regression is caught at the same time.
        np.testing.assert_array_equal(
            actual, snapshot["psi_in"][space.inflow_indices],
        )


# ─────────────────────────────────────────────────────────────────────
# Case 5 — White(xmax, α=1) × LevelSymmetric(4)
# ─────────────────────────────────────────────────────────────────────


class TestWhiteXmaxLS4Snapshot:
    """White (Lambertian) xmax at α=1.

    The cosine-weighted hemisphere average is an N-term reduction. The
    realiser's fast path returns the bare
    :class:`AngularAverageOperator` at α=1, bit-exact against the
    snapshot.
    """

    case_id = "white_xmax_LS4"

    @pytest.fixture(scope="class")
    def snapshot(self) -> np.lib.npyio.NpzFile:
        return _load_snapshot(self.case_id)

    def test_realizer_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Realized bare :class:`AngularAverageOperator` is bit-exact.

        ``assert_array_equal`` is the right gate; the snapshot was
        generated from this exact path.
        """
        quad = Quadrature.level_symmetric(sn_order=4)
        bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
        space = face_method_space(
            quad, face="xmax", faces=("xmin", "xmax", "ymin", "ymax"),
        )
        op = SNBoundaryRealizer().realize(bc, space)
        actual = op.apply(snapshot["psi_out"][space.outflow_indices])
        # ⭐ B3.4a — the SAME bit-identity shape the specular cases already
        # use: the frozen ``psi_in`` is the pre-narrowing FULL-FACE image, so
        # the narrowed law's image must be that array RESTRICTED to Γ₋,
        # bit-for-bit. The snapshot is NOT regenerated; the reference stays the
        # frozen artefact, never a re-run of the new code.
        #
        # `[M]` It holds EXACTLY here (maxdiff 0.0). It need not in general:
        # the narrowed sum runs over the restricted array instead of a
        # zero-padded full-N one, and that reduction-order change is worth ~1
        # ULP on quadratures where it bites. `level_symmetric(4)` carries ZERO
        # tangential ordinates, so nothing was mis-admitted and nothing moved.
        np.testing.assert_array_equal(
            actual, snapshot["psi_in"][space.inflow_indices],
        )


# ─────────────────────────────────────────────────────────────────────
# Case 6 — White(xmin, α=0.3) × GaussLegendre1D(8)
# ─────────────────────────────────────────────────────────────────────


class TestWhiteXminPartial03GLSnapshot:
    """White (Lambertian) xmin at α=0.3 on a 1-D GL quadrature.

    Body identical to Case 5 (cosine-weighted average), but with a
    :class:`ScaledOperator(0.3, ...)` wrap. GL is selected here
    because it exercises a 1-D adapter (different ``Σ w`` than
    spherical Lebedev / level-symmetric) — see ``numerical-bug-signatures``
    Signature 4 (quadrature-dependent constant hardcoded). This case
    is the canary against a future refactor introducing a
    ``4*pi``-hardcode regression.
    """

    case_id = "white_xmin_partial_03_GL"

    @pytest.fixture(scope="class")
    def snapshot(self) -> np.lib.npyio.NpzFile:
        return _load_snapshot(self.case_id)

    def test_realizer_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Realized ``ScaledOperator(0.3, AngularAverageOperator)`` is bit-exact."""
        quad = Quadrature.gauss_legendre(n_ordinates=8)
        bc = WhiteBoundary(axis="x", outward_sign=-1, albedo=0.3)
        space = face_method_space(quad, face="xmin")
        op = SNBoundaryRealizer().realize(bc, space)
        actual = op.apply(snapshot["psi_out"][space.outflow_indices])
        # B3.4a — frozen pre-narrowing image restricted to Γ₋; see Case 5 for
        # why the artefact is not regenerated. `[M]` maxdiff 0.0 here too: a
        # 1-D GL quadrature carries no tangential ordinates. This row
        # additionally carries the SCALED path, so an α-fold regression is
        # caught at the same time.
        np.testing.assert_array_equal(
            actual, snapshot["psi_in"][space.inflow_indices],
        )


# ─────────────────────────────────────────────────────────────────────
# Case 7 — Periodic × Lebedev(17)
# ─────────────────────────────────────────────────────────────────────


class TestPeriodicLebedev17Snapshot:
    """Periodic BC: identity body. Bit-exact."""

    case_id = "periodic_lebedev17"

    @pytest.fixture(scope="class")
    def snapshot(self) -> np.lib.npyio.NpzFile:
        return _load_snapshot(self.case_id)

    def test_realizer_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Realized :class:`PeriodicWrapOperator` is bit-exact.

        Body is a no-op copy; no FP operations, no reduction tree.
        ``assert_array_equal`` is the right gate.
        """
        quad = Quadrature.lebedev(17)
        bc = PeriodicBoundary()
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        actual = op.apply(snapshot["psi_out"])
        np.testing.assert_array_equal(actual, snapshot["psi_in"])


# ─────────────────────────────────────────────────────────────────────
# Case 8 — Mixed (0.3 Specular + 0.7 White) × LevelSymmetric(4)
# ─────────────────────────────────────────────────────────────────────


class TestMixed30Spec70WhiteLS4Snapshot:
    r"""0.3-Specular + 0.7-White mixed boundary on LS_4.

    Realises the rank-N tensor decomposition

    .. math::

        R = 0.3 \, G_{\text{refl}} + 0.7 \, G_{\text{diff}}

    via the Wave-0 operator algebra: each leaf
    :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` is realised
    independently against the SN method space, and the rank-N
    composition is expressed as a Wave-0
    :class:`~orpheus.numerics.operator.OperatorSum` of
    :class:`~orpheus.numerics.operator.ScaledOperator` wrappers.

    The composed-vs-snapshot assertion stays at ``nulp=64`` (the FP
    safety margin for the OperatorSum reduction-tree change vs the
    pre-Wave-11 snapshot if it shifts by a ULP across the renamed
    field-access path).
    """

    case_id = "mixed_30spec_70white_LS4"

    @pytest.fixture(scope="class")
    def snapshot(self) -> np.lib.npyio.NpzFile:
        return _load_snapshot(self.case_id)

    def test_realizer_apply_matches_snapshot(
        self, snapshot: np.lib.npyio.NpzFile,
    ) -> None:
        """Wave-0 ``0.3 * spec_realised + 0.7 * white_realised`` matches
        snapshot at ``nulp=64``.

        RE-POSED at **B3.2** onto Γ₊ (the narrowed leaf's domain) and against
        the frozen image restricted to Γ₋ — the same bit-identity shape as the
        pure-specular cases above. **UNBLOCKED at B3.4a**: white narrowed, so
        both summands are now typed ``Γ₊ → Γ₋`` and the rank-N composition
        claim is stateable. Both leaves are realized on the SAME face-ful
        space, which is what makes the sum well-typed — realizing one of them
        on a faceless space is no longer even possible.
        """
        quad = Quadrature.level_symmetric(sn_order=4)
        space = face_method_space(
            quad, face="xmax", faces=("xmin", "xmax", "ymin", "ymax"),
        )
        spec_realized = SNBoundaryRealizer().realize(
            ReflectiveBoundary(axis="x", albedo=1.0), space,
        )
        white_realized = SNBoundaryRealizer().realize(
            WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0), space,
        )
        composed = 0.3 * spec_realized + 0.7 * white_realized
        actual = composed.apply(snapshot["psi_out"][space.outflow_indices])
        np.testing.assert_array_almost_equal_nulp(
            actual, snapshot["psi_in"][space.inflow_indices], nulp=64,
        )
