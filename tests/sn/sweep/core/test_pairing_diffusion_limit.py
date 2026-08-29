r"""Phase 1a (#236 ST4): the (spatial ⊗ angular) diffusion-limit pairing predicate.

``@foundation`` — software-invariant tests of the pairing-validity predicate
(:func:`orpheus.sn.sweep.pairing.pair_diffusion_limit_consistent`) and the
literature-declared per-scheme / per-closure booleans.  The booleans transcribe
Larsen–Morel–Miller 1987 (spatial) / Bailey–Morel–Chang 2010 (angular); the
underlying PHYSICS is verified elsewhere
(``tests/sn/verification/mms/test_mms_ld_slab.py::test_ld_thick_diffusive_limit``
for LD, the curvilinear MMS gates for M-M).  These tests pin the transcription
(a flipped boolean is a regression) and the predicate logic — positive AND
negative coverage per ``vv-principles`` anti-pattern #11.

NO ``@verifies`` — a pairing-validity predicate is a software invariant, not an
equation discretization (vv-principles: ``foundation`` carries no theory label).
"""

from types import SimpleNamespace

import pytest

from orpheus.transport.spatial.diamond import DiamondDifference
from orpheus.transport.spatial.linear_discontinuous import LinearDiscontinuous
from orpheus.sn.sweep.pairing import pair_diffusion_limit_consistent
from orpheus.sn.angular.closure import (
    IdentityAngularClosure,
    MorelMontryAngularSweep,
    AngularClosureBase,
)
from orpheus.transport.spatial.scheme import DiscretizationSchemeBase


@pytest.mark.foundation
def test_spatial_scheme_diffusion_limit_booleans() -> None:
    """Spatial trait transcribes LMM-1987 / Larsen–Morel 1989: DD and full LD True."""
    # DD: LMM-1987 Eq. (4.24) leading-order limit.
    assert DiamondDifference.diffusion_limit_consistent is True
    # Full LD (slope source threaded, D5b-S3): Larsen–Morel 1989 II Eq. (4.16).
    assert LinearDiscontinuous.diffusion_limit_consistent is True


@pytest.mark.foundation
def test_angular_closure_beta_booleans() -> None:
    """Angular trait transcribes BMC-2010: M-M True (Eq. 42); Cartesian identity vacuously True."""
    assert MorelMontryAngularSweep.beta_first_order_consistent is True
    assert IdentityAngularClosure.beta_first_order_consistent is True


@pytest.mark.foundation
def test_base_defaults_are_conservative_opt_in() -> None:
    """A scheme/closure that declares nothing is NOT assumed consistent (opt-in default)."""
    assert DiscretizationSchemeBase.diffusion_limit_consistent is False
    assert AngularClosureBase.beta_first_order_consistent is False


@pytest.mark.foundation
def test_pair_predicate_truth_table() -> None:
    """The predicate is the conjunction — positive AND every negative branch (vv #11).

    Stubs (``SimpleNamespace``) exercise the predicate LOGIC across the full
    2×2 truth table without instantiating a mesh.  The (False, *) rows are the
    negative coverage anti-pattern #11 requires; ``Step × M-M`` will be the
    first REAL (False, True) pairing once Step is implemented (#158).
    """
    def scheme(ok: bool) -> SimpleNamespace:
        return SimpleNamespace(diffusion_limit_consistent=ok)

    def closure(ok: bool) -> SimpleNamespace:
        return SimpleNamespace(beta_first_order_consistent=ok)

    assert pair_diffusion_limit_consistent(scheme(True), closure(True)) is True
    # spatial fails (e.g. Step × M-M) → pair invalid even with a good closure
    assert pair_diffusion_limit_consistent(scheme(False), closure(True)) is False
    # angular fails (e.g. LD × step-in-angle on a curvilinear mesh) → pair invalid
    assert pair_diffusion_limit_consistent(scheme(True), closure(False)) is False
    assert pair_diffusion_limit_consistent(scheme(False), closure(False)) is False


@pytest.mark.foundation
def test_pair_predicate_on_canonical_real_pairings() -> None:
    """The real canonical pairs are valid; the predicate reads class-level traits.

    ``DiamondDifference × IdentityAngularClosure`` is the Cartesian default;
    ``LinearDiscontinuous × MorelMontryAngularSweep`` is the diffusion-limit-
    consistent curvilinear pairing (the #233 floor-lifting target).  Passing the
    closure CLASSES (rather than mesh-bound instances) suffices because the
    predicate reads only the ClassVar trait.
    """
    assert pair_diffusion_limit_consistent(
        DiamondDifference(), IdentityAngularClosure,
    ) is True
    assert pair_diffusion_limit_consistent(
        LinearDiscontinuous(), MorelMontryAngularSweep,
    ) is True


@pytest.mark.foundation
def test_pair_predicate_reads_a_real_scheme_instance() -> None:
    """The predicate honours a REAL scheme instance's trait against a stub closure."""
    assert pair_diffusion_limit_consistent(
        LinearDiscontinuous(), SimpleNamespace(beta_first_order_consistent=True),
    ) is True
    # A real LD spatial scheme cannot rescue a β-inconsistent angular closure.
    assert pair_diffusion_limit_consistent(
        LinearDiscontinuous(), SimpleNamespace(beta_first_order_consistent=False),
    ) is False
