r"""The absorbing-chain predicate: which energy classes of a closed domain can never lose a neutron.

Space-membership campaign, unit U0 (spec ``scratch/membership_covector/u0/spec_u0.md`` rows P1-P5, by the
test-architect). :func:`~orpheus.transport.absorbing_chain.absorption_free_closed_classes` reads DECLARED data only (the
present mixtures and the binding) and returns the witness, the absorption-free final classes of the group-transfer graph
on a closed domain, never a bool. It mints no basis: the left co-vector depends on the assembled operator, so each
method mints its bases from its own loss on the classes named here.

The values were validated on a prototype against the solvers (``scratch/membership_covector/spec/battery/
check_u0_predicate.py``: 11 of 11 posings classified as measured), and against a dense SVD in the method gates.
"""
from __future__ import annotations

import pytest

from orpheus.geometry.deck import MirrorDeck
from orpheus.transport.absorbing_chain import absorption_free_closed_classes as _classes_of, final_classes
from orpheus.transport.boundary_conditions import WallCondition
from orpheus.transport.boundary_conditions.binding import bind

from tests.gates.transport._singular_loss_fixtures import (
    ABSORBER, BOTH_2G, FISSILE_NO_CAPTURE, G2_ONLY, G2_UPSCATTER, NEAR_1E14, PURE_1G, mix, slab,
)

pytestmark = [pytest.mark.l0]


def absorption_free_closed_classes(mixtures, binding):
    """The predicate's classes as member tuples (the witness each row compares)."""
    return tuple(c.members for c in _classes_of(mixtures, binding))


@pytest.mark.parametrize("mixture,classes", [
    (PURE_1G, ((0,),)), (BOTH_2G, ((1,),)), (G2_ONLY, ((1,),)),
    (G2_UPSCATTER, ()), (NEAR_1E14, ()), (FISSILE_NO_CAPTURE, ()), (ABSORBER, ()),
], ids=["pure1g", "both2g", "g2_only", "upscatter", "sigma_a_1e-14", "fissile_no_capture", "absorber"])
def test_P1_the_absorption_free_final_classes_of_an_infinite_medium(mixture, classes):
    """``g2_only`` is the row that refutes "Sigma_a = 0 everywhere": one absorbing group does not save the closed class
    downstream of it. ``sigma_a_1e-14`` is regular (an exact-zero predicate; the conditioning is recorded, not judged)."""
    assert absorption_free_closed_classes([mixture], None) == classes


def test_P2_a_class_with_removal_in_one_region_is_regular():
    """Every present material must be absorption-free on the class; one absorbing region makes it transient."""
    mesh = slab(bps=(0.0, 5.0, 10.0), mats=(0, 1), n=(4, 4))
    assert absorption_free_closed_classes([ABSORBER, PURE_1G], bind(mesh, {})) == ()
    assert absorption_free_closed_classes([PURE_1G], bind(mesh, {})) == ((0,),)


@pytest.mark.parametrize("condition,closed", [
    (WallCondition.white(1.0), True), (WallCondition.specular(1.0), True),
    (WallCondition.white(0.999), False), (WallCondition.vacuum(), False),
], ids=["white1", "specular1", "white0.999", "vacuum"])
def test_P3_a_face_conserves_iff_its_return_fraction_is_exactly_one(condition, closed):
    """The abstract verb ``BoundaryCondition.return_fraction``, not ``isinstance(WallCondition)`` (diffusion's zero-flux
    condition, return fraction -1, is the diffusion gate's row); not SN's realized ``closes_axis`` either: white(1.0)
    conserves without closing an axis."""
    binding = bind(slab(0, 0), {0: condition})
    assert (absorption_free_closed_classes([PURE_1G], binding) == ((0,),)) is closed


def test_P4_one_leaking_face_opens_the_domain():
    binding = bind(slab(MirrorDeck(), 0), {0: WallCondition.vacuum()})
    assert absorption_free_closed_classes([PURE_1G], binding) == ()


def test_P5_reads_the_loss_removal_not_absorption_xs():
    """A mixture whose total is not the sum of its channels (SigT = 1, every channel 0: the SN placeholder material's
    shape) has absorption_xs = 0 and removal 1, so it is REGULAR. First red of an absorption_xs-keyed predicate: 527
    tests in 66 files flagged (spec §7, ``battery/u0_census_absorption_xs.json``)."""
    placeholder = mix([1.0], [0.0], [[0.0]])
    assert float(placeholder.absorption_xs[0]) == 0.0
    assert absorption_free_closed_classes([placeholder], None) == ()


def test_the_final_classes_of_a_chain_with_a_cycle():
    """The graph primitive on its own: 0 -> 1 <-> 2 -> 3, and 4 isolated: final classes {3} and {4}; {1, 2} is a
    communicating class that is not final."""
    import numpy as np

    t = np.zeros((5, 5), dtype=bool)
    t[0, 1] = t[1, 2] = t[2, 1] = t[2, 3] = True
    assert final_classes(t) == ((3,), (4,))
    t[3, 1] = True
    assert final_classes(t) == ((1, 2, 3), (4,))
