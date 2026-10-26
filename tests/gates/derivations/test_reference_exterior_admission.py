r"""Every reference that holds a ``(geometry, exterior)`` pair admits it through ONE door (wave 2, S5 amend; review findings 1, 6, 7, 8).

The exact-cover predicate on the exterior (``exterior.keys() == geometry.tags``) is
:func:`orpheus.specification.exterior.admit_exterior`, written once. Before the amend it ran only inside
``GeometrySpecification``, so the three reference spaces (``MomentSpace``, ``Spectrum``, ``BasisSpace``), ``Walls.of`` and
``specular_albedos`` held a raw mapping and the states the specification refuses could still be built: a missing tag failed as a
bare ``KeyError`` with no owner, an extra entry was accepted silently, a caller's later mutation of the mapping changed the
object after it had been validated.

What each row pins, with the mutation that reddens it (both described, and run in a copy at review time):

* **One door, by tag number** (``test_a_missing_or_unused_tag_is_refused_by_number_at_every_consumer``): five consumers, two
  refusals each. Mutation: delete the ``admit_exterior`` call from any one consumer, and that consumer's two cells redden (a
  ``KeyError``, and a silent build).
* **The copy is frozen** (``test_the_exterior_is_held_as_a_frozen_copy``): the caller mutates the dict it passed and the space's
  ``exterior`` does not move. Mutation: store ``self.exterior`` as passed.
* **One reader of a role's amplitudes** (``test_amplitudes_at_reads_a_deck_as_itself_and_a_tag_through_the_exterior``) and the
  two readers agree on a mirror (``test_the_two_readers_agree_on_a_mirror``). Mutation: write the mirror as ``(0.5, 0)`` in
  either reader.
* **Named vertices** (``test_the_vertices_of_the_triangle_are_named_once``).

Positive controls: the same five consumers accept the valid cover (``test_a_valid_cover_builds_at_every_consumer``), so a
refusal above is the cover's and not a broken construction (``instrument-doctrine`` X1).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.reference_body import amplitudes_at, specular_albedos
from orpheus.data.macro_xs.synthetic import make_mixture
from orpheus.derivations.continuous.characteristic.walls import Walls
from orpheus.derivations.continuous.fn_method.moment_space import MomentSpace
from orpheus.derivations.continuous.galerkin_spectral.basis_space import BasisSpace
from orpheus.derivations.continuous.singular_eigenfunction.spectrum import Spectrum
from orpheus.geometry import CoordSystem, MirrorDeck, StructuredGeometry, TranslationDeck
from orpheus.specification.exterior import WallAmplitudes, admit_exterior

_MIX = make_mixture(
    sig_t=np.array([1.0]), sig_c=np.array([0.0]), sig_f=np.array([0.3]), nu=np.array([2.0]),
    chi=np.array([1.0]), sig_s=np.array([[0.7]]),
)
_SLAB = StructuredGeometry(coord=CoordSystem.CARTESIAN, breakpoints=(0.0, 2.0), mat_ids=(0,), face_roles=(0, 0))


def _moment(exterior):
    return MomentSpace(geometry=_SLAB, exterior=exterior, materials={0: _MIX})


def _spectrum(exterior):
    return Spectrum(geometry=_SLAB, exterior=exterior, materials={0: _MIX})


def _basis(exterior):
    return BasisSpace(geometry=_SLAB, exterior=exterior, materials={0: _MIX})


def _walls(exterior):
    return Walls.of(_SLAB, exterior)


def _albedos(exterior):
    return specular_albedos(_SLAB, exterior, owner="test")


_CONSUMERS = {"MomentSpace": _moment, "Spectrum": _spectrum, "BasisSpace": _basis, "Walls.of": _walls, "specular_albedos": _albedos}


@pytest.mark.parametrize("name", _CONSUMERS)
def test_a_valid_cover_builds_at_every_consumer(name: str) -> None:
    """Positive control: tag 0 bound to the vacuum wall is the exact cover of the slab's one tag."""
    assert _CONSUMERS[name]({0: WallAmplitudes.VACUUM}) is not None


@pytest.mark.parametrize("name", _CONSUMERS)
def test_a_missing_or_unused_tag_is_refused_by_number_at_every_consumer(name: str) -> None:
    build = _CONSUMERS[name]
    with pytest.raises(ValueError, match=r"states no amplitudes for the tag\(s\) \[0\]"):
        build({})
    with pytest.raises(ValueError, match=r"states the tag\(s\) \[9\] that no face of the geometry carries"):
        build({0: WallAmplitudes.VACUUM, 9: WallAmplitudes.MIRROR})
    with pytest.raises(TypeError, match="a wall is a WallAmplitudes"):
        build({0: (0.0, 0.0)})


@pytest.mark.parametrize("name", ["MomentSpace", "Spectrum", "BasisSpace"])
def test_the_exterior_is_held_as_a_frozen_copy(name: str) -> None:
    passed = {0: WallAmplitudes.VACUUM}
    held = _CONSUMERS[name](passed)
    passed[0] = WallAmplitudes.MIRROR
    assert held.exterior[0] == WallAmplitudes.VACUUM


def test_the_owner_prefixes_the_refusal() -> None:
    with pytest.raises(ValueError, match="^SomeReader: the exterior"):
        admit_exterior({}, _SLAB, owner="SomeReader")


def test_amplitudes_at_reads_a_deck_as_itself_and_a_tag_through_the_exterior() -> None:
    exterior = {0: WallAmplitudes(0.25, 0.0)}
    assert amplitudes_at(MirrorDeck(), exterior) == WallAmplitudes.MIRROR
    assert amplitudes_at(TranslationDeck(), exterior) == WallAmplitudes.MIRROR
    assert amplitudes_at(0, exterior) == WallAmplitudes(0.25, 0.0)
    with pytest.raises(AssertionError):
        amplitudes_at("not a role", exterior)  # type: ignore[arg-type]


def test_the_two_readers_agree_on_a_mirror() -> None:
    mirror_slab = StructuredGeometry(coord=CoordSystem.CARTESIAN, breakpoints=(0.0, 2.0), mat_ids=(0,), face_roles=(MirrorDeck(), MirrorDeck()))
    walls = Walls.of(mirror_slab, {})
    assert [w.amplitudes for w in walls.walls] == [WallAmplitudes.MIRROR] * 2
    assert specular_albedos(mirror_slab, {}, owner="test") == (1.0, 1.0)


def test_the_vertices_of_the_triangle_are_named_once() -> None:
    assert WallAmplitudes.VACUUM == WallAmplitudes(0.0, 0.0)
    assert WallAmplitudes.MIRROR == WallAmplitudes(1.0, 0.0)
    assert WallAmplitudes.WHITE == WallAmplitudes(0.0, 1.0)
