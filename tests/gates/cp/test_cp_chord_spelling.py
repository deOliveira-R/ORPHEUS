r"""G4 of the unit ``insulation`` (#641 item 1): CP's half-chords are the partition's chord, spelled once.

``CPMesh._setup_radial_quadrature`` reads the outbound half of ``ConcentricPartition.chord(line).slot_length`` (regions ``0..n-1``,
transposed to ``(n, n_y)``) where it read ``orpheus.derivations.common.kernels.chord_half_lengths(edges[1:], y)``.

What this gate proves and does not. The agreement is TAUTOLOGICAL in its upstream (X4): ``_radial_chord``, the body the partition's chord
runs, forms its slots from the same ``region_half_lengths`` of the same ``half_chord`` that ``chord_half_lengths`` calls. So the chord swap has
no behavioural first red: this gate lands green at the base tip, and its first red is a test-only input, the mutation. What it pins is the
SLOT CHOICE and the TRANSPOSITION (the outbound half, regions ``0..n-1``, ``.T``): the slot layout decoded by position
(``[..., -n:]``) is the producer's convention read at the consumer. The behavioural pin of the chord is the carry row ``sphere-white_1.0``
(``test_cp_boundary_carry``), and the arm below that runs the carry file with the inbound half unreversed.

Hollow bodies. A hollow ``CPMesh`` cannot be built (#513 refuses the declared inner face, ``BoundaryRefused``), so no hollow ``k`` exists to
capture; the function-level leg discharges the question: for breakpoints ``[0.1, 0.3, 0.5]`` on the cylinder and the sphere and seven
values of ``y`` in ``[0.1, 0.5]`` the partition's chord equals ``chord_half_lengths`` bit for bit.

Mutations (in-process, by ``monkeypatch``): M-inbound-slots (``Chord.slot_length`` returns the array with the inbound half ``[..., :n]``,
unreversed, selected where the outbound was; the comparison leg reads red).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from orpheus.derivations.common.kernels import chord_half_lengths
from orpheus.geometry.chart import Chart
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.line import Line
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates.cp import _cp_s4_helpers as CP

_CHARTS = ["cylinder", "sphere"]
_COORD = {"cylinder": CoordSystem.CYLINDRICAL, "sphere": CoordSystem.SPHERICAL}


def _cp_mesh(chart: str) -> Any:
    return CP.cp_head(CP.cp_mesh(chart), {0: WallCondition.vacuum()})


@pytest.mark.l1
@pytest.mark.parametrize("chart", _CHARTS)
def test_the_cp_chords_equal_the_reference_half_chords_bit_for_bit(chart: str) -> None:
    r"""``mesh._chords`` equals ``chord_half_lengths(edges[1:], mesh._y_pts)`` (``array_equal``) and has shape ``(12, 768)`` (12 regions,
    12 panels x 64 points per panel) on the carried cylinder and sphere. The comparison is tautological upstream (module docstring).

    Red input at the base tip: none behavioural (green, by X4). Mutation: M-inbound-slots.
    """
    mesh = _cp_mesh(chart)
    expected = chord_half_lengths(mesh.mesh.edges[1:], mesh._y_pts)
    assert mesh._chords.shape == expected.shape == (12, 768)
    assert np.array_equal(mesh._chords, expected)


@pytest.mark.l1
@pytest.mark.parametrize("chart", _CHARTS)
def test_the_hollow_function_level_leg(chart: str) -> None:
    r"""A hollow body is not constructible as a ``CPMesh`` (#513), so its leg is the function: for breakpoints ``[0.1, 0.3, 0.5]`` and
    seven values of ``y`` in ``[0.1, 0.5]`` the outbound slots of the partition's chord equal ``chord_half_lengths(breakpoints[1:], y)``
    bit for bit (the equality needs ``y >= r_0``, and the y-rule starts at ``edges[0]``). Green at the base tip (X4); pins the slot choice."""
    breakpoints = (0.1, 0.3, 0.5)
    partition = ConcentricPartition(Chart(_COORD[chart]), breakpoints)
    y = np.linspace(0.1, 0.5, 7)
    n = partition.n_regions
    points = np.stack([y, np.zeros_like(y), np.zeros_like(y)], axis=-1)
    line = Line.through(points, np.broadcast_to(np.array([0.0, 1.0, 0.0]), points.shape))
    got = partition.chord(line).slot_length[..., -n:].T
    assert np.array_equal(got, chord_half_lengths(np.asarray(breakpoints[1:]), y))


def _inbound_slot_chord(n: int) -> Any:
    """``Chord.slot_length`` with the outbound slots ``[..., -n:]`` holding the inbound half ``[..., :n]``, unreversed (M-inbound-slots)."""
    from orpheus.geometry import chord as chord_module

    original = chord_module.Chord.slot_length

    def inbound(self: Any) -> Any:
        a = np.array(original.fget(self))  # type: ignore[union-attr]
        a[..., -n:] = a[..., :n]
        return a

    return property(inbound)


@pytest.mark.foundation
def test_the_slot_choice_reddens_when_the_inbound_slots_are_read(monkeypatch: pytest.MonkeyPatch) -> None:
    r"""M-inbound-slots at function level: the hollow leg's comparison is red under the mutated ``Chord.slot_length``. The positive control
    is the honest read of the same line. Green at the base tip (it touches the geometry type only)."""
    from orpheus.geometry import chord as chord_module

    breakpoints = (0.1, 0.3, 0.5)
    partition = ConcentricPartition(Chart(CoordSystem.CYLINDRICAL), breakpoints)
    y = np.linspace(0.1, 0.5, 7)
    points = np.stack([y, np.zeros_like(y), np.zeros_like(y)], axis=-1)
    line = Line.through(points, np.broadcast_to(np.array([0.0, 1.0, 0.0]), points.shape))
    n = partition.n_regions
    reference = chord_half_lengths(np.asarray(breakpoints[1:]), y)
    assert np.array_equal(partition.chord(line).slot_length[..., -n:].T, reference)
    monkeypatch.setattr(chord_module.Chord, "slot_length", _inbound_slot_chord(n))
    assert not np.array_equal(partition.chord(line).slot_length[..., -n:].T, reference)


@pytest.mark.l1
@pytest.mark.parametrize("chart", _CHARTS)
def test_the_cp_chords_are_read_from_the_partitions_chord(chart: str, monkeypatch: pytest.MonkeyPatch) -> None:
    r"""The activation count of the carve: the built ``mesh._chords`` MOVES when ``Chord.slot_length`` returns the inbound half in the
    outbound slots. Red input at the base tip, for the real reason: ``CPMesh`` still reads ``chord_half_lengths`` from the reference, so the
    mutation reaches nothing and ``_chords`` equals the reference (the gate lands red and goes green when the carve routes the path)."""
    from orpheus.geometry import chord as chord_module

    honest = _cp_mesh(chart)
    expected = chord_half_lengths(honest.mesh.edges[1:], honest._y_pts)
    monkeypatch.setattr(chord_module.Chord, "slot_length", _inbound_slot_chord(len(honest.mesh.edges) - 1))
    mutant = _cp_mesh(chart)
    assert not np.array_equal(mutant._chords, expected), "CPMesh does not read the partition's chord: the slot choice is unpinned"
