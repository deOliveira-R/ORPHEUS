r"""The population of the reference carry: bodies and the walls they state, as DATA (wave 2, S5; plan 6.5 gates 2 and 3).

Not a test module, and importing nothing from ``orpheus``: the SAME table is read by the HEAD capture script
(``_capture_reference_carry_head.py``, which spells each wall as the retired string tag ``BC(...)`` and is run once, from a
detached tree of ``main`` ``5d04c296``) and by the post-carve tests (``test_reference_carry.py`` and
``test_walls_of_reads_decks_and_amplitudes.py``, which spell the same wall as a face role plus an ``exterior`` of
``WallAmplitudes``). One table, two spellings, one fixture between them: that is the whole instrument. The fixture they compare
against is ``reference_carry_head.json``.

A wall declaration is a tuple, one per boundary point, inner first (a slab: left, right; a solid radial body: outer only):

* ``("vacuum",)`` nothing returned, :math:`(s, d) = (0, 0)`;
* ``("mirror",)`` everything returned along the reflected line, :math:`(1, 0)` (a ``MirrorDeck`` on a slab, a tag on a curved face);
* ``("partial", a)`` the fraction ``a`` returned specularly, :math:`(a, 0)`;
* ``("white", a)`` the fraction ``a`` returned isotropically, :math:`(0, a)`;
* ``("periodic",)`` a wrap to the opposite wall (a ``TranslationDeck``; Cartesian only, a curved wrap is a ``DeckError``).

The materials are the two-group Sood-registry mixtures (``PU_2_0_SL_STUB`` downscatter, ``URRB_2_0_IN`` upscatter), alternating by
interval, so the reference solves a multi-group heterogeneous body (``vv-principles`` anti-patterns 3 and 4: a one-group or
homogeneous population nulls the terms a wrong wall reading moves).
"""

from __future__ import annotations

from typing import Final

Decl = tuple  # ("vacuum",) | ("mirror",) | ("partial", a) | ("white", a) | ("periodic",)

#: ``(name, chart, breakpoints, faces)``; ``mat_ids`` alternate 0, 1, 0, ... over the intervals.
CASES: Final[tuple[tuple[str, str, tuple[float, ...], tuple[Decl, ...]], ...]] = (
    ("slab-vacuum-vacuum", "slab", (0.0, 0.4, 1.5, 2.3), (("vacuum",), ("vacuum",))),
    ("slab-mirror-vacuum", "slab", (0.0, 0.4, 1.5, 2.3), (("mirror",), ("vacuum",))),
    ("slab-mirror-mirror", "slab", (0.0, 0.4, 1.5, 2.3), (("mirror",), ("mirror",))),
    ("slab-partial-vacuum", "slab", (0.0, 0.4, 1.5, 2.3), (("partial", 0.37), ("vacuum",))),
    ("slab-partial-partial", "slab", (0.0, 0.4, 1.5, 2.3), (("partial", 0.5), ("partial", 0.25))),
    ("slab-white-vacuum", "slab", (0.0, 0.4, 1.5, 2.3), (("white", 0.6), ("vacuum",))),
    ("slab-white-mirror", "slab", (0.0, 0.4, 1.5, 2.3), (("white", 1.0), ("mirror",))),
    ("slab-periodic", "slab", (0.0, 0.4, 1.5, 2.3), (("periodic",), ("periodic",))),
    ("sphere-vacuum", "sphere", (0.0, 0.5, 1.5, 2.0), (("vacuum",),)),
    ("sphere-mirror", "sphere", (0.0, 0.5, 1.5, 2.0), (("mirror",),)),
    ("sphere-partial", "sphere", (0.0, 0.5, 1.5, 2.0), (("partial", 0.6),)),
    ("sphere-white", "sphere", (0.0, 0.5, 1.5, 2.0), (("white", 0.8),)),
    ("cylinder-vacuum", "cylinder", (0.0, 0.6, 1.0), (("vacuum",),)),
    ("cylinder-mirror", "cylinder", (0.0, 0.6, 1.0), (("mirror",),)),
    ("cylinder-white", "cylinder", (0.0, 0.6, 1.0), (("white", 1.0),)),
    ("hollow-sphere-vacuum-vacuum", "sphere", (0.4, 0.5, 1.5, 2.0), (("vacuum",), ("vacuum",))),
    ("hollow-sphere-partial-mirror", "sphere", (0.4, 0.5, 1.5, 2.0), (("partial", 0.3), ("mirror",))),
    ("hollow-sphere-white-vacuum", "sphere", (0.4, 0.5, 1.5, 2.0), (("white", 0.5), ("vacuum",))),
    ("hollow-cylinder-mirror-white", "cylinder", (0.3, 0.6, 1.0), (("mirror",), ("white", 0.7))),
)

#: The cases whose every wall is vacuum or specular, the population ``specular_albedos`` serves (no diffuse wall, no wrap).
SPECULAR_ONLY: Final[tuple[str, ...]] = tuple(
    name for name, _chart, _bp, faces in CASES if all(d[0] in ("vacuum", "mirror", "partial") for d in faces)
)

#: The albedo each declaration is read as by ``specular_albedos``, written by hand.
ALBEDO_OF: Final[dict[str, float]] = {"vacuum": 0.0, "mirror": 1.0}


def albedo(decl: Decl) -> float:
    """The specular albedo a vacuum, mirror or partial declaration states, by hand."""
    return float(decl[1]) if decl[0] == "partial" else ALBEDO_OF[decl[0]]


def amplitudes(decl: Decl) -> tuple[float, float]:
    """The cone coordinates ``(specular, diffuse)`` of a non-wrap declaration, by hand."""
    kind = decl[0]
    if kind == "vacuum":
        return (0.0, 0.0)
    if kind == "mirror":
        return (1.0, 0.0)
    if kind == "partial":
        return (float(decl[1]), 0.0)
    if kind == "white":
        return (0.0, float(decl[1]))
    raise ValueError(f"a {kind!r} wall has no amplitudes: it is a wrap")


def mat_ids(breakpoints: tuple[float, ...]) -> tuple[int, ...]:
    """Alternating 0, 1, 0, ... over the intervals."""
    return tuple(k % 2 for k in range(len(breakpoints) - 1))
