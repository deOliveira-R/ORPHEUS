r"""The NEW spelling of the reference-carry population: face roles plus an ``exterior`` of ``WallAmplitudes`` (wave 2, S5).

Not a test module. :mod:`tests.gates.derivations._reference_carry_cases` states each body's walls as data; the HEAD capture
spelled them as ``BC(...)``; this module spells the SAME declarations the way the carved tree does (rules R1 to R3 of
``.claude/plans/wave2_boundary_carve.md`` 6.0), so the two spellings are two readers of one table:

* a mirror on a Cartesian face is ``MirrorDeck()`` and a periodic one ``TranslationDeck()`` (a symmetry of the domain, no entry);
* every other declaration on a Cartesian face, and EVERY declaration on a curved face (a mirror on a curved face is a tag bound
  to the specular wall ``(1, 0)``, the meaning change of Appendix B), is an integer tag;
* tags are numbered ``0, 1, ...`` in face order, one per DISTINCT declaration (R1: equal declarations share a tag, so a slab
  with ``vacuum | vacuum`` carries the one tag ``0``), and ``exterior`` binds each tag to the amplitudes
  :func:`~tests.gates.derivations._reference_carry_cases.amplitudes` states by hand.

``orpheus`` is imported inside the functions: a module-level import would make a collection error of every gate that uses the
spelling, where each gate must record its own red line.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from tests.gates.derivations import _reference_carry_cases as cases

if TYPE_CHECKING:
    from orpheus.geometry import StructuredGeometry
    from orpheus.specification.exterior import WallAmplitudes


def roles_and_exterior(chart: str, faces: tuple[cases.Decl, ...]) -> tuple[tuple[Any, ...], "dict[int, WallAmplitudes]"]:
    """The face roles (deck or tag, inner first) and the ``exterior`` of one body's declarations."""
    from orpheus.geometry import MirrorDeck, TranslationDeck
    from orpheus.specification.exterior import WallAmplitudes

    tags: dict[tuple, int] = {}
    exterior: dict[int, WallAmplitudes] = {}
    roles: list[Any] = []
    for decl in faces:
        if chart == "slab" and decl[0] == "mirror":
            roles.append(MirrorDeck())
        elif chart == "slab" and decl[0] == "periodic":
            roles.append(TranslationDeck())
        else:
            if decl not in tags:
                tags[decl] = len(tags)
                exterior[tags[decl]] = WallAmplitudes(*cases.amplitudes(decl))
            roles.append(tags[decl])
    return tuple(roles), exterior


def geometry_and_exterior(chart: str, breakpoints: tuple[float, ...], faces: tuple[cases.Decl, ...]) -> "tuple[StructuredGeometry, dict[int, WallAmplitudes]]":
    """The carved-tree geometry of a case and the exterior its tags state."""
    from orpheus.geometry import StructuredGeometry

    roles, exterior = roles_and_exterior(chart, faces)
    mats = cases.mat_ids(breakpoints)
    if chart == "slab":
        return StructuredGeometry.slab(breakpoints, mats, left=roles[0], right=roles[1]), exterior
    build = StructuredGeometry.sphere if chart == "sphere" else StructuredGeometry.cylinder
    if len(roles) == 1:
        return build(breakpoints, mats, outer=roles[0]), exterior
    return build(breakpoints, mats, inner=roles[0], outer=roles[1]), exterior
