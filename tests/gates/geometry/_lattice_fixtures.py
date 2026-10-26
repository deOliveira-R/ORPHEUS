"""Lattice fixtures shared by the geometry and the numerics gates (one home, X4).

Geometry-pure: numpy and :mod:`orpheus.geometry` only, so a numerics test may import
it (numerics knows geometry; the reverse edge is forbidden).  Re-filed 2026-10-25 by the
layer-DAG campaign (ruling 23, the geometry -> numerics move): every helper here was
read by a test that moved to ``tests/gates/numerics/`` AND by one that stayed in
``tests/gates/geometry/``, so it could neither travel nor stay.  ``_g22b_AXES`` and
``_g22b_mirror`` were pasted character-identically into ``test_symmetry.py`` and
``test_manifold.py``; the other four lived in one of the two.  The bodies are verbatim.
"""

from __future__ import annotations

import numpy as np

from orpheus.geometry.symmetry import SubgroupOfO3
from orpheus.geometry.transformation import RigidMotion


_g22b_AXES = ("x", "y", "z")


def _g22b_mirror(axis: str) -> RigidMotion:
    return RigidMotion.reflection(normal=np.eye(3)["xyz".index(axis)])


def _every_spellable_group() -> list[SubgroupOfO3]:
    named = [
        SubgroupOfO3.Trivial,
        SubgroupOfO3.Dinfh,
        SubgroupOfO3.OctahedralOh,
        SubgroupOfO3.IcosahedralIh,
        SubgroupOfO3.SO3,
        SubgroupOfO3.O3,
    ]
    named += [SubgroupOfO3.Mirror(a) for a in _g22b_AXES]
    named += [SubgroupOfO3.SO2(a) for a in _g22b_AXES]
    named += [SubgroupOfO3.O2(a) for a in _g22b_AXES]
    named += [SubgroupOfO3.Cn(n) for n in (1, 2, 3, 4)]
    named += [SubgroupOfO3.Dnh(n) for n in (1, 2, 3, 4)]
    return named


#: Every shipped orbit space of the sphere: the six catalogue keys plus the
#: derived identity quotient.  A universal below is a universal over THIS.
_ROSTER = [
    SubgroupOfO3.O2("x"), SubgroupOfO3.O2("y"), SubgroupOfO3.O2("z"),
    SubgroupOfO3.Mirror("x"), SubgroupOfO3.Mirror("y"), SubgroupOfO3.Mirror("z"),
    SubgroupOfO3.Trivial,
]


def _generic_unit_directions(seed: int = 20260902, n: int = 9) -> np.ndarray:
    """Deterministic generic unit directions — no component zero, none related by a coordinate symmetry."""
    raw = np.random.default_rng(seed).normal(size=(n, 3))
    return raw / np.linalg.norm(raw, axis=1)[:, None]


_R4_I3 = np.eye(3)
