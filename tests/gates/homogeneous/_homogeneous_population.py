r"""The ONE homogeneous population: the eight producing mixtures the tree ships.

Shared by every homogeneous gate that runs over the population
(``test_kinf_exact_reference``, ``test_coda_anchors``,
``test_operator_spaces`` G2.1, ``test_homogeneous_population``), so there is
one list (``coding-elegance`` Pattern 2). It lived in
``test_byte_stability.py`` until that gate retired on 2026-10-01: its byte
comparison pinned the platform's LAPACK output (a macOS update moved
``geev``'s k∞ by 1 ULP on an unchanged tree), and its successor is the exact
rational reference with a derived bound, ``test_kinf_exact_reference.py``.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations import get
from orpheus.derivations.common.xs_library import get_mixture

#: The eg-bearing variant's grid (the one homogeneous eg idiom in the tree,
#: ``test_homogeneous.py`` eg-block).
_EDGES_2G = np.array([1.0e7, 1.0e3, 1.0e-3])

EXPECTED_CASES = (
    "homo_1eg", "homo_2eg", "homo_2eg_n2n", "homo_2eg_with_eg", "homo_4eg",
    "mixture_A_1g", "mixture_A_2g", "mixture_A_4g",
)


def mixture_cases() -> dict[str, Mixture]:
    """name -> producing Mixture, exhaustive over what the tree ships.

    The ONE homogeneous population (Pattern 2): gate (c)
    ``test_kinf_exact_reference``, ``test_coda_anchors`` and
    ``test_operator_spaces`` G2.1 all read it. Regions B/C/D of the xs
    library are non-producing (``[M]`` the ``is_producing`` screen at
    ``24a991ba``), so the eigenvalue entry is meaningless there.
    """
    cases: dict[str, Mixture] = {}
    for name in ("homo_1eg", "homo_2eg", "homo_4eg", "homo_2eg_n2n"):
        cases[name] = next(iter(get(name).materials.values()))
    cases["homo_2eg_with_eg"] = dataclasses.replace(
        next(iter(get("homo_2eg").materials.values())), eg=_EDGES_2G
    )
    for k in ("1g", "2g", "4g"):
        cases[f"mixture_A_{k}"] = get_mixture("A", k)
    return cases


