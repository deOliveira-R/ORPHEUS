r"""The carried problems of wave 2, as DATA (slice S2b, gate 3; plan section 6.0 "The carry fixture").

Not a test module (no ``test_`` prefix). Two readers share this file and nothing else:

* ``_capture_carry_head.py`` (deleted in WF2.4; recover it with ``git show 197549f6:tests/gates/transport/_capture_carry_head.py``) builds each problem in the SPELLING OF ``main`` ``5d04c296`` (the typed laws of
  ``orpheus.geometry.boundary``) inside a detached worktree of that commit and freezes the answers into
  ``carry_head.npz``;
* ``tests/gates/sn/operators/test_wave2_carry.py`` builds the SAME problems in the new spelling (a tag bound to a
  :class:`WallCondition`, or a deck) and requires the answers bit-equal.

So this module imports nothing from ``orpheus`` except at call time and holds no boundary vocabulary: a declaration is
a neutral tuple, ``("vacuum",)``, ``("mirror",)``, ``("specular", a)``, ``("white", a)``, ``("uniform_inflow", v)``, and
each reader translates it once. That is what keeps the two sides from drifting into different problems (a case whose
law and quadrature were written out twice, once per side, is the Pattern-2 twin the geometry harness names).

Every setting the answer depends on is a field here (the body, the cell counts, the quadrature, the solver
tolerances), because a bit-equal ``k`` is meaningful only if the solver is the same solver: nothing is left to a
default that a later change could move.

Populations (X2): 4 bodies, 6 + 5 + 5 + 2 declarations; each declaration is solved as a k-eigenvalue problem and as a
fixed-source problem with a bulk source (the slab also with a uniform inflow); the realized operator matrix of both
faces of an S16 slab is frozen for each slab declaration.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Any, TypedDict

#: The commit the fixture is captured at (section 6.0).
HEAD_COMMIT = "5d04c296"

#: The frozen HEAD answers (``_capture_carry_head.py`` writes it).
CARRY_FIXTURE = Path(__file__).with_name("carry_head.npz")

#: The seven independent ``bc_equivalence`` references (R-10: NOT re-baselined; the directory moved from
#: ``tests/gates/geometry/snapshots/`` to here in WF2.4, with the retirement of the law tests that read it).
SNAPSHOT_DIR = Path(__file__).resolve().parent.parent / "sn" / "operators" / "snapshots"

#: ``(kind, value)`` neutral declarations. ``mirror`` is the perfect mirror (HEAD ``BC.reflective``).
Declaration = tuple


@dataclasses.dataclass(frozen=True)
class Body:
    """A 2-group two-region body: ``mat_ids`` read from the centre (or the left face) outward."""

    name: str
    coord: str  # "CARTESIAN" | "CYLINDRICAL" | "SPHERICAL"
    breakpoints: tuple[float, ...]
    mat_ids: tuple[int, ...]
    cells: tuple[int, ...]
    quadrature: tuple  # ("gauss_legendre", n) | ("folded_product", n_mu, n_phi)

    @property
    def hollow(self) -> bool:
        return self.coord != "CARTESIAN" and self.breakpoints[0] > 0.0

    @property
    def faces(self) -> tuple[str, ...]:
        """The boundary faces the declaration is installed on (the tag side of the new spelling)."""
        if self.coord == "CARTESIAN":
            return ("xmin", "xmax")
        return ("xmin", "xmax") if self.hollow else ("xmax",)


#: Material 1 is the fuel ("A" of the 2-group library), material 0 the moderator ("B").
BODIES: dict[str, Body] = {
    "slab": Body("slab", "CARTESIAN", (0.0, 1.5, 3.0), (1, 0), (6, 6), ("gauss_legendre", 8)),
    "cylinder": Body("cylinder", "CYLINDRICAL", (0.0, 1.0, 2.0), (1, 0), (5, 5), ("folded_product", 4, 8)),
    "sphere": Body("sphere", "SPHERICAL", (0.0, 1.0, 2.0), (1, 0), (5, 5), ("gauss_legendre", 8)),
    "hollow_cylinder": Body("hollow_cylinder", "CYLINDRICAL", (0.5, 1.25, 2.0), (1, 0), (4, 4), ("folded_product", 4, 8)),
}

#: The OUTER wall declarations per body. On the hollow cylinder the inner face is always the reflective cavity
#: (``("mirror",)``, which SN computes today); on the slab the declaration is installed on both faces.
DECLARATIONS: dict[str, tuple[Declaration, ...]] = {
    "slab": (("vacuum",), ("mirror",), ("specular", 0.3), ("specular", 0.7), ("specular", 1.0), ("white", 0.5), ("white", 1.0)),
    "cylinder": (("vacuum",), ("mirror",), ("specular", 0.3), ("specular", 0.7)),
    "sphere": (("vacuum",), ("mirror",), ("specular", 0.3), ("specular", 0.7)),
    "hollow_cylinder": (("vacuum",), ("specular", 0.7)),
}

#: The slab fixed-source problems under a uniform inflow on the left face (the right face vacuum).
INFLOW_VALUE = 1.0

#: The S16 slab whose realized face operators are frozen.
S16_ORDINATES = 16

#: Solver settings, spelled out (nothing left to a default).
class _EigenSettings(TypedDict):
    max_outer: int
    keff_tol: float
    flux_tol: float
    max_inner: int
    inner_tol: float


class _FixedSettings(TypedDict):
    max_inner: int
    inner_tol: float


EIGEN = _EigenSettings(max_outer=500, keff_tol=1e-10, flux_tol=1e-9, max_inner=4000, inner_tol=1e-11)
FIXED = _FixedSettings(max_inner=4000, inner_tol=1e-12)


def tag(declaration: Declaration) -> str:
    """A key fragment: ``vacuum``, ``mirror``, ``specular_0.3`` ..."""
    return declaration[0] if len(declaration) == 1 else f"{declaration[0]}_{declaration[1]}"


def key(body: str, declaration: Declaration, what: str) -> str:
    """The npz key of one frozen array: ``<what>/<body>/<declaration>``."""
    return f"{what}/{body}/{tag(declaration)}"


def inflow_key(what: str) -> str:
    return f"{what}/slab/uniform_inflow_{INFLOW_VALUE}"


def face_key(declaration: Declaration, face: str) -> str:
    return f"face_matrix/slab_s16/{tag(declaration)}/{face}"


def materials() -> dict[int, Any]:
    """``{1: fuel, 0: moderator}``, the 2-group library mixtures (the data layer, unchanged by the carve)."""
    from orpheus.data.macro_xs.synthetic import get_mixture

    return {1: get_mixture("A", "2g"), 0: get_mixture("B", "2g")}


def quadrature(body: Body, ordinates: int | None = None) -> Any:
    from orpheus.numerics.quadrature import Quadrature

    kind, *params = body.quadrature
    if kind == "gauss_legendre":
        return Quadrature.gauss_legendre(n_ordinates=ordinates if ordinates is not None else params[0])
    return Quadrature.folded_product(n_mu=params[0], n_phi=params[1])


def bulk_source(n_ordinates: int, ng: int, shape: tuple[int, ...]) -> Any:
    """The per-ordinate bulk source ``(N, ng, *spatial)``: ones in group 0, zero in group 1."""
    import numpy as np

    q = np.zeros((n_ordinates, ng, *shape))
    q[:, 0] = 1.0
    return q


def declared_cases() -> list[tuple[str, Declaration]]:
    """Every ``(body, declaration)`` pair, in the order the fixture stores them."""
    return [(name, d) for name, ds in DECLARATIONS.items() for d in ds]


# ── the other methods' carried numbers (slices S3 and S4 gate them; captured here, from the one HEAD tree) ──────────

#: ``(chart, breakpoints)`` of the diffusion and collision-probability bodies: a two-region 2-group body, 6 cells each.
CHARTS: dict[str, tuple[str, tuple[float, ...]]] = {
    "slab": ("CARTESIAN", (0.0, 1.5, 3.0)),
    "cylinder": ("CYLINDRICAL", (0.0, 1.0, 2.0)),
    "sphere": ("SPHERICAL", (0.0, 1.0, 2.0)),
}
OTHER_CELLS = (6, 6)
DIFFUSION_DECLARATIONS: tuple[Declaration, ...] = (("vacuum",), ("mirror",), ("albedo", 0.5), ("zero_flux",))
#: CP: the slab's left face is a forced mirror (section 4.2), so the slab is declared mirror on the left and the
#: declaration on the right; the curved bodies carry the declaration on their one outer face.
CP_DECLARATIONS: tuple[Declaration, ...] = (("vacuum",), ("white", 1.0))
MOC_CELLS = (5, 5)
MOC_ARGS = dict(n_azi=8, n_polar=2, ray_spacing=0.1)
#: MC: the slab periodic on both faces, materials 2 (fuel) and 0 (moderator), a fixed seed.
MC_PARAMS = dict(n_neutrons=200, n_inactive=5, n_active=20, pitch=3.0, seed=12345)


def other_key(method: str, chart: str, declaration: Declaration) -> str:
    return f"keff_{method}/{chart}/{tag(declaration)}"
