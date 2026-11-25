r"""The fixtures of the bit-identity harness: five problem families at three group counts, and 21 dense-metric Gram families.

Not a test module (no ``test_`` prefix), and it imports no other test module: what the harness needs from the gates is copied
here, each copy with a comment naming its source ``file:line``. A copied literal or function is a definition the harness
owns, so a later edit of a gate cannot move the stored bytes. Nothing here is a copy of an algorithm of ``orpheus``: every
problem is built through the public constructors, and the constructor call of each is the text returned beside the object.

The extents are NON-dyadic on purpose. A power-of-two cell width makes the product order of two diagonal slots invisible in the
bytes (``x * a * b == x * b * a`` exactly when ``a`` is a power of two), so the slot-order control of
``test_snapshot_controls.py`` could not see a reversed ``FactoredMetric._walk`` on the slab and the box. A later unit that
re-meshes the harness keeps the extents non-dyadic.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

from orpheus.data.macro_xs.synthetic import get_mixture, make_mixture
from orpheus.diffusion import DiffusionMesh, DiffusionSolver
from orpheus.geometry import CoordSystem, StructuredGeometry
from orpheus.geometry.deck import DeckTransformation, MirrorDeck
from orpheus.mesh import CellEdges, CellsByCount, Mesh2D, Mesher
from orpheus.mesh.axis import roles_of_axes
from orpheus.numerics.declaration import Coord, Declaration, Spanned
from orpheus.numerics.quadrature import Quadrature
from orpheus.values.content import content_digest
from orpheus.sn.problem import SNProblem
from orpheus.transport.boundary_conditions import WallCondition

GROUPS = (1, 2, 3)
SN_KINDS = ("slab", "box", "sphere", "cylinder")

# The sizes of the proposed small configuration (the prototype ``scratch/ontology/w1/open/N00/capture.py``, ``SIZES_SMALL``).
SLAB_REGION_CELLS = (2, 2, 2)  # cells per region, for the slab, the sphere and the cylinder
SLAB_THICKNESSES = (0.5, 0.9, 0.5)  # non-dyadic; the sphere and the cylinder reuse them as radii
GAUSS_LEGENDRE_ORDER = 4  # slab and sphere
CYLINDER_FOLDED = (2, 4)  # ``Quadrature.folded_product(n_mu, n_phi)``
BOX_CELLS = (2, 2)
BOX_EXTENTS = (1.3, 0.7)  # non-dyadic
BOX_LEVEL_SYMMETRIC_ORDER = 2
DIFFUSION_CELL_EDGES = ((0.0, 0.3, 0.5), (0.5, 1.4, 3.0), (3.0, 4.1, 5.0))  # non-dyadic; region edges 0.5 | 3.0 | 5.0
SPECULAR_ALBEDO = 1.0  # the albedos of the wall conditions of the fixtures; the texts of ``snapshots/BUILDERS.txt`` read them
WHITE_ALBEDO = 0.3

# ─────────────────────────── copied from the gates ───────────────────────────

# tests/gates/diffusion/test_solver.py:56-71 (the 2-group and 3-group tables).
_SIG_T_A = np.array([0.2181, 0.7850])
_SIG_S_A = np.array([[0.1900, 0.0160], [0.0, 0.4200]])
_SIG_F_A = np.array([0.0024, 0.0489])
_SIG_T_B = np.array([0.3416, 0.9431])
_SIG_S_B = np.array([[0.1000, 0.0020], [0.0, 0.0500]])
_SIG_T_3G = np.array([0.30, 0.55, 1.10])
_SIG_S_3G = np.array([
    [0.10, 0.05, 0.01],
    [0.00, 0.30, 0.08],
    [0.00, 0.00, 0.60],
])
_SIG_F_3G = np.array([0.002, 0.010, 0.120])
_NU_3G = np.array([2.9, 2.5, 2.4])
_CHI_3G = np.array([0.7, 0.3, 0.0])


def _fuel_2g():
    # tests/gates/diffusion/test_solver.py:74
    return make_mixture(
        sig_t=_SIG_T_A,
        sig_c=_SIG_T_A - _SIG_F_A - _SIG_S_A.sum(axis=1),
        sig_f=_SIG_F_A, nu=np.array([2.54, 2.47]),
        chi=np.array([1.0, 0.0]), sig_s=_SIG_S_A,
    )


def _reflector_2g():
    # tests/gates/diffusion/test_solver.py:83
    return make_mixture(
        sig_t=_SIG_T_B,
        sig_c=_SIG_T_B - _SIG_S_B.sum(axis=1),
        sig_f=np.zeros(2), nu=np.zeros(2),
        chi=np.zeros(2), sig_s=_SIG_S_B,
    )


def _fuel_3g():
    # tests/gates/diffusion/test_solver.py:92
    return make_mixture(
        sig_t=_SIG_T_3G,
        sig_c=_SIG_T_3G - _SIG_F_3G - _SIG_S_3G.sum(axis=1),
        sig_f=_SIG_F_3G, nu=_NU_3G,
        chi=_CHI_3G, sig_s=_SIG_S_3G,
    )


def _het_materials():
    # tests/gates/diffusion/test_solver.py:117
    return {0: _fuel_2g(), 1: _reflector_2g()}


def declare_sides(*sides: Any) -> tuple[tuple[Any, ...], dict[int, Any]]:
    """The roles and the ``{tag: condition}`` of faces declared by value, in face order (``tests/gates/_exterior.py:30``)."""
    roles: list[Any] = []
    conditions: dict[int, Any] = {}
    for side in sides:
        if isinstance(side, DeckTransformation):
            roles.append(side)
            continue
        tag = next((t for t, held in conditions.items() if held == side), None)
        if tag is None:
            tag = len(conditions)
            conditions[tag] = side
        roles.append(tag)
    return tuple(roles), conditions


def conditions_for(carrier: Any, table: dict[int, Any]) -> dict[int, Any]:
    """The entries of ``table`` for the tags ``carrier``'s faces carry (``tests/gates/_exterior.py:51`` and ``:57``)."""
    roles = carrier.face_roles if hasattr(carrier, "face_roles") else roles_of_axes(tuple(carrier))
    tags = tuple(sorted({int(role) for role in roles.values() if not isinstance(role, DeckTransformation)}))
    missing = [tag for tag in tags if tag not in table]
    if missing:
        raise KeyError(f"the table declares {sorted(table)} but the carrier's faces carry the tags {list(tags)}: no condition for {missing}")
    return {tag: table[tag] for tag in tags}


# tests/gates/numerics/test_quotient_metric.py:93 and :96. The rotation runs through LAPACK, so the copied code reproduces the
# matrix only on a host with the same stamp (the first line of ``snapshots/MANIFEST.txt``).
_Q5 = np.linalg.qr(np.random.default_rng(0).standard_normal((5, 5)))[0]


def _gram(eigs):
    G = _Q5 @ np.diag(eigs) @ _Q5.T
    return (G + G.T) / 2


# ─────────────────────────── the materials ───────────────────────────


def _check_groups(ng: int) -> None:
    if ng not in GROUPS:
        raise ValueError(f"the harness holds the group counts {GROUPS}, not {ng!r}")


def _sn_mixtures(ng: int):
    """The (region 0, region 1) mixtures and the text of their construction.

    ``get_mixture`` knows ``"1g"``, ``"2g"`` and ``"4g"`` only (``orpheus/data/macro_xs/synthetic.py:300-303``), so three groups
    use ``_fuel_3g`` for both regions.
    """
    _check_groups(ng)
    if ng in (1, 2):
        return get_mixture("A", f"{ng}g"), get_mixture("B", f"{ng}g"), f"get_mixture('A'/'B', '{ng}g')"
    m = _fuel_3g()
    return m, m, "_fuel_3g() for both regions"


# ─────────────────────────── the SN problems ───────────────────────────


def _sn_table() -> dict[int, Any]:
    return {0: WallCondition.specular(SPECULAR_ALBEDO), 1: WallCondition.vacuum(), 2: WallCondition.white(WHITE_ALBEDO)}


def _regions(counts):
    return tuple(CellsByCount.uniform_volume(n) for n in counts)


def sn_problem(kind: str, ng: int) -> tuple[SNProblem, str]:
    """One SN problem of the harness and the text of its constructor call (the content of ``snapshots/BUILDERS.txt``)."""
    if kind not in SN_KINDS:
        raise ValueError(f"unknown SN fixture {kind!r}; the harness holds {SN_KINDS}")
    m0, m1, mix_note = _sn_mixtures(ng)
    table = _sn_table()
    if kind == "slab":
        geometry = StructuredGeometry.from_thicknesses(
            coord=CoordSystem.CARTESIAN, thicknesses=SLAB_THICKNESSES, mat_ids=(0, 1, 0), face_roles=(MirrorDeck(), 1))
        mesh = Mesher(geometry).partition(_regions(SLAB_REGION_CELLS)).mesh
        quad = Quadrature.gauss_legendre(n_ordinates=GAUSS_LEGENDRE_ORDER)
        qn = f"Quadrature.gauss_legendre(n_ordinates={GAUSS_LEGENDRE_ORDER})"
        mats = {0: m0, 1: m1}
        note = f"slab, thicknesses {SLAB_THICKNESSES}, faces (MirrorDeck, tag 1 vacuum)"
    elif kind in ("sphere", "cylinder"):
        coord = CoordSystem.SPHERICAL if kind == "sphere" else CoordSystem.CYLINDRICAL
        geometry = StructuredGeometry.from_thicknesses(
            coord=coord, thicknesses=SLAB_THICKNESSES, mat_ids=(0, 1, 0), face_roles=(0,))
        mesh = Mesher(geometry).partition(_regions(SLAB_REGION_CELLS)).mesh
        if kind == "sphere":
            quad = Quadrature.gauss_legendre(n_ordinates=GAUSS_LEGENDRE_ORDER)
            qn = f"Quadrature.gauss_legendre(n_ordinates={GAUSS_LEGENDRE_ORDER})"
        else:
            n_mu, n_phi = CYLINDER_FOLDED
            quad = Quadrature.folded_product(n_mu=n_mu, n_phi=n_phi)
            qn = f"Quadrature.folded_product(n_mu={n_mu}, n_phi={n_phi})"
        mats = {0: m0, 1: m1}
        # A vacuum outer face made the boundary coupling identically zero, so the outer face is specular.
        note = f"{kind}, radii {SLAB_THICKNESSES} by CellsByCount.uniform_volume, outer face tag 0 specular({SPECULAR_ALBEDO})"
    else:
        nx, ny = BOX_CELLS
        material_map = np.zeros((nx, ny), dtype=int)
        material_map[: nx // 2, :] = 2
        mesh = Mesh2D(
            edges_x=np.linspace(0.0, BOX_EXTENTS[0], nx + 1), edges_y=np.linspace(0.0, BOX_EXTENTS[1], ny + 1), mat_map=material_map,
            face_roles={"xmin": MirrorDeck(), "xmax": 1, "ymin": 0, "ymax": 2})
        quad = Quadrature.level_symmetric(sn_order=BOX_LEVEL_SYMMETRIC_ORDER)
        qn = f"Quadrature.level_symmetric(sn_order={BOX_LEVEL_SYMMETRIC_ORDER})"
        mats = {2: m0, 0: m1}
        note = (f"box {nx}x{ny} cells, extents {BOX_EXTENTS}; faces xmin MirrorDeck, xmax tag 1 vacuum, "
                f"ymin tag 0 specular({SPECULAR_ALBEDO}), ymax tag 2 white({WHITE_ALBEDO})")
    problem = SNProblem(mesh, quad, mats, scheme=None, scattering_order=0, boundary_conditions=conditions_for(mesh, table))
    text = (f"{note}; quadrature {qn}; materials {mix_note}; SNProblem(mesh, quad, mats, scheme=None, scattering_order=0, "
            f"boundary_conditions=conditions_for(mesh, {{0: specular({SPECULAR_ALBEDO}), 1: vacuum, 2: white({WHITE_ALBEDO})}}))")
    return problem, text


# ─────────────────────────── the diffusion problem ───────────────────────────


def diffusion_problem(ng: int) -> tuple[DiffusionMesh, DiffusionSolver, str]:
    """The diffusion slab at ``ng`` groups: the mesh, the solver and the text of the constructor call."""
    _check_groups(ng)
    # The body of ``_het_mesh(WallCondition.vacuum(), WallCondition.white(0.3))`` (tests/gates/diffusion/test_solver.py:101),
    # with the cell edges of this harness in place of the gate's 1 / 2 / 1 cells.
    (role_left, role_right), conditions = declare_sides(WallCondition.vacuum(), WallCondition.white(WHITE_ALBEDO))
    geometry = StructuredGeometry.slab((0.0, 0.5, 3.0, 5.0), (0, 1, 0), left=role_left, right=role_right)
    mesh = Mesher(geometry).partition(tuple(CellEdges(np.array(e)) for e in DIFFUSION_CELL_EDGES)).mesh
    if ng == 2:
        mats, note = _het_materials(), "_het_materials() (_fuel_2g, _reflector_2g)"
    elif ng == 1:
        mats, note = {0: get_mixture("A", "1g"), 1: get_mixture("B", "1g")}, "get_mixture('A'/'B', '1g')"
    else:
        m = _fuel_3g()
        mats, note = {0: m, 1: m}, "_fuel_3g() for both regions"
    diffusion_mesh = DiffusionMesh(mesh, mats, boundary_conditions=conditions)
    text = (f"the body of _het_mesh(WallCondition.vacuum(), WallCondition.white({WHITE_ALBEDO})) with cell edges {DIFFUSION_CELL_EDGES}; "
            f"materials {note}; DiffusionSolver(DiffusionMesh(mesh, materials, boundary_conditions=conditions))")
    return diffusion_mesh, DiffusionSolver(diffusion_mesh), text


# ─────────────────────────── the dense-metric Gram families ───────────────────────────


#: The four hand-written families and the nullity each is declared with (the builder's own declaration: an independent witness of the
#: frame's rank rule on those four). The other 17 families are built by a frame and declare through ``frame.frame_gram_kernel``.
HAND_KERNEL_DIMS: dict[str, int] = {"syn.rank2_3x3": 1, "syn.pou_3x3": 0, "syn.kappa1e4_8x8": 0, "syn.near_cutoff_5x5": 1}


def declared_kernel_of(family: str, gram: np.ndarray) -> tuple[Declaration, str]:
    """``(kernel, constructor text)`` the harness declares for the family's raw ``gram``.

    A frame-built family declares through the frame's rule; a hand-written one by its own dimension, spelled inline here (this
    module imports no other test module), keyed ``("test-gram-kernel", content_digest(g_sym))``.
    """
    from orpheus.numerics import frame as FR

    if family not in HAND_KERNEL_DIMS:
        return FR.frame_gram_kernel(gram), "frame_gram_kernel(gram)"
    g_sym = (gram + gram.T) / 2.0
    n, dim = g_sym.shape[0], HAND_KERNEL_DIMS[family]
    if dim == 0:
        return Coord(n, ()), f"Coord({n}, ())"
    frozen = np.array(g_sym, dtype=float)
    frozen.setflags(write=False)
    return Spanned(n, ("test-gram-kernel", content_digest(frozen)), dim), f"Spanned({n}, ('test-gram-kernel', content_digest(g_sym)), {dim})"


def frame_families() -> list[tuple[str, str, Callable[[], np.ndarray]]]:
    """The 21 families as ``(id, constructor call, builder of the raw Gram)``: 17 frame-built and 4 hand-written.

    The 17 frame-built Grams are the ones the 101-file spy of ``scratch/ontology/w1/spec/inv/spy2.pkl.rows`` saw, rebuilt
    explicitly (``syn.equi8.L4_5x5`` is one of them). The synthetic ids start with ``syn.``; four of those are hand-written.
    """
    from orpheus.geometry.manifold import SPHERE, RealSpace
    from orpheus.geometry.symmetry import SubgroupOfO3
    from orpheus.numerics import frame_faces as FF
    from orpheus.numerics.basis import OverlapBasis
    from orpheus.numerics.basis.legendre_basis import LegendreBasis
    from orpheus.numerics.frame import GalerkinFrame
    from orpheus.numerics.measure import DiscreteMeasure

    def ang(q, L):
        return lambda: np.asarray(FF.angular_frame(q(), L).discrete_gram, float)

    def equi(L, n):
        # tests/gates/numerics/test_dense_metric.py:45 (``_equispaced_legendre_gram``): the same linspace, weights 2/n and S^2/O2_x measure.
        def build():
            mu = np.linspace(-1.0, 1.0, n + 2)[1:-1]
            w = np.full(mu.size, 2.0 / mu.size)
            m = DiscreteMeasure(nodes=mu, weights=w, support=SPHERE.quotient(SubgroupOfO3.O2("x")))
            return np.asarray(FF.angular_frame(Quadrature(measure=m), L).discrete_gram, float)
        return build

    def overlap():
        # A construction of the harness (no gate holds it): a two-cell overlap basis on the line and three equal-weight nodes.
        ob = OverlapBasis(edges_per_axis=(np.array([-0.5, 0.5, 1.5]),), partition_of=RealSpace(1),
                          overlap_table=np.array([[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]]), fine=RealSpace(1))
        m = DiscreteMeasure(nodes=np.array([0.0, 0.5, 1.0]), weights=np.ones(3), support=RealSpace(1))
        return np.asarray(GalerkinFrame(ob, m).discrete_gram, float)

    def ls4_legendre():
        # A construction of the harness (no gate holds it): the Legendre basis of degree 4 on the level-symmetric order-4 measure.
        q = Quadrature.level_symmetric(4)
        return np.asarray(GalerkinFrame(LegendreBasis(L=4, axis="x"), q.measure).discrete_gram, float)

    def kappa():
        q, _ = np.linalg.qr(np.random.default_rng(7).standard_normal((8, 8)))
        g = (q * np.geomspace(1.0, 1.0 / 1e4, 8)) @ q.T
        return (g + g.T) / 2.0

    LS, FP, PR = Quadrature.level_symmetric, Quadrature.folded_product, Quadrature.product
    return [
        ("SH.Oh24.L3", "frame_faces.angular_frame(Quadrature.level_symmetric(4), 3).discrete_gram", ang(lambda: LS(4), 3)),
        ("SH.Oh24.L4", "frame_faces.angular_frame(Quadrature.level_symmetric(4), 4).discrete_gram", ang(lambda: LS(4), 4)),
        ("SH.Oh24.L5", "frame_faces.angular_frame(Quadrature.level_symmetric(4), 5).discrete_gram", ang(lambda: LS(4), 5)),
        ("SH.Dnh8.L4", "frame_faces.angular_frame(Quadrature.product(4, 8), 4).discrete_gram", ang(lambda: PR(4, 8), 4)),
        ("SH.Dnh4.L2", "frame_faces.angular_frame(Quadrature.product(4, 4), 2).discrete_gram", ang(lambda: PR(4, 4), 2)),
        ("MirrorSH.fold2x4.L2", "frame_faces.angular_frame(Quadrature.folded_product(2, 4), 2).discrete_gram", ang(lambda: FP(2, 4), 2)),
        ("MirrorSH.fold2x4.L3", "frame_faces.angular_frame(Quadrature.folded_product(2, 4), 3).discrete_gram", ang(lambda: FP(2, 4), 3)),
        ("MirrorSH.fold2x6.L2", "frame_faces.angular_frame(Quadrature.folded_product(2, 6), 2).discrete_gram", ang(lambda: FP(2, 6), 2)),
        ("MirrorSH.fold4x4.L2", "frame_faces.angular_frame(Quadrature.folded_product(4, 4), 2).discrete_gram", ang(lambda: FP(4, 4), 2)),
        ("MirrorSH.fold4x4.L4", "frame_faces.angular_frame(Quadrature.folded_product(4, 4), 4).discrete_gram", ang(lambda: FP(4, 4), 4)),
        ("Legendre.GL2.L2", "frame_faces.angular_frame(Quadrature.gauss_legendre(2), 2).discrete_gram", ang(lambda: Quadrature.gauss_legendre(2), 2)),
        ("Legendre.equi4.L2", "frame_faces.angular_frame(Quadrature(measure=DiscreteMeasure(mu=linspace(-1,1,6)[1:-1], w=2/4, S^2/O2_x)), 2)", equi(2, 4)),
        ("Legendre.equi5.L3", "same, n=5, L=3", equi(3, 5)),
        ("Legendre.equi8.L3", "same, n=8, L=3", equi(3, 8)),
        ("Legendre.LS4.L4", "GalerkinFrame(LegendreBasis(L=4, axis='x'), Quadrature.level_symmetric(4).measure).discrete_gram", ls4_legendre),
        ("Overlap.R1", "GalerkinFrame(OverlapBasis(edges (-.5,.5,1.5), table [[1,0],[.5,.5],[0,1]]), DiscreteMeasure(nodes (0,.5,1), ones, RealSpace(1))).discrete_gram", overlap),
        ("syn.rank2_3x3", "np.array([[1,1,0],[1,1,0],[0,0,2]], float)  (tests/gates/numerics/test_dense_metric.py:234)", lambda: np.array([[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 2.0]])),
        ("syn.pou_3x3", "np.array([[1.25,.25,.25],[.25,1.25,.25],[.25,.25,1.25]])  (a literal of the harness)", lambda: np.array([[1.25, 0.25, 0.25], [0.25, 1.25, 0.25], [0.25, 0.25, 1.25]])),
        ("syn.kappa1e4_8x8", "qr(default_rng(7).standard_normal((8,8))) with geomspace(1, 1e-4, 8), symmetrised  (tests/gates/numerics/test_dense_metric.py:161)", kappa),
        ("syn.equi8.L4_5x5", "frame_faces.angular_frame(Quadrature(equispaced n=8 on S^2/O2_x), 4).discrete_gram  (frame-built, counted synthetic by the memo)", equi(4, 8)),
        ("syn.near_cutoff_5x5", "_gram([1.0, .5, .3, .2, 5e-13])  (tests/gates/numerics/test_quotient_metric.py:96)", lambda: _gram([1.0, 0.5, 0.3, 0.2, 5e-13])),
    ]
