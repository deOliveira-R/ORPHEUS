r"""The SN loss kernel's gauge — the closed-form projector onto :math:`\ker(L+C-S-B)`.

Since the space-membership campaign's U3 (2026-10-31) the projector is the numerics
:class:`~orpheus.numerics.kernel_projection_gauge.KernelProjectionGauge`, whose laws are gated on a
numerics fixture in ``tests/gates/numerics/test_kernel_projection_gauge.py``; SN supplies the pair
generators. Since U0 (2026-10-31) the generators are the diamond-difference half of the hub's predicted kernel
(:func:`~orpheus.sn.operators.loss_kernel_gauge.predicted_discretization_kernel`) and the gauge acts on the FULL
field (bulk then trace, flat): a trace vector ``x`` enters as ``[0_bulk, x]`` (:func:`_on_trace`), its summands
sit at trace indices offset by the bulk size, and every law below reads the trace part. The rows below
are its production instances, and section 10 pins the rebuilt gauge bit-identical to the retired
``LossKernelGauge`` (a carry captured before the carve).

**What has to be true, and how each is checked without using the gauge as its
own oracle.** The construction is a closed form, so a gate that only asked
"does it project onto what it built?" would be vacuous. Every substantive claim
here is pinned against something the construction does not touch:

============================================  ==========================================
claim                                          independent reference
============================================  ==========================================
the dimension is right                         the combinatorial counting law
                                               (:func:`predicted_kernel_dimension`,
                                               which builds no vector), AND a dense
                                               SVD of the assembled :math:`A`
the vectors really are in :math:`\ker A`       the production matvec, per mode, ALL modes
the deviation it removes is the REAL one       a production Gauss-Seidel solve, against
                                               the analytic uniform trace
it is not a universal absorber                 the Jacobi control, which already lands on
                                               the correct member and must be left alone
============================================  ==========================================

**The dense-SVD ground truth, in full.** Four configurations re-run the SVD
live (:func:`test_the_dimension_matches_a_DENSE_SVD_of_the_assembled_operator`);
the construction was originally validated against **13**, and the rest are
recorded here because their fixtures are too large to re-solve every run and
the table is otherwise only in an untracked scratch memo. ``mine`` = the
closed form's column count, ``law`` = the combinatorial count computed by a
separate function that builds no vector — **both of these are** :math:`R`
**alone**; ``T`` = the tangential slots; ``svd`` = dense SVD of the assembled
:math:`A`, which sees the WHOLE kernel, so the table's own identity is
:math:`{\rm svd} = T + R`. ``gap`` =
:math:`\lVert(I - P_{\rm mine})P_{\rm svd}\rVert_2`, meaningful only where
:math:`T = 0` (see below).

====================================  ====================  ====  ====  ====  ====  ===========  =========
case                                  quadrature            mine   law     T   svd      max res        gap
====================================  ====================  ====  ====  ====  ====  ===========  =========
d2 (3,4) refl abs ng=2                level_symmetric(S4)     12    12     0    12     2.80e-16  1.288e-14
d2 (3,4) refl abs ng=2                level_symmetric(S2)      4     4     0     4     2.21e-16  2.220e-14
d2 (5,6) refl c=0.9 ng=3              level_symmetric(S4)     18    18     0    18     4.16e-16  2.049e-14
d2 (3,4) refl abs ng=2                product(4,4)             0     0   224   224  n/a (empty)        n/a
d2 (3,4) refl abs ng=2                product(8,8)            16    16   448   464     9.90e-16        n/a
d2 (3,4) refl abs ng=2                lebedev(11)             18    18   224   242     8.08e-16        n/a
d2 (3,4) GRADED refl abs ng=2         level_symmetric(S4)     12    12     0    12     5.37e-16  1.416e-14
d2 (3,4) x-VAC y-refl ng=2            level_symmetric(S4)      0     0     0     0     0.00e+00        n/a
d3 (2,2,2) refl abs ng=1              level_symmetric(S4)     33    33     0    33     1.50e-16  9.801e-15
d3 (2,2,3) refl abs ng=1              level_symmetric(S4)     39    39     0    39     1.39e-16  1.018e-14
d3 (2,2,2) refl abs ng=2              level_symmetric(S2)     22    22     0    22     1.49e-16  5.790e-15
d3 (2,2,2) refl abs ng=1              product(4,4)             8     8   128   136     4.00e-16        n/a
d3 (2,2,2) xy-refl z-VAC ng=1         level_symmetric(S4)     12    12     0    12     0.00e+00  1.526e-14
====================================  ====================  ====  ====  ====  ====  ===========  =========

Note the structural rows: ``product(4,4)`` at d=2 is **pure T** (224
tangential, ``R = 0``, so the construction correctly returns an EMPTY basis and
``max res`` has no modes to range over), and the mixed-BC row is a negative
control where the correct answer is 0 — an empty basis rather than manufactured
spurious modes.

⛔ **This table was WRONG when it landed at** ``f934ff57``, **and the correction
is instructive.** Every T-bearing row recorded :math:`T + R` under the ``mine``
and ``law`` columns — ``224 / 464 / 242 / 136`` where the true :math:`R` is
``0 / 16 / 18 / 8``. It contradicted its own prose three lines below (which said
``R = 0``) and ``docs/theory/methods/sn/cartesian_multid.rst`` (which was
correct). ``[M]`` re-measured 2026-08-15: **4 of 4** T-bearing rows, not the
3 first reported — the ``d3 (2,2,2) product(4,4)`` row is the one an audit
scanning only the d=2 block misses. The ``gap`` column was wrong for the same
reason and is now ``n/a`` there *by construction*: :math:`P_{\rm mine}` spans
:math:`R` and :math:`P_{\rm svd}` spans :math:`T \oplus R`, so a whole-space
projector distance between them is :math:`\approx 1` and says nothing about the
closed form. Containment — the claim that actually matters — is carried by
``max res`` on every row and by
:func:`test_EVERY_basis_vector_is_annihilated_by_the_production_matvec`.
The gate that would have caught this is
``tests/gates/sn/operators/test_loss_nullspace_reflective_box.py``'s T+R row, which
did not exist until Step 7; every ``_SINGULAR`` fixture here is
``level_symmetric``, where :math:`T = 0` and the error is unspellable.

⚠ **The load-bearing gate is**
:func:`test_each_block_gram_is_the_identity_which_is_what_earns_DIAGONAL`.
``[M]`` the *raw pair generators* — the natural thing to ship — have a Gram that
is ``0.000e+00`` off-diagonal at :math:`d = 2` and **``4.05e-01``** at
``(3,4,5)``. The d=2 reading is vacuous (:math:`\kappa(\{x,y\}) = 1`, so a
:math:`1\times1` Gram is diagonal for free), and a suite that only exercised
:math:`d = 2` would certify ``DIAGONAL`` on a basis that is 43 % off-diagonal
where it matters — producing a silently mis-normalised projector, because
:meth:`~orpheus.numerics.frame_faces.CrossGramInverse.of`'s row-sum probe equals the true
diagonal only when :math:`MR` is diagonal. **Every d=3 row in this file exists
for that reason; do not trim them to speed the suite up.**
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from orpheus.numerics.index_set import IndexSet
from orpheus.data.macro_xs.synthetic import make_mixture
from orpheus.numerics.basis import GramStructure
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.coupled_system import build_within_group_system
from orpheus.sn.operators.loss_kernel_gauge import (
    GaugeFreedom,
    predicted_discretization_kernel,
    predicted_loss_kernel,
    predicted_physical_kernel,
    _anova_dimension,
    _reflection_orbits,
    gauge_freedom,
    predicted_kernel_dimension,
)
from orpheus.sn.solver import SNSolver, _as_problem, _unwindowed_cold_start
from orpheus.mesh import AxisMesh
from orpheus.transport.spatial.linear_discontinuous import LinearDiscontinuous
from orpheus.numerics import frame_faces
from orpheus.geometry.deck import MirrorDeck
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for

#: Exterior conditions of the tags this file's meshes carry (wave 2: a face's role is a deck or an int tag; the condition a tag
#: stands for is declared at the solve, keyword ``boundary_conditions``).
_CONDITIONS = {0: WallCondition.vacuum()}


_R, _V = MirrorDeck(), 0


# ─────────────────────────────────────────────────────────────────────
# fixtures
# ─────────────────────────────────────────────────────────────────────
def _mixture(ng: int = 2, c: float = 0.0, fissile: bool = False):
    sig_t = np.linspace(0.8, 1.6, ng)
    sig_s = np.zeros((ng, ng))
    np.fill_diagonal(sig_s, c * sig_t)
    sig_f = 0.3 * sig_t if fissile else np.zeros(ng)
    chi = np.zeros(ng)
    if fissile:
        chi[0] = 1.0
    return make_mixture(
        sig_t=sig_t, sig_c=sig_t - sig_s.sum(axis=0) - sig_f,
        sig_f=sig_f, nu=np.full(ng, 2.4) if fissile else np.zeros(ng),
        chi=chi, sig_s=sig_s,
    )


def _graded(extent: float, n: int, stretch: float = 1.0) -> np.ndarray:
    """Geometric cell widths — ``stretch == 1.0`` is uniform."""
    widths = stretch ** np.arange(n)
    return np.concatenate([[0.0], np.cumsum(widths / widths.sum() * extent)])


def _mesh(cells, bcs, *, ng=2, quad=None, c=0.0, stretch=1.0,
          fissile=False, scheme=None):
    extents = (1.0, 2.0, 3.0)[: len(cells)]
    axes = tuple(
        AxisMesh(edges=_graded(e, n, stretch), role_low=lo, role_high=hi)
        for e, n, (lo, hi) in zip(extents, cells, bcs)
    )
    return _as_problem(
        axes, quad or Quadrature.level_symmetric(sn_order=4),
        {0: _mixture(ng, c, fissile)}, scheme=scheme, boundary_conditions=conditions_for(axes, _CONDITIONS),
    )


def _loss_system(problem):
    """``(system, template)`` — the PRODUCTION splitting the SI driver iterates."""
    solver = SNSolver(problem, inner_solver="source_iteration")
    system = build_within_group_system(
        problem, solver.problem.mat_xs,
    )
    return system, _unwindowed_cold_start(problem, history_depth=0)


def _apply_loss(system, template, flat: np.ndarray) -> np.ndarray:
    r""":math:`A x = (L + C - S - N_{2n} - B_a)\,x` — the record's factors,
    by role, in the loss-sign order (the same arithmetic the Jacobi value's
    ``M − ΣN`` spells)."""
    f = system.factors
    x = type(template).from_flat(flat, template)
    out = f.streaming_collision.apply(x)
    for gain in (f.scattering, f.n2n, f.boundary):
        out = out - gain.apply(x)
    return out.to_flat()


def _n_bulk(problem) -> int:
    return int(np.prod(problem.full_field_space.interior_space.shape))


def _dd_gauge(problem):
    """The diamond-difference half of the hub's predicted kernel, its projector on the full field (U0). The retired
    ``build_loss_kernel_gauge`` (retired at U0) returned this object on the trace space."""
    return predicted_discretization_kernel(problem).projector


def _on_trace(problem, gauge, verb: str, trace: np.ndarray) -> np.ndarray:
    """A trace vector through a full-field verb: embed as ``[0_bulk, x]``, act, read the trace part back. The bulk part
    of the answer must be exactly zero (the diamond-difference summands live on the trace)."""
    nb = _n_bulk(problem)
    out = np.asarray(getattr(gauge, verb)(np.concatenate([np.zeros(nb), trace])), dtype=float)
    np.testing.assert_array_equal(out[:nb], 0.0)
    return out[nb:]


def _trace_indices(problem, block) -> np.ndarray:
    """A summand's support in the TRACE chart (its gather indices are full-field ones)."""
    return np.asarray(block.gather.indices) - _n_bulk(problem)


def _embed_trace(template, trace: np.ndarray) -> np.ndarray:
    """A trace vector as a full-field flat vector with zero bulk."""
    return np.concatenate(
        [np.zeros(template.interior.values.size), trace])


#: ``(label, cells, bcs, kwargs)``. The d=3 rows are NOT optional — see the
#: module docstring.
_SINGULAR = [
    ("d2 (3,4) LS4 ng=2", (3, 4), [(_R, _R)] * 2, {}),
    ("d2 (5,6) LS4 ng=3 c=0.9", (5, 6), [(_R, _R)] * 2, {"ng": 3, "c": 0.9}),
    ("d2 (3,4) GRADED 1.7x", (3, 4), [(_R, _R)] * 2, {"stretch": 1.7}),
    ("d2 (3,4) LS2 ng=1", (3, 4), [(_R, _R)] * 2,
     {"ng": 1, "quad": Quadrature.level_symmetric(sn_order=2)}),
    ("d3 (2,2,2) LS4 ng=1", (2, 2, 2), [(_R, _R)] * 3, {"ng": 1}),
    ("d3 (2,2,3) LS4 ng=1", (2, 2, 3), [(_R, _R)] * 3, {"ng": 1}),
]

#: Small enough that a dense SVD of ``A`` is viable as a second oracle.
_DENSE_SVD_VIABLE = {
    "d2 (3,4) LS4 ng=2", "d2 (3,4) GRADED 1.7x", "d2 (3,4) LS2 ng=1",
    "d3 (2,2,2) LS4 ng=1",
}


# ─────────────────────────────────────────────────────────────────────
# 1. the dimension, against two independent oracles
# ─────────────────────────────────────────────────────────────────────
@pytest.mark.foundation
@pytest.mark.verifies("dd-null-counting-law")
@pytest.mark.parametrize("label,cells,bcs,kwargs", _SINGULAR,
                         ids=[row[0] for row in _SINGULAR])
def test_the_dimension_matches_the_combinatorial_counting_law(
        label, cells, bcs, kwargs):
    r"""``gauge.dimension`` == the law, computed by a route that builds nothing.

    :func:`predicted_kernel_dimension` walks the ANOVA combinatorics of the
    theory and returns an integer; the gauge walks the pair generators and takes
    an SVD. They share the orbit decomposition and nothing else, so agreement is
    real evidence — and it is the ONLY available oracle at production sizes,
    where the dense SVD is not viable.
    """
    mesh = _mesh(cells, bcs, **kwargs)
    gauge = _dd_gauge(mesh)
    assert gauge.dimension == predicted_kernel_dimension(mesh) > 0


@pytest.mark.foundation
@pytest.mark.verifies("dd-null-counting-law")
@pytest.mark.parametrize(
    "label,cells,bcs,kwargs",
    [row for row in _SINGULAR if row[0] in _DENSE_SVD_VIABLE],
    ids=[row[0] for row in _SINGULAR if row[0] in _DENSE_SVD_VIABLE],
)
def test_the_dimension_matches_a_DENSE_SVD_of_the_assembled_operator(
        label, cells, bcs, kwargs):
    """The numerical oracle: rank-deficiency of ``A`` itself.

    ``level_symmetric`` places no tangential ordinate, so ``T`` is empty and the
    whole kernel is the component ``R`` the gauge builds — which is what makes
    this an equality rather than an inequality. The ``product`` case, where
    ``T`` is large, is covered separately.
    """
    mesh = _mesh(cells, bcs, **kwargs)
    system, template = _loss_system(mesh)
    n_dof = template.to_flat().size

    dense = np.empty((n_dof, n_dof))
    unit = np.zeros(n_dof)
    for column in range(n_dof):
        unit[column] = 1.0
        dense[:, column] = _apply_loss(system, template, unit)
        unit[column] = 0.0
    singular = np.linalg.svd(dense, compute_uv=False)
    rank_deficiency = int(np.sum(singular < 1e-10 * singular[0]))

    assert _dd_gauge(mesh).dimension == rank_deficiency


@pytest.mark.foundation
def test_the_counting_law_reproduces_the_two_closed_form_specialisations():
    r"""``ng * N/4`` at d=2 (mesh-INDEPENDENT) and ``ng*(N/8)*(2*sum(n)-1)`` at d=3.

    Pinned as literals derived from the theory, so a change to the orbit
    decomposition cannot quietly move both the law and the construction
    together.
    """
    n_ordinates = Quadrature.level_symmetric(sn_order=4).weights.size
    for cells in [(3, 4), (5, 6), (2, 2)]:          # d=2 is mesh-independent
        mesh = _mesh(cells, [(_R, _R)] * 2, ng=2)
        assert predicted_kernel_dimension(mesh) == 2 * (n_ordinates // 4)

    for cells in [(2, 2, 2), (2, 2, 3), (3, 4, 5)]:
        mesh = _mesh(cells, [(_R, _R)] * 3, ng=1)
        assert predicted_kernel_dimension(mesh) == (
            (n_ordinates // 8) * (2 * sum(cells) - 1)
        )


@pytest.mark.foundation
def test_the_anova_dimension_is_the_separable_equation_solution_count():
    r""":math:`\kappa(\{a,b\}) = 1` and :math:`\kappa(\{a,b,c\}) = n_a+n_b+n_c-1`.

    The two specialisations the theory names, checked against the general
    formula — the two-term case says "a function of :math:`i_y` plus a function
    of :math:`i_x` vanishes ⟹ both are constant", which is why the d=2 kernel
    carries no mesh freedom at all.
    """
    for n_a, n_b in [(2, 2), (3, 4), (5, 7), (1, 9)]:
        assert _anova_dimension((n_a, n_b)) == 1
    for cells in [(2, 2, 2), (3, 4, 5), (2, 3, 7)]:
        assert _anova_dimension(cells) == sum(cells) - 1


# ─────────────────────────────────────────────────────────────────────
# 2. the vectors are in ker A — every one of them, via the production matvec
# ─────────────────────────────────────────────────────────────────────
@pytest.mark.foundation
@pytest.mark.verifies("dd-null-sawtooth")
@pytest.mark.parametrize("label,cells,bcs,kwargs", _SINGULAR,
                         ids=[row[0] for row in _SINGULAR])
def test_EVERY_basis_vector_is_annihilated_by_the_production_matvec(
        label, cells, bcs, kwargs):
    r"""``max ||A phi|| / ||phi||`` over ALL modes — not a sample (`vv` #13).

    What each fixture is here to catch, **measured by mutation** rather than
    asserted (an unwitnessed fixture is `plan-authoring` §6c):

    * **GRADED 1.7x** pins that the face area is a *function of the transverse
      index*, :math:`A_a(i_\perp) = \prod_{b \neq a} h_b(i_b)`. ``[M]``
      replacing :math:`h_b(i_b)` with :math:`h_b(0)` — so the area becomes a
      constant — reddens **exactly these two graded rows and nothing else**.
      ⚠ It does NOT pin the area's *presence*: deleting the factor outright
      reddens all six rows, because the extents ``(1.0, 2.0, 3.0)`` carry
      unequal cell counts, so :math:`h_x \neq h_y` on a uniform mesh too.
    * **LS2 ng=1 / LS4 ng=2 / ng=3 c=0.9** cover the ordinate count, the group
      count and a scattering ratio, all of which cancel out of the derivation
      and must therefore change nothing.
    * **d3 (2,2,3)** carries an ODD cell count, and that is load-bearing:
      ``[M]`` a mutation that drops the :math:`(-1)^{n_a}` far-face sign
      reddens ``(2,2,3)`` but **not** ``(2,2,2)``, where every :math:`n_a` is
      even and the dropped factor is :math:`+1` anyway (`vv` #13's
      congruence-class trap, in this suite's own fixture set).
    """
    mesh = _mesh(cells, bcs, **kwargs)
    gauge = _dd_gauge(mesh)
    system, template = _loss_system(mesh)
    n_trace = int(np.prod(mesh.angular_trace.shape))

    worst = 0.0
    modes = 0
    for block in gauge.summands:
        for mode in range(block.basis.table.shape[1]):
            trace = np.zeros(n_trace)
            trace[_trace_indices(mesh, block)] = block.basis.table[:, mode]
            full = _embed_trace(template, trace)
            residual = np.linalg.norm(_apply_loss(system, template, full))
            worst = max(worst, residual / np.linalg.norm(full))
            modes += 1
    assert modes == gauge.dimension
    assert worst < 1e-13, f"{label}: worst residual {worst:.3e} over {modes} modes"


# ─────────────────────────────────────────────────────────────────────
# 3. ⭐ the orthonormality that EARNS the DIAGONAL declaration
# ─────────────────────────────────────────────────────────────────────
@pytest.mark.foundation
@pytest.mark.parametrize("label,cells,bcs,kwargs", _SINGULAR,
                         ids=[row[0] for row in _SINGULAR])
def test_each_block_gram_is_the_identity_which_is_what_earns_DIAGONAL(
        label, cells, bcs, kwargs):
    r"""⭐ :math:`\Phi^{\mathsf T} G \Phi = I` per block — the load-bearing gate.

    :attr:`~orpheus.numerics.basis.orthonormal_span_basis.OrthonormalSpanBasis.gram_structure` declares ``DIAGONAL``, and
    :meth:`~orpheus.numerics.frame_faces.CrossGramInverse.of` computes its diagonal by the
    row-sum probe ``analysis(reconstruction(ones))``, which equals the true
    diagonal ONLY if :math:`MR` is diagonal. So the declaration is a promise
    about this matrix, and if it were false the projector would be silently
    mis-normalised rather than loudly wrong.

    ⚠ ``[M]`` the raw pair generators — what the derivation produces before the
    :math:`\sqrt G`-weighted SVD — read ``0.000e+00`` off-diagonal at d=2 and
    **``4.05e-01``** at ``(3,4,5)``. The d=2 reading is *vacuous*
    (:math:`\kappa(\{x,y\}) = 1`: one mode per orbit, so the Gram is
    :math:`1\times1`). **The d=3 rows are the only ones that can fail this.**
    """
    mesh = _mesh(cells, bcs, **kwargs)
    gauge = _dd_gauge(mesh)
    metric = np.asarray(mesh.angular_trace.inner_product_weights, dtype=float)

    saw_multimode = False
    for block in gauge.summands:
        table = block.basis.table
        assert block.basis.gram_structure is GramStructure.DIAGONAL
        gram = table.T @ (metric[_trace_indices(mesh, block)][:, None] * table)
        saw_multimode |= table.shape[1] > 1
        assert np.allclose(gram, np.eye(table.shape[1]), atol=1e-12), (
            f"{label}: summand {block.basis.label} Gram is "
            f"off the identity by {np.abs(gram - np.eye(len(gram))).max():.3e}"
        )
    if len(cells) >= 3:
        assert saw_multimode, (
            "a d=3 fixture must carry a block with >1 mode, or this gate is "
            "the vacuous d=2 reading in disguise"
        )


@pytest.mark.foundation
def test_the_frame_gram_probe_agrees_with_the_true_gram():
    """The frame's own row-sum shortcut returns all ones, as orthonormality implies.

    Closes the loop between the declaration and the machinery that consumes it:
    this asserts the value ``CrossGramInverse.of(frame)`` actually computes, not the matrix
    we believe it stands for.
    """
    from orpheus.numerics.frame import GalerkinFrame
    from orpheus.numerics.measure import DiscreteMeasure

    mesh = _mesh((2, 2, 2), [(_R, _R)] * 3, ng=1)
    gauge = _dd_gauge(mesh)
    assert gauge.summands, "fixture is no longer singular"
    metric = np.asarray(mesh.angular_trace.inner_product_weights, dtype=float)
    for block in gauge.summands:
        indices = _trace_indices(mesh, block)
        # ⛔ RE-KEYED 2026-09-02 (#429 tracker 2.2). The probe measure named
        # its own manifold (``IndexSet(label="probe")``) while the block's
        # basis eats ``index(sn_trace_orbit(...)_g0)`` — two different point
        # sets, which the frame's new G0 arrow refuses. The probe was never
        # about the LABEL: it integrates the block's own trace metric over the
        # block's own indices, so it must name the basis's manifold. Naming
        # the same object on both halves is what production does
        # (``frame.py:818`` binds one ``points``).
        frame = GalerkinFrame(
            block.basis,
            DiscreteMeasure(nodes=indices.astype(float),
                            weights=metric[indices], support=block.basis.domain),
        )
        diagonal = np.asarray(frame_faces.CrossGramInverse.of(frame).diagonal)
        np.testing.assert_allclose(diagonal, np.ones_like(diagonal), atol=1e-12)


# ─────────────────────────────────────────────────────────────────────
# 4. it is a projector, and a G-self-adjoint one
# ─────────────────────────────────────────────────────────────────────
@pytest.mark.foundation
@pytest.mark.parametrize("label,cells,bcs,kwargs", _SINGULAR,
                         ids=[row[0] for row in _SINGULAR])
def test_it_is_an_idempotent_G_self_adjoint_projector(
        label, cells, bcs, kwargs):
    r""":math:`\Pi^2 = \Pi` and :math:`\langle \Pi x, y\rangle_G =
    \langle x, \Pi y\rangle_G`.

    The defining laws of a :math:`G`-orthogonal projector, asserted as laws
    rather than inferred from the construction.
    """
    mesh = _mesh(cells, bcs, **kwargs)
    gauge = _dd_gauge(mesh)
    metric = np.asarray(mesh.angular_trace.inner_product_weights, dtype=float)
    rng = np.random.default_rng(0)
    x = rng.standard_normal(metric.size)
    y = rng.standard_normal(metric.size)

    def apply(v):
        return _on_trace(mesh, gauge, "apply", v)

    projected = apply(x)
    assert np.linalg.norm(projected) > 1e-6, "degenerate fixture: nothing projected"
    np.testing.assert_allclose(apply(projected), projected, atol=1e-12)
    assert float(np.sum(metric * apply(x) * y)) == pytest.approx(
        float(np.sum(metric * x * apply(y))), rel=1e-11
    )


@pytest.mark.foundation
@pytest.mark.rests_on(
    "tests/gates/numerics/test_kernel_projection_gauge.py::TestTheDirectSum",
)
def test_the_summands_have_DISJOINT_supports():
    """SN's builder hands disjoint supports (one per ordinate orbit and group), so the direct sum needs no
    orthogonality check here.

    Until U3 this row also asserted that ``LossKernelGauge`` REFUSED any overlap (fragment "DISJOINT"). That contract
    is inverted (the space-membership campaign's R3: overlapping summands are admitted iff mutually G-orthogonal), it
    is numerics', and it is gated there with both legs
    (``TestTheDirectSum::test_overlapping_summands_that_are_NOT_G_orthogonal_are_refused``). The disjointness that
    remains is a fact of SN's builder.
    """
    mesh = _mesh((2, 2, 2), [(_R, _R)] * 3, ng=2)
    gauge = _dd_gauge(mesh)
    assert len(gauge.summands) > 1
    seen = np.concatenate([b.gather.indices for b in gauge.summands])
    assert np.unique(seen).size == seen.size


# ─────────────────────────────────────────────────────────────────────
# 5. residual-neutrality — no convergence check may move
# ─────────────────────────────────────────────────────────────────────
@pytest.mark.foundation
@pytest.mark.verifies("sn-loss-kernel-gauge-projection")
@pytest.mark.parametrize("label,cells,bcs,kwargs", _SINGULAR,
                         ids=[row[0] for row in _SINGULAR])
def test_gauging_cannot_move_any_convergence_check(
        label, cells, bcs, kwargs):
    r""":math:`A(\psi - \Pi\psi) = A\psi`, because :math:`\Pi\psi \in \ker A`.

    Asserted, not assumed. This is what makes the gauge safe to fire at a
    converged exit: the residual, the balance projection and every other
    functional of :math:`A\psi` are bit-unaffected.
    """
    mesh = _mesh(cells, bcs, **kwargs)
    gauge = _dd_gauge(mesh)
    system, template = _loss_system(mesh)
    rng = np.random.default_rng(1)
    trace = rng.standard_normal(int(np.prod(mesh.angular_trace.shape)))
    bulk = rng.standard_normal(template.interior.values.size)

    raw = np.concatenate([bulk, trace])
    gauged = np.concatenate([bulk, _on_trace(mesh, gauge, "gauge", trace)])
    before = _apply_loss(system, template, raw)
    after = _apply_loss(system, template, gauged)
    assert np.linalg.norm(after - before) / np.linalg.norm(before) < 1e-13


# ─────────────────────────────────────────────────────────────────────
# 6. ⭐ THE FLAGSHIP — and its negative control
# ─────────────────────────────────────────────────────────────────────
def _uniform_source_fixture(cells, ng=2):
    r"""All-reflective box + uniform isotropic source, whose exact answer is flat.

    :math:`\psi = Q / (\Sigma_t \sum w)` everywhere, bulk and trace — an
    analytic reference that owes nothing to the gauge.
    """
    from orpheus.sn.solver import _build_fixed_source_rhs

    mesh = _mesh(cells, [(_R, _R)] * len(cells), ng=ng)
    system, template = _loss_system(mesh)
    total_weight = float(mesh.quad.weights.sum())
    sig_t = np.asarray(mesh.materials[0].SigT, dtype=float)
    per_group = 1.0 / (total_weight * sig_t)

    interior = np.zeros(template.interior.values.shape)
    interior[:] = 1.0 / total_weight
    source = _build_fixed_source_rhs(interior, mesh)

    n_bulk = template.interior.values.size
    exact = np.zeros(template.to_flat().size)
    block = np.zeros(template.interior.values.shape)
    for g in range(ng):
        block[:, g, ...] = per_group[g]
    exact[:n_bulk] = block.ravel()
    for slot in mesh.angular_trace.layout.faces.values():
        face = np.zeros(slot.shape)
        for g in range(ng):
            face[:, g, ...] = per_group[g]
        exact[n_bulk + slot.offset:
              n_bulk + slot.offset + slot.flat_size] = face.ravel()
    return mesh, system, template, source, exact, n_bulk


def _solve(system, mesh, template, source, schedule):
    from orpheus.numerics.iteration import SourceIteration
    from orpheus.sn.splitting import Splitting, resolve_schedule

    # The Strategy VALUE for ``schedule`` (step 2 of the consumers campaign):
    # its derived ``M`` and lagged pieces, exactly as the driver runs them.
    value = Splitting.from_schedule(system, resolve_schedule(mesh, schedule))
    iteration = SourceIteration(value.implicit.inverse(), *value.explicit,
                                max_iter=400_000, tol=1e-13)
    zero = type(template).from_flat(
        np.zeros(template.to_flat().size), template)
    solution, record = iteration.solve(source, initial_guess=zero)
    assert record.converged, f"{schedule} did not converge"
    return solution.to_flat()


@pytest.mark.foundation
@pytest.mark.parametrize("cells", [(3, 4), (3, 3)])
def test_the_gauge_recovers_the_PHYSICAL_trace_from_a_production_solve(cells):
    r"""⭐ The flagship: Gauss-Seidel's 6 % trace error collapses to round-off.

    The deviation is produced by the **production** boundary Gauss-Seidel
    splitting on a fixture whose exact answer is analytic, so neither side of
    the comparison comes from this module. ``[M]`` at ``(3,4)``:
    ``6.09e-02 -> 1.19e-13``, with **100.0000 %** of the deviation inside the
    closed-form span.
    """
    mesh, system, template, source, exact, n_bulk = _uniform_source_fixture(cells)
    gauge = _dd_gauge(mesh)
    psi = _solve(system, mesh, template, source, "gauss_seidel")

    reference = np.linalg.norm(exact[n_bulk:])
    deviation = psi[n_bulk:] - exact[n_bulk:]
    before = np.linalg.norm(deviation) / reference
    after = np.linalg.norm(_on_trace(mesh, gauge, "gauge", psi[n_bulk:]) - exact[n_bulk:]) / reference

    assert before > 1e-3, (
        f"fixture no longer exhibits the defect (deviation {before:.3e}) — if "
        f"the splitting was fixed, this gate has nothing to catch"
    )
    assert after < 1e-11, f"{before:.3e} -> {after:.3e}"
    # The whole deviation is kernel content, so the projector is the right tool.
    in_span = (np.linalg.norm(_on_trace(mesh, gauge, "apply", deviation))
               / np.linalg.norm(deviation))
    assert in_span == pytest.approx(1.0, abs=1e-6)


@pytest.mark.foundation
def test_it_is_NOT_a_universal_absorber_the_jacobi_control():
    r"""⭐ The negative control (`vv` #19) — Jacobi is already right, and stays right.

    Jacobi lands on the correct member of the solution manifold, so its residual
    deviation is pure round-off and lies essentially OUTSIDE ``ker A``. A
    projector that "improved" it would be absorbing whatever it is handed.
    ``[M]`` ``1.7741e-13 -> 1.7741e-13``, with only 0.3 % of the round-off
    deviation in span.

    Without this leg, the flagship above is satisfied by ``gauge = identity map
    onto the exact answer``, which is not a thing but is exactly the shape of
    error a green-only suite cannot see.
    """
    mesh, system, template, source, exact, n_bulk = _uniform_source_fixture((3, 4))
    gauge = _dd_gauge(mesh)
    psi = _solve(system, mesh, template, source, "jacobi")

    reference = np.linalg.norm(exact[n_bulk:])
    deviation = psi[n_bulk:] - exact[n_bulk:]
    before = np.linalg.norm(deviation) / reference
    assert before < 1e-11, "jacobi is supposed to already be correct here"

    after = np.linalg.norm(_on_trace(mesh, gauge, "gauge", psi[n_bulk:]) - exact[n_bulk:]) / reference
    assert after == pytest.approx(before, rel=1e-6), (
        "the gauge moved a solution that was already the canonical member"
    )
    in_span = (np.linalg.norm(_on_trace(mesh, gauge, "apply", deviation))
               / np.linalg.norm(deviation))
    assert in_span < 0.05, (
        f"{in_span:.4f} of pure round-off reads as kernel content — the "
        f"projector is absorbing more than it should"
    )


# ─────────────────────────────────────────────────────────────────────
# 7. the predicate — every way of NOT having a gauge
# ─────────────────────────────────────────────────────────────────────
@pytest.mark.foundation
def test_one_vacuum_axis_removes_the_freedom():
    """`[M]` (#344) a single vacuum face collapses ``dim ker A`` from 12 to 0.

    An undamped face mode with only one closed axis pair has nowhere to return
    from — the loop does not close.
    """
    mesh = _mesh((3, 4), [(_V, _V), (_R, _R)])
    verdict = gauge_freedom(mesh)
    assert not verdict.present and not verdict.undetermined
    assert "1 reflective axis pair" in verdict.because
    assert _dd_gauge(mesh).dimension == 0
    assert predicted_kernel_dimension(mesh) == 0


@pytest.mark.foundation
def test_a_damping_closure_removes_the_freedom_at_the_ROOT():
    r"""Linear-discontinuous on the IDENTICAL box has no kernel.

    This is the remedy the warning is obliged to name: the freedom is a property
    of the *closure*, not of the geometry, so changing the discretization
    removes it at the root rather than projecting it out afterwards. ``[M]``
    (#344) ``dim ker A == 0`` for LD on the same all-reflective box where
    diamond gives 12.
    """
    mesh = _mesh((3, 4), [(_R, _R)] * 2, scheme=LinearDiscontinuous())
    verdict = gauge_freedom(mesh)
    assert not verdict.present and not verdict.undetermined
    assert "DAMPS" in verdict.because
    assert mesh.reflective_axis_pairs == 2, (
        "the geometry conjunct must still be satisfied, or this row proves "
        "nothing about the CLOSURE conjunct"
    )
    assert _dd_gauge(mesh).dimension == 0


@pytest.mark.foundation
def test_one_dimension_has_no_face_mode_and_therefore_no_gauge():
    """``d = 1`` falls out of the spectral radius — no special case anywhere."""
    mesh = _mesh((5,), [(_R, _R)])
    assert not gauge_freedom(mesh).present
    assert _dd_gauge(mesh).dimension == 0


@pytest.mark.foundation
def test_an_UNDETERMINED_closure_does_NOT_gauge_and_says_why():
    """The third state: unclassified is not "no freedom".

    ``[M]`` LD at ``ndim = 3`` cannot be driven (its ``assemble_inflow_axis``
    handles ``axis in {0, d-1}`` only). The ruled behaviour is to warn loudly
    and NOT gauge — silently treating it as DAMPED would skip the gauge on a
    scheme that was never examined.
    """
    mesh = _mesh((2, 2, 2), [(_R, _R)] * 3, ng=1, scheme=LinearDiscontinuous())
    verdict = gauge_freedom(mesh)
    assert verdict.undetermined and not verdict.present
    assert "could not be classified" in verdict.because
    assert _dd_gauge(mesh).dimension == 0


@pytest.mark.foundation
def test_the_verdict_refuses_to_be_both_present_and_undetermined():
    """An unclassified closure yields no verdict about its face mode."""
    with pytest.raises(ValueError, match="cannot be both"):
        GaugeFreedom(present=True, undetermined=True, because="impossible")


@pytest.mark.foundation
def test_a_quadrature_not_closed_under_a_mirror_is_REFUSED():
    """A reflective BC needs the mirror partner to exist in the rule.

    Not a limitation of this construction — the reflected ordinate has nowhere
    to land, so the boundary condition itself is unrealisable. Refusing here is
    cheaper than producing modes that are not in ``ker A``.
    """
    from types import SimpleNamespace

    mesh = _mesh((3, 4), [(_R, _R)] * 2)
    mu_x = np.asarray(mesh.quad.mu_x, dtype=float).copy()
    mu_x[0] += 0.25                     # its mirror partner no longer exists

    # A stub carrying exactly what `_reflection_orbits` reads — a real mesh
    # cannot hold this state, which is the point: the refusal must fire before
    # anything downstream trusts the orbit decomposition.
    #
    # ⚠ `reflective_axes` is read from the MESH, not recomputed here, and that
    # is deliberate: the criterion is single-sourced on
    # `SNProblem.reflective_axes` (2026-08-15 — it used to be a twin inside this
    # module). A surrogate must honour the contract it stands in for, so it
    # forwards the real mesh's answer rather than re-deriving one; the thing
    # this stub perturbs is the QUADRATURE, and nothing else.
    broken = SimpleNamespace(
        quad=SimpleNamespace(
            mu_x=mu_x,
            mu_y=np.asarray(mesh.quad.mu_y, dtype=float),
            mu_z=np.asarray(mesh.quad.mu_z, dtype=float),
        ),
        face_labels=mesh.face_labels,
        bc=mesh.bc,
        reflective_axes=mesh.reflective_axes,
    )
    with pytest.raises(ValueError, match="not closed under the axis-0 mirror"):
        _reflection_orbits(broken)  # type: ignore[arg-type]  # deliberate stub

    # …and the unmodified rule is accepted, so the guard is not always-on.
    assert _reflection_orbits(mesh)


# ─────────────────────────────────────────────────────────────────────
# 8. the kernel is geometry-only, and the mesh caches it
# ─────────────────────────────────────────────────────────────────────
@pytest.mark.foundation
def test_the_basis_never_reads_a_CROSS_SECTION():
    """`[M]` (#344) an absorber and a fissile mixture give the identical basis.

    The Stratum-1 claim that licenses caching on the mesh and reusing one build
    across every group, outer and eigenvalue iterate. Asserted bit-exactly:
    every factor in the mode is a mesh, quadrature or boundary quantity, so
    anything less than bit-identity means a cross-section leaked in.
    """
    absorber = _dd_gauge(_mesh((3, 4), [(_R, _R)] * 2))
    fissile = _dd_gauge(
        _mesh((3, 4), [(_R, _R)] * 2, fissile=True, c=0.4))

    assert absorber.dimension == fissile.dimension > 0
    for left, right in zip(absorber.summands, fissile.summands):
        np.testing.assert_array_equal(left.gather.indices, right.gather.indices)
        np.testing.assert_array_equal(left.basis.table, right.basis.table)


@pytest.mark.foundation
def test_the_mesh_caches_the_gauge():
    """One build per mesh — the setup cost must not ride the per-solve path."""
    mesh = _mesh((3, 4), [(_R, _R)] * 2)
    assert mesh.loss_kernel_gauge is mesh.loss_kernel_gauge
    assert mesh.loss_kernel_gauge.dimension == predicted_kernel_dimension(mesh)


# ─────────────────────────────────────────────────────────────────────
# 9. component T — out of scope, and LEFT ALONE rather than mangled
# ─────────────────────────────────────────────────────────────────────
@pytest.mark.foundation
def test_the_TANGENTIAL_component_is_untouched_not_annihilated():
    r"""``T`` lies in :math:`\ker G`, so :math:`(I - \Pi)t = t` — exactly.

    ``product(4,4)`` puts **half** the trace rows on tangential ordinates, where
    the metric is exactly zero and no minimum-:math:`G`-norm representative
    exists. The right behaviour is to leave them alone, and it is what falls out:
    :math:`Gt = 0` makes :math:`t` :math:`G`-orthogonal to every mode, so the
    projection sees nothing.

    ⚠ A suite using only ``level_symmetric`` (zero tangential rows) is blind to
    this whole class.
    """
    mesh = _mesh((3, 4), [(_R, _R)] * 2,
                 quad=Quadrature.product(n_mu=4, n_phi=4))
    gauge = _dd_gauge(mesh)
    metric = np.asarray(mesh.angular_trace.inner_product_weights, dtype=float)

    tangential = metric == 0.0
    assert tangential.sum() > 0, "fixture has no tangential rows to test"

    field = np.zeros(metric.size)
    field[tangential] = np.random.default_rng(2).standard_normal(
        int(tangential.sum()))
    np.testing.assert_array_equal(_on_trace(mesh, gauge, "gauge", field), field)

    # And no block ever claims a zero-metric DOF in the first place.
    for block in gauge.summands:
        assert np.all(metric[_trace_indices(mesh, block)] > 0.0)


@pytest.mark.foundation
@pytest.mark.parametrize("label, cells, bcs, kwargs", _SINGULAR,
                         ids=[row[0] for row in _SINGULAR])
def test_each_blocks_frame_names_ONE_manifold_on_both_halves(
        label, cells, bcs, kwargs):
    r"""⭐ The production gauge's basis and measure agree on the point set.

    Every block is a :class:`~orpheus.numerics.frame.GalerkinFrame` over a
    :class:`~orpheus.numerics.basis.orthonormal_span_basis.OrthonormalSpanBasis` and a
    :class:`~orpheus.numerics.measure.DiscreteMeasure` built four lines apart
    in the summand mint (``KernelProjectionGauge.direct_sum`` since U3). Until 2026-09-01 the two named the same point
    set in **two spellings** — the measure tagged the bare label
    ``sn_trace_orbit(…)_g…`` while the basis wrapped it as ``index(…)`` — and
    tracker 2.1 pinned that divergence in
    ``tests/gates/numerics/test_basis_domain.py::test_d6``, explicitly so that
    tracker 2.0c *"must come back here and cannot resolve it by accident"*.

    ⛔ **It was resolved by accident, and this gate is why that is now
    visible.** 2.0c made the measure read ``support=basis.domain``, which
    closes the divergence by construction — and a mutation battery then
    measured the repair as **BLIND**: replacing that expression with a wrong
    ``IndexSet`` reddened **nothing tree-wide**, because ``test_d6`` asserts
    the *basis*'s half and no gate reached the production *measure*'s.

    The observable path is the block's own operator wiring: the gather's
    codomain IS the frame's measure space, one axis whose support is the
    measure's and whose label is the index set's name (S7, ruling 63). So
    this compares the two halves through production plumbing rather than by
    re-reading the constructor. `[M]` the mutation now reddens every row.

    (``vv-principles`` #17: the verdict of a battery is a TABLE, and an arm
    that reddens nothing is a claim with no witness — not a claim that holds.)
    """
    gauge = _dd_gauge(_mesh(cells, bcs, **kwargs))
    assert gauge.summands, f"{label}: no summands — the fixture is not singular"
    for block in gauge.summands:
        domain = block.basis.domain
        assert isinstance(domain, IndexSet)
        assert block.gather.codomain is not None
        axes = block.gather.codomain.axes
        assert axes is not None and axes[0].support == domain and axes[0].label == domain.name, (
            f"{label}: summand {block.basis.label} — the "
            f"basis says its functions live on {domain.name!r} while its "
            f"measure induced {axes[0].support!r}. One frame, two "
            f"manifolds."
        )


@pytest.mark.parametrize("label,cells,bcs,kwargs", _SINGULAR, ids=[row[0] for row in _SINGULAR])
def test_the_kernel_modes_are_orthogonal_to_the_constants(label, cells, bcs, kwargs) -> None:
    r"""The WITNESS of the declared zero unit (S7; the user's ruling, 2026-10-21).

    The SN builder declares (``orthogonal_to_constants=True``) its
    coefficient axis's unit zero: the kernel modes are orthogonal to the
    constants, :math:`\Phi^{\mathsf T} G\,1 = 0`. A computed unit would be
    round-off (`[M]` S7 reviews: at most 1.8e-16 of its scale on 102 blocks,
    exactly 0.0 on only 12), and the unit is identity content, so the zero is
    declared rather than computed and this row holds the claim to the shipped
    blocks: on every block of every singular fixture,
    :math:`|\Phi^{\mathsf T} G 1| \le 64\,\varepsilon\,\|\Phi\|\,\|G 1\|`.
    The general proof is owed: #636.

    Reddens under: a gauge whose modes are not orthogonal to the constants
    (e.g. a block table with a constant column appended), against which the
    declared zero would be false.
    """
    gauge = _dd_gauge(_mesh(cells, bcs, **kwargs))
    assert gauge.summands, f"{label}: no summands — the fixture is not singular"
    eps = np.finfo(float).eps
    for block in gauge.summands:
        codomain = block.gather.codomain
        assert codomain is not None
        axes = codomain.axes
        assert axes is not None
        table = block.basis.table
        g = np.ones(table.shape[0]) if axes[0].weights is None else np.asarray(axes[0].weights)
        moments = table.T @ g
        bound = 64.0 * eps * np.linalg.norm(table) * np.linalg.norm(g)
        assert np.all(np.abs(moments) <= bound), (
            f"{label}: summand {block.basis.label}: |Φᵀ G 1| = "
            f"{np.max(np.abs(moments)):.3e} exceeds {bound:.3e} — the declared zero unit is false here"
        )
        unit = block.basis.layout.unit
        assert not np.any(unit), "the declared unit is exactly zero"


# ─────────────────────────────────────────────────────────────────────
# 10. the U3 re-base: bit-identical to the retired class; the exit reads the Protocol, in the G-norm
# ─────────────────────────────────────────────────────────────────────
#: Captured by ``LossKernelGauge.for_mesh`` (retired at U3) at ``1b2a6a14``, the six ``_SINGULAR`` fixtures, one seeded
#: trace each: ``apply``, ``gauge``, ``dimension``. `[M]` array_equal on 24 of 24 keys to the ``e42fd498`` capture
#: (``scratch/membership_covector/u3/baseline_currency.py``, with a one-ULP positive control).
_CARRY = Path(__file__).with_name("gauge_carry_1b2a6a14.npz")


#: The carry's re-baselines, as DELTAS with a stated bound (the .npz is never rewritten).
_REBASELINES = Path(__file__).with_name("gauge_carry_rebaselines.json")


@pytest.mark.foundation
@pytest.mark.rests_on(
    "tests/gates/numerics/test_kernel_projection_gauge.py::TestTheProjectorLaws",
    "tests/gates/numerics/test_kernel_projection_gauge.py::TestTheDirectSum",
    "tests/gates/sn/operators/test_loss_kernel_gauge.py::test_the_dimension_matches_the_combinatorial_counting_law",
)
@pytest.mark.parametrize("label,cells,bcs,kwargs", _SINGULAR, ids=[row[0] for row in _SINGULAR])
def test_the_rebuilt_gauge_agrees_with_the_retired_class_within_the_declared_rebaseline(label, cells, bcs, kwargs):
    """The U3 carve moved the sqrt(G)-SVD verbatim and was bit-identical to the retired ``LossKernelGauge``. The U3
    review then spelled Pi in the operator algebra (the user's ruling (c), 2026-10-31), which reorders the summation:
    accepted under principled equivalence (vv-principles, "Bit-identity vs principled-equivalence"), the bound and its
    reason read from ``gauge_carry_rebaselines.json`` (`[M]` at most 0.55 eps of max|Pi x| on ``apply``, 0.52 eps of
    max|x| on ``gauge``, LS2 still bit-equal; the bound is 2 eps, two roundings per element). The dimension is exact.
    Reddens under every numerics arm that moves a value (``leak``, ``euclid``, ``gauge_sign``, ``drop_summand``)."""
    import json
    from orpheus.numerics.kernel_projection_gauge import KernelProjectionGauge
    carry = np.load(_CARRY)
    bound = json.loads(_REBASELINES.read_text())["events"][-1]["bound"]
    eps = np.finfo(float).eps
    mesh = _mesh(cells, bcs, **kwargs)
    gauge = mesh.loss_kernel_gauge   # since U0 the full-field projector: the carry reads its trace part on [0_bulk, x]
    assert isinstance(gauge, KernelProjectionGauge)
    x = carry[f"{label}|x"]
    applied, reference = _on_trace(mesh, gauge, "apply", x), carry[f"{label}|apply"]
    np.testing.assert_allclose(applied, reference, rtol=0,
                               atol=bound["apply_atol_in_eps_max_abs_apply"] * eps * np.abs(reference).max())
    np.testing.assert_allclose(_on_trace(mesh, gauge, "gauge", x), carry[f"{label}|gauge"], rtol=0,
                               atol=bound["gauge_atol_in_eps_max_abs_x"] * eps * np.abs(x).max())
    assert gauge.dimension == int(carry[f"{label}|dimension"])


def test_the_rebaseline_is_a_declared_delta_not_a_rewrite():
    """The carry is the pre-carve reference: its bytes are pinned by the delta file's measured values (a rewrite of the
    .npz would make every measured drift 0, and LS2 is the only fixture whose drift is 0)."""
    import json
    doc = json.loads(_REBASELINES.read_text())
    assert doc["carry"] == _CARRY.name
    measured = doc["events"][-1]["measured"]
    assert set(measured) == {row[0] for row in _SINGULAR}
    assert sum(m["apply_in_eps_max_abs_apply"] > 0 for m in measured.values()) == 5


def _trace_metric(mesh) -> np.ndarray:
    """The trace metric as the row's INPUT, through the space's metric leg (the hand sum below is the oracle)."""
    space = mesh.angular_trace
    return np.asarray(space.apply_metric(np.ones(tuple(space.shape))), dtype=float)


_THE_HUBS = object()


def _exit_reading(mesh, x, kernel=_THE_HUBS, bulk_seed=None):
    """The exit on a state whose trace is ``x``; its bulk is zero, or a seeded random field (so the full-field norm
    differs from the trace norm: U0 moved the correction onto the whole state)."""
    from dataclasses import replace
    from orpheus.sn.solver import _exit_gauge
    psi = _unwindowed_cold_start(mesh, history_depth=0)
    interior = psi.interior.values
    if bulk_seed is not None:
        interior = np.random.default_rng(bulk_seed).standard_normal(np.shape(interior))
    psi = replace(psi, interior=replace(psi.interior, values=interior), boundary=replace(psi.boundary, values=x))
    return _exit_gauge(psi, kernel=mesh.pencil.kernel if kernel is _THE_HUBS else kernel)


_UNIFORM_TRACE_METRIC = "d2 (3,4) LS2 ng=1"


@pytest.mark.foundation
def test_premise_the_exit_rows_carry_a_non_uniform_trace_metric():
    """Five of six carry fixtures read trace G max/min 2.48; LS2 reads 1.0 (a uniform G makes the G-norm ratio equal to
    the Euclidean one: that row is the control leg, never evidence)."""
    for label, cells, bcs, kwargs in _SINGULAR:
        g = _trace_metric(_mesh(cells, bcs, **kwargs))
        ratio = g.max() / g.min()
        assert (ratio == 1.0) if label == _UNIFORM_TRACE_METRIC else (ratio > 2.0), (label, ratio)


@pytest.mark.foundation
@pytest.mark.rests_on(
    "tests/gates/sn/operators/test_loss_kernel_gauge.py::test_the_rebuilt_gauge_agrees_with_the_retired_class_within_the_declared_rebaseline",
    "tests/gates/sn/operators/test_loss_kernel_gauge.py::test_premise_the_exit_rows_carry_a_non_uniform_trace_metric",
)
@pytest.mark.parametrize("label,cells,bcs,kwargs", _SINGULAR, ids=[row[0] for row in _SINGULAR])
def test_the_exit_reports_the_correction_in_the_G_norm(label, cells, bcs, kwargs):
    """RULED (U3, re-posed U0): ``||psi - sigma(psi)||_G / ||psi||_G`` in the FULL FIELD's norm (bulk + trace; until U0
    the trace's). Oracle: the hand sum over the full-field metric, the bulk a seeded random field, so a trace-norm
    denominator reads differently on every row. `[M]` before U3 the exit reported the Euclidean ratio, 5.7 %-18.6 % away
    on the five non-uniform fixtures; the bulk leg is the U0 tooth."""
    carry = np.load(_CARRY)
    mesh = _mesh(cells, bcs, **kwargs)
    g = np.asarray(mesh.full_field_space.metric_diagonal(), dtype=float)
    x = carry[f"{label}|x"]
    out, correction = _exit_reading(mesh, x, bulk_seed=4)
    before = np.concatenate([np.ravel(_exit_reading(mesh, x, kernel=None, bulk_seed=4)[0].interior.values), x])
    after = np.concatenate([np.ravel(out.interior.values), np.ravel(out.boundary.values)])
    d = before - after
    np.testing.assert_array_equal(d[: d.size - x.size], 0.0)   # the DD modes live on the trace: the bulk is untouched
    assert correction == pytest.approx(np.sqrt(np.sum(g * d * d) / np.sum(g * before * before)), rel=1e-14), label


@pytest.mark.foundation
@pytest.mark.rests_on(
    "tests/gates/numerics/test_kernel_projection_gauge.py::TestTheProtocolAndTheMetric",
)
def test_the_exit_reads_the_protocol_only():
    """A posing with no freedom leaves the state untouched and reports no correction (``None``, never "measured and
    zero"). `[M]` before U3 a Protocol conformer with no ``.blocks`` raised ``AttributeError: ... 'blocks'`` there;
    since U0 the exit takes the POSING's predicted kernel, and a trivial one is the no-freedom value."""
    from orpheus.numerics.fredholm import PredictedKernel
    mesh = _mesh((3, 4), [(_R, _R)] * 2)
    x = np.random.default_rng(0).standard_normal(np.shape(_unwindowed_cold_start(mesh, history_depth=0).boundary.values))
    out, correction = _exit_reading(mesh, x, kernel=PredictedKernel.trivial(mesh.pencil.lhs.domain, "no freedom"))
    assert correction is None and np.array_equal(np.asarray(out.boundary.values), x)
