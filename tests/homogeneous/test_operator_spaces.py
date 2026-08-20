r"""The homogeneous operators pose on ONE real space (campaign 1, CS1 3b).

The positive FLOOR that succeeded the four ``test_monomorphic_leaves``
strict-xfail rows (``test_model_generic_leaf_declares_a_space[C|F-2g|4g]``,
deleted at 3b), plus the refusal witnesses and the vv#19 loaded/blind
``.H`` pair. Gate ids D1–D11 refer to the CS1 battery of record
(``scratch/cs1_verification_plan.md`` §2).

⭐⭐ HOME: deliberately NOT ``test_homogeneous.py`` — that module's
``pytestmark = [l1, verifies(<19 labels>)]`` would write FALSE
equation-TESTS edges for every foundation invariant added there (the
in-tree precedent is ``tests/numerics/test_matrix_inverse_operator.py``,
hosted out of that file for exactly this reason). These are software
invariants of the posing, so ``foundation``, never ``verifies(...)``.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.homogeneous.solver import _assemble_loss_operator
from orpheus.numerics.axis import Axis, BasisKind, EnergyAxis
from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator
from orpheus.numerics.operator import (
    IncompatibleOperatorComposition,
    MissingAssembly,
)
from orpheus.numerics.space import FunctionSpace
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.operators.fission import FissionOperator
from orpheus.transport.operators.isotropic_scattering import (
    IsotropicN2N,
    IsotropicScattering,
)
from orpheus.transport.operators.multiplication_operator import (
    MultiplicationOperator,
)

pytestmark = pytest.mark.foundation

_EDGES_2G = np.array([1.0e7, 1.0e3, 1.0e-3])
_EDGES_4G = np.array([1.0e7, 1.0e5, 1.0e3, 1.0, 1.0e-3])


def _require(condition: bool, message: str) -> None:
    """A ``-O``-firing assertion (NOT a bare ``assert``)."""
    if not condition:
        pytest.fail(message)


def _mat_xs(groups: str, edges: np.ndarray | None = None):
    """The meshless carrier exactly as ``solve_homogeneous_infinite`` builds it."""
    mix = get_mixture("A", groups)
    if edges is not None:
        mix = dataclasses.replace(mix, eg=edges)
    return MaterialMesh.from_materials({0: mix}).material_xs_field()


def _fused_loss_matrix(mat_xs) -> np.ndarray:
    """The independent loss reference ``diag(Σt) − (Σs0 + 2Σ2)ᵀ`` from raw XS."""
    sig_t = mat_xs.total_cross_section[:, 0]
    sig_s0 = mat_xs.sig_s_legendre(0)[0]
    sig_2 = mat_xs.n2n_matrix(0)
    return np.diag(sig_t) - (sig_s0 + 2.0 * sig_2).T


@pytest.mark.parametrize("groups", ["2g", "4g"])
@pytest.mark.parametrize("with_eg", [False, True], ids=["synthetic", "from_grid"])
def test_every_homogeneous_operator_reports_the_same_space(
    groups: str, with_eg: bool
) -> None:
    r"""D1 ⭐ — the FLOOR: C, IsoS, IsoN2N, F, M⁻¹, K all pose on ONE space.

    Successor of the four deleted R1 xfail rows. Mirrors the production
    construction (``_assemble_loss_operator`` IS the production assembler;
    F mirrors ``solver.py``'s ``space=mat_xs.mesh.bulk_space`` spelling)
    and asserts the whole posing agrees on the carrier's axis-built
    Energy ⊗ point space: shape ``(ng, 1)``, all-NODAL (coordinate cone
    present), energy arm per the carrier's data.
    """
    edges = {False: None, True: {"2g": _EDGES_2G, "4g": _EDGES_4G}[groups]}[with_eg]
    mat_xs = _mat_xs(groups, edges)
    ng = mat_xs.mesh.ng
    space = mat_xs.mesh.bulk_space

    _require(space.shape == (ng, 1), f"bulk_space shape {space.shape} != ({ng}, 1)")
    _require(space.has_coordinate_cone is True, "the scalar bulk is all-NODAL")
    axes = space.axes
    assert axes is not None
    energy_axis = axes[0]
    _require(
        isinstance(energy_axis, EnergyAxis),
        "the first factor must be the EnergyAxis",
    )
    assert isinstance(energy_axis, EnergyAxis)
    _require(
        (energy_axis.edges is not None) == with_eg,
        f"energy arm wrong: edges={'present' if energy_axis.edges is not None else 'absent'} "
        f"for the {'eg-bearing' if with_eg else 'grid-less'} carrier",
    )
    _require(
        axes[1].weights is None,
        "the quotient point's unit volume must canonicalize to the counting "
        "weight (the normalized density convention)",
    )

    loss = _assemble_loss_operator(mat_xs)
    production = FissionOperator.from_solver_data(
        mat_xs=mat_xs, space=mat_xs.mesh.bulk_space,
    )
    inverse = MatrixInverseOperator(loss)
    K = inverse @ production
    operators = {
        "C": MultiplicationOperator.from_mesh(
            mat_xs.total_cross_section_field, mat_xs.mesh,
        ),
        "IsoS": IsotropicScattering(mat_xs, space=space),
        "IsoN2N": IsotropicN2N(mat_xs, space=space),
        "F": production,
        "loss": loss,
        "M_inv": inverse,
        "K": K,
    }
    for name, op in operators.items():
        _require(
            op.domain == space and op.codomain == space,
            f"{name}: domain={op.domain!r}, codomain={op.codomain!r} — the "
            f"posing does not agree on the carrier's bulk_space",
        )


def test_two_group_and_four_group_sum_is_REFUSED() -> None:
    """D2 — the refusal witness the space threading ACTIVATES: an ill-posed
    cross-group sum, both arms bound, dies at construction with the
    established provenance fragment (the ``_agreed_space`` pins' shared
    vocabulary — do not invent new wording)."""
    c_2g = MultiplicationOperator.from_mesh(
        _mat_xs("2g").total_cross_section_field, _mat_xs("2g").mesh,
    )
    c_4g = MultiplicationOperator.from_mesh(
        _mat_xs("4g").total_cross_section_field, _mat_xs("4g").mesh,
    )
    _require(
        c_2g.domain is not None and c_4g.domain is not None,
        "precondition lost: the chain no longer binds the degenerate carrier",
    )
    with pytest.raises(IncompatibleOperatorComposition, match="equal domains"):
        _ = c_2g + c_4g


def test_matrix_inverse_of_2g_loss_composed_with_4g_fission_is_REFUSED() -> None:
    """D3 — the product-guard witness: ``M⁻¹(2g) @ F(4g)`` dies at
    construction naming the composition law."""
    loss_2g = _assemble_loss_operator(_mat_xs("2g"))
    mat_4g = _mat_xs("4g")
    f_4g = FissionOperator.from_solver_data(
        mat_xs=mat_4g, space=mat_4g.mesh.bulk_space,
    )
    with pytest.raises(
        IncompatibleOperatorComposition, match=r"A\.domain == B\.codomain"
    ):
        _ = MatrixInverseOperator(loss_2g) @ f_4g


def test_H_is_bit_identical_to_the_pre_CS1_euclidean_transpose() -> None:
    r"""D4a ⭐ — the vv#19 NEUTRALITY leg (the LOADED leg is D4b).

    The threaded ``bulk_space`` carries identity metrics BY THE COUNTING
    THEOREM (group integrals × group averages pair without widths ⟹
    energy metric = I) and the quotient point's unit volume canonicalizes
    to the counting weight — so the threaded ``loss.H`` must stay
    BIT-identical to the pre-CS1 path. The comparison is the vv#12 direct
    form: the SAME loss built space-less (the None path — still legal
    until CS4) versus the production threaded one, ``np.array_equal``
    (an independent matmul reference would associate differently and
    fail at 1 ULP — measured; the *value* claim rides on D9's fused
    matrix at ``atol=1e-12``). ⚠ And by the SCALAR-COMMUTATOR argument
    (F2, measured): even a non-unit uniform cell volume could not move
    ``.H`` — ``G = cI`` commutes with everything — so this gate alone can
    never certify the metric plumbing; D4b's control is the loaded leg.
    """
    mat_xs = _mat_xs("2g")
    loss_threaded = _assemble_loss_operator(mat_xs)
    loss_bare = MultiplicationOperator(
        coefficient=mat_xs.total_cross_section_field,
    ) - (IsotropicScattering(mat_xs) + IsotropicN2N(mat_xs))
    x = np.array([[1.0], [2.0]])
    got = np.asarray(loss_threaded.H.apply(x))
    old = np.asarray(loss_bare.H.apply(x))
    _require(
        bool(np.array_equal(got, old)),
        f"threaded loss.H moved off the pre-CS1 Euclidean-transpose path: "
        f"{got.ravel()} != {old.ravel()} — the counting theorem promises "
        f"bit-identity",
    )


def test_H_MOVES_under_a_per_group_weighted_axis() -> None:
    r"""D4b ⭐ — the vv#19 CONTROL: a per-GROUP energy weight moves ``.H``.

    ⚠ Deliberately NON-PHYSICAL: the counting-measure theorem forbids a
    weighted energy axis on a real problem (``EnergyAxis`` REFUSES
    weights at construction), so the toy uses a generic ``Axis`` — its
    whole job is to prove the adjoint machinery actually consults the
    threaded metric (``[M]`` component 0 moves ~4.75×; the whole vector
    is asserted). Without this leg, D4a's green is indistinguishable
    from a blind gate (vv#19: only the deliberately-wrong structure
    discriminates loaded from blind).
    """
    mat_xs = _mat_xs("2g")
    w_energy = np.array([2.0, 5.0])
    weighted_space = FunctionSpace.of_axes(
        Axis("energy", (2,), weights=w_energy, kind=BasisKind.NODAL),
        Axis("spatial", (1,), kind=BasisKind.NODAL),
    )
    collision = MultiplicationOperator.from_mesh(
        mat_xs.total_cross_section_field, mat_xs.mesh, space=weighted_space,
    )
    k_iso = IsotropicScattering(mat_xs, space=weighted_space) + IsotropicN2N(
        mat_xs, space=weighted_space,
    )
    loss_weighted = collision - k_iso

    x = np.array([[1.0], [2.0]])
    a_fused = _fused_loss_matrix(mat_xs)
    w = w_energy.reshape(2, 1)
    # The Hilbert adjoint under G = diag(w) ⊗ 1: G⁻¹ Aᵀ (G x), built
    # independently from the fused matrix.
    reference = ((a_fused.T @ (w * x).ravel()).reshape(2, 1)) / w
    euclidean = (a_fused.T @ x.ravel()).reshape(2, 1)
    got = np.asarray(loss_weighted.H.apply(x))
    _require(
        bool(np.allclose(got, reference, rtol=1e-14, atol=0.0)),
        f"weighted .H disagrees with the independent G⁻¹AᵀG reference: "
        f"{got.ravel()} != {reference.ravel()}",
    )
    _require(
        not np.allclose(got, euclidean, rtol=1e-10, atol=0.0),
        "the control did not MOVE — the adjoint machinery is not "
        "consulting the threaded metric (a blind gate, vv#19)",
    )


def test_bulk_space_energy_arm_distinguishes_from_grid_from_synthetic() -> None:
    r"""D6 ⭐ — the SAME mixture with and without ``eg`` mints UNEQUAL
    spaces (F1's production witness, and B2's injectivity in production:
    same shape, different partition ⟹ different name ⟹ different space).
    """
    bare = _mat_xs("2g").mesh.bulk_space
    gridded = _mat_xs("2g", _EDGES_2G).mesh.bulk_space
    _require(bare.shape == gridded.shape, "precondition: same index set")
    _require(
        bare != gridded,
        "the energy partition was dropped from the space identity — a "
        "gridded and a grid-less problem would compose silently",
    )


def test_assemble_refusal_names_the_real_reason() -> None:
    """D8 — the homogeneous C is now space-BEARING, so the assemble
    refusal may no longer claim it is "space-anonymous"; it must name the
    real reason (no composite flat layout on a plain bulk space)."""
    mat_xs = _mat_xs("2g")
    collision = MultiplicationOperator.from_mesh(
        mat_xs.total_cross_section_field, mat_xs.mesh,
    )
    _require(collision.space is not None, "precondition lost: C is space-less")
    with pytest.raises(MissingAssembly, match="composite flat layout") as exc:
        collision.assemble()
    _require(
        "space-anonymous" not in str(exc.value),
        f"the refusal still claims space-anonymity of a space-bearing "
        f"multiplier: {exc.value}",
    )


def test_as_matrix_derives_the_basis_shape_from_the_threaded_domain() -> None:
    r"""D9 — with both production ``basis_shape=(ng, 1)`` spellings
    deleted, the bare ``as_matrix()`` derivation is LOAD-BEARING: it must
    reproduce the independent fused loss matrix.

    (⚠ battery M23: a leftover explicit spelling is value-identical, so
    "the spellings are gone" has no runtime witness — that half is a grep
    obligation on the 3b commit, not this gate's claim.)
    """
    mat_xs = _mat_xs("2g")
    loss = _assemble_loss_operator(mat_xs)
    got = loss.as_matrix()
    _require(
        bool(np.allclose(got, _fused_loss_matrix(mat_xs), rtol=0.0, atol=1e-12)),
        "bare as_matrix() (deriving (ng, 1) from the threaded domain) "
        "disagrees with the independent fused loss matrix",
    )


def test_from_solver_data_does_NOT_default_derive_a_space() -> None:
    """D11 — pin the ruled ABSENCE: a bare ``from_solver_data(mat_xs=…)``
    stays space-less.

    The wrong-family hazard: ``mat_xs`` cannot know which family's space
    the caller poses on (an SN caller threads the angular composite, the
    homogeneous caller its scalar bulk) — a silently-added default
    derivation would pick one and type-check the other family's posing
    against it. Nothing else reds if that default appears (the rows that
    would have caught it are the ones CS1 deleted), so this gate is the
    only witness. WHEN CS4's Optional→mandatory flip lands: delete this
    gate (the constructor will then REQUIRE the space).
    """
    op = FissionOperator.from_solver_data(mat_xs=_mat_xs("2g"))
    _require(
        op.domain is None and op.codomain is None,
        "from_solver_data silently derived a space — the wrong-family "
        "hazard is live",
    )
