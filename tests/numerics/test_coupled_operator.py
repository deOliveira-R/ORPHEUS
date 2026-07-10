r"""B.2a — the N-system coupled block machinery (M1–M5 gate suite).

Intrinsic-property gates for
:mod:`orpheus.numerics.coupled_system` — ``CoupledField`` (the block
vector), ``CoupledSpace`` (the direct sum + the scoped local→global offsets),
``CoupledOperator`` (the typed N×N block grid) — realizing the test-architect
spec ``coupled_operator_step4_verification.md`` rows on SYNTHETIC toys (the
machinery is semantics-agnostic; the ψ½ 2×2 is wired as instance #1 in
later Phase-B sub-steps):

* **M1** — block matvec ``y_i = Σ_j A_ij x_j`` vs the hand sum
  (``array_equal`` — the toy dense blocks add no reduction); ``None`` block =
  structural zero; 1×1 degenerate.
* **M2 (THE LOAD-BEARING ROW)** — ``assemble().as_matrix() ≡`` the
  apply-probed dense (principled-equiv ``rtol=1e-11``, NEVER ``array_equal``
  — L16 scatter-vs-apply order). ⚠ tautology trap: the probe builds basis
  ``CoupledField``\ s and runs ``.apply`` — it must NEVER route through
  ``as_matrix``/``assemble`` (which would compare assembly with itself).
  Teeth (permanent, in-process monkeypatch): transpose-placement mutant and
  dropped-block mutant both RED the gate on an asymmetric same-size grid
  (a symmetric grid is offset-transpose-blind — Mode 12).
* **M3** — per-block type safety: a space-declaring block at the wrong grid
  position is UNCONSTRUCTABLE (``IncompatibleOperatorComposition``, the
  specific wrong-column / wrong-row message + positive control); ragged /
  all-``None`` row / all-``None`` column grids rejected; wrong
  container/arity at apply; the space-less-block legacy seam skips the
  construction check and the BLOCK's own guard raises at apply (the
  wrong-component runtime raise).
* **M4 (HIGHEST-VALUE ROW)** — the Euclidean transposed-grid law
  ``(Aᵀ)_ji = (A_ij)ᵀ``; Hilbert reciprocity
  ``⟨A x, y⟩_G = ⟨x, A.H y⟩_G`` closes under the composite NONUNIFORM
  member metrics (control < 1e-12) and the M-ADJ-metric tooth (skip the
  ``G⁺/G`` conjugation in ``_AdjointOperator.apply``) REDs O(1) — the
  ERR-067 / Mode-12 reopening a Euclidean block adjoint would cause;
  non-adjointable block ⟹ eager ``MissingAdjoint``.
* **M5** — ``CoupledOperator.system_role is SystemRole.COUPLED``.

Plus the ravellable-protocol interop rider (``_is_ravellable`` /
``_zeros_like`` / ``_unravel_like`` accept a ``CoupledField``): the ERR-053
family (``restart = n_dof = to_flat().size`` GMRES sizing) closes by this
conformance, pinned here before the B.2d driver wire consumes it.

vv Mode-8 discipline: ``np.testing.assert_*`` / ``pytest.raises`` /
``pytest.fail`` only (the suite runs under ``python -O``).
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import pytest
from scipy import sparse

from orpheus.numerics.assembled_operator import SparseAssembledOperator
from orpheus.numerics.coupled_system import (
    CoupledField,
    CoupledOperator,
    CoupledSpace,
)
from orpheus.numerics.iteration import (
    _is_ravellable,
    _ravel,
    _unravel_like,
    _zeros_like,
)
from orpheus.numerics.operator import (
    IncompatibleOperatorComposition,
    LinearOperator,
    MissingAdjoint,
    MissingAssembly,
    SystemRole,
    _AdjointOperator,
)
from orpheus.numerics.space import FunctionSpace

pytestmark = pytest.mark.foundation


# ── Synthetic toys: two DISTINCT member families + carrier-typed spaces ──
#
# Two member classes with DIFFERENT shapes so offsets AND type-safety are
# both activated (the memo's "two distinct toy System types"); a same-size
# second pair (_Gamma/_Delta) for the transpose-placement tooth, where the
# swap must be shape-legal so the VALUE gate (not a shape error) is what
# reds.


class _ToyAlgebra:
    """The member contract: ravellable pair + copy + vector dunders.

    Same-class-only arithmetic — a foreign member raises, so the suite can
    assert the coupled algebra DELEGATES role law to members (the machinery
    adds no cross-member tolerance).
    """

    values: np.ndarray

    def to_flat(self) -> np.ndarray:
        return np.asarray(self.values, dtype=float).ravel()

    @classmethod
    def from_flat(cls, flat: np.ndarray, template: "_ToyAlgebra"):
        return cls(  # type: ignore[call-arg]
            values=np.asarray(flat, dtype=float).reshape(
                np.asarray(template.values).shape
            ),
        )

    def copy(self):
        return type(self)(values=np.asarray(self.values).copy())  # type: ignore[call-arg]

    def _check(self, other: object) -> "_ToyAlgebra":
        if type(other) is not type(self):
            raise TypeError(
                f"{type(self).__name__}: foreign member "
                f"{type(other).__name__}"
            )
        return other  # type: ignore[return-value]

    def __add__(self, other):
        return type(self)(values=self.values + self._check(other).values)  # type: ignore[call-arg]

    def __sub__(self, other):
        return type(self)(values=self.values - self._check(other).values)  # type: ignore[call-arg]

    def __neg__(self):
        return type(self)(values=-self.values)  # type: ignore[call-arg]

    def __mul__(self, scalar: float):
        return type(self)(values=self.values * float(scalar))  # type: ignore[call-arg]

    __rmul__ = __mul__

    def __truediv__(self, scalar: float):
        return type(self)(values=self.values / float(scalar))  # type: ignore[call-arg]


@dataclass(frozen=True)
class _AlphaField(_ToyAlgebra):
    values: np.ndarray


@dataclass(frozen=True)
class _BetaField(_ToyAlgebra):
    values: np.ndarray


@dataclass(frozen=True)
class _GammaField(_ToyAlgebra):
    values: np.ndarray


@dataclass(frozen=True)
class _DeltaField(_ToyAlgebra):
    values: np.ndarray


@dataclass(frozen=True)
class _ToySpace(FunctionSpace):
    """A carrier-typed member space: the diagonal metric on ``field.values``."""

    def apply_metric(self, x):
        return replace(x, values=self._diagonal_apply_metric(x.values))

    def apply_inverse_metric(self, x):
        return replace(x, values=self._diagonal_apply_inverse_metric(x.values))

    def inner_product(self, x, y) -> float:
        return self._diagonal_inner_product(x.values, y.values)


class _Block(LinearOperator):
    """A dense toy block: typed in/out members, adjointable, assemblable.

    The in-guard is the BLOCK-level wrong-component raise M3 exercises
    (``expected _BetaField, got _AlphaField``).
    """

    def __init__(
        self,
        matrix,
        *,
        in_cls: type,
        out_cls: type,
        domain=None,
        codomain=None,
        adjointable_: bool = True,
        assemblable_: bool = True,
    ) -> None:
        self.matrix = np.asarray(matrix, dtype=float)
        self._in_cls, self._out_cls = in_cls, out_cls
        self._domain, self._codomain = domain, codomain
        self._adjointable, self._assemblable = adjointable_, assemblable_

    @property
    def domain(self):
        return self._domain

    @property
    def codomain(self):
        return self._codomain

    def apply(self, x, /):
        if type(x) is not self._in_cls:
            raise TypeError(
                f"_Block: expected {self._in_cls.__name__}, got "
                f"{type(x).__name__}"
            )
        return self._out_cls(values=self.matrix @ x.values)

    def apply_transpose(self, y, /):
        if type(y) is not self._out_cls:
            raise TypeError(
                f"_Block transpose: expected {self._out_cls.__name__}, got "
                f"{type(y).__name__}"
            )
        return self._in_cls(values=self.matrix.T @ y.values)

    @property
    def is_adjointable(self) -> bool:
        return self._adjointable

    @property
    def is_assemblable(self) -> bool:
        return self._assemblable

    def assemble(self) -> SparseAssembledOperator:
        return SparseAssembledOperator(
            sparse.coo_array(self.matrix),
            domain=self._domain,
            codomain=self._codomain,
        )


_RNG = np.random.default_rng(20260710)

# The rectangular main fixture: alpha (3 DOF, nonuniform metric) ⊕ beta
# (2 DOF, nonuniform metric) — distinct shapes activate the offsets.
_SP_A = _ToySpace(
    name="alpha", shape=(3,),
    inner_product_weights=np.array([1.0, 2.0, 0.5]),
)
_SP_B = _ToySpace(
    name="beta", shape=(2,),
    inner_product_weights=np.array([3.0, 0.25]),
)


def _coupled_space() -> CoupledSpace:
    return CoupledSpace.from_systems((_SP_A, _SP_B))


def _grid_blocks() -> dict[str, _Block]:
    return {
        "AA": _Block(
            _RNG.random((3, 3)), in_cls=_AlphaField, out_cls=_AlphaField,
            domain=_SP_A, codomain=_SP_A,
        ),
        "AB": _Block(
            _RNG.random((3, 2)), in_cls=_BetaField, out_cls=_AlphaField,
            domain=_SP_B, codomain=_SP_A,
        ),
        "BA": _Block(
            _RNG.random((2, 3)), in_cls=_AlphaField, out_cls=_BetaField,
            domain=_SP_A, codomain=_SP_B,
        ),
        "BB": _Block(
            _RNG.random((2, 2)), in_cls=_BetaField, out_cls=_BetaField,
            domain=_SP_B, codomain=_SP_B,
        ),
    }


def _operator(blocks: dict[str, _Block] | None = None) -> CoupledOperator:
    b = blocks if blocks is not None else _grid_blocks()
    cs = _coupled_space()
    return CoupledOperator(
        [[b["AA"], b["AB"]], [b["BA"], b["BB"]]], domain=cs, codomain=cs,
    )


def _x(seed: int) -> CoupledField:
    rng = np.random.default_rng(seed)
    return CoupledField(
        systems=(
            _AlphaField(values=rng.random(3) + 0.5),
            _BetaField(values=rng.random(2) + 0.5),
        ),
    )


def _probe_dense(op: CoupledOperator, template: CoupledField) -> np.ndarray:
    """Column j = ``op.apply(e_j)`` on basis coupled fields — the BLOCK-MATVEC
    probe. Deliberately NEVER routes through ``assemble``/``as_matrix`` (the
    L16 tautology trap: probing assembly against its own densification is
    vacuous)."""
    n = template.to_flat().size
    columns = []
    for j in range(n):
        e = np.zeros(n)
        e[j] = 1.0
        columns.append(op.apply(CoupledField.from_flat(e, template)).to_flat())
    return np.column_stack(columns)


# ── CoupledField — the block vector's intrinsic laws ──────────────────


class TestCoupledFieldIntrinsic:
    def test_systems_must_be_a_nonempty_tuple(self) -> None:
        with pytest.raises(ValueError, match="at least one system"):
            CoupledField(systems=())
        with pytest.raises(TypeError, match="must be a tuple"):
            CoupledField(systems=[_x(1).systems[0]])  # type: ignore[arg-type]

    def test_member_contract_guard_names_the_offender(self) -> None:
        good = _x(2).systems[0]
        with pytest.raises(TypeError, match="system 1"):
            CoupledField(systems=(good, object()))  # type: ignore[arg-type]

    def test_arithmetic_requires_a_same_class_arity_matched_partner(self) -> None:
        x = _x(3)
        with pytest.raises(TypeError, match="same-class partner"):
            _ = x + 42  # type: ignore[operator]
        one = CoupledField(systems=(x.systems[0],))
        with pytest.raises(ValueError, match="matching system arity"):
            _ = x + one

    def test_member_wise_algebra_delegates_to_the_members(self) -> None:
        x, y = _x(4), _x(5)
        s = x + y
        np.testing.assert_array_equal(
            s.systems[0].values, x.systems[0].values + y.systems[0].values,
        )
        d = x - y
        np.testing.assert_array_equal(
            d.systems[1].values, x.systems[1].values - y.systems[1].values,
        )
        np.testing.assert_array_equal(
            (-x).systems[0].values, -x.systems[0].values,
        )
        np.testing.assert_allclose(
            (2.0 * x).systems[1].values, 2.0 * x.systems[1].values,
        )
        np.testing.assert_allclose(
            (x / 2.0).systems[0].values, 0.5 * x.systems[0].values,
        )

    def test_member_role_law_is_delegated_not_absorbed(self) -> None:
        # A partner with the member ORDER swapped is arity-legal at the
        # container but must raise from the MEMBER algebra (the machinery
        # adds no cross-member tolerance — role/type law lives on members).
        x = _x(6)
        swapped = CoupledField(systems=(x.systems[1], x.systems[0]))
        with pytest.raises(TypeError, match="foreign member"):
            _ = x + swapped

    def test_to_flat_packs_in_system_order(self) -> None:
        x = _x(7)
        np.testing.assert_array_equal(
            x.to_flat(),
            np.concatenate(
                [x.systems[0].to_flat(), x.systems[1].to_flat()],
            ),
        )

    def test_from_flat_round_trip_preserves_member_types(self) -> None:
        x = _x(8)
        back = CoupledField.from_flat(x.to_flat(), x)
        if type(back) is not CoupledField:
            pytest.fail(f"from_flat returned {type(back).__name__}")
        if type(back.systems[0]) is not _AlphaField:
            pytest.fail(f"member 0 became {type(back.systems[0]).__name__}")
        if type(back.systems[1]) is not _BetaField:
            pytest.fail(f"member 1 became {type(back.systems[1]).__name__}")
        np.testing.assert_array_equal(
            back.systems[0].values, x.systems[0].values,
        )
        np.testing.assert_array_equal(
            back.systems[1].values, x.systems[1].values,
        )

    def test_from_flat_rejects_wrong_size(self) -> None:
        x = _x(9)
        with pytest.raises(ValueError, match="does not match template total"):
            CoupledField.from_flat(x.to_flat()[:-1], x)

    def test_copy_is_value_independent(self) -> None:
        x = _x(10)
        c = x.copy()
        np.testing.assert_array_equal(
            c.systems[0].values, x.systems[0].values,
        )
        c.systems[0].values[0] += 1.0
        if x.systems[0].values[0] == c.systems[0].values[0]:
            pytest.fail("copy aliases the original member buffer")

    def test_ravellable_protocol_interop(self) -> None:
        # The ERR-053 closure rider: the Krylov boundary's duck-typed
        # helpers must accept a CoupledField, so every
        # ``n_dof = to_flat().size`` sizing site counts BOTH systems.
        x = _x(11)
        if not _is_ravellable(x):
            pytest.fail("CoupledField does not satisfy the ravellable pair")
        np.testing.assert_array_equal(_ravel(x), x.to_flat())
        z = _zeros_like(x)
        if type(z) is not CoupledField:
            pytest.fail(f"_zeros_like returned {type(z).__name__}")
        np.testing.assert_array_equal(z.to_flat(), np.zeros(x.to_flat().size))
        u = _unravel_like(x, x.to_flat())
        if type(u) is not CoupledField:
            pytest.fail(f"_unravel_like returned {type(u).__name__}")
        np.testing.assert_array_equal(u.to_flat(), x.to_flat())


# ── CoupledSpace — the direct sum + the scoped local→global offsets ───


class TestCoupledSpace:
    def test_from_systems_derives_name_and_flat_shape(self) -> None:
        cs = _coupled_space()
        if cs.name != "coupled(alpha ⊕ beta)":
            pytest.fail(f"derived name {cs.name!r}")
        if cs.shape != (5,):
            pytest.fail(f"flat direct-sum shape {cs.shape}")
        if cs.n_systems != 2:
            pytest.fail(f"n_systems {cs.n_systems}")

    def test_from_systems_rejects_empty(self) -> None:
        with pytest.raises(ValueError, match="at least one member space"):
            CoupledSpace.from_systems(())

    def test_bare_constructor_footgun_guard(self) -> None:
        bare = CoupledSpace(name="coupled(?)", shape=(5,))
        with pytest.raises(RuntimeError, match="from_systems"):
            _ = bare.n_systems

    def test_system_slices_are_the_local_to_global_offsets(self) -> None:
        # The scoped LocalToGlobalMap: the slice table IS the to_flat layout.
        cs = _coupled_space()
        if cs.system_slices != (slice(0, 3), slice(3, 5)):
            pytest.fail(f"offset table {cs.system_slices}")
        x = _x(12)
        flat = x.to_flat()
        for member, sl in zip(x.systems, cs.system_slices):
            np.testing.assert_array_equal(flat[sl], member.to_flat())

    def test_apply_metric_dispatches_member_wise(self) -> None:
        cs = _coupled_space()
        x = _x(13)
        gx = cs.apply_metric(x)
        np.testing.assert_array_equal(
            gx.systems[0].values,
            np.array([1.0, 2.0, 0.5]) * x.systems[0].values,
        )
        np.testing.assert_array_equal(
            gx.systems[1].values,
            np.array([3.0, 0.25]) * x.systems[1].values,
        )

    def test_apply_inverse_metric_round_trips(self) -> None:
        cs = _coupled_space()
        x = _x(14)
        rt = cs.apply_inverse_metric(cs.apply_metric(x))
        np.testing.assert_allclose(
            rt.systems[0].values, x.systems[0].values, rtol=1e-14,
        )
        np.testing.assert_allclose(
            rt.systems[1].values, x.systems[1].values, rtol=1e-14,
        )

    def test_inner_product_is_the_member_sum(self) -> None:
        cs = _coupled_space()
        x, y = _x(15), _x(16)
        expected = float(
            np.sum(
                np.array([1.0, 2.0, 0.5])
                * x.systems[0].values
                * y.systems[0].values
            )
            + np.sum(
                np.array([3.0, 0.25])
                * x.systems[1].values
                * y.systems[1].values
            )
        )
        np.testing.assert_allclose(cs.inner_product(x, y), expected, rtol=1e-14)

    def test_metric_arity_mismatch_raises(self) -> None:
        cs = _coupled_space()
        one = CoupledField(systems=(_x(17).systems[0],))
        with pytest.raises(ValueError, match="pairing is inconsistent"):
            cs.apply_metric(one)


# ── M1 — the block matvec ──────────────────────────────────────────────


class TestBlockMatvec:
    def test_block_matvec_matches_the_hand_sum(self) -> None:
        b = _grid_blocks()
        op = _operator(b)
        x = _x(20)
        xa, xb = x.systems[0].values, x.systems[1].values
        out = op.apply(x)
        np.testing.assert_array_equal(
            out.systems[0].values, b["AA"].matrix @ xa + b["AB"].matrix @ xb,
        )
        np.testing.assert_array_equal(
            out.systems[1].values, b["BA"].matrix @ xa + b["BB"].matrix @ xb,
        )

    def test_missing_block_is_the_structural_zero(self) -> None:
        b = _grid_blocks()
        cs = _coupled_space()
        op = CoupledOperator(
            [[b["AA"], None], [b["BA"], b["BB"]]], domain=cs, codomain=cs,
        )
        x = _x(21)
        out = op.apply(x)
        # Row A carries ONLY the AA contribution — no zero arithmetic ran.
        np.testing.assert_array_equal(
            out.systems[0].values, b["AA"].matrix @ x.systems[0].values,
        )

    def test_degenerate_single_system_grid(self) -> None:
        b = _grid_blocks()
        cs1 = CoupledSpace.from_systems((_SP_A,))
        op = CoupledOperator([[b["AA"]]], domain=cs1, codomain=cs1)
        xa = _x(22).systems[0]
        out = op.apply(CoupledField(systems=(xa,)))
        np.testing.assert_array_equal(
            out.systems[0].values, b["AA"].matrix @ xa.values,
        )

    def test_apply_rejects_wrong_container_and_arity(self) -> None:
        op = _operator()
        with pytest.raises(TypeError, match="expected a CoupledField"):
            op.apply(np.ones(5))  # type: ignore[arg-type]
        one = CoupledField(systems=(_x(23).systems[0],))
        with pytest.raises(ValueError, match="carries 1 systems"):
            op.apply(one)


# ── M2 — assemble ≡ probe (THE load-bearing gate) ──────────────────────


class TestAssembleProbe:
    def test_assemble_equals_the_apply_probe(self) -> None:
        op = _operator()
        template = _x(30)
        assembled = op.assemble().as_matrix()
        probed = _probe_dense(op, template)
        # Principled-equiv, NEVER array_equal: the sparse scatter-sum order
        # is not the block-matvec apply order (L16).
        np.testing.assert_allclose(assembled, probed, rtol=1e-11)

    def test_assembled_blocks_land_at_the_system_offsets(self) -> None:
        b = _grid_blocks()
        op = _operator(b)
        cs = _coupled_space()
        assembled = op.assemble().as_matrix()
        rows, cols = cs.system_slices, cs.system_slices
        np.testing.assert_allclose(assembled[rows[0], cols[1]], b["AB"].matrix)
        np.testing.assert_allclose(assembled[rows[1], cols[0]], b["BA"].matrix)

    def test_none_block_region_is_identically_zero(self) -> None:
        b = _grid_blocks()
        cs = _coupled_space()
        op = CoupledOperator(
            [[b["AA"], None], [b["BA"], b["BB"]]], domain=cs, codomain=cs,
        )
        assembled = op.assemble().as_matrix()
        np.testing.assert_array_equal(
            assembled[cs.system_slices[0], cs.system_slices[1]],
            np.zeros((3, 2)),
        )

    @staticmethod
    def _square_fixture() -> tuple[CoupledOperator, CoupledField, CoupledSpace]:
        """Same-size (2 ⊕ 2) ASYMMETRIC grid: the transpose-placement mutant
        is shape-legal here, so the VALUE gate (not a shape error) is what
        must red — a symmetric grid would be offset-transpose-blind (Mode 12).
        """
        sg = _ToySpace(
            name="gamma", shape=(2,),
            inner_product_weights=np.array([1.0, 4.0]),
        )
        sd = _ToySpace(
            name="delta", shape=(2,),
            inner_product_weights=np.array([0.5, 2.0]),
        )
        rng = np.random.default_rng(99)
        blocks = [
            [
                _Block(rng.random((2, 2)), in_cls=_GammaField,
                       out_cls=_GammaField, domain=sg, codomain=sg),
                _Block(rng.random((2, 2)), in_cls=_DeltaField,
                       out_cls=_GammaField, domain=sd, codomain=sg),
            ],
            [
                _Block(rng.random((2, 2)), in_cls=_GammaField,
                       out_cls=_DeltaField, domain=sg, codomain=sd),
                _Block(rng.random((2, 2)), in_cls=_DeltaField,
                       out_cls=_DeltaField, domain=sd, codomain=sd),
            ],
        ]
        cs = CoupledSpace.from_systems((sg, sd))
        op = CoupledOperator(blocks, domain=cs, codomain=cs)
        template = CoupledField(
            systems=(
                _GammaField(values=rng.random(2)),
                _DeltaField(values=rng.random(2)),
            ),
        )
        return op, template, cs

    def test_offset_swap_tooth_reds_the_probe_gate(self, monkeypatch) -> None:
        # The M2 tooth: a mutant assemble that scatters block (i, j) at the
        # TRANSPOSED placement (j, i). In-process monkeypatch (never
        # git-checkout); the asymmetric same-size grid makes the swap
        # shape-legal and value-visible O(1).
        op, template, _ = self._square_fixture()

        def _swapped_assemble(self: CoupledOperator) -> SparseAssembledOperator:
            grid = [
                [
                    None
                    if self.blocks[j][i] is None
                    else self.blocks[j][i].assemble().matrix  # type: ignore[union-attr]
                    for j in range(self.n_cols)
                ]
                for i in range(self.n_rows)
            ]
            return SparseAssembledOperator(
                sparse.block_array(grid), domain=self.domain,
                codomain=self.codomain,
            )

        monkeypatch.setattr(CoupledOperator, "assemble", _swapped_assemble)
        assembled = op.assemble().as_matrix()
        probed = _probe_dense(op, template)
        defect = float(np.abs(assembled - probed).max())
        if defect < 1e-3:
            pytest.fail(
                f"offset-swap mutant survived the assemble≡probe gate "
                f"(max diff {defect:.2e}) — the gate has no placement teeth"
            )

    def test_dropped_block_tooth_reds_the_probe_gate(self, monkeypatch) -> None:
        op, template, _ = self._square_fixture()

        def _dropping_assemble(self: CoupledOperator) -> SparseAssembledOperator:
            grid = [
                [
                    None
                    if (b is None or (i, j) == (0, 1))
                    else b.assemble().matrix  # type: ignore[union-attr]
                    for j, b in enumerate(row)
                ]
                for i, row in enumerate(self.blocks)
            ]
            return SparseAssembledOperator(
                sparse.block_array(grid), domain=self.domain,
                codomain=self.codomain,
            )

        monkeypatch.setattr(CoupledOperator, "assemble", _dropping_assemble)
        assembled = op.assemble().as_matrix()
        probed = _probe_dense(op, template)
        defect = float(np.abs(assembled - probed).max())
        if defect < 1e-3:
            pytest.fail(
                f"dropped-block mutant survived the assemble≡probe gate "
                f"(max diff {defect:.2e})"
            )

    def test_assemble_names_a_non_assemblable_block(self) -> None:
        b = _grid_blocks()
        b["AB"] = _Block(
            b["AB"].matrix, in_cls=_BetaField, out_cls=_AlphaField,
            domain=_SP_B, codomain=_SP_A, assemblable_=False,
        )
        op = _operator(b)
        if op.is_assemblable:
            pytest.fail("is_assemblable must derive False from the block")
        with pytest.raises(MissingAssembly, match=r"block \(0, 1\)"):
            op.assemble()


# ── M3 — per-block type safety (net-new teeth) ─────────────────────────


class TestTypeSafety:
    def test_wrong_column_placement_is_unconstructable(self) -> None:
        b = _grid_blocks()
        cs = _coupled_space()
        # Positive control FIRST: the correct grid constructs and applies.
        op = CoupledOperator(
            [[b["AA"], b["AB"]], [b["BA"], b["BB"]]], domain=cs, codomain=cs,
        )
        out = op.apply(_x(40))
        if type(out) is not CoupledField:
            pytest.fail("positive control failed to apply")
        # The negative: AB and BA swapped — a space-declaring block at the
        # wrong grid position must be unconstructable, with the specific
        # wrong-column message (a bare raises would false-green on any
        # downstream shape crash).
        with pytest.raises(
            IncompatibleOperatorComposition,
            match=r"block \(0, 1\).*wrong column",
        ):
            CoupledOperator(
                [[b["AA"], b["BA"]], [b["AB"], b["BB"]]],
                domain=cs, codomain=cs,
            )

    def test_wrong_row_placement_names_the_row(self) -> None:
        b = _grid_blocks()
        cs = _coupled_space()
        # Domain-legal but codomain-wrong: swap the grid ROWS. Block (0, 0)
        # becomes BA (domain alpha ✓ column 0, codomain beta ✗ row 0).
        with pytest.raises(
            IncompatibleOperatorComposition,
            match=r"block \(0, 0\).*wrong row",
        ):
            CoupledOperator(
                [[b["BA"], b["BB"]], [b["AA"], b["AB"]]],
                domain=cs, codomain=cs,
            )

    def test_grid_shape_must_match_the_spaces(self) -> None:
        b = _grid_blocks()
        cs = _coupled_space()
        with pytest.raises(ValueError, match="one row per codomain system"):
            CoupledOperator([[b["AA"], b["AB"]]], domain=cs, codomain=cs)
        with pytest.raises(ValueError, match="one column per domain system"):
            CoupledOperator(
                [[b["AA"]], [b["BA"]]], domain=cs, codomain=cs,
            )

    def test_all_none_row_and_column_rejected(self) -> None:
        b = _grid_blocks()
        cs = _coupled_space()
        with pytest.raises(ValueError, match="row 1 has no blocks"):
            CoupledOperator(
                [[b["AA"], b["AB"]], [None, None]], domain=cs, codomain=cs,
            )
        with pytest.raises(ValueError, match="column 1 has no blocks"):
            CoupledOperator(
                [[b["AA"], None], [b["BA"], None]], domain=cs, codomain=cs,
            )

    def test_space_less_block_skips_the_check_and_raises_at_apply(self) -> None:
        # The legacy seam (house convention: a None-spaced operand skips the
        # composability check). The mis-placed space-LESS block constructs —
        # and the BLOCK's own in-guard raises at apply time (the M3
        # wrong-component runtime raise, block-level).
        b = _grid_blocks()
        cs = _coupled_space()
        naked_ba = _Block(
            b["BA"].matrix, in_cls=_AlphaField, out_cls=_BetaField,
        )
        op = CoupledOperator(
            [[b["AA"], naked_ba], [b["BA"], b["BB"]]], domain=cs, codomain=cs,
        )
        with pytest.raises(TypeError, match="expected _AlphaField"):
            op.apply(_x(41))


# ── M4 — the block adjoint (Mode-12 closure; highest-value row) ────────


class TestBlockAdjoint:
    def test_transposed_grid_law(self) -> None:
        b = _grid_blocks()
        op = _operator(b)
        y = _x(50)
        t = op.apply_transpose(y)
        ya, yb = y.systems[0].values, y.systems[1].values
        np.testing.assert_array_equal(
            t.systems[0].values,
            b["AA"].matrix.T @ ya + b["BA"].matrix.T @ yb,
        )
        np.testing.assert_array_equal(
            t.systems[1].values,
            b["AB"].matrix.T @ ya + b["BB"].matrix.T @ yb,
        )

    def test_hilbert_reciprocity_closes_under_the_composite_metric(self) -> None:
        # The CONTROL leg (L18 both-legs discipline; the tooth below is the
        # mutated leg): ⟨A x, y⟩_G = ⟨x, A.H y⟩_G under the NONUNIFORM
        # member metrics — the metric conjugation is realized by the generic
        # _AdjointOperator over CoupledSpace.apply_metric/-inverse_metric.
        op = _operator()
        cs = _coupled_space()
        x, y = _x(51), _x(52)
        lhs = cs.inner_product(op.apply(x), y)
        rhs = cs.inner_product(x, op.H.apply(y))
        defect = abs(lhs - rhs) / abs(lhs)
        if defect > 1e-12:
            pytest.fail(f"reciprocity defect {defect:.2e} exceeds 1e-12")

    def test_euclidean_block_adjoint_tooth_reds(self, monkeypatch) -> None:
        # The M-ADJ-metric tooth: skip the G⁺/G conjugation in
        # _AdjointOperator.apply (the Euclidean block adjoint) — the exact
        # ERR-067 / Mode-12 reopening. Must red O(1) under the nonuniform
        # composite metric.
        op = _operator()
        cs = _coupled_space()
        x, y = _x(53), _x(54)
        monkeypatch.setattr(
            _AdjointOperator, "apply",
            lambda self, z: self.inner.apply_transpose(z),
        )
        lhs = cs.inner_product(op.apply(x), y)
        rhs = cs.inner_product(x, op.H.apply(y))
        defect = abs(lhs - rhs) / abs(lhs)
        if defect < 1e-3:
            pytest.fail(
                f"Euclidean block adjoint survived the reciprocity gate "
                f"(defect {defect:.2e}) — the metric conjugation is not "
                f"load-bearing"
            )

    def test_non_adjointable_block_refuses_eagerly(self) -> None:
        b = _grid_blocks()
        b["BA"] = _Block(
            b["BA"].matrix, in_cls=_AlphaField, out_cls=_BetaField,
            domain=_SP_A, codomain=_SP_B, adjointable_=False,
        )
        op = _operator(b)
        if op.is_adjointable:
            pytest.fail("is_adjointable must derive False from the block")
        with pytest.raises(MissingAdjoint):
            _ = op.H
        with pytest.raises(MissingAdjoint, match=r"block \(1, 0\)"):
            op.apply_transpose(_x(55))

    def test_adjoint_swaps_the_coupled_spaces(self) -> None:
        op = _operator()
        adj = op.H
        if adj.domain is not op.codomain or adj.codomain is not op.domain:
            pytest.fail("A.H must present the swapped coupled spaces")


# ── M5 — the system-role stamp ─────────────────────────────────────────


class TestSystemRoleStamp:
    def test_coupled_operator_is_stamped_coupled(self) -> None:
        if CoupledOperator.system_role is not SystemRole.COUPLED:
            pytest.fail("class stamp missing")
        if _operator().system_role is not SystemRole.COUPLED:
            pytest.fail("instance stamp missing")
