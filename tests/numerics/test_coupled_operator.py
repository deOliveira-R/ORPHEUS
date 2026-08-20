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
    CoupledSubstitutionOperator,
)
from orpheus.numerics.iteration import (
    _is_ravellable,
    _ravel,
    _unravel_like,
    _zeros_like,
)
from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator
from orpheus.numerics.operator import (
    IncompatibleOperatorComposition,
    LinearOperator,
    MatrixTooLarge,
    MissingAdjoint,
    MissingAssembly,
    NotInvertible,
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
        invertible_: bool = True,
    ) -> None:
        self.matrix = np.asarray(matrix, dtype=float)
        self._in_cls, self._out_cls = in_cls, out_cls
        self._domain, self._codomain = domain, codomain
        self._adjointable, self._assemblable = adjointable_, assemblable_
        self._invertible = invertible_

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

    # The direct-invertible face (the step-5 substitution consumes it):
    # ``solve`` inverts the forward (out → in), ``solve_transpose`` the
    # transposed forward (in → out) — with the same class guards as
    # apply/apply_transpose, so a member routed to the wrong verb raises.

    @property
    def is_invertible(self) -> bool:
        return self._invertible

    def solve(self, b):
        if type(b) is not self._out_cls:
            raise TypeError(
                f"_Block solve: expected {self._out_cls.__name__}, got "
                f"{type(b).__name__}"
            )
        return self._in_cls(values=np.linalg.solve(self.matrix, b.values))

    def solve_transpose(self, b):
        if type(b) is not self._in_cls:
            raise TypeError(
                f"_Block solve_transpose: expected {self._in_cls.__name__}, "
                f"got {type(b).__name__}"
            )
        return self._out_cls(values=np.linalg.solve(self.matrix.T, b.values))


class _NoReverseBlock(_Block):
    """A diagonal block WITHOUT the direct transpose-solve verb (the
    ``solve_transpose = None`` class attr makes ``callable(getattr(...))``
    False) — the transposed-substitution MissingAdjoint leg's fixture."""

    solve_transpose = None  # type: ignore[assignment]


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

    def test_principal_bulk_leaf_delegates_to_the_first_system(self) -> None:
        """The coupled block vector's convergence-diagnostic carrier is the
        PRINCIPAL (first) system's own ``principal_bulk_leaf`` (CS3-R —
        the same system-order convention that fixes the flat layout).
        The toy members expose none, so the duck delegation degrades to
        ``None`` (no diagnostics — the #391 member contract does not yet
        declare the property); the real chain (CoupledField → FullField →
        interior) is pinned end-to-end by
        tests/sn/operators/test_psi_half_coupling.py G-d1.7."""
        x = _x(9)
        assert x.principal_bulk_leaf is None

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


# ── Step 5 — the block SOLVE mode (S1–S4 + the space zeros seam) ───────
#
# The step-5 gate spec rows (``coupled_operator_step5_solve_verification.md``):
# S1 structural triangularity detection + route keying, S2 the substitution
# value + the order/drop/sign teeth, S3 the transposed substitution vs the
# dense ``Mᵀ``, S4 the structure-keyed ``inverse()`` factory — on the SAME
# synthetic toys as M1–M5. Every fixture grid is ASYMMETRIC (random blocks;
# a symmetric grid is order/transpose-blind — Mode 12), with diagonally-
# dominant diagonal blocks (``+3I``) so the LU/dense references are
# well-conditioned. A LOCAL rng — the module-level ``_RNG`` stream feeding
# the M1–M5 fixtures stays untouched.


_SOLVE_RNG = np.random.default_rng(20260712)

_D11 = _SOLVE_RNG.random((3, 3)) + 3.0 * np.eye(3)
_D12 = _SOLVE_RNG.random((3, 2))
_D21 = _SOLVE_RNG.random((2, 3))
_D22 = _SOLVE_RNG.random((2, 2)) + 3.0 * np.eye(2)


def _solve_blocks(*, assemblable: bool = True, **overrides) -> dict:
    blocks: dict = {
        "AA": _Block(
            _D11, in_cls=_AlphaField, out_cls=_AlphaField,
            domain=_SP_A, codomain=_SP_A, assemblable_=assemblable,
        ),
        "AB": _Block(
            _D12, in_cls=_BetaField, out_cls=_AlphaField,
            domain=_SP_B, codomain=_SP_A, assemblable_=assemblable,
        ),
        "BA": _Block(
            _D21, in_cls=_AlphaField, out_cls=_BetaField,
            domain=_SP_A, codomain=_SP_B, assemblable_=assemblable,
        ),
        "BB": _Block(
            _D22, in_cls=_BetaField, out_cls=_BetaField,
            domain=_SP_B, codomain=_SP_B, assemblable_=assemblable,
        ),
    }
    blocks.update(overrides)
    return blocks


def _zeros_field() -> CoupledField:
    return CoupledField(
        systems=(
            _AlphaField(values=np.zeros(3)),
            _BetaField(values=np.zeros(2)),
        ),
    )


def _solvable_space(*, zeros: bool = True) -> CoupledSpace:
    return CoupledSpace.from_systems(
        (_SP_A, _SP_B), zeros=_zeros_field if zeros else None,
    )


_PATTERNS: dict = {
    "upper": [["AA", "AB"], [None, "BB"]],
    "lower": [["AA", None], ["BA", "BB"]],
    "full": [["AA", "AB"], ["BA", "BB"]],
    "antidiag": [[None, "AB"], ["BA", None]],
}


def _grid(
    pattern: str, *, zeros: bool = True, assemblable: bool = True, **overrides,
) -> CoupledOperator:
    """Build the ``pattern`` ∈ {upper, lower, full, antidiag} grid — the
    ``None``-mask IS the structure the detector reads."""
    blocks = _solve_blocks(assemblable=assemblable, **overrides)
    cs = _solvable_space(zeros=zeros)
    rows = [
        [blocks[key] if key is not None else None for key in row]
        for row in _PATTERNS[pattern]
    ]
    return CoupledOperator(rows, domain=cs, codomain=cs)


def _rhs(seed: int) -> CoupledField:
    rng = np.random.default_rng(seed)
    return CoupledField(
        systems=(
            _AlphaField(values=rng.random(3) + 0.5),
            _BetaField(values=rng.random(2) + 0.5),
        ),
    )


class TestBlockSolve:
    # ── S1 — triangularity is STRUCTURAL (the None-pattern) ───────────

    def test_triangular_orientation_reads_the_none_pattern(self) -> None:
        if _grid("upper")._triangular_orientation() != "upper":
            pytest.fail("upper grid not detected")
        if _grid("lower")._triangular_orientation() != "lower":
            pytest.fail("lower grid not detected")
        if _grid("full")._triangular_orientation() is not None:
            pytest.fail(
                "a FULL grid read as triangular — the substitution would "
                "silently drop a coupling (the S1 tooth)"
            )
        if _grid("antidiag")._triangular_orientation() is not None:
            pytest.fail("a diagonal-less grid read as triangular")

    def test_is_invertible_routes(self) -> None:
        if not _grid("upper").is_invertible:
            pytest.fail("triangular substitution route not advertised")
        if not _grid("full").is_invertible:
            pytest.fail("square materializable grid not advertised")
        # Falls back to the zeros-factory probe when assembly is absent:
        bare = _grid(
            "full", zeros=True,
            AA=_Block(_D11, in_cls=_AlphaField, out_cls=_AlphaField,
                      domain=_SP_A, codomain=_SP_A, assemblable_=False),
        )
        if not bare.is_invertible:
            pytest.fail("zeros-factory materialization route not advertised")
        # No route at all: non-assemblable + zeroless space.
        no_route = _grid("full", zeros=False, assemblable=False)
        if no_route.is_invertible:
            pytest.fail("a route was advertised with no assembly, no zeros")
        with pytest.raises(NotInvertible, match="no direct route"):
            no_route.inverse()
        with pytest.raises(NotInvertible, match="no direct route"):
            no_route.solve(_rhs(1))

    # ── S2 — the substitution value + order/drop/sign teeth ───────────

    def test_upper_substitution_bit_matches_the_hand_substitution(self) -> None:
        grid, q = _grid("upper"), _rhs(11)
        x = grid.solve(q)
        blocks = _solve_blocks()
        x_b = blocks["BB"].solve(q.systems[1])
        x_a = blocks["AA"].solve(q.systems[0] - blocks["AB"].apply(x_b))
        np.testing.assert_array_equal(x.systems[1].values, x_b.values)
        np.testing.assert_array_equal(x.systems[0].values, x_a.values)

    def test_substitution_matches_the_dense_lu_reference(self) -> None:
        template = _zeros_field()
        for pattern in ("upper", "lower"):
            grid, q = _grid(pattern), _rhs(12)
            reference = np.linalg.solve(
                _probe_dense(grid, template), q.to_flat(),
            )
            np.testing.assert_allclose(
                grid.solve(q).to_flat(), reference, rtol=1e-12, atol=1e-14,
                err_msg=f"{pattern} substitution off the LU reference",
            )

    def test_substitution_drop_and_sign_teeth(self, monkeypatch) -> None:
        """The off-diagonal update ``rhs_i − A_ij·x_j`` has teeth: dropping
        the subtraction or flipping its sign REDs against the dense-LU
        reference O(1) (control leg first — L18 both-legs discipline)."""
        grid, q, template = _grid("upper"), _rhs(13), _zeros_field()
        reference = np.linalg.solve(_probe_dense(grid, template), q.to_flat())
        np.testing.assert_allclose(  # control: unmutated matches
            grid.solve(q).to_flat(), reference, rtol=1e-12, atol=1e-14,
        )

        def _mutant(sign: float):
            def _body(self, rhs, orientation, *, transpose):
                n = self.n_rows
                order = (
                    range(n - 1, -1, -1)
                    if (orientation == "upper") != transpose
                    else range(n)
                )
                solved: list = [None] * n
                for i in order:
                    acc = rhs.systems[i]
                    if sign != 0.0:
                        for j in range(n):
                            if j == i:
                                continue
                            block = self.blocks[i][j]
                            if block is None or solved[j] is None:
                                continue
                            update = block.apply(solved[j])
                            acc = (
                                acc - update if sign > 0 else acc + update
                            )
                    diagonal = self.blocks[i][i]
                    solved[i] = diagonal.solve(acc)
                return CoupledField(systems=tuple(solved))

            return _body

        for name, sign in (("drop", 0.0), ("sign-flip", -1.0)):
            with monkeypatch.context() as m:
                m.setattr(CoupledOperator, "_solve_triangular", _mutant(sign))
                mutated = grid.solve(q).to_flat()
            if np.allclose(mutated, reference, rtol=1e-6):
                pytest.fail(
                    f"the {name} mutation left the substitution on the "
                    f"reference — the off-diagonal update has no teeth"
                )

    def test_substitution_order_guard_is_loud(self, monkeypatch) -> None:
        """Visiting members in the WRONG order must hit the loud
        ordering-bug guard, never silently drop the coupling (the
        production RuntimeError is itself mutation-verified here)."""
        grid, q = _grid("upper"), _rhs(14)
        original = CoupledOperator._solve_triangular

        def _flipped(self, rhs, orientation, *, transpose):
            wrong = "lower" if orientation == "upper" else "upper"
            return original(self, rhs, wrong, transpose=transpose)

        monkeypatch.setattr(CoupledOperator, "_solve_triangular", _flipped)
        with pytest.raises(RuntimeError, match="ordering bug"):
            grid.solve(q)

    # ── S3 — the transposed substitution ──────────────────────────────

    def test_transpose_substitution_matches_the_dense_transpose(self) -> None:
        template = _zeros_field()
        for pattern in ("upper", "lower"):
            grid, b = _grid(pattern), _rhs(15)
            dense = _probe_dense(grid, template)
            if np.allclose(dense, dense.T):
                pytest.fail("fixture drift: the grid is symmetric — the "
                            "transpose gate is Mode-12-blind")
            reference = np.linalg.solve(dense.T, b.to_flat())
            np.testing.assert_allclose(
                grid.solve_transpose(b).to_flat(), reference,
                rtol=1e-12, atol=1e-14,
                err_msg=f"{pattern} transposed substitution off dense-Mᵀ",
            )
            np.testing.assert_array_equal(
                grid.inverse().apply_transpose(b).to_flat(),
                grid.solve_transpose(b).to_flat(),
                err_msg="the wrap's apply_transpose must delegate to "
                        "solve_transpose (one body)",
            )

    def test_transpose_guards_and_the_wrap_predicate(self) -> None:
        clean = _grid("upper")
        if not clean.inverse().is_adjointable:
            pytest.fail("positive control: the clean upper grid's "
                        "substitution wrap must advertise the transpose")
        # (a) a non-adjointable coupling block — the per-block guard.
        no_adj = _grid(
            "upper",
            AB=_Block(_D12, in_cls=_BetaField, out_cls=_AlphaField,
                      domain=_SP_B, codomain=_SP_A, adjointable_=False),
        )
        if no_adj.inverse().is_adjointable:
            pytest.fail("wrap advertises a transpose over a "
                        "non-adjointable coupling block")
        with pytest.raises(MissingAdjoint, match=r"coupling block \(0, 1\)"):
            no_adj.solve_transpose(_rhs(16))
        # (b) a diagonal without the direct transpose-solve verb.
        no_reverse = _grid(
            "upper",
            BB=_NoReverseBlock(_D22, in_cls=_BetaField, out_cls=_BetaField,
                               domain=_SP_B, codomain=_SP_B),
        )
        if no_reverse.inverse().is_adjointable:
            pytest.fail("wrap advertises a transpose over a diagonal "
                        "with no solve_transpose")
        with pytest.raises(MissingAdjoint, match="no solve_transpose"):
            no_reverse.solve_transpose(_rhs(17))

    # ── S4 — the structure-keyed inverse() factory ────────────────────

    def test_inverse_is_structure_keyed(self) -> None:
        upper, q = _grid("upper"), _rhs(18)
        substitution = upper.inverse()
        if type(substitution) is not CoupledSubstitutionOperator:
            pytest.fail("triangular grid must key the substitution wrap")
        if substitution.inverse() is not upper:
            pytest.fail("the involution (A⁻¹)⁻¹ is not the grid itself")
        round_trip = upper.apply(substitution.apply(q))
        np.testing.assert_allclose(
            round_trip.to_flat(), q.to_flat(), rtol=1e-12, atol=1e-14,
        )
        np.testing.assert_array_equal(  # seeded-apply: accepted and DROPPED
            substitution.apply(q).to_flat(),
            substitution.apply(q, initial_guess=_rhs(9)).to_flat(),
        )
        full = _grid("full")
        if type(full.inverse()) is not MatrixInverseOperator:
            pytest.fail("full square grid must key the materialized LU")
        with pytest.raises(NotInvertible, match="block-triangular"):
            CoupledSubstitutionOperator(full)

    def test_full_grid_lu_solve_value_and_member_types(self) -> None:
        grid, q, template = _grid("full"), _rhs(19), _zeros_field()
        dense = _probe_dense(grid, template)
        x = grid.solve(q)  # the documented one-shot LU convenience
        np.testing.assert_allclose(
            x.to_flat(), np.linalg.solve(dense, q.to_flat()),
            rtol=1e-12, atol=1e-14,
        )
        if type(x.systems[0]) is not _AlphaField or (
            type(x.systems[1]) is not _BetaField
        ):
            pytest.fail(
                "the LU route must mint the solution from the DOMAIN's "
                "zero exemplar (typed members), never return flats"
            )
        np.testing.assert_allclose(
            grid.solve_transpose(q).to_flat(),
            np.linalg.solve(dense.T, q.to_flat()),
            rtol=1e-12, atol=1e-14,
        )

    def test_dead_diagonal_falls_through_to_the_lu_route(self) -> None:
        grid = _grid(
            "upper",
            BB=_Block(_D22, in_cls=_BetaField, out_cls=_BetaField,
                      domain=_SP_B, codomain=_SP_B, invertible_=False),
        )
        if grid._substitution_ready():
            pytest.fail("a non-invertible diagonal must disqualify the "
                        "substitution route")
        if not grid.is_invertible:
            pytest.fail("the grid is still materializable — the LU route "
                        "must carry it")
        if type(grid.inverse()) is not MatrixInverseOperator:
            pytest.fail("the fall-through must key the materialized LU")
        q, template = _rhs(20), _zeros_field()
        np.testing.assert_allclose(  # and it is CORRECT, not just typed
            grid.solve(q).to_flat(),
            np.linalg.solve(_probe_dense(grid, template), q.to_flat()),
            rtol=1e-12, atol=1e-14,
        )

    def test_singular_grid_fails_loud_at_inverse_construction(self) -> None:
        grid = _grid(
            "full",
            BA=_Block(np.zeros((2, 3)), in_cls=_AlphaField,
                      out_cls=_BetaField, domain=_SP_A, codomain=_SP_B),
            BB=_Block(np.zeros((2, 2)), in_cls=_BetaField,
                      out_cls=_BetaField, domain=_SP_B, codomain=_SP_B,
                      invertible_=False),
        )
        with pytest.raises(np.linalg.LinAlgError, match="exactly singular"):
            grid.inverse()

    # ── as_matrix — the typed-carrier Op → Mat realization ────────────

    def test_as_matrix_typed_probe_route(self) -> None:
        """A non-assemblable grid materializes by TYPED basis probing
        through the domain's zero exemplar — bit-identical to the test's
        own apply-probe (same basis drive), with the base gate contract
        (basis_shape agreement + MatrixTooLarge) preserved."""
        grid = _grid("upper", assemblable=False)
        np.testing.assert_array_equal(
            grid.as_matrix(), _probe_dense(grid, _zeros_field()),
        )
        with pytest.raises(ValueError, match="contradicts"):
            grid.as_matrix(basis_shape=(7,))
        with pytest.raises(MatrixTooLarge):
            grid.as_matrix(max_dimension=2)
        zeroless = _grid("upper", zeros=False, assemblable=False)
        with pytest.raises(RuntimeError, match="zero-element factory"):
            zeroless.as_matrix()

    def test_as_matrix_assemblable_grid_routes_through_assembly(self) -> None:
        grid = _grid("upper")
        np.testing.assert_allclose(  # principled-equiv — L16 scatter order
            grid.as_matrix(), _probe_dense(grid, _zeros_field()),
            rtol=1e-11,
        )


class TestCoupledSpaceZeros:
    """The zero-element seam (:meth:`CoupledSpace.zeros`) — the
    typed-carrier materialization's template source."""

    def test_unwired_zeros_raises_the_seam_message(self) -> None:
        with pytest.raises(RuntimeError, match="zero-element factory"):
            _solvable_space(zeros=False).zeros()

    def test_wired_zeros_mints_fresh_owned_buffers(self) -> None:
        space = _solvable_space()
        first, second = space.zeros(), space.zeros()
        if first is second:
            pytest.fail("zeros() must mint a FRESH field per call")
        np.asarray(first.systems[0].values)[...] = 7.0
        np.testing.assert_array_equal(
            second.systems[0].values, np.zeros(3),
            err_msg="zeros() results share a buffer — factory semantics "
                    "broken",
        )

    def test_zeros_arity_is_checked(self) -> None:
        lopsided = CoupledSpace.from_systems(
            (_SP_A, _SP_B),
            zeros=lambda: CoupledField(
                systems=(_AlphaField(values=np.zeros(3)),),
            ),
        )
        with pytest.raises(ValueError, match="pairing is inconsistent"):
            lopsided.zeros()
