r"""The assembly functor downstream of the operator algebra (numerics_untangling S3).

Rulings 24, 60, 61, 62 of ``.claude/plans/numerics_untangling.md``; the
verification spec is ``scratch/numerics_untangling/packages/S3_gate_spec.md``.
Assembly is a homomorphism out of the free algebra:

* the functor :func:`~orpheus.numerics.assembled_operator.assemble` has one
  arm per composite of the free algebra (a combination folds left in
  insertion order, a word multiplies, a fused atom assembles its
  expansion), returns an assembled operator itself, and otherwise asks the
  operator's own ``assembled_matrix(emit)``;
* a leaf's emission is a RAW sparse matrix and only the functor wraps it;
  a domain composite (the coupled block law, the lift's embedding, the
  flattened pass-through) recurses through the ``emit`` it is handed;
* ``operator.py`` knows the flag ``is_assemblable`` and nothing else of
  assembly;
* ``as_matrix`` is the probe and ``dense`` the assembly route: two
  algorithms for one matrix, and ``MatrixInverseOperator`` factors
  ``dense``.

Fixtures: products on integer-valued leaves (exact under any summation
order); sums on a seeded non-dyadic fixture whose left fold differs from
the other associations, so the fold ORDER is visible (the dyadic law rows
of ``test_assembled_operator.py`` cannot see it).
"""

from __future__ import annotations

from orpheus.transport.fields.cross_section_field import CrossSectionField

import ast
import pathlib
import re
from typing import Any, Optional

import numpy as np
import pytest
from scipy import sparse

import orpheus.numerics.assembled_operator as asm
import orpheus.numerics.matrix_inverse_operator as mio
from orpheus.numerics.assembled_operator import (
    MissingAssembly,
    SparseAssembledOperator,
    assemble,
    dense,
)
from orpheus.numerics.elementary_operators import DiagonalOperator
from orpheus.numerics.flat_operator import FlattenedOperator
from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator
from orpheus.numerics.operator import (
    FusedAtom,
    IdentityOperator,
    LinearOperator,
    MatrixTooLarge,
    OperatorSum,
)
from orpheus.numerics.space import FunctionSpace

pytestmark = [pytest.mark.foundation]

ORPHEUS = pathlib.Path(asm.__file__).resolve().parents[1]
SUM_LAW = "[Σ cᵢwᵢ] = Σ cᵢ[wᵢ]"
WORD_LAW = "[A₁…Aₙ] = [A₁]…[Aₙ]"


# ── Fixtures ───────────────────────────────────────────────────────────


def _leaf(dense_matrix: np.ndarray) -> SparseAssembledOperator:
    return SparseAssembledOperator(sparse.csr_array(dense_matrix))


def _non_dyadic() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Three seeded 4×4 matrices, one structural zero each (`[M]` S3 spec
    finding 5: the left fold of A − B + 0.3C differs from the right
    association and from the reversed fold at 2 entries)."""
    rng = np.random.default_rng(2026)
    a, b, c = (rng.random((4, 4)) for _ in range(3))
    a[0, 1] = 0.0
    b[2, 3] = 0.0
    c[1, 1] = 0.0
    return a, b, c


def _integer(seed: int) -> np.ndarray:
    return np.random.default_rng(seed).integers(-3, 4, (4, 4)).astype(float)


class _Leaf(LinearOperator):
    """A leaf that knows its stencil: emits a raw matrix and records the
    ``emit`` it was handed."""

    def __init__(self, matrix: np.ndarray) -> None:
        self.matrix = matrix
        self.received: list[object] = []

    @property
    def domain(self) -> Optional[FunctionSpace]:
        return None

    @property
    def codomain(self) -> Optional[FunctionSpace]:
        return None

    def apply(self, x, /):
        return self.matrix @ np.asarray(x, dtype=float).ravel()

    @property
    def is_assemblable(self) -> bool:
        return True

    def assembled_matrix(self, emit, /) -> Any:
        self.received.append(emit)
        return sparse.csr_array(self.matrix)


class _Atom(FusedAtom):
    """A fused atom built as production atoms are: an expansion, an action
    of its own, the expansion's roles adopted at construction."""

    def __init__(self, expansion: LinearOperator) -> None:
        self._expansion = expansion
        self._adopt_expansion_roles()

    def expansion(self) -> LinearOperator:
        return self._expansion

    def apply(self, x, /):
        return self._expansion.apply(x)


# ── G1: the functor's arms ─────────────────────────────────────────────


class TestTheArms:
    @pytest.mark.verifies("matrix-functor-homomorphism")
    def test_g1a_the_sum_arm_is_the_left_fold_in_insertion_order(self) -> None:
        a, b, c = _non_dyadic()
        op = _leaf(a) - _leaf(b) + 0.3 * _leaf(c)
        assert isinstance(op, OperatorSum)
        assert [coef for coef, _ in op.terms] == [1.0, -1.0, 0.3], "premise: three terms in order"
        left = (a - b) + 0.3 * c
        # Premise: the fixture makes the fold order visible.
        assert not np.array_equal(left, a + (-b + 0.3 * c)), "premise: right association coincides"
        assert not np.array_equal(left, (0.3 * c - b) + a), "premise: reversed fold coincides"
        np.testing.assert_array_equal(assemble(op).matrix.toarray(), left)

    def test_g1g_a_leaf_emission_is_summed_before_its_coefficient(self) -> None:
        """The leaf arm canonicalises each emission to CSR (duplicate
        ``(row, col)`` entries summed, the FEM scatter) BEFORE the sum's
        coefficient touches it: ``[c·L] = c·(v₁ + v₂)``, the old carrier's
        rounding. Scaling the raw COO first gives ``c·v₁ + c·v₂``, which
        differs in the last bit on this seeded pair (premise asserted)."""
        rng = np.random.default_rng(31)
        pairs = rng.random((3, 2))
        v1, v2 = pairs[2]
        c = 0.3
        assert c * (v1 + v2) != c * v1 + c * v2, "premise: the two roundings coincide"

        class _ScatterLeaf(_Leaf):
            def assembled_matrix(self, emit, /) -> Any:
                return sparse.coo_array(
                    (np.array([v1, v2, 1.0]), (np.array([0, 0, 1]), np.array([0, 0, 1]))),
                    shape=(2, 2),
                )

        got = assemble(c * _ScatterLeaf(np.zeros((2, 2)))).matrix.toarray()
        assert got[0, 0] == c * (v1 + v2), (got[0, 0], c * (v1 + v2), c * v1 + c * v2)

    @pytest.mark.verifies("matrix-functor-homomorphism")
    def test_g1b_the_product_arm_multiplies_in_word_order(self) -> None:
        p, q, r = _integer(1), _integer(2), _integer(3)
        assert not np.array_equal(p @ q @ r, r @ q @ p), "premise: the word is order-sensitive"
        np.testing.assert_array_equal(
            assemble(_leaf(p) @ _leaf(q) @ _leaf(r)).matrix.toarray(), p @ q @ r,
        )

    def test_g1c_the_fused_arm_assembles_the_expansion(self) -> None:
        a, b = _integer(4), _integer(5)
        atom = _Atom(_leaf(a) - 2.0 * _leaf(b))
        assert atom.is_assemblable
        np.testing.assert_array_equal(assemble(atom).matrix.toarray(), a - 2.0 * b)
        silent = _Atom(_leaf(a) + IdentityOperator())
        assert not silent.is_assemblable
        with pytest.raises(MissingAssembly, match="IdentityOperator emits no assembly"):
            assemble(silent)
        bare = _Atom(IdentityOperator())
        with pytest.raises(
            MissingAssembly,
            match=re.escape("[_Atom] = [its expansion]: IdentityOperator emits no assembly"),
        ):
            assemble(bare)

    def test_g1h_the_fused_arm_builds_the_expansion_once(self) -> None:
        """The arm reads the atom's cached ``expanded``; repeated assembly
        never rebuilds the expansion (S3 review)."""
        a, b = _integer(4), _integer(5)
        built: list[int] = []

        class _Counting(_Atom):
            def expansion(self) -> LinearOperator:
                built.append(1)
                return super().expansion()

        atom = _Counting(_leaf(a) - 2.0 * _leaf(b))
        for _ in range(3):
            np.testing.assert_array_equal(assemble(atom).matrix.toarray(), a - 2.0 * b)
            dense(atom, basis_shape=(4,))
        assert len(built) == 1, f"expansion() built {len(built)} times"

    def test_g1d_an_assembled_operator_is_its_own_assembly(self) -> None:
        s = _leaf(_integer(6))
        assert assemble(s) is s

    def test_g1e_a_leaf_emits_raw_and_the_functor_wraps_once_with_its_ends(
        self, monkeypatch,
    ) -> None:
        space = FunctionSpace(name="g1e", shape=(4,))
        m = _integer(7)

        class _BoundLeaf(_Leaf):
            @property
            def domain(self) -> Optional[FunctionSpace]:
                return space

            @property
            def codomain(self) -> Optional[FunctionSpace]:
                return space

        leaf = _BoundLeaf(m)
        got = assemble(leaf)
        assert type(got) is SparseAssembledOperator
        np.testing.assert_array_equal(got.matrix.toarray(), m)
        assert got.domain is space and got.codomain is space
        assert len(leaf.received) == 1 and callable(leaf.received[0])
        # The handed emitter is the functor's raw recursion: a matrix, not
        # a carrier.
        raw = leaf.received[0](_leaf(m))
        assert not isinstance(raw, LinearOperator)
        np.testing.assert_array_equal(sparse.csr_array(raw).toarray(), m)

    @pytest.mark.verifies("matrix-functor-homomorphism")
    def test_g1f_nesting_recurses_through_every_arm(self) -> None:
        a, b, c, d = (_integer(s) for s in (8, 9, 10, 11))
        op = _leaf(a) @ (_leaf(b) + 2.0 * _leaf(c)) - _leaf(d)
        np.testing.assert_array_equal(
            assemble(op).matrix.toarray(), a @ (b + 2.0 * c) - d,
        )


# ── G2: the domain composites recurse through the emit they are handed ──


def _recording_emit():
    calls: list[object] = []

    def emit(part):
        calls.append(part)
        return asm._emit(part)

    return emit, calls


class TestDomainCompositesRecurseThroughTheArgument:
    def test_g2_coupled_block_law(self) -> None:
        from tests.gates.numerics.test_coupled_operator import _grid, _probe_dense, _zeros_field

        grid = _grid("full")
        emit, calls = _recording_emit()
        grid.assembled_matrix(emit)
        present = [b for row in grid.blocks for b in row if b is not None]
        assert len(calls) == len(present) == 4
        assert all(any(c is b for c in calls) for b in present)
        # A composite block (a combination) is reached through the functor.
        from tests.gates.numerics.test_coupled_operator import _solve_blocks

        blocks = _solve_blocks()
        composite = grid.__class__(
            [[blocks["AA"] + blocks["AA"], blocks["AB"]], [blocks["BA"], blocks["BB"]]],
            domain=grid.domain, codomain=grid.codomain,
        )
        np.testing.assert_allclose(
            dense(composite), _probe_dense(composite, _zeros_field()), rtol=1e-11,
        )

    def test_g2_flattened_pass_through(self) -> None:
        inner = _leaf(_integer(12))

        class _Duck:
            def to_flat(self):
                return np.zeros(4)

        emit, calls = _recording_emit()
        FlattenedOperator(inner, _Duck()).assembled_matrix(emit)
        assert len(calls) == 1 and calls[0] is inner

    def test_g2_bulk_lift_embedding(self) -> None:
        from tests.gates.transport.test_bulk_lift import _ShapedDiagonal
        from orpheus.transport.operators.lift import BulkLift

        from tests.gates.diffusion.test_operators import _config_setup

        mesh, _, _, _ = _config_setup("zero_flux/zero_flux")
        inner = _ShapedDiagonal(
            np.random.default_rng(13).random(mesh.bulk_space.shape) + 0.5,
            domain=mesh.bulk_space, codomain=mesh.bulk_space,
        )
        space = mesh.full_field_space
        emit, calls = _recording_emit()
        BulkLift(inner, domain=space, codomain=space).assembled_matrix(emit)
        assert len(calls) == 1 and calls[0] is inner


# ── G3: the refusal names the part and the law ─────────────────────────


class TestRefusal:
    def test_g3_the_deepest_part_and_the_law_kind(self) -> None:
        a, b = _leaf(_integer(14)), _leaf(_integer(15))
        with pytest.raises(MissingAssembly) as nested:
            assemble(a + b @ IdentityOperator())
        assert re.search(re.escape(WORD_LAW) + ": IdentityOperator emits no assembly", str(nested.value))
        assert SUM_LAW not in str(nested.value)
        with pytest.raises(MissingAssembly) as summed:
            assemble(a + IdentityOperator())
        assert re.search(re.escape(SUM_LAW) + ": IdentityOperator emits no assembly", str(summed.value))
        assert WORD_LAW not in str(summed.value)
        with pytest.raises(MissingAssembly, match="assemble: IdentityOperator emits no assembly"):
            assemble(IdentityOperator())
        assert issubclass(MissingAssembly, TypeError)

    def test_g3_a_domain_composite_names_its_non_emitting_part(self) -> None:
        """A domain composite's flag is its parts'; the functor runs the
        composite's own recursion, so the refusal names the part and the
        composite's law, never only the composite (qa F1). Disjointness: the
        top-level fragment ``assemble: <Composite>`` is absent."""
        from tests.gates.numerics.test_coupled_operator import _grid

        grid = _grid("full", assemblable=False)
        with pytest.raises(MissingAssembly) as coupled:
            assemble(grid)
        assert re.search(r"CoupledOperator's assembly: block \(0, 0\) \(_Block\)", str(coupled.value))
        assert "assemble: CoupledOperator" not in str(coupled.value)

        class _Duck:
            def to_flat(self):
                return np.zeros(4)

        with pytest.raises(MissingAssembly) as flat:
            assemble(FlattenedOperator(_leaf(_integer(16)) + IdentityOperator(), _Duck()))
        assert re.search(re.escape(SUM_LAW) + ": IdentityOperator emits no assembly", str(flat.value))
        assert "assemble: FlattenedOperator" not in str(flat.value)
        with pytest.raises(MissingAssembly) as bare:
            assemble(FlattenedOperator(IdentityOperator(), _Duck()))
        assert "FlattenedOperator's assembly: IdentityOperator emits no assembly" in str(bare.value)

    def test_g3_the_size_refusal_names_the_densification(self) -> None:
        a, b, c = _non_dyadic()
        op = _leaf(a) - _leaf(b) + 0.3 * _leaf(c)
        for densify in (lambda: op.as_matrix(basis_shape=(4,), max_dimension=2),
                        lambda: dense(op, basis_shape=(4,), max_dimension=2)):
            with pytest.raises(MatrixTooLarge, match="densifying OperatorSum: basis dimension 4"):
                densify()


# ── G4: who constructs the carrier, and what operator.py knows ─────────


def _emission_bodies(tree: ast.AST):
    for cls in ast.walk(tree):
        if isinstance(cls, ast.ClassDef):
            for fn in cls.body:
                if isinstance(fn, ast.FunctionDef) and fn.name == "assembled_matrix":
                    is_stub = all(  # a Protocol's declaration, not an emission body
                        isinstance(s, ast.Expr) and isinstance(s.value, ast.Constant)
                        for s in fn.body
                    )
                    if not is_stub:
                        yield cls.name, fn


def _constructs_the_carrier(fn: ast.FunctionDef) -> bool:
    return any(
        isinstance(n, ast.Call)
        and (getattr(n.func, "id", None) or getattr(n.func, "attr", None)) == "SparseAssembledOperator"
        for n in ast.walk(fn)
    )


class TestWhoConstructs:
    def test_g4a_no_emission_body_constructs_the_carrier(self) -> None:
        definers, offenders = set(), []
        for path in sorted(ORPHEUS.rglob("*.py")):
            for cls, fn in _emission_bodies(ast.parse(path.read_text(encoding="utf-8"))):
                definers.add(cls)
                if _constructs_the_carrier(fn):
                    offenders.append(f"{path.relative_to(ORPHEUS)}::{cls}")
        assert definers == {
            "SparseAssembledOperator", "CoupledOperator", "FlattenedOperator", "BulkLift",
            "MultiplicationOperator", "IsotropicTransfer", "LeakageOperator",
            "DiffusionBoundaryOperator",
        }, f"the emission population moved: {sorted(definers)}"
        assert offenders == [], f"emission bodies that wrap their own matrix: {offenders}"
        planted = ast.parse(
            "class X:\n    def assembled_matrix(self, emit, /):\n"
            "        return SparseAssembledOperator(m)\n"
        )
        assert [c for c, f in _emission_bodies(planted) if _constructs_the_carrier(f)] == ["X"]

    @pytest.mark.parametrize("subject", ["diffusion_loss", "homogeneous_loss"])
    def test_g4b_one_wrap_per_functor_call(self, subject: str, monkeypatch) -> None:
        """Transitively: a helper an emitter calls (``embed_bulk_assembly``,
        ``_bulk_assembly``) wrapping its result is invisible to G4a."""
        if subject == "diffusion_loss":
            from tests.gates.diffusion.test_operators import _config_setup, _loss

            mesh, mat_xs, template, _ = _config_setup("albedo(0.3)/reflective")
            op: LinearOperator = FlattenedOperator(_loss(mesh, mat_xs), template)
        else:
            from orpheus.derivations import get
            from orpheus.homogeneous.solver import HomogeneousProblem

            op = HomogeneousProblem(next(iter(get("homo_2eg_n2n").materials.values()))).loss
        wraps: list[int] = []
        original = SparseAssembledOperator.__init__

        def counting(self, *args, **kwargs):
            wraps.append(1)
            original(self, *args, **kwargs)

        monkeypatch.setattr(SparseAssembledOperator, "__init__", counting)
        assemble(op)
        assert len(wraps) == 1, f"{len(wraps)} carriers built for one assembly"

    def test_g4c_operator_py_knows_only_the_flag(self) -> None:
        names = {"SupportsAssembly", "assemblable", "MissingAssembly",
                 "SparseAssembledOperator", "assemble", "dense",
                 "_as_matrix_by_probing", "assembled_matrix"}

        def scan(path: pathlib.Path) -> tuple[set[str], set[str], set[str]]:
            tree = ast.parse(path.read_text(encoding="utf-8"))
            imports = {
                n.module or "" for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)
            } | {a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
            defined = {
                n.name for n in ast.walk(tree)
                if isinstance(n, (ast.FunctionDef, ast.ClassDef))
            }
            loaded = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)} | {
                n.attr for n in ast.walk(tree) if isinstance(n, ast.Attribute)
            }
            return imports, defined, loaded

        imports, defined, loaded = scan(ORPHEUS / "numerics" / "operator.py")
        assert not any("assembled_operator" in m for m in imports)
        assert not (defined | loaded) & names, sorted((defined | loaded) & names)
        assert "is_assemblable" in defined
        # Positive control: the same scan finds the functor downstream.
        _, defined_asm, _ = scan(ORPHEUS / "numerics" / "assembled_operator.py")
        assert {"assemble", "dense", "MissingAssembly", "assemblable"} <= defined_asm


# ── G5-G7: the two densifications ──────────────────────────────────────


class TestTwoAlgorithms:
    def test_g5_dense_routes_by_the_flag(self, monkeypatch) -> None:
        a, b, c = _non_dyadic()
        op = _leaf(a) - _leaf(b) + 0.3 * _leaf(c)
        functor_calls: list[object] = []
        probe_calls: list[object] = []
        original_assemble, original_probe = asm.assemble, LinearOperator.as_matrix

        def spy_assemble(x):
            functor_calls.append(x)
            return original_assemble(x)

        def spy_probe(self, *args, **kwargs):
            probe_calls.append(self)
            return original_probe(self, *args, **kwargs)

        monkeypatch.setattr(asm, "assemble", spy_assemble)
        monkeypatch.setattr(LinearOperator, "as_matrix", spy_probe)
        got = dense(op, basis_shape=(4,))
        assert functor_calls == [op] and probe_calls == []
        np.testing.assert_array_equal(got, original_assemble(op).matrix.toarray())

        functor_calls.clear()
        diagonal = DiagonalOperator(np.array([1.0, 2.0, 3.0]))
        np.testing.assert_array_equal(dense(diagonal, basis_shape=(3,)), np.diag([1.0, 2.0, 3.0]))
        assert functor_calls == [] and probe_calls == [diagonal]

        for subject, n in ((op, 4), (diagonal, 3)):
            with pytest.raises(MatrixTooLarge):
                dense(subject, basis_shape=(n,), max_dimension=2)
        with pytest.raises(ValueError):
            dense(op, basis_shape=(5,))

    def test_g6_as_matrix_is_the_probe(self, monkeypatch) -> None:
        a, b, c = _non_dyadic()
        leaves = [_leaf(a), _leaf(b), _leaf(c)]
        op = leaves[0] - leaves[1] + 0.3 * leaves[2]
        functor_calls: list[object] = []
        monkeypatch.setattr(asm, "_emit", lambda *x: functor_calls.append(x))
        applies: list[int] = []
        original_apply = SparseAssembledOperator.apply

        def counting_apply(self, x, /):
            applies.append(1)
            return original_apply(self, x)

        monkeypatch.setattr(SparseAssembledOperator, "apply", counting_apply)
        got = op.as_matrix(basis_shape=(4,))
        assert functor_calls == []
        assert len(applies) == 3 * 4, "one apply per leaf per basis vector"
        e2 = np.zeros(4)
        e2[2] = 1.0
        np.testing.assert_array_equal(got[:, 2], op.apply(e2))

    def test_g7a_probe_and_dense_two_algorithms_one_matrix(self) -> None:
        """REFERENCE row. Shared upstream, stated: both read the leaves'
        stored entries, so a wrong stored entry is outside this row (G1a's
        numpy oracle sees it); a wrong functor arm is inside it."""
        a, b, c = _non_dyadic()
        op = _leaf(a) - _leaf(b) + 0.3 * _leaf(c)
        np.testing.assert_array_equal(op.as_matrix(basis_shape=(4,)), dense(op, basis_shape=(4,)))

    def test_g7b_probe_and_dense_agree_on_every_homogeneous_loss(self) -> None:
        """`[M]` 8 of 8 mixtures at S3: the k∞ bound reads ``dense`` (the
        matrix the solver factors); this row keeps the probe and the
        assembly one matrix on the same population."""
        from tests.gates.homogeneous._homogeneous_population import mixture_cases
        from orpheus.homogeneous.solver import HomogeneousProblem

        cases = mixture_cases()
        assert len(cases) == 8, sorted(cases)
        for name, mix in cases.items():
            loss = HomogeneousProblem(mix).loss
            assert loss.is_assemblable, name
            np.testing.assert_array_equal(loss.as_matrix(), dense(loss), err_msg=name)


# ── G8-G9: MatrixInverseOperator factors dense; the coupled grid probes ──


def _lu_inputs(monkeypatch) -> list[np.ndarray]:
    seen: list[np.ndarray] = []
    original = mio.lu_factor

    def spy(a, *args, **kwargs):
        seen.append(np.array(a, copy=True))
        return original(a, *args, **kwargs)

    monkeypatch.setattr(mio, "lu_factor", spy)
    return seen


def _functor_calls(monkeypatch) -> list[LinearOperator]:
    seen: list[LinearOperator] = []
    original = asm.assemble

    def spy(op):
        seen.append(op)
        return original(op)

    monkeypatch.setattr(asm, "assemble", spy)
    return seen


class TestMatrixInverseFactorsDense:
    def test_g8a_the_factored_matrix_is_dense_on_a_route_discriminating_loss(
        self, monkeypatch,
    ) -> None:
        """The diffusion loss with its collision from Σ_a: a legitimate
        operator whose probe and assembly differ in the last bit (`[M]` S3
        spec finding 3: entries (2, 2), (6, 6)), so this VALUE row sees
        which route ``MatrixInverseOperator`` took."""
        from tests.gates.diffusion.test_operators import _config_setup
        from orpheus.diffusion.operators import DiffusionBoundaryOperator, LeakageOperator
        from orpheus.transport.operators.isotropic_transfer import IsotropicScattering
        from orpheus.transport.operators.lift import BulkLift
        from orpheus.transport.operators.multiplication_operator import MultiplicationOperator

        mesh, mat_xs, template, _ = _config_setup("zero_flux/zero_flux")
        ffs = mesh.full_field_space
        loss = (
            LeakageOperator(mesh)
            + MultiplicationOperator(CrossSectionField.absorption_of(mat_xs), domain=ffs, codomain=ffs)
            - BulkLift(IsotropicScattering.from_material_xs(mat_xs, space=mesh.bulk_space), domain=ffs, codomain=ffs)
            - DiffusionBoundaryOperator(mesh)
        )
        flat = FlattenedOperator(loss, template)
        assembled = dense(flat)
        assert not np.array_equal(flat.as_matrix(), assembled), (
            "premise lost: probe and assembly agree bitwise, the row cannot see the route"
        )
        seen = _lu_inputs(monkeypatch)
        MatrixInverseOperator(flat)
        assert len(seen) == 1
        np.testing.assert_array_equal(seen[0], assembled)

    def test_g8b_the_diffusion_resolvent_assembles(self, monkeypatch) -> None:
        from tests.gates.diffusion.test_operators import _config_setup, _loss

        mesh, mat_xs, template, _ = _config_setup("reflective/marshak_vacuum")
        flat = FlattenedOperator(_loss(mesh, mat_xs), template)
        calls = _functor_calls(monkeypatch)
        MatrixInverseOperator(flat)
        assert [op for op in calls if op is flat], calls

    def test_g8b_the_homogeneous_solver_assembles(self, monkeypatch) -> None:
        from orpheus.derivations import get
        from orpheus.homogeneous.solver import solve_homogeneous_infinite

        calls = _functor_calls(monkeypatch)
        seen = _lu_inputs(monkeypatch)
        solve_homogeneous_infinite(next(iter(get("homo_2eg_n2n").materials.values())))
        assert calls, "the homogeneous solver never reached the functor"
        assert len(seen) >= 1
        reached = list(calls)
        assert any(np.array_equal(s, dense(op)) for s in seen for op in reached), (
            "no factored matrix is the assembly of an operator the solver assembled"
        )

    def test_g8b_the_coupled_inverse_factors_dense(self, monkeypatch) -> None:
        from tests.gates.numerics.test_coupled_operator import _grid

        grid = _grid("full")
        assert not grid._substitution_ready(), "premise: the dense route is taken"
        expected = dense(grid)
        calls = _functor_calls(monkeypatch)
        seen = _lu_inputs(monkeypatch)
        inverse = grid.inverse()
        assert isinstance(inverse, MatrixInverseOperator)
        assert [op for op in calls if op is grid], calls
        assert len(seen) == 1
        np.testing.assert_array_equal(seen[0], expected)


class TestCoupledAsMatrixProbes:
    def test_g9_an_assemblable_grid_probes(self, monkeypatch) -> None:
        """ROUTE row: ``as_matrix`` is the typed probe (bit-identical to the
        test's own drive) and never calls the functor; ``dense`` is the
        block scatter, equal at the L16 scatter-order band."""
        from tests.gates.numerics.test_coupled_operator import _grid, _probe_dense, _zeros_field

        grid = _grid("full")
        assert grid.is_assemblable
        calls = _functor_calls(monkeypatch)
        probed = grid.as_matrix()
        assert calls == []
        np.testing.assert_array_equal(probed, _probe_dense(grid, _zeros_field()))
        np.testing.assert_allclose(dense(grid), probed, rtol=1e-11)


# ── G10: the flag and the functor agree (the witness of ruling F2) ─────


def _registry() -> dict[str, "Any"]:
    """Every operator kind the S3 fixtures build, keyed by name and built
    lazily (a rebinding battery arm must reach the rows: lessons L100)."""

    def diffusion():
        from tests.gates.diffusion.test_operators import _config_setup

        return _config_setup("albedo(0.3)/reflective")

    def diffusion_parts():
        from orpheus.diffusion.operators import DiffusionBoundaryOperator, LeakageOperator
        from orpheus.transport.operators.isotropic_transfer import IsotropicN2N, IsotropicScattering
        from orpheus.transport.operators.lift import BulkLift
        from orpheus.transport.operators.multiplication_operator import MultiplicationOperator

        mesh, mat_xs, template, _ = diffusion()
        ffs = mesh.full_field_space
        S = IsotropicScattering.from_material_xs(mat_xs, space=mesh.bulk_space)
        return {
            "L": LeakageOperator(mesh),
            "B": DiffusionBoundaryOperator(mesh),
            "C": MultiplicationOperator(CrossSectionField.total_of(mat_xs), domain=ffs, codomain=ffs),
            "S_plain": S,
            "N2N_plain": IsotropicN2N.from_material_xs(mat_xs, space=mesh.bulk_space),
            "S_lift": BulkLift(S, domain=ffs, codomain=ffs),
            "template": template,
            "mesh": mesh,
        }

    def homogeneous():
        from orpheus.derivations import get
        from orpheus.homogeneous.solver import HomogeneousProblem

        return HomogeneousProblem(next(iter(get("homo_2eg_n2n").materials.values())))

    def coupled(**kw):
        from tests.gates.numerics.test_coupled_operator import _grid

        return _grid("full", **kw)

    class _Duck:
        def to_flat(self):
            return np.zeros(4)

    a, b = _leaf(_integer(20)), _leaf(_integer(21))
    ident = IdentityOperator()
    diag = DiagonalOperator(np.arange(1.0, 5.0))
    builders = {
        "carrier": lambda: a,
        "sum": lambda: a - 0.3 * b,
        "word": lambda: a @ b,
        "nested": lambda: a @ (b + 2.0 * a) - b,
        "atom": lambda: _Atom(a - 2.0 * b),
        "sum_with_identity": lambda: a + ident,
        "word_with_diagonal": lambda: a @ diag,
        "atom_over_identity": lambda: _Atom(a + ident),
        "identity": lambda: ident,
        "diagonal": lambda: diag,
        "dual_of_carrier": lambda: a.dual(),
        "inverse_of_diagonal": lambda: diag.inverse(),
        "flattened": lambda: FlattenedOperator(a + b, _Duck()),
        "flattened_identity": lambda: FlattenedOperator(ident, _Duck()),
        "coupled": lambda: coupled(),
        "coupled_silent_block": lambda: coupled(assemblable=False),
        "diffusion_L": lambda: diffusion_parts()["L"],
        "diffusion_B": lambda: diffusion_parts()["B"],
        "diffusion_C": lambda: diffusion_parts()["C"],
        "diffusion_S_plain": lambda: diffusion_parts()["S_plain"],
        "diffusion_N2N_plain": lambda: diffusion_parts()["N2N_plain"],
        "diffusion_S_lift": lambda: diffusion_parts()["S_lift"],
        "diffusion_loss_flat": lambda: (lambda p: FlattenedOperator(
            p["L"] + p["C"] - p["S_lift"] - p["B"], p["template"]))(diffusion_parts()),
        "homogeneous_loss": lambda: homogeneous().loss,
        "homogeneous_production": lambda: homogeneous().production,
    }
    return builders


_EXPECTED_ASSEMBLABLE = {
    "carrier", "sum", "word", "nested", "atom", "flattened", "coupled",
    "diffusion_L", "diffusion_B", "diffusion_C", "diffusion_S_plain",
    "diffusion_N2N_plain", "diffusion_S_lift", "diffusion_loss_flat",
    "homogeneous_loss",
}


class TestTheFlagIsTheFunctorsDomain:
    def test_g10_the_registry_has_both_kinds(self) -> None:
        names = set(_registry())
        assert len(names) == 25
        assert _EXPECTED_ASSEMBLABLE < names
        assert len(names - _EXPECTED_ASSEMBLABLE) == 10

    @pytest.mark.parametrize("name", sorted(_registry()))
    def test_g10_is_assemblable_iff_assemble_succeeds(self, name: str) -> None:
        """The composites' flags live in ``operator.py`` beside the other
        capability axes; the functor lives downstream (user ruling on
        elegance F2). This row is the witness that keeps the two equal: the
        flag says True exactly when ``assemble`` returns, False exactly when
        it refuses with :class:`MissingAssembly`."""
        op = _registry()[name]()
        flag = op.is_assemblable
        assert flag == (name in _EXPECTED_ASSEMBLABLE), f"{name}: flag {flag}"
        try:
            assemble(op)
        except MissingAssembly:
            succeeded = False
        else:
            succeeded = True
        assert flag == succeeded, f"{name}: is_assemblable {flag}, assemble succeeded {succeeded}"
