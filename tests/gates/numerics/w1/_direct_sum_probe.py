r"""Constructors the N05 gates share (spec 7.6). Not a test module: every function RETURNS data and none asserts.

``python -O`` strips a bare ``assert`` of a helper, so the collected test files assert. Every name that exists only after the unit
(``DirectSumSpace``, ``bulk_and_trace``, the label constants) is read from its module AT CALL TIME, inside each function, so that on a tree without
it a gate row fails in its body (a behavioural red a reader can interpret) and not at collection.

The summands are tensor spaces of the N04 helper ``_tensor_probe`` (a non-uniform diagonal form per slot, so a coupled summand pair is visible in the
last bit); the hand rule of a summand's metric is written here from the stated weights and reads no form object of production.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from tests.gates.numerics.w1 import _tensor_probe as tp

ROOT = Path(__file__).resolve().parents[4]

#: The weights of the three summands the gates build (``A`` and ``C`` have one shape and two forms; ``B`` has another shape).
A_ARGS = (tp.ANGULAR_W, 2, tp.SPATIAL_W)
B_ARGS = ((0.9, 1.4, 0.6), 3, (0.5, 1.9))
C_ARGS = ((0.3, 1.7, 0.9, 1.1), 2, tp.SPATIAL_W)


def direct_sum() -> Any:
    from orpheus.numerics.direct_sum import DirectSumSpace

    return DirectSumSpace


def space_of(args: tuple[Any, ...]) -> Any:
    angular, ng, spatial = args
    return tp.tensor_space().of_axes(*tp.triple_factors(angular, ng, spatial))


def hand_diagonal(args: tuple[Any, ...]) -> np.ndarray:
    """The metric diagonal of a triple space written from its stated weights: angular times spatial, one on the formless energy slot, in slot order, flat C-order."""
    angular, ng, spatial = args
    wa, ws = np.asarray(angular), np.asarray(spatial)
    return (wa[:, None, None] * np.ones((1, ng, 1)) * ws[None, None, :]).reshape(-1)


def signed_space() -> Any:
    """One formed slot with the signed diagonal form ``[1, -1, 2]``: ``negative`` is 1."""
    return tp.tensor_space().of_axes(tp.density_factor(tp.diag_axis((1.0, -1.0, 2.0), "s")))


def positive_space() -> Any:
    return tp.tensor_space().of_axes(tp.density_factor(tp.diag_axis((1.0, 2.0, 3.0), "p")))


@dataclass(frozen=True)
class ToyCarrier:
    """A typed element of a sum built in a gate: the tuple of its summand elements in label order (a ``Member`` per tensor summand, a ``ToyCarrier`` per nested sum).

    It implements the carriers' structural protocol of ON05-1 and ruling V1 and nothing else: ``blocks``, ``labels`` (the summand labels the element is
    admitted by) and ``with_blocks``. It names no space (the sum asks no membership of its operand, issue #666).
    """

    blocks: tuple[Any, ...]
    labels: tuple[str, ...]

    def with_blocks(self, blocks: tuple[Any, ...]) -> "ToyCarrier":
        return ToyCarrier(tuple(blocks), self.labels)


def member_of(space: Any, values: np.ndarray) -> Any:
    from orpheus.numerics.space import Member

    return Member(values=np.asarray(values, dtype=float).reshape(space.shape), space=space)


def toy_element(space: Any, seed: int) -> ToyCarrier:
    """A ``ToyCarrier`` of seeded members on the sum ``space``, labelled by its labels, one block per summand (a summand that is itself a sum makes a nested ``ToyCarrier``)."""
    rng = np.random.default_rng(seed)
    blocks = []
    for _, s in space.summands:
        if hasattr(s, "summands"):
            blocks.append(toy_element(s, int(rng.integers(1 << 30))))
        else:
            blocks.append(member_of(s, rng.standard_normal(s.shape)))
    return ToyCarrier(tuple(blocks), tuple(space.labels))


def flat_of(x: Any) -> np.ndarray:
    """The flat coordinates of a ``ToyCarrier`` (nested ones flattened in label order) by a hand walk."""
    parts = [flat_of(b) if isinstance(b, ToyCarrier) else np.asarray(b.values).reshape(-1) for b in x.blocks]
    return np.concatenate(parts)


def modal_head_and_tail() -> tuple[Any, Any]:
    """``(head, tail)``: a modal angular head (``frame_faces.basis_space`` of the folded product (2, 4) frame at order 2) and a 3-point function-density tensor tail."""
    from orpheus.numerics import frame_faces as ff
    from orpheus.numerics.factor import Kind, Variance, VarianceChart
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.numerics.space import TensorSpace

    head = ff.basis_space(ff.angular_frame(Quadrature.folded_product(n_mu=2, n_phi=4), 2))
    tail = TensorSpace.declared("t", (3,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    return head, tail


# ─────────────────────────── ast predicates over orpheus/ ───────────────────────────


def orpheus_sources() -> list[tuple[str, str]]:
    """``(posix path relative to the root, source)`` of every module of ``orpheus/`` outside ``derivations``."""
    out = []
    for path in sorted((ROOT / "orpheus").rglob("*.py")):
        rel = path.relative_to(ROOT).as_posix()
        if "/derivations/" in rel:
            continue
        out.append((rel, path.read_text()))
    return out


def tree_sources(subdir: str) -> list[tuple[str, str]]:
    return [(p.relative_to(ROOT).as_posix(), p.read_text()) for p in sorted((ROOT / subdir).rglob("*.py"))]


def class_defs(source: str) -> list[str]:
    return [n.name for n in ast.walk(ast.parse(source)) if isinstance(n, ast.ClassDef)]


def called_names(source: str) -> list[str]:
    """The name of every called function or method (``f(...)`` and ``x.f(...)``) of a source."""
    out = []
    for n in ast.walk(ast.parse(source)):
        if isinstance(n, ast.Call):
            out.append(n.func.id if isinstance(n.func, ast.Name) else n.func.attr if isinstance(n.func, ast.Attribute) else "")
    return out


def attribute_reads(source: str, names: frozenset[str]) -> list[int]:
    """The line of every attribute node whose attribute is one of ``names``."""
    return [n.lineno for n in ast.walk(ast.parse(source)) if isinstance(n, ast.Attribute) and n.attr in names]


def accumulator_loops(source: str, names: frozenset[str]) -> list[int]:
    """The line of every ``NAME += <expr>`` (``NAME`` in ``names``) that sits inside a ``for``."""
    tree = ast.parse(source)
    hits: list[int] = []

    def visit(node: ast.AST, in_for: bool) -> None:
        for child in ast.iter_child_nodes(node):
            inside = in_for or isinstance(child, (ast.For, ast.AsyncFor))
            if in_for and isinstance(child, ast.AugAssign) and isinstance(child.op, ast.Add) and isinstance(child.target, ast.Name) and child.target.id in names:
                hits.append(child.lineno)
            visit(child, inside)

    visit(tree, False)
    return hits
