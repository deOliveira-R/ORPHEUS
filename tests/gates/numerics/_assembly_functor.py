"""The assembly functor's shared fixtures (numerics_untangling S3): integer and fused-atom leaves, the registry of
numerics operator kinds with their expected flag, and the G10 law (the flag is the functor's domain), stated once.

Read by ``tests/gates/numerics/test_assembly_functor.py`` (the numerics kinds) and by the census
``tests/gates/test_assembly_functor_census.py`` (the consumers' kinds and the union; layer DAG Phase 2, ruling 28).
Not collected (the leading underscore). The asserts are in collected test functions' frames only through this
module's callers, so they are ``-O``-stripped here: the law raises ``AssertionError`` explicitly instead.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy import sparse

from orpheus.numerics.assembled_operator import MissingAssembly, SparseAssembledOperator, assemble
from orpheus.numerics.elementary_operators import DiagonalOperator
from orpheus.numerics.flat_operator import FlattenedOperator
from orpheus.numerics.operator import FusedAtom, IdentityOperator, LinearOperator


def leaf(dense_matrix: np.ndarray) -> SparseAssembledOperator:
    return SparseAssembledOperator(sparse.csr_array(dense_matrix))


def integer(seed: int) -> np.ndarray:
    return np.random.default_rng(seed).integers(-3, 4, (4, 4)).astype(float)


class Atom(FusedAtom):
    """A fused atom built as production atoms are: an expansion, an action
    of its own, the expansion's roles adopted at construction."""

    def __init__(self, expansion: LinearOperator) -> None:
        self._expansion = expansion
        self._adopt_expansion_roles()

    def expansion(self) -> LinearOperator:
        return self._expansion

    def apply(self, x, /):
        return self._expansion.apply(x)


def registry() -> dict[str, "Any"]:
    """Every NUMERICS operator kind the S3 fixtures build, keyed by name and built lazily (a rebinding battery arm
    must reach the rows: lessons L100). The consumers' kinds (diffusion, homogeneous) are the root census's,
    ``tests/gates/test_assembly_functor_census.py`` (layer DAG Phase 2 refresh, ruling 28)."""

    def coupled(**kw):
        from tests.gates.numerics.test_coupled_operator import _grid

        return _grid("full", **kw)

    class _Duck:
        def to_flat(self):
            return np.zeros(4)

    a, b = leaf(integer(20)), leaf(integer(21))
    ident = IdentityOperator()
    diag = DiagonalOperator(np.arange(1.0, 5.0), shape=(4,), spans=(0,))
    builders = {
        "carrier": lambda: a,
        "sum": lambda: a - 0.3 * b,
        "word": lambda: a @ b,
        "nested": lambda: a @ (b + 2.0 * a) - b,
        "atom": lambda: Atom(a - 2.0 * b),
        "sum_with_identity": lambda: a + ident,
        "word_with_diagonal": lambda: a @ diag,
        "atom_over_identity": lambda: Atom(a + ident),
        "identity": lambda: ident,
        "diagonal": lambda: diag,
        "dual_of_carrier": lambda: a.dual(),
        "inverse_of_diagonal": lambda: diag.inverse(),
        "flattened": lambda: FlattenedOperator(a + b, _Duck()),
        "flattened_identity": lambda: FlattenedOperator(ident, _Duck()),
        "coupled": lambda: coupled(),
        "coupled_silent_block": lambda: coupled(assemblable=False),
    }
    return builders


EXPECTED_ASSEMBLABLE = {"carrier", "sum", "word", "nested", "atom", "flattened", "coupled"}


def assert_is_assemblable_iff_assemble_succeeds(name: str, op: "Any", expected: set[str]) -> None:
    """G10, the law ONCE: the flag ``is_assemblable`` is the functor's domain. The flag says True exactly when
    ``assemble`` returns and False exactly when it refuses with :class:`MissingAssembly`, and it equals the expected
    kind (``name in expected``). The composites' flags live in ``operator.py`` beside the other capability axes; the
    functor lives downstream (user ruling on elegance F2)."""
    flag = op.is_assemblable
    if flag != (name in expected):
        raise AssertionError(f"{name}: flag {flag}")
    try:
        assemble(op)
    except MissingAssembly:
        succeeded = False
    else:
        succeeded = True
    if flag != succeeded:
        raise AssertionError(f"{name}: is_assemblable {flag}, assemble succeeded {succeeded}")
