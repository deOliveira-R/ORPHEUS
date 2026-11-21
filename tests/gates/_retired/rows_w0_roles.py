"""The rows of unit ``roles`` (W0, specification ``ontology_w0_spec.md`` section 4, R1): the block and system role tags retire.

25 rows (13 ``Retired`` by definition site, 12 ``RetiredName`` by identifier) and one ``Kept`` entry. Every row reads RED at the base
``954c4a4f`` (the definitions of ``orpheus/numerics/operator.py`` and the stamp are there) and green once the unit's change lands.

The 12 names are ``BlockRole SystemRole _BlockRoleMeta BulkOperator FullOperator BoundaryOperator block_role system_role
_join_block_roles _join_system_roles stamp_boundary_role _adopt_expansion_roles``. ``BoundaryOperator`` is a whole identifier, so the
live ``SNBoundaryOperator``, ``DiffusionBoundaryOperator`` and the like are other identifiers and are not read by its row. A
``Retired`` row binds a definition site; a ``RetiredName`` row reads the bare word anywhere in the five trees (identifier, ``__all__``
member, ``getattr``/``setattr`` string); each name below has 0 surviving homonym (specification 4.2 R1, ``roles_count.py``).

``D2`` is kept: the ``TypeVar`` of the quoted annotation of ``LinearOperator.__matmul__`` (row 63), read only inside a string, so a
builder pruning the imports of ``operator.py`` could delete it silently; the engine asserts it is still defined.

Red input of each row: the definition or the identifier named, present in the tree today. Mutation for the reviewers (built at run
time from the row's own fields; no planted source is stored): copy ``operator.py`` alone to a temporary root at its relative path,
require the unmutated walk green, append ``class LinearOperator:`` with the body ``block_role = None``, and require the
``LinearOperator.block_role`` ``Retired`` row and the ``block_role`` ``RetiredName`` row red.
"""

from __future__ import annotations

from tests.gates._retired import Kept, Retired, RetiredName, Row

_OPERATOR = "orpheus/numerics/operator.py"

ROWS: tuple[Row, ...] = (
    # the 13 definition sites
    Retired(_OPERATOR, "BlockRole"),
    Retired(_OPERATOR, "SystemRole"),
    Retired(_OPERATOR, "_BlockRoleMeta"),
    Retired(_OPERATOR, "BulkOperator"),
    Retired(_OPERATOR, "FullOperator"),
    Retired(_OPERATOR, "BoundaryOperator"),
    Retired(_OPERATOR, "_join_block_roles"),
    Retired(_OPERATOR, "_join_system_roles"),
    Retired(_OPERATOR, "FusedAtom._adopt_expansion_roles"),
    Retired(_OPERATOR, "LinearOperator.block_role"),
    Retired(_OPERATOR, "LinearOperator.system_role"),
    Retired("orpheus/transport/boundary_conditions/realization.py", "stamp_boundary_role"),
    Retired("orpheus/sn/boundary/realized_operators.py", "RealizedBoundary.block_role"),
    # the 12 bare identifiers (each with 0 homonym)
    RetiredName("BlockRole"),
    RetiredName("SystemRole"),
    RetiredName("_BlockRoleMeta"),
    RetiredName("BulkOperator"),
    RetiredName("FullOperator"),
    RetiredName("BoundaryOperator"),
    RetiredName("block_role"),
    RetiredName("system_role"),
    RetiredName("_join_block_roles"),
    RetiredName("_join_system_roles"),
    RetiredName("stamp_boundary_role"),
    RetiredName("_adopt_expansion_roles"),
)
KEPT: tuple[Kept, ...] = (
    Kept(_OPERATOR, "D2", "row 63: the TypeVar of the quoted annotation of LinearOperator.__matmul__"),
)
