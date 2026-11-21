"""The rows of unit ``roles`` (W0, specification ``ontology_w0_spec.md`` section 4, R1): the block and system role tags retire.

25 rows (13 ``Retired`` by definition site, 12 ``RetiredName`` by identifier) and one ``Kept`` entry. Every row reads RED at the base
``954c4a4f`` (the definitions of ``orpheus/numerics/operator.py`` and the stamp are there) and green once the unit's change lands.

The 12 names are ``BlockRole SystemRole _BlockRoleMeta BulkOperator FullOperator BoundaryOperator block_role system_role
_join_block_roles _join_system_roles stamp_boundary_role _adopt_expansion_roles``. ``BoundaryOperator`` is a whole identifier, so the
live ``SNBoundaryOperator``, ``DiffusionBoundaryOperator`` and the like are other identifiers and are not read by its row. A
``Retired`` row binds a definition site; a ``RetiredName`` row reads the bare word anywhere in the five trees (identifier, ``__all__``
member, ``getattr``/``setattr`` string); each name below has 0 surviving homonym (specification 4.2 R1; the counts are the trailing comments of the rows, measured with the counter of the 4.5 list).

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
    # the 12 bare identifiers (counts: identifier tokens and string constants equal to the name over the five trees at 954c4a4f, `git ls-tree` + `tokenize` + `ast`; 0 homonym = the surviving count outside this module's own two strings per name)
    RetiredName("BlockRole"),  # 105 identifier occurrences at 954c4a4f, 0 homonym
    RetiredName("SystemRole"),  # 41 identifier occurrences at 954c4a4f, 0 homonym
    RetiredName("_BlockRoleMeta"),  # 4 identifier occurrences at 954c4a4f, 0 homonym
    RetiredName("BulkOperator"),  # 22 identifier occurrences at 954c4a4f, 0 homonym
    RetiredName("FullOperator"),  # 28 identifier occurrences at 954c4a4f, 0 homonym
    RetiredName("BoundaryOperator"),  # 30 identifier occurrences at 954c4a4f, 0 homonym
    RetiredName("block_role"),  # 91 identifier occurrences at 954c4a4f, 0 homonym
    RetiredName("system_role"),  # 35 identifier occurrences at 954c4a4f, 0 homonym
    RetiredName("_join_block_roles"),  # 3 identifier occurrences at 954c4a4f, 0 homonym
    RetiredName("_join_system_roles"),  # 14 identifier occurrences at 954c4a4f, 0 homonym
    RetiredName("stamp_boundary_role"),  # 14 identifier occurrences at 954c4a4f, 0 homonym
    RetiredName("_adopt_expansion_roles"),  # 4 identifier occurrences at 954c4a4f, 0 homonym
)
KEPT: tuple[Kept, ...] = (
    Kept(_OPERATOR, "D2", "row 63: the TypeVar of the quoted annotation of LinearOperator.__matmul__"),
)
