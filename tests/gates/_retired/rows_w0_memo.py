"""The rows of unit ``memo-digest`` (specification part B, 1.4 and 4.1 G2): the traced memo moves to ``orpheus/memo/`` and the hash
convention gets one definition, ``orpheus/values/digest.py``.

Sub-unit M (the move): 21 ``RetiredModule`` rows, one per old module path (decided B5, not 250 ``Moved`` rows): the three memo
modules and the eighteen old memo test modules (``git mv`` carries each test to ``tests.gates.memo``). Each is red while its old
file exists (the absence leg) and while any statement imports it.

Sub-unit D (the digest): 3 ``Retired`` rows. The boot's ``digest`` and the two width constants of ``content.py`` are gone once the
convention is the one definition. These three rows read green before the move (the boot file is absent at its new path) and red at
the end of M, where the boot still defines ``digest``: they belong to commit D.

The old paths are spelled from pieces here as everywhere outside this package (``tests/gates/test_memo_old_path.py`` reads the text of
every tracked file and excludes only this package).
"""

from __future__ import annotations

from tests.gates._retired import Kept, Retired, RetiredModule

_OLD = "orpheus." + "numerics."
_OLD_TESTS = "tests.gates." + "numerics."

MEMO_MODULES: tuple[str, ...] = ("traced_memo", "_traced_memo_boot", "_import_time")

MEMO_TEST_MODULES: tuple[str, ...] = (
    "test_traced_memo_apparatus", "test_traced_memo_compiled_bytes", "test_traced_memo_data", "test_traced_memo_findings",
    "test_traced_memo_listing_ignores_bytecode", "test_traced_memo_manifest", "test_traced_memo_nested_recording",
    "test_traced_memo_ownership", "test_traced_memo_process", "test_traced_memo_records_each_file", "test_traced_memo_source_index",
    "test_traced_memo_statement_pins", "test_traced_memo_template", "test_traced_memo_validation",
    "_traced_memo_api", "_traced_memo_synthetic", "_statement_pin_cases", "_statement_pin_qa_rows",
)

ROWS: tuple[RetiredModule, ...] = (
    *(RetiredModule(_OLD + name) for name in MEMO_MODULES),
    *(RetiredModule(_OLD_TESTS + name) for name in MEMO_TEST_MODULES),
)

KEPT: tuple[Kept, ...] = ()
