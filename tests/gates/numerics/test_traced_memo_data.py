r"""Step 3 of #405 P3, the rider on data files: what a generation READ is in the manifest, not only what ran.

Spec ``.claude/plans/reference_p3_spec.md`` §1.3, gate M3.17 (an open question to the user, §4 Q4, with this
gate as its evidence). The trace records code; a generator that reads a data file beside its module (a table, an
``.npz``, nuclear data) answers from bytes no AST holds, so an edit to the file would be served stale. The
recommended remedy is an audit hook on ``open`` in the generating process that pins every file it opened for
reading by its bytes. Its blind spot, named: a C library that opens a file itself (HDF5 through ``h5py``) raises
no Python audit event.
"""
from __future__ import annotations

import importlib
import json
import sys
import textwrap
import uuid

import pytest

from . import _traced_memo_api as api

pytestmark = pytest.mark.foundation

READER = '''
import json
from pathlib import Path

from {memo_module} import traced_memo

TABLE = Path(__file__).with_name("table.json")


@traced_memo
def lookup(key: str) -> float:
    return float(json.loads(TABLE.read_text())[key])
'''


@pytest.fixture
def reader(tmp_path):
    name = f"memo_data_{uuid.uuid4().hex[:12]}"
    package = tmp_path / "src" / name
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("")
    (package / "reader.py").write_text(textwrap.dedent(READER).format(memo_module=api.MODULE))
    (package / "table.json").write_text(json.dumps({"a": 1.5}))
    sys.path.insert(0, str(tmp_path / "src"))
    importlib.invalidate_caches()
    with api.cache_root(tmp_path / "cache"):
        yield importlib.import_module(f"{name}.reader"), package / "table.json"
    sys.path.remove(str(tmp_path / "src"))
    for key in [k for k in sys.modules if k.startswith(name)]:
        del sys.modules[key]


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_1_a_miss_generates_once_and_a_hit_never')
def test_m3_17_an_edit_to_a_data_file_the_generation_read_is_a_miss(reader):
    """M3.17: ``lookup("a")`` reads ``table.json``; after the table's value changes, the lookup is ``Stale``
    and the call returns the NEW value. Without a record of files read, the entry stays a ``Hit`` and serves
    1.5: the stale hit the trace cannot see (``[M]`` 2026-10-04 on the prototype without the hook: the red)."""
    module, table = reader
    assert module.lookup("a") == 1.5
    table.write_text(json.dumps({"a": 2.5}))
    assert api.verdict_kind(module.lookup.lookup("a")) == "Stale"
    assert module.lookup("a") == 2.5
