"""Gate G4 of unit ``memo-digest`` (specification part B 3.2, 4.2): the one hash convention, ``orpheus/values/digest.py``.

Kinds. The laws against ``hashlib`` spelled in this file are REFERENCE rows (an independent route to the same blake2b-256). The vector
table is a RECORD (the outputs of ``content_digest`` on the unmodified tree at ``e3025655``, frozen as literals that no refactor
can move); it is green at base by design and is the bit-identity witness for every stored key. The spelling census is a THEOREM
about the tree: ``hashlib`` is imported by ``orpheus/values/digest.py`` alone, and no digest width literal remains elsewhere.

``orpheus.values.digest`` is imported inside the test bodies, so its absence is one red per test and not a collection error. Every
census walk takes its file list (the planted legs), prints its input count, and refuses an empty list.

First reds at base (measured at the commit that adds this file): the laws, ``NAME_SIZE`` and the one-spelling census are red (the module
is absent; the census reads 0 files in ``orpheus/memo/``); the vector table is green. A scratch probe of the same predicates on the
base locations reads 6 ``hashlib.blake2b`` calls in 2 modules and 9 + 1 width literals (commit body).
"""
from __future__ import annotations

import ast
import hashlib
from pathlib import Path

import pytest

pytestmark = pytest.mark.foundation

REPO = Path(__file__).resolve().parents[3]

VECTORS = {
    "int0": "66647095f396ec024ede27f99211f8295bced12baa12de937a23ebeb994552ef",
    "int1": "8b63300be12e6a865d09db958abaae1f1d6a7c69806fd0163dad6361a2ece393",
    "float1": "8b63300be12e6a865d09db958abaae1f1d6a7c69806fd0163dad6361a2ece393",
    "negzero": "66647095f396ec024ede27f99211f8295bced12baa12de937a23ebeb994552ef",
    "true": "8b63300be12e6a865d09db958abaae1f1d6a7c69806fd0163dad6361a2ece393",
    "str": "a00d3b459b43674138a54caf9e7820263dea1960a220cb2085828737336aa0d1",
    "none": "9b5fb044c4b8983390c4e0ccf33e7e1cfd4cc200feafdb54aa9c040a87dba389",
    "tuple": "caba50112342afde2f0d5a8c5a004900b09911c8a2bd4f3e7635ccc43659bb9f",
    "arr": "90134377cc95fa4812ba116861ccda01d3fe386b33dd8ff79f20e8bc392be597",
    "arr2d": "eb9d4a47bc665abd6ff1cf377c7f6a776ef57c22f0da9342239407d7ec2b7eed",
    "frozenset": "fd063b0e7933cd430d199776ddd32af159b864b319e6ba182ee60c5ab01b1e66",
    "identity": "e65cd67fe1c3527b91d79629843f75c88a4fcd0bd4f49625f27cc548b296927b",
}


def _convention():
    """The module under test, imported at call time (an import STATEMENT: its absence is one red per gate, not a collection error)."""
    from orpheus.values import digest

    return digest


def _values() -> dict[str, object]:
    """Built INSIDE each test: a module-level ``FrozenMapping`` keeps its cached digest and would stay green under a mutation of the
    convention (``content._DIGESTS``)."""
    import numpy as np

    from orpheus.values.content import FrozenMapping

    return {"int0": 0, "int1": 1, "float1": 1.0, "negzero": -0.0, "true": True, "str": "a", "none": None,
            "tuple": (1, 2.5, "s"), "arr": np.arange(3.0), "arr2d": np.eye(2), "frozenset": frozenset({1, 2}),
            "identity": FrozenMapping({"k": 1.0, "j": (1, 2)})}


def _blobs() -> list[bytes]:
    return [b"", b"\x00", bytes(range(63)), bytes(range(64)), bytes(range(65)), bytes(1 << 20)]


# ---------------------------------------------------------------- the laws (REFERENCE: hashlib spelled here)
def test_the_convention_is_blake2b_256_as_bytes_and_as_hex() -> None:
    convention = _convention()
    assert convention.DIGEST_SIZE == 32 and convention.NAME_SIZE == 8
    for blob in _blobs():
        expected = hashlib.blake2b(blob, digest_size=32)
        assert len(convention.digest(blob)) == convention.DIGEST_SIZE
        assert convention.digest(blob) == expected.digest()
        assert convention.hexdigest(blob) == convention.digest(blob).hex() == expected.hexdigest()


def test_the_convention_is_total_on_buffers_and_refuses_text() -> None:
    convention = _convention()
    for buffer in (bytearray(b"abc"), memoryview(b"abc")):
        assert convention.digest(buffer) == hashlib.blake2b(b"abc", digest_size=32).digest()
    with pytest.raises(TypeError):
        convention.digest("abc")  # type: ignore[arg-type]


# ---------------------------------------------------------------- the vector table (RECORD) and the name prefix law
def test_the_content_digest_of_each_value_is_the_frozen_vector() -> None:
    """Red on any change of the encoder, the hash or the width: every stored key moves (8 526 manifests)."""
    from orpheus.values.content import content_digest, name_digest

    for name, value in _values().items():
        assert content_digest(value).hex() == VECTORS[name], name
        assert name_digest(value) == VECTORS[name][:16], name


def test_the_name_digest_is_a_prefix_of_the_content_digest_and_the_two_routes_agree() -> None:
    from orpheus.values.content import FrozenMapping, content_digest, name_digest

    prefix_width = _convention().NAME_SIZE
    for name, value in _values().items():
        assert name_digest(value) == content_digest(value)[:prefix_width].hex(), name
    spelled = {"k": 1.0, "j": (1, 2)}
    assert FrozenMapping(spelled).content_digest == content_digest(FrozenMapping(spelled))


# ---------------------------------------------------------------- the spelling census (THEOREM over the tree)
def hashlib_uses(source: str) -> list[int]:
    """Lines where ``source`` imports ``hashlib`` in any form or calls through one: ``import hashlib``, ``import hashlib as h``,
    ``from hashlib import blake2b as b``, ``importlib.import_module("hashlib")``, ``__import__("hashlib")``, ``h.sha256(...)``."""
    tree = ast.parse(source)
    lines: set[int] = set()
    aliases: set[str] = set()
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.split(".")[0] == "hashlib":
                    lines.add(node.lineno)
                    aliases.add(alias.asname or alias.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom) and node.module and node.module.split(".")[0] == "hashlib" and not node.level:
            lines.add(node.lineno)
            names |= {a.asname or a.name for a in node.names}
        elif isinstance(node, ast.Call) and node.args and isinstance(node.args[0], ast.Constant) and node.args[0].value == "hashlib":
            if ast.unparse(node.func).rsplit(".", 1)[-1] in ("import_module", "__import__"):
                lines.add(node.lineno)
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            func = node.func
            if (isinstance(func, ast.Attribute) and isinstance(func.value, ast.Name) and func.value.id in aliases) or (
                    isinstance(func, ast.Name) and func.id in names):
                lines.add(node.lineno)
    return sorted(lines)


def hashlib_calls(source: str) -> int:
    tree = ast.parse(source)
    aliases = {a.asname or a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names if a.name == "hashlib"}
    names = {a.asname or a.name for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.module == "hashlib" for a in n.names}
    return sum(1 for n in ast.walk(tree) if isinstance(n, ast.Call) and (
        (isinstance(n.func, ast.Attribute) and isinstance(n.func.value, ast.Name) and n.func.value.id in aliases)
        or (isinstance(n.func, ast.Name) and n.func.id in names)))


def width_literals(source: str) -> list[int]:
    """The integer literals 32 and 64 (a ``bool`` is not one), by line."""
    return sorted(n.lineno for n in ast.walk(ast.parse(source))
                  if isinstance(n, ast.Constant) and type(n.value) is int and n.value in (32, 64))


def _census_files(root: Path) -> list[Path]:
    """``orpheus/`` and ``tools/`` under ``root`` outside the closed ``orpheus/derivations/`` (stated exclusion)."""
    files = [p for tree in ("orpheus", "tools") for p in sorted((root / tree).rglob("*.py"))
             if "derivations" not in p.relative_to(root).parts[:2] and "__pycache__" not in p.parts]
    if not files:
        raise AssertionError(f"input count 0 under {root}")
    return files


def _spelling_findings(root: Path) -> dict[str, object]:
    importers = {p.relative_to(root).as_posix(): hashlib_uses(p.read_text(encoding="utf-8")) for p in _census_files(root)}
    importers = {k: v for k, v in importers.items() if v}
    calls = sum(hashlib_calls(p.read_text(encoding="utf-8")) for p in _census_files(root))
    guarded = [*sorted((root / "orpheus" / "memo").glob("*.py")), root / "orpheus" / "values" / "content.py",
               *sorted((root / "tools" / "cone_store").glob("*.py"))]
    literals = {p.relative_to(root).as_posix(): width_literals(p.read_text(encoding="utf-8")) for p in guarded if p.is_file()}
    return {"importers": sorted(importers), "calls": calls, "literals": {k: v for k, v in literals.items() if v},
            "memo_files": len(list((root / "orpheus" / "memo").glob("*.py")))}


def test_hashlib_is_imported_by_the_convention_alone_and_no_width_literal_remains() -> None:
    found = _spelling_findings(REPO)
    assert found["memo_files"] >= 4, f"input count {found['memo_files']} in orpheus/memo/ (the package and three modules)"
    assert found["importers"] == ["orpheus/values/digest.py"], found["importers"]
    assert found["calls"] == 1, f"{found['calls']} hashlib calls"
    assert found["literals"] == {}, found["literals"]
    source = (REPO / "orpheus" / "values" / "digest.py").read_text(encoding="utf-8")
    assert len(width_literals(source)) == 1
    assert sum(1 for n in ast.walk(ast.parse(source)) if isinstance(n, ast.Constant) and type(n.value) is int and n.value == 8) == 1


@pytest.mark.parametrize(("label", "source", "uses", "calls"), [
    ("import", "import hashlib\nhashlib.blake2b(b'')\n", [1, 2], 1),
    ("import-as", "import hashlib as h\nh.sha256(b'')\n", [1, 2], 1),
    ("from-import-as", "from hashlib import blake2b as b\nb(b'')\n", [1, 2], 1),
    ("import_module", "import importlib\nimportlib.import_module('hashlib')\n", [2], 0),
    ("dunder-import", "m = __import__('hashlib')\n", [1], 0),
    ("negative-name", "import hashlibx\nhashlibx.f()\nx = 'hashlib'\n", [], 0),
])
def test_the_hashlib_census_finds_each_of_five_shapes_and_reads_its_negative_as_zero(label: str, source: str, uses: list[int], calls: int) -> None:
    assert hashlib_uses(source) == uses, label
    assert hashlib_calls(source) == calls, label


def test_the_width_census_finds_a_slice_and_a_keyword_and_skips_bools_and_other_ints() -> None:
    assert width_literals("def f(d):\n    return d[:32], d[32:64], g(digest_size=32), 64\n") == [2, 2, 2, 2, 2]  # 32; 32, 64; 32; 64
    assert width_literals("x = (True, 8, 33, 2 * 3)\n") == []


def test_a_planted_spelling_in_a_copy_of_the_tree_is_read_by_the_census(tmp_path: Path) -> None:
    """The census walks a root: a copy with one planted ``hashlib`` alias and one planted slice is red (the two arms of 4.2 G4 (2))."""
    (tmp_path / "orpheus" / "memo").mkdir(parents=True)
    (tmp_path / "orpheus" / "values").mkdir()
    (tmp_path / "tools" / "cone_store").mkdir(parents=True)
    (tmp_path / "orpheus" / "values" / "digest.py").write_text("import hashlib\nhashlib.blake2b(b'')\nDIGEST_SIZE = 32\nNAME_SIZE = 8\n")
    (tmp_path / "orpheus" / "memo" / "a.py").write_text("x = 1\n")
    (tmp_path / "tools" / "cone_store" / "b.py").write_text("y = 2\n")
    clean = _spelling_findings(tmp_path)
    assert clean["importers"] == ["orpheus/values/digest.py"] and clean["literals"] == {} and clean["calls"] == 1
    (tmp_path / "orpheus" / "memo" / "a.py").write_text("import hashlib as h\nh.sha256(b'')\n")
    (tmp_path / "tools" / "cone_store" / "b.py").write_text("def f(data):\n    return data[:32]\n")
    planted = _spelling_findings(tmp_path)
    assert planted["importers"] == ["orpheus/memo/a.py", "orpheus/values/digest.py"] and planted["calls"] == 2
    assert planted["literals"] == {"tools/cone_store/b.py": [2]}
    with pytest.raises(AssertionError, match="input count 0"):
        _census_files(tmp_path / "empty")
