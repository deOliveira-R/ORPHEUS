r"""#602 round two, D: the distribution that owns a site-packages file is the one whose ``RECORD`` lists it.

Spec ``scratch/characteristic_architecture/p3/m602/round2_gate_spec.md`` §2, rows D1–D8.

The manifest pins a site-packages file by the version of the distribution that owns it. Ownership used to be
read off the file's top-level name (``importlib.metadata.packages_distributions``, first claimant
alphabetically). A namespace top such as ``sphinxcontrib`` is claimed by several distributions, so that map
pinned ``sphinxcontrib/tikz.py`` to ``sphinxcontrib-applehelp`` (``[M]`` 2026-10-15). An upgrade of
``sphinxcontrib-tikz`` would then have left valid every entry that ran its code. Ownership is now per FILE:
the distribution whose ``RECORD`` lists it, built once per checkout.

The fake site directory below is a ``tmp_path`` monkeypatched in as the only site-packages directory, holding
hand-written ``dist-info`` directories; no installed file is touched.
"""
from __future__ import annotations

import importlib
import importlib.metadata
import json
import os
import site
import sys
import uuid
from pathlib import Path

import pytest

from . import _traced_memo_api as api

pytestmark = pytest.mark.foundation

_M15 = 'tests/gates/numerics/test_traced_memo_manifest.py::test_m1_5_every_shape_of_traced_file_is_classified'
_D2 = 'tests/gates/numerics/test_traced_memo_ownership.py::test_d2_files_under_a_shared_namespace_top_belong_to_their_own_distributions'


class _FakeSite:
    """A site-packages directory of hand-written distributions, the only one ``origin`` sees."""

    def __init__(self, root: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        self.dir = Path(os.path.realpath(root)) / "site-packages"
        self.dir.mkdir(parents=True)
        monkeypatch.setattr(api.module(), "_site_dirs", lambda: (str(self.dir) + os.sep,))

    def install(self, name: str, version: str, files: dict[str, str], *, record_rows: list[str] | None = None,
                editable: bool = False) -> Path:
        info = self.dir / f"{name.replace('-', '_')}-{version}.dist-info"
        info.mkdir()
        (info / "METADATA").write_text(f"Metadata-Version: 2.1\nName: {name}\nVersion: {version}\n")
        for relative, text in files.items():
            path = self.dir / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)
        rows = record_rows if record_rows is not None else [f"{relative},sha256=x,1" for relative in files]
        (info / "RECORD").write_text("\n".join([*rows, f"{info.name}/RECORD,,"]) + "\n")
        if editable:
            (info / "direct_url.json").write_text(json.dumps({"url": "file:///x", "dir_info": {"editable": True}}))
        return info

    def origin(self, relative: str):
        return api.name("origin")(str(self.dir / relative), api.name("_InstalledFiles")())


@pytest.fixture
def fake_site(tmp_path, monkeypatch):
    return _FakeSite(tmp_path, monkeypatch)


def _kind_and_name(found) -> tuple[str, str]:
    return type(found).__name__, found[0]


def _site_files_of_a_characteristic_block() -> list[str]:
    import orpheus.derivations.continuous.characteristic  # noqa: F401  (what a block imports)

    sites = [os.path.realpath(d) + os.sep for d in site.getsitepackages()]
    files = {os.path.realpath(f) for m in list(sys.modules.values()) if isinstance(f := getattr(m, "__file__", None), str)}
    return sorted(f for f in files if f.startswith(tuple(sites)))


@pytest.mark.rests_on(_M15)
def test_d1_every_imported_site_file_is_owned_by_the_distribution_whose_record_lists_it():
    """D1, the real tree: each site-packages file a characteristic block imports belongs to a distribution whose
    file list, read by the standard library's own ``RECORD`` parser (``Distribution.files``, an independent
    route), contains it; the editable install's finder is pinned by content. The namespace witnesses
    (``sphinxcontrib/tikz.py``, ``jaraco/context/__init__.py``) belong to their own distributions, not to
    the first claimant of their top-level name. First red: ownership by the top-level name's first claimant
    (``sphinxcontrib-applehelp``)."""
    installed = api.name("_InstalledFiles")()  # one table for the whole row, as one manifest holds one

    def origin(path: str):
        return api.name("origin")(path, installed)

    site_dir = os.path.realpath(site.getsitepackages()[0]) + os.sep
    files = _site_files_of_a_characteristic_block()
    witnesses = [site_dir + w for w in ("sphinxcontrib/tikz.py", "jaraco/context/__init__.py") if os.path.isfile(site_dir + w)]
    assert len(files) > 100 and witnesses, (len(files), witnesses)
    listed: dict[str, set[str]] = {}
    wrong = []
    for path in files + witnesses:
        kind, name = _kind_and_name(origin(path))
        if kind == "Source":
            assert "__editable__" in os.path.basename(path), f"{path}: a site file pinned by content"
            continue
        assert kind == "Distribution", (path, kind)
        if name not in listed:
            listed[name] = {str(p) for p in (importlib.metadata.distribution(name).files or [])}
        if os.path.relpath(path, site_dir) not in listed[name]:
            wrong.append((os.path.relpath(path, site_dir), name))
    assert not wrong, f"{len(wrong)} of {len(files) + len(witnesses)} files pinned to a distribution that does not list them: {wrong[:5]}"
    names = {os.path.relpath(w, site_dir): origin(w)[0] for w in witnesses}
    assert names.get("sphinxcontrib/tikz.py", "sphinxcontrib-tikz") == "sphinxcontrib-tikz", names
    assert names.get("jaraco/context/__init__.py", "jaraco.context") == "jaraco.context", names


@pytest.mark.rests_on(_M15)
def test_d2_files_under_a_shared_namespace_top_belong_to_their_own_distributions(fake_site):
    """D2: two distributions share the namespace top ``ns``; each file belongs to the distribution that lists it.
    First red: ownership by the top-level name's first claimant, which gives ``ns/b.py`` to ``ns-a``."""
    fake_site.install("ns-a", "1.0", {"ns/a.py": "A = 1\n"})
    fake_site.install("ns-b", "2.0", {"ns/b.py": "B = 2\n"})
    assert _kind_and_name(fake_site.origin("ns/a.py")) == ("Distribution", "ns-a")
    assert _kind_and_name(fake_site.origin("ns/b.py")) == ("Distribution", "ns-b")


@pytest.mark.rests_on(_D2)
def test_d3_a_file_of_an_editable_install_is_pinned_by_its_content(fake_site):
    """D3: a file listed by an editable install's ``RECORD`` (``direct_url.json`` ``dir_info.editable``) is
    source, pinned by its bytes, since an editable version never moves; the same file of a regular install is
    the distribution. First red: the editable flag dropped."""
    fake_site.install("ed", "0.1", {"__editable___ed_finder.py": "X = 1\n"}, editable=True)
    fake_site.install("reg", "0.1", {"reg_finder.py": "X = 1\n"})
    assert _kind_and_name(fake_site.origin("__editable___ed_finder.py"))[0] == "Source"
    assert _kind_and_name(fake_site.origin("reg_finder.py")) == ("Distribution", "reg")


@pytest.mark.rests_on(_D2)
def test_d4_a_record_path_holding_a_comma_is_read_as_csv(fake_site):
    """D4: ``RECORD`` is CSV; a path that holds a comma is quoted, and the file it names is owned. First red: a
    ``RECORD`` row split on its first comma, which reads the path as ``"ns/a``."""
    fake_site.install("comma", "1.0", {"ns/a,b.py": "X = 1\n"}, record_rows=['"ns/a,b.py",sha256=x,1'])
    assert _kind_and_name(fake_site.origin("ns/a,b.py")) == ("Distribution", "comma")


@pytest.mark.rests_on(_D2)
def test_d5_a_site_file_listed_by_no_record_is_refused(fake_site):
    """D5 (ruled 2026-10-15): a site-packages file that no ``RECORD`` lists has no version that determines it,
    even under a top another distribution claims, so it is refused. First red: a fallback to a claimant's
    version."""
    fake_site.install("claimant", "1.0", {"ns/listed.py": "X = 1\n"})
    (fake_site.dir / "ns" / "stray.py").write_text("X = 2\n")
    with pytest.raises(api.name("Unpinnable"), match="listed by no installed distribution's RECORD"):
        fake_site.origin("ns/stray.py")
    assert _kind_and_name(fake_site.origin("ns/listed.py")) == ("Distribution", "claimant")


@pytest.mark.rests_on(_D2)
def test_d6_an_upgrade_between_two_manifests_in_one_process_is_seen(fake_site, monkeypatch):
    """D6: a distribution upgraded between two traced runs in one process (its ``dist-info`` renamed, its
    version bumped, a new file listed) is read afresh by the second manifest: ownership lives one checkout.
    First red: an ownership table cached for the process, which still names the removed ``dist-info``."""
    module = f"own6_{uuid.uuid4().hex[:8]}"
    info = fake_site.install("own6", "1.0", {f"{module}.py": "def f(x):\n    return x + 1.0\n"})
    monkeypatch.syspath_prepend(str(fake_site.dir))
    importlib.invalidate_caches()
    f = importlib.import_module(module).f
    _, first = api.trace_call(f, 1.0)
    assert ("own6", "1.0") in set(first.distributions), first.distributions
    upgraded = info.with_name(f"own6-2.0.dist-info")
    info.rename(upgraded)
    (upgraded / "METADATA").write_text("Metadata-Version: 2.1\nName: own6\nVersion: 2.0\n")
    importlib.invalidate_caches()
    _, second = api.trace_call(f, 1.0)
    assert ("own6", "2.0") in set(second.distributions), second.distributions
    sys.modules.pop(module, None)


def test_d7_the_top_level_map_is_retired(monkeypatch):
    """D7, retirement witnesses: the two routes to a file's distribution that preceded the ownership table are
    gone from the module, and a traced run's manifest never builds the top-level map. The spy is checked to
    count (one direct call is counted) before its zero is read. First red: the module before D."""
    memo = api.module()
    assert not hasattr(memo, "_distributions_by_top_level")
    assert not hasattr(memo, "_recording_distribution")
    calls: list[int] = []
    real = importlib.metadata.packages_distributions
    monkeypatch.setattr(importlib.metadata, "packages_distributions", lambda: (calls.append(1), real())[1])
    importlib.metadata.packages_distributions()
    assert calls == [1], "the spy does not count"
    calls.clear()
    import numpy as np

    _, manifest = api.trace_call(lambda: float(np.linalg.norm(np.ones(3))))
    assert any(name == "numpy" for name, _ in manifest.distributions), manifest.distributions
    assert calls == [], f"the manifest built the top-level map {len(calls)} time(s)"


def test_d8_a_distribution_pin_carries_the_installed_name_and_version():
    """D8: the pins of a traced run that executes numpy and scipy name each distribution as its ``METADATA``
    does and hold its installed version, so a validation looks up the same name. First red: a name spelled
    from the ``dist-info`` directory (``sphinxcontrib_tikz`` against ``sphinxcontrib-tikz``; D1 reddens on it
    too, these rows carry no such name)."""
    import scipy.special

    _, manifest = api.trace_call(lambda: float(scipy.special.roots_legendre(4)[1].sum()))  # Python-level code of both
    pins = set(manifest.distributions)
    assert {"numpy", "scipy"} <= {name for name, _ in pins}, pins
    for name, version in pins:
        assert importlib.metadata.metadata(name)["Name"] == name, name
        assert importlib.metadata.version(name) == version, (name, version)
    assert api.validate(manifest) == ()
