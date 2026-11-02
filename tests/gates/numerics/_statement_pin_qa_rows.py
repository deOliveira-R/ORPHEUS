"""qa's 42 attack rows on the B2 design (``scratch/verification_tree/cone_store/b2/qa_b2_design.md``), transcribed
by ``scratch/verification_tree/cone_store/b2/ta/gen_qa_rows.py`` from ``attacks/attacks.py`` (do not edit by hand;
re-run the generator). Each row: the files of a package ``pkg``, ONE edited file with its whole new text, and a run
snippet that sets ``RESULT``. ``MEMO_ROWS`` (38) run in a traced call; ``PYTEST_ROWS`` (4) run pytest on a test module
and reach only the cone store (a test module's constants)."""

MEMO_ROWS = [{'id': 'PC1',
  'title': 'selected constant edited -> must be STALE',
  'files': {'pkg/__init__.py': '', 'pkg/c.py': 'TOL = 1e-6\nOTHER = 2\n'},
  'file': 'pkg/c.py',
  'new': 'TOL = 1e-3\nOTHER = 2\n',
  'run': 'from pkg import c\nRESULT = c.TOL\n'},
 {'id': 'PC2',
  'title': 'unselected constant edited -> valid AND unchanged (the B2 benefit)',
  'files': {'pkg/__init__.py': '', 'pkg/c.py': 'TOL = 1e-6\nOTHER = 2\n'},
  'file': 'pkg/c.py',
  'new': 'TOL = 1e-6\nOTHER = 3\n',
  'run': 'from pkg import c\nRESULT = c.TOL\n'},
 {'id': 'S2c',
  'title': "registry function's default, reached by a registry lookup through a non-first-party decorator root",
  'files': {'pkg/__init__.py': '',
            'pkg/r.py': 'import functools\n'
                        'SCALE = 2\n'
                        '@functools.singledispatch\n'
                        'def f(x):\n'
                        '    return 0\n'
                        '@f.register\n'
                        'def _(x: int, s=SCALE):\n'
                        '    return x * s\n'},
  'file': 'pkg/r.py',
  'new': 'import functools\n'
         'SCALE = 3\n'
         '@functools.singledispatch\n'
         'def f(x):\n'
         '    return 0\n'
         '@f.register\n'
         'def _(x: int, s=SCALE):\n'
         '    return x * s\n',
  'run': 'from pkg import r\nRESULT = r.f(5)\n'},
 {'id': 'S2d',
  'title': 'side-effect-only import removed (the plugin module registers into a registry)',
  'files': {'pkg/__init__.py': 'from . import plugins\n',
            'pkg/reg.py': 'REGISTRY = {}\ndef register(k, v):\n    REGISTRY[k] = v\n',
            'pkg/plugins.py': "from .reg import register\nregister('a', 1)\n"},
  'file': 'pkg/__init__.py',
  'new': '',
  'run': 'import pkg\nfrom pkg.reg import REGISTRY\nRESULT = sorted(REGISTRY)\n'},
 {'id': 'S2e',
  'title': '__init_subclass__ registry: a new subclass (no Call node at class creation)',
  'files': {'pkg/__init__.py': '',
            'pkg/b.py': 'class Base:\n'
                        '    registry = []\n'
                        '    def __init_subclass__(cls):\n'
                        '        Base.registry.append(cls.__name__)\n',
            'pkg/m.py': 'from .b import Base\nclass A(Base):\n    pass\n'},
  'file': 'pkg/m.py',
  'new': 'from .b import Base\nclass A(Base):\n    pass\nclass B(Base):\n    pass\n',
  'run': 'from pkg import m\nfrom pkg.b import Base\nRESULT = len(Base.registry)\n'},
 {'id': 'S2f',
  'title': 'np.seterr bound to an unread name (third-party effect)',
  'files': {'pkg/__init__.py': '',
            'pkg/e.py': 'import numpy as np\n'
                        "_OLD = np.seterr(divide='ignore')\n"
                        'def ratio(a, b):\n'
                        '    return float(np.float64(a) / np.float64(b))\n'},
  'file': 'pkg/e.py',
  'new': 'import numpy as np\n'
         "_OLD = np.seterr(divide='raise')\n"
         'def ratio(a, b):\n'
         '    return float(np.float64(a) / np.float64(b))\n',
  'run': 'from pkg import e\nRESULT = e.ratio(1, 0)\n'},
 {'id': 'S2g',
  'title': 'warnings.filterwarnings bound (a with-free alias), DriftWarning-style filter',
  'files': {'pkg/__init__.py': '',
            'pkg/w.py': 'import warnings\n'
                        "_F = warnings.filterwarnings('ignore', category=UserWarning)\n"
                        'def f():\n'
                        "    warnings.warn('drift', UserWarning)\n"
                        '    return 1\n'},
  'file': 'pkg/w.py',
  'new': 'import warnings\n'
         "_F = warnings.filterwarnings('error', category=UserWarning)\n"
         'def f():\n'
         "    warnings.warn('drift', UserWarning)\n"
         '    return 1\n',
  'run': 'from pkg import w\nRESULT = w.f()\n'},
 {'id': 'S2h',
  'title': 'mutation through a builtin method in an assignment (dict.setdefault on a module object)',
  'files': {'pkg/__init__.py': '',
            'pkg/d.py': 'DEFAULTS = {}\n'
                        "_TOL = DEFAULTS.setdefault('tol', 1e-6)\n"
                        'def tol():\n'
                        "    return DEFAULTS['tol']\n"},
  'file': 'pkg/d.py',
  'new': "DEFAULTS = {}\n_TOL = DEFAULTS.setdefault('tol', 1e-2)\ndef tol():\n    return DEFAULTS['tol']\n",
  'run': 'from pkg import d\nRESULT = d.tol()\n'},
 {'id': 'S2i',
  'title': 'first-party call hidden behind a module-level alias root (fp classifier misses it)',
  'files': {'pkg/__init__.py': '',
            'pkg/r.py': 'REG = {}\n'
                        'class Registry:\n'
                        '    def add(self, k, v):\n'
                        '        REG[k] = v\n'
                        '        return v\n'
                        'registry = Registry()\n'
                        "_x = registry.add('k', 1)\n"},
  'file': 'pkg/r.py',
  'new': 'REG = {}\n'
         'class Registry:\n'
         '    def add(self, k, v):\n'
         '        REG[k] = v\n'
         '        return v\n'
         'registry = Registry()\n'
         "_x = registry.add('k', 2)\n",
  'run': "from pkg import r\nRESULT = r.REG['k']\n"},
 {'id': 'S2j',
  'title': 'for-loop registering at module level (subscript store nested in a compound statement)',
  'files': {'pkg/__init__.py': '', 'pkg/f.py': "REG = {}\nfor _k in ('a', 'b'):\n    REG[_k] = len(_k)\n"},
  'file': 'pkg/f.py',
  'new': "REG = {}\nfor _k in ('a', 'b', 'cc'):\n    REG[_k] = len(_k)\n",
  'run': 'from pkg import f\nRESULT = sorted(f.REG.items())\n'},
 {'id': 'S2k',
  'title': "import-time raise in an unselected statement (a 'dead' table referencing a removed def)",
  'files': {'pkg/__init__.py': '', 'pkg/t.py': "def old():\n    return 1\nTOL = 1e-6\n_LEGACY = {'old': old}\n"},
  'file': 'pkg/t.py',
  'new': "TOL = 1e-6\n_LEGACY = {'old': old}\n",
  'run': 'from pkg import t\nRESULT = t.TOL\n'},
 {'id': 'S2l',
  'title': 'import-time warning under an error filter (unselected statement now warns)',
  'files': {'pkg/__init__.py': "import warnings\nwarnings.simplefilter('error')\n",
            'pkg/t.py': 'import warnings\nTOL = 1e-6\n_UNUSED = 1\n'},
  'file': 'pkg/t.py',
  'new': "import warnings\nTOL = 1e-6\n_UNUSED = warnings.warn('drift') or 1\n",
  'run': 'from pkg import t\nRESULT = t.TOL\n'},
 {'id': 'S3a',
  'title': 'subclass in another module gains an override of an inherited class attribute',
  'files': {'pkg/__init__.py': '',
            'pkg/b.py': 'class Base:\n    tol = 1e-6\n    def get(self):\n        return self.tol\n',
            'pkg/c.py': 'from .b import Base\nclass Child(Base):\n    pass\n'},
  'file': 'pkg/c.py',
  'new': 'from .b import Base\nclass Child(Base):\n    tol = 1e-3\n',
  'run': 'from pkg.c import Child\nRESULT = Child().get()\n'},
 {'id': 'S3b',
  'title': 'dataclass gains a defaulted field (read by asdict, never by name)',
  'files': {'pkg/__init__.py': '',
            'pkg/p.py': 'import dataclasses\n@dataclasses.dataclass\nclass P:\n    a: int = 1\n    b: int = 2\n'},
  'file': 'pkg/p.py',
  'new': 'import dataclasses\n@dataclasses.dataclass\nclass P:\n    a: int = 1\n    b: int = 2\n    c: int = 3\n',
  'run': 'import dataclasses\nfrom pkg.p import P\nRESULT = sum(dataclasses.astuple(P()))\n'},
 {'id': 'S3c',
  'title': 'enum gains a member (read by iteration)',
  'files': {'pkg/__init__.py': '', 'pkg/k.py': 'import enum\nclass Kind(enum.Enum):\n    A = 1\n    B = 2\n'},
  'file': 'pkg/k.py',
  'new': 'import enum\nclass Kind(enum.Enum):\n    A = 1\n    B = 2\n    C = 3\n',
  'run': 'from pkg.k import Kind\nRESULT = [k.value for k in Kind]\n'},
 {'id': "S3c'",
  'title': 'enum member value edited (read by iteration, members never named)',
  'files': {'pkg/__init__.py': '', 'pkg/k.py': 'import enum\nclass Kind(enum.Enum):\n    A = 1\n    B = 2\n'},
  'file': 'pkg/k.py',
  'new': 'import enum\nclass Kind(enum.Enum):\n    A = 1\n    B = 5\n',
  'run': 'from pkg.k import Kind\nRESULT = [k.value for k in Kind]\n'},
 {'id': 'S3d',
  'title': 'PEP 562 module __getattr__ shadowed by a new binding',
  'files': {'pkg/__init__.py': 'def __getattr__(name):\n    return len(name)\n'},
  'file': 'pkg/__init__.py',
  'new': 'FOO = 99\ndef __getattr__(name):\n    return len(name)\n',
  'run': 'import pkg\nRESULT = pkg.FOO\n'},
 {'id': 'S3e',
  'title': 'submodule shadowed by a new binding in the package __init__',
  'files': {'pkg/__init__.py': '', 'pkg/thing.py': 'VALUE = 1\n'},
  'file': 'pkg/__init__.py',
  'new': 'class thing:\n    VALUE = 2\n',
  'run': 'from pkg import thing\nRESULT = thing.VALUE\n'},
 {'id': 'S3f',
  'title': 'unguarded read that raised at capture, succeeds after the edit (AttributeError caught by the caller)',
  'files': {'pkg/__init__.py': '', 'pkg/c.py': 'X = 1\n'},
  'file': 'pkg/c.py',
  'new': 'X = 1\nLATE = 7\n',
  'run': "from pkg import c\ntry:\n    RESULT = c.LATE\nexcept AttributeError:\n    RESULT = 'absent'\n"},
 {'id': 'S3g',
  'title': 'dataclass field reordered (positional construction, both names read)',
  'files': {'pkg/__init__.py': '',
            'pkg/p.py': 'import dataclasses\n@dataclasses.dataclass\nclass P:\n    a: int = 0\n    b: int = 0\n'},
  'file': 'pkg/p.py',
  'new': 'import dataclasses\n@dataclasses.dataclass\nclass P:\n    b: int = 0\n    a: int = 0\n',
  'run': 'from pkg.p import P\np = P(1, 2)\nRESULT = (p.a, p.b)\n'},
 {'id': 'N1',
  'title': 'operator.attrgetter with a dotted path',
  'files': {'pkg/__init__.py': '', 'pkg/o.py': 'class Inner:\n    val = 1\nclass Outer:\n    inner = Inner()\n'},
  'file': 'pkg/o.py',
  'new': 'class Inner:\n    val = 2\nclass Outer:\n    inner = Inner()\n',
  'run': "import operator\nfrom pkg.o import Outer\nRESULT = operator.attrgetter('inner.val')(Outer())\n"},
 {'id': 'N2',
  'title': 'str.format reading an attribute from C',
  'files': {'pkg/__init__.py': '', 'pkg/o.py': "class C:\n    label = 'x'\n"},
  'file': 'pkg/o.py',
  'new': "class C:\n    label = 'y'\n",
  'run': "from pkg import o\nRESULT = '{0.label}'.format(o.C())\n"},
 {'id': 'N3',
  'title': 'frozenset constant membership (co_consts recursion through tuples only)',
  'files': {'pkg/__init__.py': '', 'pkg/o.py': 'class C:\n    alpha = 1\n    beta = 2\n'},
  'file': 'pkg/o.py',
  'new': 'class C:\n    alpha = 1\n    beta = 5\n',
  'run': 'from pkg import o\n'
         "RESULT = [getattr(o.C, n) for n in sorted(n for n in vars(o.C) if n in {'alpha', 'beta', 'gamma'})]\n"},
 {'id': 'R1',
  'title': 'walrus at module level',
  'files': {'pkg/__init__.py': '', 'pkg/w.py': 'Y = (X := 3) + 1\n'},
  'file': 'pkg/w.py',
  'new': 'Y = (X := 4) + 1\n',
  'run': 'from pkg import w\nRESULT = w.X\n'},
 {'id': 'R2',
  'title': 'global rebind in a def',
  'files': {'pkg/__init__.py': '',
            'pkg/g.py': "MODE = 'a'\ndef set_mode():\n    global MODE\n    MODE = MODE + '!'\n"},
  'file': 'pkg/g.py',
  'new': "MODE = 'b'\ndef set_mode():\n    global MODE\n    MODE = MODE + '!'\n",
  'run': 'from pkg import g\ng.set_mode()\nRESULT = g.MODE\n'},
 {'id': 'R3',
  'title': 'del added',
  'files': {'pkg/__init__.py': '', 'pkg/d.py': 'X = 1\ndef f():\n    return X\n'},
  'file': 'pkg/d.py',
  'new': 'X = 1\ndel X\ndef f():\n    return X\n',
  'run': "from pkg import d\ntry:\n    RESULT = d.f()\nexcept NameError:\n    RESULT = 'gone'\n"},
 {'id': 'R4',
  'title': 'star-unpacking target',
  'files': {'pkg/__init__.py': '', 'pkg/s.py': 'A, *REST = (1, 2, 3)\n'},
  'file': 'pkg/s.py',
  'new': 'A, *REST = (1, 2, 4)\n',
  'run': 'from pkg import s\nRESULT = s.REST\n'},
 {'id': 'R5',
  'title': 'PEP 695 type alias',
  'files': {'pkg/__init__.py': '', 'pkg/t.py': 'type V = int\n'},
  'file': 'pkg/t.py',
  'new': 'type V = float\n',
  'run': 'from pkg import t\nRESULT = t.V.__value__\n'},
 {'id': 'R6',
  'title': '__slots__ change',
  'files': {'pkg/__init__.py': '',
            'pkg/s.py': "class C:\n    __slots__ = ('a',)\n    def __init__(self):\n        self.a = 1\n"},
  'file': 'pkg/s.py',
  'new': "class C:\n    __slots__ = ('a', '__dict__')\n    def __init__(self):\n        self.a = 1\n",
  'run': 'from pkg.s import C\n'
         'c = C()\n'
         'try:\n'
         "    c.z = 1; RESULT = 'ok'\n"
         'except AttributeError:\n'
         "    RESULT = 'slots'\n"},
 {'id': 'R7',
  'title': 'augmented assignment to a read name',
  'files': {'pkg/__init__.py': '', 'pkg/a.py': 'L = [1]\nL += [2]\n'},
  'file': 'pkg/a.py',
  'new': 'L = [1]\nL += [3]\n',
  'run': 'from pkg import a\nRESULT = a.L\n'},
 {'id': 'R8',
  'title': '__set_name__ descriptor (first-party call in the class body)',
  'files': {'pkg/__init__.py': '',
            'pkg/d.py': 'NAMES = []\n'
                        'class D:\n'
                        '    def __set_name__(self, owner, name):\n'
                        '        NAMES.append(name)\n'
                        'class C:\n'
                        '    x = D()\n'},
  'file': 'pkg/d.py',
  'new': 'NAMES = []\n'
         'class D:\n'
         '    def __set_name__(self, owner, name):\n'
         '        NAMES.append(name)\n'
         'class C:\n'
         '    y = D()\n',
  'run': 'from pkg import d\nRESULT = d.NAMES\n'},
 {'id': 'R9',
  'title': 'module-level lambda run after import',
  'files': {'pkg/__init__.py': '', 'pkg/l.py': 'K = 2\nSQ = lambda x: x * K\n'},
  'file': 'pkg/l.py',
  'new': 'K = 3\nSQ = lambda x: x * K\n',
  'run': 'from pkg import l\nRESULT = l.SQ(2)\n'},
 {'id': 'R10',
  'title': 're-export by assignment X = k.X across modules',
  'files': {'pkg/__init__.py': '', 'pkg/k.py': 'X = 1\n', 'pkg/r.py': 'from . import k\nX = k.X\n'},
  'file': 'pkg/k.py',
  'new': 'X = 2\n',
  'run': 'from pkg.r import X\nRESULT = X\n'},
 {'id': 'R11',
  'title': 'import a.b as c',
  'files': {'pkg/__init__.py': '',
            'pkg/sub/__init__.py': '',
            'pkg/sub/v.py': 'V = 1\n',
            'pkg/u.py': 'import pkg.sub.v as c\ndef get():\n    return c.V\n'},
  'file': 'pkg/sub/v.py',
  'new': 'V = 2\n',
  'run': 'from pkg import u\nRESULT = u.get()\n'},
 {'id': "S2c'",
  'title': 'singledispatch.register with a distinctive impl name, default edited',
  'files': {'pkg/__init__.py': '',
            'pkg/r.py': 'import functools\n'
                        'SCALE = 2\n'
                        '@functools.singledispatch\n'
                        'def f(x):\n'
                        '    return 0\n'
                        '@f.register\n'
                        'def qa_int_impl(x: int, s=SCALE):\n'
                        '    return x * s\n'},
  'file': 'pkg/r.py',
  'new': 'import functools\n'
         'SCALE = 3\n'
         '@functools.singledispatch\n'
         'def f(x):\n'
         '    return 0\n'
         '@f.register\n'
         'def qa_int_impl(x: int, s=SCALE):\n'
         '    return x * s\n',
  'run': 'from pkg import r\nRESULT = r.f(5)\n'},
 {'id': "N3'",
  'title': 'frozenset constant, no vars()',
  'files': {'pkg/__init__.py': '', 'pkg/o.py': 'class C:\n    qa_alpha = 1\n    qa_beta = 2\n'},
  'file': 'pkg/o.py',
  'new': 'class C:\n    qa_alpha = 1\n    qa_beta = 5\n',
  'run': "from pkg import o\nRESULT = sorted(getattr(o.C, n) for n in {'qa_alpha', 'qa_beta'})\n"},
 {'id': 'W1',
  'title': 'the wims.py shape: a group-boundary table read only by a module-level for loop',
  'files': {'pkg/__init__.py': '',
            'pkg/wims.py': '_LAST = [8, 11, 14]\n'
                           'CONDENSE = []\n'
                           '_start = 1\n'
                           'for _g, _last in enumerate(_LAST, start=1):\n'
                           '    CONDENSE.append((_g, _start, _last))\n'
                           '    _start = _last + 1\n'},
  'file': 'pkg/wims.py',
  'new': '_LAST = [8, 12, 14]\n'
         'CONDENSE = []\n'
         '_start = 1\n'
         'for _g, _last in enumerate(_LAST, start=1):\n'
         '    CONDENSE.append((_g, _start, _last))\n'
         '    _start = _last + 1\n',
  'run': 'from pkg.wims import CONDENSE\nRESULT = CONDENSE\n'},
 {'id': "S2k'",
  'title': 'a def deleted in module b, imported by name in module a whose binding the run never reads',
  'files': {'pkg/__init__.py': '',
            'pkg/b.py': 'def old():\n    return 1\n',
            'pkg/a.py': 'from .b import old\nTOL = 1\n'},
  'file': 'pkg/b.py',
  'new': 'def other():\n    return 1\n',
  'run': 'from pkg import a\nRESULT = a.TOL\n'}]

PYTEST_ROWS = [{'id': 'S2a',
  'title': 'pytest parametrize list: the test def is found by __dict__ iteration, never by name',
  'files': {'test_t.py': 'import pytest\n'
                         'CASES = [1, 2]\n'
                         "@pytest.mark.parametrize('x', CASES)\n"
                         'def test_x(x):\n'
                         '    assert x < 3\n'},
  'file': 'test_t.py',
  'new': "import pytest\nCASES = [1, 5]\n@pytest.mark.parametrize('x', CASES)\ndef test_x(x):\n    assert x < 3\n",
  'run': 'import pytest\n'
         "RESULT = int(pytest.main(['-q', '-p', 'no:cacheprovider', '--color=no', '-x', 'test_t.py']))\n"},
 {'id': 'S2b',
  'title': 'fixture params constant (conftest-style fixture found by argument name)',
  'files': {'test_t.py': 'import pytest\n'
                         'MESHES = [4, 8]\n'
                         '@pytest.fixture(params=MESHES)\n'
                         'def n(request):\n'
                         '    return request.param\n'
                         'def test_n(n):\n'
                         '    assert n % 4 == 0\n'},
  'file': 'test_t.py',
  'new': 'import pytest\n'
         'MESHES = [4, 9]\n'
         '@pytest.fixture(params=MESHES)\n'
         'def n(request):\n'
         '    return request.param\n'
         'def test_n(n):\n'
         '    assert n % 4 == 0\n',
  'run': 'import pytest\n'
         "RESULT = int(pytest.main(['-q', '-p', 'no:cacheprovider', '--color=no', '-x', 'test_t.py']))\n"},
 {'id': "S2b'",
  'title': 'fixture params constant, fixture with a distinctive name',
  'files': {'test_t.py': 'import pytest\n'
                         'MESHES = [4, 8]\n'
                         '@pytest.fixture(params=MESHES)\n'
                         'def qa_cells_zz(request):\n'
                         '    return request.param\n'
                         'def test_n(qa_cells_zz):\n'
                         '    assert qa_cells_zz % 4 == 0\n'},
  'file': 'test_t.py',
  'new': 'import pytest\n'
         'MESHES = [4, 9]\n'
         '@pytest.fixture(params=MESHES)\n'
         'def qa_cells_zz(request):\n'
         '    return request.param\n'
         'def test_n(qa_cells_zz):\n'
         '    assert qa_cells_zz % 4 == 0\n',
  'run': 'import pytest\n'
         "RESULT = int(pytest.main(['-q', '-p', 'no:cacheprovider', '--color=no', '-x', 'test_t.py']))\n"},
 {'id': "S2b''",
  'title': 'fixture params constant, fixture name never in an assert (no rewrite string)',
  'files': {'test_t.py': 'import pytest\n'
                         'MESHES = [4, 8]\n'
                         '@pytest.fixture(params=MESHES)\n'
                         'def qa_cells_zz(request):\n'
                         '    return request.param\n'
                         'def test_n(qa_cells_zz):\n'
                         '    if qa_cells_zz % 4:\n'
                         '        raise AssertionError\n'},
  'file': 'test_t.py',
  'new': 'import pytest\n'
         'MESHES = [4, 9]\n'
         '@pytest.fixture(params=MESHES)\n'
         'def qa_cells_zz(request):\n'
         '    return request.param\n'
         'def test_n(qa_cells_zz):\n'
         '    if qa_cells_zz % 4:\n'
         '        raise AssertionError\n',
  'run': 'import pytest\n'
         "RESULT = int(pytest.main(['-q', '-p', 'no:cacheprovider', '--color=no', '-x', 'test_t.py']))\n"}]
