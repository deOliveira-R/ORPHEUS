r"""F3's ⛔ gate (#426 step 2) — a role subclass of the transfer family carries NO arithmetic.

F3's ruling: *"the kernel tier names the MATHEMATICAL OBJECT; the operator tier
names the TERM of the algebra it realises."*  Its structural half is that
``ScatteringOperator`` / ``N2NOperator`` / ``IsotropicScattering`` /
``IsotropicN2N`` are **thin** subclasses of two shared cores
(:class:`~orpheus.transport.operators.transfer.TransferOperator`,
:class:`~orpheus.transport.operators.isotropic_scattering.IsotropicTransfer`),
whose only content is two class constants — the channel the ONE tier-2 mint
reads and the P0 binding the term lifts — and the role name: NO code.  Without a gate the twin path regrows one override at a time —
exactly the shape the carve removed (`[M]` 2026-09-03, the explorer's
``twin_paths``: ``N2NOperator._interior_space`` ≡
``ScatteringOperator._interior_space`` at **0.94**, ``apply_transpose`` at
**0.85**, ``IsotropicN2N`` member-for-member identical to
``IsotropicScattering``).

**TWO filters, because neither alone is sound** (the 2026-08-19 two-filter
countermeasure):

* **AST** answers *"does this class body define arithmetic?"* — a runtime
  ``dir()`` cannot, because a thin subclass INHERITS every verb and looks
  identical to a fat one.
* **runtime ``__subclasses__``** answers *"is this the whole population?"* — an
  AST census keyed on four names is structurally blind to a fifth subclass
  someone adds later (lessons `L61h`: a class-NAME census misses subclasses;
  `[M]` on this tree an AST base-name census read 3 direct / 5 recursive where
  the runtime answer is 4 / 6).

⚠ This gate is a STRUCTURE claim, not a value claim.  It cannot see a wrong
core; the value rows for that are ``tests/transport/test_transfer_kernel.py``,
``tests/transport/test_material_field.py`` and the tier-2 equivalence family.
Say so, or an audit reads it as coverage of the arithmetic it forbids.
"""
from __future__ import annotations

import ast
import importlib
import inspect
import textwrap

import pytest

pytestmark = pytest.mark.foundation

#: role class → (module path, core class name).  The declared population; the
#: runtime ``__subclasses__`` row below asserts it is also the COMPLETE one.
_ROLES = {
    "ScatteringOperator": ("orpheus.transport.operators.scattering", "TransferOperator"),
    "N2NOperator": ("orpheus.transport.operators.n2n", "TransferOperator"),
    "IsotropicScattering": (
        "orpheus.transport.operators.isotropic_scattering", "IsotropicTransfer",
    ),
    "IsotropicN2N": (
        "orpheus.transport.operators.isotropic_scattering", "IsotropicTransfer",
    ),
}

#: The two cores, and where they live.
_CORES = {
    "TransferOperator": "orpheus.transport.operators.transfer",
    "IsotropicTransfer": "orpheus.transport.operators.isotropic_scattering",
}


def _all_subclasses(cls) -> set:
    """Every subclass at any depth (a direct walk certifies only the first level)."""
    direct = set(cls.__subclasses__())
    return direct.union(*(_all_subclasses(s) for s in direct))


def _resolve(module_path: str, name: str):
    mod = importlib.import_module(module_path)
    obj = getattr(mod, name, None)
    if obj is None:
        pytest.fail(f"{module_path}.{name} does not exist — F3's core/role split is broken")
    return obj


def _class_body(cls) -> list[ast.stmt]:
    r"""The class's own statements, from SOURCE.

    ⚠ ``textwrap.dedent`` on a method's source strips four spaces (`L75`); a
    top-level ``class`` is already at column 0, so the dedent is a no-op here
    and is kept only so a future nested role class does not silently fail to
    parse.  ⚠ ``inspect.getsource`` reads the DECORATED source, so a
    ``@dataclass`` would land in ``decorator_list``, not in the body — the thin
    predicate below is about the body alone, by design: whether the role is a
    dataclass is the CORE's decision.
    """
    tree = ast.parse(textwrap.dedent(inspect.getsource(cls)))
    node = tree.body[0]
    assert isinstance(node, ast.ClassDef)
    return node.body


def _is_thin(stmt: ast.stmt) -> tuple[bool, str]:
    """A role-class body statement is thin iff it is a docstring, a ClassVar
    annotation, or ``pass`` — a role has NO code, not even a classmethod
    (the mints live on the cores since the elegance review's R1: two mint
    bodies of one recipe was the defect this family retired)."""
    if isinstance(stmt, ast.Pass):
        return True, ""
    if isinstance(stmt, ast.Expr) and isinstance(stmt.value, ast.Constant):
        return True, ""  # docstring / bare string
    if isinstance(stmt, ast.AnnAssign):
        ann = ast.unparse(stmt.annotation)
        if ann.startswith("ClassVar"):
            return True, ""
        return False, f"dataclass FIELD `{ast.unparse(stmt.target)}: {ann}`"
    if isinstance(stmt, ast.FunctionDef):
        decos = {ast.unparse(d) for d in stmt.decorator_list}
        return False, f"method `{stmt.name}` (decorators: {sorted(decos) or 'none'})"
    if isinstance(stmt, ast.Assign):
        return False, f"class attribute `{ast.unparse(stmt.targets[0])}`"
    return False, f"statement `{ast.unparse(stmt)[:60]}`"


def test_the_thin_predicate_positive_control():
    r"""`vv` #17 — validate the FILTER before trusting any verdict it gives.

    A broken predicate and a landed carve print the same thing on the
    production rows.  This row runs the predicate on a synthetic thin/fat
    pair covering all five shapes the roles could regrow — an instance
    method, a ``@property``, a ``@classmethod`` (a mint body of its own),
    a dataclass field, a class attribute.
    """
    thin = ast.parse(
        "@dataclass\n"
        "class Thin(Core):\n"
        "    '''doc'''\n"
        "    channel: ClassVar[Callable] = Field.n2n\n"
        "    isotropic_binding: ClassVar[type[Core]] = IsoN2N\n"
    ).body[0]
    fat = ast.parse(
        "@dataclass\n"
        "class Fat(Core):\n"
        "    '''doc'''\n"
        "    frame: HarmonicFrame\n"
        "    _cache = {}\n"
        "    def apply(self, x):\n"
        "        return x\n"
        "    @property\n"
        "    def full_n2n_kernel(self):\n"
        "        return None\n"
        "    @classmethod\n"
        "    def from_solver_data(cls, *, mat_xs, space):\n"
        "        return cls(...)\n"
    ).body[0]
    assert isinstance(thin, ast.ClassDef) and isinstance(fat, ast.ClassDef)

    thin_offenders = [w for ok, w in map(_is_thin, thin.body) if not ok]
    if thin_offenders:
        pytest.fail(
            f"the predicate rejects a legitimately THIN role: {thin_offenders}"
        )
    fat_offenders = [w for ok, w in map(_is_thin, fat.body) if not ok]
    if len(fat_offenders) != 5:
        pytest.fail(
            f"the predicate caught {len(fat_offenders)} of the 5 fat shapes "
            f"(field, class attribute, method, property, classmethod): "
            f"{fat_offenders}"
        )


@pytest.mark.parametrize("role", sorted(_ROLES))
def test_the_role_class_is_a_thin_subclass_of_its_core(role):
    r"""The role class defines NOTHING but ``ClassVar``s.

    The channel constant (which ``Mixture`` channel the core's ONE mint
    reads), the P0-binding constant and the role name are the role's whole
    content; everything with a body — the mints, ``apply``,
    ``apply_transpose``, the faces, the kernel properties — belongs to the
    core, once.
    """
    module_path, core_name = _ROLES[role]
    core = _resolve(_CORES[core_name], core_name)
    cls = _resolve(module_path, role)
    assert issubclass(cls, core), (
        f"{role} is not a subclass of {core_name} — F3's shared core has not "
        f"landed"
    )
    offenders = []
    for stmt in _class_body(cls):
        ok, why = _is_thin(stmt)
        if not ok:
            offenders.append(f"line {stmt.lineno}: {why}")
    if offenders:
        pytest.fail(
            f"{role} carries arithmetic of its own — a twin path regrowing one "
            f"override at a time (F3's ⛔ gate):\n  " + "\n  ".join(offenders)
        )


@pytest.mark.parametrize("core_name", sorted(_CORES))
def test_the_core_has_exactly_the_declared_roles(core_name):
    r"""Runtime population check — the AST filter's blind spot.

    An AST census keyed on four names cannot see a FIFTH subclass; this row
    can.  Assert set EQUALITY (not ``⊆``), so both a new role and a retired one
    are a red with a name.

    ⚠ ``__subclasses__`` only sees classes whose module has been IMPORTED —
    every role module is imported first, so the row does not depend on
    collection order. And it is walked RECURSIVELY: a ``class FastN2N(N2NOperator)``
    carrying an override is the cheapest way to regrow the twin one level
    down, and a direct-only walk is blind to it (the elegance review's S6).
    """
    for module_path, _ in _ROLES.values():
        importlib.import_module(module_path)
    core = _resolve(_CORES[core_name], core_name)
    declared = {r for r, (_, c) in _ROLES.items() if c == core_name}
    got = {c.__name__ for c in _all_subclasses(core)}
    assert got == declared, (
        f"{core_name}'s subclasses are {sorted(got)}; the declared roles are "
        f"{sorted(declared)}.  A new role owes this list an entry AND owes the "
        f"thin-body row above; a retired one owes its removal."
    )
