#!/usr/bin/env python
r"""Resolve Sphinx Python-domain cross-references against the LIVE namespace.

Why this exists — no Sphinx-side gate can do it
------------------------------------------------

An unresolved Python-domain xref (``:func:`` / ``:class:`` / ``:meth:`` /
``:mod:`` / ``:attr:``) renders as **plain text** and draws no ``-W`` warning.
The usual next thought is "run a nitpicky build" — but ``-n`` cannot help
either, and the reason is structural rather than a config gap:

**Sphinx can only diagnose what it RENDERS.** A docstring in a module that is
never ``automodule``'d is never read, so its roles never become nodes and
nothing can warn about them. Measured 2026-08-03 on this tree: the doc source
carries **45** live ``automodule`` directives against **263** modules in
``orpheus/`` — so ~83 % of the package is invisible to Sphinx at every
severity, and ``tests/`` (which carries many xrefs) is never rendered at all.

This checker sidesteps rendering entirely: it reads prose with ``ast`` /
``tokenize`` (``.py``) or verbatim (``.rst``) and resolves each target by
**importing** it. Render coverage is irrelevant, so it reaches every module,
``tests/``, and the doc source.

Honest scope note — ``.rst`` is a WEAKER case for this tool than ``.py``.
A page under ``docs/`` *is* rendered by definition, so a nitpicky build genuinely
could flag its dead Python-domain roles; the structural argument above applies
only to docstrings. ``.rst`` is covered here anyway for three practical reasons:
one gate instead of two, no ``nitpick_ignore`` curation (this checker only
judges targets it can decide, so it does not flood on ``int``/``ndarray``/
annotation noise), and the same per-site adjudication discipline for both trees.

⚠ **This tool is deliberately NOT the only instrument, and must not become
it.** Its companion is the project's knowledge graph (``nexus
dead-references``, wired into SessionStart). The two ask *different* questions
and — measured 2026-08-10 — **each one's false answer sits in the other's true
column**:

===============================  ==================================  ===============================
instrument                       asks                                structurally blind to
===============================  ==================================  ===============================
this checker (import)            does the name resolve at runtime?    attributes only ``tokenize``
                                                                     can see; see below
nexus (graph)                    does the ref reach an indexed node?  runtime-created attributes
                                                                     (``__getattr__``, decorators)
===============================  ==================================  ===============================

So the **set difference between the two is the triage queue**, and that only
works while they resolve independently. Do NOT "unify" this tool onto the graph
to remove the duplication: sharing a resolution source would give both the same
blind spots and silently destroy the cross-catching. (``vv-principles`` L11 —
structural independence — applied to instruments rather than to references.)

What it decides, and what it declines
-------------------------------------

A target is judged only when it can be *decided*. Two shapes are decidable:

**Absolute** — rooted in something the interpreter can produce: a builtin, or
an importable top-level module (:func:`_root_is_resolvable`). This is *computed*,
not curated. It used to be a hand-written root tuple, and a curated list does not
announce what it is missing: that one omitted ``mpmath`` — the dependency the
whole Peierls reference family rests on — along with ``functools`` /
``dataclasses`` / ``typing`` and every builtin, which between them were **86 of
the 90** roles the tool had been dismissing as foreign.

**Relative, resolved in the prose's own namespace.** ``:meth:`Quadrature.product```
is resolved by Sphinx against the current module; this tool derives that
context from the prose's LOCATION, with zero configuration — the module from
the file path, the class from the AST — and tries the same order Sphinx's
Python domain does (``module.Class.target``, then ``module.target``, then the
target as absolute). This works because **Python's ``import`` statement already
IS the namespace declaration**: ``from orpheus.sn.mesh import SNMesh`` at module
top means ``getattr(module, "SNMesh")`` succeeds, so a docstring saying
``:class:`SNMesh``` names something that module can genuinely see.

Everything else is DECLINED rather than reported, and the declining is what
keeps the gate credible: if the *head* of a relative target is not resolvable
in context either (``:exc:`NotImplementedError``` — a builtin; ``:class:`mpmath.mpf```
— not imported there), the tool says nothing. Measured 2026-08-10: relative
resolution decided **2965** further roles and honestly declined **2925**.

An ``.rst`` page has no module context (this project carries **zero**
``currentmodule`` directives), so a relative target on a page stays
undecidable. That is a real residual hole, not an oversight: ``peierls.rst``
carries ``:func:`peierls_geometry.build_volume_kernel``` naming a module
retired in ``bda76faf``, and no import-based gate without a page-level module
context can see it.

A target is reported ONLY when the failure is a genuinely missing module or
attribute. If the module exists but raises on import, that is reported
separately as a tooling problem, never as a dead reference — the distinction
matters because the two need opposite fixes.

Usage
-----

    .venv/bin/python tools/check_docstring_xrefs.py            # orpheus + tests + docs
    .venv/bin/python tools/check_docstring_xrefs.py orpheus    # narrower

Exit code is 1 when any dead reference is found, so it can gate CI.
"""

from __future__ import annotations

import argparse
import ast
import builtins
import enum
import functools
import importlib
import importlib.util
import inspect
import io
import pathlib
import re
import sys
import textwrap
import tokenize
from collections import defaultdict
from dataclasses import dataclass, field

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent

# Running a script puts the SCRIPT's directory on `sys.path`, not the working
# directory — so this tool started life able to import `orpheus` (pip-installed
# editable) and unable to import `tests` (not installed at all). Every
# `tests.*` target then failed to resolve and was reported DEAD.
#
# `[M]` 2026-08-18: 49 of 49 dead targets in `docs/` were `tests.*`, and the
# sampled ones all EXIST. With the repo root on the path the same scan reports
# **0** in `docs/` and **5** tree-wide. The gate had been red on all three
# roots — issue #302's "71 dead sites" is largely this, not dead references.
#
# ⚠ This is the failure mode the tool's own `DECLINED` outcome exists to
# prevent, arrived at from the other side: not "declined read as alive", but
# "unimportable read as missing". A gate that cries wolf gets ignored, which
# costs more than no gate.
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

#: ``:role:`target``` — the Python-domain roles that render silently when dead.
ROLE_PATTERN = re.compile(r":(func|class|meth|mod|attr|exc|data|obj):`([^`]+)`")

#: A role target is a dotted Python identifier path; nothing else is decidable.
#: Checking this POSITIVELY (rather than blacklisting ``(``/``[``/space) is what
#: makes whitespace-collapsing safe — a signature fragment or a prose phrase
#: fails the pattern instead of being silently mangled into a plausible path.
DOTTED_PATH = re.compile(r"[A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*\Z")

#: ⛔ RETIRED 2026-08-10 — was a hand-curated tuple of "roots decidable by
#: import" (``orpheus``, ``tests``, ``tools``, ``examples``, ``numpy``,
#: ``scipy``, ``sympy``, ``pytest``, ``matplotlib``). It was a hardcoded
#: enumeration standing in for a property that :func:`_root_is_resolvable`
#: computes in one line — and the list silently omitted **mpmath**, the
#: dependency the entire Peierls reference family is built on (28 sites), plus
#: ``functools`` / ``dataclasses`` / ``typing`` and every builtin. `[M]` those
#: omissions were 86 of the 90 roles this tool had been reporting as
#: "foreign, therefore nobody's problem".
#:
#: The lesson is the generic one: a curated list approximating a decidable
#: predicate does not announce what it is missing. The predicate does.


@functools.lru_cache(maxsize=None)
def _root_is_resolvable(root: str) -> bool:
    """True if ``root`` names something the interpreter can produce.

    Two ways: a **builtin** (``ValueError``, ``property``) or an **importable
    top-level module** (``orpheus``, ``mpmath``, ``functools``). Between them
    they are the whole decidable universe, so this replaces the curated root
    list rather than extending it.

    ``find_spec`` locates without importing, so this costs a stat walk per
    distinct root and never executes third-party module bodies.

    ⚠ It answers only about the ROOT, and an import **alias** is not a root:
    ``np`` in ``:func:`np.array_equal``` names nothing importable, which is
    exactly why Sphinx cannot resolve that role either.

    ⛔ But do NOT read that as "so this tool catches it". It does not, and the
    reason is the honest limit of resolving by import. On an ``.rst`` page the
    alias makes the target *relative*, and a page has no namespace, so it is
    DECLINED. In a ``.py`` docstring it is worse: the module very likely did
    ``import numpy as np``, so ``getattr(module, "np")`` **succeeds** and the
    target reads ALIVE — correct for this tool's question ("does the name
    resolve?") and wrong for Sphinx's ("does the role link?"). An alias in a
    role is a real defect that belongs to the *graph* instrument, which resolves
    against registered targets and has no node for ``np`` (`[M]` 2026-08-10:
    766 ``external`` nodes, none of them ``np.*``). Another instance of the
    standing rule that the two instruments' blind spots are complementary.
    """
    if hasattr(builtins, root):
        return True
    try:
        return importlib.util.find_spec(root) is not None
    except (ImportError, ValueError, ModuleNotFoundError, AttributeError, TypeError):
        return False

#: Targets that are legitimately unresolvable — keep EMPTY if at all possible,
#: and always with a comment saying why. An allowlist is where a gate goes to die.
ALLOWLIST: frozenset[str] = frozenset()

#: A ``#:`` run documents the statement BELOW it (Sphinx's attribute-comment
#: form). The prefix is stripped exactly as Sphinx strips it, so a role split
#: across two ``#:`` lines rejoins into one target.
ATTRIBUTE_COMMENT = re.compile(r"^\s*#:\s?(.*)$")


@dataclass(frozen=True, slots=True)
class ProseBlock:
    """A run of prose, plus the namespaces it was WRITTEN IN.

    ``namespaces`` are the dotted prefixes a *relative* role in this block
    resolves against, innermost first. For an ``.rst`` page the tuple is empty
    — a page has no module context — so only absolute targets are decidable
    there.
    """

    lineno: int
    text: str
    namespaces: tuple[str, ...] = field(default=())


def extract_target(raw: str) -> str | None:
    """Normalise a role body to a bare dotted path, or ``None`` if not one.

    Handles the two Sphinx spellings that wrap a target (``~shortened.path``
    and ``Display Text <real.target>``) and the one the *formatter* creates:
    RST joins a paragraph's lines before parsing inline markup, so a role the
    79-column convention wrapped is ONE target that happens to contain a
    newline. Rejoining it is not a normalisation of a defect — it is reading
    the target the way Sphinx reads it.

    ⚠ Only a **line break** is joined, never a bare space, and the distinction
    is load-bearing rather than pedantic. A wrap always contains a newline; a
    space on one line does not. Collapsing all whitespace turns the prose
    phrase ``not a dotted path at all`` into ``notadottedpathatall`` — a
    perfectly valid single identifier that sails through :data:`DOTTED_PATH`,
    so the positive whitelist does NOT catch it. (Measured 2026-08-10 by this
    function's own negative-leg test, after its first docstring claimed the
    whitelist made the collapse safe. It did not.) In outcome such a target
    would be declined for having a non-local head — but relying on that is
    relying on a *second* check to cover a mangling this one should not
    perform, and a two-word phrase that happens to collapse onto a real name
    would be read as a reference to it.
    """
    raw = raw.strip()
    if "<" in raw and raw.endswith(">"):
        raw = raw[raw.index("<") + 1 : -1].strip()
    raw = re.sub(r"\s*\n\s*", "", raw).lstrip("~").lstrip(".")
    return raw if DOTTED_PATH.match(raw) else None


def resolve(dotted: str) -> tuple[bool, str | None]:
    """Resolve ``dotted`` by importing the longest importable prefix and
    ``getattr``-ing the remainder.

    Returns ``(resolved, problem)``. ``problem`` is ``None`` on success,
    ``"import-error:<Type>"`` when a module EXISTS but fails to import (a
    tooling problem, not a dead ref), ``"empty-namespace-package"`` for a
    directory masquerading as a module, or ``"missing"`` when the module or
    attribute genuinely is not there.
    """
    parts = dotted.split(".")
    if hasattr(builtins, parts[0]):
        # `ValueError`, `dict.get`, `property.setter` — the interpreter carries
        # these, no module hosts them, and the import loop below would report
        # every one of them missing.
        obj: object = getattr(builtins, parts[0])
        for attribute in parts[1:]:
            try:
                obj = getattr(obj, attribute)
            except AttributeError:
                return False, "missing"
        return True, None
    for split in range(len(parts), 0, -1):
        module_name = ".".join(parts[:split])
        try:
            obj = importlib.import_module(module_name)
        except ModuleNotFoundError:
            continue  # try a shorter prefix — this segment may be an attribute
        except Exception as exc:  # exists, but importing it blew up
            return False, f"import-error:{type(exc).__name__}"
        remainder = parts[split:]
        if not remainder and _is_empty_namespace_package(obj):
            return False, "empty-namespace-package"
        for attribute in remainder:
            try:
                obj = getattr(obj, attribute)
            except AttributeError:
                if _is_annotated_attribute(obj, attribute):
                    # Declared but not bound on the class — a dataclass field or
                    # a bare ``x: T`` / ``ClassVar[T]``. ``getattr`` on the CLASS
                    # raises (only instances carry it), yet autodoc documents it
                    # and the xref resolves. Treating this as dead would delete
                    # live cross-references — the inverse of the tool's purpose.
                    return True, None
                if isinstance(obj, type) and attribute in _self_attributes(obj):
                    return True, None
                return False, "missing"
        return True, None
    return False, "missing"


def _is_empty_namespace_package(module: object) -> bool:
    """True for a PEP 420 namespace package that contains no Python module.

    ``importlib`` materialises a namespace package from a bare DIRECTORY, so
    ``import orpheus.foo.bar`` succeeds for **any** path under ``orpheus/`` —
    including a directory holding nothing but a ``README.md``. Such a target
    names no code, Sphinx renders nothing for it, and reporting it alive is a
    false ALIVE: it is the one failure mode where importing is *less*
    informative than reading the tree.

    Measured 2026-08-10: the five ``orpheus/derivations/continuous/`` reserved
    placeholders (``spectral_resolvent``, ``pn_method``, ``escape_probability``,
    ``spectral_collocation``, ``spn_method``) are README-only directories that
    all imported clean, and the knowledge graph was right about all five while
    this tool was wrong about all five.

    A namespace package WITH content stays alive — ``tests/sn/`` has no
    ``__init__.py`` yet carries hundreds of modules, so the test is "does any
    Python module live under it", not "is it a namespace package".
    """
    if getattr(module, "__file__", None) is not None:
        return False
    paths = getattr(module, "__path__", None)
    if not paths:
        return False
    return not any(any(pathlib.Path(p).rglob("*.py")) for p in paths)


def _is_annotated_attribute(owner: object, attribute: str) -> bool:
    """True if ``attribute`` is declared as an annotation anywhere in ``owner``'s MRO.

    The MRO walk matters: a dataclass field declared on a base class is absent
    from the subclass's own ``__annotations__``, so an own-class-only check
    would still report inherited fields as dead.
    """
    if not isinstance(owner, type):
        return False
    for klass in getattr(owner, "__mro__", ()):
        if attribute in getattr(klass, "__annotations__", {}):
            return True
    return False


@functools.lru_cache(maxsize=None)
def _self_attributes(klass: type) -> frozenset[str]:
    """Attribute names assigned to ``self`` anywhere in ``klass``'s MRO source.

    The complement of :func:`_is_annotated_attribute`, and the reason it is
    needed: an instance attribute is **invisible to ``getattr`` on the class**,
    so deciding such a target by import alone is a false DEAD — the failure
    direction that makes a tool delete correct cross-references.

    Two shapes, and both are needed because they fail for different reasons:

    ``self.x: T = ...``
        **PEP 526 records this nowhere at runtime.** The annotation is
        function-local and discarded, so the name is absent from ``getattr``
        on the class *and* from every ``__annotations__`` in the MRO, while
        autodoc reads it from source and the xref resolves.
        Worked (2026-08-10): ``orpheus/numerics/face_layout.py:89`` cites
        ``:attr:`SNMesh.bc``` and ``SNMesh`` sets ``self.bc: dict[...] = ...``.
        That citation lives inside a ``#:`` block, so widening the scanned
        surface WITHOUT this function would have made the gate's very first
        output a false red.

    ``self.x = ...``
        Plain assignment, no annotation — invisible to the clause above too.
        `[M]` 2026-08-10: ``SNMesh.mesh`` is set unannotated by
        ``MaterialMesh._init_data`` and reads as a live ``Mesh2D`` on any
        instance, yet ``resolve`` called it *missing*. It was latent rather
        than active only because every citation of it happened to be
        unqualified; the first person to qualify one would have been handed a
        false red on a correct reference. Found while adjudicating #346 W1,
        where a sub-agent's tree check contradicted my class-level probe.

    Any method is scanned, not only ``__init__``: autodoc documents an instance
    attribute wherever it is first assigned.
    """
    names: set[str] = set()
    for owner in getattr(klass, "__mro__", ()):
        try:
            tree = ast.parse(textwrap.dedent(inspect.getsource(owner)))
        except (OSError, TypeError, SyntaxError, IndentationError):
            continue  # C-implemented, dynamically built, or otherwise sourceless
        for node in ast.walk(tree):
            if isinstance(node, ast.AnnAssign):
                targets = [node.target]
            elif isinstance(node, ast.Assign):
                targets = list(node.targets)
            else:
                continue
            for target in targets:
                if (
                    isinstance(target, ast.Attribute)
                    and isinstance(target.value, ast.Name)
                    and target.value.id == "self"
                ):
                    names.add(target.attr)
    return frozenset(names)


def module_name_of(path: pathlib.Path) -> str | None:
    """The dotted module name a ``.py`` file under the repo root implements."""
    if path.suffix != ".py":
        return None
    parts = list(path.relative_to(REPO_ROOT).parts)
    if parts[-1] == "__init__.py":
        parts = parts[:-1]
    else:
        parts[-1] = parts[-1].removesuffix(".py")
    return ".".join(parts) or None


def iter_text_blocks(path: pathlib.Path):
    """Yield a :class:`ProseBlock` for every run of prose worth scanning.

    For ``.py`` that is each module/class/function docstring **and** each
    ``#:`` attribute-comment run, reported at the prose's OWN line so a hit
    resolves to the text rather than to the enclosing ``def``. For ``.rst`` the
    whole file is one block: a doc page is prose end to end.
    """
    if path.suffix != ".py":
        yield ProseBlock(1, path.read_text(encoding="utf-8", errors="replace"))
        return

    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source)
    module = module_name_of(path)
    blocks: list[ProseBlock] = []
    class_scopes: list[tuple[int, int, str]] = []  # (first line, last line, qualname)

    def namespaces_for(qualname: str) -> tuple[str, ...]:
        """Innermost-first, the order Sphinx's Python domain tries them."""
        if module is None:
            return ()
        return (f"{module}.{qualname}", module) if qualname else (module,)

    def collect(node: ast.AST, qualname: str) -> None:
        docstring = _docstring_of(node)
        if docstring is not None:
            lineno, text = docstring
            blocks.append(ProseBlock(lineno, text, namespaces_for(qualname)))
        for child in ast.iter_child_nodes(node):
            if isinstance(child, ast.ClassDef):
                nested = f"{qualname}.{child.name}" if qualname else child.name
                class_scopes.append(
                    (child.lineno, child.end_lineno or child.lineno, nested)
                )
                collect(child, nested)
            elif isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
                # A function does not extend the namespace: a role in a method's
                # docstring resolves against the CLASS, not against the method.
                collect(child, qualname)

    collect(tree, "")

    # ``#:`` runs are COMMENTS, so ``ast`` discards them entirely — they need
    # the token stream. Scoping them needs the AST's class line ranges, which
    # is why this comes after the walk.
    class_scopes.sort(key=lambda scope: scope[1] - scope[0])  # narrowest first

    def namespace_at(lineno: int) -> tuple[str, ...]:
        for start, end, qualname in class_scopes:
            if start <= lineno <= end:
                return namespaces_for(qualname)
        return namespaces_for("")

    blocks.extend(_attribute_comment_blocks(source, namespace_at))
    yield from blocks


def _docstring_of(node: ast.AST) -> tuple[int, str] | None:
    """``(lineno, text)`` of ``node``'s docstring, or ``None`` if it has none.

    Reported at the docstring literal's OWN line so a hit resolves to the prose
    rather than to the enclosing ``def``.
    """
    body = getattr(node, "body", None)
    if not body:
        return None
    first = body[0]
    if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant):
        text = first.value.value
        if isinstance(text, str):
            return first.lineno, text
    return None


def _attribute_comment_blocks(source: str, namespace_at):
    """Yield one :class:`ProseBlock` per maximal run of ``#:`` comment lines.

    A run is joined with newlines and reported at its FIRST line, so a role
    wrapped across two ``#:`` lines is one target (:func:`extract_target`
    collapses the whitespace) at the line a reader would look for it.

    ``tokenize`` rather than a line regex, so a ``#:`` inside a string literal
    is not mistaken for an attribute comment.
    """
    try:
        tokens = list(tokenize.generate_tokens(io.StringIO(source).readline))
    except (tokenize.TokenError, IndentationError, SyntaxError):
        return
    run: list[str] = []
    start = 0
    previous_line = -2
    for token in tokens:
        if token.type != tokenize.COMMENT:
            continue
        match = ATTRIBUTE_COMMENT.match(token.line.rstrip("\n"))
        if match is None:
            continue
        lineno = token.start[0]
        if lineno != previous_line + 1 and run:
            yield ProseBlock(start, "\n".join(run), namespace_at(start))
            run = []
        if not run:
            start = lineno
        run.append(match.group(1))
        previous_line = lineno
    if run:
        yield ProseBlock(start, "\n".join(run), namespace_at(start))


def candidate_paths(
    target: str, namespaces: tuple[str, ...], role: str | None = None
) -> tuple[str, ...]:
    """Every absolute path ``target`` could name, innermost namespace first.

    An absolute target means itself and nothing else. A relative one means
    whatever the prose's namespaces make of it — and with no namespaces (an
    ``.rst`` page) it means nothing this tool can decide.

    "Absolute" is decided by :func:`_root_is_resolvable`, not by a curated list
    (see the retirement note on ``DECIDABLE_ROOTS``).

    ⚠ One guard the curated list did not need. A **single-segment** target that
    resolves *as a module* is almost never a module reference: ``:class:`array```
    means a local class, and stdlib has an ``array`` module, so treating it as
    absolute would report ALIVE for a role Sphinx leaves broken. A bare module
    name is legitimate only under ``:mod:``, so that is the one role allowed to
    take this path. Builtins are exempt — ``:exc:`ValueError``` is single-segment
    and genuinely absolute.
    """
    root = target.split(".")[0]
    bare_module_guess = "." not in target and role != "mod" and not hasattr(builtins, root)
    if _root_is_resolvable(root) and not bare_module_guess:
        return (target,)
    return tuple(f"{namespace}.{target}" for namespace in namespaces)


class Outcome(enum.Enum):
    """The four things that can be true of one role occurrence.

    Spelling these as a closed set rather than as ``(bool, bool, str | None)``
    is what makes "declined" impossible to confuse with "alive" — the single
    most dangerous conflation available to this tool, because both are silent.
    """

    DECLINED = "declined"        # not ours to judge; says nothing, counts nothing
    ALIVE = "alive"
    DEAD = "dead"
    UNIMPORTABLE = "unimportable"  # module exists but blew up — a TOOLING problem


@dataclass(frozen=True, slots=True)
class Judgement:
    """An :class:`Outcome` plus the evidence for it."""

    outcome: Outcome
    spelling: str | None = None  # the path actually judged (innermost candidate)
    problem: str | None = None  # why, when DEAD or UNIMPORTABLE


def judge(
    target: str,
    namespaces: tuple[str, ...],
    lookup=None,
    role: str | None = None,
) -> Judgement:
    """Decide ``target``, read in ``namespaces``, or decline to.

    Extracted from the scan loop so the gate's teeth are provable in-process:
    a subprocess-only decision path can be shown to exit 0, which is exactly
    what a scanner that silently found nothing also does (``vv-principles``
    #17 — verify the instrument on a known positive).

    ``lookup`` defaults to :func:`resolve` and exists so a caller can supply a
    memo across many calls.
    """
    lookup = lookup or resolve
    candidates = candidate_paths(target, namespaces, role)
    if not candidates:
        return Judgement(Outcome.DECLINED)  # relative, no namespace to read it in
    if any(lookup(candidate)[0] for candidate in candidates):
        return Judgement(Outcome.ALIVE, candidates[0])
    problems = [lookup(candidate)[1] for candidate in candidates]
    broken = next((p for p in problems if p and p.startswith("import-error")), None)
    if broken is not None:
        return Judgement(Outcome.UNIMPORTABLE, candidates[0], broken)
    # A relative target whose HEAD is not in the namespace either is a builtin,
    # an external, or prose — not ours to judge. Declining here is what keeps
    # the gate credible; reporting it is how a gate earns the right to be ignored.
    head = target.split(".")[0]
    if not any(lookup(c)[0] for c in candidate_paths(head, namespaces, role)):
        return Judgement(Outcome.DECLINED)
    return Judgement(Outcome.DEAD, candidates[0], problems[0] or "missing")


def source_files(root: pathlib.Path):
    """Every scannable file under ``root``, skipping build output.

    ``docs/_build/**`` holds stale copies of the very pages being checked —
    scanning them reports long-fixed refs as live findings.
    """
    for pattern in ("*.py", "*.rst"):
        for path in sorted(root.rglob(pattern)):
            if "_build" not in path.parts:
                yield path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("roots", nargs="*", default=None,
                        help="directories to scan (default: orpheus tests docs)")
    parser.add_argument("--quiet", action="store_true",
                        help="summary only, no per-site listing")
    args = parser.parse_args()
    roots = args.roots or ["orpheus", "tests", "docs"]

    dead: dict[str, list[tuple[str, int, str]]] = defaultdict(list)
    diagnosis: dict[str, str] = {}
    broken_imports: dict[str, str] = {}
    cache: dict[str, tuple[bool, str | None]] = {}
    n_files = n_roles = n_decided = 0

    def lookup(dotted: str) -> tuple[bool, str | None]:
        if dotted not in cache:
            cache[dotted] = resolve(dotted)
        return cache[dotted]

    for top in roots:
        for path in source_files(REPO_ROOT / top):
            try:
                blocks = list(iter_text_blocks(path))
            except (SyntaxError, UnicodeDecodeError):
                continue
            n_files += 1
            for block in blocks:
                for match in ROLE_PATTERN.finditer(block.text):
                    n_roles += 1
                    role, raw = match.group(1), match.group(2)
                    target = extract_target(raw)
                    if target is None or target in ALLOWLIST:
                        continue
                    verdict = judge(target, block.namespaces, lookup, role)
                    if verdict.outcome is Outcome.DECLINED:
                        continue
                    if verdict.outcome is Outcome.UNIMPORTABLE:
                        assert verdict.spelling and verdict.problem
                        broken_imports[verdict.spelling] = verdict.problem
                        continue
                    n_decided += 1
                    if verdict.outcome is Outcome.ALIVE:
                        continue
                    assert verdict.spelling and verdict.problem
                    diagnosis[verdict.spelling] = verdict.problem
                    rel = str(path.relative_to(REPO_ROOT))
                    lineno = block.lineno + block.text[: match.start()].count("\n")
                    dead[verdict.spelling].append((rel, lineno, role))

    n_sites = sum(len(v) for v in dead.values())
    print(f"files scanned          : {n_files}")
    print(f"xref roles found       : {n_roles}")
    print(f"decidable (gated)      : {n_decided}")
    print(f"DEAD TARGETS           : {len(dead)}  across {n_sites} sites")
    if broken_imports:
        print(f"unimportable (TOOLING, not dead refs): {len(broken_imports)}")
        for target, problem in sorted(broken_imports.items()):
            print(f"    {target} -> {problem}")

    if dead and not args.quiet:
        print()
        for target, sites in sorted(dead.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            print(f"  {target}   [{len(sites)} site(s)]  {diagnosis[target]}")
            for rel, lineno, role in sites:
                print(f"      {rel}:{lineno}  (:{role}:)")

    return 1 if dead else 0


if __name__ == "__main__":
    raise SystemExit(main())
