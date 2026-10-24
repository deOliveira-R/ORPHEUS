r"""The three method tables of slice S4 together: the domain name, the torus refusal, the state-cone door, the retirement (wave 2, WF2.3).

Spec: ``.claude/plans/wave2_boundary_carve.md`` sections 2.8, 4.2, 4.3 (deltas 8 and 9), 6.4 gates 6 and 7. The per-method tables are
gated in ``tests/gates/cp/test_cp_realizer_table.py``, ``tests/gates/moc/test_moc_boundary_table.py`` and
``tests/gates/mc/test_mc_boundary_table.py``, and cell by cell in the shared pin (``test_admission_matrix.py``); this file holds what
spans two or three of them.

The ladder
==========

::

    test_the_scanner_sees_what_it_must                                 foundation: the AST census reads planted members
        ^
    test_the_string_registries_and_the_typed_law_reads_are_gone        structural: 0 of 8 retired names and 0 BC imports in cp, moc, mc
        ^
    test_moc_names_its_domain_approximation_and_mc_refuses_its_torus   gate 6: MoC A(547) with the name of 4.2; MC F(TOPOLOGY, 547) for every weight
    test_zero_flux_is_state_cone_in_cp_moc_mc                          gate 7: the unregistered door, STATE_CONE, no issue; the diffusion control
    test_zero_flux_is_refused_at_construction_in_cp_moc_mc             gate 7: the same refusal through bind and realize_boundary

Gate 6, and why both halves are one gate
========================================

The same declaration, a specular wall of weight 1 on the solid cylinder's outer tag, is the SAME wall for MoC and MC (the declaration does
not know the method), and the two methods compute it on the same shape, an equal-area square cell, with different carriers: MoC links
specular tracks across the square (a specular return on a different SHAPE: a named DOMAIN approximation, ``Approximated``), MC wraps the
square periodically (a TRANSLATION, a torus: no one-dimensional declaration describes it, ``Refused(TOPOLOGY)``). The gate reads both
cells from their tables on one declaration so that a drift of either one is seen against the other (the contrast is the claim).

Gate 7, the door
================

``ZeroFluxCondition`` is diffusion's method-only condition (``return_fraction == -1``). CP, MoC and MC register no entry for it, and their
tables are on the non-negative-flux cone, so the unregistered door classes it as ``Refused(STATE_CONE)`` (no issue: the cone is theory, not
plumbing; section 4.3 delta 8: three cells flip axis from the measured grid, which counted them plumbing). Diffusion's table, on the signed
cone, admits it, and that is the control: the door is the table's, not the condition's.

First reds at the base tip ``d4da4328``: the retirement row finds the registries and the ``BC`` imports (the census at the base reads the
counts in its message); the two gate rows raise ``ImportError`` at the accessor of their first table. Mutations for the reviewers:
M-moc-exact (MoC's specular admitted as ``Realized(True)``), M-mc-translation-of-the-cylinder (MC's specular admitted on the cylinder),
M-door-cone (a table on the signed cone: ``ZeroFluxCondition`` is then plumbing, not state cone).
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from tests.gates.transport import _wave2_s4_helpers as H

_THIS = "tests/gates/transport/test_wave2_s4_tables.py"
_ROOT = Path(__file__).resolve().parents[3] / "orpheus"
_PACKAGES = ("cp", "moc", "mc")

#: The names the carve retires from the three packages in the same commit as their replacements (section 2.8).
_RETIRED = (
    "BC_REGISTRY", "bc_kind", "_resolve_bc", "_refuse_a_law_cp_drops", "_cp_bc_white", "_cp_bc_vacuum", "_mc_bc_periodic", "_moc_bc_reflective",
    "face_laws", "outer_law",
)


def retired_names(source: str) -> list[str]:
    """The retired identifiers a module defines, assigns, reads or imports, and every import of the typed-law vocabulary (``BC``, and
    anything from ``orpheus.geometry.boundary``). A docstring or comment naming one is not a use and is not found (AST, not text)."""
    found: list[str] = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Name) and node.id in _RETIRED:
            found.append(f"{node.lineno}: name {node.id}")
        elif isinstance(node, ast.Attribute) and node.attr in _RETIRED:
            found.append(f"{node.lineno}: attribute {node.attr}")
        elif isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name in _RETIRED:
            found.append(f"{node.lineno}: def {node.name}")
        elif isinstance(node, ast.ImportFrom):
            names = [a.name for a in node.names]
            if (node.module or "") == "orpheus.geometry" and "BC" in names:
                found.append(f"{node.lineno}: import BC")
            if (node.module or "").startswith("orpheus.geometry.boundary"):
                found.append(f"{node.lineno}: import from {node.module}")
    return found


@pytest.mark.foundation
def test_the_scanner_sees_what_it_must() -> None:
    r"""The positive controls of the census (``instrument-doctrine`` X1): a planted module with one member of each shape (a class
    attribute named ``BC_REGISTRY``, an assignment to ``self.bc_kind``, a retired ``def``, ``from orpheus.geometry import BC``, an
    import from ``orpheus.geometry.boundary``) reads exactly those, a docstring naming ``BC_REGISTRY`` reads none, and a clean module
    reads none. Passes today; reddens on a scanner that reads nothing."""
    planted = (
        "from orpheus.geometry import BC, CoordSystem\n"
        "from orpheus.geometry.boundary import ReflectiveBoundary\n"
        "class M:\n"
        "    BC_REGISTRY = {}\n"
        "    def __init__(self):\n"
        "        self.bc_kind = 'x'\n"
        "def _cp_bc_white(): pass\n"
    )
    found = retired_names(planted)
    assert len(found) == 5, found
    assert retired_names('"""the BC_REGISTRY of the old days"""\nx = 1\n') == []
    assert retired_names("from orpheus.geometry import CoordSystem\n") == []


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_scanner_sees_what_it_must")
def test_the_string_registries_and_the_typed_law_reads_are_gone() -> None:
    r"""The three packages ``orpheus/cp``, ``orpheus/moc`` and ``orpheus/mc`` define, read and import none of :data:`_RETIRED` and import
    neither ``BC`` nor anything of ``orpheus.geometry.boundary``: the registries, their factories, ``bc_kind`` (0 readers) and the reads
    of the retired ``face_laws`` and ``outer_law`` are replaced in the same commit as the tables (section 2.8, the retirement rule).

    Population (X2): every ``.py`` file under the three packages, parsed with ``ast`` (printed in the message); positive control: the
    planted module of the row above. Red input at the base tip: the census reads the registries and the BC imports (see the message).
    Mutation: ``bc_kind`` assigned again in ``MCMesh``.
    """
    files = sorted(p for package in _PACKAGES for p in (_ROOT / package).rglob("*.py"))
    assert len(files) >= 6, f"the census read only {len(files)} files: {files}"
    findings = {str(p.relative_to(_ROOT)): hits for p in files if (hits := retired_names(p.read_text()))}
    assert not findings, f"retired names still used in {len(findings)} of {len(files)} files: {findings}"


@pytest.mark.l1
def test_moc_names_its_domain_approximation_and_mc_refuses_its_torus() -> None:
    r"""Gate 6 of 6.4, read on one declaration. MoC's cell for ``specular(1.0)`` on a solid cylinder's outer site is ``Approximated`` with
    EXACTLY the name ``"an equal-area square cell with reflecting flat faces, not the declared cylinder"`` and issue 547; MoC's
    ``specular(0.5)`` is ``Refused(PLUMBING, #654)``. MC's cell for ``specular(w)`` there is ``Refused(TOPOLOGY, ..., #547)`` for EVERY
    weight of ``0.1, 0.3, 0.5, 0.7, 1.0``, the reason naming the torus, and the contrast is the claim: the same declaration is a named
    approximation under MoC and a refusal under MC. At a Cartesian site MC's specular is plumbing, #654 (the torus is a cylinder's
    problem).

    Red input at the base tip: ``ImportError`` (neither table exists; MoC's registry has no admission and MC's constructs and runs as
    periodic). Mutations: M-moc-exact, M-mc-translation-of-the-cylinder.
    """
    B = H.bc()
    W = B.WallCondition
    where = H.site("RAD_OUTER")
    moc = H.moc_table().admit(W.specular(1.0), where)
    assert isinstance(moc, B.Approximated), f"MoC specular(1.0) is {moc!r}"
    assert H.summary(moc) == ("A", (H.MOC_DOMAIN_NAME,), (547,)), H.summary(moc)
    assert H.summary(H.moc_table().admit(W.specular(0.5), where)) == ("F", "PLUMBING", 654)
    for weight in (0.1, 0.3, 0.5, 0.7, 1.0):
        mc = H.mc_table().admit(W.specular(weight), where)
        assert isinstance(mc, B.Refused) and mc.axis is B.RefusalAxis.TOPOLOGY and mc.issue == 547, f"MC specular({weight}) is {mc!r}"
        assert H.MC_TORUS_FRAGMENT in mc.reason, mc.reason
    assert H.summary(H.mc_table().admit(W.specular(1.0), H.site("CART"))) == ("F", "PLUMBING", 654)


@pytest.mark.l1
def test_zero_flux_is_state_cone_in_cp_moc_mc() -> None:
    r"""Gate 7 of 6.4, the table leg. ``ZeroFluxCondition()`` is ``Refused(STATE_CONE)`` (no issue) at every site of the three tables, through
    the unregistered door (their tables register no entry for it and sit on the non-negative-flux cone). The control: diffusion's table
    (the signed cone, where the condition is registered) admits it, ``Realized(True)``, so the door is the table's and not the condition's.

    Population: CP at 4 sites, MoC at 1, MC at 2, plus the control. Red input at the base tip: ``ImportError`` (the tables do not exist; on
    the HEAD tree the registry refused ``zero_flux`` with ``ValueError``, a plumbing-class message). Mutation: M-door-cone (a table on
    the signed cone reads plumbing, #654).
    """
    zero = H.zero_flux()
    cells = (
        (H.cp_table(), ("CART_XMIN", "CART_XMAX", "RAD_OUTER", "RAD_INNER")),
        (H.moc_table(), ("RAD_OUTER",)),
        (H.mc_table(), ("CART", "RAD_OUTER")),
    )
    for table, sites in cells:
        for label in sites:
            assert H.summary(table.admit(zero, H.site(label))) == ("F", "STATE_CONE", None), f"{table.method} @ {label}"
    from orpheus.diffusion.boundary_realizer import DIFFUSION_REALIZERS

    assert H.summary(DIFFUSION_REALIZERS.admit(zero, H.site("CART"))) == ("R", True), "the control: diffusion registers zero flux"


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_zero_flux_is_state_cone_in_cp_moc_mc")
def test_zero_flux_is_refused_at_construction_in_cp_moc_mc() -> None:
    r"""Gate 7 of 6.4, the construction leg: a mesh whose tag is bound to ``ZeroFluxCondition()`` is refused at the head of CP, MoC and MC
    with ``BoundaryRefused`` carrying the state-cone admission (the same cells, read through ``bind`` and ``realize_boundary``).

    Red input at the base tip: ``ImportError``. Mutations: a head that realizes before it admits; M-door-cone.
    """
    B = H.bc()
    zero = H.zero_flux()
    constructions = (
        lambda: H.cp_head(H.cp_mesh("cylinder"), {0: zero}),
        lambda: H.moc_head(H.moc_mesh1d(), {0: zero}),
        lambda: H.mc_head(H.mc_slab_mesh1d(left=0, right=0), 3.0, {0: zero}),
    )
    for build in constructions:
        with pytest.raises(B.BoundaryRefused) as info:
            build()
        assert info.value.admission.axis is B.RefusalAxis.STATE_CONE
