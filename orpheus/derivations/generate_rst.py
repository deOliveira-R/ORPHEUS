"""Generate RST fragments for Sphinx documentation.

Imports all verification cases from the derivation modules and writes
RST files with LaTeX equations and results tables to docs/_generated/.

Run before sphinx-build:
    python -m derivations.generate_rst
"""

# ---------------------------------------------------------------------------
# DERIVATION-FRAGMENT-GENERATOR (#231) -- FIRST INSTANCE SHIPPED (P6 #281)
#
# The documentation methodology (.claude/plans/sn_book_architecture.md,
# "Doc-generation methodology") is that the algebra-of-record GENERATES the
# documentation: the doc's math is rendered straight from the verified SymPy so
# it can never drift -- NOT hand-transcribed, and NEVER merely referenced.
#
# First shipped instance (user-ruled at the P6 open, 2026-07-26):
# generate_homogenization_collapse_rules() renders the adjoint-weighted
# collapse-rule display equations from
# orpheus.derivations.common.homogenization.collapse_rules() -- each one
# PROOF-WELDED inside its theorem to the builder the production code mirrors
# (display == builder at the concrete grids), so the rendered math, the
# verified math, and the implemented math cannot drift apart. The fragment is
# ``.. include::``-ed by docs/theory/foundations/frame.rst.
#
# REMAINING SCOPE (still manual-leveraged): the discrete.sn.* derivation
# fragments (2x2 cell systems, Schur eliminations, cell-balance solves) --
# extend with the same pattern as chapters stabilise. Worked manual example:
# docs/theory/foundations/discretization.rst's LD Schur, hand-expressed from
# derive_d1_reduction_to_production. Greppable tag: DERIVATION-FRAGMENT-GENERATOR.
# ---------------------------------------------------------------------------

from __future__ import annotations

from pathlib import Path

from .reference_values import all_cases, by_method


OUTPUT_DIR = Path(__file__).resolve().parent.parent.parent / "docs" / "_generated"

# Display names for method and geometry codes
_METHOD_LABELS = {
    "homo": "Homogeneous",
    "cp": "Collision Probability",
    "sn": "Discrete Ordinates",
    "moc": "Method of Characteristics",
    "mc": "Monte Carlo",
    "dif": "Diffusion",
}

_GEOM_LABELS = {
    "--": "Homogeneous (0D)",
    "slab": "Slab 1D",
    "cyl1D": "Cylindrical 1D",
    "sph1D": "Spherical 1D",
}

# Sort order for method-grouped display
_METHOD_ORDER = ["homo", "sn", "cp", "moc", "mc", "dif"]


def _write(name: str, content: str) -> None:
    """Write an RST fragment to the output directory."""
    path = OUTPUT_DIR / name
    path.write_text(content)


def _method_fragment(method: str, geometry: str | None = None) -> str:
    """Build an RST fragment for all cases of a given method (and optional geometry)."""
    cases = by_method(method)
    if geometry is not None:
        cases = [c for c in cases if c.geometry == geometry]
    parts = []
    for c in sorted(cases, key=lambda c: c.n_groups):
        parts.append(f"**{c.name}**: {c.description}")
        parts.append("")
        parts.append(c.latex)
        parts.append("")
    return "\n".join(parts) + "\n"


def generate_verification_table() -> None:
    """Master table of all verification cases as a sortable HTML table."""
    def sort_key(c):
        m = _METHOD_ORDER.index(c.method) if c.method in _METHOD_ORDER else 99
        return (m, c.n_groups, c.geometry)

    cases = sorted(all_cases(), key=sort_key)

    rows = []
    for c in cases:
        method_label = _METHOD_LABELS.get(c.method, c.method)
        geom_label = _GEOM_LABELS.get(c.geometry, c.geometry)
        tol_label = c.tolerance or "&mdash;"
        rows.append(
            f"      <tr>"
            f"<td><code>{c.name}</code></td>"
            f"<td>{method_label}</td>"
            f"<td>{geom_label}</td>"
            f"<td>{c.n_groups}</td>"
            f"<td>{c.n_regions}</td>"
            f"<td>{tol_label}</td>"
            f'<td>{c.k_inf:.10f}</td>'
            f"</tr>"
        )

    html = "\n".join([
        '.. raw:: html',
        '',
        '   <table class="sortable docutils align-default">',
        '     <thead>',
        '       <tr>',
        '         <th>Name</th>',
        '         <th>Method</th>',
        '         <th>Geometry</th>',
        '         <th>Groups</th>',
        '         <th>Regions</th>',
        '         <th>Tolerance</th>',
        '         <th><em>k</em><sub>&infin;</sub></th>',
        '       </tr>',
        '     </thead>',
        '     <tbody>',
        *rows,
        '     </tbody>',
        '   </table>',
    ])

    _write("verification_table.rst", html + "\n")


_COLLAPSE_RULE_CONTEXT = {
    # key -> (display title, one-line context; the theorem name is the proof)
    "vector": (
        "Vector channels (T1)",
        r"The response vectors :math:`\Sigma_c,\Sigma_L,\Sigma_f` collapse "
        r"with the bilinear PAIR weight :math:`\varphi^*\!\odot\varphi` — "
        r"the unique rule zeroing the P0 reaction worth "
        r"(``derive_vector_channel_rule``).",
    ),
    "angular_sigma_t": (
        "The collision channel (T1b)",
        r"The pencil's collision term pairs ANGULARLY: :math:`\Sigma_t` "
        r"collapses with :math:`\rho_{i,g} = \sum_n w_n\psi^*\psi` — unique, "
        r"and identical to the scalar pair on isotropic shapes "
        r"(``derive_angular_sigma_t_rule``; user-ruled implemented).",
    ),
    "matrix_per_pair": (
        "Matrix channels (T2)",
        r"Every :math:`(g'\!\to\!g)` entry collapses with the per-pair "
        r"sink-adjoint × source-flux weight — the source-product broadcast "
        r"provably does NOT zero the worth "
        r"(``derive_matrix_channel_rule``; B&G (6.136) is the classical "
        r"statement).",
    ),
    "fission_nsf_mixed_fold": (
        "Fission production (T3)",
        r"The mixed-fold rule — numerator folded by the fine emission "
        r"importance :math:`\iota_i=\sum_g\varphi^*_{i,g}\chi_{i,g}`, "
        r"denominator by the collapsed :math:`\tilde\iota_i` — zeroes the "
        r"TOTAL fission worth for ANY simplex :math:`\chi_R` "
        r"(``derive_fission_factored_rule``).",
    ),
    "fission_chi_canonical": (
        "Fission emission (T3)",
        r"The canonical :math:`\chi_R` — the adjoint-weighted-emission "
        r"convex average, weights :math:`q_i = V_i\,\iota_i\,p_i` — a "
        r"simplex by construction, production-weighted at flat "
        r":math:`\varphi^*`.",
    ),
}


def generate_homogenization_collapse_rules() -> None:
    """Render the P6 adjoint-weighted collapse rules from the algebra of record.

    The display equations come from
    :func:`orpheus.derivations.common.homogenization.collapse_rules` — each
    proof-welded to its builder inside the theorem suite — so the fragment
    regenerates on every build and can never drift from the verified rules.
    """
    import sympy as sp

    from orpheus.derivations.common.homogenization import collapse_rules

    parts = [
        ".. GENERATED by orpheus.derivations.generate_rst — do NOT edit.",
        ".. Source of truth: orpheus/derivations/common/homogenization.py",
        ".. (each display equation is proof-welded to its builder in the",
        "..  theorem suite; the pin is tests/derivations/",
        "..  test_homogenization_rules.py).",
        "",
    ]
    rules = collapse_rules()
    for key, eq in rules.items():
        title, context = _COLLAPSE_RULE_CONTEXT[key]
        parts.append(f".. rubric:: {title}")
        parts.append("")
        parts.append(context)
        parts.append("")
        parts.append(".. math::")
        parts.append("")
        parts.append(f"   {sp.latex(eq)}")
        parts.append("")
    _write("homogenization_collapse_rules.inc.rst", "\n".join(parts))


def main() -> None:
    """Generate all RST fragments."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    generate_verification_table()
    generate_homogenization_collapse_rules()
    _write("homogeneous_derivation.rst", _method_fragment("homo"))
    _write("sn_derivation.rst", _method_fragment("sn"))
    _write("cp_slab_derivation.rst", _method_fragment("cp", "slab"))
    _write("cp_cylinder_derivation.rst", _method_fragment("cp", "cyl1D"))
    _write("moc_derivation.rst", _method_fragment("moc"))
    _write("mc_derivation.rst", _method_fragment("mc"))
    _write("diffusion_derivation.rst", _method_fragment("dif"))

    n_files = len(list(OUTPUT_DIR.glob("*.rst")))
    print(f"Generated {n_files} RST fragments in {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
