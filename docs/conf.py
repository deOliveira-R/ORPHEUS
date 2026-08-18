# ORPHEUS Documentation — Sphinx Configuration
# ==============================================

import sys
from pathlib import Path

# Add project root to Python path (orpheus package lives there)
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

# -- Project information -----------------------------------------------

project = 'ORPHEUS'
copyright = '2026, Rodrigo de Oliveira'
author = 'Rodrigo de Oliveira'
release = '0.1'

# -- General configuration ---------------------------------------------

extensions = [
    'sphinx.ext.autodoc',       # Pull docstrings from Python code
    'sphinx.ext.mathjax',       # LaTeX math rendering
    'sphinx.ext.viewcode',      # [source] links to highlighted code
    'sphinx.ext.intersphinx',   # Cross-reference external docs (numpy, scipy)
    'sphinx.ext.napoleon',      # Google/NumPy-style docstrings
    'sphinx.ext.todo',          # .. todo:: method-implementer stubs (algebra-of-record)
    'matplotlib.sphinxext.plot_directive',  # .. plot:: for auto-generated figures
    'sphinxcontrib.nexus',                  # Knowledge graph extraction
    'sphinx_design',                        # dropdowns/cards (machine header, gotchas) — #231
    'sphinxcontrib.bibtex',                 # central refs.bib -> citation nodes — #231
]

# sphinxcontrib-bibtex: single source of truth for citations (Zotero is upstream ->
# manual Better-BibTeX export -> refs.bib). Migrating pages replace docutils
# ``.. [Key]`` definitions with ``:cite:`Key```; the remaining pages get a mechanical
# pass afterward (#231 D). Keyed by the existing citation labels for a drop-in swap.
bibtex_bibfiles = ['refs.bib']

# Bibliography labels = the BibTeX key itself ([LewisMiller1984]), so the
# post-migration inline rendering is character-identical to the historical
# docutils bracket citations — the drop-in-swap promise above extends to the
# READER, not just the author. pybtex's stock styles would relabel every
# citation numerically ([1]) and silently break the author-year signal the
# corpus prose leans on.
from pybtex.plugin import register_plugin
from pybtex.style.formatting.plain import Style as PlainStyle
from pybtex.style.labels import BaseLabelStyle


class _KeyLabelStyle(BaseLabelStyle):
    def format_labels(self, sorted_entries):
        for entry in sorted_entries:
            yield entry.key


class _KeyLabelPlainStyle(PlainStyle):
    default_label_style = _KeyLabelStyle


register_plugin('pybtex.style.formatting', 'keylabel', _KeyLabelPlainStyle)
bibtex_default_style = 'keylabel'

# Render `.. todo::` blocks in the output. These are the
# method-implementer's "Archivist expansion needed" stubs per the
# `algebra-of-record` skill's Sphinx-stub vs rich-narrative
# discipline — they MUST be visible so a future session sees which
# theory sections still await rich-narrative expansion.
todo_include_todos = True

templates_path = ['_templates']
# '**/*.inc.rst' — fragments meant only for `.. include::`. Without this they
# ALSO build as standalone documents, so each renders twice and, being in no
# toctree, warns as an orphan. `.. include::` reads them as raw text and is
# unaffected by exclusion.
exclude_patterns = ['_build', '_generated', '**/*.inc.rst', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_js_files = ['sortable.js']
html_css_files = ['sortable.css']

# -- Options for autodoc -----------------------------------------------

# -- Options for Nexus knowledge graph ------------------------------------
#
# Settings live in ``.nexus/config.toml`` at the repository root, NOT here.
# The CLI and the MCP server cannot read conf.py, so a setting declared here
# is invisible to two of the three surfaces that need it — which is how the
# graph's path came to be declared in three places at once.
#
# The rationale for each value moved with it; read that file.

# -- Options for autodoc -----------------------------------------------

autodoc_member_order = 'bysource'
autodoc_typehints = 'description'

# -- Options for napoleon ---------------------------------------------
#
# Render NumPy-style ``Attributes`` sections as ``:ivar:`` info fields
# rather than standalone ``.. attribute::`` directives. Without this
# flag, dataclass fields are double-documented (once by autodoc and
# once by napoleon), causing "duplicate object description" warnings.

napoleon_use_ivar = True

# -- Auto-generate documentation fragments ----------------------------
#
# Each entry runs `python -m <module>` before Sphinx collects sources
# so generated pages always match their live producer — a hand-run
# generator is a frozen-results hazard. The shared failure handler
# uses `sphinx.util.logging` (Sphinx 9.1 removed ``app.warn``; the
# old call would itself AttributeError on the failure path, masking
# the real error — flagged by the archivist at Wave T T.5.2 close-out,
# commit `40e9ccc`).

# Each entry is (module, label, event). Nearly all run at
# `builder-inited`, before Sphinx collects sources, because what they
# write is INCLUDED in the build. One runs at `build-finished` instead:
# it reads the knowledge graph, and nexus writes `graph.db` at
# build-finished, so a `builder-inited` hook would silently read the
# PREVIOUS build's graph.

_GENERATORS = [
    # docs/theory/verification/matrix.rst from the pytest test
    # registry. Closes ORPHEUS issue #79.
    ("tools.verification.generate_matrix", "verification matrix",
     "builder-inited"),
    # One docs/theory/_<pkg>_capability_matrix.inc.rst per package in
    # `orpheus.derivations.continuous` exposing
    # `cases.py:capability_rows()` (auto-discovered; replaced the
    # per-method hooks for peierls_nystrom and fn_method).
    ("tools.verification.generate_capability_matrices",
     "capability matrices", "builder-inited"),
    # The docs/_generated/ fragments — the verification-case table
    # included by docs/theory/verification/summary.rst and the
    # per-method derivation tables — from the continuous-reference
    # registry.
    ("orpheus.derivations.generate_rst", "reference tables",
     "builder-inited"),
    # The vv-principles skill's error-catalogue INDEX, from the graph's
    # `.. error-entry::` nodes and their `catches` edges. Not included
    # by any page — it is `!cat`-injected into the skill — so nothing
    # needs it during the build, which is what lets it run late enough
    # to see a fresh graph.
    ("tools.verification.generate_error_index", "error-catalogue index",
     "build-finished"),
]


def _make_regenerator(module, label):
    def _regenerate(app, *_event_args):
        import subprocess
        from sphinx.util.logging import getLogger
        try:
            subprocess.run(
                [sys.executable, "-m", module],
                cwd=project_root,
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            getLogger(__name__).warning(
                f"{label} regeneration failed: {e.stderr}"
            )
    return _regenerate


def setup(app):
    for module, label, event in _GENERATORS:
        app.connect(event, _make_regenerator(module, label))

# -- Options for mathjax -----------------------------------------------

mathjax3_config = {
    'tex': {
        'macros': {
            # Macros that take a subscript argument to avoid double-subscript errors
            'Sigt': [r'\Sigma_{\mathrm{t},#1}', 1, ''],
            'Sigs': [r'\Sigma_{\mathrm{s},#1}', 1, ''],
            'Siga': [r'\Sigma_{\mathrm{a},#1}', 1, ''],
            'Sigf': [r'\Sigma_{\mathrm{f},#1}', 1, ''],
            'nSigf': [r'\nu\Sigma_{\mathrm{f},#1}', 1, ''],
            'keff': r'k_{\mathrm{eff}}',
            'kinf': r'k_{\infty}',
        },
    },
}

# -- Options for plot directive ----------------------------------------

plot_include_source = True
plot_html_show_source_link = False
plot_formats = ['png']
plot_rcparams = {
    'figure.figsize': (8, 5),
    'font.size': 11,
}

# -- Intersphinx mapping -----------------------------------------------

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
}
