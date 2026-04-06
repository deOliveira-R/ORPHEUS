Rebuild the Sphinx documentation and open it in the browser.

1. Run `python -m derivations.generate_rst` from the repository root to regenerate RST fragments from SymPy derivations
2. Run `python -m sphinx -b html docs docs/_build/html` to build the HTML documentation
3. Check for warnings or errors in the build output
4. If the build succeeds, run `open docs/_build/html/index.html` to open in the browser
5. If there are warnings, report them and suggest fixes
6. If there are errors, fix them and rebuild
