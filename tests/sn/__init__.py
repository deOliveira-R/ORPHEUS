"""Marks ``tests/sn`` as a package.

⛔ Not cosmetic, and not optional. Without it this directory is a NON-package
containing packages (`tests/sn/operators/`, `tests/sn/regression/`, … all carry
`__init__.py`), so pytest's rootdir walk stops here and every module below is
imported under a TRUNCATED name — `regression.test_dd_regression` rather than
`tests.sn.regression.test_dd_regression`.

That truncation is invisible to the suite (collection and imports work either
way) and fatal to per-test attribution: `coverage`'s `dynamic_context =
test_function` records the module name, the knowledge graph derives its node id
from the file PATH, and the two then disagree for the whole SN tree — the
largest in the project at 3330 tests.

`[M]` 2026-08-19, measured on the wide-capture slices: `exercised_by` was
**0** with `unknown_context` at 696 (`tests/sn/operators`) and 373
(`tests/sn/sweep`), against 0 unknown for `tests/geometry` and
`tests/homogeneous`. With this file the context reads
`tests.sn.regression.test_dd_regression.test_dd_regression` and binds. The
historical `qa_quadrature_cov` run carries the same signature — captured with
contexts, `exercised_by` absent — so this has been silently true of every SN
capture.
"""
