"""Every binder of ``Withdrawal`` binds the reference package's one class (#405 P2 step 5, R5w.2).

Specified by the test-architect (2026-10-03, ``.claude/plans/reference_p2_spec.md``
§1.5). The row lived in ``tests/gates/reference/test_published.py`` until the
layer-DAG re-filing (Phase 1b, 2026-10-26, ``.claude/plans/layer_dag.md``): it is
a census of one layer's type over its consumers in the layers above it (the
derivations' lock and the withdrawn Peierls family, the test harness's registry
and the root conftest), and ruling 28 homes such a census at the ``tests/gates/``
root. Its binders are import statements (ruling 27): they were module names
handed to ``importlib.import_module`` until then. The row it rests on, the one
definition of the class (R5w.1), stays in the reference tests.
"""

from __future__ import annotations

import pytest

from orpheus.reference import withdrawal
from tests.gates._content_identity_helpers import require

pytestmark = pytest.mark.foundation

_ONE_DEFINITION = (
    "tests/gates/reference/test_published.py::test_r5w_1_withdrawal_is_defined_once_in_the_reference_package"
)


@pytest.mark.rests_on(_ONE_DEFINITION)
def test_r5w_2_every_binding_is_the_one_object() -> None:
    """Every module that binds the name ``Withdrawal`` (the lock, the
    withdrawn Peierls family, the test registry, the root conftest) binds the
    reference package's class, by identity: one type, so a reference
    certificate's ``Withdrawal`` standing and P0's lock agree by construction.
    The lock's own machinery (``withdrawn_generator``, ``GeneratorWithdrawn``)
    stays in derivations."""
    import tests.conftest as root_conftest
    from orpheus.derivations.common import withdrawal as lock
    from orpheus.derivations.continuous import peierls_nystrom
    from tests._harness import registry

    one = withdrawal.Withdrawal
    for binder in (lock, peierls_nystrom, registry, root_conftest):
        bound = getattr(binder, "Withdrawal", None)
        require(bound is one, f"{binder.__name__}.Withdrawal is {bound!r}, not {one!r}")
    require(hasattr(lock, "withdrawn_generator") and hasattr(lock, "GeneratorWithdrawn"), "the lock moved with the type")
