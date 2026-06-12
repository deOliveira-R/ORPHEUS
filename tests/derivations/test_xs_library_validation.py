"""Positive + negative gates for :func:`xs_library.validate_all`.

Issue #228 (vv-principles Mode 8): ``validate_all`` used to be a bare
``assert`` in a production module — stripped to a NO-OP under the
canonical ``python -O`` invocation, leaving the hardcoded verification
XS tables with NO live consistency gate (``Mixture.__init__`` does not
re-validate :math:`\\Sigma_t = \\Sigma_c + \\Sigma_f + \\sum_g
\\Sigma_{s,g}`). It now raises :class:`ValueError` explicitly.

Per vv-principles anti-pattern #11, a contract-validation function
needs BOTH directions:

* **positive** — the shipped tables MUST pass (and the call must
  survive ``-O``, which this collected test exercises since the suite
  runs under ``-O``);
* **negative** — a deliberately-inconsistent table MUST raise. Without
  this leg the positive test cannot distinguish "validated" from
  "validation silently skipped" — exactly the Mode-8 failure this
  guards against.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common import xs_library


pytestmark = pytest.mark.foundation


def test_shipped_tables_pass_validation():
    """The hardcoded XS tables satisfy Σ_t = Σ_c + Σ_f + ΣΣ_s."""
    xs_library.validate_all()  # MUST NOT raise


def test_inconsistent_sig_t_raises(monkeypatch):
    """A broken σ_t entry MUST raise ValueError — even under ``-O``."""
    broken = {
        "broken-region": {
            "2eg": {
                "sig_t": np.array([1.0, 2.0]),
                "sig_c": np.array([0.1, 0.2]),
                "sig_f": np.array([0.0, 0.0]),
                # Row sums 0.5 / 0.5 → sig_t should be 0.6 / 0.7, not 1 / 2.
                "sig_s": np.array([[0.3, 0.2], [0.1, 0.4]]),
            },
        },
    }
    monkeypatch.setattr(xs_library, "XS", {**xs_library.XS, **broken})
    with pytest.raises(ValueError, match="XS inconsistency"):
        xs_library.validate_all()
