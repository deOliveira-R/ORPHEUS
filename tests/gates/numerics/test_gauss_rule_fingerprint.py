r"""Gate (b): the PLATFORM-INDEPENDENCE WITNESS — a frozen byte fingerprint of the Gauss–Legendre rules.

**What this is.** A PRODUCER pin (``vv-principles``, "a producer's TOLERANCE
pin can never warn its consumers' BIT-IDENTITY pins"): about 100 frozen
artefacts in 15 groups downstream (the SN regression snapshots, the
``--capture-baseline`` arrays, the walk-matvec baselines, ...) are bytes fed by
a Gauss–Legendre rule. When the rule's bytes move, all of them red at once
with messages about fluxes; this gate reds FIRST and by name, saying that the
rule moved and the consumers need no investigation of their own.

**It must be green on every platform.** The digests are those of the
CORRECTLY ROUNDED rules (gate (a),
``tests/gates/numerics/test_gauss_rules_correctly_rounded.py``), computed from
that gate's 80-digit reference and NEVER from production. A correctly rounded
real is unique, so the digests hold on Linux CI (OpenBLAS) and on macOS
(Accelerate) alike, and on any future LAPACK: a red here on one platform and
green on another means production has stopped being correctly rounded there.
``[M]`` 2026-10-01: red at HEAD ``0a5a23fa`` for n in {4, 8, 16, 20, 40, 64}
(Golub–Welsch ``eigh``); n = 2 is green at HEAD because symmetry + mass
imposition fix both numbers exactly.

**The population.** LEGENDRE is the only family production reaches (every
product or folded angular quadrature carries GL(n_mu) through
``gauss_legendre_on_mu``); n = 2, 4, 8, 16, 64 are the brief's orders and
20, 40 are the orders of the SN regression snapshots. Each order has its own
digest, so a red names the order.

**The digest.** ``sha256(nodes.astype('<f8').tobytes() + weights.astype('<f8').tobytes())``,
ascending nodes, little-endian float64 whatever the host byte order.

**Re-capture is forbidden by construction.** The second row recomputes each
digest from gate (a)'s reference; a literal re-pinned from a non-correctly-
rounded production (the move this gate exists to prevent) reds there.
"""

from __future__ import annotations

import hashlib

import numpy as np
import pytest

from orpheus.numerics.generating_measure import LEGENDRE
from orpheus.numerics.quadrature.rules_1d import gauss_legendre_on_mu

from .test_gauss_rules_correctly_rounded import _rule_legendre, reference_rule

_HERE = "tests/gates/numerics/test_gauss_rule_fingerprint.py"
_GATE_A = "tests/gates/numerics/test_gauss_rules_correctly_rounded.py"

#: sha256 of the correctly rounded GL-n rule (gate (a)'s reference, 80 digits,
#: rounded once). Computed 2026-10-01; NEVER regenerate from production.
_FROZEN_GL_DIGEST: dict[int, str] = {
    2: "5a6acf9e0b0b6dc44a60ab078a8be45dbe41c0d303877f4441880b7a434e0408",
    4: "dfa897b41394cc0b2d3187ee491145cea7a9e04c33872acf198715e9a5b2d988",
    8: "cc553b3d933bf2e27d8333f30542767a3d3347ecf6b89ed83264916a29d1fb95",
    16: "279caf968c379e74b42f1cb8d5e800fb553000351b377a66e4968693d5b64713",
    20: "1734f39d54076f331aca615e486976aa5c0675275880f9a0bf5266df936627aa",
    40: "494247905a27ab12116e4a1ccb85c3150e124373229c19415542d193003fc548",
    64: "a6d554a65b2125786f951eff31920ba4f29ce78efcd6c9049df46e803b9121b0",
}


def _digest(nodes: np.ndarray, weights: np.ndarray) -> str:
    payload = np.asarray(nodes).astype("<f8").tobytes() + np.asarray(weights).astype("<f8").tobytes()
    return hashlib.sha256(payload).hexdigest()


@pytest.mark.foundation
@pytest.mark.parametrize("n", sorted(_FROZEN_GL_DIGEST))
def test_frozen_digest_is_the_correctly_rounded_rule(n: int) -> None:
    """The literal is the digest of gate (a)'s reference, not of any production run."""
    nodes, weights = reference_rule(_rule_legendre, n)
    if _digest(nodes, weights) != _FROZEN_GL_DIGEST[n]:
        pytest.fail(
            f"GL-{n}: the frozen digest is not the correctly rounded rule's — "
            "the literal was re-pinned from somewhere else; restore it from the reference"
        )


@pytest.mark.l1
@pytest.mark.rests_on(
    f"{_HERE}::test_frozen_digest_is_the_correctly_rounded_rule",
    f"{_GATE_A}::test_rule_is_the_correctly_rounded_rule",
)
@pytest.mark.parametrize("n", sorted(_FROZEN_GL_DIGEST))
@pytest.mark.parametrize(
    "entry",
    [LEGENDRE.gauss, gauss_legendre_on_mu],
    ids=["LEGENDRE.gauss", "gauss_legendre_on_mu"],
)
def test_gauss_legendre_bytes_are_platform_independent(entry, n: int) -> None:
    """The producer's bytes (and the production entry's) equal the frozen digest.

    A red here means the Gauss–Legendre rule's bytes moved on THIS platform:
    every GL-fed snapshot downstream will red with it, and none of them is a
    separate finding. Re-capturing those snapshots is NOT the remedy; making
    the rule correctly rounded again is.
    """
    rule = entry(n)
    live = _digest(rule.nodes, rule.weights)
    if live != _FROZEN_GL_DIGEST[n]:
        pytest.fail(
            f"GL-{n} bytes moved on this platform ({live[:16]}... != "
            f"{_FROZEN_GL_DIGEST[n][:16]}...): the rule is no longer the correctly "
            "rounded one — every GL-fed snapshot downstream will red with it"
        )
