r"""#257 S5 — Spec A: the ``Functional`` category is strictly distinct.

The intrinsic-property gate the user mandates (``feedback_test_intrinsic_properties``):
a math-bearing type ships a test of its DEFINING property. The defining
property of the ``Functional`` category (grand-report §5.6 suffix law) is
that it is **structurally NOT a** :class:`~orpheus.numerics.operator.LinearOperator`.
A Functional maps a field → scalar (or field → scalar-field) via
``evaluate(x) -> float | V``; a LinearOperator maps a field → field via
``apply(x) -> V`` plus ``capabilities`` (and optionally ``solve`` /
``apply_transpose`` / ``H``). The two Protocols share no method.

Teeth (anti-pattern #11 — positive AND negative, and avoid the L11
self-referential trap):

* POSITIVE — a ``ReactionRateFunctional`` (and an estimator wrapper)
  IS a ``Functional``: runtime ``isinstance`` + it has ``evaluate``.
* NEGATIVE — a ``Functional`` is NOT a ``LinearOperator``
  (no ``apply`` / ``capabilities``); CONVERSELY a ``LinearOperator``
  (``RankOneOperator``, ``MultiplicationOperator``) is NOT a
  ``Functional`` (no ``evaluate``). Both directions are required —
  a one-directional test would pass even if ``Functional`` were an
  alias of ``LinearOperator`` with an extra method.

* DISCRIMINATOR (Frame-4 style — prove strict distinctness, not a
  re-label). Two foils answer "is this category a useless alias?":
    - If ``Functional`` were just ``Vector`` (the carrier), a bare
      ``np.ndarray`` would satisfy it. The discriminator asserts a
      Vector that lacks ``evaluate`` is NOT a Functional → the
      category adds the ``evaluate`` suffix law, it is not ``Vector``.
    - If ``Functional`` were just ``LinearOperator`` (rename), then a
      LinearOperator would satisfy it. The discriminator asserts a
      LinearOperator is NOT a Functional → the category is not an
      ``apply``-Protocol in disguise.

vv claim layer (1.5 gate): every row here is a CATEGORY / structural
claim (Protocol membership). Zero eigenvalue claims, zero MMS, zero
convergence-order. No structurally-independent NUMERICAL reference is
needed because no number is asserted — the reference is the Protocol
definitions themselves (positive + negative + discriminator make the
gate non-self-referential).

vv Mode-8: structural assertions route through ``require`` (a function
call, fires under ``python -O``), NEVER a bare ``assert``.

``foundation`` — a software invariant on the type surface, no
theory-page ``:label:``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.functional import InnerProductFunctional
from orpheus.numerics.operator import (
    IdentityOperator,
    LinearOperator,
    RankOneOperator,
    outer,
)
from orpheus.numerics.vector import Vector
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator

from tests.transport._functional_helpers import (
    asymmetric_nu_sigma_f,
    asymmetric_phi,
    build_production_rate_functional,
    cartesian_2d_mesh,
    cross_section_field,
    require,
    require_functional,
)

pytestmark = pytest.mark.foundation


def _make_functional_and_operators():
    """Build one ReactionRateFunctional + two LinearOperators on one mesh."""
    sn = cartesian_2d_mesh(nx=5, ny=3, ng=2)
    nu_sf = asymmetric_nu_sigma_f(ng=2, spatial_shape=sn.spatial_shape)
    func = build_production_rate_functional(cross_section_field(nu_sf, sn))

    # Two concrete LinearOperators to foil against.
    chi = np.ones((2, *sn.spatial_shape))
    rank_one = outer(chi, InnerProductFunctional(nu_sf, axis=0))
    cf = cross_section_field(nu_sf, sn)
    mult = MultiplicationOperator(coefficient=cf, domain=cf.space, codomain=cf.space)
    return func, rank_one, mult, sn


# ═══════════════════════════════════════════════════════════════════════
# POSITIVE — the SUT IS a Functional.
# ═══════════════════════════════════════════════════════════════════════


class TestFunctionalPositive:
    def test_production_rate_functional_is_a_functional(self):
        """``ReactionRateFunctional`` satisfies the ``Functional`` Protocol."""
        Functional = require_functional()
        func, _, _, _ = _make_functional_and_operators()
        require(
            isinstance(func, Functional),
            f"ReactionRateFunctional must satisfy the Functional Protocol; "
            f"isinstance returned False for {type(func).__name__}.",
        )

    def test_production_rate_functional_has_evaluate(self):
        """The suffix law surface: ``evaluate`` is present and callable."""
        func, _, _, _ = _make_functional_and_operators()
        require(
            callable(getattr(func, "evaluate", None)),
            "ReactionRateFunctional must expose a callable `evaluate` "
            "(the §5.6 suffix law field→scalar/scalar-field map).",
        )


# ═══════════════════════════════════════════════════════════════════════
# NEGATIVE — the category does not collapse into LinearOperator.
# ═══════════════════════════════════════════════════════════════════════


class TestFunctionalIsNotLinearOperator:
    def test_functional_is_not_a_linear_operator(self):
        """A Functional must NOT satisfy the LinearOperator contract.

        The defining distinction: a Functional has ``evaluate``, NOT
        ``apply`` + ``capabilities``. The runtime-checkable
        ``LinearOperator`` Protocol requires ``apply`` and
        ``capabilities`` — a genuine Functional has neither, so this
        isinstance MUST be False. If it is True, the category was built
        as a LinearOperator subtype (a category leak).
        """
        func, _, _, _ = _make_functional_and_operators()
        require(
            not isinstance(func, LinearOperator),
            "ReactionRateFunctional must NOT satisfy the LinearOperator "
            "Protocol (a Functional has `evaluate`, not `apply`+`capabilities`). "
            "isinstance returned True — the category leaked into LinearOperator.",
        )

    def test_functional_lacks_apply(self):
        """A Functional has no ``apply`` — it is not an operator."""
        func, _, _, _ = _make_functional_and_operators()
        require(
            not hasattr(func, "apply"),
            "A Functional must NOT carry `apply` (that is the operator "
            "surface). Found `apply` on ReactionRateFunctional.",
        )

    def test_functional_lacks_capabilities(self):
        """A Functional has no ``capabilities`` frozenset."""
        func, _, _, _ = _make_functional_and_operators()
        require(
            not hasattr(func, "capabilities"),
            "A Functional must NOT carry a `capabilities` frozenset (that "
            "is the operator capability-set surface).",
        )


class TestLinearOperatorIsNotFunctional:
    def test_rank_one_operator_is_not_a_functional(self):
        """The fission ``RankOneOperator`` is NOT a Functional (no ``evaluate``).

        This is the converse leg: it proves the Protocols partition the
        objects (the operator the SUT will eventually be composed with
        does not accidentally satisfy the new category).
        """
        Functional = require_functional()
        _, rank_one, _, _ = _make_functional_and_operators()
        require(
            not isinstance(rank_one, Functional),
            "RankOneOperator (a LinearOperator) must NOT satisfy the "
            "Functional Protocol — it has `apply`, not `evaluate`.",
        )

    def test_multiplication_operator_is_not_a_functional(self):
        """The S3 ``MultiplicationOperator`` is NOT a Functional."""
        Functional = require_functional()
        _, _, mult, _ = _make_functional_and_operators()
        require(
            not isinstance(mult, Functional),
            "MultiplicationOperator (a LinearOperator) must NOT satisfy "
            "the Functional Protocol — it has `apply`, not `evaluate`.",
        )

    def test_identity_operator_is_not_a_functional(self):
        """A pure-protocol LinearOperator (no SN fixture) is not a Functional."""
        Functional = require_functional()
        require(
            not isinstance(IdentityOperator(), Functional),
            "IdentityOperator (a LinearOperator) must NOT satisfy the "
            "Functional Protocol.",
        )


# ═══════════════════════════════════════════════════════════════════════
# DISCRIMINATOR — prove strict distinctness (Frame-4 style), not a rename.
# ═══════════════════════════════════════════════════════════════════════


class TestCategoryIsNotAnAlias:
    def test_functional_is_not_just_vector(self):
        """If ``Functional`` were ``Vector``, a bare ndarray would satisfy it.

        A ``np.ndarray`` IS a ``Vector`` (test_vector_protocol pins this)
        but has no ``evaluate``. Asserting it is NOT a ``Functional``
        proves the category adds the suffix law — it is strictly finer
        than ``Vector``, not an alias of it.
        """
        Functional = require_functional()
        arr = np.zeros(3)
        require(
            isinstance(arr, Vector),
            "Precondition: a bare ndarray is a Vector (sanity of the foil).",
        )
        require(
            not isinstance(arr, Functional),
            "A bare ndarray (a Vector with no `evaluate`) must NOT satisfy "
            "Functional — if it does, Functional is an alias of Vector, "
            "not the §5.6 suffix-law category.",
        )

    def test_functional_is_not_just_linear_operator(self):
        """If ``Functional`` were ``LinearOperator``, an operator would satisfy it.

        Already covered object-by-object above; this row states the
        discriminator intent explicitly: a generic LinearOperator
        (``IdentityOperator``) is a clean LinearOperator with no
        ``evaluate``. It must NOT be a Functional — proving the category
        is not a rename of the operator Protocol.
        """
        Functional = require_functional()
        require(
            isinstance(IdentityOperator(), LinearOperator),
            "Precondition: IdentityOperator is a LinearOperator (foil sanity).",
        )
        require(
            not isinstance(IdentityOperator(), Functional),
            "A LinearOperator must NOT satisfy Functional — if it does, "
            "Functional is a rename of LinearOperator, not a new category.",
        )

    def test_evaluate_and_apply_are_disjoint_surfaces(self):
        """The two category surfaces share no method name.

        ``Functional`` speaks ``evaluate``; ``LinearOperator`` speaks
        ``apply`` (+ ``capabilities``). A type satisfying ONE must not be
        forced to carry the OTHER's surface. The SUT has ``evaluate`` and
        not ``apply``; the operators have ``apply`` and not ``evaluate``.
        This is the symmetric statement of the partition.
        """
        func, rank_one, mult, _ = _make_functional_and_operators()
        require(
            hasattr(func, "evaluate") and not hasattr(func, "apply"),
            "Functional surface: has `evaluate`, lacks `apply`.",
        )
        for op in (rank_one, mult):
            require(
                hasattr(op, "apply") and not hasattr(op, "evaluate"),
                f"Operator surface for {type(op).__name__}: has `apply`, "
                f"must lack `evaluate`.",
            )
