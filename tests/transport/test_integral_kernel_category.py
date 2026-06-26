r"""#257 S6 — Spec A: the ``IntegralKernelOperator`` category is a strict refinement.

The intrinsic-property gate the user mandates
(``feedback_test_intrinsic_properties``): a math-bearing type ships a test
of its DEFINING property. The defining property of the §5.6 *Kernel*
category (grand-report suffix law) is:

    an ``IntegralKernelOperator`` exposes a ``kernel`` (a nonlocal,
    measure-integrating action); a LOCAL operator (``MultiplicationOperator``,
    ``IdentityOperator``) and a ``Functional`` (field → scalar) do NOT.

Crucially the Kernel is a **refinement of** ``LinearOperator`` (it still
has ``apply`` + ``capabilities``), NOT disjoint from it — UNLIKE the
``Functional`` (S5), which shares no member with ``LinearOperator``. So
the partition is asymmetric: a Kernel IS a LinearOperator, but only SOME
LinearOperators are Kernels (those exposing ``kernel``).

Teeth (anti-pattern #11 — positive AND negative both directions, and
avoid the L11 self-referential trap):

* POSITIVE — ``FissionOperator`` and ``ScatteringOperator`` ARE
  ``IntegralKernelOperator`` s: runtime ``isinstance`` + each has a
  ``kernel`` member returning a ``LinearOperator``. They are ALSO still
  ``LinearOperator`` s (refinement, not disjoint).
* NEGATIVE both directions — a ``MultiplicationOperator`` (S3b,
  local/diagonal) is NOT an ``IntegralKernelOperator`` (no ``kernel``);
  a ``Functional`` / ``ReactionRateFunctional`` (S5, field→scalar) is
  NOT (no ``kernel``, no ``apply``). The Functional NEGATIVE is the
  sharpest: it confirms the new category did NOT accidentally admit the
  field→scalar maps.
* DISCRIMINATOR (Frame-4 style — prove a strict refinement, not a useless
  alias of ``LinearOperator``): an ``IdentityOperator`` / a
  ``MultiplicationOperator`` IS a ``LinearOperator`` but exposes no
  ``kernel`` → it must NOT satisfy ``IntegralKernelOperator``. If it does,
  the Protocol is just ``LinearOperator`` renamed (the refinement adds
  nothing).

The ``@runtime_checkable`` caveat (S5 lesson): a runtime_checkable
isinstance only checks member PRESENCE, so it can be fooled by an
INCIDENTAL same-named attribute. Every negative gate therefore also
asserts the ``kernel``-attribute ABSENCE directly (defense-in-depth),
and the positive gate asserts the ``kernel`` returns a genuine
``LinearOperator`` (not merely that the attribute is present).

vv claim layer (1.5 gate): every row is a CATEGORY / structural claim
(Protocol membership). Zero eigenvalue / MMS / convergence-order claims —
the reference is the Protocol definitions themselves (positive + negative
+ discriminator make the gate non-self-referential, L11).

vv Mode-8: structural assertions route through ``require`` (a function
call, fires under ``python -O``), NEVER a bare ``assert``.

``foundation`` — a software invariant on the type surface, no
theory-page ``:label:``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture, make_mixture
from orpheus.geometry import Mesh2D
from orpheus.numerics.operator import (
    IdentityOperator,
    LinearOperator,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.solver import SNSolver
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator

from tests.transport._integral_kernel_helpers import (
    require,
    require_integral_kernel_operator,
    require_production_rate_property,
    require_scattering_kernel_property,
)
from tests.transport._functional_helpers import (
    asymmetric_nu_sigma_f,
    cartesian_2d_mesh,
    cross_section_field,
)

pytestmark = pytest.mark.foundation


# ═══════════════════════════════════════════════════════════════════════
# Fixtures — a fission operator + an anisotropic scattering operator, both
# carrying real SN data so they EXPOSE the kernel/production_rate surface.
# (placeholder_materials has zero fission/scatter; the category gate needs
# operators that genuinely satisfy the Protocol once landed.)
# ═══════════════════════════════════════════════════════════════════════


def _uniform_2d(nx, ny, delta, mat_map):
    return Mesh2D(
        edges_x=np.linspace(0, nx * delta, nx + 1),
        edges_y=np.linspace(0, ny * delta, ny + 1),
        mat_map=np.asarray(mat_map, dtype=int),
    )


@pytest.fixture
def fission_op():
    """The §5.6 Kernel instance ``F`` — 2G fissile, heterogeneous."""
    fuel = get_mixture("A", "2g")
    mod = get_mixture("B", "2g")
    nx, ny = 6, 4
    mat = np.zeros((nx, ny), dtype=int)
    mat[:3, :] = 2
    mat[3:, :] = 0
    mesh = _uniform_2d(nx, ny, 0.2, mat)
    quad = Quadrature.lebedev(order=17)
    sn_mesh = SNMesh(mesh, quad, {2: fuel, 0: mod})
    return SNSolver(sn_mesh).fission_op


@pytest.fixture
def scattering_op():
    """The §5.6 Kernel instance ``S`` — 2G P1 aniso + asymmetric SigS + n2n."""
    from scipy.sparse import csr_matrix

    p0 = np.array([[0.38, 0.10], [0.05, 0.90]])
    p1 = np.array([[0.02, 0.01], [0.00, 0.04]])
    mix = make_mixture(
        sig_t=np.array([0.5, 1.0]),
        sig_c=np.array([0.01, 0.02]),
        sig_f=np.array([0.01, 0.08]),
        nu=np.array([2.5, 2.5]),
        chi=np.array([1.0, 0.0]),
        sig_s=p0,
    )
    mix.SigS = [csr_matrix(p0), csr_matrix(p1)]
    mix.Sig2 = csr_matrix(np.array([[0.0, 0.03], [0.01, 0.0]]))
    nx, ny = 3, 2
    mesh = _uniform_2d(nx, ny, 0.4, np.zeros((nx, ny), dtype=int))
    quad = Quadrature.lebedev(order=17)
    return SNSolver(
        SNMesh(mesh, quad, {0: mix}), scattering_order=1,
    ).scattering_op


@pytest.fixture
def production_rate_functional():
    """A bare ``ReactionRateFunctional`` foil (field→scalar, no kernel).

    The retired ``ProductionRateFunctional`` is superseded by the
    production-rate :class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`;
    it is still a bare ``Functional`` (no ``kernel``, no ``apply``), so it
    foils the ``IntegralKernelOperator`` negative gates identically.
    """
    from tests.transport._functional_helpers import (
        build_production_rate_functional,
    )

    sn = cartesian_2d_mesh(nx=5, ny=3, ng=2)
    nu_sf = asymmetric_nu_sigma_f(ng=2, spatial_shape=sn.spatial_shape)
    return build_production_rate_functional(cross_section_field(nu_sf, sn))


@pytest.fixture
def multiplication_operator():
    """A bare S3b ``MultiplicationOperator`` foil (local/diagonal, no kernel)."""
    sn = cartesian_2d_mesh(nx=5, ny=3, ng=2)
    nu_sf = asymmetric_nu_sigma_f(ng=2, spatial_shape=sn.spatial_shape)
    return MultiplicationOperator(coefficient=cross_section_field(nu_sf, sn))


# ═══════════════════════════════════════════════════════════════════════
# POSITIVE — the kernel-bearing operators ARE IntegralKernelOperators.
# ═══════════════════════════════════════════════════════════════════════


class TestIntegralKernelPositive:
    def test_fission_operator_is_an_integral_kernel_operator(self, fission_op):
        """``FissionOperator`` satisfies the ``IntegralKernelOperator`` Protocol."""
        IKO = require_integral_kernel_operator()
        # production_rate is the NEW S6 member that completes the Protocol;
        # skip cleanly if not yet landed (PRE-IMPL).
        require_production_rate_property(fission_op)
        require(
            isinstance(fission_op, IKO),
            "FissionOperator must satisfy the IntegralKernelOperator Protocol; "
            f"isinstance returned False for {type(fission_op).__name__}.",
        )

    def test_scattering_operator_is_an_integral_kernel_operator(
        self, scattering_op,
    ):
        """``ScatteringOperator`` satisfies the ``IntegralKernelOperator`` Protocol."""
        IKO = require_integral_kernel_operator()
        require_scattering_kernel_property(scattering_op)
        require(
            isinstance(scattering_op, IKO),
            "ScatteringOperator must satisfy the IntegralKernelOperator "
            f"Protocol; isinstance returned False for "
            f"{type(scattering_op).__name__}.",
        )

    def test_fission_kernel_member_returns_a_linear_operator(self, fission_op):
        """The ``kernel`` member returns a genuine ``LinearOperator``.

        Defense-in-depth against the ``@runtime_checkable`` member-PRESENCE
        loophole: the Protocol's defining payload is that ``kernel`` IS a
        ``LinearOperator`` (an ``apply``-bearing object), not merely that
        the attribute exists.
        """
        kernel = fission_op.kernel
        require(
            isinstance(kernel, LinearOperator),
            "FissionOperator.kernel must return a LinearOperator (an "
            f"apply-bearing operator); got {type(kernel).__name__}.",
        )

    def test_scattering_kernel_member_returns_a_linear_operator(
        self, scattering_op,
    ):
        """``ScatteringOperator.kernel`` returns a genuine ``LinearOperator``."""
        kernel = require_scattering_kernel_property(scattering_op)
        require(
            isinstance(kernel, LinearOperator),
            "ScatteringOperator.kernel must return a LinearOperator; "
            f"got {type(kernel).__name__}.",
        )

    def test_integral_kernel_is_still_a_linear_operator(self, fission_op):
        """A Kernel is a REFINEMENT of LinearOperator — still has ``apply``.

        Unlike the Functional (S5, disjoint from LinearOperator), the
        IntegralKernelOperator is an ``apply``-bearing operator that ALSO
        exposes ``kernel``. This row pins the refinement relationship: the
        SUT remains a LinearOperator (the matvec arms are untouched in S6).
        """
        require_production_rate_property(fission_op)
        require(
            isinstance(fission_op, LinearOperator),
            "FissionOperator (an IntegralKernelOperator) must STILL satisfy "
            "the LinearOperator Protocol — the Kernel is a refinement, not a "
            "disjoint sibling. Found it is no longer a LinearOperator.",
        )
        require(
            hasattr(fission_op, "apply") and hasattr(fission_op, "capabilities"),
            "An IntegralKernelOperator must carry the operator surface "
            "(`apply` + `capabilities`) — it is a LinearOperator refinement.",
        )


# ═══════════════════════════════════════════════════════════════════════
# NEGATIVE — a local operator / a functional is NOT an IntegralKernelOperator.
# Both directions, plus the direct attribute-absence (runtime_checkable
# defense-in-depth).
# ═══════════════════════════════════════════════════════════════════════


class TestLocalOperatorIsNotIntegralKernel:
    def test_multiplication_operator_is_not_an_integral_kernel(
        self, multiplication_operator,
    ):
        """A ``MultiplicationOperator`` (local/diagonal) is NOT a Kernel.

        It is a LINEAR operator (has ``apply``) but its action is LOCAL —
        pointwise multiply, no integration against a measure. It exposes
        no ``kernel``, so it must NOT satisfy the Protocol.
        """
        IKO = require_integral_kernel_operator()
        require(
            not isinstance(multiplication_operator, IKO),
            "MultiplicationOperator (a LOCAL/diagonal operator) must NOT "
            "satisfy IntegralKernelOperator — it has no `kernel` (nonlocal "
            "measure-integration). isinstance returned True: the category "
            "leaked into the local-multiplier algebra.",
        )

    def test_multiplication_operator_lacks_kernel(self, multiplication_operator):
        """Defense-in-depth: the local operator has no ``kernel`` attribute.

        A ``@runtime_checkable`` isinstance only checks member PRESENCE, so
        it could be fooled by an incidental same-named attribute. This row
        asserts the ABSENCE directly — independent of the Protocol's
        isinstance machinery.
        """
        require(
            not hasattr(multiplication_operator, "kernel"),
            "MultiplicationOperator must NOT carry a `kernel` member (that "
            "is the IntegralKernelOperator surface). Found `kernel` on it.",
        )

    def test_identity_operator_is_not_an_integral_kernel(self):
        """A pure-protocol LinearOperator (``IdentityOperator``) is not a Kernel."""
        IKO = require_integral_kernel_operator()
        require(
            not isinstance(IdentityOperator(), IKO),
            "IdentityOperator (a LinearOperator) must NOT satisfy "
            "IntegralKernelOperator — it has no `kernel`.",
        )

    def test_identity_operator_lacks_kernel(self):
        """Defense-in-depth: ``IdentityOperator`` has no ``kernel`` attribute."""
        require(
            not hasattr(IdentityOperator(), "kernel"),
            "IdentityOperator must NOT carry a `kernel` member.",
        )


class TestFunctionalIsNotIntegralKernel:
    def test_production_rate_functional_is_not_an_integral_kernel(
        self, production_rate_functional,
    ):
        """A ``Functional`` (field→scalar) is NOT an ``IntegralKernelOperator``.

        The sharpest negative: the new Kernel category must NOT admit the
        S5 field→scalar maps. A Functional has neither ``kernel`` nor
        ``apply`` — it is a different suffix of the taxonomy.
        """
        IKO = require_integral_kernel_operator()
        require(
            not isinstance(production_rate_functional, IKO),
            "ReactionRateFunctional (a Functional, field→scalar) must NOT "
            "satisfy IntegralKernelOperator — it has no `kernel` and no "
            "`apply`. isinstance returned True: the Kernel category leaked "
            "into the Functional category.",
        )

    def test_production_rate_functional_lacks_kernel_and_apply(
        self, production_rate_functional,
    ):
        """Defense-in-depth: the Functional carries neither surface member."""
        require(
            not hasattr(production_rate_functional, "kernel"),
            "A Functional must NOT carry a `kernel` member.",
        )
        require(
            not hasattr(production_rate_functional, "apply"),
            "A Functional must NOT carry `apply` (the operator surface).",
        )


# ═══════════════════════════════════════════════════════════════════════
# DISCRIMINATOR — prove a strict refinement (Frame-4), not a useless alias
# of LinearOperator.
# ═══════════════════════════════════════════════════════════════════════


class TestKernelIsNotJustLinearOperator:
    def test_a_linear_operator_without_kernel_is_not_an_integral_kernel(self):
        """If ``IntegralKernelOperator`` were ``LinearOperator``, this fails.

        The discriminator: an ``IdentityOperator`` is a clean
        ``LinearOperator`` with NO ``kernel``. If it satisfied
        ``IntegralKernelOperator``, the Protocol would be a rename of
        ``LinearOperator`` (the ``kernel`` refinement adds nothing). It
        must NOT.
        """
        IKO = require_integral_kernel_operator()
        ident = IdentityOperator()
        require(
            isinstance(ident, LinearOperator),
            "Precondition: IdentityOperator is a LinearOperator (foil sanity).",
        )
        require(
            not isinstance(ident, IKO),
            "A LinearOperator that exposes no `kernel` (IdentityOperator) "
            "must NOT satisfy IntegralKernelOperator — else the Protocol is "
            "an alias of LinearOperator, not a refinement adding `kernel`.",
        )

    def test_kernel_refines_but_does_not_collapse_into_functional(
        self, fission_op, production_rate_functional,
    ):
        """The three categories partition asymmetrically.

        * Kernel ⊂ LinearOperator (has ``apply`` + ``kernel``).
        * Functional disjoint from LinearOperator (has ``evaluate``, not
          ``apply``).
        So a Kernel has ``apply`` AND ``kernel`` and NOT ``evaluate``-only;
        a Functional has ``evaluate`` and NOT ``kernel``/``apply``. This
        states the partition symmetrically across the two SUT objects.
        """
        require_production_rate_property(fission_op)
        require(
            hasattr(fission_op, "apply") and hasattr(fission_op, "kernel"),
            "Kernel surface: the IntegralKernelOperator has `apply` + `kernel`.",
        )
        require(
            hasattr(production_rate_functional, "evaluate")
            and not hasattr(production_rate_functional, "kernel"),
            "Functional surface: has `evaluate`, lacks `kernel` — the two "
            "categories do not collapse.",
        )
