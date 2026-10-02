r"""The two lifts on the continuous sphere, in SymPy (#405 P1 step 6, S6.10).

Branch 1 of the lift a per-region table takes into phase space
(:mod:`orpheus.derivations.common.angular_measure`): a source rate through the
section E (``R(E Q) = Q``), a detector through the retraction's adjoint R†
(``⟨R†Σ, ψ⟩ = ⟨Σ, Rψ⟩``). The mass of the measure, 4π, is derived by
integration and never typed; the AST row holds the module to that.
"""

from __future__ import annotations

import ast
import inspect

import pytest
import sympy as sp

from orpheus.derivations.common import angular_measure
from orpheus.derivations.common.angular_measure import derive_adjoint_identity, derive_section_identity, mass

pytestmark = pytest.mark.foundation


def test_s6_10_a_the_section_keeps_the_rate() -> None:
    result = derive_section_identity()
    if not result["pass"]:
        pytest.fail(f"R(E Q) != Q: residuals {result['residuals']}")


def test_s6_10_b_the_pullback_is_the_adjoint_of_the_retraction() -> None:
    """And the wrong arrow is measured, not assumed: lifting the detector by
    the section E costs a factor 1/m = 1/(4π) in the response."""
    result = derive_adjoint_identity()
    if not result["pass"]:
        pytest.fail(f"⟨R†Σ, ψ⟩ = {result['lhs']} != ⟨Σ, Rψ⟩ = {result['rhs']}")
    if sp.simplify(result["ratio_if_lifted_by_the_section"] - 1 / mass()) != 0:
        pytest.fail(f"the section-lifted detector's ratio is {result['ratio_if_lifted_by_the_section']}, not 1/m")


def test_s6_10_the_mass_is_derived_and_is_the_sphere() -> None:
    """m = ∫ dΩ over the declared domain is 4π (the comparison constant is
    written only here, in the test, as the independent value)."""
    if sp.simplify(mass() - 4 * sp.pi) != 0:
        pytest.fail(f"the derived mass is {mass()}, not the sphere's 4π")


def test_s6_10_c_pi_is_written_only_in_the_domain() -> None:
    """The module spells ``pi`` exactly once, inside the ``SPHERE`` domain (the
    azimuth's range); a typed mass (``4*pi``, ``1/(4*pi)``) would add a site."""
    tree = ast.parse(inspect.getsource(angular_measure))
    sites = [node for node in ast.walk(tree) if isinstance(node, ast.Attribute) and node.attr == "pi"]
    if len(sites) == 0:
        pytest.fail("activation: the AST pass found no `pi` at all, not even in SPHERE")
    sphere = next(
        node for node in ast.walk(tree)
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "SPHERE" for t in node.targets)
    )
    inside = {id(node) for node in ast.walk(sphere)}
    outside = [node.lineno for node in sites if id(node) not in inside]
    if outside:
        pytest.fail(f"`pi` written outside the SPHERE domain at lines {outside}")
