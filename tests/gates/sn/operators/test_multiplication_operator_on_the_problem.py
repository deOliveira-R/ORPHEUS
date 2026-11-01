r"""Foundation / L0 tests for :class:`MultiplicationOperator` (#257 S3b).

The §5.7 promotion ``M[f]`` of a
:class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`
to a diagonal :class:`~orpheus.numerics.operator.LinearOperator`. Two
test families:

1. **Broadcast oracle** — the generalised N-D
   :class:`~orpheus.numerics.elementary_operators.DiagonalOperator` engine on the
   σ_t leading-ordinate form reduces EXACTLY to the legacy
   ``sigma[None] * psi`` at the VALUES level. The discriminating regime
   is a 2-D mesh carrier ``(N_ord, ng, nx, ny)`` with **nx ≠ ny** and
   **ng = 2** — the only shape where a transposed broadcast axis would
   mis-scale rather than silently agree (variable-swap failure mode #2).
   Primary assertion :func:`numpy.testing.assert_array_equal` (it IS the
   same operation); a ``nulp=1`` fallback is documented only against a
   future re-association (per the test-architect spec, gate item 1).

2. **Multiplier-algebra law-suite** (the intrinsic-properties directive)
   — the laws that make "a coefficient IS a diagonal operator" literally
   true: :math:`M[1]=I`, :math:`M[0]=0` (codomain-aware: collision's zero
   is a SOURCE, not a flux), linearity :math:`M[af+bg]=aM[f]+bM[g]`
   (tested on ≥2G ASYMMETRIC HETEROGENEOUS coefficients — anti-pattern
   #3/#4: a 1G / homogeneous coefficient is blind here), self-adjointness
   :math:`M[f]^*=M[f]`, the spectrum law
   :math:`\mathrm{spec}(M[f])=\mathrm{ess\,range}(f)` (``is_invertible``
   iff ``min|f|>0``), and the homomorphism :math:`M[f]M[g]=M[fg]` tested at
   the VALUES level via the numerics engine on the raw product array
   (``σ·σ'`` has units ``cm⁻²``, so the homomorphism is a numerics-engine
   property — built directly, NOT as a ``CrossSectionField``).

The resolvent-assembly leg (``L + M[σ] → StreamingCollisionOperator``, the
sweep/matvec matching the legacy resolvent and the analytical
``k_∞ = νΣf/Σa`` / streaming-equilibrium references) is verified by the
EXISTING structurally-independent gates that stay green:
``tests/gates/sn/verification/analytical/test_kinf_homogeneous.py`` and
``tests/gates/sn/operators/test_streaming_collision_operator.py`` (the
``test_si_carve_recovers_analytical_kinf`` + the
``..._recovers_q_over_sigma_composite`` streaming-equilibrium leg). This
file adds an L0 bare-engine streaming-equilibrium values check so the
broadcast itself is pinned without a full solve.

``foundation`` — software invariants (the multiplier-algebra laws), not
a theory-page equation claim.

References
----------

* ``.claude/plans/issue_257_coefficient_field_promotion.md`` — S3b.
* ``.claude/agent-memory/test-architect/issue_257_s3_multiplication_operator_verification.md``
  — the 5 gate verdicts + failure-mode exposure.
* Grand Report v3 §5.5–5.7 (the operator promotion).

Moved from tests/gates/transport/test_multiplication_operator.py by layer DAG Phase 4 (2026-11-01): these rows are
the multiplier's behaviour ON THE SN PROBLEM: from_mesh's default ends are the problem's cached full_field_space (by
identity), and the plain L + C on the problem's streaming operator equals the SN sweep atom (ERR-119's catcher on
the LD rows). The class names are kept, so the ids keep their tails; the multiplier's own laws stay in transport,
posed on the raw constructors.
"""

from __future__ import annotations


import numpy as np
import pytest

from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.problem import SNProblem
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from tests.gates.transport._materials import placeholder_materials
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for

#: Exterior conditions of the tags this file's meshes carry (wave 2: a face's role is a deck or an int tag; the condition a tag
#: stands for is declared at the solve, keyword ``boundary_conditions``).
_CONDITIONS = {0: WallCondition.vacuum()}


pytestmark = pytest.mark.foundation


# vv Mode-8: the canonical ORPHEUS invocation is ``python -O``, which
# strips bare ``assert`` to a NO-OP. These structural gates (isinstance /
# capability-set / membership) MUST fire under -O, so they route through
# ``pytest.fail`` (a function call) rather than a bare ``assert`` (a
# false-green tripwire under -O).
def _require(condition: bool, message: str) -> None:
    """A -O-firing assertion (NOT a bare ``assert``)."""
    if not condition:
        pytest.fail(message)


# ═══════════════════════════════════════════════════════════════════════
# Fixtures — a 1-D slab and the discriminating 2-D nx≠ny carrier.
# ═══════════════════════════════════════════════════════════════════════


def _slab_mesh(nx: int = 4, ng: int = 2) -> SNProblem:
    mesh = Mesher(StructuredGeometry.from_homogeneous(1.0, 0)).partition(
        CellsByCount.uniform_width(nx),
    ).mesh
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNProblem(mesh, quad, placeholder_materials(ng=ng), boundary_conditions=conditions_for(mesh, _CONDITIONS))


def _positive_sigma(problem: SNProblem, ng: int = 2, seed: int = 11) -> np.ndarray:
    """Heterogeneous positive σ ``(ng, *spatial)``, bounded away from 0."""
    rng = np.random.default_rng(seed)
    return 0.3 + 0.5 * rng.random((ng, *problem.spatial_shape))


# ═══════════════════════════════════════════════════════════════════════
# 1. Broadcast oracle — the generalised engine ≡ legacy sigma[None]*psi.
# ═══════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════
# 2. The multiplier-algebra law-suite (intrinsic-properties directive).
# ═══════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════
# 3. L0 streaming-equilibrium — the bare-engine values check.
#
# The full resolvent-assembly leg (L + M[σ] → StreamingCollisionOperator, sweep
# matching kinf / streaming-equilibrium) is the EXISTING gate set that
# stays green:
#   tests/gates/sn/verification/analytical/test_kinf_homogeneous.py
#   tests/gates/sn/operators/test_streaming_collision_operator.py::
#     TestStreamingCollisionSolveBridgeRegression::{test_si_carve_recovers_analytical_kinf,
#       test_fixed_source_homogeneous_reflective_recovers_q_over_sigma_composite}
# This pins the BARE multiplier algebra without a solve.
# ═══════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════
# 3. #261 — the composition-guard space metadata.
#
# The gap the collapse of ``CollisionOperator`` folded UP into the base:
# a multiplier may carry the optional composite ``FunctionSpace`` and so
# JOIN the W-D ``OperatorSum`` composition guard (``L + C``). The retired
# subclass added NOTHING the base lacked once the base gained ``space`` —
# these gates pin that the base now does the job (the collapse's premise).
# ═══════════════════════════════════════════════════════════════════════


class TestSpaceMetadataAndGuardJoin:
    """``space`` → ``domain``/``codomain`` → the multiplier joins the guard."""


    def test_from_mesh_defaults_space_from_mesh(self):
        """``from_mesh(σ, problem)`` defaults ``space`` to the mesh's
        ``full_field_space`` — the faithful drop-in for the retired
        ``CollisionOperator(problem, σ)``, which reached the same composite
        space through ``problem.full_field_space`` (#261)."""
        sn = _slab_mesh()
        M = MultiplicationOperator.from_mesh(_positive_sigma(sn), sn)
        _require(
            M.domain is sn.full_field_space,
            "from_mesh must default both ends to mesh.full_field_space",
        )

    def test_from_mesh_chain_short_circuits_on_full_field_space(self):
        """D7 (CS1, path gate; re-derived at CS4b S2) — on a method mesh the
        OPERATOR's space resolves to ``full_field_space`` by identity, never
        to the scalar bulk.

        The chain order is load-bearing: ``SNProblem.bulk_space`` is the
        scalar ``(ng, *spatial)`` bulk, NOT the angular composite, so a
        flipped chain would silently re-space every SN multiplier. Until
        CS4b this was pinned with a poison on ``bulk_space`` ("the chain
        NEVER consults it") — that instrument died when the space-source
        flip made ``bulk_space`` the shared arm of every carrier mint (the
        COEFFICIENT field's space is now correctly sourced from it, and
        the energy-arm single source routes ``angular_bulk_space`` through
        it too). The claim that survives is the OUTCOME, pinned by
        identity on both legs: the operator's space is the composite, and
        the coefficient's is the carrier's scalar bulk.
        """
        sn = _slab_mesh()
        sigma = _positive_sigma(sn)

        M = MultiplicationOperator.from_mesh(sigma, sn)
        _require(
            M.domain is sn.full_field_space,
            "the chain must resolve the composite by identity, first arm",
        )
        _require(
            M.domain is not sn.bulk_space,
            "the operator's ends must never be the scalar bulk",
        )
        _require(
            M.coefficient.space is sn.bulk_space,
            "the coefficient's space is the carrier's cached scalar bulk "
            "(CS4b: one mint per carrier)",
        )


# ═══════════════════════════════════════════════════════════════════════
# 4. The meshless bare-ndarray apply arm (#276 — homogeneous C − K_iso).
#
# The singledispatch arm that lets ``M[σ]`` act on a bare ``(ng, *spatial)``
# group block (NO ordinate axis, NO mesh), so the collision operator
# composes on a mesh-free MaterialMesh and realizes as a dense matrix via
# apply-to-basis-vectors. The load-bearing gate is CROSS-ARM CONSISTENCY:
# the bare arm must agree with the validated FullField-engine arm per
# ordinate (structurally independent — raw ``coeff*x`` vs the engine's
# leading-axis broadcast). The single-cell homogeneous config is
# axis-/group-blind, so these gates run on the discriminating 2-D nx≠ny
# carrier (a wrong axis or a group-drop is observable here, NOT meshless).
# ═══════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════
# G-C1 — the tier-2 classmethod is provably the ctor (vv#28 closure)
# ═══════════════════════════════════════════════════════════════════════


# The LD rows catch ERR-119 (σ broadcast over the leading ordinate axis only);
# the DD rows are its controls, green under the re-dropped bug, so they carry
# no marker.
# The LD rows catch ERR-119 (each reds under the leading-axis-only broadcast);
# the DD rows are the controls (green under it). They are separate test
# functions because the catches marker must be a DECORATOR for the graph to
# read it (a pytest.param mark is invisible to the error index).
_LD_CASES = ("ld_slab_2g", "ld_2d_2g")
_DD_CASES = ("dd_slab_2g", "dd_cart2d_2g")


class TestTheCoefficientBroadcastsOverTheAxesItDoesNotSpan:
    r"""``M[σ]`` is constant over every bulk axis σ does not span (S2, 2026-10-20).

    A cross section spans the energy and spatial axes, so on a
    linear-discontinuous bulk ``(N, ng, *spatial, n_moments)`` the multiply
    is :math:`(σψ)_{n,g,\mathbf c,k} = σ_{g,\mathbf c}\,ψ_{n,g,\mathbf c,k}`:
    the cell-constant σ multiplies each spatial moment alike. Until S2 the
    engine broadcast over the LEADING ordinate axis only, so ``M[σ]`` refused
    an LD carrier outright (``DiagonalOperator(broadcast_axes=(0,)) expects a
    rank-4 coefficient for a 5-D carrier``). Nothing reached it while ``L + C``
    fused into the sweep atom, whose own walk applies σ; once ``+`` stopped
    fusing (numerics_untangling.md ruling 52) the plain ``L + C`` applies C by
    itself. The DD rows are the controls: their bulk has no moment axis.
    """

    @staticmethod
    def _cases():
        from tests.gates.sn.operators.test_g_adjoint_reciprocity import (
            _make_cart2d,
            _make_ld_2d,
            _make_ld_slab,
            _make_slab,
        )

        return {
            "ld_slab_2g": _make_ld_slab,
            "ld_2d_2g": _make_ld_2d,
            "dd_slab_2g": _make_slab,
            "dd_cart2d_2g": _make_cart2d,
        }


    @pytest.mark.catches("ERR-119")
    @pytest.mark.parametrize("case", _LD_CASES)
    def test_the_plain_sum_applies_like_the_sweep_atom(self, case):
        """The LD rows of :meth:`_plain_sum_applies_like_the_sweep_atom`."""
        self._plain_sum_applies_like_the_sweep_atom(case)

    @pytest.mark.parametrize("case", _DD_CASES)
    def test_the_plain_sum_applies_like_the_sweep_atom_dd_control(self, case):
        """The DD controls of :meth:`_plain_sum_applies_like_the_sweep_atom`."""
        self._plain_sum_applies_like_the_sweep_atom(case)

    def _plain_sum_applies_like_the_sweep_atom(self, case):
        r"""``(L + C)ψ``, the plain combination, equals the sweep atom's own
        ``M(σ)ψ`` to ``8ε`` relative. `[M]` 2026-10-20: ≤ 2.3e-16 over three
        seeds on each of these four fixtures; the two differ only by the
        affine-in-σ association ``L ψ + σψ`` against the walk's one fused
        emission. The atom is the independent side: its σ enters through the
        loss representation's walk, never through this multiplier."""
        from orpheus.sn.operators.streaming import StreamingCollisionOperator, StreamingOperator
        from tests.gates.sn.operators.test_g_adjoint_reciprocity import _random_composite

        sn, sig = self._cases()[case](ng=2)
        L = StreamingOperator.pose(sn)
        C = MultiplicationOperator.from_mesh(sig, sn)
        psi = _random_composite(sn, np.random.default_rng(7))
        plain = np.asarray((L + C).apply(psi).interior.values)
        fused = np.asarray(StreamingCollisionOperator(L, C).apply(psi).interior.values)
        scale = float(np.max(np.abs(fused)))
        assert scale > 0.0
        assert float(np.max(np.abs(plain - fused))) <= 8.0 * np.finfo(float).eps * scale


# ── moved here by the layer DAG Phase 2 refresh (specs_2r §5) ──


