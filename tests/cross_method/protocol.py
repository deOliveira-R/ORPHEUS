r"""Cross-method test protocol — shared schema for cross-solver
regression cases.

Three pieces:

1. :class:`ScalarResult` — what every adapter returns.
2. :class:`SolverAdapter` — the ``solve(case) -> ScalarResult``
   protocol every solver wraps to.
3. :class:`CrossMethodCase` — case spec carrying registry case
   reference, per-solver tolerances, claim layer, and pillar.

The protocol is intentionally minimal. It does NOT replace
:class:`~orpheus.derivations.continuous.sood_registry.la13511.La13511Case`
— it CONSUMES it. Each :class:`CrossMethodCase` references an
existing registry case (or carries inline truth for cases that
don't yet have a registry entry, e.g. NM 1980 reflected slab) and
adds:

* the cross-method-specific scalar tag (``"k_eff"``,
  ``"a_critical_mfp"``, ``"R_critical_mfp"``,
  ``"tau_critical_mfp"``);
* the per-solver tolerance map;
* the claim-layer / pillar tags from
  :doc:`/skills/vv-principles`.

When a second project (e.g. spectral_resolvent or singular_eigenfunction)
reaches the same point, this module gets promoted into
:mod:`orpheus.derivations.registry` per the wave3 plan
(``.claude/plans/wave3/architecture.md``). Until then it lives in
the test tree.

References
----------

* :doc:`/testing/cross_method` — architecture page.
* ``.claude/scratch/cross_method_test_protocol_assessment.md`` — the
  Phase-1 architectural assessment that produced this protocol.
* :doc:`/skills/vv-principles` — claim-layer + pillar discipline.
* :doc:`/skills/algebra-of-record` — structural-independence ladder.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal, Mapping, Protocol, runtime_checkable

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.geometry_template import MeshTemplate

ClaimLayer = Literal["eigenvalue", "flux-shape", "convergence-order"]
"""The claim layer per :doc:`/skills/vv-principles` §"Hierarchical
claim taxonomy". Eigenvalue and critical-dimension claims live here;
flux-shape and convergence-order are reserved for tests that pin
those quantities specifically.
"""

Pillar = Literal["closed-form", "MMS", "semi-analytical", "ancillary"]
"""The verification pillar per :doc:`/skills/vv-principles` §"The
three pillars of verification". ``"ancillary"`` covers L4 cross-
implementation references that are NOT pillars but are useful as
cross-checks for a chain whose L1 backing comes from elsewhere.
"""

ScalarTag = Literal[
    "k_eff",
    "k_inf",
    "a_critical_mfp",
    "R_critical_mfp",
    "tau_critical_mfp",
]
"""The scalar quantity each adapter returns. The protocol's
agreement check compares two :class:`ScalarResult` instances with
the same ``tag``.
"""


@dataclass(frozen=True)
class ScalarResult:
    """The per-solver scalar output of a cross-method comparison.

    Attributes
    ----------
    tag : ScalarTag
        Which scalar this result represents. The cross-method
        agreement check requires both results carry the same tag.
    value : float
        The scalar value (e.g. ``k_eff = 1.0``,
        ``a_critical_mfp = 0.93772556``).
    solver_name : str
        The adapter's ``name`` field — for error messages.
    metadata : Mapping[str, Any]
        Solver-specific diagnostics: ``iterations``, ``converged``,
        ``determinant_residual``, etc. Not used by the agreement
        check; surfaced in failure messages.
    """

    tag: ScalarTag
    value: float
    solver_name: str
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def agrees_with(self, other: ScalarResult, tol: float) -> bool:
        """Absolute-tolerance agreement; raises on tag mismatch."""
        if self.tag != other.tag:
            raise ValueError(
                f"ScalarResult tag mismatch: {self.tag!r} vs {other.tag!r}"
            )
        return abs(self.value - other.value) <= tol

    def diff(self, other: ScalarResult) -> float:
        """Absolute difference; raises on tag mismatch."""
        if self.tag != other.tag:
            raise ValueError(
                f"ScalarResult tag mismatch: {self.tag!r} vs {other.tag!r}"
            )
        return abs(self.value - other.value)


@runtime_checkable
class SolverAdapter(Protocol):
    """Protocol every solver wraps to.

    Implementations live in :mod:`tests.cross_method.adapters`. Each
    adapter handles unit conversions and parameter selection
    internally; the test code only sees ``solve(case) -> ScalarResult``.

    Attributes
    ----------
    name : str
        Stable identifier for tolerance lookups in
        :attr:`CrossMethodCase.tolerances`. Convention: lowercase
        snake_case, e.g. ``"fn_slab"``, ``"trajectory_resolvent_slab"``,
        ``"fn_reflected_slab"``.
    method : str
        The continuous method family (``"fn_method"``,
        ``"trajectory_resolvent"``, ``"singular_eigenfunction"``,
        ...). Used by the agreement-matrix renderer to group adapters.
    geometry : str
        Geometry the adapter handles (``"slab"``, ``"sphere-1d"``,
        ``"reflected-slab"``, ``"sphere-mr-fixed-source"``, ...).
    """

    name: str
    method: str
    geometry: str

    def solve(self, case: CrossMethodCase) -> ScalarResult:
        """Solve the case and return the scalar result.

        The adapter is responsible for:
        - extracting the right cross sections from
          ``case.registry_case.materials`` (or the inline
          ``case.materials`` if no registry case);
        - selecting the right solver parameters (n_modes for fn_method,
          n_r/n_mu/n_traj for trajectory_resolvent) for the
          requested ``case.tolerances`` floor;
        - performing any unit conversions (mfp ↔ cm, half-thickness
          ↔ full slab, etc.);
        - returning a :class:`ScalarResult` with the right ``tag``.

        Adapters MUST NOT consult ``case.tolerances`` to alter the
        answer — only to scale internal numerical parameters (more
        modes for a tighter requested tolerance).
        """
        ...


@dataclass(frozen=True)
class CrossMethodCase:
    """A single physical case carrying enough metadata for any
    cross-method comparison gate.

    A case provides its physical parameters either via a
    ``registry_case`` (the registry-backed path — typically a
    :class:`~orpheus.derivations.continuous.sood_registry.la13511.La13511Case`)
    or via inline ``materials`` + ``mesh_template`` fields (for cases
    NOT in the registry — closed-sphere k_inf, MMS, custom
    configurations). Exactly ONE of these two paths must be populated;
    :meth:`__post_init__` enforces the invariant.

    Attributes
    ----------
    case_id : str
        Stable identifier (typically matches the registry case ID,
        with extra suffix for parametric variants). Used as the
        pytest test ID.
    description : str
        One-line human-readable description.
    registry_case : Any | None
        Reference into the upstream case registry (e.g. a
        :class:`~orpheus.derivations.continuous.sood_registry.la13511.La13511Case`).
        ``None`` for cases that don't yet have a registry entry —
        those carry inline ``materials`` / ``mesh_template`` instead.
    materials : dict[int, Mixture] | None
        Optional inline cross sections for cases without a registry
        entry, keyed by integer material ID (matches the production-
        solver convention). MUST be paired with a non-None
        :attr:`mesh_template`. Defaults to ``None`` for registry-
        backed cases.
    mesh_template : MeshTemplate | None
        Optional inline geometry recipe for cases without a registry
        entry. MUST be paired with non-None :attr:`materials`.
        Defaults to ``None`` for registry-backed cases. The
        :class:`~orpheus.derivations.common.geometry_template.MeshTemplate`
        carries ``(geometry_kind, critical_dimension_{mfp,cm},
        n_groups, mat_id, bc_left, bc_right)`` — enough for adapters
        to dispatch all three families uniformly.
    geometry : str
        Geometry tag matching the adapter's geometry field
        (``"slab"``, ``"sphere-1d"``, ``"reflected-slab"``, ...).
    truth_tag : ScalarTag
        Which scalar this case's truth value represents.
    truth_value : float
        The reference value to compare against. For bare-critical
        cases this is typically ``a_critical_mfp`` or
        ``R_critical_mfp``; for k-eigenvalue cases it is ``1.0``
        (criticality) or ``k_inf`` (closed configuration).
    truth_source : str
        Primary literature citation (e.g. "KLL 1974 Table I via
        Sood LA-13511 Table 4", "Grandjean-Siewert 1979 Table XI",
        "NM 1980 Table 2"). Carried verbatim into error messages
        so a failing test names the source it disagrees with.
    pillar : Pillar
        Verification pillar of the truth value. ``"closed-form"``
        for KLL/Sood transcendental dispersion roots;
        ``"semi-analytical"`` for arbitrary-precision integral
        evaluations (e.g. PS-1982 Peierls integral); ``"MMS"`` for
        manufactured solutions. ``"ancillary"`` is reserved.
    claim_layer : ClaimLayer
        Per :doc:`/skills/vv-principles` §"Hierarchical claim
        taxonomy". Cross-method bare-critical cases are eigenvalue
        claims (``k_eff = 1.0`` or ``critical_dimension``).
    tolerances : Mapping[str, float]
        Per-adapter absolute tolerance: ``{adapter.name: tol}``. An
        adapter not listed here will fail the agreement check by
        default — case authors must explicitly opt every adapter
        in by providing a tolerance.
    notes : str
        Free-form notes (convention drift, known issues, references
        to specific lines of the citation, etc.). Surfaced in error
        messages.
    """

    case_id: str
    description: str
    registry_case: Any | None
    geometry: str
    truth_tag: ScalarTag
    truth_value: float
    truth_source: str
    pillar: Pillar
    claim_layer: ClaimLayer
    tolerances: Mapping[str, float]
    notes: str = ""
    materials: dict[int, Mixture] | None = None
    mesh_template: MeshTemplate | None = None

    def __post_init__(self) -> None:
        """Validate the provisioning paths for materials + geometry.

        A case may carry physical parameters via three paths,
        depending on its kind:

        1. **Registry-backed**: ``registry_case`` is set, both
           ``materials`` and ``mesh_template`` are ``None``. The
           adapter reads XS + geometry off
           ``case.registry_case.materials`` /
           ``case.registry_case.mesh_template``.
        2. **Inline**: ``materials`` and ``mesh_template`` are both
           set, ``registry_case`` is ``None``. The adapter reads
           XS + geometry directly off the case.
        3. **Notes-only** (legacy / awaiting Step 5): all three are
           ``None``; the adapter parses parameters from
           ``case.notes`` via ``_parse_notes_kv``. Today only the
           reflected-slab cases use this path; multi-region
           ``MeshTemplate`` will retire it.

        Validation rules:

        * ``materials`` and ``mesh_template`` MUST be set together
          (a partial inline path is a programming error).
        * ``registry_case`` and inline ``materials``/``mesh_template``
          MUST NOT both be set (one provisioning path per case).
        """
        has_registry = self.registry_case is not None
        has_inline = self.materials is not None and self.mesh_template is not None
        partial_inline = (
            (self.materials is None) ^ (self.mesh_template is None)
        )

        if partial_inline:
            raise ValueError(
                f"CrossMethodCase {self.case_id!r}: 'materials' and "
                f"'mesh_template' must be provided together "
                f"(materials={'set' if self.materials else 'None'}, "
                f"mesh_template="
                f"{'set' if self.mesh_template else 'None'})."
            )
        if has_registry and has_inline:
            raise ValueError(
                f"CrossMethodCase {self.case_id!r}: cannot have BOTH a "
                f"registry_case and inline materials/mesh_template — "
                f"pick one provisioning path."
            )

    def tolerance_for(self, adapter: SolverAdapter | str) -> float:
        """Return the absolute tolerance for ``adapter`` against
        this case's truth.

        Parameters
        ----------
        adapter : SolverAdapter | str
            Either an adapter instance (uses ``adapter.name``) or
            the adapter name directly.

        Raises
        ------
        KeyError
            If the adapter is not opted in for this case.
        """
        name = adapter if isinstance(adapter, str) else adapter.name
        if name not in self.tolerances:
            raise KeyError(
                f"CrossMethodCase {self.case_id!r}: no tolerance "
                f"declared for adapter {name!r}. "
                f"Available: {list(self.tolerances)}"
            )
        return self.tolerances[name]


def agreement_tolerance(
    case: CrossMethodCase,
    a: SolverAdapter | str,
    b: SolverAdapter | str,
) -> float:
    """The cross-method agreement tolerance is the **larger** of the
    two adapter-vs-truth tolerances.

    Both adapters reach the truth to within their respective
    floors; their pairwise agreement therefore cannot be tighter
    than the looser of the two. Tight pairwise tolerance is the
    classic reference-contamination anti-pattern (per
    :doc:`/skills/vv-principles` §6 "AI failure modes #6 reference
    contamination") — an unrealistic agreement floor falsely
    suggests the methods are bit-equal when they are merely both
    correct to ~1e-5.
    """
    return max(
        case.tolerance_for(a),
        case.tolerance_for(b),
    )
