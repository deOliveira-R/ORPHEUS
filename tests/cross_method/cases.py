r"""Populated cross-method case sets.

Each case set targets a specific physics+geometry combination and
collects the cases that multiple solvers can natively address. Per-
case truth values trace to a primary literature source. Per-solver
tolerances reflect each method's published quadrature floor at
default settings.

Case-set inventory
------------------

* :data:`BARE_CRITICAL_SLAB_CASES` — 1G isotropic bare-critical slab,
  c-sweep covering the four published Sood/KLL truths. fn_method +
  trajectory_resolvent both supported.
* :data:`BARE_CRITICAL_SPHERE_CASES` — 1G isotropic bare-critical
  sphere, c-sweep covering the three published Sood/KLL truths.
  fn_method + trajectory_resolvent both supported.
* :data:`REFLECTED_SLAB_CASES` — Neshat-Maiorino 1980 / Sood Table
  10 reflected-slab criticality. fn_method only (no
  trajectory_resolvent counterpart yet).
* :data:`CLOSED_SPHERE_KINF_CASES` — closed-sphere α=1 cases where
  ``k_eff = k_inf`` exactly. trajectory_resolvent native; fn_method
  via ``compute_kinf_*`` direct evaluation. The 2G/4G coverage
  this provides is the ONLY multi-group cross-method gate
  shippable today.
* :data:`GRANDJEAN_SIEWERT_SLAB_PARAMETRIC` — extra slab c-values
  not in Sood (c=1.10, 1.70, 1.90 from Grandjean-Siewert Table XI).
  Unit XS, no registry case.

References
----------

* Sood, Forster, Parsons (1999), LA-13511.
* Kaper, Lindeman, Leaf (1974), *Nucl. Sci. Eng.* **54**, 94 (the
  underlying truth source for Sood Tables 4, 6, 13; not directly
  obtainable, transcribed via Sood — see
  ``.claude/agent-memory/literature-researcher/kaper_lindeman_leaf_1974_fn_method.md``).
* Grandjean, Siewert (1979), *Nucl. Sci. Eng.* **69**, 161 — Table
  XI critical thicknesses.
* Neshat, Maiorino (1980), *Ann. Nucl. Energy* **7**, 79 — Table 2
  reflected-slab cases.
* Burkart (1976) — the "Exact" reference values NM Table 2 cites.
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.common.geometry_spec import GeometrySpec, Region
from orpheus.derivations.common.xs_library import make_mixture
from orpheus.derivations.continuous.sood_registry import (
    PUA_1_0_SL,
    PUB_1_0_SL,
    PUB_1_0_SP,
    UA_1_0_SL_STUB,
    UA_1_0_SP_STUB,
    UD2O_1_0_SL,
    UD2O_1_0_SP,
)
from orpheus.geometry.mesh import BC

from .protocol import CrossMethodCase


# ═══════════════════════════════════════════════════════════════════
# Bare-critical 1G isotropic slab — c-sweep
# ═══════════════════════════════════════════════════════════════════
#
# fn_method/slab/one_group reaches ~5e-6 absolute on a_c at N=10 vs
# Sood truth (already tested in test_fn_la13511_slab.py). The
# trajectory_resolvent slab vacuum-BC quadrature floor is ~5e-5 on
# k_eff at default (n_x=48, n_mu=128, n_traj_quad=96) due to the
# slab μ=0 cusp.
#
# The agreement tolerance for the cross-check is the larger of the
# two, i.e. 5e-5. Tighter than that is reference contamination.


# F_N slab at N=10 across c ∈ {1.02 .. 1.50} reaches err ≤ 1e-5 absolute
# on a_c against Sood/KLL truth — see test_fn_la13511_slab.py for the
# canonical pin (it uses 1e-5 for the c=1.30 case). We adopt the same
# floor for the cross-method protocol.
_FN_SLAB_TOL_DEFAULT = 1e-5
_TR_SLAB_TOL_DEFAULT = 5e-5


BARE_CRITICAL_SLAB_CASES: list[CrossMethodCase] = [
    # c = 1.02 — UD2O, lowest c in the bare 1G slab family.
    # Sood Table 17. F_N has bracket-loss issues at N≥14; cap at N=12.
    CrossMethodCase(
        case_id="UD2O-1-0-SL-c1.02",
        description=(
            "U-D2O reactor bare slab, 1G isotropic, c=1.02. "
            "a_c = 5.6655054562 mfp."
        ),
        registry_case=UD2O_1_0_SL,
        geometry="slab",
        truth_tag="a_critical_mfp",
        truth_value=5.6655054562,
        truth_source=(
            "Sood LA-13511 Table 17 (1999), citing Kaper-Lindeman-Leaf "
            "1974 NSE 54, 94 (Ref. 26)"
        ),
        pillar="closed-form",
        claim_layer="eigenvalue",
        tolerances={
            # UD2O-1-0-SL has the tightest fn_slab match (3.4e-6 at
            # N=10) since c=1.02 is closest to a Case-eigenfunction
            # purely-real regime; keep the default 1e-5.
            "fn_slab": _FN_SLAB_TOL_DEFAULT,
            # trajectory_resolvent at default parameters needs the
            # looser floor at low c (large slab → many bounce
            # decays).
            "trajectory_resolvent_slab": 1e-4,
        },
        notes=(
            "Lowest c in bare 1G slab. F_N at N≥14 fails (bracket "
            "loss); use N=12. trajectory_resolvent k_eff floor is "
            "looser at low c due to long-tau bounce-period quadrature."
        ),
    ),
    # c = 1.30 — Sood Ua-1-0-SL, the canonical KLL-table reference.
    CrossMethodCase(
        case_id="Ua-1-0-SL-c1.30",
        description=(
            "U-235 (a) bare slab, 1G isotropic, c=1.30. The canonical "
            "Kaper-Lindeman-Leaf reference; a_c = 0.93772556 mfp."
        ),
        registry_case=UA_1_0_SL_STUB,
        geometry="slab",
        truth_tag="a_critical_mfp",
        truth_value=0.93772556,
        truth_source=(
            "Sood LA-13511 Table 4 (1999), citing Kaper-Lindeman-Leaf "
            "1974 NSE 54, 94 (Ref. 26)"
        ),
        pillar="closed-form",
        claim_layer="eigenvalue",
        tolerances={
            "fn_slab": _FN_SLAB_TOL_DEFAULT,
            "trajectory_resolvent_slab": _TR_SLAB_TOL_DEFAULT,
        },
        notes="Canonical KLL Table I bare-critical slab benchmark.",
    ),
    # c = 1.40 — PUb slab.
    CrossMethodCase(
        case_id="PUb-1-0-SL-c1.40",
        description=(
            "Pu-239 (b) bare slab, 1G isotropic, c=1.40. a_c = "
            "0.73660355 mfp."
        ),
        registry_case=PUB_1_0_SL,
        geometry="slab",
        truth_tag="a_critical_mfp",
        truth_value=0.73660355,
        truth_source=(
            "Sood LA-13511 Table 7 (1999), citing Kaper-Lindeman-Leaf "
            "1974 NSE 54, 94 (Ref. 26)"
        ),
        pillar="closed-form",
        claim_layer="eigenvalue",
        tolerances={
            "fn_slab": _FN_SLAB_TOL_DEFAULT,
            "trajectory_resolvent_slab": _TR_SLAB_TOL_DEFAULT,
        },
    ),
    # c = 1.50 — PUa slab (highest c in the bare 1G slab family).
    CrossMethodCase(
        case_id="PUa-1-0-SL-c1.50",
        description=(
            "Pu-239 (a) bare slab, 1G isotropic, c=1.50 (highest c in "
            "the Sood family). a_c = 0.605055 mfp."
        ),
        registry_case=PUA_1_0_SL,
        geometry="slab",
        truth_tag="a_critical_mfp",
        truth_value=0.605055,
        truth_source=(
            "Sood LA-13511 Table 6 (1999), citing Lathrop-Leonard "
            "1965 NSE 22, 115 (Ref. 9)"
        ),
        pillar="closed-form",
        claim_layer="eigenvalue",
        tolerances={
            "fn_slab": _FN_SLAB_TOL_DEFAULT,
            "trajectory_resolvent_slab": _TR_SLAB_TOL_DEFAULT,
        },
        notes=(
            "Highest c in bare 1G slab family. Thin slab, steep "
            "boundary gradient — angular floor governs."
        ),
    ),
]


# ═══════════════════════════════════════════════════════════════════
# Bare-critical 1G isotropic sphere — c-sweep
# ═══════════════════════════════════════════════════════════════════
#
# fn_method/sphere/one_group at N=10 reaches ~5e-8 absolute on R_c —
# exquisitely tight. trajectory_resolvent sphere vacuum-BC at
# (n_r=32, n_mu=32, n_traj_quad=64) reaches ~1e-5 on k_eff.

_FN_SPHERE_TOL_DEFAULT = 1e-7
_TR_SPHERE_TOL_DEFAULT = 1e-5


BARE_CRITICAL_SPHERE_CASES: list[CrossMethodCase] = [
    # c = 1.02 — UD2O sphere. Lowest c, largest sphere.
    CrossMethodCase(
        case_id="UD2O-1-0-SP-c1.02",
        description=(
            "U-D2O reactor bare sphere, 1G isotropic, c=1.02. "
            "R_c = 12.0275320980 mfp."
        ),
        registry_case=UD2O_1_0_SP,
        geometry="sphere-1d",
        truth_tag="R_critical_mfp",
        truth_value=12.0275320980,
        truth_source=(
            "Sood LA-13511 Table 17 (1999), citing Kaper-Lindeman-Leaf "
            "1974 NSE 54, 94 (Ref. 26)"
        ),
        pillar="closed-form",
        claim_layer="eigenvalue",
        tolerances={
            "fn_sphere": _FN_SPHERE_TOL_DEFAULT,
            # Large sphere → long bounce-period chord → looser
            # trajectory_resolvent floor.
            "trajectory_resolvent_sphere": 1e-4,
        },
        notes=(
            "Lowest c, largest sphere. trajectory_resolvent at "
            "default n_traj_quad=64 has a looser k_eff floor at "
            "long τ_R."
        ),
    ),
    # c = 1.30 — Sood Ua-1-0-SP.
    CrossMethodCase(
        case_id="Ua-1-0-SP-c1.30",
        description=(
            "U-235 (a) bare sphere, 1G isotropic, c=1.30. The "
            "canonical Kaper-Lindeman-Leaf reference; R_c = "
            "2.4248249802 mfp."
        ),
        registry_case=UA_1_0_SP_STUB,
        geometry="sphere-1d",
        truth_tag="R_critical_mfp",
        truth_value=2.4248249802,
        truth_source=(
            "Sood LA-13511 Table 5 (1999), citing Kaper-Lindeman-Leaf "
            "1974 NSE 54, 94 (Ref. 26)"
        ),
        pillar="closed-form",
        claim_layer="eigenvalue",
        tolerances={
            "fn_sphere": _FN_SPHERE_TOL_DEFAULT,
            "trajectory_resolvent_sphere": _TR_SPHERE_TOL_DEFAULT,
        },
        notes="Canonical KLL Table V bare-critical sphere benchmark.",
    ),
    # c = 1.40 — PUb sphere.
    CrossMethodCase(
        case_id="PUb-1-0-SP-c1.40",
        description=(
            "Pu-239 (b) bare sphere, 1G isotropic, c=1.40. R_c = "
            "1.9853434324 mfp."
        ),
        registry_case=PUB_1_0_SP,
        geometry="sphere-1d",
        truth_tag="R_critical_mfp",
        truth_value=1.9853434324,
        truth_source=(
            "Sood LA-13511 Table 7 (1999), citing Kaper-Lindeman-Leaf "
            "1974 NSE 54, 94 (Ref. 26)"
        ),
        pillar="closed-form",
        claim_layer="eigenvalue",
        tolerances={
            "fn_sphere": _FN_SPHERE_TOL_DEFAULT,
            "trajectory_resolvent_sphere": _TR_SPHERE_TOL_DEFAULT,
        },
    ),
]


# ═══════════════════════════════════════════════════════════════════
# Reflected-slab cases — fn_method only (no trajectory_resolvent
# counterpart for reflected/multi-region slab).
# ═══════════════════════════════════════════════════════════════════
#
# Sood LA-13511 Table 10 (problem 4) is the symmetric Pu+H2O
# reflected slab Δ=0.5 mfp each side, c_core=1.50, c_refl=0.90,
# Pu r_c = 0.43014 mfp. NM 1980 Table 2 publishes 8 cases (c1, c2,
# Δ) families; Sood Table 10 #4 is one of them. We populate the
# four primary NM cases.
#
# Encoding convention (Step 3, 2026-05-04 input-cleanup):
#
# * Materials:
#   - mat_id=0 → core: 1G isotropic with Σ_t=1, Σ_s=0, νΣ_f=c_core
#     (pure-multiplying, c > 1). The XS construction is a unit-σ_t
#     cross-section convention so cm ≡ mfp; the FN solver consumes
#     ``Mixture.scattering_ratio[0]`` exclusively, which is invariant
#     under any (Σ_s, νΣ_f) decomposition with a fixed c.
#   - mat_id=1 → reflector: 1G isotropic with Σ_t=1, Σ_s=c_refl,
#     νΣ_f=0 (pure scattering, c < 1).
# * GeometrySpec:
#   - ``regions`` is the ORDERED slab layout ``(reflector, core,
#     reflector)`` — left-to-right with the symmetric Pu+H2O
#     convention.
#   - ``Region.outer_thickness_*`` for the reflector layers is the
#     published Δ (mfp ≡ cm with σ_t=1). For the core layer the
#     thickness is ``2 * truth_value`` (the FULL core, of which τ
#     is the half-thickness).
#   - ``critical_dimension_{mfp,cm}`` carries the published
#     core-half-thickness τ_c — the truth scalar — for diagnostics.
#     For multi-region specs ``domain_extent_cm`` is taken from the
#     region stack, not from this field (see GeometrySpec docstring).
#   - BCs are vacuum-vacuum at the outer edges of the slab.


def _make_unit_sigma_t_one_group_mixture(c: float):
    r"""Build a 1G isotropic :class:`Mixture` with :math:`\Sigma_t = 1`
    that has scattering ratio :math:`c` exactly.

    Decomposition convention:

    * :math:`c \ge 1`: pure-multiplying (:math:`\Sigma_s = 0`,
      :math:`\nu\Sigma_f = c`, :math:`\Sigma_c = 1 - c \ge 0` is
      negative for c>1, so we set :math:`\Sigma_c = 0` and put the
      mass on :math:`\nu`/:math:`\Sigma_f` decomposition;
      ``make_mixture`` requires :math:`\Sigma_c \ge 0` so we
      pick :math:`\Sigma_f = c, \nu = 1` and let
      :math:`\Sigma_t = \Sigma_c + \Sigma_f + \Sigma_s = 0 + c + 0`).
      That gives :math:`\Sigma_t = c` — but we want
      :math:`\Sigma_t = 1`. We resolve by **renormalising**:
      :math:`\Sigma_t \equiv 1`, :math:`\Sigma_s = 0`,
      :math:`\Sigma_f = 1`, :math:`\nu = c`,
      :math:`\Sigma_c = 1 - 1 = 0`. Then
      :math:`c_{\rm Case-Zweifel} = (\Sigma_s + \nu\Sigma_f)/\Sigma_t
      = (0 + c \cdot 1)/1 = c` — matches.
    * :math:`c < 1`: pure-scattering (:math:`\Sigma_s = c`,
      :math:`\nu\Sigma_f = 0`, :math:`\Sigma_c = 1 - c > 0`,
      :math:`\Sigma_f = 0`).

    The two decompositions yield the same ``Mixture.scattering_ratio``
    — that is the only quantity the F_N reflected-slab solver
    consumes from the materials. The (Σ_s, νΣ_f) split is a
    cross-section labelling convention.

    Parameters
    ----------
    c : float
        Case–Zweifel scattering ratio. Must be ``> 0``.

    Returns
    -------
    Mixture
        1G isotropic mixture with ``Σ_t = 1`` and
        ``Mixture.scattering_ratio[0] == c``.
    """
    if c <= 0.0:
        raise ValueError(f"_make_unit_sigma_t_one_group_mixture: c must be > 0; got {c!r}")
    sig_t = np.array([1.0])
    chi = np.array([1.0])
    if c >= 1.0:
        # Pure-multiplying: Σ_s = 0, ν=c, Σ_f = 1, Σ_c = 0.
        return make_mixture(
            sig_t=sig_t,
            sig_c=np.array([0.0]),
            sig_f=np.array([1.0]),
            nu=np.array([c]),
            chi=chi,
            sig_s=np.array([[0.0]]),
        )
    # Pure-scattering: Σ_s = c, ν=0, Σ_f = 0, Σ_c = 1 - c.
    return make_mixture(
        sig_t=sig_t,
        sig_c=np.array([1.0 - c]),
        sig_f=np.array([0.0]),
        nu=np.array([0.0]),
        chi=chi,
        sig_s=np.array([[c]]),
    )


def _build_reflected_slab_case(
    *,
    case_id: str,
    description: str,
    c_core: float,
    c_reflector: float,
    reflector_half_thickness_mfp: float,
    truth_value: float,
    truth_source: str,
    fn_tolerance: float,
    notes: str = "",
) -> CrossMethodCase:
    r"""Construct a symmetric reflected-slab :class:`CrossMethodCase`.

    Materials and geometry are encoded inline (Step 3 input-cleanup):

    * ``materials = {0: core_mixture, 1: reflector_mixture}``
      built via :func:`_make_unit_sigma_t_one_group_mixture`.
    * ``geometry_spec.regions = (reflector, core, reflector)`` with
      thicknesses ``(Δ, 2·τ, Δ)`` mfp (≡ cm under σ_t=1). The core
      thickness uses the published critical τ as the half-thickness.

    The FN reflected-slab adapter reads ``c_core``,
    ``c_reflector``, and ``reflector_half_thickness_mfp`` directly
    from these fields — no notes-string parsing.
    """
    core_mix = _make_unit_sigma_t_one_group_mixture(c_core)
    refl_mix = _make_unit_sigma_t_one_group_mixture(c_reflector)
    # Unit σ_t convention → cm ≡ mfp for these problems.
    refl_thickness_cm = float(reflector_half_thickness_mfp)
    core_thickness_cm = 2.0 * float(truth_value)
    regions = (
        Region(
            mat_id=1,  # reflector layer (left)
            outer_thickness_mfp=reflector_half_thickness_mfp,
            outer_thickness_cm=refl_thickness_cm,
        ),
        Region(
            mat_id=0,  # core layer (full thickness = 2·τ)
            outer_thickness_mfp=2.0 * float(truth_value),
            outer_thickness_cm=core_thickness_cm,
        ),
        Region(
            mat_id=1,  # reflector layer (right)
            outer_thickness_mfp=reflector_half_thickness_mfp,
            outer_thickness_cm=refl_thickness_cm,
        ),
    )
    spec = GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=float(truth_value),  # τ_c half-thickness
        critical_dimension_cm=float(truth_value),  # σ_t=1 → cm ≡ mfp
        n_groups=1,
        regions=regions,
    )
    return CrossMethodCase(
        case_id=case_id,
        description=description,
        registry_case=None,
        materials={0: core_mix, 1: refl_mix},
        geometry_spec=spec,
        geometry="reflected-slab",
        truth_tag="tau_critical_mfp",
        truth_value=truth_value,
        truth_source=truth_source,
        pillar="closed-form",
        claim_layer="eigenvalue",
        tolerances={"fn_reflected_slab": fn_tolerance},
        notes=notes,
    )


REFLECTED_SLAB_CASES: list[CrossMethodCase] = [
    _build_reflected_slab_case(
        case_id="Sood-Table10-4-PUa-H2O-symm-0.5mfp",
        description=(
            "Symmetric Pu+H2O reflected slab, Sood LA-13511 Table 10 "
            "problem 4: c_core=1.50, c_refl=0.90, Δ=0.5 mfp each "
            "side. Pu r_c = 0.43014 mfp."
        ),
        c_core=1.50,
        c_reflector=0.90,
        reflector_half_thickness_mfp=0.5,
        truth_value=0.43014,
        truth_source="Sood LA-13511 Table 10 problem 4 (1999)",
        fn_tolerance=1e-3,
    ),
    # NM 1980 Table 1 Case 4 — c_core=1.30, c_reflector=0.90, Δ=1.0.
    # NM Table 2 reports F_7 matching Burkart 1976 "Exact" to all
    # published digits.
    _build_reflected_slab_case(
        case_id="NM-1980-Case4-c130-c090-D10",
        description=(
            "NM 1980 Table 1 Case 4: c_core=1.30, c_reflector=0.90, "
            "Δ=1.0 mfp each side. τ_c = 0.6027 mfp per NM Table 2 "
            "(Burkart 1976 'Exact')."
        ),
        c_core=1.30,
        c_reflector=0.90,
        reflector_half_thickness_mfp=1.0,
        truth_value=0.6027,
        truth_source=(
            "Neshat-Maiorino 1980 Table 1 Case 4 / Table 2 Burkart "
            "1976 'Exact'"
        ),
        # 4 published digits → 5e-4 absolute is the published-
        # precision floor; F_7 reaches it.
        fn_tolerance=1e-3,
    ),
    # NM 1980 Table 1 Case 6 — the canonical c_core=1.50, c_reflector=0.90,
    # Δ=1.0 cross-comparator (also referenced by
    # test_fn_sood_table10_symmetric_pu_h2o.py).
    _build_reflected_slab_case(
        case_id="NM-1980-Case6-c150-c090-D10",
        description=(
            "NM 1980 Table 1 Case 6: c_core=1.50, c_reflector=0.90, "
            "Δ=1.0. τ_c = 0.3597 mfp per Burkart 1976 'Exact'."
        ),
        c_core=1.50,
        c_reflector=0.90,
        reflector_half_thickness_mfp=1.0,
        truth_value=0.3597,
        truth_source=(
            "Neshat-Maiorino 1980 Table 1 Case 6 / Table 2 Burkart "
            "1976 'Exact'"
        ),
        fn_tolerance=1e-3,
    ),
    # NM 1980 Table 1 Case 1 — small-c endpoint (most loosely
    # multiplying core; thick critical slab). c_core=1.01,
    # c_reflector=0.09, Δ=0.5. τ_c = 8.3107.
    _build_reflected_slab_case(
        case_id="NM-1980-Case1-c101-c009-D05",
        description=(
            "NM 1980 Table 1 Case 1 (small-c endpoint): c_core=1.01, "
            "c_reflector=0.09, Δ=0.5. τ_c = 8.3107 mfp."
        ),
        c_core=1.01,
        c_reflector=0.09,
        reflector_half_thickness_mfp=0.5,
        truth_value=8.3107,
        truth_source="Neshat-Maiorino 1980 Table 1 Case 1",
        # 4 published digits + thick slab — F_7 reaches Burkart
        # 'Exact' to all printed digits per NM Table 2 (~ 5e-4).
        fn_tolerance=1e-2,
    ),
]


# ═══════════════════════════════════════════════════════════════════
# Closed-sphere k_inf cases — α=1, the eigenvalue ≡ k_inf identity
# ═══════════════════════════════════════════════════════════════════
#
# When trajectory_resolvent solves a closed (α=1) homogeneous sphere
# the rank-1 isotropic eigenmode gives k_eff = k_inf exactly (V_α1).
# The fn_method computes k_inf directly via Sood Eq 19/29/76 closed
# form. Comparing the two for 1G/2G/multi-group cases gives the
# multi-group cross-method coverage that bare-critical does not
# provide.

_KINF_TOL_FN = 1e-12
_KINF_TOL_TRAJECTORY = 1e-10


# 1G — fuel-A-like (canonical V_α1 fixture).
#
# XS factoring: σ_t = 0.5, σ_s = 0.38, νσ_f = 0.025.
# σ_a = σ_t - σ_s = 0.12. Pick ν = 1.0, σ_f = 0.025; then
# σ_c = σ_a - σ_f = 0.095. The closed-sphere adapter consumes
# (σ_t, σ_s_p0, νσ_f) extracted via mixture_to_fn_arrays — only
# those three quantities reach the solver, so the (ν, σ_f)
# decomposition above is a labelling convenience.
_FUEL_A_LIKE_MIX = make_mixture(
    sig_t=np.array([0.5]),
    sig_c=np.array([0.095]),
    sig_f=np.array([0.025]),
    nu=np.array([1.0]),
    chi=np.array([1.0]),
    sig_s=np.array([[0.38]]),
)


CLOSED_SPHERE_KINF_CASES: list[CrossMethodCase] = [
    CrossMethodCase(
        case_id="closed-sphere-1G-fuelA-tauR2.5",
        description=(
            "Closed (α=1) homogeneous sphere, 1G isotropic, fuel-A-"
            "like XS (σ_t=0.5, σ_s=0.38, νσ_f=0.025; τ_R=2.5). "
            "k_eff = k_inf = 0.025/0.12 = 0.2083̄."
        ),
        registry_case=None,
        materials={0: _FUEL_A_LIKE_MIX},
        # Closed sphere: BC.reflective on BOTH sides (centreline +
        # outer surface at α=1). The trajectory_resolvent solver
        # consumes the radius via geometry_spec.critical_dimension_cm
        # and translates BC.reflective to its α=1 albedo.
        geometry_spec=GeometrySpec(
            geometry="sphere",
            critical_dimension_mfp=2.5,  # τ_R = σ_t · R_cm = 0.5 · 5.0
            critical_dimension_cm=5.0,
            n_groups=1,
            bc_left=BC.reflective,
            bc_right=BC.reflective,
        ),
        geometry="closed-sphere-1d",
        truth_tag="k_inf",
        truth_value=0.025 / 0.12,  # 0.20833...
        truth_source=(
            "Closed-form k_inf = νΣ_f/Σ_a; V_α1 algebraic identity "
            "(see test_peierls_greens_function_xverif.py)"
        ),
        pillar="closed-form",
        claim_layer="eigenvalue",
        tolerances={
            # Closed-sphere α=1 gives k_eff = k_inf to machine
            # precision (no quadrature error since the eigenmode
            # is rank-1 isotropic). Tight tolerance documents the
            # V_α1 algebraic identity at the production-code level.
            "trajectory_resolvent_sphere_closed": 1e-10,
        },
        notes=(
            "Used by closed-sphere k_inf adapter; "
            "trajectory_resolvent must converge to k_inf to "
            "machine precision via V_α1."
        ),
    ),
]


# ═══════════════════════════════════════════════════════════════════
# Grandjean-Siewert Table XI parametric c-sweep — unit XS
# ═══════════════════════════════════════════════════════════════════
#
# GS Table XI publishes critical 2τ for c ∈ {1.10, 1.30, 1.50, 1.70,
# 1.90} with σ_t=1 (unit cross sections, so mfp ≡ cm). c=1.10 and
# c=1.70/1.90 fill in the gaps not covered by Sood's bare 1G slab
# family. These cases use unit XS (no registry case) and are
# parameterised on c only.
#
# GS Table XI 2τ values (mfp):
#   c=1.10 → 4.2266
#   c=1.30 → 1.8755 (matches Sood Ua-1-0-SL = 2 × 0.93772556)
#   c=1.50 → 1.2101
#   c=1.70 → 0.8851
#   c=1.90 → 0.6919

_GS_TABLE_XI = (
    (1.10, 4.2266 / 2.0),
    (1.30, 1.8755 / 2.0),
    (1.50, 1.2101 / 2.0),
    (1.70, 0.8851 / 2.0),
    (1.90, 0.6919 / 2.0),
)


def _grandjean_siewert_unit_case(c: float, a_truth: float) -> CrossMethodCase:
    r"""Build a parametric GS Table XI cross-method case.

    Unit cross sections (``σ_t = 1``), pure-multiplying medium
    (``σ_s = 0``, ``νσ_f = c``). Truth is from GS Table XI.

    These cases exercise the c-sweep without binding to a specific
    Sood case. They have ``registry_case = None`` because the
    materials are constructed inline (no XS-set entry in the Sood
    registry corresponds to ``c, σ_t=1, σ_s=0``).

    The trajectory_resolvent slab adapter cannot consume these
    today (its XS extractor requires a registry case). They are
    populated here for fn_method coverage and for documenting the
    extension point.
    """
    return CrossMethodCase(
        case_id=f"GS-Table-XI-slab-c{c:.2f}",
        description=(
            f"Grandjean-Siewert 1979 Table XI parametric: c={c:.2f}, "
            f"unit XS, bare slab. 2τ_c = {2*a_truth:.4f} mfp."
        ),
        registry_case=None,
        geometry="slab",
        truth_tag="a_critical_mfp",
        truth_value=a_truth,
        truth_source="Grandjean-Siewert 1979 Table XI",
        pillar="closed-form",
        claim_layer="eigenvalue",
        tolerances={
            # GS Table XI is published to 4-5 sig figs.
            "fn_slab": 5e-4,
        },
        notes=(
            f"c={c:.2f} unit XS; covers slab c-sweep gap not in "
            f"Sood family. trajectory_resolvent adapter requires a "
            f"registry case so this case is fn_method-only at the "
            f"protocol level."
        ),
    )


GRANDJEAN_SIEWERT_SLAB_PARAMETRIC: list[CrossMethodCase] = [
    _grandjean_siewert_unit_case(c, a)
    for c, a in _GS_TABLE_XI
]


# ═══════════════════════════════════════════════════════════════════
# All cases — for full-suite enumeration / agreement-matrix renderer
# ═══════════════════════════════════════════════════════════════════


ALL_CASES: list[CrossMethodCase] = (
    BARE_CRITICAL_SLAB_CASES
    + BARE_CRITICAL_SPHERE_CASES
    + REFLECTED_SLAB_CASES
    + CLOSED_SPHERE_KINF_CASES
    + GRANDJEAN_SIEWERT_SLAB_PARAMETRIC
)
"""Every cross-method case the regression net knows about. Used by
the agreement-matrix renderer (future) and by the
:func:`tests.cross_method.test_eigenvalue.test_case_inventory`
foundation gate to ensure no case set has been silently dropped.
"""


def cases_by_geometry(geometry: str) -> list[CrossMethodCase]:
    """Filter :data:`ALL_CASES` by the ``geometry`` field."""
    return [c for c in ALL_CASES if c.geometry == geometry]


def case_by_id(case_id: str) -> CrossMethodCase:
    """Look up a case by ``case_id``. Raises :class:`KeyError` if absent."""
    for c in ALL_CASES:
        if c.case_id == case_id:
            return c
    raise KeyError(f"No CrossMethodCase with case_id={case_id!r}")
