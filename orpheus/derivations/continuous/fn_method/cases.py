r"""Static metadata for the F_N method capability matrix.

Companion to :mod:`orpheus.derivations.continuous.peierls_nystrom.cases`
— exposes :func:`capability_rows`, the metadata-only function the
capability-matrix infrastructure discovers. The discovering agent is
the meta-generator :mod:`tools.verification.generate_capability_matrices`,
which walks every package under ``orpheus.derivations.continuous``,
calls whatever ``cases.capability_rows()`` it finds, and writes
``docs/theory/references/_<package>_capability_matrix.inc.rst``. (It
superseded the per-method generators — ``generate_fn_method_matrix``
and its siblings — which were deleted at ``045afeca``.)

The fn_method package ships three solver families anchored on the
Sood/Forster/Parsons LA-13511 (1999) benchmark catalogue and the
Neshat-Maiorino 1980 reflected-slab F_N paper:

1. **Multi-group :math:`k_\infty`** — :func:`...multi_group.k_inf.compute_kinf_*`
   for 1G / 2G no-upscatter / 2G general / mG (Eq 76). No spatial
   discretisation; pure rational-algebra in cross sections. Anchored
   on LA-13511 Tables 21+38 reference :math:`k_\infty` values.
2. **Bare-critical slab/sphere F_N** — :func:`...slab.one_group.solve_fn_slab_bare_critical`
   (Siewert-Benoist 1979 / Grandjean-Siewert 1979) and
   :func:`...sphere.one_group.solve_fn_sphere_bare_critical`
   (Siewert-Thomas 1986). 1G isotropic, vacuum BC. Anchored on
   LA-13511 Tables 1+7 / 5+7 critical dimensions.
3. **Reflected-slab F_N** — :func:`...slab.reflected.solve_fn_slab_reflected_critical`
   (Neshat-Maiorino 1980, *Ann. Nucl. Energy* 7, 79–81). 1G isotropic,
   finite reflector both sides. Anchored on NM Table 2 Burkart 1976
   "Exact" critical core half-thicknesses; not a registry case
   (parametrised over ``c1``/``c2``/``Δ``, no LA-13511 anchor).

The cylinder F_N stub (:mod:`...cylinder`) is registered in the
LA-13511 wide-slice STUBS but not yet activated — Mitsis-style
Wiener-Hopf is documented as non-convergent for bare cylinder
(Westfall-Metcalf 1972). Cylinder critical dimensions ship via
:mod:`orpheus.derivations.continuous.singular_eigenfunction.cylinder`
on a different mathematical pillar.

This function does **not** call any solver — it is safe to invoke at
Sphinx build time without paying any solver cost.
"""
from __future__ import annotations


def capability_rows() -> list[dict[str, object]]:
    """Return one row dict per shipped F_N method reference.

    Schema (codified in the meta-generator docstring):

    Required keys
        - ``name`` — registry / human-readable case name
        - ``geometry`` — ``"infinite" | "slab" | "sphere-1d" | "cylinder-1d"``
        - ``n_groups`` — energy-group count
        - ``n_regions`` — spatial-region count (``0`` for ``infinite``,
          ``1`` for bare, ``3`` for symmetric reflector + core)
        - ``bc`` — boundary-condition descriptor (RST-ready)
        - ``status`` — production-readiness label (``"shipped"``,
          ``"shipped (anisotropic-pending)"``, ``"stub"``, ...)

    Optional keys
        - ``accuracy`` — accuracy-class string for the matrix column
        - ``scattering_order`` — Legendre-order (``0`` isotropic, ``1`` linear-anisotropic)
        - ``multiplying`` — ``True`` for ``c > 1`` cases

    The matrix infrastructure auto-detects which optional columns are
    present and renders only those that appear in at least one row.
    """
    rows: list[dict[str, object]] = []

    # -----------------------------------------------------------------
    # Family 1: multi-group k_inf — pure rational-algebra in XS.
    # -----------------------------------------------------------------
    # 19 wide-slice + 2 first-slice = 21 catalogued LA-13511 cases.
    # Activated by ``compute_kinf_*`` in ``multi_group/k_inf.py``.
    # See ``tests/derivations/test_sood_registry_wide_kinf.py``.

    rows.append({
        "name": "k_inf — 1G isotropic (Sood Eq 19)",
        "geometry": "infinite",
        "n_groups": 1,
        "n_regions": 0,
        "bc": "—",
        "accuracy": "exact (rational algebra)",
        "scattering_order": 0,
        "multiplying": True,
        "status": (
            "shipped — 12 LA-13511 1G k_inf cases "
            "(PUa/PUb/Ua/Ub/Uc/Ud/UD2O/Ue/PU-1-1/UD2Oa-1-1/UD2Ob-1-1/UD2Oc-1-1)"
        ),
    })
    rows.append({
        "name": "k_inf — 2G no-upscatter (Sood Eq 29)",
        "geometry": "infinite",
        "n_groups": 2,
        "n_regions": 0,
        "bc": "—",
        "accuracy": "exact (rational algebra)",
        "scattering_order": 0,
        "multiplying": True,
        "status": (
            "shipped — 6 LA-13511 cases "
            "(PU-2-0-IN, U-2-0-IN, UAl-2-0-IN, URRa-2-0-IN, URRd-2-0-IN, UD2O-2-0-IN)"
        ),
    })
    rows.append({
        "name": "k_inf — 2G general (Sood Eq 28)",
        "geometry": "infinite",
        "n_groups": 2,
        "n_regions": 0,
        "bc": "—",
        "accuracy": "exact (2x2 dominant eigenvalue)",
        "scattering_order": 0,
        "multiplying": True,
        "status": (
            "shipped — 2 LA-13511 with-upscatter cases "
            "(URRb-2-0-IN, URRc-2-0-IN)"
        ),
    })
    # Both rows realise the SAME mG family (``compute_kinf_mg``, Sood
    # Eq 76) at different group counts. They carried the identical label
    # ``k_inf — mG general (Sood Eq 76)`` until 2026-08-09 (#345): only
    # the ``n_g`` column told them apart, which is off-pattern against
    # every sibling above (``1G isotropic``, ``2G no-upscatter``,
    # ``2G general``) and makes ``name`` unusable as a key. Naming the
    # order restores both.
    rows.append({
        "name": "k_inf — 3G general (Sood Eq 76)",
        "geometry": "infinite",
        "n_groups": 3,
        "n_regions": 0,
        "bc": "—",
        "accuracy": "exact (GxG dominant eigenvalue)",
        "scattering_order": 0,
        "multiplying": True,
        "status": "shipped — LA-13511 URR-3-0-IN (3G)",
    })
    rows.append({
        "name": "k_inf — 6G general (Sood Eq 76)",
        "geometry": "infinite",
        "n_groups": 6,
        "n_regions": 0,
        "bc": "—",
        "accuracy": "exact (GxG dominant eigenvalue)",
        "scattering_order": 0,
        "multiplying": True,
        "status": "shipped — LA-13511 URR-6-0-IN (6G)",
    })

    # -----------------------------------------------------------------
    # Family 2: bare-critical slab F_N (Siewert-Benoist 1979 +
    # Grandjean-Siewert 1979) — 1G isotropic, vacuum BC.
    # -----------------------------------------------------------------

    rows.append({
        "name": "F_N slab bare-critical (Siewert-Benoist + Grandjean-Siewert 1979)",
        "geometry": "slab",
        "n_groups": 1,
        "n_regions": 1,
        "bc": "vacuum (symmetric)",
        "accuracy": "≤ 5e-6 abs at N=10–12 vs Sood truth",
        "scattering_order": 0,
        "multiplying": True,
        "status": (
            "shipped — 4 LA-13511 1G bare slab cases "
            "(Ua-1-0-SL, PUa-1-0-SL, PUb-1-0-SL, UD2O-1-0-SL)"
        ),
    })

    # -----------------------------------------------------------------
    # Family 3: bare-critical sphere F_N (Siewert-Thomas 1986).
    # -----------------------------------------------------------------

    rows.append({
        "name": "F_N sphere bare-critical (Siewert-Thomas 1986)",
        "geometry": "sphere-1d",
        "n_groups": 1,
        "n_regions": 1,
        "bc": "vacuum",
        "accuracy": "≤ 5e-8 abs at N=10 vs Sood truth",
        "scattering_order": 0,
        "multiplying": True,
        "status": (
            "shipped — 3 LA-13511 1G bare sphere cases "
            "(Ua-1-0-SP, PUb-1-0-SP, UD2O-1-0-SP)"
        ),
    })

    # -----------------------------------------------------------------
    # Family 4: reflected-slab F_N (Neshat-Maiorino 1980).
    # -----------------------------------------------------------------

    rows.append({
        "name": "F_N reflected-slab critical (Neshat-Maiorino 1980)",
        "geometry": "slab",
        "n_groups": 1,
        "n_regions": 3,
        "bc": "vacuum on outer reflector face; symmetric",
        "accuracy": "≤ 5e-5 abs at F_7 vs NM Table 2 Burkart 1976 \"Exact\"",
        "scattering_order": 0,
        "multiplying": True,
        "status": (
            "shipped — 8 NM Table 2 cases "
            "(c_1∈{1.02,1.10,1.30,1.50}, c_2∈{0.10..0.90}, Δ∈{0.5,1,2,5})"
        ),
    })

    # -----------------------------------------------------------------
    # Family 5: bare-critical cylinder F_N — STUB (Mitsis-style WH does
    # not converge for cylinder; Westfall-Metcalf 1972 documents this).
    # The cylinder critical dimension ships via singular_eigenfunction
    # on a different mathematical pillar (Mitsis-WM Fredholm method).
    # -----------------------------------------------------------------

    rows.append({
        "name": "F_N cylinder bare-critical",
        "geometry": "cylinder-1d",
        "n_groups": 1,
        "n_regions": 1,
        "bc": "vacuum",
        "accuracy": "—",
        "scattering_order": 0,
        "multiplying": True,
        "status": (
            "stub — Mitsis Wiener-Hopf not convergent for bare "
            "cylinder (WM-72); see ``singular_eigenfunction.cylinder`` "
            "for the alternative pillar"
        ),
    })

    return rows
