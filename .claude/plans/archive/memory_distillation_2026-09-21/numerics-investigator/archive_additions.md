# Archive additions — numerics-investigator, 2026-09-21

Six new COLD files. Every one is a **lossless verbatim move** in the standard's `git mv` shape:
a `sed -n 'A,Bp'` extraction of a contiguous block of the pre-distillation `lessons.md`, with a
5-line provenance header prepended. No transcription, no paraphrase.

The materialised files are already written, ready to `cp` into
`.claude/agent-memory/numerics-investigator/_archive/`:

```
scratch/_memory_distillation/numerics-investigator/archive/
```

## Why these blocks went cold

Two reasons, and the second is the one worth noting. (1) They are the measured war stories the
standard sends cold — numbers, fixture names, per-probe tables. (2) **Five of the eight records
they point at live in `scratch/` and are UNTRACKED** (`git ls-files --error-unmatch`,
2026-09-21): `issue_344_null_space_structure.md`, `issue_344_kernel_basis.md`,
`adjoint_gram_ownership_audit.md`, `tau_under_ld_dip_analysis.md`, `rcond_rederivation.md`.
For those four campaigns the archive copy is the only version-controlled record of the
measurements, so moving the block cold *adds* durability rather than spending it.

## The six files, with the command that reproduces each

`$L` is the pre-distillation `.claude/agent-memory/numerics-investigator/lessons.md`
(1 228 lines, `wc -l`).

| target file under `_archive/` | source lines | digest successor |
|---|---|---|
| `issue_344_singularity_kernel_and_gauge.md` | 803–1012 | L23, L24 |
| `curvilinear_tau_ld_and_gram_ownership.md` | 1015–1165 | L25, L26 |
| `rcond_threshold_rederivation.md` | 1168–1228 | L27 |
| `krylov_composite_restart_and_stale_references.md` | 394–439, then 749–800 | L15, L22 |
| `rate_spectrum_certificate_and_angular_axis.md` | 582–635, then 638–748 | L19, L20, L21 |
| `curvilinear_seed_metric_and_ordering.md` | 354–393, then 440–579 | L14, L16, L17, L18 |

```sh
L=.claude/agent-memory/numerics-investigator/lessons.md
A=.claude/agent-memory/numerics-investigator/_archive

sed -n '803,1012p'  $L   # -> $A/issue_344_singularity_kernel_and_gauge.md
sed -n '1015,1165p' $L   # -> $A/curvilinear_tau_ld_and_gram_ownership.md
sed -n '1168,1228p' $L   # -> $A/rcond_threshold_rederivation.md
sed -n '394,439p;749,800p'  $L   # -> $A/krylov_composite_restart_and_stale_references.md
sed -n '582,635p;638,748p'  $L   # -> $A/rate_spectrum_certificate_and_angular_axis.md
sed -n '354,393p;440,579p'  $L   # -> $A/curvilinear_seed_metric_and_ordering.md
```

## Losslessness, verified

Each extraction was `diff`-ed against the corresponding slice of the materialised file (header
offset accounted for); all six returned empty. Re-runnable:

```sh
D=scratch/_memory_distillation/numerics-investigator/archive
diff <(sed -n '803,1012p' $L) <(sed -n '7,216p' $D/issue_344_singularity_kernel_and_gauge.md)
diff <(sed -n '1015,1165p' $L) <(sed -n '7,157p' $D/curvilinear_tau_ld_and_gram_ownership.md)
diff <(sed -n '1168,1228p' $L) <(sed -n '6,66p'  $D/rcond_threshold_rederivation.md)
diff <(sed -n '394,439p' $L) <(sed -n '6,51p'   $D/krylov_composite_restart_and_stale_references.md)
diff <(sed -n '749,800p' $L) <(sed -n '53,104p' $D/krylov_composite_restart_and_stale_references.md)
diff <(sed -n '582,635p' $L) <(sed -n '8,61p'   $D/rate_spectrum_certificate_and_angular_axis.md)
diff <(sed -n '638,748p' $L) <(sed -n '63,173p' $D/rate_spectrum_certificate_and_angular_axis.md)
diff <(sed -n '354,393p' $L) <(sed -n '8,47p'   $D/curvilinear_seed_metric_and_ordering.md)
diff <(sed -n '440,579p' $L) <(sed -n '49,188p' $D/curvilinear_seed_metric_and_ordering.md)
```

Total moved: **904 lines** written, of which **865** are the verbatim source blocks and 39 are
the six provenance headers.

## Blast radius before applying

`_archive/` files are Nexus graph nodes. These six are NEW files, so nothing can reference them
yet and no `dead_references` exposure is created by adding them. The exposure runs the other
way: the digest's `_archive/<name>.md` pointers must match the file names above exactly, and
the digest as proposed does.
