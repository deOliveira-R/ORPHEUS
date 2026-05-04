#!/bin/bash
# Replace <!-- COMMIT-HASH --> in 16 issue comments + #140 body with the
# Peierls docs cleanup commit hash.
set -euo pipefail

TIP=18a852b
RANGE=742d3b0..18a852b
REPLACEMENT='[`'"$TIP"'`](https://github.com/deOliveira-R/ORPHEUS/commit/'"$TIP"') (cleanup tip; full range `'"$RANGE"'`, 25 commits)'
SCRATCH=/workspaces/ORPHEUS/.claude/scratch/peierls-cleanup

# Map: <comment-id> <local-draft-file>
declare -a MAP=(
  "4348744341 094_bickley_retirement.md"
  "4348744491 100_sphere_white_bc_retraction.md"
  "4348744624 112_marshak_empirical_scan.md"
  "4348744788 113_failed_scheme_forensic.md"
  "4348745026 117_moment_form_archived.md"
  "4348745180 119_F5_villarino_stammler.md"
  "4348745355 122_lambert_marshak_gauge.md"
  "4348745506 123_L19_stability_protocol.md"
  "4348745701 129_planar_limit.md"
  "4348745846 131_probe_cascade.md"
  "4348746014 132_class_b_mr_mg_falsification.md"
  "4348746147 133_phase5_continuous_mu_retreat.md"
  "4348746336 133_section_22_9_audit.md"
  "4348746496 134_cp_chord_migration.md"
  "4348746667 135_volume_kernel_migration.md"
  "4348746836 136_residual_leggauss.md"
)

for entry in "${MAP[@]}"; do
  cid="${entry%% *}"
  file="${entry#* }"
  in="$SCRATCH/comments/$file"
  out="$SCRATCH/comments/_filled_${file}"
  python -c "
import sys
src = open(sys.argv[1]).read()
src = src.replace('<!-- COMMIT-HASH -->', sys.argv[2])
open(sys.argv[3], 'w').write(src)
" "$in" "$REPLACEMENT" "$out"
  gh api -X PATCH "repos/deOliveira-R/ORPHEUS/issues/comments/$cid" \
    -f body="$(cat "$out")" \
    --jq '.html_url' || echo "  (PATCH failed for comment $cid)"
done

# Issue #140 body — strip the leading 4 metadata HTML comments before posting
in="$SCRATCH/comments/NEW_mg_residual_benchmark.md"
out="$SCRATCH/comments/_filled_NEW_mg_residual_benchmark.md"
python -c "
import sys
src = open(sys.argv[1]).read()
src = src.replace('<!-- COMMIT-HASH -->', sys.argv[2])
# Strip the 4 leading metadata HTML comments (they say 'New issue body — to be filed' etc.)
lines = src.split('\n')
while lines and lines[0].strip().startswith('<!--') and 'COMMIT-HASH' not in lines[0]:
    lines.pop(0)
while lines and lines[0].strip() == '':
    lines.pop(0)
open(sys.argv[3], 'w').write('\n'.join(lines))
" "$in" "$REPLACEMENT" "$out"

gh issue edit 140 --body-file "$out" || echo "  (edit failed for #140)"

echo "DONE."
