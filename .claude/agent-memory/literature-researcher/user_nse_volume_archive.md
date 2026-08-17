---
name: user-nse-volume-archive
description: The user's COMPLETE NSE run lives at /Users/rodrigo/Downloads/NSE as per-volume zip bundles — a far richer Tier-0 surface than scratch/literature; grep the zip listings by title before ever reporting an NSE paper as not-local
metadata:
  type: reference
---

# Tier 0 is BIGGER than `scratch/literature/` — the full NSE run is on disk

`/Users/rodrigo/Downloads/NSE/` holds the user's complete Nuclear Science &
Engineering archive as **per-volume zip bundles**:

```
Vol_<NNN>(<issue>)_Nuclear Science and Engineering.zip
  → NN-<Article Title Truncated To ~95 chars>.pdf
```

(Also present: `Taylor_&_Francis_Articles(25Apr2026)-*.zip` staging bundles —
the archive was bulk-downloaded from T&F, which is NSE's publisher.)

**Why this matters:** `scratch/literature/` holds only ~70 hand-picked PDFs.
The brief's standing claim that *"the user has ALL Nuclear Science &
Engineering volumes locally"* refers to THIS directory, not to
`scratch/literature/`. Reporting an NSE paper as "not in local folder" without
checking here is a false negative.

## The lookup (read-only, no extraction needed to search)

```bash
cd /Users/rodrigo/Downloads/NSE && for z in *.zip; do
  unzip -l "$z" 2>/dev/null | sed "s|^|$z :: |"
done 2>/dev/null | grep -i -E "<title fragment>"
```

⚠ **Filenames are TRUNCATED** (~95 chars), so a grep for a word late in a long
title silently misses. Grep a **distinctive early fragment**; e.g. `curvilinear`
missed `…Calculation in Curvili.pdf`. Same failure family as the
shortest-distinctive-fragment rule for message strings.

## Acquisition into the working folder

```bash
unzip -j -o "Vol_182(4)_Nuclear Science and Engineering.zip" "*<fragment>*" -d /tmp/nseext
cp "/tmp/nseext/<file>.pdf" "scratch/literature/<Author>(<Year>)<Full Title>.pdf"
.venv/bin/python tools/ocr_literature.py --glob '<Author>*'
```

Rename to the folder's `Author(Year)Title.pdf` convention on the way in — the
archive's `NN-Title` names lose authorship and are truncated. OCR is ~$0.01–0.02
per paper; cache-idempotent.

**Worked 2026-08-12 (Q68):** found Chaland & Samba 2016 NSE 182 and
Wu-Xie-Fischer 1999 NSE 133 this way after CrossRef named them; both added to
`scratch/literature/` and OCR'd. Chaland–Samba overturned a negative I had
already written. See [[lathrop-2000-angular-scheme-taxonomy]].

## Companion: pick the right database to NAME the paper first

⛔ **OpenAlex free-text is useless for `S_N` topics** — `SN` matches *supernova*;
7 formulations returned WMAP / Type Ia / MHD / sea-ice results. Use
**ISSN-scoped CrossRef `journal_search`** (NSE `0029-5639`, ANE `0306-4549`,
JCP `0021-9991`, TTSP `0041-1450`, JQSRT `0022-4073`), then look the title up in
this archive. Also: **Semantic Scholar `get_citations` on a seed DOI** is the
best net for "what came after paper X" — it found every post-2010 follow-up that
title search missed. Semantic Scholar `429`-rate-limits after ~2 rapid calls.
