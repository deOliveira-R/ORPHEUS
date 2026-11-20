# Memory distillation, 2026-11-16: the brief every owner receives

Workflow: harness update, memory distillation. Phase: the owner PROPOSES; the
orchestrator evaluates and applies. You are distilling YOUR OWN memory,
`.claude/agent-memory/<you>/`. You change nothing there in this dispatch: every
output is a proposal file under `scratch/memory_distillation/<you>/`.

## Why, with the measurements

- Your `MEMORY.md` index loads whole into every dispatch of you. Your
  `lessons.md` is opened in a fraction of your dispatches (`[M]` 2026-11-16, a
  count over this project's stored sub-agent transcripts, a lower bound: qa 47
  of 179, elegance-enforcer 36 of 145, test-architect 77 of 126, explorer 97 of
  123, archivist 68 of 112, cross-domain-attacker 33 of 40,
  numerics-investigator 17 of 30). Both are read in full when read, at about
  2.6 characters per token.
- The digests were cut to at most 400 lines on 2026-09-22 and have regrown
  (525 to 800 lines for six owners). Several indexes carry hooks of 40 words
  where 15 decide relevance.

## The standard (the user's, 2026-06-21, sharpened 2026-08-03 and 2026-09-22)

1. The test for one entry: *what mistake was made, and what did I learn that
   changed how I work?* An entry that answers it stays. Campaign narrative (a
   wave, a phase, an issue number as the subject) is noise: the name points into
   git and brings nothing forward.
2. Compress only where it is lossless (cut prose, keep the behaviour and the
   pointer). Where it would be lossy, distil: find the one lesson that drives
   several entries and write THAT, keeping every distinct `check:`.
3. An entry that a rule, a preloaded skill or your own definition already
   carries is not memory (the workflows rule, invariant 6). It retires, naming
   the file and the clause. Verify before claiming: grep the carrier for the
   entry's most distinctive token, and read the entry's LAST paragraph, where a
   body states what still holds (`retirement-audit` item 19). An entry carried
   only in part keeps the uncarried part.
4. The index holds hooks of at most 15 words, enough to decide whether to open
   the file, never the content. A cold file stays cold: no file is deleted or
   renamed by this pass, and no lesson ID is renumbered (IDs are cited from
   other files).
5. A lesson you apply on essentially every task is not memory either: propose
   it for your definition (one standing sentence, no war story). A correction
   that generalises across agents is a rule or skill candidate: name the target
   file. Both go to the user as candidates; be conservative, the always-on
   block is over its target.

## What you write, under `scratch/memory_distillation/<you>/`

- `dispositions.md`: one row per entry of your current `lessons.md` and one per
  line or block of your current `MEMORY.md`. Columns: the entry's ID or first
  words; one disposition of KEEP, TIGHTEN (lossless), MERGE (into which ID),
  COVERED (file and clause, with the token you grepped and its hit count),
  ARCHIVE (the war-story part moves to the cold archive), PROMOTE (the proposed
  definition sentence), UPLIFT (the target rule or skill and the proposed
  sentence); one line of reason. End with the counts per disposition and the
  line and character counts before and after.
- `lessons.md` and `MEMORY.md`: the proposed files, complete, ready to replace
  the current ones.
- `archive_append.md`: the exact text to append to your cold archive, if any
  (text that leaves the digest and exists nowhere else).
- `candidates.md`: the PROMOTE and UPLIFT rows in full, each with its evidence.

## Bounds

- Read your whole `lessons.md` and `MEMORY.md`; open an archive entry or a
  topic file only to decide a row. Read the rules in `.claude/rules/`, your
  preloaded skills in `.claude/skills/<name>/SKILL.md` and your definition in
  `.claude/agents/<you>/AGENT.md` as the carriers.
- No size target overrides rule 2: a digest that is honestly 450 lines is
  right, a 300-line one that dropped a `check:` is wrong.
- Edit nothing outside `scratch/memory_distillation/<you>/`. Spawn no agent.
- Return at most 200 words: the counts, the three rows you were least sure
  of, and a `NEEDS:` block.

## Settled by the pilot (qa, 2026-11-16)

- Every row of a "retired into the rules" table in your digest is an entry:
  re-verify its carrier claim (the pilot found four false ones).
- A clause carried in paraphrase, with 0 hits for the entry's token: quote the
  carrier's sentence in the row beside the 0.
- A path-scoped rule (front matter `paths:`) counts as a carrier only for work
  on those paths; say so in the row.
- A recipe naming a symbol or a path that no longer exists: correct it against
  the tree and mark it `[R]` if you did not run it; never copy it forward stale.
- The carrier census is a Python script (`re` + `pathlib`) with one positive
  control printed before any zero is believed; a zsh loop over an unsplit
  variable reads 0 silently.
- If your proposed index or digest cites text that `archive_append.md` adds,
  say so in `NEEDS:`.
- An honest result may be small: the pilot's digest went from 734 to 689 lines
  and its index from 85 to 29.
