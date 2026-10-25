export const meta = {
  name: 'docs-pass',
  description: 'Per section: a Sonnet archivist drafts in a worktree the orchestrator cut from the named base, an Opus qa verifies read-only against the code, a Sonnet archivist amends; the orchestrator integrates and builds Sphinx',
  whenToUse: 'A theory page, or a set of sections, made present-tense false by a carve. The orchestrator first creates one worktree per section on its own branch FROM args.base (never origin/main: the carve is usually unmerged). Plan: .claude/plans/workflow_orchestration.md',
  phases: [
    { title: 'Draft', detail: 'executor archivist, one per section, in its pre-created worktree', model: 'sonnet' },
    { title: 'Verify', detail: 'reviewer qa, read-only: re-runs each [M], reads the code behind each present-tense claim', model: 'opus' },
    { title: 'Amend', detail: 'executor archivist applies the findings (only when findings remain)', model: 'sonnet' },
  ],
}

// args = { campaign, base, rules, sections: [{ key, worktree, branch, brief }] }
//   base      the commit every worktree was cut from (the positive control below checks it)
//   worktree  an absolute path the orchestrator created with `git worktree add -b <branch> <path> <base>`
//   rules     the shared rules block: tense rule, label rule, file ownership, generated files, Python path, escalation
// Why no `isolation: 'worktree'`: it branches from the default branch, so on an unmerged carve every unit
// documents the WRONG tree (WF2.5, 2026-10-24, caught before launch). Why an Amend stage: in WF2.5 qa raised
// 108 findings over 14 units and every unit needed one (plan wave2_boundary_carve.md section 11.8).
const A = args || {}
const FINDINGS = {
  type: 'object',
  properties: {
    verdict: { type: 'string', enum: ['sound', 'amend'] },
    claims_checked: { type: 'string' },
    findings: { type: 'array', items: { type: 'object', properties: {
      severity: { type: 'string', enum: ['high', 'medium', 'low'] },
      where: { type: 'string' }, finding: { type: 'string' }, evidence: { type: 'string' }, fix: { type: 'string' },
    }, required: ['severity', 'where', 'finding', 'evidence', 'fix'] } },
  },
  required: ['verdict', 'claims_checked', 'findings'],
}
const CONTROL = (s) => `Work only in ${s.worktree} (branch ${s.branch}). Positive control before any edit: \`git -C ${s.worktree} log -1 --format=%h\` prints ${A.base} (or a commit of yours on top of it), and a probe run from OUTSIDE the repository with PYTHONPATH=${s.worktree} prints an \`orpheus.__file__\` inside ${s.worktree}. If either fails, stop and return the reading in NEEDS.`
const COMMIT = `Commit: stage by explicit path after reading \`git status --porcelain\`; \`git commit -F -\` from a quoted heredoc; Conventional Commits; the Co-Authored-By trailer names the model that actually wrote the commit (yours), followed by the session trailer in the rules.`
const results = await pipeline(
  A.sections || [],
  (s) => agent(
    `Workflow: ${A.campaign}, phase Draft (docs, section ${s.key}). ${CONTROL(s)}\n${s.brief}\n${A.rules || ''}\nThe articulation rule binds: full derivations, the design rationale, what was tried and failed, the numerical evidence with [M] and its command; a present-tense-false sentence is a bug. A code block you publish is RUN, and its printed output is JUDGED against the physics before it is pasted (WF2.5 pasted a defect's k as a result). Do not run Sphinx (the orchestrator does) and never hand-edit generated files. ${COMMIT} Return at most 400 words: the commit, the sections touched, every [M] with its command, the label counts before and after, and NEEDS.`,
    { label: `draft:${s.key}`, phase: 'Draft', model: 'sonnet', agentType: 'archivist' },
  ),
  (draft, s) => draft && agent(
    `Workflow: ${A.campaign}, phase Verify (docs, section ${s.key}). READ-ONLY on ${s.worktree}; probes go in a fresh directory named for this review. The archivist's report:\n${draft}\nRead the diff \`git -C ${s.worktree} diff ${A.base}\`. For every [M] claim re-run the command and compare; for every present-tense claim about code, read the code; for every published code block, run it and ask whether its output is physically right; check the tense of every surviving mention of a retired name, and that no label was removed. ${A.rules || ''}\nReturn the findings object; claims_checked states k of N per kind.`,
    { label: `verify:${s.key}`, phase: 'Verify', model: 'opus', agentType: 'qa', schema: FINDINGS },
  ).then((v) => ({ draft, v })),
  (r, s) => {
    if (!r) return null
    if (!r.v || r.v.verdict === 'sound' || !r.v.findings.length) return { section: s.key, draft: r.draft, verify: r.v, amend: null }
    return agent(
      `Workflow: ${A.campaign}, phase Amend (docs, section ${s.key}). ${CONTROL(s)}\nqa's findings on your draft:\n${JSON.stringify(r.v.findings)}\nApply each finding inside your scope, or decline it with the structural reason; a fix outside your scope goes to NEEDS with the exact text. ${A.rules || ''}\n${COMMIT} Return at most 400 words: one disposition per finding, the commit, the label counts, and NEEDS.`,
      { label: `amend:${s.key}`, phase: 'Amend', model: 'sonnet', agentType: 'archivist' },
    ).then((a) => ({ section: s.key, draft: r.draft, verify: r.v, amend: a }))
  },
)
return results.filter(Boolean)
