export const meta = {
  name: 'docs-pass',
  description: 'Per section: a Sonnet archivist drafts in one worktree per section, a Fable qa re-runs every [M] claim; the orchestrator integrates and builds Sphinx',
  whenToUse: 'A theory page, or a set of sections, made present-tense false by a carve. Plan: .claude/plans/workflow_orchestration.md',
  phases: [
    { title: 'Draft', detail: 'executor archivist, one per section', model: 'sonnet' },
    { title: 'Verify', detail: 'advisor qa re-runs each measured claim', model: 'fable' },
  ],
}

// args = { campaign, main, sections: [{ key, branch, brief: '<page, section anchors, what changed, the code to read>' }], rules }
const A = args || {}
const results = await pipeline(
  A.sections || [],
  (s) => agent(
    `Workflow: ${A.campaign}, phase Draft (docs, section ${s.key}). ${s.brief}\n${A.rules || ''}\nWork only in your worktree (print pwd and branch first; rename the branch to ${s.branch}). The articulation rule binds: full derivations, the design rationale, what was tried and failed, the numerical evidence with [M]; a present-tense-false sentence is a bug. Do not run Sphinx (the orchestrator does) and never hand-edit generated files. Commit (stage by path; the trailers in the rules). Return at most 300 words: the worktree, the branch, the commit, the sections touched, every [M] you wrote with its command.`,
    { label: `docs:${s.key}`, phase: 'Draft', isolation: 'worktree', model: 'sonnet', agentType: 'archivist' },
  ),
  (draft, s) => draft && agent(
    `Workflow: ${A.campaign}, phase Verify (docs, section ${s.key}). READ-ONLY. The archivist's report:\n${draft}\nFor every [M] claim in its diff, re-run the command on a copy and compare; for every present-tense claim about code, read the code. Return at most 300 words: the claims checked (k of N), each mismatch with its file:line and the correct text.`,
    { label: `verify:${s.key}`, phase: 'Verify', model: 'fable', agentType: 'qa' },
  ).then((v) => ({ section: s.key, draft, verify: v })),
)
return results.filter(Boolean)
