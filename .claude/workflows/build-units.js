export const meta = {
  name: 'build-units',
  description: 'Per unit, in its own worktree: Sonnet gates (red) then build, Fable qa + elegance in parallel, Sonnet amend, Fable verify',
  whenToUse: 'Additive or file-disjoint units of a planned change, built in parallel (layer DAG waves; any build-review-amend fan-out). Plan: .claude/plans/workflow_orchestration.md',
  phases: [
    { title: 'Gates', detail: 'executor test-architect writes the done-when gates, red at the base tip', model: 'sonnet' },
    { title: 'Build', detail: 'executor implementer turns the gates green', model: 'sonnet' },
    { title: 'Review', detail: 'advisor qa and elegance-enforcer in parallel', model: 'fable' },
    { title: 'Amend', detail: 'executor applies the inline-fixable findings', model: 'sonnet' },
    { title: 'Verify', detail: 'advisor checks each finding was answered', model: 'fable' },
  ],
}

// args = {
//   campaign: 'layer DAG wave 2',          // named in every brief
//   main: '/abs/path/to/main/tree',         // read-only reference tree (scratch, plans)
//   plan: '.claude/plans/<plan>.md#<section>',
//   rules: '<the shared rules block: git, tests, residue spellings, rulings that bind>',
//   units: [{ key, branch, brief, gates?: '<gate spec, if the test-architect already wrote it>',
//             worktree?: '<abs path of an EXISTING worktree to continue in (sequential slices on one branch); no new worktree is made>' }],
//   skip_gates?: false,                     // true when every unit's brief carries its gates
// }
const A = args || {}
if (!A.units || !A.units.length) throw new Error('build-units: args.units is empty')
const MAIN = A.main
const PY = `${MAIN}/.venv/bin/python`

const RULES = `
${A.rules || ''}
Invariants (.claude/plans/workflow_orchestration.md):
- Work ONLY in your worktree; print \`pwd\` and \`git rev-parse --abbrev-ref HEAD\` first. You may READ ${MAIN} (plans, scratch); never edit it.
- Commit on your worktree's branch; stage by explicit path; read \`git status --porcelain\` first; \`git commit -F -\` from a quoted heredoc, Conventional Commits, ending with:
Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01BAyyA444HWpLA91BrSQ3LD
- Tests: \`${PY} -O -m pytest -p no:cacheprovider -q --color=no <paths>\` from your worktree root, TARGETED, serially; prove the import resolves to your tree (\`${PY} -c "import orpheus; print(orpheus.__file__)"\` from the root; else set PYTHONPATH and say so). No full suite, no Sphinx, no pyright claims (vacuous in worktrees, #640): the orchestrator runs them at integration.
- Never hand-edit generated files (CLAUDE.md, docs/theory/verification/matrix.rst, .claude/rules/*, .claude/skills/*).
- A question the brief does not settle goes to NEEDS; do the independent work; never message the orchestrator mid-run.
- Every claim with its population and instrument; a zero only beside a positive control.`

const BUILT = {
  type: 'object',
  properties: {
    worktree: { type: 'string' }, branch: { type: 'string' },
    commits: { type: 'array', items: { type: 'string' } },
    report: { type: 'string' }, needs: { type: 'string' },
  },
  required: ['worktree', 'branch', 'commits', 'report', 'needs'],
}
const FINDINGS = {
  type: 'object',
  properties: {
    verdict: { type: 'string', enum: ['sound', 'amend'] },
    findings: { type: 'array', items: { type: 'object', properties: {
      severity: { type: 'string', enum: ['high', 'medium', 'low'] },
      where: { type: 'string' }, finding: { type: 'string' }, evidence: { type: 'string' },
      fix: { type: 'string' }, inline_fixable: { type: 'boolean' },
    }, required: ['severity', 'where', 'finding', 'evidence', 'fix', 'inline_fixable'] } },
    rejected: { type: 'string' },
  },
  required: ['verdict', 'findings', 'rejected'],
}
const VERIFIED = {
  type: 'object',
  properties: {
    all_answered: { type: 'boolean' },
    unanswered: { type: 'array', items: { type: 'string' } },
    notes: { type: 'string' },
  },
  required: ['all_answered', 'unanswered', 'notes'],
}

const head = (u, phase) => `Workflow: ${A.campaign}, phase ${phase}, unit ${u.key}. Plan: ${A.plan}.`

const results = await pipeline(
  A.units,
  // Gates + Build in ONE isolated worktree: the test-architect writes the gates and records their red,
  // then the implementer (a second agent, same worktree by path) turns them green.
  async (u) => {
    const into = u.worktree ? `Work ONLY in the existing worktree ${u.worktree} (cd there first; print pwd and the branch, which must be ${u.branch}); it continues a branch an earlier unit left, so build on its tip. ` : ''
    const fresh = u.worktree ? {} : { isolation: 'worktree' }
    const gate = A.skip_gates ? null : await agent(
      `${head(u, 'Gates')} ${into}You are the test-architect. Write the unit's done-when gates as tests (and nothing else of the change), commit them, and record each gate's RED at the base tip (the command and the failing line). Each gate names the input that reddens it (X1); a structural repair gets a mutation you describe for the reviewers. Do not make them green.\n\nThe unit:\n${u.brief}\n${u.gates ? '\nThe gate specification:\n' + u.gates : ''}\n${RULES}\n${u.worktree ? '' : 'Rename your branch to ' + u.branch + ' (\`git branch -m\`). '}Return at most 300 words through the schema.`,
      { label: `gates:${u.key}`, phase: 'Gates', ...fresh, model: 'sonnet', agentType: 'test-architect', schema: BUILT },
    )
    if (!A.skip_gates && (!gate || !gate.worktree)) return { unit: u.key, gate, built: null }
    const where = gate ? `Work in the existing worktree ${gate.worktree} (cd there first), on branch ${gate.branch}, where the test-architect committed the gates (${gate.commits.join(', ')}); its report:\n${gate.report}` : (u.worktree ? into : `You are in a fresh worktree; rename its branch to ${u.branch}.`)
    const built = await agent(
      `${head(u, 'Build')} You are the implementer. ${where}\n\nTurn the gates green by making the change below, with no shim and no residue; do not weaken a gate (if one is wrong, say so in NEEDS and leave it red).\n\nThe unit:\n${u.brief}\n${RULES}\nReturn at most 500 words through the schema: commits (\`git log --oneline main..HEAD\`), what changed where, each gate red-then-green with its command, the tests run with counts, the residue census (instrument, population, control), what you did not do and why.`,
      gate ? { label: `build:${u.key}`, phase: 'Build', model: 'sonnet', schema: BUILT }
           : { label: `build:${u.key}`, phase: 'Build', ...fresh, model: 'sonnet', schema: BUILT },
    )
    return { unit: u.key, gate, built }
  },
  // Review: two advisors, independent, read-only on a copy.
  async (r, u) => {
    if (!r.built || !r.built.commits || !r.built.commits.length) return { ...r, reviews: [] }
    const ctx = `${head(u, 'Review')} The unit was built in worktree ${r.built.worktree}, branch ${r.built.branch}; review \`git -C ${r.built.worktree} diff main...${r.built.branch}\`. The unit:\n${u.brief}\n\nThe gates' report:\n${r.gate ? r.gate.report : '(gates in the brief)'}\n\nThe implementer's report:\n${r.built.report}\n\nREAD-ONLY: mutate only a copy (\`cp -a\` to a temp dir), print orpheus.__file__ before any probe, targeted tests only with PYTHONPATH set to the copy. First pass adversarial: how would I break this; how would I make it 100x better. Return through the schema, at most 500 words; each finding with file:line, evidence ([M] with the command) and a concrete fix; rejected candidates with one-line reasons.`
    const reviews = await parallel([
      () => agent(`${ctx}\nYour lens (qa): residue the change missed (every spelling class in the rules), whether each gate truly reddens (mutate in the copy), a surviving gate made vacuous or demoted (retirement-audit D.14-D.16), behaviour preserved (digests, memo, layer census), present-tense prose made false.`, { label: `qa:${u.key}`, phase: 'Review', model: 'fable', agentType: 'qa', schema: FINDINGS }),
      () => agent(`${ctx}\nYour lens (coding-elegance): each object in its right place with its right shape; twin paths; a repeated conditional that is a missing type; stringly-typed dispatch; an illegal state still representable; names; the code read aloud against the mathematics.`, { label: `elegance:${u.key}`, phase: 'Review', model: 'fable', agentType: 'elegance-enforcer', schema: FINDINGS }),
    ])
    return { ...r, reviews: reviews.filter(Boolean) }
  },
  // Amend (executor) then Verify (advisor), only when a reviewer found something inline-fixable.
  async (r, u) => {
    const todo = (r.reviews || []).flatMap((v) => v.findings || []).filter((f) => f.inline_fixable)
    if (!r.built || !todo.length) return { ...r, amended: null, verified: null }
    const amended = await agent(
      `${head(u, 'Amend')} Work ONLY in worktree ${r.built.worktree} on branch ${r.built.branch} (cd there; print pwd and branch). Apply each finding below as stated or with a better fix you justify; for one you reject, give the structural reason. Re-run the targeted tests; one commit listing each finding and its disposition.\n\nFindings:\n${JSON.stringify(todo, null, 1)}\n\nThe unit:\n${u.brief}\n${RULES}\nReturn at most 400 words: the commit, each finding's disposition, tests with counts, NEEDS.`,
      { label: `amend:${u.key}`, phase: 'Amend', model: 'sonnet' },
    )
    if (!amended) return { ...r, amended: null, verified: null }
    const verified = await agent(
      `${head(u, 'Verify')} READ-ONLY. In worktree ${r.built.worktree} (branch ${r.built.branch}), read the amend commit (the newest) and check, finding by finding, that each was answered: fixed as stated, fixed better, or rejected with a structural reason that holds. Findings:\n${JSON.stringify(todo, null, 1)}\n\nThe amender's report:\n${amended}\n\nReturn through the schema, at most 250 words.`,
      { label: `verify:${u.key}`, phase: 'Verify', model: 'fable', schema: VERIFIED },
    )
    return { ...r, amended, verified }
  },
)

return results.map((r) => r && ({
  unit: r.unit,
  worktree: r.built && r.built.worktree,
  branch: r.built && r.built.branch,
  commits: r.built && r.built.commits,
  gates: r.gate && r.gate.report,
  build: r.built && r.built.report,
  needs: [r.gate && r.gate.needs, r.built && r.built.needs].filter(Boolean).join('\n'),
  reviews: r.reviews,
  amended: r.amended,
  verified: r.verified,
}))
