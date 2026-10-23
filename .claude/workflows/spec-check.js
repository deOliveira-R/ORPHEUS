export const meta = {
  name: 'spec-check',
  description: 'Sonnet measures the open items in parallel, a Sonnet drafter writes one specification from them, then an Opus elegance cold read and an Opus structural check',
  whenToUse: 'Turning ruled design decisions into one executable specification before gates are written. Plan: .claude/plans/workflow_orchestration.md',
  phases: [
    { title: 'Measure', detail: 'executors, one per open item', model: 'sonnet' },
    { title: 'Draft', detail: 'executor writes the specification', model: 'sonnet' },
    { title: 'Check', detail: 'reviewer elegance cold read and reviewer structural check', model: 'opus' },
  ],
}

// args = {
//   campaign, out: '/abs/dir for measurement files', spec: '/abs/path of the specification to write',
//   inputs: '<files the specification is built from + the rulings that bind>',
//   measure: [{ key, question: '<the open item, with its control>' }],
//   contents: '<what the specification must fix>',
// }
const A = args || {}
const measured = await parallel((A.measure || []).map((m) => () => agent(
  `Workflow: ${A.campaign}, phase Measure (${m.key}). READ-ONLY on tracked files; write under ${A.out}/. ${A.inputs}\n\nMeasure: ${m.question}\nInstrument-doctrine: the instrument, the population, a positive control; AST for membership, Python re for text. Write ${A.out}/measure_${m.key}.md; return at most 300 words: the answer with its [M], and what it decides.`,
  { label: `measure:${m.key}`, phase: 'Measure', model: 'sonnet' },
)))

const draft = await agent(
  `Workflow: ${A.campaign}, phase Draft. Write ONE specification to ${A.spec} (the only tracked file you write; do not commit). Inputs: ${A.inputs}\n\nThe measurements (files under ${A.out}/measure_*.md):\n${measured.filter(Boolean).join('\n---\n')}\n\nIt must fix: ${A.contents}\nWhere the inputs disagree, the rulings decide; where no ruling decides, choose, give the structural reason, and list the choice under "Decisions for the orchestrator". Every claim about the tree [M] with its instrument. Return at most 300 words: the outline, the decisions you took, NEEDS.`,
  { label: 'draft:spec', phase: 'Draft', model: 'sonnet' },
)

const checks = await parallel([
  () => agent(`Workflow: ${A.campaign}, phase Check (elegance, COLD read: you have not seen the discussion). READ-ONLY. Read ${A.spec} alone first, then its inputs. ${A.inputs}\nCan every step be executed from the text alone? Twin paths, a repeated conditional that is a missing type, an illegal state representable, a name that collides, a verb on the wrong object. Write ${A.out}/check_elegance.md; return at most 400 words, findings ranked with the fix.`, { label: 'check:elegance', phase: 'Check', model: 'opus', agentType: 'elegance-enforcer' }),
  () => agent(`Workflow: ${A.campaign}, phase Check (structure). READ-ONLY. Read ${A.spec} and its inputs. ${A.inputs}\nDoes each type model the native mathematical structure; is every law stated (adjoint, composition, invariance); does each gate red on a state in the tree; is the migration order leaves first with no orphan? Write ${A.out}/check_structure.md; return at most 400 words, findings ranked, an UNEXPLORED block.`, { label: 'check:structure', phase: 'Check', model: 'opus', agentType: 'cross-domain-attacker' }),
])

return { measured: measured.filter(Boolean), draft, checks: checks.filter(Boolean) }
