export const meta = {
  name: 'spec-check',
  description: 'Sonnet measures the open items in parallel, a Sonnet drafter writes one specification from them, then an Opus elegance cold read and an Opus structural check',
  whenToUse: 'Turning ruled design decisions into one executable specification before gates are written. Plan: .claude/plans/workflow_orchestration.md',
  phases: [
    { title: 'Measure', detail: 'executors, one per open item', model: 'sonnet' },
    { title: 'Draft', detail: 'executor writes the specification', model: 'sonnet' },
    { title: 'Check', detail: 'reviewer elegance cold read and reviewer structural check, of the whole text or of each unit of the final text', model: 'opus' },
  ],
}

// args = {
//   campaign, out: '/abs/dir for measurement files', spec: '/abs/path of the specification to write',
//   inputs: '<files the specification is built from + the rulings that bind>',
//   measure: [{ key, question: '<the open item, with its control>' }],
//   contents: '<what the specification must fix>',
//   units?: [{ key, sections: '<the unit's sections, e.g. "2.5 and 5.4 (with decision 49 of section 9)">' }],
//   check_only?: false,   // true: no Measure and no Draft; the Check phase reads the text at args.spec as it stands
// }
// The check that pays is per unit, on the text the builders will receive ([M] 2026-11-15, .claude/plans/harness_eval.md, the
// specification rung): of 8 defects per unit that the code review later traced to the specification, one run of each lens
// on ONE unit's sections of the final draft named 5.3 and 6.4 together, where the whole-text check of the first draft had
// caught none of the 14 mandated ones and three were introduced by amendments nobody re-read. So: pass `units`, and after
// the LAST amendment of the text run this template again with `check_only: true` before the gates stage.
const A = args || {}
const UNITS = A.units && A.units.length ? A.units : [null]
const measured = A.check_only ? [] : await parallel((A.measure || []).map((m) => () => agent(
  `Workflow: ${A.campaign}, phase Measure (${m.key}). READ-ONLY on tracked files; write under ${A.out}/. ${A.inputs}\n\nMeasure: ${m.question}\nInstrument-doctrine: the instrument, the population, a positive control; AST for membership, Python re for text. Write ${A.out}/measure_${m.key}.md; return at most 300 words: the answer with its [M], and what it decides.`,
  { label: `measure:${m.key}`, phase: 'Measure', model: 'sonnet', agentType: 'executor' },
)))

const draft = A.check_only ? '(check only: the text at args.spec as it stands)' : await agent(
  `Workflow: ${A.campaign}, phase Draft. Write ONE specification to ${A.spec} (the only tracked file you write; do not commit). Inputs: ${A.inputs}\n\nThe measurements (files under ${A.out}/measure_*.md):\n${measured.filter(Boolean).join('\n---\n')}\n\nIt must fix: ${A.contents}\nWhere the inputs disagree, the rulings decide; where no ruling decides, choose, give the structural reason, and list the choice under "Decisions for the orchestrator". Every claim about the tree [M] with its instrument. Return at most 300 words: the outline, the decisions you took, NEEDS.`,
  { label: 'draft:spec', phase: 'Draft', model: 'sonnet', agentType: 'executor' },
)

const scope = (u) => u ? `Your scope is unit ${u.key}: sections ${u.sections}, and whatever those sections cite; the other units have their own readers. ` : ''
const file = (name, u) => `${A.out}/check_${name}${u ? '_' + u.key : ''}.md`
const checks = await parallel(UNITS.flatMap((u) => [
  () => agent(`Workflow: ${A.campaign}, phase Check (elegance, COLD read: you have not seen the discussion). READ-ONLY. ${scope(u)}Read ${A.spec} alone first, then its inputs. ${A.inputs}\nCan every step be executed from the text alone? Twin paths, a repeated conditional that is a missing type, an illegal state representable, a name that collides, a verb on the wrong object. Write ${file('elegance', u)}: findings ranked, each with the line it concerns (quoted) and the fix; return at most 400 words.`, { label: `check:elegance${u ? ':' + u.key : ''}`, phase: 'Check', model: 'opus', agentType: 'elegance-enforcer' }),
  () => agent(`Workflow: ${A.campaign}, phase Check (structure). READ-ONLY. ${scope(u)}Read ${A.spec} and its inputs. ${A.inputs}\nDoes each type model the native mathematical structure; is every law stated (adjoint, composition, invariance); does each gate red on a state in the tree; is the migration order leaves first with no orphan? Write ${file('structure', u)}: findings ranked, each with the line it concerns (quoted) and the fix, an UNEXPLORED block; return at most 400 words.`, { label: `check:structure${u ? ':' + u.key : ''}`, phase: 'Check', model: 'opus', agentType: 'cross-domain-attacker' }),
]))

return { measured: measured.filter(Boolean), draft, checks: checks.filter(Boolean) }
