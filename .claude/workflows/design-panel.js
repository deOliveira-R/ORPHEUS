export const meta = {
  name: 'design-panel',
  description: 'Sonnet drafts from distinct angles in parallel, then Opus attackers each read every draft; the orchestrator synthesizes',
  whenToUse: 'An ontology being searched: two or more independent designs, then adversarial attacks, before the user rules. Plan: .claude/plans/workflow_orchestration.md',
  phases: [
    { title: 'Draft', detail: 'executor designers, one per angle', model: 'sonnet' },
    { title: 'Attack', detail: 'reviewer attackers, one per lens, reading every draft', model: 'opus' },
  ],
}

// args = {
//   campaign, out: '/abs/dir for the files', context: '<files to read + the rulings that bind, not re-litigated>',
//   angles: [{ key, angle: '<where this designer starts from>' }],          // 2-3
//   lenses: [{ key, agentType: 'cross-domain-attacker' | 'elegance-enforcer' | 'qa', lens: '<what to attack>' }],
// }
const A = args || {}
if (!A.angles || A.angles.length < 2) throw new Error('design-panel: give at least two angles')
const READONLY = `READ-ONLY on tracked files; write ONLY under ${A.out}/.`
const DELIVER = `Write your design to the named file (markdown, maximal precision: every type with its fields and invariant, every module path, every verb with its signature, the algebra in symbols, the migration order leaves first with the dependency table, the gates that make each illegal state unrepresentable with its first red, what you left out with the structural reason; every claim about the tree [M] with its instrument). Return at most 400 words: the path, the design in brief, the three decisions you were least sure of, NEEDS.`

const drafts = await parallel(A.angles.map((a) => () => agent(
  `Workflow: ${A.campaign}, phase Draft (design ${a.key}). ${A.context}\n${READONLY}\n\nYour angle: ${a.angle}\n\n${DELIVER}\nFile: ${A.out}/design_${a.key}.md`,
  { label: `draft:${a.key}`, phase: 'Draft', model: 'sonnet' },
)))
const files = A.angles.map((a) => `${A.out}/design_${a.key}.md`).join(', ')
log(`${drafts.filter(Boolean).length} of ${A.angles.length} drafts returned`)

const attacks = await parallel((A.lenses || []).map((l) => () => agent(
  `Workflow: ${A.campaign}, phase Attack (${l.key}). First pass adversarial: how would I break this; how would I make it 100x better. Read every design: ${files}, and the context they cite. ${A.context}\n${READONLY}\n\nYour lens: ${l.lens}\nWhere the designs disagree, say which is right and why. Write ${A.out}/attack_${l.key}.md; return at most 400 words: findings ranked, the winner on each disagreement, the decisions that need the user, an UNEXPLORED block.`,
  { label: `attack:${l.key}`, phase: 'Attack', model: 'opus', agentType: l.agentType },
)))

return { drafts: drafts.filter(Boolean), attacks: attacks.filter(Boolean), files }
