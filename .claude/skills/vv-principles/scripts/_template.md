<!-- Template for a vv-principles case study.
     Copy this file when adding a new entry. File name pattern:
     err<NNN>_<one_or_two_word_principle>.md     for ERR-NNN bugs
     issue<NN>_<one_or_two_word_principle>.md    for issue-tracked failures

     Drop the "<…>" placeholders and fill in the four bullets.
     The bar: every entry MUST answer bullet 4 sharply. If you cannot
     answer "what evidence class would have caught this," do not write
     the file — the case is catalog material, not skill material. -->

# <ERR-NNN | Issue #NN> — <epistemic-failure principle in five words>

**The bug.** <one to three sentences. What was wrong, in plain language.>

**Evidence that existed before it was caught.** <what tests passed;
what made the bug look like correctness; what conventional standard of
"verified" was satisfied; how many people / methods / cross-checks
endorsed the result before falsification.>

**Why that evidence didn't catch it.** <the epistemic mechanism. Be
specific: which structural property (independence, ladder level,
quadrature stability, BC fidelity, mode coverage) the existing tests
shared with the bug or failed to exercise. This bullet is about
*why the test apparatus was blind*, not about the bug's algebraic
form.>

**What evidence class would have caught it.** <THE LOAD-BEARING
BULLET. Use the redirect pattern: "X is necessary, NEVER sufficient
— instead, require Y." Name the specific evidence class. Cite the
relevant section of vv-principles SKILL.md or reference.md (e.g.
"§4.3 MMS rule 2"; "§1 structural-vs-procedural independence";
"§4.5 reference contamination"). If the bullet ends up vague, the
case study is not yet ready — sharpen until the prescription is
operational.>

**References.** <ERR-NNN catalog entry; the @pytest.mark.catches test
that gates against regression; the originating issue / commit / PR;
relevant agent memory or lessons.md entries.>
