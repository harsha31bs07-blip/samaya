# Plan to the grand finale (from 2026-10-13)

The idea submission for PS 26119 went in on 29 September; screening results are expected in late
October or November, and the grand finale (a 36-hour software sprint) is expected in December. No
official finale date is out yet: when it is, shift the calendar in section 6 and keep the order.

This plan does not replace `docs/POST_SUBMISSION_PLAN.md` (W0–W8, with their methods, tests and
gates). It updates that plan's status after 7–9 October, reorders the work by what decides a
finale, and adds:
- how the work is split between Claude and Codex, both working in parallel (sections 4 and 5);
- preparation for the sprint itself: the judges' rounds, a likely on-site change request, and a
  rehearsal.

The rules in `CLAUDE.md` hold for every item, whoever writes the code: differential tests, a
planted bug that the test catches, the debug, release and asan presets (and cuda for GPU code)
with zero warnings, every number through `src/verify/`, tolerances never loosened.

## 1. What wins this finale

Expect 4–5 teams on PS 26119 at the nodal centre. MRPL's people and domain experts judge over
several rounds in the 36 hours. They will look for, in this order:

1. **Does it solve MRPL's problem?** Their models are planning LPs, blending with quality
   recursion (PIMS-style SLP), crude and berth scheduling MILPs, and utilities. A judge will ask
   "can I run *my* model, and what do I get back?"
2. **The problem statement's explicit asks:** GPU acceleration, "thousands to millions of
   variables", numerical robustness, and a comparison with an established solver.
3. **Trust:** answers that are checked. The verifier is our strongest point that no other team
   is likely to have; keep it at the front.
4. **A live demo that does not fail,** and visible progress between rounds.
5. **Defending the internals:** every member can explain the sparse LU, the simplex, the
   barrier, PDLP, branch-and-cut and the verifier, and say honestly where we are behind.

So the order of the work is: scale and GPU (one headline run), refinery fit (ranging, SLP,
adoption through Python), MILP strength, then speed and polish.

## 2. Status on 13 October

Done since the plan of 7 October (branch `claude/kind-hamilton-jr2b6u`, on GitHub):
- **W1 crossover:** `--crossover` from barrier or PDLP; Netlib 93/93 verified by barrier +
  crossover and by PDLP + crossover (b4858b5, 4ac6de6).
- **W2 GPU residuals on the device and tiled SpMV:** on the 467k × 1M refinery LP, 76,160 →
  238,784 iterations in 600 s (3.1×), still not converged to 1e-8 (0d1b3c4, 3f63fb2).
- **Time limits** hold through presolve, PDLP setup and fallbacks (22b767b).
- **W3 step 0 attribution:** presolve is HiGHS's largest single gain; SCIP loses 7 of 9 without
  cuts; conflict analysis carries nothing (f11eb1e).

On local branches, not yet merged:
- `claude/amd`: AMD approximate degrees, barrier ordering 16.2 s → 1.9 s on the 400 × 400
  transportation LP (W4 item 4).
- `claude/pdlp-logs` (contains `claude/amd`): PDLP's per-column certificate lines moved to log
  level 3.

Open from POST_SUBMISSION_PLAN: the large LP to convergence, the Mittelmann subset, W3 C–G, W4,
W5, W6, W7 and W8.

## 3. Workstreams, by priority

Effort is in working days. Each item says who writes it: **[Claude]**, **[Codex tier/effort]**
from a spec Claude writes, or **[machine]** for unattended runs on the laptop at night. Claude
reviews and integrates every Codex result.

### P0. Setup this week (0.5 day)

1. **One source of truth** [Claude]. Edits happen in `C:\dev\samaya\mrpl-pro`, the Windows
   clone. WSL `~/mrpl-pro` stays the place for the Linux presets, CUDA and benchmarks: add the
   Windows clone as a remote there and fetch branches from it instead of committing in both
   places. Retire the WSL worktrees as their branches are settled (four have uncommitted
   changes; look at them first).
2. **Merge `claude/pdlp-logs`** (AMD ordering + log levels) [Claude], after a review
   [Codex sol/high] and the presets plus cuda on the merged tree.
3. **W0 housekeeping** from POST_SUBMISSION_PLAN, including the slide fix "Express / CEPLEX"
   [Codex luna/low].

### P1. Headline scale run (machine 2–3 nights, 1.5 days of Claude's work)

The problem statement asks for millions of variables on the GPU. One verified run above 10⁶
columns, quoted on one slide, answers it.

1. **The 467k × 1M refinery LP with `--lp-method pdlp --gpu --crossover`** [machine]: PDLP now
   stops at 1e-4 and the simplex finishes. One run at 3,600 s, alone on the machine.
2. **If the crossover is too slow there:** the same LP by barrier (AMD ordering now merged) and
   by the dual simplex [machine], so at least one method gives a verified optimum.
3. **The Mittelmann subset** (`bench/mittelmann_lp.test`) [machine]: dual simplex, barrier,
   PDLP on CPU and GPU, HiGHS; transfer times stated separately.
4. **Write `docs/results/large_lp.md`** from these runs [Claude], including where the CPU still
   wins and the runs that did not finish.

**Gate:** at least one verified solve above 10⁶ columns, or a plain statement in the results
of why not and what the best verified size is.

### P2. Refinery fit (about 2 weeks)

1. **Sensitivity ranging** (W5.3), 2 days.
   - The method, the API and the finite-difference test plan [Claude spec].
   - The ranging code and the tests [Codex sol/high]; it is textbook basis arithmetic.
   - The CLI flag `--ranging`, the JSON fields, Studio and the Excel export, together with the
     column bounds the tank capacity lines need [Codex sol/medium].
2. **Python adoption layer**, 1 day [Codex sol/medium]. `python/samaya_pulp.py` and
   `python/samaya_pyomo.py`: solver plugins that write MPS, call the `samaya` CLI and read its
   JSON back, including the verifier's status. Demo line: "your existing Pyomo model, one line
   changed".
3. **SLP / distributive recursion** (W5.1), 5–6 days.
   - Design, one second opinion on it [Codex astra/high, read-only], then the loop, the trust
     region and the final nonlinear check [Claude].
   - Haverly pooling instances with their known optima, and the generated multi-crude blending
     case [Codex sol/medium].
4. **QP blending case** (W5.2) [Codex sol/medium] and **Studio tank capacity lines and
   utilization** (W5.4) [Codex sol/medium].
5. **Realistic cases** (W5.5): `cases/mrpl.py` at full-year planning and a 4-week, 20-cargo
   schedule [Codex sol/medium]; Claude checks the models are sensible.

### P3. MILP strength (about 2.5 weeks)

In the order of the attribution results:
1. **C: presolve aggregator** (implied-free columns, doubleton equations, substitution with a
   fill limit), 2–3 days [Claude]: postsolve mistakes give wrong answers, so it stays with the
   orchestrator. Design opinion and branch review [Codex astra/high].
2. **SCIP's separators one at a time** [machine, 1 night], then zero-half cuts (F) or stronger
   c-MIR aggregation, whichever they point to [Codex sol/high from Claude's spec].
3. **D: clique table** in propagation and probing, implied-bound cuts, 1.5 days
   [Codex sol/high from Claude's spec; Claude integrates and runs the A/B].
4. **E: inference history** in branching, 0.5 day [Codex sol/medium; A/B by Claude].
5. **MIPLIB 60 s and 600 s re-run** after C and after D [machine]; `docs/results/comparison.md`
   [Claude].

G (symmetry) [Claude] and A4 (shifting, zi-round) [Codex sol/high] only if C and D have landed
by 16 November.

**Gate:** IMPROVEMENT_PLAN §4 A/B gates; no wrong answer anywhere. Aim at 600 s: 17 → 23 of 62.

### P4. Speed, parallelism, distribution (fill the gaps)

- **W4 simplex speed** [Claude]: profile first; hyper-sparse FTRAN/BTRAN; faster
  infeasibility proofs. Netlib 0.37 s → 0.25 s.
- **W6 parallel table** [machine, 2 nights]: MIPLIB 600 s at 1, 4 and 8 threads.
- **W7 Windows:** `samaya-setup.exe` per-user installer [Codex sol/medium]; `.mps.gz` input with
  our own inflate and zlib-fixture tests [Codex sol/high]; WDDM measurements [machine]; release
  v0.2.0 after the owner's OK [Claude].

### P5. Finale preparation (from 1 December, or 10 days before the announced date)

1. **Feature freeze** a week before the finale; fixes only.
2. **Final results** on one machine, clock checked, all in `docs/results/` [machine, Claude].
3. **Finale kit:**
   - a 5-minute scripted Studio demo on an MRPL-type case: plan, ranging, verifier panel, the
     GPU toggle, the Excel export [Claude];
   - a 2-minute backup demo from the CLI, and a recorded video of the full demo in case the
     venue's machine or network fails;
   - `cases/demo.sh` with the QP blending case, an SLP pooling case and the large GPU LP
     [Codex luna/medium];
   - the deck and the video updated with the new numbers [Claude]; the robustness report as
     backup.
4. **"Defend the internals" sheets**, one page per component (W8.4): first drafts from the code
   [Codex sol/medium, read-only], checked and finished [Claude]. Every member presents each one
   once before the finale; a mock Q&A with the hardest questions: "why not just use HiGHS", "how
   do you know the answer is right", "why is the GPU slower on small models", "what does MRPL
   need to adopt this".
5. **Sprint backlog:** 8–10 small, visible, testable features for the 36 hours, each with a
   one-paragraph spec written beforehand, so each judging round sees something new. Candidates:
   a solution pool, MIQP cuts, a ranging chart in Studio, a "what changed between two plans"
   view, an MPS export of a Studio case, a REST endpoint around the CLI. Check the finale rules
   on what may be prepared in advance, and build only what they allow.
6. **Change-request drill:** rehearse three likely on-site requests against the clock: "solve
   this MPS file we give you", "add this constraint type to the case", "show this against
   another solver". Know beforehand which files each would touch.
7. **One full rehearsal** of the 36-hour format in the last week of November: 12 hours of work
   with three mock judging rounds.

## 4. How Claude and Codex share the work

Codex is used like a second engineer, not on a quota: any piece of work goes to whichever agent
suits it, decided task by task, and both work at the same time wherever the files allow.

### Who takes what

| Usually Claude | Usually Codex |
|---|---|
| Decisions and specs: method, contract, files, tests, planted bugs, how to verify | Implementation from a spec, including numerics: ranging, clique table, cuts, heuristics, `.mps.gz` inflate |
| Work where a subtle mistake gives wrong answers: presolve aggregator and postsolve, SLP loop, simplex internals | Plumbing and UI: CLI, JSON, Studio, Excel, installer |
| Integration: review of every Codex diff, the Linux debug, release, asan and cuda presets, commits | Python plugins, case generators, test instances, demo scripts, tests where coverage is thin |
| Benchmarks and `docs/results/`; slides, demo, internals sheets (final versions) | Reviews of Claude's branches, design second opinions (astra for the hard ones), research in large files, first drafts |

Tier and effort follow `C:\dev\CLAUDE.md`: luna for mechanical work, sol for most
implementation and reviews, astra for the hardest problems and high-stakes reviews.

### Working rules

- **In parallel:** Codex writes in its own worktree (`-Worktree`) on files Claude is not
  touching (section 5), so both can work at once; several Codex tasks can run side by side on
  disjoint files.
- **Specs** live in `C:\dev\.agents\specs\` and are self-contained: the goal, the files in
  scope, the rules of `CLAUDE.md`, the acceptance criteria, how to verify, and the build recipe
  for Codex's sandbox (`samaya-build-recipe.md`: NMake, Release, Studio off). Codex never
  commits or pushes.
- **When Codex finds a problem in a spec,** it stops and reports; Claude fixes the contract
  rather than letting a test be weakened (as with ranging's dual degeneracy on 13 October).
- **Every Codex result** is reviewed and integrated by Claude, with the four Linux presets and
  the planted bugs, before it is committed.
- **The clock rule:** no builds or agent runs on the laptop while a timed benchmark runs.
- **Limits:** if either agent's 5-hour or weekly limit runs out, the other carries on; nothing
  waits on one agent.

## 5. Who owns which files each week

Codex writes in a worktree while Claude works in the main clone, so the two must not touch the
same files in the same week. The calendar below keeps them apart: Claude's area first, then
Codex's.

| Week | Claude's files | Codex's files |
|---|---|---|
| 1 | `src/linalg/ordering.cpp` merge, `docs/results/` | new `src/lp/ranging.{hpp,cpp}`, `tests/test_ranging.cpp`, new `python/` |
| 2 | `src/presolve/` (aggregator) | `apps/cli`, `apps/studio`, `apps/ui`, `cases/report.py` (Excel) |
| 3 | `src/presolve/`, `docs/results/` | new `src/mip/clique.{hpp,cpp}`, `src/mip/probing.cpp`, their tests; `cases/mrpl.py` |
| 4 | new `src/slp/` | `cases/` (Haverly, blending, QP blending), `apps/studio` + `apps/ui` (tank lines), branching score in `src/mip/branch_and_bound.cpp` |
| 5 | `src/slp/`, `docs/` | `src/mip/cuts.{hpp,cpp}` (zero-half or c-MIR) |
| 6 | `src/lp/simplex.*`, `src/linalg/basis_factor.*` | `src/io/mps_reader.cpp` (`.mps.gz`), `tools/windows/` (installer), `src/mip/heuristics.cpp` |
| 7 | `docs/`, release | `cases/demo.sh`, internals drafts in `docs/sih/` |

The exact paths are checked against `docs/ARCHITECTURE.md` §6 when each spec is written.

## 6. Calendar (IST; shift with the official finale date)

| Week | Dates | Claude | Codex | Machine (nights) |
|---|---|---|---|---|
| 1 | 13–19 Oct | P0 setup, merge `claude/pdlp-logs`; ranging spec; Python plugin spec | `pdlp-logs` review (sol/high); W0 (luna/low); ranging (sol/high); PuLP/Pyomo plugins (sol/medium); tests where thin (luna/medium). | 467k LP, PDLP + crossover; Mittelmann subset |
| 2 | 20–26 Oct | `large_lp.md`; presolve aggregator design and code | aggregator design opinion (astra/high); ranging export + column bounds to CLI/JSON/Studio/Excel (sol/medium); review of ranging integration (sol/high). | barrier and dual simplex on the large LP if needed; SCIP separators |
| 3 | 27 Oct – 2 Nov | aggregator finished, A/B; clique table spec | clique table (sol/high); aggregator review (astra/high); case generator scale-up (sol/medium); luna fixes. | MIPLIB 60/600 s after C |
| 4 | 3–9 Nov | clique A/B and merge; SLP design, then the loop | SLP design opinion (astra/high); Haverly + blending + QP blending cases (sol/medium); Studio tank lines (sol/medium); inference history (sol/medium). | MIPLIB after D |
| 5 | 10–16 Nov | SLP finished; cut spec from the separator results | zero-half or c-MIR (sol/high); SLP review (sol/high); tests where thin (luna). | cases re-run; Maros–Mészáros |
| 6 | 17–23 Nov | W4 profile and hyper-sparse solves; G if allowed | `.mps.gz` (sol/high); installer (sol/medium); shifting/zi-round if allowed (sol/high); simplex review (sol/high). | W6 parallel table |
| 7 | 24–30 Nov | final results; v0.2.0; deck, video; rehearsal | release review (astra/high); `demo.sh` (luna/medium); internals drafts (sol/medium, read-only); luna fixes. | final runs: Netlib, cases, MIPLIB, QP, Mittelmann, large LP |
| 8+ | from 1 Dec | freeze; internals sheets, mock Q&A, demo drills | fixes only (luna/sol). | – |

The Codex column is a starting split; items move between the columns whenever that is faster.

**Critical path:** P1 (the scale run) and P2.1–P2.3 (ranging, Python, SLP). MILP work (P3) is
the largest single improvement but is not needed for a demo; it moves if P1 or P2 runs late.

If the screening results are late, the order stays the same. If the team is not shortlisted,
P2 and P3 are still the useful core for the solver.

## 7. Decisions for the owner

1. **The finale date and the team's roles,** once announced: who presents, who demos, who
   answers internals questions.
2. **Comparison with a commercial solver:** Gurobi and CPLEX have free academic licences, but
   their terms may restrict publishing benchmark results. Read the terms first; otherwise
   position against Mittelmann's public benchmarks through HiGHS.
3. **Code signing** for v0.2.0 (cost), as in POST_SUBMISSION_PLAN.
4. **The four WSL worktrees with uncommitted changes** (`mrpl-pro-act`, `-c3`, `-rens`,
   `-rens2`): keep or drop.
5. **Whether SLP comes before the MILP items C and D** if both cannot fit by mid-November. This
   plan puts the refinery fit first because it is what MRPL's judges will ask about.
