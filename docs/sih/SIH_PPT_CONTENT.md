# SIH 2026 presentation: content for the standard template

Problem statement 26119 (MRPL). Team solution: **samaya**.

This file is for the person building the slides. It follows the standard SIH idea-presentation
template: 6 slides, in the order and with the headings below.

**Before you start:**
- **Download the official 2026 template** from the SIH portal (sih.gov.in) and put this
  content into it.
- If the 2026 template differs from the 6 slides here, keep the official one and move the
  content to the matching headings.
- Do not change the title-slide layout or remove the SIH logos.
- Submit as PDF if the portal asks for PDF.

**Rules for the content:**
- Keep every number exactly as written; each comes from a measured run in the repository
  (`docs/results/`).
- Items marked **[FILL]** need the team's details or final results.
- Short bullets, big numbers, one diagram per slide. The judges read fast.

---

## Slide 1: Title page (template fields)

| Field | Content |
|---|---|
| Problem Statement ID | 26119 |
| Problem Statement Title | **[FILL: copy exactly from the SIH portal]**. Ours is about an indigenous, GPU-accelerated optimization (LP/MILP) solver for MRPL. |
| Theme | **[FILL: from the portal]** |
| PS Category | Software |
| Team ID | **[FILL]** |
| Team Name | **[FILL]** |

Idea title to use on slide 2: **"samaya: an indigenous, verified, GPU-accelerated
optimization solver for refinery planning and scheduling"**.
(samaya = Sovereign Accelerated Mathematical Algorithms for Yield & Allocation.)

---

## Slide 2: Idea title and proposed solution

**Heading:** samaya: our own optimization solver, built from scratch

**Proposed solution (bullets for the slide):**
- Refineries decide every day which crude to buy, how to run units, how to blend products and
  when to receive ships. These are **linear and mixed-integer optimization problems (LP/MILP)**.
- Today they depend on foreign solvers (CPLEX, Gurobi) or open-source ones built abroad.
- We built **samaya, a complete LP/MILP solver written from scratch in C++20**. No external
  solver or linear-algebra code is inside; HiGHS, SCIP, CBC and GLPK are used only to compare
  results.
- It reads standard model files (MPS), so existing refinery models run unchanged. It has a
  command-line tool and C/C++ APIs.

**How it addresses the problem:**
- It solves the three MRPL-type decisions we modelled:
  - refinery planning (LP);
  - crude receipt scheduling (MILP);
  - captive power and steam unit commitment (MILP).
- **Daily re-planning:** it starts from yesterday's plan and keeps it if it is still feasible,
  so the plan does not change without reason.

**Innovation and uniqueness:**
- **Every answer is independently verified.** A separate checker recomputes each result from
  the original model with mathematical certificates (optimality, infeasibility proofs) before
  it is reported. It is on by default.
- **GPU acceleration** for very large LPs with PDLP, the first-order method behind Google's PDLP
  and NVIDIA's cuOpt. **[FILL: GPU teammate's one-line result]**
- **Modern MILP techniques implemented by us:**
  - Feasibility Jump (2023 research);
  - cutting planes inside the search tree;
  - multi-threaded tree search.

**Visual:** three refinery icons (planning, scheduling, utilities) → "samaya" box → "verified
plan" with a check mark.

---

## Slide 3: Technical approach

**Technologies:**
- **Language and build:** C++20, CMake; runs on Linux and Windows (WSL).
- **GPU:** CUDA on an NVIDIA RTX 3060, for the PDLP LP solver.
- **Benchmarks and case generators:** Python.
- **No external solver libraries.**

**Methodology (draw as a flow chart):**
```
Refinery model (MPS file)
      │
      ▼
Presolve: shrinks the model
      │
      ├──► LP:   dual simplex with our own sparse LU factorization
      │          (or GPU PDLP for very large LPs)
      │
      └──► MILP: branch-and-cut
                 • cutting planes (Gomory, MIR, covers, flow covers), also in the tree
                 • heuristics (Feasibility Jump, feasibility pump, diving, RINS/RENS)
                 • reliability branching, restarts, parallel search
      │
      ▼
Independent verifier (certificates) ──► Result / plan report
```

**Built and working today (for the speaker):**
- An LP engine: dual simplex, sparse LU with updates, presolve, scaling.
- A MILP engine: branch-and-cut with the techniques above.
- The verifier.
- The case studies with a readable plan report.
- The benchmark harness against HiGHS, SCIP, CBC and GLPK.

**Visual:** the flow chart above, as boxes and arrows. Optionally, a screenshot of a crude
schedule from the plan report (`cases/mrpl.py report`).

---

## Slide 4: Feasibility and viability

**Feasibility, already shown by working results** (one laptop, Intel Core i5-13450HX, same
settings for all solvers):

| Test set | samaya | HiGHS | SCIP | CBC | GLPK |
|---|---|---|---|---|---|
| Netlib LP, 93 models: solved | **93/93** | 93/93 | 93/93 | 93/93 | 93/93 |
| Netlib LP: time (geomean) | 0.36 s | 0.19 s | 0.47 s | 0.20 s | 0.42 s |
| Refinery cases + generated MILPs, 16 models: solved | **16/16** | 16/16 | 16/16 | 16/16 | 11/16 |
| Refinery cases: time (geomean) | **3.43 s (fastest)** | 5.18 s | 3.59 s | 4.02 s | 20.8 s |
| MIPLIB 2017, 62 hard MILPs, 60 s: solved | 13 (was 7 before 26 Sep) | 18 | 19 | 8 | – |

Source: `docs/results/comparison.md`, laptop section (it also has the earlier cloud-machine run,
where samaya was third on the refinery cases). The MIPLIB run gave every solver the same real
time (a laptop clock quirk, corrected and explained there). samaya solves 4 MIPLIB instances
that neither HiGHS nor SCIP solves in 60 s. **[FILL]** Add the 600 s MIPLIB row if it arrives
in time.

- **Correctness:** no wrong answer on any instance; every result is cross-checked against the
  other solvers and published optimal values.
- On Netlib, **28 of 29 infeasible models are proven infeasible with a certificate**.
- samaya solves **4 MIPLIB instances that neither HiGHS nor SCIP solves in 60 s**
  (markshare_4_0, neos5, mas76, pk1).

**Challenges and risks → how we handle them:**

| Challenge | Strategy |
|---|---|
| Hard MILPs: mature solvers have 10–15 years of techniques | We add the proven techniques one at a time and measure each on fixed test sets. Since 26 September: conflict analysis, cut-set flow covers on network models (p200x1188c from unsolved to 3 s), clique cuts; MIPLIB solved at 60 s went from 7 to 13. |
| Numerical errors (wrong answers) | Independent verifier; automatic re-solve with tighter tolerances; tests against a reference solver; deliberately planted bugs must be caught by the tests. |
| Very large LPs | GPU PDLP (first-order method), with the simplex polishing the result. |
| Adoption | Standard MPS input, a CLI and APIs; it runs on ordinary hardware. |

**Viability:**
- It runs on a laptop and has no licence cost.
- The code is modular: QP (blending and pooling) and larger models fit into the same
  structure.

---

## Slide 5: Impact and benefits

**Impact on MRPL and the refining sector:**
- **A sovereign solver core:** no dependence on foreign licences for planning and scheduling
  decisions.
- **Trust:** every plan comes with a verified check that it meets all constraints (tanks,
  units, specifications, demand).
- **Daily re-planning** from yesterday's plan, with stable plans that don't change without
  reason.

**Benefits:**
- **Economic:**
  - No per-seat solver licences.
  - Better crude slates and schedules. Example from our case study: the LP reports the value of
    each unit's capacity; one more kt a week of reformer capacity is worth about 1,800 lakh
    rupees in the small planning case, which informs debottlenecking decisions.
- **Strategic:** Atmanirbhar Bharat, an indigenous core technology for the energy sector,
  reusable by other PSUs (refineries, power, logistics).
- **Environmental:** optimized unit commitment and fuel use in the captive power and steam
  case, with less fuel wasted.
- **Academic:** an open, documented codebase for teaching and research in optimization.

**Visual:** 3–4 icons with one line each (licence-free, verified, faster decisions, greener
operation).

---

## Slide 6: Research and references

1. **Dual simplex:**
   - Koberstein, A., *The dual simplex method, techniques for a fast and stable
     implementation*, PhD thesis, 2005.
   - Forrest, J. and Goldfarb, D., "Steepest-edge simplex algorithms for linear programming",
     Math. Programming, 1992.
2. **MILP:** Achterberg, T., *Constraint Integer Programming*, PhD thesis, TU Berlin, 2007
   (SCIP's design).
3. **Cutting planes:** Marchand, H. and Wolsey, L., "Aggregated mixed-integer rounding",
   Operations Research, 2001.
4. **Feasibility Jump:** Luteberget, B. and Sartor, G., "Feasibility Jump: an LP-free
   Lagrangian MIP heuristic", Math. Programming Computation, 2023.
5. **Conflict analysis:** Witzig, J., Berthold, T. and Heinz, S., "Experiments with conflict
   analysis in mixed integer programming", CPAIOR, 2017.
6. **PDLP:** Applegate, D. et al., "Practical large-scale linear programming using
   primal-dual hybrid gradient", NeurIPS, 2021.
7. **Benchmarks:** MIPLIB 2017 (miplib.zib.de); the Netlib LP collection.
8. **Compared solvers:** HiGHS (highs.dev), SCIP 10 (scipopt.org), CBC (COIN-OR), GLPK.
9. **Project repository:** **[FILL: GitHub link, if the team shares it]**

---

## Speaker notes: questions the jury may ask

- **"Is it really from scratch?"** Yes. No external solver or factorization library is inside.
  HiGHS, SCIP, CBC and GLPK are only used in benchmark scripts for comparison.
- **"Why is it slower than HiGHS on MIPLIB?"** MILP speed comes from decades of techniques. We
  built the core ones in weeks, and each new one closes part of the gap on fixed test sets.
  On refinery-type models we are already competitive.
- **"How do you know the answers are right?"** Through the independent verifier with
  certificates, plus cross-checks against four other solvers and published optima: zero
  disagreements.
- **"Did you choose tests that make you look good?"** No. The test sets were fixed before
  running, and every result is reported, including the ones where we lose.
