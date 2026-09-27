# Convex QP: Maros–Mészáros (138 models), samaya against HiGHS

## Setup

- **Solvers:** samaya (commit 457101d, interior point on our own sparse LDLᵀ) and HiGHS 1.15.1
  (highspy, its QP solver).
- **Settings:** 300 s per model, one thread each, 4 models at a time.
- **Machine:** the Dell G15 laptop, WSL2. Its clock was checked before the run: factor 1.000.
- **Date:** 27 September 2026, 15:03–15:48 IST.
- **Instances:** the 138 QPS files of the Maros–Mészáros repository (the
  `optimizers/maros-meszaros-mirror` copy), fetched with `bench/fetch_instances.sh maros` and
  run with `bench/compare.sh maros`.
- **Verification:** every samaya answer reported optimal passed `verify_qp_optimality`: the
  reduced costs d = c + Qx − A'y are recomputed from the model in extended precision, with 1e-6
  tolerances.
- **Raw data:** `~/compare-partial-1215/maros/` on the laptop (to be pushed with the other
  laptop results).

**One deviation from the plan:** HiGHS did not return on QFORPLAN within its 300 s limit (it
ran for 24 minutes), and its harness shard was stopped. The shard's 21 finished models are
taken from its log, and the remaining 12 were rerun unchanged (the samaya result for QFORPLAN is
from the first run). QFORPLAN counts as unsolved for HiGHS.

## Results

| | samaya | HiGHS 1.15 |
|---|---|---|
| Solved to optimality | **106 / 138** | 104 / 138 |
| Shifted geometric mean, all 138 (unsolved at 300 s) | **13.1 s** | 16.5 s |
| Shifted geometric mean, the 80 both solve | **0.20 s** | 1.82 s |
| Faster on (of the 80 both solve) | **43** | 37 |
| Solved only by this solver | 26 | 24 |
| Wrong "optimal" answers | **0** | **2** |

- **Where samaya wins:** the families an interior-point method suits: AUG2D/AUG3D (0.03–0.07 s
  vs 18–67 s), CONT-100 to CONT-300, STCQP1 (0.16 s vs 65.7 s), DTOC3, UBH1, the HUESTIS pair,
  MOSARQP1/2.
- **Where HiGHS wins:** the large CVXQP*_L; the LISWET family (samaya's interior point stalls
  on its slow steps); and the degenerate, LP-like Q-versions of Netlib models (QSHIP*,
  QSCFXM*, QSCRS8). On those, samaya's iterate and its polished point both fail the strict
  verifier, and samaya reports `numerical_error` rather than an unverified answer; a crossover
  would finish them.
- **Where HiGHS returns no answer:** `solve_error` on 12 models and no status (`not_set`) on 7.

**Wrong answers (HiGHS).** On two models HiGHS reports "optimal" with a wrong objective. SCIP
(PySCIPOpt, run separately) confirms samaya's values:

| Model | samaya (verified) | SCIP | HiGHS |
|---|---:|---:|---:|
| DPKLO1 | 0.370096217 | 0.370095867 | 0.712522210 |
| QBORE3D | 3100.20080176 | 3100.20080176 | 3102.13879676 |

(SCIP's DPKLO1 value differs from samaya's by 3.5e-7, SCIP's default gap tolerance.)

**Non-convex model.** VALUES is refused by samaya as not convex: its Q has a negative
eigenvalue, −1.27e-5, computed independently with numpy. HiGHS solves it; its answer is not
checked here.

## Reproduce

```sh
bench/fetch_instances.sh maros
SAMAYA_CLOCK_FACTOR=$(python3 bench/clock_factor.py 120) bash bench/compare.sh maros 4
```
