# Netlib results

Generated with `bench/fetch_instances.sh netlib netlib-infeas` and
`bench/harness.py <set> --baseline highspy --time-limit 120` on a 4-core Intel Xeon
(2.8 GHz), release build, one thread. HiGHS 1.15.1 (highspy) is the comparison baseline
only; it is not linked into samaya. Objectives agree when the relative difference is
below 1e-6. Every samaya result below was checked by the independent verifier.

## Feasible set (93 instances)

samaya solved 93/93, shifted geometric mean 0.315 s (HiGHS 93/93, 0.167 s; shift 10 s, as in the Mittelmann benchmarks).

| instance | samaya status | samaya objective | samaya s | HiGHS status | HiGHS objective | HiGHS s |
|---|---|---:|---:|---|---:|---:|
| 25fv47 | optimal | 5501.845888 | 0.181 | optimal | 5501.845888 | 0.157 |
| 80bau3b | optimal | 987224.1924 | 0.239 | optimal | 987224.1924 | 0.127 |
| adlittle | optimal | 225494.9632 | 0.001 | optimal | 225494.9632 | 0.001 |
| afiro | optimal | -464.7531429 | 0.000 | optimal | -464.7531429 | 0.001 |
| agg | optimal | -35991767.29 | 0.004 | optimal | -35991767.29 | 0.004 |
| agg2 | optimal | -20239252.36 | 0.004 | optimal | -20239252.36 | 0.006 |
| agg3 | optimal | 10312115.94 | 0.005 | optimal | 10312115.94 | 0.006 |
| bandm | optimal | -158.6280185 | 0.012 | optimal | -158.6280185 | 0.009 |
| beaconfd | optimal | 33592.48581 | 0.002 | optimal | 33592.48581 | 0.002 |
| blend | optimal | -30.81214985 | 0.001 | optimal | -30.81214985 | 0.002 |
| bnl1 | optimal | 1977.629562 | 0.034 | optimal | 1977.629562 | 0.025 |
| bnl2 | optimal | 1811.23654 | 0.123 | optimal | 1811.23654 | 0.053 |
| boeing1 | optimal | -335.2135675 | 0.012 | optimal | -335.2135675 | 0.011 |
| boeing2 | optimal | -315.018728 | 0.002 | optimal | -315.018728 | 0.003 |
| bore3d | optimal | 1373.080394 | 0.003 | optimal | 1373.080394 | 0.002 |
| brandy | optimal | 1518.509896 | 0.005 | optimal | 1518.509896 | 0.005 |
| capri | optimal | 2690.012914 | 0.004 | optimal | 2690.012914 | 0.004 |
| cycle | optimal | -5.226393025 | 0.132 | optimal | -5.226393025 | 0.249 |
| czprob | optimal | 2185196.699 | 0.049 | optimal | 2185196.699 | 0.027 |
| d2q06c | optimal | 122784.2108 | 1.234 | optimal | 122784.2108 | 0.835 |
| d6cube | optimal | 315.4916667 | 0.059 | optimal | 315.4916667 | 0.127 |
| degen2 | optimal | -1435.178 | 0.018 | optimal | -1435.178 | 0.014 |
| degen3 | optimal | -987.294 | 0.246 | optimal | -987.294 | 0.198 |
| dfl001 | optimal | 11266396.05 | 15.240 | optimal | 11266396.05 | 6.974 |
| e226 | optimal | -11.63892907 | 0.007 | optimal | -11.63892907 | 0.007 |
| etamacro | optimal | -755.7152333 | 0.012 | optimal | -755.7152333 | 0.009 |
| fffff800 | optimal | 555679.5648 | 0.024 | optimal | 555679.5648 | 0.011 |
| finnis | optimal | 172791.0656 | 0.006 | optimal | 172791.0656 | 0.005 |
| fit1d | optimal | -9146.378092 | 0.006 | optimal | -9146.378092 | 0.011 |
| fit1p | optimal | 9146.378092 | 0.046 | optimal | 9146.378092 | 0.040 |
| fit2d | optimal | -68464.29329 | 0.252 | optimal | -68464.29329 | 0.159 |
| fit2p | optimal | 68464.29329 | 1.554 | optimal | 68464.29329 | 1.079 |
| forplan | optimal | -664.2189613 | 0.003 | optimal | -664.2189613 | 0.010 |
| ganges | optimal | -109585.7361 | 0.035 | optimal | -109585.7361 | 0.013 |
| gfrd-pnc | optimal | 6902236 | 0.010 | optimal | 6902236 | 0.006 |
| greenbea | optimal | -72555248.13 | 1.158 | optimal | -72555248.13 | 0.238 |
| greenbeb | optimal | -4302260.261 | 1.381 | optimal | -4302260.261 | 0.431 |
| grow15 | optimal | -106870941.3 | 0.045 | optimal | -106870941.3 | 0.040 |
| grow22 | optimal | -160834336.5 | 0.103 | optimal | -160834336.5 | 0.083 |
| grow7 | optimal | -47787811.81 | 0.010 | optimal | -47787811.81 | 0.012 |
| israel | optimal | -896644.8219 | 0.002 | optimal | -896644.8219 | 0.003 |
| kb2 | optimal | -1749.90013 | 0.000 | optimal | -1749.90013 | 0.001 |
| lotfi | optimal | -25.26470606 | 0.003 | optimal | -25.26470606 | 0.002 |
| maros | optimal | -58063.7437 | 0.087 | optimal | -58063.7437 | 0.041 |
| maros-r7 | optimal | 1497185.166 | 5.792 | optimal | 1497185.166 | 0.771 |
| modszk1 | optimal | 320.6197291 | 0.022 | optimal | 320.6197291 | 0.014 |
| nesm | optimal | 14076036.49 | 0.184 | optimal | 14076036.49 | 0.150 |
| perold | optimal | -9380.755278 | 0.063 | optimal | -9380.755278 | 0.057 |
| pilot | optimal | -557.4897062 | 1.450 | optimal | -557.4897293 | 1.081 |
| pilot.ja | optimal | -6113.136466 | 0.117 | optimal | -6113.136466 | 0.079 |
| pilot.we | optimal | -2720107.533 | 0.189 | optimal | -2720107.533 | 0.139 |
| pilot4 | optimal | -2581.139259 | 0.034 | optimal | -2581.139259 | 0.028 |
| pilot87 | optimal | 301.7103475 | 7.494 | optimal | 301.7103473 | 3.872 |
| pilotnov | optimal | -4497.276188 | 0.059 | optimal | -4497.276188 | 0.151 |
| recipe | optimal | -266.616 | 0.000 | optimal | -266.616 | 0.001 |
| sc105 | optimal | -52.20206121 | 0.001 | optimal | -52.20206121 | 0.001 |
| sc205 | optimal | -52.20206121 | 0.003 | optimal | -52.20206121 | 0.002 |
| sc50a | optimal | -64.57507706 | 0.000 | optimal | -64.57507706 | 0.001 |
| sc50b | optimal | -70 | 0.000 | optimal | -70 | 0.001 |
| scagr25 | optimal | -14753433.06 | 0.009 | optimal | -14753433.06 | 0.006 |
| scagr7 | optimal | -2331389.824 | 0.001 | optimal | -2331389.824 | 0.002 |
| scfxm1 | optimal | 18416.75903 | 0.007 | optimal | 18416.75903 | 0.008 |
| scfxm2 | optimal | 36660.26156 | 0.022 | optimal | 36660.26156 | 0.018 |
| scfxm3 | optimal | 54901.25455 | 0.042 | optimal | 54901.25455 | 0.029 |
| scorpion | optimal | 1878.124823 | 0.004 | optimal | 1878.124823 | 0.003 |
| scrs8 | optimal | 904.2969538 | 0.011 | optimal | 904.2969538 | 0.008 |
| scsd1 | optimal | 8.666666674 | 0.002 | optimal | 8.666666674 | 0.002 |
| scsd6 | optimal | 50.50000008 | 0.005 | optimal | 50.50000008 | 0.007 |
| scsd8 | optimal | 904.9999999 | 0.024 | optimal | 904.9999999 | 0.038 |
| sctap1 | optimal | 1412.25 | 0.004 | optimal | 1412.25 | 0.005 |
| sctap2 | optimal | 1724.807143 | 0.016 | optimal | 1724.807143 | 0.012 |
| sctap3 | optimal | 1424 | 0.026 | optimal | 1424 | 0.015 |
| seba | optimal | 15711.6 | 0.006 | optimal | 15711.6 | 0.004 |
| share1b | optimal | -76589.31858 | 0.002 | optimal | -76589.31858 | 0.003 |
| share2b | optimal | -415.7322407 | 0.001 | optimal | -415.7322407 | 0.002 |
| shell | optimal | 1208825346 | 0.010 | optimal | 1208825346 | 0.012 |
| ship04l | optimal | 1793324.538 | 0.005 | optimal | 1793324.538 | 0.007 |
| ship04s | optimal | 1798714.7 | 0.004 | optimal | 1798714.7 | 0.005 |
| ship08l | optimal | 1909055.211 | 0.015 | optimal | 1909055.211 | 0.017 |
| ship08s | optimal | 1920098.211 | 0.011 | optimal | 1920098.211 | 0.009 |
| ship12l | optimal | 1470187.919 | 0.030 | optimal | 1470187.919 | 0.024 |
| ship12s | optimal | 1489236.134 | 0.025 | optimal | 1489236.134 | 0.011 |
| sierra | optimal | 15394362.18 | 0.015 | optimal | 15394362.18 | 0.014 |
| stair | optimal | -251.2669512 | 0.016 | optimal | -251.2669512 | 0.013 |
| standata | optimal | 1257.6995 | 0.001 | optimal | 1257.6995 | 0.004 |
| standgub | optimal | 1257.6995 | 0.001 | optimal | 1257.6995 | 0.004 |
| standmps | optimal | 1406.0175 | 0.003 | optimal | 1406.0175 | 0.005 |
| stocfor1 | optimal | -41131.97622 | 0.001 | optimal | -41131.97622 | 0.001 |
| stocfor2 | optimal | -39024.40854 | 0.141 | optimal | -39024.40854 | 0.025 |
| tuff | optimal | 0.2921477651 | 0.005 | optimal | 0.2921477651 | 0.007 |
| vtp.base | optimal | 129831.4625 | 0.002 | optimal | 129831.4625 | 0.001 |
| wood1p | optimal | 1.442902412 | 0.018 | optimal | 1.442902412 | 0.068 |
| woodw | optimal | 1.304476333 | 0.120 | optimal | 1.304476333 | 0.064 |

## Barrier and crossover (laptop, 7 October 2026)

`bench/compare.sh netlib-barrier`, `netlib-crossover` and `netlib` on the Dell G15 laptop
(Intel Core i5-13450HX, WSL2), alone on the machine, clock checked (factor 1.000), 4 instances at
a time, one thread each, build b5cc218 (branch `claude/crossover`). Every answer below passed
the verifier, and every objective equals HiGHS's to 1e-6.

| Method | Solved | Shifted geomean |
|---|---|---|
| samaya, dual simplex (default) | 93/93 | 0.37 s |
| samaya, `--lp-method barrier` | 93/93 | 0.75 s |
| samaya, `--lp-method barrier --crossover` | 93/93 | 0.77 s |
| HiGHS 1.15.1 | 93/93 | 0.21 s |
| CBC 2.10.11 (Clp) | 93/93 | 0.23 s |
| GLPK 5.0 | 93/93 | 0.45 s |
| SCIP 10.0 (SoPlex) | 93/93 | 0.50 s |

- With `--crossover`, 73 of the 93 answers are basic solutions from the crossover basis
  (median 110 simplex iterations after it). The crossover costs about 0.02 s of geomean over
  the barrier alone.
- The other 20 are the models where the barrier stops at its iteration limit (80bau3b, agg,
  agg2, agg3, bnl2, d2q06c, d6cube, dfl001, etamacro, finnis, greenbea, greenbeb, modszk1, nesm,
  pilot, pilot.ja, pilot.we, pilot87, recipe, share1b); the dual simplex solves them from the
  slack basis. Starting it from the crossover basis of the unconverged iterate was tried and made
  most of them slower (pilot87 11,109 to 115,446 iterations; dfl001 failed), so it is not done.
- On Netlib the barrier is about twice as slow as the dual simplex; it is there for large and
  dense models, and the crossover makes its answers basic.

## Infeasible set (29 instances)

samaya solved 28/29, shifted geometric mean 0.103 s (HiGHS 29/29, 0.008 s; shift 10 s, as in the Mittelmann benchmarks).

| instance | samaya status | samaya objective | samaya s | HiGHS status | HiGHS objective | HiGHS s |
|---|---|---:|---:|---|---:|---:|
| bgdbg1 | infeasible |  | 0.000 | infeasible |  | 0.000 |
| bgetam | infeasible |  | 0.001 | infeasible |  | 0.003 |
| bgindy | infeasible |  | 0.006 | infeasible |  | 0.012 |
| bgprtr | infeasible |  | 0.000 | infeasible |  | 0.000 |
| box1 | infeasible |  | 0.000 | infeasible |  | 0.001 |
| ceria3d | infeasible |  | 0.483 | infeasible |  | 0.006 |
| chemcom | infeasible |  | 0.000 | infeasible |  | 0.002 |
| cplex1 | infeasible |  | 0.088 | infeasible |  | 0.118 |
| cplex2 | numerical_error |  | 0.014 | infeasible |  | 0.005 |
| ex72a | infeasible |  | 0.001 | infeasible |  | 0.001 |
| ex73a | infeasible |  | 0.001 | infeasible |  | 0.001 |
| forest6 | infeasible |  | 0.000 | infeasible |  | 0.001 |
| galenet | infeasible |  | 0.000 | infeasible |  | 0.000 |
| gosh | infeasible |  | 2.339 | infeasible |  | 0.021 |
| gran | infeasible |  | 0.080 | infeasible |  | 0.003 |
| greenbea | infeasible |  | 0.109 | infeasible |  | 0.005 |
| itest2 | infeasible |  | 0.000 | infeasible |  | 0.000 |
| itest6 | infeasible |  | 0.000 | infeasible |  | 0.000 |
| klein1 | infeasible |  | 0.001 | infeasible |  | 0.002 |
| klein2 | infeasible |  | 0.007 | infeasible |  | 0.008 |
| klein3 | infeasible |  | 0.031 | infeasible |  | 0.026 |
| mondou2 | infeasible |  | 0.002 | infeasible |  | 0.001 |
| pang | infeasible |  | 0.007 | infeasible |  | 0.007 |
| pilot4i | infeasible |  | 0.021 | infeasible |  | 0.001 |
| qual | infeasible |  | 0.007 | infeasible |  | 0.006 |
| reactor | infeasible |  | 0.001 | infeasible |  | 0.001 |
| refinery | infeasible |  | 0.005 | infeasible |  | 0.004 |
| vol1 | infeasible |  | 0.007 | infeasible |  | 0.009 |
| woodinfe | infeasible |  | 0.000 | infeasible |  | 0.000 |

## Notes

* `cplex2` is infeasible by less than the feasibility tolerance: the Farkas certificate
  samaya finds proves infeasibility only by about 1e-9, below the 1e-7 tolerance, so the
  verifier does not accept it and samaya reports `numerical_error` instead of claiming a
  result it cannot prove. HiGHS reports it as infeasible.
* `qap8`, `qap12`, `qap15` (generated, not stored as data on netlib) and `truss` (Fortran
  source) are not part of the downloadable set.
