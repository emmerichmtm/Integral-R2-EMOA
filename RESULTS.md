# Pilot IR2-EMOA results

These results are **software-verification pilots**, not the final EMO 2027 experiment.

Settings: DTLZ1 and DTLZ2; 2 and 3 objectives; population 30; 500 function evaluations; 3 common-start runs; SBX + polynomial mutation. Primary measures are Delta_p (p=2) and mean convergence distance to a dense known Pareto-front reference set. Hypervolume is not used as an evaluation measure.

| Problem | m | Algorithm | Median Delta_p | IQR Delta_p | Median convergence | IQR convergence | Median runtime (s) |
|---|---:|---|---:|---:|---:|---:|---:|
| DTLZ1 | 2 | IR2-EMOA | 19.7174 | 16.0400 | 18.3006 | 13.1053 | 2.199 |
| DTLZ1 | 2 | finite R2-EMOA | 37.3407 | 8.6445 | 23.9890 | 7.6579 | 2.288 |
| DTLZ1 | 2 | SMS-EMOA | 22.4962 | 15.1913 | 21.3223 | 5.4395 | 2.195 |
| DTLZ1 | 3 | IR2-EMOA | 56.8925 | 27.0299 | 46.5300 | 18.1391 | 2.158 |
| DTLZ1 | 3 | finite R2-EMOA | 42.3898 | 3.2421 | 33.7213 | 4.0929 | 2.176 |
| DTLZ1 | 3 | SMS-EMOA | 81.1614 | 20.9246 | 61.5049 | 20.7782 | 2.100 |
| DTLZ2 | 2 | IR2-EMOA | 0.09060 | 0.02288 | 0.04754 | 0.01261 | 2.085 |
| DTLZ2 | 2 | finite R2-EMOA | 0.08811 | 0.02098 | 0.04602 | 0.02336 | 2.206 |
| DTLZ2 | 2 | SMS-EMOA | 0.10049 | 0.03512 | 0.04605 | 0.00724 | 2.111 |
| DTLZ2 | 3 | IR2-EMOA | 0.22760 | 0.02928 | 0.15419 | 0.01072 | 2.227 |
| DTLZ2 | 3 | finite R2-EMOA | 0.21887 | 0.02043 | 0.13227 | 0.03892 | 2.248 |
| DTLZ2 | 3 | SMS-EMOA | 0.18908 | 0.01000 | 0.12345 | 0.01363 | 2.168 |

## Interpretation

The pilot is only large enough to validate code paths and data handling. DTLZ2 shows comparable 2D behavior among the three methods and a better pilot Delta_p for SMS-EMOA in 3D. DTLZ1 remains far from convergence and highly variable at 500 evaluations; no substantive algorithmic conclusion should be drawn from those rows.

The main experimental campaign should therefore use a much larger evaluation budget and at least the planned 30 independent runs before statistical testing.

## Exact-ideal 2D single shots

The following deterministic single shots were rerun after removing the old strict-ideal requirement.  Settings are population 100, 10,000 evaluations, and seed 2027.  The exact ideal point is passed directly to IR2-EMOA.

| Problem | Exact ideal point | Delta_p | Convergence | Final nondominated size |
|---|---|---:|---:|---:|
| ZDT1 | (0, 0) | 0.008408 | 0.006272 | 100 |
| ZDT2 | (0, 0) | 0.011097 | 0.007000 | 100 |
| ZDT3 | (0, -0.7733690123) | 0.098782 | 0.002498 | 100 |
| DTLZ1 | (0, 0) | 0.022007 | 0.021915 | 100 |
| DTLZ2 | (0, 0) | 0.009135 | 0.000072 | 100 |
| DTLZ7 | (0, 2.3070043655) | 0.007219 | 0.001557 | 100 |

These are illustrative single runs, not inferential results.  For ZDT3 and DTLZ7 the exact ideal point is componentwise and is not itself generally feasible.  The PNGs are in `figures/`; raw final objective vectors and metadata are in `results/single_shot_exact_ideal/`.
