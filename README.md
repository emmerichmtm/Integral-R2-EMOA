# IR2-EMOA

Reference implementation and pilot experiments for **IR2-EMOA: An Evolutionary Multiobjective Optimization Algorithm Based on the Integral R2 Indicator**.

The repository is intentionally focused on two and three objectives and on the steady-state indicator-based comparison with SMS-EMOA and finite-weight R2-EMOA.

## Core convention

For minimization, losses are measured from an ideal point `z*` and may be **zero**. Under perspective mapping,

```text
loss_j = 0  ->  reciprocal_j = +infinity
```

The code uses an actual IEEE `+inf`; it never substitutes epsilon or a large finite number. Degenerate cells (`lo == hi`, including `[inf, inf]`) are assigned zero before the 4-/8-corner antiderivative is evaluated. For 3D, nondegenerate reciprocal cells are integrated with the Jacobian density `(x+y+z)^(-4)`.

## Algorithms

- `ir2`: proposed Integral R2 EMOA, steady-state, nondominated sorting, least Integral R2 deletion loss in the worst front.
- `r2`: finite-weight R2-EMOA in the same steady-state architecture, using a simplex lattice of weight vectors.
- `sms`: SMS-style steady-state deletion by hypervolume contribution. For the pilot code, the worst front is affinely normalized and a dystopian reference `(1.1,...,1.1)` is used.

Variation is SBX (`pc=0.9`, `eta_c=15`) plus polynomial mutation (`pm=1/n`, `eta_m=20`), matching the classic R2-EMOA settings.

## Pilot results

The checked-in pilot uses DTLZ1 and DTLZ2 in 2D and 3D, 3 common-start runs, population 30, and 500 evaluations. This is a software/regression study, **not yet the final EMO 2027 experiment**. Publication runs should use the frozen larger budget and number of repetitions.

Primary quality measures are `Delta_p` (averaged Hausdorff distance, p=2) and mean convergence distance to a dense known Pareto-front reference set. Hypervolume is not used as a performance metric.

Reproduce:

```bash
python -m unittest discover -s tests -v
python run_experiments.py --runs 3 --evaluations 500
python analyze_results.py
```

## Files

- `integral_r2.py` — exact 2D/3D deletion losses, including ideal-boundary `+inf` handling.
- `ir2_emoa.py` — IR2-EMOA, finite-weight R2-EMOA comparator, and SMS-style comparator.
- `benchmarks.py` — DTLZ1/2/4.
- `metrics.py` — Delta_p and convergence.
- `run_experiments.py` — reproducible pilot experiment runner.
- `analyze_results.py` — median/IQR summary.
- `results/` — raw and summarized pilot results.

## Status

The 3D Integral R2 routine here is the transparent exact O(n^2) z-slab/x-sweep implementation with correct `+inf` boundary semantics. The separate contribution-calculator repository contains the research fast sweep. The final paper repository should vendor or pin the same verified fast implementation once the degenerate/infinite-coordinate rules are merged there.
