# Git update notes

This folder is ready to become the IR2-EMOA research repository or to be copied into an existing one.

## Main changes

1. Added exact 2D/3D Integral R2 deletion losses with the corrected ideal-boundary convention:
   - loss `0` maps to actual `+inf`;
   - degenerate cells, including `[inf, inf]`, return zero before corner evaluation;
   - 3D cells use the Jacobian density `(x+y+z)^(-4)` and the eight-corner antiderivative.
2. Added the steady-state IR2-EMOA implementation.
3. Added finite-weight R2-EMOA and SMS-style controlled comparators.
4. Added DTLZ1/2/4, Delta_p, convergence, and reproducible experiment scripts.
5. Added checked-in pilot raw/summary CSV results and a manuscript source with the pilot table.

## Suggested commit

```bash
git add .
git commit -m "Add IR2-EMOA pilot, exact infinity boundary handling, and reproducible results"
git push
```

The checked-in results are explicitly marked as pilot/software-verification results. Do not present them as the final EMO 2027 benchmark.

## Exact-ideal single-shot rerun

- Added ZDT1, ZDT2, ZDT3 and DTLZ7 benchmark definitions.
- `ideal()` now returns the exact componentwise benchmark ideal; no artificial strict-dominance offset is used.
- ZDT3 uses `z*_2=-0.7733690123266405`; DTLZ7 (2D) uses `z*_2=2.3070043655015775`.
- Repeated the six 2D IR2-EMOA single shots for ZDT1/2/3 and DTLZ1/2/7 with population 100, 10,000 evaluations, seed 2027.
- Added the six PNGs and raw final objective vectors.
- Vectorized nondominated sorting and reused surviving ranks after deletion; this changes runtime only, not the steady-state selection rule.

- Replotted DTLZ1 3D with the true Pareto background drawn first so the approximation points stay in the foreground.
- Increased the DTLZ7 3D single-shot budget to 60,000 evaluations and updated the recorded figure and metadata accordingly.

- Added a six-plot gallery to the LaTeX manuscript for the current exact-ideal single-shot figures (ZDT1, ZDT2, ZDT3, DTLZ1, DTLZ2, DTLZ7).
