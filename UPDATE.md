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
