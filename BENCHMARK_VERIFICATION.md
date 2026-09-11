# Benchmark verification

The benchmark result set in `results/benchmark_10runs_dim5/` was checked before packaging.

- 6 problems: ZDT1, ZDT2, ZDT3, DTLZ1, DTLZ2, DTLZ7.
- 3 algorithms per problem: IR2-EMOA, finite R2-EMOA, SMS-EMOA.
- 10 paired runs per problem/algorithm.
- 180 run-level records in `raw.csv` and `raw.json`.
- 10,000 evaluations, population 50, decision dimension 5.
- Every recorded final approximation has 50 nondominated points.
- Within each problem/run, all three algorithms use the same variation seed and common initial population.
- Exact variation seeds and initialization seeds are saved in `seeds.csv`.
- Initialization seed = variation seed + 17.
- Table reporting is median +/- sample standard deviation across the ten runs.
- `summary.csv` retains median, IQR, mean, and sample standard deviation for audit/reanalysis.

The benchmark script now embeds the exact saved seeds explicitly, so rerunning the script does not depend on a hidden seed-generation convention.
