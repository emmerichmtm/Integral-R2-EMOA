# IR2-EMOA EMO 2027 benchmark release

This folder contains the updated paper source and the complete 10-run benchmark used in the current draft.

Benchmark settings: ZDT1/2/3 in 2D; DTLZ1/2/7 in 3D; decision-space dimension 5; population 50; 10,000 evaluations; 10 paired runs per algorithm; algorithms IR2-EMOA, finite R2-EMOA, and SMS-EMOA. Reported performance indicators are exact Integral R2, hypervolume, and Delta_p.

The six diagnostic gallery PNGs are stored in the same directory as `IR2_EMOA_EMO2027.tex`, so the LaTeX source contains no `figures/` path prefixes.

Benchmark scripts and complete run-wise data are under `results/benchmark_10runs_dim5/`.

Current manuscript revision: IR2-EMOA is stated explicitly as a deliberate SMS-EMOA-style $\mu+1\to\mu$ steady-state design. The detailed bi-objective least-contributor derivation/algorithm has been removed and attributed to the R2 v2 paper. A remark explains that exact fixed-cardinality Integral R2 subset selection is efficient in 2D but NP-hard in 3D, citing arXiv:2606.23365 and arXiv:2606.26591. The benchmark table now reports mean +/- sample standard deviation over the ten paired runs; median/IQR summaries remain in `summary.csv`.


## Reproducibility update

The complete benchmark contains 180 run-level results (6 problems x 3 algorithms x 10 paired runs). Exact variation and initialization seeds are saved in `results/benchmark_10runs_dim5/seeds.csv`, and `benchmark_10runs_dim5.py` embeds the same seed values explicitly. The manuscript table reports median +/- sample standard deviation; `summary.csv` also retains means and IQRs.
