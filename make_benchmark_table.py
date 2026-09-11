"""Generate the compact LaTeX benchmark table as median +/- sample SD."""
from pathlib import Path
import csv

IN = Path('results/benchmark_10runs_dim5/summary.csv')
OUT = Path('results/benchmark_10runs_dim5/benchmark_table.tex')
ORDER = ['ir2', 'r2', 'sms']
LABEL = {'ir2': 'IR2-EMOA', 'r2': 'finite R2-EMOA', 'sms': 'SMS-EMOA'}
PROBS = ['ZDT1', 'ZDT2', 'ZDT3', 'DTLZ1', 'DTLZ2', 'DTLZ7']


def fmt(v):
    x=float(v)
    if abs(x) < 1e-4 and x != 0.0:
        return f'{x:.2e}'
    return f'{x:.6f}'


def pm(median, sd):
    return f'{fmt(median)}\\pm{fmt(sd)}'


rows = list(csv.DictReader(IN.open()))
by = {(r['problem'], r['algorithm']): r for r in rows}
lines = [
    r'\begin{table}[t]',
    r'\centering',
    r'\caption{Ten-run benchmark results as median $\pm$ sample standard deviation with five decision variables, population size $50$, and $10{,}000$ evaluations. Smaller is better for IR2 and $\Delta_p$; larger is better for HV. Exact paired seeds and all run-wise values are provided in the repository.}',
    r'\label{tab:benchmark10}',
    r'\scriptsize',
    r'\setlength{\tabcolsep}{2.5pt}',
    r'\begin{tabular}{llccc}',
    r'\toprule',
    r'Problem & Algorithm & IR2 $\downarrow$ & HV $\uparrow$ & $\Delta_p\downarrow$\\',
    r'\midrule',
]
for pi, p in enumerate(PROBS):
    S = [by[(p, a)] for a in ORDER]
    best_ir2 = min(range(3), key=lambda i: float(S[i]['integral_r2_median']))
    best_hv = max(range(3), key=lambda i: float(S[i]['hypervolume_median']))
    best_dp = min(range(3), key=lambda i: float(S[i]['delta_p_median']))
    for i, a in enumerate(ORDER):
        r = by[(p, a)]
        vals = [
            pm(r['integral_r2_median'], r['integral_r2_std']),
            pm(r['hypervolume_median'], r['hypervolume_std']),
            pm(r['delta_p_median'], r['delta_p_std']),
        ]
        for j, best in enumerate((best_ir2, best_hv, best_dp)):
            if i == best:
                vals[j] = r'\mathbf{' + vals[j] + '}'
        pp = p if i == 0 else ''
        lines.append(f'{pp} & {LABEL[a]} & ${vals[0]}$ & ${vals[1]}$ & ${vals[2]}$\\\\')
    if pi < len(PROBS) - 1:
        lines.append(r'\addlinespace[1pt]')
lines += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
OUT.write_text('\n'.join(lines) + '\n')
print(OUT)
