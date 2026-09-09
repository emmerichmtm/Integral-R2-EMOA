from __future__ import annotations
import csv,statistics
from collections import defaultdict
from pathlib import Path

IN='results/pilot_raw.csv'; OUT='results/pilot_summary.csv'
rows=list(csv.DictReader(open(IN)))
g=defaultdict(list)
for r in rows: g[(r['problem'],int(r['objectives']),r['algorithm'])].append(r)
fields=['problem','objectives','algorithm','runs','delta_p_median','delta_p_iqr','convergence_median','convergence_iqr','runtime_median_s']
out=[]
def iqr(v):
    q=statistics.quantiles(sorted(v),n=4,method='inclusive'); return q[2]-q[0]
for k,rs in sorted(g.items()):
    d=[float(x['delta_p']) for x in rs]; c=[float(x['convergence']) for x in rs]; t=[float(x['runtime_s']) for x in rs]
    out.append(dict(problem=k[0],objectives=k[1],algorithm=k[2],runs=len(rs),
        delta_p_median=statistics.median(d),delta_p_iqr=iqr(d),convergence_median=statistics.median(c),convergence_iqr=iqr(c),runtime_median_s=statistics.median(t)))
Path(OUT).parent.mkdir(exist_ok=True)
with open(OUT,'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(out)
for r in out: print(r)
