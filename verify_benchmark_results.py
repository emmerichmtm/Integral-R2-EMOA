"""Verify completeness, seed pairing, and summary statistics of the frozen benchmark."""
from pathlib import Path
import csv, math
import numpy as np

ROOT=Path('results/benchmark_10runs_dim5')
raw=list(csv.DictReader((ROOT/'raw.csv').open()))
summary=list(csv.DictReader((ROOT/'summary.csv').open()))
assert len(raw)==180, len(raw)
problems=['ZDT1','ZDT2','ZDT3','DTLZ1','DTLZ2','DTLZ7']
algs=['ir2','r2','sms']
for p in problems:
    for run in range(10):
        rr=[r for r in raw if r['problem']==p and int(r['run'])==run]
        assert len(rr)==3, (p,run,len(rr))
        assert len({r['seed'] for r in rr})==1, (p,run,'variation seed mismatch')
        assert len({r['initialization_seed'] for r in rr})==1, (p,run,'initialization seed mismatch')
        assert all(int(r['nondominated_size'])==50 for r in rr), (p,run,'ND size')
for s in summary:
    rr=[r for r in raw if r['problem']==s['problem'] and r['algorithm']==s['algorithm']]
    assert len(rr)==10
    for key in ['integral_r2','hypervolume','delta_p']:
        a=np.array([float(r[key]) for r in rr])
        med=float(np.median(a)); sd=float(np.std(a,ddof=1))
        assert math.isclose(med,float(s[f'{key}_median']),rel_tol=1e-12,abs_tol=1e-14)
        assert math.isclose(sd,float(s[f'{key}_std']),rel_tol=1e-12,abs_tol=1e-14)
print('OK: 180 records; 18 groups x 10 runs; paired seeds; median/SD summary verified.')
