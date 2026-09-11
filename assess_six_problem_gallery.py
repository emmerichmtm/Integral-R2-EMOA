"""Assess the six diagnostic IR2-EMOA approximation sets.

Metrics: Integral R2, hypervolume, and Delta_p.
Each problem uses a deterministic 10,000-point Pareto reference set.
"""
from __future__ import annotations
import csv, json
from pathlib import Path
import numpy as np
from benchmarks import ideal, pareto_reference
from performance_indicators import ir2_value, hypervolume, delta_p

CASES=[('ZDT1',2,30000),('ZDT2',2,30000),('ZDT3',2,30000),
       ('DTLZ1',3,30000),('DTLZ2',3,30000),('DTLZ7',3,60000)]


def nondominated(A):
    A=np.asarray(A,float); keep=np.ones(len(A),bool)
    for i in range(len(A)):
        if not keep[i]: continue
        dom=np.all(A<=A[i],axis=1)&np.any(A<A[i],axis=1)
        if np.any(dom): keep[i]=False
    return A[keep]


def fixed_hv_reference(R):
    lo=R.min(axis=0); nad=R.max(axis=0); span=np.maximum(nad-lo,1e-12)
    return nad + 0.1*span


def main():
    rows=[]
    result_dir=Path('results/single_shot_exact_ideal_30000')
    for problem,m,budget in CASES:
        f=result_dir/f'{problem.lower()}_{m}d_ir2_seed2027.csv'
        A=np.loadtxt(f,delimiter=',',skiprows=1)
        if A.ndim==1: A=A[None,:]
        A=nondominated(A)
        R=pareto_reference(problem,m,10000)
        # For generators that return slightly more/less than requested, use a
        # deterministic evenly spaced subsample of exactly 10,000 when possible.
        if len(R)>10000:
            idx=np.linspace(0,len(R)-1,10000).round().astype(int); R=R[idx]
        z=np.asarray(ideal(problem,m),float)
        href=fixed_hv_reference(R)
        rows.append(dict(problem=problem,objectives=m,evaluations=budget,
                         n_reference=int(len(R)),integral_r2=ir2_value(A,z),
                         hypervolume=hypervolume(A,href),delta_p=delta_p(A,R,2),
                         ideal=z.tolist(),hv_reference=href.tolist()))
    Path('results').mkdir(exist_ok=True)
    with open('results/six_problem_indicator_assessment.csv','w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    Path('results/six_problem_indicator_assessment.json').write_text(json.dumps(rows,indent=2)+'\n')
    for r in rows: print(r)

if __name__=='__main__': main()
