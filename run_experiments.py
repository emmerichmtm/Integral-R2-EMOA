from __future__ import annotations
import argparse,csv,time
from pathlib import Path
import numpy as np
from benchmarks import n_var, pareto_reference
from ir2_emoa import run
from metrics import delta_p, convergence

FIELDS=['problem','objectives','algorithm','run','evaluations','population','delta_p','convergence','runtime_s','nondominated_size']

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--runs',type=int,default=3)
    ap.add_argument('--evaluations',type=int,default=500)
    ap.add_argument('--pop2',type=int,default=30)
    ap.add_argument('--pop3',type=int,default=30)
    ap.add_argument('--problems',nargs='+',default=['DTLZ1','DTLZ2'])
    ap.add_argument('--objectives',nargs='+',type=int,default=[2,3])
    ap.add_argument('--output',default='results/pilot_raw.csv')
    args=ap.parse_args()
    out=Path(args.output); out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',newline='') as f:
      w=csv.DictWriter(f,fieldnames=FIELDS); w.writeheader(); f.flush()
      for prob in args.problems:
       for m in args.objectives:
        ref=pareto_reference(prob,m,2500); pop=args.pop2 if m==2 else args.pop3
        for r in range(args.runs):
          init_rng=np.random.default_rng(100000*m+1000*r+sum(map(ord,prob)))
          init=init_rng.random((pop,n_var(prob,m)))
          for alg in ('ir2','r2','sms'):
            seed=200000*m+1000*r+sum(map(ord,prob))
            t=time.perf_counter(); _,F=run(prob,m,alg,pop,args.evaluations,seed,init); dt=time.perf_counter()-t
            row=dict(problem=prob,objectives=m,algorithm=alg,run=r,evaluations=args.evaluations,population=pop,
                     delta_p=delta_p(F,ref),convergence=convergence(F,ref),runtime_s=dt,nondominated_size=len(F))
            print(row,flush=True); w.writerow(row); f.flush()
if __name__=='__main__': main()
