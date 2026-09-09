from __future__ import annotations
import argparse, json
from pathlib import Path
from make_single_shot_figures import _run_case

ap=argparse.ArgumentParser()
ap.add_argument('problem')
ap.add_argument('m',type=int)
ap.add_argument('--population',type=int,default=60)
ap.add_argument('--evaluations',type=int,default=30000)
ap.add_argument('--seed',type=int,default=2027)
ap.add_argument('--figdir',default='figures')
ap.add_argument('--resultdir',default='results/single_shot_exact_ideal_30000')
args=ap.parse_args()
row=_run_case((args.problem,args.m,args.population,args.evaluations,args.seed,args.figdir,args.resultdir))
out=Path(args.resultdir)/f'{args.problem.lower()}_{args.m}d_metadata.json'
out.write_text(json.dumps(row,indent=2)+'\n')
print(json.dumps(row),flush=True)
