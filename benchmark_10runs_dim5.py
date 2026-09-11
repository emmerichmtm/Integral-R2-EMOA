"""Ten-run benchmark for IR2-EMOA, finite-weight R2-EMOA, and SMS-EMOA.

Design frozen for the EMO 2027 working paper:
  * problems: ZDT1, ZDT2, ZDT3 (2 objectives), DTLZ1, DTLZ2, DTLZ7 (3 objectives)
  * decision-space dimension: 5 for every problem
  * population: 50
  * evaluations: 10,000 per run
  * runs: 10 paired/common-start runs
  * algorithms: IR2-EMOA, finite-weight R2-EMOA, SMS-EMOA
  * assessment: exact Integral R2, hypervolume, and Delta_p (p=2)

The final approximation is evaluated against a deterministic 10,000-point
Pareto-front reference set. Hypervolume uses one fixed reference point per
problem: the sampled PF nadir shifted outward by 10% of its ideal-to-nadir
range. The same reference is used for all algorithms and all runs.

IR2 selection and evaluation use the exact low-dimensional implementation
based on perspective mapping distributed with this repository and associated
with https://github.com/emmerichmtm/IntegralR2ByPerspectiveMapping .
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from benchmarks import ideal, pareto_reference
from ir2_emoa import run
from performance_indicators import delta_p, hypervolume, ir2_value

CASES = (
    ("ZDT1", 2), ("ZDT2", 2), ("ZDT3", 2),
    ("DTLZ1", 3), ("DTLZ2", 3), ("DTLZ7", 3),
)
ALGORITHMS = ("ir2", "r2", "sms")
ALG_LABEL = {"ir2": "IR2-EMOA", "r2": "finite R2-EMOA", "sms": "SMS-EMOA"}


def exact_reference(problem: str, m: int, n: int = 10000) -> np.ndarray:
    """Return a deterministic reference set with exactly n points."""
    R = np.asarray(pareto_reference(problem, m, n), float)
    if len(R) == n:
        return R
    if len(R) > n:
        # deterministic coverage-preserving subsample
        idx = np.linspace(0, len(R)-1, n).round().astype(int)
        return R[idx]
    # If a generator returns fewer than requested, regenerate with a larger
    # target until a deterministic n-point subsample is available.
    target = max(n + 1, int(math.ceil(1.25*n)))
    for _ in range(8):
        R = np.asarray(pareto_reference(problem, m, target), float)
        if len(R) >= n:
            idx = np.linspace(0, len(R)-1, n).round().astype(int)
            return R[idx]
        target = int(math.ceil(1.5*target))
    raise RuntimeError(f"Could not generate {n} reference points for {problem}, m={m}")


def hv_reference(problem: str, m: int, R: np.ndarray) -> np.ndarray:
    """Fixed slightly shifted PF nadir used only for performance HV."""
    z = np.asarray(ideal(problem, m), float)
    nad = np.max(R, axis=0)
    span = np.maximum(nad - z, 1e-12)
    return nad + 0.1 * span


EXPLICIT_SEEDS = {
    'ZDT1': [3304432, 3305441, 3306450, 3307459, 3308468, 3309477, 3310486, 3311495, 3312504, 3313513],
    'ZDT2': [3312351, 3313360, 3314369, 3315378, 3316387, 3317396, 3318405, 3319414, 3320423, 3321432],
    'ZDT3': [3320270, 3321279, 3322288, 3323297, 3324306, 3325315, 3326324, 3327333, 3328342, 3329351],
    'DTLZ1': [3906276, 3907285, 3908294, 3909303, 3910312, 3911321, 3912330, 3913339, 3914348, 3915357],
    'DTLZ2': [3914195, 3915204, 3916213, 3917222, 3918231, 3919240, 3920249, 3921258, 3922267, 3923276],
    'DTLZ7': [3953790, 3954799, 3955808, 3956817, 3957826, 3958835, 3959844, 3960853, 3961862, 3962871],
}

def _seed_code(problem: str, run_index: int) -> int:
    return EXPLICIT_SEEDS[problem][run_index]


def _job(task):
    problem, m, run_index, evaluations, population, dimension, ref_list, hvref_list, outdir_s = task
    outdir = Path(outdir_s)
    outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir / f"{problem.lower()}_run{run_index:02d}.json"
    if outfile.exists():
        return json.loads(outfile.read_text())

    seed = _seed_code(problem, run_index)
    initialization_seed = seed + 17
    init_rng = np.random.default_rng(initialization_seed)
    initial_X = init_rng.random((population, dimension))
    R = np.asarray(ref_list, float)
    href = np.asarray(hvref_list, float)
    z = np.asarray(ideal(problem, m), float)

    rows = []
    # Same initial population and same variation RNG seed across algorithms.
    for alg in ALGORITHMS:
        _, F = run(problem, m, alg, population, evaluations, seed, initial_X=initial_X, n_variables=dimension)
        F = np.asarray(F, float)
        row = {
            "problem": problem,
            "objectives": m,
            "decision_dimension": dimension,
            "algorithm": alg,
            "algorithm_label": ALG_LABEL[alg],
            "run": run_index,
            "seed": seed,
            "initialization_seed": initialization_seed,
            "evaluations": evaluations,
            "population": population,
            "nondominated_size": int(len(F)),
            "integral_r2": float(ir2_value(F, z)),
            "hypervolume": float(hypervolume(F, href)),
            "delta_p": float(delta_p(F, R, p=2)),
            "ideal": [float(v) for v in z],
            "hv_reference": [float(v) for v in href],
            "finite_r2_weights": (101 if (alg == "r2" and m == 2) else 120 if alg == "r2" else None),
        }
        rows.append(row)
    outfile.write_text(json.dumps(rows, indent=2) + "\n")
    return rows


def _iqr(values):
    a = np.asarray(values, float)
    return float(np.quantile(a, 0.75) - np.quantile(a, 0.25))


def summarize(rows):
    summary=[]
    for problem,m in CASES:
        for alg in ALGORITHMS:
            S=[r for r in rows if r["problem"]==problem and r["algorithm"]==alg]
            if not S: continue
            d={
                "problem": problem,
                "objectives": m,
                "decision_dimension": S[0]["decision_dimension"],
                "algorithm": alg,
                "algorithm_label": ALG_LABEL[alg],
                "runs": len(S),
                "evaluations": S[0]["evaluations"],
                "population": S[0]["population"],
            }
            for key in ("integral_r2", "hypervolume", "delta_p"):
                vals=[x[key] for x in S]
                d[f"{key}_median"] = float(np.median(vals))
                d[f"{key}_iqr"] = _iqr(vals)
                d[f"{key}_mean"] = float(np.mean(vals))
                d[f"{key}_std"] = float(np.std(vals, ddof=1)) if len(vals)>1 else 0.0
            summary.append(d)
    return summary


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--runs", type=int, default=10)
    ap.add_argument("--evaluations", type=int, default=10000)
    ap.add_argument("--population", type=int, default=50)
    ap.add_argument("--dimension", type=int, default=5)
    ap.add_argument("--reference-points", type=int, default=10000)
    ap.add_argument("--jobs", type=int, default=min(5, os.cpu_count() or 1))
    ap.add_argument("--outdir", default="results/benchmark_10runs_dim5")
    args=ap.parse_args()
    if args.dimension < 3:
        raise ValueError("dimension must be at least 3 for the three-objective DTLZ cases")

    outdir=Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    refs={}; hrefs={}
    for problem,m in CASES:
        R=exact_reference(problem,m,args.reference_points)
        refs[(problem,m)]=R
        hrefs[(problem,m)]=hv_reference(problem,m,R)
        np.savetxt(outdir/f"reference_{problem.lower()}_{m}d.csv",R,delimiter=",",
                   header=",".join(f"f{i+1}" for i in range(m)),comments="")

    tasks=[]
    for problem,m in CASES:
        for r in range(args.runs):
            tasks.append((problem,m,r,args.evaluations,args.population,args.dimension,
                          refs[(problem,m)].tolist(),hrefs[(problem,m)].tolist(),str(outdir/"partial")))

    nested=[]
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs={ex.submit(_job,t): (t[0],t[2]) for t in tasks}
        for fut in as_completed(futs):
            problem,r=futs[fut]
            result=fut.result(); nested.append(result)
            print(f"completed {problem} run {r}: " + ", ".join(
                f"{x['algorithm_label']} Delta_p={x['delta_p']:.6g}" for x in result), flush=True)

    rows=[x for group in nested for x in group]
    rows.sort(key=lambda x:(x["problem"],x["run"],ALGORITHMS.index(x["algorithm"])))
    fields=[
        "problem","objectives","decision_dimension","algorithm","algorithm_label","run","seed","initialization_seed",
        "evaluations","population","nondominated_size","integral_r2","hypervolume","delta_p",
        "ideal","hv_reference","finite_r2_weights"
    ]
    with (outdir/"raw.csv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields); w.writeheader(); w.writerows(rows)
    (outdir/"raw.json").write_text(json.dumps(rows,indent=2)+"\n")

    summary=summarize(rows)
    sfields=list(summary[0].keys())
    with (outdir/"summary.csv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=sfields); w.writeheader(); w.writerows(summary)
    (outdir/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")

    with (outdir/"seeds.csv").open("w", newline="") as fh:
        w=csv.writer(fh)
        w.writerow(["problem","run","variation_seed","initialization_seed"])
        for problem,m in CASES:
            for run_index, seed in enumerate(EXPLICIT_SEEDS[problem][:args.runs]):
                w.writerow([problem, run_index, seed, seed+17])

    metadata={
        "problems":[{"problem":p,"objectives":m} for p,m in CASES],
        "algorithms":[ALG_LABEL[a] for a in ALGORITHMS],
        "runs":args.runs,"evaluations":args.evaluations,"population":args.population,
        "decision_dimension":args.dimension,"reference_points":args.reference_points,
        "paired_common_start":True,
        "seed_file":"seeds.csv",
        "explicit_variation_seeds":EXPLICIT_SEEDS,
        "initialization_seed_rule":"variation_seed + 17",
        "finite_r2_weights":{"2d":101,"3d":120},
        "performance_indicators":["Integral R2","hypervolume","Delta_p (p=2)"],
        "hv_reference_rule":"sampled PF nadir + 0.1*(sampled PF nadir - exact ideal)",
        "integral_r2_source":"https://github.com/emmerichmtm/IntegralR2ByPerspectiveMapping",
        "algorithm_source":"https://github.com/emmerichmtm/Integral-R2-EMOA",
    }
    (outdir/"metadata.json").write_text(json.dumps(metadata,indent=2)+"\n")
    print("wrote", outdir/"raw.csv")
    print("wrote", outdir/"summary.csv")

if __name__=="__main__":
    main()
