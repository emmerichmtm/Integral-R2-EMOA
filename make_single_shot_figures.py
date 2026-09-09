"""Reproduce IR2-EMOA single-shot figures with exact ideal points.

ZDT1/2/3 are bi-objective. DTLZ1/2/7 are tri-objective.
The population size is 60. The default budget is 30,000 function
 evaluations per run, except DTLZ7 which uses 60,000 by default. No artificial strict-dominance shift is applied to the
ideal point; zero loss is valid and maps to reciprocal +infinity.

The 3-D plots explicitly use the same perspective projection and Matplotlib
view as the earlier figures: elevation 30 degrees and azimuth -60 degrees.
This revision also increases the visibility of the approximation points by
using larger markers, black edges, and a light gray reference front.
"""
from __future__ import annotations
import argparse, csv, json, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from benchmarks import ideal, pareto_reference
from ir2_emoa import run
from metrics import delta_p, convergence

CASES = (
    ('ZDT1', 2), ('ZDT2', 2), ('ZDT3', 2),
    ('DTLZ1', 3), ('DTLZ2', 3), ('DTLZ7', 3),
)

VIEW_ELEV = 30.0
VIEW_AZIM = -60.0


def _run_case(task):
    problem, m, population, evaluations, seed, figdir_s, resultdir_s = task
    figdir = Path(figdir_s); resultdir = Path(resultdir_s)
    figdir.mkdir(parents=True, exist_ok=True)
    resultdir.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()
    z = np.asarray(ideal(problem, m), float)
    _, F = run(problem, m, 'ir2', population, evaluations, seed)
    ref = pareto_reference(problem, m, 20001 if m == 2 else 4096)

    csv_path = resultdir / f'{problem.lower()}_{m}d_ir2_seed{seed}.csv'
    with csv_path.open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow([f'f{i+1}' for i in range(m)])
        w.writerows(F.tolist())

    if m == 2:
        plt.figure(figsize=(6.2, 4.8))
        plt.scatter(ref[:,0], ref[:,1], s=5, c='0.78', alpha=0.65, label='Pareto front')
        plt.scatter(F[:,0], F[:,1], s=42, c='#1f77b4', edgecolors='k', linewidths=0.45, label='IR2-EMOA', zorder=3)
        plt.scatter([z[0]], [z[1]], s=145, marker='*', c='gold', edgecolors='k', linewidths=0.9, label=r'exact ideal $z^\star$', zorder=4)
        plt.xlabel(r'$f_1$'); plt.ylabel(r'$f_2$')
        plt.title(f'{problem}: IR2-EMOA single shot')
        plt.legend(); plt.tight_layout()
    else:
        fig = plt.figure(figsize=(7.0, 5.6))
        # Disable automatic depth-based artist ordering so the true PF can be
        # drawn as a background and the approximation set can remain visible
        # in the foreground even when it lies just behind the front in the
        # current camera projection.
        ax = fig.add_subplot(111, projection='3d', computed_zorder=False)
        ax.scatter(ref[:,0], ref[:,1], ref[:,2], s=3, c='0.84', alpha=0.10,
                   depthshade=False, zorder=1, label='Pareto front')
        ax.scatter(F[:,0], F[:,1], F[:,2], s=42, c='#1f77b4', edgecolors='k',
                   linewidths=0.5, depthshade=False, zorder=50, label='IR2-EMOA')
        ax.scatter([z[0]], [z[1]], [z[2]], s=190, marker='*', c='gold',
                   edgecolors='k', linewidths=0.9, depthshade=False, zorder=60,
                   label=r'exact ideal $z^\star$')
        ax.set_proj_type('persp')
        ax.view_init(elev=VIEW_ELEV, azim=VIEW_AZIM)
        ax.set_xlabel(r'$f_1$'); ax.set_ylabel(r'$f_2$'); ax.set_zlabel(r'$f_3$')
        ax.set_title(f'{problem}: IR2-EMOA single shot')
        ax.legend(loc='upper left')
        plt.tight_layout()

    png = figdir / f'single_shot_{problem.lower()}_{m}d_exact_ideal_{evaluations}evals.png'
    plt.savefig(png, dpi=240, bbox_inches='tight'); plt.close()

    # Validate the defining requirement for an ideal point on the plotted PF.
    if not np.all(ref >= z[None,:] - 2e-12):
        raise AssertionError(f'{problem}: ideal point does not dominate PF reference sample')

    return dict(
        problem=problem, objectives=m, population=population,
        evaluations=evaluations, seed=seed,
        ideal=[float(v) for v in z],
        nondominated_size=int(len(F)),
        delta_p=float(delta_p(F,ref)),
        convergence=float(convergence(F,ref)),
        runtime_s=float(time.perf_counter()-t0),
        figure=str(png), objectives_csv=str(csv_path),
        view=(dict(projection='persp', elevation=VIEW_ELEV, azimuth=VIEW_AZIM)
              if m == 3 else None),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--population', type=int, default=60)
    ap.add_argument('--evaluations', type=int, default=30000)
    ap.add_argument('--dtlz7-evaluations', type=int, default=60000,
                    help='override evaluation budget for the DTLZ7 3D single-shot figure')
    ap.add_argument('--seed', type=int, default=2027)
    ap.add_argument('--jobs', type=int, default=1,
                    help='independent cases to run in parallel')
    ap.add_argument('--figdir', default='figures')
    ap.add_argument('--resultdir', default='results/single_shot_exact_ideal_30000')
    args = ap.parse_args()

    Path(args.figdir).mkdir(parents=True, exist_ok=True)
    resultdir = Path(args.resultdir); resultdir.mkdir(parents=True, exist_ok=True)
    tasks=[]
    for p, m in CASES:
        evals = args.dtlz7_evaluations if (p == 'DTLZ7' and m == 3) else args.evaluations
        tasks.append((p, m, args.population, evals, args.seed, args.figdir, args.resultdir))

    rows=[]
    if args.jobs == 1:
        for task in tasks:
            row=_run_case(task); rows.append(row); print(row, flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futures={ex.submit(_run_case,t):t[:2] for t in tasks}
            for fut in as_completed(futures):
                row=fut.result(); rows.append(row); print(row, flush=True)

    rows.sort(key=lambda r:(r['objectives'],r['problem']))
    (resultdir/'metadata.json').write_text(json.dumps(rows,indent=2)+'\n')


if __name__=='__main__':
    main()
