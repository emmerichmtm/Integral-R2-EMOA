"""Reproduce single-shot IR2-EMOA figures with exact ideal points.

ZDT1/2/3 are drawn in 2D. DTLZ1/2/7 are drawn in 3D.
No artificial strict-dominance shift is used for the ideal point.
A zero loss is valid and the Integral R2 code maps its reciprocal to +infinity.
"""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from benchmarks import ideal, pareto_reference
from ir2_emoa import run
from metrics import delta_p, convergence

CASES = (
    ('ZDT1', 2), ('ZDT2', 2), ('ZDT3', 2),
    ('DTLZ1', 3), ('DTLZ2', 3), ('DTLZ7', 3),
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--population', type=int, default=60)
    ap.add_argument('--evaluations', type=int, default=3000)
    ap.add_argument('--seed', type=int, default=2027)
    ap.add_argument('--figdir', default='figures')
    ap.add_argument('--resultdir', default='results/single_shot_exact_ideal_60')
    args = ap.parse_args()

    figdir = Path(args.figdir); figdir.mkdir(parents=True, exist_ok=True)
    resultdir = Path(args.resultdir); resultdir.mkdir(parents=True, exist_ok=True)
    metadata = []

    for problem, m in CASES:
        z = np.asarray(ideal(problem, m), float)
        _, F = run(problem, m, 'ir2', args.population, args.evaluations, args.seed)
        ref = pareto_reference(problem, m, 20001 if m == 2 else 4000)

        csv_path = resultdir / f'{problem.lower()}_{m}d_ir2_seed{args.seed}.csv'
        with csv_path.open('w', newline='') as fh:
            w = csv.writer(fh)
            w.writerow([f'f{i+1}' for i in range(m)])
            w.writerows(F.tolist())

        if m == 2:
            plt.figure(figsize=(6.2, 4.8))
            plt.scatter(ref[:,0], ref[:,1], s=3, label='Pareto front')
            plt.scatter(F[:,0], F[:,1], s=24, label='IR2-EMOA')
            plt.scatter([z[0]], [z[1]], s=95, marker='*', label=r'exact ideal $z^\star$')
            plt.xlabel(r'$f_1$'); plt.ylabel(r'$f_2$')
            plt.title(f'{problem}: IR2-EMOA single shot')
            plt.legend(); plt.tight_layout()
        else:
            fig = plt.figure(figsize=(7.0, 5.6))
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(ref[:,0], ref[:,1], ref[:,2], s=3, label='Pareto front')
            ax.scatter(F[:,0], F[:,1], F[:,2], s=18, label='IR2-EMOA')
            ax.scatter([z[0]], [z[1]], [z[2]], s=120, marker='*', label=r'exact ideal $z^\star$')
            ax.set_xlabel(r'$f_1$'); ax.set_ylabel(r'$f_2$'); ax.set_zlabel(r'$f_3$')
            ax.set_title(f'{problem}: IR2-EMOA single shot')
            ax.legend(loc='upper left')
            plt.tight_layout()
        png = figdir / f'single_shot_{problem.lower()}_{m}d_exact_ideal_60pts.png'
        plt.savefig(png, dpi=240, bbox_inches='tight'); plt.close()

        row = dict(problem=problem, objectives=m, population=args.population,
                   evaluations=args.evaluations, seed=args.seed,
                   ideal=[float(v) for v in z],
                   nondominated_size=int(len(F)),
                   delta_p=float(delta_p(F, ref)),
                   convergence=float(convergence(F, ref)),
                   figure=str(png), objectives_csv=str(csv_path))
        metadata.append(row)
        print(row, flush=True)

    (resultdir / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')

if __name__ == '__main__':
    main()
