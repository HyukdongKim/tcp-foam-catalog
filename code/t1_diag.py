"""Diagnose the best stored candidate of a T1 search: relax at several resolutions and print the
area of each polygon class relative to the mean face area (a collapsing triangle shows ratio ~0).
usage: python code/t1_diag.py PATH_TO_t1_jsonl(.gz) r1 [r2 ...]"""
import gzip, json, os, sys, time, warnings
import numpy as np
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t1search import deser
from trifoam import TriFoam
from relax import Mesh, relax_box

fn = sys.argv[1]; rs = [int(x) for x in sys.argv[2:]]
op = gzip.open if fn.endswith('.gz') else open
R = [json.loads(l) for l in op(fn, 'rt')]
b = min([x for x in R if 'tets' in x], key=lambda x: x['cost'])
f = TriFoam(np.array(b['X']), np.array(b['H']), deser(b['tets']))
cells = f.cells()
sizes = []
for i, c in enumerate(cells):
    for fc in c['faces']:
        j, n = fc['j'], tuple(int(v) for v in fc['n'])
        if j < i or (j == i and n < (0, 0, 0)):
            continue
        sizes.append(len(fc['verts']))
sizes = np.array(sizes)
for r in rs:
    t = time.time(); m = Mesh(f.X, f.H, cells, r=r); c = relax_box(m, gtol=1e-10, tol=1e-10, mu0=300, maxit=40)
    a, bb, cc = m._abs(m.P); ar = .5 * np.linalg.norm(np.cross(bb - a, cc - a), axis=1)
    FA = np.bincount(m.tface, weights=ar, minlength=m.nfaces); mean = FA.mean()
    print(f"r={r} cost={c:.6f} t={time.time()-t:.0f}s", flush=True)
    for k in range(3, 9):
        if (sizes == k).any():
            print(f"   {k}-gon: n={int((sizes==k).sum())} area/mean min={FA[sizes==k].min()/mean:.4f} max={FA[sizes==k].max()/mean:.4f}", flush=True)
