"""r=2,3 relaxation of extra TCP records with the settings of catalog.py.
usage: python code/shoemaker_run.py [workers=2] [input=data/tcp_shoemaker.json] [output=rerun/shoemaker_phases.jsonl]
Archived outputs: results/shoemaker_phases.jsonl (data/tcp_shoemaker.json), results/cif_phases.jsonl (data/tcp_cif.json)."""
import json, os, sys, time, warnings
import numpy as np
from multiprocessing import Pool
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)


def one(rec, r, mu0):
    from tcpbuild import build_cells
    from relax import Mesh, relax_box
    S = [set((j, tuple(n)) for j, n in s) for s in rec['S']]; H = np.array(rec['H'])
    X, cells = build_cells(np.array(rec['F']), H, np.array(rec['Z']), S)
    m = Mesh(X, H, cells, r=r); c = float(relax_box(m, gtol=1e-9, tol=1e-9, mu0=mu0, maxit=30))
    a, b, cc = m._abs(m.P); cr = np.cross(b - a, cc - a); ar = .5 * np.linalg.norm(cr, axis=1)
    bad = np.einsum('ij,ij->i', cr, m.face_normals[m.tface]) < 0; A, V = m.energy(m.P, False)
    o = dict(cost=c, flip_frac=float(ar[bad].sum() / ar.sum()), vol_err=float(np.abs(V / m.Vt - 1).max()))
    o['ok'] = bool(np.isfinite(c) and o['flip_frac'] < 1e-6 and o['vol_err'] < 1e-5)
    return o


def job(a):
    rec, r = a; t = time.time(); o = one(rec, r, 300.)
    if not o['ok']:
        o2 = one(rec, r, 1.)
        if o2['ok'] or o2['cost'] < o['cost']:
            o = o2
    return dict(name=rec['name'], N=rec['N'], meanZ=rec['meanZ'], r=r, t=time.time() - t, **o)


if __name__ == '__main__':
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    os.chdir(os.path.join(HERE, '..')); os.makedirs('rerun', exist_ok=True)
    src = sys.argv[2] if len(sys.argv) > 2 else 'data/tcp_shoemaker.json'
    dst = sys.argv[3] if len(sys.argv) > 3 else 'rerun/shoemaker_phases.jsonl'
    recs = sorted(json.load(open(src)), key=lambda x: x['N'])
    with Pool(nw) as p, open(dst, 'a') as f:
        for res in p.imap_unordered(job, [(x, r) for r in (2, 3) for x in recs]):
            f.write(json.dumps(res) + '\n'); f.flush(); print(res, flush=True)
