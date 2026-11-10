"""T2-type moves (cell number changes by one): insert one cell centre at a foam vertex, at a foam-edge
centroid or at a face centre, or delete one cell; re-derive the Voronoi topology by periodic Delaunay and
relax at r = 1 with all N +/- 1 cells at equal volume (box free).
usage: python code/t2_probe.py START REPS WORKERS  -> rerun/t2p_<tag>.jsonl
Follow-up T1 rearrangement of a stored state: python code/t1_deep.py file:rerun/t2p_<tag>.jsonl:LINE 1x1x1 3 10"""
import json, os, sys, warnings
import numpy as np
from multiprocessing import Pool
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from trifoam import from_tcp, from_delaunay
from t1_deep import wl_hash
from t1search import evaluate, ser


def candidates(f0):
    Hinv = np.linalg.inv(f0.H)
    pts = [('vertex', np.mean([f0.pos(v) for v in T], 0)) for T in f0.tets]
    pts += [('edge', np.mean([f0.pos(v) for v in c], 0)) for c in f0.tri]
    pts += [('face', np.mean([f0.pos(v) for v in c], 0)) for c in f0.edge]
    F0 = np.mod(f0.X @ Hinv, 1)
    out = [(f'insert-{k}', np.vstack([F0, np.mod(p @ Hinv, 1)])) for k, p in pts]
    out += [(f'delete-{i}', np.delete(F0, i, 0)) for i in range(f0.N)]
    return out


if __name__ == '__main__':
    start, reps, nw = sys.argv[1], tuple(int(x) for x in sys.argv[2].split('x')), int(sys.argv[3])
    os.chdir(os.path.join(HERE, '..')); os.makedirs('rerun', exist_ok=True)
    R = {d['name']: d for d in json.load(open('data/tcp_all.json'))}
    f0 = from_tcp(R[start]).supercell(reps)
    seen, jobs = set(), []
    for kind, F in candidates(f0):
        try:
            g = from_delaunay(F, f0.H)
        except Exception:
            continue
        if g.check() is not None:
            continue
        h = wl_hash(g)
        if h not in seen:
            seen.add(h); jobs.append((kind, g))
    tag = f"{start.split(':')[-1]}_{'x'.join(map(str, reps))}"
    print(f'[{tag}] N={f0.N}; {len(jobs)} distinct T2-type topologies', flush=True)
    with Pool(nw) as p:
        base = p.apply(evaluate, ((f0.X, f0.H, sorted(f0.tets), 1),))['cost']
        res = p.map(evaluate, [(g.X, g.H, sorted(g.tets), 1) for _, g in jobs])
    with open(f'rerun/t2p_{tag}.jsonl', 'w') as out:
        for (kind, g), r in zip(jobs, res):
            rec = dict(kind=kind, N=g.N, cost=r['cost'], ok=r['ok'], faces=g.face_hist(), nf=g.nfaces_hist(), base=base)
            if r['ok']:
                rec.update(tets=ser(g.tets), X=r['Xn'], H=r['Hn'])
            out.write(json.dumps(rec) + '\n')
    for (kind, g), r in sorted(zip(jobs, res), key=lambda t: t[1]['cost'])[:5]:
        print(f"  {r['cost']:.6f}  defect energy {(r['cost'] - base) * g.N:+.3f}  {kind}  N={g.N}", flush=True)
