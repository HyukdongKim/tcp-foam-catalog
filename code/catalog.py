"""High-resolution catalogue of equal-volume foam costs for every structure in tcp_all.json.
usage: python catalog.py [workers=10]   -> rerun/catalog.jsonl  (resumable)
Phase A: r=2 for all.  Phase B: r=3 for all that were valid at r=2.
Phase C: r=4 for valid structures with r=3 cost < 5.310.
Each job relaxes with mu0=300; if the result is not valid it is retried with mu0=1.
valid = flipped-triangle area fraction < 1e-6 and max volume error < 1e-5."""
import sys, os, json, time, signal, warnings
import numpy as np
from multiprocessing import Pool
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OUT = 'rerun/catalog.jsonl'
LIMIT = {2: 6 * 3600, 3: 12 * 3600, 4: 20 * 3600}


def _alarm(s, f):
    raise TimeoutError()


def one(rec, r, mu0):
    from tcpbuild import build_cells
    from relax import Mesh, relax_box
    S = [set((j, tuple(n)) for j, n in s) for s in rec['S']]
    H = np.array(rec['H'])
    X, cells = build_cells(np.array(rec['F']), H, np.array(rec['Z']), S)
    m = Mesh(X, H, cells, r=r)
    c = float(relax_box(m, gtol=1e-9, tol=1e-9, mu0=mu0, maxit=30))
    a, b, cc = m._abs(m.P)
    cr = np.cross(b - a, cc - a); ar = 0.5 * np.linalg.norm(cr, axis=1)
    bad = np.einsum('ij,ij->i', cr, m.face_normals[m.tface]) < 0
    A, V = m.energy(m.P, False)
    out = dict(cost=c, flip_frac=float(ar[bad].sum() / ar.sum()),
               vol_err=float(np.abs(V / m.Vt - 1).max()), tris=int(len(m.tid)), Hrel=m.H.tolist())
    out['ok'] = bool(np.isfinite(c) and out['flip_frac'] < 1e-6 and out['vol_err'] < 1e-5)
    return out


def job(a):
    rec, r = a
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(LIMIT[r])
    res = dict(name=rec['name'], N=rec['N'], meanZ=rec['meanZ'], r=r); t = time.time()
    try:
        o = one(rec, r, 300.0); o['mu0'] = 300.0
        if not o['ok']:
            o2 = one(rec, r, 1.0); o2['mu0'] = 1.0
            if o2['ok'] or o2['cost'] < o['cost']:
                o = o2
        res.update(o)
    except Exception as e:
        res['error'] = repr(e)[:200]; res['ok'] = False
    signal.alarm(0); res['t'] = time.time() - t
    return res


def load():
    if not os.path.exists(OUT):
        return {}
    return {(d['name'], d['r']): d for d in map(json.loads, open(OUT))}


def phase(recs, r, nw):
    done = load()
    todo = [x for x in recs if (x['name'], r) not in done]
    todo.sort(key=lambda x: -x['N'])
    print(f'[r={r}] {len(todo)} structures to relax', flush=True)
    if not todo:
        return
    with Pool(nw, maxtasksperchild=1) as pool, open(OUT, 'a') as f:
        for k, res in enumerate(pool.imap_unordered(job, [(x, r) for x in todo]), 1):
            f.write(json.dumps(res) + '\n'); f.flush()
            tag = f"{res['cost']:.6f} ok={res['ok']}" if 'cost' in res else res.get('error')
            print(f"[r={r}] {k}/{len(todo)} {res['name']:24s} N={res['N']:4d} {tag} ({res['t']:.0f}s)", flush=True)


if __name__ == '__main__':
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    os.chdir(os.path.join(HERE, '..')); os.makedirs('rerun', exist_ok=True)
    recs = json.load(open('data/tcp_all.json'))
    phase(recs, 2, nw)
    d = load()
    phase([x for x in recs if d.get((x['name'], 2), {}).get('ok')], 3, nw)
    d = load()
    phase([x for x in recs if d.get((x['name'], 3), {}).get('ok') and d[(x['name'], 3)]['cost'] < 5.310], 4, nw)
    print('[완료] 130종 고해상도 계산이 끝났습니다.', flush=True)
