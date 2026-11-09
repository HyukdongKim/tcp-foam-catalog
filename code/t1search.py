"""Beam search over foam topologies by T1 moves (Pachner 2-3 / 3-2 on the dual triangulation),
judged by curved equal-volume relaxed cost at resolution r (free box).
usage: python t1search.py START REPS DEPTH BEAM WORKERS [r=1]
  START = name in tcp_all.json (e.g. alloy:Cr3Si_A15, n20:Phys26, n20:N56); REPS = e.g. 2x1x1
Writes rerun/t1_<tag>.jsonl (one line per evaluated topology; resumable is not needed)."""
import sys, os, json, time, signal, warnings
import numpy as np
from multiprocessing import Pool
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from trifoam import TriFoam, from_tcp, relaxed_centroids

EVAL_LIMIT = 1800


def _alarm(s, f):
    raise TimeoutError()


def evaluate(args):
    X, H, tets, r = args
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(EVAL_LIMIT)
    from relax import Mesh, relax_box
    t = time.time()
    try:
        f = TriFoam(X, H, []); f.tets = set(tets); f._index()
        best = None
        for mu0 in (300.0, 1.0):
            m = Mesh(f.X, f.H, f.cells(), r=r)
            c = float(relax_box(m, gtol=1e-9, tol=1e-9, mu0=mu0, maxit=30))
            a, b, cc = m._abs(m.P)
            cr = np.cross(b - a, cc - a); ar = 0.5 * np.linalg.norm(cr, axis=1)
            bad = np.einsum('ij,ij->i', cr, m.face_normals[m.tface]) < 0
            A, V = m.energy(m.P, False)
            o = dict(cost=c, flip=float(ar[bad].sum() / ar.sum()), verr=float(np.abs(V / m.Vt - 1).max()))
            o['ok'] = bool(np.isfinite(c) and o['flip'] < 1e-6 and o['verr'] < 1e-5)
            if o['ok']:
                o['Xn'] = relaxed_centroids(m).tolist(); o['Hn'] = m.H.tolist()
            if best is None or (o['ok'] and not best['ok']) or (o['ok'] == best['ok'] and o['cost'] < best['cost']):
                best = o
            if o['ok']:
                break
    except Exception as e:
        best = dict(ok=False, error=repr(e)[:200], cost=float('inf'))
    signal.alarm(0)
    best['t'] = time.time() - t
    return best


def ser(tets):
    return [[[v[0], list(v[1])] for v in T] for T in sorted(tets)]


def deser(L):
    return [tuple((v[0], tuple(v[1])) for v in T) for T in L]


def main():
    start, reps, depth, beam, nw = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    r = int(sys.argv[6]) if len(sys.argv) > 6 else 1
    reps = tuple(int(x) for x in reps.split('x'))
    os.chdir(os.path.join(HERE, '..')); os.makedirs('rerun', exist_ok=True)
    tag = f"{start.split(':')[-1]}_{'x'.join(map(str, reps))}_r{r}"
    out = open(f'rerun/t1_{tag}.jsonl', 'a')
    R = {d['name']: d for d in json.load(open('data/tcp_all.json'))}
    f0 = from_tcp(R[start]).supercell(reps)
    assert f0.check() is None
    pool = Pool(nw, maxtasksperchild=20)
    res0 = pool.apply(evaluate, ((f0.X, f0.H, sorted(f0.tets), r),))
    base = res0['cost']
    print(f'[{tag}] start N={f0.N} cost={base:.6f} ok={res0["ok"]}', flush=True)
    rec = dict(tag=tag, depth=0, path=[], cost=base, ok=res0['ok'], faces=f0.face_hist(), nf=f0.nfaces_hist(), sig=str(f0.signature()))
    out.write(json.dumps(rec) + '\n'); out.flush()
    seen = {f0.signature()}
    level = [(f0, res0, [])]
    best_overall = (base, [])
    for d in range(1, depth + 1):
        cands = []
        for f, res, path in level:
            Xp = np.array(res['Xn']) if res.get('Xn') else f.X
            Hp = np.array(res['Hn']) if res.get('Hn') else f.H
            for mv in f.moves():
                g = f.apply(mv)
                if g.check() is not None:
                    continue
                s = g.signature()
                if s in seen:
                    continue
                seen.add(s)
                g.X = Xp; g.H = Hp
                cands.append((g, path + [mv[0]]))
        print(f'[{tag}] depth {d}: {len(cands)} new topologies to relax', flush=True)
        if not cands:
            break
        results = []
        t0 = time.time()
        for k, res in enumerate(pool.imap(evaluate, [(g.X, g.H, sorted(g.tets), r) for g, _ in cands]), 1):
            g, path = cands[k - 1]
            rec = dict(tag=tag, depth=d, path=path, cost=res['cost'], ok=res['ok'], flip=res.get('flip'),
                       faces=g.face_hist(), nf=g.nfaces_hist(), sig=str(g.signature()), t=res['t'])
            if res['ok'] and res['cost'] < base + 0.004:
                rec['tets'] = ser(g.tets); rec['X'] = res['Xn']; rec['H'] = res['Hn']
            out.write(json.dumps(rec) + '\n'); out.flush()
            results.append((g, res, path))
            if res['ok'] and res['cost'] < best_overall[0]:
                best_overall = (res['cost'], path)
                print(f'[{tag}] *** NEW BEST {res["cost"]:.6f} (start {base:.6f}) faces={g.face_hist()}', flush=True)
            if k % 10 == 0 or k == len(cands):
                okc = [x[1]['cost'] for x in results if x[1]['ok']]
                print(f'[{tag}] depth {d}: {k}/{len(cands)}  best this depth {min(okc) if okc else float("nan"):.6f}  '
                      f'(Δ vs start {min(okc)-base if okc else float("nan"):+.5f})  {time.time()-t0:.0f}s', flush=True)
        good = sorted([x for x in results if x[1]['ok']], key=lambda x: x[1]['cost'])
        level = good[:beam]
    print(f'[{tag}] done. start {base:.6f}, best {best_overall[0]:.6f} via {best_overall[1]}', flush=True)


if __name__ == '__main__':
    main()
