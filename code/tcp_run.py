"""Relax every TCP structure in tcp_all.json (fixed topology, equal volumes).
usage: python tcp_run.py [workers=10]
Phase 1: r=1 for all.  Phase 2: r=2 for those with r1 < 5.300.
Results: rerun/tcp_r1.jsonl, rerun/tcp_r2.jsonl  (WP: r1 5.29043, r2 5.28884)"""
import sys, os, json, time, signal, warnings
import numpy as np
from multiprocessing import Pool

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
LIMIT = {1: 4 * 3600, 2: 10 * 3600}


class Timeout(Exception):
    pass


def _alarm(s, f):
    raise Timeout()


def job(args):
    rec, r = args
    from tcpbuild import build_cells
    from relax import Mesh, relax_box
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(LIMIT[r])
    t = time.time()
    out = dict(name=rec['name'], N=rec['N'], meanZ=rec['meanZ'], r=r)
    try:
        F = np.array(rec['F']); H = np.array(rec['H']); Z = np.array(rec['Z'])
        S = [set((j, tuple(n)) for j, n in s) for s in rec['S']]
        X, cells = build_cells(F, H, Z, S)
        m = Mesh(X, H, cells, r=r)
        c = relax_box(m, gtol=1e-8 if r == 1 else 1e-9, tol=1e-8 if r == 1 else 1e-9, mu0=300, maxit=15)
        v = float(m.valid())
        out.update(cost=float(c), valid=v, ok=bool(np.isfinite(c) and v > 0), tris=int(len(m.tid)),
                   Hrel=m.H.tolist())
    except Timeout:
        out.update(error='timeout')
    except Exception as e:
        out.update(error=repr(e)[:200])
    finally:
        signal.alarm(0)
    out['t'] = time.time() - t
    return out


def done(path):
    if not os.path.exists(path):
        return {}
    return {d['name']: d for d in map(json.loads, open(path))}


def run_phase(recs, r, nw, path):
    prev = done(path)
    todo = [x for x in recs if x['name'] not in prev]
    todo.sort(key=lambda x: -x['N'])  # biggest first for load balance
    print(f'[r={r}] {len(todo)} structures to relax ({len(prev)} already done)', flush=True)
    with Pool(nw, maxtasksperchild=1) as pool, open(path, 'a') as f:
        for k, res in enumerate(pool.imap_unordered(job, [(x, r) for x in todo]), 1):
            f.write(json.dumps(res) + '\n'); f.flush()
            tag = f"{res['cost']:.5f} valid={res['valid']:.2f}" if 'cost' in res else res.get('error')
            print(f"[r={r}] {k}/{len(todo)} {res['name']:24s} N={res['N']:4d} {tag} ({res['t']:.0f}s)", flush=True)
    return done(path)


if __name__ == '__main__':
    for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        os.environ[v] = '1'
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    os.chdir(os.path.join(HERE, '..')); os.makedirs('rerun', exist_ok=True)
    recs = json.load(open('data/tcp_all.json'))
    r1 = run_phase(recs, 1, nw, 'rerun/tcp_r1.jsonl')
    sel = [x for x in recs if x['name'] in r1 and r1[x['name']].get('ok') and r1[x['name']]['cost'] < 5.300]
    run_phase(sel, 2, nw, 'rerun/tcp_r2.jsonl')
    print('[완료] 모든 단계가 끝났습니다.', flush=True)
