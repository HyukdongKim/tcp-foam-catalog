"""(a) r=3 check of the top TCP structures; (b) retry failed relaxations with a gentler constraint ramp.
usage: python tcp_verify.py   -> rerun/tcp_verify.jsonl"""
import sys, os, json, time, signal, warnings
import numpy as np
from multiprocessing import Pool
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
TOP = ['n20:Phys27', 'n20:N56', 'n20:Phys26', 'n20:Phys25', 'n20:E11']
RETRY = ['n21:Struct19', 'n21:Struct20', 'n21:Struct37', 'n20:N60', 'n21:Struct10',
         'n20:N13', 'n20:N15', 'n21:Struct27', 'n20:E11']


def _alarm(s, f):
    raise TimeoutError()


def job(a):
    rec, r, mu0, tag = a
    from tcpbuild import build_cells
    from relax import Mesh, relax_box
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(12 * 3600)
    out = dict(name=rec['name'], N=rec['N'], meanZ=rec['meanZ'], r=r, mu0=mu0, task=tag); t = time.time()
    try:
        S = [set((j, tuple(n)) for j, n in s) for s in rec['S']]
        X, cells = build_cells(np.array(rec['F']), np.array(rec['H']), np.array(rec['Z']), S)
        m = Mesh(X, np.array(rec['H']), cells, r=r)
        c = relax_box(m, gtol=1e-9, tol=1e-9, mu0=mu0, maxit=30)
        out.update(cost=float(c), valid=float(m.valid()), tris=int(len(m.tid)))
    except Exception as e:
        out.update(error=repr(e)[:200])
    signal.alarm(0); out['t'] = time.time() - t
    return out


if __name__ == '__main__':
    os.chdir(os.path.join(HERE, '..')); os.makedirs('rerun', exist_ok=True)
    recs = {x['name']: x for x in json.load(open('data/tcp_all.json'))}
    jobs = [(recs[n], 3, 300.0, 'r3') for n in TOP]
    jobs += [(recs[n], r, mu, 'retry') for n in RETRY for r, mu in ((1, 1.0), (2, 1.0))]
    jobs.sort(key=lambda j: -(j[0]['N'] * j[1] ** 2))
    with Pool(10, maxtasksperchild=1) as pool, open('rerun/tcp_verify.jsonl', 'a') as f:
        for k, res in enumerate(pool.imap_unordered(job, jobs), 1):
            f.write(json.dumps(res) + '\n'); f.flush()
            tag = f"{res['cost']:.5f} valid={res['valid']:.2f}" if 'cost' in res else res.get('error')
            print(f"{k}/{len(jobs)} [{res['task']} r={res['r']}] {res['name']:16s} {tag} ({res['t']:.0f}s)", flush=True)
    print('[완료] 검증이 끝났습니다.', flush=True)
