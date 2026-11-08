"""Rigorous checks of the N56 vs sigma vs WP ranking.
usage: python rigor.py   -> rerun/rigor.jsonl
Checks
  conv   : graph-built relaxation at r = 1..4 (convergence / extrapolation)
  voro   : independent construction from an equal-volume power diagram (r = 1, 2)
  perturb: relax (r=2), randomly perturb every mesh point + the box, relax again
  plateau: at r=3, film angles along triple lines (ideal 120 deg) and volume error
"""
import sys, os, json, time, signal, warnings
import numpy as np
from multiprocessing import Pool
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
NAMES = {'WP': 'n20:Phys27', 'sigma': 'n20:Phys26', 'N56': 'n20:N56'}


def _alarm(s, f):
    raise TimeoutError()


def graph_mesh(rec, r):
    from tcpbuild import build_cells
    from relax import Mesh
    S = [set((j, tuple(n)) for j, n in s) for s in rec['S']]
    H = np.array(rec['H'])
    X, cells = build_cells(np.array(rec['F']), H, np.array(rec['Z']), S)
    return Mesh(X, H, cells, r=r)


def voro_mesh(rec, r):
    from power import equalize
    from relax import Mesh
    H = np.array(rec['H']); X = np.mod(np.array(rec['F']), 1) @ H
    w, cells = equalize(X, H, iters=80, tol=1e-12)
    nf = np.array([len(c['faces']) for c in cells])
    if not np.all(nf == np.array(rec['Z'])):
        raise ValueError('power diagram topology differs from graph')
    return Mesh(X, H, cells, r=r)


def relax(m, tol=1e-9):
    from relax import relax_box
    return float(relax_box(m, gtol=tol, tol=tol, mu0=300, maxit=30))


def flip_frac(m):
    """Area fraction of triangles whose orientation opposes their film's normal."""
    a, b, c = m._abs(m.P)
    cr = np.cross(b - a, c - a)
    ar = 0.5 * np.linalg.norm(cr, axis=1)
    bad = np.einsum('ij,ij->i', cr, m.face_normals[m.tface]) < 0
    return float(ar[bad].sum() / ar.sum())


def vol_err(m):
    A, V = m.energy(m.P, False)
    return float(np.abs(V / m.Vt - 1).max())


def plateau(m):
    """Angles between films meeting along mesh edges shared by 3 different faces."""
    a, b, c = m._abs(m.P)
    tri = np.stack([a, b, c], 1)
    edges = {}
    for t in range(len(m.tid)):
        for e in ((0, 1), (1, 2), (2, 0)):
            k1, k2 = m.tid[t, e[0]], m.tid[t, e[1]]
            s1, s2 = m.tsh[t, e[0]], m.tsh[t, e[1]]
            if (k1, tuple(s1)) > (k2, tuple(s2)):
                k1, k2, s1, s2 = k2, k1, s2, s1
            key = (k1, k2, tuple(np.round(s2 - s1).astype(int)))
            opp = 3 - e[0] - e[1]
            edges.setdefault(key, []).append((t, e, opp))
    ang = []
    for key, lst in edges.items():
        if len(lst) != 3 or len({m.tface[t] for t, _, _ in lst}) != 3:
            continue
        t0, e0, _ = lst[0]
        p = tri[t0, e0[0]]; q = tri[t0, e0[1]]
        u = (q - p) / np.linalg.norm(q - p)
        dirs = []
        for t, e, opp in lst:
            # bring this triangle's copy of the edge onto (p, q)
            pe = tri[t, e[0]]; qe = tri[t, e[1]]
            same = np.linalg.norm(pe - p) + np.linalg.norm(qe - q) < np.linalg.norm(pe - q) + np.linalg.norm(qe - p)
            base = pe if same else qe
            v = tri[t, opp] - base
            v = v - (v @ u) * u
            dirs.append(v / np.linalg.norm(v))
        for i in range(3):
            for j in range(i + 1, 3):
                ang.append(np.degrees(np.arccos(np.clip(dirs[i] @ dirs[j], -1, 1))))
    ang = np.array(ang)
    return dict(n=int(len(ang)), mean=float(ang.mean()), rms_dev=float(np.sqrt(((ang - 120) ** 2).mean())),
                max_dev=float(np.abs(ang - 120).max()))


def job(a):
    kind, label, r, seed, rec = a
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(20 * 3600)
    out = dict(kind=kind, label=label, r=r, seed=seed); t = time.time()
    try:
        if kind == 'conv':
            m = graph_mesh(rec, r); out['cost'] = relax(m); out['vol_err'] = vol_err(m); out['valid'] = float(m.valid()); out['flip_frac'] = flip_frac(m)
        elif kind == 'voro':
            m = voro_mesh(rec, r); out['cost'] = relax(m); out['vol_err'] = vol_err(m); out['valid'] = float(m.valid()); out['flip_frac'] = flip_frac(m)
        elif kind == 'perturb':
            m = graph_mesh(rec, r); c0 = relax(m)
            rng = np.random.default_rng(seed)
            s = 0.02 * m.Vt ** (1 / 3)
            m.P = m.P + s * rng.standard_normal(m.P.shape)
            Hn = m.H @ (np.eye(3) + 0.01 * rng.standard_normal((3, 3)))
            Hi = np.linalg.inv(m.H)
            Pf, Xf, tcsf = m.P @ Hi, m.X @ Hi, m.tcs @ Hi
            m.H = Hn; m.Hinv = np.linalg.inv(Hn); m.P = Pf @ Hn; m.X = Xf @ Hn
            m.tshc = m.tsh @ Hn; m.tcs = tcsf @ Hn
            m.ri = m.X[m.tcell[:, 0]]; m.rj = m.X[m.tcell[:, 1]] + m.tcs
            m.Vt = abs(np.linalg.det(Hn)) / m.N
            out['perturbed_cost'] = float(m.cost())
            out['cost0'] = c0; out['cost'] = relax(m); out['valid'] = float(m.valid()); out['flip_frac'] = flip_frac(m)
        elif kind == 'plateau':
            m = graph_mesh(rec, r); out['cost'] = relax(m); out['vol_err'] = vol_err(m)
            out['valid'] = float(m.valid()); out['flip_frac'] = flip_frac(m); out.update(plateau(m))
    except Exception as e:
        out['error'] = repr(e)[:200]
    signal.alarm(0); out['t'] = time.time() - t
    return out


if __name__ == '__main__':
    os.chdir(os.path.join(HERE, '..')); os.makedirs('rerun', exist_ok=True)
    allr = {x['name']: x for x in json.load(open('data/tcp_all.json'))}
    R = {k: allr[v] for k, v in NAMES.items()}
    jobs = []
    for k in R:
        jobs += [('conv', k, 4, 0, R[k])]
        jobs += [('voro', k, r, 0, R[k]) for r in (1, 2)]
        jobs += [('perturb', k, 2, s, R[k]) for s in (1, 2)]
        jobs += [('plateau', k, 3, 0, R[k])]
    jobs.sort(key=lambda j: -(j[4]['N'] * j[2] ** 2))
    with Pool(10, maxtasksperchild=1) as pool, open('rerun/rigor.jsonl', 'a') as f:
        for k, res in enumerate(pool.imap_unordered(job, jobs), 1):
            f.write(json.dumps(res) + '\n'); f.flush()
            print(f"{k}/{len(jobs)} {res['kind']:8s} {res['label']:6s} r={res['r']} "
                  f"{res.get('cost', res.get('error'))} ({res['t']:.0f}s)", flush=True)
    print('[완료] 엄밀성 검증이 끝났습니다.', flush=True)
