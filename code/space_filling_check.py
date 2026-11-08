"""Space-filling check: every point of the box must lie inside exactly one cell (winding number)."""
import sys, os, json, numpy as np, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rigor
from relax import relax_box
a = {x['name']: x for x in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'tcp_all.json')))}


def solid(P, T):  # P (M,3) points, T (K,3,3) triangles -> sum of signed solid angles / 4pi per point
    w = np.zeros(len(P))
    for s in range(0, len(T), 400):
        A = T[None, s:s+400, 0, :] - P[:, None, :]; B = T[None, s:s+400, 1, :] - P[:, None, :]; C = T[None, s:s+400, 2, :] - P[:, None, :]
        a_, b_, c_ = (np.linalg.norm(v, axis=2) for v in (A, B, C))
        num = np.einsum('mkj,mkj->mk', A, np.cross(B, C))
        den = a_ * b_ * c_ + np.einsum('mkj,mkj->mk', A, B) * c_ + np.einsum('mkj,mkj->mk', A, C) * b_ + np.einsum('mkj,mkj->mk', B, C) * a_
        w += (2 * np.arctan2(num, den)).sum(1)
    return w / (4 * np.pi)


for n in sys.argv[1:]:
    m = rigor.graph_mesh(a[n], 1); c = float(relax_box(m, gtol=1e-9, tol=1e-9, mu0=300, maxit=30))
    A, B, C = m._abs(m.P); T = np.stack([A, B, C], 1)
    rng = np.random.default_rng(0); P = rng.random((1500, 3)) @ m.H
    shifts = np.array(list(itertools.product((-1, 0, 1), repeat=3))) @ m.H
    count = np.zeros(len(P))
    for i in range(m.N):
        mi = m.tcell[:, 0] == i; mj = m.tcell[:, 1] == i
        Ti = np.concatenate([T[mi], (T[mj] - m.tcs[mj][:, None, :])[:, ::-1]])
        for sh in shifts:
            count += solid(P - sh, Ti)
    k = np.round(count).astype(int)
    print(f"{n:14s} cost={c:.5f}  points in exactly one cell: {np.mean(k == 1)*100:.1f}%   counts seen: {dict(zip(*np.unique(k, return_counts=True)))}", flush=True)
