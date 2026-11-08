"""Periodic power (Laguerre) diagrams with equal-volume cells.

Lattice H: 3x3 array, rows are lattice vectors. Sites X: (N,3) cartesian.
Cell i = { x : |x-xi|^2 - wi <= |x - xj - nH|^2 - wj  for all j, n }.
"""
import numpy as np
from scipy.spatial import HalfspaceIntersection, ConvexHull
from scipy.optimize import linprog
import itertools

TOL = 1e-8


def shifts(H, nmax=2):
    ns = np.array(list(itertools.product(range(-nmax, nmax + 1), repeat=3)))
    return ns, ns @ H


def _poly_area(P, normal):
    c = P.mean(axis=0)
    n = normal / np.linalg.norm(normal)
    u = P[0] - c
    if np.linalg.norm(u) < 1e-14:
        u = P[1] - c
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    ang = np.arctan2((P - c) @ v, (P - c) @ u)
    order = np.argsort(ang)
    Q = P[order]
    s = np.zeros(3)
    for k in range(len(Q)):
        s += np.cross(Q[k] - c, Q[(k + 1) % len(Q)] - c)
    return 0.5 * abs(s @ n), order


def cell(i, X, w, H, NS, SV, rc):
    """Return dict with volume, area, faces for cell i (coordinates absolute)."""
    xi = X[i]
    D = (X[None, :, :] + SV[:, None, :]) - xi  # (S,N,3)
    dist = np.linalg.norm(D, axis=2)
    S_idx, J_idx = np.nonzero((dist < rc) & (dist > 1e-12))
    d = D[S_idx, J_idx]
    dd = dist[S_idx, J_idx]
    rhs = dd ** 2 - w[J_idx] + w[i]           # 2 d.y <= rhs
    A = 2 * d
    hs = np.hstack([A, -rhs[:, None]])
    if np.all(rhs > 1e-9):
        ip = np.zeros(3)
    else:  # Chebyshev center
        nrm = np.linalg.norm(A, axis=1)
        res = linprog(np.r_[0, 0, 0, -1], A_ub=np.hstack([A, nrm[:, None]]), b_ub=rhs,
                      bounds=[(None, None)] * 3 + [(0, None)])
        if res.status != 0 or res.x[3] < 1e-9:
            return None
        ip = res.x[:3]
    hsi = HalfspaceIntersection(hs, ip)
    V = hsi.intersections
    hull = ConvexHull(V)
    faces = []
    nrm = np.linalg.norm(A, axis=1)
    resid = np.abs(V @ A.T - rhs[None, :]) / nrm[None, :]  # (nv, nh)
    scale = max(1.0, np.abs(V).max())
    for k in range(len(rhs)):
        on = np.nonzero(resid[:, k] < 1e-7 * scale)[0]
        if len(on) < 3:
            continue
        P = V[on]
        # dedupe
        Pu = [P[0]]
        for p in P[1:]:
            if min(np.linalg.norm(p - q) for q in Pu) > 1e-9 * scale:
                Pu.append(p)
        P = np.array(Pu)
        if len(P) < 3:
            continue
        a, order = _poly_area(P, A[k])
        if a < 1e-12:
            continue
        faces.append(dict(j=int(J_idx[k]), n=NS[S_idx[k]].copy(), area=a, d=dd[k],
                          verts=P[order] + xi))
    rmax = np.linalg.norm(V - ip, axis=1).max()
    pc = V[hull.simplices]  # (k,3,3)
    o = V.mean(axis=0)
    tv = np.abs(np.einsum('ij,ij->i', pc[:, 0] - o, np.cross(pc[:, 1] - o, pc[:, 2] - o))) / 6
    cen = ((pc.sum(1) + o) / 4 * tv[:, None]).sum(0) / tv.sum() + xi
    return dict(volume=hull.volume, centroid=cen, area=sum(f['area'] for f in faces),
                faces=faces, rmax=rmax, hull_area=hull.area)


def diagram(X, w, H, nmax=2, rc=None):
    N = len(X)
    vol = abs(np.linalg.det(H))
    if rc is None:
        rc = 3.2 * (vol / N) ** (1 / 3)
    NS, SV = shifts(H, nmax)
    cells = []
    for i in range(N):
        c = cell(i, X, w, H, NS, SV, rc)
        if c is None:
            return None
        # validity: all neighbours that could matter are within rc
        if 2 * c['rmax'] + np.sqrt(max(0, w.max() - w.min())) > rc:
            return diagram(X, w, H, nmax + 1, rc * 1.5) if nmax < 4 else None
        cells.append(c)
    return cells


def equalize(X, H, w=None, iters=60, tol=1e-10, verbose=False):
    """Find weights making all cells equal volume (Newton)."""
    N = len(X)
    w = np.zeros(N) if w is None else w.copy()
    Vt = abs(np.linalg.det(H)) / N
    for it in range(iters):
        cells = diagram(X, w, H)
        if cells is None:
            raise RuntimeError('degenerate diagram')
        Vs = np.array([c['volume'] for c in cells])
        err = np.abs(Vs - Vt).max() / Vt
        if verbose:
            print(it, err)
        if err < tol:
            return w, cells
        J = np.zeros((N, N))
        for i, c in enumerate(cells):
            for f in c['faces']:
                if f['j'] == i:
                    continue
                g = f['area'] / (2 * f['d'])
                J[i, i] += g
                J[i, f['j']] -= g
        dw = np.linalg.lstsq(J, Vt - Vs, rcond=None)[0]
        dw -= dw.mean()
        # damped step
        step = 1.0
        while step > 1e-4:
            wn = w + step * dw
            cn = diagram(X, wn, H)
            if cn is not None:
                Vn = np.array([c['volume'] for c in cn])
                if np.abs(Vn - Vt).max() / Vt < err * (1 - 0.1 * step) or err < 1e-6:
                    break
            step *= 0.5
        w = wn
    return w, cells


def cost(cells, H):
    N = len(cells)
    Vt = abs(np.linalg.det(H)) / N
    return np.mean([c['area'] for c in cells]) / Vt ** (2 / 3)


def face_stats(cells):
    from collections import Counter
    fc = Counter()
    nf = []
    for c in cells:
        nf.append(len(c['faces']))
        for f in c['faces']:
            fc[len(f['verts'])] += 1
    return Counter(nf), fc
