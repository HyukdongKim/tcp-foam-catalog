"""Build a TCP foam directly from its neighbour graph (topology fixed).
Cells = atoms; face(i,j) = ring of common neighbours of i and j;
foam vertex = centroid of the tetrahedron (i, j, k_m, k_m+1)."""
import itertools
import numpy as np


def neighbours(F, H, Z):
    X = F @ H
    N = len(F)
    ns = np.array(list(itertools.product((-1, 0, 1), repeat=3)))
    nb = []
    for i in range(N):
        cand = []
        for j in range(N):
            for n in ns:
                if j == i and not n.any():
                    continue
                d = np.linalg.norm(X[j] + n @ H - X[i])
                cand.append((d, j, tuple(n)))
        cand.sort()
        nb.append([(j, n) for _, j, n in cand[:Z[i]]])
    # symmetry check
    sets = [set(l) for l in nb]
    ok = all((i, tuple(-np.array(n))) in sets[j] for i in range(N) for j, n in nb[i])
    return nb, sets, ok


def build_cells(F, H, Z, sets=None):
    if sets is None:
        F = np.mod(F, 1)
        nb, sets, ok = neighbours(F, H, Z)
        if not ok:
            raise ValueError('neighbour graph not symmetric')
    else:
        nb = [sorted(s) for s in sets]
        if any(len(s) != z for s, z in zip(sets, Z)):
            raise ValueError('edge sets do not match coordination')
    X = F @ H
    N = len(F)

    def pos(j, n):
        return X[j] + np.array(n) @ H

    def adj(a, b):  # a,b = (atom, shift) absolute; adjacent?
        (ja, na), (jb, nb_) = a, b
        rel = tuple(np.array(nb_) - np.array(na))
        return (jb, rel) in sets[ja]

    cells = []
    for i in range(N):
        faces = []
        for j, n in nb[i]:
            pj = pos(j, n)
            A = (i, (0, 0, 0)); B = (j, n)
            ring = [(k, m) for k, m in nb[i] if (k, m) != B and adj(B, (k, m))]
            if len(ring) < 3:
                raise ValueError('bad ring')
            u = pj - X[i]; u /= np.linalg.norm(u)
            mid = 0.5 * (X[i] + pj)
            e1 = np.cross(u, [1, 0, 0] if abs(u[0]) < .9 else [0, 1, 0]); e1 /= np.linalg.norm(e1)
            e2 = np.cross(u, e1)
            P = np.array([pos(*r) - mid for r in ring])
            order = np.argsort(np.arctan2(P @ e2, P @ e1))
            ring = [ring[o] for o in order]
            # consecutive ring members must be adjacent (TCP)
            for a in range(len(ring)):
                if not adj(ring[a], ring[(a + 1) % len(ring)]):
                    raise ValueError('ring not a cycle')
            V = np.array([(X[i] + pj + pos(*ring[a]) + pos(*ring[(a + 1) % len(ring)])) / 4
                          for a in range(len(ring))])
            faces.append(dict(j=j, n=np.array(n), verts=V, area=0.0,
                              d=np.linalg.norm(pj - X[i])))
        cells.append(dict(faces=faces))
    return X, cells
