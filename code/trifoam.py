"""Periodic foam represented by its dual triangulation.
Cells = vertices (i, shift); foam faces = edges; foam edges = triangles; foam vertices = tetrahedra.
Foam T1 moves = Pachner 2-3 (edge -> new triangular face) and 3-2 (triangular face -> edge).
Topologies are not restricted to TCP or to power diagrams.

A vertex is (i, (n1,n2,n3)). A tetrahedron is stored once per lattice-translation class
in canonical form (see canon)."""
import itertools
import numpy as np
from collections import defaultdict, Counter

Z3 = (0, 0, 0)


def _sub(v, s):
    return (v[0], (v[1][0] - s[0], v[1][1] - s[1], v[1][2] - s[2]))


def _add(v, s):
    return (v[0], (v[1][0] + s[0], v[1][1] + s[1], v[1][2] + s[2]))


def canon(simplex):
    """Canonical representative of the translation class of a simplex (tuple of vertices).
    Returns (canonical tuple, shift t) with canonical = simplex - t."""
    best = None
    for v in simplex:
        t = v[1]
        c = tuple(sorted(_sub(w, t) for w in simplex))
        if best is None or c < best[0]:
            best = (c, t)
    return best


class TriFoam:
    def __init__(self, X, H, tets):
        """X: (N,3) cartesian cell reference points; H: rows = lattice vectors;
        tets: iterable of 4-tuples of vertices (any translates)."""
        self.X = np.asarray(X, float)
        self.H = np.asarray(H, float)
        self.N = len(self.X)
        self.tets = set(canon(t)[0] for t in tets)
        self._index()

    # ------------------------------------------------------------------ indices
    def _index(self):
        self.tri = defaultdict(list)   # canonical triangle -> list of apex vertices (canonical frame)
        self.edge = defaultdict(list)  # canonical edge -> list of (k, l) pairs (canonical frame)
        for T in self.tets:
            for o in range(4):
                tri = tuple(T[q] for q in range(4) if q != o)
                c, t = canon(tri)
                self.tri[c].append(_sub(T[o], t))
            for a, b in itertools.combinations(range(4), 2):
                e = (T[a], T[b])
                c, t = canon(e)
                rest = [_sub(T[q], t) for q in range(4) if q not in (a, b)]
                self.edge[c].append(tuple(rest))

    def check(self):
        """Return None if a valid closed 3-manifold triangulation (every triangle in 2 tets,
        every edge link a single cycle, no repeated vertex); else an error string."""
        for T in self.tets:
            if len(set(T)) < 4:
                return 'degenerate tet'
        for c, ap in self.tri.items():
            if len(ap) != 2:
                return f'triangle with {len(ap)} tets'
            if ap[0] == ap[1]:
                return 'two tets identical'
        for c, pairs in self.edge.items():
            if len(pairs) < 3:
                return 'edge degree < 3'
            if self._cycle(pairs) is None:
                return 'edge link not a cycle'
        # vertex links: Euler characteristic of each vertex link must be 2 (sphere)
        for i in range(self.N):
            V, E, F = set(), set(), 0
            for T in self.tets:
                for v in T:
                    if v[0] != i:
                        continue
                    L = [_sub(w, v[1]) for w in T if w != v]
                    F += 1
                    for w in L:
                        V.add(w)
                    for a, b in itertools.combinations(L, 2):
                        E.add(tuple(sorted((a, b))))
            if len(V) - len(E) + F != 2:
                return f'vertex link of cell {i} not a sphere (chi={len(V)-len(E)+F})'
        return None

    @staticmethod
    def _cycle(pairs):
        adj = defaultdict(list)
        for k, l in pairs:
            adj[k].append(l); adj[l].append(k)
        if any(len(v) != 2 for v in adj.values()):
            return None
        start = pairs[0][0]
        ring = [start]; prev = None; cur = start
        while True:
            nx = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
            if nx == start:
                break
            ring.append(nx); prev, cur = cur, nx
            if len(ring) > len(pairs):
                return None
        return ring if len(ring) == len(pairs) else None

    # ------------------------------------------------------------------ moves
    def moves(self):
        """All candidate moves: ('23', triangle) and ('32', edge)."""
        out = []
        for c, ap in self.tri.items():
            d, e = ap
            ce, _ = canon((d, e))
            if ce not in self.edge and d != e:
                out.append(('23', c))
        for c, pairs in self.edge.items():
            if len(pairs) == 3:
                ring = self._cycle(pairs)
                if ring is not None:
                    ct, _ = canon(tuple(ring))
                    if ct not in self.tri:
                        out.append(('32', c))
        return out

    def apply(self, move):
        kind, c = move
        tets = set(self.tets)
        if kind == '23':
            a, b, cc = c
            d, e = self.tri[c]
            tets.discard(canon((a, b, cc, d))[0]); tets.discard(canon((a, b, cc, e))[0])
            new = [(a, b, d, e), (b, cc, d, e), (cc, a, d, e)]
        else:
            d, e = c
            ring = self._cycle(self.edge[c])
            for k in range(3):
                tets.discard(canon((d, e, ring[k], ring[(k + 1) % 3]))[0])
            new = [tuple(ring) + (d,), tuple(ring) + (e,)]
        for t in new:
            tets.add(canon(t)[0])
        f = TriFoam(self.X, self.H, [])
        f.tets = tets
        f._index()
        return f

    # ------------------------------------------------------------------ geometry
    def pos(self, v):
        return self.X[v[0]] + np.array(v[1]) @ self.H

    def cells(self):
        """Cells in the format used by relax.Mesh."""
        per = defaultdict(lambda: defaultdict(list))  # cell i -> neighbour (j,n) -> list of pairs
        for T in self.tets:
            for v in T:
                Tl = [_sub(w, v[1]) for w in T]
                vi = _sub(v, v[1])
                others = [w for w in Tl if w != vi]
                for w in others:
                    rest = tuple(x for x in others if x != w)
                    per[vi[0]][w].append(rest)
        out = []
        for i in range(self.N):
            faces = []
            Xi = self.X[i]
            for w, pairs in per[i].items():
                ring = self._cycle(pairs)
                if ring is None:
                    raise ValueError('bad edge link')
                pw = self.pos(w)
                V = np.array([(Xi + pw + self.pos(ring[a]) + self.pos(ring[(a + 1) % len(ring)])) / 4
                              for a in range(len(ring))])
                faces.append(dict(j=w[0], n=np.array(w[1]), verts=V, area=0.0))
            out.append(dict(faces=faces))
        return out

    # ------------------------------------------------------------------ descriptors
    def signature(self):
        """Topological fingerprint: per-cell (face count, sorted polygon sizes) multiset,
        normalised by gcd. Not a full isomorphism test, but separates most cases."""
        deg = {c: len(p) for c, p in self.edge.items()}
        cellsig = []
        for i in range(self.N):
            sizes = []
            for c, d in deg.items():
                for v in c:
                    if v[0] == i:
                        sizes.append(d)
            # each self-edge (i,i+n) counted twice -> fine, consistent
            cellsig.append(tuple(sorted(Counter(sizes).items())))
        cnt = Counter(cellsig)
        from math import gcd
        from functools import reduce
        g = reduce(gcd, cnt.values())
        return tuple(sorted((k, v // g) for k, v in cnt.items()))

    def face_hist(self):
        return dict(sorted(Counter(len(p) for p in self.edge.values()).items()))

    def nfaces_hist(self):
        nf = Counter()
        for c in self.edge:
            for v in c:
                nf[v[0]] += 1
        return dict(sorted(Counter(nf[i] for i in range(self.N)).items()))

    def supercell(self, reps):
        reps = tuple(reps)
        cells = list(itertools.product(*[range(r) for r in reps]))
        idx = {(i, c): k for k, (c, i) in enumerate(itertools.product(cells, range(self.N)))}
        X = np.array([self.X[i] + np.array(c) @ self.H for c in cells for i in range(self.N)])
        H = np.diag(reps).astype(float) @ self.H

        def lift(v, c0):
            g = np.array(v[1]) + np.array(c0)
            q = np.floor_divide(g, reps); r = g - q * np.array(reps)
            return (idx[(v[0], tuple(int(x) for x in r))], tuple(int(x) for x in q))
        tets = [tuple(lift(v, c0) for v in T) for T in self.tets for c0 in cells]
        return TriFoam(X, H, tets)


def from_tcp(rec):
    """Build from a tcp_all.json record (neighbour sets S): tets = 4-cliques."""
    H = np.array(rec['H'], float)
    X = np.mod(np.array(rec['F'], float), 1) @ H
    S = [set((j, tuple(n)) for j, n in s) for s in rec['S']]

    def adj(a, b):
        return _sub(b, a[1]) in S[a[0]] if a[0] < len(S) else False
    tets = set()
    for i in range(len(S)):
        vi = (i, Z3)
        nb = sorted(S[i])
        for a, b, c in itertools.combinations(nb, 3):
            if adj(a, b) and adj(b, c) and adj(a, c):
                tets.add(canon((vi, a, b, c))[0])
    return TriFoam(X, H, tets)


def relaxed_centroids(m):
    """Volume centroids of cells from a relaxed relax.Mesh (cartesian, cell i's frame)."""
    a, b, c = m._abs(m.P)
    N = m.N
    out = np.zeros((N, 3)); vol = np.zeros(N)
    for side, ref in ((0, m.ri), (1, m.rj)):
        ai, bi, ci = a - ref, b - ref, c - ref
        v = np.einsum('ij,ij->i', ai, np.cross(bi, ci)) / 6 * (1 if side == 0 else -1)
        cen = ref + (ai + bi + ci) / 4
        if side == 1:
            cen = cen - m.tcs  # bring back to frame of cell j at zero shift
        idx = m.tcell[:, side]
        np.add.at(vol, idx, v)
        np.add.at(out, idx, v[:, None] * cen)
    return out / vol[:, None]


def from_delaunay(F, H):
    """Dual triangulation = periodic Delaunay of the points (foam = Voronoi topology)."""
    from scipy.spatial import Delaunay
    F = np.mod(np.asarray(F, float), 1); H = np.asarray(H, float); N = len(F)
    sh = np.array(list(itertools.product((-1, 0, 1), repeat=3)))
    P = np.concatenate([(F + s) @ H for s in sh]); lab = [(i, tuple(int(x) for x in s)) for s in sh for i in range(N)]
    D = Delaunay(P)
    tets = set()
    for simp in D.simplices:
        vs = [lab[k] for k in simp]
        if any(v[1] == Z3 for v in vs):
            tets.add(canon(tuple(vs))[0])
    return TriFoam(F @ H, H, tets)
