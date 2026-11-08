"""Curved relaxation of a periodic foam: minimise total film area subject to
equal cell volumes. Faces are triangulated; all mesh points move freely.
Augmented Lagrangian + L-BFGS. Self-written stand-in for Surface Evolver."""
import numpy as np
from scipy.optimize import minimize
import itertools


class Mesh:
    def __init__(self, X, H, cells, r=2):
        self.X = np.asarray(X, float)
        self.H = np.asarray(H, float)
        self.Hinv = np.linalg.inv(self.H)
        self.N = len(X)
        self.keys = {}
        self.P = []
        tid, tsh, tcell, tcs, tface = [], [], [], [], []
        self.face_normals = []
        nface = 0
        for i, c in enumerate(cells):
            for f in c['faces']:
                j, n = f['j'], tuple(int(v) for v in f['n'])
                if j < i or (j == i and n < (0, 0, 0)):
                    continue
                V = f['verts']
                cen = V.mean(axis=0)
                d = self.X[j] + np.array(n) @ self.H - self.X[i]
                nrm = np.cross(V[1] - cen, V[0] - cen) * 0 + np.cross(V[0] - cen, V[1] - cen)
                if nrm @ d < 0:
                    V = V[::-1]
                fn = d / np.linalg.norm(d)
                self.face_normals.append(fn)
                m = len(V)
                for k in range(m):
                    A, B = V[k], V[(k + 1) % m]
                    # sub-triangulate fan triangle (cen, A, B)
                    def pt(a, b):
                        if a + b == r:  # on boundary edge: symmetric formula
                            t = b / r
                            p = (1 - t) * A + t * B if a >= b else A + t * (B - A)
                            return self._get(p)
                        return self._get(cen + (a / r) * (A - cen) + (b / r) * (B - cen))
                    grid = {}
                    for a in range(r + 1):
                        for b in range(r + 1 - a):
                            grid[a, b] = pt(a, b)
                    for a in range(r):
                        for b in range(r - a):
                            tris = [(grid[a, b], grid[a + 1, b], grid[a, b + 1])]
                            if a + b + 2 <= r:
                                tris.append((grid[a + 1, b], grid[a + 1, b + 1], grid[a, b + 1]))
                            for tri in tris:
                                tid.append([t[0] for t in tri])
                                tsh.append([t[1] for t in tri])
                                tcell.append([i, j])
                                tcs.append(n)
                                tface.append(nface)
                nface += 1
        self.P = np.array(self.P)
        self.tid = np.array(tid)
        self.tsh = np.array(tsh, float)  # (T,3,3) integer lattice shifts
        self.tshc = self.tsh @ self.H    # cartesian shifts
        self.tcell = np.array(tcell)
        self.tcs = np.array(tcs, float) @ self.H
        self.tface = np.array(tface)
        self.face_normals = np.array(self.face_normals)
        self.nfaces = nface
        # orientation fix per triangle relative to face normal
        a, b, c = self._abs(self.P)
        cr = np.cross(b - a, c - a)
        flip = np.einsum('ij,ij->i', cr, self.face_normals[self.tface]) < 0
        self.tid[flip] = self.tid[flip][:, [0, 2, 1]]
        self.tsh[flip] = self.tsh[flip][:, [0, 2, 1]]
        self.tshc = self.tsh @ self.H
        self.ri = self.X[self.tcell[:, 0]]
        self.rj = self.X[self.tcell[:, 1]] + self.tcs
        self.Vt = abs(np.linalg.det(self.H)) / self.N

    def _get(self, p, tol=1e-5, g=1e-3):
        """Merge points equal modulo the lattice within fractional tol."""
        f = p @ self.Hinv
        fw = np.mod(f, 1.0)
        b = np.floor(fw / g).astype(np.int64)
        nb = int(round(1 / g))
        found = None
        for d in itertools.product((-1, 0, 1), repeat=3):
            key = tuple(int(v) for v in np.mod(b + d, nb))
            for k in self.keys.get(key, ()):
                q = np.asarray(self.P[k]) @ self.Hinv
                if np.abs(((fw - q + .5) % 1) - .5).max() < tol:
                    found = k; break
            if found is not None:
                break
        if found is None:
            found = len(self.P)
            self.P.append(fw @ self.H)
            self.keys.setdefault(tuple(int(v) for v in b % nb), []).append(found)
        canon = np.asarray(self.P[found]) @ self.Hinv
        sh = np.round(f - canon).astype(int)
        return found, tuple(sh)

    def _abs(self, P):
        Q = P[self.tid] + self.tshc
        return Q[:, 0], Q[:, 1], Q[:, 2]

    def energy(self, P, grad=True):
        a, b, c = self._abs(P)
        cr = np.cross(b - a, c - a)
        nn = np.linalg.norm(cr, axis=1)
        area = 0.5 * nn.sum()
        # volumes
        ai, bi, ci = a - self.ri, b - self.ri, c - self.ri
        aj, bj, cj = a - self.rj, b - self.rj, c - self.rj
        di = np.einsum('ij,ij->i', ai, np.cross(bi, ci)) / 6
        dj = np.einsum('ij,ij->i', aj, np.cross(bj, cj)) / 6
        V = np.zeros(self.N)
        np.add.at(V, self.tcell[:, 0], di)
        np.add.at(V, self.tcell[:, 1], -dj)
        if not grad:
            return area, V
        nh = cr / np.maximum(nn, 1e-300)[:, None]
        gA = 0.5 * np.stack([np.cross(c - b, nh), np.cross(a - c, nh), np.cross(b - a, nh)], 1) * -1
        # check sign: dA/da = 0.5 * nh x (b - c)? handled by unit test
        gVi = np.stack([np.cross(bi, ci), np.cross(ci, ai), np.cross(ai, bi)], 1) / 6
        gVj = -np.stack([np.cross(bj, cj), np.cross(cj, aj), np.cross(aj, bj)], 1) / 6
        return area, V, gA, gVi, gVj

    def scatter(self, gT):
        g = np.zeros_like(self.P)
        np.add.at(g, self.tid.ravel(), gT.reshape(-1, 3))
        return g

    def relax(self, maxit=40, mu=None, tol=1e-10, verbose=False):
        P = self.P.copy()
        lam = np.zeros(self.N)
        mu = 10.0 / self.Vt ** (4 / 3) if mu is None else mu
        for outer in range(maxit):
            def f(x):
                Pm = x.reshape(-1, 3)
                A, V, gA, gVi, gVj = self.energy(Pm)
                dv = V - self.Vt
                coef = lam + mu * dv
                E = A + lam @ dv + 0.5 * mu * dv @ dv
                gT = gA + coef[self.tcell[:, 0]][:, None, None] * gVi + coef[self.tcell[:, 1]][:, None, None] * gVj
                return E, self.scatter(gT).ravel()
            res = minimize(f, P.ravel(), jac=True, method='L-BFGS-B',
                           options=dict(maxiter=3000, maxcor=30, ftol=1e-15, gtol=1e-10))
            P = res.x.reshape(-1, 3)
            A, V = self.energy(P, grad=False)
            dv = V - self.Vt
            lam = lam + mu * dv
            verr = np.abs(dv).max() / self.Vt
            if verbose:
                print(outer, A, verr, res.nit)
            if verr < tol and outer > 2:
                break
            if outer > 3 and verr > 1e-6:
                mu *= 3
        self.P = P
        return self.cost(P)

    def cost(self, P=None):
        P = self.P if P is None else P
        A, V = self.energy(P, grad=False)
        return 2 * A / self.N / self.Vt ** (2 / 3)

    def valid(self, P=None):
        """No triangle flipped relative to its flat face normal."""
        P = self.P if P is None else P
        a, b, c = self._abs(P)
        cr = np.cross(b - a, c - a)
        cosn = np.einsum('ij,ij->i', cr, self.face_normals[self.tface]) / np.linalg.norm(cr, axis=1)
        return cosn.min()


def relax_box(m, maxit=40, tol=1e-10, verbose=False, free_box=True, gtol=1e-11, mu0=10.0, maxiter=4000):
    """Relax mesh m including the lattice H (all 9 components). Dimensionless:
    minimise 2A/(N s^2) with V_i = s^3, s^3 = det(H)/N."""
    N = m.N
    Pf = m.P @ m.Hinv
    Xf = m.X @ m.Hinv
    tshf = m.tsh                      # (T,3,3) fractional shifts
    rjf_shift = m.tcs @ m.Hinv       # fractional cell shifts (T,3)
    H0 = m.H.copy()
    nP = Pf.size

    def unpack(x):
        Pf_ = x[:nP].reshape(-1, 3)
        H_ = x[nP:].reshape(3, 3) if free_box else H0
        return Pf_, H_

    def evaluate(x, lam, mu):
        Pf_, H_ = unpack(x)
        F = Pf_[m.tid] + tshf                     # (T,3,3) fractional
        Q = F @ H_
        a, b, c = Q[:, 0], Q[:, 1], Q[:, 2]
        ri = Xf[m.tcell[:, 0]] @ H_
        rj = (Xf[m.tcell[:, 1]] + rjf_shift) @ H_
        cr = np.cross(b - a, c - a)
        nn = np.linalg.norm(cr, axis=1)
        A = 0.5 * nn.sum()
        nh = cr / nn[:, None]
        gA = -0.5 * np.stack([np.cross(c - b, nh), np.cross(a - c, nh), np.cross(b - a, nh)], 1)
        ai, bi, ci = a - ri, b - ri, c - ri
        aj, bj, cj = a - rj, b - rj, c - rj
        di = np.einsum('ij,ij->i', ai, np.cross(bi, ci)) / 6
        dj = np.einsum('ij,ij->i', aj, np.cross(bj, cj)) / 6
        V = np.zeros(N)
        np.add.at(V, m.tcell[:, 0], di)
        np.add.at(V, m.tcell[:, 1], -dj)
        gVi = np.stack([np.cross(bi, ci), np.cross(ci, ai), np.cross(ai, bi)], 1) / 6
        gVj = -np.stack([np.cross(bj, cj), np.cross(cj, aj), np.cross(aj, bj)], 1) / 6
        detH = np.linalg.det(H_)
        s3 = abs(detH) / N
        s2 = s3 ** (2 / 3)
        E0 = 2 * A / N / s2
        g = V / s3 - 1
        coef = lam + mu * g
        E = E0 + lam @ g + 0.5 * mu * g @ g
        # dE/dq (cartesian per triangle vertex)
        ci_ = coef[m.tcell[:, 0]] / s3
        cj_ = coef[m.tcell[:, 1]] / s3
        GQ = (2 / N / s2) * gA + ci_[:, None, None] * gVi + cj_[:, None, None] * gVj
        Gri = -ci_[:, None] * gVi.sum(1)
        Grj = -cj_[:, None] * gVj.sum(1)
        gPf = np.zeros_like(Pf_)
        np.add.at(gPf, m.tid.ravel(), (GQ @ H_.T).reshape(-1, 3))
        out = gPf.ravel()
        if free_box:
            gH = np.einsum('tki,tkj->ij', F, GQ)
            gH += np.einsum('ti,tj->ij', Xf[m.tcell[:, 0]], Gri)
            gH += np.einsum('ti,tj->ij', Xf[m.tcell[:, 1]] + rjf_shift, Grj)
            # derivative through s3 (det H): dE/dlogdet * Hinv^T
            dE_dlogs3 = -(2 / 3) * E0 - (coef * (V / s3)).sum()
            gH += dE_dlogs3 * np.linalg.inv(H_).T
            out = np.concatenate([out, gH.ravel()])
        return E, out, E0, g

    x = np.concatenate([Pf.ravel(), H0.ravel()]) if free_box else Pf.ravel()
    lam = np.zeros(N)
    mu = mu0
    for outer in range(maxit):
        res = minimize(lambda y: evaluate(y, lam, mu)[:2], x, jac=True, method='L-BFGS-B',
                       options=dict(maxiter=maxiter, maxcor=30, ftol=1e-15, gtol=gtol))
        x = res.x
        _, _, E0, g = evaluate(x, lam, mu)
        lam = lam + mu * g
        verr = np.abs(g).max()
        if verbose:
            print(outer, E0, verr, res.nit)
        if verr < tol and outer > 2:
            break
        if outer > 3 and verr > 1e-6:
            mu *= 3
    Pf_, H_ = unpack(x)
    m.H = H_.copy(); m.Hinv = np.linalg.inv(H_)
    m.P = Pf_ @ H_
    m.tshc = m.tsh @ H_
    m.tcs = rjf_shift @ H_
    m.X = Xf @ H_
    m.ri = m.X[m.tcell[:, 0]]
    m.rj = m.X[m.tcell[:, 1]] + m.tcs
    m.Vt = abs(np.linalg.det(H_)) / N
    return E0
