"""Write a Surface Evolver TORUS_FILLED datafile for a TCP foam (flat polyhedral start).
usage: python export_fe.py tcp_all.json NAME Hbox.json out.fe"""
import sys, json
import numpy as np
from tcpbuild import build_cells

src, name, hfile, outp = sys.argv[1:5]
FIXV = '--fixv' in sys.argv
rec = {x['name']: x for x in json.load(open(src))}[name]
F = np.array(rec['F']); H = np.array(rec['H']); Z = np.array(rec['Z'])
S = [set((j, tuple(n)) for j, n in s) for s in rec['S']]
if '--voro' in sys.argv:
    # flat equal-volume power diagram: every cell starts at exactly volume 1
    from power import equalize
    X = np.mod(F, 1) @ H
    w, cells = equalize(X, H, iters=80, tol=1e-13)
    if not all(len(c['faces']) == z for c, z in zip(cells, Z)):
        raise ValueError('power diagram topology differs from the TCP graph')
else:
    X, cells = build_cells(F, H, Z, S)
Hbox = np.array(json.load(open(hfile))[name]) if hfile != '-' else H
Hinv = np.linalg.inv(H)
N = len(F)

verts = []          # canonical fractional coords


def vid(p):
    f = p @ Hinv
    fw = np.mod(f, 1.0)
    for k, q in enumerate(verts):
        d = ((fw - q + .5) % 1) - .5
        if np.abs(d).max() < 1e-5:
            return k, tuple(int(v) for v in np.round(f - q))
    verts.append(fw)
    return len(verts) - 1, tuple(int(v) for v in np.round(f - fw))


edges = {}          # (a, b, wrap) -> id


def eid(a, sa, b, sb):
    w = tuple(np.subtract(sb, sa))
    if (a, b, w) in edges:
        return edges[(a, b, w)]
    if (b, a, tuple(-np.array(w))) in edges:
        return -edges[(b, a, tuple(-np.array(w)))]
    edges[(a, b, w)] = len(edges) + 1
    return edges[(a, b, w)]


faces, body = [], [[] for _ in range(N)]
for i, c in enumerate(cells):
    for f in c['faces']:
        j, n = f['j'], tuple(int(v) for v in f['n'])
        if j < i or (j == i and n < (0, 0, 0)):
            continue
        V = f['verts']
        d = X[j] + np.array(n) @ H - X[i]
        cen = V.mean(0)
        nrm = sum(np.cross(V[k] - cen, V[(k + 1) % len(V)] - cen) for k in range(len(V)))
        if nrm @ d < 0:
            V = V[::-1]
        ids = [vid(p) for p in V]
        loop = [eid(*ids[k], *ids[(k + 1) % len(ids)]) for k in range(len(ids))]
        faces.append(loop)
        fidx = len(faces)
        body[i].append(fidx)       # normal points out of cell i
        body[j].append(-fidx)      # and into cell j

sym = {0: '*', 1: '+', -1: '-'}
L = [f'// {name}: {N} cells, TCP foam, flat start (topology from Dutour Sikiric et al. data)',
     'TORUS_FILLED', 'PERIODS'] + [' '.join(f'{v:.15f}' for v in row) for row in Hbox]
L += ['', 'vertices']
# vertex 1 fixed: removes the 3 rigid-translation null modes of a periodic foam (needed for the Hessian)
L += [f'{k + 1} ' + ' '.join(f'{v:.15f}' for v in (q @ Hbox)) + (' fixed' if FIXV and k == 0 else '') for k, q in enumerate(verts)]
L += ['', 'edges']
for (a, b, w), k in sorted(edges.items(), key=lambda kv: kv[1]):
    if any(abs(x) > 1 for x in w):
        raise ValueError('edge wraps more than one period')
    L.append(f'{k} {a + 1} {b + 1} ' + ' '.join(sym[x] for x in w))
L += ['', 'faces'] + [f'{k + 1} ' + ' '.join(str(e) for e in lp) for k, lp in enumerate(faces)]
# all bodies volume 1 (TORUS_FILLED itself drops the one redundant constraint)
L += ['', 'bodies'] + [f'{i + 1} ' + ' '.join(str(f) for f in b) + ' volume 1' for i, b in enumerate(body)]
open(outp, 'w').write('\n'.join(L) + '\n')
print(f'{name}: {len(verts)} vertices, {len(edges)} edges, {len(faces)} faces, {N} bodies -> {outp}')
