"""Build TCP foam records from crystal structures (CIF) in data/source_cif/ -> data/tcp_cif.json.
Every crystallographic site is a foam cell centre. The topology is the Voronoi topology, obtained from the
periodic Delaunay triangulation (trifoam.from_delaunay); a structure is accepted only if every cell has
12/14/15/16 faces and every face is a pentagon or hexagon (TCP). R-3 (hexagonal setting) is reduced to the
primitive rhombohedral cell. Duplicates of data/tcp_all.json are reported by Weisfeiler-Lehman hash.
Requires gemmi.  usage: python code/cif_build.py"""
import collections, json, os, sys, warnings
import numpy as np
import gemmi
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from trifoam import from_delaunay, from_tcp
from t1_deep import wl_hash

SOURCES = [('cif:R-phase', '2310299.cif'), ('cif:P-phase', '2310312.cif'), ('cif:P-phase_b', '2311049.cif')]
HEX_TO_RHOMB = np.array([[2 / 3, 1 / 3, 1 / 3], [-1 / 3, 1 / 3, 1 / 3], [-1 / 3, -2 / 3, 1 / 3]])


def dedup(F, tol=1e-4):
    out = []
    for p in np.mod(F, 1):
        if not any(np.abs(((p - q + .5) % 1) - .5).max() < tol for q in out):
            out.append(p)
    return np.array(out)


def record(name, fn):
    st = gemmi.read_small_structure(fn)
    F = dedup([[s.fract.x, s.fract.y, s.fract.z] for s in st.get_all_unit_cell_sites()])
    H = np.array(st.cell.orth.mat.tolist()).T          # rows = lattice vectors
    if st.spacegroup_hm.replace(' ', '').startswith('R-3:H'):
        H = HEX_TO_RHOMB @ H
        F = dedup([np.linalg.solve(HEX_TO_RHOMB.T, p) for p in F])
    f = from_delaunay(F, H)
    assert f.check() is None
    deg = collections.Counter()
    for c in f.edge:
        for v in c:
            deg[v[0]] += 1
    Z = [deg[i] for i in range(len(F))]
    if not (set(Z) <= {12, 14, 15, 16} and set(f.face_hist()) <= {5, 6}):
        raise ValueError(f'{fn}: Voronoi topology is not TCP')
    S = [set() for _ in F]
    for u, v in f.edge:
        n = tuple(int(x) for x in np.subtract(v[1], u[1]))
        S[u[0]].add((v[0], n)); S[v[0]].add((u[0], tuple(-x for x in n)))
    return dict(name=name, group=st.spacegroup_hm, source='COD ' + os.path.basename(fn)[:-4], N=len(F), meanZ=float(np.mean(Z)),
                F=F.tolist(), H=H.tolist(), Z=Z, S=[sorted([int(j), list(n)] for j, n in s) for s in S])


if __name__ == '__main__':
    os.chdir(os.path.join(HERE, '..'))
    T = json.load(open('data/tcp_all.json'))
    known = collections.defaultdict(list)
    for t in T:
        known[(t['N'], wl_hash(from_tcp(t)))].append(t['name'])
    out, seen = [], {}
    for name, fn in SOURCES:
        r = record(name, 'data/source_cif/' + fn)
        h = wl_hash(from_tcp(r))
        comp = dict(sorted(collections.Counter(r['Z']).items()))
        dup = known.get((r['N'], h), []) + ([seen[h]] if h in seen else [])
        print(f"{name:16s} {r['source']}  {r['group']}  N={r['N']}  mean Z={r['meanZ']:.4f}  {comp}  duplicate of: {dup or 'none'}")
        if not dup:
            out.append(r); seen[h] = name
    json.dump(out, open('data/tcp_cif.json', 'w'))
