"""Build the Shoemaker pseudo-hexagonal phases F, J, K (data/source_cgd/*-phase.cgd, from the
Dutour Sikiric space-fullerene dataset, SpecialStructure/) into data/tcp_shoemaker.json
(same record format as data/tcp_all.json) and check the dual triangulation."""
import json, os, sys, collections
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cgd import read, expand_graph
from tcpbuild import neighbours, build_cells
from trifoam import from_tcp

os.chdir(os.path.join(HERE, '..'))
out = []
for ph in ('F', 'J', 'K'):
    rec = read(f'data/source_cgd/{ph}-phase.cgd')[0]
    F, H, Z, group, sets = expand_graph(rec)
    F = np.mod(F, 1)
    if sets is None:
        _, sets, ok = neighbours(F, H, Z)
        assert ok, 'neighbour graph not symmetric'
    r = dict(name=f'shoemaker:{ph}-phase', group=group, N=len(F), meanZ=float(np.mean(Z)), F=F.tolist(), H=H.tolist(),
             Z=[int(z) for z in Z], S=[sorted([int(j), [int(v) for v in n]] for j, n in s) for s in sets])
    build_cells(F, H, Z, [set((j, tuple(n)) for j, n in s) for s in r['S']])
    f = from_tcp(r)
    assert f.check() is None
    print(r['name'], group, 'N', r['N'], 'mean Z %.4f' % r['meanZ'], dict(sorted(collections.Counter(r['Z']).items())), 'faces', f.face_hist())
    out.append(r)
json.dump(out, open('data/tcp_shoemaker.json', 'w'))
