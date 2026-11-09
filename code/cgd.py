"""Read Systre .cgd crystal nets (Dutour Sikiric space-fullerene data):
expand atoms (= foam cell centres) by space-group symmetry."""
import numpy as np
import gemmi


def lattice(a, b, c, al, be, ga):
    al, be, ga = np.radians([al, be, ga])
    ax = [a, 0, 0]
    bx = [b * np.cos(ga), b * np.sin(ga), 0]
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    cz = np.sqrt(max(c * c - cx * cx - cy * cy, 0))
    return np.array([ax, bx, [cx, cy, cz]])


def _sg(name):
    for cand in ((name + ':R', name + ':2', name + ':H', name) if name.startswith('R') else (name + ':2', name + ':H', name)):  # Systre: origin choice 2, hexagonal axes
        try:
            sg = gemmi.find_spacegroup_by_name(cand)
            if sg is not None:
                return sg
        except Exception:
            pass
    raise ValueError('unknown group ' + name)


def _num(x):
    """Accept fractional coordinates written as '1/3'."""
    if '/' in x:
        a, b = x.split('/'); return float(a) / float(b)
    return float(x)


def read(path):
    out, cur = [], None
    for line in open(path):
        t = line.split()
        if not t or t[0].startswith('#'):
            continue
        t[0] = t[0].lower()
        if t[0] == 'crystal':
            cur = dict(atoms=[], edges=[], name='?')
            out.append(cur)
        elif t[0] == 'name':
            cur['name'] = t[1]
        elif t[0] == 'cell':
            cur['cell'] = [float(x) for x in t[1:7]]
        elif t[0] == 'group':
            cur['group'] = t[1]
        elif t[0] == 'end':
            pass
        elif t[0] == 'atom':
            cur['atoms'].append((int(t[2]), [_num(x) for x in t[3:6]]))
            cur['edges'].append([])
        elif t[0] == 'edge':
            cur['edges'][-1].append([_num(x) for x in t[2:5]])
    return out


def expand(rec, tol=1e-3):
    sg = _sg(rec['group'])
    ops = sg.operations()
    F, Z = [], []
    for z, p in rec['atoms']:
        for op in ops:
            q = np.mod(np.array(op.apply_to_xyz(p)), 1.0)
            if not any(np.linalg.norm(((q - r + .5) % 1) - .5) < tol for r in F):
                F.append(q); Z.append(z)
    H = lattice(*rec['cell'])
    return np.array(F), H, np.array(Z), sg.xhm()


def expand_graph(rec, tol=1e-3):
    """Atoms plus explicit neighbour sets [(j, shift)] from the cgd edge list.
    Without edges (alloy files) returns sets=None (nearest-neighbour graph is used)."""
    F, H, Z, name = expand(rec, tol)
    if not any(rec['edges']):
        return F, H, Z, name, None
    ops = _sg(rec['group']).operations()

    def locate(q):
        q = np.array(q)
        for i, r in enumerate(F):
            d = q - r
            if np.abs(d - np.round(d)).max() < tol:
                return i, tuple(int(v) for v in np.round(d))
        raise ValueError('edge endpoint not an atom')

    sets = [set() for _ in F]
    for (z, p), E in zip(rec['atoms'], rec['edges']):
        for op in ops:
            i, si = locate(op.apply_to_xyz(p))
            for q in E:
                j, sj = locate(op.apply_to_xyz(q))
                n = tuple(np.subtract(sj, si))
                sets[i].add((j, n)); sets[j].add((i, tuple(-np.array(n))))
    return F, H, Z, name, sets
