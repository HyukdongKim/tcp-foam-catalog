"""Candidate site sets (fractional coordinates) for known structures."""
import numpy as np


def cubic(a=1.0):
    return np.eye(3) * a


def hexagonal(a, c):
    return np.array([[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]])


def tetragonal(a, c):
    return np.diag([a, a, c])


def expand(gen, ops):
    """Apply symmetry ops (callables on frac coords) and dedupe mod 1."""
    pts = []
    for g in gen:
        for op in ops:
            p = np.mod(np.array(op(*g), float), 1.0)
            if not any(np.linalg.norm(((p - q + 0.5) % 1) - 0.5) < 1e-6 for q in pts):
                pts.append(p)
    return np.array(pts)


def kelvin():  # BCC
    return np.array([[0, 0, 0], [.5, .5, .5]]), cubic()


def wp():  # A15 / Pm-3n
    F = np.array([[0, 0, 0], [.5, .5, .5],
                  [.25, 0, .5], [.75, 0, .5],
                  [.5, .25, 0], [.5, .75, 0],
                  [0, .5, .25], [0, .5, .75]])
    return F, cubic()


def fcc_ops():
    tr = [(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)]
    return tr


def c15():  # MgCu2 Laves, Fd-3m: 8a (Mg) + 16d (Cu)  -> 24 sites
    # origin choice 2: 8a at 1/8,1/8,1/8 ; 16d at 1/2,1/2,1/2
    base8 = [(1/8, 1/8, 1/8), (7/8, 3/8, 3/8)]
    base16 = [(1/2, 1/2, 1/2), (1/4, 3/4, 0), (3/4, 0, 1/4), (0, 1/4, 3/4)]
    pts = []
    for b in base8 + base16:
        for t in fcc_ops():
            pts.append(np.mod(np.add(b, t), 1))
    return np.array(pts), cubic()


def c14(c_over_a=np.sqrt(8 / 3) * 1.0):  # MgZn2 P6_3/mmc, 12 sites
    x = -1 / 6  # ideal 6h x ~ -1/6 (i.e. 0.83)
    z4f = 1 / 16 * 1.0 + 0.0
    z4f = 0.0625
    F = [
        # 2a Zn
        (0, 0, 0), (0, 0, .5),
        # 4f Mg
        (1/3, 2/3, z4f), (2/3, 1/3, z4f + .5), (2/3, 1/3, -z4f), (1/3, 2/3, .5 - z4f),
        # 6h Zn
        (x, 2 * x, .25), (-2 * x, -x, .25), (x, -x, .25),
        (-x, -2 * x, .75), (2 * x, x, .75), (-x, x, .75),
    ]
    # c/a for C14 with unit sizes ~1.633 ; lattice: a, c
    a = 1.0
    return np.mod(np.array(F), 1), hexagonal(a, c_over_a * a)


def z_phase():  # Zr4Al3, P6/mmm, 7 sites
    F = [(0, 0, 0),
         (1/3, 2/3, .5), (2/3, 1/3, .5),
         (.5, 0, 0), (0, .5, 0), (.5, .5, 0),
         ]
    # Zr4Al3: Al 1a(0,0,0)? use common TCP description: 1a, 2d, 3f + ... (7 = 1+2+... )
    # Z phase cell: 7 atoms: 1a(0,0,0) Z14? , 2e(0,0,z~0.31)?? -> use Frank-Kasper
    F = [(0, 0, 0),
         (1/3, 2/3, .5), (2/3, 1/3, .5),
         (.5, 0, .5), (0, .5, .5), (.5, .5, .5),
         (0, 0, .5)]
    return np.mod(np.array(F, float), 1), hexagonal(1.0, 1.0)


def sigma():  # sigma phase P4_2/mnm, 30 sites (ideal coords)
    # Wyckoff 2a, 4f(x=.4), 8i(x=.13,y=.46),8i(.74,.07),8j(x=.18,z=.25)
    def ops_P42mnm(x, y, z):
        return [(x, y, z), (-x, -y, z), (-y + .5, x + .5, z + .5), (y + .5, -x + .5, z + .5),
                (-x + .5, y + .5, -z + .5), (x + .5, -y + .5, -z + .5), (y, x, -z), (-y, -x, -z),
                (-x, -y, -z), (x, y, -z), (y + .5, -x + .5, -z + .5), (-y + .5, x + .5, -z + .5),
                (x + .5, -y + .5, z + .5), (-x + .5, y + .5, z + .5), (-y, -x, z), (y, x, z)]
    gens = [(0, 0, 0), (.3981, .3981, 0), (.4632, .1316, 0), (.7376, .0653, 0), (.1823, .1823, .2524)]
    pts = []
    for g in gens:
        for p in ops_P42mnm(*g):
            p = np.mod(np.array(p, float), 1)
            if not any(np.linalg.norm(((p - q + .5) % 1) - .5) < 1e-4 for q in pts):
                pts.append(p)
    return np.array(pts), tetragonal(1.0, 0.52)


def frac_to_cart(F, H):
    return np.asarray(F) @ H


def normalize_volume(F, H):
    """Scale lattice so each cell has unit volume."""
    N = len(F)
    s = (N / abs(np.linalg.det(H))) ** (1 / 3)
    return F, H * s
