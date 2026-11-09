"""Deep exhaustive T1 enumeration (parallel BFS), triangle-free end states relaxed.
Deduplication by a Weisfeiler-Lehman hash of the quotient graph (cells, edges labelled by polygon size),
which is finer than TriFoam.signature(); symmetric/translated copies always merge.
Triangle-free states (no 3-gon faces) occur only at odd depth from a TCP start.
usage: python t1_deep.py START REPS DEPTH WORKERS [r=1] [maxlevel=3000000]
writes rerun/t1d_<tag>.jsonl (relaxed triangle-free states) and rerun/t1d_<tag>_stable.pkl."""
import sys, os, json, time, warnings, pickle, zlib
from collections import defaultdict, Counter
import numpy as np
from multiprocessing import Pool
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from trifoam import from_tcp, TriFoam


def wl_hash(f, rounds=None):
    N = f.N
    adj = defaultdict(list)
    nf = Counter()
    for (u, v), pairs in f.edge.items():
        d = len(pairs)
        adj[u[0]].append((d, v[0]))
        adj[v[0]].append((d, u[0]))
        nf[u[0]] += 1; nf[v[0]] += 1
    col = [hash((nf[i], tuple(sorted(d for d, _ in adj[i])))) for i in range(N)]
    for _ in range(rounds or min(N, 8)):
        new = [hash((col[i], tuple(sorted((d, col[j]) for d, j in adj[i])))) for i in range(N)]
        if len(set(new)) == len(set(col)):
            col = new; break
        col = new
    return hash((len(f.tets), tuple(sorted(col))))


G = {}


def pack(tets):
    return zlib.compress(pickle.dumps(sorted(tets), protocol=4), 6)


def unpack(b):
    return set(pickle.loads(zlib.decompress(b)))


def init(X, H):
    G['X'] = X; G['H'] = H


def expand(args):
    """Expand a chunk of parent states. last=True: return only valid triangle-free children."""
    chunk, last = args
    out = {}
    for blob, path in chunk:
        tets = unpack(blob)
        f = TriFoam(G['X'], G['H'], []); f.tets = tets; f._index()
        for mv in f.moves():
            g = f.apply(mv)
            ntri = sum(1 for p in g.edge.values() if len(p) == 3)
            if last and ntri:
                continue
            h = wl_hash(g)
            if h in out:
                continue
            if g.check() is not None:
                continue
            out[h] = (pack(g.tets), path + [mv[0]], ntri)
    return out


def main():
    start, reps, depth, nw = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    r = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    maxlevel = int(sys.argv[6]) if len(sys.argv) > 6 else 3000000
    reps = tuple(int(x) for x in reps.split('x'))
    os.chdir(os.path.join(HERE, '..')); os.makedirs('rerun', exist_ok=True)
    tag = f"{start.split(':')[-1]}_{'x'.join(map(str, reps))}_r{r}"
    R = {d['name']: d for d in json.load(open('data/tcp_all.json'))}
    f0 = from_tcp(R[start]).supercell(reps)
    seen = {wl_hash(f0)}
    level = [(pack(f0.tets), [])]
    stable = {}
    t0 = time.time()
    pool = Pool(nw, initializer=init, initargs=(f0.X, f0.H))
    for d in range(1, depth + 1):
        last = d == depth
        csz = max(1, min(200, len(level) // (nw * 8) + 1))
        chunks = [(level[k:k + csz], last) for k in range(0, len(level), csz)]
        nxt = []
        done = 0
        for out in pool.imap_unordered(expand, chunks):
            done += 1
            for h, (tets, path, ntri) in out.items():
                if h in seen:
                    continue
                seen.add(h)
                if ntri == 0:
                    stable[h] = (unpack(tets), path, d)
                if not last:
                    nxt.append((tets, path))
            if done % max(1, len(chunks) // 20) == 0:
                print(f'[{tag}] depth {d}: chunk {done}/{len(chunks)}  new {len(nxt)}  triangle-free {len(stable)}  ({time.time()-t0:.0f}s)', flush=True)
        print(f'[{tag}] depth {d} done: {len(nxt)} states kept, triangle-free total {len(stable)}  ({time.time()-t0:.0f}s)', flush=True)
        level = nxt
        if len(level) > maxlevel:
            print(f'[{tag}] level too large ({len(level)} > {maxlevel}); stopping enumeration', flush=True)
            break
    pickle.dump((f0.X, f0.H, stable), open(f'rerun/t1d_{tag}_stable.pkl', 'wb'))
    from t1search import evaluate, ser
    out = open(f'rerun/t1d_{tag}.jsonl', 'a')
    base = pool.apply(evaluate, ((f0.X, f0.H, sorted(f0.tets), r),))['cost']
    out.write(json.dumps(dict(tag=tag, depth=0, path=[], cost=base, ok=True)) + '\n')
    items = list(stable.values())
    print(f'[{tag}] start cost {base:.6f}; relaxing {len(items)} triangle-free topologies', flush=True)
    best = None
    for k, res in enumerate(pool.imap_unordered(_eval, [(f0.X, f0.H, tets, r, path, d) for tets, path, d in items], chunksize=1), 1):
        rec = dict(tag=tag, depth=res['d'], path=res['path'], cost=res['cost'], ok=res['ok'], flip=res.get('flip'),
                   faces=res['faces'], nf=res['nf'], t=res['t'])
        if res['ok'] and res['cost'] < base + 0.02:
            rec['tets'] = res['tets']; rec['X'] = res['Xn']; rec['H'] = res['Hn']
        out.write(json.dumps(rec) + '\n'); out.flush()
        if res['ok'] and (best is None or res['cost'] < best):
            best = res['cost']
            flag = '  <<< 시작점보다 낮음!' if best < base else ''
            print(f'[{tag}] new best {best:.6f} (Δ {best-base:+.5f}) faces={res["faces"]} nf={res["nf"]}{flag}', flush=True)
        if k % 50 == 0 or k == len(items):
            print(f'[{tag}] {k}/{len(items)} relaxed  ({time.time()-t0:.0f}s)', flush=True)
    print(f'[{tag}] done. start {base:.6f} best {best}', flush=True)


def _eval(a):
    X, H, tets, r, path, d = a
    from t1search import evaluate, ser
    res = evaluate((X, H, sorted(tets), r))
    f = TriFoam(X, H, []); f.tets = tets; f._index()
    res.update(path=path, d=d, faces=f.face_hist(), nf=f.nfaces_hist(), tets=ser(tets))
    return res


if __name__ == '__main__':
    main()
