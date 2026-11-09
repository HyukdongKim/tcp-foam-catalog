"""Exhaustive (signature-deduplicated) BFS over T1 move sequences, WITHOUT relaxing intermediates;
only triangle-free end states (no 3-sided faces, which collapse back under relaxation) are relaxed.
usage: python t1_stable.py START REPS DEPTH WORKERS [r=1] [maxstates=200000]
Writes rerun/t1s_<tag>.jsonl"""
import sys, os, json, time, warnings
import numpy as np
from multiprocessing import Pool
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from trifoam import from_tcp, TriFoam
from t1search import evaluate, ser


def main():
    start, reps, depth, nw = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    r = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    maxst = int(sys.argv[6]) if len(sys.argv) > 6 else 200000
    reps = tuple(int(x) for x in reps.split('x'))
    os.chdir(os.path.join(HERE, '..')); os.makedirs('rerun', exist_ok=True)
    tag = f"{start.split(':')[-1]}_{'x'.join(map(str, reps))}_r{r}"
    R = {d['name']: d for d in json.load(open('data/tcp_all.json'))}
    f0 = from_tcp(R[start]).supercell(reps)
    seen = {f0.signature()}
    level = [(f0.tets, [])]
    stable = []
    t0 = time.time()
    for d in range(1, depth + 1):
        nxt = []
        for tets, path in level:
            f = TriFoam(f0.X, f0.H, []); f.tets = tets; f._index()
            for mv in f.moves():
                g = f.apply(mv)
                s = g.signature()
                if s in seen:
                    continue
                seen.add(s)
                if g.check() is not None:
                    continue
                p = path + [mv[0]]
                tri_free = 3 not in g.face_hist()
                if d < depth:
                    nxt.append((g.tets, p))
                if tri_free:
                    stable.append((g, p, d))
            if len(seen) > maxst:
                break
        print(f'[{tag}] depth {d}: {len(nxt)} topologies, triangle-free so far {len(stable)}  ({time.time()-t0:.0f}s)', flush=True)
        level = nxt
        if len(seen) > maxst:
            print(f'[{tag}] state limit reached', flush=True); break
    out = open(f'rerun/t1s_{tag}.jsonl', 'a')
    pool = Pool(nw, maxtasksperchild=20)
    base = pool.apply(evaluate, ((f0.X, f0.H, sorted(f0.tets), r),))['cost']
    out.write(json.dumps(dict(tag=tag, depth=0, path=[], cost=base, ok=True, faces=f0.face_hist(), nf=f0.nfaces_hist())) + '\n')
    print(f'[{tag}] start cost {base:.6f}; relaxing {len(stable)} triangle-free topologies', flush=True)
    best = None
    for k, res in enumerate(pool.imap(evaluate, [(g.X, g.H, sorted(g.tets), r) for g, _, _ in stable]), 1):
        g, p, d = stable[k - 1]
        rec = dict(tag=tag, depth=d, path=p, cost=res['cost'], ok=res['ok'], flip=res.get('flip'),
                   faces=g.face_hist(), nf=g.nfaces_hist(), t=res['t'])
        if res['ok']:
            rec['tets'] = ser(g.tets); rec['X'] = res['Xn']; rec['H'] = res['Hn']
            if best is None or res['cost'] < best:
                best = res['cost']
                flag = '  <<< 시작점보다 낮음!' if best < base else ''
                print(f'[{tag}] new best {best:.6f} (Δ {best-base:+.5f}) faces={g.face_hist()} nf={g.nfaces_hist()}{flag}', flush=True)
        out.write(json.dumps(rec) + '\n'); out.flush()
        if k % 20 == 0 or k == len(stable):
            print(f'[{tag}] {k}/{len(stable)} relaxed', flush=True)
    print(f'[{tag}] done. start {base:.6f} best {best}', flush=True)


if __name__ == '__main__':
    main()
