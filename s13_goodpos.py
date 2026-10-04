"""s13_goodpos.py -- goodness of a pair of free rows on D versus their relative position on D (first candidate).

For joined instances (x1 ~ y1 at (q1,c1)): for all unordered pairs {z,z'} of free rows on D, record whether the pair is
good (row cycle through c1 avoids q1), whether it is interleaving (exactly one of x1,y1 strictly between), and the cyclic
distance min(d, ell-d) along D.  Checks the adjacency fact: (z, pi(z)) is never good.
usage: ./jmfix ... | python3 s13_goodpos.py n [--marks K] [--seed s]
"""
import sys, random
from collections import defaultdict
from s8_real_bias import read_squares
from pf_cl import cycles_of
from s13_amtwo import rowcycle_cols, colcycle_rows
from s13_hazard import pos_table, first_pair, candidates

SPECIAL = (0, 1)


def main():
    n = int(sys.argv[1]); args = sys.argv[2:]
    def opt(name, default, conv=int):
        return conv(args[args.index(name) + 1]) if name in args else default
    K = opt('--marks', 3); seed = opt('--seed', 1)
    rng = random.Random(seed)
    sig = None; marks = []
    byint = defaultdict(lambda: [0, 0]); bydist = defaultdict(lambda: [0, 0]); adj = [0, 0]; nsq = 0; joined = 0
    for L in read_squares(n, sys.stdin.buffer):
        nsq += 1
        if sig is None:
            sig = L[1][:]
            for cyc in cycles_of(sig):
                m = len(cyc)
                if m < 4: continue
                for ji, j in enumerate(cyc):
                    for alpha in range(2, m - 1):
                        marks.append((j, cyc[(ji + alpha) % m]))
        pos = pos_table(L, n)
        for (p, pp) in rng.sample(marks, min(K, len(marks))):
            x1, y1 = first_pair(L, n, pos, p, pp)
            nu, cands = candidates(L, n, pos, p, pp, x1, y1, 1)
            if not cands: continue
            q, c = cands[0]
            D = colcycle_rows(L, n, q, c, x1)
            if y1 not in D: continue
            joined += 1
            ell = len(D); m0 = D.index(y1); excl = set(SPECIAL) | {x1, y1}
            for i in range(ell):
                for j in range(i + 1, ell):
                    z, zp = D[i], D[j]
                    if z in excl or zp in excl: continue
                    good = q not in rowcycle_cols(L, n, z, zp, c)
                    inter = (i < m0 < j)            # x1 at position 0, y1 at m0: interleaving iff exactly one of i,j in (0,m0)
                    d = j - i; d = min(d, ell - d)
                    byint[inter][0] += good; byint[inter][1] += 1
                    bydist[min(d, 12)][0] += good; bydist[min(d, 12)][1] += 1
                    if d == 1: adj[0] += good; adj[1] += 1
    print(f'n={n} squares={nsq} joined instances={joined}')
    print(f'  good fraction: interleaving {byint[True][0]/byint[True][1]:.3f} ({byint[True][1]} pairs), '
          f'same-arc {byint[False][0]/byint[False][1]:.3f} ({byint[False][1]} pairs)')
    print(f'  adjacent pairs (z, pi(z)): good {adj[0]}/{adj[1]}')
    print('  good fraction by cyclic distance along D (12 = >=12):',
          {d: round(v[0] / v[1], 3) for d, v in sorted(bydist.items()) if v[1]})


if __name__ == '__main__':
    main()
