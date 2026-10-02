"""s13_interleave.py -- the blocked (interleaving) fraction for CGW's moves on a generic row pair in X.

For a ladder pair (x,y) of a uniform square with a true mark (p,p'): a CGW flip at a column pair
(q,c) (q,c in different rho_{x,y}-cycles) for the rows x,y is legal in X iff one of the two
(q,c)-column cycles Z_x, Z_y (through x, through y; after the free switch if x,y share a cycle)
meets {1,2} in 0 or 2 rows.  It is BLOCKED iff rows 1,2 and rows x,y are INTERLEAVED in the
column-cycle structure of (q,c): (A-state) 1 in Z_x and 2 in Z_y or vice versa; (B-state)
x,1,y,2 alternate on the common cycle.  Interleaving is invariant under the switches of both
pairs.  This script measures, for joining pairs (q in Q_p, c not in Q_p) and for all pairs:
  P[interleaved], split by A/B state and by whether 1,2 lie in different sigma-cycles,
and the distribution of |Q_p| (length of the rho_{x,y}-cycle through the mark column).
usage: ./jm n nsamp thin seed | python3 s13_interleave.py n [--inst I] [--seed s]
"""
import sys, random
from collections import defaultdict
from verify_identity import sigma_of
from pf_cl import cycles_of
from pf_master import pi_P
from s8_real_bias import read_squares


def main():
    n = int(sys.argv[1]); args = sys.argv[2:]
    def opt(name, default):
        return int(args[args.index(name) + 1]) if name in args else default
    ninst = opt("--inst", 3); seed = opt("--seed", 23)
    rng = random.Random(seed)
    nsq = 0; N = 0
    stat = defaultdict(lambda: [0, 0])     # key -> [pairs, interleaved]
    qlen = defaultdict(int)
    for L in read_squares(n, sys.stdin.buffer):
        nsq += 1
        sig = sigma_of(L); pairs = []
        col = {L[0][c]: c for c in range(n)}
        for cyc in cycles_of(sig):
            m = len(cyc)
            if m < 4: continue
            for ji, jsym in enumerate(cyc):
                for alpha in range(2, m - 1): pairs.append((col[jsym], col[cyc[(ji + alpha) % m]]))
        if not pairs: continue
        pos = [[0] * n for _ in range(n)]
        for r in range(n):
            for q in range(n): pos[r][L[r][q]] = q
        colpos = [[0] * n for _ in range(n)]
        for r in range(n):
            for q in range(n): colpos[q][L[r][q]] = r
        # sigma-cycle id of each column (sigma as column permutation)
        sigcol = [0] * n
        for ci, cyc in enumerate(cycles_of(sig)):
            for jsym in cyc: sigcol[col[jsym]] = ci
        for (p, pp) in rng.sample(pairs, min(ninst, len(pairs))):
            pi = pi_P(L, p, pp); inv = [0] * n
            for r in range(n): inv[pi[r]] = r
            c1 = 0; r = 0
            while True:
                r = pi[r]; c1 += 1
                if r == 0: break
            together = False; r = 0; d12 = 0
            while True:
                r = pi[r]; d12 += 1
                if r == 1: together = True; break
                if r == 0: break
            if together: K0 = min(d12, c1 - d12)
            else:
                c2 = 0; r = 1
                while True:
                    r = pi[r]; c2 += 1
                    if r == 1: break
                K0 = min(c1, c2)
            if K0 < 2: continue
            k = rng.randrange(1, K0); x, y = 0, 1
            for _ in range(k): x, y = inv[x], inv[y]
            rho = [pos[y][L[x][q]] for q in range(n)]
            # rho-cycle ids
            rcyc = [-1] * n; nc = 0
            for q in range(n):
                if rcyc[q] >= 0: continue
                r_ = q
                while rcyc[r_] < 0: rcyc[r_] = nc; r_ = rho[r_]
                nc += 1
            Qp = sum(1 for q in range(n) if rcyc[q] == rcyc[p]); qlen[min(Qp, 11)] += 1
            N += 1
            # all column pairs with q,c in different rho-cycles
            for q in range(n):
                for c in range(q + 1, n):
                    if rcyc[q] == rcyc[c]: continue
                    # column cycle structure of (q,c): position of rows x,y,0,1
                    def cyc_through(r0):
                        Z = [r0]; r_ = colpos[c][L[r0][q]]
                        while r_ != r0: Z.append(r_); r_ = colpos[c][L[r_][q]]
                        return Z
                    Zx = cyc_through(x); sx = set(Zx)
                    if y in sx:
                        # B-state: alternation of x,0,y,1 on the cycle
                        if 0 in sx and 1 in sx:
                            ix = 0; iy = Zx.index(y); i0 = Zx.index(0); i1 = Zx.index(1)
                            inter = ((i0 < iy) != (i1 < iy))     # exactly one of 0,1 on the arc x->y
                        else: inter = False
                        state = 'B'
                    else:
                        Zy = set(cyc_through(y))
                        inter = ((0 in sx and 1 in Zy) or (1 in sx and 0 in Zy))
                        state = 'A'
                    join = (rcyc[q] == rcyc[p]) != (rcyc[c] == rcyc[p])   # joining pair for Q_p
                    diffsig = (sigcol[q] != sigcol[c])
                    for key in (('all', state), ('all', '*'), ('join' if join else 'nonjoin', '*'),
                                ('sig-diff' if diffsig else 'sig-same', '*')):
                        s = stat[key]; s[0] += 1; s[1] += inter
    print(f"n={n}: squares={nsq} ladder-pair instances={N}")
    print("  |Q_p| distribution (2..10, 11+):", [qlen[l] for l in range(2, 12)], " expected ~N/n each =", N / n)
    print("  P[interleaved] (rows 1,2 vs rows x,y in the (q,c)-column cycles; q,c in different rho_{x,y}-cycles):")
    for key in sorted(stat):
        a, b = stat[key]
        print(f"    {key[0]:8s} {key[1]:2s}: pairs={a:9d}  P[interleaved]={b/max(1,a):.4f}")


if __name__ == "__main__":
    main()
