"""s13_trade.py -- the trade chain on real squares (session 13, sec:s13tradechain).

From uniform squares (jm.c) and marks P={j,jp}: run the symmetric chain
  with prob 1/2: pick rows x!=y in 3..n; if (x,y) flippable (j, jp in
                 different cycles of rho_{x,y}) trade the rho-cycle through jp;
  with prob 1/2: pick row x in 3..n; if C_x avoids rows 1,2, turn C_x.
Rows 1,2 (indices 0,1) are never touched.  Record the bit eps = [u~v in pi_P]
after every step, and the number of trades performed so far.  Report
  * flippable fraction, and flippable fraction split by the frame relation
    of (x,y): adjacent (pi(x)=y or pi(y)=x), same cycle non-adjacent, different cycles;
  * autocorrelation E[eps_s eps_{s+t}] as a function of the number of trades m
    between the two times, compared with <eps,P^m eps>_N (exact spectral formula,
    N=n) -- the pure walk with iid transpositions.
usage: ./jm n nsamp thin seed | python3 s13_trade.py n [--inst I] [--steps T] [--seed s]
"""
import sys, random, math
import numpy as np
from verify_identity import sigma_of
from pf_cl import cycles_of
from pf_master import pi_P


def spectral(N):
    b = np.arange(1, N, dtype=np.float64)[:, None]; j = np.arange(0, N, dtype=np.float64)[None, :]
    a = N - b - j; ok = a >= b
    c = np.where(ok, (a - b + 1) / ((a + j + 1) * (b + j)), 0.0)
    r = np.where(ok, (a * (a - 1) + b * (b - 3) - j * (j + 3)) / (N * (N - 1)), 0.0)
    nu = np.where(ok, 2 * c * c * (1 - r), 0.0)
    return r[ok], nu[ok]


class Square:
    def __init__(self, L, j, jp):
        self.n = len(L); self.L = [list(r) for r in L]; self.j = j; self.jp = jp
        self.pos = [[0] * self.n for _ in range(self.n)]  # pos[r][s] = column of symbol s in row r
        for r in range(self.n):
            for q in range(self.n): self.pos[r][self.L[r][q]] = q
        self.rebuild_frame()

    def rebuild_frame(self):
        self.pi = pi_P(self.L, self.j, self.jp)
        cyc = cycles_of(self.pi)
        self.cid = [0] * self.n
        for i, c in enumerate(cyc):
            for x in c: self.cid[x] = i
        self.cycles = cyc

    def rho_cycle_through(self, x, y, q0):
        """cycle of rho_{x,y} (columns) through column q0: q -> column where row y holds L[x][q]"""
        cyc = [q0]; q = self.pos[y][self.L[x][q0]]
        while q != q0:
            cyc.append(q); q = self.pos[y][self.L[x][q]]
        return cyc

    def flippable(self, x, y):
        return self.jp not in set(self.rho_cycle_through(x, y, self.j))

    def trade(self, x, y):
        Q = self.rho_cycle_through(x, y, self.jp)
        for q in Q:
            a, b = self.L[x][q], self.L[y][q]
            self.L[x][q], self.L[y][q] = b, a
            self.pos[x][b] = q; self.pos[y][a] = q
        self.rebuild_frame()

    def turn(self, c):
        for r in c:
            a, b = self.L[r][self.j], self.L[r][self.jp]
            self.L[r][self.j], self.L[r][self.jp] = b, a
            self.pos[r][b] = self.j; self.pos[r][a] = self.jp
        self.rebuild_frame()

    def eps(self):
        return 1 if self.cid[0] == self.cid[1] else -1

    def relation(self, x, y):
        if self.pi[x] == y or self.pi[y] == x: return 'adj'
        return 'same' if self.cid[x] == self.cid[y] else 'diff'


def is_latin(L):
    n = len(L)
    return all(sorted(r) == list(range(n)) for r in L) and all(sorted(L[r][q] for r in range(n)) == list(range(n)) for q in range(n))


def main():
    n = int(sys.argv[1]); args = sys.argv[2:]
    def opt(name, default):
        return int(args[args.index(name) + 1]) if name in args else default
    ninst = opt("--inst", 2); steps = opt("--steps", 4000); seed = opt("--seed", 5)
    rng = random.Random(seed)
    from s8_real_bias import read_squares
    flip = {'adj': [0, 0], 'same': [0, 0], 'diff': [0, 0]}
    maxlag = 64
    acc = np.zeros(maxlag + 1); cnt = np.zeros(maxlag + 1)
    ninst_done = 0; nsq = 0; latin_checks = 0
    for L in read_squares(n, sys.stdin.buffer):
        nsq += 1
        sig = sigma_of(L); pairs = []
        for cyc in cycles_of(sig):
            m = len(cyc)
            if m < 4: continue
            for ji, j in enumerate(cyc):
                for alpha in range(2, m - 1): pairs.append((j, cyc[(ji + alpha) % m]))
        if not pairs: continue
        for (j, jp) in rng.sample(pairs, min(ninst, len(pairs))):
            S = Square(L, j, jp)
            sig0 = sigma_of(S.L)
            eps_hist = []; trades_hist = []; ntr = 0
            for t in range(steps):
                if rng.random() < 0.5:
                    x, y = rng.sample(range(2, n), 2)
                    rel = S.relation(x, y); f = S.flippable(x, y)
                    flip[rel][0] += 1; flip[rel][1] += f
                    if f: S.trade(x, y); ntr += 1
                else:
                    x = rng.randrange(2, n); c = S.cycles[S.cid[x]]
                    if S.cid[0] != S.cid[x] and S.cid[1] != S.cid[x]: S.turn(c)
                eps_hist.append(S.eps()); trades_hist.append(ntr)
            # invariants
            assert sigma_of(S.L) == sig0, "sigma changed!"
            if latin_checks < 5: assert is_latin(S.L); latin_checks += 1
            e = np.array(eps_hist); tr = np.array(trades_hist)
            # autocorrelation by number of trades between times (subsample start times)
            burn = steps // 4
            for s0 in range(burn, steps - 1, 7):
                for s1 in range(s0, min(steps, s0 + 400)):
                    m = tr[s1] - tr[s0]
                    if m > maxlag: break
                    acc[m] += e[s0] * e[s1]; cnt[m] += 1
            ninst_done += 1
    print(f"n={n}: squares={nsq} instances={ninst_done} steps/instance={steps}; sigma preserved, Latin checked")
    for rel in ('adj', 'same', 'diff'):
        a, f = flip[rel]
        print(f"  flippable fraction | {rel:4s}: {f/max(1,a):.4f}  ({a} pairs)")
    r, nu = spectral(n)
    # baseline: uniform derangement frame, iid uniform transpositions avoiding rows 0,1
    bacc = np.zeros(maxlag + 1); bcnt = np.zeros(maxlag + 1)
    for it in range(ninst_done):
        while True:
            p = list(range(n)); rng.shuffle(p)
            if all(p[i] != i for i in range(n)) and p[0] != 1 and p[1] != 0: break   # real frames are derangements with pi(1)!=2, pi(2)!=1
        def same(p):
            x = 0
            while True:
                x = p[x]
                if x == 1: return 1
                if x == 0: return -1
        hist = []
        for t in range(3 * maxlag + steps // 4):
            x, y = rng.sample(range(2, n), 2)
            p = [p[q] for q in range(n)]
            p[x], p[y] = p[y], p[x]   # (x y) o p : swap images? we want pi -> (x y) o pi: relabel values x<->y
            p = [y if v == x else x if v == y else v for v in p]
            p[x], p[y] = p[y], p[x]   # undo the first swap (we only want left multiplication)
            hist.append(same(p))
        h = np.array(hist)
        for s0 in range(0, len(h) - maxlag, 3):
            for m in range(maxlag + 1):
                bacc[m] += h[s0] * h[s0 + m]; bcnt[m] += 1
    print("  trades m : E[eps eps] chain     baseline (uniform derangement frame, pi(u)!=v, pi(v)!=u; iid transp. avoiding u,v)   exact walk (all transp.)")
    for m in [0, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]:
        if cnt[m] == 0: continue
        print(f"  {m:3d} : {acc[m]/cnt[m]:+.4f}        {bacc[m]/max(1,bcnt[m]):+.4f}        {(nu*r**m).sum():+.4f}")


if __name__ == "__main__":
    main()
