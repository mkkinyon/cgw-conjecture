"""s13_cleanid.py -- the marked-pair identity for a generic row pair (x,y) inside the space with rows 1,2 FIXED.

For rows (x,y) = (3,4) and a mark {q,q'} with q' = rho_{x,y}^a(q) on an m-cycle (2 <= a <= m-2), the flip at an
A-pair (turn of the column cycle of (q,q') through y) is legal in S(R) (rows 1,2 fixed) iff that cycle avoids
rows 1,2 ("clean"); the switch (trade of rows x,y along a row cycle) is always legal.  So the exact identity
available with rows 1,2 fixed is the clean part only.  This script measures, on uniform squares, how the clean
restriction interacts with the A/B status: P[A], P[clean], P[A | clean], P[clean | A], P[clean | B], the length
distribution of the cycle through y, and the same conditioned on classes of rho_{1,2} (LONG = longest cycle > n/2).
usage: ./jm n nsamp thin seed | python3 s13_cleanid.py n
"""
import sys
from collections import Counter, defaultdict
from s8_real_bias import read_squares


def rho_of(L, x, y):
    n = len(L); pos = {L[y][q]: q for q in range(n)}
    return [pos[L[x][q]] for q in range(n)]


def cycles_of(perm):
    n = len(perm); seen = [False] * n; out = []
    for s in range(n):
        if seen[s]: continue
        c = []; q = s
        while not seen[q]: seen[q] = True; c.append(q); q = perm[q]
        out.append(c)
    return out


def main():
    n = int(sys.argv[1]); x, y = 2, 3
    st = defaultdict(Counter); lens = defaultdict(Counter); nsq = 0
    for L in read_squares(n, sys.stdin.buffer):
        nsq += 1
        colpos = [[0] * n for _ in range(n)]
        for r in range(n):
            for q in range(n): colpos[q][L[r][q]] = r
        r12 = cycles_of(rho_of(L, 0, 1)); long12 = max(len(c) for c in r12) > n / 2
        cls = 'LONG12' if long12 else 'short12'
        rho = rho_of(L, x, y)
        for cyc in cycles_of(rho):
            m = len(cyc)
            if m < 4: continue
            for i, q in enumerate(cyc):
                for a in range(2, m - 1):
                    qp = cyc[(i + a) % m]
                    # column cycle of (q,qp) through y: r -> row of column qp holding L(r,q)
                    Z = [y]; r = colpos[qp][L[y][q]]
                    while r != y: Z.append(r); r = colpos[qp][L[r][q]]
                    A = x not in Z
                    clean = (0 not in Z) and (1 not in Z)
                    for k in ('all', cls):
                        st[k]['n'] += 1; st[k]['A'] += A; st[k]['clean'] += clean
                        st[k]['A&clean'] += A and clean; st[k]['B&clean'] += (not A) and clean
                        lens[k][min(len(Z), n)] += 1
    print(f"n={n}: squares={nsq}, rows (x,y)=(3,4), marks on rho_xy-cycles of length>=4 at distance 2..m-2")
    for k in ('all', 'LONG12', 'short12'):
        c = st[k]; N = c['n']
        if not N: continue
        pA = c['A'] / N; pc = c['clean'] / N
        print(f"  [{k:8s}] N={N:7d}  P[A]={pA:.4f}  P[clean]={pc:.4f}  P[A|clean]={c['A&clean']/max(1,c['clean']):.4f}"
              f"  P[clean|A]={c['A&clean']/max(1,c['A']):.4f}  P[clean|B]={c['B&clean']/max(1,N-c['A']):.4f}"
              f"  P[B&clean]/P[A&clean]={c['B&clean']/max(1,c['A&clean']):.4f}")
    lc = lens['all']; N = sum(lc.values())
    print("  length of the column cycle through y: P[len<=k] for k = 2,3,5,10,n/4,n/2:",
          " ".join(f"{sum(v for l, v in lc.items() if l <= k)/N:.3f}" for k in (2, 3, 5, 10, n // 4, n // 2)))
    print("  (uniform on {2..n} would give:", " ".join(f"{(k-1)/(n-1):.3f}" for k in (2, 3, 5, 10, n // 4, n // 2)), ")")
    # P[clean | len = l] ~ (1 - l/n)^2 ?  report by length bins
    print("  P[clean | len in bin] vs ((n-len)/(n-1))*((n-len-1)/(n-2)) at bin midpoints: (not tracked; see lengths)")


if __name__ == "__main__":
    main()
