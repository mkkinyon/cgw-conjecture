"""s13_rowpairs.py -- decorrelation of the cycle types of two DISJOINT row pairs of a uniform Latin square.

Question (2) of the plan: is  P[rho_{x,y} has a macroscopic cycle | class of rho_{1,2}]  bounded below uniformly
in the class?  For each sampled square we take the pair (1,2) and the disjoint pairs (3,4),(5,6),... and record
for every row pair: kappa (number of cycles), Lmax (longest cycle), N2 (2-cycles), and the indicator
LONG = [Lmax > n/2]  (P = ln 2 for a uniform permutation).  We report
  * P[LONG(x,y)] unconditionally and conditioned on classes of (1,2): kappa(1,2) = 1,2,3,4+; Lmax(1,2) in bins;
    LONG(1,2) / not; N2(1,2) = 0 / >= 1;
  * the same with (1,2) and (x,y) exchanged (symmetric by exchangeability; a check on the sampler);
  * Pearson correlations of kappa, Lmax, log Lmax across the two pairs, with standard errors from per-square blocks;
  * P[LONG(x,y) | all cycles of (1,2) have length <= n/5]  if that class is populated.
usage: ./jm n nsamp thin seed | python3 s13_rowpairs.py n [--pairs P]
"""
import sys, math
from collections import defaultdict
from s8_real_bias import read_squares


def cyc_stats(L, x, y):
    n = len(L)
    pos = {L[y][q]: q for q in range(n)}
    rho = [pos[L[x][q]] for q in range(n)]
    seen = [False] * n; lens = []
    for s in range(n):
        if seen[s]: continue
        c = 0; q = s
        while not seen[q]: seen[q] = True; q = rho[q]; c += 1
        lens.append(c)
    return len(lens), max(lens), sum(1 for l in lens if l == 2)


def main():
    n = int(sys.argv[1]); args = sys.argv[2:]
    P = int(args[args.index('--pairs') + 1]) if '--pairs' in args else (n - 2) // 2
    recs = []   # (square index, stats12, statsxy)
    nsq = 0
    for L in read_squares(n, sys.stdin.buffer):
        s12 = cyc_stats(L, 0, 1)
        for i in range(P):
            x, y = 2 + 2 * i, 3 + 2 * i
            if y >= n: break
            recs.append((nsq, s12, cyc_stats(L, x, y)))
        nsq += 1
    N = len(recs)
    long_ = lambda s: s[1] > n / 2
    def cond(name, f12, fxy=long_):
        sel = [r for r in recs if f12(r[1])]
        if not sel: print(f"  {name:38s}  (empty)"); return
        p = sum(fxy(r[2]) for r in sel) / len(sel)
        # block standard error over squares
        bysq = defaultdict(list)
        for r in sel: bysq[r[0]].append(fxy(r[2]))
        means = [sum(v) / len(v) for v in bysq.values()]
        se = (sum((m - p) ** 2 for m in means) / max(1, len(means) - 1)) ** 0.5 / len(means) ** 0.5
        print(f"  {name:38s}  P[LONG(x,y)] = {p:.4f} +- {se:.4f}   (n_pairs={len(sel)}, squares={len(bysq)})")
    print(f"n={n}: squares={nsq} row-pair records={N}  (LONG = longest cycle > n/2; ln2 = {math.log(2):.4f})")
    cond("unconditional", lambda s: True)
    for k in (1, 2, 3):
        cond(f"kappa(1,2) = {k}", lambda s, k=k: s[0] == k)
    cond("kappa(1,2) >= 4", lambda s: s[0] >= 4)
    cond("kappa(1,2) >= 6", lambda s: s[0] >= 6)
    cond("LONG(1,2)", long_)
    cond("not LONG(1,2)  (Lmax <= n/2)", lambda s: not long_(s))
    cond("Lmax(1,2) <= n/3", lambda s: s[1] <= n / 3)
    cond("Lmax(1,2) <= n/4", lambda s: s[1] <= n / 4)
    cond("Lmax(1,2) <= n/5", lambda s: s[1] <= n / 5)
    cond("Lmax(1,2) >= 0.9 n", lambda s: s[1] >= 0.9 * n)
    cond("N2(1,2) = 0", lambda s: s[2] == 0)
    cond("N2(1,2) >= 1", lambda s: s[2] >= 1)
    cond("N2(1,2) >= 2", lambda s: s[2] >= 2)
    print("  --- exchanged roles (should agree within error) ---")
    swap = [(r[0], r[2], r[1]) for r in recs]
    recs_saved = recs[:]
    recs[:] = swap
    cond("not LONG(x,y) -> P[LONG(1,2)]", lambda s: not long_(s))
    cond("kappa(x,y) >= 4 -> P[LONG(1,2)]", lambda s: s[0] >= 4)
    recs[:] = recs_saved
    # correlations
    def corr(f):
        a = [f(r[1]) for r in recs]; b = [f(r[2]) for r in recs]
        ma = sum(a) / N; mb = sum(b) / N
        cov = sum((u - ma) * (v - mb) for u, v in zip(a, b)) / N
        va = sum((u - ma) ** 2 for u in a) / N; vb = sum((v - mb) ** 2 for v in b) / N
        return cov / (va * vb) ** 0.5
    print(f"  corr(kappa(1,2), kappa(x,y)) = {corr(lambda s: s[0]):+.4f}")
    print(f"  corr(Lmax(1,2), Lmax(x,y))   = {corr(lambda s: s[1]):+.4f}")
    print(f"  corr(LONG(1,2), LONG(x,y))   = {corr(lambda s: float(long_(s))):+.4f}")
    print(f"  corr(N2(1,2), N2(x,y))       = {corr(lambda s: s[2]):+.4f}")
    print(f"  (naive s.e. of a correlation on {N} records: {1/N**0.5:.4f}; records within a square are correlated)")
    # kappa joint table
    print("  P[kappa(x,y)=k' | kappa(1,2)=k]  (rows k=1..5, cols k'=1..6):")
    tab = defaultdict(lambda: defaultdict(int))
    for r in recs: tab[min(r[1][0], 6)][min(r[2][0], 6)] += 1
    for k in range(1, 6):
        tot = sum(tab[k].values())
        if tot: print(f"    k={k} (n={tot:6d}): " + " ".join(f"{tab[k][kp]/tot:.3f}" for kp in range(1, 7)))
    tot = N
    print(f"    all (n={tot:6d}): " + " ".join(f"{sum(tab[k][kp] for k in tab)/tot:.3f}" for kp in range(1, 7)))


if __name__ == "__main__":
    main()
