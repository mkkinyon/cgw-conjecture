"""s13_wtail.py -- the witness lemma (hyp:witness) on completions of a fixed rectangle: lower tails.

Reads squares (n*n bytes) from stdin (./jmfix output, rows 0,1 = (id, sigma)).  For sampled marks (p,p') of
class (m, alpha):
 (i)  LADDER: pairs (x_k,y_k) = (pi^-k(0), pi^-k(1)), 1 <= k < K0.  For each pair: |Q_p| (rho_{x,y}-cycle through p)
      and whether it avoids p'.  A witness at threshold l' is |Q_p| <= l' and p' not in Q_p.  Reports, for each l',
      P[no witness among the first K pairs | K0-1 >= K] against the independent prediction (1-w1)^K, the hazard
      P[witness at k | none before k], and P[Flip = empty | K0-1 >= K] against 2^-K.
 (ii) DIAGONALS: on B-instances with l = d12 <= n/4 and |C1| >= n/2 + l, the diagonals D_c, 0 <= c <= n/2 - l,
      each the l-1 pairs (x_{c+i}, y_i), i = 1..l-1 (x-row c steps deeper).  G = number of diagonals containing a
      witness (threshold l'), Y = number of witness pairs.  Reports per l: number of instances, mean per-pair witness
      rate, E[G]/#diag, P[G = 0] and the lower tail P[G < G_mean/2] against the binomial prediction with the same
      per-diagonal rate (diagonals share their y-rows, so correlations are possible).
usage: ./jmfix ... | python3 s13_wtail.py n [--marks M] [--seed s] [--classes m:alpha,...] [--lmax Lmax]
"""
import sys, random, math
from collections import defaultdict, Counter
from s8_real_bias import read_squares
from pf_cl import cycles_of

THR = (4, 8, 16, 24)


def main():
    n = int(sys.argv[1]); args = sys.argv[2:]
    def opt(name, default, conv=int):
        return conv(args[args.index(name) + 1]) if name in args else default
    M = opt('--marks', 6); seed = opt('--seed', 7); LMAX = opt('--lmax', n // 4)
    classes = opt('--classes', None, str)
    classes = None if classes is None else set(tuple(int(t) for t in c.split(':')) for c in classes.split(','))
    rng = random.Random(seed)
    sig = None; marks = []
    # ladder accumulators (pooled over classes; per-class first-pair rates also kept)
    NK = Counter()                      # NK[K] = instances with K0-1 >= K
    nowit = {t: Counter() for t in THR}  # nowit[t][K] = instances with K0-1 >= K and no witness among first K
    noflip = Counter()
    hz_num = {t: Counter() for t in THR}; hz_den = {t: Counter() for t in THR}
    w1 = {t: [0, 0] for t in THR}
    # diagonal accumulators
    diag = defaultdict(lambda: {t: [] for t in THR})   # diag[l][t] = list of (ndiag, npairs, G, Y)
    nsq = 0; ninst = 0
    for L in read_squares(n, sys.stdin.buffer):
        nsq += 1
        if sig is None:
            assert L[0] == list(range(n)); sig = L[1][:]
            for cyc in cycles_of(sig):
                m = len(cyc)
                if m < 4: continue
                for ji, j in enumerate(cyc):
                    for alpha in range(2, m - 1):
                        if classes is None or (m, alpha) in classes:
                            marks.append((j, cyc[(ji + alpha) % m]))
            if not marks: print('no admissible marks'); return
        pos = [[0] * n for _ in range(n)]
        for r in range(n):
            for q in range(n): pos[r][L[r][q]] = q
        for (p, pp) in rng.sample(marks, min(M, len(marks))):
            colpp = [0] * n
            for r in range(n): colpp[L[r][pp]] = r
            pi = [colpp[L[r][p]] for r in range(n)]
            inv = [0] * n
            for r in range(n): inv[pi[r]] = r
            c1 = 0; r = 0; d12 = None
            while True:
                r = pi[r]; c1 += 1
                if r == 1: d12 = c1
                if r == 0: break
            together = d12 is not None
            if together:
                d21 = c1 - d12; K0 = min(d12, d21)
            else:
                c2 = 0; r = 1
                while True:
                    r = pi[r]; c2 += 1
                    if r == 1: break
                K0 = min(c1, c2)
            if K0 < 2: continue
            ninst += 1
            # chains of preimages
            xs = [0]; ys = [1]
            def qp_info(x, y):
                """(|Q_p|, avoids p') for rho_{x,y}."""
                rowy = pos[y]; ln = 1; q = rowy[L[x][p]]; hit = False
                while q != p:
                    if q == pp: hit = True
                    q = rowy[L[x][q]]; ln += 1
                return ln, (not hit)
            # (i) ladder
            x, y = 0, 1; first_wit = {t: None for t in THR}; first_flip = None
            for k in range(1, K0):
                x, y = inv[x], inv[y]
                ln, av = qp_info(x, y)
                if av and first_flip is None: first_flip = k
                for t in THR:
                    if av and ln <= t and first_wit[t] is None: first_wit[t] = k
                    if k == 1:
                        w1[t][0] += (av and ln <= t); w1[t][1] += 1
            Kmax = K0 - 1
            for K in range(1, Kmax + 1):
                NK[K] += 1
                if first_flip is None or first_flip > K: noflip[K] += 1
                for t in THR:
                    if first_wit[t] is None or first_wit[t] > K: nowit[t][K] += 1
            for t in THR:
                fw = first_wit[t]
                for k in range(1, Kmax + 1):
                    if fw is not None and k > fw: break
                    hz_den[t][k] += 1
                    if fw == k: hz_num[t][k] += 1
            # (ii) diagonals
            if together and 2 <= d12 <= LMAX and c1 >= n // 2 + d12:
                l = d12; ncs = n // 2 - l + 1
                # x-chain to depth ncs-1 + l-1, y-chain to depth l-1
                xch = [0]; r = 0
                for _ in range(ncs + l - 1): r = inv[r]; xch.append(r)
                ych = [1]; r = 1
                for _ in range(l - 1): r = inv[r]; ych.append(r)
                G = {t: 0 for t in THR}; Y = {t: 0 for t in THR}
                for c in range(ncs):
                    has = {t: False for t in THR}
                    for i in range(1, l):
                        ln, av = qp_info(xch[c + i], ych[i])
                        if av:
                            for t in THR:
                                if ln <= t: has[t] = True; Y[t] += 1
                    for t in THR: G[t] += has[t]
                for t in THR: diag[l][t].append((ncs, ncs * (l - 1), G[t], Y[t]))
    print(f'n={n}: squares={nsq} instances={ninst} type={sorted((len(c) for c in cycles_of(sig)), reverse=True)} marks/square={M} classes={sorted(classes) if classes else "all"}')
    print('(i) LADDER.  first-pair witness rate w1:', {t: round(w1[t][0] / max(1, w1[t][1]), 4) for t in THR})
    print('   K : N(K0-1>=K) | P[Flip=0] vs 2^-K | for l\'=4,8,16,24: P[no witness in first K] vs (1-w1)^K')
    for K in (1, 2, 3, 4, 6, 8, 10, 12, 15, 20, 25, 30, 40):
        if NK[K] < 200: continue
        s = f'   {K:3d}: {NK[K]:7d} | {noflip[K]/NK[K]:.5f} vs {2.0**-K:.5f} |'
        for t in THR:
            r = w1[t][0] / max(1, w1[t][1])
            s += f'  {nowit[t][K]/NK[K]:.4f} vs {(1-r)**K:.4f}'
        print(s)
    print('   hazard P[witness at k | none before k], k=1..12:')
    for t in THR:
        print(f'     l\'={t:2d}:', ' '.join(f'{hz_num[t][k]/hz_den[t][k]:.3f}' if hz_den[t][k] >= 300 else '  -  ' for k in range(1, 13)), f'  (den at k=12: {hz_den[t][12]})')
    print('(ii) DIAGONALS (B, l=d12, |C1|>=n/2+l).  per l: N; for each l\': per-pair rate, E[G]/#diag, P[G=0] vs binomial, P[G<E[G]/2] vs binomial')
    for l in sorted(diag):
        rows = diag[l][THR[0]]; N = len(rows)
        if N < 100: continue
        s = f'   l={l:2d} N={N:5d} #diag={rows[0][0]:3d}:'
        for t in THR:
            rows = diag[l][t]
            pr = sum(Y for _, np_, G, Y in rows) / sum(np_ for _, np_, G, Y in rows)
            gd = sum(G for *_, G, Y in rows) / sum(nd for nd, *_ in rows)
            nd = rows[0][0]
            p0 = sum(1 for *_, G, Y in rows if G == 0) / N
            b0 = (1 - gd) ** nd
            half = gd * nd / 2
            plow = sum(1 for *_, G, Y in rows if G < half) / N
            blow = sum(math.comb(nd, g) * gd ** g * (1 - gd) ** (nd - g) for g in range(0, nd + 1) if g < half)
            s += f' | l\'={t}: {pr:.3f} {gd:.3f} P0={p0:.4f}/{b0:.4f} low={plow:.3f}/{blow:.3f}'
        print(s)


if __name__ == '__main__':
    main()
