"""s13_f1.py -- (F1) statistics on real squares (session 13, sec:s13comp7).

For uniform squares (jm.c) and true marks P={j,jp} (jp = sigma^alpha(j), 2<=alpha<=m-2):
  frame pi, cycles C1 (row 0), C2 (row 1), bit = [C1 == C2].
  bit-changing pairs: together -> separating pairs (x on arc 0->1, y on arc 1->0);
                      apart    -> C1 x C2 (rows other than 0,1).
  flippable(x,y) <=> a_x = L[x][jp] and b_x = L[x][j] lie in different cycles of
                     phi_{x,y}: s -> L[x][pos_y(s)]  (equivalently j not on the
                     rho_{x,y}-cycle through jp).
Records, per instance: state, m = #bit-changing pairs, f = #flippable among them,
and pairwise products of indicators for pairs sharing a row / disjoint, plus the
flippable fraction of the other (non-bit-changing) non-adjacent pairs inside
C1 u C2 and of pairs meeting a free cycle.
usage: ./jm n nsamp thin seed | python3 s13_f1.py n [--inst I] [--seed s] [--csv f]
"""
import sys, random, math
from collections import defaultdict
from verify_identity import sigma_of
from pf_cl import cycles_of
from pf_master import pi_P
from s8_real_bias import read_squares


def main():
    n = int(sys.argv[1]); args = sys.argv[2:]
    def opt(name, default):
        return int(args[args.index(name) + 1]) if name in args else default
    ninst = opt("--inst", 3); seed = opt("--seed", 7)
    csvf = open(args[args.index("--csv") + 1], "w") if "--csv" in args else None
    if csvf: csvf.write("state,l1,l2,m,f,f_other,m_other,f_free,m_free\n")
    rng = random.Random(seed)
    nsq = 0; ninst_done = 0
    # aggregates
    tot = defaultdict(lambda: [0, 0, 0, 0])      # state -> [sum m, sum f, count, count f==0]
    binsum = defaultdict(lambda: [0, 0.0, 0.0, 0.0, 0])  # (state, mbin) -> [count, sum f/m, sum (f-m/2)^2, sum m/4, count f==0]
    pair_corr = defaultdict(lambda: [0.0, 0])  # (state, kind) -> [sum (I-1/2)(I'-1/2), count]
    other = defaultdict(lambda: [0, 0]); free = defaultdict(lambda: [0, 0])
    for L in read_squares(n, sys.stdin.buffer):
        nsq += 1
        sig = sigma_of(L); pairs = []
        col = {L[0][c]: c for c in range(n)}   # symbol -> column in row 0 (jm squares are NOT normalised)
        for cyc in cycles_of(sig):
            mlen = len(cyc)
            if mlen < 4: continue
            for ji, jsym in enumerate(cyc):
                for alpha in range(2, mlen - 1): pairs.append((col[jsym], col[cyc[(ji + alpha) % mlen]]))
        if not pairs: continue
        pos = [[0] * n for _ in range(n)]
        for r in range(n):
            for q in range(n): pos[r][L[r][q]] = q
        for (j, jp) in rng.sample(pairs, min(ninst, len(pairs))):
            pi = pi_P(L, j, jp)
            cyc = cycles_of(pi); cid = [0] * n
            for i, c in enumerate(cyc):
                for x in c: cid[x] = i
            def flippable(x, y):
                ax, bx = L[x][jp], L[x][j]
                s = ax
                while True:
                    s = L[x][pos[y][s]]
                    if s == bx: return False
                    if s == ax: return True
            together = cid[0] == cid[1]
            if together:
                C = cyc[cid[0]]; k = len(C); i0 = C.index(0); i1 = C.index(1)
                arcA = [C[(i0 + t) % k] for t in range(1, (i1 - i0) % k)]   # strictly between 0 and 1
                arcB = [C[(i1 + t) % k] for t in range(1, (i0 - i1) % k)]   # strictly between 1 and 0
                bc = [(x, y) for x in arcA for y in arcB]
                l1, l2 = len(arcA), len(arcB)
                others = [(x, y) for x in arcA for y in arcA if x < y] + [(x, y) for x in arcB for y in arcB if x < y]
                state = "tog"
            else:
                A = [x for x in cyc[cid[0]] if x != 0]; B = [x for x in cyc[cid[1]] if x != 1]
                bc = [(x, y) for x in A for y in B]; l1, l2 = len(A), len(B)
                others = [(x, y) for x in A for y in A if x < y] + [(x, y) for x in B for y in B if x < y]
                state = "apt"
            others = [(x, y) for (x, y) in others if pi[x] != y and pi[y] != x]
            ind = {e: int(flippable(*e)) for e in bc}
            m = len(bc); f = sum(ind.values())
            t = tot[state]; t[0] += m; t[1] += f; t[2] += 1; t[3] += (f == 0)
            mb = 0 if m == 0 else min(8, int(math.log2(m)) + 1)
            b = binsum[(state, mb)]; b[0] += 1
            if m: b[1] += f / m; b[2] += (f - m / 2) ** 2; b[3] += m / 4; b[4] += (f == 0)
            # pairwise correlations (subsample)
            es = list(bc)
            for _ in range(min(400, len(es) * (len(es) - 1) // 2)):
                e1, e2 = rng.sample(es, 2)
                kind = "share" if set(e1) & set(e2) else "disj"
                pc = pair_corr[(state, kind)]; pc[0] += (ind[e1] - 0.5) * (ind[e2] - 0.5); pc[1] += 1
            fo = sum(flippable(*e) for e in others); other[state][0] += len(others); other[state][1] += fo
            # pairs meeting a free cycle
            freerows = [x for x in range(2, n) if cid[x] != cid[0] and cid[x] != cid[1]]
            ff = mf = 0
            if freerows:
                for _ in range(40):
                    x = rng.choice(freerows); y = rng.randrange(2, n)
                    if y == x or pi[x] == y or pi[y] == x: continue
                    mf += 1; ff += flippable(x, y)
            free[state][0] += mf; free[state][1] += ff
            if csvf: csvf.write(f"{state},{l1},{l2},{m},{f},{fo},{len(others)},{ff},{mf}\n")
            ninst_done += 1
    if csvf: csvf.close()
    print(f"n={n}: squares={nsq} instances={ninst_done}")
    for state in ("tog", "apt"):
        t = tot[state]
        if not t[2]: continue
        print(f" {state}: instances={t[2]}  mean #bit-changing pairs m={t[0]/t[2]:.1f}  overall flippable fraction={t[1]/max(1,t[0]):.4f}  P[f=0]={t[3]/t[2]:.4f}")
        o = other[state]; fr = free[state]
        print(f"      other non-adjacent pairs inside C1uC2: {o[1]/max(1,o[0]):.4f} ({o[0]});  pairs meeting a free cycle: {fr[1]/max(1,fr[0]):.4f} ({fr[0]})")
        for kind in ("share", "disj"):
            pc = pair_corr[(state, kind)]
            if pc[1]: print(f"      corr of indicators, {kind:5s} pairs: {4*pc[0]/pc[1]:+.4f}  ({pc[1]} samples)")
        print("      m-bin      count   mean f/m   var(f)/(m/4)   P[f=0]   2^-m(pred, bin mid)")
        for mb in range(0, 9):
            b = binsum.get((state, mb))
            if not b: continue
            lo = 0 if mb == 0 else 2 ** (mb - 1); hi = 0 if mb == 0 else 2 ** mb - 1
            print(f"      [{lo:3d},{hi:3d}]  {b[0]:6d}   {b[1]/b[0] if lo else float('nan'):.4f}    {b[2]/max(1e-9,b[3]):.3f}        {b[4]/b[0]:.4f}   {2**-((lo+hi)/2) if lo else 1:.4f}")


if __name__ == "__main__":
    main()
