"""s13_witness.py -- scoping the witness lemma (hyp:witness) on real squares.

witness(x,y; l) = [the rho_{x,y}-cycle through column j has length <= l and avoids jp].
For each instance (uniform square, true mark) record
  (i)  ladder: K0, G_l = # ladder pairs with a witness (l = 5, 10, 20)   -> Var(G | K0) vs binomial;
  (ii) star:   for y1 = pi^-1(2) and all x in C1 \\ {1} (apart) or on the arc 2->1 (together):
               S_l = # x with witness(x, y1)                            -> Var(S | size) vs binomial;
  (iii) pairwise products of witness indicators for disjoint ladder pairs vs product of marginals.
usage: ./jm n nsamp thin seed | python3 s13_witness.py n [--inst I] [--seed s]
"""
import sys, random, math
from collections import defaultdict
from verify_identity import sigma_of
from pf_cl import cycles_of
from pf_master import pi_P
from s8_real_bias import read_squares

LS = (5, 10, 20)


def main():
    n = int(sys.argv[1]); args = sys.argv[2:]
    def opt(name, default):
        return int(args[args.index(name) + 1]) if name in args else default
    ninst = opt("--inst", 4); seed = opt("--seed", 13)
    rng = random.Random(seed)
    nsq = 0; N = 0
    lad = {l: defaultdict(lambda: [0, 0.0, 0.0]) for l in LS}   # l -> K0bin -> [count, sum G, sum G^2]
    star = {l: defaultdict(lambda: [0, 0.0, 0.0]) for l in LS}  # l -> sizebin -> [count, sum S, sum S^2]
    marg = {l: [0, 0] for l in LS}                               # single-pair witness frequency (ladder pairs)
    pair = {l: [0, 0.0] for l in LS}                             # products of indicators, disjoint ladder pairs
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
        for (j, jp) in rng.sample(pairs, min(ninst, len(pairs))):
            pi = pi_P(L, j, jp); inv = [0] * n
            for r in range(n): inv[pi[r]] = r
            def jlen(x, y):
                q = j; ln = 0; hasjp = False
                while True:
                    q = pos[y][L[x][q]]; ln += 1
                    if q == jp: hasjp = True
                    if q == j: return (n + 1) if hasjp else ln
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
            # (i) ladder
            x, y = 0, 1; lens = []
            for k in range(1, K0):
                x, y = inv[x], inv[y]; lens.append(jlen(x, y))
            for l in LS:
                ind = [1 if ln <= l else 0 for ln in lens]
                G = sum(ind); b = lad[l][min(K0 // 8, 6)]; b[0] += 1; b[1] += G; b[2] += G * G
                marg[l][0] += len(ind); marg[l][1] += G
                for _ in range(min(30, len(ind) * (len(ind) - 1) // 2)):
                    a, c = rng.sample(range(len(ind)), 2)
                    pair[l][0] += 1; pair[l][1] += ind[a] * ind[c]
            # (ii) star around y1 = pi^-1(2): x over the rows preceding 1 (C1 \ {1} apart; arc 2->1 together)
            y1 = inv[1]; xs = []; x = 0
            while True:
                x = inv[x]
                if x == 1 or x == 0: break
                xs.append(x)
            if xs:
                lens = [jlen(x, y1) for x in xs]
                for l in LS:
                    S = sum(1 for ln in lens if ln <= l); b = star[l][min(len(xs) // 8, 6)]; b[0] += 1; b[1] += S; b[2] += S * S
            N += 1
    print(f"n={n}: squares={nsq} instances={N}")
    for l in LS:
        p = marg[l][1] / max(1, marg[l][0])
        print(f" l={l:2d}: single-pair witness freq p={p:.4f}  ((l-1)/(n-1)={(l-1)/(n-1):.4f});  "
              f"disjoint ladder pairs E[W W']/p^2 = {pair[l][1]/max(1,pair[l][0])/p**2:.3f} ({pair[l][0]} samples)")
        print("   ladder G:  K0-bin   count   mean G   var G   binomial var   ratio")
        for kb in sorted(lad[l]):
            c, s1, s2 = lad[l][kb]
            if c < 20: continue
            mean = s1 / c; var = s2 / c - mean * mean; bv = mean * (1 - p)
            print(f"              {8*kb:3d}+    {c:5d}   {mean:6.3f}  {var:6.3f}   {bv:6.3f}       {var/max(bv,1e-9):.2f}")
        print("   star  S:  size-bin count   mean S   var S   binomial var   ratio")
        for kb in sorted(star[l]):
            c, s1, s2 = star[l][kb]
            if c < 20: continue
            mean = s1 / c; var = s2 / c - mean * mean; bv = mean * (1 - p)
            print(f"              {8*kb:3d}+    {c:5d}   {mean:6.3f}  {var:6.3f}   {bv:6.3f}       {var/max(bv,1e-9):.2f}")


if __name__ == "__main__":
    main()
