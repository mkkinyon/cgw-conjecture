"""s13_ladder.py -- the ladder identity on real squares (session 13, sec:s13ladder).

Ladder pairs: (x_k, y_k) = (pi^-k(1), pi^-k(2)), 1 <= k < K0, K0 = min(|C1|,|C2|)
(apart) or min(d12, d21) (together).  The j'-trades at these pairs commute, preserve
R and the mark, and flip the bit iff the pair is flippable, so in every orbit of the
group they generate the bit is exactly fair unless NO ladder pair is flippable:
    P[B|R] - P[A|R] = P[B, Flip=0 | R] - P[A, Flip=0 | R].
This script measures P[Flip=0], the law of K0, and compares with E[2^-(K0-1)].
usage: ./jm n nsamp thin seed | python3 s13_ladder.py n [--inst I] [--seed s]
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
    ninst = opt("--inst", 4); seed = opt("--seed", 9)
    rng = random.Random(seed)
    nsq = 0; N = 0; nB = nA = 0; noneB = noneA = 0; pred = 0.0
    k0hist = defaultdict(lambda: [0, 0]); firstflip = defaultdict(int)
    coinpos = [0, 0]  # ladder coins: count, flippable (all k < K0)
    lenhist = defaultdict(int); big = max(4, n // 4); nbig = 0
    wit = {2: 0, 3: 0, 5: 0, 10: 0, 20: 0, 50: 0, 100: 0}
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
            def jcycle(x, y):
                """length of the rho_{x,y}-cycle through column j, and whether it contains jp"""
                q = j; ln = 0; hasjp = False
                while True:
                    q = pos[y][L[x][q]]; ln += 1
                    if q == jp: hasjp = True
                    if q == j: return ln, hasjp
            def flippable(x, y):
                return not jcycle(x, y)[1]
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
            x, y = 0, 1; anyflip = 0; ff = 0; minlen = n + 1
            for k in range(1, K0):
                x, y = inv[x], inv[y]
                ln, hasjp = jcycle(x, y); f = not hasjp; coinpos[0] += 1; coinpos[1] += f
                if f and not anyflip: anyflip = k
                if f and ln < minlen: minlen = ln
                lenhist[min(ln, n)] += 1
            if K0 >= big:
                nbig += 1
                for l in wit:
                    if minlen <= l: wit[l] += 1
            N += 1; pred += 2.0 ** (-(K0 - 1))
            h = k0hist[min(K0, 40)]; h[0] += 1; h[1] += (anyflip == 0)
            if anyflip: firstflip[anyflip] += 1
            if together:
                nB += 1; noneB += (anyflip == 0)
            else:
                nA += 1; noneA += (anyflip == 0)
    print(f"n={n}: squares={nsq} instances={N}")
    print(f"  P[B]-P[A] = {(nB-nA)/N:+.4f} (+-{1/math.sqrt(N):.4f});  P[B,Flip=0]-P[A,Flip=0] = {(noneB-noneA)/N:+.4f};  "
          f"P[Flip=0] = {(noneB+noneA)/N:.4f};  E[2^-(K0-1)] = {pred/N:.4f}")
    print(f"  ladder coins flippable: {coinpos[1]/max(1,coinpos[0]):.4f} of {coinpos[0]}")
    print("  K0 : count  P[Flip=0 | K0]  2^-(K0-1)")
    for k in sorted(k0hist):
        c, z = k0hist[k]
        print(f"  {k:3d}{'+' if k == 40 else ' '}: {c:5d}   {z/c:.4f}         {2.0**(-(k-1)):.4f}")
    tot = sum(lenhist.values())
    print("  ladder j-cycle length <= l (fraction of pairs): " + ", ".join(f"{l}:{sum(v for k,v in lenhist.items() if k<=l)/tot:.4f}" for l in (1,2,3,5,10,20,50,100) if l <= n))
    print(f"  instances with K0 >= {big}: {nbig}; of these, fraction having a flippable ladder pair whose j-cycle has length <= l: "
          + ", ".join(f"{l}:{wit[l]/max(1,nbig):.4f}" for l in sorted(wit) if l <= n))
    print("  first flippable ladder index k: " + ", ".join(f"{k}:{firstflip[k]}" for k in sorted(firstflip)[:12]))


if __name__ == "__main__":
    main()
