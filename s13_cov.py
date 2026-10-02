"""s13_cov.py -- concentration input (gamma) of hyp:orbit: covariance of the orbit weights ybar_k across a sub-ladder.

For each instance (uniform square, true mark), the sub-ladder I = first min(K0-1, n/8) ladder pairs, a FIXED matching of
the free columns (consecutive free columns), and for each k in I the clean coordinates T_k (column cycles through x_k of
length <= LMAX meeting {1,2} in 0 or 2 rows and containing no other row of I u {y_k}).  ybar_k = E_eta[fraction of T_k
separated by {p,p'} in rho_k(eta)] (NEPS samples).  Reports: mean of |I|, E[ybar], Var(sum_k ybar_k) vs sum_k Var(ybar_k)
(the ratio is 1 + (|I|-1) x average correlation), and the fraction of instances with sum_k ybar_k < kappa |I| for a few kappa.
usage: ./jm n nsamp thin seed | python3 s13_cov.py n [--inst I] [--seed s] [--lmax L] [--neps E]
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
    ninst = opt("--inst", 2); seed = opt("--seed", 29); LMAX = opt("--lmax", 6); NEPS = opt("--neps", 32)
    rng = random.Random(seed)
    nsq = 0; N = 0
    S1 = 0.0; S2 = 0.0; sumI = 0; sumvar = 0.0; sumy = 0.0; sumy2 = 0.0; npairs = 0
    below = defaultdict(int); kappas = (0.05, 0.1, 0.15)
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
            m = min(K0 - 1, n // 8)
            if m < 2: continue
            lad = []; x, y = 0, 1
            for k in range(m): x, y = inv[x], inv[y]; lad.append((x, y))
            ladrows = set(r for xy in lad for r in xy)
            free = [q for q in range(n) if q not in (p, pp)]
            matching = [(free[a], free[a + 1]) for a in range(0, len(free) - 1, 2)]
            ys = []
            for (x, y) in lad:
                rho = [pos[y][L[x][q]] for q in range(n)]
                T = []
                for (q, c) in matching:
                    Z = [x]; r = colpos[c][L[x][q]]
                    while r != x and len(Z) <= LMAX: Z.append(r); r = colpos[c][L[r][q]]
                    if r != x: continue
                    m12 = (0 in Z) + (1 in Z)
                    if m12 == 1: continue
                    if any(z in ladrows for z in Z[1:]): continue
                    T.append((q, c))
                if not T: ys.append(0.0); continue
                def sides(rh):
                    side = [0] * n
                    Qp = [p]; q = rh[p]
                    while q != p: Qp.append(q); q = rh[q]
                    if pp in Qp:
                        i = Qp.index(pp)
                        for q in Qp[1:i]: side[q] = 1
                        for q in Qp[i + 1:]: side[q] = 2
                        return side
                    for q in Qp[1:]: side[q] = 1
                    q = rh[pp]
                    while q != pp: side[q] = 2; q = rh[q]
                    return side
                tot = 0.0
                for _ in range(NEPS):
                    rh = list(rho)
                    for (q, c) in T:
                        if rng.random() < 0.5: rh[q], rh[c] = rh[c], rh[q]
                    sd = sides(rh)
                    tot += sum(1 for (q, c) in T if sd[q] and sd[c] and sd[q] != sd[c]) / len(T)
                ys.append(tot / NEPS)
            Y = sum(ys); N += 1; sumI += m; S1 += Y; S2 += Y * Y
            for v in ys: sumy += v; sumy2 += v * v; npairs += 1
            for kap in kappas:
                if Y < kap * m: below[kap] += 1
    EY = S1 / N; VarY = S2 / N - EY * EY
    Ey = sumy / npairs; Vary = sumy2 / npairs - Ey * Ey; mI = sumI / N
    print(f"n={n}: squares={nsq} instances={N} LMAX={LMAX} NEPS={NEPS}  mean |I|={mI:.2f}")
    print(f"  per pair: E[ybar]={Ey:.4f}  Var(ybar)={Vary:.4f}")
    print(f"  per instance: E[Y]={EY:.3f}  Var(Y)={VarY:.3f}  vs independent sum |I| Var(ybar) = {mI*Vary:.3f}  (ratio {VarY/max(1e-9,mI*Vary):.2f}; includes the variation of |I|)")
    for kap in kappas: print(f"  P[Y < {kap}|I|] = {below[kap]/N:.4f}")


if __name__ == "__main__":
    main()
