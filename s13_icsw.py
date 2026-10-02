"""s13_icsw.py -- intercalate switches on a ladder row (session 13, sec:crossspace).

For a ladder pair (x,y) of a uniform square with a true mark (p,p') = (j,jp):
  Q_p = rho_{x,y}-cycle through p.  Crossed: p' in Q_p, arcs P1 (p->p') and P2 (p'->p)
  of interior lengths d1-1, d2-1.  Parallel: Q_p, Q_p' separate, interior sizes l-1, b-1.
An intercalate on rows (x,z), columns (q,c) with q,c not in {p,p'} and z not in {1,2}
is a legal move preserving rows 1,2 and columns p,p'; it right-multiplies rho_{x,y}
by (q c): it SEPARATES p,p' iff crossed and q in P1, c in P2 (or vice versa), and
MERGES iff parallel and q in Q_p\\p, c in Q_p'\\p'.  The switch is an involution, so
   E[I_sep ; crossed] = E[I_merge ; parallel]  exactly (conditionally on the frame and rows 1,2).
This script measures I_sep, I_merge against (d1-1)(d2-1)/n and (l-1)(b-1)/n.
usage: ./jm n nsamp thin seed | python3 s13_icsw.py n [--inst I] [--seed s]
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
    ninst = opt("--inst", 4); seed = opt("--seed", 17)
    rng = random.Random(seed)
    nsq = 0; N = 0
    acc = {'cr': [0, 0.0, 0.0, 0.0], 'par': [0, 0.0, 0.0, 0.0]}  # count, sum I, sum product/n, sum I*n/product
    bins = defaultdict(lambda: [0, 0.0, 0.0])  # (state, product-bin) -> count, sum I, sum product/n
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
        colpos = [[0] * n for _ in range(n)]   # colpos[q][s] = row holding s in column q
        for r in range(n):
            for q in range(n): colpos[q][L[r][q]] = r
        for (p, pp) in rng.sample(pairs, min(ninst, len(pairs))):
            pi = pi_P(L, p, pp); inv = [0] * n
            for r in range(n): inv[pi[r]] = r
            # ladder pairs
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
            # rho_{x,y}-cycle through p
            def rho(q): return pos[y][L[x][q]]
            Qp = [p]; q = rho(p)
            while q != p: Qp.append(q); q = rho(q)
            if pp in Qp:
                i = Qp.index(pp); P1 = Qp[1:i]; P2 = Qp[i + 1:]
                A, B = P1, P2; state = 'cr'
            else:
                Qpp = [pp]; q = rho(pp)
                while q != pp: Qpp.append(q); q = rho(q)
                A, B = Qp[1:], Qpp[1:]; state = 'par'
            # intercalates of row x with q in A, c in B: z = row holding L[x][q] in column c; need L[z][q] == L[x][c], z not in {0,1}
            I = 0
            Bset = set(B)
            for q in A:
                for c in B:
                    z = colpos[c][L[x][q]]
                    if z != x and z not in (0, 1) and L[z][q] == L[x][c]: I += 1
            prod = len(A) * len(B)
            a = acc[state]; a[0] += 1; a[1] += I; a[2] += prod / n
            if prod: a[3] += I * n / prod
            pb = 0 if prod == 0 else min(6, int(prod * 16 / (n * n)))   # bins of product/n^2 in 1/16 steps
            b = bins[(state, pb)]; b[0] += 1; b[1] += I; b[2] += prod / n
            N += 1
    print(f"n={n}: squares={nsq} ladder-pair instances={N}")
    for state in ('cr', 'par'):
        a = acc[state]
        if a[0]: print(f"  {state}: instances={a[0]}  E[I]={a[1]/a[0]:.3f}  E[product/n]={a[2]/a[0]:.3f}  ratio={a[1]/max(1e-9,a[2]):.3f}")
    print("  product/n^2 bin : state  count  E[I]  E[product/n]  ratio")
    for key in sorted(bins):
        c, sI, sp = bins[key]
        if c >= 10: print(f"   [{key[1]/16:.3f},{(key[1]+1)/16:.3f}) : {key[0]:3s} {c:6d}  {sI/c:7.3f}  {sp/c:7.3f}   {sI/max(1e-9,sp):.3f}")


if __name__ == "__main__":
    main()
