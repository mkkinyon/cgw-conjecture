"""s13_nutail.py -- R4: the switching for the nu-tail at the first ladder pair, checked on samples.

Event E(start, end, k): the rho-path from start in {p,p'} reaches end in {p,p'} at step k (first return/arrival),
k >= 2; nu <= M is the union over k <= M.  Forward move: for a row u notin {0,1,x1,y1} such that the rho_{x1,u}-cycle
W through w_{k-1} avoids S = {p} u {w_0,...,w_{k-2}}, trade rows x1,u along W.  Claims checked: the result is a Latin
square with rows 0,1 unchanged (so in X) and the same first pair; and the move is inverted from (L', k, start, end)
alone: the rho'-path from start is w_0..w_{k-1} (unchanged), s = L'(y1, end), u = row of s in column w_{k-1} of L',
W = rho'_{x1,u}-cycle through w_{k-1}, trade back.  Also records f(L)/n = fraction of admissible u.
usage: ./jmfix ... | python3 s13_nutail.py n [--marks M] [--kmax K] [--seed s]
"""
import sys, random
from collections import Counter, defaultdict
from s8_real_bias import read_squares
from pf_cl import cycles_of


def main():
    n = int(sys.argv[1]); args = sys.argv[2:]
    def opt(name, default):
        return int(args[args.index(name) + 1]) if name in args else default
    M = opt('--marks', 4); KMAX = opt('--kmax', 6); seed = opt('--seed', 5)
    rng = random.Random(seed)
    sig = None; marks = []; st = Counter(); fhist = defaultdict(list); nsq = 0
    for L in read_squares(n, sys.stdin.buffer):
        nsq += 1
        if sig is None:
            sig = L[1][:]
            for cyc in cycles_of(sig):
                m = len(cyc)
                if m < 4: continue
                for ji, j in enumerate(cyc):
                    for alpha in range(2, m - 1): marks.append((j, cyc[(ji + alpha) % m]))
        pos = [[0] * n for _ in range(n)]
        for r in range(n):
            for q in range(n): pos[r][L[r][q]] = q
        for (p, pp) in rng.sample(marks, min(M, len(marks))):
            colpp = [0] * n
            for r in range(n): colpp[L[r][pp]] = r
            pi = [colpp[L[r][p]] for r in range(n)]
            inv = [0] * n
            for r in range(n): inv[pi[r]] = r
            x1, y1 = inv[0], inv[1]
            rho = [pos[y1][L[x1][q]] for q in range(n)]
            # the four events: first arrival at an end in {p,pp} from a start, at step k <= KMAX
            for start in (p, pp):
                path = [start]; q = start
                for k in range(1, KMAX + 1):
                    q = rho[q]; path.append(q)
                    if q in (p, pp):
                        end = q
                        if k >= 2: check(L, n, pos, p, pp, x1, y1, rho, path, k, start, end, st, fhist)
                        break
    print(f'n={n}: squares={nsq}  marks/square={M}  kmax={KMAX}')
    print('  sources (events with k<=kmax):', st['sources'], ' moves checked:', st['moves'])
    print('  Latin & rows 0,1 fixed:', st['latin'], ' first pair unchanged:', st['pair'], ' inverted exactly:', st['inverted'])
    print('  f(L)/n by k (mean, min) and 2^-(k):')
    for k in sorted(fhist):
        v = sorted(fhist[k]); N=len(v); print(f'    k={k}: N={N} mean={sum(v)/N:.3f} min={v[0]:.3f} 1%={v[N//100]:.3f} 10%={v[N//10]:.3f}  P[f<n/20]={sum(1 for t in v if t<0.05)/N:.4f}  ref 2^-k={2.0**-k:.3f}')


def cycle_through(L, n, x, u, c):
    """rho_{x,u}-cycle through column c: rho_{x,u}(q) = column of L[x][q] in row u."""
    posu = {L[u][q]: q for q in range(n)}
    W = [c]; q = posu[L[x][c]]
    while q != c: W.append(q); q = posu[L[x][q]]
    return W


def check(L, n, pos, p, pp, x1, y1, rho, path, k, start, end, st, fhist):
    st['sources'] += 1
    S = set([p] + path[:k - 1])           # {p} u {w_0..w_{k-2}}
    wk1 = path[k - 1]
    adm = []
    for u in range(2, n):
        if u in (x1, y1): continue
        W = cycle_through(L, n, x1, u, wk1)
        if S.isdisjoint(W): adm.append((u, W))
    fhist[k].append(len(adm) / n)
    # check up to 3 admissible moves per source
    for (u, W) in adm[:3]:
        st['moves'] += 1
        Lp = [row[:] for row in L]
        for w in W: Lp[x1][w], Lp[u][w] = Lp[u][w], Lp[x1][w]
        ok = all(sorted(row) == list(range(n)) for row in Lp) and all(sorted(Lp[r][c] for r in range(n)) == list(range(n)) for c in range(n)) and Lp[0] == L[0] and Lp[1] == L[1]
        st['latin'] += ok
        colpp = [0] * n
        for r in range(n): colpp[Lp[r][pp]] = r
        pi2 = [colpp[Lp[r][p]] for r in range(n)]
        inv2 = [0] * n
        for r in range(n): inv2[pi2[r]] = r
        st['pair'] += (inv2[0] == x1 and inv2[1] == y1)
        # inversion from (Lp, k, start, end)
        posy = {Lp[y1][q]: q for q in range(n)}
        rho2 = [posy[Lp[x1][q]] for q in range(n)]
        q = start; path2 = [start]
        for _ in range(k - 1): q = rho2[q]; path2.append(q)
        s = Lp[y1][end]
        u2 = next(r for r in range(n) if Lp[r][path2[k - 1]] == s)
        W2 = cycle_through(Lp, n, x1, u2, path2[k - 1])
        Lq = [row[:] for row in Lp]
        for w in W2: Lq[x1][w], Lq[u2][w] = Lq[u2][w], Lq[x1][w]
        st['inverted'] += (Lq == L and u2 == u and set(W2) == set(W))


if __name__ == '__main__':
    main()
