"""s13_fixstats.py -- statistics of marked instances on completions of a FIXED 2 x n rectangle (rows 0,1 = (id, sigma)).

Reads squares (n*n bytes each) from stdin (output of ./jmfix).  For every square and every sampled mark
(p, p') with p' = sigma^alpha(p) in a sigma-cycle of length m (2 <= alpha <= m-2), computes the frame pi
(L[pi(r)][p'] = L[r][p]), B = [rows 0,1 in a common pi-cycle], K0, the ladder (x_k,y_k) = (pi^-k(0), pi^-k(1)),
k = 1..K0-1, flippability bits (p, p' in different cycles of rho_{x,y}), and for the first ladder pair the
length of the rho-cycle through p and whether it avoids p' (the witness event at threshold l').
usage: ./jmfix ... | python3 s13_fixstats.py n [--marks M] [--seed s] [--classes m:alpha,m:alpha,...]
Reports per class (m, alpha): N, P[B], qbar = P[B]/P[A], P[K0>=2], P[(x1,y1) flippable] (all, |B, |A),
P[Flip_ladder = empty | K0 = k], witness frequencies at pair 1 (cycle through p of length <= L avoiding p').
"""
import sys, random
from collections import defaultdict, Counter
from s8_real_bias import read_squares
from pf_cl import cycles_of


def main():
    n = int(sys.argv[1]); args = sys.argv[2:]
    def opt(name, default, conv=int):
        return conv(args[args.index(name) + 1]) if name in args else default
    M = opt('--marks', 8); seed = opt('--seed', 7); TC = opt('--cand', 0)
    classes = opt('--classes', None, str)
    classes = None if classes is None else set(tuple(int(t) for t in c.split(':')) for c in classes.split(','))
    rng = random.Random(seed)
    Lthr = (4, 8, 16, 32)
    S = defaultdict(Counter)   # per class
    nsq = 0; sig = None; marks_by_class = None
    for L in read_squares(n, sys.stdin.buffer):
        nsq += 1
        if sig is None:
            assert L[0] == list(range(n))
            sig = L[1][:]
            marks_by_class = defaultdict(list)
            for cyc in cycles_of(sig):
                m = len(cyc)
                if m < 4: continue
                for ji, j in enumerate(cyc):
                    for alpha in range(2, m - 1):
                        if classes is None or (m, alpha) in classes:
                            marks_by_class[(m, alpha)].append((j, cyc[(ji + alpha) % m]))
            if not marks_by_class:
                print('no admissible mark classes'); return
        assert L[0] == list(range(n)) and L[1] == sig, 'rows 0,1 changed!'
        pos = [[0] * n for _ in range(n)]
        for r in range(n):
            for q in range(n): pos[r][L[r][q]] = q
        for cls, marks in marks_by_class.items():
            st = S[cls]
            for (p, pp) in rng.sample(marks, min(M, len(marks))):
                # frame
                colpp = [0] * n
                for r in range(n): colpp[L[r][pp]] = r
                pi = [colpp[L[r][p]] for r in range(n)]
                inv = [0] * n
                for r in range(n): inv[pi[r]] = r
                # cycle of 0, is 1 on it, distances
                c1 = 0; r = 0; d12 = None
                while True:
                    r = pi[r]; c1 += 1
                    if r == 1: d12 = c1
                    if r == 0: break
                together = d12 is not None
                if together: K0 = min(d12, c1 - d12)
                else:
                    c2 = 0; r = 1
                    while True:
                        r = pi[r]; c2 += 1
                        if r == 1: break
                    K0 = min(c1, c2)
                st['N'] += 1; st['B'] += together
                if K0 < 2: continue
                st['K0>=2'] += 1; st[('K0', min(K0, 10))] += 1
                x, y = inv[0], inv[1]
                anyflip = False; first = None
                for k in range(1, K0):
                    rho = [pos[y][L[x][q]] for q in range(n)]
                    # cycle through p
                    ln = 1; q = rho[p]; hit = False
                    while q != p:
                        if q == pp: hit = True
                        q = rho[q]; ln += 1
                    flip = not hit
                    if k == 1:
                        first = (flip, ln)
                    anyflip = anyflip or flip
                    x, y = inv[x], inv[y]
                flip1, ln1 = first
                if TC:
                    # adaptive candidates at the first ladder pair: (q_t,c_t) = (rho^t(p), rho^t(p')), t < mu;
                    # good_t: x1, y1 in different cycles of the column permutation tau_{q,c} (tau(r) = row of L[r][q]
                    # in column c) and the cycle through x1 or the one through y1 meets {0,1} in 0 or 2 rows.
                    x1, y1 = inv[0], inv[1]
                    rho = [pos[y1][L[x1][q]] for q in range(n)]
                    # mu: crossed -> min arc length; parallel -> min cycle length
                    d = 1; q = rho[p]
                    while q != p and q != pp: q = rho[q]; d += 1
                    if q == pp:
                        e = 1; q = rho[pp]
                        while q != p: q = rho[q]; e += 1
                    else:
                        e = 1; q = rho[pp]
                        while q != pp: q = rho[q]; e += 1
                    mu = min(d, e)
                    st['mu>=2'] += (mu >= 2); st['musum'] += mu
                    q, c = p, pp; tgood = None; g1 = None
                    for t in range(1, min(TC, mu - 1) + 1):
                        q, c = rho[q], rho[c]
                        colc = [0] * n
                        for r in range(n): colc[L[r][c]] = r
                        def cyc_through(r0):
                            Z = [r0]; r = colc[L[r0][q]]
                            while r != r0: Z.append(r); r = colc[L[r][q]]
                            return Z
                        Zx = cyc_through(x1)
                        if y1 in Zx: ok = False
                        else:
                            Zy = cyc_through(y1)
                            ok = ((0 in Zx) + (1 in Zx)) != 1 or ((0 in Zy) + (1 in Zy)) != 1
                        if t == 1: g1 = ok
                        if ok: tgood = t; break
                    st['G1'] += bool(g1); st['anygood'] += tgood is not None
                    st['flip&good'] += flip1 and tgood is not None; st['noflip&good'] += (not flip1) and tgood is not None
                    st['cand>=TC'] += (mu - 1 >= TC)
                st['flip1'] += flip1; st['flip1&B'] += flip1 and together; st['flip1&A'] += flip1 and not together
                st['B&K0>=2'] += together
                st[('noflip', min(K0, 10))] += (not anyflip)
                st['len1sum'] += ln1
                for Lt in Lthr:
                    st[('wit', Lt)] += (ln1 <= Lt and flip1)
                    st[('short', Lt)] += (ln1 <= Lt)
    print(f'n={n}: squares={nsq}  sigma cycle type={sorted((len(c) for c in cycles_of(sig)), reverse=True)}  marks/class/square={M}')
    for cls in sorted(S):
        st = S[cls]; N = st['N']; B = st['B']; A = N - B
        print(f'class m={cls[0]} alpha={cls[1]}: N={N}  P[B]={B/N:.4f}  qbar=P[B]/P[A]={B/max(1,A):.4f}  P[K0>=2]={st["K0>=2"]/N:.4f}')
        K2 = st['K0>=2']; BK = st['B&K0>=2']; AK = K2 - BK
        print(f'   P[(x1,y1) flippable | K0>=2] = {st["flip1"]/max(1,K2):.4f}   | B: {st["flip1&B"]/max(1,BK):.4f}  | A: {st["flip1&A"]/max(1,AK):.4f}   mean |cycle through p| at pair 1 = {st["len1sum"]/max(1,K2):.2f}')
        print('   P[Flip_ladder = empty | K0=k] (k: N, p):', {k: (st[('K0', k)], round(st[('noflip', k)] / st[('K0', k)], 4)) for k in range(2, 11) if st[('K0', k)]})
        print('   witness at pair 1, P[len<=L and avoids p\'] / P[len<=L]:', {Lt: (round(st[('wit', Lt)] / max(1, K2), 4), round(st[('short', Lt)] / max(1, K2), 4)) for Lt in Lthr})
        if TC:
            print(f'   adaptive candidates at pair 1 (t<=min({TC},mu-1)): P[mu>=2]={st["mu>=2"]/K2:.4f} E[mu]={st["musum"]/K2:.2f} P[mu-1>={TC}]={st["cand>=TC"]/K2:.4f}  P[G_1]={st["G1"]/K2:.4f}  P[some good t]={st["anygood"]/K2:.4f}')
            print(f'      involution check: P[flip1 & good]={st["flip&good"]/K2:.4f}  P[not flip1 & good]={st["noflip&good"]/K2:.4f}   bound P[flip1] >= P[good]/2 = {st["anygood"]/K2/2:.4f}')


if __name__ == '__main__':
    main()
