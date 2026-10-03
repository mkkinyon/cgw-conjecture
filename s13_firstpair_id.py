"""s13_firstpair_id.py -- the first-pair identity  P[B] - P[A] = E[(1_B - 1_A)(1_{F^c} - 1_F); G^c]
(F = (x1,y1) flippable, G = some candidate t < mu is good; see prop:firstpair and section (k)), and its pieces
   P[B & F] = P[A & F]                    (trade T_1)
   P[B & F & G] = P[B & F^c & G],  P[A & F & G] = P[A & F^c & G]     (the turn at the least good t)
on squares read from stdin (n*n bytes each; rows 0,1 = (id, sigma)): exact if stdin is the enumeration of all
completions, sampled if it is ./jmfix output.  All marks of each class are used unless --marks M is given.
usage: python3 s13_firstpair_id.py n [--marks M] [--seed s] [--classes m:alpha,...] < squares.bin
"""
import sys, random
from collections import defaultdict, Counter
from s8_real_bias import read_squares
from pf_cl import cycles_of


def analyse(L, n, pos, p, pp):
    colpp = [0] * n
    for r in range(n): colpp[L[r][pp]] = r
    pi = [colpp[L[r][p]] for r in range(n)]
    inv = [0] * n
    for r in range(n): inv[pi[r]] = r
    r = 0; together = False
    while True:
        r = pi[r]
        if r == 1: together = True
        if r == 0: break
    x1, y1 = inv[0], inv[1]
    assert x1 >= 2 and y1 >= 2
    rho = [pos[y1][L[x1][q]] for q in range(n)]
    d = 1; q = rho[p]
    while q != p and q != pp: q = rho[q]; d += 1
    crossed = (q == pp)
    if crossed:
        e = 1; q = rho[pp]
        while q != p: q = rho[q]; e += 1
    else:
        e = 1; q = rho[pp]
        while q != pp: q = rho[q]; e += 1
    mu = min(d, e)
    F = not crossed
    q, c = p, pp; good = False
    for t in range(1, mu):
        q, c = rho[q], rho[c]
        colc = [0] * n
        for r in range(n): colc[L[r][c]] = r
        Z = [x1]; r = colc[L[x1][q]]
        while r != x1: Z.append(r); r = colc[L[r][q]]
        if y1 in Z: continue
        Zy = [y1]; r = colc[L[y1][q]]
        while r != y1: Zy.append(r); r = colc[L[r][q]]
        if ((0 in Z) + (1 in Z)) != 1 or ((0 in Zy) + (1 in Zy)) != 1:
            good = True; break
    return together, F, good, mu


def main():
    n = int(sys.argv[1]); args = sys.argv[2:]
    def opt(name, default, conv=int):
        return conv(args[args.index(name) + 1]) if name in args else default
    M = opt('--marks', 0); seed = opt('--seed', 3)
    classes = opt('--classes', None, str)
    classes = None if classes is None else set(tuple(int(t) for t in c.split(':')) for c in classes.split(','))
    rng = random.Random(seed)
    sig = None; marks_by_class = None; S = defaultdict(Counter); nsq = 0
    for L in read_squares(n, sys.stdin.buffer):
        nsq += 1
        if sig is None:
            assert L[0] == list(range(n)); sig = L[1][:]
            marks_by_class = defaultdict(list)
            for cyc in cycles_of(sig):
                m = len(cyc)
                if m < 4: continue
                for ji, j in enumerate(cyc):
                    for alpha in range(2, m - 1):
                        if classes is None or (m, alpha) in classes:
                            marks_by_class[(m, alpha)].append((j, cyc[(ji + alpha) % m]))
        pos = [[0] * n for _ in range(n)]
        for r in range(n):
            for q in range(n): pos[r][L[r][q]] = q
        for cls, marks in marks_by_class.items():
            use = marks if M == 0 else rng.sample(marks, min(M, len(marks)))
            st = S[cls]
            for (p, pp) in use:
                B, F, G, mu = analyse(L, n, pos, p, pp)
                st['N'] += 1
                st[(B, F, G)] += 1
                st[('mu', min(mu, 12))] += 1
                if not G: st[('Gc_mu', min(mu, 12))] += 1
    print(f'n={n}: squares={nsq} type={sorted((len(c) for c in cycles_of(sig)), reverse=True)} marks={"all" if M == 0 else M}')
    for cls in sorted(S):
        st = S[cls]; N = st['N']
        P = lambda *keys: sum(st[k] for k in keys) / N
        B = lambda F, G: P((True, F, G)); A = lambda F, G: P((False, F, G))
        lhs = (B(True, True) + B(True, False) + B(False, True) + B(False, False)) - (A(True, True) + A(True, False) + A(False, True) + A(False, False))
        rhs = (B(False, False) - B(True, False)) - (A(False, False) - A(True, False))
        PGc = B(True, False) + B(False, False) + A(True, False) + A(False, False)
        print(f' class m={cls[0]} alpha={cls[1]}: N={N}')
        print(f'   P[B]-P[A] = {lhs:+.6f}   RHS E[(1_B-1_A)(1_Fc-1_F);Gc] = {rhs:+.6f}   P[G^c] = {PGc:.6f}   bound |P[B]-1/2| <= P[G^c]/2: {abs(lhs)/2:.6f} <= {PGc/2:.6f}')
        print(f'   T_1:  P[B&F] = {B(True,True)+B(True,False):.6f}  P[A&F] = {A(True,True)+A(True,False):.6f}')
        print(f'   turn: P[B&F&G] = {B(True,True):.6f}  P[B&Fc&G] = {B(False,True):.6f}   P[A&F&G] = {A(True,True):.6f}  P[A&Fc&G] = {A(False,True):.6f}')
        print('   mu histogram:', {k: round(st[("mu", k)] / N, 4) for k in range(2, 13)})
        print('   P[G^c | mu=k] vs 0.59^(k-1):', {k: (round(st[('Gc_mu', k)] / st[('mu', k)], 3), round(0.59 ** (k - 1), 3)) for k in range(2, 13) if st[('mu', k)] >= 100})


if __name__ == '__main__':
    main()
