"""s13_adaptive.py -- the orbit method with ADAPTIVE coordinates: separation is automatic.

For a ladder pair (x,y) with mark (p,p') and rho = rho_{x,y}, the column pair (q,c) = (rho^t(p), rho^t(p')) is
automatically separated by {p,p'} (different arcs when crossed, different cycles when parallel) as soon as
t < min(arc lengths) resp. t < min(|Q_p|,|Q_p'|); and this coordinate is invariant under its own turn (the paths
p -> q and p' -> c are untouched) and under the turns of the other ladder pairs' coordinates.  So the genericity
input (beta) of hyp:orbit disappears: a ladder pair contributes a toggling coordinate as soon as ONE of its
candidates t = 1..T has a short CLEAN column cycle through x (length <= LMAX, meets {1,2} in 0 or 2 rows, no other
family row), and the chosen coordinates of different pairs do not collide with other pairs' candidate columns.
With U the set of such pairs,  P[Flip_I = 0 | orbit] <= 2^{-|U|}.

This script measures, on true-marked instances: per pair, the number of candidates, whether a good candidate exists
(t0), |U|; the bound E[2^{-|U|}] against P[Flip=0] (which must hold -- a consistency check of the theory); and it
verifies directly, for pairs in U, that turning Z_k toggles the flippability of (x_k,y_k), leaves the other pairs'
bits and candidate columns unchanged, keeps the type of rows 1,2 and the frame, and that the k-th candidate data
(t0, C_k) are unchanged after the turn (orbit invariance).
usage: ./jm n nsamp thin seed | python3 s13_adaptive.py n [--inst I] [--seed s] [--lmax L] [--T T] [--fam F]
  (family = first min(K0-1, F) ladder pairs; default F = n//8)
"""
import sys, random
from collections import defaultdict, Counter
from verify_identity import sigma_of
from pf_cl import cycles_of
from pf_master import pi_P
from s8_real_bias import read_squares


def rho_of(L, pos, x, y):
    return [pos[y][L[x][q]] for q in range(len(L))]


def flippable(rho, p, pp):
    q = rho[p]
    while q != p:
        if q == pp: return False
        q = rho[q]
    return True


def arcs(rho, p, pp):
    """returns (crossed, a, b): crossed with arc lengths a = d(p->p'), b = d(p'->p), or parallel with a=|Q_p|, b=|Q_p'|."""
    d = 1; q = rho[p]
    while q != p and q != pp: q = rho[q]; d += 1
    if q == pp:
        e = 1; q = rho[pp]
        while q != p: q = rho[q]; e += 1
        return True, d, e
    e = 1; q = rho[pp]
    while q != pp: q = rho[q]; e += 1
    return False, d, e


def col_cycle(L, colpos, q, c, x):
    Z = [x]; r = colpos[c][L[x][q]]
    while r != x: Z.append(r); r = colpos[c][L[r][q]]
    return Z


def turn(L, q, c, rows):
    M = [row[:] for row in L]
    for r in rows: M[r][q], M[r][c] = M[r][c], M[r][q]
    return M


def analyse(L, p, pp, lad, T, LMAX):
    """per ladder pair: candidate count, t0, C_k, chosen (q,c), Z.  Returns list of dicts and U."""
    n = len(L)
    pos = [[0] * n for _ in range(n)]
    for r in range(n):
        for q in range(n): pos[r][L[r][q]] = q
    colpos = [[0] * n for _ in range(n)]
    for r in range(n):
        for q in range(n): colpos[q][L[r][q]] = r
    fam = set(r for xy in lad for r in xy)
    info = []
    for (x, y) in lad:
        rho = rho_of(L, pos, x, y)
        crossed, a, b = arcs(rho, p, pp)
        ncand = min(T, min(a, b) - 1)
        d = {'x': x, 'y': y, 'bit': flippable(rho, p, pp), 'ncand': ncand, 't0': None, 'C': set(), 'Z': None, 'qc': None}
        q, c = p, pp
        for t in range(1, ncand + 1):
            q, c = rho[q], rho[c]
            d['C'].add(q); d['C'].add(c)
            Z = col_cycle(L, colpos, q, c, x)
            if len(Z) > LMAX: continue
            m12 = (0 in Z) + (1 in Z)
            if m12 == 1: continue
            if any(r in fam for r in Z[1:]): continue
            d['t0'] = t; d['Z'] = Z; d['qc'] = (q, c); break
        info.append(d)
    U = []
    for k, d in enumerate(info):
        if d['t0'] is None: continue
        if any(d['qc'][0] in info[l]['C'] or d['qc'][1] in info[l]['C'] for l in range(len(info)) if l != k): continue
        U.append(k)
    return info, U


def main():
    n = int(sys.argv[1]); args = sys.argv[2:]
    def opt(name, default):
        return int(args[args.index(name) + 1]) if name in args else default
    ninst = opt("--inst", 3); seed = opt("--seed", 23); LMAX = opt("--lmax", 8); T = opt("--T", 8); F = opt("--fam", n // 8)
    rng = random.Random(seed)
    nsq = 0; N = 0; st = Counter(); sumU = 0; sum2U = 0.0; flip0 = 0; byK = defaultdict(lambda: [0, 0, 0.0])
    Uhist = Counter(); ncand_hist = Counter(); t0_hist = Counter(); checks = Counter()
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
            if K0 < 2: continue
            m = min(K0 - 1, F)
            lad = []; x, y = 0, 1
            for _ in range(m): x, y = inv[x], inv[y]; lad.append((x, y))
            info, U = analyse(L, p, pp, lad, T, LMAX)
            N += 1; sumU += len(U); sum2U += 2.0 ** (-len(U)); Uhist[len(U)] += 1
            anyflip = any(d['bit'] for d in info)
            flip0 += (not anyflip)
            byK[min(K0, 12)][0] += 1; byK[min(K0, 12)][1] += (not anyflip); byK[min(K0, 12)][2] += 2.0 ** (-len(U))
            for d in info:
                ncand_hist[d['ncand']] += 1; st['pairs'] += 1; st['pairs with t0'] += d['t0'] is not None
                if d['t0'] is not None: t0_hist[d['t0']] += 1
            st['pairs in U'] += len(U)
            # direct verification for one pair in U: turn Z_k, check toggling and invariances
            if U:
                k = rng.choice(U); d = info[k]; q, c = d['qc']
                M = turn(L, q, c, d['Z'])
                assert all(sorted(row) == list(range(n)) for row in M)
                assert sorted(sig) == sorted(sigma_of(M)) and sorted(len(cy) for cy in cycles_of(sig)) == sorted(len(cy) for cy in cycles_of(sigma_of(M)))
                assert pi_P(M, p, pp) == pi
                info2, U2 = analyse(M, p, pp, lad, T, LMAX)
                checks['bit k toggled'] += (info2[k]['bit'] != d['bit'])
                checks['other bits unchanged'] += all(info2[l]['bit'] == info[l]['bit'] for l in range(len(info)) if l != k)
                checks['candidate data unchanged (all pairs)'] += all(info2[l]['t0'] == info[l]['t0'] and info2[l]['C'] == info[l]['C'] and info2[l]['qc'] == info[l]['qc'] and info2[l]['ncand'] == info[l]['ncand'] for l in range(len(info)))
                checks['U unchanged'] += (U2 == U)
                checks['turn checked'] += 1
    print(f"n={n}: squares={nsq} instances={N}  T={T} LMAX={LMAX} family=first min(K0-1,{F}) pairs")
    print(f"  pairs={st['pairs']}  with a good candidate: {st['pairs with t0']/max(1,st['pairs']):.4f}  in U: {st['pairs in U']/max(1,st['pairs']):.4f}")
    print(f"  E|U| = {sumU/N:.3f}   P[Flip=0] = {flip0/N:.4f}   E[2^-|U|] = {sum2U/N:.4f}   (bound must hold: P[Flip=0] <= E[2^-|U|])")
    print("  |U| histogram:", dict(sorted(Uhist.items())))
    print("  candidates per pair (min(T, min arc)-1):", dict(sorted(ncand_hist.items())))
    print("  t0 histogram:", dict(sorted(t0_hist.items())))
    print("  by K0 (K0: instances, P[Flip=0], E[2^-|U|]):", {k: (v[0], round(v[1] / v[0], 3), round(v[2] / v[0], 3)) for k, v in sorted(byK.items())})
    print("  direct checks:", dict(checks))


if __name__ == "__main__":
    main()
