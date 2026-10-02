"""s13_cgwcase3.py -- independent check of CGW08 Lemma 3.12, splitting case 3.

CGW's cross-switch at a B-pair (j,j') with omega^alpha(j)=j' (Def 3.9, first form, min symbol in
{e1,e4}): rows 1,2 on the omega-arc from j to j' and columns j,j' on the delta-arc from row 2 to
row 1 (the arc that balances the corner cells).  Rows 1,2 part (from the proof of Lemma 3.10): row 1 at omega^k j (1<=k<alpha) gets the old
row-2 symbol there, row 1 at j' gets e3 = L(2,j); row 2 at omega^k j gets the old row-1 symbol,
row 2 at j gets e2 = L(1,j').  Columns j,j' part: rows delta^k(2), 1<=k<b (delta = delta_{j,j'},
delta^b(2) = 1) swap their entries at j and j'.  (Second form: the same with the roles of
(j, row 1) and (j', row 2) exchanged.)  We verify that the result is a Latin square of the same
type (Lemma 3.10), then apply the splitting procedure's case 3: at a pair (j,j') that is B with
B-neighbour (omega j, omega j'), cross-switch at the neighbour and backflip at (j,j'), and report
the resulting type.  CGW claim it lies in S(mu) (the (alpha+beta)-cycle split into alpha, beta).
usage: ./jm n nsamp thin seed | python3 s13_cgwcase3.py n
"""
import sys
from collections import Counter
from s8_real_bias import read_squares


def omega_of(L):
    n = len(L); col1 = {L[0][c]: c for c in range(n)}
    return [col1[L[1][c]] for c in range(n)]          # 1 o omega(j) = 2 o j


def delta_of(L, j, jp):
    n = len(L); row = {L[r][jp]: r for r in range(n)}
    return [row[L[r][j]] for r in range(n)]           # delta(i) o j' = i o j


def cycles(perm):
    n = len(perm); seen = [False] * n; out = []
    for s in range(n):
        if seen[s]: continue
        c = []; r = s
        while not seen[r]: seen[r] = True; c.append(r); r = perm[r]
        out.append(c)
    return out


def ctype(L): return tuple(sorted(len(c) for c in cycles(omega_of(L))))


def is_latin(L):
    n = len(L)
    return all(sorted(r) == list(range(n)) for r in L) and all(sorted(L[r][c] for r in range(n)) == list(range(n)) for c in range(n))


def col_cycle_through(L, j, jp, r0):
    d = delta_of(L, j, jp); Z = [r0]; r = d[r0]
    while r != r0: Z.append(r); r = d[r]
    return Z


def turn(L, j, jp, rows):
    M = [list(r) for r in L]
    for r in rows: M[r][j], M[r][jp] = M[r][jp], M[r][j]
    return M


def is_A(L, j, jp): return 1 not in col_cycle_through(L, j, jp, 0)


def cross_switch(L, j, jp):
    """CGW Def 3.9 at the B-pair (j,jp) with omega^alpha(j) = jp (rows 1,2 = indices 0,1)."""
    om = omega_of(L); n = len(L)
    e1, e2, e3, e4 = L[0][j], L[0][jp], L[1][j], L[1][jp]
    if min(e1, e2, e3, e4) in (e2, e3):               # second form: exchange roles
        return cross_switch_form(L, jp, j, 1, 0)
    return cross_switch_form(L, j, jp, 0, 1)


def cross_switch_form(L, j, jp, r1, r2):
    """first form with distinguished rows (r1,r2) and omega-arc from j to jp in the direction
    'r1 o omega(c) = r2 o c'."""
    n = len(L); M = [list(r) for r in L]
    col = {L[r1][c]: c for c in range(n)}; om = [col[L[r2][c]] for c in range(n)]
    arc = []; c = om[j]
    while c != jp: arc.append(c); c = om[c]
    for c in arc: M[r1][c], M[r2][c] = L[r2][c], L[r1][c]
    M[r1][jp] = L[r2][j]; M[r2][j] = L[r1][jp]
    # delta-arc from r2 to r1 in columns j,jp (interior rows swap); delta(i) o jp = i o j
    row = {L[r][jp]: r for r in range(n)}; d = [row[L[r][j]] for r in range(n)]
    r = d[r2]
    while r != r1: M[r][j], M[r][jp] = L[r][jp], L[r][j]; r = d[r]
    return M


def main():
    n = int(sys.argv[1]); stats = Counter(); examples = 0
    for L in read_squares(n, sys.stdin.buffer):
        om = omega_of(L); lam = ctype(L)
        for cyc in cycles(om):
            m = len(cyc)
            if m < 6: continue
            for alpha in range(3, m - 2):
                beta = m - alpha
                for i, j in enumerate(cyc):
                    jp = cyc[(i + alpha) % m]
                    if is_A(L, j, jp): continue
                    J, Jp = om[j], om[jp]
                    if is_A(L, J, Jp): continue
                    # case 3: cross-switch at the neighbour, then backflip at (j,jp)
                    L1 = cross_switch(L, J, Jp)
                    ok1 = is_latin(L1) and ctype(L1) == lam
                    stats['cross-switch Latin & type preserved' if ok1 else 'cross-switch FAILS'] += 1
                    if not ok1: continue
                    om1 = omega_of(L1)
                    dist = 1; c = om1[j]
                    while c != jp: c = om1[c]; dist += 1
                    stats[f'distance j->jp after cross-switch at neighbour: {"2" if dist == 2 else ("alpha" if dist == alpha else "other")}'] += 1
                    if not is_A(L1, j, jp):
                        stats['(j,jp) still B after cross-switch'] += 1; continue
                    L2 = turn(L1, j, jp, col_cycle_through(L1, j, jp, 1))     # backflip: cycle through row 2
                    mu = sorted(list(lam)); mu.remove(m); mu = tuple(sorted(mu + [alpha, beta]))
                    t2 = ctype(L2)
                    stats['backflip lands in S(mu)' if t2 == mu else 'backflip lands in WRONG type'] += 1
                    if t2 != mu and examples < 3:
                        examples += 1
                        print(f"example: n={n} lambda={lam} (alpha,beta)=({alpha},{beta}) -> type after case 3 = {t2}, expected {mu}")
    for k, v in sorted(stats.items()): print(f"{v:7d}  {k}")


if __name__ == "__main__":
    main()
