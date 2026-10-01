"""Exact verification of the B-pair identity extracted from the CGW switching
multigraph:

  |S(lambda)| / |S(mu)|  =  (gamma(lambda)/gamma(mu)) * (1 + qbar)/2,
  equivalently  qbar = 2 C(lambda)/C(mu) - 1,

where mu is obtained from lambda by splitting a part of size alpha+beta into
alpha, beta, and qbar is the mean of

  q(L; pair) = 1_B({w'^{-1}(j), w'^{-1}(j')} in L') + 1_B({w'(j), w'(j')} in L')

over uniform L in S(mu, empty), uniform (cycle pair, column pair (j,j')),
where L' = flip(switch-if-B(L)) as in CGW Section 3.

We enumerate ALL Latin squares of order n with first row = identity (symbol
relabelling makes this WLOG; the min-symbol tie-break in the switch is averaged
exactly with weights alpha/(alpha+beta), beta/(alpha+beta)).

Everything exact (integer counts, Fractions at the end).
"""
from fractions import Fraction
from itertools import permutations
import sys

sys.setrecursionlimit(10000)


def latin_squares_first_row_id(n):
    """Yield Latin squares (tuple of rows, each a tuple) with first row id."""
    rows = [tuple(range(n))]
    colmask = [1 << j for j in range(n)]  # symbols used per column
    full = (1 << n) - 1

    def rec(r):
        if r == n:
            yield tuple(rows)
            return
        row = [0] * n

        def fill(c, used):
            if c == n:
                rows.append(tuple(row))
                yield from rec(r + 1)
                rows.pop()
                return
            avail = full & ~used & ~colmask[c]
            while avail:
                b = avail & (-avail)
                avail ^= b
                s = b.bit_length() - 1
                row[c] = s
                colmask[c] |= b
                yield from fill(c + 1, used | b)
                colmask[c] ^= b

        yield from fill(0, 0)

    yield from rec(1)


def sigma_of(L):
    """Row permutation sigma_{1,2}: sigma(L[0][j]) = L[1][j]."""
    n = len(L)
    sig = [0] * n
    for j in range(n):
        sig[L[0][j]] = L[1][j]
    return sig


def cycles_of_perm(p):
    n = len(p)
    seen = [False] * n
    cyc = []
    for i in range(n):
        if not seen[i]:
            c = []
            x = i
            while not seen[x]:
                seen[x] = True
                c.append(x)
                x = p[x]
            cyc.append(c)
    return cyc


def row_cycle_columns(L, cyc_symbols):
    """Columns of the row cycle (rows 1,2) corresponding to a sigma-cycle given
    by its symbols: columns j with L[0][j] in the symbol set."""
    s = set(cyc_symbols)
    return [j for j in range(len(L)) if L[0][j] in s]


def delta_same_cycle(L, j, jp, r1=0, r2=1):
    """Is row r1 in the same column-cycle as row r2 for column pair (j, jp)?
    delta(i) = i' where L[i'][jp] = L[i][j]."""
    n = len(L)
    pos_jp = [0] * n
    for i in range(n):
        pos_jp[L[i][jp]] = i
    i = r1
    while True:
        i = pos_jp[L[i][j]]
        if i == r1:
            return False
        if i == r2:
            return True


def col_cycle_rows(L, j, jp, r):
    """Rows of the column cycle of pair (j,jp) containing row r."""
    n = len(L)
    pos_jp = [0] * n
    pos_j = [0] * n
    for i in range(n):
        pos_jp[L[i][jp]] = i
        pos_j[L[i][j]] = i
    out = [r]
    i = pos_jp[L[r][j]]
    while i != r:
        out.append(i)
        i = pos_jp[L[i][j]]
    return out


def do_switch(L, cols):
    """Trade the row cycle occupying `cols`: swap rows 0,1 entries there."""
    L2 = [list(r) for r in L]
    for c in cols:
        L2[0][c], L2[1][c] = L2[1][c], L2[0][c]
    return L2


def do_flip(L, j, jp):
    """Flip wrt A-pair {j,jp}: trade the column cycle through row index 1."""
    rows = col_cycle_rows(L, j, jp, 1)
    L2 = [list(r) for r in L]
    for i in rows:
        L2[i][j], L2[i][jp] = L2[i][jp], L2[i][j]
    return L2


def omega_map(L):
    """w(c) with row1[w(c)] = row2[c]  (CGW: 1 o w(j) = 2 o j)."""
    n = len(L)
    inv0 = [0] * n
    for c in range(n):
        inv0[L[0][c]] = c
    return [inv0[L[1][c]] for c in range(n)]


def q_events(L, j, jp, alpha, beta):
    """Given L with {j,jp} an A-pair joining an alpha-cycle (through j) and
    beta-cycle (through jp), flip and evaluate the two B-events."""
    Lp = do_flip(L, j, jp)
    w = omega_map(Lp)
    n = len(Lp)
    winv = [0] * n
    for c in range(n):
        winv[w[c]] = c
    e2 = delta_same_cycle(Lp, w[j], w[jp])       # {w'(j), w'(j')}
    e1 = delta_same_cycle(Lp, winv[j], winv[jp])  # {w'^{-1}(j), w'^{-1}(j')}
    # sanity: merged cycle length
    sig = sigma_of(Lp)
    cyc = cycles_of_perm(sig)
    lens = sorted(len(c) for c in cyc)
    return e1, e2, lens


def run(n, mu, alpha, beta, lam_expected, report_each=False):
    """Compute qbar for the join alpha+beta from mu; also P(B1),P(B2)."""
    mu = tuple(sorted(mu))
    # accumulators: [count, sum q, sum e1, sum e2] for A-instances;
    # for B-instances two variants (trade alpha-cycle / beta-cycle)
    A_cnt = A_q = A_e1 = A_e2 = 0
    B_cnt = 0
    BX_q = BX_e1 = BX_e2 = 0
    BY_q = BY_e1 = BY_e2 = 0
    n_sq = 0
    n_mu = 0
    for L in latin_squares_first_row_id(n):
        n_sq += 1
        sig = sigma_of(L)
        cyc = cycles_of_perm(sig)
        typ = tuple(sorted(len(c) for c in cyc))
        if typ != mu:
            continue
        n_mu += 1
        acyc = [c for c in cyc if len(c) == alpha]
        bcyc = [c for c in cyc if len(c) == beta]
        if alpha == beta:
            pairs = [(acyc[i], acyc[k]) for i in range(len(acyc))
                     for k in range(i + 1, len(acyc))]
        else:
            pairs = [(ca, cb) for ca in acyc for cb in bcyc]
        for (ca, cb) in pairs:
            cols_a = row_cycle_columns(L, ca)
            cols_b = row_cycle_columns(L, cb)
            for j in cols_a:
                for jp in cols_b:
                    isB = delta_same_cycle(L, j, jp)
                    if not isB:
                        e1, e2, lens = q_events(L, j, jp, alpha, beta)
                        exp_lens = sorted(list(mu))
                        exp_lens.remove(alpha)
                        exp_lens.remove(beta)
                        exp_lens.append(alpha + beta)
                        assert lens == sorted(exp_lens), (lens, exp_lens)
                        A_cnt += 1
                        A_q += e1 + e2
                        A_e1 += e1
                        A_e2 += e2
                    else:
                        B_cnt += 1
                        for variant, cols in (("X", cols_a), ("Y", cols_b)):
                            L2 = do_switch(L, cols)
                            # after switch, {j,jp} must be an A-pair
                            assert not delta_same_cycle(L2, j, jp)
                            e1, e2, lens = q_events(L2, j, jp, alpha, beta)
                            if variant == "X":
                                BX_q += e1 + e2
                                BX_e1 += e1
                                BX_e2 += e2
                            else:
                                BY_q += e1 + e2
                                BY_e1 += e1
                                BY_e2 += e2
    tot = A_cnt + B_cnt
    wA = Fraction(alpha, alpha + beta)
    wB = Fraction(beta, alpha + beta)
    qbar = (Fraction(A_q) + wA * BX_q + wB * BY_q) / tot
    p1 = (Fraction(A_e1) + wA * BX_e1 + wB * BY_e1) / tot
    p2 = (Fraction(A_e2) + wA * BX_e2 + wB * BY_e2) / tot
    print(f"n={n} mu={mu} join {alpha}+{beta}:")
    print(f"  squares(first row id) = {n_sq}, in S(mu): {n_mu}, instances: {tot} (A: {A_cnt}, B: {B_cnt})")
    print(f"  qbar = {qbar} = {float(qbar):.6f}   expected {lam_expected} = {float(lam_expected):.6f}"
          f"   {'MATCH' if qbar == lam_expected else 'MISMATCH'}")
    print(f"  P(B1) = {p1} = {float(p1):.6f},  P(B2) = {p2} = {float(p2):.6f}")
    return qbar


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "5"
    if which == "5":
        run(5, (2, 3), 2, 3, Fraction(2))
    elif which == "6a":
        run(6, (2, 4), 2, 4, Fraction(10, 11))
    elif which == "6b":
        run(6, (3, 3), 3, 3, Fraction(3, 4))
    elif which == "6c":
        run(6, (2, 2, 2), 2, 2, Fraction(4, 7))
