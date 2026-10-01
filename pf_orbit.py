"""Splitting-side (post-flip) instance space and the CGW cross-switch.

Instance space X = X(n, lambda, alpha, beta):
  (L', c, j) with L' in S(lambda) (first row id), c a merged (alpha+beta)-cycle
  of sigma(L'), j a column of c, marked pair {j, j'} with j' = omega^alpha(j)
  (so the pair is at omega-distances (alpha, beta)).

Events (as in the notes, sec:identity):
  B1: {w^{-1}(j), w^{-1}(j')} is a B-pair of L'
  B2: {w(j),     w(j')}     is a B-pair of L'
  MB: the marked pair {j,j'} itself is a B-pair of L'.

Cross-switch (CGW08 Definition 3.9; the paper is in the tarball as
CGW_2008.pdf, Def 3.9 on pp. 292-293): defined at a B-pair {j,j'} lying in
one row cycle, j' = omega^alpha(j), alpha,beta >= 2.  Trades one
row-1/row-2 arc plus one column-j/column-j' delta-path.  Lemma 3.10
(verbatim family): L' stays in S(lambda) and every pair
{w^k(j), w^k'(j')} with INDEPENDENT 1 <= k < alpha, 1 <= k' < beta
toggles A/B status.  Both shifted pairs are in this family:
P_{+1} = (k,k')=(1,1) and P_{-1} = (alpha-1, beta-1) -- so their
deterministic toggling (observed below at n=5) is a theorem.
Equal-shift pairs {w^k(j), w^k(j')} with k >= min(alpha,beta) are NOT
covered, matching the partial toggle rates observed.  Corollary 3.11
(sic; involution).

Everything exact; complete enumeration at n=5,6.
"""
from fractions import Fraction
from collections import Counter
import sys

from verify_identity import (latin_squares_first_row_id, sigma_of,
                             cycles_of_perm, row_cycle_columns,
                             delta_same_cycle, omega_map)


def col_delta(L, j, jp):
    """delta_{j,jp}: row map i -> i' with L[i'][jp] = L[i][j]."""
    n = len(L)
    pos_jp = [0] * n
    for i in range(n):
        pos_jp[L[i][jp]] = i
    return [pos_jp[L[i][j]] for i in range(n)]


def cross_switch(L, j, jp, alpha, beta):
    """CGW Def 3.9 at the B-pair {j,jp}, jp = omega^alpha(j).
    Returns L2 (list of lists) or raises AssertionError on precondition."""
    n = len(L)
    w = omega_map(L)
    # arcs: from j, omega^alpha reaches jp; from jp, omega^beta reaches j
    arc_j = []  # omega^k(j), 1<=k<alpha  (interior of j->jp arc)
    c = j
    for _ in range(1, alpha):
        c = w[c]
        arc_j.append(c)
    assert w[c] == jp if alpha >= 1 else True
    arc_jp = []  # omega^k(jp), 1<=k<beta (interior of jp->j arc)
    c = jp
    for _ in range(1, beta):
        c = w[c]
        arc_jp.append(c)
    assert w[c] == j
    e1, e2, e3, e4 = L[0][j], L[0][jp], L[1][j], L[1][jp]
    assert len({e1, e2, e3, e4}) == 4
    d = col_delta(L, j, jp)
    # delta-path row 1 -> row 2 (0-indexed 0 -> 1) and row 2 -> row 1
    a = 1
    i = d[0]
    path_a = []  # delta^k(1), 1<=k<a (interior rows of 1->2 path)
    while i != 1:
        path_a.append(i)
        i = d[i]
        a += 1
    b = 1
    i = d[1]
    path_b = []  # delta^k(2), 1<=k<b
    while i != 0:
        path_b.append(i)
        i = d[i]
        b += 1
    L2 = [list(r) for r in L]
    if min(e1, e2, e3, e4) in (e1, e4):
        # case 1: trade (2,j)<->(1,jp) endpoints, arc_j rows-swap,
        # column-swap on the interior of the delta path row2 -> row1
        L2[1][j], L2[0][jp] = e2, e3
        for cc in arc_j:
            L2[0][cc], L2[1][cc] = L2[1][cc], L2[0][cc]
        for rr in path_b:
            L2[rr][j], L2[rr][jp] = L2[rr][jp], L2[rr][j]
    else:
        # case 2: trade (1,j)<->(2,jp) endpoints, arc_jp rows-swap,
        # column-swap on the interior of the delta path row1 -> row2
        L2[0][j], L2[1][jp] = e4, e1
        for cc in arc_jp:
            L2[0][cc], L2[1][cc] = L2[1][cc], L2[0][cc]
        for rr in path_a:
            L2[rr][j], L2[rr][jp] = L2[rr][jp], L2[rr][j]
    return L2


def is_latin(L):
    n = len(L)
    for r in L:
        if len(set(r)) != n:
            return False
    for c in range(n):
        if len(set(L[i][c] for i in range(n))) != n:
            return False
    return True


def typ_of(L):
    return tuple(sorted(len(c) for c in cycles_of_perm(sigma_of(L))))


def instances(L, lam, alpha, beta):
    """Yield (j, jp) marked instances: jp = omega^alpha(j), j in a merged
    (alpha+beta)-cycle.  For alpha != beta each unordered pair once; for
    alpha == beta each unordered pair arises twice (j and jp) -- keep both
    (uniform over ordered marks)."""
    sig = sigma_of(L)
    w = omega_map(L)
    for cyc in cycles_of_perm(sig):
        if len(cyc) != alpha + beta:
            continue
        cols = row_cycle_columns(L, cyc)
        for j in cols:
            jp = j
            for _ in range(alpha):
                jp = w[jp]
            yield j, jp


def shifted_pairs(L, j, jp):
    w = omega_map(L)
    n = len(L)
    winv = [0] * n
    for c in range(n):
        winv[w[c]] = c
    return (winv[j], winv[jp]), (w[j], w[jp])


def run(n, lam, alpha, beta, maxsq=None, do_cross=True):
    lam = tuple(sorted(lam))
    tot = 0
    MB = 0          # marked pair B
    B1 = B2 = 0
    B1_MB = B2_MB = 0     # joint with marked B
    B1_MA = B2_MA = 0     # joint with marked A
    # cross-switch property counters (over marked-B instances)
    cs_fail = 0
    cs_latin = cs_type = cs_invol = cs_markB = 0
    toggle = Counter()    # k -> #toggles of pair P_k among cross-switched
    toggle_tot = 0
    nsq = 0
    for L in latin_squares_first_row_id(n):
        if typ_of(L) != lam:
            continue
        nsq += 1
        if maxsq and nsq > maxsq:
            break
        w = omega_map(L)
        for (j, jp) in instances(L, lam, alpha, beta):
            tot += 1
            mb = delta_same_cycle(L, j, jp)
            p1, p2 = shifted_pairs(L, j, jp)
            b1 = delta_same_cycle(L, *p1)
            b2 = delta_same_cycle(L, *p2)
            MB += mb
            B1 += b1
            B2 += b2
            if mb:
                B1_MB += b1
                B2_MB += b2
            else:
                B1_MA += b1
                B2_MA += b2
            if do_cross and mb:
                toggle_tot += 1
                L2 = cross_switch(L, j, jp, alpha, beta)
                ok_lat = is_latin(L2)
                ok_typ = ok_lat and (typ_of(L2) == lam)
                cs_latin += ok_lat
                cs_type += ok_typ
                if not ok_typ:
                    cs_fail += 1
                    continue
                # involution?
                L3 = cross_switch(L2, j, jp, alpha, beta)
                cs_invol += (tuple(map(tuple, L3)) == tuple(map(tuple, L)))
                # marked pair still B in L2?
                cs_markB += delta_same_cycle(L2, j, jp)
                # toggle profile over translates P_k, k=1..alpha+beta-1
                u, up = j, jp
                for k in range(1, alpha + beta):
                    u, up = w[u], w[up]  # NOTE: w of L (should equal w of L2
                    # on the cycle? not necessarily -- use L's w for labeling)
                    s_old = delta_same_cycle(L, u, up)
                    s_new = delta_same_cycle(L2, u, up)
                    if s_old != s_new:
                        toggle[k] += 1
    print(f"n={n} lam={lam} (alpha,beta)=({alpha},{beta}): "
          f"squares={nsq} instances={tot}")
    if tot == 0:
        return
    print(f"  P(marked B) = {Fraction(MB, tot)} = {MB/tot:.6f}")
    print(f"  P(B1) = {Fraction(B1, tot)} = {B1/tot:.6f}   "
          f"P(B2) = {Fraction(B2, tot)} = {B2/tot:.6f}")
    if MB:
        print(f"  P(B1|MB) = {Fraction(B1_MB, MB)} = {B1_MB/MB:.6f}   "
              f"P(B2|MB) = {Fraction(B2_MB, MB)} = {B2_MB/MB:.6f}")
    MA = tot - MB
    if MA:
        print(f"  P(B1|MA) = {Fraction(B1_MA, MA)} = {B1_MA/MA:.6f}   "
              f"P(B2|MA) = {Fraction(B2_MA, MA)} = {B2_MA/MA:.6f}")
    if do_cross and toggle_tot:
        print(f"  cross-switch on {toggle_tot} marked-B instances: "
              f"latin {cs_latin}/{toggle_tot}, type {cs_type}/{toggle_tot}, "
              f"involution {cs_invol}/{cs_type}, markB-preserved "
              f"{cs_markB}/{cs_type}, hard-fail {cs_fail}")
        print("  toggle profile (k: toggles/opportunities):")
        for k in range(1, alpha + beta):
            print(f"    k={k}: {toggle[k]}/{cs_type}"
                  f"  ({toggle[k]/max(cs_type,1):.4f})")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "5"
    if which == "5":
        run(5, (5,), 2, 3)
    elif which == "6":
        run(6, (6,), 2, 4)
    elif which == "633":
        run(6, (6,), 3, 3)
    elif which == "624":
        run(6, (2, 4), 2, 2)
