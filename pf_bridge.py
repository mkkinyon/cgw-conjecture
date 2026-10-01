"""Exact bridge between the joining instance space J (on S(mu)) and the
splitting instance space X (on S(lambda)).

Quantities, all by complete enumeration:
  J_A, J_B : joins (L, {j,j'}) with the spanning pair A / B in L
  X_A, X_B : marked instances (L', {j,j'}) with the marked pair A / B in L'
  counting identity: qbar = 2|X|/J - 1  (CGW multigraph balance)

Maps tested:
  unflip: x in X_A  ->  (L, {j,j'}) = flip at the marked A-pair (splits);
          record the status of {j,j'} in L, and whether join(unflip(x)) == x.
  join:   (L,{j,j'}) B-join -> switch (both variants) then flip -> x';
          record marked status of x', multiplicity of hits.
          A-join -> flip -> x'; same records.
"""
from fractions import Fraction
from collections import Counter
import sys

from verify_identity import (latin_squares_first_row_id, sigma_of,
                             cycles_of_perm, row_cycle_columns,
                             delta_same_cycle, omega_map, do_switch, do_flip)
from pf_orbit import typ_of, instances


def canon(L):
    return tuple(map(tuple, L))


def renorm(L):
    """Relabel symbols so row 0 becomes the identity."""
    n = len(L)
    rho = [0] * n
    for c in range(n):
        rho[L[0][c]] = c
    return [[rho[L[i][c]] for c in range(n)] for i in range(n)]


def wdist(L, j, jp):
    """omega-distance from j to jp (steps of omega), or None."""
    w = omega_map(L)
    c = j
    for k in range(1, len(L) + 1):
        c = w[c]
        if c == jp:
            return k
    return None


def run(n, mu, alpha, beta):
    mu = tuple(sorted(mu))
    lam = list(mu)
    lam.remove(alpha)
    lam.remove(beta)
    lam.append(alpha + beta)
    lam = tuple(sorted(lam))

    # ---------- X side ----------
    X_A = X_B = 0
    x_index = {}          # (canon(L'), frozenset({j,jp})) -> marked B?
    xa_unflip_status = Counter()   # status of {j,jp} in unflipped L
    xa_roundtrip = Counter()       # join(unflip(x)) == x ?
    squares_lam = []
    for L in latin_squares_first_row_id(n):
        t = typ_of(L)
        if t == lam:
            squares_lam.append(canon(L))
        elif t == mu:
            pass
    # second pass for mu handled below; store lam squares now
    for Lc in squares_lam:
        L = [list(r) for r in Lc]
        for (j, jp) in instances(L, lam, alpha, beta):
            mb = delta_same_cycle(L, j, jp)
            key = (Lc, frozenset((j, jp)))
            if key in x_index:
                # alpha==beta: same unordered pair marked twice; skip dup
                continue
            x_index[key] = mb
            if mb:
                X_B += 1
            else:
                X_A += 1
                # unflip: flip at the marked A-pair -> splits
                L0 = do_flip(L, j, jp)
                t0 = typ_of(L0)
                assert t0 == mu, (t0, mu)
                st = delta_same_cycle(L0, j, jp)
                xa_unflip_status["B" if st else "A"] += 1
                if not st:
                    # A-join: direct flip back; do we return to x?
                    L1 = do_flip(L0, j, jp)
                    xa_roundtrip[canon(L1) == Lc] += 1

    # ---------- J side ----------
    J_A = J_B = 0
    hitsA = Counter()   # x-key -> #hits by A-join flips
    hitsB = Counter()   # x-key -> #hits by B-join switch(min-variant)+flip
    ja_marked = Counter()   # marked status of x' produced by A-joins
    jb_marked = Counter()   # same for B-joins (deterministic min-symbol switch)
    jb_marked_both = Counter()  # over both switch variants
    for L in latin_squares_first_row_id(n):
        if typ_of(L) != mu:
            continue
        sig = sigma_of(L)
        cyc = cycles_of_perm(sig)
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
                        J_A += 1
                        Lp = renorm(do_flip(L, j, jp))
                        assert typ_of(Lp) == lam
                        key = (canon(Lp), frozenset((j, jp)))
                        hitsA[key] += 1
                        ja_marked[("B" if delta_same_cycle(Lp, j, jp)
                                   else "A", wdist(Lp, j, jp))] += 1
                    else:
                        J_B += 1
                        # min-symbol rule: trade the row cycle containing the
                        # minimal symbol among the two cycles' symbols
                        min_a = min(ca)
                        min_b = min(cb)
                        cols = cols_a if min_a < min_b else cols_b
                        L2 = do_switch(L, cols)
                        assert not delta_same_cycle(L2, j, jp)
                        Lp = renorm(do_flip(L2, j, jp))
                        assert typ_of(Lp) == lam
                        key = (canon(Lp), frozenset((j, jp)))
                        hitsB[key] += 1
                        jb_marked[("B" if delta_same_cycle(Lp, j, jp)
                                   else "A", wdist(Lp, j, jp))] += 1
                        for cols2 in (cols_a, cols_b):
                            L2v = do_switch(L, cols2)
                            Lpv = renorm(do_flip(L2v, j, jp))
                            jb_marked_both[
                                ("B" if delta_same_cycle(Lpv, j, jp)
                                 else "A", wdist(Lpv, j, jp))] += 1

    J = J_A + J_B
    X = X_A + X_B
    print(f"n={n} mu={mu} -> lam={lam} (alpha,beta)=({alpha},{beta})")
    print(f"  J = {J} (A: {J_A}, B: {J_B});  |X| = {X} (A: {X_A}, B: {X_B})")
    if J:
        print(f"  counting qbar = 2|X|/J - 1 = {Fraction(2*X, J)-1} "
              f"= {2*X/J-1:.6f}")
        print(f"  J_A - X_A = {J_A - X_A};   J_B - X_B = {J_B - X_B}")
    print(f"  unflip status of marked-A instances: {dict(xa_unflip_status)}")
    print(f"  A-join roundtrip flip(flip)==id: {dict(xa_roundtrip)}")
    print(f"  marked status of A-join images: {dict(ja_marked)}")
    print(f"  marked status of B-join images (min-rule): {dict(jb_marked)}")
    print(f"  marked status of B-join images (both variants): "
          f"{dict(jb_marked_both)}")
    # multiplicity structure
    ca = Counter(hitsA.values())
    cb = Counter(hitsB.values())
    print(f"  A-join hit multiplicities on X: {dict(ca)}  "
          f"(distinct images: {len(hitsA)})")
    print(f"  B-join hit multiplicities on X: {dict(cb)}  "
          f"(distinct images: {len(hitsB)})")
    both = set(hitsA) & set(hitsB)
    print(f"  images hit by both A- and B-joins: {len(both)}")
    # orphans: X instances never hit
    allhit = set(hitsA) | set(hitsB)
    orphan_A = orphan_B = 0
    for key, mb in x_index.items():
        if key not in allhit:
            if mb:
                orphan_B += 1
            else:
                orphan_A += 1
    print(f"  orphans (never hit): A: {orphan_A}, B: {orphan_B}")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "5"
    if which == "5":
        run(5, (2, 3), 2, 3)
    elif which == "6":
        run(6, (2, 4), 2, 4)
    elif which == "624":
        run(6, (2, 2, 2), 2, 2)
