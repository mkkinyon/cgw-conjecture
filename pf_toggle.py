"""Grounding statistics for the re-randomisation programme (notes sec:kps):

For every L in S(lambda) (n=6, lambda=(6)) and every marked pair
P = {j, sigma^alpha(j)}:
  - U(x) = # intercalates with both rows in completion rows (>=2) meeting
    EXACTLY ONE of the two marked columns ("useful" switches);
  - F(x) = # of those whose switch FLIPS the marked-pair status
    (exact: perform the switch, recompute status; the switch avoids rows
    0,1 so it stays in S(lambda) and preserves the marking);
  - the joint law of (status, U, F): the EQ imbalance should concentrate
    on instances with F = 0.
"""
from collections import Counter
import sys

from verify_identity import (latin_squares_first_row_id, sigma_of,
                             cycles_of_perm, delta_same_cycle)
from pf_orbit import typ_of, instances


def intercalates(L):
    n = len(L)
    out = []
    for r in range(n):
        for rp in range(r + 1, n):
            for c in range(n):
                for cp in range(c + 1, n):
                    if L[r][c] == L[rp][cp] and L[r][cp] == L[rp][c]:
                        out.append((r, rp, c, cp))
    return out


def run(n, lam, alpha, beta):
    lam = tuple(sorted(lam))
    stat = Counter()   # (status, U, F) -> count
    tot = 0
    for L in latin_squares_first_row_id(n):
        if typ_of(L) != lam:
            continue
        ics = [ic for ic in intercalates(L) if ic[0] >= 2]  # rows >= 2
        seen = set()
        for (j, jp) in instances(L, lam, alpha, beta):
            key = frozenset((j, jp))
            if key in seen:
                continue
            seen.add(key)
            tot += 1
            st = delta_same_cycle(L, j, jp)
            useful = [ic for ic in ics
                      if (ic[2] in (j, jp)) != (ic[3] in (j, jp))]
            F = 0
            for (r, rp, c, cp) in useful:
                L2 = [list(row) for row in L]
                L2[r][c], L2[r][cp] = L2[r][cp], L2[r][c]
                L2[rp][c], L2[rp][cp] = L2[rp][cp], L2[rp][c]
                st2 = delta_same_cycle(L2, j, jp)
                if st2 != st:
                    F += 1
            stat[(st, len(useful), F)] += 1
    print(f"n={n} lam={lam} (alpha,beta)=({alpha},{beta}) instances={tot}")
    # marginals
    byU = Counter()
    byF = Counter()
    bal_by_F = Counter()
    for (st, U, F), c in stat.items():
        byU[U] += c
        byF[F] += c
        bal_by_F[(F, st)] += c
    print("  #useful intercalates U:", dict(sorted(byU.items())))
    print("  #flipping F:", dict(sorted(byF.items())))
    print("  status balance by F (F: B-count/A-count, B-A imbalance):")
    Fs = sorted(set(f for (f, s) in bal_by_F))
    totB = tot_imb = 0
    for f in Fs:
        b = bal_by_F.get((f, True), 0)
        a = bal_by_F.get((f, False), 0)
        print(f"    F={f}: B={b} A={a}  B-A={b-a}")
        tot_imb += b - a
        totB += b
    print(f"  total B-A imbalance: {tot_imb}  (= X_B - X_A)")
    f0b = bal_by_F.get((0, True), 0)
    f0a = bal_by_F.get((0, False), 0)
    print(f"  imbalance carried by F=0 stratum: {f0b - f0a}")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "6"
    if which == "6":
        run(6, (6,), 2, 4)
    elif which == "633":
        run(6, (6,), 3, 3)
    elif which == "5":
        run(5, (5,), 2, 3)
