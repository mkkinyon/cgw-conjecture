"""End-to-end grounding of the master formula (Theorem thm:master) on
REAL Latin squares at n=6.

For every L in S(lambda) with first row id, every marked pair
P = {j, sigma^alpha(j)}, and a greedy maximal family of pairwise
cell-disjoint intercalates in rows >= 2 (0-indexed; i.e. avoiding rows
1,2 in 1-indexed terms), classified by column incidence with P:

  - j-only  -> right chord (copy 1)
  - j'-only -> left chord  (copy 0)
  - both    -> conjugating coordinate (both lifts, tied)
  - neither -> must not change the bit (verified too)

check for EVERY activation subset that the rank-event prediction of
the master formula (doubling + bridges + bordered rank over the
pi_P(L)-cycle structure) equals the actual A/B status of P in the
switched square.

Usage: python3 pf_master.py 6 [maxfam]
"""
import sys
from collections import Counter

from verify_identity import (latin_squares_first_row_id, sigma_of,
                             delta_same_cycle)
from pf_orbit import typ_of, instances
from pf_toggle import intercalates
from pf_cl import (cycles_of, interlacement_data, in_colspace,
                   multi_cycle_augment, same_cycle)


def pi_P(L, j, jp):
    n = len(L)
    # pi(r) = row r' with L[r'][jp] == L[r][j]
    row_of = [0] * n
    for r in range(n):
        row_of[L[r][jp]] = r
    return [row_of[L[r][j]] for r in range(n)]


def switch_ic(L, ic):
    r, rp, c, cp = ic
    L2 = [list(row) for row in L]
    L2[r][c], L2[r][cp] = L2[r][cp], L2[r][c]
    L2[rp][c], L2[rp][cp] = L2[rp][cp], L2[rp][c]
    return L2


def run(n, maxfam=10, stride=1):
    total = bad = 0
    inst = 0
    fam_sizes = Counter()
    for nsq, L in enumerate(latin_squares_first_row_id(n)):
        if nsq % stride:
            continue
        lam = typ_of(L)
        seen = set()
        for (j, jp) in instances(L, lam, None, None) if False else []:
            pass
        # enumerate marked pairs over all merges (alpha from 2..m-2)
        sig = sigma_of(L)
        # find cycles of sigma
        cyc = cycles_of(sig)
        for c in cyc:
            m = len(c)
            if m < 4:
                continue
            for alpha in range(2, m - 1):
                for j in c:
                    jp = j
                    for _ in range(alpha):
                        jp = sig[jp]
                    key = frozenset((j, jp))
                    if key in seen or j == jp:
                        continue
                    seen.add(key)
                    inst += 1
                    # greedy cell-disjoint intercalate family, rows >= 2
                    ics = [ic for ic in intercalates(L) if ic[0] >= 2]
                    fam = []
                    used = set()
                    for (r, rp, cc, ccp) in ics:
                        cells = {(r, cc), (r, ccp), (rp, cc), (rp, ccp)}
                        if cells & used:
                            continue
                        # classify
                        meets = (cc in (j, jp)) + (ccp in (j, jp))
                        fam.append((r, rp, cc, ccp))
                        used |= cells
                        if len(fam) >= maxfam:
                            break
                    fam_sizes[len(fam)] += 1
                    # build master-formula data from pi0 = pi_P(L)
                    pi0 = pi_P(L, j, jp)
                    pcyc = cycles_of(pi0)
                    # doubled circle
                    lab = {}
                    dcyc = []
                    for cyc0 in pcyc:
                        sc = []
                        for x in cyc0:
                            sc.append((x, 0))
                            sc.append((x, 1))
                        dcyc.append(sc)
                    for sc in dcyc:
                        for q in sc:
                            lab[q] = len(lab)
                    icyc = [[lab[q] for q in sc] for sc in dcyc]
                    circle, bridges, tot = multi_cycle_augment(
                        icyc, len(lab))
                    # per-coordinate lifted chords
                    coord_chords = []
                    for (r, rp, cc, ccp) in fam:
                        mj = (cc == j) + (ccp == j)
                        mjp = (cc == jp) + (ccp == jp)
                        ch = []
                        if mj:   # right side -> copy 1
                            ch.append((lab[(r, 1)], lab[(rp, 1)]))
                        if mjp:  # left side -> copy 0
                            ch.append((lab[(r, 0)], lab[(rp, 0)]))
                        coord_chords.append(ch)  # may be [] (neither)
                    k = len(fam)
                    u2, v2 = lab[(0, 0)], lab[(1, 0)]
                    for Sb in range(1 << k):
                        act = list(bridges)
                        L2 = L
                        for i in range(k):
                            if (Sb >> i) & 1:
                                act.extend(coord_chords[i])
                                L2 = switch_ic(L2, fam[i])
                        rows, x = interlacement_data(circle, act, u2, v2)
                        pred = in_colspace(rows, len(act), x)
                        truth = delta_same_cycle(L2, j, jp)
                        total += 1
                        if pred != truth:
                            bad += 1
                            if bad <= 3:
                                print("FAIL", L, (j, jp), fam, Sb)
    print(f"n={n}: instances={inst} subset-checks={total} "
          f"failures={bad}")
    print("  family sizes:", dict(sorted(fam_sizes.items())))


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    mf = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    st = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    run(n, mf, st)
