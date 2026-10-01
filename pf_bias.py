"""Bias laboratory for the orbit-mixing core (session 6).

By the verified master formula (pf_cl.py), the conditional law of the
marked bit given an orbit is

    P_S(bit) = P_{S uniform}[ x_{S u B} in colspace(A_{S u B}) ]

for a FIXED interlacement graph G on ground set [k] u B (bridges frozen
in).  Equivalently, with e the virtual u-v chord adjoined,

    bias := P(bit) - 1/2 = (1/2) E_S[ null(G[S u B u e]) - null(G[S u B]) ].

This file studies bias exactly (exhaustive over S) for structured and
random chord geometries:

  1. pure ladder (r pairwise-interlacing separating chords): bias = 0
     for r >= 1 (bit = parity) -- verify.
  2. ladder + one anti-parallel separating chord: bias = -1/4 -- verify.
  3. nested separating family: bias = 2^-k - 1/2 -- verify.
  4. random geometric ensemble ("generic position"): u,v on a circle,
     k disjoint chords with uniform random endpoints; distribution of
     |bias| vs k.  Also conditional on both u-v arcs being long/short.
  5. multi-cycle bases with frozen bridges.
"""
import itertools
import random
import sys
from collections import Counter

from pf_cl import (cycles_of, same_cycle, right_mult_transpositions,
                   f2_rank, in_colspace, interlacement_data,
                   multi_cycle_augment, circle_positions)


def bias_exact(circle, chords, frozen, u, v):
    """Exhaustive bias over uniform subsets of `chords`; `frozen` always
    active."""
    k = len(chords)
    tot = 0
    for Sb in range(1 << k):
        act = list(frozen) + [chords[i] for i in range(k) if (Sb >> i) & 1]
        rows, x = interlacement_data(circle, act, u, v)
        if in_colspace(rows, len(act), x):
            tot += 1
    return tot / (1 << k) - 0.5


def bias_mc(circle, chords, frozen, u, v, samples, rng):
    k = len(chords)
    tot = 0
    for _ in range(samples):
        act = list(frozen) + [chords[i] for i in range(k)
                              if rng.random() < 0.5]
        rows, x = interlacement_data(circle, act, u, v)
        if in_colspace(rows, len(act), x):
            tot += 1
    return tot / samples - 0.5


# --------------------------------------------------- structured examples

def structured():
    N = 64
    circle = list(range(N))
    u, v = 0, N // 2
    # 1. ladder: a_i in (u,v) increasing, b_i in (v,u) increasing
    for r in (1, 2, 3, 5, 8):
        chords = [(1 + i, N // 2 + 1 + i) for i in range(r)]
        b = bias_exact(circle, chords, [], u, v)
        print(f"ladder r={r}: bias = {b:+.6f}")
    # 2. ladder + anti-parallel separating chord
    for r in (2, 4, 6):
        chords = [(2 + i, N // 2 + 2 + i) for i in range(r)]
        chords.append((1, N - 1))  # cuts off u, crosses nothing
        b = bias_exact(circle, chords, [], u, v)
        print(f"ladder r={r} + uncrossed cut chord: bias = {b:+.6f}")
    # 3. nested separating family (pairwise non-interlacing)
    for r in (2, 4, 6):
        chords = [(1 + i, N - 1 - i) for i in range(r)]
        b = bias_exact(circle, chords, [], u, v)
        print(f"nested r={r}: bias = {b:+.6f}  (predict {2**-r - 0.5:+.6f})")


# ---------------------------------------------- random geometric ensemble

def random_geometric(kvals=(4, 6, 8, 10, 12), trials=400, N=200, seed=5,
                     antipodal=True):
    rng = random.Random(seed)
    print(f"random geometric ensemble, N={N}, antipodal={antipodal}, "
          f"{trials} instances per k")
    for k in kvals:
        biases = []
        for t in range(trials):
            circle = list(range(N))
            if antipodal:
                u, v = 0, N // 2
            else:
                u, v = 0, rng.randrange(1, N)
            pool = [p for p in range(N) if p not in (u, v)]
            rng.shuffle(pool)
            chords = [(pool[2 * i], pool[2 * i + 1]) for i in range(k)]
            biases.append(bias_exact(circle, chords, [], u, v))
        ab = sorted(abs(b) for b in biases)
        mean = sum(biases) / len(biases)
        meanabs = sum(ab) / len(ab)
        q90 = ab[int(0.9 * len(ab))]
        mx = ab[-1]
        print(f"  k={k:3d}: mean bias {mean:+.5f}  E|bias| {meanabs:.5f}  "
              f"q90 {q90:.5f}  max {mx:.5f}")


def random_geometric_short_arc(kvals=(6, 10), trials=300, N=200, seed=6,
                               arclen=3):
    """u,v close together (short arc): the hard case flagged in the
    handoff.  Chord endpoints uniform: rarely inside the short arc."""
    rng = random.Random(seed)
    print(f"short-arc ensemble: arc(u->v) length {arclen}, N={N}")
    for k in kvals:
        biases = []
        for t in range(trials):
            circle = list(range(N))
            u, v = 0, arclen
            pool = [p for p in range(N) if p not in (u, v)]
            rng.shuffle(pool)
            chords = [(pool[2 * i], pool[2 * i + 1]) for i in range(k)]
            biases.append(bias_exact(circle, chords, [], u, v))
        ab = sorted(abs(b) for b in biases)
        mean = sum(biases) / len(biases)
        meanabs = sum(ab) / len(ab)
        mx = ab[-1]
        print(f"  k={k:3d}: mean bias {mean:+.5f}  E|bias| {meanabs:.5f}  "
              f"max {mx:.5f}")


def random_multicycle(kvals=(6, 10), trials=300, n=40, seed=7):
    """pi0 random permutation of [n] (random cycle type), u,v random,
    chords random disjoint pairs; frozen bridges."""
    rng = random.Random(seed)
    print(f"multi-cycle ensemble, n={n}")
    for k in kvals:
        biases = []
        samecyc = 0
        for t in range(trials):
            perm = list(range(n))
            rng.shuffle(perm)
            cyc = cycles_of(perm)
            pts = list(range(n))
            rng.shuffle(pts)
            u, v = pts[0], pts[1]
            rest = pts[2:]
            chords = [(rest[2 * i], rest[2 * i + 1]) for i in range(k)]
            circle, bridges, tot = multi_cycle_augment(cyc, n)
            biases.append(bias_exact(circle, chords, bridges, u, v))
            if same_cycle(perm, u, v):
                samecyc += 1
        ab = sorted(abs(b) for b in biases)
        mean = sum(biases) / len(biases)
        meanabs = sum(ab) / len(ab)
        mx = ab[-1]
        print(f"  k={k:3d}: mean bias {mean:+.5f}  E|bias| {meanabs:.5f}  "
              f"max {mx:.5f}   (u~v in pi0: {samecyc}/{trials})")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "struct"):
        structured()
    if which in ("all", "geom"):
        random_geometric()
    if which in ("all", "short"):
        random_geometric_short_arc()
    if which in ("all", "multi"):
        random_multicycle()
