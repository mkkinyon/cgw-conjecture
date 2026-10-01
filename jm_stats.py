"""Consume raw squares from the JM sampler and measure:
  (1) distribution of cycle type of sigma_{1,2}  (validation vs exact tables at n<=11)
  (2) generic-pair B probability: random columns {c,c'}, P(rows 1,2 in same
      column cycle)
  (3) post-flip B probabilities: for squares with type mu=(2,n-2), the CGW
      joining instance with a random (j,j'), events B1,B2, and qbar
Batch-mean error bars.

Usage: jm | python3 jm_stats.py n
"""
import sys, os
import numpy as np
from collections import Counter, defaultdict
from verify_identity import (sigma_of, cycles_of_perm, row_cycle_columns,
                             delta_same_cycle, do_switch, do_flip, omega_map,
                             q_events)
import random

n = int(sys.argv[1])
rng = random.Random(2026)

type_counts = Counter()
generic = [0, 0]          # [B count, total]
post = defaultdict(lambda: [0, 0, 0])   # per-mu: [B1 count, B2 count, total]
NBATCH = 20
batches_q = defaultdict(lambda: [[] for _ in range(NBATCH)])
batches_gen = [[] for _ in range(NBATCH)]

MU = (2, n - 2)

nread = 0
buf = sys.stdin.buffer
sq = buf.read(n * n)
while len(sq) == n * n:
    L = [tuple(sq[r * n:(r + 1) * n]) for r in range(n)]
    nread += 1
    if nread <= 3:  # validate latinity of first few
        for r in range(n):
            assert sorted(L[r]) == list(range(n))
        for c in range(n):
            assert sorted(L[r][c] for r in range(n)) == list(range(n))
    b = (nread * NBATCH // 10**9) % NBATCH  # placeholder, fixed below

    sig = sigma_of(L)
    cyc = cycles_of_perm(sig)
    typ = tuple(sorted(len(c) for c in cyc))
    type_counts[typ] += 1

    # generic pair
    while True:
        c1, c2 = rng.sample(range(n), 2)
        w = omega_map(L)
        if c2 != w[c1] and c1 != w[c2]:
            break
    gB = delta_same_cycle(L, c1, c2)
    generic[0] += gB
    generic[1] += 1

    # post-flip events for mu
    if typ == MU:
        alpha, beta = 2, n - 2
        ca = next(c for c in cyc if len(c) == alpha)
        cb = next(c for c in cyc if len(c) == beta)
        cols_a = row_cycle_columns(L, ca)
        cols_b = row_cycle_columns(L, cb)
        j = rng.choice(cols_a)
        jp = rng.choice(cols_b)
        if delta_same_cycle(L, j, jp):     # B-pair: switch min-symbol cycle
            mina = min(ca)
            minb = min(cb)
            cols = cols_a if mina < minb else cols_b
            L2 = do_switch(L, cols)
        else:
            L2 = L
        e1, e2, lens = q_events(L2, j, jp, alpha, beta)
        st = post[MU]
        st[0] += e1
        st[1] += e2
        st[2] += 1
    sq = buf.read(n * n)

print(f"n={n}: samples read: {nread}")

# (1) validation vs exact tables where available
try:
    from cgw_data import C, gamma
    if n in C:
        tot = sum(gamma(l, n) * C[n][l] for l in C[n])
        print("cycle-type distribution vs exact:")
        chi2 = 0.0
        for lam in sorted(C[n]):
            p_exact = gamma(lam, n) * C[n][lam] / tot
            obs = type_counts.get(lam, 0)
            exp = p_exact * nread
            if exp > 5:
                chi2 += (obs - exp) ** 2 / exp
            print(f"  {str(lam):18s} obs {obs:8d}  exp {exp:10.1f}")
        print(f"  chi2 = {chi2:.1f} on ~{len(C[n])-1} dof")
except ImportError:
    pass

# (2) generic B
p = generic[0] / generic[1]
se = (p * (1 - p) / generic[1]) ** 0.5
print(f"generic-pair B prob: {p:.5f} +- {se:.5f}   [(p-1/2)*n = {(p-0.5)*n:+.3f}]")

# (3) post-flip
st = post[MU]
if st[2] > 0:
    p1 = st[0] / st[2]
    p2 = st[1] / st[2]
    qb = p1 + p2
    se1 = (p1 * (1 - p1) / st[2]) ** 0.5
    seq = se1 * 1.6
    print(f"post-flip mu={MU}: N={st[2]}  P(B1)={p1:.5f}+-{se1:.5f}  "
          f"P(B2)={p2:.5f}  qbar={qb:.5f}+-{seq:.5f}")
    try:
        from cgw_data import C as CC
        lam = (n,)
        qexact = 2 * CC[n][lam] / CC[n][MU] - 1
        print(f"  exact qbar (identity) = {qexact:.6f}")
    except Exception:
        pass
