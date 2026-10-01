#!/usr/bin/env python3
"""Monte Carlo check of the block-count law P(B_t=m) of uniform split-merge from one
block against the Harer-Zagier law P(B_t = m) = [z^{t+1}] (2 artanh z)^m / (2 m!).
usage: python3 s10_blocksim.py tmax samples"""
import sys, random
from fractions import Fraction as Fr

def hz_law(t):
    # coefficients of artanh z = sum z^{2i+1}/(2i+1)
    N = t + 2
    at = [Fr(0)] * N
    for i in range(0, N, 1):
        if i % 2 == 1: at[i] = Fr(1, i)
    out = {}
    # (2 artanh z)^m / (2 m!)
    poly = [Fr(1)] + [Fr(0)] * (N - 1)
    fact = 1
    for m in range(1, t + 2):
        new = [Fr(0)] * N
        for i in range(N):
            if poly[i] == 0: continue
            for j in range(N - i):
                new[i + j] += poly[i] * 2 * at[j]
        poly = new; fact *= m
        c = poly[t + 1] / (2 * fact)
        if c: out[m] = c
    return out

def simulate(tmax, S):
    counts = {t: {} for t in range(tmax + 1)}
    for _ in range(S):
        blocks = [1.0]
        for t in range(1, tmax + 1):
            # two uniform cuts: pick blocks size-biased
            c1 = random.random(); c2 = random.random()
            acc = 0.0; i1 = i2 = None
            for i, b in enumerate(blocks):
                if i1 is None and c1 < acc + b: i1 = i
                if i2 is None and c2 < acc + b: i2 = i
                acc += b
            if i1 is None: i1 = len(blocks) - 1
            if i2 is None: i2 = len(blocks) - 1
            if i1 == i2:
                b = blocks[i1]; m = abs(random.random() - random.random()) * b
                blocks[i1] = b - m; blocks.append(m)
            else:
                b = blocks[i1] + blocks[i2]
                for i in sorted([i1, i2], reverse=True): del blocks[i]
                blocks.append(b)
            k = len(blocks); counts[t][k] = counts[t].get(k, 0) + 1
    return counts

if __name__ == "__main__":
    tmax = int(sys.argv[1]); S = int(sys.argv[2])
    counts = simulate(tmax, S)
    for t in range(1, tmax + 1):
        law = hz_law(t)
        e2 = sum(float(p) * 2 ** m for m, p in law.items()); e3 = sum(float(p) * 3 ** m for m, p in law.items())
        print(f"t={t}: HZ E[2^B]={e2:.3f} (2+2t={2+2*t}), E[3^B]={e3:.3f} (3+4t+2t^2={3+4*t+2*t*t})")
        for m in sorted(law):
            mc = counts[t].get(m, 0) / S
            print(f"   m={m}: HZ {float(law[m]):.5f}   MC {mc:.5f}")
