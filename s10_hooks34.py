#!/usr/bin/env python3
"""s10_hooks34.py -- hook inner products <chi_{(n-k,1^k)}, sum_i c_i^r> for r=2,3,4
(brute force over partitions, exact), used to guess closed forms; and the exact
E[sum p_i^r(m)] along the chain from one block at finite n.
usage: python3 s10_hooks34.py n M"""
import sys
from fractions import Fraction as Fr
sys.path.insert(0, '.')
from s9_hooks import partitions, hook_chars, zmu

def inner_r(n, r):
    ip = [Fr(0)] * n
    for mu in partitions(n):
        ch = hook_chars(mu, n); F = sum(x ** r for x in mu); z = zmu(mu)
        for k in range(n): ip[k] += Fr(ch[k] * F, z)
    return ip

def E_sum(n, M, ip, r):
    out = []
    for m in range(M + 1):
        s = Fr(0)
        for k in range(n):
            s += (-1) ** k * Fr(n - 1 - 2 * k, n - 1) ** m * ip[k]
        out.append(s / n ** r)
    return out

if __name__ == '__main__':
    n = int(sys.argv[1]); M = int(sys.argv[2])
    for r in (2, 3, 4):
        ip = inner_r(n, r)
        print(f"r={r} n={n}: (-1)^k <chi_k, sum c^{r}> for k=0..n-1:")
        print("   ", [str((-1) ** k * ip[k]) for k in range(n)])
        vals = E_sum(n, M, ip, r)
        print("    E[sum p^%d(m)] m=0..%d:" % (r, M), ", ".join(f"{float(v):.6f}" for v in vals))
