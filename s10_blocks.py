#!/usr/bin/env python3
"""s10_blocks.py -- exact pgf of the block count B_t of the uniform split-merge
chain started from ONE block (continuum limit of #cycles of (n-cycle)*(t random
transpositions)), via hook characters:

    E[x^{B_t}] = sum_{k=0}^{n-1} (-1)^k (1-2k/(n-1))^t  s_{(n-k,1^k)}(1^x),
    s_{(n-k,1^k)}(1^x) = [x(x+1)...(x+n-k-1)] [(x-1)(x-2)...(x-k)] / ( n (n-k-1)! k! ).

For integer x >= 1 only k <= x-1 contribute, and the n->infinity limit is a
polynomial in t of degree x-1:   E[2^B] = 2+2t,  E[3^B] = 3+4t+2t^2, ...
usage: python3 s10_blocks.py n tmax [xmax]
"""
import sys
from fractions import Fraction as Fr
from math import factorial

def rising(x, m):
    p = Fr(1)
    for i in range(m): p *= (x + i)
    return p

def falling(x, m):
    p = Fr(1)
    for i in range(m): p *= (x - i)
    return p

def hook_schur_ones(n, k, x):
    """s_{(n-k,1^k)}(1^x) for integer/rational x (hook-content formula)."""
    return rising(Fr(x), n - k) * falling(Fr(x) - 1, k) / (n * factorial(n - k - 1) * factorial(k))

def pgf(n, t, x):
    tot = Fr(0)
    kmax = n - 1
    if float(x) == int(x) and int(x) >= 1: kmax = min(n - 1, int(x) - 1)
    for k in range(kmax + 1):
        r = Fr(n - 1 - 2 * k, n - 1)
        tot += (-1) ** k * r ** t * hook_schur_ones(n, k, x)
    return tot

def continuum_poly(x, tmax):
    """Continuum polynomial in t for integer x, obtained exactly by taking n large
    enough that the rational function of n is evaluated at several n and the
    limit read off by fitting a polynomial in 1/n (degree <= x-1+tmax)."""
    # exact limit: use n as a symbol via evaluation at many n and Richardson? Simpler:
    # the finite-n expression is a polynomial in 1/(n-1) of bounded degree; evaluate at
    # D+1 values and Lagrange-interpolate at 1/(n-1)=0.
    D = x - 1 + tmax + 2
    ns = [200 + 50 * i for i in range(D + 1)]
    out = []
    for t in range(tmax + 1):
        pts = [(Fr(1, n - 1), pgf(n, t, x)) for n in ns]
        # Lagrange at 0
        val = Fr(0)
        for i, (xi, yi) in enumerate(pts):
            L = Fr(1)
            for j, (xj, _) in enumerate(pts):
                if i != j: L *= (0 - xj) / (xi - xj)
            val += yi * L
        out.append(val)
    return out

if __name__ == "__main__":
    n = int(sys.argv[1]); tmax = int(sys.argv[2]); xmax = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    for x in range(2, xmax + 1):
        vals = continuum_poly(x, tmax)
        print(f"x={x}: continuum E[x^B_t] for t=0..{tmax}: " + ", ".join(str(v) for v in vals))
        # finite n values
        print(f"     finite n={n}: " + ", ".join(f"{float(pgf(n,t,x)):.6f}" for t in range(tmax + 1)))
