"""Exact joint moments of (X12, X13, X23) in a uniform 3 x n Latin rectangle.

Normalisation: rows (id, sigma, tau); the uniform measure weights each pair
(sigma, tau) equally, sigma in D_n, tau discordant with id and sigma.
  X12 = #2-cycles of sigma,
  X13 = #2-cycles of tau,
  X23 = #2-cycles of tau o sigma^{-1}.

We compute exactly, using cycle-type sums with gamma weights and punctured
counts:
  E[X12], E[X13], E[X23]  (the last two must agree with each other by row
  exchangeability; X12 differs from them only through... in fact all three are
  exchangeable -- a nontrivial code check),
  E[X12 X13], E[X13 X23], E[X12 X23],
  E[(X13)_2] (second factorial moment).

Brute-force certification at n = 7.
"""
from fractions import Fraction
from math import factorial
import itertools, sys
from puncture import make_sigma
from cgw_data import gamma  # gamma(lam, n)
from layer4_production import fast_R


def partitions_no_ones(n, min_part=2):
    if n == 0:
        yield ()
        return
    for p in range(min_part, n + 1):
        if n - p == 0 or n - p >= p:
            for rest in partitions_no_ones(n - p, p):
                yield tuple(sorted((p,) + rest))


def ctx(n, lam):
    sigma = make_sigma(n, lam)
    sig_inv = [0] * n
    for i in range(n):
        sig_inv[sigma[i]] = i
    return sigma, sig_inv


def E_X13_given(n, sigma, sig_inv):
    """sum over admissible pairs {i,j} of R({(i,j),(j,i)})  (tau 2-cycles)."""
    tot = 0
    for i in range(n):
        for j in range(i + 1, n):
            if j in (i, sigma[i]) or i in (j, sigma[j]):
                continue
            tot += fast_R(n, sigma, sig_inv, ((i, j), (j, i)))
    return tot


def E_X23_given(n, sigma, sig_inv):
    """tau o sigma^{-1} has 2-cycle {x,y} <=> tau(sigma^{-1}x)=y, tau(sigma^{-1}y)=x.
    Cells ((si_x, y), (si_y, x)) with si_x = sigma^{-1}(x).  Admissible iff cells
    avoid B2 and form a matching."""
    tot = 0
    for x in range(n):
        for y in range(x + 1, n):
            ix, iy = sig_inv[x], sig_inv[y]
            # cells (ix, y), (iy, x); avoid B2: y not in {ix, sigma[ix]=x} etc.
            if y == ix or y == x:      # y == sigma[ix] means y == x
                continue
            if x == iy or x == y:
                continue
            # matching condition: rows ix != iy (true), cols y != x (true)
            tot += fast_R(n, sigma, sig_inv, ((ix, y), (iy, x)))
    return tot


def E2_X13_given(n, sigma, sig_inv):
    """sum over ordered pairs of DISTINCT tau-2-cycles => (X13)_2 numerator /2...
    We return sum over unordered pairs of disjoint admissible {i,j},{k,l}:
    contributes E[C(X13,2)]."""
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            if j == sigma[i] or i == sigma[j]:
                continue
            pairs.append((i, j))
    tot = 0
    for a in range(len(pairs)):
        i, j = pairs[a]
        for b in range(a + 1, len(pairs)):
            k, l = pairs[b]
            if len({i, j, k, l}) < 4:
                continue
            tot += fast_R(n, sigma, sig_inv, ((i, j), (j, i), (k, l), (l, k)))
    return tot


def E_X13X23_given(n, sigma, sig_inv):
    """sum over admissible (tau-2-cycle {i,j}) x (tausigma^{-1}-2-cycle {x,y}),
    cells forming a matching of size <= 4 (cells may coincide!).
    Coinciding cells are possible: (i,j)=(six, y) etc. -- handle by building
    the union of the cell set and checking matching validity."""
    tot = 0
    tpairs = []
    for i in range(n):
        for j in range(i + 1, n):
            if j == sigma[i] or i == sigma[j]:
                continue
            tpairs.append(((i, j), (j, i)))
    spairs = []
    for x in range(n):
        for y in range(x + 1, n):
            ix, iy = sig_inv[x], sig_inv[y]
            if y == ix or x == iy:
                continue
            spairs.append(((ix, y), (iy, x)))
    for A in tpairs:
        for B in spairs:
            cells = set(A) | set(B)
            rows = [c[0] for c in cells]
            cols = [c[1] for c in cells]
            if len(set(rows)) != len(cells) or len(set(cols)) != len(cells):
                continue  # not a matching: tau cannot contain both
            tot += fast_R(n, sigma, sig_inv, tuple(sorted(cells)))
    return tot


def moments(n, verbose=True):
    from collections import Counter
    Z = Fraction(0)
    acc = {k: Fraction(0) for k in
           ("X12", "X13", "X23", "X12X13", "X12X23", "CX13_2", "X13X23")}
    for lam in partitions_no_ones(n):
        g = gamma(lam, n)
        sigma, sig_inv = ctx(n, lam)
        N3 = fast_R(n, sigma, sig_inv, ())
        if N3 == 0:
            continue
        lam2 = Counter(lam).get(2, 0)
        Z += g * N3
        e13 = E_X13_given(n, sigma, sig_inv)
        e23 = E_X23_given(n, sigma, sig_inv)
        acc["X12"] += g * N3 * lam2
        acc["X13"] += g * e13
        acc["X23"] += g * e23
        acc["X12X13"] += g * e13 * lam2
        acc["X12X23"] += g * e23 * lam2
        acc["CX13_2"] += g * E2_X13_given(n, sigma, sig_inv)
        acc["X13X23"] += g * E_X13X23_given(n, sigma, sig_inv)
    out = {k: v / Z for k, v in acc.items()}
    if verbose:
        print(f"n={n}:")
        print(f"  E[X12] = {float(out['X12']):.6f}   E[X13] = {float(out['X13']):.6f}"
              f"   E[X23] = {float(out['X23']):.6f}   (Poisson: 0.5)")
        print(f"  E[X12X13] = {float(out['X12X13']):.6f}  E[X12X23] = {float(out['X12X23']):.6f}"
              f"  E[X13X23] = {float(out['X13X23']):.6f}   (indep: 0.25)")
        print(f"  E[C(X13,2)] = {float(out['CX13_2']):.6f}   (Poisson: 0.125)")
        c1213 = out['X12X13'] - out['X12'] * out['X13']
        c1323 = out['X13X23'] - out['X13'] * out['X23']
        print(f"  Cov(X12,X13) = {float(c1213):+.6f}   Cov(X13,X23) = {float(c1323):+.6f}")
    return out


def brute(n):
    """Brute force at small n over all (sigma, tau)."""
    from collections import Counter
    perms = list(itertools.permutations(range(n)))
    ders = [p for p in perms if all(p[i] != i for i in range(n))]
    acc = Counter()
    Z = 0
    for sigma in ders:
        sig_inv = [0] * n
        for i in range(n):
            sig_inv[sigma[i]] = i
        for tau in ders:
            if any(tau[i] == sigma[i] for i in range(n)):
                continue
            Z += 1
            def two_cycles(p):
                return sum(1 for i in range(n) for j in range(i + 1, n)
                           if p[i] == j and p[j] == i)
            x12 = two_cycles(sigma)
            x13 = two_cycles(tau)
            ts = [tau[sig_inv[x]] for x in range(n)]
            x23 = two_cycles(ts)
            acc["X12"] += x12; acc["X13"] += x13; acc["X23"] += x23
            acc["X12X13"] += x12 * x13; acc["X12X23"] += x12 * x23
            acc["X13X23"] += x13 * x23
            acc["CX13_2"] += x13 * (x13 - 1) // 2
    return {k: Fraction(v, Z) for k, v in acc.items()}


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "brute7":
        b = brute(7)
        m = moments(7, verbose=False)
        ok = all(b[k] == m[k] for k in b)
        print("n=7 brute vs toolkit:", "MATCH" if ok else "MISMATCH")
        for k in sorted(b):
            print(f"  {k}: brute {b[k]} toolkit {m[k]}")
    else:
        for n in [int(a) for a in sys.argv[1:]] or [8, 10, 12]:
            moments(n)
