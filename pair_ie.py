"""The pair inclusion-exclusion identity:

    N4(sigma) := #{4 x n Latin rectangles with rows 1,2 = (id, sigma)}
              = #{(tau,eta): both avoid B2, tau(i) != eta(i) for all i}
              = sum_j (-1)^j  sum_{M in M_j(complement of B2)}  R(M;sigma)^2

(inclusion-exclusion over the set of common cells of tau and eta; the pairs
containing a given common j-set M number R(M)^2 by independence).

This file: (1) verify the identity at n = 6, 7 by brute force;
(2) exact T_j = sum_M R(M)^2 for j <= 3 at general n for sigma of type
(n) and (2, n-2), by direct enumeration of puncture positions.
"""
from math import factorial
from puncture import punctured_count, make_sigma
import itertools, sys


def brute_N4(sigma):
    n = len(sigma)
    disc = [t for t in itertools.permutations(range(n))
            if all(t[i] != i and t[i] != sigma[i] for i in range(n))]
    return sum(1 for t in disc for e in disc
               if all(t[i] != e[i] for i in range(n)))


def all_matchings_sum(sigma):
    """sum_j (-1)^j sum_M R(M)^2 over ALL matchings M of the complement."""
    n = len(sigma)
    cells = [(i, c) for i in range(n) for c in range(n)
             if c != i and c != sigma[i]]
    total = 0

    def rec(start, M, usedr, usedc):
        nonlocal total
        R = punctured_count(sigma, M)
        total += (-1) ** len(M) * R * R
        for k in range(start, len(cells)):
            i, c = cells[k]
            if not (usedr >> i) & 1 and not (usedc >> c) & 1:
                M.append((i, c))
                rec(k + 1, M, usedr | (1 << i), usedc | (1 << c))
                M.pop()

    rec(0, [], 0, 0)
    return total


def Tj(sigma, j):
    """T_j = sum over j-matchings M of complement of R(M)^2 (exact)."""
    n = len(sigma)
    cells = [(i, c) for i in range(n) for c in range(n)
             if c != i and c != sigma[i]]
    if j == 0:
        R = punctured_count(sigma, [])
        return R * R
    tot = 0
    if j == 1:
        for cell in cells:
            R = punctured_count(sigma, [cell])
            tot += R * R
        return tot
    if j == 2:
        for a in range(len(cells)):
            i1, c1 = cells[a]
            for b in range(a + 1, len(cells)):
                i2, c2 = cells[b]
                if i1 != i2 and c1 != c2:
                    R = punctured_count(sigma, [cells[a], cells[b]])
                    tot += R * R
        return tot
    if j == 3:
        for a in range(len(cells)):
            i1, c1 = cells[a]
            for b in range(a + 1, len(cells)):
                i2, c2 = cells[b]
                if i1 == i2 or c1 == c2:
                    continue
                for d in range(b + 1, len(cells)):
                    i3, c3 = cells[d]
                    if i3 in (i1, i2) or c3 in (c1, c2):
                        continue
                    R = punctured_count(sigma, [cells[a], cells[b], cells[d]])
                    tot += R * R
        return tot
    raise ValueError


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "verify"
    if which == "verify":
        for n, lam in [(6, (6,)), (6, (2, 4)), (7, (7,)), (7, (2, 5))]:
            sigma = make_sigma(n, lam)
            b = brute_N4(sigma)
            s = all_matchings_sum(sigma)
            print(f"n={n} lam={lam}: brute N4 = {b}, pair-IE sum = {s}, "
                  f"{'MATCH' if b == s else 'MISMATCH'}")
