"""Layer analysis, k=3: exact computation of
    R3(sigma) = per(J - I - P_sigma)
= number of permutations discordant with both the identity and sigma
= number of choices of a third row given first two rows (id, sigma).

For sigma a derangement with cycle type lambda, the forbidden bipartite graph
B = I \cup P_sigma is a disjoint union of even cycles C_{2l}, one per l-cycle.
Then per(J - B) = sum_k (-1)^k m_k(B) (n-k)!  where m_k(B) = #k-matchings of B,
and the matching generating polynomial factorises over cycles:
    M_{C_m}(x) = sum_k  m/(m-k) * C(m-k, k) x^k.

We compute R3 exactly for all cycle types (or a family of types) and analyse
the dependence of log R3 on the cycle type.
"""
from fractions import Fraction
from math import factorial, comb, log
from functools import lru_cache
import itertools


@lru_cache(maxsize=None)
def cycle_matching_poly(m):
    """Matching polynomial coefficients [m_0, m_1, ...] of the cycle C_m (m>=3),
    and for m=2 interpret as a doubled edge? Not needed: cycles here have length
    2*l with l>=2, so m>=4."""
    assert m >= 3
    return tuple(m * comb(m - k, k) // (m - k) for k in range(m // 2 + 1))


def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def matchings_of_union(lam):
    """m_k for B = disjoint union of C_{2l}, l in lam."""
    poly = [1]
    for l in lam:
        poly = poly_mul(poly, list(cycle_matching_poly(2 * l)))
    return poly


def R3(lam):
    n = sum(lam)
    mk = matchings_of_union(lam)
    return sum((-1) ** k * mk[k] * factorial(n - k) for k in range(len(mk)))


def partitions_no_ones(n, min_part=2):
    if n == 0:
        yield ()
        return
    for p in range(min_part, n + 1):
        if n - p == 0 or n - p >= p:
            for rest in partitions_no_ones(n - p, p):
                yield (p,) + rest


def brute_R3(sigma):
    n = len(sigma)
    cnt = 0
    for tau in itertools.permutations(range(n)):
        if all(tau[i] != i and tau[i] != sigma[i] for i in range(n)):
            cnt += 1
    return cnt


if __name__ == "__main__":
    # sanity: brute force for n=6,7
    import random
    for n in [6, 7]:
        for lam in partitions_no_ones(n):
            # build a permutation with cycle type lam
            sigma = list(range(n))
            pos = 0
            for l in sorted(lam):
                cyc = list(range(pos, pos + l))
                for i in range(l):
                    sigma[cyc[i]] = cyc[(i + 1) % l]
                pos += l
            assert brute_R3(sigma) == R3(tuple(sorted(lam))), (n, lam)
    print("brute-force check OK for n=6,7")

    # per-2-cycle effect at the k=3 layer:
    # compare (2, n-2) vs (n): ratio - 1, scaled by n^4
    print("\n l=2 effect:  r = R3((2,n-2))/R3((n,)) - 1")
    for n in range(6, 41, 2):
        r = Fraction(R3((2, n - 2)), R3((n,))) - 1
        print(f"  n={n:3d}  r={float(r):+.6e}   r*n^4={float(r) * n**4:+.6f}")
    print("\n l=3 effect:  r = R3((3,n-3))/R3((n,)) - 1")
    for n in range(8, 41, 2):
        r = Fraction(R3((3, n - 3)), R3((n,))) - 1
        print(f"  n={n:3d}  r={float(r):+.6e}   r*n^6={float(r) * n**6:+.6f}")
    print("\n additivity check at n=30: resid of log R3 additive model")
    n = 30
    lams = list(partitions_no_ones(n))
    print(f"  {len(lams)} cycle types at n={n}")
    # fit log R3 = c0 + sum_l h(l) lam_l  with gauge h(n)=0, using numpy lstsq
    import numpy as np
    from collections import Counter
    import mpmath
    mpmath.mp.dps = 60
    parts = sorted({p for lam in lams for p in lam if p != n})
    A, y = [], []
    base = mpmath.log(mpmath.mpf(R3((n,))))
    for lam in lams:
        c = Counter(lam)
        A.append([1.0] + [float(c.get(p, 0)) for p in parts])
        y.append(float(mpmath.log(mpmath.mpf(R3(lam))) - base))
    A = np.array(A); y = np.array(y)
    sol, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ sol
    h = {p: sol[1 + i] for i, p in enumerate(parts)}
    print(f"  max|resid| = {np.abs(resid).max():.3e}")
    for p in [2, 3, 4, 5, 6]:
        print(f"  h({p}) = {h[p]:+.6e}   h({p})*n^(2*{p}) = {h[p] * float(n)**(2*p):+.4f}")
