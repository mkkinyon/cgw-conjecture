"""The TRIPLE inclusion-exclusion identity (layer 5 / TODO 3).

    N5(sigma) := #{5 x n Latin rectangles with rows 1,2 = (id, sigma)}
              = #{(tau,eta,theta): all avoid B2, pairwise cell-disjoint}

Derivation.  For each of the three pairs, expand the disjointness indicator
by inclusion-exclusion over subsets of the pair's coincidence-cell set:

  N5 = sum over (M1,M2,M3) of (-1)^{|M1|+|M2|+|M3|}
         R(M1 u M2) R(M1 u M3) R(M2 u M3)                      (M-form)

where M1 carries tau^eta coincidences, M2 tau^theta, M3 eta^theta, each Mi
is a matching of the complement of B2, and COMPATIBILITY is required:
each pairwise union M_i u M_j must itself be a matching (shared cells
allowed; if two M's meet the same row/col they must do so in the SAME cell).

KEY RESUMMATION (this file, verified): pairwise compatibility is equivalent
to U = M1 u M2 u M3 being a matching.  Each cell u in U carries a nonempty
label set T(u) <= {1,2,3} (which Mi contain it).  The three R-arguments
A = M1uM2 (tau's cells), B = M1uM3 (eta's), C = M2uM3 (theta's) depend only
on the membership pattern; summing the 7 label sets per cell gives per-cell
patterns   (in A,B only): weight -1   [T={1}]
           (in A,C only): weight -1   [T={2}]
           (in B,C only): weight -1   [T={3}]
           (in A,B,C):    weight +2   [T in {12,13,23} minus {123}: +3-1]
Hence the exact resummed form

  N5 = sum over matchings U of  sum over partitions of U into
       (Uab, Uac, Ubc, Uabc) of  (-1)^{|Uab|+|Uac|+|Ubc|} 2^{|Uabc|}
       R(Uab u Uac u Uabc) R(Uab u Ubc u Uabc) R(Uac u Ubc u Uabc). (U-form)

Every coincidence cell lies in >= 2 of the three constraint sets A,B,C.

Verification: brute force vs U-form (and, at n=5, vs the raw M-form too).
Usage: python3 triple_ie.py verify         (n=5, n=6 both types; ~1 min)
       python3 triple_ie.py verify7        (adds n=7 type (7) and (2,5))
"""
from math import factorial
from puncture import punctured_count, make_sigma
import itertools, sys

_Rcache = {}


def R(sigma, cells):
    key = (id(sigma), cells)
    v = _Rcache.get(key)
    if v is None:
        v = punctured_count(sigma, list(cells))
        _Rcache[key] = v
    return v


def complement_cells(sigma):
    n = len(sigma)
    return [(i, c) for i in range(n) for c in range(n)
            if c != i and c != sigma[i]]


def all_matchings(sigma):
    """All matchings (as tuples, cells in listing order) of the complement."""
    cells = complement_cells(sigma)
    out = []

    def rec(start, M, usedr, usedc):
        out.append(tuple(M))
        for k in range(start, len(cells)):
            i, c = cells[k]
            if not (usedr >> i) & 1 and not (usedc >> c) & 1:
                M.append((i, c))
                rec(k + 1, M, usedr | (1 << i), usedc | (1 << c))
                M.pop()

    rec(0, [], 0, 0)
    return out


def u_form_sum(sigma, maxorder=None, report_by_order=False):
    """The resummed identity.  Total order = |M1|+|M2|+|M3| =
    #2-label cells * 1 ... careful: order(U-partition) = |Uab|+|Uac|+|Ubc|
    + 3|Uabc| counted with multiplicity: |M1|+|M2|+|M3| = sum |T(u)| =
    (#2-cells)*1? No: T={1} is ONE M containing it -> contributes 1.
    Pattern ab <-> T={1}: order 1.  Pattern abc <-> |T| in {2,2,2,3}.
    We track order = |M1|+|M2|+|M3| only approximately for abc cells
    (mixed 2/3); for truncation purposes we use ord(ab/ac/bc)=1,
    ord(abc)=2 (the minimal representative).  maxorder=None: full sum."""
    total = 0
    by_order = {}
    for U in all_matchings(sigma):
        k = len(U)
        if maxorder is not None and k > maxorder:
            # each cell contributes >= 1 to the order
            continue
        # iterate over assignments of each cell to one of 4 classes
        for assign in itertools.product(range(4), repeat=k):
            # classes 0=ab,1=ac,2=bc,3=abc
            order = sum(1 if a < 3 else 2 for a in assign)
            if maxorder is not None and order > maxorder:
                continue
            A = tuple(u for u, a in zip(U, assign) if a in (0, 1, 3))
            B = tuple(u for u, a in zip(U, assign) if a in (0, 2, 3))
            C = tuple(u for u, a in zip(U, assign) if a in (1, 2, 3))
            n2 = sum(1 for a in assign if a < 3)
            n3 = k - n2
            w = (-1) ** n2 * 2 ** n3
            term = w * R(sigma, A) * R(sigma, B) * R(sigma, C)
            total += term
            by_order[order] = by_order.get(order, 0) + term
    if report_by_order:
        return total, by_order
    return total


def m_form_sum(sigma):
    """The raw M-form (slow; use only at n=5): enumerate compatible triples
    (M1,M2,M3) directly."""
    Ms = all_matchings(sigma)
    total = 0

    def compatible(Ma, Mb):
        rows = {}
        cols = {}
        for (i, c) in Ma:
            rows[i] = c
            cols[c] = i
        for (i, c) in Mb:
            if rows.get(i, c) != c:
                return False
            if cols.get(c, i) != i:
                return False
        return True

    for M1 in Ms:
        for M2 in Ms:
            if not compatible(M1, M2):
                continue
            U12 = tuple(sorted(set(M1) | set(M2)))
            for M3 in Ms:
                if not (compatible(M1, M3) and compatible(M2, M3)):
                    continue
                U13 = tuple(sorted(set(M1) | set(M3)))
                U23 = tuple(sorted(set(M2) | set(M3)))
                s = (-1) ** (len(M1) + len(M2) + len(M3))
                total += s * R(sigma, tuple(sorted(U12))) * \
                    R(sigma, U13) * R(sigma, U23)
    return total


def brute_N5(sigma):
    n = len(sigma)
    disc = [t for t in itertools.permutations(range(n))
            if all(t[i] != i and t[i] != sigma[i] for i in range(n))]
    masks = []
    for t in disc:
        m = 0
        for i in range(n):
            m |= 1 << (i * n + t[i])
        masks.append(m)
    N = len(masks)
    total = 0
    for a in range(N):
        ma = masks[a]
        for b in range(a + 1, N):
            mb = masks[b]
            if ma & mb:
                continue
            mab = ma | mb
            for c in range(b + 1, N):
                if not (mab & masks[c]):
                    total += 1
    return 6 * total  # ordered triples


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "verify"
    cases = [(5, (5,)), (6, (6,)), (6, (2, 4))]
    if which == "verify7":
        cases += [(7, (7,)), (7, (2, 5))]
    for n, lam in cases:
        sigma = make_sigma(n, lam)
        b = brute_N5(sigma)
        u = u_form_sum(sigma)
        line = f"n={n} lam={lam}: brute N5 = {b}, U-form = {u}"
        ok = (b == u)
        if n == 5:
            m = m_form_sum(sigma)
            line += f", M-form = {m}"
            ok = ok and (m == b)
        print(line + ("  MATCH" if ok else "  MISMATCH"))
