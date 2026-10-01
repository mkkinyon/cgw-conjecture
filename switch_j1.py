"""Exact j=1 bracket structure: Br({m}) for every admissible cell m,
as exact rationals, to find the organizing principle of the signed sum.

Also: exact signed sums S_J = sum_{j<=J} (-1)^j sum_{M in Mj} R0(M) Br(M)
for J = 0,1 (and 2 at small n), compared against exact 2*BrG (brute force
via per-permanent A(tau) at n=8,9).
"""
from fractions import Fraction
from math import factorial
import itertools, sys

from switch import (path_mp, cycle_mp, poly_mul, menage, board_count,
                    make_boards)


def per_ryser(A):
    """Permanent of 0/1 matrix A (list of row bitmasks? use simple Ryser)."""
    n = len(A)
    # A: list of lists
    total = 0
    for S in range(1, 1 << n):
        prod = 1
        for i in range(n):
            s = 0
            for j in range(n):
                if S >> j & 1:
                    s += A[i][j]
            prod *= s
            if prod == 0:
                break
        if prod:
            k = bin(S).count('1')
            total += (-1) ** (n - k) * prod
    return total


def brute_BrG(n):
    """Exact BrG = G(d1)+G(d2)-G(e1)-G(e2) via A(tau) permanents."""
    B2, B2p, b0, Dsig, Dsigp = make_boards(n)
    d1, d2 = Dsig
    e1, e2 = Dsigp
    good = [t for t in itertools.permutations(range(n))
            if all((i, t[i]) not in b0 for i in range(n))]
    print(f"  n={n}: |A0|={len(good)}")
    def A(t):
        mat = [[0 if ((i, c) in b0 or t[i] == c) else 1 for c in range(n)]
               for i in range(n)]
        return per_ryser(mat)
    G = {d: 0 for d in (d1, d2, e1, e2)}
    for t in good:
        a = None
        for d in (d1, d2, e1, e2):
            if t[d[0]] == d[1]:
                if a is None:
                    a = A(t)
                G[d] += a
    return G[d1] + G[d2] - G[e1] - G[e2], len(good)


def all_cells(n, b0):
    return [(i, c) for i in range(n) for c in range(n) if (i, c) not in b0]


def bracket(b0, n, M, Dsig, Dsigp):
    d1, d2 = Dsig
    e1, e2 = Dsigp
    def pin(d):
        if d in M:
            return 0  # incompatible (same cell pinned twice) -- handle outside
        return board_count(b0, n, tuple(M) + (d,))
    return pin(d1) + pin(d2) - pin(e1) - pin(e2)


def j1_scan(n):
    """Print Br({m}) * Mn / (R0({m}) * M_{n-4}) for all cells m, grouped."""
    B2, B2p, b0, Dsig, Dsigp = make_boards(n)
    Mn4 = menage(n - 4)
    groups = {}
    for m in all_cells(n, b0):
        br = bracket(b0, n, (m,), Dsig, Dsigp)
        r0 = board_count(b0, n, (m,))
        rat = Fraction(br) / (Fraction(r0) * Fraction(Mn4, menage(n)))
        i, c = m
        # classify cell: which lines does it touch?
        tag = []
        if i == 0: tag.append('r0')       # P1 interior row
        if c == 1: tag.append('c1')       # P1 interior col
        if i in (1, n - 1): tag.append('switchrow')
        if c in (0, 2): tag.append('switchcol')
        # omega-distance along sigma for context
        tag = ','.join(tag) if tag else 'generic'
        groups.setdefault(tag, []).append((float(rat), m))
    for tag in sorted(groups):
        vals = sorted(groups[tag])
        lo, hi = vals[0], vals[-1]
        mean = sum(v for v, _ in vals) / len(vals)
        print(f"  [{tag}] ({len(vals)} cells) min={lo[0]:.4f}@{lo[1]} "
              f"max={hi[0]:.4f}@{hi[1]} mean={mean:.4f}")
    return groups


def signed_sums(n, Jmax=2):
    """Exact signed sums sum_j (-1)^j sum_M R0(M) Br(M) vs brute 2BrG-ish.
    Note: M ranges over ALL matchings avoiding b0 (including through switch
    cells d/e themselves: then some pins are incompatible -> pin=0 unless
    handled).  We handle M containing a switch cell separately: if d in M,
    the G(d)-IE term is R0(M) (tau pin absorbed)."""
    B2, B2p, b0, Dsig, Dsigp = make_boards(n)
    d1, d2 = Dsig
    e1, e2 = Dsigp
    cells = all_cells(n, b0)
    total_by_j = {}
    def G_term(d, M):
        # IE term for G(d): #(tau >= M u {d}) * #(eta >= M): need M u {d} matching.
        Ms = set(M)
        if d in Ms:
            t = board_count(b0, n, tuple(M))
        else:
            if any(x[0] == d[0] or x[1] == d[1] for x in M):
                return 0
            t = board_count(b0, n, tuple(M) + (d,))
        return t * board_count(b0, n, tuple(M))
    from itertools import combinations
    for j in range(Jmax + 1):
        s = 0
        for Mc in combinations(cells, j):
            rows = set(x[0] for x in Mc)
            cols = set(x[1] for x in Mc)
            if len(rows) != j or len(cols) != j:
                continue
            term = (G_term(d1, Mc) + G_term(d2, Mc)
                    - G_term(e1, Mc) - G_term(e2, Mc))
            s += (-1) ** j * term
        total_by_j[j] = s
        print(f"  n={n} j={j}: signed term={s}  running={sum(total_by_j.values())}")
    return total_by_j


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'scan':
        for n in [int(a) for a in sys.argv[2:]]:
            print(f"j=1 scan at n={n} (ratio = Br*Mn/(R0*M_(n-4))):")
            j1_scan(n)
    elif cmd == 'signed':
        n = int(sys.argv[2]); J = int(sys.argv[3]) if len(sys.argv) > 3 else 2
        signed_sums(n, J)
    elif cmd == 'brg':
        for n in [int(a) for a in sys.argv[2:]]:
            brg, sz = brute_BrG(n)
            print(f"  n={n}: exact BrG={brg}")


def anatomy(n):
    """Exact anatomy of BrG: masses m+/-, conditional means E+-[A], and
    comparison with M_{n-4} * E[A]."""
    B2, B2p, b0, Dsig, Dsigp = make_boards(n)
    d1, d2 = Dsig; e1, e2 = Dsigp
    good = [t for t in itertools.permutations(range(n))
            if all((i, t[i]) not in b0 for i in range(n))]
    N0 = len(good)
    def A(t):
        mat = [[0 if ((i, c) in b0 or t[i] == c) else 1 for c in range(n)]
               for i in range(n)]
        return per_ryser(mat)
    sumA = 0; Gp = 0; Gm = 0; mp = 0; mm = 0
    for t in good:
        a = A(t)
        sumA += a
        wp = (1 if t[d1[0]] == d1[1] else 0) + (1 if t[d2[0]] == d2[1] else 0)
        wm = (1 if t[e1[0]] == e1[1] else 0) + (1 if t[e2[0]] == e2[1] else 0)
        Gp += wp * a; Gm += wm * a; mp += wp; mm += wm
    EA = Fraction(sumA, N0)
    EpA = Fraction(Gp, mp); EmA = Fraction(Gm, mm)
    BrG = Gp - Gm
    Mn4 = menage(n - 4)
    print(f"n={n}: N0={N0}  m+={mp} m-={mm}  m+-m-={mp-mm}  M_(n-4)={Mn4}  "
          f"{'OK' if mp-mm == Mn4 else 'MISMATCH'}")
    print(f"  E[A]={float(EA):.6f}  E+[A]={float(EpA):.6f}  E-[A]={float(EmA):.6f}")
    print(f"  (E+-E-)/E = {float((EpA-EmA)/EA):.3e}   (E- - E)/E = {float((EmA-EA)/EA):.3e}")
    print(f"  BrG={BrG}   M_(n-4)*E[A]={float(Mn4*EA):.2f}   ratio={float(Fraction(BrG)/(Mn4*EA)):.6f}")
    print(f"  N4_0=sum A={sumA}")
