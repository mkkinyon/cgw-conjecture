"""Layer-5 effect via the TRIPLE-IE series (TODO 3), types (n) vs (2,n-2).

Truncated triple-IE (see triple_ie.py for the identity and its verification):

  N5 = sum over matchings U of the complement of B2, and label maps
       T: U -> nonempty subsets of {1,2,3},  of
       (-1)^{sum |T(u)|}  R(A) R(B) R(C),
  A = {u: T(u) n {1,2} != 0}, B = {u: T(u) n {1,3} != 0},
  C = {u: T(u) n {2,3} != 0}.

Total order q = |M1|+|M2|+|M3| = sum_u |T(u)|.  S(Q) = partial sum over
q <= Q.  This file computes S(Q), Q <= 3, exactly (integer arithmetic) for
both cycle types, forms the truncated layer-5 effects

  r5(Q) = [S_b(Q)/N4_b] / [S_a(Q)/N4_a] - 1     (N4 exact for n <= 11,
                                                 else pair-IE J=3)

and validates against exact N5 at n = 6, 7 (triple_ie.py brute force).

Usage: python3 layer5_production.py 8 9 10 11 12 13 14
       python3 layer5_production.py validate     (n=6,7 vs brute force)
"""
from fractions import Fraction
from math import factorial, comb, log
import sys, time, math
from puncture import make_sigma
from layer4_production import fast_R, make_ctx, Ts, EXACT

# nonempty subsets of {1,2,3} as bitmasks 1..7, |T| = popcount
SUBSETS = [(t, bin(t).count('1')) for t in range(1, 8)]
# membership: A iff T n {1,2} (bits 0,1), B iff T n {1,3} (bits 0,2),
# C iff T n {2,3} (bits 1,2)
MEMB = {t: (bool(t & 0b011), bool(t & 0b101), bool(t & 0b110))
        for t in range(1, 8)}


def S_orders(n, lam, Q=3):
    """Exact by-order sums s[q], q = 0..Q, of the triple-IE series."""
    sigma, sig_inv, cells = make_ctx(n, lam)
    Rc = {}

    def R(sub):
        key = sub
        v = Rc.get(key)
        if v is None:
            v = fast_R(n, sigma, sig_inv, sub)
            Rc[key] = v
        return v

    s = [0] * (Q + 1)
    N0 = R(())
    s[0] = N0 ** 3

    # enumerate matchings U of size 1..Q
    L = len(cells)

    def labels(U):
        """yield (q, A, B, C) over label maps with sum|T| <= Q."""
        k = len(U)
        stack = [(0, 0, [], [], [])]
        while stack:
            idx, q, A, B, C = stack.pop()
            if idx == k:
                yield q, tuple(A), tuple(B), tuple(C)
                continue
            u = U[idx]
            for t, sz in SUBSETS:
                q2 = q + sz
                if q2 + (k - idx - 1) > Q:
                    continue
                a, b, c = MEMB[t]
                stack.append((idx + 1, q2,
                              A + [u] if a else A,
                              B + [u] if b else B,
                              C + [u] if c else C))

    def visit(U):
        for q, A, B, C in labels(U):
            s[q] += (-1) ** q * R(A) * R(B) * R(C)

    # sizes 1..Q (each cell contributes >= 1 to q)
    for a in range(L):
        ia, ca = cells[a]
        visit((cells[a],))
        if Q >= 2:
            for b in range(a + 1, L):
                ib, cb = cells[b]
                if ib == ia or cb == ca:
                    continue
                visit((cells[a], cells[b]))
                if Q >= 3:
                    for d in range(b + 1, L):
                        id_, cd = cells[d]
                        if id_ in (ia, ib) or cd in (ca, cb):
                            continue
                        visit((cells[a], cells[b], cells[d]))
    return s


def N4_exact_or_trunc(n, which):
    """which = 0 for type (n), 1 for (2, n-2)."""
    if n in EXACT:
        return EXACT[n][1 + 2 * which], True
    lam = (n,) if which == 0 else (2, n - 2)
    T = Ts(n, lam, 3)
    return sum((-1) ** j * T[j] for j in range(4)), False


def run(n, Q=3):
    t0 = time.time()
    sa = S_orders(n, (n,), Q)
    sb = S_orders(n, (2, n - 2), Q)
    n4a, exact4a = N4_exact_or_trunc(n, 0)
    n4b, _ = N4_exact_or_trunc(n, 1)
    print(f"n={n}  [{time.time()-t0:.0f}s]  N4 {'exact' if exact4a else 'pair-IE J=3'}")
    for q in range(Q + 1):
        Sa = sum(sa[:q + 1])
        Sb = sum(sb[:q + 1])
        eff = Fraction(Sb * n4a, Sa * n4b) - 1
        print(f"  S(Q={q}): a={Sa}  b={Sb}  r5(Q) = {float(eff):+.6e}"
              f"   r5*(n)_4 = {float(eff) * n*(n-1)*(n-2)*(n-3):+.4f}")
    return sa, sb


def validate():
    import triple_ie
    for n, lam_pairs, exacts in [
            (6, ((6,), (2, 4)), None),
            (7, ((7,), (2, 5)), None)]:
        for lam in lam_pairs:
            sigma = make_sigma(n, lam)
            exact = triple_ie.brute_N5(sigma)
            s = S_orders(n, lam, 3)
            partial = [sum(s[:q + 1]) for q in range(4)]
            print(f"n={n} lam={lam}: exact N5 = {exact}")
            for q, P in enumerate(partial):
                print(f"   S(Q={q}) = {P}   rel err {P/exact - 1:+.4e}")


if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "validate":
        validate()
    else:
        ns = [int(a) for a in args] or [8, 9, 10, 11]
        for n in ns:
            run(n)
