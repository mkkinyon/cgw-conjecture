"""Switch calculus for the FIFTH layer (TODO 3): r5 = r3 (1 + O(1/n)).

Same adjacent representatives and pins as switch.py (layer 4):
sigma = (0 1 ... n-1), sigma' = (0 1)(2 3 ... n-1); common board0 =
two paths; pins d1=(1,2), d2=(n-1,0) (sigma's), e1=(1,0), e2=(n-1,2)
(sigma's').  All 0-indexed.

N5(sigma) = # ordered triples (tau,eta,theta) of B2-avoiding permutations,
pairwise cell-disjoint.  IE over the two extra cells {d1,d2} on board0:

  N5(sigma) = P0 - 3 G5(d1) - 3 G5(d2) + 3 G25(d1,d2) + 6 H5(d1,d2)

where all counts are over triples avoiding board0, pairwise disjoint:
  P0        = # triples                       (no pin)
  G5(z)     = # triples with tau through z
  G25(z,z') = # triples with tau through BOTH z, z'
  H5(z,z')  = # triples with tau through z and eta through z'
(3 = choice of the pinned row; 6 = ordered choice of two pinned rows.)

  DeltaN5 := N5(sigma') - N5(sigma)
           = 3 BrG5 + 3 [G25(e1,e2) - G25(d1,d2)] + 6 [H5(e1,e2) - H5(d1,d2)]
  BrG5    := G5(d1) + G5(d2) - G5(e1) - G5(e2).

STRUCTURE (the bridge to the layer-4 machinery): expanding G5(z) by the
triple-IE (triple_ie.py) with the pin in tau's factor,

  G5(z) = sum over compatible (M1,M2,M3) of (-1)^{|M1|+|M2|+|M3|}
          R0(A u {z}) R0(B) R0(C),   A = M1uM2, B = M1uM3, C = M2uM3,

so BrG5 = sum (-1)^... B_A . R0(B) R0(C)  with B_A the SAME four-term
layer-4 bracket in the A-slot.  The A = empty stratum (M1 = M2 = empty,
M3 free) resums by the PAIR-IE to an exact leading term:

  BrG5|_{A=empty} = M_{n-4} * P0^(4),   P0^(4) = # ordered disjoint pairs
                                                 avoiding board0.

Checks:
  check1: the IE decomposition and DeltaN5 identity, brute force n=6,7.
  check2: BrG5 vs M_{n-4} P0^(4) (the A=empty leading term), n=6,7[,8].
Usage: python3 switch5.py 1     ;  python3 switch5.py 2
"""
import itertools, sys
from switch import make_boards, board_count, menage


def good_perms(n, B):
    return [t for t in itertools.permutations(range(n))
            if all((i, t[i]) not in B for i in range(n))]


def disjoint(t, u):
    return all(a != b for a, b in zip(t, u))


def triples_count(rows, pred=None):
    """# ordered pairwise-disjoint triples (tau,eta,theta), tau restricted
    by pred (or all)."""
    taus = [t for t in rows if pred is None or pred(t)]
    tot = 0
    for t in taus:
        L = [u for u in rows if disjoint(t, u)]
        for i, u in enumerate(L):
            for v in L:
                if v is not u and disjoint(u, v):
                    tot += 1
    return tot


def check1(ns=(6, 7)):
    for n in ns:
        B2, B2p, b0, Dsig, Dsigp = make_boards(n)
        d1, d2 = Dsig
        e1, e2 = Dsigp
        g0 = good_perms(n, b0)

        def N5_of(B):
            g = [t for t in g0 if all((i, t[i]) not in B for i in range(n))]
            tot = 0
            for t in g:
                L = [u for u in g if disjoint(t, u)]
                for u in L:
                    for v in L:
                        if v is not u and disjoint(u, v):
                            tot += 1
            return tot

        def G5(z):
            return triples_count(g0, lambda t: t[z[0]] == z[1])

        def G25(z, zp):
            return triples_count(g0, lambda t: t[z[0]] == z[1] and t[zp[0]] == zp[1])

        def H5(z, zp):
            tot = 0
            for t in g0:
                if t[z[0]] != z[1]:
                    continue
                L = [u for u in g0 if disjoint(t, u)]
                for u in L:
                    if u[zp[0]] != zp[1]:
                        continue
                    for v in L:
                        if v is not u and disjoint(u, v):
                            tot += 1
            return tot

        N5a, N5b = N5_of(B2), N5_of(B2p)
        P0 = triples_count(g0)
        # direct IE for N5(sigma)
        ie_a = P0 - 3 * G5(d1) - 3 * G5(d2) + 3 * G25(d1, d2) + 6 * H5(d1, d2)
        ie_b = P0 - 3 * G5(e1) - 3 * G5(e2) + 3 * G25(e1, e2) + 6 * H5(e1, e2)
        BrG5 = G5(d1) + G5(d2) - G5(e1) - G5(e2)
        dG2 = G25(e1, e2) - G25(d1, d2)
        dH = H5(e1, e2) - H5(d1, d2)
        lhs = N5b - N5a
        rhs = 3 * BrG5 + 3 * dG2 + 6 * dH
        print(f"n={n}: N5={N5a},{N5b}  IE: {'PASS' if (ie_a, ie_b) == (N5a, N5b) else 'FAIL'}"
              f"  DeltaN5={lhs} = 3BrG5+3dG2+6dH={rhs} "
              f"{'PASS' if lhs == rhs else 'FAIL'}")
        print(f"   pieces: 3BrG5={3*BrG5}  3dG2={3*dG2}  6dH={6*dH}")
        assert (ie_a, ie_b) == (N5a, N5b) and lhs == rhs


def check2(ns=(6, 7)):
    for n in ns:
        B2, B2p, b0, Dsig, Dsigp = make_boards(n)
        d1, d2 = Dsig
        e1, e2 = Dsigp
        g0 = good_perms(n, b0)

        def G5(z):
            return triples_count(g0, lambda t: t[z[0]] == z[1])

        BrG5 = G5(d1) + G5(d2) - G5(e1) - G5(e2)
        P04 = sum(1 for t in g0 for u in g0 if disjoint(t, u))
        P05 = triples_count(g0)
        N0 = len(g0)
        exact_lead = menage(n - 4) * P04   # the exact A=empty stratum
        # asymptotic normalisation: generic strata resum P04 -> P05/N0
        pred = menage(n - 4) * P05
        print(f"n={n}: BrG5={BrG5}  A=empty stratum M_(n-4)*P0^(4)={exact_lead}")
        if pred:
            print(f"   BrG5*N0/(M_(n-4)*P0^(5)) = {BrG5 * N0 / pred:.4f}"
                  f"   [layer-4 analogue was ~0.87-0.89 at n=8..11]")


def check3(ns=(6, 7)):
    """Triple-IE (U-form) evaluated ON THE BOARD board0 equals P0^(5), and
    the A=empty spectator weight w(empty) equals P0^(4) (pair-IE on
    board0).  These are identities (i),(ii) of the notes' eq:BrG5strata."""
    import itertools as it
    for n in ns:
        B2, B2p, b0, Dsig, Dsigp = make_boards(n)
        g0 = good_perms(n, b0)
        N0 = len(g0)
        P04 = sum(1 for t in g0 for u in g0 if disjoint(t, u))
        P05 = triples_count(g0)
        cells = [(i, c) for i in range(n) for c in range(n) if (i, c) not in b0]
        # enumerate matchings of the complement of board0
        Ms = []

        def rec(start, M, ur, uc):
            Ms.append(tuple(M))
            for k in range(start, len(cells)):
                i, c = cells[k]
                if not (ur >> i) & 1 and not (uc >> c) & 1:
                    M.append(cells[k])
                    rec(k + 1, M, ur | 1 << i, uc | 1 << c)
                    M.pop()

        rec(0, [], 0, 0)
        Rc = {}

        def R(sub):
            v = Rc.get(sub)
            if v is None:
                v = board_count(b0, n, sub)
                Rc[sub] = v
            return v

        # U-form sum
        tot = 0
        for U in Ms:
            k = len(U)
            for assign in it.product(range(4), repeat=k):
                A = tuple(u for u, a in zip(U, assign) if a in (0, 1, 3))
                B = tuple(u for u, a in zip(U, assign) if a in (0, 2, 3))
                C = tuple(u for u, a in zip(U, assign) if a in (1, 2, 3))
                n2 = sum(1 for a in assign if a < 3)
                tot += (-1) ** n2 * 2 ** (k - n2) * R(A) * R(B) * R(C)
        # pair-IE on board0 (w(empty))
        wempty = sum((-1) ** len(M) * R(M) ** 2 for M in Ms)
        print(f"n={n}: U-form on board0 = {tot}  P0^(5) = {P05}  "
              f"{'PASS' if tot == P05 else 'FAIL'};  w(empty) = {wempty}  "
              f"P0^(4) = {P04}  {'PASS' if wempty == P04 else 'FAIL'}")
        assert tot == P05 and wempty == P04


if __name__ == "__main__":
    args = sys.argv[1:] or ['1', '2']
    for a in args:
        if a == '1':
            check1()
        elif a == '2':
            check2()
        elif a == '3':
            check3()
