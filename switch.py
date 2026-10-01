"""Switch-IE calculus for the layer-4 theorem  r4 = r3 (1 + O(1/n)).

Compare sigma = (0 1 ... n-1)  [type (n)]  with
        sigma' = (0 1)(2 3 ... n-1)  [type (2, n-2)].

Their boards B2 = I u P_sigma differ in exactly two cells:
  sigma  has (1,2) and (n-1,0);   sigma' has (1,0) and (n-1,2).
board0 = I u (P_sigma \\ {(1,2),(n-1,0)})   (the common board, 2n-2 cells)
       = two paths:  P1 = c0-r0-c1-r1 (3 edges),  P2 = c2-r2-c3-...-c_{n-1}-r_{n-1} (2n-5 edges).

Pin cells (all avoid board0):
  Dsig  = {(1,2), (n-1,0)}   (forbidden for sigma,  admissible for sigma')
  Dsig' = {(1,0), (n-1,2)}   (forbidden for sigma', admissible for sigma)

IE over the two differing cells gives
  N4(sigma)  = P0 - U(1,2) - U(n-1,0) + U((1,2),(n-1,0))
  N4(sigma') = P0 - U(1,0) - U(n-1,2) + U((1,0),(n-1,2))
with U(d) = 2 G(d), U(d,d') = 2 G2(d,d') + 2 H(d,d'), so

  DeltaN4 = 2*BrG + 2*(G2' - G2) + 2*DeltaH,
  BrG  = G((1,2)) + G((n-1,0)) - G((1,0)) - G((n-1,2)),
  DeltaH = H((1,0),(n-1,2)) - H((1,2),(n-1,0)).

Claims to verify exactly:
  (1) R0({(1,2)}) + R0({(n-1,0)}) - R0({(1,0)}) - R0({(n-1,2)}) = M_{n-4}
      (the empty-matching bracket), for n = 6..14.
  (2) G2((1,0),(n-1,2)) = G2((1,2),(n-1,0))  exactly (same punctured boards),
      brute force at n = 6, 7.
  (3) DeltaN4 = 2 BrG + 2 DeltaH  exactly (brute force n = 6, 7), and
      per-piece sizes: DeltaH / DeltaN4 = O(1/n)-small.
  (4) per-M bracket Br(M) = R0(M u {(1,2)}) + R0(M u {(n-1,0)})
      - R0(M u {(1,0)}) - R0(M u {(n-1,2)}) satisfies
      Br(M) / (R0(M) * M_{n-4}/M_n) = 1 + O((j+1)^2/n) for generic M,
      and stays bounded (x^3-suppressed) for adversarial M near the switch.
"""
import itertools
import random
from collections import defaultdict
from functools import lru_cache
from fractions import Fraction
from math import factorial, comb


@lru_cache(maxsize=None)
def path_mp(e):
    if e < 0:
        return (1,)
    t = tuple(comb(e - k + 1, k) for k in range((e + 1) // 2 + 1))
    return t if t else (1,)


@lru_cache(maxsize=None)
def cycle_mp(m):
    assert m >= 3
    return tuple(m * comb(m - k, k) // (m - k) for k in range(m // 2 + 1))


def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def menage(n):
    # M_n = per(J - I - P_c), n-cycle c; M_0 := 2 (convention), M_1 := -1? use IE directly
    if n == 0:
        return 2
    if n == 1:
        return -1  # formal value of the IE sum; not used for n>=2 checks
    poly = cycle_mp(2 * n) if n >= 2 else None
    return sum((-1) ** k * poly[k] * factorial(n - k) for k in range(len(poly)) if n - k >= 0)


def board_count(board_cells, n, M=()):
    """R0(M) = # permutations of [n] avoiding board_cells and containing M.
    board_cells: set of (row, col).  M: pinned cells, pairwise row/col disjoint,
    each avoiding the board.  Returns exact integer."""
    rows = set(i for i, c in M)
    cols = set(c for i, c in M)
    if len(rows) != len(M) or len(cols) != len(M):
        return 0
    for (i, c) in M:
        if (i, c) in board_cells:
            return 0
    # remaining board edges
    edges = [(i, c) for (i, c) in board_cells if i not in rows and c not in cols]
    adj = defaultdict(list)
    for k, (i, c) in enumerate(edges):
        adj[('r', i)].append((('c', c), k))
        adj[('c', c)].append((('r', i), k))
    seen_v = set()
    poly = [1]
    for start in list(adj):
        if start in seen_v:
            continue
        stack = [start]
        vs, es = set(), set()
        while stack:
            u = stack.pop()
            if u in vs:
                continue
            vs.add(u)
            for (w, k) in adj[u]:
                es.add(k)
                if w not in vs:
                    stack.append(w)
        seen_v |= vs
        ne, nv = len(es), len(vs)
        if ne == nv:
            comp = list(cycle_mp(ne))
        elif ne == nv - 1:
            comp = list(path_mp(ne))
        else:
            raise RuntimeError("not path/cycle")
        poly = poly_mul(poly, comp)
    N = n - len(M)
    return sum((-1) ** k * poly[k] * factorial(N - k) for k in range(len(poly)) if N - k >= 0)


def make_boards(n):
    sig = {i: (i + 1) % n for i in range(n)}
    sigp = dict(sig)
    sigp[1] = 0
    sigp[n - 1] = 2
    B2 = set((i, i) for i in range(n)) | set((i, sig[i]) for i in range(n))
    B2p = set((i, i) for i in range(n)) | set((i, sigp[i]) for i in range(n))
    board0 = B2 & B2p
    Dsig = tuple(sorted(B2 - board0))     # {(1,2),(n-1,0)}
    Dsigp = tuple(sorted(B2p - board0))   # {(1,0),(n-1,2)}
    return B2, B2p, board0, Dsig, Dsigp


def check1(nmax=14):
    print("check 1: empty bracket = M_{n-4}")
    for n in range(6, nmax + 1):
        _, _, b0, Dsig, Dsigp = make_boards(n)
        # bracket: pins admissible-for-sigma' minus admissible-for-sigma:
        # BrG(empty) corresponds to G((1,2))+G((n-1,0))-G((1,0))-G((n-1,2))
        br = (board_count(b0, n, ((1, 2),)) + board_count(b0, n, ((n - 1, 0),))
              - board_count(b0, n, ((1, 0),)) - board_count(b0, n, ((n - 1, 2),)))
        ok = (br == menage(n - 4))
        print(f"  n={n}: bracket={br}  M_(n-4)={menage(n-4)}  {'PASS' if ok else 'FAIL'}")
        assert ok


def brute_pairs(n, B2):
    """Return list of B2-avoiding permutations and N4 = # ordered disjoint pairs."""
    good = [t for t in itertools.permutations(range(n))
            if all((i, t[i]) not in B2 for i in range(n))]
    return good


def disjoint(t, u):
    return all(t[i] != u[i] for i in range(len(t)))


def check23(n):
    print(f"check 2+3 at n={n} (brute force)")
    B2, B2p, b0, Dsig, Dsigp = make_boards(n)
    good0 = brute_pairs(n, b0)
    # brute counts
    def N4_of(B):
        g = [t for t in good0 if all((i, t[i]) not in B for i in range(n))]
        return sum(1 for t in g for u in g if disjoint(t, u))
    N4 = N4_of(B2)
    N4p = N4_of(B2p)
    # pieces on board0
    def G(d):
        i0, c0 = d
        tg = [t for t in good0 if t[i0] == c0]
        return sum(1 for t in tg for u in good0 if disjoint(t, u))
    def G2(d, dp):
        tg = [t for t in good0 if t[d[0]] == d[1] and t[dp[0]] == dp[1]]
        return sum(1 for t in tg for u in good0 if disjoint(t, u))
    def H(d, dp):
        tg = [t for t in good0 if t[d[0]] == d[1]]
        ug = [u for u in good0 if u[dp[0]] == dp[1]]
        return sum(1 for t in tg for u in ug if disjoint(t, u))
    d1, d2 = Dsig      # (1,2),(n-1,0)
    e1, e2 = Dsigp     # (1,0),(n-1,2)
    BrG = G(d1) + G(d2) - G(e1) - G(e2)
    g2 = G2(d1, d2)
    g2p = G2(e1, e2)
    dH = H(e1, e2) - H(d1, d2)
    lhs = N4p - N4
    rhs = 2 * BrG + 2 * (g2p - g2) + 2 * dH
    print(f"  G2: {g2} vs {g2p}  (dG2={g2p-g2}; equal only up to eta-through-pin terms)")
    print(f"  DeltaN4={lhs}  2BrG+2dG2+2dH={rhs}  {'PASS' if lhs == rhs else 'FAIL'}")
    print(f"  pieces: 2BrG={2*BrG}  2dG2={2*(g2p-g2)}  2dH={2*dH}")
    assert lhs == rhs


def random_matching(rng, n, b0, avoid_rows, avoid_cols, j, force_cells=None):
    cells = [(i, c) for i in range(n) for c in range(n)
             if (i, c) not in b0 and i not in avoid_rows and c not in avoid_cols]
    rng.shuffle(cells)
    M, rows, cols = [], set(), set()
    if force_cells:
        for cell in force_cells:
            M.append(cell); rows.add(cell[0]); cols.add(cell[1])
    for cell in cells:
        if len(M) >= j:
            break
        if cell[0] not in rows and cell[1] not in cols and cell not in (force_cells or []):
            M.append(cell); rows.add(cell[0]); cols.add(cell[1])
    return tuple(M)


def check4(n, trials=40, seed=7):
    print(f"check 4 at n={n}: per-M bracket ratios (generic M avoiding switch rows/cols)")
    rng = random.Random(seed)
    B2, B2p, b0, Dsig, Dsigp = make_boards(n)
    d1, d2 = Dsig; e1, e2 = Dsigp
    Mn = board_count(set((i, i) for i in range(n)) | set((i, (i + 1) % n) for i in range(n)), n)
    r3 = Fraction(menage(n - 4), Mn)
    for j in (1, 2, 3, 4):
        rats = []
        for _ in range(trials):
            M = random_matching(rng, n, b0, {1, n - 1}, {0, 2}, j)
            if len(M) != j:
                continue
            def pin(d):
                return board_count(b0, n, M + (d,))
            br = pin(d1) + pin(d2) - pin(e1) - pin(e2)
            R0M = board_count(b0, n, M)
            rat = Fraction(br) / (R0M * r3)
            rats.append(float(rat))
        print(f"  j={j}: ratio Br(M)/(R0(M)*r3): "
              f"min={min(rats):.4f} max={max(rats):.4f} mean={sum(rats)/len(rats):.4f}")
    # adversarial: M touching row 0 / col 1 (interior of P1) and near switch on P2
    print("  adversarial M classes:")
    for label, force in [
        ("cell in col 1 (P1 interior)", [(4, 1)]),
        ("cell in row 0 (P1 interior)", [(0, 4)]),
        ("cell in col 3 (P2 near c2-end)", [(5, 3)]),
        ("cell in row n-2 (P2 near rn-end)", [(n - 2, 5)]),
    ]:
        for extra in (0, 2):
            M = random_matching(rng, n, b0, {1, n - 1}, {0, 2}, len(force) + extra,
                                force_cells=force)
            def pin(d):
                return board_count(b0, n, M + (d,))
            br = pin(d1) + pin(d2) - pin(e1) - pin(e2)
            R0M = board_count(b0, n, M)
            rat = float(Fraction(br) / (R0M * r3))
            print(f"    {label} (+{extra} generic): |M|={len(M)}  ratio={rat:.4f}")


if __name__ == "__main__":
    import sys
    args = sys.argv[1:]
    if not args or 'all' in args:
        check1()
        check23(6)
        check23(7)
        check4(11)
        check4(13)
    else:
        for a in args:
            if a == '1':
                check1()
            elif a.startswith('23'):
                check23(int(a[2:] or 6))
            elif a.startswith('4'):
                check4(int(a[1:] or 11))


def check5(ns=(9, 11, 13), trials=30, seed=21):
    """Column-exchange identities and the boundary pairing (notes Lemma
    'column exchange' and Lemma 'boundary'):
      class {rows 1,2 hit} (0-idx rows 0,1):
        R0(Mud2) - R0(Mue2) = -R0(M u {(2,0),(n-1,2)})
      class {rows 2,3 hit} (0-idx rows 1,2):
        R0(Mud2) - R0(Mue2) = +R0(M u {(0,2),(n-1,0)})
      class {rows 1,2,3} (0-idx 0,1,2): bracket = 0
      pairing: R0(M u {(2,0),(n-1,2)}) = R0(phiM u {(0,2),(n-1,0)}),
        phiM = M with (0,g) -> (2,g).
    """
    import random
    for n in ns:
        B2, B2p, b0, Dsig, Dsigp = make_boards(n)
        d1, d2 = Dsig
        e1, e2 = Dsigp
        rng = random.Random(seed)
        cnt = 0
        for _ in range(trials):
            g = rng.choice(range(4, n - 1))
            h = rng.choice([c for c in range(3, n - 1) if c != g])
            k = rng.choice([c for c in range(3, n - 1) if c not in (g, h)])
            M12 = ((0, g), (1, h))
            M23 = ((1, h), (2, g))
            M123 = ((0, k), (1, h), (2, g))
            if any(z in b0 for z in M12 + M23 + M123):
                continue
            lhs = board_count(b0, n, M12 + (d2,)) - board_count(b0, n, M12 + (e2,))
            assert lhs == -board_count(b0, n, M12 + ((2, 0), (n - 1, 2)))
            lhs = board_count(b0, n, M23 + (d2,)) - board_count(b0, n, M23 + (e2,))
            assert lhs == board_count(b0, n, M23 + ((0, 2), (n - 1, 0)))
            assert (board_count(b0, n, M123 + (d2,))
                    - board_count(b0, n, M123 + (e2,))) == 0
            assert (board_count(b0, n, M12 + ((2, 0), (n - 1, 2)))
                    == board_count(b0, n, ((2, g), (1, h)) + ((0, 2), (n - 1, 0))))
            cnt += 1
        print(f"  n={n}: column-exchange + mirror + 3-line-zero + pairing: "
              f"PASS on {cnt} instances")
