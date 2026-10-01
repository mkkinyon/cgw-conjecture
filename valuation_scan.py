"""Exhaustive x-valuation scan of the bracket polynomial
B_M = poly(M u d1) + poly(M u d2) - poly(M u e1) - poly(M u e2)
over all matchings M with |M| <= 2 compatible with all four pins
(and a sample of |M| = 3), at n = 9, 10.

Claims tested:
  (V1) val(B_M) >= 2 always;
  (V2) val(B_M) = 2 only when M makes >= 2 "window hits", where a window
       hit is a cell touching {row 2, row n-2, col 3, col n-1} (the
       end-adjacent lines of P2) or {row 0, col 1} (P1 interior)?? --
       we EMPIRICALLY DISCOVER the v=2 classes and print them.
Baseline valuation of a single poly is 0; the empty-M bracket has val 3
(= -x^3 c-form).  Note polynomials here are exact integer coeff lists.
"""
import itertools
from collections import defaultdict

from switch import path_mp, cycle_mp, poly_mul, make_boards


def poly_board(board_cells, n, M):
    """Matching polynomial (coeff list) of board minus rows/cols of M."""
    rows = set(i for i, c in M)
    cols = set(c for i, c in M)
    assert len(rows) == len(M) and len(cols) == len(M)
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
        comp = list(cycle_mp(ne)) if ne == nv else list(path_mp(ne))
        poly = poly_mul(poly, comp)
    return poly


def padd(a, b, s=1):
    out = [0] * max(len(a), len(b))
    for i, v in enumerate(a):
        out[i] += v
    for i, v in enumerate(b):
        out[i] += s * v
    return out


def valuation(p):
    for i, v in enumerate(p):
        if v != 0:
            return i
    return None  # zero polynomial


def scan(n, jmax=2):
    B2, B2p, b0, Dsig, Dsigp = make_boards(n)
    d1, d2 = Dsig
    e1, e2 = Dsigp
    badrows = {1, n - 1}
    badcols = {0, 2}
    cells = [(i, c) for i in range(n) for c in range(n)
             if (i, c) not in b0 and i not in badrows and c not in badcols]
    def bracket_poly(M):
        P = poly_board(b0, n, M + (d1,))
        P = padd(P, poly_board(b0, n, M + (d2,)))
        P = padd(P, poly_board(b0, n, M + (e1,)), -1)
        P = padd(P, poly_board(b0, n, M + (e2,)), -1)
        return P
    stats = defaultdict(int)
    v2examples = []
    minval = 99
    for j in range(0, jmax + 1):
        for Mc in itertools.combinations(cells, j):
            rows = set(x[0] for x in Mc)
            cols = set(x[1] for x in Mc)
            if len(rows) != j or len(cols) != j:
                continue
            v = valuation(bracket_poly(tuple(Mc)))
            key = ('zero' if v is None else v, j)
            stats[key] += 1
            if v is not None and v < minval:
                minval = v
            if v == 2 and len(v2examples) < 40:
                v2examples.append(Mc)
    print(f"n={n}: min valuation over all |M|<={jmax}: {minval}")
    for key in sorted(stats, key=str):
        print(f"  val={key[0]} j={key[1]}: {stats[key]} matchings")
    print("  v=2 examples:", v2examples[:12])
    # classify v=2 examples by which lines they touch
    lines = defaultdict(int)
    for Mc in v2examples:
        tags = []
        for (i, c) in Mc:
            t = []
            if i == 0: t.append('r0')
            if c == 1: t.append('c1')
            if i == 2: t.append('r2')
            if c == 3: t.append('c3')
            if i == n - 2: t.append(f'r{n-2}=rn-2')
            if c == n - 1: t.append(f'c{n-1}=cn-1')
            tags.append('+'.join(t) if t else 'generic')
        lines[tuple(sorted(tags))] += 1
    print("  v=2 classes (line-touch tags):")
    for k, v in sorted(lines.items(), key=str):
        print(f"    {k}: {v}")


if __name__ == '__main__':
    scan(9, 2)
    scan(10, 2)
