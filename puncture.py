"""Touchard-Riordan puncture calculus.

Board B2 = I u P_sigma (forbidden cells for any further row), sigma a
derangement of [n].  As a bipartite graph on rows x columns it is a disjoint
union of cycles C_{2l}, one per l-cycle of sigma.

A *puncture set* is a partial matching M of the complement board (cells
(i,c) with c not in {i, sigma(i)}, no two sharing a row/column).  Define

    R(M; sigma) = #{ permutations tau of [n] : tau avoids B2, tau contains M }.

Removing the rows and columns of M from B2 leaves a board whose bipartite
graph is a disjoint union of paths and cycles; with mk = matching numbers of
that graph and N = n - |M|,

    R(M; sigma) = sum_k (-1)^k mk[k] (N-k)! .

This module computes R exactly for explicit sigma and M, via component
decomposition, with matching polynomials of paths/cycles by recurrence.
"""
from functools import lru_cache
from math import factorial, comb


@lru_cache(maxsize=None)
def path_mp(e):
    """Matching numbers of a path with e edges (e+1 vertices): m_k = C(e-k+1, k)."""
    if e < 0:
        return (1,)
    return tuple(comb(e - k + 1, k) for k in range(e // 2 + 2) if comb(e - k + 1, k) > 0) or (1,)


@lru_cache(maxsize=None)
def cycle_mp(m):
    """Matching numbers of a cycle with m vertices/edges (m >= 2):
    m_k = m/(m-k) C(m-k,k).  (For m=2 this is the multigraph 2-cycle; our
    boards never produce it since B2 is simple: minimal cycle length 4.)"""
    assert m >= 3
    return tuple(m * comb(m - k, k) // (m - k) for k in range(m // 2 + 1))


def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def board_components(sigma, rows_removed, cols_removed):
    """Component structure of B2 restricted to remaining rows/cols.

    Returns list of ('cycle', m) / ('path', e) with m = #edges of the cycle,
    e = #edges of the path (paths with 0 edges are dropped).

    B2's edges: for each remaining row i: (i -> col i) if col i remains,
    and (i -> col sigma(i)) if col sigma(i) remains.
    Vertices: remaining rows and remaining cols that appear in >= 1 edge...
    (isolated vertices contribute nothing to matchings).
    """
    n = len(sigma)
    RR = [i for i in range(n) if i not in rows_removed]
    edges = []
    for i in RR:
        if i not in cols_removed:
            edges.append(('r%d' % i, 'c%d' % i))
        if sigma[i] not in cols_removed:
            edges.append(('r%d' % i, 'c%d' % sigma[i]))
    # build adjacency
    from collections import defaultdict
    adj = defaultdict(list)
    for k, (u, v) in enumerate(edges):
        adj[u].append((v, k))
        adj[v].append((u, k))
    seen_e = [False] * len(edges)
    seen_v = set()
    comps = []
    for start in list(adj):
        if start in seen_v:
            continue
        # BFS collecting component
        stack = [start]
        vs = set()
        es = set()
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
            comps.append(('cycle', ne))
        elif ne == nv - 1:
            comps.append(('path', ne))
        else:
            raise RuntimeError("component not path/cycle")
    return comps


def punctured_count(sigma, M):
    """R(M; sigma): M a list of cells (i, c), pairwise row/col-disjoint,
    each avoiding B2."""
    n = len(sigma)
    rows = set(i for i, c in M)
    cols = set(c for i, c in M)
    assert len(rows) == len(M) and len(cols) == len(M)
    for i, c in M:
        assert c != i and c != sigma[i], "cell in B2"
    comps = board_components(sigma, rows, cols)
    poly = [1]
    for typ, m in comps:
        poly = poly_mul(poly, list(cycle_mp(m)) if typ == 'cycle' else list(path_mp(m)))
    N = n - len(M)
    return sum((-1) ** k * poly[k] * factorial(N - k) for k in range(len(poly)) if N - k >= 0)


def make_sigma(n, lam):
    """A derangement with cycle type lam (sum = n)."""
    sigma = [0] * n
    pos = 0
    for l in lam:
        for i in range(l):
            sigma[pos + i] = pos + (i + 1) % l
        pos += l
    return sigma


if __name__ == "__main__":
    import itertools
    # brute-force verification at n = 6, 7
    for n, lam in [(6, (6,)), (6, (2, 4)), (6, (3, 3)), (6, (2, 2, 2)),
                   (7, (7,)), (7, (2, 5)), (7, (2, 2, 3)), (7, (3, 4))]:
        sigma = make_sigma(n, lam)
        allperms = list(itertools.permutations(range(n)))
        disc = [t for t in allperms if all(t[i] != i and t[i] != sigma[i] for i in range(n))]
        # R(empty)
        assert punctured_count(sigma, []) == len(disc)
        # single punctures and pairs
        import random
        rng = random.Random(1)
        cells = [(i, c) for i in range(n) for c in range(n) if c != i and c != sigma[i]]
        for _ in range(60):
            j = rng.choice([1, 2, 3])
            M = []
            rows, cols = set(), set()
            for cell in rng.sample(cells, len(cells)):
                if cell[0] not in rows and cell[1] not in cols:
                    M.append(cell); rows.add(cell[0]); cols.add(cell[1])
                    if len(M) == j:
                        break
            truth = sum(1 for t in disc if all(t[i] == c for i, c in M))
            assert punctured_count(sigma, M) == truth, (n, lam, M, truth)
        print(f"n={n} lam={lam}: punctured counts verified (empty + 60 random punctures)")
