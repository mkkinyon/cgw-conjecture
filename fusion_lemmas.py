"""Symbolic + exact-integer verification of the path-fusion toolkit
(Lemma W) and the bracket structure (Lemma B, generic case).

Conventions: p_e = matching polynomial (generating m_k x^k) of the path
with e edges;  p_0 = 1, p_{-1} := 1, p_{-2} := 0, and for e <= -3 the
Chebyshev continuation p_e = -(-x)^{e+2} p_{-e-4}  (from d_{-m} = -(-x)^{-m} d_m,
p_e = d_{e+2}).

Lemma W(ii):  p_a p_b - p_{a-1} p_{b+1} = -(-x)^{b+2} p_{a-b-3}.
Consequences:
  W_C := p_2 p_{C-1} - p_1 p_C = -(-x)^{C+1} p_{-C} = ... = -x^3 p_{C-4}
     (check both closed forms)
  V_C := p_2 p_C - p_3 p_{C-1} = -x^4 p_{C-5}
Lemma B generic (j-independent, P1 intact, both cuts interior to P2):
  R0(M u d1) - R0(M u e1) = eval[ -x^3 p_{C-4} p_D  * middle ]
  R0(M u d2) - R0(M u e2) = eval[ -x^4 p_{C-5} p_{D-1} * middle ]
where C = edge-length of the piece of punctured-P2 containing c2,
      D = edge-length of the piece containing r_{n-1}.
We verify the two pin-difference identities as EXACT INTEGER identities
via board_count vs direct evaluation of the predicted reduced boards.
"""
import sympy as sp
import random
from math import factorial

from switch import path_mp, board_count, make_boards, menage

x = sp.symbols('x')


def p(e):
    """Path matching polynomial with Chebyshev continuation for e < -2."""
    if e >= -1:
        if e == -1:
            return sp.Integer(1)
        return sum(sp.binomial(e - k + 1, k) * x ** k for k in range(e // 2 + 2))
    if e == -2:
        return sp.Integer(0)
    return sp.expand(-(-x) ** (e + 2) * p(-e - 4))


def check_W2():
    print("Lemma W(ii): p_a p_b - p_{a-1} p_{b+1} = -(-x)^{b+2} p_{a-b-3}")
    ok = True
    for a in range(0, 12):
        for b in range(0, 12):
            lhs = sp.expand(p(a) * p(b) - p(a - 1) * p(b + 1))
            rhs = sp.expand(-(-x) ** (b + 2) * p(a - b - 3))
            if sp.expand(lhs - rhs) != 0:
                print(f"  FAIL a={a} b={b}")
                ok = False
    print("  PASS (a,b in 0..11)" if ok else "  FAILED")
    # W_C and V_C
    for C in range(1, 12):
        WC = sp.expand(p(2) * p(C - 1) - p(1) * p(C))
        assert sp.expand(WC + x ** 3 * p(C - 4)) == 0, ("W_C fail", C)
        VC = sp.expand(p(2) * p(C) - p(3) * p(C - 1))
        assert sp.expand(VC + x ** 4 * p(C - 5)) == 0, ("V_C fail", C)
    print("  W_C = -x^3 p_{C-4},  V_C = -x^4 p_{C-5}: PASS (C=1..11)")


def eval_pieces(pieces, N):
    """eval_N of a product of path polynomials given by edge lengths,
    allowing a global x-power shift: pieces = (shift, [e1,e2,...])."""
    shift, es = pieces
    poly = [1]
    from switch import poly_mul
    for e in es:
        poly = poly_mul(poly, list(path_mp(e)))
    # multiply by x^shift with sign (+x)^shift -- handle sign outside
    return sum((-1) ** (k + shift) * poly[k] * factorial(N - k - shift)
               for k in range(len(poly)) if N - k - shift >= 0)


def p2_positions(n):
    """P2 vertex order: position t=0..2n-5: even t = c_{2+t//2}, odd t = r_{2+t//2}."""
    names = []
    for t in range(2 * n - 4):
        k = 2 + t // 2
        names.append(('c', k) if t % 2 == 0 else ('r', k))
    return names


def check_pin_diffs(n, trials=200, seed=3):
    """Exact integer check of the generic pin-difference identities."""
    rng = random.Random(seed)
    B2, B2p, b0, Dsig, Dsigp = make_boards(n)
    d1, d2 = Dsig  # (1,2),(n-1,0)
    e1, e2 = Dsigp # (1,0),(n-1,2)
    pos = p2_positions(n)
    idx = {v: t for t, v in enumerate(pos)}
    npass = 0
    for _ in range(trials):
        # random generic M: j cells, rows in 3..n-3, cols in 4..n-2 (interior of P2,
        # away from P1 (rows 0,1 cols 0,1,2? P1 = rows 0,1 cols 0,1; keep rows>=3, cols>=4)
        j = rng.choice([1, 2, 3, 4])
        rows = rng.sample(range(3, n - 2), j)
        cols = rng.sample(range(4, n - 1), j)
        # build a matching avoiding b0: cell (i,c), c != i, i+1
        rng.shuffle(cols)
        M = []
        okM = True
        for i, c in zip(rows, cols):
            if c == i or c == (i + 1) % n:
                okM = False
                break
            M.append((i, c))
        if not okM:
            continue
        M = tuple(M)
        # compute punctured-P2 piece structure: deleted positions
        dead = sorted(idx[('r', i)] for i, c in M) + sorted(idx[('c', c)] for i, c in M)
        dead = sorted(dead)
        # pieces of P2 as maximal runs of alive positions
        alive_runs = []
        prev = -1
        for t in dead + [2 * n - 4]:
            run = t - prev - 1
            alive_runs.append(run)      # number of alive vertices in run
            prev = t
        # piece edge lengths = vertices-1 (>=0 vertices)
        # first run contains c2 (position 0) iff dead[0] > 0; last contains r_{n-1}
        pieces_edges = [r - 1 for r in alive_runs if r >= 1]
        Cv = alive_runs[0]          # vertices in c2-piece
        Dv = alive_runs[-1]         # vertices in r_{n-1}-piece
        C = Cv - 1                  # edge length
        D = Dv - 1
        if C < 5 or D < 3:
            continue  # generic case only here
        middle = [r - 1 for r in alive_runs[1:-1] if r >= 2] # drop isolated vertices (0 edges fine but no effect)
        middle = [e for e in middle if e >= 1]
        # identity 1: R0(M u d1) - R0(M u e1) = eval[-x^3 p_{C-4} p_D * middle * ...]
        lhs1 = board_count(b0, n, M + (d1,)) - board_count(b0, n, M + (e1,))
        # reduced board: pieces p_{C-4}, p_D, middle; P1 contributes p2? NO:
        # in the difference, P1\r1 (=p2) multiplied W_C: full poly:
        #   p2 p_{C-1} pD - p1 p_C pD = W_C pD = -x^3 p_{C-4} p_D  (P1-side absorbed)
        # eval at N = n - j - 1 (one pin) with x^3 shift:
        N1 = n - len(M) - 1
        rhs1 = -eval_pieces((3, [C - 4, D] + middle), N1)
        # identity 2: R0(M u d2) - R0(M u e2) = eval[-x^4 p_{C-5} p_{D-1} * middle]
        lhs2 = board_count(b0, n, M + (d2,)) - board_count(b0, n, M + (e2,))
        rhs2 = -eval_pieces((4, [C - 5, D - 1] + middle), N1)
        assert lhs1 == rhs1, (n, M, C, D, lhs1, rhs1)
        assert lhs2 == rhs2, (n, M, C, D, lhs2, rhs2)
        npass += 1
    print(f"  n={n}: pin-difference identities verified exactly on {npass} random generic M")


if __name__ == '__main__':
    check_W2()
    for n in (9, 11, 13):
        check_pin_diffs(n)
