"""s13_cgw_min2.py -- CGW's upper bound with constant 2 instead of 3/2 for splits with min(alpha,beta) = 2.

Claim.  Let alpha = 2 <= beta, lambda -> mu the split of an m-cycle (m = beta + 2) into 2 + beta.  Every splitting
instance (L, {j,j'}) of CGW's procedure (cases 1-3, plus the switch edge 5) produces a square in S(mu) (for alpha = 2
case 3 lands in S(mu) in both forms of the cross-switch), so each lambda-square has exactly 2 (mu_m+1) m valid edges.
On the mu side, every edge into M at the joining pair {j,j'} is recovered from (M, j, j') by one of at most FOUR
inverse recipes (applied to M if {j,j'} is an A-pair of M, to switch(M) if it is a B-pair):
   (1) L = flip(M at {j,j'})                                   [case 1, and the duplicate case-2 edge]
   (3g) L = cross-switch(L1, {w1 j, w1 j'}) with L1 = flip(M)   [case 3, good form: nothing reversed]
   (3w) L = cross-switch(L1, {w1(j'), w1^{-1}(j')})  (j' the column with w^2 j = j' in L; the other ordering
       gives the pair of (3g) again)                           [case 3, wrong form: the long arc was reversed]
so deg_{G_S}(M) <= 4 mu_alpha mu_beta alpha beta and |S(lambda)| (mu_m+1) m / (|S(mu)| mu_2 mu_beta 2 beta) <= 2.
This script verifies the recipe on sampled lambda-squares: for every splitting instance it computes the target M and
checks that (L, instance) is recovered by the recipe from (M, j, j') -- i.e. the recipe is complete -- and that each
recipe candidate, when it is a valid preimage, is unique.  (Exact G_S degrees at n=5 are 12/18/24 = 2,3,4 x alpha beta.)
usage: ./jm n nsamp thin seed | python3 s13_cgw_min2.py n
"""
import sys
from collections import Counter
from s8_real_bias import read_squares
from s13_cgwcase3 import omega_of, cycles, ctype, col_cycle_through, turn, is_A, cross_switch
from s13_cgw_multigraph import switch


def key(L): return tuple(map(tuple, L))


def split_target(L, j, jp):
    """CGW splitting at the unordered pair {j,jp}: returns (M, case)."""
    om = omega_of(L)
    if is_A(L, j, jp): return turn(L, j, jp, col_cycle_through(L, j, jp, 1)), 1
    if is_A(L, om[j], om[jp]): return turn(L, om[j], om[jp], col_cycle_through(L, om[j], om[jp], 1)), 2
    L1 = cross_switch(L, om[j], om[jp]); return turn(L1, j, jp, col_cycle_through(L1, j, jp, 1)), 3


def recipe(M, j, jp, beta):
    """candidate preimages (L, pair) of M at the joining pair {j,jp}: list of (L, (a,b)) where the splitting at {a,b} should give M."""
    out = []
    if not is_A(M, j, jp): return out            # B-pair: edges arrive via the switch; handled by the caller on switch(M)
    L1 = turn(M, j, jp, col_cycle_through(M, j, jp, 1))       # flip = inverse of backflip
    w1 = omega_of(L1)
    out.append((L1, (j, jp)))                                   # (1): case 1 at {j,jp} (and case 2 at the predecessor pair)
    inv = [0] * len(M)
    for c in range(len(M)): inv[w1[c]] = c
    out.append((L1, (inv[j], inv[jp])))                         # (1'): case 2 at {w^-1 j, w^-1 j'} (same square, other instance)
    J, Jp = w1[j], w1[jp]
    if not is_A(L1, J, Jp): out.append((cross_switch(L1, J, Jp), (j, jp)))          # (3g)
    for (a, b) in ((j, jp), (jp, j)):
        Jp2 = w1[b]; c = b
        for _ in range(beta + 1): c = w1[c]
        J2 = c                                                  # J = w1^{beta+1}(j') = w1^{-1}(j')
        if J2 != Jp2 and not is_A(L1, J2, Jp2): out.append((cross_switch(L1, J2, Jp2), (j, jp)))   # (3w)
    return out


def main():
    n = int(sys.argv[1]); st = Counter(); shown = 0
    for L in read_squares(n, sys.stdin.buffer):
        om = omega_of(L); lam = ctype(L)
        for cyc in cycles(om):
            m = len(cyc)
            if m < 5: continue
            beta = m - 2
            for i, j in enumerate(cyc):
                jp = cyc[(i + 2) % m]                       # alpha = 2
                M, case = split_target(L, j, jp)
                mu = sorted(lam); mu.remove(m); mu = tuple(sorted(mu + [2, beta]))
                assert ctype(M) == mu, (case, ctype(M), mu)
                st[f'case {case}: lands in S(mu)'] += 1
                # recover the direct edge (L, {j,jp}) -> M from (M, j, jp)
                if case == 2: jj, jjp = om[j], om[jp]           # the edge arrives at the backflipped pair {w j, w j'}
                else: jj, jjp = j, jp
                cands = recipe(M, jj, jjp, beta)
                hits = [(Lc, pr) for (Lc, pr) in cands if key(Lc) == key(L) and set(pr) == {j, jp}]
                st[f'case {case}: direct edge recovered by recipe'] += bool(hits)
                if not hits and shown < 3:
                    shown += 1; print("NOT recovered:", lam, "case", case, "j,jp", j, jp)
                # verify every candidate that claims to be a preimage really is one (soundness), and count valid ones
                seen = set(); valid = 0
                for (Lc, (a, b)) in cands:
                    if ctype(Lc) != lam: continue
                    try: Mc, cc = split_target(Lc, a, b)
                    except Exception: continue
                    if key(Mc) == key(M) and (key(Lc), frozenset((a, b))) not in seen:
                        seen.add((key(Lc), frozenset((a, b)))); valid += 1
                st[f'valid candidates at this joining pair: {valid}'] += 1
                # the switch edge (L, switch(M)): recovered from switch(M) by applying the recipe to switch(switch(M)) = M
                Ms = switch(M, j, jp)
                assert key(switch(Ms, j, jp)) == key(M)
                # recover the switch edge (L, Ms) from Ms and its attributed pair: {j,jp} for cases 1,3 (B in Ms),
                # {w j, w jp} for case 2 (B in Ms; the switched pair {j,jp} is A in Ms)
                if case != 2:
                    assert not is_A(Ms, j, jp)
                    Mp = switch(Ms, j, jp); assert key(Mp) == key(M)
                    c2 = recipe(Mp, j, jp, beta)
                    st[f'case {case}: switch edge recovered'] += any(key(Lc) == key(L) and set(pr) == {j, jp} for Lc, pr in c2)
                else:
                    J, Jp = om[j], om[jp]
                    assert not is_A(Ms, J, Jp) and is_A(Ms, j, jp)
                    ws = omega_of(Ms); n_ = len(Ms)
                    def cyc_through(c):
                        out = [c]; x = ws[c]
                        while x != c: out.append(x); x = ws[x]
                        return out
                    CJ, CJp = cyc_through(J), cyc_through(Jp)
                    mJ = min(min(Ms[0][c], Ms[1][c]) for c in CJ); mJp = min(min(Ms[0][c], Ms[1][c]) for c in CJp)
                    inv = [0] * n_
                    for c in range(n_): inv[ws[c]] = c
                    jj = ws[J] if mJ < mJp else inv[J]        # the traded cycle is inverted in Ms
                    jjp = ws[Jp] if mJp < mJ else inv[Jp]
                    Mp = switch(Ms, jj, jjp)
                    Lc = turn(Mp, J, Jp, col_cycle_through(Mp, J, Jp, 1)) if is_A(Mp, J, Jp) else None
                    st['case 2: switch edge recovered'] += (Lc is not None and key(Lc) == key(L) and {jj, jjp} == {j, jp})
                    if not (Lc is not None and key(Lc) == key(L)) and shown < 6:
                        shown += 1; print("case-2 switch edge NOT recovered", (j, jp), (J, Jp), (jj, jjp))
    for k, v in sorted(st.items()): print(f"{v:7d}  {k}")


if __name__ == "__main__":
    main()
