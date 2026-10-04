"""s13_dyadic.py -- check of the insert/cut switching for cycle and segment lengths at the square-selected pair (q_1,c_1).

For completions of a fixed rectangle (rows 0,1), marks (p,pp), first pair (x1,y1), candidate 1 = (q1,c1), restricted to
instances with x1 ~ y1 on the (q1,c1)-cycle D (the 'joined' part of E_1):
  ell = |D|; a, a' = sizes of the two x1-y1 arcs (a + a' = ell); g = free rows on (q1,c1)-cycles avoiding rows 0,1 and D.
Insert move: z free on D, z' free on an inactive cycle C; if the row cycle of (z,z') through c1 avoids q1, row-switch; else
first column-switch C (columns q1,c1) then row-switch.  Checked per move: result Latin, rows 0,1,x1,y1 unchanged, D' = D u C,
inverse (row-switch back, then column-switch the cycle through z' if 2-step) recovers L.  Preimage count of L': number of
free pairs (z,z') on D' whose split leaves an x1y1-piece of size in [k,2k] (must be <= 6k^2; counted exactly).
Reports: distribution of ell/n, a/n, min-arc/n conditioned on joined, P[ell in [k,2k] | joined] per dyadic k against the
bound 400k/n and the permutation value ~3k^2/n^2, mean g/(n-ell) per ell-bin (Claim 3'' predicts >= 1/31).
usage: ./jmfix ... | python3 s13_dyadic.py n [--marks K] [--seed s] [--moves m]
"""
import sys, random
from collections import Counter, defaultdict
from s8_real_bias import read_squares
from pf_cl import cycles_of
from s13_amtwo import rowcycle_cols, colcycle_rows, rowswitch, colswitch, is_latin
from s13_hazard import pos_table, first_pair, candidates

SPECIAL = (0, 1)


def all_cycles(L, n, q, c):
    seen = [False] * n; cyc = []
    for r in range(n):
        if seen[r]: continue
        Z = colcycle_rows(L, n, q, c, r)
        for x in Z: seen[x] = True
        cyc.append(Z)
    return cyc


def main():
    n = int(sys.argv[1]); args = sys.argv[2:]
    def opt(name, default, conv=int):
        return conv(args[args.index(name) + 1]) if name in args else default
    K = opt('--marks', 3); seed = opt('--seed', 1); nmoves = opt('--moves', 4)
    rng = random.Random(seed)
    sig = None; marks = []
    ells = []; arcs = []; minarcs = []; gfrac = defaultdict(list); st = Counter(); nsq = 0
    for L in read_squares(n, sys.stdin.buffer):
        nsq += 1
        if sig is None:
            sig = L[1][:]
            for cyc in cycles_of(sig):
                m = len(cyc)
                if m < 4: continue
                for ji, j in enumerate(cyc):
                    for alpha in range(2, m - 1):
                        marks.append((j, cyc[(ji + alpha) % m]))
        pos = pos_table(L, n)
        for (p, pp) in rng.sample(marks, min(K, len(marks))):
            x1, y1 = first_pair(L, n, pos, p, pp)
            nu, cands = candidates(L, n, pos, p, pp, x1, y1, 1)
            if not cands: continue
            q, c = cands[0]
            D = colcycle_rows(L, n, q, c, x1)
            st['instances'] += 1
            if y1 not in D: st['separated'] += 1; continue
            st['joined'] += 1
            ell = len(D); m0 = D.index(y1); a = m0; ap = ell - m0
            ells.append(ell); arcs.append(a); minarcs.append(min(a, ap))
            excl = set(SPECIAL) | {x1, y1}
            cyc = all_cycles(L, n, q, c)
            inactive = [Z for Z in cyc if not (set(Z) & excl)]
            g = sum(len(Z) for Z in inactive)
            b = min(int(10 * ell / n), 9)
            gfrac[b].append(g / max(1, n - ell))
            # insert moves
            Dfree = [r for r in D if r not in excl]
            if not Dfree or not inactive: continue
            for _ in range(nmoves):
                z = rng.choice(Dfree); C = rng.choice(inactive); zp = rng.choice(C)
                cols = rowcycle_cols(L, n, z, zp, c)
                two = q in cols
                L1 = L
                if two:
                    L1, _ = colswitch(L, n, q, c, zp)
                    cols1 = rowcycle_cols(L1, n, z, zp, c)
                    assert q not in cols1, 'Lemma 2.3(1) failed'
                Lp, _ = rowswitch(L1, n, z, zp, c)
                st['moves'] += 1; st['two_step'] += two
                ok = is_latin(Lp, n) and all(Lp[r] == L[r] for r in excl)
                Dp = colcycle_rows(Lp, n, q, c, x1)
                ok = ok and set(Dp) == set(D) | set(C) and y1 in Dp
                # cyclic order of marked rows preserved and C inserted after z
                marked = [r for r in D if r in excl]; markedp = [r for r in Dp if r in excl]
                i0 = markedp.index(marked[0]); ok = ok and markedp[i0:] + markedp[:i0] == marked
                # inverse
                Lb, _ = rowswitch(Lp, n, z, zp, c)
                if two:
                    Lb, _ = colswitch(Lb, n, q, c, zp)
                ok = ok and Lb == L
                st['moves_ok'] += ok
                # preimage count of Lp for k = ell (pieces of size in [ell, 2 ell] containing x1,y1): count split pairs
                k = ell
                Dpf = Dp
                cnt = 0; cnt_good = 0
                m = len(Dpf)
                for i in range(m):
                    for j in range(i + 1, m):
                        # split at pair (Dpf[i], Dpf[j]) in list order: S = strictly between i and j plus Dpf[j]; P = rest
                        P = Dpf[:i + 1] + Dpf[j + 1:]
                        if x1 in P and y1 in P and k <= len(P) <= 2 * k and Dpf[i] not in excl and Dpf[j] not in excl:
                            cnt += 1
                            if q not in rowcycle_cols(Lp, n, Dpf[i], Dpf[j], c): cnt_good += 1
                st['preimage_max'] = max(st['preimage_max'], cnt)
                st['preimage_bound_ok'] += (cnt <= 3 * k * k)
                st['preimage_checked'] += 1
            # transfer moves (part (ii)): cut S from the other arc at a good same-arc pair (w,w') inside a special-row-free
            # stretch, then insert S into the x1-y1 arc at (z, z'=w')
            arc1 = D[1:m0]                       # strictly between x1 and y1 in list order
            arc2 = D[m0 + 1:]                    # strictly between y1 and x1
            stretches = []; cur = []
            for r in arc2:
                if r in excl: stretches.append(cur); cur = []
                else: cur.append(r)
            stretches.append(cur)
            B = max(stretches, key=len) if stretches else []
            z_choices = [r for r in arc1 if r not in excl]
            if len(B) >= 4 and z_choices:
                for _ in range(nmoves):
                    i, j = sorted(rng.sample(range(len(B)), 2))
                    w, wp = B[i], B[j]
                    if q in rowcycle_cols(L, n, w, wp, c): st['transfer_badpair'] += 1; continue
                    z = rng.choice(z_choices)
                    L1, _ = rowswitch(L, n, w, wp, c)            # cut
                    S = colcycle_rows(L1, n, q, c, wp)
                    ok = set(S) == set(B[i + 1:j + 1]) and not (set(S) & excl)
                    ok = ok and colcycle_rows(L1, n, q, c, x1) == [r for r in D if r not in set(S)]
                    two = q in rowcycle_cols(L1, n, z, wp, c)
                    L2 = L1
                    if two: L2, _ = colswitch(L1, n, q, c, wp)
                    Lp, _ = rowswitch(L2, n, z, wp, c)           # insert at (z, z'=w')
                    Dp = colcycle_rows(Lp, n, q, c, x1)
                    ok = ok and is_latin(Lp, n) and all(Lp[r] == L[r] for r in excl) and len(Dp) == ell and y1 in Dp
                    m0p = Dp.index(y1)
                    ok = ok and m0p == m0 + len(S)
                    marked = [r for r in D if r in excl]; markedp = [r for r in Dp if r in excl]
                    i0 = markedp.index(marked[0]); ok = ok and markedp[i0:] + markedp[:i0] == marked
                    # S spliced in after z, ending at w'
                    iz = Dp.index(z); ok = ok and set(Dp[iz + 1:iz + 1 + len(S)]) == set(S) and Dp[iz + len(S)] == wp
                    # inverse: row switch (z, w'), column switch through w' if two-step, then row switch (w, w')
                    Lb, _ = rowswitch(Lp, n, z, wp, c)
                    if two: Lb, _ = colswitch(Lb, n, q, c, wp)
                    Lb, _ = rowswitch(Lb, n, w, wp, c)
                    ok = ok and Lb == L
                    st['transfer'] += 1; st['transfer_two'] += two; st['transfer_ok'] += ok
    N = len(ells)
    print(f'n={n} squares={nsq} instances={st["instances"]} joined={st["joined"]} separated={st["separated"]}')
    print(f'  insert moves: {st["moves"]} (two-step {st["two_step"]}), all checks passed: {st["moves_ok"]}/{st["moves"]}; '
          f'preimage counts N <= 3k^2: {st["preimage_bound_ok"]}/{st["preimage_checked"]} (max N {st["preimage_max"]})')
    print(f'  transfer moves: {st["transfer"]} (two-step insert {st["transfer_two"]}; bad pairs skipped {st["transfer_badpair"]}), '
          f'all checks passed: {st["transfer_ok"]}/{st["transfer"]}')
    print('  conditioned on joined: P[ell in [k,2k]], bound 400k/n, permutation value 3k^2/n^2 (dyadic k)')
    k = 2
    while 2 * k <= n:
        pk = sum(1 for e in ells if k <= e < 2 * k) / N
        pa = sum(1 for e in minarcs if k <= e < 2 * k) / N
        print(f'    k={k}: P[ell in [k,2k))={pk:.4f}  P[min arc in [k,2k))={pa:.4f}   400k/n={400*k/n:.2f}   3k^2/n^2={3*k*k/n/n:.4f}')
        k *= 2
    print(f'  P[ell <= n/10 | joined]={sum(1 for e in ells if e <= n/10)/N:.3f}  P[min arc <= n/20 | joined]={sum(1 for e in minarcs if e <= n/20)/N:.3f}'
          f'  mean ell/n={sum(ells)/N/n:.3f}  mean min arc/n={sum(minarcs)/N/n:.3f}')
    print('  mean g/(n-ell) by ell/n bin (Claim 3\'\' lower bound 1/30):', {b/10: round(sum(v)/len(v), 3) for b, v in sorted(gfrac.items())})


if __name__ == '__main__':
    main()
