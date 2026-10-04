"""s13_amtwo.py -- numerical checks for notes section (q): the Allsop-Morris Theorem 3.1 machinery with two complete
rows (rows 0,1 fixed; Q empty, so E = {0,1}).  On squares read from stdin (jmfix output), for random (r, r', c) with
r, r' notin {0,1}, run the AM recursion; on squares in X_3 check:
 (a) structure facts: gamma_L = (c2,c3)-cycle through r1 avoids r2,r3; {r2,r3} is a (c2,c3)-cycle;
 (b) Claim 2': restricted cross-switches zeta(c2,c3,x_u,x_v) within an arc: Latin, rows 0,1 unchanged, structure cells
     unchanged, row set of D unchanged, list reversed between x_u,x_v, goodness of an interleaved same-arc pair toggles;
     inverse cross-switch returns L;
 (c) Claim 3'': row switch rho(x_i,x_j,c3) for a good same-arc pair: Latin, rows 0,1 and structure unchanged, D splits
     into (reference piece containing all special rows of D and exactly one of x_i,x_j) and (the other rows, no special
     row); other cycles unchanged; switching back returns L; g increases;
 (d) Claim 4': for free x on an inactive cycle, the AM move (possibly preceded by the column-cycle switch through x)
     lands in Y_2 minus X_3 with the same r_l,c_l,s_l, and the inverse recovers L;
 (e) statistics: g(L)/f.
usage: ./jmfix n nsamp thin seed < init.bin | python3 s13_amtwo.py n [--triples T] [--seed s]
"""
import sys, random
from collections import Counter
from s8_real_bias import read_squares

SPECIAL = (0, 1)


def is_latin(L, n):
    return all(sorted(row) == list(range(n)) for row in L) and all(sorted(L[r][c] for r in range(n)) == list(range(n)) for c in range(n))


def row_of(L, n, c, s):
    for r in range(n):
        if L[r][c] == s: return r
    raise ValueError


def col_of(L, n, r, s):
    return L[r].index(s)


def rowcycle_cols(L, n, r, rp, c):
    """columns of the row cycle of rows r, rp through column c (orientation: next col has L[r][.] = L[rp][cur])."""
    cols = [c]; cur = c
    while True:
        cur = col_of(L, n, r, L[rp][cur])
        if cur == c: return cols
        cols.append(cur)


def colcycle_rows(L, n, c, cp, r):
    """rows of the (c,cp)-column cycle through r, in list order: next row has L[.][c] = L[cur][cp]."""
    rows = [r]; cur = r
    while True:
        cur = row_of(L, n, c, L[cur][cp])
        if cur == r: return rows
        rows.append(cur)


def rowswitch(L, n, r, rp, c):
    cols = rowcycle_cols(L, n, r, rp, c)
    Lp = [row[:] for row in L]
    for y in cols: Lp[r][y], Lp[rp][y] = L[rp][y], L[r][y]
    return Lp, cols


def colswitch(L, n, c, cp, r):
    rows = colcycle_rows(L, n, c, cp, r)
    Lp = [row[:] for row in L]
    for x in rows: Lp[x][c], Lp[x][cp] = L[x][cp], L[x][c]
    return Lp, rows


def cross_switch(L, n, c, cp, r, rp):
    """zeta_L(c, cp, r, rp) per Allsop-Morris 2.3; returns new square or None if undefined."""
    rows = colcycle_rows(L, n, c, cp, r)
    if rp not in rows: return None
    m = rows.index(rp)              # rp = rows[m]; d_gamma = m
    if m <= 1: return None          # need 1 < d
    rcols = rowcycle_cols(L, n, r, rp, c)
    if cp not in rcols: return None  # need d_rho(r, rp, c -> cp) finite
    spp = L[rp][c]; cpp = col_of(L, n, r, spp)        # c'' with L[r,c''] = L[rp,c]
    # partial row cycle from c'' to cp along rho(r, rp, c'')
    pcols = [cpp]; cur = cpp
    while True:
        cur = col_of(L, n, r, L[rp][cur])
        if cur == cp: break
        pcols.append(cur)
        if len(pcols) > n: return None
    Lp = [row[:] for row in L]
    for x in rows[1:m]:                # rows strictly between r and rp along the partial column cycle
        Lp[x][c], Lp[x][cp] = L[x][cp], L[x][c]
    for y in pcols:
        Lp[r][y], Lp[rp][y] = L[rp][y], L[r][y]
    Lp[r][cp] = spp
    Lp[rp][c] = L[r][cp]
    return Lp


def eta_type(L, n, r, rp, c):
    s, sp = L[r][c], L[rp][c]
    cp = col_of(L, n, r, sp)
    rpp = row_of(L, n, cp, s)
    return 1 if rpp == rp else 0


def recursion(L, n, r, rp, c, E, CQ, jmax=3):
    """AM recursion; returns (status, data) with status in {'X',k} meaning L in X_k, or ('Y', ...) if in Y_jmax, or None."""
    rr = [r, rp]; cc = [c]; ss = [L[r][c]]
    cc.append(col_of(L, n, r, L[rp][c]))          # c1: L[r0,c1] = L[r1,c0]
    for j in range(1, jmax + 1):
        sj = L[rr[j]][cc[j]]; ss.append(sj)        # s_j
        if sj == ss[j - 1]:
            return ('X', j, rr, cc, ss)
        rnext = row_of(L, n, cc[j], ss[j - 1])
        cnext = col_of(L, n, rr[j], ss[j - 1])
        if rnext in E or rnext in rr: return None
        if cnext in CQ or cnext in cc: return None
        rr.append(rnext); cc.append(cnext)
    return ('Y', jmax, rr, cc, ss)


def structure_cells(rr, cc):
    cells = set()
    for l in range(3):
        cells |= {(rr[l], cc[l]), (rr[l + 1], cc[l]), (rr[l], cc[l + 1])}
    cells.add((rr[3], cc[3]))
    return cells


def cycles_c2c3(L, n, c2, c3):
    seen = set(); cyc = []
    for r in range(n):
        if r in seen: continue
        rows = colcycle_rows(L, n, c2, c3, r)
        seen |= set(rows); cyc.append(rows)
    return cyc


def active_structure(L, n, rr, cc, E):
    """returns dict: cycles (list of lists, rotated to start at reference), active flags, arcs, free rows, g."""
    c2, c3 = cc[2], cc[3]
    excluded = set(E) | set(rr)
    cyc = cycles_c2c3(L, n, c2, c3)
    acts = []
    for rows in cyc:
        S = set(rows)
        if rr[1] in S: ref = rr[1]
        elif 0 in S: ref = 0
        elif 1 in S: ref = 1
        else: continue
        k = rows.index(ref); lst = rows[k:] + rows[:k]      # list from the reference row
        # arcs: split free positions by special rows (other than reference) -- positions 1..m-1
        arcs = []; cur = []
        for pos in range(1, len(lst)):
            x = lst[pos]
            if x in SPECIAL:
                if cur: arcs.append(cur)
                cur = []
            elif x not in excluded:
                cur.append(pos)
        if cur: arcs.append(cur)
        acts.append({'list': lst, 'ref': ref, 'arcs': arcs, 'set': S})
    free = [x for x in range(n) if x not in excluded]
    Sfree = set()
    for a in acts: Sfree |= (a['set'] - excluded)
    return acts, free, len(free) - len(Sfree)


def good(L, n, xi, xj, c2, c3):
    return c2 not in rowcycle_cols(L, n, xi, xj, c3)


def main():
    n = int(sys.argv[1]); args = sys.argv[2:]
    T = int(args[args.index('--triples') + 1]) if '--triples' in args else 400
    seed = int(args[args.index('--seed') + 1]) if '--seed' in args else 1
    rng = random.Random(seed)
    E = set(SPECIAL); CQ = set()
    st = Counter(); gf = []; agg_good = Counter(); agg_N = Counter()
    for L in read_squares(n, sys.stdin.buffer):
        st['squares'] += 1
        for _ in range(T):
            r, rp = rng.sample(range(2, n), 2); c = rng.randrange(n)
            res = recursion(L, n, r, rp, c, E, CQ)
            if res is None: st['notY'] += 1; continue
            if res[0] == 'Y': st['Y3'] += 1; continue
            _, j, rr, cc, ss = res
            st[f'X{j}'] += 1
            if j != 3: continue
            c2, c3 = cc[2], cc[3]
            Qc = structure_cells(rr, cc)
            # (a) structure facts
            gam = colcycle_rows(L, n, c2, c3, rr[1])
            st['a_gamma_avoids_r2r3'] += (rr[2] not in gam and rr[3] not in gam)
            st['a_intercalate'] += (set(colcycle_rows(L, n, c2, c3, rr[2])) == {rr[2], rr[3]})
            acts, free, g = active_structure(L, n, rr, cc, E)
            f = len(free); gf.append(g / f)
            st['a_eta_type1'] += eta_type(L, n, rr[2], rr[3], c2)
            # (b) Claim 2': pick an arc with >= 4 free rows, an interleaved same-arc pair
            for a in acts:
                lst = a['list']
                for arc in a['arcs']:
                    if len(arc) < 4: continue
                    # positions in arc (increasing); choose i<j and (u,v) interleaved within arc
                    for _try in range(3):
                        i, j2 = sorted(rng.sample(arc, 2))
                        inside = [p for p in arc if i < p < j2]; outside = [p for p in arc if p < i or p > j2]
                        if not inside or not outside: continue
                        u, v = sorted((rng.choice(inside), rng.choice(outside)))
                        xi, xj, xu, xv = lst[i], lst[j2], lst[u], lst[v]
                        if good(L, n, xu, xv, c2, c3): st['b_uv_good_skip'] += 1; continue
                        st['b_tried'] += 1
                        Lp = cross_switch(L, n, c2, c3, xu, xv)
                        if Lp is None: st['b_undefined'] += 1; continue
                        ok = is_latin(Lp, n) and Lp[0] == L[0] and Lp[1] == L[1] and all(Lp[x][y] == L[x][y] for (x, y) in Qc)
                        st['b_legal'] += ok
                        newlist = colcycle_rows(Lp, n, c2, c3, a['ref'])
                        st['b_rowset_same'] += (set(newlist) == a['set'])
                        exp = lst[:u + 1] + lst[u + 1:v][::-1] + lst[v:]
                        st['b_reversed'] += (newlist == exp)
                        st['b_toggle'] += (good(Lp, n, xi, xj, c2, c3) != good(L, n, xi, xj, c2, c3))
                        Lq = cross_switch(Lp, n, c2, c3, xu, xv)
                        st['b_inverse'] += (Lq == L)
                        break
            # (b2) Claim 2' aggregate count: good same-arc pairs per arc, keyed by arc size
            for a in acts:
                lst = a['list']
                for arc in a['arcs']:
                    aa = len(arc)
                    if aa < 4: continue
                    ng = sum(1 for p in arc for q in arc if p < q and good(L, n, lst[p], lst[q], c2, c3))
                    agg_good[aa] += ng; agg_N[aa] += 1
                    st['b2_arcs_with_no_good_pair'] += (ng == 0)
            # (c) Claim 3'': good same-arc pair row switch
            for a in acts:
                lst = a['list']
                for arc in a['arcs']:
                    if len(arc) < 2: continue
                    pairs = [(p, q) for p in arc for q in arc if p < q and good(L, n, lst[p], lst[q], c2, c3)]
                    if not pairs: st['c_nogood'] += 1; continue
                    i, j2 = rng.choice(pairs); xi, xj = lst[i], lst[j2]
                    st['c_tried'] += 1
                    Lp, cols = rowswitch(L, n, xi, xj, c3)
                    ok = is_latin(Lp, n) and Lp[0] == L[0] and Lp[1] == L[1] and all(Lp[x][y] == L[x][y] for (x, y) in Qc) and c2 not in cols
                    st['c_legal'] += ok
                    refpiece = set(colcycle_rows(Lp, n, c2, c3, a['ref']))
                    specials_on_D = {x for x in a['set'] if x in SPECIAL}
                    between = set(lst[i + 1:j2])
                    other = (a['set'] - refpiece)
                    st['c_split_ok'] += (specials_on_D <= refpiece and len({xi, xj} & refpiece) == 1 and other - {xi, xj} == between - {xi, xj} and not (other & set(SPECIAL)))
                    st['c_xi_leaves'] += (xi not in refpiece)
                    # other cycles unchanged
                    oth = [set(cy) for cy in cycles_c2c3(L, n, c2, c3) if set(cy) != a['set']]
                    oth2 = [set(cy) for cy in cycles_c2c3(Lp, n, c2, c3)]
                    st['c_others_same'] += all(o in oth2 for o in oth)
                    Lq, _ = rowswitch(Lp, n, xi, xj, c3)
                    st['c_inverse'] += (Lq == L)
                    acts2, free2, g2 = active_structure(Lp, n, rr, cc, E)
                    Sfree_new = (between | {xi, xj}) - refpiece
                    st['c_g_increase_exact'] += (g2 == g + len([y for y in Sfree_new if y not in set(E) | set(rr)]))
                    res2 = recursion(Lp, n, rr[0], rr[1], cc[0], E, CQ)
                    st['c_stays_X3'] += (res2 is not None and res2[0] == 'X' and res2[1] == 3 and res2[2] == rr and res2[3] == cc)
                    break
            # (d) Claim 4': free x on an inactive cycle
            Sfree = set()
            for a in acts: Sfree |= a['set']
            cand = [x for x in free if x not in Sfree]
            if cand:
                x = rng.choice(cand); st['d_tried'] += 1
                cols = rowcycle_cols(L, n, rr[3], x, c3)
                if c2 in cols:
                    Lpp, rows = colswitch(L, n, c2, c3, x)
                    st['d_case2'] += 1
                    okc = is_latin(Lpp, n) and Lpp[0] == L[0] and Lpp[1] == L[1] and all(Lpp[a_][b_] == L[a_][b_] for (a_, b_) in Qc) and not (set(rows) & (set(SPECIAL) | {rr[1], rr[2], rr[3]}))
                    st['d_colswitch_legal'] += okc
                    res2 = recursion(Lpp, n, rr[0], rr[1], cc[0], E, CQ)
                    st['d_Lpp_in_X3'] += (res2 is not None and res2[0] == 'X' and res2[1] == 3 and res2[2] == rr and res2[3] == cc)
                    cols2 = rowcycle_cols(Lpp, n, rr[3], x, c3)
                    st['d_rho_avoids_c2'] += (c2 not in cols2)
                    base = Lpp
                else:
                    st['d_case1'] += 1; base = L
                Lp, _ = rowswitch(base, n, rr[3], x, c3)
                res3 = recursion(Lp, n, rr[0], rr[1], cc[0], E, CQ, jmax=2)
                inY2notX3 = (res3 is not None and res3[0] == 'Y' and res3[2][:4] == rr and res3[3][:4] == cc and res3[4][:3] == ss[:3] and Lp[rr[3]][cc[3]] != ss[2])
                st['d_lands_Y2_minus_X3'] += inY2notX3
                # inverse: r'' = row with Lp[r'', c3] = s2
                rpp = row_of(Lp, n, c3, ss[2]); st['d_rpp_is_x'] += (rpp == x)
                Lq, _ = rowswitch(Lp, n, rr[3], rpp, c3)
                if c2 in cols:
                    Lq, _ = colswitch(Lq, n, c2, c3, rpp)
                st['d_inverse'] += (Lq == L)
    print(f'n={n} squares={st["squares"]} triples/square={T}')
    for k in sorted(st):
        if k != 'squares': print(f'  {k}: {st[k]}')
    print('  Claim 2\' aggregate (arc size a: arcs N_a, good pairs, required N_a*a(a-1)(a-3)/(9(a-2)), ok?):')
    for aa in sorted(agg_N):
        req = agg_N[aa] * aa * (aa - 1) * (aa - 3) / (9 * (aa - 2))
        print(f'    a={aa}: N_a={agg_N[aa]} good={agg_good[aa]} required={req:.1f} {"OK" if agg_good[aa] >= req else "FAIL"}')
    if gf:
        gf.sort(); N = len(gf)
        print(f'  g/f over X3 instances: N={N} mean={sum(gf)/N:.3f} min={gf[0]:.3f} 10%={gf[N//10]:.3f} median={gf[N//2]:.3f}  (bound 1/30={1/30:.3f})')


if __name__ == '__main__':
    main()
