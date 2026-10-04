"""s13_hazard.py -- what a chain-of-toggles argument for the hazard hypothesis would need, measured.

For sampled completions of a fixed rectangle (rows 0,1 fixed), marks (p,pp) of a class, first ladder pair (x1,y1),
rho = rho_{x1,y1}, nu, candidates t < nu with columns (q_t, c_t) = (rho^t p, rho^t pp) and status
  good_t  :=  x1, y1 on different (q_t,c_t)-column cycles  AND  the cycle through x1 or through y1 meets {0,1} in 0 or 2 rows.
For j = 1..M with candidates 1..j all bad, the 'toggle' for candidate j is a row switch rho(z,z') of two free rows z,z'
(not 0,1,x1,y1) on D_j := the (q_j,c_j)-cycle through x1 (which contains y1), with y1 strictly between z and z' in the
cyclic order from x1, whose row cycle through c_j avoids q_j (AM 'good pair').  After the toggle x1,y1 are on different
(q_j,c_j)-cycles (AM Lemma 2.4(1)).  Measured per (L,j):
  d        = number of toggles,
  d_ok     = number of toggles that make candidate j GOOD (separation + admissibility),
  d_pres   = number of toggles that are ok AND leave the statuses of candidates 1..j-1 unchanged (all still bad),
  flips    = for ok toggles, how many of the candidates 1..j-1 flipped to good,
  len      = length of the row cycle of (z,z') through c_j (the footprint).
Also the empirical hazard P[j good | 1..j-1 bad, nu>j] and prod(1-hazard).
Backward side of the switching identity P[E_j]/P[F_j] = E[b_j|F_j]/E[f_j|E_j]: for each instance in F_j (first good
candidate j), b_j = number of 'merges' = row switches rho(z,z',c_j) with z free on the (q_j,c_j)-cycle through x1, z' free
on the cycle through y1, footprint avoiding q_j, whose result has candidates 1..j bad; also the number of merges and of
free pairs.  --maxtog T: at most T toggles/merges per instance are tested (uniform subsample) and counts are reweighted by
(total/T); per-toggle statistics (flip distribution, bin rates, footprint means) are toggle-weighted.
usage: ./jmfix ... | python3 s13_hazard.py n [--marks K] [--M M] [--seed s] [--classes m:alpha,...] [--maxtog T]
"""
import sys, random
from collections import Counter, defaultdict
from s8_real_bias import read_squares
from pf_cl import cycles_of
from s13_amtwo import rowcycle_cols, colcycle_rows, rowswitch

SPECIAL = (0, 1)


def pos_table(L, n):
    pos = [[0] * n for _ in range(n)]
    for r in range(n):
        for q in range(n): pos[r][L[r][q]] = q
    return pos


def first_pair(L, n, pos, p, pp):
    colpp = [0] * n
    for r in range(n): colpp[L[r][pp]] = r
    pi = [colpp[L[r][p]] for r in range(n)]
    inv = [0] * n
    for r in range(n): inv[pi[r]] = r
    return inv[0], inv[1]


def candidates(L, n, pos, p, pp, x1, y1, M):
    """returns nu and the list of (q_t, c_t) for t = 1..min(nu-1, M)."""
    rho = [pos[y1][L[x1][q]] for q in range(n)]
    d = 1; q = rho[p]
    while q != p and q != pp: q = rho[q]; d += 1
    crossed = (q == pp)
    e = 1; q = rho[pp]
    if crossed:
        while q != p: q = rho[q]; e += 1
    else:
        while q != pp: q = rho[q]; e += 1
    nu = min(d, e)
    cands = []; q, c = p, pp
    for t in range(1, min(nu, M + 1)):
        q, c = rho[q], rho[c]; cands.append((q, c))
    return nu, cands


def status(L, n, x1, y1, q, c):
    """good iff x1,y1 on different (q,c)-cycles and (cycle through x1 or y1) meets SPECIAL in 0 or 2 rows."""
    Z = colcycle_rows(L, n, q, c, x1)
    if y1 in Z: return False, Z
    Zy = colcycle_rows(L, n, q, c, y1)
    adm = (len(set(Z) & set(SPECIAL)) != 1) or (len(set(Zy) & set(SPECIAL)) != 1)
    return adm, Z   # Z does not contain y1 here: bad-but-separated (inadmissible) case


def main():
    n = int(sys.argv[1]); args = sys.argv[2:]
    def opt(name, default, conv=int):
        return conv(args[args.index(name) + 1]) if name in args else default
    K = opt('--marks', 3); M = opt('--M', 6); seed = opt('--seed', 1); maxtog = opt('--maxtog', 40)
    classes = opt('--classes', None, str)
    classes = None if classes is None else set(tuple(int(t) for t in c.split(':')) for c in classes.split(','))
    rng = random.Random(seed)
    sig = None; marks = []
    st = Counter(); allbad = Counter(); pres_by_j = defaultdict(list); ok_by_j = defaultdict(list)
    flips_by_j = defaultdict(Counter); lens_by_j = defaultdict(list); dcount_by_j = defaultdict(list); binflip = defaultdict(lambda: [0, 0, 0]); bwd_by_j = defaultdict(list); bwdm_by_j = defaultdict(list); zz_by_j = defaultdict(list); fwd_all = defaultdict(list)
    nsq = 0
    for L in read_squares(n, sys.stdin.buffer):
        nsq += 1
        if sig is None:
            sig = L[1][:]
            for cyc in cycles_of(sig):
                m = len(cyc)
                if m < 4: continue
                for ji, j in enumerate(cyc):
                    for alpha in range(2, m - 1):
                        if classes is None or (m, alpha) in classes:
                            marks.append((j, cyc[(ji + alpha) % m]))
        pos = pos_table(L, n)
        for (p, pp) in rng.sample(marks, min(K, len(marks))):
            x1, y1 = first_pair(L, n, pos, p, pp)
            nu, cands = candidates(L, n, pos, p, pp, x1, y1, M)
            st['instances'] += 1
            stats = [status(L, n, x1, y1, q, c) for (q, c) in cands]
            # empirical hazard
            prefix_bad = True
            for j, (good, Z) in enumerate(stats, start=1):
                if not prefix_bad: break
                allbad[(j, 'cond')] += 1
                if good: allbad[(j, 'good')] += 1; prefix_bad = False
                else: allbad[(j, 'bad')] += 1
            # toggles for candidate j given 1..j-1 bad and j bad
            prefix_bad = True
            for j, (good, Z) in enumerate(stats, start=1):
                if not prefix_bad: break
                if good: prefix_bad = False; break
                q, c = cands[j - 1]
                if y1 not in Z:
                    st[('sepbad', j)] += 1; fwd_all[j].append(0.0); continue   # separated but inadmissible: no straddling toggle
                st[('joinbad', j)] += 1
                # D_j = cycle through x1 (contains y1); list from x1
                lst = Z; m0 = lst.index(y1)
                excl = set(SPECIAL) | {x1, y1}
                A = [r for r in lst[1:m0] if r not in excl]; B = [r for r in lst[m0 + 1:] if r not in excl]
                toggles = [(z, zp) for z in A for zp in B]
                tot = len(toggles)
                rng.shuffle(toggles); toggles = toggles[:maxtog]
                w = (tot / len(toggles)) if toggles else 0   # weight of each sampled toggle
                d = 0; d_ok = 0; d_pres = 0
                for (z, zp) in toggles:
                    cols = rowcycle_cols(L, n, z, zp, c)
                    if q in cols: continue          # not a good pair
                    d += 1
                    Lp, _ = rowswitch(L, n, z, zp, c)
                    g2, _ = status(Lp, n, x1, y1, q, c)
                    sep = y1 not in colcycle_rows(Lp, n, q, c, x1)
                    st['sep_ok'] += sep; st['sep_tried'] += 1
                    if not g2: continue
                    d_ok += 1
                    lens_by_j[j].append((len(cols), w))
                    # statuses of earlier candidates (same columns: rows x1,y1 untouched so cands unchanged)
                    nflip = 0
                    for s in range(1, j):
                        qs, cs = cands[s - 1]
                        gs, _ = status(Lp, n, x1, y1, qs, cs)
                        if gs: nflip += 1
                    flips_by_j[j][nflip] += w
                    if j > 1:
                        b = min(int(len(cols) / n * 10), 9); binflip[b][0] += nflip * w; binflip[b][1] += (j - 1) * w; binflip[b][2] += 1
                    if nflip == 0: d_pres += 1
                # scale counts back to the full toggle set
                dcount_by_j[j].append(d * w); ok_by_j[j].append(d_ok * w); pres_by_j[j].append(d_pres * w); fwd_all[j].append(d_pres * w)
            # backward count b_j for L in F_j (1..j-1 bad, j good): merging switches (z on cycle of x1, z' on cycle of y1,
            # free, footprint through c_j avoiding q_j) whose result has candidates 1..j bad
            prefix_bad = True
            for j, (good, Z) in enumerate(stats, start=1):
                if not good: continue
                q, c = cands[j - 1]
                Zy = colcycle_rows(L, n, q, c, y1)
                excl = set(SPECIAL) | {x1, y1}
                A = [r for r in Z if r not in excl]; B = [r for r in Zy if r not in excl]
                merges = [(z, zp) for z in A for zp in B]
                tot = len(merges); rng.shuffle(merges); merges = merges[:maxtog]
                w = (tot / len(merges)) if merges else 0
                b = 0; bm = 0
                for (z, zp) in merges:
                    cols = rowcycle_cols(L, n, z, zp, c)
                    if q in cols: continue
                    bm += 1
                    Lp, _ = rowswitch(L, n, z, zp, c)
                    gj, Zp = status(Lp, n, x1, y1, q, c)
                    st['merge_joined'] += (y1 in Zp); st['merge_tried'] += 1
                    if gj: continue
                    if all(not status(Lp, n, x1, y1, *cands[s - 1])[0] for s in range(1, j)): b += 1
                bwd_by_j[j].append(b * w); bwdm_by_j[j].append(bm * w); zz_by_j[j].append(tot)
                break   # only the first good candidate: L is in F_j for exactly this j
    print(f'n={n} squares={nsq} marks/square={K} M={M} seed={seed} maxtog={maxtog} instances={st["instances"]}  separation check: {st["sep_ok"]}/{st["sep_tried"]}')
    print('  empirical hazard: j, N(1..j-1 bad), P[j good | 1..j-1 bad], P[1..j all bad]')
    cum = 1.0
    for j in range(1, M + 1):
        N = allbad[(j, 'cond')]
        if N == 0: break
        pg = allbad[(j, 'good')] / N; cum *= (1 - pg)
        print(f'    j={j}: N={N} P[good|prefix bad]={pg:.3f}  P[all bad up to j]={cum:.4f}')
    print('  bad candidates split: ' + ', '.join(f'j={j}: joined {st[("joinbad",j)]} / separated-inadmissible {st[("sepbad",j)]}' for j in range(1, M+1) if st[("joinbad",j)]+st[("sepbad",j)]))
    print('  toggles for candidate j given 1..j bad: mean #good-pair toggles d, mean #ok (j becomes good), mean #preserving (1..j-1 still bad),')
    print('  fraction preserving among ok, distribution of #earlier candidates flipped, mean footprint length/n')
    for j in range(1, M + 1):
        if not dcount_by_j[j]: break
        N = len(dcount_by_j[j]); d = sum(dcount_by_j[j]) / N; dok = sum(ok_by_j[j]) / N; dp = sum(pres_by_j[j]) / N
        fl = flips_by_j[j]; tf = sum(fl.values())
        dist = {k: round(v / tf, 3) for k, v in sorted(fl.items())} if tf else {}
        ml = sum(l * ww for l, ww in lens_by_j[j]) / sum(ww for _, ww in lens_by_j[j]) / n if lens_by_j[j] else float('nan')
        print(f'    j={j}: N={N} d={d:.1f} ok={dok:.1f} pres={dp:.1f} pres/ok={dp/dok if dok else float("nan"):.3f} flips={dist} footprint/n={ml:.3f}')
    _tail(binflip, n)
    print(f'  merging switches: joined x1,y1 in {st["merge_joined"]}/{st["merge_tried"]} cases')
    print('  identity check: j, E[f_j|E_j] (all of E_j), E[b_j|F_j], mean #merging switches |F_j, mean |Z||Z\'| |F_j, ratio E[b|F]/E[f|E], N(E_j)/N(F_j)')
    for j in range(1, M + 1):
        if not fwd_all[j] or not bwd_by_j[j]: break
        fE = sum(fwd_all[j]) / len(fwd_all[j]); bF = sum(bwd_by_j[j]) / len(bwd_by_j[j]); bm = sum(bwdm_by_j[j]) / len(bwdm_by_j[j]); zz = sum(zz_by_j[j]) / len(zz_by_j[j])
        NE = allbad[(j, 'bad')]; NF = allbad[(j, 'good')]
        print(f'    j={j}: E[f|E]={fE:.1f} E[b|F]={bF:.1f} merging={bm:.1f} |Z||Z\'|={zz:.1f} b/n^2={bF/n/n:.4f} ratio={bF/fE if fE else float("nan"):.3f} N(E)/N(F)={NE/NF if NF else float("nan"):.3f}')


def _tail(binflip, n):
    print('  flip rate per earlier candidate by footprint length bin (len/n in [b/10,(b+1)/10)): bin: rate (#ok toggles)')
    print('    ' + '  '.join(f'{b/10:.1f}: {v[0]/v[1]:.3f} ({v[2]})' for b, v in sorted(binflip.items()) if v[1]))


if __name__ == '__main__':
    main()
