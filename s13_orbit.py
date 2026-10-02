"""s13_orbit.py -- the orbit method for a ladder pair (session 13, sec:orbit).

For a ladder pair (x,y) of a uniform square with a true mark (p,p'):
  a column pair (q,c), q,c not in {p,p'}, gives the column cycle Z of (q,c) through x
  (rows x -> row holding L(x,q) in column c -> ...).  Turning Z is LEGAL in X iff
  Z meets {row 1, row 2} in 0 or 2 rows (type and mark preserved); it is CLEAN iff
  legal and y not in Z (then rho_{x,y} -> rho_{x,y} o (q c) and no other ladder row
  of a short sub-ladder is touched when Z is short).  (q,c) is CROSS-PATH for (x,y)
  iff q,c are separated by {p,p'} in rho_{x,y} (different arcs if crossed, different
  cycles Q_p, Q_p' if parallel); the turn toggles the bit iff cross-path.
Records, per instance:
  (i)   length distribution of Z over all (q,c); fraction legal / clean by length;
  (ii)  among clean pairs of length <= LMAX: fraction cross-path, vs. the fraction of
        all pairs that are separated by {p,p'} (the "genericity" (Pi));
  (iii) orbit test: T = greedy matching of clean pairs with length <= LMAX;
        rho(eps) = rho_0 o prod_{eps_Z=1} (q_Z c_Z); estimate P_eps[parallel] by
        NEPS samples; report the mean over instances, split by the state at eps=0,
        and a histogram of the per-instance orbit means.
usage: ./jm n nsamp thin seed | python3 s13_orbit.py n [--inst I] [--seed s] [--lmax L] [--neps E]
"""
import sys, random
from collections import defaultdict
from verify_identity import sigma_of
from pf_cl import cycles_of
from pf_master import pi_P
from s8_real_bias import read_squares


def main():
    n = int(sys.argv[1]); args = sys.argv[2:]
    def opt(name, default):
        return int(args[args.index(name) + 1]) if name in args else default
    ninst = opt("--inst", 3); seed = opt("--seed", 19); LMAX = opt("--lmax", 4); NEPS = opt("--neps", 64)
    FIXED = "--fixed" in args   # use the fixed matching of consecutive free columns (orbit-invariant)
    ybar = [0, 0.0, 0]          # instances, sum of orbit weight ybar, instances with empty T
    rng = random.Random(seed)
    nsq = 0; N = 0
    lenhist = defaultdict(lambda: [0, 0, 0])      # length -> [all, legal, clean]
    cross = defaultdict(lambda: [0, 0])           # length -> [clean count, clean & cross-path]
    sep_all = [0, 0]                              # all pairs: [count, separated]
    orb = {'cr': [0, 0.0, 0.0], 'par': [0, 0.0, 0.0]}   # state at eps=0 -> [inst, sum P[par], sum |T|]
    orbhist = defaultdict(int)                    # bin of per-instance P_eps[par] (10 bins)
    notT = 0
    for L in read_squares(n, sys.stdin.buffer):
        nsq += 1
        sig = sigma_of(L); pairs = []
        col = {L[0][c]: c for c in range(n)}
        for cyc in cycles_of(sig):
            m = len(cyc)
            if m < 4: continue
            for ji, jsym in enumerate(cyc):
                for alpha in range(2, m - 1): pairs.append((col[jsym], col[cyc[(ji + alpha) % m]]))
        if not pairs: continue
        pos = [[0] * n for _ in range(n)]
        for r in range(n):
            for q in range(n): pos[r][L[r][q]] = q
        colpos = [[0] * n for _ in range(n)]
        for r in range(n):
            for q in range(n): colpos[q][L[r][q]] = r
        for (p, pp) in rng.sample(pairs, min(ninst, len(pairs))):
            pi = pi_P(L, p, pp); inv = [0] * n
            for r in range(n): inv[pi[r]] = r
            c1 = 0; r = 0
            while True:
                r = pi[r]; c1 += 1
                if r == 0: break
            together = False; r = 0; d12 = 0
            while True:
                r = pi[r]; d12 += 1
                if r == 1: together = True; break
                if r == 0: break
            if together: K0 = min(d12, c1 - d12)
            else:
                c2 = 0; r = 1
                while True:
                    r = pi[r]; c2 += 1
                    if r == 1: break
                K0 = min(c1, c2)
            if K0 < 2: continue
            k = rng.randrange(1, K0); x, y = 0, 1
            for _ in range(k): x, y = inv[x], inv[y]
            rho = [pos[y][L[x][q]] for q in range(n)]
            # side labels for separation: side[q] in {1,2} for q in Q_p u Q_p' \ {p,p'} (arcs or cycles), 0 otherwise
            def sides(rh):
                side = [0] * n
                Qp = [p]; q = rh[p]
                while q != p: Qp.append(q); q = rh[q]
                if pp in Qp:
                    i = Qp.index(pp)
                    for q in Qp[1:i]: side[q] = 1
                    for q in Qp[i + 1:]: side[q] = 2
                    return side, True
                for q in Qp[1:]: side[q] = 1
                q = rh[pp]
                while q != pp: side[q] = 2; q = rh[q]
                return side, False
            side0, crossed0 = sides(rho)
            free = [q for q in range(n) if q not in (p, pp)]
            clean_pairs = []
            for a in range(len(free)):
                q = free[a]
                for b in range(a + 1, len(free)):
                    c = free[b]
                    # column cycle of (q,c) through x
                    Z = [x]; r = colpos[c][L[x][q]]
                    while r != x: Z.append(r); r = colpos[c][L[r][q]]
                    ln = len(Z)
                    m12 = (0 in Z) + (1 in Z)
                    legal = (m12 != 1); clean = legal and (y not in Z)
                    h = lenhist[min(ln, 12)]; h[0] += 1; h[1] += legal; h[2] += clean
                    sp = (side0[q] and side0[c] and side0[q] != side0[c])
                    sep_all[0] += 1; sep_all[1] += sp
                    if clean and ln <= LMAX:
                        cc = cross[ln]; cc[0] += 1; cc[1] += sp
                        clean_pairs.append((q, c))
            # matching: fixed (consecutive free columns) or greedy among clean short pairs
            if FIXED:
                cs = set(clean_pairs); T = []
                for a in range(0, len(free) - 1, 2):
                    q, c = free[a], free[a + 1]
                    if (q, c) in cs or (c, q) in cs: T.append((q, c))
            else:
                rng.shuffle(clean_pairs); used = set(); T = []
                for (q, c) in clean_pairs:
                    if q in used or c in used: continue
                    used.add(q); used.add(c); T.append((q, c))
            # orbit weight ybar = E_eta[fraction of T separated in rho(eta)]
            ybar[0] += 1
            if not T: ybar[2] += 1
            else:
                tot = 0
                for _ in range(NEPS):
                    rh = list(rho)
                    for (q, c) in T:
                        if rng.random() < 0.5: rh[q], rh[c] = rh[c], rh[q]
                    sd, _ = sides(rh)
                    tot += sum(1 for (q, c) in T if sd[q] and sd[c] and sd[q] != sd[c]) / len(T)
                ybar[1] += tot / NEPS
            state = 'cr' if crossed0 else 'par'
            o = orb[state]; o[0] += 1; o[2] += len(T)
            if not T:
                notT += 1; pm = 0.0 if crossed0 else 1.0
            else:
                cnt = 0
                for _ in range(NEPS):
                    rh = list(rho)
                    for (q, c) in T:
                        if rng.random() < 0.5: rh[q], rh[c] = rh[c], rh[q]   # rho o (q c)
                    q = rh[p]; par = True
                    while q != p:
                        if q == pp: par = False; break
                        q = rh[q]
                    cnt += par
                pm = cnt / NEPS
            o[1] += pm; orbhist[min(9, int(pm * 10))] += 1
            N += 1
    print(f"n={n}: squares={nsq} ladder-pair instances={N}  LMAX={LMAX} NEPS={NEPS}")
    print("  column-cycle length through x:  len   count   P[legal]  P[clean]")
    for ln in sorted(lenhist):
        a, lg, cl = lenhist[ln]
        print(f"      {ln:3d}{'+' if ln == 12 else ' '}  {a:8d}   {lg/a:.3f}    {cl/a:.3f}")
    print(f"  all pairs: P[separated by {{p,p'}}] = {sep_all[1]/max(1,sep_all[0]):.4f}")
    print("  clean pairs by length: len  count  P[cross-path | clean]")
    for ln in sorted(cross):
        a, s = cross[ln]
        print(f"      {ln:3d}  {a:7d}   {s/max(1,a):.4f}")
    print("  orbit test (random sub-matching of clean short pairs):")
    for st in ('cr', 'par'):
        o = orb[st]
        if o[0]: print(f"    state at eps=0 = {st:3s}: inst={o[0]:5d}  mean |T|={o[2]/o[0]:.2f}  mean P_eps[parallel]={o[1]/o[0]:.4f}")
    print(f"    instances with empty T: {notT}")
    print(f"  orbit weight ({'fixed' if FIXED else 'greedy'} matching): E[ybar] = {ybar[1]/max(1,ybar[0]):.4f} over {ybar[0]} instances, empty T in {ybar[2]}")
    print("    histogram of per-instance P_eps[parallel] (bins of 0.1):", [orbhist[b] for b in range(10)])


if __name__ == "__main__":
    main()
