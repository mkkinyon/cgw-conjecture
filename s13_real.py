"""s13_real.py -- the transfer (II) tested on real squares (session 13).

For uniform Latin squares (jm.c stream) and sampled marks P={j,jp}
(jp = sigma^alpha(j), alpha>=2, in a sigma-cycle of length >=4):
  frame  pi0 = pi_P(L)  (permutation of rows; u=row 0, v=row 1),
  pool   = randomized greedy cell-disjoint family of intercalates in rows>=2
           meeting column j (R-side), jp (M-side) or both (conjugating).
  bit(w) = [u ~ v in  M_w pi0 R_w]  (thm:master; VERIFIED here against the
           switched square with --verify).
Measured per instance (and for a MATCHED MODEL instance: uniform derangement
frame, uniform disjoint chords with the same side counts):
  * eps_S eps_T for independent uniform S,T   (4 E[bias^2])
  * level test: eps(w) eps(w^D) with |D|=M non-conjugating coords, M=1,2,3
  * frame: #cycles, u~v, |C_u|,|C_v|
  * pool genericity: #chords with both ends in one pi0-cycle vs null,
    and the relative position of r' along the cycle of r (4 bins).
Reference: the exact spectral prediction E_M <eps,P^M eps>_{N=n},
M ~ Bin(k',1/2) (iid transpositions, so only up to O(k^2/n)).

usage: ./jm n nsamp thin seed | python3 s13_real.py n [--inst I] [--kcap K]
          [--pairs P] [--verify V] [--seed s] [--csv file]
"""
import sys, random, math
from collections import defaultdict
from verify_identity import sigma_of, delta_same_cycle
from pf_cl import cycles_of
from pf_master import pi_P, switch_ic
from s8_real_bias import read_squares, intercalates_through


def cyc_info(p):
    """cycles_of plus index maps: cycle id and position of each point."""
    cyc = cycles_of(p)
    cid = [0] * len(p); pos = [0] * len(p)
    for i, c in enumerate(cyc):
        for t, x in enumerate(c):
            cid[x] = i; pos[x] = t
    return cyc, cid, pos


def same_cycle(p, u, v):
    x = u
    while True:
        x = p[x]
        if x == v: return True
        if x == u: return False


def bit_of(pi0, R, M, w):
    """[u~v in M_w pi0 R_w], w a set of coordinate indices; R,M: lists of
    (r,rp) or None per coordinate."""
    n = len(pi0)
    Rw = list(range(n)); Mw = list(range(n))
    for i in w:
        if R[i] is not None:
            a, b = R[i]; Rw[a], Rw[b] = b, a
        if M[i] is not None:
            a, b = M[i]; Mw[a], Mw[b] = b, a
    p = [Mw[pi0[Rw[x]]] for x in range(n)]
    return same_cycle(p, 0, 1)


def build_pool(L, j, jp, kcap, rng):
    ics = intercalates_through(L, (j, jp))
    rng.shuffle(ics)
    fam, used = [], set()
    for (r, rp, cc, ccp) in ics:
        cells = {(r, cc), (r, ccp), (rp, cc), (rp, ccp)}
        if cells & used: continue
        fam.append((r, rp, cc, ccp)); used |= cells
        if len(fam) >= kcap: break
    R, M = [], []
    for (r, rp, cc, ccp) in fam:
        R.append((r, rp) if j in (cc, ccp) else None)
        M.append((r, rp) if jp in (cc, ccp) else None)
    return fam, R, M


def model_pool(n, kR, kM, kB, rng):
    """uniform disjoint chords on rows 2..n-1, per side; conjugating ones
    are shared chords on both sides."""
    def disjoint_pairs(cnt, avoid):
        rows = [r for r in range(2, n) if r not in avoid]
        rng.shuffle(rows)
        return [(rows[2 * i], rows[2 * i + 1]) for i in range(cnt)]
    both = disjoint_pairs(kB, set())
    usedB = {x for p in both for x in p}
    rr = disjoint_pairs(kR, usedB)
    mm = disjoint_pairs(kM, usedB)
    R = [p for p in both] + rr + [None] * kM
    M = [p for p in both] + [None] * kR + mm
    return R, M


def random_derangement(n, rng):
    while True:
        p = list(range(n)); rng.shuffle(p)
        if all(p[i] != i for i in range(n)): return p


def stats_of(pi0, R, M, rng, npairs, nlevel):
    """returns dict of measured quantities for one (frame, pool)."""
    k = len(R)
    nonconj = [i for i in range(k) if not (R[i] is not None and M[i] is not None)]
    out = {}
    # eps_S eps_T
    acc = 0.0
    for _ in range(npairs):
        S = {i for i in range(k) if rng.random() < 0.5}
        T = {i for i in range(k) if rng.random() < 0.5}
        eS = 1 if bit_of(pi0, R, M, S) else -1
        eT = 1 if bit_of(pi0, R, M, T) else -1
        acc += eS * eT
    out["ee"] = acc / npairs
    # level tests
    for Mlev in (1, 2, 3):
        if len(nonconj) < Mlev: out["lev%d" % Mlev] = None; continue
        acc = 0.0
        for _ in range(nlevel):
            w = {i for i in range(k) if rng.random() < 0.5}
            D = set(rng.sample(nonconj, Mlev))
            e0 = 1 if bit_of(pi0, R, M, w) else -1
            e1 = 1 if bit_of(pi0, R, M, w ^ D) else -1
            acc += e0 * e1
        out["lev%d" % Mlev] = acc / nlevel
    # frame
    cyc, cid, pos = cyc_info(pi0)
    out["ncyc"] = len(cyc); out["tog"] = int(cid[0] == cid[1])
    out["cu"] = len(cyc[cid[0]]); out["cv"] = len(cyc[cid[1]])
    # pool genericity (non-conjugating chords, one end each)
    same = 0.0; expect = 0.0; bins = [0] * 4; nch = 0
    for i in nonconj:
        ch = R[i] if R[i] is not None else M[i]
        r, rp = ch
        c = cyc[cid[r]]
        avail = len(c) - 1 - sum(1 for x in (0, 1) if cid[x] == cid[r])
        expect += avail / (len(pi0) - 3)
        nch += 1
        if cid[rp] == cid[r]:
            same += 1
            d = (pos[rp] - pos[r]) % len(c)
            bins[min(3, int(4 * d / len(c)))] += 1
    out["same"] = same; out["same_exp"] = expect; out["bins"] = bins; out["nch"] = nch
    return out


def main():
    n = int(sys.argv[1]); args = sys.argv[2:]
    def opt(name, default):
        return int(args[args.index(name) + 1]) if name in args else default
    ninst = opt("--inst", 4); kcap = opt("--kcap", 10 ** 6); npairs = opt("--pairs", 8)
    nver = opt("--verify", 0); seed = opt("--seed", 11); nlevel = opt("--nlevel", 8)
    csvf = open(args[args.index("--csv") + 1], "w") if "--csv" in args else None
    if csvf: csvf.write("kind,n,mcyc,kR,kM,kB,ee,lev1,lev2,lev3,ncyc,tog,cu,cv,same,same_exp,b0,b1,b2,b3,nch\n")
    rng = random.Random(seed)
    agg = defaultdict(lambda: defaultdict(float)); cnt = defaultdict(int)
    nsq = 0; vbad = vtot = 0
    for L in read_squares(n, sys.stdin.buffer):
        nsq += 1
        sig = sigma_of(L)
        pairs = []
        for cyc in cycles_of(sig):
            m = len(cyc)
            if m < 4: continue
            for ji, j in enumerate(cyc):
                for alpha in range(2, m - 1):
                    pairs.append((j, cyc[(ji + alpha) % m], m))
        if not pairs: continue
        for (j, jp, mcyc) in rng.sample(pairs, min(ninst, len(pairs))):
            fam, R, M = build_pool(L, j, jp, kcap, rng)
            k = len(fam)
            if k == 0: continue
            pi0 = pi_P(L, j, jp)
            if vtot < nver:
                for _ in range(10):
                    w = {i for i in range(k) if rng.random() < 0.5}
                    L2 = L
                    for i in w: L2 = switch_ic(L2, fam[i])
                    vbad += int(bit_of(pi0, R, M, w) != delta_same_cycle(L2, j, jp)); vtot += 1
            kB = sum(1 for i in range(k) if R[i] is not None and M[i] is not None)
            kR = sum(1 for i in range(k) if R[i] is not None) - kB
            kM = sum(1 for i in range(k) if M[i] is not None) - kB
            real = stats_of(pi0, R, M, rng, npairs, nlevel)
            Rm, Mm = model_pool(n, kR, kM, kB, rng)
            model = stats_of(random_derangement(n, rng), Rm, Mm, rng, npairs, nlevel)
            for kind, st in (("real", real), ("model", model)):
                if csvf:
                    csvf.write(f"{kind},{n},{mcyc},{kR},{kM},{kB},{st['ee']:.4f},"
                               + ",".join("" if st['lev%d' % m] is None else f"{st['lev%d' % m]:.4f}" for m in (1, 2, 3))
                               + f",{st['ncyc']},{st['tog']},{st['cu']},{st['cv']},{st['same']},{st['same_exp']:.3f},"
                               + ",".join(str(b) for b in st['bins']) + f",{st['nch']}\n")
                a = agg[kind]; cnt[kind] += 1
                a["k"] += k; a["kp"] += kR + kM; a["ee"] += st["ee"]; a["ee2"] += st["ee"] ** 2
                for m in (1, 2, 3):
                    if st["lev%d" % m] is not None:
                        a["lev%d" % m] += st["lev%d" % m]; a["nlev%d" % m] += 1
                a["ncyc"] += st["ncyc"]; a["tog"] += st["tog"]
                a["same"] += st["same"]; a["same_exp"] += st["same_exp"]; a["nch"] += st["nch"]
                for b in range(4): a["b%d" % b] += st["bins"][b]
    if csvf: csvf.close()
    print(f"n={n}: squares={nsq}, instances={cnt['real']}" + (f", VERIFY {vbad}/{vtot} mismatches" if vtot else ""))
    for kind in ("real", "model"):
        a = agg[kind]; c = cnt[kind]
        if c == 0: continue
        ee = a["ee"] / c; se = math.sqrt(max(a["ee2"] / c - ee ** 2, 0) / c)
        kp = a["kp"] / c
        print(f"  {kind:5s}: mean k={a['k']/c:.1f} (k'={kp:.1f})  E[eps_S eps_T]={ee:.4f}({se:.4f})  "
              f"E[bias^2]*k'={ee*kp/4:.3f}   levels M=1,2,3: "
              + ", ".join(f"{a['lev%d'%m]/a['nlev%d'%m]:.4f}" if a['nlev%d'%m] else "-" for m in (1, 2, 3)))
        print(f"         frame: E#cycles={a['ncyc']/c:.2f}  P(u~v)={a['tog']/c:.3f}   "
              f"pool: same-cycle chords {a['same']:.0f} vs null {a['same_exp']:.1f} (of {a['nch']:.0f}); "
              f"position bins {[int(a['b%d'%b]) for b in range(4)]}")


if __name__ == "__main__":
    main()
