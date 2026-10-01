"""TODO 7c: real-ensemble grounding of the annealed bias at n=30..100.

Reads uniform Latin squares from the Jacobson-Matthews sampler (jm.c,
raw bytes), and for sampled merge instances (marked pair P = {j,jp} =
{j, sigma^alpha(j)}, alpha >= 2, inside a sigma-cycle of length >= 4):

  - extracts a randomized greedy maximal family of pairwise
    cell-disjoint intercalates in rows >= 2 meeting column j or jp
    (the useful coordinates: j-only -> copy-1 chord, jp-only ->
    copy-0 chord, both -> two tied lifts),
  - builds the master-formula data (doubled circle + bridges) exactly
    as in pf_master.py (END-TO-END VERIFIED at n=5,6),
  - computes the orbit bias exactly (k <= kexh) or by MC over subsets,
  - diagnostics: #separating lifted chords, #cut chords, r=0 orbits.

With --verify, additionally checks the formula's bit against the
actual switched square on random subsets (ground truth at large n).

Usage: ./jm n nsamp thin seed | python3 s7_real_bias.py n
           [--kcap K] [--inst I] [--ssamp S] [--verify V] [--seed s]
"""
import sys
import random
from collections import defaultdict

from verify_identity import sigma_of, delta_same_cycle
from pf_cl import (cycles_of, interlacement_data, in_colspace,
                   multi_cycle_augment)
from pf_master import pi_P, switch_ic


def read_squares(n, stream):
    need = n * n
    while True:
        buf = stream.read(need)
        if len(buf) < need:
            return
        yield [[buf[r * n + c] for c in range(n)] for r in range(n)]


def intercalates_through(L, cols):
    """intercalates (r, rp, c, cp) with rows >= 2 meeting a column in
    `cols`; O(n^2) per anchor column."""
    n = len(L)
    out = set()
    for c in cols:
        for cp in range(n):
            if cp == c:
                continue
            pairmap = {}
            for r in range(2, n):
                pairmap[(L[r][c], L[r][cp])] = r
            for r in range(2, n):
                a, b = L[r][c], L[r][cp]
                rp = pairmap.get((b, a))
                if rp is not None and rp > r:
                    cc, ccp = min(c, cp), max(c, cp)
                    out.add((r, rp, cc, ccp))
    return sorted(out)


def build_instance(L, j, jp, kcap, rng):
    """greedy family + master-formula data; returns dict or None."""
    ics = intercalates_through(L, (j, jp))
    rng.shuffle(ics)
    fam, used = [], set()
    for (r, rp, cc, ccp) in ics:
        cells = {(r, cc), (r, ccp), (rp, cc), (rp, ccp)}
        if cells & used:
            continue
        fam.append((r, rp, cc, ccp))
        used |= cells
        if len(fam) >= kcap:
            break
    pi0 = pi_P(L, j, jp)
    pcyc = cycles_of(pi0)
    lab = {}
    dcyc = []
    for cyc0 in pcyc:
        sc = []
        for x in cyc0:
            sc.append((x, 0))
            sc.append((x, 1))
        dcyc.append(sc)
    for sc in dcyc:
        for q in sc:
            lab[q] = len(lab)
    icyc = [[lab[q] for q in sc] for sc in dcyc]
    circle, bridges, tot = multi_cycle_augment(icyc, len(lab))
    coord_chords = []
    for (r, rp, cc, ccp) in fam:
        ch = []
        if cc == j or ccp == j:
            ch.append((lab[(r, 1)], lab[(rp, 1)]))
        if cc == jp or ccp == jp:
            ch.append((lab[(r, 0)], lab[(rp, 0)]))
        coord_chords.append(ch)
    return dict(fam=fam, circle=circle, bridges=bridges,
                chords=coord_chords, u=lab[(0, 0)], v=lab[(1, 0)],
                ncyc=len(pcyc))


def bias_of(inst, kexh, ssamp, rng):
    k = len(inst["chords"])
    circle, bridges = inst["circle"], inst["bridges"]
    u, v = inst["u"], inst["v"]
    if k <= kexh:
        hits = 0
        for Sb in range(1 << k):
            act = list(bridges)
            for i in range(k):
                if (Sb >> i) & 1:
                    act.extend(inst["chords"][i])
            rows, x = interlacement_data(circle, act, u, v)
            hits += in_colspace(rows, len(act), x)
        return hits / (1 << k) - 0.5, (1 << k)
    hits = 0
    for _ in range(ssamp):
        act = list(bridges)
        for i in range(k):
            if rng.random() < 0.5:
                act.extend(inst["chords"][i])
        rows, x = interlacement_data(circle, act, u, v)
        hits += in_colspace(rows, len(act), x)
    return hits / ssamp - 0.5, ssamp


def diagnostics(inst):
    """with ALL coordinates active: #separating chords, #cut chords."""
    act = list(inst["bridges"])
    for ch in inst["chords"]:
        act.extend(ch)
    rows, x = interlacement_data(inst["circle"], act, inst["u"], inst["v"])
    nb = len(inst["bridges"])
    nsep = ncut = 0
    for i in range(nb, len(act)):
        if (x >> i) & 1:
            nsep += 1
            if rows[i] == 0:
                ncut += 1
    return nsep, ncut


def verify_instance(L, j, jp, inst, rng, nsub=25):
    k = len(inst["chords"])
    bad = 0
    for _ in range(nsub):
        Sb = rng.randrange(1 << k)
        act = list(inst["bridges"])
        L2 = L
        for i in range(k):
            if (Sb >> i) & 1:
                act.extend(inst["chords"][i])
                L2 = switch_ic(L2, inst["fam"][i])
        rows, x = interlacement_data(inst["circle"], act,
                                     inst["u"], inst["v"])
        pred = in_colspace(rows, len(act), x)
        if pred != delta_same_cycle(L2, j, jp):
            bad += 1
    return bad


def main():
    n = int(sys.argv[1])
    args = sys.argv[2:]

    def opt(name, default):
        return int(args[args.index(name) + 1]) if name in args else default
    kcap = opt("--kcap", 18)
    ninst = opt("--inst", 2)
    ssamp = opt("--ssamp", 1500)
    nver = opt("--verify", 0)
    kexh = opt("--kexh", 14)
    seed = opt("--seed", 7)
    csvf = None
    if "--csv" in args:
        csvf = open(args[args.index("--csv") + 1], "w")
        csvf.write("n,k,bias,nsep,ncut,ncyc\n")
    rng = random.Random(seed)

    stats = defaultdict(lambda: [0, 0.0, 0.0, 0.0])  # kbin -> n,E|b|,mean,max
    ktot = 0
    r0 = 0
    cutn = 0
    nsq = ninst_done = 0
    vbad = vtot = 0
    kcounts = defaultdict(int)
    for L in read_squares(n, sys.stdin.buffer):
        nsq += 1
        sig = sigma_of(L)
        pairs = []
        for cyc in cycles_of(sig):
            m = len(cyc)
            if m < 4:
                continue
            for ji, j in enumerate(cyc):
                for alpha in range(2, m - 1):
                    pairs.append((j, cyc[(ji + alpha) % m]))
        if not pairs:
            continue
        for (j, jp) in rng.sample(pairs, min(ninst, len(pairs))):
            inst = build_instance(L, j, jp, kcap, rng)
            k = len(inst["chords"])
            kcounts[k] += 1
            if k == 0:
                continue
            if nver and ninst_done < nver:
                vbad += verify_instance(L, j, jp, inst, rng)
                vtot += 25
            b, _ = bias_of(inst, kexh, ssamp, rng)
            nsep, ncut = diagnostics(inst)
            if csvf:
                csvf.write(f"{n},{k},{b:.6f},{nsep},{ncut},"
                           f"{inst['ncyc']}\n")
            if nsep == 0:
                r0 += 1
            if ncut:
                cutn += 1
            ninst_done += 1
            ktot += k
            kb = k
            st = stats[kb]
            st[0] += 1
            st[1] += abs(b)
            st[2] += b
            st[3] = max(st[3], abs(b))
    print(f"n={n}: squares={nsq} instances={ninst_done} "
          f"mean k={ktot / max(1, ninst_done):.2f}  "
          f"r=0 (bias=+1/2) frac={r0 / max(1, ninst_done):.4f}  "
          f"cut-chord frac={cutn / max(1, ninst_done):.4f}")
    if vtot:
        print(f"  VERIFY: {vbad}/{vtot} mismatches")
    print("  k :  count  E|bias|  mean(signed)  max|b|   "
          "[iid 0.35k^-1/2]")
    for kb in sorted(stats):
        c, sa, sg, mx = stats[kb]
        print(f"  {kb:3d}: {c:6d}  {sa / c:.4f}   {sg / c:+.4f}     "
              f"{mx:.4f}   [{0.35 * kb ** -0.5:.4f}]")


if __name__ == "__main__":
    main()
