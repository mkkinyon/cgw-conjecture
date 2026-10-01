"""Cohn-Lempel framework for the orbit-mixing lemma (session 6).

The target bit of the KPS-transfer programme (notes sec:kps item 4) is

    bit(S) = [ u ~ v in the same cycle of  pi0 * prod_{i in S} tau_i ]

with tau_i pairwise disjoint transpositions avoiding {u,v}, S a uniform
subset.  This file verifies, instance-by-instance:

  (CL)   Cohn-Lempel: for pi0 a single cycle and disjoint chords,
         #cycles(pi0 * tau_S) = nullity_F2(A_S) + 1,
         A_S = interlacement matrix of the activated chords.

  (BIT)  Bordered-rank bit formula (single cycle):
         u ~ v  in pi0*tau_S   <=>   x_S in colspace(A_S)
         where x_i = [chord i separates u from v on the circle].

  (MC)   Multi-cycle reduction: passenger construction
         gamma = (C1 s1 C2 t1 s2 C3 t2 ... s_{p-1} Cp t_{p-1}),
         bridges B = {(s_i,t_i)}: contraction of gamma*tau_B on the
         fresh points is pi0, so
         u ~ v in pi0*tau_S  <=>  xbar_{S u B} in colspace(Abar_{S u B})
         on the augmented circle with frozen bridge coordinates.

  (SPLIT) Point-splitting: ordered product pi0 * R * M where R, M are
         products of internally-disjoint chords but a point may be
         shared between one R-chord and one M-chord: split each shared
         point r into adjacent copies r',r'' on the circle (R-chord
         attaches to the LATER copy r'', M-chord to r'... both orders
         tested) and verify the disjoint-chord formula on the split
         circle computes the bit of the ordered product.

All checks are exhaustive over S in {0,1}^k for many random instances.
"""
import itertools
import random
import sys


# ---------------------------------------------------------------- utils

def cycles_of(perm):
    n = len(perm)
    seen = [False] * n
    out = []
    for s in range(n):
        if seen[s]:
            continue
        c = []
        x = s
        while not seen[x]:
            seen[x] = True
            c.append(x)
            x = perm[x]
        out.append(c)
    return out


def same_cycle(perm, u, v):
    x = perm[u]
    while x != u:
        if x == v:
            return True
        x = perm[x]
    return u == v


def right_mult_transpositions(perm, chords):
    """perm * prod tau  (apply tau first): (pi tau)(x) = pi(tau(x))."""
    p = list(perm)
    for (a, b) in chords:  # disjoint => order irrelevant
        p[a], p[b] = p[b], p[a]
    # careful: (pi*tau)(x) = pi(tau(x)); building: result[x] = perm[tau(x)]
    tau = list(range(len(perm)))
    for (a, b) in chords:
        tau[a], tau[b] = tau[b], tau[a]
    return [perm[tau[x]] for x in range(len(perm))]


def compose(p, q):
    """(p q)(x) = p(q(x))."""
    return [p[q[x]] for x in range(len(q))]


def transposition(n, a, b):
    t = list(range(n))
    t[a], t[b] = t[b], t[a]
    return t


# ------------------------------------------------------- F2 linear algebra

def f2_rank(rows, ncols):
    rows = [r for r in rows]
    rank = 0
    pivcol = 0
    m = len(rows)
    i = 0
    for col in range(ncols):
        piv = None
        for r in range(rank, m):
            if (rows[r] >> col) & 1:
                piv = r
                break
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for r in range(m):
            if r != rank and ((rows[r] >> col) & 1):
                rows[r] ^= rows[rank]
        rank += 1
    return rank


def in_colspace(mat_rows, ncols, vec):
    """mat is symmetric (rows as ints over ncols bits); vec int.
    x in colspace(A)  <=>  rank([A|x]) == rank(A)  (A symmetric)."""
    rA = f2_rank(list(mat_rows), ncols)
    aug = [mat_rows[i] | (((vec >> i) & 1) << ncols) for i in range(len(mat_rows))]
    rAx = f2_rank(aug, ncols + 1)
    return rA == rAx


# ------------------------------------------------ circle / chord machinery

def circle_positions(circle):
    return {p: i for i, p in enumerate(circle)}


def interlace(pos, c1, c2):
    """chords c1={a,b}, c2={c,d} on circle given by pos: interlace?"""
    a, b = sorted((pos[c1[0]], pos[c1[1]]))
    c, d = sorted((pos[c2[0]], pos[c2[1]]))
    return (a < c < b < d) or (c < a < d < b)


def separates(pos, chord, u, v):
    a, b = sorted((pos[chord[0]], pos[chord[1]]))
    pu, pv = pos[u], pos[v]
    inu = a < pu <= b if False else (a < pu <= b)
    # point strictly inside arc (a,b]?  Use open intervals on positions:
    inu = a < pu < b
    inv = a < pv < b
    return inu != inv


def interlacement_data(circle, chords, u, v):
    """A as list of int bitrows, x as int, on the given chord order."""
    pos = circle_positions(circle)
    k = len(chords)
    rows = [0] * k
    x = 0
    for i in range(k):
        for j in range(i + 1, k):
            if interlace(pos, chords[i], chords[j]):
                rows[i] |= (1 << j)
                rows[j] |= (1 << i)
        if separates(pos, chords[i], u, v):
            x |= (1 << i)
    return rows, x


def bit_formula(circle, chords, frozen, free, S, u, v):
    """bit for activated set = frozen + [free[i] for i in S]."""
    act = list(frozen) + [free[i] for i in S]
    rows, x = interlacement_data(circle, act, u, v)
    k = len(act)
    xv = x
    return in_colspace(rows, k, xv)


# ----------------------------------------------------------- verification

def check_CL_and_BIT(trials=200, nmax=14, kmax=6, seed=1):
    rng = random.Random(seed)
    bad = 0
    for t in range(trials):
        n = rng.randrange(6, nmax + 1)
        circle = list(range(n))
        rng.shuffle(circle)
        pi0 = [0] * n
        for i in range(n):
            pi0[circle[i]] = circle[(i + 1) % n]
        # pick u,v and disjoint chords avoiding u,v
        pts = list(range(n))
        rng.shuffle(pts)
        u, v = pts[0], pts[1]
        rest = pts[2:]
        k = rng.randrange(1, min(kmax, len(rest) // 2) + 1)
        chords = [(rest[2 * i], rest[2 * i + 1]) for i in range(k)]
        rowsA, x = interlacement_data(circle, chords, u, v)
        for Sbits in range(1 << k):
            S = [i for i in range(k) if (Sbits >> i) & 1]
            act = [chords[i] for i in S]
            p = right_mult_transpositions(pi0, act)
            # CL
            ncyc = len(cycles_of(p))
            sub = [0] * len(S)
            for a, i in enumerate(S):
                for b, j in enumerate(S):
                    if (rowsA[i] >> j) & 1:
                        sub[a] |= (1 << b)
            nullity = len(S) - f2_rank(list(sub), len(S))
            if ncyc != nullity + 1:
                bad += 1
                print("CL FAIL", n, circle, chords, S, ncyc, nullity)
            # BIT
            xs = 0
            for a, i in enumerate(S):
                if (x >> i) & 1:
                    xs |= (1 << a)
            pred = in_colspace(sub, len(S), xs)
            if pred != same_cycle(p, u, v):
                bad += 1
                print("BIT FAIL", n, circle, chords, S, u, v)
    return bad


def multi_cycle_augment(pi0_cycles, npts):
    """Return (circle, bridges, total_points).  Fresh points get labels
    npts, npts+1, ...  circle = C1 s1 C2 t1 s2 C3 t2 ... s_{p-1} Cp t_{p-1}."""
    p = len(pi0_cycles)
    circle = list(pi0_cycles[0])
    bridges = []
    nxt = npts
    for i in range(1, p):
        s, tpt = nxt, nxt + 1
        nxt += 2
        circle.append(s)
        circle.extend(pi0_cycles[i])
        circle.append(tpt)
        bridges.append((s, tpt))
    return circle, bridges, nxt


def check_MC(trials=200, nmax=13, kmax=5, seed=2):
    rng = random.Random(seed)
    bad = 0
    for t in range(trials):
        n = rng.randrange(6, nmax + 1)
        perm = list(range(n))
        rng.shuffle(perm)  # random pi0 (any cycle type)
        pi0 = perm
        cyc = cycles_of(pi0)
        pts = list(range(n))
        rng.shuffle(pts)
        u, v = pts[0], pts[1]
        rest = pts[2:]
        k = rng.randrange(1, min(kmax, len(rest) // 2) + 1)
        chords = [(rest[2 * i], rest[2 * i + 1]) for i in range(k)]
        circle, bridges, tot = multi_cycle_augment(cyc, n)
        for Sbits in range(1 << k):
            S = [i for i in range(k) if (Sbits >> i) & 1]
            act = [chords[i] for i in S]
            p = right_mult_transpositions(pi0, act)
            truth = same_cycle(p, u, v)
            pred = bit_formula(circle, chords, bridges, chords,
                               S, u, v) if False else None
            # frozen bridges + activated subset:
            actall = list(bridges) + act
            rows, x = interlacement_data(circle, actall, u, v)
            pred = in_colspace(rows, len(actall), x)
            if pred != truth:
                bad += 1
                print("MC FAIL", n, pi0, chords, S, u, v)
    return bad


def check_SPLIT(trials=300, nmax=12, seed=3):
    """Ordered product pi0 * R * M, R/M chords may share a point.

    Model: bit = [u~v in pi0 * rho * mu] where rho = prod of R-chords
    (disjoint among themselves), mu = prod of M-chords (disjoint among
    themselves), each R-chord and M-chord may share at most one point.
    Split each shared point r into r' (first copy) r'' (second copy)
    adjacent on the circle; attach the EARLIER factor (rho, applied
    first in x -> pi0(rho(mu(x))) reading?? we fix the convention by
    testing both) to one copy and the later to the other; verify against
    direct computation, for all activation subsets of R u M.
    """
    rng = random.Random(seed)
    bad = 0
    conventions = [0, 1]
    conv_fail = [0, 0]
    for t in range(trials):
        n = rng.randrange(7, nmax + 1)
        perm = list(range(n))
        rng.shuffle(perm)
        pi0 = perm
        cyc = cycles_of(pi0)
        pts = list(range(n))
        rng.shuffle(pts)
        u, v = pts[0], pts[1]
        rest = pts[2:]
        if len(rest) < 4:
            continue
        # build R-chords and M-chords with one shared point
        r1 = (rest[0], rest[1])
        m1 = (rest[1], rest[2])  # shares rest[1] with r1
        extra = []
        if len(rest) >= 5:
            extra = [(rest[3], rest[4])]  # an M-chord disjoint from all
        Rch = [r1]
        Mch = [m1] + extra
        for Sb in range(1 << (len(Rch) + len(Mch))):
            Ra = [Rch[i] for i in range(len(Rch)) if (Sb >> i) & 1]
            Ma = [Mch[i] for i in range(len(Mch))
                  if (Sb >> (len(Rch) + i)) & 1]
            # direct: pi = pi0 * rho * mu   (apply mu, then rho, then pi0)
            nn = len(pi0)
            rho = list(range(nn))
            for (a, b) in Ra:
                rho[a], rho[b] = rho[b], rho[a]
            mu = list(range(nn))
            for (a, b) in Ma:
                mu[a], mu[b] = mu[b], mu[a]
            p = compose(pi0, compose(rho, mu))
            truth = same_cycle(p, u, v)
            for conv in conventions:
                # split circle: every point r gets copies (r,0),(r,1)
                # adjacent in pi0-cycle order; chord from factor rho
                # attaches to copy conv, chords from mu to copy 1-conv.
                circle = []
                for c in cyc:
                    for x in c:
                        circle.append((x, 0))
                        circle.append((x, 1))
                # rebuild as multi-cycle augment over split cycles
                scyc = []
                for c in cyc:
                    sc = []
                    for x in c:
                        sc.append((x, 0))
                        sc.append((x, 1))
                    scyc.append(sc)
                # relabel to ints
                lab = {}
                for sc in scyc:
                    for q in sc:
                        lab[q] = len(lab)
                icyc = [[lab[q] for q in sc] for sc in scyc]
                circle2, bridges, tot = multi_cycle_augment(
                    icyc, len(lab))
                chR = [(lab[(a, conv)], lab[(b, conv)]) for (a, b) in Ra]
                chM = [(lab[(a, 1 - conv)], lab[(b, 1 - conv)])
                       for (a, b) in Ma]
                actall = list(bridges) + chR + chM
                rows, x = interlacement_data(circle2, actall,
                                             lab[(u, 0)], lab[(v, 0)])
                pred = in_colspace(rows, len(actall), x)
                if pred != truth:
                    conv_fail[conv] += 1
    print("SPLIT conv failures (conv0, conv1):", conv_fail)
    return conv_fail


def check_SPLIT2(trials=150, nmax=12, seed=4):
    """Random two-sided configurations with paired (conjugation)
    coordinates.  Product pi0 * rho * mu; coordinates:
      - R-coords: one chord in rho
      - M-coords: one chord in mu
      - C-coords (conjugation/both-column intercalates): SAME chord in
        rho and in mu (tau pi tau), activated together.
    R-chords pairwise disjoint; M-chords pairwise disjoint; sharing of
    single points across sides allowed and handled by point-splitting
    with the verified convention (earlier-applied factor -> copy 0).
    Exhaustive over all activation vectors."""
    rng = random.Random(seed)
    bad = 0
    tested = 0
    for t in range(trials):
        n = rng.randrange(8, nmax + 1)
        perm = list(range(n))
        rng.shuffle(perm)
        pi0 = perm
        cyc = cycles_of(pi0)
        pts = list(range(n))
        rng.shuffle(pts)
        u, v = pts[0], pts[1]
        rest = pts[2:]
        # random R-side disjoint chords, M-side disjoint chords,
        # C-side (both) chords disjoint from both sides.
        rng.shuffle(rest)
        nC = rng.randrange(0, 2)
        Cch = []
        idx = 0
        for _ in range(nC):
            if idx + 1 < len(rest):
                Cch.append((rest[idx], rest[idx + 1]))
                idx += 2
        pool = rest[idx:]
        # R-chords from pool
        nR = rng.randrange(0, 3)
        Rch = []
        rpool = list(pool)
        rng.shuffle(rpool)
        i = 0
        for _ in range(nR):
            if i + 1 < len(rpool):
                Rch.append((rpool[i], rpool[i + 1]))
                i += 2
        # M-chords from pool (can reuse points used by R => sharing)
        mpool = list(pool)
        rng.shuffle(mpool)
        nM = rng.randrange(0, 3)
        Mch = []
        used_m = set()
        j = 0
        while len(Mch) < nM and j + 1 < len(mpool):
            a, b = mpool[j], mpool[j + 1]
            j += 2
            if a in used_m or b in used_m:
                continue
            Mch.append((a, b))
            used_m.update((a, b))
        k = len(Rch) + len(Mch) + len(Cch)
        if k == 0:
            continue
        # split-circle setup (copies (x,0),(x,1)) once per instance
        lab = {}
        scyc = []
        for c in cyc:
            sc = []
            for x in c:
                sc.append((x, 0))
                sc.append((x, 1))
            scyc.append(sc)
        for sc in scyc:
            for q in sc:
                lab[q] = len(lab)
        icyc = [[lab[q] for q in sc] for sc in scyc]
        circle2, bridges, tot = multi_cycle_augment(icyc, len(lab))
        # convention (verified): p = pi0(rho(mu(x))): mu -> copy 0,
        # rho -> copy 1.
        def mchord(ch):
            return (lab[(ch[0], 0)], lab[(ch[1], 0)])

        def rchord(ch):
            return (lab[(ch[0], 1)], lab[(ch[1], 1)])

        for Sb in range(1 << k):
            Ra = [Rch[i2] for i2 in range(len(Rch)) if (Sb >> i2) & 1]
            Ma = [Mch[i2] for i2 in range(len(Mch))
                  if (Sb >> (len(Rch) + i2)) & 1]
            Ca = [Cch[i2] for i2 in range(len(Cch))
                  if (Sb >> (len(Rch) + len(Mch) + i2)) & 1]
            nn = len(pi0)
            rho = list(range(nn))
            for (a, b) in Ra + Ca:
                rho[a], rho[b] = rho[b], rho[a]
            mu = list(range(nn))
            for (a, b) in Ma + Ca:
                mu[a], mu[b] = mu[b], mu[a]
            p = compose(pi0, compose(rho, mu))
            truth = same_cycle(p, u, v)
            actall = (list(bridges) + [rchord(c) for c in Ra + Ca]
                      + [mchord(c) for c in Ma + Ca])
            rows, x = interlacement_data(circle2, actall,
                                         lab[(u, 0)], lab[(v, 0)])
            pred = in_colspace(rows, len(actall), x)
            tested += 1
            if pred != truth:
                bad += 1
                if bad < 5:
                    print("SPLIT2 FAIL", n, pi0, Rch, Mch, Cch, Sb, u, v)
    print(f"SPLIT2: {tested} cases, failures = {bad}")
    return bad


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "cl"):
        b = check_CL_and_BIT()
        print("CL+BIT: failures =", b)
    if which in ("all", "mc"):
        b = check_MC()
        print("MC: failures =", b)
    if which in ("all", "split"):
        check_SPLIT()
    if which in ("all", "split2"):
        check_SPLIT2()
