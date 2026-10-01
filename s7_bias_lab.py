"""Session-7 bias laboratory: which structural conditions cap |bias|?

Built on the verified master-formula machinery (pf_cl.py).  Setting:
u,v marked points on a circle, k disjoint chords, uniform activation
subset S; bit(S) = [x_S in col(A_S)]; bias = P(bit) - 1/2.

Experiments:

  pendant   E1: the pendant-star REFUTATION of local crossing conditions.
            r nested separating chords, each crossed by exactly m private
            non-separating pendants (and by nothing else).  Blocks are
            independent =>  P(bit) = (1 - 2^{-(m+1)})^r, so
            bias -> -1/2 for EVERY fixed m.  Hence "every separating
            chord is crossed by >= m others" does NOT imply small bias,
            for any constant m.  Verified here against the exact rank
            computation on an explicit circle realization.

  exhaust   E2: exhaustive classification.  ALL configurations with
            k <= kmax chords: circle = 2k+2 labelled points, u = 0,
            v = p (all p), chords = all perfect matchings of the rest.
            For each: exact bias + structural predicates.  Output: a
            CSV (s7_exhaust_k{k}.csv) and a summary of max|bias| per
            predicate class.  This is the ground-truth table for
            guessing the right genericity condition (TODO 7a(i)).

  toggle    E3: per-coordinate toggle probabilities in the generic
            ensemble; checks the toggle bound
            |bias| <= (1/2) min_i P(i does not toggle)
            and measures how min_i P(no toggle) scales with k.
"""
import csv
import itertools
import random
import sys

from pf_cl import (interlacement_data, in_colspace, f2_rank)


# ----------------------------------------------------------- exact bias

def bias_exact(circle, chords, u, v):
    k = len(chords)
    tot = 0
    for Sb in range(1 << k):
        act = [chords[i] for i in range(k) if (Sb >> i) & 1]
        rows, x = interlacement_data(circle, act, u, v)
        if in_colspace(rows, len(act), x):
            tot += 1
    return tot / (1 << k) - 0.5


def bit_of(circle, act, u, v):
    rows, x = interlacement_data(circle, act, u, v)
    return in_colspace(rows, len(act), x)


# ------------------------------------------------------------ E1 pendant

def pendant(rvals=(1, 2, 3, 4), mvals=(1, 2, 3)):
    print("E1: pendant-star construction "
          "(predict bias = (1-2^-(m+1))^r - 1/2)")
    for m in mvals:
        for r in rvals:
            # circle layout: block i occupies slot band around a_i
            # a_i = (2m+4)*i + m + 1 on the u-side, b_i mirrored on the
            # far side; u = 0 = slot 'gap' start, v = mid.
            span = (2 * m + 4) * r + 4
            N = 2 * span
            u, v = 0, span
            chords = []
            sep_idx = []
            for i in range(r):
                a = (2 * m + 4) * i + m + 2
                b = N - 1 - ((2 * m + 4) * i + m + 2)
                sep_idx.append(len(chords))
                chords.append((a, b))
                for p in range(1, m + 1):
                    chords.append((a - p, a + p))
            b = bias_exact(list(range(N)), chords, u, v)
            pred = (1 - 2 ** (-(m + 1))) ** r - 0.5
            flag = "OK" if abs(b - pred) < 1e-12 else "MISMATCH"
            print(f"  m={m} r={r} k={len(chords)}: bias={b:+.6f} "
                  f"predict={pred:+.6f}  {flag}")


# ---------------------------------------------------------- E2 exhaustive

def perfect_matchings(pts):
    if not pts:
        yield []
        return
    a = pts[0]
    for i in range(1, len(pts)):
        b = pts[i]
        rest = pts[1:i] + pts[i + 1:]
        for m in perfect_matchings(rest):
            yield [(a, b)] + m


def components(adj_rows, verts):
    """connected components of the graph restricted to verts
    (adj_rows: int bitmask rows over all chords)."""
    verts = list(verts)
    seen = set()
    ncomp = 0
    vs = set(verts)
    for s in verts:
        if s in seen:
            continue
        ncomp += 1
        stack = [s]
        seen.add(s)
        while stack:
            x = stack.pop()
            for y in vs:
                if y not in seen and ((adj_rows[x] >> y) & 1):
                    seen.add(y)
                    stack.append(y)
    return ncomp


def exhaust(k, csvout=True):
    N = 2 * k + 2
    circle = list(range(N))
    u = 0
    rows_out = []
    npos = 0
    for pv in range(1, N):
        v = pv
        pts = [p for p in range(N) if p not in (u, v)]
        for mat in perfect_matchings(pts):
            chords = mat
            kk = len(chords)
            rowsA, x = interlacement_data(circle, chords, u, v)
            sep = [i for i in range(kk) if (x >> i) & 1]
            r = len(sep)
            # predicates
            deg_all = [bin(rowsA[i]).count("1") for i in range(kk)]
            deg_sep = [sum(1 for j in sep if j != i and
                           ((rowsA[i] >> j) & 1)) for i in sep]
            mins_all = min((deg_all[i] for i in sep), default=-1)
            mins_sep = min(deg_sep, default=-1)
            conn_sep = (components(rowsA, sep) <= 1)
            # rank of A restricted to sep
            sub = []
            for a_, i in enumerate(sep):
                rrow = 0
                for b_, j in enumerate(sep):
                    if (rowsA[i] >> j) & 1:
                        rrow |= (1 << b_)
                sub.append(rrow)
            rk_sep = f2_rank(sub, r)
            b = bias_exact(circle, chords, u, v)
            npos += 1
            rows_out.append((pv, r, mins_all, mins_sep,
                             int(conn_sep), rk_sep, b))
    if csvout:
        with open(f"s7_exhaust_k{k}.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["vpos", "r", "min_sepcross_all",
                        "min_sepcross_sep", "conn_sep", "rank_sep",
                        "bias"])
            w.writerows(rows_out)
    # summaries
    print(f"E2: k={k}: {npos} configurations")

    def summarize(name, keyf):
        groups = {}
        for row in rows_out:
            key = keyf(row)
            g = groups.setdefault(key, [0, 0.0, 0.0, 0])
            g[0] += 1
            g[1] += abs(row[-1])
            g[2] = max(g[2], abs(row[-1]))
            g[3] += row[-1]
        print(f"  by {name}: (count, E|bias|, max|bias|, mean bias)")
        for key in sorted(groups):
            c, s, mx, sg = groups[key]
            print(f"    {key}: n={c:6d}  E|b|={s / c:.4f}  "
                  f"max|b|={mx:.4f}  mean={sg / c:+.4f}")

    summarize("r (num separating)", lambda t: t[1])
    summarize("(r, min sep-cross by ALL)", lambda t: (t[1], t[2]))
    summarize("(r, min sep-cross by SEP)", lambda t: (t[1], t[3]))
    summarize("(r, conn_sep)", lambda t: (t[1], t[4]))
    summarize("(r, rank_sep)", lambda t: (t[1], t[5]))
    summarize("(r, r - rank_sep)", lambda t: (t[1], t[1] - t[5]))


# ------------------------------------------------------------- E3 toggle

def toggle_stats(kvals=(4, 6, 8, 10, 12), trials=200, N=400, seed=11):
    rng = random.Random(seed)
    print("E3: toggle bound  |bias| <= (1/2) min_i P(no toggle)")
    print("  k   E|bias|   E[minNoTog/2]   bound viol.  corr-ish")
    for k in kvals:
        viol = 0
        sb = 0.0
        sm = 0.0
        for t in range(trials):
            pts = list(range(N))
            u, v = 0, N // 2
            pool = [p for p in pts if p not in (u, v)]
            rng.shuffle(pool)
            chords = [(pool[2 * i], pool[2 * i + 1]) for i in range(k)]
            circle = list(range(N))
            # bits for all subsets
            bits = [0] * (1 << k)
            for Sb in range(1 << k):
                act = [chords[i] for i in range(k) if (Sb >> i) & 1]
                bits[Sb] = 1 if bit_of(circle, act, u, v) else 0
            bias = sum(bits) / (1 << k) - 0.5
            # per-coordinate no-toggle probability
            minnt = 1.0
            for i in range(k):
                nt = 0
                cnt = 0
                for Sb in range(1 << k):
                    if (Sb >> i) & 1:
                        continue
                    cnt += 1
                    if bits[Sb] == bits[Sb | (1 << i)]:
                        nt += 1
                minnt = min(minnt, nt / cnt)
            if abs(bias) > minnt / 2 + 1e-12:
                viol += 1
            sb += abs(bias)
            sm += minnt / 2
        print(f"  {k:3d}  {sb / trials:.5f}   {sm / trials:.5f}"
              f"      {viol}")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "pendant"):
        pendant()
    if which in ("all", "exhaust"):
        kmax = int(sys.argv[2]) if len(sys.argv) > 2 else 4
        for k in range(2, kmax + 1):
            exhaust(k)
    if which == "exhaust1":
        exhaust(int(sys.argv[2]))
    if which in ("all", "toggle"):
        toggle_stats()
