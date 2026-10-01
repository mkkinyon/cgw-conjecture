"""Session-8 harvest: (a) n=200 real-ensemble summary + k-binned
E[b^2]k with s7 data; (b) decoupling tables (E[k|mcyc], P(k small));
(c) psibar exponent fit from s8_psibar.log."""
import csv, math, re, sys
from collections import defaultdict


def kbin(fname, bins):
    rows = [r for r in csv.DictReader(open(fname))]
    out = []
    byk = defaultdict(list)
    for r in rows:
        if int(r['k']) > 0:
            byk[int(r['k'])].append(float(r['bias']))
    for lo, hi in bins:
        bs = [b for k in range(lo, hi + 1) for b in byk.get(k, [])]
        ks = [k for k in range(lo, hi + 1) for _ in byk.get(k, [])]
        if len(bs) < 25:
            continue
        kbar = sum(ks) / len(ks)
        eb2 = sum(b * b for b in bs) / len(bs)
        se = (sum((b * b - eb2) ** 2 for b in bs) / len(bs)) ** .5 \
            / len(bs) ** .5
        eab = sum(abs(b) for b in bs) / len(bs)
        out.append((lo, hi, kbar, len(bs), eb2 * kbar, se * kbar,
                    eab))
    return rows, out


def real_summary(fname, n):
    rows = [r for r in csv.DictReader(open(fname))]
    rows = [r for r in rows if int(r['k']) > 0]
    if not rows:
        print(f"n={n}: no rows")
        return
    N = len(rows)
    kb = sum(int(r['k']) for r in rows) / N
    eab = sum(abs(float(r['bias'])) for r in rows) / N
    sgn = sum(float(r['bias']) for r in rows) / N
    r0 = sum(1 for r in rows if int(r['nsep']) == 0) / N
    cut = sum(1 for r in rows if int(r['ncut']) > 0) / N
    eb2k = sum(float(r['bias']) ** 2 for r in rows) / N * kb
    print(f"n={n}: inst={N} kbar={kb:.1f} E|b|={eab:.3f} "
          f"signed={sgn:+.3f} P(r=0)={r0:.3f} cut={cut:.3f} "
          f"Eb2*k={eb2k:.2f}  atom*0.5={r0/2:+.3f}")


def decouple(fname):
    rows = [r for r in csv.DictReader(open(fname))]
    bym = defaultdict(list)
    for r in rows:
        bym[int(r['mcyc'])].append(int(r['k']))
    tot = [k for v in bym.values() for k in v]
    print(f"  {fname}: {len(tot)} inst, E[k]={sum(tot)/len(tot):.2f}, "
          f"P(k=0)={sum(1 for k in tot if k==0)/len(tot):.4f}")
    # correlation coefficient k vs mcyc
    ms = [m for m in bym for _ in bym[m]]
    ks = [k for m in bym for k in bym[m]]
    mm = sum(ms) / len(ms); km = sum(ks) / len(ks)
    cov = sum((a - mm) * (b - km) for a, b in zip(ms, ks)) / len(ms)
    sm = (sum((a - mm) ** 2 for a in ms) / len(ms)) ** .5
    sk = (sum((b - km) ** 2 for b in ks) / len(ks)) ** .5
    print(f"    corr(k, mcyc) = {cov/(sm*sk):+.4f} "
          f"(se ~ {1/len(ms)**.5:.4f})")
    for m in sorted(bym):
        v = bym[m]
        if len(v) < 40:
            continue
        print(f"    mcyc={m:3d} N={len(v):5d} E[k]={sum(v)/len(v):5.2f}"
              f" P(k<=4)={sum(1 for k in v if k<=4)/len(v):.4f}")


def psifit():
    pts = []
    import glob as _g
    lines = []
    for f in sorted(_g.glob('s8_psibar*.log')):
        lines += list(open(f))
    for line in lines:
        m = re.match(r"SM psibar d=\s*(\d+) .*: psi=([0-9.e-]+) ", line)
        if m:
            pts.append((int(m.group(1)), float(m.group(2))))
    pts.sort()
    print("psibar points:", [(d, f"{p:.3e}", f"{p*d:.3f}")
                             for d, p in pts])
    if len(pts) >= 2:
        (d1, p1), (d2, p2) = pts[0], pts[-1]
        print(f"  slope over [{d1},{d2}]: "
              f"{-math.log(p2/p1)/math.log(d2/d1):.3f}")


if __name__ == "__main__":
    import os
    print("== real-ensemble summaries ==")
    for f, n in [("s7_inst30.csv", 30), ("s7_inst50.csv", 50),
                 ("s7_inst100.csv", 100), ("s8_inst200.csv", 200)]:
        if os.path.exists(f) and os.path.getsize(f) > 100:
            real_summary(f, n)
    print("== k-binned Eb2*k, n=200 ==")
    if os.path.exists("s8_inst200.csv") and \
            os.path.getsize("s8_inst200.csv") > 100:
        _, out = kbin("s8_inst200.csv",
                      [(1, 8), (9, 16), (17, 25), (26, 26)])
        for lo, hi, kbar, N, e, se, eab in out:
            print(f"  k[{lo},{hi}] kbar={kbar:.1f} N={N} "
                  f"Eb2*k={e:.3f}(±{se:.3f}) E|b|={eab:.3f}")
    print("== decoupling ==")
    for f in ["s8_dec30.csv", "s8_dec50.csv"]:
        if os.path.exists(f) and os.path.getsize(f) > 100:
            decouple(f)
    print("== psibar ==")
    psifit()
