"""Structural analysis of the CGW completion counts C(lambda).

Questions:
  (1) How far is C(lambda) from constant?  (max/min - 1, coefficient of variation,
      dTV(P_n, Q_n) -- reproduce CGW's numbers as a check.)
  (2) Is log C(lambda) additive over cycles:  log C = c0 + sum_l h(l)*lambda_l,
      possibly plus a sign term eps*sgn(lambda)?
  (3) How do the fitted h(l) scale with n?
"""
from mpmath import mp, mpf, log, fabs
import itertools
from cgw_data import C, gamma, sgn
import numpy as np

mp.dps = 60


def logC(n):
    return {lam: log(mpf(C[n][lam])) for lam in C[n]}


def counts_vec(lam, n):
    """lambda_l counts for l = 2..n."""
    from collections import Counter
    c = Counter(lam)
    return {l: c.get(l, 0) for l in range(2, n + 1)}


def dTV_and_spread(n):
    tot = sum(gamma(l, n) * C[n][l] for l in C[n])
    totg = sum(gamma(l, n) for l in C[n])  # = number of derangements D_n
    dtv = mpf(0)
    for lam in C[n]:
        P = mpf(gamma(lam, n) * C[n][lam]) / mpf(tot)
        Q = mpf(gamma(lam, n)) / mpf(totg)
        dtv += fabs(P - Q)
    dtv /= 2
    vals = [mpf(C[n][l]) for l in C[n]]
    mx, mn = max(vals), min(vals)
    # gamma-weighted mean and sd of C
    mean = mpf(tot) / mpf(totg)
    var = sum(mpf(gamma(l, n)) * (mpf(C[n][l]) - mean) ** 2 for l in C[n]) / mpf(totg)
    return dtv, mx / mn - 1, (var ** mpf(0.5)) / mean


def fit(n, use_sign=True, gauge_part=None):
    """Least-squares fit of log C(lam) = c0 + sum_l h(l) lam_l (+ eps*sgn).

    Gauge: h(l) -> h(l) + t*l changes predictions by t*n (absorbed in c0),
    so one h is redundant; we fix h(gauge_part) = 0 (default: largest part = n).
    Returns dict h, c0, eps, residuals per partition."""
    if gauge_part is None:
        gauge_part = n
    lams = sorted(C[n].keys())
    parts = sorted({p for lam in lams for p in lam if p != gauge_part})
    lc = logC(n)
    # scale: subtract log C of the n-cycle partition for conditioning
    base = lc[(n,)] if (n,) in lc else list(lc.values())[0]
    rows, ys = [], []
    for lam in lams:
        cv = counts_vec(lam, n)
        row = [1.0] + [float(cv[p]) for p in parts]
        if use_sign:
            row.append(float(sgn(lam, n)))
        rows.append(row)
        ys.append(float(lc[lam] - base))
    A = np.array(rows)
    y = np.array(ys)
    sol, res, rank, sv = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ sol
    resid = {lam: float(y[i] - pred[i]) for i, lam in enumerate(lams)}
    h = {p: sol[1 + i] for i, p in enumerate(parts)}
    eps = sol[-1] if use_sign else 0.0
    dof = len(lams) - A.shape[1]
    return dict(h=h, c0=sol[0], eps=eps, resid=resid, dof=dof, rank=rank, ncols=A.shape[1])


if __name__ == "__main__":
    print("== global statistics (check against CGW Section 5) ==")
    for n in sorted(C):
        dtv, spread, cov = dTV_and_spread(n)
        print(f"n={n:2d}  dTV={float(dtv):.6g}  max/min-1={float(spread):.6g}  cv={float(cov):.6g}")

    print("\n== per-2-cycle increments along the chain (2^k, n-2k) ==")
    for n in sorted(C):
        chain = []
        k = 0
        while True:
            lam = tuple([2] * k + ([n - 2 * k] if n - 2 * k >= 2 else []))
            lam = tuple(sorted(lam))
            if sum(lam) != n or lam not in C[n]:
                break
            chain.append(lam)
            k += 1
        for a, b in zip(chain, chain[1:]):
            r = mpf(C[n][b]) / mpf(C[n][a]) - 1
            print(f"n={n:2d}  {a} -> {b}: ratio-1 = {float(r):+.6e}   *n^3 = {float(r)*n**3:+.4f}")

    print("\n== additive fits ==")
    for n in [8, 9, 10, 11]:
        for use_sign in [False, True]:
            f = fit(n, use_sign=use_sign)
            print(f"\nn={n}, sign_term={use_sign}  (dof={f['dof']}, rank={f['rank']}/{f['ncols']})")
            hs = "  ".join(f"h({p})={f['h'][p]:+.3e}" for p in sorted(f['h']))
            print("   " + hs)
            print(f"   eps={f['eps']:+.3e}  c0={f['c0']:+.3e}")
            mr = max(abs(v) for v in f['resid'].values())
            print(f"   max|resid| = {mr:.3e}")
            for lam, r in sorted(f['resid'].items(), key=lambda kv: -abs(kv[1]))[:4]:
                print(f"     resid {lam}: {r:+.3e}")
