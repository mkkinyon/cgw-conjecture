"""s13_spectral.py -- exact spectral measure of the same-cycle sign under the
random-transposition walk (session 13).

Claim: for eps(sigma) = 2*1{u ~ v in sigma} - 1 on S_N (u != v fixed), with
P f(sigma) = E_tau f(sigma tau) (tau uniform transposition),
    ||Pi_lambda eps||^2 = 2 c_lambda^2 (1 - r_lambda)        (lambda != (N))
where r_lambda = chi_lambda(tau)/d_lambda = 2 sum(contents)/(N(N-1)) and
c_lambda = <c, chi_lambda> = (d_lambda/N!) * d/dx|_{x=1} prod_{boxes}(x+content),
c(sigma) = number of cycles.  c_lambda != 0 only for lambda = (a,b,1^j).
Hence  psi_pi(d) = ||P^d eps||^2 = sum_lambda r^{2d} nu_lambda,
       E_pi[Phi_m] = <eps, P^m eps> = sum_lambda r^m nu_lambda.

usage: python3 s13_spectral.py brute N        -- brute-force check, N <= 7
       python3 s13_spectral.py table N dmax   -- psi_pi(d) and E_pi[Phi_m]
       python3 s13_spectral.py tail N         -- spectral mass near r=1 and r=-1
"""
import sys, itertools
from fractions import Fraction
from math import factorial

def partitions_abj(N):
    """shapes (a,b,1^j), a>=b>=1, a+b+j=N, plus the trivial (N)."""
    out = [(N,)]
    for b in range(1, N):
        for a in range(b, N - b + 1):
            j = N - a - b
            if j < 0: continue
            if j >= 1 and b < 1: continue
            lam = (a, b) + (1,) * j
            out.append(lam)
    return out

def hook_dim(lam):
    N = sum(lam); conj = [sum(1 for x in lam if x > j) for j in range(lam[0])]
    H = 1
    for i, li in enumerate(lam):
        for j in range(li):
            H *= (li - j - 1) + (conj[j] - i - 1) + 1
    return factorial(N) // H

def contents(lam):
    return [j - i for i, li in enumerate(lam) for j in range(li)]

def spectral_measure(N):
    """returns list of (lambda, r_lambda, nu_lambda) with exact Fractions."""
    res = []
    for lam in partitions_abj(N):
        if lam == (N,): continue
        cs = contents(lam); d = hook_dim(lam)
        # derivative at x=1 of prod (x + c): sum over boxes of prod over others
        # exactly one box has content -1 (box (2,1)) for these shapes
        prod_others = 1
        for c in cs:
            if c != -1: prod_others *= (1 + c)
        # if more than one -1 appears the derivative vanishes; our shapes have one
        assert cs.count(-1) == 1
        chat = Fraction(d * prod_others, factorial(N))
        r = Fraction(2 * sum(cs), N * (N - 1))
        nu = 2 * chat * chat * (1 - r)
        res.append((lam, r, nu))
    return res

def brute(N):
    import random
    perms = list(itertools.permutations(range(N)))
    def cyc(p):
        seen = [False]*N; c = 0
        for i in range(N):
            if not seen[i]:
                c += 1; x = i
                while not seen[x]: seen[x] = True; x = p[x]
        return c
    def same(p, u, v):
        x = u
        while True:
            x = p[x]
            if x == v: return True
            if x == u: return False
    u, v = 0, 1
    eps = {p: (1 if same(p, u, v) else -1) for p in perms}
    trans = [(i, j) for i in range(N) for j in range(i+1, N)]
    def mult(p, t):  # p composed with transposition t on the right: (p t)(x) = p(t(x))
        q = list(p); i, j = t; q[i], q[j] = p[j], p[i]; return tuple(q)
    # P eps
    f = dict(eps)
    out = []
    for d in range(0, 5):
        l2 = sum(f[p]**2 for p in perms) / len(perms)
        inner = sum(eps[p] * f[p] for p in perms) / len(perms)
        out.append((d, l2, inner))
        f = {p: sum(f[mult(p, t)] for t in trans) / len(trans) for p in perms}
    sm = spectral_measure(N)
    print("N=%d: sum nu = %s (should be 1)" % (N, sum(nu for _, _, nu in sm)))
    for d, l2, inner in out:
        pred_l2 = float(sum(nu * r**(2*d) for _, r, nu in sm))
        pred_in = float(sum(nu * r**d for _, r, nu in sm))
        print("  d=%d  ||P^d eps||^2 brute=%.6f exact=%.6f   <eps,P^d eps> brute=%.6f exact=%.6f"
              % (d, l2, pred_l2, inner, pred_in))

def table(N, dmax):
    sm = spectral_measure(N)
    print("N=%d, %d shapes, sum nu = %.12f" % (N, len(sm), float(sum(nu for _,_,nu in sm))))
    import math
    for m in [1,2,3,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512]:
        if m > dmax: break
        phim = sum(float(nu) * float(r)**m for _, r, nu in sm)
        print("  m=%4d  E_pi[Phi_m]=%+.6f   m*E=%.4f" % (m, phim, m*phim))

def tail(N):
    sm = spectral_measure(N)
    sm.sort(key=lambda t: -float(t[1]))
    print("top of spectrum (near r=1):")
    acc = 0
    for lam, r, nu in sm[:12]:
        acc += float(nu); print("   %-24s r=%.6f nu=%.3e  cum=%.4f" % (str(lam), float(r), float(nu), acc))
    print("bottom (near r=-1):")
    acc = 0
    for lam, r, nu in sorted(sm, key=lambda t: float(t[1]))[:6]:
        acc += float(nu); print("   %-24s r=%.6f nu=%.3e  cum=%.4f" % (str(lam), float(r), float(nu), acc))
    # mass within delta of 1
    for delta in [0.5, 0.2, 0.1, 0.05, 0.02, 0.01]:
        m1 = sum(float(nu) for _, r, nu in sm if float(r) >= 1 - delta)
        m2 = sum(float(nu) for _, r, nu in sm if float(r) <= -1 + delta)
        print("  delta=%.3f  nu[1-delta,1]=%.5f (ratio/delta=%.3f)   nu[-1,-1+delta]=%.3e" % (delta, m1, m1/delta, m2))

if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "brute": brute(int(sys.argv[2]))
    elif cmd == "table": table(int(sys.argv[2]), int(sys.argv[3]))
    elif cmd == "tail": tail(int(sys.argv[2]))
