"""Production layer-4 analysis via the pair-IE series, types (n) vs (2,n-2).

Outputs per n:
  exact T_0..T_3 for both types, t_j = T_j/T_0,
  truncated layer-4 effects eff(J), J=0..3,
  decomposition of r4 - r3 = Delta log Phi via partial sums,
  comparison with exact values (n <= 11).
"""
from fractions import Fraction
from math import factorial, comb, log
from functools import lru_cache
import sys, time
from puncture import make_sigma

@lru_cache(maxsize=None)
def path_mp(e):
    return tuple(comb(e - k + 1, k) for k in range((e + 1) // 2 + 1))

@lru_cache(maxsize=None)
def cycle_mp(m):
    return tuple(m * comb(m - k, k) // (m - k) for k in range(m // 2 + 1))

@lru_cache(maxsize=None)
def count_from_sig(sig, N):
    poly = [1]
    for typ, m in sig:
        q = cycle_mp(m) if typ == 1 else path_mp(m)
        out = [0] * (len(poly) + len(q) - 1)
        for i, x in enumerate(poly):
            if x:
                for j, y in enumerate(q):
                    out[i + j] += x * y
        poly = out
    return sum((-1) ** k * poly[k] * factorial(N - k) for k in range(min(len(poly), N + 1)))

def make_ctx(n, lam):
    sigma = make_sigma(n, lam)
    sig_inv = [0] * n
    for i in range(n):
        sig_inv[sigma[i]] = i
    cells = [(i, c) for i in range(n) for c in range(n) if c != i and c != sigma[i]]
    return sigma, sig_inv, cells

def fast_R(n, sigma, sig_inv, M):
    rows_rm = frozenset(i for i, c in M)
    cols_rm = frozenset(c for i, c in M)
    seen = [False] * (2 * n)
    comps = []
    for s in range(2 * n):
        if seen[s]:
            continue
        if (s < n and s in rows_rm) or (s >= n and (s - n) in cols_rm):
            seen[s] = True
            continue
        stack = [s]
        seen[s] = True
        nv = 0
        ecount = 0
        while stack:
            u = stack.pop()
            nv += 1
            if u < n:
                nb = []
                if u not in cols_rm:
                    nb.append(n + u)
                if sigma[u] not in cols_rm:
                    nb.append(n + sigma[u])
            else:
                c = u - n
                nb = []
                if c not in rows_rm:
                    nb.append(c)
                if sig_inv[c] not in rows_rm:
                    nb.append(sig_inv[c])
            ecount += len(nb)
            for w in nb:
                if not seen[w]:
                    seen[w] = True
                    stack.append(w)
        ne = ecount // 2
        if ne:
            comps.append((1 if ne == nv else 0, ne))
    return count_from_sig(tuple(sorted(comps)), n - len(M))

def Ts(n, lam, jmax=3):
    sigma, sig_inv, cells = make_ctx(n, lam)
    T = [0] * (jmax + 1)
    T[0] = fast_R(n, sigma, sig_inv, ()) ** 2
    if jmax >= 1:
        for cell in cells:
            T[1] += fast_R(n, sigma, sig_inv, (cell,)) ** 2
    if jmax >= 2:
        L = len(cells)
        for a in range(L):
            ia, ca = cells[a]
            for b in range(a + 1, L):
                ib, cb = cells[b]
                if ia != ib and ca != cb:
                    T[2] += fast_R(n, sigma, sig_inv, (cells[a], cells[b])) ** 2
    if jmax >= 3:
        L = len(cells)
        for a in range(L):
            ia, ca = cells[a]
            for b in range(a + 1, L):
                ib, cb = cells[b]
                if ia == ib or ca == cb:
                    continue
                for d in range(b + 1, L):
                    id_, cd = cells[d]
                    if id_ in (ia, ib) or cd in (ca, cb):
                        continue
                    T[3] += fast_R(n, sigma, sig_inv,
                                   (cells[a], cells[b], cells[d])) ** 2
    return T

EXACT = {
    8: (4738, 5955474, 4740, 5960976),
    9: (43387, 522417852, 43400, 522738160),
    10: (439792, 55553633070, 439872, 55574137280),
    11: (4890741, 7058091360392, 4891320, 7059777493824),
}

def run(n, jmax=3):
    t0 = time.time()
    Ta = Ts(n, (n,), jmax)       # type 0
    Tb = Ts(n, (2, n - 2), jmax)  # type 1
    N3a = int(Ta[0] ** 0.5 + .5)
    # use exact sqrt via integer isqrt
    import math
    N3a = math.isqrt(Ta[0]); N3b = math.isqrt(Tb[0])
    assert N3a * N3a == Ta[0] and N3b * N3b == Tb[0]
    r3 = Fraction(N3b, N3a) - 1
    print(f"n={n}  [{time.time()-t0:.0f}s]  r3 = {float(r3):.6e}")
    effs = []
    for J in range(jmax + 1):
        Sa = sum((-1) ** j * Ta[j] for j in range(J + 1))
        Sb = sum((-1) ** j * Tb[j] for j in range(J + 1))
        eff = Fraction(Sb * N3a, Sa * N3b) - 1
        effs.append(eff)
        line = f"  eff(J={J}) = {float(eff):+.6e}"
        if n in EXACT:
            n30, n40, n31, n41 = EXACT[n]
            true_eff = Fraction(n41 * n30, n40 * n31) - 1
            line += f"   exact = {float(true_eff):+.6e}  ratio {float(eff/true_eff):.4f}"
        print(line)
    # decomposition of r4 - r3 (Delta log Phi) via truncations
    for J in range(1, jmax + 1):
        Pa = sum((-1) ** j * Fraction(Ta[j], Ta[0]) for j in range(J + 1))
        Pb = sum((-1) ** j * Fraction(Tb[j], Tb[0]) for j in range(J + 1))
        rat = Pb / Pa
        line = (f"  Dlog Phi_J (J={J}) = {log(float(rat)):+.3e}" if rat > 0
                else f"  Dlog Phi_J (J={J}) = n/a (truncation sign change)")
        if n in EXACT:
            n30, n40, n31, n41 = EXACT[n]
            Phia = Fraction(n40, n30 * n30)
            Phib = Fraction(n41, n31 * n31)
            line += f"   exact DlogPhi = {log(float(Phib)) - log(float(Phia)):+.3e}"
        print(line)
    return effs

if __name__ == "__main__":
    ns = [int(a) for a in sys.argv[1:]] or [8, 9, 10, 11, 12, 13, 14]
    for n in ns:
        run(n)
