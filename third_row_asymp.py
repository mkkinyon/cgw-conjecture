"""Pin down the asymptotics of the third-row-layer effects:
    rho_2(n) = R3((2,n-2))/R3((n,)) - 1   (conjectured ~ a2/n^4 * (1+b2/n+...))
    rho_3(n) = R3((3,n-3))/R3((n,)) - 1   (conjectured ~ a3/n^6 * ...)
computed EXACTLY with integer arithmetic for n up to N, then Richardson-
extrapolated in powers of 1/n.
"""
from math import comb, factorial
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 80


def cycle_mp(m):
    return [m * comb(m - k, k) // (m - k) for k in range(m // 2 + 1)]


def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def R3_from_poly(mk, n):
    # sum (-1)^k m_k (n-k)!  -- compute with running factorial
    tot = 0
    f = factorial(n)
    for k, m in enumerate(mk):
        if k > 0:
            f //= (n - k + 1)
        tot += (-1) ** k * m * f
    return tot


def rho(n, ell):
    base = cycle_mp(2 * n)
    split = mul(cycle_mp(2 * ell), cycle_mp(2 * (n - ell)))
    Rb = R3_from_poly(base, n)
    Rs = R3_from_poly(split, n)
    return Fraction(Rs, Rb) - 1


def richardson(ns, vals, order):
    """Fit vals(n) = a0 + a1/n + ... + a_{order}/n^order using the last
    order+1 points, solve exactly-ish with mpmath lu_solve."""
    pts = list(zip(ns, vals))[-(order + 1):]
    A = mp.matrix([[mp.mpf(1) / mp.mpf(n) ** j for j in range(order + 1)] for n, _ in pts])
    b = mp.matrix([v for _, v in pts])
    sol = mp.lu_solve(A, b)
    return sol


if __name__ == "__main__":
    ns = list(range(100, 401, 50))
    for ell, power in [(2, 4), (3, 6)]:
        vals = []
        for n in ns:
            r = rho(n, ell)
            v = mp.mpf(r.numerator) / mp.mpf(r.denominator) * mp.mpf(n) ** power
            vals.append(v)
            print(f"ell={ell} n={n:4d}  rho*n^{power} = {mp.nstr(v, 12)}")
        for order in [2, 3, 4, 5]:
            sol = richardson(ns, vals, order)
            print(f"  ell={ell}: extrapolation order {order}: "
                  + ", ".join(mp.nstr(sol[j], 8) for j in range(min(order + 1, 4))))
