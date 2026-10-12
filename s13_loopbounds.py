"""s13_loopbounds.py -- the finite computations quoted in the sharpness remark of the loops note.

For composite n and each divisor a of n with 1 < a < n (b = n/a), the imprimitive union bound of the proof has the term
    n! * N(a) * (b!)^n * L(a)^{b^2} / L(n),
where N(a) is the number of partitions of [n] into b blocks of size a.  Two versions of the sum over a are computed:
  A ('crude'): N(a) replaced by n!, log L(a) by the upper bound and log L(n) by the lower bound of eq. (1) of the note
     (log L(a) <= a^2 log a - 2a^2 + a f(a), f(a) = log^2(a)/2 + 2 log a + 3; log L(n) >= n^2 log n - 2n^2); n!, (b!)^n exact.
  B ('exact'): N(a) = n!/((a!)^b b!), L(a) bounded by Bregman's product prod_{k<=a} (k!)^{a/k}, L(n) by van der Waerden's
     (n!)^{2n}/n^{n^2}.
The script reports, for each version, the largest composite n at which the sum is >= 1 (so the sum is < 1 for every larger
composite n up to NMAX) and the exponents at n = 30, 100, 1000 in units of n^2 log 2.
Analytic tail: for every n >= NMAX each term of version A (hence of B, which is termwise smaller) is less than 1/n by the two
inequalities used in the proof of the theorem,
  a >= sqrt n:  log term <= -n^2 log 2 + n^{3/2}(log^2(n)/2 + 3 log n + 3) + 2 n log n,
  a <  sqrt n:  log term <= -(n^2/4) log n + (5/2) n^2 + 2 n log n,
and the script checks that both right-hand sides are < -log n at n = NMAX and that both are decreasing in n beyond it (their
quotients by n^2 are decreasing), so the sums are < 1 for every composite n > NMAX as well.
usage: python3 s13_loopbounds.py [NMAX]
"""
import sys
from math import lgamma, log, sqrt

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 25000


def lf(k):                       # log k!
    return lgamma(k + 1)


def f(a):
    return 0.5 * log(a) ** 2 + 2 * log(a) + 3


# cumulative Bregman sums S[a] = sum_{k<=a} log(k!)/k, so that log prod_{k<=a}(k!)^{a/k} = a S[a]
S = [0.0] * (NMAX + 1)
for k in range(1, NMAX + 1):
    S[k] = S[k - 1] + lf(k) / k


def divisors(n):
    ds = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            ds.append(d)
            if d * d != n: ds.append(n // d)
        d += 1
    return sorted(ds)


def sums(n):
    ds = divisors(n)
    if not ds: return None
    logA = []; logB = []
    for a in ds:
        b = n // a
        A = lf(n) + lf(n) + n * lf(b) + b * b * (a * a * log(a) - 2 * a * a + a * f(a)) - (n * n * log(n) - 2 * n * n)
        B = lf(n) + (lf(n) - b * lf(a) - lf(b)) + n * lf(b) + b * b * (a * S[a]) - (2 * n * lf(n) - n * n * log(n))
        logA.append(A); logB.append(B)
    def logsum(xs):
        m = max(xs)
        return m + log(sum(2.718281828459045 ** (x - m) for x in xs))
    return logsum(logA), logsum(logB)


lastA = lastB = None
report = {}
for n in range(4, NMAX + 1):
    r = sums(n)
    if r is None: continue
    A, B = r
    if A >= 0: lastA = n
    if B >= 0: lastB = n
    if n in (30, 100, 1000): report[n] = (A / (n * n * log(2)), B / (n * n * log(2)))
print(f'NMAX={NMAX}')
print(f'version A (eq.(1), crude n! partitions): sum >= 1 last at composite n={lastA}; '
      f'so < 1 for every composite n >= {lastA + 1 if lastA else "?"} up to NMAX')
print(f'version B (exact vdW/Bregman, exact partition count): sum >= 1 last at composite n={lastB}; '
      f'so < 1 for every composite n >= {lastB + 1 if lastB else "?"} up to NMAX')
for n, (a, b) in sorted(report.items()):
    print(f'  n={n}: log(sum)/(n^2 log 2): A={a:.4f}  B={b:.4f}')
# analytic tail
def tail1(n): return -n * n * log(2) + n ** 1.5 * (0.5 * log(n) ** 2 + 3 * log(n) + 3) + 2 * n * log(n)
def tail2(n): return -(n * n / 4) * log(n) + 2.5 * n * n + 2 * n * log(n)
n = NMAX
print(f'analytic tail at n=NMAX: (a>=sqrt n) bound/(n^2)={tail1(n)/n/n:.5f}, (a<sqrt n) bound/(n^2)={tail2(n)/n/n:.5f}; '
      f'need both < -log(n)/n^2 = {-log(n)/n/n:.2e}: {"OK" if tail1(n) < -log(n) and tail2(n) < -log(n) else "FAIL"}')
# monotonicity beyond NMAX: the quotients by n^2 are -log2 + n^{-1/2}(log^2 n/2 + 3 log n + 3) + 2 log n / n and
# -(log n)/4 + 5/2 + 2 log n / n; both decreasing for n >= NMAX (checked on a grid up to 10^9, and both terms with n^{-1/2}
# or 1/n are decreasing once log n >= 6).
ok = True; prev1 = tail1(n) / n / n; prev2 = tail2(n) / n / n
m = n
while m < 10 ** 9:
    m = int(m * 1.1) + 1
    t1 = tail1(m) / m / m; t2 = tail2(m) / m / m
    ok = ok and t1 < prev1 and t2 < prev2
    prev1, prev2 = t1, t2
print(f'quotients decreasing on a geometric grid from NMAX to 1e9: {ok}')
# smallest n from which the tail bounds alone suffice
nn = 1000
while not (tail1(nn) < -log(nn) and tail2(nn) < -log(nn)): nn += 1
print(f'tail bounds alone give sum < 1 for all n >= {nn} (given monotonicity)')
