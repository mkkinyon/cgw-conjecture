"""Truncated pair-IE analysis of the layer-4 effect.

For sigma of types (n) and (2, n-2):
  T_j(sigma), j = 0..3  (exact),
  partial sums S_J = sum_{j<=J} (-1)^j T_j,
  the truncated layer-4 effect
     eff(J) = [S_J((2,n-2)) * N3((n))] / [S_J((n)) * N3((2,n-2))] - 1,
compared (n <= 11) with the exact effect from layer4.c data.
"""
from math import factorial
from puncture import punctured_count, make_sigma
from pair_ie import Tj
import sys, time

EXACT = {  # n: (N3_0, N4_0, N3_1, N4_1)   0 = (n,), 1 = (2, n-2)
    8: (4738, 5955474, 4740, 5960976),
    9: (43387, 522417852, 43400, 522738160),
    10: (439792, 55553633070, 439872, 55574137280),
    11: (4890741, 7058091360392, 4891320, 7059777493824),
}

def analyse(n, jmax=3):
    out = {}
    for tag, lam in (("0", (n,)), ("1", (2, n - 2))):
        sigma = make_sigma(n, lam)
        T = [Tj(sigma, j) for j in range(jmax + 1)]
        out[tag] = T
    from fractions import Fraction
    T0 = out["0"]; T1 = out["1"]
    N3_0 = punctured_count(make_sigma(n, (n,)), [])
    N3_1 = punctured_count(make_sigma(n, (2, n - 2)), [])
    print(f"n={n}:")
    print(f"  N3: {N3_0} / {N3_1}   t_j ratios type0: " +
          ", ".join(f"{float(Fraction(T0[j], T0[0])):.5f}" for j in range(1, jmax+1)))
    for J in range(jmax + 1):
        S0 = sum((-1) ** j * T0[j] for j in range(J + 1))
        S1 = sum((-1) ** j * T1[j] for j in range(J + 1))
        eff = Fraction(S1 * N3_0, S0 * N3_1) - 1
        line = f"  J={J}: eff = {float(eff):+.6e}"
        if n in EXACT:
            n30, n40, n31, n41 = EXACT[n]
            true_eff = Fraction(n41 * n30, n40 * n31) - 1
            line += f"   exact = {float(true_eff):+.6e}   (J-est/exact = {float(eff/true_eff):.4f})"
        print(line)
    if n in EXACT:
        n30, n40, n31, n41 = EXACT[n]
        # also check absolute truncation: S_J vs true N4
        S3_0 = sum((-1) ** j * T0[j] for j in range(jmax + 1))
        print(f"  S3(type0)/N4(type0) = {S3_0 / n40:.6f}")

if __name__ == "__main__":
    ns = [int(a) for a in sys.argv[1:]] or [8, 9, 10, 11]
    for n in ns:
        t0 = time.time()
        analyse(n)
        print(f"  [{time.time()-t0:.1f}s]")
