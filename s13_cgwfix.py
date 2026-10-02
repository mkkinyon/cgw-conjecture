"""s13_cgwfix.py -- test of a repair of CGW case 3 when alpha is even.
For (j,j') B at distance alpha (alpha even) on an (alpha+beta)-cycle, let R = (omega^{alpha/2} j, omega^{3alpha/2} j):
the pair symmetric about j' at distance alpha.  If R is A: backflip at R (own edge, type mu).  (requires alpha < 2 beta so that the arc of length alpha around j' avoids j.)  If R is B:
cross-switch at R (its arc has j' at the midpoint, j outside), which toggles (j,j') and keeps its distance
alpha; then backflip at (j,j').  Check the resulting type in every instance.
usage: ./jm n nsamp thin seed | python3 s13_cgwfix.py n
"""
import sys
from collections import Counter
from s8_real_bias import read_squares
from s13_cgwcase3 import omega_of, cycles, ctype, is_latin, col_cycle_through, turn, is_A, cross_switch

def main():
    n = int(sys.argv[1]); st = Counter()
    for L in read_squares(n, sys.stdin.buffer):
        om = omega_of(L); lam = ctype(L)
        for cyc in cycles(om):
            m = len(cyc)
            if m < 5: continue
            for alpha in range(2, m - 1):
                if alpha % 2: continue
                beta = m - alpha
                if 3 * alpha // 2 >= m and False: pass
                mu = sorted(list(lam)); mu.remove(m); mu = tuple(sorted(mu + [alpha, beta]))
                for i, j in enumerate(cyc):
                    jp = cyc[(i + alpha) % m]
                    if is_A(L, j, jp): continue
                    r1, r2 = cyc[(i + alpha // 2) % m], cyc[(i + 3 * alpha // 2) % m]
                    if 3 * alpha // 2 >= m: st['not covered: alpha >= 2 beta (arc around jp reaches j)'] += 1; continue
                    if is_A(L, r1, r2):
                        L2 = turn(L, r1, r2, col_cycle_through(L, r1, r2, 1)); key = 'R is A: backflip at R'
                    else:
                        L1 = cross_switch(L, r1, r2)
                        if not (is_latin(L1) and ctype(L1) == lam): st['cross-switch FAILS'] += 1; continue
                        if not is_A(L1, j, jp): st['R is B: (j,jp) NOT toggled'] += 1; continue
                        L2 = turn(L1, j, jp, col_cycle_through(L1, j, jp, 1)); key = 'R is B: cross-switch at R, backflip at (j,jp)'
                    st[key + (' -> type mu' if ctype(L2) == mu else ' -> WRONG type')] += 1
    for k, v in sorted(st.items()): print(f"{v:7d}  {k}")

if __name__ == "__main__": main()
