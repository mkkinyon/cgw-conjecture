"""Remark 'reduced squares' check: enumerate all Latin squares with first
row = id at n=5,6; verify that (a) the type law of row 2 (= P_n), (b) the
type law of the row with first entry 2 (= second row of the reduced square
/ loop), and (c) the law of row 2 conditioned on sigma(1)=2, coincide as
exact rationals."""
import itertools
from collections import Counter
from fractions import Fraction

def cyc_type(p):
    n=len(p); seen=[False]*n; t=[]
    for i in range(n):
        if not seen[i]:
            l=0; j=i
            while not seen[j]: seen[j]=True; j=p[j]; l+=1
            t.append(l)
    return tuple(sorted(t))

def all_squares(n):
    def rec(rows, colsets):
        if len(rows)==n:
            yield rows; return
        for p in itertools.permutations(range(n)):
            if all(p[c] not in colsets[c] for c in range(n)):
                yield from rec(rows+[p], [colsets[c]|{p[c]} for c in range(n)])
    yield from rec([tuple(range(n))], [{i} for i in range(n)])

if __name__ == '__main__':
    for n in (5,6):
        P=Counter(); Ploop=Counter(); Pcond=Counter(); tot=0
        for sq in all_squares(n):
            tot+=1
            r2=sq[1]; P[cyc_type(r2)]+=1
            if r2[0]==1: Pcond[cyc_type(r2)]+=1
            loop_row=[r for r in sq[1:] if r[0]==1][0]
            Ploop[cyc_type(loop_row)]+=1
        ok=True
        for t in sorted(P):
            a=Fraction(P[t],tot); b=Fraction(Ploop[t],tot)
            c=Fraction(Pcond[t],sum(Pcond.values()))
            ok &= (a==b==c)
            print(f"  n={n} type {t}: P_n={a} loop={b} cond={c} "
                  f"{'MATCH' if a==b==c else 'DIFFER'}")
        assert ok
        print(f"n={n}: all three laws coincide exactly ({tot} squares)")
