"""s13_cgw_altpair.py -- case 3 of CGW splitting: after the published cross-switch at {omega j, omega j'}, the pairs {j, omega^2 j} and
{j', omega^2 j'} both toggle status and both sit at separation (alpha, beta) in omega'; whichever was a B-pair in L gives a valid
splitting edge (backflip there lands in S(mu)).  Reports toggling, separation and type checks and the fraction of case-3 instances
repaired.  usage: ./jm n nsamp thin seed | python3 s13_cgw_altpair.py n"""
import sys; 
from s13_cgwcase3 import *
from s8_real_bias import read_squares
from collections import Counter
n=int(sys.argv[1]); st=Counter()
for L in read_squares(n, sys.stdin.buffer):
    om=omega_of(L); lam=ctype(L)
    for cyc in cycles(om):
        m=len(cyc)
        if m<6: continue
        for alpha in range(3,m-2):
            beta=m-alpha
            for i,j in enumerate(cyc):
                jp=cyc[(i+alpha)%m]
                if is_A(L,j,jp) or is_A(L,om[j],om[jp]): continue
                mu=sorted(lam); mu.remove(m); mu=tuple(sorted(mu+[alpha,beta]))
                e1,e2,e3,e4=L[0][om[j]],L[0][om[jp]],L[1][om[j]],L[1][om[jp]]
                first = min(e1,e2,e3,e4) in (e1,e4)
                L1=cross_switch(L,om[j],om[jp]); om1=omega_of(L1)
                # candidate pairs: first form -> {j, om^2 j} and {jp, om^2 jp}; second form -> {j, om^-2 j}, {jp, om^-2 jp}
                cands=[(j,om[om[j]]),(jp,om[om[jp]])]
                st['case3 instances']+=1; ok_any=False
                for (a,c) in cands:
                    wasB = not is_A(L,a,c)
                    nowA = is_A(L1,a,c)
                    st['toggled' if wasB==nowA else 'NOT toggled']+=1
                    d=1; x=om1[a]
                    while x!=c: x=om1[x]; d+=1
                    st[f'omega1-distance a->c is alpha' if d==alpha else f'distance other']+=1
                    if nowA:
                        L2=turn(L1,a,c,col_cycle_through(L1,a,c,1))
                        st['backflip at candidate lands in S(mu)' if ctype(L2)==mu else 'backflip WRONG']+=1
                        ok_any=True
                st['instance repaired (some candidate was B in L)' if ok_any else 'instance NOT repaired (both candidates A in L)']+=1
for k,v in sorted(st.items()): print(f'{v:7d}  {k}')
