#!/usr/bin/env python3
"""MC check: law of the size-biased block mass at time m (chain from one block) vs
rho_m(x) = 1-(1-2x)^m on (0,1), atom (1+(-1)^m)/(2(m+1)) at x=1; and E[sum p^r] vs
1/r + int_0^1 (1-x^{r-1})(1-2x)^m dx.   usage: python3 s10_sizebiased.py m samples"""
import sys, random
from s10_blocksim import simulate  # not used; own loop below
def run(m, S):
    bins = 10; hist=[0]*bins; atom=0; mom={2:0.0,3:0.0,4:0.0}
    for _ in range(S):
        blocks=[1.0]
        for t in range(m):
            c1=random.random(); c2=random.random(); acc=0.0; i1=i2=None
            for i,b in enumerate(blocks):
                if i1 is None and c1<acc+b: i1=i
                if i2 is None and c2<acc+b: i2=i
                acc+=b
            if i1 is None: i1=len(blocks)-1
            if i2 is None: i2=len(blocks)-1
            if i1==i2:
                b=blocks[i1]; mm=abs(random.random()-random.random())*b; blocks[i1]=b-mm; blocks.append(mm)
            else:
                b=blocks[i1]+blocks[i2]
                for i in sorted([i1,i2],reverse=True): del blocks[i]
                blocks.append(b)
        # size-biased pick
        c=random.random(); acc=0.0; x=blocks[-1]
        for b in blocks:
            if c<acc+b: x=b; break
            acc+=b
        if x>=1.0-1e-12: atom+=1
        else: hist[min(int(x*bins),bins-1)]+=1
        for r in mom: mom[r]+=sum(b**r for b in blocks)
    return hist, atom, mom
if __name__=='__main__':
    m=int(sys.argv[1]); S=int(sys.argv[2]); hist,atom,mom=run(m,S)
    from math import comb
    def integ(f,N=20000):
        return sum(f((i+0.5)/N) for i in range(N))/N
    print(f"m={m}: atom MC {atom/S:.5f}  exact {(1+(-1)**m)/(2*(m+1)):.5f}")
    for i in range(10):
        lo,hi=i/10,(i+1)/10
        ex=integ(lambda x: (1-(1-2*(lo+x*(hi-lo)))**m))*(hi-lo)
        print(f"  x in [{lo:.1f},{hi:.1f}): MC {hist[i]/S:.5f}  exact {ex:.5f}")
    for r in (2,3,4):
        ex=1/r+integ(lambda x:(1-x**(r-1))*(1-2*x)**m)
        print(f"  E[sum p^{r}]: MC {mom[r]/S:.5f}  exact {ex:.5f}")
