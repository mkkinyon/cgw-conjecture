"""s13_cgw_multigraph.py -- exact comparison of CGW08's multigraphs G_S (splitting, p.294) and G_J (joining, p.295)
for lambda=(n) -> mu=(alpha, n-alpha), F = empty, over ALL Latin squares of order n (unnormalised: the switch and the
cross-switch choose their form by the smallest symbol, so symbol relabelling is not a symmetry of the procedures).
The proof of Lemma 3.12 asserts G_S = G_J and that every mu-square has degree in [alpha*beta, 3*alpha*beta].
usage: python3 s13_cgw_multigraph.py n alpha        (n=4: (4)->(2,2) agree; n=5: (5)->(2,3) differ; ~6 min at n=5)
switch = CGW Def 3.6 (trade the row cycle with the smaller symbol); joining_edges = rules 1-4 of p.295."""
import sys
sys.path.insert(0,'/home/claude/cgw'); 
from collections import Counter
from s13_cgwcase3 import omega_of, cycles, ctype, col_cycle_through, turn, is_A, cross_switch
def switch(L, j, jp):
    # Def 3.6: trade the row cycle (rows 1,2) through whichever of j, jp has the smaller symbol-min
    n = len(L); om = omega_of(L)
    def rowcycle(c):
        cyc=[c]; x=om[c]
        while x!=c: cyc.append(x); x=om[x]
        return cyc
    Q1 = rowcycle(j); Q2 = rowcycle(jp)
    m1 = min(min(L[0][c],L[1][c]) for c in Q1); m2 = min(min(L[0][c],L[1][c]) for c in Q2)
    Q = Q1 if m1 < m2 else Q2
    M=[list(r) for r in L]
    for c in Q: M[0][c],M[1][c]=L[1][c],L[0][c]
    return M

def joining_edges(L2, j, jp):
    """CGW joining rules 1-4 from L2 at the pair {j,jp} (j,jp on different cycles).  Returns list of squares."""
    out=[]
    L2b = L2
    if not is_A(L2, j, jp): L2b = switch(L2, j, jp)
    assert is_A(L2b, j, jp)
    Lp = turn(L2b, j, jp, col_cycle_through(L2b, j, jp, 1))   # flip: cycle through row 2
    out.append(Lp)
    omp = omega_of(Lp)
    inv = [0]*len(L2); 
    for c in range(len(L2)): inv[omp[c]] = c
    if not is_A(Lp, inv[j], inv[jp]): out.append(Lp)          # rule 3
    if not is_A(Lp, omp[j], omp[jp]): out.append(cross_switch(Lp, omp[j], omp[jp]))  # rule 4
    return out


def splitting_edges(L, alpha, beta):
    om=omega_of(L); out=[]
    for cyc in cycles(om):
        m=len(cyc)
        if m!=alpha+beta: continue
        pairs=set()
        for i,j in enumerate(cyc):
            jp=cyc[(i+alpha)%m]; pairs.add(frozenset((j,jp)))
        for pr in pairs:
            j,jp=tuple(pr)
            if is_A(L,j,jp): L2=turn(L,j,jp,col_cycle_through(L,j,jp,1)); case=1
            elif is_A(L,om[j],om[jp]): L2=turn(L,om[j],om[jp],col_cycle_through(L,om[j],om[jp],1)); case=2
            else:
                L1=cross_switch(L,om[j],om[jp]); L2=turn(L1,j,jp,col_cycle_through(L1,j,jp,1)); case=3
            out.append((L2,case)); out.append((switch(L2,j,jp),5))
    return out

def main():
    n=int(sys.argv[1]); alpha=int(sys.argv[2]); beta=n-alpha
    lam=(n,); mu=tuple(sorted((alpha,beta)))
    # enumerate all Latin squares with first row identity
    squares=[]
    def rec(L, r):
        if r==n: squares.append([row[:] for row in L]); return
        used_col=[set(L[i][c] for i in range(r)) for c in range(n)]
        def fill(c, row):
            if c==n:
                L.append(row[:]); rec(L, r+1); L.pop(); return
            for s in range(n):
                if s not in row and s not in used_col[c]:
                    row.append(s); fill(c+1,row); row.pop()
        fill(0,[])
    rec([list(range(n))],1)
    from itertools import permutations
    allsq=[]
    for L in squares:
        for perm in permutations(range(n)):
            allsq.append([[perm[s] for s in row] for row in L])
    squares=allsq
    print("squares (all, unnormalised):",len(squares))
    key=lambda L: tuple(map(tuple,L))
    Sl=[L for L in squares if ctype(L)==lam]; Sm=[L for L in squares if ctype(L)==mu]
    print("|S(lam)|=",len(Sl),"|S(mu)|=",len(Sm), "ratio C(lam)/C(mu)*gamma... |S(lam)|/|S(mu)| =", len(Sl)/len(Sm))
    GS=Counter(); GJ=Counter()
    for L in Sl:
        for (M,case) in splitting_edges(L,alpha,beta):
            assert ctype(M)==mu, (ctype(M),case)
            GS[(key(L),key(M))]+=1
    for M in Sm:
        om=omega_of(M); cyc=cycles(om)
        ca=[c for c in cyc if len(c)==alpha]; cb=[c for c in cyc if len(c)==beta]
        for ia,A in enumerate(ca):
            for ib,B in enumerate(cb):
                if A is B: continue
                if alpha==beta and ib<ia: continue
                for j in A:
                    for jp in B:
                        for Lp in joining_edges(M,j,jp):
                            assert ctype(Lp)==lam
                            GJ[(key(Lp),key(M))]+=1
    print("edges G_S:",sum(GS.values()),"G_J:",sum(GJ.values()))
    print("G_S == G_J as multisets:", GS==GJ)
    only_S=sum((GS-GJ).values()); only_J=sum((GJ-GS).values())
    print("edges in G_S not in G_J (with multiplicity):",only_S," in G_J not in G_S:",only_J)
    degS=Counter(); degJ=Counter()
    for (a,b),v in GS.items(): degS[b]+=v
    for (a,b),v in GJ.items(): degJ[b]+=v
    print("G_S degrees of mu-squares:",Counter(degS[key(M)] for M in Sm))
    print("G_J degrees of mu-squares:",Counter(degJ[key(M)] for M in Sm))
    dl=Counter()
    for (a,b),v in GS.items(): dl[a]+=v
    print("G_S degrees of lam-squares:",Counter(dl[key(L)] for L in Sl))
    print("sum degS:",sum(degS.values()),"distinct targets in G_S:",len(degS),"distinct mu keys:",len(set(key(M) for M in Sm)))
    tS=set(b for (a,b) in GS); tM=set(key(M) for M in Sm)
    print("targets of G_S that are mu-squares:",len(tS&tM),"targets not in S(mu):",len(tS-tM))
    bad=[b for b in tS-tM][:1]
    if bad:
        B=[list(r) for r in bad[0]]; print("example target:",B,"type",ctype(B),"first row",B[0])

    bycase=Counter(); tot=Counter()
    for L in Sl:
        for (M,case) in splitting_edges(L,alpha,beta):
            tot[case]+=1
            if GJ[(key(L),key(M))]==0: bycase[case]+=1
    print("G_S edges by case:",dict(tot)," of which with multiplicity 0 in G_J:",dict(bycase))


if __name__ == '__main__':
    main()
