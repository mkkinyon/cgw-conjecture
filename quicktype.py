import sys
from collections import Counter
n=int(sys.argv[1])
cnt = Counter()
buf = sys.stdin.buffer
while True:
    sq = buf.read(n*n)
    if len(sq)<n*n: break
    L=[sq[r*n:(r+1)*n] for r in range(n)]
    for r in range(n): assert sorted(L[r])==list(range(n))
    for c in range(n): assert sorted(L[r][c] for r in range(n))==list(range(n))
    sig=[0]*n
    for j in range(n): sig[L[0][j]]=L[1][j]
    seen=[False]*n; lens=[]
    for i in range(n):
    	if not seen[i]:
            l=0; x=i
            while not seen[x]: seen[x]=True; x=sig[x]; l+=1
            lens.append(l)
    cnt[tuple(sorted(lens))]+=1
tot=sum(cnt.values())
for k,v in sorted(cnt.items()):
    print(k, v, f"{v/tot:.5f}")
